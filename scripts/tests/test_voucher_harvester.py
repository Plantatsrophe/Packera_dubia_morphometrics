#!/usr/bin/env python3
"""
===============================================================================
Unit tests for botanical voucher harvesting, authority stratification,
Darwin Core metadata normalization, and atomic persistence.
===============================================================================
"""

import sys
import tempfile
import unittest

from scripts.tests._r_session import rscript_available, run_warm_r_cmd
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
from PIL import Image

from scripts.core.config import (
    DEFAULT_MIN_MEGAPIXELS,
    EXCLUDED_WESTERN_STATES,
    WESTERN_LONGITUDE_THRESHOLD,
)
from scripts.core.taxonomy_authority import EXPORT_COLUMNS
from scripts.core.harvester import (
    VoucherHarvester,
    export_curated_table,
    extract_high_res_image_url,
    generate_exsiccatae_key,
    harvest_taxa_occurrences,
    haversine_distance_meters,
    infer_regional_group,
    is_excluded_western_region,
    normalize_collector,
    normalize_event_date,
    normalize_record_number,
    optimize_herbarium_image_url,
    parse_determiner_tier,
    sanitize_filename,
    spatially_balance_and_thin_vouchers,
    stratify_duplicates,
    validate_image_quality,
)


class TestVoucherHarvesterGeographicFiltering(unittest.TestCase):
    """Test suite for western region exclusion and geographic filtering."""

    def test_western_states_exclusion_full_names(self):
        """Verify that all full names of states farther west than TX/OK are excluded."""
        western_names = [
            "Colorado", "New Mexico", "Wyoming", "Montana", "Utah",
            "Idaho", "Arizona", "Nevada", "Washington", "Oregon",
            "California", "Alaska", "Hawaii",
        ]
        for state in western_names:
            with self.subTest(state=state):
                self.assertTrue(
                    is_excluded_western_region(state),
                    f"State '{state}' should be excluded.",
                )

    def test_western_states_exclusion_abbreviations(self):
        """Verify that 2-letter postal codes for western states are excluded."""
        western_codes = [
            "CO", "NM", "WY", "MT", "UT", "ID", "AZ", "NV", "WA", "OR", "CA", "AK", "HI",
        ]
        for code in western_codes:
            with self.subTest(code=code):
                self.assertTrue(
                    is_excluded_western_region(code),
                    f"Code '{code}' should be excluded.",
                )

    def test_eastern_and_plains_states_retained(self):
        """Verify that Texas, Oklahoma, Great Plains, and Eastern states are NOT excluded."""
        retained_states = [
            "Texas", "TX",
            "Oklahoma", "OK",
            "Kansas", "KS",
            "Nebraska", "NE",
            "South Dakota", "SD",
            "North Dakota", "ND",
            "North Carolina", "NC", "North Carolina (State)",
            "Virginia", "VA",
            "Georgia", "GA",
            "Wisconsin", "WI", "Wisconsin (State)",
            "Florida", "FL",
            "Louisiana", "LA",
            "Arkansas", "AR",
            "Missouri", "MO",
            "Illinois", "IL",
            "Michigan", "MI",
            "Minnesota", "MN",
            "Pennsylvania", "PA",
            "New York", "NY",
            "Maine", "ME",
            "Ontario",
        ]
        for state in retained_states:
            with self.subTest(state=state):
                self.assertFalse(
                    is_excluded_western_region(state),
                    f"State '{state}' should be retained (not excluded).",
                )

    def test_washington_dc_disambiguation(self):
        """Verify that Washington, D.C. is retained while Washington State is excluded."""
        self.assertTrue(is_excluded_western_region("Washington"))
        self.assertTrue(is_excluded_western_region("WA"))
        self.assertFalse(is_excluded_western_region("Washington, D.C."))
        self.assertFalse(is_excluded_western_region("Washington D.C."))
        self.assertFalse(is_excluded_western_region("Washington DC"))
        self.assertFalse(is_excluded_western_region("District of Columbia"))
        self.assertFalse(is_excluded_western_region("DC"))

    def test_coordinate_fallback_filtering(self):
        """Verify longitude threshold fallback when stateProvince is missing/unrecorded."""
        self.assertTrue(
            is_excluded_western_region(None, lat=39.5, lon=-108.5),
            "Unrecorded state with longitude < -106.65 should be excluded.",
        )
        self.assertFalse(
            is_excluded_western_region(None, lat=35.9, lon=-79.0),
            "Unrecorded state with Eastern longitude should be retained.",
        )
        self.assertFalse(
            is_excluded_western_region(None, lat=31.5, lon=-98.0),
            "Unrecorded state in central Texas longitude should be retained.",
        )

    def test_harvester_ingestion_drops_western_records_and_retains_dwc_fields(self):
        """Verify harvester normalizes raw DwC temporal fields and excludes western records."""
        mock_results = [
            {
                "key": 1001,
                "scientificName": "Packera paupercula",
                "stateProvince": "Colorado",
                "decimalLatitude": 39.18,
                "decimalLongitude": -106.05,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2020, "month": 6, "day": 15,
                "media": [{"identifier": "https://example.com/colo1.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "COLO001",
                "institutionCode": "COLO",
            },
            {
                "key": 1002,
                "scientificName": "Packera paupercula",
                "stateProvince": "North Carolina",
                "decimalLatitude": 35.90,
                "decimalLongitude": -79.05,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2020, "month": 5, "day": 10,
                "eventDate": "2020-05-10",
                "media": [{"identifier": "https://example.com/ncu1.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "NCU001",
                "institutionCode": "NCU",
                "identifiedBy": "Debra Trock",
            },
            {
                "key": 1003,
                "scientificName": "Packera paupercula",
                "stateProvince": "Texas",
                "decimalLatitude": 30.50,
                "decimalLongitude": -97.50,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2021, "month": 4, "day": 20,
                "media": [{"identifier": "https://example.com/tex1.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "TEX001",
                "institutionCode": "TEX",
            },
            {
                "key": 1004,
                "scientificName": "Packera paupercula",
                "stateProvince": "Wyoming",
                "decimalLatitude": 43.00,
                "decimalLongitude": -108.00,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2022, "month": 7, "day": 4,
                "media": [{"identifier": "https://example.com/rm1.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "RM001",
                "institutionCode": "RM",
            },
        ]

        with patch("scripts.core.harvester.occ.search") as mock_occ_search:
            mock_occ_search.return_value = {
                "results": mock_results,
                "count": len(mock_results),
            }

            df = harvest_taxa_occurrences(
                taxa_list=["Packera paupercula"],
                max_records_per_taxon=10,
                exclude_western=True,
            )

            # Only NCU001 (NC) and TEX001 (TX) should be retained
            self.assertEqual(len(df), 2)
            retained_catalogs = set(df["catalogNumber"])
            self.assertIn("NCU001", retained_catalogs)
            self.assertIn("TEX001", retained_catalogs)
            self.assertNotIn("COLO001", retained_catalogs)
            self.assertNotIn("RM001", retained_catalogs)

            # Check raw Darwin Core temporal fields are present
            self.assertIn("year", df.columns)
            self.assertIn("month", df.columns)
            self.assertIn("day", df.columns)
            self.assertIn("eventDate", df.columns)

            # Verify derived calendar day-of-year is decoupled from Step 01 curated export
            self.assertNotIn("doy", df.columns)

            # Check Determiner Tier Stratification (NCU001 was identified by Debra Trock -> Tier_1_Gold)
            ncu_rec = df[df["catalogNumber"] == "NCU001"].iloc[0]
            self.assertEqual(ncu_rec["determiner_tier"], "Tier_1_Gold")


class TestDeterminerAuthorityStratification(unittest.TestCase):
    """Test suite for 3-tier determiner authority classification."""

    def test_tier_1_specialist(self):
        """Specialist annotations (Trock, Barkley, Kowal) map to Tier_1_Gold."""
        tier, t_clean, determiner = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw="D.K. Trock",
            recorded_by_raw="John Doe",
            history_raw=None,
            institution_code_raw="NCU",
            locality_raw="Granite outcrop near Wake Forest",
            habitat_raw="Flatrock",
        )
        self.assertEqual(tier, "Tier_1_Gold")
        self.assertEqual(determiner, "D.K. Trock")

    def test_tier_1_type_specimen(self):
        """Holotypes and isotypes map to Tier_1_Gold regardless of determiner."""
        tier, t_clean, determiner = parse_determiner_tier(
            type_status_raw="Holotype",
            identified_by_raw=None,
            recorded_by_raw="A. Radford",
            history_raw=None,
            institution_code_raw="NCU",
            locality_raw="Rich mountain slope",
            habitat_raw="Cove forest",
        )
        self.assertEqual(tier, "Tier_1_Gold")
        self.assertEqual(t_clean, "Holotype")

    def test_tier_2_major_herbarium_rich_locality(self):
        """Major herbarium with rich locality and verified determiner maps to Tier_2_Silver."""
        tier, t_clean, determiner = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw="Botanist Jane",
            recorded_by_raw="Jane Collector",
            history_raw=None,
            institution_code_raw="NCU",
            locality_raw="Roadside embankment 5 miles north of Chapel Hill",
            habitat_raw="Granitic glade",
        )
        self.assertEqual(tier, "Tier_2_Silver")

    def test_tier_3_bronze_unverified(self):
        """General collectors or unknown herbaria without rich ecology map to Tier_3_Bronze."""
        tier, t_clean, determiner = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw=None,
            recorded_by_raw="Unknown Collector",
            history_raw=None,
            institution_code_raw="XYZ_HERB",
            locality_raw="Roadside",
            habitat_raw=None,
        )
        self.assertEqual(tier, "Tier_3_Bronze")


class TestImageQualityAndUrlOptimization(unittest.TestCase):
    """Test suite for URL rewriting, media ranking, and optical image quality validation."""

    def test_optimize_smithsonian_nmnh_urls(self):
        """Verify that dimension clamps (&h=2000) are removed from Smithsonian NMNH URLs."""
        raw_url = "https://collections.nmnh.si.edu/media/?i=11709418&h=2000"
        optimized = optimize_herbarium_image_url(raw_url)
        self.assertEqual(optimized, "https://collections.nmnh.si.edu/media/?i=11709418")

    def test_optimize_symbiota_urls(self):
        """Verify that web/thumbnail paths and suffixes are upgraded to orig/large."""
        web_url = "https://media01.symbiota.org/media/seinet/sernec/NCU/web/NCU001_web.jpg"
        optimized = optimize_herbarium_image_url(web_url)
        self.assertIn("/orig/", optimized)
        self.assertIn("_lg.jpg", optimized)

        tn_url = "https://media01.symbiota.org/media/seinet/sernec/NCU/tn/NCU001_tn.jpg"
        optimized_tn = optimize_herbarium_image_url(tn_url)
        self.assertIn("/orig/", optimized_tn)
        self.assertIn("_lg.jpg", optimized_tn)

    def test_optimize_iiif_urls(self):
        """Verify that IIIF URLs are rewritten to request max/full resolution."""
        iiif_url = "https://images.herbarium.org/iiif/2/NCU001/full/!1000,1000/0/default.jpg"
        optimized = optimize_herbarium_image_url(iiif_url)
        self.assertIn("/full/max/", optimized)

    def test_extract_high_res_image_url_media_ranking(self):
        """Verify that extract_high_res_image_url selects original/large image even when thumbnail is listed first."""
        media_list = [
            {
                "type": "StillImage",
                "format": "image/jpeg",
                "identifier": "https://media.symbiota.org/NCU/thumbnails/NCU001_tn.jpg",
            },
            {
                "type": "StillImage",
                "format": "image/jpeg",
                "identifier": "https://media.symbiota.org/NCU/original/NCU001_lg.jpg",
            },
        ]
        selected_url = extract_high_res_image_url(media_list)
        self.assertIn("_lg.jpg", selected_url)
        self.assertNotIn("_tn.jpg", selected_url)

    def test_validate_image_quality_high_vs_low_resolution(self):
        """Verify that validate_image_quality accepts high-res sheets and rejects low-res images."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            # 1. Low-res image (800 x 1000 = 0.8 MP)
            low_res_img = tmp_path / "low_res.jpg"
            img_small = Image.new("RGB", (800, 1000), color=(128, 128, 128))
            img_small.save(low_res_img, "JPEG")

            is_valid, metrics = validate_image_quality(low_res_img, min_megapixels=8.0, min_file_size_kb=1.0)
            self.assertFalse(is_valid)
            self.assertEqual(metrics["reason"], "low_resolution")
            self.assertLess(metrics["megapixels"], 8.0)

            # 2. High-res image (3000 x 4000 = 12.0 MP)
            high_res_img = tmp_path / "high_res.jpg"
            img_large = Image.new("RGB", (3000, 4000), color=(128, 128, 128))
            img_large.save(high_res_img, "JPEG")

            is_valid_high, metrics_high = validate_image_quality(high_res_img, min_megapixels=8.0, min_file_size_kb=1.0)
            self.assertTrue(is_valid_high)
            self.assertEqual(metrics_high["megapixels"], 12.0)


class TestAtomicTableExport(unittest.TestCase):
    """Test suite for atomic CSV table persistence."""

    def test_atomic_export_creates_file_with_dwc_columns(self):
        """Verify export_curated_table creates destination CSV with DwC headers."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = Path(tmpdir) / "curated_test.csv"
            sample_df = pd.DataFrame([
                {
                    "catalogNumber": "NCU001",
                    "institutionCode": "NCU",
                    "scientificName": "Packera dubia",
                    "species_raw": "Packera dubia",
                    "identifiedBy": "D.K. Trock",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "None",
                    "county": "Wake",
                    "stateProvince": "North Carolina",
                    "decimalLatitude": 35.8,
                    "decimalLongitude": -78.6,
                    "coordinateUncertainty": 500.0,
                    "year": 1998,
                    "month": 5,
                    "day": 12,
                    "eventDate": "1998-05-12",
                    "regional_group": "Piedmont_Granite_Flatrocks",
                    "image_path": "data/raw_vouchers/NCU001.jpg",
                }
            ])

            exported_path = export_curated_table(sample_df, out_file)
            self.assertTrue(exported_path.exists())
            read_df = pd.read_csv(exported_path)
            self.assertEqual(len(read_df), 1)
            self.assertEqual(read_df.iloc[0]["catalogNumber"], "NCU001")
            self.assertEqual(read_df.iloc[0]["year"], 1998)
            self.assertEqual(read_df.iloc[0]["determiner_tier"], "Tier_1_Gold")

            # Validate Darwin Core contract headers
            expected_contract_headers = [
                "catalogNumber",
                "scientificName",
                "decimalLatitude",
                "decimalLongitude",
                "eventDate",
                "identifiedBy",
                "determiner_tier",
                "image_path",
            ]
            for col in expected_contract_headers:
                self.assertIn(col, read_df.columns, f"Missing required contract header '{col}'")

            self.assertIn("species_raw", read_df.columns)
            for purged_col in ("latitude", "longitude"):
                self.assertNotIn(purged_col, read_df.columns)
            self.assertEqual(read_df.iloc[0]["scientificName"], "Packera dubia")
            self.assertEqual(read_df.iloc[0]["species_raw"], "Packera dubia")
            self.assertEqual(read_df.iloc[0]["identifiedBy"], "D.K. Trock")
            self.assertAlmostEqual(read_df.iloc[0]["decimalLatitude"], 35.8)
            self.assertAlmostEqual(read_df.iloc[0]["decimalLongitude"], -78.6)

class TestVoucherHarvesterCLI(unittest.TestCase):
    """Test suite for 01_voucher_harvester.py CLI parser and main entrypoint."""

    def test_cli_parser_defaults(self):
        """Verify default CLI arguments match PipelineConfig specifications."""
        import importlib
        mod_cli = importlib.import_module("scripts.data_prep.01_voucher_harvester")
        parser = mod_cli.build_cli_parser()
        args = parser.parse_args([])

        self.assertIsInstance(args.taxa, list)
        self.assertEqual(args.taxa, [])
        self.assertFalse(hasattr(args, "all_eastern_taxa"))
        with self.assertRaises(SystemExit):
            parser.parse_args(["--all-eastern-taxa"])
        self.assertEqual(args.max_records, 5000)
        self.assertEqual(args.min_megapixels, 8.0)
        self.assertEqual(args.min_file_size_kb, 500.0)
        self.assertTrue(args.exclude_western)
        self.assertFalse(args.download_images)

    def test_cli_parser_custom_args(self):
        """Verify custom CLI arguments are correctly parsed."""
        import importlib
        mod_cli = importlib.import_module("scripts.data_prep.01_voucher_harvester")
        parser = mod_cli.build_cli_parser()
        args = parser.parse_args([
            "--taxa", "Packera dubia", "Packera anonyma",
            "--max-records", "50",
            "--out-dir", "/tmp/custom_vouchers",
            "--output-csv", "/tmp/custom_curated.csv",
            "--download-images",
            "--check-sharpness",
            "--force",
        ])

        self.assertEqual(args.taxa, ["Packera dubia", "Packera anonyma"])
        self.assertEqual(args.max_records, 50)
        self.assertEqual(args.out_dir, "/tmp/custom_vouchers")
        self.assertEqual(args.output_csv, "/tmp/custom_curated.csv")
        self.assertTrue(args.download_images)
        self.assertTrue(args.check_sharpness)
        self.assertTrue(args.force)

    def test_cli_main_invocation(self):
        """Verify main() entrypoint constructs VoucherHarvester and triggers run()."""
        import importlib
        mod_cli = importlib.import_module("scripts.data_prep.01_voucher_harvester")

        with patch.object(sys, "argv", ["01_voucher_harvester.py", "--max-records", "5", "--no-exclude-western", "--force"]), \
             patch.object(mod_cli, "VoucherHarvester") as mock_harvester_cls:
            mock_instance = MagicMock()
            mock_harvester_cls.return_value = mock_instance

            mod_cli.main()

            mock_harvester_cls.assert_called_once()
            _, kwargs = mock_harvester_cls.call_args
            self.assertEqual(kwargs["max_records_per_taxon"], 5)
            self.assertFalse(kwargs["exclude_western"])
            self.assertTrue(kwargs["force"])
            mock_instance.run.assert_called_once_with(download_images=False)


class TestVoucherHarvesterResumptionLogic(unittest.TestCase):
    """Test suite for download resumption, atomic downloads, and CSV merging."""

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_download_single_image_resumption_skips_existing(self):
        """Verify that an existing non-empty voucher image is retained and skipped when force=False."""
        import asyncio
        from scripts.core.harvester import download_single_image

        dest_file = self.dir_path / "NCU000123.jpg"
        dest_file.write_bytes(b"\xff\xd8\xff" + b"A" * 2048)

        mock_session = MagicMock()
        mock_logger = MagicMock()
        semaphore = asyncio.Semaphore(1)

        success, was_skipped = asyncio.run(
            download_single_image(
                session=mock_session,
                image_url="https://example.com/sheet.jpg",
                destination_path=dest_file,
                semaphore=semaphore,
                force=False,
                logger=mock_logger,
            )
        )

        self.assertTrue(success)
        self.assertTrue(was_skipped)
        mock_session.get.assert_not_called()
        mock_logger.info.assert_called()
        self.assertIn("Skipping download for existing voucher", mock_logger.info.call_args[0][0])

    def test_download_single_image_force_redownloads(self):
        """Verify that when force=True, existing files are re-downloaded and replaced atomically."""
        import asyncio
        from scripts.core.harvester import download_single_image

        dest_file = self.dir_path / "NCU000456.jpg"
        dest_file.write_bytes(b"old_data")

        new_image_data = b"\xff\xd8\xff" + b"X" * 2048

        # Mock aiohttp response context manager
        mock_response = MagicMock()
        mock_response.status = 200
        mock_response.headers = {"Content-Type": "image/jpeg"}

        async def mock_read():
            return new_image_data
        mock_response.read = mock_read

        class MockGetContext:
            async def __aenter__(self):
                return mock_response
            async def __aexit__(self, exc_type, exc_val, exc_tb):
                pass

        mock_session = MagicMock()
        mock_session.get.return_value = MockGetContext()
        semaphore = asyncio.Semaphore(1)

        success, was_skipped = asyncio.run(
            download_single_image(
                session=mock_session,
                image_url="https://example.com/sheet.jpg",
                destination_path=dest_file,
                semaphore=semaphore,
                force=True,
            )
        )

        self.assertTrue(success)
        self.assertFalse(was_skipped)
        mock_session.get.assert_called_once()
        self.assertEqual(dest_file.read_bytes(), new_image_data)

    def test_download_all_voucher_images_counters(self):
        """Verify download_all_voucher_images aggregates progress counters [Processed: X | Skipped: Y | Total: Z]."""
        import asyncio
        from scripts.core.harvester import download_all_voucher_images

        img1 = self.dir_path / "NCU_EXISTING.jpg"
        img1.write_bytes(b"\xff\xd8\xff" + b"Z" * 1500)
        img2 = self.dir_path / "NCU_PENDING.jpg"

        records = [
            ("https://example.com/img1.jpg", img1),
            ("https://example.com/img2.jpg", img2),
        ]

        mock_logger = MagicMock()
        with patch("scripts.core.harvester.download_single_image", return_value=(True, False)):
            stats = asyncio.run(
                download_all_voucher_images(
                    records_to_download=records,
                    concurrency_limit=2,
                    force=False,
                    logger=mock_logger,
                )
            )

        self.assertEqual(stats["total"], 2)
        self.assertEqual(stats["skipped"], 1)
        self.assertEqual(stats["success"], 1)
        mock_logger.info.assert_called()

    def test_export_curated_table_merges_without_duplicates(self):
        """Verify export_curated_table merges new records with existing CSV without duplicate catalogNumber rows."""
        csv_path = self.dir_path / "curated_vouchers.csv"

        df_initial = pd.DataFrame([
            {"catalogNumber": "NCU001", "scientificName": "Packera dubia", "year": 2020},
            {"catalogNumber": "NCU002", "scientificName": "Packera dubia", "year": 2021},
        ])
        export_curated_table(df_initial, csv_path)
        self.assertEqual(len(pd.read_csv(csv_path)), 2)

        # Second harvest includes NCU002 (duplicate) and NCU003 (new)
        df_new = pd.DataFrame([
            {"catalogNumber": "NCU002", "scientificName": "Packera dubia", "year": 2022},
            {"catalogNumber": "NCU003", "scientificName": "Packera dubia", "year": 2023},
        ])
        export_curated_table(df_new, csv_path)

        df_merged = pd.read_csv(csv_path)
        self.assertEqual(len(df_merged), 3, "Merged table must contain exactly 3 unique vouchers.")
        self.assertListEqual(list(df_merged["catalogNumber"]), ["NCU001", "NCU002", "NCU003"])
        # Should keep last updated value for NCU002
        ncu002_year = df_merged[df_merged["catalogNumber"] == "NCU002"]["year"].iloc[0]
        self.assertEqual(ncu002_year, 2022)


class TestExsiccataeFingerprintingAndStratification(unittest.TestCase):
    """Test suite for exsiccatae collection event fingerprinting and multi-tier duplicate stratification."""

    def test_normalize_collector(self):
        """Verify collector name normalization: lowercase, strip initials, extract primary surname."""
        test_cases = [
            ("J. Brandon Fuller", "fuller"),
            ("J.B. Fuller", "fuller"),
            ("Fuller, J. Brandon", "fuller"),
            ("Fuller, J. B.", "fuller"),
            ("Fuller", "fuller"),
            ("A. Cronquist", "cronquist"),
            ("R. R. Kowal & J. B. Fuller", "kowal"),
            ("Robert R. Kowal with J. Brandon Fuller", "kowal"),
            ("Radford, Ahles & Bell", "radford"),
            ("J. B.", "jb"),
            ("Unknown", "unknown"),
            ("", "unknown"),
            (None, "unknown"),
        ]
        for raw_name, expected in test_cases:
            with self.subTest(raw_name=raw_name):
                self.assertEqual(normalize_collector(raw_name), expected)

    def test_normalize_record_number(self):
        """Verify collection number normalization: strip prefixes/suffixes, extract numeric component."""
        test_cases = [
            ("#1042", "1042"),
            ("1042", "1042"),
            ("No. 1042", "1042"),
            ("No. 1042b", "1042b"),
            ("1042-A", "1042a"),
            ("2022-014", "2022-014"),
            ("2022-015", "2022-015"),
            ("Coll #5821", "5821"),
            ("s.n.", ""),
            ("s. n.", ""),
            ("sn", ""),
            ("unnumbered", ""),
            ("", ""),
            (None, ""),
        ]
        for raw_num, expected in test_cases:
            with self.subTest(raw_num=raw_num):
                self.assertEqual(normalize_record_number(raw_num), expected)

    def test_normalize_event_date(self):
        """Verify eventDate normalization to standard YYYY-MM-DD format."""
        test_cases = [
            ({"event_date_raw": "1984-05-12"}, "1984-05-12"),
            ({"event_date_raw": "1984-05-12T00:00:00Z"}, "1984-05-12"),
            ({"event_date_raw": "1984/05/12"}, "1984-05-12"),
            ({"year": 1984, "month": 5, "day": 12}, "1984-05-12"),
            ({"year": 1984}, "1984-00-00"),
            ({"event_date_raw": "1984"}, "1984-00-00"),
            ({}, "unknown_date"),
        ]
        for kwargs, expected in test_cases:
            with self.subTest(kwargs=kwargs):
                self.assertEqual(normalize_event_date(**kwargs), expected)

    def test_generate_exsiccatae_key(self):
        """Verify construction of standardized signature string f'{norm_collector}_{norm_number}_{norm_date}'."""
        key = generate_exsiccatae_key(
            collector="J. Brandon Fuller",
            number="#1042",
            date="1984-05-12",
        )
        self.assertEqual(key, "fuller_1042_1984-05-12")

    def test_haversine_distance_meters(self):
        """Verify haversine distance calculates accurate geodesic distances."""
        # Exact same point
        self.assertAlmostEqual(haversine_distance_meters(35.9123, -79.0512, 35.9123, -79.0512), 0.0, places=1)
        # Shift latitude by ~0.0003 deg (~33 m < 50 m)
        dist_near = haversine_distance_meters(35.9123, -79.0512, 35.9126, -79.0512)
        self.assertTrue(25.0 < dist_near < 40.0)
        # Shift latitude by ~0.001 deg (~111 m > 50 m)
        dist_far = haversine_distance_meters(35.9123, -79.0512, 35.9133, -79.0512)
        self.assertTrue(100.0 < dist_far < 125.0)

    def test_stratify_duplicates_numbered_exsiccatae(self):
        """Verify grouping of duplicate sheets sharing collector, collection number, and date."""
        df_exsiccatae = pd.DataFrame([
            {
                "catalogNumber": "NCU001",
                "institutionCode": "NCU",
                "recordedBy": "J. Brandon Fuller",
                "recordNumber": "1042",
                "eventDate": "1984-05-12",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
            {
                "catalogNumber": "GA002",
                "institutionCode": "GA",
                "recordedBy": "J.B. Fuller",
                "recordNumber": "#1042",
                "eventDate": "1984-05-12",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "None",
            },
            {
                "catalogNumber": "US003",
                "institutionCode": "US",
                "recordedBy": "Fuller, J. Brandon",
                "recordNumber": "No. 1042",
                "eventDate": "1984-05-12",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "None",
            },
        ])
        strat_df = stratify_duplicates(df_exsiccatae)

        # All 3 records must share identical exsiccatae_key
        self.assertEqual(len(strat_df["exsiccatae_key"].unique()), 1)
        self.assertEqual(strat_df["exsiccatae_key"].iloc[0], "fuller_1042_1984-05-12")

        # duplicate_count must equal 3 for all records
        self.assertTrue((strat_df["duplicate_count"] == 3).all())

        # Exactly 1 primary duplicate, remaining are False
        self.assertEqual(strat_df["is_primary_duplicate"].sum(), 1)
        # Priority rule: NCU beats GA and US (Tier 2 vs Tier 3, and NCU is Home Herbarium)
        ncu_row = strat_df[strat_df["catalogNumber"] == "NCU001"].iloc[0]
        self.assertTrue(ncu_row["is_primary_duplicate"])
        self.assertFalse(strat_df[strat_df["catalogNumber"] == "GA002"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(strat_df[strat_df["catalogNumber"] == "US003"]["is_primary_duplicate"].iloc[0])

    def test_stratify_duplicates_coordinate_proximity_fallback(self):
        """Verify unnumbered collections (< 50m) on identical dates are clustered as duplicates."""
        df_proximity = pd.DataFrame([
            {
                "catalogNumber": "SHEET_A",
                "institutionCode": "NCU",
                "recordedBy": "J. Brandon Fuller",
                "recordNumber": "s.n.",
                "eventDate": "1984-05-12",
                "decimalLatitude": 35.91230,
                "decimalLongitude": -79.05120,
                "determiner_tier": "Tier_2_Silver",
            },
            {
                "catalogNumber": "SHEET_B",
                "institutionCode": "GA",
                "recordedBy": "Fuller",
                "recordNumber": "",
                "eventDate": "1984-05-12",
                "decimalLatitude": 35.91245,  # ~17 m away (< 50 m)
                "decimalLongitude": -79.05125,
                "determiner_tier": "Tier_3_Bronze",
            },
            {
                "catalogNumber": "SHEET_C",
                "institutionCode": "MO",
                "recordedBy": "Fuller",
                "recordNumber": "s.n.",
                "eventDate": "1984-05-12",
                "decimalLatitude": 36.50000,  # ~65 km away (> 50 m)
                "decimalLongitude": -79.05120,
                "determiner_tier": "Tier_3_Bronze",
            },
        ])
        strat_df = stratify_duplicates(df_proximity, max_distance_meters=50.0)

        # SHEET_A and SHEET_B should share the same exsiccatae cluster key
        key_a = strat_df[strat_df["catalogNumber"] == "SHEET_A"]["exsiccatae_key"].iloc[0]
        key_b = strat_df[strat_df["catalogNumber"] == "SHEET_B"]["exsiccatae_key"].iloc[0]
        key_c = strat_df[strat_df["catalogNumber"] == "SHEET_C"]["exsiccatae_key"].iloc[0]

        self.assertEqual(key_a, key_b, "SHEET_A and SHEET_B must be clustered together within 50m proximity.")
        self.assertNotEqual(key_a, key_c, "SHEET_C (>50m away) must not be clustered with SHEET_A/B.")

        # Duplicate counts
        self.assertEqual(strat_df[strat_df["catalogNumber"] == "SHEET_A"]["duplicate_count"].iloc[0], 2)
        self.assertEqual(strat_df[strat_df["catalogNumber"] == "SHEET_B"]["duplicate_count"].iloc[0], 2)
        self.assertEqual(strat_df[strat_df["catalogNumber"] == "SHEET_C"]["duplicate_count"].iloc[0], 1)

        # Primary duplicate flag: NCU (SHEET_A) wins over GA (SHEET_B)
        self.assertTrue(strat_df[strat_df["catalogNumber"] == "SHEET_A"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(strat_df[strat_df["catalogNumber"] == "SHEET_B"]["is_primary_duplicate"].iloc[0])
        # SHEET_C is singleton, so it is primary for its own group
        self.assertTrue(strat_df[strat_df["catalogNumber"] == "SHEET_C"]["is_primary_duplicate"].iloc[0])

    def test_institutional_priority_rule(self):
        """
        Verify institutional hierarchy:
        (1) Verified Type Status -> (2) Home Herbarium (NCU) -> (3) Tier 1 Specialist -> (4) Image resolution.
        """
        # Scenario 1: Verified Type status (GA) beats NCU (non-type)
        df_type_test = pd.DataFrame([
            {
                "catalogNumber": "NCU_NORM",
                "institutionCode": "NCU",
                "recordedBy": "Fuller",
                "recordNumber": "100",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "None",
            },
            {
                "catalogNumber": "GA_ISOTYPE",
                "institutionCode": "GA",
                "recordedBy": "Fuller",
                "recordNumber": "100",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "Isotype",
            },
        ])
        res1 = stratify_duplicates(df_type_test)
        self.assertTrue(res1[res1["catalogNumber"] == "GA_ISOTYPE"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(res1[res1["catalogNumber"] == "NCU_NORM"]["is_primary_duplicate"].iloc[0])

        # Scenario 2: Home Herbarium NCU beats non-NCU when neither is a type
        df_ncu_test = pd.DataFrame([
            {
                "catalogNumber": "NCU_SHEET",
                "institutionCode": "NCU",
                "recordedBy": "Fuller",
                "recordNumber": "200",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
            {
                "catalogNumber": "MO_SHEET",
                "institutionCode": "MO",
                "recordedBy": "Fuller",
                "recordNumber": "200",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
        ])
        res2 = stratify_duplicates(df_ncu_test)
        self.assertTrue(res2[res2["catalogNumber"] == "NCU_SHEET"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(res2[res2["catalogNumber"] == "MO_SHEET"]["is_primary_duplicate"].iloc[0])

        # Scenario 3: Tier 1 Specialist determination beats Tier 2/3 when neither is NCU and neither is Type
        df_tier1_test = pd.DataFrame([
            {
                "catalogNumber": "GA_TIER1",
                "institutionCode": "GA",
                "recordedBy": "Fuller",
                "recordNumber": "300",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "None",
            },
            {
                "catalogNumber": "MO_TIER3",
                "institutionCode": "MO",
                "recordedBy": "Fuller",
                "recordNumber": "300",
                "eventDate": "1990-06-01",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "None",
            },
        ])
        res3 = stratify_duplicates(df_tier1_test)
        self.assertTrue(res3[res3["catalogNumber"] == "GA_TIER1"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(res3[res3["catalogNumber"] == "MO_TIER3"]["is_primary_duplicate"].iloc[0])


class TestSpatialBalancingAndThinning(unittest.TestCase):
    """Test suite for Milestone 1/4 spatial grid thinning, state/region quotas, and inverse-density weights."""

    def test_spatially_balance_and_thin_vouchers_nc_cluster_vs_ar_spread(self):
        """
        Verify that 200 NC vouchers clustered in 2 counties and 15 AR vouchers across 10 counties
        are thinned by spatial_thin_km=15.0, max_per_state=25 to at most 25 NC vouchers while
        retaining all Tier_1_Gold vouchers and computing valid spatial_balance_weight in [0.20, 3.00].
        """
        rng = np.random.default_rng(42)
        records = []

        # 200 NC vouchers clustered across 2 counties (Wake & Franklin), including 6 Tier_1_Gold anchors
        # and 10 secondary duplicate sheets
        gold_catalog_numbers = set()
        for i in range(200):
            is_gold = i < 6
            cat = f"NC_VOUCHER_{i:03d}"
            if is_gold:
                gold_catalog_numbers.add(cat)
            county = "Wake" if i % 2 == 0 else "Franklin"
            # Spread across ~30 grid cells in Wake/Franklin so grid thinning leaves >25 before max_per_state=25 caps at 25
            base_lat = 35.78 + (i % 10) * 0.16
            base_lon = -78.64 + ((i // 10) % 6) * 0.18
            records.append({
                "catalogNumber": cat,
                "institutionCode": "NCU",
                "scientificName": "Packera anonyma",
                "species_raw": "Packera anonyma",
                "determiner_tier": "Tier_1_Gold" if is_gold else ("Tier_2_Silver" if i % 5 == 0 else "Tier_3_Bronze"),
                "type_status": "Holotype" if i == 0 else "None",
                "county": county,
                "stateProvince": "NC",
                "decimalLatitude": float(base_lat + rng.normal(0.0, 0.005)),
                "decimalLongitude": float(base_lon + rng.normal(0.0, 0.005)),
                "coordinateUncertainty": 50.0 if is_gold else 500.0,
                "regional_group": "Piedmont_Granite_Flatrocks",
                "is_primary_duplicate": False if i >= 190 else True,
            })

        # 15 AR vouchers spread across 10 counties, separated by >25 km each
        for j in range(15):
            cat = f"AR_VOUCHER_{j:03d}"
            if j == 0:
                gold_catalog_numbers.add(cat)
            records.append({
                "catalogNumber": cat,
                "institutionCode": "UARK",
                "scientificName": "Packera anonyma",
                "species_raw": "Packera anonyma",
                "determiner_tier": "Tier_1_Gold" if j == 0 else "Tier_2_Silver",
                "type_status": "None",
                "county": f"AR_County_{j % 10}",
                "stateProvince": "AR",
                "decimalLatitude": float(34.20 + (j % 5) * 0.35),
                "decimalLongitude": float(-94.10 + (j // 5) * 0.45),
                "coordinateUncertainty": 100.0,
                "regional_group": "Interior_Prairie_Midwest",
                "is_primary_duplicate": True,
            })

        df_in = pd.DataFrame(records)
        manifest_df = pd.DataFrame([
            {"catalogNumber": "NC_VOUCHER_010", "leaf_id": 1, "efa_qc_pass": True},
            {"catalogNumber": "NC_VOUCHER_010", "leaf_id": 2, "efa_qc_pass": True},
        ])

        thinned = spatially_balance_and_thin_vouchers(
            df_in,
            spatial_thin_km=15.0,
            max_per_state=25,
            primary_duplicates_only=True,
            preserve_all_gold=True,
            manifest_df=manifest_df,
        )

        # 1. Secondary duplicates must be excluded
        self.assertTrue(thinned["is_primary_duplicate"].all())
        for i in range(190, 200):
            self.assertNotIn(f"NC_VOUCHER_{i:03d}", set(thinned["catalogNumber"]))

        # 2. NC vouchers capped at max_per_state=25
        nc_retained = thinned[thinned["stateProvince"] == "NC"]
        ar_retained = thinned[thinned["stateProvince"] == "AR"]
        self.assertLessEqual(len(nc_retained), 25)
        self.assertEqual(len(nc_retained), 25)
        self.assertEqual(len(ar_retained), 15)

        # 3. Both NC counties represented evenly and 100% of Tier_1_Gold vouchers retained
        self.assertEqual(set(nc_retained["county"].unique()), {"Wake", "Franklin"})
        self.assertTrue(gold_catalog_numbers.issubset(set(thinned["catalogNumber"])))

        # 4. Continuous spatial_balance_weight bounded in [0.20, 3.00] and normalized around 1.0
        self.assertIn("spatial_balance_weight", thinned.columns)
        weights = thinned["spatial_balance_weight"].to_numpy(dtype=float)
        self.assertTrue(np.all(np.isfinite(weights)))
        self.assertTrue(np.all((weights >= 0.20) & (weights <= 3.00)))
        self.assertAlmostEqual(float(np.mean(weights)), 1.0, delta=0.15)

        # Under-represented, spatially dispersed AR vouchers should have higher mean weight than dense NC vouchers
        mean_w_ar = float(ar_retained["spatial_balance_weight"].mean())
        mean_w_nc = float(nc_retained["spatial_balance_weight"].mean())
        self.assertGreater(mean_w_ar, mean_w_nc)

    def test_max_per_region_and_cli_arguments(self):
        """Verify max_per_region quota balancing and CLI argument parsing in 01_voucher_harvester.py."""
        import importlib
        cli_mod = importlib.import_module("scripts.data_prep.01_voucher_harvester")
        parser = cli_mod.build_cli_parser()
        args = parser.parse_args([
            "--spatial-thin-km", "15.0",
            "--max-per-state", "50",
            "--max-per-region", "150",
        ])
        self.assertAlmostEqual(args.spatial_thin_km, 15.0)
        self.assertEqual(args.max_per_state, 50)
        self.assertEqual(args.max_per_region, 150)


class TestEPALevel3EcoregionAssignment(unittest.TestCase):
    """Test suite for official US EPA Level III Ecoregion polygon assignment & county-containment fallback."""

    def test_point_in_polygon_epa_level3_ecoregion_assignment(self):
        """Verify georeferenced points map strictly by polygon intersection across Atlantic, Gulf, TX/AR, and Interior ecoregions."""
        cases = [
            # (label, lat, lon, expected_ecoregion)
            ("Chapel Hill, NC (Piedmont)", 35.9132, -79.0558, "Piedmont"),
            ("Wilmington, NC (Middle Atlantic Coastal Plain)", 34.2257, -77.9447, "Middle_Atlantic_Coastal_Plain"),
            ("Dublin, GA (Southeastern Plains)", 32.5404, -82.9037, "Southeastern_Plains"),
            ("Gainesville, FL (Southern Coastal Plain)", 29.6516, -82.3248, "Southern_Coastal_Plain"),
            ("Nacogdoches, TX (South Central Plains)", 31.6035, -94.6555, "South_Central_Plains"),
            ("Camden, AR (South Central Plains)", 33.5846, -92.8343, "South_Central_Plains"),
            ("College Station, TX (East Central Texas Plains)", 30.6280, -96.3344, "East_Central_Texas_Plains"),
            ("Greenville, MS (Mississippi Alluvial Plain)", 33.4101, -91.0618, "Mississippi_Alluvial_Plain"),
            ("Asheville, NC (Blue Ridge)", 35.5951, -82.5515, "Blue_Ridge"),
            ("Knoxville, TN (Ridge and Valley)", 35.9606, -83.9207, "Ridge_and_Valley"),
            ("Sewanee, TN (Southwestern Appalachians)", 35.2031, -85.9211, "Southwestern_Appalachians"),
            ("Nashville, TN (Interior Plateau)", 36.1627, -86.7816, "Interior_Plateau"),
            ("Harrison, AR (Ozark Highlands)", 36.2298, -93.1077, "Ozark_Highlands"),
            ("Hot Springs, AR (Ouachita Mountains)", 34.5037, -93.0552, "Ouachita_Mountains"),
        ]
        for desc, lat, lon, expected in cases:
            with self.subTest(desc=desc):
                # Also pass misleading habitat/state text to prove text heuristics and state lumping are ignored
                got = infer_regional_group(
                    lat=lat,
                    lon=lon,
                    state_province="Florida",
                    habitat="sandhills granite flatrock",
                    locality="Blue Ridge mountains",
                )
                self.assertEqual(got, expected, f"Failed for {desc}: got {got!r}, expected {expected!r}")
                self.assertNotEqual(got, "Other_US")

        # Offshore / Null Island points outside all US EPA Level III polygons must return None (NA), never 'Other_US'
        self.assertIsNone(infer_regional_group(0.0, 0.0, state_province="Tennessee", county="Cheatham"))
        self.assertIsNone(infer_regional_group(35.16, -22.28, state_province="South Carolina", county="Greenville"))

    def test_county_containment_fallback_and_straddling_county_na(self):
        """Verify ungeoreferenced single-ecoregion counties resolve via centroid while straddling/missing counties return None (NA)."""
        # Single-ecoregion (>= 95% share) counties
        self.assertEqual(
            infer_regional_group(None, None, state_province="North Carolina", county="Orange County"),
            "Piedmont",
        )
        self.assertEqual(
            infer_regional_group(np.nan, np.nan, state_province="NC", county="New Hanover"),
            "Middle_Atlantic_Coastal_Plain",
        )
        self.assertEqual(
            infer_regional_group(None, None, state_province="TX", county="Nacogdoches Co."),
            "South_Central_Plains",
        )
        self.assertEqual(
            infer_regional_group(None, None, state_province="Arkansas", county="Baxter County"),
            "Ozark_Highlands",
        )
        self.assertEqual(
            infer_regional_group(None, None, state_province="AR", county="Garland"),
            "Ouachita_Mountains",
        )

        # Multi-ecoregion straddling counties (Fall Line & Blue Ridge Escarpment) must return None (NA)
        self.assertIsNone(infer_regional_group(None, None, state_province="North Carolina", county="Wake County"))
        self.assertIsNone(infer_regional_group(None, None, state_province="South Carolina", county="Richland"))
        self.assertIsNone(infer_regional_group(None, None, state_province="Georgia", county="Bibb County"))
        self.assertIsNone(infer_regional_group(None, None, state_province="North Carolina", county="Burke"))

        # Missing, ambiguous, or unmatched counties must return None (NA), never 'Other_US'
        self.assertIsNone(infer_regional_group(None, None, state_province="North Carolina", county=None))
        self.assertIsNone(infer_regional_group(None, None, state_province="Texas", county=""))
        self.assertIsNone(infer_regional_group(None, None, state_province="Tennessee", county="Sevier or Blount"))
        self.assertIsNone(infer_regional_group(None, None, state_province="Ontario", county="Bruce"))

    def test_deterministic_offline_fallback_and_curated_vouchers_integrity(self):
        """Verify deterministic offline polygon mode and zero 'Other_US' in curated_vouchers.csv."""
        with tempfile.TemporaryDirectory() as tmpdir:
            missing_eco = Path(tmpdir) / "nonexistent_eco.geojson"
            missing_cnty = Path(tmpdir) / "nonexistent_cnty.geojson"
            self.assertEqual(
                infer_regional_group(35.91, -79.05, ecoregion_path=missing_eco, county_path=missing_cnty),
                "Piedmont",
            )
            self.assertEqual(
                infer_regional_group(31.60, -94.65, ecoregion_path=missing_eco, county_path=missing_cnty),
                "South_Central_Plains",
            )
            self.assertEqual(
                infer_regional_group(None, None, state_province="NC", county="Orange", ecoregion_path=missing_eco, county_path=missing_cnty),
                "Piedmont",
            )
            self.assertIsNone(
                infer_regional_group(None, None, state_province="NC", county="Wake", ecoregion_path=missing_eco, county_path=missing_cnty),
            )
            # Offline fallback must never write synthetic files to disk
            self.assertFalse(missing_eco.exists())
            self.assertFalse(missing_cnty.exists())

        prod_csv = Path(__file__).resolve().parents[2] / "data" / "tables" / "curated_vouchers.csv"
        if prod_csv.exists():
            df = pd.read_csv(prod_csv, low_memory=False)
            vals = set(df["regional_group"].dropna().astype(str))
            self.assertNotIn("Other_US", vals)
            self.assertNotIn("Coastal_Plain_Sandhills", vals)
            self.assertNotIn("Piedmont_Granite_Flatrocks", vals)
            self.assertIn("Piedmont", vals)
            self.assertIn("Middle_Atlantic_Coastal_Plain", vals)
            self.assertIn("South_Central_Plains", vals)


class TestDynamicAuthorityEngineAndSynonymExpansion(unittest.TestCase):
    """Test suite for the 3-Layer Dynamic Authority Engine, GBIF synonym expansion, and geographic scopes."""

    def test_extract_protologue_authors_basionym_and_combination(self):
        """Verify Layer 2 extracts basionym/combination authors and IPNI expansions from scientificName strings."""
        from scripts.core.harvester import extract_protologue_authors

        dubia_authors = extract_protologue_authors(["Packera dubia (Spreng.) Trock & Mabb."])
        self.assertTrue(
            {"spreng", "sprengel", "trock", "mabb", "mabberley"}.issubset(dubia_authors),
            f"Expected Spreng./Trock/Mabb. tokens in {dubia_authors}",
        )

        paupercula_authors = extract_protologue_authors(["Packera paupercula (Michx.) Á.Löve & D.Löve"])
        self.assertTrue(
            {"michx", "michaux", "love", "löve"}.issubset(paupercula_authors),
            f"Expected Michx./Löve tokens in {paupercula_authors}",
        )

        western_authors = extract_protologue_authors([
            "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve",
            "Packera aurea var. gracilis (Pursh) Wight",
        ])
        self.assertTrue({"greene", "weber", "love", "pursh", "wight"}.issubset(western_authors))
        self.assertNotIn("gracilis", western_authors)
        self.assertNotIn("var", western_authors)

        # Verify Layer 2 promotes nomenclatural author annotations to Tier_1_Gold
        tier, _, _ = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw="D. J. Mabberley",
            recorded_by_raw="Field Collector",
            history_raw=None,
            institution_code_raw="XYZ_LOCAL",
            locality_raw="Sandy slope",
            habitat_raw=None,
            dynamic_specialists=dubia_authors,
        )
        self.assertEqual(tier, "Tier_1_Gold")

    def test_layer1_expanded_north_american_monographers(self):
        """Verify expanded genus-wide Packera/Senecioneae monographers classify as Tier_1_Gold."""
        monographers = [
            "C. C. Freeman",
            "A. Cronquist",
            "J. M. Greenman",
            "W. A. Weber",
            "Á. Löve",
            "R. Spellenberg",
            "R. W. Pippen",
            "L. J. Uttal",
            "M. L. Fernald",
            "P. A. Rydberg",
            "A. A. Heller",
            "R. Kral",
            "R. K. Godfrey",
            "B. A. Sorrie",
            "G. L. Nesom",
            "J. L. Strother",
            "B. L. Turner",
            "J. L. Golden",
        ]
        for det_name in monographers:
            with self.subTest(determiner=det_name):
                tier, _, _ = parse_determiner_tier(
                    type_status_raw=None,
                    identified_by_raw=det_name,
                    recorded_by_raw="Collector",
                    history_raw=None,
                    institution_code_raw="SMALL_HERB",
                    locality_raw=" Meadow",
                    habitat_raw=None,
                )
                self.assertEqual(tier, "Tier_1_Gold", f"{det_name} should map to Tier_1_Gold")

    def test_discover_corpus_authorities_promotes_multi_institution_specialist_not_local_collector(self):
        """
        Verify Layer 3 promotes a determiner with >= 5 annotations across >= 3 institutions
        (or type annotation) to Tier_1_Gold while keeping single-institution collectors in Silver/Bronze.
        """
        from scripts.core.harvester import apply_corpus_authority_tiers, discover_corpus_authorities

        rows = []
        # 1. Loan specialist: 5 sheets of Packera streptanthifolia across 3 institutions (COLO, RM, UTC)
        loan_insts = ["COLO", "RM", "UTC", "COLO", "RM"]
        for idx, inst in enumerate(loan_insts):
            rows.append({
                "catalogNumber": f"LOAN_{idx:02d}",
                "institutionCode": inst,
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "Dr. Elena V. Monographia, 2019",
                "recordedBy": f"Collector {idx}",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Subalpine rocky meadow at 2800 m elevation",
                "habitat": "Granitic talus",
            })

        # 2. Single-institution local collector/determiner: 6 sheets all at NCU (only 1 institution)
        for idx in range(6):
            rows.append({
                "catalogNumber": f"LOCAL_NCU_{idx:02d}",
                "institutionCode": "NCU",
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "M. Localcollector",
                "recordedBy": "M. Localcollector",
                "type_status": "None",
                "determiner_tier": "Tier_2_Silver",
                "locality": "Roadside clearing 4 miles west of town center",
                "habitat": "Open woods",
            })

        # 3. Single-institution collector at unknown herbarium without rich locality -> remains Bronze
        for idx in range(6):
            rows.append({
                "catalogNumber": f"LOCAL_BRONZE_{idx:02d}",
                "institutionCode": "TINY_COLL",
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "P. Singleherb",
                "recordedBy": "P. Singleherb",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Woods",
                "habitat": "",
            })

        # 4. Generic non-person strings across 4 institutions must NEVER be promoted
        for idx, inst in enumerate(["NCU", "MO", "NY", "US", "GH"]):
            rows.append({
                "catalogNumber": f"GENERIC_{idx:02d}",
                "institutionCode": inst,
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": "Herbarium Staff",
                "recordedBy": "Unknown",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Detailed locality over ten characters long",
                "habitat": "Rocky slope",
            })

        # 5. Dynamic regional research herbarium ('REG_UNIV') with 12 curated sheets & >=70% determiners
        for idx in range(12):
            rows.append({
                "catalogNumber": f"REG_UNIV_{idx:02d}",
                "institutionCode": "REG_UNIV",
                "scientificName": "Packera streptanthifolia",
                "species_raw": "Packera streptanthifolia",
                "identifiedBy": f"Curator Person{idx}",
                "recordedBy": "Regional Collector",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Limestone bluff 6 km north of river crossing",
                "habitat": "Calcareous outcrop",
            })

        df = pd.DataFrame(rows)
        dyn_specs, dyn_herbs = discover_corpus_authorities(df)

        self.assertIn("monographia", dyn_specs)
        self.assertNotIn("localcollector", dyn_specs)
        self.assertNotIn("singleherb", dyn_specs)
        self.assertNotIn("staff", dyn_specs)
        self.assertNotIn("herbarium", dyn_specs)
        self.assertIn("REG_UNIV", dyn_herbs)

        df_scored = apply_corpus_authority_tiers(df, dynamic_specialists=dyn_specs, dynamic_herbaria=dyn_herbs)
        tier_by_cat = dict(zip(df_scored["catalogNumber"], df_scored["determiner_tier"]))

        self.assertEqual(tier_by_cat["LOAN_00"], "Tier_1_Gold")
        self.assertEqual(tier_by_cat["LOAN_04"], "Tier_1_Gold")
        self.assertEqual(tier_by_cat["LOCAL_NCU_00"], "Tier_2_Silver")
        self.assertEqual(tier_by_cat["LOCAL_BRONZE_00"], "Tier_3_Bronze")
        self.assertEqual(tier_by_cat["GENERIC_00"], "Tier_3_Bronze")
        self.assertEqual(tier_by_cat["REG_UNIV_00"], "Tier_2_Silver")

    def test_expand_taxon_synonyms_gbif_and_dynamic_standardization(self):
        """Verify dynamic GBIF Backbone synonym expansion (Packera/Senecio/Cineraria) and reverse standardization."""
        from scripts.core.harvester import (
            _GBIF_SYNONYM_CACHE,
            expand_taxon_synonyms_gbif,
            standardize_packera_taxon,
        )

        mock_species = MagicMock()
        mock_species.name_backbone.return_value = {
            "usageKey": 999123,
            "acceptedUsageKey": 999123,
            "canonicalName": "Packera neomexicana",
            "scientificName": "Packera neomexicana (A.Gray) W.A.Weber & Á.Löve",
        }
        mock_species.name_usage.return_value = {
            "results": [
                {
                    "canonicalName": "Senecio neomexicanus",
                    "scientificName": "Senecio neomexicanus A.Gray",
                },
                {
                    "canonicalName": "Cineraria neomexicana",
                    "scientificName": "Cineraria neomexicana Hort.",
                },
            ]
        }

        _GBIF_SYNONYM_CACHE.pop("packera neomexicana", None)
        with patch("scripts.core.harvester.species", mock_species):
            expanded = expand_taxon_synonyms_gbif("Packera neomexicana")

        self.assertIn("Packera neomexicana", expanded)
        self.assertIn("Senecio neomexicanus", expanded)
        self.assertIn("Cineraria neomexicana", expanded)

        # Harvested synonyms must map back to canonical Packera <epithet>
        self.assertEqual(standardize_packera_taxon("Senecio neomexicanus A.Gray"), "Packera neomexicana")
        self.assertEqual(standardize_packera_taxon("Cineraria neomexicana"), "Packera neomexicana")
        self.assertEqual(standardize_packera_taxon("Packera cana (Hook.) W.A.Weber & Á.Löve"), "Packera cana")

        # Unseen historical Senecio synonym resolves dynamically via GBIF Backbone lookup
        mock_strep = MagicMock()
        mock_strep.name_backbone.return_value = {
            "usageKey": 888456,
            "acceptedUsageKey": 888000,
            "canonicalName": "Senecio streptanthifolius",
            "species": "Packera streptanthifolia",
            "accepted": "Packera streptanthifolia (Greene) W.A.Weber & Á.Löve",
        }
        mock_strep.name_usage.return_value = {"results": []}
        _GBIF_SYNONYM_CACHE.pop("senecio streptanthifolius", None)
        with patch("scripts.core.harvester.species", mock_strep):
            self.assertEqual(
                standardize_packera_taxon("Senecio streptanthifolius Greene"),
                "Packera streptanthifolia",
            )

    def test_geographic_scope_north_america_and_include_western_retain_western_vouchers(self):
        """Verify north_america / include_western retain < -106.65 western records while eastern_na excludes them."""
        mock_occ_records = [
            {
                "key": 5001,
                "scientificName": "Packera cana (Hook.) W.A.Weber & Á.Löve",
                "stateProvince": "Colorado",
                "decimalLatitude": 39.74,
                "decimalLongitude": -108.50,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2019, "month": 6, "day": 18,
                "media": [{"identifier": "https://example.com/colo_cana.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "COLO_WEST_01",
                "institutionCode": "COLO",
                "identifiedBy": "W. A. Weber",
            },
            {
                "key": 5002,
                "scientificName": "Packera breweri (Burtt Davy) W.A.Weber & Á.Löve",
                "stateProvince": "California",
                "decimalLatitude": 36.80,
                "decimalLongitude": -120.40,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2020, "month": 5, "day": 12,
                "media": [{"identifier": "https://example.com/uc_breweri.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "UC_WEST_02",
                "institutionCode": "UC",
                "identifiedBy": "J. L. Strother",
            },
            {
                "key": 5003,
                "scientificName": "Packera paupercula (Michx.) Á.Löve & D.Löve",
                "stateProvince": "North Carolina",
                "decimalLatitude": 35.91,
                "decimalLongitude": -79.05,
                "coordinateUncertaintyInMeters": 100.0,
                "year": 2021, "month": 5, "day": 3,
                "media": [{"identifier": "https://example.com/ncu_east.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "NCU_EAST_03",
                "institutionCode": "NCU",
                "identifiedBy": "A. S. Weakley",
            },
        ]

        def _fake_occ_search(**kwargs):
            country = kwargs.get("country")
            if country in (None, "US"):
                return {"results": mock_occ_records, "count": len(mock_occ_records)}
            return {"results": [], "count": 0}

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_search):
            # 1. Default eastern_na excludes Colorado & California records (< -106.65)
            df_east = harvest_taxa_occurrences(
                taxa_list=["Packera cana"],
                max_records_per_taxon=10,
            )
            self.assertEqual(set(df_east["catalogNumber"]), {"NCU_EAST_03"})

            # 2. geographic_scope='north_america' retains Colorado & California records and queries US/CA/MX
            harvester_na = VoucherHarvester(
                taxa=["Packera cana"],
                max_records_per_taxon=10,
                geographic_scope="north_america",
            )
            self.assertFalse(harvester_na.exclude_western)
            self.assertEqual(harvester_na.countries, ["US", "CA", "MX"])
            df_na = harvester_na.harvest()
            self.assertEqual(
                set(df_na["catalogNumber"]),
                {"COLO_WEST_01", "UC_WEST_02", "NCU_EAST_03"},
            )

            # 3. include_western=True also disables western exclusion
            df_inc_west = harvest_taxa_occurrences(
                taxa_list=["Packera cana"],
                max_records_per_taxon=10,
                include_western=True,
            )
            self.assertEqual(
                set(df_inc_west["catalogNumber"]),
                {"COLO_WEST_01", "UC_WEST_02", "NCU_EAST_03"},
            )

    def test_select_cascading_outgroup_anchors_and_home_herbarium(self):
        """Verify preflight outgroup anchor cascade (Gold -> Silver -> Bronze) and configurable home_herbarium."""
        from scripts.core.harvester import select_cascading_outgroup_anchors

        df_pool = pd.DataFrame([
            # Outgroup A has NO Gold vouchers, only 2 Silver and 2 Bronze
            {
                "catalogNumber": "OG_A_BRONZE_1",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_3_Bronze",
                "image_path": "data/raw_vouchers/OG_A_BRONZE_1.jpg",
                "is_primary_duplicate": True,
            },
            {
                "catalogNumber": "OG_A_SILVER_1",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_2_Silver",
                "image_path": "data/raw_vouchers/OG_A_SILVER_1.jpg",
                "is_primary_duplicate": True,
            },
            {
                "catalogNumber": "OG_A_SILVER_2",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_2_Silver",
                "image_path": "data/raw_vouchers/OG_A_SILVER_2.jpg",
                "is_primary_duplicate": True,
            },
            {
                "catalogNumber": "OG_A_BRONZE_2",
                "species_standardized": "Packera fendleri",
                "determiner_tier": "Tier_3_Bronze",
                "image_path": "data/raw_vouchers/OG_A_BRONZE_2.jpg",
                "is_primary_duplicate": True,
            },
        ])
        out_dfs, counts, tiers = select_cascading_outgroup_anchors(
            df_pool,
            outgroup_taxa=["Packera fendleri"],
            max_anchors=3,
        )
        self.assertEqual(counts["Packera fendleri"], 3)
        selected_cats = out_dfs[0]["catalogNumber"].tolist()
        # Both Silver sheets must rank ahead of Bronze fallback
        self.assertEqual(selected_cats[:2], ["OG_A_SILVER_1", "OG_A_SILVER_2"])
        self.assertEqual(tiers["Packera fendleri"]["Tier_2_Silver"], 2)
        self.assertEqual(tiers["Packera fendleri"]["Tier_3_Bronze"], 1)

        # Verify configurable home_herbarium override in stratify_duplicates
        df_dups = pd.DataFrame([
            {
                "catalogNumber": "NCU_DUP",
                "institutionCode": "NCU",
                "recordedBy": "W. A. Weber",
                "recordNumber": "4401",
                "eventDate": "1988-07-10",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
            {
                "catalogNumber": "COLO_DUP",
                "institutionCode": "COLO",
                "recordedBy": "W. A. Weber",
                "recordNumber": "4401",
                "eventDate": "1988-07-10",
                "determiner_tier": "Tier_2_Silver",
                "type_status": "None",
            },
        ])
        strat_colo = stratify_duplicates(df_dups, home_herbarium="COLO")
        self.assertTrue(strat_colo[strat_colo["catalogNumber"] == "COLO_DUP"]["is_primary_duplicate"].iloc[0])
        self.assertFalse(strat_colo[strat_colo["catalogNumber"] == "NCU_DUP"]["is_primary_duplicate"].iloc[0])

    def test_cli_geographic_scope_and_run_pipeline_shell_syntax(self):
        """Verify new CLI flags in 01_voucher_harvester.py and shell syntax/synonym generalization in run_pipeline.sh."""
        import importlib
        import subprocess

        mod_cli = importlib.import_module("scripts.data_prep.01_voucher_harvester")
        parser = mod_cli.build_cli_parser()

        defaults = parser.parse_args([])
        self.assertEqual(defaults.geographic_scope, "eastern_na")
        self.assertFalse(defaults.include_western)
        self.assertTrue(defaults.expand_synonyms)
        self.assertTrue(defaults.exclude_western)

        na_args = parser.parse_args([
            "--geographic-scope", "north_america",
            "--countries", "US", "CA", "MX",
            "--include-western",
        ])
        self.assertEqual(na_args.geographic_scope, "north_america")
        self.assertTrue(na_args.include_western)
        self.assertEqual(na_args.countries, ["US", "CA", "MX"])

        run_sh = Path(__file__).resolve().parents[2] / "run_pipeline.sh"
        sh_text = run_sh.read_text(encoding="utf-8")
        self.assertNotIn('if [[ "${FOCAL_TAXON}" == "Packera dubia" ]]', sh_text)
        self.assertNotIn("--all-eastern-taxa", sh_text)
        self.assertNotIn("ALL_EASTERN_TAXA", sh_text)
        self.assertNotIn("mod05.TARGET_TAXA", sh_text)
        self.assertNotIn('"Packera dubia" "Packera anonyma" "Packera paupercula" "Packera plattensis"', sh_text)
        self.assertIn("Packera_Taxonomic_Revision_Summary.md", sh_text)
        self.assertIn("fig06a_integrative_taxonomic_synthesis.png", sh_text)
        self.assertIn("select_cascading_outgroup_anchors", sh_text)
        self.assertIn("--geographic-scope", sh_text)
        self.assertIn("--include-western", sh_text)

        proc = subprocess.run(["bash", "-n", str(run_sh)], capture_output=True, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"bash -n run_pipeline.sh failed: {proc.stderr}")

    def test_genus_wide_packera_harvest_when_taxa_empty(self):
        """Running VoucherHarvester with taxa=[] dynamically discovers and harvests genus-wide Packera species."""
        import scripts.core.config as core_cfg
        import scripts.core.harvester as core_harv
        from scripts.core.harvester import discover_genus_packera_species_gbif

        self.assertEqual(core_cfg.DEFAULT_TARGET_TAXA, [])
        cfg = core_cfg.PipelineConfig.from_yaml()
        self.assertEqual(cfg.taxa.target_species, [])
        self.assertEqual(cfg.taxa.synonyms, {})
        self.assertEqual(cfg.taxa.outgroups, [])
        self.assertFalse(hasattr(core_harv, "EASTERN_NA_PACKERA_TAXA"))
        self.assertFalse(hasattr(core_harv, "EASTERN_NA_GBIF_QUERY_TAXA"))

        # 1. Dynamic genus-wide species discovery via GBIF Species Backbone + name_lookup/name_usage
        mock_species = MagicMock()
        mock_species.name_backbone.return_value = {
            "usageKey": 3109225,
            "canonicalName": "Packera",
            "rank": "GENUS",
            "matchType": "EXACT",
        }
        mock_species.name_lookup.return_value = {
            "results": [
                {"canonicalName": "Packera bolanderi", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
                {"canonicalName": "Packera cardamine", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
            ]
        }
        mock_species.name_usage.return_value = {
            "results": [
                {"canonicalName": "Packera franciscana", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
                {"canonicalName": "Senecio bolanderi", "rank": "SPECIES", "taxonomicStatus": "SYNONYM"},
                {"canonicalName": "Packera dubitabilis", "rank": "SPECIES", "taxonomicStatus": "DOUBTFUL"},
            ]
        }

        queried_names = []

        def _fake_occ_by_species(**kwargs):
            if "typeStatus" in kwargs:
                return {"results": [], "count": 0}
            sci = kwargs.get("scientificName")
            queried_names.append(sci)
            sci_ord = len(dict.fromkeys(queried_names))
            recs = [
                {
                    "key": sci_ord * 10000 + idx + 1,
                    "scientificName": f"{sci} Auth.",
                    "stateProvince": "North Carolina",
                    "decimalLatitude": 35.80 + idx * 0.05,
                    "decimalLongitude": -79.50 - idx * 0.05,
                    "coordinateUncertaintyInMeters": 100.0,
                    "year": 2021,
                    "month": 5,
                    "day": 14 + idx,
                    "recordedBy": "R. Collector",
                    "recordNumber": str(sci_ord * 1000 + idx),
                    "media": [{"identifier": f"https://example.com/{sci.replace(' ', '_')}_{idx}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                    "catalogNumber": f"CAT_{sci.replace(' ', '_').upper()}_{idx}",
                    "institutionCode": "TEX",
                }
                for idx in range(3)
            ]
            return {
                "results": recs,
                "count": len(recs),
            }

        with patch("scripts.core.harvester.species", mock_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_by_species):
            discovered = discover_genus_packera_species_gbif()
            self.assertEqual(discovered, ["Packera bolanderi", "Packera cardamine", "Packera franciscana"])

            harvester = VoucherHarvester(taxa=[], max_records_per_taxon=5, expand_synonyms=False)
            self.assertEqual(harvester.taxa, ["Packera bolanderi", "Packera cardamine", "Packera franciscana"])
            df_genus = harvester.harvest()

        self.assertEqual(list(dict.fromkeys(queried_names)), ["Packera bolanderi", "Packera cardamine", "Packera franciscana"])
        self.assertEqual(len(df_genus), 9)
        self.assertTrue(all(str(p).endswith((".webp", ".jpg")) for p in df_genus["image_path"]))
        for unrequested in ("Packera dubia", "Packera anonyma", "Packera paupercula", "Packera plattensis"):
            self.assertNotIn(unrequested, queried_names)

        # 2. Direct genus-level GBIF occurrence query (taxa=["Packera"]) filters indets and resolves Senecio synonyms
        mock_direct_species = MagicMock()

        def _backbone_dispatch(**kwargs):
            if kwargs.get("rank") == "genus":
                return {"usageKey": 3109225, "canonicalName": "Packera", "matchType": "EXACT"}
            if "cardamine" in str(kwargs.get("name", "")).lower():
                return {
                    "usageKey": 777111,
                    "canonicalName": "Senecio cardamine",
                    "species": "Packera cardamine",
                    "accepted": "Packera cardamine (Greene) W.A.Weber & Á.Löve",
                }
            return {}

        mock_direct_species.name_backbone.side_effect = _backbone_dispatch
        mock_direct_species.name_usage.return_value = {"results": []}

        direct_search_calls = []

        def _fake_direct_genus_occ(**kwargs):
            direct_search_calls.append(dict(kwargs))
            return {
                "results": [
                    {
                        "key": 9001,
                        "scientificName": "Packera hartiana (A.Heller) W.A.Weber & Á.Löve",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.1,
                        "decimalLongitude": -95.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 4,
                        "day": 20,
                        "media": [{"identifier": "https://example.com/hartiana.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "TEX_HART_01",
                        "institutionCode": "TEX",
                    },
                    {
                        "key": 9002,
                        "scientificName": "Senecio cardamine Greene",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.3,
                        "decimalLongitude": -95.3,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 4,
                        "day": 22,
                        "media": [{"identifier": "https://example.com/cardamine.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "TEX_CARD_02",
                        "institutionCode": "TEX",
                    },
                    {
                        "key": 9003,
                        "scientificName": "Packera sp.",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.4,
                        "decimalLongitude": -95.4,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 4,
                        "day": 25,
                        "media": [{"identifier": "https://example.com/indet.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "TEX_INDET_03",
                        "institutionCode": "TEX",
                    },
                ],
                "count": 3,
            }

        core_harv._DYNAMIC_SYNONYM_TO_CANONICAL.pop("senecio cardamine", None)
        core_harv._GBIF_SYNONYM_CACHE.pop("senecio cardamine", None)
        with patch("scripts.core.harvester.species", mock_direct_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_direct_genus_occ):
            harvester_direct = VoucherHarvester(taxa=["Packera"], max_records_per_taxon=10)
            df_direct = harvester_direct.harvest()

        self.assertEqual(direct_search_calls[0].get("taxonKey"), 3109225)
        self.assertEqual(set(df_direct["catalogNumber"]), {"TEX_HART_01", "TEX_CARD_02"})
        card_row = df_direct[df_direct["catalogNumber"] == "TEX_CARD_02"].iloc[0]
        self.assertEqual(card_row["scientificName"], "Packera cardamine")
        self.assertEqual(card_row["species_raw"], "Senecio cardamine Greene")

        # 3. Fallback when taxa=[] and species discovery returns empty: queries genus-level Packera directly
        mock_direct_species.name_lookup.return_value = {"results": []}
        direct_search_calls.clear()
        with patch("scripts.core.harvester.species", mock_direct_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_direct_genus_occ):
            harvester_fallback = VoucherHarvester(taxa=[], max_records_per_taxon=10)
            self.assertEqual(harvester_fallback.taxa, ["Packera"])
            df_fb = harvester_fallback.harvest()
        self.assertEqual(direct_search_calls[0].get("taxonKey"), 3109225)
        self.assertEqual(set(df_fb["catalogNumber"]), {"TEX_HART_01", "TEX_CARD_02"})

        # 4. Dynamic genus-wide standardize_packera_taxon edge cases (hybrid markers & open nomenclature)
        self.assertEqual(core_harv.standardize_packera_taxon("Packera ×memmingeri (Britton) Weakley"), "Packera memmingeri")
        self.assertEqual(core_harv.standardize_packera_taxon("Packera x memmingeri"), "Packera memmingeri")
        self.assertIsNone(core_harv.standardize_packera_taxon("Packera x"))
        self.assertIsNone(core_harv.standardize_packera_taxon("Packera ×"))
        self.assertIsNone(core_harv.standardize_packera_taxon("Packera sp."))

    def test_explicit_taxa_harvests_strictly_requested_species(self):
        """Running VoucherHarvester with explicit taxa harvests strictly those species without injecting Packera dubia."""
        queried_taxa = []

        def _fake_occ(**kwargs):
            if "typeStatus" in kwargs:
                return {"results": [], "count": 0}
            sci = kwargs.get("scientificName")
            queried_taxa.append(sci)
            return {
                "results": [
                    {
                        "key": 7101 if "aurea" in str(sci) else 7102,
                        "scientificName": sci,
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.91,
                        "decimalLongitude": -79.05,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2022,
                        "month": 4,
                        "day": 18,
                        "media": [{"identifier": f"https://example.com/{sci.replace(' ', '_')}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": f"NCU_{sci.split()[-1].upper()}_01",
                        "institutionCode": "NCU",
                    }
                ],
                "count": 1,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ):
            harvester = VoucherHarvester(
                taxa=["Packera aurea", "Packera obovata"],
                max_records_per_taxon=5,
                expand_synonyms=False,
            )
            self.assertEqual(harvester.taxa, ["Packera aurea", "Packera obovata"])
            df = harvester.harvest()

        self.assertEqual(list(dict.fromkeys(queried_taxa)), ["Packera aurea", "Packera obovata"])
        self.assertEqual(set(df["scientificName"]), {"Packera aurea", "Packera obovata"})
        for forbidden in ("Packera dubia", "Packera anonyma", "Packera paupercula", "Packera plattensis"):
            self.assertNotIn(forbidden, queried_taxa)
            self.assertNotIn(forbidden, set(df["scientificName"]))

        # Verify main.run_harvest omits --taxa when none supplied and passes strictly explicit --taxa
        import argparse
        import main as root_main

        cfg_default = root_main.PipelineConfig.from_yaml()
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_csv = Path(tmp_d) / "curated_vouchers.csv"
            tmp_csv.write_text("catalogNumber,scientificName\nNCU1,Packera aurea\n", encoding="utf-8")
            base_ns = dict(
                run_id=None,
                out_dir=Path(tmp_d) / "raw",
                output_csv=tmp_csv,
                max_records=10,
                max_uncertainty=1000.0,
                min_megapixels=8.0,
                min_file_size_kb=500.0,
                concurrency=4,
                download_images=False,
                force=False,
                verbose=False,
            )

            with patch.object(root_main.subprocess, "run") as mock_sub_empty:
                root_main.run_harvest(argparse.Namespace(taxa=None, **base_ns), cfg_default)
                cmd_empty = mock_sub_empty.call_args[0][0]
                self.assertNotIn("--taxa", cmd_empty)
                self.assertNotIn("--taxon", cmd_empty)

            with patch.object(root_main.subprocess, "run") as mock_sub_explicit:
                root_main.run_harvest(
                    argparse.Namespace(taxa=["Packera aurea", "Packera obovata"], **base_ns),
                    cfg_default,
                )
                cmd_explicit = mock_sub_explicit.call_args[0][0]
                self.assertTrue("--taxa" in cmd_explicit or "--taxon" in cmd_explicit)
                taxa_flag = "--taxa" if "--taxa" in cmd_explicit else "--taxon"
                taxa_idx = cmd_explicit.index(taxa_flag)
                self.assertEqual(cmd_explicit[taxa_idx + 1 : taxa_idx + 3], ["Packera aurea", "Packera obovata"])
                self.assertNotIn("Packera dubia", cmd_explicit)

            main_parser = root_main.build_parser()
            parsed_harvest = main_parser.parse_args(["harvest", "--taxon", "Packera aurea", "Packera obovata"])
            self.assertEqual(parsed_harvest.taxa, ["Packera aurea", "Packera obovata"])
            parsed_all_empty = main_parser.parse_args(["run-all"])
            self.assertIsNone(parsed_all_empty.taxa)

    def test_cross_stage_cli_table_handoffs_in_main_and_run_pipeline(self):
        """Verify main.py and run_pipeline.sh wire all cross-stage CLI table handoffs across Steps 04, 06, 07, and 08."""
        import os
        import subprocess
        import main as root_main

        def _arg_value(cmd_list, flag_name):
            self.assertIn(flag_name, cmd_list, f"Missing CLI flag {flag_name!r} in command: {cmd_list}")
            idx = cmd_list.index(flag_name)
            self.assertLess(idx + 1, len(cmd_list), f"Flag {flag_name!r} has no value in command: {cmd_list}")
            return cmd_list[idx + 1]

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_run = Path(tmp_d) / "run_handoff_audit"
            tables_dir = tmp_run / "data" / "tables"
            figures_dir = tmp_run / "outputs" / "figures"
            contours_dir = tmp_run / "data" / "contours"
            tables_dir.mkdir(parents=True, exist_ok=True)
            figures_dir.mkdir(parents=True, exist_ok=True)
            contours_dir.mkdir(parents=True, exist_ok=True)
            (contours_dir / "NCU001_leaf1.csv").write_text("x,y\n0,0\n1,1\n", encoding="utf-8")

            cfg_scoped = root_main.PipelineConfig.from_yaml().with_run_dir(tmp_run)
            parser = root_main.build_parser()

            # 1. Verify Step 03 -> Step 04 handoff in main.run_morphometrics()
            morph_args = parser.parse_args([
                "morphometrics",
                "--focal-taxon", "Packera aurea",
                "--intraspecific-mode",
            ])
            with patch.dict(os.environ, {"PACKERA_RUN_DIR": ""}, clear=False), \
                 patch.object(root_main.shutil, "which", return_value="/usr/bin/Rscript"), \
                 patch.object(root_main, "verify_file_exists"), \
                 patch.object(root_main, "verify_dir_has_files"), \
                 patch.object(root_main.subprocess, "run") as mock_morph_sub:
                root_main.run_morphometrics(morph_args, cfg_scoped)
                self.assertEqual(mock_morph_sub.call_count, 2)
                cmd_step04 = mock_morph_sub.call_args_list[1][0][0]

            self.assertTrue(any("04_gmm_morphotools.R" in tok for tok in cmd_step04))
            self.assertEqual(_arg_value(cmd_step04, "--reproductive"), str(tables_dir / "voucher_reproductive_metrics.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--allometry-report"), str(tables_dir / "allometry_permanova_summary.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--allometry-pairwise"), str(tables_dir / "allometry_pairwise_slopes.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--allometry-plot"), str(figures_dir / "fig04a_foliar_allometry_audit.png"))
            self.assertEqual(_arg_value(cmd_step04, "--output-batch-audit"), str(tables_dir / "batch_variance_audit.csv"))
            self.assertEqual(_arg_value(cmd_step04, "--focal-taxon"), "Packera aurea")
            self.assertIn("--intraspecific-mode", cmd_step04)

            # 2. Verify Step 06 and Step 07 handoffs in main.run_synthesis()
            synth_args = parser.parse_args([
                "synthesis",
                "--focal-taxon", "Packera aurea",
                "--intraspecific-mode",
            ])
            with patch.dict(os.environ, {"PACKERA_RUN_DIR": ""}, clear=False), \
                 patch.object(root_main.shutil, "which", return_value="/usr/bin/Rscript"), \
                 patch.object(root_main, "verify_file_exists"), \
                 patch.object(root_main.subprocess, "run") as mock_synth_sub:
                root_main.run_synthesis(synth_args, cfg_scoped)
                self.assertEqual(mock_synth_sub.call_count, 2)
                cmd_step06 = mock_synth_sub.call_args_list[0][0][0]
                cmd_step07 = mock_synth_sub.call_args_list[1][0][0]

            # Step 05 assertions
            self.assertTrue(any("05_multimodal_spatial_rf.R" in tok for tok in cmd_step06))
            reports_dir = tmp_run / "outputs" / "reports"
            self.assertEqual(_arg_value(cmd_step06, "--pheno-states"), str(tables_dir / "voucher_phenological_states.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--ssurgo-table"), str(tables_dir / "ssurgo_edaphic_validation.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--contingency-table"), str(reports_dir / "micro_edaphic_outcrop_contingency.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--pheno-anomalies"), str(tables_dir / "phenological_anomalies.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--pheno-summary"), str(reports_dir / "phenological_anomaly_summary.csv"))
            self.assertEqual(_arg_value(cmd_step06, "--pheno-plot"), str(figures_dir / "fig05b_hopkins_circular_phenology.png"))
            self.assertEqual(_arg_value(cmd_step06, "--focal-taxon"), "Packera aurea")
            self.assertIn("--intraspecific-mode", cmd_step06)

            # Unified Step 06 synthesis + complex delimitation assertions
            self.assertTrue(any("06_triage_dashboard_synthesis.R" in tok for tok in cmd_step07))
            self.assertEqual(_arg_value(cmd_step07, "--harmonics"), str(tables_dir / "leaf_efa_harmonics.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-queue"), str(tables_dir / "triage_queue.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-html"), str(reports_dir / "Triage_Dashboard_Interactive.html"))
            self.assertEqual(_arg_value(cmd_step07, "--output-complexes"), str(reports_dir / "empirical_taxon_complexes.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-affinity"), str(reports_dir / "pairwise_taxon_affinity_matrix.csv"))
            self.assertEqual(_arg_value(cmd_step07, "--output-neighbors"), str(reports_dir / "focal_auto_neighbors.txt"))
            self.assertEqual(_arg_value(cmd_step07, "--output-complex-plot"), str(figures_dir / "fig06b_taxon_complex_delimitation_network.png"))
            self.assertEqual(_arg_value(cmd_step07, "--focal-taxon"), "Packera aurea")
            self.assertIn("--intraspecific-mode", cmd_step07)

        # 3. Verify run_pipeline.sh unified 6-stage ordering (Step 05 -> unified Step 06 with complex discovery)
        run_sh = Path(__file__).resolve().parents[2] / "run_pipeline.sh"
        sh_text = run_sh.read_text(encoding="utf-8")

        pos_step05 = sh_text.index("scripts/analysis/05_multimodal_spatial_rf.R")
        pos_step06 = sh_text.index("scripts/analysis/06_triage_dashboard_synthesis.R")
        self.assertLess(pos_step05, pos_step06, "Unified Step 06 must run AFTER Step 05 completes in run_pipeline.sh")

        step06_block = sh_text[pos_step06:]
        self.assertIn('--multimodal-flags "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv"', step06_block)
        self.assertIn('--output-complexes "${RUN_DIR}/outputs/reports/empirical_taxon_complexes.csv"', step06_block)
        self.assertIn('--output-complex-plot "${RUN_DIR}/outputs/figures/fig06b_taxon_complex_delimitation_network.png"', step06_block)

        proc = subprocess.run(["bash", "-n", str(run_sh)], capture_output=True, text=True, check=False)
        self.assertEqual(proc.returncode, 0, f"bash -n run_pipeline.sh failed: {proc.stderr}")


class TestPhase1AuthorityHygieneAndExsiccataeConflict(unittest.TestCase):
    """Tests for collector-as-determiner leakage prevention, corpus authority duplicate screening, and exsiccatae conflict auditing."""

    def test_collector_as_determiner_leakage_prevented_in_parse_determiner_tier(self):
        """Specialist in recordedBy must NOT promote a generalist/student identifiedBy sheet to Tier_1_Gold unless identifiedBy is blank or references collector."""
        # 1. Specialist collected, student determined at local herbarium -> Tier_3_Bronze (NOT Tier_1_Gold)
        tier_student_local, _, det_raw_1 = parse_determiner_tier(
            type_status_raw="",
            identified_by_raw="J. Undergrad Student",
            recorded_by_raw="T.M. Barkley",
            history_raw="",
            institution_code_raw="UNCC",
            locality_raw="Granite outcrop near Concord",
            habitat_raw="Shallow sandy soil",
        )
        self.assertEqual(tier_student_local, "Tier_3_Bronze")
        self.assertEqual(det_raw_1, "J. Undergrad Student")

        # 2. Specialist collected, generalist determined at major research herbarium (NCU) -> Tier_2_Silver (NOT Tier_1_Gold)
        tier_generalist_ncu, _, _ = parse_determiner_tier(
            type_status_raw="",
            identified_by_raw="J. Generalist",
            recorded_by_raw="T.M. Barkley",
            history_raw="",
            institution_code_raw="NCU",
            locality_raw="Wake Co., granite flatrock 3 km E of Rolesville",
            habitat_raw="Granitic outcrop seepage margin",
            dynamic_specialists={"mahoney", "kowal"},
        )
        self.assertEqual(tier_generalist_ncu, "Tier_2_Silver")

        # 3. Dynamic corpus specialist in recordedBy, student in identifiedBy -> NOT Tier_1_Gold
        tier_dyn_leak, _, _ = parse_determiner_tier(
            type_status_raw="",
            identified_by_raw="Field Intern",
            recorded_by_raw="A.M. Mahoney",
            history_raw="",
            institution_code_raw="UNCC",
            locality_raw="Roadside meadow",
            habitat_raw="",
            dynamic_specialists={"mahoney"},
        )
        self.assertNotEqual(tier_dyn_leak, "Tier_1_Gold")

        # 4. Specialist collected, identifiedBy is blank or explicitly references collector -> Tier_1_Gold
        for id_by in ["", "collector", "id. by collector", "same"]:
            with self.subTest(identified_by=id_by):
                tier_fb, _, _ = parse_determiner_tier(
                    type_status_raw="",
                    identified_by_raw=id_by,
                    recorded_by_raw="T.M. Barkley",
                    history_raw="",
                    institution_code_raw="UNCC",
                    locality_raw="Prairie slope",
                    habitat_raw="",
                )
                self.assertEqual(tier_fb, "Tier_1_Gold")

    def test_discover_corpus_authorities_disqualifies_high_duplicate_disagreement(self):
        """Candidate corpus specialists with >25% cross-duplicate taxonomic disagreement across shared exsiccatae_key groups are disqualified."""
        from scripts.core.harvester import discover_corpus_authorities

        rows = []
        insts = ["NCU", "NY", "MO", "GH", "US", "BRIT"]
        # Consistent determiner: 6 determinations across 6 herbaria, 3 duplicate groups all concordant (0% disagreement)
        for i in range(6):
            grp_num = 100 + (i // 2)
            rows.append({
                "catalogNumber": f"CONS-{i:02d}",
                "institutionCode": insts[i],
                "scientificName": "Packera obovata",
                "identifiedBy": "E.R. Consistent",
                "recordedBy": "R. Collector",
                "recordNumber": str(grp_num),
                "eventDate": "2020-05-10",
            })

        # Noisy determiner: 6 determinations of Packera aurea across 6 herbaria,
        # participating in 3 duplicate groups where 2 of 3 (66.7% > 25%) have duplicate sheets determined as Packera obovata
        for i in range(6):
            grp_num = 200 + (i % 3)
            rows.append({
                "catalogNumber": f"NOISY-{i:02d}",
                "institutionCode": insts[i],
                "scientificName": "Packera aurea",
                "identifiedBy": "B.L. Noisy",
                "recordedBy": "R. Collector",
                "recordNumber": str(grp_num),
                "eventDate": "2021-04-18",
            })
        # Add conflicting duplicate partners for groups 200 and 201 (2 of 3 duplicate groups conflicted = 66.7%)
        rows.append({
            "catalogNumber": "DUP-PARTNER-200",
            "institutionCode": "TENN",
            "scientificName": "Packera obovata",
            "identifiedBy": "T.M. Barkley",
            "recordedBy": "R. Collector",
            "recordNumber": "200",
            "eventDate": "2021-04-18",
        })
        rows.append({
            "catalogNumber": "DUP-PARTNER-201",
            "institutionCode": "DUKE",
            "scientificName": "Packera obovata",
            "identifiedBy": "T.M. Barkley",
            "recordedBy": "R. Collector",
            "recordNumber": "201",
            "eventDate": "2021-04-18",
        })

        dyn_specs, _ = discover_corpus_authorities(pd.DataFrame(rows))
        self.assertIn("consistent", dyn_specs)
        self.assertNotIn("noisy", dyn_specs)

    def test_stratify_duplicates_conflict_elevation_sqlite_persistence_and_26_col_csv(self):
        """Conflicting Tier_3_Bronze home-herbarium primary duplicate elevates Tier_1_Gold duplicate determination, logs to SQLite, and preserves 26-column CSV."""
        import sqlite3
        from scripts.core.harvester import EXPORT_COLUMNS

        df_in = pd.DataFrame([
            {
                "catalogNumber": "NCU-DUP-1",
                "institutionCode": "NCU",
                "scientificName": "Packera aurea",
                "species_raw": "Packera aurea",
                "recordedBy": "J. B. Fuller",
                "recordNumber": "1042",
                "identifiedBy": "Field Student",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "",
                "county": "Wake",
                "stateProvince": "North Carolina",
                "decimalLatitude": 35.85,
                "decimalLongitude": -78.55,
                "coordinateUncertainty": 50.0,
                "year": 2022,
                "month": 5,
                "day": 1,
                "eventDate": "2022-05-01",
                "regional_group": "Piedmont",
                "image_path": "",
            },
            {
                "catalogNumber": "MO-DUP-2",
                "institutionCode": "MO",
                "scientificName": "Packera obovata",
                "species_raw": "Packera obovata",
                "recordedBy": "J. B. Fuller",
                "recordNumber": "1042",
                "identifiedBy": "T.M. Barkley",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "",
                "county": "Wake",
                "stateProvince": "North Carolina",
                "decimalLatitude": 35.85,
                "decimalLongitude": -78.55,
                "coordinateUncertainty": 50.0,
                "year": 2022,
                "month": 5,
                "day": 1,
                "eventDate": "2022-05-01",
                "regional_group": "Piedmont",
                "image_path": "",
            },
        ])

        strat_df = stratify_duplicates(df_in, home_herbarium="NCU")
        self.assertTrue((strat_df["exsiccatae_gold_taxon"] == "Packera obovata").all())
        mo_gold = strat_df[strat_df["catalogNumber"] == "MO-DUP-2"].iloc[0]
        self.assertFalse(bool(mo_gold["exsiccatae_taxonomic_conflict"]))

        prim = strat_df[strat_df["catalogNumber"] == "NCU-DUP-1"].iloc[0]
        self.assertTrue(bool(prim["exsiccatae_taxonomic_conflict"]))
        self.assertTrue(bool(prim["is_primary_duplicate"]))
        self.assertEqual(prim["scientificName"], "Packera obovata")
        self.assertEqual(prim["determiner_tier"], "Tier_2_Silver")
        # Verbatim original taxon field and exsiccatae provenance remain intact for auditability
        self.assertEqual(prim["species_raw"], "Packera aurea")
        self.assertIn("[exsiccatae Tier_1_Gold duplicate MO-DUP-2]", prim["identifiedBy"])

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            out_csv = tmp_p / "curated_vouchers.csv"
            db_path = tmp_p / "test_packera.sqlite"
            export_curated_table(strat_df, out_csv, home_herbarium="NCU", db_path=db_path)

            csv_df = pd.read_csv(out_csv)
            from scripts.core.harvester import EXPORT_COLUMNS
            self.assertEqual(list(csv_df.columns), EXPORT_COLUMNS)
            self.assertEqual(len(csv_df.columns), 23)

            conn = sqlite3.connect(str(db_path))
            try:
                row_db = conn.execute(
                    "SELECT scientificName, determiner_tier, exsiccatae_taxonomic_conflict, exsiccatae_gold_taxon "
                    "FROM vouchers WHERE catalogNumber = 'NCU-DUP-1'"
                ).fetchone()
            finally:
                conn.close()
            self.assertIsNotNone(row_db)
            self.assertEqual(row_db[0], "Packera obovata")
            self.assertEqual(row_db[1], "Tier_2_Silver")
            self.assertEqual(int(row_db[2]), 1)
            self.assertEqual(row_db[3], "Packera obovata")

    def test_classify_nomenclatural_type_rank_and_cascading_outgroup_type_priority(self):
        """Verify classify_nomenclatural_type_rank(), classify_nomenclatural_type_subrank(), type-annotator leakage guard, and cascading outgroup type priority."""
        from scripts.core.harvester import (
            classify_nomenclatural_type_rank,
            classify_nomenclatural_type_subrank,
            discover_corpus_authorities,
            select_cascading_outgroup_anchors,
        )

        self.assertEqual(classify_nomenclatural_type_rank("Holotype of Senecio dubius")[0], 2)
        self.assertEqual(classify_nomenclatural_type_rank("Isolectotype")[0], 2)
        self.assertEqual(classify_nomenclatural_type_rank("Paratype")[0], 1)
        self.assertEqual(classify_nomenclatural_type_rank("Syntype")[0], 1)
        self.assertEqual(classify_nomenclatural_type_rank("Not a type"), (0, "None"))
        self.assertEqual(classify_nomenclatural_type_rank(None), (0, "None"))

        # Verify ICN Art. 9.6 sub-ranking: Holotype > Isotype > Syntype > Paratype > Non-type
        self.assertGreater(
            classify_nomenclatural_type_subrank("Holotype"),
            classify_nomenclatural_type_subrank("Isotype"),
        )
        self.assertGreater(
            classify_nomenclatural_type_subrank("Syntype"),
            classify_nomenclatural_type_subrank("Paratype"),
        )
        self.assertEqual(classify_nomenclatural_type_subrank("Not a type"), 0.0)

        # Verify Layer 3 type-annotator leakage guard: a single-sheet annotator on 1 type specimen
        # at 1 herbarium is NOT promoted to dynamic_specialists unless multi-institution criteria are met
        leak_rows = [
            {
                "catalogNumber": "TYPE-SINGLE-1",
                "institutionCode": "NCU",
                "scientificName": "Packera aurea",
                "identifiedBy": "Q.Z. Singleannotator",
                "recordedBy": "A. Collector",
                "recordNumber": "10",
                "type_status": "Holotype",
                "eventDate": "1950-05-01",
            },
            {
                "catalogNumber": "TYPE-MULTI-1",
                "institutionCode": "NCU",
                "scientificName": "Packera aurea",
                "identifiedBy": "M.K. Multiannotator",
                "recordedBy": "B. Collector",
                "recordNumber": "20",
                "type_status": "Holotype",
                "eventDate": "1960-05-01",
            },
            {
                "catalogNumber": "NONTYPE-MULTI-2",
                "institutionCode": "NY",
                "scientificName": "Packera aurea",
                "identifiedBy": "M.K. Multiannotator",
                "recordedBy": "C. Collector",
                "recordNumber": "21",
                "type_status": "None",
                "eventDate": "1962-05-01",
            },
        ]
        dyn_specs_leak, _ = discover_corpus_authorities(pd.DataFrame(leak_rows))
        self.assertNotIn("singleannotator", dyn_specs_leak)
        self.assertIn("multiannotator", dyn_specs_leak)

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            # Create local images where non-type Gold sheets have LARGER dimensions than the Holotype/Paratype
            holo_img = tmp_p / "holo.jpg"
            para_img = tmp_p / "para.jpg"
            Image.new("RGB", (2000, 2000), color=(120, 120, 120)).save(holo_img)
            Image.new("RGB", (2100, 2100), color=(120, 120, 120)).save(para_img)

            rows = [
                {
                    "catalogNumber": "HOLO-01",
                    "scientificName": "Packera aurea",
                    "species_raw": "Packera aurea",
                    "species_standardized": "Packera aurea",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "Holotype",
                    "is_primary_duplicate": True,
                    "image_path": str(holo_img),
                },
                {
                    "catalogNumber": "PARA-01",
                    "scientificName": "Packera aurea",
                    "species_raw": "Packera aurea",
                    "species_standardized": "Packera aurea",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "Paratype",
                    "is_primary_duplicate": True,
                    "image_path": str(para_img),
                },
            ]
            for idx in range(5):
                hi_img = tmp_p / f"gold_{idx}.jpg"
                Image.new("RGB", (3500, 3500), color=(120, 120, 120)).save(hi_img)
                rows.append({
                    "catalogNumber": f"GOLD-NONTYPE-{idx:02d}",
                    "scientificName": "Packera aurea",
                    "species_raw": "Packera aurea",
                    "species_standardized": "Packera aurea",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "None",
                    "is_primary_duplicate": True,
                    "image_path": str(hi_img),
                })

            df_in = pd.DataFrame(rows)
            outgroup_dfs, counts, breakdown = select_cascading_outgroup_anchors(
                df_in,
                outgroup_taxa=["Packera aurea"],
                max_anchors=3,
                workspace_dir=tmp_p,
            )
            self.assertEqual(counts["Packera aurea"], 3)
            selected_df = outgroup_dfs[0]
            self.assertEqual(len(selected_df), 3)
            sel_cats = list(selected_df["catalogNumber"])
            self.assertEqual(sel_cats[0], "HOLO-01")
            self.assertEqual(sel_cats[1], "PARA-01")
            self.assertEqual(breakdown["Packera aurea"]["Type_Specimens"], 2)
            self.assertEqual(breakdown["Packera aurea"]["Tier_0_Primary_Type"], 1)
            self.assertEqual(breakdown["Packera aurea"]["Tier_0_Secondary_Type"], 1)

    def test_canadian_western_provinces_and_ungeoreferenced_nontype_exclusion(self):
        """Verify Canadian western/northern provinces are excluded and un-georeferenced non-types without state are excluded while types are retained."""
        canadian_western = [
            "BC", "British Columbia", "AB", "Alberta", "SK", "Saskatchewan",
            "YT", "Yukon", "Yukon Territory", "NT", "Northwest Territories", "NU", "Nunavut",
        ]
        for prov in canadian_western:
            with self.subTest(province=prov):
                self.assertTrue(
                    is_excluded_western_region(prov, exclude_western=True),
                    f"Canadian western/northern province/territory '{prov}' should be excluded.",
                )

        # Un-georeferenced non-type lacking stateProvince is excluded when exclude_western=True
        self.assertTrue(
            is_excluded_western_region(None, lat=None, lon=None, exclude_western=True, is_type_specimen=False)
        )
        self.assertTrue(
            is_excluded_western_region("", lat=np.nan, lon=np.nan, exclude_western=True, is_type_specimen=False)
        )
        # Un-georeferenced type specimen without western state/lon is retained
        self.assertFalse(
            is_excluded_western_region(None, lat=None, lon=None, exclude_western=True, is_type_specimen=True)
        )
        # Un-georeferenced non-type WITH valid eastern state/province is retained
        self.assertFalse(
            is_excluded_western_region("North Carolina", lat=None, lon=None, exclude_western=True, is_type_specimen=False)
        )
        self.assertFalse(
            is_excluded_western_region("Ontario", lat=None, lon=None, exclude_western=True, is_type_specimen=False)
        )

    def test_genus_wide_eastern_na_residency_filter_and_type_rescue(self):
        """Verify genus-wide eastern_na residency check filters out western stray species (< 3 eastern records) and rescues un-georeferenced types."""
        mock_species = MagicMock()
        mock_species.name_backbone.return_value = {
            "usageKey": 3109225,
            "canonicalName": "Packera",
            "rank": "GENUS",
            "matchType": "EXACT",
        }
        mock_species.name_lookup.return_value = {
            "results": [
                {"canonicalName": "Packera aurea", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
                {"canonicalName": "Packera bolanderi", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
            ]
        }
        mock_species.name_usage.return_value = {"results": []}

        def _fake_occ_residency(**kwargs):
            sci = kwargs.get("scientificName")
            if kwargs.get("typeStatus") == "*":
                if sci == "Packera aurea":
                    # 19th-century Holotype lacking decimalLongitude
                    return {
                        "results": [
                            {
                                "key": 88001,
                                "scientificName": "Packera aurea (L.) Á.Löve & D.Löve",
                                "typeStatus": "Holotype",
                                "stateProvince": "Pennsylvania",
                                "decimalLatitude": None,
                                "decimalLongitude": None,
                                "year": 1814,
                                "recordedBy": "F. Pursh",
                                "recordNumber": "1",
                                "media": [{"identifier": "https://example.com/aurea_holo.jpg", "type": "StillImage", "format": "image/jpeg"}],
                                "catalogNumber": "PH_AUREA_HOLOTYPE",
                                "institutionCode": "PH",
                            }
                        ],
                        "count": 1,
                    }
                return {"results": [], "count": 0}

            ctry = kwargs.get("country")
            if sci == "Packera aurea":
                if ctry == "CA":
                    recs = [
                        {
                            "key": 8150,
                            "scientificName": "Packera aurea",
                            "stateProvince": "Ontario",
                            "decimalLatitude": 44.2,
                            "decimalLongitude": -79.5,
                            "coordinateUncertaintyInMeters": 50.0,
                            "year": 2020,
                            "month": 6,
                            "day": 12,
                            "recordedBy": "C. Collector",
                            "recordNumber": "701",
                            "media": [{"identifier": "https://example.com/aurea_on.jpg", "type": "StillImage", "format": "image/jpeg"}],
                            "catalogNumber": "CAN_AUREA_ON",
                            "institutionCode": "CAN",
                        }
                    ]
                    return {"results": recs, "count": len(recs)}
                # 2 US eastern records + 1 ON CA record = 3 Eastern NA records (>= 3 residency threshold),
                # plus 1 western Montana conspecific retained under full continental range harvesting!
                recs = [
                    {
                        "key": 8100 + i,
                        "scientificName": "Packera aurea",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.8 + i * 0.1,
                        "decimalLongitude": -79.0 - i * 0.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 5,
                        "day": 10 + i,
                        "recordedBy": "A. Collector",
                        "recordNumber": str(500 + i),
                        "media": [{"identifier": f"https://example.com/aurea_{i}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": f"NCU_AUREA_{i}",
                        "institutionCode": "NCU",
                    }
                    for i in range(2)
                ] + [
                    {
                        "key": 8199,
                        "scientificName": "Packera aurea",
                        "stateProvince": "Montana",
                        "decimalLatitude": 46.8,
                        "decimalLongitude": -110.5,
                        "coordinateUncertaintyInMeters": 100.0,
                        "year": 2019,
                        "month": 7,
                        "day": 4,
                        "recordedBy": "W. Collector",
                        "recordNumber": "902",
                        "media": [{"identifier": "https://example.com/aurea_mt.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "MONT_AUREA_WEST",
                        "institutionCode": "MONT",
                    }
                ]
                return {"results": recs, "count": len(recs)}
            else:
                if ctry == "CA":
                    return {"results": [], "count": 0}
                # Western species with only 2 eastern stray records (< 3 threshold) + 2 Oregon records -> filtered out
                recs = [
                    {
                        "key": 8200 + i,
                        "scientificName": "Packera bolanderi",
                        "stateProvince": "Texas",
                        "decimalLatitude": 31.0 + i * 0.1,
                        "decimalLongitude": -97.0 - i * 0.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 5,
                        "day": 10 + i,
                        "recordedBy": "B. Collector",
                        "recordNumber": str(600 + i),
                        "media": [{"identifier": f"https://example.com/bolanderi_{i}.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": f"TEX_BOLANDERI_{i}",
                        "institutionCode": "TEX",
                    }
                    for i in range(2)
                ] + [
                    {
                        "key": 8250,
                        "scientificName": "Packera bolanderi",
                        "stateProvince": "Oregon",
                        "decimalLatitude": 44.0,
                        "decimalLongitude": -123.1,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2020,
                        "month": 6,
                        "day": 1,
                        "recordedBy": "P. Collector",
                        "recordNumber": "650",
                        "media": [{"identifier": "https://example.com/bolanderi_or.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "OSC_BOLANDERI_OR",
                        "institutionCode": "OSC",
                    }
                ]
                return {"results": recs, "count": len(recs)}

        with patch("scripts.core.harvester.species", mock_species), \
             patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_residency):
            # Default eastern_na drops Montana (-110.5W) and non-cis-Mississippi P. bolanderi while keeping cis-Mississippi + Holotype
            harvester = VoucherHarvester(taxa=[], max_records_per_taxon=10, expand_synonyms=False)
            df = harvester.harvest()
            self.assertEqual(set(df["species_standardized"]), {"Packera aurea"})
            self.assertNotIn("Packera bolanderi", set(df["species_standardized"]))
            self.assertIn("PH_AUREA_HOLOTYPE", set(df["catalogNumber"]))
            self.assertIn("CAN_AUREA_ON", set(df["catalogNumber"]))
            self.assertNotIn("MONT_AUREA_WEST", set(df["catalogNumber"]))

            # Explicit harvest_full_continental_range=True opt-in still retains Montana conspecific
            harvester_full = VoucherHarvester(
                taxa=[],
                max_records_per_taxon=10,
                expand_synonyms=False,
                harvest_full_continental_range=True,
            )
            df_full = harvester_full.harvest()
            self.assertIn("MONT_AUREA_WEST", set(df_full["catalogNumber"]))

    def test_harvest_deduplicates_secondary_exsiccatae_and_migrates_jpg_to_webp(self):
        """Verify harvest() deduplicates secondary exsiccatae before download and migrate_existing_vouchers_to_webp converts .jpg to .webp."""
        from scripts.core.harvester import migrate_existing_vouchers_to_webp, resolve_voucher_image_path

        def _fake_occ_dups(**kwargs):
            if "typeStatus" in kwargs:
                return {"results": [], "count": 0}
            # Two duplicate sheets of the exact same collection event (Fuller #777 on 2021-05-01)
            return {
                "results": [
                    {
                        "key": 9101,
                        "scientificName": "Packera dubia",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.0,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2021,
                        "month": 5,
                        "day": 1,
                        "recordedBy": "J. B. Fuller",
                        "recordNumber": "777",
                        "media": [{"identifier": "https://example.com/ncu_dup.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "NCU_PRIMARY_777",
                        "institutionCode": "NCU",
                    },
                    {
                        "key": 9102,
                        "scientificName": "Packera dubia",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.0,
                        "coordinateUncertaintyInMeters": 50.0,
                        "year": 2021,
                        "month": 5,
                        "day": 1,
                        "recordedBy": "J. Brandon Fuller",
                        "recordNumber": "#777",
                        "media": [{"identifier": "https://example.com/ga_dup.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "GA_SECONDARY_777",
                        "institutionCode": "GA",
                    },
                ],
                "count": 2,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_dups):
            harvester = VoucherHarvester(taxa=["Packera dubia"], max_records_per_taxon=10, expand_synonyms=False)
            df_dedup = harvester.harvest()

        self.assertEqual(len(df_dedup), 1)
        self.assertEqual(df_dedup.iloc[0]["catalogNumber"], "NCU_PRIMARY_777")
        self.assertEqual(df_dedup.iloc[0]["image_path"], "data/raw_vouchers/NCU_PRIMARY_777.webp")

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            raw_dir = tmp_p / "data" / "raw_vouchers"
            raw_dir.mkdir(parents=True, exist_ok=True)
            jpg_file = raw_dir / "NCU_LEGACY_001.jpg"
            Image.new("RGB", (600, 800), color=(90, 140, 100)).save(jpg_file, format="JPEG")

            csv_file = tmp_p / "curated_vouchers.csv"
            pd.DataFrame([
                {
                    "catalogNumber": "NCU_LEGACY_001",
                    "image_path": "data/raw_vouchers/NCU_LEGACY_001.jpg",
                }
            ]).to_csv(csv_file, index=False)

            migrated_n = migrate_existing_vouchers_to_webp(raw_dir, curated_csv=csv_file)
            self.assertEqual(migrated_n, 1)
            self.assertFalse(jpg_file.exists())
            webp_file = raw_dir / "NCU_LEGACY_001.webp"
            self.assertTrue(webp_file.is_file())
            self.assertGreater(webp_file.stat().st_size, 0)

            updated_csv = pd.read_csv(csv_file)
            self.assertEqual(updated_csv.iloc[0]["image_path"], "data/raw_vouchers/NCU_LEGACY_001.webp")
            resolved = resolve_voucher_image_path("data/raw_vouchers/NCU_LEGACY_001.jpg", workspace_dir=tmp_p)
            self.assertEqual(resolved, webp_file)

    def test_nomenclatural_adjustment_isolation_and_institutional_batch_demotion(self):
        """
        Verify that:
          1. `identifiedBy = 'Nomenclatural Adjustment'` is not promoted to `Tier_1_Gold`
             even when `recordedBy` contains a specialist surname.
          2. Undated institutional batch determinations (`>= 25` sheets across `>= 10`
             collectors at one herbarium) are demoted to `Tier_2_Silver` for historical
             collectors while personal collections (`recordedBy == identifiedBy`) remain
             `Tier_1_Gold`.
        """
        from scripts.core.taxonomy_authority import (
            apply_corpus_authority_tiers,
            parse_determiner_tier,
        )

        tier_adj, _, _ = parse_determiner_tier(
            type_status_raw="None",
            identified_by_raw="Nomenclatural Adjustment",
            recorded_by_raw="D. K. Trock",
            history_raw="",
            institution_code_raw="NCU",
            locality_raw="Orange County, Chapel Hill, rich mesic slope",
            habitat_raw="Deciduous forest",
            dynamic_specialists={"trock", "weakley"},
            dynamic_herbaria={"NCU"},
        )
        self.assertNotEqual(tier_adj, "Tier_1_Gold")

        batch_rows = [
            {
                "catalogNumber": "NCU_NOM_ADJ_01",
                "institutionCode": "NCU",
                "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                "species_raw": "Packera dubia",
                "type_status": "None",
                "identifiedBy": "Nomenclatural Adjustment",
                "recordedBy": "A. S. Weakley",
                "dateIdentified": "",
                "locality": "Orange County, Chapel Hill, rich mesic slope",
                "habitat": "Deciduous forest",
            }
        ]
        for idx in range(26):
            batch_rows.append(
                {
                    "catalogNumber": f"NCU_BATCH_HIST_{idx:02d}",
                    "institutionCode": "NCU",
                    "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species_raw": "Packera dubia",
                    "type_status": "None",
                    "identifiedBy": "D. K. Trock",
                    "recordedBy": f"J. HistoricalCollector{idx:02d}",
                    "dateIdentified": "",
                    "locality": "Orange County, North Carolina, alluvial woodland",
                    "habitat": "Bottomland forest",
                }
            )
        batch_rows.append(
            {
                "catalogNumber": "NCU_BATCH_PERSONAL_01",
                "institutionCode": "NCU",
                "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                "species_raw": "Packera dubia",
                "type_status": "None",
                "identifiedBy": "D. K. Trock",
                "recordedBy": "D. K. Trock",
                "dateIdentified": "",
                "locality": "Orange County, North Carolina, alluvial woodland",
                "habitat": "Bottomland forest",
            }
        )

        df_scored = apply_corpus_authority_tiers(pd.DataFrame(batch_rows))
        tier_by_cat = dict(zip(df_scored["catalogNumber"], df_scored["determiner_tier"]))

        self.assertNotEqual(tier_by_cat["NCU_NOM_ADJ_01"], "Tier_1_Gold")
        for idx in range(26):
            self.assertEqual(tier_by_cat[f"NCU_BATCH_HIST_{idx:02d}"], "Tier_2_Silver")
        self.assertEqual(tier_by_cat["NCU_BATCH_PERSONAL_01"], "Tier_1_Gold")

    def test_north_american_packera_specialists_and_monographer_loan_protection(self):
        """
        Verify that:
          1. Newly added North American Packera specialists are classified as Tier_1_Gold.
          2. apply_corpus_authority_tiers() preserves Tier_1_Gold for a specialist
             annotating >= 30 sheets across multiple species or herbaria even when
             dateIdentified is missing, while demoting an undated single-herbarium,
             single-species batch transfer (>= 25 sheets, >= 10 collectors, >= 85%
             in 1 species at 1 institution where recordedBy != identifiedBy) to Tier_2_Silver.
        """
        import inspect
        import scripts.core.taxonomy_authority as tax_auth
        from scripts.core.taxonomy_authority import (
            apply_corpus_authority_tiers,
            parse_determiner_tier,
        )

        classify_fn = getattr(tax_auth, "classify_determiner_tier", None)

        def _eval_tier(det_str: str) -> str:
            if callable(classify_fn):
                sig = inspect.signature(classify_fn)
                if "identified_by_raw" in sig.parameters:
                    res = classify_fn(
                        type_status_raw="None",
                        identified_by_raw=det_str,
                        recorded_by_raw="J. Collector",
                        history_raw="",
                        institution_code_raw="NCU",
                        locality_raw="Orange County, Chapel Hill, rich mesic slope",
                        habitat_raw="Deciduous forest",
                    )
                else:
                    res = classify_fn(det_str)
                return res[0] if isinstance(res, tuple) else str(res)
            tier_out, _, _ = parse_determiner_tier(
                type_status_raw="None",
                identified_by_raw=det_str,
                recorded_by_raw="J. Collector",
                history_raw="",
                institution_code_raw="NCU",
                locality_raw="Orange County, Chapel Hill, rich mesic slope",
                habitat_raw="Deciduous forest",
            )
            return tier_out

        specialists = [
            "A.M. Mahoney",
            "R.R. Kowal",
            "D.K. Trock",
            "T.M. Barkley 2002",
            "L.J. Uttal",
            "J.F. Bain",
            "J.L. Golden",
            "C.C. Freeman",
            "G.L. Nesom",
            "E.L. Greene",
            "E.R. Moore",
            "A.S. Weakley",
            "D.E. Boufford",
        ]
        for spec_name in specialists:
            with self.subTest(specialist=spec_name):
                self.assertEqual(_eval_tier(spec_name), "Tier_1_Gold")

        rows = []
        loan_species = [
            "Packera paupercula",
            "Packera plattensis",
            "Packera aurea",
            "Packera obovata",
        ]
        loan_herbaria = ["NCU", "WIS", "MO", "NY"]
        for idx in range(32):
            sp = loan_species[idx % len(loan_species)]
            inst = loan_herbaria[idx % len(loan_herbaria)]
            rows.append(
                {
                    "catalogNumber": f"LOAN_MAHONEY_{idx:02d}",
                    "institutionCode": inst,
                    "scientificName": sp,
                    "species_raw": sp,
                    "type_status": "None",
                    "identifiedBy": "A.M. Mahoney",
                    "recordedBy": f"M. LoanCollector{idx:02d}",
                    "dateIdentified": "",
                    "locality": "Calcareous bluff along river valley",
                    "habitat": "Limestone glade",
                }
            )

        for idx in range(28):
            rows.append(
                {
                    "catalogNumber": f"SINGLE_BATCH_{idx:02d}",
                    "institutionCode": "TENN",
                    "scientificName": "Packera anonyma (Wood) W.A. Weber & A. Love",
                    "species_raw": "Packera anonyma",
                    "type_status": "None",
                    "identifiedBy": "T.M. Barkley",
                    "recordedBy": f"H. BatchCollector{idx:02d}",
                    "dateIdentified": "",
                    "locality": "Knox County, Tennessee, roadside clay bank",
                    "habitat": "Disturbed meadow",
                }
            )

        df_out = apply_corpus_authority_tiers(pd.DataFrame(rows))
        tier_map = dict(zip(df_out["catalogNumber"], df_out["determiner_tier"]))

        for idx in range(32):
            self.assertEqual(tier_map[f"LOAN_MAHONEY_{idx:02d}"], "Tier_1_Gold")
        for idx in range(28):
            self.assertEqual(tier_map[f"SINGLE_BATCH_{idx:02d}"], "Tier_2_Silver")



    def test_usda_plants_symbol_resolution_csv_parsing_and_sqlite_caching(self):
        """Verify USDA PLANTS synonym symbol resolution (Packera dubia -> PATO4), CSV parsing, Excluded filtering, and SQLite caching."""
        from scripts.core.usda_plants_ranges import (
            fetch_or_load_usda_taxon_range,
            parse_usda_distribution_csv,
            reset_usda_network_state,
            resolve_usda_symbols_for_taxon,
        )

        raw_csv = (
            '"Symbol","Scientific Name","Common Name","State","County","FIPS Code","Native Status"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","NC","Orange","37135","L48 (N)"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","NC","Wake","37183","L48 (N)"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","DE","","10","L48 (N)"\n'
            '"PATO4","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","NY","Essex","36031","L48 (X) Excluded"\n'
        )
        parsed = parse_usda_distribution_csv(raw_csv, default_symbol="PATO4")
        self.assertEqual(parsed["counties"], {"north carolina|orange", "north carolina|wake"})
        self.assertEqual(parsed["state_only_states"], {"delaware"})
        self.assertNotIn("new york", parsed["documented_states"])

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            cache_dir = tmp_p / "usda_cache"
            cache_dir.mkdir(parents=True, exist_ok=True)
            db_file = tmp_p / "usda_test.sqlite"

            # 1. Resolution via local plantlst.txt
            plantlst_sample = (
                '"Symbol","Synonym Symbol","Scientific Name with Author","Common Name","Family"\n'
                '"PATO4","","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","Asteraceae"\n'
                '"PATO4","SETO2","Senecio tomentosus Michx.","","Asteraceae"\n'
                '"PASE16","","Packera serpenticola Boufford, Kartesz, S.H. Shi & Renchao Zhou","serpentine ragwort","Asteraceae"\n'
            )
            (cache_dir / "plantlst.txt").write_text(plantlst_sample, encoding="utf-8")
            mock_dubia_sp = MagicMock()
            mock_dubia_sp.name_backbone.return_value = {
                "matchType": "EXACT",
                "status": "ACCEPTED",
                "usageKey": 5387757,
                "canonicalName": "Packera dubia",
                "species": "Packera dubia",
            }
            mock_dubia_sp.name_usage.return_value = {
                "results": [
                    {"canonicalName": "Packera tomentosa", "scientificName": "Packera tomentosa (Michx.) C. Jeffrey", "taxonomicStatus": "SYNONYM"},
                    {"canonicalName": "Senecio tomentosus", "scientificName": "Senecio tomentosus Michx.", "taxonomicStatus": "SYNONYM"},
                ]
            }
            with patch("scripts.core.harvester.species", mock_dubia_sp):
                syms = resolve_usda_symbols_for_taxon(
                    "Packera dubia",
                    cache_dir=cache_dir,
                    allow_network=False,
                )
            self.assertEqual(syms, ["PATO4"])

            # 2. Fetch via mocked HTTP endpoint, verify CSV + SQLite persistence, then reload offline from SQLite
            reset_usda_network_state()
            with patch("scripts.core.harvester.species", mock_dubia_sp), \
                 patch("scripts.core.usda_plants_ranges._http_get_text", return_value=raw_csv) as mock_get:
                r_live = fetch_or_load_usda_taxon_range(
                    canonical_packera="Packera dubia",
                    cache_dir=cache_dir,
                    db_path=db_file,
                    workspace_dir=tmp_p,
                    allow_network=True,
                )
                self.assertEqual(mock_get.call_count, 1)
                self.assertEqual(r_live["symbols"], ["PATO4"])
                self.assertIn("north carolina|orange", r_live["counties"])
                self.assertIn("delaware", r_live["state_only_states"])
                self.assertTrue((cache_dir / "distributions" / "PATO4.csv").is_file())

            # Delete CSV cache and verify offline reload directly from SQLite usda_county_ranges
            (cache_dir / "distributions" / "PATO4.csv").unlink()
            r_sql = fetch_or_load_usda_taxon_range(
                canonical_packera="Packera dubia",
                cache_dir=cache_dir,
                db_path=db_file,
                workspace_dir=tmp_p,
                allow_network=False,
            )
            self.assertEqual(r_sql["symbols"], ["PATO4"])
            self.assertEqual(r_sql["counties"], {"north carolina|orange", "north carolina|wake"})
            self.assertEqual(r_sql["state_only_states"], {"delaware"})

            # 3. Verify USDA PLANTS 3.0 "State and Province,County Name,Distribution Status" CSV schema
            #    and two-step GET /api/PlantProfile -> POST /api/PlantProfile/getDownloadDistributionDocumentation flow
            v3_csv = (
                "Symbol,Scientific Name,Common Name,State and Province,County Name,FIPS Code,Distribution Status\n"
                "PATO4,Packera tomentosa (Michx.) C. Jeffrey,woolly ragwort,Alabama,Autauga County,01001,Native\n"
                "PATO4,Packera tomentosa (Michx.) C. Jeffrey,woolly ragwort,North Carolina,Durham County,37063,Native\n"
                "PATO4,Packera tomentosa (Michx.) C. Jeffrey,woolly ragwort,Michigan,Washtenaw County,26161,Excluded\n"
            )
            parsed_v3 = parse_usda_distribution_csv(v3_csv, default_symbol="PATO4")
            self.assertEqual(parsed_v3["counties"], {"alabama|autauga", "north carolina|durham"})

            db_v3 = tmp_p / "test_usda_v3.sqlite"
            cache_v3 = tmp_p / "usda_cache_v3"
            (cache_v3 / "plantlst.txt").parent.mkdir(parents=True, exist_ok=True)
            (cache_v3 / "plantlst.txt").write_text(plantlst_sample, encoding="utf-8")
            reset_usda_network_state()
            with patch("scripts.core.harvester.species", mock_dubia_sp), patch(
                "scripts.core.usda_plants_ranges._http_get_text",
                return_value='{"Id": 38530, "Symbol": "PATO4"}',
            ) as m_prof, patch(
                "scripts.core.usda_plants_ranges._http_post_text",
                return_value=v3_csv,
            ) as m_post:
                r_v3 = fetch_or_load_usda_taxon_range(
                    canonical_packera="Packera dubia",
                    cache_dir=cache_v3,
                    db_path=db_v3,
                    workspace_dir=tmp_p,
                    allow_network=True,
                )
                self.assertEqual(m_prof.call_count, 1)
                self.assertEqual(m_post.call_count, 1)
                post_payload_sent = m_post.call_args.kwargs.get("payload") or m_post.call_args[0][1]
                self.assertEqual(post_payload_sent["MasterId"], 38530)
                self.assertEqual(post_payload_sent["Symbol"], "PATO4")
                self.assertEqual(r_v3["counties"], {"alabama|autauga", "north carolina|durham"})

    def test_buffered_usda_county_filter_exemptions_and_harvester_integration(self):
        """Verify 1-ring/geodesic county buffer, Type & Gold exemptions, out-of-range misidentification removal, and 26-column CSV invariant."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            usda_dir = tmp_p / "usda_plants"
            (usda_dir / "distributions").mkdir(parents=True, exist_ok=True)
            db_file = tmp_p / "harvester_usda.sqlite"
            out_csv = tmp_p / "curated_vouchers.csv"

            # Synthetic 3-county GeoJSON: Orange (documented), Durham (adjacent 1-ring neighbor ~18 km), and Essex NY (~850 km away)
            counties_geojson = {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"state_name": "north carolina", "county_norm": "orange", "centroid_lon": -79.10, "centroid_lat": 36.05},
                        "geometry": {"type": "Polygon", "coordinates": [[[-79.25, 35.90], [-78.95, 35.90], [-78.95, 36.20], [-79.25, 36.20], [-79.25, 35.90]]]},
                    },
                    {
                        "type": "Feature",
                        "properties": {"state_name": "north carolina", "county_norm": "durham", "centroid_lon": -78.85, "centroid_lat": 36.05},
                        "geometry": {"type": "Polygon", "coordinates": [[[-78.95, 35.90], [-78.70, 35.90], [-78.70, 36.20], [-78.95, 36.20], [-78.95, 35.90]]]},
                    },
                    {
                        "type": "Feature",
                        "properties": {"state_name": "new york", "county_norm": "essex", "centroid_lon": -73.90, "centroid_lat": 44.10},
                        "geometry": {"type": "Polygon", "coordinates": [[[-74.20, 43.80], [-73.50, 43.80], [-73.50, 44.40], [-74.20, 44.40], [-74.20, 43.80]]]},
                    },
                ],
            }
            cnty_geojson_path = tmp_p / "us_counties_test.geojson"
            import json
            cnty_geojson_path.write_text(json.dumps(counties_geojson), encoding="utf-8")

            (usda_dir / "plantlst.txt").write_text(
                '"Symbol","Synonym Symbol","Scientific Name with Author","Common Name","Family"\n'
                '"PATO4","","Packera tomentosa (Michx.) C. Jeffrey","woolly ragwort","Asteraceae"\n',
                encoding="utf-8",
            )
            (usda_dir / "distributions" / "PATO4.csv").write_text(
                '"Symbol","Scientific Name","Common Name","State","County","FIPS Code","Native Status"\n'
                '"PATO4","Packera tomentosa","woolly ragwort","NC","Orange","37135","L48 (N)"\n',
                encoding="utf-8",
            )

            mock_gbif_results = [
                # 1. Documented in-range county (Orange Co., NC) -> RETAIN (in_range)
                {
                    "catalogNumber": "NCU-INRANGE-01",
                    "institutionCode": "NCU",
                    "scientificName": "Packera dubia",
                    "recordedBy": "A. Collector",
                    "recordNumber": "101",
                    "identifiedBy": "Local Botanist",
                    "stateProvince": "North Carolina",
                    "county": "Orange",
                    "decimalLatitude": 36.05,
                    "decimalLongitude": -79.10,
                    "coordinateUncertaintyInMeters": 100.0,
                    "eventDate": "2021-04-20",
                    "year": 2021, "month": 4, "day": 20,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img1.jpg"}],
                },
                # 2. Adjacent county not in USDA checklist ("Swiss-cheese" gap: Durham Co., NC) -> RETAIN (buffered_adjacent)
                {
                    "catalogNumber": "DUKE-ADJ-02",
                    "institutionCode": "DUKE",
                    "scientificName": "Packera dubia",
                    "recordedBy": "B. Collector",
                    "recordNumber": "102",
                    "identifiedBy": "Local Botanist",
                    "stateProvince": "North Carolina",
                    "county": "Durham",
                    "decimalLatitude": 36.05,
                    "decimalLongitude": -78.85,
                    "coordinateUncertaintyInMeters": 100.0,
                    "eventDate": "2021-04-22",
                    "year": 2021, "month": 4, "day": 22,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img2.jpg"}],
                },
                # 3. Far out-of-range misidentification (Essex Co., NY ~850 km away, Bronze/Silver) -> DROP
                {
                    "catalogNumber": "NY-MISID-03",
                    "institutionCode": "NY",
                    "scientificName": "Packera dubia",
                    "recordedBy": "C. Collector",
                    "recordNumber": "103",
                    "identifiedBy": "Student Determiner",
                    "stateProvince": "New York",
                    "county": "Essex",
                    "decimalLatitude": 44.10,
                    "decimalLongitude": -73.90,
                    "coordinateUncertaintyInMeters": 250.0,
                    "eventDate": "2020-06-15",
                    "year": 2020, "month": 6, "day": 15,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img3.jpg"}],
                },
                # 4. Out-of-range Tier_1_Gold specialist determination (A.S. Weakley) -> EXEMPT RETAIN
                {
                    "catalogNumber": "NCU-GOLD-EXT-04",
                    "institutionCode": "NCU",
                    "scientificName": "Packera dubia",
                    "recordedBy": "D. Collector",
                    "recordNumber": "104",
                    "identifiedBy": "A.S. Weakley",
                    "stateProvince": "New York",
                    "county": "Essex",
                    "decimalLatitude": 44.12,
                    "decimalLongitude": -73.92,
                    "coordinateUncertaintyInMeters": 100.0,
                    "eventDate": "2022-05-30",
                    "year": 2022, "month": 5, "day": 30,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img4.jpg"}],
                },
                # 5. Out-of-range Holotype -> EXEMPT RETAIN
                {
                    "catalogNumber": "US-TYPE-05",
                    "institutionCode": "US",
                    "scientificName": "Packera dubia",
                    "typeStatus": "Holotype",
                    "recordedBy": "A. Michaux",
                    "recordNumber": "1",
                    "identifiedBy": "Unknown",
                    "stateProvince": "New York",
                    "county": "Essex",
                    "decimalLatitude": 44.05,
                    "decimalLongitude": -73.88,
                    "coordinateUncertaintyInMeters": 500.0,
                    "eventDate": "1803-05-10",
                    "year": 1803, "month": 5, "day": 10,
                    "media": [{"type": "StillImage", "identifier": "https://example.org/img5.jpg"}],
                },
            ]

            mock_occ = MagicMock()
            mock_occ.search.side_effect = lambda **kw: (
                {"results": [], "endOfRecords": True}
                if "typeStatus" in kw
                else {"results": mock_gbif_results, "endOfRecords": True}
            )

            mock_dubia_sp2 = MagicMock()
            mock_dubia_sp2.name_backbone.return_value = {
                "matchType": "EXACT",
                "status": "ACCEPTED",
                "usageKey": 5387757,
                "canonicalName": "Packera dubia",
                "species": "Packera dubia",
            }
            mock_dubia_sp2.name_usage.return_value = {
                "results": [
                    {"canonicalName": "Packera tomentosa", "scientificName": "Packera tomentosa (Michx.) C. Jeffrey", "taxonomicStatus": "SYNONYM"},
                ]
            }
            with patch("scripts.core.harvester.occ", mock_occ), patch("scripts.core.harvester.species", mock_dubia_sp2):
                harv = VoucherHarvester(
                    taxa=["Packera dubia"],
                    expand_synonyms=False,
                    output_csv=out_csv,
                    raw_dir=tmp_p / "raw",
                    workspace_dir=tmp_p,
                    db_path=db_file,
                    usda_county_filter=True,
                    usda_county_buffer_km=50.0,
                    usda_exempt_gold=True,
                    usda_cache_dir=usda_dir,
                    usda_county_geojson_path=cnty_geojson_path,
                )
                df_res = harv.harvest()
                harv.export(df_res)

            retained_cats = set(df_res["catalogNumber"].tolist())
            self.assertEqual(
                retained_cats,
                {"NCU-INRANGE-01", "DUKE-ADJ-02", "NCU-GOLD-EXT-04", "US-TYPE-05"},
            )
            self.assertNotIn("NY-MISID-03", retained_cats)
            self.assertEqual(harv.last_usda_filter_stats["dropped_out_of_range"], 1)
            self.assertEqual(harv.last_usda_filter_stats["buffered_adjacent"], 1)
            self.assertEqual(harv.last_usda_filter_stats["exempt_retained"], 2)

            # Verify SQLite-only USDA diagnostic columns never leak into curated_vouchers.csv
            csv_df = pd.read_csv(out_csv)
            from scripts.core.harvester import EXPORT_COLUMNS
            self.assertEqual(list(csv_df.columns[: len(EXPORT_COLUMNS)]), EXPORT_COLUMNS)
            self.assertNotIn("usda_symbol", csv_df.columns)
            self.assertNotIn("usda_county_status", csv_df.columns)

            import sqlite3
            conn = sqlite3.connect(str(db_file))
            try:
                rows_sql = dict(
                    conn.execute(
                        "SELECT catalogNumber, usda_county_status FROM vouchers;"
                    ).fetchall()
                )
                cnt_ranges = conn.execute(
                    "SELECT COUNT(*) FROM usda_county_ranges WHERE canonical_packera = 'Packera dubia';"
                ).fetchone()[0]
            finally:
                conn.close()
            self.assertEqual(rows_sql["NCU-INRANGE-01"], "in_range")
            self.assertEqual(rows_sql["DUKE-ADJ-02"], "buffered_adjacent")
            self.assertEqual(rows_sql["NCU-GOLD-EXT-04"], "out_of_range")
            self.assertEqual(rows_sql["US-TYPE-05"], "out_of_range")
            self.assertEqual(int(cnt_ranges), 1)



class TestStep1SystematicsAndDownstreamAuthorityAudit(unittest.TestCase):
    """Comprehensive regression tests for the Step 1 botanical systematics & downstream authority-tier review."""

    def test_authority_engine_checklist_blocking_surname_qualification_and_self_id_guard(self):
        """Verify flora/checklist blocking, initial-qualified common surnames, regional floristician Silver routing, and Layer 3 self-ID exclusion."""
        from scripts.core.taxonomy_authority import (
            parse_determiner_tier,
            discover_corpus_authorities,
        )

        # 1. Flora / checklist strings must NOT promote to Tier_1_Gold or Tier_2_Silver
        for flora_str in ("Weakley 2015", "Flora of Virginia (Weakley)", "FNA", "Kartesz", "BONAP", "sensu T.M. Barkley"):
            tier, _, _ = parse_determiner_tier(
                type_status_raw=None,
                identified_by_raw=flora_str,
                recorded_by_raw="Field Collector",
                history_raw=None,
                institution_code_raw="NCU",
                locality_raw="Rich calcareous seepage slope over dolostone along rocky bluff, elev. 420 m",
                habitat_raw="Calcareous seepage",
            )
            self.assertEqual(tier, "Tier_3_Bronze", f"Flora citation '{flora_str}' leaked above Bronze: {tier}")

        # 2. Unrelated collectors/annotators sharing common surnames (Fuller, Turner, Freeman, Weber) do NOT get Tier_1_Gold
        for non_spec in ("S. Fuller", "M. Turner", "J. Freeman", "K. Weber"):
            tier, _, _ = parse_determiner_tier(
                type_status_raw=None,
                identified_by_raw=non_spec,
                recorded_by_raw=non_spec,
                history_raw=None,
                institution_code_raw="UNKNOWN_HERB",
                locality_raw="Roadside",
                habitat_raw=None,
            )
            self.assertEqual(tier, "Tier_3_Bronze", f"Non-specialist '{non_spec}' falsely promoted to {tier}")

        # While initial-qualified core Packera specialists DO get Tier_1_Gold
        for real_spec in ("T.M. Barkley", "Debra K. Trock", "R.R. Kowal", "J.F. Bain", "B.L. Turner", "C.C. Freeman", "B. Fuller"):
            tier, _, _ = parse_determiner_tier(
                type_status_raw=None,
                identified_by_raw=real_spec,
                recorded_by_raw="Field Collector",
                history_raw=None,
                institution_code_raw="UNKNOWN_HERB",
                locality_raw="Roadside",
                habitat_raw=None,
            )
            self.assertEqual(tier, "Tier_1_Gold", f"Core specialist '{real_spec}' failed Tier_1_Gold: {tier}")

        # 3. Verbatim identification history does NOT override a later explicit generalist identifiedBy
        tier_reannot, _, _ = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw="Undergraduate Botany Class",
            recorded_by_raw="Field Collector",
            history_raw="Det. T.M. Barkley, 1978",
            institution_code_raw="NCU",
            locality_raw="Calcareous cliff ledge along river",
            habitat_raw="Limestone bluff",
        )
        self.assertEqual(tier_reannot, "Tier_3_Bronze")

        # 4. Regional floristician field collection with blank identifiedBy routes to Tier_2_Silver (not blind Tier_1_Gold)
        tier_reg_collector, _, _ = parse_determiner_tier(
            type_status_raw=None,
            identified_by_raw="",
            recorded_by_raw="Julian Campbell",
            history_raw=None,
            institution_code_raw="EKY",
            locality_raw="Limestone glade",
            habitat_raw=None,
        )
        self.assertEqual(tier_reg_collector, "Tier_2_Silver")

        # 5. Layer 3 self-identification across distributed duplicates (identifiedBy == recordedBy) must NOT promote collector
        dup_rows = [
            {
                "catalogNumber": f"DUP_{inst}_{i}",
                "institutionCode": inst,
                "scientificName": "Packera aurea" if i % 2 == 0 else "Packera anonyma",
                "identifiedBy": "H.E. Generalcollector",
                "recordedBy": "H.E. Generalcollector",
                "recordNumber": str(100 + i),
                "type_status": "None",
                "eventDate": "2018-05-10",
            }
            for i, inst in enumerate(["NCU", "NY", "MO", "US", "GH", "DUKE"])
        ]
        dyn_specs, _ = discover_corpus_authorities(pd.DataFrame(dup_rows))
        self.assertNotIn("generalcollector", dyn_specs)

    def test_no_packera_to_packera_backbone_lumping_and_infraspecific_persistence(self):
        """Verify GBIF Backbone only auto-maps historical generic transfers (Senecio/Cineraria -> Packera) and saves infraspecific_taxon."""
        import sqlite3
        from scripts.core import harvester as ch
        from scripts.core.harvester import VoucherHarvester, export_curated_table

        ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop("packera pseudaurea", None)
        ch._GBIF_SYNONYM_CACHE.pop("packera pseudaurea", None)
        ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop("senecio aureus", None)
        ch._GBIF_SYNONYM_CACHE.pop("senecio aureus", None)

        mock_sp = MagicMock()
        def _fake_bb(**kwargs):
            nm = str(kwargs.get("name", "")).lower()
            if "pseudaurea" in nm:
                return {"canonicalName": "Packera pseudaurea", "species": "Packera aurea", "accepted": "Packera aurea"}
            if "senecio aureus" in nm:
                return {"canonicalName": "Senecio aureus", "species": "Packera aurea", "accepted": "Packera aurea"}
            return {}
        mock_sp.name_backbone.side_effect = _fake_bb

        with patch("scripts.core.harvester.species", mock_sp):
            # Heterotypic Packera -> Packera backbone lump must NOT occur automatically
            self.assertEqual(ch.standardize_packera_taxon("Packera pseudaurea var. semicordata"), "Packera pseudaurea")
            # Historical Senecio -> Packera generic transfer DOES resolve cleanly
            self.assertEqual(ch.standardize_packera_taxon("Senecio aureus L."), "Packera aurea")

        def _fake_occ_var(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            return {
                "results": [
                    {
                        "key": 55101,
                        "scientificName": "Packera pseudaurea (Rydb.) W.A.Weber & Á.Löve var. semicordata (Mack. & Bush) Trock & T.M.Barkley",
                        "infraspecificEpithet": "semicordata",
                        "taxonRank": "VARIETY",
                        "stateProvince": "Wisconsin",
                        "decimalLatitude": 43.1,
                        "decimalLongitude": -89.4,
                        "year": 1995,
                        "month": 5,
                        "recordedBy": "T.M. Barkley",
                        "identifiedBy": "Debra K. Trock",
                        "recordNumber": "2022-014",
                        "media": [{"identifier": "https://example.com/semicordata.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "12345",
                        "institutionCode": "WIS",
                    },
                    {
                        "key": 55102,
                        "scientificName": "Packera pseudaurea var. semicordata",
                        "stateProvince": "Minnesota",
                        "decimalLatitude": 44.9,
                        "decimalLongitude": -93.2,
                        "eventDate": "1988-01-01",
                        "year": 1988,
                        "month": 1,
                        "day": 1,
                        "recordedBy": "T.M. Barkley",
                        "identifiedBy": "Debra K. Trock",
                        "recordNumber": "2022-015",
                        "media": [{"identifier": "https://example.com/semicordata2.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "12345",
                        "institutionCode": "MIN",
                    },
                ],
                "count": 2,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_var):
            harv = VoucherHarvester(taxa=["Packera pseudaurea"], expand_synonyms=False)
            df_v = harv.harvest()

        self.assertEqual(len(df_v), 2)
        # Bare numeric catalogNumbers 12345 from WIS and MIN are namespaced instead of colliding!
        self.assertEqual(set(df_v["catalogNumber"]), {"WIS_12345", "MIN_12345"})
        self.assertEqual(set(df_v["infraspecific_taxon"]), {"var. semicordata"})
        # Missing coordinateUncertaintyInMeters stays NaN/None instead of fabricated 5000.0 m
        self.assertTrue(df_v["coordinateUncertainty"].isna().all())
        # Partial year+month (1995-05) retained, and 1988-01-01 placeholder masked to year-only (month/day NaN)
        min_row = df_v[df_v["catalogNumber"] == "MIN_12345"].iloc[0]
        self.assertTrue(pd.isna(min_row["month"]) and pd.isna(min_row["day"]))

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            db_p = tmp_p / "var_test.sqlite"
            csv_p = tmp_p / "curated_vouchers.csv"
            export_curated_table(df_v, csv_p, db_path=db_p)
            conn = sqlite3.connect(str(db_p))
            try:
                rows_infra = dict(conn.execute("SELECT catalogNumber, infraspecific_taxon FROM vouchers").fetchall())
            finally:
                conn.close()
            self.assertEqual(rows_infra["WIS_12345"], "var. semicordata")
            self.assertEqual(rows_infra["MIN_12345"], "var. semicordata")

    def test_coordinate_cleaner_cultivated_filter_anonyma_usda_chesapeake_ecoregion_and_media_rescue(self):
        """Verify CoordinateCleaner rejection, cultivated filtering, Packera anonyma USDA symbol, Chesapeake VA ecoregion, and secondary duplicate media rescue."""
        from scripts.core.harvester import VoucherHarvester
        from scripts.core.usda_plants_ranges import resolve_usda_symbols_for_taxon
        from scripts.core.epa_ecoregions import _HISTORIC_COUNTY_CONSOLIDATIONS, _normalize_us_county_name, infer_regional_group

        # 1. Packera anonyma dynamic USDA symbol resolution (resolves via GBIF synonym Senecio smallii + algorithmic ICN generic transfer)
        with tempfile.TemporaryDirectory() as tmp_usda:
            p_lst = Path(tmp_usda) / "plantlst.txt"
            p_lst.write_text(
                '"Symbol","Synonym Symbol","Scientific Name with Author","Common Name","Family"\n'
                '"PAAN6","","Senecio smallii Britton","Small\'s ragwort","Asteraceae"\n'
                '"PACR11","","Senecio crawfordii (Britton) G.W. Douglas & Ruyle-Douglas","Crawford\'s ragwort","Asteraceae"\n',
                encoding="utf-8",
            )
            syms_craw = resolve_usda_symbols_for_taxon("Packera crawfordii", cache_dir=Path(tmp_usda), allow_network=False)
            self.assertEqual(syms_craw, ["PACR11"])
            mock_anon_sp = MagicMock()
            mock_anon_sp.name_backbone.return_value = {
                "matchType": "EXACT",
                "status": "ACCEPTED",
                "usageKey": 3109240,
                "canonicalName": "Packera anonyma",
                "species": "Packera anonyma",
            }
            mock_anon_sp.name_usage.return_value = {
                "results": [
                    {"canonicalName": "Senecio smallii", "scientificName": "Senecio smallii Britton", "taxonomicStatus": "SYNONYM"},
                ]
            }
            with patch("scripts.core.harvester.species", mock_anon_sp):
                syms_anon = resolve_usda_symbols_for_taxon("Packera anonyma", cache_dir=Path(tmp_usda), allow_network=False)
            self.assertEqual(syms_anon, ["PAAN6"])

        # 2. Chesapeake VA county name consolidation fix (not overwritten to 'norfolk'; South Norfolk -> 'chesapeake')
        self.assertEqual(_normalize_us_county_name("Chesapeake", "virginia"), "chesapeake")
        self.assertEqual(_normalize_us_county_name("South Norfolk", "virginia"), "chesapeake")
        # Eastern Canadian CEC Level III ecoregion assignment (even when CONUS-only us_eco_l3.geojson is cached on disk!)
        import json
        with tempfile.TemporaryDirectory() as tmp_eco_d:
            conus_only_geojson = Path(tmp_eco_d) / "us_eco_l3.geojson"
            conus_only_geojson.write_text(json.dumps({
                "type": "FeatureCollection",
                "features": [{
                    "type": "Feature",
                    "properties": {"US_L3CODE": "45", "US_L3NAME": "Piedmont", "ecoregion_label": "Piedmont"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [[[-80.0, 35.0], [-78.5, 35.0], [-78.5, 36.5], [-80.0, 36.5], [-80.0, 35.0]]],
                    },
                }],
            }))
            with patch("scripts.core.epa_ecoregions.DEFAULT_EPA_L3_GEOJSON", conus_only_geojson),                  patch("scripts.core.epa_ecoregions._ensure_epa_level3_cache_files", return_value=False):
                eco_ottawa = infer_regional_group(lat=45.42, lon=-75.69, state_province="Ontario")
                eco_erie = infer_regional_group(lat=42.98, lon=-81.25, state_province="Ontario")
                eco_nc = infer_regional_group(lat=35.91, lon=-79.05, state_province="North Carolina")
            self.assertEqual(eco_ottawa, "Eastern_Great_Lakes_Lowlands")
            self.assertEqual(eco_erie, "Lake_Erie_Lowland")
            self.assertEqual(eco_nc, "Piedmont")

        # 3. CoordinateCleaner + Cultivated occurrence filtering during harvest()
        def _fake_occ_hygiene(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            return {
                "results": [
                    # Valid NC sheet
                    {
                        "key": 66001,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.05,
                        "year": 2015,
                        "recordedBy": "R. Collector",
                        "recordNumber": "101",
                        "media": [{"identifier": "https://example.com/valid.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "NCU_VALID_1",
                        "institutionCode": "NCU",
                    },
                    # (0, 0) Gulf of Guinea error -> rejected
                    {
                        "key": 66002,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 0.0,
                        "decimalLongitude": 0.0,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/zero.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_ZERO_2",
                        "institutionCode": "NCU",
                    },
                    # Sign-inverted positive longitude in NA -> rejected
                    {
                        "key": 66003,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": 79.05,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/poslon.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_POSLON_3",
                        "institutionCode": "NCU",
                    },
                    # State-coordinate mismatch (says North Carolina, coordinates in Utah) -> rejected
                    {
                        "key": 66004,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 39.5,
                        "decimalLongitude": -111.5,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/utah.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_STATEMISMATCH_4",
                        "institutionCode": "NCU",
                    },
                    # Cultivated botanical garden sheet -> completely excluded
                    {
                        "key": 66005,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.9,
                        "decimalLongitude": -79.03,
                        "locality": "Cultivated in Sarah P. Duke Gardens ornamental bed",
                        "year": 2018,
                        "media": [{"identifier": "https://example.com/cult.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_CULT_5",
                        "institutionCode": "DUKE",
                    },
                    # (0, 0) error AND missing stateProvince -> completely excluded (cannot place biogeographically)
                    {
                        "key": 66006,
                        "scientificName": "Packera anonyma",
                        "stateProvince": "",
                        "decimalLatitude": 0.0,
                        "decimalLongitude": 0.0,
                        "year": 2015,
                        "media": [{"identifier": "https://example.com/zeronostate.jpg", "type": "StillImage", "format": "image/jpeg"}],
                        "catalogNumber": "BAD_ZERO_NOSTATE_6",
                        "institutionCode": "NCU",
                    },
                ],
                "count": 6,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_fake_occ_hygiene):
            harv_h = VoucherHarvester(taxa=["Packera anonyma"], expand_synonyms=False)
            df_h = harv_h.harvest().set_index("catalogNumber")
        self.assertNotIn("BAD_CULT_5", df_h.index)
        self.assertNotIn("BAD_ZERO_NOSTATE_6", df_h.index)
        self.assertAlmostEqual(float(df_h.loc["NCU_VALID_1", "decimalLatitude"]), 35.9, places=4)
        self.assertAlmostEqual(float(df_h.loc["NCU_VALID_1", "decimalLongitude"]), -79.05, places=4)
        self.assertEqual(df_h.loc["NCU_VALID_1", "georef_qc_status"], "valid")
        for bad_cat, exp_qc in [
            ("BAD_ZERO_2", "zero_coord_nullified"),
            ("BAD_POSLON_3", "hemisphere_sign_mismatch_nullified"),
            ("BAD_STATEMISMATCH_4", "state_polygon_mismatch_nullified"),
        ]:
            self.assertIn(bad_cat, df_h.index)
            self.assertTrue(pd.isna(df_h.loc[bad_cat, "decimalLatitude"]))
            self.assertTrue(pd.isna(df_h.loc[bad_cat, "decimalLongitude"]))
            self.assertEqual(df_h.loc[bad_cat, "georef_qc_status"], exp_qc)

        # 4. Secondary duplicate media rescue when primary URL fails download/QC
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            good_src = tmp_p / "secondary_good.jpg"
            Image.new("RGB", (3200, 4000), color=(110, 130, 105)).save(good_src, format="JPEG", quality=95)

            harv_m = VoucherHarvester(
                taxa=["Packera anonyma"],
                workspace_dir=tmp_p,
                min_megapixels=4.0,
                min_file_size_kb=5,
                expand_synonyms=False,
            )
            harv_m.workspace_dir = tmp_p
            harv_m._secondary_duplicate_pool = pd.DataFrame([
                {
                    "catalogNumber": "MO_SEC_01",
                    "institutionCode": "MO",
                    "scientificName": "Packera anonyma",
                    "recordedBy": "T.M. Barkley",
                    "recordNumber": "500",
                    "eventDate": "1980-05-10",
                    "determiner_tier": "Tier_1_Gold",
                    "exsiccatae_key": "barkley_500_1980-05-10",
                    "image_path": "raw_vouchers/MO_SEC_01.webp",
                    "_image_url": "https://example.com/secondary_good.jpg",
                }
            ])
            df_prim = pd.DataFrame([
                {
                    "catalogNumber": "NCU_PRIM_01",
                    "institutionCode": "NCU",
                    "scientificName": "Packera anonyma",
                    "recordedBy": "T.M. Barkley",
                    "recordNumber": "500",
                    "eventDate": "1980-05-10",
                    "determiner_tier": "Tier_1_Gold",
                    "exsiccatae_key": "barkley_500_1980-05-10",
                    "image_path": "raw_vouchers/NCU_PRIM_01.webp",
                    "_image_url": "https://example.com/broken_primary.jpg",
                }
            ])

            async def _fake_dl_all(records_to_download, **kwargs):
                is_rescue = any("secondary_good" in str(u) for u, _ in records_to_download)
                if not is_rescue:
                    return {"success": 0, "skipped": 0, "failed": 1}
                dest = Path(records_to_download[0][1])
                dest.parent.mkdir(parents=True, exist_ok=True)
                Image.open(good_src).save(dest, format="WEBP", quality=92)
                return {"success": 1, "skipped": 0, "failed": 0}

            with patch("scripts.core.harvester.download_all_voucher_images", side_effect=_fake_dl_all):
                df_rescued, rescue_stats = harv_m.download_and_validate_media(df_prim)
            self.assertEqual(len(df_rescued), 1)
            self.assertEqual(df_rescued.iloc[0]["catalogNumber"], "MO_SEC_01")
            self.assertTrue(bool(df_rescued.iloc[0]["is_primary_duplicate"]))
            self.assertEqual(rescue_stats["quality_rejected"], 1)
            self.assertTrue((tmp_p / "raw_vouchers" / "MO_SEC_01.webp").exists())

    def test_zero_hardcoded_synonyms_and_dynamic_packera_tomentosa_dubia_unification(self):
        """Verify zero hardcoded taxon names in source and 100% dynamic unification of Packera tomentosa into Packera dubia."""
        import subprocess
        import main as root_main
        from scripts.core import harvester as ch
        from scripts.core.harvester import VoucherHarvester

        project_root = Path(__file__).resolve().parents[2]
        harv_src = (project_root / "scripts/core/harvester.py").read_text(encoding="utf-8")
        usda_src = (project_root / "scripts/core/usda_plants_ranges.py").read_text(encoding="utf-8")
        aw_r_src = (project_root / "scripts/morphometrics/authority_weighting_utils.R").read_text(encoding="utf-8")
        self.assertNotIn('"packera tomentosa"', harv_src.lower())
        self.assertNotIn('"senecio tomentosus"', harv_src.lower())
        self.assertNotIn("fallback_map", usda_src)
        self.assertNotIn('"packera dubia":', usda_src.lower())
        self.assertNotIn('grepl("dubia|tomentosa', aw_r_src)

        for k in ("packera dubia", "packera tomentosa", "senecio tomentosus", "cineraria dubia", "cineraria tomentosa", "senecio dubius"):
            ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop(k, None)
            ch._GBIF_SYNONYM_CACHE.pop(k, None)
            ch._GBIF_BACKBONE_SYNONYM_FLAG.pop(k, None)

        mock_gbif_sp = MagicMock()
        def _gbif_bb(**kwargs):
            nm = str(kwargs.get("name") or kwargs.get("scientificName") or "").strip().lower()
            if "packera dubia" in nm:
                return {
                    "matchType": "EXACT",
                    "status": "ACCEPTED",
                    "usageKey": 5387757,
                    "canonicalName": "Packera dubia",
                    "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                }
            if "packera tomentosa" in nm:
                return {
                    "matchType": "EXACT",
                    "status": "SYNONYM",
                    "usageKey": 3109265,
                    "acceptedUsageKey": 5387757,
                    "canonicalName": "Packera tomentosa",
                    "scientificName": "Packera tomentosa (Michx.) C. Jeffrey",
                    "accepted": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                }
            return {"matchType": "NONE"}

        def _gbif_nu(**kwargs):
            if kwargs.get("data") == "synonyms":
                return {
                    "results": [
                        {"canonicalName": "Packera tomentosa", "scientificName": "Packera tomentosa (Michx.) C. Jeffrey", "taxonomicStatus": "SYNONYM"},
                        {"canonicalName": "Senecio tomentosus", "scientificName": "Senecio tomentosus Michx.", "taxonomicStatus": "SYNONYM"},
                        {"canonicalName": "Cineraria dubia", "scientificName": "Cineraria dubia Spreng.", "taxonomicStatus": "HOMOTYPIC_SYNONYM"},
                        {"canonicalName": "Senecio dubius", "scientificName": "Senecio dubius (Spreng.)", "taxonomicStatus": "HOMOTYPIC_SYNONYM"},
                        {"canonicalName": "Packera aurea", "scientificName": "Packera aurea auct. non (L.)", "taxonomicStatus": "MISAPPLIED"},
                    ]
                }
            return {}

        mock_gbif_sp.name_backbone.side_effect = _gbif_bb
        mock_gbif_sp.name_usage.side_effect = _gbif_nu

        def _gbif_occ(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            sci_q = str(kwargs.get("scientificName") or "")
            recs = []
            if sci_q == "Packera dubia":
                recs.append({
                    "key": 90101,
                    "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "taxonomicStatus": "ACCEPTED",
                    "species": "Packera dubia",
                    "stateProvince": "North Carolina",
                    "decimalLatitude": 35.91,
                    "decimalLongitude": -79.05,
                    "year": 2021, "month": 4, "day": 18,
                    "recordedBy": "B. Fuller", "recordNumber": "401",
                    "identifiedBy": "A.S. Weakley",
                    "media": [{"identifier": "https://example.org/dubia1.jpg", "type": "StillImage"}],
                    "catalogNumber": "NCU_DUBIA_01",
                    "institutionCode": "NCU",
                })
            elif sci_q == "Packera tomentosa":
                recs.append({
                    "key": 90102,
                    "scientificName": "Packera tomentosa (Michx.) C. Jeffrey",
                    "taxonomicStatus": "SYNONYM",
                    "acceptedScientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                    "stateProvince": "South Carolina",
                    "decimalLatitude": 34.00,
                    "decimalLongitude": -81.03,
                    "year": 1998, "month": 4, "day": 12,
                    "recordedBy": "J.B. Nelson", "recordNumber": "19200",
                    "identifiedBy": "T.M. Barkley",
                    "media": [{"identifier": "https://example.org/tomentosa1.jpg", "type": "StillImage"}],
                    "catalogNumber": "USCH_TOM_02",
                    "institutionCode": "USCH",
                })
            elif sci_q == "Senecio tomentosus":
                recs.append({
                    "key": 90103,
                    "scientificName": "Senecio tomentosus Michx.",
                    "taxonomicStatus": "SYNONYM",
                    "acceptedScientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                    "species": "Packera dubia",
                    "stateProvince": "Virginia",
                    "decimalLatitude": 36.85,
                    "decimalLongitude": -76.29,
                    "year": 1968, "month": 5, "day": 3,
                    "recordedBy": "A.E. Radford", "recordNumber": "45012",
                    "identifiedBy": "A.E. Radford",
                    "media": [{"identifier": "https://example.org/sentom1.jpg", "type": "StillImage"}],
                    "catalogNumber": "NCU_SENTOM_03",
                    "institutionCode": "NCU",
                })
            return {"results": recs, "count": len(recs)}

        # Case A: Querying Packera dubia expands and combines Packera tomentosa + Senecio tomentosus into Packera dubia
        with patch("scripts.core.harvester.species", mock_gbif_sp), \
             patch("scripts.core.harvester.occ.search", side_effect=_gbif_occ):
            harv_dubia = VoucherHarvester(taxa=["Packera dubia"], expand_synonyms=True)
            self.assertIn("Packera tomentosa", harv_dubia.taxa)
            self.assertIn("Senecio tomentosus", harv_dubia.taxa)
            self.assertNotIn("Packera aurea", harv_dubia.taxa)  # MISAPPLIED rejected!
            df_combined = harv_dubia.harvest()

        self.assertEqual(set(df_combined["catalogNumber"]), {"NCU_DUBIA_01", "USCH_TOM_02", "NCU_SENTOM_03"})
        self.assertEqual(set(df_combined["species_standardized"]), {"Packera dubia"})
        tom_row = df_combined[df_combined["catalogNumber"] == "USCH_TOM_02"].iloc[0]
        self.assertEqual(tom_row["scientificName"], "Packera dubia")
        self.assertEqual(tom_row["species_raw"], "Packera tomentosa (Michx.) C. Jeffrey")

        # Case B: Querying unaccepted Packera tomentosa directly resolves canonical to Packera dubia and combines both
        for k in ("packera dubia", "packera tomentosa", "senecio tomentosus", "cineraria dubia", "senecio dubius"):
            ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop(k, None)
            ch._GBIF_SYNONYM_CACHE.pop(k, None)
            ch._GBIF_BACKBONE_SYNONYM_FLAG.pop(k, None)

        with patch("scripts.core.harvester.species", mock_gbif_sp), \
             patch("scripts.core.harvester.occ.search", side_effect=_gbif_occ):
            harv_tom = VoucherHarvester(taxa=["Packera tomentosa"], expand_synonyms=True)
            df_from_tom = harv_tom.harvest()
            self.assertEqual(set(df_from_tom["catalogNumber"]), {"NCU_DUBIA_01", "USCH_TOM_02", "NCU_SENTOM_03"})
            self.assertEqual(set(df_from_tom["species_standardized"]), {"Packera dubia"})

        # Case C: Occurrence-level SYNONYM status combines Packera tomentosa into Packera dubia even when expand_synonyms=False
        for k in ("packera dubia", "packera tomentosa", "senecio tomentosus"):
            ch._DYNAMIC_SYNONYM_TO_CANONICAL.pop(k, None)
            ch._GBIF_SYNONYM_CACHE.pop(k, None)
            ch._GBIF_BACKBONE_SYNONYM_FLAG.pop(k, None)

        def _mixed_occ(**kwargs):
            if "typeStatus" in kwargs or kwargs.get("country") == "CA":
                return {"results": [], "count": 0}
            return {
                "results": [
                    {
                        "key": 90201,
                        "scientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                        "taxonomicStatus": "ACCEPTED",
                        "species": "Packera dubia",
                        "stateProvince": "North Carolina",
                        "decimalLatitude": 35.91, "decimalLongitude": -79.05,
                        "year": 2021, "recordedBy": "B. Fuller", "recordNumber": "401",
                        "media": [{"identifier": "https://example.org/d.jpg", "type": "StillImage"}],
                        "catalogNumber": "NCU_D_01", "institutionCode": "NCU",
                    },
                    {
                        "key": 90202,
                        "scientificName": "Packera tomentosa (Michx.) C. Jeffrey",
                        "taxonomicStatus": "SYNONYM",
                        "acceptedScientificName": "Packera dubia (Spreng.) Trock & Mabb.",
                        "species": "Packera dubia",
                        "stateProvince": "South Carolina",
                        "decimalLatitude": 34.00, "decimalLongitude": -81.03,
                        "year": 1998, "recordedBy": "J.B. Nelson", "recordNumber": "19200",
                        "media": [{"identifier": "https://example.org/t.jpg", "type": "StillImage"}],
                        "catalogNumber": "USCH_T_02", "institutionCode": "USCH",
                    },
                ],
                "count": 2,
            }

        with patch("scripts.core.harvester.occ.search", side_effect=_mixed_occ):
            harv_occ_syn = VoucherHarvester(taxa=["Packera dubia"], expand_synonyms=False)
            df_occ_syn = harv_occ_syn.harvest()
            self.assertEqual(set(df_occ_syn["species_standardized"]), {"Packera dubia"})

        # Verify main.py forwards the 4 Step 1 flags when explicitly specified
        parser = root_main.build_parser()
        args_flags = parser.parse_args([
            "harvest",
            "--no-harvest-full-continental-range",
            "--spatial-thin-km", "12.5",
            "--max-per-state", "40",
            "--max-per-region", "80",
        ])
        cfg_stub = {
            "paths": {"raw_vouchers_dir": "/tmp/r", "curated_vouchers_csv": "/tmp/c.csv"},
            "harvesting": {
                "max_records_per_taxon": 100, "max_uncertainty_meters": 5000,
                "min_megapixels": 8.0, "min_file_size_kb": 500, "download_concurrency": 5,
            },
        }
        with patch("subprocess.run") as m_sub, patch("main.verify_file_exists"):
            root_main.run_harvest(args_flags, cfg_stub)
            sent_cmd = m_sub.call_args[0][0]
            self.assertIn("--no-harvest-full-continental-range", sent_cmd)
            self.assertIn("--spatial-thin-km", sent_cmd)
            self.assertIn("12.5", sent_cmd)
            self.assertIn("--max-per-state", sent_cmd)
            self.assertIn("40", sent_cmd)
            self.assertIn("--max-per-region", sent_cmd)
            self.assertIn("80", sent_cmd)

        # Verify R .load_authority_weight_config reads nested thresholds$authority_weights from YAML
        self.assertIn("cfg$thresholds$authority_weights", aw_r_src)
        self.assertIn("dyn_equiv", aw_r_src)
        if rscript_available():
            with tempfile.TemporaryDirectory() as tmp_r_d:
                custom_yaml = Path(tmp_r_d) / "cfg.yaml"
                custom_yaml.write_text(
                    "thresholds:\n  authority_weights:\n    tier_0_type: 3.5\n    silver_max_pool_share: 0.35\n",
                    encoding="utf-8",
                )
                r_out = run_warm_r_cmd(
                    [
                        "Rscript", "-e",
                        f"source('scripts/morphometrics/authority_weighting_utils.R'); "
                        f"rm(list=ls(envir=.AUTH_WEIGHT_CFG_ENV), envir=.AUTH_WEIGHT_CFG_ENV); "
                        f"on.exit(rm(list=ls(envir=.AUTH_WEIGHT_CFG_ENV), envir=.AUTH_WEIGHT_CFG_ENV), add=TRUE); "
                        f"w <- .load_authority_weight_config({repr(str(custom_yaml))}); "
                        f"stopifnot(abs(w$tier_0_type - 3.5) < 1e-6, abs(w$silver_max_pool_share - 0.35) < 1e-6); "
                        f"cat('R_CFG_OK\\n')",
                    ],
                    stage_key="step03",
                    cwd=project_root,
                )
                self.assertEqual(r_out.returncode, 0, r_out.stderr)
                self.assertIn("R_CFG_OK", r_out.stdout)


    def test_cis_mississippi_disjunct_filter_and_prune_eastern_cohort(self):
        """Verify Cis-Mississippi filter keeps P. dubia / P. obovata disjuncts up to -100W, drops P. cana & western P. paupercula, and prunes SQLite + CSVs."""
        import sqlite3
        from scripts.core.db import get_db_connection
        from scripts.data_prep.prune_eastern_cohort import prune_workspace_cohort

        with tempfile.TemporaryDirectory() as tmp_d:
            ws = Path(tmp_d)
            tables_dir = ws / "data" / "tables"
            tables_dir.mkdir(parents=True, exist_ok=True)
            db_p = ws / "data" / "packera_pipeline.sqlite"

            conn = get_db_connection(db_path=db_p, init_schema=True, workspace_root=ws)
            rows = [
                # 1. Packera dubia cis-Mississippi (NC) -> KEEP
                ("NCU_DUBIA_NC", "NCU", "Packera dubia", "Cineraria dubia", "North Carolina", 35.9, -79.0, "Tier_1_Gold", ""),
                # 2. Packera dubia trans-Mississippi disjunct in East Texas (-95.2W >= -100W) -> KEEP
                ("TEX_DUBIA_ETX", "TEX", "Packera dubia", "Senecio tomentosus", "Texas", 31.2, -95.2, "Tier_2_Silver", ""),
                # 3. Packera obovata Edwards Plateau TX disjunct (-98.6W >= -100W) -> KEEP
                ("BRIT_OBOVATA_TX", "BRIT", "Packera obovata", "Packera obovata", "Texas", 30.1, -98.6, "Tier_1_Gold", ""),
                # 4. Packera obovata West Texas stray (-103.9W < -100W) -> DROP
                ("SRSC_OBOVATA_WTX", "SRSC", "Packera obovata", "Packera obovata", "Texas", 30.6, -103.9, "Tier_3_Bronze", ""),
                # 5. Packera obovata Western Type Rescue (-104.5W Holotype) -> KEEP (nomenclatural type)
                ("NY_OBOVATA_TYPE", "NY", "Packera obovata", "Packera obovata", "New Mexico", 35.1, -104.5, "Tier_1_Gold", "Holotype"),
                # 6. Packera paupercula cis-Mississippi (NC) & Driftless IA (-91.5W) -> KEEP
                ("NCU_PAUP_NC", "NCU", "Packera paupercula", "Packera paupercula", "North Carolina", 36.1, -81.2, "Tier_1_Gold", ""),
                ("ISC_PAUP_IA", "ISC", "Packera paupercula", "Packera paupercula", "Iowa", 42.9, -91.5, "Tier_2_Silver", ""),
                # 7. Packera paupercula Western / Great Plains (SD -103.5W, TX -98.0W, BC -120.5W) -> DROP
                ("RM_PAUP_SD", "RM", "Packera paupercula", "Packera paupercula", "South Dakota", 43.8, -103.5, "Tier_2_Silver", ""),
                ("TEX_PAUP_TX", "TEX", "Packera paupercula", "Packera paupercula", "Texas", 30.4, -98.0, "Tier_3_Bronze", ""),
                ("UBC_PAUP_BC", "UBC", "Packera paupercula", "Packera paupercula", "British Columbia", 50.2, -120.5, "Tier_2_Silver", ""),
                # 8. Packera cana (Cordilleran/Plains with 1 stray IL sheet + 5 western sheets) -> DROP ENTIRE SPECIES
                ("ILL_CANA_STRAY", "ILL", "Packera cana", "Packera cana", "Illinois", 40.1, -89.2, "Tier_3_Bronze", ""),
                ("COLO_CANA_1", "COLO", "Packera cana", "Packera cana", "Colorado", 40.0, -105.3, "Tier_1_Gold", ""),
                ("COLO_CANA_2", "COLO", "Packera cana", "Packera cana", "Wyoming", 43.1, -108.2, "Tier_1_Gold", ""),
                ("MONT_CANA_3", "MONT", "Packera cana", "Packera cana", "Montana", 46.5, -112.0, "Tier_2_Silver", ""),
            ]
            with conn:
                for cat, inst, sn, sr, st, lat, lon, tier, tstat in rows:
                    conn.execute(
                        """
                        INSERT INTO vouchers (
                            catalogNumber, institutionCode, scientificName, species_raw, species_standardized,
                            stateProvince, decimalLatitude, decimalLongitude, determiner_tier, type_status
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                        """,
                        (cat, inst, sn, sr, sn, st, lat, lon, tier, tstat),
                    )
                    conn.execute(
                        "INSERT INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 0);",
                        (cat,),
                    )
                    conn.execute(
                        "INSERT INTO extracted_leaves (leaf_id, catalogNumber, plant_individual_id, leaf_zone, contour_xy_json) VALUES (?, ?, 0, 'basal', '[]');",
                        (f"{cat}_L1", cat),
                    )
            v_df = pd.read_sql_query("SELECT * FROM vouchers;", conn)
            v_df.to_csv(tables_dir / "curated_vouchers.csv", index=False)
            conn.close()

            res = prune_workspace_cohort(workspace_root=ws, db_path=db_p, vacuum=False, dry_run=False)
            expected_kept = {
                "NCU_DUBIA_NC",
                "TEX_DUBIA_ETX",
                "BRIT_OBOVATA_TX",
                "NY_OBOVATA_TYPE",
                "NCU_PAUP_NC",
                "ISC_PAUP_IA",
            }
            self.assertEqual(res["keep_catalog_count"], len(expected_kept))
            self.assertIn("Packera cana", res["dropped_species"])

            conn_chk = sqlite3.connect(str(db_p))
            try:
                rem_cats = {r[0] for r in conn_chk.execute("SELECT catalogNumber FROM vouchers;").fetchall()}
                rem_leaves = {r[0] for r in conn_chk.execute("SELECT catalogNumber FROM extracted_leaves;").fetchall()}
                self.assertEqual(rem_cats, expected_kept)
                self.assertEqual(rem_leaves, expected_kept)
            finally:
                conn_chk.close()

            pruned_csv = pd.read_csv(tables_dir / "curated_vouchers.csv")
            self.assertEqual(set(pruned_csv["catalogNumber"]), expected_kept)
            # Also verify missing extracted_leaf_manifest.csv was reconstructed automatically from SQLite!
            self.assertTrue((tables_dir / "extracted_leaf_manifest.csv").is_file())
            rebuilt_leaves = pd.read_csv(tables_dir / "extracted_leaf_manifest.csv")
            self.assertEqual(set(rebuilt_leaves["catalogNumber"]), expected_kept)


if __name__ == "__main__":
    unittest.main()
