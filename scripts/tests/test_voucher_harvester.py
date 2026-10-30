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
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
from PIL import Image

from scripts.core.config import (
    DEFAULT_MIN_MEGAPIXELS,
    DEFAULT_TARGET_TAXA,
    EXCLUDED_WESTERN_STATES,
    WESTERN_LONGITUDE_THRESHOLD,
)
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

            # Verify circular phenology features are DECOUPLED and not present
            self.assertNotIn("pheno_sin", df.columns)
            self.assertNotIn("pheno_cos", df.columns)
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
                    "species_raw": "Packera dubia",
                    "determiner_raw": "D.K. Trock",
                    "determiner_tier": "Tier_1_Gold",
                    "type_status": "None",
                    "county": "Wake",
                    "stateProvince": "North Carolina",
                    "latitude": 35.8,
                    "longitude": -78.6,
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

            # Validate backwards-compatible legacy column preservation
            legacy_headers = ["species_raw", "determiner_raw", "latitude", "longitude"]
            for col in legacy_headers:
                self.assertIn(col, read_df.columns, f"Missing backwards-compatible header '{col}'")

            # Check value consistency between DwC and legacy aliases
            self.assertEqual(read_df.iloc[0]["scientificName"], "Packera dubia")
            self.assertEqual(read_df.iloc[0]["species_raw"], "Packera dubia")
            self.assertEqual(read_df.iloc[0]["identifiedBy"], "D.K. Trock")
            self.assertEqual(read_df.iloc[0]["determiner_raw"], "D.K. Trock")
            self.assertAlmostEqual(read_df.iloc[0]["decimalLatitude"], 35.8)
            self.assertAlmostEqual(read_df.iloc[0]["latitude"], 35.8)
            self.assertAlmostEqual(read_df.iloc[0]["decimalLongitude"], -78.6)
            self.assertAlmostEqual(read_df.iloc[0]["longitude"], -78.6)

    def test_curated_vouchers_table_contract(self):
        """Verify production curated_vouchers.csv adheres to Darwin Core contract headers."""
        prod_csv = Path(__file__).resolve().parents[2] / "data" / "tables" / "curated_vouchers.csv"
        if not prod_csv.exists():
            self.skipTest("Production curated_vouchers.csv not present.")
        df = pd.read_csv(prod_csv)
        contract_cols = [
            "catalogNumber",
            "scientificName",
            "decimalLatitude",
            "decimalLongitude",
            "eventDate",
            "identifiedBy",
            "determiner_tier",
            "image_path",
        ]
        for col in contract_cols:
            self.assertIn(col, df.columns, f"Contract column '{col}' missing from curated_vouchers.csv")

        legacy_cols = ["species_raw", "determiner_raw", "latitude", "longitude"]
        for col in legacy_cols:
            self.assertIn(col, df.columns, f"Legacy column '{col}' missing from curated_vouchers.csv")


class TestVoucherHarvesterCLI(unittest.TestCase):
    """Test suite for 01_voucher_harvester.py CLI parser and main entrypoint."""

    def test_cli_parser_defaults(self):
        """Verify default CLI arguments match PipelineConfig specifications."""
        import importlib
        mod_cli = importlib.import_module("scripts.data_prep.01_voucher_harvester")
        parser = mod_cli.build_cli_parser()
        args = parser.parse_args([])

        self.assertIsInstance(args.taxa, list)
        self.assertGreater(len(args.taxa), 0)
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
            ("No. 1042b", "1042"),
            ("1042-A", "1042"),
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
                "recordNumber": "1042b",
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
                "latitude": 35.91230,
                "longitude": -79.05120,
                "determiner_tier": "Tier_2_Silver",
            },
            {
                "catalogNumber": "SHEET_B",
                "institutionCode": "GA",
                "recordedBy": "Fuller",
                "recordNumber": "",
                "eventDate": "1984-05-12",
                "latitude": 35.91245,  # ~17 m away (< 50 m)
                "longitude": -79.05125,
                "determiner_tier": "Tier_3_Bronze",
            },
            {
                "catalogNumber": "SHEET_C",
                "institutionCode": "MO",
                "recordedBy": "Fuller",
                "recordNumber": "s.n.",
                "eventDate": "1984-05-12",
                "latitude": 36.50000,  # ~65 km away (> 50 m)
                "longitude": -79.05120,
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
            {"catalogNumber": "NC_VOUCHER_010", "leaf_id": 1, "assigned_tier": "Tier_1_Direct"},
            {"catalogNumber": "NC_VOUCHER_010", "leaf_id": 2, "assigned_tier": "Tier_1_Direct"},
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


if __name__ == "__main__":
    unittest.main()

