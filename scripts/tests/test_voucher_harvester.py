#!/usr/bin/env python3
"""
===============================================================================
Unit tests for botanical voucher harvesting, authority stratification,
Darwin Core metadata normalization, and atomic persistence.
===============================================================================
"""

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
    harvest_taxa_occurrences,
    infer_regional_group,
    is_excluded_western_region,
    optimize_herbarium_image_url,
    parse_determiner_tier,
    sanitize_filename,
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


if __name__ == "__main__":
    unittest.main()
