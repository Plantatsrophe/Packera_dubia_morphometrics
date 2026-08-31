#!/usr/bin/env python3
"""
Unit tests for botanical voucher harvesting and geographic exclusion filters.
Tests the exclusion of western US states (farther west than Texas and Oklahoma)
and verifies standard botanical ingestion filtering.
"""

import unittest
from unittest.mock import patch, MagicMock
from pathlib import Path
import pandas as pd
import numpy as np

import tempfile
from PIL import Image

from scripts.core.config import (
    EXCLUDED_WESTERN_STATES,
    WESTERN_LONGITUDE_THRESHOLD,
    DEFAULT_TARGET_TAXA,
    DEFAULT_MIN_MEGAPIXELS,
)
from scripts.core.harvester_utils import (
    is_excluded_western_region,
    sanitize_filename,
    parse_determiner_tier,
    calculate_circular_phenology,
    infer_regional_group,
    optimize_herbarium_image_url,
    extract_high_res_image_url,
    validate_image_quality,
)
from scripts.core.harvester import harvest_taxa_occurrences


class TestVoucherHarvesterGeographicFiltering(unittest.TestCase):
    """Test suite for western region exclusion and geographic filtering."""

    def test_western_states_exclusion_full_names(self):
        """Verify that all full names of states farther west than TX/OK are excluded."""
        western_names = [
            "Colorado", "New Mexico", "Wyoming", "Montana", "Utah",
            "Idaho", "Arizona", "Nevada", "Washington", "Oregon",
            "California", "Alaska", "Hawaii"
        ]
        for state in western_names:
            with self.subTest(state=state):
                self.assertTrue(
                    is_excluded_western_region(state),
                    f"State '{state}' should be excluded."
                )

    def test_western_states_exclusion_abbreviations(self):
        """Verify that 2-letter postal codes for western states are excluded."""
        western_codes = [
            "CO", "NM", "WY", "MT", "UT", "ID", "AZ", "NV", "WA", "OR", "CA", "AK", "HI"
        ]
        for code in western_codes:
            with self.subTest(code=code):
                self.assertTrue(
                    is_excluded_western_region(code),
                    f"Code '{code}' should be excluded."
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
                    f"State '{state}' should be retained (not excluded)."
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
        # Specimen in Colorado/Utah region with missing state
        self.assertTrue(
            is_excluded_western_region(None, lat=39.5, lon=-108.5),
            "Unrecorded state with longitude < -106.65 should be excluded."
        )
        # Specimen in North Carolina with missing state
        self.assertFalse(
            is_excluded_western_region(None, lat=35.9, lon=-79.0),
            "Unrecorded state with Eastern longitude should be retained."
        )
        # Specimen in Texas with missing state
        self.assertFalse(
            is_excluded_western_region(None, lat=31.5, lon=-98.0),
            "Unrecorded state in central Texas longitude should be retained."
        )

    def test_harvester_ingestion_drops_western_records(self):
        """Verify harvest_taxa_occurrences filters out western records from mock API response."""
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
                "media": [{"identifier": "https://example.com/ncu1.jpg", "type": "StillImage", "format": "image/jpeg"}],
                "catalogNumber": "NCU001",
                "institutionCode": "NCU",
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
                "count": len(mock_results)
            }

            # Harvest with exclude_western=True (default)
            df = harvest_taxa_occurrences(
                taxa_list=["Packera paupercula"],
                max_records_per_taxon=10,
                exclude_western=True
            )

            # Only NCU001 (NC) and TEX001 (TX) should be retained
            self.assertEqual(len(df), 2)
            retained_catalogs = set(df["catalogNumber"])
            self.assertIn("NCU001", retained_catalogs)
            self.assertIn("TEX001", retained_catalogs)
            self.assertNotIn("COLO001", retained_catalogs)
            self.assertNotIn("RM001", retained_catalogs)


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
                "identifier": "https://media.symbiota.org/NCU/thumbnails/NCU001_tn.jpg"
            },
            {
                "type": "StillImage",
                "format": "image/jpeg",
                "identifier": "https://media.symbiota.org/NCU/original/NCU001_lg.jpg"
            },
        ]
        selected_url = extract_high_res_image_url(media_list)
        self.assertIn("_lg.jpg", selected_url)
        self.assertNotIn("_tn.jpg", selected_url)

    def test_validate_image_quality_high_vs_low_resolution(self):
        """Verify that validate_image_quality accepts high-res sheets and rejects low-res images."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)

            # 1. Create small low-res image (800 x 1000 = 0.8 MP)
            low_res_img = tmp_path / "low_res.jpg"
            img_small = Image.new("RGB", (800, 1000), color=(128, 128, 128))
            img_small.save(low_res_img, "JPEG")

            is_valid, metrics = validate_image_quality(low_res_img, min_megapixels=8.0, min_file_size_kb=1.0)
            self.assertFalse(is_valid)
            self.assertEqual(metrics["reason"], "low_resolution")
            self.assertLess(metrics["megapixels"], 8.0)

            # 2. Create high-res image (3000 x 4000 = 12.0 MP)
            high_res_img = tmp_path / "high_res.jpg"
            img_large = Image.new("RGB", (3000, 4000), color=(128, 128, 128))
            img_large.save(high_res_img, "JPEG")

            is_valid_high, metrics_high = validate_image_quality(high_res_img, min_megapixels=8.0, min_file_size_kb=1.0)
            self.assertTrue(is_valid_high)
            self.assertEqual(metrics_high["megapixels"], 12.0)


if __name__ == "__main__":
    unittest.main()
