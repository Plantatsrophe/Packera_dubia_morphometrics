#!/usr/bin/env python3
"""
Step 01 statistical, taxonomic, and biogeographic robustness audit tests.
"""

from __future__ import annotations

import importlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from PIL import Image, ImageFilter

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core import taxonomy_authority as ta
from scripts.core.config import PipelineConfig
from scripts.core.db import get_db_connection
from scripts.core.harvester import VoucherHarvester, validate_image_quality
from scripts.core.step01_biogeography_utils import (
    evaluate_transcontinental_ecoregion_gate,
    normalize_packera_orthography,
    resolve_eastern_resident_whitelist,
    resolve_synonym_genera,
)
from scripts.tests._r_session import rscript_available, run_warm_r_cmd


def _mk_occ(key: int, cat: str, inst: str, state: str, lat: float, lon: float, det: str = "T.M. Barkley") -> dict:
    return {
        "key": key,
        "scientificName": "Packera dubia",
        "stateProvince": state,
        "decimalLatitude": lat,
        "decimalLongitude": lon,
        "coordinateUncertaintyInMeters": 50.0,
        "year": 2020,
        "month": 5,
        "day": 10,
        "recordedBy": f"Collector {key}",
        "recordNumber": str(key),
        "identifiedBy": det,
        "media": [{"identifier": f"https://example.org/{cat}.jpg", "type": "StillImage", "format": "image/jpeg"}],
        "catalogNumber": cat,
        "institutionCode": inst,
    }


class TestStep01RobustnessAudit(unittest.TestCase):
    """Focused verification for Step 01 pagination, authority, blur, ecoregion, and weight invariants."""

    def test_still_image_and_two_pass_country_quota_allocation(self) -> None:
        """Verify mediaType='StillImage', fair US/CA Pass-1 allocation, and CA->US capacity rollover."""
        us_large = [_mk_occ(1000 + i, f"US_{i:02d}", "NCU", "North Carolina", 35.9, -79.0 - i * 0.01) for i in range(12)]
        ca_large = [_mk_occ(2000 + i, f"CA_{i:02d}", "TRT", "Ontario", 44.0, -79.5 - i * 0.01) for i in range(12)]
        seen_calls = []

        def _occ_balanced(**kw):
            seen_calls.append(dict(kw))
            if "typeStatus" in kw:
                return {"results": [], "count": 0}
            pool = us_large if kw.get("country") == "US" else ca_large
            off = int(kw.get("offset", 0))
            lim = int(kw.get("limit", 300))
            return {"results": pool[off:off + lim], "count": len(pool)}

        with patch("scripts.core.harvester.occ.search", side_effect=_occ_balanced):
            harv = VoucherHarvester(taxa=["Packera dubia"], max_records_per_taxon=10, expand_synonyms=False)
            df_bal = harv.harvest()

        main_calls = [c for c in seen_calls if "typeStatus" not in c]
        self.assertTrue(main_calls and all(c.get("mediaType") == "StillImage" for c in main_calls))
        self.assertEqual(len(df_bal), 10)
        us_n = int(df_bal["catalogNumber"].str.startswith("US_").sum())
        ca_n = int(df_bal["catalogNumber"].str.startswith("CA_").sum())
        self.assertEqual((us_n, ca_n), (5, 5))

        # Roll unused CA capacity back to US when CA holds only 2 records
        ca_small = ca_large[:2]
        seen_calls.clear()

        def _occ_rollover(**kw):
            seen_calls.append(dict(kw))
            if "typeStatus" in kw:
                return {"results": [], "count": 0}
            pool = us_large if kw.get("country") == "US" else ca_small
            off = int(kw.get("offset", 0))
            lim = int(kw.get("limit", 300))
            return {"results": pool[off:off + lim], "count": len(pool)}

        with patch("scripts.core.harvester.occ.search", side_effect=_occ_rollover):
            harv2 = VoucherHarvester(taxa=["Packera dubia"], max_records_per_taxon=10, expand_synonyms=False)
            df_roll = harv2.harvest()

        self.assertEqual(len(df_roll), 10)
        self.assertEqual(int(df_roll["catalogNumber"].str.startswith("US_").sum()), 8)
        self.assertEqual(int(df_roll["catalogNumber"].str.startswith("CA_").sum()), 2)

    def test_authority_decircularization_no_self_gold_no_protologue_and_silver_cap(self) -> None:
        """Verify investigator self-ID, protologue removal, and corpus-specialist Tier_2_Silver cap."""
        self.assertFalse(hasattr(ta, "extract_" + "protologue_authors"))

        for det_name in ("J. Brandon Fuller", "B. Fuller", "K. Sprengel", "F. Michaux"):
            tier_ncu, _, _ = ta.parse_determiner_tier(
                None, det_name, "Field Collector", None, "NCU", "Granite flatrock 3 km E of Wake Forest", None
            )
            tier_unk, _, _ = ta.parse_determiner_tier(
                None, det_name, "Field Collector", None, "UNKNOWN_HERB", "Roadside", None
            )
            self.assertEqual(tier_ncu, "Tier_2_Silver")
            self.assertEqual(tier_unk, "Tier_3_Bronze")

        rows = [
            {
                "catalogNumber": f"CORP_{i}",
                "institutionCode": f"INST_{i}",
                "scientificName": "Packera anonyma",
                "species_raw": "Packera anonyma",
                "identifiedBy": "Q. Specialist",
                "recordedBy": f"Collector {i}",
                "type_status": "None",
                "determiner_tier": "Tier_3_Bronze",
                "locality": "Granite outcrop slope",
                "habitat": "Granite",
            }
            for i in range(6)
        ]
        df_c = pd.DataFrame(rows)
        specs, herbs = ta.discover_corpus_authorities(df_c)
        self.assertIn("specialist", specs)
        scored = ta.apply_corpus_authority_tiers(df_c, dynamic_specialists=specs, dynamic_herbaria=herbs)
        self.assertTrue((scored["determiner_tier"] == "Tier_2_Silver").all())

    def test_stratify_duplicates_gold_first_and_zero_silent_relabeling(self) -> None:
        """Verify Tier_1_Gold beats home_herbarium and conflicting sheets remain unmutated."""
        df_dups = pd.DataFrame([
            {
                "catalogNumber": "NCU_BRZ",
                "institutionCode": "NCU",
                "scientificName": "Packera aurea",
                "species_raw": "Packera aurea",
                "species_standardized": "Packera aurea",
                "recordedBy": "Collector A",
                "recordNumber": "901",
                "eventDate": "2020-05-12",
                "identifiedBy": "Generalist Bot",
                "determiner_tier": "Tier_3_Bronze",
                "type_status": "None",
            },
            {
                "catalogNumber": "MO_GLD",
                "institutionCode": "MO",
                "scientificName": "Packera obovata",
                "species_raw": "Packera obovata",
                "species_standardized": "Packera obovata",
                "recordedBy": "Collector A",
                "recordNumber": "901",
                "eventDate": "2020-05-12",
                "identifiedBy": "T.M. Barkley",
                "determiner_tier": "Tier_1_Gold",
                "type_status": "None",
            },
        ])
        strat = ta.stratify_duplicates(df_dups, home_herbarium="NCU").set_index("catalogNumber")
        self.assertTrue(bool(strat.loc["MO_GLD", "is_primary_duplicate"]))
        self.assertFalse(bool(strat.loc["NCU_BRZ", "is_primary_duplicate"]))
        self.assertTrue(bool(strat.loc["NCU_BRZ", "exsiccatae_taxonomic_conflict"]))
        self.assertEqual(strat.loc["NCU_BRZ", "exsiccatae_gold_taxon"], "Packera obovata")
        self.assertEqual(strat.loc["NCU_BRZ", "scientificName"], "Packera aurea")
        self.assertEqual(strat.loc["NCU_BRZ", "determiner_tier"], "Tier_3_Bronze")
        self.assertEqual(strat.loc["NCU_BRZ", "identifiedBy"], "Generalist Bot")

    def test_indumentum_safe_tile_p90_sharpness_gate(self) -> None:
        """Verify tile-p90 Laplacian accepts crisp label/margin edges on smooth sheets and rejects blurred scans."""
        cfg = PipelineConfig.from_yaml()
        self.assertTrue(cfg.harvesting.check_sharpness)
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            arr = np.full((3000, 4000, 3), 145, dtype=np.uint8)
            # Crisp checkerboard label/margin patch occupying ~16% of central ROI
            block = np.zeros((750, 950, 3), dtype=np.uint8)
            block[::4, ::4] = 255
            arr[1400:2150, 1800:2750] = block
            tom_p = tmp_p / "tomentose_crisp.png"
            Image.fromarray(arr, mode="RGB").save(tom_p, format="PNG")
            ok_tom, m_tom = validate_image_quality(
                tom_p, min_megapixels=8.0, min_file_size_kb=1.0, check_sharpness=True, min_sharpness=80.0
            )
            self.assertTrue(ok_tom, f"Expected tomentose + crisp edge scan to pass: {m_tom}")
            self.assertGreaterEqual(float(m_tom["sharpness"]), 80.0)

            blur_p = tmp_p / "blurred.png"
            Image.fromarray(arr, mode="RGB").filter(ImageFilter.GaussianBlur(radius=6.0)).save(blur_p, format="PNG")
            ok_blur, m_blur = validate_image_quality(
                blur_p, min_megapixels=8.0, min_file_size_kb=1.0, check_sharpness=True, min_sharpness=80.0
            )
            self.assertFalse(ok_blur)
            self.assertEqual(m_blur["reason"], "blurry_or_upscaled")

    def test_ecoregion_bounded_transcontinental_gate_and_sqlite_whitelist(self) -> None:
        """Verify northern MN / Boreal Shield P. paupercula passes while Plains/Cordilleran records drop."""
        self.assertEqual(resolve_synonym_genera(), {"Packera", "Senecio", "Cineraria"})
        self.assertEqual(normalize_packera_orthography("Packera millefolia"), "Packera millefolium")
        self.assertEqual(normalize_packera_orthography("Packera ×memmingeri"), "Packera memmingeri")

        drop_mn, reason_mn = evaluate_transcontinental_ecoregion_gate(
            "Minnesota", lat=48.5, lon=-96.2, regional_group="Boreal_Shield"
        )
        self.assertFalse(drop_mn)
        self.assertEqual(reason_mn, "boreal_ecoregion_extension")

        drop_sd, reason_sd = evaluate_transcontinental_ecoregion_gate(
            "South Dakota", lat=43.8, lon=-103.5, regional_group="Northwestern_Great_Plains"
        )
        self.assertTrue(drop_sd)
        self.assertEqual(reason_sd, "excluded_plains_cordilleran")

        df_coh = pd.DataFrame([
            {"catalogNumber": "NC_P", "scientificName": "Packera paupercula", "stateProvince": "North Carolina", "decimalLatitude": 36.1, "decimalLongitude": -81.2},
            {"catalogNumber": "MN_B", "scientificName": "Packera paupercula", "stateProvince": "Minnesota", "decimalLatitude": 48.5, "decimalLongitude": -96.2, "regional_group": "Northern_Lakes_and_Forests"},
            {"catalogNumber": "SD_P", "scientificName": "Packera paupercula", "stateProvince": "South Dakota", "decimalLatitude": 43.8, "decimalLongitude": -103.5, "regional_group": "Northwestern_Great_Plains"},
            {"catalogNumber": "TX_P", "scientificName": "Packera paupercula", "stateProvince": "Texas", "decimalLatitude": 30.4, "decimalLongitude": -98.0, "regional_group": "Edwards_Plateau"},
            {"catalogNumber": "BC_P", "scientificName": "Packera paupercula", "stateProvince": "British Columbia", "decimalLatitude": 50.2, "decimalLongitude": -120.5},
        ])
        keep_m, stats = ta.evaluate_eastern_na_cohort_filter(df_coh, lambda x: str(x).strip())
        self.assertEqual(set(df_coh.loc[keep_m, "catalogNumber"]), {"NC_P", "MN_B"})
        self.assertGreater(stats["biogeographic_boundary_status"].get("boreal_ecoregion_extension", 0), 0)

        with tempfile.TemporaryDirectory() as tmp_d:
            db_p = Path(tmp_d) / "test.sqlite"
            conn = get_db_connection(db_path=db_p, init_schema=True)
            with conn:
                conn.execute(
                    "INSERT INTO taxon_backbone_profiles (species_standardized, tdwg_native_codes_json) VALUES (?, ?);",
                    ("Packera eurycephala", '["NCA"]'),
                )
                conn.execute(
                    "INSERT INTO natureserve_subnational_ranges (species_standardized, country_code, subnation_code, is_native) VALUES (?, 'US', 'VA', 1);",
                    ("Packera rosei",),
                )
            conn.close()
            wl = resolve_eastern_resident_whitelist(db_path=db_p)
            self.assertIn("Packera eurycephala", wl)
            self.assertIn("Packera rosei", wl)
            self.assertIn("Packera millefolium", wl)

    def test_lossless_step01_and_active_anchor_weights_python_and_r(self) -> None:
        """Verify neutral Step 01 weights, CLI flag removal, and centroid/Bronze exclusion in Python & R."""
        cli_mod = importlib.import_module("scripts.data_prep.01_voucher_harvester")
        parser = cli_mod.build_cli_parser()
        for removed_flag in ("--spatial-thin-km", "--max-per-state", "--max-per-region"):
            with self.assertRaises(SystemExit):
                parser.parse_args([removed_flag, "10"])

        recs = [
            _mk_occ(3001, "NCU_G1", "NCU", "North Carolina", 35.90, -79.05),
            _mk_occ(3002, "NCU_G2", "NCU", "North Carolina", 35.91, -79.06),
            _mk_occ(3003, "AR_G3", "UARK", "Arkansas", 35.50, -93.80),
        ]
        with patch("scripts.core.harvester.occ.search", side_effect=lambda **kw: {"results": [] if ("typeStatus" in kw or kw.get("country") == "CA") else recs, "count": 3}):
            harv = VoucherHarvester(taxa=["Packera dubia"], max_records_per_taxon=10, expand_synonyms=False)
            df_h = harv.harvest()
        self.assertEqual(len(df_h), 3)
        self.assertTrue((df_h["spatial_balance_weight"] == 1.0).all())

        df_tier_b = pd.DataFrame([
            {"catalogNumber": "A1", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.90, "decimalLongitude": -79.05, "is_centroid_or_inst_coord": False, "authority_weight": 2.0},
            {"catalogNumber": "A2", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.91, "decimalLongitude": -79.06, "is_centroid_or_inst_coord": False, "authority_weight": 2.0},
            {"catalogNumber": "A3", "species_standardized": "Packera dubia", "determiner_tier": "Tier_2_Silver", "decimalLatitude": 35.50, "decimalLongitude": -93.80, "is_centroid_or_inst_coord": False, "authority_weight": 0.75},
            {"catalogNumber": "CENT", "species_standardized": "Packera dubia", "determiner_tier": "Tier_1_Gold", "decimalLatitude": 35.50, "decimalLongitude": -93.80, "is_centroid_or_inst_coord": True, "authority_weight": 2.0},
            {"catalogNumber": "BRZ", "species_standardized": "Packera dubia", "determiner_tier": "Tier_3_Bronze", "decimalLatitude": 35.50, "decimalLongitude": -93.80, "is_centroid_or_inst_coord": False, "authority_weight": 0.0},
        ])
        py_w = ta.spatially_balance_and_thin_vouchers(df_tier_b, primary_duplicates_only=False).set_index("catalogNumber")
        self.assertEqual(float(py_w.loc["CENT", "spatial_balance_weight"]), 1.0)
        self.assertEqual(float(py_w.loc["BRZ", "spatial_balance_weight"]), 1.0)
        self.assertGreater(float(py_w.loc["A3", "spatial_balance_weight"]), float(py_w.loc["A1", "spatial_balance_weight"]))

        if rscript_available():
            r_script = (
                "source('scripts/morphometrics/pub_figures_utils.R'); "
                "df <- data.frame("
                "species_standardized=rep('Packera dubia',5), "
                "decimalLatitude=c(35.90,35.91,35.50,35.50,35.50), "
                "decimalLongitude=c(-79.05,-79.06,-93.80,-93.80,-93.80), "
                "is_centroid_or_inst_coord=c(FALSE,FALSE,FALSE,TRUE,FALSE), "
                "authority_weight=c(2.0,2.0,0.75,2.0,0.0), "
                "spatial_balance_weight=rep(1.0,5), stringsAsFactors=FALSE); "
                "w <- compute_spatial_balance_weights(df); "
                "stopifnot(abs(w[4]-1.0)<1e-6, abs(w[5]-1.0)<1e-6, w[3]>w[1]); "
                "cat('R_SPATIAL_OK\\n')"
            )
            r_out = run_warm_r_cmd(["Rscript", "-e", r_script], stage_key="step03", cwd=PROJECT_ROOT)
            self.assertEqual(r_out.returncode, 0, r_out.stderr)
            self.assertIn("R_SPATIAL_OK", r_out.stdout)


if __name__ == "__main__":
    unittest.main()

