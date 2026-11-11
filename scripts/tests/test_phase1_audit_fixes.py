#!/usr/bin/env python3
"""
===============================================================================
Script: test_phase1_audit_fixes.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Regression test suite verifying all Phase 1 Scientific, Statistical, and
    Software Audit remediation patches across Steps 01-07.
===============================================================================
"""

from __future__ import annotations

import importlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import main as root_main
from scripts.core.config import PipelineConfig
from scripts.core.epa_ecoregions import infer_regional_group
from scripts.core.harvester import VoucherHarvester
from scripts.vision import lm3_data_loader
from scripts.vision.configure_leafmachine3 import NATIVE_LM3_SETTINGS_PATH, load_native_lm3_config
from scripts.vision.lm3_compat import DBSCAN


class TestPhase1AuditFixes(unittest.TestCase):
    """Verify all Phase 1 codebase bug fixes and statistical anti-leakage safeguards."""

    def test_config_and_cli_readonly_lm3_and_usda_defaults(self) -> None:
        """Verify Phase 1 defaults: legacy Python gatekeeper purged, usda_county_filter=False, use_cf_predicted_by_mp=False, min_flowering_doy=30, permutations=999."""
        cfg = PipelineConfig.from_yaml()
        self.assertFalse(cfg.harvesting.usda_county_filter)
        self.assertFalse(cfg.segmentation.use_cf_predicted_by_mp)
        self.assertEqual(cfg.macroecology.phenology.min_flowering_doy, 30)
        self.assertFalse(any(p.name.endswith("gatekeeper.py") for p in (PROJECT_ROOT / "scripts/vision").glob("*.py")))

        parser = root_main.build_parser()
        synth_args = parser.parse_args(["synthesis"])
        self.assertEqual(synth_args.permutations, 999)

        step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        with patch.object(sys, "argv", ["02_segment_and_extract.py"]):
            p02_def = step02.parse_args()
            self.assertEqual(p02_def.score_thresh, cfg.segmentation.score_thresh)

    def test_step02_readonly_lm3_sqlite_persistence_and_step03_qc_contract(self) -> None:
        """Verify Step 02 persists 100% of LM3 leaves + native QC flags to SQLite and Step 03 filters clumps in R."""
        import cv2
        from scripts.core.db import get_db_connection, load_step02_records_from_db, record_voucher_segmentation_transaction

        step02 = importlib.import_module("scripts.pipeline.02_segment_and_extract")
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            dummy_w = tmp_p / "dummy.onnx"
            dummy_w.write_bytes(b"0")
            v_csv = tmp_p / "vouchers.csv"
            img_p = tmp_p / "NCU_RO_01.jpg"
            img_p.write_bytes(b"\xff\xd8\xff" + b"0" * 200)
            pd.DataFrame([{"catalogNumber": "NCU_RO_01", "image_path": str(img_p), "scientificName": "Packera aurea"}]).to_csv(
                v_csv, index=False
            )
            db_file = tmp_p / "readonly.sqlite"
            with patch.object(step02, "LM3LeafSegmenterEngine"), patch.object(step02, "LM3PlantDetector"):
                pipe = step02.SegmentAndExtractPipeline(
                    vouchers_csv=v_csv,
                    model_weights=dummy_w,
                    stage1_weights=dummy_w,
                    output_dir=tmp_p / "out",
                    device="cpu",
                    db_path=db_file,
                )
            pipe.ruler_calibrations["NCU_RO_01"] = 12.0

            sheet = np.full((600, 800, 3), 232, dtype=np.uint8)
            clean_leaf = np.zeros((600, 800), dtype=np.uint8)
            cv2.ellipse(clean_leaf, (200, 300), (45, 95), 0, 0, 360, 255, -1)

            clump_leaf = np.zeros((600, 800), dtype=np.uint8)
            pts_clump = np.array(
                [[520, 420], [516, 360], [514, 295], [520, 230], [548, 252], [608, 268], [616, 330], [582, 385], [532, 418]],
                dtype=np.int32,
            )
            cv2.fillPoly(clump_leaf, [pts_clump], 255)

            with patch("cv2.imread", return_value=sheet):
                with patch.object(
                    pipe.model_engine,
                    "predict",
                    return_value=[
                        ((180, 140, 420, 260), clean_leaf, 0.96, 0),
                        ((210, 500, 440, 630), clump_leaf, 0.91, 0),
                    ],
                ):
                    recs, fail = pipe.process_voucher("NCU_RO_01", img_p)

            self.assertIsNone(fail)
            self.assertEqual(len(recs), 2, "Read-only Step 02 must persist both clean and clump LM3 leaves.")
            gates_vals = {bool(r["lm3_gates_pass"]) for r in recs}
            self.assertEqual(gates_vals, {True, False})

            conn = get_db_connection(db_path=db_file, init_schema=True)
            record_voucher_segmentation_transaction(
                conn,
                catalog_number="NCU_RO_01",
                image_path=str(img_p),
                extracted_records=recs,
                failure_record=None,
                output_dir=pipe.output_dir,
            )
            rows = conn.execute(
                "SELECT leaf_id, bilateral_dice, si_a, lm3_gates_pass, n_lamina_components, scale_source FROM extracted_leaves WHERE catalogNumber = 'NCU_RO_01' ORDER BY leaf_id"
            ).fetchall()
            self.assertEqual(len(rows), 2)
            self.assertEqual({int(r["lm3_gates_pass"]) for r in rows}, {0, 1})
            loaded_ext, _, _, _ = load_step02_records_from_db(conn, ["NCU_RO_01"])
            conn.close()
            self.assertEqual(len(loaded_ext["NCU_RO_01"]), 2)

        s_db_r = (PROJECT_ROOT / "scripts/core/db_utils.R").read_text(encoding="utf-8")
        s_03_r = (PROJECT_ROOT / "scripts/morphometrics/03_fourier_extractor.R").read_text(encoding="utf-8")
        self.assertIn("evaluate_step03_lm3_efa_qc", s_db_r)
        self.assertIn("evaluate_step03_lm3_efa_qc", s_03_r)

    def test_audit_only_usda_county_tagging_retains_disjuncts_by_default(self) -> None:
        """Verify usda_county_filter=False retains out-of-range vouchers while persisting usda_county_status in SQLite."""
        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            usda_dir = tmp_p / "usda_plants"
            (usda_dir / "distributions").mkdir(parents=True, exist_ok=True)
            db_file = tmp_p / "audit_usda.sqlite"
            out_csv = tmp_p / "curated_vouchers.csv"

            cnty_geojson_path = tmp_p / "counties.geojson"
            cnty_geojson_path.write_text(
                json.dumps({
                    "type": "FeatureCollection",
                    "features": [
                        {
                            "type": "Feature",
                            "properties": {"state_name": "north carolina", "county_norm": "orange", "centroid_lon": -79.10, "centroid_lat": 36.05},
                            "geometry": {"type": "Polygon", "coordinates": [[[-79.25, 35.90], [-78.95, 35.90], [-78.95, 36.20], [-79.25, 36.20], [-79.25, 35.90]]]},
                        },
                        {
                            "type": "Feature",
                            "properties": {"state_name": "new york", "county_norm": "essex", "centroid_lon": -73.90, "centroid_lat": 44.10},
                            "geometry": {"type": "Polygon", "coordinates": [[[-74.20, 43.80], [-73.50, 43.80], [-73.50, 44.40], [-74.20, 44.40], [-74.20, 43.80]]]},
                        },
                    ],
                }),
                encoding="utf-8",
            )
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

            mock_occ = MagicMock()
            mock_occ.search.side_effect = lambda **kw: (
                {"results": [], "endOfRecords": True}
                if "typeStatus" in kw
                else {
                    "results": [
                        {
                            "catalogNumber": "NCU-IN-01", "institutionCode": "NCU", "scientificName": "Packera dubia",
                            "recordedBy": "Coll A", "recordNumber": "1", "identifiedBy": "Botanist",
                            "stateProvince": "North Carolina", "county": "Orange",
                            "decimalLatitude": 36.05, "decimalLongitude": -79.10, "coordinateUncertaintyInMeters": 50.0,
                            "eventDate": "2022-04-15", "year": 2022, "month": 4, "day": 15,
                            "media": [{"type": "StillImage", "identifier": "https://example.org/1.jpg"}],
                        },
                        {
                            "catalogNumber": "NY-DISJUNCT-02", "institutionCode": "NY", "scientificName": "Packera dubia",
                            "recordedBy": "Coll B", "recordNumber": "2", "identifiedBy": "Student",
                            "stateProvince": "New York", "county": "Essex",
                            "decimalLatitude": 44.10, "decimalLongitude": -73.90, "coordinateUncertaintyInMeters": 150.0,
                            "eventDate": "2022-05-15", "year": 2022, "month": 5, "day": 15,
                            "media": [{"type": "StillImage", "identifier": "https://example.org/2.jpg"}],
                        },
                    ],
                    "endOfRecords": True,
                }
            )

            mock_sp = MagicMock()
            mock_sp.name_backbone.return_value = {
                "matchType": "EXACT",
                "status": "ACCEPTED",
                "usageKey": 5387757,
                "canonicalName": "Packera dubia",
                "species": "Packera dubia",
            }
            mock_sp.name_usage.return_value = {
                "results": [
                    {"canonicalName": "Packera tomentosa", "scientificName": "Packera tomentosa (Michx.) C. Jeffrey", "taxonomicStatus": "SYNONYM"},
                    {"canonicalName": "Senecio tomentosus", "scientificName": "Senecio tomentosus Michx.", "taxonomicStatus": "SYNONYM"},
                ]
            }
            with patch("scripts.core.harvester.occ", mock_occ), patch("scripts.core.harvester.species", mock_sp):
                harv = VoucherHarvester(
                    taxa=["Packera dubia"],
                    expand_synonyms=False,
                    output_csv=out_csv,
                    raw_dir=tmp_p / "raw",
                    workspace_dir=tmp_p,
                    db_path=db_file,
                    usda_cache_dir=usda_dir,
                    usda_county_geojson_path=cnty_geojson_path,
                )
                self.assertFalse(harv.usda_county_filter)
                df_res = harv.harvest()
                harv.export(df_res)

            self.assertEqual(set(df_res["catalogNumber"]), {"NCU-IN-01", "NY-DISJUNCT-02"})
            import sqlite3
            conn = sqlite3.connect(str(db_file))
            try:
                rows = dict(conn.execute("SELECT catalogNumber, usda_county_status FROM vouchers;").fetchall())
            finally:
                conn.close()
            self.assertEqual(rows["NCU-IN-01"], "in_range")
            self.assertEqual(rows["NY-DISJUNCT-02"], "out_of_range")

    def test_scale_source_separation_and_baxter_county_ozarks(self) -> None:
        """Verify scale_source distinguishes ruler_calibrated vs mp_estimated, and Baxter Co. AR maps to Ozark_Highlands."""
        self.assertEqual(
            infer_regional_group(None, None, state_province="Arkansas", county="Baxter County"),
            "Ozark_Highlands",
        )
        self.assertTrue(NATIVE_LM3_SETTINGS_PATH.exists())
        lm3_cfg = load_native_lm3_config()
        self.assertEqual(int(lm3_cfg["modules"]["ect"]["num_dirs"]), 360)

        with tempfile.TemporaryDirectory() as tmp_d:
            rdir = Path(tmp_d) / "Data" / "Ruler"
            rdir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([
                {
                    "filename": "NCU_RUL__R__1", "ruler_success": True,
                    "conversion_mean": 120.0, "ruler_class_confidence": 95.0,
                    "conversion_determination": "measured_from_ruler",
                },
                {
                    "filename": "NCU_MP__R__1", "ruler_success": True,
                    "conversion_mean": 100.0, "ruler_class_confidence": 85.0,
                    "conversion_determination": "predicted_from_megapixels",
                },
            ]).to_csv(rdir / "ruler_summary.csv", index=False)
            recs = lm3_data_loader.load_ruler_calibration_records(Path(tmp_d))
            self.assertEqual(recs["NCU_RUL"]["scale_source"], "ruler_calibrated")
            self.assertTrue(recs["NCU_RUL"]["is_empirical_ruler"])
            self.assertEqual(recs["NCU_MP"]["scale_source"], "mp_estimated")
            self.assertFalse(recs["NCU_MP"]["is_empirical_ruler"])

    def test_numpy_dbscan_fallback_clusters_and_outliers(self) -> None:
        """Verify lm3_compat DBSCAN separates spatial clusters and marks isolated noise as -1."""
        pts = np.array([
            [10.0, 10.0], [11.0, 10.5], [9.5, 10.2],
            [100.0, 100.0], [101.0, 99.5], [99.8, 100.4],
            [500.0, 500.0],
        ])
        labels = DBSCAN(eps=5.0, min_samples=2).fit_predict(pts)
        self.assertEqual(len(set(labels[:3])), 1)
        self.assertEqual(len(set(labels[3:6])), 1)
        self.assertNotEqual(labels[0], labels[3])
        self.assertEqual(int(labels[6]), -1)

    def test_r_phase1_statistical_validity_contracts(self) -> None:
        """Verify Steps 03-07 R scripts enforce Phase 1 anti-leakage, log10(CS), uncertainty, IPP, and circular p-value contracts."""
        r03 = (PROJECT_ROOT / "scripts/morphometrics/03_fourier_extractor.R").read_text(encoding="utf-8")
        r04 = (PROJECT_ROOT / "scripts/morphometrics/04_gmm_morphotools.R").read_text(encoding="utf-8")
        raw = (PROJECT_ROOT / "scripts/morphometrics/authority_weighting_utils.R").read_text(encoding="utf-8")
        rman = (PROJECT_ROOT / "scripts/morphometrics/pub_ecoregion_manova_utils.R").read_text(encoding="utf-8")
        red = (PROJECT_ROOT / "scripts/analysis/edaphic_fine_grain_utils.R").read_text(encoding="utf-8")
        r05 = (PROJECT_ROOT / "scripts/analysis/05_multimodal_spatial_rf.R").read_text(encoding="utf-8")
        rph = (PROJECT_ROOT / "scripts/analysis/step05_phenology_ecotypes.R").read_text(encoding="utf-8")
        rsdm = (PROJECT_ROOT / "scripts/analysis/step05_sdm_maxent.R").read_text(encoding="utf-8")
        r07 = (PROJECT_ROOT / "scripts/analysis/07_discover_taxon_complexes.R").read_text(encoding="utf-8")

        # Step 03: base-10 log10(centsize), centroid_size_px separated, NA when uncalibrated
        self.assertIn("log10(centsize)", r03)
        self.assertIn("centroid_size_px", r03)
        self.assertIn("scale_source", r03)

        # Step 04: zero pre-CV deme shrinkage, zero soft geo prior in CDA, zero u_cut <= 0.70 boundary trimming, include_gmm_clusters = !is_intra_mode
        self.assertNotIn("apply_hierarchical_population_deme_shrinkage(closed_df", r04)
        self.assertIn("assign_spatial_gathering_demes(closed_df", r04)
        self.assertNotIn("quantile(unc_v, 0.70", r04)
        self.assertIn("include_gmm_clusters = !is_intra_mode", r04)
        self.assertIn("include_gmm_clusters = TRUE", rman)
        self.assertIn("build_taxon_blind_gathering_deme_folds", raw)
        self.assertNotIn("geo_ll_mat[i, ]", raw)

        # Steps 05-07: <= 2000 m edaphic gate, ILR texture, include_pheno = FALSE, Hopkins lon+year, cv.glmnet IPP, Watson-Williams circular p-values
        self.assertIn("max_edaphic_uncertainty_m = 2000.0", red)
        self.assertIn("ilr_texture_1", red)
        self.assertIn("include_pheno = FALSE", r05)
        self.assertIn("b_lon", rph)
        self.assertIn("b_year", rph)
        self.assertNotIn("rep(pres_idx, times = rep_p)", rsdm)
        self.assertIn("glmnet::cv.glmnet", rsdm)
        self.assertNotIn("exp(-(delta_circ / 18.0))", r07)
        self.assertIn(".circular_watson_williams_or_perm_p", r07)

    def test_phase1_nine_audit_regressions(self) -> None:
        """Verify fixes for aiohttp counter tallying, Step 03 CLI/OTU/QC contracts, ba_cm2, and Step 06 multi-plant leaf_count coalescing."""
        import asyncio
        from scripts.core.harvester import download_all_voucher_images

        with tempfile.TemporaryDirectory() as tmp_d:
            tp = Path(tmp_d)
            img_exist = tp / "NCU_EXIST.jpg"
            img_exist.write_bytes(b"\xff\xd8\xff" + b"Z" * 1500)
            img_pend = tp / "NCU_PEND.jpg"
            records = [
                ("https://example.com/1.jpg", img_exist),
                ("https://example.com/2.jpg", img_pend),
            ]

            class _DummySessionCM:
                async def __aenter__(self):
                    return MagicMock()

                async def __aexit__(self, exc_type, exc, tb):
                    return False

            fake_aiohttp = MagicMock()
            fake_aiohttp.TCPConnector.return_value = MagicMock()
            fake_aiohttp.ClientSession.return_value = _DummySessionCM()

            with patch("scripts.core.harvester.aiohttp", fake_aiohttp):
                with patch("scripts.core.harvester.download_single_image", return_value=(True, False)):
                    stats_aio = asyncio.run(
                        download_all_voucher_images(records_to_download=records, concurrency_limit=2, force=False)
                    )
            self.assertEqual(stats_aio["total"], 2)
            self.assertEqual(stats_aio["skipped"], 1)
            self.assertEqual(stats_aio["success"], 1)

        r03 = (PROJECT_ROOT / "scripts/morphometrics/03_fourier_extractor.R").read_text(encoding="utf-8")
        rfpu = (PROJECT_ROOT / "scripts/morphometrics/fourier_parallel_utils.R").read_text(encoding="utf-8")
        rdlu = (PROJECT_ROOT / "scripts/morphometrics/dissected_leaf_utils.R").read_text(encoding="utf-8")
        rdb = (PROJECT_ROOT / "scripts/core/db_utils.R").read_text(encoding="utf-8")
        r06 = (PROJECT_ROOT / "scripts/analysis/06_triage_dashboard_synthesis.R").read_text(encoding="utf-8")

        self.assertIn('make_option(c("-r", "--report-out")', r03)
        self.assertIn('make_option(c("-g", "--plot-out")', r03)
        self.assertNotIn("ucs_v", rfpu)
        self.assertNotIn("ucs_m", rfpu)
        self.assertIn('pull_num("archetype_score")', rfpu)
        self.assertIn("spec_df$archetype_score <- res$arch", rfpu)
        self.assertIn("log10(res$cs)", rfpu)
        self.assertIn('area_mm2 <- if ("area_mm2" %in% names(fac_df))', rdlu)
        self.assertIn('if ("efa_qc_pass" %in% names(df))', rdb)
        self.assertIn('if ("qc_flag_reason" %in% names(df))', rdb)
        self.assertIn("dplyr::coalesce(n_lvs_raw, leaf_cnt_raw, bc_sum, bas_cnt_raw)", r06)




class TestSteps1To4ZeroTrustVariableContracts(unittest.TestCase):
    """Zero-trust regression suite verifying Steps 1-4 variable, geographic, anatomical, and schema contracts."""

    def test_steps1_to_4_zero_trust_variable_and_schema_contracts(self):
        import cv2
        import numpy as np
        import pandas as pd
        from scripts.core.config import PipelineConfig
        from scripts.core.db import get_db_connection
        from scripts.core.db_schema import PACKERA_SCHEMA_USER_VERSION, init_db_schema
        from scripts.core.epa_ecoregions import DEFAULT_CEC_NA_L3_GEOJSON, infer_regional_group, _normalize_us_county_name
        from scripts.core.usda_plants_ranges import _normalize_us_state_name, load_us_county_spatial_graph, resolve_voucher_us_county
        from scripts.core.taxonomy_authority import stratify_duplicates
        import importlib
        compute_capitulum_metrics = importlib.import_module('scripts.pipeline.02_segment_and_extract').compute_capitulum_metrics

        root = Path(__file__).resolve().parents[2]
        cfg = PipelineConfig.from_yaml(root / "config" / "config.yaml")

        # 1. Step 1: western_longitude_threshold == -104.05, Canadian CEC ecoregions & provinces, VA independent cities
        self.assertAlmostEqual(float(cfg.harvesting.western_longitude_threshold), -104.05, places=4)
        self.assertIn("na_cec_eco_l3.geojson", str(DEFAULT_CEC_NA_L3_GEOJSON))
        self.assertEqual(_normalize_us_state_name("ON"), "ontario")
        self.assertEqual(_normalize_us_state_name("Québec"), "quebec")
        self.assertEqual(_normalize_us_county_name("Fairfax City", "virginia"), "fairfax city")
        self.assertEqual(_normalize_us_county_name("Fairfax County", "virginia"), "fairfax")
        self.assertEqual(_normalize_us_county_name("City of Baltimore", "maryland"), "baltimore city")
        self.assertEqual(_normalize_us_county_name("St. Louis City", "missouri"), "st louis city")
        self.assertEqual(_normalize_us_county_name("Chesapeake", "virginia"), "chesapeake")
        self.assertEqual(_normalize_us_county_name("Norfolk", "virginia"), "norfolk")

        sg_ca = load_us_county_spatial_graph()
        st_ca, co_ca, is_non_us = resolve_voucher_us_county(45.42, -75.69, "Ontario", "Ottawa", sg_ca)
        self.assertEqual(st_ca, "ontario")
        self.assertFalse(is_non_us)
        eco_on = infer_regional_group(43.70, -79.40, state_province="Ontario", county="Toronto")
        self.assertTrue(bool(eco_on) and eco_on != "Unmapped_Ecoregion")

        # Missing pre-GPS coordinateUncertainty receives neutral 1000 m weight (never 5000 m penalty)
        dup_df = pd.DataFrame([
            {"catalogNumber": "DUP_NO_UNC", "institutionCode": "GH", "scientificName": "Packera aurea",
             "recordedBy": "M.L. Fernald", "recordNumber": "401", "eventDate": "1912-06-10",
             "determiner_tier": "Tier_1_Gold", "decimalLatitude": 45.1, "decimalLongitude": -68.9,
             "coordinateUncertainty": np.nan},
            {"catalogNumber": "DUP_2500M", "institutionCode": "NY", "scientificName": "Packera aurea",
             "recordedBy": "M.L. Fernald", "recordNumber": "401", "eventDate": "1912-06-10",
             "determiner_tier": "Tier_1_Gold", "decimalLatitude": 45.1, "decimalLongitude": -68.9,
             "coordinateUncertainty": 2500.0},
        ])
        strat = stratify_duplicates(dup_df, home_herbarium="NCU")
        prim_cat = strat.loc[strat["is_primary_duplicate"], "catalogNumber"].iloc[0]
        self.assertEqual(prim_cat, "DUP_NO_UNC", "Neutral 1000 m missing uncertainty must outrank coarse 2500 m sheet.")

        # 2. Step 2: Single-head capitulum HSV phyllary isolation & [0.25, 4.5] bounds
        crop_bgr = np.full((240, 200, 3), 245, dtype=np.uint8)  # pale herbarium sheet background
        mask_cap = np.zeros((240, 200), dtype=np.uint8)
        # Green herbaceous involucre (H ~ 60 px, W ~ 50 px) at bottom + white capillary pappus plume at top (H ~ 80 px, W ~ 110 px)
        cv2.ellipse(crop_bgr, (100, 155), (25, 30), 0, 0, 360, (45, 95, 55), -1)   # BGR olive-green phyllaries
        cv2.ellipse(crop_bgr, (100, 85), (55, 42), 0, 0, 360, (238, 240, 242), -1)  # BGR bright white pappus
        cv2.ellipse(mask_cap, (100, 155), (25, 30), 0, 0, 360, 255, -1)
        cv2.ellipse(mask_cap, (100, 85), (55, 42), 0, 0, 360, 255, -1)
        cap_m = compute_capitulum_metrics(mask_cap, (0, 0, 240, 200), pixels_per_mm=10.0, crop_bgr=crop_bgr, class_id=2)
        self.assertIsNotNone(cap_m)
        # Isolated phyllary cup width should reflect ~50 px (~5.0 mm), not the 110 px (11.0 mm) white pappus flare!
        self.assertLess(cap_m["involucre_width_mm"], 7.5)
        self.assertGreaterEqual(cap_m["AR_inv"], 0.90)
        # Pressed campanulate head with AR_inv = 0.55 survives new [0.25, 4.5] bounds
        mask_camp = np.zeros((200, 200), dtype=np.uint8)
        cv2.ellipse(mask_camp, (100, 100), (55, 30), 0, 0, 360, 255, -1)
        camp_m = compute_capitulum_metrics(mask_camp, (40, 30, 160, 170), pixels_per_mm=10.0, class_id=4)
        self.assertIsNotNone(camp_m)
        self.assertLess(camp_m["AR_inv"], 0.70)
        self.assertGreaterEqual(camp_m["AR_inv"], 0.25)

        # 3. SQLite schema user_version == 20261114 & zero twin columns across tables and v_integrative_triage_master
        self.assertEqual(PACKERA_SCHEMA_USER_VERSION, 20261114)
        conn = get_db_connection(":memory:", init_schema=True)
        try:
            init_db_schema(conn)
            v_cols = {r[1] for r in conn.execute("PRAGMA table_info(vouchers)").fetchall()}
            el_cols = {r[1] for r in conn.execute("PRAGMA table_info(extracted_leaves)").fetchall()}
            otu_cols = {r[1] for r in conn.execute("PRAGMA table_info(otu_morphometrics)").fetchall()}
            view_cols = {r[1] for r in conn.execute("PRAGMA table_info(v_integrative_triage_master)").fetchall()}
            for twin in ("latitude", "longitude"):
                self.assertNotIn(twin, v_cols)
                self.assertNotIn(twin, otu_cols)
                self.assertNotIn(twin, view_cols)
            self.assertIn("within_rosette_efa_dispersion", otu_cols)
            self.assertIn("morphometric_outlier_flag", otu_cols)
            self.assertIn("misidentification_flag", otu_cols)
            self.assertIn("morphometric_outlier_flag", view_cols)
            self.assertIn("misidentification_flag", view_cols)
            self.assertIn("area_mm2", el_cols)
            self.assertIn("aspect_ratio", el_cols)
            for k in range(1, 25):
                self.assertNotIn(f"B{k}", otu_cols)
                self.assertNotIn(f"C{k}", otu_cols)
        finally:
            conn.close()

        # 4. Step 3 & Step 4 R anti-circularity & leakage-free contracts
        aw_src = (root / "scripts" / "morphometrics" / "authority_weighting_utils.R").read_text(encoding="utf-8")
        pm_src = (root / "scripts" / "morphometrics" / "plasticity_mitigation_utils.R").read_text(encoding="utf-8")
        pf_src = (root / "scripts" / "morphometrics" / "pub_figures_utils.R").read_text(encoding="utf-8")
        s04_src = (root / "scripts" / "morphometrics" / "04_gmm_morphotools.R").read_text(encoding="utf-8")
        self.assertIn('(G >= 2L) && (target_modality != "morph")', aw_src)
        self.assertIn("allow_supervised_f = FALSE", aw_src)
        self.assertIn("allow_supervised_f = TRUE", aw_src)
        self.assertNotIn("decimalLatitude", pm_src)
        self.assertNotIn("regional_group", pm_src)
        self.assertIn("kern_fine + 0.35 * kern_reg", pf_src)
        self.assertIn('c("AR_inv", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio")', s04_src)


if __name__ == "__main__":
    unittest.main()
