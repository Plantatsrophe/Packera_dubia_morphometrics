#!/usr/bin/env python3
"""
Step 02 scientific hygiene and vision invariants test suite:
  1. Shape-decoupled anatomical leaf_zone classification.
  2. Tightened cluster_plant_individuals rosette reach and mask bridging.
  3. Millimeter-calibrated CIELAB lamina photometry invariance and delta b*.
  4. Verified archival ruler vs. unboxed margin tick-lattice provenance.
  5. Corymbose *_MANY head counting and nodding/aged-ochre capitulum phenology.
  6. Harmonized EFA/ECT sinus-bridge opening and save_debug_crops gating.
"""

from __future__ import annotations

import importlib
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

SegmentAndExtractPipeline = importlib.import_module(
    "scripts.pipeline.02_segment_and_extract"
).SegmentAndExtractPipeline
from scripts.vision.capitulum_phenology_classifier import (
    HEAD_ANTHESIS,
    HEAD_BUD,
    classify_single_capitulum,
)
from scripts.vision.leaf_geometry_utils import export_standardized_contour
from scripts.vision.lm3_inference_engines import (
    DetectedInstance,
    LM3ScaleCalibrator,
    _ensure_ect_module_loaded,
    compute_leaf_ect_descriptor,
)
from scripts.vision.organ_crown_metrics import (
    classify_leaf_zone_anatomical,
    cluster_plant_individuals,
)
from scripts.vision.step02_photometry_corymb_utils import (
    compute_calibrated_blade_photometry,
    harmonize_upright_mask_sinus_bridges,
)


def _draw_mm_ruler_strip(
    canvas: np.ndarray,
    x_start: int,
    y_top: int,
    px_per_mm: float,
    n_mm: int,
    thickness: int = 2,
) -> None:
    """Draw a synthetic 1-5-10 mm hierarchical metric tick strip onto a sheet or crop."""
    for mm in range(n_mm):
        x = int(round(x_start + mm * px_per_mm))
        if 0 <= x < canvas.shape[1]:
            tick_h = 38 if mm % 10 == 0 else (26 if mm % 5 == 0 else 16)
            cv2.line(canvas, (x, y_top), (x, min(canvas.shape[0] - 1, y_top + tick_h)), (20, 20, 20), thickness)


def _make_flowering_head_crop(corolla_hsv: tuple[int, int, int]) -> np.ndarray:
    """Build a 120x120 RGB capitulum crop with an erect green cup and specified upper disc HSV."""
    img = np.full((120, 120, 3), 240, dtype=np.uint8)
    cv2.rectangle(img, (40, 60), (80, 108), color=(35, 95, 30), thickness=-1)
    hsv_px = np.uint8([[list(corolla_hsv)]])
    rgb_col = tuple(int(c) for c in cv2.cvtColor(hsv_px, cv2.COLOR_HSV2RGB)[0, 0])
    cv2.circle(img, (60, 38), 28, color=rgb_col, thickness=-1)
    return img


def _run_mocked_voucher(
    tmp_dir: Path,
    catalog_number: str,
    sheet: np.ndarray,
    stage1_dets: list,
    save_debug_crops: bool = False,
) -> tuple[object, list[Path]]:
    """Run SegmentAndExtractPipeline.process_voucher with mocked Stage 1/2 predictions."""
    img_p = tmp_dir / f"{catalog_number}.jpg"
    img_p.write_bytes(b"\xff\xd8\xff" + b"0" * 128)
    v_csv = tmp_dir / f"{catalog_number}_vouchers.csv"
    pd.DataFrame([{"catalogNumber": catalog_number, "image_path": str(img_p)}]).to_csv(v_csv, index=False)

    pipe = SegmentAndExtractPipeline(
        vouchers_csv=v_csv,
        model_weights=Path("models/dummy.pth"),
        output_dir=tmp_dir / "data",
        device="cpu",
        force=True,
        save_debug_crops=save_debug_crops,
    )

    def _fake_s2(crop_img: np.ndarray):
        ch, cw = crop_img.shape[:2]
        m = np.zeros((ch, cw), dtype=np.uint8)
        cv2.ellipse(m, (cw // 2, ch // 2), (int(cw * 0.26), int(ch * 0.40)), 0, 0, 360, 255, -1)
        return [((8, 8, ch - 8, cw - 8), m, 0.96, 0)]

    try:
        with patch("cv2.imread", return_value=sheet):
            with patch.object(pipe.stage1_detector, "predict", return_value=stage1_dets):
                with patch.object(pipe.model_engine, "predict", side_effect=_fake_s2):
                    bundle = pipe.process_voucher(catalog_number, img_p, _return_bundle=True)
        patches = list(pipe.capitula_patches_dir.glob("*.jpg"))
        return bundle, patches
    finally:
        pipe.shutdown()


class TestStep02VisionHygiene(unittest.TestCase):
    """Focused verification of the 6 Step 02 scientific and visual hygiene invariants."""

    def test_shape_decoupled_leaf_zone_anatomical(self) -> None:
        """Linear rosette leaves classify basal while subpetiolate lower-cauline and unrooted stem leaves classify cauline."""
        caps, roots = [(200, 480, 250, 520)], [(2000, 480, 2080, 520)]
        leaves = [
            {
                "leaf_id": 1, "plant_individual_id": 0,
                "p_petiole_tip": (495.0, 1995.0), "p_base_sheet": (480.0, 1980.0), "p_apex_sheet": (415.0, 1915.0),
                "blade_length_px": 92.0, "aspect_ratio": 4.6, "area_px": 450.0,
                "has_petiole_mask": True, "petiole_length_px": 21.0,
            },
            {
                "leaf_id": 2, "plant_individual_id": 0,
                "p_petiole_tip": (505.0, 1995.0), "p_base_sheet": (520.0, 1980.0), "p_apex_sheet": (585.0, 1915.0),
                "blade_length_px": 92.0, "aspect_ratio": 4.6, "area_px": 450.0,
                "has_petiole_mask": True, "petiole_length_px": 21.0,
            },
            {
                "leaf_id": 3, "plant_individual_id": 0,
                "p_petiole_tip": (500.0, 2010.0), "p_base_sheet": (500.0, 2030.0), "p_apex_sheet": (500.0, 2120.0),
                "blade_length_px": 90.0, "aspect_ratio": 4.6, "area_px": 450.0,
                "has_petiole_mask": True, "petiole_length_px": 20.0,
            },
            {
                "leaf_id": 4, "plant_individual_id": 0,
                "p_petiole_tip": (500.0, 1432.0), "p_base_sheet": (518.0, 1420.0), "p_apex_sheet": (590.0, 1365.0),
                "blade_length_px": 90.0, "aspect_ratio": 3.4, "area_px": 780.0,
                "has_petiole_mask": True, "petiole_length_px": 22.0,
            },
        ]
        zones = classify_leaf_zone_anatomical(
            leaves, crowns_by_pid={0: (500.0, 2000.0)}, capitula_detections=caps, root_boxes=roots
        )
        for lid in (1, 2, 3):
            self.assertEqual(zones[lid]["leaf_zone"], "basal")
        self.assertEqual(zones[4]["leaf_zone"], "cauline")

        sessile_stem = [
            {
                "leaf_id": 10 + i, "plant_individual_id": 0,
                "p_base_sheet": (500.0, float(y_pos)), "p_apex_sheet": (565.0, float(y_pos - 60)),
                "blade_length_px": 88.0, "aspect_ratio": 2.8, "has_petiole_mask": False,
            }
            for i, y_pos in enumerate((950, 720, 500))
        ]
        zones_sessile = classify_leaf_zone_anatomical(sessile_stem, capitula_detections=caps, root_boxes=[])
        for lid in (10, 11, 12):
            self.assertEqual(zones_sessile[lid]["leaf_zone"], "cauline")

    def test_tightened_cluster_plant_individuals_and_mask_bridge(self) -> None:
        """Isolated packet leaves 120 mm away stay distinct unless bridged by specimen_mask."""
        dummy = np.ones((12, 12), dtype=np.uint8)
        b1 = DetectedInstance("PKT1", 1, (2400, 400, 2500, 500), dummy, 0.95, 0)
        b2 = DetectedInstance("PKT1", 2, (2420, 460, 2520, 560), dummy, 0.94, 0)
        pkt_cal = DetectedInstance("PKT1", 3, (2410, 1630, 2510, 1730), dummy, 0.91, 0)

        out_cal = cluster_plant_individuals([b1, b2, pkt_cal], 2600, 3000, pixels_per_mm=10.0)
        self.assertEqual([inst.plant_individual_id for inst in out_cal], [0, 0, 1])

        pkt_uncal = DetectedInstance("PKT2", 3, (2410, 1270, 2510, 1370), dummy, 0.91, 0)
        out_uncal = cluster_plant_individuals([b1, b2, pkt_uncal], 2600, 3000, pixels_per_mm=None)
        self.assertEqual([inst.plant_individual_id for inst in out_uncal], [0, 0, 1])

        spec_mask = np.zeros((3000, 2600), dtype=np.uint8)
        cv2.line(spec_mask, (480, 2460), (1680, 2460), 255, thickness=16)
        out_bridged = cluster_plant_individuals(
            [b1, b2, pkt_cal], 2600, 3000, specimen_mask=spec_mask, pixels_per_mm=10.0
        )
        self.assertEqual([inst.plant_individual_id for inst in out_bridged], [0, 0, 0])

    def test_physical_mm_photometry_invariance_and_b_delta(self) -> None:
        """Calibrated photometry stays invariant across 10 vs 25 px/mm and detects tomentum plus yellowing."""
        rng = np.random.default_rng(7)
        sheet_10 = np.full((320, 480, 3), (236, 236, 236), dtype=np.uint8)
        mask_10 = np.zeros((160, 160), dtype=np.uint8)
        cv2.circle(mask_10, (80, 80), 62, 255, -1)

        sheet_glab = sheet_10.copy()
        glab_tex = np.clip(64.0 + rng.normal(0.0, 1.8, size=(160, 160, 1)), 40, 90).astype(np.uint8)
        glab_bgr = np.concatenate([glab_tex - 16, glab_tex + 8, glab_tex + 4], axis=2)
        sheet_glab[80:240, 40:200][mask_10 > 0] = glab_bgr[mask_10 > 0]

        sheet_tom = sheet_10.copy()
        wool_hi = rng.normal(0.0, 28.0, size=(100, 100)).astype(np.float32)
        wool_10 = cv2.resize(wool_hi, (160, 160), interpolation=cv2.INTER_AREA)
        tom_lum = np.clip(170.0 + wool_10, 92.0, 224.0).astype(np.uint8)
        sheet_tom[80:240, 40:200][mask_10 > 0] = np.stack([tom_lum, tom_lum, tom_lum], axis=2)[mask_10 > 0]

        sheet_brown = sheet_10.copy()
        sheet_brown[80:240, 40:200][mask_10 > 0] = (55, 138, 178)

        m_glab_10 = compute_calibrated_blade_photometry(sheet_glab, mask_10, (40, 80), pixels_per_mm=10.0)
        m_tom_10 = compute_calibrated_blade_photometry(sheet_tom, mask_10, (40, 80), pixels_per_mm=10.0)
        m_brown_10 = compute_calibrated_blade_photometry(sheet_brown, mask_10, (40, 80), pixels_per_mm=10.0)

        self.assertGreater(m_tom_10["blade_cie_l_norm"], m_glab_10["blade_cie_l_norm"] + 0.25)
        self.assertGreater(m_tom_10["blade_cie_l_sd"], m_glab_10["blade_cie_l_sd"] * 2.5)
        self.assertGreaterEqual(m_brown_10["blade_cie_b_delta"], 8.0)

        mask_25 = cv2.resize(mask_10, None, fx=2.5, fy=2.5, interpolation=cv2.INTER_NEAREST)
        sheet_glab_25 = cv2.resize(sheet_glab, None, fx=2.5, fy=2.5, interpolation=cv2.INTER_LINEAR)
        sheet_tom_25 = cv2.resize(sheet_tom, None, fx=2.5, fy=2.5, interpolation=cv2.INTER_LINEAR)
        m_glab_25 = compute_calibrated_blade_photometry(sheet_glab_25, mask_25, (100, 200), pixels_per_mm=25.0)
        m_tom_25 = compute_calibrated_blade_photometry(sheet_tom_25, mask_25, (100, 200), pixels_per_mm=25.0)

        self.assertLessEqual(abs(m_glab_10["blade_cie_l_norm"] - m_glab_25["blade_cie_l_norm"]), 0.05)
        self.assertLessEqual(abs(m_tom_10["blade_cie_l_norm"] - m_tom_25["blade_cie_l_norm"]), 0.05)
        max_tom_sd = max(m_tom_10["blade_cie_l_sd"], m_tom_25["blade_cie_l_sd"], 1e-6)
        self.assertLessEqual(abs(m_tom_10["blade_cie_l_sd"] - m_tom_25["blade_cie_l_sd"]) / max_tom_sd, 0.25)

    def test_ruler_calibration_unboxed_vs_verified_crop(self) -> None:
        """Unboxed margin tick lattices stay uncalibrated while verified ruler boxes calibrate cropped sheets."""
        cal = LM3ScaleCalibrator(device="cpu", providers=[], use_cf_predicted_by_mp=False)
        cal.archival_session = None

        sheet_unboxed = np.full((1600, 2000, 3), 240, dtype=np.uint8)
        cv2.rectangle(sheet_unboxed, (100, 30), (1800, 130), (225, 225, 225), -1)
        _draw_mm_ruler_strip(sheet_unboxed, x_start=120, y_top=30, px_per_mm=5.0, n_mm=300, thickness=2)

        res_unb = cal.calibrate_sheet(sheet_unboxed, catalog_number="NCU_UNBOXED", use_cf_predicted_by_mp=False)
        self.assertEqual(res_unb["cf_source"], "tick_lattice_unboxed")
        self.assertEqual(res_unb["scale_source"], "tick_lattice_unboxed")
        self.assertFalse(res_unb["is_calibrated_mm"])
        self.assertAlmostEqual(float(res_unb["pixels_per_mm"]), 5.0, delta=0.15)

        sheet_crop = np.full((3400, 4000, 3), 240, dtype=np.uint8)
        mp_est = (cal.MP_PRIOR_SLOPE_PX_PER_CM * math.sqrt((3400.0 * 4000.0) / 1e6)) / 10.0
        true_ppm = 0.60 * mp_est
        ruler_bbox = (40, 200, 180, 1100)
        cv2.rectangle(sheet_crop, (300, 60), (960, 150), (225, 225, 225), -1)
        _draw_mm_ruler_strip(sheet_crop, x_start=320, y_top=60, px_per_mm=true_ppm, n_mm=95, thickness=2)

        cal.archival_session, cal.ruler_classifier, cal.ruler_lattice = object(), None, None
        with patch.object(cal, "_detect_ruler_boxes", return_value=[(ruler_bbox, 0.90)]):
            res_cal = cal.calibrate_sheet(sheet_crop, catalog_number="NCU_VERIFIED", use_cf_predicted_by_mp=False)

        self.assertTrue(res_cal["is_calibrated_mm"])
        self.assertEqual(res_cal["scale_source"], "ruler_calibrated")
        self.assertAlmostEqual(float(res_cal["pixels_per_mm"]), true_ppm, delta=0.20)

    def test_corymb_many_head_counting_and_aged_nodding_phenology(self) -> None:
        """Flower_MANY clusters increment counts without fake involucre dimensions, and nodding/ochre heads score anthesis."""
        sheet = np.full((1000, 800, 3), 235, dtype=np.uint8)
        sheet[60:220, 260:460] = (30, 175, 215)
        s1_dets = [((320, 220, 560, 420), 0.95, 0), ((60, 260, 220, 460), 0.91, 5)]
        with tempfile.TemporaryDirectory() as tmp_d:
            bundle, _ = _run_mocked_voucher(Path(tmp_d), "NCU_MANY", sheet, s1_dets, save_debug_crops=False)

        self.assertIsNotNone(bundle.repro_record)
        self.assertIsNotNone(bundle.pheno_record)
        self.assertGreaterEqual(int(bundle.repro_record["capitula_count"]), 2)
        self.assertGreaterEqual(int(bundle.pheno_record["total_capitula"]), 2)
        self.assertTrue(pd.isna(bundle.repro_record["involucre_height_mm"]))
        self.assertTrue(pd.isna(bundle.repro_record["involucre_width_mm"]))

        nodding_rgb = cv2.rotate(_make_flowering_head_crop((15, 153, 217)), cv2.ROTATE_180)
        self.assertEqual(classify_single_capitulum(nodding_rgb)["head_state"], HEAD_ANTHESIS)

        ochre_rgb = _make_flowering_head_crop((14, 61, 179))
        res_ochre = classify_single_capitulum(ochre_rgb, class_id=4)
        self.assertEqual(res_ochre["head_state"], HEAD_ANTHESIS)
        self.assertTrue(res_ochre["has_aged_ochre_corollas"])
        res_bud = classify_single_capitulum(ochre_rgb, class_id=6)
        self.assertEqual(res_bud["head_state"], HEAD_BUD)
        self.assertFalse(res_bud["has_aged_ochre_corollas"])

    def test_sinus_bridge_efa_ect_parity_and_debug_crop_gate(self) -> None:
        """Harmonized sinus opening feeds 120-vertex contours and 32x360 ECT, while save_debug_crops=False skips crop JPEGs."""
        mask_u8 = np.zeros((160, 160), dtype=np.uint8)
        cv2.ellipse(mask_u8, (80, 80), (45, 58), 0, 0, 360, 255, -1)
        mask_u8[68:92, 82:128] = 0
        mask_u8[68:92, 118:120] = 255

        cleaned_u8, evf = harmonize_upright_mask_sinus_bridges(mask_u8, pixels_per_mm=10.0)
        self.assertLess(evf, 0.01)

        with tempfile.TemporaryDirectory() as tmp_d:
            tmp_p = Path(tmp_d)
            csv_path, coords = export_standardized_contour(
                cleaned_u8, "TEST_SINUS", 1, tmp_p / "data",
                return_coords=True, already_oriented=True, bridge_already_opened=True,
            )
            self.assertTrue(Path(csv_path).is_file())
            self.assertEqual(len(coords), 120)
            arr_xy = np.asarray(coords, dtype=np.float64)
            self.assertTrue(np.any((arr_xy[:, 0] >= 80.0) & (arr_xy[:, 0] <= 116.0) & (arr_xy[:, 1] >= 68.0) & (arr_xy[:, 1] <= 92.0)))

            _ensure_ect_module_loaded()
            if "ect" in sys.modules:
                ect_clean = compute_leaf_ect_descriptor(cleaned_u8, num_dirs=360, num_thresh=32, already_oriented=True)
                ect_orig = compute_leaf_ect_descriptor(mask_u8, num_dirs=360, num_thresh=32, already_oriented=True)
                self.assertIsNotNone(ect_clean)
                self.assertEqual(ect_clean.shape, (32, 360))
                self.assertFalse(np.allclose(ect_clean, ect_orig))

            sheet = np.full((1000, 800, 3), 235, dtype=np.uint8)
            sheet[60:180, 310:410] = (20, 200, 240)
            s1_dets = [((320, 220, 560, 420), 0.95, 0), ((60, 310, 180, 410), 0.92, 4)]
            _, patches = _run_mocked_voucher(tmp_p, "NCU_NOCROP", sheet, s1_dets, save_debug_crops=False)
            self.assertEqual(patches, [])


if __name__ == "__main__":
    unittest.main()
