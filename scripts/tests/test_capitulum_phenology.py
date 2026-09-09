#!/usr/bin/env python3
"""
===============================================================================
Script: test_capitulum_phenology.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Comprehensive unit and integration test suite validating:
    1. Synthetic Head Image Fixtures:
       - Fixture A: Fresh Flowering Head (Anthesis; H/W=1.2 green cup, H=30°, S=0.60, V=0.85 yellow disc)
       - Fixture B: Fruiting Head (Brown reflexed base H/W=0.70, fibrous dome L*=85, S=0.08 with edge noise)
       - Fixture C: Compact Green Bud (Tight green oval, no yellow or white pixels)
    2. Voucher-Level Botanical Determinate Cyme Rule:
       - Case 1 (Mixed Cyme): 4 Anthesis + 2 Fruit -> "anthesis"
       - Case 2 (All Fruit Cyme): 0 Anthesis + 5 Fruit -> "fruit"
       - Case 3 (Sterile Rosette): 0 capitula -> "sterile"
    3. Biomarker extractions (pappus plume, yellow corollas, involucre posture).
    4. Guardrails for empty/degenerate crops and latency budget (< 0.05s / crop).

Execution:
    python -m unittest scripts/tests/test_capitulum_phenology.py
===============================================================================
"""

from __future__ import annotations

import sys
import time
import unittest
from pathlib import Path
from typing import Any, Dict, List

import cv2
import numpy as np
import pandas as pd

# Ensure repository root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.capitulum_phenology_classifier import (
    HEAD_ANTHESIS,
    HEAD_BUD,
    HEAD_FRUIT,
    HEAD_UNKNOWN,
    VOUCHER_ANTHESIS,
    VOUCHER_BUD,
    VOUCHER_FRUIT,
    VOUCHER_STERILE,
    VOUCHER_UNKNOWN,
    classify_single_capitulum,
    classify_voucher_phenology,
    extract_biomarker_involucre_posture,
    extract_biomarker_pappus_plume,
    extract_biomarker_yellow_corollas,
)


def create_fixture_a_fresh_flowering_head(h: int = 120, w: int = 120) -> np.ndarray:
    """
    Fixture A: Fresh Flowering Head (Anthesis).

    Synthesizes an in-memory RGB NumPy array with:
    - Neutral mounting paper background (RGB: 240, 240, 240).
    - Green lower involucre cup with aspect ratio H/W = 1.2 (e.g., W=40, H=48 px).
    - Saturated golden-yellow upper disc (H = 30°, S = 0.60, V = 0.85).
    """
    img = np.full((h, w, 3), 240, dtype=np.uint8)

    # Green lower cup (H/W = 1.2)
    cup_w = max(4, int(round(w * 0.333)))
    cup_h = int(round(cup_w * 1.2))
    x1 = (w - cup_w) // 2
    x2 = x1 + cup_w
    y1 = int(round(h * 0.50))
    y2 = min(h - 2, y1 + cup_h)
    cv2.rectangle(img, (x1, y1), (x2, y2), color=(35, 95, 30), thickness=-1)

    # Saturated yellow upper disc: H = 30°, S = 0.60, V = 0.85
    # In OpenCV HSV: H_cv2 = 30/2 = 15, S_cv2 = round(0.60 * 255) = 153, V_cv2 = round(0.85 * 255) = 217
    hsv_col = np.uint8([[[15, int(round(0.60 * 255.0)), int(round(0.85 * 255.0))]]])
    rgb_col = tuple(int(c) for c in cv2.cvtColor(hsv_col, cv2.COLOR_HSV2RGB)[0, 0])
    disc_radius = max(3, int(round(min(h, w) * 0.233)))
    cv2.circle(img, (w // 2, int(round(h * 0.32))), disc_radius, color=rgb_col, thickness=-1)

    return img


def create_fixture_b_fruiting_head(h: int = 120, w: int = 120) -> np.ndarray:
    """
    Fixture B: Fruiting Head with Expanded Pappus Plume.

    Synthesizes an in-memory RGB NumPy array with:
    - Neutral mounting paper background (RGB: 240, 240, 240).
    - Brown reflexed base with aspect ratio H/W = 0.70 (e.g., W=50, H=35 px) and flared/reflexed phyllaries.
    - Bright white, fibrous upper dome (L* = 85, S = 0.08) with synthetic edge noise.
    """
    img = np.full((h, w, 3), 240, dtype=np.uint8)

    # Brown reflexed base with H/W = 0.70 (e.g. W=50, H=35 px)
    base_w = max(6, int(round(w * 0.417)))
    target_h = int(round(base_w * 0.70))
    tooth_h = max(2, int(round(target_h * 0.15)))
    rect_h = target_h - tooth_h

    bx1 = (w - base_w) // 2
    bx2 = bx1 + base_w
    by1 = int(round(h * 0.54))
    by2 = by1 + rect_h
    cv2.rectangle(img, (bx1, by1), (bx2, by2), color=(70, 48, 25), thickness=-1)

    # Add reflexed jagged teeth along basal perimeter to enforce boundary roughness >= 1.15
    for i in range(bx1 + 2, bx2 - 4, 4):
        tooth = np.array([[i, by2], [i + 2, by2 + tooth_h], [i + 4, by2]], dtype=np.int32)
        cv2.fillPoly(img, [tooth], color=(65, 42, 20))

    # Bright white, fibrous upper dome (L* = 85, S = 0.08) with synthetic edge noise
    center_x = w / 2.0
    rad_x = w * 0.30
    rad_y = h * 0.33
    y_start = max(1, int(round(h * 0.08)))
    y_end = min(h - 1, int(round(h * 0.44)))
    x_start = max(1, int(round(w * 0.20)))
    x_end = min(w - 1, int(round(w * 0.80)))

    for y in range(y_start, y_end + 1):
        for x in range(x_start, x_end + 1):
            if ((x - center_x) / rad_x) ** 2 + ((y - h * 0.42) / rad_y) ** 2 <= 1.0:
                noise = 25 if (x + y) % 2 == 0 else -25
                v_val = np.clip(int(round(0.85 * 255.0)) + noise, 0, 255)
                hsv_px = np.uint8([[[0, int(round(0.08 * 255.0)), v_val]]])
                rgb_px = cv2.cvtColor(hsv_px, cv2.COLOR_HSV2RGB)[0, 0]
                img[y, x] = rgb_px

    return img


def create_fixture_c_compact_green_bud(h: int = 120, w: int = 120) -> np.ndarray:
    """
    Fixture C: Compact Green Bud.

    Synthesizes an in-memory RGB NumPy array with:
    - Neutral mounting paper background (RGB: 240, 240, 240).
    - Solid green, tight, small oval without yellow or white pixels.
    """
    img = np.full((h, w, 3), 240, dtype=np.uint8)
    ax_x = max(2, int(round(w * 0.133)))
    ax_y = max(4, int(round(h * 0.217)))
    cv2.ellipse(img, (w // 2, h // 2), (ax_x, ax_y), 0, 0, 360, color=(35, 90, 30), thickness=-1)
    return img


# Backwards compatibility aliases
create_synthetic_anthesis_capitulum = create_fixture_a_fresh_flowering_head
create_synthetic_fruiting_capitulum = create_fixture_b_fruiting_head
create_synthetic_bud_capitulum = create_fixture_c_compact_green_bud


class TestCapitulumPhenologyClassifier(unittest.TestCase):
    """Test suite for head-level phenological biomarkers, fixtures, and voucher aggregation."""

    # -------------------------------------------------------------------------
    # 1. Synthetic Head Image Fixtures Tests
    # -------------------------------------------------------------------------
    def test_fixture_a_fresh_flowering_head(self) -> None:
        """Fixture A: Fresh Flowering Head (Anthesis) asserts classify_single_capitulum returns HEAD_ANTHESIS."""
        crop = create_fixture_a_fresh_flowering_head(120, 120)
        self.assertEqual(crop.shape, (120, 120, 3))

        res = classify_single_capitulum(crop)
        self.assertEqual(res["head_state"], HEAD_ANTHESIS)
        self.assertEqual(res, HEAD_ANTHESIS)
        self.assertTrue(res["has_fresh_yellow_corollas"])
        self.assertTrue(res["erect_involucre"])
        self.assertFalse(res["has_pappus_plume"])

    def test_fixture_b_fruiting_head_expanded_pappus_plume(self) -> None:
        """Fixture B: Fruiting Head with Expanded Pappus Plume asserts classify_single_capitulum returns HEAD_FRUIT."""
        crop = create_fixture_b_fruiting_head(120, 120)
        self.assertEqual(crop.shape, (120, 120, 3))

        res = classify_single_capitulum(crop)
        self.assertEqual(res["head_state"], HEAD_FRUIT)
        self.assertEqual(res, HEAD_FRUIT)
        self.assertTrue(res["has_pappus_plume"])
        self.assertFalse(res["has_fresh_yellow_corollas"])

    def test_fixture_c_compact_green_bud(self) -> None:
        """Fixture C: Compact Green Bud asserts classify_single_capitulum returns HEAD_BUD."""
        crop = create_fixture_c_compact_green_bud(120, 120)
        self.assertEqual(crop.shape, (120, 120, 3))

        res = classify_single_capitulum(crop)
        self.assertEqual(res["head_state"], HEAD_BUD)
        self.assertEqual(res, HEAD_BUD)
        self.assertTrue(res["erect_involucre"])
        self.assertFalse(res["has_pappus_plume"])
        self.assertFalse(res["has_fresh_yellow_corollas"])

    # -------------------------------------------------------------------------
    # 2. Voucher-Level Cyme Classification Tests
    # -------------------------------------------------------------------------
    def test_voucher_case_1_mixed_cyme(self) -> None:
        """Case 1 (Mixed Cyme): 4 Anthesis heads + 2 Fruit heads -> voucher state is 'anthesis'."""
        anthesis_heads = [classify_single_capitulum(create_fixture_a_fresh_flowering_head()) for _ in range(4)]
        fruit_heads = [classify_single_capitulum(create_fixture_b_fruiting_head()) for _ in range(2)]
        all_heads = anthesis_heads + fruit_heads

        res = classify_voucher_phenology("NCU_TEST_MIXED", all_heads)
        self.assertEqual(res["voucher_phenological_state"], "anthesis")
        self.assertEqual(res, "anthesis")
        self.assertEqual(res["voucher_phenological_state"], VOUCHER_ANTHESIS)
        self.assertEqual(res["n_anthesis"], 4)
        self.assertEqual(res["n_fruit"], 2)
        self.assertEqual(res["total_capitula"], 6)

    def test_voucher_case_2_all_fruit_cyme(self) -> None:
        """Case 2 (All Fruit Cyme): 0 Anthesis heads + 5 Fruit heads -> voucher state is 'fruit'."""
        fruit_heads = [classify_single_capitulum(create_fixture_b_fruiting_head()) for _ in range(5)]

        res = classify_voucher_phenology("NCU_TEST_ALL_FRUIT", fruit_heads)
        self.assertEqual(res["voucher_phenological_state"], "fruit")
        self.assertEqual(res, "fruit")
        self.assertEqual(res["voucher_phenological_state"], VOUCHER_FRUIT)
        self.assertEqual(res["n_anthesis"], 0)
        self.assertEqual(res["n_fruit"], 5)
        self.assertEqual(res["total_capitula"], 5)

    def test_voucher_case_3_sterile_rosette(self) -> None:
        """Case 3 (Sterile Rosette): 0 capitula -> voucher state is 'sterile'."""
        res = classify_voucher_phenology("NCU_TEST_STERILE", [])
        self.assertEqual(res["voucher_phenological_state"], "sterile")
        self.assertEqual(res, "sterile")
        self.assertEqual(res["voucher_phenological_state"], VOUCHER_STERILE)
        self.assertEqual(res["total_capitula"], 0)
        self.assertEqual(res["n_anthesis"], 0)
        self.assertEqual(res["n_fruit"], 0)
        self.assertEqual(res["n_bud"], 0)
        self.assertEqual(res["dominant_pappus_ratio"], 0.0)

    # -------------------------------------------------------------------------
    # 3. Biomarker Extractions & Profile Logic Tests
    # -------------------------------------------------------------------------
    def test_biomarker_pappus_plume_detection(self) -> None:
        """Validates that fibrous white bristles trigger plume detection while flat paper does not."""
        fruiting_crop = create_fixture_b_fruiting_head(120, 120)
        has_plume, plume_ratio, diag = extract_biomarker_pappus_plume(fruiting_crop)

        self.assertTrue(has_plume, f"Expected pappus plume on fruiting crop, got {diag}")
        self.assertGreater(plume_ratio, 0.18, "Plume ratio must exceed 18%")
        self.assertGreater(diag["laplacian_var"], 20.0, "Fibrous plume must exhibit high-frequency texture")

        # Flat white background mounting sheet (no fibrous texture)
        flat_white_paper = np.full((120, 120, 3), 245, dtype=np.uint8)
        has_paper_plume, paper_ratio, paper_diag = extract_biomarker_pappus_plume(flat_white_paper)

        self.assertFalse(has_paper_plume, "Flat white paper must not trigger pappus plume false positive")
        self.assertLess(paper_diag["laplacian_var"], 10.0, "Flat paper should have near-zero Laplacian variance")

    def test_biomarker_yellow_corollas_detection(self) -> None:
        """Validates isolation of the 18°-45° saturated yellow corolla band."""
        anthesis_crop = create_fixture_a_fresh_flowering_head(120, 120)
        has_yellow, yellow_ratio, diag = extract_biomarker_yellow_corollas(anthesis_crop)

        self.assertTrue(has_yellow, f"Anthesis crop must detect yellow corollas: {diag}")
        self.assertGreater(yellow_ratio, 0.12, "Yellow corolla ratio must exceed 12%")

        bud_crop = create_fixture_c_compact_green_bud(120, 120)
        has_bud_yellow, bud_yellow_ratio, bud_diag = extract_biomarker_yellow_corollas(bud_crop)
        self.assertFalse(has_bud_yellow, f"Vegetative bud must not exhibit fresh yellow corollas: {bud_diag}")
        self.assertLess(bud_yellow_ratio, 0.05)

    def test_biomarker_involucre_posture(self) -> None:
        """Validates profile classification into erect cylindrical vs. reflexed/flared."""
        bud_crop = create_fixture_c_compact_green_bud(120, 120)
        erect, reflexed, diag = extract_biomarker_involucre_posture(bud_crop)

        self.assertTrue(erect, f"Bud involucre must be erect: {diag}")
        self.assertFalse(reflexed, f"Bud involucre should not be reflexed: {diag}")
        self.assertGreaterEqual(diag["involucre_aspect_ratio"], 0.8)

    def test_head_level_decision_logic(self) -> None:
        """Validates classify_single_capitulum on standard developmental stages."""
        fruiting_crop = create_fixture_b_fruiting_head(120, 120)
        res_fruit = classify_single_capitulum(fruiting_crop)
        self.assertEqual(res_fruit["head_state"], HEAD_FRUIT)
        self.assertTrue(res_fruit["has_pappus_plume"])

        anthesis_crop = create_fixture_a_fresh_flowering_head(120, 120)
        res_anthesis = classify_single_capitulum(anthesis_crop)
        self.assertEqual(res_anthesis["head_state"], HEAD_ANTHESIS)
        self.assertTrue(res_anthesis["has_fresh_yellow_corollas"])

        bud_crop = create_fixture_c_compact_green_bud(120, 120)
        res_bud = classify_single_capitulum(bud_crop)
        self.assertEqual(res_bud["head_state"], HEAD_BUD)
        self.assertTrue(res_bud["erect_involucre"])

    def test_degenerate_and_empty_crops(self) -> None:
        """Validates guardrails against zero-division on empty or sub-pixel crops."""
        empty_crop = np.zeros((0, 0, 3), dtype=np.uint8)
        res_empty = classify_single_capitulum(empty_crop)
        self.assertEqual(res_empty["head_state"], HEAD_UNKNOWN)
        self.assertEqual(res_empty["white_fibrous_ratio"], 0.0)

        tiny_crop = np.full((2, 2, 3), 200, dtype=np.uint8)
        res_tiny = classify_single_capitulum(tiny_crop)
        self.assertEqual(res_tiny["head_state"], HEAD_UNKNOWN)

    def test_voucher_cyme_heterogeneity_rules(self) -> None:
        """Validates Botanical Determinate Cyme Rule across various voucher configurations."""
        cat_num = "NCU00123456"

        # Case 1: At least one anthesis head -> VOUCHER_ANTHESIS
        heads_anthesis = [
            {"head_state": HEAD_ANTHESIS, "white_fibrous_ratio": 0.02},
            {"head_state": HEAD_FRUIT, "white_fibrous_ratio": 0.35},
            {"head_state": HEAD_BUD, "white_fibrous_ratio": 0.01},
        ]
        res1 = classify_voucher_phenology(cat_num, heads_anthesis)
        self.assertEqual(res1["voucher_phenological_state"], VOUCHER_ANTHESIS)
        self.assertEqual(res1["n_anthesis"], 1)
        self.assertEqual(res1["n_fruit"], 1)
        self.assertEqual(res1["n_bud"], 1)
        self.assertEqual(res1["total_capitula"], 3)
        self.assertAlmostEqual(res1["dominant_pappus_ratio"], 0.35, places=4)

        # Case 2: Fruit present, but 0 anthesis -> VOUCHER_FRUIT
        heads_fruit = [
            {"head_state": HEAD_FRUIT, "white_fibrous_ratio": 0.42},
            {"head_state": HEAD_FRUIT, "white_fibrous_ratio": 0.28},
        ]
        res2 = classify_voucher_phenology(cat_num, heads_fruit)
        self.assertEqual(res2["voucher_phenological_state"], VOUCHER_FRUIT)
        self.assertEqual(res2["n_anthesis"], 0)
        self.assertEqual(res2["n_fruit"], 2)
        self.assertAlmostEqual(res2["dominant_pappus_ratio"], 0.42, places=4)

        # Case 3: Only buds -> VOUCHER_BUD
        heads_bud = [
            {"head_state": HEAD_BUD, "white_fibrous_ratio": 0.01},
            {"head_state": HEAD_BUD, "white_fibrous_ratio": 0.02},
        ]
        res3 = classify_voucher_phenology(cat_num, heads_bud)
        self.assertEqual(res3["voucher_phenological_state"], VOUCHER_BUD)
        self.assertEqual(res3["n_bud"], 2)

        # Case 4: No capitula (vegetative rosette) -> VOUCHER_STERILE
        res4 = classify_voucher_phenology(cat_num, [])
        self.assertEqual(res4["voucher_phenological_state"], VOUCHER_STERILE)
        self.assertEqual(res4["total_capitula"], 0)
        self.assertEqual(res4["dominant_pappus_ratio"], 0.0)

    def test_execution_latency_constraint(self) -> None:
        """Validates execution constraint: < 0.05 seconds (50 ms) per crop."""
        sample_crop = create_fixture_b_fruiting_head(200, 150)

        # Warm up
        classify_single_capitulum(sample_crop)

        n_trials = 50
        start_time = time.perf_counter()
        for _ in range(n_trials):
            classify_single_capitulum(sample_crop)
        elapsed = time.perf_counter() - start_time
        avg_time = elapsed / n_trials

        self.assertLess(
            avg_time,
            0.05,
            f"Execution time must be < 0.05s per crop, achieved {avg_time:.4f}s",
        )

    def test_manifest_schema_and_dataframe(self) -> None:
        """Validates that output manifest contains exactly the required columns."""
        records = [
            classify_voucher_phenology("NCU001", [
                {"head_state": HEAD_ANTHESIS, "white_fibrous_ratio": 0.05}
            ]),
            classify_voucher_phenology("NCU002", []),
        ]
        df = pd.DataFrame(records)
        expected_cols = [
            "catalogNumber",
            "total_capitula",
            "n_anthesis",
            "n_fruit",
            "n_bud",
            "voucher_phenological_state",
            "dominant_pappus_ratio",
        ]
        self.assertListEqual(list(df.columns), expected_cols)
        self.assertEqual(len(df), 2)


if __name__ == "__main__":
    unittest.main()
