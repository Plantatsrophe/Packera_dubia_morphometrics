#!/usr/bin/env python3
"""
===============================================================================
Script: test_capitulum_phenology.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration test suite validating:
    1. Biomarker A: Capillary pappus plume extraction (white fibrous bristles vs.
       flat mounting sheet paper vs. vegetative structures).
    2. Biomarker B: Active yellow corollas/ligules extraction (anthesis indicator).
    3. Biomarker C: Involucre posture profiling (erect vs. reflexed phyllaries).
    4. Head-level decision logic (HEAD_ANTHESIS, HEAD_FRUIT, HEAD_BUD).
    5. Voucher-level Botanical Determinate Cyme Rule (heterogeneity resolution).
    6. Performance constraint verification (< 0.05 seconds per crop).
    7. Schema compliance of exported manifest.

Execution:
    python -m unittest scripts/tests/test_capitulum_phenology.py
===============================================================================
"""

from __future__ import annotations

import sys
import time
import unittest
from pathlib import Path
from typing import Dict, List

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
    classify_single_capitulum,
    classify_voucher_phenology,
    extract_biomarker_involucre_posture,
    extract_biomarker_pappus_plume,
    extract_biomarker_yellow_corollas,
)


def create_synthetic_fruiting_capitulum(h: int = 120, w: int = 80) -> np.ndarray:
    """
    Synthesizes an RGB crop of a fruiting capitulum:
    - Upper 45%: white fibrous pappus plume with fine bristle texture.
    - Lower 55%: brownish/greenish reflexed or flaring involucre.
    """
    img = np.full((h, w, 3), 240, dtype=np.uint8)  # White/neutral paper background

    # Upper 45%: White fibrous bristles on neutral paper
    upper_h = int(h * 0.45)
    for y in range(upper_h):
        for x in range(10, w - 10):
            # Alternating fine fibrous strands (high L*, low S, high texture)
            if (x + y) % 3 == 0 or (x * 2 + y) % 5 == 0:
                img[y, x] = [245, 245, 240]
            else:
                img[y, x] = [200, 195, 185]

    # Lower 55%: Greenish/brown reflexed involucre
    cv2.ellipse(
        img,
        center=(w // 2, int(h * 0.75)),
        axes=(w // 2 - 8, int((h - upper_h) * 0.35)),
        angle=0,
        startAngle=0,
        endAngle=180,
        color=(60, 80, 40),
        thickness=-1,
    )
    return img


def create_synthetic_anthesis_capitulum(h: int = 120, w: int = 80) -> np.ndarray:
    """
    Synthesizes an RGB crop of an active flowering capitulum:
    - Upper & middle: saturated golden-yellow corollas (disc + ray florets).
    - Lower: erect green cylindrical involucre.
    """
    img = np.full((h, w, 3), 235, dtype=np.uint8)

    # Golden yellow disc & ray ligules in upper/middle region (H_deg ~ 38°, in 18°-45° band)
    cv2.circle(img, (w // 2, int(h * 0.35)), radius=28, color=(255, 165, 0), thickness=-1)
    # Add radiating golden ray petals
    for angle in range(0, 360, 30):
        rad = np.deg2rad(angle)
        x2 = int(w // 2 + 35 * np.cos(rad))
        y2 = int(h * 0.35 + 25 * np.sin(rad))
        cv2.line(img, (w // 2, int(h * 0.35)), (x2, y2), color=(255, 150, 0), thickness=4)

    # Erect green cylindrical involucre in lower half (H_inv >= W_inv)
    cv2.rectangle(img, (w // 2 - 16, int(h * 0.45)), (w // 2 + 16, h - 5), color=(40, 95, 30), thickness=-1)
    return img


def create_synthetic_bud_capitulum(h: int = 120, w: int = 80) -> np.ndarray:
    """
    Synthesizes an RGB crop of a compact vegetative bud:
    - Compact, erect green tapered involucre.
    - No yellow ligules and no white fibrous bristles.
    """
    img = np.full((h, w, 3), 240, dtype=np.uint8)
    # Erect tapered green bud: H_inv > W_inv
    pts = np.array([
        [w // 2, int(h * 0.30)],
        [w // 2 + 15, int(h * 0.60)],
        [w // 2 + 12, h - 10],
        [w // 2 - 12, h - 10],
        [w // 2 - 15, int(h * 0.60)],
    ], dtype=np.int32)
    cv2.fillPoly(img, [pts], color=(35, 90, 35))
    return img


class TestCapitulumPhenologyClassifier(unittest.TestCase):
    """Test suite for head-level phenological biomarkers and voucher aggregation."""

    def test_biomarker_pappus_plume_detection(self) -> None:
        """Validates that fibrous white bristles trigger plume detection while flat paper does not."""
        fruiting_crop = create_synthetic_fruiting_capitulum(120, 80)
        has_plume, plume_ratio, diag = extract_biomarker_pappus_plume(fruiting_crop)

        self.assertTrue(has_plume, f"Expected pappus plume on fruiting crop, got {diag}")
        self.assertGreater(plume_ratio, 0.18, "Plume ratio must exceed 18%")
        self.assertGreater(diag["laplacian_var"], 20.0, "Fibrous plume must exhibit high-frequency texture")

        # Flat white background mounting sheet (no fibrous texture)
        flat_white_paper = np.full((120, 80, 3), 245, dtype=np.uint8)
        has_paper_plume, paper_ratio, paper_diag = extract_biomarker_pappus_plume(flat_white_paper)

        self.assertFalse(has_paper_plume, "Flat white paper must not trigger pappus plume false positive")
        self.assertLess(paper_diag["laplacian_var"], 10.0, "Flat paper should have near-zero Laplacian variance")

    def test_biomarker_yellow_corollas_detection(self) -> None:
        """Validates isolation of the 18°-45° saturated yellow corolla band."""
        anthesis_crop = create_synthetic_anthesis_capitulum(120, 80)
        has_yellow, yellow_ratio, diag = extract_biomarker_yellow_corollas(anthesis_crop)

        self.assertTrue(has_yellow, f"Anthesis crop must detect yellow corollas: {diag}")
        self.assertGreater(yellow_ratio, 0.12, "Yellow corolla ratio must exceed 12%")

        bud_crop = create_synthetic_bud_capitulum(120, 80)
        has_bud_yellow, bud_yellow_ratio, bud_diag = extract_biomarker_yellow_corollas(bud_crop)
        self.assertFalse(has_bud_yellow, f"Vegetative bud must not exhibit fresh yellow corollas: {bud_diag}")
        self.assertLess(bud_yellow_ratio, 0.05)

    def test_biomarker_involucre_posture(self) -> None:
        """Validates profile classification into erect cylindrical vs. reflexed/flared."""
        bud_crop = create_synthetic_bud_capitulum(120, 80)
        erect, reflexed, diag = extract_biomarker_involucre_posture(bud_crop)

        self.assertTrue(erect, f"Bud involucre must be erect: {diag}")
        self.assertFalse(reflexed, f"Bud involucre should not be reflexed: {diag}")
        self.assertGreaterEqual(diag["involucre_aspect_ratio"], 0.8)

    def test_head_level_decision_logic(self) -> None:
        """Validates classify_single_capitulum on standard developmental stages."""
        fruiting_crop = create_synthetic_fruiting_capitulum(120, 80)
        res_fruit = classify_single_capitulum(fruiting_crop)
        self.assertEqual(res_fruit["head_state"], HEAD_FRUIT)
        self.assertTrue(res_fruit["has_pappus_plume"])

        anthesis_crop = create_synthetic_anthesis_capitulum(120, 80)
        res_anthesis = classify_single_capitulum(anthesis_crop)
        self.assertEqual(res_anthesis["head_state"], HEAD_ANTHESIS)
        self.assertTrue(res_anthesis["has_fresh_yellow_corollas"])

        bud_crop = create_synthetic_bud_capitulum(120, 80)
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
        sample_crop = create_synthetic_fruiting_capitulum(200, 150)

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
