#!/usr/bin/env python3
"""
===============================================================================
Script: test_geometry_gatekeeper_enhancements.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit test suite verifying geometric gatekeeper enhancements:
      1. Straight-chord fold detection (detect_folded_leaf) and bilateral synthesis (Tier 2).
      2. Irregular and diagonal fold rejection (Failed QC).
      3. Convexity defect inspection for botanical dissection (is_botanical_dissection) (Tier 1).
      4. Asymmetric foreign leaf occlusion rejection from botanical dissection.
      5. Diagnostic helper routing (classify_and_route_leaf and LeafRoutingResult).
===============================================================================
"""

from __future__ import annotations

import math
import unittest
from pathlib import Path
from typing import Tuple

import cv2
import numpy as np

import sys
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.lm2_geometry_utils import (
    classify_and_route_leaf,
    compute_geometric_metrics,
    detect_folded_leaf,
    is_botanical_dissection,
    reflect_folded_mask,
    resample_normalized_contour,
    LeafRoutingResult,
)


def create_synthetic_folded_leaf(
    canvas_size: Tuple[int, int] = (400, 400),
    center: Tuple[int, int] = (200, 200),
    axes: Tuple[int, int] = (50, 150),
    angle_deg: float = 0.0
) -> np.ndarray:
    """
    Creates a synthetic folded leaf: an ellipse bisected along its longitudinal midrib,
    yielding one perfectly straight boundary (chord) and one curved margin.
    """
    mask = np.zeros(canvas_size, dtype=np.uint8)
    # Draw half-ellipse from angle 270 to 90 (or 0 to 180)
    cv2.ellipse(mask, center, axes, angle_deg, -90, 90, 255, -1)
    return mask


def create_synthetic_dissected_leaf(
    canvas_size: Tuple[int, int] = (400, 400),
    center: Tuple[int, int] = (200, 200),
    length: int = 240,
    max_width: int = 120
) -> np.ndarray:
    """
    Creates a synthetic pinnatifid/lyrate leaf with multiple (>= 4) bilateral sinuses
    and lateral lobes alternating along the blade margin.
    Solidity falls cleanly into the 0.48 <= Solidity < 0.72 botanical dissection window.
    """
    mask = np.zeros(canvas_size, dtype=np.uint8)
    cx, cy = center
    half_l = length // 2

    # Draw main elliptical terminal lobe
    cv2.ellipse(mask, (cx, cy - half_l + 50), (max_width // 3, 50), 0, 0, 360, 255, -1)
    # Draw longitudinal rachis/midrib stem
    cv2.line(mask, (cx, cy - half_l), (cx, cy + half_l), 255, thickness=18)

    # Draw 4 pairs of bilateral lateral lobes
    lobe_y_offsets = [-30, 10, 50, 90]
    for i, y_off in enumerate(lobe_y_offsets):
        lobe_w = int(max_width * (0.45 - 0.05 * i))
        lobe_h = 16
        # Left lobe
        cv2.ellipse(mask, (cx - lobe_w, cy + y_off), (lobe_w, lobe_h), -15, 0, 360, 255, -1)
        # Right lobe
        cv2.ellipse(mask, (cx + lobe_w, cy + y_off), (lobe_w, lobe_h), 15, 0, 360, 255, -1)

    # Smooth the synthetic outline slightly
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    return mask


def create_synthetic_occluded_leaf(
    canvas_size: Tuple[int, int] = (400, 400),
    center: Tuple[int, int] = (200, 200),
    axes: Tuple[int, int] = (70, 150)
) -> np.ndarray:
    """
    Creates an elliptical leaf with a single foreign leaf cutout / occlusion
    on only one lateral margin.
    """
    mask = np.zeros(canvas_size, dtype=np.uint8)
    cv2.ellipse(mask, center, axes, 0, 0, 360, 255, -1)
    # Foreign occlusion cutout on left side only
    cv2.circle(mask, (center[0] - 60, center[1] + 20), 45, 0, -1)
    return mask


class TestGeometryGatekeeperEnhancements(unittest.TestCase):
    """Test cases for folded leaf detection, dissection inspection, and routing."""

    def test_detect_folded_leaf_straight_midrib(self) -> None:
        """Verify that a leaf folded along its midrib is detected and reflected."""
        # Vertical half-ellipse: straight chord along x = 200
        mask = create_synthetic_folded_leaf(canvas_size=(400, 400), center=(200, 200), axes=(40, 140), angle_deg=0.0)

        # Before reflection, solidity is high (> 0.88) because half-ellipse is convex!
        ucs, solidity, _, _, _ = compute_geometric_metrics(mask)
        self.assertGreaterEqual(solidity, 0.85, "Half-ellipse has high solidity despite being folded.")

        res = detect_folded_leaf(mask)
        self.assertTrue(res["is_folded"], "detect_folded_leaf must identify straight chord fold.")
        self.assertFalse(res["is_irregular_fold"])
        self.assertIsNotNone(res["straight_side"])
        self.assertLessEqual(res["wl_ratio"], 0.45)
        self.assertIsNotNone(res["synthesized_mask"])

        # Check synthesized mask
        synth = res["synthesized_mask"]
        synth_ucs, synth_solidity, _, _, _ = compute_geometric_metrics(synth)
        self.assertGreater(np.count_nonzero(synth), np.count_nonzero(mask) * 1.5,
                           "Synthesized mask must reconstruct bilateral hemi-blade.")
        self.assertGreaterEqual(synth_solidity, 0.72)

    def test_classify_and_route_folded_leaf_to_tier2(self) -> None:
        """Verify that a folded leaf routes to Tier 2 (not Tier 1) despite high solidity."""
        mask = create_synthetic_folded_leaf(canvas_size=(400, 400), center=(200, 200), axes=(35, 130), angle_deg=0.0)

        result: LeafRoutingResult = classify_and_route_leaf(mask)
        self.assertTrue(result.is_folded)
        self.assertEqual(result.tier, "Tier 2", "Folded leaf must route to Tier 2 after bilateral reflection.")
        self.assertEqual(result.assigned_tier, "tier2")
        self.assertIsNotNone(result.processed_mask)
        self.assertIsNotNone(result.contour_coords)

    def test_unfolded_whole_leaf_routes_to_tier1(self) -> None:
        """Verify that an intact full leaf (not folded) routes directly to Tier 1."""
        mask = np.zeros((400, 400), dtype=np.uint8)
        # Full ellipse with W/L > 0.45
        cv2.ellipse(mask, (200, 200), (80, 140), 0, 0, 360, 255, -1)

        result: LeafRoutingResult = classify_and_route_leaf(mask)
        self.assertFalse(result.is_folded)
        self.assertFalse(result.is_irregular_fold)
        self.assertEqual(result.tier, "Tier 1")
        self.assertEqual(result.assigned_tier, "tier1")

    def test_irregular_diagonal_fold_routes_to_failed_qc(self) -> None:
        """Verify that a leaf folded diagonally or irregularly routes to Failed QC."""
        # Create a narrow mask with a prominent diagonal fold crease
        mask = np.zeros((400, 400), dtype=np.uint8)
        # Draw a polygon representing a diagonally folded blade
        pts = np.array([
            [190, 60],   # Apex
            [230, 180],  # Diagonal crease endpoint 1
            [200, 320],  # Base
            [160, 200],  # Curved side point
        ], dtype=np.int32)
        cv2.fillPoly(mask, [pts], 255)

        fold_res = detect_folded_leaf(mask)
        self.assertTrue(fold_res["is_irregular_fold"], "Diagonal crease must trigger is_irregular_fold.")
        self.assertFalse(fold_res["is_folded"])

        result: LeafRoutingResult = classify_and_route_leaf(mask)
        self.assertEqual(result.tier, "Failed QC")
        self.assertEqual(result.assigned_tier, "rejected")
        self.assertTrue(result.is_irregular_fold)

    def test_botanical_dissection_detection_and_routing(self) -> None:
        """Verify that a dissected basal leaf with bilateral sinuses passes to Tier 1."""
        mask = create_synthetic_dissected_leaf(canvas_size=(400, 400), center=(200, 200), length=260, max_width=130)

        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(mask)
        # Confirm that solidity falls in the typical dissection window [0.48, 0.72)
        self.assertGreaterEqual(solidity, 0.48, f"Dissected leaf solidity {solidity} should be >= 0.48")
        self.assertLess(solidity, 0.72, f"Dissected leaf solidity {solidity} should be < 0.72")

        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        main_cnt = max(contours, key=cv2.contourArea)

        is_dissected = is_botanical_dissection(main_cnt, solidity, p_apex, p_base)
        self.assertTrue(is_dissected, "is_botanical_dissection must recognize bilateral sinuses.")

        # Classify and route
        result: LeafRoutingResult = classify_and_route_leaf(mask)
        self.assertIn(result.tier, ["Tier 1", "Tier 1 (Dissected)"], "Botanical dissection leaf must be rescued and routed to Tier 1.")
        self.assertEqual(result.assigned_tier, "tier1")

    def test_foreign_occlusion_rejected_from_dissection(self) -> None:
        """Verify that asymmetric unilateral occlusion is NOT classified as botanical dissection."""
        mask = create_synthetic_occluded_leaf(canvas_size=(400, 400), center=(200, 200), axes=(65, 145))

        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(mask)
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        main_cnt = max(contours, key=cv2.contourArea)

        # Unilateral defect: should NOT be botanical dissection
        is_dissected = is_botanical_dissection(main_cnt, solidity, p_apex, p_base)
        self.assertFalse(is_dissected, "Unilateral occlusion must not be misclassified as botanical dissection.")

    def test_empty_and_degenerate_mask_handling(self) -> None:
        """Verify numerical stability on empty or tiny masks."""
        empty_mask = np.zeros((100, 100), dtype=np.uint8)
        result = classify_and_route_leaf(empty_mask)
        self.assertEqual(result.tier, "Failed QC")
        self.assertEqual(result.assigned_tier, "rejected")

        # Tiny dot
        tiny_mask = np.zeros((50, 50), dtype=np.uint8)
        tiny_mask[24:26, 24:26] = 255
        res_tiny = classify_and_route_leaf(tiny_mask)
        self.assertEqual(res_tiny.tier, "Failed QC")

    def test_resample_normalized_contour(self) -> None:
        """Verify resampled contour normalization coordinates."""
        mask = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(mask, (100, 100), 40, 255, -1)

        coords = resample_normalized_contour(mask, num_points=120)
        self.assertIsNotNone(coords)
        self.assertEqual(coords.shape, (120, 2))
        self.assertGreaterEqual(np.min(coords), 0.0)
        self.assertLessEqual(np.max(coords), 1.0)


if __name__ == "__main__":
    unittest.main()
