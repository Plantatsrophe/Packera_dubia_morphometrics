#!/usr/bin/env python3
"""
scripts/tests/test_geometry_edge_cases.py
=========================================
Unit test suite verifying edge-case geometric morphologies in Packera dubia:
  1. Fixture 1: Normal Unfolded Oval Leaf (W/L ~ 0.65, Solidity ~ 0.88, crenate teeth)
  2. Fixture 2: Leaf Folded in Half Along Midrib (Solidity > 0.85, W/L ~ 0.32, reflected ~2x area)
  3. Fixture 3: Highly Dissected / Lyrate Leaf (Ovate terminal blade, 4 bilateral sinuses, Solidity ~ 0.58)
  4. Fixture 4: Damaged Leaf with Foreign Leaf Occlusion (Solidity ~ 0.60, non-botanical unilateral defect)
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path
from typing import Tuple

import cv2
import numpy as np

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.lm2_geometry_utils import (
    classify_and_route_leaf,
    compute_geometric_metrics,
    detect_folded_leaf,
    is_botanical_dissection,
    trim_petiole_tail,
    FoldDetectionResult,
    LeafRoutingResult,
)


def create_fixture_1_unfolded_oval_leaf(
    canvas_size: Tuple[int, int] = (400, 400),
    center: Tuple[int, int] = (200, 200),
    length: int = 240,
    wl_ratio: float = 0.65,
) -> np.ndarray:
    """
    Fixture 1: Normal Unfolded Oval Leaf.
    Generates a synthetic binary mask of a standard ovate leaf with small crenate teeth
    (W/L = 0.65, Solidity ~ 0.88).
    """
    ry = length / 2.0
    rx = (length * wl_ratio) / 2.0
    thetas = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    num_teeth = 28
    amp = 9.0
    tooth = - amp * (np.sin(num_teeth * thetas) ** 2)

    y_norm = np.sin(thetas)
    shape = 1.0 + 0.06 * y_norm

    r_x = (rx + tooth) * shape
    r_y = ry + tooth * (ry / rx)
    xs = center[0] + r_x * np.cos(thetas)
    ys = center[1] + r_y * np.sin(thetas)
    pts = np.column_stack([xs, ys]).astype(np.int32)

    mask = np.zeros(canvas_size, dtype=np.uint8)
    cv2.fillPoly(mask, [pts], 255)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    return mask


def create_fixture_2_folded_midrib_leaf(
    f1_mask: np.ndarray,
    split_x: int = 200,
) -> np.ndarray:
    """
    Fixture 2: Leaf Folded in Half Along Midrib.
    Slices Fixture 1 in half vertically along the midline x = split_x, keeping one half.
    Has Solidity > 0.85 and W/L ~ 0.32.
    """
    f2 = f1_mask.copy()
    f2[:, :split_x] = 0
    return f2


def create_fixture_3_dissected_lyrate_leaf(
    canvas_size: Tuple[int, int] = (400, 400),
    center: Tuple[int, int] = (200, 200),
    length: int = 260,
    max_width: int = 130,
) -> np.ndarray:
    """
    Fixture 3: Highly Dissected / Lyrate Leaf.
    Generates a synthetic leaf with an ovate terminal blade and 4 deep, alternating
    bilateral sinuses (Solidity ~ 0.58).
    """
    mask = np.zeros(canvas_size, dtype=np.uint8)
    cx, cy = center
    half_l = length // 2

    # Draw ovate terminal blade
    cv2.ellipse(mask, (cx, cy - half_l + 55), (max_width // 3, 45), 0, 0, 360, 255, -1)
    # Draw central longitudinal midrib / rachis
    cv2.line(mask, (cx, cy - half_l), (cx, cy + half_l), 255, thickness=18)

    # 4 deep, alternating bilateral sinuses / lobes
    lobes = [
        (-1, -30, 46, 20),
        (1, 2, 48, 20),
        (-1, 38, 44, 20),
        (1, 74, 38, 20),
    ]
    for side, y_off, lobe_w, lh in lobes:
        lobe_cx = cx + side * (lobe_w + 5)
        cv2.ellipse(mask, (lobe_cx, cy + y_off), (lobe_w, lh), side * 15, 0, 360, 255, -1)

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    return mask


def create_fixture_4_occluded_damaged_leaf(
    f1_mask: np.ndarray,
    center: Tuple[int, int] = (200, 200),
) -> np.ndarray:
    """
    Fixture 4: Damaged Leaf with Foreign Leaf Occlusion.
    Takes Fixture 1 and subtracts an irregular, deep elliptical bite on one side (Solidity ~ 0.60).
    """
    mask = f1_mask.copy()
    cx, cy = center
    cv2.ellipse(mask, (cx - 50, cy - 10), (58, 70), -20, 0, 360, 0, -1)
    return mask


def create_fixture_5_leaf_with_petiole_stalk(
    f1_mask: np.ndarray,
    stalk_length: int = 45,
    stalk_width: int = 10,
) -> np.ndarray:
    """
    Fixture 5: Leaf with Narrow Slender Petiole Stalk.
    Takes Fixture 1 (ovate leaf terminating at y ~ 320) and appends a narrow linear
    petiole stalk (W = 10 px, length = 45 px) extending below the lamina down to y = 365.
    """
    f5_mask = f1_mask.copy()
    half_w = stalk_width // 2
    cv2.rectangle(f5_mask, (200 - half_w, 315), (200 + half_w, 315 + stalk_length), 255, -1)
    return f5_mask


def create_fixture_6_cordate_leaf(
    canvas_size: Tuple[int, int] = (400, 400),
    center: Tuple[int, int] = (200, 200),
    length: int = 240,
    max_width: int = 160,
) -> np.ndarray:
    """
    Fixture 6: Broad Cordate Leaf with Basal Lobes.
    Generates a synthetic cordate leaf mask where the blade base is broad (W_base > 50 px)
    with bilateral lobes, testing the cordate/truncate protection guard.
    """
    mask = np.zeros(canvas_size, dtype=np.uint8)
    cx, cy = center
    thetas = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    r_x = (max_width / 2.0) * np.sin(thetas)
    r_y = - (length / 2.0) * (np.cos(thetas) - 0.2 * np.sin(thetas) ** 2)
    lobe_notch = 15.0 * np.exp(-((thetas - np.pi) ** 2) / 0.15)
    r_y += lobe_notch
    xs = cx + r_x
    ys = cy + r_y
    pts = np.column_stack([xs, ys]).astype(np.int32)
    cv2.fillPoly(mask, [pts], 255)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)
    return mask


class TestGeometryEdgeCases(unittest.TestCase):
    """Synthetic unit tests verifying folded and dissected edge-case morphologies."""

    def setUp(self) -> None:
        self.f1_mask = create_fixture_1_unfolded_oval_leaf()

    def test_fixture_1_normal_unfolded_oval_leaf(self) -> None:
        """Fixture 1: Normal unfolded leaf (W/L ~ 0.65, Solidity ~ 0.88) routes to Tier 1."""
        mask = self.f1_mask
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(mask)
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        self.assertTrue(len(contours) > 0, "Fixture 1 mask must produce a valid external contour.")
        cnt = max(contours, key=cv2.contourArea)

        rect = cv2.minAreaRect(cnt.astype(np.float32))
        wl_ratio = float(min(rect[1]) / max(rect[1]))

        # Geometric parameter validations
        self.assertAlmostEqual(wl_ratio, 0.65, delta=0.05, msg="W/L ratio must be approximately 0.65.")
        self.assertAlmostEqual(solidity, 0.88, delta=0.05, msg="Solidity must be approximately 0.88.")

        # Fold detection assertion: must return False
        fold_res = detect_folded_leaf(mask)
        self.assertFalse(detect_folded_leaf(mask), "detect_folded_leaf(mask) must return False for unfolded leaf.")
        self.assertFalse(fold_res["is_folded"], "fold_res['is_folded'] must be False.")
        self.assertFalse(fold_res["is_irregular_fold"], "fold_res['is_irregular_fold'] must be False.")

        # Routing assertion: must route to Tier 1
        routing = classify_and_route_leaf(mask)
        self.assertEqual(routing.tier, "Tier 1", "Unfolded normal leaf must route to Tier 1.")
        self.assertEqual(routing.assigned_tier, "tier1")
        self.assertFalse(routing.is_folded)
        self.assertFalse(routing.is_dissected)

    def test_fixture_2_folded_in_half_along_midrib(self) -> None:
        """Fixture 2: Leaf folded along midrib (Solidity > 0.85, W/L ~ 0.32) reflects and routes to Tier 2."""
        f2_mask = create_fixture_2_folded_midrib_leaf(self.f1_mask, split_x=200)
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(f2_mask)
        contours, _ = cv2.findContours(f2_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        self.assertTrue(len(contours) > 0, "Fixture 2 mask must produce a valid external contour.")
        cnt = max(contours, key=cv2.contourArea)

        rect = cv2.minAreaRect(cnt.astype(np.float32))
        wl_ratio = float(min(rect[1]) / max(rect[1]))

        # Half-leaf geometric properties
        self.assertGreater(solidity, 0.85, f"Folded half-leaf solidity {solidity:.3f} must be > 0.85.")
        self.assertAlmostEqual(wl_ratio, 0.32, delta=0.05, msg="Folded half-leaf W/L must be approximately 0.32.")

        # Fold detection assertion: must return True
        fold_res = detect_folded_leaf(f2_mask)
        self.assertTrue(detect_folded_leaf(f2_mask), "detect_folded_leaf(mask) must return True for folded leaf.")
        self.assertTrue(fold_res["is_folded"], "fold_res['is_folded'] must be True.")
        self.assertFalse(fold_res["is_irregular_fold"], "fold_res['is_irregular_fold'] must be False.")
        self.assertIsNotNone(fold_res["straight_side"], "Folded leaf must identify straight side.")

        # Reflected synthesized mask area assertion: ~2x input area
        synth_mask = fold_res["synthesized_mask"]
        self.assertIsNotNone(synth_mask, "Synthesized mask must be generated for folded leaf.")
        orig_area = float(np.count_nonzero(f2_mask))
        synth_area = float(np.count_nonzero(synth_mask))
        area_ratio = synth_area / orig_area
        self.assertAlmostEqual(
            area_ratio, 2.0, delta=0.25,
            msg=f"Synthesized reflected mask must restore bilateral area (~2x input area, got {area_ratio:.2f}x)."
        )

        # Routing assertion: must route to Tier 2
        routing = classify_and_route_leaf(f2_mask)
        self.assertEqual(routing.tier, "Tier 2", "Folded leaf must route to Tier 2.")
        self.assertEqual(routing.assigned_tier, "tier2")
        self.assertTrue(routing.is_folded)
        self.assertIsNotNone(routing.processed_mask)

    def test_fixture_3_highly_dissected_lyrate_leaf(self) -> None:
        """Fixture 3: Dissected lyrate leaf (Solidity ~ 0.58) accepted into Tier 1 rather than rejected."""
        f3_mask = create_fixture_3_dissected_lyrate_leaf()
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(f3_mask)
        contours, _ = cv2.findContours(f3_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        self.assertTrue(len(contours) > 0, "Fixture 3 mask must produce a valid external contour.")
        cnt = max(contours, key=cv2.contourArea)

        # Dissected solidity validation
        self.assertAlmostEqual(solidity, 0.58, delta=0.07, msg=f"Solidity {solidity:.3f} must be approximately 0.58.")

        # Dissection check assertion: is_botanical_dissection(contour) returns True
        is_dissected = is_botanical_dissection(cnt)
        self.assertTrue(is_dissected, "is_botanical_dissection(contour) must return True for lyrate leaf.")

        # Routing assertion: accepts leaf into Tier 1 rather than rejecting for low solidity
        routing = classify_and_route_leaf(f3_mask)
        self.assertIn(
            routing.tier, ["Tier 1", "Tier 1 (Dissected)"],
            f"Lyrate dissected leaf must be accepted into Tier 1, got {routing.tier}."
        )
        self.assertEqual(routing.assigned_tier, "tier1")
        self.assertTrue(routing.is_dissected)

    def test_fixture_4_damaged_leaf_with_foreign_occlusion(self) -> None:
        """Fixture 4: Occluded/damaged leaf (Solidity ~ 0.60) rejected from dissection and fails Tier 1."""
        f4_mask = create_fixture_4_occluded_damaged_leaf(self.f1_mask)
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(f4_mask)
        contours, _ = cv2.findContours(f4_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        self.assertTrue(len(contours) > 0, "Fixture 4 mask must produce a valid external contour.")
        cnt = max(contours, key=cv2.contourArea)

        # Low solidity from unilateral defect
        self.assertAlmostEqual(solidity, 0.60, delta=0.07, msg=f"Solidity {solidity:.3f} must be approximately 0.60.")

        # Dissection check assertion: is_botanical_dissection(contour) returns False
        is_dissected = is_botanical_dissection(cnt)
        self.assertFalse(
            is_dissected,
            "is_botanical_dissection(contour) must return False for unilateral foreign occlusion."
        )

        # Routing assertion: must NOT pass Tier 1 (routes to Tier 2 reflection if clean half, or fails QC)
        routing = classify_and_route_leaf(f4_mask)
        self.assertNotIn(
            routing.tier, ["Tier 1", "Tier 1 (Dissected)"],
            "Occluded/damaged leaf must NOT pass Tier 1."
        )
        self.assertNotEqual(routing.assigned_tier, "tier1")
        self.assertIn(
            routing.tier, ["Tier 2", "Failed QC"],
            f"Occluded leaf must route to Tier 2 reflection or Failed QC, got {routing.tier}."
        )

    def test_fixture_5_petiole_stalk_trimming_and_restoration(self) -> None:
        """Fixture 5: Slender petiole stalk trimmed at junction, restoring blade metrics and Tier 1 routing."""
        f5_mask = create_fixture_5_leaf_with_petiole_stalk(self.f1_mask, stalk_length=45, stalk_width=10)
        orig_area = np.count_nonzero(self.f1_mask)
        stalk_area = np.count_nonzero(f5_mask)
        self.assertGreater(stalk_area, orig_area, "Fixture 5 mask must include attached stalk pixels.")

        # Metric degradation before trimming
        ucs_stalk, sol_stalk, _, _, _ = compute_geometric_metrics(f5_mask)
        self.assertLess(sol_stalk, 0.89, f"Untrimmed stalk should degrade solidity (got {sol_stalk:.3f}).")

        # Trim petiole tail
        trimmed_mask, meta = trim_petiole_tail(f5_mask, return_metadata=True)
        self.assertTrue(meta["is_trimmed"], "trim_petiole_tail must detect and trim petiole stalk.")
        self.assertEqual(meta["reason"], "WIDTH_INFLECTION_JUNCTION_TRIMMED")
        self.assertGreaterEqual(meta["trim_distance_px"], 35.0, "Trim distance must span the narrow stalk length.")
        self.assertIsNotNone(meta["junction_pt"])

        # Metric restoration after trimming
        trimmed_area = np.count_nonzero(trimmed_mask)
        area_diff = abs(trimmed_area - orig_area)
        self.assertLessEqual(area_diff, 15, f"Trimmed mask must match intact blade area within 15 px (diff={area_diff}).")

        ucs_trimmed, sol_trimmed, _, _, _ = compute_geometric_metrics(trimmed_mask)
        self.assertGreaterEqual(sol_trimmed, 0.90, f"Trimmed mask solidity ({sol_trimmed:.3f}) must be restored to >= 0.90.")
        self.assertGreaterEqual(ucs_trimmed, 0.95, f"Trimmed mask UCS ({ucs_trimmed:.3f}) must be restored to >= 0.95.")

        # Tier routing integration assertion: classify_and_route_leaf automatically trims petiole and routes to Tier 1
        routing = classify_and_route_leaf(f5_mask)
        self.assertEqual(routing.tier, "Tier 1", "Leaf with trimmed petiole tail must route to Tier 1.")
        self.assertEqual(routing.assigned_tier, "tier1")
        self.assertFalse(routing.is_folded)
        self.assertTrue(routing.metadata.get("petiole_trim_info", {}).get("is_trimmed", False))

    def test_fixture_6_cordate_blade_protection(self) -> None:
        """Fixture 6: Broad cordate blade base is guarded and 100% protected against false clipping."""
        f6_mask = create_fixture_6_cordate_leaf()
        orig_area = np.count_nonzero(f6_mask)

        trimmed_mask, meta = trim_petiole_tail(f6_mask, return_metadata=True)
        self.assertFalse(meta["is_trimmed"], "Cordate leaf base must NOT be clipped.")
        self.assertEqual(meta["reason"], "CORDATE_OR_BROAD_BASE_PROTECTED")
        self.assertEqual(np.count_nonzero(trimmed_mask), orig_area, "Cordate leaf area must be 100% preserved.")

    def test_clean_unfolded_leaf_left_unclipped(self) -> None:
        """Verify that a cleanly clipped blade without petiole stalk passes fallback check and remains unclipped."""
        f1_mask = self.f1_mask
        orig_area = np.count_nonzero(f1_mask)

        trimmed_mask, meta = trim_petiole_tail(f1_mask, return_metadata=True)
        self.assertFalse(meta["is_trimmed"], "Cleanly clipped blade must NOT be clipped.")
        self.assertEqual(meta["reason"], "CLEAN_BLADE_UNCLIPPED")
        self.assertEqual(np.count_nonzero(trimmed_mask), orig_area, "Intact blade area must be 100% preserved.")

    def test_fallback_morphological_linear_protrusion(self) -> None:
        """Verify that fallback morphological check catches and trims narrow linear protrusions (aspect ratio > 3:1)."""
        # Create an oval blade with a slender protrusion (length 25 px, width 6 px, aspect ratio ~ 4:1)
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 130), (50, 90), 0, 0, 360, 255, -1)
        # Add narrow protrusion at basal tip (y from 220 to 245, width 6)
        cv2.rectangle(mask, (147, 220), (153, 245), 255, -1)

        trimmed_mask, meta = trim_petiole_tail(mask, return_metadata=True)
        self.assertTrue(meta["is_trimmed"], "Narrow linear protrusion must be trimmed.")
        self.assertIn(
            meta["reason"],
            ["WIDTH_INFLECTION_JUNCTION_TRIMMED", "FALLBACK_MORPHOLOGICAL_PROTRUSION_TRIMMED"]
        )
        self.assertLess(np.count_nonzero(trimmed_mask), np.count_nonzero(mask))


if __name__ == "__main__":
    unittest.main()
