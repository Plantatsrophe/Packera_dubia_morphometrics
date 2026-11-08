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

from scripts.vision.leaf_geometry_utils import (
    _prepare_lamina_mask,
    approximate_petiole_and_scape_metrics,
    compute_classical_shape_descriptors,
    compute_geometric_metrics,
    estimate_rosette_crown_point,
    is_botanical_dissection,
    resample_normalized_contour,
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
    """Synthetic unit tests verifying geometric descriptors, lamina preparation, and organ metrics."""

    def setUp(self) -> None:
        self.f1_mask = create_fixture_1_unfolded_oval_leaf()

    def test_fixture_1_normal_unfolded_oval_leaf(self) -> None:
        """Fixture 1: Normal unfolded leaf (W/L ~ 0.65, Solidity ~ 0.88) computes expected metrics."""
        mask = self.f1_mask
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(mask)
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        self.assertTrue(len(contours) > 0, "Fixture 1 mask must produce a valid external contour.")
        cnt = max(contours, key=cv2.contourArea)

        rect = cv2.minAreaRect(cnt.astype(np.float32))
        wl_ratio = float(min(rect[1]) / max(rect[1]))

        self.assertAlmostEqual(wl_ratio, 0.65, delta=0.05, msg="W/L ratio must be approximately 0.65.")
        self.assertAlmostEqual(solidity, 0.88, delta=0.05, msg="Solidity must be approximately 0.88.")
        self.assertGreaterEqual(ucs, 0.90)

    def test_fixture_3_highly_dissected_lyrate_leaf(self) -> None:
        """Fixture 3: Dissected lyrate leaf (Solidity ~ 0.58) recognized by is_botanical_dissection."""
        f3_mask = create_fixture_3_dissected_lyrate_leaf()
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(f3_mask)
        contours, _ = cv2.findContours(f3_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        self.assertTrue(len(contours) > 0, "Fixture 3 mask must produce a valid external contour.")
        cnt = max(contours, key=cv2.contourArea)

        self.assertAlmostEqual(solidity, 0.58, delta=0.07, msg=f"Solidity {solidity:.3f} must be approximately 0.58.")
        is_dissected = is_botanical_dissection(cnt)
        self.assertTrue(is_dissected, "is_botanical_dissection(contour) must return True for lyrate leaf.")

    def test_fixture_4_damaged_leaf_with_foreign_occlusion(self) -> None:
        """Fixture 4: Occluded/damaged leaf (Solidity ~ 0.60) rejected from botanical dissection."""
        f4_mask = create_fixture_4_occluded_damaged_leaf(self.f1_mask)
        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(f4_mask)
        contours, _ = cv2.findContours(f4_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        self.assertTrue(len(contours) > 0, "Fixture 4 mask must produce a valid external contour.")
        cnt = max(contours, key=cv2.contourArea)

        self.assertAlmostEqual(solidity, 0.60, delta=0.07, msg=f"Solidity {solidity:.3f} must be approximately 0.60.")
        is_dissected = is_botanical_dissection(cnt)
        self.assertFalse(
            is_dissected,
            "is_botanical_dissection(contour) must return False for unilateral foreign occlusion."
        )

    def test_fixture_5_petiole_stalk_trimming_and_restoration(self) -> None:
        """Fixture 5: Explicit LM3 petiole_mask subtraction and hole_mask subtraction in _prepare_lamina_mask."""
        f5_mask = create_fixture_5_leaf_with_petiole_stalk(self.f1_mask, stalk_length=45, stalk_width=10)
        orig_area = np.count_nonzero(self.f1_mask)
        stalk_area = np.count_nonzero(f5_mask)
        self.assertGreater(stalk_area, orig_area, "Fixture 5 mask must include attached stalk pixels.")

        petiole_mask = np.zeros_like(f5_mask)
        cv2.rectangle(petiole_mask, (195, 315), (205, 360), 255, -1)

        cleaned = _prepare_lamina_mask(f5_mask, petiole_mask=petiole_mask, hole_mask=None)
        self.assertEqual(int(np.count_nonzero(cleaned[325:, :])), 0)
        cleaned_area = np.count_nonzero(cleaned)
        area_diff = abs(cleaned_area - orig_area)
        self.assertLessEqual(area_diff, 25, f"Cleaned mask must match intact blade area within 25 px (diff={area_diff}).")

        ucs_clean, sol_clean, _, _, _ = compute_geometric_metrics(cleaned)
        self.assertGreaterEqual(sol_clean, 0.85)
        self.assertGreaterEqual(ucs_clean, 0.90)

    def test_fixture_6_cordate_blade_protection(self) -> None:
        """Fixture 6: Broad cordate blade base is 100% preserved and unclipped in _prepare_lamina_mask and resample_normalized_contour."""
        f6_mask = create_fixture_6_cordate_leaf()
        orig_area = np.count_nonzero(f6_mask)

        cleaned = _prepare_lamina_mask(f6_mask)
        self.assertEqual(np.count_nonzero(cleaned), orig_area, "Cordate leaf area must be 100% preserved.")

        pts = resample_normalized_contour(cleaned, num_points=120)
        self.assertEqual(pts.shape, (120, 2))

    def test_clean_unfolded_leaf_left_unclipped(self) -> None:
        """Verify that a cleanly segmented blade without petiole_mask remains 100% unclipped."""
        f1_mask = self.f1_mask
        orig_area = np.count_nonzero(f1_mask)

        cleaned = _prepare_lamina_mask(f1_mask)
        self.assertEqual(np.count_nonzero(cleaned), orig_area, "Intact blade area must be 100% preserved.")

        pts = resample_normalized_contour(cleaned, num_points=120)
        self.assertEqual(pts.shape, (120, 2))

    def test_fallback_morphological_linear_protrusion(self) -> None:
        """Verify that explicit LM3 petiole_mask removes narrow linear petiole protrusions cleanly."""
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 130), (50, 90), 0, 0, 360, 255, -1)
        cv2.rectangle(mask, (147, 220), (153, 245), 255, -1)

        petiole_mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.rectangle(petiole_mask, (147, 220), (153, 245), 255, -1)

        cleaned = _prepare_lamina_mask(mask, petiole_mask=petiole_mask)
        self.assertLess(np.count_nonzero(cleaned), np.count_nonzero(mask))
        self.assertEqual(int(np.count_nonzero(cleaned[222:, :])), 0)

    def test_estimate_rosette_crown_point_midrib_intersection_and_root_anchor(self) -> None:
        """Verify estimate_rosette_crown_point() intersects inward basal midrib rays and anchors on Class 9 roots."""
        basal_leaves = [
            {"leaf_id": 1, "plant_individual_id": 0, "p_apex": (120.0, 190.0), "p_base": (160.0, 270.0)},
            {"leaf_id": 2, "plant_individual_id": 0, "p_apex": (200.0, 170.0), "p_base": (200.0, 260.0)},
            {"leaf_id": 3, "plant_individual_id": 0, "p_apex": (280.0, 190.0), "p_base": (240.0, 270.0)},
        ]
        crown_pt = estimate_rosette_crown_point(basal_leaves)
        self.assertAlmostEqual(crown_pt[0], 200.0, delta=1.5)
        self.assertAlmostEqual(crown_pt[1], 350.0, delta=1.5)

        single_leaf = [{"leaf_id": 1, "plant_individual_id": 0, "p_apex": (200.0, 180.0), "p_base": (200.0, 280.0)}]
        crown_anchored = estimate_rosette_crown_point(
            single_leaf,
            root_boxes=[((340, 190, 380, 210), 0.92, 9)],
        )
        self.assertAlmostEqual(crown_anchored[0], 200.0, delta=2.0)
        self.assertGreater(crown_anchored[1], 310.0, "Stage 1 Class 9 root box must anchor p_crown below p_base.")

        multi_leaves = basal_leaves + [
            {"leaf_id": 4, "plant_individual_id": 1, "p_apex": (520.0, 190.0), "p_base": (560.0, 270.0)},
            {"leaf_id": 5, "plant_individual_id": 1, "p_apex": (680.0, 190.0), "p_base": (640.0, 270.0)},
        ]
        by_pid = estimate_rosette_crown_point(multi_leaves, return_by_plant=True)
        self.assertIn(0, by_pid)
        self.assertIn(1, by_pid)
        self.assertAlmostEqual(by_pid[0][0], 200.0, delta=1.5)
        self.assertAlmostEqual(by_pid[1][0], 600.0, delta=1.5)

    def test_approximate_petiole_and_scape_metrics_straight_and_bent_stem(self) -> None:
        """Verify per-leaf petiole_length_px = max(geodesic, ||p_base - p_crown||_2), 75th-%ile aggregation, and bent scape height."""
        leaves = [
            {
                "leaf_id": 1,
                "plant_individual_id": 0,
                "p_base": (200.0, 740.0),
                "p_apex": (200.0, 620.0),
                "geodesic_petiole_skeleton_px": 35.0,
            },
            {
                "leaf_id": 2,
                "plant_individual_id": 0,
                "p_base": (200.0, 760.0),
                "p_apex": (200.0, 660.0),
                "geodesic_petiole_skeleton_px": 90.0,
            },
        ]
        res_straight = approximate_petiole_and_scape_metrics(
            basal_leaves=leaves,
            p_crown=(200.0, 800.0),
            px_per_cm=100.0,
            p_corymb=(200.0, 200.0),
        )
        self.assertEqual(len(res_straight["leaf_metrics"]), 2)
        m1, m2 = res_straight["leaf_metrics"]
        self.assertAlmostEqual(m1["petiole_length_px"], 60.0, places=3)
        self.assertAlmostEqual(m1["blade_length_px"], 120.0, places=3)
        self.assertAlmostEqual(m1["petiole_to_blade_ratio"], 0.50, places=3)
        self.assertAlmostEqual(m1["petiole_length_mm"], 6.0, places=3)
        self.assertAlmostEqual(m2["petiole_length_px"], 90.0, places=3)
        self.assertAlmostEqual(m2["petiole_to_blade_ratio"], 0.90, places=3)
        self.assertAlmostEqual(m2["petiole_length_mm"], 9.0, places=3)

        self.assertAlmostEqual(res_straight["petiole_length_mm"], 8.25, places=2)
        self.assertAlmostEqual(res_straight["petiole_to_blade_ratio"], 0.80, places=3)
        self.assertFalse(res_straight["is_bent_stem"])
        self.assertAlmostEqual(res_straight["plant_height_cm"], 6.0, places=3)
        self.assertGreater(res_straight["scape_to_basal_leaf_ratio"], 2.5)

        res_bent = approximate_petiole_and_scape_metrics(
            basal_leaves=leaves,
            p_crown=(200.0, 800.0),
            px_per_cm=100.0,
            p_corymb=(500.0, 500.0),
            upper_cauline_detections=[(200.0, 100.0)],
        )
        self.assertTrue(res_bent["is_bent_stem"])
        self.assertAlmostEqual(res_bent["scape_length_px"], 1200.0, places=2)
        self.assertAlmostEqual(res_bent["plant_height_cm"], 12.0, places=2)
        self.assertGreater(res_bent["plant_height_cm"], res_straight["plant_height_cm"])

    def test_resample_normalized_contour(self) -> None:
        """Verify resampled contour normalization coordinates."""
        mask = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(mask, (100, 100), 40, 255, -1)

        coords = resample_normalized_contour(mask, num_points=120)
        self.assertIsNotNone(coords)
        self.assertEqual(coords.shape, (120, 2))
        self.assertGreaterEqual(np.min(coords), 0.0)
        self.assertLessEqual(np.max(coords), 1.0)

    def test_compute_classical_shape_descriptors_solidity_circularity_extent(self) -> None:
        """Verify compute_classical_shape_descriptors returns accurate Solidity, Circularity, and Extent."""
        mask = np.zeros((300, 300), dtype=np.uint8)
        cv2.ellipse(mask, (150, 150), (50, 90), 0, 0, 360, 255, -1)

        desc_tight = compute_classical_shape_descriptors(mask, bbox=(60, 100, 240, 200))
        self.assertGreaterEqual(desc_tight["solidity"], 0.95)
        self.assertGreaterEqual(desc_tight["circularity"], 0.65)
        self.assertGreaterEqual(desc_tight["extent"], 0.70)

        desc_partial = compute_classical_shape_descriptors(mask, bbox=(0, 0, 300, 300))
        self.assertLess(desc_partial["extent"], 0.35)

    def test_rotation_invariant_extent_preserves_oversized_bbox_penalty(self) -> None:
        """Verify 45-degree rotated leaves achieve high oriented Extent while oversized external bboxes still fail."""
        mask_rot = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(mask_rot, (200, 200), (36, 130), 45, 0, 360, 255, -1)
        bx, by, bw, bh = cv2.boundingRect(mask_rot)
        tight_bbox = (by, bx, by + bh, bx + bw)

        desc_rot = compute_classical_shape_descriptors(mask_rot, bbox=tight_bbox)
        self.assertGreaterEqual(desc_rot["extent"], 0.72)

        desc_oversized = compute_classical_shape_descriptors(mask_rot, bbox=(0, 0, 400, 400))
        self.assertLess(desc_oversized["extent"], 0.35)


if __name__ == "__main__":
    unittest.main()
