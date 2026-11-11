#!/usr/bin/env python3
"""
scripts/tests/test_geometry_edge_cases.py
=========================================
Unit tests verifying LeafMachine3 shape metric delegation, contour resampling,
rosette crown ray intersection, and bent-stem scape measurements.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.leaf_geometry_utils import (
    approximate_petiole_and_scape_metrics,
    compute_lm3_shape_metrics,
    estimate_rosette_crown_point,
    resample_normalized_contour,
)


class TestGeometryEdgeCases(unittest.TestCase):
    """Unit tests verifying LM3 shape metrics, contour resampling, and organ-crown geometry."""

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
        # Direct petiole mask (35.0 px) takes strict priority over crown_fallback_px (60.0 px)
        self.assertAlmostEqual(m1["petiole_length_px"], 35.0, places=3)
        self.assertAlmostEqual(m1["blade_length_px"], 120.0, places=3)
        self.assertAlmostEqual(m1["petiole_to_blade_ratio"], 35.0 / 120.0, places=3)
        self.assertAlmostEqual(m1["petiole_length_mm"], 3.5, places=3)
        self.assertAlmostEqual(m2["petiole_length_px"], 90.0, places=3)
        self.assertAlmostEqual(m2["petiole_to_blade_ratio"], 0.90, places=3)
        self.assertAlmostEqual(m2["petiole_length_mm"], 9.0, places=3)

        self.assertAlmostEqual(res_straight["petiole_length_mm"], 7.625, places=2)
        self.assertAlmostEqual(res_straight["petiole_to_blade_ratio"], 0.75 * 0.90 + 0.25 * (35.0 / 120.0), places=3)
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

        leaves_with_cauline = leaves + [
            {
                "leaf_id": 3,
                "plant_individual_id": 0,
                "leaf_zone": "cauline",
                "p_base": (200.0, 320.0),
                "p_apex": (200.0, 240.0),
                "geodesic_petiole_skeleton_px": 0.0,
            }
        ]
        res_caul = approximate_petiole_and_scape_metrics(
            basal_leaves=leaves_with_cauline,
            p_crown=(200.0, 800.0),
            px_per_cm=100.0,
            p_corymb=(200.0, 200.0),
        )
        self.assertEqual(len(res_caul["leaf_metrics"]), 3)
        m_caul = [r for r in res_caul["leaf_metrics"] if r["leaf_id"] == 3][0]
        self.assertAlmostEqual(m_caul["petiole_length_mm"], 0.0, places=3)
        self.assertAlmostEqual(m_caul["blade_length_mm"], 8.0, places=3)
        self.assertAlmostEqual(res_caul["petiole_length_mm"], 7.625, places=2)

        # Verify crown_fallback_px activates when direct petiole skeleton is absent (None)
        res_fb = approximate_petiole_and_scape_metrics(
            basal_leaves=[{"leaf_id": 1, "plant_individual_id": 0, "p_base": (200.0, 740.0), "p_apex": (200.0, 620.0), "geodesic_petiole_skeleton_px": None}],
            p_crown=(200.0, 800.0),
            px_per_cm=100.0,
            p_corymb=(200.0, 200.0),
        )
        self.assertAlmostEqual(res_fb["leaf_metrics"][0]["petiole_length_px"], 60.0, places=3)

    def test_resample_normalized_contour(self) -> None:
        """Verify resampled contour normalization coordinates."""
        mask = np.zeros((200, 200), dtype=np.uint8)
        cv2.circle(mask, (100, 100), 40, 255, -1)

        coords = resample_normalized_contour(mask, num_points=120)
        self.assertIsNotNone(coords)
        self.assertEqual(coords.shape, (120, 2))
        self.assertGreaterEqual(np.min(coords), 0.0)
        self.assertLessEqual(np.max(coords), 1.0)

    def test_compute_lm3_shape_metrics_delegates_to_polygon_morphology(self) -> None:
        """Verify compute_lm3_shape_metrics delegates solidity, circularity, and oriented extent to LeafMachine3."""
        mask_rot = np.zeros((400, 400), dtype=np.uint8)
        cv2.ellipse(mask_rot, (200, 200), (36, 130), 45, 0, 360, 255, -1)
        bx, by, bw, bh = cv2.boundingRect(mask_rot)
        tight_bbox = (by, bx, by + bh, bx + bw)

        desc_rot = compute_lm3_shape_metrics(mask_rot, bbox=tight_bbox)
        self.assertGreaterEqual(desc_rot["solidity"], 0.95)
        self.assertGreaterEqual(desc_rot["circularity"], 0.50)
        self.assertGreaterEqual(desc_rot["extent"], 0.72)

        desc_oversized = compute_lm3_shape_metrics(mask_rot, bbox=(0, 0, 400, 400))
        self.assertLess(desc_oversized["extent"], 0.35)


if __name__ == "__main__":
    unittest.main()
