"""
Unit Tests for Botanical Organ Topology Classifier
==================================================
Tests medial axis skeletonization, 3x3 neighbor convolution, branch/endpoint identification,
tortuosity computation, and spatial prior decision logic for differentiating
root_rhizome, cauline_stem, and leaf_petiole.
"""

from __future__ import annotations

import unittest
import numpy as np
import cv2

from scripts.core.botanical_topology_classifier import (
    ROOT_RHIZOME,
    CAULINE_STEM,
    LEAF_PETIOLE,
    NEIGHBOR_KERNEL,
    analyze_skeleton_topology,
    calculate_mask_centroid,
    calculate_skeleton_path_and_chord_length,
    classify_elongated_botanical_organ,
    compute_topological_summary,
    extract_medial_axis_skeleton,
    generate_synthetic_stem_mask,
    generate_synthetic_root_mask,
    generate_synthetic_petiole_with_blade_mask,
)


class TestBotanicalTopologyClassifier(unittest.TestCase):
    """Test suite for botanical organ topology classification and morphometrics."""

    def test_neighbor_convolution_kernel_values(self):
        """
        Verify exact convolution response for isolated points (10), endpoints (11),
        straight line points (12), 3-way junctions (13), and 4-way cross junctions (14).
        """
        kernel_f32 = NEIGHBOR_KERNEL.astype(np.float32)
        # Test 1: Isolated point
        grid_iso = np.zeros((5, 5), dtype=np.float32)
        grid_iso[2, 2] = 1.0
        conv_iso = np.round(cv2.filter2D(grid_iso, cv2.CV_32F, kernel_f32, borderType=cv2.BORDER_CONSTANT)).astype(int)
        self.assertEqual(conv_iso[2, 2], 10)

        # Test 2: Endpoint (1 neighbor)
        grid_end = np.zeros((5, 5), dtype=np.float32)
        grid_end[2, 2] = 1.0
        grid_end[2, 3] = 1.0
        conv_end = np.round(cv2.filter2D(grid_end, cv2.CV_32F, kernel_f32, borderType=cv2.BORDER_CONSTANT)).astype(int)
        self.assertEqual(conv_end[2, 2], 11)
        self.assertEqual(conv_end[2, 3], 11)

        # Test 3: Segment point (2 neighbors)
        grid_seg = np.zeros((5, 5), dtype=np.float32)
        grid_seg[2, 1] = 1.0
        grid_seg[2, 2] = 1.0
        grid_seg[2, 3] = 1.0
        conv_seg = np.round(cv2.filter2D(grid_seg, cv2.CV_32F, kernel_f32, borderType=cv2.BORDER_CONSTANT)).astype(int)
        self.assertEqual(conv_seg[2, 2], 12)
        self.assertEqual(conv_seg[2, 1], 11)
        self.assertEqual(conv_seg[2, 3], 11)

        # Test 4: 3-way T-junction (3 neighbors)
        grid_junc = np.zeros((5, 5), dtype=np.float32)
        grid_junc[2, 1] = 1.0
        grid_junc[2, 2] = 1.0
        grid_junc[2, 3] = 1.0
        grid_junc[3, 2] = 1.0
        conv_junc = np.round(cv2.filter2D(grid_junc, cv2.CV_32F, kernel_f32, borderType=cv2.BORDER_CONSTANT)).astype(int)
        self.assertEqual(conv_junc[2, 2], 13)

        # Test 5: 4-way cross junction (4 neighbors)
        grid_cross = np.zeros((5, 5), dtype=np.float32)
        grid_cross[2, 2] = 1.0
        grid_cross[1, 2] = 1.0
        grid_cross[3, 2] = 1.0
        grid_cross[2, 1] = 1.0
        grid_cross[2, 3] = 1.0
        conv_cross = np.round(cv2.filter2D(grid_cross, cv2.CV_32F, kernel_f32, borderType=cv2.BORDER_CONSTANT)).astype(int)
        self.assertEqual(conv_cross[2, 2], 14)

    def test_classify_has_connected_blade(self):
        """Rule 1: If has_connected_blade is True, must immediately return leaf_petiole."""
        # Even with high branches and low on sheet (like roots), blade attachment forces leaf_petiole
        mask = generate_synthetic_root_mask(shape=(1000, 1000), origin=(500, 800))
        result = classify_elongated_botanical_organ(
            mask=mask,
            full_sheet_height=1000,
            y_centroid=800,
            has_connected_blade=True
        )
        self.assertEqual(result, LEAF_PETIOLE)

    def test_classify_linear_stem(self):
        """
        Rule 5.2: Continuous vertical axis (norm_y < 0.70, branch points <= 2)
        must classify as cauline_stem.
        """
        # Stem spanning upper/middle sheet from y=100 to y=500 on a 1000px sheet (norm_y = 0.30)
        stem_mask = generate_synthetic_stem_mask(
            shape=(1000, 1000),
            start_pt=(500, 100),
            end_pt=(500, 500),
            thickness=8,
            waviness_amplitude=0.0
        )
        result = classify_elongated_botanical_organ(
            mask=stem_mask,
            full_sheet_height=1000,
            has_connected_blade=False
        )
        self.assertEqual(result, CAULINE_STEM)

    def test_classify_curved_stem(self):
        """A wavy cauline stem with minimal branching in upper sheet should classify as cauline_stem."""
        wavy_stem_mask = generate_synthetic_stem_mask(
            shape=(1000, 1000),
            start_pt=(500, 150),
            end_pt=(500, 600),
            thickness=6,
            waviness_amplitude=25.0,
            waviness_frequency=0.03
        )
        result = classify_elongated_botanical_organ(
            mask=wavy_stem_mask,
            full_sheet_height=1000,
            has_connected_blade=False
        )
        self.assertEqual(result, CAULINE_STEM)

    def test_classify_branched_root_rhizome(self):
        """
        Rule 5.1: Subterranean highly branched organ (norm_y > 0.65, branch points >= 3)
        must classify as root_rhizome.
        """
        root_mask = generate_synthetic_root_mask(
            shape=(1000, 1000),
            origin=(500, 720),
            num_branches=4,
            max_depth=3,
            branch_length=50,
            thickness=5
        )
        # Verify centroid is > 0.65
        _, y_c = calculate_mask_centroid(root_mask)
        self.assertGreater(y_c / 1000.0, 0.65)

        # Verify skeleton topology has at least 3 branch points
        skel = extract_medial_axis_skeleton(root_mask)
        _, _, _, num_branches, _, _ = analyze_skeleton_topology(skel)
        self.assertGreaterEqual(num_branches, 3)

        result = classify_elongated_botanical_organ(
            mask=root_mask,
            full_sheet_height=1000,
            has_connected_blade=False
        )
        self.assertEqual(result, ROOT_RHIZOME)

    def test_fallback_classification_logic(self):
        """
        Fallback logic:
        - If norm_y > 0.50 (and doesn't meet root branching >=3 or stem norm_y < 0.70), fallback to leaf_petiole.
        - If norm_y <= 0.50, fallback to cauline_stem.
        """
        # Create an unbranched horizontal petiole stalk near base (e.g. norm_y = 0.75, branches = 0)
        # norm_y = 0.75 is >= 0.70 (fails stem rule), branches = 0 < 3 (fails root rule)
        # Fallback norm_y > 0.50 -> leaf_petiole
        petiole_mask = np.zeros((1000, 1000), dtype=np.uint8)
        cv2.line(petiole_mask, (400, 750), (600, 750), 255, 4)

        result = classify_elongated_botanical_organ(
            mask=petiole_mask,
            full_sheet_height=1000,
            has_connected_blade=False
        )
        self.assertEqual(result, LEAF_PETIOLE)

        # A synthetic short stalk in upper sheet with norm_y = 0.40
        upper_mask = np.zeros((1000, 1000), dtype=np.uint8)
        cv2.line(upper_mask, (450, 400), (550, 400), 255, 4)
        result_upper = classify_elongated_botanical_organ(
            mask=upper_mask,
            full_sheet_height=1000,
            has_connected_blade=False
        )
        self.assertEqual(result_upper, CAULINE_STEM)

    def test_path_length_and_tortuosity_straight_vs_wavy(self):
        """
        Tortuosity of a straight line should be approximately 1.0,
        while a wavy / sinusoidal line must have tortuosity significantly > 1.0.
        """
        straight_mask = generate_synthetic_stem_mask(
            shape=(1000, 1000),
            start_pt=(500, 100),
            end_pt=(500, 500),
            thickness=4,
            waviness_amplitude=0.0
        )
        skel_straight = extract_medial_axis_skeleton(straight_mask)
        _, _, end_straight, _, _, _ = analyze_skeleton_topology(skel_straight)
        path_s, chord_s, tau_straight = calculate_skeleton_path_and_chord_length(skel_straight, end_straight)

        self.assertAlmostEqual(tau_straight, 1.0, delta=0.05)
        self.assertAlmostEqual(chord_s, 400.0, delta=5.0)

        wavy_mask = generate_synthetic_stem_mask(
            shape=(1000, 1000),
            start_pt=(500, 100),
            end_pt=(500, 500),
            thickness=4,
            waviness_amplitude=50.0,
            waviness_frequency=0.04
        )
        skel_wavy = extract_medial_axis_skeleton(wavy_mask)
        _, _, end_wavy, _, _, _ = analyze_skeleton_topology(skel_wavy)
        path_w, chord_w, tau_wavy = calculate_skeleton_path_and_chord_length(skel_wavy, end_wavy)

        self.assertGreater(tau_wavy, 1.15)
        self.assertGreater(path_w, chord_w)

    def test_compute_topological_summary_export(self):
        """
        Verify that compute_topological_summary returns a comprehensive dictionary/dataclass
        with expected fields and valid scaling for morphometric tables.
        """
        stem_mask = generate_synthetic_stem_mask(
            shape=(1000, 1000),
            start_pt=(500, 100),
            end_pt=(500, 400),
            thickness=6
        )
        pixel_size_mm = 0.05  # 50 microns per pixel
        metrics = compute_topological_summary(
            mask=stem_mask,
            full_sheet_height=1000,
            pixel_size_mm=pixel_size_mm,
            has_connected_blade=False
        )

        self.assertEqual(metrics.predicted_class, CAULINE_STEM)
        self.assertGreater(metrics.skeleton_pixel_count, 280)
        self.assertGreater(metrics.mask_area_px, 1000)
        self.assertIsNotNone(metrics.path_length_mm)
        self.assertAlmostEqual(metrics.path_length_mm, metrics.path_length_px * pixel_size_mm, places=4)
        self.assertAlmostEqual(metrics.chord_length_mm, metrics.chord_length_px * pixel_size_mm, places=4)
        self.assertAlmostEqual(metrics.mask_area_mm2, metrics.mask_area_px * (pixel_size_mm ** 2), places=4)

        data_dict = metrics.to_dict()
        self.assertIn("predicted_class", data_dict)
        self.assertIn("tortuosity", data_dict)
        self.assertIn("num_branch_points", data_dict)
        self.assertIn("num_end_points", data_dict)
        self.assertIn("norm_y", data_dict)

    def test_edge_cases(self):
        """Test empty masks, single pixel masks, and boundary coordinates."""
        empty_mask = np.zeros((100, 100), dtype=np.uint8)
        skel = extract_medial_axis_skeleton(empty_mask)
        self.assertEqual(np.count_nonzero(skel), 0)

        metrics_empty = compute_topological_summary(empty_mask, full_sheet_height=1000)
        self.assertEqual(metrics_empty.skeleton_pixel_count, 0)
        self.assertEqual(metrics_empty.num_branch_points, 0)
        self.assertEqual(metrics_empty.tortuosity, 1.0)

        # Single pixel mask
        single_px = np.zeros((100, 100), dtype=np.uint8)
        single_px[50, 50] = 1
        metrics_single = compute_topological_summary(single_px, full_sheet_height=100)
        self.assertEqual(metrics_single.skeleton_pixel_count, 1)
        self.assertEqual(metrics_single.tortuosity, 1.0)


if __name__ == "__main__":
    unittest.main()
