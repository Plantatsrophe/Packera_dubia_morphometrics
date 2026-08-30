#!/usr/bin/env python3
"""
===============================================================================
Script: test_postprocess_lm2_routing.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit and integration tests for Step 02:
    LeafMachine2 post-processing, DBSCAN spatial clustering, geometric
    gatekeeping (solidity, convexity, pose), and 4-tier morphometric routing.
===============================================================================
"""

import unittest
from pathlib import Path
import numpy as np
import cv2

PROJECT_ROOT = Path(__file__).resolve().parents[2]

from scripts.vision.geometric_gatekeeper import compute_geometric_metrics_and_pose
from scripts.vision.spatial_clustering import cluster_voucher_plants_dbscan
from scripts.vision.lm2_data_loader import LeafCandidate


class TestLM2PostProcessingAndRouting(unittest.TestCase):
    """Test suite for Step 02 vision gatekeeping, clustering, and routing."""

    def test_geometric_gatekeeper_pristine_mask(self):
        """Verify that a solid elliptic contour achieves high solidity (Tier 1)."""
        # Create a synthetic solid ellipse mask (height 100, width 40)
        mask = np.zeros((200, 200), dtype=np.uint8)
        cv2.ellipse(mask, (100, 100), (20, 50), 0, 0, 360, 255, -1)

        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics_and_pose(mask)
        self.assertGreaterEqual(solidity, 0.72, "Solid ellipse should pass solidity threshold.")
        self.assertGreater(ucs, 0.0)
        self.assertIsInstance(angle_deg, float)

    def test_spatial_clustering_dbscan(self):
        """Verify that leaf candidates separated spatially are clustered into distinct plants."""
        cand1 = LeafCandidate(catalog_number="TEST1", leaf_id=1, bbox=(10, 10, 30, 30))
        cand2 = LeafCandidate(catalog_number="TEST1", leaf_id=2, bbox=(15, 15, 35, 35))
        cand3 = LeafCandidate(catalog_number="TEST1", leaf_id=3, bbox=(800, 800, 820, 820))
        cand4 = LeafCandidate(catalog_number="TEST1", leaf_id=4, bbox=(810, 810, 830, 830))

        candidates = [cand1, cand2, cand3, cand4]
        cluster_voucher_plants_dbscan(candidates, sheet_width=1000, sheet_height=1000)

        self.assertEqual(cand1.plant_individual_id, cand2.plant_individual_id, "Candidates in close proximity must share plant_individual_id.")
        self.assertEqual(cand3.plant_individual_id, cand4.plant_individual_id, "Candidates in close proximity must share plant_individual_id.")
        self.assertNotEqual(cand1.plant_individual_id, cand3.plant_individual_id, "Distant candidates must have different plant_individual_ids.")

    def test_entrypoint_exists(self):
        """Verify that the standardized step 02 script exists."""
        script_path = PROJECT_ROOT / "scripts" / "vision" / "02_postprocess_lm2_routing.py"
        self.assertTrue(script_path.exists(), "02_postprocess_lm2_routing.py must exist.")


if __name__ == "__main__":
    unittest.main()
