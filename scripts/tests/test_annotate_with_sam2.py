#!/usr/bin/env python3
"""
Unit test for PrecisionSAM2Annotator in scripts/data_prep/annotate_with_sam2.py
Verifies polygon extraction, label mapping, and binary mask generation.
"""

import tempfile
import unittest
from pathlib import Path
import numpy as np
import cv2

from scripts.data_prep.annotate_with_sam2 import PrecisionSAM2Annotator, CLASS_NAMES

class TestPrecisionSAM2Annotator(unittest.TestCase):
    def test_class_names(self):
        self.assertIn("basal_leaf_whole", CLASS_NAMES)
        self.assertIn("basal_leaf_partial", CLASS_NAMES)
        self.assertIn("cauline_leaf", CLASS_NAMES)
        self.assertEqual(CLASS_NAMES[0], "basal_leaf_whole")
        self.assertEqual(CLASS_NAMES[1], "basal_leaf_partial")
        self.assertEqual(CLASS_NAMES[2], "cauline_leaf")

    def test_save_and_mask_export(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_p = Path(tmp_dir)
            img_dir = tmp_p / "vouchers"
            out_dir = tmp_p / "annotations"
            img_dir.mkdir(parents=True, exist_ok=True)
            
            # Create a dummy image
            dummy_img_path = img_dir / "VOUCHER_TEST001.jpg"
            dummy_img = np.zeros((100, 100, 3), dtype=np.uint8)
            cv2.imwrite(str(dummy_img_path), dummy_img)

            # Initialize mock instance without running model inference
            annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
            annotator.project_root = tmp_p
            annotator.images_dir = img_dir
            annotator.output_dir = out_dir
            annotator.masks_dir = out_dir / "masks"
            annotator.output_dir.mkdir(parents=True, exist_ok=True)
            annotator.masks_dir.mkdir(parents=True, exist_ok=True)
            annotator.image_files = [dummy_img_path]
            annotator.current_idx = 0
            annotator.orig_h = 100
            annotator.orig_w = 100

            # Mock saved instances for requested labels
            mask1 = np.zeros((100, 100), dtype=bool)
            mask1[10:50, 10:50] = True
            
            annotator.saved_instances = [
                {
                    "class_id": 0,
                    "label": "basal_leaf_whole",
                    "polygon": [0.1, 0.1, 0.5, 0.1, 0.5, 0.5, 0.1, 0.5],
                    "binary_mask": mask1
                },
                {
                    "class_id": 1,
                    "label": "basal_leaf_partial",
                    "polygon": [0.2, 0.2, 0.4, 0.2, 0.4, 0.4, 0.2, 0.4],
                    "binary_mask": mask1
                },
                {
                    "class_id": 2,
                    "label": "cauline_leaf",
                    "polygon": [0.6, 0.6, 0.8, 0.6, 0.8, 0.8, 0.6, 0.8],
                    "binary_mask": mask1
                }
            ]

            annotator.save_current_sheet()

            # Verify YOLO polygon text file exists
            txt_file = out_dir / "VOUCHER_TEST001.txt"
            self.assertTrue(txt_file.exists())
            with open(txt_file, "r") as f:
                lines = f.readlines()
                self.assertEqual(len(lines), 3)
                self.assertTrue(lines[0].startswith("0 "))
                self.assertTrue(lines[1].startswith("1 "))
                self.assertTrue(lines[2].startswith("2 "))

            # Verify binary pixel mask files exist and are tagged explicitly
            m0 = out_dir / "masks" / "VOUCHER_TEST001_inst00_basal_leaf_whole.png"
            m1 = out_dir / "masks" / "VOUCHER_TEST001_inst01_basal_leaf_partial.png"
            m2 = out_dir / "masks" / "VOUCHER_TEST001_inst02_cauline_leaf.png"

            self.assertTrue(m0.exists())
            self.assertTrue(m1.exists())
            self.assertTrue(m2.exists())

            # Read back binary mask image to confirm binary uint8 format (0/255)
            img_m0 = cv2.imread(str(m0), cv2.IMREAD_GRAYSCALE)
            self.assertEqual(img_m0.shape, (100, 100))
            self.assertEqual(set(np.unique(img_m0)), {0, 255})

    def test_polygon_bounding_box_and_interior_selection(self):
        from scripts.data_prep.sam2_geometry import polygon_to_bounding_box, polygon_interior_point, rasterize_lasso_polygon
        
        # Test 1: Bounding box from vertices
        poly = [(10, 20), (50, 15), (75, 60), (30, 80), (10, 40)]
        bbox = polygon_to_bounding_box(poly)
        self.assertIsNotNone(bbox)
        self.assertEqual(bbox, (10, 15, 75, 80))

        # Test 2: Interior point inside mask
        interior = polygon_interior_point(poly, img_h=100, img_w=100)
        self.assertIsNotNone(interior)
        ix, iy = interior
        # Verify interior point falls within bounding box
        self.assertTrue(10 <= ix <= 75)
        self.assertTrue(15 <= iy <= 80)
        # Verify interior point falls inside rasterized polygon mask
        mask = rasterize_lasso_polygon(poly, 100, 100)
        self.assertEqual(mask[int(iy), int(ix)], 255)

        # Test 3: Annotator finalize_polygon_selection integration
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.orig_h = 100
        annotator.orig_w = 100
        annotator.polygon_points = [(10, 20), (50, 20), (50, 60), (10, 60)]
        annotator.candidate_mask = None
        annotator.point_coords = []
        annotator.point_labels = []
        annotator.box_prompt = None
        annotator.predictor = None  # Mock no predictor to test fallback

        annotator.finalize_polygon_selection()

        # Bounding box prompt should be set
        self.assertEqual(annotator.box_prompt, [10.0, 20.0, 50.0, 60.0])
        # Candidate mask should be rasterized polygon
        self.assertIsNotNone(annotator.candidate_mask)
        self.assertGreater(np.count_nonzero(annotator.candidate_mask), 0)
        self.assertEqual(annotator.polygon_points, [])


if __name__ == "__main__":
    unittest.main()

