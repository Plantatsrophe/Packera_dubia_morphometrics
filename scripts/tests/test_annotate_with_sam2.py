#!/usr/bin/env python3
"""
Unit test for PrecisionSAM2Annotator in scripts/annotation_and_training/annotate_with_sam2.py
and helper utilities in scripts/annotation_and_training/sam2_annotator_utils.py.
Verifies polygon extraction, label mapping, binary mask generation, COCO export, and filename parsing.
"""

import json
import tempfile
import unittest
from pathlib import Path
import cv2
import numpy as np

from scripts.annotation_and_training.annotate_with_sam2 import PrecisionSAM2Annotator, CLASS_NAMES
from scripts.annotation_and_training.sam2_annotator_utils import (
    convert_masks_to_coco_dataset,
    export_coco_annotations,
    parse_mask_filename,
    polygon_interior_point,
    polygon_to_bounding_box,
    rasterize_lasso_polygon,
)


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
            coco_out = tmp_p / "coco" / "annotations_packera_train.json"
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
            annotator.coco_output = coco_out
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

            # Verify COCO output was created and synchronized
            self.assertTrue(coco_out.exists())
            with open(coco_out, "r") as f:
                coco_data = json.load(f)
            self.assertIn("images", coco_data)
            self.assertIn("annotations", coco_data)
            self.assertIn("categories", coco_data)

    def test_polygon_bounding_box_and_interior_selection(self):
        poly = [(10, 20), (50, 15), (75, 60), (30, 80), (10, 40)]
        bbox = polygon_to_bounding_box(poly)
        self.assertIsNotNone(bbox)
        self.assertEqual(bbox, (10, 15, 75, 80))

        interior = polygon_interior_point(poly, img_h=100, img_w=100)
        self.assertIsNotNone(interior)
        ix, iy = interior
        self.assertTrue(10 <= ix <= 75)
        self.assertTrue(15 <= iy <= 80)

        mask = rasterize_lasso_polygon(poly, 100, 100)
        self.assertEqual(mask[int(iy), int(ix)], 255)

        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.orig_h = 100
        annotator.orig_w = 100
        annotator.polygon_points = [(10, 20), (50, 20), (50, 60), (10, 60)]
        annotator.candidate_mask = None
        annotator.point_coords = []
        annotator.point_labels = []
        annotator.box_prompt = None
        annotator.predictor = None

        annotator.finalize_polygon_selection()

        self.assertEqual(annotator.box_prompt, [10.0, 20.0, 50.0, 60.0])
        self.assertIsNotNone(annotator.candidate_mask)
        self.assertGreater(np.count_nonzero(annotator.candidate_mask), 0)
        self.assertEqual(annotator.polygon_points, [])

    def test_parse_mask_filename(self):
        cat1, lbl1, inst1 = parse_mask_filename("000331814_inst00_basal_leaf_partial.png")
        self.assertEqual(cat1, "000331814")
        self.assertEqual(lbl1, "basal_leaf_partial")
        self.assertEqual(inst1, 0)

        cat2, lbl2, inst2 = parse_mask_filename("NCU00001234_basal_leaf_whole_instance_1.png")
        self.assertEqual(cat2, "NCU00001234")
        self.assertEqual(lbl2, "basal_leaf_whole")
        self.assertEqual(inst2, 1)

    def test_export_coco_annotations(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_p = Path(tmp_dir)
            masks_dir = tmp_p / "masks"
            masks_dir.mkdir(parents=True, exist_ok=True)
            coco_path = tmp_p / "annotations_packera_train.json"

            # Create dummy mask for ideal_leaf (basal_leaf_whole)
            mask_img = np.zeros((200, 200), dtype=np.uint8)
            cv2.rectangle(mask_img, (30, 30), (100, 100), 255, -1)
            cv2.imwrite(str(masks_dir / "VOUCHER001_inst00_basal_leaf_whole.png"), mask_img)

            # Create dummy mask for partial_leaf (basal_leaf_partial)
            mask_img2 = np.zeros((200, 200), dtype=np.uint8)
            cv2.circle(mask_img2, (150, 150), 30, 255, -1)
            cv2.imwrite(str(masks_dir / "VOUCHER001_inst01_basal_leaf_partial.png"), mask_img2)

            coco_doc = export_coco_annotations(
                masks_dir=masks_dir,
                output_coco_path=coco_path,
                min_area_px=10.0,
            )

            self.assertTrue(coco_path.exists())
            self.assertEqual(len(coco_doc["images"]), 1)
            self.assertEqual(len(coco_doc["annotations"]), 2)
            # Check categories: ideal_leaf (1) and partial_leaf (2)
            cat_ids = {ann["category_id"] for ann in coco_doc["annotations"]}
            self.assertEqual(cat_ids, {1, 2})


if __name__ == "__main__":
    unittest.main()
