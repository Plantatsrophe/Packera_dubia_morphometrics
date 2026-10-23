#!/usr/bin/env python3
"""
Unit tests for SAM 2 native-DPI single-leaf crop COCO export in
`scripts/annotation_and_training/sam2_annotator_utils.py` and
`scripts/annotation_and_training/annotate_with_sam2.py`.

Verifies:
  1. `export_coco_from_yolo_and_masks()` exports native-DPI cropped single-leaf
     patches by default (`full_sheet=False`) for `basal_leaf_whole` (`ideal_leaf`,
     category_id=1) and `basal_leaf_partial` (`partial_leaf`, category_id=2).
  2. Tight polygon bounding boxes `(x_min, y_min, w, h)` are expanded isotropically
     by 12% padding on each side (`padding_ratio=0.12`) and clamped to voucher
     sheet boundaries.
  3. Cropped RGB patches are saved as `{catalogNumber}_inst{idx}.jpg` at native DPI
     and match the COCO `images` entry dimensions `(crop_w, crop_h)`.
  4. Translated segmentation polygon vertices `(x - x0, y - y0)` and local bounding
     boxes lie strictly within `[0, crop_w] x [0, crop_h]`.
  5. Optional `--full-sheet` (`full_sheet=True`) mode preserves full-sheet COCO export.
"""

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import cv2
import numpy as np

from scripts.annotation_and_training.annotate_with_sam2 import (
    PrecisionSAM2Annotator,
    parse_args as parse_annotator_args,
)
from scripts.annotation_and_training.sam2_annotator_utils import (
    DEFAULT_CROP_PADDING_RATIO,
    DEFAULT_LEAF_CROPS_DIR,
    convert_masks_to_coco_dataset,
    export_coco_annotations,
    export_coco_from_yolo_and_masks,
    mask_to_polygons,
)


class TestSAM2AnnotatorLeafCropExport(unittest.TestCase):
    """Tests for LeafMachine2 Stage 2 PointRend native-DPI cropped leaf patch COCO export."""

    def test_default_constants(self):
        self.assertAlmostEqual(DEFAULT_CROP_PADDING_RATIO, 0.12, places=4)
        self.assertEqual(DEFAULT_LEAF_CROPS_DIR, Path("data/annotations/lm2_leaf_crops"))

    def test_export_coco_from_yolo_and_masks_native_dpi_crops_and_local_coords(self):
        """
        Verifies 12% isotropic padding, native-DPI RGB patch saving, and local
        coordinate frame translation strictly within [0, crop_w] x [0, crop_h].
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_p = Path(tmp_dir)
            images_dir = tmp_p / "raw_vouchers"
            raw_ann_dir = tmp_p / "raw_annotations"
            masks_dir = raw_ann_dir / "masks"
            crops_dir = tmp_p / "annotations" / "lm2_leaf_crops"
            coco_path = tmp_p / "annotations" / "annotations_packera_train.json"

            images_dir.mkdir(parents=True, exist_ok=True)
            masks_dir.mkdir(parents=True, exist_ok=True)

            sheet_h, sheet_w = 800, 600
            voucher_id = "NCU00099887"
            sheet_bgr = np.zeros((sheet_h, sheet_w, 3), dtype=np.uint8)
            # Distinct RGB pattern in leaf regions to verify native-DPI pixel-exact cropping
            sheet_bgr[100:260, 120:240] = (34, 139, 34)
            sheet_bgr[400:550, 300:480] = (50, 205, 50)
            cv2.imwrite(str(images_dir / f"{voucher_id}.jpg"), sheet_bgr)

            # Instance 0: basal_leaf_whole (ideal_leaf, category_id=1)
            mask0 = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
            cv2.rectangle(mask0, (120, 100), (220, 250), 255, -1)
            cv2.imwrite(str(masks_dir / f"{voucher_id}_inst00_basal_leaf_whole.png"), mask0)

            # Instance 1: basal_leaf_partial (partial_leaf, category_id=2)
            mask1 = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
            cv2.rectangle(mask1, (300, 400), (450, 520), 255, -1)
            cv2.imwrite(str(masks_dir / f"{voucher_id}_inst01_basal_leaf_partial.png"), mask1)

            # Instance 2: cauline_leaf (omitted class, should not be exported)
            mask2 = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
            cv2.rectangle(mask2, (50, 600), (150, 700), 255, -1)
            cv2.imwrite(str(masks_dir / f"{voucher_id}_inst02_cauline_leaf.png"), mask2)

            coco_doc = export_coco_from_yolo_and_masks(
                masks_dir=masks_dir,
                output_coco_path=coco_path,
                images_dir=images_dir,
                crops_dir=crops_dir,
                min_area_px=20.0,
            )

            self.assertTrue(coco_path.exists())
            self.assertEqual(len(coco_doc["images"]), 2)
            self.assertEqual(len(coco_doc["annotations"]), 2)

            # Check expected 12% isotropic padding on Instance 0:
            # Tight polygon bbox from rectangle (120, 100) to (220, 250) -> w=100, h=150
            # pad_x = 12.0 -> x0 = 108, x1 = 232 -> crop_w = 124
            # pad_y = 18.0 -> y0 = 82,  y1 = 268 -> crop_h = 186
            img0 = coco_doc["images"][0]
            self.assertEqual(img0["file_name"], f"{voucher_id}_inst0.jpg")
            self.assertEqual(img0["width"], 124)
            self.assertEqual(img0["height"], 186)

            crop0_path = crops_dir / f"{voucher_id}_inst0.jpg"
            self.assertTrue(crop0_path.exists())
            crop0_disk = cv2.imread(str(crop0_path))
            self.assertIsNotNone(crop0_disk)
            self.assertEqual(crop0_disk.shape, (186, 124, 3))

            # Verify every annotation's polygon vertices and bbox lie strictly within [0, crop_w] x [0, crop_h]
            img_map = {im["id"]: im for im in coco_doc["images"]}
            for ann in coco_doc["annotations"]:
                im_info = img_map[ann["image_id"]]
                crop_w = im_info["width"]
                crop_h = im_info["height"]

                crop_file = crops_dir / im_info["file_name"]
                self.assertTrue(crop_file.exists())
                disk_im = cv2.imread(str(crop_file))
                self.assertEqual((disk_im.shape[1], disk_im.shape[0]), (crop_w, crop_h))

                bx, by, bw, bh = ann["bbox"]
                self.assertGreaterEqual(bx, 0.0)
                self.assertGreaterEqual(by, 0.0)
                self.assertGreater(bw, 0.0)
                self.assertGreater(bh, 0.0)
                self.assertLessEqual(bx + bw, float(crop_w))
                self.assertLessEqual(by + bh, float(crop_h))

                for poly in ann["segmentation"]:
                    xs = poly[0::2]
                    ys = poly[1::2]
                    for x in xs:
                        self.assertGreaterEqual(x, 0.0)
                        self.assertLessEqual(x, float(crop_w))
                    for y in ys:
                        self.assertGreaterEqual(y, 0.0)
                        self.assertLessEqual(y, float(crop_h))

                # Also verify local polygon coordinates match (x_sheet - x0, y_sheet - y0)
                if im_info["file_name"].endswith("_inst0.jpg"):
                    self.assertEqual(ann["category_id"], 1)
                    # x_min=120 in sheet space -> 120 - 108 = 12 in local patch space
                    # y_min=100 in sheet space -> 100 - 82 = 18 in local patch space
                    self.assertAlmostEqual(bx, 12.0, delta=1.0)
                    self.assertAlmostEqual(by, 18.0, delta=1.0)
                elif im_info["file_name"].endswith("_inst1.jpg"):
                    self.assertEqual(ann["category_id"], 2)

    def test_crop_boundary_clamping_near_sheet_edges(self):
        """
        Verifies that when a leaf mask touches or approaches the herbarium sheet border,
        the 12% expanded crop window clamps to [0, sheet_w] x [0, sheet_h] and local
        polygon coordinates remain strictly within [0, crop_w] x [0, crop_h].
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_p = Path(tmp_dir)
            images_dir = tmp_p / "vouchers"
            masks_dir = tmp_p / "masks"
            crops_dir = tmp_p / "crops"
            images_dir.mkdir()
            masks_dir.mkdir()

            sheet_h, sheet_w = 300, 400
            sheet_img = np.full((sheet_h, sheet_w, 3), 200, dtype=np.uint8)
            cv2.imwrite(str(images_dir / "EDGE_SHEET.jpg"), sheet_img)

            # Leaf right against top-left corner (2, 2) to (80, 60)
            mask_tl = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
            cv2.rectangle(mask_tl, (2, 2), (80, 60), 255, -1)
            cv2.imwrite(str(masks_dir / "EDGE_SHEET_inst00_basal_leaf_whole.png"), mask_tl)

            # Leaf right against bottom-right corner (320, 240) to (399, 299)
            mask_br = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
            cv2.rectangle(mask_br, (320, 240), (398, 298), 255, -1)
            cv2.imwrite(str(masks_dir / "EDGE_SHEET_inst01_basal_leaf_partial.png"), mask_br)

            coco_doc = export_coco_from_yolo_and_masks(
                masks_dir=masks_dir,
                images_dir=images_dir,
                crops_dir=crops_dir,
                min_area_px=10.0,
            )

            self.assertEqual(len(coco_doc["images"]), 2)
            img_map = {im["id"]: im for im in coco_doc["images"]}

            for ann in coco_doc["annotations"]:
                im_info = img_map[ann["image_id"]]
                crop_w, crop_h = im_info["width"], im_info["height"]
                self.assertLessEqual(crop_w, sheet_w)
                self.assertLessEqual(crop_h, sheet_h)
                for poly in ann["segmentation"]:
                    for x in poly[0::2]:
                        self.assertTrue(0.0 <= x <= float(crop_w))
                    for y in poly[1::2]:
                        self.assertTrue(0.0 <= y <= float(crop_h))

    def test_yolo_polygon_fallback_and_full_sheet_flag(self):
        """
        Verifies that `export_coco_from_yolo_and_masks` can also export from normalized
        YOLO polygon `.txt` files when PNG masks are absent, and honors `--full-sheet`.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_p = Path(tmp_dir)
            images_dir = tmp_p / "vouchers"
            yolo_dir = tmp_p / "raw_annotations"
            crops_dir = tmp_p / "crops"
            images_dir.mkdir()
            yolo_dir.mkdir()

            sheet_h, sheet_w = 500, 400
            cv2.imwrite(
                str(images_dir / "V_YOLO.jpg"),
                np.ones((sheet_h, sheet_w, 3), dtype=np.uint8) * 128,
            )

            # Class 0 (basal_leaf_whole): normalized rectangle [0.2, 0.2] to [0.5, 0.6]
            yolo_txt = yolo_dir / "V_YOLO.txt"
            yolo_txt.write_text("0 0.2 0.2 0.5 0.2 0.5 0.6 0.2 0.6\n", encoding="utf-8")

            coco_crop = export_coco_from_yolo_and_masks(
                yolo_dir=yolo_dir,
                images_dir=images_dir,
                crops_dir=crops_dir,
                min_area_px=10.0,
                full_sheet=False,
            )
            self.assertEqual(len(coco_crop["images"]), 1)
            self.assertEqual(coco_crop["images"][0]["file_name"], "V_YOLO_inst0.jpg")
            self.assertLess(coco_crop["images"][0]["width"], sheet_w)
            self.assertLess(coco_crop["images"][0]["height"], sheet_h)
            self.assertTrue((crops_dir / "V_YOLO_inst0.jpg").exists())

            coco_full = export_coco_from_yolo_and_masks(
                yolo_dir=yolo_dir,
                images_dir=images_dir,
                min_area_px=10.0,
                full_sheet=True,
            )
            self.assertEqual(len(coco_full["images"]), 1)
            self.assertEqual(coco_full["images"][0]["file_name"], "V_YOLO.jpg")
            self.assertEqual(coco_full["images"][0]["width"], sheet_w)
            self.assertEqual(coco_full["images"][0]["height"], sheet_h)


if __name__ == "__main__":
    unittest.main()
