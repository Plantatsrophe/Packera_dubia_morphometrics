"""
Tests for Precision SAM 2 Botanical Annotator
==============================================
Verifies window configurations, multi-modal exclusion point handling,
bounding box prompting, freehand lasso conversion, knife cut geometry,
YOLO polygon extraction, zoom/pan clamping, hard negative background saving,
voucher skipping, and HUD mode state transitions.
"""

import sys
import unittest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:
    from scripts.data_prep.annotate_with_sam2 import (
        PrecisionSAM2Annotator,
        CLASS_NAMES,
        CLASS_COLORS,
        get_project_root
    )
except ImportError:
    from scripts.annotate_with_sam2 import (
        PrecisionSAM2Annotator,
        CLASS_NAMES,
        CLASS_COLORS,
        get_project_root
    )


class TestPrecisionSAM2Annotator(unittest.TestCase):
    def setUp(self):
        self.temp_dir = Path(tempfile.mkdtemp())
        self.vouchers_dir = self.temp_dir / "data" / "raw_vouchers"
        self.vouchers_dir.mkdir(parents=True, exist_ok=True)
        self.annotations_dir = self.temp_dir / "data" / "raw_annotations"
        self.annotations_dir.mkdir(parents=True, exist_ok=True)

        # Create dummy test voucher image
        self.test_img = np.zeros((1000, 800, 3), dtype=np.uint8)
        cv2.ellipse(self.test_img, (400, 300), (100, 200), 0, 0, 360, (0, 180, 0), -1)
        cv2.ellipse(self.test_img, (400, 700), (80, 150), 0, 0, 360, (30, 60, 120), -1)
        self.test_img_path = self.vouchers_dir / "NCU00099999.jpg"
        cv2.imwrite(str(self.test_img_path), self.test_img)

        with patch("scripts.data_prep.annotate_with_sam2.build_sam2") as mock_build_sam2, \
             patch("scripts.data_prep.annotate_with_sam2.SAM2ImagePredictor") as mock_predictor_cls:
            
            mock_predictor = MagicMock()
            mock_predictor.predict.return_value = (
                np.ones((1, 1000, 800), dtype=bool),
                np.array([0.99]),
                None
            )
            mock_predictor_cls.return_value = mock_predictor

            self.annotator = PrecisionSAM2Annotator(
                images_dir=self.vouchers_dir,
                output_dir=self.annotations_dir,
                checkpoint_path="dummy.pt",
                config_path="dummy.yaml",
                window_w=800,
                window_h=600
            )
            self.annotator.current_img = self.test_img
            self.annotator.orig_h, self.annotator.orig_w = 1000, 800

    def tearDown(self):
        shutil.rmtree(self.temp_dir, ignore_errors=True)

    def test_project_root_resolution(self):
        root = get_project_root()
        self.assertTrue(root.exists())
        self.assertTrue((root / "scripts").exists())

    def test_annotator_initialization(self):
        self.assertEqual(len(self.annotator.image_files), 1)
        self.assertEqual(self.annotator.image_files[0].stem, "NCU00099999")
        self.assertFalse(self.annotator.exclusion_mode)
        self.assertFalse(self.annotator.lasso_mode)
        self.assertEqual(self.annotator.zoom_level, 1.0)
        self.assertEqual(self.annotator.pan_offset, [0, 0])

    def test_auto_resume_unannotated_voucher(self):
        # Create second voucher image in test directory
        img2_path = self.vouchers_dir / "NCU00088888.jpg"
        cv2.imwrite(str(img2_path), self.test_img)

        # Create annotation file for the first voucher (NCU00099999.txt)
        (self.annotations_dir / "NCU00099999.txt").write_text("0 0.1 0.1 0.2 0.1 0.2 0.2 0.1 0.2\n")

        with patch("scripts.data_prep.annotate_with_sam2.build_sam2") as mock_build_sam2, \
             patch("scripts.data_prep.annotate_with_sam2.SAM2ImagePredictor") as mock_predictor_cls:
            annotator = PrecisionSAM2Annotator(
                images_dir=self.vouchers_dir,
                output_dir=self.annotations_dir,
                checkpoint_path="dummy.pt",
                config_path="dummy.yaml"
            )
            # Should auto-advance to index of NCU00088888 because NCU00099999 is already annotated
            self.assertEqual(annotator.image_files[annotator.current_idx].stem, "NCU00088888")

    def test_coordinate_mapping(self):
        # 1.0x Zoom at window 800x600 mapping to 800x1000 image
        orig_x, orig_y = self.annotator.get_orig_coords(400, 300, 800, 600)
        self.assertEqual(orig_x, 400)
        self.assertEqual(orig_y, 500)

        # 2.0x Zoom centered with pan offset
        self.annotator.zoom_level = 2.0
        self.annotator.pan_offset = [200, 250]
        orig_x, orig_y = self.annotator.get_orig_coords(400, 300, 800, 600)
        self.assertEqual(orig_x, 400)
        self.assertEqual(orig_y, 500)

    def test_zoom_and_clamping(self):
        # Zoom in by 2.0x anchored at (400, 500)
        self.annotator.zoom(2.0, center_x=400, center_y=500)
        self.assertEqual(self.annotator.zoom_level, 2.0)
        self.assertEqual(self.annotator.pan_offset, [200, 250])

        # Zoom exceeding 6.0x should clamp to 6.0x
        self.annotator.zoom(10.0, center_x=400, center_y=500)
        self.assertEqual(self.annotator.zoom_level, 6.0)

        # Zoom out below 1.0x should clamp to 1.0x and reset pan
        self.annotator.zoom(0.01)
        self.assertEqual(self.annotator.zoom_level, 1.0)
        self.assertEqual(self.annotator.pan_offset, [0, 0])

    def test_pan_controls(self):
        # Pan at 1.0x zoom automatically triggers 1.5x magnification
        self.annotator.zoom_level = 1.0
        self.annotator.pan_offset = [0, 0]
        self.annotator.pan(dy_frac=0.20)
        self.assertGreater(self.annotator.zoom_level, 1.0)
        self.assertGreater(self.annotator.pan_offset[1], 0)


    def test_multi_modal_exclusion_points(self):
        window_dims = (800, 600)

        # 1. Standard Left Click in INCLUDE mode -> Label 1 (Positive)
        self.annotator.exclusion_mode = False
        self.annotator.mouse_callback(cv2.EVENT_LBUTTONDOWN, 400, 300, 0, window_dims)
        self.assertEqual(len(self.annotator.prompt_points), 1)
        self.assertEqual(self.annotator.prompt_labels[-1], 1)

        # 2. Native Right-Click -> Label 0 (Negative Exclusion)
        self.annotator.mouse_callback(cv2.EVENT_RBUTTONDOWN, 400, 450, 0, window_dims)
        self.assertEqual(len(self.annotator.prompt_points), 2)
        self.assertEqual(self.annotator.prompt_labels[-1], 0)

        # 3. Alt + Left-Click -> Label 0 (Negative Exclusion)
        self.annotator.mouse_callback(cv2.EVENT_LBUTTONDOWN, 400, 500, cv2.EVENT_FLAG_ALTKEY, window_dims)
        self.assertEqual(len(self.annotator.prompt_points), 3)
        self.assertEqual(self.annotator.prompt_labels[-1], 0)

        # 4. Toggle Exclusion Mode ('e') -> Left Click produces Label 0
        self.annotator.exclusion_mode = True
        self.annotator.mouse_callback(cv2.EVENT_LBUTTONDOWN, 400, 550, 0, window_dims)
        self.assertEqual(len(self.annotator.prompt_points), 4)
        self.assertEqual(self.annotator.prompt_labels[-1], 0)

    def test_freehand_lasso_tool(self):
        window_dims = (800, 600)
        self.annotator.lasso_mode = True
        self.annotator.exclusion_mode = False

        # Start Lasso Drag
        self.annotator.mouse_callback(cv2.EVENT_LBUTTONDOWN, 100, 100, 0, window_dims)
        self.assertTrue(self.annotator.is_drawing_lasso)

        # Mouse Moves along triangle path
        self.annotator.mouse_callback(cv2.EVENT_MOUSEMOVE, 200, 100, 0, window_dims)
        self.annotator.mouse_callback(cv2.EVENT_MOUSEMOVE, 200, 200, 0, window_dims)
        self.annotator.mouse_callback(cv2.EVENT_MOUSEMOVE, 100, 200, 0, window_dims)

        # Release Mouse -> Fill polygon
        self.annotator.mouse_callback(cv2.EVENT_LBUTTONUP, 100, 100, 0, window_dims)
        self.assertFalse(self.annotator.is_drawing_lasso)
        self.assertIsNotNone(self.annotator.current_candidate_mask)
        self.assertTrue(self.annotator.current_candidate_mask.dtype == bool)
        # Inside of the triangle (mapped coordinates) should be True
        poly = self.annotator.mask_to_yolo_polygon(self.annotator.current_candidate_mask)
        self.assertGreaterEqual(len(poly), 6)

    def test_bounding_box_prompt(self):
        window_dims = (800, 600)

        # Shift + Left Click Down
        self.annotator.mouse_callback(cv2.EVENT_LBUTTONDOWN, 100, 100, cv2.EVENT_FLAG_SHIFTKEY, window_dims)
        self.assertTrue(self.annotator.is_drawing_box)

        # Mouse Move
        self.annotator.mouse_callback(cv2.EVENT_MOUSEMOVE, 500, 400, cv2.EVENT_FLAG_SHIFTKEY, window_dims)
        self.assertTrue(self.annotator.is_drawing_box)

        # Left Click Up
        self.annotator.mouse_callback(cv2.EVENT_LBUTTONUP, 500, 400, 0, window_dims)
        self.assertFalse(self.annotator.is_drawing_box)
        self.assertIsNotNone(self.annotator.prompt_box)
        x1, y1, x2, y2 = self.annotator.prompt_box
        self.assertTrue(x1 < x2 and y1 < y2)
        self.assertEqual(x1, 100)
        self.assertEqual(x2, 500)

    def test_interactive_knife_cutter(self):
        # Create a candidate mask with two connected components (petiole and root)
        candidate_mask = np.zeros((1000, 800), dtype=bool)
        candidate_mask[200:400, 350:450] = True  # Component 1 (Petiole)
        candidate_mask[400:600, 350:450] = True  # Connected to root
        self.annotator.current_candidate_mask = candidate_mask
        self.annotator.prompt_points = [[400, 300]]  # Point inside Component 1
        self.annotator.prompt_labels = [1]

        # Apply knife cut across y = 400
        self.annotator.apply_knife_cut((300, 400), (500, 400))

        self.assertEqual(len(self.annotator.cut_lines), 1)
        # Check that cut carved zeros across y=400
        self.assertFalse(self.annotator.current_candidate_mask[400, 400])
        # Check that component containing (400, 300) is preserved
        self.assertTrue(self.annotator.current_candidate_mask[300, 400])

    def test_mask_to_yolo_polygon(self):
        mask = np.zeros((1000, 800), dtype=bool)
        mask[100:300, 100:300] = True
        poly = self.annotator.mask_to_yolo_polygon(mask)

        self.assertGreaterEqual(len(poly), 6)
        # Coordinates must be normalized between 0.0 and 1.0
        for v in poly:
            self.assertTrue(0.0 <= v <= 1.0)

    def test_save_current_sheet_with_instances(self):
        self.annotator.saved_instances = [
            {"class_id": 0, "polygon": [0.1, 0.1, 0.2, 0.1, 0.2, 0.2, 0.1, 0.2]},
            {"class_id": 1, "polygon": [0.3, 0.3, 0.4, 0.3, 0.4, 0.4, 0.3, 0.4]},
            {"class_id": 4, "polygon": [0.5, 0.5, 0.6, 0.5, 0.6, 0.6, 0.5, 0.6]}
        ]
        self.annotator.save_current_sheet()

        label_file = self.annotator.output_dir / "NCU00099999.txt"
        self.assertTrue(label_file.exists())

        with open(label_file, "r") as f:
            lines = [line.strip() for line in f.readlines() if line.strip()]

        self.assertEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("0 "))
        self.assertTrue(lines[1].startswith("1 "))
        self.assertTrue(lines[2].startswith("4 "))

    def test_save_current_sheet_hard_negative_empty_sample(self):
        # 0 instances saved -> exports empty .txt hard negative
        self.annotator.saved_instances = []
        self.annotator.save_current_sheet()

        label_file = self.annotator.output_dir / "NCU00099999.txt"
        self.assertTrue(label_file.exists())
        self.assertEqual(label_file.stat().st_size, 0)

    def test_hud_rendering_and_mode_display(self):
        self.annotator.exclusion_mode = False
        self.annotator.lasso_mode = False
        display_inc = self.annotator.render_display(800, 600)
        self.assertEqual(display_inc.shape, (600, 800, 3))

        self.annotator.exclusion_mode = True
        self.annotator.lasso_mode = False
        display_exc = self.annotator.render_display(800, 600)
        self.assertEqual(display_exc.shape, (600, 800, 3))

    def test_load_existing_annotations(self):
        # Write dummy YOLO annotation file
        label_file = self.annotator.output_dir / "NCU00099999.txt"
        with open(label_file, "w") as f:
            f.write("0 0.1 0.1 0.2 0.1 0.2 0.2 0.1 0.2\n")
            f.write("2 0.3 0.3 0.4 0.3 0.4 0.4 0.3 0.4\n")

        self.annotator.load_existing_annotations("NCU00099999")
        self.assertEqual(len(self.annotator.saved_instances), 2)
        self.assertEqual(self.annotator.saved_instances[0]["class_id"], 0)
        self.assertEqual(self.annotator.saved_instances[1]["class_id"], 2)


if __name__ == "__main__":
    unittest.main()

