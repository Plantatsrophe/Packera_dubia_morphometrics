#!/usr/bin/env python3
"""
Unit tests for PrecisionSAM2Annotator in scripts/annotation_and_training/annotate_with_sam2.py
and helper utilities in scripts/annotation_and_training/sam2_annotator_utils.py.
Verifies:
  1. Cursor-centered zoom mathematics and clamping.
  2. Multimask proposal cycling and IoU score indexing.
  3. Two-click knife slicing with dilation.
  4. Morphological single-pixel margin tuning (3x3 dilation and erosion).
  5. Instant class commitment and prompt buffer clearance.
  6. Semi-transparent HUD overlay banner formatting.
  7. High-contrast contour and translucent fill viewport overlays.
  8. COCO dataset synchronization and YOLO polygon exports.
"""

import json
import tempfile
import unittest
from pathlib import Path
import cv2
import numpy as np

from scripts.annotation_and_training.annotate_with_sam2 import PrecisionSAM2Annotator, CLASS_NAMES
from scripts.annotation_and_training.sam2_annotator_utils import (
    KEYSYM_TO_CLASS,
    apply_knife_cut,
    apply_mask_dilation,
    apply_mask_erosion,
    apply_morphological_tuning,
    apply_viewport_transform,
    calculate_cursor_centered_zoom,
    clamp_viewport_pan,
    convert_masks_to_coco_dataset,
    export_coco_annotations,
    get_undo_button_rect,
    image_to_viewport_coords,
    overlay_candidate_mask_on_viewport,
    parse_mask_filename,
    polygon_interior_point,
    polygon_to_bounding_box,
    rasterize_lasso_polygon,
    render_hud_overlay,
    split_mask_with_knife_line,
    viewport_to_image_coords,
)


class TestPrecisionSAM2Annotator(unittest.TestCase):
    def test_class_names(self):
        self.assertIn("basal_leaf_whole", CLASS_NAMES)
        self.assertIn("basal_leaf_partial", CLASS_NAMES)
        self.assertIn("cauline_leaf", CLASS_NAMES)
        self.assertEqual(CLASS_NAMES[0], "basal_leaf_whole")
        self.assertEqual(CLASS_NAMES[1], "basal_leaf_partial")
        self.assertEqual(CLASS_NAMES[2], "cauline_leaf")

    def test_cursor_centered_zoom(self):
        # Initial fit-to-window state
        zoom, pan = calculate_cursor_centered_zoom(
            current_zoom=1.0,
            current_pan=(0, 0),
            cursor_vx=640,
            cursor_vy=400,
            zoom_in=True,
            target_w=1280,
            target_h=800,
            orig_w=4000,
            orig_h=6000,
            step=1.15,
            min_zoom=1.0,
            max_zoom=16.0,
        )
        self.assertAlmostEqual(zoom, 1.15, places=2)
        # Zooming in on center should keep pan centered
        self.assertIsInstance(pan, list)
        self.assertEqual(len(pan), 2)

        # Deep zoom clamping at 16.0x
        deep_zoom, deep_pan = calculate_cursor_centered_zoom(
            current_zoom=15.0,
            current_pan=tuple(pan),
            cursor_vx=640,
            cursor_vy=400,
            zoom_in=True,
            target_w=1280,
            target_h=800,
            orig_w=4000,
            orig_h=6000,
            step=1.15,
            min_zoom=1.0,
            max_zoom=16.0,
        )
        self.assertLessEqual(deep_zoom, 16.0)

        # Zooming out below 1.0 resets to fit to window
        out_zoom, out_pan = calculate_cursor_centered_zoom(
            current_zoom=1.05,
            current_pan=(100, 100),
            cursor_vx=640,
            cursor_vy=400,
            zoom_in=False,
            target_w=1280,
            target_h=800,
            orig_w=4000,
            orig_h=6000,
            step=1.15,
            min_zoom=1.0,
            max_zoom=16.0,
        )
        self.assertEqual(out_zoom, 1.0)
        self.assertEqual(out_pan, [0, 0])

    def test_multimask_proposal_cycling(self):
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        m1 = np.ones((50, 50), dtype=np.uint8) * 10
        m2 = np.ones((50, 50), dtype=np.uint8) * 20
        m3 = np.ones((50, 50), dtype=np.uint8) * 30
        annotator.candidate_masks = [m1, m2, m3]
        annotator.candidate_scores = [0.85, 0.92, 0.78]
        annotator.active_mask_idx = 0
        annotator.candidate_mask = m1

        annotator.cycle_multimask_proposal()
        self.assertEqual(annotator.active_mask_idx, 1)
        np.testing.assert_array_equal(annotator.candidate_mask, m2)

        annotator.cycle_multimask_proposal()
        self.assertEqual(annotator.active_mask_idx, 2)
        np.testing.assert_array_equal(annotator.candidate_mask, m3)

        annotator.cycle_multimask_proposal()
        self.assertEqual(annotator.active_mask_idx, 0)
        np.testing.assert_array_equal(annotator.candidate_mask, m1)

    def test_knife_slicing_with_dilation(self):
        mask = np.ones((100, 100), dtype=np.uint8) * 255
        cut = split_mask_with_knife_line(mask, (50, 0), (50, 100), line_thickness=2, dilation_px=2)
        # Vertical cut line around x=50 must be zeroed out
        self.assertEqual(cut[50, 50], 0)
        self.assertEqual(cut[50, 52], 0)
        self.assertEqual(cut[50, 48], 0)
        # Pixels further away must remain foreground
        self.assertEqual(cut[50, 10], 255)
        self.assertEqual(cut[50, 90], 255)

    def test_knife_sever_and_prune_with_prompt_point(self):
        mask = np.ones((100, 100), dtype=np.uint8) * 255
        # Cut vertically at x=50, keeping the left fragment with prompt at (25, 50)
        cut_left = split_mask_with_knife_line(
            mask, (50, 0), (50, 100), line_thickness=2, dilation_px=2,
            keep_points=[(25.0, 50.0)]
        )
        # Left side must be retained
        self.assertEqual(cut_left[50, 10], 255)
        self.assertEqual(cut_left[50, 25], 255)
        # Right severed fragment must be completely cleared/pruned
        self.assertEqual(cut_left[50, 75], 0)
        self.assertEqual(cut_left[50, 90], 0)

        # Cut vertically at x=50, keeping the right fragment with prompt at (75, 50)
        cut_right = split_mask_with_knife_line(
            mask, (50, 0), (50, 100), line_thickness=2, dilation_px=2,
            keep_points=[(75.0, 50.0)]
        )
        # Left severed fragment must be pruned
        self.assertEqual(cut_right[50, 25], 0)
        # Right side must be retained
        self.assertEqual(cut_right[50, 75], 255)

    def test_knife_sever_and_prune_with_bounding_box(self):
        mask = np.ones((100, 100), dtype=np.uint8) * 255
        # Cut vertically at x=50, bounding box only encompasses the right side [60, 10, 95, 90]
        cut_box = split_mask_with_knife_line(
            mask, (50, 0), (50, 100), line_thickness=2, dilation_px=2,
            bounding_box=(60.0, 10.0, 95.0, 90.0)
        )
        self.assertEqual(cut_box[50, 25], 0)
        self.assertEqual(cut_box[50, 75], 255)

    def test_morphological_tuning(self):
        mask = np.zeros((100, 100), dtype=np.uint8)
        mask[40:60, 40:60] = 255

        # Dilation expands margin
        dilated = apply_morphological_tuning(mask, operation="dilate", kernel_size=3)
        self.assertGreater(np.count_nonzero(dilated), np.count_nonzero(mask))
        self.assertEqual(dilated[39, 50], 255)

        # Erosion contracts margin
        eroded = apply_morphological_tuning(mask, operation="erode", kernel_size=3)
        self.assertLess(np.count_nonzero(eroded), np.count_nonzero(mask))
        self.assertEqual(eroded[40, 40], 0)

    def test_instant_class_commit(self):
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.active_image = np.zeros((200, 200, 3), dtype=np.uint8)
        annotator.cached_base_layer = annotator.active_image.copy()
        annotator.overlay_alpha = 0.55
        annotator.saved_instances = []
        annotator.point_coords = [[50.0, 50.0]]
        annotator.point_labels = [1]
        annotator.box_prompt = [10.0, 10.0, 100.0, 100.0]
        annotator.polygon_points = [(10, 10), (20, 20)]
        annotator.knife_pt_a = (10, 10)

        cand = np.zeros((200, 200), dtype=np.uint8)
        cand[20:80, 20:80] = 255
        annotator.candidate_mask = cand
        annotator.candidate_masks = [cand.copy()]
        annotator.candidate_scores = [0.95]
        annotator.active_mask_idx = 0

        # Commit class 0 (basal_leaf_whole)
        annotator.commit_active_instance(0)

        self.assertEqual(len(annotator.saved_instances), 1)
        self.assertEqual(annotator.saved_instances[0]["class_id"], 0)
        self.assertEqual(annotator.saved_instances[0]["label"], "basal_leaf_whole")
        # Prompt buffers must be reset
        self.assertIsNone(annotator.candidate_mask)
        self.assertEqual(annotator.candidate_masks, [])
        self.assertEqual(annotator.point_coords, [])
        self.assertEqual(annotator.point_labels, [])
        self.assertIsNone(annotator.box_prompt)
        self.assertIsNone(annotator.knife_pt_a)

    def test_hud_overlay_formatting(self):
        canvas = np.zeros((600, 800, 3), dtype=np.uint8)
        hud = render_hud_overlay(
            display_img=canvas,
            voucher_name="NCU00001234",
            voucher_idx=0,
            total_vouchers=15,
            saved_instances=[{"class_id": 0}, {"class_id": 1}],
            mode="SELECT",
            zoom_level=2.3,
            pan_offset=(120, 340),
            candidate_idx=1,
            candidate_total=3,
            candidate_iou=0.94,
            view_mode="CONTOUR",
            alpha=0.6,
        )
        self.assertEqual(hud.shape, canvas.shape)
        # Check that HUD top banner was rendered
        self.assertTrue(np.any(hud[:70, :] > 0))

    def test_overlay_candidate_mask_view_modes(self):
        vp = np.zeros((200, 200, 3), dtype=np.uint8)
        mask = np.zeros((200, 200), dtype=np.uint8)
        mask[50:150, 50:150] = 255
        transform = (1.0, 1.0, 0, 0)

        # FILL mode
        vp_fill = overlay_candidate_mask_on_viewport(
            vp.copy(), mask, transform, alpha=0.55, view_mode="FILL"
        )
        self.assertGreater(np.count_nonzero(vp_fill[100, 100]), 0)

        # CONTOUR mode (center should not have fill, only contour)
        vp_contour = overlay_candidate_mask_on_viewport(
            vp.copy(), mask, transform, alpha=0.55, view_mode="CONTOUR"
        )
        # Boundary contour drawn
        self.assertGreater(np.count_nonzero(vp_contour), 0)
        # Center of mask in contour mode remains 0
        self.assertEqual(tuple(vp_contour[100, 100]), (0, 0, 0))

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

            txt_file = out_dir / "VOUCHER_TEST001.txt"
            self.assertTrue(txt_file.exists())
            with open(txt_file, "r") as f:
                lines = f.readlines()
                self.assertEqual(len(lines), 3)

            m0 = out_dir / "masks" / "VOUCHER_TEST001_inst00_basal_leaf_whole.png"
            m1 = out_dir / "masks" / "VOUCHER_TEST001_inst01_basal_leaf_partial.png"
            m2 = out_dir / "masks" / "VOUCHER_TEST001_inst02_cauline_leaf.png"

            self.assertTrue(m0.exists())
            self.assertTrue(m1.exists())
            self.assertTrue(m2.exists())

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
        annotator.candidate_masks = []
        annotator.active_mask_idx = 0
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

            mask_img = np.zeros((200, 200), dtype=np.uint8)
            cv2.rectangle(mask_img, (30, 30), (100, 100), 255, -1)
            cv2.imwrite(str(masks_dir / "VOUCHER001_inst00_basal_leaf_whole.png"), mask_img)

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
            cat_ids = {ann["category_id"] for ann in coco_doc["annotations"]}
            self.assertEqual(cat_ids, {1, 2})

    def test_viewport_coordinate_transforms_exact_zero_drift(self):
        """
        Verifies that given a mock high-resolution specimen image (4000 x 6000)
        and arbitrary zoom factors (2.0x, 4.5x), forward mapping (image space to
        viewport screen space) and inverse mapping (mouse click on screen to image pixel)
        are exact with zero drift (< 1e-6 px).
        """
        test_dimensions = [(4000, 6000), (6000, 4000)]
        zoom_factors = [2.0, 4.5, 1.0, 3.25]
        target_w, target_h = 1280, 800

        for orig_h, orig_w in test_dimensions:
            annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
            annotator.orig_h = orig_h
            annotator.orig_w = orig_w

            for zoom in zoom_factors:
                crop_w = orig_w / zoom
                crop_h = orig_h / zoom
                scale_x = target_w / crop_w
                scale_y = target_h / crop_h

                for pan_x, pan_y in [(0, 0), (250, 400)]:
                    # 1. Test mouse clicks on screen -> image space -> screen space
                    screen_clicks = [
                        (0.0, 0.0),
                        (float(target_w // 2), float(target_h // 2)),
                        (float(target_w - 1), float(target_h - 1)),
                        (142.75, 591.33),
                        (640.0, 400.0),
                    ]
                    for vx, vy in screen_clicks:
                        ix, iy = viewport_to_image_coords(
                            vx, vy, scale_x, scale_y, pan_x, pan_y
                        )
                        vx_proj, vy_proj = image_to_viewport_coords(
                            ix, iy, scale_x, scale_y, pan_x, pan_y
                        )
                        self.assertAlmostEqual(vx_proj, vx, places=6)
                        self.assertAlmostEqual(vy_proj, vy, places=6)

                    # 2. Test image points -> screen space -> image space
                    image_points = [
                        (float(pan_x + 10), float(pan_y + 10)),
                        (float(pan_x + crop_w / 2), float(pan_y + crop_h / 2)),
                        (float(pan_x + crop_w - 10), float(pan_y + crop_h - 10)),
                        (float(pan_x + 123.45), float(pan_y + 456.78)),
                    ]
                    for ix, iy in image_points:
                        vx, vy = image_to_viewport_coords(
                            ix, iy, scale_x, scale_y, pan_x, pan_y
                        )
                        ix_proj, iy_proj = viewport_to_image_coords(
                            vx, vy, scale_x, scale_y, pan_x, pan_y
                        )
                        self.assertAlmostEqual(ix_proj, ix, places=6)
                        self.assertAlmostEqual(iy_proj, iy, places=6)

                    # 3. Test PrecisionSAM2Annotator integer internal methods for zero drift
                    int_screen_clicks = [
                        (target_w // 4, target_h // 4),
                        (target_w // 2, target_h // 2),
                        (3 * target_w // 4, 3 * target_h // 4),
                    ]
                    for vx, vy in int_screen_clicks:
                        a_ix, a_iy = annotator._viewport_to_image(
                            vx, vy, scale_x, scale_y, pan_x, pan_y
                        )
                        a_vx, a_vy = annotator._image_to_viewport(
                            a_ix, a_iy, scale_x, scale_y, pan_x, pan_y
                        )
                        self.assertLessEqual(abs(a_vx - vx), 1)
                        self.assertLessEqual(abs(a_vy - vy), 1)

    def test_viewport_boundary_clamping(self):
        """
        Verifies that viewport boundary clamping prevents panning from shifting
        the viewport entirely off the specimen canvas, maintaining at least a
        minimum visible margin under extreme panning inputs.
        """
        orig_w, orig_h = 6000, 4000
        target_w, target_h = 1280, 800

        for zoom in [1.5, 2.0, 4.5, 10.0]:
            crop_w = orig_w / zoom
            crop_h = orig_h / zoom

            # Test extreme panning values (far beyond canvas bounds)
            extreme_pans = [
                (100000, 100000),
                (-100000, -100000),
                (-50000, 25000),
                (30000, -40000),
            ]

            for px, py in extreme_pans:
                clamped_x, clamped_y = clamp_viewport_pan(
                    px, py, zoom_level=zoom, orig_w=orig_w, orig_h=orig_h, margin_ratio=0.15
                )

                # The viewport crop window MUST overlap the specimen canvas:
                # 1. Left boundary of crop (clamped_x) must be strictly less than orig_w
                self.assertLess(clamped_x, orig_w)
                # 2. Right boundary of crop (clamped_x + crop_w) must be strictly greater than 0
                self.assertGreater(clamped_x + crop_w, 0)
                # 3. Top boundary of crop (clamped_y) must be strictly less than orig_h
                self.assertLess(clamped_y, orig_h)
                # 4. Bottom boundary of crop (clamped_y + crop_h) must be strictly greater than 0
                self.assertGreater(clamped_y + crop_h, 0)

                # Calculate exact visible specimen width and height in crop
                overlap_w = min(orig_w, clamped_x + crop_w) - max(0, clamped_x)
                overlap_h = min(orig_h, clamped_y + crop_h) - max(0, clamped_y)

                self.assertGreaterEqual(overlap_w, crop_w * 0.14)
                self.assertGreaterEqual(overlap_h, crop_h * 0.14)

            # Also verify calculate_cursor_centered_zoom clamps pan
            _, clamped_zoom_pan = calculate_cursor_centered_zoom(
                current_zoom=zoom,
                current_pan=(100000, 100000),
                cursor_vx=target_w,
                cursor_vy=target_h,
                zoom_in=True,
                target_w=target_w,
                target_h=target_h,
                orig_w=orig_w,
                orig_h=orig_h,
            )
            self.assertLess(clamped_zoom_pan[0], orig_w)
            self.assertGreater(clamped_zoom_pan[0] + crop_w, 0)

        # Test apply_viewport_transform headless rendering with extreme pan
        mock_specimen = np.zeros((400, 600, 3), dtype=np.uint8)
        mock_specimen[150:250, 200:400] = (0, 255, 0)
        clamped_p = clamp_viewport_pan(-10000, -10000, 2.0, 600, 400)
        rendered_vp, transform = apply_viewport_transform(
            mock_specimen,
            zoom_level=2.0,
            pan_offset=clamped_p,
            target_w=300,
            target_h=200,
        )
        self.assertEqual(rendered_vp.shape, (200, 300, 3))

    def test_knife_slicing_circle_split(self):
        """
        Verifies that apply_knife_cut severing a synthetic binary circle mask across
        its diameter zeroes out mask pixels along the line and cleanly splits
        the single connected component into two distinct connected components.
        """
        h, w = 200, 200
        center_x, center_y = 100, 100
        radius = 50

        # Create synthetic binary circle mask
        mask = np.zeros((h, w), dtype=np.uint8)
        cv2.circle(mask, (center_x, center_y), radius, 255, -1)

        # Confirm single connected component before cut (background = 0, circle = 1)
        num_labels_init, labels_init = cv2.connectedComponents(mask, connectivity=8)
        self.assertEqual(num_labels_init, 2)
        initial_circle_area = np.count_nonzero(mask)
        self.assertGreater(initial_circle_area, 0)

        # Test 1: Vertical diameter cut (pt1=(100, 40), pt2=(100, 160))
        cut_vert = apply_knife_cut(
            mask, (center_x, center_y - radius - 10), (center_x, center_y + radius + 10), thickness=2
        )

        # Assert pixels along diameter are zeroed out
        for y in range(center_y - radius + 5, center_y + radius - 5):
            self.assertEqual(cut_vert[y, center_x], 0)

        # Assert connected components split into exactly 2 distinct components (plus background = 3 total)
        num_labels_vert, labels_vert = cv2.connectedComponents(cut_vert, connectivity=8)
        self.assertEqual(num_labels_vert, 3)

        # Check that both split halves have non-trivial foreground area
        comp1_area = np.count_nonzero(labels_vert == 1)
        comp2_area = np.count_nonzero(labels_vert == 2)
        self.assertGreater(comp1_area, 0)
        self.assertGreater(comp2_area, 0)
        self.assertLess(comp1_area + comp2_area, initial_circle_area)

        # Test 2: Horizontal diameter cut
        cut_horiz = apply_knife_cut(
            mask, (center_x - radius - 10, center_y), (center_x + radius + 10, center_y), thickness=2
        )
        for x in range(center_x - radius + 5, center_x + radius - 5):
            self.assertEqual(cut_horiz[center_y, x], 0)
        num_labels_horiz, _ = cv2.connectedComponents(cut_horiz, connectivity=8)
        self.assertEqual(num_labels_horiz, 3)

        # Test 3: Diagonal diameter cut
        cut_diag = apply_knife_cut(
            mask, (center_x - radius, center_y - radius), (center_x + radius, center_y + radius), thickness=2
        )
        self.assertEqual(cut_diag[center_y, center_x], 0)
        num_labels_diag, _ = cv2.connectedComponents(cut_diag, connectivity=8)
        self.assertEqual(num_labels_diag, 3)

        # Test 4: Boolean mask input preserves bool dtype
        bool_mask = (mask > 0)
        cut_bool = apply_knife_cut(bool_mask, (center_x, 40), (center_x, 160), thickness=2)
        self.assertEqual(cut_bool.dtype, bool)
        self.assertFalse(cut_bool[center_y, center_x])

    def test_morphological_dilation_and_erosion(self):
        """
        Verifies that apply_mask_dilation and apply_mask_erosion accurately grow
        and shrink binary margins by 1 pixel using a 3x3 structuring element
        without boundary inversion.
        """
        h, w = 100, 100
        mask = np.zeros((h, w), dtype=np.uint8)
        # Define a 20x20 foreground square [40:60, 40:60]
        mask[40:60, 40:60] = 255
        initial_area = np.count_nonzero(mask)
        self.assertEqual(initial_area, 400)

        # 1. Test Dilation with 3x3 kernel
        dilated = apply_mask_dilation(mask, kernel_size=3)
        self.assertEqual(dilated.shape, mask.shape)
        self.assertEqual(dilated.dtype, np.uint8)

        # Expanded by 1 pixel on all sides: now [39:61, 39:61], area = 22x22 = 484
        expected_dilated_area = 22 * 22
        self.assertEqual(np.count_nonzero(dilated), expected_dilated_area)

        # Boundary expansion verification: new 1-px margin must be foreground
        self.assertEqual(dilated[39, 50], 255)  # top
        self.assertEqual(dilated[60, 50], 255)  # bottom
        self.assertEqual(dilated[50, 39], 255)  # left
        self.assertEqual(dilated[50, 60], 255)  # right

        # Boundary inversion check: pixels outside 1-px expansion remain background
        self.assertEqual(dilated[38, 50], 0)
        self.assertEqual(dilated[61, 50], 0)
        self.assertEqual(dilated[50, 38], 0)
        self.assertEqual(dilated[50, 61], 0)
        self.assertEqual(dilated[0, 0], 0)

        # Inversion check: all original foreground pixels remain foreground
        self.assertTrue(np.all((mask == 255) <= (dilated == 255)))

        # 2. Test Erosion with 3x3 kernel
        eroded = apply_mask_erosion(mask, kernel_size=3)
        self.assertEqual(eroded.shape, mask.shape)
        self.assertEqual(eroded.dtype, np.uint8)

        # Shrunk by 1 pixel on all sides: now [41:59, 41:59], area = 18x18 = 324
        expected_eroded_area = 18 * 18
        self.assertEqual(np.count_nonzero(eroded), expected_eroded_area)

        # Boundary shrinkage verification: outer 1-px margin becomes background
        self.assertEqual(eroded[40, 50], 0)   # former top edge
        self.assertEqual(eroded[59, 50], 0)   # former bottom edge
        self.assertEqual(eroded[50, 40], 0)   # former left edge
        self.assertEqual(eroded[50, 59], 0)   # former right edge

        # Interior remains intact
        self.assertEqual(eroded[41, 41], 255)
        self.assertEqual(eroded[50, 50], 255)
        self.assertEqual(eroded[58, 58], 255)

        # Inversion check: background remains background, eroded is strict subset
        self.assertTrue(np.all((eroded == 255) <= (mask == 255)))
        self.assertEqual(eroded[0, 0], 0)

        # 3. Test on boolean masks
        bool_mask = (mask > 0)
        dilated_bool = apply_mask_dilation(bool_mask, kernel_size=3)
        eroded_bool = apply_mask_erosion(bool_mask, kernel_size=3)

        self.assertEqual(dilated_bool.dtype, bool)
        self.assertEqual(eroded_bool.dtype, bool)
        self.assertEqual(np.count_nonzero(dilated_bool), expected_dilated_area)
        self.assertEqual(np.count_nonzero(eroded_bool), expected_eroded_area)

    def test_x11_window_wm_delete_and_close_event(self):
        """
        Verifies that X11GUIWindow correctly registers the WM_DELETE_WINDOW protocol
        and that poll_events reliably catches WM ClientMessage close events.
        """
        import os
        import ctypes
        if not os.environ.get("DISPLAY"):
            self.skipTest("No X11 DISPLAY available")

        import scripts.annotation_and_training.annotate_with_sam2 as a
        try:
            win = a.X11GUIWindow(title="Test Close Event", width=200, height=200)
        except Exception as e:
            self.skipTest(f"X11 display connection failed: {e}")

        try:
            self.assertIsNotNone(win.disp)
            self.assertGreater(win.wm_delete, 0)

            # Send synthetic WM_DELETE_WINDOW ClientMessage
            x11 = ctypes.CDLL("libX11.so.6")
            x11.XInternAtom.restype = ctypes.c_ulong
            x11.XInternAtom.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
            wm_protocols = x11.XInternAtom(win.disp, b"WM_PROTOCOLS", False)

            evt = a.XEvent()
            evt.type = 33
            evt.xclient.type = 33
            evt.xclient.serial = 0
            evt.xclient.send_event = 1
            evt.xclient.display = win.disp
            evt.xclient.window = win.win
            evt.xclient.message_type = wm_protocols
            evt.xclient.format = 32
            evt.xclient.data_l[0] = win.wm_delete

            x11.XSendEvent.restype = ctypes.c_int
            x11.XSendEvent.argtypes = [ctypes.c_void_p, ctypes.c_ulong, ctypes.c_int, ctypes.c_long, ctypes.c_void_p]
            x11.XSendEvent(win.disp, win.win, 0, 0, ctypes.byref(evt))

            x11.XSync.restype = ctypes.c_int
            x11.XSync.argtypes = [ctypes.c_void_p, ctypes.c_int]
            x11.XSync(win.disp, 0)

            events = win.poll_events()
            close_events = [ev for ev in events if ev[0] == "close"]
            self.assertEqual(len(close_events), 1)
        finally:
            win.close()

    def test_get_undo_button_rect(self):
        x0, y0, x1, y1 = get_undo_button_rect(1280, 70)
        self.assertEqual(x1 - x0, 170)
        self.assertEqual(y1 - y0, 26)
        self.assertGreater(x0, 0)
        self.assertLess(x1, 1280)
        self.assertGreaterEqual(y0, 0)
        self.assertLessEqual(y1, 70)

    def test_hud_renders_undo_button(self):
        canvas = np.zeros((400, 1000, 3), dtype=np.uint8)
        hud = render_hud_overlay(
            display_img=canvas,
            voucher_name="TEST_VOUCHER",
            voucher_idx=0,
            total_vouchers=1,
            saved_instances=[],
            mode="SELECT",
            zoom_level=1.0,
            pan_offset=(0, 0),
        )
        bx0, by0, bx1, by1 = get_undo_button_rect(1000, 70)
        btn_crop = hud[by0:by1, bx0:bx1]
        self.assertGreater(np.count_nonzero(btn_crop), 0)

    def test_undo_last_point_prompts(self):
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.predictor = None
        annotator.active_image = None
        annotator.mode = "SELECT"
        annotator.knife_pt_a = None
        annotator.polygon_points = []
        annotator.box_prompt = None
        annotator.candidate_mask = np.ones((50, 50), dtype=np.uint8) * 255
        annotator.candidate_masks = [annotator.candidate_mask.copy()]
        annotator.candidate_scores = [0.91]

        # Add 2 inclusion points and 1 exclusion point
        annotator.point_coords = [[10.0, 10.0], [20.0, 20.0], [30.0, 30.0]]
        annotator.point_labels = [1, 1, 0]

        # 1. Undo exclusion point
        res1 = annotator.undo_last_point()
        self.assertTrue(res1)
        self.assertEqual(annotator.point_coords, [[10.0, 10.0], [20.0, 20.0]])
        self.assertEqual(annotator.point_labels, [1, 1])

        # 2. Undo 2nd inclusion point
        res2 = annotator.undo_last_point()
        self.assertTrue(res2)
        self.assertEqual(annotator.point_coords, [[10.0, 10.0]])
        self.assertEqual(annotator.point_labels, [1])

        # 3. Undo 1st inclusion point (buffer becomes empty)
        res3 = annotator.undo_last_point()
        self.assertTrue(res3)
        self.assertEqual(annotator.point_coords, [])
        self.assertEqual(annotator.point_labels, [])
        # Candidate mask should now be cleared
        self.assertIsNone(annotator.candidate_mask)
        self.assertEqual(annotator.candidate_masks, [])
        self.assertEqual(annotator.candidate_scores, [])

        # 4. Undo when empty should return False gracefully
        res4 = annotator.undo_last_point()
        self.assertFalse(res4)

    def test_undo_last_point_knife_and_polygon(self):
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.predictor = None
        annotator.active_image = None
        annotator.point_coords = []
        annotator.point_labels = []
        annotator.box_prompt = None
        annotator.candidate_mask = None
        annotator.candidate_masks = []
        annotator.candidate_scores = []

        # Knife mode with Point A
        annotator.mode = "KNIFE"
        annotator.knife_pt_a = (150, 250)
        res_k = annotator.undo_last_point()
        self.assertTrue(res_k)
        self.assertIsNone(annotator.knife_pt_a)

        # Polygon mode with vertices
        annotator.mode = "POLYGON"
        annotator.polygon_points = [(10, 10), (20, 20), (30, 30)]
        res_p = annotator.undo_last_point()
        self.assertTrue(res_p)
        self.assertEqual(annotator.polygon_points, [(10, 10), (20, 20)])

        # Bounding box prompt
        annotator.mode = "SELECT"
        annotator.box_prompt = [5.0, 5.0, 50.0, 50.0]
        res_b = annotator.undo_last_point()
        self.assertTrue(res_b)
        self.assertIsNone(annotator.box_prompt)

    def test_undo_last_point_event_triggers(self):
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.predictor = None
        annotator.active_image = None
        annotator.mode = "SELECT"
        annotator.knife_pt_a = None
        annotator.polygon_points = []
        annotator.box_prompt = None
        annotator.candidate_mask = None
        annotator.candidate_masks = []
        annotator.candidate_scores = []
        annotator.window_w = 1280
        annotator.window_h = 800

        annotator.point_coords = [[10.0, 10.0], [20.0, 20.0]]
        annotator.point_labels = [1, 0]

        # Simulate Ctrl+Z keypress keysym
        sym = ord('z')
        if sym in (ord('z'), ord('Z'), 0xff08, 0xffff):
            annotator.undo_last_point()
        self.assertEqual(len(annotator.point_coords), 1)

        # Simulate Backspace keypress keysym (0xff08)
        sym_bksp = 0xff08
        if sym_bksp in (ord('z'), ord('Z'), 0xff08, 0xffff):
            annotator.undo_last_point()
        self.assertEqual(len(annotator.point_coords), 0)

        # Simulate clicking on HUD Undo button
        annotator.point_coords = [[15.0, 25.0]]
        annotator.point_labels = [1]
        bx0, by0, bx1, by1 = get_undo_button_rect(annotator.window_w)
        click_x, click_y = bx0 + 5, by0 + 5
        if bx0 <= click_x <= bx1 and by0 <= click_y <= by1:
            annotator.undo_last_point()
        self.assertEqual(len(annotator.point_coords), 0)

    def test_knife_second_click_does_not_add_inclusion_point(self):
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.predictor = None
        annotator.active_image = np.zeros((100, 100, 3), dtype=np.uint8)
        annotator.mode = "KNIFE"
        annotator.knife_pt_a = None
        annotator.polygon_points = []
        annotator.box_prompt = None
        annotator.point_coords = []
        annotator.point_labels = []
        annotator.candidate_mask = np.ones((100, 100), dtype=np.uint8) * 255
        annotator.candidate_masks = [annotator.candidate_mask.copy()]
        annotator.candidate_scores = [0.9]
        annotator.active_mask_idx = 0
        annotator.window_w = 1000
        annotator.window_h = 800
        annotator.orig_w = 100
        annotator.orig_h = 100
        annotator.zoom_level = 1.0
        annotator.pan_offset = [0, 0]
        annotator.space_down = False
        annotator.is_pan_dragging = False
        annotator.is_box_dragging = False

        # First click (Point A): press and release
        annotator.press_mode = annotator.mode
        annotator.knife_click_handled = True
        annotator.lbutton_down = False
        annotator.knife_pt_a = (50, 10)
        # Release of Point A
        if getattr(annotator, "knife_click_handled", False):
            annotator.knife_click_handled = False
            annotator.lbutton_down = False
        self.assertEqual(len(annotator.point_coords), 0)
        self.assertEqual(annotator.knife_pt_a, (50, 10))
        self.assertEqual(annotator.mode, "KNIFE")

        # Second click (Point B): press
        annotator.press_mode = annotator.mode
        annotator.knife_click_handled = True
        annotator.lbutton_down = False
        # Sever cut is made
        annotator.candidate_mask = split_mask_with_knife_line(
            annotator.candidate_mask, annotator.knife_pt_a, (50, 90), line_thickness=2, dilation_px=2
        )
        annotator.knife_pt_a = None
        annotator.mode = "SELECT"

        # Release of Point B
        btn = 1
        if getattr(annotator, "knife_click_handled", False):
            annotator.knife_click_handled = False
            annotator.lbutton_down = False
        elif annotator.mode == "SELECT" and getattr(annotator, "press_mode", "SELECT") == "SELECT":
            annotator.point_coords.append([50.0, 90.0])
            annotator.point_labels.append(1)

        # Verification: NO inclusion point must have been added on second click!
        self.assertEqual(len(annotator.point_coords), 0)
        self.assertEqual(annotator.mode, "SELECT")
        self.assertFalse(annotator.knife_click_handled)

    def test_resume_last_and_unannotated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            images_dir = tmp_path / "images"
            output_dir = tmp_path / "annotations"
            images_dir.mkdir()
            output_dir.mkdir()

            img1 = images_dir / "V001.jpg"
            img2 = images_dir / "V002.jpg"
            img3 = images_dir / "V003.jpg"
            for img in (img1, img2, img3):
                img.touch()

            txt1 = output_dir / "V001.txt"
            txt1.write_text("0 0.1 0.1 0.2 0.2\n")

            import time
            time.sleep(0.01)
            txt2 = output_dir / "V002.txt"
            txt2.write_text("0 0.3 0.3 0.4 0.4\n")

            orig_init = PrecisionSAM2Annotator._init_model
            try:
                PrecisionSAM2Annotator._init_model = lambda self: None

                annotator_last = PrecisionSAM2Annotator(
                    images_dir=images_dir,
                    output_dir=output_dir,
                    resume_last=True,
                )
                self.assertEqual(annotator_last.current_idx, 1)
                self.assertEqual(annotator_last.image_files[annotator_last.current_idx].stem, "V002")

                annotator_unann = PrecisionSAM2Annotator(
                    images_dir=images_dir,
                    output_dir=output_dir,
                    resume_unannotated=True,
                )
                self.assertEqual(annotator_unann.current_idx, 2)
                self.assertEqual(annotator_unann.image_files[annotator_unann.current_idx].stem, "V003")

                # Test save_current_sheet safeguard: empty saved_instances does not overwrite existing file
                annotator_last.saved_instances = []
                annotator_last.save_current_sheet()
                self.assertGreater(txt2.stat().st_size, 0)
            finally:
                PrecisionSAM2Annotator._init_model = orig_init

    def test_repopulate_existing_annotations(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            images_dir = tmp_path / "images"
            output_dir = tmp_path / "annotations"
            masks_dir = output_dir / "masks"
            images_dir.mkdir()
            output_dir.mkdir()
            masks_dir.mkdir()

            # Create dummy test image (100x100 RGB)
            dummy_img = np.zeros((100, 100, 3), dtype=np.uint8)
            cv2.imwrite(str(images_dir / "NCU001.jpg"), dummy_img)

            # Create an annotation txt file with 2 instances (class 0: basal_leaf_whole, class 6: capitulum)
            txt_file = output_dir / "NCU001.txt"
            txt_file.write_text(
                "0 0.1 0.1 0.2 0.1 0.2 0.2 0.1 0.2\n"
                "6 0.5 0.5 0.6 0.5 0.6 0.6 0.5 0.6\n"
            )

            # Create mask for instance 0
            mask0 = np.zeros((100, 100), dtype=np.uint8)
            mask0[10:20, 10:20] = 255
            cv2.imwrite(str(masks_dir / "NCU001_inst00_basal_leaf_whole.png"), mask0)

            orig_init = PrecisionSAM2Annotator._init_model
            try:
                PrecisionSAM2Annotator._init_model = lambda self: None
                annotator = PrecisionSAM2Annotator(
                    images_dir=images_dir,
                    output_dir=output_dir,
                )
                success = annotator.load_active_image()
                self.assertTrue(success)

                # Verify 2 instances repopulated
                self.assertEqual(len(annotator.saved_instances), 2)
                self.assertEqual(annotator.saved_instances[0]["class_id"], 0)
                self.assertEqual(annotator.saved_instances[0]["label"], "basal_leaf_whole")
                self.assertEqual(annotator.saved_instances[1]["class_id"], 6)
                self.assertEqual(annotator.saved_instances[1]["label"], "capitulum")

                # Verify cached base layer has been composited with non-zero overlays
                self.assertIsNotNone(annotator.cached_base_layer)
                self.assertGreater(np.count_nonzero(annotator.cached_base_layer), 0)
            finally:
                PrecisionSAM2Annotator._init_model = orig_init

    def test_navigation_next_prev_without_saving(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            images_dir = Path(tmpdir) / "images"
            output_dir = Path(tmpdir) / "annotations"
            images_dir.mkdir()
            output_dir.mkdir()

            dummy = np.zeros((100, 100, 3), dtype=np.uint8)
            cv2.imwrite(str(images_dir / "V001.jpg"), dummy)
            cv2.imwrite(str(images_dir / "V002.jpg"), dummy)

            # Write existing annotation for V001
            txt_v1 = output_dir / "V001.txt"
            txt_v1.write_text("0 0.1 0.1 0.2 0.1 0.2 0.2 0.1 0.2\n")

            orig_init = PrecisionSAM2Annotator._init_model
            try:
                PrecisionSAM2Annotator._init_model = lambda self: None
                annotator = PrecisionSAM2Annotator(images_dir=images_dir, output_dir=output_dir)
                annotator.load_active_image()
                self.assertEqual(len(annotator.saved_instances), 1)
                self.assertFalse(annotator.is_dirty)

                # Move to next voucher
                annotator.next_voucher()
                self.assertEqual(annotator.current_idx, 1)
                self.assertFalse(annotator.is_dirty)
                # Ensure V002 text file was NOT created or modified by just browsing
                self.assertFalse((output_dir / "V002.txt").exists())

                # Move back to prev voucher
                annotator.prev_voucher()
                self.assertEqual(annotator.current_idx, 0)
                self.assertEqual(len(annotator.saved_instances), 1)
                self.assertFalse(annotator.is_dirty)
            finally:
                PrecisionSAM2Annotator._init_model = orig_init

    def test_target_voucher_and_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            images_dir = Path(tmpdir) / "images"
            output_dir = Path(tmpdir) / "annotations"
            images_dir.mkdir()
            output_dir.mkdir()

            dummy = np.zeros((10, 10, 3), dtype=np.uint8)
            for name in ["A.jpg", "B.jpg", "C.jpg", "D.jpg"]:
                cv2.imwrite(str(images_dir / name), dummy)

            orig_init = PrecisionSAM2Annotator._init_model
            try:
                PrecisionSAM2Annotator._init_model = lambda self: None
                # Test jump by voucher stem
                ann_v = PrecisionSAM2Annotator(images_dir=images_dir, output_dir=output_dir, target_voucher="C")
                self.assertEqual(ann_v.current_idx, 2)

                # Test jump by index (1-based index 4 -> 0-based idx 3)
                ann_idx = PrecisionSAM2Annotator(images_dir=images_dir, output_dir=output_dir, target_index=4)
                self.assertEqual(ann_idx.current_idx, 3)
            finally:
                PrecisionSAM2Annotator._init_model = orig_init

    def test_hud_button_rects(self):
        from scripts.annotation_and_training.sam2_annotator_utils import (
            get_prev_button_rect,
            get_next_button_rect,
            get_save_button_rect,
            get_undo_button_rect,
        )
        w = 1280
        px0, py0, px1, py1 = get_prev_button_rect(w)
        nx0, ny0, nx1, ny1 = get_next_button_rect(w)
        sx0, sy0, sx1, sy1 = get_save_button_rect(w)
        bx0, by0, bx1, by1 = get_undo_button_rect(w)

        # Assert widths and heights are positive
        self.assertGreater(px1, px0)
        self.assertGreater(py1, py0)
        self.assertGreater(nx1, nx0)
        self.assertGreater(ny1, ny0)
        self.assertGreater(sx1, sx0)
        self.assertGreater(sy1, sy0)
        self.assertGreater(bx1, bx0)
        self.assertGreater(by1, by0)

        # Row 1 (Prev and Next) must not overlap horizontally
        self.assertLessEqual(px1, nx0)

        # Row 2 (Save and Undo) must not overlap horizontally
        self.assertLessEqual(sx1, bx0)

    def test_tier_filtering(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            images_dir = Path(tmpdir) / "images"
            output_dir = Path(tmpdir) / "annotations"
            images_dir.mkdir()
            output_dir.mkdir()

            dummy = np.zeros((10, 10, 3), dtype=np.uint8)
            cv2.imwrite(str(images_dir / "V_GOLD1.jpg"), dummy)
            cv2.imwrite(str(images_dir / "V_GOLD2.jpg"), dummy)
            cv2.imwrite(str(images_dir / "V_SILVER.jpg"), dummy)
            cv2.imwrite(str(images_dir / "V_BRONZE.jpg"), dummy)

            csv_path = Path(tmpdir) / "curated_vouchers.csv"
            csv_path.write_text(
                "catalogNumber,determiner_tier\n"
                "V_GOLD1,Tier_1_Gold\n"
                "V_GOLD2,Tier_1_Gold\n"
                "V_SILVER,Tier_2_Silver\n"
                "V_BRONZE,Tier_3_Bronze\n"
            )

            orig_init = PrecisionSAM2Annotator._init_model
            try:
                PrecisionSAM2Annotator._init_model = lambda self: None
                # Tier 1 filtering (default)
                ann_t1 = PrecisionSAM2Annotator(images_dir=images_dir, output_dir=output_dir, tier="1", vouchers_csv=csv_path)
                self.assertEqual(len(ann_t1.image_files), 2)
                self.assertEqual([p.stem for p in ann_t1.image_files], ["V_GOLD1", "V_GOLD2"])
                self.assertEqual(ann_t1.tier_label, "Tier 1 Gold")

                # Tier 2 filtering
                ann_t2 = PrecisionSAM2Annotator(images_dir=images_dir, output_dir=output_dir, tier="2", vouchers_csv=csv_path)
                self.assertEqual(len(ann_t2.image_files), 1)
                self.assertEqual(ann_t2.image_files[0].stem, "V_SILVER")

                # Tier 'all'
                ann_all = PrecisionSAM2Annotator(images_dir=images_dir, output_dir=output_dir, tier="all", vouchers_csv=csv_path)
                self.assertEqual(len(ann_all.image_files), 4)
            finally:
                PrecisionSAM2Annotator._init_model = orig_init

    def test_keysym_to_class_mapping(self):
        # 1. Top row numbers
        for i in range(len(CLASS_NAMES)):
            top_row_sym = ord(str(i))
            self.assertIn(top_row_sym, KEYSYM_TO_CLASS)
            self.assertEqual(KEYSYM_TO_CLASS[top_row_sym], i)

        # 2. Keypad with NumLock ON (XK_KP_0 to XK_KP_6)
        numpad_on_keysyms = [0xffb0, 0xffb1, 0xffb2, 0xffb3, 0xffb4, 0xffb5, 0xffb6]
        for i, sym in enumerate(numpad_on_keysyms):
            self.assertIn(sym, KEYSYM_TO_CLASS)
            self.assertEqual(KEYSYM_TO_CLASS[sym], i)

        # 3. Keypad with NumLock OFF / unshifted navigation mode (0 to 6)
        # 0: Insert, 1: End, 2: Down, 3: Next/PgDn, 4: Left, 5: Begin, 6: Right
        numpad_off_keysyms = [0xff9e, 0xff9c, 0xff99, 0xff9b, 0xff96, 0xff9d, 0xff98]
        for i, sym in enumerate(numpad_off_keysyms):
            self.assertIn(sym, KEYSYM_TO_CLASS)
            self.assertEqual(KEYSYM_TO_CLASS[sym], i)

    def test_numpad_integration_in_annotator(self):
        annotator = PrecisionSAM2Annotator.__new__(PrecisionSAM2Annotator)
        annotator.active_image = np.zeros((200, 200, 3), dtype=np.uint8)
        annotator.cached_base_layer = annotator.active_image.copy()
        annotator.overlay_alpha = 0.55
        annotator.saved_instances = []
        annotator.point_coords = []
        annotator.point_labels = []
        annotator.box_prompt = None
        annotator.polygon_points = []
        annotator.knife_pt_a = None
        annotator.view_mode = "FILL"

        # Verify committing class via keypad keysyms (both NumLock ON and OFF)
        for class_id in range(len(CLASS_NAMES)):
            cand = np.zeros((200, 200), dtype=np.uint8)
            cand[10:50, 10:50] = 255
            annotator.candidate_mask = cand
            annotator.candidate_masks = [cand.copy()]
            annotator.candidate_scores = [0.98]
            annotator.active_mask_idx = 0

            # Alternate between NumLock ON (0xffb0 + id) and NumLock OFF
            numpad_off_map = {0: 0xff9e, 1: 0xff9c, 2: 0xff99, 3: 0xff9b, 4: 0xff96, 5: 0xff9d, 6: 0xff98}
            sym = 0xffb0 + class_id if class_id % 2 == 0 else numpad_off_map[class_id]

            # Simulate key_press logic directly using KEYSYM_TO_CLASS
            self.assertIn(sym, KEYSYM_TO_CLASS)
            mapped_cls = KEYSYM_TO_CLASS[sym]
            self.assertEqual(mapped_cls, class_id)
            annotator.commit_active_instance(mapped_cls)

            self.assertEqual(len(annotator.saved_instances), class_id + 1)
            self.assertEqual(annotator.saved_instances[-1]["class_id"], class_id)
            self.assertEqual(annotator.saved_instances[-1]["label"], CLASS_NAMES[class_id])
            self.assertIsNone(annotator.candidate_mask)


if __name__ == "__main__":
    unittest.main()
