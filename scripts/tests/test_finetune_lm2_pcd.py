#!/usr/bin/env python3
"""
Unit tests for LeafMachine2 Stage 2 PointRend Leaf Segmenter fine-tuning in
`scripts/annotation_and_training/finetune_lm2_pcd.py`.

Verifies:
  1. `build_pointrend_cfg()` resolves `base_weights_path` by searching for
     LeafMachine2's pre-trained `leaf_seg_rcnn_v2_1` checkpoint
     (`LeafMachine2/leafmachine2/segmentation/models/leaf_seg_rcnn_v2_1/model_final.pth`
     or `models/leaf_seg_rcnn_v2_1.pth`) before falling back to the Detectron2
     model zoo.
  2. `build_pointrend_cfg()` configures `cfg.INPUT.MIN_SIZE_TRAIN = (512, 640, 768, 896, 1024)`
     and `cfg.INPUT.MAX_SIZE_TRAIN = 1333` for single-leaf crops.
  3. `freeze_backbone_stages(model, freeze_stages=2)` freezes early ResNet backbone
     stages (`stem`, `res2`) while leaving deeper stages trainable.
  4. `save_checkpoint_atomically()` atomically exports fine-tuned weights to
     `models/lm2_packera_pointrend_finetuned.pth` and aliases
     `models/lm2_packera_pcd_finetuned.pth` for backward compatibility.
  5. End-to-end compatibility with cropped leaf COCO annotations from
     `export_coco_from_yolo_and_masks()`, verifying patch dimensions `(crop_w, crop_h)`
     and polygon bounds `[0, crop_w] x [0, crop_h]`.
"""

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import cv2
import numpy as np

from scripts.annotation_and_training.finetune_lm2_pcd import (
    DEFAULT_MAX_SIZE_TRAIN,
    DEFAULT_MIN_SIZE_TRAIN,
    DEFAULT_POINTREND_OUTPUT_WEIGHTS,
    DEFAULT_THING_CLASSES,
    LEGACY_PCD_ALIAS_WEIGHTS,
    build_pointrend_cfg,
    freeze_backbone_stages,
    resolve_pointrend_base_weights,
    save_checkpoint_atomically,
)
from scripts.annotation_and_training.sam2_annotator_utils import (
    export_coco_from_yolo_and_masks,
)


class _DummyParam:
    def __init__(self) -> None:
        self.requires_grad = True


class _DummyLayer:
    def __init__(self, n_params: int = 2) -> None:
        self._params = [_DummyParam() for _ in range(n_params)]

    def parameters(self):
        return self._params


class TestFineTuneLM2PointRend(unittest.TestCase):
    """Tests for LeafMachine2 Stage 2 PointRend Leaf Segmenter fine-tuning configuration."""

    def test_resolve_pointrend_base_weights_priority(self):
        """
        Verifies weight resolution priority:
          1. LeafMachine2/leafmachine2/segmentation/models/leaf_seg_rcnn_v2_1/model_final.pth
          2. models/leaf_seg_rcnn_v2_1.pth
          3. Detectron2 model zoo fallback when local checkpoints are absent.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)

            # 1. Neither local file exists -> falls back to Detectron2 model zoo
            fallback_w = resolve_pointrend_base_weights(
                base_weights_path=None,
                search_roots=[root],
            )
            self.assertTrue(
                "COCO-InstanceSegmentation/mask_rcnn_R_50_FPN_3x" in fallback_w
                or "detectron2" in fallback_w.lower()
            )

            # 2. Create `models/leaf_seg_rcnn_v2_1.pth` -> resolved
            models_ckpt = root / "models" / "leaf_seg_rcnn_v2_1.pth"
            models_ckpt.parent.mkdir(parents=True, exist_ok=True)
            models_ckpt.write_bytes(b"mock_pointrend_v2_1_models")

            resolved_models = resolve_pointrend_base_weights(
                base_weights_path=None,
                search_roots=[root],
            )
            self.assertEqual(Path(resolved_models), models_ckpt)

            # 3. Create `LeafMachine2/leafmachine2/segmentation/models/leaf_seg_rcnn_v2_1/model_final.pth`
            #    -> takes precedence over `models/leaf_seg_rcnn_v2_1.pth`
            lm2_ckpt = (
                root
                / "LeafMachine2"
                / "leafmachine2"
                / "segmentation"
                / "models"
                / "leaf_seg_rcnn_v2_1"
                / "model_final.pth"
            )
            lm2_ckpt.parent.mkdir(parents=True, exist_ok=True)
            lm2_ckpt.write_bytes(b"mock_pointrend_v2_1_lm2")

            resolved_lm2 = resolve_pointrend_base_weights(
                base_weights_path=None,
                search_roots=[root],
            )
            self.assertEqual(Path(resolved_lm2), lm2_ckpt)

    def test_build_pointrend_cfg_crop_scales_and_weights(self):
        """
        Verifies that `build_pointrend_cfg` sets single-leaf crop input scales
        `MIN_SIZE_TRAIN = (512, 640, 768, 896, 1024)` and `MAX_SIZE_TRAIN = 1333`
        and resolves `leaf_seg_rcnn_v2_1` weights.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            out_dir = root / "training_out"
            lm2_ckpt = (
                root
                / "LeafMachine2"
                / "leafmachine2"
                / "segmentation"
                / "models"
                / "leaf_seg_rcnn_v2_1"
                / "model_final.pth"
            )
            lm2_ckpt.parent.mkdir(parents=True, exist_ok=True)
            lm2_ckpt.write_bytes(b"leaf_seg_rcnn_v2_1_weights")

            cfg = build_pointrend_cfg(
                train_dataset_name="packera_pointrend_leaf_train",
                val_dataset_name="packera_pointrend_leaf_val",
                output_dir=out_dir,
                base_weights_path=None,
                num_classes=len(DEFAULT_THING_CLASSES),
                base_lr=0.0001,
                max_iters=800,
                batch_size=2,
                device="cpu",
                search_roots=[root],
            )

            self.assertIsNotNone(cfg)
            self.assertEqual(tuple(cfg.INPUT.MIN_SIZE_TRAIN), (512, 640, 768, 896, 1024))
            self.assertEqual(int(cfg.INPUT.MAX_SIZE_TRAIN), 1333)
            self.assertEqual(Path(cfg.MODEL.WEIGHTS), lm2_ckpt)
            self.assertEqual(cfg.MODEL.ROI_HEADS.NUM_CLASSES, 2)
            self.assertEqual(cfg.MODEL.POINT_HEAD.NUM_CLASSES, 2)
            self.assertEqual(cfg.DATASETS.TRAIN, ("packera_pointrend_leaf_train",))
            self.assertEqual(cfg.DATASETS.TEST, ("packera_pointrend_leaf_val",))

    def test_freeze_backbone_stages_preserves_later_stages(self):
        """
        Verifies that `freeze_backbone_stages(model, freeze_stages=2)` freezes
        `stem` and `res2` while keeping `res3` and `res4` trainable.
        """
        bottom_up = SimpleNamespace(
            stem=_DummyLayer(2),
            res2=_DummyLayer(3),
            res3=_DummyLayer(3),
            res4=_DummyLayer(3),
        )
        model = SimpleNamespace(backbone=SimpleNamespace(bottom_up=bottom_up))

        freeze_backbone_stages(model, freeze_stages=2)

        for p in bottom_up.stem.parameters():
            self.assertFalse(p.requires_grad)
        for p in bottom_up.res2.parameters():
            self.assertFalse(p.requires_grad)
        for p in bottom_up.res3.parameters():
            self.assertTrue(p.requires_grad)
        for p in bottom_up.res4.parameters():
            self.assertTrue(p.requires_grad)

    def test_atomic_checkpoint_saving_and_backward_compatible_alias(self):
        """
        Verifies atomic checkpoint export to `models/lm2_packera_pointrend_finetuned.pth`
        and backward-compatible alias `models/lm2_packera_pcd_finetuned.pth`.
        """
        self.assertEqual(
            DEFAULT_POINTREND_OUTPUT_WEIGHTS,
            Path("models/lm2_packera_pointrend_finetuned.pth"),
        )
        self.assertEqual(
            LEGACY_PCD_ALIAS_WEIGHTS,
            Path("models/lm2_packera_pcd_finetuned.pth"),
        )

        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            src_ckpt = root / "output" / "model_final.pth"
            src_ckpt.parent.mkdir(parents=True, exist_ok=True)
            payload = b"PACKERA_STAGE2_POINTREND_CHECKPOINT_BYTES"
            src_ckpt.write_bytes(payload)

            primary_out = root / "models" / "lm2_packera_pointrend_finetuned.pth"
            alias_out = root / "models" / "lm2_packera_pcd_finetuned.pth"

            saved_primary, saved_alias = save_checkpoint_atomically(
                src_checkpoint=src_ckpt,
                output_weights=primary_out,
                alias_weights=alias_out,
            )

            self.assertEqual(saved_primary, primary_out)
            self.assertEqual(saved_alias, alias_out)
            self.assertTrue(primary_out.exists())
            self.assertTrue(alias_out.exists())
            self.assertEqual(primary_out.read_bytes(), payload)
            self.assertEqual(alias_out.read_bytes(), payload)

    def test_cropped_patch_coco_annotations_match_pointrend_training_contract(self):
        """
        End-to-end contract test: verifies that COCO annotations exported by
        `export_coco_from_yolo_and_masks()` match cropped patch dimensions `(crop_w, crop_h)`
        and have polygon coordinates strictly inside `[0, crop_w] x [0, crop_h]`.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            root = Path(tmp_dir)
            images_dir = root / "raw_vouchers"
            masks_dir = root / "raw_annotations" / "masks"
            crops_dir = root / "annotations" / "lm2_leaf_crops"
            coco_json = root / "annotations" / "annotations_packera_train.json"
            images_dir.mkdir(parents=True)
            masks_dir.mkdir(parents=True)

            sheet = np.full((600, 500, 3), 180, dtype=np.uint8)
            cv2.imwrite(str(images_dir / "SPECIMEN_01.jpg"), sheet)

            mask = np.zeros((600, 500), dtype=np.uint8)
            cv2.ellipse(mask, (250, 300), (60, 90), 15, 0, 360, 255, -1)
            cv2.imwrite(str(masks_dir / "SPECIMEN_01_inst00_basal_leaf_whole.png"), mask)

            coco_data = export_coco_from_yolo_and_masks(
                masks_dir=masks_dir,
                output_coco_path=coco_json,
                images_dir=images_dir,
                crops_dir=crops_dir,
                min_area_px=20.0,
            )

            self.assertEqual(len(coco_data["images"]), 1)
            self.assertEqual(len(coco_data["annotations"]), 1)

            img_entry = coco_data["images"][0]
            ann_entry = coco_data["annotations"][0]
            crop_w, crop_h = img_entry["width"], img_entry["height"]

            crop_file = crops_dir / img_entry["file_name"]
            self.assertTrue(crop_file.exists())
            crop_img = cv2.imread(str(crop_file))
            self.assertEqual((crop_img.shape[1], crop_img.shape[0]), (crop_w, crop_h))

            for poly in ann_entry["segmentation"]:
                for x in poly[0::2]:
                    self.assertGreaterEqual(x, 0.0)
                    self.assertLessEqual(x, float(crop_w))
                for y in poly[1::2]:
                    self.assertGreaterEqual(y, 0.0)
                    self.assertLessEqual(y, float(crop_h))


if __name__ == "__main__":
    unittest.main()
