#!/usr/bin/env python3
"""
===============================================================================
Script: test_configure_leafmachine3.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit test suite for scripts/vision/configure_leafmachine3.py:
      - Verifies loading native `LeafMachine3/LM3_settings.yaml` via `load_native_lm3_config`
      - Verifies `resolve_lm3_models_dir`, `resolve_lm3_weights`, and `resolve_pcd_weights`
        across both LeafMachine3 modelhub stage layouts and flat `models/lm3/` directories
      - Verifies `_load_onnx_session` integration with `resolve_lm3_weights`
===============================================================================
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.configure_leafmachine3 import (
    NATIVE_LM3_SETTINGS_PATH,
    load_native_lm3_config,
    resolve_lm3_models_dir,
    resolve_lm3_weights,
    resolve_pcd_weights,
)


class TestConfigureLeafMachine3(unittest.TestCase):
    """Test suite for native LeafMachine3 settings loading and ONNX model path resolution."""

    def test_resolve_lm3_models_dir(self) -> None:
        """Verify resolve_lm3_models_dir returns a valid absolute Path to the LM3 models directory."""
        models_dir = resolve_lm3_models_dir()
        self.assertIsInstance(models_dir, Path)
        self.assertTrue(models_dir.is_absolute())

        with tempfile.TemporaryDirectory() as tmpdir:
            custom_dir = Path(tmpdir) / "custom_lm3_models"
            custom_dir.mkdir(parents=True, exist_ok=True)
            custom_resolved = resolve_lm3_models_dir(custom_dir)
            self.assertTrue(custom_resolved.is_absolute())

    def test_load_native_lm3_config(self) -> None:
        """Verify load_native_lm3_config reads the canonical LeafMachine3/LM3_settings.yaml directly."""
        self.assertTrue(NATIVE_LM3_SETTINGS_PATH.exists())
        cfg = load_native_lm3_config()
        self.assertIsInstance(cfg, dict)
        self.assertEqual(int(cfg["ingest"]["max_working_dim"]), 3200)
        self.assertEqual(int(cfg["ingest"]["jpg_quality"]), 100)
        modules = cfg["modules"]
        self.assertIs(bool(modules["leaf_segmenter"]["retina_masks"]), True)
        self.assertIs(bool(modules["specimen_segmenter"]["paperclean"]), True)
        self.assertIs(bool(cfg["report"]["masks"]["subtract_holes"]), True)
        self.assertEqual(int(modules["ect"]["num_dirs"]), 360)
        self.assertIs(bool(modules["momocs"]["include_petiole"]), False)
        self.assertIs(bool(modules["momocs"]["oriented"]), True)

    def test_resolve_lm3_weights_modelhub_and_flat_layouts(self) -> None:
        """Verify resolve_lm3_weights resolves both LeafMachine3 modelhub stage layouts and flat models/lm3 files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            stage_onnx = root / "LeafMachine3" / "models" / "leaf_segmenter" / "model.onnx"
            stage_onnx.parent.mkdir(parents=True, exist_ok=True)
            stage_onnx.write_bytes(b"ONNX_LEAF_SEG")

            unet_onnx = (
                root
                / "LeafMachine3"
                / "models"
                / "specimen_segmenter"
                / "unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"
                / "model.onnx"
            )
            unet_onnx.parent.mkdir(parents=True, exist_ok=True)
            unet_onnx.write_bytes(b"ONNX_UNETPP")

            pcd_onnx = root / "LeafMachine3" / "models" / "plant_detector" / "model.onnx"
            pcd_onnx.parent.mkdir(parents=True, exist_ok=True)
            pcd_onnx.write_bytes(b"ONNX_PCD")

            resolved_leaf = resolve_lm3_weights(
                "lm3_leaf_segmenter__yolo26x_seg_1024.onnx",
                explicit_path="models/lm3/lm3_leaf_segmenter__yolo26x_seg_1024.onnx",
                project_root=root,
            )
            self.assertEqual(resolved_leaf, stage_onnx.resolve())

            resolved_unet = resolve_lm3_weights(
                "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain.onnx",
                project_root=root,
            )
            self.assertEqual(resolved_unet, unet_onnx.resolve())

            pcd_name, pcd_resolved = resolve_pcd_weights(project_root=root)
            self.assertEqual(pcd_name, "lm3_plant_detector__yolo26x_det_1280")
            self.assertEqual(pcd_resolved, pcd_onnx.resolve())

    def test_load_onnx_session_handles_posixpath_resolution_and_dummy_weights(self) -> None:
        """Verify _load_onnx_session unpacks single-Path resolve_lm3_weights returns and stays mock-ready on dummy paths."""
        from scripts.vision.lm3_inference_engines import _load_onnx_session

        with tempfile.TemporaryDirectory() as tmpdir:
            dummy_p = Path(tmpdir) / "dummy.pth"
            dummy_p.write_bytes(b"dummy")
            w_path, m_name, provs, sess = _load_onnx_session(
                weights_path=dummy_p,
                default_model_name="lm3_leaf_segmenter__yolo26x_seg_1024",
                device="cpu",
            )
            self.assertIsInstance(w_path, Path)
            self.assertEqual(m_name, "dummy")
            self.assertEqual(provs, ["CPUExecutionProvider"])
            self.assertIsNone(sess)

            with patch(
                "scripts.vision.lm3_inference_engines.resolve_lm3_weights",
                return_value=Path(tmpdir) / "LeafMachine3" / "models" / "plant_detector" / "model.onnx",
            ):
                w_none, m_none, provs_none, sess_none = _load_onnx_session(
                    weights_path=None,
                    default_model_name="lm3_plant_detector__yolo26x_det_1280",
                    device="cpu",
                )
                self.assertEqual(m_none, "lm3_plant_detector__yolo26x_det_1280")
                self.assertTrue(str(w_none).endswith("plant_detector/model.onnx"))
                self.assertEqual(provs_none, ["CPUExecutionProvider"])
                self.assertIsNone(sess_none)


if __name__ == "__main__":
    unittest.main()
