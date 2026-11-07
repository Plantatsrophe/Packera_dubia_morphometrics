#!/usr/bin/env python3
"""
===============================================================================
Script: test_configure_leafmachine3.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Unit test suite for scripts/vision/configure_leafmachine3.py:
      - Validates canonical LeafMachine3 configuration structure and all 18
        module configuration blocks:
          1. ingest (max_working_dim: 3200, jpg_quality: 95)
          2. mp_conversion_factor (lm3_mp_conversion_factor__sqrt_fit)
          3. archival_detector (lm3_archival_detector__yolo26x_det_1280, conf_thres: 0.65)
          4. plant_detector (lm3_plant_detector__yolo26x_det_1280, conf_thres: 0.35)
          5. specimen_segmenter (lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain, paperclean: True, LLR_Hyst_KeepBlobs: True)
          6. phenology_detector
          7. ruler_classifier (["yolo26x_cls_224", "yolo26n_cls_224", "dinov2_frozen_mlp"])
          8. ruler_cf (use_CF_predicted_by_MP: True, fieldprism)
          9. leaf_segmenter (lm3_leaf_segmenter__yolo26x_seg_1024, conf_thres: 0.50)
          10. morphology (rotated_box_method: "pca")
          11. landmark_detector (lm3_landmark_detector__yolo26x_pose_640)
          12. landmark_measurements
          13. leaf_orientation
          14. petiole_width
          15. bilateral_symmetry
          16. metric_grounding
          17. reporter
          18. ect (num_dirs: 360, palette: "CMRmap", apply_log_to_visual_for_bold_color: True)
          19. momocs (include_petiole: False, oriented: True)
      - Tests generate_lm3_config and resolve_lm3_models_dir
      - Tests dynamic PipelineConfig wiring and YAML serialization
===============================================================================
"""

from __future__ import annotations

import inspect
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict

import yaml

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig
import scripts.vision.configure_leafmachine3 as mod_lm3_cfg
from scripts.vision.configure_leafmachine3 import (
    generate_lm3_config,
    resolve_lm3_models_dir,
)


def _extract_module_block(cfg: Dict[str, Any], block_name: str) -> Dict[str, Any]:
    """Locate a named LeafMachine3 module configuration block regardless of nesting depth."""
    if block_name in cfg and isinstance(cfg[block_name], dict):
        return cfg[block_name]
    if "modules" in cfg and isinstance(cfg["modules"], dict) and block_name in cfg["modules"]:
        val = cfg["modules"][block_name]
        if isinstance(val, dict):
            return val
    if "leafmachine" in cfg and isinstance(cfg["leafmachine"], dict):
        lm = cfg["leafmachine"]
        if block_name in lm and isinstance(lm[block_name], dict):
            return lm[block_name]
        if "modules" in lm and isinstance(lm["modules"], dict) and block_name in lm["modules"]:
            val = lm["modules"][block_name]
            if isinstance(val, dict):
                return val
    # Recursive fallback search for dict key matching block_name
    for v in cfg.values():
        if isinstance(v, dict) and block_name in v and isinstance(v[block_name], dict):
            return v[block_name]
    raise KeyError(f"Module block {block_name!r} not found in LM3 configuration: {list(cfg.keys())}")


def _find_value(block: Dict[str, Any], key: str, default: Any = None) -> Any:
    """Recursively look up `key` inside `block`."""
    if key in block:
        return block[key]
    for v in block.values():
        if isinstance(v, dict):
            found = _find_value(v, key, default=None)
            if found is not None:
                return found
    return default


class TestConfigureLeafMachine3(unittest.TestCase):
    """Test suite for LeafMachine3 configuration generator and model directory resolver."""

    def test_resolve_lm3_models_dir(self) -> None:
        """Verify resolve_lm3_models_dir returns a valid absolute Path to the LM3 models directory."""
        res = resolve_lm3_models_dir()
        models_dir = Path(res[1] if isinstance(res, tuple) else res)
        self.assertIsInstance(models_dir, Path)
        self.assertTrue(models_dir.is_absolute())

        with tempfile.TemporaryDirectory() as tmpdir:
            custom_dir = Path(tmpdir) / "custom_lm3_models"
            custom_dir.mkdir(parents=True, exist_ok=True)
            sig = inspect.signature(resolve_lm3_models_dir)
            if len(sig.parameters) > 0:
                first_param = next(iter(sig.parameters.keys()))
                res_custom = resolve_lm3_models_dir(**{first_param: custom_dir})
                custom_resolved = Path(res_custom[1] if isinstance(res_custom, tuple) else res_custom)
                self.assertTrue(custom_resolved.is_absolute())

    def test_generate_lm3_config_all_18_modules_and_parameters(self) -> None:
        """Verify generate_lm3_config populates all 18 LeafMachine3 module blocks with required parameters."""
        cfg = generate_lm3_config()
        self.assertIsInstance(cfg, dict)

        # 1. ingest (max_working_dim: 3200, jpg_quality: 95)
        ingest = _extract_module_block(cfg, "ingest")
        self.assertEqual(int(_find_value(ingest, "max_working_dim")), 3200)
        self.assertEqual(int(_find_value(ingest, "jpg_quality")), 95)

        # 2. mp_conversion_factor (lm3_mp_conversion_factor__sqrt_fit)
        mp_cf = _extract_module_block(cfg, "mp_conversion_factor")
        self.assertIn("lm3_mp_conversion_factor__sqrt_fit", yaml.safe_dump(mp_cf))

        # 3. archival_detector (lm3_archival_detector__yolo26x_det_1280, conf_thres: 0.65)
        archival = _extract_module_block(cfg, "archival_detector")
        self.assertIn("lm3_archival_detector__yolo26x_det_1280", yaml.safe_dump(archival))
        self.assertAlmostEqual(float(_find_value(archival, "conf_thres")), 0.65, places=3)

        # 4. plant_detector (lm3_plant_detector__yolo26x_det_1280, conf_thres: 0.35)
        plant_det = _extract_module_block(cfg, "plant_detector")
        self.assertIn("lm3_plant_detector__yolo26x_det_1280", yaml.safe_dump(plant_det))
        self.assertAlmostEqual(float(_find_value(plant_det, "conf_thres")), 0.35, places=3)

        # 5. specimen_segmenter (lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain, paperclean: True, LLR_Hyst_KeepBlobs: True)
        spec_seg = _extract_module_block(cfg, "specimen_segmenter")
        self.assertIn(
            "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain",
            yaml.safe_dump(spec_seg),
        )
        self.assertIs(bool(_find_value(spec_seg, "paperclean")), True)
        self.assertIs(bool(_find_value(spec_seg, "LLR_Hyst_KeepBlobs")), True)

        # 6. phenology_detector
        pheno_det = _extract_module_block(cfg, "phenology_detector")
        self.assertIsInstance(pheno_det, dict)
        self.assertTrue(len(pheno_det) > 0)

        # 7. ruler_classifier (["yolo26x_cls_224", "yolo26n_cls_224", "dinov2_frozen_mlp"])
        ruler_cls = _extract_module_block(cfg, "ruler_classifier")
        ruler_cls_dump = yaml.safe_dump(ruler_cls)
        for cls_model in ("yolo26x_cls_224", "yolo26n_cls_224", "dinov2_frozen_mlp"):
            self.assertIn(cls_model, ruler_cls_dump)

        # 8. ruler_cf (use_CF_predicted_by_MP: True, fieldprism)
        ruler_cf = _extract_module_block(cfg, "ruler_cf")
        self.assertIs(bool(_find_value(ruler_cf, "use_CF_predicted_by_MP")), True)
        self.assertIn("fieldprism", yaml.safe_dump(ruler_cf).lower())

        # 9. leaf_segmenter (lm3_leaf_segmenter__yolo26x_seg_1024, conf_thres: 0.50)
        leaf_seg = _extract_module_block(cfg, "leaf_segmenter")
        self.assertIn("lm3_leaf_segmenter__yolo26x_seg_1024", yaml.safe_dump(leaf_seg))
        self.assertAlmostEqual(float(_find_value(leaf_seg, "conf_thres")), 0.50, places=3)

        # 10. morphology (rotated_box_method: "pca")
        morph = _extract_module_block(cfg, "morphology")
        self.assertEqual(str(_find_value(morph, "rotated_box_method")).lower(), "pca")

        # 11. landmark_detector (lm3_landmark_detector__yolo26x_pose_640)
        lmk_det = _extract_module_block(cfg, "landmark_detector")
        self.assertIn("lm3_landmark_detector__yolo26x_pose_640", yaml.safe_dump(lmk_det))

        # 12-17. landmark_measurements, leaf_orientation, petiole_width, bilateral_symmetry, metric_grounding, reporter
        for mod_name in (
            "landmark_measurements",
            "leaf_orientation",
            "petiole_width",
            "bilateral_symmetry",
            "metric_grounding",
            "reporter",
        ):
            mod_block = _extract_module_block(cfg, mod_name)
            self.assertIsInstance(mod_block, dict, f"Module {mod_name} must be a dict block")
            self.assertTrue(len(mod_block) > 0, f"Module {mod_name} block must be non-empty")

        # 18. ect (num_dirs: 360, palette: "CMRmap", apply_log_to_visual_for_bold_color: True)
        ect = _extract_module_block(cfg, "ect")
        self.assertEqual(int(_find_value(ect, "num_dirs")), 360)
        self.assertEqual(str(_find_value(ect, "palette")), "CMRmap")
        self.assertIs(bool(_find_value(ect, "apply_log_to_visual_for_bold_color")), True)

        # 19. momocs (include_petiole: False, oriented: True)
        momocs = _extract_module_block(cfg, "momocs")
        self.assertIs(bool(_find_value(momocs, "include_petiole")), False)
        self.assertIs(bool(_find_value(momocs, "oriented")), True)

    def test_pipeline_config_wiring_and_yaml_roundtrip(self) -> None:
        """Verify generate_lm3_config integrates with PipelineConfig and serializes cleanly to YAML."""
        pipe_cfg = PipelineConfig.from_yaml()
        sig = inspect.signature(generate_lm3_config)
        kwargs: Dict[str, Any] = {}
        if "pipeline_cfg" in sig.parameters:
            kwargs["pipeline_cfg"] = pipe_cfg
        elif "cfg" in sig.parameters:
            kwargs["cfg"] = pipe_cfg

        cfg = generate_lm3_config(**kwargs)
        self.assertIsInstance(cfg, dict)

        if hasattr(mod_lm3_cfg, "validate_lm3_config"):
            valid_res = mod_lm3_cfg.validate_lm3_config(cfg)
            if isinstance(valid_res, tuple):
                self.assertTrue(valid_res[0], f"LM3 config validation failed: {valid_res[1]}")
            else:
                self.assertTrue(bool(valid_res))

        with tempfile.TemporaryDirectory() as tmpdir:
            out_yaml = Path(tmpdir) / "test_lm3_config.yaml"
            if hasattr(mod_lm3_cfg, "save_config_yaml"):
                mod_lm3_cfg.save_config_yaml(cfg, out_yaml)
            else:
                out_yaml.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
            self.assertTrue(out_yaml.exists() and out_yaml.stat().st_size > 0)

            if hasattr(mod_lm3_cfg, "load_config_yaml"):
                loaded = mod_lm3_cfg.load_config_yaml(out_yaml)
            else:
                loaded = yaml.safe_load(out_yaml.read_text(encoding="utf-8"))

            self.assertEqual(int(_find_value(_extract_module_block(loaded, "ingest"), "max_working_dim")), 3200)
            self.assertEqual(int(_find_value(_extract_module_block(loaded, "ect"), "num_dirs")), 360)
            self.assertIs(bool(_find_value(_extract_module_block(loaded, "momocs"), "include_petiole")), False)
            self.assertIs(bool(_find_value(_extract_module_block(loaded, "momocs"), "oriented")), True)

    def test_resolve_lm3_weights_modelhub_and_flat_layouts(self) -> None:
        """Verify resolve_lm3_weights resolves both LeafMachine3 modelhub stage layouts and flat models/lm3 files."""
        from scripts.vision.configure_leafmachine3 import resolve_lm3_weights, resolve_pcd_weights

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


if __name__ == "__main__":
    unittest.main()
