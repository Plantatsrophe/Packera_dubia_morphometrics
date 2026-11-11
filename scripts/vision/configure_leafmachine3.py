#!/usr/bin/env python3
"""
===============================================================================
Script: configure_leafmachine3.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    ONNX checkpoint path resolver and native LeafMachine3 configuration loader.
    Resolves model checkpoints across `LeafMachine3/models/<stage>/` modelhub
    and `models/lm3/` directories and loads native `LeafMachine3/LM3_settings.yaml`.
===============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

import importlib.util

from scripts.core.config import PROJECT_ROOT, PipelineConfig


def _resolve_default_lm3_settings_path(project_root: Path = PROJECT_ROOT) -> Path:
    """Resolve LM3_settings.yaml across local LeafMachine3 clone, installed leafmachine3 package, and bundled config fallback."""
    candidates: List[Path] = [project_root / "LeafMachine3" / "LM3_settings.yaml"]
    try:
        spec = importlib.util.find_spec("leafmachine3")
        if spec is not None and spec.origin:
            pkg_dir = Path(spec.origin).resolve().parent
            candidates.extend([
                pkg_dir / "LM3_settings.yaml",
                pkg_dir.parent / "LM3_settings.yaml",
            ])
    except Exception:
        pass
    candidates.append(project_root / "config" / "LM3_settings.yaml")
    for cand in candidates:
        if cand.is_file():
            return cand.resolve()
    return candidates[0].resolve()


NATIVE_LM3_SETTINGS_PATH: Path = _resolve_default_lm3_settings_path()
DEFAULT_PLANT_DETECTOR_MODEL: str = "lm3_plant_detector__yolo26x_det_1280"

_LM3_MODELHUB_SUBPATHS: Dict[str, List[str]] = {
    "lm3_plant_detector__yolo26x_det_1280": [
        "plant_detector/model.onnx",
        "plant_detector/yolo26x_det_1280/model.onnx",
        "plant_detector/best.pt",
        "plant_detector/model.torchscript",
    ],
    "lm3_leaf_segmenter__yolo26x_seg_1024": [
        "leaf_segmenter/model.onnx",
        "leaf_segmenter/yolo26x_seg_1024/model.onnx",
        "leaf_segmenter/best.pt",
        "leaf_segmenter/model.torchscript",
    ],
    "lm3_landmark_detector__yolo26x_pose_640_canonical": [
        "landmark_detector/model.onnx",
        "landmark_detector/yolo26x_pose_640_canonical/model.onnx",
        "landmark_detector/best.pt",
        "landmark_detector/model.torchscript",
    ],
    "lm3_landmark_detector__yolo26x_pose_640": [
        "landmark_detector/model.onnx",
        "landmark_detector/yolo26x_pose_640_canonical/model.onnx",
        "landmark_detector/best.pt",
        "landmark_detector/model.torchscript",
    ],
    "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain": [
        "specimen_segmenter/unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain/model.onnx",
        "specimen_segmenter/model.onnx",
        "specimen_segmenter/unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain/best.pt",
        "specimen_segmenter/unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain/model.torchscript.pt",
    ],
    "lm3_archival_detector__yolo26x_det_1280": [
        "archival_detector/model.onnx",
        "archival_detector/yolo26x_det_1280/model.onnx",
        "archival_detector/best.pt",
        "archival_detector/model.torchscript",
    ],
    "lm3_ruler_classifier__dinov3_vits16_tile_four_896_mlp": [
        "ruler_classifier/model.onnx",
        "ruler_classifier/dinov3_vits16_tile_four_896_mlp/model.onnx",
        "ruler_classifier/ckpt_best.pt",
    ],
}


def load_native_lm3_config(
    settings_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Load the canonical LeafMachine3 configuration across local clone, installed package, or bundled fallback."""
    target = Path(settings_path).resolve() if settings_path is not None else _resolve_default_lm3_settings_path()
    if not target.exists():
        raise FileNotFoundError(f"Native LeafMachine3 settings file not found: {target}")
    with open(target, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        raise ValueError(f"Expected YAML mapping in {target}, got {type(data).__name__}")
    return data


def resolve_lm3_models_dir(
    pipeline_cfg: Optional[Union[PipelineConfig, Dict[str, Any], str, Path]] = None,
    project_root: Optional[Path] = None,
) -> Path:
    """
    Resolve the local directory containing LeafMachine3 ONNX models across
    PipelineConfig (`paths.models_dir`, `models.models_dir`, `segmentation.models_dir`),
    `LeafMachine3/models`, and `models/lm3`.
    """
    root = (Path(project_root) if project_root else PROJECT_ROOT).resolve()
    raw_candidates: List[Union[str, Path]] = []

    if isinstance(pipeline_cfg, (str, Path)):
        raw_candidates.append(pipeline_cfg)
    elif isinstance(pipeline_cfg, PipelineConfig):
        if getattr(pipeline_cfg.models, "models_dir", None):
            raw_candidates.append(pipeline_cfg.models.models_dir)
        if getattr(pipeline_cfg.paths, "models_dir", None):
            raw_candidates.append(pipeline_cfg.paths.models_dir)
        if getattr(pipeline_cfg.segmentation, "models_dir", None):
            raw_candidates.append(pipeline_cfg.segmentation.models_dir)
        if getattr(pipeline_cfg.paths, "lm3_root", None):
            raw_candidates.append(Path(pipeline_cfg.paths.lm3_root) / "models")
    elif isinstance(pipeline_cfg, dict):
        for sec in ("models", "paths", "segmentation"):
            sec_dict = pipeline_cfg.get(sec)
            if isinstance(sec_dict, dict):
                if sec_dict.get("models_dir"):
                    raw_candidates.append(sec_dict["models_dir"])
                if sec_dict.get("lm3_root"):
                    raw_candidates.append(Path(sec_dict["lm3_root"]) / "models")

    raw_candidates.extend([
        root / "models" / "lm3",
        root / "LeafMachine3" / "models",
        root / "models",
    ])

    resolved_candidates: List[Path] = []
    for cand in raw_candidates:
        p = Path(cand)
        p_abs = p.resolve() if p.is_absolute() else (root / p).resolve()
        if p_abs not in resolved_candidates:
            resolved_candidates.append(p_abs)

    for p_abs in resolved_candidates:
        if p_abs.is_dir() and any(p_abs.glob("**/*.onnx")):
            return p_abs
    for p_abs in resolved_candidates:
        if p_abs.exists():
            return p_abs

    return resolved_candidates[0] if resolved_candidates else (root / "models" / "lm3").resolve()


def _candidate_relpaths_for_stem(stem: str) -> List[str]:
    """Build relative candidate paths inside an LM3 models root for a given model stem."""
    rels: List[str] = []
    for ext in (".onnx", ".pt", ".torchscript", ".torchscript.pt"):
        rels.append(f"{stem}{ext}")
    rels.append(f"{stem}/model.onnx")
    for sub in _LM3_MODELHUB_SUBPATHS.get(stem, []):
        if sub not in rels:
            rels.append(sub)
    if stem.startswith("lm3_") and "__" in stem:
        stage, model_key = stem[4:].split("__", 1)
        for sub in (
            f"{stage}/model.onnx",
            f"{stage}/{model_key}/model.onnx",
            f"{stage}/best.pt",
            f"{stage}/{model_key}/best.pt",
            f"{stage}/model.torchscript",
            f"{stage}/{model_key}/model.torchscript",
            f"{stage}/{model_key}/model.torchscript.pt",
        ):
            if sub not in rels:
                rels.append(sub)
    return rels


def resolve_lm3_weights(
    default_filename: str,
    explicit_path: Optional[Union[str, Path]] = None,
    project_root: Optional[Path] = None,
) -> Path:
    """
    Resolve a LeafMachine3 ONNX/PyTorch model checkpoint path across flat `models/lm3/`
    and `LeafMachine3/models/<stage>/` modelhub directory layouts.
    """
    root = (Path(project_root) if project_root else PROJECT_ROOT).resolve()
    explicit_resolved: Optional[Path] = None
    if explicit_path is not None and str(explicit_path).strip():
        p = Path(explicit_path)
        explicit_resolved = p.resolve() if p.is_absolute() else (root / p).resolve()
        if explicit_resolved.exists():
            return explicit_resolved

    search_roots = [
        root / "models" / "lm3",
        root / "LeafMachine3" / "models",
        root / "LeafMachine3" / "leafmachine3" / "models",
        root / "models",
    ]
    stems_to_check: List[str] = []
    if explicit_resolved is not None:
        exp_stem = explicit_resolved.stem
        if exp_stem.startswith("lm3_") or exp_stem == Path(default_filename).stem:
            stems_to_check.append(exp_stem)
        else:
            return explicit_resolved

    def_stem = Path(default_filename).stem
    if def_stem not in stems_to_check:
        stems_to_check.append(def_stem)

    for stem in stems_to_check:
        rel_paths = _candidate_relpaths_for_stem(stem)
        for base_dir in search_roots:
            for rel in rel_paths:
                cand = base_dir / rel
                if cand.exists() and cand.is_file():
                    return cand.resolve()

    if explicit_resolved is not None:
        return explicit_resolved
    models_dir = resolve_lm3_models_dir(project_root=root)
    onnx_name = def_stem if def_stem.endswith(".onnx") else f"{def_stem}.onnx"
    return (models_dir / onnx_name).resolve()


def resolve_pcd_weights(
    weights_path: Optional[Union[str, Path]] = None,
    project_root: Optional[Path] = None,
) -> Tuple[str, Path]:
    """Resolve the LeafMachine3 Plant Detector ONNX model name and path."""
    root = (Path(project_root) if project_root else PROJECT_ROOT).resolve()
    if weights_path:
        p = Path(weights_path)
        target_name = p.stem if p.suffix else p.name
    else:
        target_name = DEFAULT_PLANT_DETECTOR_MODEL
    resolved = resolve_lm3_weights(
        DEFAULT_PLANT_DETECTOR_MODEL,
        explicit_path=weights_path,
        project_root=root,
    )
    return target_name, resolved


__all__ = [
    "NATIVE_LM3_SETTINGS_PATH",
    "DEFAULT_PLANT_DETECTOR_MODEL",
    "load_native_lm3_config",
    "resolve_lm3_models_dir",
    "resolve_lm3_weights",
    "resolve_pcd_weights",
]
