#!/usr/bin/env python3
"""
===============================================================================
Script: configure_leafmachine3.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Consolidated, standalone configuration generator and validation engine for
    LeafMachine3 (LM3). Dynamically injects parameters from the centralized
    PipelineConfig (config/config.yaml) and synthesizes fully formed
    LM3_settings.local.yaml configuration files covering the 18-module
    LeafMachine3 ONNX pipeline for Packera morphometrics.
===============================================================================
"""

from __future__ import annotations

import argparse
import copy
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

from scripts.core.config import PROJECT_ROOT, PipelineConfig

# Canonical LeafMachine3 ONNX model identifiers and defaults
DEFAULT_MP_CONVERSION_MODEL = "lm3_mp_conversion_factor__sqrt_fit"
DEFAULT_ARCHIVAL_DETECTOR_MODEL = "lm3_archival_detector__yolo26x_det_1280"
DEFAULT_PLANT_DETECTOR_MODEL = "lm3_plant_detector__yolo26x_det_1280"
DEFAULT_SPECIMEN_SEGMENTER_MODEL = (
    "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"
)
DEFAULT_LEAF_SEGMENTER_MODEL = "lm3_leaf_segmenter__yolo26x_seg_1024"
DEFAULT_LANDMARK_DETECTOR_MODEL = "lm3_landmark_detector__yolo26x_pose_640"
DEFAULT_RULER_CLASSIFIER_ENSEMBLE: List[str] = [
    "yolo26x_cls_224",
    "yolo26n_cls_224",
    "dinov2_frozen_mlp",
]

DEFAULT_MAX_WORKING_DIM = 3200
DEFAULT_JPG_QUALITY = 95
DEFAULT_ARCHIVAL_CONF_THRES = 0.65
DEFAULT_PLANT_CONF_THRES = 0.35
DEFAULT_LEAF_SEG_CONF_THRES = 0.50
DEFAULT_ECT_NUM_DIRS = 360
DEFAULT_BATCH_SIZE = 4
DEFAULT_NUM_WORKERS = 4

LM3_18_MODULES: Tuple[str, ...] = (
    "mp_conversion_factor",
    "archival_detector",
    "plant_detector",
    "specimen_segmenter",
    "phenology_detector",
    "ruler_classifier",
    "ruler_cf",
    "leaf_segmenter",
    "morphology",
    "landmark_detector",
    "landmark_measurements",
    "leaf_orientation",
    "petiole_width",
    "bilateral_symmetry",
    "metric_grounding",
    "reporter",
    "ect",
    "momocs",
)

BASE_LM3_TEMPLATE_YAML = """
project:
  run_name: packera_lm3_run
  input_dir: data/raw_vouchers
  output_dir: LM3_Project/output
  dir_images_local: data/raw_vouchers
  dir_output: LM3_Project/output
  lm3_root: LeafMachine3
  models_dir: models/lm3
  device: cuda
  batch_size: 4
  num_workers: 4
ingest:
  max_working_dim: 3200
  jpg_quality: 95
modules:
  mp_conversion_factor:
    enabled: true
    model: "lm3_mp_conversion_factor__sqrt_fit"
  archival_detector:
    enabled: true
    model: "lm3_archival_detector__yolo26x_det_1280"
    conf_thres: 0.65
  plant_detector:
    enabled: true
    model: "lm3_plant_detector__yolo26x_det_1280"
    conf_thres: 0.35
  specimen_segmenter:
    enabled: true
    model: "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"
    paperclean: true
    LLR_Hyst_KeepBlobs: true
  phenology_detector:
    enabled: true
  ruler_classifier:
    enabled: true
    ensemble:
      - "yolo26x_cls_224"
      - "yolo26n_cls_224"
      - "dinov2_frozen_mlp"
  ruler_cf:
    enabled: true
    use_CF_predicted_by_MP: true
    fieldprism:
      enabled: true
      peer_tol: 0.03
      anchor_tol: 0.03
      allow_single_marker: true
  leaf_segmenter:
    enabled: true
    model: "lm3_leaf_segmenter__yolo26x_seg_1024"
    conf_thres: 0.50
  morphology:
    enabled: true
    rotated_box_method: "pca"
  landmark_detector:
    enabled: true
    model: "lm3_landmark_detector__yolo26x_pose_640"
  landmark_measurements:
    enabled: true
  leaf_orientation:
    enabled: true
  petiole_width:
    enabled: true
  bilateral_symmetry:
    enabled: true
  metric_grounding:
    enabled: true
  reporter:
    enabled: true
  ect:
    enabled: true
    num_dirs: 360
    palette: "CMRmap"
    apply_log_to_visual_for_bold_color: true
  momocs:
    enabled: true
    include_petiole: false
    oriented: true
"""

DEFAULT_PACKERA_LM3_CONFIG_TEMPLATE: Dict[str, Any] = yaml.safe_load(BASE_LM3_TEMPLATE_YAML)
DEFAULT_PACKERA_CONFIG_TEMPLATE: Dict[str, Any] = DEFAULT_PACKERA_LM3_CONFIG_TEMPLATE

logger = logging.getLogger("LM3_Configurator")


def get_default_lm3_template(base_yaml_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """Returns a deep copy of the canonical LeafMachine3 18-module template."""
    if base_yaml_path and Path(base_yaml_path).exists():
        with open(base_yaml_path, "r", encoding="utf-8") as f:
            loaded = yaml.safe_load(f)
        if isinstance(loaded, dict) and "modules" in loaded:
            return copy.deepcopy(loaded)
    return copy.deepcopy(DEFAULT_PACKERA_LM3_CONFIG_TEMPLATE)


def resolve_lm3_models_dir(
    pipeline_cfg: Optional[Union[PipelineConfig, Dict[str, Any]]] = None,
    project_root: Optional[Path] = None,
) -> Path:
    """
    Resolves the local directory containing LeafMachine3 ONNX models, checking
    PipelineConfig (`paths.models_dir`, `models.models_dir`, `segmentation.models_dir`),
    `LeafMachine3/models`, and `models/lm3`.
    """
    root = (Path(project_root) if project_root else PROJECT_ROOT).resolve()
    raw_candidates: List[Union[str, Path]] = []

    if isinstance(pipeline_cfg, PipelineConfig):
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


def resolve_pcd_weights(
    weights_path: Optional[Union[str, Path]] = None,
    project_root: Optional[Path] = None,
) -> Tuple[str, Path]:
    """Resolves the LeafMachine3 Plant Detector ONNX model name and path."""
    root = (Path(project_root) if project_root else PROJECT_ROOT).resolve()
    if weights_path:
        p = Path(weights_path)
        if p.is_absolute() and p.exists():
            return p.stem, p.resolve()
        if (root / p).exists():
            resolved = (root / p).resolve()
            return resolved.stem, resolved
        target_name = p.stem if p.suffix else p.name
    else:
        target_name = DEFAULT_PLANT_DETECTOR_MODEL

    models_dir = resolve_lm3_models_dir(project_root=root)
    onnx_name = target_name if target_name.endswith(".onnx") else f"{target_name}.onnx"
    candidates = [
        models_dir / onnx_name,
        models_dir / target_name / "model.onnx",
        root / "LeafMachine3" / "models" / onnx_name,
        root / "LeafMachine3" / "models" / target_name / "model.onnx",
        root / "models" / "lm3" / onnx_name,
    ]
    for c in candidates:
        if c.exists():
            return target_name, c.resolve()
    return target_name, (models_dir / onnx_name).resolve()


def _extract_cfg_attr(
    pipe_cfg: Union[PipelineConfig, Dict[str, Any]],
    section: str,
    key: str,
    default: Any,
) -> Any:
    """Safely extracts a configuration parameter from either PipelineConfig or a raw dict."""
    if isinstance(pipe_cfg, PipelineConfig):
        sec_obj = getattr(pipe_cfg, section, None)
        if sec_obj is not None and hasattr(sec_obj, key):
            val = getattr(sec_obj, key)
            if val is not None:
                return val
        raw_sec = pipe_cfg.get(section, {}) if hasattr(pipe_cfg, "get") else {}
        if isinstance(raw_sec, dict) and key in raw_sec and raw_sec[key] is not None:
            return raw_sec[key]
    elif isinstance(pipe_cfg, dict):
        raw_sec = pipe_cfg.get(section, {})
        if isinstance(raw_sec, dict) and key in raw_sec and raw_sec[key] is not None:
            return raw_sec[key]
        if key in pipe_cfg and pipe_cfg[key] is not None:
            return pipe_cfg[key]
    return default


def generate_lm3_config(
    pipeline_cfg: Optional[Union[PipelineConfig, Dict[str, Any]]] = None,
    input_dir: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    run_name: str = "packera_lm3_run",
    batch_size: int = DEFAULT_BATCH_SIZE,
    num_workers: int = DEFAULT_NUM_WORKERS,
    config_out_path: Optional[Union[str, Path]] = None,
    images_dir: Optional[Union[str, Path]] = None,
    device: Optional[str] = None,
    base_yaml_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Synthesizes and validates a LeafMachine3 configuration dictionary covering all
    18 pipeline modules from PipelineConfig (or dict), optionally writing
    `LM3_settings.local.yaml` to `config_out_path`.
    """
    pipe_cfg: Union[PipelineConfig, Dict[str, Any]] = (
        pipeline_cfg if pipeline_cfg is not None else PipelineConfig.from_yaml()
    )

    cfg = get_default_lm3_template(base_yaml_path=base_yaml_path)

    eff_input_dir = input_dir if input_dir is not None else images_dir
    if eff_input_dir is None:
        eff_input_dir = _extract_cfg_attr(
            pipe_cfg, "paths", "raw_vouchers_dir", PROJECT_ROOT / "data" / "raw_vouchers"
        )
    eff_output_dir = (
        output_dir
        if output_dir is not None
        else (PROJECT_ROOT / "LM3_Project" / "output")
    )

    in_path = Path(eff_input_dir)
    in_resolved = in_path.resolve() if in_path.is_absolute() else (PROJECT_ROOT / in_path).resolve()
    out_path = Path(eff_output_dir)
    out_resolved = out_path.resolve() if out_path.is_absolute() else (PROJECT_ROOT / out_path).resolve()

    lm3_root_raw = _extract_cfg_attr(
        pipe_cfg,
        "paths",
        "lm3_root",
        _extract_cfg_attr(pipe_cfg, "segmentation", "lm3_root", PROJECT_ROOT / "LeafMachine3"),
    )
    lm3_root_path = Path(lm3_root_raw)
    lm3_root_resolved = (
        lm3_root_path.resolve()
        if lm3_root_path.is_absolute()
        else (PROJECT_ROOT / lm3_root_path).resolve()
    )
    models_dir_resolved = resolve_lm3_models_dir(pipe_cfg)

    resolved_device = str(
        device
        or _extract_cfg_attr(pipe_cfg, "segmentation", "device", "cuda")
        or "cuda"
    )
    resolved_batch_size = max(1, int(batch_size if batch_size is not None else DEFAULT_BATCH_SIZE))
    resolved_num_workers = max(1, int(num_workers if num_workers is not None else DEFAULT_NUM_WORKERS))

    max_working_dim = int(
        _extract_cfg_attr(
            pipe_cfg,
            "segmentation",
            "max_working_dim",
            _extract_cfg_attr(pipe_cfg, "leaf_segmentation", "max_working_dim", DEFAULT_MAX_WORKING_DIM),
        )
    )
    paperclean = bool(
        _extract_cfg_attr(
            pipe_cfg,
            "segmentation",
            "paperclean",
            _extract_cfg_attr(pipe_cfg, "leaf_segmentation", "paperclean", True),
        )
    )
    llr_hyst_keep_blobs = bool(
        _extract_cfg_attr(
            pipe_cfg,
            "segmentation",
            "llr_hyst_keep_blobs",
            _extract_cfg_attr(pipe_cfg, "leaf_segmentation", "llr_hyst_keep_blobs", True),
        )
    )
    use_cf_predicted_by_mp = bool(
        _extract_cfg_attr(
            pipe_cfg,
            "segmentation",
            "use_cf_predicted_by_mp",
            _extract_cfg_attr(pipe_cfg, "leaf_segmentation", "use_cf_predicted_by_mp", True),
        )
    )
    ect_num_dirs = int(
        _extract_cfg_attr(
            pipe_cfg,
            "segmentation",
            "ect_num_dirs",
            _extract_cfg_attr(pipe_cfg, "leaf_segmentation", "ect_num_dirs", DEFAULT_ECT_NUM_DIRS),
        )
    )

    mp_model = str(
        _extract_cfg_attr(
            pipe_cfg,
            "models",
            "mp_conversion_model",
            _extract_cfg_attr(pipe_cfg, "segmentation", "mp_conversion_model", DEFAULT_MP_CONVERSION_MODEL),
        )
    )
    archival_model = str(
        _extract_cfg_attr(
            pipe_cfg,
            "models",
            "archival_detector_model",
            _extract_cfg_attr(
                pipe_cfg, "segmentation", "archival_detector_model", DEFAULT_ARCHIVAL_DETECTOR_MODEL
            ),
        )
    )
    plant_model = str(
        _extract_cfg_attr(
            pipe_cfg,
            "models",
            "plant_detector_model",
            _extract_cfg_attr(pipe_cfg, "segmentation", "plant_detector_model", DEFAULT_PLANT_DETECTOR_MODEL),
        )
    )
    specimen_seg_model = str(
        _extract_cfg_attr(
            pipe_cfg,
            "models",
            "specimen_segmenter_model",
            _extract_cfg_attr(
                pipe_cfg, "segmentation", "specimen_segmenter_model", DEFAULT_SPECIMEN_SEGMENTER_MODEL
            ),
        )
    )
    leaf_seg_model = str(
        _extract_cfg_attr(
            pipe_cfg,
            "models",
            "leaf_segmenter_model",
            _extract_cfg_attr(pipe_cfg, "segmentation", "leaf_segmenter_model", DEFAULT_LEAF_SEGMENTER_MODEL),
        )
    )
    landmark_model = str(
        _extract_cfg_attr(
            pipe_cfg,
            "models",
            "landmark_detector_model",
            _extract_cfg_attr(
                pipe_cfg, "segmentation", "landmark_detector_model", DEFAULT_LANDMARK_DETECTOR_MODEL
            ),
        )
    )

    archival_conf = float(
        _extract_cfg_attr(pipe_cfg, "segmentation", "archival_conf_thres", DEFAULT_ARCHIVAL_CONF_THRES)
    )
    plant_conf = float(
        _extract_cfg_attr(
            pipe_cfg,
            "segmentation",
            "score_thresh",
            _extract_cfg_attr(
                pipe_cfg, "plant_component_detector", "minimum_confidence_threshold", DEFAULT_PLANT_CONF_THRES
            ),
        )
    )
    leaf_seg_conf = float(
        _extract_cfg_attr(
            pipe_cfg, "segmentation", "minimum_confidence_threshold", DEFAULT_LEAF_SEG_CONF_THRES
        )
    )
    # Keep canonical LM3 defaults when config carries legacy higher threshold
    if abs(leaf_seg_conf - 0.65) < 1e-6:
        leaf_seg_conf = DEFAULT_LEAF_SEG_CONF_THRES

    # 1. Top-level & project settings
    cfg["run_name"] = str(run_name)
    cfg["input_dir"] = str(in_resolved)
    cfg["output_dir"] = str(out_resolved)
    cfg["models_dir"] = str(models_dir_resolved)
    cfg["lm3_root"] = str(lm3_root_resolved)
    cfg["device"] = resolved_device
    cfg["batch_size"] = resolved_batch_size
    cfg["num_workers"] = resolved_num_workers

    proj = cfg.setdefault("project", {})
    proj["run_name"] = str(run_name)
    proj["input_dir"] = str(in_resolved)
    proj["output_dir"] = str(out_resolved)
    proj["dir_images_local"] = str(in_resolved)
    proj["dir_output"] = str(out_resolved)
    proj["lm3_root"] = str(lm3_root_resolved)
    proj["models_dir"] = str(models_dir_resolved)
    proj["device"] = resolved_device
    proj["batch_size"] = resolved_batch_size
    proj["num_workers"] = resolved_num_workers

    # 2. Ingest settings
    ingest = cfg.setdefault("ingest", {})
    ingest["max_working_dim"] = max_working_dim
    ingest["jpg_quality"] = DEFAULT_JPG_QUALITY

    # 3. 18-module LeafMachine3 pipeline configuration
    modules = cfg.setdefault("modules", {})
    modules["mp_conversion_factor"] = {
        "enabled": True,
        "model": mp_model,
    }
    modules["archival_detector"] = {
        "enabled": True,
        "model": archival_model,
        "conf_thres": archival_conf,
    }
    modules["plant_detector"] = {
        "enabled": True,
        "model": plant_model,
        "conf_thres": plant_conf,
    }
    modules["specimen_segmenter"] = {
        "enabled": True,
        "model": specimen_seg_model,
        "paperclean": paperclean,
        "LLR_Hyst_KeepBlobs": llr_hyst_keep_blobs,
    }
    modules["phenology_detector"] = {
        "enabled": True,
    }
    modules["ruler_classifier"] = {
        "enabled": True,
        "ensemble": list(DEFAULT_RULER_CLASSIFIER_ENSEMBLE),
    }
    modules["ruler_cf"] = {
        "enabled": True,
        "use_CF_predicted_by_MP": use_cf_predicted_by_mp,
        "fieldprism": {
            "enabled": True,
            "peer_tol": 0.03,
            "anchor_tol": 0.03,
            "allow_single_marker": True,
        },
    }
    modules["leaf_segmenter"] = {
        "enabled": True,
        "model": leaf_seg_model,
        "conf_thres": leaf_seg_conf,
    }
    modules["morphology"] = {
        "enabled": True,
        "rotated_box_method": "pca",
    }
    modules["landmark_detector"] = {
        "enabled": True,
        "model": landmark_model,
    }
    modules["landmark_measurements"] = {
        "enabled": True,
    }
    modules["leaf_orientation"] = {
        "enabled": True,
    }
    modules["petiole_width"] = {
        "enabled": True,
    }
    modules["bilateral_symmetry"] = {
        "enabled": True,
    }
    modules["metric_grounding"] = {
        "enabled": True,
    }
    modules["reporter"] = {
        "enabled": True,
    }
    modules["ect"] = {
        "enabled": True,
        "num_dirs": ect_num_dirs,
        "palette": "CMRmap",
        "apply_log_to_visual_for_bold_color": True,
    }
    modules["momocs"] = {
        "enabled": True,
        "include_petiole": False,
        "oriented": True,
    }

    valid, err = validate_lm3_config(cfg)
    if not valid:
        raise ValueError(f"Invalid LeafMachine3 configuration: {err}")

    if config_out_path is not None:
        save_config_yaml(cfg, config_out_path)

    return cfg


# Alias for callers using the high-performance generator name
generate_high_performance_config = generate_lm3_config


def validate_lm3_config(config_dict: Dict[str, Any]) -> Tuple[bool, Optional[str]]:
    """Validates structural integrity and data types of a LeafMachine3 configuration dictionary."""
    if not isinstance(config_dict, dict):
        return False, "Configuration must be a dictionary"

    ingest = config_dict.get("ingest")
    if not isinstance(ingest, dict):
        return False, "Missing 'ingest' section"
    max_dim = ingest.get("max_working_dim", 0)
    if not isinstance(max_dim, int) or max_dim < 256:
        return False, f"Invalid ingest.max_working_dim ({max_dim}): must be int >= 256"
    jpg_q = ingest.get("jpg_quality", 0)
    if not isinstance(jpg_q, int) or not (1 <= jpg_q <= 100):
        return False, f"Invalid ingest.jpg_quality ({jpg_q}): must be in [1, 100]"

    modules = config_dict.get("modules")
    if not isinstance(modules, dict):
        return False, "Missing 'modules' section"

    for mod_name in LM3_18_MODULES:
        if mod_name not in modules or not isinstance(modules[mod_name], dict):
            return False, f"Missing required LM3 module configuration: 'modules.{mod_name}'"
        if "enabled" not in modules[mod_name]:
            return False, f"Module 'modules.{mod_name}' is missing 'enabled' boolean"

    for det_key in ("archival_detector", "plant_detector", "leaf_segmenter"):
        conf = float(modules[det_key].get("conf_thres", 0.0))
        if conf <= 0.0 or conf > 1.0:
            return False, f"Invalid modules.{det_key}.conf_thres ({conf}): must be in (0, 1]"

    ens = modules["ruler_classifier"].get("ensemble", [])
    if not isinstance(ens, list) or len(ens) == 0:
        return False, "modules.ruler_classifier.ensemble must be a non-empty list"

    ect_dirs = int(modules["ect"].get("num_dirs", 0))
    if ect_dirs <= 0:
        return False, f"Invalid modules.ect.num_dirs ({ect_dirs}): must be > 0"

    return True, None


def save_config_yaml(config_dict: Dict[str, Any], output_path: Union[str, Path]) -> Path:
    """Serializes configuration dictionary to a YAML file, ensuring parent dirs exist."""
    p = Path(output_path)
    p = p if p.is_absolute() else (PROJECT_ROOT / p).resolve()
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
    logger.info("Saved LM3 configuration YAML to %s", p)
    return p


def load_config_yaml(yaml_path: Union[str, Path]) -> Dict[str, Any]:
    """Loads a YAML configuration file into a Python dictionary."""
    p = Path(yaml_path)
    p = p if p.is_absolute() else (PROJECT_ROOT / p).resolve()
    if not p.exists():
        raise FileNotFoundError(f"Configuration file not found: {p}")
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def update_main_lm3_yaml(
    base_yaml_path: Union[str, Path] = "LeafMachine3/LM3_settings.local.yaml",
    project_yaml_path: Union[str, Path] = "LM3_settings.local.yaml",
    **override_kwargs: Any,
) -> Path:
    """
    Constructs LeafMachine3 configuration and writes both the project-level
    `LM3_settings.local.yaml` and optional `LeafMachine3/LM3_settings.local.yaml`.
    """
    cfg = generate_lm3_config(**override_kwargs)
    saved = save_config_yaml(cfg, project_yaml_path)

    base_p = Path(base_yaml_path)
    base_p = base_p if base_p.is_absolute() else (PROJECT_ROOT / base_p).resolve()
    if base_p.parent.exists() and base_p != saved:
        save_config_yaml(cfg, base_p)
        logger.info("Synchronized LeafMachine3 settings at %s", base_p)
    return saved


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments, with dynamic defaults from PipelineConfig."""
    pipe_cfg = PipelineConfig.from_yaml()

    parser = argparse.ArgumentParser(
        description="Programmatically configure LeafMachine3 (LM3_settings.local.yaml) for Packera morphometrics."
    )
    parser.add_argument(
        "-o",
        "--output-yaml",
        type=Path,
        default=PROJECT_ROOT / "LM3_settings.local.yaml",
        help="Destination path for the generated LeafMachine3 YAML.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=pipe_cfg.paths.raw_vouchers_dir,
        help="Directory containing input herbarium voucher images.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "LM3_Project" / "output",
        help="Destination directory for LeafMachine3 run outputs.",
    )
    parser.add_argument(
        "--run-name",
        type=str,
        default="packera_lm3_run",
        help="LeafMachine3 run identifier.",
    )
    parser.add_argument(
        "--update-main-config",
        action="store_true",
        help="Synchronize LeafMachine3/LM3_settings.local.yaml in addition to output YAML.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Batch size for LM3 ONNX inference (default: {DEFAULT_BATCH_SIZE}).",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=DEFAULT_NUM_WORKERS,
        help=f"Worker processes for LM3 pipeline stages (default: {DEFAULT_NUM_WORKERS}).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dump generated YAML to stdout without writing to disk.",
    )
    return parser.parse_args()


def main() -> None:
    """Main execution entrypoint."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    args = parse_args()

    cfg = generate_lm3_config(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        run_name=args.run_name,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
    )

    if args.dry_run:
        print(yaml.dump(cfg, default_flow_style=False, sort_keys=False))
        logger.info("Dry-run inspection complete. No files written.")
        return

    if args.update_main_config:
        update_main_lm3_yaml(
            project_yaml_path=args.output_yaml,
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            run_name=args.run_name,
            batch_size=args.batch_size,
            num_workers=args.num_workers,
        )
    else:
        save_config_yaml(cfg, args.output_yaml)


__all__ = [
    "BASE_LM3_TEMPLATE_YAML",
    "DEFAULT_PACKERA_LM3_CONFIG_TEMPLATE",
    "DEFAULT_PACKERA_CONFIG_TEMPLATE",
    "LM3_18_MODULES",
    "DEFAULT_MP_CONVERSION_MODEL",
    "DEFAULT_ARCHIVAL_DETECTOR_MODEL",
    "DEFAULT_PLANT_DETECTOR_MODEL",
    "DEFAULT_SPECIMEN_SEGMENTER_MODEL",
    "DEFAULT_LEAF_SEGMENTER_MODEL",
    "DEFAULT_LANDMARK_DETECTOR_MODEL",
    "DEFAULT_RULER_CLASSIFIER_ENSEMBLE",
    "DEFAULT_MAX_WORKING_DIM",
    "DEFAULT_JPG_QUALITY",
    "DEFAULT_ARCHIVAL_CONF_THRES",
    "DEFAULT_PLANT_CONF_THRES",
    "DEFAULT_LEAF_SEG_CONF_THRES",
    "DEFAULT_ECT_NUM_DIRS",
    "DEFAULT_BATCH_SIZE",
    "DEFAULT_NUM_WORKERS",
    "get_default_lm3_template",
    "resolve_lm3_models_dir",
    "resolve_pcd_weights",
    "generate_lm3_config",
    "generate_high_performance_config",
    "validate_lm3_config",
    "save_config_yaml",
    "load_config_yaml",
    "update_main_lm3_yaml",
    "parse_args",
    "main",
]

if __name__ == "__main__":
    main()
