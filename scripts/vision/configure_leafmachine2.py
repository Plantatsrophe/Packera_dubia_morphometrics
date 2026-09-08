#!/usr/bin/env python3
"""
===============================================================================
Script: configure_leafmachine2.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Consolidated, standalone configuration generator and validation engine for
    LeafMachine2 (LM2). Dynamically injects high-performance parameters from
    the centralized PipelineConfig (config/config.yaml) and synthesizes
    fully formed YAML configuration files for Packera morphometrics.
===============================================================================
"""

from __future__ import annotations

import argparse
import copy
import logging
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Union

import yaml

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

from scripts.core.config import PROJECT_ROOT, PipelineConfig

# Default Packera hyperparameter constants calibrated for dense rosettes
DEFAULT_PCD_WEIGHTS_NAME = "lm2_packera_pcd_finetuned.pth"
DEFAULT_SEG_MODEL_NAME = "Packera_LeafPriority"
DEFAULT_PCD_CONFIDENCE = 0.72
DEFAULT_NMS_THRESH = 0.75
DEFAULT_SUBDIVISION_STEPS = 5
DEFAULT_MIN_LEAF_AREA = 500
DEFAULT_SEG_CONFIDENCE = 0.70
DEFAULT_EFD_ORDER = 40
DEFAULT_BATCH_SIZE = 50
DEFAULT_NUM_WORKERS = 8

# Canonical LeafMachine2 base template matching full schema
BASE_LM2_TEMPLATE_YAML = """
leafmachine:
  archival_component_detector:
    detector_iteration: PREP_final
    detector_type: Archival_Detector
    detector_version: PREP_final
    detector_weights: best.pt
    do_save_prediction_overlay_images: true
    ignore_objects_for_overlay: []
    minimum_confidence_threshold: 0.7
    ruler_binary_detector: model_scripted_resnet_720_withCompression.pt
    ruler_detector: ruler_classifier_38classes_v-1.pt
  armature_component_detector:
    detector_iteration: ARM_A_1000
    detector_type: Armature_Detector
    detector_version: ARM_A_1000
    detector_weights: best.pt
    do_save_prediction_overlay_images: true
    ignore_objects_for_overlay: []
    minimum_confidence_threshold: 0.5
  cropped_components:
    binarize_labels: false
    binarize_labels_skeletonize: false
    do_save_cropped_annotations: false
    save_cropped_annotations: [label]
    save_per_annotation_class: false
    save_per_image: false
  data:
    do_apply_conversion_factor: true
    include_darwin_core_data_from_combined_file: false
    save_individual_csv_files_landmarks: true
    save_individual_csv_files_measurements: true
    save_individual_csv_files_rulers: true
    save_individual_efd_files: false
    save_json_measurements: false
    save_json_rulers: false
  do:
    check_for_corrupt_images_make_vertical: false
    check_for_illegal_filenames: false
    run_leaf_processing: true
  landmark_detector:
    detector_iteration: Landmarks
    detector_type: Landmark_Detector_YOLO
    detector_version: Landmarks
    detector_weights: best.pt
    do_save_QC_images: false
    do_save_final_images: false
    do_save_prediction_overlay_images: false
    do_show_QC_images: false
    do_show_final_images: false
    ignore_objects_for_overlay: []
    landmark_partial_leaves: false
    landmark_whole_leaves: true
    minimum_confidence_threshold: 0.02
    use_existing_landmark_detections: null
  landmark_detector_armature:
    detector_iteration: Landmarks_Arm_A_200
    detector_type: Landmark_Detector_YOLO
    detector_version: Landmarks_Arm_A_200
    detector_weights: last.pt
    do_save_QC_images: false
    do_save_final_images: false
    do_save_prediction_overlay_images: false
    do_show_QC_images: false
    do_show_final_images: false
    ignore_objects_for_overlay: []
    minimum_confidence_threshold: 0.06
    upscale_factor: 10
    use_existing_landmark_detections: null
  leaf_segmentation:
    calculate_elliptic_fourier_descriptors: true
    detector_version: uniform_spaced_oriented_traces_mid15_pet5_clean_640_flipidx_pt2
    elliptic_fourier_descriptor_order: 40
    find_minimum_bounding_box: true
    generate_overlay: true
    keep_only_best_one_leaf_one_petiole: true
    minimum_confidence_threshold: 0.7
    NMS_thresh: 0.75
    min_leaf_area: 500
    pointrend_subdivision_steps: 5
    overlay_background_color: black
    overlay_dpi: 300
    overlay_line_width: 1
    save_each_segmentation_overlay_image: true
    save_full_image_masks_color: true
    save_individual_overlay_images: true
    save_keypoint_overlay: true
    save_masks_color: true
    save_oriented_images: true
    save_oriented_mask: true
    save_rgb_cropped_images: true
    save_segmentation_overlay_images_to_pdf: false
    save_simple_txt: true
    segment_partial_leaves: false
    segment_whole_leaves: true
    segmentation_model: Packera_LeafPriority
    use_efds_for_png_masks: false
  logging:
    log_level: null
  modules:
    armature: false
    specimen_crop: false
  overlay:
    alpha_transparency_archival: 0.3
    alpha_transparency_plant: 0.0
    alpha_transparency_seg_partial_leaf: 0.3
    alpha_transparency_seg_whole_leaf: 0.4
    ignore_archival_detections_classes: []
    ignore_landmark_classes: []
    ignore_plant_detections_classes: [leaf_whole]
    line_width_archival: 12
    line_width_efd: 12
    line_width_plant: 12
    line_width_seg: 12
    overlay_background_color: black
    overlay_dpi: 300
    save_overlay_to_jpgs: true
    save_overlay_to_pdf: false
    show_archival_detections: true
    show_landmarks: true
    show_plant_detections: true
    show_segmentations: true
  plant_component_detector:
    detector_iteration: PLANT_LeafPriority
    detector_type: Plant_Detector
    detector_version: PLANT_LeafPriority
    detector_weights: LeafPriority.pt
    do_save_prediction_overlay_images: true
    ignore_objects_for_overlay: []
    minimum_confidence_threshold: 0.72
    PCD_confidence: 0.72
  print:
    optional_warnings: true
    verbose: true
  project:
    device: cuda
    GBIF_mode: all
    accept_only_ideal_leaves: true
    treat_leaflet_as_leaf_whole: false
    batch_size: 50
    censor_archival_components: false
    dir_images_local: LM2_Project/Data/images
    dir_images_subset: ''
    dir_output: LM2_Project/Data/output
    hide_archival_components: [ruler, barcode, label, colorcard, map, photo, weights]
    image_location: local
    minimum_total_reproductive_counts: 0
    n_images_per_species: 1000000
    num_workers: 1
    num_workers_ruler: 8
    num_workers_seg: 8
    overlay_background_color: black
    overlay_dpi: 300
    path_combined_csv_local: null
    path_images_csv_local: null
    path_occurrence_csv_local: null
    process_subset_of_images: false
    replacement_color: '#FFFFFF'
    run_name: Packera_dubia_LM2
    species_list: ''
    use_CF_predictor: true
    auto_cache_annotations: true
    regenerate_annotations: false
    use_existing_archival_component_detections: LM2_Project/Data/precomputed_annotations/Archival_Components/labels
    use_existing_plant_component_detections: LM2_Project/Data/precomputed_annotations/Plant_Components/labels
  ruler_detection:
    detect_ruler_type: true
    minimum_confidence_threshold: 0.5
    ruler_binary_detector: model_scripted_resnet_720_withCompression.pt
    ruler_detector: ruler_classifier_38classes_v-1.pt
    save_ruler_processed: false
    save_ruler_validation: false
    save_ruler_validation_summary: true
"""

DEFAULT_PACKERA_CONFIG_TEMPLATE: Dict[str, Any] = yaml.safe_load(BASE_LM2_TEMPLATE_YAML)

logger = logging.getLogger("LM2_Configurator")


def get_default_lm2_template(base_yaml_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """Returns a copy of the canonical LeafMachine2 template, optionally from a base YAML."""
    if base_yaml_path and Path(base_yaml_path).exists():
        with open(base_yaml_path, "r", encoding="utf-8") as f:
            loaded = yaml.safe_load(f)
        if isinstance(loaded, dict) and "leafmachine" in loaded:
            return copy.deepcopy(loaded)
    return copy.deepcopy(DEFAULT_PACKERA_CONFIG_TEMPLATE)


def resolve_pcd_weights(
    weights_path: Optional[Union[str, Path]] = None,
    project_root: Optional[Path] = None,
) -> Tuple[str, Path]:
    """Resolves the weights name and path to the Plant Component Detector weights file."""
    root = project_root or PROJECT_ROOT
    target_name = Path(weights_path).name if weights_path else DEFAULT_PCD_WEIGHTS_NAME

    if weights_path:
        p = Path(weights_path)
        if p.is_absolute() and p.exists():
            return p.name, p
        if (root / p).exists():
            return p.name, (root / p).resolve()

    candidates = [
        root / "models" / target_name,
        root / "models" / "lm2_packera_pcd_finetuned.pth",
        root / "models" / "Packera_LeafPriority.pth",
        root / "models" / "LeafPriority.pth",
        root / "models" / "LeafPriority.pt",
        root / "LeafMachine2" / "leafmachine2" / "component_detector" / "models" / target_name,
    ]
    for c in candidates:
        if c.exists():
            return c.name, c.resolve()

    return target_name, (root / "models" / target_name).resolve()


def generate_high_performance_config(
    pipeline_cfg: Optional[PipelineConfig] = None,
    images_dir: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    run_name: str = "Packera_dubia_LM2",
    pcd_weights: Optional[Union[str, Path]] = None,
    pcd_confidence: Optional[float] = None,
    nms_thresh: float = DEFAULT_NMS_THRESH,
    pointrend_subdivision_steps: int = DEFAULT_SUBDIVISION_STEPS,
    min_leaf_area: int = DEFAULT_MIN_LEAF_AREA,
    seg_confidence: float = DEFAULT_SEG_CONFIDENCE,
    batch_size: Optional[int] = None,
    num_workers: Optional[int] = None,
    device: Optional[str] = None,
    use_pointrend: bool = True,
    base_yaml_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Generates a high-performance LeafMachine2 configuration dictionary, dynamically
    wiring parameters from PipelineConfig (config/config.yaml).
    """
    pipe_cfg = pipeline_cfg or PipelineConfig.from_yaml()

    # Dynamic defaults from PipelineConfig
    resolved_pcd_conf = pcd_confidence if pcd_confidence is not None else float(pipe_cfg.thresholds.pcd_conf)
    resolved_device = device or pipe_cfg.segmentation.device or "cuda"
    resolved_batch_size = batch_size if batch_size is not None else DEFAULT_BATCH_SIZE
    resolved_num_workers = num_workers if num_workers is not None else DEFAULT_NUM_WORKERS

    # Model weights resolution
    target_weights = pcd_weights or pipe_cfg.paths.model_weights or pipe_cfg.segmentation.model_weights
    weights_name, weights_path = resolve_pcd_weights(target_weights)

    # Base template
    cfg = get_default_lm2_template(base_yaml_path=base_yaml_path)
    lm = cfg.setdefault("leafmachine", {})

    # 1. Project-level configuration
    root = PROJECT_ROOT
    staged_images = Path(images_dir) if images_dir else root / "LM2_Project" / "Data" / "images"
    out_dir = Path(output_dir) if output_dir else root / "LM2_Project" / "Data" / "output"
    archival_det = root / "LM2_Project" / "Data" / "precomputed_annotations" / "Archival_Components" / "labels"
    plant_det = root / "LM2_Project" / "Data" / "precomputed_annotations" / "Plant_Components" / "labels"

    proj = lm.setdefault("project", {})
    proj["dir_images_local"] = str(staged_images.resolve())
    proj["dir_output"] = str(out_dir.resolve())
    proj["run_name"] = run_name
    proj["device"] = resolved_device
    proj["batch_size"] = resolved_batch_size
    proj["num_workers_ruler"] = resolved_num_workers
    proj["num_workers_seg"] = resolved_num_workers
    proj["use_segmentation"] = True
    proj["use_existing_archival_component_detections"] = str(archival_det.resolve())
    proj["use_existing_plant_component_detections"] = str(plant_det.resolve())

    # 2. Plant Component Detector (PCD)
    pcd = lm.setdefault("plant_component_detector", {})
    pcd["detector_weights"] = weights_name
    pcd["minimum_confidence_threshold"] = resolved_pcd_conf
    pcd["PCD_confidence"] = resolved_pcd_conf
    pcd["iou_threshold"] = nms_thresh
    pcd["batch_size"] = resolved_batch_size
    pcd["num_workers"] = resolved_num_workers
    pcd["device"] = resolved_device

    # 3. Leaf Segmentation (PointRend)
    seg = lm.setdefault("leaf_segmentation", {})
    seg["segmentation_model"] = DEFAULT_SEG_MODEL_NAME
    seg["segmentation_type"] = "Detectron2_PointRend"
    seg["use_pointrend"] = use_pointrend
    seg["minimum_confidence_threshold"] = seg_confidence
    seg["NMS_thresh"] = nms_thresh
    seg["pointrend_subdivision_steps"] = pointrend_subdivision_steps
    seg["min_leaf_area"] = min_leaf_area
    seg["minimum_leaf_area_px"] = min_leaf_area
    seg["calculate_elliptic_fourier_descriptors"] = True
    seg["elliptic_fourier_descriptor_order"] = DEFAULT_EFD_ORDER
    seg["efd_order"] = DEFAULT_EFD_ORDER
    seg["device"] = resolved_device
    seg["num_workers_seg"] = resolved_num_workers

    # 4. Ruler Detection
    ruler = lm.setdefault("ruler_detection", {})
    ruler["num_workers_ruler"] = resolved_num_workers

    return cfg


def validate_lm2_config(config_dict: Dict[str, Any]) -> Tuple[bool, Optional[str]]:
    """Validates structural integrity and data types of an LM2 configuration dictionary."""
    if not isinstance(config_dict, dict) or "leafmachine" not in config_dict:
        return False, "Missing top-level 'leafmachine' key"

    lm = config_dict["leafmachine"]
    if "project" not in lm:
        return False, "Missing 'project' subsection"

    pcd = lm.get("plant_component_detector", {})
    pcd_conf = pcd.get("minimum_confidence_threshold", 0.0)
    if pcd_conf <= 0.0 or pcd_conf > 1.0:
        return False, f"Invalid PCD confidence threshold ({pcd_conf}): must be in (0, 1]"

    seg = lm.get("leaf_segmentation", {})
    seg_conf = seg.get("minimum_confidence_threshold", 0.0)
    if seg_conf <= 0.0 or seg_conf > 1.0:
        return False, f"Invalid segmentation confidence threshold ({seg_conf}): must be in (0, 1]"

    return True, None


def save_config_yaml(config_dict: Dict[str, Any], output_path: Union[str, Path]) -> Path:
    """Serializes configuration dictionary to a YAML file, ensuring parent dirs exist."""
    p = Path(output_path)
    p = p if p.is_absolute() else (PROJECT_ROOT / p).resolve()
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
    logger.info("Saved LM2 configuration YAML to %s", p)
    return p


def load_config_yaml(yaml_path: Union[str, Path]) -> Dict[str, Any]:
    """Loads a YAML configuration file into a Python dictionary."""
    p = Path(yaml_path)
    p = p if p.is_absolute() else (PROJECT_ROOT / p).resolve()
    if not p.exists():
        raise FileNotFoundError(f"Configuration file not found: {p}")
    with open(p, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def update_main_lm2_yaml(
    base_yaml_path: Union[str, Path] = "LeafMachine2/LeafMachine2.yaml",
    project_yaml_path: Union[str, Path] = "LM2_Project/configs/lm2_packera_highperf.yaml",
    **override_kwargs,
) -> None:
    """
    Constructs high-performance configuration and synchronizes both base LeafMachine2.yaml
    and the project configuration file.
    """
    cfg = generate_high_performance_config(**override_kwargs)
    valid, err = validate_lm2_config(cfg)
    if not valid:
        raise ValueError(f"Generated configuration is invalid: {err}")

    save_config_yaml(cfg, project_yaml_path)

    base_p = Path(base_yaml_path)
    base_p = base_p if base_p.is_absolute() else (PROJECT_ROOT / base_p).resolve()
    if base_p.parent.exists():
        save_config_yaml(cfg, base_p)
        logger.info("Synchronized main base config at %s", base_p)


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments, with dynamic defaults from PipelineConfig."""
    pipe_cfg = PipelineConfig.from_yaml()

    parser = argparse.ArgumentParser(
        description="Programmatically configure LeafMachine2 for Packera dubia morphometrics."
    )
    parser.add_argument(
        "-o",
        "--output-yaml",
        type=Path,
        default=PROJECT_ROOT / "LM2_Project" / "configs" / "lm2_packera_highperf.yaml",
        help="Destination path for the generated LeafMachine2 YAML.",
    )
    parser.add_argument(
        "--update-main-config",
        action="store_true",
        help="Synchronize LeafMachine2/LeafMachine2.yaml in addition to output YAML.",
    )
    parser.add_argument(
        "--pcd-weights",
        type=str,
        default=None,
        help="Path or filename of PCD weights (defaults to models/lm2_packera_pcd_finetuned.pth).",
    )
    parser.add_argument(
        "--pcd-confidence",
        type=float,
        default=float(pipe_cfg.thresholds.pcd_conf),
        help=f"PCD detection confidence threshold (default: {pipe_cfg.thresholds.pcd_conf}).",
    )
    parser.add_argument(
        "--nms-thresh",
        type=float,
        default=DEFAULT_NMS_THRESH,
        help=f"IoU / NMS threshold (default: {DEFAULT_NMS_THRESH}).",
    )
    parser.add_argument(
        "--pointrend-subdivision-steps",
        type=int,
        default=DEFAULT_SUBDIVISION_STEPS,
        help=f"PointRend subdivision steps (default: {DEFAULT_SUBDIVISION_STEPS}).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Batch size for PCD / LM2 processing (default: {DEFAULT_BATCH_SIZE}).",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=DEFAULT_NUM_WORKERS,
        help=f"Worker processes for ruler/segmentation (default: {DEFAULT_NUM_WORKERS}).",
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

    cfg = generate_high_performance_config(
        pcd_weights=args.pcd_weights,
        pcd_confidence=args.pcd_confidence,
        nms_thresh=args.nms_thresh,
        pointrend_subdivision_steps=args.pointrend_subdivision_steps,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
    )

    valid, err = validate_lm2_config(cfg)
    if not valid:
        logger.error("Configuration validation failed: %s", err)
        sys.exit(1)

    if args.dry_run:
        print(yaml.dump(cfg, default_flow_style=False, sort_keys=False))
        logger.info("Dry-run inspection complete. No files written.")
        return

    if args.update_main_config:
        update_main_lm2_yaml(project_yaml_path=args.output_yaml)
    else:
        save_config_yaml(cfg, args.output_yaml)


__all__ = [
    "BASE_LM2_TEMPLATE_YAML",
    "DEFAULT_PACKERA_CONFIG_TEMPLATE",
    "DEFAULT_PCD_WEIGHTS_NAME",
    "DEFAULT_SEG_MODEL_NAME",
    "DEFAULT_PCD_CONFIDENCE",
    "DEFAULT_NMS_THRESH",
    "DEFAULT_SUBDIVISION_STEPS",
    "DEFAULT_MIN_LEAF_AREA",
    "DEFAULT_SEG_CONFIDENCE",
    "DEFAULT_EFD_ORDER",
    "DEFAULT_BATCH_SIZE",
    "DEFAULT_NUM_WORKERS",
    "get_default_lm2_template",
    "resolve_pcd_weights",
    "generate_high_performance_config",
    "validate_lm2_config",
    "save_config_yaml",
    "load_config_yaml",
    "update_main_lm2_yaml",
    "parse_args",
    "main",
]

if __name__ == "__main__":
    main()
