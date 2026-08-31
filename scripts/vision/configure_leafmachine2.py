#!/usr/bin/env python3
"""
===============================================================================
Script: configure_leafmachine2.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Programmatically creates, validates, and updates LeafMachine2 (LM2)
    configuration files and Python dictionaries optimized for high-performance
    compute hardware (48GB+ RAM, 8GB+ VRAM GPU) and tailored for the dense,
    overlapping, tomentose basal rosettes of the Packera dubia complex.

    Key Specialized Defaults Configured:
      - Plant Component Detector (PCD) Weights: 'Packera_LeafPriority.pth'
          Custom fine-tuned PointRend checkpoint trained on curated Packera vouchers.
      - Component Detection Thresholds (Overlapping Rosettes):
          PCD minimum confidence threshold: 0.20 (captures occluded/partial leaves)
          PCD_confidence: 0.20
      - Non-Maximum Suppression (NMS) Overlap Threshold: 0.75
          Permits heavy bounding box overlap without prematurely suppressing
          densely clustered basal leaves in rosette formations.
      - PointRend Subdivision Steps: 5 (Maximum Resolution)
          Forces PointRend to execute recursive 5-step point sampling,
          guaranteeing high-fidelity boundary delineation along complex,
          tomentose, crenate, and lyrately-lobed leaf margins.
      - Minimum Leaf Area: 500 px (retains fragmented rosette blades)
      - Batch size: 50 (optimal throughput for 8GB+ VRAM without CUDA OOM)
      - Parallel workers: 8 (num_workers_ruler: 8, num_workers_seg: 8)
      - Compute device: 'cuda' (GPU acceleration)
      - Integrated paths:
          dir_images_local: LM2_Project/Data/images
          dir_output: LM2_Project/Data/output
          run_name: Packera_dubia_LM2

Usage:
    # 1. Update both LeafMachine2.yaml and LM2_Project/configs/lm2_packera_highperf.yaml:
    python scripts/vision/configure_leafmachine2.py --update-main-config

    # 2. Generate customized YAML configuration:
    python scripts/vision/configure_leafmachine2.py \\
        -o LM2_Project/configs/custom_packera.yaml \\
        --pcd-weights Packera_LeafPriority.pth \\
        --pcd-confidence 0.20 \\
        --nms-thresh 0.75 \\
        --pointrend-subdivision-steps 5

    # 3. Dry-run inspection without writing files:
    python scripts/vision/configure_leafmachine2.py --dry-run

    # 4. Programmatic Python API:
    from scripts.vision.configure_leafmachine2 import generate_high_performance_config, save_config_yaml
    cfg = generate_high_performance_config(pcd_confidence=0.20, nms_thresh=0.75, pointrend_subdivision_steps=5)
    save_config_yaml(cfg, "LM2_Project/configs/lm2_packera_highperf.yaml")
===============================================================================
"""

import argparse
import copy
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import yaml

# Determine project root and LeafMachine2 root paths
PROJECT_ROOT = Path(__file__).resolve().parents[2]
LM2_ROOT = PROJECT_ROOT / "LeafMachine2"
LM2_PROJECT_DIR = PROJECT_ROOT / "LM2_Project"

DEFAULT_BASE_YAML = LM2_ROOT / "LeafMachine2.yaml"
DEFAULT_OUTPUT_YAML = LM2_PROJECT_DIR / "configs" / "lm2_packera_highperf.yaml"
DEFAULT_IMAGES_DIR = LM2_PROJECT_DIR / "Data" / "images"
DEFAULT_OUTPUT_DIR = LM2_PROJECT_DIR / "Data" / "output"

# Packera-specific defaults for dense rosettes & PointRend
DEFAULT_PCD_WEIGHTS_NAME = "LeafPriority.pt"
DEFAULT_SEG_MODEL_NAME = "Packera_LeafPriority"
DEFAULT_PCD_CONFIDENCE = 0.20
DEFAULT_NMS_THRESH = 0.75
DEFAULT_SUBDIVISION_STEPS = 5
DEFAULT_MIN_LEAF_AREA = 500
DEFAULT_SEG_CONFIDENCE = 0.70
DEFAULT_EFD_ORDER = 40


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure structured console logging."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("LM2_Configurator")


def resolve_pcd_weights(
    weights_name_or_path: Union[str, Path],
    logger: Optional[logging.Logger] = None,
) -> Tuple[str, Optional[Path]]:
    """
    Resolves the location of the Plant Component Detector weights file.
    Searches known model directories in priority order:
      1. Absolute or relative path passed directly
      2. LeafMachine2 embedded component detector weights
      3. models/LeafPriority.pt
      4. LM2_Project/Data/output/pcd_training/LeafPriority.pt
    """
    target_name = Path(weights_name_or_path).name
    p = Path(weights_name_or_path)

    candidates: List[Tuple[str, Path]] = []
    if p.is_absolute() and p.exists():
        candidates.append(("Direct Absolute Path", p))
    elif (PROJECT_ROOT / p).exists():
        candidates.append(("Project Relative Path", (PROJECT_ROOT / p).resolve()))

    candidates.extend([
        ("LM2 Component Detector Directory", LM2_ROOT / "leafmachine2" / "component_detector" / "runs" / "train" / "Plant_Detector" / "PLANT_LeafPriority" / "PLANT_LeafPriority" / "weights" / target_name),
        ("Project Models Directory", PROJECT_ROOT / "models" / target_name),
        ("PCD Training Output Directory", LM2_PROJECT_DIR / "Data" / "output" / "pcd_training" / target_name),
        ("Project Models Checkpoints", PROJECT_ROOT / "models" / "checkpoints" / target_name),
    ])

    for desc, cand in candidates:
        if cand.exists() and cand.is_file() and cand.stat().st_size > 1024:
            if logger:
                logger.debug(f"Resolved PCD weights via [{desc}]: {cand}")
            return target_name, cand

    if logger:
        logger.warning(f"Could not locate physical weights file for '{weights_name_or_path}' on disk.")
    return target_name, None


def sync_weights_to_lm2_environment(
    weights_name: str = DEFAULT_PCD_WEIGHTS_NAME,
    seg_model_name: str = DEFAULT_SEG_MODEL_NAME,
    logger: Optional[logging.Logger] = None,
) -> None:
    """
    Ensures that the PCD weights and trained PointRend segmentation artifacts
    (cfg_output.yaml, metadata.json, model_final.pth) are mirrored into
    LeafMachine2's component detector run directory and segmentation model tree.
    """
    # 1. Mirror PCD weights into LeafMachine2 component detector weights directory
    _, resolved_pcd = resolve_pcd_weights(weights_name, logger=logger)
    target_lm2_pcd_dir = (
        LM2_ROOT
        / "leafmachine2"
        / "component_detector"
        / "runs"
        / "train"
        / "Plant_Detector"
        / "PLANT_LeafPriority"
        / "PLANT_LeafPriority"
        / "weights"
    )
    target_lm2_pcd_dir.mkdir(parents=True, exist_ok=True)
    if resolved_pcd and resolved_pcd.exists():
        target_pcd_weights = target_lm2_pcd_dir / weights_name
        try:
            if not target_pcd_weights.exists() or target_pcd_weights.stat().st_size != resolved_pcd.stat().st_size:
                shutil.copyfile(resolved_pcd, target_pcd_weights)
                if logger:
                    logger.info(f"Mirrored PCD weights to LM2 detector tree: {target_pcd_weights}")
        except Exception as err:
            if logger:
                logger.debug(f"PCD weights sync notice: {err}")

    # 2. Mirror PointRend segmentation artifacts into LeafMachine2 segmentation models tree
    target_lm2_seg_dir = LM2_ROOT / "leafmachine2" / "segmentation" / "models" / seg_model_name
    target_lm2_seg_dir.mkdir(parents=True, exist_ok=True)

    pcd_training_dir = LM2_PROJECT_DIR / "Data" / "output" / "pcd_training"
    if pcd_training_dir.exists():
        for fname in ["cfg_output.yaml", "metadata.json", "model_final.pth"]:
            src_file = pcd_training_dir / fname
            if not src_file.exists() and fname == "model_final.pth":
                src_file = pcd_training_dir / "Packera_LeafPriority.pth"
            if src_file.exists():
                dst_file = target_lm2_seg_dir / fname
                try:
                    if not dst_file.exists() or dst_file.stat().st_size != src_file.stat().st_size or dst_file.stat().st_mtime < src_file.stat().st_mtime:
                        shutil.copyfile(src_file, dst_file)
                        if logger:
                            logger.info(f"Mirrored PointRend artifact {fname} to {target_lm2_seg_dir}")
                except Exception as err:
                    if logger:
                        logger.debug(f"Segmentation artifact sync notice: {err}")


def get_default_lm2_template() -> Dict[str, Any]:
    """
    Returns a complete, canonical LeafMachine2 configuration template dictionary
    matching the LeafMachine2 specification, pre-configured with Packera-specific
    thresholds and PointRend parameters.
    """
    return {
        "leafmachine": {
            "project": {
                "run_name": "Packera_dubia_LM2",
                "dir_images_local": str(DEFAULT_IMAGES_DIR),
                "dir_output": str(DEFAULT_OUTPUT_DIR),
                "dir_images_subset": "",
                "image_location": "local",
                "process_subset_of_images": False,
                "n_images_per_species": 1000000,
                "species_list": "",
                "GBIF_mode": "all",
                "batch_size": 50,
                "num_workers": 1,
                "num_workers_ruler": 8,
                "num_workers_seg": 8,
                "device": "cuda",
                "use_CF_predictor": True,
                "accept_only_ideal_leaves": True,
                "treat_leaflet_as_leaf_whole": False,
                "censor_archival_components": False,
                "hide_archival_components": [
                    "ruler",
                    "barcode",
                    "label",
                    "colorcard",
                    "map",
                    "photo",
                    "weights",
                ],
                "replacement_color": "#FFFFFF",
                "overlay_dpi": 300,
                "overlay_background_color": "black",
                "path_combined_csv_local": None,
                "path_images_csv_local": None,
                "path_occurrence_csv_local": None,
                "minimum_total_reproductive_counts": 0,
                "auto_cache_annotations": False,
                "regenerate_annotations": True,
                "use_existing_archival_component_detections": None,
                "use_existing_plant_component_detections": None,
            },
            "plant_component_detector": {
                "detector_type": "Plant_Detector",
                "detector_version": "PLANT_LeafPriority",
                "detector_iteration": "PLANT_LeafPriority",
                "detector_weights": DEFAULT_PCD_WEIGHTS_NAME,
                "minimum_confidence_threshold": DEFAULT_PCD_CONFIDENCE,
                "PCD_confidence": DEFAULT_PCD_CONFIDENCE,
                "do_save_prediction_overlay_images": True,
                "ignore_objects_for_overlay": [],
            },
            "archival_component_detector": {
                "detector_type": "Archival_Detector",
                "detector_version": "PREP_final",
                "detector_iteration": "PREP_final",
                "detector_weights": "best.pt",
                "minimum_confidence_threshold": 0.7,
                "ruler_binary_detector": "model_scripted_resnet_720_withCompression.pt",
                "ruler_detector": "ruler_classifier_38classes_v-1.pt",
                "do_save_prediction_overlay_images": True,
                "ignore_objects_for_overlay": [],
            },
            "armature_component_detector": {
                "detector_type": "Armature_Detector",
                "detector_version": "ARM_A_1000",
                "detector_iteration": "ARM_A_1000",
                "detector_weights": "best.pt",
                "minimum_confidence_threshold": 0.5,
                "do_save_prediction_overlay_images": True,
                "ignore_objects_for_overlay": [],
            },
            "landmark_detector": {
                "landmark_whole_leaves": True,
                "landmark_partial_leaves": False,
                "detector_type": "Landmark_Detector_YOLO",
                "detector_version": "Landmarks",
                "detector_iteration": "Landmarks_V2",
                "detector_weights": "best.pt",
                "minimum_confidence_threshold": 0.02,
                "do_save_prediction_overlay_images": True,
                "ignore_objects_for_overlay": [],
                "use_existing_landmark_detections": None,
                "do_show_QC_images": False,
                "do_save_QC_images": True,
                "do_show_final_images": False,
                "do_save_final_images": True,
            },
            "landmark_detector_armature": {
                "upscale_factor": 10,
                "detector_type": "Landmark_Detector_YOLO",
                "detector_version": "Landmarks_Arm_A_200",
                "detector_iteration": "Landmarks_Arm_A_200",
                "detector_weights": "last.pt",
                "minimum_confidence_threshold": 0.06,
                "do_save_prediction_overlay_images": True,
                "ignore_objects_for_overlay": [],
                "use_existing_landmark_detections": None,
                "do_show_QC_images": True,
                "do_save_QC_images": True,
                "do_show_final_images": True,
                "do_save_final_images": True,
            },
            "leaf_segmentation": {
                "segment_whole_leaves": True,
                "segment_partial_leaves": False,
                "segmentation_model": DEFAULT_SEG_MODEL_NAME,
                "detector_version": "uniform_spaced_oriented_traces_mid15_pet5_clean_640_flipidx_pt2",
                "minimum_confidence_threshold": DEFAULT_SEG_CONFIDENCE,
                "NMS_thresh": DEFAULT_NMS_THRESH,
                "min_leaf_area": DEFAULT_MIN_LEAF_AREA,
                "pointrend_subdivision_steps": DEFAULT_SUBDIVISION_STEPS,
                "keep_only_best_one_leaf_one_petiole": True,
                "find_minimum_bounding_box": True,
                "calculate_elliptic_fourier_descriptors": True,
                "elliptic_fourier_descriptor_order": DEFAULT_EFD_ORDER,
                "use_efds_for_png_masks": False,
                "save_masks_color": True,
                "save_full_image_masks_color": True,
                "save_rgb_cropped_images": True,
                "save_oriented_images": True,
                "save_oriented_mask": True,
                "save_keypoint_overlay": True,
                "save_simple_txt": True,
                "generate_overlay": True,
                "save_individual_overlay_images": True,
                "save_each_segmentation_overlay_image": True,
                "save_segmentation_overlay_images_to_pdf": False,
                "overlay_dpi": 300,
                "overlay_background_color": "black",
                "overlay_line_width": 1,
            },
            "ruler_detection": {
                "detect_ruler_type": True,
                "minimum_confidence_threshold": 0.5,
                "ruler_binary_detector": "model_scripted_resnet_720_withCompression.pt",
                "ruler_detector": "ruler_classifier_38classes_v-1.pt",
                "save_ruler_processed": False,
                "save_ruler_validation": False,
                "save_ruler_validation_summary": True,
            },
            "modules": {
                "specimen_crop": False,
                "armature": False,
            },
            "cropped_components": {
                "do_save_cropped_annotations": False,
                "save_per_image": False,
                "save_per_annotation_class": False,
                "save_cropped_annotations": ["label"],
                "binarize_labels": False,
                "binarize_labels_skeletonize": False,
            },
            "data": {
                "save_json_measurements": False,
                "save_individual_csv_files_measurements": True,
                "save_json_rulers": False,
                "save_individual_csv_files_rulers": True,
                "save_individual_csv_files_landmarks": True,
                "save_individual_efd_files": False,
                "include_darwin_core_data_from_combined_file": False,
                "do_apply_conversion_factor": True,
            },
            "do": {
                "check_for_corrupt_images_make_vertical": False,
                "check_for_illegal_filenames": False,
                "run_leaf_processing": True,
            },
            "overlay": {
                "save_overlay_to_pdf": False,
                "save_overlay_to_jpgs": True,
                "overlay_dpi": 300,
                "overlay_background_color": "black",
                "show_archival_detections": True,
                "show_plant_detections": True,
                "show_segmentations": True,
                "show_landmarks": True,
                "ignore_archival_detections_classes": [],
                "ignore_plant_detections_classes": ["leaf_whole"],
                "ignore_landmark_classes": [],
                "line_width_archival": 12,
                "line_width_plant": 12,
                "line_width_seg": 12,
                "line_width_efd": 12,
                "alpha_transparency_archival": 0.3,
                "alpha_transparency_plant": 0.0,
                "alpha_transparency_seg_whole_leaf": 0.4,
                "alpha_transparency_seg_partial_leaf": 0.3,
            },
            "print": {
                "verbose": True,
                "optional_warnings": True,
            },
            "logging": {
                "log_level": None,
            },
        }
    }


def load_config_yaml(yaml_path: Union[str, Path]) -> Dict[str, Any]:
    """Load an existing YAML configuration file."""
    path = Path(yaml_path)
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    return cfg


def save_config_yaml(cfg: Dict[str, Any], output_path: Union[str, Path]) -> Path:
    """Save configuration dictionary to a YAML file, creating parent directories if needed."""
    out_file = Path(output_path).resolve()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(out_file, "w", encoding="utf-8") as f:
            yaml.dump(cfg, f, default_flow_style=False, sort_keys=False, indent=2)
    except PermissionError as e:
        raise PermissionError(
            f"Permission denied writing to '{out_file}'. If updating system-protected files, "
            f"run with elevated permissions (e.g. sudo) or specify an output path inside LM2_Project/configs/."
        ) from e
    return out_file


def generate_high_performance_config(
    base_config: Optional[Union[str, Path, Dict[str, Any]]] = None,
    batch_size: int = 50,
    num_workers: int = 1,
    num_workers_ruler: int = 8,
    num_workers_seg: int = 8,
    device: str = "cuda",
    use_leaf_priority: bool = True,
    pcd_weights: str = DEFAULT_PCD_WEIGHTS_NAME,
    pcd_confidence: float = DEFAULT_PCD_CONFIDENCE,
    seg_model: str = DEFAULT_SEG_MODEL_NAME,
    nms_thresh: float = DEFAULT_NMS_THRESH,
    pointrend_subdivision_steps: int = DEFAULT_SUBDIVISION_STEPS,
    min_leaf_area: int = DEFAULT_MIN_LEAF_AREA,
    seg_confidence: float = DEFAULT_SEG_CONFIDENCE,
    images_dir: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    run_name: str = "Packera_dubia_LM2",
    save_overlays: bool = True,
    calculate_efds: bool = True,
    efd_order: int = DEFAULT_EFD_ORDER,
) -> Dict[str, Any]:
    """
    Generate or modify a LeafMachine2 configuration dictionary tailored for
    high-performance compute environments and dense Packera rosette morphology.

    Args:
        base_config: Path to existing YAML file, a dict, or None (uses default template).
        batch_size: Processing batch size (default: 50).
        num_workers: Parallel workers for component detection (default: 1 for single GPU).
        num_workers_ruler: Parallel workers for ruler conversion (default: 8).
        num_workers_seg: Parallel workers for leaf segmentation (default: 8).
        device: Compute hardware device ('cuda' or 'cpu', default: 'cuda').
        use_leaf_priority: Set PCD to 'LeafPriority' / custom weights (default: True).
        pcd_weights: Name or path to PCD YOLO weights (default: 'LeafPriority.pt').
        pcd_confidence: Component confidence threshold (default: 0.20 for dense rosettes).
        seg_model: PointRend segmentation model directory name (default: 'Packera_LeafPriority').
        nms_thresh: Non-Maximum Suppression overlap threshold (default: 0.75).
        pointrend_subdivision_steps: PointRend subdivision steps (default: 5 for complex margins).
        min_leaf_area: Minimum leaf pixel area to retain (default: 500).
        seg_confidence: Minimum confidence threshold for Leaf Segmentation (default: 0.70).
        images_dir: Directory containing input herbarium images / symlinks.
        output_dir: Directory to save LM2 output masks and measurements.
        run_name: Experiment or run identifier.
        save_overlays: Whether to generate and save JPG overlay visualizations.
        calculate_efds: Whether to compute Elliptic Fourier Descriptors in LM2.
        efd_order: Number of harmonics for EFD extraction (default: 40).

    Returns:
        Dict[str, Any]: Optimized LeafMachine2 configuration dictionary.
    """
    if base_config is None:
        cfg = get_default_lm2_template()
    elif isinstance(base_config, (str, Path)):
        if Path(base_config).exists():
            cfg = load_config_yaml(base_config)
        else:
            cfg = get_default_lm2_template()
    elif isinstance(base_config, dict):
        cfg = copy.deepcopy(base_config)
    else:
        cfg = get_default_lm2_template()

    # Ensure 'leafmachine' top-level namespace
    if "leafmachine" not in cfg:
        cfg = {"leafmachine": cfg}

    lm = cfg["leafmachine"]

    # 1. Project & Hardware Performance Settings
    if "project" not in lm:
        lm["project"] = {}
    proj = lm["project"]

    proj["batch_size"] = int(batch_size)
    proj["num_workers"] = int(num_workers)
    proj["num_workers_ruler"] = int(num_workers_ruler)
    proj["num_workers_seg"] = int(num_workers_seg)
    proj["device"] = str(device)
    proj["run_name"] = str(run_name)

    if images_dir is not None:
        proj["dir_images_local"] = str(Path(images_dir).resolve())
    elif "dir_images_local" not in proj or not proj["dir_images_local"]:
        proj["dir_images_local"] = str(DEFAULT_IMAGES_DIR.resolve())

    if output_dir is not None:
        proj["dir_output"] = str(Path(output_dir).resolve())
    elif "dir_output" not in proj or not proj["dir_output"]:
        proj["dir_output"] = str(DEFAULT_OUTPUT_DIR.resolve())

    # Fast local processing defaults
    proj["image_location"] = "local"
    proj["use_CF_predictor"] = True
    proj["accept_only_ideal_leaves"] = True
    proj["treat_leaflet_as_leaf_whole"] = False
    proj["auto_cache_annotations"] = False
    proj["regenerate_annotations"] = True

    # 2. Plant Component Detector: Packera_LeafPriority (PCD Version 2.2 Fine-Tuned)
    if "plant_component_detector" not in lm:
        lm["plant_component_detector"] = {}
    pcd = lm["plant_component_detector"]

    resolved_weights_name, _ = resolve_pcd_weights(pcd_weights)

    if use_leaf_priority:
        pcd["detector_type"] = "Plant_Detector"
        pcd["detector_version"] = "PLANT_LeafPriority"
        pcd["detector_iteration"] = "PLANT_LeafPriority"
        pcd["detector_weights"] = str(resolved_weights_name)

    # Lowered confidence threshold for dense, overlapping basal rosettes
    pcd["minimum_confidence_threshold"] = float(pcd_confidence)
    pcd["PCD_confidence"] = float(pcd_confidence)
    pcd["do_save_prediction_overlay_images"] = bool(save_overlays)

    # 3. Leaf Segmentation & PointRend Boundary Settings
    if "leaf_segmentation" not in lm:
        lm["leaf_segmentation"] = {}
    seg = lm["leaf_segmentation"]

    seg["segment_whole_leaves"] = True
    seg["segment_partial_leaves"] = False
    seg["segmentation_model"] = str(seg_model)
    seg["minimum_confidence_threshold"] = float(seg_confidence)

    # Modified thresholds for dense overlapping rosettes & tomentose margins
    seg["NMS_thresh"] = float(nms_thresh)
    seg["min_leaf_area"] = int(min_leaf_area)
    seg["pointrend_subdivision_steps"] = int(pointrend_subdivision_steps)

    seg["calculate_elliptic_fourier_descriptors"] = bool(calculate_efds)
    seg["elliptic_fourier_descriptor_order"] = int(efd_order)
    seg["find_minimum_bounding_box"] = True
    seg["keep_only_best_one_leaf_one_petiole"] = True

    seg["save_masks_color"] = True
    seg["save_full_image_masks_color"] = True
    seg["save_rgb_cropped_images"] = True
    seg["save_oriented_images"] = True
    seg["save_oriented_mask"] = True
    seg["save_keypoint_overlay"] = True
    seg["save_simple_txt"] = True
    seg["generate_overlay"] = True
    seg["save_individual_overlay_images"] = True
    seg["save_each_segmentation_overlay_image"] = True

    # 4. Processing Flags
    if "do" not in lm:
        lm["do"] = {}
    lm["do"]["run_leaf_processing"] = True
    lm["do"]["check_for_corrupt_images_make_vertical"] = False
    lm["do"]["check_for_illegal_filenames"] = False

    # 5. Overlays and Visualization
    if "overlay" not in lm:
        lm["overlay"] = {}
    lm["overlay"]["save_overlay_to_jpgs"] = bool(save_overlays)
    lm["overlay"]["save_overlay_to_pdf"] = False
    lm["overlay"]["overlay_dpi"] = 300
    lm["overlay"]["show_plant_detections"] = True
    lm["overlay"]["show_segmentations"] = True
    lm["overlay"]["show_landmarks"] = True
    lm["overlay"]["show_archival_detections"] = True

    # 6. Data Exports
    if "data" not in lm:
        lm["data"] = {}
    lm["data"]["save_individual_csv_files_measurements"] = True
    lm["data"]["save_individual_csv_files_rulers"] = True
    lm["data"]["save_individual_csv_files_landmarks"] = True
    lm["data"]["do_apply_conversion_factor"] = True

    # 7. Logging and Prints
    if "print" not in lm:
        lm["print"] = {}
    lm["print"]["verbose"] = True
    lm["print"]["optional_warnings"] = True

    return cfg


def print_config_summary(cfg: Dict[str, Any], logger: logging.Logger) -> None:
    """Log a human-readable summary of the key configuration parameters."""
    lm = cfg.get("leafmachine", {})
    proj = lm.get("project", {})
    pcd = lm.get("plant_component_detector", {})
    seg = lm.get("leaf_segmentation", {})

    logger.info("=================================================================")
    logger.info(" LeafMachine2 Packera High-Performance Configuration Summary")
    logger.info("=================================================================")
    logger.info(f"  Run Name:                       {proj.get('run_name')}")
    logger.info(f"  Compute Device:                 {proj.get('device')}")
    logger.info(f"  Batch Size:                     {proj.get('batch_size')}")
    logger.info(f"  Parallel Workers (General/PCD): {proj.get('num_workers')}")
    logger.info(f"  Parallel Workers (Ruler):       {proj.get('num_workers_ruler')}")
    logger.info(f"  Parallel Workers (Seg):         {proj.get('num_workers_seg')}")
    logger.info("  ---------------------------------------------------------------")
    logger.info(f"  PCD Detector Type:              {pcd.get('detector_type')}")
    logger.info(f"  PCD Detector Version:           {pcd.get('detector_version')}")
    logger.info(f"  PCD Detector Iteration:         {pcd.get('detector_iteration')}")
    logger.info(f"  PCD Weights File:               {pcd.get('detector_weights')}")
    logger.info(f"  PCD Confidence Threshold:       {pcd.get('minimum_confidence_threshold')} (Dense Rosette Optimized)")
    logger.info("  ---------------------------------------------------------------")
    logger.info(f"  Segmentation Model:             {seg.get('segmentation_model')}")
    logger.info(f"  NMS Overlap Threshold:          {seg.get('NMS_thresh')} (Overlapping Rosette Optimized)")
    logger.info(f"  PointRend Subdivision Steps:    {seg.get('pointrend_subdivision_steps')} (Max Resolution Margin Fidelity)")
    logger.info(f"  Minimum Leaf Pixel Area:        {seg.get('min_leaf_area')} px")
    logger.info(f"  Compute EFD Descriptors:        {seg.get('calculate_elliptic_fourier_descriptors')} (Order {seg.get('elliptic_fourier_descriptor_order')})")
    logger.info("  ---------------------------------------------------------------")
    logger.info(f"  Input Images Directory:         {proj.get('dir_images_local')}")
    logger.info(f"  Output Directory:               {proj.get('dir_output')}")
    logger.info(f"  Segment Whole Leaves:           {seg.get('segment_whole_leaves')}")
    logger.info(f"  Save Mask Visualizations:       {seg.get('save_masks_color')}")
    logger.info("=================================================================")


def run_leafmachine2(
    config: Optional[Union[str, Path, Dict[str, Any]]] = None,
    lm2_home: Optional[Union[str, Path]] = None,
    **config_kwargs,
) -> Any:
    """
    Programmatic entry point to run LeafMachine2 pipeline using the high-performance
    configuration dictionary or YAML path.

    Args:
        config: Path to YAML or dict with 'leafmachine' config. If None, generates one.
        lm2_home: Path to LeafMachine2 repository root.
        **config_kwargs: Passed to generate_high_performance_config() if config is None.

    Returns:
        Result from LeafMachine2 pipeline execution.
    """
    if config is None:
        cfg = generate_high_performance_config(**config_kwargs)
    elif isinstance(config, (str, Path)):
        cfg = load_config_yaml(config)
    elif isinstance(config, dict):
        cfg = config
    else:
        cfg = generate_high_performance_config(**config_kwargs)

    lm2_path = Path(lm2_home) if lm2_home else LM2_ROOT
    if str(lm2_path) not in sys.path:
        sys.path.insert(0, str(lm2_path))

    from LeafMachine2.leafmachine2.machine.LeafMachine2 import LeafMachine2

    # Instantiate and run
    engine = LeafMachine2(cfg)
    return engine.run()


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate, validate, and deploy LeafMachine2 configuration for Packera dubia morphometrics.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-c", "--base-config",
        type=Path,
        default=DEFAULT_BASE_YAML,
        help="Path to base LeafMachine2 YAML configuration template.",
    )
    parser.add_argument(
        "-o", "--output-config",
        type=Path,
        default=DEFAULT_OUTPUT_YAML,
        help="Path to write the generated YAML configuration.",
    )
    parser.add_argument(
        "--update-main-config",
        action="store_true",
        help="Also overwrite LeafMachine2/LeafMachine2.yaml with the generated configuration.",
    )
    parser.add_argument(
        "--pcd-weights",
        type=str,
        default=DEFAULT_PCD_WEIGHTS_NAME,
        help=f"Plant Component Detector YOLOv5 weights filename or path (default: '{DEFAULT_PCD_WEIGHTS_NAME}').",
    )
    parser.add_argument(
        "--seg-model",
        type=str,
        default=DEFAULT_SEG_MODEL_NAME,
        help=f"Detectron2 PointRend segmentation model directory name (default: '{DEFAULT_SEG_MODEL_NAME}').",
    )
    parser.add_argument(
        "--pcd-confidence",
        type=float,
        default=DEFAULT_PCD_CONFIDENCE,
        help=f"Plant Component Detector confidence threshold (default: {DEFAULT_PCD_CONFIDENCE}).",
    )
    parser.add_argument(
        "--nms-thresh",
        type=float,
        default=DEFAULT_NMS_THRESH,
        help=f"Non-Maximum Suppression (NMS) overlap threshold (default: {DEFAULT_NMS_THRESH}).",
    )
    parser.add_argument(
        "--pointrend-subdivision-steps",
        type=int,
        default=DEFAULT_SUBDIVISION_STEPS,
        help=f"PointRend subdivision steps for margin extraction (default: {DEFAULT_SUBDIVISION_STEPS}).",
    )
    parser.add_argument(
        "--min-leaf-area",
        type=int,
        default=DEFAULT_MIN_LEAF_AREA,
        help=f"Minimum leaf pixel area to retain (default: {DEFAULT_MIN_LEAF_AREA}).",
    )
    parser.add_argument(
        "--seg-confidence",
        type=float,
        default=DEFAULT_SEG_CONFIDENCE,
        help=f"Leaf segmentation confidence threshold (default: {DEFAULT_SEG_CONFIDENCE}).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Batch size for model inference (default: 50).",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=1,
        help="Number of parallel worker processes for component detector (default: 1 for single GPU).",
    )
    parser.add_argument(
        "--num-workers-ruler",
        type=int,
        default=8,
        help="Number of parallel worker processes for ruler conversion (default: 8).",
    )
    parser.add_argument(
        "--num-workers-seg",
        type=int,
        default=8,
        help="Number of parallel worker processes for leaf segmentation (default: 8).",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda",
        choices=["cuda", "cpu"],
        help="Compute device for PyTorch inference (default: 'cuda').",
    )
    parser.add_argument(
        "--images-dir",
        type=Path,
        default=DEFAULT_IMAGES_DIR,
        help=f"Path to input image dataset directory (default: {DEFAULT_IMAGES_DIR}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Path to output directory for LM2 results (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--run-name",
        type=str,
        default="Packera_dubia_LM2",
        help="Identifier name for this processing run (default: 'Packera_dubia_LM2').",
    )
    parser.add_argument(
        "--no-leaf-priority",
        dest="use_leaf_priority",
        action="store_false",
        help="Do not enforce LeafPriority detector (use base config detector).",
    )
    parser.add_argument(
        "--run-lm2",
        action="store_true",
        help="Execute LeafMachine2 immediately using the generated configuration.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Generate and validate configuration in memory without writing to disk.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose DEBUG logging.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    logger = setup_logging(verbose=args.verbose)

    logger.info("Configuring LeafMachine2 for Packera dubia rosette extraction...")

    # Ensure weights are synced across LM2 environment directories
    sync_weights_to_lm2_environment(weights_name=args.pcd_weights, seg_model_name=args.seg_model, logger=logger)

    cfg = generate_high_performance_config(
        base_config=args.base_config,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        num_workers_ruler=args.num_workers_ruler,
        num_workers_seg=args.num_workers_seg,
        device=args.device,
        use_leaf_priority=args.use_leaf_priority,
        pcd_weights=args.pcd_weights,
        pcd_confidence=args.pcd_confidence,
        seg_model=args.seg_model,
        nms_thresh=args.nms_thresh,
        pointrend_subdivision_steps=args.pointrend_subdivision_steps,
        min_leaf_area=args.min_leaf_area,
        seg_confidence=args.seg_confidence,
        images_dir=args.images_dir,
        output_dir=args.output_dir,
        run_name=args.run_name,
    )

    print_config_summary(cfg, logger)

    if args.dry_run:
        logger.info("[DRY-RUN] Configuration generated and validated successfully. File write skipped.")
        return

    # Write target project configuration (e.g. LM2_Project/configs/lm2_packera_highperf.yaml)
    out_path = save_config_yaml(cfg, args.output_config)
    logger.info(f"Saved optimized configuration to: {out_path}")

    # Programmatically update main LeafMachine2/LeafMachine2.yaml
    if args.update_main_config:
        main_yaml = LM2_ROOT / "LeafMachine2.yaml"
        try:
            save_config_yaml(cfg, main_yaml)
            logger.info(f"Overrode default LeafMachine2 configuration: {main_yaml}")
        except PermissionError:
            logger.warning(
                f"Could not overwrite '{main_yaml}' due to permissions. "
                f"Use sudo to update it directly, or pass -c '{out_path}' to LeafMachine2."
            )

    if args.run_lm2:
        logger.info("Launching LeafMachine2 execution engine...")
        run_leafmachine2(config=cfg)


if __name__ == "__main__":
    main()
