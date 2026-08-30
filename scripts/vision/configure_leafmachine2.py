#!/usr/bin/env python3
"""
===============================================================================
Script: configure_leafmachine2.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Programmatically creates, validates, and updates LeafMachine2 (LM2)
    configuration files and Python dictionaries optimized for high-performance
    compute hardware (48GB+ RAM, 8GB+ VRAM GPU).

    Key Defaults Configured:
      - Batch size: 50 (optimal throughput for 8GB+ VRAM without CUDA OOM)
      - Parallel workers: 8 (num_workers, num_workers_ruler, num_workers_seg)
      - Compute device: 'cuda' (GPU acceleration)
      - Plant Component Detector (PCD): 'LeafPriority' (Version 2.2 detector)
          detector_type: Plant_Detector
          detector_version: PLANT_LeafPriority
          detector_iteration: PLANT_LeafPriority
          detector_weights: LeafPriority.pt
      - Integrated paths:
          dir_images_local: LM2_Project/Data/images
          dir_output: LM2_Project/Data/output
          run_name: Packera_dubia_LM2

Usage:
    # Generate optimized YAML in LM2_Project/configs/
    python scripts/vision/configure_leafmachine2.py

    # Update LeafMachine2 default config in-place
    python scripts/vision/configure_leafmachine2.py --update-main-config

    # Custom output path and batch size
    python scripts/vision/configure_leafmachine2.py -o LM2_Project/configs/custom.yaml --batch-size 50 --num-workers 8

    # Programmatic Python API:
    from scripts.vision.configure_leafmachine2 import generate_high_performance_config, save_config_yaml
    cfg = generate_high_performance_config(batch_size=50, num_workers=8)
    save_config_yaml(cfg, "LM2_Project/configs/lm2_highperf.yaml")
===============================================================================
"""

import argparse
import copy
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Union
import yaml

# Determine project root and LeafMachine2 root paths
PROJECT_ROOT = Path(__file__).resolve().parents[2]
LM2_ROOT = PROJECT_ROOT / "LeafMachine2"
LM2_PROJECT_DIR = PROJECT_ROOT / "LM2_Project"

DEFAULT_BASE_YAML = LM2_ROOT / "LeafMachine2.yaml"
DEFAULT_OUTPUT_YAML = LM2_PROJECT_DIR / "configs" / "lm2_packera_highperf.yaml"
DEFAULT_IMAGES_DIR = LM2_PROJECT_DIR / "Data" / "images"
DEFAULT_OUTPUT_DIR = LM2_PROJECT_DIR / "Data" / "output"


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure structured console logging."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("LM2_Configurator")


def get_default_lm2_template() -> Dict[str, Any]:
    """
    Returns a complete, canonical LeafMachine2 configuration template dictionary
    matching the LeafMachine2 specification.
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
                "num_workers": 8,
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
                "use_existing_archival_component_detections": None,
                "use_existing_plant_component_detections": None,
            },
            "plant_component_detector": {
                "detector_type": "Plant_Detector",
                "detector_version": "PLANT_LeafPriority",
                "detector_iteration": "PLANT_LeafPriority",
                "detector_weights": "LeafPriority.pt",
                "minimum_confidence_threshold": 0.5,
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
                "detector_iteration": "Landmarks",
                "detector_weights": "best.pt",
                "minimum_confidence_threshold": 0.02,
                "do_save_prediction_overlay_images": False,
                "ignore_objects_for_overlay": [],
                "use_existing_landmark_detections": None,
                "do_show_QC_images": False,
                "do_save_QC_images": False,
                "do_show_final_images": False,
                "do_save_final_images": False,
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
                "segmentation_model": "Group3_Dataset_100000_Iter_1176PTS_512Batch_smooth_l1_LR00025_BGR",
                "detector_version": "uniform_spaced_oriented_traces_mid15_pet5_clean_640_flipidx_pt2",
                "minimum_confidence_threshold": 0.7,
                "keep_only_best_one_leaf_one_petiole": True,
                "find_minimum_bounding_box": True,
                "calculate_elliptic_fourier_descriptors": False,
                "elliptic_fourier_descriptor_order": 40,
                "use_efds_for_png_masks": False,
                "save_masks_color": True,
                "save_full_image_masks_color": True,
                "save_rgb_cropped_images": True,
                "save_oriented_images": True,
                "save_oriented_mask": True,
                "save_keypoint_overlay": True,
                "save_simple_txt": True,
                "generate_overlay": False,
                "save_individual_overlay_images": False,
                "save_each_segmentation_overlay_image": False,
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
                "save_individual_csv_files_measurements": False,
                "save_json_rulers": False,
                "save_individual_csv_files_rulers": False,
                "save_individual_csv_files_landmarks": False,
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


def run_leafmachine2(
    config: Optional[Union[str, Path, Dict[str, Any]]] = None,
    lm2_home: Optional[Union[str, Path]] = None,
    **config_kwargs,
) -> Any:
    """
    Programmatically launch LeafMachine2 with high-performance configuration.

    Args:
        config: Path to YAML file or configuration dictionary. If None, generates
                a high-performance configuration automatically.
        lm2_home: Path to LeafMachine2 package root directory.
        **config_kwargs: Additional configuration parameters passed to
                         generate_high_performance_config().

    Returns:
        Result from leafmachine2 machine execution.
    """
    home = Path(lm2_home or LM2_ROOT).resolve()
    if not home.exists():
        raise FileNotFoundError(f"LeafMachine2 root directory not found at: {home}")

    # Add LeafMachine2 and internal directories to Python path
    if str(home) not in sys.path:
        sys.path.insert(0, str(home))

    try:
        from leafmachine2.machine.machine import machine
    except ImportError as e:
        raise ImportError(
            f"Could not import leafmachine2. Ensure the LM2 virtual environment is activated. Error: {e}"
        ) from e

    if isinstance(config, (str, Path)):
        cfg_file_path = str(Path(config).resolve())
        cfg_test = None
    elif isinstance(config, dict):
        cfg_file_path = None
        cfg_test = config
    else:
        cfg_dict = generate_high_performance_config(**config_kwargs)
        cfg_file_path = None
        cfg_test = cfg_dict

    return machine(cfg_file_path=cfg_file_path, dir_home=str(home), cfg_test=cfg_test)



def generate_high_performance_config(
    base_config: Optional[Union[str, Path, Dict[str, Any]]] = None,
    batch_size: int = 50,
    num_workers: int = 8,
    device: str = "cuda",
    use_leaf_priority: bool = True,
    images_dir: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    run_name: str = "Packera_dubia_LM2",
    pcd_confidence: float = 0.5,
    seg_confidence: float = 0.7,
    save_overlays: bool = True,
) -> Dict[str, Any]:
    """
    Generate or modify a LeafMachine2 configuration dictionary optimized for
    high-performance Ubuntu compute environments (48GB+ RAM, 8GB+ VRAM).

    Args:
        base_config: Path to existing YAML file, or a dict, or None (uses default template).
        batch_size: Processing batch size (default: 50).
        num_workers: Parallel workers for processing, ruler, and seg (default: 8).
        device: Compute hardware device ('cuda' or 'cpu', default: 'cuda').
        use_leaf_priority: Set PCD to 'LeafPriority' (Version 2.2) (default: True).
        images_dir: Directory containing input herbarium images / symlinks.
        output_dir: Directory to save LM2 output masks and measurements.
        run_name: Experiment or run identifier.
        pcd_confidence: Minimum confidence threshold for Plant Component Detector.
        seg_confidence: Minimum confidence threshold for Leaf Segmentation.
        save_overlays: Whether to generate and save JPG overlay visualizations.

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

    # 1. Project & Hardware Performance Settings (48GB+ RAM, 8GB+ VRAM)
    if "project" not in lm:
        lm["project"] = {}
    proj = lm["project"]

    proj["batch_size"] = int(batch_size)
    proj["num_workers"] = int(num_workers)
    proj["num_workers_ruler"] = int(num_workers)
    proj["num_workers_seg"] = int(num_workers)
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

    # 2. Plant Component Detector: LeafPriority (PCD Version 2.2)
    if "plant_component_detector" not in lm:
        lm["plant_component_detector"] = {}
    pcd = lm["plant_component_detector"]

    if use_leaf_priority:
        pcd["detector_type"] = "Plant_Detector"
        pcd["detector_version"] = "PLANT_LeafPriority"
        pcd["detector_iteration"] = "PLANT_LeafPriority"
        pcd["detector_weights"] = "LeafPriority.pt"
    pcd["minimum_confidence_threshold"] = float(pcd_confidence)
    pcd["do_save_prediction_overlay_images"] = bool(save_overlays)

    # 3. Leaf Segmentation Settings
    if "leaf_segmentation" not in lm:
        lm["leaf_segmentation"] = {}
    seg = lm["leaf_segmentation"]
    seg["segment_whole_leaves"] = True
    seg["segment_partial_leaves"] = False
    seg["minimum_confidence_threshold"] = float(seg_confidence)
    seg["save_masks_color"] = True
    seg["save_full_image_masks_color"] = True
    seg["save_rgb_cropped_images"] = True
    seg["save_oriented_images"] = True
    seg["save_oriented_mask"] = True
    seg["save_keypoint_overlay"] = True
    seg["save_simple_txt"] = True

    # 4. Processing Flags (Do not re-rotate or mutate filenames if already sanitized)
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

    # 6. Logging and Prints
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
    logger.info(" LeafMachine2 High-Performance Configuration Summary")
    logger.info("=================================================================")
    logger.info(f"  Run Name:                   {proj.get('run_name')}")
    logger.info(f"  Compute Device:             {proj.get('device')}")
    logger.info(f"  Batch Size:                 {proj.get('batch_size')}")
    logger.info(f"  Parallel Workers (General): {proj.get('num_workers')}")
    logger.info(f"  Parallel Workers (Ruler):   {proj.get('num_workers_ruler')}")
    logger.info(f"  Parallel Workers (Seg):     {proj.get('num_workers_seg')}")
    logger.info("  ---------------------------------------------------------------")
    logger.info(f"  PCD Detector Type:          {pcd.get('detector_type')}")
    logger.info(f"  PCD Detector Version:       {pcd.get('detector_version')}")
    logger.info(f"  PCD Detector Iteration:     {pcd.get('detector_iteration')}")
    logger.info(f"  PCD Weights File:           {pcd.get('detector_weights')}")
    logger.info(f"  PCD Confidence Threshold:   {pcd.get('minimum_confidence_threshold')}")
    logger.info("  ---------------------------------------------------------------")
    logger.info(f"  Input Images Directory:     {proj.get('dir_images_local')}")
    logger.info(f"  Output Directory:           {proj.get('dir_output')}")
    logger.info(f"  Segment Whole Leaves:       {seg.get('segment_whole_leaves')}")
    logger.info(f"  Save Mask Visualizations:   {seg.get('save_masks_color')}")
    logger.info("=================================================================")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="LeafMachine2 Configuration Generator & Optimizer (High Performance)."
    )
    parser.add_argument(
        "--base-config", "-b",
        type=Path,
        default=DEFAULT_BASE_YAML if DEFAULT_BASE_YAML.exists() else None,
        help="Path to source LeafMachine2 YAML config to modify (default: LeafMachine2/LeafMachine2.yaml).",
    )
    parser.add_argument(
        "--output-config", "-o",
        type=Path,
        default=DEFAULT_OUTPUT_YAML,
        help=f"Target path for generated YAML configuration (default: {DEFAULT_OUTPUT_YAML}).",
    )
    parser.add_argument(
        "--update-main-config",
        action="store_true",
        help="Also write changes directly to LeafMachine2/LeafMachine2.yaml.",
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
        default=8,
        help="Number of parallel worker processes/threads (default: 8).",
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
        "--pcd-confidence",
        type=float,
        default=0.5,
        help="Plant Component Detector confidence threshold (default: 0.5).",
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

    logger.info("Generating optimized LeafMachine2 configuration...")

    cfg = generate_high_performance_config(
        base_config=args.base_config,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        device=args.device,
        use_leaf_priority=args.use_leaf_priority,
        images_dir=args.images_dir,
        output_dir=args.output_dir,
        run_name=args.run_name,
        pcd_confidence=args.pcd_confidence,
    )

    print_config_summary(cfg, logger)

    if args.dry_run:
        logger.info("[DRY-RUN] Configuration generated successfully. File write skipped.")
        return

    # Write target config
    out_path = save_config_yaml(cfg, args.output_config)
    logger.info(f"Saved optimized configuration to: {out_path}")

    # Optionally update main LeafMachine2.yaml
    if args.update_main_config:
        main_yaml = LM2_ROOT / "LeafMachine2.yaml"
        try:
            save_config_yaml(cfg, main_yaml)
            logger.info(f"Updated default LeafMachine2 configuration: {main_yaml}")
        except PermissionError:
            logger.warning(
                f"Could not overwrite root-owned '{main_yaml}'. "
                f"Use sudo to update it directly, or pass -c '{out_path}' to LeafMachine2."
            )

    if args.run_lm2:
        logger.info("Launching LeafMachine2 execution engine...")
        run_leafmachine2(config=cfg)


if __name__ == "__main__":
    main()

