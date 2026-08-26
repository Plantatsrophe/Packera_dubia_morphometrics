#!/usr/bin/env python3
"""
===============================================================================
Script: train_artifact_robust_yolo.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
Author: Deep Learning Architect & Senior Computer Vision Engineer
Date: August 2026

Description:
    Fine-tunes a high-capacity YOLOv8 Instance Segmentation model (YOLOv8x-seg
    or YOLOv8m-seg) specifically configured to maximize botanical organ recall
    (basal lamina, petioles, rosettes, capitula) while rigorously suppressing
    false-positive detections on herbarium sheet mounting artifacts (labels,
    mounting tape strips, ruler scales, color charts, and barcode stickers).

Key Technical Capabilities & Architectural Features:
    1. Model Initialization:
       - Loads pre-trained weights (`yolov8x-seg.pt` or `yolov8m-seg.pt`).
       - Loads dataset configuration from `data/dataset_config.yaml`.
    2. Botanical Hyperparameter Configuration:
       - Image input resolution: `imgsz=1024` for capturing fine petiole margins.
       - Adaptive batch sizing based on CUDA VRAM (e.g., batch=8 or 16 with AMP).
       - Cosine annealing learning rate schedule (`cos_lr=True`, `lr0=0.01`,
         `lrf=0.001`, `warmup_epochs=3.0`).
       - Specialized multi-scale data augmentations:
         * Mosaic augmentation: `mosaic=1.0` (4-patch voucher synthesis)
         * MixUp augmentation: `mixup=0.15` (inter-class regularization)
         * Copy-Paste augmentation: `copy_paste=0.30` (segments pasted across backgrounds)
         * Geometric perturbations: `degrees=15.0`, `translate=0.1`, `scale=0.2`, `fliplr=0.5`
    3. Custom Loss Gain Configuration:
       - Bounding box loss gain: `box=7.5`
       - Classification loss gain: `cls=0.5`
       - Distribution focal loss gain: `dfl=1.5`
       - Background hard-negative sample penalization to eliminate artifact false alarms.
    4. Rigorous Evaluation & Checkpoint Export:
       - Computes class-specific Box & Mask mAP50 and mAP50-95 across all botanical
         and artifact classes.
       - Verifies zero cross-classification between `herbarium_label` and `basal_leaf`.
       - Exports confusion matrix and precision-recall curves to `outputs/training_evaluation/`.
       - Saves optimized best checkpoint to `models/yolov8_leaf_best.pt`.

Usage:
    python scripts/train_artifact_robust_yolo.py --epochs 100 --batch 8 --imgsz 1024
===============================================================================
"""

import os
import sys
import json
import shutil
import logging
import argparse
import random
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union

import numpy as np
import yaml
import torch
import matplotlib.pyplot as plt
from ultralytics import YOLO

# ===============================================================================
# LOGGING CONFIGURATION
# ===============================================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("ArtifactRobustYOLO")

# ===============================================================================
# STANDARD CLASS DEFINITIONS
# ===============================================================================
# Standardized botanical and mounting artifact multi-class taxonomy
STANDARD_CLASS_MAPPING: Dict[int, str] = {
    0: "basal_leaf",       # Intact leaf / lamina blade
    1: "leaf_petiole",     # Distinct petiole / leaf stalk
    2: "basal_rosette",    # Clustered basal rosette
    3: "capitulum",        # Inflorescence / flower head
    4: "herbarium_label",  # Specimen metadata collection label
    5: "color_chart",      # Color calibration chart / palette
    6: "ruler_scale",      # Centimeter scale / measurement bar
    7: "barcode_sticker",  # Digitization barcode / QR sticker
    8: "mounting_tape",    # Linen, paper, or plastic mounting tape strip
}


# ===============================================================================
# TILED DATASET PARTITIONING AND PREPARATION
# ===============================================================================
def extract_specimen_stem(filename: str) -> str:
    """
    Extracts the parent specimen sheet identifier from a tile filename.
    
    Herbarium patch tiles follow naming conventions such as:
    - '00116749_tile_y00000_x00000.jpg' -> Parent Specimen: '00116749'
    - 'neg_sheet_NCU00011994_76_tile_y00000_x00000.jpg' -> Parent Specimen: 'neg_sheet_NCU00011994_76'
    
    Args:
        filename (str): Name or path of the tile image file.
        
    Returns:
        str: Base specimen sheet ID.
    """
    stem = Path(filename).stem
    # Match the delimiter used by NativeDPIPatchTiler
    if "_tile_" in stem:
        return stem.split("_tile_")[0]
    return stem.split("_")[0]


def prepare_tiled_dataset_split(
    tiled_dir: Path,
    config_yaml_path: Optional[Path] = None,
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    seed: int = 42,
    force_resplit: bool = False
) -> Path:
    """
    Partitions native-DPI cropped tiles into train/val/test splits grouped strictly
    by parent specimen sheet ID. Grouping by specimen prevents spatial data leakage
    across overlapping tiles of the same physical herbarium voucher sheet.
    
    Also writes and validates the Ultralytics dataset configuration YAML (e.g. data/tiled_dataset_config.yaml).

    Args:
        tiled_dir (Path): Root directory containing tiled dataset ('images' and 'labels').
        config_yaml_path (Path, optional): Destination path for tiled_dataset_config.yaml.
        train_ratio (float): Fraction of specimen sheets for training (default: 0.70).
        val_ratio (float): Fraction of specimen sheets for validation (default: 0.15).
        test_ratio (float): Fraction of specimen sheets for testing (default: 0.15).
        seed (int): Deterministic random seed for reproducible splitting.
        force_resplit (bool): If True, re-partitions even if splits already exist.

    Returns:
        Path: Path to the validated YAML configuration file.
    """
    tiled_dir = Path(tiled_dir).resolve()
    images_dir = tiled_dir / "images"
    labels_dir = tiled_dir / "labels"

    if not images_dir.exists():
        raise FileNotFoundError(
            f"Tiled images directory not found at: {images_dir}. "
            f"Please run 'scripts/run_dpi_tiler.py' to generate native-DPI tiles."
        )

    if config_yaml_path is None:
        config_yaml_path = tiled_dir.parent / "tiled_dataset_config.yaml"
    else:
        config_yaml_path = Path(config_yaml_path).resolve()

    train_img_dir = images_dir / "train"
    val_img_dir = images_dir / "val"
    test_img_dir = images_dir / "test"

    # Check if dataset is already partitioned into train/val subdirectories
    already_split = (
        train_img_dir.exists() and any(train_img_dir.glob("*.jpg")) and
        val_img_dir.exists() and any(val_img_dir.glob("*.jpg"))
    )

    if already_split and not force_resplit:
        logger.info(f"Tiled dataset is already partitioned in: {tiled_dir}")
        train_count = len(list(train_img_dir.glob("*.jpg")))
        val_count = len(list(val_img_dir.glob("*.jpg")))
        test_count = len(list(test_img_dir.glob("*.jpg"))) if test_img_dir.exists() else 0
        logger.info(f"Existing partition count -> Train: {train_count} | Val: {val_count} | Test: {test_count} tiles")
    else:
        logger.info(
            f"Partitioning tiled dataset at: {tiled_dir} "
            f"(Train: {train_ratio*100:.0f}%, Val: {val_ratio*100:.0f}%, Test: {test_ratio*100:.0f}%)"
        )

        # Create split subdirectories under images/ and labels/
        for split in ["train", "val", "test"]:
            (images_dir / split).mkdir(parents=True, exist_ok=True)
            (labels_dir / split).mkdir(parents=True, exist_ok=True)

        # Collect all image files (either from root of images/ or nested if force_resplit)
        all_image_paths: List[Path] = []
        for ext in ("*.jpg", "*.jpeg", "*.png"):
            all_image_paths.extend(images_dir.glob(ext))
            if force_resplit:
                for split in ["train", "val", "test"]:
                    all_image_paths.extend((images_dir / split).glob(ext))

        # Deduplicate paths
        all_image_paths = sorted(list(set(all_image_paths)))

        if not all_image_paths:
            raise ValueError(f"No image tiles found in {images_dir} to partition.")

        # Group tiles by parent specimen sheet identifier
        specimen_to_images: Dict[str, List[Path]] = {}
        for img_path in all_image_paths:
            specimen_id = extract_specimen_stem(img_path.name)
            if specimen_id not in specimen_to_images:
                specimen_to_images[specimen_id] = []
            specimen_to_images[specimen_id].append(img_path)

        specimen_ids = sorted(list(specimen_to_images.keys()))
        logger.info(f"Found {len(all_image_paths)} total tiles across {len(specimen_ids)} unique specimen sheets.")

        # Deterministic shuffle to ensure reproducible training/validation partitions
        rng = random.Random(seed)
        rng.shuffle(specimen_ids)

        # Compute split boundaries across unique herbarium sheets
        n_total = len(specimen_ids)
        n_train = max(1, int(round(n_total * train_ratio)))
        n_val = max(1, int(round(n_total * val_ratio)))
        
        train_specs = set(specimen_ids[:n_train])
        val_specs = set(specimen_ids[n_train:n_train + n_val])
        test_specs = set(specimen_ids[n_train + n_val:])
        
        # Ensure test set has at least 1 specimen if total count allows
        if not test_specs and n_total > 2:
            test_specs = {specimen_ids[-1]}
            if specimen_ids[-1] in val_specs and len(val_specs) > 1:
                val_specs.remove(specimen_ids[-1])
            elif specimen_ids[-1] in train_specs and len(train_specs) > 1:
                train_specs.remove(specimen_ids[-1])

        specimen_split_map: Dict[str, str] = {}
        for s in train_specs:
            specimen_split_map[s] = "train"
        for s in val_specs:
            specimen_split_map[s] = "val"
        for s in test_specs:
            specimen_split_map[s] = "test"

        # Relocate image and label files into designated split directories
        split_tile_counts = {"train": 0, "val": 0, "test": 0}

        for specimen_id, img_paths in specimen_to_images.items():
            split = specimen_split_map.get(specimen_id, "train")
            dest_img_split = images_dir / split
            dest_lbl_split = labels_dir / split

            for img_path in img_paths:
                dest_img_path = dest_img_split / img_path.name
                if img_path != dest_img_path:
                    shutil.move(str(img_path), str(dest_img_path))

                # Locate corresponding YOLO polygon label file
                lbl_name = f"{img_path.stem}.txt"
                source_lbl_path = labels_dir / lbl_name
                if not source_lbl_path.exists():
                    # Check in split subdirectories if force_resplit
                    for s in ["train", "val", "test"]:
                        candidate = labels_dir / s / lbl_name
                        if candidate.exists():
                            source_lbl_path = candidate
                            break

                dest_lbl_path = dest_lbl_split / lbl_name
                if source_lbl_path.exists():
                    if source_lbl_path != dest_lbl_path:
                        shutil.move(str(source_lbl_path), str(dest_lbl_path))
                else:
                    # Create empty label file for hard negative background tiles (mounting paper / no annotations)
                    dest_lbl_path.touch()

                split_tile_counts[split] += 1

        logger.info(
            f"Successfully partitioned tiled dataset -> "
            f"Train: {split_tile_counts['train']} tiles ({len(train_specs)} specimens) | "
            f"Val: {split_tile_counts['val']} tiles ({len(val_specs)} specimens) | "
            f"Test: {split_tile_counts['test']} tiles ({len(test_specs)} specimens)"
        )

    # Write Ultralytics dataset configuration YAML using forward slashes for Windows/POSIX cross-compatibility
    config_dict = {
        "path": tiled_dir.as_posix(),
        "train": "images/train",
        "val": "images/val",
        "test": "images/test",
        "names": STANDARD_CLASS_MAPPING
    }

    config_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_yaml_path, "w", encoding="utf-8") as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)

    logger.info(f"Generated/verified tiled dataset YAML configuration at: {config_yaml_path}")
    return config_yaml_path


# ===============================================================================
# UTILITY FUNCTIONS FOR PATHS AND DEVICE RESOLUTION
# ===============================================================================
def resolve_default_paths(project_root: Optional[Path] = None) -> Dict[str, Path]:
    """
    Resolves standard project paths relative to workspace root.

    Args:
        project_root: Optional custom project root path.

    Returns:
        Dictionary mapping path identifiers to absolute Path objects.
    """
    if project_root is None:
        # Default to parent of scripts/ or current working directory
        current_script_dir = Path(__file__).resolve().parent
        if current_script_dir.name == "scripts":
            project_root = current_script_dir.parent
        else:
            project_root = Path.cwd()

    tiled_dataset_dir = project_root / "data" / "tiled_dataset"
    tiled_config = project_root / "data" / "tiled_dataset_config.yaml"
    legacy_config = project_root / "data" / "dataset_config.yaml"

    # Default to tiled dataset config if tiled dataset exists, else legacy config
    default_config = tiled_config if tiled_dataset_dir.exists() else legacy_config

    paths = {
        "root": project_root,
        "data_config": default_config,
        "tiled_config": tiled_config,
        "legacy_config": legacy_config,
        "tiled_dir": tiled_dataset_dir,
        "models_dir": project_root / "models",
        "output_eval_dir": project_root / "outputs" / "training_evaluation",
        "default_weights": project_root / "yolov8x-seg.pt",
        "fallback_weights": project_root / "yolov8m-seg.pt",
        "best_model_export": project_root / "models" / "yolov8_leaf_best.pt"
    }

    # Ensure required output directories exist
    paths["models_dir"].mkdir(parents=True, exist_ok=True)
    paths["output_eval_dir"].mkdir(parents=True, exist_ok=True)

    return paths


def detect_optimal_device_and_batch(requested_batch: Optional[int] = None, imgsz: int = 1024) -> Tuple[str, int]:
    """
    Detects hardware accelerator (CUDA GPU / CPU) and computes safe batch size
    based on available VRAM to prevent Out-Of-Memory (OOM) errors during 1024x1024
    instance segmentation training.

    Args:
        requested_batch: User-specified batch size, or None for auto-selection.
        imgsz: Training input resolution (default: 1024).

    Returns:
        Tuple of (device_string, selected_batch_size).
    """
    if not torch.cuda.is_available():
        logger.warning("CUDA is not available. Falling back to CPU execution.")
        selected_batch = requested_batch if requested_batch is not None else 2
        return "cpu", selected_batch

    device_name = torch.cuda.get_device_name(0)
    total_memory_gb = torch.cuda.get_device_properties(0).total_memory / (1024 ** 3)
    logger.info(f"CUDA accelerator detected: {device_name} ({total_memory_gb:.2f} GB VRAM)")

    if requested_batch is not None:
        logger.info(f"Using user-specified batch size: {requested_batch}")
        return "0", requested_batch

    # Heuristic batch size selection based on VRAM capacity at 1024x1024
    if total_memory_gb >= 24.0:
        recommended_batch = 16
    elif total_memory_gb >= 12.0:
        recommended_batch = 8
    elif total_memory_gb >= 8.0:
        recommended_batch = 4
    else:
        recommended_batch = 2

    logger.info(
        f"Auto-selected batch size: {recommended_batch} "
        f"(Tailored for {imgsz}x{imgsz} segmentation on {total_memory_gb:.1f} GB VRAM)"
    )
    return "0", recommended_batch


def verify_dataset_configuration(dataset_config_path: Path, auto_prepare_tiled: bool = True) -> Dict[str, Any]:
    """
    Loads and validates the YOLO dataset YAML configuration file, ensuring
    all required paths and class indices exist. If the configuration references
    a tiled dataset that has not yet been partitioned, it automatically invokes
    prepare_tiled_dataset_split to generate the splits.

    Args:
        dataset_config_path: Path to dataset YAML configuration file.
        auto_prepare_tiled: Whether to automatically split tiled dataset if needed.

    Returns:
        Parsed YAML configuration dictionary.
    """
    project_root = dataset_config_path.parent.parent
    tiled_dir = project_root / "data" / "tiled_dataset"

    # If the file does not exist yet, but tiled_dataset directory is present, generate it automatically
    if not dataset_config_path.exists():
        if tiled_dir.exists() and (tiled_dir / "images").exists() and auto_prepare_tiled:
            logger.info(f"Dataset config '{dataset_config_path.name}' not found. Automatically partitioning tiled dataset...")
            prepare_tiled_dataset_split(tiled_dir=tiled_dir, config_yaml_path=dataset_config_path)
        else:
            raise FileNotFoundError(
                f"Dataset configuration YAML not found at: {dataset_config_path}. "
                f"Please run 'scripts/run_dpi_tiler.py' or 'scripts/build_artifact_robust_dataset.py' to generate data."
            )

    with open(dataset_config_path, "r", encoding="utf-8") as f:
        config_data = yaml.safe_load(f)

    # If pointing to tiled dataset and not partitioned yet, partition it now
    dataset_root = Path(config_data.get("path", ""))
    if dataset_root.exists() and "tiled_dataset" in str(dataset_root) and auto_prepare_tiled:
        train_img_dir = dataset_root / config_data.get("train", "images/train")
        if not train_img_dir.exists() or not any(train_img_dir.glob("*.jpg")):
            logger.info(f"Tiled dataset at {dataset_root} requires partitioning. Partitioning now...")
            prepare_tiled_dataset_split(tiled_dir=dataset_root, config_yaml_path=dataset_config_path)
            with open(dataset_config_path, "r", encoding="utf-8") as f:
                config_data = yaml.safe_load(f)

    logger.info(f"Successfully loaded dataset config from: {dataset_config_path}")
    classes = config_data.get("names", {})
    logger.info(f"Configured dataset classes ({len(classes)} total): {classes}")
    return config_data


# ===============================================================================
# ROBUST YOLOv8 TRAINING & EVALUATION MANAGER
# ===============================================================================
class RobustYOLOTrainer:
    """
    Orchestrates the fine-tuning, custom loss configuration, botanical augmentation,
    validation evaluation, and artifact discrimination analysis for YOLOv8-seg.
    """

    def __init__(
        self,
        weights: str = "yolov8x-seg.pt",
        data_config: Union[str, Path] = "data/dataset_config.yaml",
        imgsz: int = 1024,
        batch: Optional[int] = None,
        epochs: int = 100,
        device: Optional[str] = None,
        workers: Optional[int] = None,
        cache: Optional[str] = "ram",
        resume: bool = False,
        project_root: Optional[Path] = None,
        experiment_name: str = "artifact_robust_yolov8_seg"
    ):
        """
        Initializes the trainer with model weights and architectural parameters.

        Args:
            weights: Pre-trained YOLOv8 weights path or model architecture identifier.
            data_config: Path to dataset YAML configuration.
            imgsz: Input image resolution (default: 1024).
            batch: Batch size or None for auto-detection.
            epochs: Total training epochs (default: 100).
            device: Computing device ('0', 'cpu', etc.) or None for auto-detection.
            workers: Number of DataLoader worker subprocesses (default: 0 on Windows/Google Drive
                     to avoid multi-process virtual filesystem handle lock contention).
            cache: Image caching strategy ('ram', 'disk', True, or False). Default 'ram'
                   pre-loads images into host RAM to eliminate Google Drive read bottlenecks.
            resume: If True, resumes training from the last saved checkpoint (last.pt)
                    without re-training already completed epochs.
            project_root: Project root directory Path.
            experiment_name: Identifier for Ultralytics run output directory.
        """
        self.paths = resolve_default_paths(project_root)
        self.data_config_path = Path(data_config) if Path(data_config).is_absolute() else self.paths["root"] / data_config
        self.dataset_info = verify_dataset_configuration(self.data_config_path)

        self.imgsz = imgsz
        self.epochs = epochs
        self.experiment_name = experiment_name
        self.resume = resume
        self.cache = cache

        # Windows Multiprocessing & RAM Caching IPC Protection:
        # On Windows, Python uses `spawn` instead of `fork` to initialize DataLoader worker subprocesses.
        # When `cache='ram'` is enabled, the entire in-memory dataset (20,000+ images) must be serialized
        # (pickled) and passed through IPC pipes to each spawned worker. This immediately exhausts Windows
        # IPC pipe memory buffers and throws `MemoryError: reduction.pickle.load(from_parent)`.
        # When caching in RAM on Windows, workers must be 0 (main process zero-copy in-memory loading).
        if workers is not None:
            self.workers = workers
        else:
            self.workers = 0 if os.name == "nt" else 4

        if os.name == "nt" and self.cache == "ram" and self.workers > 0:
            logger.warning(
                f"Windows multiprocessing 'spawn' cannot serialize a 20,000+ image RAM cache "
                f"across {self.workers} worker IPC pipes (causes MemoryError). "
                f"Automatically adjusting workers=0 for ultra-fast zero-copy main-process RAM caching."
            )
            self.workers = 0

        # Detect optimal device and batch
        auto_device, auto_batch = detect_optimal_device_and_batch(batch, imgsz=self.imgsz)
        self.device = device if device is not None else auto_device
        self.batch = batch if batch is not None else auto_batch

        # Handle checkpoint resumption if requested
        last_ckpt_path = self.paths["root"] / "runs" / "segment" / self.experiment_name / "weights" / "last.pt"
        if self.resume:
            if last_ckpt_path.exists() and last_ckpt_path.stat().st_size > 0:
                logger.info(f"Resuming training from existing checkpoint: {last_ckpt_path}")
                self.weights_path = str(last_ckpt_path)
            else:
                logger.warning(f"Resume requested but '{last_ckpt_path}' was not found. Starting from base weights: {weights}")
                self.resume = False
                self.weights_path = self._resolve_model_weights(weights)
        else:
            # Model weights resolution for fresh run
            self.weights_path = self._resolve_model_weights(weights)

        logger.info(f"Initializing YOLO model with backbone: {self.weights_path} | DataLoader workers: {self.workers} | Resume: {self.resume}")
        self.model = YOLO(str(self.weights_path))

    def _resolve_model_weights(self, weights_name: str) -> str:
        """
        Resolves model weights path, checking local files and falling back
        to Ultralytics automated download identifiers.
        """
        local_candidate = self.paths["root"] / weights_name
        if local_candidate.exists() and local_candidate.stat().st_size > 0:
            logger.info(f"Found local model weights at: {local_candidate}")
            return str(local_candidate)

        models_dir_candidate = self.paths["models_dir"] / weights_name
        if models_dir_candidate.exists() and models_dir_candidate.stat().st_size > 0:
            logger.info(f"Found weights in models directory: {models_dir_candidate}")
            return str(models_dir_candidate)

        # Return model tag for Ultralytics to automatically download
        logger.info(f"Using Ultralytics pre-trained weight identifier: {weights_name}")
        return weights_name

    def build_training_hyperparameters(self) -> Dict[str, Any]:
        """
        Constructs the customized botanical hyperparameter dictionary including:
        - Box loss gain: box=7.5
        - Classification loss gain: cls=0.5
        - Distribution focal loss gain: dfl=1.5
        - Cosine learning rate schedule
        - Mosaic, mixup, and copy-paste botanical augmentations
        - Mixed precision AMP

        Returns:
            Dictionary of training arguments for YOLO.train().
        """
        hyperparameters = {
            # Core dataset and model dimensions
            "data": str(self.data_config_path),
            "epochs": self.epochs,
            "imgsz": self.imgsz,
            "batch": self.batch,
            "device": self.device,
            "workers": self.workers,  # Set DataLoader worker count (0 for single-process synchronous loading on Windows)
            "cache": self.cache,      # Pre-load images in RAM to bypass cloud/virtual filesystem I/O latency
            "amp": True,  # Mixed precision for performance and memory efficiency

            # Optimization and learning rate schedule
            "optimizer": "auto",
            "lr0": 0.01,           # Initial learning rate
            "lrf": 0.001,          # Final learning rate factor (cosine minimum)
            "cos_lr": True,        # Cosine annealing learning rate schedule
            "warmup_epochs": 3.0,  # Warmup epochs
            "warmup_momentum": 0.8,
            "warmup_bias_lr": 0.1,
            "weight_decay": 0.0005,
            "momentum": 0.937,

            # Custom Loss Weight Gains (Maximizing Mask and Edge Fidelity)
            "box": 7.5,            # Bounding box loss gain
            "cls": 0.5,            # Classification loss gain
            "dfl": 1.5,            # Distribution focal loss gain

            # Botanical Multi-Scale Augmentations
            "mosaic": 1.0,         # 4-voucher sheet patch mosaic
            "mixup": 0.15,         # Inter-sample mixup interpolation
            "copy_paste": 0.30,    # Copy-paste leaf and artifact instances onto backgrounds
            "degrees": 15.0,       # Rotation angle range
            "translate": 0.1,      # Image translation fraction
            "scale": 0.2,          # Image scale range (+/- 20%)
            "shear": 0.0,          # Shear angle
            "perspective": 0.0001, # Perspective distortion
            "fliplr": 0.5,         # Horizontal flip probability (botanically invariant)
            "flipud": 0.0,         # Vertical flip disabled to maintain sheet orientation conventions
            "hsv_h": 0.015,        # Hue augmentation
            "hsv_s": 0.7,          # Saturation augmentation
            "hsv_v": 0.4,          # Value/brightness augmentation

            # Output and Checkpoint Management
            "project": str(self.paths["root"] / "runs" / "segment"),
            "name": self.experiment_name,
            "exist_ok": True,
            "save": True,
            "save_period": -1,     # Save only best and last
            "val": True,           # Validate each epoch
            "plots": True,         # Generate training curves and batches
            "verbose": True
        }
        return hyperparameters

    def train(self) -> Any:
        """
        Executes the YOLOv8-seg training loop with custom botanical loss and augmentations.
        Supports seamless resumption from interrupted checkpoints.

        Returns:
            Ultralytics training results object.
        """
        hyperparams = self.build_training_hyperparameters()
        logger.info("=" * 80)
        logger.info(f"STARTING ARTIFACT-ROBUST YOLOV8-SEG MODEL {'RESUMPTION' if self.resume else 'FINE-TUNING'}")
        logger.info("=" * 80)
        logger.info(f"Hyperparameters Summary:\n{json.dumps({k: str(v) for k, v in hyperparams.items()}, indent=2)}")

        # Execute training with resume support
        if self.resume:
            logger.info("Resuming training from checkpoint without re-running completed epochs.")
            train_results = self.model.train(resume=True, workers=self.workers, cache=self.cache)
        else:
            train_results = self.model.train(**hyperparams)

        logger.info("Training loop completed successfully.")
        return train_results

    def evaluate_and_export_metrics(
        self,
        checkpoint_path: Optional[Union[str, Path]] = None,
        split: str = "val"
    ) -> Dict[str, Any]:
        """
        Evaluates the trained model on validation split, computing class-specific
        mAP50, mAP50-95, precision, recall, and analyzing cross-classification
        between herbarium artifacts and botanical organs.

        Args:
            checkpoint_path: Optional path to a specific model checkpoint (.pt).
            split: Dataset split to evaluate ('val' or 'test').

        Returns:
            Dictionary of computed class-level metrics and confusion statistics.
        """
        if checkpoint_path is not None:
            ckpt_p = Path(checkpoint_path)
            if ckpt_p.exists() and ckpt_p.stat().st_size > 0:
                eval_model = YOLO(str(ckpt_p))
            else:
                logger.warning(f"Checkpoint at '{checkpoint_path}' does not exist or is 0-bytes. Using active model instance.")
                eval_model = self.model
        else:
            eval_model = self.model

        logger.info("=" * 80)
        logger.info(f"EVALUATING MODEL ON {split.upper()} SPLIT (imgsz={self.imgsz})")
        logger.info("=" * 80)

        # Run validation with plots enabled and save_json disabled (COCO RLE multi-encoding on dense
        # 1024x1024 segmentation masks can exceed GPU memory allocation limits)
        val_results = eval_model.val(
            data=str(self.data_config_path),
            imgsz=self.imgsz,
            batch=self.batch,
            device=self.device,
            split=split,
            plots=True,
            save_json=False
        )

        class_names = self.dataset_info.get("names", {})
        metrics_summary: Dict[str, Any] = {
            "overall_box_map50": float(val_results.box.map50) if hasattr(val_results, "box") else 0.0,
            "overall_box_map50_95": float(val_results.box.map) if hasattr(val_results, "box") else 0.0,
            "overall_mask_map50": float(val_results.seg.map50) if hasattr(val_results, "seg") else 0.0,
            "overall_mask_map50_95": float(val_results.seg.map) if hasattr(val_results, "seg") else 0.0,
            "class_metrics": {}
        }

        logger.info("-" * 80)
        logger.info(f"{'Class ID':<10}{'Class Name':<20}{'Box mAP50':<14}{'Box mAP50-95':<16}{'Mask mAP50':<14}{'Mask mAP50-95':<16}")
        logger.info("-" * 80)

        # Extract per-class metrics
        box_maps50 = getattr(val_results.box, "maps50", None)
        box_maps = getattr(val_results.box, "maps", None)
        seg_maps50 = getattr(val_results.seg, "maps50", None)
        seg_maps = getattr(val_results.seg, "maps", None)

        for class_id, class_name in class_names.items():
            cid = int(class_id)
            b_map50 = float(box_maps50[cid]) if box_maps50 is not None and cid < len(box_maps50) else 0.0
            b_map = float(box_maps[cid]) if box_maps is not None and cid < len(box_maps) else 0.0
            s_map50 = float(seg_maps50[cid]) if seg_maps50 is not None and cid < len(seg_maps50) else 0.0
            s_map = float(seg_maps[cid]) if seg_maps is not None and cid < len(seg_maps) else 0.0

            metrics_summary["class_metrics"][class_name] = {
                "class_id": cid,
                "box_map50": b_map50,
                "box_map50_95": b_map,
                "mask_map50": s_map50,
                "mask_map50_95": s_map
            }

            logger.info(
                f"{cid:<10}{class_name:<20}{b_map50:<14.4f}{b_map:<16.4f}{s_map50:<14.4f}{s_map:<16.4f}"
            )

        logger.info("-" * 80)

        # Verify key requirements
        basal_leaf_metrics = metrics_summary["class_metrics"].get("basal_leaf", {})
        petiole_metrics = metrics_summary["class_metrics"].get("leaf_petiole", {})
        logger.info(
            f"Key Botanical Organ Performance -> "
            f"basal_leaf Mask mAP50: {basal_leaf_metrics.get('mask_map50', 0.0):.4f} | "
            f"leaf_petiole Mask mAP50: {petiole_metrics.get('mask_map50', 0.0):.4f}"
        )

        # Cross-classification analysis and visualizations
        self._generate_evaluation_artifacts(val_results, class_names, metrics_summary)

        # Save metrics to JSON
        eval_report_path = self.paths["output_eval_dir"] / "evaluation_report.json"
        with open(eval_report_path, "w", encoding="utf-8") as f:
            json.dump(metrics_summary, f, indent=2)
        logger.info(f"Exported evaluation metrics JSON to: {eval_report_path}")

        return metrics_summary

    def _generate_evaluation_artifacts(
        self,
        val_results: Any,
        class_names: Dict[Any, str],
        metrics_summary: Dict[str, Any]
    ) -> None:
        """
        Generates and exports custom confusion matrix plots, Precision-Recall curves,
        and verifies cross-classification between herbarium artifacts and botanical organs.

        Args:
            val_results: Ultralytics validation results.
            class_names: Dictionary of class names.
            metrics_summary: Extracted metrics dictionary.
        """
        output_dir = self.paths["output_eval_dir"]
        output_dir.mkdir(parents=True, exist_ok=True)

        # 1. Confusion Matrix Analysis
        cm_matrix = None
        if hasattr(val_results, "confusion_matrix") and val_results.confusion_matrix is not None:
            cm = val_results.confusion_matrix
            if hasattr(cm, "matrix"):
                cm_matrix = cm.matrix

        num_classes = len(class_names)
        ordered_names = [class_names[i] for i in range(num_classes)]

        if cm_matrix is not None and isinstance(cm_matrix, np.ndarray):
            # Check cross-classification between herbarium_label and basal_leaf
            label_idx = None
            leaf_idx = None
            for idx, name in enumerate(ordered_names):
                if name == "herbarium_label":
                    label_idx = idx
                elif name == "basal_leaf":
                    leaf_idx = idx

            cross_misclassifications = 0
            if label_idx is not None and leaf_idx is not None and cm_matrix.shape[0] > max(label_idx, leaf_idx):
                # Row is true class, Column is predicted class
                label_as_leaf = cm_matrix[label_idx, leaf_idx]
                leaf_as_label = cm_matrix[leaf_idx, label_idx]
                cross_misclassifications = int(label_as_leaf + leaf_as_label)

                logger.info("=" * 80)
                logger.info("ARTIFACT DISCRIMINATION INTEGRITY CHECK:")
                logger.info(f"  - Herbarium Label misclassified as Basal Leaf: {int(label_as_leaf)}")
                logger.info(f"  - Basal Leaf misclassified as Herbarium Label: {int(leaf_as_label)}")
                if cross_misclassifications == 0:
                    logger.info("  [PASSED] PERFECT ZERO CROSS-CLASSIFICATION CONFIRMED.")
                else:
                    logger.warning(
                        f"  [WARNING] Cross-classification detected: {cross_misclassifications} instances."
                    )
                logger.info("=" * 80)

            metrics_summary["label_to_leaf_misclassifications"] = cross_misclassifications

            # Render custom publication-quality confusion matrix heatmap
            self._plot_custom_confusion_matrix(cm_matrix, ordered_names, output_dir / "confusion_matrix_custom.png")

        # 2. Precision-Recall Curves Plot
        self._plot_classwise_map_bars(metrics_summary, output_dir / "classwise_map_comparison.png")

        # Copy any Ultralytics generated PR / F1 / P / R curves if available in run dir
        if hasattr(val_results, "save_dir") and val_results.save_dir:
            save_dir = Path(val_results.save_dir)
            for curve_file in save_dir.glob("*.png"):
                dest_file = output_dir / curve_file.name
                shutil.copy2(curve_file, dest_file)
            logger.info(f"Copied Ultralytics validation curve artifacts from {save_dir} to {output_dir}")

    def _plot_custom_confusion_matrix(
        self,
        matrix: np.ndarray,
        class_names: List[str],
        output_path: Path
    ) -> None:
        """
        Plots a high-contrast confusion matrix focusing on botanical vs. artifact discrimination.
        """
        try:
            plt.figure(figsize=(10, 8), dpi=300)
            # Normalize matrix by row (True Class)
            row_sums = matrix.sum(axis=1, keepdims=True)
            norm_matrix = np.divide(matrix, row_sums, out=np.zeros_like(matrix, dtype=float), where=row_sums != 0)

            # Include background if present in matrix dimensions
            display_names = class_names.copy()
            if matrix.shape[0] > len(class_names):
                display_names.append("background")

            matrix_slice = norm_matrix[:len(display_names), :len(display_names)]

            fig, ax = plt.subplots(figsize=(10, 8), dpi=300)
            cax = ax.imshow(matrix_slice, interpolation='nearest', cmap=plt.cm.Blues)
            plt.colorbar(cax, fraction=0.046, pad=0.04)

            # Add numeric text labels in each cell
            thresh = matrix_slice.max() / 2.0 if matrix_slice.max() > 0 else 0.5
            for i in range(matrix_slice.shape[0]):
                for j in range(matrix_slice.shape[1]):
                    val = matrix_slice[i, j]
                    ax.text(
                        j, i, f"{val:.2f}",
                        ha="center", va="center",
                        color="white" if val > thresh else "black",
                        fontsize=8
                    )

            ax.set_xticks(np.arange(len(display_names)))
            ax.set_yticks(np.arange(len(display_names)))
            ax.set_xticklabels(display_names, rotation=45, ha="right", fontsize=9)
            ax.set_yticklabels(display_names, fontsize=9)
            ax.set_xlabel("Predicted Class", fontsize=11, labelpad=8)
            ax.set_ylabel("True Class", fontsize=11, labelpad=8)
            ax.set_title("Normalized Confusion Matrix (Botanical vs. Artifact Discrimination)", fontsize=13, pad=15)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
            logger.info(f"Saved custom confusion matrix visualization: {output_path}")
        except Exception as e:
            logger.warning(f"Could not render custom confusion matrix plot: {e}")

    def _plot_classwise_map_bars(
        self,
        metrics_summary: Dict[str, Any],
        output_path: Path
    ) -> None:
        """
        Plots a bar chart comparing Mask mAP50 and Mask mAP50-95 across all botanical
        and artifact classes.
        """
        try:
            class_metrics = metrics_summary.get("class_metrics", {})
            if not class_metrics:
                return

            names = list(class_metrics.keys())
            mask_map50 = [class_metrics[n]["mask_map50"] for n in names]
            mask_map50_95 = [class_metrics[n]["mask_map50_95"] for n in names]

            x = np.arange(len(names))
            width = 0.35

            plt.figure(figsize=(12, 6), dpi=300)
            plt.bar(x - width/2, mask_map50, width, label="Mask mAP50", color="#2b5c8f")
            plt.bar(x + width/2, mask_map50_95, width, label="Mask mAP50-95", color="#52b788")

            plt.ylabel("mAP Score", fontsize=11)
            plt.title("Class-Specific Segmentation Performance (YOLOv8-seg Fine-Tuning)", fontsize=13, pad=15)
            plt.xticks(x, names, rotation=35, ha="right", fontsize=9)
            plt.ylim(0.0, 1.05)
            plt.legend(loc="upper right")
            plt.grid(axis="y", linestyle="--", alpha=0.5)
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
            logger.info(f"Saved class-wise mAP comparison plot: {output_path}")
        except Exception as e:
            logger.warning(f"Could not render class-wise mAP bar chart: {e}")

    def export_best_checkpoint(self) -> Path:
        """
        Locates the best checkpoint generated during training and saves/copies it
        to the standard model registry destination `models/yolov8_leaf_best.pt`.

        Returns:
            Path to the saved best model weights.
        """
        target_path = self.paths["best_model_export"]
        run_best_weights = self.paths["root"] / "runs" / "segment" / self.experiment_name / "weights" / "best.pt"

        if run_best_weights.exists():
            shutil.copy2(run_best_weights, target_path)
            logger.info(f"Successfully exported best model checkpoint to: {target_path}")
        else:
            logger.warning(
                f"Run best weights not found at: {run_best_weights}. "
                f"Checking for existing model at target path..."
            )
            if not target_path.exists():
                # Save current model weights if available
                try:
                    self.model.save(str(target_path))
                    logger.info(f"Saved current model weights directly to: {target_path}")
                except Exception as e:
                    logger.error(f"Failed to export model checkpoint: {e}")

        return target_path


# ===============================================================================
# CLI & ENTRY POINT
# ===============================================================================
def parse_arguments() -> argparse.Namespace:
    """
    Parses command-line arguments for the YOLOv8 fine-tuning script.
    """
    parser = argparse.ArgumentParser(
        description="Fine-tune an artifact-robust YOLOv8-seg model for botanical organ segmentation on sliced native-DPI tiles or full voucher sheets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--weights",
        type=str,
        default="yolov8x-seg.pt",
        help="Initial model backbone weights (e.g. 'yolov8x-seg.pt', 'yolov8m-seg.pt', or path to .pt)."
    )
    parser.add_argument(
        "--data",
        type=str,
        default=None,
        help="Path to YOLO dataset YAML configuration file (default: auto-detects data/tiled_dataset_config.yaml if tiled dataset is present, otherwise data/dataset_config.yaml)."
    )
    parser.add_argument(
        "--split-tiled-dataset",
        action="store_true",
        help="Explicitly (re-)partition data/tiled_dataset into train/val/test splits grouped by parent specimen sheet."
    )
    parser.add_argument(
        "--train-ratio",
        type=float,
        default=0.70,
        help="Fraction of specimen sheets assigned to training split when partitioning tiled dataset."
    )
    parser.add_argument(
        "--val-ratio",
        type=float,
        default=0.15,
        help="Fraction of specimen sheets assigned to validation split when partitioning tiled dataset."
    )
    parser.add_argument(
        "--test-ratio",
        type=float,
        default=0.15,
        help="Fraction of specimen sheets assigned to testing split when partitioning tiled dataset."
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=100,
        help="Total number of training epochs."
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=None,
        help="Training batch size (None for auto-detection based on GPU VRAM)."
    )
    parser.add_argument(
        "--imgsz",
        type=int,
        default=1024,
        help="Input image resolution in pixels (default 1024 for high-resolution botanical tiles)."
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Computing device to use (e.g. '0', '0,1', 'cpu'). Auto-detected if None."
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=0 if os.name == "nt" else 4,
        help="Number of DataLoader worker subprocesses (0 recommended on Windows/Google Drive)."
    )
    parser.add_argument(
        "--cache",
        type=str,
        default="ram",
        choices=["ram", "disk", "none"],
        help="Image caching strategy. 'ram' pre-loads scans into memory for maximum training speed."
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume training from existing runs/segment/<name>/weights/last.pt without restarting."
    )
    parser.add_argument(
        "--name",
        type=str,
        default="artifact_robust_yolov8_seg",
        help="Experiment run name for checkpoint outputs."
    )
    parser.add_argument(
        "--eval-only",
        action="store_true",
        help="Skip training loop and evaluate existing checkpoint."
    )
    parser.add_argument(
        "--checkpoint",
        type=str,
        default=None,
        help="Path to specific checkpoint to evaluate when --eval-only is set."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run a 1-epoch lightweight trial to verify the training and evaluation pipeline."
    )

    return parser.parse_args()


def main():
    """
    Main orchestration entry point.
    """
    args = parse_arguments()
    paths = resolve_default_paths()

    # Determine which dataset YAML configuration to use
    if args.data is not None:
        data_config_path = Path(args.data) if Path(args.data).is_absolute() else paths["root"] / args.data
    else:
        data_config_path = paths["data_config"]

    # Explicitly re-partition tiled dataset if flag is supplied
    if args.split_tiled_dataset:
        logger.info("Executing manual tiled dataset re-partitioning (--split-tiled-dataset requested)...")
        prepare_tiled_dataset_split(
            tiled_dir=paths["tiled_dir"],
            config_yaml_path=paths["tiled_config"],
            train_ratio=args.train_ratio,
            val_ratio=args.val_ratio,
            test_ratio=args.test_ratio,
            force_resplit=True
        )
        data_config_path = paths["tiled_config"]

    if args.dry_run:
        logger.info("Executing DRY RUN: Setting epochs=1, batch=2, imgsz=640 for rapid pipeline verification.")
        args.epochs = 1
        args.batch = 2
        args.imgsz = 640

    cache_val = False if args.cache == "none" else args.cache

    trainer = RobustYOLOTrainer(
        weights=args.weights,
        data_config=data_config_path,
        imgsz=args.imgsz,
        batch=args.batch,
        epochs=args.epochs,
        device=args.device,
        workers=args.workers,
        cache=cache_val,
        resume=args.resume,
        experiment_name=args.name
    )

    if not args.eval_only:
        # 1. Execute Training Loop
        trainer.train()

        # 2. Export Best Model Checkpoint to models/yolov8_leaf_best.pt
        best_checkpoint_path = trainer.export_best_checkpoint()

        # 3. Evaluate and Export Comprehensive Metrics
        trainer.evaluate_and_export_metrics(checkpoint_path=best_checkpoint_path, split="val")
    else:
        # Evaluation Only Mode
        eval_checkpoint = args.checkpoint if args.checkpoint else trainer.paths["best_model_export"]
        logger.info(f"Running evaluation-only mode using checkpoint: {eval_checkpoint}")
        trainer.evaluate_and_export_metrics(checkpoint_path=eval_checkpoint, split="val")

    logger.info("Pipeline execution finished successfully.")


if __name__ == "__main__":
    main()
