import os
import shutil
import random
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any

from .config import STANDARD_CLASS_MAPPING

logger = logging.getLogger("ArtifactRobustYOLO")

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
    
    Also writes and validates the Ultralytics dataset configuration YAML.

    Args:
        tiled_dir (Path): Root directory containing tiled dataset ('images' and 'labels').
        config_yaml_path (Path, optional): Destination path for config yaml.
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

        # Collect all image files
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
                    # Create empty label file for hard negative background tiles
                    dest_lbl_path.touch()

                split_tile_counts[split] += 1

        logger.info(
            f"Successfully partitioned tiled dataset -> "
            f"Train: {split_tile_counts['train']} tiles ({len(train_specs)} specimens) | "
            f"Val: {split_tile_counts['val']} tiles ({len(val_specs)} specimens) | "
            f"Test: {split_tile_counts['test']} tiles ({len(test_specs)} specimens)"
        )

    # Write Ultralytics dataset configuration YAML
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
