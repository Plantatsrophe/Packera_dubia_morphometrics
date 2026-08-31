#!/usr/bin/env python3
"""
===============================================================================
Script: prepare_lm2_dataset.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Dual-purpose Pre-Processor and Dataset Generator for LeafMachine2 (LM2).
    
    Key Capabilities:
      1. SAM 2 Binary Mask Ingestion & COCO 1.0 JSON Dataset Generation:
         - Ingests high-resolution binary instance segmentation masks from SAM 2.
         - Uses OpenCV (cv2.findContours) to extract vector polygon boundaries.
         - Converts contours into standardized COCO 1.0 polygon annotations.
         - Implements a strict class mapping dictionary to align with LeafMachine2's
           native Plant Component Detector (PCD) taxonomy:
             * 'basal_leaf_whole'   -> 'ideal_leaf' (Category ID 1)
             * 'basal_leaf_partial' -> 'partial_leaf' (Category ID 2)
         - Intentionally OMITS 'cauline_leaf' (and non-target structures) from the
           final JSON so that the Detectron2 loss function treats cauline leaves
           entirely as background noise during PCD fine-tuning/training.
         - Supports optional train/val splitting (--val-split) and background negative mining.

      2. Herbarium Sheet Image Staging & Symlink Management:
         - Creates the standardized LM2_Project directory structure (Data/images, Data/output, Data/annotations).
         - Sets up safe relative or absolute symlinks to raw voucher images.
         - Generates an asset audit manifest CSV.

Usage Examples:
    # 1. Full Pipeline: Stage images and convert SAM 2 masks to COCO JSON:
    python scripts/data_prep/prepare_lm2_dataset.py \\
        --masks-dir data/raw_annotations/masks \\
        --coco-output LM2_Project/Data/annotations/coco_pcd_packera.json

    # 2. Only Convert SAM 2 Masks to COCO JSON (with 80/20 train/val split):
    python scripts/data_prep/prepare_lm2_dataset.py \\
        --skip-symlinks \\
        --masks-dir data/raw_annotations/masks \\
        --coco-output data/raw_annotations/coco_pcd_packera.json \\
        --val-split 0.20

    # 3. Only Stage Images and Symlinks:
    python scripts/data_prep/prepare_lm2_dataset.py \\
        --skip-coco \\
        --input-dirs data/raw_vouchers \\
        --lm2-root LM2_Project
===============================================================================
"""

import argparse
import csv
import json
import logging
import os
import random
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

try:
    import cv2
    import numpy as np
except ImportError:
    # Graceful fallback if invoked outside active virtualenv
    cv2 = None
    np = None

# ===============================================================================
# TAXONOMIC CLASS MAPPING & COCO SCHEMA CONSTANTS
# ===============================================================================

# Strict class mapping from SAM 2 botanical annotations to LeafMachine2 PCD taxonomy
PCD_CLASS_MAPPING: Dict[str, str] = {
    "basal_leaf_whole": "ideal_leaf",
    "basal_leaf_partial": "partial_leaf",
}

# Categories intentionally omitted so Detectron2 treats them as background noise
OMITTED_CLASSES: Set[str] = {
    "cauline_leaf",
    "cauline_stem",
    "root_rhizome",
    "basal_rosette_clump",
    "capitulum",
}

# Standardized COCO 1.0 categories for LeafMachine2 Plant Component Detector
PCD_COCO_CATEGORIES: List[Dict[str, Any]] = [
    {
        "id": 1,
        "name": "ideal_leaf",
        "supercategory": "plant_component",
    },
    {
        "id": 2,
        "name": "partial_leaf",
        "supercategory": "plant_component",
    },
]

PCD_CATEGORY_ID_MAP: Dict[str, int] = {
    "ideal_leaf": 1,
    "partial_leaf": 2,
}

# Supported raster image extensions for herbarium processing
VALID_IMAGE_EXTENSIONS: Set[str] = {
    ".jpg",
    ".jpeg",
    ".png",
    ".tif",
    ".tiff",
    ".webp",
    ".bmp",
}


def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure structured logging output."""
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    return logging.getLogger("LM2_DatasetPrep")


def get_project_root() -> Path:
    """Dynamically resolves the workspace root directory."""
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "data").exists() or (parent / "LM2_Project").exists() or (parent / ".git").exists():
            return parent
    return current.parents[2] if len(current.parents) > 2 else current.parent


def verify_image_file(file_path: Path) -> bool:
    """
    Perform a lightweight header verification using PIL to ensure the file
    is a readable, valid image format.
    """
    try:
        from PIL import Image
        with Image.open(file_path) as img:
            img.verify()
        return True
    except Exception:
        return False


def create_lm2_directories(
    lm2_root: Path,
    images_dir_name: str = "images",
    output_dir_name: str = "output",
    configs_dir_name: str = "configs",
    annotations_dir_name: str = "annotations",
    dry_run: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Tuple[Path, Path, Path, Path]:
    """
    Create the standardized LM2_Project directory structure.
    
    Structure:
      <lm2_root>/
        Data/
          <images_dir_name>/       <- Symlinks to raw herbarium images
          <output_dir_name>/       <- Destination for LM2 analysis output
          <annotations_dir_name>/  <- COCO 1.0 JSON datasets for PCD training
        <configs_dir_name>/        <- Project-specific configuration files
    """
    data_dir = lm2_root / "Data"
    images_dir = data_dir / images_dir_name
    output_dir = data_dir / output_dir_name
    annotations_dir = data_dir / annotations_dir_name
    configs_dir = lm2_root / configs_dir_name

    dirs_to_create = [lm2_root, data_dir, images_dir, output_dir, annotations_dir, configs_dir]

    for d in dirs_to_create:
        if not dry_run:
            d.mkdir(parents=True, exist_ok=True)
        if logger:
            logger.debug(f"Ensured directory exists: {d}")

    return images_dir, output_dir, annotations_dir, configs_dir


def collect_raw_images(
    input_dirs: List[Path],
    verify_images: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[Path], List[Tuple[Path, str]]]:
    """
    Scan provided input directories and collect valid herbarium sheet images.
    Returns:
      (valid_images, filtered_out_files)
    """
    valid_images: List[Path] = []
    filtered_out: List[Tuple[Path, str]] = []

    for input_dir in input_dirs:
        if not input_dir.exists():
            if logger:
                logger.warning(f"Input directory does not exist: {input_dir}")
            continue

        if logger:
            logger.info(f"Scanning image directory: {input_dir.resolve()}")

        for entry in input_dir.iterdir():
            if not entry.is_file():
                continue

            if entry.name.startswith("."):
                filtered_out.append((entry, "Hidden file"))
                continue

            try:
                if entry.stat().st_size == 0:
                    filtered_out.append((entry, "Zero-byte file"))
                    continue
            except OSError as e:
                filtered_out.append((entry, f"Stat error: {e}"))
                continue

            ext = entry.suffix.lower()
            if ext not in VALID_IMAGE_EXTENSIONS:
                filtered_out.append((entry, f"Unsupported extension: '{ext}'"))
                continue

            if verify_images:
                if not verify_image_file(entry):
                    filtered_out.append((entry, "Corrupt or unreadable image header"))
                    continue

            valid_images.append(entry.resolve())

    return valid_images, filtered_out


def create_symlinks(
    source_images: List[Path],
    target_dir: Path,
    relative: bool = False,
    overwrite: bool = False,
    limit: Optional[int] = None,
    dry_run: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[Dict[str, str]], int, int, int]:
    """
    Safely create symbolic links pointing from target_dir to source_images.

    Returns:
      (manifest_records, created_count, skipped_count, error_count)
    """
    manifest: List[Dict[str, str]] = []
    created_count = 0
    skipped_count = 0
    error_count = 0

    if limit and limit > 0:
        source_images = source_images[:limit]
        if logger:
            logger.info(f"Limiting symlinks to first {limit} images.")

    seen_names: Set[str] = set()

    for src_path in source_images:
        dest_filename = src_path.name

        if dest_filename in seen_names:
            if logger:
                logger.warning(f"Filename collision detected for '{dest_filename}'. Skipping duplicate.")
            skipped_count += 1
            continue

        seen_names.add(dest_filename)
        dest_link_path = target_dir / dest_filename

        if relative:
            try:
                link_target = os.path.relpath(src_path, target_dir)
            except ValueError:
                link_target = str(src_path)
        else:
            link_target = str(src_path)

        link_status = "PENDING"
        
        if dest_link_path.is_symlink() or dest_link_path.exists():
            if overwrite:
                if not dry_run:
                    try:
                        dest_link_path.unlink()
                        os.symlink(link_target, dest_link_path)
                        link_status = "OVERWRITTEN"
                        created_count += 1
                    except OSError as err:
                        if logger:
                            logger.error(f"Failed to overwrite symlink {dest_link_path}: {err}")
                        link_status = f"ERROR: {err}"
                        error_count += 1
                else:
                    link_status = "DRY_RUN_OVERWRITE"
                    created_count += 1
            else:
                link_status = "SKIPPED_EXISTS"
                skipped_count += 1
        else:
            if not dry_run:
                try:
                    os.symlink(link_target, dest_link_path)
                    link_status = "CREATED"
                    created_count += 1
                except OSError as err:
                    if logger:
                        logger.error(f"Failed to create symlink {dest_link_path} -> {link_target}: {err}")
                    link_status = f"ERROR: {err}"
                    error_count += 1
            else:
                link_status = "DRY_RUN_CREATED"
                created_count += 1

        manifest.append({
            "filename": dest_filename,
            "source_path": str(src_path),
            "symlink_path": str(dest_link_path),
            "link_target": str(link_target),
            "file_size_bytes": str(src_path.stat().st_size) if src_path.exists() else "0",
            "status": link_status,
        })

    return manifest, created_count, skipped_count, error_count


def write_manifest_csv(manifest: List[Dict[str, str]], output_csv: Path, logger: Optional[logging.Logger] = None):
    """Write the symlink operation manifest to a CSV file."""
    if not manifest:
        return

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["filename", "source_path", "symlink_path", "link_target", "file_size_bytes", "status"]

    with open(output_csv, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest)

    if logger:
        logger.info(f"Wrote asset symlink manifest to: {output_csv.resolve()}")


# ===============================================================================
# SAM 2 MASK INGESTION & COCO 1.0 JSON GENERATION LOGIC
# ===============================================================================

def parse_mask_filename(mask_filename: str) -> Optional[Tuple[str, int, str]]:
    """
    Parse a SAM 2 mask filename into (voucher_stem, instance_index, label_name).
    Expected pattern: '<voucher_id>_inst<idx>_<label_name>.png'
    """
    match = re.match(r"^(.+)_inst(\d+)_(.+)\.(png|jpg|jpeg|tif|tiff|bmp)$", mask_filename, re.IGNORECASE)
    if match:
        voucher_id = match.group(1)
        instance_idx = int(match.group(2))
        label_name = match.group(3).strip()
        return voucher_id, instance_idx, label_name
    return None


def extract_contours_from_mask(
    mask_path: Path,
    min_contour_area: float = 1.0,
    simplify_tolerance: float = 0.0,
    logger: Optional[logging.Logger] = None,
) -> Tuple[List[List[float]], float, List[float], Tuple[int, int]]:
    """
    Load a binary segmentation mask, extract outer boundary contours using OpenCV,
    and convert them to flattened COCO polygon coordinates [x1, y1, x2, y2, ...].
    
    Returns:
        (polygons, total_area, bbox_xywh, (height, width))
    """
    if cv2 is None or np is None:
        raise ImportError("OpenCV (cv2) and NumPy (numpy) are required for mask contour processing.")

    mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
    if mask is None:
        raise ValueError(f"Unable to read mask image: {mask_path}")

    h, w = mask.shape[:2]
    binary = (mask > 127).astype(np.uint8)

    # Find external contours
    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    polygons: List[List[float]] = []
    total_area: float = 0.0
    all_points: List[np.ndarray] = []

    for cnt in contours:
        if len(cnt) < 3:
            continue

        if simplify_tolerance > 0.0:
            cnt = cv2.approxPolyDP(cnt, epsilon=simplify_tolerance, closed=True)
            if len(cnt) < 3:
                continue

        area = float(cv2.contourArea(cnt))
        if area < min_contour_area:
            continue

        poly = cnt.reshape(-1, 2).astype(float).flatten().tolist()
        if len(poly) >= 6:  # At least 3 (x, y) coordinate pairs
            polygons.append(poly)
            total_area += area
            all_points.append(cnt.reshape(-1, 2))

    if all_points:
        pts_concat = np.concatenate(all_points, axis=0)
        min_x = float(pts_concat[:, 0].min())
        min_y = float(pts_concat[:, 1].min())
        max_x = float(pts_concat[:, 0].max())
        max_y = float(pts_concat[:, 1].max())
        bbox = [
            round(min_x, 2),
            round(min_y, 2),
            round(max_x - min_x, 2),
            round(max_y - min_y, 2),
        ]
    else:
        bbox = [0.0, 0.0, 0.0, 0.0]

    return polygons, round(total_area, 2), bbox, (h, w)


def build_coco_dataset_from_masks(
    masks_dir: Path,
    image_search_dirs: List[Path],
    min_contour_area: float = 1.0,
    simplify_tolerance: float = 0.0,
    include_unannotated_images: bool = False,
    logger: Optional[logging.Logger] = None,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Ingest SAM 2 binary masks, extract polygon contours with OpenCV, apply strict PCD
    class mapping, omit 'cauline_leaf' background noise, and construct a complete COCO 1.0 JSON.

    Returns:
        (coco_dict, statistics_summary_dict)
    """
    if not masks_dir.exists():
        raise FileNotFoundError(f"Masks directory not found: {masks_dir.resolve()}")

    if logger:
        logger.info(f"Ingesting SAM 2 binary masks from: {masks_dir.resolve()}")

    mask_files = sorted(list(masks_dir.glob("*.png")) + list(masks_dir.glob("*.jpg")) + list(masks_dir.glob("*.tif")))
    if not mask_files:
        if logger:
            logger.warning(f"No mask image files found in {masks_dir.resolve()}")
        return {}, {}

    if logger:
        logger.info(f"Discovered {len(mask_files)} total mask files.")

    # Index available herbarium voucher sheet images across search directories
    image_index: Dict[str, Path] = {}
    for sdir in image_search_dirs:
        if sdir.exists():
            for img_p in sdir.glob("*.*"):
                if img_p.suffix.lower() in VALID_IMAGE_EXTENSIONS:
                    if img_p.stem not in image_index:
                        image_index[img_p.stem] = img_p

    if logger:
        logger.info(f"Indexed {len(image_index)} voucher sheet images across search paths.")

    # Group masks by voucher catalog number
    voucher_masks: Dict[str, List[Tuple[int, str, Path]]] = defaultdict(list)
    unmatched_masks: List[Path] = []
    raw_label_counts: Counter = Counter()

    for mf in mask_files:
        parsed = parse_mask_filename(mf.name)
        if parsed:
            voucher_id, inst_idx, label = parsed
            voucher_masks[voucher_id].append((inst_idx, label, mf))
            raw_label_counts[label] += 1
        else:
            unmatched_masks.append(mf)

    if unmatched_masks and logger:
        logger.warning(f"Found {len(unmatched_masks)} masks with unrecognized naming convention.")

    # Build COCO 1.0 Data Structures
    images: List[Dict[str, Any]] = []
    annotations: List[Dict[str, Any]] = []
    
    img_id_counter = 1
    ann_id_counter = 1

    stats = {
        "total_masks_ingested": len(mask_files),
        "unique_vouchers_annotated": len(voucher_masks),
        "target_classes_mapped": Counter(),
        "background_classes_omitted": Counter(),
        "total_annotations_created": 0,
        "degenerate_masks_skipped": 0,
        "total_images_in_dataset": 0,
    }

    # Iterate through each annotated voucher
    for voucher_id, inst_list in sorted(voucher_masks.items()):
        inst_list.sort(key=lambda x: x[0])  # Sort by instance index

        # Determine native sheet dimensions
        resolved_img = image_index.get(voucher_id)
        if resolved_img and resolved_img.exists() and cv2 is not None:
            # Quick dimension read from raster header or shape
            img_sample = cv2.imread(str(resolved_img))
            if img_sample is not None:
                sheet_h, sheet_w = img_sample.shape[:2]
                file_name = resolved_img.name
            else:
                # Fallback to mask dimension
                _, _, _, (sheet_h, sheet_w) = extract_contours_from_mask(inst_list[0][2], min_contour_area=1.0)
                file_name = f"{voucher_id}.jpg"
        else:
            # Fallback to mask shape
            _, _, _, (sheet_h, sheet_w) = extract_contours_from_mask(inst_list[0][2], min_contour_area=1.0)
            file_name = f"{voucher_id}.jpg"

        image_record = {
            "id": img_id_counter,
            "file_name": file_name,
            "width": int(sheet_w),
            "height": int(sheet_h),
            "license": 1,
            "date_captured": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        }
        images.append(image_record)
        current_img_id = img_id_counter
        img_id_counter += 1

        # Process instances for this sheet
        for inst_idx, label_name, mask_path in inst_list:
            # Strict Taxonomy Enforcement:
            # If class is 'cauline_leaf' or any non-target structure, intentionally omit from JSON
            if label_name not in PCD_CLASS_MAPPING:
                stats["background_classes_omitted"][label_name] += 1
                if logger:
                    logger.debug(
                        f"Voucher {voucher_id} inst {inst_idx:02d}: OMITTING '{label_name}' "
                        f"(treated as background noise for Detectron2)."
                    )
                continue

            target_class = PCD_CLASS_MAPPING[label_name]
            cat_id = PCD_CATEGORY_ID_MAP[target_class]

            polygons, area, bbox, _ = extract_contours_from_mask(
                mask_path,
                min_contour_area=min_contour_area,
                simplify_tolerance=simplify_tolerance,
                logger=logger,
            )

            if not polygons:
                stats["degenerate_masks_skipped"] += 1
                if logger:
                    logger.debug(f"Skipped degenerate/empty mask: {mask_path.name}")
                continue

            ann_record = {
                "id": ann_id_counter,
                "image_id": current_img_id,
                "category_id": cat_id,
                "segmentation": polygons,
                "area": area,
                "bbox": bbox,
                "iscrowd": 0,
            }
            annotations.append(ann_record)
            stats["target_classes_mapped"][target_class] += 1
            ann_id_counter += 1

    # Optional: Include unannotated images as pure negative background images
    if include_unannotated_images:
        annotated_stems = set(voucher_masks.keys())
        for stem, img_path in sorted(image_index.items()):
            if stem not in annotated_stems:
                img_sample = cv2.imread(str(img_path)) if cv2 is not None else None
                if img_sample is not None:
                    h, w = img_sample.shape[:2]
                else:
                    h, w = 4000, 3000
                images.append({
                    "id": img_id_counter,
                    "file_name": img_path.name,
                    "width": int(w),
                    "height": int(h),
                    "license": 1,
                    "date_captured": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                })
                img_id_counter += 1

    stats["total_annotations_created"] = len(annotations)
    stats["total_images_in_dataset"] = len(images)

    coco_dataset = {
        "info": {
            "description": "LeafMachine2 Plant Component Detector (PCD) Packera Dataset (SAM 2 Derived)",
            "url": "https://github.com/Plantatsrophe/Packera_dubia_morphometrics",
            "version": "1.0",
            "year": datetime.now().year,
            "contributor": "University of North Carolina at Chapel Hill Herbarium (NCU)",
            "date_created": datetime.now().strftime("%Y-%m-%d"),
        },
        "licenses": [
            {
                "id": 1,
                "name": "Research / Educational Non-Commercial Use",
                "url": "",
            }
        ],
        "images": images,
        "annotations": annotations,
        "categories": PCD_COCO_CATEGORIES,
    }

    return coco_dataset, stats


def split_coco_dataset(
    coco_dataset: Dict[str, Any],
    val_split: float = 0.20,
    seed: int = 42,
    logger: Optional[logging.Logger] = None,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Split a COCO dataset into training and validation sets by image ID.
    Ensures that all annotations for an image stay together in the same split.
    """
    random.seed(seed)
    all_images = list(coco_dataset.get("images", []))
    all_annotations = list(coco_dataset.get("annotations", []))
    categories = list(coco_dataset.get("categories", []))
    info = dict(coco_dataset.get("info", {}))
    licenses = list(coco_dataset.get("licenses", []))

    # Shuffle image list
    shuffled_images = list(all_images)
    random.shuffle(shuffled_images)

    val_count = int(round(len(shuffled_images) * val_split))
    val_count = max(1, val_count) if (val_split > 0 and len(shuffled_images) > 1) else 0

    val_images = shuffled_images[:val_count]
    train_images = shuffled_images[val_count:]

    val_img_ids = {img["id"] for img in val_images}
    train_img_ids = {img["id"] for img in train_images}

    train_annotations = [ann for ann in all_annotations if ann["image_id"] in train_img_ids]
    val_annotations = [ann for ann in all_annotations if ann["image_id"] in val_img_ids]

    train_coco = {
        "info": {**info, "description": f"{info.get('description', '')} [Train Split]"},
        "licenses": licenses,
        "images": train_images,
        "annotations": train_annotations,
        "categories": categories,
    }

    val_coco = {
        "info": {**info, "description": f"{info.get('description', '')} [Validation Split]"},
        "licenses": licenses,
        "images": val_images,
        "annotations": val_annotations,
        "categories": categories,
    }

    if logger:
        logger.info(
            f"COCO Dataset Split ({100*(1-val_split):.0f}% Train / {100*val_split:.0f}% Val): "
            f"Train={len(train_images)} imgs ({len(train_annotations)} anns), "
            f"Val={len(val_images)} imgs ({len(val_annotations)} anns)"
        )

    return train_coco, val_coco


def export_coco_json(coco_data: Dict[str, Any], output_path: Path, logger: Optional[logging.Logger] = None):
    """Serialize and write COCO JSON file with standard formatting."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(coco_data, f, indent=2)
    if logger:
        logger.info(f"Successfully exported COCO 1.0 JSON to: {output_path.resolve()}")


def log_dataset_summary(stats: Dict[str, Any], logger: logging.Logger):
    """Print comprehensive summary report of the COCO conversion."""
    logger.info("=================================================================")
    logger.info("       LEAFMACHINE2 PCD / SAM 2 COCO DATASET SUMMARY             ")
    logger.info("=================================================================")
    logger.info(f" Total SAM 2 Mask Files Ingested    : {stats.get('total_masks_ingested', 0):,}")
    logger.info(f" Unique Voucher Sheets Annotated    : {stats.get('unique_vouchers_annotated', 0):,}")
    logger.info(f" Total Herbarium Sheets in Dataset  : {stats.get('total_images_in_dataset', 0):,}")
    logger.info(f" Total COCO Annotations Generated   : {stats.get('total_annotations_created', 0):,}")
    
    logger.info("")
    logger.info("--- INCLUDED PCD TARGET CATEGORIES (Detectron2 Supervised) ---")
    target_classes = stats.get("target_classes_mapped", {})
    if target_classes:
        for cat_name, count in target_classes.items():
            cat_id = PCD_CATEGORY_ID_MAP.get(cat_name, "?")
            logger.info(f"  * [Category ID {cat_id}] {cat_name:<16} : {count:>5} instances")
    else:
        logger.info("  (None mapped)")

    logger.info("")
    logger.info("--- INTENTIONALLY OMITTED CLASSES (Detectron2 Background Noise) ---")
    omitted_classes = stats.get("background_classes_omitted", {})
    if omitted_classes:
        for cat_name, count in omitted_classes.items():
            logger.info(f"  * {cat_name:<22} : {count:>5} instances (Omitted from JSON)")
    else:
        logger.info("  (None omitted)")

    if stats.get("degenerate_masks_skipped", 0) > 0:
        logger.info(f" Degenerate / Empty Masks Skipped   : {stats.get('degenerate_masks_skipped', 0)}")
    logger.info("=================================================================")


# ===============================================================================
# CLI ARGUMENT PARSER & ENTRYPOINT
# ===============================================================================

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments with backwards-compatible aliases."""
    project_root = get_project_root()
    default_raw_dir = project_root / "data" / "raw_vouchers"
    default_lm2_dir = project_root / "LM2_Project"
    default_masks_dir = project_root / "data" / "raw_annotations" / "masks"
    default_coco_output = project_root / "LM2_Project" / "Data" / "annotations" / "coco_pcd_packera.json"

    parser = argparse.ArgumentParser(
        description="LeafMachine2 Dataset Pre-Processor & SAM 2 to COCO 1.0 JSON Dataset Converter.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Group 1: Image Staging & Symlink Arguments
    staging_group = parser.add_argument_group("Herbarium Image Staging & Symlinks")
    staging_group.add_argument(
        "--input-dirs",
        "--image-src",
        "-i",
        dest="input_dirs",
        nargs="+",
        type=Path,
        default=[default_raw_dir],
        help="One or more directories containing raw herbarium voucher sheet images.",
    )
    staging_group.add_argument(
        "--lm2-root",
        "--output",
        "-o",
        dest="lm2_root",
        type=Path,
        default=default_lm2_dir,
        help="Root directory for LeafMachine2 project (LM2_Project).",
    )
    staging_group.add_argument(
        "--images-subdir",
        type=str,
        default="images",
        help="Subdirectory inside LM2_Project/Data where image symlinks will reside.",
    )
    staging_group.add_argument(
        "--output-subdir",
        type=str,
        default="output",
        help="Subdirectory inside LM2_Project/Data for LM2 analysis output.",
    )
    staging_group.add_argument(
        "--annotations-subdir",
        type=str,
        default="annotations",
        help="Subdirectory inside LM2_Project/Data where COCO annotations will reside.",
    )
    staging_group.add_argument(
        "--relative",
        action="store_true",
        help="Create relative symlinks instead of absolute symlinks.",
    )
    staging_group.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing symlinks or files at destination.",
    )
    staging_group.add_argument(
        "--verify-images",
        action="store_true",
        help="Perform PIL image header verification on each file before symlinking.",
    )
    staging_group.add_argument(
        "--limit",
        "-n",
        type=int,
        default=None,
        help="Limit number of images to symlink (useful for testing).",
    )
    staging_group.add_argument(
        "--skip-symlinks",
        action="store_true",
        help="Skip image symlink creation and only perform SAM 2 to COCO conversion.",
    )

    # Group 2: SAM 2 Binary Mask & COCO Conversion Arguments
    coco_group = parser.add_argument_group("SAM 2 Mask Ingestion & COCO 1.0 Conversion")
    coco_group.add_argument(
        "--masks-dir",
        "-m",
        type=Path,
        default=default_masks_dir,
        help="Directory containing binary PNG masks exported by SAM 2.",
    )
    coco_group.add_argument(
        "--coco-output",
        "--output-coco",
        dest="coco_output",
        type=Path,
        default=default_coco_output,
        help="Output destination path for the exported COCO 1.0 JSON dataset.",
    )
    coco_group.add_argument(
        "--val-split",
        type=float,
        default=0.0,
        help="Optional validation split ratio (e.g. 0.20 for 80/20 train/val split).",
    )
    coco_group.add_argument(
        "--train-coco-output",
        type=Path,
        default=None,
        help="Custom output path for train split COCO JSON (used if --val-split > 0).",
    )
    coco_group.add_argument(
        "--val-coco-output",
        type=Path,
        default=None,
        help="Custom output path for validation split COCO JSON (used if --val-split > 0).",
    )
    coco_group.add_argument(
        "--min-contour-area",
        type=float,
        default=1.0,
        help="Minimum contour area in pixels to retain (filters out tiny island noise).",
    )
    coco_group.add_argument(
        "--simplify-tolerance",
        type=float,
        default=0.0,
        help="Contour approximation epsilon for Douglas-Peucker simplification (0.0 = exact boundaries).",
    )
    coco_group.add_argument(
        "--include-unannotated-images",
        action="store_true",
        help="Include unannotated raw voucher images as negative background samples in COCO.",
    )
    coco_group.add_argument(
        "--skip-coco",
        action="store_true",
        help="Skip SAM 2 mask to COCO conversion and only create symlinks.",
    )
    coco_group.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible train/val splits.",
    )

    # Group 3: General Execution Flags
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Simulate the workflow without writing files or creating symlinks.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose debug logging.",
    )

    return parser.parse_args()


def main():
    args = parse_arguments()
    logger = setup_logging(verbose=args.verbose)

    logger.info("=================================================================")
    logger.info(" LeafMachine2 Dataset Pre-Processor & SAM 2 COCO Builder")
    logger.info("=================================================================")
    logger.info(f"Target LM2 Project Root : {args.lm2_root.resolve()}")
    logger.info(f"Dry Run Mode            : {args.dry_run}")

    # Step 1: Ensure directory structure
    images_dir, output_dir, annotations_dir, configs_dir = create_lm2_directories(
        lm2_root=args.lm2_root,
        images_dir_name=args.images_subdir,
        output_dir_name=args.output_subdir,
        configs_dir_name="configs",
        annotations_dir_name=args.annotations_subdir,
        dry_run=args.dry_run,
        logger=logger,
    )

    # Step 2: Handle Herbarium Sheet Image Staging & Symlinks
    if not args.skip_symlinks:
        logger.info("\n--- Phase 1: Herbarium Sheet Image Staging & Symlink Creation ---")
        valid_images, filtered_out = collect_raw_images(
            input_dirs=args.input_dirs,
            verify_images=args.verify_images,
            logger=logger,
        )

        logger.info(f"Discovered {len(valid_images)} valid herbarium sheet image(s).")
        if filtered_out:
            logger.info(f"Filtered out {len(filtered_out)} non-image / invalid file(s).")

        if valid_images:
            manifest, created, skipped, errors = create_symlinks(
                source_images=valid_images,
                target_dir=images_dir,
                relative=args.relative,
                overwrite=args.overwrite,
                limit=args.limit,
                dry_run=args.dry_run,
                logger=logger,
            )

            manifest_path = args.lm2_root / "Data" / "symlink_manifest.csv"
            if not args.dry_run:
                write_manifest_csv(manifest, manifest_path, logger=logger)

            logger.info(
                f"Symlinks: Created/Updated={created}, Skipped={skipped}, Errors={errors}"
            )
        else:
            logger.warning("No valid images found for symlinking.")
    else:
        logger.info("Skipping image symlink creation as requested (--skip-symlinks).")

    # Step 3: Handle SAM 2 Binary Mask Ingestion & COCO 1.0 Dataset Generation
    if not args.skip_coco:
        logger.info("\n--- Phase 2: SAM 2 Binary Mask Ingestion & COCO 1.0 Dataset Generation ---")
        search_dirs = [images_dir] + args.input_dirs + [get_project_root() / "data" / "raw_vouchers_quarantine"]
        
        coco_dataset, stats = build_coco_dataset_from_masks(
            masks_dir=args.masks_dir,
            image_search_dirs=search_dirs,
            min_contour_area=args.min_contour_area,
            simplify_tolerance=args.simplify_tolerance,
            include_unannotated_images=args.include_unannotated_images,
            logger=logger,
        )

        if coco_dataset and not args.dry_run:
            if args.val_split > 0.0:
                train_coco, val_coco = split_coco_dataset(
                    coco_dataset,
                    val_split=args.val_split,
                    seed=args.seed,
                    logger=logger,
                )

                train_path = (
                    args.train_coco_output
                    if args.train_coco_output
                    else args.coco_output.parent / f"{args.coco_output.stem}_train{args.coco_output.suffix}"
                )
                val_path = (
                    args.val_coco_output
                    if args.val_coco_output
                    else args.coco_output.parent / f"{args.coco_output.stem}_val{args.coco_output.suffix}"
                )

                export_coco_json(train_coco, train_path, logger=logger)
                export_coco_json(val_coco, val_path, logger=logger)
            else:
                export_coco_json(coco_dataset, args.coco_output, logger=logger)

            log_dataset_summary(stats, logger=logger)
        elif args.dry_run:
            logger.info("[Dry Run] Simulated COCO dataset generation without writing to disk.")
            if stats:
                log_dataset_summary(stats, logger=logger)
    else:
        logger.info("Skipping SAM 2 to COCO conversion as requested (--skip-coco).")

    logger.info("\nDataset preparation workflow completed successfully.")


if __name__ == "__main__":
    main()
