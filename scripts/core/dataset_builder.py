
import os
import sys
import shutil
import cv2
import math
import numpy as np
import random
import yaml
import json
import glob
import logging
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from pathlib import Path

# Common imports for the project
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE, DEFAULT_RAW_DIR, DEFAULT_CURATED_CSV, DEFAULT_OUTPUT_DIR, DEFAULT_CONFIG_PATH, DEFAULT_QC_DIR
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

import argparse
from scripts.core.artifact_harvester import ArtifactHarvester
from scripts.core.augmentation import extract_hard_negative_background_crop, SyntheticOcclusionAugmenter
from scripts.core.botanical_annotations import extract_botanical_annotations
from scripts.core.dataset_utils import stratify_and_partition_dataset, render_qc_verification_overlay, export_ultralytics_dataset_yaml

logger = setup_logging()

def build_artifact_robust_dataset(
    raw_vouchers_dir: Path = DEFAULT_RAW_DIR,
    annotations_dir: Optional[Path] = None,
    curated_csv_path: Path = DEFAULT_CURATED_CSV,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    config_yaml_path: Path = DEFAULT_CONFIG_PATH,
    qc_output_dir: Path = DEFAULT_QC_DIR,
    negative_ratio: float = 0.09,
    augment_prob: float = 0.75,
    copy_paste_prob: Optional[float] = None,
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    min_instance_area: float = 150.0,
    num_qc_plots: int = 25,
    num_workers: int = 4,
    limit: Optional[int] = None,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Main orchestration routine for building the artifact-robust dataset.
    """
    if copy_paste_prob is not None:
        augment_prob = copy_paste_prob

    logger.info("Initializing Artifact-Robust Botanical Dataset Builder...")
    logger.info(f"YOLO Schema: {CLASS_NAMES}")
    logger.info(f"Target Hard Negative Proportion: {negative_ratio * 100:.1f}%")
    logger.info(f"Augmentation Probability: {augment_prob * 100:.1f}% | Min Instance Area: {min_instance_area} px")

    # Set seeds for complete reproducibility
    random.seed(seed)
    np.random.seed(seed)

    # 0. Purge existing dataset directory to prevent ghost label contamination
    if output_dir.exists():
        logger.warning(f"Purging existing dataset output directory to prevent stale label contamination: {output_dir}")
        import shutil
        shutil.rmtree(output_dir)

    # 1. Discover voucher records
    vouchers: List[Dict[str, Any]] = []
    if curated_csv_path.exists():
        df = pd.read_csv(curated_csv_path)
        logger.info(f"Loaded {len(df)} records from curated vouchers table: {curated_csv_path}")
        for _, row in df.iterrows():
            cat_num = str(row.get("catalogNumber", ""))
            img_path_str = str(row.get("image_path", ""))
            img_path = DEFAULT_WORKSPACE / img_path_str if not os.path.isabs(img_path_str) else Path(img_path_str)
            if not img_path.exists():
                img_path = raw_vouchers_dir / f"{cat_num}.jpg"
            if img_path.exists():
                vouchers.append({
                    "catalogNumber": cat_num,
                    "institutionCode": row.get("institutionCode", "NCU"),
                    "determiner_tier": row.get("determiner_tier", "Tier_1_Gold"),
                    "image_path": img_path
                })
    else:
        logger.warning(f"Curated CSV not found at {curated_csv_path}. Scanning raw directory directly...")
        img_files = sorted(glob.glob(str(raw_vouchers_dir / "*.jpg")))
        for p_str in img_files:
            p = Path(p_str)
            cat_num = p.stem
            vouchers.append({
                "catalogNumber": cat_num,
                "institutionCode": cat_num[:3] if cat_num[:3].isalpha() else "NCU",
                "determiner_tier": "Tier_1_Gold",
                "image_path": p
            })

    if limit and limit > 0:
        vouchers = vouchers[:limit]
        logger.info(f"Limiting processing to first {limit} vouchers as requested.")

    if not vouchers:
        logger.error(f"No valid voucher images located in {raw_vouchers_dir}!")
        return {"status": "error", "message": "No voucher images found"}

    # 2. Stratified Partitioning
    partitions = stratify_and_partition_dataset(
        vouchers, train_ratio=train_ratio, val_ratio=val_ratio, test_ratio=test_ratio, rng_seed=seed
    )

    # 3. Create YOLO directory hierarchy
    for split in ["train", "val", "test"]:
        (output_dir / "images" / split).mkdir(parents=True, exist_ok=True)
        (output_dir / "labels" / split).mkdir(parents=True, exist_ok=True)

    # 4. Initialize Harvester & Augmenter
    harvester = ArtifactHarvester(rng_seed=seed)
    augmenter = SyntheticOcclusionAugmenter(harvester=harvester, rng_seed=seed)

    # First pass: Ingest clean artifacts from vouchers across all partitions
    logger.info("Pass 1/2: Harvesting artifact patches across voucher collection...")
    for v in tqdm(vouchers, desc="Harvesting Artifacts"):
        img = cv2.imread(str(v["image_path"]))
        if img is not None:
            harvester.detect_and_extract_sheet_artifacts(img, catalog_number=v["catalogNumber"])

    logger.info(
        f"Artifact Bank Populated: " +
        ", ".join([f"{k}: {len(v)}" for k, v in harvester.artifact_bank.items()])
    )

    # Second pass: Generate images and labels with copy-paste augmentations and hard negatives
    logger.info("Pass 2/2: Building YOLO partitions, annotations, and negative injections...")
    manifest_records: List[Dict[str, Any]] = []
    qc_count = 0

    total_stats = {
        "train_samples": 0, "val_samples": 0, "test_samples": 0,
        "negative_samples": 0, "total_instances": 0, "class_counts": {c: 0 for c in CLASS_NAMES}
    }

    for split_name, split_vouchers in partitions.items():
        logger.info(f"Processing split '{split_name}' ({len(split_vouchers)} base vouchers)...")

        # Calculate number of hard negative background tiles to generate for this split
        num_negatives = max(1, int(round(len(split_vouchers) * negative_ratio)))

        for idx, v in enumerate(tqdm(split_vouchers, desc=f"Split [{split_name}]")):
            cat_num = v["catalogNumber"]
            img_path = v["image_path"]
            img = cv2.imread(str(img_path))
            if img is None:
                continue

            h, w = img.shape[:2]

            # 1. Extract natural sheet artifacts
            art_anns = harvester.detect_and_extract_sheet_artifacts(img, catalog_number=cat_num)

            # 2. Extract botanical organs
            bot_anns = extract_botanical_annotations(
                img, artifact_anns=art_anns, min_instance_area=min_instance_area
            )

            all_anns = art_anns + bot_anns

            # 3. Apply Synthetic Copy-Paste Occlusion Augmentations (primarily in train split)
            if split_name == "train":
                aug_img, final_anns = augmenter.apply_copy_paste_augmentation(
                    img, all_anns, paste_probability=augment_prob, max_pastes_per_image=3
                )
            else:
                aug_img = img
                final_anns = all_anns

            # 4. Save Image and YOLO Label .txt
            out_img_name = f"{cat_num}.jpg"
            out_txt_name = f"{cat_num}.txt"
            dest_img_path = output_dir / "images" / split_name / out_img_name
            dest_txt_path = output_dir / "labels" / split_name / out_txt_name

            cv2.imwrite(str(dest_img_path), aug_img)

            # Write YOLO segmentation lines
            with open(dest_txt_path, "w", encoding="utf-8") as f_lbl:
                for ann in final_anns:
                    line = ann.to_yolo_seg_line(img_w=w, img_h=h)
                    f_lbl.write(f"{line}\n")
                    total_stats["class_counts"][ann.class_name] += 1
                    total_stats["total_instances"] += 1

            total_stats[f"{split_name}_samples"] += 1

            manifest_records.append({
                "catalogNumber": cat_num,
                "split": split_name,
                "is_negative": False,
                "num_instances": len(final_anns),
                "image_path": str(dest_img_path),
                "label_path": str(dest_txt_path)
            })

            # Render QC verification overlay for top N samples
            if qc_count < num_qc_plots:
                qc_path = qc_output_dir / f"qc_{split_name}_{cat_num}.jpg"
                render_qc_verification_overlay(aug_img, final_anns, catalog_number=cat_num, output_path=qc_path)
                qc_count += 1

        # 5. Hard Negative Injection (Pure background sheet regions with empty .txt label)
        logger.info(f"Injecting {num_negatives} hard negative background sheets into split '{split_name}'...")
        donor_pool = split_vouchers if split_vouchers else vouchers
        if not donor_pool:
            continue

        for neg_idx in range(num_negatives):
            # Select background donor voucher safely from available pool
            donor_v = donor_pool[neg_idx % len(donor_pool)]
            donor_img = cv2.imread(str(donor_v["image_path"]))
            if donor_img is None:
                continue

            neg_crop = extract_hard_negative_background_crop(
                donor_img, all_annotations=[], crop_size=(1024, 1024)
            )
            if neg_crop is None:
                continue

            neg_name = f"neg_sheet_{donor_v['catalogNumber']}_{neg_idx:02d}"
            neg_img_path = output_dir / "images" / split_name / f"{neg_name}.jpg"
            neg_txt_path = output_dir / "labels" / split_name / f"{neg_name}.txt"

            cv2.imwrite(str(neg_img_path), neg_crop)
            # Create an explicitly EMPTY .txt annotation file (0 lines)
            with open(neg_txt_path, "w", encoding="utf-8") as f_neg:
                pass  # Empty file signals true background negative to YOLO

            total_stats["negative_samples"] += 1
            total_stats[f"{split_name}_samples"] += 1

            manifest_records.append({
                "catalogNumber": neg_name,
                "split": split_name,
                "is_negative": True,
                "num_instances": 0,
                "image_path": str(neg_img_path),
                "label_path": str(neg_txt_path)
            })

    # 6. Export Ultralytics dataset YAML
    export_ultralytics_dataset_yaml(output_dir, config_yaml_path)

    # 7. Export dataset manifest table & summary report
    manifest_csv_path = DEFAULT_WORKSPACE / "data" / "tables" / "dataset_manifest.csv"
    manifest_csv_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(manifest_records).to_csv(manifest_csv_path, index=False)

    summary_report_path = DEFAULT_WORKSPACE / "outputs" / "reports" / "dataset_build_summary.json"
    summary_report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_report_path, "w", encoding="utf-8") as f_rep:
        json.dump({
            "total_samples": len(manifest_records),
            "stats": total_stats,
            "config_yaml": str(config_yaml_path),
            "qc_plots_dir": str(qc_output_dir)
        }, f_rep, indent=2)

    logger.info("=" * 78)
    logger.info("ARTIFACT-ROBUST DATASET BUILD COMPLETED SUCCESSFULLY")
    logger.info(f"Total Samples Generated: {len(manifest_records)}")
    logger.info(f"Train: {total_stats['train_samples']} | Val: {total_stats['val_samples']} | Test: {total_stats['test_samples']}")
    logger.info(f"Hard Negative Background Sheets: {total_stats['negative_samples']} ({total_stats['negative_samples']/max(1, len(manifest_records))*100:.1f}%)")
    logger.info(f"Total Instances Annotated: {total_stats['total_instances']}")
    logger.info(f"Class Distribution: {total_stats['class_counts']}")
    logger.info(f"Config YAML: {config_yaml_path}")
    logger.info(f"QC Overlay Figures: {qc_output_dir}")
    logger.info("=" * 78)

    return {
        "status": "success",
        "total_samples": len(manifest_records),
        "stats": total_stats,
        "config_yaml": str(config_yaml_path),
        "qc_dir": str(qc_output_dir)
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Construct an artifact-robust multi-class YOLO segmentation dataset for botanical herbarium specimens."
    )
    parser.add_argument(
        "--raw-images-dir", "--raw-dir", dest="raw_dir", type=Path, default=DEFAULT_RAW_DIR,
        help="Path to directory containing raw voucher JPG images."
    )
    parser.add_argument(
        "--annotations-dir", type=Path, default=None,
        help="Optional directory containing pre-computed raw annotations."
    )
    parser.add_argument(
        "--curated-csv", type=Path, default=DEFAULT_CURATED_CSV,
        help="Path to curated_vouchers.csv containing metadata for stratification."
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR,
        help="Destination directory for the formatted YOLO dataset."
    )
    parser.add_argument(
        "--config-yaml", type=Path, default=DEFAULT_CONFIG_PATH,
        help="Path to export the Ultralytics dataset_config.yaml."
    )
    parser.add_argument(
        "--qc-dir", type=Path, default=DEFAULT_QC_DIR,
        help="Directory to save QC overlay visualization figures."
    )
    parser.add_argument(
        "--negative-ratio", type=float, default=0.09,
        help="Proportion of dataset comprising pure background negative sheets (default: 0.09 / ~9%%)."
    )
    parser.add_argument(
        "--augment-prob", type=float, default=0.75,
        help="Probability of applying synthetic copy-paste augmentations to training samples."
    )
    parser.add_argument(
        "--copy-paste-prob", type=float, default=None,
        help="Alias for copy-paste augmentation probability."
    )
    parser.add_argument(
        "--train-ratio", type=float, default=0.70,
        help="Fraction of dataset allocated to training partition (default 0.70)."
    )
    parser.add_argument(
        "--val-ratio", type=float, default=0.15,
        help="Fraction of dataset allocated to validation partition (default 0.15)."
    )
    parser.add_argument(
        "--test-ratio", type=float, default=0.15,
        help="Fraction of dataset allocated to test partition (default 0.15)."
    )
    parser.add_argument(
        "--min-instance-area", type=float, default=150.0,
        help="Minimum bounding/contour area in square pixels for a botanical instance (default 150)."
    )
    parser.add_argument(
        "--num-qc-plots", type=int, default=25,
        help="Number of verification QC overlay images to generate in qc-dir."
    )
    parser.add_argument(
        "--num-workers", type=int, default=4,
        help="Number of processing threads/workers (default 4)."
    )
    parser.add_argument(
        "--limit", type=int, default=1500,
        help="Optional limit on the number of vouchers to process."
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for deterministic stratification and augmentations."
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable detailed debug logging."
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    build_artifact_robust_dataset(
        raw_vouchers_dir=Path(args.input_dir),
        curated_csv_path=Path(args.curated_csv),
        output_dir=Path(args.output_dir),
        limit=args.limit,
        num_workers=args.concurrency
    )
