"""
scripts/core/dataset_builder.py
===============================
Artifact-Robust Botanical Dataset Builder module.
Constructs a verified multi-class instance segmentation YOLO dataset for
herbarium specimen phenotyping (Packera dubia complex).
"""

from __future__ import annotations

import argparse
import datetime
import json
import logging
import os
import random
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import cv2
import numpy as np
import pandas as pd
import yaml

from scripts.core.botanical_annotations import (
    BOTANICAL_CLASS_MAP,
    BOTANICAL_CLASSES,
    CLASS_COLORS_BGR,
    CLASS_NAMES,
    BotanicalAnnotation,
    VerifiedBotanicalLabelParser,
)
from scripts.core.config import (
    DEFAULT_CONFIG_PATH,
    DEFAULT_CURATED_CSV,
    DEFAULT_OUTPUT_DIR,
    DEFAULT_QC_DIR,
    DEFAULT_RAW_DIR,
    DEFAULT_WORKSPACE,
)
from scripts.core.dataset_utils import (
    render_qc_verification_overlay,
    stratify_and_partition_dataset,
)
from scripts.core.logger import setup_logging

logger = setup_logging()

DEFAULT_ANNOTATIONS_DIR = DEFAULT_WORKSPACE / "data" / "raw_annotations"

__all__ = [
    "build_artifact_robust_dataset",
    "extract_paper_background_crop",
    "run_synthetic_test_suite",
    "parse_args",
    "BOTANICAL_CLASSES",
    "BOTANICAL_CLASS_MAP",
    "CLASS_NAMES",
    "CLASS_COLORS_BGR",
    "VerifiedBotanicalLabelParser",
    "BotanicalAnnotation",
]


def extract_paper_background_crop(
    image_bgr: np.ndarray,
    annotations: List[BotanicalAnnotation],
    crop_size: Tuple[int, int] = (1024, 1024),
    rng: Optional[random.Random] = None
) -> Optional[np.ndarray]:
    """
    Extracts a pure background paper crop containing 0 plant organs or annotations.
    """
    if rng is None:
        rng = random.Random()

    h, w = image_bgr.shape[:2]
    cw, ch = crop_size
    if w <= cw or h <= ch:
        return None

    occupied = np.zeros((h, w), dtype=np.uint8)
    for ann in annotations:
        if len(ann.polygon) >= 3:
            pts = ann.polygon.astype(np.int32)
            cv2.fillPoly(occupied, [pts], 255)
        else:
            x1, y1, x2, y2 = [int(v) for v in ann.bbox]
            cv2.rectangle(occupied, (x1, y1), (x2, y2), 255, -1)

    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (25, 25))
    occupied = cv2.dilate(occupied, kernel)

    for _ in range(30):
        rx = rng.randint(0, w - cw)
        ry = rng.randint(0, h - ch)
        roi_mask = occupied[ry:ry + ch, rx:rx + cw]
        if cv2.countNonZero(roi_mask) == 0:
            crop = image_bgr[ry:ry + ch, rx:rx + cw].copy()
            return crop

    rx = min(w - cw, max(0, w - cw - 20))
    ry = 20
    return image_bgr[ry:ry + ch, rx:rx + cw].copy()


def build_artifact_robust_dataset(
    raw_vouchers_dir: Path = DEFAULT_RAW_DIR,
    annotations_dir: Optional[Path] = DEFAULT_ANNOTATIONS_DIR,
    curated_csv_path: Path = DEFAULT_CURATED_CSV,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    config_yaml_path: Path = DEFAULT_CONFIG_PATH,
    qc_output_dir: Path = DEFAULT_QC_DIR,
    negative_ratio: float = 0.09,
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    min_instance_area: float = 50.0,
    num_qc_plots: int = 25,
    limit: Optional[int] = None,
    seed: int = 42
) -> Dict[str, Any]:
    """
    Main orchestration routine for building the 7-class botanical YOLO dataset.
    """
    logger.info("=" * 78)
    logger.info("STARTING BOTANICAL ARTIFACT-ROBUST DATASET BUILD")
    logger.info(f"Target 7-Class Ontology: {BOTANICAL_CLASSES}")
    logger.info(f"Splits: Train={train_ratio*100:.0f}% | Val={val_ratio*100:.0f}% | Test={test_ratio*100:.0f}%")
    logger.info(f"Hard Negative Background Ratio: {negative_ratio*100:.1f}%")
    logger.info("=" * 78)

    random.seed(seed)
    np.random.seed(seed)
    rng = random.Random(seed)

    raw_vouchers_dir = Path(raw_vouchers_dir)
    curated_csv_path = Path(curated_csv_path)
    output_dir = Path(output_dir)
    config_yaml_path = Path(config_yaml_path)
    qc_output_dir = Path(qc_output_dir)
    if annotations_dir:
        annotations_dir = Path(annotations_dir)

    # 1. Clean & Prepare YOLO Dataset Output Directory
    if output_dir.exists():
        logger.warning(f"Purging existing dataset output directory: {output_dir}")
        shutil.rmtree(output_dir)

    for split in ["train", "val", "test"]:
        (output_dir / "images" / split).mkdir(parents=True, exist_ok=True)
        (output_dir / "labels" / split).mkdir(parents=True, exist_ok=True)

    # 2. Discover Voucher Images & Records
    vouchers: List[Dict[str, Any]] = []
    if curated_csv_path.exists():
        df = pd.read_csv(curated_csv_path)
        logger.info(f"Loaded {len(df)} metadata records from {curated_csv_path}")
        for _, row in df.iterrows():
            cat_num = str(row.get("catalogNumber", "")).strip()
            if not cat_num:
                continue
            img_path_str = str(row.get("image_path", ""))
            img_path = DEFAULT_WORKSPACE / img_path_str if not os.path.isabs(img_path_str) else Path(img_path_str)
            if not img_path.exists():
                img_path = raw_vouchers_dir / f"{cat_num}.jpg"
            if img_path.exists():
                vouchers.append({
                    "catalogNumber": cat_num,
                    "institutionCode": str(row.get("institutionCode", cat_num[:3] if cat_num[:3].isalpha() else "NCU")),
                    "image_path": img_path
                })
    else:
        logger.info(f"Curated CSV not found at {curated_csv_path}. Scanning {raw_vouchers_dir} directly...")
        if raw_vouchers_dir.exists():
            for p in sorted(raw_vouchers_dir.glob("*.jpg")):
                cat_num = p.stem
                inst = cat_num[:3] if cat_num[:3].isalpha() else "NCU"
                vouchers.append({
                    "catalogNumber": cat_num,
                    "institutionCode": inst,
                    "image_path": p
                })

    if limit and limit > 0:
        vouchers = vouchers[:limit]
        logger.info(f"Applied limit: processing {len(vouchers)} vouchers.")

    if not vouchers:
        logger.error(f"No valid voucher images found in {raw_vouchers_dir}")
        return {"status": "error", "message": "No voucher images found"}

    # 3. Initialize Verified Label Parser
    parser = VerifiedBotanicalLabelParser(min_instance_area=min_instance_area)

    # 4. Stratified Partitioning (70/15/15)
    partitions = stratify_and_partition_dataset(
        vouchers, train_ratio=train_ratio, val_ratio=val_ratio, test_ratio=test_ratio, rng_seed=seed
    )

    manifest_records: List[Dict[str, Any]] = []
    qc_count = 0
    stats = {
        "train_samples": 0, "val_samples": 0, "test_samples": 0,
        "negative_samples": 0, "total_instances": 0,
        "class_counts": {name: 0 for name in CLASS_NAMES}
    }

    # 5. Process Vouchers by Partition
    for split_name, split_vouchers in partitions.items():
        logger.info(f"Building split '{split_name}' with {len(split_vouchers)} specimen sheets...")

        for v in split_vouchers:
            cat_num = v["catalogNumber"]
            img_path = v["image_path"]
            img = cv2.imread(str(img_path))
            if img is None:
                continue

            h, w = img.shape[:2]
            annotations: List[BotanicalAnnotation] = []

            if annotations_dir and annotations_dir.exists():
                txt_candidates = [
                    annotations_dir / f"{cat_num}.txt",
                    annotations_dir / f"{cat_num}.yolo.txt",
                    annotations_dir / split_name / f"{cat_num}.txt",
                ]
                for txt_p in txt_candidates:
                    if txt_p.exists():
                        annotations.extend(parser.parse_yolo_txt(txt_p, img_w=w, img_h=h))
                        break

                json_candidates = [
                    annotations_dir / f"{cat_num}.json",
                    annotations_dir / split_name / f"{cat_num}.json",
                ]
                for json_p in json_candidates:
                    if json_p.exists():
                        annotations.extend(parser.parse_json_file(json_p, img_w=w, img_h=h))
                        break

            dest_img_path = output_dir / "images" / split_name / f"{cat_num}.jpg"
            dest_txt_path = output_dir / "labels" / split_name / f"{cat_num}.txt"

            cv2.imwrite(str(dest_img_path), img)

            with open(dest_txt_path, "w", encoding="utf-8") as f_lbl:
                for ann in annotations:
                    line = ann.to_yolo_seg_line(img_w=w, img_h=h)
                    if line:
                        f_lbl.write(f"{line}\n")
                        stats["class_counts"][ann.class_name] += 1
                        stats["total_instances"] += 1

            stats[f"{split_name}_samples"] += 1
            manifest_records.append({
                "catalogNumber": cat_num,
                "split": split_name,
                "is_negative": False,
                "num_instances": len(annotations),
                "image_path": str(dest_img_path),
                "label_path": str(dest_txt_path)
            })

            if qc_count < num_qc_plots and len(annotations) > 0:
                qc_img_path = qc_output_dir / f"qc_{split_name}_{cat_num}.jpg"
                render_qc_verification_overlay(
                    image_bgr=img,
                    annotations=annotations,
                    catalog_number=cat_num,
                    output_path=qc_img_path
                )
                qc_count += 1

        # 6. Hard Negative Background Paper Injections (~9%)
        num_negatives = max(1, int(round(len(split_vouchers) * negative_ratio)))
        logger.info(f"Injecting {num_negatives} hard negative background sheet tiles into split '{split_name}'...")

        for neg_idx in range(num_negatives):
            donor_v = split_vouchers[neg_idx % len(split_vouchers)]
            donor_img = cv2.imread(str(donor_v["image_path"]))
            if donor_img is None:
                continue

            neg_crop = extract_paper_background_crop(
                image_bgr=donor_img,
                annotations=[],
                crop_size=(1024, 1024),
                rng=rng
            )
            if neg_crop is None:
                continue

            neg_name = f"neg_paper_{donor_v['catalogNumber']}_{neg_idx:02d}"
            neg_img_path = output_dir / "images" / split_name / f"{neg_name}.jpg"
            neg_txt_path = output_dir / "labels" / split_name / f"{neg_name}.txt"

            cv2.imwrite(str(neg_img_path), neg_crop)
            with open(neg_txt_path, "w", encoding="utf-8") as f_neg:
                pass

            stats["negative_samples"] += 1
            stats[f"{split_name}_samples"] += 1
            manifest_records.append({
                "catalogNumber": neg_name,
                "split": split_name,
                "is_negative": True,
                "num_instances": 0,
                "image_path": str(neg_img_path),
                "label_path": str(neg_txt_path)
            })

    # 7. Write Ultralytics Dataset YAML
    config_dict = {
        "path": str(output_dir.resolve()).replace("\\", "/"),
        "train": "images/train",
        "val": "images/val",
        "test": "images/test",
        "names": {idx: name for idx, name in BOTANICAL_CLASSES.items()}
    }
    config_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_yaml_path, "w", encoding="utf-8") as f_yaml:
        yaml.dump(config_dict, f_yaml, default_flow_style=False, sort_keys=False)

    logger.info(f"Exported Ultralytics dataset configuration to: {config_yaml_path}")

    # 8. Export Manifest CSV & Summary JSON
    manifest_csv_path = DEFAULT_WORKSPACE / "data" / "tables" / "dataset_manifest.csv"
    manifest_csv_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(manifest_records).to_csv(manifest_csv_path, index=False)

    summary_json_path = DEFAULT_WORKSPACE / "outputs" / "reports" / "dataset_build_summary.json"
    summary_json_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_json_path, "w", encoding="utf-8") as f_rep:
        json.dump({
            "timestamp": datetime.datetime.now().isoformat(),
            "total_samples": len(manifest_records),
            "stats": stats,
            "classes": BOTANICAL_CLASSES,
            "config_yaml": str(config_yaml_path),
            "qc_plots_dir": str(qc_output_dir)
        }, f_rep, indent=2)

    logger.info("=" * 78)
    logger.info("BOTANICAL DATASET BUILD COMPLETED SUCCESSFULLY")
    logger.info(f"Total Dataset Samples: {len(manifest_records)}")
    logger.info(f"Train: {stats['train_samples']} | Val: {stats['val_samples']} | Test: {stats['test_samples']}")
    logger.info(f"Hard Negative Background Sheets: {stats['negative_samples']} ({stats['negative_samples']/max(1, len(manifest_records))*100:.1f}%)")
    logger.info(f"Total Instances Ingested: {stats['total_instances']}")
    logger.info(f"Class Breakdown: {stats['class_counts']}")
    logger.info("=" * 78)

    return {
        "status": "success",
        "total_samples": len(manifest_records),
        "stats": stats,
        "config_yaml": str(config_yaml_path),
        "qc_dir": str(qc_output_dir)
    }


def run_synthetic_test_suite() -> bool:
    """
    Generates synthetic botanical annotations and voucher sheets to test
    and verify the label parser, 7-class ontology, stratification, negative injection,
    and QC overlay generation.
    """
    logger.info("=" * 78)
    logger.info("RUNNING SYNTHETIC BOTANICAL DATASET VERIFICATION SUITE (--test)")
    logger.info("=" * 78)

    test_root = DEFAULT_WORKSPACE / "outputs" / "test_synthetic_dataset"
    test_raw_vouchers = test_root / "raw_vouchers"
    test_raw_annotations = test_root / "raw_annotations"
    test_output_dataset = test_root / "yolo_dataset"
    test_config_yaml = test_root / "dataset_config.yaml"
    test_qc_dir = test_root / "qc_overlays"

    if test_root.exists():
        shutil.rmtree(test_root)

    test_raw_vouchers.mkdir(parents=True, exist_ok=True)
    test_raw_annotations.mkdir(parents=True, exist_ok=True)

    test_vouchers = [
        {"cat": "NCU00099901", "inst": "NCU"},
        {"cat": "NCU00099902", "inst": "NCU"},
        {"cat": "GA00088801", "inst": "GA"},
        {"cat": "GA00088802", "inst": "GA"},
        {"cat": "US00077701", "inst": "US"},
        {"cat": "NY00066601", "inst": "NY"},
        {"cat": "BRIT00055501", "inst": "BRIT"},
    ]

    img_w, img_h = 1200, 1600

    for idx, v in enumerate(test_vouchers):
        cat = v["cat"]
        paper_bg = np.full((img_h, img_w, 3), (240, 245, 248), dtype=np.uint8)
        noise = np.random.normal(0, 4, (img_h, img_w, 3)).astype(np.float32)
        sheet_img = np.clip(paper_bg.astype(np.float32) + noise, 0, 255).astype(np.uint8)

        cv2.ellipse(sheet_img, (400, 1200), (80, 150), 30, 0, 360, (30, 120, 40), -1)
        cv2.line(sheet_img, (400, 1350), (450, 1480), (40, 140, 60), 6)
        cv2.ellipse(sheet_img, (580, 800), (40, 90), -45, 0, 360, (35, 130, 45), -1)
        cv2.line(sheet_img, (600, 400), (600, 1400), (50, 100, 50), 10)
        cv2.circle(sheet_img, (500, 1520), 45, (40, 60, 100), -1)
        cv2.circle(sheet_img, (450, 1400), 70, (20, 90, 30), -1)
        cv2.circle(sheet_img, (600, 380), 35, (40, 180, 220), -1)

        img_save_path = test_raw_vouchers / f"{cat}.jpg"
        cv2.imwrite(str(img_save_path), sheet_img)

        if idx % 2 == 0:
            txt_path = test_raw_annotations / f"{cat}.txt"
            with open(txt_path, "w", encoding="utf-8") as f_txt:
                f_txt.write("0 0.28 0.70 0.38 0.68 0.40 0.82 0.30 0.84\n")
                f_txt.write("leaf_petiole 0.33 0.84 0.35 0.84 0.38 0.93 0.36 0.93\n")
                f_txt.write("cauline_leaf 0.45 0.48 0.52 0.46 0.50 0.55 0.44 0.54\n")
                f_txt.write("3 0.49 0.25 0.51 0.25 0.51 0.88 0.49 0.88\n")
                f_txt.write("root_rhizome 0.38 0.92 0.45 0.92 0.45 0.98 0.38 0.98\n")
                f_txt.write("basal_rosette_clump 0.33 0.83 0.42 0.83 0.42 0.92 0.33 0.92\n")
                f_txt.write("capitulum 0.47 0.21 0.53 0.21 0.53 0.27 0.47 0.27\n")
        else:
            json_path = test_raw_annotations / f"{cat}.json"
            labelme_data = {
                "version": "5.0.1",
                "flags": {},
                "shapes": [
                    {"label": "basal_leaf_blade", "points": [[340, 1120], [460, 1090], [480, 1310], [360, 1340]], "shape_type": "polygon"},
                    {"label": "leaf_petiole", "points": [[395, 1345], [420, 1345], [455, 1490], [430, 1490]], "shape_type": "polygon"},
                    {"label": "cauline_leaf", "points": [[540, 770], [620, 740], [600, 880], [530, 860]], "shape_type": "polygon"},
                    {"label": "cauline_stem", "points": [[590, 400], [610, 400], [610, 1400], [590, 1400]], "shape_type": "polygon"},
                    {"label": "root_rhizome", "points": [[455, 1475], [545, 1475], [545, 1565], [455, 1565]], "shape_type": "polygon"},
                    {"label": "basal_rosette_clump", "points": [[380, 1330], [520, 1330], [520, 1470], [380, 1470]], "shape_type": "polygon"},
                    {"label": "capitulum", "points": [[565, 345], [635, 345], [635, 415], [565, 415]], "shape_type": "polygon"}
                ],
                "imagePath": f"{cat}.jpg",
                "imageHeight": img_h,
                "imageWidth": img_w
            }
            with open(json_path, "w", encoding="utf-8") as f_json:
                json.dump(labelme_data, f_json, indent=2)

    result = build_artifact_robust_dataset(
        raw_vouchers_dir=test_raw_vouchers,
        annotations_dir=test_raw_annotations,
        curated_csv_path=test_root / "non_existent.csv",
        output_dir=test_output_dataset,
        config_yaml_path=test_config_yaml,
        qc_output_dir=test_qc_dir,
        negative_ratio=0.09,
        train_ratio=0.70,
        val_ratio=0.15,
        test_ratio=0.15,
        min_instance_area=10.0,
        num_qc_plots=5,
        seed=42
    )

    if result.get("status") == "success":
        logger.info("SYNTHETIC TEST SUITE PASSED ALL VERIFICATIONS.")
        return True
    else:
        logger.error("SYNTHETIC TEST SUITE FAILED.")
        return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Construct Verified Artifact-Robust 7-Class YOLO Dataset for Packera dubia"
    )
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR, help="Path to raw voucher images")
    parser.add_argument("--annotations-dir", type=Path, default=DEFAULT_ANNOTATIONS_DIR, help="Path to raw annotations")
    parser.add_argument("--curated-csv", type=Path, default=DEFAULT_CURATED_CSV, help="Path to curated voucher metadata CSV")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Output YOLO dataset root directory")
    parser.add_argument("--config-yaml", type=Path, default=DEFAULT_CONFIG_PATH, help="Path to output dataset YAML")
    parser.add_argument("--qc-dir", type=Path, default=DEFAULT_QC_DIR, help="Path to save QC overlay images")
    parser.add_argument("--negative-ratio", type=float, default=0.09, help="Hard negative sheet ratio")
    parser.add_argument("--train-ratio", type=float, default=0.70, help="Train split ratio")
    parser.add_argument("--val-ratio", type=float, default=0.15, help="Val split ratio")
    parser.add_argument("--test-ratio", type=float, default=0.15, help="Test split ratio")
    parser.add_argument("--min-instance-area", type=float, default=50.0, help="Minimum polygon area in px^2")
    parser.add_argument("--num-qc-plots", type=int, default=25, help="Number of QC overlay figures to render")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of vouchers to process")
    parser.add_argument("--seed", type=int, default=42, help="Deterministic random seed")
    parser.add_argument("--test", action="store_true", help="Run synthetic verification test suite and exit")
    return parser.parse_args()
