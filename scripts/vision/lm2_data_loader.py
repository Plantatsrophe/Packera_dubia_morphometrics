#!/usr/bin/env python3
"""
scripts/vision/lm2_data_loader.py
=================================
Data ingestion and metadata linking module for LeafMachine2 (LM2) outputs.
Fast ingestion of ruler scale calibrations, leaf crops, masks, and Darwin Core metadata.
"""

from __future__ import annotations

import csv
import logging
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import cv2
import numpy as np
import pandas as pd

logger = logging.getLogger("LM2_PostProcessing")


@dataclass
class LeafCandidate:
    """Represents an extracted leaf instance candidate from LeafMachine2."""
    catalog_number: str
    leaf_id: int
    crop_image_path: Optional[Path] = None
    mask_image_path: Optional[Path] = None
    full_mask_path: Optional[Path] = None
    bbox: Tuple[int, int, int, int] = (0, 0, 0, 0)  # ymin, xmin, ymax, xmax on full sheet
    sheet_width: int = 3000
    sheet_height: int = 4000
    mask_array: Optional[np.ndarray] = None
    scale_mm_per_px: float = 0.04233  # default ~600 DPI fallback (1 inch / 600 px * 25.4)
    plant_individual_id: int = 0
    assigned_tier: str = "tier3"
    ucs_score: float = 0.0
    solidity: float = 0.0
    midrib_angle_deg: float = 0.0
    length_mm: float = 0.0
    width_mm: float = 0.0
    area_mm2: float = 0.0
    saved_mask_path: str = ""


def parse_bbox_from_filename(filename: str) -> Optional[Tuple[int, int, int, int]]:
    """Extracts (ymin, xmin, ymax, xmax) from LM2 filename format."""
    match = re.search(r"__(\d+)-(\d+)-(\d+)-(\d+)", filename)
    if match:
        return int(match.group(1)), int(match.group(2)), int(match.group(3)), int(match.group(4))
    return None


def load_ruler_calibrations(lm2_dir: Path) -> Dict[str, float]:
    """
    Fast loader for ruler scale calibrations across LM2 output directories.
    Reads master combined CSVs first, then batch files if available.
    """
    calibrations: Dict[str, float] = {}

    # 1. First check consolidated ruler CSV files
    master_ruler_files = list(lm2_dir.glob("**/Data/Ruler/*.csv")) + \
                         list(lm2_dir.glob("**/Data/Batch/Ruler/*.csv"))

    if master_ruler_files:
        for rf in master_ruler_files:
            try:
                df = pd.read_csv(rf)
                if df.empty:
                    continue
                for _, row in df.iterrows():
                    fn = str(row.get("filename", "")).strip()
                    if not fn or fn.lower() == "nan":
                        continue
                    cat = Path(fn).stem.split("__")[0]
                    success = str(row.get("ruler_success", "False")).strip().lower() == "true"
                    conv_mean = row.get("conversion_mean", np.nan)
                    pred_conv = row.get("predicted_conversion_factor_cm", np.nan)

                    conversion_val = None
                    try:
                        if pd.notna(conv_mean) and float(conv_mean) > 0:
                            conversion_val = float(conv_mean)
                        elif pd.notna(pred_conv) and float(pred_conv) > 0:
                            conversion_val = float(pred_conv)
                    except (ValueError, TypeError):
                        pass

                    if success and conversion_val and conversion_val > 0:
                        calibrations[cat] = 10.0 / conversion_val
                    elif conversion_val and conversion_val > 10.0 and cat not in calibrations:
                        calibrations[cat] = 10.0 / conversion_val
            except Exception as e:
                logger.debug(f"Error reading master ruler {rf}: {e}")

    logger.info(f"Loaded {len(calibrations)} specimen ruler scale calibrations from LM2.")
    return calibrations


def discover_lm2_candidates(
    lm2_dir: Path,
    curated_vouchers: pd.DataFrame,
    ruler_calibs: Dict[str, float],
    raw_images_dir: Optional[Path] = None
) -> Dict[str, List[LeafCandidate]]:
    """
    Discovers all candidate basal leaves from the latest/aggregated LM2 runs.
    Groups candidates by catalogNumber.
    """
    candidates_by_voucher: Dict[str, List[LeafCandidate]] = defaultdict(list)

    # Search for crop and mask directories
    leaf_crop_dirs = list(lm2_dir.glob("**/Plant_Components/Leaves_Whole")) + \
                     list(lm2_dir.glob("**/Plant_Components/Leaves_Partial")) + \
                     list(lm2_dir.glob("**/Keypoints/Oriented_Cropped_Leaves")) + \
                     list(lm2_dir.glob("**/crops/leaves"))

    mask_dirs = list(lm2_dir.glob("**/Plant_Components/Segmentation_Whole_Leaves")) + \
                list(lm2_dir.glob("**/Plant_Components/Segmentation_Masks_Color_Whole_Leaves")) + \
                list(lm2_dir.glob("**/Keypoints/Oriented_Masks")) + \
                list(lm2_dir.glob("**/crops/masks"))

    mask_lookup: Dict[str, Path] = {}
    for md in mask_dirs:
        for mf in md.glob("*.*"):
            if mf.suffix.lower() in [".jpg", ".jpeg", ".png", ".bmp", ".tif"]:
                mask_lookup[mf.stem] = mf

    leaf_counts: Dict[str, int] = defaultdict(int)

    for lcd in leaf_crop_dirs:
        for crop_path in lcd.glob("*.*"):
            if crop_path.suffix.lower() not in [".jpg", ".jpeg", ".png"]:
                continue

            stem = crop_path.stem
            cat_part = stem.split("__")[0].strip()
            bbox = parse_bbox_from_filename(stem) or (0, 0, 0, 0)
            matched_mask = mask_lookup.get(stem)

            scale = ruler_calibs.get(cat_part, 0.04233)
            leaf_counts[cat_part] += 1

            candidate = LeafCandidate(
                catalog_number=cat_part,
                leaf_id=leaf_counts[cat_part],
                crop_image_path=crop_path,
                mask_image_path=matched_mask,
                bbox=bbox,
                scale_mm_per_px=scale
            )
            candidates_by_voucher[cat_part].append(candidate)

    logger.info(f"Discovered leaf candidates for {len(candidates_by_voucher)} voucher sheets ({sum(len(v) for v in candidates_by_voucher.values())} total leaves).")
    return candidates_by_voucher


def load_and_preprocess_mask(candidate: LeafCandidate) -> Optional[np.ndarray]:
    """Loads and binarizes the leaf mask. Returns 8-bit binary mask (0 or 255)."""
    mask: Optional[np.ndarray] = None

    if candidate.mask_image_path and candidate.mask_image_path.exists():
        img = cv2.imread(str(candidate.mask_image_path), cv2.IMREAD_GRAYSCALE)
        if img is not None:
            _, mask = cv2.threshold(img, 10, 255, cv2.THRESH_BINARY)

    if mask is None and candidate.crop_image_path and candidate.crop_image_path.exists():
        img = cv2.imread(str(candidate.crop_image_path))
        if img is not None:
            gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
            blurred = cv2.GaussianBlur(gray, (5, 5), 0)
            _, mask = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

    if mask is not None:
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if contours:
            largest = max(contours, key=cv2.contourArea)
            clean_mask = np.zeros_like(mask)
            cv2.drawContours(clean_mask, [largest], -1, 255, thickness=cv2.FILLED)
            mask = clean_mask

    return mask
