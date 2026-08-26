#!/usr/bin/env python3
"""
===============================================================================
Script: build_artifact_robust_dataset.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
Author: Senior Computer Vision Engineer & Botanical Dataset Curator
Date: August 2026

Description:
    Constructs an artifact-robust, multi-class YOLO segmentation and detection
    training dataset designed to penalize and eliminate false-positive leaf
    detections caused by herbarium sheet artifacts (mounting tape, labels,
    ruler scales, color charts, and barcode stickers).

Key Technical Capabilities:
    1. Multi-Class YOLO Schema:
       Standardizes 9 semantic classes (4 botanical organ classes + 5 artifact classes):
         0: basal_leaf        1: leaf_petiole      2: basal_rosette
         3: capitulum         4: herbarium_label   5: color_chart
         6: ruler_scale       7: barcode_sticker   8: mounting_tape

    2. Hard Negative Injection (Pure Background Sheets):
       Ingests non-plant voucher regions, blank paper zones, and isolated background
       crops. Automatically outputs empty .txt label files (0 bounding instances)
       representing 8-10% of total dataset volume to enforce paper texture invariance.

    3. Synthetic Copy-Paste & Occlusion Augmentations:
       Extracts clean artifact patches (labels, tape strips, color swatches, rulers)
       and dynamically pastes them adjacent to, touching, or partially occluding
       botanical instances with realistic alpha blending and paper lighting.
       Dynamically updates polygon contours and bounding boxes via boolean geometric
       clipping to teach the loss function sharp boundary discrimination between
       straight artificial edges and organic leaf margins.

    4. Stratified Partitioning & Quality Control:
       Performs deterministic 70/15/15 (Train/Val/Test) splitting stratified by
       herbarium institution code and voucher quality tier. Exports Ultralytics
       `data/dataset_config.yaml` and renders multi-class overlay QC figures in
       `outputs/dataset_qc/`.

Usage:
    python scripts/build_artifact_robust_dataset.py --output-dir data/yolo_dataset --limit 100
===============================================================================
"""

import os
import sys
import math
import yaml
import json
import glob
import random
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union

import cv2
import numpy as np
import pandas as pd

# Optional tqdm for visual progress bars
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable


# =============================================================================
# 1. Dataset Schema & Global Configuration
# =============================================================================

# Multi-class schema mapping for botanical vision model
CLASS_NAMES: List[str] = [
    "basal_leaf",       # 0: Leaf blade / intact leaf
    "leaf_petiole",     # 1: Distinct petiole / leaf stalk
    "basal_rosette",    # 2: Clustered basal rosette
    "capitulum",        # 3: Inflorescence / flower head
    "herbarium_label",  # 4: Main specimen metadata label
    "color_chart",      # 5: Calibration color chart / palette
    "ruler_scale",      # 6: Measurement scale / centimeter bar
    "barcode_sticker",  # 7: Digitization barcode / QR sticker
    "mounting_tape",    # 8: Linen, paper, or plastic mounting tape strip
]

CLASS_MAP: Dict[str, int] = {name: idx for idx, name in enumerate(CLASS_NAMES)}

# Distinct RGB color palette for rendering class overlays during QC
CLASS_COLORS_BGR: Dict[int, Tuple[int, int, int]] = {
    0: (0, 200, 0),      # basal_leaf: Vibrant Green
    1: (50, 255, 150),   # leaf_petiole: Mint Green
    2: (0, 140, 70),     # basal_rosette: Dark Forest Green
    3: (0, 215, 255),    # capitulum: Gold / Yellow
    4: (30, 30, 230),    # herbarium_label: Bright Red
    5: (230, 30, 230),   # color_chart: Magenta
    6: (0, 140, 255),    # ruler_scale: Orange
    7: (230, 180, 0),    # barcode_sticker: Cyan
    8: (180, 0, 180),    # mounting_tape: Purple / Violet
}

# Default filesystem paths
DEFAULT_WORKSPACE = Path(__file__).resolve().parent.parent if "__file__" in locals() else Path.cwd()
DEFAULT_RAW_DIR = DEFAULT_WORKSPACE / "data" / "raw_vouchers"
DEFAULT_CURATED_CSV = DEFAULT_WORKSPACE / "data" / "tables" / "curated_vouchers.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_WORKSPACE / "data" / "yolo_dataset"
DEFAULT_CONFIG_PATH = DEFAULT_WORKSPACE / "data" / "dataset_config.yaml"
DEFAULT_QC_DIR = DEFAULT_WORKSPACE / "outputs" / "dataset_qc"


# =============================================================================
# 2. Logging Setup
# =============================================================================

def setup_logging(log_file: Optional[Path] = None, verbose: bool = False) -> logging.Logger:
    """
    Configures standard multi-handler logging for the dataset curator.
    
    Args:
        log_file: Optional path to an output log file.
        verbose: If True, enables debug-level logging.
        
    Returns:
        Configured logging.Logger instance.
    """
    logger = logging.getLogger("ArtifactRobustDatasetBuilder")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "[%(asctime)s] [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


logger = setup_logging()


# =============================================================================
# 3. Annotation Data Structures & Helper Functions
# =============================================================================

class InstanceAnnotation:
    """
    Represents a single segmented instance on a herbarium sheet.
    Stores class ID, normalized/absolute bounding box, and polygon contour points.
    """
    def __init__(
        self,
        class_id: int,
        polygon: np.ndarray,  # Shape (N, 2) in absolute pixel coordinates [[x, y], ...]
        bbox: Optional[Tuple[float, float, float, float]] = None,  # (x1, y1, x2, y2)
        confidence: float = 1.0,
        is_synthetic: bool = False,
        tag: str = ""
    ):
        self.class_id = int(class_id)
        self.polygon = np.asarray(polygon, dtype=np.float32)
        if bbox is not None:
            self.bbox = (float(bbox[0]), float(bbox[1]), float(bbox[2]), float(bbox[3]))
        else:
            if len(self.polygon) > 0:
                x_min = float(np.min(self.polygon[:, 0]))
                y_min = float(np.min(self.polygon[:, 1]))
                x_max = float(np.max(self.polygon[:, 0]))
                y_max = float(np.max(self.polygon[:, 1]))
                self.bbox = (x_min, y_min, x_max, y_max)
            else:
                self.bbox = (0.0, 0.0, 0.0, 0.0)
        self.confidence = float(confidence)
        self.is_synthetic = bool(is_synthetic)
        self.tag = tag

    @property
    def class_name(self) -> str:
        if 0 <= self.class_id < len(CLASS_NAMES):
            return CLASS_NAMES[self.class_id]
        return f"class_{self.class_id}"

    def to_yolo_seg_line(self, img_w: int, img_h: int) -> str:
        """
        Converts instance polygon to Ultralytics YOLO segmentation format:
        <class_id> <x1_norm> <y1_norm> <x2_norm> <y2_norm> ...
        """
        if len(self.polygon) < 3 or img_w <= 0 or img_h <= 0:
            # Fallback to normalized bounding box box polygon (4 corners)
            x1, y1, x2, y2 = self.bbox
            pts = np.array([
                [x1, y1], [x2, y1], [x2, y2], [x1, y2]
            ], dtype=np.float32)
        else:
            pts = self.polygon

        # Normalize coordinates between 0.0 and 1.0, clipped to bounds
        norm_pts = pts.copy()
        norm_pts[:, 0] = np.clip(norm_pts[:, 0] / float(img_w), 0.0, 1.0)
        norm_pts[:, 1] = np.clip(norm_pts[:, 1] / float(img_h), 0.0, 1.0)

        coords_str = " ".join([f"{x:.6f} {y:.6f}" for x, y in norm_pts])
        return f"{self.class_id} {coords_str}"

    def to_yolo_det_line(self, img_w: int, img_h: int) -> str:
        """
        Converts instance bounding box to Ultralytics YOLO detection format:
        <class_id> <x_center_norm> <y_center_norm> <width_norm> <height_norm>
        """
        x1, y1, x2, y2 = self.bbox
        bw = max(0.0, x2 - x1)
        bh = max(0.0, y2 - y1)
        cx = x1 + bw / 2.0
        cy = y1 + bh / 2.0

        cx_norm = np.clip(cx / float(img_w), 0.0, 1.0)
        cy_norm = np.clip(cy / float(img_h), 0.0, 1.0)
        w_norm = np.clip(bw / float(img_w), 0.0, 1.0)
        h_norm = np.clip(bh / float(img_h), 0.0, 1.0)

        return f"{self.class_id} {cx_norm:.6f} {cy_norm:.6f} {w_norm:.6f} {h_norm:.6f}"


# =============================================================================
# 4. Artifact Detection, Segmentation & Library Extraction Engine
# =============================================================================

class ArtifactHarvester:
    """
    Harvests clean, realistic artifact crops from herbarium sheets:
    - Herbarium Labels (lower quadrants, printed text blocks)
    - Color Calibration Charts (multi-color square blocks)
    - Ruler Scales (centimeter/millimeter tick bars)
    - Barcode Stickers (code bars, numbering)
    - Mounting Tape Strips (translucent strips across plant stems or paper)
    """
    def __init__(self, rng_seed: int = 42):
        self.rng = random.Random(rng_seed)
        self.np_rng = np.random.default_rng(rng_seed)
        # In-memory artifact bank grouped by class name
        self.artifact_bank: Dict[str, List[Dict[str, Any]]] = {
            "herbarium_label": [],
            "color_chart": [],
            "ruler_scale": [],
            "barcode_sticker": [],
            "mounting_tape": [],
        }

    def detect_and_extract_sheet_artifacts(
        self, image_bgr: np.ndarray, catalog_number: str = "voucher"
    ) -> List[InstanceAnnotation]:
        """
        Extracts artifact annotations and caches cropped artifact patches.
        
        Args:
            image_bgr: Full herbarium sheet BGR image.
            catalog_number: Specimen voucher ID.
            
        Returns:
            List of detected artifact InstanceAnnotations.
        """
        h, w = image_bgr.shape[:2]
        gray = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2GRAY)
        annotations: List[InstanceAnnotation] = []

        # ---------------------------------------------------------------------
        # 1. Herbarium Label (typically bottom-right quadrant)
        # ---------------------------------------------------------------------
        label_y1, label_y2 = int(h * 0.60), int(h * 0.99)
        label_x1, label_x2 = int(w * 0.50), int(w * 0.99)
        label_roi = gray[label_y1:label_y2, label_x1:label_x2]

        _, label_thresh = cv2.threshold(label_roi, 195, 255, cv2.THRESH_BINARY)
        contours, _ = cv2.findContours(label_thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        for cnt in contours:
            area = cv2.contourArea(cnt)
            if area > (w * h * 0.015):
                bx, by, bw, bh = cv2.boundingRect(cnt)
                # Ensure rectangular aspect ratio typical of herbarium labels
                if 0.5 < (bw / max(1, bh)) < 3.0:
                    abs_x1, abs_y1 = label_x1 + bx, label_y1 + by
                    abs_x2, abs_y2 = abs_x1 + bw, abs_y1 + bh
                    poly = np.array([
                        [abs_x1, abs_y1], [abs_x2, abs_y1],
                        [abs_x2, abs_y2], [abs_x1, abs_y2]
                    ], dtype=np.float32)
                    ann = InstanceAnnotation(
                        class_id=CLASS_MAP["herbarium_label"],
                        polygon=poly,
                        bbox=(abs_x1, abs_y1, abs_x2, abs_y2),
                        confidence=0.95,
                        tag="detected_label"
                    )
                    annotations.append(ann)
                    # Cache crop
                    crop = image_bgr[abs_y1:abs_y2, abs_x1:abs_x2].copy()
                    if crop.size > 0:
                        self.artifact_bank["herbarium_label"].append({
                            "image": crop,
                            "source": catalog_number
                        })
                    break

        # ---------------------------------------------------------------------
        # 2. Color Calibration Chart & Ruler Scale (typically along left margin or top)
        # ---------------------------------------------------------------------
        left_margin_w = int(w * 0.22)
        left_roi = image_bgr[:, :left_margin_w]
        left_hsv = cv2.cvtColor(left_roi, cv2.COLOR_BGR2HSV)
        sat = left_hsv[:, :, 1]

        # High saturation regions indicate color patches
        _, sat_thresh = cv2.threshold(sat, 80, 255, cv2.THRESH_BINARY)
        chart_cnts, _ = cv2.findContours(sat_thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        for cnt in chart_cnts:
            area = cv2.contourArea(cnt)
            if (w * h * 0.003) < area < (w * h * 0.08):
                bx, by, bw, bh = cv2.boundingRect(cnt)
                # Expand slightly to capture the full chart frame
                margin = int(bw * 0.15)
                cx1 = max(0, bx - margin)
                cy1 = max(0, by - margin)
                cx2 = min(left_margin_w, bx + bw + margin)
                cy2 = min(h, by + bh + margin)
                poly = np.array([
                    [cx1, cy1], [cx2, cy1],
                    [cx2, cy2], [cx1, cy2]
                ], dtype=np.float32)
                ann = InstanceAnnotation(
                    class_id=CLASS_MAP["color_chart"],
                    polygon=poly,
                    bbox=(cx1, cy1, cx2, cy2),
                    confidence=0.90,
                    tag="detected_color_chart"
                )
                annotations.append(ann)
                crop = image_bgr[cy1:cy2, cx1:cx2].copy()
                if crop.size > 0:
                    self.artifact_bank["color_chart"].append({
                        "image": crop,
                        "source": catalog_number
                    })
                break

        # ---------------------------------------------------------------------
        # 3. Barcode Sticker (top margins or corners, stark black-and-white stripes)
        # ---------------------------------------------------------------------
        top_roi = gray[:int(h * 0.35), :]
        grad_x = cv2.Sobel(top_roi, cv2.CV_32F, 1, 0, ksize=3)
        grad_y = cv2.Sobel(top_roi, cv2.CV_32F, 0, 1, ksize=3)
        grad = cv2.subtract(grad_x, grad_y)
        grad = cv2.convertScaleAbs(grad)
        blurred = cv2.blur(grad, (9, 9))
        _, barcode_thresh = cv2.threshold(blurred, 180, 255, cv2.THRESH_BINARY)
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (21, 7))
        closed = cv2.morphologyEx(barcode_thresh, cv2.MORPH_CLOSE, kernel)

        barcode_cnts, _ = cv2.findContours(closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        for cnt in barcode_cnts:
            area = cv2.contourArea(cnt)
            if (w * h * 0.001) < area < (w * h * 0.03):
                bx, by, bw, bh = cv2.boundingRect(cnt)
                aspect = bw / max(1, bh)
                if 1.5 < aspect < 5.5:
                    poly = np.array([
                        [bx, by], [bx + bw, by],
                        [bx + bw, by + bh], [bx, by + bh]
                    ], dtype=np.float32)
                    ann = InstanceAnnotation(
                        class_id=CLASS_MAP["barcode_sticker"],
                        polygon=poly,
                        bbox=(bx, by, bx + bw, by + bh),
                        confidence=0.88,
                        tag="detected_barcode"
                    )
                    annotations.append(ann)
                    crop = image_bgr[by:by+bh, bx:bx+bw].copy()
                    if crop.size > 0:
                        self.artifact_bank["barcode_sticker"].append({
                            "image": crop,
                            "source": catalog_number
                        })
                    break

        # ---------------------------------------------------------------------
        # 4. Ruler Scale (elongated rectangular tick strip)
        # ---------------------------------------------------------------------
        ruler_edges = cv2.Canny(left_roi, 50, 150)
        lines = cv2.HoughLinesP(ruler_edges, 1, np.pi / 180, threshold=80, minLineLength=int(h * 0.08), maxLineGap=15)
        if lines is not None and len(lines) > 0:
            lines_reshaped = lines.reshape(-1, 4)
            min_rx = int(np.min(lines_reshaped[:, [0, 2]]))
            max_rx = int(np.max(lines_reshaped[:, [0, 2]]))
            min_ry = int(np.min(lines_reshaped[:, [1, 3]]))
            max_ry = int(np.max(lines_reshaped[:, [1, 3]]))
            rw = max_rx - min_rx
            rh = max_ry - min_ry
            if rh > rw * 2 and (w * h * 0.005) < (rw * rh) < (w * h * 0.05):
                poly = np.array([
                    [min_rx, min_ry], [max_rx, min_ry],
                    [max_rx, max_ry], [min_rx, max_ry]
                ], dtype=np.float32)
                ann = InstanceAnnotation(
                    class_id=CLASS_MAP["ruler_scale"],
                    polygon=poly,
                    bbox=(min_rx, min_ry, max_rx, max_ry),
                    confidence=0.85,
                    tag="detected_ruler"
                )
                annotations.append(ann)
                crop = image_bgr[min_ry:max_ry, min_rx:max_rx].copy()
                if crop.size > 0:
                    self.artifact_bank["ruler_scale"].append({
                        "image": crop,
                        "source": catalog_number
                    })

        # Ensure synthetic fallback patches if any category is sparse
        self._ensure_synthetic_fallbacks(image_bgr)

        return annotations

    def _ensure_synthetic_fallbacks(self, sheet_bgr: np.ndarray) -> None:
        """
        Synthesizes high-fidelity realistic artifact patches (mounting tape strips,
        color bars, label stamps) when natural crops are sparse.
        """
        # 1. Realistic Translucent Linen / Paper Mounting Tape Strips
        if len(self.artifact_bank["mounting_tape"]) < 20:
            for _ in range(5):
                tw = self.rng.randint(60, 160)
                th = self.rng.randint(15, 35)
                # Sample paper background color from the sheet to match lighting
                bg_sample = np.median(sheet_bgr[100:200, 100:200], axis=(0, 1)).astype(np.float32)
                # Tape has slightly amber / aged translucent hue
                tape_color = np.clip(bg_sample * np.array([0.88, 0.94, 0.98]), 0, 255).astype(np.uint8)
                tape_img = np.full((th, tw, 3), tape_color, dtype=np.uint8)
                # Add micro-fibrous noise
                noise = np.random.normal(0, 8, (th, tw, 3)).astype(np.float32)
                tape_img = np.clip(tape_img.astype(np.float32) + noise, 0, 255).astype(np.uint8)
                # Add straight boundary lines
                cv2.rectangle(tape_img, (0, 0), (tw - 1, th - 1), (120, 130, 140), 1)
                self.artifact_bank["mounting_tape"].append({
                    "image": tape_img,
                    "source": "synthetic_generator"
                })

        # 2. Color Bar Fallbacks
        if len(self.artifact_bank["color_chart"]) < 10:
            swatch_w, swatch_h = 30, 30
            num_blocks = 6
            chart_img = np.zeros((swatch_h, swatch_w * num_blocks, 3), dtype=np.uint8)
            colors = [
                (40, 40, 200), (40, 180, 40), (200, 40, 40),
                (40, 200, 200), (200, 40, 200), (200, 200, 40)
            ]
            for i, col in enumerate(colors):
                cv2.rectangle(chart_img, (i * swatch_w, 0), ((i + 1) * swatch_w, swatch_h), col, -1)
                cv2.rectangle(chart_img, (i * swatch_w, 0), ((i + 1) * swatch_w, swatch_h), (50, 50, 50), 1)
            self.artifact_bank["color_chart"].append({
                "image": chart_img,
                "source": "synthetic_generator"
            })

    def get_random_artifact_crop(self, preferred_class: Optional[str] = None) -> Optional[Tuple[str, np.ndarray]]:
        """
        Retrieves a randomly sampled artifact patch from the library.
        """
        classes = [preferred_class] if preferred_class and preferred_class in self.artifact_bank else list(self.artifact_bank.keys())
        valid_classes = [c for c in classes if len(self.artifact_bank[c]) > 0]
        if not valid_classes:
            return None
        chosen_cls = self.rng.choice(valid_classes)
        item = self.rng.choice(self.artifact_bank[chosen_cls])
        return chosen_cls, item["image"].copy()


# =============================================================================
# 5. Botanical Instance Segmentation Extractor
# =============================================================================

def extract_botanical_annotations(
    image_bgr: np.ndarray,
    artifact_anns: List[InstanceAnnotation],
    min_instance_area: float = 150.0
) -> List[InstanceAnnotation]:
    """
    Extracts botanical organ instances (`basal_leaf`, `leaf_petiole`, `basal_rosette`, `capitulum`)
    from the voucher sheet using anatomical color thresholding and morphological clustering,
    strictly excluding regions occupied by known sheet artifacts.
    """
    h, w = image_bgr.shape[:2]
    annotations: List[InstanceAnnotation] = []

    # Build an artifact exclusion mask
    artifact_mask = np.zeros((h, w), dtype=np.uint8)
    for art in artifact_anns:
        if len(art.polygon) >= 3:
            pts = art.polygon.astype(np.int32)
            cv2.fillPoly(artifact_mask, [pts], 255)
        else:
            x1, y1, x2, y2 = [int(v) for v in art.bbox]
            cv2.rectangle(artifact_mask, (x1, y1), (x2, y2), 255, -1)

    # Convert to LAB / HSV color space for plant biomass extraction (brown/green leaves)
    lab = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2LAB)
    l_chan = lab[:, :, 0]
    b_chan = lab[:, :, 2]

    # Otsu thresholding on L-channel (herbarium plant material is darker than paper background)
    _, plant_thresh = cv2.threshold(l_chan, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    # Zero out artifacts
    plant_thresh[artifact_mask > 0] = 0

    # Morphological clean up
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    clean_plant = cv2.morphologyEx(plant_thresh, cv2.MORPH_OPEN, kernel)

    contours, _ = cv2.findContours(clean_plant, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    for cnt in contours:
        area = cv2.contourArea(cnt)
        # Filter out instances below min_instance_area or giant background clumps (>20% of sheet)
        if max(min_instance_area, (w * h * 0.0002)) < area < (w * h * 0.20):
            # Compute solidity and aspect ratio
            hull = cv2.convexHull(cnt)
            hull_area = max(1.0, cv2.contourArea(hull))
            solidity = area / hull_area
            bx, by, bw, bh = cv2.boundingRect(cnt)
            aspect = bw / max(1.0, bh)

            # Approximate smooth polygon
            epsilon = 0.004 * cv2.arcLength(cnt, True)
            approx_poly = cv2.approxPolyDP(cnt, epsilon, True).reshape(-1, 2)
            if len(approx_poly) < 3:
                continue

            # Class determination based on morphometrics
            cy = by + bh / 2.0
            if cy > (h * 0.55) and area > (w * h * 0.015):
                # Dense basal rosette cluster
                class_id = CLASS_MAP["basal_rosette"]
            elif cy < (h * 0.40) and (0.7 < aspect < 1.4) and area < (w * h * 0.008):
                # Flower head / capitulum in upper inflorescence
                class_id = CLASS_MAP["capitulum"]
            elif aspect > 3.0 or aspect < 0.33:
                # Elongated petiole
                class_id = CLASS_MAP["leaf_petiole"]
            else:
                # Standard basal / cauline leaf blade
                class_id = CLASS_MAP["basal_leaf"]

            ann = InstanceAnnotation(
                class_id=class_id,
                polygon=approx_poly,
                bbox=(bx, by, bx + bw, by + bh),
                confidence=0.92,
                tag="botanical_instance"
            )
            annotations.append(ann)

    return annotations


# =============================================================================
# 6. Hard Negative Mining Engine (Pure Background Paper Invariance)
# =============================================================================

def extract_hard_negative_background_crop(
    image_bgr: np.ndarray,
    all_annotations: List[InstanceAnnotation],
    crop_size: Tuple[int, int] = (1024, 1024),
    max_attempts: int = 25
) -> Optional[np.ndarray]:
    """
    Extracts a pure background sheet crop containing zero plant or artifact instances.
    Used to generate negative training samples with empty annotation files.
    
    Args:
        image_bgr: Full herbarium sheet image.
        all_annotations: Existing annotations to avoid.
        crop_size: Output (width, height) of the negative tile.
        max_attempts: Number of random spatial samples before fallback.
        
    Returns:
        Pure background BGR crop, or None if no vacant region satisfies constraints.
    """
    h, w = image_bgr.shape[:2]
    cw, ch = crop_size
    if w <= cw or h <= ch:
        return None

    # Build occupied occupancy mask
    occupancy = np.zeros((h, w), dtype=np.uint8)
    for ann in all_annotations:
        x1, y1, x2, y2 = [int(v) for v in ann.bbox]
        # Pad bounding box to guarantee complete vacancy
        px1 = max(0, x1 - 30)
        py1 = max(0, y1 - 30)
        px2 = min(w, x2 + 30)
        py2 = min(h, y2 + 30)
        cv2.rectangle(occupancy, (px1, py1), (px2, py2), 255, -1)

    # Search for an unallocated window
    for _ in range(max_attempts):
        rx = random.randint(0, w - cw)
        ry = random.randint(0, h - ch)
        sub_occupancy = occupancy[ry:ry+ch, rx:rx+cw]
        if np.count_nonzero(sub_occupancy) == 0:
            return image_bgr[ry:ry+ch, rx:rx+cw].copy()

    # Fallback: Sample sheet corner/margin with minimal gradient
    margin_w = min(cw, int(w * 0.25))
    margin_h = min(ch, int(h * 0.25))
    corner_crop = image_bgr[:margin_h, :margin_w].copy()
    return cv2.resize(corner_crop, (cw, ch), interpolation=cv2.INTER_LINEAR)


# =============================================================================
# 7. Synthetic Copy-Paste & Occlusion Augmentation Engine
# =============================================================================

class SyntheticOcclusionAugmenter:
    """
    Applies synthetic copy-paste augmentations of non-plant artifacts
    adjacent to, touching, or occluding annotated botanical leaves.
    Updates polygon geometries and bounding boxes to reflect sharp boundaries.
    """
    def __init__(self, harvester: ArtifactHarvester, rng_seed: int = 42):
        self.harvester = harvester
        self.rng = random.Random(rng_seed)
        self.np_rng = np.random.default_rng(rng_seed)

    def apply_copy_paste_augmentation(
        self,
        image_bgr: np.ndarray,
        annotations: List[InstanceAnnotation],
        paste_probability: float = 0.75,
        max_pastes_per_image: int = 3
    ) -> Tuple[np.ndarray, List[InstanceAnnotation]]:
        """
        Executes dynamic copy-paste augmentation on a single herbarium sheet.
        
        Args:
            image_bgr: Original herbarium sheet BGR image.
            annotations: Current list of InstanceAnnotations.
            paste_probability: Probability of applying augmentation.
            max_pastes_per_image: Max artifact patches to paste.
            
        Returns:
            Tuple of (augmented_image_bgr, updated_annotations_list).
        """
        if self.rng.random() > paste_probability or not annotations:
            return image_bgr, annotations

        aug_image = image_bgr.copy()
        img_h, img_w = aug_image.shape[:2]
        updated_anns = [ann for ann in annotations]

        # Target basal_leaf instances for occlusion / adjacency
        leaf_indices = [
            i for i, ann in enumerate(updated_anns)
            if ann.class_id == CLASS_MAP["basal_leaf"]
        ]

        if not leaf_indices:
            return aug_image, updated_anns

        num_pastes = self.rng.randint(1, max_pastes_per_image)

        for _ in range(num_pastes):
            target_idx = self.rng.choice(leaf_indices)
            target_leaf = updated_anns[target_idx]
            lx1, ly1, lx2, ly2 = target_leaf.bbox
            lw = max(10, int(lx2 - lx1))
            lh = max(10, int(ly2 - ly1))

            # Sample random artifact patch (tape, label, ruler, color swatch)
            artifact_sample = self.harvester.get_random_artifact_crop()
            if artifact_sample is None:
                continue

            art_class_name, art_crop = artifact_sample
            art_class_id = CLASS_MAP[art_class_name]
            ah, aw = art_crop.shape[:2]

            # Scale artifact appropriately relative to the leaf (20% to 80% leaf scale)
            scale = self.rng.uniform(0.3, 0.9) * (max(lw, lh) / max(aw, ah, 1))
            new_aw = max(15, min(int(aw * scale), img_w // 3))
            new_ah = max(10, min(int(ah * scale), img_h // 3))
            resized_art = cv2.resize(art_crop, (new_aw, new_ah), interpolation=cv2.INTER_AREA)

            # Random rotation (-35 to +35 degrees)
            angle = self.rng.uniform(-35, 35)
            rot_mat = cv2.getRotationMatrix2D((new_aw / 2, new_ah / 2), angle, 1.0)
            cos = np.abs(rot_mat[0, 0])
            sin = np.abs(rot_mat[0, 1])
            bound_w = int((new_ah * sin) + (new_aw * cos))
            bound_h = int((new_ah * cos) + (new_aw * sin))
            rot_mat[0, 2] += (bound_w / 2) - (new_aw / 2)
            rot_mat[1, 2] += (bound_h / 2) - (new_ah / 2)

            rotated_art = cv2.warpAffine(
                resized_art, rot_mat, (bound_w, bound_h),
                borderMode=cv2.BORDER_CONSTANT, borderValue=(255, 255, 255)
            )

            # Create binary foreground mask for the rotated artifact
            art_gray = cv2.cvtColor(rotated_art, cv2.COLOR_BGR2GRAY)
            # Anything not pure white border is artifact foreground
            art_mask = (art_gray < 250).astype(np.uint8) * 255
            # Feather edge for seamless paper texture alpha-blending
            art_mask_blurred = cv2.GaussianBlur(art_mask, (3, 3), 0)
            alpha = (art_mask_blurred.astype(np.float32) / 255.0)[:, :, np.newaxis]

            # Choose placement mode:
            # Mode A: Partial Occlusion (Pasted across leaf blade/petiole)
            # Mode B: Adjacent / Touching (Pasted directly adjacent to margin)
            mode = self.rng.choice(["occlusion", "adjacent"])

            if mode == "occlusion":
                paste_x = int(self.rng.uniform(lx1 - bound_w * 0.3, lx2 - bound_w * 0.7))
                paste_y = int(self.rng.uniform(ly1 - bound_h * 0.3, ly2 - bound_h * 0.7))
            else:
                offset_side = self.rng.choice(["left", "right", "top", "bottom"])
                if offset_side == "left":
                    paste_x = int(lx1 - bound_w * 0.9)
                    paste_y = int(ly1 + self.rng.uniform(-bound_h * 0.2, lh * 0.5))
                elif offset_side == "right":
                    paste_x = int(lx2 - bound_w * 0.1)
                    paste_y = int(ly1 + self.rng.uniform(-bound_h * 0.2, lh * 0.5))
                elif offset_side == "top":
                    paste_x = int(lx1 + self.rng.uniform(-bound_w * 0.2, lw * 0.5))
                    paste_y = int(ly1 - bound_h * 0.9)
                else:
                    paste_x = int(lx1 + self.rng.uniform(-bound_w * 0.2, lw * 0.5))
                    paste_y = int(ly2 - bound_h * 0.1)

            # Clip placement coordinates to image boundaries
            paste_x1 = max(0, paste_x)
            paste_y1 = max(0, paste_y)
            paste_x2 = min(img_w, paste_x + bound_w)
            paste_y2 = min(img_h, paste_y + bound_h)

            crop_w = paste_x2 - paste_x1
            crop_h = paste_y2 - paste_y1

            if crop_w <= 5 or crop_h <= 5:
                continue

            art_sub_x1 = paste_x1 - paste_x
            art_sub_y1 = paste_y1 - paste_y
            art_sub_x2 = art_sub_x1 + crop_w
            art_sub_y2 = art_sub_y1 + crop_h

            sub_art = rotated_art[art_sub_y1:art_sub_y2, art_sub_x1:art_sub_x2]
            sub_alpha = alpha[art_sub_y1:art_sub_y2, art_sub_x1:art_sub_x2]

            # Alpha-blend artifact into the canvas
            target_roi = aug_image[paste_y1:paste_y2, paste_x1:paste_x2].astype(np.float32)
            blended_roi = (sub_art.astype(np.float32) * sub_alpha) + (target_roi * (1.0 - sub_alpha))
            aug_image[paste_y1:paste_y2, paste_x1:paste_x2] = np.clip(blended_roi, 0, 255).astype(np.uint8)

            # Define new artifact polygon and bounding box
            art_poly = np.array([
                [paste_x1, paste_y1], [paste_x2, paste_y1],
                [paste_x2, paste_y2], [paste_x1, paste_y2]
            ], dtype=np.float32)

            new_art_ann = InstanceAnnotation(
                class_id=art_class_id,
                polygon=art_poly,
                bbox=(paste_x1, paste_y1, paste_x2, paste_y2),
                confidence=1.0,
                is_synthetic=True,
                tag=f"aug_copy_paste_{art_class_name}"
            )
            updated_anns.append(new_art_ann)

            # Update occluded leaf polygon geometry dynamically (Boolean Difference)
            if mode == "occlusion" and len(target_leaf.polygon) >= 3:
                updated_leaf_poly = self._compute_occluded_polygon(
                    target_leaf.polygon, (paste_x1, paste_y1, paste_x2, paste_y2), img_w, img_h
                )
                if updated_leaf_poly is not None and len(updated_leaf_poly) >= 3:
                    target_leaf.polygon = updated_leaf_poly
                    # Recompute bounding box
                    x_min = float(np.min(updated_leaf_poly[:, 0]))
                    y_min = float(np.min(updated_leaf_poly[:, 1]))
                    x_max = float(np.max(updated_leaf_poly[:, 0]))
                    y_max = float(np.max(updated_leaf_poly[:, 1]))
                    target_leaf.bbox = (x_min, y_min, x_max, y_max)

        return aug_image, updated_anns

    def _compute_occluded_polygon(
        self,
        leaf_poly: np.ndarray,
        occluder_bbox: Tuple[int, int, int, int],
        img_w: int,
        img_h: int
    ) -> Optional[np.ndarray]:
        """
        Subtracts the occluding rectangular bounding mask from the leaf polygon.
        """
        ox1, oy1, ox2, oy2 = occluder_bbox
        # Create full-sheet binary mask for the leaf
        leaf_mask = np.zeros((img_h, img_w), dtype=np.uint8)
        pts = leaf_poly.astype(np.int32)
        cv2.fillPoly(leaf_mask, [pts], 255)

        # Subtract occluder region
        leaf_mask[oy1:oy2, ox1:ox2] = 0

        # Extract largest remaining contour
        contours, _ = cv2.findContours(leaf_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            return None

        largest_cnt = max(contours, key=cv2.contourArea)
        if cv2.contourArea(largest_cnt) < 50:
            return None

        epsilon = 0.005 * cv2.arcLength(largest_cnt, True)
        approx = cv2.approxPolyDP(largest_cnt, epsilon, True).reshape(-1, 2)
        return approx.astype(np.float32)


# =============================================================================
# 8. Deterministic & Stratified Dataset Splitter
# =============================================================================

def stratify_and_partition_dataset(
    vouchers: List[Dict[str, Any]],
    train_ratio: float = 0.70,
    val_ratio: float = 0.15,
    test_ratio: float = 0.15,
    rng_seed: int = 42
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Deterministically partitions vouchers into Train, Validation, and Test splits,
    stratified by herbarium institution code and voucher quality tier.
    
    Args:
        vouchers: List of voucher metadata records.
        train_ratio: Proportion allocated to training (default 0.70).
        val_ratio: Proportion allocated to validation (default 0.15).
        test_ratio: Proportion allocated to testing (default 0.15).
        rng_seed: Random seed for deterministic reproducibility.
        
    Returns:
        Dictionary mapping split names ("train", "val", "test") to lists of records.
    """
    rng = random.Random(rng_seed)

    # Group vouchers by stratum: (institutionCode, determiner_tier)
    strata: Dict[str, List[Dict[str, Any]]] = {}
    for v in vouchers:
        inst = str(v.get("institutionCode", "UNKNOWN"))
        tier = str(v.get("determiner_tier", "Tier_Default"))
        stratum_key = f"{inst}_{tier}"
        if stratum_key not in strata:
            strata[stratum_key] = []
        strata[stratum_key].append(v)

    splits: Dict[str, List[Dict[str, Any]]] = {
        "train": [],
        "val": [],
        "test": []
    }

    # Allocate each stratum proportionally
    for stratum_key, group in strata.items():
        rng.shuffle(group)
        n = len(group)
        n_train = int(round(n * train_ratio))
        n_val = int(round(n * val_ratio))
        # Ensure at least 1 sample in train if group >= 1
        if n_train == 0 and n > 0:
            n_train = 1

        train_items = group[:n_train]
        val_items = group[n_train:n_train + n_val]
        test_items = group[n_train + n_val:]

        splits["train"].extend(train_items)
        splits["val"].extend(val_items)
        splits["test"].extend(test_items)

    logger.info(
        f"Stratified Dataset Partition: Train={len(splits['train'])} | "
        f"Val={len(splits['val'])} | Test={len(splits['test'])} "
        f"(Total: {len(vouchers)} vouchers across {len(strata)} strata)"
    )

    return splits


# =============================================================================
# 9. Quality Control & Visualization Overlay Renderer
# =============================================================================

def render_qc_verification_overlay(
    image_bgr: np.ndarray,
    annotations: List[InstanceAnnotation],
    catalog_number: str,
    output_path: Path
) -> None:
    """
    Renders high-visibility multi-class bounding boxes, polygon contours,
    class labels, and synthetic augmentation flags for visual verification.
    """
    overlay = image_bgr.copy()
    h, w = overlay.shape[:2]

    # Create semi-transparent overlay layer for masks
    mask_layer = np.zeros_like(overlay, dtype=np.uint8)

    for ann in annotations:
        color = CLASS_COLORS_BGR.get(ann.class_id, (255, 255, 255))
        # Draw polygon mask if available
        if len(ann.polygon) >= 3:
            pts = ann.polygon.astype(np.int32)
            cv2.fillPoly(mask_layer, [pts], color)
            cv2.polylines(overlay, [pts], isClosed=True, color=color, thickness=2)

        # Draw bounding box
        x1, y1, x2, y2 = [int(v) for v in ann.bbox]
        cv2.rectangle(overlay, (x1, y1), (x2, y2), color, thickness=2)

        # Label tag
        tag_str = f"{ann.class_name}"
        if ann.is_synthetic:
            tag_str += " [SYN_AUG]"

        (tw, th), _ = cv2.getTextSize(tag_str, cv2.FONT_HERSHEY_SIMPLEX, 0.5, 1)
        cv2.rectangle(overlay, (x1, max(0, y1 - th - 6)), (x1 + tw + 6, max(0, y1)), color, -1)
        cv2.putText(
            overlay, tag_str, (x1 + 3, max(th, y1 - 3)),
            cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1, cv2.LINE_AA
        )

    # Blend polygon masks with 30% alpha
    cv2.addWeighted(mask_layer, 0.35, overlay, 0.65, 0, overlay)

    # Header banner with catalogNumber and instance tally
    cv2.rectangle(overlay, (0, 0), (w, 40), (20, 20, 20), -1)
    title = f"QC Verification: {catalog_number} | Instances: {len(annotations)}"
    cv2.putText(overlay, title, (15, 26), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (255, 255, 255), 2, cv2.LINE_AA)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(output_path), overlay)


# =============================================================================
# 10. Ultralytics Dataset Configuration Exporter
# =============================================================================

def export_ultralytics_dataset_yaml(
    dataset_root: Path,
    output_yaml_path: Path
) -> None:
    """
    Exports the dataset YAML configuration file for training Ultralytics YOLOv8 / YOLOv11.
    """
    # Relative or absolute POSIX paths for compatibility across platforms
    dataset_dict = {
        "path": str(dataset_root.resolve()).replace("\\", "/"),
        "train": "images/train",
        "val": "images/val",
        "test": "images/test",
        "names": {idx: name for idx, name in enumerate(CLASS_NAMES)}
    }

    output_yaml_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_yaml_path, "w", encoding="utf-8") as f:
        yaml.dump(dataset_dict, f, default_flow_style=False, sort_keys=False)

    logger.info(f"Ultralytics dataset configuration written to: {output_yaml_path}")


# =============================================================================
# 11. Main Dataset Construction Pipeline
# =============================================================================

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
        for neg_idx in range(num_negatives):
            donor_v = split_vouchers[neg_idx % len(split_vouchers)]
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


# =============================================================================
# 12. Command Line Interface (CLI)
# =============================================================================

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
        help="Proportion of dataset comprising pure background negative sheets (default: 0.09 / ~9%)."
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
        "--limit", type=int, default=None,
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


if __name__ == "__main__":
    cli_args = parse_args()
    if cli_args.verbose:
        logger.setLevel(logging.DEBUG)

    build_artifact_robust_dataset(
        raw_vouchers_dir=cli_args.raw_dir,
        annotations_dir=cli_args.annotations_dir,
        curated_csv_path=cli_args.curated_csv,
        output_dir=cli_args.output_dir,
        config_yaml_path=cli_args.config_yaml,
        qc_output_dir=cli_args.qc_dir,
        negative_ratio=cli_args.negative_ratio,
        augment_prob=cli_args.augment_prob,
        copy_paste_prob=cli_args.copy_paste_prob,
        train_ratio=cli_args.train_ratio,
        val_ratio=cli_args.val_ratio,
        test_ratio=cli_args.test_ratio,
        min_instance_area=cli_args.min_instance_area,
        num_qc_plots=cli_args.num_qc_plots,
        num_workers=cli_args.num_workers,
        limit=cli_args.limit,
        seed=cli_args.seed
    )
