
import os
import sys
import shutil
import cv2
import math
import numpy as np
import random
import yaml
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


