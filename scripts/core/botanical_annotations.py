
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


