#!/usr/bin/env python3
"""
===============================================================================
Script: 02_segment_and_extract.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Production Runner for Step 02:
    Direct voucher image segmentation and leaf extraction pipeline for
    downstream Elliptic Fourier Analysis (EFA).
    
    Key Capabilities:
      1. Direct Ingestion: Reads curated_vouchers.csv directly without intermediate
         staging or symlink copying.
      2. PointRend Segmentation: Direct inference using the fine-tuned checkpoint
         (models/lm2_packera_pcd_finetuned.pth).
      3. Decoupled Asynchronous Scale Detection: Hough-transform ruler detection
         runs in parallel; scale failure never aborts silhouette extraction.
      4. In-Sheet Centroid Clustering: Groups organs into plant individuals
         without separate DBSCAN micro-scripts.
      5. 2-Path Botanical Leaf Extraction:
           - Tier 1 (Direct Pristine): UCS >= 0.85 & Solidity >= 0.72.
           - Tier 2 (Hemi-Blade Bilateral Reflection): Midrib axis alignment,
             half-blade defect gating, and symmetrical reflection.
           - Failed QC Routing: Heavily clumped rosettes or unsegmentable leaves
             routed to data/tables/failed_qc_vouchers.csv.
      6. Standardized Contour Export: Saves continuous 2D (x, y) coordinates
         to data/contours/{catalogNumber}_leaf{id}.csv and binary silhouettes
         to data/masks/.

Usage:
    python scripts/pipeline/02_segment_and_extract.py \\
        --vouchers data/tables/curated_vouchers.csv \\
        --model-weights models/lm2_packera_pcd_finetuned.pth \\
        --output-dir data/ \\
        --device cuda
===============================================================================
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# LeafMachine2 and Detectron2 paths
LM2_DIR = PROJECT_ROOT / "LeafMachine2" / "leafmachine2"
if LM2_DIR.exists() and str(LM2_DIR) not in sys.path:
    sys.path.insert(0, str(LM2_DIR))

DETECTRON2_PATH = LM2_DIR / "segmentation" / "detectron2"
if DETECTRON2_PATH.exists() and str(DETECTRON2_PATH) not in sys.path:
    sys.path.insert(0, str(DETECTRON2_PATH))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("SegmentAndExtract")


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class DetectedInstance:
    """Represents a segmented botanical instance on a herbarium sheet."""
    catalog_number: str
    leaf_id: int
    bbox: Tuple[int, int, int, int]  # (ymin, xmin, ymax, xmax)
    mask: np.ndarray                 # Binary uint8 mask (H, W)
    score: float
    class_id: int
    plant_individual_id: int = 0
    assigned_tier: str = "unassigned"
    ucs_score: float = 0.0
    solidity: float = 0.0
    midrib_angle_deg: float = 0.0
    pixels_per_mm: Optional[float] = None
    mask_path: Optional[str] = None
    contour_path: Optional[str] = None
    rejection_reason: Optional[str] = None


# =============================================================================
# Decoupled Asynchronous Ruler Detector (Hough Transform)
# =============================================================================

def detect_ruler_scale_hough(image: np.ndarray) -> Optional[float]:
    """
    Asynchronous Hough-transform ruler detection.
    Scans margin regions of the herbarium sheet for periodic tick marks.
    Returns pixels_per_mm or None if not confidently detected.
    """
    if image is None or image.size == 0:
        return None

    h, w = image.shape[:2]
    # Rulers are typically along sheet borders (margins)
    margins = [
        image[int(h * 0.75):, :],               # Bottom
        image[:int(h * 0.25), :],               # Top
        image[:, :int(w * 0.20)],               # Left
        image[:, int(w * 0.80):],               # Right
    ]

    candidate_scales: List[float] = []

    for region in margins:
        if region.size == 0:
            continue
        gray = cv2.cvtColor(region, cv2.COLOR_BGR2GRAY) if len(region.shape) == 3 else region
        
        # Downscale if extremely large for rapid Hough transform
        rh, rw = gray.shape[:2]
        max_dim = max(rh, rw)
        scale_factor = 1.0
        if max_dim > 1500:
            scale_factor = 1500.0 / max_dim
            gray = cv2.resize(gray, (int(rw * scale_factor), int(rh * scale_factor)), interpolation=cv2.INTER_AREA)

        # Bilateral / Gaussian blur + Canny edge detection
        blurred = cv2.GaussianBlur(gray, (3, 3), 0)
        edges = cv2.Canny(blurred, 30, 100)

        # Detect line segments
        min_line_len = max(int(12 * scale_factor), 10)
        lines = cv2.HoughLinesP(edges, rho=1, theta=np.pi / 180, threshold=25,
                                minLineLength=min_line_len, maxLineGap=4)
        if lines is None or len(lines) < 8:
            continue

        # Extract line angles and lengths
        line_data = []
        for line in lines:
            x1, y1, x2, y2 = line[0]
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy)
            angle = math.degrees(math.atan2(dy, dx)) % 180.0
            line_data.append((x1, y1, x2, y2, angle, length))

        # Check for dominant parallel angle (horizontal or vertical ticks)
        angles = np.array([ld[4] for ld in line_data])
        hist, bin_edges = np.histogram(angles, bins=18, range=(0, 180))
        dominant_bin = np.argmax(hist)
        if hist[dominant_bin] < 6:
            continue

        dom_angle_low = bin_edges[dominant_bin]
        dom_angle_high = bin_edges[dominant_bin + 1]

        # Filter lines matching dominant tick orientation
        matching_lines = [ld for ld in line_data if dom_angle_low <= ld[4] <= dom_angle_high]
        if len(matching_lines) < 6:
            continue

        # Project lines onto orthogonal axis to determine inter-tick pitch
        dom_angle = (dom_angle_low + dom_angle_high) / 2.0
        ortho_rad = math.radians(dom_angle + 90.0)
        proj_coords = []
        for ld in matching_lines:
            mid_x, mid_y = (ld[0] + ld[2]) / 2.0, (ld[1] + ld[3]) / 2.0
            p = mid_x * math.cos(ortho_rad) + mid_y * math.sin(ortho_rad)
            proj_coords.append(p)

        proj_coords = np.sort(proj_coords)
        tick_positions = []
        for p in proj_coords:
            if not tick_positions or abs(p - tick_positions[-1]) > 5.0:
                tick_positions.append(p)

        if len(tick_positions) >= 5:
            deltas = np.diff(tick_positions)
            med_delta = float(np.median(deltas))
            mad = float(np.median(np.abs(deltas - med_delta)))
            # Consistent periodic tick spacing condition
            if mad < 0.25 * med_delta and (8.0 <= med_delta / scale_factor <= 150.0):
                px_per_mm = med_delta / scale_factor
                candidate_scales.append(px_per_mm)

    if candidate_scales:
        return float(np.median(candidate_scales))
    return None


# =============================================================================
# Geometric Quality Gatekeeper & Pose Estimation
# =============================================================================

def compute_geometric_metrics(
    mask: np.ndarray
) -> Tuple[float, float, float, Tuple[int, int], Tuple[int, int]]:
    """
    Computes:
      - ucs: Unoccluded Completeness Score = Area_mask / (0.70 * Area_min_rect)
      - solidity: Area_mask / Area_convex_hull
      - midrib_angle_deg: Angle of longitudinal midrib axis from apex to base
      - p_apex: Estimated apex coordinate (x, y)
      - p_base: Estimated base/petiole coordinate (x, y)
    """
    area_mask = float(np.count_nonzero(mask))
    if area_mask < 30:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return 0.0, 0.0, 0.0, (0, 0), (0, 0)

    main_contour = max(contours, key=cv2.contourArea)

    # 1. Convex Hull & Solidity
    hull = cv2.convexHull(main_contour)
    area_hull = float(cv2.contourArea(hull))
    solidity = (area_mask / area_hull) if area_hull > 0 else 0.0
    solidity = min(max(solidity, 0.0), 1.0)

    # 2. Expected Area via Minimum Area Rect (Oriented Bounding Box)
    rect = cv2.minAreaRect(main_contour)
    rect_w, rect_h = rect[1]
    area_rect = max(rect_w * rect_h, 1.0)
    expected_area = 0.70 * area_rect
    ucs = min(max(area_mask / expected_area, 0.0), 1.0)

    # 3. Longitudinal Midrib Axis via cv2.fitLine
    line_params = cv2.fitLine(main_contour, cv2.DIST_L2, 0, 0.01, 0.01)
    vx, vy, x0, y0 = [float(v[0]) for v in line_params]
    angle_rad = math.atan2(vy, vx)
    angle_deg = math.degrees(angle_rad)

    # Project contour points along longitudinal midrib axis
    pts = main_contour.reshape(-1, 2).astype(np.float32)
    projections = (pts[:, 0] - x0) * vx + (pts[:, 1] - y0) * vy

    min_idx = int(np.argmin(projections))
    max_idx = int(np.argmax(projections))

    p_apex = (int(pts[min_idx, 0]), int(pts[min_idx, 1]))
    p_base = (int(pts[max_idx, 0]), int(pts[max_idx, 1]))

    return ucs, solidity, angle_deg, p_apex, p_base


# =============================================================================
# 2-Path Extraction: Tier 1 (Pristine) & Tier 2 (Bilateral Reflection)
# =============================================================================

def extract_tier1_pristine(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    output_dir: Path
) -> str:
    """
    Tier 1: Direct Pristine Silhouette extraction.
    Saves binary mask directly to data/masks/.
    """
    mask_dir = output_dir / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    tier1_dir = mask_dir / "tier1_intact"
    tier1_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{catalog_number}_leaf{leaf_id}.png"
    main_save_path = mask_dir / filename
    tier_save_path = tier1_dir / filename

    bin_mask = (mask > 0).astype(np.uint8) * 255
    cv2.imwrite(str(main_save_path), bin_mask)
    cv2.imwrite(str(tier_save_path), bin_mask)
    return str(main_save_path)


def extract_tier2_reflected(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    p_apex: Tuple[int, int],
    p_base: Tuple[int, int],
    output_dir: Path
) -> Tuple[Optional[str], Optional[np.ndarray]]:
    """
    Tier 2: Hemi-Blade Bilateral Symmetry Reflection.
    Aligns mask along midrib axis, assesses margin defect profile of upper vs
    lower half, reflects the clean half across the midrib axis symmetrically,
    and applies morphological closing to produce a watertight closed silhouette.
    """
    h, w = mask.shape[:2]
    if h < 15 or w < 15:
        return None, None

    dx = p_base[0] - p_apex[0]
    dy = p_base[1] - p_apex[1]
    angle = math.degrees(math.atan2(dy, dx))

    cx, cy = w / 2.0, h / 2.0
    rot_mat = cv2.getRotationMatrix2D((cx, cy), -angle, 1.0)
    aligned = cv2.warpAffine(mask, rot_mat, (w, h), flags=cv2.INTER_NEAREST)

    y_indices, _ = np.where(aligned > 0)
    if len(y_indices) < 30:
        return None, None

    y_midrib = int(np.median(y_indices))

    upper_half = np.zeros_like(aligned)
    lower_half = np.zeros_like(aligned)
    upper_half[:y_midrib, :] = aligned[:y_midrib, :]
    lower_half[y_midrib:, :] = aligned[y_midrib:, :]

    # Vectorized margin profile detection across columns
    upper_has_fg = upper_half > 0
    lower_has_fg = lower_half > 0

    col_has_upper = np.any(upper_has_fg, axis=0)
    col_has_lower = np.any(lower_has_fg, axis=0)

    # Profile arrays
    upper_profile = np.full(w, y_midrib, dtype=np.int32)
    lower_profile = np.full(w, y_midrib, dtype=np.int32)

    for x in range(w):
        if col_has_upper[x]:
            upper_profile[x] = np.min(np.where(upper_has_fg[:, x])[0])
        if col_has_lower[x]:
            lower_profile[x] = np.max(np.where(lower_has_fg[:, x])[0])

    # Assess margin defect score (jumps/discontinuities)
    upper_diff = np.abs(np.diff(upper_profile))
    lower_diff = np.abs(np.diff(lower_profile))
    upper_defects = float(np.sum(upper_diff > 8) * 10.0 + np.std(upper_diff))
    lower_defects = float(np.sum(lower_diff > 8) * 10.0 + np.std(lower_diff))

    threshold_defect = 35.0
    if upper_defects < threshold_defect and upper_defects <= lower_defects:
        selected_half = "upper"
    elif lower_defects < threshold_defect:
        selected_half = "lower"
    else:
        return None, None

    # Vectorized reflection across midrib
    reflected = np.zeros_like(aligned)
    if selected_half == "upper":
        reflected[:y_midrib, :] = upper_half[:y_midrib, :]
        for y in range(y_midrib):
            target_y = y_midrib + (y_midrib - y)
            if 0 <= target_y < h:
                reflected[target_y, :] = upper_half[y, :]
    else:
        reflected[y_midrib:, :] = lower_half[y_midrib:, :]
        for y in range(y_midrib, h):
            target_y = y_midrib - (y - y_midrib)
            if 0 <= target_y < h:
                reflected[target_y, :] = lower_half[y, :]

    # Morphological closing to seal midrib seam
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    reflected = cv2.morphologyEx(reflected, cv2.MORPH_CLOSE, kernel)
    reflected_bin = (reflected > 127).astype(np.uint8) * 255

    # Verify reflected silhouette quality
    ref_ucs, ref_solidity, _, _, _ = compute_geometric_metrics(reflected_bin)
    if ref_solidity < 0.68 or ref_ucs < 0.50:
        return None, None

    mask_dir = output_dir / "masks"
    mask_dir.mkdir(parents=True, exist_ok=True)
    tier2_dir = mask_dir / "tier2_reflected"
    tier2_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{catalog_number}_leaf{leaf_id}.png"
    main_save_path = mask_dir / filename
    tier_save_path = tier2_dir / filename

    cv2.imwrite(str(main_save_path), reflected_bin)
    cv2.imwrite(str(tier_save_path), reflected_bin)
    return str(main_save_path), reflected_bin


# =============================================================================
# Standardized Coordinate Contour Exporter
# =============================================================================

def export_standardized_contour(
    mask: np.ndarray,
    catalog_number: str,
    leaf_id: int,
    output_dir: Path,
    num_points: int = 120
) -> Optional[str]:
    """
    Extracts vectorized 2D (x, y) boundary coordinates from a binary silhouette mask,
    normalizes coordinates (x_norm, y_norm) for downstream Elliptic Fourier Analysis,
    and saves to data/contours/{catalogNumber}_leaf{id}.csv.
    """
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        return None

    largest_cnt = max(contours, key=cv2.contourArea).reshape(-1, 2)
    if len(largest_cnt) < 10:
        return None

    # Resample smoothly along perimeter to standardized point count
    pts = largest_cnt.astype(np.float64)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])

    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    valid = dt > 1e-7
    if np.sum(valid) < 5:
        return None

    pts = pts[np.concatenate([[True], valid])]
    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    cum_dist = np.concatenate([[0.0], np.cumsum(dt)])
    total_len = cum_dist[-1]
    if total_len <= 0:
        return None

    target_distances = np.linspace(0.0, total_len, num_points, endpoint=False)
    x_resamp = np.interp(target_distances, cum_dist, pts[:, 0])
    y_resamp = np.interp(target_distances, cum_dist, pts[:, 1])

    min_x, max_x = np.min(x_resamp), np.max(x_resamp)
    min_y, max_y = np.min(y_resamp), np.max(y_resamp)
    span_x = max(max_x - min_x, 1e-6)
    span_y = max(max_y - min_y, 1e-6)

    contour_df = pd.DataFrame({
        "catalogNumber": catalog_number,
        "leaf_id": leaf_id,
        "point_index": np.arange(num_points),
        "x": np.round(x_resamp, 2),
        "y": np.round(y_resamp, 2),
        "x_norm": np.round((x_resamp - min_x) / span_x, 6),
        "y_norm": np.round((y_resamp - min_y) / span_y, 6),
    })

    contour_dir = output_dir / "contours"
    contour_dir.mkdir(parents=True, exist_ok=True)
    out_csv = contour_dir / f"{catalog_number}_leaf{leaf_id}.csv"
    contour_df.to_csv(out_csv, index=False)
    return str(out_csv)


# =============================================================================
# In-Sheet Centroid Clustering for Multi-Plant Specimen Sheets
# =============================================================================

def cluster_plant_individuals(
    instances: List[DetectedInstance],
    sheet_width: int,
    sheet_height: int,
    eps_ratio: float = 0.15
) -> List[DetectedInstance]:
    """
    Vectorized centroid clustering across detected organs on a specimen sheet.
    Assigns plant_individual_id to avoid trait averaging across distinct plants.
    """
    if not instances:
        return instances

    if len(instances) == 1:
        instances[0].plant_individual_id = 0
        return instances

    centroids = []
    for inst in instances:
        ymin, xmin, ymax, xmax = inst.bbox
        cx = (xmin + xmax) / 2.0
        cy = (ymin + ymax) / 2.0
        centroids.append([cx, cy])

    coords = np.array(centroids, dtype=np.float32)
    eps = eps_ratio * float(max(sheet_width, sheet_height))

    clustering = DBSCAN(eps=eps, min_samples=1).fit(coords)
    labels = clustering.labels_

    unique_labels = sorted(list(set(labels)))
    label_map = {lbl: idx for idx, lbl in enumerate(unique_labels)}

    for inst, lbl in zip(instances, labels):
        inst.plant_individual_id = label_map.get(lbl, 0)

    return instances


# =============================================================================
# Direct PointRend Model Inference Engine
# =============================================================================

class PointRendInferenceEngine:
    """Encapsulates Detectron2 PointRend model loading and inference."""

    def __init__(
        self,
        weights_path: Union[str, Path],
        device: str = "cuda",
        score_thresh: float = 0.50
    ) -> None:
        self.weights_path = Path(weights_path)
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.score_thresh = score_thresh
        self.predictor = None
        self.model = None
        self.cfg = None
        self.aug = None
        self._load_model()

    def _load_model(self) -> None:
        """Loads Detectron2 PointRend model from weights checkpoint."""
        try:
            from detectron2.config import get_cfg
            from detectron2.projects import point_rend
            from detectron2.checkpoint import DetectionCheckpointer
            from detectron2.modeling import build_model
            import detectron2.data.transforms as T
        except ImportError as e:
            logger.warning(f"Detectron2 not available: {e}. Running in simulation/fallback mode.")
            return

        cfg = get_cfg()
        point_rend.add_pointrend_config(cfg)

        # Locate accompanying configuration YAML
        cfg_candidates = [
            self.weights_path.parent / "cfg_output.yaml",
            PROJECT_ROOT / "LeafMachine2" / "leafmachine2" / "segmentation" / "models" / "Packera_LeafPriority" / "cfg_output.yaml",
            PROJECT_ROOT / "models" / "cfg_output.yaml",
        ]
        cfg_yaml = None
        for cand in cfg_candidates:
            if cand.exists():
                cfg_yaml = cand
                break

        if cfg_yaml:
            cfg.merge_from_file(str(cfg_yaml))
        else:
            logger.warning("No cfg_output.yaml found; using default PointRend config.")

        cfg.MODEL.WEIGHTS = str(self.weights_path.resolve())
        cfg.MODEL.ROI_HEADS.SCORE_THRESH_TEST = self.score_thresh
        cfg.MODEL.DEVICE = self.device

        model = build_model(cfg)
        model.eval()
        checkpointer = DetectionCheckpointer(model)
        checkpointer.load(cfg.MODEL.WEIGHTS)

        self.cfg = cfg
        self.model = model
        self.aug = T.ResizeShortestEdge([cfg.INPUT.MIN_SIZE_TEST, cfg.INPUT.MIN_SIZE_TEST], cfg.INPUT.MAX_SIZE_TEST)
        logger.info(f"Loaded PointRend model from {self.weights_path} on {self.device}")

    def predict(self, image_bgr: np.ndarray) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """
        Runs inference on BGR image.
        Returns list of (bbox, mask, score, class_id).
        """
        if self.model is None:
            return []

        import torch
        height, width = image_bgr.shape[:2]
        with torch.no_grad():
            img_transformed = self.aug.get_transform(image_bgr).apply_image(image_bgr)
            tensor_img = torch.as_tensor(img_transformed.astype("float32").transpose(2, 0, 1))
            inputs = [{"image": tensor_img, "height": height, "width": width}]
            preds = self.model(inputs)[0]

        instances = preds["instances"].to("cpu")
        results = []

        if not instances.has("pred_masks"):
            return results

        scores = instances.scores.numpy()
        classes = instances.pred_classes.numpy()
        boxes = instances.pred_boxes.tensor.numpy()
        masks = instances.pred_masks.numpy()

        for i in range(len(scores)):
            if scores[i] < self.score_thresh:
                continue
            box = boxes[i]
            ymin, xmin, ymax, xmax = int(box[1]), int(box[0]), int(box[3]), int(box[2])
            mask_uint8 = (masks[i] > 0).astype(np.uint8) * 255
            results.append(((ymin, xmin, ymax, xmax), mask_uint8, float(scores[i]), int(classes[i])))

        return results


def torch_cuda_available() -> bool:
    """Helper to safely check CUDA availability."""
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


# =============================================================================
# Core Pipeline Execution Runner
# =============================================================================

class SegmentAndExtractPipeline:
    """Unified runner for leaf segmentation, 2-path extraction, and manifest logging."""

    def __init__(
        self,
        vouchers_csv: Union[str, Path],
        model_weights: Union[str, Path],
        output_dir: Union[str, Path],
        device: str = "cuda",
        min_solidity: float = 0.72,
        min_ucs: float = 0.85,
        score_thresh: float = 0.40,
        num_ruler_workers: int = 4
    ) -> None:
        self.vouchers_csv = Path(vouchers_csv)
        self.model_weights = Path(model_weights)
        self.output_dir = Path(output_dir)
        self.device = device
        self.min_solidity = min_solidity
        self.min_ucs = min_ucs
        self.score_thresh = score_thresh
        self.num_ruler_workers = num_ruler_workers

        self.masks_dir = self.output_dir / "masks"
        self.contours_dir = self.output_dir / "contours"
        self.tables_dir = self.output_dir / "tables"

        for d in [self.masks_dir, self.contours_dir, self.tables_dir]:
            d.mkdir(parents=True, exist_ok=True)

        self.ruler_executor = concurrent.futures.ThreadPoolExecutor(max_workers=self.num_ruler_workers)
        self.model_engine = PointRendInferenceEngine(
            weights_path=self.model_weights,
            device=self.device,
            score_thresh=self.score_thresh
        )

    def shutdown(self) -> None:
        """Cleans up background thread pools."""
        self.ruler_executor.shutdown(wait=False)

    def process_voucher(
        self,
        catalog_number: str,
        image_path: Path
    ) -> Tuple[List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Processes a single herbarium voucher specimen sheet."""
        if not image_path.exists():
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "image_not_found",
                "details": f"File does not exist at {image_path}"
            }
            return [], failure_record

        img_bgr = cv2.imread(str(image_path))
        if img_bgr is None:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "corrupted_image",
                "details": "Failed to decode image with OpenCV"
            }
            return [], failure_record

        sheet_h, sheet_w = img_bgr.shape[:2]

        # 1. Asynchronously launch decoupled Hough ruler detection
        ruler_future = self.ruler_executor.submit(detect_ruler_scale_hough, img_bgr)

        # 2. Run PointRend instance segmentation
        raw_detections = self.model_engine.predict(img_bgr)

        # Collect asynchronous scale result (with graceful fallback)
        try:
            pixels_per_mm = ruler_future.result(timeout=5.0)
        except Exception as e:
            logger.debug(f"Ruler detection timeout or exception for {catalog_number}: {e}")
            pixels_per_mm = None

        if not raw_detections:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "no_leaves_detected",
                "details": "Model produced 0 candidate leaf segmentations"
            }
            return [], failure_record

        # Build detected instances
        instances: List[DetectedInstance] = []
        for idx, (bbox, mask, score, class_id) in enumerate(raw_detections, 1):
            instances.append(DetectedInstance(
                catalog_number=catalog_number,
                leaf_id=idx,
                bbox=bbox,
                mask=mask,
                score=score,
                class_id=class_id,
                pixels_per_mm=pixels_per_mm
            ))

        # 3. Vectorized Centroid Spatial Clustering
        instances = cluster_plant_individuals(instances, sheet_w, sheet_h)

        extracted_records: List[Dict[str, Any]] = []
        rejected_instances: List[DetectedInstance] = []

        # 4. 2-Path Botanical Leaf Extraction Logic
        for inst in instances:
            ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(inst.mask)
            inst.ucs_score = ucs
            inst.solidity = solidity
            inst.midrib_angle_deg = angle_deg

            # Path A: Tier 1 Direct Pristine
            if ucs >= self.min_ucs and solidity >= self.min_solidity:
                inst.assigned_tier = "tier1"
                mask_path = extract_tier1_pristine(inst.mask, catalog_number, inst.leaf_id, self.output_dir)
                contour_path = export_standardized_contour(inst.mask, catalog_number, inst.leaf_id, self.output_dir)
                inst.mask_path = mask_path
                inst.contour_path = contour_path
            else:
                # Path B: Tier 2 Hemi-Blade Bilateral Reflection
                mask_path, reflected_mask = extract_tier2_reflected(
                    inst.mask, catalog_number, inst.leaf_id, p_apex, p_base, self.output_dir
                )
                if mask_path and reflected_mask is not None:
                    inst.assigned_tier = "tier2"
                    contour_path = export_standardized_contour(
                        reflected_mask, catalog_number, inst.leaf_id, self.output_dir
                    )
                    inst.mask_path = mask_path
                    inst.contour_path = contour_path
                else:
                    inst.assigned_tier = "rejected"
                    inst.rejection_reason = "clumped_rosette_or_severe_occlusion"
                    rejected_instances.append(inst)
                    continue

            # Record successfully extracted leaf
            extracted_records.append({
                "catalogNumber": catalog_number,
                "plant_individual_id": inst.plant_individual_id,
                "leaf_id": inst.leaf_id,
                "assigned_tier": inst.assigned_tier,
                "ucs_score": round(inst.ucs_score, 4),
                "solidity": round(inst.solidity, 4),
                "midrib_angle_deg": round(inst.midrib_angle_deg, 2),
                "pixels_per_mm": round(inst.pixels_per_mm, 4) if inst.pixels_per_mm else np.nan,
                "mask_path": inst.mask_path,
                "contour_path": inst.contour_path
            })

        # Route vouchers with 0 valid silhouettes to failed manifest
        failure_record = None
        if not extracted_records:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "clumped_rosette_occlusion",
                "details": f"{len(rejected_instances)} leaves detected but all rejected by geometric gatekeeping"
            }

        return extracted_records, failure_record

    def run(self, limit: Optional[int] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Executes full extraction pipeline across curated vouchers."""
        logger.info("==================================================================")
        logger.info("Starting Step 02: Direct Segmentation & Botanical Leaf Extraction")
        logger.info(f"Vouchers Metadata: {self.vouchers_csv}")
        logger.info(f"PointRend Model: {self.model_weights}")
        logger.info(f"Output Directory: {self.output_dir}")
        logger.info(f"Device: {self.device}")
        logger.info(f"Gatekeeper Thresholds: Min Solidity={self.min_solidity}, Min UCS={self.min_ucs}")
        logger.info("==================================================================")

        if not self.vouchers_csv.exists():
            raise FileNotFoundError(f"Curated vouchers table not found at {self.vouchers_csv}")

        df_vouchers = pd.read_csv(self.vouchers_csv)
        logger.info(f"Loaded {len(df_vouchers)} voucher records.")

        if "catalogNumber" not in df_vouchers.columns or "image_path" not in df_vouchers.columns:
            raise KeyError("Vouchers CSV must contain 'catalogNumber' and 'image_path' columns.")

        if limit and limit > 0:
            df_vouchers = df_vouchers.head(limit)
            logger.info(f"Applied execution limit: processing first {limit} vouchers.")

        all_extracted_records: List[Dict[str, Any]] = []
        all_failed_records: List[Dict[str, Any]] = []

        total_vouchers = len(df_vouchers)
        for i, (_, row) in enumerate(df_vouchers.iterrows(), 1):
            cat_num = str(row["catalogNumber"]).strip()
            raw_img_path = Path(str(row["image_path"]).strip())
            if not raw_img_path.is_absolute():
                raw_img_path = PROJECT_ROOT / raw_img_path

            extracted, failure = self.process_voucher(cat_num, raw_img_path)
            all_extracted_records.extend(extracted)
            if failure:
                all_failed_records.append(failure)

            if i % 25 == 0 or i == total_vouchers:
                logger.info(f"Progress: {i}/{total_vouchers} vouchers processed "
                            f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC).")

        # Save manifests
        extracted_df = pd.DataFrame(all_extracted_records)
        if extracted_df.empty:
            extracted_df = pd.DataFrame(columns=[
                "catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier",
                "ucs_score", "solidity", "midrib_angle_deg", "pixels_per_mm",
                "mask_path", "contour_path"
            ])

        failed_df = pd.DataFrame(all_failed_records)
        if failed_df.empty:
            failed_df = pd.DataFrame(columns=[
                "catalogNumber", "image_path", "failure_reason", "details"
            ])

        manifest_path = self.tables_dir / "extracted_leaves_manifest.csv"
        failed_qc_path = self.tables_dir / "failed_qc_vouchers.csv"

        extracted_df.to_csv(manifest_path, index=False)
        failed_df.to_csv(failed_qc_path, index=False)

        logger.info("==================================================================")
        logger.info(f"Extraction completed!")
        logger.info(f"  Extracted Leaves: {len(extracted_df)} -> {manifest_path}")
        logger.info(f"  Failed QC Vouchers: {len(failed_df)} -> {failed_qc_path}")
        if not extracted_df.empty and "assigned_tier" in extracted_df.columns:
            tier_dist = extracted_df["assigned_tier"].value_counts().to_dict()
            for tier, count in tier_dist.items():
                pct = (count / len(extracted_df)) * 100
                logger.info(f"    {tier.upper()}: {count} ({pct:.1f}%)")
        logger.info("==================================================================")

        self.shutdown()
        return extracted_df, failed_df


# =============================================================================
# CLI Interface
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Direct PointRend segmentation, 2-path leaf extraction, and contour coordinate export."
    )
    parser.add_argument(
        "--vouchers",
        type=Path,
        default=PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv",
        help="Path to curated vouchers CSV containing catalogNumber and image_path."
    )
    parser.add_argument(
        "--model-weights",
        type=Path,
        default=PROJECT_ROOT / "models" / "lm2_packera_pcd_finetuned.pth",
        help="Path to fine-tuned PointRend model checkpoint."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_ROOT / "data",
        help="Root output directory for masks, contours, and manifest tables."
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cuda" if torch_cuda_available() else "cpu",
        choices=["cuda", "cpu"],
        help="Compute device for inference (default: cuda if available, else cpu)."
    )
    parser.add_argument(
        "--min-solidity",
        type=float,
        default=0.72,
        help="Minimum solidity for Tier 1 direct pristine leaf extraction (default: 0.72)."
    )
    parser.add_argument(
        "--min-ucs",
        type=float,
        default=0.85,
        help="Minimum Unoccluded Completeness Score for Tier 1 (default: 0.85)."
    )
    parser.add_argument(
        "--score-thresh",
        type=float,
        default=0.40,
        help="Detection confidence score threshold (default: 0.40)."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of vouchers to process."
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pipeline = SegmentAndExtractPipeline(
        vouchers_csv=args.vouchers,
        model_weights=args.model_weights,
        output_dir=args.output_dir,
        device=args.device,
        min_solidity=args.min_solidity,
        min_ucs=args.min_ucs,
        score_thresh=args.score_thresh
    )
    pipeline.run(limit=args.limit)


if __name__ == "__main__":
    main()
