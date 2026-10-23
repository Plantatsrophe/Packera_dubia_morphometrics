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
import tempfile
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

from scripts.core.config import PipelineConfig
from scripts.vision.configure_leafmachine2 import (
    generate_high_performance_config,
    resolve_pcd_weights,
)
from scripts.vision.lm2_data_loader import (
    discover_lm2_candidates,
    LeafCandidate,
    LM2_LEAF_CLASS_IDS,
    LM2_PCD_CLASSES,
    LM2_REPRODUCTIVE_CLASS_IDS,
    load_precomputed_pcd_labels,
    load_ruler_calibrations,
    parse_bbox_from_filename,
)
from scripts.vision.lm2_geometry_utils import (
    classify_and_route_leaf,
    cluster_plant_individuals,
    compute_geometric_metrics,
    crop_dense_rosette_tier4,
    detect_folded_leaf,
    detect_ruler_scale_hough,
    export_standardized_contour,
    extract_tier1_pristine,
    extract_tier2_reflected,
    generate_sam2_prompts_from_lm2_mask,
    homologize_contour_starting_point,
    is_botanical_dissection,
    is_topologically_connected_mask,
    LeafRoutingResult,
    normalize_assigned_tier,
    SAM2PromptBundle,
    TIER_1_DIRECT,
    TIER_2_REFLECTED,
    trim_petiole_tail,
)
from scripts.vision.capitulum_phenology_classifier import (
    classify_single_capitulum,
    classify_voucher_phenology,
)


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
# Data Structures & Two-Stage Crop Reprojection Utilities
# =============================================================================

@dataclass
class DetectedInstance:
    """Represents a segmented botanical instance on a herbarium sheet."""
    catalog_number: str
    leaf_id: int
    bbox: Tuple[int, int, int, int]  # (ymin, xmin, ymax, xmax) in full-sheet coordinates
    mask: np.ndarray                 # Binary uint8 mask (H, W)
    score: float
    class_id: int
    plant_individual_id: int = 0
    assigned_tier: str = "unassigned"
    ucs_score: float = 0.0
    solidity: float = 0.0
    midrib_angle_deg: float = 0.0
    pixels_per_mm: Optional[float] = None
    is_calibrated_mm: bool = False
    crop_offset: Tuple[int, int] = (0, 0)  # (x0, y0) crop offset on full sheet
    mask_path: Optional[str] = None
    contour_path: Optional[str] = None
    rejection_reason: Optional[str] = None


def crop_leaf_bbox_with_padding(
    image_bgr: np.ndarray,
    bbox: Tuple[int, int, int, int],
    pad_ratio: float = 0.12,
) -> Tuple[np.ndarray, Tuple[int, int, int, int]]:
    """
    Crops a candidate leaf bounding box (ymin, xmin, ymax, xmax) from the native-DPI
    voucher sheet with isotropic padding (default 12%), clamped to sheet boundaries.
    Returns (crop_bgr, (crop_y0, crop_x0, crop_y1, crop_x1)).
    """
    sheet_h, sheet_w = image_bgr.shape[:2]
    ymin, xmin, ymax, xmax = [int(v) for v in bbox]
    box_h = max(0, ymax - ymin)
    box_w = max(0, xmax - xmin)
    pad_h = int(round(box_h * pad_ratio))
    pad_w = int(round(box_w * pad_ratio))

    crop_y0 = max(0, ymin - pad_h)
    crop_x0 = max(0, xmin - pad_w)
    crop_y1 = min(sheet_h, ymax + pad_h)
    crop_x1 = min(sheet_w, xmax + pad_w)

    if crop_y1 <= crop_y0 or crop_x1 <= crop_x0:
        return np.zeros((0, 0, 3), dtype=image_bgr.dtype), (crop_y0, crop_x0, crop_y0, crop_x0)

    return image_bgr[crop_y0:crop_y1, crop_x0:crop_x1].copy(), (crop_y0, crop_x0, crop_y1, crop_x1)


def reproject_crop_mask_to_sheet(
    crop_mask: np.ndarray,
    crop_box: Tuple[int, int, int, int],
    sheet_shape: Tuple[int, int],
) -> Tuple[np.ndarray, Tuple[int, int, int, int]]:
    """
    Projects a binary instance mask (crop_h, crop_w) from a padded Stage 1 crop
    back into full-sheet coordinates (sheet_h, sheet_w) and returns
    (sheet_mask, full_sheet_bbox).
    """
    sheet_h, sheet_w = int(sheet_shape[0]), int(sheet_shape[1])
    crop_y0, crop_x0, crop_y1, crop_x1 = [int(v) for v in crop_box]
    crop_y0 = max(0, min(sheet_h, crop_y0))
    crop_x0 = max(0, min(sheet_w, crop_x0))
    crop_y1 = max(crop_y0, min(sheet_h, crop_y1))
    crop_x1 = max(crop_x0, min(sheet_w, crop_x1))

    target_h = crop_y1 - crop_y0
    target_w = crop_x1 - crop_x0
    sheet_mask = np.zeros((sheet_h, sheet_w), dtype=np.uint8)
    if target_h <= 0 or target_w <= 0 or crop_mask is None or crop_mask.size == 0:
        return sheet_mask, (crop_y0, crop_x0, crop_y1, crop_x1)

    mask_bin = (crop_mask > 0).astype(np.uint8) * 255
    if mask_bin.shape[0] != target_h or mask_bin.shape[1] != target_w:
        mask_bin = cv2.resize(mask_bin, (target_w, target_h), interpolation=cv2.INTER_NEAREST)

    sheet_mask[crop_y0:crop_y1, crop_x0:crop_x1] = mask_bin

    ys, xs = np.where(mask_bin > 0)
    if len(ys) > 0 and len(xs) > 0:
        ymin_s = crop_y0 + int(np.min(ys))
        ymax_s = crop_y0 + int(np.max(ys)) + 1
        xmin_s = crop_x0 + int(np.min(xs))
        xmax_s = crop_x0 + int(np.max(xs)) + 1
        full_bbox = (ymin_s, xmin_s, min(sheet_h, ymax_s), min(sheet_w, xmax_s))
    else:
        full_bbox = (crop_y0, crop_x0, crop_y1, crop_x1)

    return sheet_mask, full_bbox


# =============================================================================
# Capitulum Macro-Reproductive Metric Extraction
# =============================================================================

def compute_capitulum_metrics(
    mask: Optional[np.ndarray],
    bbox: Tuple[int, int, int, int],
    pixels_per_mm: Optional[float] = None,
    min_ar: float = 0.7,
    max_ar: float = 2.0,
) -> Optional[Dict[str, Any]]:
    """
    Computes involucre height (H_inv), width (W_inv), and aspect ratio (AR_inv = H/W)
    for a detected reproductive / capitulum instance.

    Filters for valid cylindrical capitula where min_ar <= H/W <= max_ar (default 0.7 to 2.0).
    Dimensions are computed via oriented bounding box (cv2.minAreaRect), resolving the longitudinal
    axis (peduncle-to-apex height) versus transverse diameter (involucre width) based on dominant vertical
    alignment on the herbarium sheet.
    """
    h_inv_px = 0.0
    w_inv_px = 0.0

    contours = []
    if mask is not None and isinstance(mask, np.ndarray) and mask.size > 0:
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if contours and cv2.contourArea(max(contours, key=cv2.contourArea)) >= 10:
        main_contour = max(contours, key=cv2.contourArea)
        rect = cv2.minAreaRect(main_contour)
        box = cv2.boxPoints(rect)
        v1 = box[1] - box[0]
        v2 = box[2] - box[1]
        len1 = float(np.linalg.norm(v1))
        len2 = float(np.linalg.norm(v2))

        # Herbarium specimens are mounted with flowering stalks ascending vertically (y-axis).
        # Resolve the longitudinal axis (height) as the edge with greater vertical projection (|Δy|/len).
        norm1 = max(len1, 1e-6)
        norm2 = max(len2, 1e-6)
        vert1 = abs(float(v1[1])) / norm1
        vert2 = abs(float(v2[1])) / norm2

        if vert1 >= vert2:
            h_inv_px = len1
            w_inv_px = len2
        else:
            h_inv_px = len2
            w_inv_px = len1
    else:
        # Fallback to axis-aligned bounding box
        ymin, xmin, ymax, xmax = bbox
        h_inv_px = float(max(ymax - ymin, 1))
        w_inv_px = float(max(xmax - xmin, 1))

    if w_inv_px <= 1e-6:
        return None

    aspect_ratio = h_inv_px / w_inv_px
    if not (min_ar <= aspect_ratio <= max_ar):
        return None

    h_inv_mm = (h_inv_px / pixels_per_mm) if (pixels_per_mm and pixels_per_mm > 0) else np.nan
    w_inv_mm = (w_inv_px / pixels_per_mm) if (pixels_per_mm and pixels_per_mm > 0) else np.nan

    return {
        "involucre_height_px": round(h_inv_px, 2),
        "involucre_width_px": round(w_inv_px, 2),
        "involucre_height_mm": round(h_inv_mm, 3) if not np.isnan(h_inv_mm) else np.nan,
        "involucre_width_mm": round(w_inv_mm, 3) if not np.isnan(w_inv_mm) else np.nan,
        "capitulum_aspect_ratio": round(aspect_ratio, 4),
    }


# =============================================================================
# Stage 1 LeafMachine2 Component Detector (LeafPriority.pt, 11 Classes)
# =============================================================================

class LM2Stage1Detector:
    """
    Encapsulates LeafMachine2's Stage 1 YOLOv5 Plant Component Detector (LeafPriority.pt).
    11 Canonical Classes:
      0: leaf_whole, 1: leaf_partial, 2: leaflet,
      3: seed_fruit_one, 4: seed_fruit_many, 5: flower_one, 6: flower_many, 7: bud,
      8: specimen, 9: roots, 10: wood.
    """

    CLASSES: Dict[int, str] = LM2_PCD_CLASSES
    LEAF_CLASS_IDS: Tuple[int, ...] = LM2_LEAF_CLASS_IDS
    REPRODUCTIVE_CLASS_IDS: Tuple[int, ...] = LM2_REPRODUCTIVE_CLASS_IDS

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        conf_thres: float = 0.25,
        iou_thres: float = 0.45,
        img_size: int = 1280,
        lm2_project_dir: Optional[Path] = None,
    ) -> None:
        if weights_path is not None:
            self.weights_path = Path(weights_path)
            self.model_name = self.weights_path.stem
        else:
            self.model_name, self.weights_path = resolve_pcd_weights("LeafPriority.pt")
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.conf_thres = conf_thres
        self.iou_thres = iou_thres
        self.img_size = img_size
        self.lm2_project_dir = Path(lm2_project_dir) if lm2_project_dir else None
        self.model = None
        self.stride = 32
        self._load_model()

    def _load_model(self) -> None:
        """Loads YOLOv5 LeafPriority.pt model if available in current environment."""
        if not self.weights_path.exists():
            logger.debug(f"Stage 1 YOLOv5 weights not found at {self.weights_path}; using precomputed/fallback mode.")
            return

        pcd_dir = LM2_DIR / "component_detector"
        if not pcd_dir.exists():
            return

        try:
            import torch
            if str(pcd_dir) not in sys.path:
                sys.path.insert(0, str(pcd_dir))
            from models.common import DetectMultiBackend
            from utils.torch_utils import select_device

            dev = select_device(self.device if self.device != "cuda" else "0")
            self.model = DetectMultiBackend(str(self.weights_path), device=dev, dnn=False)
            self.stride = int(max(getattr(self.model, "stride", 32), 32))
            logger.info(f"Loaded Stage 1 LM2 YOLOv5 detector ({self.model_name}) from {self.weights_path}")
        except Exception as e:
            logger.debug(f"Stage 1 YOLOv5 runtime not loaded in current interpreter ({e}).")
            self.model = None

    def predict(
        self,
        image_bgr: np.ndarray,
        catalog_number: Optional[str] = None,
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        """
        Runs Stage 1 full-sheet organ detection on `image_bgr`.
        Returns list of ((ymin, xmin, ymax, xmax), score, class_id) in 11-class LM2 taxonomy.
        """
        sheet_h, sheet_w = image_bgr.shape[:2]

        # 1. Check precomputed LM2 Plant_Components YOLO labels if available
        if self.lm2_project_dir and catalog_number:
            label_candidates = [
                self.lm2_project_dir / "Data" / "precomputed_annotations" / "Plant_Components" / "labels" / f"{catalog_number}.txt",
                self.lm2_project_dir / "Plant_Components" / "labels" / f"{catalog_number}.txt",
            ]
            for lf in label_candidates:
                if lf.exists():
                    dets = load_precomputed_pcd_labels(
                        lf, sheet_w, sheet_h, min_confidence=self.conf_thres
                    )
                    if dets:
                        return dets

        # 2. Run live YOLOv5 LeafPriority.pt inference if model is loaded
        if self.model is None:
            return []

        try:
            import torch
            from utils.augmentations import letterbox
            from utils.general import non_max_suppression, scale_coords

            img = letterbox(image_bgr, self.img_size, stride=self.stride, auto=self.model.pt)[0]
            img = img.transpose((2, 0, 1))[::-1]
            img = np.ascontiguousarray(img)

            im = torch.from_numpy(img).to(self.model.device)
            im = im.float() / 255.0
            if len(im.shape) == 3:
                im = im[None]

            with torch.no_grad():
                pred = self.model(im, augment=False, visualize=False)
                pred = non_max_suppression(
                    pred,
                    conf_thres=self.conf_thres,
                    iou_thres=self.iou_thres,
                    max_det=1000,
                )

            results: List[Tuple[Tuple[int, int, int, int], float, int]] = []
            for det in pred:
                if len(det):
                    det[:, :4] = scale_coords(im.shape[2:], det[:, :4], image_bgr.shape).round()
                    for *xyxy, conf, cls in reversed(det):
                        xmin, ymin, xmax, ymax = (
                            max(0, int(xyxy[0].item())),
                            max(0, int(xyxy[1].item())),
                            min(sheet_w, int(xyxy[2].item())),
                            min(sheet_h, int(xyxy[3].item())),
                        )
                        if ymax > ymin and xmax > xmin:
                            results.append(((ymin, xmin, ymax, xmax), float(conf.item()), int(cls.item())))
            return results
        except Exception as e:
            logger.debug(f"Stage 1 YOLOv5 inference error: {e}")
            return []


# =============================================================================
# Stage 2 Cropped PointRend Model Inference Engine
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
        if self.weights_path.exists():
            checkpointer = DetectionCheckpointer(model)
            checkpointer.load(cfg.MODEL.WEIGHTS)
        else:
            logger.warning(f"Weights file not found at {self.weights_path}, running with uninitialized weights for testing.")

        self.cfg = cfg
        self.model = model
        self.aug = T.ResizeShortestEdge([cfg.INPUT.MIN_SIZE_TEST, cfg.INPUT.MIN_SIZE_TEST], cfg.INPUT.MAX_SIZE_TEST)
        logger.info(f"Loaded PointRend model from {self.weights_path} on {self.device}")

    def predict(self, image_bgr: np.ndarray) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """
        Runs inference on BGR image (cropped leaf patch or full sheet).
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


def is_stock_pointrend_weights(weights_path: Union[str, Path]) -> bool:
    """
    Returns True when using stock un-finetuned LeafMachine2 PointRend weights
    (e.g. model_final.pth from leaf_seg_rcnn_v2_1) rather than a Packera-finetuned checkpoint.
    """
    name = Path(weights_path).name.lower()
    if not name:
        return True
    return "finetuned" not in name and "packera" not in name


# =============================================================================
# Automated LM2-Guided SAM 2 Boundary Refinement Engine
# =============================================================================

class SAM2BoundaryRefiner:
    """
    Automated zero-shot SAM 2 boundary refiner guided by Stage 2 PointRend initial masks (m_init).

    Given a 12%-padded native-DPI leaf crop and initial binary mask `m_init`:
      1. Generates positive medial/lobe prompts inside eroded `m_init`, negative corner/concavity
         prompts outside dilated `m_init`, crop bounding box, and `(1, 256, 256)` low-res logits.
      2. Queries `SAM2ImagePredictor.predict()` deterministically (`torch.no_grad()`, fixed seeds).
      3. Validates topological connectedness (`is_topologically_connected_mask`); falls back cleanly
         to `m_init` if SAM 2 weights are unavailable or if the refined mask violates connectedness.
    """

    def __init__(
        self,
        checkpoint_path: Optional[Union[str, Path]] = None,
        model_cfg: str = "configs/sam2.1/sam2.1_hiera_l.yaml",
        device: str = "cuda",
        predictor: Optional[Any] = None,
        seed: int = 42,
    ) -> None:
        self.checkpoint_path = (
            Path(checkpoint_path)
            if checkpoint_path is not None
            else (PROJECT_ROOT / "models" / "sam2.1_hiera_large.pt")
        )
        self.model_cfg = model_cfg
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.seed = int(seed)
        self.predictor = predictor
        if self.predictor is None:
            self._load_predictor()

    def _load_predictor(self) -> None:
        """Loads SAM2ImagePredictor if checkpoint and sam2 package are available."""
        if not self.checkpoint_path.exists():
            logger.debug(
                f"SAM 2 checkpoint not found at {self.checkpoint_path}; "
                "SAM2BoundaryRefiner will fall back cleanly to PointRend m_init."
            )
            return

        try:
            import torch
            from sam2.build_sam import build_sam2
            from sam2.sam2_image_predictor import SAM2ImagePredictor

            torch.manual_seed(self.seed)
            if torch.cuda.is_available():
                torch.cuda.manual_seed_all(self.seed)

            sam2_model = build_sam2(
                self.model_cfg,
                str(self.checkpoint_path),
                device=self.device,
            )
            sam2_model.eval()
            self.predictor = SAM2ImagePredictor(sam2_model)
            logger.info(f"Loaded SAM 2 boundary refiner from {self.checkpoint_path} on {self.device}")
        except Exception as e:
            logger.debug(f"SAM 2 runtime unavailable ({e}); falling back to PointRend m_init.")
            self.predictor = None

    @property
    def is_available(self) -> bool:
        return self.predictor is not None

    def refine_crop_mask(
        self,
        crop_bgr: np.ndarray,
        m_init: np.ndarray,
    ) -> np.ndarray:
        """
        Refines an initial PointRend binary mask `m_init` on `crop_bgr` using automated
        SAM 2 prompting. Falls back cleanly to `m_init` if SAM 2 is unavailable or if
        the refined candidate violates topological connectedness.
        """
        if m_init is None or m_init.size == 0:
            return np.zeros((0, 0), dtype=np.uint8)

        fallback_mask = (m_init > 0).astype(np.uint8) * 255
        if np.count_nonzero(fallback_mask) == 0 or crop_bgr is None or crop_bgr.size == 0:
            return fallback_mask

        if self.predictor is None:
            return fallback_mask

        prompt_bundle = generate_sam2_prompts_from_lm2_mask(fallback_mask)
        if prompt_bundle["point_coords"].shape[0] == 0 or prompt_bundle["positive_points"].shape[0] == 0:
            return fallback_mask

        try:
            np.random.seed(self.seed)
            try:
                import torch
                torch.manual_seed(self.seed)
                grad_ctx = torch.no_grad()
            except ImportError:
                from contextlib import nullcontext
                grad_ctx = nullcontext()

            crop_rgb = (
                cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGB)
                if (crop_bgr.ndim == 3 and crop_bgr.shape[2] == 3)
                else crop_bgr
            )

            with grad_ctx:
                self.predictor.set_image(crop_rgb)
                masks, scores, _ = self.predictor.predict(
                    point_coords=prompt_bundle["point_coords"],
                    point_labels=prompt_bundle["point_labels"],
                    box=prompt_bundle["box"],
                    mask_input=prompt_bundle["mask_input"],
                    multimask_output=True,
                )

            if masks is None or len(masks) == 0:
                return fallback_mask

            scores_arr = np.asarray(scores, dtype=np.float32).reshape(-1)
            order = np.argsort(scores_arr)[::-1]
            target_h, target_w = fallback_mask.shape[:2]

            for idx in order:
                cand = np.asarray(masks[int(idx)])
                if cand.ndim == 3:
                    cand = cand[0]
                cand_u8 = (cand > 0.0).astype(np.uint8) * 255
                if cand_u8.shape[0] != target_h or cand_u8.shape[1] != target_w:
                    cand_u8 = cv2.resize(cand_u8, (target_w, target_h), interpolation=cv2.INTER_NEAREST)

                if is_topologically_connected_mask(cand_u8, reference_mask=fallback_mask):
                    return cand_u8

            logger.debug("SAM 2 candidate masks violated topological connectedness; falling back to m_init.")
            return fallback_mask
        except Exception as e:
            logger.debug(f"SAM 2 boundary refinement failed ({e}); falling back to m_init.")
            return fallback_mask


def refine_mask_with_sam2(
    crop_bgr: np.ndarray,
    m_init: np.ndarray,
    predictor: Optional[Any] = None,
    checkpoint_path: Optional[Union[str, Path]] = None,
    model_cfg: str = "configs/sam2.1/sam2.1_hiera_l.yaml",
    device: str = "cuda",
) -> np.ndarray:
    """Functional helper for automated LM2 -> SAM 2 mask boundary refinement."""
    refiner = SAM2BoundaryRefiner(
        checkpoint_path=checkpoint_path,
        model_cfg=model_cfg,
        device=device,
        predictor=predictor,
    )
    return refiner.refine_crop_mask(crop_bgr, m_init)


# =============================================================================
# Core Pipeline Execution Runner
# =============================================================================

class SegmentAndExtractPipeline:
    """Unified two-stage LM2 detection, cropped PointRend segmentation, and extraction runner."""

    def __init__(
        self,
        vouchers_csv: Union[str, Path],
        model_weights: Union[str, Path],
        output_dir: Union[str, Path],
        device: str = "cuda",
        min_solidity: Optional[float] = None,
        min_ucs: Optional[float] = None,
        score_thresh: float = 0.40,
        num_ruler_workers: int = 4,
        force: bool = False,
        config: Optional[PipelineConfig] = None,
        lm2_project_dir: Optional[Union[str, Path]] = None,
        stage1_weights: Optional[Union[str, Path]] = None,
        ruler_confidence_thresh: float = 0.50,
        sam2_refine: bool = False,
        sam2_checkpoint: Optional[Union[str, Path]] = None,
        sam2_config: str = "configs/sam2.1/sam2.1_hiera_l.yaml",
        sam2_predictor: Optional[Any] = None,
    ) -> None:
        self.vouchers_csv = Path(vouchers_csv)
        self.model_weights = Path(model_weights)
        self.output_dir = Path(output_dir).resolve()
        if self.output_dir.name != "data":
            self.output_dir = self.output_dir / "data"
        self.device = device
        self.config = config or PipelineConfig.from_yaml()

        if min_solidity is not None:
            self.min_solidity = min_solidity
        else:
            self.min_solidity = self.config.thresholds.min_solidity

        if min_ucs is not None:
            self.min_ucs = min_ucs
        else:
            self.min_ucs = self.config.thresholds.min_ucs

        dissection_cfg = getattr(self.config.thresholds, "dissection", None) or getattr(self.config.thresholds, "solidity", None)
        if dissection_cfg is not None:
            self.min_dissected = getattr(dissection_cfg, "min_solidity_dissected", getattr(dissection_cfg, "min_dissected", 0.50))
        else:
            self.min_dissected = 0.50

        self.score_thresh = score_thresh
        self.num_ruler_workers = num_ruler_workers
        self.force = force
        self.ruler_confidence_thresh = ruler_confidence_thresh

        # Resolve LM2 project directory and build canonical LM2 runtime config
        if lm2_project_dir is not None:
            self.lm2_project_dir: Optional[Path] = Path(lm2_project_dir)
        else:
            cand_lm2_dir = self.output_dir.parent / "LM2_Project"
            self.lm2_project_dir = cand_lm2_dir if cand_lm2_dir.exists() else None

        self.lm2_config = generate_high_performance_config(
            pipeline_cfg=self.config,
            output_dir=str(self.lm2_project_dir / "Data" / "output") if self.lm2_project_dir else None,
        )

        # Output artifact directories
        self.masks_dir = self.output_dir / "masks"
        self.contours_dir = self.output_dir / "contours"
        self.tables_dir = self.output_dir / "tables"
        self.cropped_patches_dir = self.output_dir / "cropped_patches"
        self.basal_leaves_raw_dir = self.cropped_patches_dir / "basal_leaves_raw"
        self.capitula_patches_dir = self.cropped_patches_dir / "capitula"
        self.rosettes_dense_dir = self.cropped_patches_dir / "rosettes_dense"

        for d in [
            self.masks_dir,
            self.contours_dir,
            self.tables_dir,
            self.basal_leaves_raw_dir,
            self.capitula_patches_dir,
            self.rosettes_dense_dir,
        ]:
            d.mkdir(parents=True, exist_ok=True)

        # Preload LM2 38-class ruler scale calibrations (keeping Hough as secondary fallback)
        self.ruler_calibrations: Dict[str, float] = {}
        if self.lm2_project_dir and self.lm2_project_dir.exists():
            self.ruler_calibrations.update(
                load_ruler_calibrations(self.lm2_project_dir, min_confidence=self.ruler_confidence_thresh)
            )
        if self.output_dir.exists():
            self.ruler_calibrations.update(
                load_ruler_calibrations(self.output_dir, min_confidence=self.ruler_confidence_thresh)
            )

        self.ruler_executor = concurrent.futures.ThreadPoolExecutor(max_workers=self.num_ruler_workers)

        # Stage 1 Full-Sheet Organ Detector (LeafPriority.pt, 11 classes)
        self.stage1_detector = LM2Stage1Detector(
            weights_path=stage1_weights,
            device=self.device,
            conf_thres=min(self.score_thresh, 0.25),
            lm2_project_dir=self.lm2_project_dir,
        )

        # Stage 2 Cropped Leaf PointRend Segmentation Engine
        self.model_engine = PointRendInferenceEngine(
            weights_path=self.model_weights,
            device=self.device,
            score_thresh=self.score_thresh
        )
        self.is_stock_pointrend = is_stock_pointrend_weights(self.model_weights)
        self.sam2_refine = bool(sam2_refine)
        self.sam2_checkpoint = (
            Path(sam2_checkpoint)
            if sam2_checkpoint is not None
            else (PROJECT_ROOT / "models" / "sam2.1_hiera_large.pt")
        )
        self.sam2_config = sam2_config
        self.sam2_refiner = SAM2BoundaryRefiner(
            checkpoint_path=self.sam2_checkpoint,
            model_cfg=self.sam2_config,
            device=self.device,
            predictor=sam2_predictor,
        )
        self.voucher_repro_records: Dict[str, Dict[str, Any]] = {}
        self.voucher_pheno_records: Dict[str, Dict[str, Any]] = {}

        # Preload catalogNumber -> scientificName mapping if available
        self.vouchers_taxa: Dict[str, str] = {}
        if self.vouchers_csv.exists() and self.vouchers_csv.stat().st_size > 0:
            try:
                vdf = pd.read_csv(self.vouchers_csv)
                if "catalogNumber" in vdf.columns:
                    for _, r in vdf.iterrows():
                        cat = str(r["catalogNumber"]).strip()
                        tax = str(r.get("scientificName", r.get("species_raw", ""))).strip()
                        if tax:
                            self.vouchers_taxa[cat] = tax
            except Exception:
                pass

    def _should_refine_with_sam2(self, m_init: np.ndarray) -> bool:
        """
        Determines whether to invoke automated LM2-guided SAM 2 boundary refinement
        on an initial PointRend leaf mask `m_init`.
        Triggered when `--sam2-refine` is enabled or when using stock un-finetuned
        PointRend weights on dissected leaves.
        """
        if self.sam2_refiner is None or not self.sam2_refiner.is_available:
            return False
        if self.sam2_refine:
            return True
        if not self.is_stock_pointrend:
            return False
        if m_init is None or m_init.size == 0 or np.count_nonzero(m_init) == 0:
            return False
        return bool(is_botanical_dissection((m_init > 0).astype(np.uint8) * 255))

    def shutdown(self) -> None:
        """Cleans up background thread pools."""
        self.ruler_executor.shutdown(wait=False)

    def is_voucher_completed(self, catalog_number: str) -> bool:
        """
        Checks whether output artifacts already exist for a voucher specimen.
        Verifies:
          1. Standardized contour CSV: data/contours/{catalogNumber}_leaf*.csv (size > 0)
          2. Binary silhouette mask: data/masks/{catalogNumber}_leaf*.png (size > 0)
             or within tier1_pristine / tier2_reflected subdirectories.
        """
        contours = list(self.contours_dir.glob(f"{catalog_number}_leaf*.csv"))
        if contours and any(c.is_file() and c.stat().st_size > 0 for c in contours):
            return True

        masks = list(self.masks_dir.glob(f"{catalog_number}_leaf*.png"))
        if masks and any(m.is_file() and m.stat().st_size > 0 for m in masks):
            return True

        for sub in ["tier1_pristine", "tier1_intact", "tier2_reflected"]:
            sub_dir = self.masks_dir / sub
            if sub_dir.exists():
                sub_masks = list(sub_dir.glob(f"{catalog_number}_leaf*.png"))
                if sub_masks and any(m.is_file() and m.stat().st_size > 0 for m in sub_masks):
                    return True

        return False

    @staticmethod
    def _atomic_save_csv(df: pd.DataFrame, target_path: Path) -> None:
        """Atomically persists DataFrame to CSV via temporary file replacement."""
        target_path.parent.mkdir(parents=True, exist_ok=True)
        temp_file = tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            dir=target_path.parent,
            suffix=".tmp",
            encoding="utf-8",
        )
        temp_path = Path(temp_file.name)
        try:
            df.to_csv(temp_file, index=False, encoding="utf-8")
            temp_file.flush()
            os.fsync(temp_file.fileno())
            temp_file.close()
            temp_path.replace(target_path)
        except Exception:
            if temp_path.exists():
                try:
                    temp_path.unlink()
                except Exception:
                    pass
            raise

    def process_voucher(
        self,
        catalog_number: str,
        image_path: Path,
        taxon: Optional[str] = None
    ) -> Tuple[List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Processes a single herbarium voucher specimen sheet using Two-Stage LM2 + PointRend."""
        if taxon is None:
            taxon = self.vouchers_taxa.get(catalog_number, None)

        if not image_path.exists():
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "image_not_found",
                "details": f"File does not exist at {image_path}"
            }
            self.voucher_pheno_records[catalog_number] = classify_voucher_phenology(catalog_number, [])
            return [], failure_record

        img_bgr = cv2.imread(str(image_path))
        if img_bgr is None:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "corrupted_image",
                "details": "Failed to decode image with OpenCV"
            }
            self.voucher_pheno_records[catalog_number] = classify_voucher_phenology(catalog_number, [])
            return [], failure_record

        sheet_h, sheet_w = img_bgr.shape[:2]

        # 1. Ruler Scale Calibration: Prefer Stage 1 LM2 Ruler Calibration, fallback to Hough
        lm2_px_per_mm = self.ruler_calibrations.get(catalog_number)
        ruler_future = None
        if lm2_px_per_mm is not None and not np.isnan(lm2_px_per_mm) and lm2_px_per_mm > 0:
            pixels_per_mm: Optional[float] = (
                float(1.0 / lm2_px_per_mm) if lm2_px_per_mm < 1.0 else float(lm2_px_per_mm)
            )
        else:
            pixels_per_mm = None
            ruler_future = self.ruler_executor.submit(detect_ruler_scale_hough, img_bgr)

        # 2. Stage 1 Full-Sheet Organ Detection (LeafPriority.pt 11-class detector)
        stage1_detections = self.stage1_detector.predict(img_bgr, catalog_number=catalog_number)

        capitulum_detections: List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]] = []
        leaf_detections: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int, Tuple[int, int]]] = []

        if stage1_detections:
            # Stage 1 -> Stage 2 Two-Stage Flow
            for s1_item in stage1_detections:
                if len(s1_item) == 4:
                    s1_bbox, s1_mask, s1_score, s1_cls = s1_item
                else:
                    s1_bbox, s1_score, s1_cls = s1_item
                    s1_mask = None

                if s1_cls in LM2_REPRODUCTIVE_CLASS_IDS:
                    capitulum_detections.append((s1_bbox, s1_mask, float(s1_score), int(s1_cls)))
                elif s1_cls in LM2_LEAF_CLASS_IDS:
                    # Crop candidate leaf box with 12% padding and pass to Stage 2 PointRend
                    crop_bgr, crop_box = crop_leaf_bbox_with_padding(img_bgr, s1_bbox, pad_ratio=0.12)
                    if crop_bgr.size == 0:
                        continue
                    crop_preds = self.model_engine.predict(crop_bgr)
                    for pr_bbox, pr_mask, pr_score, _pr_cls in crop_preds:
                        m_refined = pr_mask
                        if self._should_refine_with_sam2(pr_mask):
                            m_refined = self.sam2_refiner.refine_crop_mask(crop_bgr, pr_mask)
                        sheet_mask, full_bbox = reproject_crop_mask_to_sheet(
                            m_refined, crop_box, (sheet_h, sheet_w)
                        )
                        leaf_detections.append((
                            full_bbox,
                            sheet_mask,
                            float(pr_score),
                            int(s1_cls),
                            (crop_box[1], crop_box[0]),
                        ))
        else:
            # Fallback when Stage 1 detector is uninitialized (e.g., direct PointRend unit tests)
            raw_detections = self.model_engine.predict(img_bgr)
            for bbox, mask, score, cls_id in raw_detections:
                if cls_id in LM2_REPRODUCTIVE_CLASS_IDS:
                    capitulum_detections.append((bbox, mask, float(score), int(cls_id)))
                elif cls_id in LM2_LEAF_CLASS_IDS:
                    sheet_mask = mask
                    full_bbox = bbox
                    if self._should_refine_with_sam2(mask):
                        crop_bgr, crop_box = crop_leaf_bbox_with_padding(img_bgr, bbox, pad_ratio=0.12)
                        cy0, cx0, cy1, cx1 = crop_box
                        if crop_bgr.size > 0 and cy1 > cy0 and cx1 > cx0:
                            crop_m_init = mask[cy0:cy1, cx0:cx1]
                            crop_m_refined = self.sam2_refiner.refine_crop_mask(crop_bgr, crop_m_init)
                            sheet_mask, full_bbox = reproject_crop_mask_to_sheet(
                                crop_m_refined, crop_box, (sheet_h, sheet_w)
                            )
                    leaf_detections.append((
                        full_bbox,
                        sheet_mask,
                        float(score),
                        int(cls_id),
                        (int(full_bbox[1]), int(full_bbox[0])),
                    ))

        # Collect asynchronous Hough scale fallback if LM2 ruler calibration was absent
        if ruler_future is not None:
            try:
                pixels_per_mm = ruler_future.result(timeout=5.0)
            except Exception as e:
                logger.debug(f"Ruler detection timeout or exception for {catalog_number}: {e}")
                pixels_per_mm = None

        is_calibrated_mm = bool(
            pixels_per_mm is not None and not np.isnan(pixels_per_mm) and float(pixels_per_mm) > 0.0
        )

        if not leaf_detections and not capitulum_detections:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "no_leaves_detected",
                "details": "Model produced 0 candidate leaf segmentations"
            }
            self.voucher_repro_records[catalog_number] = {
                "catalogNumber": catalog_number,
                "capitula_count": 0,
                "involucre_height_px": np.nan,
                "involucre_width_px": np.nan,
                "involucre_height_mm": np.nan,
                "involucre_width_mm": np.nan,
                "capitulum_aspect_ratio": np.nan,
            }
            self.voucher_pheno_records[catalog_number] = classify_voucher_phenology(catalog_number, [])
            return [], failure_record

        # 3. Route Reproductive Detections (seed_fruit_one/many, flower_one/many, bud) to Phenology Classifier
        valid_capitula: List[Dict[str, Any]] = []
        head_pheno_results: List[Dict[str, Any]] = []
        for cap_idx, (bbox, mask, score, class_id) in enumerate(capitulum_detections, 1):
            metrics = compute_capitulum_metrics(mask, bbox, pixels_per_mm)
            if metrics is not None:
                metrics["score"] = score
                metrics["class_id"] = class_id
                valid_capitula.append(metrics)

            ymin, xmin, ymax, xmax = [int(v) for v in bbox]
            y0, y1 = max(0, ymin), min(sheet_h, max(0, ymax))
            x0, x1 = max(0, xmin), min(sheet_w, max(0, xmax))
            crop_bgr = img_bgr[y0:y1, x0:x1]
            if crop_bgr.size > 0:
                cap_patch_path = self.capitula_patches_dir / f"{catalog_number}_cap{cap_idx}.jpg"
                self.capitula_patches_dir.mkdir(parents=True, exist_ok=True)
                cv2.imwrite(str(cap_patch_path), crop_bgr)

                if crop_bgr.shape[0] >= 5 and crop_bgr.shape[1] >= 5:
                    crop_rgb = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGB)
                    head_pheno = classify_single_capitulum(crop_rgb)
                    head_pheno_results.append(head_pheno)

        self.voucher_pheno_records[catalog_number] = classify_voucher_phenology(
            catalog_number, head_pheno_results
        )

        capitula_count = len(valid_capitula)
        if capitula_count > 0:
            med_h_px = float(np.median([c["involucre_height_px"] for c in valid_capitula]))
            med_w_px = float(np.median([c["involucre_width_px"] for c in valid_capitula]))
            h_mm_list = [c["involucre_height_mm"] for c in valid_capitula if not np.isnan(c["involucre_height_mm"])]
            w_mm_list = [c["involucre_width_mm"] for c in valid_capitula if not np.isnan(c["involucre_width_mm"])]
            med_h_mm = float(np.median(h_mm_list)) if h_mm_list else np.nan
            med_w_mm = float(np.median(w_mm_list)) if w_mm_list else np.nan
            med_ar = float(np.median([c["capitulum_aspect_ratio"] for c in valid_capitula]))
        else:
            med_h_px = np.nan
            med_w_px = np.nan
            med_h_mm = np.nan
            med_w_mm = np.nan
            med_ar = np.nan

        repro_summary = {
            "catalogNumber": catalog_number,
            "capitula_count": capitula_count,
            "involucre_height_px": round(med_h_px, 2) if not np.isnan(med_h_px) else np.nan,
            "involucre_width_px": round(med_w_px, 2) if not np.isnan(med_w_px) else np.nan,
            "involucre_height_mm": round(med_h_mm, 3) if not np.isnan(med_h_mm) else np.nan,
            "involucre_width_mm": round(med_w_mm, 3) if not np.isnan(med_w_mm) else np.nan,
            "capitulum_aspect_ratio": round(med_ar, 4) if not np.isnan(med_ar) else np.nan,
        }
        self.voucher_repro_records[catalog_number] = repro_summary

        if not leaf_detections:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "no_leaves_detected",
                "details": f"Model produced 0 candidate leaf segmentations ({capitula_count} capitula detected)"
            }
            return [], failure_record

        # Build detected leaf instances
        instances: List[DetectedInstance] = []
        for idx, (bbox, mask, score, class_id, crop_offset) in enumerate(leaf_detections, 1):
            instances.append(DetectedInstance(
                catalog_number=catalog_number,
                leaf_id=idx,
                bbox=bbox,
                mask=mask,
                score=score,
                class_id=class_id,
                pixels_per_mm=pixels_per_mm,
                is_calibrated_mm=is_calibrated_mm,
                crop_offset=crop_offset,
            ))

        # 4. Vectorized Centroid Spatial Clustering (DBSCAN)
        instances = cluster_plant_individuals(instances, sheet_w, sheet_h)

        # 5. Tier 4 Dense Rosette Patch Generation per plant_individual_id
        clusters_by_pid: Dict[int, List[DetectedInstance]] = defaultdict(list)
        for inst in instances:
            clusters_by_pid[inst.plant_individual_id].append(inst)
        for pid, pid_instances in sorted(clusters_by_pid.items()):
            crop_dense_rosette_tier4(
                raw_image_path=image_path,
                plant_candidates=pid_instances,
                plant_id=pid,
                output_dir=self.output_dir,
                raw_img=img_bgr,
                catalog_number=catalog_number,
            )

        extracted_records: List[Dict[str, Any]] = []
        rejected_instances: List[DetectedInstance] = []

        # 6. Geometric Quality Gatekeeper, 2-Path Extraction & Basal Leaf RGB Patch Saving
        for inst in instances:
            routing_res = classify_and_route_leaf(
                inst.mask,
                config=self.config,
                min_ucs=self.min_ucs,
                min_solidity=self.min_solidity,
                taxon=taxon,
            )
            inst.ucs_score = routing_res.ucs
            inst.solidity = routing_res.solidity
            inst.midrib_angle_deg = float(routing_res.metadata.get("midrib_angle_deg", 0.0))

            # Diagnostic classification flags
            is_folded = bool(routing_res.is_folded)
            is_dissected = bool(routing_res.is_dissected)
            reflection_applied = False
            tier: str = "Failed QC"
            assigned_tier: str = "rejected"

            # Helper to extract, homologize, and export standardized contour
            def _export_homologized(mask_proc: np.ndarray, is_reflected: bool) -> Optional[str]:
                p_base_pt = routing_res.metadata.get("p_base")
                if is_reflected or p_base_pt is None or p_base_pt == (0, 0):
                    _, _, _, _, p_base_pt = compute_geometric_metrics(mask_proc)

                cnts, _ = cv2.findContours(mask_proc, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
                homologized_c = None
                if cnts:
                    raw_c = max(cnts, key=cv2.contourArea).reshape(-1, 2)
                    homologized_c = homologize_contour_starting_point(raw_c, p_base_pt)

                return export_standardized_contour(
                    mask_proc,
                    catalog_number,
                    inst.leaf_id,
                    self.output_dir,
                    petiole_attachment_pt=p_base_pt,
                    contour=homologized_c,
                )

            # Routing Decision Tree
            # 1. Check for fold: If is_folded == True -> synthesize full blade via reflection -> route to Tier 2
            if is_folded and routing_res.assigned_tier == "tier2":
                tier = "Tier 2"
                assigned_tier = TIER_2_REFLECTED
                reflection_applied = True
                out_base = (self.output_dir / "data") if (self.output_dir.name != "data" and (self.output_dir / "data").is_dir()) else self.output_dir
                mask_dir = out_base / "masks"
                tier2_dir = mask_dir / "tier2_reflected"
                mask_dir.mkdir(parents=True, exist_ok=True)
                tier2_dir.mkdir(parents=True, exist_ok=True)

                filename = f"{catalog_number}_leaf{inst.leaf_id}.png"
                main_save_path = mask_dir / filename
                tier_save_path = tier2_dir / filename

                cv2.imwrite(str(main_save_path), routing_res.processed_mask)
                cv2.imwrite(str(tier_save_path), routing_res.processed_mask)

                contour_path = _export_homologized(routing_res.processed_mask, is_reflected=True)
                inst.mask_path = str(main_save_path)
                inst.contour_path = contour_path

            # 2. Check for pristine whole leaf: If solidity >= min_solidity and ucs >= 0.85 -> route to Tier 1
            elif routing_res.solidity >= self.min_solidity and routing_res.ucs >= self.min_ucs and routing_res.assigned_tier == "tier1" and not is_dissected:
                tier = "Tier 1"
                assigned_tier = TIER_1_DIRECT
                reflection_applied = False
                mask_path = extract_tier1_pristine(routing_res.processed_mask, catalog_number, inst.leaf_id, self.output_dir)
                contour_path = _export_homologized(routing_res.processed_mask, is_reflected=False)
                inst.mask_path = mask_path
                inst.contour_path = contour_path

            # 3. Check for botanical dissection: If solidity < min_solidity but is_botanical_dissection == True -> route to Tier 1 (Dissected)
            elif (routing_res.solidity < self.min_solidity) and is_dissected and routing_res.solidity >= self.min_dissected and routing_res.assigned_tier == "tier1":
                tier = "Tier 1 (Dissected)"
                assigned_tier = TIER_1_DIRECT
                is_dissected = True
                reflection_applied = False
                mask_path = extract_tier1_pristine(routing_res.processed_mask, catalog_number, inst.leaf_id, self.output_dir)
                contour_path = _export_homologized(routing_res.processed_mask, is_reflected=False)
                inst.mask_path = mask_path
                inst.contour_path = contour_path

            # 4. Check for partial occlusion with intact half-blade -> route to Tier 2
            elif routing_res.assigned_tier == "tier2":
                tier = "Tier 2"
                assigned_tier = TIER_2_REFLECTED
                reflection_applied = True
                out_base = (self.output_dir / "data") if (self.output_dir.name != "data" and (self.output_dir / "data").is_dir()) else self.output_dir
                mask_dir = out_base / "masks"
                tier2_dir = mask_dir / "tier2_reflected"
                mask_dir.mkdir(parents=True, exist_ok=True)
                tier2_dir.mkdir(parents=True, exist_ok=True)

                filename = f"{catalog_number}_leaf{inst.leaf_id}.png"
                main_save_path = mask_dir / filename
                tier_save_path = tier2_dir / filename

                cv2.imwrite(str(main_save_path), routing_res.processed_mask)
                cv2.imwrite(str(tier_save_path), routing_res.processed_mask)

                contour_path = _export_homologized(routing_res.processed_mask, is_reflected=True)
                inst.mask_path = str(main_save_path)
                inst.contour_path = contour_path

            # 5. Else -> log rejection reason (FAILED_SOLIDITY, IRREGULAR_FOLD, UNRESOLVED_CLUMP) and route to Failed QC
            else:
                if routing_res.is_irregular_fold:
                    rejection_reason = "IRREGULAR_FOLD"
                elif routing_res.rejection_reason in ["FAILED_SOLIDITY", "IRREGULAR_FOLD", "UNRESOLVED_CLUMP"]:
                    rejection_reason = routing_res.rejection_reason
                elif routing_res.solidity < self.min_dissected or routing_res.ucs < 0.60:
                    rejection_reason = "UNRESOLVED_CLUMP"
                elif routing_res.solidity < self.min_solidity:
                    rejection_reason = "FAILED_SOLIDITY"
                else:
                    rejection_reason = "UNRESOLVED_CLUMP"

                inst.assigned_tier = "rejected"
                inst.rejection_reason = rejection_reason
                logger.debug(
                    f"Voucher {catalog_number} leaf {inst.leaf_id} rejected by gatekeeper: {rejection_reason} "
                    f"(solidity={inst.solidity:.3f}, ucs={inst.ucs_score:.3f})"
                )
                rejected_instances.append(inst)
                continue

            inst.assigned_tier = assigned_tier

            # Save extracted basal leaf RGB crop to data/cropped_patches/basal_leaves_raw/{catalogNumber}_leaf{leaf_id}.jpg
            ly0, lx0, ly1, lx1 = [int(v) for v in inst.bbox]
            ly0, ly1 = max(0, ly0), min(sheet_h, max(0, ly1))
            lx0, lx1 = max(0, lx0), min(sheet_w, max(0, lx1))
            if ly1 <= ly0 or lx1 <= lx0:
                ys_m, xs_m = np.where(inst.mask > 0)
                if len(ys_m) > 0 and len(xs_m) > 0:
                    ly0, ly1 = int(np.min(ys_m)), min(sheet_h, int(np.max(ys_m)) + 1)
                    lx0, lx1 = int(np.min(xs_m)), min(sheet_w, int(np.max(xs_m)) + 1)
            leaf_rgb_crop = img_bgr[ly0:ly1, lx0:lx1]
            if leaf_rgb_crop.size > 0:
                self.basal_leaves_raw_dir.mkdir(parents=True, exist_ok=True)
                leaf_patch_path = self.basal_leaves_raw_dir / f"{catalog_number}_leaf{inst.leaf_id}.jpg"
                cv2.imwrite(str(leaf_patch_path), leaf_rgb_crop)

            # Record successfully extracted leaf with diagnostic classification flags
            extracted_records.append({
                "catalogNumber": catalog_number,
                "leaf_id": inst.leaf_id,
                "tier": tier,
                "solidity": round(float(inst.solidity), 4),
                "ucs": round(float(inst.ucs_score), 4),
                "is_folded": is_folded,
                "is_dissected": is_dissected,
                "reflection_applied": reflection_applied,
                "contour_path": inst.contour_path,
                "plant_individual_id": inst.plant_individual_id,
                "assigned_tier": inst.assigned_tier,
                "ucs_score": round(float(inst.ucs_score), 4),
                "midrib_angle_deg": round(float(inst.midrib_angle_deg), 2),
                "pixels_per_mm": round(float(inst.pixels_per_mm), 4) if inst.pixels_per_mm else np.nan,
                "is_calibrated_mm": bool(inst.is_calibrated_mm),
                "mask_path": inst.mask_path,
                "capitula_count": capitula_count,
                "involucre_height_px": repro_summary["involucre_height_px"],
                "involucre_width_px": repro_summary["involucre_width_px"],
                "involucre_height_mm": repro_summary["involucre_height_mm"],
                "involucre_width_mm": repro_summary["involucre_width_mm"],
                "capitulum_aspect_ratio": repro_summary["capitulum_aspect_ratio"],
            })

        # Route vouchers with 0 valid silhouettes to failed manifest
        failure_record = None
        if not extracted_records:
            primary_reason = "UNRESOLVED_CLUMP"
            if rejected_instances:
                reasons = [inst.rejection_reason for inst in rejected_instances if inst.rejection_reason]
                if reasons:
                    from collections import Counter
                    primary_reason = Counter(reasons).most_common(1)[0][0]
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": primary_reason,
                "details": f"{len(rejected_instances)} leaves detected but all rejected by geometric gatekeeping"
            }

        return extracted_records, failure_record

    def run(self, limit: Optional[int] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Executes full extraction pipeline across curated vouchers with resumption support."""
        logger.info("==================================================================")
        logger.info("Starting Step 02: Direct Segmentation & Botanical Leaf Extraction")
        logger.info(f"Vouchers Metadata: {self.vouchers_csv}")
        logger.info(f"PointRend Model: {self.model_weights}")
        logger.info(f"Output Directory: {self.output_dir}")
        logger.info(f"Device: {self.device}")
        logger.info(f"Force Overwrite: {self.force}")
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

        manifest_path = self.tables_dir / "extracted_leaf_manifest.csv"
        legacy_manifest_path = self.tables_dir / "extracted_leaves_manifest.csv"
        failed_qc_path = self.tables_dir / "failed_qc_vouchers.csv"
        repro_manifest_path = self.tables_dir / "voucher_reproductive_metrics.csv"
        pheno_manifest_path = self.tables_dir / "voucher_phenological_states.csv"

        existing_extracted_by_cat: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        existing_failed_by_cat: Dict[str, Dict[str, Any]] = {}
        existing_repro_by_cat: Dict[str, Dict[str, Any]] = {}
        existing_pheno_by_cat: Dict[str, Dict[str, Any]] = {}

        if not self.force:
            for path in [manifest_path, legacy_manifest_path]:
                if path.exists() and path.stat().st_size > 0:
                    try:
                        prev_manifest = pd.read_csv(path)
                        if not prev_manifest.empty and "catalogNumber" in prev_manifest.columns:
                            for rec in prev_manifest.to_dict(orient="records"):
                                existing_extracted_by_cat[str(rec["catalogNumber"]).strip()].append(rec)
                        break  # Found a valid manifest, don't need to check the legacy one
                    except Exception as e:
                        logger.warning(f"Could not load previous manifest at {path}: {e}")

            if failed_qc_path.exists() and failed_qc_path.stat().st_size > 0:
                try:
                    prev_failed = pd.read_csv(failed_qc_path)
                    if not prev_failed.empty and "catalogNumber" in prev_failed.columns:
                        for rec in prev_failed.to_dict(orient="records"):
                            existing_failed_by_cat[str(rec["catalogNumber"]).strip()] = rec
                except Exception as e:
                    logger.warning(f"Could not load previous failed QC table at {failed_qc_path}: {e}")

            if repro_manifest_path.exists() and repro_manifest_path.stat().st_size > 0:
                try:
                    prev_repro = pd.read_csv(repro_manifest_path)
                    if not prev_repro.empty and "catalogNumber" in prev_repro.columns:
                        for rec in prev_repro.to_dict(orient="records"):
                            existing_repro_by_cat[str(rec["catalogNumber"]).strip()] = rec
                except Exception as e:
                    logger.warning(f"Could not load previous reproductive table at {repro_manifest_path}: {e}")

            if pheno_manifest_path.exists() and pheno_manifest_path.stat().st_size > 0:
                try:
                    prev_pheno = pd.read_csv(pheno_manifest_path)
                    if not prev_pheno.empty and "catalogNumber" in prev_pheno.columns:
                        for rec in prev_pheno.to_dict(orient="records"):
                            existing_pheno_by_cat[str(rec["catalogNumber"]).strip()] = rec
                except Exception as e:
                    logger.warning(f"Could not load previous phenology table at {pheno_manifest_path}: {e}")

        all_extracted_records: List[Dict[str, Any]] = []
        all_failed_records: List[Dict[str, Any]] = []

        total_vouchers = len(df_vouchers)
        processed_count = 0
        skipped_count = 0

        for i, (_, row) in enumerate(df_vouchers.iterrows(), 1):
            cat_num = str(row["catalogNumber"]).strip()
            raw_img_path = Path(str(row["image_path"]).strip())
            if not raw_img_path.is_absolute():
                raw_img_path = PROJECT_ROOT / raw_img_path

            # Check if specimen already processed and resumption is enabled
            if not self.force and self.is_voucher_completed(cat_num):
                skipped_count += 1
                logger.info(
                    f"Skipping completed voucher: {cat_num} (artifacts exist in contours/masks). "
                    f"[Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}]"
                )
                if cat_num in existing_extracted_by_cat:
                    all_extracted_records.extend(existing_extracted_by_cat[cat_num])
                elif cat_num in existing_failed_by_cat:
                    all_failed_records.append(existing_failed_by_cat[cat_num])

                if cat_num in existing_repro_by_cat:
                    self.voucher_repro_records[cat_num] = existing_repro_by_cat[cat_num]
                if cat_num in existing_pheno_by_cat:
                    self.voucher_pheno_records[cat_num] = existing_pheno_by_cat[cat_num]
                continue

            processed_count += 1
            taxon = str(row.get("scientificName", row.get("species_raw", ""))).strip()
            extracted, failure = self.process_voucher(cat_num, raw_img_path, taxon=taxon)
            all_extracted_records.extend(extracted)
            if failure:
                all_failed_records.append(failure)

            if processed_count % 25 == 0 or (processed_count + skipped_count) == total_vouchers:
                logger.info(
                    f"Progress: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}] "
                    f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC)."
                )

        logger.info(
            f"Extraction Summary: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}]"
        )

        # Build manifest tables with diagnostic classification schema
        manifest_cols = [
            "catalogNumber", "leaf_id", "tier", "solidity", "ucs",
            "is_folded", "is_dissected", "reflection_applied", "contour_path"
        ]
        extracted_df = pd.DataFrame(all_extracted_records)
        if extracted_df.empty:
            extracted_df = pd.DataFrame(columns=manifest_cols + [
                "plant_individual_id", "assigned_tier", "ucs_score",
                "midrib_angle_deg", "pixels_per_mm", "is_calibrated_mm", "mask_path",
                "capitula_count", "involucre_height_px", "involucre_width_px",
                "involucre_height_mm", "involucre_width_mm", "capitulum_aspect_ratio"
            ])
        else:
            # Normalize assigned_tier to canonical contract ("Tier_1_Direct", "Tier_2_Reflected")
            if "assigned_tier" in extracted_df.columns:
                extracted_df["assigned_tier"] = extracted_df["assigned_tier"].apply(normalize_assigned_tier)
            # Backfill any missing diagnostic fields for legacy resumption entries
            if "tier" not in extracted_df.columns:
                extracted_df["tier"] = extracted_df.get("assigned_tier", TIER_1_DIRECT).apply(
                    lambda x: "Tier 2" if "2" in str(x).lower() else "Tier 1"
                )
            if "ucs" not in extracted_df.columns:
                extracted_df["ucs"] = extracted_df.get("ucs_score", np.nan)
            if "is_folded" not in extracted_df.columns:
                extracted_df["is_folded"] = False
            if "is_dissected" not in extracted_df.columns:
                extracted_df["is_dissected"] = False
            if "reflection_applied" not in extracted_df.columns:
                extracted_df["reflection_applied"] = extracted_df.get("assigned_tier", TIER_1_DIRECT).apply(
                    lambda x: True if "2" in str(x).lower() else False
                )
            if "is_calibrated_mm" not in extracted_df.columns:
                extracted_df["is_calibrated_mm"] = extracted_df.get("pixels_per_mm", pd.Series(dtype=float)).apply(
                    lambda x: bool(pd.notna(x) and float(x) > 0)
                )
            for c in manifest_cols:
                if c not in extracted_df.columns:
                    extracted_df[c] = np.nan

            ordered_cols = manifest_cols + [c for c in extracted_df.columns if c not in manifest_cols]
            extracted_df = extracted_df[ordered_cols]

        failed_df = pd.DataFrame(all_failed_records)
        if failed_df.empty:
            failed_df = pd.DataFrame(columns=[
                "catalogNumber", "image_path", "failure_reason", "details"
            ])

        repro_df = pd.DataFrame(list(self.voucher_repro_records.values()))
        if repro_df.empty:
            repro_df = pd.DataFrame(columns=[
                "catalogNumber", "capitula_count", "involucre_height_px",
                "involucre_width_px", "involucre_height_mm", "involucre_width_mm",
                "capitulum_aspect_ratio"
            ])
        else:
            repro_cols = [
                "catalogNumber", "capitula_count", "involucre_height_px",
                "involucre_width_px", "involucre_height_mm", "involucre_width_mm",
                "capitulum_aspect_ratio"
            ]
            repro_df = repro_df[[c for c in repro_cols if c in repro_df.columns]]

        pheno_df = pd.DataFrame(list(self.voucher_pheno_records.values()))
        if pheno_df.empty:
            pheno_df = pd.DataFrame(columns=[
                "catalogNumber", "total_capitula", "n_anthesis", "n_fruit",
                "n_bud", "voucher_phenological_state", "dominant_pappus_ratio"
            ])
        else:
            pheno_cols = [
                "catalogNumber", "total_capitula", "n_anthesis", "n_fruit",
                "n_bud", "voucher_phenological_state", "dominant_pappus_ratio"
            ]
            pheno_df = pheno_df[[c for c in pheno_cols if c in pheno_df.columns]]

        # Atomically save manifests via temporary files
        self._atomic_save_csv(extracted_df, manifest_path)
        self._atomic_save_csv(extracted_df, legacy_manifest_path)
        self._atomic_save_csv(failed_df, failed_qc_path)
        self._atomic_save_csv(repro_df, repro_manifest_path)
        self._atomic_save_csv(pheno_df, pheno_manifest_path)

        logger.info("==================================================================")
        logger.info(f"Extraction completed!")
        logger.info(f"  Extracted Leaves: {len(extracted_df)} -> {manifest_path}")
        logger.info(f"  Legacy Manifest: {len(extracted_df)} -> {legacy_manifest_path}")
        logger.info(f"  Failed QC Vouchers: {len(failed_df)} -> {failed_qc_path}")
        logger.info(f"  Reproductive Metrics: {len(repro_df)} -> {repro_manifest_path}")
        logger.info(f"  Phenological States: {len(pheno_df)} -> {pheno_manifest_path}")
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
    cfg = PipelineConfig.from_yaml()
    parser = argparse.ArgumentParser(
        description="Direct PointRend segmentation, 2-path leaf extraction, and contour coordinate export."
    )
    parser.add_argument(
        "--vouchers",
        type=Path,
        default=cfg.paths.curated_vouchers_csv,
        help="Path to curated vouchers CSV containing catalogNumber and image_path."
    )
    parser.add_argument(
        "--model-weights",
        type=Path,
        default=cfg.segmentation.model_weights,
        help="Path to fine-tuned PointRend model checkpoint."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=cfg.paths.workspace_root / "data",
        help="Root output directory for masks, contours, and manifest tables."
    )
    parser.add_argument(
        "--device",
        type=str,
        default=cfg.segmentation.device if torch_cuda_available() else "cpu",
        choices=["cuda", "cpu"],
        help="Compute device for inference (default: cuda if available, else cpu)."
    )
    parser.add_argument(
        "--min-solidity",
        type=float,
        default=cfg.thresholds.min_solidity,
        help=f"Minimum solidity for Tier 1 direct pristine leaf extraction (default: {cfg.thresholds.min_solidity})."
    )
    parser.add_argument(
        "--min-ucs",
        type=float,
        default=cfg.thresholds.min_ucs,
        help=f"Minimum Unoccluded Completeness Score for Tier 1 (default: {cfg.thresholds.min_ucs})."
    )
    parser.add_argument(
        "--score-thresh",
        type=float,
        default=cfg.segmentation.score_thresh,
        help=f"Detection confidence score threshold (default: {cfg.segmentation.score_thresh})."
    )
    parser.add_argument(
        "--sam2-refine",
        dest="sam2_refine",
        action="store_true",
        default=False,
        help="Enable automated LM2-guided SAM 2 sub-pixel boundary refinement on Stage 2 leaf crops.",
    )
    parser.add_argument(
        "--sam2-checkpoint",
        type=Path,
        default=cfg.paths.workspace_root / "models" / "sam2.1_hiera_large.pt",
        help="Path to SAM 2.1 checkpoint for automated boundary refinement.",
    )
    parser.add_argument(
        "--sam2-config",
        type=str,
        default="configs/sam2.1/sam2.1_hiera_l.yaml",
        help="SAM 2 model config path for automated boundary refinement.",
    )
    parser.add_argument(
        "--force",
        "--overwrite",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-segmentation and overwrite existing contour CSVs and masks.",
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
        score_thresh=args.score_thresh,
        force=args.force,
        sam2_refine=args.sam2_refine,
        sam2_checkpoint=args.sam2_checkpoint,
        sam2_config=args.sam2_config,
    )
    pipeline.run(limit=args.limit)


if __name__ == "__main__":
    main()
