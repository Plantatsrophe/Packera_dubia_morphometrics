#!/usr/bin/env python3
"""
===============================================================================
Module: lm2_inference_engines.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Stage 1 YOLOv5/PCD organ detection, Stage 2 cropped PointRend instance
    segmentation, SAM 2 boundary refinement, and two-stage crop-to-sheet
    mask reprojection engines for Step 02.
===============================================================================
"""

from __future__ import annotations

import gc
import logging
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.configure_leafmachine2 import resolve_pcd_weights
from scripts.vision.lm2_data_loader import (
    LM2_LEAF_CLASS_IDS,
    LM2_PCD_CLASSES,
    LM2_REPRODUCTIVE_CLASS_IDS,
    load_precomputed_pcd_labels,
)
from scripts.vision.lm2_geometry_utils import (
    bridge_dissected_rachis_components,
    generate_sam2_prompts_from_lm2_mask,
    is_topologically_connected_mask,
)

LM2_DIR = PROJECT_ROOT / "LeafMachine2" / "leafmachine2"
if LM2_DIR.exists() and str(LM2_DIR) not in sys.path:
    sys.path.insert(0, str(LM2_DIR))

DETECTRON2_PATH = LM2_DIR / "segmentation" / "detectron2"
if DETECTRON2_PATH.exists() and str(DETECTRON2_PATH) not in sys.path:
    sys.path.insert(0, str(DETECTRON2_PATH))

logger = logging.getLogger("SegmentAndExtract")

__all__ = [
    "DetectedInstance",
    "_VoucherFinalizeBundle",
    "crop_leaf_bbox_with_padding",
    "reproject_crop_mask_to_sheet",
    "compute_capitulum_metrics",
    "LM2Stage1Detector",
    "PointRendInferenceEngine",
    "torch_cuda_available",
    "is_stock_pointrend_weights",
    "SAM2BoundaryRefiner",
    "refine_mask_with_sam2",
    "_extract_native_lm2_petiole_anchor",
    "_merge_and_bridge_dissected_crop_mask",
]


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
    mask_origin: Tuple[int, int] = (0, 0)  # (x0, y0) origin of mask array on full sheet
    already_trimmed: bool = False
    trim_meta: Optional[Dict[str, Any]] = None
    pcd_bbox: Optional[Tuple[int, int, int, int]] = None
    circularity: float = 0.0
    extent: float = 0.0
    blade_area_cm2: float = 0.0
    is_severed_lobe: bool = False
    mask_path: Optional[str] = None
    contour_path: Optional[str] = None
    rejection_reason: Optional[str] = None


@dataclass
class _VoucherFinalizeBundle:
    """Thread-safe per-voucher post-processing result bundle for wave-pipelined execution."""
    catalog_number: str
    image_path: str
    extracted_records: List[Dict[str, Any]]
    failure_record: Optional[Dict[str, Any]]
    repro_record: Optional[Dict[str, Any]] = None
    pheno_record: Optional[Dict[str, Any]] = None
    crowns_by_pid: Dict[int, Tuple[float, float]] = field(default_factory=dict)
    rosette_centroids_by_pid: Dict[int, Tuple[float, float]] = field(default_factory=dict)
    by_plant_petiole: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    contour_coords_by_leaf: Dict[int, List[List[float]]] = field(default_factory=dict)
    telemetry_deltas: Dict[str, int] = field(default_factory=dict)


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

    _MODEL_CACHE: Dict[Tuple[str, str], Tuple[Any, int]] = {}

    def _load_model(self) -> None:
        """Loads YOLOv5 LeafPriority.pt model if available in current environment."""
        if not self.weights_path.exists():
            logger.debug(f"Stage 1 YOLOv5 weights not found at {self.weights_path}; using precomputed/fallback mode.")
            return

        pcd_dir = LM2_DIR / "component_detector"
        if not pcd_dir.exists():
            return

        cache_key = (str(self.weights_path.resolve()), str(self.device))
        if cache_key in LM2Stage1Detector._MODEL_CACHE:
            self.model, self.stride = LM2Stage1Detector._MODEL_CACHE[cache_key]
            return

        try:
            import torch
            if str(pcd_dir) not in sys.path:
                sys.path.insert(0, str(pcd_dir))
            from models.common import DetectMultiBackend
            from utils.torch_utils import select_device

            if torch.cuda.is_available():
                torch.backends.cudnn.benchmark = False
                torch.backends.cuda.matmul.allow_tf32 = True
                torch.backends.cudnn.allow_tf32 = True
            dev = select_device(self.device if self.device != "cuda" else "0")
            use_fp16 = bool(getattr(dev, "type", "cpu") != "cpu" and torch.cuda.is_available())
            self.model = DetectMultiBackend(str(self.weights_path), device=dev, dnn=False, fp16=use_fp16)
            self.stride = int(max(getattr(self.model, "stride", 32), 32))
            LM2Stage1Detector._MODEL_CACHE[cache_key] = (self.model, self.stride)
            logger.info(
                f"Loaded Stage 1 LM2 YOLOv5 detector ({self.model_name}, fp16={getattr(self.model, 'fp16', False)}) "
                f"from {self.weights_path}"
            )
        except Exception as e:
            logger.debug(f"Stage 1 YOLOv5 runtime not loaded in current interpreter ({e}).")
            self.model = None

    def prepare_letterbox(self, image_bgr: np.ndarray) -> Optional[np.ndarray]:
        """Precomputes the resized letterbox BGR image on a background CPU worker thread."""
        if self.model is None or image_bgr is None or image_bgr.size == 0:
            return None
        try:
            from utils.augmentations import letterbox
            return letterbox(image_bgr, self.img_size, stride=self.stride, auto=self.model.pt)[0]
        except Exception:
            return None

    def predict(
        self,
        image_bgr: np.ndarray,
        catalog_number: Optional[str] = None,
        use_cache: bool = False,
        precomputed_letterbox: Optional[np.ndarray] = None,
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        """
        Runs Stage 1 full-sheet organ detection on `image_bgr`.
        Returns list of ((ymin, xmin, ymax, xmax), score, class_id) in 11-class LM2 taxonomy.
        """
        sheet_h, sheet_w = image_bgr.shape[:2]

        # 1. Check precomputed LM2 Plant_Components YOLO labels only when use_cache=True and live model is unavailable
        if use_cache and self.model is None and self.lm2_project_dir and catalog_number:
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

            if precomputed_letterbox is not None:
                img = precomputed_letterbox
            else:
                img = letterbox(image_bgr, self.img_size, stride=self.stride, auto=self.model.pt)[0]

            target_dtype = torch.float16 if getattr(self.model, "fp16", False) else torch.float32
            im = torch.from_numpy(np.ascontiguousarray(img)).to(self.model.device, non_blocking=True)
            im = im.permute(2, 0, 1).flip(0).contiguous().to(dtype=target_dtype).div_(255.0)
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

    def predict_batch(
        self,
        items: Sequence[Tuple[np.ndarray, Optional[str], Optional[np.ndarray]]],
        use_cache: bool = False,
        max_batch_size: int = 16,
    ) -> List[List[Tuple[Tuple[int, int, int, int], float, int]]]:
        """
        Runs batched FP16/FP32 Stage 1 YOLOv5 organ detection across a wave of herbarium sheets,
        grouping precomputed letterbox arrays by exact (H, W) shape so letterbox padding
        remains 100% identical to single-sheet `predict()`.
        """
        if not items:
            return []
        if (
            self.model is None
            or "predict" in self.__dict__
            or getattr(self.predict, "__func__", None) is not LM2Stage1Detector.predict
            or len(items) <= 1
        ):
            return [
                self.predict(img, catalog_number=cat, use_cache=use_cache, precomputed_letterbox=lb)
                for img, cat, lb in items
            ]

        try:
            import torch
            from utils.augmentations import letterbox
            from utils.general import non_max_suppression, scale_coords

            out_results: List[List[Tuple[Tuple[int, int, int, int], float, int]]] = [[] for _ in items]
            buckets: Dict[Tuple[int, int], List[Tuple[int, np.ndarray, np.ndarray]]] = defaultdict(list)
            for idx, (img_bgr, _cat, pre_lb) in enumerate(items):
                if img_bgr is None or img_bgr.size == 0:
                    continue
                lb = (
                    pre_lb
                    if pre_lb is not None
                    else letterbox(img_bgr, self.img_size, stride=self.stride, auto=self.model.pt)[0]
                )
                buckets[(int(lb.shape[0]), int(lb.shape[1]))].append((idx, img_bgr, lb))

            target_dtype = torch.float16 if getattr(self.model, "fp16", False) else torch.float32
            step = max(1, int(max_batch_size))
            for _shape_hw, group in buckets.items():
                for start in range(0, len(group), step):
                    chunk = group[start : start + step]
                    if len(chunk) == 1:
                        c_idx, c_img, c_lb = chunk[0]
                        out_results[c_idx] = self.predict(
                            c_img, use_cache=use_cache, precomputed_letterbox=c_lb
                        )
                        continue
                    stacked = np.stack([c[2] for c in chunk], axis=0)
                    im = torch.from_numpy(stacked).to(self.model.device, non_blocking=True)
                    im = im.permute(0, 3, 1, 2).flip(1).contiguous().to(dtype=target_dtype).div_(255.0)
                    with torch.no_grad():
                        pred = self.model(im, augment=False, visualize=False)
                        preds_nms = non_max_suppression(
                            pred,
                            conf_thres=self.conf_thres,
                            iou_thres=self.iou_thres,
                            max_det=1000,
                        )
                    for (c_idx, c_img, _), det in zip(chunk, preds_nms):
                        sh, sw = c_img.shape[:2]
                        img_dets: List[Tuple[Tuple[int, int, int, int], float, int]] = []
                        if len(det):
                            det[:, :4] = scale_coords(im.shape[2:], det[:, :4], c_img.shape).round()
                            for *xyxy, conf, cls in reversed(det):
                                xmin, ymin, xmax, ymax = (
                                    max(0, int(xyxy[0].item())),
                                    max(0, int(xyxy[1].item())),
                                    min(sw, int(xyxy[2].item())),
                                    min(sh, int(xyxy[3].item())),
                                )
                                if ymax > ymin and xmax > xmin:
                                    img_dets.append(((ymin, xmin, ymax, xmax), float(conf.item()), int(cls.item())))
                        out_results[c_idx] = img_dets
                    del im, pred, preds_nms
            return out_results
        except Exception as e:
            logger.debug(f"Batched Stage 1 YOLOv5 fallback to sequential predict ({e})")
            return [
                self.predict(img, catalog_number=cat, use_cache=use_cache, precomputed_letterbox=lb)
                for img, cat, lb in items
            ]


# =============================================================================
# Stage 2 Cropped PointRend Model Inference Engine
# =============================================================================

class PointRendInferenceEngine:
    """Encapsulates Detectron2 PointRend model loading and inference."""

    def __init__(
        self,
        weights_path: Union[str, Path],
        device: str = "cuda",
        score_thresh: float = 0.50,
        subdivision_steps: int = 6,
        subdivision_num_points: int = 2304,
    ) -> None:
        self.weights_path = Path(weights_path)
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.score_thresh = score_thresh
        self.subdivision_steps = int(subdivision_steps)
        self.subdivision_num_points = int(subdivision_num_points)
        self.predictor = None
        self.model = None
        self.cfg = None
        self.aug = None
        self.use_amp = False
        self._load_model()

    def _load_model(self) -> None:
        """Loads Detectron2 PointRend model from weights checkpoint."""
        try:
            import torch
            from detectron2.config import get_cfg
            from detectron2.projects import point_rend
            from detectron2.checkpoint import DetectionCheckpointer
            from detectron2.modeling import build_model
            import detectron2.data.transforms as T
        except ImportError as e:
            logger.warning(f"Detectron2 not available: {e}. Running in simulation/fallback mode.")
            return

        if self.device == "cuda" and torch.cuda.is_available():
            torch.backends.cudnn.benchmark = False
            torch.backends.cuda.matmul.allow_tf32 = True
            torch.backends.cudnn.allow_tf32 = True
            self.use_amp = True

        cfg = get_cfg()
        point_rend.add_pointrend_config(cfg)

        # Locate accompanying configuration YAML
        cfg_candidates = [
            self.weights_path.parent / "cfg_output.yaml",
            PROJECT_ROOT / "LeafMachine2" / "leafmachine2" / "segmentation" / "models" / "Group3_Dataset_100000_Iter_1176PTS_512Batch_smooth_l1_LR00025_BGR" / "cfg_output.yaml",
            PROJECT_ROOT / "LeafMachine2" / "leafmachine2" / "segmentation" / "models" / "leaf_seg_rcnn_v2_1" / "cfg_output.yaml",
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
        if hasattr(cfg.MODEL, "POINT_HEAD"):
            cfg.MODEL.POINT_HEAD.SUBDIVISION_STEPS = int(self.subdivision_steps)
            cfg.MODEL.POINT_HEAD.SUBDIVISION_NUM_POINTS = int(self.subdivision_num_points)

        model = build_model(cfg)
        model.eval()
        if self.weights_path.exists():
            checkpointer = DetectionCheckpointer(model)
            checkpointer.load(cfg.MODEL.WEIGHTS)
        else:
            logger.warning(f"Weights file not found at {self.weights_path}, running with uninitialized weights for testing.")

        # Keep PointRendMaskHead sub-pixel `F.grid_sample` / `point_sample` subdivision in exact FP32
        # while ResNet-101 + FPN + RPN + Box Head execute in FP16 Tensor Core autocast.
        if self.use_amp and hasattr(model, "roi_heads") and getattr(model.roi_heads, "mask_head", None) is not None:
            _orig_mask_forward = model.roi_heads.mask_head.forward

            def _fp32_mask_forward(features: Any, instances: Any) -> Any:
                with torch.amp.autocast("cuda", enabled=False):
                    if isinstance(features, dict):
                        fp32_feats = {
                            k: (v.float() if torch.is_tensor(v) and v.is_floating_point() and v.dtype != torch.float32 else v)
                            for k, v in features.items()
                        }
                    elif torch.is_tensor(features) and features.is_floating_point() and features.dtype != torch.float32:
                        fp32_feats = features.float()
                    else:
                        fp32_feats = features
                    return _orig_mask_forward(fp32_feats, instances)

            model.roi_heads.mask_head.forward = _fp32_mask_forward

        self.cfg = cfg
        self.model = model
        self.aug = T.ResizeShortestEdge([cfg.INPUT.MIN_SIZE_TEST, cfg.INPUT.MIN_SIZE_TEST], cfg.INPUT.MAX_SIZE_TEST)
        logger.info(f"Loaded PointRend model from {self.weights_path} on {self.device} (amp_fp16={self.use_amp})")

    def prepare_transformed_crop(self, image_bgr: np.ndarray) -> Optional[Tuple[np.ndarray, int, int]]:
        """
        Thread-safe CPU helper that applies Detectron2 `ResizeShortestEdge` and CHW
        transposition on worker threads so GPU `predict_batch` transfers compact `uint8` tensors.
        """
        if self.aug is None or image_bgr is None or image_bgr.size == 0:
            return None
        h, w = image_bgr.shape[:2]
        img_t = self.aug.get_transform(image_bgr).apply_image(image_bgr)
        chw_u8 = np.ascontiguousarray(np.transpose(img_t, (2, 0, 1)), dtype=np.uint8)
        return chw_u8, h, w

    def _extract_instances(self, preds: Dict[str, Any]) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """Extracts (bbox, mask_uint8, score, class_id) from a single Detectron2 prediction dict."""
        import torch
        instances = preds["instances"]
        if not instances.has("pred_masks") or len(instances) == 0:
            return []
        keep = instances.scores >= self.score_thresh
        if not bool(keep.any()):
            return []
        instances = instances[keep]
        if torch.is_tensor(instances.pred_masks):
            masks_u8 = (instances.pred_masks.to(dtype=torch.uint8) * 255).cpu().numpy()
            scores = instances.scores.cpu().numpy()
            classes = instances.pred_classes.cpu().numpy()
            boxes = instances.pred_boxes.tensor.cpu().numpy()
        else:
            instances = instances.to("cpu")
            scores = instances.scores.numpy()
            classes = instances.pred_classes.numpy()
            boxes = instances.pred_boxes.tensor.numpy()
            masks_u8 = (instances.pred_masks.numpy() > 0).astype(np.uint8) * 255

        results: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]] = []
        for i in range(len(scores)):
            box = boxes[i]
            ymin, xmin, ymax, xmax = int(box[1]), int(box[0]), int(box[3]), int(box[2])
            results.append(((ymin, xmin, ymax, xmax), masks_u8[i], float(scores[i]), int(classes[i])))
        return results

    def predict(
        self,
        image_bgr: np.ndarray,
        pretransformed: Optional[Tuple[np.ndarray, int, int]] = None,
    ) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """
        Runs inference on BGR image (cropped leaf patch or full sheet).
        Returns list of (bbox, mask, score, class_id).
        """
        if self.model is None:
            return []

        import torch
        with torch.no_grad():
            if pretransformed is not None:
                chw_u8, height, width = pretransformed
                tensor_img = torch.from_numpy(chw_u8).to(self.device, non_blocking=True).float()
            else:
                height, width = image_bgr.shape[:2]
                img_transformed = self.aug.get_transform(image_bgr).apply_image(image_bgr)
                chw_u8 = np.ascontiguousarray(np.transpose(img_transformed, (2, 0, 1)), dtype=np.uint8)
                tensor_img = torch.from_numpy(chw_u8).to(self.device, non_blocking=True).float()
            inputs = [{"image": tensor_img, "height": height, "width": width}]
            if self.use_amp and self.device == "cuda":
                with torch.amp.autocast("cuda", dtype=torch.float16):
                    preds = self.model(inputs)[0]
            else:
                preds = self.model(inputs)[0]

        return self._extract_instances(preds)

    def predict_batch(
        self,
        images_bgr: List[np.ndarray],
        batch_size: int = 10,
        pretransformed: Optional[Sequence[Optional[Tuple[np.ndarray, int, int]]]] = None,
    ) -> List[List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]]:
        """
        Runs batched Stage 2 PointRend inference across multiple candidate leaf crops
        (including cross-voucher crop waves on high-VRAM GPUs).
        Sorts crops by orientation and aspect ratio before batching so Detectron2's
        `ImageList` minimizes zero-padding and VRAM footprint, restoring original crop order.
        Automatically delegates to `self.predict(img)` if `predict` is mocked in unit tests
        or if `len(images_bgr) <= 1`.
        """
        if not images_bgr:
            return []
        has_pre = pretransformed is not None and len(pretransformed) == len(images_bgr)
        if (
            self.model is None
            or "predict" in self.__dict__
            or getattr(self.predict, "__func__", None) is not PointRendInferenceEngine.predict
            or len(images_bgr) <= 1
        ):
            if (
                has_pre
                and len(images_bgr) == 1
                and "predict" not in self.__dict__
                and getattr(self.predict, "__func__", None) is PointRendInferenceEngine.predict
            ):
                return [self.predict(images_bgr[0], pretransformed=pretransformed[0])]  # type: ignore[index]
            return [self.predict(img) for img in images_bgr]

        import torch
        n_imgs = len(images_bgr)
        all_results: List[List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]] = [
            [] for _ in range(n_imgs)
        ]
        eff_bs = max(1, int(batch_size))

        def _aspect_sort_key(idx: int) -> Tuple[int, float]:
            pre_it = pretransformed[idx] if has_pre else None  # type: ignore[index]
            if pre_it is not None:
                _, h_i, w_i = pre_it
            else:
                h_i, w_i = images_bgr[idx].shape[:2]
            is_land = 1 if w_i > h_i else 0
            return (is_land, round(float(h_i) / float(max(w_i, 1)), 2))

        sorted_indices = sorted(range(n_imgs), key=_aspect_sort_key)
        pending_chunks: List[List[int]] = [
            sorted_indices[s : s + eff_bs] for s in range(0, n_imgs, eff_bs)
        ]

        while pending_chunks:
            chunk_idx_list = pending_chunks.pop(0)
            if len(chunk_idx_list) == 1:
                only_i = chunk_idx_list[0]
                pre_only = pretransformed[only_i] if has_pre else None  # type: ignore[index]
                all_results[only_i] = self.predict(images_bgr[only_i], pretransformed=pre_only)
                continue
            try:
                with torch.no_grad():
                    inputs = []
                    for orig_i in chunk_idx_list:
                        img_bgr = images_bgr[orig_i]
                        pre_item = pretransformed[orig_i] if has_pre else None  # type: ignore[index]
                        if pre_item is not None:
                            chw_u8, h, w = pre_item
                        else:
                            h, w = img_bgr.shape[:2]
                            img_t = self.aug.get_transform(img_bgr).apply_image(img_bgr)
                            chw_u8 = np.ascontiguousarray(np.transpose(img_t, (2, 0, 1)), dtype=np.uint8)
                        t_img = torch.from_numpy(chw_u8).to(self.device, non_blocking=True).float()
                        inputs.append({"image": t_img, "height": h, "width": w})
                    if self.use_amp and self.device == "cuda":
                        with torch.amp.autocast("cuda", dtype=torch.float16):
                            preds_list = self.model(inputs)
                    else:
                        preds_list = self.model(inputs)
                for orig_i, preds in zip(chunk_idx_list, preds_list):
                    all_results[orig_i] = self._extract_instances(preds)
                del inputs, preds_list
            except Exception as e:
                logger.debug(f"Batched PointRend halving backoff from chunk={len(chunk_idx_list)} ({e})")
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                mid = max(1, len(chunk_idx_list) // 2)
                pending_chunks = [chunk_idx_list[:mid], chunk_idx_list[mid:]] + pending_chunks
        return all_results


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


def _extract_native_lm2_petiole_anchor(
    leaf_mask: np.ndarray,
    petiole_mask: Optional[np.ndarray],
) -> Optional[Dict[str, Any]]:
    """
    Uses LeafMachine2's native Stage 2 PointRend `petiole` (`class_id == 1`) mask alongside
    the raw `leaf` (`class_id == 0`) mask to locate the anatomical petiole-lamina junction
    `p_base_crop`, distal lamina apex `p_apex_crop`, and native petiole pixel span without
    clipping or modifying `leaf_mask`.
    """
    if (
        petiole_mask is None
        or leaf_mask is None
        or leaf_mask.size == 0
        or petiole_mask.size == 0
        or petiole_mask.shape[:2] != leaf_mask.shape[:2]
    ):
        return None
    if int(np.count_nonzero(petiole_mask)) == 0 or int(np.count_nonzero(leaf_mask)) == 0:
        return None

    bx_l, by_l, bw_l, bh_l = cv2.boundingRect(leaf_mask)
    bx_p, by_p, bw_p, bh_p = cv2.boundingRect(petiole_mask)
    if bw_l == 0 or bh_l == 0 or bw_p == 0 or bh_p == 0:
        return None
    cnts_l, _ = cv2.findContours(
        leaf_mask[by_l : by_l + bh_l, bx_l : bx_l + bw_l],
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_NONE,
        offset=(bx_l, by_l),
    )
    if not cnts_l:
        return None
    leaf_pts = max(cnts_l, key=cv2.contourArea).reshape(-1, 2)
    if len(leaf_pts) == 0:
        return None

    h_full, w_full = leaf_mask.shape[:2]
    ux0 = max(0, min(bx_l, bx_p) - 2)
    uy0 = max(0, min(by_l, by_p) - 2)
    ux1 = min(w_full, max(bx_l + bw_l, bx_p + bw_p) + 2)
    uy1 = min(h_full, max(by_l + bh_l, by_p + bh_p) + 2)
    inv_pet_sub = (petiole_mask[uy0:uy1, ux0:ux1] == 0).astype(np.uint8)
    dt_pet = cv2.distanceTransform(inv_pet_sub, cv2.DIST_L2, 3)
    dists_to_pet = dt_pet[leaf_pts[:, 1] - uy0, leaf_pts[:, 0] - ux0]
    min_idx = int(np.argmin(dists_to_pet))
    min_dist = float(dists_to_pet[min_idx])
    if min_dist > max(28.0, 0.30 * float(max(bw_l, bh_l))):
        return None

    p_base_crop = (int(leaf_pts[min_idx, 0]), int(leaf_pts[min_idx, 1]))
    dists_from_base = np.hypot(leaf_pts[:, 0] - p_base_crop[0], leaf_pts[:, 1] - p_base_crop[1])
    apex_idx = int(np.argmax(dists_from_base))
    p_apex_crop = (int(leaf_pts[apex_idx, 0]), int(leaf_pts[apex_idx, 1]))

    cnts_p, _ = cv2.findContours(
        petiole_mask[by_p : by_p + bh_p, bx_p : bx_p + bw_p],
        cv2.RETR_EXTERNAL,
        cv2.CHAIN_APPROX_NONE,
        offset=(bx_p, by_p),
    )
    if cnts_p:
        pet_pts = np.vstack([c.reshape(-1, 2) for c in cnts_p])
        pet_span_dists = np.hypot(pet_pts[:, 0] - p_base_crop[0], pet_pts[:, 1] - p_base_crop[1])
        petiole_len_px = float(np.max(pet_span_dists)) if len(pet_span_dists) > 0 else 0.0
    else:
        petiole_len_px = 0.0

    return {
        "is_trimmed": False,
        "native_lm2_petiole": True,
        "trim_distance_px": round(petiole_len_px, 2),
        "p_base_crop": p_base_crop,
        "p_apex_crop": p_apex_crop,
    }


def _merge_and_bridge_dissected_crop_mask(
    best_mask: np.ndarray,
    leaf_cands: Sequence[Tuple[float, int, Any, np.ndarray]],
    petiole_mask: Optional[np.ndarray] = None,
    crop_bgr: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Reconstructs a unified dissected leaf silhouette when PointRend predicts lateral pinna
    pairs as either:
      (a) disconnected foreground components within `best_mask`, OR
      (b) secondary `Class 0` (`leaf`) predictions in `leaf_cands[1:]` that are spatially
          disjoint from `best_mask` and anchored to the same native `petiole_mask` (`Class 1`)
          rachis stalk.
    Fast no-op (`0` pixels modified) when `best_mask` is a single connected component and
    no co-predicted pinna lobes share the petiole/rachis.
    """
    if best_mask is None or best_mask.size == 0:
        return best_mask

    merged = best_mask
    best_area = int(np.count_nonzero(best_mask))
    if len(leaf_cands) > 1 and petiole_mask is not None and petiole_mask.shape[:2] == best_mask.shape[:2]:
        if int(np.count_nonzero(petiole_mask)) > 0 and best_area > 0:
            k_pet = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9))
            pet_dil = cv2.dilate((petiole_mask > 0).astype(np.uint8) * 255, k_pet, iterations=1) > 0
            best_bool = best_mask > 0
            if np.any(best_bool & pet_dil):
                merged = best_mask.copy()
                for c_score, c_area, _c_bbox, c_mask in leaf_cands[1:]:
                    if float(c_score) < 0.45 or c_mask is None or c_mask.shape[:2] != best_mask.shape[:2]:
                        continue
                    if c_area < max(16, int(0.025 * best_area)) or c_area > int(0.85 * best_area):
                        continue
                    c_bool = c_mask > 0
                    inter = int(np.count_nonzero(best_bool & c_bool))
                    if float(inter) / float(max(1, min(c_area, best_area))) >= 0.10:
                        continue
                    if np.any(c_bool & pet_dil):
                        merged = np.maximum(merged, (c_bool.astype(np.uint8) * 255))

    return bridge_dissected_rachis_components(merged, petiole_mask=petiole_mask, crop_bgr=crop_bgr)
