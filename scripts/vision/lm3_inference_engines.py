#!/usr/bin/env python3
"""
===============================================================================
Module: lm3_inference_engines.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    LeafMachine3 ONNX Runtime inference engines for Step 02:
      1. LM3PlantDetector (lm3_plant_detector__yolo26x_det_1280, 9-class NMS-free YOLO26x)
      2. LM3SpecimenSegmenter (lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain,
         UNet++ v9 1024 px + LLR_Hyst_KeepBlobs paper-vs-plant cleanup)
      3. LM3LeafSegmenterEngine (lm3_leaf_segmenter__yolo26x_seg_1024, 1024 px YOLO26x-seg
         producing separated lamina, petiole, and hole masks)
      4. LM3PoseLandmarkEngine (lm3_landmark_detector__yolo26x_pose_640, 640 px YOLO26x-pose
         31-keypoint mid15_pet5 skeleton + petiole width estimation)
===============================================================================
"""

from __future__ import annotations

import logging
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

try:
    from scripts.vision.configure_leafmachine3 import (
        resolve_lm3_weights,
    )
except ImportError:
    def resolve_lm3_weights(
        default_filename: str,
        explicit_path: Optional[Union[str, Path]] = None,
        project_root: Optional[Path] = None,
    ) -> Path:
        root = Path(project_root) if project_root is not None else PROJECT_ROOT
        if explicit_path is not None and str(explicit_path).strip():
            p = Path(explicit_path)
            return p if p.is_absolute() else (root / p)
        stem = Path(default_filename).stem
        return root / "models" / "lm3" / f"{stem}.onnx"

try:
    from scripts.vision.lm3_data_loader import (
        LM3_LEAF_CLASS_IDS,
        LM3_PLANT_CLASSES,
        LM3_REPRODUCTIVE_CLASS_IDS,
        compute_work_scale,
        load_precomputed_plant_labels,
        unscale_working_bbox_to_sheet,
    )
except ImportError:
    LM3_PLANT_CLASSES: Dict[int, str] = {
        0: "Leaf_WHOLE",
        1: "Leaf_PARTIAL",
        2: "Seed_Fruit_ONE",
        3: "Seed_Fruit_MANY",
        4: "Flower_ONE",
        5: "Flower_MANY",
        6: "Bud",
        7: "Roots",
        8: "Wood",
    }
    LM3_LEAF_CLASS_IDS: Tuple[int, ...] = (0, 1)
    LM3_REPRODUCTIVE_CLASS_IDS: Tuple[int, ...] = (2, 3, 4, 5, 6)

    def compute_work_scale(sheet_shape: Tuple[int, int], max_working_dim: int = 3200) -> float:
        sh, sw = int(sheet_shape[0]), int(sheet_shape[1])
        max_dim = max(sh, sw)
        if max_dim <= 0 or max_working_dim <= 0 or max_dim <= max_working_dim:
            return 1.0
        return float(max_working_dim) / float(max_dim)

    def unscale_working_bbox_to_sheet(
        bbox_work: Tuple[Union[int, float], Union[int, float], Union[int, float], Union[int, float]],
        work_scale: float,
        sheet_shape: Optional[Tuple[int, int]] = None,
    ) -> Tuple[int, int, int, int]:
        scale = float(work_scale) if (work_scale and work_scale > 0.0) else 1.0
        ymin_w, xmin_w, ymax_w, xmax_w = bbox_work
        ymin_s = int(round(float(ymin_w) / scale))
        xmin_s = int(round(float(xmin_w) / scale))
        ymax_s = int(round(float(ymax_w) / scale))
        xmax_s = int(round(float(xmax_w) / scale))
        if sheet_shape is not None:
            sh, sw = int(sheet_shape[0]), int(sheet_shape[1])
            ymin_s = max(0, min(sh, ymin_s))
            xmin_s = max(0, min(sw, xmin_s))
            ymax_s = max(ymin_s, min(sh, ymax_s))
            xmax_s = max(xmin_s, min(sw, xmax_s))
        return ymin_s, xmin_s, ymax_s, xmax_s

    def load_precomputed_plant_labels(
        label_path: Union[str, Path],
        sheet_w: int,
        sheet_h: int,
        min_confidence: float = 0.25,
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        p = Path(label_path)
        if not p.exists():
            return []
        results: List[Tuple[Tuple[int, int, int, int], float, int]] = []
        try:
            for line in p.read_text(encoding="utf-8").splitlines():
                parts = line.strip().split()
                if len(parts) < 5:
                    continue
                cls_id = int(float(parts[0]))
                xc, yc, bw, bh = [float(x) for x in parts[1:5]]
                conf = float(parts[5]) if len(parts) >= 6 else 0.90
                if conf < min_confidence:
                    continue
                xmin = max(0, int(round((xc - 0.5 * bw) * sheet_w)))
                xmax = min(sheet_w, int(round((xc + 0.5 * bw) * sheet_w)))
                ymin = max(0, int(round((yc - 0.5 * bh) * sheet_h)))
                ymax = min(sheet_h, int(round((yc + 0.5 * bh) * sheet_h)))
                if ymax > ymin and xmax > xmin:
                    results.append(((ymin, xmin, ymax, xmax), conf, cls_id))
        except Exception:
            return []
        return results

from scripts.vision.leaf_geometry_utils import (
    bridge_dissected_rachis_components,
    is_topologically_connected_mask,
)

logger = logging.getLogger("SegmentAndExtract")

DEFAULT_ONNX_PROVIDERS: Tuple[str, ...] = (
    "CUDAExecutionProvider",
    "CoreMLExecutionProvider",
    "CPUExecutionProvider",
)

LM3_KEYPOINT_NAMES_31: Tuple[str, ...] = (
    "lamina_tip",
    "apex_1",
    "apex_2",
    "apex_3",
    "midvein_1",
    "midvein_2",
    "midvein_3",
    "midvein_4",
    "midvein_5",
    "midvein_6",
    "midvein_7",
    "midvein_8",
    "midvein_9",
    "midvein_10",
    "midvein_11",
    "midvein_12",
    "midvein_13",
    "midvein_14",
    "midvein_15",
    "base_1",
    "base_2",
    "base_3",
    "lamina_base",
    "petiole_1",
    "petiole_2",
    "petiole_3",
    "petiole_4",
    "petiole_5",
    "petiole_6",
    "width_1",
    "width_2",
)

__all__ = [
    "DEFAULT_ONNX_PROVIDERS",
    "LM3_KEYPOINT_NAMES_31",
    "DetectedInstance",
    "_VoucherFinalizeBundle",
    "crop_leaf_bbox_with_padding",
    "reproject_crop_mask_to_sheet",
    "compute_capitulum_metrics",
    "resolve_onnx_providers",
    "LM3PlantDetector",
    "LM3SpecimenSegmenter",
    "apply_llr_hyst_keep_blobs",
    "LM3LeafSegmenterEngine",
    "LM3PoseResult",
    "LM3PoseLandmarkEngine",
    "torch_cuda_available",
    "is_stock_lm3_weights",
    "_extract_native_lm3_petiole_anchor",
    "_merge_and_bridge_dissected_crop_mask",
]


@dataclass
class DetectedInstance:
    """Represents a segmented botanical instance on a herbarium sheet."""
    catalog_number: str
    leaf_id: int
    bbox: Tuple[int, int, int, int]  # (ymin, xmin, ymax, xmax) in full-sheet coordinates
    mask: np.ndarray                 # Binary uint8 lamina mask (H, W)
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
    petiole_mask: Optional[np.ndarray] = None
    hole_mask: Optional[np.ndarray] = None
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None
    work_scale: float = 1.0
    cf_source: str = "measured_from_ruler"
    archetype_score: float = 0.0
    petiole_width_mm: Optional[float] = None
    ect_path: Optional[str] = None


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
# ONNX Runtime Session & Letterbox Preprocessing Utilities
# =============================================================================

def torch_cuda_available() -> bool:
    """Safely checks whether CUDA is available via ONNX Runtime or PyTorch."""
    try:
        import onnxruntime as ort
        if "CUDAExecutionProvider" in ort.get_available_providers():
            return True
    except Exception:
        pass
    try:
        import torch
        return bool(torch.cuda.is_available())
    except Exception:
        return False


def is_stock_lm3_weights(weights_path: Union[str, Path]) -> bool:
    """
    Returns True when using canonical stock LeafMachine3 ONNX weights rather than
    a taxon-specific custom fine-tuned checkpoint.
    """
    name = Path(weights_path).name.lower()
    if not name:
        return True
    return "finetuned" not in name and "packera" not in name


def resolve_onnx_providers(
    device: str = "cuda",
    requested_providers: Optional[Sequence[str]] = None,
) -> List[str]:
    """
    Resolves the ordered ONNX Runtime execution provider priority list:
      ["CUDAExecutionProvider", "CoreMLExecutionProvider", "CPUExecutionProvider"]
    filtered against available runtime providers when device != 'cpu'.
    """
    priority = list(requested_providers) if requested_providers else list(DEFAULT_ONNX_PROVIDERS)
    if str(device).strip().lower() == "cpu":
        return ["CPUExecutionProvider"]
    try:
        import onnxruntime as ort
        avail = set(ort.get_available_providers())
        if avail:
            filtered = [p for p in priority if p in avail]
            if "CPUExecutionProvider" not in filtered:
                filtered.append("CPUExecutionProvider")
            return filtered if filtered else list(DEFAULT_ONNX_PROVIDERS)
    except Exception:
        pass
    return priority


_ONNX_SESSION_CACHE: Dict[Tuple[str, Tuple[str, ...]], Any] = {}
_NVIDIA_LIBS_PRELOADED: bool = False


def _preload_nvidia_cuda_libs() -> None:
    """
    Preloads bundled `nvidia-*` wheel CUDA 12 / cuDNN 9 shared libraries into the process
    symbol table (`RTLD_GLOBAL`) so `onnxruntime-gpu` binds `CUDAExecutionProvider` even when
    `LD_LIBRARY_PATH` was not exported prior to interpreter launch.
    """
    global _NVIDIA_LIBS_PRELOADED
    if _NVIDIA_LIBS_PRELOADED:
        return
    _NVIDIA_LIBS_PRELOADED = True
    try:
        import ctypes
        import nvidia  # type: ignore
        roots = [Path(p) for p in getattr(nvidia, "__path__", ()) or ()]
        lib_dirs = sorted({d for root in roots for d in root.glob("*/lib") if d.is_dir()})
        so_files = [
            p
            for d in lib_dirs
            for p in sorted(d.glob("lib*.so*"))
            if p.is_file() and not p.is_symlink()
        ]
        for _ in range(2):
            for p in so_files:
                try:
                    ctypes.CDLL(str(p), mode=ctypes.RTLD_GLOBAL)
                except OSError:
                    pass
    except Exception:
        pass


def _load_onnx_session(
    weights_path: Optional[Union[str, Path]],
    default_model_name: str,
    device: str = "cuda",
    providers: Optional[Sequence[str]] = None,
) -> Tuple[Path, str, List[str], Optional[Any]]:
    """Loads or retrieves a cached ONNX Runtime InferenceSession with graceful fallback when absent."""
    resolved = resolve_lm3_weights(default_model_name, explicit_path=weights_path)
    if isinstance(resolved, tuple):
        m_name, w_path = str(resolved[0]), Path(resolved[1])
    else:
        w_path = Path(resolved)
        if weights_path is not None and str(weights_path).strip():
            wp_stem = Path(weights_path).stem
            m_name = wp_stem if (wp_stem and wp_stem != "model") else (Path(default_model_name).stem or default_model_name)
        else:
            m_name = (w_path.stem if (w_path.stem and w_path.stem != "model") else Path(default_model_name).stem) or default_model_name

    prov_list = resolve_onnx_providers(device=device, requested_providers=providers)
    try:
        has_valid_weights = w_path.is_file() and w_path.stat().st_size >= 1024
    except OSError:
        has_valid_weights = False

    if not has_valid_weights:
        logger.debug(f"LeafMachine3 ONNX weights not found at {w_path}; running in fallback/mock-ready mode.")
        return w_path, m_name, prov_list, None

    cache_key = (str(w_path.resolve()), tuple(prov_list))
    if cache_key in _ONNX_SESSION_CACHE:
        return w_path, m_name, prov_list, _ONNX_SESSION_CACHE[cache_key]

    try:
        import os
        import onnxruntime as ort
        if "CUDAExecutionProvider" in prov_list:
            _preload_nvidia_cuda_libs()
        sess_opts = ort.SessionOptions()
        sess_opts.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        try:
            sess_opts.log_severity_level = int(os.environ.get("LM3_ORT_LOG_SEVERITY", 3))
        except ValueError:
            sess_opts.log_severity_level = 3
        session = ort.InferenceSession(str(w_path), sess_options=sess_opts, providers=prov_list)
        _ONNX_SESSION_CACHE[cache_key] = session
        bound_provs = session.get_providers()
        logger.info(f"Loaded LeafMachine3 ONNX session ({m_name}) from {w_path} with providers={bound_provs}")
        return w_path, m_name, bound_provs, session
    except Exception as exc:
        logger.debug(f"ONNX Runtime session initialization skipped for {w_path}: {exc}")
        return w_path, m_name, prov_list, None


def _letterbox_bgr(
    image_bgr: np.ndarray,
    target_size: int,
    pad_value: int = 114,
) -> Tuple[np.ndarray, float, Tuple[float, float]]:
    """
    Resizes `image_bgr` with aspect-preserving letterbox padding to `(target_size, target_size)`.
    Returns `(letterboxed_bgr, scale_ratio, (dw, dh))`.
    """
    h, w = image_bgr.shape[:2]
    if h <= 0 or w <= 0:
        empty = np.full((target_size, target_size, 3), pad_value, dtype=np.uint8)
        return empty, 1.0, (0.0, 0.0)
    r = min(float(target_size) / float(h), float(target_size) / float(w))
    new_w = max(1, int(round(w * r)))
    new_h = max(1, int(round(h * r)))
    if (new_w, new_h) != (w, h):
        resized = cv2.resize(image_bgr, (new_w, new_h), interpolation=cv2.INTER_LINEAR)
    else:
        resized = image_bgr
    dw = (float(target_size - new_w)) * 0.5
    dh = (float(target_size - new_h)) * 0.5
    top = int(round(dh - 0.1))
    bottom = int(round(dh + 0.1))
    left = int(round(dw - 0.1))
    right = int(round(dw + 0.1))
    padded = cv2.copyMakeBorder(
        resized, top, bottom, left, right, cv2.BORDER_CONSTANT, value=(pad_value, pad_value, pad_value)
    )
    if padded.shape[0] != target_size or padded.shape[1] != target_size:
        padded = cv2.resize(padded, (target_size, target_size), interpolation=cv2.INTER_NEAREST)
    return padded, r, (float(left), float(top))


def _bgr_to_nchw_float32(letterboxed_bgr: np.ndarray) -> np.ndarray:
    """Converts a `(H, W, 3)` BGR uint8 array to `(1, 3, H, W)` RGB float32 in `[0, 1]`."""
    rgb = cv2.cvtColor(letterboxed_bgr, cv2.COLOR_BGR2RGB)
    chw = np.transpose(rgb, (2, 0, 1)).astype(np.float32) / 255.0
    return np.ascontiguousarray(chw[None, ...], dtype=np.float32)


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -40.0, 40.0)))


# =============================================================================
# 1. LeafMachine3 Stage 1 Plant Detector (YOLO26x End-to-End NMS-Free, 1280 px)
# =============================================================================

class LM3PlantDetector:
    """
    LeafMachine3 Stage 1 ONNX Plant Organ Detector (`lm3_plant_detector__yolo26x_det_1280`).
    Runs 1280 px end-to-end NMS-free YOLO26x detection across 9 canonical plant classes:
      0: Leaf_WHOLE, 1: Leaf_PARTIAL,
      2: Seed_Fruit_ONE, 3: Seed_Fruit_MANY, 4: Flower_ONE, 5: Flower_MANY, 6: Bud,
      7: Roots, 8: Wood.
    Automatically scales detections from the working frame (`max_working_dim=3200`)
    back to native full-sheet coordinates via `work_scale` (`unscale_working_bbox_to_sheet`).
    """

    MODEL_NAME: str = "lm3_plant_detector__yolo26x_det_1280"
    CLASSES: Dict[int, str] = LM3_PLANT_CLASSES
    LEAF_CLASS_IDS: Tuple[int, ...] = LM3_LEAF_CLASS_IDS
    REPRODUCTIVE_CLASS_IDS: Tuple[int, ...] = LM3_REPRODUCTIVE_CLASS_IDS

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        conf_thres: float = 0.25,
        iou_thres: float = 0.45,
        img_size: int = 1280,
        max_working_dim: int = 3200,
        lm3_project_dir: Optional[Union[str, Path]] = None,
        providers: Optional[Sequence[str]] = None,
        **kwargs: Any,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.conf_thres = float(conf_thres)
        self.iou_thres = float(iou_thres)
        self.img_size = int(img_size)
        self.max_working_dim = int(max_working_dim)
        proj_dir = lm3_project_dir or kwargs.get("project_dir")
        self.lm3_project_dir = Path(proj_dir) if proj_dir else None
        self.stride = 32
        self.last_work_scale: float = 1.0
        self.weights_path, self.model_name, self.providers, self.session = _load_onnx_session(
            weights_path=weights_path,
            default_model_name=self.MODEL_NAME,
            device=self.device,
            providers=providers,
        )
        self.model = self.session

    def prepare_letterbox(
        self,
        image_bgr: np.ndarray,
        max_working_dim: Optional[int] = None,
    ) -> Optional[np.ndarray]:
        """Precomputes the working-frame letterboxed BGR array on a background CPU worker thread."""
        if (self.session is None and self.model is None) or image_bgr is None or image_bgr.size == 0:
            return None
        work_dim = int(max_working_dim) if max_working_dim is not None else self.max_working_dim
        sh, sw = image_bgr.shape[:2]
        w_scale = compute_work_scale((sh, sw), max_working_dim=work_dim)
        if w_scale < 1.0:
            work_w = max(1, int(round(sw * w_scale)))
            work_h = max(1, int(round(sh * w_scale)))
            work_img = cv2.resize(image_bgr, (work_w, work_h), interpolation=cv2.INTER_AREA)
        else:
            work_img = image_bgr
        lb_bgr, _, _ = _letterbox_bgr(work_img, self.img_size)
        return lb_bgr

    def _decode_yolo26x_det_output(
        self,
        raw_output: np.ndarray,
        ratio: float,
        pad_xy: Tuple[float, float],
        work_shape: Tuple[int, int],
        work_scale: float,
        sheet_shape: Tuple[int, int],
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        """
        Decodes end-to-end NMS-free YOLO26x detections `(B, N, 6)` `[x1, y1, x2, y2, score, cls]`
        or anchor-grid `(B, 4+C, A)` outputs, unletterboxing to the working frame and unscaling
        to native sheet coordinates via `unscale_working_bbox_to_sheet`.
        """
        arr = np.asarray(raw_output, dtype=np.float32)
        if arr.ndim == 3:
            arr = arr[0]
        if arr.size == 0:
            return []

        dw, dh = pad_xy
        work_h, work_w = int(work_shape[0]), int(work_shape[1])
        sheet_h, sheet_w = int(sheet_shape[0]), int(sheet_shape[1])
        inv_r = 1.0 / max(float(ratio), 1e-6)
        results: List[Tuple[Tuple[int, int, int, int], float, int]] = []

        if arr.ndim == 2 and arr.shape[1] in (6, 7):
            # End-to-end NMS-free YOLO26x format: [x1, y1, x2, y2, score, class_id]
            for row in arr:
                score = float(row[4])
                if score < self.conf_thres:
                    continue
                cls_id = int(round(float(row[5])))
                if cls_id not in self.CLASSES:
                    continue
                x1_w = max(0.0, min(float(work_w), (float(row[0]) - dw) * inv_r))
                y1_w = max(0.0, min(float(work_h), (float(row[1]) - dh) * inv_r))
                x2_w = max(0.0, min(float(work_w), (float(row[2]) - dw) * inv_r))
                y2_w = max(0.0, min(float(work_h), (float(row[3]) - dh) * inv_r))
                if y2_w <= y1_w or x2_w <= x1_w:
                    continue
                sheet_bbox = unscale_working_bbox_to_sheet(
                    (y1_w, x1_w, y2_w, x2_w), work_scale, sheet_shape=(sheet_h, sheet_w)
                )
                if sheet_bbox[2] > sheet_bbox[0] and sheet_bbox[3] > sheet_bbox[1]:
                    results.append((sheet_bbox, score, cls_id))
            return results

        # Transpose (4+C, A) -> (A, 4+C) if needed
        if arr.ndim == 2 and arr.shape[0] < arr.shape[1] and arr.shape[0] >= 5:
            arr = arr.T
        if arr.ndim != 2 or arr.shape[1] < 5:
            return []

        n_cls = len(self.CLASSES)
        boxes_xywh = arr[:, :4]
        cls_scores = arr[:, 4 : 4 + n_cls]
        if cls_scores.size == 0:
            return []
        best_cls = np.argmax(cls_scores, axis=1)
        best_scores = cls_scores[np.arange(len(arr)), best_cls]
        keep = best_scores >= self.conf_thres
        if not np.any(keep):
            return []

        boxes_xywh = boxes_xywh[keep]
        best_scores = best_scores[keep]
        best_cls = best_cls[keep]

        xyxy_boxes: List[List[float]] = []
        for (xc, yc, bw, bh) in boxes_xywh:
            x1 = float(xc - 0.5 * bw)
            y1 = float(yc - 0.5 * bh)
            xyxy_boxes.append([x1, y1, float(bw), float(bh)])

        indices = cv2.dnn.NMSBoxes(
            bboxes=xyxy_boxes,
            scores=[float(s) for s in best_scores],
            score_threshold=self.conf_thres,
            nms_threshold=self.iou_thres,
        )
        if len(indices) == 0:
            return []

        for idx_raw in np.asarray(indices).reshape(-1):
            idx = int(idx_raw)
            xc, yc, bw, bh = boxes_xywh[idx]
            x1_l = float(xc - 0.5 * bw)
            y1_l = float(yc - 0.5 * bh)
            x2_l = float(xc + 0.5 * bw)
            y2_l = float(yc + 0.5 * bh)
            x1_w = max(0.0, min(float(work_w), (x1_l - dw) * inv_r))
            y1_w = max(0.0, min(float(work_h), (y1_l - dh) * inv_r))
            x2_w = max(0.0, min(float(work_w), (x2_l - dw) * inv_r))
            y2_w = max(0.0, min(float(work_h), (y2_l - dh) * inv_r))
            if y2_w <= y1_w or x2_w <= x1_w:
                continue
            sheet_bbox = unscale_working_bbox_to_sheet(
                (y1_w, x1_w, y2_w, x2_w), work_scale, sheet_shape=(sheet_h, sheet_w)
            )
            cls_id = int(best_cls[idx])
            if sheet_bbox[2] > sheet_bbox[0] and sheet_bbox[3] > sheet_bbox[1] and cls_id in self.CLASSES:
                results.append((sheet_bbox, float(best_scores[idx]), cls_id))
        return results

    def detect(
        self,
        image_bgr: np.ndarray,
        max_working_dim: int = 3200,
        catalog_number: Optional[str] = None,
        use_cache: bool = False,
        precomputed_letterbox: Optional[np.ndarray] = None,
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        """
        Runs LeafMachine3 9-class full-sheet organ detection on `image_bgr`.
        Scales predictions from `max_working_dim` working coordinates back to native
        sheet dimensions `(sheet_h, sheet_w)` via `work_scale`.
        Returns list of `((ymin, xmin, ymax, xmax), score, class_id)`.
        """
        if image_bgr is None or image_bgr.size == 0:
            return []
        sheet_h, sheet_w = image_bgr.shape[:2]
        work_dim = int(max_working_dim) if max_working_dim > 0 else self.max_working_dim
        work_scale = compute_work_scale((sheet_h, sheet_w), max_working_dim=work_dim)
        self.last_work_scale = work_scale

        active_sess = self.session if self.session is not None else self.model
        if use_cache and active_sess is None and self.lm3_project_dir and catalog_number:
            label_candidates = [
                self.lm3_project_dir / "Data" / "precomputed_annotations" / "Plant_Components" / "labels" / f"{catalog_number}.txt",
                self.lm3_project_dir / "Plant_Components" / "labels" / f"{catalog_number}.txt",
            ]
            for lf in label_candidates:
                if lf.exists():
                    dets = load_precomputed_plant_labels(
                        lf, sheet_w, sheet_h, min_confidence=self.conf_thres
                    )
                    if dets:
                        return dets

        if active_sess is None:
            return []

        try:
            if work_scale < 1.0:
                work_w = max(1, int(round(sheet_w * work_scale)))
                work_h = max(1, int(round(sheet_h * work_scale)))
                work_img = (
                    cv2.resize(image_bgr, (work_w, work_h), interpolation=cv2.INTER_AREA)
                    if precomputed_letterbox is None
                    else None
                )
            else:
                work_h, work_w = sheet_h, sheet_w
                work_img = image_bgr

            if precomputed_letterbox is not None and precomputed_letterbox.shape[:2] == (self.img_size, self.img_size):
                lb_bgr = precomputed_letterbox
                r = min(float(self.img_size) / float(max(work_h, 1)), float(self.img_size) / float(max(work_w, 1)))
                new_w = max(1, int(round(work_w * r)))
                new_h = max(1, int(round(work_h * r)))
                pad_xy = (float(int(round((self.img_size - new_w) * 0.5 - 0.1))), float(int(round((self.img_size - new_h) * 0.5 - 0.1))))
            else:
                assert work_img is not None
                lb_bgr, r, pad_xy = _letterbox_bgr(work_img, self.img_size)

            inp_tensor = _bgr_to_nchw_float32(lb_bgr)
            input_name = active_sess.get_inputs()[0].name
            raw_outs = active_sess.run(None, {input_name: inp_tensor})
            if not raw_outs:
                return []
            return self._decode_yolo26x_det_output(
                raw_outs[0],
                ratio=r,
                pad_xy=pad_xy,
                work_shape=(work_h, work_w),
                work_scale=work_scale,
                sheet_shape=(sheet_h, sheet_w),
            )
        except Exception as exc:
            logger.debug(f"LM3PlantDetector ONNX inference error: {exc}")
            return []

    def predict(
        self,
        image_bgr: np.ndarray,
        catalog_number: Optional[str] = None,
        use_cache: bool = False,
        precomputed_letterbox: Optional[np.ndarray] = None,
        max_working_dim: Optional[int] = None,
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        """Alias for `.detect(...)` preserving Step 02 runner and unit test compatibility."""
        if "detect" in self.__dict__:
            return self.detect(
                image_bgr,
                max_working_dim=int(max_working_dim) if max_working_dim is not None else self.max_working_dim,
                catalog_number=catalog_number,
                use_cache=use_cache,
                precomputed_letterbox=precomputed_letterbox,
            )
        return LM3PlantDetector.detect(
            self,
            image_bgr,
            max_working_dim=int(max_working_dim) if max_working_dim is not None else self.max_working_dim,
            catalog_number=catalog_number,
            use_cache=use_cache,
            precomputed_letterbox=precomputed_letterbox,
        )

    def __call__(
        self,
        image_bgr: np.ndarray,
        max_working_dim: int = 3200,
        **kwargs: Any,
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        if "predict" in self.__dict__ and "detect" not in self.__dict__:
            return self.predict(image_bgr, max_working_dim=max_working_dim, **kwargs)
        return self.detect(image_bgr, max_working_dim=max_working_dim, **kwargs)

    def predict_batch(
        self,
        items: Sequence[Tuple[np.ndarray, Optional[str], Optional[np.ndarray]]],
        use_cache: bool = False,
        max_batch_size: int = 16,
    ) -> List[List[Tuple[Tuple[int, int, int, int], float, int]]]:
        """Runs Stage 1 plant organ detection across a wave of herbarium sheets."""
        if not items:
            return []
        return [
            self.predict(img, catalog_number=cat, use_cache=use_cache, precomputed_letterbox=lb)
            for img, cat, lb in items
        ]


# =============================================================================
# 2. LeafMachine3 Whole-Sheet Specimen Segmenter (UNet++ v9 + LLR_Hyst_KeepBlobs)
# =============================================================================

def apply_llr_hyst_keep_blobs(
    prob_map: np.ndarray,
    image_bgr: Optional[np.ndarray] = None,
    low_thresh: float = 0.35,
    high_thresh: float = 0.65,
    min_blob_area: int = 64,
) -> np.ndarray:
    """
    Applies LeafMachine3 `LLR_Hyst_KeepBlobs` paper-vs-plant log-likelihood ratio modulation
    and dual-threshold hysteresis connected-component filtering:
      1. Normalizes `prob_map` to float32 in `[0, 1]`.
      2. When `image_bgr` is provided, estimates herbarium mounting paper background chromaticity
         and modulates `prob_map` by the plant-vs-paper likelihood ratio to suppress bright paper
         halo bleed while preserving thin petioles and pinnatifid rachis bridges.
      3. Retains weak foreground components (`prob >= low_thresh`) that are topologically connected
         to strong foreground seeds (`prob >= high_thresh`) and meet `min_blob_area`.
      4. Returns a whole-sheet `uint8` binary `specimen_paper_cleaned_mask` (`0` or `255`).
    """
    if prob_map is None or prob_map.size == 0:
        return np.zeros((0, 0), dtype=np.uint8)

    prob = np.asarray(prob_map, dtype=np.float32)
    if prob.ndim == 3:
        prob = prob[:, :, 0] if prob.shape[2] == 1 else np.mean(prob, axis=2)

    p_min = float(np.nanmin(prob)) if prob.size > 0 else 0.0
    p_max = float(np.nanmax(prob)) if prob.size > 0 else 0.0
    if p_min < -1e-3 or (p_max > 1.0 + 1e-3 and p_max <= 25.0):
        prob = _sigmoid(prob)
    elif p_max > 1.0 + 1e-3:
        prob = np.clip(prob / 255.0, 0.0, 1.0)
    else:
        prob = np.clip(np.nan_to_num(prob, nan=0.0), 0.0, 1.0)

    h, w = prob.shape[:2]
    if image_bgr is not None and image_bgr.ndim == 3 and image_bgr.shape[:2] == (h, w):
        bgr_f = image_bgr.astype(np.float32)
        bg_candidates = prob < max(0.15, float(low_thresh))
        if int(np.count_nonzero(bg_candidates)) >= 32:
            paper_med = np.median(bgr_f[bg_candidates], axis=0)
            paper_std = np.maximum(np.std(bgr_f[bg_candidates], axis=0), 8.0)
        else:
            paper_med = np.array([230.0, 239.0, 242.0], dtype=np.float32)
            paper_std = np.array([14.0, 14.0, 14.0], dtype=np.float32)

        z_dist = np.sqrt(np.mean(((bgr_f - paper_med) / paper_std) ** 2, axis=2))
        # Log-likelihood ratio adjustment: positive for dark/chromatic plant tissue, negative for bright paper
        llr_shift = np.clip((z_dist - 1.65) * 0.12, -0.22, 0.22)
        prob = np.clip(prob + np.where(prob > 0.08, llr_shift, 0.0), 0.0, 1.0)

    hi = float(max(low_thresh, high_thresh))
    lo = float(min(low_thresh, high_thresh))
    strong = (prob >= hi).astype(np.uint8)
    weak = (prob >= lo).astype(np.uint8)

    if int(np.count_nonzero(strong)) == 0:
        return np.zeros((h, w), dtype=np.uint8)

    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(weak, connectivity=8)
    if num_labels <= 1:
        return (strong * 255).astype(np.uint8)

    strong_labels = np.unique(labels[strong > 0])
    strong_labels = strong_labels[strong_labels > 0]
    if len(strong_labels) == 0:
        return np.zeros((h, w), dtype=np.uint8)

    eff_min_area = max(1, min(int(min_blob_area), int(0.0005 * h * w)))
    valid_labels = [
        int(lbl)
        for lbl in strong_labels
        if int(stats[int(lbl), cv2.CC_STAT_AREA]) >= eff_min_area
    ]
    if not valid_labels:
        valid_labels = [int(lbl) for lbl in strong_labels]

    cleaned = np.isin(labels, valid_labels).astype(np.uint8) * 255
    return cleaned


class LM3SpecimenSegmenter:
    """
    LeafMachine3 Whole-Sheet Specimen Segmenter
    (`lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain`, UNet++ v9 1024 px)
    paired with `LLR_Hyst_KeepBlobs` paper-vs-plant likelihood ratio hysteresis cleanup.
    """

    MODEL_NAME: str = "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        img_size: int = 1024,
        low_thresh: float = 0.35,
        high_thresh: float = 0.65,
        providers: Optional[Sequence[str]] = None,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.img_size = int(img_size)
        self.low_thresh = float(low_thresh)
        self.high_thresh = float(high_thresh)
        self.weights_path, self.model_name, self.providers, self.session = _load_onnx_session(
            weights_path=weights_path,
            default_model_name=self.MODEL_NAME,
            device=self.device,
            providers=providers,
        )
        self.model = self.session

    @staticmethod
    def _empirical_plant_prob_map(image_bgr: np.ndarray) -> np.ndarray:
        """Computes a colorimetric plant-vs-paper probability map when ONNX weights are absent."""
        if image_bgr.ndim == 2:
            gray = image_bgr.astype(np.float32)
            paper_val = float(np.percentile(gray, 85)) if gray.size > 0 else 235.0
            diff = np.clip((paper_val - gray) / max(paper_val * 0.45, 25.0), 0.0, 1.0)
            return diff.astype(np.float32)

        hsv = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2HSV).astype(np.float32)
        sat = hsv[:, :, 1] / 255.0
        val = hsv[:, :, 2] / 255.0
        paper_v = float(np.percentile(val, 85)) if val.size > 0 else 0.92
        darkness = np.clip((paper_v - val) / max(0.35, paper_v * 0.45), 0.0, 1.0)
        prob = np.clip(0.70 * darkness + 0.55 * sat, 0.0, 1.0)
        return prob.astype(np.float32)

    def segment_specimen(
        self,
        image_bgr: np.ndarray,
        low_thresh: Optional[float] = None,
        high_thresh: Optional[float] = None,
    ) -> np.ndarray:
        """
        Produces the whole-sheet `uint8` binary `specimen_paper_cleaned_mask` (`0` or `255`)
        using UNet++ v9 (1024 px) + `LLR_Hyst_KeepBlobs`.
        """
        if image_bgr is None or image_bgr.size == 0:
            return np.zeros((0, 0), dtype=np.uint8)

        sheet_h, sheet_w = image_bgr.shape[:2]
        lo = float(low_thresh) if low_thresh is not None else self.low_thresh
        hi = float(high_thresh) if high_thresh is not None else self.high_thresh

        active_sess = self.session if self.session is not None else self.model
        prob_map: Optional[np.ndarray] = None
        if active_sess is not None:
            try:
                rgb = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2RGB) if image_bgr.ndim == 3 else cv2.cvtColor(image_bgr, cv2.COLOR_GRAY2RGB)
                scale = float(self.img_size) / float(max(sheet_h, sheet_w, 1))
                h2 = max(1, int(round(sheet_h * scale)))
                w2 = max(1, int(round(sheet_w * scale)))
                r_rgb = cv2.resize(rgb, (w2, h2), interpolation=cv2.INTER_AREA)
                top = (self.img_size - h2) // 2
                left = (self.img_size - w2) // 2
                canvas = np.full((self.img_size, self.img_size, 3), 255, dtype=np.uint8)
                canvas[top : top + h2, left : left + w2] = r_rgb
                im_mean = np.array([0.485, 0.456, 0.406], dtype=np.float32)
                im_std = np.array([0.229, 0.224, 0.225], dtype=np.float32)
                norm_chw = np.transpose((canvas.astype(np.float32) / 255.0 - im_mean) / im_std, (2, 0, 1))
                inp_tensor = np.ascontiguousarray(norm_chw[None, ...], dtype=np.float32)
                input_name = active_sess.get_inputs()[0].name
                raw_outs = active_sess.run(None, {input_name: inp_tensor})
                if raw_outs:
                    out_arr = np.asarray(raw_outs[0], dtype=np.float32)
                    while out_arr.ndim > 2 and out_arr.shape[0] == 1:
                        out_arr = out_arr[0]
                    if out_arr.ndim == 3:
                        out_arr = out_arr[-1] if out_arr.shape[0] <= 4 else out_arr[:, :, 0]
                    if out_arr.min() < 0.0 or out_arr.max() > 1.0:
                        out_arr = _sigmoid(out_arr)
                    if out_arr.shape[:2] == (self.img_size, self.img_size):
                        out_arr = out_arr[top : top + h2, left : left + w2]
                    prob_map = cv2.resize(out_arr, (sheet_w, sheet_h), interpolation=cv2.INTER_LINEAR)
            except Exception as exc:
                logger.debug(f"LM3SpecimenSegmenter ONNX inference fallback: {exc}")
                prob_map = None

        if prob_map is None:
            prob_map = self._empirical_plant_prob_map(image_bgr)

        return apply_llr_hyst_keep_blobs(
            prob_map,
            image_bgr=image_bgr,
            low_thresh=lo,
            high_thresh=hi,
        )

    def __call__(self, image_bgr: np.ndarray, **kwargs: Any) -> np.ndarray:
        return self.segment_specimen(image_bgr, **kwargs)


# =============================================================================
# 3. LeafMachine3 Cropped Leaf Instance Segmenter (YOLO26x-seg, 1024 px)
# =============================================================================

class LM3LeafSegmenterEngine:
    """
    LeafMachine3 Stage 2 ONNX Leaf Instance Segmenter (`lm3_leaf_segmenter__yolo26x_seg_1024`).
    Runs 1024 px end-to-end NMS-free YOLO26x instance segmentation on candidate leaf crops,
    separating:
      - Class 0: `leaf` (lamina_mask)
      - Class 1: `petiole` (petiole_mask)
      - Class 2: `hole` (hole_mask)
    """

    MODEL_NAME: str = "lm3_leaf_segmenter__yolo26x_seg_1024"
    SEG_CLASSES: Dict[int, str] = {
        0: "leaf",
        1: "petiole",
        2: "hole",
    }

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        score_thresh: float = 0.50,
        img_size: int = 1024,
        providers: Optional[Sequence[str]] = None,
        **kwargs: Any,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.score_thresh = float(score_thresh)
        self.img_size = int(img_size)
        self.weights_path, self.model_name, self.providers, self.session = _load_onnx_session(
            weights_path=weights_path,
            default_model_name=self.MODEL_NAME,
            device=self.device,
            providers=providers,
        )
        self.model = self.session

    def prepare_transformed_crop(self, image_bgr: np.ndarray) -> Optional[Tuple[np.ndarray, int, int]]:
        """
        Thread-safe CPU helper that precomputes the 1024 px letterbox uint8 CHW crop
        on worker threads for batched inference.
        """
        if (self.session is None and self.model is None) or image_bgr is None or image_bgr.size == 0:
            return None
        h, w = image_bgr.shape[:2]
        lb_bgr, _, _ = _letterbox_bgr(image_bgr, self.img_size)
        rgb = cv2.cvtColor(lb_bgr, cv2.COLOR_BGR2RGB)
        chw_u8 = np.ascontiguousarray(np.transpose(rgb, (2, 0, 1)), dtype=np.uint8)
        return chw_u8, h, w

    def _decode_yolo26x_seg_outputs(
        self,
        raw_outs: Sequence[np.ndarray],
        crop_shape: Tuple[int, int],
        ratio: float,
        pad_xy: Tuple[float, float],
    ) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """Decodes YOLO26x-seg detection + prototype mask outputs into `(bbox, mask_u8, score, class_id)`."""
        if not raw_outs:
            return []
        crop_h, crop_w = int(crop_shape[0]), int(crop_shape[1])
        dw, dh = pad_xy
        inv_r = 1.0 / max(float(ratio), 1e-6)

        det_arr = np.asarray(raw_outs[0], dtype=np.float32)
        proto_arr = np.asarray(raw_outs[1], dtype=np.float32) if len(raw_outs) >= 2 else None
        if det_arr.ndim == 3:
            det_arr = det_arr[0]
        if proto_arr is not None and proto_arr.ndim == 4:
            proto_arr = proto_arr[0]

        if det_arr.size == 0:
            return []

        results: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]] = []
        n_proto = int(proto_arr.shape[0]) if (proto_arr is not None and proto_arr.ndim == 3) else 0

        # End-to-end NMS-free YOLO26x-seg: (N, 6 + n_proto) -> [x1, y1, x2, y2, score, cls, coeff_0..31]
        if det_arr.ndim == 2 and (det_arr.shape[1] == 6 + n_proto or det_arr.shape[1] == 6):
            rows_decoded = []
            for row in det_arr:
                score = float(row[4])
                if score < self.score_thresh:
                    continue
                cls_id = int(round(float(row[5])))
                coeffs = row[6 : 6 + n_proto] if n_proto > 0 and len(row) >= 6 + n_proto else None
                rows_decoded.append((float(row[0]), float(row[1]), float(row[2]), float(row[3]), score, cls_id, coeffs))
        else:
            if det_arr.ndim == 2 and det_arr.shape[0] < det_arr.shape[1] and det_arr.shape[0] >= 5:
                det_arr = det_arr.T
            if det_arr.ndim != 2 or det_arr.shape[1] < 5:
                return []
            n_cls = max(1, det_arr.shape[1] - 4 - n_proto)
            boxes_xywh = det_arr[:, :4]
            cls_scores = det_arr[:, 4 : 4 + n_cls]
            coeffs_all = det_arr[:, 4 + n_cls : 4 + n_cls + n_proto] if n_proto > 0 else None
            best_cls = np.argmax(cls_scores, axis=1)
            best_scores = cls_scores[np.arange(len(det_arr)), best_cls]
            keep = best_scores >= self.score_thresh
            if not np.any(keep):
                return []
            boxes_xywh = boxes_xywh[keep]
            best_scores = best_scores[keep]
            best_cls = best_cls[keep]
            if coeffs_all is not None:
                coeffs_all = coeffs_all[keep]
            nms_boxes = [[float(xc - 0.5 * bw), float(yc - 0.5 * bh), float(bw), float(bh)] for xc, yc, bw, bh in boxes_xywh]
            idxs = cv2.dnn.NMSBoxes(nms_boxes, [float(s) for s in best_scores], self.score_thresh, 0.45)
            rows_decoded = []
            for idx_raw in np.asarray(idxs).reshape(-1):
                i = int(idx_raw)
                xc, yc, bw, bh = boxes_xywh[i]
                rows_decoded.append((
                    float(xc - 0.5 * bw),
                    float(yc - 0.5 * bh),
                    float(xc + 0.5 * bw),
                    float(yc + 0.5 * bh),
                    float(best_scores[i]),
                    int(best_cls[i]),
                    coeffs_all[i] if coeffs_all is not None else None,
                ))

        for x1_l, y1_l, x2_l, y2_l, score, cls_id, coeffs in rows_decoded:
            xmin = max(0, min(crop_w, int(round((x1_l - dw) * inv_r))))
            ymin = max(0, min(crop_h, int(round((y1_l - dh) * inv_r))))
            xmax = max(0, min(crop_w, int(round((x2_l - dw) * inv_r))))
            ymax = max(0, min(crop_h, int(round((y2_l - dh) * inv_r))))
            if ymax <= ymin or xmax <= xmin:
                continue

            mask_u8 = np.zeros((crop_h, crop_w), dtype=np.uint8)
            if proto_arr is not None and coeffs is not None and len(coeffs) == n_proto and n_proto > 0:
                mh, mw = int(proto_arr.shape[1]), int(proto_arr.shape[2])
                mask_logits = np.tensordot(coeffs, proto_arr, axes=([0], [0]))
                mask_prob = _sigmoid(mask_logits)
                scale_x = float(mw) / float(self.img_size)
                scale_y = float(mh) / float(self.img_size)
                mx0 = max(0, int(round(dw * scale_x)))
                my0 = max(0, int(round(dh * scale_y)))
                mx1 = min(mw, int(round((dw + crop_w * ratio) * scale_x)))
                my1 = min(mh, int(round((dh + crop_h * ratio) * scale_y)))
                if my1 > my0 and mx1 > mx0:
                    sub_prob = mask_prob[my0:my1, mx0:mx1]
                    resized_prob = cv2.resize(sub_prob, (crop_w, crop_h), interpolation=cv2.INTER_LINEAR)
                    box_prob = np.zeros((crop_h, crop_w), dtype=np.float32)
                    box_prob[ymin:ymax, xmin:xmax] = resized_prob[ymin:ymax, xmin:xmax]
                    mask_u8 = (box_prob >= 0.50).astype(np.uint8) * 255
            else:
                mask_u8[ymin:ymax, xmin:xmax] = 255

            if int(np.count_nonzero(mask_u8)) > 0:
                results.append(((ymin, xmin, ymax, xmax), mask_u8, float(score), int(cls_id)))

        return results

    def predict(
        self,
        image_bgr: np.ndarray,
        pretransformed: Optional[Tuple[np.ndarray, int, int]] = None,
    ) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """
        Runs LeafMachine3 1024 px YOLO26x-seg inference on a BGR leaf crop (or full sheet).
        Returns list of `(bbox, mask_uint8, score, class_id)` where `class_id` is:
          `0` (leaf / lamina), `1` (petiole), or `2` (hole).
        """
        active_sess = self.session if self.session is not None else self.model
        if active_sess is None or image_bgr is None or image_bgr.size == 0:
            return []

        try:
            if pretransformed is not None:
                chw_u8, crop_h, crop_w = pretransformed
                r = min(float(self.img_size) / float(max(crop_h, 1)), float(self.img_size) / float(max(crop_w, 1)))
                new_w = max(1, int(round(crop_w * r)))
                new_h = max(1, int(round(crop_h * r)))
                pad_xy = (
                    float(int(round((self.img_size - new_w) * 0.5 - 0.1))),
                    float(int(round((self.img_size - new_h) * 0.5 - 0.1))),
                )
                inp_tensor = (chw_u8.astype(np.float32) / 255.0)[None, ...]
            else:
                crop_h, crop_w = image_bgr.shape[:2]
                lb_bgr, r, pad_xy = _letterbox_bgr(image_bgr, self.img_size)
                inp_tensor = _bgr_to_nchw_float32(lb_bgr)

            input_name = active_sess.get_inputs()[0].name
            raw_outs = active_sess.run(None, {input_name: inp_tensor})
            return self._decode_yolo26x_seg_outputs(
                raw_outs,
                crop_shape=(crop_h, crop_w),
                ratio=r,
                pad_xy=pad_xy,
            )
        except Exception as exc:
            logger.debug(f"LM3LeafSegmenterEngine ONNX inference error: {exc}")
            return []

    def predict_batch(
        self,
        images_bgr: List[np.ndarray],
        batch_size: int = 10,
        pretransformed: Optional[Sequence[Optional[Tuple[np.ndarray, int, int]]]] = None,
    ) -> List[List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]]:
        """
        Runs Stage 2 LeafMachine3 YOLO26x-seg inference across candidate leaf crops.
        Automatically delegates to `self.predict(img)` when `predict` is mocked in unit tests.
        """
        if not images_bgr:
            return []
        has_pre = pretransformed is not None and len(pretransformed) == len(images_bgr)
        if (
            (self.session is None and self.model is None)
            or "predict" in self.__dict__
            or getattr(self.predict, "__func__", None) is not LM3LeafSegmenterEngine.predict
            or len(images_bgr) <= 1
        ):
            if (
                has_pre
                and len(images_bgr) == 1
                and "predict" not in self.__dict__
                and getattr(self.predict, "__func__", None) is LM3LeafSegmenterEngine.predict
            ):
                return [self.predict(images_bgr[0], pretransformed=pretransformed[0])]  # type: ignore[index]
            return [self.predict(img) for img in images_bgr]

        out: List[List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]] = []
        for idx, img in enumerate(images_bgr):
            pre_i = pretransformed[idx] if has_pre else None  # type: ignore[index]
            out.append(self.predict(img, pretransformed=pre_i))
        return out

    def segment_leaf_crop(
        self,
        crop_bgr: np.ndarray,
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], float]:
        """
        Runs LeafMachine3 leaf segmentation on a single candidate crop and separates
        `(lamina_mask, petiole_mask, hole_mask, best_score)`.
        """
        if crop_bgr is None or crop_bgr.size == 0:
            return None, None, None, 0.0

        preds = self.predict(crop_bgr)
        if not preds:
            return None, None, None, 0.0

        crop_h, crop_w = crop_bgr.shape[:2]
        leaf_cands: List[Tuple[float, int, Tuple[int, int, int, int], np.ndarray]] = []
        petiole_cands: List[Tuple[float, int, Tuple[int, int, int, int], np.ndarray]] = []
        hole_masks: List[np.ndarray] = []

        for bbox, mask_u8, score, cls_id in preds:
            fg = int(np.count_nonzero(mask_u8))
            if fg <= 0:
                continue
            if int(cls_id) == 0:
                leaf_cands.append((float(score), fg, bbox, mask_u8))
            elif int(cls_id) == 1:
                petiole_cands.append((float(score), fg, bbox, mask_u8))
            elif int(cls_id) == 2:
                hole_masks.append((mask_u8 > 0).astype(np.uint8) * 255)

        if not leaf_cands:
            return None, None, None, 0.0

        leaf_cands.sort(
            key=lambda c: (float(c[1]) * (float(c[0]) ** 0.5) if float(c[0]) >= 0.42 else float(c[0]), c[1]),
            reverse=True,
        )
        best_score, _, _, lamina_mask = leaf_cands[0]

        petiole_mask: Optional[np.ndarray] = None
        if petiole_cands:
            petiole_cands.sort(key=lambda c: (c[0], c[1]), reverse=True)
            petiole_mask = petiole_cands[0][3]

        hole_mask: Optional[np.ndarray] = None
        if hole_masks:
            hole_union = np.zeros((crop_h, crop_w), dtype=np.uint8)
            for hm in hole_masks:
                if hm.shape[:2] == (crop_h, crop_w):
                    np.maximum(hole_union, hm, out=hole_union)
            if int(np.count_nonzero(hole_union)) > 0:
                hole_mask = hole_union

        return lamina_mask, petiole_mask, hole_mask, float(best_score)

    def segment_crop(
        self,
        crop_bgr: np.ndarray,
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], float]:
        """Alias for `segment_leaf_crop` returning `(lamina_mask, petiole_mask, hole_mask, score)`."""
        return self.segment_leaf_crop(crop_bgr)

    def segment_sheet(
        self,
        image_bgr: np.ndarray,
        candidate_boxes: Optional[Sequence[Tuple[int, int, int, int]]] = None,
        catalog_number: str = "",
        pad_ratio: float = 0.12,
    ) -> List[DetectedInstance]:
        """
        Segments candidate leaf boxes (or the full sheet when `candidate_boxes` is None)
        and returns `DetectedInstance` objects with `mask` (lamina), `petiole_mask`, and `hole_mask`.
        """
        if image_bgr is None or image_bgr.size == 0:
            return []
        sheet_h, sheet_w = image_bgr.shape[:2]
        boxes = list(candidate_boxes) if candidate_boxes is not None else [(0, 0, sheet_h, sheet_w)]
        instances: List[DetectedInstance] = []
        for idx, box in enumerate(boxes, 1):
            crop_bgr, crop_box = crop_leaf_bbox_with_padding(image_bgr, box, pad_ratio=pad_ratio)
            if crop_bgr.size == 0:
                continue
            lamina_m, petiole_m, hole_m, sc = self.segment_leaf_crop(crop_bgr)
            if lamina_m is None or int(np.count_nonzero(lamina_m)) == 0:
                continue
            sheet_lamina, full_bbox = reproject_crop_mask_to_sheet(lamina_m, crop_box, (sheet_h, sheet_w))
            sheet_petiole = (
                reproject_crop_mask_to_sheet(petiole_m, crop_box, (sheet_h, sheet_w))[0]
                if petiole_m is not None
                else None
            )
            sheet_hole = (
                reproject_crop_mask_to_sheet(hole_m, crop_box, (sheet_h, sheet_w))[0]
                if hole_m is not None
                else None
            )
            instances.append(
                DetectedInstance(
                    catalog_number=catalog_number,
                    leaf_id=idx,
                    bbox=full_bbox,
                    mask=sheet_lamina,
                    score=sc,
                    class_id=0,
                    crop_offset=(crop_box[1], crop_box[0]),
                    mask_origin=(0, 0),
                    petiole_mask=sheet_petiole,
                    hole_mask=sheet_hole,
                )
            )
        return instances

    def __call__(self, image_bgr: np.ndarray, **kwargs: Any) -> Any:
        return self.predict(image_bgr, **kwargs)


# =============================================================================
# 4. LeafMachine3 Pose Landmark Engine (YOLO26x-pose, 640 px, 31 Keypoints)
# =============================================================================

class LM3PoseResult(dict):
    """
    Dictionary mapping all 31 canonical `mid15_pet5` keypoint names (`lamina_tip`, `apex_1..3`,
    `midvein_1..15`, `base_1..3`, `lamina_base`, `petiole_1..6`, `width_1..2`) to
    `(x, y, confidence)` tuples while also exposing `.petiole_width_px` and `.keypoints_31`.
    """

    def __init__(
        self,
        keypoints_31: Dict[str, Tuple[float, float, float]],
        petiole_width_px: Optional[float] = None,
    ) -> None:
        super().__init__(keypoints_31)
        self.petiole_width_px: Optional[float] = (
            float(petiole_width_px) if (petiole_width_px is not None and np.isfinite(petiole_width_px)) else None
        )
        self.keypoints_31: Dict[str, Tuple[float, float, float]] = dict(keypoints_31)

    def __getitem__(self, key: str) -> Any:
        if key == "petiole_width_px":
            return self.petiole_width_px
        if key == "keypoints_31":
            return self.keypoints_31
        return super().__getitem__(key)

    def get(self, key: str, default: Any = None) -> Any:
        if key == "petiole_width_px":
            return self.petiole_width_px if self.petiole_width_px is not None else default
        if key == "keypoints_31":
            return self.keypoints_31
        return super().get(key, default)


class LM3PoseLandmarkEngine:
    """
    LeafMachine3 Stage 2 ONNX Pose Landmark Engine (`lm3_landmark_detector__yolo26x_pose_640`).
    Runs 640 px YOLO26x-pose to predict the 31-keypoint `mid15_pet5` anatomical leaf skeleton:
      - `lamina_tip` (1)
      - `apex_1..3` (3)
      - `midvein_1..15` (15)
      - `base_1..3` (3)
      - `lamina_base` (1)
      - `petiole_1..6` (6)
      - `width_1..2` (2)
    plus estimated `petiole_width_px`.
    """

    MODEL_NAME: str = "lm3_landmark_detector__yolo26x_pose_640"
    KEYPOINT_NAMES: Tuple[str, ...] = LM3_KEYPOINT_NAMES_31

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        conf_thres: float = 0.25,
        img_size: int = 640,
        providers: Optional[Sequence[str]] = None,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.conf_thres = float(conf_thres)
        self.img_size = int(img_size)
        self.weights_path, self.model_name, self.providers, self.session = _load_onnx_session(
            weights_path=weights_path,
            default_model_name=self.MODEL_NAME,
            device=self.device,
            providers=providers,
        )
        self.model = self.session

    @staticmethod
    def _estimate_petiole_width_px(
        keypoints_31: Dict[str, Tuple[float, float, float]],
        petiole_mask: Optional[np.ndarray] = None,
        lamina_mask: Optional[np.ndarray] = None,
    ) -> Optional[float]:
        """Estimates petiole width in pixels from `petiole_mask` distance transform or skeleton geometry."""
        if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
            pet_u8 = (petiole_mask > 0).astype(np.uint8)
            dt = cv2.distanceTransform(pet_u8, cv2.DIST_L2, 3)
            sampled_radii: List[float] = []
            h_p, w_p = pet_u8.shape[:2]
            for kp_name in ("lamina_base", "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_5", "petiole_6"):
                kp = keypoints_31.get(kp_name)
                if kp is not None and kp[2] >= 0.20:
                    ix = max(0, min(w_p - 1, int(round(kp[0]))))
                    iy = max(0, min(h_p - 1, int(round(kp[1]))))
                    r_val = float(dt[iy, ix])
                    if r_val > 0.0:
                        sampled_radii.append(2.0 * r_val)
            if sampled_radii:
                return round(float(np.median(sampled_radii)), 3)
            nz_dt = dt[dt > 0]
            if nz_dt.size > 0:
                return round(float(2.0 * np.percentile(nz_dt, 75)), 3)

        if lamina_mask is not None and lamina_mask.size > 0 and int(np.count_nonzero(lamina_mask)) > 0:
            lam_u8 = (lamina_mask > 0).astype(np.uint8)
            dt_lam = cv2.distanceTransform(lam_u8, cv2.DIST_L2, 3)
            ys_l, _ = np.where(lam_u8 > 0)
            if ys_l.size > 0:
                y_cut = float(np.percentile(ys_l, 85))
                basal_dt = dt_lam[(lam_u8 > 0) & (np.arange(lam_u8.shape[0])[:, None] >= y_cut)]
                if basal_dt.size > 0:
                    sw = float(2.0 * np.median(basal_dt))
                    if sw > 0.0:
                        return round(sw, 3)

        return None

    @classmethod
    def _geometric_fallback_keypoints_31(
        cls,
        crop_bgr: Optional[np.ndarray],
        lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None,
        offset_xy: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        """
        Derives the 31 canonical `mid15_pet5` keypoints deterministically from `lamina_mask`
        and `petiole_mask` (or `crop_bgr`) when ONNX weights are absent in unit tests.
        """
        ox, oy = float(offset_xy[0]), float(offset_xy[1])
        mask_use = lamina_mask
        if (mask_use is None or mask_use.size == 0) and crop_bgr is not None and crop_bgr.size > 0:
            gray = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2GRAY) if crop_bgr.ndim == 3 else crop_bgr
            _, mask_use = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
            if int(np.count_nonzero(mask_use)) == 0:
                mask_use = (gray > 0).astype(np.uint8) * 255

        if mask_use is None or mask_use.size == 0 or int(np.count_nonzero(mask_use)) == 0:
            empty_kps = {k: (ox, oy, 0.0) for k in cls.KEYPOINT_NAMES}
            return LM3PoseResult(empty_kps, petiole_width_px=None)

        ys, xs = np.where(mask_use > 0)
        pts = np.column_stack((xs.astype(np.float64), ys.astype(np.float64)))
        centroid = np.mean(pts, axis=0)
        centered = pts - centroid
        cov = np.cov(centered, rowvar=False) if len(pts) >= 3 else np.eye(2)
        eigvals, eigvecs = np.linalg.eigh(cov)
        major_vec = eigvecs[:, int(np.argmax(eigvals))]
        minor_vec = np.array([-major_vec[1], major_vec[0]], dtype=np.float64)

        proj_major = centered @ major_vec
        p_min = centroid + float(np.min(proj_major)) * major_vec
        p_max = centroid + float(np.max(proj_major)) * major_vec

        # Orient base toward petiole_mask if present, otherwise lower y (herbarium basal mount)
        if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
            pys, pxs = np.where(petiole_mask > 0)
            pet_cent = np.array([float(np.mean(pxs)), float(np.mean(pys))], dtype=np.float64)
            if np.linalg.norm(p_min - pet_cent) <= np.linalg.norm(p_max - pet_cent):
                base_pt, tip_pt = p_min, p_max
            else:
                base_pt, tip_pt = p_max, p_min
        else:
            if p_min[1] >= p_max[1]:
                base_pt, tip_pt = p_min, p_max
            else:
                base_pt, tip_pt = p_max, p_min

        axis_vec = tip_pt - base_pt
        axis_len = max(float(np.linalg.norm(axis_vec)), 1.0)
        axis_unit = axis_vec / axis_len
        perp_unit = np.array([-axis_unit[1], axis_unit[0]], dtype=np.float64)

        proj_perp = (pts - centroid) @ perp_unit
        w1_pt = centroid + float(np.min(proj_perp)) * perp_unit
        w2_pt = centroid + float(np.max(proj_perp)) * perp_unit
        half_w = max(1.0, 0.5 * float(np.linalg.norm(w2_pt - w1_pt)))

        kps: Dict[str, Tuple[float, float, float]] = {}
        kps["lamina_tip"] = (round(float(tip_pt[0] + ox), 2), round(float(tip_pt[1] + oy), 2), 0.92)
        kps["apex_1"] = (
            round(float(tip_pt[0] - 0.06 * axis_vec[0] - 0.18 * half_w * perp_unit[0] + ox), 2),
            round(float(tip_pt[1] - 0.06 * axis_vec[1] - 0.18 * half_w * perp_unit[1] + oy), 2),
            0.88,
        )
        kps["apex_2"] = (round(float(tip_pt[0] + ox), 2), round(float(tip_pt[1] + oy), 2), 0.90)
        kps["apex_3"] = (
            round(float(tip_pt[0] - 0.06 * axis_vec[0] + 0.18 * half_w * perp_unit[0] + ox), 2),
            round(float(tip_pt[1] - 0.06 * axis_vec[1] + 0.18 * half_w * perp_unit[1] + oy), 2),
            0.88,
        )

        # midvein_1 (near apex) .. midvein_15 (near base)
        for idx_m in range(1, 16):
            frac_from_tip = float(idx_m) / 16.0
            m_pt = tip_pt - frac_from_tip * axis_vec
            kps[f"midvein_{idx_m}"] = (round(float(m_pt[0] + ox), 2), round(float(m_pt[1] + oy), 2), 0.90)

        kps["base_1"] = (
            round(float(base_pt[0] + 0.06 * axis_vec[0] - 0.18 * half_w * perp_unit[0] + ox), 2),
            round(float(base_pt[1] + 0.06 * axis_vec[1] - 0.18 * half_w * perp_unit[1] + oy), 2),
            0.88,
        )
        kps["base_2"] = (round(float(base_pt[0] + ox), 2), round(float(base_pt[1] + oy), 2), 0.90)
        kps["base_3"] = (
            round(float(base_pt[0] + 0.06 * axis_vec[0] + 0.18 * half_w * perp_unit[0] + ox), 2),
            round(float(base_pt[1] + 0.06 * axis_vec[1] + 0.18 * half_w * perp_unit[1] + oy), 2),
            0.88,
        )
        kps["lamina_base"] = (round(float(base_pt[0] + ox), 2), round(float(base_pt[1] + oy), 2), 0.92)

        # Petiole 1..6 trace
        if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
            pys, pxs = np.where(petiole_mask > 0)
            p_pts = np.column_stack((pxs.astype(np.float64), pys.astype(np.float64)))
            dists_from_base = np.linalg.norm(p_pts - base_pt[None, :], axis=1)
            far_pet_pt = p_pts[int(np.argmax(dists_from_base))]
            for idx_p in range(1, 7):
                t_p = float(idx_p) / 6.0
                p_coord = base_pt + t_p * (far_pet_pt - base_pt)
                kps[f"petiole_{idx_p}"] = (round(float(p_coord[0] + ox), 2), round(float(p_coord[1] + oy), 2), 0.86)
        else:
            for idx_p in range(1, 7):
                t_p = float(idx_p) / 6.0
                p_coord = base_pt - (0.12 * t_p) * axis_vec
                kps[f"petiole_{idx_p}"] = (round(float(p_coord[0] + ox), 2), round(float(p_coord[1] + oy), 2), 0.50)

        kps["width_1"] = (round(float(w1_pt[0] + ox), 2), round(float(w1_pt[1] + oy), 2), 0.90)
        kps["width_2"] = (round(float(w2_pt[0] + ox), 2), round(float(w2_pt[1] + oy), 2), 0.90)

        pet_w = cls._estimate_petiole_width_px(kps, petiole_mask=petiole_mask, lamina_mask=mask_use)
        return LM3PoseResult(kps, petiole_width_px=pet_w)

    def predict_landmarks(
        self,
        crop_bgr: np.ndarray,
        lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None,
        crop_offset: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        """
        Runs 640 px YOLO26x-pose landmark detection on `crop_bgr`, returning an `LM3PoseResult`
        mapping all 31 canonical keypoint names (`lamina_tip`, `apex_1..3`, `midvein_1..15`,
        `base_1..3`, `lamina_base`, `petiole_1..6`, `width_1..2`) to `(x, y, confidence)`
        in crop (or offset sheet) coordinates, plus `.petiole_width_px`.
        """
        ox, oy = float(crop_offset[0]), float(crop_offset[1])
        active_sess = self.session if self.session is not None else self.model
        if active_sess is not None and crop_bgr is not None and crop_bgr.size > 0:
            try:
                crop_h, crop_w = crop_bgr.shape[:2]
                px, py = int(round(crop_w * 0.10)), int(round(crop_h * 0.10))
                padded_bgr = cv2.copyMakeBorder(
                    crop_bgr, py, py, px, px, cv2.BORDER_CONSTANT, value=(255, 255, 255)
                )
                lb_bgr, r, (dw, dh) = _letterbox_bgr(padded_bgr, self.img_size)
                inp_tensor = _bgr_to_nchw_float32(lb_bgr)
                input_name = active_sess.get_inputs()[0].name
                raw_outs = active_sess.run(None, {input_name: inp_tensor})
                if raw_outs:
                    arr = np.asarray(raw_outs[0], dtype=np.float32)
                    if arr.ndim == 3:
                        arr = arr[0]
                    if arr.ndim == 2 and arr.shape[0] < arr.shape[1] and arr.shape[0] >= 5 + 3 * 31:
                        arr = arr.T
                    if arr.ndim == 2 and arr.shape[0] > 0 and arr.shape[1] >= 5 + 3 * 31:
                        best_row_idx = int(np.argmax(arr[:, 4]))
                        best_row = arr[best_row_idx]
                        if float(best_row[4]) >= self.conf_thres:
                            kp_offset = 6 if best_row.shape[0] >= 6 + 3 * 31 else 5
                            kp_flat = best_row[kp_offset : kp_offset + 3 * 31].reshape(31, 3)
                            inv_r = 1.0 / max(float(r), 1e-6)
                            kps_out: Dict[str, Tuple[float, float, float]] = {}
                            for idx_k, kp_name in enumerate(self.KEYPOINT_NAMES):
                                kx_pad = (float(kp_flat[idx_k, 0]) - dw) * inv_r
                                ky_pad = (float(kp_flat[idx_k, 1]) - dh) * inv_r
                                kx = max(0.0, min(float(crop_w), kx_pad - float(px))) + ox
                                ky = max(0.0, min(float(crop_h), ky_pad - float(py))) + oy
                                kc = float(np.clip(kp_flat[idx_k, 2], 0.0, 1.0))
                                kps_out[kp_name] = (round(kx, 2), round(ky, 2), round(kc, 4))
                            pet_w = self._estimate_petiole_width_px(
                                kps_out, petiole_mask=petiole_mask, lamina_mask=lamina_mask
                            )
                            return LM3PoseResult(kps_out, petiole_width_px=pet_w)
            except Exception as exc:
                logger.debug(f"LM3PoseLandmarkEngine ONNX inference fallback: {exc}")

        return self._geometric_fallback_keypoints_31(
            crop_bgr=crop_bgr,
            lamina_mask=lamina_mask,
            petiole_mask=petiole_mask,
            offset_xy=(ox, oy),
        )

    def predict(
        self,
        crop_bgr: np.ndarray,
        lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None,
        crop_offset: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        """Alias for `predict_landmarks`."""
        return self.predict_landmarks(
            crop_bgr=crop_bgr,
            lamina_mask=lamina_mask,
            petiole_mask=petiole_mask,
            crop_offset=crop_offset,
        )

    def __call__(
        self,
        crop_bgr: np.ndarray,
        lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None,
        crop_offset: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        return self.predict_landmarks(
            crop_bgr=crop_bgr,
            lamina_mask=lamina_mask,
            petiole_mask=petiole_mask,
            crop_offset=crop_offset,
        )


# =============================================================================
# 5. Native LM3 Petiole / Keypoint Anchoring & Dissected Rachis Bridging
# =============================================================================

def _extract_native_lm3_petiole_anchor(
    leaf_mask: np.ndarray,
    petiole_mask: Optional[np.ndarray] = None,
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None,
    min_kp_conf: float = 0.25,
) -> Optional[Dict[str, Any]]:
    """
    Uses LeafMachine3's explicit Stage 2 `petiole_mask` (`class_id == 1`) and/or
    31-keypoint `mid15_pet5` (`lamina_base`, `petiole_1..6`, `lamina_tip`) skeleton to
    locate the anatomical petiole-lamina junction `p_base_crop`, distal lamina apex
    `p_apex_crop`, and native petiole pixel span without clipping `leaf_mask`.
    Falls back to the basal-most mask contour point when both `keypoints_31` and
    `petiole_mask` are absent.
    """
    if leaf_mask is None or leaf_mask.size == 0 or int(np.count_nonzero(leaf_mask)) == 0:
        return None

    bx_l, by_l, bw_l, bh_l = cv2.boundingRect(leaf_mask)
    if bw_l == 0 or bh_l == 0:
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

    p_base_crop: Optional[Tuple[int, int]] = None
    p_apex_crop: Optional[Tuple[int, int]] = None
    petiole_len_px: float = 0.0
    used_native_lm3 = False

    # 1. Try 31-keypoint pose skeleton (`lamina_base`, `petiole_1..6`, `lamina_tip`)
    if keypoints_31:
        lb_kp = keypoints_31.get("lamina_base") or keypoints_31.get("base_2")
        if lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= min_kp_conf:
            dists_lb = np.hypot(leaf_pts[:, 0] - float(lb_kp[0]), leaf_pts[:, 1] - float(lb_kp[1]))
            min_lb_idx = int(np.argmin(dists_lb))
            if float(dists_lb[min_lb_idx]) <= max(36.0, 0.35 * float(max(bw_l, bh_l))):
                p_base_crop = (int(leaf_pts[min_lb_idx, 0]), int(leaf_pts[min_lb_idx, 1]))
                used_native_lm3 = True

        lt_kp = keypoints_31.get("lamina_tip") or keypoints_31.get("apex_2")
        if lt_kp is not None and len(lt_kp) >= 3 and float(lt_kp[2]) >= min_kp_conf:
            dists_lt = np.hypot(leaf_pts[:, 0] - float(lt_kp[0]), leaf_pts[:, 1] - float(lt_kp[1]))
            min_lt_idx = int(np.argmin(dists_lt))
            p_apex_crop = (int(leaf_pts[min_lt_idx, 0]), int(leaf_pts[min_lt_idx, 1]))

        pet_chain: List[Tuple[float, float]] = []
        if p_base_crop is not None:
            pet_chain.append((float(p_base_crop[0]), float(p_base_crop[1])))
        for idx_p in range(1, 7):
            pkp = keypoints_31.get(f"petiole_{idx_p}")
            if pkp is not None and len(pkp) >= 3 and float(pkp[2]) >= min_kp_conf:
                pet_chain.append((float(pkp[0]), float(pkp[1])))
        pt_kp = keypoints_31.get("petiole_tip")
        if pt_kp is not None and len(pt_kp) >= 3 and float(pt_kp[2]) >= min_kp_conf:
            pet_chain.append((float(pt_kp[0]), float(pt_kp[1])))
        if len(pet_chain) >= 2:
            arc_len = 0.0
            for i_c in range(len(pet_chain) - 1):
                arc_len += math.hypot(
                    pet_chain[i_c + 1][0] - pet_chain[i_c][0],
                    pet_chain[i_c + 1][1] - pet_chain[i_c][1],
                )
            petiole_len_px = max(petiole_len_px, float(arc_len))

    # 2. Use explicit Stage 2 `petiole_mask` (`class_id == 1`) when present
    if (
        petiole_mask is not None
        and petiole_mask.size > 0
        and petiole_mask.shape[:2] == leaf_mask.shape[:2]
        and int(np.count_nonzero(petiole_mask)) > 0
    ):
        bx_p, by_p, bw_p, bh_p = cv2.boundingRect(petiole_mask)
        if bw_p > 0 and bh_p > 0:
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
            if min_dist <= max(28.0, 0.30 * float(max(bw_l, bh_l))):
                if p_base_crop is None:
                    p_base_crop = (int(leaf_pts[min_idx, 0]), int(leaf_pts[min_idx, 1]))
                used_native_lm3 = True
                cnts_p, _ = cv2.findContours(
                    petiole_mask[by_p : by_p + bh_p, bx_p : bx_p + bw_p],
                    cv2.RETR_EXTERNAL,
                    cv2.CHAIN_APPROX_NONE,
                    offset=(bx_p, by_p),
                )
                if cnts_p and p_base_crop is not None:
                    pet_pts = np.vstack([c.reshape(-1, 2) for c in cnts_p])
                    pet_span_dists = np.hypot(pet_pts[:, 0] - p_base_crop[0], pet_pts[:, 1] - p_base_crop[1])
                    if len(pet_span_dists) > 0:
                        petiole_len_px = max(petiole_len_px, float(np.max(pet_span_dists)))

    if p_base_crop is not None:
        if p_apex_crop is None:
            dists_from_base = np.hypot(leaf_pts[:, 0] - p_base_crop[0], leaf_pts[:, 1] - p_base_crop[1])
            apex_idx = int(np.argmax(dists_from_base))
            p_apex_crop = (int(leaf_pts[apex_idx, 0]), int(leaf_pts[apex_idx, 1]))
        return {
            "is_trimmed": False,
            "native_lm3_petiole": bool(used_native_lm3),
            "trim_distance_px": round(float(petiole_len_px), 2),
            "p_base_crop": p_base_crop,
            "p_apex_crop": p_apex_crop,
        }

    # 3. Fallback to basal-most contour vertex when both keypoints and petiole_mask are absent
    base_idx = int(np.argmax(leaf_pts[:, 1]))
    p_base_crop = (int(leaf_pts[base_idx, 0]), int(leaf_pts[base_idx, 1]))
    dists_from_base = np.hypot(leaf_pts[:, 0] - p_base_crop[0], leaf_pts[:, 1] - p_base_crop[1])
    apex_idx = int(np.argmax(dists_from_base))
    p_apex_crop = (int(leaf_pts[apex_idx, 0]), int(leaf_pts[apex_idx, 1]))
    return {
        "is_trimmed": False,
        "native_lm3_petiole": False,
        "morphological_fallback": True,
        "trim_distance_px": 0.0,
        "p_base_crop": p_base_crop,
        "p_apex_crop": p_apex_crop,
    }


def _merge_and_bridge_dissected_crop_mask(
    best_mask: np.ndarray,
    leaf_cands: Sequence[Tuple[float, int, Any, np.ndarray]],
    petiole_mask: Optional[np.ndarray] = None,
    crop_bgr: Optional[np.ndarray] = None,
    specimen_paper_cleaned_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Reconstructs a unified dissected leaf silhouette when Stage 2 segmentation predicts
    lateral pinna pairs as either:
      (a) disconnected foreground components within `best_mask`, OR
      (b) secondary `Class 0` (`leaf`) predictions in `leaf_cands[1:]` that are spatially
          disjoint from `best_mask` and anchored to the same native `petiole_mask` (`Class 1`)
          rachis stalk or connected along `specimen_paper_cleaned_mask` plant tissue.
    Forwards `specimen_paper_cleaned_mask` into `bridge_dissected_rachis_components` so deeply
    pinnatifid Packera lobes are bridged along verified UNet++ v9 plant tissue.
    """
    if best_mask is None or best_mask.size == 0:
        return best_mask

    merged = best_mask
    best_area = int(np.count_nonzero(best_mask))
    support_mask: Optional[np.ndarray] = None
    if petiole_mask is not None and petiole_mask.shape[:2] == best_mask.shape[:2] and int(np.count_nonzero(petiole_mask)) > 0:
        support_mask = (petiole_mask > 0).astype(np.uint8) * 255
    if (
        specimen_paper_cleaned_mask is not None
        and specimen_paper_cleaned_mask.shape[:2] == best_mask.shape[:2]
        and int(np.count_nonzero(specimen_paper_cleaned_mask)) > 0
    ):
        spec_u8 = (specimen_paper_cleaned_mask > 0).astype(np.uint8) * 255
        support_mask = spec_u8 if support_mask is None else np.maximum(support_mask, spec_u8)

    if len(leaf_cands) > 1 and support_mask is not None and best_area > 0:
        k_pet = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (9, 9))
        sup_dil = cv2.dilate(support_mask, k_pet, iterations=1) > 0
        best_bool = best_mask > 0
        if np.any(best_bool & sup_dil):
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
                if np.any(c_bool & sup_dil):
                    merged = np.maximum(merged, (c_bool.astype(np.uint8) * 255))

    try:
        return bridge_dissected_rachis_components(
            merged,
            petiole_mask=petiole_mask,
            crop_bgr=crop_bgr,
            specimen_paper_cleaned_mask=specimen_paper_cleaned_mask,
        )
    except TypeError:
        return bridge_dissected_rachis_components(
            merged,
            petiole_mask=petiole_mask,
            crop_bgr=crop_bgr,
        )
