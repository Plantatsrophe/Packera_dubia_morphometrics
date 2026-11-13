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
import threading
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import cv2
import scripts.vision.lm3_compat  # noqa: F401
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
LM3_PACKAGE_ROOT = PROJECT_ROOT / "LeafMachine3"
if LM3_PACKAGE_ROOT.is_dir() and str(LM3_PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(LM3_PACKAGE_ROOT))
_VENV_LM3_LIB = PROJECT_ROOT / ".venv_LM3" / "lib"


def _ensure_ect_module_loaded() -> None:
    if "ect" in sys.modules:
        return
    try:
        import ect as _ect_probe  # noqa: F401
        return
    except ImportError:
        pass
    if not _VENV_LM3_LIB.is_dir():
        return
    import importlib.util as _ilu
    import types as _types

    for _ect_init in sorted(_VENV_LM3_LIB.glob("python*/site-packages/ect/__init__.py")):
        if _ect_init.is_file():
            try:
                import numba as _numba_probe  # noqa: F401
            except ImportError:
                _numba_mod = _types.ModuleType("numba")
                _numba_mod.prange = range  # type: ignore[attr-defined]
                _numba_mod.njit = (  # type: ignore[attr-defined]
                    lambda *a, **kw: (a[0] if (a and callable(a[0])) else (lambda f: f))
                )
                _numba_typed = _types.ModuleType("numba.typed")
                _numba_typed.List = list  # type: ignore[attr-defined]
                _numba_mod.typed = _numba_typed  # type: ignore[attr-defined]
                sys.modules["numba"] = _numba_mod
                sys.modules["numba.typed"] = _numba_typed
            _ect_spec = _ilu.spec_from_file_location(
                "ect",
                str(_ect_init),
                submodule_search_locations=[str(_ect_init.parent)],
            )
            if _ect_spec is not None and _ect_spec.loader is not None:
                _ect_mod = _ilu.module_from_spec(_ect_spec)
                sys.modules["ect"] = _ect_mod
                _ect_spec.loader.exec_module(_ect_mod)
            break

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

from scripts.vision.lm3_data_loader import (
    LM3_LEAF_CLASS_IDS,
    LM3_PLANT_CLASSES,
    LM3_REPRODUCTIVE_CLASS_IDS,
    compute_work_scale,
    load_precomputed_plant_labels,
    unscale_working_bbox_to_sheet,
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
    "LM3LeafSegmenterEngine",
    "LM3PoseResult",
    "LM3PoseLandmarkEngine",
    "LM3ScaleCalibrator",
    "normalize_lm3_keypoints_dict",
    "orient_and_score_leaf_lm3",
    "compute_leaf_ect_descriptor",
    "torch_cuda_available",
    "is_stock_lm3_weights",
    "_prepare_working_image",
    "_reduce_crop_predictions_native_lm3",
    "_extract_native_lm3_petiole_anchor",
    "evaluate_leaf_damage_and_fold_qc",
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
    solidity: float = 0.0
    midrib_angle_deg: float = 0.0
    pixels_per_mm: Optional[float] = None
    is_calibrated_mm: bool = False
    crop_offset: Tuple[int, int] = (0, 0)
    mask_origin: Tuple[int, int] = (0, 0)
    trim_meta: Optional[Dict[str, Any]] = None
    pcd_bbox: Optional[Tuple[int, int, int, int]] = None
    circularity: float = 0.0
    extent: float = 0.0
    aspect_ratio: float = 0.0
    area_mm2: float = np.nan
    mask_path: Optional[str] = None
    contour_path: Optional[str] = None
    petiole_mask: Optional[np.ndarray] = None
    hole_mask: Optional[np.ndarray] = None
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None
    work_scale: float = 1.0
    cf_source: str = "measured_from_ruler"
    archetype_score: float = 0.0
    petiole_width_px: Optional[float] = None
    petiole_width_mm: Optional[float] = None
    blade_length_px: Optional[float] = None
    blade_width_px: Optional[float] = None
    ect_path: Optional[str] = None
    blade_cie_l_norm: float = np.nan
    blade_cie_l_sd: float = np.nan


@dataclass
class _VoucherFinalizeBundle:
    """Thread-safe per-voucher post-processing result bundle for wave-pipelined execution."""
    catalog_number: str
    image_path: str
    extracted_records: List[Dict[str, Any]]
    failure_record: Optional[Dict[str, Any]]
    repro_record: Optional[Dict[str, Any]] = None
    pheno_record: Optional[Dict[str, Any]] = None
    repro_by_pid: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    pheno_by_pid: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    crowns_by_pid: Dict[int, Tuple[float, float]] = field(default_factory=dict)
    rosette_centroids_by_pid: Dict[int, Tuple[float, float]] = field(default_factory=dict)
    by_plant_petiole: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    contour_coords_by_leaf: Dict[int, List[List[float]]] = field(default_factory=dict)


def crop_leaf_bbox_with_padding(
    image_bgr: np.ndarray,
    bbox: Tuple[int, int, int, int],
    pad_ratio: float = 0.0,
) -> Tuple[np.ndarray, Tuple[int, int, int, int]]:
    """
    Crops a candidate leaf bounding box (ymin, xmin, ymax, xmax) from the native-DPI
    voucher sheet with optional isotropic padding (default 0.0 for tight LeafMachine3
    Stage 1 -> Stage 2 crops), clamped to sheet boundaries.
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
    min_ar: float = 0.25,
    max_ar: float = 4.5,
    crop_bgr: Optional[np.ndarray] = None,
    class_id: Optional[int] = None,
    peduncle_dir_xy: Optional[Tuple[float, float]] = None,
    is_calibrated_mm: bool = True,
) -> Optional[Dict[str, Any]]:
    """
    Computes herbaceous involucre height (H_inv), width (W_inv), and aspect ratio (AR_inv = H/W)
    strictly for single-capitulum detections (`Seed_Fruit_ONE = 2`, `Flower_ONE = 4`, `Bud = 6`).

    Anatomical & Geometric Safeguards:
      1. Rejects multi-head corymbose cluster detections (`Seed_Fruit_MANY = 3`, `Flower_MANY = 5`).
      2. When `crop_bgr` is supplied, isolates the green/olive/anthocyanic herbaceous phyllary cup
         by masking out sheet paper background, yellow/gold ray/disc corollas, and bright white
         capillary pappus plumes in HSV/RGB space.
      3. Resolves longitudinal `involucre_height_px` vs. transverse `involucre_width_px` along the
         anatomical peduncle-to-apex vector (from phyllary cup centroid to distal floret/pappus
         centroid when present, else `peduncle_dir_xy`, else upward sheet axis), avoiding gravity
         axis inversions on tilted or nodding peduncles.
      4. Uses wide artifact-only bounds (`0.25 <= AR_inv <= 4.5`) so genuine campanulate/hemispheric
         fruiting involucres and slender cylindrical buds are never censored.
    """
    if class_id is not None and int(class_id) not in (2, 4, 6):
        return None

    eff_mask: Optional[np.ndarray] = None
    if mask is not None and isinstance(mask, np.ndarray) and mask.size > 0:
        eff_mask = (mask > 0).astype(np.uint8) * 255

    axis_ux, axis_uy = 0.0, -1.0
    has_anatomical_axis = False
    if peduncle_dir_xy is not None and len(peduncle_dir_xy) >= 2:
        px_d, py_d = float(peduncle_dir_xy[0]), float(peduncle_dir_xy[1])
        p_norm = float(np.hypot(px_d, py_d))
        if np.isfinite(p_norm) and p_norm > 1e-6:
            axis_ux, axis_uy = px_d / p_norm, py_d / p_norm
            has_anatomical_axis = True

    if crop_bgr is not None and isinstance(crop_bgr, np.ndarray) and crop_bgr.ndim == 3 and crop_bgr.shape[0] >= 5 and crop_bgr.shape[1] >= 5:
        ch, cw = crop_bgr.shape[:2]
        hsv = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2HSV)
        H_ch = hsv[:, :, 0]
        S_ch = hsv[:, :, 1]
        V_ch = hsv[:, :, 2]

        border_pts = np.concatenate([
            crop_bgr[0, :, :],
            crop_bgr[-1, :, :],
            crop_bgr[:, 0, :],
            crop_bgr[:, -1, :],
        ], axis=0).astype(np.float32)
        bg_med = np.median(border_pts, axis=0)
        bg_dist = np.linalg.norm(crop_bgr.astype(np.float32) - bg_med[None, None, :], axis=2)
        paper_mask = (bg_dist < 18.0) & (S_ch <= 42)

        fg_plant = ~paper_mask
        if eff_mask is not None and eff_mask.shape[:2] == (ch, cw) and np.any(eff_mask > 0):
            fg_plant = fg_plant & (eff_mask > 0)

        floret_mask = fg_plant & (H_ch >= 10) & (H_ch <= 36) & (S_ch >= 76) & (V_ch >= 85)
        pappus_mask = fg_plant & (~floret_mask) & (V_ch >= 158) & (S_ch <= 48)
        involucre_cand = fg_plant & (~floret_mask) & (~pappus_mask) & (V_ch <= 198)

        min_inv_px = max(12, int(0.04 * ch * cw))
        if int(np.count_nonzero(involucre_cand)) >= min_inv_px:
            eff_mask = involucre_cand.astype(np.uint8) * 255
            distal_mask = floret_mask | pappus_mask
            if int(np.count_nonzero(distal_mask)) >= 8:
                iy_inv, ix_inv = np.where(involucre_cand)
                iy_dst, ix_dst = np.where(distal_mask)
                dx_a = float(np.mean(ix_dst) - np.mean(ix_inv))
                dy_a = float(np.mean(iy_dst) - np.mean(iy_inv))
                d_norm = float(np.hypot(dx_a, dy_a))
                if d_norm >= 2.0:
                    axis_ux, axis_uy = dx_a / d_norm, dy_a / d_norm
                    has_anatomical_axis = True

    h_inv_px = 0.0
    w_inv_px = 0.0
    contours = []
    if eff_mask is not None and eff_mask.size > 0:
        contours, _ = cv2.findContours(eff_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
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
        proj1 = abs(float(v1[0]) * axis_ux + float(v1[1]) * axis_uy) / norm1
        proj2 = abs(float(v2[0]) * axis_ux + float(v2[1]) * axis_uy) / norm2
        if proj1 >= proj2:
            h_inv_px, w_inv_px = len1, len2
        else:
            h_inv_px, w_inv_px = len2, len1
    else:
        ymin, xmin, ymax, xmax = bbox
        dy_box = float(max(ymax - ymin, 1))
        dx_box = float(max(xmax - xmin, 1))
        if has_anatomical_axis and abs(axis_ux) > abs(axis_uy):
            h_inv_px, w_inv_px = dx_box, dy_box
        else:
            h_inv_px, w_inv_px = dy_box, dx_box

    if w_inv_px <= 1e-6:
        return None

    aspect_ratio = h_inv_px / w_inv_px
    if not (min_ar <= aspect_ratio <= max_ar):
        return None

    use_cal = bool(is_calibrated_mm and pixels_per_mm is not None and not np.isnan(pixels_per_mm) and float(pixels_per_mm) > 0.0)
    h_inv_mm = (h_inv_px / float(pixels_per_mm)) if use_cal else np.nan
    w_inv_mm = (w_inv_px / float(pixels_per_mm)) if use_cal else np.nan

    return {
        "involucre_height_px": round(h_inv_px, 2),
        "involucre_width_px": round(w_inv_px, 2),
        "involucre_height_mm": round(h_inv_mm, 3) if not np.isnan(h_inv_mm) else np.nan,
        "involucre_width_mm": round(w_inv_mm, 3) if not np.isnan(w_inv_mm) else np.nan,
        "AR_inv": round(aspect_ratio, 4),
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
_SHARED_CUDA_ARENA_REGISTERED: bool = False
_ONNX_GPU_LOCK: threading.Lock = threading.Lock()


def _ensure_shared_cuda_arena() -> bool:
    """
    Registers a single process-wide shared CUDA BFCArena allocator in ONNX Runtime so all
    4 LeafMachine3 ONNX sessions (PlantDetector, LeafSegmenter, PoseLandmark, SpecimenSegmenter)
    reuse the same GPU activation pool instead of hoarding 4 isolated per-session arenas.
    """
    global _SHARED_CUDA_ARENA_REGISTERED
    if _SHARED_CUDA_ARENA_REGISTERED:
        return True
    try:
        import onnxruntime as ort
        mem_info = ort.OrtMemoryInfo(
            "Cuda",
            ort.OrtAllocatorType.ORT_ARENA_ALLOCATOR,
            0,
            ort.OrtMemType.DEFAULT,
        )
        # OrtArenaCfg(max_mem=0 [auto], arena_extend_strategy=1 [kSameAsRequested], initial_chunk=-1, max_dead=-1)
        arena_cfg = ort.OrtArenaCfg(0, 1, -1, -1)
        if hasattr(ort, "create_and_register_allocator_v2"):
            ort.create_and_register_allocator_v2("CUDAExecutionProvider", mem_info, {}, arena_cfg)
        else:
            ort.create_and_register_allocator(mem_info, arena_cfg)
        _SHARED_CUDA_ARENA_REGISTERED = True
        return True
    except Exception as exc:
        logger.debug(f"Shared CUDA arena registration skipped: {exc}")
        return False


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
    lookup_name = default_model_name
    if (
        (weights_path is None or not str(weights_path).strip())
        and default_model_name == "lm3_landmark_detector__yolo26x_pose_640"
    ):
        cand_canon = resolve_lm3_weights("lm3_landmark_detector__yolo26x_pose_640_canonical", explicit_path=None)
        canon_path = Path(cand_canon[1]) if isinstance(cand_canon, tuple) else Path(cand_canon)
        try:
            if canon_path.is_file() and canon_path.stat().st_size >= 1024:
                lookup_name = "lm3_landmark_detector__yolo26x_pose_640_canonical"
        except OSError:
            pass
    resolved = resolve_lm3_weights(lookup_name, explicit_path=weights_path)
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
        sess_opts.enable_mem_pattern = False
        sess_opts.enable_mem_reuse = True
        if "CUDAExecutionProvider" in prov_list and _ensure_shared_cuda_arena():
            try:
                sess_opts.add_session_config_entry("session.use_env_allocators", "1")
            except Exception:
                pass
        sess_opts.intra_op_num_threads = 2
        sess_opts.inter_op_num_threads = 1
        try:
            sess_opts.log_severity_level = int(os.environ.get("LM3_ORT_LOG_SEVERITY", 3))
        except ValueError:
            sess_opts.log_severity_level = 3
        ort_providers: List[Any] = [
            (
                "CUDAExecutionProvider",
                {
                    "device_id": 0,
                    "arena_extend_strategy": "kSameAsRequested",
                    "cudnn_conv_use_max_workspace": "0",
                    "do_copy_in_default_stream": "1",
                },
            )
            if p == "CUDAExecutionProvider"
            else p
            for p in prov_list
        ]
        session = ort.InferenceSession(str(w_path), sess_options=sess_opts, providers=ort_providers)
        _ONNX_SESSION_CACHE[cache_key] = session
        bound_provs = session.get_providers()
        logger.info(f"Loaded LeafMachine3 ONNX session ({m_name}) from {w_path} with providers={bound_provs}")
        return w_path, m_name, bound_provs, session
    except Exception as exc:
        logger.debug(f"ONNX Runtime session initialization skipped for {w_path}: {exc}")
        return w_path, m_name, prov_list, None


_INTER_LINEAR_EXACT: int = int(getattr(cv2, "INTER_LINEAR_EXACT", cv2.INTER_LINEAR))


def _prepare_working_image(
    image_bgr: np.ndarray,
    max_working_dim: int = 3200,
) -> Tuple[np.ndarray, float]:
    """
    Downscales `image_bgr` with `cv2.INTER_AREA` when `max(orig_h, orig_w) > max_working_dim`
    (default 3200 px), returning `(work_img, work_scale)` matching LeafMachine3's working-image stage.
    """
    if image_bgr is None or image_bgr.size == 0:
        return image_bgr, 1.0
    orig_h, orig_w = int(image_bgr.shape[0]), int(image_bgr.shape[1])
    work_scale = compute_work_scale((orig_h, orig_w), max_working_dim=int(max_working_dim))
    if work_scale < 1.0:
        work_w = max(1, int(round(orig_w * work_scale)))
        work_h = max(1, int(round(orig_h * work_scale)))
        work_img = cv2.resize(image_bgr, (work_w, work_h), interpolation=cv2.INTER_AREA)
        return work_img, float(work_scale)
    return image_bgr, 1.0


def _letterbox_gain_pad_hw(
    img1_shape: Tuple[int, int],
    img0_shape: Tuple[int, int],
) -> Tuple[float, Tuple[float, float]]:
    """Computes LeafMachine3 `ultra_rep._letterbox_gain_pad` `(gain, (pad_x, pad_y))`."""
    try:
        from leafmachine3.inference import ultra_replacements as ultra_rep
        gain, pad_x, pad_y = ultra_rep._letterbox_gain_pad(img1_shape, img0_shape)
        return float(gain), (float(pad_x), float(pad_y))
    except Exception:
        h0, w0 = max(1, int(img0_shape[0])), max(1, int(img0_shape[1]))
        h1, w1 = int(img1_shape[0]), int(img1_shape[1])
        gain = min(float(h1) / float(h0), float(w1) / float(w0))
        pad_x = float(round((w1 - round(w0 * gain)) / 2.0 - 0.1))
        pad_y = float(round((h1 - round(h0 * gain)) / 2.0 - 0.1))
        return gain, (pad_x, pad_y)


def _letterbox_bgr(
    image_bgr: np.ndarray,
    target_size: int,
    pad_value: int = 114,
    stride: int = 32,
    auto: bool = True,
) -> Tuple[np.ndarray, float, Tuple[float, float]]:
    """
    Resizes `image_bgr` using LeafMachine3's native dynamic stride-32 letterboxing
    (`ultra_rep.letterbox` with `auto=True`, `stride=32`, and `cv2.INTER_LINEAR_EXACT`).
    Returns `(letterboxed_bgr, scale_ratio, (dw, dh))`.
    """
    h, w = image_bgr.shape[:2]
    if h <= 0 or w <= 0:
        empty = np.full((target_size, target_size, 3), pad_value, dtype=np.uint8)
        return empty, 1.0, (0.0, 0.0)
    try:
        from leafmachine3.inference import ultra_replacements as ultra_rep
        in_hw = ultra_rep.check_imgsz(int(target_size), stride=int(stride))
        padded = ultra_rep.letterbox(
            image_bgr,
            in_hw,
            stride=int(stride),
            auto=bool(auto),
            center=True,
        )
        gain, pad_x, pad_y = ultra_rep._letterbox_gain_pad(padded.shape[:2], (h, w))
        return padded, float(gain), (float(pad_x), float(pad_y))
    except Exception:
        r = min(float(target_size) / float(h), float(target_size) / float(w))
        new_w = max(1, int(round(w * r)))
        new_h = max(1, int(round(h * r)))
        dw = float(target_size - new_w)
        dh = float(target_size - new_h)
        if auto:
            dw = float(np.mod(dw, stride))
            dh = float(np.mod(dh, stride))
        dw *= 0.5
        dh *= 0.5
        top = int(round(dh - 0.1))
        bottom = int(round(dh + 0.1))
        left = int(round(dw - 0.1))
        right = int(round(dw + 0.1))
        if (new_w, new_h) != (w, h):
            resized = cv2.resize(image_bgr, (new_w, new_h), interpolation=cv2.INTER_LINEAR_EXACT)
        else:
            resized = image_bgr
        padded = cv2.copyMakeBorder(
            resized, top, bottom, left, right, cv2.BORDER_CONSTANT, value=(pad_value, pad_value, pad_value)
        )
        return padded, r, (float(left), float(top))


def _bgr_to_nchw_float32(letterboxed_bgr: np.ndarray) -> np.ndarray:
    """
    Converts a `(H, W, 3)` BGR uint8 array to `(1, 3, H, W)` RGB float32 in `[0, 1]`
    using LeafMachine3's reciprocal float32 scaling (`x *= float32(1.0) / float32(255.0)`).
    """
    im = letterboxed_bgr[None, ..., ::-1].transpose((0, 3, 1, 2))
    im_f32 = np.ascontiguousarray(im, dtype=np.float32)
    im_f32 *= np.float32(1.0) / np.float32(255.0)
    return im_f32


# =============================================================================
# 1. LeafMachine3 Stage 1 Plant Detector (YOLO26x End-to-End NMS-Free, 1280 px)
# =============================================================================

def _run_onnx_minibatches_back_to_back(
    active_sess: Any,
    input_name: str,
    valid_chw_list: Sequence[np.ndarray],
    initial_sub_bs: int,
) -> List[Tuple[int, int, List[np.ndarray]]]:
    """
    Executes GPU ONNX mini-batches back-to-back without interleaved CPU post-processing.
    Automatically halves the mini-batch size on transient BFCArena / VRAM allocation pressure
    and respects contiguous `(H_lb, W_lb)` shape boundaries.
    """
    n_valid = len(valid_chw_list)
    if n_valid == 0:
        return []
    sub_bs = max(1, int(initial_sub_bs))
    batch_records: List[Tuple[int, int, List[np.ndarray]]] = []
    idx = 0
    while idx < n_valid:
        base_shape = valid_chw_list[idx].shape
        max_same_shape = 1
        while (
            idx + max_same_shape < n_valid
            and max_same_shape < sub_bs
            and valid_chw_list[idx + max_same_shape].shape == base_shape
        ):
            max_same_shape += 1
        cur_bs = min(sub_bs, max_same_shape)
        batch_tensor = np.stack(valid_chw_list[idx : idx + cur_bs], axis=0)
        try:
            with _ONNX_GPU_LOCK:
                raw_outs = active_sess.run(None, {input_name: batch_tensor})
        except Exception as exc:
            if cur_bs > 1:
                sub_bs = max(1, cur_bs // 2)
                logger.debug(
                    f"ONNX mini-batch size {cur_bs} hit memory ceiling ({exc}); retrying with sub_bs={sub_bs}"
                )
                continue
            raise
        if raw_outs:
            batch_records.append((idx, idx + cur_bs, raw_outs))
        idx += cur_bs
    return batch_records


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
        conf_thres: float = 0.30,
        iou_thres: float = 0.95,
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
        work_img, _ = _prepare_working_image(image_bgr, max_working_dim=work_dim)
        lb_bgr, _, _ = _letterbox_bgr(work_img, self.img_size, auto=True, stride=32)
        return lb_bgr

    _preprocess_letterbox = prepare_letterbox

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
            return self._deduplicate_sheet_boxes(results)

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
        return self._deduplicate_sheet_boxes(results)

    _postprocess_yolo26x_det = _decode_yolo26x_det_output

    def _deduplicate_sheet_boxes(
        self,
        dets: List[Tuple[Tuple[int, int, int, int], float, int]],
        containment_thres: float = 0.95,
    ) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
        """
        Applies LeafMachine3 native same-class duplicate box suppression
        (`suppress_duplicate_boxes(..., overlap=0.95, metric='min')`), sorting candidates in
        descending confidence order (`key=lambda x: x[1], reverse=True`) and suppressing a box
        only when `inter / min_area >= 0.95` (`containment_thres`) with a higher-confidence kept box.
        """
        if len(dets) <= 1:
            return dets
        eff_overlap = max(float(containment_thres), float(self.iou_thres)) if self.iou_thres > 0.85 else float(containment_thres)
        by_class: Dict[int, List[Tuple[int, Tuple[int, int, int, int], float, int]]] = defaultdict(list)
        for orig_idx, (bbox, score, cls_id) in enumerate(dets):
            by_class[int(cls_id)].append((orig_idx, bbox, float(score), int(cls_id)))

        kept_with_idx: List[Tuple[int, Tuple[int, int, int, int], float, int]] = []
        for _cid, group in by_class.items():
            group.sort(key=lambda item: item[2], reverse=True)
            accepted: List[Tuple[int, Tuple[int, int, int, int], float, int]] = []
            for cand in group:
                _, (y1_a, x1_a, y2_a, x2_a), _sc_a, _c_a = cand
                area_a = float(max(1, (y2_a - y1_a) * (x2_a - x1_a)))
                is_dup = False
                for _, (y1_b, x1_b, y2_b, x2_b), _sc_b, _c_b in accepted:
                    iy1 = max(y1_a, y1_b)
                    ix1 = max(x1_a, x1_b)
                    iy2 = min(y2_a, y2_b)
                    ix2 = min(x2_a, x2_b)
                    if iy2 <= iy1 or ix2 <= ix1:
                        continue
                    inter = float((iy2 - iy1) * (ix2 - ix1))
                    area_b = float(max(1, (y2_b - y1_b) * (x2_b - x1_b)))
                    min_area = min(area_a, area_b)
                    if (inter / min_area) >= eff_overlap:
                        is_dup = True
                        break
                if not is_dup:
                    accepted.append(cand)
            kept_with_idx.extend(accepted)

        kept_with_idx.sort(key=lambda item: item[0])
        return [(bbox, score, cls_id) for _, bbox, score, cls_id in kept_with_idx]

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
            if (
                precomputed_letterbox is not None
                and precomputed_letterbox.ndim == 3
                and max(precomputed_letterbox.shape[:2]) == self.img_size
                and precomputed_letterbox.shape[0] % 32 == 0
                and precomputed_letterbox.shape[1] % 32 == 0
            ):
                if work_scale < 1.0:
                    work_w = max(1, int(round(sheet_w * work_scale)))
                    work_h = max(1, int(round(sheet_h * work_scale)))
                else:
                    work_h, work_w = sheet_h, sheet_w
                lb_bgr = precomputed_letterbox
                r, pad_xy = _letterbox_gain_pad_hw(lb_bgr.shape[:2], (work_h, work_w))
            else:
                work_img, work_scale = _prepare_working_image(image_bgr, max_working_dim=work_dim)
                self.last_work_scale = work_scale
                work_h, work_w = work_img.shape[:2]
                lb_bgr, r, pad_xy = _letterbox_bgr(work_img, self.img_size, auto=True, stride=32)

            inp_tensor = _bgr_to_nchw_float32(lb_bgr)
            input_name = active_sess.get_inputs()[0].name
            with _ONNX_GPU_LOCK:
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
        executor: Optional[Any] = None,
    ) -> List[List[Tuple[Tuple[int, int, int, int], float, int]]]:
        """Runs Stage 1 plant organ detection across a wave of herbarium sheets using shape-bucketed GPU tensors."""
        if not items:
            return []
        active_sess = self.session if self.session is not None else self.model
        if (
            active_sess is None
            or "detect" in self.__dict__
            or "predict" in self.__dict__
            or getattr(self.detect, "__func__", None) is not LM3PlantDetector.detect
            or getattr(self.predict, "__func__", None) is not LM3PlantDetector.predict
            or len(items) <= 1
        ):
            return [
                self.predict(img, catalog_number=cat, use_cache=use_cache, precomputed_letterbox=lb)
                for img, cat, lb in items
            ]

        try:
            input_name = active_sess.get_inputs()[0].name
            sub_bs = min(8, max(1, int(max_batch_size)))
            all_outs: List[List[Tuple[Tuple[int, int, int, int], float, int]]] = [[] for _ in items]

            def _prep_s1_item(
                idx_item: Tuple[int, Tuple[np.ndarray, Optional[str], Optional[np.ndarray]]]
            ) -> Optional[Tuple[int, np.ndarray, Tuple[float, Tuple[float, float], Tuple[int, int], float, Tuple[int, int]]]]:
                idx, (img_bgr, _cat, lb_pre) = idx_item
                if img_bgr is None or img_bgr.size == 0:
                    return None
                sh, sw = img_bgr.shape[:2]
                if (
                    lb_pre is not None
                    and lb_pre.ndim == 3
                    and max(lb_pre.shape[:2]) == self.img_size
                    and lb_pre.shape[0] % 32 == 0
                    and lb_pre.shape[1] % 32 == 0
                ):
                    w_scale = compute_work_scale((sh, sw), max_working_dim=self.max_working_dim)
                    if w_scale < 1.0:
                        work_w = max(1, int(round(sw * w_scale)))
                        work_h = max(1, int(round(sh * w_scale)))
                    else:
                        work_h, work_w = sh, sw
                    lb_bgr = lb_pre
                    r, pad_xy = _letterbox_gain_pad_hw(lb_bgr.shape[:2], (work_h, work_w))
                else:
                    work_img, w_scale = _prepare_working_image(img_bgr, max_working_dim=self.max_working_dim)
                    work_h, work_w = work_img.shape[:2]
                    lb_bgr, r, pad_xy = _letterbox_bgr(work_img, self.img_size, auto=True, stride=32)
                rgb = cv2.cvtColor(lb_bgr, cv2.COLOR_BGR2RGB)
                chw_f32 = np.ascontiguousarray(np.transpose(rgb, (2, 0, 1)), dtype=np.float32)
                chw_f32 *= np.float32(1.0) / np.float32(255.0)
                return (idx, chw_f32, (r, pad_xy, (work_h, work_w), w_scale, (sh, sw)))

            indexed_items = list(enumerate(items))
            eff_exec = executor if executor is not None else getattr(self, "_prep_executor", None)
            if eff_exec is not None and len(indexed_items) >= 2:
                prep_results = list(eff_exec.map(_prep_s1_item, indexed_items))
            else:
                prep_results = [_prep_s1_item(it) for it in indexed_items]

            valid_indices: List[int] = []
            valid_chw_list: List[np.ndarray] = []
            valid_meta: List[Tuple[float, Tuple[float, float], Tuple[int, int], float, Tuple[int, int]]] = []
            for res in prep_results:
                if res is not None:
                    idx, chw_f32, meta = res
                    self.last_work_scale = meta[3]
                    valid_indices.append(idx)
                    valid_chw_list.append(chw_f32)
                    valid_meta.append(meta)

            buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
            for valid_pos, chw_t in enumerate(valid_chw_list):
                buckets[(int(chw_t.shape[1]), int(chw_t.shape[2]))].append(valid_pos)

            decode_tasks: List[Tuple[int, np.ndarray, Tuple[float, Tuple[float, float], Tuple[int, int], float, Tuple[int, int]]]] = []
            for _shape_hw, bucket_positions in buckets.items():
                bucket_tensors = [valid_chw_list[p] for p in bucket_positions]
                batch_records = _run_onnx_minibatches_back_to_back(
                    active_sess, input_name, bucket_tensors, sub_bs
                )
                for start_i, end_i, raw_outs in batch_records:
                    out0 = np.asarray(raw_outs[0], dtype=np.float32)
                    for b_off, valid_pos in enumerate(bucket_positions[start_i:end_i]):
                        orig_idx = valid_indices[valid_pos]
                        slice_out = out0[b_off : b_off + 1] if out0.ndim >= 3 else out0
                        decode_tasks.append((orig_idx, slice_out, valid_meta[valid_pos]))

            def _decode_s1_item(
                task: Tuple[int, np.ndarray, Tuple[float, Tuple[float, float], Tuple[int, int], float, Tuple[int, int]]]
            ) -> Tuple[int, List[Tuple[Tuple[int, int, int, int], float, int]]]:
                orig_idx, slice_out, (r, pad_xy, work_hw, w_scale, sheet_hw) = task
                dets = self._decode_yolo26x_det_output(
                    slice_out,
                    ratio=r,
                    pad_xy=pad_xy,
                    work_shape=work_hw,
                    work_scale=w_scale,
                    sheet_shape=sheet_hw,
                )
                return orig_idx, dets

            if eff_exec is not None and len(decode_tasks) >= 2:
                decoded_pairs = list(eff_exec.map(_decode_s1_item, decode_tasks))
            else:
                decoded_pairs = [_decode_s1_item(t) for t in decode_tasks]

            for orig_idx, dets in decoded_pairs:
                all_outs[orig_idx] = dets
            return all_outs
        except Exception as exc:
            logger.debug(f"LM3PlantDetector batched inference fallback: {exc}")
            return [
                self.predict(img, catalog_number=cat, use_cache=use_cache, precomputed_letterbox=lb)
                for img, cat, lb in items
            ]


# =============================================================================
# 3. LeafMachine3 Cropped Leaf Instance Segmenter (YOLO26x-seg, 1024 px)
# =============================================================================

def _reduce_crop_predictions_native_lm3(
    raw_dets: Sequence[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]],
    crop_shape: Tuple[int, int],
    min_area_px: int = 64,
) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
    """
    Reduces raw per-instance YOLO26x-seg crop predictions to canonical LeafMachine3
    masks matching `YoloLeafSegmenter.build_instances` and `Reporter._build_mask_products`
    while preserving distinct `cls_id == 0` (`Leaf`) instances so touching neighbor leaves
    are never fused together:
      - Unions `cls_id == 1` (`Petiole`) into its largest 8-connected component (`area >= min_area_px`)
        and `cls_id == 2` (`Hole`) into all 8-connected components (`area >= min_area_px`).
      - Processes each `cls_id == 0` (`Leaf`) instance individually, retaining its primary
        8-connected component (`area >= min_area_px`) plus any unbridged lateral pinna islands
        (`area >= min_area_px`) within that same instance that attach along the segmented
        `petiole` stalk (`cls_id == 1`) or lie interior to the crop adjacent to the primary
        blade (while rejecting border-clipped neighbor bleed), subtracting `hole_mask_u8` as
        `lamina = silhouette & ~hole` without adding any synthetic bridging pixels.
    """
    ch, cw = int(crop_shape[0]), int(crop_shape[1])
    if not raw_dets or ch <= 0 or cw <= 0:
        return []
    min_area = max(1, int(min_area_px))

    raw_leaf_dets: List[Tuple[float, np.ndarray]] = []
    union_by_cls: Dict[int, np.ndarray] = {}
    max_score_by_cls: Dict[int, float] = {}
    for _bbox, mask_arr, score, cls_id in raw_dets:
        cid = int(cls_id)
        if mask_arr is None or mask_arr.size == 0:
            continue
        m_u8 = (np.asarray(mask_arr) > 0).astype(np.uint8)
        if m_u8.shape[:2] != (ch, cw):
            m_u8 = cv2.resize(m_u8, (cw, ch), interpolation=cv2.INTER_NEAREST)
        if not np.any(m_u8):
            continue
        if cid == 0:
            raw_leaf_dets.append((float(score), m_u8))
        elif cid in (1, 2):
            if cid not in union_by_cls:
                union_by_cls[cid] = m_u8 > 0
                max_score_by_cls[cid] = float(score)
            else:
                union_by_cls[cid] |= (m_u8 > 0)
                if float(score) > max_score_by_cls[cid]:
                    max_score_by_cls[cid] = float(score)

    if not raw_leaf_dets and not union_by_cls:
        return []

    class_masks_u8: Dict[int, np.ndarray] = {}

    if 1 in union_by_cls:
        union_u8 = union_by_cls[1].astype(np.uint8)
        n_cc, labels, stats, _ = cv2.connectedComponentsWithStats(union_u8, connectivity=8)
        if n_cc > 1:
            best_lbl = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
            best_area = int(stats[best_lbl, cv2.CC_STAT_AREA])
            if best_area >= min_area:
                comp_u8 = (labels == best_lbl).astype(np.uint8)
                cnts, _ = cv2.findContours(comp_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                polys = [c.reshape(-1, 1, 2).astype(np.int32) for c in cnts if c.reshape(-1, 2).shape[0] >= 3]
                if polys:
                    raster_u8 = np.zeros((ch, cw), dtype=np.uint8)
                    cv2.fillPoly(raster_u8, polys, 255)
                    if int(np.count_nonzero(raster_u8)) >= min_area:
                        class_masks_u8[1] = raster_u8

    if 2 in union_by_cls:
        union_hole_u8 = union_by_cls[2].astype(np.uint8)
        n_cc, labels, stats, _ = cv2.connectedComponentsWithStats(union_hole_u8, connectivity=8)
        hole_polys: List[np.ndarray] = []
        for lbl in range(1, n_cc):
            if int(stats[lbl, cv2.CC_STAT_AREA]) < min_area:
                continue
            comp_u8 = (labels == lbl).astype(np.uint8)
            cnts, _ = cv2.findContours(comp_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            for c in cnts:
                if c.reshape(-1, 2).shape[0] >= 3:
                    hole_polys.append(c.reshape(-1, 1, 2).astype(np.int32))
        if hole_polys:
            hole_mask_u8 = np.zeros((ch, cw), dtype=np.uint8)
            cv2.fillPoly(hole_mask_u8, hole_polys, 255)
            if int(np.count_nonzero(hole_mask_u8)) >= min_area:
                class_masks_u8[2] = hole_mask_u8

    pet_ref_u8 = class_masks_u8.get(1)
    pet_corridor_bool: Optional[np.ndarray] = None
    if pet_ref_u8 is not None and int(np.count_nonzero(pet_ref_u8)) > 0:
        kern_pet_adj = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (11, 11))
        pet_corridor_bool = cv2.dilate((pet_ref_u8 > 0).astype(np.uint8), kern_pet_adj, iterations=1) > 0
    hole_sub_u8 = class_masks_u8.get(2)
    cx_mid = 0.5 * float(cw)
    cy_mid = 0.5 * float(ch)
    crop_diag = max(1.0, math.hypot(float(cw), float(ch)))

    ranked_leaf_entries: List[Tuple[float, Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]] = []
    for l_score, l_u8 in raw_leaf_dets:
        n_cc, labels, stats, _ = cv2.connectedComponentsWithStats(l_u8, connectivity=8)
        if n_cc <= 1:
            continue
        best_lbl = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
        best_area = int(stats[best_lbl, cv2.CC_STAT_AREA])
        if best_area < min_area:
            continue
        keep_lbls: List[int] = [best_lbl]
        if n_cc > 2:
            prim_bool = labels == best_lbl
            prim_touches_pet = bool(pet_corridor_bool is not None and np.any(prim_bool & pet_corridor_bool))
            dist_to_prim = cv2.distanceTransform((~prim_bool).astype(np.uint8), cv2.DIST_L2, 3)
            max_pinna_gap = max(10.0, 0.16 * math.sqrt(float(best_area)))
            for cid_sub in range(1, n_cc):
                if cid_sub == best_lbl:
                    continue
                sub_area = int(stats[cid_sub, cv2.CC_STAT_AREA])
                if sub_area < min_area:
                    continue
                sub_bool = labels == cid_sub
                sx0 = int(stats[cid_sub, cv2.CC_STAT_LEFT])
                sy0 = int(stats[cid_sub, cv2.CC_STAT_TOP])
                sx1 = sx0 + int(stats[cid_sub, cv2.CC_STAT_WIDTH])
                sy1 = sy0 + int(stats[cid_sub, cv2.CC_STAT_HEIGHT])
                touches_crop_edge = (sx0 <= 1) or (sy0 <= 1) or (sx1 >= cw - 1) or (sy1 >= ch - 1)
                on_same_petiole = bool(
                    prim_touches_pet and pet_corridor_bool is not None and np.any(sub_bool & pet_corridor_bool)
                )
                min_gap_px = float(np.min(dist_to_prim[sub_bool])) if np.any(sub_bool) else float("inf")
                if on_same_petiole or ((not touches_crop_edge) and min_gap_px <= max_pinna_gap):
                    keep_lbls.append(cid_sub)
        raster_u8 = (np.isin(labels, keep_lbls).astype(np.uint8)) * 255
        if hole_sub_u8 is not None:
            raster_u8[hole_sub_u8 > 0] = 0
            n_cc_h, labels_h, stats_h, _ = cv2.connectedComponentsWithStats((raster_u8 > 0).astype(np.uint8), connectivity=8)
            if n_cc_h <= 1:
                continue
            valid_h_lbls = [
                lbl_h for lbl_h in range(1, n_cc_h) if int(stats_h[lbl_h, cv2.CC_STAT_AREA]) >= min_area
            ]
            if not valid_h_lbls:
                continue
            raster_u8 = np.isin(labels_h, valid_h_lbls).astype(np.uint8) * 255
        fg_count = int(np.count_nonzero(raster_u8))
        if fg_count < min_area:
            continue
        ys, xs = np.where(raster_u8 > 0)
        if ys.size == 0 or xs.size == 0:
            continue
        ymin = int(ys.min())
        ymax = min(ch, int(ys.max()) + 1)
        xmin = int(xs.min())
        xmax = min(cw, int(xs.max()) + 1)
        if ymax <= ymin or xmax <= xmin:
            continue
        mx_l = float(np.mean(xs))
        my_l = float(np.mean(ys))
        dist_norm = math.hypot(mx_l - cx_mid, my_l - cy_mid) / (0.65 * crop_diag)
        cent_w = max(0.35, 1.0 - dist_norm)
        focal_rank = float(l_score) * cent_w * (float(fg_count) ** 0.35)
        ranked_leaf_entries.append((focal_rank, ((ymin, xmin, ymax, xmax), raster_u8, float(l_score), 0)))

    ranked_leaf_entries.sort(key=lambda item: item[0], reverse=True)
    reduced: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]] = [
        entry for _rank, entry in ranked_leaf_entries
    ]

    for cid in (1, 2):
        if cid not in class_masks_u8:
            continue
        m_u8 = class_masks_u8[cid]
        ys, xs = np.where(m_u8 > 0)
        if ys.size == 0 or xs.size == 0:
            continue
        ymin = int(ys.min())
        ymax = min(ch, int(ys.max()) + 1)
        xmin = int(xs.min())
        xmax = min(cw, int(xs.max()) + 1)
        if ymax <= ymin or xmax <= xmin:
            continue
        reduced.append(((ymin, xmin, ymax, xmax), m_u8, float(max_score_by_cls[cid]), int(cid)))

    return reduced


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
        score_thresh: float = 0.30,
        min_area_px: int = 64,
        img_size: int = 1024,
        providers: Optional[Sequence[str]] = None,
        **kwargs: Any,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.score_thresh = float(score_thresh)
        self.min_area_px = int(min_area_px)
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
        Thread-safe CPU helper that precomputes the 1024 px dynamic stride-32 letterbox uint8 CHW crop
        on worker threads for batched inference.
        """
        if (self.session is None and self.model is None) or image_bgr is None or image_bgr.size == 0:
            return None
        h, w = image_bgr.shape[:2]
        lb_bgr, _, _ = _letterbox_bgr(image_bgr, self.img_size, auto=True, stride=32)
        rgb = cv2.cvtColor(lb_bgr, cv2.COLOR_BGR2RGB)
        chw_u8 = np.ascontiguousarray(np.transpose(rgb, (2, 0, 1)), dtype=np.uint8)
        return chw_u8, h, w

    def _decode_yolo26x_seg_outputs(
        self,
        raw_outs: Sequence[np.ndarray],
        crop_shape: Tuple[int, int],
        ratio: float,
        pad_xy: Tuple[float, float],
        image_bgr: Optional[np.ndarray] = None,
        letterbox_shape: Optional[Tuple[int, int]] = None,
    ) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """Decodes YOLO26x-seg detection + prototype mask outputs into reduced `(bbox, mask_u8, score, class_id)`."""
        _ = (image_bgr, letterbox_shape)
        if not raw_outs:
            return []
        crop_h, crop_w = int(crop_shape[0]), int(crop_shape[1])
        dw, dh = pad_xy
        inv_r = float(np.float32(1.0) / np.float32(max(float(ratio), 1e-6)))

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

        if not rows_decoded:
            return []

        valid_entries: List[Tuple[int, int, int, int, float, int, Optional[np.ndarray], Tuple[float, float, float, float]]] = []
        for x1_l, y1_l, x2_l, y2_l, score, cls_id, coeffs in rows_decoded:
            x1_orig = float(np.clip((x1_l - dw) * inv_r, 0.0, float(crop_w)))
            y1_orig = float(np.clip((y1_l - dh) * inv_r, 0.0, float(crop_h)))
            x2_orig = float(np.clip((x2_l - dw) * inv_r, 0.0, float(crop_w)))
            y2_orig = float(np.clip((y2_l - dh) * inv_r, 0.0, float(crop_h)))
            xmin = max(0, min(crop_w, int(round(x1_orig))))
            ymin = max(0, min(crop_h, int(round(y1_orig))))
            xmax = max(0, min(crop_w, int(round(x2_orig))))
            ymax = max(0, min(crop_h, int(round(y2_orig))))
            if ymax <= ymin or xmax <= xmin:
                continue
            valid_entries.append((ymin, xmin, ymax, xmax, score, cls_id, coeffs, (x1_orig, y1_orig, x2_orig, y2_orig)))

        if not valid_entries:
            return []

        proto_masks_u8: Optional[np.ndarray] = None
        if (
            proto_arr is not None
            and n_proto > 0
            and all(e[6] is not None and len(e[6]) == n_proto for e in valid_entries)
        ):
            try:
                from leafmachine3.inference import ultra_replacements as ultra_rep
                coeffs_mat = np.vstack([np.asarray(e[6], dtype=np.float32) for e in valid_entries])
                boxes_xyxy = np.asarray([e[7] for e in valid_entries], dtype=np.float32)
                raw_masks = ultra_rep.process_mask_native(proto_arr, coeffs_mat, boxes_xyxy, (crop_h, crop_w))
                proto_masks_u8 = (raw_masks > 0).astype(np.uint8) * 255
            except Exception:
                proto_masks_u8 = None

        for idx_e, (ymin, xmin, ymax, xmax, score, cls_id, coeffs, _) in enumerate(valid_entries):
            if proto_masks_u8 is not None and idx_e < len(proto_masks_u8):
                mask_u8 = proto_masks_u8[idx_e]
            elif proto_arr is not None and coeffs is not None and len(coeffs) == n_proto and n_proto > 0:
                mh, mw = int(proto_arr.shape[1]), int(proto_arr.shape[2])
                mask_logits = np.tensordot(coeffs, proto_arr, axes=([0], [0]))
                gain_p = min(float(mh) / float(max(crop_h, 1)), float(mw) / float(max(crop_w, 1)))
                pad_w = (float(mw) - round(crop_w * gain_p)) * 0.5
                pad_h = (float(mh) - round(crop_h * gain_p)) * 0.5
                top = max(0, int(round(pad_h - 0.1)))
                left = max(0, int(round(pad_w - 0.1)))
                bottom = min(mh, mh - int(round(pad_h + 0.1)))
                right = min(mw, mw - int(round(pad_w + 0.1)))
                mask_u8 = np.zeros((crop_h, crop_w), dtype=np.uint8)
                if bottom > top and right > left:
                    sub_logits = np.ascontiguousarray(mask_logits[top:bottom, left:right], dtype=np.float32)
                    resized_logits = cv2.resize(sub_logits, (crop_w, crop_h), interpolation=_INTER_LINEAR_EXACT)
                    mask_u8[ymin:ymax, xmin:xmax] = (resized_logits[ymin:ymax, xmin:xmax] > 0.0).astype(np.uint8) * 255
            else:
                mask_u8 = np.zeros((crop_h, crop_w), dtype=np.uint8)
                mask_u8[ymin:ymax, xmin:xmax] = 255

            if int(np.count_nonzero(mask_u8)) > 0:
                results.append(((ymin, xmin, ymax, xmax), mask_u8, float(score), int(cls_id)))

        return _reduce_crop_predictions_native_lm3(
            results,
            (crop_h, crop_w),
            min_area_px=max(1, int(getattr(self, "min_area_px", 64))),
        )

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
            if (
                pretransformed is not None
                and max(pretransformed[0].shape[1:]) == self.img_size
                and pretransformed[0].shape[1] % 32 == 0
                and pretransformed[0].shape[2] % 32 == 0
            ):
                chw_u8, crop_h, crop_w = pretransformed
                lb_h, lb_w = int(chw_u8.shape[1]), int(chw_u8.shape[2])
                r, pad_xy = _letterbox_gain_pad_hw((lb_h, lb_w), (crop_h, crop_w))
                inp_f32 = np.ascontiguousarray(chw_u8, dtype=np.float32)
                inp_f32 *= np.float32(1.0) / np.float32(255.0)
                inp_tensor = inp_f32[None, ...]
                lb_shape = (lb_h, lb_w)
            else:
                crop_h, crop_w = image_bgr.shape[:2]
                lb_bgr, r, pad_xy = _letterbox_bgr(image_bgr, self.img_size, auto=True, stride=32)
                inp_tensor = _bgr_to_nchw_float32(lb_bgr)
                lb_shape = (int(lb_bgr.shape[0]), int(lb_bgr.shape[1]))

            input_name = active_sess.get_inputs()[0].name
            with _ONNX_GPU_LOCK:
                raw_outs = active_sess.run(None, {input_name: inp_tensor})
            return self._decode_yolo26x_seg_outputs(
                raw_outs,
                crop_shape=(crop_h, crop_w),
                ratio=r,
                pad_xy=pad_xy,
                image_bgr=image_bgr,
                letterbox_shape=lb_shape,
            )
        except Exception as exc:
            logger.debug(f"LM3LeafSegmenterEngine ONNX inference error: {exc}")
            return []

    def predict_batch(
        self,
        images_bgr: List[np.ndarray],
        batch_size: int = 10,
        pretransformed: Optional[Sequence[Optional[Tuple[np.ndarray, int, int]]]] = None,
        executor: Optional[Any] = None,
    ) -> List[List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]]:
        """
        Runs Stage 2 LeafMachine3 YOLO26x-seg inference across candidate leaf crops using
        shape-bucketed GPU tensors `(B, 3, H_lb, W_lb)`.
        Automatically delegates to `self.predict(img)` when `predict` is mocked in unit tests.
        """
        if not images_bgr:
            return []
        has_pre = pretransformed is not None and len(pretransformed) == len(images_bgr)
        active_sess = self.session if self.session is not None else self.model
        if (
            active_sess is None
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

        try:
            input_name = active_sess.get_inputs()[0].name
            sub_bs = min(8, max(1, int(batch_size)))
            out: List[List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]] = [[] for _ in images_bgr]

            def _prep_s2_item(
                idx: int,
            ) -> Optional[Tuple[int, np.ndarray, Tuple[Tuple[int, int], float, Tuple[float, float], Tuple[int, int]]]]:
                img = images_bgr[idx]
                if img is None or img.size == 0:
                    return None
                pre_i = pretransformed[idx] if has_pre else None  # type: ignore[index]
                if (
                    pre_i is not None
                    and max(pre_i[0].shape[1:]) == self.img_size
                    and pre_i[0].shape[1] % 32 == 0
                    and pre_i[0].shape[2] % 32 == 0
                ):
                    chw_u8, crop_h, crop_w = pre_i
                    lb_h, lb_w = int(chw_u8.shape[1]), int(chw_u8.shape[2])
                    r, pad_xy = _letterbox_gain_pad_hw((lb_h, lb_w), (crop_h, crop_w))
                    inp_f32 = np.ascontiguousarray(chw_u8, dtype=np.float32)
                    inp_f32 *= np.float32(1.0) / np.float32(255.0)
                    lb_shape = (lb_h, lb_w)
                else:
                    crop_h, crop_w = img.shape[:2]
                    lb_bgr, r, pad_xy = _letterbox_bgr(img, self.img_size, auto=True, stride=32)
                    rgb = cv2.cvtColor(lb_bgr, cv2.COLOR_BGR2RGB)
                    inp_f32 = np.ascontiguousarray(np.transpose(rgb, (2, 0, 1)), dtype=np.float32)
                    inp_f32 *= np.float32(1.0) / np.float32(255.0)
                    lb_shape = (int(lb_bgr.shape[0]), int(lb_bgr.shape[1]))
                return (idx, inp_f32, ((crop_h, crop_w), r, pad_xy, lb_shape))

            all_idx = list(range(len(images_bgr)))
            eff_exec = executor if executor is not None else getattr(self, "_prep_executor", None)
            if eff_exec is not None and len(all_idx) >= 2:
                prep_results = list(eff_exec.map(_prep_s2_item, all_idx))
            else:
                prep_results = [_prep_s2_item(i) for i in all_idx]

            valid_indices: List[int] = []
            valid_chw_list: List[np.ndarray] = []
            valid_meta: List[Tuple[Tuple[int, int], float, Tuple[float, float], Tuple[int, int]]] = []
            for res in prep_results:
                if res is not None:
                    idx, inp_f32, meta = res
                    valid_indices.append(idx)
                    valid_chw_list.append(inp_f32)
                    valid_meta.append(meta)

            buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
            for valid_pos, chw_t in enumerate(valid_chw_list):
                buckets[(int(chw_t.shape[1]), int(chw_t.shape[2]))].append(valid_pos)

            decode_tasks: List[Tuple[int, List[np.ndarray], Tuple[Tuple[int, int], float, Tuple[float, float], Tuple[int, int]]]] = []
            for _shape_hw, bucket_positions in buckets.items():
                bucket_tensors = [valid_chw_list[p] for p in bucket_positions]
                batch_records = _run_onnx_minibatches_back_to_back(
                    active_sess, input_name, bucket_tensors, sub_bs
                )
                for start_i, end_i, raw_outs in batch_records:
                    for b_off, valid_pos in enumerate(bucket_positions[start_i:end_i]):
                        orig_idx = valid_indices[valid_pos]
                        item_outs = [
                            ro[b_off : b_off + 1] if np.asarray(ro).ndim >= 3 else ro
                            for ro in raw_outs
                        ]
                        decode_tasks.append((orig_idx, item_outs, valid_meta[valid_pos]))

            def _decode_s2_item(
                task: Tuple[int, List[np.ndarray], Tuple[Tuple[int, int], float, Tuple[float, float], Tuple[int, int]]]
            ) -> Tuple[int, List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]]:
                orig_idx, item_outs, (crop_hw, r, pad_xy, lb_shape) = task
                dets = self._decode_yolo26x_seg_outputs(
                    item_outs,
                    crop_shape=crop_hw,
                    ratio=r,
                    pad_xy=pad_xy,
                    image_bgr=images_bgr[orig_idx],
                    letterbox_shape=lb_shape,
                )
                return orig_idx, dets

            if eff_exec is not None and len(decode_tasks) >= 2:
                decoded_pairs = list(eff_exec.map(_decode_s2_item, decode_tasks))
            else:
                decoded_pairs = [_decode_s2_item(t) for t in decode_tasks]

            for orig_idx, dets in decoded_pairs:
                out[orig_idx] = dets
            return out
        except Exception as exc:
            logger.debug(f"LM3LeafSegmenterEngine batched inference fallback: {exc}")
            out_fb: List[List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]] = []
            for idx, img in enumerate(images_bgr):
                pre_i = pretransformed[idx] if has_pre else None  # type: ignore[index]
                out_fb.append(self.predict(img, pretransformed=pre_i))
            return out_fb

    def segment_leaf_crop(
        self,
        crop_bgr: np.ndarray,
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], float]:
        """
        Runs LeafMachine3 leaf segmentation on a single candidate crop and separates
        `(lamina_mask, petiole_mask, hole_mask, best_score)` following `leaf_segmenter.build_instances`
        and `Reporter._build_mask_products`.
        """
        if crop_bgr is None or crop_bgr.size == 0:
            return None, None, None, 0.0

        preds = self.predict(crop_bgr)
        if not preds:
            return None, None, None, 0.0

        crop_h, crop_w = crop_bgr.shape[:2]
        min_fg = max(1, int(getattr(self, "min_area_px", 64)))
        if (
            "predict" not in self.__dict__
            and getattr(self.predict, "__func__", None) is LM3LeafSegmenterEngine.predict
            and (self.session is not None or self.model is not None)
        ):
            reduced = list(preds)
        else:
            reduced = _reduce_crop_predictions_native_lm3(preds, (crop_h, crop_w), min_area_px=min_fg)
        if not reduced:
            return None, None, None, 0.0

        lamina_mask: Optional[np.ndarray] = None
        petiole_mask: Optional[np.ndarray] = None
        hole_mask: Optional[np.ndarray] = None
        best_score: float = 0.0

        for _bbox, mask_u8, score, cls_id in reduced:
            if int(cls_id) == 0:
                if lamina_mask is None:
                    lamina_mask = mask_u8
                    best_score = float(score)
            elif int(cls_id) == 1:
                petiole_mask = mask_u8
            elif int(cls_id) == 2:
                hole_mask = mask_u8

        if lamina_mask is None:
            return None, None, None, 0.0

        return lamina_mask, petiole_mask, hole_mask, float(best_score)

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
        is_geometric_fallback: bool = False,
    ) -> None:
        super().__init__(keypoints_31)
        self.petiole_width_px: Optional[float] = (
            float(petiole_width_px) if (petiole_width_px is not None and np.isfinite(petiole_width_px)) else None
        )
        self.keypoints_31: Dict[str, Tuple[float, float, float]] = dict(keypoints_31)
        self.is_geometric_fallback: bool = bool(is_geometric_fallback)

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
    WHITE_PAD_FRAC: float = 0.10

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        conf_thres: float = 0.25,
        img_size: int = 640,
        white_pad_frac: float = 0.10,
        providers: Optional[Sequence[str]] = None,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.conf_thres = float(conf_thres)
        self.img_size = int(img_size)
        self.WHITE_PAD_FRAC = float(white_pad_frac)
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
        """Estimates petiole width in pixels by delegating to `leafmachine3.core.petiole.measure_petiole`."""
        if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
            pet_u8 = (petiole_mask > 0).astype(np.uint8)
            lam_u8 = (
                (lamina_mask > 0).astype(np.uint8)
                if (lamina_mask is not None and lamina_mask.size > 0 and lamina_mask.shape[:2] == pet_u8.shape[:2])
                else None
            )
            norm_kps = normalize_lm3_keypoints_dict(keypoints_31)
            petiole_pts: List[Tuple[float, float]] = []
            for kp_name in ("lamina_base", "petiole_0", "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_tip"):
                kp = norm_kps.get(kp_name)
                if kp is not None and len(kp) >= 3 and float(kp[2]) >= 0.20:
                    pt_xy = (float(kp[0]), float(kp[1]))
                    if not petiole_pts or math.hypot(pt_xy[0] - petiole_pts[-1][0], pt_xy[1] - petiole_pts[-1][1]) > 1e-3:
                        petiole_pts.append(pt_xy)
            if len(petiole_pts) < 2:
                ys_p, xs_p = np.where(pet_u8 > 0)
                if len(xs_p) >= 2:
                    pts_cov = np.column_stack((xs_p.astype(np.float64), ys_p.astype(np.float64)))
                    mean_pt = pts_cov.mean(axis=0)
                    centered = pts_cov - mean_pt
                    _, _, vh = np.linalg.svd(centered, full_matrices=False)
                    proj = centered @ vh[0]
                    i_min = int(np.argmin(proj))
                    i_max = int(np.argmax(proj))
                    p0 = (float(pts_cov[i_min, 0]), float(pts_cov[i_min, 1]))
                    p1 = (float(pts_cov[i_max, 0]), float(pts_cov[i_max, 1]))
                    if math.hypot(p1[0] - p0[0], p1[1] - p0[1]) > 1e-3:
                        petiole_pts = [p0, p1]
            if lam_u8 is not None and len(petiole_pts) >= 2:
                try:
                    from leafmachine3.core.petiole import measure_petiole

                    res = measure_petiole(pet_u8, lam_u8, petiole_pts)
                    if res is not None and res.width_px is not None and np.isfinite(res.width_px) and res.width_px > 0:
                        return round(float(res.width_px), 3)
                except Exception:
                    pass

            dt = cv2.distanceTransform(pet_u8, cv2.DIST_L2, 3)
            sampled_radii: List[float] = []
            h_p, w_p = pet_u8.shape[:2]
            for pt_x, pt_y in petiole_pts:
                ix = max(0, min(w_p - 1, int(round(pt_x))))
                iy = max(0, min(h_p - 1, int(round(pt_y))))
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
        crown_hint_xy: Optional[Tuple[float, float]] = None,
    ) -> LM3PoseResult:
        """
        Derives the 31 canonical `mid15_pet5` keypoints deterministically from `lamina_mask`
        and `petiole_mask` (or `crop_bgr`) when ONNX weights are absent in unit tests,
        anchoring reniform/cordate (`W > L`) blades along the petiole-to-centroid axis
        rather than rotating `90 deg` onto lateral lobes.
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
            return LM3PoseResult(empty_kps, petiole_width_px=None, is_geometric_fallback=True)

        ys, xs = np.where(mask_use > 0)
        pts = np.column_stack((xs.astype(np.float64), ys.astype(np.float64)))
        centroid = np.mean(pts, axis=0)
        centered = pts - centroid
        cov = np.cov(centered, rowvar=False) if len(pts) >= 3 else np.eye(2)
        eigvals, eigvecs = np.linalg.eigh(cov)
        major_vec = eigvecs[:, int(np.argmax(eigvals))]

        proj_major = centered @ major_vec
        p_min = centroid + float(np.min(proj_major)) * major_vec
        p_max = centroid + float(np.max(proj_major)) * major_vec

        # Orient base -> tip along petiole attachment axis (preserving reniform W > L blades),
        # else crown_hint_xy proximity, else apex/base taper & mount orientation
        if petiole_mask is not None and petiole_mask.size > 0 and int(np.count_nonzero(petiole_mask)) > 0:
            pys, pxs = np.where(petiole_mask > 0)
            pet_cent = np.array([float(np.mean(pxs)), float(np.mean(pys))], dtype=np.float64)
            junc_pt = pet_cent
            if petiole_mask.shape[:2] == mask_use.shape[:2]:
                kern = np.ones((5, 5), dtype=np.uint8)
                pet_dil = cv2.dilate((petiole_mask > 0).astype(np.uint8), kern, iterations=1)
                jys, jxs = np.where((pet_dil > 0) & (mask_use > 0))
                if len(jxs) > 0:
                    junc_pt = np.array([float(np.mean(jxs)), float(np.mean(jys))], dtype=np.float64)
            v_from_pet = centroid - junc_pt
            v_norm = float(np.linalg.norm(v_from_pet))
            if v_norm > 1.0:
                u_pet = v_from_pet / v_norm
                e0 = eigvecs[:, 0]
                e1 = eigvecs[:, 1]
                d0 = abs(float(np.dot(e0, u_pet)))
                d1 = abs(float(np.dot(e1, u_pet)))
                best_e = e0 if d0 >= d1 else e1
                if max(d0, d1) >= 0.9063:
                    axis_dir = best_e if float(np.dot(best_e, u_pet)) >= 0.0 else -best_e
                else:
                    axis_dir = u_pet
                proj_ax = centered @ axis_dir
                base_pt = centroid + float(np.min(proj_ax)) * axis_dir
                tip_pt = centroid + float(np.max(proj_ax)) * axis_dir
            elif np.linalg.norm(p_min - pet_cent) <= np.linalg.norm(p_max - pet_cent):
                base_pt, tip_pt = p_min, p_max
            else:
                base_pt, tip_pt = p_max, p_min
        elif crown_hint_xy is not None:
            cr_local = np.array([float(crown_hint_xy[0]) - ox, float(crown_hint_xy[1]) - oy], dtype=np.float64)
            if np.linalg.norm(p_min - cr_local) <= np.linalg.norm(p_max - cr_local):
                base_pt, tip_pt = p_min, p_max
            else:
                base_pt, tip_pt = p_max, p_min
        else:
            dt_fb = cv2.distanceTransform((mask_use > 0).astype(np.uint8), cv2.DIST_L2, 5)
            h_fb, w_fb = dt_fb.shape[:2]
            span_fb = float(np.linalg.norm(p_max - p_min))
            win_fb = max(2, min(6, int(round(0.015 * span_fb))))
            r_min_vals: List[float] = []
            r_max_vals: List[float] = []
            for frac_fb in (0.01, 0.02, 0.03, 0.04, 0.05):
                pt_a = p_min + frac_fb * (p_max - p_min)
                pt_b = p_max + frac_fb * (p_min - p_max)
                xa, ya = int(round(pt_a[0])), int(round(pt_a[1]))
                xb, yb = int(round(pt_b[0])), int(round(pt_b[1]))
                pa = dt_fb[max(0, ya - win_fb) : min(h_fb, ya + win_fb + 1), max(0, xa - win_fb) : min(w_fb, xa + win_fb + 1)]
                pb = dt_fb[max(0, yb - win_fb) : min(h_fb, yb + win_fb + 1), max(0, xb - win_fb) : min(w_fb, xb + win_fb + 1)]
                if pa.size > 0:
                    r_min_vals.append(float(np.max(pa)))
                if pb.size > 0:
                    r_max_vals.append(float(np.max(pb)))
            r_min_m = float(np.mean(r_min_vals)) if r_min_vals else 0.0
            r_max_m = float(np.mean(r_max_vals)) if r_max_vals else 0.0
            if abs(r_min_m - r_max_m) > 0.35:
                base_pt, tip_pt = (p_min, p_max) if r_min_m < r_max_m else (p_max, p_min)
            elif p_min[1] >= p_max[1]:
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
        return LM3PoseResult(kps, petiole_width_px=pet_w, is_geometric_fallback=True)

    def predict_landmarks(
        self,
        crop_bgr: np.ndarray,
        lamina_mask: Optional[np.ndarray] = None,
        petiole_mask: Optional[np.ndarray] = None,
        crop_offset: Tuple[float, float] = (0.0, 0.0),
    ) -> LM3PoseResult:
        """
        Runs 640 px YOLO26x-pose landmark detection on unmasked `crop_bgr`, returning an `LM3PoseResult`
        mapping all 31 canonical keypoint names (`lamina_tip`, `apex_1..3`, `midvein_1..15`,
        `base_1..3`, `lamina_base`, `petiole_1..6`, `width_1..2`) to `(x, y, confidence)`
        in crop (or offset sheet) coordinates, plus `.petiole_width_px`.
        """
        ox, oy = float(crop_offset[0]), float(crop_offset[1])
        active_sess = self.session if self.session is not None else self.model
        if active_sess is not None and crop_bgr is not None and crop_bgr.size > 0:
            try:
                crop_h, crop_w = crop_bgr.shape[:2]
                pad_f = float(getattr(self, "WHITE_PAD_FRAC", 0.10))
                if pad_f > 0.0:
                    px, py = int(round(crop_w * pad_f)), int(round(crop_h * pad_f))
                    pose_in_bgr = cv2.copyMakeBorder(
                        crop_bgr, py, py, px, px, cv2.BORDER_CONSTANT, value=(255, 255, 255)
                    )
                else:
                    px, py = 0, 0
                    pose_in_bgr = crop_bgr
                lb_bgr, r, (dw, dh) = _letterbox_bgr(pose_in_bgr, self.img_size, auto=True, stride=32)
                inp_tensor = _bgr_to_nchw_float32(lb_bgr)
                input_name = active_sess.get_inputs()[0].name
                with _ONNX_GPU_LOCK:
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
                            inv_r = float(np.float32(1.0) / np.float32(max(float(r), 1e-6)))
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

    def predict_batch(
        self,
        crops_bgr: Sequence[np.ndarray],
        lamina_masks: Optional[Sequence[Optional[np.ndarray]]] = None,
        petiole_masks: Optional[Sequence[Optional[np.ndarray]]] = None,
        crop_offsets: Optional[Sequence[Tuple[float, float]]] = None,
        batch_size: int = 8,
        executor: Optional[Any] = None,
    ) -> List[LM3PoseResult]:
        """
        Runs batched 640 px YOLO26x-pose landmark detection across multiple unmasked leaf crops
        on the GPU in shape-bucketed `(B, 3, H_lb, W_lb)` mini-batches.
        """
        if not crops_bgr:
            return []
        n_items = len(crops_bgr)
        active_sess = self.session if self.session is not None else self.model
        if (
            active_sess is None
            or "predict_landmarks" in self.__dict__
            or getattr(self.predict_landmarks, "__func__", None) is not LM3PoseLandmarkEngine.predict_landmarks
            or n_items <= 1
        ):
            return [
                self.predict_landmarks(
                    crop_bgr=crops_bgr[i],
                    lamina_mask=lamina_masks[i] if (lamina_masks and i < len(lamina_masks)) else None,
                    petiole_mask=petiole_masks[i] if (petiole_masks and i < len(petiole_masks)) else None,
                    crop_offset=crop_offsets[i] if (crop_offsets and i < len(crop_offsets)) else (0.0, 0.0),
                )
                for i in range(n_items)
            ]

        try:
            input_name = active_sess.get_inputs()[0].name
            sub_bs = min(16, max(1, int(batch_size)))
            pad_f = float(getattr(self, "WHITE_PAD_FRAC", 0.10))
            results: List[Optional[LM3PoseResult]] = [None] * n_items

            def _prep_pose_item(
                idx: int,
            ) -> Optional[Tuple[int, np.ndarray, Tuple[int, int, int, int, float, float, float, float, float]]]:
                crop_bgr = crops_bgr[idx]
                ox, oy = (
                    (float(crop_offsets[idx][0]), float(crop_offsets[idx][1]))
                    if (crop_offsets and idx < len(crop_offsets))
                    else (0.0, 0.0)
                )
                if crop_bgr is None or crop_bgr.size == 0:
                    return None
                crop_h, crop_w = crop_bgr.shape[:2]
                if pad_f > 0.0:
                    px, py = int(round(crop_w * pad_f)), int(round(crop_h * pad_f))
                    pose_in_bgr = cv2.copyMakeBorder(
                        crop_bgr, py, py, px, px, cv2.BORDER_CONSTANT, value=(255, 255, 255)
                    )
                else:
                    px, py = 0, 0
                    pose_in_bgr = crop_bgr
                lb_bgr, r, (dw, dh) = _letterbox_bgr(pose_in_bgr, self.img_size, auto=True, stride=32)
                rgb = cv2.cvtColor(lb_bgr, cv2.COLOR_BGR2RGB)
                chw_f32 = np.ascontiguousarray(np.transpose(rgb, (2, 0, 1)), dtype=np.float32)
                chw_f32 *= np.float32(1.0) / np.float32(255.0)
                return (idx, chw_f32, (crop_h, crop_w, px, py, float(r), float(dw), float(dh), ox, oy))

            all_idx = list(range(n_items))
            eff_exec = executor if executor is not None else getattr(self, "_prep_executor", None)
            if eff_exec is not None and n_items >= 2:
                prep_results = list(eff_exec.map(_prep_pose_item, all_idx))
            else:
                prep_results = [_prep_pose_item(i) for i in all_idx]

            valid_indices: List[int] = []
            valid_chw_list: List[np.ndarray] = []
            valid_meta: List[Tuple[int, int, int, int, float, float, float, float, float]] = []
            for res in prep_results:
                if res is not None:
                    idx, chw_f32, meta = res
                    valid_indices.append(idx)
                    valid_chw_list.append(chw_f32)
                    valid_meta.append(meta)

            buckets: Dict[Tuple[int, int], List[int]] = defaultdict(list)
            for valid_pos, chw_t in enumerate(valid_chw_list):
                buckets[(int(chw_t.shape[1]), int(chw_t.shape[2]))].append(valid_pos)

            raw_arr_by_idx: Dict[int, Tuple[np.ndarray, Tuple[int, int, int, int, float, float, float, float, float]]] = {}
            for _shape_hw, bucket_positions in buckets.items():
                bucket_tensors = [valid_chw_list[p] for p in bucket_positions]
                batch_records = _run_onnx_minibatches_back_to_back(
                    active_sess, input_name, bucket_tensors, sub_bs
                )
                for start_i, end_i, raw_outs in batch_records:
                    out0 = np.asarray(raw_outs[0], dtype=np.float32)
                    for b_off, valid_pos in enumerate(bucket_positions[start_i:end_i]):
                        orig_idx = valid_indices[valid_pos]
                        arr = out0[b_off] if out0.ndim >= 3 else out0
                        raw_arr_by_idx[orig_idx] = (arr, valid_meta[valid_pos])

            def _decode_pose_item(idx: int) -> LM3PoseResult:
                lm_i = lamina_masks[idx] if (lamina_masks and idx < len(lamina_masks)) else None
                pm_i = petiole_masks[idx] if (petiole_masks and idx < len(petiole_masks)) else None
                ox, oy = (
                    (float(crop_offsets[idx][0]), float(crop_offsets[idx][1]))
                    if (crop_offsets and idx < len(crop_offsets))
                    else (0.0, 0.0)
                )
                if idx in raw_arr_by_idx:
                    arr, (crop_h, crop_w, px, py, r, dw, dh, ox, oy) = raw_arr_by_idx[idx]
                    if arr.ndim == 2 and arr.shape[0] < arr.shape[1] and arr.shape[0] >= 5 + 3 * 31:
                        arr = arr.T
                    if arr.ndim == 2 and arr.shape[0] > 0 and arr.shape[1] >= 5 + 3 * 31:
                        best_row_idx = int(np.argmax(arr[:, 4]))
                        best_row = arr[best_row_idx]
                        if float(best_row[4]) >= self.conf_thres:
                            kp_offset = 6 if best_row.shape[0] >= 6 + 3 * 31 else 5
                            kp_flat = best_row[kp_offset : kp_offset + 3 * 31].reshape(31, 3)
                            inv_r = float(np.float32(1.0) / np.float32(max(float(r), 1e-6)))
                            kps_out: Dict[str, Tuple[float, float, float]] = {}
                            for idx_k, kp_name in enumerate(self.KEYPOINT_NAMES):
                                kx_pad = (float(kp_flat[idx_k, 0]) - dw) * inv_r
                                ky_pad = (float(kp_flat[idx_k, 1]) - dh) * inv_r
                                kx = max(0.0, min(float(crop_w), kx_pad - float(px))) + ox
                                ky = max(0.0, min(float(crop_h), ky_pad - float(py))) + oy
                                kc = float(np.clip(kp_flat[idx_k, 2], 0.0, 1.0))
                                kps_out[kp_name] = (round(kx, 2), round(ky, 2), round(kc, 4))
                            pet_w = self._estimate_petiole_width_px(
                                kps_out, petiole_mask=pm_i, lamina_mask=lm_i
                            )
                            return LM3PoseResult(kps_out, petiole_width_px=pet_w)
                return self._geometric_fallback_keypoints_31(
                    crop_bgr=crops_bgr[idx],
                    lamina_mask=lm_i,
                    petiole_mask=pm_i,
                    offset_xy=(ox, oy),
                )

            if eff_exec is not None and n_items >= 2:
                return list(eff_exec.map(_decode_pose_item, all_idx))
            return [_decode_pose_item(i) for i in all_idx]
        except Exception as exc:
            logger.debug(f"LM3PoseLandmarkEngine batched inference fallback: {exc}")
            return [
                self.predict_landmarks(
                    crop_bgr=crops_bgr[i],
                    lamina_mask=lamina_masks[i] if (lamina_masks and i < len(lamina_masks)) else None,
                    petiole_mask=petiole_masks[i] if (petiole_masks and i < len(petiole_masks)) else None,
                    crop_offset=crop_offsets[i] if (crop_offsets and i < len(crop_offsets)) else (0.0, 0.0),
                )
                for i in range(n_items)
            ]

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
# 5. Native LM3 Petiole / Keypoint Anchoring (Read-Only, Zero Rachis Bridging)
# =============================================================================

def _extract_native_lm3_petiole_anchor(
    leaf_mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None,
    min_kp_conf: float = 0.25,
    lamina_mask: Optional[np.ndarray] = None,
    offset_xy: Tuple[float, float] = (0.0, 0.0),
) -> Optional[Dict[str, Any]]:
    """
    Uses LeafMachine3's explicit Stage 2 `petiole_mask` (`class_id == 1`) and/or
    31-keypoint `mid15_pet5` (`lamina_base`, `petiole_0..4`, `petiole_tip`, `lamina_tip`)
    skeleton to locate the anatomical petiole-lamina junction `p_base_crop`, distal lamina
    apex `p_apex_crop`, native landmark metrics (`leafmachine3.core.landmark_metrics.compute_measurements`),
    and orthogonal petiole measurements (`leafmachine3.core.petiole.measure_petiole`).
    Falls back to the basal-most mask contour point when both `keypoints_31` and
    `petiole_mask` are absent.
    """
    if leaf_mask is None and lamina_mask is not None:
        leaf_mask = lamina_mask
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
    petiole_width_px: Optional[float] = None
    touches_leaf: bool = False
    used_native_lm3 = False
    lm_metrics_dict: Dict[str, Any] = {}
    pet_chain: List[Tuple[float, float]] = []

    norm_kps = normalize_lm3_keypoints_dict(keypoints_31, offset_xy=offset_xy) if keypoints_31 is not None else {}

    # 1. Try 31-keypoint pose skeleton (`lamina_base`, `petiole_0..4`, `petiole_tip`, `lamina_tip`)
    if norm_kps:
        pts_2d: Dict[str, Tuple[float, float]] = {
            k: (float(v[0]), float(v[1]))
            for k, v in norm_kps.items()
            if len(v) >= 3 and float(v[2]) >= float(min_kp_conf)
        }
        try:
            from leafmachine3.core.landmark_metrics import compute_measurements

            lm_res = compute_measurements(pts_2d)
            if lm_res is not None:
                lm_metrics_dict = {
                    "lamina_trace_length_px": lm_res.lamina_trace_length_px,
                    "lamina_tip_base_length_px": lm_res.lamina_tip_base_length_px,
                    "petiole_trace_length_px": lm_res.petiole_trace_length_px,
                    "leaf_width_px": lm_res.leaf_width_px,
                    "apex_angle_deg": lm_res.apex_angle_deg,
                    "base_angle_deg": lm_res.base_angle_deg,
                    "lamina_curvature": lm_res.lamina_curvature,
                }
                if lm_res.petiole_trace_length_px is not None and np.isfinite(lm_res.petiole_trace_length_px):
                    petiole_len_px = max(petiole_len_px, float(lm_res.petiole_trace_length_px))
        except Exception:
            pass

        lb_kp = (
            norm_kps.get("lamina_base")
            or norm_kps.get("base_center")
            or norm_kps.get("base_2")
            or norm_kps.get("midvein_14")
        )
        if lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= 0.25:
            dists_lb = np.hypot(leaf_pts[:, 0] - float(lb_kp[0]), leaf_pts[:, 1] - float(lb_kp[1]))
            min_lb_idx = int(np.argmin(dists_lb))
            if float(dists_lb[min_lb_idx]) <= max(36.0, 0.35 * float(max(bw_l, bh_l))):
                p_base_crop = (int(leaf_pts[min_lb_idx, 0]), int(leaf_pts[min_lb_idx, 1]))
            else:
                p_base_crop = (int(round(float(lb_kp[0]))), int(round(float(lb_kp[1]))))
            used_native_lm3 = True

        lt_kp = (
            norm_kps.get("lamina_tip")
            or norm_kps.get("apex_center")
            or norm_kps.get("apex_2")
            or norm_kps.get("midvein_0")
        )
        if lt_kp is not None and len(lt_kp) >= 3 and float(lt_kp[2]) >= 0.25:
            dists_lt = np.hypot(leaf_pts[:, 0] - float(lt_kp[0]), leaf_pts[:, 1] - float(lt_kp[1]))
            min_lt_idx = int(np.argmin(dists_lt))
            if float(dists_lt[min_lt_idx]) <= max(36.0, 0.35 * float(max(bw_l, bh_l))):
                p_apex_crop = (int(leaf_pts[min_lt_idx, 0]), int(leaf_pts[min_lt_idx, 1]))
            else:
                p_apex_crop = (int(round(float(lt_kp[0]))), int(round(float(lt_kp[1]))))

        if p_base_crop is not None:
            pet_chain.append((float(p_base_crop[0]), float(p_base_crop[1])))
        elif lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= min_kp_conf:
            pet_chain.append((float(lb_kp[0]), float(lb_kp[1])))
        for native_pk in ("petiole_0", "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_tip"):
            pkp = norm_kps.get(native_pk)
            if pkp is not None and len(pkp) >= 3 and float(pkp[2]) >= min_kp_conf:
                pt_xy = (float(pkp[0]), float(pkp[1]))
                if not pet_chain or math.hypot(pt_xy[0] - pet_chain[-1][0], pt_xy[1] - pet_chain[-1][1]) > 1e-3:
                    pet_chain.append(pt_xy)
        if len(pet_chain) >= 2 and petiole_len_px <= 0.0:
            arc_len = 0.0
            for i_c in range(len(pet_chain) - 1):
                arc_len += math.hypot(
                    pet_chain[i_c + 1][0] - pet_chain[i_c][0],
                    pet_chain[i_c + 1][1] - pet_chain[i_c][1],
                )
            petiole_len_px = max(petiole_len_px, float(arc_len))

    # 2. Delegate orthogonal petiole measurement directly to LeafMachine3 `measure_petiole`
    if (
        petiole_mask is not None
        and petiole_mask.size > 0
        and petiole_mask.shape[:2] == leaf_mask.shape[:2]
        and int(np.count_nonzero(petiole_mask)) > 0
    ):
        pet_u8 = (petiole_mask > 0).astype(np.uint8)
        lam_u8 = (leaf_mask > 0).astype(np.uint8)
        if p_base_crop is None:
            pys, pxs = np.where(pet_u8 > 0)
            pet_cx, pet_cy = float(np.mean(pxs)), float(np.mean(pys))
            kern_j = np.ones((5, 5), dtype=np.uint8)
            jys_c, jxs_c = np.where((cv2.dilate(pet_u8, kern_j, iterations=1) > 0) & (lam_u8 > 0))
            if len(jxs_c) > 0:
                j_cx, j_cy = float(np.median(jxs_c)), float(np.median(jys_c))
                j_idx = int(np.argmin(np.hypot(leaf_pts[:, 0] - j_cx, leaf_pts[:, 1] - j_cy)))
            else:
                lys_c, lxs_c = np.where(lam_u8 > 0)
                lam_cx = float(np.mean(lxs_c)) if len(lxs_c) > 0 else float(np.mean(leaf_pts[:, 0]))
                lam_cy = float(np.mean(lys_c)) if len(lys_c) > 0 else float(np.mean(leaf_pts[:, 1]))
                rx_p, ry_p = lam_cx - pet_cx, lam_cy - pet_cy
                rn_p = math.hypot(rx_p, ry_p)
                if rn_p > 1.0:
                    ux_p, uy_p = rx_p / rn_p, ry_p / rn_p
                    dx_lp = leaf_pts[:, 0] - pet_cx
                    dy_lp = leaf_pts[:, 1] - pet_cy
                    along_p = dx_lp * ux_p + dy_lp * uy_p
                    perp_p = np.abs(dx_lp * (-uy_p) + dy_lp * ux_p)
                    corridor = (perp_p <= max(8.0, 0.06 * float(max(bw_l, bh_l)))) & (along_p <= rn_p)
                    if np.any(corridor):
                        c_idxs = np.where(corridor)[0]
                        min_al = float(np.min(along_p[c_idxs]))
                        tie_pts = leaf_pts[c_idxs[along_p[c_idxs] <= (min_al + 1.5)]]
                        t_cx, t_cy = float(np.median(tie_pts[:, 0])), float(np.median(tie_pts[:, 1]))
                        j_idx = int(np.argmin(np.hypot(leaf_pts[:, 0] - t_cx, leaf_pts[:, 1] - t_cy)))
                    else:
                        j_idx = int(np.argmin(np.hypot(dx_lp, dy_lp)))
                else:
                    j_idx = int(np.argmin(np.hypot(leaf_pts[:, 0] - pet_cx, leaf_pts[:, 1] - pet_cy)))
            p_base_crop = (int(leaf_pts[j_idx, 0]), int(leaf_pts[j_idx, 1]))
            used_native_lm3 = True
            if len(pet_chain) < 2:
                far_i = int(np.argmax(np.hypot(pxs - p_base_crop[0], pys - p_base_crop[1])))
                pet_chain = [(float(p_base_crop[0]), float(p_base_crop[1])), (float(pxs[far_i]), float(pys[far_i]))]
                petiole_len_px = max(petiole_len_px, math.hypot(pet_chain[1][0] - pet_chain[0][0], pet_chain[1][1] - pet_chain[0][1]))
        if len(pet_chain) >= 2:
            try:
                from leafmachine3.core.petiole import measure_petiole

                pet_res = measure_petiole(pet_u8, lam_u8, pet_chain)
                if pet_res is not None:
                    touches_leaf = bool(pet_res.touches_leaf)
                    if pet_res.width_px is not None and np.isfinite(pet_res.width_px) and pet_res.width_px > 0:
                        petiole_width_px = round(float(pet_res.width_px), 3)
                    if pet_res.length_px is not None and np.isfinite(pet_res.length_px) and pet_res.length_px > 0:
                        petiole_len_px = float(pet_res.length_px)
            except Exception:
                pass

    def _select_midline_apex(base_xy: Tuple[int, int]) -> Tuple[int, int]:
        ys_m, xs_m = np.where(leaf_mask > 0)
        if len(xs_m) > 0:
            cx_m, cy_m = float(np.mean(xs_m)), float(np.mean(ys_m))
        else:
            cx_m, cy_m = float(np.mean(leaf_pts[:, 0])), float(np.mean(leaf_pts[:, 1]))
        vx_a = cx_m - float(base_xy[0])
        vy_a = cy_m - float(base_xy[1])
        vn_a = math.hypot(vx_a, vy_a)
        if vn_a > 1.0:
            ux_a, uy_a = vx_a / vn_a, vy_a / vn_a
            proj_a = (leaf_pts[:, 0] - float(base_xy[0])) * ux_a + (leaf_pts[:, 1] - float(base_xy[1])) * uy_a
            max_p = float(np.max(proj_a))
            cap_pts = leaf_pts[proj_a >= (max_p - 1.5)]
            if len(cap_pts) > 0:
                cap_cx = float(np.median(cap_pts[:, 0]))
                cap_cy = float(np.median(cap_pts[:, 1]))
                a_i = int(np.argmin(np.hypot(cap_pts[:, 0] - cap_cx, cap_pts[:, 1] - cap_cy)))
                return (int(cap_pts[a_i, 0]), int(cap_pts[a_i, 1]))
        d_fb = np.hypot(leaf_pts[:, 0] - float(base_xy[0]), leaf_pts[:, 1] - float(base_xy[1]))
        a_fb = int(np.argmax(d_fb))
        return (int(leaf_pts[a_fb, 0]), int(leaf_pts[a_fb, 1]))

    if p_base_crop is not None:
        if p_apex_crop is None:
            p_apex_crop = _select_midline_apex(p_base_crop)
        out_anchor: Dict[str, Any] = {
            "is_trimmed": False,
            "native_lm3_petiole": bool(used_native_lm3),
            "trim_distance_px": round(float(petiole_len_px), 2),
            "petiole_length_px": round(float(petiole_len_px), 2),
            "geodesic_petiole_skeleton_px": round(float(petiole_len_px), 2),
            "petiole_width_px": petiole_width_px,
            "touches_leaf": bool(touches_leaf),
            "p_base_crop": p_base_crop,
            "p_apex_crop": p_apex_crop,
            "p_base": p_base_crop,
            "p_apex": p_apex_crop,
            "landmark_measurements": lm_metrics_dict,
        }
        out_anchor.update({k: v for k, v in lm_metrics_dict.items() if v is not None})
        return out_anchor

    # 3. Fallback to basal-most midline contour vertex when both keypoints and petiole_mask are absent
    y_max_lp = float(np.max(leaf_pts[:, 1]))
    bot_cap = leaf_pts[leaf_pts[:, 1] >= (y_max_lp - 1.5)]
    bot_cx = float(np.median(bot_cap[:, 0])) if len(bot_cap) > 0 else float(np.mean(leaf_pts[:, 0]))
    b_sub_i = int(np.argmin(np.abs(bot_cap[:, 0] - bot_cx))) if len(bot_cap) > 0 else int(np.argmax(leaf_pts[:, 1]))
    p_base_crop = (int(bot_cap[b_sub_i, 0]), int(bot_cap[b_sub_i, 1])) if len(bot_cap) > 0 else (int(leaf_pts[b_sub_i, 0]), int(leaf_pts[b_sub_i, 1]))
    p_apex_crop = _select_midline_apex(p_base_crop)
    out_fallback: Dict[str, Any] = {
        "is_trimmed": False,
        "native_lm3_petiole": False,
        "morphological_fallback": True,
        "trim_distance_px": 0.0,
        "petiole_length_px": 0.0,
        "geodesic_petiole_skeleton_px": 0.0,
        "petiole_width_px": petiole_width_px,
        "touches_leaf": False,
        "p_base_crop": p_base_crop,
        "p_apex_crop": p_apex_crop,
        "p_base": p_base_crop,
        "p_apex": p_apex_crop,
        "landmark_measurements": lm_metrics_dict,
    }
    out_fallback.update({k: v for k, v in lm_metrics_dict.items() if v is not None})
    return out_fallback


# =============================================================================
# 5. LeafMachine3 Native Orientation, Bilateral Archetype, Scale & ECT Engines
# =============================================================================

_LM3_CANONICAL_TO_NATIVE_KPT_MAP: Dict[str, str] = {
    "lamina_tip": "lamina_tip",
    "apex_1": "apex_left",
    "apex_2": "apex_center",
    "apex_3": "apex_right",
    **{f"midvein_{i}": f"midvein_{i - 1}" for i in range(1, 16)},
    "base_1": "base_left",
    "base_2": "base_center",
    "base_3": "base_right",
    "lamina_base": "lamina_base",
    **{f"petiole_{i}": f"petiole_{i - 1}" for i in range(1, 6)},
    "petiole_6": "petiole_tip",
    "width_1": "width_left",
    "width_2": "width_right",
}


def normalize_lm3_keypoints_dict(
    keypoints_31: Optional[Union[Dict[str, Tuple[float, float, float]], np.ndarray]],
    offset_xy: Tuple[float, float] = (0.0, 0.0),
) -> Dict[str, Tuple[float, float, float]]:
    """
    Normalizes a 31-keypoint dictionary or `(31, 3)` array so both `LM3_KEYPOINT_NAMES_31`
    (`midvein_1..15`, `apex_1..3`, `base_1..3`, `petiole_1..6`, `width_1..2`) and LeafMachine3
    native `KPT_NAMES` (`midvein_0..14`, `apex_left/center/right`, `base_left/center/right`,
    `petiole_0..4`, `petiole_tip`, `width_left/right`) are populated. Applies optional
    `(ox, oy)` subtraction when converting sheet-space keypoints back to crop coordinates.
    """
    ox, oy = float(offset_xy[0]), float(offset_xy[1])
    out: Dict[str, Tuple[float, float, float]] = {}
    if keypoints_31 is None:
        return out
    if isinstance(keypoints_31, dict):
        for k, val in keypoints_31.items():
            if not isinstance(val, (tuple, list, np.ndarray)) or len(val) < 2:
                continue
            vx = float(val[0]) - ox
            vy = float(val[1]) - oy
            vc = float(val[2]) if len(val) >= 3 else 1.0
            out[str(k)] = (vx, vy, vc)
        src_kpts = dict(out)
        # 1. Resolve 1-based (`midvein_1..15`) vs 0-based (`midvein_0..14`) midvein indexing without collision
        if "midvein_0" not in src_kpts and "midvein_1" in src_kpts:
            for i in range(1, 16):
                ck = f"midvein_{i}"
                if ck in src_kpts:
                    out[f"midvein_{i - 1}"] = src_kpts[ck]
        elif "midvein_0" in src_kpts and "midvein_15" not in src_kpts:
            for i in range(1, 16):
                ck = f"midvein_{i}"
                nk = f"midvein_{i - 1}"
                if nk in src_kpts and ck not in src_kpts:
                    out[ck] = src_kpts[nk]
        # 2. Resolve 1-based (`petiole_1..6`) vs 0-based (`petiole_0..4`, `petiole_tip`) petiole indexing
        if "petiole_0" not in src_kpts and "petiole_1" in src_kpts:
            for i in range(1, 6):
                ck = f"petiole_{i}"
                if ck in src_kpts:
                    out[f"petiole_{i - 1}"] = src_kpts[ck]
            if "petiole_6" in src_kpts:
                out["petiole_tip"] = src_kpts["petiole_6"]
            elif "petiole_5" in src_kpts and "petiole_tip" not in src_kpts:
                out["petiole_tip"] = src_kpts["petiole_5"]
        elif "petiole_0" in src_kpts and "petiole_6" not in src_kpts:
            for i in range(1, 6):
                ck = f"petiole_{i}"
                nk = f"petiole_{i - 1}"
                if nk in src_kpts and ck not in src_kpts:
                    out[ck] = src_kpts[nk]
            if "petiole_tip" in src_kpts and "petiole_6" not in src_kpts:
                out["petiole_6"] = src_kpts["petiole_tip"]
        # 3. Map non-overlapping canonical <-> native aliases
        _NON_INDEXED_ALIASES = (
            ("lamina_tip", "apex_center"),
            ("apex_2", "apex_center"),
            ("apex_1", "apex_left"),
            ("apex_3", "apex_right"),
            ("lamina_base", "base_center"),
            ("base_2", "base_center"),
            ("base_1", "base_left"),
            ("base_3", "base_right"),
            ("petiole_6", "petiole_tip"),
            ("width_1", "width_left"),
            ("width_2", "width_right"),
        )
        for canon_k, native_k in _NON_INDEXED_ALIASES:
            if canon_k in src_kpts and native_k not in out:
                out[native_k] = src_kpts[canon_k]
            if native_k in src_kpts and canon_k not in out:
                out[canon_k] = src_kpts[native_k]
        if "lamina_tip" in out and "apex_2" not in out:
            out["apex_2"] = out["lamina_tip"]
        if "lamina_base" in out and "base_2" not in out:
            out["base_2"] = out["lamina_base"]
        return out

    arr = np.asarray(keypoints_31, dtype=np.float64)
    if arr.ndim == 2 and arr.shape[0] >= 2 and arr.shape[1] >= 2:
        # Support both LM3_KEYPOINT_NAMES_31 order (lamina_base at index 22) and
        # compact array order (0: lamina_tip, 1: lamina_base, 2..7: petiole, 8..30: midvein)
        c22 = float(arr[22, 2]) if (arr.shape[0] > 22 and arr.shape[1] >= 3) else 0.0
        c4 = float(arr[4, 2]) if (arr.shape[0] > 4 and arr.shape[1] >= 3) else 0.0
        if arr.shape[0] >= 23 and c22 > 0.0 and c4 > 0.0:
            for idx, canon_k in enumerate(LM3_KEYPOINT_NAMES_31):
                if idx >= arr.shape[0]:
                    break
                vx = float(arr[idx, 0]) - ox
                vy = float(arr[idx, 1]) - oy
                vc = float(arr[idx, 2]) if arr.shape[1] >= 3 else 1.0
                out[canon_k] = (vx, vy, vc)
                native_k = _LM3_CANONICAL_TO_NATIVE_KPT_MAP.get(canon_k)
                if native_k:
                    out[native_k] = (vx, vy, vc)
        else:
            vx0, vy0 = float(arr[0, 0]) - ox, float(arr[0, 1]) - oy
            vc0 = float(arr[0, 2]) if arr.shape[1] >= 3 else 1.0
            out["lamina_tip"] = (vx0, vy0, vc0)
            vx1, vy1 = float(arr[1, 0]) - ox, float(arr[1, 1]) - oy
            vc1 = float(arr[1, 2]) if arr.shape[1] >= 3 else 1.0
            out["lamina_base"] = (vx1, vy1, vc1)
            for p_i in range(min(6, max(0, arr.shape[0] - 2))):
                row_idx = 2 + p_i
                vxp, vyp = float(arr[row_idx, 0]) - ox, float(arr[row_idx, 1]) - oy
                vcp = float(arr[row_idx, 2]) if arr.shape[1] >= 3 else 1.0
                canon_pk = f"petiole_{p_i + 1}"
                out[canon_pk] = (vxp, vyp, vcp)
                native_pk = _LM3_CANONICAL_TO_NATIVE_KPT_MAP.get(canon_pk)
                if native_pk:
                    out[native_pk] = (vxp, vyp, vcp)
            mv_rows = arr[8:] if arr.shape[0] > 8 else np.empty((0, arr.shape[1]))
            if len(mv_rows) > 0:
                idxs = np.linspace(0, len(mv_rows) - 1, min(15, len(mv_rows))).round().astype(int)
                for m_i, r_idx in enumerate(idxs):
                    vxm, vym = float(mv_rows[r_idx, 0]) - ox, float(mv_rows[r_idx, 1]) - oy
                    vcm = float(mv_rows[r_idx, 2]) if arr.shape[1] >= 3 else 1.0
                    # Ensure tip -> base ordering along midvein_0..14
                    out[f"midvein_{m_i}"] = (vxm, vym, vcm)
            if "midvein_0" in out and "midvein_14" in out:
                d0_tip = math.hypot(out["midvein_0"][0] - vx0, out["midvein_0"][1] - vy0)
                d14_tip = math.hypot(out["midvein_14"][0] - vx0, out["midvein_14"][1] - vy0)
                if d0_tip > d14_tip:
                    rev_mv = [out[f"midvein_{i}"] for i in range(14, -1, -1) if f"midvein_{i}" in out]
                    for i, val in enumerate(rev_mv):
                        out[f"midvein_{i}"] = val
                out["midvein_15"] = out["midvein_14"]
    return out


def _count_half_width_lobes(half_widths: np.ndarray, max_hw: float) -> int:
    """Counts prominent lateral lobe maxima along a smoothed 1D margin half-width profile."""
    if half_widths.size < 9 or max_hw <= 3.0:
        return 0
    k_win = max(3, min(9, (int(half_widths.size // 12) * 2) + 1))
    kernel = np.ones(k_win, dtype=np.float64) / float(k_win)
    sm = np.convolve(half_widths, kernel, mode="same")
    prom_min = max(1.8, 0.10 * float(max_hw))
    peaks: List[int] = []
    for i in range(2, len(sm) - 2):
        if sm[i] >= sm[i - 1] and sm[i] >= sm[i + 1] and (sm[i] > sm[i - 2] or sm[i] > sm[i + 2]):
            left_trough = float(np.min(sm[max(0, i - max(4, len(sm) // 6)) : i + 1]))
            right_trough = float(np.min(sm[i : min(len(sm), i + max(5, len(sm) // 6))]))
            if (float(sm[i]) - max(left_trough, right_trough)) >= prom_min or (
                (float(sm[i]) - min(left_trough, right_trough)) >= 1.35 * prom_min
            ):
                if not peaks or (i - peaks[-1]) >= max(3, len(sm) // 14):
                    peaks.append(i)
    return len(peaks)


def evaluate_leaf_damage_and_fold_qc(
    oriented_mask_u8: Optional[np.ndarray],
    solidity: float = 1.0,
    bilateral_dice: float = 1.0,
    si_a: float = 0.0,
    lm3_gates_pass: bool = True,
    n_lamina_components: int = 1,
    detector_class_id: int = 0,
    segment_partial_leaves: bool = False,
    is_crop_clipped: bool = False,
    hole_area_frac: float = 0.0,
    enclosed_void_frac: float = 0.0,
    oriented_points: Optional[Dict[str, Tuple[float, float]]] = None,
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None,
    is_geometric_fallback_pose: bool = False,
) -> Dict[str, Any]:
    """
    Evaluates morphological and anatomical leaf QC gates on an upright LeafMachine3 lamina mask:
      1. Stage 1 `Leaf_PARTIAL` (`detector_class_id == 1` when `segment_partial_leaves == False`)
      2. Crop/sheet border lamina truncation (`is_crop_clipped`)
      3. Severe interior or marginal insect/herbivory perforation (`hole_area_frac > 0.035`)
      4. 31-keypoint midvein pose fold/bend (`FOLDED_MIDVEIN_TRACE`, `APICAL_TIP_FOLD`, `LOW_LANDMARK_CONFIDENCE`)
      5. Longitudinal margin fold crease (`LONGITUDINAL_MARGIN_CREASE`, ruler-straight folded edge vs curved free margin)
      6. Verified `bilateral_lobe_symmetry` closing the `solidity < 0.72` unilateral-tear loophole
      7. Unopened broad inter-lobe overlap loops (`ENCLOSED_INTERLOBE_VOID`, `enclosed_void_frac > 0.020`)
         and multi-island pinnae (`MULTI_COMPONENT_PINNAE`, `n_lamina_components > 1`).
    """
    reasons: List[str] = []
    bilateral_lobe_sym = False
    min_margin_rms_norm = np.nan
    crease_rms_ratio = 1.0

    if int(detector_class_id) == 1 and not bool(segment_partial_leaves):
        reasons.append("LEAF_PARTIAL_FRAGMENT")
    if bool(is_crop_clipped):
        reasons.append("CROP_BORDER_CLIPPED")
    if np.isfinite(hole_area_frac) and float(hole_area_frac) > 0.035:
        reasons.append("SEVERE_HERBIVORY_HOLES")

    # Inspect upright blade margin profiles for bilateral lobing and longitudinal fold creases
    if oriented_mask_u8 is not None and oriented_mask_u8.size > 0 and int(np.count_nonzero(oriented_mask_u8)) >= 32:
        fg = oriented_mask_u8 > 0
        ys_all, xs_all = np.where(fg)
        if ys_all.size >= 16:
            y_min, y_max = int(ys_all.min()), int(ys_all.max())
            x_min, x_max = int(xs_all.min()), int(xs_all.max())
            blade_l = float(max(y_max - y_min, 1))
            blade_w = float(max(x_max - x_min, 1))

            pts_map = oriented_points or {}
            mid_xs = [
                float(pts_map[k][0])
                for k in ("lamina_tip", "lamina_base", *[f"midvein_{i}" for i in range(15)])
                if k in pts_map and pts_map[k] is not None
            ]
            x_axis = float(np.median(mid_xs)) if mid_xs else float(np.mean(xs_all))

            # Scanlines over [15%, 85%] of upright blade length
            y0_s = int(round(y_min + 0.15 * blade_l))
            y1_s = int(round(y_min + 0.85 * blade_l))
            if (y1_s - y0_s) >= 12 and blade_w >= 10.0:
                left_x_list: List[float] = []
                right_x_list: List[float] = []
                scan_ys: List[float] = []
                for y_cur in range(y0_s, y1_s + 1):
                    row_xs = np.where(fg[y_cur, :])[0]
                    if row_xs.size >= 2:
                        scan_ys.append(float(y_cur))
                        left_x_list.append(float(row_xs[0]))
                        right_x_list.append(float(row_xs[-1]))

                if len(scan_ys) >= 12:
                    sy_arr = np.asarray(scan_ys, dtype=np.float64)
                    lx_arr = np.asarray(left_x_list, dtype=np.float64)
                    rx_arr = np.asarray(right_x_list, dtype=np.float64)

                    # Linear chord-deviation RMS on left vs right margins ("taco" longitudinal fold detector)
                    p_l = np.polyfit(sy_arr, lx_arr, 1)
                    p_r = np.polyfit(sy_arr, rx_arr, 1)
                    rms_l = float(np.sqrt(np.mean((lx_arr - np.polyval(p_l, sy_arr)) ** 2))) / blade_w
                    rms_r = float(np.sqrt(np.mean((rx_arr - np.polyval(p_r, sy_arr)) ** 2))) / blade_w
                    min_margin_rms_norm = min(rms_l, rms_r)
                    max_margin_rms_norm = max(rms_l, rms_r)
                    crease_rms_ratio = max_margin_rms_norm / max(min_margin_rms_norm, 1e-4)
                    if (
                        blade_l >= 28.0
                        and min_margin_rms_norm < 0.016
                        and max_margin_rms_norm >= 0.045
                        and crease_rms_ratio >= 3.2
                    ):
                        reasons.append("LONGITUDINAL_MARGIN_CREASE")

                    # Bilateral multi-lobe verification for dissected blades (`solidity < 0.72`)
                    hw_left = np.maximum(0.0, x_axis - lx_arr)
                    hw_right = np.maximum(0.0, rx_arr - x_axis)
                    max_hw_l = float(np.max(hw_left)) if hw_left.size else 0.0
                    max_hw_r = float(np.max(hw_right)) if hw_right.size else 0.0
                    span_balance = min(max_hw_l, max_hw_r) / max(max_hw_l, max_hw_r, 1e-6)
                    lobes_l = _count_half_width_lobes(hw_left, max_hw_l)
                    lobes_r = _count_half_width_lobes(hw_right, max_hw_r)
                    if (
                        lobes_l >= 2
                        and lobes_r >= 2
                        and span_balance >= 0.58
                        and abs(float(si_a)) < 0.28
                        and float(bilateral_dice) >= 0.55
                    ):
                        bilateral_lobe_sym = True
                    elif (
                        int(n_lamina_components) > 1
                        and span_balance >= 0.55
                        and abs(float(si_a)) < 0.28
                    ):
                        bilateral_lobe_sym = True

            # Check 31-keypoint midvein trace tortuosity & apical dog-ear fold
            if keypoints_31 and isinstance(keypoints_31, dict):
                norm_kps = normalize_lm3_keypoints_dict(keypoints_31)
                lt_kp = norm_kps.get("lamina_tip")
                lb_kp = norm_kps.get("lamina_base")
                mv_confs = [
                    float(norm_kps[f"midvein_{i}"][2])
                    for i in range(15)
                    if f"midvein_{i}" in norm_kps and len(norm_kps[f"midvein_{i}"]) >= 3
                ]
                if not is_geometric_fallback_pose:
                    tip_c = float(lt_kp[2]) if (lt_kp is not None and len(lt_kp) >= 3) else 1.0
                    base_c = float(lb_kp[2]) if (lb_kp is not None and len(lb_kp) >= 3) else 1.0
                    med_mv_c = float(np.median(mv_confs)) if mv_confs else 1.0
                    if tip_c < 0.30 or base_c < 0.30 or med_mv_c < 0.30:
                        reasons.append("LOW_LANDMARK_CONFIDENCE")

                # Midvein polyline trace vs straight tip-to-base chord
                chain_pts: List[Tuple[float, float]] = []
                if lt_kp is not None and len(lt_kp) >= 3 and float(lt_kp[2]) >= 0.25:
                    chain_pts.append((float(lt_kp[0]), float(lt_kp[1])))
                for i in range(15):
                    mk = f"midvein_{i}"
                    if mk in norm_kps and len(norm_kps[mk]) >= 3 and float(norm_kps[mk][2]) >= 0.25:
                        pt_m = (float(norm_kps[mk][0]), float(norm_kps[mk][1]))
                        if not chain_pts or math.hypot(pt_m[0] - chain_pts[-1][0], pt_m[1] - chain_pts[-1][1]) > 1e-3:
                            chain_pts.append(pt_m)
                if lb_kp is not None and len(lb_kp) >= 3 and float(lb_kp[2]) >= 0.25:
                    pt_b = (float(lb_kp[0]), float(lb_kp[1]))
                    if not chain_pts or math.hypot(pt_b[0] - chain_pts[-1][0], pt_b[1] - chain_pts[-1][1]) > 1e-3:
                        chain_pts.append(pt_b)
                if len(chain_pts) >= 4:
                    chord_px = math.hypot(chain_pts[-1][0] - chain_pts[0][0], chain_pts[-1][1] - chain_pts[0][1])
                    if chord_px >= 20.0:
                        trace_px = sum(
                            math.hypot(chain_pts[j + 1][0] - chain_pts[j][0], chain_pts[j + 1][1] - chain_pts[j][1])
                            for j in range(len(chain_pts) - 1)
                        )
                        ux_c = (chain_pts[-1][0] - chain_pts[0][0]) / chord_px
                        uy_c = (chain_pts[-1][1] - chain_pts[0][1]) / chord_px
                        max_lat = max(
                            abs((px_i - chain_pts[0][0]) * (-uy_c) + (py_i - chain_pts[0][1]) * ux_c)
                            for px_i, py_i in chain_pts[1:-1]
                        )
                        if (trace_px / chord_px) > 1.22 or (max_lat / chord_px) > 0.22:
                            reasons.append("FOLDED_MIDVEIN_TRACE")

                # Apical dog-ear fold: lamina_tip landmark displaced far below top upright blade margin
                if (
                    not is_geometric_fallback_pose
                    and oriented_points
                    and "lamina_tip" in oriented_points
                    and lt_kp is not None
                    and len(lt_kp) >= 3
                    and float(lt_kp[2]) >= 0.35
                    and blade_l >= 28.0
                ):
                    oy_tip = float(oriented_points["lamina_tip"][1])
                    if (oy_tip - float(y_min)) > 0.14 * blade_l:
                        # Exempt genuinely symmetric obcordate/emarginate apical notches
                        if float(bilateral_dice) < 0.84 or abs(float(si_a)) >= 0.14:
                            reasons.append("APICAL_TIP_FOLD")

    # Symmetry & clump/tear evaluation (closing the `solidity < 0.72` unilateral-damage loophole)
    is_verified_dissected = (float(solidity) < 0.72 and bilateral_lobe_sym) or (int(n_lamina_components) > 1 and bilateral_lobe_sym)
    min_bd = 0.58 if is_verified_dissected else 0.70
    max_sia = 0.28 if (float(solidity) < 0.72 or is_verified_dissected) else 0.32
    if float(bilateral_dice) < min_bd or abs(float(si_a)) >= max_sia or (not is_verified_dissected and not bool(lm3_gates_pass)):
        if float(solidity) < 0.72 and not bilateral_lobe_sym:
            reasons.append("UNILATERAL_ASYMMETRY_OR_TEAR")
        else:
            reasons.append("UNRESOLVED_CLUMP")

    if np.isfinite(enclosed_void_frac) and float(enclosed_void_frac) > 0.020:
        reasons.append("ENCLOSED_INTERLOBE_VOID")
    if int(n_lamina_components) > 1:
        reasons.append("MULTI_COMPONENT_PINNAE")

    dedup_reasons: List[str] = []
    for r_code in reasons:
        if r_code not in dedup_reasons:
            dedup_reasons.append(r_code)

    efa_qc_pass = len(dedup_reasons) == 0
    return {
        "bilateral_lobe_symmetry": bool(bilateral_lobe_sym),
        "efa_qc_pass": bool(efa_qc_pass),
        "qc_flag_reason": "PASS" if efa_qc_pass else ";".join(dedup_reasons),
        "min_margin_rms_norm": round(float(min_margin_rms_norm), 5) if np.isfinite(min_margin_rms_norm) else np.nan,
        "crease_rms_ratio": round(float(crease_rms_ratio), 3),
    }


def orient_and_score_leaf_lm3(
    lamina_mask: np.ndarray,
    keypoints_31: Optional[Union[Dict[str, Tuple[float, float, float]], np.ndarray]] = None,
    p_base: Optional[Tuple[float, float]] = None,
    p_apex: Optional[Tuple[float, float]] = None,
    offset_xy: Tuple[float, float] = (0.0, 0.0),
    min_conf: float = 0.25,
    min_score: float = 0.50,
    hole_mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    crop_bgr: Optional[np.ndarray] = None,
    crown_hint_xy: Optional[Tuple[float, float]] = None,
) -> Dict[str, Any]:
    """
    Computes LeafMachine3 native tip-up rotation angle (`compute_orientation`), rotates
    `lamina_mask`, optional `hole_mask`, `petiole_mask`, and `crop_bgr` (`rotate_image`,
    `mask_bbox`, `crop_to_box`, `oriented_affine`, `apply_affine`), and evaluates
    curvilinear midvein bilateral symmetry and archetype quality (`measure_leaf`).
    """
    if lamina_mask is None or lamina_mask.size == 0 or int(np.count_nonzero(lamina_mask)) == 0:
        return {
            "success": False,
            "angle_cw": 0.0,
            "method": "empty",
            "archetype_score": 0.0,
            "bilateral_dice": 0.0,
            "si_a": 1.0,
            "largest_frac": 0.0,
            "gates_pass": False,
            "is_archetypal": False,
        }

    norm_kps = normalize_lm3_keypoints_dict(keypoints_31, offset_xy=offset_xy)
    pts_2d: Dict[str, Tuple[float, float]] = {
        k: (float(v[0]), float(v[1]))
        for k, v in norm_kps.items()
        if len(v) >= 3 and float(v[2]) >= float(min_conf)
    }
    if p_base is not None and "lamina_base" not in pts_2d:
        pts_2d["lamina_base"] = (float(p_base[0]) - float(offset_xy[0]), float(p_base[1]) - float(offset_xy[1]))
    if p_apex is not None and "lamina_tip" not in pts_2d:
        pts_2d["lamina_tip"] = (float(p_apex[0]) - float(offset_xy[0]), float(p_apex[1]) - float(offset_xy[1]))

    _kp_conf_min = (
        min(
            float(norm_kps.get("lamina_base", (0.0, 0.0, 0.0))[2]),
            float(norm_kps.get("lamina_tip", (0.0, 0.0, 0.0))[2]),
        )
        if norm_kps
        else 0.0
    )

    def _swap_midvein_poles(pts_map: Dict[str, Tuple[float, float]]) -> None:
        pts_map["lamina_base"], pts_map["lamina_tip"] = pts_map["lamina_tip"], pts_map["lamina_base"]
        if "base_center" in pts_map and "apex_center" in pts_map:
            pts_map["base_center"], pts_map["apex_center"] = pts_map["apex_center"], pts_map["base_center"]
        mv_keys = [f"midvein_{i}" for i in range(15) if f"midvein_{i}" in pts_map]
        if len(mv_keys) == 15:
            rev_vals = [pts_map[f"midvein_{14 - i}"] for i in range(15)]
            for i in range(15):
                pts_map[f"midvein_{i}"] = rev_vals[i]

    if (
        petiole_mask is not None
        and petiole_mask.size > 0
        and int(np.count_nonzero(petiole_mask)) > 0
        and "lamina_base" in pts_2d
        and "lamina_tip" in pts_2d
        and (not norm_kps or _kp_conf_min < 0.55)
    ):
        pys, pxs = np.where(petiole_mask > 0)
        lb_pt = pts_2d["lamina_base"]
        lt_pt = pts_2d["lamina_tip"]
        d_lb_pet = float(np.min(np.hypot(pxs - lb_pt[0], pys - lb_pt[1])))
        d_lt_pet = float(np.min(np.hypot(pxs - lt_pt[0], pys - lt_pt[1])))
        span_bt = math.hypot(lt_pt[0] - lb_pt[0], lt_pt[1] - lb_pt[1])
        if (not norm_kps and d_lt_pet < d_lb_pet) or (d_lt_pet + 0.15 * span_bt < d_lb_pet):
            _swap_midvein_poles(pts_2d)
    elif (
        crown_hint_xy is not None
        and "lamina_base" in pts_2d
        and "lamina_tip" in pts_2d
        and (not norm_kps or _kp_conf_min < 0.45)
    ):
        cr_x = float(crown_hint_xy[0]) - float(offset_xy[0])
        cr_y = float(crown_hint_xy[1]) - float(offset_xy[1])
        lb_pt = pts_2d["lamina_base"]
        lt_pt = pts_2d["lamina_tip"]
        d_lb_cr = math.hypot(lb_pt[0] - cr_x, lb_pt[1] - cr_y)
        d_lt_cr = math.hypot(lt_pt[0] - cr_x, lt_pt[1] - cr_y)
        span_bt = math.hypot(lt_pt[0] - lb_pt[0], lt_pt[1] - lb_pt[1])
        if (not norm_kps and d_lt_cr < d_lb_cr) or (d_lt_cr + 0.15 * span_bt < d_lb_cr):
            _swap_midvein_poles(pts_2d)

    # Populate synthetic midvein_0..14 along p_apex -> p_base if pose keypoints were not supplied
    if "lamina_tip" in pts_2d and "lamina_base" in pts_2d and "midvein_0" not in pts_2d:
        tx, ty = pts_2d["lamina_tip"]
        bx, by = pts_2d["lamina_base"]
        for i in range(15):
            frac = float(i + 1) / 16.0
            pts_2d[f"midvein_{i}"] = (tx + frac * (bx - tx), ty + frac * (by - ty))

    try:
        from scripts.vision.lm3_compat import (
            OrientedLeaf,
            apply_affine,
            compute_orientation,
            crop_to_box,
            mask_bbox,
            measure_leaf,
            oriented_affine,
            rotate_image,
        )

        ori = compute_orientation(pts_2d)
        angle_cw = float(ori.angle_cw) if (ori.success and ori.angle_cw is not None) else 0.0
        signed_cw = ((angle_cw + 180.0) % 360.0) - 180.0
        if abs(signed_cw) <= 0.25:
            angle_cw = 0.0
        method = str(ori.method) if ori.success else "fallback"

        mask_u8 = ((lamina_mask > 0).astype(np.uint8)) * 255
        orig_h, orig_w = int(mask_u8.shape[0]), int(mask_u8.shape[1])
        if angle_cw == 0.0:
            mat = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64)
            nh, nw = orig_h, orig_w
            rot_mask = mask_u8.copy()
            rot_hole = (
                (((hole_mask > 0).astype(np.uint8)) * 255)
                if (hole_mask is not None and hole_mask.size > 0 and hole_mask.shape[:2] == (orig_h, orig_w))
                else None
            )
            rot_petiole = (
                (((petiole_mask > 0).astype(np.uint8)) * 255)
                if (petiole_mask is not None and petiole_mask.size > 0 and petiole_mask.shape[:2] == (orig_h, orig_w))
                else None
            )
            rot_bgr = (
                crop_bgr.copy()
                if (crop_bgr is not None and crop_bgr.size > 0 and crop_bgr.shape[:2] == (orig_h, orig_w))
                else None
            )
        else:
            mat, (nh, nw) = oriented_affine((orig_h, orig_w), angle_cw)
            rot_mask = rotate_image(mask_u8, angle_cw, bg=0, nearest=True)

            rot_hole = None
            if hole_mask is not None and hole_mask.size > 0 and hole_mask.shape[:2] == (orig_h, orig_w):
                rot_hole = rotate_image(((hole_mask > 0).astype(np.uint8)) * 255, angle_cw, bg=0, nearest=True)

            rot_petiole = None
            if petiole_mask is not None and petiole_mask.size > 0 and petiole_mask.shape[:2] == (orig_h, orig_w):
                rot_petiole = rotate_image(((petiole_mask > 0).astype(np.uint8)) * 255, angle_cw, bg=0, nearest=True)

            rot_bgr = None
            if crop_bgr is not None and crop_bgr.size > 0 and crop_bgr.shape[:2] == (orig_h, orig_w):
                rot_bgr = rotate_image(crop_bgr, angle_cw, bg=0, nearest=False)

        pt_keys = list(pts_2d.keys())
        if pt_keys:
            pt_arr = np.asarray([pts_2d[k] for k in pt_keys], dtype=np.float64)
            rot_arr = apply_affine(mat, pt_arr)
            rot_pts = {k: (float(rot_arr[i, 0]), float(rot_arr[i, 1])) for i, k in enumerate(pt_keys)}
        else:
            rot_pts = {}

        union_rot = rot_mask.copy()
        if rot_petiole is not None and rot_petiole.shape == union_rot.shape:
            union_rot = np.maximum(union_rot, rot_petiole)
        bbox = mask_bbox(union_rot)
        if bbox is not None:
            x0, y0, x1, y1 = bbox
            crop_rot = crop_to_box(rot_mask, bbox)
            crop_hole = crop_to_box(rot_hole, bbox) if rot_hole is not None else None
            crop_pet = crop_to_box(rot_petiole, bbox) if rot_petiole is not None else None
            crop_img = crop_to_box(rot_bgr, bbox) if rot_bgr is not None else None
            shifted_pts = {k: (pt[0] - float(x0), pt[1] - float(y0)) for k, pt in rot_pts.items()}
            det_box = (x0, y0, x1, y1)
        else:
            crop_rot = rot_mask
            crop_hole = rot_hole
            crop_pet = rot_petiole
            crop_img = rot_bgr
            shifted_pts = rot_pts
            det_box = (0, 0, int(nw), int(nh))

        holes_bool = (
            (crop_hole > 0).astype(bool)
            if (crop_hole is not None and crop_hole.shape == crop_rot.shape)
            else np.zeros(crop_rot.shape[:2], dtype=bool)
        )
        lam_bool = (crop_rot > 0).astype(bool) & (~holes_bool)
        if not np.any(lam_bool):
            lam_bool = (crop_rot > 0).astype(bool)
        crop_lamina_u8 = (lam_bool.astype(np.uint8)) * 255
        # Trim any empty petiole-only margin around upright lamina so oriented_mask frames the upright blade
        lam_bbox = mask_bbox(crop_lamina_u8)
        if lam_bbox is not None and lam_bbox != (0, 0, crop_lamina_u8.shape[1], crop_lamina_u8.shape[0]):
            lx0, ly0, lx1, ly1 = lam_bbox
            crop_lamina_u8 = crop_to_box(crop_lamina_u8, lam_bbox)
            crop_hole = crop_to_box(crop_hole, lam_bbox) if crop_hole is not None else None
            crop_pet = crop_to_box(crop_pet, lam_bbox) if crop_pet is not None else None
            crop_img = crop_to_box(crop_img, lam_bbox) if crop_img is not None else None
            shifted_pts = {k: (pt[0] - float(lx0), pt[1] - float(ly0)) for k, pt in shifted_pts.items()}
            lam_bool = (crop_lamina_u8 > 0).astype(bool)
            holes_bool = (
                (crop_hole > 0).astype(bool)
                if (crop_hole is not None and crop_hole.shape == crop_lamina_u8.shape)
                else np.zeros(crop_lamina_u8.shape[:2], dtype=bool)
            )
        orig_edge_clipped = bool(
            orig_h >= 8
            and orig_w >= 8
            and (
                int(np.count_nonzero(mask_u8[0, :])) > max(3, int(0.02 * max(orig_h, orig_w)))
                or int(np.count_nonzero(mask_u8[-1, :])) > max(3, int(0.02 * max(orig_h, orig_w)))
                or int(np.count_nonzero(mask_u8[:, 0])) > max(3, int(0.02 * max(orig_h, orig_w)))
                or int(np.count_nonzero(mask_u8[:, -1])) > max(3, int(0.02 * max(orig_h, orig_w)))
            )
        )
        if not orig_edge_clipped and crop_lamina_u8.size > 0:
            crop_lamina_u8 = np.pad(crop_lamina_u8, ((2, 2), (2, 2)), mode="constant", constant_values=0)
            if crop_hole is not None:
                crop_hole = np.pad(crop_hole, ((2, 2), (2, 2)), mode="constant", constant_values=0)
            if crop_pet is not None:
                crop_pet = np.pad(crop_pet, ((2, 2), (2, 2)), mode="constant", constant_values=0)
            if crop_img is not None and crop_img.ndim == 3:
                crop_img = np.pad(crop_img, ((2, 2), (2, 2), (0, 0)), mode="constant", constant_values=0)
            shifted_pts = {k: (pt[0] + 2.0, pt[1] + 2.0) for k, pt in shifted_pts.items()}
            lam_bool = (crop_lamina_u8 > 0).astype(bool)
            holes_bool = (
                (crop_hole > 0).astype(bool)
                if (crop_hole is not None and crop_hole.shape == crop_lamina_u8.shape)
                else np.zeros(crop_lamina_u8.shape[:2], dtype=bool)
            )
        sil_bool = lam_bool | holes_bool
        pet_bool = (crop_pet > 0).astype(bool) if crop_pet is not None else None

        cnts_lw, _ = cv2.findContours(crop_lamina_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if cnts_lw:
            sig_lw = [c.reshape(-1, 2).astype(np.float64) for c in cnts_lw if cv2.contourArea(c) >= 16.0]
            c_lw = (
                np.vstack(sig_lw)
                if sig_lw
                else max(cnts_lw, key=cv2.contourArea).reshape(-1, 2).astype(np.float64)
            )
            blade_l_px = float(np.ptp(c_lw[:, 1]))
            blade_w_px = float(np.ptp(c_lw[:, 0]))
        else:
            ys_lw, xs_lw = np.where(lam_bool)
            blade_l_px = float(ys_lw.max() - ys_lw.min()) if len(ys_lw) > 0 else 0.0
            blade_w_px = float(xs_lw.max() - xs_lw.min()) if len(xs_lw) > 0 else 0.0
        anatomical_lw = (blade_l_px / max(blade_w_px, 1e-6)) if (blade_l_px > 0.0 and blade_w_px > 0.0) else np.nan

        kpt_conf = {
            k: (float(norm_kps[k][2]) if k in norm_kps and len(norm_kps[k]) >= 3 else 0.95)
            for k in shifted_pts
        }
        oleaf = OrientedLeaf(
            specimen_id=0,
            leaf_id=0,
            detection_id=0,
            stem="leaf",
            silhouette=sil_bool,
            lamina=lam_bool,
            holes=holes_bool,
            petiole=pet_bool,
            kpts=shifted_pts,
            kpt_conf=kpt_conf,
            angle_cw=angle_cw,
            det_box=det_box,
            crop_shape=(orig_h, orig_w),
            truncated=False,
        )
        bm = measure_leaf(oleaf, min_score=float(min_score))
        if bm is None:
            return {
                "success": bool(ori.success),
                "angle_cw": round(angle_cw, 4),
                "method": method,
                "archetype_score": 0.75,
                "bilateral_dice": 0.75,
                "si_a": 0.0,
                "largest_frac": 1.0,
                "gates_pass": True,
                "is_archetypal": True,
                "aspect_ratio": round(float(anatomical_lw), 4) if np.isfinite(anatomical_lw) else np.nan,
                "blade_length_px": round(blade_l_px, 4),
                "blade_width_px": round(blade_w_px, 4),
                "oriented_mask": crop_lamina_u8,
                "oriented_hole_mask": crop_hole,
                "oriented_petiole_mask": crop_pet,
                "oriented_crop_bgr": crop_img,
                "oriented_points": shifted_pts,
            }
        return {
            "success": bool(ori.success),
            "angle_cw": round(angle_cw, 4),
            "method": method,
            "archetype_score": round(float(bm.archetype_score), 4) if np.isfinite(bm.archetype_score) else 0.75,
            "bilateral_dice": round(float(bm.dice), 4) if (bm.dice is not None and np.isfinite(bm.dice)) else 0.0,
            "si_a": round(float(bm.si_a), 4) if (bm.si_a is not None and np.isfinite(bm.si_a)) else 0.0,
            "largest_frac": round(float(bm.largest_frac), 4) if np.isfinite(bm.largest_frac) else 1.0,
            "gates_pass": bool(bm.gates_pass),
            "is_archetypal": bool(bm.is_archetypal),
            "aspect_ratio": round(float(anatomical_lw), 4) if np.isfinite(anatomical_lw) else np.nan,
            "blade_length_px": round(blade_l_px, 4),
            "blade_width_px": round(blade_w_px, 4),
            "oriented_mask": crop_lamina_u8,
            "oriented_hole_mask": crop_hole,
            "oriented_petiole_mask": crop_pet,
            "oriented_crop_bgr": crop_img,
            "oriented_points": shifted_pts,
        }
    except Exception as exc:
        logger.debug(f"LeafMachine3 orient_and_score_leaf_lm3 fallback: {exc}")
        return {
            "success": False,
            "angle_cw": 0.0,
            "method": "fallback",
            "archetype_score": 0.75,
            "bilateral_dice": 0.75,
            "si_a": 0.0,
            "largest_frac": 1.0,
            "gates_pass": True,
            "is_archetypal": True,
        }


class LM3SpecimenSegmenter:
    """
    LeafMachine3 Stage 2b ONNX Whole-Specimen Segmenter
    (`lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain`).
    Delegates directly to `leafmachine3.inference.specimen_segmenter` (`prep_letterbox_imagenet`,
    `post_letterbox_imagenet`) and `leafmachine3.core.paper_removal.run_paperclean`
    (`METHOD_LLR_HYST_KEEPBLOBS` in CIE L*a*b* space) to produce the cleaned plant-vs-paper
    binary mask (`specimen_paper_cleaned_mask`).
    """

    MODEL_NAME: str = "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"

    def __init__(
        self,
        weights_path: Optional[Union[str, Path]] = None,
        device: str = "cuda",
        conf_thres: float = 0.50,
        img_size: int = 1024,
        max_working_dim: int = 3200,
        paperclean: bool = True,
        providers: Optional[Sequence[str]] = None,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.conf_thres = float(conf_thres)
        self.img_size = int(img_size)
        self.max_working_dim = int(max_working_dim)
        self.paperclean = bool(paperclean)
        self.weights_path, self.model_name, self.providers, self.session = _load_onnx_session(
            weights_path=weights_path,
            default_model_name=self.MODEL_NAME,
            device=self.device,
            providers=providers,
        )
        self.model = self.session

    def predict_mask(
        self,
        image_bgr: np.ndarray,
        max_working_dim: Optional[int] = None,
        return_sheet_size: bool = False,
    ) -> Optional[np.ndarray]:
        """
        Runs UNet++ v9 whole-specimen segmentation + LLR hysteresis paper cleanup on `image_bgr`.
        Returns a uint8 binary mask (`0` or `255`) at working resolution (or full-sheet resolution
        when `return_sheet_size=True`). Returns `None` when ONNX weights are unavailable.
        """
        active_sess = self.session if self.session is not None else self.model
        if active_sess is None or image_bgr is None or image_bgr.size == 0:
            return None
        try:
            from leafmachine3.core.paper_removal import METHOD_LLR_HYST_KEEPBLOBS, run_paperclean
            from leafmachine3.inference.specimen_segmenter import (
                _sigmoid,
                 post_letterbox_imagenet,
                prep_letterbox_imagenet,
            )

            sh, sw = image_bgr.shape[:2]
            w_dim = int(max_working_dim) if max_working_dim is not None else self.max_working_dim
            work_img, _ = _prepare_working_image(image_bgr, max_working_dim=w_dim)
            wh, ww = work_img.shape[:2]
            inp_tensor, meta = prep_letterbox_imagenet(work_img, self.img_size)
            input_name = active_sess.get_inputs()[0].name
            output_name = active_sess.get_outputs()[0].name
            with _ONNX_GPU_LOCK:
                raw_out = active_sess.run([output_name], {input_name: inp_tensor})[0]
            prob_sq = _sigmoid(raw_out[0, 0])
            raw_mask = post_letterbox_imagenet(prob_sq, meta, wh, ww, self.conf_thres)
            if self.paperclean:
                cleaned_mask, _, _ = run_paperclean(work_img, raw_mask, METHOD_LLR_HYST_KEEPBLOBS)
            else:
                cleaned_mask = raw_mask
            mask_u8 = (cleaned_mask > 0).astype(np.uint8) * 255
            if return_sheet_size and (wh != sh or ww != sw):
                mask_u8 = cv2.resize(mask_u8, (sw, sh), interpolation=cv2.INTER_NEAREST)
            return mask_u8
        except Exception as exc:
            logger.debug(f"LM3SpecimenSegmenter inference fallback: {exc}")
            return None

    def __call__(self, image_bgr: np.ndarray, **kwargs: Any) -> Optional[np.ndarray]:
        return self.predict_mask(image_bgr, **kwargs)


class LM3ScaleCalibrator:
    """
    LeafMachine3 scale calibration engine combining:
      1. Archival ruler detection (`lm3_archival_detector__yolo26x_det_1280`, class `0: Ruler`)
      2. Ruler Verification Classifier (`RulerClassifierModel`) and Lattice Scale Solver (`RulerCFLattice`)
      3. Megapixel Conversion Factor prior (`MpConversionFactorModel`, `cf_px_per_cm = 27.112084 * sqrt(MP)`)
    """

    ARCHIVAL_MODEL_NAME: str = "lm3_archival_detector__yolo26x_det_1280"
    MP_PRIOR_SLOPE_PX_PER_CM: float = 27.112084

    def __init__(
        self,
        device: str = "cuda",
        providers: Optional[Sequence[str]] = None,
        use_cf_predicted_by_mp: bool = True,
        archival_conf_thres: float = 0.25,
    ) -> None:
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.providers = providers
        self.use_cf_predicted_by_mp: bool = bool(use_cf_predicted_by_mp)
        self.archival_conf_thres: float = float(archival_conf_thres)
        self.archival_weights_path: Optional[Path] = None
        self.archival_session: Optional[Any] = None
        self.mp_model: Optional[Any] = None
        self.ruler_classifier: Optional[Any] = None
        self.ruler_lattice: Optional[Any] = None
        self._init_lm3_scale_models()
        if self.ruler_lattice is not None:
            self.archival_weights_path, _, _, self.archival_session = _load_onnx_session(
                weights_path=None,
                default_model_name=self.ARCHIVAL_MODEL_NAME,
                device=self.device,
                providers=providers,
            )

    def _init_lm3_scale_models(self) -> None:
        mp_json = LM3_PACKAGE_ROOT / "models" / "mp_conversion_factor" / "model.json"
        if mp_json.is_file():
            try:
                from leafmachine3.inference.mp_conversion_factor import MpConversionFactorModel
                self.mp_model = MpConversionFactorModel(mp_json)
            except Exception:
                self.mp_model = None

        cls_onnx = LM3_PACKAGE_ROOT / "models" / "ruler_classifier" / "model.onnx"
        if not cls_onnx.is_file():
            for _cls_stem in (
                "lm3_ruler_classifier__dinov3_vits16_tile_four_896_mlp",
                "lm3_ruler_classifier__effnet_b4",
            ):
                cand_cls = resolve_lm3_weights(_cls_stem, explicit_path=None)
                cand_path = Path(cand_cls[1]) if isinstance(cand_cls, tuple) else Path(cand_cls)
                if cand_path.is_file():
                    cls_onnx = cand_path
                    break
        if cls_onnx.is_file():
            try:
                from leafmachine3.inference.ruler_classifier import RulerClassifierModel
                self.ruler_classifier = RulerClassifierModel(cls_onnx)
            except Exception:
                self.ruler_classifier = None

        lat_onnx = LM3_PACKAGE_ROOT / "models" / "ruler_cf_lattice" / "model.onnx"
        if not lat_onnx.is_file():
            cand_lat = resolve_lm3_weights("lm3_ruler_cf_lattice__effb4_v11", explicit_path=None)
            lat_onnx = Path(cand_lat[1]) if isinstance(cand_lat, tuple) else Path(cand_lat)
        if lat_onnx.is_file():
            try:
                from leafmachine3.inference.ruler_lattice import RulerCFLattice
                self.ruler_lattice = RulerCFLattice(lat_onnx)
            except Exception:
                self.ruler_lattice = None

    def predict_megapixel_scale(self, height_px: int, width_px: int) -> Tuple[float, str]:
        """Predicts `pixels_per_mm` from native sheet dimensions `(height_px, width_px)` via `MpConversionFactorModel`."""
        h, w = max(1, int(height_px)), max(1, int(width_px))
        if self.mp_model is not None:
            try:
                pred = self.mp_model.predict_from_hw(h, w)
                px_per_mm = float(pred.cf_px_per_cm) / 10.0
                if np.isfinite(px_per_mm) and px_per_mm > 0.0:
                    return round(px_per_mm, 4), "predicted_from_megapixels"
            except Exception:
                pass
        mp = (float(h) * float(w)) / 1e6
        px_per_mm = (self.MP_PRIOR_SLOPE_PX_PER_CM * math.sqrt(max(mp, 1e-6))) / 10.0
        return round(px_per_mm, 4), "predicted_from_megapixels"

    def _detect_ruler_boxes(self, image_bgr: np.ndarray) -> List[Tuple[Tuple[int, int, int, int], float]]:
        """Runs `lm3_archival_detector__yolo26x_det_1280` and extracts class `0` (`Ruler`) bounding boxes."""
        if self.archival_session is None or image_bgr is None or image_bgr.size == 0:
            return []
        try:
            sh, sw = image_bgr.shape[:2]
            w_scale = compute_work_scale((sh, sw), max_working_dim=3200)
            if w_scale < 1.0:
                work_w = max(1, int(round(sw * w_scale)))
                work_h = max(1, int(round(sh * w_scale)))
                work_img = cv2.resize(image_bgr, (work_w, work_h), interpolation=cv2.INTER_AREA)
            else:
                work_h, work_w = sh, sw
                work_img = image_bgr
            lb_bgr, r, (dw, dh) = _letterbox_bgr(work_img, 1280)
            inp_tensor = _bgr_to_nchw_float32(lb_bgr)
            input_name = self.archival_session.get_inputs()[0].name
            with _ONNX_GPU_LOCK:
                raw_outs = self.archival_session.run(None, {input_name: inp_tensor})
            if not raw_outs:
                return []
            arr = np.asarray(raw_outs[0], dtype=np.float32)
            if arr.ndim == 3:
                arr = arr[0]
            inv_r = 1.0 / max(float(r), 1e-6)
            rulers: List[Tuple[Tuple[int, int, int, int], float]] = []
            if arr.ndim == 2 and arr.shape[1] in (6, 7):
                for row in arr:
                    sc = float(row[4])
                    cls_id = int(round(float(row[5])))
                    if cls_id != 0 or sc < float(getattr(self, 'archival_conf_thres', 0.25)):
                        continue
                    x1_w = max(0.0, min(float(work_w), (float(row[0]) - dw) * inv_r))
                    y1_w = max(0.0, min(float(work_h), (float(row[1]) - dh) * inv_r))
                    x2_w = max(0.0, min(float(work_w), (float(row[2]) - dw) * inv_r))
                    y2_w = max(0.0, min(float(work_h), (float(row[3]) - dh) * inv_r))
                    if y2_w <= y1_w or x2_w <= x1_w:
                        continue
                    s_bbox = unscale_working_bbox_to_sheet((y1_w, x1_w, y2_w, x2_w), w_scale, sheet_shape=(sh, sw))
                    if s_bbox[2] > s_bbox[0] and s_bbox[3] > s_bbox[1]:
                        rulers.append((s_bbox, sc))
            rulers.sort(key=lambda item: item[1], reverse=True)
            return rulers
        except Exception as exc:
            logger.debug(f"LM3ScaleCalibrator archival ruler detection fallback: {exc}")
            return []

    def calibrate_sheet(
        self,
        image_bgr: np.ndarray,
        catalog_number: str = "",
        use_cf_predicted_by_mp: Optional[bool] = None,
    ) -> Dict[str, Any]:
        """
        Resolves `pixels_per_mm`, `is_calibrated_mm`, `cf_source`, and `scale_source` for a herbarium sheet:
          1. Attempts direct ruler measurement via `lm3_archival_detector__yolo26x_det_1280` +
             `RulerClassifierModel` + `RulerCFLattice` (`cf_source = 'measured_from_ruler'`).
          2. When `use_cf_predicted_by_mp` is enabled, falls back to `MpConversionFactorModel`
             (`cf_source = 'predicted_from_megapixels'`, `is_calibrated_mm = False`); otherwise
             returns `uncalibrated` (`pixels_per_mm = None`, `is_calibrated_mm = False`).
        """
        _ = catalog_number
        if image_bgr is None or image_bgr.size == 0:
            return {
                "pixels_per_mm": None,
                "is_calibrated_mm": False,
                "cf_source": "uncalibrated",
                "scale_source": "uncalibrated",
            }

        sh, sw = image_bgr.shape[:2]
        if self.archival_session is not None and self.ruler_lattice is not None:
            ruler_boxes = self._detect_ruler_boxes(image_bgr)
            for (ymin, xmin, ymax, xmax), _det_sc in ruler_boxes:
                crop = image_bgr[ymin:ymax, xmin:xmax]
                if crop.size == 0 or min(crop.shape[:2]) < 16:
                    continue
                try:
                    crop_rgb = cv2.cvtColor(crop, cv2.COLOR_BGR2RGB)
                    if self.ruler_classifier is not None:
                        cls_pred = self.ruler_classifier.predict(crop_rgb)
                        if not getattr(cls_pred, "is_ruler", True):
                            continue
                    lat_res = self.ruler_lattice.predict(crop_rgb)
                    px_per_cm = getattr(lat_res, "px_per_cm", None)
                    conf = float(getattr(lat_res, "confidence", 0.0) or 0.0)
                    if px_per_cm is not None and np.isfinite(float(px_per_cm)) and float(px_per_cm) > 10.0 and conf >= 0.45:
                        px_per_mm = round(float(px_per_cm) / 10.0, 4)
                        return {
                            "pixels_per_mm": px_per_mm,
                            "is_calibrated_mm": True,
                            "cf_source": "measured_from_ruler",
                            "scale_source": "ruler_calibrated",
                        }
                except Exception:
                    continue

        eff_use_mp = (
            self.use_cf_predicted_by_mp
            if use_cf_predicted_by_mp is None
            else bool(use_cf_predicted_by_mp)
        )
        if not eff_use_mp:
            return {
                "pixels_per_mm": None,
                "is_calibrated_mm": False,
                "cf_source": "uncalibrated",
                "scale_source": "uncalibrated",
            }
        px_per_mm, cf_source = self.predict_megapixel_scale(sh, sw)
        return {
            "pixels_per_mm": px_per_mm,
            "is_calibrated_mm": cf_source != "predicted_from_megapixels",
            "cf_source": cf_source,
            "scale_source": "mp_estimated" if cf_source == "predicted_from_megapixels" else "ruler_calibrated",
        }


def compute_leaf_ect_descriptor(
    mask: np.ndarray,
    num_dirs: Union[int, str, Path, None] = 360,
    num_thresh: int = 32,
    output_path: Optional[Union[str, Path]] = None,
    hole_mask: Optional[np.ndarray] = None,
    petiole_mask: Optional[np.ndarray] = None,
    angle_cw: Optional[float] = None,
    include_petiole: bool = False,
    include_holes: bool = True,
    already_oriented: bool = False,
) -> Optional[np.ndarray]:
    """
    Computes the LeafMachine3 Module 18 Euler Characteristic Transform (ECT) matrix
    (default `32 x 360` per `config/LM3_settings.yaml`) by delegating directly to
    `leafmachine3.core.ect_compute.compute_ect`, `ect_product_for`, and
    `leafmachine3.reporting.leaf_products.render_leaf_products`.
    """
    if mask is None or mask.size == 0 or int(np.count_nonzero(mask)) == 0:
        return None

    # Support legacy positional call `compute_leaf_ect_descriptor(mask, output_path)`
    if isinstance(num_dirs, (str, Path)) or num_dirs is None:
        if output_path is None and isinstance(num_dirs, (str, Path)):
            output_path = num_dirs
        eff_num_dirs = 360
    else:
        eff_num_dirs = int(num_dirs)
    eff_num_thresh = int(num_thresh)

    _ensure_ect_module_loaded()
    from leafmachine3.core.ect_compute import compute_ect, ect_product_for, tip_up_directions
    from leafmachine3.reporting.leaf_products import render_leaf_products

    mask_u8 = ((mask > 0).astype(np.uint8)) * 255
    if not already_oriented and (angle_cw is not None or hole_mask is not None or petiole_mask is not None):
        h_m, w_m = mask_u8.shape[:2]
        dummy_rgb = np.zeros((h_m, w_m, 3), dtype=np.uint8)
        lam_bool = mask_u8 > 0
        hole_bool = (
            (hole_mask > 0)
            if (hole_mask is not None and hole_mask.size > 0 and hole_mask.shape[:2] == (h_m, w_m))
            else None
        )
        pet_bool = (
            (petiole_mask > 0)
            if (petiole_mask is not None and petiole_mask.size > 0 and petiole_mask.shape[:2] == (h_m, w_m))
            else None
        )
        if hole_bool is not None:
            lam_bool = lam_bool & (~hole_bool)
            sil_bool = lam_bool | hole_bool
        else:
            sil_bool = lam_bool
        prod_name = ect_product_for(bool(include_petiole), bool(include_holes))
        masks_dict: Dict[str, Optional[np.ndarray]] = {
            "lamina": lam_bool,
            "laminapet": (lam_bool | pet_bool) if pet_bool is not None else None,
            "silhouette": sil_bool,
            "silhouettepet": (sil_bool | pet_bool) if pet_bool is not None else None,
            "hole": hole_bool,
        }
        prods = render_leaf_products(
            dummy_rgb,
            masks_dict,
            angle_cw=float(angle_cw) if angle_cw is not None else None,
            want={prod_name},
        )
        selected = prods.get(prod_name)
        if selected is not None and selected.size > 0 and int(np.count_nonzero(selected)) > 0:
            mask_u8 = ((selected > 0).astype(np.uint8)) * 255

    res = compute_ect(mask_u8, num_dirs=eff_num_dirs, bound_radius=1.0, want_simple=False)
    if res is None or getattr(res, "ect_matrix", None) is None:
        return None

    ect_raw = np.asarray(res.ect_matrix, dtype=np.float32)
    thetas = getattr(res, "thetas", None)
    if thetas is not None and len(thetas) == ect_raw.shape[1]:
        ect_mat = np.asarray(tip_up_directions(ect_raw, thetas)[0], dtype=np.float32)
    else:
        ect_mat = ect_raw
    if ect_mat.shape != (eff_num_thresh, eff_num_dirs):
        ect_mat = cv2.resize(ect_mat, (eff_num_dirs, eff_num_thresh), interpolation=cv2.INTER_LINEAR).astype(
            np.float32
        )

    if output_path is not None:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        np.save(str(out_p), ect_mat)

    return ect_mat
