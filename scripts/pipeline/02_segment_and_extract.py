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
      2. LeafMachine3 ONNX Segmentation: Direct inference using canonical LeafMachine3
         YOLO26x plant/leaf/pose engines (`models/lm3/lm3_leaf_segmenter__yolo26x_seg_1024.onnx`)
         and UNet++ v9 specimen segmentation.
      3. LeafMachine3 Scale Calibration: Resolves metric scale from LeafMachine3's
         Module 7 Ruler Lattice Solver and Module 1 Megapixel Conversion Prior.
      4. In-Sheet Organ & Rosette Crown Clustering: Groups organs into plant
         individuals and estimates per-individual root/crown hubs via weighted
         least-squares petiole/midrib ray intersections.
      5. Read-Only LeafMachine3 Ingestion & Native QC Persistence:
           - Persists segmented LeafMachine3 leaves with native bilateral symmetry
             (`bilateral_dice`, `si_a`, `lm3_gates_pass`), unbridged focal pinna
             component counts (`n_lamina_components`), and empirical ruler vs.
             megapixel-prior scale provenance (`scale_source`, `is_calibrated_mm`).
           - Defers closed-contour EFA QC gating to Step 03 while routing sheets
             with zero extractable leaves to data/tables/failed_qc_vouchers.csv.
      6. Standardized Upright Asset Export: Saves upright base-to-tip binary
         lamina masks to data/masks/{catalogNumber}_leaf{id}.png, 120-vertex
         tip-up contours to data/contours/{catalogNumber}_leaf{id}.csv, and
         32 x 360 Euler Characteristic Transform matrices to data/ect/.

Usage:
    python scripts/pipeline/02_segment_and_extract.py \
        --vouchers data/tables/curated_vouchers.csv \
        --model-weights models/lm3/lm3_leaf_segmenter__yolo26x_seg_1024.onnx \
        --output-dir data/ \
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
import re
import shutil
import sys
import tempfile
import threading
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from scripts.vision.lm3_compat import DBSCAN

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import (
    PipelineConfig,
    resolve_cpu_threads,
    resolve_gpu_batch_size,
)
from scripts.core.db import (
    get_completed_qc_vouchers,
    get_db_connection,
    is_canonical_or_explicit_db_context,
    load_step02_records_from_db,
    purge_segments_and_runs,
    record_voucher_segmentation_transaction,
    upsert_vouchers_records,
)
import yaml

from scripts.vision.configure_leafmachine3 import load_native_lm3_config
from scripts.vision.lm3_data_loader import (
    RulerCalibrationMap,
    compute_area_mm2,
    DEFAULT_HERBARIUM_DPI,
    DEFAULT_PX_PER_CM,
    LM3_LEAF_CLASS_IDS,
    LM3_PLANT_CLASSES,
    LM3_REPRODUCTIVE_CLASS_IDS,
    load_ruler_calibrations,
    normalize_pixels_per_mm,
    resolve_px_per_cm,
)
from scripts.vision.leaf_geometry_utils import (
    approximate_petiole_and_scape_metrics,
    cluster_plant_individuals,
    compute_lm3_shape_metrics,
    estimate_rosette_crown_point,
    export_standardized_contour,
    save_leaf_mask_png,
)
from scripts.vision.capitulum_phenology_classifier import (
    classify_single_capitulum,
    classify_voucher_phenology,
)
from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical

from scripts.vision.lm3_compat import (
    _components,
    _largest_component,
    _outer_polygons,
    compute_orientation,
    ect_product_for,
    polygon_mask,
    render_leaf_products,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("SegmentAndExtract")


# =============================================================================
# Data Structures & Two-Stage Crop Reprojection Utilities
# =============================================================================

from scripts.vision.lm3_inference_engines import (
    DetectedInstance,
    LM3LeafSegmenterEngine,
    LM3PlantDetector,
    LM3PoseLandmarkEngine,
    LM3ScaleCalibrator,
    LM3SpecimenSegmenter,
    _VoucherFinalizeBundle,
    _extract_native_lm3_petiole_anchor,
    compute_capitulum_metrics,
    compute_leaf_ect_descriptor,
    crop_leaf_bbox_with_padding,
    is_stock_lm3_weights,
    orient_and_score_leaf_lm3,
    reproject_crop_mask_to_sheet,
    torch_cuda_available,
)


def _build_native_lm3_crop_masks(
    crop_preds: Sequence[Any],
    crop_shape: Tuple[int, int],
    is_high_dpi: bool = True,
    min_leaf_area: int = 0,
    specimen_paper_cleaned_mask: Optional[np.ndarray] = None,
    subtract_holes: bool = True,
) -> Tuple[Optional[np.ndarray], float, Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Constructs `best_mask` (`cls == 0`), `petiole_crop_mask` (`cls == 1`), and `hole_crop_mask` (`cls == 2`)
    delegating component isolation to `leafmachine3.inference.leaf_segmenter._largest_component`,
    `_components`, `_outer_polygons`, and `leafmachine3.core.imaging.polygon_mask`, while isolating
    the focal leaf from overlapping neighbor leaves and retaining strictly unbridged focal pinnae:
      - Rank `cls == 0` (`Leaf`) candidate instance masks by confidence, crop centrality, and area to
        select the primary focal leaf instance (`leaf_primary = _largest_component(primary_inst)`),
        preventing touching neighbor leaf instances from fusing into the focal leaf.
      - Retain all connected pinna components (`>= min_cc_area`) belonging to the primary focal
        `cls == 0` (`Leaf`) instance strictly unbridged (zero synthetic morphological dilation or
        cross-instance bridging), preserving exact LeafMachine3 lamina boundaries.
      - Rasterize outer polygons of `leaf_mask` via `_outer_polygons` + `polygon_mask` as `silhouette`.
      - Union all `cls == 1` (`Petiole`) candidate masks -> keep `_largest_component` (`area >= 64` px).
      - Union all `cls == 2` (`Hole`) candidate masks -> keep `_components(min_area=64)` -> rasterize outer
        polygons via `_outer_polygons` + `polygon_mask` as `hole_crop_mask`.
      - When `subtract_holes=True`, compute `best_mask` as `silhouette & ~hole_crop_mask`, keeping
        `petiole_crop_mask` separate from `best_mask`.
    """
    _ = specimen_paper_cleaned_mask
    crop_h, crop_w = int(crop_shape[0]), int(crop_shape[1])
    if crop_h <= 0 or crop_w <= 0 or not crop_preds:
        return None, 0.0, None, None

    leaf_scores: List[float] = []
    leaf_instances: List[Tuple[float, int, np.ndarray]] = []
    pet_union = np.zeros((crop_h, crop_w), dtype=np.uint8)
    hole_union = np.zeros((crop_h, crop_w), dtype=np.uint8)
    has_pet = False
    has_hole = False

    for _pr_bbox, pr_mask, pr_score, _pr_cls in crop_preds:
        if pr_mask is None or not isinstance(pr_mask, np.ndarray) or pr_mask.size == 0:
            continue
        if pr_mask.shape[:2] != (crop_h, crop_w):
            m_u8 = cv2.resize(
                (pr_mask > 0).astype(np.uint8) * 255,
                (crop_w, crop_h),
                interpolation=cv2.INTER_NEAREST,
            )
        else:
            m_u8 = (pr_mask > 0).astype(np.uint8) * 255
        fg_area = int(np.count_nonzero(m_u8))
        if fg_area <= 0:
            continue
        cls_i = int(_pr_cls)
        if cls_i == 0:
            if min_leaf_area > 0 and fg_area < min_leaf_area:
                continue
            sc_f = float(pr_score)
            leaf_scores.append(sc_f)
            leaf_instances.append((sc_f, fg_area, m_u8 > 0))
        elif cls_i == 1:
            has_pet = True
            np.maximum(pet_union, m_u8, out=pet_union)
        elif cls_i == 2:
            has_hole = True
            np.maximum(hole_union, m_u8, out=hole_union)

    if not leaf_instances:
        return None, 0.0, None, None

    best_score = float(max(leaf_scores))
    is_full_crop = bool(is_high_dpi and (crop_h * crop_w) >= 4096)
    min_cc_area = 64 if is_full_crop else max(1, int(min_leaf_area))

    petiole_crop_mask: Optional[np.ndarray] = None
    if has_pet and np.any(pet_union > 0):
        min_pet_cc = 64 if is_full_crop else 1
        pet_primary = _largest_component(pet_union > 0)
        if int(pet_primary.sum()) >= min_pet_cc:
            pet_sil = np.zeros((crop_h, crop_w), dtype=bool)
            for poly_p in _outer_polygons(pet_primary):
                pet_sil |= polygon_mask(poly_p, (crop_h, crop_w))
            pet_sil |= pet_primary
            petiole_crop_mask = np.where(pet_sil, 255, 0).astype(np.uint8)

    if len(leaf_instances) == 1:
        primary_inst_bool = leaf_instances[0][2]
    else:
        cx_mid = 0.5 * float(crop_w)
        cy_mid = 0.5 * float(crop_h)
        crop_diag = max(1.0, math.hypot(float(crop_w), float(crop_h)))
        best_focal_idx = 0
        best_focal_rank = -1.0
        for idx_li, (sc_li, area_li, mask_li) in enumerate(leaf_instances):
            ys_li, xs_li = np.where(mask_li)
            if len(xs_li) == 0:
                continue
            mx_li = float(np.mean(xs_li))
            my_li = float(np.mean(ys_li))
            dist_norm = math.hypot(mx_li - cx_mid, my_li - cy_mid) / (0.65 * crop_diag)
            cent_w = max(0.35, 1.0 - dist_norm)
            rank_li = float(sc_li) * cent_w * (float(area_li) ** 0.35)
            if rank_li > best_focal_rank:
                best_focal_rank = rank_li
                best_focal_idx = idx_li
        primary_inst_bool = leaf_instances[best_focal_idx][2]

    leaf_primary = _largest_component(primary_inst_bool)
    if int(leaf_primary.sum()) < min_cc_area:
        if is_full_crop:
            return None, 0.0, None, None
        if int(leaf_primary.sum()) <= 0:
            return None, 0.0, None, None

    # Preserve multi-component pinnae belonging to the focal LeafMachine3 leaf instance untouched
    raw_primary_comps = _components(primary_inst_bool, min_area=min_cc_area)
    if isinstance(raw_primary_comps, list) and raw_primary_comps:
        leaf_mask_bool = np.zeros((crop_h, crop_w), dtype=bool)
        for c_m in raw_primary_comps:
            leaf_mask_bool |= np.asarray(c_m, dtype=bool)
    elif isinstance(raw_primary_comps, np.ndarray) and np.any(raw_primary_comps > 0):
        leaf_mask_bool = np.asarray(raw_primary_comps > 0, dtype=bool)
    else:
        leaf_mask_bool = leaf_primary.copy()

    silhouette_bool = np.zeros((crop_h, crop_w), dtype=bool)
    for poly_l in _outer_polygons(leaf_mask_bool):
        silhouette_bool |= polygon_mask(poly_l, (crop_h, crop_w))
    silhouette_bool |= leaf_mask_bool
    silhouette = np.where(silhouette_bool, 255, 0).astype(np.uint8)

    hole_crop_mask = np.zeros((crop_h, crop_w), dtype=np.uint8)
    if has_hole and np.any(hole_union > 0):
        min_hole_cc = 64 if is_full_crop else 1
        raw_hole_comps = _components(hole_union > 0, min_area=min_hole_cc)
        hole_comp_list = (
            [np.asarray(hc, dtype=bool) for hc in raw_hole_comps if np.any(hc)]
            if isinstance(raw_hole_comps, list)
            else ([np.asarray(raw_hole_comps > 0, dtype=bool)] if np.any(raw_hole_comps) else [])
        )
        hole_bool = np.zeros((crop_h, crop_w), dtype=bool)
        for hc_mask in hole_comp_list:
            for poly_h in _outer_polygons(hc_mask):
                hole_bool |= polygon_mask(poly_h, (crop_h, crop_w))
            hole_bool |= hc_mask
        hole_crop_mask = np.where(hole_bool, 255, 0).astype(np.uint8)

    if subtract_holes:
        best_mask = np.where((silhouette > 0) & (hole_crop_mask == 0), 255, 0).astype(np.uint8)
        if not np.any(best_mask > 0):
            best_mask = silhouette
    else:
        best_mask = silhouette

    hole_out: Optional[np.ndarray] = hole_crop_mask if np.any(hole_crop_mask > 0) else None
    return best_mask, best_score, petiole_crop_mask, hole_out


# =============================================================================
# Core Pipeline Execution Runner
# =============================================================================

class SegmentAndExtractPipeline:
    """Unified LeafMachine3 ONNX detection, specimen cleanup, leaf segmentation, pose landmarking, and extraction runner."""

    def __init__(
        self,
        vouchers_csv: Union[str, Path],
        model_weights: Union[str, Path],
        output_dir: Union[str, Path],
        device: str = "cuda",
        score_thresh: float = 0.35,
        force: bool = False,
        config: Optional[PipelineConfig] = None,
        lm3_project_dir: Optional[Union[str, Path]] = None,
        stage1_weights: Optional[Union[str, Path]] = None,
        ruler_confidence_thresh: float = 0.50,
        db_path: Optional[Union[str, Path]] = None,
        threads: Optional[int] = None,
        gpu_batch_size: Optional[int] = None,
        wave_size: Optional[int] = None,
        prefetch_workers: Optional[int] = None,
        prep_workers: Optional[int] = None,
        post_workers: Optional[int] = None,
    ) -> None:
        self.db_path = db_path
        self.vouchers_csv = Path(vouchers_csv)
        self.model_weights = Path(model_weights)
        self.output_dir = Path(output_dir).resolve()
        if self.output_dir.name != "data":
            self.output_dir = self.output_dir / "data"
        self.device = device
        self.config = config or PipelineConfig.from_yaml()

        self.score_thresh = score_thresh
        self.threads: Optional[int] = resolve_cpu_threads(threads, self.config)
        if self.threads is not None and self.device == "cpu":
            try:
                import torch as _torch_cpu
                _torch_cpu.set_num_threads(max(1, int(self.threads)))
            except Exception:
                pass
        self.gpu_batch_size: int = resolve_gpu_batch_size(
            gpu_batch_size,
            self.config,
            default_cuda=10,
            default_cpu=8,
            use_cuda=(self.device == "cuda"),
        )
        _has_explicit_gpu_bs = (
            gpu_batch_size is not None
            or bool(os.environ.get("PACKERA_GPU_BATCH_SIZE", "").strip())
            or self.config.compute.gpu_batch_size is not None
        )
        self.force = force
        self.ruler_confidence_thresh = ruler_confidence_thresh

        # Resolve LM3 project directory
        self.lm3_project_dir: Optional[Path] = Path(lm3_project_dir) if lm3_project_dir is not None else None

        # Output artifact directories
        self.masks_dir = self.output_dir / "masks"
        self.contours_dir = self.output_dir / "contours"
        self.ect_dir = self.output_dir / "ect"
        self.tables_dir = self.output_dir / "tables"
        self.cropped_patches_dir = self.output_dir / "cropped_patches"
        self.capitula_patches_dir = self.cropped_patches_dir / "capitula"

        for d in [
            self.masks_dir,
            self.contours_dir,
            self.ect_dir,
            self.tables_dir,
            self.capitula_patches_dir,
        ]:
            d.mkdir(parents=True, exist_ok=True)

        # If running inside a run directory and curated_vouchers.csv is not yet present in run_dir/data/tables,
        # snapshot the shared data/tables/curated_vouchers.csv into the run directory.
        if not self.vouchers_csv.exists() and self.vouchers_csv.name == "curated_vouchers.csv":
            shared_vouchers = PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv"
            if shared_vouchers.exists() and shared_vouchers.resolve() != self.vouchers_csv.resolve():
                self.vouchers_csv.parent.mkdir(parents=True, exist_ok=True)
                tmp_v = self.vouchers_csv.with_suffix(".csv.tmp")
                shutil.copy2(shared_vouchers, tmp_v)
                tmp_v.replace(self.vouchers_csv)

        # Load canonical LeafMachine3 runtime settings (native LM3 defaults overlaid with config/LM3_settings.yaml)
        try:
            self.lm3_settings: Dict[str, Any] = load_native_lm3_config(project_root=PROJECT_ROOT)
        except Exception:
            self.lm3_settings = {}
        _lm3_cfg_path = PROJECT_ROOT / "config" / "LM3_settings.yaml"
        if _lm3_cfg_path.is_file():
            try:
                _loaded_lm3 = yaml.safe_load(_lm3_cfg_path.read_text(encoding="utf-8"))
                if isinstance(_loaded_lm3, dict):
                    for _k_top, _v_top in _loaded_lm3.items():
                        if isinstance(_v_top, dict) and isinstance(self.lm3_settings.get(_k_top), dict):
                            _merged_sub = dict(self.lm3_settings[_k_top])
                            for _k_sub, _v_sub in _v_top.items():
                                if isinstance(_v_sub, dict) and isinstance(_merged_sub.get(_k_sub), dict):
                                    _merged_sub[_k_sub] = {**_merged_sub[_k_sub], **_v_sub}
                                else:
                                    _merged_sub[_k_sub] = _v_sub
                            self.lm3_settings[_k_top] = _merged_sub
                        else:
                            self.lm3_settings[_k_top] = _v_top
            except Exception:
                pass
        _mods_cfg = self.lm3_settings.get("modules") if isinstance(self.lm3_settings.get("modules"), dict) else {}

        def _mod_section(name: str, legacy_alias: Optional[str] = None) -> Dict[str, Any]:
            sec = _mods_cfg.get(name)
            if isinstance(sec, dict):
                return sec
            if legacy_alias and isinstance(_mods_cfg.get(legacy_alias), dict):
                return _mods_cfg[legacy_alias]
            top_sec = self.lm3_settings.get(name)
            if isinstance(top_sec, dict):
                return top_sec
            if legacy_alias and isinstance(self.lm3_settings.get(legacy_alias), dict):
                return self.lm3_settings[legacy_alias]
            return {}

        _ect_cfg = _mod_section("ect")
        _seg_cfg = getattr(self.config, "segmentation", None)
        _cfg_ect_dirs = int(getattr(_seg_cfg, "ect_num_dirs", 0) or 0) if _seg_cfg is not None else 0
        self.ect_enabled: bool = bool(_ect_cfg.get("enabled", True))
        self.ect_num_dirs: int = max(8, _cfg_ect_dirs if _cfg_ect_dirs > 0 else int(_ect_cfg.get("num_dirs", 360) or 360))
        self.ect_num_thetas: int = max(
            8,
            int(_ect_cfg.get("num_thresh", _ect_cfg.get("num_thetas", 32)) or 32),
        )
        self.ect_include_petiole: bool = bool(_ect_cfg.get("include_petiole", False))
        self.ect_include_holes: bool = bool(_ect_cfg.get("include_holes", True))

        _pd_cfg = _mod_section("plant_detector")
        _ad_cfg = _mod_section("archival_detector")
        _ss_cfg = _mod_section("specimen_segmenter")
        _ls_cfg = _mod_section("leaf_segmenter")
        _kp_cfg = _mod_section("landmark_detector", legacy_alias="pose_detector")

        _s1_conf_cfg = float(_pd_cfg.get("conf_thres", 0.25) or 0.25)
        _s1_iou_cfg = float(_pd_cfg.get("iou_thres", 0.95) or 0.95)
        _s1_img_size = max(320, int(_pd_cfg.get("img_size", 1280) or 1280))
        _s1_work_dim = max(
            640,
            int(getattr(_seg_cfg, "max_working_dim", 0) or _pd_cfg.get("max_working_dim", 3200) or 3200),
        )
        _ad_conf_cfg = float(_ad_cfg.get("conf_thres", 0.25) or 0.25)
        _s2_min_area = max(1, int(_ls_cfg.get("min_area_px", 64) or 64))
        _s2_img_size = max(320, int(_ls_cfg.get("img_size", 1024) or 1024))
        _kp_conf_cfg = float(_kp_cfg.get("conf_thres", 0.25) or 0.25)
        _kp_img_size = max(256, int(_kp_cfg.get("img_size", 640) or 640))
        _kp_pad_frac = float(_kp_cfg.get("white_pad_frac", 0.10) if _kp_cfg.get("white_pad_frac") is not None else 0.10)
        _ss_conf_cfg = float(_ss_cfg.get("conf_thres", 0.50) or 0.50)
        _ss_img_size = max(320, int(_ss_cfg.get("img_size", 1024) or 1024))
        _ss_work_dim = max(
            640,
            int(getattr(_seg_cfg, "specimen_max_working_dim", 0) or _ss_cfg.get("max_working_dim", 3200) or 3200),
        )
        self.paperclean: bool = bool(
            getattr(_seg_cfg, "paperclean", _ss_cfg.get("paperclean", True))
            if _seg_cfg is not None
            else _ss_cfg.get("paperclean", True)
        )
        self.subtract_holes: bool = bool(getattr(_seg_cfg, "subtract_holes", True)) if _seg_cfg is not None else True
        self.do_apply_landmarks: bool = bool(getattr(_seg_cfg, "do_apply_landmarks", True)) if _seg_cfg is not None else True
        self.use_cf_predicted_by_mp: bool = bool(
            getattr(_seg_cfg, "use_cf_predicted_by_mp", False)
        ) if _seg_cfg is not None else False

        # Preload LM3 ruler scale calibrations (preserving empirical vs. megapixel-prior provenance)
        self.ruler_calibrations: Dict[str, float] = RulerCalibrationMap()
        if self.lm3_project_dir and self.lm3_project_dir.exists():
            self.ruler_calibrations.update(
                load_ruler_calibrations(self.lm3_project_dir, min_confidence=self.ruler_confidence_thresh)
            )
        if self.output_dir.exists():
            self.ruler_calibrations.update(
                load_ruler_calibrations(self.output_dir, min_confidence=self.ruler_confidence_thresh)
            )

        try:
            _is_dummy_test_weights = (
                not self.model_weights.stem.lower().startswith("lm3_")
                and (not self.model_weights.is_file() or self.model_weights.stat().st_size < 1024)
            )
        except OSError:
            _is_dummy_test_weights = True
        self._is_dummy_test_weights = bool(_is_dummy_test_weights)
        _aux_default_weights = self.model_weights if _is_dummy_test_weights else None

        # 1. Stage 1 Full-Sheet Plant Organ Detector (YOLO26x 1280px, 9 classes)
        self.stage1_detector = LM3PlantDetector(
            weights_path=stage1_weights if stage1_weights is not None else _aux_default_weights,
            device=self.device,
            conf_thres=min(self.score_thresh, _s1_conf_cfg),
            iou_thres=_s1_iou_cfg,
            img_size=_s1_img_size,
            max_working_dim=_s1_work_dim,
            lm3_project_dir=self.lm3_project_dir,
        )
        self.plant_detector = self.stage1_detector

        # 2. Stage 2 Cropped Leaf Instance Segmenter (YOLO26x-seg 1024px: Leaf, Petiole, Hole)
        self.model_engine = LM3LeafSegmenterEngine(
            weights_path=self.model_weights,
            device=self.device,
            score_thresh=self.score_thresh,
            min_area_px=_s2_min_area,
            img_size=_s2_img_size,
        )
        self.leaf_segmenter = self.model_engine

        # 3. Stage 2 Cropped Leaf Pose Landmark Detector (YOLO26x-pose 640px, 31 keypoints)
        self.pose_engine = LM3PoseLandmarkEngine(
            weights_path=_aux_default_weights,
            device=self.device,
            conf_thres=_kp_conf_cfg,
            img_size=_kp_img_size,
            white_pad_frac=_kp_pad_frac,
        )

        # 4. LeafMachine3 Archival Ruler Detector + Lattice Solver + Megapixel Scale Calibrator
        self.scale_calibrator = LM3ScaleCalibrator(
            device=self.device,
            providers=([] if _is_dummy_test_weights else None),
            use_cf_predicted_by_mp=self.use_cf_predicted_by_mp,
            archival_conf_thres=_ad_conf_cfg,
        )

        # 5. LeafMachine3 Full-Sheet Specimen Segmenter (UNet++ v9 + CIE L*a*b* paper-clean)
        self.specimen_segmenter = LM3SpecimenSegmenter(
            weights_path=_aux_default_weights,
            device=self.device,
            conf_thres=_ss_conf_cfg,
            img_size=_ss_img_size,
            max_working_dim=_ss_work_dim,
            paperclean=self.paperclean,
            providers=([] if _is_dummy_test_weights else None),
        )

        self.is_stock_lm3 = is_stock_lm3_weights(self.model_weights)
        self.voucher_repro_records: Dict[str, Dict[str, Any]] = {}
        self.voucher_pheno_records: Dict[str, Dict[str, Any]] = {}
        self.voucher_repro_by_pid: Dict[str, Dict[int, Dict[str, Any]]] = {}
        self.voucher_pheno_by_pid: Dict[str, Dict[int, Dict[str, Any]]] = {}
        self.voucher_crowns_by_pid: Dict[str, Dict[int, Tuple[float, float]]] = {}
        self.voucher_rosette_centroids_by_pid: Dict[str, Dict[int, Tuple[float, float]]] = {}
        self.voucher_petiole_by_pid: Dict[str, Dict[int, Dict[str, Any]]] = {}
        self.voucher_contour_coords_by_leaf: Dict[str, Dict[int, List[List[float]]]] = {}
        self.stage2_batch_size: int = max(1, int(self.gpu_batch_size))
        env_wave_size = os.environ.get("PACKERA_WAVE_SIZE", "").strip()
        if wave_size is not None:
            self.stage1_wave_size: int = max(1, int(wave_size))
        elif env_wave_size:
            self.stage1_wave_size = max(1, int(env_wave_size))
        elif _has_explicit_gpu_bs:
            self.stage1_wave_size = max(1, min(64, int(self.gpu_batch_size)))
        else:
            self.stage1_wave_size = 16

        env_prefetch_w = os.environ.get("PACKERA_PREFETCH_WORKERS", "").strip()
        env_prep_w = os.environ.get("PACKERA_PREP_WORKERS", "").strip()
        env_post_w = os.environ.get("PACKERA_POST_WORKERS", "").strip()
        self.prefetch_workers_override: Optional[int] = (
            max(1, int(prefetch_workers))
            if prefetch_workers is not None
            else (max(1, int(env_prefetch_w)) if env_prefetch_w else None)
        )
        self.prep_workers_override: Optional[int] = (
            max(1, int(prep_workers))
            if prep_workers is not None
            else (max(1, int(env_prep_w)) if env_prep_w else None)
        )
        self.post_workers_override: Optional[int] = (
            max(1, int(post_workers))
            if post_workers is not None
            else (max(1, int(env_post_w)) if env_post_w else None)
        )

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

    def _apply_finalize_bundle(self, bundle: _VoucherFinalizeBundle) -> None:
        """Merges a completed per-voucher bundle into pipeline state on the coordinator thread."""
        cat = bundle.catalog_number
        if bundle.repro_record is not None:
            self.voucher_repro_records[cat] = bundle.repro_record
        if bundle.pheno_record is not None:
            self.voucher_pheno_records[cat] = bundle.pheno_record
        if bundle.repro_by_pid:
            self.voucher_repro_by_pid[cat] = bundle.repro_by_pid
        if bundle.pheno_by_pid:
            self.voucher_pheno_by_pid[cat] = bundle.pheno_by_pid
        if bundle.crowns_by_pid:
            self.voucher_crowns_by_pid[cat] = bundle.crowns_by_pid
        if bundle.rosette_centroids_by_pid:
            self.voucher_rosette_centroids_by_pid[cat] = bundle.rosette_centroids_by_pid
        if bundle.by_plant_petiole:
            self.voucher_petiole_by_pid[cat] = bundle.by_plant_petiole
        if bundle.contour_coords_by_leaf:
            self.voucher_contour_coords_by_leaf[cat] = bundle.contour_coords_by_leaf

    def shutdown(self) -> None:
        """Cleans up background thread pools."""
        pass

    def _index_completed_vouchers(self) -> Set[str]:
        """Builds an O(1) lookup set of catalogNumbers that already have mask or contour files."""
        completed: Set[str] = set()
        scan_dirs = [self.contours_dir, self.masks_dir]
        for d in scan_dirs:
            if not d.exists():
                continue
            try:
                with os.scandir(d) as it:
                    for entry in it:
                        if not entry.is_file():
                            continue
                        name = entry.name
                        if "_leaf" in name and (name.endswith(".csv") or name.endswith(".png")):
                            try:
                                if entry.stat().st_size > 0:
                                    completed.add(name.split("_leaf")[0])
                            except OSError:
                                pass
            except OSError:
                pass
        return completed

    def is_voucher_completed(self, catalog_number: str) -> bool:
        """
        Checks whether output artifacts already exist for a voucher specimen.
        Verifies:
          1. Standardized contour CSV: data/contours/{catalogNumber}_leaf*.csv (size > 0)
          2. Binary silhouette mask: data/masks/{catalogNumber}_leaf*.png (size > 0)
             in data/masks/.
        """
        contours = list(self.contours_dir.glob(f"{catalog_number}_leaf*.csv"))
        if contours and any(c.is_file() and c.stat().st_size > 0 for c in contours):
            return True

        masks = list(self.masks_dir.glob(f"{catalog_number}_leaf*.png"))
        return bool(masks and any(m.is_file() and m.stat().st_size > 0 for m in masks))

    def _purge_voucher_disk_artifacts(self, catalog_number: str) -> None:
        """
        Purges all prior mask, contour, and capitula crop files for `catalog_number`
        before re-segmenting under `--force`, preventing orphan files (e.g., `_leaf5.png`
        when a voucher now yields 3 leaves) from lingering.
        """
        cat = str(catalog_number).strip()
        if not cat:
            return
        roots = [self.output_dir]
        if (
            self.force
            and os.environ.get("PACKERA_DEMO_MODE") != "1"
            and is_canonical_or_explicit_db_context(self.output_dir)
        ):
            canon_data = PROJECT_ROOT / "data"
            if canon_data.exists() and canon_data.resolve() != self.output_dir.resolve():
                roots.append(canon_data)

        for base_dir in roots:
            dir_patterns = [
                (base_dir / "contours", f"{cat}_leaf*.csv"),
                (base_dir / "ect", f"{cat}_leaf*_ect.npy"),
                (base_dir / "masks", f"{cat}_leaf*.png"),
                (base_dir / "cropped_patches", f"{cat}_p*.*"),
                (base_dir / "cropped_patches" / "capitula", f"{cat}_cap*.jpg"),
            ]
            for d_path, pat in dir_patterns:
                if not d_path.exists():
                    continue
                for f_old in d_path.glob(pat):
                    try:
                        if f_old.is_file() or f_old.is_symlink():
                            f_old.unlink(missing_ok=True)
                    except OSError:
                        pass

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

    def _unpack_stage1_candidates(
        self,
        catalog_number: str,
        img_bgr: np.ndarray,
        stage1_detections: Sequence[Any],
        pixels_per_mm: Optional[float],
    ) -> Tuple[
        List[Tuple[Tuple[int, int, int, int], float, int, np.ndarray, Tuple[int, int, int, int]]],
        List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]],
        List[Tuple[int, int, int, int]],
        List[Tuple[float, float]],
        bool,
        int,
        bool,
        float,
        float,
    ]:
        """
        Thread-safe CPU Stage 1 -> Stage 2 candidate crop extractor.
        Returns:
          (candidate_s1_crops, capitulum_detections, root_detections, upper_cauline_pts,
           rosette_bbox_fallback_used, rejected_leaflet_delta, s1_is_calibrated_mm,
           s1_px_per_cm, s1_eff_min_cm2)
        """
        sheet_h, sheet_w = img_bgr.shape[:2]

        s1_is_calibrated_mm = bool(
            pixels_per_mm is not None and not np.isnan(pixels_per_mm) and float(pixels_per_mm) > 0.0
        )
        s1_px_per_cm = (
            resolve_px_per_cm(
                pixels_per_mm=pixels_per_mm,
                default_dpi=DEFAULT_HERBARIUM_DPI,
                sheet_shape=(sheet_h, sheet_w),
            )
            if s1_is_calibrated_mm
            else 0.0
        )

        candidate_s1_crops: List[Tuple[Tuple[int, int, int, int], float, int, np.ndarray, Tuple[int, int, int, int]]] = []
        capitulum_detections: List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]] = []
        root_detections: List[Tuple[int, int, int, int]] = []
        upper_cauline_pts: List[Tuple[float, float]] = []
        s1_plant_boxes: List[Tuple[int, int, int, int]] = []
        rosette_bbox_fallback_used = False
        rejected_leaflet_delta = 0
        crop_pad = 0.0
        seg_partial_cfg = bool(getattr(self.config.segmentation, "segment_partial_leaves", False))

        # Native LeafMachine3 Two-Stage Handshake:
        #   Stage 1 (`lm3_plant_detector__yolo26x_det_1280` 9-class YOLO26x) separates
        #   `0: Leaf_WHOLE`, `1: Leaf_PARTIAL`, `2..6: reproductive` (`Seed_Fruit_ONE`,
        #   `Seed_Fruit_MANY`, `Flower_ONE`, `Flower_MANY`, `Bud`), `7: Roots`, and `8: Wood`.
        #   Unpack `Leaf_WHOLE` and rescued `Leaf_PARTIAL` boxes directly to Stage 2
        #   `LM3LeafSegmenterEngine` and `LM3PoseLandmarkEngine`.
        c0_items: List[Tuple[Tuple[int, int, int, int], float, int]] = []
        c1_items: List[Tuple[Tuple[int, int, int, int], float, int]] = []
        for s1_item in stage1_detections:
            if len(s1_item) == 4:
                s1_bbox, s1_mask, s1_score, s1_cls = s1_item
            else:
                s1_bbox, s1_score, s1_cls = s1_item
                s1_mask = None

            if s1_cls in LM3_REPRODUCTIVE_CLASS_IDS:
                capitulum_detections.append((s1_bbox, s1_mask, float(s1_score), int(s1_cls)))
                s1_plant_boxes.append(s1_bbox)
            elif int(s1_cls) in (7, 9):
                root_detections.append(s1_bbox)
                s1_plant_boxes.append(s1_bbox)
            elif int(s1_cls) == 0:
                s1_plant_boxes.append(s1_bbox)
                c0_items.append((s1_bbox, float(s1_score), 0))
            elif int(s1_cls) == 1:
                s1_plant_boxes.append(s1_bbox)
                c1_items.append((s1_bbox, float(s1_score), 1))

        selected_s1_leaves: List[Tuple[Tuple[int, int, int, int], float, int]] = list(c0_items)
        if not c0_items:
            selected_s1_leaves.extend(c1_items)
        elif seg_partial_cfg:
            for c1_box, c1_sc, c1_cl in c1_items:
                area_c1 = float(max(1, (c1_box[2] - c1_box[0]) * (c1_box[3] - c1_box[1])))
                overlaps_c0 = False
                for c0_box, _c0_sc, _c0_cl in c0_items:
                    area_c0 = float(max(1, (c0_box[2] - c0_box[0]) * (c0_box[3] - c0_box[1])))
                    iy0, ix0 = max(c1_box[0], c0_box[0]), max(c1_box[1], c0_box[1])
                    iy1, ix1 = min(c1_box[2], c0_box[2]), min(c1_box[3], c0_box[3])
                    if iy1 > iy0 and ix1 > ix0:
                        inter_c = float((iy1 - iy0) * (ix1 - ix0))
                        if (inter_c / max(min(area_c1, area_c0), 1.0)) >= 0.20:
                            overlaps_c0 = True
                            break
                if not overlaps_c0:
                    selected_s1_leaves.append((c1_box, c1_sc, c1_cl))

        # Confidence-ordered nested-box deduplication:
        #   1. Already-kept box `k` suppresses candidate `b` if symmetric `IoU >= 0.70` OR `k` encloses `b`
        #      (`area_k >= 0.85 * area_b` and `inter / max(area_b, 1.0) >= 0.90`).
        #   2. Conversely, if candidate `b` is a whole-leaf box (`b_cls == 0`) that encloses an already-kept
        #      smaller sub-lobe box `k` (`k_cls == 0`, `area_b >= 1.35 * area_k` and `inter / max(area_k, 1.0) >= 0.90`
        #      and `b_score >= 0.50`), replace `k` in `deduped_s1` with the enclosing parent leaf box `b`.
        if len(selected_s1_leaves) >= 2:
            selected_s1_leaves.sort(key=lambda x: (int(x[2]) == 0, float(x[1])), reverse=True)
            deduped_s1: List[Tuple[Tuple[int, int, int, int], float, int]] = []
            for b_box, b_score, b_cls in selected_s1_leaves:
                area_b = float(max(1, (b_box[2] - b_box[0]) * (b_box[3] - b_box[1])))
                suppressed = False
                replaced_indices: List[int] = []
                for idx_k, (k_box, _k_score, k_cls) in enumerate(deduped_s1):
                    area_k = float(max(1, (k_box[2] - k_box[0]) * (k_box[3] - k_box[1])))
                    iy0, ix0 = max(b_box[0], k_box[0]), max(b_box[1], k_box[1])
                    iy1, ix1 = min(b_box[2], k_box[2]), min(b_box[3], k_box[3])
                    if iy1 > iy0 and ix1 > ix0:
                        inter = float((iy1 - iy0) * (ix1 - ix0))
                        iou = inter / max(area_b + area_k - inter, 1.0)
                        if (
                            int(b_cls) == 0
                            and int(k_cls) == 0
                            and area_b >= 1.35 * area_k
                            and (inter / max(area_k, 1.0)) >= 0.90
                            and float(b_score) >= 0.50
                        ):
                            replaced_indices.append(idx_k)
                        elif iou >= 0.70 or (
                            area_k >= 0.85 * area_b and (inter / max(area_b, 1.0)) >= 0.90
                        ):
                            suppressed = True
                            break
                if suppressed:
                    continue
                if replaced_indices:
                    first_idx = replaced_indices[0]
                    deduped_s1[first_idx] = (b_box, b_score, b_cls)
                    for rem_idx in reversed(replaced_indices[1:]):
                        deduped_s1.pop(rem_idx)
                else:
                    deduped_s1.append((b_box, b_score, b_cls))
            selected_s1_leaves = deduped_s1

        for s1_bbox, s1_score, s1_cls in selected_s1_leaves:
            crop_bgr, crop_box = crop_leaf_bbox_with_padding(img_bgr, s1_bbox, pad_ratio=crop_pad)
            if crop_bgr.size > 0:
                candidate_s1_crops.append((s1_bbox, float(s1_score), int(s1_cls), crop_bgr, crop_box))
        s1_eff_min_cm2 = 0.0

        return (
            candidate_s1_crops,
            capitulum_detections,
            root_detections,
            upper_cauline_pts,
            rosette_bbox_fallback_used,
            rejected_leaflet_delta,
            s1_is_calibrated_mm,
            s1_px_per_cm,
            s1_eff_min_cm2,
        )

    def process_voucher(
        self,
        catalog_number: str,
        image_path: Path,
        taxon: Optional[str] = None,
        preloaded_img: Optional[np.ndarray] = None,
        precomputed_s1_letterbox: Optional[np.ndarray] = None,
        precomputed_s1_detections: Optional[List[Any]] = None,
        precomputed_s2_payload: Optional[Tuple[Any, ...]] = None,
        precomputed_single_stage_detections: Optional[Union[List[Any], Tuple[List[Any], bool]]] = None,
        _return_bundle: bool = False,
    ) -> Any:
        """Processes a single herbarium voucher specimen sheet using LeafMachine3 ONNX engines."""
        if taxon is None:
            taxon = self.vouchers_taxa.get(catalog_number, None)

        def _finish(
            ext_recs: List[Dict[str, Any]],
            fail_rec: Optional[Dict[str, Any]],
            repro_rec: Optional[Dict[str, Any]] = None,
            pheno_rec: Optional[Dict[str, Any]] = None,
            crowns_out: Optional[Dict[int, Tuple[float, float]]] = None,
            centroids_out: Optional[Dict[int, Tuple[float, float]]] = None,
            petiole_out: Optional[Dict[int, Dict[str, Any]]] = None,
            contours_out: Optional[Dict[int, List[List[float]]]] = None,
            repro_pid_out: Optional[Dict[int, Dict[str, Any]]] = None,
            pheno_pid_out: Optional[Dict[int, Dict[str, Any]]] = None,
        ) -> Any:
            bundle = _VoucherFinalizeBundle(
                catalog_number=catalog_number,
                image_path=str(image_path),
                extracted_records=ext_recs,
                failure_record=fail_rec,
                repro_record=repro_rec,
                pheno_record=pheno_rec,
                repro_by_pid=repro_pid_out or {},
                pheno_by_pid=pheno_pid_out or {},
                crowns_by_pid=crowns_out or {},
                rosette_centroids_by_pid=centroids_out or {},
                by_plant_petiole=petiole_out or {},
                contour_coords_by_leaf=contours_out or {},
            )
            if _return_bundle:
                return bundle
            self._apply_finalize_bundle(bundle)
            return bundle.extracted_records, bundle.failure_record

        if preloaded_img is not None:
            img_bgr = preloaded_img
        else:
            if not image_path.exists():
                for _ext in (".webp", ".jpg", ".jpeg", ".png"):
                    _cand_p = image_path.with_suffix(_ext)
                    if _cand_p.exists():
                        image_path = _cand_p
                        break
            if not image_path.exists():
                failure_record = {
                    "catalogNumber": catalog_number,
                    "image_path": str(image_path),
                    "failure_reason": "image_not_found",
                    "details": f"File does not exist at {image_path}"
                }
                return _finish([], failure_record, pheno_rec=classify_voucher_phenology(catalog_number, []))

            img_bgr = cv2.imread(str(image_path))
            if img_bgr is None:
                failure_record = {
                    "catalogNumber": catalog_number,
                    "image_path": str(image_path),
                    "failure_reason": "corrupted_image",
                    "details": "Failed to decode image with OpenCV"
                }
                return _finish([], failure_record, pheno_rec=classify_voucher_phenology(catalog_number, []))

        if self.force:
            self._purge_voucher_disk_artifacts(catalog_number)

        sheet_h, sheet_w = img_bgr.shape[:2]
        scale_f = max(sheet_h, sheet_w) / 5600.0 if max(sheet_h, sheet_w) > 1000 else 0.0
        is_high_dpi = max(sheet_h, sheet_w) > 1000 and scale_f > 0.4

        # 1. Ruler Scale Calibration: LeafMachine3 Module 7 Ruler Lattice Solver & Optional Module 1 Megapixel Prior
        lm3_px_per_mm = self.ruler_calibrations.get(catalog_number)
        _cal_sources = getattr(self.ruler_calibrations, "sources", None)
        if lm3_px_per_mm is not None and not np.isnan(lm3_px_per_mm) and lm3_px_per_mm > 0:
            pixels_per_mm: Optional[float] = normalize_pixels_per_mm(lm3_px_per_mm)
            cf_source = (
                str(_cal_sources.get(catalog_number, "measured_from_ruler"))
                if isinstance(_cal_sources, dict)
                else "measured_from_ruler"
            )
        elif (
            not getattr(self, "_is_dummy_test_weights", True)
            and getattr(self, "scale_calibrator", None) is not None
        ):
            cal_res = self.scale_calibrator.calibrate_sheet(
                img_bgr,
                catalog_number=catalog_number,
                use_cf_predicted_by_mp=bool(getattr(self, "use_cf_predicted_by_mp", False)),
            )
            pixels_per_mm = normalize_pixels_per_mm(cal_res.get("pixels_per_mm"))
            cf_source = str(cal_res.get("cf_source", "uncalibrated"))
            if pixels_per_mm is not None and float(pixels_per_mm) > 0.0:
                self.ruler_calibrations[catalog_number] = float(pixels_per_mm)
                if isinstance(_cal_sources, dict):
                    _cal_sources[catalog_number] = cf_source
        else:
            pixels_per_mm = None
            cf_source = "uncalibrated"

        # 2. Stage 1 Full-Sheet Plant Organ Detection (YOLO26x 9-class detector)
        use_s1_cache = bool(self.lm3_project_dir is not None and not self.force)
        direct_s2_mock = (
            "predict" in getattr(self.model_engine, "__dict__", {})
            and "predict" not in getattr(self.stage1_detector, "__dict__", {})
            and "predict_batch" not in getattr(self.model_engine, "__dict__", {})
        )
        if precomputed_s1_detections is not None:
            stage1_detections = precomputed_s1_detections
        elif direct_s2_mock:
            stage1_detections = []
        elif precomputed_s1_letterbox is not None:
            stage1_detections = self.stage1_detector.predict(
                img_bgr,
                catalog_number=catalog_number,
                use_cache=use_s1_cache,
                precomputed_letterbox=precomputed_s1_letterbox,
            )
        else:
            stage1_detections = self.stage1_detector.predict(
                img_bgr, catalog_number=catalog_number, use_cache=use_s1_cache
            )

        capitulum_detections: List[Tuple[Tuple[int, int, int, int], Optional[np.ndarray], float, int]] = []
        root_detections: List[Tuple[int, int, int, int]] = []
        upper_cauline_pts: List[Tuple[float, float]] = []
        rosette_bbox_fallback_used: bool = False
        # Each tuple: (full_bbox, mask, score, cls_id, crop_offset, mask_origin, already_trimmed, trim_meta, m_solidity, fg_area)
        leaf_detections: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int, Tuple[int, int], Tuple[int, int], bool, Optional[Dict[str, Any]], float, int]] = []

        if stage1_detections:
            all_crop_pose_preds: List[Optional[Any]] = []
            if precomputed_s2_payload is not None:
                if len(precomputed_s2_payload) >= 11:
                    (
                        candidate_s1_crops,
                        capitulum_detections,
                        root_detections,
                        upper_cauline_pts,
                        rosette_bbox_fallback_used,
                        _rej_lflt,
                        s1_is_calibrated_mm,
                        s1_px_per_cm,
                        s1_eff_min_cm2,
                        all_crop_preds,
                        all_crop_pose_preds,
                    ) = precomputed_s2_payload[:11]
                else:
                    (
                        candidate_s1_crops,
                        capitulum_detections,
                        root_detections,
                        upper_cauline_pts,
                        rosette_bbox_fallback_used,
                        _rej_lflt,
                        s1_is_calibrated_mm,
                        s1_px_per_cm,
                        s1_eff_min_cm2,
                        all_crop_preds,
                    ) = precomputed_s2_payload[:10]
            else:
                (
                    candidate_s1_crops,
                    capitulum_detections,
                    root_detections,
                    upper_cauline_pts,
                    rosette_bbox_fallback_used,
                    _rej_lflt,
                    s1_is_calibrated_mm,
                    s1_px_per_cm,
                    s1_eff_min_cm2,
                ) = self._unpack_stage1_candidates(
                    catalog_number=catalog_number,
                    img_bgr=img_bgr,
                    stage1_detections=stage1_detections,
                    pixels_per_mm=pixels_per_mm,
                )
                all_crop_preds = (
                    self.model_engine.predict_batch([c[3] for c in candidate_s1_crops])
                    if candidate_s1_crops
                    else []
                )
                if (
                    candidate_s1_crops
                    and all_crop_preds
                    and getattr(self, "pose_engine", None) is not None
                    and getattr(self.pose_engine, "model", None) is not None
                    and "predict_landmarks" not in getattr(self.pose_engine, "__dict__", {})
                    and hasattr(self.pose_engine, "predict_batch")
                ):
                    all_crop_pose_preds = [None] * len(candidate_s1_crops)
                    p_idxs: List[int] = []
                    p_crops: List[np.ndarray] = []
                    p_lamina_masks: List[Optional[np.ndarray]] = []
                    p_petiole_masks: List[Optional[np.ndarray]] = []
                    for idx_c, (c_item, c_preds) in enumerate(zip(candidate_s1_crops, all_crop_preds)):
                        if not c_preds:
                            continue
                        has_leaf_cand = any(
                            int(p_cls) == 0 and np.any(pm_u8 > 0)
                            for _pb, pm_u8, _p_sc, p_cls in c_preds
                        )
                        if has_leaf_cand:
                            _bm_c, _, _pm_c, _ = _build_native_lm3_crop_masks(
                                c_preds,
                                c_item[3].shape[:2],
                                is_high_dpi=is_high_dpi,
                                min_leaf_area=0,
                            )
                            p_idxs.append(idx_c)
                            p_crops.append(c_item[3])
                            p_lamina_masks.append(_bm_c)
                            p_petiole_masks.append(_pm_c)
                    if p_crops:
                        for idx_c, p_res in zip(
                            p_idxs,
                            self.pose_engine.predict_batch(
                                p_crops,
                                lamina_masks=p_lamina_masks,
                                petiole_masks=p_petiole_masks,
                                batch_size=max(8, self.stage2_batch_size),
                            ),
                        ):
                            all_crop_pose_preds[idx_c] = p_res

            def _reproject_crop_mask_to_local_win(
                cmask: Optional[np.ndarray],
                crop_box_yx: Tuple[int, int, int, int],
                win_box_yx: Tuple[int, int, int, int],
            ) -> Optional[np.ndarray]:
                if cmask is None:
                    return None
                cy0_c, cx0_c, cy1_c, cx1_c = crop_box_yx
                wy0_w, wx0_w, wy1_w, wx1_w = win_box_yx
                wh = wy1_w - wy0_w
                ww = wx1_w - wx0_w
                if wh <= 0 or ww <= 0:
                    return None
                out_m = np.zeros((wh, ww), dtype=np.uint8)
                ch_m, cw_m = cmask.shape[:2]
                iy0 = max(cy0_c, wy0_w)
                ix0 = max(cx0_c, wx0_w)
                iy1 = min(cy0_c + ch_m, wy1_w)
                ix1 = min(cx0_c + cw_m, wx1_w)
                if iy1 > iy0 and ix1 > ix0:
                    out_m[iy0 - wy0_w : iy1 - wy0_w, ix0 - wx0_w : ix1 - wx0_w] = cmask[
                        iy0 - cy0_c : iy1 - cy0_c, ix0 - cx0_c : ix1 - cx0_c
                    ]
                return out_m if np.any(out_m > 0) else None

            def _shift_crop_keypoints_to_local_win(
                kpts: Optional[Any],
                crop_box_yx: Tuple[int, int, int, int],
                win_box_yx: Tuple[int, int, int, int],
            ) -> Optional[Dict[str, Tuple[float, float, float]]]:
                if not kpts or not isinstance(kpts, dict):
                    return None
                cy0_c, cx0_c, _, _ = crop_box_yx
                wy0_w, wx0_w, _, _ = win_box_yx
                dx_w = float(cx0_c - wx0_w)
                dy_w = float(cy0_c - wy0_w)
                shifted: Dict[str, Tuple[float, float, float]] = {}
                for k_name, vals in kpts.items():
                    if isinstance(vals, (tuple, list, np.ndarray)) and len(vals) >= 2:
                        conf_v = float(vals[2]) if len(vals) >= 3 else 1.0
                        shifted[str(k_name)] = (float(vals[0]) + dx_w, float(vals[1]) + dy_w, conf_v)
                return shifted or None

            # Stage 2 LeafMachine3 YOLO26x-seg (`0: Leaf`, `1: Petiole`, `2: Hole`) + YOLO26x-pose (`keypoints_31`)
            for crop_idx_s2, ((s1_bbox, s1_score, s1_cls, crop_bgr, crop_box), crop_preds) in enumerate(zip(candidate_s1_crops, all_crop_preds)):
                if not crop_preds:
                    continue

                best_mask, best_score, best_petiole_mask, best_hole_mask = _build_native_lm3_crop_masks(
                    crop_preds,
                    crop_bgr.shape[:2],
                    is_high_dpi=is_high_dpi,
                    min_leaf_area=0,
                    subtract_holes=bool(getattr(self, "subtract_holes", True)),
                )
                if best_mask is None:
                    continue

                keypoints_31_crop: Optional[Dict[str, Tuple[float, float, float]]] = None
                petiole_width_px_crop: Optional[float] = None
                pre_pose_crop = all_crop_pose_preds[crop_idx_s2] if crop_idx_s2 < len(all_crop_pose_preds) else None
                if pre_pose_crop is not None and bool(getattr(self, "do_apply_landmarks", True)):
                    keypoints_31_crop = getattr(pre_pose_crop, "keypoints_31", None) or None
                    petiole_width_px_crop = getattr(pre_pose_crop, "petiole_width_px", None)
                elif (
                    bool(getattr(self, "do_apply_landmarks", True))
                    and getattr(self, "pose_engine", None) is not None
                    and (
                        getattr(self.pose_engine, "model", None) is not None
                        or "predict_landmarks" in getattr(self.pose_engine, "__dict__", {})
                    )
                ):
                    try:
                        pose_res = self.pose_engine.predict_landmarks(
                            crop_bgr, lamina_mask=best_mask, petiole_mask=best_petiole_mask
                        )
                        keypoints_31_crop = pose_res.keypoints_31 or None
                        petiole_width_px_crop = pose_res.petiole_width_px
                    except Exception:
                        keypoints_31_crop = None

                cy0, cx0, cy1, cx1 = crop_box
                m_raw = best_mask
                fg_area = int(np.count_nonzero(m_raw))
                if fg_area <= 0:
                    continue

                bx_c, by_c, bw_c, bh_c = cv2.boundingRect(m_raw)
                if bw_c == 0 or bh_c == 0:
                    continue
                sub_raw = m_raw[by_c : by_c + bh_c, bx_c : bx_c + bw_c]
                cnts_m, _ = cv2.findContours(sub_raw, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                if not cnts_m:
                    continue
                main_m_cnt = max(cnts_m, key=cv2.contourArea)
                hull_m_area = float(cv2.contourArea(cv2.convexHull(main_m_cnt)))
                m_solidity = float(fg_area) / hull_m_area if hull_m_area > 0 else 1.0

                ymin_s = max(0, cy0 + by_c)
                ymax_s = min(sheet_h, cy0 + by_c + bh_c)
                xmin_s = max(0, cx0 + bx_c)
                xmax_s = min(sheet_w, cx0 + bx_c + bw_c)
                if ymax_s <= ymin_s or xmax_s <= xmin_s:
                    continue
                full_bbox = (ymin_s, xmin_s, ymax_s, xmax_s)

                margin = max(24, int(0.35 * max(bw_c, bh_c)))
                win_y0 = max(0, ymin_s - margin)
                win_x0 = max(0, xmin_s - margin)
                win_y1 = min(sheet_h, ymax_s + margin)
                win_x1 = min(sheet_w, xmax_s + margin)
                local_mask = np.zeros((win_y1 - win_y0, win_x1 - win_x0), dtype=np.uint8)
                dst_y0 = ymin_s - win_y0
                dst_x0 = xmin_s - win_x0
                eff_h = min(bh_c, win_y1 - ymin_s)
                eff_w = min(bw_c, win_x1 - xmin_s)
                local_mask[dst_y0 : dst_y0 + eff_h, dst_x0 : dst_x0 + eff_w] = sub_raw[:eff_h, :eff_w]

                native_pet_meta = _extract_native_lm3_petiole_anchor(
                    m_raw, best_petiole_mask, keypoints_31=keypoints_31_crop
                )
                if native_pet_meta is not None:
                    p_bc = native_pet_meta["p_base_crop"]
                    p_ac = native_pet_meta["p_apex_crop"]
                    native_pet_meta["p_base"] = (
                        int(cx0 + p_bc[0] - win_x0),
                        int(cy0 + p_bc[1] - win_y0),
                    )
                    native_pet_meta["p_apex"] = (
                        int(cx0 + p_ac[0] - win_x0),
                        int(cy0 + p_ac[1] - win_y0),
                    )
                else:
                    native_pet_meta = {}
                native_pet_meta["petiole_mask"] = _reproject_crop_mask_to_local_win(
                    best_petiole_mask, (cy0, cx0, cy1, cx1), (win_y0, win_x0, win_y1, win_x1)
                )
                native_pet_meta["hole_mask"] = _reproject_crop_mask_to_local_win(
                    best_hole_mask, (cy0, cx0, cy1, cx1), (win_y0, win_x0, win_y1, win_x1)
                )
                native_pet_meta["keypoints_31"] = _shift_crop_keypoints_to_local_win(
                    keypoints_31_crop, (cy0, cx0, cy1, cx1), (win_y0, win_x0, win_y1, win_x1)
                )
                native_pet_meta["petiole_width_px"] = petiole_width_px_crop

                leaf_detections.append((
                    full_bbox,
                    local_mask,
                    float(best_score),
                    int(s1_cls),
                    (crop_box[1], crop_box[0]),
                    (win_x0, win_y0),
                    False,
                    native_pet_meta,
                    m_solidity,
                    fg_area,
                ))
        elif not is_high_dpi or self.stage1_detector.model is None or "predict" in self.model_engine.__dict__:
            # Fallback when Stage 1 detector is uninitialized (e.g., direct single-stage unit tests)
            if precomputed_single_stage_detections is not None:
                raw_detections = (
                    precomputed_single_stage_detections[0]
                    if isinstance(precomputed_single_stage_detections, tuple)
                    else precomputed_single_stage_detections
                )
            else:
                raw_detections = self.model_engine.predict(img_bgr)
            raw_s2_leaves: List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]] = []
            raw_s2_petioles: List[Tuple[Tuple[int, int, int, int], np.ndarray, float]] = []
            for bbox, mask, score, cls_id in raw_detections:
                if cls_id in LM3_REPRODUCTIVE_CLASS_IDS:
                    capitulum_detections.append((bbox, mask, float(score), int(cls_id)))
                elif int(cls_id) in (7, 9):
                    root_detections.append(bbox)
                elif int(cls_id) == 0:
                    raw_s2_leaves.append((bbox, mask, float(score), int(cls_id)))
                elif int(cls_id) == 1:
                    raw_s2_petioles.append((bbox, mask, float(score)))

            for bbox, mask, score, cls_id in raw_s2_leaves:
                sheet_mask = mask
                full_bbox = bbox
                fg_area = int(np.count_nonzero(sheet_mask))
                if fg_area <= 0:
                    continue

                native_pet_meta: Optional[Dict[str, Any]] = None
                if raw_s2_petioles and sheet_mask.shape[:2] == (sheet_h, sheet_w):
                    ly0, lx0, ly1, lx1 = [int(v) for v in full_bbox]
                    prox = max(32, int(0.35 * max(ly1 - ly0, lx1 - lx0, 1)))
                    nearby_pets = [
                        (pscore, pbbox, pmask)
                        for pbbox, pmask, pscore in raw_s2_petioles
                        if pmask is not None
                        and pmask.shape[:2] == sheet_mask.shape[:2]
                        and not (
                            int(pbbox[2]) < ly0 - prox
                            or int(pbbox[0]) > ly1 + prox
                            or int(pbbox[3]) < lx0 - prox
                            or int(pbbox[1]) > lx1 + prox
                        )
                    ]
                    if nearby_pets:
                        nearby_pets.sort(key=lambda item: item[0], reverse=True)
                        _, pbbox_best, pmask_best = nearby_pets[0]
                        jy0 = max(0, min(ly0, int(pbbox_best[0])) - 8)
                        jx0 = max(0, min(lx0, int(pbbox_best[1])) - 8)
                        jy1 = min(sheet_h, max(ly1, int(pbbox_best[2])) + 8)
                        jx1 = min(sheet_w, max(lx1, int(pbbox_best[3])) + 8)
                        if jy1 > jy0 and jx1 > jx0:
                            native_pet_meta = _extract_native_lm3_petiole_anchor(
                                sheet_mask[jy0:jy1, jx0:jx1],
                                pmask_best[jy0:jy1, jx0:jx1],
                            )
                            if native_pet_meta is not None:
                                p_bc = native_pet_meta["p_base_crop"]
                                p_ac = native_pet_meta["p_apex_crop"]
                                native_pet_meta["p_base"] = (int(jx0 + p_bc[0]), int(jy0 + p_bc[1]))
                                native_pet_meta["p_apex"] = (int(jx0 + p_ac[0]), int(jy0 + p_ac[1]))

                leaf_detections.append((
                    full_bbox,
                    sheet_mask,
                    float(score),
                    int(cls_id),
                    (int(full_bbox[1]), int(full_bbox[0])),
                    (0, 0),
                    False,
                    native_pet_meta,
                    0.0,
                    fg_area,
                ))

        is_calibrated_mm = bool(
            pixels_per_mm is not None
            and not np.isnan(pixels_per_mm)
            and float(pixels_per_mm) > 0.0
            and str(cf_source).strip().lower() not in {"predicted_from_megapixels", "mp_estimated", "uncalibrated"}
        )

        if not leaf_detections and not capitulum_detections:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "no_leaves_detected",
                "details": "Model produced 0 candidate leaf segmentations",
                "pixels_per_mm": round(float(pixels_per_mm), 4) if is_calibrated_mm else np.nan,
                "is_calibrated_mm": bool(is_calibrated_mm),
            }
            empty_repro = {
                "catalogNumber": catalog_number,
                "plant_individual_id": 0,
                "capitula_count": 0,
                "involucre_height_px": np.nan,
                "involucre_width_px": np.nan,
                "involucre_height_mm": np.nan,
                "involucre_width_mm": np.nan,
                "AR_inv": np.nan,
                "plant_height_cm": np.nan,
                "scape_to_basal_leaf_ratio": np.nan,
                "is_bent_stem": False,
            }
            _pheno_cfg = getattr(getattr(self.config, "thresholds", None), "capitulum_phenology", None)
            _min_anth_heads = int(getattr(_pheno_cfg, "min_anthesis_heads_for_flowering", 1) or 1)
            empty_pheno = classify_voucher_phenology(
                catalog_number,
                [],
                min_anthesis_heads_for_flowering=_min_anth_heads,
            )
            return _finish([], failure_record, repro_rec=empty_repro, pheno_rec=empty_pheno)

        # 3. Route Reproductive Detections (seed_fruit_one/many, flower_one/many, bud) to Phenology Classifier
        _pheno_cfg = getattr(getattr(self.config, "thresholds", None), "capitulum_phenology", None)
        _min_pappus_lum = float(getattr(_pheno_cfg, "min_pappus_luminance", 75.0) if _pheno_cfg is not None else 75.0)
        _max_pappus_sat = float(getattr(_pheno_cfg, "max_pappus_saturation", 0.18) if _pheno_cfg is not None else 0.18)
        _min_pappus_area = float(getattr(_pheno_cfg, "min_pappus_area_ratio", 0.18) if _pheno_cfg is not None else 0.18)
        _min_floret_sat = float(getattr(_pheno_cfg, "min_floret_yellow_saturation", 0.32) if _pheno_cfg is not None else 0.32)
        _min_floret_rat = float(getattr(_pheno_cfg, "min_floret_yellow_ratio", 0.12) if _pheno_cfg is not None else 0.12)
        _min_anth_heads = int(getattr(_pheno_cfg, "min_anthesis_heads_for_flowering", 1) or 1) if _pheno_cfg is not None else 1
        valid_capitula: List[Dict[str, Any]] = []
        head_pheno_results: List[Dict[str, Any]] = []
        raw_cap_items: List[Dict[str, Any]] = []
        for cap_idx, (bbox, mask, score, class_id) in enumerate(capitulum_detections, 1):
            ymin, xmin, ymax, xmax = [int(v) for v in bbox]
            y0, y1 = max(0, ymin), min(sheet_h, max(0, ymax))
            x0, x1 = max(0, xmin), min(sheet_w, max(0, xmax))
            crop_bgr = img_bgr[y0:y1, x0:x1]
            local_cap_mask: Optional[np.ndarray] = None
            if mask is not None and isinstance(mask, np.ndarray) and mask.size > 0:
                if mask.shape[:2] == (sheet_h, sheet_w):
                    local_cap_mask = mask[y0:y1, x0:x1]
                else:
                    local_cap_mask = mask

            head_pheno_item: Optional[Dict[str, Any]] = None
            if crop_bgr.size > 0:
                cap_patch_path = self.capitula_patches_dir / f"{catalog_number}_cap{cap_idx}.jpg"
                try:
                    cv2.imwrite(str(cap_patch_path), crop_bgr, [int(cv2.IMWRITE_JPEG_QUALITY), 85])
                except TypeError:
                    cv2.imwrite(str(cap_patch_path), crop_bgr)

                if crop_bgr.shape[0] >= 5 and crop_bgr.shape[1] >= 5:
                    crop_rgb = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGB)
                    head_pheno_item = classify_single_capitulum(
                        crop_rgb,
                        min_pappus_luminance=_min_pappus_lum,
                        max_pappus_saturation=_max_pappus_sat,
                        min_pappus_area_ratio=_min_pappus_area,
                        min_floret_yellow_saturation=_min_floret_sat,
                        min_floret_yellow_ratio=_min_floret_rat,
                    )
                    head_pheno_results.append(head_pheno_item)

            metrics = compute_capitulum_metrics(
                local_cap_mask if local_cap_mask is not None else mask,
                bbox,
                pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
                crop_bgr=crop_bgr if crop_bgr.size > 0 else None,
                class_id=int(class_id),
                is_calibrated_mm=is_calibrated_mm,
            )
            if metrics is not None:
                metrics["score"] = score
                metrics["class_id"] = int(class_id)
                valid_capitula.append(metrics)

            raw_cap_items.append({
                "bbox": bbox,
                "center_xy": (0.5 * float(xmin + xmax), 0.5 * float(ymin + ymax)),
                "mask": local_cap_mask if local_cap_mask is not None else mask,
                "crop_bgr": crop_bgr if crop_bgr.size > 0 else None,
                "score": float(score),
                "class_id": int(class_id),
                "metrics": metrics,
                "head_pheno": head_pheno_item,
            })

        def _summarize_repro_subset(cap_metrics_sub: List[Dict[str, Any]], pid_val: int = 0) -> Dict[str, Any]:
            c_cnt = len(cap_metrics_sub)
            if c_cnt > 0:
                mh_px = float(np.median([c["involucre_height_px"] for c in cap_metrics_sub]))
                mw_px = float(np.median([c["involucre_width_px"] for c in cap_metrics_sub]))
                hmm_l = [c["involucre_height_mm"] for c in cap_metrics_sub if not np.isnan(c["involucre_height_mm"])]
                wmm_l = [c["involucre_width_mm"] for c in cap_metrics_sub if not np.isnan(c["involucre_width_mm"])]
                mh_mm = float(np.median(hmm_l)) if hmm_l else np.nan
                mw_mm = float(np.median(wmm_l)) if wmm_l else np.nan
                m_ar = float(np.median([c["AR_inv"] for c in cap_metrics_sub]))
            else:
                mh_px = mw_px = mh_mm = mw_mm = m_ar = np.nan
            return {
                "catalogNumber": catalog_number,
                "plant_individual_id": int(pid_val),
                "capitula_count": int(c_cnt),
                "involucre_height_px": round(mh_px, 2) if not np.isnan(mh_px) else np.nan,
                "involucre_width_px": round(mw_px, 2) if not np.isnan(mw_px) else np.nan,
                "involucre_height_mm": round(mh_mm, 3) if not np.isnan(mh_mm) else np.nan,
                "involucre_width_mm": round(mw_mm, 3) if not np.isnan(mw_mm) else np.nan,
                "AR_inv": round(m_ar, 4) if not np.isnan(m_ar) else np.nan,
                "plant_height_cm": np.nan,
                "scape_to_basal_leaf_ratio": np.nan,
                "is_bent_stem": False,
            }

        _pheno_rec = dict(classify_voucher_phenology(
            catalog_number,
            head_pheno_results,
            min_anthesis_heads_for_flowering=_min_anth_heads,
        ))
        _pheno_rec["plant_individual_id"] = 0
        capitula_count = len(valid_capitula)
        repro_summary = _summarize_repro_subset(valid_capitula, pid_val=0)

        if not leaf_detections:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "no_leaves_detected",
                "details": f"Model produced 0 candidate leaf segmentations ({capitula_count} capitula detected)",
                "pixels_per_mm": round(float(pixels_per_mm), 4) if is_calibrated_mm else np.nan,
                "is_calibrated_mm": bool(is_calibrated_mm),
            }
            return _finish([], failure_record, repro_rec=repro_summary, pheno_rec=_pheno_rec)

        # Build detected leaf instances
        instances: List[DetectedInstance] = []
        for idx, (bbox, mask, score, class_id, crop_offset, mask_origin, already_trimmed, trim_meta, _, _) in enumerate(leaf_detections, 1):
            _pet_mask = trim_meta.get("petiole_mask") if isinstance(trim_meta, dict) else None
            _hole_mask = trim_meta.get("hole_mask") if isinstance(trim_meta, dict) else None
            _kpts_31 = trim_meta.get("keypoints_31") if isinstance(trim_meta, dict) else None
            _pw_px_raw = trim_meta.get("petiole_width_px") if isinstance(trim_meta, dict) else None
            _pw_px = (
                float(_pw_px_raw)
                if (_pw_px_raw is not None and np.isfinite(float(_pw_px_raw)) and float(_pw_px_raw) > 0.0)
                else None
            )
            _pw_mm = (
                float(_pw_px) / float(pixels_per_mm)
                if (_pw_px is not None and is_calibrated_mm and pixels_per_mm and float(pixels_per_mm) > 0.0)
                else None
            )
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
                mask_origin=mask_origin,
                trim_meta=trim_meta,
                petiole_mask=_pet_mask,
                hole_mask=_hole_mask,
                keypoints_31=_kpts_31,
                petiole_width_px=_pw_px,
                petiole_width_mm=_pw_mm,
            ))

        # 4. Vectorized Centroid Spatial Clustering (DBSCAN)
        instances = cluster_plant_individuals(instances, sheet_w, sheet_h)

        # 5. Basal Rosette Zone Statistics per plant_individual_id
        def _inst_sheet_geom(inst_obj: DetectedInstance) -> Tuple[float, float, float, float, float]:
            ox_i, oy_i = inst_obj.mask_origin
            m_u8 = inst_obj.mask if inst_obj.mask.dtype == np.uint8 else ((inst_obj.mask > 0).astype(np.uint8) * 255)
            bx_i, by_i, bw_i, bh_i = cv2.boundingRect(m_u8)
            if bw_i > 0 and bh_i > 0:
                sub_i = m_u8[by_i : by_i + bh_i, bx_i : bx_i + bw_i]
                mom_i = cv2.moments(sub_i, binaryImage=True)
                area_i = float(mom_i["m00"])
            else:
                sub_i = None
                mom_i = None
                area_i = 0.0
            if sub_i is not None and mom_i is not None and area_i > 0.0:
                cx_i = float(mom_i["m10"] / area_i) + float(bx_i + ox_i)
                cy_i = float(mom_i["m01"] / area_i) + float(by_i + oy_i)
                cnts_i, _ = cv2.findContours(sub_i, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                if cnts_i:
                    main_c = max(cnts_i, key=cv2.contourArea)
                    rw, rh = cv2.minAreaRect(main_c)[1]
                    ar_i = float(max(rw, rh)) / max(float(min(rw, rh)), 1.0)
                    hull_a = float(cv2.contourArea(cv2.convexHull(main_c)))
                    sol_i = float(area_i / hull_a) if hull_a > 0 else 1.0
                else:
                    ar_i = float(max(bw_i, bh_i)) / max(float(min(bw_i, bh_i)), 1.0)
                    sol_i = 1.0
            else:
                cx_i = 0.5 * float(inst_obj.bbox[1] + inst_obj.bbox[3])
                cy_i = 0.5 * float(inst_obj.bbox[0] + inst_obj.bbox[2])
                bh_b = max(1.0, float(inst_obj.bbox[2] - inst_obj.bbox[0]))
                bw_b = max(1.0, float(inst_obj.bbox[3] - inst_obj.bbox[1]))
                area_i = bh_b * bw_b
                ar_i = max(bh_b, bw_b) / max(min(bh_b, bw_b), 1.0)
                sol_i = 1.0
            return cx_i, cy_i, area_i, ar_i, sol_i

        inst_geom_map: Dict[int, Tuple[float, float, float, float, float]] = {
            inst.leaf_id: _inst_sheet_geom(inst) for inst in instances
        }
        _lower_frac_cfg = float(
            getattr(getattr(self.config, "segmentation", None), "rosette_plant_bbox_lower_fraction", 0.60) or 0.60
        )
        _basal_y_quantile = float(np.clip(1.0 - _lower_frac_cfg, 0.05, 0.90))
        all_cx = np.array([inst_geom_map[i.leaf_id][0] for i in instances], dtype=np.float64)
        all_cy = np.array([inst_geom_map[i.leaf_id][1] for i in instances], dtype=np.float64)
        all_area = np.array([inst_geom_map[i.leaf_id][2] for i in instances], dtype=np.float64)
        sheet_b_mask = (
            all_cy >= float(np.quantile(all_cy, _basal_y_quantile))
            if len(all_cy) > 0
            else np.array([], dtype=bool)
        )
        if len(sheet_b_mask) > 0 and not np.any(sheet_b_mask):
            sheet_b_mask = np.ones(len(all_cy), dtype=bool)
        sheet_basal_cx = float(np.median(all_cx[sheet_b_mask])) if len(all_cx) > 0 else (0.5 * sheet_w)
        sheet_basal_cy = float(np.median(all_cy[sheet_b_mask])) if len(all_cy) > 0 else (0.75 * sheet_h)
        sheet_med_area = float(np.median(all_area[sheet_b_mask])) if len(all_area) > 0 else 1.0
        sheet_radii = np.hypot(all_cx[sheet_b_mask] - sheet_basal_cx, all_cy[sheet_b_mask] - sheet_basal_cy) if len(all_cx) > 0 else np.array([0.0])
        sheet_rosette_rad = max(float(np.median(sheet_radii)), 0.04 * float(min(sheet_h, sheet_w))) if len(instances) >= 2 else 0.0

        cap_cy: Optional[float] = (
            float(np.median([0.5 * float(b[0] + b[2]) for b, _, _, _ in capitulum_detections]))
            if capitulum_detections
            else None
        )

        clusters_by_pid: Dict[int, List[DetectedInstance]] = defaultdict(list)
        for inst in instances:
            clusters_by_pid[inst.plant_individual_id].append(inst)
        rosette_centroids_by_pid: Dict[int, Optional[Tuple[float, float]]] = {}
        basal_stats_by_pid: Dict[int, Tuple[float, float, float, float]] = {}

        for pid, pid_instances in sorted(clusters_by_pid.items()):
            p_cx = np.array([inst_geom_map[i.leaf_id][0] for i in pid_instances], dtype=np.float64)
            p_cy = np.array([inst_geom_map[i.leaf_id][1] for i in pid_instances], dtype=np.float64)
            p_area = np.array([inst_geom_map[i.leaf_id][2] for i in pid_instances], dtype=np.float64)
            # Lower configured vertical quantile of detected leaf centroids (in image coords, larger y = lower on sheet)
            b_mask = p_cy >= float(np.quantile(p_cy, _basal_y_quantile))
            if not np.any(b_mask):
                b_mask = np.ones(len(pid_instances), dtype=bool)
            basal_cx_med = float(np.median(p_cx[b_mask]))
            basal_cy_med = float(np.median(p_cy[b_mask]))
            med_area = float(np.median(p_area[b_mask]))
            radii = np.hypot(p_cx[b_mask] - basal_cx_med, p_cy[b_mask] - basal_cy_med)
            med_rosette_radius = (
                max(float(np.median(radii)), 0.04 * float(min(sheet_h, sheet_w)))
                if len(pid_instances) >= 2
                else 0.0
            )
            if (
                len(pid_instances) == 1
                and len(instances) >= 2
                and (
                    p_cy[0] < (sheet_basal_cy - 0.15 * float(sheet_h))
                    or p_area[0] < (0.55 * sheet_med_area)
                )
            ):
                basal_stats_by_pid[pid] = (sheet_basal_cx, sheet_basal_cy, sheet_med_area, sheet_rosette_rad)
            else:
                basal_stats_by_pid[pid] = (basal_cx_med, basal_cy_med, med_area, med_rosette_radius)

            if len(pid_instances) >= 2:
                rosette_centroids_by_pid[pid] = (basal_cx_med, basal_cy_med)
            else:
                rosette_centroids_by_pid[pid] = None

        extracted_records: List[Dict[str, Any]] = []
        candidate_extractions: List[
            Tuple[
                DetectedInstance,
                Optional[Tuple[int, int]],
                Optional[Tuple[int, int]],
                Optional[float],
                np.ndarray,
                Dict[str, Tuple[float, float]],
                bool,
            ]
        ] = []

        # 6. Read-Only LeafMachine3 Orientation, Native Morphology Readout, and Persistence
        px_per_cm_sheet = (
            resolve_px_per_cm(
                pixels_per_mm=pixels_per_mm,
                default_dpi=DEFAULT_HERBARIUM_DPI,
                sheet_shape=(sheet_h, sheet_w),
            )
            if is_calibrated_mm
            else None
        )

        _pid_ymins: Dict[int, List[float]] = {}
        _pid_ymaxs: Dict[int, List[float]] = {}
        _pid_xs: Dict[int, List[float]] = {}
        for _inst_b in instances:
            _p_id = int(_inst_b.plant_individual_id)
            _pid_ymins.setdefault(_p_id, []).append(float(_inst_b.bbox[0]))
            _pid_ymaxs.setdefault(_p_id, []).append(float(_inst_b.bbox[2]))
            _pid_xs.setdefault(_p_id, []).extend([float(_inst_b.bbox[1]), float(_inst_b.bbox[3])])

        if len(_pid_ymins) == 1:
            _sole_pid = next(iter(_pid_ymins.keys()))
            for _cap_b, _, _, _ in capitulum_detections:
                _pid_ymins[_sole_pid].append(float(_cap_b[0]))
                _pid_ymaxs[_sole_pid].append(float(_cap_b[2]))
            for _root_b in root_detections:
                _pid_ymins[_sole_pid].append(float(_root_b[0]))
                _pid_ymaxs[_sole_pid].append(float(_root_b[2]))
            for _uc_x, _uc_y in upper_cauline_pts:
                _pid_ymins[_sole_pid].append(float(_uc_y))
                _pid_ymaxs[_sole_pid].append(float(_uc_y))
        elif len(_pid_ymins) > 1:
            _pid_centers_xy = {
                _k: (float(np.median(_pid_xs[_k])), float(np.median(_pid_ymins[_k] + _pid_ymaxs[_k])))
                for _k in _pid_ymins.keys()
            }
            for _cap_b, _, _, _ in capitulum_detections:
                _cx_b = 0.5 * (float(_cap_b[1]) + float(_cap_b[3]))
                _cy_b = 0.5 * (float(_cap_b[0]) + float(_cap_b[2]))
                _best_k = min(
                    _pid_centers_xy.keys(),
                    key=lambda k: math.hypot(_cx_b - _pid_centers_xy[k][0], _cy_b - _pid_centers_xy[k][1]),
                )
                if math.hypot(_cx_b - _pid_centers_xy[_best_k][0], _cy_b - _pid_centers_xy[_best_k][1]) <= 0.35 * max(sheet_h, sheet_w):
                    _pid_ymins[_best_k].append(float(_cap_b[0]))
                    _pid_ymaxs[_best_k].append(float(_cap_b[2]))
            for _root_b in root_detections:
                _rx_b = 0.5 * (float(_root_b[1]) + float(_root_b[3]))
                _ry_b = 0.5 * (float(_root_b[0]) + float(_root_b[2]))
                _best_k = min(
                    _pid_centers_xy.keys(),
                    key=lambda k: math.hypot(_rx_b - _pid_centers_xy[k][0], _ry_b - _pid_centers_xy[k][1]),
                )
                if math.hypot(_rx_b - _pid_centers_xy[_best_k][0], _ry_b - _pid_centers_xy[_best_k][1]) <= 0.35 * max(sheet_h, sheet_w):
                    _pid_ymins[_best_k].append(float(_root_b[0]))
                    _pid_ymaxs[_best_k].append(float(_root_b[2]))

        inst_cy_rel_plant_bbox: Dict[int, float] = {}
        for _inst_b in instances:
            _p_id = int(_inst_b.plant_individual_id)
            _ymin_p = float(min(_pid_ymins.get(_p_id, [0.0])))
            _ymax_p = float(max(_pid_ymaxs.get(_p_id, [float(sheet_h)])))
            _span_p = max(1.0, _ymax_p - _ymin_p)
            inst_cy_rel_plant_bbox[_inst_b.leaf_id] = max(
                0.0,
                min(1.0, (inst_geom_map[_inst_b.leaf_id][1] - _ymin_p) / _span_p),
            )

        _raw_px_by_leaf_id: Dict[int, float] = {}
        for _inst_pre in instances:
            _rpx = float(np.count_nonzero(_inst_pre.mask))
            _raw_px_by_leaf_id[_inst_pre.leaf_id] = _rpx
            _inst_pre.area_mm2 = compute_area_mm2(_rpx, pixels_per_mm=pixels_per_mm, is_calibrated_mm=is_calibrated_mm)

        # Pre-compute per-individual crown / root / rosette hub hints for 180-deg tie-breaking
        _crown_hints_sheet: Dict[int, Optional[Tuple[float, float]]] = {}
        for _pid_h, _rc_h in rosette_centroids_by_pid.items():
            _best_root_pt: Optional[Tuple[float, float]] = None
            if root_detections:
                _ref_xy = _rc_h if _rc_h is not None else (
                    basal_stats_by_pid[_pid_h][0],
                    basal_stats_by_pid[_pid_h][1],
                )
                _r_cands = [
                    (0.5 * float(rb[1] + rb[3]), 0.5 * float(rb[0] + rb[2]))
                    for rb in root_detections
                ]
                _nearest_r = min(_r_cands, key=lambda rxy: math.hypot(rxy[0] - _ref_xy[0], rxy[1] - _ref_xy[1]))
                if len(rosette_centroids_by_pid) == 1 or math.hypot(
                    _nearest_r[0] - _ref_xy[0], _nearest_r[1] - _ref_xy[1]
                ) <= 0.35 * max(sheet_h, sheet_w):
                    _best_root_pt = _nearest_r
            _crown_hints_sheet[_pid_h] = _best_root_pt if _best_root_pt is not None else _rc_h

        for inst in instances:
            if int(np.count_nonzero(inst.mask)) <= 0:
                continue
            _ox_l, _oy_l = inst.mask_origin
            _pid_hint_s = _crown_hints_sheet.get(int(inst.plant_individual_id))
            _crown_local: Optional[Tuple[float, float]] = (
                (float(_pid_hint_s[0]) - float(_ox_l), float(_pid_hint_s[1]) - float(_oy_l))
                if _pid_hint_s is not None
                else None
            )
            p_base_use: Optional[Tuple[int, int]] = None
            p_apex_use: Optional[Tuple[int, int]] = None
            if (
                isinstance(inst.keypoints_31, dict)
                and "lamina_base" in inst.keypoints_31
                and "lamina_tip" in inst.keypoints_31
                and len(inst.keypoints_31["lamina_base"]) >= 2
                and len(inst.keypoints_31["lamina_tip"]) >= 2
            ):
                _lb_inst = inst.keypoints_31["lamina_base"]
                _lt_inst = inst.keypoints_31["lamina_tip"]
                p_base_use = (int(round(float(_lb_inst[0]))), int(round(float(_lb_inst[1]))))
                p_apex_use = (int(round(float(_lt_inst[0]))), int(round(float(_lt_inst[1]))))
            elif isinstance(inst.trim_meta, dict):
                if inst.trim_meta.get("p_base") is not None:
                    p_base_use = inst.trim_meta["p_base"]
                if inst.trim_meta.get("p_apex") is not None:
                    p_apex_use = inst.trim_meta["p_apex"]

            if p_base_use is None or p_apex_use is None:
                _m_u8 = (inst.mask > 0).astype(np.uint8) * 255
                _cnts_fb, _ = cv2.findContours(_m_u8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                if _cnts_fb:
                    _pts_fb = max(_cnts_fb, key=cv2.contourArea).reshape(-1, 2)
                    if p_base_use is None:
                        if inst.petiole_mask is not None and int(np.count_nonzero(inst.petiole_mask)) > 0:
                            _pys, _pxs = np.where(inst.petiole_mask > 0)
                            _pcx, _pcy = float(np.mean(_pxs)), float(np.mean(_pys))
                            _bi = int(np.argmin(np.hypot(_pts_fb[:, 0] - _pcx, _pts_fb[:, 1] - _pcy)))
                        elif _crown_local is not None:
                            _bi = int(np.argmin(np.hypot(_pts_fb[:, 0] - _crown_local[0], _pts_fb[:, 1] - _crown_local[1])))
                        else:
                            _bi = int(np.argmax(_pts_fb[:, 1]))
                        p_base_use = (int(_pts_fb[_bi, 0]), int(_pts_fb[_bi, 1]))
                    _ai = int(np.argmax(np.hypot(_pts_fb[:, 0] - p_base_use[0], _pts_fb[:, 1] - p_base_use[1])))
                    p_apex_use = p_apex_use or (int(_pts_fb[_ai, 0]), int(_pts_fb[_ai, 1]))
                else:
                    p_base_use = p_base_use or (0, 0)
                    p_apex_use = p_apex_use or (0, 0)

            inst.cf_source = cf_source

            lm3_bm_raw = orient_and_score_leaf_lm3(
                inst.mask,
                keypoints_31=inst.keypoints_31,
                p_base=p_base_use,
                p_apex=p_apex_use,
                offset_xy=(0, 0),
                hole_mask=inst.hole_mask,
                petiole_mask=inst.petiole_mask,
                crown_hint_xy=_crown_local,
            )
            angle_cw_val = float(lm3_bm_raw.get("angle_cw", 0.0)) if lm3_bm_raw.get("success") else None
            # Record sheet-space midrib angle reconciled with LM3 orientation (angle_cw rotates midrib onto -90 deg)
            if angle_cw_val is not None:
                _sheet_deg = -90.0 - float(angle_cw_val)
                inst.midrib_angle_deg = float(((_sheet_deg + 180.0) % 360.0) - 180.0)
            else:
                inst.midrib_angle_deg = float(
                    math.degrees(math.atan2(float(p_apex_use[1] - p_base_use[1]), float(p_apex_use[0] - p_base_use[0])))
                )

            lm3_morph = compute_lm3_shape_metrics(inst.mask, angle_cw=angle_cw_val)
            inst.solidity = float(lm3_morph["solidity"])
            inst.circularity = float(lm3_morph["circularity"])
            inst.extent = float(lm3_morph["extent"])
            _bm_lw = lm3_bm_raw.get("aspect_ratio", np.nan)
            _bm_bl_px = lm3_bm_raw.get("blade_length_px", None)
            _bm_bw_px = lm3_bm_raw.get("blade_width_px", None)
            if lm3_bm_raw.get("success") and _bm_lw is not None and np.isfinite(float(_bm_lw)) and float(_bm_lw) > 0.0:
                inst.aspect_ratio = float(_bm_lw)
                if _bm_bl_px is not None and np.isfinite(float(_bm_bl_px)) and float(_bm_bl_px) > 0.0:
                    inst.blade_length_px = float(_bm_bl_px)
                if _bm_bw_px is not None and np.isfinite(float(_bm_bw_px)) and float(_bm_bw_px) > 0.0:
                    inst.blade_width_px = float(_bm_bw_px)
            else:
                inst.aspect_ratio = float(lm3_morph["aspect_ratio"])
            inst.archetype_score = float(lm3_bm_raw.get("archetype_score", 0.0))
            inst.bilateral_dice = float(lm3_bm_raw.get("bilateral_dice", 1.0))  # type: ignore[attr-defined]
            inst.si_a = float(lm3_bm_raw.get("si_a", 0.0))  # type: ignore[attr-defined]
            inst.lm3_gates_pass = bool(lm3_bm_raw.get("gates_pass", True))  # type: ignore[attr-defined]

            _ori_m = lm3_bm_raw.get("oriented_mask")
            if isinstance(_ori_m, np.ndarray) and _ori_m.size > 0 and int(np.count_nonzero(_ori_m)) > 0:
                _persisted_mask_u8 = (_ori_m > 0).astype(np.uint8) * 255
                _ori_pts_map: Dict[str, Tuple[float, float]] = dict(lm3_bm_raw.get("oriented_points") or {})
                _is_upright_mask = True
            else:
                _persisted_mask_u8 = (inst.mask > 0).astype(np.uint8) * 255
                if (
                    bool(getattr(self, "subtract_holes", True))
                    and inst.hole_mask is not None
                    and isinstance(inst.hole_mask, np.ndarray)
                    and inst.hole_mask.shape[:2] == _persisted_mask_u8.shape[:2]
                ):
                    _persisted_mask_u8 = np.where((_persisted_mask_u8 > 0) & (inst.hole_mask == 0), 255, 0).astype(np.uint8)
                _ori_pts_map = {}
                _is_upright_mask = False

            _n_cc_lbls, _, _cc_stats, _ = cv2.connectedComponentsWithStats(
                (inst.mask > 0).astype(np.uint8), connectivity=8
            )
            _n_sig_comps = (
                sum(1 for cid in range(1, _n_cc_lbls) if int(_cc_stats[cid, cv2.CC_STAT_AREA]) >= 32)
                if _n_cc_lbls > 1
                else 1
            )
            inst.n_lamina_components = max(1, _n_sig_comps)  # type: ignore[attr-defined]
            candidate_extractions.append(
                (inst, p_base_use, p_apex_use, angle_cw_val, _persisted_mask_u8, _ori_pts_map, _is_upright_mask)
            )

        pose_inputs: List[Dict[str, Any]] = []
        for _inst, _pb, _pa, _ang_cw, _pmask, _opts, _is_up in candidate_extractions:
            _ox, _oy = _inst.mask_origin
            _pb_use = _pb or (0, 0)
            _pa_use = _pa or (0, 0)
            _pb_sheet = (float(_pb_use[0] + _ox), float(_pb_use[1] + _oy))
            _pa_sheet = (float(_pa_use[0] + _ox), float(_pa_use[1] + _oy))
            _trim_info = _inst.trim_meta or {}
            _geod_px = float(
                _trim_info.get("petiole_length_px", 0.0)
                or _trim_info.get("geodesic_petiole_skeleton_px", 0.0)
                or _trim_info.get("trim_distance_px", 0.0)
                or 0.0
            )
            pose_inputs.append({
                "leaf_id": _inst.leaf_id,
                "plant_individual_id": _inst.plant_individual_id,
                "p_base": _pb_sheet,
                "p_apex": _pa_sheet,
                "p_base_sheet": _pb_sheet,
                "p_apex_sheet": _pa_sheet,
                "mask_origin": (0.0, 0.0),
                "bbox": _inst.bbox,
                "mask": _inst.mask,
                "area_px": _raw_px_by_leaf_id.get(_inst.leaf_id, float(np.count_nonzero(_inst.mask))),
                "area_mm2": _inst.area_mm2,
                "blade_length_px": _inst.blade_length_px,
                "blade_width_px": _inst.blade_width_px,
                "petiole_width_px": _inst.petiole_width_px,
                "petiole_width_mm": _inst.petiole_width_mm,
                "aspect_ratio": _inst.aspect_ratio,
                "solidity": _inst.solidity,
                "geodesic_petiole_px": _geod_px,
                "geodesic_petiole_skeleton_px": _geod_px,
                "cy_rel_plant_bbox": inst_cy_rel_plant_bbox.get(_inst.leaf_id, np.nan),
                "keypoints_31": _inst.keypoints_31,
                "petiole_mask": _inst.petiole_mask,
                "local_mask_origin": (_ox, _oy),
            })

        crowns_by_pid = estimate_rosette_crown_point(
            pose_inputs,
            root_boxes=root_detections,
            rosette_centroids=rosette_centroids_by_pid,
        )
        _spec_mask_sheet: Optional[np.ndarray] = None
        if (
            pose_inputs
            and len(pose_inputs) >= 2
            and is_high_dpi
            and getattr(self, "specimen_segmenter", None) is not None
            and getattr(self.specimen_segmenter, "model", None) is not None
        ):
            try:
                _spec_mask_sheet = self.specimen_segmenter.predict_mask(img_bgr, return_sheet_size=True)
            except Exception:
                _spec_mask_sheet = None

        # Partition capitula across plant_individual_id so sterile individuals never inherit another plant's scape/involucre
        unique_leaf_pids = sorted({int(p["plant_individual_id"]) for p in pose_inputs})
        caps_by_pid_boxes: Dict[int, List[Tuple[int, int, int, int]]] = {pid: [] for pid in unique_leaf_pids}
        cap_items_by_pid: Dict[int, List[Dict[str, Any]]] = {pid: [] for pid in unique_leaf_pids}

        if len(unique_leaf_pids) <= 1:
            only_pid = unique_leaf_pids[0] if unique_leaf_pids else 0
            cr_only = crowns_by_pid.get(only_pid) if isinstance(crowns_by_pid, dict) else None
            refreshed_valid: List[Dict[str, Any]] = []
            for c_it in raw_cap_items:
                if cr_only is not None:
                    ped_dir = (c_it["center_xy"][0] - float(cr_only[0]), c_it["center_xy"][1] - float(cr_only[1]))
                    m_ref = compute_capitulum_metrics(
                        c_it["mask"],
                        c_it["bbox"],
                        pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
                        crop_bgr=c_it["crop_bgr"],
                        class_id=c_it["class_id"],
                        peduncle_dir_xy=ped_dir,
                        is_calibrated_mm=is_calibrated_mm,
                    )
                    if m_ref is not None:
                        m_ref["score"] = c_it["score"]
                        m_ref["class_id"] = c_it["class_id"]
                    c_it["metrics"] = m_ref
                if c_it["metrics"] is not None:
                    refreshed_valid.append(c_it["metrics"])
                caps_by_pid_boxes.setdefault(only_pid, []).append(c_it["bbox"])
                cap_items_by_pid.setdefault(only_pid, []).append(c_it)
            valid_capitula = refreshed_valid
            capitula_count = len(valid_capitula)
            repro_summary = _summarize_repro_subset(valid_capitula, pid_val=only_pid)
            _pheno_rec["plant_individual_id"] = int(only_pid)
        else:
            cc_labels_sheet: Optional[np.ndarray] = None
            pid_cc_sets: Dict[int, Set[int]] = {pid: set() for pid in unique_leaf_pids}
            if isinstance(_spec_mask_sheet, np.ndarray) and _spec_mask_sheet.shape[:2] == (sheet_h, sheet_w) and np.any(_spec_mask_sheet > 0):
                _, cc_labels_sheet = cv2.connectedComponents((_spec_mask_sheet > 0).astype(np.uint8), connectivity=8)
                for p_in in pose_inputs:
                    px_i = int(np.clip(round(p_in["p_base"][0]), 0, sheet_w - 1))
                    py_i = int(np.clip(round(p_in["p_base"][1]), 0, sheet_h - 1))
                    lbl_p = int(cc_labels_sheet[py_i, px_i])
                    if lbl_p > 0:
                        pid_cc_sets[int(p_in["plant_individual_id"])].add(lbl_p)
            max_link_dist = 0.55 * float(max(sheet_h, sheet_w))
            refreshed_valid = []
            for c_it in raw_cap_items:
                cx_c, cy_c = c_it["center_xy"]
                assigned_pid: Optional[int] = None
                if cc_labels_sheet is not None:
                    cxi = int(np.clip(round(cx_c), 0, sheet_w - 1))
                    cyi = int(np.clip(round(cy_c), 0, sheet_h - 1))
                    lbl_c = int(cc_labels_sheet[cyi, cxi])
                    if lbl_c > 0:
                        matching_pids = [pid for pid, s_lbls in pid_cc_sets.items() if lbl_c in s_lbls]
                        if len(matching_pids) == 1:
                            assigned_pid = matching_pids[0]
                if assigned_pid is None:
                    best_d = float("inf")
                    for pid in unique_leaf_pids:
                        pid_pts = [p["p_base"] for p in pose_inputs if int(p["plant_individual_id"]) == pid] + [
                            p["p_apex"] for p in pose_inputs if int(p["plant_individual_id"]) == pid
                        ]
                        cr_p = crowns_by_pid.get(pid) if isinstance(crowns_by_pid, dict) else None
                        if cr_p is not None:
                            pid_pts.append((float(cr_p[0]), float(cr_p[1])))
                        if not pid_pts:
                            continue
                        min_leaf_d = min(math.hypot(cx_c - px_p, cy_c - py_p) for px_p, py_p in pid_pts)
                        # Penalize assigning an upper-sheet corymb to a rosette situated far laterally or above the capitulum
                        below_pen = 1.35 if (cr_p is not None and cy_c > float(cr_p[1]) + 40.0) else 1.0
                        eff_d = min_leaf_d * below_pen
                        if eff_d < best_d:
                            best_d = eff_d
                            assigned_pid = pid
                    if best_d > max_link_dist:
                        assigned_pid = None
                if assigned_pid is not None:
                    cr_ass = crowns_by_pid.get(assigned_pid) if isinstance(crowns_by_pid, dict) else None
                    if cr_ass is not None:
                        ped_dir = (cx_c - float(cr_ass[0]), cy_c - float(cr_ass[1]))
                        m_ref = compute_capitulum_metrics(
                            c_it["mask"],
                            c_it["bbox"],
                            pixels_per_mm=pixels_per_mm if is_calibrated_mm else None,
                            crop_bgr=c_it["crop_bgr"],
                            class_id=c_it["class_id"],
                            peduncle_dir_xy=ped_dir,
                            is_calibrated_mm=is_calibrated_mm,
                        )
                        if m_ref is not None:
                            m_ref["score"] = c_it["score"]
                            m_ref["class_id"] = c_it["class_id"]
                        c_it["metrics"] = m_ref
                    caps_by_pid_boxes[assigned_pid].append(c_it["bbox"])
                    cap_items_by_pid[assigned_pid].append(c_it)
                if c_it["metrics"] is not None:
                    refreshed_valid.append(c_it["metrics"])
            valid_capitula = refreshed_valid
            capitula_count = len(valid_capitula)
            repro_summary = _summarize_repro_subset(valid_capitula, pid_val=0)

        petiole_scape_res = approximate_petiole_and_scape_metrics(
            basal_leaves=pose_inputs,
            p_crown=crowns_by_pid,
            px_per_cm=px_per_cm_sheet,
            is_calibrated_mm=is_calibrated_mm,
            upper_cauline_detections=upper_cauline_pts,
            root_boxes=root_detections,
            capitula_detections=[b for b, _, _, _ in capitulum_detections],
            capitula_by_pid=caps_by_pid_boxes,
        )
        repro_summary["plant_height_cm"] = petiole_scape_res.get("plant_height_cm", np.nan)
        repro_summary["scape_to_basal_leaf_ratio"] = petiole_scape_res.get("scape_to_basal_leaf_ratio", np.nan)
        repro_summary["is_bent_stem"] = bool(petiole_scape_res.get("is_bent_stem", False))

        leaf_zones_by_id = classify_leaf_zone_anatomical(
            pose_inputs,
            crowns_by_pid=crowns_by_pid if isinstance(crowns_by_pid, dict) else None,
            rosette_centroids_by_pid=rosette_centroids_by_pid,
            basal_stats_by_pid=basal_stats_by_pid,
            capitula_detections=[b for b, _, _, _ in capitulum_detections],
            root_boxes=root_detections,
            is_bent_stem=bool(repro_summary.get("is_bent_stem", False)),
            sheet_shape=(sheet_h, sheet_w),
            specimen_mask=_spec_mask_sheet,
        )
        has_cauline_zone = False
        for p_in in pose_inputs:
            _lz = str((leaf_zones_by_id.get(int(p_in["leaf_id"])) or {}).get("leaf_zone", "basal") or "basal").lower()
            p_in["leaf_zone"] = _lz
            if _lz == "cauline":
                has_cauline_zone = True
        if has_cauline_zone:
            petiole_scape_res = approximate_petiole_and_scape_metrics(
                basal_leaves=pose_inputs,
                p_crown=crowns_by_pid,
                px_per_cm=px_per_cm_sheet,
                is_calibrated_mm=is_calibrated_mm,
                upper_cauline_detections=upper_cauline_pts,
                root_boxes=root_detections,
                capitula_detections=[b for b, _, _, _ in capitulum_detections],
                capitula_by_pid=caps_by_pid_boxes,
            )
            if np.isfinite(petiole_scape_res.get("scape_to_basal_leaf_ratio", np.nan)):
                repro_summary["scape_to_basal_leaf_ratio"] = petiole_scape_res.get("scape_to_basal_leaf_ratio", np.nan)
        by_plant_petiole = petiole_scape_res.get("by_plant", {}) or {}
        per_leaf_petiole: Dict[int, Dict[str, Any]] = {
            int(r["leaf_id"]): r
            for r in (petiole_scape_res.get("leaf_metrics") or [])
            if isinstance(r, dict) and r.get("leaf_id") is not None
        }
        _centroids_out = {
            int(k): (float(v[0]), float(v[1]))
            for k, v in rosette_centroids_by_pid.items()
            if v is not None
        }
        _crowns_out: Dict[int, Tuple[float, float]] = (
            {
                int(k): (float(v[0]), float(v[1]))
                for k, v in crowns_by_pid.items()
                if isinstance(v, (tuple, list)) and len(v) >= 2
            }
            if isinstance(crowns_by_pid, dict)
            else {}
        )
        _petiole_out: Dict[int, Dict[str, Any]] = (
            {
                int(k): dict(v)
                for k, v in by_plant_petiole.items()
                if isinstance(v, dict)
            }
            if isinstance(by_plant_petiole, dict)
            else {}
        )
        _repro_by_pid_out: Dict[int, Dict[str, Any]] = {}
        _pheno_by_pid_out: Dict[int, Dict[str, Any]] = {}
        for _pid_r in unique_leaf_pids:
            _p_caps = cap_items_by_pid.get(_pid_r) or []
            _p_valid_m = [it["metrics"] for it in _p_caps if it.get("metrics") is not None]
            _p_heads_ph = [it["head_pheno"] for it in _p_caps if it.get("head_pheno") is not None]
            _r_sub = _summarize_repro_subset(_p_valid_m, pid_val=_pid_r)
            _pet_sub = _petiole_out.get(_pid_r, {})
            if len(unique_leaf_pids) <= 1:
                _r_sub["plant_height_cm"] = repro_summary["plant_height_cm"]
                _r_sub["scape_to_basal_leaf_ratio"] = repro_summary["scape_to_basal_leaf_ratio"]
                _r_sub["is_bent_stem"] = repro_summary["is_bent_stem"]
            else:
                _r_sub["plant_height_cm"] = _pet_sub.get("plant_height_cm", np.nan)
                _r_sub["scape_to_basal_leaf_ratio"] = _pet_sub.get("scape_to_basal_leaf_ratio", np.nan)
                _r_sub["is_bent_stem"] = bool(_pet_sub.get("is_bent_stem", False))
            _repro_by_pid_out[_pid_r] = _r_sub
            _ph_sub = dict(classify_voucher_phenology(
                catalog_number,
                _p_heads_ph,
                min_anthesis_heads_for_flowering=_min_anth_heads,
            ))
            _ph_sub["plant_individual_id"] = int(_pid_r)
            _pheno_by_pid_out[_pid_r] = _ph_sub

        _contours_by_leaf: Dict[int, List[List[float]]] = {}

        for inst, p_base_pt, p_apex_pt, angle_cw, persisted_mask_u8, ori_pts_map, is_upright_mask in candidate_extractions:
            inst.mask_path = save_leaf_mask_png(persisted_mask_u8, catalog_number, inst.leaf_id, self.output_dir)
            if is_upright_mask:
                _ori_base = ori_pts_map.get("lamina_base") or ori_pts_map.get("base_center")
                _ori_apex = ori_pts_map.get("lamina_tip") or ori_pts_map.get("apex_center")
                c_res = export_standardized_contour(
                    persisted_mask_u8,
                    catalog_number,
                    inst.leaf_id,
                    self.output_dir,
                    petiole_attachment_pt=_ori_base,
                    p_apex_pt=_ori_apex,
                    return_coords=True,
                    already_oriented=True,
                )
            else:
                c_res = export_standardized_contour(
                    persisted_mask_u8,
                    catalog_number,
                    inst.leaf_id,
                    self.output_dir,
                    petiole_attachment_pt=p_base_pt,
                    p_apex_pt=p_apex_pt,
                    hole_mask=inst.hole_mask,
                    return_coords=True,
                    keypoints_31=inst.keypoints_31,
                )
            if isinstance(c_res, tuple):
                inst.contour_path, c_xy = c_res
                if c_xy:
                    _contours_by_leaf[int(inst.leaf_id)] = c_xy
            else:
                inst.contour_path = c_res

            if getattr(self, "ect_dir", None) is not None and getattr(self, "ect_enabled", True):
                ect_out = self.ect_dir / f"{catalog_number}_leaf{inst.leaf_id}_ect.npy"
                try:
                    ect_mat = compute_leaf_ect_descriptor(
                        persisted_mask_u8,
                        num_dirs=getattr(self, "ect_num_dirs", 360),
                        num_thresh=getattr(self, "ect_num_thetas", 32),
                        output_path=ect_out,
                        hole_mask=(None if is_upright_mask else inst.hole_mask),
                        petiole_mask=(None if is_upright_mask else inst.petiole_mask),
                        angle_cw=(None if is_upright_mask else angle_cw),
                        include_petiole=bool(getattr(self, "ect_include_petiole", False)),
                        include_holes=bool(getattr(self, "ect_include_holes", True)),
                        already_oriented=is_upright_mask,
                    )
                except TypeError:
                    ect_mat = compute_leaf_ect_descriptor(persisted_mask_u8, output_path=ect_out)
                if ect_mat is not None:
                    inst.ect_path = str(ect_out)

            pid_pet = by_plant_petiole.get(int(inst.plant_individual_id), petiole_scape_res)
            leaf_pet = per_leaf_petiole.get(int(inst.leaf_id), {})
            _inst_area_px_val = _raw_px_by_leaf_id.get(inst.leaf_id, np.nan)
            _lz_info = leaf_zones_by_id.get(int(inst.leaf_id), {}) if isinstance(leaf_zones_by_id, dict) else {}
            _leaf_zone_str = str(_lz_info.get("leaf_zone", "basal") or "basal")
            _leaf_zone_conf = _lz_info.get("leaf_zone_confidence", np.nan)
            _dist_crown_norm = _lz_info.get("dist_from_crown_norm", np.nan)
            _pet_len_px = leaf_pet.get("petiole_length_px", np.nan)
            _blade_len_px = leaf_pet.get("blade_length_px", np.nan)
            _pet_len_mm = leaf_pet.get("petiole_length_mm", np.nan)
            _blade_len_mm = leaf_pet.get("blade_length_mm", np.nan)
            _pet_ratio = leaf_pet.get("petiole_to_blade_ratio", np.nan)
            _pid_repro_rec = _repro_by_pid_out.get(int(inst.plant_individual_id), repro_summary)
            _pet_w_px = (
                round(float(inst.petiole_width_px), 4)
                if inst.petiole_width_px is not None and np.isfinite(inst.petiole_width_px)
                else leaf_pet.get("petiole_width_px", np.nan)
            )
            _pet_w_mm = (
                round(float(inst.petiole_width_mm), 4)
                if inst.petiole_width_mm is not None and np.isfinite(inst.petiole_width_mm)
                else leaf_pet.get("petiole_width_mm", np.nan)
            )
            _blade_w_px = (
                round(float(inst.blade_width_px), 4)
                if inst.blade_width_px is not None and np.isfinite(inst.blade_width_px) and inst.blade_width_px > 0.0
                else (
                    round(float(_blade_len_px) / float(inst.aspect_ratio), 4)
                    if (
                        _blade_len_px is not None
                        and np.isfinite(float(_blade_len_px))
                        and float(_blade_len_px) > 0.0
                        and np.isfinite(float(inst.aspect_ratio))
                        and float(inst.aspect_ratio) > 1e-6
                    )
                    else np.nan
                )
            )
            _blade_w_mm = (
                round(float(_blade_w_px) / float(inst.pixels_per_mm), 3)
                if (
                    np.isfinite(_blade_w_px)
                    and inst.is_calibrated_mm
                    and inst.pixels_per_mm is not None
                    and float(inst.pixels_per_mm) > 0.0
                )
                else np.nan
            )

            extracted_records.append({
                "catalogNumber": catalog_number,
                "plant_individual_id": inst.plant_individual_id,
                "leaf_id": inst.leaf_id,
                "solidity": round(float(inst.solidity), 4),
                "circularity": round(float(inst.circularity), 4) if np.isfinite(inst.circularity) else np.nan,
                "extent": round(float(inst.extent), 4) if np.isfinite(inst.extent) else np.nan,
                "aspect_ratio": round(float(inst.aspect_ratio), 4) if np.isfinite(inst.aspect_ratio) else np.nan,
                "area_px": round(float(_inst_area_px_val), 1) if np.isfinite(_inst_area_px_val) else np.nan,
                "area_mm2": round(float(inst.area_mm2), 2) if np.isfinite(inst.area_mm2) else np.nan,
                "archetype_score": round(float(inst.archetype_score), 4),
                "bilateral_dice": round(float(getattr(inst, "bilateral_dice", 1.0)), 4),
                "si_a": round(float(getattr(inst, "si_a", 0.0)), 4),
                "lm3_gates_pass": bool(getattr(inst, "lm3_gates_pass", True)),
                "n_lamina_components": int(getattr(inst, "n_lamina_components", 1)),
                "midrib_angle_deg": round(float(inst.midrib_angle_deg), 2),
                "pixels_per_mm": round(float(inst.pixels_per_mm), 4) if inst.pixels_per_mm else np.nan,
                "is_calibrated_mm": bool(inst.is_calibrated_mm),
                "cf_source": inst.cf_source,
                "scale_source": (
                    "mp_estimated"
                    if (
                        inst.pixels_per_mm
                        and str(inst.cf_source).strip().lower() in {"predicted_from_megapixels", "mp_estimated"}
                    )
                    else ("ruler_calibrated" if (inst.pixels_per_mm and inst.is_calibrated_mm) else "uncalibrated")
                ),
                "mask_path": inst.mask_path,
                "mask_is_upright": bool(is_upright_mask),
                "contour_path": inst.contour_path,
                "ect_path": inst.ect_path,
                "petiole_length_px": _pet_len_px,
                "petiole_width_px": _pet_w_px,
                "blade_length_px": _blade_len_px,
                "blade_width_px": _blade_w_px,
                "petiole_length_mm": _pet_len_mm,
                "petiole_width_mm": _pet_w_mm,
                "blade_length_mm": _blade_len_mm,
                "blade_width_mm": _blade_w_mm,
                "petiole_to_blade_ratio": _pet_ratio,
                "capitula_count": _pid_repro_rec["capitula_count"],
                "involucre_height_px": _pid_repro_rec["involucre_height_px"],
                "involucre_width_px": _pid_repro_rec["involucre_width_px"],
                "involucre_height_mm": _pid_repro_rec["involucre_height_mm"],
                "involucre_width_mm": _pid_repro_rec["involucre_width_mm"],
                "AR_inv": _pid_repro_rec["AR_inv"],
                "leaf_zone": _leaf_zone_str,
                "leaf_zone_confidence": round(float(_leaf_zone_conf), 4) if _leaf_zone_conf is not None and np.isfinite(_leaf_zone_conf) else np.nan,
                "dist_from_crown_norm": round(float(_dist_crown_norm), 4) if _dist_crown_norm is not None and np.isfinite(_dist_crown_norm) else np.nan,
                "cy_rel_plant_bbox": round(float(inst_cy_rel_plant_bbox.get(inst.leaf_id, np.nan)), 4),
                "rosette_bbox_fallback_used": bool(rosette_bbox_fallback_used),
            })

        failure_record = None
        if not extracted_records:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "no_leaves_detected",
                "details": "LeafMachine3 produced 0 non-empty leaf masks",
                "pixels_per_mm": round(float(pixels_per_mm), 4) if is_calibrated_mm else np.nan,
                "is_calibrated_mm": bool(is_calibrated_mm),
            }

        return _finish(
            extracted_records,
            failure_record,
            repro_rec=repro_summary,
            pheno_rec=_pheno_rec,
            crowns_out=_crowns_out,
            centroids_out=_centroids_out,
            petiole_out=_petiole_out,
            contours_out=_contours_by_leaf,
            repro_pid_out=_repro_by_pid_out,
            pheno_pid_out=_pheno_by_pid_out,
        )

    def _flush_manifests(
        self,
        all_extracted_records: List[Dict[str, Any]],
        all_failed_records: List[Dict[str, Any]],
        manifest_path: Path,
        failed_qc_path: Path,
        repro_manifest_path: Path,
        pheno_manifest_path: Path,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Builds and atomically persists all 4 canonical Phase 2 output tables."""
        manifest_cols = [
            "catalogNumber", "plant_individual_id", "leaf_id", "solidity",
            "circularity", "extent", "aspect_ratio", "area_px", "area_mm2",
            "archetype_score", "bilateral_dice", "si_a", "lm3_gates_pass",
            "n_lamina_components", "scale_source", "pixels_per_mm", "is_calibrated_mm",
            "mask_path", "mask_is_upright", "contour_path", "ect_path", "leaf_zone",
        ]
        extracted_df = pd.DataFrame(all_extracted_records)
        if extracted_df.empty:
            extracted_df = pd.DataFrame(columns=manifest_cols + [
                "midrib_angle_deg", "petiole_length_px", "blade_length_px",
                "petiole_length_mm", "blade_length_mm",
                "petiole_to_blade_ratio", "capitula_count", "involucre_height_px",
                "involucre_width_px", "involucre_height_mm", "involucre_width_mm",
                "AR_inv",
            ])
        else:
            if "scale_source" not in extracted_df.columns:
                extracted_df["scale_source"] = extracted_df.get("pixels_per_mm", pd.Series(dtype=float)).apply(
                    lambda x: "ruler_calibrated" if (pd.notna(x) and float(x) > 0) else "uncalibrated"
                )
            if "is_calibrated_mm" not in extracted_df.columns:
                extracted_df["is_calibrated_mm"] = [
                    bool(pd.notna(px) and float(px) > 0 and str(src).strip().lower() not in {"mp_estimated", "predicted_from_megapixels", "uncalibrated"})
                    for px, src in zip(
                        extracted_df.get("pixels_per_mm", pd.Series([np.nan] * len(extracted_df))),
                        extracted_df["scale_source"],
                    )
                ]
            if "plant_individual_id" not in extracted_df.columns:
                extracted_df["plant_individual_id"] = 0
            for pet_col in ("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio"):
                if pet_col not in extracted_df.columns:
                    extracted_df[pet_col] = np.nan
                else:
                    extracted_df[pet_col] = pd.to_numeric(extracted_df[pet_col], errors="coerce")
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

        repro_cols = [
            "catalogNumber", "plant_individual_id", "capitula_count", "involucre_height_px",
            "involucre_width_px", "involucre_height_mm", "involucre_width_mm",
            "AR_inv", "plant_height_cm",
            "scape_to_basal_leaf_ratio", "is_bent_stem"
        ]
        repro_rows_all: List[Dict[str, Any]] = []
        for c_k, base_rep in self.voucher_repro_records.items():
            pid_map_r = self.voucher_repro_by_pid.get(c_k)
            if pid_map_r and len(pid_map_r) > 1:
                for p_k in sorted(pid_map_r.keys()):
                    repro_rows_all.append(pid_map_r[p_k])
            elif pid_map_r and len(pid_map_r) == 1:
                repro_rows_all.append(next(iter(pid_map_r.values())))
            else:
                repro_rows_all.append(base_rep)

        repro_df = pd.DataFrame(repro_rows_all)
        if repro_df.empty:
            repro_df = pd.DataFrame(columns=repro_cols)
        else:
            if "plant_individual_id" not in repro_df.columns:
                repro_df["plant_individual_id"] = 0
            if "plant_height_cm" not in repro_df.columns:
                repro_df["plant_height_cm"] = np.nan
            if "scape_to_basal_leaf_ratio" not in repro_df.columns:
                repro_df["scape_to_basal_leaf_ratio"] = np.nan
            if "is_bent_stem" not in repro_df.columns:
                repro_df["is_bent_stem"] = False
            repro_df = repro_df[[c for c in repro_cols if c in repro_df.columns]]

        pheno_cols = [
            "catalogNumber", "plant_individual_id", "total_capitula", "n_anthesis", "n_fruiting",
            "n_bud", "voucher_phenological_state", "dominant_pappus_ratio"
        ]
        pheno_rows_all: List[Dict[str, Any]] = []
        for c_k, base_phe in self.voucher_pheno_records.items():
            pid_map_p = self.voucher_pheno_by_pid.get(c_k)
            if pid_map_p and len(pid_map_p) > 1:
                for p_k in sorted(pid_map_p.keys()):
                    pheno_rows_all.append(pid_map_p[p_k])
            elif pid_map_p and len(pid_map_p) == 1:
                pheno_rows_all.append(next(iter(pid_map_p.values())))
            else:
                pheno_rows_all.append(base_phe)

        pheno_df = pd.DataFrame(pheno_rows_all)
        if pheno_df.empty:
            pheno_df = pd.DataFrame(columns=pheno_cols)
        else:
            if "plant_individual_id" not in pheno_df.columns:
                pheno_df["plant_individual_id"] = 0
            pheno_df = pheno_df[[c for c in pheno_cols if c in pheno_df.columns]]

        self._atomic_save_csv(extracted_df, manifest_path)
        self._atomic_save_csv(failed_df, failed_qc_path)
        self._atomic_save_csv(repro_df, repro_manifest_path)
        self._atomic_save_csv(pheno_df, pheno_manifest_path)
        return extracted_df, failed_df, repro_df, pheno_df

    def run(self, limit: Optional[int] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Executes full extraction pipeline across curated vouchers with resumption support."""
        logger.info("==================================================================")
        logger.info("Starting Step 02: Direct Segmentation & Botanical Leaf Extraction")
        logger.info(f"Vouchers Metadata: {self.vouchers_csv}")
        logger.info(f"LeafMachine3 Model: {self.model_weights}")
        logger.info(f"Output Directory: {self.output_dir}")
        logger.info(f"Device: {self.device}")
        logger.info(f"Force Overwrite: {self.force}")
        logger.info(
            "Segmentation Mode: Read-Only LeafMachine3 Ingestion Adapter "
            "(100% LM3 masks + native QC & organismal indexing saved to SQLite; EFA filtering in Step 03)"
        )
        logger.info("==================================================================")

        if not self.vouchers_csv.exists():
            raise FileNotFoundError(f"Curated vouchers table not found at {self.vouchers_csv}")

        df_vouchers = pd.read_csv(self.vouchers_csv)
        logger.info(f"Loaded {len(df_vouchers)} voucher records.")

        if "catalogNumber" not in df_vouchers.columns or "image_path" not in df_vouchers.columns:
            raise KeyError("Vouchers CSV must contain 'catalogNumber' and 'image_path' columns.")

        # Automatic one-time migration of existing uncompressed .jpg/.jpeg in raw_vouchers/ to .webp
        try:
            from PIL import Image as _PILImage

            _webp_q = int(getattr(getattr(self.config, "harvesting", None), "webp_quality", 92) or 92)
            _raw_vouchers_dir = self.output_dir / "raw_vouchers"
            if _raw_vouchers_dir.is_dir():
                _migrated_map: Dict[str, str] = {}
                for _jpg_f in list(_raw_vouchers_dir.glob("*.jpg")) + list(_raw_vouchers_dir.glob("*.jpeg")):
                    _webp_f = _jpg_f.with_suffix(".webp")
                    if not _webp_f.exists() or _webp_f.stat().st_size == 0:
                        try:
                            with _PILImage.open(_jpg_f) as _im:
                                _im.convert("RGB").save(_webp_f, format="WEBP", quality=_webp_q, method=4)
                            if _webp_f.exists() and _webp_f.stat().st_size > 0:
                                _jpg_f.unlink()
                                _migrated_map[_jpg_f.name] = _webp_f.name
                        except Exception as _mig_exc:
                            logger.debug(f"Skipped JPEG->WebP migration for {_jpg_f}: {_mig_exc}")
                    elif _webp_f.stat().st_size > 0:
                        try:
                            _jpg_f.unlink()
                            _migrated_map[_jpg_f.name] = _webp_f.name
                        except Exception as _unlink_exc:
                            logger.debug(f"Could not remove migrated JPEG {_jpg_f}: {_unlink_exc}")
                if _migrated_map and "image_path" in df_vouchers.columns:
                    df_vouchers["image_path"] = df_vouchers["image_path"].apply(
                        lambda p: str(Path(str(p)).with_name(_migrated_map[Path(str(p)).name]))
                        if pd.notna(p) and Path(str(p)).name in _migrated_map
                        else p
                    )
        except Exception as _outer_mig_exc:
            logger.debug(f"Raw voucher JPEG->WebP migration check skipped: {_outer_mig_exc}")

        if limit and limit > 0:
            df_vouchers = df_vouchers.head(limit)
            logger.info(f"Applied execution limit: processing first {limit} vouchers.")

        manifest_path = self.tables_dir / "extracted_leaf_manifest.csv"
        failed_qc_path = self.tables_dir / "failed_qc_vouchers.csv"
        repro_manifest_path = self.tables_dir / "voucher_reproductive_metrics.csv"
        pheno_manifest_path = self.tables_dir / "voucher_phenological_states.csv"

        existing_extracted_by_cat: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
        existing_failed_by_cat: Dict[str, Dict[str, Any]] = {}
        existing_repro_by_cat: Dict[str, Dict[str, Any]] = {}
        existing_pheno_by_cat: Dict[str, Dict[str, Any]] = {}

        use_db = (self.db_path is not None) or is_canonical_or_explicit_db_context(self.output_dir)
        db_conn = None
        if use_db:
            try:
                db_conn = get_db_connection(db_path=self.db_path, init_schema=True)
                if not df_vouchers.empty:
                    upsert_vouchers_records(db_conn, df_vouchers.to_dict(orient="records"))
            except Exception as db_init_exc:
                logger.warning(f"SQLite initialization skipped ({db_init_exc}); falling back to CSV manifests.")
                db_conn = None

        completed_in_db: Set[str] = set()
        if not self.force:
            if db_conn is not None:
                try:
                    target_cats = df_vouchers["catalogNumber"].astype(str).str.strip().tolist()
                    completed_in_db = get_completed_qc_vouchers(db_conn, catalog_filter=target_cats)
                    if completed_in_db:
                        db_ext, db_fail, db_rep, db_phe = load_step02_records_from_db(db_conn, list(completed_in_db))
                        for k_cat, recs in db_ext.items():
                            existing_extracted_by_cat[k_cat].extend(recs)
                        existing_failed_by_cat.update(db_fail)
                        existing_repro_by_cat.update(db_rep)
                        existing_pheno_by_cat.update(db_phe)
                except Exception as db_load_exc:
                    logger.warning(f"Could not query completed vouchers from SQLite: {db_load_exc}")

            if manifest_path.exists() and manifest_path.stat().st_size > 0:
                try:
                    prev_manifest = pd.read_csv(manifest_path)
                    if not prev_manifest.empty and "catalogNumber" in prev_manifest.columns:
                        for rec in prev_manifest.to_dict(orient="records"):
                            c_key = str(rec["catalogNumber"]).strip()
                            if c_key not in existing_extracted_by_cat:
                                existing_extracted_by_cat[c_key].append(rec)
                except Exception as e:
                    logger.warning(f"Could not load previous manifest at {manifest_path}: {e}")

            if failed_qc_path.exists() and failed_qc_path.stat().st_size > 0:
                try:
                    prev_failed = pd.read_csv(failed_qc_path)
                    if not prev_failed.empty and "catalogNumber" in prev_failed.columns:
                        for rec in prev_failed.to_dict(orient="records"):
                            existing_failed_by_cat.setdefault(str(rec["catalogNumber"]).strip(), rec)
                except Exception as e:
                    logger.warning(f"Could not load previous failed QC table at {failed_qc_path}: {e}")

            if repro_manifest_path.exists() and repro_manifest_path.stat().st_size > 0:
                try:
                    prev_repro = pd.read_csv(repro_manifest_path)
                    if not prev_repro.empty and "catalogNumber" in prev_repro.columns:
                        for rec in prev_repro.to_dict(orient="records"):
                            existing_repro_by_cat.setdefault(str(rec["catalogNumber"]).strip(), rec)
                except Exception as e:
                    logger.warning(f"Could not load previous reproductive table at {repro_manifest_path}: {e}")

            if pheno_manifest_path.exists() and pheno_manifest_path.stat().st_size > 0:
                try:
                    prev_pheno = pd.read_csv(pheno_manifest_path)
                    if not prev_pheno.empty and "catalogNumber" in prev_pheno.columns:
                        for rec in prev_pheno.to_dict(orient="records"):
                            existing_pheno_by_cat.setdefault(str(rec["catalogNumber"]).strip(), rec)
                except Exception as e:
                    logger.warning(f"Could not load previous phenology table at {pheno_manifest_path}: {e}")

        all_extracted_records: List[Dict[str, Any]] = []
        all_failed_records: List[Dict[str, Any]] = []

        total_vouchers = len(df_vouchers)
        processed_count = 0
        skipped_count = 0

        if not self.force:
            disk_scan = self._index_completed_vouchers() if len(completed_in_db) < total_vouchers else set()
            completed_on_disk = completed_in_db | disk_scan
        else:
            completed_on_disk = set()
        custom_process_voucher = "process_voucher" in self.__dict__
        custom_is_completed = (
            "is_voucher_completed" in self.__dict__
            or getattr(self.is_voucher_completed, "__func__", None) is not SegmentAndExtractPipeline.is_voucher_completed
        )
        custom_hooks = (
            custom_process_voucher
            or "predict" in getattr(self.stage1_detector, "__dict__", {})
            or "predict_batch" in getattr(self.stage1_detector, "__dict__", {})
            or "predict" in getattr(self.model_engine, "__dict__", {})
            or "predict_batch" in getattr(self.model_engine, "__dict__", {})
        )
        use_wave_engine = (not custom_hooks) and (total_vouchers > 4)

        def _prefetch_voucher_data(cat: str, img_p: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
            if not img_p.exists():
                for _ext in (".webp", ".jpg", ".jpeg", ".png"):
                    _cand_p = img_p.with_suffix(_ext)
                    if _cand_p.exists():
                        img_p = _cand_p
                        break
            if not img_p.exists():
                return None, None
            img = cv2.imread(str(img_p))
            if img is None:
                return None, None
            s1_lb: Optional[np.ndarray] = None
            if (
                hasattr(self.stage1_detector, "prepare_letterbox")
                and getattr(self.stage1_detector, "model", None) is not None
                and (self.force or cat not in getattr(self.stage1_detector, "cached_detections", {}))
            ):
                try:
                    s1_lb = self.stage1_detector.prepare_letterbox(img)
                except Exception:
                    s1_lb = None
            return img, s1_lb

        def _should_skip_voucher(cat_num: str) -> bool:
            if self.force:
                return False
            has_disk_artifacts = (
                (cat_num in completed_on_disk)
                or (self.is_voucher_completed(cat_num) if custom_is_completed else False)
            )
            has_manifest_entry = (cat_num in existing_extracted_by_cat) or (cat_num in existing_failed_by_cat)
            is_prev_failed_only = (
                (cat_num in existing_failed_by_cat)
                and (cat_num not in existing_extracted_by_cat)
                and (total_vouchers > 20)
            )
            return bool(
                (
                    has_disk_artifacts
                    and (
                        has_manifest_entry
                        or custom_process_voucher
                        or (not manifest_path.exists() and total_vouchers <= 20)
                    )
                )
                or is_prev_failed_only
            )

        def _retire_skipped_voucher(cat_num: str) -> None:
            nonlocal skipped_count
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

        voucher_rows = list(df_vouchers.iterrows())
        prev_cv2_threads: Optional[int] = None
        prefetch_pool: Optional[concurrent.futures.ThreadPoolExecutor] = None
        prep_pool: Optional[concurrent.futures.ThreadPoolExecutor] = None
        post_pool: Optional[concurrent.futures.ThreadPoolExecutor] = None
        retire_pool: Optional[concurrent.futures.ThreadPoolExecutor] = None

        try:
            if use_wave_engine:
                if self.threads is not None:
                    t_val = max(1, int(self.threads))
                    if t_val == 1:
                        prefetch_workers = 1
                        prep_workers = 1
                        post_workers = 1
                    else:
                        prefetch_workers = max(4, min(32, t_val // 2))
                        prep_workers = max(4, min(48, (t_val * 3) // 4))
                        post_workers = max(4, min(64, t_val))
                else:
                    n_cpu = os.cpu_count() or 32
                    # On 32-core / 64-thread SMT CPUs, scale worker pools to full physical core count (32) by default
                    eff_cores = min(32, max(4, n_cpu // 2)) if n_cpu >= 32 else max(4, (n_cpu * 3) // 4)
                    prefetch_workers = max(4, min(24, eff_cores // 2))
                    prep_workers = max(4, min(32, (eff_cores * 3) // 4))
                    post_workers = max(4, min(64, eff_cores))
                if self.prefetch_workers_override is not None:
                    prefetch_workers = max(1, int(self.prefetch_workers_override))
                if self.prep_workers_override is not None:
                    prep_workers = max(1, int(self.prep_workers_override))
                if self.post_workers_override is not None:
                    post_workers = max(1, int(self.post_workers_override))
                wave_size = max(1, int(self.stage1_wave_size))
                prefetch_window = max(wave_size * 3, prefetch_workers * 2)

                try:
                    prev_cv2_threads = cv2.getNumThreads()
                    cv2.setNumThreads(1 if (self.threads is not None and int(self.threads) == 1) else 2)
                except Exception:
                    prev_cv2_threads = None

                # Build deterministic ordered retirement plan separating completed skips from active GPU/CPU work
                ordered_plan: List[Tuple[str, Any]] = []
                active_tasks: List[Tuple[int, str, Path, str]] = []
                seen_active_cats: Set[str] = set()
                for _, row in voucher_rows:
                    cat_num = str(row["catalogNumber"]).strip()
                    raw_img_path = Path(str(row["image_path"]).strip())
                    if not raw_img_path.is_absolute():
                        raw_img_path = PROJECT_ROOT / raw_img_path
                    if not raw_img_path.exists():
                        for _ext in (".webp", ".jpg", ".jpeg", ".png"):
                            _cand_p = raw_img_path.with_suffix(_ext)
                            if _cand_p.exists():
                                raw_img_path = _cand_p
                                break
                    if _should_skip_voucher(cat_num) or (not self.force and cat_num in seen_active_cats):
                        ordered_plan.append(("skip", cat_num))
                    else:
                        a_idx = len(active_tasks)
                        taxon = str(row.get("scientificName", row.get("species_raw", ""))).strip()
                        active_tasks.append((a_idx, cat_num, raw_img_path, taxon))
                        ordered_plan.append(("active", a_idx))
                        seen_active_cats.add(cat_num)

                plan_cursor = 0
                completed_bundles: Dict[int, _VoucherFinalizeBundle] = {}

                def _drain_ordered_plan() -> None:
                    nonlocal plan_cursor, processed_count
                    while plan_cursor < len(ordered_plan):
                        kind, payload = ordered_plan[plan_cursor]
                        if kind == "skip":
                            _retire_skipped_voucher(payload)
                            plan_cursor += 1
                            if (processed_count + skipped_count) == total_vouchers and processed_count > 0:
                                logger.info(
                                    f"Progress: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}] "
                                    f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC)."
                                )
                        else:
                            if payload not in completed_bundles:
                                break
                            bundle = completed_bundles.pop(payload)
                            plan_cursor += 1
                            processed_count += 1
                            self._apply_finalize_bundle(bundle)
                            if bundle.extracted_records:
                                completed_on_disk.add(bundle.catalog_number)
                            all_extracted_records.extend(bundle.extracted_records)
                            if bundle.failure_record:
                                all_failed_records.append(bundle.failure_record)
                            if processed_count % 25 == 0 or (processed_count + skipped_count) == total_vouchers:
                                logger.info(
                                    f"Progress: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}] "
                                    f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC)."
                                )

                # Flush any leading skipped vouchers before spinning up thread pools
                _drain_ordered_plan()

                if active_tasks:
                    logger.info(
                        f"3-Stage Wave Engine active: {len(active_tasks)} vouchers queued "
                        f"(wave_size={wave_size}, s2_batch={self.stage2_batch_size}, "
                        f"prefetch_workers={prefetch_workers}, prep_workers={prep_workers}, post_workers={post_workers})."
                    )
                    prefetch_pool = concurrent.futures.ThreadPoolExecutor(max_workers=prefetch_workers)
                    prep_pool = concurrent.futures.ThreadPoolExecutor(max_workers=prep_workers)
                    post_pool = concurrent.futures.ThreadPoolExecutor(max_workers=post_workers)
                    retire_pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
                    self.stage1_detector._prep_executor = prep_pool  # type: ignore[attr-defined]
                    self.model_engine._prep_executor = prep_pool  # type: ignore[attr-defined]
                    if getattr(self, "pose_engine", None) is not None:
                        self.pose_engine._prep_executor = prep_pool  # type: ignore[attr-defined]
                    prefetch_futures: Dict[int, concurrent.futures.Future] = {}
                    retire_futures: List[concurrent.futures.Future] = []
                    use_s1_cache = bool(self.lm3_project_dir is not None and not self.force)

                    def _retire_wave_batch(
                        oldest_wave: List[Tuple[int, str, Path, concurrent.futures.Future]]
                    ) -> None:
                        wave_bundles: List[_VoucherFinalizeBundle] = []
                        for a_idx_w, cat_w, path_w, fut_w in oldest_wave:
                            try:
                                bundle = fut_w.result()
                            except Exception as exc_w:
                                logger.warning(f"Unexpected error processing voucher {cat_w} ({path_w}): {exc_w}")
                                bundle = _VoucherFinalizeBundle(
                                    catalog_number=cat_w,
                                    image_path=str(path_w),
                                    extracted_records=[],
                                    failure_record={
                                        "catalogNumber": cat_w,
                                        "image_path": str(path_w),
                                        "failure_reason": "SEGMENTATION_RUNTIME_ERROR",
                                        "details": f"Runtime exception ({type(exc_w).__name__}) during voucher extraction: {exc_w}",
                                    },
                                )
                            completed_bundles[a_idx_w] = bundle
                            wave_bundles.append(bundle)

                        _drain_ordered_plan()

                        if db_conn is not None and wave_bundles:
                            try:
                                with db_conn:
                                    for idx_b, b_item in enumerate(wave_bundles):
                                        sp_name = f"sp_v_{idx_b}"
                                        try:
                                            db_conn.execute(f"SAVEPOINT {sp_name};")
                                            record_voucher_segmentation_transaction(
                                                conn=db_conn,
                                                catalog_number=b_item.catalog_number,
                                                extracted_records=b_item.extracted_records,
                                                failure_record=b_item.failure_record,
                                                repro_record=b_item.repro_record,
                                                pheno_record=b_item.pheno_record,
                                                repro_by_pid=b_item.repro_by_pid,
                                                pheno_by_pid=b_item.pheno_by_pid,
                                                crowns_by_pid=b_item.crowns_by_pid,
                                                rosette_centroids_by_pid=b_item.rosette_centroids_by_pid,
                                                by_plant_petiole=b_item.by_plant_petiole,
                                                contour_coords_by_leaf=b_item.contour_coords_by_leaf,
                                                output_dir=self.output_dir,
                                                image_path=b_item.image_path,
                                                commit=False,
                                            )
                                            db_conn.execute(f"RELEASE SAVEPOINT {sp_name};")
                                        except Exception as v_db_exc:
                                            try:
                                                db_conn.execute(f"ROLLBACK TO SAVEPOINT {sp_name};")
                                                db_conn.execute(f"RELEASE SAVEPOINT {sp_name};")
                                            except Exception:
                                                pass
                                            logger.warning(
                                                f"Per-voucher SQLite transaction failed for {b_item.catalog_number}: {v_db_exc}"
                                            )
                            except Exception as wave_db_exc:
                                logger.warning(f"Wave SQLite WAL commit failed: {wave_db_exc}")

                    n_active = len(active_tasks)
                    for w_start in range(0, n_active, wave_size):
                        w_end = min(n_active, w_start + wave_size)
                        # Stage A: Keep bounded lookahead window of high-DPI sheets prefetched in RAM
                        la_limit = min(n_active, w_end + prefetch_window)
                        for la_idx in range(w_start, la_limit):
                            if la_idx not in prefetch_futures:
                                _, la_cat, la_path, _ = active_tasks[la_idx]
                                prefetch_futures[la_idx] = prefetch_pool.submit(_prefetch_voucher_data, la_cat, la_path)

                        wave_slice = active_tasks[w_start:w_end]
                        wave_prefetched: List[Tuple[int, str, Path, str, Optional[np.ndarray], Optional[np.ndarray]]] = []
                        valid_indices: List[int] = []
                        valid_s1_inputs: List[Tuple[np.ndarray, Optional[str], Optional[np.ndarray]]] = []
                        for local_i, (a_idx, cat_num, raw_img_path, taxon) in enumerate(wave_slice):
                            pre_img, pre_s1_lb = prefetch_futures.pop(a_idx).result()
                            wave_prefetched.append((a_idx, cat_num, raw_img_path, taxon, pre_img, pre_s1_lb))
                            if pre_img is not None:
                                valid_indices.append(local_i)
                                valid_s1_inputs.append((pre_img, cat_num, pre_s1_lb))

                        # Stage B1: Batched Stage 1 LeafMachine3 YOLO26x ONNX organ detection on GPU thread
                        s1_dets_by_local: Dict[int, List[Any]] = {}
                        if valid_s1_inputs:
                            s1_batch_out = self.stage1_detector.predict_batch(
                                valid_s1_inputs,
                                use_cache=use_s1_cache,
                                max_batch_size=wave_size,
                            )
                            for local_i, dets in zip(valid_indices, s1_batch_out):
                                s1_dets_by_local[local_i] = dets

                        # Stage B2: Parallel CPU crop extraction + LM3 resize transform across Threadripper pool
                        def _prep_s2_voucher(local_i: int) -> Tuple[int, Optional[Tuple[Any, ...]], List[Optional[Tuple[np.ndarray, int, int]]], bool]:
                            _, cat_v, _, _, img_v, _ = wave_prefetched[local_i]
                            dets_v = s1_dets_by_local.get(local_i, [])
                            assert img_v is not None
                            sh_v, sw_v = img_v.shape[:2]
                            sf_v = max(sh_v, sw_v) / 5600.0 if max(sh_v, sw_v) > 1000 else 0.0
                            hd_v = max(sh_v, sw_v) > 1000 and sf_v > 0.4
                            if not dets_v:
                                needs_ss = (not hd_v) or (self.stage1_detector.model is None)
                                ss_pre = (
                                    [self.model_engine.prepare_transformed_crop(img_v)]
                                    if (needs_ss and self.model_engine.model is not None)
                                    else []
                                )
                                return local_i, None, ss_pre, needs_ss

                            lm3_cal_v = self.ruler_calibrations.get(cat_v)
                            ppm_v: Optional[float] = normalize_pixels_per_mm(lm3_cal_v)
                            unpacked = self._unpack_stage1_candidates(
                                catalog_number=cat_v,
                                img_bgr=img_v,
                                stage1_detections=dets_v,
                                pixels_per_mm=ppm_v,
                            )
                            cand_crops = unpacked[0]
                            pre_crops = [
                                self.model_engine.prepare_transformed_crop(c[3])
                                for c in cand_crops
                            ]
                            return local_i, unpacked, pre_crops, False

                        unpacked_by_local: Dict[int, Tuple[Any, ...]] = {}
                        precrops_by_local: Dict[int, List[Optional[Tuple[np.ndarray, int, int]]]] = {}
                        ss_Fallback_locals: List[int] = []
                        if len(valid_indices) >= 2:
                            for loc_i, unp, pre_c, needs_ss in prep_pool.map(_prep_s2_voucher, valid_indices):
                                if unp is not None:
                                    unpacked_by_local[loc_i] = unp
                                    precrops_by_local[loc_i] = pre_c
                                elif needs_ss:
                                    ss_Fallback_locals.append(loc_i)
                                    precrops_by_local[loc_i] = pre_c
                        else:
                            for loc_i in valid_indices:
                                _, unp, pre_c, needs_ss = _prep_s2_voucher(loc_i)
                                if unp is not None:
                                    unpacked_by_local[loc_i] = unp
                                    precrops_by_local[loc_i] = pre_c
                                elif needs_ss:
                                    ss_Fallback_locals.append(loc_i)
                                    precrops_by_local[loc_i] = pre_c

                        # Stage B3: Cross-voucher batched Stage 2 LeafMachine3 YOLO26x-seg on GPU thread
                        flat_crops: List[np.ndarray] = []
                        flat_pre: List[Optional[Tuple[np.ndarray, int, int]]] = []
                        slice_ranges: List[Tuple[int, int, int]] = []
                        for loc_i in valid_indices:
                            if loc_i not in unpacked_by_local:
                                continue
                            cands = unpacked_by_local[loc_i][0]
                            if not cands:
                                slice_ranges.append((loc_i, 0, 0))
                                continue
                            s_idx = len(flat_crops)
                            for c_item in cands:
                                flat_crops.append(c_item[3])
                            flat_pre.extend(precrops_by_local.get(loc_i, [None] * len(cands)))
                            slice_ranges.append((loc_i, s_idx, len(flat_crops)))

                        flat_preds = (
                            self.model_engine.predict_batch(
                                flat_crops,
                                batch_size=self.stage2_batch_size,
                                pretransformed=flat_pre,
                            )
                            if flat_crops
                            else []
                        )

                        flat_pose_preds: List[Optional[Any]] = [None] * len(flat_crops)
                        if (
                            flat_crops
                            and flat_preds
                            and bool(getattr(self, "do_apply_landmarks", True))
                            and getattr(self, "pose_engine", None) is not None
                            and getattr(self.pose_engine, "model", None) is not None
                            and "predict_landmarks" not in getattr(self.pose_engine, "__dict__", {})
                            and hasattr(self.pose_engine, "predict_batch")
                        ):
                            pose_indices: List[int] = []
                            pose_crops: List[np.ndarray] = []
                            pose_lamina_masks: List[Optional[np.ndarray]] = []
                            pose_petiole_masks: List[Optional[np.ndarray]] = []
                            for idx_f, (c_bgr, c_preds) in enumerate(zip(flat_crops, flat_preds)):
                                if not c_preds:
                                    continue
                                has_leaf_f = any(
                                    int(p_cls) == 0 and np.any(pm_u8 > 0)
                                    for _bb, pm_u8, _p_sc, p_cls in c_preds
                                )
                                if has_leaf_f:
                                    _bm_f, _, _pm_f, _ = _build_native_lm3_crop_masks(
                                        c_preds,
                                        c_bgr.shape[:2],
                                        is_high_dpi=True,
                                        min_leaf_area=0,
                                        subtract_holes=bool(getattr(self, "subtract_holes", True)),
                                    )
                                    pose_indices.append(idx_f)
                                    pose_crops.append(c_bgr)
                                    pose_lamina_masks.append(_bm_f)
                                    pose_petiole_masks.append(_pm_f)
                            if pose_crops:
                                pose_batch_out = self.pose_engine.predict_batch(
                                    pose_crops,
                                    lamina_masks=pose_lamina_masks,
                                    petiole_masks=pose_petiole_masks,
                                    batch_size=max(8, self.stage2_batch_size),
                                )
                                for idx_f, p_res in zip(pose_indices, pose_batch_out):
                                    flat_pose_preds[idx_f] = p_res

                        s2_payload_by_local: Dict[int, Tuple[Any, ...]] = {}
                        for loc_i, s_idx, e_idx in slice_ranges:
                            unp = unpacked_by_local[loc_i]
                            v_preds = list(flat_preds[s_idx:e_idx]) if e_idx > s_idx else []
                            v_pose_preds = list(flat_pose_preds[s_idx:e_idx]) if e_idx > s_idx else []
                            s2_payload_by_local[loc_i] = (*unp, v_preds, v_pose_preds)

                        ss_payload_by_local: Dict[int, List[Any]] = {}
                        if ss_Fallback_locals:
                            ss_imgs = [wave_prefetched[loc_i][4] for loc_i in ss_Fallback_locals]  # type: ignore[misc]
                            ss_pre_list = [
                                (precrops_by_local[loc_i][0] if precrops_by_local.get(loc_i) else None)
                                for loc_i in ss_Fallback_locals
                            ]
                            ss_preds_batch = self.model_engine.predict_batch(
                                ss_imgs,  # type: ignore[arg-type]
                                batch_size=self.stage2_batch_size,
                                pretransformed=ss_pre_list,
                            )
                            for loc_i, ss_raw in zip(ss_Fallback_locals, ss_preds_batch):
                                ss_payload_by_local[loc_i] = ss_raw

                        # Stage C: Submit parallel CPU contour/mask/rosette post-processing while GPU moves to next wave
                        wave_post_tasks: List[Tuple[int, str, Path, concurrent.futures.Future]] = []
                        for local_i, (a_idx, cat_num, raw_img_path, taxon, pre_img, _pre_lb) in enumerate(wave_prefetched):
                            fut_post = post_pool.submit(
                                self.process_voucher,
                                cat_num,
                                raw_img_path,
                                taxon,
                                pre_img,
                                None,
                                s1_dets_by_local.get(local_i, []) if pre_img is not None else None,
                                s2_payload_by_local.get(local_i),
                                ss_payload_by_local.get(local_i),
                                True,
                            )
                            wave_post_tasks.append((a_idx, cat_num, raw_img_path, fut_post))

                        retire_futures.append(retire_pool.submit(_retire_wave_batch, wave_post_tasks))
                        while len(retire_futures) >= 4:
                            retire_futures.pop(0).result()

                    while retire_futures:
                        retire_futures.pop(0).result()
                _drain_ordered_plan()
            else:
                for _, row in voucher_rows:
                    cat_num = str(row["catalogNumber"]).strip()
                    raw_img_path = Path(str(row["image_path"]).strip())
                    if not raw_img_path.is_absolute():
                        raw_img_path = PROJECT_ROOT / raw_img_path
                    if not raw_img_path.exists():
                        for _ext in (".webp", ".jpg", ".jpeg", ".png"):
                            _cand_p = raw_img_path.with_suffix(_ext)
                            if _cand_p.exists():
                                raw_img_path = _cand_p
                                break

                    if _should_skip_voucher(cat_num):
                        _retire_skipped_voucher(cat_num)
                        continue

                    processed_count += 1
                    taxon = str(row.get("scientificName", row.get("species_raw", ""))).strip()
                    try:
                        extracted, failure = self.process_voucher(cat_num, raw_img_path, taxon=taxon)
                    except Exception as exc:
                        logger.warning(f"Unexpected error processing voucher {cat_num} ({raw_img_path}): {exc}")
                        extracted = []
                        failure = {
                            "catalogNumber": cat_num,
                            "image_path": str(raw_img_path),
                            "failure_reason": "SEGMENTATION_RUNTIME_ERROR",
                            "details": f"Runtime exception ({type(exc).__name__}) during voucher extraction: {exc}",
                        }

                    if extracted:
                        completed_on_disk.add(cat_num)
                    all_extracted_records.extend(extracted)
                    if failure:
                        all_failed_records.append(failure)

                    if db_conn is not None:
                        try:
                            record_voucher_segmentation_transaction(
                                conn=db_conn,
                                catalog_number=cat_num,
                                extracted_records=extracted,
                                failure_record=failure,
                                repro_record=self.voucher_repro_records.get(cat_num),
                                pheno_record=self.voucher_pheno_records.get(cat_num),
                                repro_by_pid=self.voucher_repro_by_pid.get(cat_num),
                                pheno_by_pid=self.voucher_pheno_by_pid.get(cat_num),
                                crowns_by_pid=self.voucher_crowns_by_pid.get(cat_num),
                                rosette_centroids_by_pid=self.voucher_rosette_centroids_by_pid.get(cat_num),
                                by_plant_petiole=self.voucher_petiole_by_pid.get(cat_num),
                                contour_coords_by_leaf=self.voucher_contour_coords_by_leaf.get(cat_num),
                                output_dir=self.output_dir,
                                image_path=str(raw_img_path),
                            )
                        except Exception as db_tx_exc:
                            logger.warning(f"Per-voucher SQLite transaction failed for {cat_num}: {db_tx_exc}")

                    if processed_count % 25 == 0 or (processed_count + skipped_count) == total_vouchers:
                        logger.info(
                            f"Progress: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}] "
                            f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC)."
                        )
        finally:
            self.stage1_detector._prep_executor = None  # type: ignore[attr-defined]
            self.model_engine._prep_executor = None  # type: ignore[attr-defined]
            if getattr(self, "pose_engine", None) is not None:
                self.pose_engine._prep_executor = None  # type: ignore[attr-defined]
            if retire_pool is not None:
                retire_pool.shutdown(wait=True)
            if post_pool is not None:
                post_pool.shutdown(wait=True)
            if prep_pool is not None:
                prep_pool.shutdown(wait=False, cancel_futures=True)
            if prefetch_pool is not None:
                prefetch_pool.shutdown(wait=False, cancel_futures=True)
            if prev_cv2_threads is not None:
                try:
                    cv2.setNumThreads(prev_cv2_threads)
                except Exception:
                    pass
            if db_conn is not None:
                try:
                    db_conn.close()
                except Exception:
                    pass

        logger.info(
            f"Extraction Summary: [Processed: {processed_count} | Skipped: {skipped_count} | Total: {total_vouchers}]"
        )

        extracted_df, failed_df, repro_df, pheno_df = self._flush_manifests(
            all_extracted_records,
            all_failed_records,
            manifest_path,
            failed_qc_path,
            repro_manifest_path,
            pheno_manifest_path,
        )

        logger.info("==================================================================")
        logger.info(f"Extraction completed!")
        logger.info(f"  Extracted Leaves: {len(extracted_df)} -> {manifest_path}")
        logger.info(f"  Failed QC Vouchers: {len(failed_df)} -> {failed_qc_path}")
        logger.info(f"  Reproductive Metrics: {len(repro_df)} -> {repro_manifest_path}")
        logger.info(f"  Phenological States: {len(pheno_df)} -> {pheno_manifest_path}")

        logger.info("==================================================================")

        if len(extracted_df) > 0:
            try:
                env_run_dir = os.environ.get("PACKERA_RUN_DIR", "").strip()
                if env_run_dir:
                    figures_dir = Path(env_run_dir).resolve() / "outputs" / "figures"
                elif self.output_dir.name == "data":
                    figures_dir = self.output_dir.parent / "outputs" / "figures"
                else:
                    figures_dir = self.output_dir / "outputs" / "figures"
                generate_step02_publication_figures(
                    extracted_df=extracted_df,
                    failed_df=failed_df,
                    repro_df=repro_df,
                    pheno_df=pheno_df,
                    figures_dir=figures_dir,
                )
            except Exception as fig_exc:
                logger.warning(f"Step 02 publication figure generation skipped due to non-fatal error: {fig_exc}")

        self.shutdown()
        return extracted_df, failed_df


Step02SegmentAndExtract = SegmentAndExtractPipeline

from scripts.vision.step02_figures import (
    _atomic_save_matplotlib_fig,
    _discover_downstream_table,
    _is_valid_smooth_closed_mask,
    _parse_bool_val,
    _resolve_asset_path,
    _standardize_packera_species_py,
    generate_step02_publication_figures,
    regenerate_step02_figures_from_tables,
)


# =============================================================================
# CLI Interface
# =============================================================================

def parse_args() -> argparse.Namespace:
    cfg = PipelineConfig.from_yaml()
    parser = argparse.ArgumentParser(
        description="Direct LeafMachine3 ONNX segmentation, upright leaf/ECT extraction, and 120-vertex contour export."
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
        help="Path to LeafMachine3 YOLO26x leaf segmenter ONNX model checkpoint."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=cfg.paths.workspace_root / "data",
        help="Root output directory for masks, contours, and manifest tables."
    )
    parser.add_argument(
        "--lm3-dir",
        dest="lm3_project_dir",
        type=Path,
        default=None,
        help="Optional LeafMachine3 project output directory for cached stage-1 detections and ruler calibrations.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=cfg.segmentation.device if torch_cuda_available() else "cpu",
        choices=["cuda", "cpu"],
        help="Compute device for inference (default: cuda if available, else cpu)."
    )
    parser.add_argument(
        "--score-thresh",
        type=float,
        default=cfg.segmentation.score_thresh,
        help=f"Detection confidence score threshold (default: {cfg.segmentation.score_thresh})."
    )
    parser.add_argument(
        "-j",
        "--threads",
        "--workers",
        dest="threads",
        type=int,
        default=None,
        help="Maximum CPU worker threads for image prefetching, tensor prep, and post-processing (default: auto-scaled for 32-core/64-thread CPUs).",
    )
    parser.add_argument(
        "--gpu-batch-size",
        dest="gpu_batch_size",
        type=int,
        default=None,
        help="GPU inference batch size for Stage 2 LeafMachine3 leaf crops and Stage 1 wave scaling (default: 10 on CUDA, 8 on CPU).",
    )
    parser.add_argument(
        "--batch-size",
        dest="gpu_batch_size",
        type=int,
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--wave-size",
        dest="wave_size",
        type=int,
        default=None,
        help="Number of herbarium sheets processed per cross-voucher GPU wave (default: 16 or --gpu-batch-size).",
    )
    parser.add_argument(
        "--prefetch-workers",
        dest="prefetch_workers",
        type=int,
        default=None,
        help="Optional override for Stage A disk image prefetch worker thread count.",
    )
    parser.add_argument(
        "--prep-workers",
        dest="prep_workers",
        type=int,
        default=None,
        help="Optional override for Stage B2/B3 CPU crop prep and prototype mask decode worker thread count.",
    )
    parser.add_argument(
        "--post-workers",
        dest="post_workers",
        type=int,
        default=None,
        help="Optional override for Stage C contour/ECT/rosette post-processing worker thread count.",
    )
    parser.add_argument(
        "--purge-segments",
        "--purge-old-data",
        dest="purge_segments",
        action="store_true",
        default=False,
        help="Purge all previous Step 02 to 07 segments, masks, contours, run directories, and SQLite segmentation tables before running.",
    )
    parser.add_argument(
        "--force",
        dest="force",
        action="store_true",
        default=False,
        help="Force re-segmentation and overwrite existing contour CSVs and masks.",
    )
    parser.add_argument(
        "--overwrite",
        dest="force",
        action="store_true",
        default=argparse.SUPPRESS,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of vouchers to process."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="Optional SQLite database path override (default: data/packera_pipeline.sqlite).",
    )
    parser.add_argument(
        "--regenerate-figures-only",
        dest="regenerate_figures_only",
        action="store_true",
        default=False,
        help="Reload Step 02 to 04 tables and regenerate post-QC Step 02 publication figures without running segmentation.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if getattr(args, "purge_segments", False):
        purge_segments_and_runs(
            workspace_root=PROJECT_ROOT,
            db_path=getattr(args, "db", None),
            output_dir=getattr(args, "output_dir", None),
            purge_runs=False,
            vacuum=True,
        )
    if getattr(args, "regenerate_figures_only", False):
        regenerate_step02_figures_from_tables(
            output_dir=args.output_dir,
            vouchers_csv=args.vouchers,
        )
        return
    pipeline = SegmentAndExtractPipeline(
        vouchers_csv=args.vouchers,
        model_weights=args.model_weights,
        output_dir=args.output_dir,
        device=args.device,
        score_thresh=args.score_thresh,
        force=bool(args.force or getattr(args, "purge_segments", False)),
        lm3_project_dir=getattr(args, "lm3_project_dir", None),
        db_path=getattr(args, "db", None),
        threads=getattr(args, "threads", None),
        gpu_batch_size=getattr(args, "gpu_batch_size", None),
        wave_size=getattr(args, "wave_size", None),
        prefetch_workers=getattr(args, "prefetch_workers", None),
        prep_workers=getattr(args, "prep_workers", None),
        post_workers=getattr(args, "post_workers", None),
    )
    pipeline.run(limit=args.limit)


if __name__ == "__main__":
    main()
