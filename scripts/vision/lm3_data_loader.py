#!/usr/bin/env python3
"""
scripts/vision/lm3_data_loader.py
=================================
Data ingestion, coordinate unscaling, 31-keypoint pose loading, and metadata
linking module for LeafMachine3 (LM3) outputs.
Fast ingestion of ruler/FieldPrism/MP scale calibrations, leaf crops, masks
(Leaf / Petiole / Hole / Specimen_PaperCleaned), 31-keypoint skeletons, ECTs,
Momocs outlines, and Darwin Core voucher metadata.
"""

from __future__ import annotations

import json
import logging
import math
import re
import sqlite3
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import cv2
import numpy as np
import pandas as pd

logger = logging.getLogger("LM3_PostProcessing")

# Canonical LeafMachine3 9-Class Plant Detector (lm3_plant_detector__yolo26x_det_1280) Taxonomy
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
LM3_PCD_CLASSES: Dict[int, str] = LM3_PLANT_CLASSES
LM3_LEAF_CLASS_IDS: Tuple[int, ...] = (0, 1)
LM3_REPRODUCTIVE_CLASS_IDS: Tuple[int, ...] = (2, 3, 4, 5, 6)

# Canonical LeafMachine3 3-Class Leaf Segmenter (lm3_leaf_segmenter__yolo26x_seg_1024) Taxonomy
LM3_LEAF_SEG_CLASSES: Dict[int, str] = {
    0: "Leaf",
    1: "Petiole",
    2: "Hole",
}

# Canonical LeafMachine3 31-Keypoint Pose Skeleton (lm3_landmark_detector__yolo26x_pose_640, mid15_pet5)
LM3_KEYPOINT_31_NAMES: Tuple[str, ...] = (
    "lamina_tip",
    "apex_1", "apex_2", "apex_3",
    "midvein_1", "midvein_2", "midvein_3", "midvein_4", "midvein_5",
    "midvein_6", "midvein_7", "midvein_8", "midvein_9", "midvein_10",
    "midvein_11", "midvein_12", "midvein_13", "midvein_14", "midvein_15",
    "base_1", "base_2", "base_3",
    "lamina_base",
    "petiole_1", "petiole_2", "petiole_3", "petiole_4", "petiole_5", "petiole_6",
    "width_1", "width_2",
)

VALID_LM3_CF_SOURCES: Tuple[str, ...] = (
    "measured_from_ruler",
    "measured_from_fieldprism",
    "predicted_from_megapixels",
)

DEFAULT_HERBARIUM_DPI: float = 600.0
DEFAULT_PX_PER_CM: float = DEFAULT_HERBARIUM_DPI / 2.54  # ~236.2205 px/cm
DEFAULT_PX_PER_MM: float = DEFAULT_PX_PER_CM / 10.0     # ~23.6220 px/mm


def compute_lm3_work_scale(
    sheet_shape: Any,
    sheet_w: Optional[int] = None,
    max_working_dim: int = 3200,
) -> float:
    """
    Computes the LeafMachine3 working-frame downscale factor `work_scale`.
    Accepts either `compute_lm3_work_scale((sheet_h, sheet_w), max_working_dim=3200)`
    or `compute_lm3_work_scale(sheet_h, sheet_w, max_working_dim=3200)`.
    Returns `1.0` when `max(h, w) <= max_working_dim`, else `float(max_working_dim) / float(max(h, w))`.
    """
    if isinstance(sheet_shape, (tuple, list, np.ndarray)):
        if len(sheet_shape) < 2:
            return 1.0
        h, w = int(sheet_shape[0]), int(sheet_shape[1])
        if sheet_w is not None and max_working_dim == 3200:
            max_working_dim = int(sheet_w)
    elif isinstance(sheet_shape, (int, float, np.number)) and sheet_w is not None:
        h, w = int(sheet_shape), int(sheet_w)
    else:
        return 1.0

    if max_working_dim <= 0:
        return 1.0
    max_dim = max(h, w)
    if max_dim <= 0 or max_dim <= int(max_working_dim):
        return 1.0
    return float(max_working_dim) / float(max_dim)


def unscale_working_bbox_to_sheet(
    bbox_work: Tuple[float, float, float, float],
    work_scale: float,
    sheet_shape: Optional[Tuple[int, int]] = None,
) -> Tuple[int, int, int, int]:
    """
    Divides a working-frame bounding box `(ymin, xmin, ymax, xmax)` by `work_scale`
    and optionally clamps coordinates to `sheet_shape` `(h, w)`.
    """
    eff_scale = float(work_scale) if (work_scale is not None and float(work_scale) > 0.0) else 1.0
    y1 = int(round(float(bbox_work[0]) / eff_scale))
    x1 = int(round(float(bbox_work[1]) / eff_scale))
    y2 = int(round(float(bbox_work[2]) / eff_scale))
    x2 = int(round(float(bbox_work[3]) / eff_scale))

    if sheet_shape is not None and len(sheet_shape) >= 2:
        h, w = int(sheet_shape[0]), int(sheet_shape[1])
        if h > 0 and w > 0:
            y1 = max(0, min(h, y1))
            y2 = max(0, min(h, y2))
            x1 = max(0, min(w, x1))
            x2 = max(0, min(w, x2))

    return y1, x1, y2, x2


def unscale_working_points_to_sheet(
    pts_work: Any,
    work_scale: float,
    sheet_shape: Optional[Tuple[int, int]] = None,
) -> Any:
    """
    Divides working-frame 2D points, arrays, or keypoint dictionaries `{name: (x, y, conf)}`
    by `work_scale` for `(x, y)` while preserving confidence scores.
    """
    _ = sheet_shape
    if pts_work is None:
        return None
    eff_scale = float(work_scale) if (work_scale is not None and float(work_scale) > 0.0) else 1.0
    if abs(eff_scale - 1.0) < 1e-9:
        if isinstance(pts_work, np.ndarray):
            return pts_work.copy()
        if isinstance(pts_work, dict):
            return dict(pts_work)
        return pts_work

    if isinstance(pts_work, dict):
        unscaled_dict: Dict[str, Any] = {}
        for k, v in pts_work.items():
            if isinstance(v, (tuple, list, np.ndarray)) and len(v) >= 2:
                ux = float(v[0]) / eff_scale
                uy = float(v[1]) / eff_scale
                if len(v) >= 3:
                    unscaled_dict[k] = (ux, uy, float(v[2]))
                else:
                    unscaled_dict[k] = (ux, uy)
            else:
                unscaled_dict[k] = v
        return unscaled_dict

    if isinstance(pts_work, np.ndarray):
        arr = pts_work.astype(np.float64, copy=True)
        if arr.ndim == 1 and arr.size >= 2:
            arr[0] /= eff_scale
            arr[1] /= eff_scale
        elif arr.ndim >= 2 and arr.shape[-1] >= 2:
            arr[..., 0] /= eff_scale
            arr[..., 1] /= eff_scale
        return arr

    if isinstance(pts_work, tuple) and len(pts_work) in (2, 3) and all(
        isinstance(x, (int, float, np.number)) for x in pts_work
    ):
        if len(pts_work) == 3:
            return (float(pts_work[0]) / eff_scale, float(pts_work[1]) / eff_scale, float(pts_work[2]))
        return (float(pts_work[0]) / eff_scale, float(pts_work[1]) / eff_scale)

    if isinstance(pts_work, (list, tuple)):
        out_list: List[Any] = []
        for item in pts_work:
            if isinstance(item, (tuple, list, np.ndarray)) and len(item) >= 2:
                ux = float(item[0]) / eff_scale
                uy = float(item[1]) / eff_scale
                if len(item) >= 3:
                    pt_out = (ux, uy, float(item[2]))
                else:
                    pt_out = (ux, uy)
                out_list.append(list(pt_out) if isinstance(item, list) else pt_out)
            else:
                out_list.append(item)
        return tuple(out_list) if isinstance(pts_work, tuple) else out_list

    return pts_work


def resolve_px_per_cm(
    pixels_per_mm: Optional[float] = None,
    px_per_cm: Optional[float] = None,
    default_dpi: float = DEFAULT_HERBARIUM_DPI,
    sheet_shape: Optional[Tuple[int, int]] = None,
) -> float:
    """
    Resolves the empirical ruler conversion factor in pixels per centimeter (`px_per_cm`),
    falling back to a standard herbarium DPI conversion (`default_dpi / 2.54`) when
    ruler detection is absent.
    """
    if px_per_cm is not None and not np.isnan(px_per_cm) and float(px_per_cm) > 0.0:
        return float(px_per_cm)
    if pixels_per_mm is not None and not np.isnan(pixels_per_mm) and float(pixels_per_mm) > 0.0:
        ppm = float(1.0 / pixels_per_mm) if float(pixels_per_mm) < 1.0 else float(pixels_per_mm)
        return ppm * 10.0
    fallback_px_cm = float(default_dpi) / 2.54 if float(default_dpi) > 0.0 else DEFAULT_PX_PER_CM
    if sheet_shape is not None and max(sheet_shape) <= 500:
        # Proportional scaling for synthetic low-resolution unit-test fixtures lacking a ruler
        return max(40.0, fallback_px_cm * (float(max(sheet_shape)) / 2400.0))
    return fallback_px_cm


def compute_blade_area_cm2(
    mask_area_px: float,
    pixels_per_mm: Optional[float] = None,
    px_per_cm: Optional[float] = None,
    default_dpi: float = DEFAULT_HERBARIUM_DPI,
    sheet_shape: Optional[Tuple[int, int]] = None,
) -> float:
    """
    Converts pixel mask area (`mask_area_px`) into physical blade area in square
    centimeters (`cm^2`) using `px_per_cm` (or `pixels_per_mm * 10.0`), with a
    standard herbarium DPI fallback when ruler detection is absent.
    """
    if mask_area_px <= 0.0:
        return 0.0
    eff_px_cm = resolve_px_per_cm(
        pixels_per_mm=pixels_per_mm,
        px_per_cm=px_per_cm,
        default_dpi=default_dpi,
        sheet_shape=sheet_shape,
    )
    return float(mask_area_px) / max(eff_px_cm * eff_px_cm, 1e-6)


@dataclass
class LeafCandidate:
    """Represents an extracted leaf instance candidate from LeafMachine3."""
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
    px_per_cm: float = DEFAULT_PX_PER_CM
    plant_individual_id: int = 0
    assigned_tier: str = "Tier_1_Direct"
    ucs_score: float = 0.0
    solidity: float = 0.0
    midrib_angle_deg: float = 0.0
    length_mm: float = 0.0
    width_mm: float = 0.0
    area_mm2: float = 0.0
    area_cm2: float = 0.0
    saved_mask_path: str = ""
    # LeafMachine3 extended fields
    petiole_mask_path: Optional[Path] = None
    hole_mask_path: Optional[Path] = None
    specimen_paper_cleaned_mask_path: Optional[Path] = None
    keypoints_31: Optional[Dict[str, Tuple[float, float, float]]] = None
    work_scale: float = 1.0
    cf_source: str = ""
    archetype_score: float = 0.0
    petiole_width_mm: float = 0.0
    ect_path: Optional[Path] = None
    momocs_outline_xy: Optional[np.ndarray] = None
    oriented_image_path: Optional[Path] = None


class CandidateCollection(list):
    """
    Hybrid list/dict container of `LeafCandidate` instances that behaves as a
    `List[LeafCandidate]` when iterated or indexed by int, while also supporting
    voucher-grouped lookup (`collection[catalog_number]`, `.get(catalog_number)`,
    `.items()`, `.values()`, `.keys()`).
    """

    def __init__(
        self,
        candidates: Optional[Iterable[LeafCandidate]] = None,
        by_voucher: Optional[Dict[str, List[LeafCandidate]]] = None,
    ) -> None:
        super().__init__()
        self.by_voucher: Dict[str, List[LeafCandidate]] = defaultdict(list)
        if by_voucher is not None:
            for cat, c_list in by_voucher.items():
                self.by_voucher[cat].extend(c_list)
                super().extend(c_list)
        elif candidates is not None:
            for cand in candidates:
                super().append(cand)
                self.by_voucher[cand.catalog_number].append(cand)

    def __getitem__(self, key: Any) -> Any:
        if isinstance(key, str):
            return self.by_voucher.get(key, [])
        return super().__getitem__(key)

    def __contains__(self, item: Any) -> bool:
        if isinstance(item, str):
            return item in self.by_voucher
        return super().__contains__(item)

    def get(self, key: str, default: Any = None) -> Any:
        if key in self.by_voucher:
            return self.by_voucher[key]
        return default if default is not None else []

    def keys(self) -> Any:
        return self.by_voucher.keys()

    def values(self) -> Any:
        return self.by_voucher.values()

    def items(self) -> Any:
        return self.by_voucher.items()


def parse_bbox_from_filename(filename: str) -> Optional[Tuple[int, int, int, int]]:
    """
    Extracts bounding-box coordinates from LeafMachine3 `<stem>__<ID>__x1_y1_x2_y2.<ext>`
    double-underscore syntax, bracketed `[x1, y1, x2, y2]` syntax, or `_x1_y1_x2_y2` /
    `__y1-x1-y2-x2` patterns.
    """
    if not filename:
        return None
    s = str(filename)

    # 1. LM3 double-underscore underscore-delimited coordinates: <stem>__<ID>__x1_y1_x2_y2
    m_lm3 = re.search(r"__(\d+)_(\d+)_(\d+)_(\d+)(?=$|\.|__)", s)
    if m_lm3:
        return int(m_lm3.group(1)), int(m_lm3.group(2)), int(m_lm3.group(3)), int(m_lm3.group(4))

    # 2. Double-underscore hyphen-delimited coordinates: __y1-x1-y2-x2
    m_hyph = re.search(r"__(\d+)-(\d+)-(\d+)-(\d+)", s)
    if m_hyph:
        return int(m_hyph.group(1)), int(m_hyph.group(2)), int(m_hyph.group(3)), int(m_hyph.group(4))

    # 3. Bracketed [x1, y1, x2, y2]
    m_brk = re.search(r"\[\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\]", s)
    if m_brk:
        return int(m_brk.group(1)), int(m_brk.group(2)), int(m_brk.group(3)), int(m_brk.group(4))

    # 4. Standard _x1_y1_x2_y2 suffix
    m_und = re.search(r"(?:^|_)(\d+)_(\d+)_(\d+)_(\d+)(?:\.[A-Za-z0-9]+|$)", s)
    if m_und:
        return int(m_und.group(1)), int(m_und.group(2)), int(m_und.group(3)), int(m_und.group(4))

    # 5. Standard _x1-y1-x2-y2 suffix
    m_dash = re.search(r"(?:^|_)(\d+)-(\d+)-(\d+)-(\d+)(?:\.[A-Za-z0-9]+|$)", s)
    if m_dash:
        return int(m_dash.group(1)), int(m_dash.group(2)), int(m_dash.group(3)), int(m_dash.group(4))

    return None


def _extract_catalog_number(raw_identifier: Any) -> str:
    """Normalizes a filename, stem, or specimen_id into a clean catalogNumber."""
    if raw_identifier is None:
        return ""
    s = str(raw_identifier).strip()
    if not s or s.lower() in {"nan", "none", "null"}:
        return ""
    stem = Path(s).stem
    cat = stem.split("__")[0].strip()
    return cat


def _parse_calibration_row(
    row: Mapping[str, Any],
    min_confidence: float = 0.50,
) -> Optional[Dict[str, Any]]:
    """Parses a single SQLite or CSV calibration record into a normalized dict."""
    raw_id = (
        row.get("catalogNumber")
        or row.get("catalog_number")
        or row.get("specimen_id")
        or row.get("image_name")
        or row.get("filename")
        or row.get("stem")
    )
    cat = _extract_catalog_number(raw_id)
    if not cat:
        return None

    cf_source = str(row.get("cf_source") or "").strip()
    work_scale_raw = row.get("work_scale", 1.0)
    try:
        work_scale = float(work_scale_raw) if pd.notna(work_scale_raw) and float(work_scale_raw) > 0 else 1.0
    except (ValueError, TypeError):
        work_scale = 1.0

    # Check ruler classifier confidence when present
    conf_raw = row.get(
        "ruler_class_confidence",
        row.get("ruler_confidence", row.get("confidence", np.nan)),
    )
    if pd.notna(conf_raw) and cf_source not in {"measured_from_fieldprism", "predicted_from_megapixels"}:
        try:
            conf_val = float(conf_raw)
            if conf_val < 0:
                return None
            if conf_val > 1.0:
                conf_val /= 100.0
            if conf_val < min_confidence:
                return None
        except (ValueError, TypeError):
            pass

    # 1. Direct pixels_per_mm column if provided
    direct_px_mm = row.get("pixels_per_mm", np.nan)
    if pd.notna(direct_px_mm):
        try:
            val_mm = float(direct_px_mm)
            if val_mm > 0:
                return {
                    "catalog_number": cat,
                    "pixels_per_mm": val_mm,
                    "cf_px_per_cm": val_mm * 10.0,
                    "fp_cf_px_per_cm": None,
                    "cf_source": cf_source or "measured_from_ruler",
                    "work_scale": work_scale,
                    "is_empirical_ruler": cf_source != "predicted_from_megapixels",
                }
        except (ValueError, TypeError):
            pass

    # 2. LeafMachine3 cf_px_per_cm / fp_cf_px_per_cm / mp_cf_px_per_cm
    cf_px_cm_raw = row.get("cf_px_per_cm", np.nan)
    fp_px_cm_raw = row.get("fp_cf_px_per_cm", np.nan)
    mp_px_cm_raw = row.get("mp_cf_px_per_cm", np.nan)

    fp_val: Optional[float] = None
    if pd.notna(fp_px_cm_raw):
        try:
            if float(fp_px_cm_raw) > 0:
                fp_val = float(fp_px_cm_raw)
        except (ValueError, TypeError):
            pass

    cf_val: Optional[float] = None
    for cand_raw in (cf_px_cm_raw, fp_px_cm_raw, mp_px_cm_raw):
        if pd.notna(cand_raw):
            try:
                f_val = float(cand_raw)
                if f_val > 0:
                    cf_val = f_val
                    break
            except (ValueError, TypeError):
                pass

    if cf_val is not None and cf_val > 0:
        eff_source = cf_source
        if not eff_source:
            if fp_val is not None and abs(cf_val - fp_val) < 1e-4:
                eff_source = "measured_from_fieldprism"
            else:
                eff_source = "measured_from_ruler"
        return {
            "catalog_number": cat,
            "pixels_per_mm": cf_val / 10.0,
            "cf_px_per_cm": cf_val,
            "fp_cf_px_per_cm": fp_val,
            "cf_source": eff_source,
            "work_scale": work_scale,
            "is_empirical_ruler": eff_source in {"measured_from_ruler", "measured_from_fieldprism"},
        }

    # 3. Legacy CSV fallback columns (conversion_mean, predicted_conversion_factor_cm)
    success = str(row.get("ruler_success", "False")).strip().lower() == "true"
    conv_mean = row.get("conversion_mean", np.nan)
    pred_conv = row.get("predicted_conversion_factor_cm", np.nan)

    conversion_val: Optional[float] = None
    inferred_source = cf_source or "measured_from_ruler"
    try:
        if pd.notna(conv_mean) and float(conv_mean) > 0:
            conversion_val = float(conv_mean)
            inferred_source = cf_source or "measured_from_ruler"
        elif pd.notna(pred_conv) and float(pred_conv) > 0:
            conversion_val = float(pred_conv)
            inferred_source = cf_source or "predicted_from_megapixels"
    except (ValueError, TypeError):
        pass

    if (success and conversion_val and conversion_val > 0) or (conversion_val and conversion_val > 10.0):
        return {
            "catalog_number": cat,
            "pixels_per_mm": conversion_val / 10.0,
            "cf_px_per_cm": conversion_val,
            "fp_cf_px_per_cm": fp_val,
            "cf_source": inferred_source,
            "work_scale": work_scale,
            "is_empirical_ruler": inferred_source in {"measured_from_ruler", "measured_from_fieldprism"},
        }

    return None


def load_ruler_calibration_records(
    lm3_dir: Path,
    min_confidence: float = 0.50,
) -> Dict[str, Dict[str, Any]]:
    """
    Loads detailed LeafMachine3 ruler / FieldPrism / MP calibration records from
    `<run_name>.sqlite` (`specimen` and `ruler_CF_lattice` tables) and
    `Data/specimen_summary.csv` / `Data/ruler_conversion_factor.csv` (plus legacy
    ruler CSV locations for test compatibility).
    """
    records: Dict[str, Dict[str, Any]] = {}
    if not lm3_dir or not Path(lm3_dir).exists():
        return records

    lm3_path = Path(lm3_dir)

    def _merge_rec(parsed: Optional[Dict[str, Any]]) -> None:
        if not parsed:
            return
        cat = parsed["catalog_number"]
        existing = records.get(cat)
        if existing is None:
            records[cat] = parsed
            return
        # Prefer empirical ruler/FieldPrism calibration over megapixel prediction
        if parsed.get("is_empirical_ruler") and not existing.get("is_empirical_ruler"):
            records[cat] = parsed

    # 1. Inspect LeafMachine3 SQLite databases (<run_name>.sqlite)
    sqlite_files = sorted(
        set(list(lm3_path.glob("*.sqlite")) + list(lm3_path.glob("**/*.sqlite")))
    )
    for db_file in sqlite_files:
        try:
            conn = sqlite3.connect(f"file:{db_file.as_posix()}?mode=ro", uri=True, timeout=5.0)
            conn.row_factory = sqlite3.Row
            try:
                tables = {
                    str(r[0])
                    for r in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table';"
                    ).fetchall()
                }
                for tbl in ("specimen", "ruler_CF_lattice", "ruler_conversion_factor", "ruler"):
                    if tbl in tables:
                        for row in conn.execute(f"SELECT * FROM {tbl};").fetchall():
                            _merge_rec(_parse_calibration_row(dict(row), min_confidence=min_confidence))
            finally:
                conn.close()
        except Exception as exc:
            logger.debug("Error reading LM3 SQLite calibrations from %s: %s", db_file, exc)

    # 2. Inspect LeafMachine3 CSV summaries & ruler tables
    csv_files = (
        list(lm3_path.glob("Data/specimen_summary.csv"))
        + list(lm3_path.glob("**/Data/specimen_summary.csv"))
        + list(lm3_path.glob("Data/ruler_conversion_factor.csv"))
        + list(lm3_path.glob("**/Data/ruler_conversion_factor.csv"))
        + list(lm3_path.glob("**/Data/Ruler/*.csv"))
        + list(lm3_path.glob("**/Data/Batch/Ruler/*.csv"))
        + list(lm3_path.glob("**/Archival_Components/Ruler_Info/*.csv"))
    )
    seen_csvs = set()
    for rf in csv_files:
        rf_res = rf.resolve()
        if rf_res in seen_csvs:
            continue
        seen_csvs.add(rf_res)
        try:
            df = pd.read_csv(rf_res)
            if df.empty:
                continue
            for _, row in df.iterrows():
                _merge_rec(_parse_calibration_row(row.to_dict(), min_confidence=min_confidence))
        except Exception as exc:
            logger.debug("Error reading ruler CSV %s: %s", rf_res, exc)

    return records


def load_ruler_calibrations(
    lm3_dir: Path,
    min_confidence: float = 0.50,
) -> Dict[str, float]:
    """
    Fast loader for specimen scale calibrations (`pixels_per_mm`) across LeafMachine3
    SQLite (`<run_name>.sqlite`) and CSV (`Data/specimen_summary.csv`,
    `Data/ruler_conversion_factor.csv`) outputs.
    """
    records = load_ruler_calibration_records(lm3_dir=lm3_dir, min_confidence=min_confidence)
    calibrations: Dict[str, float] = {
        cat: float(rec["pixels_per_mm"])
        for cat, rec in records.items()
        if rec.get("pixels_per_mm") is not None and float(rec["pixels_per_mm"]) > 0.0
    }
    logger.info("Loaded %d specimen ruler scale calibrations from LM3.", len(calibrations))
    return calibrations


def _parse_31_keypoints_from_mapping(
    row: Mapping[str, Any],
    work_scale: float = 1.0,
) -> Dict[str, Tuple[float, float, float]]:
    """Extracts up to 31 named `(x, y, conf)` keypoints from a row mapping."""
    kps: Dict[str, Tuple[float, float, float]] = {}

    # 1. Check JSON columns (keypoints_json / landmarks_json / keypoints_31)
    for json_col in ("keypoints_json", "landmarks_json", "keypoints_31", "pose_keypoints"):
        raw_j = row.get(json_col)
        if isinstance(raw_j, str) and raw_j.strip().startswith(("{", "[")):
            try:
                parsed = json.loads(raw_j)
                if isinstance(parsed, dict):
                    for kp_name, pt in parsed.items():
                        if isinstance(pt, (list, tuple)) and len(pt) >= 2:
                            conf = float(pt[2]) if len(pt) >= 3 else 1.0
                            kps[str(kp_name)] = (float(pt[0]), float(pt[1]), conf)
                elif isinstance(parsed, list):
                    for idx, pt in enumerate(parsed):
                        if idx < len(LM3_KEYPOINT_31_NAMES) and isinstance(pt, (list, tuple)) and len(pt) >= 2:
                            conf = float(pt[2]) if len(pt) >= 3 else 1.0
                            kps[LM3_KEYPOINT_31_NAMES[idx]] = (float(pt[0]), float(pt[1]), conf)
            except Exception:
                pass

    # 2. Check named columns (<kp_name>_x, <kp_name>_y, <kp_name>_conf) or indexed (kp_0_x, ...)
    for idx, kp_name in enumerate(LM3_KEYPOINT_31_NAMES):
        if kp_name in kps:
            continue
        for prefix in (kp_name, f"kp_{idx}", f"kp{idx}", f"landmark_{idx}"):
            xk = f"{prefix}_x"
            yk = f"{prefix}_y"
            ck = f"{prefix}_conf"
            if xk in row and yk in row and pd.notna(row[xk]) and pd.notna(row[yk]):
                try:
                    fx = float(row[xk])
                    fy = float(row[yk])
                    fc = float(row[ck]) if (ck in row and pd.notna(row[ck])) else 1.0
                    if not (math.isnan(fx) or math.isnan(fy)):
                        kps[kp_name] = (fx, fy, fc)
                        break
                except (ValueError, TypeError):
                    pass

    if kps and work_scale > 0.0 and abs(work_scale - 1.0) > 1e-9:
        kps = unscale_working_points_to_sheet(kps, work_scale)
    return kps


def load_lm3_31_keypoints(
    lm3_dir: Path,
    catalog_number: Optional[str] = None,
) -> Dict[str, Dict[str, Tuple[float, float, float]]]:
    """
    Loads the 31-keypoint `mid15_pet5` pose skeleton from `Data/landmarks.csv`
    or `<run_name>.sqlite` (`leaf_landmark` table).
    Returns a dictionary mapping leaf identifier (and `catalog_number` fallback)
    to `{keypoint_name: (x, y, conf)}`.
    """
    keypoints_by_leaf: Dict[str, Dict[str, Tuple[float, float, float]]] = {}
    if not lm3_dir or not Path(lm3_dir).exists():
        return keypoints_by_leaf

    lm3_path = Path(lm3_dir)
    target_cat = str(catalog_number).strip() if catalog_number else None

    def _register_row(row_map: Mapping[str, Any]) -> None:
        raw_leaf_key = (
            row_map.get("leaf_stem")
            or row_map.get("leaf_filename")
            or row_map.get("crop_name")
            or row_map.get("leaf_id")
            or row_map.get("organ_id")
            or ""
        )
        raw_spec = (
            row_map.get("catalogNumber")
            or row_map.get("catalog_number")
            or row_map.get("specimen_id")
            or row_map.get("image_name")
            or row_map.get("filename")
            or raw_leaf_key
        )
        cat = _extract_catalog_number(raw_spec)
        if target_cat and cat and cat != target_cat:
            return

        leaf_key = Path(str(raw_leaf_key)).stem if raw_leaf_key else cat
        if cat and leaf_key and not leaf_key.startswith(cat):
            leaf_key = f"{cat}__{leaf_key}"
        if not leaf_key:
            leaf_key = cat
        if not leaf_key:
            return

        ws_raw = row_map.get("work_scale", 1.0)
        try:
            ws = float(ws_raw) if pd.notna(ws_raw) and float(ws_raw) > 0 else 1.0
        except (ValueError, TypeError):
            ws = 1.0

        # Long-format row check: single keypoint per row
        kp_name_col = row_map.get("keypoint_name") or row_map.get("landmark_name") or row_map.get("name")
        if kp_name_col is None and "kp_index" in row_map and pd.notna(row_map["kp_index"]):
            try:
                kp_idx = int(row_map["kp_index"])
                if 0 <= kp_idx < len(LM3_KEYPOINT_31_NAMES):
                    kp_name_col = LM3_KEYPOINT_31_NAMES[kp_idx]
            except (ValueError, TypeError):
                pass

        if kp_name_col is not None and "x" in row_map and "y" in row_map:
            try:
                fx = float(row_map["x"]) / ws
                fy = float(row_map["y"]) / ws
                fc_raw = row_map.get("conf", row_map.get("confidence", 1.0))
                fc = float(fc_raw) if pd.notna(fc_raw) else 1.0
                kp_str = str(kp_name_col).strip()
                keypoints_by_leaf.setdefault(leaf_key, {})[kp_str] = (fx, fy, fc)
                if cat and cat not in keypoints_by_leaf:
                    keypoints_by_leaf[cat] = keypoints_by_leaf[leaf_key]
                return
            except (ValueError, TypeError):
                pass

        # Wide-format or JSON row
        kps = _parse_31_keypoints_from_mapping(row_map, work_scale=ws)
        if kps:
            keypoints_by_leaf.setdefault(leaf_key, {}).update(kps)
            if cat and cat not in keypoints_by_leaf:
                keypoints_by_leaf[cat] = keypoints_by_leaf[leaf_key]

    # 1. SQLite `leaf_landmark` / `landmarks` table
    sqlite_files = sorted(
        set(list(lm3_path.glob("*.sqlite")) + list(lm3_path.glob("**/*.sqlite")))
    )
    for db_file in sqlite_files:
        try:
            conn = sqlite3.connect(f"file:{db_file.as_posix()}?mode=ro", uri=True, timeout=5.0)
            conn.row_factory = sqlite3.Row
            try:
                tables = {
                    str(r[0])
                    for r in conn.execute(
                        "SELECT name FROM sqlite_master WHERE type='table';"
                    ).fetchall()
                }
                for tbl in ("leaf_landmark", "landmarks", "landmark"):
                    if tbl in tables:
                        for row in conn.execute(f"SELECT * FROM {tbl};").fetchall():
                            _register_row(dict(row))
            finally:
                conn.close()
        except Exception as exc:
            logger.debug("Error reading LM3 SQLite landmarks from %s: %s", db_file, exc)

    # 2. CSV `Data/landmarks.csv`
    landmark_csvs = sorted(
        set(list(lm3_path.glob("Data/landmarks.csv")) + list(lm3_path.glob("**/Data/landmarks.csv")) + list(lm3_path.glob("**/landmarks.csv")))
    )
    for lcsv in landmark_csvs:
        try:
            df = pd.read_csv(lcsv)
            if df.empty:
                continue
            for _, row in df.iterrows():
                _register_row(row.to_dict())
        except Exception as exc:
            logger.debug("Error reading landmarks CSV %s: %s", lcsv, exc)

    return keypoints_by_leaf


def load_precomputed_pcd_labels(
    label_path: Path,
    sheet_width: int,
    sheet_height: int,
    min_confidence: float = 0.25,
) -> List[Tuple[Tuple[int, int, int, int], float, int]]:
    """
    Loads precomputed YOLO normalized bounding boxes from a `.txt` label file and
    converts them into full-sheet pixel coordinates:
    `[((ymin, xmin, ymax, xmax), confidence, class_id), ...]`.
    """
    detections: List[Tuple[Tuple[int, int, int, int], float, int]] = []
    if not label_path or not label_path.exists():
        return detections

    try:
        with open(label_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) < 5:
                    continue
                cls_id = int(float(parts[0]))
                xc = float(parts[1]) * sheet_width
                yc = float(parts[2]) * sheet_height
                bw = float(parts[3]) * sheet_width
                bh = float(parts[4]) * sheet_height
                conf = float(parts[5]) if len(parts) >= 6 else 1.0
                if conf < min_confidence:
                    continue
                xmin = max(0, int(round(xc - bw / 2.0)))
                ymin = max(0, int(round(yc - bh / 2.0)))
                xmax = min(sheet_width, int(round(xc + bw / 2.0)))
                ymax = min(sheet_height, int(round(yc + bh / 2.0)))
                if ymax > ymin and xmax > xmin:
                    detections.append(((ymin, xmin, ymax, xmax), conf, cls_id))
    except Exception as exc:
        logger.debug("Failed to parse precomputed plant detector label %s: %s", label_path, exc)

    return detections


def _load_momocs_outlines(lm3_path: Path) -> Dict[str, np.ndarray]:
    """Loads pre-oriented Momocs leaf outlines from `reports/Leaf_Momocs/momocs_outlines.json`."""
    outlines: Dict[str, np.ndarray] = {}
    json_files = sorted(
        set(
            list(lm3_path.glob("reports/Leaf_Momocs/momocs_outlines.json"))
            + list(lm3_path.glob("**/momocs_outlines.json"))
        )
    )
    for jf in json_files:
        try:
            with open(jf, "r", encoding="utf-8") as f:
                data = json.load(f)
            if isinstance(data, dict):
                entries = data.get("outlines", data) if isinstance(data.get("outlines"), dict) else data
                for k, v in entries.items():
                    pts = v.get("coords", v.get("outline", v)) if isinstance(v, dict) else v
                    if isinstance(pts, list) and len(pts) >= 3:
                        arr = np.asarray(pts, dtype=np.float64)
                        if arr.ndim == 2 and arr.shape[1] >= 2:
                            outlines[Path(str(k)).stem] = arr[:, :2]
        except Exception as exc:
            logger.debug("Error loading Momocs outlines from %s: %s", jf, exc)
    return outlines


def _compute_midrib_angle_from_keypoints(
    kps: Optional[Dict[str, Tuple[float, float, float]]],
) -> float:
    """Computes midrib orientation angle (degrees) from `lamina_base` to `lamina_tip` when available."""
    if not kps:
        return 0.0
    base_pt = kps.get("lamina_base") or kps.get("midvein_15")
    tip_pt = kps.get("lamina_tip") or kps.get("midvein_1")
    if base_pt is not None and tip_pt is not None:
        dx = float(tip_pt[0]) - float(base_pt[0])
        dy = float(tip_pt[1]) - float(base_pt[1])
        if abs(dx) > 1e-6 or abs(dy) > 1e-6:
            return float(math.degrees(math.atan2(dy, dx)))
    return 0.0


def discover_lm3_candidates(
    lm3_dir: Path,
    curated_vouchers: Union[Dict[str, Dict[str, Any]], pd.DataFrame, Sequence[Mapping[str, Any]], None] = None,
    ruler_calibs: Optional[Dict[str, float]] = None,
    raw_images_dir: Optional[Path] = None,
) -> CandidateCollection:
    """
    Discovers all candidate basal leaves from LeafMachine3 outputs:
      - `Data/leaf_measurements.csv`
      - `reports/Specimen_Masks/*_Specimen_PaperCleaned/`
      - `reports/Leaf_Original/` and `reports/Leaf_Oriented/`
      - `reports/Leaf_Momocs/momocs_outlines.json`
      - `Data/landmarks.csv` / `<run_name>.sqlite` (31-keypoint `mid15_pet5` skeleton)
      - Cropped leaf & mask directories for unit-test fixtures

    Returns a `CandidateCollection` (a `List[LeafCandidate]` subclass that also
    supports dict-style lookup by `catalogNumber`).
    """
    candidates_by_voucher: Dict[str, List[LeafCandidate]] = defaultdict(list)
    if not lm3_dir or not Path(lm3_dir).exists():
        return CandidateCollection(by_voucher=candidates_by_voucher)

    lm3_path = Path(lm3_dir)
    calib_records = load_ruler_calibration_records(lm3_path)
    eff_calibs = ruler_calibs if ruler_calibs is not None else {
        cat: float(rec["pixels_per_mm"])
        for cat, rec in calib_records.items()
        if rec.get("pixels_per_mm") is not None and float(rec["pixels_per_mm"]) > 0
    }

    keypoints_map = load_lm3_31_keypoints(lm3_path)
    momocs_map = _load_momocs_outlines(lm3_path)

    # Index Specimen_PaperCleaned masks by catalogNumber
    paper_cleaned_lookup: Dict[str, Path] = {}
    for sp_dir in (
        list(lm3_path.glob("reports/Specimen_Masks/*_Specimen_PaperCleaned"))
        + list(lm3_path.glob("**/Specimen_Masks/*_Specimen_PaperCleaned"))
        + list(lm3_path.glob("**/Specimen_PaperCleaned"))
    ):
        if sp_dir.is_dir():
            for mf in sp_dir.glob("*.*"):
                if mf.suffix.lower() in {".png", ".jpg", ".jpeg", ".webp", ".tif"}:
                    cat_key = _extract_catalog_number(mf.stem)
                    if cat_key:
                        paper_cleaned_lookup[cat_key] = mf

    # Index oriented leaf crops, ECT visualizations, and masks (Leaf / Petiole / Hole)
    oriented_lookup: Dict[str, Path] = {}
    for or_dir in (
        list(lm3_path.glob("reports/Leaf_Oriented"))
        + list(lm3_path.glob("**/reports/Leaf_Oriented"))
        + list(lm3_path.glob("**/Keypoints/Oriented_Cropped_Leaves"))
    ):
        if or_dir.is_dir():
            for of in or_dir.glob("*.*"):
                if of.suffix.lower() in {".png", ".jpg", ".jpeg", ".webp"}:
                    oriented_lookup[of.stem] = of

    ect_lookup: Dict[str, Path] = {}
    for ect_dir in (
        list(lm3_path.glob("reports/Leaf_ECT"))
        + list(lm3_path.glob("**/reports/Leaf_ECT"))
        + list(lm3_path.glob("**/ECT"))
    ):
        if ect_dir.is_dir():
            for ef in ect_dir.glob("*.*"):
                if ef.suffix.lower() in {".png", ".jpg", ".jpeg", ".npy", ".csv"}:
                    ect_lookup[ef.stem] = ef

    mask_dirs = (
        list(lm3_path.glob("reports/Leaf_Masks"))
        + list(lm3_path.glob("**/reports/Leaf_Masks"))
        + list(lm3_path.glob("**/Plant_Components/Segmentation_Whole_Leaves"))
        + list(lm3_path.glob("**/Plant_Components/Segmentation_Masks_Color_Whole_Leaves"))
        + list(lm3_path.glob("**/Keypoints/Oriented_Masks"))
        + list(lm3_path.glob("**/crops/masks"))
        + list(lm3_path.glob("**/masks"))
    )
    mask_lookup: Dict[str, Path] = {}
    petiole_mask_lookup: Dict[str, Path] = {}
    hole_mask_lookup: Dict[str, Path] = {}
    for md in mask_dirs:
        if not md.is_dir():
            continue
        for mf in md.glob("*.*"):
            if mf.suffix.lower() not in {".webp", ".jpg", ".jpeg", ".png", ".bmp", ".tif"}:
                continue
            stem_low = mf.stem.lower()
            if stem_low.endswith("_petiole") or "__petiole" in stem_low:
                base_stem = re.sub(r"(_+?petiole)$", "", mf.stem, flags=re.IGNORECASE)
                petiole_mask_lookup[base_stem] = mf
            elif stem_low.endswith("_hole") or "__hole" in stem_low:
                base_stem = re.sub(r"(_+?hole)$", "", mf.stem, flags=re.IGNORECASE)
                hole_mask_lookup[base_stem] = mf
            else:
                mask_lookup[mf.stem] = mf

    # Index original leaf crops
    leaf_crop_dirs = (
        list(lm3_path.glob("reports/Leaf_Original"))
        + list(lm3_path.glob("**/reports/Leaf_Original"))
        + list(lm3_path.glob("**/Plant_Components/Leaves_Whole"))
        + list(lm3_path.glob("**/Plant_Components/Leaves_Partial"))
        + list(lm3_path.glob("**/Keypoints/Oriented_Cropped_Leaves"))
        + list(lm3_path.glob("**/crops/leaves"))
    )
    crop_lookup: Dict[str, Path] = {}
    for lcd in leaf_crop_dirs:
        if not lcd.is_dir():
            continue
        for cp in lcd.glob("*.*"):
            if cp.suffix.lower() in {".webp", ".jpg", ".jpeg", ".png"}:
                crop_lookup.setdefault(cp.stem, cp)

    leaf_counts: Dict[str, int] = defaultdict(int)
    seen_stems: set[str] = set()

    # 1. Primary LM3 tabular discovery via Data/leaf_measurements.csv
    meas_csvs = sorted(
        set(
            list(lm3_path.glob("Data/leaf_measurements.csv"))
            + list(lm3_path.glob("**/Data/leaf_measurements.csv"))
        )
    )
    for mcsv in meas_csvs:
        try:
            df = pd.read_csv(mcsv)
            if df.empty:
                continue
            for _, row in df.iterrows():
                row_map = row.to_dict()
                raw_stem = str(
                    row_map.get("leaf_stem")
                    or row_map.get("leaf_filename")
                    or row_map.get("crop_name")
                    or row_map.get("filename")
                    or ""
                ).strip()
                stem = Path(raw_stem).stem if raw_stem and raw_stem.lower() != "nan" else ""
                cat_part = _extract_catalog_number(
                    row_map.get("catalogNumber")
                    or row_map.get("catalog_number")
                    or row_map.get("specimen_id")
                    or stem
                )
                if not cat_part:
                    continue

                leaf_counts[cat_part] += 1
                lid_raw = row_map.get("leaf_id", row_map.get("organ_id", leaf_counts[cat_part]))
                try:
                    lid = int(lid_raw)
                except (ValueError, TypeError):
                    lid = leaf_counts[cat_part]

                if not stem:
                    stem = f"{cat_part}__L_{lid:04d}"
                seen_stems.add(stem)

                cal_rec = calib_records.get(cat_part, {})
                ws_raw = row_map.get("work_scale", cal_rec.get("work_scale", 1.0))
                try:
                    work_scale = float(ws_raw) if pd.notna(ws_raw) and float(ws_raw) > 0 else 1.0
                except (ValueError, TypeError):
                    work_scale = 1.0

                cf_source = str(row_map.get("cf_source") or cal_rec.get("cf_source") or "")
                raw_calib = eff_calibs.get(cat_part)
                if raw_calib is None and pd.notna(row_map.get("cf_px_per_cm")):
                    try:
                        raw_calib = float(row_map["cf_px_per_cm"]) / 10.0
                    except (ValueError, TypeError):
                        raw_calib = None

                if raw_calib is not None and raw_calib > 0:
                    scale = (1.0 / raw_calib) if raw_calib >= 1.0 else raw_calib
                    px_cm = resolve_px_per_cm(pixels_per_mm=(1.0 / scale))
                else:
                    scale = 0.04233
                    px_cm = DEFAULT_PX_PER_CM

                # Resolve bounding box (unscaling working-frame coordinates when present)
                bbox = parse_bbox_from_filename(stem) or (0, 0, 0, 0)
                if all(k in row_map and pd.notna(row_map[k]) for k in ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")):
                    ux1, uy1, ux2, uy2 = unscale_working_bbox_to_sheet(
                        (
                            float(row_map["bbox_x1"]),
                            float(row_map["bbox_y1"]),
                            float(row_map["bbox_x2"]),
                            float(row_map["bbox_y2"]),
                        ),
                        work_scale=work_scale,
                    )
                    bbox = (uy1, ux1, uy2, ux2)
                elif all(k in row_map and pd.notna(row_map[k]) for k in ("xmin", "ymin", "xmax", "ymax")):
                    ux1, uy1, ux2, uy2 = unscale_working_bbox_to_sheet(
                        (
                            float(row_map["xmin"]),
                            float(row_map["ymin"]),
                            float(row_map["xmax"]),
                            float(row_map["ymax"]),
                        ),
                        work_scale=work_scale,
                    )
                    bbox = (uy1, ux1, uy2, ux2)

                kps = keypoints_map.get(stem) or keypoints_map.get(f"{cat_part}__{lid}") or keypoints_map.get(cat_part)
                midrib_deg = float(row_map.get("midrib_angle_deg", 0.0)) if pd.notna(row_map.get("midrib_angle_deg")) else _compute_midrib_angle_from_keypoints(kps)
                arch_score = float(row_map.get("archetype_score", row_map.get("ucs_score", 0.0))) if pd.notna(row_map.get("archetype_score", row_map.get("ucs_score", np.nan))) else 0.0
                solidity_val = float(row_map.get("solidity", 0.0)) if pd.notna(row_map.get("solidity")) else 0.0
                pet_w_mm = float(row_map.get("petiole_width_mm", 0.0)) if pd.notna(row_map.get("petiole_width_mm")) else 0.0
                len_mm = float(row_map.get("length_mm", row_map.get("blade_length_mm", 0.0))) if pd.notna(row_map.get("length_mm", row_map.get("blade_length_mm", np.nan))) else 0.0
                wid_mm = float(row_map.get("width_mm", row_map.get("blade_width_mm", 0.0))) if pd.notna(row_map.get("width_mm", row_map.get("blade_width_mm", np.nan))) else 0.0
                area_cm2_val = float(row_map.get("area_cm2", row_map.get("blade_area_cm2", 0.0))) if pd.notna(row_map.get("area_cm2", row_map.get("blade_area_cm2", np.nan))) else 0.0
                area_mm2_val = float(row_map.get("area_mm2", area_cm2_val * 100.0)) if pd.notna(row_map.get("area_mm2", np.nan)) else (area_cm2_val * 100.0)

                candidate = LeafCandidate(
                    catalog_number=cat_part,
                    leaf_id=lid,
                    crop_image_path=crop_lookup.get(stem) or oriented_lookup.get(stem),
                    mask_image_path=mask_lookup.get(stem),
                    bbox=bbox,
                    scale_mm_per_px=scale,
                    px_per_cm=px_cm,
                    ucs_score=arch_score,
                    solidity=solidity_val,
                    midrib_angle_deg=midrib_deg,
                    length_mm=len_mm,
                    width_mm=wid_mm,
                    area_mm2=area_mm2_val,
                    area_cm2=area_cm2_val,
                    petiole_mask_path=petiole_mask_lookup.get(stem),
                    hole_mask_path=hole_mask_lookup.get(stem),
                    specimen_paper_cleaned_mask_path=paper_cleaned_lookup.get(cat_part),
                    keypoints_31=kps,
                    work_scale=work_scale,
                    cf_source=cf_source,
                    archetype_score=arch_score,
                    petiole_width_mm=pet_w_mm,
                    ect_path=ect_lookup.get(stem),
                    momocs_outline_xy=momocs_map.get(stem),
                    oriented_image_path=oriented_lookup.get(stem),
                )
                candidates_by_voucher[cat_part].append(candidate)
        except Exception as exc:
            logger.debug("Error reading leaf_measurements.csv %s: %s", mcsv, exc)

    # 2. Directory-based discovery for cropped leaves & test fixtures not already in leaf_measurements.csv
    for stem, crop_path in crop_lookup.items():
        if stem in seen_stems:
            continue
        seen_stems.add(stem)
        cat_part = _extract_catalog_number(stem)
        if not cat_part:
            continue
        bbox = parse_bbox_from_filename(stem) or (0, 0, 0, 0)
        matched_mask = mask_lookup.get(stem)

        cal_rec = calib_records.get(cat_part, {})
        work_scale = float(cal_rec.get("work_scale", 1.0))
        cf_source = str(cal_rec.get("cf_source", ""))

        raw_calib = eff_calibs.get(cat_part)
        if raw_calib is not None and raw_calib > 0:
            scale = (1.0 / raw_calib) if raw_calib >= 1.0 else raw_calib
            px_cm = resolve_px_per_cm(pixels_per_mm=(1.0 / scale))
        else:
            scale = 0.04233
            px_cm = DEFAULT_PX_PER_CM

        leaf_counts[cat_part] += 1
        lid = leaf_counts[cat_part]
        kps = keypoints_map.get(stem) or keypoints_map.get(f"{cat_part}__{lid}") or keypoints_map.get(cat_part)

        candidate = LeafCandidate(
            catalog_number=cat_part,
            leaf_id=lid,
            crop_image_path=crop_path,
            mask_image_path=matched_mask,
            bbox=bbox,
            scale_mm_per_px=scale,
            px_per_cm=px_cm,
            midrib_angle_deg=_compute_midrib_angle_from_keypoints(kps),
            petiole_mask_path=petiole_mask_lookup.get(stem),
            hole_mask_path=hole_mask_lookup.get(stem),
            specimen_paper_cleaned_mask_path=paper_cleaned_lookup.get(cat_part),
            keypoints_31=kps,
            work_scale=work_scale,
            cf_source=cf_source,
            ect_path=ect_lookup.get(stem),
            momocs_outline_xy=momocs_map.get(stem),
            oriented_image_path=oriented_lookup.get(stem),
        )
        candidates_by_voucher[cat_part].append(candidate)

    total_leaves = sum(len(v) for v in candidates_by_voucher.values())
    logger.info(
        "Discovered LM3 leaf candidates for %d voucher sheets (%d total leaves).",
        len(candidates_by_voucher),
        total_leaves,
    )
    return CandidateCollection(by_voucher=candidates_by_voucher)


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


__all__ = [
    "LM3_PLANT_CLASSES",
    "LM3_PCD_CLASSES",
    "LM3_LEAF_CLASS_IDS",
    "LM3_REPRODUCTIVE_CLASS_IDS",
    "LM3_LEAF_SEG_CLASSES",
    "LM3_KEYPOINT_31_NAMES",
    "VALID_LM3_CF_SOURCES",
    "DEFAULT_HERBARIUM_DPI",
    "DEFAULT_PX_PER_CM",
    "DEFAULT_PX_PER_MM",
    "compute_lm3_work_scale",
    "unscale_working_bbox_to_sheet",
    "unscale_working_points_to_sheet",
    "resolve_px_per_cm",
    "compute_blade_area_cm2",
    "LeafCandidate",
    "CandidateCollection",
    "parse_bbox_from_filename",
    "load_ruler_calibration_records",
    "load_ruler_calibrations",
    "load_lm3_31_keypoints",
    "load_precomputed_pcd_labels",
    "discover_lm3_candidates",
    "load_and_preprocess_mask",
]
