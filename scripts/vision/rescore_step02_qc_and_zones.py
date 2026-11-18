#!/usr/bin/env python3
"""
===============================================================================
Script: rescore_step02_qc_and_zones.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Fast CPU post-Step-02 re-scoring utility that refreshes Tier A leaf QC flags
    (`detector_class_id`, `is_crop_clipped`, `hole_area_frac`, `enclosed_void_frac`,
    `bilateral_lobe_symmetry`, `efa_qc_pass`, `qc_flag_reason`), void-aware
    `solidity` / `circularity`, hairline sinus bridge contours, and anatomical
    `basal` vs. `cauline` leaf zones across `extracted_leaves` in SQLite and
    `data/tables/extracted_leaf_manifest.csv` without re-running GPU inference.
===============================================================================
"""

from __future__ import annotations

import argparse
import inspect
import json
import logging
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig
from scripts.core.db import get_db_connection
from scripts.core.io_utils import atomic_write_text
from scripts.vision.leaf_geometry_utils import (
    compute_lm3_shape_metrics,
    export_standardized_contour,
)
from scripts.vision.lm3_data_loader import _CALIBRATED_SCALE_SOURCES
from scripts.vision.lm3_inference_engines import (
    LM3ScaleCalibrator,
    evaluate_leaf_damage_and_fold_qc,
    solve_ruler_tick_lattice_px_per_mm,
)
from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical
from scripts.vision.step02_photometry_corymb_utils import (
    compute_calibrated_blade_photometry,
    estimate_sheet_paper_cie_l_ref,
    harmonize_upright_mask_sinus_bridges,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("RescoreStep02")


def _resolve_mask_file(raw_mask_path: Any, data_dir: Path) -> Optional[Path]:
    s = str(raw_mask_path or "").strip()
    if not s or s.lower() in {"nan", "none", "na"}:
        return None
    p = Path(s)
    if p.is_file():
        return p
    for cand in (
        data_dir / "masks" / p.name,
        PROJECT_ROOT / s,
        PROJECT_ROOT / "data" / "masks" / p.name,
    ):
        if cand.is_file():
            return cand
    return None


def _resolve_voucher_image_file(raw_image_path: Any, catalog_number: str, data_dir: Path) -> Optional[Path]:
    s = str(raw_image_path or "").strip()
    candidates: List[Path] = []
    if s and s.lower() not in {"nan", "none", "na"}:
        p = Path(s)
        candidates.extend([p, PROJECT_ROOT / s, data_dir / "raw_vouchers" / p.name])
    cat_safe = str(catalog_number or "").strip()
    if cat_safe:
        for ext in (".webp", ".jpg", ".jpeg", ".png"):
            candidates.append(data_dir / "raw_vouchers" / f"{cat_safe}{ext}")
            candidates.append(PROJECT_ROOT / "data" / "raw_vouchers" / f"{cat_safe}{ext}")
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def _eval_single_leaf_qc_worker(task: Tuple[Dict[str, Any], str, bool]) -> Dict[str, Any]:
    """Evaluate hairline sinus bridges, void-aware solidity/circ, and EFA QC gates for 1 leaf."""
    row, data_dir_str, seg_partial = task
    eff_data_dir = Path(data_dir_str)
    m_file = _resolve_mask_file(row.get("mask_path"), eff_data_dir)
    mask_u8: Optional[np.ndarray] = None
    if m_file is not None:
        m_raw = cv2.imread(str(m_file), cv2.IMREAD_GRAYSCALE)
        if m_raw is not None and np.count_nonzero(m_raw) > 0:
            mask_u8 = ((m_raw > 0).astype(np.uint8)) * 255

    def _fval(key: str, default: float) -> float:
        v = row.get(key)
        return float(v) if v is not None and pd.notna(v) else default

    def _ival(key: str, default: int) -> int:
        v = row.get(key)
        return int(v) if v is not None and pd.notna(v) else default

    sol_val = _fval("solidity", 0.85)
    circ_val = _fval("circularity", float("nan"))
    bd_val = _fval("bilateral_dice", 0.85)
    sia_val = _fval("si_a", 0.0)
    gp_val = bool(_ival("lm3_gates_pass", 1))
    ncc_val = _ival("n_lamina_components", 1)
    dcid_val = _ival("detector_class_id", 0)
    haf_val = _fval("hole_area_frac", 0.0)
    clipped_val = bool(_ival("is_crop_clipped", 0))
    evf_val = 0.0
    new_cjson: Optional[str] = None
    cal_src_set = {str(s).strip().lower() for s in _CALIBRATED_SCALE_SOURCES}
    ss_row = str(row.get("scale_source") or "").strip().lower()
    raw_px_row = row.get("pixels_per_mm")
    if raw_px_row is None or pd.isna(raw_px_row):
        raw_px_row = row.get("voucher_px_per_mm")
    is_cal_row = (
        (_ival("is_calibrated_mm", 0) == 1 or _ival("voucher_is_cal", 0) == 1)
        and (not ss_row or ss_row in cal_src_set)
    )
    px_mm_eval: Optional[float] = (
        float(raw_px_row)
        if (
            is_cal_row
            and raw_px_row is not None
            and pd.notna(raw_px_row)
            and math.isfinite(float(raw_px_row))
            and float(raw_px_row) > 0.0
        )
        else None
    )

    if mask_u8 is not None:
        hm, wm = mask_u8.shape[:2]
        edge_allow = max(6, int(round(0.06 * max(hm, wm))))
        if (
            int(np.count_nonzero(mask_u8[0, :])) > edge_allow
            or int(np.count_nonzero(mask_u8[-1, :])) > edge_allow
            or int(np.count_nonzero(mask_u8[:, 0])) > edge_allow
            or int(np.count_nonzero(mask_u8[:, -1])) > edge_allow
        ):
            clipped_val = True

        orig_u8 = mask_u8
        cleaned_u8, evf_val = harmonize_upright_mask_sinus_bridges(
            orig_u8,
            pixels_per_mm=px_mm_eval,
        )
        was_opened = not np.array_equal(cleaned_u8, orig_u8)
        mask_u8 = cleaned_u8
        shape_m = compute_lm3_shape_metrics(mask_u8)
        if shape_m["area_px"] > 0.0:
            sol_val = round(float(shape_m["solidity"]), 4)
            circ_val = round(float(shape_m["circularity"]), 4)
        if was_opened:
            c_res = export_standardized_contour(
                mask_u8,
                str(row["catalogNumber"]),
                int(row["leaf_id"]),
                eff_data_dir,
                return_coords=True,
                already_oriented=bool(_ival("mask_is_upright", 1) or 1),
                bridge_already_opened=True,
            )
            if isinstance(c_res, tuple) and c_res[1]:
                new_cjson = json.dumps(c_res[1], separators=(",", ":"))

    qc_res = evaluate_leaf_damage_and_fold_qc(
        oriented_mask_u8=mask_u8,
        solidity=sol_val,
        bilateral_dice=bd_val,
        si_a=sia_val,
        lm3_gates_pass=gp_val,
        n_lamina_components=ncc_val,
        detector_class_id=dcid_val,
        segment_partial_leaves=seg_partial,
        is_crop_clipped=clipped_val,
        hole_area_frac=haf_val,
        enclosed_void_frac=evf_val,
        is_geometric_fallback_pose=True,
    )
    return {
        "solidity": sol_val,
        "circularity": circ_val,
        "enclosed_void_frac": round(float(evf_val), 4),
        "hole_area_frac": round(float(haf_val), 4),
        "is_crop_clipped": 1 if clipped_val else 0,
        "detector_class_id": dcid_val,
        "bilateral_lobe_symmetry": 1 if qc_res["bilateral_lobe_symmetry"] else 0,
        "efa_qc_pass": 1 if qc_res["efa_qc_pass"] else 0,
        "qc_flag_reason": str(qc_res["qc_flag_reason"]),
        "contour_xy_json": new_cjson,
        "min_rms": float(qc_res["min_margin_rms_norm"]),
    }


def _calibrate_voucher_ruler_worker(task: Tuple[str, str, Optional[float]]) -> Optional[Tuple[str, float]]:
    """Run CPU unboxed margin metric tick-lattice scan for one uncalibrated voucher sheet."""
    cat, img_path_str, mp_prior_mm = task
    img = cv2.imread(img_path_str, cv2.IMREAD_COLOR)
    if img is None or img.size == 0:
        return None
    if "verified_ruler_crop" in inspect.signature(solve_ruler_tick_lattice_px_per_mm).parameters:
        px_per_mm = solve_ruler_tick_lattice_px_per_mm(
            img, mp_prior_px_per_mm=mp_prior_mm, verified_ruler_crop=False
        )
    else:
        px_per_mm = solve_ruler_tick_lattice_px_per_mm(img, mp_prior_px_per_mm=mp_prior_mm)
    if px_per_mm is None or not np.isfinite(float(px_per_mm)) or float(px_per_mm) <= 1.0:
        return None
    return (cat, round(float(px_per_mm), 4))


def _resolve_crop_file(raw_crop_path: Any, data_dir: Path) -> Optional[Path]:
    s = str(raw_crop_path or "").strip()
    if not s or s.lower() in {"nan", "none", "na"}:
        return None
    p = Path(s)
    for cand in (p, PROJECT_ROOT / s, data_dir / "crops" / p.name, PROJECT_ROOT / "data" / "crops" / p.name):
        if cand.is_file():
            return cand
    return None


def rescore_tier_a_leaves(
    db_path: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    manifest_csv: Optional[Path] = None,
    segment_partial_leaves: Optional[bool] = None,
    recompute_photometry: bool = False,
    workers: Optional[int] = None,
) -> Dict[str, int]:
    """Re-evaluates leaf QC gates, hairline sinus bridges, ruler scales, and basal/cauline zones."""
    cfg = PipelineConfig.from_yaml()
    eff_data_dir = (Path(data_dir) if data_dir is not None else (PROJECT_ROOT / "data")).resolve()
    eff_manifest = (
        Path(manifest_csv).resolve()
        if manifest_csv is not None
        else (eff_data_dir / "tables" / "extracted_leaf_manifest.csv")
    )
    seg_partial = (
        bool(segment_partial_leaves)
        if segment_partial_leaves is not None
        else bool(getattr(cfg.segmentation, "segment_partial_leaves", False))
    )

    conn = get_db_connection(db_path=db_path, init_schema=True)
    try:
        el_df = pd.read_sql_query(
            """
            SELECT
                el.*,
                pi.crown_x, pi.crown_y,
                pi.rosette_centroid_x, pi.rosette_centroid_y,
                pi.plant_height_px, pi.is_bent_stem,
                v.image_path, v.pixels_per_mm AS voucher_px_per_mm,
                v.is_calibrated_mm AS voucher_is_cal,
                COALESCE(vrp.capitula_count, vrp0.capitula_count, 0) AS eff_capitula_count
            FROM extracted_leaves el
            LEFT JOIN vouchers v
                ON v.catalogNumber = el.catalogNumber
            LEFT JOIN plant_individuals pi
                ON pi.catalogNumber = el.catalogNumber
               AND pi.plant_individual_id = el.plant_individual_id
            LEFT JOIN voucher_reproductive_phenology vrp
                ON vrp.catalogNumber = el.catalogNumber
               AND vrp.plant_individual_id = el.plant_individual_id
            LEFT JOIN voucher_reproductive_phenology vrp0
                ON vrp0.catalogNumber = el.catalogNumber
               AND vrp0.plant_individual_id = 0
            ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;
            """,
            conn,
        )
        if el_df.empty:
            logger.info("No rows in extracted_leaves to re-score.")
            return {
                "total": 0,
                "efa_qc_pass": 0,
                "flagged": 0,
                "basal": 0,
                "cauline": 0,
                "ruler_calibrated_vouchers": 0,
                "tick_lattice_unboxed_vouchers": 0,
                "photometry_recomputed": 0,
            }

        n_rows = len(el_df)
        n_workers = min(32, max(1, int(workers) if workers is not None else (os.cpu_count() or 8)))
        cal_src_set = {str(s).strip().lower() for s in _CALIBRATED_SCALE_SOURCES}
        if "scale_source" not in el_df.columns:
            el_df["scale_source"] = ""
        if "is_calibrated_mm" not in el_df.columns:
            el_df["is_calibrated_mm"] = 0
        if "pixels_per_mm" not in el_df.columns:
            el_df["pixels_per_mm"] = np.nan

        # CPU ruler-calibration backfill (LM3ScaleCalibrator archival YOLO26x + verified tick lattice,
        # falling back to unboxed outer-margin tick-lattice scan) for uncalibrated vouchers
        ruler_tasks: List[Tuple[str, str, Optional[float]]] = []
        seen_cats: set[str] = set()
        for _, rv in el_df[["catalogNumber", "image_path", "voucher_px_per_mm", "voucher_is_cal"]].drop_duplicates("catalogNumber").iterrows():
            cat_s = str(rv["catalogNumber"]).strip()
            if not cat_s or cat_s in seen_cats:
                continue
            seen_cats.add(cat_s)
            is_cal_v = int(rv["voucher_is_cal"]) if pd.notna(rv.get("voucher_is_cal")) else 0
            cat_mask = el_df["catalogNumber"].astype(str).str.strip() == cat_s
            l_ss = el_df.loc[cat_mask, "scale_source"].fillna("").astype(str).str.strip().str.lower()
            l_cal = pd.to_numeric(el_df.loc[cat_mask, "is_calibrated_mm"], errors="coerce").fillna(0).astype(int)
            leaf_cal_ok = bool(((l_cal == 1) & (l_ss.isin(cal_src_set) | (l_ss == ""))).any())
            if is_cal_v != 0 and leaf_cal_ok:
                continue
            img_f = _resolve_voucher_image_file(rv.get("image_path"), cat_s, eff_data_dir)
            if img_f is None:
                continue
            v_px = float(rv["voucher_px_per_mm"]) if pd.notna(rv.get("voucher_px_per_mm")) else None
            ruler_tasks.append((cat_s, str(img_f), v_px))

        ruler_cal_scales: Dict[str, float] = {}
        unboxed_scales: Dict[str, float] = {}
        if ruler_tasks:
            calibrator: Optional[LM3ScaleCalibrator] = None
            try:
                calibrator = LM3ScaleCalibrator(device="cpu", use_cf_predicted_by_mp=False)
            except Exception as cal_exc:
                logger.warning(f"LM3ScaleCalibrator CPU init skipped ({cal_exc}); using unboxed tick-lattice scan.")
                calibrator = None

            if calibrator is not None:
                def _cal_sheet_task(
                    task: Tuple[str, str, Optional[float]],
                ) -> Tuple[str, str, Optional[float], bool, Tuple[str, str, Optional[float]]]:
                    cat_t, img_p_t, _v_px_t = task
                    img_c = cv2.imread(img_p_t, cv2.IMREAD_COLOR)
                    if img_c is None or img_c.size == 0:
                        return (cat_t, "uncalibrated", None, False, task)
                    try:
                        res_c = calibrator.calibrate_sheet(
                            img_c, catalog_number=cat_t, use_cf_predicted_by_mp=False
                        )
                        ss_c = str(res_c.get("scale_source") or "").strip().lower()
                        cal_c = bool(res_c.get("is_calibrated_mm"))
                        px_c = res_c.get("pixels_per_mm")
                        if px_c is not None and np.isfinite(float(px_c)) and float(px_c) > 1.0:
                            return (cat_t, ss_c, round(float(px_c), 4), cal_c, task)
                        return (cat_t, ss_c, None, False, task)
                    except Exception:
                        return (cat_t, "error", None, False, task)

                if (
                    len(ruler_tasks) >= 8
                    and n_workers > 1
                    and getattr(calibrator, "archival_session", None) is None
                    and "calibrate_sheet" not in getattr(calibrator, "__dict__", {})
                    and "_detect_ruler_boxes" not in getattr(calibrator, "__dict__", {})
                    and getattr(type(calibrator), "calibrate_sheet", None) is LM3ScaleCalibrator.calibrate_sheet
                ):
                    with ProcessPoolExecutor(max_workers=n_workers) as pool:
                        hits = list(pool.map(_calibrate_voucher_ruler_worker, ruler_tasks, chunksize=8))
                    for h in hits:
                        if h is not None:
                            unboxed_scales[h[0]] = h[1]
                else:
                    if len(ruler_tasks) >= 8 and n_workers > 1:
                        with ThreadPoolExecutor(max_workers=min(4, n_workers)) as tpool:
                            cal_hits = list(tpool.map(_cal_sheet_task, ruler_tasks))
                    else:
                        cal_hits = [_cal_sheet_task(t) for t in ruler_tasks]
                    fallback_tasks: List[Tuple[str, str, Optional[float]]] = []
                    for cat_r, ss_r, px_r, is_cal_r, orig_task in cal_hits:
                        if ss_r == "ruler_calibrated" and is_cal_r and px_r is not None:
                            ruler_cal_scales[cat_r] = px_r
                        elif ss_r == "tick_lattice_unboxed" and px_r is not None:
                            unboxed_scales[cat_r] = px_r
                        elif ss_r == "error":
                            fallback_tasks.append(orig_task)
                    for fb_t in fallback_tasks:
                        hit_fb = _calibrate_voucher_ruler_worker(fb_t)
                        if hit_fb is not None:
                            unboxed_scales[hit_fb[0]] = hit_fb[1]
            else:
                if len(ruler_tasks) >= 8 and n_workers > 1:
                    with ProcessPoolExecutor(max_workers=n_workers) as pool:
                        hits = list(pool.map(_calibrate_voucher_ruler_worker, ruler_tasks, chunksize=8))
                else:
                    hits = [_calibrate_voucher_ruler_worker(t) for t in ruler_tasks]
                for h in hits:
                    if h is not None:
                        unboxed_scales[h[0]] = h[1]
            logger.info(
                f"CPU scale backfill: ruler_calibrated={len(ruler_cal_scales)}, "
                f"tick_lattice_unboxed={len(unboxed_scales)} / {len(ruler_tasks)} uncalibrated vouchers."
            )

        if ruler_cal_scales:
            cat_series = el_df["catalogNumber"].astype(str).str.strip()
            rc_mask = cat_series.isin(ruler_cal_scales)
            for idx_c in el_df.index[rc_mask]:
                cat_c = str(el_df.at[idx_c, "catalogNumber"]).strip()
                new_px = float(ruler_cal_scales[cat_c])
                el_df.at[idx_c, "pixels_per_mm"] = new_px
                el_df.at[idx_c, "is_calibrated_mm"] = 1
                el_df.at[idx_c, "scale_source"] = "ruler_calibrated"
                if "voucher_px_per_mm" in el_df.columns:
                    el_df.at[idx_c, "voucher_px_per_mm"] = new_px
                if "voucher_is_cal" in el_df.columns:
                    el_df.at[idx_c, "voucher_is_cal"] = 1
                apx = pd.to_numeric(el_df.at[idx_c, "area_px"], errors="coerce") if "area_px" in el_df.columns else np.nan
                if pd.notna(apx) and float(apx) > 0.0:
                    el_df.at[idx_c, "area_mm2"] = round(float(apx) / (new_px ** 2), 2)
                blpx = pd.to_numeric(el_df.at[idx_c, "blade_length_px"], errors="coerce") if "blade_length_px" in el_df.columns else np.nan
                if pd.notna(blpx) and float(blpx) > 0.0:
                    el_df.at[idx_c, "blade_length_mm"] = round(float(blpx) / new_px, 3)
                bwpx = pd.to_numeric(el_df.at[idx_c, "blade_width_px"], errors="coerce") if "blade_width_px" in el_df.columns else np.nan
                if (pd.isna(bwpx) or float(bwpx) <= 0.0) and pd.notna(blpx) and float(blpx) > 0.0 and "aspect_ratio" in el_df.columns:
                    ar_v = pd.to_numeric(el_df.at[idx_c, "aspect_ratio"], errors="coerce")
                    if pd.notna(ar_v) and float(ar_v) > 1e-6:
                        bwpx = float(blpx) / float(ar_v)
                if pd.notna(bwpx) and float(bwpx) > 0.0:
                    el_df.at[idx_c, "blade_width_mm"] = round(float(bwpx) / new_px, 3)
                plpx = pd.to_numeric(el_df.at[idx_c, "petiole_length_px"], errors="coerce") if "petiole_length_px" in el_df.columns else np.nan
                if pd.notna(plpx) and float(plpx) >= 0.0:
                    el_df.at[idx_c, "petiole_length_mm"] = round(float(plpx) / new_px, 3)
                pwpx = pd.to_numeric(el_df.at[idx_c, "petiole_width_px"], errors="coerce") if "petiole_width_px" in el_df.columns else np.nan
                if pd.notna(pwpx) and float(pwpx) > 0.0:
                    el_df.at[idx_c, "petiole_width_mm"] = round(float(pwpx) / new_px, 4)

        if unboxed_scales:
            cat_series = el_df["catalogNumber"].astype(str).str.strip()
            unb_mask = cat_series.isin(unboxed_scales) & (~cat_series.isin(ruler_cal_scales))
            for idx_c in el_df.index[unb_mask]:
                cat_c = str(el_df.at[idx_c, "catalogNumber"]).strip()
                new_px = float(unboxed_scales[cat_c])
                el_df.at[idx_c, "pixels_per_mm"] = new_px
                el_df.at[idx_c, "is_calibrated_mm"] = 0
                el_df.at[idx_c, "scale_source"] = "tick_lattice_unboxed"

        # Normalize scale_source across calibrated, unboxed tick-lattice, and mp_estimated rows
        is_cal_arr = pd.to_numeric(el_df["is_calibrated_mm"], errors="coerce").fillna(0).astype(int)
        px_arr = pd.to_numeric(el_df["pixels_per_mm"], errors="coerce")
        ss_cur = el_df["scale_source"].fillna("").astype(str).str.strip()
        ss_low = ss_cur.str.lower()
        is_cal_arr = (
            (is_cal_arr == 1)
            & px_arr.notna()
            & (px_arr > 0.0)
            & (ss_low.isin(cal_src_set) | (ss_low == ""))
        ).astype(int)
        el_df["is_calibrated_mm"] = is_cal_arr
        el_df.loc[(is_cal_arr == 1) & (~ss_low.isin(cal_src_set)), "scale_source"] = "ruler_calibrated"
        el_df.loc[ss_low == "tick_lattice_unboxed", "scale_source"] = "tick_lattice_unboxed"
        el_df.loc[ss_low == "tick_lattice_unboxed", "is_calibrated_mm"] = 0
        mp_fix = (
            (el_df["is_calibrated_mm"] == 0)
            & px_arr.notna()
            & (px_arr > 0.0)
            & (~ss_low.isin(cal_src_set | {"mp_estimated", "tick_lattice_unboxed"}))
        )
        el_df.loc[mp_fix, "scale_source"] = "mp_estimated"

        qc_cols = [
            "catalogNumber", "leaf_id", "mask_path", "solidity", "circularity",
            "bilateral_dice", "si_a", "lm3_gates_pass", "n_lamina_components",
            "detector_class_id", "hole_area_frac", "is_crop_clipped", "mask_is_upright",
            "pixels_per_mm", "is_calibrated_mm", "scale_source", "voucher_px_per_mm", "voucher_is_cal",
        ]
        avail_qc_cols = [c for c in qc_cols if c in el_df.columns]
        qc_tasks = [
            (rec, str(eff_data_dir), seg_partial)
            for rec in el_df[avail_qc_cols].to_dict(orient="records")
        ]
        if n_rows >= 64 and n_workers > 1:
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                qc_results = list(pool.map(_eval_single_leaf_qc_worker, qc_tasks, chunksize=64))
        else:
            qc_results = [_eval_single_leaf_qc_worker(t) for t in qc_tasks]

        el_df["solidity"] = [r["solidity"] for r in qc_results]
        el_df["circularity"] = [r["circularity"] for r in qc_results]
        el_df["enclosed_void_frac"] = [r["enclosed_void_frac"] for r in qc_results]
        el_df["hole_area_frac"] = [r["hole_area_frac"] for r in qc_results]
        el_df["is_crop_clipped"] = [r["is_crop_clipped"] for r in qc_results]
        el_df["detector_class_id"] = [r["detector_class_id"] for r in qc_results]
        el_df["bilateral_lobe_symmetry"] = [r["bilateral_lobe_symmetry"] for r in qc_results]
        el_df["efa_qc_pass"] = [r["efa_qc_pass"] for r in qc_results]
        el_df["qc_flag_reason"] = [r["qc_flag_reason"] for r in qc_results]
        el_df["_min_rms"] = [r["min_rms"] for r in qc_results]
        upd_contour_json: List[Optional[str]] = [r["contour_xy_json"] for r in qc_results]

        # Re-evaluate basal vs. cauline leaf zones per (catalogNumber, plant_individual_id)
        crown_backfills: List[Tuple[float, float, str, int]] = []
        for (_cat, _pid), grp in el_df.groupby(["catalogNumber", "plant_individual_id"], sort=False):
            if len(grp) == 0:
                continue
            r0 = grp.iloc[0]
            cx_cr = float(r0["crown_x"]) if pd.notna(r0.get("crown_x")) else None
            cy_cr = float(r0["crown_y"]) if pd.notna(r0.get("crown_y")) else None
            if cx_cr is None or cy_cr is None:
                rx = float(r0["rosette_centroid_x"]) if pd.notna(r0.get("rosette_centroid_x")) else None
                ry = float(r0["rosette_centroid_y"]) if pd.notna(r0.get("rosette_centroid_y")) else None
                if rx is not None and ry is not None:
                    cx_cr, cy_cr = rx, ry
                    crown_backfills.append((cx_cr, cy_cr, str(_cat), int(_pid)))

            has_caps = int(r0.get("eff_capitula_count") or 0) > 0
            med_bld = float(grp["blade_length_px"].median()) if grp["blade_length_px"].notna().any() else 90.0
            p_height = (
                float(r0["plant_height_px"])
                if pd.notna(r0.get("plant_height_px")) and float(r0["plant_height_px"]) > 20.0
                else (max(320.0, 3.8 * max(med_bld, 40.0)) if len(grp) >= 2 else None)
            )
            cx_eval = cx_cr if cx_cr is not None else 1500.0
            cy_eval = cy_cr if cy_cr is not None else 2200.0

            # Reconstruct anatomical frame whenever >= 2 leaves exist (flowering OR sterile)
            if len(grp) >= 2 and p_height is not None:
                leaf_recs: List[Dict[str, Any]] = []
                n_prox_pet = 0
                for _, r_l in grp.iterrows():
                    lid_i = int(r_l["leaf_id"])
                    bld_px = float(r_l["blade_length_px"]) if pd.notna(r_l.get("blade_length_px")) and float(r_l["blade_length_px"]) > 0 else 80.0
                    bld_w_px = float(r_l["blade_width_px"]) if pd.notna(r_l.get("blade_width_px")) and float(r_l["blade_width_px"]) > 0 else np.nan
                    pet_px = float(r_l["petiole_length_px"]) if pd.notna(r_l.get("petiole_length_px")) and float(r_l["petiole_length_px"]) >= 0 else 0.0
                    pet_w_px = float(r_l["petiole_width_px"]) if pd.notna(r_l.get("petiole_width_px")) and float(r_l["petiole_width_px"]) >= 0 else np.nan
                    pbr_v = float(r_l["petiole_to_blade_ratio"]) if pd.notna(r_l.get("petiole_to_blade_ratio")) else (pet_px / max(bld_px, 1.0))
                    dc_prev = float(r_l["dist_from_crown_norm"]) if pd.notna(r_l.get("dist_from_crown_norm")) else np.nan
                    cy_rel = float(r_l["cy_rel_plant_bbox"]) if pd.notna(r_l.get("cy_rel_plant_bbox")) else np.nan
                    if np.isfinite(cy_rel):
                        s_est = max(0.0, min(1.0, 1.0 - cy_rel))
                    elif np.isfinite(dc_prev):
                        s_est = max(0.0, min(1.0, dc_prev))
                    else:
                        s_est = 0.0
                    y_c = cy_eval - s_est * p_height
                    y_base = min(cy_eval, y_c + 0.5 * bld_px)
                    y_apex = y_c - 0.5 * bld_px
                    y_tip = min(cy_eval, y_base + min(pet_px, 0.45 * p_height))
                    has_pet_m = bool(pet_px > 0.0 or pbr_v >= 0.12)
                    if has_pet_m and s_est <= 0.45:
                        n_prox_pet += 1
                    leaf_recs.append({
                        "leaf_id": lid_i,
                        "plant_individual_id": int(_pid),
                        "p_base_sheet": (cx_eval, y_base),
                        "p_apex_sheet": (cx_eval, y_apex),
                        "p_petiole_tip": (cx_eval, y_tip),
                        "blade_length_px": bld_px,
                        "blade_width_px": bld_w_px,
                        "petiole_length_px": pet_px,
                        "petiole_width_px": pet_w_px,
                        "petiole_to_blade_ratio": pbr_v,
                        "has_petiole_mask": has_pet_m,
                        "cy_rel_plant_bbox": cy_rel if np.isfinite(cy_rel) else (1.0 - s_est),
                        "bilateral_dice": float(r_l["bilateral_dice"]) if pd.notna(r_l.get("bilateral_dice")) else np.nan,
                        "si_a": float(r_l["si_a"]) if pd.notna(r_l.get("si_a")) else np.nan,
                        "lm3_gates_pass": bool(int(r_l["lm3_gates_pass"]) != 0) if pd.notna(r_l.get("lm3_gates_pass")) else True,
                    })
                has_rosette_root = bool((cx_cr is not None and cy_cr is not None and n_prox_pet >= 1) or n_prox_pet >= 2)
                cap_box = (
                    [(cy_eval - p_height - 20.0, cx_eval - 20.0, cy_eval - p_height + 20.0, cx_eval + 20.0)]
                    if has_caps
                    else []
                )
                root_box = (
                    [(cy_eval - 20.0, cx_eval - 20.0, cy_eval + 20.0, cx_eval + 20.0)]
                    if has_rosette_root
                    else []
                )
                z_res = classify_leaf_zone_anatomical(
                    leaf_recs,
                    crowns_by_pid={int(_pid): (cx_eval, cy_eval)},
                    capitula_detections=cap_box,
                    root_boxes=root_box,
                    is_bent_stem=bool(int(r0.get("is_bent_stem") or 0)),
                )
                for _, r_l in grp.iterrows():
                    lid_i = int(r_l["leaf_id"])
                    if lid_i in z_res:
                        mask_r = (el_df["catalogNumber"] == _cat) & (el_df["plant_individual_id"] == _pid) & (el_df["leaf_id"] == lid_i)
                        el_df.loc[mask_r, "leaf_zone"] = z_res[lid_i]["leaf_zone"]
                        el_df.loc[mask_r, "leaf_zone_confidence"] = z_res[lid_i]["leaf_zone_confidence"]
                        el_df.loc[mask_r, "dist_from_crown_norm"] = z_res[lid_i]["dist_from_crown_norm"]

            # Within-rosette sibling aspect-ratio half-fold screen on basal leaves
            grp_after = el_df.loc[grp.index]
            bas_grp = grp_after[grp_after["leaf_zone"].astype(str).str.lower() == "basal"]
            if len(bas_grp) >= 3:
                clean_b = bas_grp[bas_grp["efa_qc_pass"] == 1]
                ref_b = clean_b if len(clean_b) >= 2 else bas_grp
                med_ar = float(ref_b["aspect_ratio"].median())
                med_area = float(ref_b["area_px"].median())
                for idx_b, rb in bas_grp.iterrows():
                    ar_b = float(rb["aspect_ratio"]) if pd.notna(rb.get("aspect_ratio")) else med_ar
                    area_b = float(rb["area_px"]) if pd.notna(rb.get("area_px")) else med_area
                    rms_b = float(rb["_min_rms"])
                    bd_b = float(rb["bilateral_dice"]) if pd.notna(rb.get("bilateral_dice")) else 1.0
                    if (
                        ar_b > 1.85 * max(med_ar, 0.5)
                        and area_b < 0.65 * max(med_area, 1.0)
                        and ((np.isfinite(rms_b) and rms_b < 0.028) or bd_b < 0.78)
                    ):
                        toks = [t for t in str(el_df.at[idx_b, "qc_flag_reason"]).split(";") if t and t != "PASS"]
                        if "LONGITUDINAL_MARGIN_CREASE" not in toks:
                            toks.append("LONGITUDINAL_MARGIN_CREASE")
                        el_df.at[idx_b, "efa_qc_pass"] = 0
                        el_df.at[idx_b, "qc_flag_reason"] = ";".join(toks)

        for _pcol in (
            "blade_cie_l_norm",
            "blade_cie_l_sd",
            "blade_cie_l_sheet_delta",
            "blade_cie_a_delta",
            "blade_cie_b_delta",
        ):
            if _pcol not in el_df.columns:
                el_df[_pcol] = np.nan

        upd_cie_l_norm: List[Optional[float]] = [None] * n_rows
        upd_cie_l_sd: List[Optional[float]] = [None] * n_rows
        upd_cie_l_sheet_delta: List[Optional[float]] = [None] * n_rows
        upd_cie_a_delta: List[Optional[float]] = [None] * n_rows
        upd_cie_b_delta: List[Optional[float]] = [None] * n_rows
        n_photo_upd = 0
        n_sheet_delta_rec = 0
        null_sheet_delta_mask = pd.to_numeric(el_df["blade_cie_l_sheet_delta"], errors="coerce").isna()
        if recompute_photometry or bool(null_sheet_delta_mask.any()):
            by_cat_photo: Dict[str, List[Tuple[int, Any, Any, Any, Optional[float], bool, Optional[float]]]] = {}
            for pos_p, (_idx_p, row_p) in enumerate(el_df.iterrows()):
                if not recompute_photometry and not bool(null_sheet_delta_mask.iloc[pos_p]):
                    continue
                cat_p = str(row_p["catalogNumber"]).strip()
                is_cal_p = int(row_p.get("is_calibrated_mm") or 0) == 1
                px_raw_p = row_p.get("pixels_per_mm")
                px_mm_p = (
                    float(px_raw_p)
                    if (
                        is_cal_p
                        and px_raw_p is not None
                        and pd.notna(px_raw_p)
                        and math.isfinite(float(px_raw_p))
                        and float(px_raw_p) > 0.0
                    )
                    else None
                )
                is_up_p = True if pd.isna(row_p.get("mask_is_upright")) else bool(int(row_p.get("mask_is_upright") or 0))
                ex_ln_raw = row_p.get("blade_cie_l_norm")
                ex_ln_p = (
                    float(ex_ln_raw)
                    if (ex_ln_raw is not None and pd.notna(ex_ln_raw) and math.isfinite(float(ex_ln_raw)))
                    else None
                )
                by_cat_photo.setdefault(cat_p, []).append((
                    pos_p,
                    row_p.get("image_path"),
                    row_p.get("mask_path"),
                    row_p.get("crop_path"),
                    px_mm_p,
                    is_up_p,
                    ex_ln_p,
                ))

            def _photo_sheet_worker(
                item: Tuple[str, List[Tuple[int, Any, Any, Any, Optional[float], bool, Optional[float]]]],
            ) -> List[Tuple[int, Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]]:
                cat_w, recs_w = item
                img_f = _resolve_voucher_image_file(recs_w[0][1], cat_w, eff_data_dir)
                img_bgr = cv2.imread(str(img_f), cv2.IMREAD_COLOR) if img_f is not None else None
                sh_hw = img_bgr.shape[:2] if (img_bgr is not None and img_bgr.ndim == 3 and img_bgr.size > 0) else None
                paper_ref: Optional[float] = estimate_sheet_paper_cie_l_ref(img_bgr) if sh_hw is not None else None
                out_w: List[
                    Tuple[int, Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]
                ] = []
                for pos_i, _, m_raw, c_raw, px_mm_i, is_up_i, ex_ln_i in recs_w:
                    m_f = _resolve_mask_file(m_raw, eff_data_dir)
                    m_u8 = cv2.imread(str(m_f), cv2.IMREAD_GRAYSCALE) if m_f is not None else None
                    if m_u8 is None or int(np.count_nonzero(m_u8)) == 0:
                        out_w.append((pos_i, None, None, None, None, None))
                        continue
                    mh_mw = m_u8.shape[:2]
                    c_f = _resolve_crop_file(c_raw, eff_data_dir)
                    c_bgr = cv2.imread(str(c_f), cv2.IMREAD_COLOR) if c_f is not None else None
                    res_p: Optional[Dict[str, float]] = None
                    try:
                        if c_bgr is not None and c_bgr.ndim == 3 and c_bgr.size > 0:
                            m_use = (
                                m_u8
                                if c_bgr.shape[:2] == mh_mw
                                else cv2.resize(m_u8, (c_bgr.shape[1], c_bgr.shape[0]), interpolation=cv2.INTER_NEAREST)
                            )
                            res_p = compute_calibrated_blade_photometry(
                                c_bgr, m_use, mask_origin=(0, 0), paper_l_ref=paper_ref, pixels_per_mm=px_mm_i
                            )
                        elif img_bgr is not None and sh_hw is not None:
                            if mh_mw == sh_hw:
                                res_p = compute_calibrated_blade_photometry(
                                    img_bgr, m_u8, mask_origin=(0, 0), paper_l_ref=paper_ref, pixels_per_mm=px_mm_i
                                )
                            else:
                                ar_m = float(mh_mw[1]) / float(max(mh_mw[0], 1))
                                ar_s = float(sh_hw[1]) / float(max(sh_hw[0], 1))
                                is_scaled_sheet_mask = (not is_up_i) or (
                                    abs(ar_m - ar_s) <= 0.08 * max(ar_s, 1e-6)
                                    and min(mh_mw) >= 0.25 * min(sh_hw)
                                )
                                if is_scaled_sheet_mask:
                                    m_res = cv2.resize(m_u8, (sh_hw[1], sh_hw[0]), interpolation=cv2.INTER_NEAREST)
                                    res_p = compute_calibrated_blade_photometry(
                                        img_bgr, m_res, mask_origin=(0, 0), paper_l_ref=paper_ref, pixels_per_mm=px_mm_i
                                    )
                                elif ex_ln_i is not None and 0.0 <= ex_ln_i <= 1.25 and paper_ref is not None:
                                    y_p = ((paper_ref + 16.0) / 116.0) ** 3 if paper_ref > 8.0 else paper_ref / 903.3
                                    y_rat = ((ex_ln_i + 0.16) / 1.16) ** 3 if ex_ln_i > 0.08 else (ex_ln_i / 9.03296296)
                                    y_bld = max(0.0, y_rat * max(y_p, 1e-4))
                                    med_l_rec = (
                                        116.0 * (y_bld ** (1.0 / 3.0)) - 16.0
                                        if y_bld > (6.0 / 29.0) ** 3
                                        else 903.3 * y_bld
                                    )
                                    lsd_rec = round(float(med_l_rec - paper_ref), 4)
                                    out_w.append((pos_i, None, None, lsd_rec if math.isfinite(lsd_rec) else None, None, None))
                                    continue
                    except Exception:
                        res_p = None
                    if not isinstance(res_p, dict):
                        out_w.append((pos_i, None, None, None, None, None))
                        continue

                    def _fin_opt(k_name: str) -> Optional[float]:
                        raw_v = res_p.get(k_name)
                        if raw_v is None:
                            return None
                        fv = float(raw_v)
                        return fv if math.isfinite(fv) else None

                    out_w.append((
                        pos_i,
                        _fin_opt("blade_cie_l_norm"),
                        _fin_opt("blade_cie_l_sd"),
                        _fin_opt("blade_cie_l_sheet_delta"),
                        _fin_opt("blade_cie_a_delta"),
                        _fin_opt("blade_cie_b_delta"),
                    ))
                return out_w

            photo_items = list(by_cat_photo.items())
            if len(photo_items) >= 4 and n_workers > 1:
                with ThreadPoolExecutor(max_workers=n_workers) as tpool:
                    photo_batches = list(tpool.map(_photo_sheet_worker, photo_items))
            else:
                photo_batches = [_photo_sheet_worker(it) for it in photo_items]
            for batch in photo_batches:
                for pos_i, ln_v, ls_v, lsd_v, a_v, b_v in batch:
                    did_upd = False
                    if ln_v is not None:
                        upd_cie_l_norm[pos_i] = ln_v
                        el_df.iat[pos_i, el_df.columns.get_loc("blade_cie_l_norm")] = ln_v
                        did_upd = True
                    if ls_v is not None:
                        upd_cie_l_sd[pos_i] = ls_v
                        el_df.iat[pos_i, el_df.columns.get_loc("blade_cie_l_sd")] = ls_v
                        did_upd = True
                    if lsd_v is not None:
                        upd_cie_l_sheet_delta[pos_i] = lsd_v
                        el_df.iat[pos_i, el_df.columns.get_loc("blade_cie_l_sheet_delta")] = lsd_v
                        did_upd = True
                    if a_v is not None:
                        upd_cie_a_delta[pos_i] = a_v
                        el_df.iat[pos_i, el_df.columns.get_loc("blade_cie_a_delta")] = a_v
                        did_upd = True
                    if b_v is not None:
                        upd_cie_b_delta[pos_i] = b_v
                        el_df.iat[pos_i, el_df.columns.get_loc("blade_cie_b_delta")] = b_v
                        did_upd = True
                    if did_upd:
                        if any(v is not None for v in (ln_v, ls_v, a_v, b_v)):
                            n_photo_upd += 1
                        else:
                            n_sheet_delta_rec += 1
            logger.info(
                f"CPU photometry recompute updated {n_photo_upd} / {n_rows} leaves where crop or sheet-matched mask existed "
                f"(+ {n_sheet_delta_rec} algebraic sheet-delta backfills)."
            )

        # Persist updated QC, ruler scale, and leaf_zone columns atomically to SQLite
        def _opt_float(val: Any) -> Optional[float]:
            if val is None or pd.isna(val):
                return None
            fv = float(val)
            return fv if math.isfinite(fv) else None

        sql_rows: List[Tuple[Any, ...]] = []
        for pos_i, (_, r_u) in enumerate(el_df.iterrows()):
            cjson_new = upd_contour_json[pos_i]
            sql_rows.append((
                _opt_float(r_u["solidity"]),
                _opt_float(r_u["circularity"]),
                int(r_u["detector_class_id"]),
                int(r_u["is_crop_clipped"]),
                float(r_u["hole_area_frac"]),
                float(r_u["enclosed_void_frac"]),
                int(r_u["bilateral_lobe_symmetry"]),
                int(r_u["efa_qc_pass"]),
                str(r_u["qc_flag_reason"]),
                str(r_u["leaf_zone"] or "basal"),
                _opt_float(r_u.get("leaf_zone_confidence")),
                _opt_float(r_u.get("dist_from_crown_norm")),
                _opt_float(r_u.get("pixels_per_mm")),
                int(r_u.get("is_calibrated_mm") or 0),
                str(r_u.get("scale_source") or ""),
                _opt_float(r_u.get("area_mm2")),
                _opt_float(r_u.get("blade_length_mm")),
                _opt_float(r_u.get("blade_width_mm")),
                _opt_float(r_u.get("petiole_length_mm")),
                _opt_float(r_u.get("petiole_width_mm")),
                cjson_new,
                _opt_float(upd_cie_l_norm[pos_i]),
                _opt_float(upd_cie_l_sd[pos_i]),
                _opt_float(upd_cie_l_sheet_delta[pos_i]),
                _opt_float(upd_cie_a_delta[pos_i]),
                _opt_float(upd_cie_b_delta[pos_i]),
                str(r_u["catalogNumber"]),
                int(r_u["plant_individual_id"]),
                int(r_u["leaf_id"]),
            ))

        with conn:
            conn.executemany(
                """
                UPDATE extracted_leaves
                SET solidity = COALESCE(?, solidity),
                    circularity = COALESCE(?, circularity),
                    detector_class_id = ?,
                    is_crop_clipped = ?,
                    hole_area_frac = ?,
                    enclosed_void_frac = ?,
                    bilateral_lobe_symmetry = ?,
                    efa_qc_pass = ?,
                    qc_flag_reason = ?,
                    leaf_zone = ?,
                    leaf_zone_confidence = COALESCE(?, leaf_zone_confidence),
                    dist_from_crown_norm = COALESCE(?, dist_from_crown_norm),
                    pixels_per_mm = COALESCE(?, pixels_per_mm),
                    is_calibrated_mm = ?,
                    scale_source = CASE WHEN ? != '' THEN ? ELSE scale_source END,
                    area_mm2 = COALESCE(?, area_mm2),
                    blade_length_mm = COALESCE(?, blade_length_mm),
                    blade_width_mm = COALESCE(?, blade_width_mm),
                    petiole_length_mm = COALESCE(?, petiole_length_mm),
                    petiole_width_mm = COALESCE(?, petiole_width_mm),
                    contour_xy_json = COALESCE(?, contour_xy_json),
                    blade_cie_l_norm = COALESCE(?, blade_cie_l_norm),
                    blade_cie_l_sd = COALESCE(?, blade_cie_l_sd),
                    blade_cie_l_sheet_delta = COALESCE(?, blade_cie_l_sheet_delta),
                    blade_cie_a_delta = COALESCE(?, blade_cie_a_delta),
                    blade_cie_b_delta = COALESCE(?, blade_cie_b_delta)
                WHERE catalogNumber = ? AND plant_individual_id = ? AND leaf_id = ?;
                """,
                [
                    r[:14] + (r[14], r[14]) + r[15:]
                    for r in sql_rows
                ],
            )
            if crown_backfills:
                conn.executemany(
                    """
                    UPDATE plant_individuals
                    SET crown_x = COALESCE(crown_x, ?),
                        crown_y = COALESCE(crown_y, ?)
                    WHERE catalogNumber = ? AND plant_individual_id = ?;
                    """,
                    crown_backfills,
                )
            if ruler_cal_scales:
                v_cal_rows = [(float(px), str(cat)) for cat, px in ruler_cal_scales.items()]
                conn.executemany(
                    """
                    UPDATE vouchers
                    SET pixels_per_mm = ?, is_calibrated_mm = 1
                    WHERE catalogNumber = ?;
                    """,
                    v_cal_rows,
                )
            if unboxed_scales:
                v_rows = [
                    (float(px), str(cat))
                    for cat, px in unboxed_scales.items()
                    if cat not in ruler_cal_scales
                ]
                if v_rows:
                    conn.executemany(
                        """
                        UPDATE vouchers
                        SET pixels_per_mm = ?, is_calibrated_mm = 0
                        WHERE catalogNumber = ?;
                        """,
                        v_rows,
                    )

        # Synchronize companion CSV snapshots when present
        tables_dir = eff_manifest.parent
        if ruler_cal_scales or unboxed_scales:
            curated_csv = tables_dir / "curated_vouchers.csv"
            if curated_csv.is_file() and curated_csv.stat().st_size > 0:
                c_df = pd.read_csv(curated_csv, low_memory=False)
                if "catalogNumber" in c_df.columns:
                    for ci in c_df.index:
                        c_key = str(c_df.at[ci, "catalogNumber"]).strip()
                        if c_key in ruler_cal_scales:
                            px_new = float(ruler_cal_scales[c_key])
                            if "pixels_per_mm" in c_df.columns:
                                c_df.at[ci, "pixels_per_mm"] = px_new
                            if "is_calibrated_mm" in c_df.columns:
                                c_df.at[ci, "is_calibrated_mm"] = True
                            if "scale_source" in c_df.columns:
                                c_df.at[ci, "scale_source"] = "ruler_calibrated"
                        elif c_key in unboxed_scales:
                            px_new = float(unboxed_scales[c_key])
                            if "pixels_per_mm" in c_df.columns:
                                c_df.at[ci, "pixels_per_mm"] = px_new
                            if "is_calibrated_mm" in c_df.columns:
                                c_df.at[ci, "is_calibrated_mm"] = False
                            if "scale_source" in c_df.columns:
                                c_df.at[ci, "scale_source"] = "tick_lattice_unboxed"
                    atomic_write_text(curated_csv, c_df.to_csv(index=False))

        if eff_manifest.is_file() and eff_manifest.stat().st_size > 0:
            m_df = pd.read_csv(eff_manifest, low_memory=False)
            if {"catalogNumber", "leaf_id"}.issubset(m_df.columns):
                if "plant_individual_id" not in m_df.columns:
                    m_df["plant_individual_id"] = 0
                key_to_idx: Dict[Tuple[str, int, int], int] = {}
                for r_i, mr in m_df.iterrows():
                    k_m = (
                        str(mr["catalogNumber"]).strip(),
                        int(mr["plant_individual_id"]) if pd.notna(mr["plant_individual_id"]) else 0,
                        int(mr["leaf_id"]),
                    )
                    key_to_idx[k_m] = int(r_i)
                sync_cols = (
                    "solidity",
                    "circularity",
                    "detector_class_id",
                    "is_crop_clipped",
                    "hole_area_frac",
                    "enclosed_void_frac",
                    "bilateral_lobe_symmetry",
                    "efa_qc_pass",
                    "qc_flag_reason",
                    "leaf_zone",
                    "leaf_zone_confidence",
                    "dist_from_crown_norm",
                    "pixels_per_mm",
                    "is_calibrated_mm",
                    "scale_source",
                    "area_mm2",
                    "blade_length_mm",
                    "blade_width_mm",
                    "petiole_length_mm",
                    "petiole_width_mm",
                    "blade_cie_l_norm",
                    "blade_cie_l_sd",
                    "blade_cie_l_sheet_delta",
                    "blade_cie_a_delta",
                    "blade_cie_b_delta",
                )
                for col_s in sync_cols:
                    if col_s not in m_df.columns:
                        if col_s in {"qc_flag_reason", "leaf_zone", "scale_source"}:
                            m_df[col_s] = ""
                        elif col_s in {"is_crop_clipped", "bilateral_lobe_symmetry", "efa_qc_pass", "is_calibrated_mm"}:
                            m_df[col_s] = col_s == "efa_qc_pass"
                        elif col_s == "detector_class_id":
                            m_df[col_s] = 0
                        else:
                            m_df[col_s] = np.nan
                    elif col_s not in {
                        "qc_flag_reason", "leaf_zone", "scale_source",
                        "is_crop_clipped", "bilateral_lobe_symmetry", "efa_qc_pass",
                        "is_calibrated_mm", "detector_class_id",
                    }:
                        m_df[col_s] = pd.to_numeric(m_df[col_s], errors="coerce")
                for _, r_u in el_df.iterrows():
                    k_u = (str(r_u["catalogNumber"]).strip(), int(r_u["plant_individual_id"]), int(r_u["leaf_id"]))
                    if k_u in key_to_idx:
                        mi = key_to_idx[k_u]
                        for col_s in sync_cols:
                            if col_s in {"is_crop_clipped", "bilateral_lobe_symmetry", "efa_qc_pass", "is_calibrated_mm"}:
                                m_df.at[mi, col_s] = bool(r_u[col_s])
                            elif col_s == "detector_class_id":
                                m_df.at[mi, col_s] = int(r_u[col_s])
                            elif col_s in {"qc_flag_reason", "leaf_zone", "scale_source"}:
                                m_df.at[mi, col_s] = str(r_u.get(col_s) or "")
                            else:
                                if col_s in r_u:
                                    m_df.at[mi, col_s] = r_u[col_s]
                atomic_write_text(eff_manifest, m_df.to_csv(index=False))
                logger.info(f"Synchronized updated QC and zone columns to {eff_manifest}")

        n_pass = int((el_df["efa_qc_pass"] == 1).sum())
        n_flag = int(n_rows - n_pass)
        n_bas = int((el_df["leaf_zone"].astype(str).str.lower() == "basal").sum())
        n_caul = int((el_df["leaf_zone"].astype(str).str.lower() == "cauline").sum())
        logger.info(
            f"Re-scored {n_rows} leaves: efa_qc_pass={n_pass}, flagged={n_flag}, "
            f"basal={n_bas}, cauline={n_caul}, ruler_calibrated={len(ruler_cal_scales)}, "
            f"tick_lattice_unboxed={len(unboxed_scales)}."
        )
        return {
            "total": n_rows,
            "efa_qc_pass": n_pass,
            "flagged": n_flag,
            "basal": n_bas,
            "cauline": n_caul,
            "ruler_calibrated_vouchers": len(ruler_cal_scales),
            "tick_lattice_unboxed_vouchers": len(unboxed_scales),
            "photometry_recomputed": n_photo_upd,
            "sheet_delta_recovered": n_sheet_delta_rec,
        }
    finally:
        conn.close()


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Re-score Step 02 leaf QC gates and anatomical basal/cauline zones in SQLite + CSV without GPU inference."
    )
    parser.add_argument("--db", type=Path, default=None, help="Optional SQLite database path (default: data/plantmatics.sqlite).")
    parser.add_argument("--data-dir", type=Path, default=PROJECT_ROOT / "data", help="Tier A data directory containing masks/, contours/, and tables/.")
    parser.add_argument("--manifest", type=Path, default=None, help="Optional extracted_leaf_manifest.csv path override.")
    parser.add_argument("--allow-partial-leaves", action="store_true", default=False, help="Allow detector_class_id == 1 (Leaf_PARTIAL) to pass QC.")
    parser.add_argument("--recompute-photometry", action="store_true", default=False, help="Recompute blade_cie_l_norm/sd where crops or sheet-matched masks exist.")
    parser.add_argument("--workers", type=int, default=None, help="Worker thread/process count cap (default: auto up to 32).")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    rescore_tier_a_leaves(
        db_path=args.db,
        data_dir=args.data_dir,
        manifest_csv=args.manifest,
        segment_partial_leaves=True if args.allow_partial_leaves else None,
        recompute_photometry=bool(args.recompute_photometry),
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
