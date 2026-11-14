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
import json
import logging
import math
import os
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
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
from scripts.vision.leaf_geometry_utils import (
    compute_lm3_shape_metrics,
    export_standardized_contour,
    open_hairline_sinus_bridges_and_measure_voids,
)
from scripts.vision.lm3_inference_engines import (
    evaluate_leaf_damage_and_fold_qc,
    solve_ruler_tick_lattice_px_per_mm,
)
from scripts.vision.organ_crown_metrics import classify_leaf_zone_anatomical

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


def _atomic_write_csv(df: pd.DataFrame, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(target.parent), prefix=f".{target.stem}_", suffix=".csv.tmp")
    tmp_p = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as fh:
            df.to_csv(fh, index=False)
            fh.flush()
            os.fsync(fh.fileno())
        tmp_p.replace(target)
    except Exception:
        if tmp_p.exists():
            tmp_p.unlink(missing_ok=True)
        raise


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

        _, evf_val, n_opened = open_hairline_sinus_bridges_and_measure_voids(mask_u8)
        shape_m = compute_lm3_shape_metrics(mask_u8)
        if shape_m["area_px"] > 0.0:
            sol_val = round(float(shape_m["solidity"]), 4)
            circ_val = round(float(shape_m["circularity"]), 4)
        if n_opened > 0:
            c_res = export_standardized_contour(
                mask_u8,
                str(row["catalogNumber"]),
                int(row["leaf_id"]),
                eff_data_dir,
                return_coords=True,
                already_oriented=bool(_ival("mask_is_upright", 1) or 1),
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
    """Run CPU margin metric tick-lattice calibration for one uncalibrated voucher sheet."""
    cat, img_path_str, mp_prior_mm = task
    img = cv2.imread(img_path_str, cv2.IMREAD_COLOR)
    if img is None or img.size == 0:
        return None
    px_per_mm = solve_ruler_tick_lattice_px_per_mm(img, mp_prior_px_per_mm=mp_prior_mm)
    if px_per_mm is None or not np.isfinite(float(px_per_mm)) or float(px_per_mm) <= 1.0:
        return None
    return (cat, round(float(px_per_mm), 4))


def _is_inverted_mount_group(grp: pd.DataFrame) -> bool:
    """Detect upside-down mounted shoots where large petiolate leaves sit near top (cy_rel < 0.35)."""
    cy_v = pd.to_numeric(grp.get("cy_rel_plant_bbox"), errors="coerce")
    area_v = pd.to_numeric(grp.get("area_px"), errors="coerce")
    pbr_v = pd.to_numeric(grp.get("petiole_to_blade_ratio"), errors="coerce").fillna(0.0)
    valid = cy_v.notna() & area_v.notna() & (area_v > 0)
    if int(valid.sum()) < 3:
        return False
    cy_f = cy_v[valid].to_numpy(dtype=float)
    ar_f = area_v[valid].to_numpy(dtype=float)
    pb_f = pbr_v[valid].to_numpy(dtype=float)
    max_area = float(np.max(ar_f))
    top_large_pet = np.any((cy_f <= 0.35) & (ar_f >= 0.70 * max_area) & (pb_f >= 0.22))
    bot_small_sess = np.any((cy_f >= 0.55) & ((ar_f <= 0.55 * max_area) | (pb_f <= 0.15)))
    if top_large_pet and bot_small_sess:
        return True
    med_area = float(np.median(ar_f))
    top_half_cy = cy_f[ar_f >= med_area]
    bot_half_cy = cy_f[ar_f < med_area]
    if top_half_cy.size >= 1 and bot_half_cy.size >= 1:
        return bool(float(np.median(top_half_cy)) < 0.38 and float(np.median(bot_half_cy)) > 0.62)
    return False


def rescore_tier_a_leaves(
    db_path: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    manifest_csv: Optional[Path] = None,
    segment_partial_leaves: Optional[bool] = None,
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
            return {"total": 0, "efa_qc_pass": 0, "flagged": 0, "basal": 0, "cauline": 0}

        n_rows = len(el_df)
        qc_cols = [
            "catalogNumber", "leaf_id", "mask_path", "solidity", "circularity",
            "bilateral_dice", "si_a", "lm3_gates_pass", "n_lamina_components",
            "detector_class_id", "hole_area_frac", "is_crop_clipped", "mask_is_upright",
        ]
        avail_qc_cols = [c for c in qc_cols if c in el_df.columns]
        qc_tasks = [
            (rec, str(eff_data_dir), seg_partial)
            for rec in el_df[avail_qc_cols].to_dict(orient="records")
        ]
        n_workers = min(32, max(1, os.cpu_count() or 8))
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

        # CPU metric margin tick-lattice calibration for uncalibrated vouchers
        ruler_tasks: List[Tuple[str, str, Optional[float]]] = []
        seen_cats: set[str] = set()
        for _, rv in el_df[["catalogNumber", "image_path", "voucher_px_per_mm", "voucher_is_cal"]].drop_duplicates("catalogNumber").iterrows():
            cat_s = str(rv["catalogNumber"]).strip()
            if not cat_s or cat_s in seen_cats:
                continue
            seen_cats.add(cat_s)
            is_cal_v = int(rv["voucher_is_cal"]) if pd.notna(rv.get("voucher_is_cal")) else 0
            if is_cal_v != 0:
                continue
            img_f = _resolve_voucher_image_file(rv.get("image_path"), cat_s, eff_data_dir)
            if img_f is None:
                continue
            v_px = float(rv["voucher_px_per_mm"]) if pd.notna(rv.get("voucher_px_per_mm")) else None
            ruler_tasks.append((cat_s, str(img_f), v_px))

        calibrated_scales: Dict[str, float] = {}
        if ruler_tasks:
            if len(ruler_tasks) >= 8 and n_workers > 1:
                with ProcessPoolExecutor(max_workers=n_workers) as pool:
                    hits = list(pool.map(_calibrate_voucher_ruler_worker, ruler_tasks, chunksize=8))
            else:
                hits = [_calibrate_voucher_ruler_worker(t) for t in ruler_tasks]
            for h in hits:
                if h is not None:
                    calibrated_scales[h[0]] = h[1]
            logger.info(
                f"CPU margin ruler calibration: upgraded {len(calibrated_scales)} / {len(ruler_tasks)} uncalibrated vouchers."
            )

        if "scale_source" not in el_df.columns:
            el_df["scale_source"] = ""
        if "is_calibrated_mm" not in el_df.columns:
            el_df["is_calibrated_mm"] = 0
        if "pixels_per_mm" not in el_df.columns:
            el_df["pixels_per_mm"] = np.nan
        if calibrated_scales:
            cat_series = el_df["catalogNumber"].astype(str).str.strip()
            cal_mask = cat_series.isin(calibrated_scales)
            for idx_c in el_df.index[cal_mask]:
                cat_c = str(el_df.at[idx_c, "catalogNumber"]).strip()
                new_px = float(calibrated_scales[cat_c])
                el_df.at[idx_c, "pixels_per_mm"] = new_px
                el_df.at[idx_c, "is_calibrated_mm"] = 1
                el_df.at[idx_c, "scale_source"] = "ruler_calibrated"
                if pd.notna(el_df.at[idx_c, "area_px"]):
                    el_df.at[idx_c, "area_mm2"] = round(float(el_df.at[idx_c, "area_px"]) / (new_px * new_px), 4)
                for px_col, mm_col in (
                    ("blade_length_px", "blade_length_mm"),
                    ("blade_width_px", "blade_width_mm"),
                    ("petiole_length_px", "petiole_length_mm"),
                    ("petiole_width_px", "petiole_width_mm"),
                ):
                    if px_col in el_df.columns and mm_col in el_df.columns and pd.notna(el_df.at[idx_c, px_col]):
                        el_df.at[idx_c, mm_col] = round(float(el_df.at[idx_c, px_col]) / new_px, 4)

        # Normalize scale_source across calibrated vs mp_estimated rows
        is_cal_arr = pd.to_numeric(el_df["is_calibrated_mm"], errors="coerce").fillna(0).astype(int)
        px_arr = pd.to_numeric(el_df["pixels_per_mm"], errors="coerce")
        ss_cur = el_df["scale_source"].fillna("").astype(str).str.strip()
        el_df.loc[is_cal_arr == 1, "scale_source"] = "ruler_calibrated"
        mp_fix = (is_cal_arr == 0) & px_arr.notna() & (px_arr > 0.0) & (~ss_cur.isin(["ruler_calibrated", "mp_estimated"]))
        el_df.loc[mp_fix, "scale_source"] = "mp_estimated"

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
                inverted_mount = _is_inverted_mount_group(grp)
                leaf_recs: List[Dict[str, Any]] = []
                for idx_df, r_l in grp.iterrows():
                    lid_i = int(r_l["leaf_id"])
                    bld_px = float(r_l["blade_length_px"]) if pd.notna(r_l.get("blade_length_px")) and float(r_l["blade_length_px"]) > 0 else 80.0
                    bld_w_px = float(r_l["blade_width_px"]) if pd.notna(r_l.get("blade_width_px")) and float(r_l["blade_width_px"]) > 0 else np.nan
                    pet_px = float(r_l["petiole_length_px"]) if pd.notna(r_l.get("petiole_length_px")) and float(r_l["petiole_length_px"]) >= 0 else 0.0
                    pet_w_px = float(r_l["petiole_width_px"]) if pd.notna(r_l.get("petiole_width_px")) and float(r_l["petiole_width_px"]) >= 0 else np.nan
                    pbr_v = float(r_l["petiole_to_blade_ratio"]) if pd.notna(r_l.get("petiole_to_blade_ratio")) else (pet_px / max(bld_px, 1.0))
                    dc_prev = float(r_l["dist_from_crown_norm"]) if pd.notna(r_l.get("dist_from_crown_norm")) else np.nan
                    cy_rel = float(r_l["cy_rel_plant_bbox"]) if pd.notna(r_l.get("cy_rel_plant_bbox")) else np.nan
                    if np.isfinite(cy_rel):
                        s_est = max(0.0, min(1.0, cy_rel if inverted_mount else (1.0 - cy_rel)))
                        eff_cy_rel = 1.0 - s_est
                    elif np.isfinite(dc_prev):
                        s_est = max(0.0, min(1.0, dc_prev))
                        eff_cy_rel = 1.0 - s_est
                    else:
                        s_est = 0.0
                        eff_cy_rel = 1.0
                    y_tip = cy_eval - s_est * p_height
                    y_base = y_tip - min(pet_px, 0.35 * p_height)
                    area_val = float(r_l["area_px"]) if pd.notna(r_l.get("area_px")) and float(r_l["area_px"]) > 0 else 0.0
                    ar_val = float(r_l["aspect_ratio"]) if pd.notna(r_l.get("aspect_ratio")) and float(r_l["aspect_ratio"]) > 0 else 1.5
                    sol_l = float(r_l["solidity"]) if pd.notna(r_l.get("solidity")) else 0.85
                    leaf_recs.append({
                        "leaf_id": lid_i,
                        "plant_individual_id": int(_pid),
                        "p_base_sheet": (cx_eval, y_base),
                        "p_apex_sheet": (cx_eval, y_base - bld_px),
                        "p_petiole_tip": (cx_eval, y_tip),
                        "blade_length_px": bld_px,
                        "blade_width_px": bld_w_px,
                        "petiole_length_px": pet_px,
                        "petiole_width_px": pet_w_px,
                        "petiole_to_blade_ratio": pbr_v,
                        "has_petiole_mask": bool(pbr_v >= 0.12),
                        "area_px": area_val,
                        "aspect_ratio": ar_val,
                        "solidity": sol_l,
                        "cy_rel_plant_bbox": eff_cy_rel,
                    })
                cap_box = (
                    [(cy_eval - p_height - 20.0, cx_eval - 20.0, cy_eval - p_height + 20.0, cx_eval + 20.0)]
                    if has_caps
                    else []
                )
                z_res = classify_leaf_zone_anatomical(
                    leaf_recs,
                    crowns_by_pid={int(_pid): (cx_eval, cy_eval)},
                    capitula_detections=cap_box,
                    is_bent_stem=bool(int(r0.get("is_bent_stem") or 0)),
                )
                for _, r_l in grp.iterrows():
                    lid_i = int(r_l["leaf_id"])
                    if lid_i in z_res:
                        el_df.loc[(el_df["catalogNumber"] == _cat) & (el_df["plant_individual_id"] == _pid) & (el_df["leaf_id"] == lid_i), "leaf_zone"] = z_res[lid_i]["leaf_zone"]
                        el_df.loc[(el_df["catalogNumber"] == _cat) & (el_df["plant_individual_id"] == _pid) & (el_df["leaf_id"] == lid_i), "leaf_zone_confidence"] = z_res[lid_i]["leaf_zone_confidence"]
                        el_df.loc[(el_df["catalogNumber"] == _cat) & (el_df["plant_individual_id"] == _pid) & (el_df["leaf_id"] == lid_i), "dist_from_crown_norm"] = z_res[lid_i]["dist_from_crown_norm"]

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
                    contour_xy_json = COALESCE(?, contour_xy_json)
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
            if calibrated_scales:
                v_rows = [(float(px), str(cat)) for cat, px in calibrated_scales.items()]
                conn.executemany(
                    """
                    UPDATE vouchers
                    SET pixels_per_mm = ?, is_calibrated_mm = 1
                    WHERE catalogNumber = ?;
                    """,
                    v_rows,
                )
                conn.executemany(
                    """
                    UPDATE plant_individuals
                    SET plant_height_mm = CASE WHEN plant_height_px IS NOT NULL THEN ROUND(plant_height_px / ?, 4) ELSE plant_height_mm END,
                        plant_height_cm = CASE WHEN plant_height_px IS NOT NULL THEN ROUND(plant_height_px / (? * 10.0), 4) ELSE plant_height_cm END,
                        petiole_length_mm = CASE WHEN petiole_length_px IS NOT NULL THEN ROUND(petiole_length_px / ?, 4) ELSE petiole_length_mm END,
                        petiole_width_mm = CASE WHEN petiole_width_px IS NOT NULL THEN ROUND(petiole_width_px / ?, 4) ELSE petiole_width_mm END,
                        blade_length_mm = CASE WHEN blade_length_px IS NOT NULL THEN ROUND(blade_length_px / ?, 4) ELSE blade_length_mm END
                    WHERE catalogNumber = ?;
                    """,
                    [(float(px), float(px), float(px), float(px), float(px), str(cat)) for cat, px in calibrated_scales.items()],
                )
                conn.executemany(
                    """
                    UPDATE voucher_reproductive_phenology
                    SET involucre_height_mm = CASE WHEN involucre_height_px IS NOT NULL THEN ROUND(involucre_height_px / ?, 4) ELSE involucre_height_mm END,
                        involucre_width_mm = CASE WHEN involucre_width_px IS NOT NULL THEN ROUND(involucre_width_px / ?, 4) ELSE involucre_width_mm END
                    WHERE catalogNumber = ?;
                    """,
                    [(float(px), float(px), str(cat)) for cat, px in calibrated_scales.items()],
                )

        # Synchronize companion CSV snapshots when present
        tables_dir = eff_manifest.parent
        if calibrated_scales:
            curated_csv = tables_dir / "curated_vouchers.csv"
            if curated_csv.is_file() and curated_csv.stat().st_size > 0:
                c_df = pd.read_csv(curated_csv, low_memory=False)
                if "catalogNumber" in c_df.columns:
                    for ci in c_df.index:
                        c_key = str(c_df.at[ci, "catalogNumber"]).strip()
                        if c_key in calibrated_scales:
                            px_new = float(calibrated_scales[c_key])
                            if "pixels_per_mm" in c_df.columns:
                                c_df.at[ci, "pixels_per_mm"] = px_new
                            if "is_calibrated_mm" in c_df.columns:
                                c_df.at[ci, "is_calibrated_mm"] = True
                            if "scale_source" in c_df.columns:
                                c_df.at[ci, "scale_source"] = "ruler_calibrated"
                    _atomic_write_csv(c_df, curated_csv)
            repro_csv = tables_dir / "voucher_reproductive_metrics.csv"
            if repro_csv.is_file() and repro_csv.stat().st_size > 0:
                vr_df = pd.read_csv(repro_csv, low_memory=False)
                if "catalogNumber" in vr_df.columns:
                    for vi in vr_df.index:
                        c_key = str(vr_df.at[vi, "catalogNumber"]).strip()
                        if c_key in calibrated_scales:
                            px_new = float(calibrated_scales[c_key])
                            for px_col, mm_col in (
                                ("involucre_height_px", "involucre_height_mm"),
                                ("involucre_width_px", "involucre_width_mm"),
                            ):
                                if px_col in vr_df.columns and mm_col in vr_df.columns and pd.notna(vr_df.at[vi, px_col]):
                                    vr_df.at[vi, mm_col] = round(float(vr_df.at[vi, px_col]) / px_new, 4)
                    _atomic_write_csv(vr_df, repro_csv)

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
                _atomic_write_csv(m_df, eff_manifest)
                logger.info(f"Synchronized updated QC and zone columns to {eff_manifest}")

        n_pass = int((el_df["efa_qc_pass"] == 1).sum())
        n_flag = int(n_rows - n_pass)
        n_bas = int((el_df["leaf_zone"].astype(str).str.lower() == "basal").sum())
        n_caul = int((el_df["leaf_zone"].astype(str).str.lower() == "cauline").sum())
        logger.info(
            f"Re-scored {n_rows} leaves: efa_qc_pass={n_pass}, flagged={n_flag}, "
            f"basal={n_bas}, cauline={n_caul}, ruler_upgraded={len(calibrated_scales)}."
        )
        return {
            "total": n_rows,
            "efa_qc_pass": n_pass,
            "flagged": n_flag,
            "basal": n_bas,
            "cauline": n_caul,
            "ruler_calibrated_vouchers": len(calibrated_scales),
        }
    finally:
        conn.close()


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Re-score Step 02 leaf QC gates and anatomical basal/cauline zones in SQLite + CSV without GPU inference."
    )
    parser.add_argument("--db", type=Path, default=None, help="Optional SQLite database path (default: data/packera_pipeline.sqlite).")
    parser.add_argument("--data-dir", type=Path, default=PROJECT_ROOT / "data", help="Tier A data directory containing masks/, contours/, and tables/.")
    parser.add_argument("--manifest", type=Path, default=None, help="Optional extracted_leaf_manifest.csv path override.")
    parser.add_argument("--allow-partial-leaves", action="store_true", default=False, help="Allow detector_class_id == 1 (Leaf_PARTIAL) to pass QC.")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    rescore_tier_a_leaves(
        db_path=args.db,
        data_dir=args.data_dir,
        manifest_csv=args.manifest,
        segment_partial_leaves=True if args.allow_partial_leaves else None,
    )


if __name__ == "__main__":
    main()
