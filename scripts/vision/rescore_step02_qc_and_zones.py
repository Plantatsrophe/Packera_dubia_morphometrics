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
from scripts.vision.lm3_inference_engines import evaluate_leaf_damage_and_fold_qc
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


def rescore_tier_a_leaves(
    db_path: Optional[Path] = None,
    data_dir: Optional[Path] = None,
    manifest_csv: Optional[Path] = None,
    segment_partial_leaves: Optional[bool] = None,
) -> Dict[str, int]:
    """Re-evaluates leaf QC gates, hairline sinus bridges, and basal/cauline zones in SQLite + CSV."""
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
                COALESCE(vrp.capitula_count, vrp0.capitula_count, 0) AS eff_capitula_count
            FROM extracted_leaves el
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
        upd_solidity: List[float] = []
        upd_circ: List[float] = []
        upd_evf: List[float] = []
        upd_haf: List[float] = []
        upd_clipped: List[int] = []
        upd_dcid: List[int] = []
        upd_bls: List[int] = []
        upd_eqc: List[int] = []
        upd_qfr: List[str] = []
        upd_contour_json: List[Optional[str]] = []
        min_rms_list: List[float] = []

        for _, row in el_df.iterrows():
            m_file = _resolve_mask_file(row.get("mask_path"), eff_data_dir)
            mask_u8: Optional[np.ndarray] = None
            if m_file is not None:
                m_raw = cv2.imread(str(m_file), cv2.IMREAD_GRAYSCALE)
                if m_raw is not None and np.count_nonzero(m_raw) > 0:
                    mask_u8 = ((m_raw > 0).astype(np.uint8)) * 255

            sol_val = float(row["solidity"]) if pd.notna(row.get("solidity")) else 0.85
            circ_val = float(row["circularity"]) if pd.notna(row.get("circularity")) else np.nan
            bd_val = float(row["bilateral_dice"]) if pd.notna(row.get("bilateral_dice")) else 0.85
            sia_val = float(row["si_a"]) if pd.notna(row.get("si_a")) else 0.0
            gp_val = bool(int(row["lm3_gates_pass"])) if pd.notna(row.get("lm3_gates_pass")) else True
            ncc_val = int(row["n_lamina_components"]) if pd.notna(row.get("n_lamina_components")) else 1
            dcid_val = int(row["detector_class_id"]) if ("detector_class_id" in row and pd.notna(row.get("detector_class_id"))) else 0
            haf_val = float(row["hole_area_frac"]) if ("hole_area_frac" in row and pd.notna(row.get("hole_area_frac"))) else 0.0
            clipped_val = bool(int(row["is_crop_clipped"])) if ("is_crop_clipped" in row and pd.notna(row.get("is_crop_clipped"))) else False
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
                        already_oriented=bool(int(row.get("mask_is_upright", 1) or 1)),
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

            upd_solidity.append(sol_val)
            upd_circ.append(circ_val)
            upd_evf.append(round(float(evf_val), 4))
            upd_haf.append(round(float(haf_val), 4))
            upd_clipped.append(1 if clipped_val else 0)
            upd_dcid.append(dcid_val)
            upd_bls.append(1 if qc_res["bilateral_lobe_symmetry"] else 0)
            upd_eqc.append(1 if qc_res["efa_qc_pass"] else 0)
            upd_qfr.append(str(qc_res["qc_flag_reason"]))
            upd_contour_json.append(new_cjson)
            min_rms_list.append(float(qc_res["min_margin_rms_norm"]))

        el_df["solidity"] = upd_solidity
        el_df["circularity"] = upd_circ
        el_df["enclosed_void_frac"] = upd_evf
        el_df["hole_area_frac"] = upd_haf
        el_df["is_crop_clipped"] = upd_clipped
        el_df["detector_class_id"] = upd_dcid
        el_df["bilateral_lobe_symmetry"] = upd_bls
        el_df["efa_qc_pass"] = upd_eqc
        el_df["qc_flag_reason"] = upd_qfr
        el_df["_min_rms"] = min_rms_list

        # Re-evaluate basal vs. cauline leaf zones per (catalogNumber, plant_individual_id)
        for (_cat, _pid), grp in el_df.groupby(["catalogNumber", "plant_individual_id"], sort=False):
            if len(grp) == 0:
                continue
            r0 = grp.iloc[0]
            cx_cr = float(r0["crown_x"]) if pd.notna(r0.get("crown_x")) else None
            cy_cr = float(r0["crown_y"]) if pd.notna(r0.get("crown_y")) else None
            has_caps = int(r0.get("eff_capitula_count") or 0) > 0
            p_height = float(r0["plant_height_px"]) if pd.notna(r0.get("plant_height_px")) and float(r0["plant_height_px"]) > 20.0 else None

            # Reconstruct vertical stem frame along crown -> corymb when stored leaf bboxes are absent
            if cx_cr is not None and cy_cr is not None and has_caps and p_height is not None:
                leaf_recs: List[Dict[str, Any]] = []
                for idx_df, r_l in grp.iterrows():
                    lid_i = int(r_l["leaf_id"])
                    bld_px = float(r_l["blade_length_px"]) if pd.notna(r_l.get("blade_length_px")) and float(r_l["blade_length_px"]) > 0 else 80.0
                    pet_px = float(r_l["petiole_length_px"]) if pd.notna(r_l.get("petiole_length_px")) and float(r_l["petiole_length_px"]) >= 0 else 0.0
                    pbr_v = float(r_l["petiole_to_blade_ratio"]) if pd.notna(r_l.get("petiole_to_blade_ratio")) else (pet_px / max(bld_px, 1.0))
                    dc_prev = float(r_l["dist_from_crown_norm"]) if pd.notna(r_l.get("dist_from_crown_norm")) else np.nan
                    cy_rel = float(r_l["cy_rel_plant_bbox"]) if pd.notna(r_l.get("cy_rel_plant_bbox")) else np.nan
                    if np.isfinite(cy_rel):
                        s_est = max(0.0, min(1.0, 1.0 - cy_rel))
                    elif np.isfinite(dc_prev):
                        s_est = max(0.0, min(1.0, dc_prev))
                    else:
                        s_est = 0.0
                    y_tip = cy_cr - s_est * p_height
                    y_base = y_tip - min(pet_px, 0.35 * p_height)
                    leaf_recs.append({
                        "leaf_id": lid_i,
                        "plant_individual_id": int(_pid),
                        "p_base_sheet": (cx_cr, y_base),
                        "p_apex_sheet": (cx_cr, y_base - bld_px),
                        "p_petiole_tip": (cx_cr, y_tip),
                        "blade_length_px": bld_px,
                        "petiole_length_px": pet_px,
                        "petiole_to_blade_ratio": pbr_v,
                        "has_petiole_mask": bool(pbr_v >= 0.12),
                    })
                cap_box = [(cy_cr - p_height - 20.0, cx_cr - 20.0, cy_cr - p_height + 20.0, cx_cr + 20.0)]
                z_res = classify_leaf_zone_anatomical(
                    leaf_recs,
                    crowns_by_pid={int(_pid): (cx_cr, cy_cr)},
                    capitula_detections=cap_box,
                    is_bent_stem=bool(int(r0.get("is_bent_stem") or 0)),
                )
                for idx_df, r_l in grp.iterrows():
                    lid_i = int(r_l["leaf_id"])
                    if lid_i in z_res:
                        el_df.at[idx_df, "leaf_zone"] = z_res[lid_i]["leaf_zone"]
                        el_df.at[idx_df, "leaf_zone_confidence"] = z_res[lid_i]["leaf_zone_confidence"]
                        el_df.at[idx_df, "dist_from_crown_norm"] = z_res[lid_i]["dist_from_crown_norm"]

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

        # Persist updated QC and leaf_zone columns atomically to SQLite extracted_leaves
        sql_rows: List[Tuple[Any, ...]] = []
        for pos_i, (_, r_u) in enumerate(el_df.iterrows()):
            cjson_new = upd_contour_json[pos_i]
            sql_rows.append((
                float(r_u["solidity"]) if pd.notna(r_u["solidity"]) else None,
                float(r_u["circularity"]) if pd.notna(r_u["circularity"]) else None,
                int(r_u["detector_class_id"]),
                int(r_u["is_crop_clipped"]),
                float(r_u["hole_area_frac"]),
                float(r_u["enclosed_void_frac"]),
                int(r_u["bilateral_lobe_symmetry"]),
                int(r_u["efa_qc_pass"]),
                str(r_u["qc_flag_reason"]),
                str(r_u["leaf_zone"] or "basal"),
                float(r_u["leaf_zone_confidence"]) if pd.notna(r_u.get("leaf_zone_confidence")) else None,
                float(r_u["dist_from_crown_norm"]) if pd.notna(r_u.get("dist_from_crown_norm")) else None,
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
                    contour_xy_json = COALESCE(?, contour_xy_json)
                WHERE catalogNumber = ? AND plant_individual_id = ? AND leaf_id = ?;
                """,
                sql_rows,
            )

        # Synchronize companion extracted_leaf_manifest.csv when present
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
                )
                for col_s in sync_cols:
                    if col_s not in m_df.columns:
                        if col_s in {"qc_flag_reason", "leaf_zone"}:
                            m_df[col_s] = ""
                        elif col_s in {"is_crop_clipped", "bilateral_lobe_symmetry", "efa_qc_pass"}:
                            m_df[col_s] = col_s == "efa_qc_pass"
                        elif col_s == "detector_class_id":
                            m_df[col_s] = 0
                        else:
                            m_df[col_s] = np.nan
                for _, r_u in el_df.iterrows():
                    k_u = (str(r_u["catalogNumber"]).strip(), int(r_u["plant_individual_id"]), int(r_u["leaf_id"]))
                    if k_u in key_to_idx:
                        mi = key_to_idx[k_u]
                        m_df.at[mi, "solidity"] = r_u["solidity"]
                        m_df.at[mi, "circularity"] = r_u["circularity"]
                        m_df.at[mi, "detector_class_id"] = int(r_u["detector_class_id"])
                        m_df.at[mi, "is_crop_clipped"] = bool(r_u["is_crop_clipped"])
                        m_df.at[mi, "hole_area_frac"] = r_u["hole_area_frac"]
                        m_df.at[mi, "enclosed_void_frac"] = r_u["enclosed_void_frac"]
                        m_df.at[mi, "bilateral_lobe_symmetry"] = bool(r_u["bilateral_lobe_symmetry"])
                        m_df.at[mi, "efa_qc_pass"] = bool(r_u["efa_qc_pass"])
                        m_df.at[mi, "qc_flag_reason"] = str(r_u["qc_flag_reason"])
                        m_df.at[mi, "leaf_zone"] = str(r_u["leaf_zone"])
                        m_df.at[mi, "leaf_zone_confidence"] = r_u["leaf_zone_confidence"]
                        m_df.at[mi, "dist_from_crown_norm"] = r_u["dist_from_crown_norm"]
                _atomic_write_csv(m_df, eff_manifest)
                logger.info(f"Synchronized updated QC and zone columns to {eff_manifest}")

        n_pass = int((el_df["efa_qc_pass"] == 1).sum())
        n_flag = int(n_rows - n_pass)
        n_bas = int((el_df["leaf_zone"].astype(str).str.lower() == "basal").sum())
        n_caul = int((el_df["leaf_zone"].astype(str).str.lower() == "cauline").sum())
        logger.info(
            f"Re-scored {n_rows} leaves: efa_qc_pass={n_pass}, flagged={n_flag}, basal={n_bas}, cauline={n_caul}."
        )
        return {
            "total": n_rows,
            "efa_qc_pass": n_pass,
            "flagged": n_flag,
            "basal": n_bas,
            "cauline": n_caul,
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
