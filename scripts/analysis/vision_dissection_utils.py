"""Dissection-neutral basal leaf linking and token-weight utilities for Vision XAI."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd


def _safe_pid_val(val: Any, default: int = 0) -> int:
    s = "" if (val is None or pd.isna(val)) else str(val).strip()
    try:
        return default if (not s or s.lower() in {"nan", "na", "none", "null", "<na>"}) else int(float(s))
    except Exception:
        return default


def _safe_float(val: Any, default: float = float("nan")) -> float:
    try:
        if val is None or pd.isna(val):
            return default
        return float(val)
    except Exception:
        return default


def filter_basal_manifest_rows_py(
    man_df: pd.DataFrame,
    tier_col: str,
    voucher_map: Optional[Dict[str, Any]] = None,
    min_cy_plant: float = 0.40,
    min_cy_sheet: float = 0.42,
    fallback_band: float = 0.15,
) -> pd.DataFrame:
    """Filter manifest Tier-1 rows to homologous basal lamina organs with lowest-lamina fallback."""
    if man_df is None or man_df.empty or not tier_col or tier_col not in man_df.columns:
        return pd.DataFrame()
    t1 = man_df[man_df[tier_col].astype(str).str.match(r"^tier[_\s]*1", case=False, na=False)].copy()
    if t1.empty or "catalogNumber" not in t1.columns:
        return pd.DataFrame()
    t1["_cat_s"] = t1["catalogNumber"].astype(str).str.strip()
    if voucher_map:
        t1 = t1[t1["_cat_s"].isin(voucher_map)].copy()
    if t1.empty:
        return t1
    has_pid = "plant_individual_id" in t1.columns
    t1["_pid_i"] = [_safe_pid_val(v, 0) for v in t1["plant_individual_id"]] if has_pid else 0

    n = len(t1)
    sol = np.array([_safe_float(v) for v in t1["solidity"]], dtype=np.float64) if "solidity" in t1.columns else np.full(n, np.nan)
    circ = np.array([_safe_float(v) for v in t1["circularity"]], dtype=np.float64) if "circularity" in t1.columns else np.full(n, np.nan)
    # Without 120-vertex contour coordinates in manifest-only mode, reject only extreme thin petiole/stem ribbons
    is_ribbon = np.isfinite(sol) & np.isfinite(circ) & (sol < 0.65) & (circ < 0.055)

    ba = np.array([_safe_float(v) for v in t1["blade_area_cm2"]], dtype=np.float64) if "blade_area_cm2" in t1.columns else np.full(n, np.nan)
    min_ba = (
        np.array([_safe_float(v) for v in t1["effective_min_blade_area_cm2"]], dtype=np.float64)
        if "effective_min_blade_area_cm2" in t1.columns
        else np.full(n, np.nan)
    )
    area_ok = np.where(np.isfinite(ba) & np.isfinite(min_ba) & (min_ba > 0.0), ba >= min_ba, True)

    cy_p = np.array([_safe_float(v) for v in t1["cy_rel_plant_bbox"]], dtype=np.float64) if "cy_rel_plant_bbox" in t1.columns else np.full(n, np.nan)
    cy_s = np.array([_safe_float(v) for v in t1["cy_rel_sheet"]], dtype=np.float64) if "cy_rel_sheet" in t1.columns else np.full(n, np.nan)

    keep = np.zeros(n, dtype=bool)
    grp_series = list(zip(t1["_cat_s"].tolist(), t1["_pid_i"].tolist()))
    idx_by_grp: Dict[Tuple[str, int], List[int]] = {}
    for pos_i, gk in enumerate(grp_series):
        idx_by_grp.setdefault(gk, []).append(pos_i)

    for g_pos in idx_by_grp.values():
        g_arr = np.asarray(g_pos, dtype=int)
        lam_ok = g_arr[(~is_ribbon[g_arr]) & area_ok[g_arr]]
        if lam_ok.size == 0:
            lam_ok = g_arr[~is_ribbon[g_arr]]
        if lam_ok.size == 0:
            continue
        has_p = np.isfinite(cy_p[lam_ok])
        has_s = (~has_p) & np.isfinite(cy_s[lam_ok])
        if not np.any(has_p) and not np.any(has_s):
            keep[lam_ok] = True
            continue
        in_primary = np.where(has_p, cy_p[lam_ok] >= min_cy_plant, np.where(has_s, cy_s[lam_ok] >= min_cy_sheet, True))
        if np.any(in_primary):
            keep[lam_ok[in_primary]] = True
        else:
            if np.any(has_p):
                c_max = float(np.nanmax(cy_p[lam_ok[has_p]]))
                fb_hit = has_p & (cy_p[lam_ok] >= (c_max - fallback_band))
            else:
                c_max = float(np.nanmax(cy_s[lam_ok[has_s]]))
                fb_hit = has_s & (cy_s[lam_ok] >= (c_max - fallback_band))
            if np.any(fb_hit):
                keep[lam_ok[fb_hit]] = True

    if not np.any(keep):
        keep[:] = True
    return t1.iloc[np.where(keep)[0]].drop(columns=["_cat_s", "_pid_i"], errors="ignore")


def compute_dissection_neutral_leaf_weight(mrow: Any) -> Tuple[float, float]:
    """Return `(quality_weight, solidity)` using `ucs_score` without penalizing low-solidity dissected blades."""
    sol = _safe_float(mrow.get("solidity", 0.85), default=0.85)
    if not np.isfinite(sol):
        sol = 0.85
    ucs_raw = mrow.get("ucs_score", mrow.get("ucs", 0.90))
    ucs = _safe_float(ucs_raw, default=0.90)
    if not np.isfinite(ucs) or ucs <= 0.0:
        ucs = 0.90
    lqw = _safe_float(mrow.get("leaf_quality_weight", np.nan))
    weight = lqw if (np.isfinite(lqw) and lqw > 0.0) else ucs
    return max(1e-3, float(weight)), float(sol)


def order_basal_leaves_by_median_solidity(
    leaf_items: Sequence[Tuple[float, Path, Path, float]],
) -> List[Tuple[float, Path, Path]]:
    """Order basal leaves so `item[0]` represents the specimen's median basal dissection & highest UCS."""
    if not leaf_items:
        return []
    if len(leaf_items) == 1:
        it0 = leaf_items[0]
        return [(float(it0[0]), it0[1], it0[2])]
    med_sol = float(np.median([float(it[3]) for it in leaf_items]))
    ordered = sorted(
        leaf_items,
        key=lambda it: (abs(float(it[3]) - med_sol) / max(float(it[0]), 0.15), -float(it[0]), str(it[1])),
    )
    return [(float(it[0]), it[1], it[2]) for it in ordered]


def compute_weighted_patch_set_sha256(
    patch_paths: Sequence[str],
    weights: Sequence[float],
    file_sha_fn: Callable[[Any], Optional[str]],
) -> Optional[str]:
    """Compute cache fingerprint across patch files and normalized multi-leaf weights."""
    ps = [str(p) for p in patch_paths if str(p)]
    if not ps:
        return None
    if len(ps) == 1:
        return file_sha_fn(ps[0])
    w_arr = np.array([max(1e-6, float(weights[i])) if i < len(weights) else 1.0 for i in range(len(ps))], dtype=np.float64)
    w_norm = w_arr / max(float(np.sum(w_arr)), 1e-12)
    payload = "|".join(f"{file_sha_fn(p) or ''}:{w_norm[i]:.4f}" for i, p in enumerate(ps))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()
