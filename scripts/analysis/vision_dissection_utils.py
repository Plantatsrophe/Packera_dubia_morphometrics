"""Dissection-neutral basal leaf linking and token-weight utilities for Vision XAI."""

from __future__ import annotations

import hashlib
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from PIL import Image


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


def resolve_specimen_keys_fast(
    grouped_by_key: Dict[Tuple[str, int], Any],
    grouped_by_cat: Dict[str, Any],
    tier1_by_key: Dict[Tuple[str, int], Any],
    tier1_by_cat: Dict[str, Any],
    man_pids_by_cat: Dict[str, set],
    voucher_pair_map: Dict[Tuple[str, int], Dict[str, Any]],
    voucher_map: Dict[str, Dict[str, Any]],
    has_v_pid: bool,
    parse_pid_fn: Callable[[Any, int], int],
) -> Tuple[List[Tuple[str, int]], set, set]:
    """Resolve `(catalogNumber, plant_individual_id)` keys and catalog lookup sets in O(N) time."""
    pids_by_cat: Dict[str, set] = {}
    grp_key_cats: set = set()
    for ck, p in grouped_by_key:
        grp_key_cats.add(ck)
        pids_by_cat.setdefault(ck, set()).add(p)
    t1_key_cats: set = set()
    for ck, p in tier1_by_key:
        t1_key_cats.add(ck)
        pids_by_cat.setdefault(ck, set()).add(p)
    for ck, pset in man_pids_by_cat.items():
        if pset:
            pids_by_cat.setdefault(ck, set()).update(pset)
    if has_v_pid:
        for ck, p in voucher_pair_map:
            pids_by_cat.setdefault(ck, set()).add(p)
    cats_seen = list(
        dict.fromkeys(
            [k[0] for k in grouped_by_key]
            + list(grouped_by_cat)
            + [k[0] for k in tier1_by_key]
            + list(tier1_by_cat)
            + list(man_pids_by_cat)
        )
    )
    all_keys: List[Tuple[str, int]] = []
    for c in cats_seen:
        pset = pids_by_cat.get(c)
        pids = sorted(pset) if pset else [parse_pid_fn(voucher_map.get(c, {}).get("plant_individual_id", 0), 0)]
        for pid in pids:
            all_keys.append((c, pid))
    return all_keys, t1_key_cats, grp_key_cats


def batch_classify_surface_orientations(
    prim_pairs: Sequence[Tuple[str, str]],
    load_mask_fn: Callable[[Optional[Path], Tuple[int, int]], Optional[np.ndarray]],
    classify_fn: Callable[..., Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Decode each primary RGB crop once and classify foliar surface orientation across threads."""
    if not prim_pairs:
        return []

    def _eval_one(pair: Tuple[str, str]) -> Dict[str, Any]:
        prim_patch, prim_mask = pair
        pp = Path(str(prim_patch)) if prim_patch else None
        if pp is None or not pp.is_file():
            return classify_fn(prim_patch, plant_mask=None)
        with Image.open(pp) as im_raw:
            img_arr = np.array(im_raw.convert("RGB"), dtype=np.uint8)
        mp = Path(str(prim_mask)) if (prim_mask and str(prim_mask).strip()) else None
        p_mask_arr = load_mask_fn(mp, img_arr.shape[:2]) if (mp is not None and mp.is_file() and img_arr.size > 0) else None
        return classify_fn(img_arr, plant_mask=p_mask_arr)

    if len(prim_pairs) >= 24:
        n_workers = min(16, os.cpu_count() or 4, len(prim_pairs))
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            return list(pool.map(_eval_one, prim_pairs))
    return [_eval_one(pair) for pair in prim_pairs]


def compute_batch_patch_set_sha256(
    spec_order: Sequence[Tuple[str, int]],
    patches_by_key: Dict[Tuple[str, int], List[str]],
    weights_by_key: Dict[Tuple[str, int], List[float]],
    file_sha_fn: Callable[[Any], Optional[str]],
) -> Dict[Tuple[str, int], Optional[str]]:
    """Hash unique patch files in parallel and assemble per-specimen weighted SHA-256 keys."""
    unique_paths = list(dict.fromkeys(str(p) for gk in spec_order for p in patches_by_key.get(gk, []) if str(p)))
    if len(unique_paths) >= 24:
        n_workers = min(16, os.cpu_count() or 4, len(unique_paths))
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            sha_lookup = dict(zip(unique_paths, pool.map(file_sha_fn, unique_paths)))
    else:
        sha_lookup = {p: file_sha_fn(p) for p in unique_paths}
    return {
        gk: compute_weighted_patch_set_sha256(patches_by_key.get(gk, []), weights_by_key.get(gk, []), sha_lookup.get)
        for gk in spec_order
    }


def pool_patch_tokens_and_fg_ratio_fast(
    ptoks: torch.Tensor,
    fg_grids: torch.Tensor,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute float32 foreground-masked L2-normalized pooled tokens and FG energy ratio on device."""
    toks = ptoks.detach().to(dtype=torch.float32)
    b, n_patches, _ = toks.shape
    m_t = fg_grids.detach().to(device=toks.device, dtype=torch.float32)
    mask_flat = (
        m_t.reshape(b, -1)
        if (m_t.ndim >= 3 or (m_t.ndim == 2 and m_t.shape == (b, n_patches)))
        else m_t.reshape(1, -1).expand(b, -1)
    )
    if mask_flat.shape[1] != n_patches:
        mask_flat = F.interpolate(mask_flat.unsqueeze(1), size=n_patches, mode="nearest")[:, 0]
    weights = torch.clamp(mask_flat, min=0.0)
    w_sums = weights.sum(dim=1, keepdim=True)
    zero_rows = w_sums[:, 0] <= 1e-8
    if bool(torch.any(zero_rows)):
        weights = weights.clone()
        weights[zero_rows, :] = 1.0
        w_sums = weights.sum(dim=1, keepdim=True)
    pooled = (toks * weights.unsqueeze(-1)).sum(dim=1) / torch.clamp(w_sums, min=1e-8)
    pooled_l2 = pooled / torch.clamp(torch.linalg.norm(pooled, dim=1, keepdim=True), min=1e-8)
    tnorms = torch.linalg.norm(toks, dim=-1)
    fg_bin = (mask_flat > 0.0).to(dtype=torch.float32)
    b_fg = torch.clamp((tnorms * fg_bin).sum(dim=1) / torch.clamp(tnorms.sum(dim=1), min=1e-8), 0.0, 1.0)
    return pooled_l2.cpu().numpy().astype(np.float32), b_fg.cpu().numpy()
