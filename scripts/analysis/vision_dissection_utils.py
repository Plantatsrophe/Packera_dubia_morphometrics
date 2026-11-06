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

from scripts.core.config import resolve_cpu_threads, resolve_gpu_batch_size


def _resolve_dissection_workers(n_items: int, max_cap: int = 16) -> int:
    req_t = resolve_cpu_threads()
    if req_t is not None:
        return max(1, min(int(req_t), max(1, n_items)))
    return min(max_cap, os.cpu_count() or 4, max(1, n_items))


def resolve_step05_loader_params(
    batch_size: Optional[int],
    threads: Optional[int],
    use_cuda: bool,
    n_records: int,
) -> Tuple[int, int]:
    """Resolve DataLoader num_workers and effective GPU/CPU batch size for Step 05 DINOv2."""
    req_t = resolve_cpu_threads(threads)
    if req_t is not None:
        os.environ["PACKERA_MAX_WORKERS"] = str(req_t)
        try:
            torch.set_num_threads(max(1, int(req_t)))
        except Exception:
            pass
        nw = 0 if req_t <= 1 else (min(req_t, 16) if (use_cuda and n_records >= 16) else 0)
    else:
        nw = (min(8, os.cpu_count() or 4) if n_records >= 64 else 0) if use_cuda else 0
    eff_bs = resolve_gpu_batch_size(
        batch_size,
        default_cuda=64,
        default_cpu=32,
        use_cuda=use_cuda,
    )
    return nw, eff_bs


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
        n_workers = _resolve_dissection_workers(len(prim_pairs))
        if n_workers > 1:
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
    n_workers = _resolve_dissection_workers(len(unique_paths)) if len(unique_paths) >= 24 else 1
    if len(unique_paths) >= 24 and n_workers > 1:
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


def _dir_has_image_files(d: Path) -> bool:
    if not d.is_dir():
        return False
    try:
        return any(p.is_file() and p.suffix.lower() in {".jpg", ".png", ".jpeg"} for p in d.iterdir())
    except OSError:
        return False


def resolve_run_scoped_patch_and_manifest_inputs(
    data_root: Path,
    dense_dir: Path,
    basal_dir: Path,
    tier1_mask_dir: Path,
    manifest_csv: Optional[Path],
    known_cats: Optional[set] = None,
    demo: bool = False,
) -> Tuple[Path, Path, Path, Optional[pd.DataFrame]]:
    """
    Resolve run-scoped `runs/<run_id>/data/...` patch/mask directories and leaf manifest DataFrame.

    When a derived run (`--start-step >= 3`) has empty run-local `cropped_patches/` or lacks a
    run-local `extracted_leaf_manifest.csv`, transparently falls back to canonical `<workspace>/data/`
    patch/mask directories and queries `extracted_leaves` from canonical SQLite (`packera_pipeline.sqlite`).
    """
    dr = Path(data_root)
    is_runs_layout = dr.name == "data" and dr.parent.parent.name == "runs"
    ws_data = (dr.parent.parent.parent / "data") if is_runs_layout else None

    if ws_data is not None and ws_data.is_dir():
        if not _dir_has_image_files(dense_dir) and _dir_has_image_files(ws_data / "cropped_patches" / "rosettes_dense"):
            dense_dir = ws_data / "cropped_patches" / "rosettes_dense"
        if not _dir_has_image_files(basal_dir) and _dir_has_image_files(ws_data / "cropped_patches" / "basal_leaves_raw"):
            basal_dir = ws_data / "cropped_patches" / "basal_leaves_raw"
        if not _dir_has_image_files(tier1_mask_dir) and _dir_has_image_files(ws_data / "masks" / "tier1_intact"):
            tier1_mask_dir = ws_data / "masks" / "tier1_intact"
        if (manifest_csv is None or not Path(manifest_csv).is_file()) and (ws_data / "tables" / "extracted_leaf_manifest.csv").is_file():
            manifest_csv = ws_data / "tables" / "extracted_leaf_manifest.csv"

    man_df: Optional[pd.DataFrame] = None
    if manifest_csv is not None and Path(manifest_csv).is_file():
        man_df = pd.read_csv(manifest_csv)
    elif (not demo) and ws_data is not None and (ws_data / "packera_pipeline.sqlite").is_file():
        try:
            from scripts.core.db import get_db_connection
            conn = get_db_connection(workspace_root=ws_data.parent, init_schema=False)
            try:
                man_df = pd.read_sql_query(
                    "SELECT catalogNumber, plant_individual_id, leaf_id, assigned_tier, ucs_score, "
                    "solidity, circularity, blade_area_cm2, mask_path, crop_path FROM extracted_leaves",
                    conn,
                )
                if known_cats and not man_df.empty:
                    man_df = man_df[man_df["catalogNumber"].astype(str).str.strip().isin(known_cats)].reset_index(drop=True)
            finally:
                conn.close()
        except Exception:
            man_df = None

    return dense_dir, basal_dir, tier1_mask_dir, man_df


def _is_truthy_val(val: Any) -> bool:
    if val is None or isinstance(val, (list, tuple, dict, set, np.ndarray, pd.Series)):
        return False
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    try:
        if bool(pd.isna(val)):
            return False
    except Exception:
        pass
    return (float(val) != 0.0) if isinstance(val, (int, float, np.number)) else (str(val).strip().lower() in {"true", "1", "1.0", "t", "yes"})


def extract_row_type_info(d: Any) -> Tuple[str, bool, int]:
    """Return `(cleaned_type_status, is_type_specimen, type_rank)` for a voucher record dict/row."""
    import re
    from scripts.core.harvester import classify_nomenclatural_type_rank

    raw_ts = d.get("type_status", d.get("typeStatus", ""))
    ts_str = "" if (raw_ts is None or pd.isna(raw_ts)) else str(raw_ts).strip()
    if ts_str.lower() in {"nan", "none", "null", "na"}:
        ts_str = ""
    t_rank, t_lbl = classify_nomenclatural_type_rank(ts_str)
    exp_b = _is_truthy_val(d.get("is_type_specimen", False))
    if t_rank == 0 and exp_b and not re.search(r"not a type|non-type|false type", ts_str, flags=re.IGNORECASE):
        t_rank, t_lbl = 1, (ts_str or "TYPE")
    return (ts_str if (t_rank >= 1 and ts_str) else (t_lbl if t_rank >= 1 else "")), bool(t_rank >= 1), int(t_rank)


def compute_type_visual_control(
    features: Optional[np.ndarray],
    labels: np.ndarray,
    class_names: List[str],
    records_df: Optional[pd.DataFrame] = None,
    pred_probs: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """Compute 768-d DINOv2 type prototype cosine similarity (`vision_type_cosine_sim`), nearest type taxon, and type margin."""
    y = np.asarray(labels, dtype=int)
    n = len(y)
    if n == 0:
        return {
            "vision_type_cosine_sim": np.array([], dtype=np.float64),
            "vision_nearest_type_taxon": [],
            "vision_type_margin": np.array([], dtype=np.float64),
        }
    X_raw = (
        np.asarray(features, dtype=np.float64)
        if (features is not None and len(features) == n and np.asarray(features).size > 0)
        else (
            np.asarray(pred_probs, dtype=np.float64)
            if (pred_probs is not None and len(pred_probs) == n)
            else np.ones((n, 2), dtype=np.float64)
        )
    )
    if X_raw.ndim == 1:
        X_raw = X_raw[:, None]
    X = X_raw / np.maximum(np.linalg.norm(X_raw, axis=1, keepdims=True), 1e-8)

    type_ranks = np.zeros(n, dtype=int)
    det_tiers = ["Tier_3_Bronze"] * n
    if records_df is not None and len(records_df) == n:
        for idx_r, r_dict in enumerate(records_df.to_dict("records")):
            _, _, rk = extract_row_type_info(r_dict)
            type_ranks[idx_r] = rk
            det_tiers[idx_r] = str(r_dict.get("determiner_tier", "Tier_3_Bronze") or "Tier_3_Bronze").strip()
    is_type_m = type_ranks >= 1
    is_gold_m = (
        np.array([t in {"Tier_0_Type_Anchor", "Tier_0_Primary_Type", "Tier_0_Secondary_Type", "Tier_1_Gold"} for t in det_tiers], dtype=bool)
        | is_type_m
    )

    u_cls = sorted(set(int(c) for c in np.unique(y) if 0 <= int(c) < max(1, len(class_names))))
    if not u_cls:
        u_cls = [0]
    prototypes = []
    proto_taxa = []
    for c in u_cls:
        c_mask = y == c
        t_mask = c_mask & is_type_m
        if np.any(t_mask):
            t_idx = np.where(t_mask)[0]
            w_t = np.where(type_ranks[t_idx] == 2, 2.0, 1.0).astype(np.float64)
            proto_type = np.average(X[t_idx], axis=0, weights=w_t)
            proto_type = proto_type / max(float(np.linalg.norm(proto_type)), 1e-8)
            if len(t_idx) == 1:
                g_non_type_mask = c_mask & is_gold_m & (~is_type_m)
                if int(np.sum(g_non_type_mask)) >= 3:
                    X_g = X[g_non_type_mask]
                    med_g_idx = int(np.argmax(np.mean(X_g @ X_g.T, axis=1)))
                    med_g = X_g[med_g_idx] / max(float(np.linalg.norm(X_g[med_g_idx])), 1e-8)
                    proto = 0.80 * proto_type + 0.20 * med_g
                    proto = proto / max(float(np.linalg.norm(proto)), 1e-8)
                else:
                    proto = proto_type
            else:
                proto = proto_type
        else:
            g_mask = c_mask & is_gold_m
            sub_X = X[g_mask] if np.any(g_mask) else (X[c_mask] if np.any(c_mask) else X)
            med = np.median(sub_X, axis=0)
            med_n = med / max(float(np.linalg.norm(med)), 1e-8)
            if len(sub_X) >= 3:
                best_i = int(np.argmax(sub_X @ med_n))
                proto = 0.5 * med_n + 0.5 * sub_X[best_i]
            else:
                proto = med
        proto = proto / max(float(np.linalg.norm(proto)), 1e-8)
        prototypes.append(proto)
        proto_taxa.append(class_names[c] if c < len(class_names) else f"Class_{c}")

    P_mat = np.vstack(prototypes)
    sim_mat = np.clip(X @ P_mat.T, 0.0, 1.0)
    c_to_col = {c: col_i for col_i, c in enumerate(u_cls)}
    sim_given = np.array([sim_mat[i, c_to_col.get(int(y[i]), 0)] for i in range(n)], dtype=np.float64)
    best_col = np.argmax(sim_mat, axis=1)
    nearest_taxa = [proto_taxa[int(ci)] for ci in best_col]
    if sim_mat.shape[1] >= 2:
        s_sorted = np.sort(sim_mat, axis=1)[:, ::-1]
        margin = np.clip(s_sorted[:, 0] - s_sorted[:, 1], 0.0, 1.0)
    else:
        margin = np.clip(sim_given, 0.0, 1.0)
    return {
        "vision_type_cosine_sim": np.round(sim_given, 4),
        "vision_nearest_type_taxon": nearest_taxa,
        "vision_type_margin": np.round(margin, 4),
    }


def fit_calibrated_oof_fold_classifier(
    X_train: np.ndarray,
    y_train: np.ndarray,
    w_train: np.ndarray,
    min_cls: int,
    random_state: int = 42,
) -> Any:
    """Fit cross-validated L2 logistic regression and propagate optimal `C_` into sigmoid Platt calibration with bounded sample weights."""
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.linear_model import LogisticRegression, LogisticRegressionCV

    w_arr = np.asarray(w_train, dtype=np.float64)
    w_mean = max(float(np.mean(w_arr)) if len(w_arr) > 0 else 1.0, 1e-8)
    w_norm = np.clip(w_arr / w_mean, 0.35, 3.0)
    if min_cls >= 2:
        base_clf: Any = LogisticRegressionCV(
            Cs=10,
            cv=min(3, min_cls),
            max_iter=2000,
            class_weight="balanced",
            random_state=random_state,
            scoring="neg_log_loss",
        )
    else:
        base_clf = LogisticRegression(
            max_iter=2000, C=0.25, class_weight="balanced", random_state=random_state
        )
    base_clf.fit(X_train, y_train, sample_weight=w_norm)
    c_attr = getattr(base_clf, "C_", None)
    best_c = (
        float(np.clip(np.median(np.asarray(c_attr, dtype=np.float64)), 0.01, 10.0))
        if c_attr is not None
        else 0.25
    )
    if min_cls >= 3:
        try:
            cal_est = LogisticRegression(
                max_iter=2000, C=best_c, class_weight="balanced", random_state=random_state
            )
            clf = CalibratedClassifierCV(
                estimator=cal_est, method="sigmoid", cv=min(3, min_cls)
            )
            clf.fit(X_train, y_train, sample_weight=w_norm)
            return clf
        except Exception:
            return base_clf
    return base_clf


def overlay_chromatic_type_anchors(
    ax: Any,
    df: pd.DataFrame,
    lum_vals: np.ndarray,
    sat_vals: np.ndarray,
) -> None:
    """Overplot nomenclatural type specimens in foliar chromatic space using hollow diamond markers."""
    if df is None or len(df) == 0 or "is_type_specimen" not in df.columns:
        return
    type_mask = np.array(
        [
            bool(v)
            if isinstance(v, (bool, np.bool_))
            else str(v).strip().lower() in {"true", "1", "yes", "t"}
            for v in df["is_type_specimen"]
        ],
        dtype=bool,
    )
    if np.any(type_mask):
        ax.scatter(
            lum_vals[type_mask],
            sat_vals[type_mask],
            facecolors="none",
            edgecolors="#111111",
            marker="D",
            s=78,
            linewidths=1.6,
            zorder=6,
            label="Nomenclatural Type",
        )



