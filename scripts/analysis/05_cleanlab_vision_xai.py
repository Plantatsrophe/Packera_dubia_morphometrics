#!/usr/bin/env python3
"""Phase 5 Indumentum-Blind Deep Vision & Explainable AI Pipeline (05_cleanlab_vision_xai.py)."""

from __future__ import annotations
import argparse, hashlib, logging, os, re, sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, numpy as np, pandas as pd, torch, torch.nn.functional as F
from PIL import Image, ImageFilter
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, Dataset
from torchvision import transforms
try: import cv2
except ImportError: cv2 = None
torch.manual_seed(42); np.random.seed(42)
if torch.cuda.is_available(): torch.cuda.manual_seed_all(42)
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path: sys.path.insert(0, str(_script_root))
from scripts.core.db import compute_file_sha256, fetch_cached_dinov2_features, get_db_connection, is_canonical_or_explicit_db_context, persist_vision_xai_audit, store_dinov2_features_cache

FOREGROUND_ATTENTION_MIN_RATIO: float = 0.65
DYNAMIC_OKABE_ITO_HEX: List[str] = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#882255", "#44AA99", "#332288", "#117733"]
__all__ = ["FOREGROUND_ATTENTION_MIN_RATIO", "DYNAMIC_OKABE_ITO_HEX", "_get_dynamic_okabe_ito_map", "standardize_packera_taxon", "neutralize_mounting_paper", "apply_tier1_leaf_mask_to_rgb", "classify_foliar_surface_orientation", "preprocess_macro_rosette_crop", "compute_foreground_tissue_mask", "foreground_masked_patch_pooling", "RosettePatchDataset", "load_and_link_rosette_patches", "extract_dinov2_embeddings", "compute_out_of_fold_probabilities", "compute_feature_space_ood_scores", "run_confident_learning_audit", "compute_foreground_attention_ratio", "blend_heatmap_on_image", "compute_vit_patch_saliency", "generate_gradcam_panel", "generate_step05_publication_figures", "main"]

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", handlers=[logging.StreamHandler(sys.stdout)])
logger = logging.getLogger("Cleanlab_Vision_XAI")
def _get_dynamic_okabe_ito_map(taxa_list: List[str]) -> Dict[str, str]:
    """Assigns distinct colorblind-safe hex colors deterministically across any active Packera species set."""
    return {str(tx): DYNAMIC_OKABE_ITO_HEX[i % len(DYNAMIC_OKABE_ITO_HEX)] for i, tx in enumerate(taxa_list)}
def standardize_packera_taxon(species_str: Optional[str]) -> str:
    """Standardize Packera scientific names and historical basionyms via scripts.core.harvester."""
    s_orig = "" if (not species_str or pd.isna(species_str)) else str(species_str).strip()
    if not s_orig or s_orig.lower() in {"none", "nan", "null", "na", "unknown"}: return "Unknown"
    try:
        from scripts.core.harvester import standardize_packera_taxon as _harvester_std
        std = _harvester_std(s_orig)
        if std and str(std).strip() and str(std).strip().lower() not in {"none", "nan", "null", "na", "unknown", "unassigned_packera"}: return str(std).strip()
    except Exception: pass
    m = re.match(r"^packera\s*(?:([x×])(?=\s)|([×]))?\s*([a-z-]+)", s_orig.lower()); epithet = m.group(3).lstrip("×-") if m else ""
    return f"Packera {'× ' if (m.group(1) or m.group(2)) else ''}{epithet}" if (m and epithet) else s_orig.split("(")[0].strip()
def _is_valid_packera_taxon(taxon_str: str) -> bool:
    return bool(taxon_str and taxon_str not in {"Unknown", "Unassigned_Packera"} and re.match(r"^(?:Packera\s+(?:[x×]\s*)?[a-z-]+|Morphotype_[A-Za-z0-9_]+|Ecotype_[A-Za-z0-9_]+)$", str(taxon_str)))
def _to_rgb_uint8(image: Image.Image | np.ndarray) -> np.ndarray:
    arr = np.array(image.convert("RGB"), dtype=np.uint8) if isinstance(image, Image.Image) else np.array(image, copy=True, dtype=np.uint8)
    return np.stack([arr] * 3, axis=-1) if arr.ndim == 2 else (arr[:, :, :3] if (arr.ndim == 3 and arr.shape[2] == 4) else arr)
def neutralize_mounting_paper(image: Image.Image | np.ndarray, neutral_color: Tuple[int, int, int] = (128, 128, 128)) -> Image.Image | np.ndarray:
    """Neutralizes herbarium mounting paper background to uniform neutral gray RGB(128, 128, 128)."""
    was_pil, img_arr = isinstance(image, Image.Image), _to_rgb_uint8(image)
    if img_arr.size == 0: return Image.fromarray(img_arr) if was_pil else img_arr
    if cv2 is not None:
        _, plant_mask = cv2.threshold(cv2.GaussianBlur(cv2.cvtColor(img_arr, cv2.COLOR_RGB2GRAY), (5, 5), 0), 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        bg_mask = cv2.morphologyEx(plant_mask, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))) == 0
    else:
        gray = 0.2989 * img_arr[:, :, 0] + 0.5870 * img_arr[:, :, 1] + 0.1140 * img_arr[:, :, 2]; bg_mask = gray > np.percentile(gray, 60)
    neutralized = img_arr.copy(); neutralized[bg_mask] = neutral_color
    return Image.fromarray(neutralized) if was_pil else neutralized
def _load_aligned_binary_mask(mask_path: Optional[Path], target_hw: Tuple[int, int], orig_hw: Optional[Tuple[int, int]] = None) -> Optional[np.ndarray]:
    """Loads a binary leaf mask and strips zero-border window padding when paired with a tight leaf RGB crop."""
    if mask_path is None or not str(mask_path).strip() or not Path(str(mask_path).strip()).is_file(): return None
    try: mp = Path(str(mask_path).strip()); raw_m = cv2.imread(str(mp), cv2.IMREAD_GRAYSCALE) if cv2 is not None else np.array(Image.open(mp).convert("L"))
    except Exception: return None
    if raw_m is None or raw_m.size == 0: return None
    bin_m, (ref_h, ref_w), (mh, mw) = (raw_m > 127).astype(np.uint8), (orig_hw if orig_hw is not None else target_hw), raw_m.shape[:2]
    if (mh, mw) != (ref_h, ref_w) and mh > ref_h and mw > ref_w and np.any(bin_m > 0):
        ys, xs = np.where(bin_m > 0); fy0, fy1, fx0, fx1 = int(ys.min()), int(ys.max()) + 1, int(xs.min()), int(xs.max()) + 1; fh, fw = fy1 - fy0, fx1 - fx0
        if abs(fh - ref_h) <= max(8, int(0.12 * ref_h)) and abs(fw - ref_w) <= max(8, int(0.12 * ref_w)) and (abs(fh - ref_h) + abs(fw - ref_w)) < 0.5 * (abs(mh - ref_h) + abs(mw - ref_w)): bin_m = bin_m[fy0:fy1, fx0:fx1]
    th, tw = int(target_hw[0]), int(target_hw[1])
    if bin_m.shape[:2] != (th, tw): bin_m = (cv2.resize(bin_m, (tw, th), interpolation=cv2.INTER_NEAREST) > 0).astype(np.uint8) if cv2 is not None else (np.array(Image.fromarray(bin_m * 255).resize((tw, th), Image.NEAREST)) > 127).astype(np.uint8)
    return bin_m
def apply_tier1_leaf_mask_to_rgb(rgb_img: np.ndarray | Image.Image, mask_path: Optional[Path], pad_px: int = 5) -> np.ndarray:
    """Applies a dilated binary mask to an RGB crop, replacing background with [242, 239, 230]."""
    out = _to_rgb_uint8(rgb_img); bin_m = _load_aligned_binary_mask(mask_path, out.shape[:2]) if out.size > 0 else None
    if bin_m is None: return out
    p, bin_u8 = max(0, int(pad_px)), bin_m * 255
    dilated = cv2.dilate(bin_u8, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * p + 1, 2 * p + 1)), iterations=1) if (cv2 is not None and p > 0) else bin_u8
    out[dilated == 0] = np.array([242, 239, 230], dtype=out.dtype)
    return out
def preprocess_macro_rosette_crop(image: Image.Image | np.ndarray, target_size: Tuple[int, int] = (224, 224), sigma: float = 1.0) -> np.ndarray:
    """Applies Gaussian anti-aliasing (sigma=1.0) and downsamples strictly to 224x224."""
    arr, (tw, th) = _to_rgb_uint8(image), (int(target_size[0]), int(target_size[1]))
    if arr.size == 0: return np.zeros((th, tw, 3), dtype=np.uint8)
    if sigma > 0.0:
        ksize = max(3, int(2 * round(2.0 * sigma) + 1)); arr = cv2.GaussianBlur(arr, (ksize, ksize), sigmaX=float(sigma), sigmaY=float(sigma)) if cv2 is not None else np.array(Image.fromarray(arr).filter(ImageFilter.GaussianBlur(radius=float(sigma))))
    if arr.shape[:2] != (th, tw):
        interp = cv2.INTER_AREA if (arr.shape[0] >= th and arr.shape[1] >= tw) else cv2.INTER_LINEAR; arr = cv2.resize(arr, (tw, th), interpolation=interp) if cv2 is not None else np.array(Image.fromarray(arr).resize((tw, th), Image.BILINEAR))
    return arr.astype(np.uint8)
def compute_foreground_tissue_mask(rgb_img: Image.Image | np.ndarray, mask_path: Optional[Path] = None, grid_size: Tuple[int, int] = (16, 16), orig_hw: Optional[Tuple[int, int]] = None) -> Tuple[np.ndarray, np.ndarray]:
    """Computes binary plant tissue mask at image resolution and downsampled 16x16 patch grid."""
    arr, (gh, gw) = _to_rgb_uint8(rgb_img), (int(grid_size[0]), int(grid_size[1])); h, w = arr.shape[:2]
    if h == 0 or w == 0: return np.ones((224, 224), dtype=np.uint8), np.ones((gh, gw), dtype=np.float32)
    fg_mask = _load_aligned_binary_mask(mask_path, (h, w), orig_hw=orig_hw)
    if fg_mask is None:
        is_neutral = ((np.abs(arr[:, :, 0].astype(int) - 128) <= 2) & (np.abs(arr[:, :, 1].astype(int) - 128) <= 2) & (np.abs(arr[:, :, 2].astype(int) - 128) <= 2)) | ((np.abs(arr[:, :, 0].astype(int) - 242) <= 3) & (np.abs(arr[:, :, 1].astype(int) - 239) <= 3) & (np.abs(arr[:, :, 2].astype(int) - 230) <= 3))
        if 0.02 <= float(np.mean(is_neutral)) < 0.98: fg_mask = (~is_neutral).astype(np.uint8)
        elif cv2 is not None:
            lab = cv2.cvtColor(arr, cv2.COLOR_RGB2LAB); _, otsu_l = cv2.threshold(cv2.GaussianBlur(lab[:, :, 0], (5, 5), 0), 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
            chroma = np.hypot(lab[:, :, 1].astype(np.float32) - 128.0, lab[:, :, 2].astype(np.float32) - 128.0)
            combined = np.where((otsu_l > 0) | (chroma > max(8.0, float(np.percentile(chroma, 55)))), 255, 0).astype(np.uint8)
            fg_mask = (cv2.morphologyEx(combined, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))) > 0).astype(np.uint8)
        else: fg_mask = ((0.2989 * arr[:, :, 0] + 0.5870 * arr[:, :, 1] + 0.1140 * arr[:, :, 2]) <= np.percentile(arr.mean(axis=2), 60)).astype(np.uint8)
    if int(fg_mask.sum()) == 0: fg_mask = np.ones((h, w), dtype=np.uint8)
    patch_frac = cv2.resize(fg_mask.astype(np.float32), (gw, gh), interpolation=cv2.INTER_AREA) if cv2 is not None else F.adaptive_avg_pool2d(torch.as_tensor(fg_mask, dtype=torch.float32).unsqueeze(0).unsqueeze(0), (gh, gw))[0, 0].numpy()
    fg_grid = (patch_frac >= 0.15).astype(np.float32)
    return fg_mask.astype(np.uint8), (np.where(patch_frac > 0.0, patch_frac, 1.0).astype(np.float32) if float(fg_grid.sum()) <= 0.0 else fg_grid)

def foreground_masked_patch_pooling(patch_tokens: torch.Tensor | np.ndarray, patch_mask_16x16: torch.Tensor | np.ndarray) -> np.ndarray:
    """Performs Foreground-Masked Average Pooling over spatial patch tokens and L2-normalizes."""
    tokens_np = patch_tokens.detach().cpu().numpy() if isinstance(patch_tokens, torch.Tensor) else np.asarray(patch_tokens, dtype=np.float32)
    mask_np = patch_mask_16x16.detach().cpu().numpy() if isinstance(patch_mask_16x16, torch.Tensor) else np.asarray(patch_mask_16x16, dtype=np.float32)
    single_input = tokens_np.ndim == 2; tokens_np = tokens_np[None, :, :] if single_input else tokens_np; b, n_patches, _ = tokens_np.shape
    mask_flat = mask_np.reshape(b, -1) if (mask_np.ndim >= 3 or (mask_np.ndim == 2 and mask_np.shape == (b, n_patches))) else np.broadcast_to(mask_np.reshape(1, -1), (b, mask_np.size))
    if mask_flat.shape[1] != n_patches: mask_flat = F.interpolate(torch.as_tensor(mask_flat, dtype=torch.float32).unsqueeze(1), size=n_patches, mode="nearest")[:, 0].numpy()
    weights = np.maximum(mask_flat, 0.0).astype(np.float32); weight_sums = np.sum(weights, axis=1, keepdims=True)
    if np.any(weight_sums[:, 0] <= 1e-8): weights[weight_sums[:, 0] <= 1e-8, :] = 1.0; weight_sums = np.sum(weights, axis=1, keepdims=True)
    pooled = np.sum(tokens_np * weights[:, :, None], axis=1) / np.maximum(weight_sums, 1e-8)
    pooled_l2 = (pooled / np.maximum(np.linalg.norm(pooled, axis=1, keepdims=True), 1e-8)).astype(np.float32)
    return pooled_l2[0] if single_input else pooled_l2

def classify_foliar_surface_orientation(image: Image.Image | np.ndarray | str | Path, plant_mask: Optional[np.ndarray] = None, neutral_color: Tuple[int, int, int] = (128, 128, 128), luminance_threshold: float = 68.0, saturation_threshold: float = 0.22, brightness_cutoff: float = 75.0, max_dimension: int = 512) -> Dict[str, Any]:
    """Classifies foliar surface orientation (abaxial vs. adaxial) in CIE LAB & HSV space."""
    default_res = {"surface_orientation": "adaxial", "mean_luminance": 0.0, "mean_saturation": 0.0, "high_brightness_ratio": 0.0}
    if isinstance(image, (str, Path)):
        if not Path(image).exists(): return default_res
        img_arr = np.array(Image.open(image).convert("RGB"))
    else: img_arr = _to_rgb_uint8(image)
    if img_arr.size == 0: return default_res
    if max(img_arr.shape[:2]) > max_dimension:
        scale = float(max_dimension) / max(img_arr.shape[:2]); nw, nh = max(1, int(img_arr.shape[1] * scale)), max(1, int(img_arr.shape[0] * scale))
        img_arr = cv2.resize(img_arr, (nw, nh), interpolation=cv2.INTER_AREA) if cv2 is not None else np.array(Image.fromarray(img_arr).resize((nw, nh), Image.BILINEAR))
        if plant_mask is not None: plant_mask = (cv2.resize(plant_mask.astype(np.uint8), (nw, nh), interpolation=cv2.INTER_NEAREST) > 0) if cv2 is not None else (np.array(Image.fromarray(plant_mask.astype(np.uint8)).resize((nw, nh), Image.NEAREST)) > 0)
    if plant_mask is not None: mask = plant_mask > 0
    else:
        is_neutral = ((img_arr[:, :, 0] == neutral_color[0]) & (img_arr[:, :, 1] == neutral_color[1]) & (img_arr[:, :, 2] == neutral_color[2])) | ((img_arr[:, :, 0] == 242) & (img_arr[:, :, 1] == 239) & (img_arr[:, :, 2] == 230))
        mask = (~is_neutral) if (np.mean(is_neutral) >= 0.02 and np.sum(~is_neutral) >= 5) else (compute_foreground_tissue_mask(img_arr)[0] > 0)
    if mask.sum() < 5: return default_res
    if cv2 is not None:
        L_vals = cv2.cvtColor(img_arr, cv2.COLOR_RGB2LAB)[mask][:, 0].astype(np.float64) * (100.0 / 255.0)
        S_vals = cv2.cvtColor(img_arr, cv2.COLOR_RGB2HSV)[mask][:, 1].astype(np.float64) / 255.0
    else:
        rgb = img_arr[mask].astype(np.float64); Y = 0.2126 * (rgb[:, 0] / 255.0) + 0.7152 * (rgb[:, 1] / 255.0) + 0.0722 * (rgb[:, 2] / 255.0)
        L_vals = np.where(Y > 0.008856, 116.0 * (Y ** (1.0 / 3.0)) - 16.0, 903.3 * Y); cmax, cmin = np.max(rgb, axis=1), np.min(rgb, axis=1); S_vals = np.where(cmax > 0, (cmax - cmin) / cmax, 0.0)
    mean_L, mean_S = float(np.mean(L_vals)), float(np.mean(S_vals))
    return {"surface_orientation": "abaxial" if (mean_L > luminance_threshold and mean_S < saturation_threshold) else "adaxial", "mean_luminance": round(mean_L, 4), "mean_saturation": round(mean_S, 4), "high_brightness_ratio": round(float(np.mean(L_vals > brightness_cutoff)), 4)}

class RosettePatchDataset(Dataset):
    """PyTorch Dataset loading 224x224 Gaussian low-pass leaf/rosette crops and 16x16 tissue masks."""
    def __init__(self, records: List[Dict[str, Any]], transform: Optional[transforms.Compose] = None, neutralize_background: bool = True, sigma: float = 1.0, demo: bool = False):
        self.records = records
        self.transform = transform or transforms.Compose([transforms.ToTensor(), transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])])
        self.neutralize_background, self.sigma, self.demo = neutralize_background, sigma, demo
    def __len__(self) -> int: return len(self.records)
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, int, str]:
        rec = self.records[idx]; image_path, mask_path = Path(str(rec["patch_path"])), rec.get("mask_path")
        has_mask = bool(mask_path and str(mask_path).strip() and Path(str(mask_path)).is_file())
        if image_path.exists():
            raw_arr = np.array(Image.open(image_path).convert("RGB")); orig_hw = raw_arr.shape[:2]
            raw_arr = apply_tier1_leaf_mask_to_rgb(raw_arr, Path(str(mask_path)), pad_px=5) if has_mask else (np.array(neutralize_mounting_paper(raw_arr)) if self.neutralize_background else raw_arr)
            macro_arr = preprocess_macro_rosette_crop(raw_arr, target_size=(224, 224), sigma=self.sigma)
            _, fg_grid = compute_foreground_tissue_mask(macro_arr, mask_path=Path(str(mask_path)) if has_mask else None, grid_size=(16, 16), orig_hw=orig_hw)
        elif self.demo: macro_arr, fg_grid = np.full((224, 224, 3), 128, dtype=np.uint8), np.ones((16, 16), dtype=np.float32)
        else: raise RuntimeError(f"Rosette patch image missing at '{image_path}'. Run Phase 2 or pass `--demo`.")
        return self.transform(Image.fromarray(macro_arr)), torch.as_tensor(fg_grid, dtype=torch.float32), int(rec.get("label_idx", 0)), str(rec["catalogNumber"])

def _parse_pid(val: Any, default: int = 0) -> int:
    s = "" if (val is None or pd.isna(val)) else str(val).strip()
    try: return default if (not s or s.lower() in {"nan", "na", "none", "null", "<na>"}) else int(float(s))
    except Exception: return default
def _is_truthy(val: Any) -> bool:
    if val is None or isinstance(val, (list, tuple, dict, set, np.ndarray, pd.Series)): return False
    if isinstance(val, (bool, np.bool_)): return bool(val)
    try:
        if bool(pd.isna(val)): return False
    except Exception: pass
    return (float(val) != 0.0) if isinstance(val, (int, float, np.number)) else (str(val).strip().lower() in {"true", "1", "1.0", "t", "yes"})
def _row_promoted_bronze(d: Any) -> bool:
    return bool(_is_truthy(d.get("is_promoted_bronze")) or _is_truthy(d.get("is_promoted_anchor")) or str(d.get("anchor_tier_source", "") or "").strip() == "Tier_3_Bronze_Promoted" or str(d.get("determiner_tier", "") or "").strip() == "Tier_3_Bronze_Promoted")
def _merge_promoted_bronze_flags(df: pd.DataFrame, morph_csv: Optional[Path] = None) -> pd.DataFrame:
    if df.empty: return df.assign(is_promoted_bronze=pd.Series(dtype=bool)) if "is_promoted_bronze" not in df.columns else df
    out = df.copy(); out["is_promoted_bronze"] = pd.Series([_row_promoted_bronze(r) for r in out.to_dict("records")], index=out.index, dtype=bool)
    if morph_csv is None or not Path(morph_csv).is_file(): return out
    m_df = pd.read_csv(morph_csv)
    if m_df.empty or "catalogNumber" not in m_df.columns or not ({"is_promoted_bronze", "is_promoted_anchor", "anchor_tier_source"} & set(m_df.columns)): return out
    m_sub, has_m_pid = pd.DataFrame({"catalogNumber": m_df["catalogNumber"].astype(str).str.strip(), "_mp": [_row_promoted_bronze(r) for r in m_df.to_dict("records")]}), ("plant_individual_id" in m_df.columns)
    if has_m_pid: m_sub["plant_individual_id"] = [_parse_pid(v, 0) for v in m_df["plant_individual_id"]]
    j_cols = ["catalogNumber", "plant_individual_id"] if has_m_pid else ["catalogNumber"]
    l_df = pd.DataFrame({"catalogNumber": out["catalogNumber"].astype(str).str.strip(), **({"plant_individual_id": ([_parse_pid(v, 0) for v in out["plant_individual_id"]] if "plant_individual_id" in out.columns else 0)} if has_m_pid else {})}, index=out.index)
    out["is_promoted_bronze"] = (out["is_promoted_bronze"].to_numpy(dtype=bool) | l_df.merge(m_sub.drop_duplicates(subset=j_cols, keep="first"), on=j_cols, how="left")["_mp"].fillna(False).astype(bool).to_numpy()).astype(bool)
    return out
def _extract_catalog_number(stem: str, known_catalogs: Optional[set] = None) -> str:
    cleaned = re.sub(r"_(?:p\d+_)?(?:rosette|leaf|patch|crop).*$", "", stem, flags=re.IGNORECASE)
    return (stem if stem in known_catalogs else cleaned) if (known_catalogs and (stem in known_catalogs or cleaned in known_catalogs)) else (cleaned or stem.split("_")[0])
def _extract_pid_from_stem(stem: str, known_catalogs: Optional[set] = None) -> Optional[int]:
    m = None if (known_catalogs and stem in known_catalogs) else re.search(r"_p(\d+)_(?:rosette|leaf|patch|crop)", stem, flags=re.IGNORECASE)
    return int(m.group(1)) if m else None
def load_and_link_rosette_patches(rosette_dir: Path, vouchers_csv: Path, demo: bool = False, manifest_csv: Optional[Path] = None, morphometrics_csv: Optional[Path] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Links vouchers per `(catalogNumber, plant_individual_id)` to verified `Tier_1_Direct` basal leaf crops/masks (or `rosettes_dense` fallbacks) and maps taxa to `0..K-1`."""
    rosette_dir, vouchers_csv = Path(rosette_dir), Path(vouchers_csv)
    vouchers_df = pd.read_csv(vouchers_csv) if vouchers_csv.exists() else pd.DataFrame()
    if not vouchers_df.empty and "is_primary_duplicate" in vouchers_df.columns:
        vouchers_df = vouchers_df.loc[vouchers_df["is_primary_duplicate"].map(lambda v: v is True or str(v).strip().lower() in {"true", "1", "t", "yes"} or pd.isna(v))].reset_index(drop=True)
    has_v_pid, voucher_map, voucher_pair_map = ((not vouchers_df.empty) and ("plant_individual_id" in vouchers_df.columns)), {}, {}
    if not vouchers_df.empty and "catalogNumber" in vouchers_df.columns:
        for _, row in vouchers_df.iterrows():
            c_str, r_dict = str(row["catalogNumber"]).strip(), row.to_dict(); voucher_map.setdefault(c_str, r_dict)
            if has_v_pid: voucher_pair_map[(c_str, _parse_pid(row.get("plant_individual_id"), 0))] = r_dict
    known_cats = set(voucher_map.keys()) if voucher_map else None
    cropped_root = rosette_dir.parent if rosette_dir.name in {"rosettes_dense", "basal_leaves_raw"} else rosette_dir; data_root = cropped_root.parent
    dense_dir, basal_dir = (rosette_dir if rosette_dir.name == "rosettes_dense" else (cropped_root / "rosettes_dense")), (rosette_dir if rosette_dir.name == "basal_leaves_raw" else (cropped_root / "basal_leaves_raw"))
    tier1_mask_dir = next((d for d in (data_root / "masks" / "tier1_intact", vouchers_csv.parent.parent / "masks" / "tier1_intact") if d.is_dir()), data_root / "masks" / "tier1_intact")
    if manifest_csv is None: manifest_csv = next((p for p in (vouchers_csv.parent / "extracted_leaf_manifest.csv", data_root / "tables" / "extracted_leaf_manifest.csv") if p.is_file()), None)
    cand_morph = Path(morphometrics_csv) if (morphometrics_csv is not None and Path(morphometrics_csv).is_file()) else next((p for p in (vouchers_csv.parent / "morphometric_flags.csv", data_root / "tables" / "morphometric_flags.csv") if p.is_file()), None)
    tier1_by_key, tier1_by_cat, man_pids_by_cat, man_has_pid, seen_basal = {}, {}, {}, False, set()
    if manifest_csv is not None and Path(manifest_csv).is_file():
        man_df = pd.read_csv(manifest_csv); man_has_pid = "plant_individual_id" in man_df.columns
        tier_col = "assigned_tier" if "assigned_tier" in man_df.columns else ("tier" if "tier" in man_df.columns else "")
        if not man_df.empty and "catalogNumber" in man_df.columns:
            if man_has_pid:
                for c_r, p_r in zip(man_df["catalogNumber"], man_df["plant_individual_id"]):
                    c_s = str(c_r).strip()
                    if c_s and (not voucher_map or c_s in voucher_map): man_pids_by_cat.setdefault(c_s, set()).add(_parse_pid(p_r, 0))
            if tier_col:
                for _, mrow in man_df[man_df[tier_col].astype(str).str.match(r"^tier[_\s]*1", case=False, na=False)].iterrows():
                    c_num = str(mrow["catalogNumber"]).strip()
                    if voucher_map and c_num not in voucher_map: continue
                    m_pid, lid_raw = (_parse_pid(mrow.get("plant_individual_id"), 0) if man_has_pid else 0), mrow.get("leaf_id", "")
                    lid = str(int(float(lid_raw))) if (pd.notna(lid_raw) and str(lid_raw).strip().replace(".", "", 1).isdigit()) else str(lid_raw).strip()
                    stems = (f"{c_num}_p{m_pid}_leaf{lid}", f"{c_num}_leaf{lid}") if man_has_pid else (f"{c_num}_leaf{lid}",)
                    cp = next((basal_dir / f"{st}{ext}" for st in stems for ext in (".jpg", ".png", ".jpeg") if (basal_dir / f"{st}{ext}").is_file()), None)
                    mp_raw = str(mrow.get("mask_path", "") or "").strip()
                    mp_cand = ([Path(mp_raw), data_root.parent / mp_raw, _script_root / mp_raw] if (mp_raw and mp_raw.lower() != "nan") else []) + [tier1_mask_dir / f"{st}.png" for st in stems]
                    mp = next((cand for cand in mp_cand if cand.is_file()), None)
                    if cp is not None and mp is not None:
                        sol = float(mrow.get("solidity", 0.85) if pd.notna(mrow.get("solidity", 0.85)) else 0.85)
                        ucs = float(mrow.get("ucs_score", mrow.get("ucs", 0.90)) if pd.notna(mrow.get("ucs_score", mrow.get("ucs", 0.90))) else 0.90)
                        qw = float(mrow.get("leaf_quality_weight", sol * ucs) if pd.notna(mrow.get("leaf_quality_weight", np.nan)) else (sol * ucs))
                        seen_basal.add(cp.resolve()); (tier1_by_key.setdefault((c_num, m_pid), []) if man_has_pid else tier1_by_cat.setdefault(c_num, [])).append((qw, cp, mp))
    if basal_dir.is_dir() and tier1_mask_dir.is_dir():
        for ext in ("*.jpg", "*.png", "*.jpeg"):
            for bp in sorted(basal_dir.glob(ext)):
                if bp.stem.lower().endswith("_mask") or bp.resolve() in seen_basal: continue
                c_num, s_pid, mp = _extract_catalog_number(bp.stem, known_cats), _extract_pid_from_stem(bp.stem, known_cats), tier1_mask_dir / f"{bp.stem}.png"
                if (not voucher_map or c_num in voucher_map) and mp.is_file():
                    (tier1_by_key.setdefault((c_num, s_pid), []) if s_pid is not None else tier1_by_cat.setdefault(c_num, [])).append((0.75, bp, mp))
    for d_t1 in (tier1_by_key, tier1_by_cat):
        for k_t1 in d_t1: d_t1[k_t1].sort(key=lambda it: (-it[0], str(it[1])))
    dense_files: List[Path] = []
    for search_dir, use_rglob in ((dense_dir, False), (rosette_dir, False), (rosette_dir, True)):
        if not dense_files and search_dir.exists():
            for ext in ("*.jpg", "*.png", "*.jpeg"):
                for fp in sorted(search_dir.rglob(ext) if use_rglob else search_dir.glob(ext)):
                    if not fp.stem.lower().endswith("_mask") and "capitula" not in fp.parts and (fp.parent != basal_dir or not (tier1_by_cat or tier1_by_key)): dense_files.append(fp)
    grouped_by_key, grouped_by_cat = {}, {}
    for p in dense_files:
        cat_num, s_pid = _extract_catalog_number(p.stem, known_cats), _extract_pid_from_stem(p.stem, known_cats)
        if not voucher_map or cat_num in voucher_map: (grouped_by_key.setdefault((cat_num, s_pid), []) if s_pid is not None else grouped_by_cat.setdefault(cat_num, [])).append(p)
    def _wh(fp: Path) -> Tuple[int, int]:
        with Image.open(fp) as _im: return _im.size
    cats_seen = list(dict.fromkeys([k[0] for k in grouped_by_key] + list(grouped_by_cat) + [k[0] for k in tier1_by_key] + list(tier1_by_cat) + list(man_pids_by_cat)))
    all_keys = [(c, pid) for c in cats_seen for pid in (sorted({p for (ck, p) in list(grouped_by_key) + list(tier1_by_key) if ck == c} | man_pids_by_cat.get(c, set()) | ({p for (ck, p) in voucher_pair_map if ck == c} if has_v_pid else set())) or [_parse_pid(voucher_map.get(c, {}).get("plant_individual_id", 0), 0)])]
    records: List[Dict[str, Any]] = []
    if all_keys and (grouped_by_key or grouped_by_cat or tier1_by_key or tier1_by_cat):
        for cat_num, pid in all_keys:
            meta = voucher_pair_map.get((cat_num, pid), voucher_map.get(cat_num, {}))
            t1_list = list(tier1_by_key.get((cat_num, pid), [])) if man_has_pid else (list(tier1_by_key.get((cat_num, pid), [])) + (tier1_by_cat.get(cat_num, []) if (not any(c == cat_num for (c, _) in tier1_by_key) or pid == 0) else []))
            p_list = list(grouped_by_key.get((cat_num, pid), [])) + (grouped_by_cat.get(cat_num, []) if (not any(c == cat_num for (c, _) in grouped_by_key) or pid == 0) else [])
            if not t1_list and not p_list: continue
            raw_sp = meta.get("species_raw", meta.get("species_standardized", meta.get("scientificName", meta.get("species", "Unknown")))); std_taxon = standardize_packera_taxon(raw_sp)
            if not _is_valid_packera_taxon(std_taxon): continue
            ros_masks = []
            for rp in p_list:
                cm = next((c for c in (rp.with_name(f"{rp.stem}_mask.png"), dense_dir / f"{cat_num}_p{pid}_rosette_mask.png", tier1_mask_dir / f"{rp.stem}_mask.png", tier1_mask_dir / f"{cat_num}_rosette_mask.png") if c.is_file()), None)
                if cm is None and tier1_mask_dir.is_dir() and rp.is_file():
                    rwh = _wh(rp); cm = next((lm for pat in (f"{cat_num}_p{pid}_*leaf*.png", f"{cat_num}_leaf*.png") for lm in sorted(tier1_mask_dir.glob(pat)) if _wh(lm) == rwh), None)
                ros_masks.append(str(cm) if cm is not None else "")
            t1_patches, t1_masks, t1_weights = [str(x[1]) for x in t1_list], [str(x[2]) for x in t1_list], [float(x[0]) for x in t1_list]
            ros_patches, ros_masks = ([str(x) for x in p_list], ros_masks) if p_list else (list(t1_patches), list(t1_masks))
            prim_patch, prim_mask = (t1_patches[0] if t1_patches else ros_patches[0]), (t1_masks[0] if t1_masks else (ros_masks[0] if ros_masks else ""))
            p_mask_arr = _load_aligned_binary_mask(Path(prim_mask), _wh(Path(prim_patch))[::-1]) if (prim_mask and Path(prim_mask).is_file() and Path(prim_patch).is_file()) else None
            records.append({"catalogNumber": cat_num, "plant_individual_id": int(pid), "patch_path": prim_patch, "mask_path": prim_mask, "rosette_patches": ros_patches, "rosette_masks": ros_masks, "tier1_leaf_patches": t1_patches, "tier1_leaf_masks": t1_masks, "tier1_leaf_weights": t1_weights, "n_rosette_patches": len(ros_patches), "taxon": std_taxon, "species_raw": raw_sp, "label_idx": 0, "institutionCode": str(meta.get("institutionCode", "UNKNOWN") or "UNKNOWN").strip(), "determiner_tier": str(meta.get("determiner_tier", "Tier_3_Bronze") or "Tier_3_Bronze"), "is_promoted_bronze": _row_promoted_bronze(meta), "spatial_balance_weight": float(meta.get("spatial_balance_weight", 1.0) if pd.notna(meta.get("spatial_balance_weight", 1.0)) else 1.0), **classify_foliar_surface_orientation(prim_patch, plant_mask=p_mask_arr)})
    elif not demo: raise RuntimeError(f"No whole-rosette or basal-leaf patch images found in '{dense_dir}', '{basal_dir}', or '{rosette_dir}'.")
    elif not vouchers_df.empty:
        seen_demo = set()
        for _, row in vouchers_df.iterrows():
            raw_sp = str(row.get("species_raw", row.get("species_standardized", row.get("scientificName", row.get("species", "Unknown"))))); std_taxon = standardize_packera_taxon(raw_sp)
            cat_num = str(row.get("catalogNumber", "")).strip()
            if _is_valid_packera_taxon(std_taxon) and cat_num:
                pids_demo = [_parse_pid(row.get("plant_individual_id", 0), 0)] if has_v_pid else (sorted(man_pids_by_cat.get(cat_num, set())) or [0])
                for pid in pids_demo:
                    if (cat_num, pid) in seen_demo: continue
                    seen_demo.add((cat_num, pid)); synth_p = str(dense_dir / f"{cat_num}_p{pid}_rosette.jpg")
                    records.append({"catalogNumber": cat_num, "plant_individual_id": int(pid), "patch_path": synth_p, "mask_path": "", "rosette_patches": [synth_p], "rosette_masks": [""], "tier1_leaf_patches": [], "tier1_leaf_masks": [], "tier1_leaf_weights": [], "n_rosette_patches": 1, "taxon": std_taxon, "species_raw": raw_sp, "label_idx": 0, "institutionCode": str(row.get("institutionCode", "UNKNOWN") or "UNKNOWN").strip(), "determiner_tier": str(row.get("determiner_tier", "Tier_3_Bronze") or "Tier_3_Bronze"), "is_promoted_bronze": _row_promoted_bronze(row), "spatial_balance_weight": 1.0, "surface_orientation": "adaxial", "mean_luminance": 0.0, "mean_saturation": 0.0, "high_brightness_ratio": 0.0})
    class_to_idx = {taxon: idx for idx, taxon in enumerate(sorted({r["taxon"] for r in records if _is_valid_packera_taxon(r["taxon"])}))}
    for r in records: r["label_idx"] = class_to_idx.get(r["taxon"], 0)
    out_df = pd.DataFrame(records) if records else pd.DataFrame(columns=["catalogNumber", "plant_individual_id", "is_promoted_bronze"])
    if not out_df.empty: out_df["plant_individual_id"] = out_df["plant_individual_id"].astype(int)
    return _merge_promoted_bronze_flags(out_df, cand_morph), class_to_idx

def _aggregate_voucher_embeddings(raw_feats: np.ndarray, expanded_meta: List[Any], raw_fg_ratios: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, List[Any], Dict[Any, float]]:
    if len(raw_feats) == 0: return np.empty((0, 768), dtype=np.float32), np.array([], dtype=int), [], {}
    voucher_groups, voucher_order = {}, []
    fg_list = np.asarray(raw_fg_ratios, dtype=np.float64) if (raw_fg_ratios is not None and len(raw_fg_ratios) == len(raw_feats)) else np.full(len(raw_feats), 0.85, dtype=np.float64)
    for feat_vec, meta_item, fg_val in zip(raw_feats, expanded_meta, fg_list):
        g_key, lbl_idx, weight = ((str(meta_item[0]).strip(), _parse_pid(meta_item[1], 0)), meta_item[2], meta_item[3]) if len(meta_item) == 4 else (((str(meta_item[0][0]).strip(), _parse_pid(meta_item[0][1], 0)) if isinstance(meta_item[0], tuple) and len(meta_item[0]) == 2 else meta_item[0]), meta_item[1], meta_item[2])
        if g_key not in voucher_groups: voucher_groups[g_key] = []; voucher_order.append(g_key)
        voucher_groups[g_key].append((feat_vec, max(1e-6, float(weight)), int(lbl_idx), float(fg_val)))
    agg_feats, agg_labels, fg_map = [], [], {}
    for g_key in voucher_order:
        items = voucher_groups[g_key]; nwts = np.array([it[1] for it in items], dtype=np.float64); nwts = nwts / np.sum(nwts)
        vec = items[0][0].astype(np.float32) if len(items) == 1 else np.sum(np.vstack([it[0] for it in items]) * nwts[:, None], axis=0).astype(np.float32)
        agg_feats.append((vec / max(float(np.linalg.norm(vec)), 1e-8)).astype(np.float32)); agg_labels.append(items[0][2])
        fg_map[g_key] = round(float(np.sum(np.array([it[3] for it in items], dtype=np.float64) * nwts)), 4)
    return np.vstack(agg_feats), np.array(agg_labels, dtype=int), voucher_order, fg_map

def extract_dinov2_embeddings(records: List[Dict[str, Any]], model_name: str = "dinov2_vitb14_reg", device: str = "cuda" if torch.cuda.is_available() else "cpu", batch_size: int = 32, neutralize_background: bool = True, demo: bool = False, return_fg_ratios: bool = False, db_path: Optional[Path] = None, sigma: float = 1.0) -> Any:
    """Extracts 768-d DINOv2 embeddings & foreground patch-energy ratios per `(catalogNumber, plant_individual_id)` with SQLite feature caching."""
    torch.manual_seed(42); np.random.seed(42)
    if not records:
        empty = (np.empty((0, 768), dtype=np.float32), np.array([], dtype=int), [])
        return (*empty, np.array([], dtype=np.float64)) if return_fg_ratios else empty
    def _expand(recs: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Any], Dict[Tuple[str, int], List[str]]]:
        exp_r, exp_m, patches_by_k = [], [], {}
        for rec in recs:
            cat_num, pid, lbl_idx = str(rec["catalogNumber"]).strip(), _parse_pid(rec.get("plant_individual_id", 0), 0), int(rec.get("label_idx", 0))
            t1_p, ros_p = rec.get("tier1_leaf_patches"), rec.get("rosette_patches")
            if isinstance(t1_p, list) and t1_p: p_list, m_list, w_list = t1_p, (rec.get("tier1_leaf_masks") or []), (rec.get("tier1_leaf_weights") or [])
            elif isinstance(ros_p, list) and ros_p: p_list, m_list, w_list = ros_p, (rec.get("rosette_masks") or []), []
            else: p_list, m_list, w_list = [rec.get("patch_path", "")], [rec.get("mask_path", "")], []
            patches_by_k.setdefault((cat_num, pid), []).extend(str(p) for p in p_list)
            for idx_p, p_str in enumerate(p_list):
                m_str = m_list[idx_p] if idx_p < len(m_list) else rec.get("mask_path", "")
                exp_r.append({**rec, "catalogNumber": cat_num, "plant_individual_id": pid, "label_idx": lbl_idx, "patch_path": p_str, "mask_path": m_str})
                exp_m.append(((cat_num, pid), lbl_idx, float(w_list[idx_p]) if idx_p < len(w_list) else 1.0))
        return exp_r, exp_m, patches_by_k
    def _finalize(agg_f: np.ndarray, agg_l: np.ndarray, v_ord: List[Tuple[str, int]], fg_map: Dict[Tuple[str, int], float]) -> Any:
        for rec in records:
            gk = (str(rec["catalogNumber"]).strip(), _parse_pid(rec.get("plant_individual_id", 0), 0))
            if gk in fg_map: rec["foreground_attention_ratio"] = fg_map[gk]
        has_multi_pid = any(p != 0 for _, p in v_ord) or len({c for c, _ in v_ord}) < len(v_ord)
        out_ord = list(v_ord) if has_multi_pid else [c for c, _ in v_ord]
        return (agg_f, agg_l, out_ord, np.array([fg_map.get(k, 0.85) for k in v_ord], dtype=np.float64)) if return_fg_ratios else (agg_f, agg_l, out_ord)
    all_exp_r, all_exp_m, patches_by_key = _expand(records)
    if not any(Path(str(r["patch_path"])).exists() for r in all_exp_r[:5]):
        if not demo: raise RuntimeError("No physical patch images found in `data/cropped_patches/`. Run Phase 2 or pass `--demo`.")
        feats = np.random.randn(len(all_exp_r), 768).astype(np.float32)
        for i, r in enumerate(all_exp_r):
            c_idx = int(r.get("label_idx", 0)) % 15; feats[i, c_idx * 50 : (c_idx + 1) * 50] += 3.0
        synth_fg = np.array([float(r.get("foreground_attention_ratio", 0.85)) for r in all_exp_r], dtype=np.float64)
        return _finalize(*_aggregate_voucher_embeddings(feats, all_exp_m, synth_fg))
    spec_order, spec_lbl, spec_rec = [], {}, {}
    for rec in records:
        gk = (str(rec["catalogNumber"]).strip(), _parse_pid(rec.get("plant_individual_id", 0), 0))
        if gk not in spec_lbl: spec_order.append(gk); spec_lbl[gk] = int(rec.get("label_idx", 0)); spec_rec[gk] = rec
    use_cache, db_conn, cache_hits, sha_by_key = ((not demo) and (db_path is not None or is_canonical_or_explicit_db_context(all_exp_r[0].get("patch_path")))), None, {}, {}
    if use_cache:
        try:
            db_conn = get_db_connection(db_path=db_path, init_schema=True)
            sha_by_key = {gk: (compute_file_sha256(ps[0]) if len(ps) == 1 else hashlib.sha256("|".join(compute_file_sha256(p) or "" for p in ps).encode()).hexdigest()) for gk in spec_order for ps in (patches_by_key.get(gk, []),)}
            cache_hits = fetch_cached_dinov2_features(db_conn, [(c, p, sha_by_key.get((c, p))) for c, p in spec_order], model_name, sigma)
        except Exception: db_conn, cache_hits = None, {}
    uncached_recs, computed_feats, computed_fg = [rec for rec in records if (str(rec["catalogNumber"]).strip(), _parse_pid(rec.get("plant_individual_id", 0), 0)) not in cache_hits], {}, {}
    try:
        if uncached_recs:
            exp_r, exp_m, _ = _expand(uncached_recs)
            loader = DataLoader(RosettePatchDataset(exp_r, neutralize_background=neutralize_background, sigma=sigma, demo=demo), batch_size=batch_size, shuffle=False, num_workers=0)
            model, last_err = None, None
            for candidate_name in (model_name, "dinov2_vitb14"):
                try: model = torch.hub.load("facebookresearch/dinov2", candidate_name); break
                except Exception as exc: last_err = exc
            if model is None:
                if not demo: raise RuntimeError(f"Could not load DINOv2 model weights ('{model_name}'): {last_err}") from last_err
                model = _DemoViTFeatureExtractor()
            model, all_feats, all_fg_ratios = (model.to(device).eval() if hasattr(model, "to") else model), [], []
            with torch.no_grad():
                for tensors, fg_grids, _labels, _cats in loader:
                    tensors = tensors.to(device)
                    if hasattr(model, "forward_features"):
                        ptoks = model.forward_features(tensors)["x_norm_patchtokens"]
                        pooled, tnorms = foreground_masked_patch_pooling(ptoks, fg_grids), torch.norm(ptoks.detach().cpu().float(), dim=-1)
                        fg_flat = (fg_grids.detach().cpu().float().reshape(tnorms.shape[0], -1) > 0).float()
                        if fg_flat.shape[1] != tnorms.shape[1]: fg_flat = (F.interpolate(fg_flat.unsqueeze(1), size=tnorms.shape[1], mode="nearest")[:, 0] > 0).float()
                        b_fg = torch.clamp((tnorms * fg_flat).sum(dim=1) / torch.clamp(tnorms.sum(dim=1), min=1e-8), 0.0, 1.0).numpy()
                    else:
                        out = model(tensors).cpu().numpy()
                        pooled = (out / np.maximum(np.linalg.norm(out, axis=1, keepdims=True), 1e-8)).astype(np.float32)
                        b_fg = (fg_grids.detach().cpu().float().reshape(len(pooled), -1) > 0).float().mean(dim=1).numpy()
                    all_feats.append(pooled); all_fg_ratios.append(b_fg)
            new_f, _, new_ord, new_fg = _aggregate_voucher_embeddings(np.vstack(all_feats) if all_feats else np.empty((0, 768), dtype=np.float32), exp_m, np.concatenate(all_fg_ratios) if all_fg_ratios else None)
            cache_rows = []
            for idx_n, gk in enumerate(new_ord):
                computed_feats[gk], computed_fg[gk], r_meta = new_f[idx_n], new_fg.get(gk, 0.85), spec_rec.get(gk, {})
                cache_rows.append({"catalogNumber": gk[0], "plant_individual_id": gk[1], "patch_sha256": sha_by_key.get(gk), "embedding": new_f[idx_n].tolist(), "foreground_attention_ratio": computed_fg[gk], "is_artifact_driven": computed_fg[gk] < FOREGROUND_ATTENTION_MIN_RATIO, "surface_orientation": r_meta.get("surface_orientation"), "mean_luminance": r_meta.get("mean_luminance"), "mean_saturation": r_meta.get("mean_saturation"), "high_brightness_ratio": r_meta.get("high_brightness_ratio")})
            if db_conn is not None and not demo and cache_rows: store_dinov2_features_cache(db_conn, cache_rows, model_name, sigma, demo=demo)
    finally:
        if db_conn is not None: db_conn.close()
    final_f = np.vstack([np.asarray(cache_hits[gk]["embedding"], dtype=np.float32) if gk in cache_hits else computed_feats[gk] for gk in spec_order])
    final_fg = {gk: (float(cache_hits[gk]["foreground_attention_ratio"]) if (gk in cache_hits and cache_hits[gk].get("foreground_attention_ratio") is not None) else computed_fg.get(gk, 0.85)) for gk in spec_order}
    return _finalize(final_f, np.array([spec_lbl[gk] for gk in spec_order], dtype=int), spec_order, final_fg)
def _resolve_adaptive_anchors(labels: np.ndarray, determiner_tiers: Optional[Any] = None, is_promoted_bronze: Optional[Any] = None, min_gs: int = 8) -> Tuple[np.ndarray, np.ndarray]:
    n = len(labels); promo = np.array([_is_truthy(v) for v in is_promoted_bronze], dtype=bool) if (is_promoted_bronze is not None and len(is_promoted_bronze) == n) else np.zeros(n, dtype=bool)
    if determiner_tiers is None or len(determiner_tiers) != n: return np.ones(n, dtype=bool), np.ones(n, dtype=np.float64)
    t_strs = [str(t).strip() for t in determiner_tiers]; gs_mask = np.array([t in {"Tier_1_Gold", "Tier_2_Silver"} for t in t_strs], dtype=bool)
    promo_cand, anchor_mask = (promo | np.array([t == "Tier_3_Bronze_Promoted" for t in t_strs], dtype=bool)), gs_mask.copy()
    weights = np.array([{"Tier_1_Gold": 3.0, "Tier_2_Silver": 1.5, "Tier_3_Bronze_Promoted": 1.0, "Tier_3_Bronze": 0.5}.get(t, 0.5) for t in t_strs], dtype=np.float64)
    for c in np.unique(labels):
        cls_m = labels == c
        if int(np.sum(gs_mask & cls_m)) < min_gs: use_p = cls_m & promo_cand & (~gs_mask); anchor_mask[use_p] = True; weights[use_p] = 1.0
    return anchor_mask, weights
def compute_out_of_fold_probabilities(features: np.ndarray, labels: np.ndarray, n_splits: int = 5, random_state: int = 42, surface_orientations: Optional[np.ndarray | List[str] | pd.Series] = None, stratify_by_orientation: bool = False, determiner_tiers: Optional[np.ndarray | List[str] | pd.Series] = None, spatial_balance_weights: Optional[np.ndarray | List[float] | pd.Series] = None, institution_codes: Optional[np.ndarray | List[str] | pd.Series] = None, is_promoted_bronze: Optional[np.ndarray | List[bool] | pd.Series] = None) -> Tuple[np.ndarray, float, float]:
    """Computes calibrated out-of-fold probabilities anchored on Gold/Silver (and Gold-poor promoted Bronze) reference specimens."""
    X = np.asarray(features, dtype=np.float64); X, labels = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-8), np.asarray(labels, dtype=int)
    anchor_mask, sample_weights = _resolve_adaptive_anchors(labels, determiner_tiers, is_promoted_bronze)
    if spatial_balance_weights is not None and len(spatial_balance_weights) == len(sample_weights):
        sp_w = np.asarray(spatial_balance_weights, dtype=np.float64); sample_weights = sample_weights * np.where(np.isfinite(sp_w) & (sp_w > 0.0), np.clip(sp_w, 0.20, 3.00), 1.0)
    unique_labels, label_counts = np.unique(labels, return_counts=True); num_classes = max(2, (int(np.max(labels)) + 1) if len(labels) > 0 else 2)
    pred_probs = np.zeros((len(labels), num_classes), dtype=np.float64)
    if len(unique_labels) < 2 or int(np.min(label_counts)) < 2:
        if len(X) == 0: return pred_probs, 1.0, 1.0
        centroid = np.median(X[anchor_mask] if np.any(anchor_mask) else X, axis=0); cos_dist = 1.0 - np.clip(X @ (centroid / max(float(np.linalg.norm(centroid)), 1e-8)), -1.0, 1.0)
        q95 = max(float(np.percentile(cos_dist, 95)) if len(cos_dist) >= 5 else float(np.max(cos_dist)) + 1e-6, 1e-6)
        is_outlier = (cos_dist > q95) & ((np.asarray([str(o).lower() for o in surface_orientations]) != "abaxial") if surface_orientations is not None else True)
        p_self = np.where(is_outlier, np.clip(0.10 * np.exp(-(cos_dist - q95) / q95), 0.01, 0.12), np.clip(1.0 - 0.50 * (cos_dist / q95), 0.50, 0.999))
        for i, lbl in enumerate(labels): c0 = int(lbl) % num_classes; pred_probs[i, c0], pred_probs[i, (c0 + 1) % num_classes] = p_self[i], 1.0 - p_self[i]
        return pred_probs, float(np.mean(~is_outlier)), float(np.mean(~is_outlier))
    u_anc, c_anc = np.unique(labels[anchor_mask], return_counts=True)
    if len(u_anc) < len(unique_labels) or int(np.min(c_anc)) < 2: anchor_mask, u_anc, c_anc = np.ones(len(labels), dtype=bool), unique_labels, label_counts
    anchor_idx, bronze_idx, eff_splits = np.where(anchor_mask)[0], np.where(~anchor_mask)[0], max(2, min(int(n_splits), int(np.min(c_anc))))
    X_anc, y_anc, w_anc = X[anchor_idx], labels[anchor_idx], sample_weights[anchor_idx]; strat_anc, cv_splits = y_anc, None
    if stratify_by_orientation and surface_orientations is not None:
        compound = np.array([f"{l}_{o}" for l, o in zip(y_anc, np.asarray(surface_orientations)[anchor_idx])])
        if all(c >= eff_splits for c in np.unique(compound, return_counts=True)[1]): strat_anc = compound
    if institution_codes is not None and len(institution_codes) == len(labels):
        groups = np.asarray([str(g).strip() if pd.notna(g) and str(g).strip() and str(g).strip().upper() != "UNKNOWN" else f"SINGLETON_{idx}" for idx, g in enumerate(np.asarray(institution_codes)[anchor_idx])])
        if len(np.unique(groups)) >= eff_splits:
            try:
                cand = list(StratifiedGroupKFold(n_splits=eff_splits, shuffle=True, random_state=random_state).split(X_anc, strat_anc, groups))
                if all(len(vi) > 0 and len(np.unique(y_anc[ti])) == len(unique_labels) for ti, vi in cand): cv_splits = cand
            except Exception: cv_splits = None
    if cv_splits is None: cv_splits = list(StratifiedKFold(n_splits=eff_splits, shuffle=True, random_state=random_state).split(X_anc, strat_anc))
    bronze_accum, n_valid_folds, y_true_anc, y_pred_anc = np.zeros((len(bronze_idx), num_classes), dtype=np.float64), 0, [], []
    for train_sub, val_sub in cv_splits:
        scaler = StandardScaler(); X_train, X_val, y_train, y_val = scaler.fit_transform(X_anc[train_sub]), scaler.transform(X_anc[val_sub]), y_anc[train_sub], y_anc[val_sub]
        w_train, min_cls = w_anc[train_sub] / max(float(np.mean(w_anc[train_sub])), 1e-8), int(np.min(np.unique(y_train, return_counts=True)[1]))
        lrcv_kw: Dict[str, Any] = {"use_legacy_attributes": False, "l1_ratios": (0.0,), "scoring": "accuracy"} if "use_legacy_attributes" in LogisticRegressionCV.__init__.__code__.co_varnames else {}
        base_clf: Any = LogisticRegressionCV(Cs=np.logspace(-2, 1, 6), cv=max(2, min(3, min_cls)), class_weight="balanced", max_iter=2000, random_state=random_state, **lrcv_kw) if min_cls >= 2 else LogisticRegression(max_iter=2000, C=0.25, class_weight="balanced", random_state=random_state)
        base_clf.fit(X_train, y_train, sample_weight=w_train); active_clf: Any = base_clf
        if min_cls >= 3:
            try:
                cal_clf = CalibratedClassifierCV(estimator=LogisticRegression(max_iter=2000, C=0.25, class_weight="balanced", random_state=random_state), method="sigmoid", cv=min(3, min_cls))
                cal_clf.fit(X_train, y_train, sample_weight=w_train); active_clf = cal_clf
            except Exception: active_clf = base_clf
        val_probs, classes_arr = active_clf.predict_proba(X_val), np.asarray(active_clf.classes_, dtype=int)
        for col_i, c_idx in enumerate(classes_arr): pred_probs[anchor_idx[val_sub], c_idx] = val_probs[:, col_i]
        if len(bronze_idx) > 0:
            brz_probs = active_clf.predict_proba(scaler.transform(X[bronze_idx])); n_valid_folds += 1
            for col_i, c_idx in enumerate(classes_arr): bronze_accum[:, c_idx] += brz_probs[:, col_i]
        y_true_anc.extend(y_val); y_pred_anc.extend(classes_arr[np.argmax(val_probs, axis=1)])
    if len(bronze_idx) > 0 and n_valid_folds > 0:
        brz_mean = bronze_accum / float(n_valid_folds); pred_probs[bronze_idx, :] = brz_mean / np.maximum(brz_mean.sum(axis=1, keepdims=True), 1e-8)
    pred_probs = pred_probs / np.maximum(pred_probs.sum(axis=1, keepdims=True), 1e-8)
    return pred_probs, float(accuracy_score(y_true_anc, y_pred_anc) if y_true_anc else 1.0), float(f1_score(y_true_anc, y_pred_anc, average="weighted") if y_true_anc else 1.0)
def compute_feature_space_ood_scores(features: np.ndarray, labels: np.ndarray, pred_probs: Optional[np.ndarray] = None, determiner_tiers: Optional[Any] = None, is_promoted_bronze: Optional[np.ndarray | List[bool] | pd.Series] = None) -> np.ndarray:
    """Computes feature-space Out-of-Distribution score (`vision_ood_score` in [0, 1])."""
    n = len(labels)
    if n == 0: return np.array([], dtype=np.float64)
    if features is None or len(features) != n: return np.clip(1.0 - np.max(pred_probs, axis=1), 0.0, 1.0) if (pred_probs is not None and len(pred_probs) == n) else np.zeros(n, dtype=np.float64)
    X = np.asarray(features, dtype=np.float64); X, y = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-8), np.asarray(labels, dtype=int)
    anchor_mask, _ = _resolve_adaptive_anchors(y, determiner_tiers, is_promoted_bronze); anchor_mask = anchor_mask if np.any(anchor_mask) else np.ones(n, dtype=bool)
    centroids = [np.mean(X[anchor_mask & (y == c)], axis=0) for c in np.unique(y[anchor_mask]) if np.any(anchor_mask & (y == c))]
    centroids = [mu / max(float(np.linalg.norm(mu)), 1e-8) for mu in (centroids or [np.mean(X, axis=0)])]
    centroid_ood = np.clip(np.min(1.0 - np.clip(X @ np.vstack(centroids).T, -1.0, 1.0), axis=1), 0.0, 1.0)
    try:
        from cleanlab.outlier import OutOfDistribution
        return np.round(0.5 * centroid_ood + 0.5 * np.clip(1.0 - np.asarray(OutOfDistribution().fit_score(features=X), dtype=np.float64), 0.0, 1.0), 4)
    except Exception: return np.round(centroid_ood, 4)
def _compute_tomentum_taxon_discordance(audit_df: pd.DataFrame, error_threshold: float) -> pd.Series:
    """Organic, taxon-agnostic chromatic vestiture discordance rule in CIE L* luminance & HSV saturation space."""
    abax_mm = (audit_df["surface_orientation"].astype(str).str.lower() == "abaxial") & (audit_df["given_label"] != audit_df["predicted_label"]) & (audit_df["c_error"] > error_threshold)
    lum_s, sat_s = pd.to_numeric(audit_df.get("mean_luminance", pd.Series(np.nan, index=audit_df.index)), errors="coerce"), pd.to_numeric(audit_df.get("mean_saturation", pd.Series(np.nan, index=audit_df.index)), errors="coerce")
    valid_l, valid_s = lum_s.dropna(), sat_s.dropna()
    if (len(valid_l) == 0 or float(valid_l.max() - valid_l.min()) <= 1e-6) and (len(valid_s) == 0 or float(valid_s.max() - valid_s.min()) <= 1e-6): return abax_mm.astype(bool)
    adax = audit_df["surface_orientation"].astype(str).str.lower() != "abaxial"
    anc = (audit_df["determiner_tier"].astype(str).str.strip().isin({"Tier_1_Gold", "Tier_2_Silver", "Tier_3_Bronze_Promoted"}) | audit_df.get("is_promoted_bronze", pd.Series(False, index=audit_df.index)).astype(bool)) if "determiner_tier" in audit_df.columns else pd.Series(True, index=audit_df.index)
    med_L, med_S, q75_L = {}, {}, {}
    for tx in set(audit_df["given_label"].dropna().unique()):
        tx_m = audit_df["given_label"] == tx
        ref_l = (tx_m & adax & anc) if np.any(tx_m & adax & anc & lum_s.notna()) else ((tx_m & adax) if np.any(tx_m & adax & lum_s.notna()) else tx_m)
        ref_s = (tx_m & adax & anc) if np.any(tx_m & adax & anc & sat_s.notna()) else ((tx_m & adax) if np.any(tx_m & adax & sat_s.notna()) else tx_m)
        if np.any(ref_l & lum_s.notna()): med_L[tx], q75_L[tx] = float(np.nanmedian(lum_s[ref_l])), float(np.nanpercentile(lum_s[tx_m & lum_s.notna()], 75))
        if np.any(ref_s & sat_s.notna()): med_S[tx] = float(np.nanmedian(sat_s[ref_s]))
    flags = []
    for idx_i, is_cand in enumerate(abax_mm.to_numpy(dtype=bool)):
        if not is_cand: flags.append(False); continue
        g_tx, p_tx, l_i = audit_df.at[idx_i, "given_label"], audit_df.at[idx_i, "predicted_label"], (float(lum_s.iat[idx_i]) if pd.notna(lum_s.iat[idx_i]) else np.nan)
        cond_a = (p_tx in med_L and g_tx in med_L and med_L[p_tx] >= med_L[g_tx] + 2.0) or (p_tx in med_S and g_tx in med_S and med_S[p_tx] <= med_S[g_tx] - 0.02)
        flags.append(bool(cond_a or (np.isfinite(l_i) and g_tx in q75_L and l_i > q75_L[g_tx])))
    return pd.Series(flags, index=audit_df.index, dtype=bool)
def run_confident_learning_audit(pred_probs: np.ndarray, labels: np.ndarray, records_df: pd.DataFrame, class_names: List[str], error_threshold: float = 0.85, filter_abaxial_artifacts: bool = True, features: Optional[np.ndarray] = None, is_promoted_bronze: Optional[np.ndarray | List[bool] | pd.Series] = None) -> pd.DataFrame:
    """Runs class-conditional Confident Learning (`prune_by_noise_rate`) & OOD audit preserving `(catalogNumber, plant_individual_id)`."""
    labels, pred_probs = np.asarray(labels, dtype=int), np.asarray(pred_probs, dtype=np.float64)
    num_classes = pred_probs.shape[1] if pred_probs.ndim == 2 and pred_probs.shape[1] > 0 else max(1, len(class_names))
    class_thresholds = np.array([float(np.mean(pred_probs[labels == c_j, c_j])) if np.any(labels == c_j) else (1.0 - error_threshold) for c_j in range(num_classes)], dtype=np.float64)
    try:
        if len(np.unique(labels)) < 2 or pred_probs.shape[1] != len(np.unique(labels)): raise ValueError("Single-class or non-contiguous label set.")
        from cleanlab.filter import find_label_issues
        from cleanlab.rank import get_label_quality_scores
        quality_scores, issues = get_label_quality_scores(labels=labels, pred_probs=pred_probs), find_label_issues(labels=labels, pred_probs=pred_probs, filter_by="prune_by_noise_rate", return_indices_ranked_by="self_confidence")
    except Exception:
        predicted_classes, quality_scores = np.argmax(pred_probs, axis=1), np.array([pred_probs[i, labels[i]] for i in range(len(labels))], dtype=np.float64)
        issues = np.where((predicted_classes != labels) & ((np.max(pred_probs, axis=1) >= np.array([class_thresholds[min(int(c), num_classes - 1)] for c in predicted_classes])) | (1.0 - quality_scores > error_threshold)))[0]
    audit_df = records_df.copy().reset_index(drop=True)
    if "plant_individual_id" in audit_df.columns: audit_df["plant_individual_id"] = [_parse_pid(v, 0) for v in audit_df["plant_individual_id"]]
    else: audit_df.insert((list(audit_df.columns).index("catalogNumber") + 1) if "catalogNumber" in audit_df.columns else 0, "plant_individual_id", np.zeros(len(audit_df), dtype=int))
    audit_df["plant_individual_id"] = audit_df["plant_individual_id"].astype(int)
    if len(audit_df) > len(pred_probs) and "catalogNumber" in audit_df.columns and len(audit_df.drop_duplicates(subset=["catalogNumber", "plant_individual_id"])) == len(pred_probs):
        fg_agg = audit_df.groupby(["catalogNumber", "plant_individual_id"], sort=False)["foreground_attention_ratio"].mean() if "foreground_attention_ratio" in audit_df.columns else None
        audit_df = audit_df.drop_duplicates(subset=["catalogNumber", "plant_individual_id"], keep="first").reset_index(drop=True)
        if fg_agg is not None: audit_df["foreground_attention_ratio"] = [round(float(fg_agg.loc[(c, p)]), 4) for c, p in zip(audit_df["catalogNumber"], audit_df["plant_individual_id"])]
    num_records = len(audit_df)
    row_promo = np.array([_row_promoted_bronze(r) for r in audit_df.to_dict("records")], dtype=bool) if num_records > 0 else np.zeros(0, dtype=bool)
    arg_promo = np.array([_is_truthy(v) for v in is_promoted_bronze], dtype=bool) if (is_promoted_bronze is not None and len(is_promoted_bronze) == num_records) else np.zeros(num_records, dtype=bool)
    promo_vec = row_promo | arg_promo
    if "is_promoted_bronze" in audit_df.columns or np.any(promo_vec): audit_df["is_promoted_bronze"] = promo_vec
    given_classes, pred_classes = [class_names[l] if l < len(class_names) else "Unknown" for l in labels], [class_names[i] if i < len(class_names) else "Outlier_Morphotype" for i in np.argmax(pred_probs, axis=1)]
    conf_given = np.array([pred_probs[i, labels[i]] if i < len(pred_probs) and labels[i] < pred_probs.shape[1] else 0.0 for i in range(num_records)], dtype=np.float64)
    sorted_probs = np.sort(pred_probs, axis=1)[:, ::-1] if pred_probs.shape[1] >= 2 else np.column_stack([np.max(pred_probs, axis=1), 1.0 - np.max(pred_probs, axis=1)])
    conf_pred, second_prob, c_error = sorted_probs[:, 0], sorted_probs[:, 1], 1.0 - conf_given
    raw_sp_col = audit_df.get("taxon", audit_df.get("species_raw", audit_df.get("species", pd.Series(given_classes, index=audit_df.index))))
    audit_df["species_raw"], std_sp = audit_df.get("species_raw", audit_df.get("species", given_classes)), pd.Series(raw_sp_col, index=audit_df.index).apply(standardize_packera_taxon)
    audit_df["species_standardized"] = [sp if _is_valid_packera_taxon(sp) else (gc if _is_valid_packera_taxon(gc) else "Unassigned_Packera") for sp, gc in zip(std_sp, given_classes)]
    audit_df["given_label"], audit_df["predicted_label"], audit_df["vision_predicted_label"] = given_classes, pred_classes, pred_classes
    audit_df["confidence_given_class"], audit_df["confidence_predicted_class"] = np.round(conf_given, 4), np.round(conf_pred, 4)
    audit_df["vision_confidence"], audit_df["vision_second_prob"], audit_df["vision_margin"] = np.round(conf_pred, 4), np.round(second_prob, 4), np.round(np.clip(conf_pred - second_prob, 0.0, 1.0), 4)
    audit_df["class_self_confidence_threshold"], audit_df["label_quality_score"], audit_df["c_error"] = np.round([class_thresholds[min(int(l), num_classes - 1)] for l in labels], 4), np.round(quality_scores, 4), np.round(c_error, 4)
    audit_df["oof_probs"] = [np.round(pred_probs[i], 6).tolist() if i < len(pred_probs) else [] for i in range(num_records)]
    if "vision_ood_score" not in audit_df.columns: audit_df["vision_ood_score"] = compute_feature_space_ood_scores(features if features is not None else np.empty((0, 768)), labels, pred_probs, audit_df.get("determiner_tier", None), is_promoted_bronze=promo_vec)
    audit_df["vision_ood_score"] = np.round(pd.to_numeric(audit_df["vision_ood_score"], errors="coerce").fillna(0.0), 4); audit_df["vision_is_ood"] = (audit_df["vision_ood_score"] >= 0.65).astype(bool)
    has_fg_col = "foreground_attention_ratio" in audit_df.columns
    audit_df["foreground_attention_ratio"] = np.round(pd.to_numeric(audit_df.get("foreground_attention_ratio", pd.Series(0.85, index=audit_df.index)), errors="coerce").fillna(1.0 if has_fg_col else 0.85), 4)
    audit_df["is_artifact_driven"] = (audit_df["foreground_attention_ratio"] < FOREGROUND_ATTENTION_MIN_RATIO).astype(bool); audit_df["is_background_shortcut"] = audit_df["is_artifact_driven"].astype(bool)
    if "determiner_tier" in audit_df.columns: _, audit_df["authority_weight"] = _resolve_adaptive_anchors(labels, audit_df["determiner_tier"], is_promoted_bronze=promo_vec)
    is_issue, raw_corrupted = np.zeros(num_records, dtype=bool), (audit_df["c_error"] > error_threshold).to_numpy(dtype=bool); is_issue[issues] = True
    if "surface_orientation" not in audit_df.columns: audit_df["surface_orientation"] = "adaxial"
    if "tomentum_taxon_discordance" in audit_df.columns and audit_df["tomentum_taxon_discordance"].notna().any(): audit_df["tomentum_taxon_discordance"] = pd.Series([_is_truthy(v) for v in audit_df["tomentum_taxon_discordance"]], index=audit_df.index, dtype=bool)
    elif "mean_luminance" in audit_df.columns or "mean_saturation" in audit_df.columns: audit_df["tomentum_taxon_discordance"] = _compute_tomentum_taxon_discordance(audit_df, error_threshold)
    else: audit_df["tomentum_taxon_discordance"] = pd.Series(False, index=audit_df.index, dtype=bool)
    abax_err_cut, final_corrupted, final_issue, triage_actions, reasons = max(float(error_threshold), 0.88), [], [], [], []
    for idx_r, (_, row) in enumerate(audit_df.iterrows()):
        was_corrupted, was_cl = bool(raw_corrupted[idx_r]), bool(is_issue[idx_r]); was_flagged = was_corrupted or was_cl
        orient, fg_ratio = str(row.get("surface_orientation", "adaxial")).strip().lower(), float(row.get("foreground_attention_ratio", 1.0))
        tom_disc, conf_p, c_err_v = bool(row.get("tomentum_taxon_discordance", False)), float(row.get("confidence_predicted_class", 0.0)), float(row.get("c_error", 0.0))
        if filter_abaxial_artifacts and orient == "abaxial":
            if (not tom_disc) and (conf_p >= 0.85) and (c_err_v > abax_err_cut) and (fg_ratio >= FOREGROUND_ATTENTION_MIN_RATIO):
                final_corrupted.append(True); final_issue.append(True); act, rsn = "Prune & Queue for Annotation Triage", f"High-confidence structural reassignment on abaxial lamina (no chromatic discordance): DINOv2 predicts {row['predicted_label']} (conf: {conf_p:.2f}) vs recorded {row['given_label']}"
            else:
                final_corrupted.append(False); final_issue.append(False); act, rsn = ("Retain (Abaxial Inversion Artifact)", "Visual divergence driven by dense abaxial arachnoid tomentum (mounting variant); retained as valid voucher.") if was_flagged else ("Retain", "High label consistency across visual self-supervised embeddings")
        else:
            final_corrupted.append(was_corrupted); final_issue.append(was_flagged)
            if fg_ratio < FOREGROUND_ATTENTION_MIN_RATIO and was_flagged: act, rsn = "Downweight (Background Artifact Attribution)", f"Low foreground attention ratio ({fg_ratio:.2f} < {FOREGROUND_ATTENTION_MIN_RATIO:.2f}); prediction driven by sheet background."
            elif was_corrupted: act, rsn = "Prune & Queue for Annotation Triage", f"Deep vision DINOv2 predicts {row['predicted_label']} (conf: {conf_p:.2f}) vs recorded {row['given_label']}"
            elif was_cl: act, rsn = "Flag for Morphometric Review", f"Confident learning flags potential discordance (quality: {row['label_quality_score']:.2f})"
            else: act, rsn = "Retain", "High label consistency across visual self-supervised embeddings"
        triage_actions.append(act); reasons.append(rsn)
    audit_df["is_cleanlab_issue"], audit_df["is_label_corrupted"], audit_df["is_label_issue"] = is_issue, pd.Series(final_corrupted, index=audit_df.index, dtype=bool), pd.Series(final_issue, index=audit_df.index, dtype=bool)
    audit_df["triage_action"], audit_df["discordance_reason"] = triage_actions, reasons
    return audit_df

def compute_foreground_attention_ratio(heatmap: np.ndarray, fg_mask: np.ndarray) -> float:
    """Computes Pointing Game Foreground Attribution Ratio (Zhang et al. 2018)."""
    hm, mask = np.maximum(np.asarray(heatmap, dtype=np.float64), 0.0), np.asarray(fg_mask)
    if hm.size == 0 or mask.size == 0: return 0.0
    if mask.shape[:2] != hm.shape[:2]:
        mask = cv2.resize(mask.astype(np.uint8), (hm.shape[1], hm.shape[0]), interpolation=cv2.INTER_NEAREST) if cv2 is not None else np.array(Image.fromarray(mask.astype(np.uint8)).resize((hm.shape[1], hm.shape[0]), Image.NEAREST))
    bin_mask, total_energy = (mask > 0).astype(np.float64), float(np.sum(hm))
    return float(np.mean(bin_mask)) if total_energy <= 1e-8 else float(np.clip(float(np.sum(hm * bin_mask)) / total_energy, 0.0, 1.0))

def blend_heatmap_on_image(rgb_image: np.ndarray, heatmap: np.ndarray, alpha: float = 0.52, colormap: str = "inferno", fg_mask: Optional[np.ndarray] = None) -> np.ndarray:
    """Blends a 2D float heatmap [0.0, 1.0] onto an RGB uint8 image using inferno colormap."""
    cmap = matplotlib.colormaps[colormap] if hasattr(matplotlib, "colormaps") else plt.get_cmap(colormap)
    colored_heatmap, (h, w) = (cmap(heatmap)[:, :, :3] * 255.0).astype(np.uint8), rgb_image.shape[:2]
    if colored_heatmap.shape[:2] != (h, w): colored_heatmap = np.array(Image.fromarray(colored_heatmap).resize((w, h), Image.BILINEAR))
    blended = (rgb_image.astype(np.float32) * (1.0 - alpha) + colored_heatmap.astype(np.float32) * alpha).astype(np.uint8)
    bg_cond = (fg_mask <= 0) if (fg_mask is not None and fg_mask.shape[:2] == (h, w)) else (np.max(np.abs(rgb_image.astype(np.int16) - np.array([242, 239, 230], dtype=np.int16)), axis=-1) <= 6)
    if np.any(bg_cond) and not np.all(bg_cond): blended[bg_cond] = rgb_image[bg_cond]
    return blended

class _DemoViTFeatureExtractor(torch.nn.Module):
    """Explicit `--demo` mode patch-token feature extractor following the DINOv2 ViT contract."""
    def __init__(self, embed_dim: int = 768, patch_size: int = 14):
        super().__init__(); self.embed_dim, self.patch_size = embed_dim, patch_size
    def forward_features(self, img_tensor: torch.Tensor) -> Dict[str, torch.Tensor]:
        b, c, h, w = img_tensor.shape; h_p, w_p = max(1, h // self.patch_size), max(1, w // self.patch_size)
        patches = F.adaptive_avg_pool2d(img_tensor, (h_p, w_p)).permute(0, 2, 3, 1).reshape(b, h_p * w_p, c)
        paper_ref = torch.tensor([(242.0 / 255.0 - 0.485) / 0.229, (239.0 / 255.0 - 0.456) / 0.224, (230.0 / 255.0 - 0.406) / 0.225], device=patches.device, dtype=patches.dtype)[:c]
        amp = torch.clamp(torch.norm(patches - paper_ref.view(1, 1, -1), dim=-1, keepdim=True), min=0.02)
        patch_tokens = F.layer_norm(patches.repeat(1, 1, (self.embed_dim + c - 1) // c)[:, :, : self.embed_dim], (self.embed_dim,)) * amp
        return {"x_norm_patchtokens": patch_tokens, "x_norm_clstoken": patch_tokens.mean(dim=1)}

def compute_vit_patch_saliency(img_arr: np.ndarray, model: Any, classifier: Optional[Any] = None, target_class_idx: Optional[int] = None, device: str = "cuda" if torch.cuda.is_available() else "cpu") -> np.ndarray:
    """Projects class weights / Captum attribution onto the 16x16 ViT patch grid upsampled to 224x224."""
    h, w = img_arr.shape[:2]; macro_arr = preprocess_macro_rosette_crop(img_arr, target_size=(224, 224), sigma=1.0)
    img_tensor = transforms.Compose([transforms.ToTensor(), transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])])(Image.fromarray(macro_arr)).unsqueeze(0).to(device)
    with torch.no_grad():
        feat_dict = model.forward_features(img_tensor); patch_tokens, cls_token = feat_dict["x_norm_patchtokens"][0], feat_dict.get("x_norm_clstoken", None)
    num_patches, embed_dim = patch_tokens.shape
    h_p = max(1, int(np.round(np.sqrt(num_patches)))); w_p = max(1, num_patches // h_p); F_patches = patch_tokens[: h_p * w_p].reshape(h_p, w_p, embed_dim)
    w_c, coef_arr, classes_list = None, None, []
    if classifier is not None:
        est = classifier if hasattr(classifier, "coef_") else (getattr(classifier.calibrated_classifiers_[0], "estimator", None) if getattr(classifier, "calibrated_classifiers_", None) else None)
        if est is not None and hasattr(est, "coef_"):
            coef_arr, classes_list = np.asarray(est.coef_, dtype=np.float32), list(getattr(est, "classes_", range(est.coef_.shape[0])))
            if coef_arr.ndim == 1: coef_arr = coef_arr[None, :]
            if coef_arr.shape[0] == 1 and len(classes_list) == 2: coef_arr = np.vstack([-coef_arr, coef_arr])
    if coef_arr is not None and coef_arr.shape[0] > 0:
        if target_class_idx is None and cls_token is not None and hasattr(classifier, "predict"): target_class_idx = int(classifier.predict(cls_token.detach().cpu().numpy())[0])
        row_idx = (classes_list.index(target_class_idx) if target_class_idx in classes_list else int(target_class_idx)) if target_class_idx is not None else 0
        w_c = torch.as_tensor(coef_arr[row_idx % coef_arr.shape[0]], dtype=F_patches.dtype, device=F_patches.device)
    elif cls_token is not None: w_c = cls_token[0].to(dtype=F_patches.dtype, device=F_patches.device)
    raw_saliency = torch.einsum("hwd,d->hw", F_patches, w_c) if (w_c is not None and w_c.numel() == embed_dim) else torch.norm(F_patches, dim=-1)
    saliency_relu = torch.relu(raw_saliency) if float(torch.relu(raw_saliency).max().item()) > 1e-8 else torch.relu(raw_saliency - raw_saliency.min())
    heatmap = torch.clamp(F.interpolate(saliency_relu.unsqueeze(0).unsqueeze(0).to(dtype=torch.float32), size=(h, w), mode="bicubic", align_corners=False)[0, 0], min=0.0).detach().cpu().numpy().astype(np.float32)
    return (heatmap - heatmap.min()) / (heatmap.max() - heatmap.min() + 1e-8)

def generate_gradcam_panel(flagged_records: List[dict], output_path: Path, max_samples: int = 8, model: Optional[Any] = None, classifier: Optional[Any] = None, model_name: str = "dinov2_vitb14_reg", device: str = "cuda" if torch.cuda.is_available() else "cpu", demo: bool = False, class_names: Optional[List[str]] = None, active_taxa: Optional[List[str]] = None) -> None:
    """Renders multi-panel diagnostic figure with 224x224 ViT patch attribution, leaf-only saliency & calibrated OOF posteriors."""
    output_path = Path(output_path); output_path.parent.mkdir(parents=True, exist_ok=True); dedup_records, seen_keys = [], set()
    for r in flagged_records:
        r_key = (str(r.get("catalogNumber", r.get("patch_path", ""))), _parse_pid(r.get("plant_individual_id", 0), 0))
        if r_key not in seen_keys: seen_keys.add(r_key); dedup_records.append(r)
    samples = dedup_records[:max_samples]
    if not samples: return
    obs_taxa = {standardize_packera_taxon(str(r.get(k, ""))) for r in flagged_records for k in ("taxon", "species_standardized", "given_label", "species_raw", "predicted_taxon", "predicted_label", "vision_predicted_label") if r.get(k) and _is_valid_packera_taxon(standardize_packera_taxon(str(r.get(k, ""))))}
    active_taxa = [t for t in (active_taxa or class_names or []) if _is_valid_packera_taxon(t)] or sorted(t for t in obs_taxa if _is_valid_packera_taxon(t)) or ["Unassigned_Packera"]
    palette_map = _get_dynamic_okabe_ito_map(active_taxa)
    if model is None:
        last_err = None
        for cand in (model_name, "dinov2_vitb14"):
            try: model = torch.hub.load("facebookresearch/dinov2", cand).to(device).eval(); break
            except Exception as exc: last_err = exc
        if model is None:
            if not demo: raise RuntimeError(f"Could not load DINOv2 model ('{model_name}'): {last_err}") from last_err
            model = _DemoViTFeatureExtractor()
    n = len(samples); fig, axes = plt.subplots(n, 3, figsize=(11.8, 3.2 * n), gridspec_kw={"width_ratios": [1.0, 1.0, 1.25], "wspace": 0.28})
    if n == 1: axes = np.expand_dims(axes, 0)
    for idx, rec in enumerate(samples):
        cat_num_str, pid_val, raw_p, raw_m = str(rec.get("catalogNumber", f"V{idx + 1}")), _parse_pid(rec.get("plant_individual_id", 0), 0), rec.get("patch_path"), rec.get("mask_path")
        img_path = Path(str(raw_p)) if (raw_p and str(raw_p).strip()) else (Path("data/cropped_patches/rosettes_dense") / f"{cat_num_str}_p{pid_val}_rosette.jpg")
        mask_p, orig_hw = (Path(str(raw_m).strip()) if (raw_m and str(raw_m).strip() and Path(str(raw_m).strip()).is_file()) else None), None
        if img_path.exists():
            raw_arr = np.array(Image.open(img_path).convert("RGB")); orig_hw = raw_arr.shape[:2]
            img_arr = preprocess_macro_rosette_crop(apply_tier1_leaf_mask_to_rgb(raw_arr, mask_p, pad_px=5) if mask_p is not None else raw_arr, target_size=(224, 224), sigma=1.0)
        elif demo:
            raw_demo = np.full((224, 224, 3), fill_value=[60, 120, 40], dtype=np.uint8); img_arr = apply_tier1_leaf_mask_to_rgb(raw_demo, mask_p, pad_px=5) if mask_p is not None else raw_demo
        else: raise RuntimeError(f"Missing patch image at '{img_path}'.")
        fg_mask, _ = compute_foreground_tissue_mask(img_arr, mask_path=mask_p, grid_size=(16, 16), orig_hw=orig_hw)
        rec_taxon = standardize_packera_taxon(str(rec.get("taxon", rec.get("species_standardized", rec.get("given_label", active_taxa[0]))))); rec_taxon = rec_taxon if _is_valid_packera_taxon(rec_taxon) else active_taxa[0]
        raw_pred = rec.get("predicted_taxon", rec.get("predicted_label", rec.get("vision_predicted_label", rec_taxon)))
        pred_taxon = standardize_packera_taxon(str(raw_pred)) if raw_pred else rec_taxon; pred_taxon = pred_taxon if pred_taxon in active_taxa else (rec_taxon if rec_taxon in active_taxa else active_taxa[0])
        heatmap = compute_vit_patch_saliency(img_arr, model=model, classifier=classifier, target_class_idx=(active_taxa.index(pred_taxon) if pred_taxon in active_taxa else int(rec.get("label_idx", 0) or 0)), device=device)
        fg_ratio = compute_foreground_attention_ratio(heatmap, fg_mask); rec["foreground_attention_ratio"] = round(fg_ratio, 4)
        fg_bin = (fg_mask > 0).astype(np.uint8); heatmap_masked = heatmap * fg_bin.astype(np.float32)
        if float(heatmap_masked.max()) > 1e-8: heatmap_masked = heatmap_masked / float(heatmap_masked.max())
        blended = blend_heatmap_on_image(img_arr, heatmap_masked, alpha=0.52, colormap="inferno", fg_mask=fg_bin)
        c_err_val = float(rec.get("c_error", 1.0 - float(rec.get("label_quality_score", 0.25)))); conf_val = float(np.clip(max(c_err_val, 1.0 - c_err_val), 0.35, 0.99))
        axes[idx, 0].imshow(img_arr); axes[idx, 0].axis("off"); axes[idx, 0].set_title(f"[1] {cat_num_str}{f' [p{pid_val}]' if pid_val else ''} (224x224 Macro)\nRecorded: {rec_taxon.replace('Packera ', 'P. ')}", fontsize=8.8, fontweight="bold")
        axes[idx, 1].imshow(blended); axes[idx, 1].axis("off"); axes[idx, 1].set_title(f"[2] 16x16 ViT Saliency (FG={fg_ratio:.2f})\nPredicted: {pred_taxon.replace('Packera ', 'P. ')} (c_err={c_err_val:.2f})", fontsize=8.8, fontweight="bold")
        ax_bar, raw_oof = axes[idx, 2], rec.get("oof_probs")
        if raw_oof is not None and len(raw_oof) == len(active_taxa): probs = np.clip(np.asarray(raw_oof, dtype=float), 0.0, 1.0); probs = probs / max(float(probs.sum()), 1e-8)
        elif raw_oof is not None and class_names and len(raw_oof) == len(class_names):
            pmap = {str(cn): float(pv) for cn, pv in zip(class_names, raw_oof)}; probs = np.array([max(0.0, pmap.get(tx, 0.0)) for tx in active_taxa], dtype=float); probs = probs / max(float(probs.sum()), 1e-8)
        else:
            probs = np.array([float(np.clip(1.0 - c_err_val, 0.01, 0.99))], dtype=float) if len(active_taxa) == 1 else np.full(len(active_taxa), (1.0 - conf_val) / max(len(active_taxa) - 1, 1), dtype=float)
            if len(active_taxa) > 1 and pred_taxon in active_taxa: probs[active_taxa.index(pred_taxon)] = conf_val; probs = probs / max(probs.sum(), 1e-8)
        y_pos = np.arange(len(active_taxa)); ax_bar.barh(y_pos, probs * 100.0, color=[palette_map.get(tx, "#0072B2") for tx in active_taxa], edgecolor="#222222", height=0.58, alpha=0.88)
        ax_bar.set_yticks(y_pos); ax_bar.set_yticklabels([tx.replace("Packera ", "P. ") for tx in active_taxa], fontsize=8.2, fontstyle="italic"); ax_bar.set_xlim(0, 112); ax_bar.invert_yaxis(); ax_bar.set_title("[3] Calibrated Posteriors", fontsize=8.8, fontweight="bold")
    fig.subplots_adjust(hspace=0.38, wspace=0.24); tmp_out = output_path.with_name(f".{output_path.stem}_{os.getpid()}_tmp{output_path.suffix}")
    plt.savefig(tmp_out, dpi=300, bbox_inches="tight"); tmp_out.replace(output_path); pdf_companion = output_path.with_suffix(".pdf")
    if pdf_companion != output_path:
        tmp_pdf = output_path.with_name(f".{output_path.stem}_{os.getpid()}_tmp.pdf"); plt.savefig(tmp_pdf, format="pdf", bbox_inches="tight"); tmp_pdf.replace(pdf_companion)
    plt.close()

def generate_step05_publication_figures(audit_df: pd.DataFrame, output_figure_path: Path) -> Dict[str, Path]:
    """Generates Step 05 publication figure `fig05b_cleanlab_quality_and_tomentum.{pdf,png}`."""
    if audit_df is None or audit_df.empty: return {}
    env_run = os.environ.get("PACKERA_RUN_DIR", "").strip()
    figures_dir = (Path(env_run).resolve() / "outputs" / "figures") if env_run else Path(output_figure_path).resolve().parent
    figures_dir.mkdir(parents=True, exist_ok=True); df = audit_df.copy()
    sp_col, pred_col = ("species_standardized" if "species_standardized" in df.columns else "given_label"), ("predicted_label" if "predicted_label" in df.columns else "vision_predicted_label")
    df["taxon_clean"] = df.get(sp_col, pd.Series("Unassigned_Packera", index=df.index)).astype(str).apply(standardize_packera_taxon)
    df["pred_clean"] = df.get(pred_col, df["taxon_clean"]).astype(str).apply(standardize_packera_taxon)
    taxa_order = sorted({t for t in df["taxon_clean"] if _is_valid_packera_taxon(t)}) or sorted({t for t in df["pred_clean"] if _is_valid_packera_taxon(t)}) or ["Unassigned_Packera"]
    palette_map, short_labels = _get_dynamic_okabe_ito_map(taxa_order), [tx.replace("Packera ", "P. ") for tx in taxa_order]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.2, 5.8), dpi=300)
    conf_counts = np.array([[int(np.sum((df["taxon_clean"] == r_tx) & (df["pred_clean"] == c_tx))) for c_tx in taxa_order] for r_tx in taxa_order], dtype=int)
    ax_a.imshow((conf_counts / np.maximum(conf_counts.sum(axis=1, keepdims=True), 1)) * 100.0, cmap="Blues", vmin=0.0, vmax=100.0, aspect="auto")
    ax_a.set_xticks(range(len(taxa_order))); ax_a.set_xticklabels(short_labels, rotation=30, ha="right", fontsize=8.5, fontstyle="italic")
    ax_a.set_yticks(range(len(taxa_order))); ax_a.set_yticklabels(short_labels, fontsize=8.5, fontstyle="italic")
    ax_a.set_title("A. Cross-Species Visual Confusion Matrix", fontsize=10.5, fontweight="bold", loc="left")
    lum_vals = pd.to_numeric(df.get("mean_luminance", pd.Series(52.0, index=df.index)), errors="coerce").fillna(52.0).to_numpy(dtype=float)
    sat_vals = pd.to_numeric(df.get("mean_saturation", pd.Series(0.30, index=df.index)), errors="coerce").fillna(0.30).to_numpy(dtype=float)
    hollow_markers = ("o", "^", "s", "D", "v", "P", "X", "<", ">", "h")
    for idx_sp, sp_name in enumerate(taxa_order):
        sp_mask = (df["taxon_clean"] == sp_name).to_numpy()
        if np.any(sp_mask): ax_b.scatter(lum_vals[sp_mask], sat_vals[sp_mask], facecolors="none", edgecolors=palette_map.get(sp_name, "#0072B2"), marker=hollow_markers[idx_sp % len(hollow_markers)], s=36, alpha=0.75, label=sp_name.replace("Packera ", "P. "))
    ax_b.set_xlabel("Mean Foliar CIE L* Luminance", fontsize=9.5, fontweight="bold"); ax_b.set_ylabel("Mean Foliar HSV Saturation", fontsize=9.5, fontweight="bold")
    if len(taxa_order) > 0: ax_b.legend(loc="best", fontsize=8.0, frameon=True)
    ax_b.set_title("B. Foliar Chromatic Space & Species Ellipses", fontsize=10.5, fontweight="bold", loc="left"); fig.tight_layout()
    pdf_out, png_out = figures_dir / "fig05b_cleanlab_quality_and_tomentum.pdf", figures_dir / "fig05b_cleanlab_quality_and_tomentum.png"
    tmp_pdf, tmp_png = figures_dir / f".fig05b_{os.getpid()}_tmp.pdf", figures_dir / f".fig05b_{os.getpid()}_tmp.png"
    try:
        fig.savefig(tmp_pdf, format="pdf", bbox_inches="tight"); tmp_pdf.replace(pdf_out)
        fig.savefig(tmp_png, format="png", dpi=300, bbox_inches="tight"); tmp_png.replace(png_out)
    finally: plt.close(fig)
    return {"pdf": pdf_out, "png": png_out}
def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for Phase 5 analysis."""
    p = argparse.ArgumentParser(description="Phase 5: Indumentum-Blind DINOv2 Macro-Architecture & Cleanlab XAI Audit")
    p.add_argument("--rosette-dir", type=Path, default=Path("data/cropped_patches")); p.add_argument("--vouchers-csv", type=Path, default=Path("data/tables/curated_vouchers.csv"))
    p.add_argument("--output-csv", type=Path, default=Path("data/tables/label_noise_audit.csv")); p.add_argument("--output-figure", type=Path, default=Path("outputs/figures/GradCAM_audit_panel.png"))
    p.add_argument("--cleanlab-threshold", dest="cleanlab_threshold", type=float, default=0.85); p.add_argument("--error-threshold", dest="cleanlab_threshold", type=float)
    p.add_argument("--neutralize-background", dest="neutralize_background", action="store_true", default=True); p.add_argument("--no-neutralize-background", dest="neutralize_background", action="store_false")
    p.add_argument("--export-figures", dest="export_figures", action="store_true", default=True); p.add_argument("--no-export-figures", dest="export_figures", action="store_false")
    p.add_argument("--filter-abaxial-artifacts", dest="filter_abaxial_artifacts", action="store_true", default=True); p.add_argument("--no-filter-abaxial-artifacts", dest="filter_abaxial_artifacts", action="store_false")
    p.add_argument("--stratify-by-orientation", dest="stratify_by_orientation", action="store_true", default=False); p.add_argument("--intraspecific-mode", dest="intraspecific_mode", action="store_true", default=False)
    p.add_argument("--morphometrics-csv", type=Path, default=None); p.add_argument("--manifest-csv", type=Path, default=None)
    p.add_argument("--focal-taxon", type=str, default=""); p.add_argument("--demo", dest="demo", action="store_true", default=False); p.add_argument("--db", type=Path, default=None)
    return p.parse_args()
def _apply_intraspecific_morphotype_labels(df: pd.DataFrame, morphometrics_csv: Optional[Path], vouchers_csv: Optional[Path] = None, intraspecific_mode: bool = False) -> pd.DataFrame:
    """Remaps `taxon` and `label_idx` to Step 04 `Morphotype_K1..Kn` clusters when single-species or `--intraspecific-mode`."""
    if df.empty or (not intraspecific_mode and df["taxon"].nunique() >= 2): return df
    m_path = Path(morphometrics_csv) if morphometrics_csv is not None else (Path(vouchers_csv).parent / "morphometric_flags.csv" if vouchers_csv is not None else None)
    if m_path is None or not m_path.is_file(): return df
    m_df = pd.read_csv(m_path)
    if m_df.empty or "catalogNumber" not in m_df.columns or "gmm_cluster" not in m_df.columns: return df
    def _fmt(v: Any) -> Optional[str]:
        s = "" if (v is None or pd.isna(v)) else str(v).strip()
        if not s or s.lower() in {"nan", "na", "none"}: return None
        m_num = re.search(r"(\d+)", s); return s if s.startswith("Morphotype_") else (f"Morphotype_K{int(m_num.group(1))}" if m_num else f"Morphotype_{re.sub(r'[^A-Za-z0-9_]', '_', s)}")
    m_sub, left = m_df.copy(), df.copy(); m_sub["catalogNumber"], left["catalogNumber"] = m_sub["catalogNumber"].astype(str).str.strip(), left["catalogNumber"].astype(str).str.strip()
    m_sub["_m_lbl"] = m_sub["gmm_cluster"].map(_fmt); m_sub = m_sub[m_sub["_m_lbl"].notna()]
    join_keys = ["catalogNumber", "plant_individual_id"] if ("plant_individual_id" in m_sub.columns and "plant_individual_id" in left.columns) else ["catalogNumber"]
    if "plant_individual_id" in join_keys: m_sub["plant_individual_id"], left["plant_individual_id"] = [_parse_pid(v, 0) for v in m_sub["plant_individual_id"]], [_parse_pid(v, 0) for v in left["plant_individual_id"]]
    merged = left.merge(m_sub.drop_duplicates(subset=join_keys, keep="first")[join_keys + ["_m_lbl"]], on=join_keys, how="left")
    if merged["_m_lbl"].dropna().nunique() < 2: return df
    out = merged[merged["_m_lbl"].notna()].copy().reset_index(drop=True); out["taxon"] = out.pop("_m_lbl")
    c2i = {tx: i for i, tx in enumerate(sorted(out["taxon"].unique()))}; out["label_idx"] = [c2i[tx] for tx in out["taxon"]]
    return out
def main() -> None:
    """Orchestrates Phase 5 Indumentum-Blind Deep Vision & Confident Learning workflow."""
    args = parse_args()
    morph_csv = Path(args.morphometrics_csv) if (args.morphometrics_csv is not None and Path(args.morphometrics_csv).is_file()) else next((p for p in (args.output_csv.parent / "morphometric_flags.csv", args.vouchers_csv.parent / "morphometric_flags.csv") if p.is_file()), None)
    df, _ = load_and_link_rosette_patches(args.rosette_dir, args.vouchers_csv, demo=args.demo, manifest_csv=args.manifest_csv, morphometrics_csv=morph_csv)
    if df.empty: logger.warning("No rosette patch records were matched. Check inputs."); return
    df["plant_individual_id"] = [_parse_pid(v, 0) for v in df["plant_individual_id"]] if "plant_individual_id" in df.columns else 0
    df = _apply_intraspecific_morphotype_labels(df, morph_csv, vouchers_csv=args.vouchers_csv, intraspecific_mode=args.intraspecific_mode)
    features, labels, key_order, fg_ratios = extract_dinov2_embeddings(df.to_dict("records"), neutralize_background=args.neutralize_background, demo=args.demo, return_fg_ratios=True, db_path=args.db)
    fg_by_key = {(k if isinstance(k, tuple) else (str(k).strip(), 0)): float(fg) for k, fg in zip(key_order, fg_ratios)}
    df["foreground_attention_ratio"] = [round(fg_by_key.get((str(c).strip(), _parse_pid(p, 0)), 0.85), 4) for c, p in zip(df["catalogNumber"], df["plant_individual_id"])]
    active_taxa, promo_s = sorted(t for t in set(df["taxon"]) if _is_valid_packera_taxon(t)) or ["Unassigned_Packera"], df.get("is_promoted_bronze", None)
    pred_probs, _, _ = compute_out_of_fold_probabilities(features, labels, n_splits=5, surface_orientations=df.get("surface_orientation", None), stratify_by_orientation=args.stratify_by_orientation, determiner_tiers=df.get("determiner_tier", None), spatial_balance_weights=df.get("spatial_balance_weight", None), institution_codes=df.get("institutionCode", None), is_promoted_bronze=promo_s)
    audit_df = run_confident_learning_audit(pred_probs=pred_probs, labels=labels, records_df=df, class_names=active_taxa, error_threshold=args.cleanlab_threshold, filter_abaxial_artifacts=args.filter_abaxial_artifacts, features=features, is_promoted_bronze=promo_s)
    if args.export_figures:
        flagged = audit_df[audit_df["is_label_corrupted"]].to_dict("records") or audit_df[audit_df["is_cleanlab_issue"]].to_dict("records") or audit_df.sort_values("c_error", ascending=False).head(8).to_dict("records")
        clf = LogisticRegression(max_iter=1000, C=0.25, class_weight="balanced", random_state=42).fit(features, labels) if len(np.unique(labels)) > 1 else None
        generate_gradcam_panel(flagged, args.output_figure, classifier=clf, demo=args.demo, class_names=active_taxa, active_taxa=active_taxa)
        fg_updates = {(str(r["catalogNumber"]).strip(), _parse_pid(r.get("plant_individual_id", 0), 0)): float(r["foreground_attention_ratio"]) for r in flagged if "catalogNumber" in r and "foreground_attention_ratio" in r}
        if fg_updates:
            audit_df["foreground_attention_ratio"] = [fg_updates.get((str(c).strip(), _parse_pid(p, 0)), float(fg)) for c, p, fg in zip(audit_df["catalogNumber"], audit_df["plant_individual_id"], audit_df["foreground_attention_ratio"])]
            audit_df["is_artifact_driven"] = (audit_df["foreground_attention_ratio"] < FOREGROUND_ATTENTION_MIN_RATIO).astype(bool); audit_df["is_background_shortcut"] = audit_df["is_artifact_driven"].astype(bool)
            if args.filter_abaxial_artifacts:
                ab_art = (audit_df["surface_orientation"].astype(str).str.lower() == "abaxial") & audit_df["is_artifact_driven"] & (audit_df["is_label_corrupted"] | audit_df["is_label_issue"])
                if bool(ab_art.any()): audit_df.loc[ab_art, ["is_label_corrupted", "is_label_issue", "triage_action"]] = [False, False, "Retain (Abaxial Inversion Artifact)"]
        try: generate_step05_publication_figures(audit_df, args.output_figure)
        except Exception as fig_exc: logger.warning(f"Step 05 publication figure generation skipped: {fig_exc}")
    if not args.demo:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True); tmp_csv = args.output_csv.with_suffix(args.output_csv.suffix + ".tmp")
        export_df = audit_df.drop(columns=["oof_probs", "rosette_patches", "rosette_masks", "tier1_leaf_patches", "tier1_leaf_masks", "tier1_leaf_weights"], errors="ignore")
        with open(tmp_csv, "w", encoding="utf-8", newline="") as f: export_df.to_csv(f, index=False); f.flush(); os.fsync(f.fileno())
        tmp_csv.replace(args.output_csv)
        if args.db is not None or is_canonical_or_explicit_db_context(args.output_csv):
            try:
                conn = get_db_connection(db_path=args.db, init_schema=True)
                try: persist_vision_xai_audit(conn, os.environ.get("RUN_ID", "canonical"), export_df.to_dict("records"), demo=args.demo)
                finally: conn.close()
            except Exception as db_exc: logger.warning(f"SQLite vision_xai_audit persistence warning: {db_exc}")
if __name__ == "__main__":
    main()
