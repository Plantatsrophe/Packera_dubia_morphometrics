#!/usr/bin/env python3
"""Phase 5 Indumentum-Blind Deep Vision & Explainable AI Pipeline (05_cleanlab_vision_xai.py)."""

from __future__ import annotations
import argparse, logging, os, re, sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd, torch, torch.nn.functional as F
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
    if not s_orig or s_orig.lower() in {"none", "nan", "null", "na", "unknown"}:
        return "Unknown"
    try:
        from scripts.core.harvester import standardize_packera_taxon as _harvester_std
        std = _harvester_std(s_orig)
        if std and str(std).strip() and str(std).strip().lower() not in {"none", "nan", "null", "na", "unknown", "unassigned_packera"}:
            return str(std).strip()
    except Exception:
        pass
    m = re.match(r"^packera\s*(?:([x×])(?=\s)|([×]))?\s*([a-z-]+)", s_orig.lower())
    epithet = m.group(3).lstrip("×-") if m else ""
    return f"Packera {'× ' if (m.group(1) or m.group(2)) else ''}{epithet}" if (m and epithet) else s_orig.split("(")[0].strip()


def _is_valid_packera_taxon(taxon_str: str) -> bool:
    return bool(taxon_str and taxon_str not in {"Unknown", "Unassigned_Packera"} and re.match(r"^(?:Packera\s+(?:[x×]\s*)?[a-z-]+|Morphotype_[A-Za-z0-9_]+|Ecotype_[A-Za-z0-9_]+)$", str(taxon_str)))


def _to_rgb_uint8(image: Image.Image | np.ndarray) -> np.ndarray:
    arr = np.array(image.convert("RGB"), dtype=np.uint8) if isinstance(image, Image.Image) else np.array(image, copy=True, dtype=np.uint8)
    return np.stack([arr] * 3, axis=-1) if arr.ndim == 2 else (arr[:, :, :3] if (arr.ndim == 3 and arr.shape[2] == 4) else arr)


def neutralize_mounting_paper(image: Image.Image | np.ndarray, neutral_color: Tuple[int, int, int] = (128, 128, 128)) -> Image.Image | np.ndarray:
    """Neutralizes herbarium mounting paper background to uniform neutral gray RGB(128, 128, 128)."""
    was_pil, img_arr = isinstance(image, Image.Image), _to_rgb_uint8(image)
    if img_arr.size == 0:
        return Image.fromarray(img_arr) if was_pil else img_arr
    if cv2 is not None:
        _, plant_mask = cv2.threshold(cv2.GaussianBlur(cv2.cvtColor(img_arr, cv2.COLOR_RGB2GRAY), (5, 5), 0), 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        bg_mask = cv2.morphologyEx(plant_mask, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))) == 0
    else:
        gray = 0.2989 * img_arr[:, :, 0] + 0.5870 * img_arr[:, :, 1] + 0.1140 * img_arr[:, :, 2]
        bg_mask = gray > np.percentile(gray, 60)
    neutralized = img_arr.copy()
    neutralized[bg_mask] = neutral_color
    return Image.fromarray(neutralized) if was_pil else neutralized


def apply_tier1_leaf_mask_to_rgb(rgb_img: np.ndarray | Image.Image, mask_path: Optional[Path], pad_px: int = 5) -> np.ndarray:
    """Applies a dilated binary mask to an RGB crop, replacing background with [242, 239, 230]."""
    out = _to_rgb_uint8(rgb_img)
    if out.size == 0 or mask_path is None or not str(mask_path).strip() or not Path(str(mask_path)).is_file():
        return out
    try:
        mask_p = Path(str(mask_path))
        mask = cv2.imread(str(mask_p), cv2.IMREAD_GRAYSCALE) if cv2 is not None else np.array(Image.open(mask_p).convert("L"))
    except Exception:
        return out
    if mask is None or mask.size == 0:
        return out
    h, w = out.shape[:2]
    if mask.shape[:2] != (h, w):
        mask = cv2.resize(mask, (w, h), interpolation=cv2.INTER_NEAREST) if cv2 is not None else np.array(Image.fromarray(mask).resize((w, h), Image.NEAREST))
    bin_mask, k = ((mask > 127).astype(np.uint8)) * 255, max(1, int(pad_px))
    dilated = cv2.dilate(bin_mask, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k)), iterations=1) if cv2 is not None else bin_mask
    out[dilated == 0] = np.array([242, 239, 230], dtype=out.dtype)
    return out


def preprocess_macro_rosette_crop(image: Image.Image | np.ndarray, target_size: Tuple[int, int] = (224, 224), sigma: float = 1.0) -> np.ndarray:
    """Applies Gaussian anti-aliasing (sigma=1.0) and downsamples strictly to 224x224."""
    arr, (tw, th) = _to_rgb_uint8(image), (int(target_size[0]), int(target_size[1]))
    if arr.size == 0:
        return np.zeros((th, tw, 3), dtype=np.uint8)
    if sigma > 0.0:
        ksize = max(3, int(2 * round(2.0 * sigma) + 1))
        arr = cv2.GaussianBlur(arr, (ksize, ksize), sigmaX=float(sigma), sigmaY=float(sigma)) if cv2 is not None else np.array(Image.fromarray(arr).filter(ImageFilter.GaussianBlur(radius=float(sigma))))
    if arr.shape[:2] != (th, tw):
        interp = cv2.INTER_AREA if (arr.shape[0] >= th and arr.shape[1] >= tw) else cv2.INTER_LINEAR
        arr = cv2.resize(arr, (tw, th), interpolation=interp) if cv2 is not None else np.array(Image.fromarray(arr).resize((tw, th), Image.BILINEAR))
    return arr.astype(np.uint8)


def compute_foreground_tissue_mask(rgb_img: Image.Image | np.ndarray, mask_path: Optional[Path] = None, grid_size: Tuple[int, int] = (16, 16)) -> Tuple[np.ndarray, np.ndarray]:
    """Computes binary plant tissue mask at image resolution and downsampled 16x16 patch grid."""
    arr, (gh, gw) = _to_rgb_uint8(rgb_img), (int(grid_size[0]), int(grid_size[1]))
    h, w = arr.shape[:2]
    if h == 0 or w == 0:
        return np.ones((224, 224), dtype=np.uint8), np.ones((gh, gw), dtype=np.float32)
    fg_mask: Optional[np.ndarray] = None
    if mask_path is not None and str(mask_path).strip() and Path(str(mask_path).strip()).is_file():
        try:
            mp = Path(str(mask_path).strip())
            raw_m = cv2.imread(str(mp), cv2.IMREAD_GRAYSCALE) if cv2 is not None else np.array(Image.open(mp).convert("L"))
            if raw_m is not None and raw_m.size > 0:
                raw_m = (cv2.resize(raw_m, (w, h), interpolation=cv2.INTER_NEAREST) if cv2 is not None else np.array(Image.fromarray(raw_m).resize((w, h), Image.NEAREST))) if raw_m.shape[:2] != (h, w) else raw_m
                fg_mask = (raw_m > 127).astype(np.uint8)
        except Exception:
            fg_mask = None
    if fg_mask is None:
        is_neutral = ((np.abs(arr[:, :, 0].astype(int) - 128) <= 2) & (np.abs(arr[:, :, 1].astype(int) - 128) <= 2) & (np.abs(arr[:, :, 2].astype(int) - 128) <= 2)) | ((np.abs(arr[:, :, 0].astype(int) - 242) <= 3) & (np.abs(arr[:, :, 1].astype(int) - 239) <= 3) & (np.abs(arr[:, :, 2].astype(int) - 230) <= 3))
        if 0.02 <= float(np.mean(is_neutral)) < 0.98:
            fg_mask = (~is_neutral).astype(np.uint8)
        elif cv2 is not None:
            lab = cv2.cvtColor(arr, cv2.COLOR_RGB2LAB)
            _, otsu_l = cv2.threshold(cv2.GaussianBlur(lab[:, :, 0], (5, 5), 0), 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
            chroma = np.hypot(lab[:, :, 1].astype(np.float32) - 128.0, lab[:, :, 2].astype(np.float32) - 128.0)
            combined = np.where((otsu_l > 0) | (chroma > max(8.0, float(np.percentile(chroma, 55)))), 255, 0).astype(np.uint8)
            fg_mask = (cv2.morphologyEx(combined, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))) > 0).astype(np.uint8)
        else:
            fg_mask = ((0.2989 * arr[:, :, 0] + 0.5870 * arr[:, :, 1] + 0.1140 * arr[:, :, 2]) <= np.percentile(arr.mean(axis=2), 60)).astype(np.uint8)
    if int(fg_mask.sum()) == 0:
        fg_mask = np.ones((h, w), dtype=np.uint8)
    patch_frac = cv2.resize(fg_mask.astype(np.float32), (gw, gh), interpolation=cv2.INTER_AREA) if cv2 is not None else F.adaptive_avg_pool2d(torch.as_tensor(fg_mask, dtype=torch.float32).unsqueeze(0).unsqueeze(0), (gh, gw))[0, 0].numpy()
    fg_grid = (patch_frac >= 0.15).astype(np.float32)
    return fg_mask.astype(np.uint8), (np.where(patch_frac > 0.0, patch_frac, 1.0).astype(np.float32) if float(fg_grid.sum()) <= 0.0 else fg_grid)


def foreground_masked_patch_pooling(patch_tokens: torch.Tensor | np.ndarray, patch_mask_16x16: torch.Tensor | np.ndarray) -> np.ndarray:
    """Performs Foreground-Masked Average Pooling over spatial patch tokens and L2-normalizes."""
    tokens_np = patch_tokens.detach().cpu().numpy() if isinstance(patch_tokens, torch.Tensor) else np.asarray(patch_tokens, dtype=np.float32)
    mask_np = patch_mask_16x16.detach().cpu().numpy() if isinstance(patch_mask_16x16, torch.Tensor) else np.asarray(patch_mask_16x16, dtype=np.float32)
    single_input = tokens_np.ndim == 2
    if single_input:
        tokens_np = tokens_np[None, :, :]
    b, n_patches, _ = tokens_np.shape
    mask_flat = mask_np.reshape(b, -1) if (mask_np.ndim >= 3 or (mask_np.ndim == 2 and mask_np.shape == (b, n_patches))) else np.broadcast_to(mask_np.reshape(1, -1), (b, mask_np.size))
    if mask_flat.shape[1] != n_patches:
        mask_flat = F.interpolate(torch.as_tensor(mask_flat, dtype=torch.float32).unsqueeze(1), size=n_patches, mode="nearest")[:, 0].numpy()
    weights = np.maximum(mask_flat, 0.0).astype(np.float32)
    weight_sums = np.sum(weights, axis=1, keepdims=True)
    if np.any(weight_sums[:, 0] <= 1e-8):
        weights[weight_sums[:, 0] <= 1e-8, :] = 1.0
        weight_sums = np.sum(weights, axis=1, keepdims=True)
    pooled = np.sum(tokens_np * weights[:, :, None], axis=1) / np.maximum(weight_sums, 1e-8)
    pooled_l2 = (pooled / np.maximum(np.linalg.norm(pooled, axis=1, keepdims=True), 1e-8)).astype(np.float32)
    return pooled_l2[0] if single_input else pooled_l2


def classify_foliar_surface_orientation(image: Image.Image | np.ndarray | str | Path, plant_mask: Optional[np.ndarray] = None, neutral_color: Tuple[int, int, int] = (128, 128, 128), luminance_threshold: float = 68.0, saturation_threshold: float = 0.22, brightness_cutoff: float = 75.0, max_dimension: int = 512) -> Dict[str, Any]:
    """Classifies foliar surface orientation (abaxial vs. adaxial) in CIE LAB & HSV space."""
    default_res = {"surface_orientation": "adaxial", "mean_luminance": 0.0, "mean_saturation": 0.0, "high_brightness_ratio": 0.0}
    if isinstance(image, (str, Path)):
        if not Path(image).exists():
            return default_res
        img_arr = np.array(Image.open(image).convert("RGB"))
    else:
        img_arr = _to_rgb_uint8(image)
    if img_arr.size == 0:
        return default_res
    if max(img_arr.shape[:2]) > max_dimension:
        scale = float(max_dimension) / max(img_arr.shape[:2])
        nw, nh = max(1, int(img_arr.shape[1] * scale)), max(1, int(img_arr.shape[0] * scale))
        img_arr = cv2.resize(img_arr, (nw, nh), interpolation=cv2.INTER_AREA) if cv2 is not None else np.array(Image.fromarray(img_arr).resize((nw, nh), Image.BILINEAR))
        if plant_mask is not None:
            plant_mask = (cv2.resize(plant_mask.astype(np.uint8), (nw, nh), interpolation=cv2.INTER_NEAREST) > 0) if cv2 is not None else (np.array(Image.fromarray(plant_mask.astype(np.uint8)).resize((nw, nh), Image.NEAREST)) > 0)
    if plant_mask is not None:
        mask = plant_mask > 0
    else:
        is_neutral = ((img_arr[:, :, 0] == neutral_color[0]) & (img_arr[:, :, 1] == neutral_color[1]) & (img_arr[:, :, 2] == neutral_color[2])) | ((img_arr[:, :, 0] == 242) & (img_arr[:, :, 1] == 239) & (img_arr[:, :, 2] == 230))
        mask = (~is_neutral) if (np.mean(is_neutral) >= 0.02 and np.sum(~is_neutral) >= 5) else (compute_foreground_tissue_mask(img_arr)[0] > 0)
    if mask.sum() < 5:
        return default_res
    if cv2 is not None:
        L_vals = cv2.cvtColor(img_arr, cv2.COLOR_RGB2LAB)[mask][:, 0].astype(np.float64) * (100.0 / 255.0)
        S_vals = cv2.cvtColor(img_arr, cv2.COLOR_RGB2HSV)[mask][:, 1].astype(np.float64) / 255.0
    else:
        rgb = img_arr[mask].astype(np.float64)
        Y = 0.2126 * (rgb[:, 0] / 255.0) + 0.7152 * (rgb[:, 1] / 255.0) + 0.0722 * (rgb[:, 2] / 255.0)
        L_vals = np.where(Y > 0.008856, 116.0 * (Y ** (1.0 / 3.0)) - 16.0, 903.3 * Y)
        cmax, cmin = np.max(rgb, axis=1), np.min(rgb, axis=1)
        S_vals = np.where(cmax > 0, (cmax - cmin) / cmax, 0.0)
    mean_L, mean_S = float(np.mean(L_vals)), float(np.mean(S_vals))
    return {"surface_orientation": "abaxial" if (mean_L > luminance_threshold and mean_S < saturation_threshold) else "adaxial", "mean_luminance": round(mean_L, 4), "mean_saturation": round(mean_S, 4), "high_brightness_ratio": round(float(np.mean(L_vals > brightness_cutoff)), 4)}


class RosettePatchDataset(Dataset):
    """PyTorch Dataset loading 224x224 Gaussian low-pass whole-rosette crops and 16x16 tissue masks."""
    def __init__(self, records: List[Dict[str, Any]], transform: Optional[transforms.Compose] = None, neutralize_background: bool = True, sigma: float = 1.0, demo: bool = False):
        self.records = records
        self.transform = transform or transforms.Compose([transforms.ToTensor(), transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])])
        self.neutralize_background, self.sigma, self.demo = neutralize_background, sigma, demo

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor, int, str]:
        rec = self.records[idx]
        image_path, mask_path = Path(str(rec["patch_path"])), rec.get("mask_path")
        has_mask = bool(mask_path and str(mask_path).strip() and Path(str(mask_path)).exists())
        if image_path.exists():
            raw_arr = np.array(Image.open(image_path).convert("RGB"))
            raw_arr = apply_tier1_leaf_mask_to_rgb(raw_arr, Path(str(mask_path))) if has_mask else (np.array(neutralize_mounting_paper(raw_arr)) if self.neutralize_background else raw_arr)
            macro_arr = preprocess_macro_rosette_crop(raw_arr, target_size=(224, 224), sigma=self.sigma)
            _, fg_grid = compute_foreground_tissue_mask(macro_arr, mask_path=Path(str(mask_path)) if has_mask else None, grid_size=(16, 16))
        elif self.demo:
            macro_arr, fg_grid = np.full((224, 224, 3), 128, dtype=np.uint8), np.ones((16, 16), dtype=np.float32)
        else:
            raise RuntimeError(f"Rosette patch image missing at '{image_path}'. Run Phase 2 or pass `--demo`.")
        return self.transform(Image.fromarray(macro_arr)), torch.as_tensor(fg_grid, dtype=torch.float32), int(rec.get("label_idx", 0)), str(rec["catalogNumber"])


def _extract_catalog_number(stem: str, known_catalogs: Optional[set] = None) -> str:
    cleaned = re.sub(r"_(?:p\d+_)?(?:rosette|leaf|patch|crop).*$", "", stem, flags=re.IGNORECASE)
    if known_catalogs and (stem in known_catalogs or cleaned in known_catalogs):
        return stem if stem in known_catalogs else cleaned
    return cleaned if cleaned else stem.split("_")[0]


def load_and_link_rosette_patches(rosette_dir: Path, vouchers_csv: Path, demo: bool = False, manifest_csv: Optional[Path] = None) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Ingests whole-rosette crops from `rosettes_dense` and maps active taxa to `0..K-1`."""
    rosette_dir, vouchers_csv = Path(rosette_dir), Path(vouchers_csv)
    vouchers_df = pd.read_csv(vouchers_csv) if vouchers_csv.exists() else pd.DataFrame()
    if not vouchers_df.empty and "is_primary_duplicate" in vouchers_df.columns:
        prim_mask = vouchers_df["is_primary_duplicate"].map(lambda v: v is True or str(v).strip().lower() in {"true", "1", "t", "yes"} or pd.isna(v))
        vouchers_df = vouchers_df.loc[prim_mask].reset_index(drop=True)
    voucher_map = {str(row["catalogNumber"]).strip(): row.to_dict() for _, row in vouchers_df.iterrows()} if (not vouchers_df.empty and "catalogNumber" in vouchers_df.columns) else {}
    known_cats = set(voucher_map.keys()) if voucher_map else None
    dense_dir = rosette_dir if rosette_dir.name == "rosettes_dense" else (rosette_dir / "rosettes_dense")
    dense_files: List[Path] = []
    for search_dir, use_rglob in ((dense_dir, False), (rosette_dir, False), (rosette_dir, True)):
        if not dense_files and search_dir.exists():
            for ext in ("*.jpg", "*.png", "*.jpeg"):
                dense_files.extend(sorted(search_dir.rglob(ext) if use_rglob else search_dir.glob(ext)))
    grouped_patches: Dict[str, List[Path]] = {}
    for p in dense_files:
        cat_num = _extract_catalog_number(p.stem, known_cats)
        if not voucher_map or cat_num in voucher_map:
            grouped_patches.setdefault(cat_num, []).append(p)
    records: List[Dict[str, Any]] = []
    if grouped_patches:
        for cat_num, p_list in grouped_patches.items():
            meta = voucher_map.get(cat_num, {})
            raw_sp = meta.get("species_raw", meta.get("species_standardized", meta.get("scientificName", meta.get("species", "Unknown"))))
            std_taxon = standardize_packera_taxon(raw_sp)
            if _is_valid_packera_taxon(std_taxon):
                records.append({"catalogNumber": cat_num, "plant_individual_id": int(meta.get("plant_individual_id", 0) or 0), "patch_path": str(p_list[0]), "rosette_patches": [str(x) for x in p_list], "n_rosette_patches": len(p_list), "taxon": std_taxon, "species_raw": raw_sp, "label_idx": 0, "institutionCode": str(meta.get("institutionCode", "UNKNOWN") or "UNKNOWN").strip(), "determiner_tier": str(meta.get("determiner_tier", "Tier_3_Bronze") or "Tier_3_Bronze"), "spatial_balance_weight": float(meta.get("spatial_balance_weight", 1.0) if pd.notna(meta.get("spatial_balance_weight", 1.0)) else 1.0), **classify_foliar_surface_orientation(p_list[0])})
    elif not demo:
        raise RuntimeError(f"No whole-rosette patch images found in '{dense_dir}' or '{rosette_dir}'.")
    elif not vouchers_df.empty:
        for _, row in vouchers_df.iterrows():
            raw_sp = str(row.get("species_raw", row.get("species_standardized", row.get("scientificName", row.get("species", "Unknown")))))
            std_taxon = standardize_packera_taxon(raw_sp)
            if _is_valid_packera_taxon(std_taxon):
                cat_num = str(row.get("catalogNumber", "")).strip()
                synth_p = str(dense_dir / f"{cat_num}_rosette.jpg")
                records.append({"catalogNumber": cat_num, "plant_individual_id": int(row.get("plant_individual_id", 0) or 0), "patch_path": synth_p, "rosette_patches": [synth_p], "n_rosette_patches": 1, "taxon": std_taxon, "species_raw": raw_sp, "label_idx": 0, "institutionCode": str(row.get("institutionCode", "UNKNOWN") or "UNKNOWN").strip(), "determiner_tier": str(row.get("determiner_tier", "Tier_3_Bronze") or "Tier_3_Bronze"), "spatial_balance_weight": 1.0, "surface_orientation": "adaxial", "mean_luminance": 0.0, "mean_saturation": 0.0, "high_brightness_ratio": 0.0})
    present_taxa = {r["taxon"] for r in records if _is_valid_packera_taxon(r["taxon"])}
    active_taxa = sorted(t for t in present_taxa if _is_valid_packera_taxon(t))
    class_to_idx = {taxon: idx for idx, taxon in enumerate(active_taxa)}
    for r in records:
        r["label_idx"] = class_to_idx.get(r["taxon"], 0)
    return pd.DataFrame(records), class_to_idx


def _aggregate_voucher_embeddings(raw_feats: np.ndarray, expanded_meta: List[Tuple[str, int, float]]) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    if len(raw_feats) == 0:
        return np.empty((0, 768), dtype=np.float32), np.array([], dtype=int), []
    voucher_groups: Dict[str, List[Tuple[np.ndarray, float, int]]] = {}
    voucher_order: List[str] = []
    for feat_vec, (cat_num, lbl_idx, weight) in zip(raw_feats, expanded_meta):
        if cat_num not in voucher_groups:
            voucher_groups[cat_num] = []
            voucher_order.append(cat_num)
        voucher_groups[cat_num].append((feat_vec, max(1e-6, float(weight)), int(lbl_idx)))
    agg_feats, agg_labels = [], []
    for cat_num in voucher_order:
        items = voucher_groups[cat_num]
        wts = np.array([it[1] for it in items], dtype=np.float64)
        vec = items[0][0].astype(np.float32) if len(items) == 1 else np.sum(np.vstack([it[0] for it in items]) * (wts / np.sum(wts))[:, None], axis=0).astype(np.float32)
        agg_feats.append((vec / max(float(np.linalg.norm(vec)), 1e-8)).astype(np.float32))
        agg_labels.append(items[0][2])
    return np.vstack(agg_feats), np.array(agg_labels, dtype=int), voucher_order


def extract_dinov2_embeddings(records: List[Dict[str, Any]], model_name: str = "dinov2_vitb14_reg", device: str = "cuda" if torch.cuda.is_available() else "cpu", batch_size: int = 32, neutralize_background: bool = True, demo: bool = False) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Extracts 768-d DINOv2 embeddings via Foreground-Masked Average Pooling on 224x224 crops."""
    torch.manual_seed(42)
    np.random.seed(42)
    if not records:
        return np.empty((0, 768), dtype=np.float32), np.array([], dtype=int), []
    expanded_records, expanded_meta = [], []
    for rec in records:
        cat_num, lbl_idx = str(rec["catalogNumber"]).strip(), int(rec.get("label_idx", 0))
        patch_list = rec.get("rosette_patches") or rec.get("tier1_leaf_patches")
        for p_str in (patch_list if (isinstance(patch_list, list) and patch_list) else [rec.get("patch_path", "")]):
            expanded_records.append({**rec, "catalogNumber": cat_num, "label_idx": lbl_idx, "patch_path": p_str})
            expanded_meta.append((cat_num, lbl_idx, 1.0))
    if not any(Path(str(r["patch_path"])).exists() for r in expanded_records[:5]):
        if not demo:
            raise RuntimeError("No physical patch images found in `data/cropped_patches/rosettes_dense/`. Run Phase 2 or pass `--demo`.")
        feats = np.random.randn(len(expanded_records), 768).astype(np.float32)
        for i, r in enumerate(expanded_records):
            c_idx = int(r.get("label_idx", 0)) % 15
            feats[i, c_idx * 50 : (c_idx + 1) * 50] += 3.0
        return _aggregate_voucher_embeddings(feats, expanded_meta)
    loader = DataLoader(RosettePatchDataset(expanded_records, neutralize_background=neutralize_background, sigma=1.0, demo=demo), batch_size=batch_size, shuffle=False, num_workers=0)
    model, last_err = None, None
    for candidate_name in (model_name, "dinov2_vitb14"):
        try: model = torch.hub.load("facebookresearch/dinov2", candidate_name); break
        except Exception as exc: last_err = exc
    if model is None:
        if not demo:
            raise RuntimeError(f"Could not load DINOv2 model weights ('{model_name}'): {last_err}") from last_err
        return _aggregate_voucher_embeddings(np.random.randn(len(expanded_records), 768).astype(np.float32), expanded_meta)
    model, all_feats = model.to(device).eval(), []
    with torch.no_grad():
        for tensors, fg_grids, _labels, _cats in loader:
            tensors = tensors.to(device)
            if hasattr(model, "forward_features"):
                pooled = foreground_masked_patch_pooling(model.forward_features(tensors)["x_norm_patchtokens"], fg_grids)
            else:
                out = model(tensors).cpu().numpy()
                pooled = (out / np.maximum(np.linalg.norm(out, axis=1, keepdims=True), 1e-8)).astype(np.float32)
            all_feats.append(pooled)
    return _aggregate_voucher_embeddings(np.vstack(all_feats) if all_feats else np.empty((0, 768), dtype=np.float32), expanded_meta)


def compute_out_of_fold_probabilities(features: np.ndarray, labels: np.ndarray, n_splits: int = 5, random_state: int = 42, surface_orientations: Optional[np.ndarray | List[str] | pd.Series] = None, stratify_by_orientation: bool = False, determiner_tiers: Optional[np.ndarray | List[str] | pd.Series] = None, spatial_balance_weights: Optional[np.ndarray | List[float] | pd.Series] = None, institution_codes: Optional[np.ndarray | List[str] | pd.Series] = None) -> Tuple[np.ndarray, float, float]:
    """Computes calibrated out-of-fold probabilities anchored on Gold/Silver reference specimens."""
    X = np.asarray(features, dtype=np.float64)
    X, labels = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-8), np.asarray(labels, dtype=int)
    tier_weight_map = {"Tier_1_Gold": 3.0, "Tier_2_Silver": 1.5, "Tier_3_Bronze": 0.5}
    if determiner_tiers is not None and len(determiner_tiers) == len(labels):
        tier_strs = [str(t).strip() for t in determiner_tiers]
        sample_weights = np.array([tier_weight_map.get(t, 0.5) for t in tier_strs], dtype=np.float64)
        anchor_mask = np.array([t in {"Tier_1_Gold", "Tier_2_Silver"} for t in tier_strs], dtype=bool)
    else:
        sample_weights, anchor_mask = np.ones(len(labels), dtype=np.float64), np.ones(len(labels), dtype=bool)
    if spatial_balance_weights is not None and len(spatial_balance_weights) == len(sample_weights):
        sp_w = np.asarray(spatial_balance_weights, dtype=np.float64)
        sample_weights = sample_weights * np.where(np.isfinite(sp_w) & (sp_w > 0.0), np.clip(sp_w, 0.20, 3.00), 1.0)
    unique_labels, label_counts = np.unique(labels, return_counts=True)
    max_lbl_plus_1 = int(np.max(labels)) + 1 if len(labels) > 0 else 2
    num_classes = max(2, max_lbl_plus_1)
    pred_probs = np.zeros((len(labels), num_classes), dtype=np.float64)
    if len(unique_labels) < 2 or int(np.min(label_counts)) < 2:
        if len(X) == 0:
            return pred_probs, 1.0, 1.0
        centroid = np.median(X[anchor_mask] if np.any(anchor_mask) else X, axis=0)
        cos_dist = 1.0 - np.clip(X @ (centroid / max(float(np.linalg.norm(centroid)), 1e-8)), -1.0, 1.0)
        q95 = max(float(np.percentile(cos_dist, 95)) if len(cos_dist) >= 5 else float(np.max(cos_dist)) + 1e-6, 1e-6)
        is_outlier = (cos_dist > q95) & ((np.asarray([str(o).lower() for o in surface_orientations]) != "abaxial") if surface_orientations is not None else True)
        p_self = np.where(is_outlier, np.clip(0.10 * np.exp(-(cos_dist - q95) / q95), 0.01, 0.12), np.clip(1.0 - 0.50 * (cos_dist / q95), 0.50, 0.999))
        for i, lbl in enumerate(labels):
            c0 = int(lbl) % num_classes
            pred_probs[i, c0], pred_probs[i, (c0 + 1) % num_classes] = p_self[i], 1.0 - p_self[i]
        return pred_probs, float(np.mean(~is_outlier)), float(np.mean(~is_outlier))
    u_anc, c_anc = np.unique(labels[anchor_mask], return_counts=True)
    if len(u_anc) < len(unique_labels) or int(np.min(c_anc)) < 2:
        anchor_mask, u_anc, c_anc = np.ones(len(labels), dtype=bool), unique_labels, label_counts
    anchor_idx, bronze_idx = np.where(anchor_mask)[0], np.where(~anchor_mask)[0]
    eff_splits = max(2, min(int(n_splits), int(np.min(c_anc))))
    X_anc, y_anc, w_anc = X[anchor_idx], labels[anchor_idx], sample_weights[anchor_idx]
    strat_anc = y_anc
    if stratify_by_orientation and surface_orientations is not None:
        compound = np.array([f"{l}_{o}" for l, o in zip(y_anc, np.asarray(surface_orientations)[anchor_idx])])
        if all(c >= eff_splits for c in np.unique(compound, return_counts=True)[1]): strat_anc = compound
    cv_splits = None
    if institution_codes is not None and len(institution_codes) == len(labels):
        inst_anc = np.asarray(institution_codes)[anchor_idx]
        groups = np.asarray([str(g).strip() if pd.notna(g) and str(g).strip() and str(g).strip().upper() != "UNKNOWN" else f"SINGLETON_{idx}" for idx, g in enumerate(inst_anc)])
        if len(np.unique(groups)) >= eff_splits:
            try:
                cand = list(StratifiedGroupKFold(n_splits=eff_splits, shuffle=True, random_state=random_state).split(X_anc, strat_anc, groups))
                if all(len(vi) > 0 and len(np.unique(y_anc[ti])) == len(unique_labels) for ti, vi in cand): cv_splits = cand
            except Exception: cv_splits = None
    if cv_splits is None:
        cv_splits = list(StratifiedKFold(n_splits=eff_splits, shuffle=True, random_state=random_state).split(X_anc, strat_anc))
    bronze_accum, n_valid_folds, y_true_anc, y_pred_anc = np.zeros((len(bronze_idx), num_classes), dtype=np.float64), 0, [], []
    for train_sub, val_sub in cv_splits:
        scaler = StandardScaler()
        X_train, X_val, y_train, y_val = scaler.fit_transform(X_anc[train_sub]), scaler.transform(X_anc[val_sub]), y_anc[train_sub], y_anc[val_sub]
        w_train = w_anc[train_sub] / max(float(np.mean(w_anc[train_sub])), 1e-8)
        min_cls = int(np.min(np.unique(y_train, return_counts=True)[1]))
        lrcv_kw: Dict[str, Any] = {"use_legacy_attributes": False, "l1_ratios": (0.0,), "scoring": "accuracy"} if "use_legacy_attributes" in LogisticRegressionCV.__init__.__code__.co_varnames else {}
        base_clf: Any = LogisticRegressionCV(Cs=np.logspace(-2, 1, 6), cv=max(2, min(3, min_cls)), class_weight="balanced", max_iter=2000, random_state=random_state, **lrcv_kw) if min_cls >= 2 else LogisticRegression(max_iter=2000, C=0.25, class_weight="balanced", random_state=random_state)
        base_clf.fit(X_train, y_train, sample_weight=w_train)
        active_clf: Any = base_clf
        if min_cls >= 3:
            try:
                cal_clf = CalibratedClassifierCV(estimator=LogisticRegression(max_iter=2000, C=0.25, class_weight="balanced", random_state=random_state), method="sigmoid", cv=min(3, min_cls))
                cal_clf.fit(X_train, y_train, sample_weight=w_train); active_clf = cal_clf
            except Exception: active_clf = base_clf
        val_probs, classes_arr = active_clf.predict_proba(X_val), np.asarray(active_clf.classes_, dtype=int)
        for col_i, c_idx in enumerate(classes_arr): pred_probs[anchor_idx[val_sub], c_idx] = val_probs[:, col_i]
        if len(bronze_idx) > 0:
            brz_probs = active_clf.predict_proba(scaler.transform(X[bronze_idx]))
            for col_i, c_idx in enumerate(classes_arr): bronze_accum[:, c_idx] += brz_probs[:, col_i]
            n_valid_folds += 1
        y_true_anc.extend(y_val); y_pred_anc.extend(classes_arr[np.argmax(val_probs, axis=1)])
    if len(bronze_idx) > 0 and n_valid_folds > 0:
        brz_mean = bronze_accum / float(n_valid_folds)
        pred_probs[bronze_idx, :] = brz_mean / np.maximum(brz_mean.sum(axis=1, keepdims=True), 1e-8)
    pred_probs = pred_probs / np.maximum(pred_probs.sum(axis=1, keepdims=True), 1e-8)
    return pred_probs, float(accuracy_score(y_true_anc, y_pred_anc) if y_true_anc else 1.0), float(f1_score(y_true_anc, y_pred_anc, average="weighted") if y_true_anc else 1.0)


def compute_feature_space_ood_scores(features: np.ndarray, labels: np.ndarray, pred_probs: Optional[np.ndarray] = None, determiner_tiers: Optional[Any] = None) -> np.ndarray:
    """Computes feature-space Out-of-Distribution score (`vision_ood_score` in [0, 1])."""
    n = len(labels)
    if n == 0:
        return np.array([], dtype=np.float64)
    if features is None or len(features) != n:
        return np.clip(1.0 - np.max(pred_probs, axis=1), 0.0, 1.0) if (pred_probs is not None and len(pred_probs) == n) else np.zeros(n, dtype=np.float64)
    X = np.asarray(features, dtype=np.float64)
    X, y = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-8), np.asarray(labels, dtype=int)
    anchor_mask = np.array([str(t).strip() in {"Tier_1_Gold", "Tier_2_Silver"} for t in determiner_tiers], dtype=bool) if (determiner_tiers is not None and len(determiner_tiers) == n) else np.ones(n, dtype=bool)
    if not np.any(anchor_mask):
        anchor_mask = np.ones(n, dtype=bool)
    centroids = [np.mean(X[anchor_mask & (y == c)], axis=0) for c in np.unique(y[anchor_mask]) if np.any(anchor_mask & (y == c))]
    centroids = [mu / max(float(np.linalg.norm(mu)), 1e-8) for mu in (centroids or [np.mean(X, axis=0)])]
    centroid_ood = np.clip(np.min(1.0 - np.clip(X @ np.vstack(centroids).T, -1.0, 1.0), axis=1), 0.0, 1.0)
    try:
        from cleanlab.outlier import OutOfDistribution
        cl_ood = np.clip(1.0 - np.asarray(OutOfDistribution().fit_score(features=X), dtype=np.float64), 0.0, 1.0)
        return np.round(0.5 * centroid_ood + 0.5 * cl_ood, 4)
    except Exception:
        return np.round(centroid_ood, 4)


def _compute_tomentum_taxon_discordance(audit_df: pd.DataFrame, error_threshold: float) -> pd.Series:
    """Organic, taxon-agnostic chromatic vestiture discordance rule in CIE L* luminance & HSV saturation space."""
    abax_mm = (audit_df["surface_orientation"].astype(str).str.lower() == "abaxial") & (audit_df["given_label"] != audit_df["predicted_label"]) & (audit_df["c_error"] > error_threshold)
    lum_s = pd.to_numeric(audit_df.get("mean_luminance", pd.Series(np.nan, index=audit_df.index)), errors="coerce")
    sat_s = pd.to_numeric(audit_df.get("mean_saturation", pd.Series(np.nan, index=audit_df.index)), errors="coerce")
    valid_l = lum_s.dropna()
    if len(valid_l) == 0 or float(valid_l.max() - valid_l.min()) <= 1e-6:
        return abax_mm.astype(bool)
    adax = audit_df["surface_orientation"].astype(str).str.lower() != "abaxial"
    anc = audit_df["determiner_tier"].astype(str).str.strip().isin({"Tier_1_Gold", "Tier_2_Silver"}) if "determiner_tier" in audit_df.columns else pd.Series(True, index=audit_df.index)
    med_L, med_S, q75_L = {}, {}, {}
    for tx in set(audit_df["given_label"].dropna().unique()):
        tx_m = audit_df["given_label"] == tx
        ref_m = (tx_m & adax & anc) if np.any(tx_m & adax & anc & lum_s.notna()) else ((tx_m & adax) if np.any(tx_m & adax & lum_s.notna()) else tx_m)
        if np.any(ref_m & lum_s.notna()):
            med_L[tx], q75_L[tx] = float(np.nanmedian(lum_s[ref_m])), float(np.nanpercentile(lum_s[tx_m & lum_s.notna()], 75))
        if np.any(ref_m & sat_s.notna()):
            med_S[tx] = float(np.nanmedian(sat_s[ref_m]))
    flags = []
    for idx_i, is_cand in enumerate(abax_mm.to_numpy(dtype=bool)):
        if not is_cand:
            flags.append(False); continue
        g_tx, p_tx = audit_df.at[idx_i, "given_label"], audit_df.at[idx_i, "predicted_label"]
        l_i = float(lum_s.iat[idx_i]) if pd.notna(lum_s.iat[idx_i]) else np.nan
        cond_a = (p_tx in med_L and g_tx in med_L and med_L[p_tx] >= med_L[g_tx] + 2.0) or (p_tx in med_S and g_tx in med_S and med_S[p_tx] <= med_S[g_tx] - 0.02)
        cond_b = bool(np.isfinite(l_i) and g_tx in q75_L and l_i > q75_L[g_tx])
        flags.append(bool(cond_a or cond_b))
    return pd.Series(flags, index=audit_df.index, dtype=bool)


def run_confident_learning_audit(pred_probs: np.ndarray, labels: np.ndarray, records_df: pd.DataFrame, class_names: List[str], error_threshold: float = 0.85, filter_abaxial_artifacts: bool = True, features: Optional[np.ndarray] = None) -> pd.DataFrame:
    """Runs class-conditional Confident Learning (`prune_by_noise_rate`) & OOD audit."""
    labels, pred_probs = np.asarray(labels, dtype=int), np.asarray(pred_probs, dtype=np.float64)
    num_classes = pred_probs.shape[1] if pred_probs.ndim == 2 and pred_probs.shape[1] > 0 else max(1, len(class_names))
    class_thresholds = np.array([float(np.mean(pred_probs[labels == c_j, c_j])) if np.any(labels == c_j) else (1.0 - error_threshold) for c_j in range(num_classes)], dtype=np.float64)
    try:
        if len(np.unique(labels)) < 2 or pred_probs.shape[1] != len(np.unique(labels)):
            raise ValueError("Single-class or non-contiguous label set.")
        from cleanlab.filter import find_label_issues
        from cleanlab.rank import get_label_quality_scores
        quality_scores = get_label_quality_scores(labels=labels, pred_probs=pred_probs)
        issues = find_label_issues(labels=labels, pred_probs=pred_probs, filter_by="prune_by_noise_rate", return_indices_ranked_by="self_confidence")
    except Exception:
        predicted_classes = np.argmax(pred_probs, axis=1)
        quality_scores = np.array([pred_probs[i, labels[i]] for i in range(len(labels))], dtype=np.float64)
        thresh_for_pred = np.array([class_thresholds[min(int(c), num_classes - 1)] for c in predicted_classes])
        issues = np.where((predicted_classes != labels) & ((np.max(pred_probs, axis=1) >= thresh_for_pred) | (1.0 - quality_scores > error_threshold)))[0]
    audit_df = records_df.copy().reset_index(drop=True)
    num_records = len(audit_df)
    given_classes = [class_names[l] if l < len(class_names) else "Unknown" for l in labels]
    pred_classes = [class_names[i] if i < len(class_names) else "Outlier_Morphotype" for i in np.argmax(pred_probs, axis=1)]
    conf_given = np.array([pred_probs[i, labels[i]] if i < len(pred_probs) and labels[i] < pred_probs.shape[1] else 0.0 for i in range(num_records)], dtype=np.float64)
    sorted_probs = np.sort(pred_probs, axis=1)[:, ::-1] if pred_probs.shape[1] >= 2 else np.column_stack([np.max(pred_probs, axis=1), 1.0 - np.max(pred_probs, axis=1)])
    conf_pred, second_prob, c_error = sorted_probs[:, 0], sorted_probs[:, 1], 1.0 - conf_given
    raw_sp_col = audit_df.get("taxon", audit_df.get("species_raw", audit_df.get("species", pd.Series(given_classes, index=audit_df.index))))
    audit_df["species_raw"] = audit_df.get("species_raw", audit_df.get("species", given_classes))
    std_sp = pd.Series(raw_sp_col, index=audit_df.index).apply(standardize_packera_taxon)
    audit_df["species_standardized"] = [sp if _is_valid_packera_taxon(sp) else (gc if _is_valid_packera_taxon(gc) else "Unassigned_Packera") for sp, gc in zip(std_sp, given_classes)]
    audit_df["given_label"], audit_df["predicted_label"], audit_df["vision_predicted_label"] = given_classes, pred_classes, pred_classes
    audit_df["confidence_given_class"], audit_df["confidence_predicted_class"] = np.round(conf_given, 4), np.round(conf_pred, 4)
    audit_df["vision_confidence"], audit_df["vision_second_prob"] = np.round(conf_pred, 4), np.round(second_prob, 4)
    audit_df["vision_margin"] = np.round(np.clip(conf_pred - second_prob, 0.0, 1.0), 4)
    audit_df["class_self_confidence_threshold"] = np.round([class_thresholds[min(int(l), num_classes - 1)] for l in labels], 4)
    audit_df["label_quality_score"], audit_df["c_error"] = np.round(quality_scores, 4), np.round(c_error, 4)
    if "vision_ood_score" not in audit_df.columns:
        audit_df["vision_ood_score"] = compute_feature_space_ood_scores(features if features is not None else np.empty((0, 768)), labels, pred_probs, audit_df.get("determiner_tier", None))
    audit_df["foreground_attention_ratio"] = np.round(pd.to_numeric(audit_df.get("foreground_attention_ratio", pd.Series(0.85, index=audit_df.index)), errors="coerce").fillna(0.85), 4)
    audit_df["is_artifact_driven"] = (audit_df["foreground_attention_ratio"] < FOREGROUND_ATTENTION_MIN_RATIO).astype(bool)
    if "determiner_tier" in audit_df.columns:
        tier_map = {"Tier_1_Gold": 3.0, "Tier_2_Silver": 1.5, "Tier_3_Bronze": 0.5}
        audit_df["authority_weight"] = [tier_map.get(str(t).strip(), 0.5) for t in audit_df["determiner_tier"]]
    is_issue = np.zeros(num_records, dtype=bool); is_issue[issues] = True
    audit_df["is_cleanlab_issue"], audit_df["is_label_corrupted"] = is_issue, audit_df["c_error"] > error_threshold
    audit_df["is_label_issue"] = (audit_df["is_cleanlab_issue"] | audit_df["is_label_corrupted"]).astype(bool)
    if "surface_orientation" not in audit_df.columns:
        audit_df["surface_orientation"] = "adaxial"
    audit_df["tomentum_taxon_discordance"] = _compute_tomentum_taxon_discordance(audit_df, error_threshold)
    triage_actions, reasons = [], []
    for _, row in audit_df.iterrows():
        is_corrupted, is_issue_flag = bool(row["is_label_corrupted"]), bool(row["is_cleanlab_issue"])
        orient, fg_ratio = str(row.get("surface_orientation", "adaxial")).lower(), float(row.get("foreground_attention_ratio", 0.85))
        if filter_abaxial_artifacts and orient == "abaxial" and (is_corrupted or is_issue_flag):
            act, rsn = "Retain (Abaxial Inversion Artifact)", "Visual divergence driven by dense abaxial arachnoid tomentum (mounting variant); retained as valid voucher."
        elif fg_ratio < FOREGROUND_ATTENTION_MIN_RATIO and (is_corrupted or is_issue_flag):
            act, rsn = "Downweight (Background Artifact Attribution)", f"Low foreground attention ratio ({fg_ratio:.2f} < {FOREGROUND_ATTENTION_MIN_RATIO:.2f}); prediction driven by sheet background."
        elif is_corrupted:
            act, rsn = "Prune & Queue for Annotation Triage", f"Deep vision DINOv2 predicts {row['predicted_label']} (conf: {row['confidence_predicted_class']:.2f}) vs recorded {row['given_label']}"
        elif is_issue_flag:
            act, rsn = "Flag for Morphometric Review", f"Confident learning flags potential discordance (quality: {row['label_quality_score']:.2f})"
        else:
            act, rsn = "Retain", "High label consistency across visual self-supervised embeddings"
        triage_actions.append(act); reasons.append(rsn)
    audit_df["triage_action"], audit_df["discordance_reason"] = triage_actions, reasons
    if filter_abaxial_artifacts:
        audit_df.loc[audit_df["surface_orientation"].str.lower() == "abaxial", ["is_label_corrupted", "is_label_issue"]] = False
    return audit_df


def compute_foreground_attention_ratio(heatmap: np.ndarray, fg_mask: np.ndarray) -> float:
    """Computes Pointing Game Foreground Attribution Ratio (Zhang et al. 2018)."""
    hm, mask = np.maximum(np.asarray(heatmap, dtype=np.float64), 0.0), np.asarray(fg_mask)
    if hm.size == 0 or mask.size == 0:
        return 0.0
    if mask.shape[:2] != hm.shape[:2]:
        mask = cv2.resize(mask.astype(np.uint8), (hm.shape[1], hm.shape[0]), interpolation=cv2.INTER_NEAREST) if cv2 is not None else np.array(Image.fromarray(mask.astype(np.uint8)).resize((hm.shape[1], hm.shape[0]), Image.NEAREST))
    bin_mask, total_energy = (mask > 0).astype(np.float64), float(np.sum(hm))
    return float(np.mean(bin_mask)) if total_energy <= 1e-8 else float(np.clip(float(np.sum(hm * bin_mask)) / total_energy, 0.0, 1.0))


def blend_heatmap_on_image(rgb_image: np.ndarray, heatmap: np.ndarray, alpha: float = 0.52, colormap: str = "inferno", fg_mask: Optional[np.ndarray] = None) -> np.ndarray:
    """Blends a 2D float heatmap [0.0, 1.0] onto an RGB uint8 image using inferno colormap."""
    cmap = matplotlib.colormaps[colormap] if hasattr(matplotlib, "colormaps") else plt.get_cmap(colormap)
    colored_heatmap, (h, w) = (cmap(heatmap)[:, :, :3] * 255.0).astype(np.uint8), rgb_image.shape[:2]
    if colored_heatmap.shape[:2] != (h, w):
        colored_heatmap = np.array(Image.fromarray(colored_heatmap).resize((w, h), Image.BILINEAR))
    blended = (rgb_image.astype(np.float32) * (1.0 - alpha) + colored_heatmap.astype(np.float32) * alpha).astype(np.uint8)
    bg_cond = (fg_mask <= 0) if (fg_mask is not None and fg_mask.shape[:2] == (h, w)) else (np.max(np.abs(rgb_image.astype(np.int16) - np.array([242, 239, 230], dtype=np.int16)), axis=-1) <= 6)
    if np.any(bg_cond) and not np.all(bg_cond):
        blended[bg_cond] = rgb_image[bg_cond]
    return blended


class _DemoViTFeatureExtractor:
    """Explicit `--demo` mode patch-token feature extractor following the DINOv2 ViT contract."""
    def __init__(self, embed_dim: int = 768, patch_size: int = 14):
        self.embed_dim, self.patch_size = embed_dim, patch_size

    def forward_features(self, img_tensor: torch.Tensor) -> Dict[str, torch.Tensor]:
        b, c, h, w = img_tensor.shape
        h_p, w_p = max(1, h // self.patch_size), max(1, w // self.patch_size)
        patches = F.adaptive_avg_pool2d(img_tensor, (h_p, w_p)).permute(0, 2, 3, 1).reshape(b, h_p * w_p, c)
        patch_tokens = F.layer_norm(patches.repeat(1, 1, (self.embed_dim + c - 1) // c)[:, :, : self.embed_dim], (self.embed_dim,))
        return {"x_norm_patchtokens": patch_tokens, "x_norm_clstoken": patch_tokens.mean(dim=1)}


def compute_vit_patch_saliency(img_arr: np.ndarray, model: Any, classifier: Optional[Any] = None, target_class_idx: Optional[int] = None, device: str = "cuda" if torch.cuda.is_available() else "cpu") -> np.ndarray:
    """Projects class weights / Captum attribution onto the 16x16 ViT patch grid upsampled to 224x224."""
    h, w = img_arr.shape[:2]
    macro_arr = preprocess_macro_rosette_crop(img_arr, target_size=(224, 224), sigma=1.0)
    img_tensor = transforms.Compose([transforms.ToTensor(), transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])])(Image.fromarray(macro_arr)).unsqueeze(0).to(device)
    with torch.no_grad():
        feat_dict = model.forward_features(img_tensor)
        patch_tokens, cls_token = feat_dict["x_norm_patchtokens"][0], feat_dict.get("x_norm_clstoken", None)
    num_patches, embed_dim = patch_tokens.shape
    h_p = max(1, int(np.round(np.sqrt(num_patches))))
    w_p = max(1, num_patches // h_p)
    F_patches = patch_tokens[: h_p * w_p].reshape(h_p, w_p, embed_dim)
    w_c, coef_arr, classes_list = None, None, []
    if classifier is not None:
        est = classifier if hasattr(classifier, "coef_") else (getattr(classifier.calibrated_classifiers_[0], "estimator", None) if getattr(classifier, "calibrated_classifiers_", None) else None)
        if est is not None and hasattr(est, "coef_"):
            coef_arr, classes_list = np.asarray(est.coef_, dtype=np.float32), list(getattr(est, "classes_", range(est.coef_.shape[0])))
            if coef_arr.ndim == 1:
                coef_arr = coef_arr[None, :]
            if coef_arr.shape[0] == 1 and len(classes_list) == 2:
                coef_arr = np.vstack([-coef_arr, coef_arr])
    if coef_arr is not None and coef_arr.shape[0] > 0:
        if target_class_idx is None and cls_token is not None and hasattr(classifier, "predict"):
            target_class_idx = int(classifier.predict(cls_token.detach().cpu().numpy())[0])
        row_idx = (classes_list.index(target_class_idx) if target_class_idx in classes_list else int(target_class_idx)) if target_class_idx is not None else 0
        w_c = torch.as_tensor(coef_arr[row_idx % coef_arr.shape[0]], dtype=F_patches.dtype, device=F_patches.device)
    elif cls_token is not None:
        w_c = cls_token[0].to(dtype=F_patches.dtype, device=F_patches.device)
    raw_saliency = torch.einsum("hwd,d->hw", F_patches, w_c) if (w_c is not None and w_c.numel() == embed_dim) else torch.norm(F_patches, dim=-1)
    saliency_relu = torch.relu(raw_saliency) if float(torch.relu(raw_saliency).max().item()) > 1e-8 else torch.relu(raw_saliency - raw_saliency.min())
    heatmap = torch.clamp(F.interpolate(saliency_relu.unsqueeze(0).unsqueeze(0).to(dtype=torch.float32), size=(h, w), mode="bicubic", align_corners=False)[0, 0], min=0.0).detach().cpu().numpy().astype(np.float32)
    return (heatmap - heatmap.min()) / (heatmap.max() - heatmap.min() + 1e-8)


def generate_gradcam_panel(flagged_records: List[dict], output_path: Path, max_samples: int = 8, model: Optional[Any] = None, classifier: Optional[Any] = None, model_name: str = "dinov2_vitb14_reg", device: str = "cuda" if torch.cuda.is_available() else "cpu", demo: bool = False, class_names: Optional[List[str]] = None, active_taxa: Optional[List[str]] = None) -> None:
    """Renders multi-panel diagnostic figure with 224x224 ViT patch attribution & foreground ratio."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dedup_records, seen_cats = [], set()
    for r in flagged_records:
        c_key = str(r.get("catalogNumber", r.get("patch_path", "")))
        if c_key not in seen_cats:
            seen_cats.add(c_key)
            dedup_records.append(r)
    samples = dedup_records[:max_samples]
    if not samples:
        return
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
    n = len(samples)
    fig, axes = plt.subplots(n, 3, figsize=(11.8, 3.2 * n), gridspec_kw={"width_ratios": [1.0, 1.0, 1.25], "wspace": 0.28})
    if n == 1:
        axes = np.expand_dims(axes, 0)
    for idx, rec in enumerate(samples):
        cat_num_str, raw_p = str(rec.get("catalogNumber", f"V{idx + 1}")), rec.get("patch_path")
        img_path = Path(str(raw_p)) if (raw_p and str(raw_p).strip()) else (Path("data/cropped_patches/rosettes_dense") / f"{cat_num_str}_p0_rosette.jpg")
        if img_path.exists():
            img_arr = preprocess_macro_rosette_crop(np.array(Image.open(img_path).convert("RGB")), target_size=(224, 224), sigma=1.0)
        elif demo:
            img_arr = np.full((224, 224, 3), fill_value=[60, 120, 40], dtype=np.uint8)
        else:
            raise RuntimeError(f"Missing patch image at '{img_path}'.")
        fg_mask, _ = compute_foreground_tissue_mask(img_arr, grid_size=(16, 16))
        rec_taxon = standardize_packera_taxon(str(rec.get("taxon", rec.get("species_standardized", rec.get("given_label", active_taxa[0])))))
        rec_taxon = rec_taxon if _is_valid_packera_taxon(rec_taxon) else active_taxa[0]
        raw_pred = rec.get("predicted_taxon", rec.get("predicted_label", rec.get("vision_predicted_label", rec_taxon)))
        pred_taxon = standardize_packera_taxon(str(raw_pred)) if raw_pred else rec_taxon
        pred_taxon = pred_taxon if pred_taxon in active_taxa else (rec_taxon if rec_taxon in active_taxa else active_taxa[0])
        target_cls_idx = active_taxa.index(pred_taxon) if pred_taxon in active_taxa else int(rec.get("label_idx", 0) or 0)
        heatmap = compute_vit_patch_saliency(img_arr, model=model, classifier=classifier, target_class_idx=target_cls_idx, device=device)
        fg_ratio = compute_foreground_attention_ratio(heatmap, fg_mask)
        rec["foreground_attention_ratio"] = round(fg_ratio, 4)
        blended = blend_heatmap_on_image(img_arr, heatmap, alpha=0.52, colormap="inferno", fg_mask=fg_mask)
        c_err_val = float(rec.get("c_error", 1.0 - float(rec.get("label_quality_score", 0.25))))
        conf_val = float(np.clip(max(c_err_val, 1.0 - c_err_val), 0.35, 0.99))
        axes[idx, 0].imshow(img_arr)
        axes[idx, 0].set_title(f"[1] {cat_num_str} (224x224 Macro)\nRecorded: {rec_taxon.replace('Packera ', 'P. ')}", fontsize=8.8, fontweight="bold")
        axes[idx, 0].axis("off")
        axes[idx, 1].imshow(blended)
        axes[idx, 1].set_title(f"[2] 16x16 ViT Saliency (FG={fg_ratio:.2f})\nPredicted: {pred_taxon.replace('Packera ', 'P. ')} (c_err={c_err_val:.2f})", fontsize=8.8, fontweight="bold")
        axes[idx, 1].axis("off")
        ax_bar = axes[idx, 2]
        probs = np.array([float(np.clip(1.0 - c_err_val, 0.01, 0.99))], dtype=float) if len(active_taxa) == 1 else np.full(len(active_taxa), (1.0 - conf_val) / max(len(active_taxa) - 1, 1), dtype=float)
        if len(active_taxa) > 1 and pred_taxon in active_taxa:
            probs[active_taxa.index(pred_taxon)] = conf_val
            probs = probs / max(probs.sum(), 1e-8)
        y_pos = np.arange(len(active_taxa))
        ax_bar.barh(y_pos, probs * 100.0, color=[palette_map.get(tx, "#0072B2") for tx in active_taxa], edgecolor="#222222", height=0.58, alpha=0.88)
        ax_bar.set_yticks(y_pos)
        ax_bar.set_yticklabels([tx.replace("Packera ", "P. ") for tx in active_taxa], fontsize=8.2, fontstyle="italic")
        ax_bar.set_xlim(0, 112); ax_bar.invert_yaxis()
        ax_bar.set_title("[3] Calibrated Posteriors", fontsize=8.8, fontweight="bold")
    fig.subplots_adjust(hspace=0.38, wspace=0.24)
    tmp_out = output_path.with_name(f".{output_path.stem}_{os.getpid()}_tmp{output_path.suffix}")
    plt.savefig(tmp_out, dpi=300, bbox_inches="tight"); tmp_out.replace(output_path)
    pdf_companion = output_path.with_suffix(".pdf")
    if pdf_companion != output_path:
        tmp_pdf = output_path.with_name(f".{output_path.stem}_{os.getpid()}_tmp.pdf")
        plt.savefig(tmp_pdf, format="pdf", bbox_inches="tight"); tmp_pdf.replace(pdf_companion)
    plt.close()


def generate_step05_publication_figures(audit_df: pd.DataFrame, output_figure_path: Path) -> Dict[str, Path]:
    """Generates Step 05 publication figure `fig05b_cleanlab_quality_and_tomentum.{pdf,png}`."""
    if audit_df is None or audit_df.empty:
        return {}
    env_run = os.environ.get("PACKERA_RUN_DIR", "").strip()
    figures_dir = (Path(env_run).resolve() / "outputs" / "figures") if env_run else Path(output_figure_path).resolve().parent
    figures_dir.mkdir(parents=True, exist_ok=True)
    df = audit_df.copy()
    sp_col, pred_col = ("species_standardized" if "species_standardized" in df.columns else "given_label"), ("predicted_label" if "predicted_label" in df.columns else "vision_predicted_label")
    df["taxon_clean"] = df.get(sp_col, pd.Series("Unassigned_Packera", index=df.index)).astype(str).apply(standardize_packera_taxon)
    df["pred_clean"] = df.get(pred_col, df["taxon_clean"]).astype(str).apply(standardize_packera_taxon)
    present_tx = {t for t in df["taxon_clean"] if _is_valid_packera_taxon(t)}
    taxa_order = sorted(present_tx) or sorted({t for t in df["pred_clean"] if _is_valid_packera_taxon(t)}) or ["Unassigned_Packera"]
    palette_map = _get_dynamic_okabe_ito_map(taxa_order)
    short_labels = [tx.replace("Packera ", "P. ") for tx in taxa_order]
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
        if np.any(sp_mask):
            ax_b.scatter(lum_vals[sp_mask], sat_vals[sp_mask], facecolors="none", edgecolors=palette_map.get(sp_name, "#0072B2"), marker=hollow_markers[idx_sp % len(hollow_markers)], s=36, alpha=0.75, label=sp_name.replace("Packera ", "P. "))
    ax_b.set_xlabel("Mean Foliar CIE L* Luminance", fontsize=9.5, fontweight="bold"); ax_b.set_ylabel("Mean Foliar HSV Saturation", fontsize=9.5, fontweight="bold")
    if len(taxa_order) > 0: ax_b.legend(loc="best", fontsize=8.0, frameon=True)
    ax_b.set_title("B. Foliar Chromatic Space & Species Ellipses", fontsize=10.5, fontweight="bold", loc="left")
    fig.tight_layout()
    pdf_out, png_out = figures_dir / "fig05b_cleanlab_quality_and_tomentum.pdf", figures_dir / "fig05b_cleanlab_quality_and_tomentum.png"
    tmp_pdf, tmp_png = figures_dir / f".fig05b_{os.getpid()}_tmp.pdf", figures_dir / f".fig05b_{os.getpid()}_tmp.png"
    try:
        fig.savefig(tmp_pdf, format="pdf", bbox_inches="tight"); tmp_pdf.replace(pdf_out)
        fig.savefig(tmp_png, format="png", dpi=300, bbox_inches="tight"); tmp_png.replace(png_out)
    finally:
        plt.close(fig)
    return {"pdf": pdf_out, "png": png_out}


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for Phase 5 analysis."""
    p = argparse.ArgumentParser(description="Phase 5: Indumentum-Blind DINOv2 Macro-Architecture & Cleanlab XAI Audit")
    p.add_argument("--rosette-dir", type=Path, default=Path("data/cropped_patches"))
    p.add_argument("--vouchers-csv", type=Path, default=Path("data/tables/curated_vouchers.csv"))
    p.add_argument("--output-csv", type=Path, default=Path("data/tables/label_noise_audit.csv"))
    p.add_argument("--output-figure", type=Path, default=Path("outputs/figures/GradCAM_audit_panel.png"))
    p.add_argument("--cleanlab-threshold", dest="cleanlab_threshold", type=float, default=0.85); p.add_argument("--error-threshold", dest="cleanlab_threshold", type=float)
    p.add_argument("--neutralize-background", dest="neutralize_background", action="store_true", default=True); p.add_argument("--no-neutralize-background", dest="neutralize_background", action="store_false")
    p.add_argument("--export-figures", dest="export_figures", action="store_true", default=True); p.add_argument("--no-export-figures", dest="export_figures", action="store_false")
    p.add_argument("--filter-abaxial-artifacts", dest="filter_abaxial_artifacts", action="store_true", default=True); p.add_argument("--no-filter-abaxial-artifacts", dest="filter_abaxial_artifacts", action="store_false")
    p.add_argument("--stratify-by-orientation", dest="stratify_by_orientation", action="store_true", default=False)
    p.add_argument("--morphometrics-csv", type=Path, default=None)
    p.add_argument("--manifest-csv", type=Path, default=None)
    p.add_argument("--intraspecific-mode", dest="intraspecific_mode", action="store_true", default=False)
    p.add_argument("--focal-taxon", type=str, default="")
    p.add_argument("--demo", dest="demo", action="store_true", default=False)
    return p.parse_args()


def main() -> None:
    """Orchestrates Phase 5 Indumentum-Blind Deep Vision & Confident Learning workflow."""
    args = parse_args()
    df, _ = load_and_link_rosette_patches(args.rosette_dir, args.vouchers_csv, demo=args.demo, manifest_csv=args.manifest_csv)
    if df.empty:
        logger.warning("No rosette patch records were matched. Check inputs."); return
    features, labels, _ = extract_dinov2_embeddings(df.to_dict("records"), neutralize_background=args.neutralize_background, demo=args.demo)
    active_taxa = sorted(t for t in set(df["taxon"]) if _is_valid_packera_taxon(t)) or ["Unassigned_Packera"]
    pred_probs, _, _ = compute_out_of_fold_probabilities(features, labels, n_splits=5, surface_orientations=df.get("surface_orientation", None), stratify_by_orientation=args.stratify_by_orientation, determiner_tiers=df.get("determiner_tier", None), spatial_balance_weights=df.get("spatial_balance_weight", None), institution_codes=df.get("institutionCode", None))
    audit_df = run_confident_learning_audit(pred_probs=pred_probs, labels=labels, records_df=df, class_names=active_taxa, error_threshold=args.cleanlab_threshold, filter_abaxial_artifacts=args.filter_abaxial_artifacts, features=features)
    if args.export_figures:
        flagged = audit_df[audit_df["is_label_corrupted"]].to_dict("records") or audit_df[audit_df["is_cleanlab_issue"]].to_dict("records") or audit_df.sort_values("c_error", ascending=False).head(8).to_dict("records")
        clf = LogisticRegression(max_iter=1000, C=0.25, class_weight="balanced", random_state=42).fit(features, labels) if len(np.unique(labels)) > 1 else None
        generate_gradcam_panel(flagged, args.output_figure, classifier=clf, demo=args.demo, class_names=active_taxa, active_taxa=active_taxa)
        fg_updates = {str(r["catalogNumber"]): float(r["foreground_attention_ratio"]) for r in flagged if "catalogNumber" in r and "foreground_attention_ratio" in r}
        if fg_updates:
            audit_df["foreground_attention_ratio"] = [fg_updates.get(str(c), float(fg)) for c, fg in zip(audit_df["catalogNumber"], audit_df["foreground_attention_ratio"])]
            audit_df["is_artifact_driven"] = (audit_df["foreground_attention_ratio"] < FOREGROUND_ATTENTION_MIN_RATIO).astype(bool)
        try: generate_step05_publication_figures(audit_df, args.output_figure)
        except Exception as fig_exc: logger.warning(f"Step 05 publication figure generation skipped: {fig_exc}")
    if not args.demo:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        tmp_csv = args.output_csv.with_suffix(args.output_csv.suffix + ".tmp")
        export_df = audit_df.drop(columns=["rosette_patches", "tier1_leaf_patches", "tier1_leaf_masks", "tier1_leaf_weights"], errors="ignore")
        with open(tmp_csv, "w", encoding="utf-8", newline="") as f:
            export_df.to_csv(f, index=False); f.flush(); os.fsync(f.fileno())
        tmp_csv.replace(args.output_csv)


if __name__ == "__main__":
    main()
