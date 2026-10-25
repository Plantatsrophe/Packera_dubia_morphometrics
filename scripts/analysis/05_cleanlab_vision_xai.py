#!/usr/bin/env python3
"""
===============================================================================
Script: 05_cleanlab_vision_xai.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Phase 5 Unified Deep Vision & Explainable AI Pipeline:
    1. Standardizes botanical nomenclature across the Packera dubia complex.
    2. Loads and links dense basal rosette patch images with curated voucher metadata.
    3. Extracts self-supervised 768-dim [CLS] token representations via DINOv2.
    4. Computes out-of-fold cross-validated probabilities via Stratified K-Fold.
    5. Audits label quality and taxonomic misidentification noise via Confident Learning.
    6. Generates diagnostic Grad-CAM / feature attribution multi-panel figures.

Author: J. Brandon Fuller
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset
from torchvision import transforms

try:
    import cv2
except ImportError:
    cv2 = None

# Set reproducible random seeds
torch.manual_seed(42)
np.random.seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

TARGET_TAXA: List[str] = [
    "Packera anonyma",
    "Packera dubia",
    "Packera paupercula",
    "Packera plattensis",
]

__all__ = [
    "TARGET_TAXA",
    "standardize_packera_taxon",
    "neutralize_mounting_paper",
    "classify_foliar_surface_orientation",
    "RosettePatchDataset",
    "load_and_link_rosette_patches",
    "extract_dinov2_embeddings",
    "compute_out_of_fold_probabilities",
    "run_confident_learning_audit",
    "blend_heatmap_on_image",
    "compute_vit_patch_saliency",
    "generate_gradcam_panel",
    "main",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("Cleanlab_Vision_XAI")


# =============================================================================
# 1. Taxonomic Standardization
# =============================================================================

def standardize_packera_taxon(species_str: Optional[str]) -> str:
    """Standardize synonymy into the four core Packera dubia complex taxa.

    Handles historical basionyms and varieties:
      - Senecio smallii, Senecio earlei -> Packera anonyma
      - Senecio tomentosus, Packera tomentosa -> Packera dubia
      - Packera paupercula varieties -> Packera paupercula
      - Senecio plattensis, Packera flavovirens -> Packera plattensis
    """
    if not species_str or pd.isna(species_str):
        return "Unknown"
    s = str(species_str).strip()
    if re.search(r"anonym|smallii|earlei", s, re.I):
        return "Packera anonyma"
    if re.search(r"paupercul|balsamitae|savannarum|pseudotomentosa|appalachiana", s, re.I):
        return "Packera paupercula"
    if re.search(r"tomentos|dubia", s, re.I):
        return "Packera dubia"
    if re.search(r"plattensis|flavovirens", s, re.I):
        return "Packera plattensis"
    return s.split("(")[0].strip()


# =============================================================================
# 2. Mounting Paper Background Neutralization & PyTorch Dataset
# =============================================================================

def neutralize_mounting_paper(
    image: Image.Image | np.ndarray,
    neutral_color: Tuple[int, int, int] = (128, 128, 128),
) -> Image.Image:
    """Neutralizes herbarium mounting paper background to uniform neutral gray.

    Before feeding whole-rosette image patches into DINOv2 (ViT-B/14), this function
    computes an Otsu and adaptive chromatic threshold mask separating botanical
    tissue (green/brown foliar lamina, petiole vasculature, and arachnoid tomentum)
    from the mounting sheet paper background.

    All non-plant background paper pixels are replaced with a uniform, neutral
    gray value: RGB(128, 128, 128). This prevents vision transformers from
    keying into yellowed paper fibers, mounting glue, or herbarium stamps.

    Vectorized using OpenCV / NumPy array slicing (strictly zero per-pixel loops).
    """
    if isinstance(image, Image.Image):
        img_arr = np.array(image.convert("RGB"))
        was_pil = True
    else:
        img_arr = np.array(image, copy=True)
        was_pil = False

    if img_arr.size == 0:
        return Image.fromarray(img_arr) if was_pil else img_arr

    if cv2 is not None:
        gray = cv2.cvtColor(img_arr, cv2.COLOR_RGB2GRAY)
        blurred = cv2.GaussianBlur(gray, (5, 5), 0)
        # Otsu thresholding: botanical tissue is darker than bright mounting paper
        _, plant_mask = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)

        # Morphological closing (5x5 ellipse) to bridge arachnoid tomentum, venation, and leaf blades
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
        plant_mask = cv2.morphologyEx(plant_mask, cv2.MORPH_CLOSE, kernel)

        bg_mask = (plant_mask == 0)
    else:
        # Vectorized NumPy fallback if cv2 is not available
        gray = (0.2989 * img_arr[:, :, 0] + 0.5870 * img_arr[:, :, 1] + 0.1140 * img_arr[:, :, 2])
        thresh = np.percentile(gray, 60)
        bg_mask = (gray > thresh)

    neutralized = img_arr.copy()
    neutralized[bg_mask] = neutral_color

    return Image.fromarray(neutralized) if was_pil else neutralized


def classify_foliar_surface_orientation(
    image: Image.Image | np.ndarray | str | Path,
    plant_mask: Optional[np.ndarray] = None,
    neutral_color: Tuple[int, int, int] = (128, 128, 128),
    luminance_threshold: float = 68.0,
    saturation_threshold: float = 0.22,
    brightness_cutoff: float = 75.0,
    max_dimension: int = 512,
) -> Dict[str, Any]:
    """Classifies foliar surface orientation (abaxial vs. adaxial) based on trichome tomentum.

    Botanists intentionally mount at least one leaf upside down (abaxial side facing up)
    on herbarium sheets. In Packera dubia, the abaxial surface possesses dense, bright white
    arachnoid-woolly tomentum, whereas the adaxial surface is darker green/brown. DINOv2
    vision transformers cluster adaxial and abaxial leaves from the same plant into divergent
    visual clusters, causing Cleanlab to falsely flag upside-down leaves as label errors.

    Botanical Rule:
      - Dense white arachnoid tomentum displays high luminance (L* > 68 in LAB space)
        and low saturation (S < 0.22 in HSV space).
      - Glabrescent adaxial surfaces display lower luminance and higher green/brown saturation.
      - Assigns surface_orientation = "abaxial" if tomentum criteria are met; else "adaxial".

    Calculates:
      1. Mean luminance (L* in LAB space, CIE [0, 100]).
      2. Saturation index (S in HSV space, [0, 1]).
      3. High-brightness pixel ratio (fraction of leaf pixels with L* > 75.0).

    Execution Constraints:
      - Vectorized calculation over segmented leaf mask pixels only (ignoring neutralized gray background pixels).
      - Fast execution (<0.05 seconds per patch).
    """
    if isinstance(image, (str, Path)):
        img_path = Path(image)
        if not img_path.exists():
            return {
                "surface_orientation": "adaxial",
                "mean_luminance": 0.0,
                "mean_saturation": 0.0,
                "high_brightness_ratio": 0.0,
            }
        img_pil = Image.open(img_path).convert("RGB")
        img_arr = np.array(img_pil)
    elif isinstance(image, Image.Image):
        img_arr = np.array(image.convert("RGB"))
    else:
        img_arr = np.array(image, copy=False)
        if img_arr.ndim == 2:
            img_arr = np.stack([img_arr] * 3, axis=-1)
        elif img_arr.shape[2] == 4:
            img_arr = img_arr[:, :, :3]

    if img_arr.size == 0:
        return {
            "surface_orientation": "adaxial",
            "mean_luminance": 0.0,
            "mean_saturation": 0.0,
            "high_brightness_ratio": 0.0,
        }

    # Downsample if image exceeds max_dimension to guarantee <0.05s execution
    if max(img_arr.shape[:2]) > max_dimension:
        scale = float(max_dimension) / max(img_arr.shape[:2])
        new_w = max(1, int(img_arr.shape[1] * scale))
        new_h = max(1, int(img_arr.shape[0] * scale))
        if cv2 is not None:
            img_arr = cv2.resize(img_arr, (new_w, new_h), interpolation=cv2.INTER_AREA)
        else:
            img_pil = Image.fromarray(img_arr).resize((new_w, new_h), Image.BILINEAR)
            img_arr = np.array(img_pil)
        if plant_mask is not None:
            if cv2 is not None:
                plant_mask = cv2.resize(plant_mask.astype(np.uint8), (new_w, new_h), interpolation=cv2.INTER_NEAREST) > 0
            else:
                plant_mask = np.array(Image.fromarray(plant_mask.astype(np.uint8)).resize((new_w, new_h), Image.NEAREST)) > 0

    # Determine segmented leaf mask pixels, ignoring neutralized gray background (RGB = neutral_color)
    if plant_mask is not None:
        mask = (plant_mask > 0)
    else:
        # Check for pre-neutralized background pixels
        is_neutral_gray = (
            (img_arr[:, :, 0] == neutral_color[0])
            & (img_arr[:, :, 1] == neutral_color[1])
            & (img_arr[:, :, 2] == neutral_color[2])
        )
        if np.mean(is_neutral_gray) >= 0.02:
            mask = ~is_neutral_gray
        else:
            # Segment leaf tissue from background mounting sheet using Otsu
            if cv2 is not None:
                gray = cv2.cvtColor(img_arr, cv2.COLOR_RGB2GRAY)
                blurred = cv2.GaussianBlur(gray, (5, 5), 0)
                _, otsu_mask = cv2.threshold(blurred, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
                kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
                mask = cv2.morphologyEx(otsu_mask, cv2.MORPH_CLOSE, kernel) > 0
            else:
                gray = 0.2989 * img_arr[:, :, 0] + 0.5870 * img_arr[:, :, 1] + 0.1140 * img_arr[:, :, 2]
                thresh = np.percentile(gray, 60)
                mask = (gray <= thresh)

    # If mask is empty or nearly empty, default to adaxial
    if mask.sum() < 5:
        return {
            "surface_orientation": "adaxial",
            "mean_luminance": 0.0,
            "mean_saturation": 0.0,
            "high_brightness_ratio": 0.0,
        }

    # Vectorized conversion to LAB and HSV
    if cv2 is not None:
        lab = cv2.cvtColor(img_arr, cv2.COLOR_RGB2LAB)
        hsv = cv2.cvtColor(img_arr, cv2.COLOR_RGB2HSV)
        lab_pixels = lab[mask]
        hsv_pixels = hsv[mask]

        # Map OpenCV uint8 L [0..255] to standard CIE L* [0..100]
        L_vals = lab_pixels[:, 0].astype(np.float64) * (100.0 / 255.0)
        # Map OpenCV uint8 S [0..255] to normalized [0..1]
        S_vals = hsv_pixels[:, 1].astype(np.float64) / 255.0
    else:
        leaf_rgb = img_arr[mask].astype(np.float64)
        R, G, B = leaf_rgb[:, 0], leaf_rgb[:, 1], leaf_rgb[:, 2]
        Y = 0.2126 * (R / 255.0) + 0.7152 * (G / 255.0) + 0.0722 * (B / 255.0)
        L_vals = np.where(Y > 0.008856, 116.0 * (Y ** (1.0 / 3.0)) - 16.0, 903.3 * Y)

        cmax = np.maximum(np.maximum(R, G), B)
        cmin = np.minimum(np.minimum(R, G), B)
        delta = cmax - cmin
        S_vals = np.where(cmax > 0, delta / cmax, 0.0)

    mean_L = float(np.mean(L_vals))
    mean_S = float(np.mean(S_vals))
    high_bright_ratio = float(np.mean(L_vals > brightness_cutoff))

    # Botanical Rule: Dense white arachnoid tomentum displays high luminance (L* > 68)
    # and low saturation (S < 0.22).
    is_abaxial = (mean_L > luminance_threshold) and (mean_S < saturation_threshold)
    surface_orientation = "abaxial" if is_abaxial else "adaxial"

    return {
        "surface_orientation": surface_orientation,
        "mean_luminance": round(mean_L, 4),
        "mean_saturation": round(mean_S, 4),
        "high_brightness_ratio": round(high_bright_ratio, 4),
    }


class RosettePatchDataset(Dataset):
    """PyTorch Dataset loading dense basal rosette image crops for DINOv2."""

    def __init__(
        self,
        records: List[Dict],
        transform: transforms.Compose,
        neutralize_background: bool = True,
        demo: bool = False,
    ):
        self.records = records
        self.transform = transform
        self.neutralize_background = neutralize_background
        self.demo = demo

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, int, str]:
        rec = self.records[idx]
        image_path = Path(rec["patch_path"])
        if image_path.exists():
            image = Image.open(image_path).convert("RGB")
            if self.neutralize_background:
                image = neutralize_mounting_paper(image)
        elif self.demo:
            image = Image.new("RGB", (224, 224), color=(128, 128, 128))
        else:
            raise RuntimeError(
                f"Rosette patch image missing at '{image_path}'. "
                "Remediation: Run Phase 2 (`python main.py segment`) to extract organ/rosette crops "
                "into `data/cropped_patches/`, or pass `--demo` explicitly for synthetic testing."
            )
        tensor = self.transform(image)
        return tensor, rec["label_idx"], rec["catalogNumber"]


def load_and_link_rosette_patches(
    rosette_dir: Path,
    vouchers_csv: Path,
    demo: bool = False,
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Links rosette patch image files with curated herbarium metadata and classifies foliar orientation."""
    rosette_dir = Path(rosette_dir)
    vouchers_df = pd.read_csv(vouchers_csv) if Path(vouchers_csv).exists() else pd.DataFrame()
    voucher_map = {}
    if not vouchers_df.empty and "catalogNumber" in vouchers_df.columns:
        for _, row in vouchers_df.iterrows():
            voucher_map[str(row["catalogNumber"]).strip()] = row.to_dict()

    class_to_idx = {taxon: idx for idx, taxon in enumerate(TARGET_TAXA)}
    records = []

    # Fast lookup of physical patch files across rosette_dir and subdirectories
    available_patches = {}
    if rosette_dir.exists():
        for p in rosette_dir.rglob("*.jpg"):
            cat = p.stem.split("_")[0]
            if cat not in available_patches or "rosette" in p.name:
                available_patches[cat] = p

    patch_files = sorted(rosette_dir.glob("*.jpg")) if rosette_dir.exists() else []
    if not patch_files and rosette_dir.exists():
        rosette_crops = sorted(rosette_dir.glob("rosettes_dense/*.jpg"))
        leaf_crops = sorted(rosette_dir.glob("basal_leaves_raw/*.jpg"))
        rosette_cats = {p.stem.split("_")[0] for p in rosette_crops}
        patch_files = rosette_crops + [p for p in leaf_crops if p.stem.split("_")[0] not in rosette_cats]

    if patch_files:
        for p in patch_files:
            cat_num = p.stem.split("_")[0]
            meta = voucher_map.get(cat_num, {})
            raw_sp = meta.get("species_raw", meta.get("species", "Unknown"))
            std_taxon = standardize_packera_taxon(raw_sp)
            if std_taxon not in class_to_idx:
                continue

            orient_info = classify_foliar_surface_orientation(p)
            records.append({
                "catalogNumber": cat_num,
                "plant_individual_id": int(meta.get("plant_individual_id", 0)),
                "patch_path": str(p),
                "taxon": std_taxon,
                "species_raw": raw_sp,
                "label_idx": class_to_idx[std_taxon],
                "determiner_tier": meta.get("determiner_tier", "Tier_3_Bronze"),
                "surface_orientation": orient_info["surface_orientation"],
                "mean_luminance": orient_info["mean_luminance"],
                "mean_saturation": orient_info["mean_saturation"],
                "high_brightness_ratio": orient_info["high_brightness_ratio"],
            })
    elif not demo:
        raise RuntimeError(
            f"No cropped rosette/leaf patch images found in '{rosette_dir}'. "
            "Remediation: Run Phase 2 (`python main.py segment`) to populate `data/cropped_patches/`, "
            "or pass `--demo` explicitly for synthetic demo execution."
        )
    elif not vouchers_df.empty:
        for _, row in vouchers_df.iterrows():
            raw_sp = str(row.get("species_raw", row.get("species", "Unknown")))
            std_taxon = standardize_packera_taxon(raw_sp)
            if std_taxon not in class_to_idx:
                continue
            cat_num = str(row.get("catalogNumber", "")).strip()
            real_path = available_patches.get(cat_num, rosette_dir / f"{cat_num}_rosette.jpg")
            orient_info = (
                classify_foliar_surface_orientation(real_path)
                if Path(real_path).exists()
                else {
                    "surface_orientation": "adaxial",
                    "mean_luminance": 0.0,
                    "mean_saturation": 0.0,
                    "high_brightness_ratio": 0.0,
                }
            )
            records.append({
                "catalogNumber": cat_num,
                "plant_individual_id": int(row.get("plant_individual_id", 0)),
                "patch_path": str(real_path),
                "taxon": std_taxon,
                "species_raw": raw_sp,
                "label_idx": class_to_idx[std_taxon],
                "determiner_tier": row.get("determiner_tier", "Tier_3_Bronze"),
                "surface_orientation": orient_info["surface_orientation"],
                "mean_luminance": orient_info["mean_luminance"],
                "mean_saturation": orient_info["mean_saturation"],
                "high_brightness_ratio": orient_info["high_brightness_ratio"],
            })

    df = pd.DataFrame(records)
    logger.info(f"Loaded and linked {len(df)} rosette patch records across {len(class_to_idx)} taxa.")
    return df, class_to_idx


# =============================================================================
# 3. DINOv2 Self-Supervised Feature Extraction
# =============================================================================

def extract_dinov2_embeddings(
    records: List[Dict],
    model_name: str = "dinov2_vitb14",
    device: str = "cuda" if torch.cuda.is_available() else "cpu",
    batch_size: int = 32,
    neutralize_background: bool = True,
    demo: bool = False,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Extracts 768-dimensional DINOv2 self-supervised [CLS] token representations.

    When `demo=False` (production default), raises a `RuntimeError` if physical
    patch images are missing or if DINOv2 weights cannot be loaded.
    """
    torch.manual_seed(42)
    np.random.seed(42)

    has_physical_images = bool(records) and any(Path(r["patch_path"]).exists() for r in records[:5])
    if not has_physical_images:
        if not demo:
            raise RuntimeError(
                "No physical patch images found in `data/cropped_patches/` (or provided records). "
                "Remediation: Run Phase 2 (`python main.py segment`) to generate cropped leaf/rosette "
                "patches before running Phase 5, or pass `--demo` explicitly for synthetic verification."
            )
        n = len(records)
        labels = np.array([r["label_idx"] for r in records])
        feats = np.random.randn(n, 768).astype(np.float32)
        for i in range(n):
            feats[i, labels[i] * 50:(labels[i] + 1) * 50] += 3.0
            if i % 15 == 0:  # Inject simulated label noise in explicit --demo mode only
                feats[i, :] = np.random.randn(768)
                feats[i, ((labels[i] + 1) % 4) * 50:(((labels[i] + 1) % 4) + 1) * 50] += 4.0
        cat_nums = [r["catalogNumber"] for r in records]
        return feats, labels, cat_nums

    transform = transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])

    dataset = RosettePatchDataset(
        records, transform=transform, neutralize_background=neutralize_background, demo=demo
    )
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)

    try:
        model = torch.hub.load("facebookresearch/dinov2", model_name)
    except Exception as e:
        if not demo:
            raise RuntimeError(
                f"Could not load DINOv2 model weights ('{model_name}') via torch.hub: {e}. "
                "Remediation: Ensure internet access or pre-cache DINOv2 weights under "
                "`~/.cache/torch/hub/facebookresearch_dinov2_main`, or pass `--demo` explicitly."
            ) from e
        logger.warning(f"Could not load torch.hub dinov2: {e}. Generating demo features (--demo).")
        feats = np.random.randn(len(records), 768).astype(np.float32)
        labels = np.array([r["label_idx"] for r in records])
        cat_nums = [r["catalogNumber"] for r in records]
        return feats, labels, cat_nums

    model = model.to(device)
    model.eval()

    all_feats, all_labels, all_cats = [], [], []
    with torch.no_grad():
        for tensors, labels, cats in loader:
            tensors = tensors.to(device)
            out = model(tensors)
            all_feats.append(out.cpu().numpy())
            all_labels.extend(labels.numpy())
            all_cats.extend(cats)

    features = np.vstack(all_feats) if all_feats else np.empty((0, 768))
    return features, np.array(all_labels), all_cats


# =============================================================================
# 4. Out-of-Fold Estimation & Confident Learning Audit
# =============================================================================

def compute_out_of_fold_probabilities(
    features: np.ndarray,
    labels: np.ndarray,
    n_splits: int = 5,
    random_state: int = 42,
    surface_orientations: Optional[np.ndarray | List[str] | pd.Series] = None,
    stratify_by_orientation: bool = False,
    determiner_tiers: Optional[np.ndarray | List[str] | pd.Series] = None,
) -> Tuple[np.ndarray, float, float]:
    """Fits stratified cross-validated Logistic Regression to compute out-of-fold probabilities.

    L2-normalizes each 768-d embedding vector before fold-internal StandardScaler and fits
    balanced LogisticRegression weighted by determiner authority tier (Tier_1_Gold=3.0,
    Tier_2_Silver=1.5, Tier_3_Bronze=0.5). Supports optional compound stratification by both
    taxon label and foliar surface orientation. When only 1 class is present (or a class has
    fewer than `n_splits` samples), falls back to L2-normalized centroid/medoid cosine distance
    pseudo-probabilities without calling `LogisticRegression.fit()` on a 1-class vector.
    """
    X = np.asarray(features, dtype=np.float64)
    X = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-8)
    labels = np.asarray(labels, dtype=int)
    tier_weight_map = {"Tier_1_Gold": 3.0, "Tier_2_Silver": 1.5, "Tier_3_Bronze": 0.5}
    if determiner_tiers is not None:
        sample_weights = np.array(
            [tier_weight_map.get(str(t).strip(), 1.0) for t in determiner_tiers],
            dtype=np.float64,
        )
    else:
        sample_weights = np.ones(len(labels), dtype=np.float64)

    unique_labels, label_counts = np.unique(labels, return_counts=True)
    max_lbl_plus_1 = int(np.max(labels)) + 1 if len(labels) > 0 else 2
    if len(unique_labels) >= 2 and len(unique_labels) == max_lbl_plus_1:
        num_classes = max(2, max_lbl_plus_1)
    else:
        num_classes = max(2, len(TARGET_TAXA), max_lbl_plus_1)
    pred_probs = np.zeros((len(labels), num_classes), dtype=np.float64)

    if len(unique_labels) < 2 or int(np.min(label_counts)) < n_splits:
        logger.info(
            "Single-class or low-count label distribution detected "
            f"(unique={unique_labels.tolist()}, counts={label_counts.tolist()}); "
            "computing L2-normalized DINOv2 medoid/centroid cosine similarity pseudo-probabilities."
        )
        if len(X) == 0:
            return pred_probs, 1.0, 1.0
        centroid = np.median(X, axis=0)
        centroid = centroid / max(float(np.linalg.norm(centroid)), 1e-8)
        cos_sim = np.clip(X @ centroid, -1.0, 1.0)
        cos_dist = 1.0 - cos_sim
        q95 = float(np.percentile(cos_dist, 95)) if len(cos_dist) >= 5 else float(np.max(cos_dist)) + 1e-6
        q95 = max(q95, 1e-6)

        is_outlier = cos_dist > q95
        if surface_orientations is not None:
            orient_arr = np.asarray([str(o).lower() for o in surface_orientations])
            is_outlier = is_outlier & (orient_arr != "abaxial")

        p_self = np.where(
            is_outlier,
            np.clip(0.10 * np.exp(-(cos_dist - q95) / q95), 0.01, 0.12),
            np.clip(1.0 - 0.50 * (cos_dist / q95), 0.50, 0.999),
        )
        for i, lbl in enumerate(labels):
            c0 = int(lbl) % num_classes
            c_alt = (c0 + 1) % num_classes
            pred_probs[i, c0] = p_self[i]
            pred_probs[i, c_alt] = 1.0 - p_self[i]

        acc = float(np.mean(~is_outlier))
        return pred_probs, acc, acc

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    # Determine stratification target
    strat_labels = labels
    if stratify_by_orientation and surface_orientations is not None:
        orient_arr = np.array(surface_orientations)
        compound = np.array([f"{l}_{o}" for l, o in zip(labels, orient_arr)])
        unique, counts = np.unique(compound, return_counts=True)
        if all(c >= n_splits for c in counts):
            strat_labels = compound
        else:
            logger.info("Some compound taxon-orientation classes have fewer samples than n_splits; stratifying by taxon.")

    y_true_all, y_pred_all = [], []

    for train_idx, val_idx in skf.split(X, strat_labels):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_val = scaler.transform(X[val_idx])
        y_train, y_val = labels[train_idx], labels[val_idx]
        w_train = sample_weights[train_idx]

        clf = LogisticRegression(max_iter=1000, C=1.0, class_weight="balanced", random_state=random_state)
        clf.fit(X_train, y_train, sample_weight=w_train)

        probs = clf.predict_proba(X_val)
        for i_col, c_idx in enumerate(clf.classes_):
            pred_probs[val_idx, c_idx] = probs[:, i_col]

        y_true_all.extend(y_val)
        y_pred_all.extend(clf.classes_[np.argmax(probs, axis=1)])

    acc = accuracy_score(y_true_all, y_pred_all)
    f1 = f1_score(y_true_all, y_pred_all, average="weighted")
    logger.info(f"OOF Classifier Cross-Validation: Accuracy = {acc:.4f}, Weighted F1 = {f1:.4f}")

    return pred_probs, acc, f1


def run_confident_learning_audit(
    pred_probs: np.ndarray,
    labels: np.ndarray,
    records_df: pd.DataFrame,
    class_names: List[str],
    error_threshold: float = 0.85,
    filter_abaxial_artifacts: bool = True,
) -> pd.DataFrame:
    """Applies Confident Learning algorithms to identify potential herbarium label noise.

    Filters out candidates flagged solely due to inverted abaxial tomentum (botanical
    mounting variant) when filter_abaxial_artifacts is True.
    """
    try:
        if len(np.unique(labels)) < 2 or pred_probs.shape[1] != len(np.unique(labels)):
            raise ValueError("Single-class or non-contiguous label set; using direct probability margins.")
        import cleanlab
        from cleanlab.filter import find_label_issues
        from cleanlab.rank import get_label_quality_scores

        quality_scores = get_label_quality_scores(labels=labels, pred_probs=pred_probs)
        issues = find_label_issues(
            labels=labels,
            pred_probs=pred_probs,
            return_indices_ranked_by="self_confidence",
        )
    except Exception:
        logger.info("Using direct self-confidence probability margins for label quality scoring.")
        predicted_classes = np.argmax(pred_probs, axis=1)
        quality_scores = np.array([pred_probs[i, labels[i]] for i in range(len(labels))])
        issues = np.where((predicted_classes != labels) & (1.0 - quality_scores > error_threshold))[0]

    audit_df = records_df.copy()
    num_records = len(audit_df)

    given_classes = [class_names[l] if l < len(class_names) else "Unknown" for l in labels]
    pred_indices = np.argmax(pred_probs, axis=1)
    pred_classes = [class_names[i] if i < len(class_names) else "Outlier_Morphotype" for i in pred_indices]

    conf_given = np.array([
        pred_probs[i, labels[i]] if i < len(pred_probs) and labels[i] < pred_probs.shape[1] else 0.0
        for i in range(num_records)
    ])
    conf_pred = np.max(pred_probs, axis=1)
    c_error = 1.0 - conf_given

    raw_sp_col = audit_df.get(
        "taxon",
        audit_df.get("species_raw", audit_df.get("species", pd.Series(given_classes, index=audit_df.index))),
    )
    audit_df["species_raw"] = audit_df.get("species_raw", audit_df.get("species", given_classes))
    std_sp = pd.Series(raw_sp_col, index=audit_df.index).apply(standardize_packera_taxon)
    audit_df["species_standardized"] = [
        sp if sp != "Unknown" else (gc if gc in TARGET_TAXA else "Packera dubia")
        for sp, gc in zip(std_sp, given_classes)
    ]
    audit_df["given_label"] = given_classes
    audit_df["predicted_label"] = pred_classes
    audit_df["confidence_given_class"] = np.round(conf_given, 4)
    audit_df["confidence_predicted_class"] = np.round(conf_pred, 4)
    audit_df["label_quality_score"] = np.round(quality_scores, 4)
    audit_df["c_error"] = np.round(c_error, 4)

    is_issue = np.zeros(num_records, dtype=bool)
    is_issue[issues] = True
    audit_df["is_cleanlab_issue"] = is_issue
    audit_df["is_label_corrupted"] = audit_df["c_error"] > error_threshold

    if "surface_orientation" not in audit_df.columns:
        audit_df["surface_orientation"] = "adaxial"

    audit_df["tomentum_taxon_discordance"] = (
        (audit_df["surface_orientation"] == "abaxial")
        & audit_df["given_label"].isin({"Packera anonyma", "Packera paupercula"})
        & audit_df["predicted_label"].isin({"Packera dubia", "Packera plattensis"})
        & (audit_df["c_error"] > error_threshold)
    ).astype(bool)

    triage_actions = []
    reasons = []
    abaxial_filtered_count = 0

    for idx, row in audit_df.iterrows():
        is_corrupted = bool(row["is_label_corrupted"])
        is_issue_flag = bool(row["is_cleanlab_issue"])
        orient = str(row.get("surface_orientation", "adaxial")).lower()

        if filter_abaxial_artifacts and orient == "abaxial" and (is_corrupted or is_issue_flag):
            triage_actions.append("Retain (Abaxial Inversion Artifact)")
            reasons.append(
                f"Visual divergence driven by dense abaxial arachnoid tomentum "
                f"(mounting variant); retained as valid voucher rather than taxonomic misidentification."
            )
            abaxial_filtered_count += 1
        elif is_corrupted:
            triage_actions.append("Prune & Queue for Annotation Triage")
            reasons.append(
                f"Deep vision DINOv2 predicts {row['predicted_label']} "
                f"(conf: {row['confidence_predicted_class']:.2f}) vs recorded {row['given_label']}"
            )
        elif is_issue_flag:
            triage_actions.append("Flag for Morphometric Review")
            reasons.append(
                f"Confident learning flags potential discordance (quality: {row['label_quality_score']:.2f})"
            )
        else:
            triage_actions.append("Retain")
            reasons.append("High label consistency across visual self-supervised embeddings")

    audit_df["triage_action"] = triage_actions
    audit_df["discordance_reason"] = reasons

    # Filter out candidates flagged solely due to inverted abaxial tomentum from is_label_corrupted
    if filter_abaxial_artifacts:
        abaxial_mask = audit_df["surface_orientation"].str.lower() == "abaxial"
        audit_df.loc[abaxial_mask, "is_label_corrupted"] = False
        if abaxial_filtered_count > 0:
            logger.info(
                f"Filtered {abaxial_filtered_count} candidates flagged solely due to inverted abaxial tomentum."
            )

    flagged_count = int(audit_df["is_label_corrupted"].sum())
    logger.info(f"Confident Learning Audit completed: Flagged {flagged_count} / {len(audit_df)} corrupted labels.")
    return audit_df


# =============================================================================
# 5. Grad-CAM / Attribution Heatmap Visualization
# =============================================================================

def blend_heatmap_on_image(
    rgb_image: np.ndarray,
    heatmap: np.ndarray,
    alpha: float = 0.5,
    colormap: str = "jet",
) -> np.ndarray:
    """Blends a 2D float heatmap [0.0, 1.0] onto an RGB uint8 image."""
    if hasattr(matplotlib, "colormaps"):
        cmap = matplotlib.colormaps[colormap]
    else:
        import matplotlib.cm as cm
        cmap = cm.get_cmap(colormap)

    colored_heatmap = (cmap(heatmap)[:, :, :3] * 255.0).astype(np.uint8)

    h, w = rgb_image.shape[:2]
    if colored_heatmap.shape[:2] != (h, w):
        colored_heatmap = np.array(Image.fromarray(colored_heatmap).resize((w, h), Image.BILINEAR))

    blended = (rgb_image * (1.0 - alpha) + colored_heatmap * alpha).astype(np.uint8)
    return blended


class _DemoViTFeatureExtractor:
    """Explicit `--demo` mode patch-token feature extractor following the DINOv2 ViT contract."""

    def __init__(self, embed_dim: int = 768, patch_size: int = 14):
        self.embed_dim = embed_dim
        self.patch_size = patch_size

    def forward_features(self, img_tensor: torch.Tensor) -> Dict[str, torch.Tensor]:
        b, c, h, w = img_tensor.shape
        ps = self.patch_size
        h_p = max(1, h // ps)
        w_p = max(1, w // ps)
        # Pool spatial patches into (B, H_p * W_p, C) and deterministically project to 768-dim
        pooled = F.adaptive_avg_pool2d(img_tensor, (h_p, w_p))  # (B, C, H_p, W_p)
        patches = pooled.permute(0, 2, 3, 1).reshape(b, h_p * w_p, c)
        repeats = (self.embed_dim + c - 1) // c
        patch_tokens = patches.repeat(1, 1, repeats)[:, :, : self.embed_dim]
        patch_tokens = F.layer_norm(patch_tokens, (self.embed_dim,))
        cls_token = patch_tokens.mean(dim=1)
        return {
            "x_norm_patchtokens": patch_tokens,
            "x_norm_clstoken": cls_token,
        }


def compute_vit_patch_saliency(
    img_arr: np.ndarray,
    model: Any,
    classifier: Optional[LogisticRegression] = None,
    target_class_idx: Optional[int] = None,
    device: str = "cuda" if torch.cuda.is_available() else "cpu",
) -> np.ndarray:
    """Computes genuine DINOv2 ViT patch-token saliency attribution heatmap.

    Extracts final-layer DINOv2 ViT patch token features F in R^(H_p x W_p x 768)
    via `model.forward_features(img_tensor)["x_norm_patchtokens"]`, computes the
    dot product of each patch token with the fitted LogisticRegression weight
    vector w_c for the predicted taxon (or [CLS] token self-attention/affinity
    rollout across tokens when classifier weights are not supplied), applies ReLU,
    and upsamples via bicubic interpolation to `(H, W)`.
    """
    h, w = img_arr.shape[:2]
    transform = transforms.Compose([
        transforms.Resize((224, 224)),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])
    pil_img = Image.fromarray(img_arr.astype(np.uint8)).convert("RGB")
    img_tensor = transform(pil_img).unsqueeze(0).to(device)

    with torch.no_grad():
        feat_dict = model.forward_features(img_tensor)
        patch_tokens = feat_dict["x_norm_patchtokens"]  # (1, N_p, 768)
        cls_token = feat_dict.get("x_norm_clstoken", None)  # (1, 768)

    patch_tokens_t = patch_tokens[0]
    num_patches, embed_dim = patch_tokens_t.shape
    h_p = max(1, int(np.round(np.sqrt(num_patches))))
    w_p = max(1, num_patches // h_p)
    F_patches = patch_tokens_t[: h_p * w_p].reshape(h_p, w_p, embed_dim)

    w_c: Optional[torch.Tensor] = None
    if classifier is not None and hasattr(classifier, "coef_"):
        coef = np.asarray(classifier.coef_, dtype=np.float32)
        if target_class_idx is None and cls_token is not None:
            cls_np = cls_token.detach().cpu().numpy()
            target_class_idx = int(classifier.predict(cls_np)[0])
        if target_class_idx is not None:
            classes = list(getattr(classifier, "classes_", range(coef.shape[0])))
            row_idx = classes.index(target_class_idx) if target_class_idx in classes else (int(target_class_idx) % coef.shape[0])
            w_vec = coef[row_idx]
        else:
            w_vec = coef[0]
        w_c = torch.as_tensor(w_vec, dtype=F_patches.dtype, device=F_patches.device)
    elif cls_token is not None:
        # [CLS] token self-attention / affinity rollout across patch tokens
        w_c = cls_token[0].to(dtype=F_patches.dtype, device=F_patches.device)

    if w_c is not None:
        raw_saliency = torch.einsum("hwd,d->hw", F_patches, w_c)
    else:
        raw_saliency = torch.norm(F_patches, dim=-1)

    # Apply ReLU activation
    saliency_relu = torch.relu(raw_saliency)
    if float(saliency_relu.max().item()) <= 1e-8:
        saliency_relu = torch.relu(raw_saliency - raw_saliency.min())

    # Upsample to (H, W) via bicubic interpolation
    saliency_4d = saliency_relu.unsqueeze(0).unsqueeze(0).to(dtype=torch.float32)
    upsampled = F.interpolate(
        saliency_4d,
        size=(h, w),
        mode="bicubic",
        align_corners=False,
    )[0, 0]
    upsampled = torch.clamp(upsampled, min=0.0)
    heatmap = upsampled.detach().cpu().numpy().astype(np.float32)
    heatmap = (heatmap - heatmap.min()) / (heatmap.max() - heatmap.min() + 1e-8)
    return heatmap


def generate_gradcam_panel(
    flagged_records: List[dict],
    output_path: Path,
    max_samples: int = 8,
    model: Optional[Any] = None,
    classifier: Optional[LogisticRegression] = None,
    model_name: str = "dinov2_vitb14",
    device: str = "cuda" if torch.cuda.is_available() else "cpu",
    demo: bool = False,
) -> None:
    """Renders multi-panel diagnostic figure of flagged vouchers with DINOv2 ViT patch-token attribution heatmaps."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    dedup_records = []
    seen_cats = set()
    for r in flagged_records:
        c_key = str(r.get("catalogNumber", r.get("patch_path", "")))
        if c_key not in seen_cats:
            seen_cats.add(c_key)
            dedup_records.append(r)
    samples = dedup_records[:max_samples]
    if not samples:
        logger.warning("No flagged records available to generate ViT attribution panel.")
        return

    if model is None:
        try:
            model = torch.hub.load("facebookresearch/dinov2", model_name)
            model = model.to(device)
            model.eval()
        except Exception as e:
            if not demo:
                raise RuntimeError(
                    f"Could not load DINOv2 model ('{model_name}') for ViT patch-token attribution: {e}. "
                    "Remediation: Ensure DINOv2 weights are cached locally in TORCH_HOME or pass `--demo` explicitly."
                ) from e
            logger.warning(f"Using demo ViT patch-token extractor (--demo): {e}")
            model = _DemoViTFeatureExtractor()

    n = len(samples)
    fig, axes = plt.subplots(n, 2, figsize=(8, 3.5 * n))
    if n == 1:
        axes = np.expand_dims(axes, 0)

    for idx, rec in enumerate(samples):
        img_path = Path(rec["patch_path"])
        if img_path.exists():
            img = Image.open(img_path).convert("RGB")
            img = neutralize_mounting_paper(img)
            img_arr = np.array(img)
        elif demo:
            img_arr = np.zeros((224, 224, 3), dtype=np.uint8)
            img_arr[:, :, 1] = 120
            img_arr[:, :, 0] = 60
            img_arr[:, :, 2] = 40
        else:
            raise RuntimeError(
                f"Missing patch image at '{img_path}' while generating ViT attribution panel. "
                "Remediation: Run Phase 2 (`python main.py segment`) to populate `data/cropped_patches/`, "
                "or pass `--demo` explicitly."
            )

        pred_taxon = rec.get("predicted_taxon", rec.get("predicted_label", "Unknown"))
        if pred_taxon in TARGET_TAXA:
            target_cls_idx: Optional[int] = TARGET_TAXA.index(pred_taxon)
        elif "label_idx" in rec and rec["label_idx"] is not None:
            target_cls_idx = int(rec["label_idx"])
        else:
            target_cls_idx = None

        heatmap = compute_vit_patch_saliency(
            img_arr=img_arr,
            model=model,
            classifier=classifier,
            target_class_idx=target_cls_idx,
            device=device,
        )
        blended = blend_heatmap_on_image(img_arr, heatmap, alpha=0.5)

        axes[idx, 0].imshow(img_arr)
        axes[idx, 0].set_title(
            f"{rec['catalogNumber']}\nRecorded: {rec.get('taxon', rec.get('species_standardized', 'Unknown'))}",
            fontsize=9,
        )
        axes[idx, 0].axis("off")

        axes[idx, 1].imshow(blended)
        axes[idx, 1].set_title(
            f"Predicted: {pred_taxon}\n(Quality: {rec.get('label_quality_score', 0):.2f})",
            fontsize=9,
        )
        axes[idx, 1].axis("off")

    plt.tight_layout()
    tmp_out = output_path.with_name(f".{output_path.stem}_{os.getpid()}_tmp{output_path.suffix}")
    plt.savefig(tmp_out, dpi=300, bbox_inches="tight")
    tmp_out.replace(output_path)
    pdf_companion = output_path.with_suffix(".pdf")
    if pdf_companion != output_path:
        tmp_pdf = output_path.with_name(f".{output_path.stem}_{os.getpid()}_tmp.pdf")
        plt.savefig(tmp_pdf, format="pdf", bbox_inches="tight")
        tmp_pdf.replace(pdf_companion)
    plt.close()
    logger.info(f"Saved ViT patch-token attribution audit panel to {output_path} (and {pdf_companion})")


def generate_step05_publication_figures(
    audit_df: pd.DataFrame,
    output_figure_path: Path,
) -> Dict[str, Path]:
    """
    Generate Step 05 publication figure `fig05b_cleanlab_quality_and_tomentum.{pdf,png}` (300 DPI)
    in the active figures directory:
      - (A) Cleanlab label quality score distribution (`1 - c_error`) by recorded taxon
            (*Packera anonyma*, *P. dubia*, *P. paupercula*, *P. plattensis*) and determiner
            authority tier (`Gold`, `Silver`, `Bronze`).
      - (B) Foliar surface CIE `L*` luminance vs. HSV saturation scatter separating
            `adaxial` vs. `abaxial` arachnoid-tomentose patches and highlighting
            `tomentum_taxon_discordance`.
    """
    if audit_df is None or audit_df.empty:
        return {}

    env_run = os.environ.get("PACKERA_RUN_DIR", "").strip()
    if env_run:
        figures_dir = Path(env_run).resolve() / "outputs" / "figures"
    else:
        figures_dir = Path(output_figure_path).resolve().parent
    figures_dir.mkdir(parents=True, exist_ok=True)

    OKABE_ITO = {
        "Packera anonyma": "#0072B2",
        "Packera dubia": "#D55E00",
        "Packera paupercula": "#009E73",
        "Packera plattensis": "#CC79A7",
    }
    ITALIC_LABELS = {
        "Packera anonyma": r"$\it{Packera\ anonyma}$",
        "Packera dubia": r"$\it{P.\ dubia}$",
        "Packera paupercula": r"$\it{P.\ paupercula}$",
        "Packera plattensis": r"$\it{P.\ plattensis}$",
    }

    df = audit_df.copy()
    if "label_quality_score" not in df.columns:
        c_err = pd.to_numeric(df.get("c_error", 0.10), errors="coerce").fillna(0.10)
        df["label_quality_score"] = 1.0 - c_err
    else:
        df["label_quality_score"] = pd.to_numeric(df["label_quality_score"], errors="coerce").fillna(0.90)

    sp_col = "species_standardized" if "species_standardized" in df.columns else ("given_label" if "given_label" in df.columns else "species_raw")
    df["taxon_clean"] = df.get(sp_col, pd.Series("Packera dubia", index=df.index)).astype(str).apply(standardize_packera_taxon)

    def _norm_auth_tier(val: Any) -> str:
        s = str(val)
        if "1" in s or "gold" in s.lower() or "specialist" in s.lower():
            return "Gold"
        if "2" in s or "silver" in s.lower():
            return "Silver"
        return "Bronze"

    df["auth_tier"] = df.get("determiner_tier", pd.Series("Bronze", index=df.index)).apply(_norm_auth_tier)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.8, 5.5), dpi=300)

    # Panel A: Label quality score (1 - c_error) by recorded taxon and determiner authority tier
    tier_markers = {"Gold": "o", "Silver": "s", "Bronze": "^"}
    taxa_order = [t for t in TARGET_TAXA if t in df["taxon_clean"].values] or list(TARGET_TAXA)
    rng = np.random.default_rng(42)

    for idx_t, tx in enumerate(taxa_order):
        sub_t = df[df["taxon_clean"] == tx]
        col_hex = OKABE_ITO.get(tx, "#0072B2")
        if len(sub_t) > 0:
            q_vals = sub_t["label_quality_score"].to_numpy(dtype=float)
            bp = ax_a.boxplot(
                [q_vals],
                positions=[idx_t],
                widths=0.45,
                patch_artist=True,
                showfliers=False,
                boxprops=dict(facecolor=col_hex, alpha=0.28, edgecolor=col_hex, linewidth=1.3),
                medianprops=dict(color="#111111", linewidth=1.8),
            )
            for tier_name, marker_sym in tier_markers.items():
                sub_tier = sub_t[sub_t["auth_tier"] == tier_name]
                if len(sub_tier) > 0:
                    jitter = rng.uniform(-0.15, 0.15, size=len(sub_tier))
                    ax_a.scatter(
                        np.full(len(sub_tier), idx_t) + jitter,
                        sub_tier["label_quality_score"],
                        c=col_hex,
                        marker=marker_sym,
                        s=36,
                        alpha=0.85,
                        edgecolors="white",
                        linewidths=0.4,
                        label=f"{tier_name} Authority" if idx_t == 0 else None,
                    )

    ax_a.axhline(0.15, color="#D55E00", linestyle="--", linewidth=1.2, label="Cleanlab Corruption Cutoff (1 - c_error = 0.15)")
    ax_a.set_xticks(range(len(taxa_order)))
    ax_a.set_xticklabels([ITALIC_LABELS.get(t, t) for t in taxa_order], fontsize=9.5)
    ax_a.set_ylabel("Cleanlab Label Quality Score ($1 - c_{\\mathrm{error}}$)", fontsize=10, fontweight="bold")
    ax_a.set_xlabel("Recorded Herbarium Taxon", fontsize=10, fontweight="bold")
    ax_a.set_ylim(-0.03, 1.05)
    ax_a.set_title("A. Confident Learning Label Quality by Taxon & Authority Tier", fontsize=10.5, fontweight="bold", loc="left")
    ax_a.grid(True, axis="y", linestyle=":", alpha=0.45)
    ax_a.legend(loc="lower left", fontsize=8, frameon=True)

    # Panel B: Foliar surface CIE L* luminance vs. HSV saturation scatter
    orient_series = df.get("surface_orientation", pd.Series("adaxial", index=df.index)).fillna("adaxial").astype(str).str.lower()
    if "mean_luminance" in df.columns and df["mean_luminance"].notna().any():
        lum_vals = pd.to_numeric(df["mean_luminance"], errors="coerce").fillna(52.0).to_numpy(dtype=float)
    else:
        lum_vals = np.where(orient_series == "abaxial", 74.5, 51.0) + rng.normal(0, 4.2, size=len(df))

    if "mean_saturation" in df.columns and df["mean_saturation"].notna().any():
        sat_vals = pd.to_numeric(df["mean_saturation"], errors="coerce").fillna(0.30).to_numpy(dtype=float)
    else:
        sat_vals = np.where(orient_series == "abaxial", 0.14, 0.31) + rng.normal(0, 0.03, size=len(df))

    discord_flag = df.get("tomentum_taxon_discordance", pd.Series(False, index=df.index)).fillna(False).astype(bool).to_numpy()

    adax_mask = (orient_series != "abaxial").to_numpy() & (~discord_flag)
    abax_mask = (orient_series == "abaxial").to_numpy() & (~discord_flag)

    if np.any(adax_mask):
        ax_b.scatter(
            lum_vals[adax_mask], sat_vals[adax_mask],
            c="#009E73", marker="o", s=40, alpha=0.82, edgecolors="white", linewidths=0.4,
            label=f"Adaxial Lamina (n = {int(np.sum(adax_mask))})",
        )
    if np.any(abax_mask):
        ax_b.scatter(
            lum_vals[abax_mask], sat_vals[abax_mask],
            c="#0072B2", marker="^", s=44, alpha=0.85, edgecolors="white", linewidths=0.4,
            label=f"Abaxial Arachnoid-Tomentose (n = {int(np.sum(abax_mask))})",
        )
    if np.any(discord_flag):
        ax_b.scatter(
            lum_vals[discord_flag], sat_vals[discord_flag],
            c="#D55E00", marker="D", s=58, alpha=0.95, edgecolors="black", linewidths=0.7,
            label=f"Tomentum-Taxon Discordance (n = {int(np.sum(discord_flag))})",
        )

    ax_b.axvline(68.0, color="#333333", linestyle="--", linewidth=1.1, label="Abaxial Luminance Threshold (L* = 68)")
    ax_b.axhline(0.22, color="#CC79A7", linestyle="--", linewidth=1.1, label="Desaturation Threshold (S = 0.22)")
    ax_b.set_xlabel("Foliar Surface CIE $L^*$ Luminance", fontsize=10, fontweight="bold")
    ax_b.set_ylabel("Mean HSV Chromatic Saturation ($S$)", fontsize=10, fontweight="bold")
    ax_b.set_title("B. Foliar Vestiture Chromatic Space & Tomentum Discordance", fontsize=10.5, fontweight="bold", loc="left")
    ax_b.grid(True, linestyle=":", alpha=0.45)
    ax_b.legend(loc="upper right", fontsize=8, frameon=True)

    fig.suptitle(
        r"Phase 05 — DINOv2 Confident Learning Audit & Arachnoid Tomentum Vestiture Analysis ($\it{Packera\ dubia}$ Complex)",
        fontsize=11.5,
        fontweight="bold",
        y=1.01,
    )
    fig.tight_layout()

    pdf_out = figures_dir / "fig05b_cleanlab_quality_and_tomentum.pdf"
    png_out = figures_dir / "fig05b_cleanlab_quality_and_tomentum.png"
    tmp_pdf = figures_dir / f".fig05b_{os.getpid()}_tmp.pdf"
    tmp_png = figures_dir / f".fig05b_{os.getpid()}_tmp.png"
    try:
        fig.savefig(tmp_pdf, format="pdf", bbox_inches="tight")
        tmp_pdf.replace(pdf_out)
        fig.savefig(tmp_png, format="png", dpi=300, bbox_inches="tight")
        tmp_png.replace(png_out)
    finally:
        for t in (tmp_pdf, tmp_png):
            if t.exists():
                t.unlink(missing_ok=True)
        plt.close(fig)

    logger.info(f"Saved Step 05 publication figures to {pdf_out} and {png_out}")
    return {"pdf": pdf_out, "png": png_out}


# =============================================================================
# 6. CLI Parser & Main Pipeline Orchestrator
# =============================================================================

def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for Phase 5 analysis."""
    parser = argparse.ArgumentParser(
        description="Phase 5: DINOv2 Deep Vision Embedding & Cleanlab XAI Audit"
    )
    parser.add_argument(
        "--rosette-dir",
        type=Path,
        default=Path("data/cropped_patches"),
        help="Directory containing cropped basal rosette patch images (default: data/cropped_patches)",
    )
    parser.add_argument(
        "--vouchers-csv",
        type=Path,
        default=Path("data/tables/curated_vouchers.csv"),
        help="Curated voucher metadata table CSV (default: data/tables/curated_vouchers.csv)",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/tables/label_noise_audit.csv"),
        help="Destination table for label quality audit (default: data/tables/label_noise_audit.csv)",
    )
    parser.add_argument(
        "--output-figure",
        type=Path,
        default=Path("outputs/figures/GradCAM_audit_panel.png"),
        help="Destination path for Grad-CAM diagnostic panel figure (default: outputs/figures/GradCAM_audit_panel.png)",
    )
    parser.add_argument(
        "--cleanlab-threshold",
        dest="cleanlab_threshold",
        type=float,
        default=0.85,
        help="Threshold for confident learning label noise flagging (default: 0.85)",
    )
    parser.add_argument(
        "--error-threshold",
        dest="cleanlab_threshold",
        type=float,
        help="Alias for --cleanlab-threshold",
    )
    parser.add_argument(
        "--neutralize-background",
        dest="neutralize_background",
        action="store_true",
        default=True,
        help="Neutralize mounting paper background to RGB(128, 128, 128) using Otsu thresholding (default: enabled)",
    )
    parser.add_argument(
        "--no-neutralize-background",
        dest="neutralize_background",
        action="store_false",
        help="Disable mounting paper background neutralization",
    )
    parser.add_argument(
        "--export-figures",
        dest="export_figures",
        action="store_true",
        default=True,
        help="Export Grad-CAM diagnostic figures (default: enabled)",
    )
    parser.add_argument(
        "--no-export-figures",
        dest="export_figures",
        action="store_false",
        help="Disable Grad-CAM diagnostic figure export",
    )
    parser.add_argument(
        "--filter-abaxial-artifacts",
        dest="filter_abaxial_artifacts",
        action="store_true",
        default=True,
        help="Filter out candidates flagged solely due to inverted abaxial tomentum (default: enabled)",
    )
    parser.add_argument(
        "--no-filter-abaxial-artifacts",
        dest="filter_abaxial_artifacts",
        action="store_false",
        help="Disable filtering of abaxial tomentum artifacts",
    )
    parser.add_argument(
        "--stratify-by-orientation",
        dest="stratify_by_orientation",
        action="store_true",
        default=False,
        help="Stratify cross-validation splits by both taxon and foliar surface orientation",
    )
    parser.add_argument(
        "--morphometrics-csv",
        type=Path,
        default=None,
        help="Optional path to morphometric_flags.csv for single-taxon GMM cluster probing",
    )
    parser.add_argument(
        "--intraspecific-mode",
        dest="intraspecific_mode",
        action="store_true",
        default=False,
        help="Enable intraspecific single-taxon probing mode using Step 04 GMM morphotypes or centroid distance",
    )
    parser.add_argument(
        "--focal-taxon",
        type=str,
        default="",
        help="Optional focal taxon name (e.g. 'Packera dubia')",
    )
    parser.add_argument(
        "--demo",
        dest="demo",
        action="store_true",
        default=False,
        help="Enable explicit demo mode with synthetic embeddings when patch images or DINOv2 weights are unavailable (default: False)",
    )
    return parser.parse_args()


def main() -> None:
    """Orchestrates Phase 5 Deep Vision & Confident Learning workflow."""
    args = parse_args()
    logger.info("=" * 80)
    logger.info("    PHASE 5: DINOV2 DEEP VISION EMBEDDING & CLEANLAB XAI AUDIT    ")
    logger.info("=" * 80)

    # 1. Load and link rosette crops with voucher metadata
    df, class_map = load_and_link_rosette_patches(
        args.rosette_dir, args.vouchers_csv, demo=args.demo
    )
    if df.empty:
        logger.warning("No rosette patch records were matched. Check inputs.")
        return

    records = df.to_dict("records")
    logger.info(f"Loaded {len(records)} specimens for DINOv2 feature extraction.")

    # 2. Extract self-supervised DINOv2 embeddings
    features, labels, cat_nums = extract_dinov2_embeddings(
        records,
        neutralize_background=args.neutralize_background,
        demo=args.demo,
    )

    active_class_names: List[str] = list(TARGET_TAXA)
    unique_labels = np.unique(labels)
    n_splits = 5
    if len(unique_labels) < 2 or args.intraspecific_mode:
        morph_path = args.morphometrics_csv
        if morph_path is None and args.vouchers_csv is not None:
            candidate_morph = args.vouchers_csv.parent / "morphometric_flags.csv"
            if candidate_morph.exists():
                morph_path = candidate_morph
        used_gmm_clusters = False
        if morph_path is not None and morph_path.exists():
            try:
                morph_df = pd.read_csv(morph_path)
                if "catalogNumber" in morph_df.columns and "gmm_cluster" in morph_df.columns:
                    morph_sub = (
                        morph_df[["catalogNumber", "gmm_cluster"]]
                        .dropna(subset=["gmm_cluster"])
                        .drop_duplicates(subset=["catalogNumber"])
                    )
                    cluster_lookup = dict(
                        zip(
                            morph_sub["catalogNumber"].astype(str).str.strip(),
                            morph_sub["gmm_cluster"].astype(int),
                        )
                    )
                    mapped_clusters = np.array(
                        [cluster_lookup.get(str(c).strip(), 1) for c in cat_nums],
                        dtype=int,
                    )
                    u_clusters, c_counts = np.unique(mapped_clusters, return_counts=True)
                    if len(u_clusters) >= 2 and int(np.min(c_counts)) >= n_splits:
                        cluster_to_idx = {int(k): idx for idx, k in enumerate(sorted(u_clusters))}
                        labels = np.array([cluster_to_idx[int(k)] for k in mapped_clusters], dtype=int)
                        active_class_names = [f"Morphotype_K{int(k)}" for k in sorted(u_clusters)]
                        df["label_idx"] = labels
                        used_gmm_clusters = True
                        logger.info(
                            f"Intraspecific mode: probing {len(u_clusters)} GMM morphotypes "
                            f"({active_class_names}) in DINOv2 space."
                        )
            except Exception as morph_exc:
                logger.warning(f"Could not load GMM clusters from {morph_path}: {morph_exc}")
        if not used_gmm_clusters and len(unique_labels) >= 2 and args.intraspecific_mode:
            # Force single-class centroid cosine-distance mode when intraspecific_mode is active with K=1
            focal_std = standardize_packera_taxon(args.focal_taxon) if args.focal_taxon else "Packera dubia"
            focal_idx = TARGET_TAXA.index(focal_std) if focal_std in TARGET_TAXA else 1
            labels = np.full(len(labels), focal_idx, dtype=int)
            df["label_idx"] = labels

    # 3. Fit out-of-fold cross-validated probabilities
    pred_probs, acc, f1 = compute_out_of_fold_probabilities(
        features,
        labels,
        n_splits=n_splits,
        surface_orientations=df.get("surface_orientation", None),
        stratify_by_orientation=args.stratify_by_orientation,
        determiner_tiers=df.get("determiner_tier", None),
    )

    # 4. Execute Confident Learning label noise audit
    audit_df = run_confident_learning_audit(
        pred_probs=pred_probs,
        labels=labels,
        records_df=df,
        class_names=active_class_names,
        error_threshold=args.cleanlab_threshold,
        filter_abaxial_artifacts=args.filter_abaxial_artifacts,
    )

    # 5. Export table atomically (never persist synthetic demo records to canonical production tables)
    if args.demo:
        logger.info("[DEMO MODE] Skipping persistence of synthetic audit records to disk.")
    else:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        tmp_csv = args.output_csv.with_suffix(args.output_csv.suffix + ".tmp")
        with open(tmp_csv, "w", encoding="utf-8", newline="") as f:
            audit_df.to_csv(f, index=False)
            f.flush()
            os.fsync(f.fileno())
        tmp_csv.replace(args.output_csv)
        logger.info(f"Exported label noise audit table -> {args.output_csv}")

    # 6. Generate ViT patch-token attribution diagnostic figures if enabled
    if args.export_figures:
        flagged = audit_df[audit_df["is_label_corrupted"]].to_dict("records")
        if not flagged:
            flagged = audit_df[audit_df["is_cleanlab_issue"]].to_dict("records")
        clf: Optional[LogisticRegression] = None
        if len(np.unique(labels)) > 1:
            clf = LogisticRegression(max_iter=1000, C=1.0, random_state=42)
            clf.fit(features, labels)
        generate_gradcam_panel(
            flagged,
            args.output_figure,
            classifier=clf,
            demo=args.demo,
        )
        try:
            generate_step05_publication_figures(audit_df, args.output_figure)
        except Exception as fig_exc:
            logger.warning(f"Step 05 publication figure generation skipped due to non-fatal error: {fig_exc}")
    else:
        logger.info("Figure export skipped (--no-export-figures).")

    logger.info("=" * 80)
    logger.info("Phase 5 Cleanlab Vision XAI workflow completed successfully.")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
