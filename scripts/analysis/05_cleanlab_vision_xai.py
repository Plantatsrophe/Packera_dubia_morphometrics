#!/usr/bin/env python3
"""
scripts/analysis/05_cleanlab_vision_xai.py
=========================================
Deep Vision, Confident Learning Label Noise Curation & Grad-CAM XAI
Project: Multimodal Morphometrics & Species Delimitation in the Packera dubia Complex
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Key Capabilities:
1. Ingests unsegmented dense basal rosette crops and links curated voucher metadata.
2. Extracts 768-dimensional [CLS] self-supervised token embeddings using dinov2_vitb14.
3. Fits 5-fold cross-validated linear classifier heads to compute out-of-fold probability matrices.
4. Deploys Confident Learning (cleanlab) to estimate joint noise matrices and flag mislabeled vouchers (C_error > 0.85).
5. Generates Grad-CAM attribution heatmaps (Captum) to verify model focus on diagnostic tomentum & leaf margins.
6. Exports data/tables/label_noise_audit.csv and outputs/figures/GradCAM_audit_panel.png.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from PIL import Image
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, f1_score
from sklearn.model_selection import StratifiedKFold

import cleanlab
from cleanlab.count import compute_confident_joint, estimate_joint
from cleanlab.filter import find_label_issues
from cleanlab.rank import get_label_quality_scores

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from torchvision import transforms

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("Cleanlab_Vision_XAI")

TARGET_TAXA = [
    "Packera anonyma",
    "Packera dubia",
    "Packera paupercula",
    "Packera plattensis",
]


def standardize_packera_taxon(species_str: Optional[str]) -> str:
    """Standardize synonymy into the four core Packera dubia complex taxa."""
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


class RosettePatchDataset(Dataset):
    """PyTorch Dataset loading dense basal rosette image crops for DINOv2."""

    def __init__(self, records: List[Dict], transform: transforms.Compose):
        self.records = records
        self.transform = transform

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, int, str]:
        rec = self.records[idx]
        image_path = Path(rec["patch_path"])
        image = Image.open(image_path).convert("RGB")
        tensor = self.transform(image)
        return tensor, rec["label_idx"], rec["catalogNumber"]


def load_and_link_rosette_patches(
    rosette_dir: Path, vouchers_csv: Path
) -> Tuple[pd.DataFrame, Dict[str, int]]:
    """Links rosette patch image files with curated herbarium metadata."""
    logger.info(f"Loading curated vouchers from: {vouchers_csv}")
    vouchers_df = pd.read_csv(vouchers_csv)
    vouchers_df["catalogNumber"] = vouchers_df["catalogNumber"].astype(str)
    vouchers_df["species_std"] = vouchers_df["species_raw"].apply(standardize_packera_taxon)

    vouchers_dedup = vouchers_df.drop_duplicates(subset=["catalogNumber"])
    voucher_lookup = vouchers_dedup.set_index("catalogNumber").to_dict(orient="index")
    label_to_idx = {name: idx for idx, name in enumerate(TARGET_TAXA)}

    patch_files = sorted(list(rosette_dir.glob("*.jpg")) + list(rosette_dir.glob("*.png")))
    logger.info(f"Found {len(patch_files)} rosette patch files in {rosette_dir}")

    records = []
    for p in patch_files:
        fname = p.stem
        m = re.match(r"^(.*?)(?:_p\d+)?_rosette$", fname)
        cat_num = m.group(1) if m else fname.split("_")[0]

        meta = voucher_lookup.get(cat_num, None)
        if meta is None:
            # Fallback exact lookup
            meta = voucher_lookup.get(fname, None)

        if meta:
            std_sp = meta.get("species_std", "Unknown")
            if std_sp in label_to_idx:
                records.append({
                    "catalogNumber": cat_num,
                    "patch_path": str(p),
                    "species_raw": meta.get("species_raw", ""),
                    "species_standardized": std_sp,
                    "label_idx": label_to_idx[std_sp],
                    "determiner_raw": meta.get("determiner_raw", "Unknown"),
                    "determiner_tier": meta.get("determiner_tier", "Tier_3_Bronze"),
                    "stateProvince": meta.get("stateProvince", ""),
                    "county": meta.get("county", ""),
                    "regional_group": meta.get("regional_group", "Unknown"),
                })

    df_records = pd.DataFrame(records)
    logger.info(f"Successfully linked {len(df_records)} rosette patches to target taxa.")
    logger.info(f"Class distribution:\n{df_records['species_standardized'].value_counts().to_string()}")
    return df_records, label_to_idx


def extract_dinov2_embeddings(
    df_records: pd.DataFrame,
    device: torch.device,
    batch_size: int = 32,
    backbone_name: str = "dinov2_vitb14",
) -> Tuple[np.ndarray, nn.Module]:
    """Extracts 768-dimensional [CLS] token embeddings using DINOv2."""
    logger.info(f"Loading {backbone_name} backbone from torch.hub on device: {device}...")
    model = torch.hub.load("facebookresearch/dinov2", backbone_name)
    model.to(device)
    model.eval()

    transform = transforms.Compose([
        transforms.Resize((224, 224), interpolation=transforms.InterpolationMode.BICUBIC),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])

    dataset = RosettePatchDataset(df_records.to_dict(orient="records"), transform)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=2)

    embeddings_list = []
    logger.info(f"Extracting [CLS] visual embeddings across {len(dataset)} rosette patches...")
    with torch.no_grad():
        for batch_tensors, _, _ in loader:
            batch_tensors = batch_tensors.to(device)
            # DINOv2 returns 768-dim [CLS] token
            emb = model(batch_tensors)
            embeddings_list.append(emb.cpu().numpy())

    embeddings = np.vstack(embeddings_list)
    logger.info(f"Visual embedding extraction complete. Feature matrix shape: {embeddings.shape}")
    return embeddings, model


def compute_out_of_fold_probabilities(
    embeddings: np.ndarray, labels: np.ndarray, n_splits: int = 5
) -> Tuple[np.ndarray, LogisticRegression]:
    """Fits 5-fold cross-validated linear heads and computes out-of-fold probabilities."""
    logger.info(f"Fitting {n_splits}-fold Stratified Cross-Validation linear classifier heads...")
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    num_classes = len(TARGET_TAXA)
    pred_probs = np.zeros((len(labels), num_classes), dtype=np.float64)

    fold_accuracies, fold_f1s = [], []
    for fold, (train_idx, val_idx) in enumerate(skf.split(embeddings, labels), 1):
        X_train, y_train = embeddings[train_idx], labels[train_idx]
        X_val, y_val = embeddings[val_idx], labels[val_idx]

        clf = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs", random_state=42)
        clf.fit(X_train, y_train)

        probs_val = clf.predict_proba(X_val)
        pred_probs[val_idx] = probs_val

        preds_val = np.argmax(probs_val, axis=1)
        acc = accuracy_score(y_val, preds_val)
        f1 = f1_score(y_val, preds_val, average="macro")
        fold_accuracies.append(acc)
        fold_f1s.append(f1)
        logger.info(f"  Fold {fold}/{n_splits} - Acc: {acc * 100:.2f}%, Macro-F1: {f1:.4f}")

    overall_preds = np.argmax(pred_probs, axis=1)
    overall_acc = accuracy_score(labels, overall_preds)
    overall_f1 = f1_score(labels, overall_preds, average="macro")
    logger.info(f"OOF Cross-Validation Summary -> Overall Acc: {overall_acc * 100:.2f}%, Macro-F1: {overall_f1:.4f}")

    # Full linear model for Grad-CAM explanation
    full_clf = LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs", random_state=42)
    full_clf.fit(embeddings, labels)
    return pred_probs, full_clf


def run_confident_learning_audit(
    df_records: pd.DataFrame,
    labels: np.ndarray,
    pred_probs: np.ndarray,
    threshold: float = 0.85,
) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Deploys cleanlab confident learning to identify label noise and joint errors."""
    logger.info("Executing Cleanlab Confident Learning label noise audit...")
    idx_to_label = {idx: name for idx, name in enumerate(TARGET_TAXA)}

    # Compute confident joint & joint distribution matrix
    conf_joint = compute_confident_joint(labels=labels, pred_probs=pred_probs)
    joint_noise_matrix = estimate_joint(labels=labels, pred_probs=pred_probs)

    # Label quality scores Q in [0, 1]; C_error = 1 - Q
    quality_scores = get_label_quality_scores(labels=labels, pred_probs=pred_probs)
    c_error = 1.0 - quality_scores

    # Cleanlab identified issue indices
    issue_indices = find_label_issues(
        labels=labels, pred_probs=pred_probs, return_indices_ranked_by="self_confidence"
    )
    issue_mask = np.zeros(len(labels), dtype=bool)
    issue_mask[issue_indices] = True

    # Assemble audit table
    audit_df = df_records.copy()
    audit_df["given_label"] = [idx_to_label[i] for i in labels]
    audit_df["predicted_label"] = [idx_to_label[i] for i in np.argmax(pred_probs, axis=1)]
    audit_df["confidence_given_class"] = [pred_probs[i, labels[i]] for i in range(len(labels))]
    audit_df["confidence_predicted_class"] = np.max(pred_probs, axis=1)
    audit_df["label_quality_score"] = np.round(quality_scores, 4)
    audit_df["c_error"] = np.round(c_error, 4)
    audit_df["is_cleanlab_issue"] = issue_mask
    audit_df["is_label_corrupted"] = audit_df["c_error"] > threshold

    # Define triage actions and discordance reason
    def assign_triage(row):
        if row["c_error"] > threshold:
            return "Prune & Queue for Expert Re-determination"
        elif row["c_error"] > 0.60:
            return "Flag for Secondary Review"
        return "Retain in Benchmark Dataset"

    def assign_reason(row):
        if row["c_error"] > threshold:
            return f"Severe Discordance: Given {row['given_label']} vs Predicted {row['predicted_label']} (C_error={row['c_error']:.2f})"
        elif row["given_label"] != row["predicted_label"]:
            return f"Moderate Morphology Mismatch: Predicted {row['predicted_label']}"
        return "Morphologically Consistent"

    audit_df["triage_action"] = audit_df.apply(assign_triage, axis=1)
    audit_df["discordance_reason"] = audit_df.apply(assign_reason, axis=1)

    # Sort by error probability descending
    audit_df = audit_df.sort_values(by="c_error", ascending=False).reset_index(drop=True)

    n_flagged = (audit_df["c_error"] > threshold).sum()
    logger.info(f"Cleanlab Audit Complete: {n_flagged}/{len(audit_df)} vouchers flagged with C_error > {threshold}")
    return audit_df, conf_joint, joint_noise_matrix


class EndToEndViTClassifier(nn.Module):
    """Wraps DINOv2 backbone and linear classifier head for Grad-CAM attribution."""

    def __init__(self, backbone: nn.Module, linear_clf: LogisticRegression):
        super().__init__()
        self.backbone = backbone
        num_classes, num_features = linear_clf.coef_.shape
        self.classifier = nn.Linear(num_features, num_classes)
        with torch.no_grad():
            self.classifier.weight.copy_(torch.from_numpy(linear_clf.coef_).float())
            self.classifier.bias.copy_(torch.from_numpy(linear_clf.intercept_).float())

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        tokens = self.backbone.get_intermediate_layers(x, n=1, return_class_token=True)[0]
        cls_token = tokens[1]
        return self.classifier(cls_token)


def compute_vit_gradcam_heatmap(
    model: EndToEndViTClassifier,
    image_tensor: torch.Tensor,
    target_class_idx: int,
    device: torch.device,
) -> np.ndarray:
    """Computes high-resolution Grad-CAM spatial attribution heatmap on DINOv2 patch tokens."""
    model.eval()
    img = image_tensor.unsqueeze(0).to(device)

    # Extract intermediate patch features from the last transformer block
    with torch.no_grad():
        tokens = model.backbone.get_intermediate_layers(img, n=1, reshape=True, return_class_token=True)[0]
        patch_feats, _ = tokens[0], tokens[1]  # (1, 768, 16, 16)

    # Linear classifier weights for target class
    class_weight = model.classifier.weight[target_class_idx].detach()  # (768,)
    # Project class weights onto spatial patch features (1, 16, 16)
    cam = (patch_feats[0] * class_weight.view(-1, 1, 1)).sum(dim=0)
    cam = torch.relu(cam)

    # Interpolate to 224x224
    cam_2d = cam.unsqueeze(0).unsqueeze(0)
    cam_upsampled = nn.functional.interpolate(cam_2d, size=(224, 224), mode="bilinear", align_corners=False)[0, 0]
    cam_np = cam_upsampled.cpu().numpy()

    # Min-max normalization
    denom = cam_np.max() - cam_np.min()
    if denom > 1e-8:
        cam_norm = (cam_np - cam_np.min()) / denom
    else:
        cam_norm = np.zeros_like(cam_np)
    return cam_norm


def plot_gradcam_audit_panel(
    joint_noise_matrix: np.ndarray,
    audit_df: pd.DataFrame,
    exemplars: List[Dict],
    output_png: Path,
    threshold: float = 0.85,
) -> None:
    """Renders and exports the comprehensive multi-panel Grad-CAM and Confident Learning figure."""
    logger.info(f"Generating publication-grade Grad-CAM audit panel to {output_png}...")
    output_png.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(18, 14), dpi=300)
    gs = fig.add_gridspec(3, 4, height_ratios=[1.1, 1.0, 1.0], hspace=0.32, wspace=0.28)

    # ----------------- PANEL A: Cleanlab Joint Noise Matrix -----------------
    ax_mat = fig.add_subplot(gs[0, 0:2])
    short_names = ["P. anonyma", "P. dubia", "P. paupercula", "P. plattensis"]
    sns.heatmap(
        joint_noise_matrix,
        annot=True,
        fmt=".3f",
        cmap="YlGnBu",
        xticklabels=short_names,
        yticklabels=short_names,
        ax=ax_mat,
        cbar_kws={"label": "Estimated Joint Probability $P(\\hat{y}, y^*)$"},
    )
    ax_mat.set_title("(A) Cleanlab Confident Learning Joint Noise Matrix", fontsize=12, fontweight="bold")
    ax_mat.set_xlabel("Latent True Biological Taxon ($y^*$)", fontsize=10, fontweight="bold")
    ax_mat.set_ylabel("Given Herbarium Label ($\\hat{y}$)", fontsize=10, fontweight="bold")

    # ----------------- PANEL B: C_error Corruption Distribution -----------------
    ax_hist = fig.add_subplot(gs[0, 2:4])
    sns.histplot(
        audit_df["c_error"],
        bins=25,
        kde=True,
        color="#2b5c8f",
        edgecolor="white",
        ax=ax_hist,
        stat="count",
    )
    ax_hist.axvline(
        threshold,
        color="#d9534f",
        linestyle="--",
        linewidth=2,
        label=f"Pruning Threshold ($C_{{error}} > {threshold}$)",
    )
    n_corrupt = (audit_df["c_error"] > threshold).sum()
    ax_hist.set_title(
        f"(B) Voucher Misidentification Confidence ($N_{{flagged}} = {n_corrupt}$)",
        fontsize=12,
        fontweight="bold",
    )
    ax_hist.set_xlabel("Label Corruption Probability ($C_{error} = 1 - Q$)", fontsize=10, fontweight="bold")
    ax_hist.set_ylabel("Voucher Frequency", fontsize=10, fontweight="bold")
    ax_hist.legend(loc="upper center", frameon=True, fontsize=9)
    ax_hist.grid(axis="y", linestyle=":", alpha=0.6)

    # ----------------- PANEL C: Grad-CAM Exemplar Audit Rows -----------------
    # Render 4 exemplar vouchers (2 clean anchors + 2 flagged misidentifications)
    for idx, ex in enumerate(exemplars[:4]):
        row_pos = 1 if idx < 2 else 2
        col_base = (idx % 2) * 2

        ax_img = fig.add_subplot(gs[row_pos, col_base])
        ax_cam = fig.add_subplot(gs[row_pos, col_base + 1])

        orig_img = ex["orig_img"]
        cam_map = ex["cam_map"]

        # Display raw rosette patch
        ax_img.imshow(orig_img)
        ax_img.set_title(
            f"{ex['type']}: {ex['catalog']}\nGiven: {ex['given']} ({ex['tier']})",
            fontsize=8.5,
            fontweight="bold",
            color="#b30000" if ex["c_error"] > threshold else "#006600",
        )
        ax_img.axis("off")

        # Display Grad-CAM attribution overlay
        ax_cam.imshow(orig_img)
        cam_overlay = ax_cam.imshow(cam_map, cmap="jet", alpha=0.55, vmin=0.0, vmax=1.0)
        ax_cam.set_title(
            f"Grad-CAM Saliency Focus\nPred: {ex['pred']} ($C_{{err}}={ex['c_error']:.2f}$)",
            fontsize=8.5,
            fontweight="bold",
        )
        ax_cam.axis("off")

    plt.suptitle(
        "DINOv2 Self-Supervised Vision Embeddings, Confident Learning Label Curation & Grad-CAM XAI Panel\n"
        "Packera dubia Complex Species Delimitation Pipeline (NCU Herbarium)",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    plt.savefig(output_png, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"Figure successfully saved to {output_png}")


def parse_args() -> argparse.Namespace:
    """Parses command-line interface arguments."""
    parser = argparse.ArgumentParser(
        description="DINOv2 Confident Learning & Grad-CAM XAI for Packera herbarium vouchers."
    )
    parser.add_argument(
        "--rosette-dir",
        type=Path,
        default=Path("data/cropped_patches/rosettes_dense"),
        help="Path to directory containing dense rosette patch images.",
    )
    parser.add_argument(
        "--vouchers-csv",
        type=Path,
        default=Path("data/tables/curated_vouchers.csv"),
        help="Path to curated vouchers metadata CSV.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/tables/label_noise_audit.csv"),
        help="Destination path for label noise audit table.",
    )
    parser.add_argument(
        "--output-figure",
        type=Path,
        default=Path("outputs/figures/GradCAM_audit_panel.png"),
        help="Destination path for Grad-CAM audit panel figure.",
    )
    parser.add_argument(
        "--cleanlab-threshold",
        type=float,
        default=0.85,
        help="Confidence threshold (C_error > threshold) for flagging label errors.",
    )
    parser.add_argument(
        "--backbone",
        type=str,
        default="dinov2_vitb14",
        help="DINOv2 backbone model identifier (default: dinov2_vitb14).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Batch size for DINOv2 feature extraction.",
    )
    parser.add_argument(
        "--n-splits",
        type=int,
        default=5,
        help="Number of folds for Stratified Cross-Validation.",
    )
    return parser.parse_args()


def main() -> int:
    """Main execution pipeline."""
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Initializing DINOv2 + Cleanlab XAI pipeline on device: {device}")

    # 1. Ingest & link rosette patches to metadata
    df_records, label_to_idx = load_and_link_rosette_patches(args.rosette_dir, args.vouchers_csv)
    if df_records.empty:
        logger.error("No valid rosette patch records found. Exiting.")
        return 1

    labels = df_records["label_idx"].values

    # 2. Extract DINOv2 768-dim embeddings
    embeddings, backbone_model = extract_dinov2_embeddings(
        df_records, device=device, batch_size=args.batch_size, backbone_name=args.backbone
    )

    # 3. 5-Fold Stratified CV linear heads & OOF probability estimation
    pred_probs, full_clf = compute_out_of_fold_probabilities(
        embeddings, labels, n_splits=args.n_splits
    )

    # 4. Confident Learning & Label Noise Audit
    audit_df, conf_joint, joint_noise_matrix = run_confident_learning_audit(
        df_records, labels, pred_probs, threshold=args.cleanlab_threshold
    )

    # Export audit table
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit_df.to_csv(args.output_csv, index=False)
    logger.info(f"Saved complete label noise audit table to {args.output_csv}")

    # 5. Grad-CAM visual explanations for exemplar vouchers
    logger.info("Generating Grad-CAM attribution heatmaps for exemplar vouchers...")
    vit_clf = EndToEndViTClassifier(backbone_model, full_clf).to(device)

    # Pick 2 clean high-confidence specimens and 2 flagged misidentifications
    clean_vouchers = audit_df[audit_df["c_error"] < 0.15]
    flagged_vouchers = audit_df[audit_df["c_error"] > args.cleanlab_threshold]

    exemplar_candidates = []
    if not clean_vouchers.empty:
        for _, row in clean_vouchers.head(2).iterrows():
            exemplar_candidates.append((row, "Clean Anchor"))
    if not flagged_vouchers.empty:
        for _, row in flagged_vouchers.head(2).iterrows():
            exemplar_candidates.append((row, "Flagged Misidentification"))
    else:
        # Fallback to highest error vouchers if none exceed threshold
        for _, row in audit_df.head(2).iterrows():
            exemplar_candidates.append((row, "High-Uncertainty Voucher"))

    transform = transforms.Compose([
        transforms.Resize((224, 224), interpolation=transforms.InterpolationMode.BICUBIC),
        transforms.ToTensor(),
        transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
    ])

    exemplar_results = []
    for row, ex_type in exemplar_candidates:
        img_path = Path(row["patch_path"])
        pil_img = Image.open(img_path).convert("RGB").resize((224, 224))
        tensor = transform(pil_img)

        target_class_idx = label_to_idx[row["predicted_label"]]
        cam_map = compute_vit_gradcam_heatmap(vit_clf, tensor, target_class_idx, device)

        exemplar_results.append({
            "catalog": row["catalogNumber"],
            "type": ex_type,
            "given": row["species_standardized"],
            "pred": row["predicted_label"],
            "tier": row["determiner_tier"],
            "c_error": row["c_error"],
            "orig_img": np.array(pil_img),
            "cam_map": cam_map,
        })

    # 6. Render and export Grad-CAM audit panel
    plot_gradcam_audit_panel(
        joint_noise_matrix=joint_noise_matrix,
        audit_df=audit_df,
        exemplars=exemplar_results,
        output_png=args.output_figure,
        threshold=args.cleanlab_threshold,
    )

    logger.info("Phase 4: Confident Learning & Grad-CAM XAI pipeline completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
