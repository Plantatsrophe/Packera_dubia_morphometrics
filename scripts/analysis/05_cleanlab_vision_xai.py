#!/usr/bin/env python3
"""
===============================================================================
Script: 05_cleanlab_vision_xai.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Phase 5 Entry Point: Deep Vision Feature Extraction (DINOv2), Confident
    Learning Label Noise Curation (Cleanlab), and Grad-CAM Explainable AI (XAI).
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Ensure project root is in sys.path
_script_root = Path(__file__).resolve().parents[2]
if str(_script_root) not in sys.path:
    sys.path.insert(0, str(_script_root))

from scripts.analysis.cleanlab_curator import (
    compute_out_of_fold_probabilities,
    run_confident_learning_audit,
)
from scripts.analysis.dinov2_embeddings import (
    TARGET_TAXA,
    RosettePatchDataset,
    extract_dinov2_embeddings,
    load_and_link_rosette_patches,
    standardize_packera_taxon,
)
from scripts.analysis.gradcam_visualizer import (
    blend_heatmap_on_image,
    generate_gradcam_panel,
)

__all__ = [
    "TARGET_TAXA",
    "standardize_packera_taxon",
    "RosettePatchDataset",
    "load_and_link_rosette_patches",
    "extract_dinov2_embeddings",
    "compute_out_of_fold_probabilities",
    "run_confident_learning_audit",
    "blend_heatmap_on_image",
    "generate_gradcam_panel",
]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger("Cleanlab_Vision_XAI")


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for Phase 5 analysis."""
    parser = argparse.ArgumentParser(
        description="Phase 5: DINOv2 Deep Vision Embedding & Cleanlab XAI Audit"
    )
    parser.add_argument(
        "--rosette-dir",
        type=Path,
        default=Path("data/cropped_patches"),
        help="Directory containing cropped basal rosette patch images",
    )
    parser.add_argument(
        "--vouchers-csv",
        type=Path,
        default=Path("data/tables/curated_vouchers.csv"),
        help="Curated voucher metadata table CSV",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/tables/label_noise_audit.csv"),
        help="Destination table for label quality audit",
    )
    parser.add_argument(
        "--output-figure",
        type=Path,
        default=Path("outputs/figures/GradCAM_audit_panel.png"),
        help="Destination path for Grad-CAM diagnostic panel figure",
    )
    parser.add_argument(
        "--error-threshold",
        type=float,
        default=0.85,
        help="Threshold for confident learning label noise flagging (default: 0.85)",
    )
    return parser.parse_args()


def main() -> None:
    """Orchestrates Phase 5 Deep Vision & Confident Learning workflow."""
    args = parse_args()
    logger.info("=" * 80)
    logger.info("    PHASE 5: DINOV2 DEEP VISION EMBEDDING & CLEANLAB XAI AUDIT    ")
    logger.info("=" * 80)

    # 1. Load and link rosette crops with voucher metadata
    df, class_map = load_and_link_rosette_patches(args.rosette_dir, args.vouchers_csv)
    if df.empty:
        logger.warning("No rosette patch records were matched. Check inputs.")
        return

    records = df.to_dict("records")
    logger.info(f"Loaded {len(records)} specimens for DINOv2 feature extraction.")

    # 2. Extract self-supervised DINOv2 embeddings
    features, labels, cat_nums = extract_dinov2_embeddings(records)

    # 3. Fit out-of-fold cross-validated probabilities
    pred_probs, acc, f1 = compute_out_of_fold_probabilities(features, labels)

    # 4. Execute Confident Learning label noise audit
    audit_df = run_confident_learning_audit(
        pred_probs=pred_probs,
        labels=labels,
        records_df=df,
        class_names=TARGET_TAXA,
        error_threshold=args.error_threshold,
    )

    # 5. Export table
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    audit_df.to_csv(args.output_csv, index=False)
    logger.info(f"Exported label noise audit table -> {args.output_csv}")

    # 6. Generate Grad-CAM diagnostic figures
    flagged = audit_df[audit_df["is_label_corrupted"]].to_dict("records")
    if not flagged:
        flagged = audit_df[audit_df["is_cleanlab_issue"]].to_dict("records")
    generate_gradcam_panel(flagged, args.output_figure)

    logger.info("=" * 80)
    logger.info("Phase 5 Cleanlab Vision XAI workflow completed successfully.")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
