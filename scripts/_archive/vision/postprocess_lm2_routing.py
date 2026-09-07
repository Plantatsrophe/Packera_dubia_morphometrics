#!/usr/bin/env python3
"""
===============================================================================
Script: postprocess_lm2_routing.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Main orchestrator for LeafMachine2 (LM2) output post-processing, spatial
    DBSCAN multi-plant clustering, geometric quality gatekeeping, 4-tier
    morphometric routing, and master quality control logging.

Usage:
    python scripts/vision/postprocess_lm2_routing.py \\
        --lm2-dir LM2_Project/Data/output/Packera_dubia_LM2/ \\
        --vouchers data/tables/curated_vouchers.csv \\
        --raw-images data/raw_vouchers/ \\
        --output-dir data/ \\
        --min-solidity 0.72 \\
        --min-ucs 0.85
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# Ensure repository root is in sys.path
_repo_root = str(Path(__file__).resolve().parent.parent.parent)
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

import cv2
import numpy as np
import pandas as pd

try:
    from scripts.vision.geometric_gatekeeper import compute_geometric_metrics_and_pose
    from scripts.vision.lm2_data_loader import (
        LeafCandidate,
        discover_lm2_candidates,
        load_and_preprocess_mask,
        load_ruler_calibrations,
    )
    from scripts.vision.morphometric_router import (
        crop_dense_rosette_tier4,
        route_tier1_silhouette,
        route_tier2_reflected,
        route_tier3_open_curve,
    )
    from scripts.vision.spatial_clustering import cluster_voucher_plants_dbscan
except ImportError:
    from geometric_gatekeeper import compute_geometric_metrics_and_pose
    from lm2_data_loader import (
        LeafCandidate,
        discover_lm2_candidates,
        load_and_preprocess_mask,
        load_ruler_calibrations,
    )
    from morphometric_router import (
        crop_dense_rosette_tier4,
        route_tier1_silhouette,
        route_tier2_reflected,
        route_tier3_open_curve,
    )
    from spatial_clustering import cluster_voucher_plants_dbscan

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("LM2_PostProcessing")


def process_voucher_routing(
    catalog_number: str,
    candidates: List[LeafCandidate],
    raw_images_dir: Path,
    output_dir: Path,
    min_solidity: float = 0.72,
    min_ucs: float = 0.85
) -> List[Dict[str, Any]]:
    """Executes clustering, quality gatekeeping, and 4-tier routing per voucher sheet."""
    if not candidates:
        return []

    valid_cands = [
        c for c in candidates
        if (c.crop_image_path and c.crop_image_path.exists()) or
           (c.mask_image_path and c.mask_image_path.exists())
    ]
    if not valid_cands:
        return []

    raw_img_path = raw_images_dir / f"{catalog_number}.jpg"
    sheet_w, sheet_h = 3000, 4000
    if raw_img_path.exists():
        img_temp = cv2.imread(str(raw_img_path))
        if img_temp is not None:
            sheet_h, sheet_w = img_temp.shape[:2]

    # 1. Multi-Plant DBSCAN Spatial Clustering
    candidates = cluster_voucher_plants_dbscan(valid_cands, sheet_width=sheet_w, sheet_height=sheet_h)

    # Group by plant for Tier 4 Dense Rosette Cropping
    candidates_by_plant = defaultdict(list)
    for c in candidates:
        candidates_by_plant[c.plant_individual_id].append(c)

    for plant_id, plant_cands in candidates_by_plant.items():
        crop_dense_rosette_tier4(raw_img_path, plant_cands, plant_id, output_dir)

    qc_records: List[Dict[str, Any]] = []

    # 2. Route Each Basal Leaf Candidate
    for cand in candidates:
        mask = load_and_preprocess_mask(cand)
        if mask is None or np.count_nonzero(mask) < 30:
            continue

        ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics_and_pose(mask)
        cand.ucs_score = float(ucs)
        cand.solidity = float(solidity)
        cand.midrib_angle_deg = float(angle_deg)

        area_px = float(np.count_nonzero(mask))
        cand.area_mm2 = area_px * (cand.scale_mm_per_px ** 2)
        h, w = mask.shape[:2]
        cand.length_mm = float(h) * cand.scale_mm_per_px
        cand.width_mm = float(w) * cand.scale_mm_per_px

        # Routing decision logic
        if ucs >= min_ucs and solidity >= min_solidity:
            cand.assigned_tier = "tier1"
            cand.saved_mask_path = route_tier1_silhouette(mask, cand, output_dir)
        else:
            saved_path, success = route_tier2_reflected(mask, cand, output_dir, p_apex, p_base)
            if success and saved_path:
                cand.assigned_tier = "tier2"
                cand.saved_mask_path = saved_path
            else:
                cand.assigned_tier = "tier3"
                saved_path, petiole_len, blade_w, apex_ang = route_tier3_open_curve(
                    mask, cand, output_dir, cand.scale_mm_per_px
                )
                cand.saved_mask_path = saved_path
                if blade_w > 0:
                    cand.width_mm = blade_w

        record = {
            "catalogNumber": cand.catalog_number,
            "plant_individual_id": cand.plant_individual_id,
            "leaf_id": cand.leaf_id,
            "assigned_tier": cand.assigned_tier,
            "ucs_score": round(cand.ucs_score, 4),
            "solidity": round(cand.solidity, 4),
            "midrib_angle_deg": round(cand.midrib_angle_deg, 2),
            "length_mm": round(cand.length_mm, 2),
            "width_mm": round(cand.width_mm, 2),
            "area_mm2": round(cand.area_mm2, 2),
            "mask_path": cand.saved_mask_path
        }
        qc_records.append(record)

    return qc_records


def run_postprocess_routing(
    lm2_dir: Union[str, Path],
    vouchers_path: Union[str, Path],
    raw_images_dir: Union[str, Path],
    output_dir: Union[str, Path],
    min_solidity: float = 0.72,
    min_ucs: float = 0.85,
    limit: Optional[int] = None
) -> pd.DataFrame:
    """Full pipeline runner for LeafMachine2 post-processing and routing."""
    lm2_dir = Path(lm2_dir)
    vouchers_path = Path(vouchers_path)
    raw_images_dir = Path(raw_images_dir)
    output_dir = Path(output_dir)

    logger.info("==================================================================")
    logger.info("Starting LeafMachine2 Post-Processing & 4-Tier Morphometric Routing")
    logger.info(f"LM2 Output Directory: {lm2_dir}")
    logger.info(f"Curated Vouchers Table: {vouchers_path}")
    logger.info(f"Raw Images Directory: {raw_images_dir}")
    logger.info(f"Output Directory: {output_dir}")
    logger.info(f"Thresholds: Min Solidity = {min_solidity}, Min UCS = {min_ucs}")
    logger.info("==================================================================")

    if not vouchers_path.exists():
        raise FileNotFoundError(f"Vouchers table not found at {vouchers_path}")
    curated_vouchers = pd.read_csv(vouchers_path)
    logger.info(f"Loaded {len(curated_vouchers)} curated voucher records.")

    ruler_calibs = load_ruler_calibrations(lm2_dir)

    candidates_by_voucher = discover_lm2_candidates(
        lm2_dir=lm2_dir,
        curated_vouchers=curated_vouchers,
        ruler_calibs=ruler_calibs,
        raw_images_dir=raw_images_dir
    )

    voucher_list = list(candidates_by_voucher.keys())
    if limit and limit > 0:
        voucher_list = voucher_list[:limit]

    logger.info(f"Processing {len(voucher_list)} voucher specimens...")

    all_qc_records: List[Dict[str, Any]] = []
    for i, cat in enumerate(voucher_list, 1):
        cands = candidates_by_voucher[cat]
        records = process_voucher_routing(
            catalog_number=cat,
            candidates=cands,
            raw_images_dir=raw_images_dir,
            output_dir=output_dir,
            min_solidity=min_solidity,
            min_ucs=min_ucs
        )
        all_qc_records.extend(records)
        if i % 25 == 0 or i == len(voucher_list):
            logger.info(f"Processed {i}/{len(voucher_list)} vouchers ({len(all_qc_records)} leaves routed).")

    qc_df = pd.DataFrame(all_qc_records)
    if qc_df.empty:
        qc_df = pd.DataFrame(columns=[
            "catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier",
            "ucs_score", "solidity", "midrib_angle_deg", "length_mm", "width_mm",
            "area_mm2", "mask_path"
        ])

    qc_output_path = output_dir / "tables" / "leaf_extraction_qc.csv"
    qc_output_path.parent.mkdir(parents=True, exist_ok=True)
    qc_df.to_csv(qc_output_path, index=False)
    logger.info(f"Master QC table exported to {qc_output_path} (Total records: {len(qc_df)}).")

    if not qc_df.empty and "assigned_tier" in qc_df.columns:
        tier_counts = qc_df["assigned_tier"].value_counts().to_dict()
        logger.info("=== 4-Tier Hierarchical Extraction Summary ===")
        for tier, count in tier_counts.items():
            logger.info(f"  {tier.upper()}: {count} leaves ({count / len(qc_df) * 100:.1f}%)")

    return qc_df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="LeafMachine2 post-processing, spatial clustering, and 4-tier morphometric routing."
    )
    parser.add_argument(
        "--lm2-dir",
        type=str,
        default="LM2_Project/Data/output/Packera_dubia_LM2/",
        help="Path to raw LeafMachine2 output directory or parent runs directory."
    )
    parser.add_argument(
        "--vouchers",
        type=str,
        default="data/tables/curated_vouchers.csv",
        help="Path to Darwin Core curated vouchers CSV metadata."
    )
    parser.add_argument(
        "--raw-images",
        type=str,
        default="data/raw_vouchers/",
        help="Path to directory containing full raw herbarium voucher JPEG images."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/",
        help="Root output directory for exported masks, curves, rosettes, and QC table."
    )
    parser.add_argument(
        "--min-solidity",
        type=float,
        default=0.72,
        help="Minimum solidity threshold for Tier 1 intact leaf silhouettes (default: 0.72)."
    )
    parser.add_argument(
        "--min-ucs",
        type=float,
        default=0.85,
        help="Minimum Unoccluded Completeness Score (UCS) for Tier 1 (default: 0.85)."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of voucher specimens to process."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_postprocess_routing(
        lm2_dir=args.lm2_dir,
        vouchers_path=args.vouchers,
        raw_images_dir=args.raw_images,
        output_dir=args.output_dir,
        min_solidity=args.min_solidity,
        min_ucs=args.min_ucs,
        limit=args.limit
    )
