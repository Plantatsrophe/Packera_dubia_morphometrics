#!/usr/bin/env python3
"""
===============================================================================
Script: 02_segment_and_extract.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Production Runner for Step 02:
    Direct voucher image segmentation and leaf extraction pipeline for
    downstream Elliptic Fourier Analysis (EFA).
    
    Key Capabilities:
      1. Direct Ingestion: Reads curated_vouchers.csv directly without intermediate
         staging or symlink copying.
      2. PointRend Segmentation: Direct inference using the fine-tuned checkpoint
         (models/lm2_packera_pcd_finetuned.pth).
      3. Decoupled Asynchronous Scale Detection: Hough-transform ruler detection
         runs in parallel; scale failure never aborts silhouette extraction.
      4. In-Sheet Centroid Clustering: Groups organs into plant individuals
         without separate DBSCAN micro-scripts.
      5. 2-Path Botanical Leaf Extraction:
           - Tier 1 (Direct Pristine): UCS >= 0.85 & Solidity >= 0.72.
           - Tier 2 (Hemi-Blade Bilateral Reflection): Midrib axis alignment,
             half-blade defect gating, and symmetrical reflection.
           - Failed QC Routing: Heavily clumped rosettes or unsegmentable leaves
             routed to data/tables/failed_qc_vouchers.csv.
      6. Standardized Contour Export: Saves continuous 2D (x, y) coordinates
         to data/contours/{catalogNumber}_leaf{id}.csv and binary silhouettes
         to data/masks/.

Usage:
    python scripts/pipeline/02_segment_and_extract.py \\
        --vouchers data/tables/curated_vouchers.csv \\
        --model-weights models/lm2_packera_pcd_finetuned.pth \\
        --output-dir data/ \\
        --device cuda
===============================================================================
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig
from scripts.vision.lm2_geometry_utils import (
    cluster_plant_individuals,
    compute_geometric_metrics,
    detect_ruler_scale_hough,
    export_standardized_contour,
    extract_tier1_pristine,
    extract_tier2_reflected,
)


# LeafMachine2 and Detectron2 paths
LM2_DIR = PROJECT_ROOT / "LeafMachine2" / "leafmachine2"
if LM2_DIR.exists() and str(LM2_DIR) not in sys.path:
    sys.path.insert(0, str(LM2_DIR))

DETECTRON2_PATH = LM2_DIR / "segmentation" / "detectron2"
if DETECTRON2_PATH.exists() and str(DETECTRON2_PATH) not in sys.path:
    sys.path.insert(0, str(DETECTRON2_PATH))

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("SegmentAndExtract")


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class DetectedInstance:
    """Represents a segmented botanical instance on a herbarium sheet."""
    catalog_number: str
    leaf_id: int
    bbox: Tuple[int, int, int, int]  # (ymin, xmin, ymax, xmax)
    mask: np.ndarray                 # Binary uint8 mask (H, W)
    score: float
    class_id: int
    plant_individual_id: int = 0
    assigned_tier: str = "unassigned"
    ucs_score: float = 0.0
    solidity: float = 0.0
    midrib_angle_deg: float = 0.0
    pixels_per_mm: Optional[float] = None
    mask_path: Optional[str] = None
    contour_path: Optional[str] = None
    rejection_reason: Optional[str] = None


# Note: Geometric quality gatekeeping, 2-path reflection/extraction, standardized
# contour export, DBSCAN spatial clustering, and Hough ruler detection are
# consolidated in scripts.vision.lm2_geometry_utils.


# =============================================================================
# Direct PointRend Model Inference Engine
# =============================================================================

class PointRendInferenceEngine:
    """Encapsulates Detectron2 PointRend model loading and inference."""

    def __init__(
        self,
        weights_path: Union[str, Path],
        device: str = "cuda",
        score_thresh: float = 0.50
    ) -> None:
        self.weights_path = Path(weights_path)
        self.device = device if (device == "cuda" and torch_cuda_available()) else "cpu"
        self.score_thresh = score_thresh
        self.predictor = None
        self.model = None
        self.cfg = None
        self.aug = None
        self._load_model()

    def _load_model(self) -> None:
        """Loads Detectron2 PointRend model from weights checkpoint."""
        try:
            from detectron2.config import get_cfg
            from detectron2.projects import point_rend
            from detectron2.checkpoint import DetectionCheckpointer
            from detectron2.modeling import build_model
            import detectron2.data.transforms as T
        except ImportError as e:
            logger.warning(f"Detectron2 not available: {e}. Running in simulation/fallback mode.")
            return

        cfg = get_cfg()
        point_rend.add_pointrend_config(cfg)

        # Locate accompanying configuration YAML
        cfg_candidates = [
            self.weights_path.parent / "cfg_output.yaml",
            PROJECT_ROOT / "LeafMachine2" / "leafmachine2" / "segmentation" / "models" / "Packera_LeafPriority" / "cfg_output.yaml",
            PROJECT_ROOT / "models" / "cfg_output.yaml",
        ]
        cfg_yaml = None
        for cand in cfg_candidates:
            if cand.exists():
                cfg_yaml = cand
                break

        if cfg_yaml:
            cfg.merge_from_file(str(cfg_yaml))
        else:
            logger.warning("No cfg_output.yaml found; using default PointRend config.")

        cfg.MODEL.WEIGHTS = str(self.weights_path.resolve())
        cfg.MODEL.ROI_HEADS.SCORE_THRESH_TEST = self.score_thresh
        cfg.MODEL.DEVICE = self.device

        model = build_model(cfg)
        model.eval()
        checkpointer = DetectionCheckpointer(model)
        checkpointer.load(cfg.MODEL.WEIGHTS)

        self.cfg = cfg
        self.model = model
        self.aug = T.ResizeShortestEdge([cfg.INPUT.MIN_SIZE_TEST, cfg.INPUT.MIN_SIZE_TEST], cfg.INPUT.MAX_SIZE_TEST)
        logger.info(f"Loaded PointRend model from {self.weights_path} on {self.device}")

    def predict(self, image_bgr: np.ndarray) -> List[Tuple[Tuple[int, int, int, int], np.ndarray, float, int]]:
        """
        Runs inference on BGR image.
        Returns list of (bbox, mask, score, class_id).
        """
        if self.model is None:
            return []

        import torch
        height, width = image_bgr.shape[:2]
        with torch.no_grad():
            img_transformed = self.aug.get_transform(image_bgr).apply_image(image_bgr)
            tensor_img = torch.as_tensor(img_transformed.astype("float32").transpose(2, 0, 1))
            inputs = [{"image": tensor_img, "height": height, "width": width}]
            preds = self.model(inputs)[0]

        instances = preds["instances"].to("cpu")
        results = []

        if not instances.has("pred_masks"):
            return results

        scores = instances.scores.numpy()
        classes = instances.pred_classes.numpy()
        boxes = instances.pred_boxes.tensor.numpy()
        masks = instances.pred_masks.numpy()

        for i in range(len(scores)):
            if scores[i] < self.score_thresh:
                continue
            box = boxes[i]
            ymin, xmin, ymax, xmax = int(box[1]), int(box[0]), int(box[3]), int(box[2])
            mask_uint8 = (masks[i] > 0).astype(np.uint8) * 255
            results.append(((ymin, xmin, ymax, xmax), mask_uint8, float(scores[i]), int(classes[i])))

        return results


def torch_cuda_available() -> bool:
    """Helper to safely check CUDA availability."""
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False


# =============================================================================
# Core Pipeline Execution Runner
# =============================================================================

class SegmentAndExtractPipeline:
    """Unified runner for leaf segmentation, 2-path extraction, and manifest logging."""

    def __init__(
        self,
        vouchers_csv: Union[str, Path],
        model_weights: Union[str, Path],
        output_dir: Union[str, Path],
        device: str = "cuda",
        min_solidity: float = 0.72,
        min_ucs: float = 0.85,
        score_thresh: float = 0.40,
        num_ruler_workers: int = 4
    ) -> None:
        self.vouchers_csv = Path(vouchers_csv)
        self.model_weights = Path(model_weights)
        self.output_dir = Path(output_dir)
        self.device = device
        self.min_solidity = min_solidity
        self.min_ucs = min_ucs
        self.score_thresh = score_thresh
        self.num_ruler_workers = num_ruler_workers

        self.masks_dir = self.output_dir / "masks"
        self.contours_dir = self.output_dir / "contours"
        self.tables_dir = self.output_dir / "tables"

        for d in [self.masks_dir, self.contours_dir, self.tables_dir]:
            d.mkdir(parents=True, exist_ok=True)

        self.ruler_executor = concurrent.futures.ThreadPoolExecutor(max_workers=self.num_ruler_workers)
        self.model_engine = PointRendInferenceEngine(
            weights_path=self.model_weights,
            device=self.device,
            score_thresh=self.score_thresh
        )

    def shutdown(self) -> None:
        """Cleans up background thread pools."""
        self.ruler_executor.shutdown(wait=False)

    def process_voucher(
        self,
        catalog_number: str,
        image_path: Path
    ) -> Tuple[List[Dict[str, Any]], Optional[Dict[str, Any]]]:
        """Processes a single herbarium voucher specimen sheet."""
        if not image_path.exists():
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "image_not_found",
                "details": f"File does not exist at {image_path}"
            }
            return [], failure_record

        img_bgr = cv2.imread(str(image_path))
        if img_bgr is None:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "corrupted_image",
                "details": "Failed to decode image with OpenCV"
            }
            return [], failure_record

        sheet_h, sheet_w = img_bgr.shape[:2]

        # 1. Asynchronously launch decoupled Hough ruler detection
        ruler_future = self.ruler_executor.submit(detect_ruler_scale_hough, img_bgr)

        # 2. Run PointRend instance segmentation
        raw_detections = self.model_engine.predict(img_bgr)

        # Collect asynchronous scale result (with graceful fallback)
        try:
            pixels_per_mm = ruler_future.result(timeout=5.0)
        except Exception as e:
            logger.debug(f"Ruler detection timeout or exception for {catalog_number}: {e}")
            pixels_per_mm = None

        if not raw_detections:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "no_leaves_detected",
                "details": "Model produced 0 candidate leaf segmentations"
            }
            return [], failure_record

        # Build detected instances
        instances: List[DetectedInstance] = []
        for idx, (bbox, mask, score, class_id) in enumerate(raw_detections, 1):
            instances.append(DetectedInstance(
                catalog_number=catalog_number,
                leaf_id=idx,
                bbox=bbox,
                mask=mask,
                score=score,
                class_id=class_id,
                pixels_per_mm=pixels_per_mm
            ))

        # 3. Vectorized Centroid Spatial Clustering
        instances = cluster_plant_individuals(instances, sheet_w, sheet_h)

        extracted_records: List[Dict[str, Any]] = []
        rejected_instances: List[DetectedInstance] = []

        # 4. 2-Path Botanical Leaf Extraction Logic
        for inst in instances:
            ucs, solidity, angle_deg, p_apex, p_base = compute_geometric_metrics(inst.mask)
            inst.ucs_score = ucs
            inst.solidity = solidity
            inst.midrib_angle_deg = angle_deg

            # Path A: Tier 1 Direct Pristine
            if ucs >= self.min_ucs and solidity >= self.min_solidity:
                inst.assigned_tier = "tier1"
                mask_path = extract_tier1_pristine(inst.mask, catalog_number, inst.leaf_id, self.output_dir)
                contour_path = export_standardized_contour(inst.mask, catalog_number, inst.leaf_id, self.output_dir)
                inst.mask_path = mask_path
                inst.contour_path = contour_path
            else:
                # Path B: Tier 2 Hemi-Blade Bilateral Reflection
                mask_path, reflected_mask = extract_tier2_reflected(
                    inst.mask, catalog_number, inst.leaf_id, p_apex, p_base, self.output_dir
                )
                if mask_path and reflected_mask is not None:
                    inst.assigned_tier = "tier2"
                    contour_path = export_standardized_contour(
                        reflected_mask, catalog_number, inst.leaf_id, self.output_dir
                    )
                    inst.mask_path = mask_path
                    inst.contour_path = contour_path
                else:
                    inst.assigned_tier = "rejected"
                    inst.rejection_reason = "clumped_rosette_or_severe_occlusion"
                    rejected_instances.append(inst)
                    continue

            # Record successfully extracted leaf
            extracted_records.append({
                "catalogNumber": catalog_number,
                "plant_individual_id": inst.plant_individual_id,
                "leaf_id": inst.leaf_id,
                "assigned_tier": inst.assigned_tier,
                "ucs_score": round(inst.ucs_score, 4),
                "solidity": round(inst.solidity, 4),
                "midrib_angle_deg": round(inst.midrib_angle_deg, 2),
                "pixels_per_mm": round(inst.pixels_per_mm, 4) if inst.pixels_per_mm else np.nan,
                "mask_path": inst.mask_path,
                "contour_path": inst.contour_path
            })

        # Route vouchers with 0 valid silhouettes to failed manifest
        failure_record = None
        if not extracted_records:
            failure_record = {
                "catalogNumber": catalog_number,
                "image_path": str(image_path),
                "failure_reason": "clumped_rosette_occlusion",
                "details": f"{len(rejected_instances)} leaves detected but all rejected by geometric gatekeeping"
            }

        return extracted_records, failure_record

    def run(self, limit: Optional[int] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """Executes full extraction pipeline across curated vouchers."""
        logger.info("==================================================================")
        logger.info("Starting Step 02: Direct Segmentation & Botanical Leaf Extraction")
        logger.info(f"Vouchers Metadata: {self.vouchers_csv}")
        logger.info(f"PointRend Model: {self.model_weights}")
        logger.info(f"Output Directory: {self.output_dir}")
        logger.info(f"Device: {self.device}")
        logger.info(f"Gatekeeper Thresholds: Min Solidity={self.min_solidity}, Min UCS={self.min_ucs}")
        logger.info("==================================================================")

        if not self.vouchers_csv.exists():
            raise FileNotFoundError(f"Curated vouchers table not found at {self.vouchers_csv}")

        df_vouchers = pd.read_csv(self.vouchers_csv)
        logger.info(f"Loaded {len(df_vouchers)} voucher records.")

        if "catalogNumber" not in df_vouchers.columns or "image_path" not in df_vouchers.columns:
            raise KeyError("Vouchers CSV must contain 'catalogNumber' and 'image_path' columns.")

        if limit and limit > 0:
            df_vouchers = df_vouchers.head(limit)
            logger.info(f"Applied execution limit: processing first {limit} vouchers.")

        all_extracted_records: List[Dict[str, Any]] = []
        all_failed_records: List[Dict[str, Any]] = []

        total_vouchers = len(df_vouchers)
        for i, (_, row) in enumerate(df_vouchers.iterrows(), 1):
            cat_num = str(row["catalogNumber"]).strip()
            raw_img_path = Path(str(row["image_path"]).strip())
            if not raw_img_path.is_absolute():
                raw_img_path = PROJECT_ROOT / raw_img_path

            extracted, failure = self.process_voucher(cat_num, raw_img_path)
            all_extracted_records.extend(extracted)
            if failure:
                all_failed_records.append(failure)

            if i % 25 == 0 or i == total_vouchers:
                logger.info(f"Progress: {i}/{total_vouchers} vouchers processed "
                            f"({len(all_extracted_records)} leaves extracted, {len(all_failed_records)} failed QC).")

        # Save manifests
        extracted_df = pd.DataFrame(all_extracted_records)
        if extracted_df.empty:
            extracted_df = pd.DataFrame(columns=[
                "catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier",
                "ucs_score", "solidity", "midrib_angle_deg", "pixels_per_mm",
                "mask_path", "contour_path"
            ])

        failed_df = pd.DataFrame(all_failed_records)
        if failed_df.empty:
            failed_df = pd.DataFrame(columns=[
                "catalogNumber", "image_path", "failure_reason", "details"
            ])

        manifest_path = self.tables_dir / "extracted_leaves_manifest.csv"
        failed_qc_path = self.tables_dir / "failed_qc_vouchers.csv"

        extracted_df.to_csv(manifest_path, index=False)
        failed_df.to_csv(failed_qc_path, index=False)

        logger.info("==================================================================")
        logger.info(f"Extraction completed!")
        logger.info(f"  Extracted Leaves: {len(extracted_df)} -> {manifest_path}")
        logger.info(f"  Failed QC Vouchers: {len(failed_df)} -> {failed_qc_path}")
        if not extracted_df.empty and "assigned_tier" in extracted_df.columns:
            tier_dist = extracted_df["assigned_tier"].value_counts().to_dict()
            for tier, count in tier_dist.items():
                pct = (count / len(extracted_df)) * 100
                logger.info(f"    {tier.upper()}: {count} ({pct:.1f}%)")
        logger.info("==================================================================")

        self.shutdown()
        return extracted_df, failed_df


# =============================================================================
# CLI Interface
# =============================================================================

def parse_args() -> argparse.Namespace:
    cfg = PipelineConfig.from_yaml()
    parser = argparse.ArgumentParser(
        description="Direct PointRend segmentation, 2-path leaf extraction, and contour coordinate export."
    )
    parser.add_argument(
        "--vouchers",
        type=Path,
        default=cfg.paths.curated_vouchers_csv,
        help="Path to curated vouchers CSV containing catalogNumber and image_path."
    )
    parser.add_argument(
        "--model-weights",
        type=Path,
        default=cfg.segmentation.model_weights,
        help="Path to fine-tuned PointRend model checkpoint."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=cfg.paths.workspace_root / "data",
        help="Root output directory for masks, contours, and manifest tables."
    )
    parser.add_argument(
        "--device",
        type=str,
        default=cfg.segmentation.device if torch_cuda_available() else "cpu",
        choices=["cuda", "cpu"],
        help="Compute device for inference (default: cuda if available, else cpu)."
    )
    parser.add_argument(
        "--min-solidity",
        type=float,
        default=cfg.thresholds.min_solidity,
        help=f"Minimum solidity for Tier 1 direct pristine leaf extraction (default: {cfg.thresholds.min_solidity})."
    )
    parser.add_argument(
        "--min-ucs",
        type=float,
        default=cfg.thresholds.min_ucs,
        help=f"Minimum Unoccluded Completeness Score for Tier 1 (default: {cfg.thresholds.min_ucs})."
    )
    parser.add_argument(
        "--score-thresh",
        type=float,
        default=cfg.segmentation.score_thresh,
        help=f"Detection confidence score threshold (default: {cfg.segmentation.score_thresh})."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of vouchers to process."
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pipeline = SegmentAndExtractPipeline(
        vouchers_csv=args.vouchers,
        model_weights=args.model_weights,
        output_dir=args.output_dir,
        device=args.device,
        min_solidity=args.min_solidity,
        min_ucs=args.min_ucs,
        score_thresh=args.score_thresh
    )
    pipeline.run(limit=args.limit)


if __name__ == "__main__":
    main()
