"""
scripts/core/leaf_extraction.py
===============================
Precision Botanical Leaf Extractor & Rosette Disentanglement Pipeline.
Integrates Euclidean Distance Transform (EDT) peak seeding with SAM 2 point
prompting to disentangle overlapping basal rosettes without full-sheet tiling fragmentation.

Adheres to the 5-Stage Precision Botanical Extraction Workflow:
  Stage 1: Pre-emptive metadata artifact sterilization via ArtifactFilterGatekeeper.
  Stage 2: Native-DPI sub-image cropping for basal rosette and synflorescence regions.
  Stage 3: YOLO organ instance detection & 3-point anatomical spine tracing.
  Stage 4: EDT peak seeding + SAM 2 point prompting (or watershed fallback) for rosette disentanglement.
  Stage 5: Solidity filtering (>= 0.72) & linear organ routing via botanical_topology_classifier.

4-Tier Extraction Routing:
  Tier 1: Intact leaves (Solidity >= 0.72, UCS >= 0.85) -> data/masks/tier1_intact/
  Tier 2: Hemi-blade bilateral symmetry reflection -> data/masks/tier2_reflected/
  Tier 3: Open-outline continuous margin curves -> data/masks/tier3_open_curves/
  Tier 4: Unsegmented dense rosette crops -> data/cropped_patches/rosettes_dense/
"""

from __future__ import annotations

import glob
import logging
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import pandas as pd

from scripts.core.botanical_topology_classifier import (
    CAULINE_STEM,
    LEAF_PETIOLE,
    ROOT_RHIZOME,
    classify_elongated_botanical_organ,
    compute_topological_summary,
)
from scripts.core.config import (
    CLASS_COLORS_BGR,
    CLASS_MAP,
    CLASS_NAMES,
    DEFAULT_RAW_DIR,
    DEFAULT_WORKSPACE,
)
from scripts.core.data_structures import (
    ArtifactDetection,
    FilterResult,
    GeometricMetrics,
    InstanceAnnotation,
    SpectralMetrics,
    TextureMetrics,
)
from scripts.core.gatekeeper_engine import ArtifactFilterGatekeeper
from scripts.core.leaf_cv_utils import (
    apply_hard_artifact_blanking,
    detect_native_dpi_regions,
    detect_sheet_artifacts,
    evaluate_convexity_and_solidity,
    extract_edt_point_seeds,
    filter_involucre_profiles,
    frangi_vesselness_filter_2d,
    load_sam2_predictor,
    segment_leaves_sam2_or_watershed,
    trace_3point_anatomical_spine,
)
from scripts.core.leaf_morphometrics import (
    align_mask_horizontally,
    detect_leaf_midrib_axis,
    extract_open_margin_curve_and_traits,
    reconstruct_bilateral_symmetry,
)
from scripts.core.logger import setup_logging

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

logger = setup_logging()

DEFAULT_MODEL_PATH = DEFAULT_WORKSPACE / "models" / "yolov8_leaf_best.pt"
DEFAULT_SAM2_CHECKPOINT = DEFAULT_WORKSPACE / "models" / "checkpoints" / "sam2_hiera_large.pt"
DEFAULT_SAM2_CONFIG = "sam2_hiera_l.yaml"

OUTPUT_DIRS = {
    "annotations": DEFAULT_WORKSPACE / "data" / "masks" / "annotations",
    "rosettes_dense": DEFAULT_WORKSPACE / "data" / "cropped_patches" / "rosettes_dense",
    "basal_leaves_raw": DEFAULT_WORKSPACE / "data" / "masks" / "basal_leaves_raw",
    "tier1_intact": DEFAULT_WORKSPACE / "data" / "masks" / "tier1_intact",
    "tier2_reflected": DEFAULT_WORKSPACE / "data" / "masks" / "tier2_reflected",
    "tier3_open_curves": DEFAULT_WORKSPACE / "data" / "masks" / "tier3_open_curves",
    "capitula": DEFAULT_WORKSPACE / "data" / "masks" / "capitula",
    "qc_overlays": DEFAULT_WORKSPACE / "outputs" / "extraction_qc",
    "tables": DEFAULT_WORKSPACE / "data" / "tables"
}

__all__ = [
    "process_voucher_precision",
    "run_pipeline",
    "DEFAULT_MODEL_PATH",
    "DEFAULT_SAM2_CHECKPOINT",
    "DEFAULT_SAM2_CONFIG",
    "OUTPUT_DIRS",
]


def process_voucher_precision(
    image_path: Path,
    output_dirs: Dict[str, Path],
    save_overlays: bool = True,
    model: Optional[Any] = None,
    conf_threshold: float = 0.25,
    use_sam2: bool = True,
    sam2_predictor: Optional[Any] = None,
    sam2_checkpoint: Optional[Union[str, Path]] = None,
    sam2_model_cfg: Optional[Union[str, Path]] = None,
    gatekeeper: Optional[ArtifactFilterGatekeeper] = None
) -> List[Dict[str, Any]]:
    """
    Executes the 5-Stage Precision Botanical Extraction Workflow for a single herbarium sheet scan.
    """
    catalog_number = image_path.stem
    logger.info(f"Processing voucher: {catalog_number} ({image_path.name})")

    image_bgr = cv2.imread(str(image_path))
    if image_bgr is None:
        logger.error(f"Failed to read image at {image_path}")
        return []

    h, w = image_bgr.shape[:2]
    qc_records: List[Dict[str, Any]] = []

    # Initialize gatekeeper if not passed
    if gatekeeper is None:
        gatekeeper = ArtifactFilterGatekeeper()

    # Stage 1: Pre-emptive metadata artifact sterilization
    artifacts = detect_sheet_artifacts(image_bgr)
    clean_sheet = apply_hard_artifact_blanking(image_bgr, artifacts)

    # Stage 2: Native-DPI sub-image cropping for basal rosette cluster
    rosette_bbox, cyme_bboxes = detect_native_dpi_regions(clean_sheet)

    if rosette_bbox is not None:
        rx1, ry1, rx2, ry2 = rosette_bbox
        rosette_crop_bgr = clean_sheet[ry1:ry2, rx1:rx2]
        rosette_gray = cv2.cvtColor(rosette_crop_bgr, cv2.COLOR_BGR2GRAY)
        inv_rosette_gray = 255 - rosette_gray

        _, rosette_bin = cv2.threshold(inv_rosette_gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
        rosette_bin = cv2.morphologyEx(rosette_bin, cv2.MORPH_OPEN, kernel)

        rosette_out_path = output_dirs["rosettes_dense"] / f"{catalog_number}_rosette_clump.jpg"
        cv2.imwrite(str(rosette_out_path), rosette_crop_bgr)

        # Stage 3: YOLO organ instance detection & 3-point anatomical spine tracing
        leaf_masks: List[np.ndarray] = []
        if model is not None:
            try:
                results = model.predict(rosette_crop_bgr, conf=conf_threshold, verbose=False)
                if results and len(results) > 0 and results[0].masks is not None:
                    yolo_masks = results[0].masks.data.cpu().numpy()
                    classes = results[0].boxes.cls.cpu().numpy()

                    for m_idx, cls_id in enumerate(classes):
                        cname = CLASS_NAMES[int(cls_id)] if int(cls_id) < len(CLASS_NAMES) else str(cls_id)
                        if "leaf" in cname or "blade" in cname:
                            m_raw = yolo_masks[m_idx]
                            m_resized = cv2.resize(m_raw, (rx2 - rx1, ry2 - ry1), interpolation=cv2.INTER_NEAREST)
                            m_u8 = (m_resized > 0.5).astype(np.uint8) * 255
                            if cv2.countNonZero(m_u8) >= 300:
                                leaf_masks.append(m_u8)
            except Exception as e:
                logger.warning(f"YOLO inference on rosette crop failed: {e}")

        # Stage 4: EDT peak seeding + SAM 2 point prompting for disentanglement
        if not leaf_masks:
            _, point_seeds = extract_edt_point_seeds(rosette_bin, min_distance_px=30, relative_threshold=0.25)
            logger.debug(f"Extracted {len(point_seeds)} EDT point seeds from rosette clump.")

            if point_seeds:
                leaf_masks = segment_leaves_sam2_or_watershed(
                    rosette_crop_bgr=rosette_crop_bgr,
                    rosette_binary_mask=rosette_bin,
                    point_seeds=point_seeds,
                    use_sam2=use_sam2,
                    sam2_predictor=sam2_predictor,
                    sam2_checkpoint=sam2_checkpoint,
                    sam2_model_cfg=sam2_model_cfg
                )

        extracted_leaf_count = 0
        for leaf_idx, l_mask in enumerate(leaf_masks):
            raw_leaf_filename = f"{catalog_number}_raw_leaf_{leaf_idx}.png"
            cv2.imwrite(str(output_dirs["basal_leaves_raw"] / raw_leaf_filename), l_mask)

            spine_info = trace_3point_anatomical_spine(l_mask, rosette_gray)
            ucs_score, solidity, is_clean, best_cnt, hull = evaluate_convexity_and_solidity(l_mask)

            # Stage 5: Gatekeeper verification and linear organ routing
            filter_res = gatekeeper.validate_candidate_leaf(
                candidate_patch=rosette_crop_bgr,
                candidate_mask=l_mask,
                catalog_number=catalog_number,
                patch_id=str(leaf_idx),
                candidate_class="basal_leaf_blade",
                is_rgb=False
            )

            if not filter_res.is_valid:
                logger.debug(f"Candidate {leaf_idx} rejected by gatekeeper: {filter_res.primary_rejection_reason}")
                continue

            topo_res = classify_elongated_botanical_organ(
                organ_mask=l_mask,
                bounding_box=[rx1, ry1, rx2, ry2],
                sheet_height=h,
                sheet_width=w
            )

            if topo_res.predicted_class in (CAULINE_STEM, ROOT_RHIZOME):
                logger.debug(f"Candidate {leaf_idx} routed to {topo_res.predicted_class} by topology classifier.")
                continue

            midrib_angle = detect_leaf_midrib_axis(l_mask)
            aligned_mask, _ = align_mask_horizontally(l_mask, midrib_angle)

            # 4-Tier Extraction Routing
            if is_clean:
                # Tier 1: Intact leaves
                mask_filename = f"{catalog_number}_intact_{leaf_idx}.png" if len(leaf_masks) > 1 else f"{catalog_number}_intact.png"
                mask_out_path = output_dirs["tier1_intact"] / mask_filename
                cv2.imwrite(str(mask_out_path), aligned_mask)

                qc_records.append({
                    "catalogNumber": catalog_number,
                    "assigned_tier": 1,
                    "ucs_score": round(ucs_score, 4),
                    "solidity": round(solidity, 4),
                    "midrib_angle": round(midrib_angle, 2),
                    "symmetry_reconstructed": "FALSE",
                    "mask_path": str(mask_out_path.as_posix()),
                    "leaf_index": leaf_idx,
                    "petiole_length_px": round(spine_info["petiole_length_px"], 2),
                    "lamina_length_px": round(spine_info["lamina_length_px"], 2),
                    "status": "TIER_1_INTACT_EXTRACTED"
                })
                extracted_leaf_count += 1
                logger.debug(f"[Tier 1] Intact leaf extracted: {mask_out_path.name} (UCS={ucs_score:.2f}, Solidity={solidity:.2f})")

            else:
                # Tier 2: Bilateral Symmetry Reflection
                reflected_mask, is_reflected, half_used = reconstruct_bilateral_symmetry(aligned_mask)

                if is_reflected and reflected_mask is not None:
                    mask_filename = f"{catalog_number}_reflected_{leaf_idx}.png" if len(leaf_masks) > 1 else f"{catalog_number}_reflected.png"
                    mask_out_path = output_dirs["tier2_reflected"] / mask_filename
                    cv2.imwrite(str(mask_out_path), reflected_mask)

                    qc_records.append({
                        "catalogNumber": catalog_number,
                        "assigned_tier": 2,
                        "ucs_score": round(ucs_score, 4),
                        "solidity": round(solidity, 4),
                        "midrib_angle": round(midrib_angle, 2),
                        "symmetry_reconstructed": "TRUE",
                        "mask_path": str(mask_out_path.as_posix()),
                        "leaf_index": leaf_idx,
                        "petiole_length_px": round(spine_info["petiole_length_px"], 2),
                        "lamina_length_px": round(spine_info["lamina_length_px"], 2),
                        "status": f"TIER_2_SYMMETRY_REFLECTED_{half_used.upper()}"
                    })
                    extracted_leaf_count += 1
                    logger.debug(f"[Tier 2] Reconstructed leaf: {mask_out_path.name} (UCS={ucs_score:.2f}, half={half_used})")

                else:
                    # Tier 3: Open-Outline Continuous Margin Curves
                    curve_filename = f"{catalog_number}_curve_{leaf_idx}.csv" if len(leaf_masks) > 1 else f"{catalog_number}_curve.csv"
                    curve_out_path = output_dirs["tier3_open_curves"] / curve_filename

                    curve_df, traits = extract_open_margin_curve_and_traits(aligned_mask, catalog_number)
                    curve_df.to_csv(curve_out_path, index=False)

                    qc_records.append({
                        "catalogNumber": catalog_number,
                        "assigned_tier": 3,
                        "ucs_score": round(ucs_score, 4),
                        "solidity": round(solidity, 4),
                        "midrib_angle": round(midrib_angle, 2),
                        "symmetry_reconstructed": "FALSE",
                        "mask_path": str(curve_out_path.as_posix()),
                        "leaf_index": leaf_idx,
                        "petiole_length_px": round(traits["petiole_length_px"], 2),
                        "lamina_length_px": round(spine_info["lamina_length_px"], 2),
                        "status": "TIER_3_OPEN_CURVE_EXTRACTED"
                    })
                    extracted_leaf_count += 1
                    logger.debug(f"[Tier 3] Open curve extracted: {curve_out_path.name}")

        # Tier 4: Dense Clump Fallback
        if extracted_leaf_count == 0 and rosette_crop_bgr.size > 0:
            qc_records.append({
                "catalogNumber": catalog_number,
                "assigned_tier": 4,
                "ucs_score": 0.0,
                "solidity": 0.0,
                "midrib_angle": 0.0,
                "symmetry_reconstructed": "FALSE",
                "mask_path": str(rosette_out_path.as_posix()),
                "leaf_index": -1,
                "petiole_length_px": 0.0,
                "lamina_length_px": 0.0,
                "status": "TIER_4_UNSEGMENTED_ROSETTE_CLUMP"
            })
            logger.debug(f"[Tier 4] Unsegmented dense rosette routed: {rosette_out_path.name}")

    # Stage 2 (cont.): Capitulescence / Involucre extraction
    clean_involucres = filter_involucre_profiles(clean_sheet, cyme_bboxes)
    for inv in clean_involucres:
        inv_idx = inv["involucre_index"]
        head_path = output_dirs["capitula"] / f"{catalog_number}_involucre_{inv_idx}.jpg"
        cv2.imwrite(str(head_path), inv["crop"])

    if save_overlays and qc_records:
        overlay_canvas = image_bgr.copy()
        for art_cls, color in [
            ("herbarium_label", (0, 0, 255)),
            ("color_chart", (255, 0, 255)),
            ("annotation_slip", (255, 128, 0)),
            ("ruler_scale", (0, 255, 255)),
            ("mounting_tape", (128, 128, 128))
        ]:
            for art in artifacts.get(art_cls, []):
                x1, y1, x2, y2 = art["bbox"]
                cv2.rectangle(overlay_canvas, (x1, y1), (x2, y2), color, 3)

        if rosette_bbox is not None:
            rx1, ry1, rx2, ry2 = rosette_bbox
            cv2.rectangle(overlay_canvas, (rx1, ry1), (rx2, ry2), (0, 255, 0), 4)
            cv2.putText(overlay_canvas, "basal_rosette_native_dpi", (rx1, max(ry1 - 10, 25)),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 255, 0), 2, cv2.LINE_AA)

        overlay_out = output_dirs["qc_overlays"] / f"{catalog_number}_qc_overlay.jpg"
        max_dim = 1600
        if max(h, w) > max_dim:
            scale = max_dim / max(h, w)
            overlay_small = cv2.resize(overlay_canvas, (int(w * scale), int(h * scale)))
            cv2.imwrite(str(overlay_out), overlay_small)
        else:
            cv2.imwrite(str(overlay_out), overlay_canvas)

    return qc_records


def run_pipeline(
    raw_dir: Union[str, Path] = DEFAULT_RAW_DIR,
    model_path: Optional[Union[str, Path]] = DEFAULT_MODEL_PATH,
    conf_threshold: float = 0.25,
    use_sam2: bool = True,
    sam2_checkpoint: Optional[Union[str, Path]] = DEFAULT_SAM2_CHECKPOINT,
    sam2_model_cfg: Optional[Union[str, Path]] = DEFAULT_SAM2_CONFIG,
    limit: Optional[int] = None,
    save_overlays: bool = True,
    clean: bool = False
) -> pd.DataFrame:
    """
    Batch executes the 5-Stage Precision Botanical Extraction Workflow over raw herbarium sheets.
    """
    raw_dir = Path(raw_dir)
    for d in OUTPUT_DIRS.values():
        if clean and d.exists() and d != OUTPUT_DIRS["tables"]:
            shutil.rmtree(d)
        d.mkdir(parents=True, exist_ok=True)

    image_paths: List[Path] = []
    for ext in ("*.jpg", "*.jpeg", "*.tif", "*.tiff", "*.png", "*.JPG", "*.JPEG"):
        image_paths.extend(raw_dir.glob(ext))
    image_paths = sorted(list(set(image_paths)))

    if not image_paths:
        logger.warning(f"No specimen images found in {raw_dir}")
        return pd.DataFrame()

    if limit is not None and limit > 0:
        image_paths = image_paths[:limit]

    logger.info(f"Initiating extraction across {len(image_paths)} voucher specimens...")

    model = None
    if model_path is not None:
        try:
            from ultralytics import YOLO
            m_path = Path(model_path)
            if m_path.is_file():
                model = YOLO(str(m_path))
                logger.info(f"Loaded YOLOv8-seg model from {m_path}")
        except Exception as e:
            logger.warning(f"Could not load YOLO model ({e}), continuing without YOLO gating.")

    sam2_predictor = None
    if use_sam2:
        sam2_predictor = load_sam2_predictor(
            checkpoint_path=sam2_checkpoint,
            model_cfg=sam2_model_cfg
        )
        if sam2_predictor is None:
            logger.warning("SAM 2 initialization unsuccessful; pipeline will use marker-controlled watershed.")

    gatekeeper = ArtifactFilterGatekeeper()
    all_qc_records: List[Dict[str, Any]] = []

    for img_path in tqdm(image_paths, desc="Extracting Basal Leaves", unit="sheet"):
        try:
            records = process_voucher_precision(
                image_path=img_path,
                output_dirs=OUTPUT_DIRS,
                save_overlays=save_overlays,
                model=model,
                conf_threshold=conf_threshold,
                use_sam2=use_sam2,
                sam2_predictor=sam2_predictor,
                sam2_checkpoint=sam2_checkpoint,
                sam2_model_cfg=sam2_model_cfg,
                gatekeeper=gatekeeper
            )
            all_qc_records.extend(records)
        except Exception as exc:
            logger.error(f"Error processing {img_path.name}: {exc}", exc_info=True)

    qc_df = pd.DataFrame(all_qc_records)

    req_cols = [
        "catalogNumber",
        "assigned_tier",
        "ucs_score",
        "solidity",
        "midrib_angle",
        "symmetry_reconstructed",
        "mask_path",
    ]
    other_cols = [col for col in qc_df.columns if col not in req_cols]
    final_cols = [col for col in req_cols if col in qc_df.columns] + other_cols
    qc_df = qc_df[final_cols] if not qc_df.empty else pd.DataFrame(columns=req_cols)

    qc_output_path = OUTPUT_DIRS["tables"] / "leaf_extraction_qc.csv"
    qc_df.to_csv(qc_output_path, index=False)
    logger.info(f"Exported master leaf extraction QC log to: {qc_output_path}")

    if not qc_df.empty and "assigned_tier" in qc_df.columns:
        tier_counts = qc_df["assigned_tier"].value_counts().to_dict()
        total_eval = len(qc_df)
        t1 = tier_counts.get(1, 0)
        t2 = tier_counts.get(2, 0)
        t3 = tier_counts.get(3, 0)
        t4 = tier_counts.get(4, 0)

        logger.info("-" * 80)
        logger.info("LEAF EXTRACTION SUMMARY & TIER BREAKDOWN:")
        logger.info(f"  Total Leaf Evaluations          : {total_eval}")
        logger.info(f"  Tier 1 (Direct Intact Leaf)     : {t1} ({t1/max(total_eval, 1)*100:.1f}%)")
        logger.info(f"  Tier 2 (Symmetry Reconstructed) : {t2} ({t2/max(total_eval, 1)*100:.1f}%)")
        logger.info(f"  Tier 3 (Open Margin Curves)     : {t3} ({t3/max(total_eval, 1)*100:.1f}%)")
        logger.info(f"  Tier 4 (Unsegmented Rosettes)   : {t4} ({t4/max(total_eval, 1)*100:.1f}%)")
        logger.info("-" * 80)

    return qc_df
