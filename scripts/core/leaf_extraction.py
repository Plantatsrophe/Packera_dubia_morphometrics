
import os
import sys
import logging
import math
import numpy as np
import cv2
import json
import glob
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
from pathlib import Path
from scipy import ndimage
from skimage.morphology import skeletonize

# Common imports
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE, DEFAULT_RAW_DIR
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable

import argparse
from scripts.core.leaf_cv_utils import detect_sheet_artifacts, apply_hard_artifact_blanking, detect_native_dpi_regions, frangi_vesselness_filter_2d, trace_3point_anatomical_spine, extract_edt_point_seeds, segment_leaves_sam2_or_watershed, evaluate_convexity_and_solidity, filter_involucre_profiles
from scripts.core.leaf_morphometrics import detect_leaf_midrib_axis, align_mask_horizontally, reconstruct_bilateral_symmetry, extract_open_margin_curve_and_traits

logger = setup_logging()

DEFAULT_MODEL_PATH = DEFAULT_WORKSPACE / "models" / "yolov8_leaf_best.pt"

OUTPUT_DIRS = {
    "annotations": DEFAULT_WORKSPACE / "data" / "masks" / "annotations",
    "rosettes_dense": DEFAULT_WORKSPACE / "data" / "masks" / "rosettes_dense",
    "basal_leaves_raw": DEFAULT_WORKSPACE / "data" / "masks" / "basal_leaves_raw",
    "tier1_intact": DEFAULT_WORKSPACE / "data" / "masks" / "tier1_intact",
    "tier2_reflected": DEFAULT_WORKSPACE / "data" / "masks" / "tier2_reflected",
    "tier3_open_curves": DEFAULT_WORKSPACE / "data" / "masks" / "tier3_open_curves",
    "capitula": DEFAULT_WORKSPACE / "data" / "masks" / "capitula",
    "qc_overlays": DEFAULT_WORKSPACE / "outputs" / "extraction_qc",
    "tables": DEFAULT_WORKSPACE / "data" / "tables"
}

def process_voucher_precision(
    image_path: Path,
    output_dirs: Dict[str, Path],
    save_overlays: bool = True,
    model: Optional[Any] = None,
    conf_threshold: float = 0.25
) -> List[Dict[str, Any]]:
    """
    Executes the full 5-Stage Precision Extraction Pipeline on a single voucher sheet.
    """
    catalog_number = image_path.stem
    logger.debug(f"Processing voucher: {catalog_number} ({image_path.name})")

    image_bgr = cv2.imread(str(image_path))
    if image_bgr is None:
        logger.error(f"Failed to read voucher image: {image_path}")
        return []

    h, w = image_bgr.shape[:2]
    qc_records: List[Dict[str, Any]] = []

    # Stage 1: Pre-Emptive Sheet Parsing & Hard Artifact Blanking
    artifacts = detect_sheet_artifacts(image_bgr)
    clean_sheet = apply_hard_artifact_blanking(image_bgr, artifacts, fill_color=(255, 255, 255))

    for idx, slip_art in enumerate(artifacts.get("annotation_slip", [])):
        sx1, sy1, sx2, sy2 = slip_art["bbox"]
        slip_crop = image_bgr[sy1:sy2, sx1:sx2]
        if slip_crop.size > 0:
            slip_path = output_dirs["annotations"] / f"{catalog_number}_slip_{idx}.jpg"
            cv2.imwrite(str(slip_path), slip_crop)

    # Stage 2: Native-DPI Rosette & Cyme Sub-Image Cropping
    rosette_bbox, cyme_bboxes = detect_native_dpi_regions(clean_sheet)

    if rosette_bbox is not None:
        rx1, ry1, rx2, ry2 = rosette_bbox
        rosette_crop_bgr = clean_sheet[ry1:ry2, rx1:rx2]
        rosette_crop_gray = cv2.cvtColor(rosette_crop_bgr, cv2.COLOR_BGR2GRAY)
        inv_rosette_gray = 255 - rosette_crop_gray

        if rosette_crop_bgr.size > 0:
            rosette_out_path = output_dirs["rosettes_dense"] / f"{catalog_number}_rosette.jpg"
            cv2.imwrite(str(rosette_out_path), rosette_crop_bgr)

        blurred_r = cv2.GaussianBlur(inv_rosette_gray, (5, 5), 0)
        _, rosette_binary = cv2.threshold(blurred_r, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

        # Stage 4: Distance Transform (EDT) Peak Seeding
        dist_map, point_seeds = extract_edt_point_seeds(rosette_binary, min_distance_px=30)
        leaf_masks = segment_leaves_sam2_or_watershed(rosette_crop_bgr, rosette_binary, point_seeds)

        if not leaf_masks:
            cnts_r, _ = cv2.findContours(rosette_binary, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            for cr in cnts_r:
                if cv2.contourArea(cr) > 400:
                    lm = np.zeros_like(rosette_binary)
                    cv2.drawContours(lm, [cr], -1, 255, -1)
                    leaf_masks.append(lm)

        # Stages 3 & 5: Anatomical Spines & Solidity Gatekeeper
        for leaf_idx, lmask in enumerate(leaf_masks):
            cnts_l, _ = cv2.findContours(lmask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            if not cnts_l:
                continue
            bc = max(cnts_l, key=cv2.contourArea)
            lx, ly, lw, lh = cv2.boundingRect(bc)
            leaf_patch = rosette_crop_bgr[ly:ly+lh, lx:lx+lw]
            if leaf_patch.size == 0:
                continue

            # YOLO Gatekeeper Injection
            if model is not None:
                results = model.predict(leaf_patch, conf=conf_threshold, verbose=False)
                has_leaf = False
                for r in results:
                    if r.boxes is not None and len(r.boxes.cls) > 0:
                        for cls_id in r.boxes.cls:
                            if int(cls_id.item()) == 0:  # 0 is basal_leaf
                                has_leaf = True
                                break
                if not has_leaf:
                    logger.debug(f"YOLO Gatekeeper REJECTED leaf_{leaf_idx}: No basal_leaf detected in patch.")
                    continue

            raw_patch_path = output_dirs["basal_leaves_raw"] / f"{catalog_number}_leaf_{leaf_idx}.jpg"
            cv2.imwrite(str(raw_patch_path), leaf_patch)

            leaf_roi_mask = lmask[ly:ly+lh, lx:lx+lw]
            leaf_roi_gray = rosette_crop_gray[ly:ly+lh, lx:lx+lw]
            ucs_score, solidity, is_clean, leaf_contour, conv_hull = evaluate_convexity_and_solidity(leaf_roi_mask)

            spine_info = trace_3point_anatomical_spine(leaf_roi_mask, leaf_roi_gray)

            midrib_angle, apex_pt, base_pt, _ = detect_leaf_midrib_axis(leaf_roi_mask, leaf_contour)
            aligned_mask, _, _ = align_mask_horizontally(leaf_roi_mask, apex_pt, base_pt)

            if is_clean and ucs_score >= 0.85 and solidity >= 0.72:
                mask_filename = f"{catalog_number}_leaf_{leaf_idx}.png" if len(leaf_masks) > 1 else f"{catalog_number}_leaf.png"
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
                logger.debug(f"[Tier 1] Intact leaf extracted: {mask_out_path.name} (UCS={ucs_score:.2f}, Solidity={solidity:.2f})")

            else:
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
                    logger.debug(f"[Tier 2] Reconstructed leaf: {mask_out_path.name} (UCS={ucs_score:.2f}, half={half_used})")

                else:
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
                    logger.debug(f"[Tier 3] Open curve extracted: {curve_out_path.name}")

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
            ("annotation_slip", (255, 128, 0))
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
            overlay_canvas = cv2.resize(overlay_canvas, (int(w * scale), int(h * scale)))
        cv2.imwrite(str(overlay_out), overlay_canvas)

    return qc_records


def run_pipeline(
    raw_dir: Path = DEFAULT_RAW_DIR,
    model_path: Path = DEFAULT_MODEL_PATH,
    conf_threshold: float = 0.25,
    limit: Optional[int] = None,
    save_overlays: bool = True,
    clean: bool = False
) -> pd.DataFrame:
    """
    Main orchestration routine iterating across all voucher images in
    data/raw_vouchers/ and exporting the master leaf extraction QC log.

    Args:
        raw_dir: Path to directory containing raw herbarium voucher scans.
        model_path: Path to fine-tuned YOLOv8-seg weights.
        conf_threshold: Confidence threshold for organ instance detection.
        limit: Optional cap on the number of vouchers to process (for testing).
        save_overlays: Whether to generate and save visual QC overlay images.
        clean: If True, purges all existing masks, cropped patches, overlays,
               and prior QC logs to guarantee a completely fresh run.

    Returns:
        pd.DataFrame: Master QC log dataframe containing all extraction records.
    """
    import shutil

    logger.info("=" * 80)
    logger.info("STARTING PACKERA 5-STAGE PRECISION BOTANICAL EXTRACTION PIPELINE")
    logger.info("=" * 80)

    # Optional Clean Reset: Purge existing masks, crops, and QC tables if --clean is specified
    if clean:
        logger.info("[CLEAN RESET] Purging prior masks, cropped patches, and QC logs as requested...")
        for name, dir_path in OUTPUT_DIRS.items():
            if dir_path.exists() and dir_path.is_dir():
                # Avoid deleting the parent data/tables/ folder itself; just remove leaf_extraction_qc.csv
                if name == "tables":
                    qc_file = dir_path / "leaf_extraction_qc.csv"
                    if qc_file.exists():
                        qc_file.unlink()
                else:
                    # Remove all files inside the output subdirectories
                    for f in dir_path.glob("*"):
                        if f.is_file() and f.name != "desktop.ini":
                            try:
                                f.unlink()
                            except Exception as e:
                                logger.debug(f"Could not delete {f}: {e}")
        logger.info("[CLEAN RESET] Output directories successfully wiped and reset.")

    # Ensure all output directories exist
    for dir_path in OUTPUT_DIRS.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    image_patterns = ["*.jpg", "*.jpeg", "*.png", "*.tif", "*.JPG", "*.JPEG", "*.PNG"]
    image_paths: List[Path] = []
    for pat in image_patterns:
        image_paths.extend(raw_dir.glob(pat))

    image_paths = sorted(list(set(image_paths)))
    total_images = len(image_paths)

    logger.info(f"Discovered {total_images} raw voucher images in: {raw_dir}")
    if total_images == 0:
        logger.warning(f"No voucher images found in {raw_dir}. Please run 01_voucher_harvester.py first.")
        return pd.DataFrame()

    if limit and limit > 0:
        image_paths = image_paths[:limit]
        logger.info(f"Subsetting to first {limit} vouchers as requested.")

    try:
        from ultralytics import YOLO
        logger.info(f"Loading YOLO gatekeeper from: {model_path}")
        model = YOLO(model_path)
    except Exception as e:
        logger.error(f"Failed to load YOLO model: {e}")
        logger.warning("Falling back to raw heuristics without gatekeeper!")
        model = None

    all_qc_records: List[Dict[str, Any]] = []

    for img_path in tqdm(image_paths, desc="Extracting Basal Leaves", unit="sheet"):
        try:
            records = process_voucher_precision(
                image_path=img_path,
                output_dirs=OUTPUT_DIRS,
                save_overlays=save_overlays,
                model=model,
                conf_threshold=conf_threshold
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

        logger.info("-" * 80)
        logger.info("LEAF EXTRACTION SUMMARY & TIER BREAKDOWN:")
        logger.info(f"  Total Leaf Evaluations          : {total_eval}")
        logger.info(f"  Tier 1 (Direct Intact Leaf)     : {t1} ({t1/max(total_eval, 1)*100:.1f}%)")
        logger.info(f"  Tier 2 (Symmetry Reconstructed) : {t2} ({t2/max(total_eval, 1)*100:.1f}%)")
        logger.info(f"  Tier 3 (Open Margin Curves)     : {t3} ({t3/max(total_eval, 1)*100:.1f}%)")
        logger.info("-" * 80)

    return qc_df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Packera 5-Stage Precision Botanical Extraction & Symmetry Reconstruction Pipeline"
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=DEFAULT_RAW_DIR,
        help="Path to directory containing raw voucher images (default: data/raw_vouchers/)"
    )
    parser.add_argument(
        "--model-path",
        type=Path,
        default=DEFAULT_MODEL_PATH,
        help="Path to fine-tuned YOLOv8-seg weights (default: models/yolov8_leaf_best.pt)"
    )
    parser.add_argument(
        "--conf-threshold",
        type=float,
        default=0.25,
        help="Confidence threshold for organ instance detection (default: 0.25)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of vouchers to process (useful for rapid testing / debugging)"
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Purge all prior masks, cropped patches, and QC logs before starting execution"
    )
    parser.add_argument(
        "--no-overlays",
        action="store_true",
        help="Disable generation of QC visualization overlay panels"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose DEBUG logging output"
    )
    return parser.parse_args()


