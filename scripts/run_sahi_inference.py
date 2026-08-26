#!/usr/bin/env python3
"""
===============================================================================
Script: run_sahi_inference.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
Author: High-Resolution Image Processing Specialist & Senior AI Engineer
Date: August 2026

Description:
    Dedicated standalone SAHI (Slicing Aided Hyper Inference) component for
    full-resolution herbarium sheet scans (24-50 MP).
    
    Splits high-resolution gigapixel-scale images into overlapping native-DPI
    windows, executes YOLOv8 instance segmentation per tile, and merges
    local bounding boxes and polygon masks across sheet coordinates using
    global Non-Maximum Suppression (NMS).

Key Features:
    - High-precision detection of small foliar organs without downsampling artifacts.
    - Checkpoint resume capability: avoids re-running already completed sheets.
    - Periodic flush to JSON summary for fault-tolerant long-running batches.
    - Optional visualization overlays for botanical morphological validation.

Usage:
    python scripts/run_sahi_inference.py --weights models/yolov8_leaf_best.pt --input-dir data/raw_vouchers --output-dir outputs/sahi_detections --conf 0.25 --iou 0.45
===============================================================================
"""

import os
import sys
import time
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Any, Set

# Add repository root and scripts directory to sys.path for absolute imports
root_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root_dir))
sys.path.insert(0, str(root_dir / "scripts"))

# Import SAHI inference engine from core module
from scripts.native_dpi_patch_tiler import (
    HerbariumSAHIInference,
    CLASS_NAMES
)

# Configure structured logging format
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] (%(name)s) %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("SAHIInferenceRunner")


def run_sahi_inference(
    weights: str = "models/yolov8_leaf_best.pt",
    input_dir: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    summary_output: Optional[Path] = None,
    slice_size: int = 1024,
    overlap: float = 0.20,
    conf_threshold: float = 0.25,
    iou_threshold: float = 0.45,
    visualize: bool = False,
    force: bool = False,
    limit: Optional[int] = None,
    flush_interval: int = 25
) -> List[Dict[str, Any]]:
    """
    Execute full-dataset SAHI sliced inference across high-resolution voucher sheets.

    Args:
        weights (str): Path to YOLOv8 segmentation weights (.pt file).
        input_dir (Path, optional): Directory containing raw voucher images.
                                    Defaults to root_dir / "data" / "raw_vouchers".
        output_dir (Path, optional): Directory to store detection outputs/visualizations.
                                     Defaults to root_dir / "outputs" / "sahi_detections".
        summary_output (Path, optional): Path to output JSON summary file.
                                         Defaults to root_dir / "outputs" / "sahi_inference_summary.json".
        slice_size (int): Dimension (height and width) of sliding window slices.
        overlap (float): Slicing overlap ratio between adjacent windows (0.0 to 0.5).
        conf_threshold (float): Model detection confidence threshold.
        iou_threshold (float): NMS IoU threshold for merging overlapping detections.
        visualize (bool): If True, generates visual detection overlays on specimen sheets.
        force (bool): If True, ignores previous summary checkpoints and re-infers all sheets.
        limit (int, optional): Optional limit on number of sheets to process.
        flush_interval (int): Frequency of writing intermediate results to JSON summary.

    Returns:
        List[Dict[str, Any]]: List of per-sheet SAHI detection summaries.
    """
    start_time = time.time()

    # Resolve paths
    resolved_input_dir = Path(input_dir) if input_dir else root_dir / "data" / "raw_vouchers"
    resolved_output_dir = Path(output_dir) if output_dir else root_dir / "outputs" / "sahi_detections"
    resolved_summary_path = Path(summary_output) if summary_output else root_dir / "outputs" / "sahi_inference_summary.json"

    # Resolve weights path with graceful fallback
    weights_path = Path(weights) if Path(weights).is_absolute() else root_dir / weights
    if not weights_path.exists():
        fallback_weights = root_dir / "yolov8x-seg.pt"
        if fallback_weights.exists():
            logger.warning("Specified weights '%s' not found. Falling back to default '%s'", weights_path, fallback_weights)
            weights_path = fallback_weights
        else:
            logger.warning("Model weights '%s' does not exist on disk yet. Proceeding with specified path.", weights_path)

    # Discover voucher image files
    image_extensions = ("*.jpg", "*.jpeg", "*.tif", "*.tiff", "*.png", "*.JPG", "*.JPEG")
    image_files: List[Path] = []
    if resolved_input_dir.exists():
        for ext in image_extensions:
            image_files.extend(resolved_input_dir.glob(ext))
    image_files = sorted(list(set(image_files)))

    if limit is not None and limit > 0:
        image_files = image_files[:limit]

    logger.info("=================================================================")
    logger.info("STARTING SAHI SLICED INFERENCE PIPELINE")
    logger.info("Model Weights:        %s", weights_path)
    logger.info("Input Directory:      %s", resolved_input_dir)
    logger.info("Output Directory:     %s", resolved_output_dir)
    logger.info("Summary Path:         %s", resolved_summary_path)
    logger.info("Total Sheets Found:   %d", len(image_files))
    logger.info("Slice Dimensions:     %dx%d px | Overlap: %.1f%%", slice_size, slice_size, overlap * 100)
    logger.info("Thresholds:           Confidence >= %.2f | NMS IoU <= %.2f", conf_threshold, iou_threshold)
    logger.info("Visualization:        %s", visualize)
    logger.info("Force Re-run:         %s", force)
    logger.info("=================================================================")

    if not image_files:
        logger.warning("No image files found in input directory '%s'.", resolved_input_dir)
        return []

    # Initialize SAHI inference engine
    sahi_engine = HerbariumSAHIInference(
        model_path=str(weights_path),
        slice_height=slice_size,
        slice_width=slice_size,
        overlap_height_ratio=overlap,
        overlap_width_ratio=overlap,
        confidence_threshold=conf_threshold,
        nms_iou_threshold=iou_threshold
    )

    resolved_output_dir.mkdir(parents=True, exist_ok=True)
    all_sahi_results: List[Dict[str, Any]] = []
    processed_images: Set[str] = set()

    # Check for existing summary checkpoint if not forcing complete re-run
    if resolved_summary_path.exists() and not force:
        try:
            with open(resolved_summary_path, "r", encoding="utf-8") as f:
                all_sahi_results = json.load(f)
                for r in all_sahi_results:
                    if "image_path" in r:
                        processed_images.add(Path(r["image_path"]).name)
            logger.info("Resuming SAHI: %d sheets already inferred, %d remaining",
                        len(processed_images), len(image_files) - len(processed_images))
        except Exception as err:
            logger.warning("Could not load previous SAHI summary checkpoint: %s. Starting fresh.", err)
            all_sahi_results = []
            processed_images = set()

    total_sheets = len(image_files)
    inferred_in_session = 0

    for idx, img_p in enumerate(image_files, 1):
        if img_p.name in processed_images:
            continue

        # Optional visual debugging overlay
        vis_path = (resolved_output_dir / f"{img_p.stem}_sahi_vis.jpg") if visualize else None

        try:
            sheet_result = sahi_engine.predict_sheet(img_p, visualize_output_path=vis_path)
            all_sahi_results.append(sheet_result)
            processed_images.add(img_p.name)
            inferred_in_session += 1

            if inferred_in_session % 10 == 0 or idx == total_sheets:
                logger.info("SAHI Inference progress: %d/%d sheets completed (%.1f%%)...",
                            len(processed_images), total_sheets, (len(processed_images) / total_sheets) * 100)

            # Periodically flush summary to disk for crash resilience
            if inferred_in_session % flush_interval == 0:
                resolved_summary_path.parent.mkdir(parents=True, exist_ok=True)
                with open(resolved_summary_path, "w", encoding="utf-8") as f:
                    json.dump(all_sahi_results, f, indent=2)

        except Exception as err:
            logger.error("Error during SAHI inference on sheet '%s': %s", img_p.name, err, exc_info=True)

    # Final summary flush to disk
    resolved_summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(resolved_summary_path, "w", encoding="utf-8") as f:
        json.dump(all_sahi_results, f, indent=2)

    elapsed_time = time.time() - start_time
    logger.info("=================================================================")
    logger.info("SAHI INFERENCE COMPLETE in %.2f minutes (%.2f seconds)", elapsed_time / 60.0, elapsed_time)
    logger.info("Total Sheets in Dataset:     %d", total_sheets)
    logger.info("Total Inferred Sheets:       %d", len(all_sahi_results))
    logger.info("New Sheets Inferred Today:   %d", inferred_in_session)
    logger.info("Summary File Saved:          %s", resolved_summary_path)
    logger.info("=================================================================")

    return all_sahi_results


def main():
    """Command-line interface entrypoint for SAHI sliced inference runner."""
    parser = argparse.ArgumentParser(
        description="Full-Dataset SAHI Sliced Inference for High-Resolution Herbarium Scans."
    )
    parser.add_argument(
        "--weights",
        type=str,
        default="models/yolov8_leaf_best.pt",
        help="Path to trained YOLOv8 segmentation model weights (.pt)."
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="data/raw_vouchers",
        help="Path to directory containing high-resolution raw voucher images."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="outputs/sahi_detections",
        help="Destination directory for detection masks and visualizations."
    )
    parser.add_argument(
        "--summary-output",
        type=str,
        default="outputs/sahi_inference_summary.json",
        help="Path to JSON summary output file."
    )
    parser.add_argument(
        "--slice-size",
        "--tile-size",
        dest="slice_size",
        type=int,
        default=1024,
        help="Sliding slice dimension in pixels (width and height)."
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=0.20,
        help="Sliding window overlap ratio between adjacent slices (e.g. 0.20 for 20%%)."
    )
    parser.add_argument(
        "--conf",
        type=float,
        default=0.25,
        help="Confidence threshold for instance detections."
    )
    parser.add_argument(
        "--iou",
        type=float,
        default=0.45,
        help="Non-Maximum Suppression (NMS) IoU threshold for stitching slice predictions."
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Generate visual overlays with rendered bounding boxes and masks."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-inference on all specimen sheets even if checkpoint summary exists."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of sheets to process (useful for rapid testing/debugging)."
    )
    parser.add_argument(
        "--flush-interval",
        type=int,
        default=25,
        help="Number of sheets between incremental JSON checkpoint flushes to disk."
    )

    args = parser.parse_args()

    run_sahi_inference(
        weights=args.weights,
        input_dir=Path(args.input_dir),
        output_dir=Path(args.output_dir),
        summary_output=Path(args.summary_output),
        slice_size=args.slice_size,
        overlap=args.overlap,
        conf_threshold=args.conf,
        iou_threshold=args.iou,
        visualize=args.visualize,
        force=args.force,
        limit=args.limit,
        flush_interval=args.flush_interval
    )


if __name__ == "__main__":
    main()
