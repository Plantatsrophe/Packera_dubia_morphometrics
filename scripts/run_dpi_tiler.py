#!/usr/bin/env python3
"""
===============================================================================
Script: run_dpi_tiler.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
Author: High-Resolution Image Processing Specialist & Senior AI Engineer
Date: August 2026

Description:
    Dedicated standalone high-throughput native-DPI patch tiling component.
    Processes full-resolution specimen sheet scans (24-50 MP) across multi-core
    CPU worker pools using sliding-window geometry, polygon clipping, and
    background paper filtering.

Key Features:
    - Multi-process parallelization with ProcessPoolExecutor.
    - Preserves native spatial resolution for micro-morphological feature extraction.
    - Performs dynamic boolean polygon re-projection & visible area filtering (< 15%).
    - Hard negative background retention sampling (5% retention).
    - Checkpoints execution metrics to JSON for auditability.

Usage:
    python scripts/run_dpi_tiler.py --input-dir data/raw_vouchers --labels-dir data/yolo_dataset/labels --output-dir data/tiled_dataset --num-workers 16
===============================================================================
"""

import os
import sys
import time
import json
import logging
import argparse
import concurrent.futures
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

# Add repository root and scripts directory to sys.path for absolute imports
root_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root_dir))
sys.path.insert(0, str(root_dir / "scripts"))

# Import native DPI patch tiler class from core module
from scripts.native_dpi_patch_tiler import (
    NativeDPIPatchTiler,
    CLASS_NAMES
)

# Configure structured logging format
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] (%(name)s) %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("DPITilerRunner")


# -----------------------------------------------------------------------------
# Standalone Worker Function for Multi-Process Parallel Patch Tiling
# -----------------------------------------------------------------------------
def _process_single_sheet_worker(task_args: Tuple) -> Dict[str, Any]:
    """
    Independent worker function invoked by ProcessPoolExecutor across multiple CPU cores.
    Processes a single specimen sheet:
      1. Decodes native image raster.
      2. Computes sliding-window grid slices.
      3. Reprojects YOLO polygon contours to local tile space.
      4. Discards truncated annotations (< min_area_ratio).
      5. Sub-samples empty paper background tiles (bg_keep_prob).
      6. Writes tile images and corresponding YOLO label files.

    Args:
        task_args (Tuple): Tuple containing:
            - img_p (Path): Absolute path to the specimen image.
            - lbl_p (Path or None): Path to matching ground-truth label file.
            - tile_size (int): Dimension of square tiles in pixels (e.g., 1024).
            - overlap (float): Overlap ratio between adjacent windows (e.g., 0.20).
            - min_area_ratio (float): Minimum visible area ratio to retain clipped polygon.
            - bg_keep_prob (float): Probability of keeping an empty background tile.
            - output_dir (Path): Destination root directory for tiled dataset.
            - visualize (bool): Whether to render visual debug overlay images.
            - skip_existing (bool): Whether to skip sheets if already tiled.

    Returns:
        dict: Worker telemetry metrics (tile counts, class counts, sheet status).
    """
    (
        img_p,
        lbl_p,
        tile_size,
        overlap,
        min_area_ratio,
        bg_keep_prob,
        output_dir,
        visualize,
        skip_existing
    ) = task_args

    # Instantiate lightweight local worker tiler instance
    worker_tiler = NativeDPIPatchTiler(
        tile_size=tile_size,
        overlap=overlap,
        min_area_ratio=min_area_ratio,
        bg_keep_prob=bg_keep_prob,
        output_dir=output_dir,
        visualize=visualize,
        skip_existing=skip_existing
    )

    # Execute processing for this individual sheet
    worker_tiler.process_sheet(image_path=img_p, label_path=lbl_p)
    return worker_tiler.metrics


def run_dpi_tiling(
    input_dir: Optional[Path] = None,
    labels_dir: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    summary_output: Optional[Path] = None,
    tile_size: int = 1024,
    overlap: float = 0.20,
    min_area_ratio: float = 0.15,
    bg_keep_prob: float = 0.05,
    num_workers: int = 32,
    visualize: bool = False,
    force: bool = False,
    limit: Optional[int] = None
) -> Dict[str, Any]:
    """
    Execute high-throughput multi-process native-DPI sliding-window patch tiling.

    Args:
        input_dir (Path, optional): Directory containing raw voucher images.
                                    Defaults to root_dir / "data" / "raw_vouchers".
        labels_dir (Path, optional): Directory containing ground-truth YOLO label files.
                                     Defaults to root_dir / "data" / "yolo_dataset" / "labels".
        output_dir (Path, optional): Destination folder for tiled dataset.
                                     Defaults to root_dir / "data" / "tiled_dataset".
        summary_output (Path, optional): Path to save summary metrics JSON.
                                         Defaults to root_dir / "outputs" / "tiling_summary.json".
        tile_size (int): Width and height of square tiles in pixels.
        overlap (float): Overlap ratio between adjacent windows (0.0 to 0.5).
        min_area_ratio (float): Minimum preserved instance area ratio (0.0 to 1.0).
        bg_keep_prob (float): Probability of keeping pure background tiles.
        num_workers (int): Number of parallel CPU worker processes.
        visualize (bool): If True, outputs diagnostic bounding box overlays.
        force (bool): If True, re-processes sheets even if already tiled.
        limit (int, optional): Optional limit on number of sheets to process (for debugging).

    Returns:
        dict: Aggregated summary metrics of all tiled sheets.
    """
    start_time = time.time()

    # Set default directory paths if not explicitly provided
    resolved_input_dir = Path(input_dir) if input_dir else root_dir / "data" / "raw_vouchers"
    resolved_labels_dir = Path(labels_dir) if labels_dir else root_dir / "data" / "yolo_dataset" / "labels"
    resolved_output_dir = Path(output_dir) if output_dir else root_dir / "data" / "tiled_dataset"
    resolved_summary_path = Path(summary_output) if summary_output else root_dir / "outputs" / "tiling_summary.json"

    # Auto-adjust worker process count to machine CPU capacity
    available_cpus = os.cpu_count() or 4
    effective_workers = max(1, min(available_cpus, num_workers))

    # Scan for voucher images across standard herbarium raster formats
    image_extensions = ("*.jpg", "*.jpeg", "*.tif", "*.tiff", "*.png", "*.JPG", "*.JPEG")
    image_files: List[Path] = []
    if resolved_input_dir.exists():
        for ext in image_extensions:
            image_files.extend(resolved_input_dir.glob(ext))
    image_files = sorted(list(set(image_files)))

    if limit is not None and limit > 0:
        image_files = image_files[:limit]

    logger.info("=================================================================")
    logger.info("STARTING NATIVE-DPI PATCH TILING PIPELINE")
    logger.info("Input Directory:      %s", resolved_input_dir)
    logger.info("Labels Directory:     %s", resolved_labels_dir)
    logger.info("Output Directory:     %s", resolved_output_dir)
    logger.info("Total Sheets Found:   %d", len(image_files))
    logger.info("CPU Worker Pool:      %d processes (Available: %d cores)", effective_workers, available_cpus)
    logger.info("Tile Parameters:      %dx%d px | Overlap: %.1f%% | Min Area: %.1f%% | BG Keep: %.1f%%",
                tile_size, tile_size, overlap * 100, min_area_ratio * 100, bg_keep_prob * 100)
    logger.info("Skip Existing Tiles:  %s (Force Mode: %s)", not force, force)
    logger.info("=================================================================")

    if not image_files:
        logger.warning("No image files found in input directory '%s'.", resolved_input_dir)
        return {
            "total_sheets_processed": 0,
            "total_tiles_generated": 0,
            "positive_tiles_generated": 0,
            "negative_tiles_retained": 0,
            "negative_tiles_discarded": 0,
            "class_tile_counts": {}
        }

    # Pre-index all ground-truth label files for fast O(1) lookups by stem
    label_index: Dict[str, Path] = {}
    if resolved_labels_dir.exists():
        for lp in resolved_labels_dir.rglob("*.txt"):
            label_index[lp.stem] = lp
        logger.info("Indexed %d ground-truth label files from '%s'", len(label_index), resolved_labels_dir)

    resolved_output_dir.mkdir(parents=True, exist_ok=True)
    master_metrics: Dict[str, Any] = {
        "total_sheets_processed": 0,
        "total_tiles_generated": 0,
        "positive_tiles_generated": 0,
        "negative_tiles_retained": 0,
        "negative_tiles_discarded": 0,
        "class_tile_counts": {name: 0 for name in CLASS_NAMES}
    }

    # Prepare worker tasks
    task_list = []
    for img_p in image_files:
        lbl_p = label_index.get(img_p.stem)
        task_list.append((
            img_p,
            lbl_p,
            tile_size,
            overlap,
            min_area_ratio,
            bg_keep_prob,
            resolved_output_dir,
            visualize,
            not force  # skip_existing
        ))

    # Execute concurrent patch tiling across CPU cores
    with concurrent.futures.ProcessPoolExecutor(max_workers=effective_workers) as executor:
        futures = {executor.submit(_process_single_sheet_worker, t): t[0] for t in task_list}
        completed_count = 0
        total_tasks = len(task_list)

        for future in concurrent.futures.as_completed(futures):
            completed_count += 1
            if completed_count % 50 == 0 or completed_count == total_tasks:
                pct = (completed_count / total_tasks) * 100
                logger.info("Tiling progress: %d/%d sheets processed (%.1f%%)...", completed_count, total_tasks, pct)

            try:
                sheet_metrics = future.result()
                master_metrics["total_sheets_processed"] += 1
                master_metrics["total_tiles_generated"] += sheet_metrics.get("total_tiles_generated", 0)
                master_metrics["positive_tiles_generated"] += sheet_metrics.get("positive_tiles_generated", 0)
                master_metrics["negative_tiles_retained"] += sheet_metrics.get("negative_tiles_retained", 0)
                master_metrics["negative_tiles_discarded"] += sheet_metrics.get("negative_tiles_discarded", 0)
                
                # Aggregate class counts
                for k, v in sheet_metrics.get("class_tile_counts", {}).items():
                    master_metrics["class_tile_counts"][k] = master_metrics["class_tile_counts"].get(k, 0) + v
            except Exception as err:
                sheet_name = futures[future].name
                logger.error("Error processing sheet '%s': %s", sheet_name, err, exc_info=True)

    # Save aggregated tiling metrics
    resolved_summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(resolved_summary_path, "w", encoding="utf-8") as f:
        json.dump(master_metrics, f, indent=2)

    elapsed_time = time.time() - start_time
    logger.info("=================================================================")
    logger.info("NATIVE-DPI TILING COMPLETE in %.2f minutes (%.2f seconds)", elapsed_time / 60.0, elapsed_time)
    logger.info("Total Sheets Processed:    %d", master_metrics["total_sheets_processed"])
    logger.info("Total Tiles Generated:    %d", master_metrics["total_tiles_generated"])
    logger.info("Positive Tiles Generated: %d", master_metrics["positive_tiles_generated"])
    logger.info("Negative Tiles Retained:  %d", master_metrics["negative_tiles_retained"])
    logger.info("Negative Tiles Discarded: %d", master_metrics["negative_tiles_discarded"])
    logger.info("Summary File Saved:       %s", resolved_summary_path)
    logger.info("=================================================================")

    return master_metrics


def main():
    """Command-line interface entrypoint for native-DPI sliding window patch tiler."""
    parser = argparse.ArgumentParser(
        description="High-Throughput Native-DPI Sliding-Window Patch Tiler for Herbarium Scans."
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default="data/raw_vouchers",
        help="Path to directory containing high-resolution raw voucher images."
    )
    parser.add_argument(
        "--labels-dir",
        type=str,
        default="data/yolo_dataset/labels",
        help="Path to directory containing YOLO segmentation labels (.txt)."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/tiled_dataset",
        help="Destination directory where generated patch tiles and labels are saved."
    )
    parser.add_argument(
        "--summary-output",
        type=str,
        default="outputs/tiling_summary.json",
        help="Path to JSON summary output file."
    )
    parser.add_argument(
        "--tile-size",
        type=int,
        default=1024,
        help="Tile dimension in pixels (width and height for square patch)."
    )
    parser.add_argument(
        "--overlap",
        type=float,
        default=0.20,
        help="Sliding window overlap ratio between adjacent tiles (e.g. 0.20 for 20%%)."
    )
    parser.add_argument(
        "--min-area-ratio",
        type=float,
        default=0.15,
        help="Minimum visible area ratio required to keep clipped polygon annotations."
    )
    parser.add_argument(
        "--bg-keep-prob",
        type=float,
        default=0.05,
        help="Probability of retaining empty hard-negative background paper tiles."
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=32,
        help="Number of concurrent CPU worker processes for parallel tiling."
    )
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="Generate QC visual inspection overlays with rendered polygon annotations."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-tiling of all sheets even if patches already exist in output directory."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of sheets to process (useful for rapid testing/debugging)."
    )

    args = parser.parse_args()

    run_dpi_tiling(
        input_dir=Path(args.input_dir),
        labels_dir=Path(args.labels_dir),
        output_dir=Path(args.output_dir),
        summary_output=Path(args.summary_output),
        tile_size=args.tile_size,
        overlap=args.overlap,
        min_area_ratio=args.min_area_ratio,
        bg_keep_prob=args.bg_keep_prob,
        num_workers=args.num_workers,
        visualize=args.visualize,
        force=args.force,
        limit=args.limit
    )


if __name__ == "__main__":
    main()
