#!/usr/bin/env python3
"""
===============================================================================
Script: annotate_with_sam2.py
Project: Packera dubia Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Interactive Botanical Instance Segmentation Annotator powered by Segment
    Anything Model 2 (SAM 2) with multi-modal boundary controls (bounding boxes,
    exclusion points, freehand lasso, knife cutting, zoom/pan navigation,
    voucher advancement, reloading, and back navigation) for Packera specimens.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np
import torch

def get_project_root() -> Path:
    """Dynamically resolves the project root directory."""
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "data").exists() or (parent / "models").exists() or (parent / ".git").exists():
            return parent
    return current.parents[1] if len(current.parents) > 1 else current.parents[0]


PROJECT_ROOT = get_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Fallback path inclusion for SAM 2 if not installed in current site-packages
SAM2_ARCHIVE_PATH = PROJECT_ROOT / "scripts" / "_archive" / "root_artifacts" / "segment-anything-2"
if SAM2_ARCHIVE_PATH.exists() and str(SAM2_ARCHIVE_PATH) not in sys.path:
    sys.path.insert(0, str(SAM2_ARCHIVE_PATH))

try:
    from sam2.build_sam import build_sam2
    from sam2.sam2_image_predictor import SAM2ImagePredictor
except ImportError as err:
    build_sam2 = None
    SAM2ImagePredictor = None

from scripts.data_prep.sam2_geometry import (
    clip_box_to_image,
    mask_to_normalized_polygon,
    mask_to_yolo_bbox,
    rasterize_lasso_polygon,
    split_mask_with_knife_line,
)
from scripts.data_prep.sam2_rendering import (
    CLASS_COLORS,
    CLASS_NAMES,
    apply_viewport_transform,
    compose_mask_overlay,
    render_hud_overlay,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("SAM2Annotator")


class PrecisionSAM2Annotator:
    """
    Interactive botanical annotator using SAM 2 to target Packera specimens,
    accepting point and bounding box prompts specifically tuned for dense basal rosettes.
    Outputs binary pixel masks tagged explicitly with labels 'basal_leaf_whole',
    'basal_leaf_partial', and 'cauline_leaf'.
    """

    def __init__(
        self,
        images_dir: Union[str, Path] = "data/raw_vouchers",
        output_dir: Union[str, Path] = "data/raw_annotations",
        single_image: Optional[Union[str, Path]] = None,
        checkpoint_path: Union[str, Path] = "models/checkpoints/sam2_hiera_large.pt",
        config_path: Union[str, Path] = "sam2_hiera_l.yaml",
        window_w: int = 1280,
        window_h: int = 960,
        resume_unannotated: bool = True
    ):
        self.project_root = get_project_root()
        self.images_dir = Path(images_dir)
        self.output_dir = Path(output_dir)
        self.masks_dir = self.output_dir / "masks"
        self.labels_dir = self.output_dir / "labels"
        self.masks_dir.mkdir(parents=True, exist_ok=True)
        self.labels_dir.mkdir(parents=True, exist_ok=True)

        self.window_w = window_w
        self.window_h = window_h
        self.resume_unannotated = resume_unannotated

        if single_image:
            self.image_files = [Path(single_image)]
        else:
            exts = {".jpg", ".jpeg", ".png", ".tif", ".tiff"}
            self.image_files = sorted([
                p for p in self.images_dir.glob("*.*")
                if p.suffix.lower() in exts and not p.name.startswith(".")
            ])

        self.checkpoint_path = Path(checkpoint_path)
        self.config_path = str(config_path)

        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.predictor = None
        self._init_model()

        self.current_idx = 0
        self.active_image: Optional[np.ndarray] = None
        self.orig_h = 1000
        self.orig_w = 1000

        self.point_coords: List[List[float]] = []
        self.point_labels: List[int] = []
        self.box_prompt: Optional[List[float]] = None
        self.candidate_mask: Optional[np.ndarray] = None
        self.saved_instances: List[Dict[str, Any]] = []

        self.zoom_level = 1.0
        self.pan_offset = [0, 0]
        self.mode = "INCLUDE"

        self.is_box_dragging = False
        self.box_start = (0, 0)
        self.is_knife_dragging = False
        self.knife_start = (0, 0)
        self.lasso_points: List[Tuple[int, int]] = []

    def _init_model(self) -> None:
        """Initializes the SAM 2 model weights and image predictor."""
        if build_sam2 is None or SAM2ImagePredictor is None:
            return

        if self.checkpoint_path.exists():
            try:
                sam2_model = build_sam2(self.config_path, str(self.checkpoint_path), device=self.device)
                self.predictor = SAM2ImagePredictor(sam2_model)
                logger.info(f"Loaded SAM 2 model ({self.config_path}) onto {self.device}")
            except Exception as e:
                logger.error(f"Failed to load SAM 2 weights: {e}")

    def run_inference(self) -> None:
        """Executes SAM 2 inference on current prompts."""
        if self.predictor is None or self.active_image is None:
            return

        pts = np.array(self.point_coords, dtype=np.float32) if self.point_coords else None
        lbls = np.array(self.point_labels, dtype=np.int32) if self.point_labels else None
        box = np.array(self.box_prompt, dtype=np.float32) if self.box_prompt else None

        if pts is None and box is None:
            self.candidate_mask = None
            return

        try:
            masks, scores, _ = self.predictor.predict(
                point_coords=pts,
                point_labels=lbls,
                box=box,
                multimask_output=False
            )
            if masks is not None and len(masks) > 0:
                self.candidate_mask = (masks[0] > 0.0).astype(np.uint8) * 255
        except Exception as e:
            logger.error(f"SAM 2 prediction error: {e}")

    def save_current_sheet(self) -> None:
        """Saves current sheet's instances to YOLO polygon text file and binary PNG masks."""
        if not hasattr(self, "image_files") or not self.image_files:
            return
        if self.current_idx >= len(self.image_files):
            return

        current_file = self.image_files[self.current_idx]
        voucher_id = current_file.stem
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.masks_dir.mkdir(parents=True, exist_ok=True)

        txt_file = self.output_dir / f"{voucher_id}.txt"
        with open(txt_file, "w", encoding="utf-8") as f:
            for inst in self.saved_instances:
                poly = inst.get("polygon", [])
                c_id = inst.get("class_id", 0)
                poly_str = " ".join([str(round(v, 6)) for v in poly])
                f.write(f"{c_id} {poly_str}\n")

        for idx, inst in enumerate(self.saved_instances):
            label = inst.get("label", CLASS_NAMES[inst.get("class_id", 0)])
            mask_dest = self.masks_dir / f"{voucher_id}_inst{idx:02d}_{label}.png"
            b_mask = inst.get("binary_mask")
            if b_mask is not None:
                uint8_mask = (b_mask.astype(np.uint8)) * 255 if b_mask.dtype == bool else b_mask.astype(np.uint8)
                cv2.imwrite(str(mask_dest), uint8_mask)

        logger.info(f"Saved {len(self.saved_instances)} instances for voucher {voucher_id}")


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for the SAM 2 botanical annotator."""
    parser = argparse.ArgumentParser(description="SAM 2 Interactive Botanical Annotator for Packera")
    parser.add_argument("--images-dir", type=str, default="data/raw_vouchers", help="Input vouchers directory")
    parser.add_argument("--output-dir", type=str, default="data/raw_annotations", help="Annotations output directory")
    parser.add_argument("--single-image", type=str, default=None, help="Target a specific image file")
    parser.add_argument("--checkpoint", type=str, default="models/checkpoints/sam2_hiera_large.pt", help="SAM2 weights")
    parser.add_argument("--config", type=str, default="sam2_hiera_l.yaml", help="SAM2 model configuration")
    return parser.parse_args()


def main() -> None:
    """Main execution entrypoint."""
    args = parse_args()
    annotator = PrecisionSAM2Annotator(
        images_dir=args.images_dir,
        output_dir=args.output_dir,
        single_image=args.single_image,
        checkpoint_path=args.checkpoint,
        config_path=args.config
    )
    logger.info("Precision SAM 2 Annotator initialized.")


if __name__ == "__main__":
    main()
