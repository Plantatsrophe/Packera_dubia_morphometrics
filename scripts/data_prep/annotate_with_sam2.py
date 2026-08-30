#!/usr/bin/env python3
"""
Script: scripts/data_prep/annotate_with_sam2.py
Project: Packera dubia Morphometrics Pipeline
Author: Botanical AI & Image Processing Specialist

Interactive click-to-annotate tool powered by SAM 2.
Left-Click = Positive prompt point (green)
Right-Click = Negative prompt point (red)
Keys 0-6 = Assign botanical class
's' = Save YOLO label file and advance to next voucher
'c' = Clear current candidate points
'u' = Undo last saved instance
'q' / Esc = Quit
"""

import os
import sys
import glob
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Any

# Add project root directory to sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import cv2
import numpy as np
import torch

from sam2.build_sam import build_sam2
from sam2.sam2_image_predictor import SAM2ImagePredictor

# Logging setup
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("SAM2Annotator")

# Class Schema Mapping
CLASS_NAMES = [
    "basal_leaf_blade",    # 0
    "leaf_petiole",         # 1
    "cauline_leaf",         # 2
    "cauline_stem",         # 3
    "root_rhizome",         # 4
    "basal_rosette_clump",  # 5
    "capitulum"             # 6
]

CLASS_COLORS = {
    0: (0, 220, 0),      # basal_leaf_blade: Green
    1: (100, 255, 100),  # leaf_petiole: Light Green
    2: (0, 200, 255),    # cauline_leaf: Yellow-Green
    3: (0, 140, 255),    # cauline_stem: Orange
    4: (50, 50, 180),    # root_rhizome: Brown / Red
    5: (0, 100, 50),     # basal_rosette_clump: Dark Green
    6: (0, 230, 255)     # capitulum: Yellow
}

class InteractiveSAM2Annotator:
    def __init__(
        self,
        images_dir: str = "data/raw_vouchers",
        output_dir: str = "data/raw_annotations",
        checkpoint_path: str = "models/checkpoints/sam2_hiera_large.pt",
        config_path: str = "sam2_hiera_l.yaml"
    ):
        self.images_dir = Path(images_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.image_files = sorted(
            list(self.images_dir.glob("*.jpg")) + 
            list(self.images_dir.glob("*.jpeg")) + 
            list(self.images_dir.glob("*.png"))
        )

        if not self.image_files:
            raise FileNotFoundError(f"No voucher images found in {self.images_dir}!")

        # Initialize SAM 2 Predictor
        device = "cuda" if torch.cuda.is_available() else "cpu"
        logger.info(f"Loading SAM 2 model onto device: {device}...")
        sam2_model = build_sam2(config_path, checkpoint_path, device=device)
        self.predictor = SAM2ImagePredictor(sam2_model)

        # State Variables
        self.current_idx = 0
        self.current_img = None
        self.display_img = None
        self.orig_h = 0
        self.orig_w = 0
        self.scale_factor = 1.0

        # Points for current candidate
        self.prompt_points: List[List[int]] = []
        self.prompt_labels: List[int] = [] # 1 = positive, 0 = negative
        self.current_candidate_mask = None

        # Saved instances for current sheet
        self.saved_instances: List[Dict[str, Any]] = []

    def mouse_callback(self, event, x, y, flags, param):
        # Map display coordinates back to original native-DPI coordinates
        orig_x = int(x / self.scale_factor)
        orig_y = int(y / self.scale_factor)

        if event == cv2.EVENT_LBUTTONDOWN:
            # Left click = Positive foreground prompt point
            self.prompt_points.append([orig_x, orig_y])
            self.prompt_labels.append(1)
            self.update_candidate_mask()

        elif event == cv2.EVENT_RBUTTONDOWN:
            # Right click = Negative background exclusion point
            self.prompt_points.append([orig_x, orig_y])
            self.prompt_labels.append(0)
            self.update_candidate_mask()

    def update_candidate_mask(self):
        if not self.prompt_points:
            self.current_candidate_mask = None
            return

        pts = np.array(self.prompt_points)
        lbls = np.array(self.prompt_labels)

        masks, scores, _ = self.predictor.predict(
            point_coords=pts,
            point_labels=lbls,
            multimask_output=True
        )
        # Select highest confidence mask
        self.current_candidate_mask = masks[np.argmax(scores)]

    def mask_to_yolo_polygon(self, mask: np.ndarray) -> List[float]:
        contours, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            return []
        c = max(contours, key=cv2.contourArea)
        epsilon = 0.002 * cv2.arcLength(c, True)
        approx = cv2.approxPolyDP(c, epsilon, True).reshape(-1, 2)
        if len(approx) < 3:
            return []

        normalized = []
        for pt in approx:
            normalized.extend([round(pt[0] / self.orig_w, 6), round(pt[1] / self.orig_h, 6)])
        return normalized

    def save_current_sheet_annotations(self):
        cat_num = self.image_files[self.current_idx].stem
        label_file = self.output_dir / f"{cat_num}.txt"

        with open(label_file, "w") as f:
            for inst in self.saved_instances:
                poly = inst["polygon"]
                if len(poly) >= 6:
                    line = f"{inst['class_id']} " + " ".join(f"{v:.6f}" for v in poly)
                    f.write(line + "\n")

        logger.info(f"Saved {len(self.saved_instances)} instances to {label_file}")

    def render_view(self) -> np.ndarray:
        render = self.current_img.copy()

        # Render saved instances
        for inst in self.saved_instances:
            color = CLASS_COLORS[inst["class_id"]]
            pts = np.array(inst["polygon"]).reshape(-1, 2)
            pts[:, 0] *= self.orig_w
            pts[:, 1] *= self.orig_h
            pts = pts.astype(np.int32)

            cv2.polylines(render, [pts], isClosed=True, color=color, thickness=3)
            # Semi-transparent overlay
            overlay = render.copy()
            cv2.fillPoly(overlay, [pts], color=color)
            cv2.addWeighted(overlay, 0.35, render, 0.65, 0, render)

        # Render current live candidate mask
        if self.current_candidate_mask is not None:
            mask_color = (255, 255, 0) # Cyan preview
            contours, _ = cv2.findContours(self.current_candidate_mask.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            cv2.drawContours(render, contours, -1, mask_color, 3)

        # Render prompt points
        for pt, lbl in zip(self.prompt_points, self.prompt_labels):
            color = (0, 255, 0) if lbl == 1 else (0, 0, 255)
            cv2.circle(render, (pt[0], pt[1]), 8, color, -1)
            cv2.circle(render, (pt[0], pt[1]), 10, (255, 255, 255), 2)

        # Resize for display screen
        target_display_h = 1000
        self.scale_factor = target_display_h / self.orig_h
        display_w = int(self.orig_w * self.scale_factor)
        display_h = target_display_h
        resized = cv2.resize(render, (display_w, display_h), interpolation=cv2.INTER_LINEAR)

        # Add Status HUD
        hud_text = f"[{self.current_idx + 1}/{len(self.image_files)}] {self.image_files[self.current_idx].name} | Saved: {len(self.saved_instances)}"
        cv2.putText(resized, hud_text, (20, 35), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 0, 0), 3)
        cv2.putText(resized, hud_text, (20, 35), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (255, 255, 255), 1)

        key_legend = "0:Blade | 1:Petiole | 2:CaulineLeaf | 3:Stem | 4:Root | 5:Rosette | 6:Head | 's':Save | 'c':Clear"
        cv2.putText(resized, key_legend, (20, 70), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (0, 0, 0), 3)
        cv2.putText(resized, key_legend, (20, 70), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (0, 255, 255), 1)

        return resized

    def run(self):
        cv2.namedWindow("SAM2 Botanical Voucher Annotator", cv2.WINDOW_AUTOSIZE)
        cv2.setMouseCallback("SAM2 Botanical Voucher Annotator", self.mouse_callback)

        while self.current_idx < len(self.image_files):
            img_path = self.image_files[self.current_idx]
            cat_num = img_path.stem
            label_file = self.output_dir / f"{cat_num}.txt"

            # Check if already annotated
            if label_file.exists():
                logger.info(f"Skipping {cat_num} (already annotated).")
                self.current_idx += 1
                continue

            logger.info(f"Loading {img_path.name}...")
            self.current_img = cv2.imread(str(img_path))
            self.orig_h, self.orig_w = self.current_img.shape[:2]

            # Set SAM 2 image features once per sheet
            img_rgb = cv2.cvtColor(self.current_img, cv2.COLOR_BGR2RGB)
            self.predictor.set_image(img_rgb)

            self.saved_instances = []
            self.prompt_points = []
            self.prompt_labels = []
            self.current_candidate_mask = None

            while True:
                display = self.render_view()
                cv2.imshow("SAM2 Botanical Voucher Annotator", display)
                key = cv2.waitKey(20) & 0xFF

                # Class Assignment Hotkeys (0-6)
                if ord('0') <= key <= ord('6'):
                    if self.current_candidate_mask is not None:
                        class_id = int(chr(key))
                        poly = self.mask_to_yolo_polygon(self.current_candidate_mask)
                        if poly:
                            self.saved_instances.append({
                                "class_id": class_id,
                                "polygon": poly
                            })
                            logger.info(f"Assigned instance to '{CLASS_NAMES[class_id]}'")
                        # Reset candidate prompts
                        self.prompt_points = []
                        self.prompt_labels = []
                        self.current_candidate_mask = None

                # 's' or Enter -> Save and advance
                elif key == ord('s') or key == 13:
                    self.save_current_sheet_annotations()
                    self.current_idx += 1
                    break

                # 'c' -> Clear current prompt points
                elif key == ord('c'):
                    self.prompt_points = []
                    self.prompt_labels = []
                    self.current_candidate_mask = None

                # 'u' -> Undo last saved instance
                elif key == ord('u'):
                    if self.saved_instances:
                        removed = self.saved_instances.pop()
                        logger.info(f"Undid instance: {CLASS_NAMES[removed['class_id']]}")

                # 'q' or Esc -> Exit
                elif key == ord('q') or key == 27:
                    cv2.destroyAllWindows()
                    logger.info("Exiting annotator.")
                    return

        cv2.destroyAllWindows()
        logger.info("All vouchers processed!")

if __name__ == "__main__":
    annotator = InteractiveSAM2Annotator()
    annotator.run()