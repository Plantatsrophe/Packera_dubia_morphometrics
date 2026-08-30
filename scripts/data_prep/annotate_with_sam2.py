#!/usr/bin/env python3
"""
Script: scripts/data_prep/annotate_with_sam2.py
Project: Packera dubia Morphometrics Pipeline
Description: Precision Botanical Instance Segmentation Annotator using SAM 2 with
             multi-modal exclusion point support, bounding box constraints,
             and knife-cut boundary splitting for herbarium voucher curation.

Controls & Shortcuts:
  - Left-Click: Add Point Prompt (Green positive dot in INCLUDE mode, Red negative dot in EXCLUDE mode)
  - Right-Click: Add Negative Exclusion Point (Red dot, label 0) - Qt context menus disabled
  - Alt + Left-Click: Add Negative Exclusion Point (Red dot, label 0)
  - 'e': Toggle Exclusion Mode [MODE: INCLUDE (Green)] <-> [MODE: EXCLUDE (Red)]
  - Shift + Left-Drag: Draw Bounding Box Prompt (Constrains SAM 2 search space)
  - Ctrl + Left-Drag: Draw Knife / Cut Line across mask to sever roots from petioles
  - Keys 0-6: Assign Botanical Class to Candidate Mask
      0: basal_leaf_blade, 1: leaf_petiole, 2: cauline_leaf, 3: cauline_stem,
      4: root_rhizome, 5: basal_rosette_clump, 6: capitulum
  - 'u': Undo last saved instance
  - 'c': Clear active candidate prompts
  - Enter / 'v' / Shift+S: Save voucher annotations (.txt) and advance to next sheet
  - 'w' / 'a' / 's' / 'd' (or Arrow Keys): Pan viewport Up / Left / Down / Right
  - 'z' / 'x': Zoom In / Zoom Out around viewport
  - 'q' / Esc: Quit Annotator
"""

import os
import sys
import glob
import logging
import argparse
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional, Union

import cv2
import numpy as np
import torch

from sam2.build_sam import build_sam2
from sam2.sam2_image_predictor import SAM2ImagePredictor

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("SAM2Annotator")

def get_project_root() -> Path:
    """Dynamically resolves the project root directory."""
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "data").exists() or (parent / "models").exists() or (parent / ".git").exists():
            return parent
    return current.parents[1]

PROJECT_ROOT = get_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

CLASS_NAMES: List[str] = [
    "basal_leaf_blade",    # 0
    "leaf_petiole",        # 1
    "cauline_leaf",        # 2
    "cauline_stem",        # 3
    "root_rhizome",        # 4
    "basal_rosette_clump", # 5
    "capitulum"            # 6
]

CLASS_COLORS: Dict[int, Tuple[int, int, int]] = {
    0: (0, 220, 0),      # basal_leaf_blade: Vibrant Green
    1: (100, 255, 100),  # leaf_petiole: Mint Green
    2: (0, 200, 255),    # cauline_leaf: Yellow-Green
    3: (0, 140, 255),    # cauline_stem: Orange
    4: (50, 50, 200),    # root_rhizome: Red/Brown
    5: (0, 100, 50),     # basal_rosette_clump: Dark Forest Green
    6: (0, 230, 255)     # capitulum: Yellow
}


class PrecisionSAM2Annotator:
    """
    Interactive high-resolution herbarium voucher annotator powered by SAM 2.
    Features robust OS/OpenCV window backend configuration, multi-modal exclusion
    point placement, bounding box constraints, knife-cutter partitioning, and YOLO export.
    """
    def __init__(
        self,
        images_dir: Union[str, Path] = "data/raw_vouchers",
        output_dir: Union[str, Path] = "data/raw_annotations",
        single_image: Optional[Union[str, Path]] = None,
        checkpoint_path: Union[str, Path] = "models/checkpoints/sam2_hiera_large.pt",
        config_path: Union[str, Path] = "sam2_hiera_l.yaml",
        window_w: int = 1280,
        window_h: int = 960
    ):
        self.project_root = get_project_root()
        self.images_dir = Path(images_dir) if Path(images_dir).is_absolute() else self.project_root / images_dir
        self.output_dir = Path(output_dir) if Path(output_dir).is_absolute() else self.project_root / output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.window_w = window_w
        self.window_h = window_h

        if single_image:
            img_p = Path(single_image) if Path(single_image).is_absolute() else self.project_root / single_image
            if not img_p.exists():
                raise FileNotFoundError(f"Voucher image not found: {img_p}")
            self.image_files = [img_p]
        else:
            self.image_files = sorted(
                list(self.images_dir.glob("*.jpg")) + 
                list(self.images_dir.glob("*.jpeg")) + 
                list(self.images_dir.glob("*.png")) +
                list(self.images_dir.glob("*.tif")) +
                list(self.images_dir.glob("*.tiff"))
            )

        if not self.image_files:
            raise FileNotFoundError(f"No voucher images found in {self.images_dir}!")

        # Resolve checkpoint path
        ckpt_p = Path(checkpoint_path) if Path(checkpoint_path).is_absolute() else self.project_root / checkpoint_path
        if not ckpt_p.exists():
            alt_ckpt = self.project_root / "models" / "checkpoints" / "sam2_hiera_large.pt"
            if alt_ckpt.exists():
                ckpt_p = alt_ckpt

        device = "cuda" if torch.cuda.is_available() else "cpu"
        logger.info(f"Loading SAM 2 on device: {device} from {ckpt_p}...")
        cfg_name = Path(config_path).name
        sam2_model = build_sam2(cfg_name, str(ckpt_p), device=device)
        self.predictor = SAM2ImagePredictor(sam2_model)

        self.current_idx = 0
        self.current_img: Optional[np.ndarray] = None
        self.orig_h = 0
        self.orig_w = 0

        # Interactive Annotation State
        self.prompt_points: List[List[int]] = []
        self.prompt_labels: List[int] = []  # 1 = positive foreground, 0 = negative exclusion
        self.prompt_box: Optional[List[int]] = None  # [x1, y1, x2, y2]
        self.current_candidate_mask: Optional[np.ndarray] = None
        self.saved_instances: List[Dict[str, Any]] = []

        # Multi-modal mode state
        self.exclusion_mode: bool = False  # Toggled via 'e' key

        # Mouse Dragging State
        self.is_drawing_box = False
        self.is_drawing_cut_line = False
        self.drag_start = (0, 0)
        self.drag_current = (0, 0)
        self.cut_lines: List[Tuple[Tuple[int, int], Tuple[int, int]]] = []

        # Viewport Zoom & Pan State
        self.zoom_level = 1.0
        self.pan_offset = [0, 0]  # [x, y] in original image space

    def get_orig_coords(self, display_x: int, display_y: int, window_w: int, window_h: int) -> Tuple[int, int]:
        """Maps window display coordinates back to original native-DPI coordinates."""
        if self.orig_w <= 0 or self.orig_h <= 0 or window_w <= 0 or window_h <= 0:
            return 0, 0
        view_w = max(1, int(self.orig_w / self.zoom_level))
        view_h = max(1, int(self.orig_h / self.zoom_level))
        
        orig_x = int(self.pan_offset[0] + (display_x / window_w) * view_w)
        orig_y = int(self.pan_offset[1] + (display_y / window_h) * view_h)
        
        orig_x = max(0, min(self.orig_w - 1, orig_x))
        orig_y = max(0, min(self.orig_h - 1, orig_y))
        return orig_x, orig_y

    def mouse_callback(self, event, x, y, flags, param):
        """
        HighGUI Mouse callback supporting bounding boxes, knife cut slicing,
        native right-click exclusion, Alt+Left click exclusion, and mode-based clicks.
        """
        window_w, window_h = param
        orig_x, orig_y = self.get_orig_coords(x, y, window_w, window_h)

        # 1. Bounding Box Prompt Drag (Shift + Left Click Drag)
        if (flags & cv2.EVENT_FLAG_SHIFTKEY) and event == cv2.EVENT_LBUTTONDOWN:
            self.is_drawing_box = True
            self.drag_start = (orig_x, orig_y)
            self.drag_current = (orig_x, orig_y)

        elif self.is_drawing_box and event == cv2.EVENT_MOUSEMOVE:
            self.drag_current = (orig_x, orig_y)

        elif self.is_drawing_box and event == cv2.EVENT_LBUTTONUP:
            self.is_drawing_box = False
            x1, x2 = min(self.drag_start[0], orig_x), max(self.drag_start[0], orig_x)
            y1, y2 = min(self.drag_start[1], orig_y), max(self.drag_start[1], orig_y)
            if (x2 - x1) > 10 and (y2 - y1) > 10:
                self.prompt_box = [x1, y1, x2, y2]
                logger.info(f"Bounding box prompt set: [{x1}, {y1}, {x2}, {y2}]")
                self.update_candidate_mask()

        # 2. Knife / Cut Line Drag (Ctrl + Left Click Drag OR Middle Click Drag)
        elif ((flags & cv2.EVENT_FLAG_CTRLKEY) and event == cv2.EVENT_LBUTTONDOWN) or (event == cv2.EVENT_MBUTTONDOWN):
            self.is_drawing_cut_line = True
            self.drag_start = (orig_x, orig_y)
            self.drag_current = (orig_x, orig_y)

        elif self.is_drawing_cut_line and event == cv2.EVENT_MOUSEMOVE:
            self.drag_current = (orig_x, orig_y)

        elif self.is_drawing_cut_line and (event == cv2.EVENT_LBUTTONUP or event == cv2.EVENT_MBUTTONUP):
            self.is_drawing_cut_line = False
            self.apply_knife_cut(self.drag_start, (orig_x, orig_y))

        # 3. Alt + Left-Click (Negative Exclusion Point)
        elif (flags & cv2.EVENT_FLAG_ALTKEY) and event == cv2.EVENT_LBUTTONDOWN:
            self.prompt_points.append([orig_x, orig_y])
            self.prompt_labels.append(0)  # Negative Background Exclusion
            logger.info(f"Placed Exclusion Point (Alt+Click) at ({orig_x}, {orig_y})")
            self.update_candidate_mask()

        # 4. Native Right-Click (Negative Exclusion Point)
        elif event == cv2.EVENT_RBUTTONDOWN:
            self.prompt_points.append([orig_x, orig_y])
            self.prompt_labels.append(0)  # Negative Background Exclusion
            logger.info(f"Placed Exclusion Point (Right-Click) at ({orig_x}, {orig_y})")
            self.update_candidate_mask()

        # 5. Standard Left-Click Point (Mode-Dependent)
        elif event == cv2.EVENT_LBUTTONDOWN and not (flags & cv2.EVENT_FLAG_SHIFTKEY) and not (flags & cv2.EVENT_FLAG_CTRLKEY) and not (flags & cv2.EVENT_FLAG_ALTKEY):
            if self.exclusion_mode:
                self.prompt_points.append([orig_x, orig_y])
                self.prompt_labels.append(0)  # Negative Background Exclusion
                logger.info(f"Placed Exclusion Point (Exclusion Mode) at ({orig_x}, {orig_y})")
            else:
                self.prompt_points.append([orig_x, orig_y])
                self.prompt_labels.append(1)  # Positive Foreground Point
                logger.info(f"Placed Include Point at ({orig_x}, {orig_y})")
            self.update_candidate_mask()

    def update_candidate_mask(self):
        """Runs SAM 2 predictor inference on the active points, labels, and bounding box."""
        if not self.prompt_points and self.prompt_box is None:
            self.current_candidate_mask = None
            return

        pts = np.array(self.prompt_points) if self.prompt_points else None
        lbls = np.array(self.prompt_labels) if self.prompt_labels else None
        box = np.array(self.prompt_box) if self.prompt_box is not None else None

        try:
            masks, scores, _ = self.predictor.predict(
                point_coords=pts,
                point_labels=lbls,
                box=box,
                multimask_output=True
            )
            self.current_candidate_mask = masks[np.argmax(scores)].copy()

            # Re-apply any active cut lines
            if self.cut_lines and self.current_candidate_mask is not None:
                mask_uint8 = self.current_candidate_mask.astype(np.uint8)
                for pt1, pt2 in self.cut_lines:
                    cv2.line(mask_uint8, pt1, pt2, 0, thickness=4)
                self.current_candidate_mask = (mask_uint8 > 0)
        except Exception as e:
            logger.error(f"Error predicting SAM 2 mask: {e}")

    def apply_knife_cut(self, pt1: Tuple[int, int], pt2: Tuple[int, int]):
        """
        Physically carves a 4-pixel zero-mask cut line across the candidate mask
        and retains only the connected component containing the active positive prompt point.
        """
        if self.current_candidate_mask is None:
            logger.warning("No candidate mask active to apply knife cut.")
            return

        self.cut_lines.append((pt1, pt2))
        mask_uint8 = self.current_candidate_mask.astype(np.uint8)
        cv2.line(mask_uint8, pt1, pt2, 0, thickness=4)
        self.current_candidate_mask = (mask_uint8 > 0)
        
        # Keep only the connected component containing the positive point prompt
        num_labels, labels_im = cv2.connectedComponents(self.current_candidate_mask.astype(np.uint8))
        if num_labels > 2 and self.prompt_points:
            pos_pts = [p for p, l in zip(self.prompt_points, self.prompt_labels) if l == 1]
            if pos_pts:
                target_x, target_y = pos_pts[0]
                if 0 <= target_y < labels_im.shape[0] and 0 <= target_x < labels_im.shape[1]:
                    target_label = labels_im[target_y, target_x]
                    if target_label > 0:
                        self.current_candidate_mask = (labels_im == target_label)
        logger.info(f"Knife cut applied between {pt1} and {pt2}!")

    def mask_to_yolo_polygon(self, mask: np.ndarray) -> List[float]:
        """
        Extracts normalized polygon coordinates (0.0 to 1.0) from a binary segmentation mask.
        """
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
            norm_x = max(0.0, min(1.0, round(float(pt[0]) / self.orig_w, 6)))
            norm_y = max(0.0, min(1.0, round(float(pt[1]) / self.orig_h, 6)))
            normalized.extend([norm_x, norm_y])
        return normalized

    def save_current_sheet(self):
        """Exports all saved instances into a normalized YOLO polygon .txt file."""
        cat_num = self.image_files[self.current_idx].stem
        label_file = self.output_dir / f"{cat_num}.txt"
        with open(label_file, "w", encoding="utf-8") as f:
            for inst in self.saved_instances:
                poly = inst["polygon"]
                if len(poly) >= 6:
                    line = f"{inst['class_id']} " + " ".join(f"{v:.6f}" for v in poly)
                    f.write(line + "\n")
        logger.info(f"Successfully saved {len(self.saved_instances)} instances to {label_file}")

    def render_display(self, window_w: int = 1280, window_h: int = 960) -> np.ndarray:
        """Renders the interactive viewport, saved overlays, candidate mask, prompts, and HUD."""
        if self.current_img is None:
            return np.zeros((window_h, window_w, 3), dtype=np.uint8)

        # Calculate viewport slice for Zoom/Pan
        view_w = max(1, int(self.orig_w / self.zoom_level))
        view_h = max(1, int(self.orig_h / self.zoom_level))
        
        # Clamp pan offsets
        self.pan_offset[0] = max(0, min(self.orig_w - view_w, self.pan_offset[0]))
        self.pan_offset[1] = max(0, min(self.orig_h - view_h, self.pan_offset[1]))
        
        x1, y1 = self.pan_offset[0], self.pan_offset[1]
        x2, y2 = x1 + view_w, y1 + view_h
        
        crop = self.current_img[y1:y2, x1:x2].copy()

        # Render Saved Polygons
        for inst in self.saved_instances:
            color = CLASS_COLORS.get(inst["class_id"], (0, 220, 0))
            pts = np.array(inst["polygon"]).reshape(-1, 2)
            pts[:, 0] = (pts[:, 0] * self.orig_w) - x1
            pts[:, 1] = (pts[:, 1] * self.orig_h) - y1
            pts = pts.astype(np.int32)
            cv2.polylines(crop, [pts], True, color, 3)
            overlay = crop.copy()
            cv2.fillPoly(overlay, [pts], color)
            cv2.addWeighted(overlay, 0.35, crop, 0.65, 0, crop)

        # Render Active Candidate Mask
        if self.current_candidate_mask is not None:
            mask_crop = self.current_candidate_mask[y1:y2, x1:x2]
            if np.any(mask_crop):
                contours, _ = cv2.findContours(mask_crop.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(crop, contours, -1, (255, 255, 0), 3)  # Cyan boundary
                overlay = crop.copy()
                cv2.drawContours(overlay, contours, -1, (255, 255, 0), -1)
                cv2.addWeighted(overlay, 0.20, crop, 0.80, 0, crop)

        # Render Bounding Box Prompt
        if self.prompt_box is not None:
            bx1, by1, bx2, by2 = self.prompt_box
            cv2.rectangle(crop, (bx1 - x1, by1 - y1), (bx2 - x1, by2 - y1), (255, 120, 0), 2)

        # Render In-Progress Drag Box
        if self.is_drawing_box:
            cv2.rectangle(crop, (self.drag_start[0] - x1, self.drag_start[1] - y1),
                                (self.drag_current[0] - x1, self.drag_current[1] - y1), (255, 255, 255), 2)

        # Render In-Progress Knife Cut Line
        if self.is_drawing_cut_line:
            cv2.line(crop, (self.drag_start[0] - x1, self.drag_start[1] - y1),
                           (self.drag_current[0] - x1, self.drag_current[1] - y1), (0, 0, 255), 3)

        # Render Prompt Points
        for pt, lbl in zip(self.prompt_points, self.prompt_labels):
            px, py = pt[0] - x1, pt[1] - y1
            color = (0, 255, 0) if lbl == 1 else (0, 0, 255)
            cv2.circle(crop, (px, py), 7, color, -1)
            cv2.circle(crop, (px, py), 9, (255, 255, 255), 2)

        # Resize cropped viewport to window display resolution
        display = cv2.resize(crop, (window_w, window_h), interpolation=cv2.INTER_LINEAR)

        # Top HUD Banner Background
        hud_bg = display[:105, :].copy()
        cv2.rectangle(display, (0, 0), (window_w, 105), (20, 20, 20), -1)
        cv2.addWeighted(hud_bg, 0.25, display[:105, :], 0.75, 0, display[:105, :])
        cv2.line(display, (0, 105), (window_w, 105), (80, 80, 80), 1)

        # 1. Main Status Line & Mode Indicator
        cat_num = self.image_files[self.current_idx].stem
        status_text = f"[{self.current_idx + 1}/{len(self.image_files)}] {cat_num}  |  Zoom: {self.zoom_level:.1f}x  |  Saved: {len(self.saved_instances)}"
        cv2.putText(display, status_text, (20, 32), cv2.FONT_HERSHEY_SIMPLEX, 0.70, (255, 255, 255), 2, cv2.LINE_AA)

        # Mode Indicator (Vibrant Green for Include, Bright Red for Exclude)
        if self.exclusion_mode:
            mode_text = "[MODE: EXCLUDE (Red)]"
            mode_color = (0, 0, 255)  # Bright Red
        else:
            mode_text = "[MODE: INCLUDE (Green)]"
            mode_color = (0, 255, 0)  # Vibrant Green

        mode_x = window_w - 390
        cv2.putText(display, mode_text, (mode_x, 32), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (0, 0, 0), 4, cv2.LINE_AA)
        cv2.putText(display, mode_text, (mode_x, 32), cv2.FONT_HERSHEY_SIMPLEX, 0.75, mode_color, 2, cv2.LINE_AA)

        # 2. Controls & Hotkeys Help
        help_line1 = "L-Click: Point  |  R-Click / Alt+Click: Exclude  |  'e': Toggle Mode  |  Shift+Drag: Box  |  Ctrl+Drag: Cut"
        cv2.putText(display, help_line1, (20, 62), cv2.FONT_HERSHEY_SIMPLEX, 0.48, (220, 220, 220), 1, cv2.LINE_AA)

        help_line2 = "0-6: Assign Class  |  'u': Undo  |  'c': Clear  |  Enter / 'v': Save Sheet  |  WASD: Pan  |  'z'/'x': Zoom"
        cv2.putText(display, help_line2, (20, 90), cv2.FONT_HERSHEY_SIMPLEX, 0.48, (0, 230, 255), 1, cv2.LINE_AA)

        return display

    def run(self):
        """Main annotation event loop with explicit HighGUI backend configuration."""
        window_name = "SAM 2 Precision Botanical Annotator"
        # Explicitly configure window backend to disable Qt context menu and toolbars
        cv2.namedWindow(window_name, cv2.WINDOW_NORMAL | cv2.WINDOW_GUI_NORMAL)
        cv2.resizeWindow(window_name, self.window_w, self.window_h)
        cv2.setMouseCallback(window_name, self.mouse_callback, param=(self.window_w, self.window_h))

        while self.current_idx < len(self.image_files):
            img_path = self.image_files[self.current_idx]
            cat_num = img_path.stem
            label_file = self.output_dir / f"{cat_num}.txt"

            if label_file.exists() and len(self.image_files) > 1:
                logger.info(f"Skipping {cat_num} (already annotated: {label_file.name}).")
                self.current_idx += 1
                continue

            logger.info(f"Loading and encoding {img_path.name} with SAM 2...")
            self.current_img = cv2.imread(str(img_path))
            if self.current_img is None:
                logger.error(f"Failed to read image: {img_path}")
                self.current_idx += 1
                continue

            self.orig_h, self.orig_w = self.current_img.shape[:2]

            img_rgb = cv2.cvtColor(self.current_img, cv2.COLOR_BGR2RGB)
            self.predictor.set_image(img_rgb)

            self.saved_instances = []
            self.prompt_points = []
            self.prompt_labels = []
            self.prompt_box = None
            self.current_candidate_mask = None
            self.cut_lines = []
            self.zoom_level = 1.0
            self.pan_offset = [0, 0]
            self.exclusion_mode = False

            while True:
                display = self.render_display(self.window_w, self.window_h)
                cv2.imshow(window_name, display)
                key = cv2.waitKey(20) & 0xFF

                # Class Hotkeys (0-6)
                if ord('0') <= key <= ord('6'):
                    if self.current_candidate_mask is not None:
                        class_id = int(chr(key))
                        poly = self.mask_to_yolo_polygon(self.current_candidate_mask)
                        if poly:
                            self.saved_instances.append({"class_id": class_id, "polygon": poly})
                            logger.info(f"Assigned instance -> Class {class_id}: {CLASS_NAMES[class_id]} ({len(poly)//2} vertices)")
                        else:
                            logger.warning("Candidate mask could not be converted to polygon (too small or empty).")
                        self.prompt_points = []
                        self.prompt_labels = []
                        self.prompt_box = None
                        self.current_candidate_mask = None
                        self.cut_lines = []

                # Exclusion Mode Toggle ('e')
                elif key == ord('e'):
                    self.exclusion_mode = not self.exclusion_mode
                    state_str = "EXCLUDE (Red, label 0)" if self.exclusion_mode else "INCLUDE (Green, label 1)"
                    logger.info(f"Exclusion mode toggled -> {state_str}")

                # Zoom Controls ('z' = Zoom In, 'x' = Zoom Out)
                elif key == ord('z'):
                    self.zoom_level = min(8.0, self.zoom_level * 1.3)
                elif key == ord('x'):
                    self.zoom_level = max(1.0, self.zoom_level / 1.3)

                # Pan Controls (WASD & Arrow keys)
                elif key in (ord('w'), 82, 0):  # Up
                    self.pan_offset[1] -= int(200 / self.zoom_level)
                elif key in (ord('s'), 84, 1):  # Down
                    self.pan_offset[1] += int(200 / self.zoom_level)
                elif key in (ord('a'), 81, 2):  # Left
                    self.pan_offset[0] -= int(200 / self.zoom_level)
                elif key in (ord('d'), 83, 3):  # Right
                    self.pan_offset[0] += int(200 / self.zoom_level)

                # Save and Advance (Enter or 'v' or 'S')
                elif key in (13, ord('v'), ord('S')):
                    self.save_current_sheet()
                    self.current_idx += 1
                    break

                # Clear Active Candidate Prompts ('c')
                elif key == ord('c'):
                    self.prompt_points = []
                    self.prompt_labels = []
                    self.prompt_box = None
                    self.current_candidate_mask = None
                    self.cut_lines = []
                    logger.info("Cleared active candidate prompts.")

                # Undo Last Saved Instance ('u')
                elif key == ord('u'):
                    if self.saved_instances:
                        removed = self.saved_instances.pop()
                        logger.info(f"Undid instance: Class {removed['class_id']} ({CLASS_NAMES[removed['class_id']]})")

                # Quit ('q' or Esc)
                elif key == ord('q') or key == 27:
                    logger.info("Exiting annotator session.")
                    cv2.destroyAllWindows()
                    return

        cv2.destroyAllWindows()
        logger.info("All vouchers annotated or queue completed.")


def parse_args():
    parser = argparse.ArgumentParser(description="SAM 2 Precision Botanical Annotator for Packera dubia")
    parser.add_argument("--images-dir", type=str, default="data/raw_vouchers", help="Directory containing raw voucher images")
    parser.add_argument("--output-dir", type=str, default="data/raw_annotations", help="Directory to save YOLO polygon annotations")
    parser.add_argument("--image", type=str, default=None, help="Path to single voucher image to annotate")
    parser.add_argument("--checkpoint", type=str, default="models/checkpoints/sam2_hiera_large.pt", help="SAM 2 model checkpoint path")
    parser.add_argument("--config", type=str, default="sam2_hiera_l.yaml", help="SAM 2 config filename")
    parser.add_argument("--width", type=int, default=1280, help="Display window width")
    parser.add_argument("--height", type=int, default=960, help="Display window height")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    annotator = PrecisionSAM2Annotator(
        images_dir=args.images_dir,
        output_dir=args.output_dir,
        single_image=args.image,
        checkpoint_path=args.checkpoint,
        config_path=args.config,
        window_w=args.width,
        window_h=args.height
    )
    annotator.run()