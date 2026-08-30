#!/usr/bin/env python3
"""
Script: scripts/annotate_with_sam2.py
Project: Packera dubia Morphometrics Pipeline
Description: Production-ready Interactive Botanical Instance Segmentation Annotator
             powered by Segment Anything Model 2 (SAM 2) with multi-modal boundary
             controls (bounding boxes, exclusion points, freehand lasso, knife cutting,
             multi-point spines, zoom/pan navigation, voucher advancement, reloading,
             and back navigation) for the 7-class botanical phenotyping schema.

Controls & Shortcuts:
  - Left-Click: Add Point Prompt (Green foreground point in INCLUDE mode, Red exclusion point in EXCLUDE mode)
  - Right-Click: Add Negative Exclusion Point (Red dot, label 0) - OS popup menus permanently suppressed
  - Alt + Left-Click: Add Negative Exclusion Point (Red dot, label 0)
  - 'e': Toggle Exclusion Mode [MODE: INCLUDE (Green)] <-> [MODE: EXCLUDE (Red)]
  - 'l': Toggle Freehand Lasso Mode [MODE: LASSO (Cyan)]
  - Shift + Left-Drag: Draw Bounding Box Prompt (Constrains SAM 2 search space)
  - Ctrl + Left-Drag: Draw Knife / Cut Line across mask to cleanly sever roots from petioles
  - Keys 0-6: Assign Botanical Class to Candidate Mask
      0: basal_leaf_blade, 1: leaf_petiole, 2: cauline_leaf, 3: cauline_stem,
      4: root_rhizome, 5: basal_rosette_clump, 6: capitulum
  - 'u': Undo last saved instance
  - 'c': Clear active candidate prompts, boxes, cut lines, and lasso paths
  - 's' / Enter / 'v': Save voucher annotations (.txt) and advance to next sheet
  - 'b': Go back to previous voucher in queue (reloads existing annotations)
  - 'n': Skip forward to next voucher (without saving)
  - 'w' / 'a' / 's' / 'd' (or Arrow Keys): Pan viewport Up / Left / Down / Right
  - 'z' / 'x': Zoom In / Zoom Out centered on mouse/viewport (1.0x to 6.0x)
  - 'q' / Esc: Quit Annotator cleanly
"""

import os
import sys
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
    return current.parents[0] if len(current.parents) > 0 else current


PROJECT_ROOT = get_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# 7-Class Botanical Schema
CLASS_NAMES: List[str] = [
    "basal_leaf_blade",     # 0: Expanded laminar portion of basal leaves
    "leaf_petiole",         # 1: Foliar stalk connecting caudex to basal blade
    "cauline_leaf",         # 2: Sessile, lyrately-pinnatifid stem leaves
    "cauline_stem",         # 3: Vertical flowering stalk / scape / peduncle
    "root_rhizome",         # 4: Subterranean fibrous roots and rhizomes
    "basal_rosette_clump",  # 5: Dense overlapping rosette center / crown
    "capitulum"             # 6: Inflorescence head / involucre / phyllaries
]

# High-contrast BGR colors for visualization
CLASS_COLORS: Dict[int, Tuple[int, int, int]] = {
    0: (0, 220, 0),       # basal_leaf_blade: Vibrant Green
    1: (100, 255, 100),   # leaf_petiole: Mint Green
    2: (0, 200, 255),     # cauline_leaf: Yellow-Green
    3: (0, 140, 255),     # cauline_stem: Orange
    4: (50, 50, 200),     # root_rhizome: Red/Brown
    5: (0, 100, 50),      # basal_rosette_clump: Dark Forest Green
    6: (0, 230, 255)      # capitulum: Yellow
}


class PrecisionSAM2Annotator:
    """
    Production-ready interactive botanical annotator integrating SAM 2 with
    multi-modal boundary controls, freehand lasso, knife cutting, zoom/pan navigation,
    and normalized YOLO polygon dataset generation.
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
            extensions = ("*.jpg", "*.jpeg", "*.png", "*.tif", "*.tiff", "*.JPG", "*.JPEG", "*.PNG")
            files = []
            for ext in extensions:
                files.extend(list(self.images_dir.glob(ext)))
            self.image_files = sorted(list(set(files)))

        if not self.image_files:
            raise FileNotFoundError(f"No voucher images found in {self.images_dir}!")

        # Resolve checkpoint and config paths
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

        # Multi-modal mode states
        self.exclusion_mode: bool = False  # Toggled via 'e'
        self.lasso_mode: bool = False      # Toggled via 'l'

        # Dragging & Interactive Tool State
        self.is_drawing_box: bool = False
        self.is_drawing_cut_line: bool = False
        self.is_drawing_lasso: bool = False
        self.drag_start: Tuple[int, int] = (0, 0)
        self.drag_current: Tuple[int, int] = (0, 0)
        self.cut_lines: List[Tuple[Tuple[int, int], Tuple[int, int]]] = []
        self.lasso_points: List[Tuple[int, int]] = []

        # Viewport Zoom & Pan State (1.0x to 6.0x native-DPI magnification)
        self.zoom_level: float = 1.0
        self.pan_offset: List[int] = [0, 0]  # [x, y] in original image space
        self.last_mouse_pos: Optional[Tuple[int, int]] = None

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

    def clamp_pan(self):
        """Ensures pan offsets keep the viewport strictly within original image boundaries."""
        if self.orig_w <= 0 or self.orig_h <= 0:
            return
        view_w = max(1, int(self.orig_w / self.zoom_level))
        view_h = max(1, int(self.orig_h / self.zoom_level))
        max_pan_x = max(0, self.orig_w - view_w)
        max_pan_y = max(0, self.orig_h - view_h)
        self.pan_offset[0] = max(0, min(max_pan_x, self.pan_offset[0]))
        self.pan_offset[1] = max(0, min(max_pan_y, self.pan_offset[1]))

    def zoom(self, factor: float, center_x: Optional[int] = None, center_y: Optional[int] = None):
        """
        Zooms viewport by the given factor, anchored around the specified center
        point in native image coordinates (or viewport center if None).
        Supports 1.0x to 6.0x native-DPI magnification.
        """
        old_zoom = self.zoom_level
        new_zoom = max(1.0, min(6.0, self.zoom_level * factor))
        if abs(new_zoom - old_zoom) < 1e-4:
            return

        if center_x is None or center_y is None:
            view_w = self.orig_w / old_zoom
            view_h = self.orig_h / old_zoom
            center_x = int(self.pan_offset[0] + view_w / 2.0)
            center_y = int(self.pan_offset[1] + view_h / 2.0)

        new_view_w = self.orig_w / new_zoom
        new_view_h = self.orig_h / new_zoom

        self.pan_offset[0] = int(center_x - (center_x - self.pan_offset[0]) * (old_zoom / new_zoom))
        self.pan_offset[1] = int(center_y - (center_y - self.pan_offset[1]) * (old_zoom / new_zoom))
        self.zoom_level = new_zoom
        self.clamp_pan()

    def mouse_callback(self, event, x, y, flags, param):
        """
        HighGUI Mouse callback supporting bounding boxes, knife cut slicing,
        freehand lasso, native right-click exclusion, Alt+Left click exclusion,
        and multi-point spine prompt placement.
        """
        window_w, window_h = param
        orig_x, orig_y = self.get_orig_coords(x, y, window_w, window_h)
        self.last_mouse_pos = (orig_x, orig_y)

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

        elif self.is_drawing_cut_line and (event in (cv2.EVENT_LBUTTONUP, cv2.EVENT_MBUTTONUP)):
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

        # 5. Freehand Lasso Tool (When Lasso Mode is active)
        elif self.lasso_mode and not (flags & cv2.EVENT_FLAG_SHIFTKEY) and not (flags & cv2.EVENT_FLAG_CTRLKEY) and not (flags & cv2.EVENT_FLAG_ALTKEY):
            if event == cv2.EVENT_LBUTTONDOWN:
                self.is_drawing_lasso = True
                self.lasso_points = [(orig_x, orig_y)]

            elif self.is_drawing_lasso and event == cv2.EVENT_MOUSEMOVE:
                self.lasso_points.append((orig_x, orig_y))

            elif self.is_drawing_lasso and event == cv2.EVENT_LBUTTONUP:
                self.is_drawing_lasso = False
                if len(self.lasso_points) >= 3:
                    mask = np.zeros((self.orig_h, self.orig_w), dtype=np.uint8)
                    pts = np.array(self.lasso_points, dtype=np.int32)
                    cv2.fillPoly(mask, [pts], 1)
                    self.current_candidate_mask = (mask > 0)
                    self.prompt_points = []
                    self.prompt_labels = []
                    self.prompt_box = None
                    self.cut_lines = []
                    logger.info(f"Freehand lasso contour converted to candidate mask ({len(self.lasso_points)} points)")
                self.lasso_points = []

        # 6. Standard Left-Click Point (Mode-Dependent: Positive Spine vs Exclusion)
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
        Carves a 4-pixel zero-mask cut line across the candidate mask
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

    def load_existing_annotations(self, cat_num: str):
        """Loads existing YOLO polygon annotations if previously saved."""
        label_file = self.output_dir / f"{cat_num}.txt"
        self.saved_instances = []
        if label_file.exists():
            with open(label_file, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 7:
                        class_id = int(parts[0])
                        poly = [float(v) for v in parts[1:]]
                        self.saved_instances.append({
                            "class_id": class_id,
                            "polygon": poly
                        })
            logger.info(f"Loaded {len(self.saved_instances)} existing instances for {cat_num}")

    def save_current_sheet(self):
        """
        Exports all saved instances into a normalized YOLO polygon .txt file.
        If 0 instances are saved, creates an empty .txt file to serve as a hard-negative sample.
        """
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

        self.clamp_pan()
        x1, y1 = self.pan_offset[0], self.pan_offset[1]
        x2, y2 = x1 + view_w, y1 + view_h

        crop = self.current_img[y1:y2, x1:x2].copy()

        # 1. Render Saved Polygons with high-contrast colored fills
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

        # 2. Render Active Candidate Mask (Cyan Boundary & Semi-transparent Fill)
        if self.current_candidate_mask is not None:
            mask_crop = self.current_candidate_mask[y1:y2, x1:x2]
            if np.any(mask_crop):
                contours, _ = cv2.findContours(mask_crop.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(crop, contours, -1, (255, 255, 0), 3)  # Bright Cyan
                overlay = crop.copy()
                cv2.drawContours(overlay, contours, -1, (255, 255, 0), -1)
                cv2.addWeighted(overlay, 0.25, crop, 0.75, 0, crop)

        # 3. Render Bounding Box Prompt
        if self.prompt_box is not None:
            bx1, by1, bx2, by2 = self.prompt_box
            cv2.rectangle(crop, (bx1 - x1, by1 - y1), (bx2 - x1, by2 - y1), (255, 120, 0), 2)

        # 4. Render In-Progress Drag Box
        if self.is_drawing_box:
            cv2.rectangle(
                crop,
                (self.drag_start[0] - x1, self.drag_start[1] - y1),
                (self.drag_current[0] - x1, self.drag_current[1] - y1),
                (255, 255, 255),
                2
            )

        # 5. Render In-Progress Knife Cut Line
        if self.is_drawing_cut_line:
            cv2.line(
                crop,
                (self.drag_start[0] - x1, self.drag_start[1] - y1),
                (self.drag_current[0] - x1, self.drag_current[1] - y1),
                (0, 0, 255),
                3
            )

        # 6. Render Freehand Lasso in progress
        if self.is_drawing_lasso and len(self.lasso_points) > 1:
            lasso_pts = np.array(self.lasso_points, dtype=np.int32)
            lasso_pts[:, 0] -= x1
            lasso_pts[:, 1] -= y1
            cv2.polylines(crop, [lasso_pts], False, (255, 255, 0), 2)

        # 7. Render Prompt Points (Positive = Green, Negative Exclusion = Red)
        for pt, lbl in zip(self.prompt_points, self.prompt_labels):
            px, py = pt[0] - x1, pt[1] - y1
            color = (0, 255, 0) if lbl == 1 else (0, 0, 255)
            cv2.circle(crop, (px, py), 7, color, -1)
            cv2.circle(crop, (px, py), 9, (255, 255, 255), 2)

        # Resize cropped viewport to target window display resolution
        display = cv2.resize(crop, (window_w, window_h), interpolation=cv2.INTER_LINEAR)

        # Top HUD Banner Background
        hud_h = 105
        hud_bg = display[:hud_h, :].copy()
        cv2.rectangle(display, (0, 0), (window_w, hud_h), (20, 20, 20), -1)
        cv2.addWeighted(hud_bg, 0.25, display[:hud_h, :], 0.75, 0, display[:hud_h, :])
        cv2.line(display, (0, hud_h), (window_w, hud_h), (80, 80, 80), 1)

        # 1. Main Status Line & Mode Indicator
        cat_num = self.image_files[self.current_idx].stem
        status_text = f"[{self.current_idx + 1}/{len(self.image_files)}] {cat_num}  |  Zoom: {self.zoom_level:.1f}x  |  Saved: {len(self.saved_instances)}"
        cv2.putText(display, status_text, (20, 32), cv2.FONT_HERSHEY_SIMPLEX, 0.70, (255, 255, 255), 2, cv2.LINE_AA)

        # Mode Indicator Banner
        if self.lasso_mode:
            mode_text = "[MODE: LASSO (Cyan)]"
            mode_color = (255, 255, 0)  # Bright Cyan
        elif self.exclusion_mode:
            mode_text = "[MODE: EXCLUDE (Red)]"
            mode_color = (0, 0, 255)    # Bright Red
        else:
            mode_text = "[MODE: INCLUDE (Green)]"
            mode_color = (0, 255, 0)    # Vibrant Green

        mode_x = window_w - 410
        cv2.putText(display, mode_text, (mode_x, 32), cv2.FONT_HERSHEY_SIMPLEX, 0.75, (0, 0, 0), 4, cv2.LINE_AA)
        cv2.putText(display, mode_text, (mode_x, 32), cv2.FONT_HERSHEY_SIMPLEX, 0.75, mode_color, 2, cv2.LINE_AA)

        # 2. Controls & Hotkeys Help
        help_line1 = "L-Click: Point | R-Click / Alt+Click: Exclude | 'e': Mode | Shift+Drag: Box | Ctrl+Drag: Cut | 'l': Lasso"
        cv2.putText(display, help_line1, (20, 62), cv2.FONT_HERSHEY_SIMPLEX, 0.48, (220, 220, 220), 1, cv2.LINE_AA)

        help_line2 = "0-6: Class | 'u': Undo | 'c': Clear | 's'/Enter: Save | 'b': Back | 'n': Skip | WASD: Pan | 'z'/'x': Zoom"
        cv2.putText(display, help_line2, (20, 90), cv2.FONT_HERSHEY_SIMPLEX, 0.48, (0, 230, 255), 1, cv2.LINE_AA)

        return display

    def run(self):
        """Main annotation event loop with explicit HighGUI backend configuration."""
        window_name = "SAM 2 Precision Botanical Annotator"
        # Configure window backend to disable Qt scroll-hand drag panning, context menu, and toolbars
        cv2.namedWindow(window_name, cv2.WINDOW_AUTOSIZE | cv2.WINDOW_GUI_NORMAL)
        cv2.resizeWindow(window_name, self.window_w, self.window_h)
        cv2.setMouseCallback(window_name, self.mouse_callback, param=(self.window_w, self.window_h))

        while 0 <= self.current_idx < len(self.image_files):
            img_path = self.image_files[self.current_idx]
            cat_num = img_path.stem

            logger.info(f"Loading [{self.current_idx + 1}/{len(self.image_files)}] {img_path.name}...")
            self.current_img = cv2.imread(str(img_path))
            if self.current_img is None:
                logger.error(f"Failed to read image: {img_path}")
                self.current_idx += 1
                continue

            self.orig_h, self.orig_w = self.current_img.shape[:2]

            img_rgb = cv2.cvtColor(self.current_img, cv2.COLOR_BGR2RGB)
            self.predictor.set_image(img_rgb)

            # Automatically load any previously saved annotations for this voucher
            self.load_existing_annotations(cat_num)

            self.prompt_points = []
            self.prompt_labels = []
            self.prompt_box = None
            self.current_candidate_mask = None
            self.cut_lines = []
            self.lasso_points = []
            self.is_drawing_lasso = False
            self.zoom_level = 1.0
            self.pan_offset = [0, 0]
            self.exclusion_mode = False
            self.lasso_mode = False

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
                        self.lasso_points = []
                        self.is_drawing_lasso = False

                # 'b' Key -> GO BACK TO PREVIOUS VOUCHER
                elif key == ord('b'):
                    if self.current_idx > 0:
                        self.current_idx -= 1
                        logger.info("Navigating back to previous voucher...")
                        break
                    else:
                        logger.info("Already at the first voucher in queue.")

                # 'n' Key -> SKIP FORWARD TO NEXT VOUCHER (Without saving)
                elif key == ord('n'):
                    self.current_idx += 1
                    logger.info("Skipping forward to next voucher...")
                    break

                # 's' Key / Enter / 'v' -> SAVE & ADVANCE TO NEXT VOUCHER
                elif key in (13, ord('s'), ord('v'), ord('S')):
                    self.save_current_sheet()
                    self.current_idx += 1
                    break

                # Exclusion Mode Toggle ('e')
                elif key == ord('e'):
                    self.exclusion_mode = not self.exclusion_mode
                    if self.exclusion_mode:
                        self.lasso_mode = False
                    state_str = "EXCLUDE (Red, label 0)" if self.exclusion_mode else "INCLUDE (Green, label 1)"
                    logger.info(f"Exclusion mode toggled -> {state_str}")

                # Freehand Lasso Mode Toggle ('l')
                elif key == ord('l'):
                    self.lasso_mode = not self.lasso_mode
                    if self.lasso_mode:
                        self.exclusion_mode = False
                    state_str = "LASSO (Cyan)" if self.lasso_mode else "POINT PROMPT"
                    logger.info(f"Lasso mode toggled -> {state_str}")

                # Zoom Controls ('z' = Zoom In, 'x' = Zoom Out)
                elif key == ord('z'):
                    cx = self.last_mouse_pos[0] if self.last_mouse_pos else None
                    cy = self.last_mouse_pos[1] if self.last_mouse_pos else None
                    self.zoom(1.3, center_x=cx, center_y=cy)
                elif key == ord('x'):
                    cx = self.last_mouse_pos[0] if self.last_mouse_pos else None
                    cy = self.last_mouse_pos[1] if self.last_mouse_pos else None
                    self.zoom(1.0 / 1.3, center_x=cx, center_y=cy)

                # Pan Controls (WASD & Arrow keys)
                elif key in (ord('w'), 82, 0):  # Up
                    self.pan_offset[1] -= int(200 / self.zoom_level)
                    self.clamp_pan()
                elif key in (ord('s'), 84, 1):  # Down
                    self.pan_offset[1] += int(200 / self.zoom_level)
                    self.clamp_pan()
                elif key in (ord('a'), 81, 2):  # Left
                    self.pan_offset[0] -= int(200 / self.zoom_level)
                    self.clamp_pan()
                elif key in (ord('d'), 83, 3):  # Right
                    self.pan_offset[0] += int(200 / self.zoom_level)
                    self.clamp_pan()

                # Clear Active Candidate Prompts ('c')
                elif key == ord('c'):
                    self.prompt_points = []
                    self.prompt_labels = []
                    self.prompt_box = None
                    self.current_candidate_mask = None
                    self.cut_lines = []
                    self.lasso_points = []
                    self.is_drawing_lasso = False
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
