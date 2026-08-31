"""
===============================================================================
Module: sam2_rendering.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Visual rendering, HUD overlays, viewport coordinate transformations, and
    multi-instance botanical mask composition for interactive SAM 2 annotation.
===============================================================================
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np

CLASS_NAMES: List[str] = [
    "basal_leaf_whole",     # 0: Intact basal leaf (blade + petiole)
    "basal_leaf_partial",   # 1: Incomplete, overlapping, or partial basal leaf in rosette
    "cauline_leaf",         # 2: Sessile or lyrately-pinnatifid stem leaf
    "cauline_stem",         # 3: Main vertical flowering stalk / scape
    "root_rhizome",         # 4: Fibrous subterranean roots, rhizomes, and caudex
    "basal_rosette_clump",  # 5: Dense overlapping basal rosette center / crown
    "capitulum"             # 6: Inflorescence head / involucre / phyllaries
]

CLASS_COLORS: Dict[int, Tuple[int, int, int]] = {
    0: (0, 220, 0),       # basal_leaf_whole: Vibrant Green
    1: (100, 255, 100),   # basal_leaf_partial: Mint Green
    2: (0, 200, 255),     # cauline_leaf: Bright Cyan / Yellow-Green
    3: (0, 140, 255),     # cauline_stem: Orange
    4: (50, 50, 200),     # root_rhizome: Red/Brown
    5: (0, 100, 50),      # basal_rosette_clump: Dark Forest Green
    6: (0, 230, 255)      # capitulum: Yellow
}


def render_hud_overlay(
    display_img: np.ndarray,
    voucher_name: str,
    voucher_idx: int,
    total_vouchers: int,
    saved_instances: List[Dict[str, Any]],
    mode: str,
    zoom_level: float,
    pan_offset: Tuple[int, int],
) -> np.ndarray:
    """
    Renders top status bar HUD with voucher progress, mode, and shortcut reminders.
    """
    canvas = display_img.copy()
    h, w = canvas.shape[:2]

    # Top HUD banner background
    hud_h = 70
    overlay = canvas.copy()
    cv2.rectangle(overlay, (0, 0), (w, hud_h), (25, 25, 25), -1)
    cv2.addWeighted(overlay, 0.75, canvas, 0.25, 0, canvas)

    # Header text: Voucher progress
    prog_text = f"[{voucher_idx + 1}/{total_vouchers}] Voucher: {voucher_name} | Zoom: {zoom_level:.1f}x | Pan: ({pan_offset[0]}, {pan_offset[1]})"
    cv2.putText(canvas, prog_text, (15, 25), cv2.FONT_HERSHEY_SIMPLEX, 0.65, (240, 240, 240), 2, cv2.LINE_AA)

    # Active Mode display with dynamic color
    mode_color = (0, 255, 0) if mode == "INCLUDE" else ((0, 0, 255) if mode == "EXCLUDE" else (255, 255, 0))
    mode_str = f"MODE: [{mode}]"
    cv2.putText(canvas, mode_str, (15, 55), cv2.FONT_HERSHEY_SIMPLEX, 0.6, mode_color, 2, cv2.LINE_AA)

    # Instance count summary
    inst_summary = f"Saved: {len(saved_instances)} instances | Keys 0-6: Classify | Enter: Save | n: Skip | b: Back | q: Quit"
    cv2.putText(canvas, inst_summary, (200, 55), cv2.FONT_HERSHEY_SIMPLEX, 0.52, (200, 200, 200), 1, cv2.LINE_AA)

    return canvas


def apply_viewport_transform(
    image: np.ndarray,
    zoom_level: float,
    pan_offset: Tuple[int, int],
    target_w: int,
    target_h: int
) -> Tuple[np.ndarray, Tuple[float, float, int, int]]:
    """
    Crops and rescales the image to match viewport pan and zoom levels.

    Returns:
        Tuple: (rendered_viewport_img, (scale_x, scale_y, crop_x0, crop_y0))
    """
    img_h, img_w = image.shape[:2]
    crop_w = int(img_w / zoom_level)
    crop_h = int(img_h / zoom_level)

    crop_x0 = max(0, min(img_w - crop_w, pan_offset[0]))
    crop_y0 = max(0, min(img_h - crop_h, pan_offset[1]))

    cropped = image[crop_y0:crop_y0 + crop_h, crop_x0:crop_x0 + crop_w]
    scaled = cv2.resize(cropped, (target_w, target_h), interpolation=cv2.INTER_LINEAR)

    scale_x = target_w / crop_w
    scale_y = target_h / crop_h

    return scaled, (scale_x, scale_y, crop_x0, crop_y0)


def compose_mask_overlay(
    base_image: np.ndarray,
    saved_instances: List[Dict[str, Any]],
    candidate_mask: Optional[np.ndarray] = None,
    alpha: float = 0.45
) -> np.ndarray:
    """
    Composites saved multi-class instances and candidate mask overlays onto image.
    """
    overlay = base_image.copy()

    for inst in saved_instances:
        mask = inst["mask"]
        class_id = inst["class_id"]
        color = CLASS_COLORS.get(class_id, (0, 255, 0))

        colored_mask = np.zeros_like(base_image)
        colored_mask[mask > 0] = color
        overlay[mask > 0] = cv2.addWeighted(base_image[mask > 0], 1.0 - alpha, colored_mask[mask > 0], alpha, 0)

        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(overlay, contours, -1, color, 2)

    if candidate_mask is not None and np.count_nonzero(candidate_mask) > 0:
        cyan = (255, 255, 0)
        c_overlay = np.zeros_like(base_image)
        c_overlay[candidate_mask > 0] = cyan
        overlay[candidate_mask > 0] = cv2.addWeighted(
            base_image[candidate_mask > 0], 0.4, c_overlay[candidate_mask > 0], 0.6, 0
        )
        c_contours, _ = cv2.findContours(candidate_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(overlay, c_contours, -1, (255, 255, 255), 2)

    return overlay
