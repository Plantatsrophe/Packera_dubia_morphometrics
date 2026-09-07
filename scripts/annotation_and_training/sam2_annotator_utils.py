"""
===============================================================================
Module: sam2_annotator_utils.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Consolidated helper library for SAM 2 interactive botanical annotation,
    including:
      1. Geometric and sub-pixel contour processing (knife slicing, lasso rasterization,
         polygon-to-box projection, interior pole calculation, sub-pixel polygon extraction).
      2. Native X11 HUD and viewport rendering (viewport coordinate transforms,
         cached base layers, live candidate mask composition).
      3. Standardized COCO 1.0 dataset export aligned with LeafMachine2 Plant Component
         Detector (PCD) taxonomy.
===============================================================================
"""

from __future__ import annotations

import json
import logging
import random
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import cv2
import numpy as np

logger = logging.getLogger("SAM2AnnotatorUtils")

# ===============================================================================
# 1. Botanical Taxonomic Classes & Visual Palette
# ===============================================================================

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

# Strict class mapping from SAM 2 botanical annotations to LeafMachine2 PCD taxonomy
PCD_CLASS_MAPPING: Dict[str, str] = {
    "basal_leaf_whole": "ideal_leaf",
    "basal_leaf_partial": "partial_leaf",
}

# Categories intentionally omitted so Detectron2 treats them as background noise
OMITTED_CLASSES: Set[str] = {
    "cauline_leaf",
    "cauline_stem",
    "root_rhizome",
    "basal_rosette_clump",
    "capitulum",
}

# Standardized COCO 1.0 categories for LeafMachine2 Plant Component Detector
PCD_COCO_CATEGORIES: List[Dict[str, Any]] = [
    {
        "id": 1,
        "name": "ideal_leaf",
        "supercategory": "plant_component",
    },
    {
        "id": 2,
        "name": "partial_leaf",
        "supercategory": "plant_component",
    },
]

CATEGORY_NAME_TO_ID: Dict[str, int] = {
    cat["name"]: cat["id"] for cat in PCD_COCO_CATEGORIES
}


# ===============================================================================
# 2. Geometric & Sub-Pixel Boundary Routines
# ===============================================================================

def clip_box_to_image(
    box: Tuple[int, int, int, int],
    img_w: int,
    img_h: int
) -> Tuple[int, int, int, int]:
    """Clips bounding box coordinates (x0, y0, x1, y1) to valid image boundaries."""
    x0, y0, x1, y1 = box
    x_min = max(0, min(x0, x1))
    y_min = max(0, min(y0, y1))
    x_max = min(img_w, max(x0, x1))
    y_max = min(img_h, max(y0, y1))
    return x_min, y_min, x_max, y_max


def split_mask_with_knife_line(
    binary_mask: np.ndarray,
    line_start: Tuple[int, int],
    line_end: Tuple[int, int],
    line_thickness: int = 3
) -> np.ndarray:
    """
    Sever a binary mask using a knife cut line (e.g. to isolate blade from petiole/roots).

    Args:
        binary_mask: 2D uint8 binary mask array.
        line_start: (x, y) start coordinate.
        line_end: (x, y) end coordinate.
        line_thickness: Cut line stroke width in pixels.

    Returns:
        np.ndarray: Mask with knife cut line zeroed out.
    """
    cut_mask = binary_mask.copy()
    cv2.line(cut_mask, line_start, line_end, 0, thickness=line_thickness)
    return cut_mask


def rasterize_lasso_polygon(
    lasso_points: List[Tuple[int, int]],
    img_h: int,
    img_w: int
) -> np.ndarray:
    """
    Rasterizes a sequence of freehand lasso boundary points into a filled binary mask.

    Args:
        lasso_points: List of (x, y) polygon vertices.
        img_h: Image height in pixels.
        img_w: Image width in pixels.

    Returns:
        np.ndarray: 2D uint8 binary mask (0 background, 255 filled polygon).
    """
    mask = np.zeros((img_h, img_w), dtype=np.uint8)
    if len(lasso_points) < 3:
        return mask
    pts = np.array(lasso_points, dtype=np.int32).reshape((-1, 1, 2))
    cv2.fillPoly(mask, [pts], 255)
    return mask


def polygon_to_bounding_box(
    polygon_points: List[Tuple[int, int]]
) -> Optional[Tuple[int, int, int, int]]:
    """
    Computes enclosing bounding box (min_x, min_y, max_x, max_y) from polygon vertices.
    """
    if not polygon_points:
        return None
    xs = [p[0] for p in polygon_points]
    ys = [p[1] for p in polygon_points]
    return min(xs), min(ys), max(xs), max(ys)


def polygon_interior_point(
    polygon_points: List[Tuple[int, int]],
    img_h: int,
    img_w: int
) -> Optional[Tuple[float, float]]:
    """
    Finds optimal interior point (pole of inaccessibility) inside an arbitrary polygon
    using Euclidean distance transform to guide SAM 2 point prompting.
    """
    if len(polygon_points) < 3:
        return None
    mask = rasterize_lasso_polygon(polygon_points, img_h, img_w)
    if np.count_nonzero(mask) == 0:
        return None
    dist = cv2.distanceTransform(mask, cv2.DIST_L2, 5)
    _, max_val, _, max_loc = cv2.minMaxLoc(dist)
    if max_val > 0:
        return float(max_loc[0]), float(max_loc[1])
    return float(polygon_points[0][0]), float(polygon_points[0][1])


def mask_to_yolo_bbox(
    binary_mask: np.ndarray,
    class_id: int
) -> Optional[Tuple[int, float, float, float, float]]:
    """
    Computes normalized YOLO format bounding box (class_id, x_center, y_center, width, height).
    """
    y_indices, x_indices = np.where(binary_mask > 0)
    if len(x_indices) == 0 or len(y_indices) == 0:
        return None

    h, w = binary_mask.shape[:2]
    x_min, x_max = float(np.min(x_indices)), float(np.max(x_indices))
    y_min, y_max = float(np.min(y_indices)), float(np.max(y_indices))

    box_w = (x_max - x_min + 1.0) / w
    box_h = (y_max - y_min + 1.0) / h
    x_center = (x_min + (x_max - x_min) / 2.0) / w
    y_center = (y_min + (y_max - y_min) / 2.0) / h

    return (
        class_id,
        round(x_center, 6),
        round(y_center, 6),
        round(box_w, 6),
        round(box_h, 6),
    )


def mask_to_normalized_polygon(
    binary_mask: np.ndarray,
    class_id: int,
    approx_epsilon: float = 1.0
) -> Optional[str]:
    """Converts binary mask to normalized YOLO segmentation string: 'class_id x1 y1 x2 y2 ...'"""
    contours, _ = cv2.findContours(binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None

    largest_cnt = max(contours, key=cv2.contourArea)
    if approx_epsilon > 0:
        largest_cnt = cv2.approxPolyDP(largest_cnt, approx_epsilon, True)

    if len(largest_cnt) < 3:
        return None

    h, w = binary_mask.shape[:2]
    pts = largest_cnt.reshape(-1, 2)
    norm_pts = []
    for x, y in pts:
        norm_pts.append(f"{round(x / w, 6):.6f}")
        norm_pts.append(f"{round(y / h, 6):.6f}")

    return f"{class_id} " + " ".join(norm_pts)


def mask_to_polygons(
    binary_mask: np.ndarray,
    min_area_px: float = 50.0,
    approx_epsilon: float = 1.0,
) -> Tuple[List[List[float]], float, List[float]]:
    """
    Extracts vectorized polygon contours with sub-pixel coordinate precision
    and bounding boxes from a 2D binary mask.

    Args:
        binary_mask: 2D uint8 binary array (0=background, >0=foreground).
        min_area_px: Minimum contour area threshold in pixels.
        approx_epsilon: Polygon approximation epsilon for contour point decimation.

    Returns:
        Tuple: (segmentation_polygons, total_area, [x_min, y_min, width, height])
    """
    if binary_mask is None or binary_mask.size == 0 or np.count_nonzero(binary_mask) == 0:
        return [], 0.0, [0.0, 0.0, 0.0, 0.0]

    contours, _ = cv2.findContours(
        binary_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_TC89_KCOS
    )

    polygons: List[List[float]] = []
    total_area = 0.0
    all_x: List[float] = []
    all_y: List[float] = []

    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area < min_area_px:
            continue

        if approx_epsilon > 0:
            cnt = cv2.approxPolyDP(cnt, approx_epsilon, True)

        if len(cnt) < 3:
            continue

        flattened = cnt.flatten().astype(float).tolist()
        if len(flattened) >= 6:
            rounded = [round(val, 2) for val in flattened]
            polygons.append(rounded)
            total_area += area
            all_x.extend(rounded[0::2])
            all_y.extend(rounded[1::2])

    if not polygons or not all_x:
        return [], 0.0, [0.0, 0.0, 0.0, 0.0]

    x_min = min(all_x)
    y_min = min(all_y)
    x_max = max(all_x)
    y_max = max(all_y)
    bbox = [round(x_min, 2), round(y_min, 2), round(x_max - x_min, 2), round(y_max - y_min, 2)]

    return polygons, round(total_area, 2), bbox


# ===============================================================================
# 3. Viewport & Native Rendering Routines
# ===============================================================================

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
    if mode == "POLYGON":
        mode_color = (0, 220, 255)
        mode_str = "MODE: [POLYGON BOX]"
        inst_summary = "Left-Click: Mark Points | Click Start / Enter / Right-Click: Close Box | Backspace: Undo | P: Exit"
    elif mode == "KNIFE":
        mode_color = (0, 100, 255)
        mode_str = "MODE: [KNIFE CUT]"
        inst_summary = "Drag line across mask to sever overlapping parts | K: Exit Knife"
    else:
        mode_color = (0, 255, 0)
        mode_str = "MODE: [SELECT]"
        inst_summary = f"Saved: {len(saved_instances)} | P: Polygon Box | K: Knife | C: Clear | U: Undo | 0-6: Classify | Enter: Save"

    cv2.putText(canvas, mode_str, (15, 55), cv2.FONT_HERSHEY_SIMPLEX, 0.6, mode_color, 2, cv2.LINE_AA)
    cv2.putText(canvas, inst_summary, (220, 55), cv2.FONT_HERSHEY_SIMPLEX, 0.48, (220, 220, 220), 1, cv2.LINE_AA)

    return canvas


def apply_viewport_transform(
    image: np.ndarray,
    zoom_level: float,
    pan_offset: Tuple[int, int],
    target_w: int,
    target_h: int,
    bg_color: Tuple[int, int, int] = (35, 35, 35)
) -> Tuple[np.ndarray, Tuple[float, float, int, int]]:
    """
    Crops and rescales image to match viewport pan and zoom levels,
    allowing smooth panning beyond image boundaries with neutral background padding.

    Returns:
        Tuple: (rendered_viewport_img, (scale_x, scale_y, crop_x0, crop_y0))
    """
    img_h, img_w = image.shape[:2]
    crop_w = max(10, int(img_w / max(zoom_level, 1.0)))
    crop_h = max(10, int(img_h / max(zoom_level, 1.0)))

    crop_x0 = pan_offset[0]
    crop_y0 = pan_offset[1]

    # Calculate valid source image overlap bounds
    src_x0 = max(0, crop_x0)
    src_x1 = min(img_w, crop_x0 + crop_w)
    src_y0 = max(0, crop_y0)
    src_y1 = min(img_h, crop_y0 + crop_h)

    # Destination canvas coordinates
    dst_x0 = max(0, -crop_x0) if crop_x0 < 0 else 0
    dst_y0 = max(0, -crop_y0) if crop_y0 < 0 else 0
    dst_x1 = dst_x0 + max(0, src_x1 - src_x0)
    dst_y1 = dst_y0 + max(0, src_y1 - src_y0)

    canvas = np.full((crop_h, crop_w, 3), bg_color, dtype=image.dtype)
    if src_x1 > src_x0 and src_y1 > src_y0 and dst_x1 <= crop_w and dst_y1 <= crop_h:
        canvas[dst_y0:dst_y1, dst_x0:dst_x1] = image[src_y0:src_y1, src_x0:src_x1]

    scaled = cv2.resize(canvas, (target_w, target_h), interpolation=cv2.INTER_LINEAR)
    scale_x = target_w / max(crop_w, 1)
    scale_y = target_h / max(crop_h, 1)

    return scaled, (scale_x, scale_y, crop_x0, crop_y0)


def compose_mask_overlay(
    base_image: np.ndarray,
    saved_instances: List[Dict[str, Any]],
    candidate_mask: Optional[np.ndarray] = None,
    alpha: float = 0.45
) -> np.ndarray:
    """
    Composites saved multi-class instances and candidate mask overlays onto image.
    Used for caching the pre-composited base layer.
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


def overlay_candidate_mask_on_viewport(
    viewport_bgr: np.ndarray,
    candidate_mask: Optional[np.ndarray],
    transform: Tuple[float, float, int, int],
    alpha: float = 0.55
) -> np.ndarray:
    """
    Overlays active SAM 2 candidate segmentation mask directly onto cropped
    viewport frame in viewport space for sub-millisecond rendering speed.
    """
    if candidate_mask is None or np.count_nonzero(candidate_mask) == 0:
        return viewport_bgr

    target_h, target_w = viewport_bgr.shape[:2]
    scale_x, scale_y, crop_x0, crop_y0 = transform
    img_h, img_w = candidate_mask.shape[:2]

    crop_w = max(10, int(target_w / max(scale_x, 1e-6)))
    crop_h = max(10, int(target_h / max(scale_y, 1e-6)))

    src_x0 = max(0, crop_x0)
    src_x1 = min(img_w, crop_x0 + crop_w)
    src_y0 = max(0, crop_y0)
    src_y1 = min(img_h, crop_y0 + crop_h)

    if src_x1 <= src_x0 or src_y1 <= src_y0:
        return viewport_bgr

    mask_slice = candidate_mask[src_y0:src_y1, src_x0:src_x1]
    if not np.any(mask_slice > 0):
        return viewport_bgr

    dst_x0 = max(0, -crop_x0) if crop_x0 < 0 else 0
    dst_y0 = max(0, -crop_y0) if crop_y0 < 0 else 0
    dst_x1 = dst_x0 + (src_x1 - src_x0)
    dst_y1 = dst_y0 + (src_y1 - src_y0)

    canvas_mask = np.zeros((crop_h, crop_w), dtype=np.uint8)
    canvas_mask[dst_y0:dst_y1, dst_x0:dst_x1] = mask_slice

    scaled_mask = cv2.resize(canvas_mask, (target_w, target_h), interpolation=cv2.INTER_NEAREST)
    locs = scaled_mask > 0
    if np.any(locs):
        viewport_bgr[locs] = (
            viewport_bgr[locs] * (1.0 - alpha) + np.array([255, 255, 0], dtype=np.uint8) * alpha
        ).astype(np.uint8)
        cnts, _ = cv2.findContours(scaled_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        cv2.drawContours(viewport_bgr, cnts, -1, (255, 255, 255), 2)

    return viewport_bgr


# ===============================================================================
# 4. Standardized COCO 1.0 Dataset Exporter
# ===============================================================================

def parse_mask_filename(mask_filename: str) -> Tuple[str, str, int]:
    """
    Parses catalogNumber / voucher ID, class label, and instance ID from a mask filename.

    Supports multiple established naming patterns:
      1. {voucher_id}_inst{id:02d}_{class_label}.png (e.g. '000331814_inst00_basal_leaf_partial.png')
      2. {voucher_id}_{class_label}_instance_{id}.png (e.g. 'NCU00001234_basal_leaf_whole_instance_1.png')
      3. {voucher_id}_{class_label}_{id}.png
    """
    stem = Path(mask_filename).stem

    # Pattern 1: {voucher_id}_inst{id}_{label}
    m1 = re.match(r"^(.*?)_inst(\d+)_(.+)$", stem)
    if m1:
        voucher_id = m1.group(1)
        inst_id = int(m1.group(2))
        label = m1.group(3)
        return voucher_id, label, inst_id

    # Pattern 2: {voucher_id}_{label}_instance_{id}
    m2 = re.match(r"^(.*?)(?:_instance_(\d+))?$", stem)
    if m2:
        core = m2.group(1)
        inst_id = int(m2.group(2)) if m2.group(2) else 1
    else:
        core = stem
        inst_id = 1

    all_classes = list(PCD_CLASS_MAPPING.keys()) + list(OMITTED_CLASSES)
    for known_class in all_classes:
        if f"_{known_class}" in core:
            voucher_id = core.replace(f"_{known_class}", "")
            return voucher_id, known_class, inst_id

    return core, "unknown", inst_id


def convert_masks_to_coco_dataset(
    masks_dir: Union[str, Path],
    images_dir: Optional[Union[str, Path]] = None,
    image_dim_map: Optional[Dict[str, Tuple[int, int]]] = None,
    default_width: int = 4000,
    default_height: int = 6000,
    min_area_px: float = 50.0,
    image_extension: str = ".jpg",
) -> Dict[str, Any]:
    """
    Scans a directory of binary PNG masks and compiles a COCO 1.0 JSON format dictionary.

    Args:
        masks_dir: Path to directory containing binary PNG masks.
        images_dir: Optional path to raw voucher images to extract exact image dimensions.
        image_dim_map: Dict mapping image stems to (width, height).
        default_width: Fallback width if dimensions cannot be inferred.
        default_height: Fallback height if dimensions cannot be inferred.
        min_area_px: Minimum polygon area threshold.
        image_extension: Target herbarium sheet extension (default .jpg).

    Returns:
        Dict[str, Any]: Standardized COCO 1.0 format dataset.
    """
    masks_dir = Path(masks_dir)
    if not masks_dir.exists():
        logger.warning(f"Masks directory not found: {masks_dir}")
        return {"images": [], "annotations": [], "categories": PCD_COCO_CATEGORIES}

    mask_files = sorted([f for f in masks_dir.iterdir() if f.suffix.lower() == ".png"])
    logger.info(f"Found {len(mask_files)} mask files in {masks_dir}")

    images_dir_path = Path(images_dir) if images_dir else None
    images_dict: Dict[str, Dict[str, Any]] = {}
    annotations_list: List[Dict[str, Any]] = []

    image_id_counter = 1
    annotation_id_counter = 1

    for mask_path in mask_files:
        cat_num, class_label, inst_id = parse_mask_filename(mask_path.name)

        if class_label in OMITTED_CLASSES:
            continue

        pcd_category = PCD_CLASS_MAPPING.get(class_label)
        if not pcd_category or pcd_category not in CATEGORY_NAME_TO_ID:
            continue

        category_id = CATEGORY_NAME_TO_ID[pcd_category]

        img_file_name = f"{cat_num}{image_extension}"
        if img_file_name not in images_dict:
            w, h = default_width, default_height

            # Resolve actual dimensions if image exists on disk or in dim_map
            if image_dim_map and cat_num in image_dim_map:
                w, h = image_dim_map[cat_num]
            elif images_dir_path:
                candidate_img = images_dir_path / img_file_name
                if candidate_img.exists():
                    try:
                        probe = cv2.imread(str(candidate_img))
                        if probe is not None:
                            h, w = probe.shape[:2]
                    except Exception:
                        pass

            images_dict[img_file_name] = {
                "id": image_id_counter,
                "file_name": img_file_name,
                "width": w,
                "height": h,
            }
            curr_image_id = image_id_counter
            image_id_counter += 1
        else:
            curr_image_id = images_dict[img_file_name]["id"]

        mask_img = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
        if mask_img is None:
            continue

        if images_dict[img_file_name]["width"] == default_width:
            h, w = mask_img.shape[:2]
            images_dict[img_file_name]["width"] = w
            images_dict[img_file_name]["height"] = h

        polygons, area, bbox = mask_to_polygons(mask_img, min_area_px=min_area_px)
        if not polygons or area < min_area_px:
            continue

        annotations_list.append({
            "id": annotation_id_counter,
            "image_id": curr_image_id,
            "category_id": category_id,
            "segmentation": polygons,
            "area": area,
            "bbox": bbox,
            "iscrowd": 0,
        })
        annotation_id_counter += 1

    coco_doc = {
        "info": {
            "description": "LeafMachine2 Plant Component Detector - Packera Dataset",
            "version": "1.0",
            "year": 2026,
            "contributor": "NCU Herbarium / UNC Chapel Hill",
        },
        "licenses": [],
        "images": list(images_dict.values()),
        "annotations": annotations_list,
        "categories": PCD_COCO_CATEGORIES,
    }

    logger.info(
        f"Generated COCO dataset: {len(coco_doc['images'])} images, "
        f"{len(coco_doc['annotations'])} annotations across {len(PCD_COCO_CATEGORIES)} categories."
    )
    return coco_doc


def split_coco_dataset(
    coco_data: Dict[str, Any],
    val_ratio: float = 0.20,
    seed: int = 42
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Partitions a COCO dataset into training and validation splits at the image level."""
    random.seed(seed)
    images = coco_data.get("images", [])
    annotations = coco_data.get("annotations", [])
    categories = coco_data.get("categories", PCD_COCO_CATEGORIES)
    info = coco_data.get("info", {})

    image_ids = [img["id"] for img in images]
    random.shuffle(image_ids)

    n_val = int(len(image_ids) * val_ratio)
    val_image_ids = set(image_ids[:n_val])
    train_image_ids = set(image_ids[n_val:])

    train_images = [img for img in images if img["id"] in train_image_ids]
    val_images = [img for img in images if img["id"] in val_image_ids]

    train_annotations = [ann for ann in annotations if ann["image_id"] in train_image_ids]
    val_annotations = [ann for ann in annotations if ann["image_id"] in val_image_ids]

    train_doc = {
        "info": {**info, "split": "train"},
        "licenses": coco_data.get("licenses", []),
        "images": train_images,
        "annotations": train_annotations,
        "categories": categories,
    }

    val_doc = {
        "info": {**info, "split": "val"},
        "licenses": coco_data.get("licenses", []),
        "images": val_images,
        "annotations": val_annotations,
        "categories": categories,
    }

    return train_doc, val_doc


def save_coco_json(coco_data: Dict[str, Any], output_path: Union[str, Path]) -> None:
    """Writes COCO dataset dictionary to a JSON file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(coco_data, f, indent=2)
    logger.info(f"Saved COCO dataset to {output_path}")


def export_coco_annotations(
    masks_dir: Union[str, Path],
    output_coco_path: Union[str, Path],
    images_dir: Optional[Union[str, Path]] = None,
    val_split: float = 0.0,
    min_area_px: float = 50.0,
) -> Dict[str, Any]:
    """
    High-level export helper: converts masks in `masks_dir` directly to COCO JSON format.
    """
    output_path = Path(output_coco_path)
    coco_data = convert_masks_to_coco_dataset(
        masks_dir=masks_dir,
        images_dir=images_dir,
        min_area_px=min_area_px,
    )

    if val_split > 0.0:
        train_doc, val_doc = split_coco_dataset(coco_data, val_ratio=val_split)
        stem = output_path.stem
        train_path = output_path.with_name(f"{stem}_train.json")
        val_path = output_path.with_name(f"{stem}_val.json")
        save_coco_json(train_doc, train_path)
        save_coco_json(val_doc, val_path)
        save_coco_json(coco_data, output_path)
    else:
        save_coco_json(coco_data, output_path)

    return coco_data
