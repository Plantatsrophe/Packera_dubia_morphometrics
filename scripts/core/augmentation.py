
import os
import sys
import shutil
import cv2
import math
import numpy as np
import random
import yaml
import glob
import logging
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any, Union
from pathlib import Path

# Common imports for the project
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE, DEFAULT_RAW_DIR, DEFAULT_CURATED_CSV, DEFAULT_OUTPUT_DIR, DEFAULT_CONFIG_PATH, DEFAULT_QC_DIR
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):
        return iterable


from scripts.core.artifact_harvester import ArtifactHarvester

def extract_hard_negative_background_crop(
    image_bgr: np.ndarray,
    all_annotations: List[InstanceAnnotation],
    crop_size: Tuple[int, int] = (1024, 1024),
    max_attempts: int = 25
) -> Optional[np.ndarray]:
    """
    Extracts a pure background sheet crop containing zero plant or artifact instances.
    Used to generate negative training samples with empty annotation files.
    
    Args:
        image_bgr: Full herbarium sheet image.
        all_annotations: Existing annotations to avoid.
        crop_size: Output (width, height) of the negative tile.
        max_attempts: Number of random spatial samples before fallback.
        
    Returns:
        Pure background BGR crop, or None if no vacant region satisfies constraints.
    """
    h, w = image_bgr.shape[:2]
    cw, ch = crop_size
    if w <= cw or h <= ch:
        return None

    # Build occupied occupancy mask
    occupancy = np.zeros((h, w), dtype=np.uint8)
    for ann in all_annotations:
        x1, y1, x2, y2 = [int(v) for v in ann.bbox]
        # Pad bounding box to guarantee complete vacancy
        px1 = max(0, x1 - 30)
        py1 = max(0, y1 - 30)
        px2 = min(w, x2 + 30)
        py2 = min(h, y2 + 30)
        cv2.rectangle(occupancy, (px1, py1), (px2, py2), 255, -1)

    # Search for an unallocated window
    for _ in range(max_attempts):
        rx = random.randint(0, w - cw)
        ry = random.randint(0, h - ch)
        sub_occupancy = occupancy[ry:ry+ch, rx:rx+cw]
        if np.count_nonzero(sub_occupancy) == 0:
            return image_bgr[ry:ry+ch, rx:rx+cw].copy()

    # Fallback: Sample sheet corner/margin with minimal gradient
    margin_w = min(cw, int(w * 0.25))
    margin_h = min(ch, int(h * 0.25))
    corner_crop = image_bgr[:margin_h, :margin_w].copy()
    return cv2.resize(corner_crop, (cw, ch), interpolation=cv2.INTER_LINEAR)


class SyntheticOcclusionAugmenter:
    """
    Applies synthetic copy-paste augmentations of non-plant artifacts
    adjacent to, touching, or occluding annotated botanical leaves.
    Updates polygon geometries and bounding boxes to reflect sharp boundaries.
    """
    def __init__(self, harvester: ArtifactHarvester, rng_seed: int = 42):
        self.harvester = harvester
        self.rng = random.Random(rng_seed)
        self.np_rng = np.random.default_rng(rng_seed)

    def apply_copy_paste_augmentation(
        self,
        image_bgr: np.ndarray,
        annotations: List[InstanceAnnotation],
        paste_probability: float = 0.75,
        max_pastes_per_image: int = 3
    ) -> Tuple[np.ndarray, List[InstanceAnnotation]]:
        """
        Executes dynamic copy-paste augmentation on a single herbarium sheet.
        
        Args:
            image_bgr: Original herbarium sheet BGR image.
            annotations: Current list of InstanceAnnotations.
            paste_probability: Probability of applying augmentation.
            max_pastes_per_image: Max artifact patches to paste.
            
        Returns:
            Tuple of (augmented_image_bgr, updated_annotations_list).
        """
        if self.rng.random() > paste_probability or not annotations:
            return image_bgr, annotations

        aug_image = image_bgr.copy()
        img_h, img_w = aug_image.shape[:2]
        updated_anns = [ann for ann in annotations]

        # Target basal_leaf instances for occlusion / adjacency
        leaf_indices = [
            i for i, ann in enumerate(updated_anns)
            if ann.class_id == CLASS_MAP["basal_leaf"]
        ]

        if not leaf_indices:
            return aug_image, updated_anns

        num_pastes = self.rng.randint(1, max_pastes_per_image)

        for _ in range(num_pastes):
            target_idx = self.rng.choice(leaf_indices)
            target_leaf = updated_anns[target_idx]
            lx1, ly1, lx2, ly2 = target_leaf.bbox
            lw = max(10, int(lx2 - lx1))
            lh = max(10, int(ly2 - ly1))

            # Sample random artifact patch (tape, label, ruler, color swatch)
            artifact_sample = self.harvester.get_random_artifact_crop()
            if artifact_sample is None:
                continue

            art_class_name, art_crop = artifact_sample
            art_class_id = CLASS_MAP[art_class_name]
            ah, aw = art_crop.shape[:2]

            # Scale artifact appropriately relative to the leaf (20% to 80% leaf scale)
            scale = self.rng.uniform(0.3, 0.9) * (max(lw, lh) / max(aw, ah, 1))
            new_aw = max(15, min(int(aw * scale), img_w // 3))
            new_ah = max(10, min(int(ah * scale), img_h // 3))
            resized_art = cv2.resize(art_crop, (new_aw, new_ah), interpolation=cv2.INTER_AREA)

            # Random rotation (-35 to +35 degrees)
            angle = self.rng.uniform(-35, 35)
            rot_mat = cv2.getRotationMatrix2D((new_aw / 2, new_ah / 2), angle, 1.0)
            cos = np.abs(rot_mat[0, 0])
            sin = np.abs(rot_mat[0, 1])
            bound_w = int((new_ah * sin) + (new_aw * cos))
            bound_h = int((new_ah * cos) + (new_aw * sin))
            rot_mat[0, 2] += (bound_w / 2) - (new_aw / 2)
            rot_mat[1, 2] += (bound_h / 2) - (new_ah / 2)

            rotated_art = cv2.warpAffine(
                resized_art, rot_mat, (bound_w, bound_h),
                borderMode=cv2.BORDER_CONSTANT, borderValue=(255, 255, 255)
            )

            # Create binary foreground mask for the rotated artifact
            art_gray = cv2.cvtColor(rotated_art, cv2.COLOR_BGR2GRAY)
            # Anything not pure white border is artifact foreground
            art_mask = (art_gray < 250).astype(np.uint8) * 255
            # Feather edge for seamless paper texture alpha-blending
            art_mask_blurred = cv2.GaussianBlur(art_mask, (3, 3), 0)
            alpha = (art_mask_blurred.astype(np.float32) / 255.0)[:, :, np.newaxis]

            # Choose placement mode:
            # Mode A: Partial Occlusion (Pasted across leaf blade/petiole)
            # Mode B: Adjacent / Touching (Pasted directly adjacent to margin)
            mode = self.rng.choice(["occlusion", "adjacent"])

            if mode == "occlusion":
                paste_x = int(self.rng.uniform(lx1 - bound_w * 0.3, lx2 - bound_w * 0.7))
                paste_y = int(self.rng.uniform(ly1 - bound_h * 0.3, ly2 - bound_h * 0.7))
            else:
                offset_side = self.rng.choice(["left", "right", "top", "bottom"])
                if offset_side == "left":
                    paste_x = int(lx1 - bound_w * 0.9)
                    paste_y = int(ly1 + self.rng.uniform(-bound_h * 0.2, lh * 0.5))
                elif offset_side == "right":
                    paste_x = int(lx2 - bound_w * 0.1)
                    paste_y = int(ly1 + self.rng.uniform(-bound_h * 0.2, lh * 0.5))
                elif offset_side == "top":
                    paste_x = int(lx1 + self.rng.uniform(-bound_w * 0.2, lw * 0.5))
                    paste_y = int(ly1 - bound_h * 0.9)
                else:
                    paste_x = int(lx1 + self.rng.uniform(-bound_w * 0.2, lw * 0.5))
                    paste_y = int(ly2 - bound_h * 0.1)

            # Clip placement coordinates to image boundaries
            paste_x1 = max(0, paste_x)
            paste_y1 = max(0, paste_y)
            paste_x2 = min(img_w, paste_x + bound_w)
            paste_y2 = min(img_h, paste_y + bound_h)

            crop_w = paste_x2 - paste_x1
            crop_h = paste_y2 - paste_y1

            if crop_w <= 5 or crop_h <= 5:
                continue

            art_sub_x1 = paste_x1 - paste_x
            art_sub_y1 = paste_y1 - paste_y
            art_sub_x2 = art_sub_x1 + crop_w
            art_sub_y2 = art_sub_y1 + crop_h

            sub_art = rotated_art[art_sub_y1:art_sub_y2, art_sub_x1:art_sub_x2]
            sub_alpha = alpha[art_sub_y1:art_sub_y2, art_sub_x1:art_sub_x2]

            # Alpha-blend artifact into the canvas
            target_roi = aug_image[paste_y1:paste_y2, paste_x1:paste_x2].astype(np.float32)
            blended_roi = (sub_art.astype(np.float32) * sub_alpha) + (target_roi * (1.0 - sub_alpha))
            aug_image[paste_y1:paste_y2, paste_x1:paste_x2] = np.clip(blended_roi, 0, 255).astype(np.uint8)

            # Define new artifact polygon and bounding box
            art_poly = np.array([
                [paste_x1, paste_y1], [paste_x2, paste_y1],
                [paste_x2, paste_y2], [paste_x1, paste_y2]
            ], dtype=np.float32)

            new_art_ann = InstanceAnnotation(
                class_id=art_class_id,
                polygon=art_poly,
                bbox=(paste_x1, paste_y1, paste_x2, paste_y2),
                confidence=1.0,
                is_synthetic=True,
                tag=f"aug_copy_paste_{art_class_name}"
            )
            updated_anns.append(new_art_ann)

            # Update occluded leaf polygon geometry dynamically (Boolean Difference)
            if mode == "occlusion" and len(target_leaf.polygon) >= 3:
                updated_leaf_poly = self._compute_occluded_polygon(
                    target_leaf.polygon, (paste_x1, paste_y1, paste_x2, paste_y2), img_w, img_h
                )
                if updated_leaf_poly is not None and len(updated_leaf_poly) >= 3:
                    target_leaf.polygon = updated_leaf_poly
                    # Recompute bounding box
                    x_min = float(np.min(updated_leaf_poly[:, 0]))
                    y_min = float(np.min(updated_leaf_poly[:, 1]))
                    x_max = float(np.max(updated_leaf_poly[:, 0]))
                    y_max = float(np.max(updated_leaf_poly[:, 1]))
                    target_leaf.bbox = (x_min, y_min, x_max, y_max)

        return aug_image, updated_anns

    def _compute_occluded_polygon(
        self,
        leaf_poly: np.ndarray,
        occluder_bbox: Tuple[int, int, int, int],
        img_w: int,
        img_h: int
    ) -> Optional[np.ndarray]:
        """
        Subtracts the occluding rectangular bounding mask from the leaf polygon.
        """
        ox1, oy1, ox2, oy2 = occluder_bbox
        # Create full-sheet binary mask for the leaf
        leaf_mask = np.zeros((img_h, img_w), dtype=np.uint8)
        pts = leaf_poly.astype(np.int32)
        cv2.fillPoly(leaf_mask, [pts], 255)

        # Subtract occluder region
        leaf_mask[oy1:oy2, ox1:ox2] = 0

        # Extract largest remaining contour
        contours, _ = cv2.findContours(leaf_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            return None

        largest_cnt = max(contours, key=cv2.contourArea)
        if cv2.contourArea(largest_cnt) < 50:
            return None

        epsilon = 0.005 * cv2.arcLength(largest_cnt, True)
        approx = cv2.approxPolyDP(largest_cnt, epsilon, True).reshape(-1, 2)
        return approx.astype(np.float32)


