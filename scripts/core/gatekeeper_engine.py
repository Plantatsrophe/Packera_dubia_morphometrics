#!/usr/bin/env python3
"""
===============================================================================
Module: artifact_filter_gatekeeper.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU) /
             Google Antigravity Image Processing & Pipeline Optimization
Author: Image Processing & Pipeline Optimization Engineer
Date: August 2026

Description:
    Production-ready deterministic filter module providing pre-segmentation
    blanking and post-segmentation validation for digitized herbarium vouchers.
    
    Herbarium specimen digitization pipelines frequently encounter severe
    non-botanical artifacts including:
      - Printed accession labels and annotation slips (high text density, 90° corners)
      - Mounting tape strips (high rectangularity)
      - Standardized calibration color charts (high HSV saturation)
      - Linear scale rulers (orthogonal edges, periodic tick marks)
      - Barcode stickers and institution stamps
      
    This module implements a rigorous 4-stage deterministic gatekeeper to ensure
    100% rejection or sterilization of artificial features while safely retaining
    authentic, curvilinear *Packera* leaf silhouettes for Fourier and geometric
    morphometrics.

Pipeline Stages:
    Stage 1: Pre-Emptive Layout Hard-Masking
             Hard zero-fills bounding boxes/polygons of known non-plant layout
             artifacts to background RGB [255, 255, 255] with a 10-pixel padding boundary.
    Stage 2: Post-Extraction Geometric Gatekeeper
             Applies mathematical filters:
             a. Rectangularity Filter: Area_mask / Area_minAreaRect (reject > 0.86)
             b. Corner Angle Analysis: Douglas-Peucker polygon approximation;
                rejects 4-vertex quadrilaterals with all internal angles in [80°, 100°].
             c. Solidity Filter: Area_mask / Area_convex_hull (require >= 0.72 for intact leaves).
    Stage 3: Spectral & Saturation Filter (Color Swatch Rejection)
             HSV color space inspection. Reclassifies candidates as 'color_chart'
             if high-saturation pixels (S > 0.45) exceed 15% of the mask area.
    Stage 4: Text & Edge Density Verification
             Laplacian gradient variance and high-frequency edge/stroke density
             measurement. Automatically routes printed label artifacts to
             'data/cropped_patches/annotations/' instead of leaf morphometric stores.

Usage:
    # As a Python library:
    from scripts.artifact_filter_gatekeeper import ArtifactFilterGatekeeper
    gatekeeper = ArtifactFilterGatekeeper()
    clean_sheet = gatekeeper.pre_emptive_hard_blanking(sheet_img, detections)
    result = gatekeeper.validate_candidate_leaf(patch_bgr, leaf_mask, catalog_number="NCU00412345")
    
    # Run comprehensive synthetic unit test suite:
    python scripts/artifact_filter_gatekeeper.py --test
===============================================================================
"""

import os
import sys
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE, DEFAULT_RAW_DIR, DEFAULT_CURATED_CSV, DEFAULT_OUTPUT_DIR, DEFAULT_CONFIG_PATH, DEFAULT_QC_DIR
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
import math
import json
import logging
import argparse
import unittest
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Dict, Optional, Union, Any

import cv2
import numpy as np

# Configure structured logging for production traceability
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [%(name)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("ArtifactFilterGatekeeper")


# =============================================================================
# DATA STRUCTURES & CONFIGURATION CLASSES
# =============================================================================

def to_dict(self) -> Dict[str, Any]:
        """Convert dataclass hierarchy to a JSON-serializable dictionary."""
        return asdict(self)


# =============================================================================
# GATEKEEPER IMPLEMENTATION CLASS
# =============================================================================

class ArtifactFilterGatekeeper:
    """
    Industrial-strength gatekeeper for pre-segmentation sheet sterilization
    and post-segmentation leaf candidate validation.
    """

    def __init__(
        self,
        padding_pixels: int = 10,
        background_fill_color: Tuple[int, int, int] = (255, 255, 255),
        max_rectangularity_threshold: float = 0.86,
        min_solidity_threshold: float = 0.72,
        douglas_peucker_epsilon_ratio: float = 0.02,
        orthogonal_angle_range: Tuple[float, float] = (80.0, 100.0),
        high_saturation_pixel_threshold: float = 0.45,
        max_color_swatch_saturation_ratio: float = 0.15,
        laplacian_text_variance_threshold: float = 120.0,
        canny_text_edge_density_threshold: float = 0.08,
        annotations_archive_dir: str = "data/cropped_patches/annotations"
    ):
        """
        Initialize the ArtifactFilterGatekeeper with calibrated deterministic thresholds.
        
        Args:
            padding_pixels: Margin added around detected artifact bounding boxes (default: 10 px).
            background_fill_color: Solid background color used for blanking (default: RGB [255, 255, 255]).
            max_rectangularity_threshold: Maximum allowable Area_mask / Area_minAreaRect before rejection (0.86).
            min_solidity_threshold: Minimum allowable Area_mask / Area_convex_hull for intact single leaves (0.72).
            douglas_peucker_epsilon_ratio: Epsilon factor for cv2.approxPolyDP relative to perimeter (0.02).
            orthogonal_angle_range: Angular tolerance in degrees for rectangular corner detection (80° to 100°).
            high_saturation_pixel_threshold: HSV Saturation threshold defining vibrant color pixels (0.45).
            max_color_swatch_saturation_ratio: Max fraction of high-saturation pixels allowed in leaf masks (0.15).
            laplacian_text_variance_threshold: Minimum Laplacian variance indicative of printed text lines (120.0).
            canny_text_edge_density_threshold: Minimum edge density ratio indicative of printed typography (0.08).
            annotations_archive_dir: Filesystem path to route and archive rejected text annotation slips.
        """
        self.padding_pixels = padding_pixels
        self.background_fill_color = background_fill_color
        self.max_rectangularity = max_rectangularity_threshold
        self.min_solidity = min_solidity_threshold
        self.dp_epsilon_ratio = douglas_peucker_epsilon_ratio
        self.orthogonal_angle_min, self.orthogonal_angle_max = orthogonal_angle_range
        self.high_sat_threshold = high_saturation_pixel_threshold
        self.max_color_swatch_ratio = max_color_swatch_saturation_ratio
        self.laplacian_text_var_threshold = laplacian_text_variance_threshold
        self.canny_text_edge_threshold = canny_text_edge_density_threshold
        self.annotations_archive_dir = Path(annotations_archive_dir)

        # Ensure archive destination exists
        self.annotations_archive_dir.mkdir(parents=True, exist_ok=True)
        
        # Target artifact categories subjected to hard blanking
        self.sterilization_categories = {
            "herbarium_label",
            "color_chart",
            "ruler_scale",
            "barcode_sticker",
            "institution_stamp",
            "mounting_tape",
            "envelope",
            "annotation_label",
            "packet"
        }

    # -------------------------------------------------------------------------
    # STAGE 1: PRE-EMPTIVE LAYOUT HARD-MASKING
    # -------------------------------------------------------------------------
    def pre_emptive_hard_blanking(
        self,
        sheet_image: np.ndarray,
        artifact_detections: List[Union[ArtifactDetection, Dict[str, Any], List[int], Tuple[int, int, int, int]]],
        is_rgb: bool = False
    ) -> np.ndarray:
        """
        Sterilize non-plant regions by hard zero-filling layout artifacts before segmentation.
        
        Applies a solid background fill (RGB/BGR [255, 255, 255]) with an expanded 10-pixel
        safety boundary over detected bounding boxes or polygons corresponding to labels,
        color charts, scale rulers, and barcodes.
        
        Args:
            sheet_image: High-resolution herbarium sheet image (H, W, 3) as a NumPy ndarray.
            artifact_detections: List of ArtifactDetection objects, bounding boxes, or polygons.
            is_rgb: Set True if sheet_image is in RGB channel order, False if BGR (OpenCV default).
            
        Returns:
            Sterilized sheet image with non-plant regions completely blanked to solid background.
        """
        # Create a contiguous copy to preserve original voucher integrity
        sterilized = sheet_image.copy()
        img_h, img_w = sterilized.shape[:2]

        fill_val = self.background_fill_color
        # Handle grayscale single-channel images if encountered
        if sterilized.ndim == 2:
            fill_val = 255

        pad = self.padding_pixels

        for det in artifact_detections:
            # Parse bounding box and polygon representations
            box = None
            polygon = None
            category = "herbarium_label"

            if isinstance(det, ArtifactDetection):
                box = det.box
                polygon = det.polygon
                category = det.category
            elif isinstance(det, dict):
                box = det.get("box") or det.get("bbox")
                polygon = det.get("polygon")
                category = det.get("category") or det.get("label", "herbarium_label")
            elif isinstance(det, (list, tuple)):
                if len(det) == 4 and all(isinstance(v, (int, float, np.integer)) for v in det):
                    # Direct [xmin, ymin, xmax, ymax] coordinate tuple
                    box = [int(v) for v in det]
                elif len(det) >= 3 and all(isinstance(pt, (list, tuple)) for pt in det):
                    # Direct polygon coordinate list [[x1, y1], [x2, y2], ...]
                    polygon = det

            # If polygon geometry is supplied, draw an expanded filled polygon
            if polygon is not None and len(polygon) >= 3:
                pts = np.array(polygon, dtype=np.int32).reshape((-1, 1, 2))
                # Dilate polygon by padding pixels using morphological mask dilation
                poly_mask = np.zeros((img_h, img_w), dtype=np.uint8)
                cv2.fillPoly(poly_mask, [pts], 255)
                if pad > 0:
                    kernel_size = 2 * pad + 1
                    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
                    poly_mask = cv2.dilate(poly_mask, kernel)
                sterilized[poly_mask > 0] = fill_val
                logger.debug("Sterilized polygon artifact '%s' with %d px padding", category, pad)
                continue

            # If bounding box is supplied, apply expanded rectangle blanking
            if box is not None and len(box) == 4:
                xmin, ymin, xmax, ymax = [int(round(v)) for v in box]
                
                # Expand box by padding pixels, clamping strictly to image bounds
                x1 = max(0, xmin - pad)
                y1 = max(0, ymin - pad)
                x2 = min(img_w, xmax + pad)
                y2 = min(img_h, ymax + pad)

                if x2 > x1 and y2 > y1:
                    sterilized[y1:y2, x1:x2] = fill_val
                    logger.debug(
                        "Hard-blanked bounding box [%d, %d, %d, %d] for artifact '%s'",
                        x1, y1, x2, y2, category
                    )

        return sterilized

    # -------------------------------------------------------------------------
    # STAGE 2: POST-EXTRACTION GEOMETRIC GATEKEEPER
    # -------------------------------------------------------------------------
    def compute_geometric_metrics(self, candidate_mask: np.ndarray) -> GeometricMetrics:
        """
        Compute deterministic geometric morphology metrics from a binary leaf candidate mask.
        
        Calculates:
          1. Rectangularity: Ratio of mask area to minimum bounding oriented rectangle area.
          2. Douglas-Peucker polygon approximation and 4-corner internal angle analysis.
          3. Solidity: Ratio of mask area to convex hull area.
          
        Args:
            candidate_mask: Binary single-channel uint8 mask (foreground=255/1, background=0).
            
        Returns:
            GeometricMetrics dataclass containing all morphological parameters.
        """
        # Ensure mask is binary uint8 (0 and 255)
        bin_mask = (candidate_mask > 0).astype(np.uint8) * 255
        mask_area = float(cv2.countNonZero(bin_mask))

        if mask_area < 1.0:
            return GeometricMetrics(
                mask_area=0.0,
                min_area_rect_area=0.0,
                rectangularity=0.0,
                convex_hull_area=0.0,
                solidity=0.0,
                num_approx_vertices=0,
                corner_angles_deg=[],
                is_rectangular=False,
                is_orthogonal_quad=False,
                is_valid_solidity=False
            )

        # Find external contours of candidate mask
        contours, _ = cv2.findContours(bin_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            return GeometricMetrics(
                mask_area=mask_area,
                min_area_rect_area=0.0,
                rectangularity=0.0,
                convex_hull_area=0.0,
                solidity=0.0,
                num_approx_vertices=0,
                corner_angles_deg=[],
                is_rectangular=False,
                is_orthogonal_quad=False,
                is_valid_solidity=False
            )

        # Select largest primary external contour
        cnt = max(contours, key=cv2.contourArea)
        cnt_area = float(cv2.contourArea(cnt))
        if cnt_area < 1.0:
            cnt_area = mask_area

        # 1. Minimum Bounding Oriented Rectangle & Rectangularity
        # minAreaRect returns ((center_x, center_y), (width, height), angle_deg)
        min_rect = cv2.minAreaRect(cnt)
        rect_w, rect_h = min_rect[1]
        min_area_rect_area = float(rect_w * rect_h)

        if min_area_rect_area > 0.0:
            rectangularity = float(mask_area / min_area_rect_area)
        else:
            rectangularity = 0.0

        is_rectangular = (rectangularity > self.max_rectangularity)

        # 2. Douglas-Peucker Polygon Approximation & Corner Angle Analysis
        perimeter = cv2.arcLength(cnt, True)
        epsilon = self.dp_epsilon_ratio * perimeter
        approx = cv2.approxPolyDP(cnt, epsilon, True)
        num_vertices = len(approx)

        corner_angles: List[float] = []
        is_orthogonal_quad = False

        if num_vertices == 4:
            # Reshape vertices to list of 2D points [[x0, y0], [x1, y1], [x2, y2], [x3, y3]]
            pts = approx.reshape((4, 2)).astype(np.float64)
            
            # Compute internal angle at each vertex i
            for i in range(4):
                p_prev = pts[(i - 1) % 4]
                p_curr = pts[i]
                p_next = pts[(i + 1) % 4]

                # Vector from curr to prev and curr to next
                v1 = p_prev - p_curr
                v2 = p_next - p_curr

                norm1 = np.linalg.norm(v1)
                norm2 = np.linalg.norm(v2)

                if norm1 > 1e-6 and norm2 > 1e-6:
                    dot_prod = np.dot(v1, v2)
                    cos_theta = dot_prod / (norm1 * norm2)
                    # Numerical clipping to prevent arccos domain errors
                    cos_theta = np.clip(cos_theta, -1.0, 1.0)
                    angle_deg = math.degrees(math.acos(cos_theta))
                else:
                    angle_deg = 0.0

                corner_angles.append(float(angle_deg))

            # Verify if ALL 4 corner angles fall strictly within orthogonal range [80°, 100°]
            if len(corner_angles) == 4 and all(
                self.orthogonal_angle_min <= ang <= self.orthogonal_angle_max for ang in corner_angles
            ):
                is_orthogonal_quad = True

        # 3. Convex Hull & Solidity Filter
        hull = cv2.convexHull(cnt)
        convex_hull_area = float(cv2.contourArea(hull))

        if convex_hull_area > 0.0:
            solidity = float(mask_area / convex_hull_area)
        else:
            solidity = 0.0

        is_valid_solidity = (solidity >= self.min_solidity)

        return GeometricMetrics(
            mask_area=mask_area,
            min_area_rect_area=min_area_rect_area,
            rectangularity=rectangularity,
            convex_hull_area=convex_hull_area,
            solidity=solidity,
            num_approx_vertices=num_vertices,
            corner_angles_deg=corner_angles,
            is_rectangular=is_rectangular,
            is_orthogonal_quad=is_orthogonal_quad,
            is_valid_solidity=is_valid_solidity
        )

    # -------------------------------------------------------------------------
    # STAGE 3: SPECTRAL & SATURATION FILTER (COLOR SWATCH REJECTION)
    # -------------------------------------------------------------------------
    def compute_spectral_metrics(
        self,
        candidate_patch: np.ndarray,
        candidate_mask: np.ndarray,
        is_rgb: bool = False
    ) -> SpectralMetrics:
        """
        Analyze the colorimetric and HSV saturation profile of foreground mask pixels.
        
        Authentic dried *Packera* herbarium vouchers are predominantly desaturated
        earth tones (ochres, tan, olive-brown; S < 0.30). Calibration color chart swatches
        (cyan, magenta, yellow, vibrant red/blue) exhibit high saturation (S > 0.45).
        
        Args:
            candidate_patch: Multi-channel image patch (H, W, 3) corresponding to the mask bounding box.
            candidate_mask: Binary single-channel uint8 mask (H, W) with foreground=255.
            is_rgb: Set True if candidate_patch is RGB, False if BGR.
            
        Returns:
            SpectralMetrics dataclass containing saturation and color distribution statistics.
        """
        if candidate_patch is None or candidate_patch.size == 0 or candidate_mask is None:
            return SpectralMetrics(
                mean_hue=0.0,
                mean_saturation=0.0,
                mean_value=0.0,
                high_saturation_pixel_count=0,
                high_saturation_ratio=0.0,
                is_color_swatch=False
            )

        # Convert to BGR if input is RGB for standard OpenCV HSV mapping
        if is_rgb:
            bgr_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_RGB2BGR)
        else:
            bgr_patch = candidate_patch

        # Convert patch to HSV color space (H in [0, 179], S in [0, 255], V in [0, 255])
        hsv_patch = cv2.cvtColor(bgr_patch, cv2.COLOR_BGR2HSV)

        # Extract foreground pixels under mask
        fg_indices = np.where(candidate_mask > 0)
        total_fg_pixels = len(fg_indices[0])

        if total_fg_pixels == 0:
            return SpectralMetrics(
                mean_hue=0.0,
                mean_saturation=0.0,
                mean_value=0.0,
                high_saturation_pixel_count=0,
                high_saturation_ratio=0.0,
                is_color_swatch=False
            )

        fg_hsv = hsv_patch[fg_indices]
        h_vals = fg_hsv[:, 0].astype(np.float64) * 2.0  # Scale OpenCV [0, 179] to degrees [0, 358]
        s_vals = fg_hsv[:, 1].astype(np.float64) / 255.0  # Normalize to [0.0, 1.0]
        v_vals = fg_hsv[:, 2].astype(np.float64) / 255.0  # Normalize to [0.0, 1.0]

        mean_hue = float(np.mean(h_vals))
        mean_saturation = float(np.mean(s_vals))
        mean_value = float(np.mean(v_vals))

        # Count pixels with saturation exceeding threshold (S > 0.45)
        high_sat_count = int(np.sum(s_vals > self.high_sat_threshold))
        high_sat_ratio = float(high_sat_count / total_fg_pixels)

        # Flag as color chart swatch if saturated pixels exceed 15% of total mask area
        is_color_swatch = (high_sat_ratio > self.max_color_swatch_ratio)

        return SpectralMetrics(
            mean_hue=mean_hue,
            mean_saturation=mean_saturation,
            mean_value=mean_value,
            high_saturation_pixel_count=high_sat_count,
            high_saturation_ratio=high_sat_ratio,
            is_color_swatch=is_color_swatch
        )

    # -------------------------------------------------------------------------
    # STAGE 4: TEXT & EDGE DENSITY VERIFICATION
    # -------------------------------------------------------------------------
    def compute_texture_metrics(
        self,
        candidate_patch: np.ndarray,
        candidate_mask: np.ndarray
    ) -> TextureMetrics:
        """
        Quantify interior high-frequency gradient variance and typographic edge stroke density.
        
        Labels and printed annotations possess sharp, high-contrast, bi-directional
        typographic glyph edges on white paper backgrounds, yielding high Laplacian variance
        and prominent horizontal text baseline alignments.
        
        Args:
            candidate_patch: Multi-channel or grayscale image patch (H, W, 3) or (H, W).
            candidate_mask: Binary single-channel uint8 mask (H, W).
            
        Returns:
            TextureMetrics dataclass summarizing gradient variance and text indicators.
        """
        if candidate_patch is None or candidate_patch.size == 0 or candidate_mask is None:
            return TextureMetrics(
                laplacian_variance=0.0,
                canny_edge_density=0.0,
                horizontal_stroke_density=0.0,
                vertical_stroke_density=0.0,
                is_printed_text=False
            )

        # Convert patch to single-channel grayscale
        if candidate_patch.ndim == 3:
            gray_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_BGR2GRAY)
        else:
            gray_patch = candidate_patch.copy()

        # Erode mask by 3 pixels to isolate interior texture and avoid perimeter border edges
        erode_kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5, 5))
        interior_mask = cv2.erode((candidate_mask > 0).astype(np.uint8) * 255, erode_kernel)
        interior_fg_pixels = int(cv2.countNonZero(interior_mask))

        # Fallback to full mask if candidate is too small to withstand erosion
        if interior_fg_pixels < 25:
            interior_mask = (candidate_mask > 0).astype(np.uint8) * 255
            interior_fg_pixels = int(cv2.countNonZero(interior_mask))

        if interior_fg_pixels == 0:
            return TextureMetrics(
                laplacian_variance=0.0,
                canny_edge_density=0.0,
                horizontal_stroke_density=0.0,
                vertical_stroke_density=0.0,
                is_printed_text=False
            )

        # 1. Laplacian Gradient Variance
        laplacian = cv2.Laplacian(gray_patch, cv2.CV_64F)
        interior_laplacian = laplacian[interior_mask > 0]
        laplacian_variance = float(np.var(interior_laplacian))

        # 2. Canny Edge Density
        canny_edges = cv2.Canny(gray_patch, 50, 150)
        interior_canny_count = int(np.sum((canny_edges > 0) & (interior_mask > 0)))
        canny_edge_density = float(interior_canny_count / interior_fg_pixels)

        # 3. Directional Sobel Typographic Stroke Analysis
        sobel_x = cv2.Sobel(gray_patch, cv2.CV_64F, 1, 0, ksize=3)
        sobel_y = cv2.Sobel(gray_patch, cv2.CV_64F, 0, 1, ksize=3)

        interior_sobel_x = np.abs(sobel_x[interior_mask > 0])
        interior_sobel_y = np.abs(sobel_y[interior_mask > 0])

        horizontal_stroke_density = float(np.mean(interior_sobel_y))
        vertical_stroke_density = float(np.mean(interior_sobel_x))

        # Decision rule: High Laplacian variance combined with high Canny edge density
        # or pronounced orthogonal gradient strokes indicates printed text
        is_printed_text = (
            (laplacian_variance > self.laplacian_text_var_threshold and
             canny_edge_density > self.canny_text_edge_threshold) or
            (canny_edge_density > 0.18 and laplacian_variance > 80.0)
        )

        return TextureMetrics(
            laplacian_variance=laplacian_variance,
            canny_edge_density=canny_edge_density,
            horizontal_stroke_density=horizontal_stroke_density,
            vertical_stroke_density=vertical_stroke_density,
            is_printed_text=is_printed_text
        )

    # -------------------------------------------------------------------------
    # COMPOSITE CANDIDATE VALIDATION & ANNOTATION ROUTING
    # -------------------------------------------------------------------------
    def validate_candidate_leaf(
        self,
        candidate_patch: np.ndarray,
        candidate_mask: np.ndarray,
        catalog_number: str = "UNKNOWN",
        patch_id: str = "0",
        candidate_class: str = "basal_leaf",
        is_rgb: bool = False
    ) -> FilterResult:
        """
        Execute the multi-stage deterministic gatekeeper against an extracted leaf candidate.
        
        Evaluates:
          1. Geometric Rectangularity, Corner Angles, and Solidity (Stage 2)
          2. Spectral Saturation / Color Swatch reclassification (Stage 3)
          3. Laplacian Text Variance & Typographic Edge Density (Stage 4)
          
        If a candidate is rejected due to printed text, it is routed and archived to
        'data/cropped_patches/annotations/' for OCR and metadata indexing.
        
        Args:
            candidate_patch: Crop of candidate image region (H, W, 3).
            candidate_mask: Binary candidate silhouette mask (H, W).
            catalog_number: Specimen voucher catalog number identifier (e.g. 'NCU00012345').
            patch_id: Unique identifier or index for the candidate patch.
            candidate_class: Nominal classification label (default: 'basal_leaf').
            is_rgb: Set True if candidate_patch is RGB, False if BGR.
            
        Returns:
            FilterResult containing complete diagnostic telemetry and pass/fail determination.
        """
        # 1. Compute all stage metrics
        geo_metrics = self.compute_geometric_metrics(candidate_mask)
        spec_metrics = self.compute_spectral_metrics(candidate_patch, candidate_mask, is_rgb=is_rgb)
        tex_metrics = self.compute_texture_metrics(candidate_patch, candidate_mask)

        rejection_reason = None
        reclassified_category = None
        routed_path = None

        # ---------------------------------------------------------------------
        # Evaluation Hierarchy
        # ---------------------------------------------------------------------
        
        # Check A: Spectral Reclassification (Color Swatches)
        if spec_metrics.is_color_swatch:
            rejection_reason = "REJECT_HIGH_SATURATION_COLOR_SWATCH"
            reclassified_category = "color_chart"

        # Check B: Text & Typography Verification (Stage 4)
        elif tex_metrics.is_printed_text:
            rejection_reason = "REJECT_PRINTED_TYPOGRAPHY_TEXT"
            reclassified_category = "herbarium_label"

        # Check C: Rectangularity Filter (Stage 2a)
        elif geo_metrics.is_rectangular:
            rejection_reason = f"REJECT_RECTANGULARITY_EXCEEDED ({geo_metrics.rectangularity:.3f} > {self.max_rectangularity})"
            reclassified_category = "mounting_tape_or_card"

        # Check D: Douglas-Peucker Orthogonal 4-Corner Quadrilateral (Stage 2b)
        elif geo_metrics.is_orthogonal_quad:
            rejection_reason = "REJECT_ORTHOGONAL_QUADRILATERAL_CORNERS"
            reclassified_category = "herbarium_label"

        # Check E: Solidity Gatekeeper (Stage 2c)
        elif not geo_metrics.is_valid_solidity:
            rejection_reason = f"REJECT_LOW_SOLIDITY_COMPOSITE_CLUMP ({geo_metrics.solidity:.3f} < {self.min_solidity})"
            reclassified_category = "severed_or_clumped_vegetation"

        # Determine overall validity
        is_valid = (rejection_reason is None)

        if is_valid:
            status = "VALID_LEAF"
        else:
            # Route text annotations and labels to dedicated annotation archive
            if reclassified_category in ("herbarium_label", "annotation_label") or tex_metrics.is_printed_text:
                status = "ROUTED_ANNOTATION"
                routed_path = self._route_annotation_patch(
                    candidate_patch=candidate_patch,
                    candidate_mask=candidate_mask,
                    catalog_number=catalog_number,
                    patch_id=patch_id,
                    metrics={
                        "geometric": asdict(geo_metrics),
                        "spectral": asdict(spec_metrics),
                        "texture": asdict(tex_metrics),
                        "rejection_reason": rejection_reason
                    },
                    is_rgb=is_rgb
                )
            else:
                status = "REJECTED_ARTIFACT"

        logger.info(
            "Candidate [%s / Patch %s]: Status=%s | Reason=%s | Rect=%.3f | Solid=%.3f | Sat=%.3f | LapVar=%.1f",
            catalog_number, patch_id, status, rejection_reason,
            geo_metrics.rectangularity, geo_metrics.solidity,
            spec_metrics.high_saturation_ratio, tex_metrics.laplacian_variance
        )

        return FilterResult(
            is_valid=is_valid,
            status=status,
            primary_rejection_reason=rejection_reason,
            reclassified_category=reclassified_category,
            geometric_metrics=geo_metrics,
            spectral_metrics=spec_metrics,
            texture_metrics=tex_metrics,
            routed_file_path=routed_path
        )

    def _route_annotation_patch(
        self,
        candidate_patch: np.ndarray,
        candidate_mask: np.ndarray,
        catalog_number: str,
        patch_id: str,
        metrics: Dict[str, Any],
        is_rgb: bool = False
    ) -> str:
        """
        Archive rejected text/label patches to the annotation storage directory.
        
        Args:
            candidate_patch: Image patch ndarray.
            candidate_mask: Binary mask ndarray.
            catalog_number: Specimen identifier.
            patch_id: Patch index.
            metrics: Diagnostic metrics dictionary.
            is_rgb: Boolean channel format flag.
            
        Returns:
            Absolute file path of the archived image patch.
        """
        base_name = f"{catalog_number}_patch{patch_id}_annotation"
        img_out_path = self.annotations_archive_dir / f"{base_name}.png"
        meta_out_path = self.annotations_archive_dir / f"{base_name}_meta.json"

        # Prepare patch for saving via OpenCV (BGR order)
        if candidate_patch is not None and candidate_patch.size > 0:
            if is_rgb and candidate_patch.ndim == 3:
                save_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_RGB2BGR)
            else:
                save_patch = candidate_patch
            cv2.imwrite(str(img_out_path), save_patch)

        # Write accompanying metadata JSON
        with open(meta_out_path, "w", encoding="utf-8") as f:
            json.dump({
                "catalog_number": catalog_number,
                "patch_id": patch_id,
                "archived_image": str(img_out_path),
                "metrics": metrics
            }, f, indent=2)

        logger.debug("Archived rejected annotation patch to '%s'", img_out_path)
        return str(img_out_path.resolve())


# =============================================================================
# SYNTHETIC TEST FIXTURE GENERATORS
# =============================================================================

def generate_synthetic_leaf(
    img_size: Tuple[int, int] = (256, 256),
    blade_radii: Tuple[int, int] = (60, 30),
    rotation_deg: float = 25.0,
    color_bgr: Tuple[int, int, int] = (75, 95, 105)  # Realistic low-saturation dried olive-tan (S ~ 0.28)
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate an authentic synthetic elliptical Packera leaf blade with natural petiole.
    
    Returns:
        Tuple of (patch_bgr, mask_uint8).
    """
    h, w = img_size
    patch = np.full((h, w, 3), (250, 248, 245), dtype=np.uint8)  # Herbarium paper bg
    mask = np.zeros((h, w), dtype=np.uint8)

    center = (w // 2, h // 2)
    # Draw smooth organic elliptical leaf blade
    cv2.ellipse(mask, center, blade_radii, rotation_deg, 0, 360, 255, -1)

    # Attach tapered curvilinear petiole extending from leaf base
    angle_rad = math.radians(rotation_deg)
    petiole_start = (
        int(center[0] + blade_radii[0] * 0.9 * math.cos(angle_rad)),
        int(center[1] + blade_radii[0] * 0.9 * math.sin(angle_rad))
    )
    petiole_end = (
        int(center[0] + (blade_radii[0] + 50) * math.cos(angle_rad + 0.15)),
        int(center[1] + (blade_radii[0] + 50) * math.sin(angle_rad + 0.15))
    )
    cv2.line(mask, petiole_start, petiole_end, 255, thickness=6)

    # Paint realistic plant coloration with subtle venation (desaturated earth tones S < 0.30)
    patch[mask > 0] = color_bgr
    # Subtle primary midrib line
    cv2.line(patch, center, petiole_end, (65, 80, 90), thickness=2)

    return patch, mask


def generate_synthetic_herbarium_label(
    img_size: Tuple[int, int] = (256, 256),
    box_size: Tuple[int, int] = (200, 120),
    rotation_deg: float = 0.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a synthetic rectangular herbarium label with printed typographic text lines.
    
    Returns:
        Tuple of (patch_bgr, mask_uint8).
    """
    h, w = img_size
    patch = np.full((h, w, 3), (255, 255, 255), dtype=np.uint8)
    mask = np.zeros((h, w), dtype=np.uint8)

    bw, bh = box_size
    x1 = (w - bw) // 2
    y1 = (h - bh) // 2
    x2 = x1 + bw
    y2 = y1 + bh

    # Draw rigid rectangular label mask
    cv2.rectangle(mask, (x1, y1), (x2, y2), 255, -1)
    
    # Fill label card background (slightly off-white paper)
    patch[mask > 0] = (245, 245, 240)
    
    # Draw dark label border
    cv2.rectangle(patch, (x1, y1), (x2, y2), (40, 40, 40), 2)

    # Simulate dense printed text lines (typographic glyphs)
    for row_y in range(y1 + 20, y2 - 15, 18):
        cv2.putText(
            patch,
            "PLANTS OF NORTH CAROLINA - Packera dubia",
            (x1 + 10, row_y),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.35,
            (15, 15, 15),
            1,
            cv2.LINE_AA
        )

    return patch, mask


def generate_synthetic_color_chart_swatch(
    img_size: Tuple[int, int] = (256, 256),
    box_size: Tuple[int, int] = (160, 160)
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a synthetic high-saturation calibration color chart swatch (vibrant cyan/magenta/yellow).
    
    Returns:
        Tuple of (patch_bgr, mask_uint8).
    """
    h, w = img_size
    patch = np.full((h, w, 3), (255, 255, 255), dtype=np.uint8)
    mask = np.zeros((h, w), dtype=np.uint8)

    bw, bh = box_size
    x1 = (w - bw) // 2
    y1 = (h - bh) // 2

    # Draw 4 vibrant color quadrants: Cyan, Magenta, Yellow, Saturated Red
    colors = [
        (255, 255, 0),    # Cyan in BGR (B=255, G=255, R=0)
        (255, 0, 255),    # Magenta in BGR
        (0, 255, 255),    # Yellow in BGR
        (0, 0, 255)       # Pure Red in BGR
    ]

    half_w = bw // 2
    half_h = bh // 2

    quads = [
        ((x1, y1), (x1 + half_w, y1 + half_h), colors[0]),
        ((x1 + half_w, y1), (x1 + bw, y1 + half_h), colors[1]),
        ((x1, y1 + half_h), (x1 + half_w, y1 + bh), colors[2]),
        ((x1 + half_w, y1 + half_h), (x1 + bw, y1 + bh), colors[3]),
    ]

    for (p1, p2, col) in quads:
        cv2.rectangle(mask, p1, p2, 255, -1)
        cv2.rectangle(patch, p1, p2, col, -1)

    return patch, mask


def generate_synthetic_scale_ruler(
    img_size: Tuple[int, int] = (256, 256),
    box_size: Tuple[int, int] = (40, 200)
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a synthetic linear scale calibration ruler with millimeter tick marks.
    
    Returns:
        Tuple of (patch_bgr, mask_uint8).
    """
    h, w = img_size
    patch = np.full((h, w, 3), (255, 255, 255), dtype=np.uint8)
    mask = np.zeros((h, w), dtype=np.uint8)

    bw, bh = box_size
    x1 = (w - bw) // 2
    y1 = (h - bh) // 2
    x2 = x1 + bw
    y2 = y1 + bh

    cv2.rectangle(mask, (x1, y1), (x2, y2), 255, -1)
    patch[mask > 0] = (240, 240, 240)
    cv2.rectangle(patch, (x1, y1), (x2, y2), (0, 0, 0), 2)

    # Periodic ruler tick marks (high orthogonal gradient)
    for tick_y in range(y1 + 5, y2 - 5, 8):
        tick_len = 15 if (tick_y - y1) % 40 == 0 else 8
        cv2.line(patch, (x1, tick_y), (x1 + tick_len, tick_y), (0, 0, 0), 2)

    return patch, mask


def generate_synthetic_mounting_tape(
    img_size: Tuple[int, int] = (256, 256),
    box_size: Tuple[int, int] = (180, 45),
    rotation_deg: float = 15.0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a synthetic rectangular mounting tape strip.
    
    Returns:
        Tuple of (patch_bgr, mask_uint8).
    """
    h, w = img_size
    patch = np.full((h, w, 3), (255, 255, 255), dtype=np.uint8)
    mask = np.zeros((h, w), dtype=np.uint8)

    center = (w // 2, h // 2)
    # Create rotated rectangle box contour
    rect = (center, box_size, rotation_deg)
    box_pts = cv2.boxPoints(rect).astype(np.int32)

    cv2.fillPoly(mask, [box_pts], 255)
    # Translucent yellowish/tan tape tint
    patch[mask > 0] = (190, 220, 235)

    return patch, mask


def generate_synthetic_clumped_rosette(
    img_size: Tuple[int, int] = (256, 256)
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate a synthetic multi-leaf fused rosette clump exhibiting severe concavities (solidity < 0.65).
    
    Returns:
        Tuple of (patch_bgr, mask_uint8).
    """
    h, w = img_size
    patch = np.full((h, w, 3), (250, 248, 245), dtype=np.uint8)
    mask = np.zeros((h, w), dtype=np.uint8)

    center = (w // 2, h // 2)
    # Draw narrow radiating lobes extending outward to produce deep concavities
    angles = [0.0, 72.0, 144.0, 216.0, 288.0]
    for ang in angles:
        rad = math.radians(ang)
        lobe_center = (
            int(center[0] + 55 * math.cos(rad)),
            int(center[1] + 55 * math.sin(rad))
        )
        cv2.ellipse(mask, lobe_center, (45, 12), ang, 0, 360, 255, -1)
    
    # Small central connecting hub
    cv2.circle(mask, center, 14, 255, -1)

    patch[mask > 0] = (75, 95, 105)
    return patch, mask


# =============================================================================
# SYNTHETIC TEST SUITE & VERIFICATION HARNESS
# =============================================================================

def run_synthetic_test_suite() -> bool:
    """
    Execute 100% automated verification against synthetic herbarium artifacts and authentic leaves.
    
    Verifies:
      1. Stage 1: Pre-emptive hard-blanking completely sterilizes layout artifact regions with 10px padding.
      2. Stage 2: Geometric filter rejects 100% of synthetic rectangular labels, tapes, and low-solidity clumps.
      3. Stage 3: Spectral filter rejects 100% of vibrant color charts (saturation > 0.45 on >15% of area).
      4. Stage 4: Text/edge density verification detects printed typography and routes to annotations/.
      5. End-to-End: 100% retention of authentic elliptical leaves and 100% rejection of artificial edge cases.
      
    Returns:
        True if all verification assertions pass, False otherwise.
    """
    print("\n" + "=" * 80)
    print("RUNNING PRODUCTION GATEKEEPER SYNTHETIC VERIFICATION SUITE")
    print("=" * 80)

    gatekeeper = ArtifactFilterGatekeeper(annotations_archive_dir="data/cropped_patches/annotations")

    test_passed = True
    total_tests = 0
    passed_tests = 0

    # -------------------------------------------------------------------------
    # TEST 1: Stage 1 Pre-Emptive Hard Blanking Sterilization
    # -------------------------------------------------------------------------
    total_tests += 1
    print("\n[TEST 1] Testing Stage 1 Pre-Emptive Hard Blanking (10-pixel padding boundary)...")
    
    sheet_canvas = np.full((1000, 800, 3), (120, 140, 160), dtype=np.uint8)
    # Define artifact bounding boxes
    mock_detections = [
        ArtifactDetection(box=[50, 50, 250, 150], category="herbarium_label"),
        ArtifactDetection(box=[500, 600, 700, 750], category="color_chart"),
        ArtifactDetection(box=[100, 700, 150, 950], category="ruler_scale"),
        ArtifactDetection(box=[600, 50, 750, 120], category="barcode_sticker")
    ]

    sterilized_sheet = gatekeeper.pre_emptive_hard_blanking(sheet_canvas, mock_detections, is_rgb=False)

    # Verify that regions including 10px padding are completely filled with (255, 255, 255)
    all_blanked = True
    for det in mock_detections:
        x1, y1, x2, y2 = det.box
        # Check interior and expanded padding pixels
        pad_x1 = max(0, x1 - 10)
        pad_y1 = max(0, y1 - 10)
        pad_x2 = min(800, x2 + 10)
        pad_y2 = min(1000, y2 + 10)
        
        region = sterilized_sheet[pad_y1:pad_y2, pad_x1:pad_x2]
        if not np.all(region == 255):
            all_blanked = False
            print(f"  FAILED: Artifact {det.category} was not fully sterilized to RGB [255, 255, 255]!")

    if all_blanked:
        print("  PASSED: 100% of layout artifact boxes and 10px margins completely hard-blanked.")
        passed_tests += 1
    else:
        test_passed = False

    # -------------------------------------------------------------------------
    # TEST 2: Genuine Packera Leaf Blade Validation (Should Pass)
    # -------------------------------------------------------------------------
    total_tests += 1
    print("\n[TEST 2] Testing Genuine Authentic Packera Basal Leaf Silhouette...")
    leaf_patch, leaf_mask = generate_synthetic_leaf(blade_radii=(65, 32), rotation_deg=20.0)
    
    leaf_res = gatekeeper.validate_candidate_leaf(
        leaf_patch, leaf_mask, catalog_number="NCU_SYNTHETIC_001", patch_id="leaf_01"
    )

    if leaf_res.is_valid and leaf_res.status == "VALID_LEAF":
        print(f"  PASSED: Authentic leaf retained (Rect={leaf_res.geometric_metrics.rectangularity:.3f}, "
              f"Solid={leaf_res.geometric_metrics.solidity:.3f}, Sat={leaf_res.spectral_metrics.high_saturation_ratio:.3f}).")
        passed_tests += 1
    else:
        print(f"  FAILED: Authentic leaf erroneously rejected! Reason: {leaf_res.primary_rejection_reason}")
        test_passed = False

    # -------------------------------------------------------------------------
    # TEST 3: Synthetic Herbarium Label Rejection & Annotation Routing
    # -------------------------------------------------------------------------
    total_tests += 1
    print("\n[TEST 3] Testing Synthetic Herbarium Label with Printed Typography...")
    label_patch, label_mask = generate_synthetic_herbarium_label()
    
    label_res = gatekeeper.validate_candidate_leaf(
        label_patch, label_mask, catalog_number="NCU_SYNTHETIC_002", patch_id="label_01"
    )

    if (not label_res.is_valid) and (label_res.status in ("ROUTED_ANNOTATION", "REJECTED_ARTIFACT")):
        print(f"  PASSED: Label correctly rejected and routed. Status={label_res.status}, "
              f"Reason={label_res.primary_rejection_reason}, "
              f"RoutedPath={label_res.routed_file_path}")
        passed_tests += 1
    else:
        print(f"  FAILED: Label was not rejected! is_valid={label_res.is_valid}")
        test_passed = False

    # -------------------------------------------------------------------------
    # TEST 4: Synthetic Color Calibration Chart Swatch Rejection (Stage 3)
    # -------------------------------------------------------------------------
    total_tests += 1
    print("\n[TEST 4] Testing Calibration Color Chart Swatch Rejection (HSV Saturation > 0.45)...")
    chart_patch, chart_mask = generate_synthetic_color_chart_swatch()
    
    chart_res = gatekeeper.validate_candidate_leaf(
        chart_patch, chart_mask, catalog_number="NCU_SYNTHETIC_003", patch_id="chart_01"
    )

    if (not chart_res.is_valid) and (chart_res.reclassified_category == "color_chart"):
        print(f"  PASSED: Color chart swatch rejected and reclassified. "
              f"HighSatRatio={chart_res.spectral_metrics.high_saturation_ratio:.3f} > 0.15, "
              f"Reason={chart_res.primary_rejection_reason}")
        passed_tests += 1
    else:
        print(f"  FAILED: Color swatch was not correctly reclassified! Reclass={chart_res.reclassified_category}")
        test_passed = False

    # -------------------------------------------------------------------------
    # TEST 5: Synthetic Scale Calibration Ruler Rejection
    # -------------------------------------------------------------------------
    total_tests += 1
    print("\n[TEST 5] Testing Linear Scale Ruler Rejection (Orthogonal Quadrilateral & Rectangularity)...")
    ruler_patch, ruler_mask = generate_synthetic_scale_ruler()
    
    ruler_res = gatekeeper.validate_candidate_leaf(
        ruler_patch, ruler_mask, catalog_number="NCU_SYNTHETIC_004", patch_id="ruler_01"
    )

    if not ruler_res.is_valid:
        print(f"  PASSED: Scale ruler rejected. Reason={ruler_res.primary_rejection_reason}, "
              f"Rect={ruler_res.geometric_metrics.rectangularity:.3f}")
        passed_tests += 1
    else:
        print(f"  FAILED: Scale ruler was not rejected! is_valid={ruler_res.is_valid}")
        test_passed = False

    # -------------------------------------------------------------------------
    # TEST 6: Synthetic Mounting Tape Strip Rejection (Stage 2a Rectangularity)
    # -------------------------------------------------------------------------
    total_tests += 1
    print("\n[TEST 6] Testing Mounting Tape Strip Rejection (Rectangularity > 0.86)...")
    tape_patch, tape_mask = generate_synthetic_mounting_tape(rotation_deg=18.0)
    
    tape_res = gatekeeper.validate_candidate_leaf(
        tape_patch, tape_mask, catalog_number="NCU_SYNTHETIC_005", patch_id="tape_01"
    )

    if (not tape_res.is_valid) and ("REJECT_RECTANGULARITY_EXCEEDED" in str(tape_res.primary_rejection_reason) or
                                    "REJECT_ORTHOGONAL_QUADRILATERAL" in str(tape_res.primary_rejection_reason)):
        print(f"  PASSED: Tape strip rejected by geometric filter. Reason={tape_res.primary_rejection_reason}")
        passed_tests += 1
    else:
        print(f"  FAILED: Tape strip was not rejected! is_valid={tape_res.is_valid}, Reason={tape_res.primary_rejection_reason}")
        test_passed = False

    # -------------------------------------------------------------------------
    # TEST 7: Synthetic Fused Clump Rejection (Stage 2c Solidity < 0.72)
    # -------------------------------------------------------------------------
    total_tests += 1
    print("\n[TEST 7] Testing Rosette Clump Rejection (Solidity < 0.72)...")
    clump_patch, clump_mask = generate_synthetic_clumped_rosette()
    
    clump_res = gatekeeper.validate_candidate_leaf(
        clump_patch, clump_mask, catalog_number="NCU_SYNTHETIC_006", patch_id="clump_01"
    )

    if (not clump_res.is_valid) and ("REJECT_LOW_SOLIDITY" in str(clump_res.primary_rejection_reason)):
        print(f"  PASSED: Multi-leaf rosette clump rejected. Reason={clump_res.primary_rejection_reason}, "
              f"Solidity={clump_res.geometric_metrics.solidity:.3f}")
        passed_tests += 1
    else:
        print(f"  FAILED: Rosette clump was not rejected! is_valid={clump_res.is_valid}")
        test_passed = False

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print(f"SYNTHETIC VERIFICATION SUITE COMPLETE: {passed_tests}/{total_tests} TESTS PASSED ({(passed_tests/total_tests)*100.0:.1f}%)")
    print("=" * 80 + "\n")

    return test_passed


# =============================================================================
# UNITTEST TEST CASE CLASS FOR AUTOMATED TEST RUNNERS
# =============================================================================

class TestArtifactFilterGatekeeper(unittest.TestCase):
    """
    Standard unittest.TestCase harness enabling automated test runners (pytest, unittest)
    to discover, execute, and assert production gatekeeper invariants across all 4 stages.
    """

    def setUp(self):
        """Initialize gatekeeper instance with deterministic thresholds for unit tests."""
        self.gatekeeper = ArtifactFilterGatekeeper(
            padding_pixels=10,
            background_fill_color=(255, 255, 255),
            max_rectangularity_threshold=0.86,
            min_solidity_threshold=0.72,
            douglas_peucker_epsilon_ratio=0.02,
            orthogonal_angle_range=(80.0, 100.0),
            high_saturation_pixel_threshold=0.45,
            max_color_swatch_saturation_ratio=0.15,
            laplacian_text_variance_threshold=120.0,
            canny_text_edge_density_threshold=0.08,
            annotations_archive_dir="data/cropped_patches/annotations"
        )

    def test_stage1_pre_emptive_hard_blanking(self):
        """Verify Stage 1 layout hard-masking completely zeroes out artifacts with 10px padding."""
        canvas = np.full((1000, 800, 3), (120, 140, 160), dtype=np.uint8)
        detections = [
            ArtifactDetection(box=[50, 50, 250, 150], category="herbarium_label"),
            ArtifactDetection(box=[500, 600, 700, 750], category="color_chart"),
            ArtifactDetection(box=[100, 700, 150, 950], category="ruler_scale"),
            ArtifactDetection(box=[600, 50, 750, 120], category="barcode_sticker")
        ]
        sterilized = self.gatekeeper.pre_emptive_hard_blanking(canvas, detections, is_rgb=False)
        
        for det in detections:
            x1, y1, x2, y2 = det.box
            pad_x1 = max(0, x1 - 10)
            pad_y1 = max(0, y1 - 10)
            pad_x2 = min(800, x2 + 10)
            pad_y2 = min(1000, y2 + 10)
            region = sterilized[pad_y1:pad_y2, pad_x1:pad_x2]
            # Ensure 100% of region is solid RGB [255, 255, 255]
            self.assertTrue(np.all(region == 255), f"Artifact {det.category} region was not fully blanked.")

    def test_stage2_authentic_leaf_retention(self):
        """Verify authentic elliptic Packera leaf blade is retained by all geometric and spectral filters."""
        leaf_patch, leaf_mask = generate_synthetic_leaf(blade_radii=(65, 32), rotation_deg=20.0)
        res = self.gatekeeper.validate_candidate_leaf(
            leaf_patch, leaf_mask, catalog_number="NCU_UNITTEST_001", patch_id="leaf_01"
        )
        self.assertTrue(res.is_valid, "Authentic Packera leaf must pass all gatekeeper checks.")
        self.assertEqual(res.status, "VALID_LEAF")
        self.assertIsNone(res.primary_rejection_reason)
        self.assertGreaterEqual(res.geometric_metrics.solidity, 0.72)
        self.assertLessEqual(res.geometric_metrics.rectangularity, 0.86)

    def test_stage2a_rectangularity_rejection(self):
        """Verify rectangular mounting tape is rejected by Rectangularity > 0.86 threshold."""
        tape_patch, tape_mask = generate_synthetic_mounting_tape(rotation_deg=15.0)
        res = self.gatekeeper.validate_candidate_leaf(
            tape_patch, tape_mask, catalog_number="NCU_UNITTEST_002", patch_id="tape_01"
        )
        self.assertFalse(res.is_valid, "Mounting tape must be rejected by geometric filter.")
        self.assertIn("REJECT_RECTANGULARITY_EXCEEDED", str(res.primary_rejection_reason))

    def test_stage2c_solidity_rejection(self):
        """Verify multi-leaf fused clump is rejected by Solidity < 0.72 threshold."""
        clump_patch, clump_mask = generate_synthetic_clumped_rosette()
        res = self.gatekeeper.validate_candidate_leaf(
            clump_patch, clump_mask, catalog_number="NCU_UNITTEST_003", patch_id="clump_01"
        )
        self.assertFalse(res.is_valid, "Low-solidity rosette clump must be rejected.")
        self.assertIn("REJECT_LOW_SOLIDITY", str(res.primary_rejection_reason))

    def test_stage3_spectral_saturation_reclassification(self):
        """Verify high-saturation color chart swatches are reclassified to color_chart."""
        chart_patch, chart_mask = generate_synthetic_color_chart_swatch()
        res = self.gatekeeper.validate_candidate_leaf(
            chart_patch, chart_mask, catalog_number="NCU_UNITTEST_004", patch_id="chart_01"
        )
        self.assertFalse(res.is_valid, "Vibrant color chart must be rejected.")
        self.assertEqual(res.reclassified_category, "color_chart")
        self.assertEqual(res.primary_rejection_reason, "REJECT_HIGH_SATURATION_COLOR_SWATCH")

    def test_stage4_printed_typography_rejection_and_routing(self):
        """Verify printed text labels are detected via Laplacian/Sobel gradients and routed to annotations/."""
        label_patch, label_mask = generate_synthetic_herbarium_label()
        res = self.gatekeeper.validate_candidate_leaf(
            label_patch, label_mask, catalog_number="NCU_UNITTEST_005", patch_id="label_01"
        )
        self.assertFalse(res.is_valid, "Printed text label must be rejected.")
        self.assertEqual(res.status, "ROUTED_ANNOTATION")
        self.assertIsNotNone(res.routed_file_path)
        self.assertTrue(Path(res.routed_file_path).exists(), "Routed annotation image must exist on disk.")


# =============================================================================
# CLI ENTRY POINT
# =============================================================================

def main():
    """Command-line execution interface for artifact filtering and batch validation."""
    parser = argparse.ArgumentParser(
        description="Production Artifact Filter Gatekeeper for Botanical Herbarium Morphometrics."
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Execute synthetic unit test suite verifying 100%% rejection of artifacts and retention of leaves."
    )
    parser.add_argument(
        "--archive-dir",
        type=str,
        default="data/cropped_patches/annotations",
        help="Target archive directory for routed text and label patches."
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=10,
        help="Padding boundary in pixels for Stage 1 hard blanking."
    )
    parser.add_argument(
        "--rect-threshold",
        type=float,
        default=0.86,
        help="Maximum rectangularity threshold for Stage 2a geometric filter."
    )
    parser.add_argument(
        "--solidity-threshold",
        type=float,
        default=0.72,
        help="Minimum solidity threshold for Stage 2c intact leaf filter."
    )

    args = parser.parse_args()

    if args.test or len(sys.argv) == 1:
        success = run_synthetic_test_suite()
        sys.exit(0 if success else 1)
    else:
        logger.info("Artifact Filter Gatekeeper initialized with archive dir: %s", args.archive_dir)


if __name__ == "__main__":
    main()
