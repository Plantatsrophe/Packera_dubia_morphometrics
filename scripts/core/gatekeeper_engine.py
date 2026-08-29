
import os
import sys
import logging
import math
import numpy as np
import cv2
import json
import random
import unittest
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass, asdict

# Common imports
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
from scripts.core.tiling_utils import HerbariumAnnotation

logger = setup_logging()

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


