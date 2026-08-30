"""
scripts/core/gatekeeper_engine.py
=================================
Industrial-strength gatekeeper for pre-segmentation sheet sterilization
and post-segmentation leaf candidate validation.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import cv2
import numpy as np

from scripts.core.config import CLASS_COLORS_BGR, CLASS_MAP, CLASS_NAMES
from scripts.core.data_structures import (
    ArtifactDetection,
    FilterResult,
    GeometricMetrics,
    SpectralMetrics,
    TextureMetrics,
)
from scripts.core.gatekeeper_metrics import (
    compute_geometric_metrics,
    compute_spectral_metrics,
    compute_texture_metrics,
)
from scripts.core.logger import setup_logging

logger = setup_logging()


class ArtifactFilterGatekeeper:
    """
    Industrial-strength gatekeeper for pre-segmentation sheet sterilization
    and post-segmentation leaf candidate validation.
    """

    def __init__(
        self,
        padding_pixels: int = 15,
        background_fill_color: Tuple[int, int, int] = (255, 255, 255),
        max_rectangularity_threshold: float = 0.86,
        min_solidity_threshold: float = 0.72,
        douglas_peucker_epsilon_ratio: float = 0.02,
        orthogonal_angle_range: Tuple[float, float] = (80.0, 100.0),
        high_saturation_pixel_threshold: float = 0.45,
        max_color_swatch_saturation_ratio: float = 0.15,
        paper_mean_val_threshold: float = 205.0,
        paper_max_sat_threshold: float = 35.0,
        laplacian_text_variance_threshold: float = 450.0,
        canny_text_edge_density_threshold: float = 0.15,
        annotations_archive_dir: str = "data/cropped_patches/annotations"
    ):
        """
        Initialize ArtifactFilterGatekeeper with calibrated deterministic thresholds.
        """
        self.padding_pixels = padding_pixels
        self.background_fill_color = background_fill_color
        self.max_rectangularity = max_rectangularity_threshold
        self.min_solidity = min_solidity_threshold
        self.dp_epsilon_ratio = douglas_peucker_epsilon_ratio
        self.orthogonal_angle_min, self.orthogonal_angle_max = orthogonal_angle_range
        self.high_sat_threshold = high_saturation_pixel_threshold
        self.max_color_swatch_ratio = max_color_swatch_saturation_ratio
        self.paper_mean_val_threshold = paper_mean_val_threshold
        self.paper_max_sat_threshold = paper_max_sat_threshold
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
        is_rgb: bool = False,
        padding_pixels: Optional[int] = None
    ) -> np.ndarray:
        """
        Sterilize non-plant regions by hard zero-filling layout artifacts before segmentation.
        """
        sterilized = sheet_image.copy()
        img_h, img_w = sterilized.shape[:2]

        fill_val = self.background_fill_color
        if sterilized.ndim == 2:
            fill_val = 255

        pad = self.padding_pixels if padding_pixels is None else padding_pixels

        for det in artifact_detections:
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
                    box = [int(v) for v in det]
                elif len(det) >= 3 and all(isinstance(pt, (list, tuple)) for pt in det):
                    polygon = det

            if polygon is not None and len(polygon) >= 3:
                pts = np.array(polygon, dtype=np.int32).reshape((-1, 1, 2))
                poly_mask = np.zeros((img_h, img_w), dtype=np.uint8)
                cv2.fillPoly(poly_mask, [pts], 255)
                if pad > 0:
                    kernel_size = 2 * pad + 1
                    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (kernel_size, kernel_size))
                    poly_mask = cv2.dilate(poly_mask, kernel)
                sterilized[poly_mask > 0] = fill_val
                logger.debug("Sterilized polygon artifact '%s' with %d px padding", category, pad)
                continue

            if box is not None and len(box) == 4:
                xmin, ymin, xmax, ymax = [int(round(v)) for v in box]
                x1 = max(0, xmin - pad)
                y1 = max(0, ymin - pad)
                x2 = min(img_w, xmax + pad)
                y2 = min(img_h, ymax + pad)

                if x2 > x1 and y2 > y1:
                    sterilized[y1:y2, x1:x2] = fill_val
                    logger.debug(
                        "Hard-blanked bounding box [%d, %d, %d, %d] for artifact '%s' with %d px padding",
                        x1, y1, x2, y2, category, pad
                    )

        return sterilized

    # -------------------------------------------------------------------------
    # METRIC COMPUTATION DELEGATION METHODS
    # -------------------------------------------------------------------------
    def compute_geometric_metrics(self, candidate_mask: np.ndarray) -> GeometricMetrics:
        """Compute deterministic geometric morphology metrics."""
        return compute_geometric_metrics(
            candidate_mask=candidate_mask,
            max_rectangularity_threshold=self.max_rectangularity,
            min_solidity_threshold=self.min_solidity,
            dp_epsilon_ratio=self.dp_epsilon_ratio,
            orthogonal_angle_range=(self.orthogonal_angle_min, self.orthogonal_angle_max)
        )

    def compute_spectral_metrics(
        self,
        candidate_patch: np.ndarray,
        candidate_mask: np.ndarray,
        is_rgb: bool = False
    ) -> SpectralMetrics:
        """Analyze colorimetric and HSV saturation profile."""
        return compute_spectral_metrics(
            candidate_patch=candidate_patch,
            candidate_mask=candidate_mask,
            high_sat_threshold=self.high_sat_threshold,
            max_color_swatch_ratio=self.max_color_swatch_ratio,
            is_rgb=is_rgb
        )

    def compute_texture_metrics(
        self,
        candidate_patch: np.ndarray,
        candidate_mask: np.ndarray,
        is_rgb: bool = False
    ) -> TextureMetrics:
        """Quantify interior gradient variance, typographic stroke density, and paper gating."""
        return compute_texture_metrics(
            candidate_patch=candidate_patch,
            candidate_mask=candidate_mask,
            paper_mean_val_threshold=self.paper_mean_val_threshold,
            paper_max_sat_threshold=self.paper_max_sat_threshold,
            laplacian_text_var_threshold=self.laplacian_text_var_threshold,
            canny_text_edge_threshold=self.canny_text_edge_threshold,
            is_rgb=is_rgb
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
        candidate_class: str = "basal_leaf_blade",
        is_rgb: bool = False
    ) -> FilterResult:
        """
        Execute multi-stage deterministic gatekeeper against an extracted leaf candidate.
        """
        geo_metrics = self.compute_geometric_metrics(candidate_mask)
        spec_metrics = self.compute_spectral_metrics(candidate_patch, candidate_mask, is_rgb=is_rgb)
        tex_metrics = self.compute_texture_metrics(candidate_patch, candidate_mask, is_rgb=is_rgb)

        rejection_reason = None
        reclassified_category = None
        routed_path = None

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

        is_valid = (rejection_reason is None)

        if is_valid:
            status = "VALID_LEAF"
        else:
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
        """
        base_name = f"{catalog_number}_patch{patch_id}_annotation"
        img_out_path = self.annotations_archive_dir / f"{base_name}.png"
        meta_out_path = self.annotations_archive_dir / f"{base_name}_meta.json"

        if candidate_patch is not None and candidate_patch.size > 0:
            if is_rgb and candidate_patch.ndim == 3:
                save_patch = cv2.cvtColor(candidate_patch, cv2.COLOR_RGB2BGR)
            else:
                save_patch = candidate_patch
            cv2.imwrite(str(img_out_path), save_patch)

        with open(meta_out_path, "w", encoding="utf-8") as f:
            json.dump({
                "catalog_number": catalog_number,
                "patch_id": patch_id,
                "archived_image": str(img_out_path),
                "metrics": metrics
            }, f, indent=2)

        logger.debug("Archived rejected annotation patch to '%s'", img_out_path)
        return str(img_out_path.resolve())
