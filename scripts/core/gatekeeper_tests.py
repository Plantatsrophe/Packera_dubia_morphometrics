
import os
import sys
import logging
import math
import numpy as np
import cv2
import random
import unittest
from typing import Dict, List, Tuple, Optional, Any, Union
from collections import defaultdict
from pathlib import Path
from dataclasses import dataclass

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Common imports
from scripts.core.config import CLASS_NAMES, CLASS_MAP, CLASS_COLORS_BGR, DEFAULT_WORKSPACE
from scripts.core.logger import setup_logging
from scripts.core.data_structures import ArtifactDetection, GeometricMetrics, SpectralMetrics, TextureMetrics, FilterResult, InstanceAnnotation
from scripts.core.tiling_utils import HerbariumAnnotation


from scripts.core.gatekeeper_engine import ArtifactFilterGatekeeper

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


if __name__ == '__main__':
    run_synthetic_test_suite()
