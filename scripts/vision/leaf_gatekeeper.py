#!/usr/bin/env python3
"""
===============================================================================
Module: leaf_gatekeeper.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Shape descriptor and botanical dissection utilities for Step 02:
    1. Classical dimensionless shape descriptors (Solidity, Circularity, Extent)
       delegating contour extraction and polygon morphology to LeafMachine3.
    2. Botanical dissection flagging (is_botanical_dissection) for lyrate-pinnatifid
       and bipinnatisect Packera silhouettes.
    3. Canonical routing tier normalization and adaptive per-plant minimum blade
       area threshold helpers.
===============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LM3_REPO_ROOT = PROJECT_ROOT / "LeafMachine3"
if LM3_REPO_ROOT.is_dir() and str(LM3_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(LM3_REPO_ROOT))

from leafmachine3.core.leaf_outline import extract_contour
from leafmachine3.core.morphometrics import polygon_morphology
from leafmachine3.core.orientation import length_width_from_box

__all__ = [
    "compute_classical_shape_descriptors",
    "is_botanical_dissection",
    "TIER_1_DIRECT",
    "TIER_2_REFLECTED",
    "CANONICAL_TIER_MAP",
    "normalize_assigned_tier",
    "compute_sheet_adaptive_min_blade_area_cm2",
    "LeafRoutingResult",
]


def compute_classical_shape_descriptors(
    mask: np.ndarray,
    bbox: Optional[Tuple[int, int, int, int]] = None,
    angle_cw: Optional[float] = None,
) -> Dict[str, float]:
    """
    Computes classical dimensionless shape descriptors (Krieger 2017; Russ 2011)
    by delegating contour extraction and polygon morphology to LeafMachine3:
      - Solidity = Area / Convex Hull Area
      - Circularity (Form Factor) = 4 * pi * Area / Perimeter^2
      - Extent = Mask Area / Oriented Bounding Box Area
    """
    if mask is None or mask.size == 0:
        return {
            "area_px": 0.0,
            "perimeter_px": 0.0,
            "convex_hull_area_px": 0.0,
            "bbox_area_px": 0.0,
            "solidity": 0.0,
            "circularity": 0.0,
            "extent": 0.0,
            "aspect_ratio": 0.0,
        }

    mask_u8 = (mask > 0).astype(np.uint8) * 255
    area_px = float(np.count_nonzero(mask_u8))
    if area_px <= 0.0:
        return {
            "area_px": 0.0,
            "perimeter_px": 0.0,
            "convex_hull_area_px": 0.0,
            "bbox_area_px": 0.0,
            "solidity": 0.0,
            "circularity": 0.0,
            "extent": 0.0,
            "aspect_ratio": 0.0,
        }

    bx, by, bw, bh = cv2.boundingRect(mask_u8)
    sub = mask_u8[by : by + bh, bx : bx + bw]
    main_cnt_xy = extract_contour(sub)
    contours, _ = cv2.findContours(sub, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if main_cnt_xy is None and not contours:
        return {
            "area_px": area_px,
            "perimeter_px": 0.0,
            "convex_hull_area_px": area_px,
            "bbox_area_px": float(max(bw * bh, 1)),
            "solidity": 1.0,
            "circularity": 0.0,
            "extent": min(1.0, area_px / float(max(bw * bh, 1))),
            "aspect_ratio": float(max(bw, bh)) / float(max(min(bw, bh), 1)),
        }

    main_cnt = (
        max(contours, key=cv2.contourArea)
        if contours
        else np.round(main_cnt_xy).astype(np.int32).reshape(-1, 1, 2)
    )
    if main_cnt_xy is None:
        main_cnt_xy = main_cnt.reshape(-1, 2).astype(np.float64)

    morph = polygon_morphology(main_cnt_xy, method="minarearect")
    morph_area = float(getattr(morph, "area_px", 0.0) or 0.0) if morph is not None else 0.0
    cnt_area = max(morph_area, float(cv2.contourArea(main_cnt)), area_px)

    perim_px = (
        float(getattr(morph, "perimeter_px", 0.0) or 0.0)
        if (morph is not None and float(getattr(morph, "perimeter_px", 0.0) or 0.0) > 0.0)
        else float(cv2.arcLength(main_cnt, True))
    )
    hull = cv2.convexHull(main_cnt)
    hull_area_cv = float(cv2.contourArea(hull))
    hull_area_morph = (
        float(getattr(morph, "convex_area_px", getattr(morph, "convex_hull_area", 0.0)) or 0.0)
        if morph is not None
        else 0.0
    )
    hull_area = hull_area_cv if hull_area_cv > 0.0 else hull_area_morph

    morph_sol = (
        float(getattr(morph, "solidity", getattr(morph, "convexity", 0.0)) or 0.0)
        if morph is not None
        else 0.0
    )
    solidity = (
        min(max(area_px / hull_area, 0.0), 1.0)
        if hull_area > 0.0
        else min(max(morph_sol, 0.0), 1.0)
    )

    morph_circ = float(getattr(morph, "circularity", 0.0) or 0.0) if morph is not None else 0.0
    circularity = (
        min(max((4.0 * math.pi * cnt_area) / (perim_px * perim_px), 0.0), 1.0)
        if perim_px > 0.0
        else min(max(morph_circ, 0.0), 1.0)
    )

    aspect_ratio = float(getattr(morph, "aspect_ratio", 0.0) or 0.0) if morph is not None else 0.0
    if angle_cw is not None and morph is not None and getattr(morph, "rotated_bbox", None) is not None:
        rbox = np.asarray(morph.rotated_bbox, dtype=np.float64)
        if rbox.shape == (4, 2):
            l_px, w_px = length_width_from_box(rbox, float(angle_cw))
            if l_px > 0.0 and w_px > 0.0:
                aspect_ratio = float(max(l_px, w_px)) / float(max(min(l_px, w_px), 1e-6))

    tight_axis_area = float(max(bw * bh, 1))
    rect_w, rect_h = cv2.minAreaRect(main_cnt)[1]
    oriented_area = min(tight_axis_area, float(max(rect_w * rect_h, 1.0)))

    if bbox is not None and len(bbox) == 4:
        ymin, xmin, ymax, xmax = [int(v) for v in bbox]
        box_h = max(0, ymax - ymin)
        box_w = max(0, xmax - xmin)
        ext_bbox_area = float(box_h * box_w) if (box_h > 0 and box_w > 0) else tight_axis_area
        inflation = max(1.0, ext_bbox_area / tight_axis_area)
        bbox_area = oriented_area * inflation
    else:
        bbox_area = oriented_area

    extent = min(max(area_px / max(bbox_area, 1.0), 0.0), 1.0)

    return {
        "area_px": area_px,
        "perimeter_px": perim_px,
        "convex_hull_area_px": hull_area,
        "bbox_area_px": bbox_area,
        "solidity": solidity,
        "circularity": circularity,
        "extent": extent,
        "aspect_ratio": aspect_ratio,
    }


def is_botanical_dissection(
    contour_or_mask: Optional[np.ndarray] = None,
    solidity: Optional[float] = None,
    blade_width: Optional[float] = None,
    p_apex: Optional[Tuple[float, float]] = None,
    p_base: Optional[Tuple[float, float]] = None,
    min_solidity: float = 0.18,
    max_solidity: float = 0.72,
    bipinnatisect_min_solidity: Optional[float] = None,
    min_defect_count: int = 3,
    defect_depth_ratio: float = 0.15,
    check_solidity_range: bool = True,
    *,
    contour: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None,
) -> bool:
    """
    Distinguishes true botanical dissection (lyrate-pinnatifid and bipinnatisect lobes
    in Packera taxa) from unilateral occlusion or physical damage.
    """
    active_input = contour if contour is not None else (contour_or_mask if contour_or_mask is not None else mask)
    if active_input is None or active_input.size == 0:
        return False
    contour = active_input

    if contour.ndim == 2 and contour.shape[1] > 2:
        bin_in = (contour > 0).astype(np.uint8) * 255
        cnts, _ = cv2.findContours(bin_in, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not cnts:
            return False
        cnt = max(cnts, key=cv2.contourArea).reshape(-1, 2).astype(np.int32)
    else:
        cnt = contour.reshape(-1, 2).astype(np.int32)

    if len(cnt) < 4:
        return False

    cnt_f = cnt.astype(np.float32)
    if solidity is None:
        area = float(cv2.contourArea(cnt_f))
        hull = cv2.convexHull(cnt_f)
        hull_area = float(cv2.contourArea(hull))
        solidity = float(area / hull_area) if hull_area > 0 else 0.0

    eff_min_solidity = (
        float(bipinnatisect_min_solidity)
        if bipinnatisect_min_solidity is not None
        else float(min_solidity)
    )
    if check_solidity_range:
        if not (eff_min_solidity <= solidity < max_solidity):
            return False

    rect = cv2.minAreaRect(cnt_f)
    minor_w = float(min(rect[1]))
    major_l = float(max(rect[1]))

    is_high_dpi = (
        (contour.ndim == 2 and max(contour.shape[:2]) > 500)
        or cnt_f[:, 0].max() > 500
        or cnt_f[:, 1].max() > 500
        or (solidity >= 0.48 and major_l > 240.0)
        or minor_w > 360.0
    )
    if is_high_dpi and major_l > 220.0 and solidity >= 0.55 and (minor_w / max(major_l, 1.0)) > 0.88:
        return False

    if blade_width is None or blade_width <= 0:
        blade_width = minor_w
    if blade_width <= 0:
        return False

    min_defect_depth = defect_depth_ratio * blade_width

    cnt_cv = cnt.reshape(-1, 1, 2)
    hull_indices = cv2.convexHull(cnt_cv, returnPoints=False)
    if hull_indices is None or len(hull_indices) < 3:
        return False

    try:
        defects = cv2.convexityDefects(cnt_cv, hull_indices)
    except cv2.error:
        try:
            hull_sorted = np.unique(hull_indices).reshape(-1, 1)
            if len(hull_sorted) < 3:
                return False
            defects = cv2.convexityDefects(cnt_cv, hull_sorted)
        except cv2.error:
            return False
    if defects is None or len(defects) == 0:
        return False

    def_arr = defects.reshape(-1, 4)
    depths = def_arr[:, 3].astype(np.float64) / 256.0
    farthest_indices = def_arr[:, 2].astype(np.int32)

    is_native_scan = (
        (contour.ndim == 2 and max(contour.shape[:2]) > 500)
        or cnt_f[:, 0].max() > 500
        or cnt_f[:, 1].max() > 500
        or major_l > 350.0
    )
    sig_mask = depths >= min_defect_depth
    if solidity < 0.48:
        eff_min_defect_count = max(min_defect_count, 4)
    else:
        eff_min_defect_count = (
            min(min_defect_count, 2)
            if (is_native_scan and solidity >= 0.48 and np.any(sig_mask) and float(np.max(depths[sig_mask])) >= 0.11 * blade_width)
            else min_defect_count
        )
    if np.sum(sig_mask) < eff_min_defect_count:
        return False

    sig_farthest_pts = cnt[farthest_indices[sig_mask]]

    if p_apex is None or p_base is None or solidity < 0.48:
        line_params = cv2.fitLine(cnt_f, cv2.DIST_L2, 0, 0.01, 0.01)
        vx, vy, x0, y0 = [float(v[0]) for v in line_params]
        projections = (cnt_f[:, 0] - x0) * vx + (cnt_f[:, 1] - y0) * vy
        t_min = float(np.min(projections))
        t_max = float(np.max(projections))
        p_apex = (x0 + t_min * vx, y0 + t_min * vy)
        p_base = (x0 + t_max * vx, y0 + t_max * vy)

    dx = float(p_base[0] - p_apex[0])
    dy = float(p_base[1] - p_apex[1])
    line_len = math.hypot(dx, dy)
    if line_len < 1e-6:
        return False

    cross_vals = dx * (sig_farthest_pts[:, 1] - p_apex[1]) - dy * (sig_farthest_pts[:, 0] - p_apex[0])
    tol = 0.01 * line_len * blade_width

    side_a_count = int(np.sum(cross_vals > tol))
    side_b_count = int(np.sum(cross_vals < -tol))
    total_sig = side_a_count + side_b_count

    is_dissected = (total_sig >= eff_min_defect_count) and (side_a_count >= 1) and (side_b_count >= 1)
    return bool(is_dissected)


TIER_1_DIRECT = "Tier_1_Direct"
TIER_2_REFLECTED = "Tier_2_Reflected"

CANONICAL_TIER_MAP: Dict[str, str] = {
    "tier1": TIER_1_DIRECT,
    "tier_1": TIER_1_DIRECT,
    "tier 1": TIER_1_DIRECT,
    "tier 1 (dissected)": TIER_1_DIRECT,
    "tier 1 (bipinnatisect)": TIER_1_DIRECT,
    "tier_1_direct": TIER_1_DIRECT,
    "Tier_1_Direct": TIER_1_DIRECT,
    "tier2": TIER_2_REFLECTED,
    "tier_2": TIER_2_REFLECTED,
    "tier 2": TIER_2_REFLECTED,
    "tier_2_reflected": TIER_2_REFLECTED,
    "Tier_2_Reflected": TIER_2_REFLECTED,
}


def normalize_assigned_tier(tier_value: str) -> str:
    """
    Normalizes routing tier values ('tier1', 'tier2', 'Tier 1', 'Tier 2')
    to the canonical cross-language CSV contract names ('Tier_1_Direct', 'Tier_2_Reflected').
    """
    if not tier_value:
        return "rejected"
    clean = str(tier_value).strip()
    return CANONICAL_TIER_MAP.get(clean, CANONICAL_TIER_MAP.get(clean.lower(), clean))


def compute_sheet_adaptive_min_blade_area_cm2(
    candidate_areas_cm2: Any,
    config_min_blade_area_cm2: float = 0.75,
    absolute_floor_cm2: float = 0.15,
) -> float:
    """
    Computes a cohort/sheet-adaptive physical minimum blade area threshold (cm^2)
    for calibrated herbarium sheets:
        effective_min_area_cm2 = min(config_min_blade_area_cm2, max(0.15, 0.25 * median_sheet_leaf_area_cm2))
    """
    cfg_floor = float(config_min_blade_area_cm2)
    abs_floor = float(absolute_floor_cm2)
    if candidate_areas_cm2 is None:
        return cfg_floor
    valid_areas: List[float] = []
    for val in candidate_areas_cm2:
        if val is None:
            continue
        try:
            fval = float(val)
        except (TypeError, ValueError):
            continue
        if np.isfinite(fval) and fval > 0.0:
            valid_areas.append(fval)
    if not valid_areas:
        return cfg_floor
    median_sheet_leaf_area_cm2 = float(np.median(valid_areas))
    return float(min(cfg_floor, max(abs_floor, 0.25 * median_sheet_leaf_area_cm2)))


@dataclass
class LeafRoutingResult:
    """Structured diagnostic and extraction result container."""
    tier: str
    assigned_tier: str
    solidity: float
    ucs: float
    is_folded: bool
    is_irregular_fold: bool
    is_dissected: bool
    processed_mask: np.ndarray
    contour_coords: Optional[np.ndarray] = None
    rejection_reason: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def canonical_tier(self) -> str:
        """Returns the canonical cross-language CSV tier ('Tier_1_Direct', 'Tier_2_Reflected', or 'rejected')."""
        return normalize_assigned_tier(self.assigned_tier)

    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        if key in self.metadata:
            return self.metadata[key]
        raise KeyError(key)

    def get(self, key: str, default: Any = None) -> Any:
        try:
            return self[key]
        except KeyError:
            return default
