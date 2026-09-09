#!/usr/bin/env python3
"""
===============================================================================
Script: tune_geometry_parameters.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Automated empirical calibration tool for geometric gatekeeper thresholds.
    Directly ingests ground-truth SAM 2 polygon annotations (COCO format) to
    derive data-driven cutoffs for:
      1. Straight-chord leaf fold detection (Aspect Ratio W/L, Chord Deviation)
      2. Lyrate / pinnatifid botanical dissection (Solidity, Sinus Defect Depth)

    Generates a 4-panel publication-grade diagnostic PDF and optionally updates
    thresholds in config/config.yaml while preserving all comments, taxa lists,
    and model weights.
===============================================================================
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.logger import setup_logging

logger = setup_logging(name="TuneGeometry")


# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class LeafGeometrySample:
    """Morphometric measurements extracted from a single annotated leaf contour."""
    annotation_id: int
    image_id: int
    contour_area: float
    aspect_ratio: float           # min(W, H) / max(W, H)
    solidity: float               # contourArea / convexHullArea
    chord_deviation_ratio: float  # flatter margin perpendicular distance / chord length
    sinus_depth_ratio: float      # max defect depth / bounding blade width
    chord_length: float
    blade_width: float
    is_flat_margin_candidate: bool = False
    is_lobed_candidate: bool = False


@dataclass
class MetricSummary:
    """Descriptive statistics and percentiles for a geometric metric."""
    name: str
    n_samples: int
    min_val: float
    max_val: float
    mean_val: float
    std_val: float
    percentiles: Dict[int, float] = field(default_factory=dict)

    def format_row(self) -> str:
        p = self.percentiles
        return (
            f"{self.name:<24} | {self.n_samples:>5} | {self.mean_val:>6.4f} | "
            f"{p.get(5, 0.0):>6.4f} | {p.get(10, 0.0):>6.4f} | {p.get(25, 0.0):>6.4f} | "
            f"{p.get(50, 0.0):>6.4f} | {p.get(75, 0.0):>6.4f} | {p.get(90, 0.0):>6.4f} | "
            f"{p.get(95, 0.0):>6.4f}"
        )


@dataclass
class RecommendedThresholds:
    """Empirically derived botanical gatekeeper thresholds."""
    max_folded_aspect_ratio: float
    max_chord_deviation_ratio: float
    min_solidity_dissected: float
    sinus_defect_min_depth_ratio: float

    # Raw empirical reference points
    raw_wl_p10: float
    raw_chord_flat_p95: float
    raw_solidity_lobed_p5: float
    raw_sinus_p25: float


# =============================================================================
# COCO Geometry Extractor
# =============================================================================

class COCOGeometryExtractor:
    """Ingests COCO polygon annotations and calculates empirical geometric metrics."""

    def __init__(
        self,
        coco_path: Path,
        min_contour_area: float = 150.0,
        target_category_name: Optional[str] = None,
        target_category_id: Optional[int] = None,
    ) -> None:
        self.coco_path = Path(coco_path)
        self.min_contour_area = min_contour_area
        self.target_category_name = target_category_name
        self.target_category_id = target_category_id

        if not self.coco_path.exists():
            raise FileNotFoundError(f"COCO annotation file not found: {self.coco_path}")

        with open(self.coco_path, "r", encoding="utf-8") as f:
            self.coco_data = json.load(f)

        self.categories = self.coco_data.get("categories", [])
        self.annotations = self.coco_data.get("annotations", [])
        self.images = {img["id"]: img for img in self.coco_data.get("images", [])}

        self.resolved_category_ids = self._resolve_target_categories()

    def _resolve_target_categories(self) -> List[int]:
        """Resolves target category IDs for basal leaf blade instances."""
        if self.target_category_id is not None:
            logger.info(f"Targeting explicit category ID: {self.target_category_id}")
            return [self.target_category_id]

        target_names = {
            "basal_leaf_blade",
            "ideal_leaf",
            "basal_leaf_whole",
            "basal_leaf",
        }
        if self.target_category_name:
            target_names = {self.target_category_name.strip().lower()}

        matching_ids: List[int] = []
        for cat in self.categories:
            c_name = cat.get("name", "").strip().lower()
            c_id = cat.get("id")
            if c_name in target_names or (c_id == 0 and "basal" in c_name):
                matching_ids.append(c_id)

        # Fallback logic if categories use 0-indexing or default names
        if not matching_ids:
            cat_ids = [c.get("id") for c in self.categories]
            if 0 in cat_ids:
                matching_ids = [0]
            elif 1 in cat_ids:
                matching_ids = [1]
            else:
                matching_ids = cat_ids[:1] if cat_ids else [0]

        cat_desc = [
            f"{c.get('name')} (ID: {c.get('id')})"
            for c in self.categories if c.get("id") in matching_ids
        ]
        logger.info(f"Resolved target leaf categories: {', '.join(cat_desc)}")
        return matching_ids

    def extract_samples(self) -> List[LeafGeometrySample]:
        """Iterates over annotations and computes the four geometric metrics."""
        samples: List[LeafGeometrySample] = []
        total_instances = 0
        filtered_noise = 0

        for ann in self.annotations:
            cat_id = ann.get("category_id")
            if cat_id not in self.resolved_category_ids:
                continue

            ann_id = ann.get("id", -1)
            img_id = ann.get("image_id", -1)
            segs = ann.get("segmentation", [])

            if not isinstance(segs, list):
                continue

            for seg in segs:
                total_instances += 1
                pts = np.array(seg, dtype=np.float32).reshape(-1, 2)
                if len(pts) < 4:
                    filtered_noise += 1
                    continue

                cnt = pts.astype(np.int32)
                area = float(cv2.contourArea(cnt))
                if area < self.min_contour_area:
                    filtered_noise += 1
                    continue

                sample = self._compute_sample_metrics(cnt, area, ann_id, img_id)
                if sample is not None:
                    samples.append(sample)

        logger.info(
            f"Extracted {len(samples)} valid leaf blade samples from {total_instances} "
            f"polygon instances ({filtered_noise} filtered noise fragments < {self.min_contour_area} px)."
        )
        return samples

    def _compute_sample_metrics(
        self,
        cnt: np.ndarray,
        area: float,
        ann_id: int,
        img_id: int
    ) -> Optional[LeafGeometrySample]:
        """Calculates W/L aspect ratio, solidity, chord deviation, and sinus depth."""
        cnt_f = cnt.reshape(-1, 2).astype(np.float32)

        # 1. Aspect Ratio (W/L) via minimum area bounding box
        rect = cv2.minAreaRect(cnt_f)
        w, h = float(rect[1][0]), float(rect[1][1])
        dim_min, dim_max = min(w, h), max(w, h)
        if dim_max < 1e-4:
            return None
        aspect_ratio = dim_min / dim_max
        blade_width = dim_min

        # 2. Solidity: contourArea / convexHullArea
        hull = cv2.convexHull(cnt)
        hull_area = float(cv2.contourArea(hull))
        solidity = (area / hull_area) if hull_area > 1e-4 else 1.0

        # 3. Chord Deviation: longest longitudinal chord from apex to base
        line_params = cv2.fitLine(cnt_f, cv2.DIST_L2, 0, 0.01, 0.01)
        vx, vy, x0, y0 = [float(v[0]) for v in line_params]
        projections = (cnt_f[:, 0] - x0) * vx + (cnt_f[:, 1] - y0) * vy
        min_idx = int(np.argmin(projections))
        max_idx = int(np.argmax(projections))
        p_apex = cnt_f[min_idx]
        p_base = cnt_f[max_idx]

        chord_dx = float(p_base[0] - p_apex[0])
        chord_dy = float(p_base[1] - p_apex[1])
        chord_len = math.hypot(chord_dx, chord_dy)

        chord_dev_ratio = 0.0
        if chord_len >= 5.0 and min_idx != max_idx:
            idx_1 = min(min_idx, max_idx)
            idx_2 = max(min_idx, max_idx)
            pts_arr = cnt.reshape(-1, 2)
            side_a_pts = pts_arr[idx_1 : idx_2 + 1]
            side_b_pts = np.vstack([pts_arr[idx_2:], pts_arr[: idx_1 + 1]])

            dev_a = self._calc_flatter_margin_deviation(side_a_pts)
            dev_b = self._calc_flatter_margin_deviation(side_b_pts)
            chord_dev_ratio = min(dev_a, dev_b)
        else:
            chord_dev_ratio = 0.05

        # 4. Sinus Depth Ratio via cv2.convexityDefects
        sinus_depth_ratio = 0.0
        cnt_cv = cnt.reshape(-1, 1, 2)
        hull_indices = cv2.convexHull(cnt_cv, returnPoints=False)
        if hull_indices is not None and len(hull_indices) >= 3 and blade_width > 0:
            defects = cv2.convexityDefects(cnt_cv, hull_indices)
            if defects is not None and len(defects) > 0:
                depths = defects.reshape(-1, 4)[:, 3].astype(np.float64) / 256.0
                max_defect_depth = float(np.max(depths))
                sinus_depth_ratio = max_defect_depth / blade_width

        is_flat = chord_dev_ratio <= 0.055
        is_lobed = solidity < 0.75 or sinus_depth_ratio >= 0.08

        return LeafGeometrySample(
            annotation_id=ann_id,
            image_id=img_id,
            contour_area=area,
            aspect_ratio=aspect_ratio,
            solidity=solidity,
            chord_deviation_ratio=chord_dev_ratio,
            sinus_depth_ratio=sinus_depth_ratio,
            chord_length=chord_len,
            blade_width=blade_width,
            is_flat_margin_candidate=is_flat,
            is_lobed_candidate=is_lobed,
        )

    @staticmethod
    def _calc_flatter_margin_deviation(side_pts: np.ndarray) -> float:
        """Calculates maximum perpendicular distance to chord line divided by chord length."""
        if len(side_pts) < 3:
            return 1.0
        p1 = side_pts[0]
        p2 = side_pts[-1]
        length = float(math.hypot(p2[0] - p1[0], p2[1] - p1[1]))
        if length < 1e-6:
            return 1.0

        dx = float(p2[0] - p1[0])
        dy = float(p2[1] - p1[1])
        # Vectorized perpendicular distance: |dy*x - dx*y + x2*y1 - y2*x1| / L
        cross = np.abs(dy * side_pts[:, 0] - dx * side_pts[:, 1] + p2[0] * p1[1] - p2[1] * p1[0])
        max_dist = float(np.max(cross)) / length
        return max_dist / length


# =============================================================================
# Statistical Analysis & Threshold Formulation
# =============================================================================

def compute_metric_summary(name: str, values: List[float]) -> MetricSummary:
    """Computes descriptive statistical percentiles (5th, 10th, 25th, 50th, 75th, 90th, 95th)."""
    if not values:
        return MetricSummary(name, 0, 0.0, 0.0, 0.0, 0.0)

    arr = np.array(values, dtype=np.float64)
    pcts = [5, 10, 25, 50, 75, 90, 95]
    p_dict = {p: float(np.percentile(arr, p)) for p in pcts}

    return MetricSummary(
        name=name,
        n_samples=len(arr),
        min_val=float(np.min(arr)),
        max_val=float(np.max(arr)),
        mean_val=float(np.mean(arr)),
        std_val=float(np.std(arr)),
        percentiles=p_dict,
    )


def formulate_recommended_thresholds(samples: List[LeafGeometrySample]) -> RecommendedThresholds:
    """
    Formulates recommended geometric thresholds based on botanical morphology:
      1. max_folded_aspect_ratio: 10th percentile of W/L (botanical range ~0.40–0.44).
         For unfolded basal blades, W/L is typically 0.60–0.85; a folded hemi-blade
         has half width, leading to W/L <= 0.40–0.44.
      2. max_chord_deviation_ratio: 95th percentile of flat-margin chords (~0.035–0.045).
      3. min_solidity_dissected: 5th percentile of solidity for lobed specimens (~0.48–0.52).
      4. sinus_defect_min_depth_ratio: 25th percentile of significant sinuses (~0.08–0.10).
    """
    wl_all = [s.aspect_ratio for s in samples]
    sol_all = [s.solidity for s in samples]
    chord_all = [s.chord_deviation_ratio for s in samples]
    sinus_all = [s.sinus_depth_ratio for s in samples if s.sinus_depth_ratio > 0.0]

    # Flat-margin chords (chords where flatter margin shows straight-line tendency)
    flat_chords = [s.chord_deviation_ratio for s in samples if s.is_flat_margin_candidate]
    if len(flat_chords) < 3:
        cutoff = float(np.percentile(chord_all, 25)) if chord_all else 0.05
        flat_chords = [c for c in chord_all if c <= max(cutoff, 0.055)]

    # Lobed / dissected specimens (solidity < 0.75 or lobed candidate)
    lobed_solidity = [s.solidity for s in samples if s.is_lobed_candidate]
    if len(lobed_solidity) < 3:
        cutoff_sol = float(np.percentile(sol_all, 35)) if sol_all else 0.72
        lobed_solidity = [sol for sol in sol_all if sol <= cutoff_sol]

    # Raw empirical percentiles
    raw_wl_p10 = float(np.percentile(wl_all, 10)) if wl_all else 0.44
    raw_chord_flat_p95 = float(np.percentile(flat_chords, 95)) if flat_chords else 0.035
    raw_solidity_lobed_p5 = float(np.percentile(lobed_solidity, 5)) if lobed_solidity else 0.50
    raw_sinus_p25 = float(np.percentile(sinus_all, 25)) if sinus_all else 0.08

    # Botanical calibration & guardrails:
    # 1. Folded Aspect Ratio:
    #    When dataset contains full rosettes (some naturally elongated), standard botanical
    #    calibration for Packera basal leaves places hemi-blade fold cutoff at ~0.40–0.44.
    if 0.38 <= raw_wl_p10 <= 0.46:
        rec_wl = round(raw_wl_p10, 3)
    elif raw_wl_p10 < 0.38:
        rec_wl = 0.42
    else:
        rec_wl = min(0.45, round(raw_wl_p10, 3))

    # 2. Chord Deviation:
    #    Bounded within [0.030, 0.048], typically ~0.035–0.045
    rec_chord = max(0.030, min(0.048, round(raw_chord_flat_p95, 4)))

    # 3. Min Dissected Solidity:
    #    Bounded within [0.48, 0.54], typically ~0.48–0.52
    rec_sol = max(0.48, min(0.54, round(raw_solidity_lobed_p5, 2)))

    # 4. Sinus Defect Depth:
    #    Bounded within [0.07, 0.12], typically ~0.08–0.10
    rec_sinus = max(0.07, min(0.12, round(raw_sinus_p25, 3)))

    return RecommendedThresholds(
        max_folded_aspect_ratio=rec_wl,
        max_chord_deviation_ratio=rec_chord,
        min_solidity_dissected=rec_sol,
        sinus_defect_min_depth_ratio=rec_sinus,
        raw_wl_p10=raw_wl_p10,
        raw_chord_flat_p95=raw_chord_flat_p95,
        raw_solidity_lobed_p5=raw_solidity_lobed_p5,
        raw_sinus_p25=raw_sinus_p25,
    )


# =============================================================================
# Diagnostic Visualization Plot (4-Panel PDF)
# =============================================================================

def export_diagnostic_plots(
    samples: List[LeafGeometrySample],
    recs: RecommendedThresholds,
    output_pdf: Path
) -> None:
    """
    Generates a 4-panel publication-grade diagnostic PDF with Agg backend:
      - Panel A: Histogram of W/L with recommended fold cutoff line.
      - Panel B: Histogram of Solidity with default (0.72) and dissected (0.50) cutoff lines.
      - Panel C: Density plot of chord deviation ratios.
      - Panel D: Distribution of sinus defect depths.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_pdf = Path(output_pdf)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)

    wl_vals = [s.aspect_ratio for s in samples]
    sol_vals = [s.solidity for s in samples]
    chord_vals = [s.chord_deviation_ratio for s in samples]
    sinus_vals = [s.sinus_depth_ratio for s in samples]

    fig, axs = plt.subplots(2, 2, figsize=(12, 10), dpi=300)
    fig.patch.set_facecolor("white")

    color_primary = "#2b5c8f"    # Deep botanical blue
    color_accent = "#2a7f62"     # Foliar green
    color_alert = "#c0392b"      # Warning red
    color_gold = "#d35400"       # Ochre / gold

    # -------------------------------------------------------------------------
    # Panel A: Aspect Ratio (W/L) Histogram
    # -------------------------------------------------------------------------
    ax_a = axs[0, 0]
    ax_a.hist(wl_vals, bins=20, color=color_primary, edgecolor="white", alpha=0.85, density=False)
    ax_a.axvline(
        recs.max_folded_aspect_ratio,
        color=color_alert,
        linestyle="--",
        linewidth=2.0,
        label=f"Recommended Cutoff ({recs.max_folded_aspect_ratio:.2f})",
    )
    ax_a.set_title("A. Aspect Ratio (W/L) Distribution", fontsize=12, fontweight="bold", pad=8)
    ax_a.set_xlabel("Width / Length Ratio ($W/L$)", fontsize=10)
    ax_a.set_ylabel("Annotated Leaf Count", fontsize=10)
    ax_a.grid(True, linestyle=":", alpha=0.6)
    ax_a.legend(loc="upper right", frameon=True, fontsize=9)

    # -------------------------------------------------------------------------
    # Panel B: Solidity Histogram
    # -------------------------------------------------------------------------
    ax_b = axs[0, 1]
    ax_b.hist(sol_vals, bins=20, color=color_accent, edgecolor="white", alpha=0.85, density=False)
    ax_b.axvline(
        0.72,
        color="#2c3e50",
        linestyle="-.",
        linewidth=1.8,
        label="Default Tier 1 Cutoff (0.72)",
    )
    ax_b.axvline(
        recs.min_solidity_dissected,
        color=color_alert,
        linestyle="--",
        linewidth=2.0,
        label=f"Dissected Lower Bound ({recs.min_solidity_dissected:.2f})",
    )
    ax_b.set_title("B. Leaf Solidity Distribution", fontsize=12, fontweight="bold", pad=8)
    ax_b.set_xlabel("Solidity ($Area / ConvexHullArea$)", fontsize=10)
    ax_b.set_ylabel("Annotated Leaf Count", fontsize=10)
    ax_b.grid(True, linestyle=":", alpha=0.6)
    ax_b.legend(loc="upper left", frameon=True, fontsize=9)

    # -------------------------------------------------------------------------
    # Panel C: Chord Deviation Ratio Density
    # -------------------------------------------------------------------------
    ax_c = axs[1, 0]
    ax_c.hist(
        chord_vals,
        bins=22,
        color="#34495e",
        edgecolor="white",
        alpha=0.75,
        density=True,
    )
    try:
        from scipy.stats import gaussian_kde
        kde = gaussian_kde(chord_vals)
        x_grid = np.linspace(min(chord_vals), max(chord_vals), 200)
        ax_c.plot(x_grid, kde(x_grid), color="#1abc9c", linewidth=2.0, label="Empirical Density (KDE)")
    except Exception:
        pass

    ax_c.axvline(
        recs.max_chord_deviation_ratio,
        color=color_alert,
        linestyle="--",
        linewidth=2.0,
        label=f"Straight-Chord Max Dev ({recs.max_chord_deviation_ratio:.3f})",
    )
    ax_c.set_title("C. Lateral Margin Chord Deviation Ratio", fontsize=12, fontweight="bold", pad=8)
    ax_c.set_xlabel("Max Lateral Perpendicular Deviation / Chord Length", fontsize=10)
    ax_c.set_ylabel("Probability Density", fontsize=10)
    ax_c.grid(True, linestyle=":", alpha=0.6)
    ax_c.legend(loc="upper right", frameon=True, fontsize=9)

    # -------------------------------------------------------------------------
    # Panel D: Sinus Defect Depth Ratio
    # -------------------------------------------------------------------------
    ax_d = axs[1, 1]
    ax_d.hist(sinus_vals, bins=20, color=color_gold, edgecolor="white", alpha=0.85, density=False)
    ax_d.axvline(
        recs.sinus_defect_min_depth_ratio,
        color=color_alert,
        linestyle="--",
        linewidth=2.0,
        label=f"Min Significant Sinus ({recs.sinus_defect_min_depth_ratio:.3f})",
    )
    ax_d.set_title("D. Convexity Defect Sinus Depth Ratio", fontsize=12, fontweight="bold", pad=8)
    ax_d.set_xlabel("Max Defect Depth / Bounding Blade Width", fontsize=10)
    ax_d.set_ylabel("Annotated Leaf Count", fontsize=10)
    ax_d.grid(True, linestyle=":", alpha=0.6)
    ax_d.legend(loc="upper right", frameon=True, fontsize=9)

    fig.suptitle(
        "Empirical Geometry Threshold Distributions (SAM 2 Botanical Calibration)\n"
        "Packera dubia Species Delimitation Pipeline — NCU Herbarium",
        fontsize=14,
        fontweight="bold",
        y=0.98,
    )
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_pdf, format="pdf", dpi=300)
    plt.close(fig)

    logger.info(f"Exported 4-panel diagnostic plot to: {output_pdf}")


# =============================================================================
# Automated Configuration Updater
# =============================================================================

def update_pipeline_config(
    config_path: Path,
    recs: RecommendedThresholds,
) -> bool:
    """
    Safely updates only thresholds.fold_detection and thresholds.solidity in config.yaml,
    strictly preserving existing comments, indentation, taxa lists, and paths.
    """
    config_path = Path(config_path)
    if not config_path.exists():
        logger.error(f"Target config file does not exist: {config_path}")
        return False

    with open(config_path, "r", encoding="utf-8") as f:
        content = f.read()

    pattern_chord = r"(max_chord_deviation_ratio:\s*)[0-9\.]+"
    pattern_aspect = r"(max_folded_aspect_ratio:\s*)[0-9\.]+"
    pattern_min_diss = r"(min_solidity_dissected:\s*|min_dissected:\s*)[0-9\.]+"
    pattern_sinus = r"(sinus_defect_min_depth_ratio:\s*)[0-9\.]+"

    new_content = content

    if re.search(pattern_chord, new_content):
        new_content = re.sub(pattern_chord, rf"\g<1>{recs.max_chord_deviation_ratio:.3f}", new_content)
    if re.search(pattern_aspect, new_content):
        new_content = re.sub(pattern_aspect, rf"\g<1>{recs.max_folded_aspect_ratio:.2f}", new_content)
    if re.search(pattern_min_diss, new_content):
        new_content = re.sub(pattern_min_diss, rf"\g<1>{recs.min_solidity_dissected:.2f}", new_content)

    if re.search(pattern_sinus, new_content):
        new_content = re.sub(pattern_sinus, rf"\g<1>{recs.sinus_defect_min_depth_ratio:.3f}", new_content)
    else:
        diss_match = re.search(r"([ ]+)(min_solidity_dissected:\s*[0-9\.]+|min_dissected:\s*[0-9\.]+)", new_content)
        if diss_match:
            indent = diss_match.group(1)
            line = diss_match.group(0)
            replacement = f"{line}\n{indent}sinus_defect_min_depth_ratio: {recs.sinus_defect_min_depth_ratio:.3f}"
            new_content = new_content.replace(line, replacement, 1)

    import yaml
    try:
        parsed = yaml.safe_load(new_content)
        thresh = parsed.get("thresholds", {})
        assert "fold_detection" in thresh, "Missing fold_detection in parsed config"
        assert "dissection" in thresh or "solidity" in thresh, "Missing dissection or solidity in parsed config"
    except Exception as exc:
        logger.error(f"Configuration safety check failed: YAML invalid after substitution: {exc}")
        return False

    with open(config_path, "w", encoding="utf-8") as f:
        f.write(new_content)

    logger.info(f"Successfully updated geometry thresholds in: {config_path}")
    logger.info(
        f"  fold_detection.max_chord_deviation_ratio -> {recs.max_chord_deviation_ratio:.3f}\n"
        f"  fold_detection.max_folded_aspect_ratio  -> {recs.max_folded_aspect_ratio:.2f}\n"
        f"  solidity.min_dissected                  -> {recs.min_solidity_dissected:.2f}\n"
        f"  solidity.sinus_defect_min_depth_ratio   -> {recs.sinus_defect_min_depth_ratio:.3f}"
    )
    return True


# =============================================================================
# CLI Orchestrator & Entry Point
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    """Constructs command-line argument parser for geometry tuning."""
    parser = argparse.ArgumentParser(
        description="Calibrate geometric fold detection and lyrate dissection thresholds from COCO annotations."
    )
    parser.add_argument(
        "--annotations",
        type=Path,
        default=PROJECT_ROOT / "data" / "annotations" / "packera_train_coco.json",
        help="Path to ground-truth COCO annotations JSON (default: data/annotations/packera_train_coco.json)",
    )
    parser.add_argument(
        "--output-plot",
        type=Path,
        default=PROJECT_ROOT / "outputs" / "figures" / "geometry_parameter_distributions.pdf",
        help="Destination path for 4-panel diagnostic distribution PDF.",
    )
    parser.add_argument(
        "--config-path",
        type=Path,
        default=PROJECT_ROOT / "config" / "config.yaml",
        help="Path to pipeline configuration YAML to calibrate.",
    )
    parser.add_argument(
        "--min-area",
        type=float,
        default=150.0,
        help="Minimum contour area threshold in pixels to filter noise (default: 150.0).",
    )
    parser.add_argument(
        "--category-id",
        type=int,
        default=None,
        help="Specific COCO category ID to target (default: auto-detect Class 0 / ideal_leaf).",
    )
    parser.add_argument(
        "--category-name",
        type=str,
        default=None,
        help="Specific COCO category name to target (e.g. basal_leaf_blade, ideal_leaf).",
    )
    parser.add_argument(
        "--update-config",
        action="store_true",
        default=False,
        help="Update thresholds.fold_detection and thresholds.solidity in config.yaml.",
    )
    return parser


def run_calibration(
    annotations_path: Path,
    output_plot: Path,
    config_path: Optional[Path] = None,
    update_config: bool = False,
    min_area: float = 150.0,
    category_id: Optional[int] = None,
    category_name: Optional[str] = None,
) -> RecommendedThresholds:
    """Executes the complete calibration pipeline."""
    logger.info("=== Starting Geometric Threshold Calibration ===")
    logger.info(f"Ingesting COCO annotations from: {annotations_path}")

    extractor = COCOGeometryExtractor(
        coco_path=annotations_path,
        min_contour_area=min_area,
        target_category_name=category_name,
        target_category_id=category_id,
    )
    samples = extractor.extract_samples()

    if not samples:
        logger.error("No valid leaf blade instances found meeting criteria.")
        sys.exit(1)

    s_wl = compute_metric_summary("Aspect Ratio (W/L)", [s.aspect_ratio for s in samples])
    s_sol = compute_metric_summary("Solidity", [s.solidity for s in samples])
    s_chord = compute_metric_summary("Chord Deviation", [s.chord_deviation_ratio for s in samples])
    s_sinus = compute_metric_summary("Sinus Depth Ratio", [s.sinus_depth_ratio for s in samples])

    recs = formulate_recommended_thresholds(samples)

    header = f"{'Metric':<24} | {'N':>5} | {'Mean':>6} | {'P5':>6} | {'P10':>6} | {'P25':>6} | {'P50':>6} | {'P75':>6} | {'P90':>6} | {'P95':>6}"
    divider = "-" * len(header)
    print("\n" + "=" * 80)
    print("EMPIRICAL GEOMETRIC DISTRIBUTIONS (SAM 2 ANNOTATED LEAVES)")
    print("=" * 80)
    print(header)
    print(divider)
    print(s_wl.format_row())
    print(s_sol.format_row())
    print(s_chord.format_row())
    print(s_sinus.format_row())
    print(divider)

    print("\n" + "=" * 80)
    print("RECOMMENDED BOTANICAL GATEKEEPER THRESHOLDS")
    print("=" * 80)
    print(f"  * Fold Detection - Max Aspect Ratio (W/L)     : {recs.max_folded_aspect_ratio:.2f}  (Empirical P10: {recs.raw_wl_p10:.3f})")
    print(f"  * Fold Detection - Max Chord Deviation Ratio  : {recs.max_chord_deviation_ratio:.3f} (Flat Chords P95: {recs.raw_chord_flat_p95:.4f})")
    print(f"  * Dissection Check - Min Dissected Solidity   : {recs.min_solidity_dissected:.2f}  (Lobed Leaves P5: {recs.raw_solidity_lobed_p5:.3f})")
    print(f"  * Dissection Check - Sinus Min Depth Ratio    : {recs.sinus_defect_min_depth_ratio:.3f} (Significant Sinuses P25: {recs.raw_sinus_p25:.4f})")
    print("=" * 80 + "\n")

    export_diagnostic_plots(samples, recs, output_plot)

    if update_config and config_path:
        logger.info(f"Applying recommended values to: {config_path}")
        success = update_pipeline_config(config_path, recs)
        if not success:
            logger.error("Failed to update pipeline configuration.")
            sys.exit(1)
    elif not update_config:
        logger.info("Config update skipped (pass --update-config to persist recommended values).")

    return recs


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    run_calibration(
        annotations_path=args.annotations,
        output_plot=args.output_plot,
        config_path=args.config_path,
        update_config=args.update_config,
        min_area=args.min_area,
        category_id=args.category_id,
        category_name=args.category_name,
    )


if __name__ == "__main__":
    main()
