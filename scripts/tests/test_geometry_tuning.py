#!/usr/bin/env python3
"""
===============================================================================
Script: test_geometry_tuning.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Automated unit test suite verifying:
      1. Synthetic COCO dataset fixture generation with precise morphometrics:
         - 5 normal ovate leaves (W/L ≈ 0.65, Solidity ≈ 0.85)
         - 5 folded leaves (W/L ≈ 0.30, one straight edge, ChordDev ≈ 0.0)
         - 5 dissected leaves (W/L ≈ 0.50, Solidity ≈ 0.55)
      2. Empirical geometric metric extraction and percentile distributions.
      3. Safe configuration updates preserving all pre-existing sections
         (paths, taxa, models, morphometrics).
===============================================================================
"""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from typing import Any, Dict, List, Tuple

import cv2
import numpy as np
import yaml

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
import sys
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.tune_geometry_parameters import (
    COCOGeometryExtractor,
    LeafGeometrySample,
    RecommendedThresholds,
    compute_metric_summary,
    export_diagnostic_plots,
    formulate_recommended_thresholds,
    update_pipeline_config,
)


class TestGeometryTuning(unittest.TestCase):
    """Unit test suite for geometric threshold tuning and safe config persistence."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.temp_dir.name)

        self.mock_coco_path = self.tmp_path / "mock_packera_coco.json"
        self._generate_synthetic_coco_dataset(self.mock_coco_path)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    # =========================================================================
    # Synthetic Polygon & COCO Fixture Generators
    # =========================================================================

    @staticmethod
    def _create_ovate_contour(
        cx: float,
        cy: float,
        length: float = 200.0,
        width: float = 130.0,
        crenation: float = 0.42,
        n_lobes: int = 10,
    ) -> List[float]:
        """
        Generates normal ovate leaf polygon contour.
        Target morphology: W/L ≈ 0.65, Solidity ≈ 0.85.
        """
        t = np.linspace(0, 2 * np.pi, 200, endpoint=False)
        r_base = (width / 2.0) * (1.0 + 0.15 * np.cos(t))
        x = r_base * np.sin(t)
        y = -(length / 2.0) * np.cos(t)

        cren = 1.0 - crenation * (np.sin(n_lobes * t) ** 2) * (np.sin(t) ** 2)
        x = x * cren + cx
        y = y + cy

        pts = np.column_stack([x, y]).astype(np.float32)
        return pts.flatten().tolist()

    @staticmethod
    def _create_folded_contour(
        cx: float,
        cy: float,
        length: float = 200.0,
        width: float = 60.0,
    ) -> List[float]:
        """
        Generates folded hemi-blade leaf polygon contour with one straight longitudinal edge.
        Target morphology: W/L ≈ 0.30, straight lateral edge (ChordDev ≈ 0.0).
        """
        y_straight = np.linspace(-length / 2.0, length / 2.0, 50)
        x_straight = np.zeros_like(y_straight)

        y_curved = np.linspace(length / 2.0, -length / 2.0, 50)
        theta = np.pi * (y_curved + length / 2.0) / length
        x_curved = width * np.sin(theta)

        x = np.concatenate([x_straight, x_curved]) + cx
        y = np.concatenate([y_straight, y_curved]) + cy

        pts = np.column_stack([x, y]).astype(np.float32)
        return pts.flatten().tolist()

    @staticmethod
    def _create_dissected_contour(
        cx: float,
        cy: float,
        length: float = 200.0,
        width: float = 100.0,
        rachis: float = 0.14,
        sinus_exp: float = 1.45,
        n_lobes: int = 5,
    ) -> List[float]:
        """
        Generates lyrate / pinnatifid dissected leaf polygon contour with deep sinuses.
        Target morphology: W/L ≈ 0.50, Solidity ≈ 0.55.
        """
        t = np.linspace(0, 2 * np.pi, 250, endpoint=False)
        y = -(length / 2.0) * np.cos(t)

        lobe_factor = np.abs(np.sin(n_lobes * np.pi * (y + length / 2.0) / length)) ** sinus_exp
        width_at_y = (rachis + (1.0 - rachis) * (1.0 - lobe_factor)) * (width / 2.0)

        is_apex = y > (length / 2.0 - length * 0.08)
        apex_curve = (width / 2.0) * np.sqrt(
            np.maximum(0.0, 1.0 - ((y[is_apex] - (length / 2.0 - length * 0.08)) / (length * 0.08)) ** 2)
        )
        width_at_y[is_apex] = np.maximum(width_at_y[is_apex], apex_curve)

        x = np.sign(np.sin(t)) * width_at_y + cx
        y = y + cy

        pts = np.column_stack([x, y]).astype(np.float32)
        return pts.flatten().tolist()

    def _generate_synthetic_coco_dataset(self, out_path: Path) -> None:
        """
        Constructs deterministic mock COCO dataset containing:
          - 5 normal ovate leaves (W/L ≈ 0.65, Solidity ≈ 0.85)
          - 5 folded leaves (W/L ≈ 0.30, one straight edge)
          - 5 dissected leaves (W/L ≈ 0.50, Solidity ≈ 0.55)
          - 2 tiny noise fragments (< 150 px area) to test noise rejection
        """
        annotations: List[Dict[str, Any]] = []
        ann_id = 1

        # 1. 5 Normal Ovate Leaves
        for i in range(5):
            seg = self._create_ovate_contour(
                cx=150.0 + i * 30.0,
                cy=200.0 + i * 20.0,
                length=200.0 + i * 5.0,
                width=130.0 + i * 3.25,
            )
            cnt = np.array(seg).reshape(-1, 2).astype(np.int32)
            area = float(cv2.contourArea(cnt))
            annotations.append({
                "id": ann_id,
                "image_id": 1,
                "category_id": 1,
                "segmentation": [seg],
                "area": area,
                "bbox": [150.0 + i * 30.0 - 70.0, 100.0, 140.0, 210.0],
                "iscrowd": 0,
            })
            ann_id += 1

        # 2. 5 Folded Leaves
        for i in range(5):
            seg = self._create_folded_contour(
                cx=450.0 + i * 30.0,
                cy=200.0 + i * 20.0,
                length=200.0 + i * 5.0,
                width=60.0 + i * 1.5,
            )
            cnt = np.array(seg).reshape(-1, 2).astype(np.int32)
            area = float(cv2.contourArea(cnt))
            annotations.append({
                "id": ann_id,
                "image_id": 1,
                "category_id": 1,
                "segmentation": [seg],
                "area": area,
                "bbox": [450.0 + i * 30.0, 100.0, 70.0, 210.0],
                "iscrowd": 0,
            })
            ann_id += 1

        # 3. 5 Dissected Leaves
        for i in range(5):
            seg = self._create_dissected_contour(
                cx=750.0 + i * 30.0,
                cy=200.0 + i * 20.0,
                length=200.0 + i * 5.0,
                width=100.0 + i * 2.5,
            )
            cnt = np.array(seg).reshape(-1, 2).astype(np.int32)
            area = float(cv2.contourArea(cnt))
            annotations.append({
                "id": ann_id,
                "image_id": 1,
                "category_id": 1,
                "segmentation": [seg],
                "area": area,
                "bbox": [750.0 + i * 30.0 - 55.0, 100.0, 110.0, 210.0],
                "iscrowd": 0,
            })
            ann_id += 1

        # 4. 2 Noise Fragments (area < 150 px)
        noise_a = [10.0, 10.0, 18.0, 10.0, 18.0, 18.0, 10.0, 18.0]  # area ~64 px
        noise_b = [30.0, 30.0, 39.0, 30.0, 39.0, 39.0, 30.0, 39.0]  # area ~81 px
        for noise_pts in [noise_a, noise_b]:
            annotations.append({
                "id": ann_id,
                "image_id": 1,
                "category_id": 1,
                "segmentation": [noise_pts],
                "area": 70.0,
                "bbox": [10.0, 10.0, 10.0, 10.0],
                "iscrowd": 0,
            })
            ann_id += 1

        coco_data = {
            "categories": [
                {"id": 1, "name": "ideal_leaf", "supercategory": "plant_component"},
                {"id": 2, "name": "petiole", "supercategory": "plant_component"},
            ],
            "images": [
                {"id": 1, "file_name": "mock_packera_specimen.jpg", "width": 2000, "height": 2000}
            ],
            "annotations": annotations,
        }

        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(coco_data, f, indent=2)

    # =========================================================================
    # Test Cases
    # =========================================================================

    def test_mock_coco_fixture_morphometrics(self) -> None:
        """
        Verify COCOGeometryExtractor accurately extracts 15 leaf instances and filters noise,
        confirming expected morphology across ovate, folded, and dissected classes.
        """
        extractor = COCOGeometryExtractor(
            coco_path=self.mock_coco_path,
            min_contour_area=150.0,
        )
        samples = extractor.extract_samples()

        # 15 real leaf samples, 2 noise fragments rejected (< 150 px)
        self.assertEqual(len(samples), 15)

        ovate_samples = samples[0:5]
        folded_samples = samples[5:10]
        dissected_samples = samples[10:15]

        # 1. Normal ovate leaves: W/L ≈ 0.65, Solidity ≈ 0.85
        for s in ovate_samples:
            self.assertAlmostEqual(s.aspect_ratio, 0.65, delta=0.04)
            self.assertAlmostEqual(s.solidity, 0.85, delta=0.04)

        # 2. Folded leaves: W/L ≈ 0.30, straight lateral edge (ChordDev ≈ 0.0, is_flat == True)
        for s in folded_samples:
            self.assertAlmostEqual(s.aspect_ratio, 0.30, delta=0.03)
            self.assertLessEqual(s.chord_deviation_ratio, 0.02)
            self.assertTrue(s.is_flat_margin_candidate)

        # 3. Dissected leaves: W/L ≈ 0.50, Solidity ≈ 0.55
        for s in dissected_samples:
            self.assertAlmostEqual(s.aspect_ratio, 0.50, delta=0.04)
            self.assertAlmostEqual(s.solidity, 0.55, delta=0.05)
            self.assertTrue(s.is_lobed_candidate)

    def test_metric_calculation_and_percentile_distributions(self) -> None:
        """
        Verify tune_geometry_parameters outputs accurate percentiles for aspect ratio
        and solidity on the calibrated 15-leaf mock dataset.
        """
        extractor = COCOGeometryExtractor(coco_path=self.mock_coco_path)
        samples = extractor.extract_samples()

        wl_summary = compute_metric_summary("Aspect Ratio (W/L)", [s.aspect_ratio for s in samples])
        sol_summary = compute_metric_summary("Solidity", [s.solidity for s in samples])

        self.assertEqual(wl_summary.n_samples, 15)
        self.assertEqual(sol_summary.n_samples, 15)

        # Aspect ratio percentiles: 5 folded (~0.30), 5 dissected (~0.50), 5 ovate (~0.65)
        self.assertAlmostEqual(wl_summary.percentiles[10], 0.30, delta=0.04)
        self.assertAlmostEqual(wl_summary.percentiles[50], 0.50, delta=0.04)
        self.assertAlmostEqual(wl_summary.percentiles[90], 0.65, delta=0.04)

        # Solidity percentiles: 5 dissected (~0.54), 5 ovate (~0.85), 5 folded (~0.99)
        self.assertAlmostEqual(sol_summary.percentiles[10], 0.54, delta=0.05)
        self.assertAlmostEqual(sol_summary.percentiles[50], 0.85, delta=0.05)
        self.assertAlmostEqual(sol_summary.percentiles[90], 0.99, delta=0.03)

        # Guardrail checks in formulate_recommended_thresholds
        recs = formulate_recommended_thresholds(samples)
        self.assertTrue(0.38 <= recs.max_folded_aspect_ratio <= 0.46)
        self.assertTrue(0.030 <= recs.max_chord_deviation_ratio <= 0.050)
        self.assertTrue(0.48 <= recs.min_solidity_dissected <= 0.55)
        self.assertTrue(0.07 <= recs.sinus_defect_min_depth_ratio <= 0.12)

    def test_safe_config_updating_preserves_intact_keys(self) -> None:
        """
        Create a temporary copy of config/config.yaml using tempfile, run the parameter updater,
        and assert thresholds update while paths, taxa, models, and morphometrics remain 100% unaltered.
        """
        source_config = PROJECT_ROOT / "config" / "config.yaml"
        self.assertTrue(source_config.exists(), "Source config/config.yaml must exist.")

        # Read original YAML
        with open(source_config, "r", encoding="utf-8") as f:
            orig_raw_text = f.read()
            orig_dict = yaml.safe_load(orig_raw_text)

        temp_config_path = self.tmp_path / "test_config_copy.yaml"
        shutil.copy(source_config, temp_config_path)

        # Calibrate recommendations from mock dataset
        extractor = COCOGeometryExtractor(coco_path=self.mock_coco_path)
        samples = extractor.extract_samples()
        recs = formulate_recommended_thresholds(samples)

        # Execute safe config update
        update_success = update_pipeline_config(temp_config_path, recs)
        self.assertTrue(update_success, "update_pipeline_config should return True on valid YAML.")

        # Read back updated YAML
        with open(temp_config_path, "r", encoding="utf-8") as f:
            updated_raw_text = f.read()
            updated_dict = yaml.safe_load(updated_raw_text)

        # 1. Assert thresholds block updated correctly
        updated_thresh = updated_dict.get("thresholds", {})
        self.assertIn("fold_detection", updated_thresh)
        self.assertIn("dissection", updated_thresh)

        fold_cfg = updated_thresh["fold_detection"]
        diss_cfg = updated_thresh["dissection"]

        self.assertAlmostEqual(
            fold_cfg["max_folded_aspect_ratio"],
            recs.max_folded_aspect_ratio,
            places=2,
        )
        self.assertAlmostEqual(
            fold_cfg["max_chord_deviation_ratio"],
            recs.max_chord_deviation_ratio,
            places=3,
        )
        self.assertAlmostEqual(
            diss_cfg["min_solidity_dissected"],
            recs.min_solidity_dissected,
            places=2,
        )
        self.assertAlmostEqual(
            diss_cfg["sinus_defect_min_depth_ratio"],
            recs.sinus_defect_min_depth_ratio,
            places=3,
        )

        # 2. Assert all other major keys remain 100% identical and unaltered
        critical_keys = ["paths", "taxa", "models", "morphometrics", "harvesting", "segmentation"]
        for key in critical_keys:
            self.assertIn(key, updated_dict, f"Section '{key}' must remain in updated configuration.")
            self.assertEqual(
                orig_dict[key],
                updated_dict[key],
                f"Section '{key}' was corrupted or altered during geometry threshold update!",
            )

        # 3. Assert header and inline comments are preserved
        self.assertIn("# Core Filesystem Paths", updated_raw_text)
        self.assertIn("# Empirical fold-detection parameters", updated_raw_text)
        self.assertIn("# Empirical dissection & lyrate sinus parameters", updated_raw_text)

    def test_headless_diagnostic_plot_generation(self) -> None:
        """
        Verify export_diagnostic_plots runs in headless mode without opening graphical
        windows or creating permanent files in production directories.
        """
        extractor = COCOGeometryExtractor(coco_path=self.mock_coco_path)
        samples = extractor.extract_samples()
        recs = formulate_recommended_thresholds(samples)

        temp_plot_path = self.tmp_path / "headless_diagnostic.pdf"
        export_diagnostic_plots(samples, recs, temp_plot_path)

        self.assertTrue(temp_plot_path.exists())
        self.assertGreater(temp_plot_path.stat().st_size, 2000)


if __name__ == "__main__":
    unittest.main()
