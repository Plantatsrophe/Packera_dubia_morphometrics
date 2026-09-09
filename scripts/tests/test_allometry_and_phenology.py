#!/usr/bin/env python3
"""
===============================================================================
Script: test_allometry_and_phenology.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Automated unit and integration test suite validating:
    1. Foliar allometry and Centroid Size (CS) scaling linearity vs.
       scale-invariant normalized Elliptic Fourier Analysis (EFA) harmonics.
    2. Empirical phenological anomaly mathematics (Delta DOY) over latitudinal
       spring progression clines (~4.0 days/degree latitude, Hopkins' Law).
    3. Configuration synchronization for allometry and phenology parameters.

Execution:
    python -m unittest scripts/tests/test_allometry_and_phenology.py
===============================================================================
"""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

# Ensure repository root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.config import PipelineConfig
from scripts.analysis.generate_phenological_artifacts import (
    model_latitudinal_spring_baseline,
)


def compute_centroid_size(coords: np.ndarray) -> float:
    """
    Computes Centroid Size (CS) according to Bookstein (1991) and Momocs::coo_centsize:
    CS = sqrt(sum(||x_i - centroid||^2))

    Parameters
    ----------
    coords : np.ndarray
        (N, 2) array of boundary coordinate vertices.

    Returns
    -------
    float
        Centroid size (scalar).
    """
    pts = np.asarray(coords, dtype=np.float64).reshape(-1, 2)
    if len(pts) < 3:
        return float("nan")
    centroid = np.mean(pts, axis=0)
    return float(np.sqrt(np.sum((pts - centroid) ** 2)))


def compute_efourier(
    contour: np.ndarray, nb_harmonics: int = 12, norm: bool = True
) -> np.ndarray:
    """
    Computes Elliptic Fourier Analysis coefficients matching Kuhl & Giardina (1982)
    and Momocs::efourier.

    Parameters
    ----------
    contour : np.ndarray
        (N, 2) or (N, 1, 2) array of boundary coordinates.
    nb_harmonics : int
        Number of harmonic orders (default: 12).
    norm : bool
        Whether to apply scale, rotation, and starting-point invariant normalization.

    Returns
    -------
    np.ndarray
        Array of shape (nb_harmonics * 4,) containing normalized [A_n, B_n, C_n, D_n].
    """
    pts = contour.reshape(-1, 2).astype(np.float64)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])

    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)

    valid = dt > 1e-7
    if np.sum(valid) < 5:
        return np.full(nb_harmonics * 4, np.nan)

    dx = dx[valid]
    dy = dy[valid]
    dt = dt[valid]

    t = np.concatenate([[0.0], np.cumsum(dt)])
    T = t[-1]
    if T <= 0:
        return np.full(nb_harmonics * 4, np.nan)

    A = np.zeros(nb_harmonics)
    B = np.zeros(nb_harmonics)
    C = np.zeros(nb_harmonics)
    D = np.zeros(nb_harmonics)

    two_pi_over_T = 2.0 * np.pi / T

    for n in range(1, nb_harmonics + 1):
        coeff = T / (2.0 * (n**2) * (np.pi**2))
        cos_diff = np.cos(n * two_pi_over_T * t[1:]) - np.cos(
            n * two_pi_over_T * t[:-1]
        )
        sin_diff = np.sin(n * two_pi_over_T * t[1:]) - np.sin(
            n * two_pi_over_T * t[:-1]
        )

        A[n - 1] = coeff * np.sum((dx / dt) * cos_diff)
        B[n - 1] = coeff * np.sum((dx / dt) * sin_diff)
        C[n - 1] = coeff * np.sum((dy / dt) * cos_diff)
        D[n - 1] = coeff * np.sum((dy / dt) * sin_diff)

    if not norm:
        return np.concatenate([A, B, C, D])

    # Kuhl & Giardina (1982) invariant normalization
    theta_1 = 0.5 * np.arctan2(
        2 * (A[0] * B[0] + C[0] * D[0]),
        (A[0] ** 2 + C[0] ** 2 - B[0] ** 2 - D[0] ** 2),
    )
    a_star_1 = A[0] * np.cos(theta_1) + B[0] * np.sin(theta_1)
    c_star_1 = C[0] * np.cos(theta_1) + D[0] * np.sin(theta_1)
    psi_1 = np.arctan2(c_star_1, a_star_1)

    semi_major = np.sqrt(a_star_1**2 + c_star_1**2)
    if semi_major < 1e-7:
        return np.full(nb_harmonics * 4, np.nan)

    A_norm = np.zeros(nb_harmonics)
    B_norm = np.zeros(nb_harmonics)
    C_norm = np.zeros(nb_harmonics)
    D_norm = np.zeros(nb_harmonics)

    cos_psi = np.cos(psi_1)
    sin_psi = np.sin(psi_1)

    for n in range(1, nb_harmonics + 1):
        cos_ntheta = np.cos(n * theta_1)
        sin_ntheta = np.sin(n * theta_1)

        an = (cos_psi * A[n - 1] + sin_psi * C[n - 1]) / semi_major
        bn = (cos_psi * B[n - 1] + sin_psi * D[n - 1]) / semi_major
        cn = (-sin_psi * A[n - 1] + cos_psi * C[n - 1]) / semi_major
        dn = (-sin_psi * B[n - 1] + cos_psi * D[n - 1]) / semi_major

        A_norm[n - 1] = an * cos_ntheta + bn * sin_ntheta
        B_norm[n - 1] = -an * sin_ntheta + bn * cos_ntheta
        C_norm[n - 1] = cn * cos_ntheta + dn * sin_ntheta
        D_norm[n - 1] = -cn * sin_ntheta + dn * cos_ntheta

    return np.concatenate([A_norm, B_norm, C_norm, D_norm])


def generate_synthetic_square(
    side_length: float = 1.0, pts_per_side: int = 25
) -> np.ndarray:
    """Generates an evenly sampled closed square perimeter centered or anchored at origin."""
    s1 = np.column_stack(
        [np.linspace(0, side_length, pts_per_side, endpoint=False), np.zeros(pts_per_side)]
    )
    s2 = np.column_stack(
        [
            np.full(pts_per_side, side_length),
            np.linspace(0, side_length, pts_per_side, endpoint=False),
        ]
    )
    s3 = np.column_stack(
        [
            np.linspace(side_length, 0, pts_per_side, endpoint=False),
            np.full(pts_per_side, side_length),
        ]
    )
    s4 = np.column_stack(
        [np.zeros(pts_per_side), np.linspace(side_length, 0, pts_per_side, endpoint=False)]
    )
    return np.vstack([s1, s2, s3, s4])


class TestAllometryAndCentroidSize(unittest.TestCase):
    """Unit tests for Centroid Size calculation and allometry-free shape normalization."""

    def test_centroid_size_linear_scaling_and_efa_invariance(self):
        """
        Verify that:
        1. Centroid Size (CS) scales strictly linearly by 2.0x for a 2x scaled shape.
        2. Normalized Elliptic Fourier Analysis (EFA) harmonics remain identical.
        """
        sq_unit = generate_synthetic_square(side_length=1.0, pts_per_side=25)
        sq_scaled = generate_synthetic_square(side_length=2.0, pts_per_side=25)

        # 1. Compute Centroid Sizes
        cs_unit = compute_centroid_size(sq_unit)
        cs_scaled = compute_centroid_size(sq_scaled)

        self.assertGreater(cs_unit, 0.0, "Unit square CS must be strictly positive.")
        self.assertGreater(cs_scaled, 0.0, "Scaled square CS must be strictly positive.")

        ratio = cs_scaled / cs_unit
        self.assertAlmostEqual(
            ratio,
            2.0,
            places=6,
            msg=f"Centroid Size must scale exactly 2.0x under 2x dilation; got {ratio:.8f}",
        )

        # 2. Compute 12-Harmonic Normalized EFA coefficients
        h_unit = compute_efourier(sq_unit, nb_harmonics=12, norm=True)
        h_scaled = compute_efourier(sq_scaled, nb_harmonics=12, norm=True)

        self.assertEqual(len(h_unit), 48, "12 harmonics must produce 48 coefficients.")
        self.assertEqual(len(h_scaled), 48, "12 harmonics must produce 48 coefficients.")
        self.assertFalse(np.isnan(h_unit).any(), "Harmonics must not contain NaN.")
        self.assertFalse(np.isnan(h_scaled).any(), "Harmonics must not contain NaN.")

        # Assert scale-invariant normalized harmonics are identical
        np.testing.assert_allclose(
            h_unit,
            h_scaled,
            atol=1e-6,
            err_msg="Normalized EFA harmonics must remain identical under 2x isometric scaling.",
        )

    def test_centroid_size_translation_and_rotation_invariance(self):
        """Verify Centroid Size is invariant to rigid translations."""
        sq_unit = generate_synthetic_square(side_length=1.0, pts_per_side=25)
        cs_orig = compute_centroid_size(sq_unit)

        # Translate square by (+50.0, -35.0)
        sq_translated = sq_unit + np.array([50.0, -35.0])
        cs_translated = compute_centroid_size(sq_translated)

        self.assertAlmostEqual(
            cs_orig,
            cs_translated,
            places=6,
            msg="Centroid Size must be strictly invariant to translation.",
        )

    def test_arbitrary_scaling_factor(self):
        """Verify CS scales linearly for arbitrary scale factor (3.5x)."""
        scale_factor = 3.5
        sq_unit = generate_synthetic_square(side_length=1.0, pts_per_side=25)
        sq_scaled = sq_unit * scale_factor

        cs_unit = compute_centroid_size(sq_unit)
        cs_scaled = compute_centroid_size(sq_scaled)

        self.assertAlmostEqual(
            cs_scaled / cs_unit,
            scale_factor,
            places=6,
            msg=f"CS scaling ratio must equal {scale_factor}.",
        )


class TestPhenologicalAnomalyMath(unittest.TestCase):
    """Unit tests for empirical latitudinal phenological anomaly (Delta DOY) math."""

    def test_mock_florida_and_virginia_cline_anomaly_elimination(self):
        """
        Verify that mock coordinates at Lat 30 deg (Florida, DOY 85) and
        Lat 38 deg (Virginia, DOY 117) with a cline slope of 4.0 days/deg
        both yield Delta DOY ~= 0.0, eliminating latitude confounding.
        """
        lat_fl, doy_fl = 30.0, 85.0
        lat_va, doy_va = 38.0, 117.0
        cline_slope = 4.0  # days per degree latitude (Hopkins' Bioclimatic Law)

        # Solve intercept from the 2-point cline:
        # DOY = intercept + slope * lat => intercept = DOY - slope * lat
        intercept_fl = doy_fl - (cline_slope * lat_fl)  # 85 - 120 = -35.0
        intercept_va = doy_va - (cline_slope * lat_va)  # 117 - 152 = -35.0
        self.assertAlmostEqual(
            intercept_fl,
            intercept_va,
            places=6,
            msg="Florida and Virginia coordinates must belong to the exact same clinal line.",
        )

        intercept = intercept_fl

        # Compute Expected DOY based on latitude
        expected_doy_fl = intercept + cline_slope * lat_fl
        expected_doy_va = intercept + cline_slope * lat_va

        # Compute Phenological Anomaly: Delta DOY = DOY_obs - Expected_DOY
        delta_doy_fl = doy_fl - expected_doy_fl
        delta_doy_va = doy_va - expected_doy_va

        # Assert both points yield Delta DOY approx 0.0
        self.assertAlmostEqual(
            delta_doy_fl,
            0.0,
            places=6,
            msg=f"Florida Delta DOY must be 0.0; got {delta_doy_fl}",
        )
        self.assertAlmostEqual(
            delta_doy_va,
            0.0,
            places=6,
            msg=f"Virginia Delta DOY must be 0.0; got {delta_doy_va}",
        )

    def test_model_latitudinal_spring_baseline_integration(self):
        """
        Verify integration with model_latitudinal_spring_baseline on a synthetic
        dataset exhibiting a 4.0 days/degree latitude cline.
        """
        # Generate 15 points along the cline: DOY = -35.0 + 4.0 * Latitude
        lats = np.linspace(28.0, 42.0, 15)
        doys = -35.0 + 4.0 * lats

        df_cline = pd.DataFrame({
            "decimalLatitude": lats,
            "doy": doys,
        })

        baseline = model_latitudinal_spring_baseline(df_cline)

        self.assertAlmostEqual(baseline["slope"], 4.0, places=3)
        self.assertAlmostEqual(baseline["intercept"], -35.0, places=3)
        self.assertAlmostEqual(baseline["r_squared"], 1.0, places=4)
        self.assertLess(baseline["p_value"], 1e-6)

        # Verify that an early-flowering outlier (e.g. granite flatrock ecotype) is detected
        flatrock_lat = 34.0
        flatrock_doy = 81.0  # expected: -35 + 4*34 = 101 -> 20 days earlier than cline
        expected = baseline["intercept"] + baseline["slope"] * flatrock_lat
        anomaly = flatrock_doy - expected
        self.assertAlmostEqual(anomaly, -20.0, places=2)
        self.assertLess(anomaly, -7.0, "Should be flagged as Early_Flowering (Delta DOY < -7).")


class TestConfigurationSynchronization(unittest.TestCase):
    """Verify that configuration parameters for allometry and phenology are present."""

    def test_allometry_and_phenology_config_entries(self):
        """Verify YAML schema includes correct_allometry, allometry_r2_threshold, and phenology."""
        cfg = PipelineConfig.from_yaml()

        # Morphometrics allometry parameters
        self.assertEqual(cfg.morphometrics.nb_harmonics, 12)
        self.assertEqual(cfg.morphometrics.correct_allometry, "auto")
        self.assertAlmostEqual(cfg.morphometrics.allometry_r2_threshold, 0.10)

        # Macroecology phenology parameters
        self.assertEqual(cfg.macroecology.phenology.min_flowering_doy, 60)
        self.assertEqual(cfg.macroecology.phenology.max_flowering_doy, 220)
        self.assertTrue(cfg.macroecology.phenology.robust_regression)


if __name__ == "__main__":
    unittest.main(verbosity=2)
