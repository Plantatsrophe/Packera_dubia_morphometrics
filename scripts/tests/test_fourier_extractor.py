"""
Unit tests for 03_fourier_extractor.py (Elliptic Fourier Analysis)
"""

import importlib
import numpy as np
import pytest

mod_fourier = importlib.import_module("scripts.morphometrics.03_fourier_extractor")
compute_efourier = mod_fourier.compute_efourier


def test_compute_efourier_circle():
    """Test that an ideal circle generates dominant first harmonic and near-zero higher harmonics."""
    theta = np.linspace(0, 2 * np.pi, 200, endpoint=False)
    r = 50.0
    x = r * np.cos(theta) + 100.0
    y = r * np.sin(theta) + 100.0
    circle_cnt = np.stack([x, y], axis=1).reshape(-1, 1, 2)

    coeffs = compute_efourier(circle_cnt, nb_harmonics=12, norm=True)
    assert coeffs is not None
    assert len(coeffs) == 48
    assert not np.isnan(coeffs).any()

    # For normalized circle: A1 ~ 1.0, D1 ~ 1.0 (or -1.0), higher harmonics near 0
    a_coeffs = coeffs[0:12]
    b_coeffs = coeffs[12:24]
    c_coeffs = coeffs[24:36]
    d_coeffs = coeffs[36:48]

    assert abs(a_coeffs[0]) > 0.8
    assert np.all(np.abs(a_coeffs[1:]) < 0.1)
    assert np.all(np.abs(b_coeffs[1:]) < 0.1)


def test_compute_efourier_invalid():
    """Test handling of degraded/degenerate contours."""
    degenerate_cnt = np.array([[[1.0, 1.0]], [[1.0, 1.0]]])
    coeffs = compute_efourier(degenerate_cnt, nb_harmonics=12, norm=True)
    assert np.isnan(coeffs).all()
