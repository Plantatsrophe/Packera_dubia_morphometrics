#!/usr/bin/env python3
"""
scratch/test_efa.py
===================
Tests Elliptic Fourier Analysis and 4-tier extraction on dataset masks.
"""

import math
from pathlib import Path
import cv2
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def compute_efourier(contour: np.ndarray, nb_harmonics: int = 12, norm: bool = True) -> np.ndarray:
    """
    Computes normalized Elliptic Fourier Analysis coefficients matching Kuhl & Giardina (1982)
    and Momocs::efourier.
    """
    pts = contour.reshape(-1, 2).astype(np.float64)
    if not np.allclose(pts[0], pts[-1]):
        pts = np.vstack([pts, pts[0]])

    dx = np.diff(pts[:, 0])
    dy = np.diff(pts[:, 1])
    dt = np.hypot(dx, dy)
    
    # Remove consecutive duplicate points
    valid = dt > 1e-7
    if np.sum(valid) < 5:
        return np.zeros(nb_harmonics * 4)
    
    dx = dx[valid]
    dy = dy[valid]
    dt = dt[valid]
    
    t = np.concatenate([[0.0], np.cumsum(dt)])
    T = t[-1]
    if T <= 0:
        return np.zeros(nb_harmonics * 4)

    A = np.zeros(nb_harmonics)
    B = np.zeros(nb_harmonics)
    C = np.zeros(nb_harmonics)
    D = np.zeros(nb_harmonics)

    two_pi_over_T = 2.0 * np.pi / T

    for n in range(1, nb_harmonics + 1):
        coeff = T / (2.0 * (n ** 2) * (np.pi ** 2))
        cos_diff = np.cos(n * two_pi_over_T * t[1:]) - np.cos(n * two_pi_over_T * t[:-1])
        sin_diff = np.sin(n * two_pi_over_T * t[1:]) - np.sin(n * two_pi_over_T * t[:-1])

        A[n - 1] = coeff * np.sum((dx / dt) * cos_diff)
        B[n - 1] = coeff * np.sum((dx / dt) * sin_diff)
        C[n - 1] = coeff * np.sum((dy / dt) * cos_diff)
        D[n - 1] = coeff * np.sum((dy / dt) * sin_diff)

    if not norm:
        harmonics = np.column_stack([A, B, C, D]).flatten()
        return harmonics

    # Normalization matching Momocs (Kuhl & Giardina 1982)
    A1, B1, C1, D1 = A[0], B[0], C[0], D[0]
    theta1 = 0.5 * math.atan2(2.0 * (A1 * B1 + C1 * D1), (A1**2 + C1**2 - B1**2 - D1**2))
    
    A1_star = A1 * math.cos(theta1) + B1 * math.sin(theta1)
    C1_star = C1 * math.cos(theta1) + D1 * math.sin(theta1)
    
    psi1 = math.atan2(C1_star, A1_star)
    E1 = math.hypot(A1_star, C1_star)
    if E1 < 1e-7:
        return np.zeros(nb_harmonics * 4)

    cos_psi1 = math.cos(psi1)
    sin_psi1 = math.sin(psi1)
    R_psi = np.array([[cos_psi1, sin_psi1], [-sin_psi1, cos_psi1]], dtype=np.float64)

    norm_A = np.zeros(nb_harmonics)
    norm_B = np.zeros(nb_harmonics)
    norm_C = np.zeros(nb_harmonics)
    norm_D = np.zeros(nb_harmonics)

    for n in range(1, nb_harmonics + 1):
        M_n = np.array([[A[n-1], B[n-1]], [C[n-1], D[n-1]]], dtype=np.float64)
        cos_ntheta = math.cos(n * theta1)
        sin_ntheta = math.sin(n * theta1)
        R_ntheta = np.array([[cos_ntheta, -sin_ntheta], [sin_ntheta, cos_ntheta]], dtype=np.float64)

        norm_M = (1.0 / E1) * (R_psi @ M_n @ R_ntheta)
        norm_A[n - 1] = norm_M[0, 0]
        norm_B[n - 1] = norm_M[0, 1]
        norm_C[n - 1] = norm_M[1, 0]
        norm_D[n - 1] = norm_M[1, 1]

    harmonics = np.column_stack([norm_A, norm_B, norm_C, norm_D]).flatten()
    return harmonics

def run_test():
    mask_files = list(Path("data/masks/tier1_intact").glob("*.png"))[:10]
    print(f"Testing on {len(mask_files)} mask files...")
    for mf in mask_files:
        img = cv2.imread(str(mf), cv2.IMREAD_GRAYSCALE)
        contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if contours:
            c = max(contours, key=cv2.contourArea)
            h = compute_efourier(c, 12, norm=True)
            print(f"Mask {mf.name}: A1={h[0]:.3f}, B1={h[1]:.3f}, C1={h[2]:.3f}, D1={h[3]:.3f}, A2={h[4]:.3f}")

if __name__ == "__main__":
    run_test()
