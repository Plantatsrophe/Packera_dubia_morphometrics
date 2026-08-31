#!/usr/bin/env python3
"""
===============================================================================
Script: 03_fourier_extractor.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Phase 3: Python implementation of normalized Elliptic Fourier Analysis (EFA).
    Extracts 12-harmonic Fourier coefficients from 4-tier leaf boundary contours,
    fits morphometric PCA, merges Darwin Core voucher metadata, and exports
    standardized morphospace tables.
===============================================================================
"""

from __future__ import annotations

import argparse
import logging
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import cv2
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("FourierExtractor")


def compute_efourier(
    contour: np.ndarray,
    nb_harmonics: int = 12,
    norm: bool = True
) -> np.ndarray:
    """
    Computes Kuhl & Giardina (1982) normalized Elliptic Fourier coefficients.

    Args:
        contour: (N, 1, 2) or (N, 2) array of boundary coordinate vertices.
        nb_harmonics: Number of harmonic Fourier orders (default: 12).
        norm: Whether to normalize for size, orientation, and starting point.

    Returns:
        np.ndarray: Array of shape (nb_harmonics * 4,) containing [A_n, B_n, C_n, D_n].
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
        coeff = T / (2.0 * (n ** 2) * (np.pi ** 2))
        cos_diff = np.cos(n * two_pi_over_T * t[1:]) - np.cos(n * two_pi_over_T * t[:-1])
        sin_diff = np.sin(n * two_pi_over_T * t[1:]) - np.sin(n * two_pi_over_T * t[:-1])

        A[n - 1] = coeff * np.sum((dx / dt) * cos_diff)
        B[n - 1] = coeff * np.sum((dx / dt) * sin_diff)
        C[n - 1] = coeff * np.sum((dy / dt) * cos_diff)
        D[n - 1] = coeff * np.sum((dy / dt) * sin_diff)

    if not norm:
        return np.concatenate([A, B, C, D])

    # Kuhl & Giardina (1982) invariant normalization
    theta_1 = 0.5 * np.arctan2(2 * (A[0] * B[0] + C[0] * D[0]), (A[0] ** 2 + C[0] ** 2 - B[0] ** 2 - D[0] ** 2))
    a_star_1 = A[0] * np.cos(theta_1) + B[0] * np.sin(theta_1)
    c_star_1 = C[0] * np.cos(theta_1) + D[0] * np.sin(theta_1)
    psi_1 = np.arctan2(c_star_1, a_star_1)

    semi_major = np.sqrt(a_star_1 ** 2 + c_star_1 ** 2)
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


def process_leaf_masks_directory(
    masks_dir: Path,
    vouchers_csv: Optional[Path] = None,
    nb_harmonics: int = 12
) -> pd.DataFrame:
    """
    Computes EFA harmonic profiles across all leaf masks in a directory.
    """
    masks_dir = Path(masks_dir)
    mask_files = sorted(list(masks_dir.glob("*.png")) + list(masks_dir.glob("*.jpg")))
    logger.info(f"Processing {len(mask_files)} leaf masks from {masks_dir}...")

    vouchers_df = pd.read_csv(vouchers_csv) if vouchers_csv and Path(vouchers_csv).exists() else None
    voucher_map = {}
    if vouchers_df is not None and "catalogNumber" in vouchers_df.columns:
        for _, row in vouchers_df.iterrows():
            voucher_map[str(row["catalogNumber"]).strip()] = row.to_dict()

    harm_cols = []
    for prefix in ["A", "B", "C", "D"]:
        for h in range(1, nb_harmonics + 1):
            harm_cols.append(f"{prefix}{h}")

    records = []
    for mf in mask_files:
        img = cv2.imread(str(mf), cv2.IMREAD_GRAYSCALE)
        if img is None:
            continue

        contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            continue

        largest_cnt = max(contours, key=cv2.contourArea)
        if len(largest_cnt) < 10:
            continue

        coeffs = compute_efourier(largest_cnt, nb_harmonics=nb_harmonics, norm=True)
        if np.isnan(coeffs).any():
            continue

        cat_num = mf.stem.split("_")[0]
        meta = voucher_map.get(cat_num, {})

        rec = {
            "mask_file": mf.name,
            "catalogNumber": cat_num,
            "species_raw": meta.get("species_raw", meta.get("species", "Unknown")),
            "determiner_tier": meta.get("determiner_tier", "Tier_3_Bronze"),
            "regional_group": meta.get("regional_group", "Other_US"),
        }
        for col_name, val in zip(harm_cols, coeffs):
            rec[col_name] = round(float(val), 6)

        records.append(rec)

    df = pd.DataFrame(records)
    logger.info(f"Computed EFA harmonics for {len(df)} leaves.")
    return df


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Phase 3: Normalized Elliptic Fourier Analysis")
    parser.add_argument("--masks-dir", type=Path, default=Path("data/raw_annotations/masks"))
    parser.add_argument("--vouchers-csv", type=Path, default=Path("data/tables/curated_vouchers.csv"))
    parser.add_argument("--output-csv", type=Path, default=Path("data/tables/leaf_efa_harmonics.csv"))
    parser.add_argument("--nb-harmonics", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    df = process_leaf_masks_directory(args.masks_dir, args.vouchers_csv, args.nb_harmonics)
    if not df.empty:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.output_csv, index=False)
        logger.info(f"Saved EFA harmonics table -> {args.output_csv}")


if __name__ == "__main__":
    main()
