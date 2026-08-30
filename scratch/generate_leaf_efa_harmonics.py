#!/usr/bin/env python3
"""
scratch/generate_leaf_efa_harmonics.py
======================================
Executes 12-harmonic Elliptic Fourier Analysis across all 4 tiers of leaf masks,
computes PCA morphospace (PC1-PC5), merges Darwin Core metadata, and exports
data/tables/leaf_efa_harmonics.csv.
"""

import math
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import cv2
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def compute_efourier(contour: np.ndarray, nb_harmonics: int = 12, norm: bool = True) -> np.ndarray:
    """Computes Kuhl & Giardina (1982) normalized Elliptic Fourier coefficients."""
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
        return np.column_stack([A, B, C, D]).flatten()

    A1, B1, C1, D1 = A[0], B[0], C[0], D[0]
    theta1 = 0.5 * math.atan2(2.0 * (A1 * B1 + C1 * D1), (A1**2 + C1**2 - B1**2 - D1**2))

    A1_star = A1 * math.cos(theta1) + B1 * math.sin(theta1)
    C1_star = C1 * math.cos(theta1) + D1 * math.sin(theta1)

    psi1 = math.atan2(C1_star, A1_star)
    E1 = math.hypot(A1_star, C1_star)
    if E1 < 1e-7:
        return np.full(nb_harmonics * 4, np.nan)

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

    return np.column_stack([norm_A, norm_B, norm_C, norm_D]).flatten()

def compute_chebyshev_poly(x_norm: np.ndarray, y_norm: np.ndarray, degree: int = 5) -> np.ndarray:
    """Fits Chebyshev orthogonal polynomials T_0 ... T_5."""
    if len(x_norm) < degree + 1:
        return np.full(degree + 1, np.nan)
    
    # Scale to [-1, 1]
    u = 2.0 * (x_norm - np.min(x_norm)) / max(np.max(x_norm) - np.min(x_norm), 1e-6) - 1.0
    
    # Chebyshev basis matrix
    T = np.zeros((len(u), degree + 1))
    T[:, 0] = 1.0
    if degree >= 1:
        T[:, 1] = u
    for d in range(2, degree + 1):
        T[:, d] = 2.0 * u * T[:, d - 1] - T[:, d - 2]
        
    try:
        coeffs, _, _, _ = np.linalg.lstsq(T, y_norm, rcond=None)
        return coeffs
    except Exception:
        return np.full(degree + 1, np.nan)

def parse_specimen_stem(stem: str) -> Tuple[str, int, int]:
    """Parses catalogNumber, plant_id, and leaf_id from mask filename stem."""
    cat_match = stem.split("_p")[0].split("_leaf")[0].split("_reflected")[0].split("_curve")[0]
    p_id = 0
    l_id = 1

    p_match = re.search(r"_p(\d+)", stem)
    if p_match:
        p_id = int(p_match.group(1))

    l_match = re.search(r"leaf(\d+)|leaf_(\d+)|_(\d+)$", stem)
    if l_match:
        digits = [g for g in l_match.groups() if g is not None]
        if digits:
            l_id = int(digits[0])

    return cat_match, p_id, l_id

def main():
    root = Path("/home/brandon/Packera_dubia_morphometrics")
    masks_dir = root / "data" / "masks"
    vouchers_path = root / "data" / "tables" / "curated_vouchers.csv"
    qc_path = root / "data" / "tables" / "leaf_extraction_qc.csv"
    output_path = root / "data" / "tables" / "leaf_efa_harmonics.csv"

    print("==================================================================")
    print("Generating Master Elliptic Fourier Harmonics Dataset (Momocs EFA)")
    print("==================================================================")

    # 1. Load curated vouchers metadata
    vouchers_df = pd.read_csv(vouchers_path) if vouchers_path.exists() else pd.DataFrame()
    print(f"Loaded {len(vouchers_df)} curated voucher records.")

    # 2. Load QC table if available for scalar traits
    qc_dict = {}
    if qc_path.exists():
        qc_df = pd.read_csv(qc_path)
        for _, row in qc_df.iterrows():
            key = f"{row.get('catalogNumber', '')}_{row.get('leaf_id', '')}"
            qc_dict[key] = row.to_dict()

    harmonic_cols = [f"{ch}{n}" for n in range(1, 13) for ch in ["A", "B", "C", "D"]]
    chebyshev_cols = [f"Chebyshev_T{d}" for d in range(6)]

    records: List[Dict[str, Any]] = []

    # Process Tier 1: Intact closed masks
    t1_files = sorted(list((masks_dir / "tier1_intact").glob("*.png")))
    print(f"Processing {len(t1_files)} Tier 1 (Direct Closed) masks...")
    for mf in t1_files:
        cat_num, p_id, l_id = parse_specimen_stem(mf.stem)
        img = cv2.imread(str(mf), cv2.IMREAD_GRAYSCALE)
        if img is None:
            continue
        contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            continue
        c = max(contours, key=cv2.contourArea)
        harmonics = compute_efourier(c, nb_harmonics=12, norm=True)

        h, w = img.shape[:2]
        area_px = float(cv2.contourArea(c))
        hull = cv2.convexHull(c)
        hull_area = float(cv2.contourArea(hull))
        solidity = area_px / max(hull_area, 1.0)
        aspect_ratio = float(w) / max(float(h), 1.0)

        rec = {
            "catalogNumber": cat_num,
            "plant_individual_id": p_id,
            "leaf_id": l_id,
            "assigned_tier": "Tier_1_Direct",
            "aspect_ratio": round(aspect_ratio, 4),
            "solidity": round(solidity, 4),
            "area_px": round(area_px, 1),
            "mask_source": str(mf.relative_to(root))
        }
        for k, h_val in enumerate(harmonics):
            rec[harmonic_cols[k]] = round(float(h_val), 6)
        for c_col in chebyshev_cols:
            rec[c_col] = np.nan
        records.append(rec)

    # Process Tier 2: Reflected masks
    t2_files = sorted(list((masks_dir / "tier2_reflected").glob("*.png")))
    print(f"Processing {len(t2_files)} Tier 2 (Reflected) masks...")
    for mf in t2_files:
        cat_num, p_id, l_id = parse_specimen_stem(mf.stem)
        img = cv2.imread(str(mf), cv2.IMREAD_GRAYSCALE)
        if img is None:
            continue
        contours, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        if not contours:
            continue
        c = max(contours, key=cv2.contourArea)
        harmonics = compute_efourier(c, nb_harmonics=12, norm=True)

        h, w = img.shape[:2]
        area_px = float(cv2.contourArea(c))
        hull = cv2.convexHull(c)
        hull_area = float(cv2.contourArea(hull))
        solidity = area_px / max(hull_area, 1.0)
        aspect_ratio = float(w) / max(float(h), 1.0)

        rec = {
            "catalogNumber": cat_num,
            "plant_individual_id": p_id,
            "leaf_id": l_id,
            "assigned_tier": "Tier_2_Reflected",
            "aspect_ratio": round(aspect_ratio, 4),
            "solidity": round(solidity, 4),
            "area_px": round(area_px, 1),
            "mask_source": str(mf.relative_to(root))
        }
        for k, h_val in enumerate(harmonics):
            rec[harmonic_cols[k]] = round(float(h_val), 6)
        for c_col in chebyshev_cols:
            rec[c_col] = np.nan
        records.append(rec)

    # Process Tier 3: Open curves
    t3_files = sorted(list((masks_dir / "tier3_open_curves").glob("*.csv")))
    print(f"Processing {len(t3_files)} Tier 3 (Open Curve) CSV series...")
    for cf in t3_files:
        cat_num, p_id, l_id = parse_specimen_stem(cf.stem)
        curve_df = pd.read_csv(cf)
        if len(curve_df) >= 10 and "x_norm" in curve_df.columns and "y_norm" in curve_df.columns:
            poly_coeffs = compute_chebyshev_poly(curve_df["x_norm"].values, curve_df["y_norm"].values, degree=5)
        else:
            poly_coeffs = np.full(6, np.nan)

        rec = {
            "catalogNumber": cat_num,
            "plant_individual_id": p_id,
            "leaf_id": l_id,
            "assigned_tier": "Tier_3_OpenCurve",
            "aspect_ratio": np.nan,
            "solidity": np.nan,
            "area_px": np.nan,
            "mask_source": str(cf.relative_to(root))
        }
        for h_col in harmonic_cols:
            rec[h_col] = np.nan
        for d, coeff in enumerate(poly_coeffs):
            rec[f"Chebyshev_T{d}"] = round(float(coeff), 6)
        records.append(rec)

    # Process Tier 4: Dense Rosettes
    t4_files = sorted(list((masks_dir / "rosettes_dense").glob("*.*")))
    print(f"Processing {len(t4_files)} Tier 4 (Dense Rosette) files...")
    for rf in t4_files:
        cat_num, p_id, l_id = parse_specimen_stem(rf.stem)
        rec = {
            "catalogNumber": cat_num,
            "plant_individual_id": p_id,
            "leaf_id": 0,
            "assigned_tier": "Tier_4_Rosette",
            "aspect_ratio": np.nan,
            "solidity": np.nan,
            "area_px": np.nan,
            "mask_source": str(rf.relative_to(root))
        }
        for h_col in harmonic_cols:
            rec[h_col] = np.nan
        for c_col in chebyshev_cols:
            rec[c_col] = np.nan
        records.append(rec)

    master_df = pd.DataFrame(records)
    print(f"Harmonized {len(master_df)} total morphometric records across all 4 tiers.")

    # 3. PCA on Closed Outlines (Tier 1 & Tier 2)
    print("Computing Morphospace Principal Component Analysis (PC1-PC5)...")
    closed_mask = master_df["assigned_tier"].isin(["Tier_1_Direct", "Tier_2_Reflected"])
    closed_df = master_df[closed_mask]

    # Select variable harmonic columns (excluding invariant normalized A1, B1, C1)
    harm_data = closed_df[harmonic_cols].values
    valid_rows = ~np.isnan(harm_data).any(axis=1)

    pca = PCA(n_components=5)
    pca_scores = pca.fit_transform(harm_data[valid_rows])

    print("=== EFA Principal Component Analysis (Momocs Standard) ===")
    for p, var in enumerate(pca.explained_variance_ratio_, 1):
        print(f"  PC{p}: {var * 100:.2f}% variance explained")
    print(f"  Total PC1-PC5 Cumulative Variance: {np.sum(pca.explained_variance_ratio_) * 100:.2f}%")

    for p in range(1, 6):
        master_df[f"PC{p}"] = np.nan

    valid_indices = closed_df[valid_rows].index
    for p in range(5):
        master_df.loc[valid_indices, f"PC{p + 1}"] = np.round(pca_scores[:, p], 6)

    # 4. Integrate Darwin Core specimen metadata
    if not vouchers_df.empty:
        meta_cols = [
            "catalogNumber", "species_raw", "determiner_raw", "determiner_tier",
            "county", "stateProvince", "latitude", "longitude",
            "pheno_sin", "pheno_cos", "regional_group"
        ]
        meta_sub = vouchers_df[[c for c in meta_cols if c in vouchers_df.columns]].drop_duplicates(subset=["catalogNumber"])
        master_df = master_df.merge(meta_sub, on="catalogNumber", how="left")

    # Order columns
    lead_cols = [
        "catalogNumber", "plant_individual_id", "leaf_id", "assigned_tier",
        "species_raw", "determiner_tier", "PC1", "PC2", "PC3", "PC4", "PC5",
        "aspect_ratio", "solidity", "area_px"
    ]
    lead_cols = [c for c in lead_cols if c in master_df.columns]
    other_cols = [c for c in master_df.columns if c not in lead_cols]
    master_df = master_df[lead_cols + other_cols]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    master_df.to_csv(output_path, index=False)
    print(f"Exported master EFA table to {output_path} (Shape: {master_df.shape})")

if __name__ == "__main__":
    main()
