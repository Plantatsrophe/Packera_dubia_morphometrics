#!/usr/bin/env python3
"""
===============================================================================
Module: step02_figures.py
Project: Packera dubia Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Publication-grade Step 02 segmentation QC distributions and taxon-stratified
    silhouette montage generation (vector PDF + 300-DPI PNG).
===============================================================================
"""

from __future__ import annotations

import json
import logging
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.vision.leaf_geometry_utils import (
    compute_geometric_metrics,
    is_botanical_dissection,
)
from scripts.vision.lm3_inference_engines import _ensure_ect_module_loaded

logger = logging.getLogger("SegmentAndExtract")

__all__ = [
    "_atomic_save_matplotlib_fig",
    "_standardize_packera_species_py",
    "_parse_bool_val",
    "_ASSET_ANCHOR_DIRS",
    "_resolve_asset_path",
    "_is_valid_smooth_closed_mask",
    "_discover_downstream_table",
    "generate_step02_publication_figures",
    "regenerate_step02_figures_from_tables",
    "export_taxon_ect_exemplars",
]


def _atomic_save_matplotlib_fig(fig: Any, out_base_path: Path, dpi: int = 300) -> Tuple[Path, Path]:
    """Atomically save a Matplotlib figure as both vector .pdf and 300-DPI .png."""
    out_base_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = out_base_path.with_suffix(".pdf")
    png_path = out_base_path.with_suffix(".png")
    tmp_pdf = out_base_path.with_name(f".{out_base_path.stem}_{os.getpid()}_tmp.pdf")
    tmp_png = out_base_path.with_name(f".{out_base_path.stem}_{os.getpid()}_tmp.png")
    try:
        fig.savefig(str(tmp_pdf), format="pdf", bbox_inches="tight")
        tmp_pdf.replace(pdf_path)
        fig.savefig(str(tmp_png), format="png", dpi=dpi, bbox_inches="tight")
        tmp_png.replace(png_path)
    finally:
        for t in (tmp_pdf, tmp_png):
            if t.exists():
                t.unlink(missing_ok=True)
    return pdf_path, png_path


def _standardize_packera_species_py(raw_val: Any) -> str:
    """Delegate taxonomic standardization to scripts.core.harvester.standardize_packera_taxon with generic binomial fallback."""
    if raw_val is None or pd.isna(raw_val):
        return "Unassigned_Packera"
    s_raw = str(raw_val).strip()
    if not s_raw or s_raw.lower() in {"none", "nan", "null", "na", "unknown", "unassigned_packera"}:
        return "Unassigned_Packera"
    try:
        from scripts.core.harvester import standardize_packera_taxon

        std = standardize_packera_taxon(s_raw)
        if std and str(std).strip() and str(std).strip().lower() not in {"none", "nan", "null", "na", "unknown", "unassigned_packera"}:
            return str(std).strip()
    except Exception:
        pass
    cleaned = s_raw.split("(")[0].strip()
    m_binom = re.match(r"^([A-Za-z][a-z]+)\s+(?:([x×])\s+)?([a-z][a-z-]+)\b", cleaned)
    if m_binom:
        genus = m_binom.group(1).capitalize()
        hyb = "× " if m_binom.group(2) else ""
        return f"{genus} {hyb}{m_binom.group(3).lower()}"
    if re.match(r"^(?:Morphotype|Ecotype)_[A-Za-z0-9_]+$", cleaned):
        return cleaned
    return "Unassigned_Packera"


def _parse_bool_val(val: Any, default: bool = False) -> bool:
    """Parse boolean flags exported by Python or R ('TRUE'/'FALSE', True/False, 1/0)."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return default
    if isinstance(val, (bool, np.bool_)):
        return bool(val)
    s = str(val).strip().upper()
    if s in {"TRUE", "T", "1", "1.0", "YES"}:
        return True
    if s in {"FALSE", "F", "0", "0.0", "NO"}:
        return False
    return default


_ASSET_ANCHOR_DIRS = ("masks", "contours", "cropped_patches")


def _resolve_asset_path(p_str: Any, out_dir: Path) -> Optional[Path]:
    """
    Resolve mask or contour asset path across absolute, run-scoped, and workspace-root paths.

    SQLite `extracted_leaves` / resumed Step 02 manifests may carry absolute paths into a
    previous (since rebuilt or deleted) `runs/<run_id>/data/` directory. When the literal path
    is missing, the trailing `masks/...`, `contours/...`, or `cropped_patches/...` segment is
    re-rooted onto the current run-local `data/` store and then the canonical `data/` store.
    """
    raw = str(p_str or "").strip()
    if not raw or raw.lower() in {"nan", "none", "na"}:
        return None
    p = Path(raw)
    if p.exists():
        return p
    for base in (out_dir.parent.parent, out_dir.parent, PROJECT_ROOT):
        cand = base / raw
        if cand.exists():
            return cand
    parts = p.parts
    anchor_idx = next((i for i in range(len(parts) - 1, -1, -1) if parts[i] in _ASSET_ANCHOR_DIRS), None)
    if anchor_idx is not None and anchor_idx < len(parts) - 1:
        rel_tail = Path(*parts[anchor_idx:])
        for data_root in (out_dir.parent.parent / "data", out_dir.parent / "data", PROJECT_ROOT / "data"):
            cand = data_root / rel_tail
            if cand.exists():
                return cand
    return None


def _is_valid_smooth_closed_mask(m_path_str: Any, out_dir: Path) -> bool:
    """
    Verify that a binary silhouette mask has a single smooth closed contour
    without internal holes or clipped border edges.
    Returns True if no on-disk mask path is configured (e.g., synthetic contour-only fixtures).
    """
    resolved = _resolve_asset_path(m_path_str, out_dir)
    raw_str = str(m_path_str or "").strip()
    if resolved is None:
        # If an explicit path was provided on a real dataset but is missing on disk, only allow
        # fallback when no mask file exists in synthetic test environments.
        return not bool(raw_str)
    mask_img = cv2.imread(str(resolved), cv2.IMREAD_GRAYSCALE)
    if mask_img is None:
        return False
    if mask_img.ndim == 3:
        mask_img = mask_img[:, :, 0]
    fg = (mask_img > 0).astype(np.uint8)
    total_fg = int(np.count_nonzero(fg))
    if total_fg < 50:
        return False

    h_m, w_m = fg.shape[:2]
    # 1. Reject masks with clipped border edges along the canvas boundary
    max_edge_touch = max(3, int(0.02 * max(h_m, w_m)))
    if (
        int(np.count_nonzero(fg[0, :])) > max_edge_touch
        or int(np.count_nonzero(fg[-1, :])) > max_edge_touch
        or int(np.count_nonzero(fg[:, 0])) > max_edge_touch
        or int(np.count_nonzero(fg[:, -1])) > max_edge_touch
    ):
        return False

    # 2. Require a single significant connected component (or dominant >=92% primary component)
    num_l, lbls_cc, stats_cc, _ = cv2.connectedComponentsWithStats(fg, connectivity=8)
    if num_l <= 1:
        return False
    min_sig_area = max(25, int(0.01 * total_fg))
    cc_areas = stats_cc[1:, cv2.CC_STAT_AREA]
    sig_comps = int(np.sum(cc_areas >= min_sig_area))
    if sig_comps != 1:
        best_cc = 1 + int(np.argmax(cc_areas))
        best_cc_area = int(stats_cc[best_cc, cv2.CC_STAT_AREA])
        if float(best_cc_area) / float(max(total_fg, 1)) >= 0.85:
            fg = (lbls_cc == best_cc).astype(np.uint8)
            total_fg = best_cc_area
            min_sig_area = max(25, int(0.01 * total_fg))
        else:
            return False

    # 3. Reject internal holes using 2-level contour hierarchy (RETR_CCOMP), allowing tiny
    #    natural herbivory pinholes (< 2.5% of blade area) on deeply dissected leaves
    cnts, hier = cv2.findContours(fg, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    if not cnts or hier is None:
        return False
    outer_cnts = [
        cnts[idx_c]
        for idx_c in range(len(cnts))
        if int(hier[0][idx_c][3]) == -1 and float(cv2.contourArea(cnts[idx_c])) >= min_sig_area
    ]
    if len(outer_cnts) != 1:
        return False

    main_cnt = outer_cnts[0]
    if len(main_cnt) < 12:
        return False
    area_c = float(cv2.contourArea(main_cnt))
    perim_c = float(cv2.arcLength(main_cnt, True))
    if area_c <= 0.0 or perim_c <= 0.0:
        return False
    hull_c = cv2.convexHull(main_cnt)
    hull_area_c = float(cv2.contourArea(hull_c))
    sol_c = area_c / max(hull_area_c, 1.0)

    min_hole_area = max(12.0, (0.025 if sol_c < 0.80 else 0.005) * float(total_fg))
    for idx_c, c in enumerate(cnts):
        parent_idx = int(hier[0][idx_c][3])
        if parent_idx != -1 and float(cv2.contourArea(c)) >= min_hole_area:
            return False

    # 4. Require smooth closed contour geometry
    circ = (4.0 * math.pi * area_c) / (perim_c * perim_c)
    if circ < 0.12:
        if not (
            circ >= 0.025
            and (
                sol_c < 0.78
                or is_botanical_dissection(main_cnt, solidity=sol_c, check_solidity_range=False)
            )
        ):
            return False
    return True


def _discover_downstream_table(filename: str, out_dir: Path, active_cats: Set[str]) -> Optional[Path]:
    """
    Locate a table in `out_dir.parent.parent / 'data' / 'tables'` (or run-local table dirs),
    falling back to `PROJECT_ROOT / 'data' / 'tables'` when `out_dir` resides under `PROJECT_ROOT`
    or when the table shares catalogNumbers with `active_cats`.
    """
    local_dirs = [
        out_dir.parent.parent / "data" / "tables",
        out_dir.parent / "data" / "tables",
        out_dir.parent / "tables",
        out_dir / "tables",
    ]
    for d in local_dirs:
        cand = d / filename
        if cand.exists() and cand.stat().st_size > 0:
            return cand

    root_cand = PROJECT_ROOT / "data" / "tables" / filename
    if root_cand.exists() and root_cand.stat().st_size > 0:
        try:
            out_resolved = out_dir.resolve()
            root_resolved = PROJECT_ROOT.resolve()
            if out_resolved == root_resolved or root_resolved in out_resolved.parents:
                return root_cand
        except Exception:
            pass
        if active_cats:
            try:
                tdf = pd.read_csv(root_cand, usecols=lambda c: c == "catalogNumber", dtype=str)
                if "catalogNumber" in tdf.columns:
                    tcats = set(tdf["catalogNumber"].astype(str).str.strip())
                    if tcats & active_cats:
                        return root_cand
            except Exception:
                pass
    return None


def generate_step02_publication_figures(
    extracted_df: pd.DataFrame,
    failed_df: pd.DataFrame,
    repro_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    figures_dir: Union[str, Path],
) -> Dict[str, Path]:
    """
    Generate Step 02 publication-ready QC and silhouette montage figures (vector .pdf + 300-DPI .png)
    using colorblind-safe Okabe-Ito palettes and italicized Packera taxon labels:
      1. fig02a_segmentation_qc_distributions.{pdf,png}:
         (A) Lamina Solidity vs. Blade Aspect Ratio (Length/Width) colored and shaped by Species
             (hollow markers) with shaded gatekeeper zones (Tier 1 Intact vs. Botanical Dissection Zone).
         (B) Taxon-Stratified Reproductive Trait Distributions:
             (B1) Involucre Aspect Ratio (L/W) boxplots + hollow jittered points by species.
             (B2) Anthesis Phenological State proportions across species.
      2. fig02b_silhouette_montage.{pdf,png}:
         Taxon-stratified grid (rows = active Packera species present in df_plot;
         4 representative Tier_1_Direct leaves per row spanning the 20th, 40th, 60th, and 80th
         percentiles of blade aspect ratio/dissection within q +/- 0.08 quantile windows),
         with 10 mm physical scale bar on calibrated tiles, homologized p_base attachment point,
         and midrib axis.
    """
    import matplotlib
    matplotlib.use("Agg")
    _step02_mod = sys.modules.get("scripts.pipeline.02_segment_and_extract")
    _save_fig_fn = (
        getattr(_step02_mod, "_atomic_save_matplotlib_fig", _atomic_save_matplotlib_fig)
        if _step02_mod is not None
        else _atomic_save_matplotlib_fig
    )
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.patches import Rectangle

    out_dir = Path(figures_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved_paths: Dict[str, Path] = {}

    if extracted_df is None or len(extracted_df) == 0:
        return saved_paths

    EXTENDED_OKABE_ITO = (
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7",
        "#E69F00",
        "#56B4E9",
        "#882255",
        "#44AA99",
        "#332288",
        "#117733",
    )
    HOLLOW_MARKER_CYCLE = ("o", "^", "s", "D", "v", "P", "X", "<", ">", "h")
    SPECIES_PALETTE: Dict[str, str] = {}
    SPECIES_MARKERS: Dict[str, str] = {}
    SPECIES_ITALIC: Dict[str, str] = {}

    df_plot = extracted_df.copy()
    df_plot["cat_str"] = df_plot.get("catalogNumber", pd.Series("", index=df_plot.index)).astype(str).str.strip()
    df_plot["leaf_id_int"] = pd.to_numeric(df_plot.get("leaf_id", 1), errors="coerce").fillna(1).astype(int)
    active_cats = set(df_plot["cat_str"].unique()) - {""}

    # Resolve voucher -> standardized species & determiner_tier map from curated_vouchers.csv if available
    cat_to_species: Dict[str, str] = {}
    cat_to_det_tier: Dict[str, str] = {}
    vp = _discover_downstream_table("curated_vouchers.csv", out_dir, active_cats)
    if vp is not None and vp.exists():
        try:
            vdf = pd.read_csv(vp, dtype=str)
            if "catalogNumber" in vdf.columns:
                tx_col = (
                    "species_standardized"
                    if "species_standardized" in vdf.columns
                    else ("scientificName" if "scientificName" in vdf.columns else "species_raw")
                )
                for _, vrow in vdf.iterrows():
                    c_val = str(vrow.get("catalogNumber", "")).strip()
                    if not c_val:
                        continue
                    if tx_col in vdf.columns and pd.notna(vrow.get(tx_col)):
                        cat_to_species[c_val] = _standardize_packera_species_py(vrow.get(tx_col))
                    if "determiner_tier" in vdf.columns and pd.notna(vrow.get("determiner_tier")):
                        cat_to_det_tier[c_val] = str(vrow.get("determiner_tier")).strip()
        except Exception:
            pass

    # 1. Enforce assigned_tier == "Tier_1_Direct"
    tier_series = df_plot.get("assigned_tier", pd.Series("Tier_1_Direct", index=df_plot.index)).astype(str).str.strip()
    df_plot["norm_tier"] = tier_series.apply(
        lambda t: "Tier_1_Direct"
        if t.lower() in {"tier_1_direct", "tier1", "tier 1", "tier 1 (dissected)", "tier 1 (bipinnatisect)"}
        else ("Tier_2_Reflected" if "2" in t else t)
    )
    df_plot = df_plot[df_plot["norm_tier"] == "Tier_1_Direct"].copy()

    # 2. Inner-join / filter against downstream Step 03 & Step 04 QC tables when present
    ind_efa_path = _discover_downstream_table("leaf_efa_harmonics_individual.csv", out_dir, active_cats)
    if ind_efa_path is not None and ind_efa_path.exists():
        try:
            ind_df = pd.read_csv(ind_efa_path)
            if "catalogNumber" in ind_df.columns and "leaf_id" in ind_df.columns:
                ind_df["cat_str"] = ind_df["catalogNumber"].astype(str).str.strip()
                ind_df["leaf_id_int"] = pd.to_numeric(ind_df["leaf_id"], errors="coerce").fillna(-1).astype(int)
                if "assigned_tier" in ind_df.columns:
                    ind_tier = ind_df["assigned_tier"].astype(str).str.strip().str.lower()
                    ind_df = ind_df[ind_tier.isin({"tier_1_direct", "tier1", "tier 1"})]
                if "is_contour_outlier" in ind_df.columns:
                    ind_df = ind_df[~ind_df["is_contour_outlier"].apply(_parse_bool_val)]
                ind_sol_col = pd.to_numeric(ind_df["solidity"], errors="coerce") if "solidity" in ind_df.columns else pd.Series(0.85, index=ind_df.index)
                ind_is_diss = ind_sol_col < 0.72
                if "pct_asymmetry" in ind_df.columns:
                    ind_asym = pd.to_numeric(ind_df["pct_asymmetry"], errors="coerce")
                    ind_df = ind_df[ind_asym.isna() | (ind_asym <= np.where(ind_is_diss, 0.42, 0.22))]
                    ind_is_diss = ind_is_diss.loc[ind_df.index]
                if "harmonic_roughness" in ind_df.columns:
                    ind_hr = pd.to_numeric(ind_df["harmonic_roughness"], errors="coerce")
                    ind_df = ind_df[ind_hr.isna() | (ind_hr <= np.where(ind_is_diss, 0.65, 0.22))]

                enrich_cols = [
                    c for c in (
                        "leaf_quality_weight",
                        "aspect_ratio",
                        "determiner_tier",
                        "species_standardized",
                        "is_calibrated_mm",
                        "pixels_per_mm",
                        "solidity",
                        "ucs_score",
                        "pct_asymmetry",
                        "harmonic_roughness",
                        "is_contour_outlier",
                    )
                    if c in ind_df.columns
                ]
                ind_sub = ind_df[["cat_str", "leaf_id_int"] + enrich_cols].drop_duplicates(
                    subset=["cat_str", "leaf_id_int"]
                )
                has_live_s2 = (
                    pd.to_numeric(df_plot["petiole_length_px"], errors="coerce").notna()
                    if "petiole_length_px" in df_plot.columns
                    else pd.Series(False, index=df_plot.index)
                )
                df_plot["_has_live_s2"] = has_live_s2
                df_plot = df_plot.merge(
                    ind_sub.assign(_in_ind_efa=True),
                    on=["cat_str", "leaf_id_int"],
                    how="left",
                    suffixes=("", "_efa"),
                )
                df_plot = (
                    df_plot[df_plot["_in_ind_efa"].fillna(False) | df_plot["_has_live_s2"].fillna(False)]
                    .drop(columns=["_in_ind_efa"])
                    .copy()
                )
                is_live_row = df_plot["_has_live_s2"].fillna(False)
                for col_e in enrich_cols:
                    efa_col = f"{col_e}_efa"
                    if efa_col in df_plot.columns:
                        if col_e not in df_plot.columns:
                            df_plot[col_e] = np.where(is_live_row, np.nan, df_plot[efa_col])
                        else:
                            df_plot[col_e] = np.where(
                                is_live_row,
                                df_plot[col_e],
                                df_plot[col_e].where(df_plot[col_e].notna(), df_plot[efa_col]),
                            )
                        df_plot.drop(columns=[efa_col], inplace=True)
                    elif is_live_row.any() and col_e in {"is_contour_outlier", "pct_asymmetry", "harmonic_roughness"}:
                        # Do not apply stale leaf_id QC flags from prior runs to newly segmented leaves
                        df_plot.loc[is_live_row, col_e] = np.nan
        except Exception as e_ind:
            logger.warning(f"Could not apply leaf_efa_harmonics_individual.csv filter: {e_ind}")

    if "_has_live_s2" not in df_plot.columns:
        df_plot["_has_live_s2"] = (
            pd.to_numeric(df_plot["petiole_length_px"], errors="coerce").notna()
            if "petiole_length_px" in df_plot.columns
            else False
        )

    # Also apply Step 03 leaf-level QC columns if directly present on df_plot
    if "is_contour_outlier" in df_plot.columns:
        df_plot = df_plot[~df_plot["is_contour_outlier"].apply(_parse_bool_val)].copy()
    plot_sol_s = pd.to_numeric(df_plot["solidity"], errors="coerce") if "solidity" in df_plot.columns else pd.Series(0.85, index=df_plot.index)
    plot_is_diss = plot_sol_s < 0.72
    if "pct_asymmetry" in df_plot.columns:
        asym_s = pd.to_numeric(df_plot["pct_asymmetry"], errors="coerce")
        df_plot = df_plot[asym_s.isna() | (asym_s <= np.where(plot_is_diss, 0.42, 0.22))].copy()
        plot_is_diss = plot_is_diss.loc[df_plot.index]
    if "harmonic_roughness" in df_plot.columns:
        hr_s = pd.to_numeric(df_plot["harmonic_roughness"], errors="coerce")
        df_plot = df_plot[hr_s.isna() | (hr_s <= np.where(plot_is_diss, 0.65, 0.22))].copy()

    # Filter to vouchers retained in leaf_efa_harmonics.csv when present
    otu_efa_path = _discover_downstream_table("leaf_efa_harmonics.csv", out_dir, active_cats)
    if otu_efa_path is not None and otu_efa_path.exists():
        try:
            otu_df = pd.read_csv(
                otu_efa_path,
                usecols=lambda c: c in {"catalogNumber", "is_contour_outlier"},
            )
            if "catalogNumber" in otu_df.columns:
                if "is_contour_outlier" in otu_df.columns:
                    otu_df = otu_df[~otu_df["is_contour_outlier"].apply(_parse_bool_val)]
                retained_cats = set(otu_df["catalogNumber"].astype(str).str.strip())
                df_plot = df_plot[df_plot["cat_str"].isin(retained_cats) | df_plot["_has_live_s2"].fillna(False)].copy()
        except Exception as e_otu:
            logger.warning(f"Could not apply leaf_efa_harmonics.csv filter: {e_otu}")

    # Filter out vouchers flagged as misidentified in Step 04 morphometric_flags.csv when present
    mflags_path = _discover_downstream_table("morphometric_flags.csv", out_dir, active_cats)
    misid_cats: Set[str] = set()
    valid_mflag_cats: Optional[Set[str]] = None
    if mflags_path is not None and mflags_path.exists():
        try:
            mf_df = pd.read_csv(
                mflags_path,
                usecols=lambda c: c in {"catalogNumber", "misidentification_flag"},
            ).copy()
            if "catalogNumber" in mf_df.columns:
                mf_df["cat_str"] = mf_df["catalogNumber"].astype(str).str.strip()
                if "misidentification_flag" in mf_df.columns:
                    is_misid = mf_df["misidentification_flag"].apply(_parse_bool_val)
                    misid_cats = set(mf_df.loc[is_misid, "cat_str"])
                    valid_mflag_cats = set(mf_df.loc[~is_misid, "cat_str"])
                else:
                    valid_mflag_cats = set(mf_df["cat_str"])
                df_plot = df_plot[
                    (df_plot["cat_str"].isin(valid_mflag_cats) & (~df_plot["cat_str"].isin(misid_cats)))
                    | df_plot["_has_live_s2"].fillna(False)
                ].copy()
        except Exception as e_mf:
            logger.warning(f"Could not apply morphometric_flags.csv filter: {e_mf}")
    df_plot.drop(columns=["_has_live_s2"], errors="ignore", inplace=True)

    if "misidentification_flag" in df_plot.columns:
        df_plot = df_plot[~df_plot["misidentification_flag"].apply(_parse_bool_val)].copy()

    if len(df_plot) == 0:
        return saved_paths

    if "species_standardized" in df_plot.columns and df_plot["species_standardized"].notna().any():
        df_plot["species_clean"] = df_plot["species_standardized"].apply(_standardize_packera_species_py)
    elif "scientificName" in df_plot.columns and df_plot["scientificName"].notna().any():
        df_plot["species_clean"] = df_plot["scientificName"].apply(_standardize_packera_species_py)
    elif "species_raw" in df_plot.columns and df_plot["species_raw"].notna().any():
        df_plot["species_clean"] = df_plot["species_raw"].apply(_standardize_packera_species_py)
    else:
        df_plot["species_clean"] = ""

    if cat_to_species:
        missing_sp = df_plot["species_clean"].isna() | df_plot["species_clean"].astype(str).str.strip().isin(
            {"", "Unassigned_Packera"}
        )
        if missing_sp.any():
            df_plot.loc[missing_sp, "species_clean"] = df_plot.loc[missing_sp, "cat_str"].map(
                lambda c: cat_to_species.get(c, "")
            )
        df_plot["species_clean"] = df_plot["species_clean"].replace({"": np.nan, "Unassigned_Packera": np.nan})
        if df_plot["species_clean"].isna().any():
            df_plot = df_plot[df_plot["species_clean"].notna()].copy()
    elif (df_plot["species_clean"].isna() | (df_plot["species_clean"].astype(str).str.strip() == "")).all():
        # Fallback label when synthetic unit-test fixture lacks voucher metadata
        df_plot["species_clean"] = "Unassigned_Packera"

    if len(df_plot) == 0:
        return saved_paths

    if "determiner_tier" not in df_plot.columns or df_plot["determiner_tier"].isna().all():
        df_plot["determiner_tier"] = df_plot["cat_str"].map(lambda c: cat_to_det_tier.get(c, "Tier_3_Bronze"))
    else:
        df_plot["determiner_tier"] = df_plot["determiner_tier"].fillna(
            df_plot["cat_str"].map(lambda c: cat_to_det_tier.get(c, "Tier_3_Bronze"))
        )

    df_plot["sol_val"] = (
        pd.to_numeric(df_plot["solidity"], errors="coerce").fillna(0.86)
        if "solidity" in df_plot.columns
        else pd.Series(0.86, index=df_plot.index)
    )
    ucs_col = "ucs_score" if "ucs_score" in df_plot.columns else ("ucs" if "ucs" in df_plot.columns else None)
    df_plot["ucs_val"] = (
        pd.to_numeric(df_plot[ucs_col], errors="coerce").fillna(0.90)
        if ucs_col is not None
        else pd.Series(0.90, index=df_plot.index)
    )

    if "aspect_ratio" in df_plot.columns:
        df_plot["ar_val"] = pd.to_numeric(df_plot["aspect_ratio"], errors="coerce")
    else:
        df_plot["ar_val"] = np.nan

    if df_plot["ar_val"].isna().any():
        # Pass 1: derive SVD aspect ratio from contour coordinates where available
        ar_pass1 = []
        for _, row_r in df_plot.iterrows():
            existing_ar = row_r.get("ar_val", np.nan)
            if pd.notna(existing_ar) and np.isfinite(float(existing_ar)) and float(existing_ar) > 0:
                ar_pass1.append(float(existing_ar))
                continue
            c_path_res = _resolve_asset_path(row_r.get("contour_path", ""), out_dir)
            ar_v = np.nan
            if c_path_res is not None and c_path_res.exists():
                try:
                    cdf = pd.read_csv(c_path_res)
                    if {"x", "y"}.issubset(cdf.columns) and len(cdf) >= 5:
                        pts = cdf[["x", "y"]].to_numpy(dtype=float)
                        pts = pts - pts.mean(axis=0)
                        _, s_vals, _ = np.linalg.svd(pts, full_matrices=False)
                        if len(s_vals) >= 2 and s_vals[1] > 1e-6:
                            ar_v = float(np.clip(s_vals[0] / s_vals[1], 1.05, 5.8))
                except Exception:
                    pass
            ar_pass1.append(ar_v)
        df_plot["ar_val"] = ar_pass1

        # Pass 2: impute remaining missing values from empirical species median (fallback global median or 2.15)
        if df_plot["ar_val"].isna().any():
            valid_ar_mask = df_plot["ar_val"].notna() & np.isfinite(df_plot["ar_val"]) & (df_plot["ar_val"] > 0)
            global_med_ar = float(df_plot.loc[valid_ar_mask, "ar_val"].median()) if valid_ar_mask.any() else 2.15
            sp_med_ar: Dict[str, float] = {}
            for sp_k in df_plot["species_clean"].dropna().unique():
                sp_vals = df_plot.loc[valid_ar_mask & (df_plot["species_clean"] == sp_k), "ar_val"]
                sp_med_ar[str(sp_k)] = float(sp_vals.median()) if len(sp_vals) > 0 else global_med_ar
            ar_final = []
            for _, row_r in df_plot.iterrows():
                cur_ar = row_r.get("ar_val", np.nan)
                if pd.notna(cur_ar) and np.isfinite(float(cur_ar)) and float(cur_ar) > 0:
                    ar_final.append(float(cur_ar))
                else:
                    sol_i = float(row_r["sol_val"])
                    base_ar = sp_med_ar.get(str(row_r["species_clean"]), global_med_ar)
                    ar_final.append(float(np.clip(base_ar + 1.8 * (0.90 - sol_i), 1.15, 5.2)))
            df_plot["ar_val"] = ar_final

    # Dynamically resolve active species present in df_plot and assign Okabe-Ito colors + hollow markers
    focal_std = _standardize_packera_species_py(
        os.environ.get("PACKERA_FOCAL_TAXON", "").strip() or "Packera dubia"
    )
    t1_counts = (
        df_plot[df_plot["norm_tier"] == "Tier_1_Direct"]["species_clean"]
        .dropna()
        .astype(str)
        .str.strip()
        .value_counts()
        .to_dict()
    )
    unique_sp = [str(sp).strip() for sp in df_plot["species_clean"].dropna().unique() if str(sp).strip()]
    present_sp = sorted(
        unique_sp,
        key=lambda sp: (0 if sp == focal_std else 1, -int(t1_counts.get(sp, 0)), sp),
    )
    if not present_sp:
        present_sp = ["Unassigned_Packera"]
    for idx_sp, sp_name in enumerate(present_sp):
        SPECIES_PALETTE[sp_name] = EXTENDED_OKABE_ITO[idx_sp % len(EXTENDED_OKABE_ITO)]
        SPECIES_MARKERS[sp_name] = HOLLOW_MARKER_CYCLE[idx_sp % len(HOLLOW_MARKER_CYCLE)]
        if sp_name.startswith("Packera "):
            ep_tex = sp_name[len("Packera "):].strip().replace("_", " ").replace(" ", r"\ ")
            SPECIES_ITALIC[sp_name] = (
                rf"$\it{{Packera\ {ep_tex}}}$" if idx_sp == 0 else rf"$\it{{P.\ {ep_tex}}}$"
            )
        elif re.match(r"^[A-Z][a-z]+\s+[a-z-]+$", sp_name):
            sp_tex = sp_name.replace(" ", r"\ ")
            SPECIES_ITALIC[sp_name] = rf"$\it{{{sp_tex}}}$"
        else:
            SPECIES_ITALIC[sp_name] = sp_name.replace("_", " ")

    for c_val, sp_val in zip(df_plot["cat_str"], df_plot["species_clean"]):
        if c_val and pd.notna(sp_val) and str(sp_val).strip() and c_val not in cat_to_species:
            cat_to_species[c_val] = str(sp_val).strip()

    if "leaf_zone" in df_plot.columns:
        df_plot["leaf_zone_norm"] = (
            df_plot["leaf_zone"]
            .fillna("basal")
            .astype(str)
            .str.strip()
            .str.lower()
            .map(lambda z: "cauline" if z == "cauline" else "basal")
        )
    else:
        df_plot["leaf_zone_norm"] = "basal"

    df_basal = df_plot[df_plot["leaf_zone_norm"] == "basal"].copy()
    df_cauline = df_plot[df_plot["leaf_zone_norm"] == "cauline"].copy()
    if len(df_basal) == 0:
        df_basal = df_plot.copy()
    has_cauline_zone = len(df_cauline) > 0

    alpha_pt = 0.45 if len(df_plot) > 1000 else 0.75

    # -------------------------------------------------------------------------
    # Figure 02a: Lamina Solidity vs. Aspect Ratio (Separated Basal vs. Cauline)
    #             & Taxon Reproductive Traits
    # -------------------------------------------------------------------------
    base_fig_h = 7.4 if has_cauline_zone else 5.8
    fig_a_height = max(base_fig_h, 0.34 * min(16, len(present_sp)) + (3.2 if has_cauline_zone else 2.4))
    fig_a = plt.figure(figsize=(13.0, fig_a_height), dpi=300)
    gs_a = GridSpec(1, 2, figure=fig_a, width_ratios=[1.15, 1.05], wspace=0.28)

    x_min = max(0.15, min(0.48, float(df_plot["sol_val"].min()) - 0.04))
    y_max = min(6.0, max(4.2, float(df_plot["ar_val"].max()) + 0.35))
    fg_sp = present_sp if len(present_sp) <= 6 else present_sp[:6]
    bg_sp = [] if len(present_sp) <= 6 else present_sp[6:]

    env_focal = os.environ.get("PACKERA_FOCAL_TAXON", "").strip()
    title_focal_std = _standardize_packera_species_py(env_focal) if env_focal else (present_sp[0] if len(present_sp) == 1 else "")

    def _draw_solidity_vs_ar_zone_panel(
        ax_sc: Any,
        zone_df: pd.DataFrame,
        panel_prefix: str,
        zone_title_lbl: str,
        show_xlabel: bool = True,
    ) -> None:
        ax_sc.add_patch(
            Rectangle(
                (0.72, 1.0), 1.02 - 0.72, y_max - 1.0,
                facecolor="#009E73", alpha=0.06, edgecolor="none", zorder=0,
                label="Tier 1 Intact Zone (Solidity >= 0.72)"
            )
        )
        ax_sc.add_patch(
            Rectangle(
                (0.48, 1.0), 0.72 - 0.48, y_max - 1.0,
                facecolor="#CC79A7", alpha=0.08, edgecolor="none", zorder=0,
                label="Pinnatifid Dissection Zone (0.48-0.72)"
            )
        )

        dense_zone = len(zone_df) > 1500
        fg_s = 18 if dense_zone else 38
        fg_alpha = 0.35 if dense_zone else alpha_pt
        fg_lw = 0.75 if dense_zone else 1.0
        bg_s = 12 if dense_zone else 22
        bg_alpha = 0.18 if dense_zone else 0.30

        if bg_sp:
            bg_sub = zone_df[zone_df["species_clean"].isin(bg_sp)]
            if len(bg_sub) > 0:
                ax_sc.scatter(
                    bg_sub["sol_val"],
                    bg_sub["ar_val"],
                    facecolors="none",
                    edgecolors="#B0B0B0",
                    marker="o",
                    label=f"Other Packera ({len(bg_sp)} spp., n={len(bg_sub)})",
                    alpha=bg_alpha,
                    s=bg_s,
                    linewidths=0.7,
                    zorder=2,
                )

        for sp_name in fg_sp:
            sub = zone_df[zone_df["species_clean"] == sp_name]
            if len(sub) > 0:
                col_hex = SPECIES_PALETTE[sp_name]
                mk_sym = SPECIES_MARKERS[sp_name]
                ax_sc.scatter(
                    sub["sol_val"],
                    sub["ar_val"],
                    facecolors="none",
                    edgecolors=col_hex,
                    marker=mk_sym,
                    label=f"{SPECIES_ITALIC[sp_name]} (n={len(sub)})",
                    alpha=fg_alpha,
                    s=fg_s,
                    linewidths=fg_lw,
                    zorder=4,
                )
                sub_xy = sub[["sol_val", "ar_val"]].dropna()
                if len(sub_xy) >= 10:
                    xv = sub_xy["sol_val"].to_numpy(dtype=float)
                    yv = sub_xy["ar_val"].to_numpy(dtype=float)
                    valid_m = np.isfinite(xv) & np.isfinite(yv)
                    xv, yv = xv[valid_m], yv[valid_m]
                    if len(xv) >= 10:
                        mx = float(np.mean(xv))
                        my = float(np.mean(yv))
                        cov_mat = np.cov(np.vstack([xv, yv]))
                        if cov_mat.shape == (2, 2) and np.all(np.isfinite(cov_mat)):
                            eigvals, eigvecs = np.linalg.eigh(cov_mat)
                            eigvals = np.maximum(eigvals, 1e-6)
                            chi2_val = 2.7726
                            theta_t = np.linspace(0.0, 2.0 * np.pi, 64)
                            unit_circle = np.vstack([np.cos(theta_t), np.sin(theta_t)])
                            ell_pts = eigvecs @ (np.sqrt(chi2_val * eigvals)[:, None] * unit_circle)
                            ex = mx + ell_pts[0, :]
                            ey = my + ell_pts[1, :]
                            ax_sc.plot(ex, ey, color="white", linewidth=2.8, zorder=5)
                            ax_sc.plot(ex, ey, color=col_hex, linewidth=1.6, zorder=6)
                            ax_sc.scatter(
                                [mx],
                                [my],
                                s=58,
                                facecolors=col_hex,
                                edgecolors="#111111",
                                linewidths=1.1,
                                marker=mk_sym,
                                zorder=7,
                            )

        if title_focal_std and title_focal_std != "Unassigned_Packera":
            sp_title = title_focal_std.replace("_", " ").replace(" ", r"\ ")
            panel_title = rf"{panel_prefix}. {zone_title_lbl} Solidity vs. Blade Elongation ($\it{{{sp_title}}}$)"
        elif len(present_sp) == 1 and present_sp[0] != "Unassigned_Packera":
            sp_title = present_sp[0].replace("_", " ").replace(" ", r"\ ")
            panel_title = rf"{panel_prefix}. {zone_title_lbl} Solidity vs. Blade Elongation ($\it{{{sp_title}}}$)"
        else:
            panel_title = f"{panel_prefix}. {zone_title_lbl} Solidity vs. Blade Elongation ($N = {len(present_sp)}$ Taxa)"

        ax_sc.axvline(0.72, color="#222222", linestyle="--", linewidth=1.3)
        ax_sc.axvline(0.48, color="#CC79A7", linestyle=":", linewidth=1.2)
        if show_xlabel:
            ax_sc.set_xlabel("Lamina Solidity (Area / Convex Hull Area)", fontsize=9.8, fontweight="bold")
        ax_sc.set_ylabel("Blade Aspect Ratio (L/W)", fontsize=9.8, fontweight="bold")
        ax_sc.set_title(
            panel_title,
            fontsize=10.2 if has_cauline_zone else 10.8,
            fontweight="bold",
            loc="left",
        )
        ax_sc.set_xlim(x_min, 1.02)
        ax_sc.set_ylim(1.0, y_max)
        ax_sc.grid(True, linestyle=":", alpha=0.45)
        ax_sc.legend(loc="upper left", fontsize=7.4 if has_cauline_zone else 7.8, frameon=True, framealpha=0.90)

    if has_cauline_zone:
        gs_left = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_a[0, 0], height_ratios=[1.0, 1.0], hspace=0.42)
        ax_sc_bas = fig_a.add_subplot(gs_left[0, 0])
        _draw_solidity_vs_ar_zone_panel(
            ax_sc_bas, df_basal, panel_prefix="A1", zone_title_lbl="Basal Rosette", show_xlabel=False
        )
        ax_sc_caul = fig_a.add_subplot(gs_left[1, 0])
        _draw_solidity_vs_ar_zone_panel(
            ax_sc_caul, df_cauline, panel_prefix="A2", zone_title_lbl="Cauline Stem", show_xlabel=True
        )
    else:
        ax_sc = fig_a.add_subplot(gs_a[0, 0])
        _draw_solidity_vs_ar_zone_panel(
            ax_sc, df_basal, panel_prefix="A", zone_title_lbl="Basal Rosette", show_xlabel=True
        )

    # Panel B: Taxon-Stratified Reproductive Trait Distributions (B1: Involucre AR, B2: Phenological State)
    gs_b = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_a[0, 1], height_ratios=[1.15, 1.0], hspace=0.52)
    ax_repro = fig_a.add_subplot(gs_b[0, 0])

    qc_voucher_cats = set(df_plot["cat_str"].unique())
    rep_df = repro_df.copy() if (repro_df is not None and len(repro_df) > 0) else pd.DataFrame()
    if len(rep_df) > 0 and "catalogNumber" in rep_df.columns:
        rep_df["cat_str"] = rep_df["catalogNumber"].astype(str).str.strip()
        if misid_cats:
            rep_df = rep_df[~rep_df["cat_str"].isin(misid_cats)].copy()
        if (rep_df["cat_str"].isin(qc_voucher_cats)).any():
            rep_df = rep_df[rep_df["cat_str"].isin(qc_voucher_cats)].copy()
        if "species_standardized" in rep_df.columns:
            rep_df["sp_clean"] = rep_df["species_standardized"].apply(_standardize_packera_species_py)
        elif cat_to_species:
            rep_df["sp_clean"] = rep_df["cat_str"].map(lambda c: cat_to_species.get(c, "Unassigned_Packera"))
        else:
            rep_df["sp_clean"] = "Unassigned_Packera"
        ar_col = (
            "capitulum_aspect_ratio"
            if "capitulum_aspect_ratio" in rep_df.columns
            else ("AR_inv" if "AR_inv" in rep_df.columns else "involucre_width_mm")
        )
        rep_df["trait_v"] = (
            pd.to_numeric(rep_df[ar_col], errors="coerce")
            if ar_col in rep_df.columns
            else pd.Series(np.nan, index=rep_df.index)
        )
        rep_valid = rep_df[rep_df["trait_v"].notna() & (rep_df["trait_v"] > 0)].copy()
    else:
        rep_valid = pd.DataFrame()

    horiz_repro = len(present_sp) > 5
    rng_jit = np.random.default_rng(42)
    box_vals: List[np.ndarray] = []
    box_positions: List[int] = []
    box_species: List[str] = []
    empty_positions: List[int] = []

    for idx_s, sp_name in enumerate(present_sp):
        pos_i = idx_s if horiz_repro else (idx_s + 1)
        if len(rep_valid) > 0 and (rep_valid["sp_clean"] == sp_name).any():
            vals = rep_valid.loc[rep_valid["sp_clean"] == sp_name, "trait_v"].dropna().to_numpy(dtype=float)
        else:
            vals = np.array([], dtype=float)

        if len(vals) > 0:
            box_vals.append(vals)
            box_positions.append(pos_i)
            box_species.append(sp_name)
            jit_offset = rng_jit.uniform(-0.16, 0.16, size=len(vals))
            if horiz_repro:
                ax_repro.scatter(
                    vals,
                    pos_i + jit_offset,
                    facecolors="none",
                    edgecolors=SPECIES_PALETTE[sp_name],
                    marker=SPECIES_MARKERS[sp_name],
                    alpha=0.22,
                    s=11,
                    linewidths=0.75,
                    zorder=2,
                )
            else:
                ax_repro.scatter(
                    pos_i + jit_offset,
                    vals,
                    facecolors="none",
                    edgecolors=SPECIES_PALETTE[sp_name],
                    marker=SPECIES_MARKERS[sp_name],
                    alpha=0.22,
                    s=11,
                    linewidths=0.75,
                    zorder=2,
                )
        else:
            empty_positions.append(pos_i)

    if len(box_vals) > 0:
        try:
            bp = ax_repro.boxplot(
                box_vals,
                positions=box_positions,
                orientation="horizontal" if horiz_repro else "vertical",
                widths=0.48,
                patch_artist=True,
                showfliers=False,
                zorder=5,
            )
        except TypeError:
            bp = ax_repro.boxplot(
                box_vals,
                positions=box_positions,
                vert=not horiz_repro,
                widths=0.48,
                patch_artist=True,
                showfliers=False,
                zorder=5,
            )
        for patch, med, sp_name in zip(bp["boxes"], bp["medians"], box_species):
            c_hex = SPECIES_PALETTE[sp_name]
            patch.set_facecolor("#F2F2F2")
            patch.set_alpha(0.72)
            patch.set_edgecolor(c_hex)
            patch.set_linewidth(1.4)
            patch.set_zorder(5)
            med.set_color("#111111")
            med.set_linewidth(2.2)
            med.set_zorder(6)
        for w_line in bp.get("whiskers", []) + bp.get("caps", []):
            w_line.set_zorder(5)
        all_repro_concat = np.concatenate(box_vals)
        baseline_val = max(0.45, float(np.nanmin(all_repro_concat)) * 0.92)
    else:
        baseline_val = 0.60
        if horiz_repro:
            ax_repro.set_xlim(0.5, 2.0)
        else:
            ax_repro.set_ylim(0.5, 2.0)

    for pos_i in empty_positions:
        if horiz_repro:
            ax_repro.text(
                baseline_val,
                pos_i,
                "n=0",
                color="#666666",
                style="italic",
                fontsize=7.8,
                va="center",
                ha="left",
            )
        else:
            ax_repro.text(
                pos_i,
                baseline_val,
                "n=0",
                color="#666666",
                style="italic",
                fontsize=7.8,
                ha="center",
                va="bottom",
            )

    if horiz_repro:
        y_ticks_repro = np.arange(len(present_sp))
        ax_repro.set_yticks(y_ticks_repro)
        ax_repro.set_yticklabels([SPECIES_ITALIC[sp] for sp in present_sp], fontsize=8.2)
        ax_repro.set_ylim(-0.6, len(present_sp) - 0.4)
        ax_repro.set_xlabel("Involucre Aspect Ratio (L/W)", fontsize=9.2, fontweight="bold")
        ax_repro.grid(True, axis="x", linestyle=":", alpha=0.45)
    else:
        x_ticks_repro = list(range(1, len(present_sp) + 1))
        ax_repro.set_xticks(x_ticks_repro)
        ax_repro.set_xticklabels([SPECIES_ITALIC[sp] for sp in present_sp], fontsize=8.5)
        ax_repro.set_xlim(0.4, len(present_sp) + 0.6)
        ax_repro.set_ylabel("Involucre Aspect Ratio (L/W)", fontsize=9.2, fontweight="bold")
        ax_repro.grid(True, axis="y", linestyle=":", alpha=0.45)

    ax_repro.set_title("B1. Capitulum / Involucre Proportions by Species", fontsize=10.2, fontweight="bold", loc="left")

    # Panel B2: Anthesis Phenological State Proportions across Species
    ax_pheno = fig_a.add_subplot(gs_b[1, 0])
    ph_df = pheno_df.copy() if (pheno_df is not None and len(pheno_df) > 0) else pd.DataFrame()
    if len(ph_df) > 0 and "catalogNumber" in ph_df.columns:
        ph_df["cat_str"] = ph_df["catalogNumber"].astype(str).str.strip()
        if misid_cats:
            ph_df = ph_df[~ph_df["cat_str"].isin(misid_cats)].copy()
        if (ph_df["cat_str"].isin(qc_voucher_cats)).any():
            ph_df = ph_df[ph_df["cat_str"].isin(qc_voucher_cats)].copy()
        if "species_standardized" in ph_df.columns:
            ph_df["sp_clean"] = ph_df["species_standardized"].apply(_standardize_packera_species_py)
        elif cat_to_species:
            ph_df["sp_clean"] = ph_df["cat_str"].map(lambda c: cat_to_species.get(c, "Unassigned_Packera"))
        else:
            ph_df["sp_clean"] = "Unassigned_Packera"
        st_col = (
            "voucher_phenological_state"
            if "voucher_phenological_state" in ph_df.columns
            else ("phenological_state" if "phenological_state" in ph_df.columns else ph_df.columns[-1])
        )
        ph_df["state_str"] = ph_df[st_col].dropna().astype(str)
        ph_df = ph_df[ph_df["state_str"].str.strip() != ""].copy()
    else:
        ph_df = pd.DataFrame(columns=["sp_clean", "state_str"])

    anthesis_pct: List[float] = []
    bud_pct: List[float] = []
    fruiting_pct: List[float] = []
    no_pheno_indices: List[int] = []
    for idx_p, sp_name in enumerate(present_sp):
        sub_p = ph_df[ph_df["sp_clean"] == sp_name] if len(ph_df) > 0 else pd.DataFrame()
        if len(sub_p) > 0:
            s_lower = sub_p["state_str"].astype(str)
            is_anth = s_lower.str.contains("anthesis|flower|open", case=False, na=False) & ~s_lower.str.contains(
                "pre|bud|post|fruit|senes", case=False, na=False
            )
            is_bud = s_lower.str.contains("bud|pre", case=False, na=False) & ~is_anth
            is_fruit = ~is_anth & ~is_bud
            p_a = float(is_anth.mean() * 100.0)
            p_b = float(is_bud.mean() * 100.0)
            p_f = float(is_fruit.mean() * 100.0)
            anthesis_pct.append(p_a)
            bud_pct.append(p_b)
            fruiting_pct.append(p_f)
        else:
            anthesis_pct.append(0.0)
            bud_pct.append(0.0)
            fruiting_pct.append(0.0)
            no_pheno_indices.append(idx_p)

    y_pos = np.arange(len(present_sp))
    left_bud = np.array(anthesis_pct, dtype=float)
    left_fruit = left_bud + np.array(bud_pct, dtype=float)
    ax_pheno.barh(y_pos, anthesis_pct, color="#009E73", edgecolor="white", height=0.58, alpha=0.88, label="Verified Anthesis")
    ax_pheno.barh(y_pos, bud_pct, left=left_bud, color="#56B4E9", edgecolor="white", height=0.58, alpha=0.78, label="Pre-Anthesis Bud")
    ax_pheno.barh(y_pos, fruiting_pct, left=left_fruit, color="#E69F00", edgecolor="white", height=0.58, alpha=0.78, label="Fruiting / Post-Anthesis")
    for idx_p in no_pheno_indices:
        ax_pheno.text(
            2.0,
            idx_p,
            "No pheno data",
            color="#666666",
            style="italic",
            fontsize=7.8,
            va="center",
            ha="left",
        )
    ax_pheno.set_yticks(y_pos)
    ax_pheno.set_yticklabels([SPECIES_ITALIC[sp] for sp in present_sp], fontsize=8.2 if horiz_repro else 8.5)
    ax_pheno.set_xlim(0, 105)
    ax_pheno.set_xlabel("Proportion of Vouchers (%)", fontsize=9.2, fontweight="bold")
    ax_pheno.set_title("B2. Optical Phenology State Proportions", fontsize=9.8, fontweight="bold", loc="left", pad=22)
    ax_pheno.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=3, frameon=False, fontsize=7.6)
    ax_pheno.grid(True, axis="x", linestyle=":", alpha=0.45)

    fig02a_pdf, fig02a_png = _save_fig_fn(
        fig_a, out_dir / "fig02a_segmentation_qc_distributions", dpi=300
    )
    plt.close(fig_a)
    saved_paths["fig02a_pdf"] = fig02a_pdf
    saved_paths["fig02a_png"] = fig02a_png

    # -------------------------------------------------------------------------
    # Figure 02b: Taxon-Stratified Basal Leaf Silhouette Montage
    #             (Rows = Active Species; Cols = Type Anchor / 20th, 35th/40th, 60th, 80th/85th %ile)
    # -------------------------------------------------------------------------
    pct_labels = [
        "20th %ile (Compact)",
        "40th %ile (Median-Low)",
        "60th %ile (Median-High)",
        "80th %ile (Elongate)",
    ]
    default_quantiles = (0.20, 0.40, 0.60, 0.80)
    type_anchored_quantiles = (0.35, 0.60, 0.85)

    t1_pool = df_plot[df_plot["norm_tier"] == "Tier_1_Direct"].copy()
    if len(t1_pool) == 0:
        t1_pool = df_plot.copy()

    # Ranking attributes on t1_pool (mask smoothness/closure validity is evaluated lazily per species)
    is_type_s = (
        t1_pool["is_type_specimen"].apply(_parse_bool_val)
        if "is_type_specimen" in t1_pool.columns
        else pd.Series(False, index=t1_pool.index)
    )
    det_type_s = (
        t1_pool["determiner_tier"].astype(str).str.strip().eq("Tier_0_Type_Anchor")
        if "determiner_tier" in t1_pool.columns
        else pd.Series(False, index=t1_pool.index)
    )
    t1_pool["is_type_bool"] = is_type_s | det_type_s
    t1_pool["is_gold_silver"] = (
        t1_pool["determiner_tier"].astype(str).str.strip().isin({"Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver"})
        | t1_pool["is_type_bool"]
    )
    if "is_calibrated_mm" in t1_pool.columns:
        t1_pool["is_cal_bool"] = t1_pool["is_calibrated_mm"].apply(_parse_bool_val)
    elif "pixels_per_mm" in t1_pool.columns:
        ppm_s = pd.to_numeric(t1_pool["pixels_per_mm"], errors="coerce")
        t1_pool["is_cal_bool"] = ppm_s.notna() & (ppm_s > 0)
    else:
        t1_pool["is_cal_bool"] = False

    lqw_num = (
        pd.to_numeric(t1_pool["leaf_quality_weight"], errors="coerce")
        if "leaf_quality_weight" in t1_pool.columns
        else pd.Series(np.nan, index=t1_pool.index)
    )
    fallback_qw = t1_pool["sol_val"] * t1_pool["ucs_val"]
    t1_pool["qual_weight"] = lqw_num.where(lqw_num.notna() & (lqw_num > 0), fallback_qw)

    t1_basal_pool = t1_pool[t1_pool["leaf_zone_norm"] == "basal"].copy()
    if len(t1_basal_pool) == 0:
        t1_basal_pool = t1_pool.copy()
    t1_cauline_pool = t1_pool[t1_pool["leaf_zone_norm"] == "cauline"].copy()

    def _render_montage_page(
        montage_sp_list: List[str],
        out_stem: str,
        title_suffix: str = "",
        pool_df: Optional[pd.DataFrame] = None,
        zone_label: str = "Basal",
        render_mode: str = "silhouette",
    ) -> Tuple[Path, Path]:
        active_pool = pool_df if pool_df is not None else t1_basal_pool
        n_sp_rows = max(1, len(montage_sp_list))
        fig_height = min(24.0, max(3.6, 3.0 * n_sp_rows))
        fig_b, axes_b = plt.subplots(
            n_sp_rows, 4, figsize=(12.2, fig_height), dpi=300, squeeze=False
        )

        for row_idx, sp_name in enumerate(montage_sp_list):
            sp_col = SPECIES_PALETTE[sp_name]
            # Strictly restrict to this species and leaf_zone pool; never substitute another species' or zone's leaves
            sp_all = active_pool[active_pool["species_clean"] == sp_name].copy()
            # Only real on-disk silhouettes are eligible (never render fabricated geometry)
            sp_all["asset_ok"] = pd.Series(
                [
                    (_resolve_asset_path(r.get("mask_path", ""), out_dir) is not None)
                    or (_resolve_asset_path(r.get("contour_path", ""), out_dir) is not None)
                    for _, r in sp_all.iterrows()
                ],
                index=sp_all.index,
                dtype=bool,
            )
            n_missing_assets = int((~sp_all["asset_ok"]).sum())
            if n_missing_assets:
                logger.warning(
                    f"fig02b: {n_missing_assets} / {len(sp_all)} {sp_name} Tier_1_Direct leaves have no resolvable mask/contour asset on disk."
                )
            sp_all = sp_all[sp_all["asset_ok"]].copy()
            if "area_px" in sp_all.columns:
                sp_all["area_val"] = pd.to_numeric(sp_all["area_px"], errors="coerce").fillna(1.0)
            elif "blade_area_cm2" in sp_all.columns:
                sp_all["area_val"] = pd.to_numeric(sp_all["blade_area_cm2"], errors="coerce").fillna(1.0)
            else:
                sp_all["area_val"] = 1.0

            is_cauline_zone = str(zone_label).strip().lower() == "cauline"
            if is_cauline_zone and len(sp_all) > 0:
                caul_conf = (
                    pd.to_numeric(sp_all["leaf_zone_confidence"], errors="coerce").fillna(1.0)
                    if "leaf_zone_confidence" in sp_all.columns
                    else pd.Series(1.0, index=sp_all.index)
                )
                caul_dist = (
                    pd.to_numeric(sp_all["dist_from_crown_norm"], errors="coerce").fillna(1.0)
                    if "dist_from_crown_norm" in sp_all.columns
                    else pd.Series(1.0, index=sp_all.index)
                )
                caul_pet_r = (
                    pd.to_numeric(sp_all["petiole_to_blade_ratio"], errors="coerce").fillna(1.0)
                    if "petiole_to_blade_ratio" in sp_all.columns
                    else pd.Series(0.0, index=sp_all.index)
                )
                sp_caul_strict = sp_all[
                    (caul_conf >= 0.62)
                    & (caul_dist >= 0.35)
                    & (
                        (caul_pet_r <= 0.35)
                        | ((sp_all["ar_val"] >= 3.0) & (sp_all["sol_val"] <= 0.83) & (caul_pet_r <= 0.75))
                        | ((sp_all["sol_val"] <= 0.72) & (caul_pet_r <= 1.35))
                    )
                    & ((sp_all["ar_val"] >= 1.65) | (sp_all["sol_val"] <= 0.72))
                    & (sp_all["sol_val"] <= 0.88)
                ].copy()
                if "petiole_length_px" in sp_caul_strict.columns:
                    has_live_pet = pd.to_numeric(sp_caul_strict["petiole_length_px"], errors="coerce").notna()
                    if int(has_live_pet.sum()) >= 4:
                        sp_caul_strict = sp_caul_strict[has_live_pet].copy()
                if len(sp_caul_strict) > 0:
                    sp_all = sp_caul_strict

            sp_med_area = float(sp_all["area_val"].median()) if len(sp_all) > 0 else 1.0
            if len(sp_all) >= 8 and sp_med_area > 0:
                sp_size_ok = sp_all[sp_all["area_val"] >= (0.20 if is_cauline_zone else 0.35) * sp_med_area]
                if len(sp_size_ok) >= 4:
                    sp_all = sp_size_ok.copy()

            is_diss_col = (
                sp_all["is_dissected"].astype(str).str.lower().isin(["true", "1", "t", "yes"])
                if "is_dissected" in sp_all.columns
                else (sp_all["sol_val"] < 0.72)
            )
            sp_all["eff_ucs"] = np.where(
                is_diss_col | (sp_all["sol_val"] < 0.72),
                np.minimum(
                    1.0,
                    np.maximum(
                        sp_all["ucs_val"] / np.maximum(sp_all["sol_val"], 0.35),
                        sp_all["ucs_val"],
                    ),
                ),
                sp_all["ucs_val"],
            )
            sp_diss_frac = float((is_diss_col | (sp_all["sol_val"] < 0.68)).mean()) if len(sp_all) > 0 else 0.0
            sp_med_sol = float(sp_all["sol_val"].median()) if len(sp_all) > 0 else 0.85
            is_dissected_cohort = bool(sp_med_sol < 0.76 or sp_diss_frac >= 0.30)

            if is_dissected_cohort:
                # Dissected / pinnatisect species: prioritize authentic dissected blades over undivided terminal lobes
                sp_hq = sp_all[
                    (sp_all["eff_ucs"] >= 0.78)
                    & (sp_all["sol_val"] >= 0.22)
                    & (sp_all["sol_val"] < (0.78 if is_cauline_zone else 0.65))
                    & (sp_all["area_val"] >= (0.35 if is_cauline_zone else 0.55) * sp_med_area)
                    & (is_diss_col | (sp_all["sol_val"] < (0.78 if is_cauline_zone else 0.65)))
                ].copy()
                if len(sp_hq) >= 4:
                    sp_hq["mask_valid"] = [
                        _is_valid_smooth_closed_mask(r.get("mask_path", ""), out_dir) for _, r in sp_hq.iterrows()
                    ]
                    sp_hq = sp_hq[sp_hq["mask_valid"].astype(bool)].sort_values("ar_val").copy()
                if len(sp_hq) < 4:
                    sp_hq = sp_all[
                        (sp_all["eff_ucs"] >= 0.75)
                        & (sp_all["sol_val"] >= 0.20)
                        & (sp_all["sol_val"] < (0.80 if is_cauline_zone else 0.68))
                        & (is_diss_col | (sp_all["sol_val"] < (0.80 if is_cauline_zone else 0.68)))
                    ].copy()
                    sp_hq["mask_valid"] = [
                        _is_valid_smooth_closed_mask(r.get("mask_path", ""), out_dir) for _, r in sp_hq.iterrows()
                    ]
                    sp_hq = sp_hq[sp_hq["mask_valid"].astype(bool)].sort_values("ar_val").copy()
            elif is_cauline_zone:
                sp_hq = sp_all[
                    (
                        (sp_all["eff_ucs"] >= 0.76)
                        & (sp_all["sol_val"] >= 0.20)
                        & (sp_all["sol_val"] < 0.82)
                    )
                    | (
                        (sp_all["ucs_val"] >= 0.84)
                        & (sp_all["ar_val"] >= 1.85)
                        & (sp_all["sol_val"] < 0.995)
                    )
                ].copy()
                if len(sp_hq) > 0:
                    sp_hq["mask_valid"] = [
                        _is_valid_smooth_closed_mask(r.get("mask_path", ""), out_dir) for _, r in sp_hq.iterrows()
                    ]
                    sp_hq = sp_hq[sp_hq["mask_valid"].astype(bool)].sort_values("ar_val").copy()
            else:
                sp_hq = pd.DataFrame()

            if len(sp_hq) < 4:
                sp_hq = sp_all[
                    (
                        (sp_all["ucs_val"] >= (0.82 if is_cauline_zone else 0.88))
                        & (sp_all["sol_val"] >= (0.20 if is_cauline_zone else 0.80))
                        & (sp_all["sol_val"] < 0.995)
                    )
                    | (
                        (is_dissected_cohort or is_cauline_zone)
                        & (sp_all["eff_ucs"] >= 0.76)
                        & (sp_all["sol_val"] >= 0.18)
                        & (is_diss_col | (sp_all["sol_val"] < 0.80))
                    )
                ].copy()
                sp_hq["mask_valid"] = [
                    _is_valid_smooth_closed_mask(r.get("mask_path", ""), out_dir) for _, r in sp_hq.iterrows()
                ]
                sp_hq = sp_hq[sp_hq["mask_valid"].astype(bool)].sort_values("ar_val").copy()
            if len(sp_hq) < 4:
                # Fallback within the SAME species only (e.g. small synthetic unit test fixture)
                sp_all["mask_valid"] = [
                    _is_valid_smooth_closed_mask(r.get("mask_path", ""), out_dir) for _, r in sp_all.iterrows()
                ]
                sp_hq = sp_all[sp_all["mask_valid"].astype(bool)].sort_values("ar_val").copy()
            if len(sp_hq) == 0:
                sp_hq = sp_all.sort_values("ar_val").copy()

            sp_target_sol = float(sp_hq["sol_val"].median()) if len(sp_hq) > 0 else sp_med_sol
            sp_hq_med_area = float(sp_hq["area_val"].median()) if len(sp_hq) > 0 else max(sp_med_area, 1.0)
            if is_dissected_cohort:
                sp_hq["rep_score"] = (
                    (1.0 - np.abs(sp_hq["sol_val"] - sp_target_sol))
                    * np.minimum(1.15, sp_hq["area_val"] / max(sp_hq_med_area, 1.0))
                    * sp_hq["eff_ucs"]
                )
            else:
                sp_hq["rep_score"] = sp_hq["qual_weight"] * np.minimum(
                    1.10, np.maximum(0.85, sp_hq["area_val"] / max(sp_hq_med_area, 1.0))
                )

            n_avail = len(sp_hq)
            if n_avail > 1:
                sp_hq["q_rank"] = np.linspace(0.0, 1.0, n_avail)
            elif n_avail == 1:
                sp_hq["q_rank"] = [0.5]

            selected_rows: List[Optional[pd.Series]] = []
            used_keys: Set[Tuple[str, int]] = set()
            used_row_ids: Set[Any] = set()

            type_candidates = (
                sp_hq[sp_hq["is_type_bool"].astype(bool)]
                if ("is_type_bool" in sp_hq.columns and not is_cauline_zone)
                else pd.DataFrame()
            )
            if len(type_candidates) == 0 and not is_dissected_cohort and not is_cauline_zone and "is_type_bool" in sp_all.columns:
                type_candidates = sp_all[sp_all["is_type_bool"].astype(bool)]

            if len(type_candidates) > 0:
                best_type_row = max(
                    [r for _, r in type_candidates.iterrows()],
                    key=lambda r: (
                        int(bool(r.get("is_cal_bool", False))),
                        float(r.get("rep_score", r.get("qual_weight", 0.0))),
                        float(r.get("qual_weight", 0.0)),
                    ),
                )
                used_keys.add((str(best_type_row.get("cat_str", "")), int(best_type_row.get("leaf_id_int", 1))))
                used_row_ids.add(best_type_row.name)
                selected_rows.append(best_type_row)
                row_quantiles: Tuple[float, ...] = type_anchored_quantiles
            else:
                row_quantiles = default_quantiles

            for q_target in row_quantiles:
                if n_avail == 0:
                    selected_rows.append(None)
                    continue
                # Primary window: q +/- 0.08
                win_df = sp_hq[(sp_hq["q_rank"] >= q_target - 0.08 - 1e-9) & (sp_hq["q_rank"] <= q_target + 0.08 + 1e-9)]
                win_unused = [
                    r for r_idx, r in win_df.iterrows()
                    if r_idx not in used_row_ids
                    and (str(r.get("cat_str", "")), int(r.get("leaf_id_int", 1))) not in used_keys
                ]
                if not win_unused:
                    # Secondary wider window within the same species if n_avail is small
                    win_df2 = sp_hq[(sp_hq["q_rank"] >= q_target - 0.18 - 1e-9) & (sp_hq["q_rank"] <= q_target + 0.18 + 1e-9)]
                    win_unused = [
                        r for r_idx, r in win_df2.iterrows()
                        if r_idx not in used_row_ids
                        and (str(r.get("cat_str", "")), int(r.get("leaf_id_int", 1))) not in used_keys
                    ]
                if not win_unused:
                    win_unused = [
                        sp_hq.loc[r_idx]
                        for r_idx in sp_hq.index
                        if r_idx not in used_row_ids
                        and (str(sp_hq.loc[r_idx].get("cat_str", "")), int(sp_hq.loc[r_idx].get("leaf_id_int", 1))) not in used_keys
                    ]
                if not win_unused:
                    selected_rows.append(None)
                    continue

                best_row = max(
                    win_unused,
                    key=lambda r: (
                        -abs(float(r.get("q_rank", q_target)) - q_target) if len(win_df) == 0 else 0.0,
                        int(bool(r.get("is_gold_silver", False))),
                        int(bool(r.get("is_cal_bool", False))),
                        float(r.get("rep_score", r.get("qual_weight", 0.0))),
                        -abs(float(r.get("q_rank", q_target)) - q_target),
                    ),
                )
                used_keys.add((str(best_row.get("cat_str", "")), int(best_row.get("leaf_id_int", 1))))
                used_row_ids.add(best_row.name)
                selected_rows.append(best_row)

            valid_sel = sorted(
                [r for r in selected_rows if r is not None],
                key=lambda r: float(r.get("ar_val", 0.0)),
            )
            selected_rows = valid_sel + [None] * max(0, 4 - len(valid_sel))

            for c_idx in range(4):
                ax_t = axes_b[row_idx, c_idx]
                row = selected_rows[c_idx] if c_idx < len(selected_rows) else None
                if row is None:
                    ax_t.set_facecolor("#F8F8F8")
                    ax_t.text(
                        0.5,
                        0.5,
                        f"n = {len(sp_hq)} leaves available\n(No additional distinct leaf)",
                        transform=ax_t.transAxes,
                        ha="center",
                        va="center",
                        fontsize=7.5,
                        color="#666666",
                        style="italic",
                    )
                    if row_idx == 0:
                        ax_t.set_title(pct_labels[c_idx], fontsize=7.6, fontweight="bold")
                    if c_idx == 0:
                        ax_t.set_ylabel(SPECIES_ITALIC[sp_name], fontsize=9.8, fontweight="bold", color=sp_col, labelpad=6)
                    ax_t.set_xticks([])
                    ax_t.set_yticks([])
                    for spine in ax_t.spines.values():
                        spine.set_edgecolor("#CCCCCC")
                        spine.set_linewidth(0.8)
                    continue

                cat_num = str(row.get("catalogNumber", f"V{row_idx * 4 + c_idx + 1}"))
                leaf_id = int(row.get("leaf_id", 1)) if pd.notna(row.get("leaf_id", 1)) else 1
                sol_v = float(row.get("sol_val", 0.88))
                ar_v = float(row.get("ar_val", 2.10))
                ppm_raw = row.get("pixels_per_mm", np.nan)
                ppm_val = float(ppm_raw) if pd.notna(ppm_raw) and str(ppm_raw).strip() != "" else np.nan

                mask_img: Optional[np.ndarray] = None
                m_path_res = _resolve_asset_path(row.get("mask_path", ""), out_dir)
                if m_path_res is not None and m_path_res.exists():
                    mask_img = cv2.imread(str(m_path_res), cv2.IMREAD_GRAYSCALE)
                    if mask_img is not None and mask_img.ndim == 3:
                        mask_img = mask_img[:, :, 0]
                    if mask_img is not None and np.count_nonzero(mask_img) > 0:
                        num_l, lbls, stats_cc, _ = cv2.connectedComponentsWithStats((mask_img > 0).astype(np.uint8), connectivity=8)
                        if num_l > 1:
                            best_l = 1 + int(np.argmax(stats_cc[1:, cv2.CC_STAT_AREA]))
                            mask_img = ((lbls == best_l).astype(np.uint8)) * 255

                if mask_img is not None and np.count_nonzero(mask_img) > 0:
                    from leafmachine3.core.bilateral import apply_affine, oriented_affine
                    from leafmachine3.core.imaging import crop_to_box, mask_bbox, rotate_image

                    pad = 14
                    m_u8 = (mask_img > 0).astype(np.uint8) * 255
                    # Seal any tiny internal herbivory pinholes (< 2.5% of blade area) for clean silhouette rendering
                    cnts_fill, hier_fill = cv2.findContours(m_u8, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
                    if cnts_fill and hier_fill is not None:
                        fg_cnt_px = float(np.count_nonzero(m_u8))
                        for idx_h, c_h in enumerate(cnts_fill):
                            if int(hier_fill[0][idx_h][3]) != -1 and float(cv2.contourArea(c_h)) < 0.025 * fg_cnt_px:
                                cv2.drawContours(m_u8, cnts_fill, idx_h, 255, thickness=cv2.FILLED)

                    midrib_raw = row.get("midrib_angle_deg", np.nan)
                    midrib_deg = (
                        float(midrib_raw)
                        if pd.notna(midrib_raw) and str(midrib_raw).strip() != ""
                        else -90.0
                    )
                    signed_cw = -(midrib_deg + 90.0)
                    angle_cw_m = 0.0 if abs(signed_cw) <= 0.25 else float(signed_cw % 360.0)

                    cxy_loaded: Optional[np.ndarray] = None
                    c_path_m = _resolve_asset_path(row.get("contour_path", ""), out_dir)
                    if c_path_m is not None and c_path_m.exists():
                        try:
                            cdf_m = pd.read_csv(c_path_m)
                            if {"x", "y"}.issubset(cdf_m.columns) and len(cdf_m) >= 3:
                                cxy_loaded = cdf_m[["x", "y"]].to_numpy(dtype=np.float64)
                        except Exception:
                            cxy_loaded = None

                    mat_m, _ = oriented_affine(m_u8.shape[:2], angle_cw_m)
                    crop_m = rotate_image(m_u8, angle_cw_m, bg=0, nearest=True)

                    rot_cxy: Optional[np.ndarray] = None
                    if cxy_loaded is not None and len(cxy_loaded) >= 3:
                        rot_rad = float(np.deg2rad(signed_cw))
                        cos_r, sin_r = float(np.cos(rot_rad)), float(np.sin(rot_rad))
                        c_mean = np.mean(cxy_loaded, axis=0)
                        rel_c = cxy_loaded - c_mean
                        unrot_xy = np.column_stack([
                            cos_r * rel_c[:, 0] + sin_r * rel_c[:, 1] + c_mean[0],
                            -sin_r * rel_c[:, 0] + cos_r * rel_c[:, 1] + c_mean[1],
                        ])
                        rot_cxy = apply_affine(mat_m, unrot_xy)

                    bbox_m = mask_bbox(crop_m)
                    if bbox_m is not None:
                        rx0 = max(0, int(bbox_m[0]) - pad)
                        ry0 = max(0, int(bbox_m[1]) - pad)
                        rx1 = min(crop_m.shape[1], int(bbox_m[2]) + pad)
                        ry1 = min(crop_m.shape[0], int(bbox_m[3]) + pad)
                        crop_m = crop_to_box(crop_m, (rx0, ry0, rx1, ry1), pad=0)
                        if rot_cxy is not None:
                            rot_cxy = rot_cxy - np.asarray([[float(rx0), float(ry0)]], dtype=np.float64)

                    ys_r, xs_r = np.where(crop_m > 0)
                    x_med_r = float(np.median(xs_r)) if len(xs_r) > 0 else float(crop_m.shape[1] * 0.5)
                    y_min_r = float(ys_r.min()) if len(ys_r) > 0 else 0.0
                    y_max_r = float(ys_r.max()) if len(ys_r) > 0 else float(crop_m.shape[0] - 1)
                    y_span_r = max(y_max_r - y_min_r, 1.0)
                    x_span_r = max(float(xs_r.max() - xs_r.min()) if len(xs_r) > 0 else 1.0, 1.0)

                    cnts_m, _ = cv2.findContours(crop_m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
                    if cnts_m:
                        cnt_pts = max(cnts_m, key=cv2.contourArea).reshape(-1, 2).astype(np.float64)
                    else:
                        cnt_pts = np.column_stack([xs_r, ys_r]).astype(np.float64)

                    # Use cxy[0] snapped to crop_m boundary when it sits in a recessed cordate sinus on a broad undivided blade;
                    # for dissected or elongate blades, anchor at the bottom of the upright rachis / petiole base on crop_m.
                    p_base = None
                    if (
                        rot_cxy is not None
                        and len(rot_cxy) >= 3
                        and ar_v < 1.80
                        and sol_v >= 0.80
                    ):
                        cand_bx, cand_by = float(rot_cxy[0, 0]), float(rot_cxy[0, 1])
                        if (
                            cand_by >= (y_max_r - 0.25 * y_span_r)
                            and abs(cand_bx - x_med_r) <= 0.16 * x_span_r
                        ):
                            bot_cnt = cnt_pts[
                                (cnt_pts[:, 1] >= y_max_r - 0.28 * y_span_r)
                                & (np.abs(cnt_pts[:, 0] - x_med_r) <= 0.18 * x_span_r)
                            ]
                            if len(bot_cnt) > 0:
                                d_sinus = np.hypot(bot_cnt[:, 0] - cand_bx, bot_cnt[:, 1] - cand_by)
                                s_idx = int(np.argmin(d_sinus))
                                p_base = (float(bot_cnt[s_idx, 0]), float(bot_cnt[s_idx, 1]))

                    if p_base is None:
                        bot_band = cnt_pts[cnt_pts[:, 1] >= (y_max_r - 0.06 * y_span_r)]
                        if len(bot_band) > 0:
                            x_bot_med = float(np.median(bot_band[:, 0]))
                            b_idx = int(
                                np.argmin(
                                    (y_max_r - bot_band[:, 1])
                                    + 0.75 * np.abs(bot_band[:, 0] - x_bot_med)
                                )
                            )
                            p_base = (float(bot_band[b_idx, 0]), float(bot_band[b_idx, 1]))
                        else:
                            p_base = (x_med_r, y_max_r)

                    upper_xs = xs_r[ys_r <= (y_min_r + 0.18 * y_span_r)]
                    x_top_ref = float(np.median(upper_xs)) if len(upper_xs) > 0 else x_med_r
                    top_band = cnt_pts[cnt_pts[:, 1] <= (y_min_r + 0.08 * y_span_r)]
                    if len(top_band) > 0:
                        a_idx = int(
                            np.argmin(
                                (top_band[:, 1] - y_min_r)
                                + 1.15 * np.abs(top_band[:, 0] - x_top_ref)
                            )
                        )
                        p_apex = (float(top_band[a_idx, 0]), float(top_band[a_idx, 1]))
                    else:
                        p_apex = (x_top_ref, y_min_r)

                    if render_mode == "ect_radial":
                        _ensure_ect_module_loaded()
                        from leafmachine3.core.ect_compute import compute_ect
                        from leafmachine3.reporting.ect_viz import render_radial_ect_overlay, visual_log

                        ect_res = compute_ect(crop_m, num_dirs=160, bound_radius=1.0, want_simple=True)
                        if ect_res is not None:
                            bgr_ect = render_radial_ect_overlay(
                                visual_log(ect_res.ect_matrix),
                                ect_res.thetas,
                                ect_res.thresholds,
                                ect_res.bound_radius,
                                ect_res.outline_norm,
                                cmap="viridis",
                            )
                            ax_t.imshow(bgr_ect[..., ::-1])
                        else:
                            ax_t.set_facecolor("#F8F8F8")
                    else:
                        h_m, w_m = crop_m.shape[:2]
                        canvas_rgb = np.full((h_m, w_m, 3), 248, dtype=np.uint8)
                        rgb_hex = sp_col.lstrip("#")
                        r_c, g_c, b_c = int(rgb_hex[0:2], 16), int(rgb_hex[2:4], 16), int(rgb_hex[4:6], 16)
                        canvas_rgb[crop_m > 0] = (r_c, g_c, b_c)

                        pad_bot = max(18, int(round(canvas_rgb.shape[0] * 0.16)))
                        canvas_rgb = np.pad(
                            canvas_rgb,
                            ((0, pad_bot), (0, 0), (0, 0)),
                            mode="constant",
                            constant_values=248,
                        )
                        h_c, w_c = canvas_rgb.shape[:2]

                        ax_t.imshow(canvas_rgb)
                        ax_t.plot(
                            [p_base[0], p_apex[0]],
                            [p_base[1], p_apex[1]],
                            color="#111111",
                            linestyle="--",
                            linewidth=1.3,
                        )
                        ax_t.scatter(
                            [p_base[0]],
                            [p_base[1]],
                            facecolors="none",
                            edgecolors="#E69F00",
                            linewidths=1.5,
                            s=48,
                            zorder=5,
                        )
                        # Draw 10 mm physical scale bar in dedicated bottom footer band
                        if np.isfinite(ppm_val) and ppm_val > 0:
                            bar_px = min(w_c * 0.45, 10.0 * ppm_val)
                            bx0 = w_c * 0.08
                            by0 = h_c - pad_bot * 0.45
                            ax_t.plot([bx0, bx0 + bar_px], [by0, by0], color="#111111", linewidth=2.2)
                            ax_t.text(
                                bx0 + bar_px / 2.0,
                                by0 - max(3.0, pad_bot * 0.18),
                                "10 mm",
                                fontsize=6.5,
                                fontweight="bold",
                                ha="center",
                                va="bottom",
                                color="#111111",
                            )
                else:
                    c_path_res = _resolve_asset_path(row.get("contour_path", ""), out_dir)
                    pts_xy = None
                    if c_path_res is not None and c_path_res.exists():
                        try:
                            cdf = pd.read_csv(c_path_res)
                            if {"x", "y"}.issubset(cdf.columns) and len(cdf) >= 3:
                                pts_xy = cdf[["x", "y"]].to_numpy(dtype=float)
                        except Exception:
                            pts_xy = None
                    if pts_xy is None:
                        # Never fabricate a silhouette: render an explicit placeholder tile instead
                        ax_t.set_facecolor("#F8F8F8")
                        ax_t.text(
                            0.5, 0.5, "Mask / contour asset\nunavailable on disk",
                            transform=ax_t.transAxes, ha="center", va="center",
                            fontsize=7.5, color="#666666", style="italic",
                        )
                    else:
                        p_base_pt = pts_xy[0].copy()
                        y_min_pt = float(pts_xy[:, 1].min())
                        y_span_pt = max(float(pts_xy[:, 1].max() - y_min_pt), 1.0)
                        top_band_pt = np.where(pts_xy[:, 1] <= y_min_pt + 0.12 * y_span_pt)[0]
                        if len(top_band_pt) > 0:
                            apex_idx_pt = int(
                                top_band_pt[
                                    int(
                                        np.argmin(
                                            pts_xy[top_band_pt, 1]
                                            + 0.45 * np.abs(pts_xy[top_band_pt, 0] - p_base_pt[0])
                                        )
                                    )
                                ]
                            )
                        else:
                            apex_idx_pt = int(np.argmin(pts_xy[:, 1]))
                        p_apex_pt = pts_xy[apex_idx_pt].copy()
                        if render_mode == "ect_radial":
                            _ensure_ect_module_loaded()
                            from leafmachine3.core.ect_compute import compute_ect
                            from leafmachine3.reporting.ect_viz import render_radial_ect_overlay, visual_log

                            min_xy = pts_xy.min(axis=0)
                            span_xy = max(float((pts_xy.max(axis=0) - min_xy).max()), 1e-6)
                            norm_pts = np.round(((pts_xy - min_xy) / span_xy) * 200.0 + 28.0).astype(np.int32)
                            poly_m = np.zeros((256, 256), dtype=np.uint8)
                            cv2.fillPoly(poly_m, [norm_pts], 255)
                            ect_res = compute_ect(poly_m, num_dirs=160, bound_radius=1.0, want_simple=True)
                            if ect_res is not None:
                                bgr_ect = render_radial_ect_overlay(
                                    visual_log(ect_res.ect_matrix),
                                    ect_res.thetas,
                                    ect_res.thresholds,
                                    ect_res.bound_radius,
                                    ect_res.outline_norm,
                                    cmap="viridis",
                                )
                                ax_t.imshow(bgr_ect[..., ::-1])
                            else:
                                ax_t.set_facecolor("#F8F8F8")
                        else:
                            ax_t.set_facecolor("#F8F8F8")
                            ax_t.fill(pts_xy[:, 0], pts_xy[:, 1], color=sp_col, alpha=0.82)
                            ax_t.plot(
                                [p_base_pt[0], p_apex_pt[0]],
                                [p_base_pt[1], p_apex_pt[1]],
                                color="#111111",
                                linestyle="--",
                                linewidth=1.3,
                            )
                            ax_t.scatter([p_base_pt[0]], [p_base_pt[1]], facecolors="none", edgecolors="#E69F00", linewidths=1.5, s=48, zorder=5)
                            ax_t.set_aspect("equal", adjustable="datalim")
                            ax_t.invert_yaxis()

                is_r_type = bool(row.get("is_type_bool", False))
                type_tag = " [TYPE]" if is_r_type else ""
                if row_idx == 0:
                    if c_idx == 0:
                        hdr_label = "Type Anchor / 20th %ile" if is_r_type else "20th %ile (Compact)"
                    else:
                        hdr_label = pct_labels[c_idx]
                    col_hdr = f"{hdr_label}\n"
                else:
                    col_hdr = ""
                ax_t.set_title(f"{col_hdr}{cat_num}{type_tag} (L{leaf_id}) | AR={ar_v:.2f}, Sol={sol_v:.2f}", fontsize=7.6, fontweight="bold")
                if c_idx == 0:
                    ax_t.set_ylabel(SPECIES_ITALIC[sp_name], fontsize=9.8, fontweight="bold", color=sp_col, labelpad=6)
                ax_t.set_xticks([])
                ax_t.set_yticks([])

        if render_mode == "ect_radial":
            fig_b.suptitle(
                f"Phase 02 - Taxon-Stratified {zone_label} Leaf Radial ECT Overlays (20th-80th %ile Aspect Ratio, Tip at North){title_suffix}",
                fontsize=11.5,
                fontweight="bold",
                y=0.99,
            )
        else:
            fig_b.suptitle(
                f"Phase 02 - Taxon-Stratified {zone_label} Leaf Silhouettes (20th-80th %ile Aspect Ratio), Midrib Axes & Petiole Anchors{title_suffix}",
                fontsize=11.5,
                fontweight="bold",
                y=0.99,
            )
        fig_b.tight_layout(rect=[0.01, 0, 1, 0.965])
        page_pdf, page_png = _save_fig_fn(fig_b, out_dir / out_stem, dpi=300)
        plt.close(fig_b)
        return page_pdf, page_png

    basal_sp_list = [sp for sp in present_sp if (t1_basal_pool["species_clean"] == sp).any()] or list(present_sp)
    fg_montage_sp = basal_sp_list[:8]
    fig02b_pdf, fig02b_png = _render_montage_page(
        fg_montage_sp,
        "fig02b_silhouette_montage",
        pool_df=t1_basal_pool,
        zone_label="Basal",
    )
    saved_paths["fig02b_pdf"] = fig02b_pdf
    saved_paths["fig02b_png"] = fig02b_png

    fig02d_pdf, fig02d_png = _render_montage_page(
        fg_montage_sp,
        "fig02d_basal_ect_radial_montage",
        pool_df=t1_basal_pool,
        zone_label="Basal",
        render_mode="ect_radial",
    )
    saved_paths["fig02d_pdf"] = fig02d_pdf
    saved_paths["fig02d_png"] = fig02d_png

    if len(basal_sp_list) > 8:
        supp_montage_sp = basal_sp_list[8:]
        supp_pdf, supp_png = _render_montage_page(
            supp_montage_sp,
            "fig02b_silhouette_montage_supplement",
            title_suffix=" (Supplementary Taxa)",
            pool_df=t1_basal_pool,
            zone_label="Basal",
        )
        saved_paths["fig02b_supplement_pdf"] = supp_pdf
        saved_paths["fig02b_supplement_png"] = supp_pdf.with_suffix(".png") if supp_png is None else supp_png

    if len(t1_cauline_pool) > 0:
        cauline_sp_list = [sp for sp in present_sp if (t1_cauline_pool["species_clean"] == sp).any()]
        if cauline_sp_list:
            fg_caul_sp = cauline_sp_list[:8]
            fig02c_pdf, fig02c_png = _render_montage_page(
                fg_caul_sp,
                "fig02c_cauline_silhouette_montage",
                pool_df=t1_cauline_pool,
                zone_label="Cauline",
            )
            saved_paths["fig02c_pdf"] = fig02c_pdf
            saved_paths["fig02c_png"] = fig02c_png
            fig02e_pdf, fig02e_png = _render_montage_page(
                fg_caul_sp,
                "fig02e_cauline_ect_radial_montage",
                pool_df=t1_cauline_pool,
                zone_label="Cauline",
                render_mode="ect_radial",
            )
            saved_paths["fig02e_pdf"] = fig02e_pdf
            saved_paths["fig02e_png"] = fig02e_png
            if len(cauline_sp_list) > 8:
                supp_caul_sp = cauline_sp_list[8:]
                supp_c_pdf, supp_c_png = _render_montage_page(
                    supp_caul_sp,
                    "fig02c_cauline_silhouette_montage_supplement",
                    title_suffix=" (Supplementary Taxa)",
                    pool_df=t1_cauline_pool,
                    zone_label="Cauline",
                )
                saved_paths["fig02c_supplement_pdf"] = supp_c_pdf
                saved_paths["fig02c_supplement_png"] = supp_c_png

    # Export authentic LeafMachine3 radial ECT overlays for downstream Step 03 publication figures
    try:
        export_taxon_ect_exemplars(
            extracted_df=df_plot,
            figures_dir=out_dir,
            zone="basal",
            n_dirs=160,
        )
        if len(cauline_sp_list) > 0:
            export_taxon_ect_exemplars(
                extracted_df=df_plot,
                figures_dir=out_dir,
                zone="cauline",
                n_dirs=160,
            )
    except Exception as e_ect:
        logger.warning(f"Could not export taxon ECT exemplars: {e_ect}")

    logger.info(f"Step 02 publication figures exported to: {out_dir}")
    return saved_paths


def regenerate_step02_figures_from_tables(
    output_dir: Union[str, Path],
    vouchers_csv: Optional[Union[str, Path]] = None,
) -> Dict[str, Path]:
    """
    Reload Step 02 manifest tables (`extracted_leaf_manifest.csv`, `failed_qc_vouchers.csv`,
    `voucher_reproductive_metrics.csv`, `voucher_phenological_states.csv`) from `output_dir / 'tables'`
    alongside downstream Step 03 (`leaf_efa_harmonics_individual.csv`, `leaf_efa_harmonics.csv`) and
    Step 04 (`morphometric_flags.csv`) tables, and re-render `fig02a_segmentation_qc_distributions.{pdf,png}`
    and `fig02b_silhouette_montage.{pdf,png}` into `output_dir.parent / 'outputs' / 'figures'`.
    """
    out_data_dir = Path(output_dir)
    tables_dir = out_data_dir / "tables" if (out_data_dir / "tables").exists() else out_data_dir
    figures_dir = out_data_dir.parent / "outputs" / "figures"

    def _safe_read_csv(p: Path) -> pd.DataFrame:
        if p.exists() and p.stat().st_size > 0:
            try:
                return pd.read_csv(p)
            except Exception as e_csv:
                logger.warning(f"Failed to read {p}: {e_csv}")
        return pd.DataFrame()

    extracted_df = _safe_read_csv(tables_dir / "extracted_leaf_manifest.csv")
    failed_df = _safe_read_csv(tables_dir / "failed_qc_vouchers.csv")
    repro_df = _safe_read_csv(tables_dir / "voucher_reproductive_metrics.csv")
    pheno_df = _safe_read_csv(tables_dir / "voucher_phenological_states.csv")

    if vouchers_csv is not None and Path(vouchers_csv).exists() and len(extracted_df) > 0:
        try:
            vdf = pd.read_csv(vouchers_csv, dtype=str)
            if "catalogNumber" in vdf.columns and "catalogNumber" in extracted_df.columns:
                active_vcats = set(vdf["catalogNumber"].astype(str).str.strip())
                extracted_df = extracted_df[
                    extracted_df["catalogNumber"].astype(str).str.strip().isin(active_vcats)
                ].copy()
        except Exception:
            pass

    return generate_step02_publication_figures(
        extracted_df=extracted_df,
        failed_df=failed_df,
        repro_df=repro_df,
        pheno_df=pheno_df,
        figures_dir=figures_dir,
    )


def export_taxon_ect_exemplars(
    extracted_df: Optional[pd.DataFrame] = None,
    figures_dir: Union[str, Path] = "outputs/figures",
    zone: str = "basal",
    n_dirs: int = 160,
    max_taxa: int = 9,
    vouchers_csv: Optional[Union[str, Path]] = None,
    tables_dir: Optional[Union[str, Path]] = None,
    db_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Exports authentic LeafMachine3 radial ECT overlays (Chitwood polar pcolormesh + leaf contour)
    for taxon medoid/exemplar leaves directly into `<figures_dir>/ect_overlays/<zone>_<taxon_slug>.png`,
    accompanied by `<figures_dir>/ect_overlays/manifest_<zone>.json`.
    """
    fig_dir = Path(figures_dir).resolve()
    overlays_dir = fig_dir / "ect_overlays"
    overlays_dir.mkdir(parents=True, exist_ok=True)
    zone_clean = str(zone).lower().strip()

    # Load extracted_df if missing
    if extracted_df is None or len(extracted_df) == 0:
        cand_tables = []
        if tables_dir is not None:
            cand_tables.append(Path(tables_dir) / "extracted_leaf_manifest.csv")
            cand_tables.append(Path(tables_dir) / "leaf_efa_harmonics_individual.csv")
            cand_tables.append(Path(tables_dir) / "leaf_efa_harmonics.csv")
        for base in (fig_dir.parent.parent / "data" / "tables", fig_dir.parent / "data" / "tables", PROJECT_ROOT / "data" / "tables"):
            cand_tables.append(base / "extracted_leaf_manifest.csv")
            cand_tables.append(base / "leaf_efa_harmonics_individual.csv")
            cand_tables.append(base / "leaf_efa_harmonics.csv")
        for ct in cand_tables:
            if ct.exists() and ct.stat().st_size > 0:
                try:
                    extracted_df = pd.read_csv(ct, low_memory=False)
                    break
                except Exception:
                    pass

    if extracted_df is None or len(extracted_df) == 0:
        logger.warning(f"export_taxon_ect_exemplars: No extracted leaf data found to export for zone '{zone_clean}'.")
        return {}

    df_work = extracted_df.copy()

    # Map species_standardized from vouchers if missing or uninformative
    has_sp = "species_standardized" in df_work.columns or "species_clean" in df_work.columns
    if not has_sp and vouchers_csv is not None and Path(vouchers_csv).exists():
        try:
            v_df = pd.read_csv(vouchers_csv, dtype=str)
            if "catalogNumber" in v_df.columns and "species_standardized" in v_df.columns:
                v_map = dict(zip(v_df["catalogNumber"].astype(str).str.strip(), v_df["species_standardized"]))
                df_work["species_standardized"] = df_work["catalogNumber"].astype(str).str.strip().map(v_map)
        except Exception:
            pass

    if "species_clean" not in df_work.columns:
        sp_src = df_work.get("species_standardized", df_work.get("scientificName", df_work.get("species_raw", df_work.get("taxon", ""))))
        df_work["species_clean"] = sp_src.astype(str).apply(_standardize_packera_species_py)

    # Normalize assigned_tier
    t_col = df_work.get("assigned_tier", df_work.get("tier_assignment", df_work.get("tier", "Tier_1_Direct")))
    is_t1 = t_col.astype(str).str.strip().str.lower().isin({
        "tier_1_direct", "tier1", "tier 1", "tier 1 (dissected)", "tier 1 (bipinnatisect)"
    })
    df_t1 = df_work[is_t1].copy()
    if len(df_t1) == 0:
        df_t1 = df_work.copy()

    # Filter by leaf_zone
    if "leaf_zone" in df_t1.columns:
        lz_norm = df_t1["leaf_zone"].astype(str).str.strip().str.lower()
        if zone_clean == "cauline":
            df_zone = df_t1[lz_norm == "cauline"].copy()
        else:
            df_zone = df_t1[lz_norm.isin({"basal", "", "nan", "none"})].copy()
            if len(df_zone) == 0:
                df_zone = df_t1.copy()
    else:
        df_zone = df_t1.copy()

    # Rank taxa by abundance
    sp_counts = df_zone["species_clean"].value_counts()
    valid_sp = [s for s in sp_counts.index if s not in {"Unassigned_Packera", "Unknown", ""} and pd.notna(s)]
    target_sp = valid_sp[:max_taxa]

    _ensure_ect_module_loaded()
    from leafmachine3.core.ect_compute import compute_ect
    from leafmachine3.reporting.ect_viz import render_radial_ect_overlay, visual_log

    manifest: Dict[str, Dict[str, Any]] = {}
    for sp_name in target_sp:
        sp_df = df_zone[df_zone["species_clean"] == sp_name]
        if len(sp_df) == 0:
            continue

        # Exemplar selection (medoid)
        has_ect_pcs = "ECT_PC1" in sp_df.columns and "ECT_PC2" in sp_df.columns
        if has_ect_pcs:
            pc1 = pd.to_numeric(sp_df["ECT_PC1"], errors="coerce")
            pc2 = pd.to_numeric(sp_df["ECT_PC2"], errors="coerce")
            valid_pcs = pc1.notna() & pc2.notna()
        else:
            valid_pcs = pd.Series(False, index=sp_df.index)

        best_row = None
        if valid_pcs.sum() >= 2:
            sub_pc = sp_df[valid_pcs]
            m_pc1 = sub_pc["ECT_PC1"].astype(float).mean()
            m_pc2 = sub_pc["ECT_PC2"].astype(float).mean()
            dist_pc = np.sqrt((sub_pc["ECT_PC1"].astype(float) - m_pc1)**2 + (sub_pc["ECT_PC2"].astype(float) - m_pc2)**2)
            # Prioritize type specimen if within 1.5 * median distance
            is_type_s = sub_pc["is_type_specimen"].apply(_parse_bool_val) if "is_type_specimen" in sub_pc.columns else pd.Series(False, index=sub_pc.index)
            type_hits = sub_pc[is_type_s & (dist_pc <= 1.5 * dist_pc.median())]
            if len(type_hits) > 0:
                best_row = type_hits.iloc[0]
            else:
                best_row = sub_pc.loc[dist_pc.idxmin()]
        else:
            sol_s = pd.to_numeric(sp_df.get("solidity", sp_df.get("sol_val", 0.88)), errors="coerce").fillna(0.88)
            ar_s = pd.to_numeric(sp_df.get("aspect_ratio", sp_df.get("ar_val", 2.10)), errors="coerce").fillna(2.10)
            d_geom = (sol_s - sol_s.median())**2 + (ar_s - ar_s.median())**2
            best_row = sp_df.loc[d_geom.idxmin()]

        if best_row is None:
            continue

        cat_num = str(best_row.get("catalogNumber", "Unknown"))
        leaf_id = int(best_row.get("leaf_id", 1)) if pd.notna(best_row.get("leaf_id", 1)) else 1
        taxon_slug = re.sub(r"[^A-Za-z0-9_]+", "_", sp_name.lower())
        out_png = overlays_dir / f"{zone_clean}_{taxon_slug}.png"

        bgr_rendered = None
        # Try loading existing pre-rendered radial overlay if valid
        if out_png.exists() and out_png.stat().st_size > 1000:
            bgr_rendered = cv2.imread(str(out_png))

        if bgr_rendered is None:
            mask_img = None
            m_path_res = _resolve_asset_path(best_row.get("mask_path", ""), fig_dir)
            if m_path_res is not None and m_path_res.exists():
                mask_img = cv2.imread(str(m_path_res), cv2.IMREAD_GRAYSCALE)
                if mask_img is not None and mask_img.ndim == 3:
                    mask_img = mask_img[:, :, 0]
            if mask_img is None or np.count_nonzero(mask_img) == 0:
                # Fallback: check contour coordinates
                c_path_res = _resolve_asset_path(best_row.get("contour_path", ""), fig_dir)
                if c_path_res is not None and c_path_res.exists():
                    try:
                        c_df = pd.read_csv(c_path_res)
                        pts = c_df[["x", "y"]].to_numpy(dtype=np.float32)
                        if len(pts) >= 6:
                            # Normalize into 400x400 canvas
                            min_xy = pts.min(axis=0)
                            span_xy = np.maximum(pts.max(axis=0) - min_xy, 1.0)
                            scale_f = 340.0 / float(np.max(span_xy))
                            norm_pts = ((pts - min_xy) * scale_f + 30.0).astype(np.int32)
                            poly_m = np.zeros((400, 400), dtype=np.uint8)
                            cv2.fillPoly(poly_m, [norm_pts], 255)
                            mask_img = poly_m
                    except Exception:
                        pass

            if mask_img is not None and np.count_nonzero(mask_img) > 0:
                m_u8 = (mask_img > 0).astype(np.uint8) * 255
                ect_res = compute_ect(m_u8, num_dirs=n_dirs, bound_radius=1.0, want_simple=True)
                if ect_res is not None:
                    bgr_rendered = render_radial_ect_overlay(
                        visual_log(ect_res.ect_matrix),
                        ect_res.thetas,
                        ect_res.thresholds,
                        ect_res.bound_radius,
                        ect_res.outline_norm,
                        cmap="viridis",
                    )
                    cv2.imwrite(str(out_png), bgr_rendered)

        comp_val = best_row.get("ect_lobe_complexity", best_row.get("cauline_ect_lobe_complexity", np.nan))
        comp_float = float(comp_val) if pd.notna(comp_val) else np.nan
        manifest[sp_name] = {
            "image_path": str(out_png),
            "catalogNumber": cat_num,
            "leaf_id": leaf_id,
            "taxon": sp_name,
            "taxon_slug": taxon_slug,
            "is_type": bool(_parse_bool_val(best_row.get("is_type_specimen", False))),
            "ect_lobe_complexity": round(comp_float, 4) if np.isfinite(comp_float) else None,
            "aspect_ratio": round(float(best_row.get("aspect_ratio", np.nan)), 4) if pd.notna(best_row.get("aspect_ratio")) else None,
            "solidity": round(float(best_row.get("solidity", np.nan)), 4) if pd.notna(best_row.get("solidity")) else None,
        }

    manifest_path = overlays_dir / f"manifest_{zone_clean}.json"
    with open(manifest_path, "w", encoding="utf-8") as f_m:
        json.dump(manifest, f_m, indent=2)

    logger.info(f"Exported {len(manifest)} authentic LeafMachine3 radial ECT overlays to {overlays_dir}")
    return manifest


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 02 Figures and LeafMachine3 ECT Exemplar Exporter")
    parser.add_argument("--export-ect-exemplars", action="store_true", help="Export authentic LeafMachine3 radial ECT overlays for taxon exemplars")
    parser.add_argument("--output-dir", type=str, default="outputs/figures", help="Output directory for figures/overlays")
    parser.add_argument("--zone", type=str, default="all", choices=["basal", "cauline", "all"], help="Leaf zone (basal, cauline, or all)")
    parser.add_argument("--vouchers-csv", type=str, default="", help="Path to curated_vouchers.csv")
    parser.add_argument("--manifest-csv", type=str, default="", help="Path to extracted_leaf_manifest.csv")
    parser.add_argument("--regenerate-figures-only", action="store_true", help="Regenerate Step 02 publication figures")
    args = parser.parse_args()

    if args.export_ect_exemplars:
        out_fdir = Path(args.output_dir)
        v_csv = args.vouchers_csv or None
        zones = ["basal", "cauline"] if args.zone == "all" else [args.zone]
        for z in zones:
            export_taxon_ect_exemplars(figures_dir=out_fdir, zone=z, vouchers_csv=v_csv)
    elif args.regenerate_figures_only:
        regenerate_step02_figures_from_tables(output_dir=args.output_dir, vouchers_csv=args.vouchers_csv or None)

