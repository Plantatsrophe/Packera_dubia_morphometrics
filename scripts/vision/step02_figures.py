#!/usr/bin/env python3
"""
===============================================================================
Module: step02_figures.py
Project: Packera Species Delimitation & Integrative Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Publication-grade 300-DPI PNG figure suite and native LeafMachine3 visual
    asset exporter for Step 02:
      1. fig02a_lm3_organ_and_ect_exemplars.png:
         Dynamic taxon-stratified organ montage delegating directly to native
         LeafMachine3 `render_leaf_products` (cropped RGB lamina, hole-aware
         lamina+petiole silhouettes, pose landmarks, 10 mm calibrated scale bar)
         and `render_radial_ect_overlay` polar ECT badges across Basal Rosette
         and Cauline Stem tracks (scales to any number of species/ranks).
      2. fig02b_segmentation_qc_and_phenophase.png:
         Multi-panel segmentation QC, ruler calibration, basal-vs-cauline
         solidity x aspect-ratio homology separation (unfilled hollow markers,
         reniform L/W < 1.0 support), involucral AR_inv distributions, and
         determinate cyme phenophase proportions.
      3. lm3_assets/ & ect_overlays/:
         Cached native LeafMachine3 exemplar silhouettes and polar ECT PNGs for
         reuse in Step 03, Step 06, and the offline interactive HTML dashboard.
      4. lm3_sheet_overlays/:
         Capped whole-sheet LeafMachine3 diagnostic overlays (Type specimens,
         DBSCAN multi-plant ramet sheets, and heteroblastic/flagged sheets).
===============================================================================
"""

from __future__ import annotations

import json
import logging
import math
import os
import re
import sys
import textwrap
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union

import cv2
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import scripts.vision.lm3_compat  # noqa: F401
from scripts.vision.lm3_inference_engines import _ensure_ect_module_loaded

logger = logging.getLogger("SegmentAndExtract")

__all__ = [
    "_atomic_save_matplotlib_fig",
    "_standardize_packera_species_py",
    "format_icn_taxon_mathtext",
    "get_combinatorial_taxon_style",
    "_parse_bool_val",
    "_ASSET_ANCHOR_DIRS",
    "_resolve_asset_path",
    "_is_valid_smooth_closed_mask",
    "_discover_downstream_table",
    "generate_step02_publication_figures",
    "regenerate_step02_figures_from_tables",
    "export_taxon_ect_exemplars",
    "export_taxon_lm3_exemplars",
    "export_lm3_sheet_diagnostic_overlays",
]

# 12-hue Tol + Okabe-Ito colorblind-safe categorical palette
TOL_OKABE_12: Tuple[str, ...] = (
    "#0072B2",  # Blue
    "#D55E00",  # Vermilion
    "#009E73",  # Bluish green
    "#CC79A7",  # Reddish purple
    "#E69F00",  # Orange
    "#56B4E9",  # Sky blue
    "#332288",  # Indigo
    "#882255",  # Wine
    "#44AA99",  # Teal
    "#999933",  # Olive
    "#888888",  # Slate
    "#117733",  # Forest green
)

# 5 unfilled/hollow marker geometries -> 12 x 5 = 60 unique taxon combinations
HOLLOW_MARKERS_5: Tuple[str, ...] = ("o", "^", "D", "s", "v")

ICN_ROMAN_TOKENS: Set[str] = {
    "subsp.",
    "subsp",
    "ssp.",
    "ssp",
    "var.",
    "var",
    "f.",
    "forma",
    "×",
    "x",
    "aff.",
    "aff",
    "cf.",
    "cf",
    "sp.",
    "sp",
    "s.s.",
    "s.l.",
    "ined.",
}


def _atomic_save_matplotlib_fig(fig: Any, out_base_path: Path, dpi: int = 300) -> Path:
    """Atomically save a Matplotlib figure strictly as a 300-DPI PNG and remove stale PDF siblings."""
    out_base = Path(out_base_path)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    png_path = out_base.with_suffix(".png")
    stale_pdf = out_base.with_suffix(".pdf")
    tmp_png = out_base.with_name(f".{out_base.stem}_{os.getpid()}_tmp.png")
    try:
        fig.savefig(str(tmp_png), format="png", dpi=int(dpi), bbox_inches="tight", facecolor="white")
        tmp_png.replace(png_path)
        if stale_pdf.exists():
            stale_pdf.unlink(missing_ok=True)
    finally:
        if tmp_png.exists():
            tmp_png.unlink(missing_ok=True)
    return png_path


def _standardize_packera_species_py(raw_val: Any, infraspecific_val: Any = None) -> str:
    """Standardize Packera species, nothospecies, infraspecific ranks (subsp./var./f.), or morphotype IDs."""
    if raw_val is None or pd.isna(raw_val):
        return "Unassigned_Packera"
    s_raw = str(raw_val).strip()
    if not s_raw or s_raw.lower() in {"none", "nan", "null", "na", "unknown", "unassigned_packera"}:
        return "Unassigned_Packera"

    if re.match(r"^(?:Morphotype|Ecotype|Complex)_[A-Za-z0-9_.-]+$", s_raw):
        return s_raw

    # Extract explicit infraspecific rank from raw string or companion column before binomial standardization
    infra_suffix = ""
    if infraspecific_val is not None and not pd.isna(infraspecific_val):
        inf_str = str(infraspecific_val).strip()
        if inf_str and inf_str.lower() not in {"none", "nan", "null", "na"}:
            m_inf = re.search(r"\b(subsp\.?|ssp\.?|var\.?|f\.?|forma)\s+([a-z][a-z-]+)\b", inf_str, flags=re.I)
            if m_inf:
                rk = m_inf.group(1).lower().rstrip(".")
                rk_std = "subsp." if rk in {"subsp", "ssp"} else ("var." if rk == "var" else "f.")
                infra_suffix = f" {rk_std} {m_inf.group(2).lower()}"
            elif re.match(r"^[a-z][a-z-]+$", inf_str.lower()):
                infra_suffix = f" var. {inf_str.lower()}"

    if not infra_suffix:
        m_raw_inf = re.search(r"\b(subsp\.?|ssp\.?|var\.?|f\.?|forma)\s+([a-z][a-z-]+)\b", s_raw, flags=re.I)
        if m_raw_inf:
            rk = m_raw_inf.group(1).lower().rstrip(".")
            rk_std = "subsp." if rk in {"subsp", "ssp"} else ("var." if rk == "var" else "f.")
            infra_suffix = f" {rk_std} {m_raw_inf.group(2).lower()}"

    base_std = None
    try:
        from scripts.core.harvester import standardize_packera_taxon

        std = standardize_packera_taxon(s_raw)
        if std and str(std).strip() and str(std).strip().lower() not in {
            "none",
            "nan",
            "null",
            "na",
            "unknown",
            "unassigned_packera",
        }:
            base_std = str(std).strip()
    except Exception:
        base_std = None

    if base_std is None:
        cleaned = s_raw.split("(")[0].strip()
        m_binom = re.match(r"^([A-Za-z][a-z]+)\s+(?:([x×])\s+)?([a-z][a-z-]+)\b", cleaned)
        if m_binom:
            genus = m_binom.group(1).capitalize()
            hyb = "× " if m_binom.group(2) else ""
            base_std = f"{genus} {hyb}{m_binom.group(3).lower()}"

    has_hybrid_mark = bool(re.search(r"(?:^|\s)[×x]\s+[a-z]", s_raw))
    if base_std and has_hybrid_mark and "×" not in base_std:
        b_parts = base_std.split()
        if len(b_parts) >= 2:
            base_std = f"{b_parts[0]} × {' '.join(b_parts[1:])}"

    if base_std:
        if infra_suffix and infra_suffix.split()[-1] not in base_std.split()[1:2]:
            return f"{base_std}{infra_suffix}"
        # Retain autonym if explicitly present in raw input
        if infra_suffix and re.search(r"\b(subsp\.?|ssp\.?|var\.?|f\.?)\s+", s_raw, flags=re.I):
            return f"{base_std}{infra_suffix}"
        return base_std

    return "Unassigned_Packera"


def format_icn_taxon_mathtext(label_raw: Any, abbreviate_genus: bool = False) -> str:
    """
    Format a botanical name using ICN typographic conventions in Matplotlib mathtext:
      - Genus, specific epithet, and infraspecific epithets rendered in italics.
      - Rank connectors (`subsp.`, `var.`, `f.`), hybrid symbols (`×`), open-nomenclature
        qualifiers (`aff.`, `cf.`, `s.s.`, `s.l.`), and morphotype/cluster IDs rendered upright.
    """
    s = str(label_raw or "").strip()
    if not s:
        return "Unknown"
    if re.match(r"^(?:Morphotype|Ecotype|Complex|Unassigned|Cluster|Unknown)\b", s, flags=re.I):
        return s.replace("_", " ")

    tokens = s.split()
    if len(tokens) < 2 or not re.match(r"^[A-Z][a-z]+$", tokens[0]):
        return s.replace("_", " ")

    out_parts: List[str] = []
    for idx, tok in enumerate(tokens):
        low = tok.lower()
        if idx == 0:
            g_txt = f"{tok[0]}." if abbreviate_genus else tok
            out_parts.append(rf"$\it{{{g_txt}}}$")
        elif low in ICN_ROMAN_TOKENS:
            if low in {"x", "×"}:
                out_parts.append("×")
            elif low in {"subsp", "ssp", "subsp.", "ssp."}:
                out_parts.append("subsp.")
            elif low in {"var", "var."}:
                out_parts.append("var.")
            elif low in {"f", "f.", "forma"}:
                out_parts.append("f.")
            else:
                out_parts.append(tok)
        elif re.match(r"^[a-z][a-z-]+$", tok):
            safe_tok = tok.replace("-", r"\text{-}") if "-" in tok else tok
            # Matplotlib mathtext handles plain hyphens cleanly inside \mathit
            out_parts.append(rf"$\it{{{tok}}}$")
        else:
            out_parts.append(tok)
    return " ".join(out_parts)


def get_combinatorial_taxon_style(taxa_seq: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """Map an arbitrary number of taxa (1..60+) to distinct (color, hollow_marker) pairs."""
    unique_taxa = list(dict.fromkeys(str(t) for t in taxa_seq if str(t).strip()))
    styles: Dict[str, Tuple[str, str]] = {}
    n_c = len(TOL_OKABE_12)
    n_m = len(HOLLOW_MARKERS_5)
    for idx, tx in enumerate(unique_taxa):
        col = TOL_OKABE_12[idx % n_c]
        mk = HOLLOW_MARKERS_5[(idx // n_c) % n_m] if len(unique_taxa) > 8 else HOLLOW_MARKERS_5[idx % n_m]
        styles[tx] = (col, mk)
    return styles


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


_ASSET_ANCHOR_DIRS = ("masks", "contours", "cropped_patches", "raw_vouchers")


def _resolve_asset_path(p_str: Any, out_dir: Path) -> Optional[Path]:
    """Resolve mask, contour, crop, or raw voucher scan path across run-scoped and canonical stores."""
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
    """Verify binary mask has unclipped borders, coherent lamina islands, and clean geometry."""
    resolved = _resolve_asset_path(m_path_str, out_dir)
    raw_str = str(m_path_str or "").strip()
    if resolved is None:
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
    max_edge_touch = max(3, int(0.02 * max(h_m, w_m)))
    if (
        int(np.count_nonzero(fg[0, :])) > max_edge_touch
        or int(np.count_nonzero(fg[-1, :])) > max_edge_touch
        or int(np.count_nonzero(fg[:, 0])) > max_edge_touch
        or int(np.count_nonzero(fg[:, -1])) > max_edge_touch
    ):
        return False

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
            # Allow multi-island deeply pinnatifid Packera leaves when clean and unclipped
            if sig_comps > 12:
                return False
            best_cc = 1 + int(np.argmax(cc_areas))
            fg = (lbls_cc == best_cc).astype(np.uint8)
            total_fg = int(np.count_nonzero(fg))

    cnts, hier = cv2.findContours(fg, cv2.RETR_CCOMP, cv2.CHAIN_APPROX_NONE)
    if not cnts or hier is None:
        return False
    outer_cnts = [
        cnts[idx_c]
        for idx_c in range(len(cnts))
        if int(hier[0][idx_c][3]) == -1 and float(cv2.contourArea(cnts[idx_c])) >= min(25, total_fg // 2)
    ]
    if not outer_cnts:
        return False
    main_cnt = max(outer_cnts, key=cv2.contourArea)
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
    circ = (4.0 * math.pi * area_c) / (perim_c * perim_c)
    if circ < 0.12 and not (circ >= 0.020 and sol_c < 0.78):
        return False
    return True


def _discover_downstream_table(filename: str, out_dir: Path, active_cats: Set[str]) -> Optional[Path]:
    """Locate downstream QC/morphometric CSV tables in run-scoped or workspace-root table directories."""
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


def _load_mask_or_contour_u8(row: pd.Series, out_dir: Path) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Load upright binary mask (`uint8` 0/255) and optional cropped RGB lamina patch (`uint8` HxWx3)
    using native LeafMachine3 `render_leaf_products` orientation rules.
    """
    mask_u8: Optional[np.ndarray] = None
    rgb_patch: Optional[np.ndarray] = None

    m_path = _resolve_asset_path(row.get("mask_path", ""), out_dir)
    if m_path is not None and m_path.exists():
        m_in = cv2.imread(str(m_path), cv2.IMREAD_GRAYSCALE)
        if m_in is not None and np.count_nonzero(m_in) > 0:
            mask_u8 = ((m_in > 0).astype(np.uint8)) * 255

    if mask_u8 is None:
        c_path = _resolve_asset_path(row.get("contour_path", ""), out_dir)
        if c_path is not None and c_path.exists():
            try:
                c_df = pd.read_csv(c_path)
                if {"x", "y"}.issubset(c_df.columns) and len(c_df) >= 6:
                    pts = c_df[["x", "y"]].to_numpy(dtype=np.float32)
                    min_xy = pts.min(axis=0)
                    span_xy = np.maximum(pts.max(axis=0) - min_xy, 1.0)
                    scale_f = 320.0 / float(np.max(span_xy))
                    norm_pts = ((pts - min_xy) * scale_f + 32.0).astype(np.int32)
                    canvas = np.zeros((384, 384), dtype=np.uint8)
                    cv2.fillPoly(canvas, [norm_pts], 255)
                    mask_u8 = canvas
            except Exception:
                pass

    if mask_u8 is None or int(np.count_nonzero(mask_u8)) == 0:
        return None, None

    is_upright = (
        _parse_bool_val(row.get("mask_is_upright", True))
        if ("mask_is_upright" in row and pd.notna(row.get("mask_is_upright")))
        else True
    )
    angle_cw: Optional[float] = None
    if not is_upright:
        m_ang = row.get("midrib_angle_deg", np.nan)
        if pd.notna(m_ang) and str(m_ang).strip() != "":
            scw = -(float(m_ang) + 90.0)
            if abs(scw % 360.0) > 0.25:
                angle_cw = float(scw % 360.0)

    # Attempt to extract authentic RGB crop from voucher sheet if bbox and scan exist
    img_cand = _resolve_asset_path(row.get("image_path", ""), out_dir)
    if img_cand is None:
        cat_s = str(row.get("catalogNumber", "")).strip()
        if cat_s:
            for ext in (".webp", ".jpg", ".jpeg", ".png", ".tif"):
                cand_raw = _resolve_asset_path(f"raw_vouchers/{cat_s}{ext}", out_dir)
                if cand_raw is not None:
                    img_cand = cand_raw
                    break

    if img_cand is not None and img_cand.exists():
        b_cols = ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")
        if all(c in row and pd.notna(row.get(c)) for c in b_cols):
            try:
                bgr_full = cv2.imread(str(img_cand), cv2.IMREAD_COLOR)
                if bgr_full is not None:
                    H_f, W_f = bgr_full.shape[:2]
                    x1 = max(0, min(W_f - 1, int(round(float(row["bbox_x1"])))))
                    y1 = max(0, min(H_f - 1, int(round(float(row["bbox_y1"])))))
                    x2 = max(x1 + 1, min(W_f, int(round(float(row["bbox_x2"])))))
                    y2 = max(y1 + 1, min(H_f, int(round(float(row["bbox_y2"])))))
                    crop_rgb = cv2.cvtColor(bgr_full[y1:y2, x1:x2], cv2.COLOR_BGR2RGB)
                    if crop_rgb.shape[:2] != mask_u8.shape[:2]:
                        crop_rgb = cv2.resize(crop_rgb, (mask_u8.shape[1], mask_u8.shape[0]), interpolation=cv2.INTER_AREA)
                    rgb_patch = crop_rgb
            except Exception:
                rgb_patch = None

    # Delegate orientation & cutout generation directly to LeafMachine3 render_leaf_products
    from leafmachine3.reporting.leaf_products import render_leaf_products

    base_rgb = rgb_patch if rgb_patch is not None else np.full((*mask_u8.shape[:2], 3), 42, dtype=np.uint8)
    prods = render_leaf_products(
        base_rgb,
        {"silhouette": mask_u8 > 0, "silhouettepet": mask_u8 > 0, "hole": None},
        angle_cw=angle_cw,
        want={"silhouette", "rgb_cropped"} if rgb_patch is not None else {"silhouette"},
    )
    sil_prod = prods.get("silhouette")
    if sil_prod is not None and np.count_nonzero(sil_prod) > 0:
        mask_u8 = (np.asarray(sil_prod) > 0).astype(np.uint8) * 255
    if rgb_patch is not None:
        rgb_prod = prods.get("rgb_cropped")
        if isinstance(rgb_prod, np.ndarray) and rgb_prod.ndim == 3 and rgb_prod.shape[:2] == mask_u8.shape[:2]:
            rgb_patch = rgb_prod
        elif angle_cw is not None:
            from leafmachine3.core.imaging import rotate_image

            rgb_patch = rotate_image(rgb_patch, angle_cw, bg=255, nearest=False)

    # Tightly crop around foreground with clean margin
    ys, xs = np.where(mask_u8 > 0)
    if len(xs) > 0 and len(ys) > 0:
        pad = max(10, int(0.08 * max(ys.max() - ys.min() + 1, xs.max() - xs.min() + 1)))
        y0 = max(0, int(ys.min()) - pad)
        y1 = min(mask_u8.shape[0], int(ys.max()) + pad + 1)
        x0 = max(0, int(xs.min()) - pad)
        x1 = min(mask_u8.shape[1], int(xs.max()) + pad + 1)
        mask_u8 = mask_u8[y0:y1, x0:x1]
        if rgb_patch is not None and rgb_patch.shape[:2] != mask_u8.shape[:2]:
            rgb_patch = rgb_patch[y0:y1, x0:x1]

    return mask_u8, rgb_patch


def _render_lm3_silhouette_card(
    mask_u8: np.ndarray,
    rgb_patch: Optional[np.ndarray],
    accent_hex: str = "#0072B2",
    pixels_per_mm: Optional[float] = None,
    is_calibrated: bool = False,
    card_size: int = 320,
    prefer_rgb: bool = True,
) -> np.ndarray:
    """
    Render a crisp square RGB card (white background) showing either the LeafMachine3
    masked herbarium RGB cutout or crisp botanical silhouette, overlaid with the homologized
    tip-up midrib vector, `lamina_base` attachment point, and a 10 mm calibrated scale bar.
    """
    card = np.full((card_size, card_size, 3), 255, dtype=np.uint8)
    h, w = mask_u8.shape[:2]
    if h == 0 or w == 0 or int(np.count_nonzero(mask_u8)) == 0:
        return card

    target_span = int(card_size * 0.80)
    scale = float(target_span) / float(max(h, w, 1))
    nh = max(4, int(round(h * scale)))
    nw = max(4, int(round(w * scale)))
    m_res = cv2.resize(mask_u8, (nw, nh), interpolation=cv2.INTER_NEAREST)
    ox = (card_size - nw) // 2
    oy = (card_size - nh) // 2

    acc_rgb = tuple(int(accent_hex.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
    fg = m_res > 0

    if prefer_rgb and rgb_patch is not None and rgb_patch.shape[:2] == (h, w):
        rgb_res = cv2.resize(rgb_patch, (nw, nh), interpolation=cv2.INTER_AREA)
        patch_canvas = card[oy : oy + nh, ox : ox + nw]
        patch_canvas[fg] = rgb_res[fg]
    else:
        tint = tuple(int(0.22 * c + 0.78 * 235) for c in acc_rgb)
        patch_canvas = card[oy : oy + nh, ox : ox + nw]
        patch_canvas[fg] = tint

    # Crisp LeafMachine3 contour outline
    cnts, _ = cv2.findContours((fg.astype(np.uint8) * 255), cv2.RETR_CCOMP, cv2.CHAIN_APPROX_SIMPLE)
    if cnts:
        shifted = [c + np.array([[[ox, oy]]], dtype=np.int32) for c in cnts]
        cv2.drawContours(card, shifted, -1, acc_rgb, thickness=2, lineType=cv2.LINE_AA)

    # Tip-up midrib axis and basal lamina attachment point (p_base at bottom pole of upright mask)
    ys, xs = np.where(fg)
    if len(ys) >= 10:
        y_top, y_bot = int(ys.min()), int(ys.max())
        top_band = xs[ys <= y_top + max(2, int(0.06 * (y_bot - y_top + 1)))]
        bot_band = xs[ys >= y_bot - max(2, int(0.06 * (y_bot - y_top + 1)))]
        x_tip = int(round(float(np.median(top_band)))) + ox if len(top_band) else (ox + nw // 2)
        x_base = int(round(float(np.median(bot_band)))) + ox if len(bot_band) else (ox + nw // 2)
        y_tip = y_top + oy
        y_base = y_bot + oy

        # Draw dashed midrib line and hollow p_base / p_apex markers
        n_seg = 14
        for s_i in range(0, n_seg, 2):
            t0, t1 = s_i / n_seg, min(1.0, (s_i + 1) / n_seg)
            p0 = (int(round(x_base + t0 * (x_tip - x_base))), int(round(y_base + t0 * (y_tip - y_base))))
            p1 = (int(round(x_base + t1 * (x_tip - x_base))), int(round(y_base + t1 * (y_tip - y_base))))
            cv2.line(card, p0, p1, (75, 85, 99), 1, cv2.LINE_AA)
        cv2.circle(card, (x_base, y_base), 5, (213, 94, 0), 2, cv2.LINE_AA)
        cv2.circle(card, (x_tip, y_tip), 4, acc_rgb, 1, cv2.LINE_AA)

    # Physical 10 mm scale bar on ruler-calibrated specimens
    if is_calibrated and pixels_per_mm is not None and np.isfinite(pixels_per_mm) and float(pixels_per_mm) > 0:
        bar_px = int(round(10.0 * float(pixels_per_mm) * scale))
        bar_mm = 10
        if bar_px > int(card_size * 0.45):
            bar_px = int(round(5.0 * float(pixels_per_mm) * scale))
            bar_mm = 5
        if 12 <= bar_px <= int(card_size * 0.60):
            bx0 = card_size - bar_px - 14
            bx1 = card_size - 14
            by = card_size - 14
            cv2.line(card, (bx0, by), (bx1, by), (30, 30, 30), 2, cv2.LINE_AA)
            cv2.line(card, (bx0, by - 3), (bx0, by + 3), (30, 30, 30), 1, cv2.LINE_AA)
            cv2.line(card, (bx1, by - 3), (bx1, by + 3), (30, 30, 30), 1, cv2.LINE_AA)
            cv2.putText(
                card,
                f"{bar_mm} mm",
                (max(6, bx0 - 2), by - 6),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.36,
                (30, 30, 30),
                1,
                cv2.LINE_AA,
            )

    return card


def _render_lm3_ect_rgb(mask_u8: np.ndarray, n_dirs: int = 160, size: int = 320) -> Optional[np.ndarray]:
    """Delegate polar radial ECT + contour rendering directly to LeafMachine3 `render_radial_ect_overlay`."""
    _ensure_ect_module_loaded()
    from leafmachine3.core.ect_compute import compute_ect
    from leafmachine3.reporting.ect_viz import render_radial_ect_overlay, visual_log

    ect_res = compute_ect(mask_u8, num_dirs=int(n_dirs), bound_radius=1.0, want_simple=True)
    if ect_res is None or getattr(ect_res, "ect_matrix", None) is None:
        return None
    bgr = render_radial_ect_overlay(
        visual_log(ect_res.ect_matrix),
        ect_res.thetas,
        ect_res.thresholds,
        ect_res.bound_radius,
        getattr(ect_res, "outline_norm", None),
        cmap="viridis",
    )
    if bgr is None or bgr.size == 0:
        return None
    rgb = cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)
    if rgb.shape[:2] != (size, size):
        rgb = cv2.resize(rgb, (size, size), interpolation=cv2.INTER_AREA)
    return rgb


def export_taxon_ect_exemplars(
    extracted_df: Optional[pd.DataFrame] = None,
    figures_dir: Union[str, Path] = "outputs/figures",
    zone: str = "basal",
    n_dirs: int = 160,
    max_taxa: Optional[int] = None,
    vouchers_csv: Optional[Union[str, Path]] = None,
    tables_dir: Optional[Union[str, Path]] = None,
    db_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Dict[str, Any]]:
    """
    Export authentic LeafMachine3 radial ECT overlays (`ect_overlays/<zone>_<slug>.png`)
    and native LeafMachine3 silhouette/RGB exemplar badges (`lm3_assets/<zone>_<slug>_*.png`)
    for ALL active taxa/ranks (unlimited by default when `max_taxa=None`).
    """
    fig_dir = Path(figures_dir).resolve()
    overlays_dir = fig_dir / "ect_overlays"
    lm3_assets_dir = fig_dir / "lm3_assets"
    overlays_dir.mkdir(parents=True, exist_ok=True)
    lm3_assets_dir.mkdir(parents=True, exist_ok=True)
    zone_clean = str(zone).lower().strip()

    if extracted_df is None or len(extracted_df) == 0:
        cand_tables: List[Path] = []
        if tables_dir is not None:
            t_p = Path(tables_dir)
            cand_tables.extend([t_p / "extracted_leaf_manifest.csv", t_p / "leaf_efa_harmonics_individual.csv", t_p / "leaf_efa_harmonics.csv"])
        for base in (fig_dir.parent.parent / "data" / "tables", fig_dir.parent / "data" / "tables", PROJECT_ROOT / "data" / "tables"):
            cand_tables.extend([base / "extracted_leaf_manifest.csv", base / "leaf_efa_harmonics_individual.csv", base / "leaf_efa_harmonics.csv"])
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
    if vouchers_csv is not None and Path(vouchers_csv).exists() and "catalogNumber" in df_work.columns:
        try:
            v_df = pd.read_csv(vouchers_csv, dtype=str)
            if "catalogNumber" in v_df.columns:
                for col_v in ("species_standardized", "scientificName", "infraspecific_taxon", "image_path"):
                    if col_v in v_df.columns and col_v not in df_work.columns:
                        v_map = dict(zip(v_df["catalogNumber"].astype(str).str.strip(), v_df[col_v]))
                        df_work[col_v] = df_work["catalogNumber"].astype(str).str.strip().map(v_map)
        except Exception:
            pass

    if "species_clean" not in df_work.columns:
        sp_col = next(
            (c for c in ("species_standardized", "scientificName", "species_name", "species_raw", "taxon") if c in df_work.columns),
            None,
        )
        inf_col = "infraspecific_taxon" if "infraspecific_taxon" in df_work.columns else None
        if sp_col is not None:
            if inf_col is not None:
                df_work["species_clean"] = [
                    _standardize_packera_species_py(s, inf) for s, inf in zip(df_work[sp_col], df_work[inf_col])
                ]
            else:
                df_work["species_clean"] = df_work[sp_col].apply(_standardize_packera_species_py)
        else:
            df_work["species_clean"] = "Unassigned_Packera"

    if "leaf_zone" in df_work.columns:
        lz_norm = df_work["leaf_zone"].astype(str).str.strip().str.lower()
        if zone_clean == "cauline":
            df_zone = df_work[lz_norm == "cauline"].copy()
        else:
            df_zone = df_work[lz_norm.isin({"basal", "", "nan", "none"})].copy()
            if len(df_zone) == 0:
                df_zone = df_work.copy()
    else:
        df_zone = df_work.iloc[0:0].copy() if zone_clean == "cauline" else df_work.copy()

    if len(df_zone) == 0:
        manifest_path = overlays_dir / f"manifest_{zone_clean}.json"
        with open(manifest_path, "w", encoding="utf-8") as f_m:
            json.dump({}, f_m, indent=2)
        logger.info(f"Exported 0 authentic LeafMachine3 radial ECT overlays to {overlays_dir}")
        return {}

    sp_counts = df_zone["species_clean"].value_counts()
    valid_sp = [s for s in sp_counts.index if s not in {"Unassigned_Packera", "Unknown", ""} and pd.notna(s)]
    target_sp = valid_sp[: int(max_taxa)] if (max_taxa is not None and int(max_taxa) > 0) else valid_sp
    style_map = get_combinatorial_taxon_style(target_sp)

    manifest: Dict[str, Dict[str, Any]] = {}
    for sp_name in target_sp:
        sp_df = df_zone[df_zone["species_clean"] == sp_name]
        if len(sp_df) == 0:
            continue

        # Filter to clean closed/unclipped masks when available
        valid_mask_flags = [
            _is_valid_smooth_closed_mask(r.get("mask_path", ""), fig_dir) for _, r in sp_df.iterrows()
        ]
        sp_clean_df = sp_df[valid_mask_flags] if any(valid_mask_flags) else sp_df

        has_ect_pcs = "ECT_PC1" in sp_clean_df.columns and "ECT_PC2" in sp_clean_df.columns
        best_row = None
        if has_ect_pcs:
            pc1 = pd.to_numeric(sp_clean_df["ECT_PC1"], errors="coerce")
            pc2 = pd.to_numeric(sp_clean_df["ECT_PC2"], errors="coerce")
            valid_pcs = pc1.notna() & pc2.notna()
            if valid_pcs.sum() >= 2:
                sub_pc = sp_clean_df[valid_pcs]
                dist_pc = np.sqrt(
                    (pc1[valid_pcs] - float(pc1[valid_pcs].median())) ** 2
                    + (pc2[valid_pcs] - float(pc2[valid_pcs].median())) ** 2
                )
                is_type_s = (
                    sub_pc["is_type_specimen"].apply(_parse_bool_val)
                    if "is_type_specimen" in sub_pc.columns
                    else pd.Series(False, index=sub_pc.index)
                )
                type_hits = sub_pc[is_type_s & (dist_pc <= 1.5 * float(dist_pc.median() + 1e-6))]
                best_row = type_hits.iloc[0] if len(type_hits) > 0 else sub_pc.loc[dist_pc.idxmin()]

        if best_row is None:
            sol_s = pd.to_numeric(sp_clean_df.get("solidity", 0.88), errors="coerce").fillna(0.88)
            ar_s = pd.to_numeric(sp_clean_df.get("aspect_ratio", 2.0), errors="coerce").fillna(2.0)
            d_geom = (sol_s - float(sol_s.median())) ** 2 + (ar_s - float(ar_s.median())) ** 2
            best_row = sp_clean_df.loc[d_geom.idxmin()]

        cat_num = str(best_row.get("catalogNumber", "Unknown"))
        leaf_id = int(best_row.get("leaf_id", 1)) if pd.notna(best_row.get("leaf_id", 1)) else 1
        taxon_slug = re.sub(r"[^A-Za-z0-9_]+", "_", sp_name.lower()).strip("_")
        out_ect_png = overlays_dir / f"{zone_clean}_{taxon_slug}.png"
        out_sil_png = lm3_assets_dir / f"{zone_clean}_{taxon_slug}_silhouette.png"
        out_rgb_png = lm3_assets_dir / f"{zone_clean}_{taxon_slug}_cutout.png"

        mask_u8, rgb_patch = _load_mask_or_contour_u8(best_row, fig_dir)
        wrote_ect = False
        if mask_u8 is not None and int(np.count_nonzero(mask_u8)) > 0:
            ect_rgb = _render_lm3_ect_rgb(mask_u8, n_dirs=n_dirs, size=320)
            if ect_rgb is not None:
                cv2.imwrite(str(out_ect_png), cv2.cvtColor(ect_rgb, cv2.COLOR_RGB2BGR))
                wrote_ect = True

            acc_col = style_map.get(sp_name, ("#0072B2", "o"))[0]
            px_mm = pd.to_numeric(best_row.get("pixels_per_mm", np.nan), errors="coerce")
            is_cal = _parse_bool_val(best_row.get("is_calibrated_mm", False)) and str(
                best_row.get("scale_source", "")
            ).strip().lower() not in {"mp_estimated", "predicted_from_megapixels", "uncalibrated"}
            sil_card = _render_lm3_silhouette_card(
                mask_u8,
                None,
                accent_hex=acc_col,
                pixels_per_mm=float(px_mm) if pd.notna(px_mm) else None,
                is_calibrated=is_cal,
                card_size=320,
                prefer_rgb=False,
            )
            cv2.imwrite(str(out_sil_png), cv2.cvtColor(sil_card, cv2.COLOR_RGB2BGR))
            if rgb_patch is not None:
                cut_card = _render_lm3_silhouette_card(
                    mask_u8,
                    rgb_patch,
                    accent_hex=acc_col,
                    pixels_per_mm=float(px_mm) if pd.notna(px_mm) else None,
                    is_calibrated=is_cal,
                    card_size=320,
                    prefer_rgb=True,
                )
                cv2.imwrite(str(out_rgb_png), cv2.cvtColor(cut_card, cv2.COLOR_RGB2BGR))
        elif out_ect_png.exists() and out_ect_png.stat().st_size > 500:
            wrote_ect = True

        if not wrote_ect:
            continue

        comp_val = best_row.get("ect_lobe_complexity", best_row.get("cauline_ect_lobe_complexity", np.nan))
        comp_float = float(comp_val) if pd.notna(comp_val) else np.nan
        manifest[sp_name] = {
            "image_path": str(out_ect_png),
            "silhouette_path": str(out_sil_png) if out_sil_png.exists() else None,
            "cutout_path": str(out_rgb_png) if out_rgb_png.exists() else None,
            "catalogNumber": cat_num,
            "leaf_id": leaf_id,
            "taxon": sp_name,
            "taxon_slug": taxon_slug,
            "zone": zone_clean,
            "is_type": bool(_parse_bool_val(best_row.get("is_type_specimen", False))),
            "ect_lobe_complexity": round(comp_float, 4) if np.isfinite(comp_float) else None,
            "aspect_ratio": round(float(best_row.get("aspect_ratio", np.nan)), 4)
            if pd.notna(best_row.get("aspect_ratio"))
            else None,
            "solidity": round(float(best_row.get("solidity", np.nan)), 4)
            if pd.notna(best_row.get("solidity"))
            else None,
        }

    manifest_path = overlays_dir / f"manifest_{zone_clean}.json"
    with open(manifest_path, "w", encoding="utf-8") as f_m:
        json.dump(manifest, f_m, indent=2)
    lm3_manifest_path = lm3_assets_dir / f"manifest_{zone_clean}.json"
    with open(lm3_manifest_path, "w", encoding="utf-8") as f_lm:
        json.dump(manifest, f_lm, indent=2)

    logger.info(f"Exported {len(manifest)} authentic LeafMachine3 radial ECT overlays to {overlays_dir}")
    return manifest


export_taxon_lm3_exemplars = export_taxon_ect_exemplars


def export_lm3_sheet_diagnostic_overlays(
    extracted_df: pd.DataFrame,
    repro_df: Optional[pd.DataFrame],
    pheno_df: Optional[pd.DataFrame],
    figures_dir: Union[str, Path],
    vouchers_df: Optional[pd.DataFrame] = None,
    max_sheets: int = 24,
) -> List[Path]:
    """
    Export capped whole-sheet LeafMachine3 diagnostic overlays (`lm3_sheet_overlays/*.png`)
    at 300 DPI for priority sheets (Nomenclatural Types, DBSCAN multi-plant ramets K >= 2,
    basal+cauline heteroblastic sheets, and ruler-calibrated exemplars) whenever full scans exist.
    """
    out_dir = Path(figures_dir)
    sheet_dir = out_dir / "lm3_sheet_overlays"
    if extracted_df is None or len(extracted_df) == 0 or "catalogNumber" not in extracted_df.columns:
        return []

    df = extracted_df.copy()
    df["cat_clean"] = df["catalogNumber"].astype(str).str.strip()

    v_map_img: Dict[str, str] = {}
    v_map_tier: Dict[str, str] = {}
    v_map_type: Dict[str, bool] = {}
    if vouchers_df is not None and len(vouchers_df) > 0 and "catalogNumber" in vouchers_df.columns:
        for _, vr in vouchers_df.iterrows():
            c_k = str(vr.get("catalogNumber", "")).strip()
            if not c_k:
                continue
            if pd.notna(vr.get("image_path")):
                v_map_img[c_k] = str(vr.get("image_path")).strip()
            if pd.notna(vr.get("determiner_tier")):
                v_map_tier[c_k] = str(vr.get("determiner_tier")).strip()
            v_map_type[c_k] = _parse_bool_val(vr.get("is_type_specimen", False)) or (
                str(vr.get("type_status", "")).strip().lower() not in {"", "none", "nan", "na", "not a type"}
            )

    pheno_map: Dict[str, str] = {}
    if pheno_df is not None and len(pheno_df) > 0 and "catalogNumber" in pheno_df.columns:
        st_col = next(
            (c for c in ("voucher_phenological_state", "phenological_state", "pheno_state") if c in pheno_df.columns),
            None,
        )
        if st_col is not None:
            for _, pr in pheno_df.iterrows():
                pheno_map[str(pr["catalogNumber"]).strip()] = str(pr.get(st_col, "")).strip()

    # Score candidate sheets so Types, multi-plant DBSCAN ramets, and heteroblastic sheets come first
    candidates: List[Tuple[int, str]] = []
    for cat_id, grp in df.groupby("cat_clean", sort=False):
        if not cat_id:
            continue
        img_p = _resolve_asset_path(v_map_img.get(cat_id, grp.iloc[0].get("image_path", "")), out_dir)
        if img_p is None:
            for ext in (".webp", ".jpg", ".jpeg", ".png", ".tif"):
                cand_raw = _resolve_asset_path(f"raw_vouchers/{cat_id}{ext}", out_dir)
                if cand_raw is not None:
                    img_p = cand_raw
                    break
        if img_p is None or not img_p.exists():
            continue

        is_t = v_map_type.get(cat_id, False) or any(
            _parse_bool_val(x) for x in grp.get("is_type_specimen", pd.Series(dtype=object))
        )
        pids = set(
            int(x)
            for x in pd.to_numeric(grp.get("plant_individual_id", pd.Series([1])), errors="coerce").dropna()
            if int(x) > 0
        )
        zones = set(str(z).strip().lower() for z in grp.get("leaf_zone", pd.Series(dtype=object)).dropna())
        is_cal = any(_parse_bool_val(x) for x in grp.get("is_calibrated_mm", pd.Series(dtype=object)))
        prio = (100 if is_t else 0) + (50 if len(pids) >= 2 else 0) + (25 if {"basal", "cauline"}.issubset(zones) else 0) + (10 if is_cal else 0) + min(len(grp), 9)
        candidates.append((prio, cat_id))

    if not candidates:
        return []

    candidates.sort(key=lambda item: (-item[0], item[1]))
    selected_cats = [cat_id for _, cat_id in candidates[: int(max_sheets)]]
    sheet_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    for cat_id in selected_cats:
        grp = df[df["cat_clean"] == cat_id]
        img_p = _resolve_asset_path(v_map_img.get(cat_id, grp.iloc[0].get("image_path", "")), out_dir)
        if img_p is None:
            for ext in (".webp", ".jpg", ".jpeg", ".png", ".tif"):
                cand_raw = _resolve_asset_path(f"raw_vouchers/{cat_id}{ext}", out_dir)
                if cand_raw is not None:
                    img_p = cand_raw
                    break
        if img_p is None or not img_p.exists():
            continue
        bgr = cv2.imread(str(img_p), cv2.IMREAD_COLOR)
        if bgr is None or bgr.size == 0:
            continue

        H_orig, W_orig = bgr.shape[:2]
        max_dim = 1800
        disp_scale = min(1.0, float(max_dim) / float(max(H_orig, W_orig, 1)))
        if disp_scale < 1.0:
            canvas = cv2.resize(
                bgr,
                (max(1, int(round(W_orig * disp_scale))), max(1, int(round(H_orig * disp_scale)))),
                interpolation=cv2.INTER_AREA,
            )
        else:
            canvas = bgr.copy()

        # Draw per-ramet envelopes when bbox coordinates exist
        b_cols = ("bbox_x1", "bbox_y1", "bbox_x2", "bbox_y2")
        has_bboxes = all(c in grp.columns for c in b_cols)
        ramet_boxes: Dict[int, List[Tuple[int, int, int, int]]] = {}

        for _, lrow in grp.iterrows():
            if not has_bboxes or any(pd.isna(lrow.get(c)) for c in b_cols):
                continue
            x1 = int(round(float(lrow["bbox_x1"]) * disp_scale))
            y1 = int(round(float(lrow["bbox_y1"]) * disp_scale))
            x2 = int(round(float(lrow["bbox_x2"]) * disp_scale))
            y2 = int(round(float(lrow["bbox_y2"]) * disp_scale))
            lz = str(lrow.get("leaf_zone", "basal")).strip().lower()
            gates_ok = _parse_bool_val(lrow.get("lm3_gates_pass", True), default=True)
            # BGR colors: Basal=#B27200 (blue), Cauline=#739E00 (green), Flagged=#888888
            col_bgr = (178, 114, 0) if lz != "cauline" else (115, 158, 0)
            if not gates_ok:
                col_bgr = (130, 130, 130)
            cv2.rectangle(canvas, (x1, y1), (x2, y2), col_bgr, 2, cv2.LINE_AA)
            pid = int(pd.to_numeric(lrow.get("plant_individual_id", 1), errors="coerce") or 1)
            ramet_boxes.setdefault(pid, []).append((x1, y1, x2, y2))
            lbl = f"L{int(lrow.get('leaf_id', 1))}:{lz[:3].upper()}"
            cv2.putText(canvas, lbl, (x1 + 3, max(14, y1 - 4)), cv2.FONT_HERSHEY_SIMPLEX, 0.42, col_bgr, 1, cv2.LINE_AA)

        if len(ramet_boxes) >= 2:
            for pid, boxes in sorted(ramet_boxes.items()):
                rx1 = max(2, min(b[0] for b in boxes) - 10)
                ry1 = max(2, min(b[1] for b in boxes) - 24)
                rx2 = min(canvas.shape[1] - 2, max(b[2] for b in boxes) + 10)
                ry2 = min(canvas.shape[0] - 2, max(b[3] for b in boxes) + 10)
                cv2.rectangle(canvas, (rx1, ry1), (rx2, ry2), (0, 94, 213), 2, cv2.LINE_AA)
                tag_txt = f"Ramet #{pid}"
                (tw, th), _ = cv2.getTextSize(tag_txt, cv2.FONT_HERSHEY_SIMPLEX, 0.48, 1)
                by_top = max(4, ry1 - th - 10)
                cv2.rectangle(canvas, (rx1, by_top), (rx1 + tw + 12, by_top + th + 8), (0, 94, 213), -1)
                cv2.putText(
                    canvas,
                    tag_txt,
                    (rx1 + 6, by_top + th + 2),
                    cv2.FONT_HERSHEY_SIMPLEX,
                    0.48,
                    (255, 255, 255),
                    1,
                    cv2.LINE_AA,
                )

        # Two-line non-clipping header banner summarizing voucher metadata & LM3 calibration
        sp_raw = grp.iloc[0].get("species_standardized", grp.iloc[0].get("scientificName", ""))
        sp_lbl = _standardize_packera_species_py(sp_raw, grp.iloc[0].get("infraspecific_taxon", None)).replace("×", "x")
        px_mm = pd.to_numeric(grp.iloc[0].get("pixels_per_mm", np.nan), errors="coerce")
        sc_src = str(grp.iloc[0].get("scale_source", "uncalibrated")).strip()
        cal_txt = f"{float(px_mm):.2f} px/mm ({sc_src})" if pd.notna(px_mm) and float(px_mm) > 0 else sc_src
        ph_txt = pheno_map.get(cat_id, "unscored")
        tier_txt = v_map_tier.get(cat_id, str(grp.iloc[0].get("determiner_tier", "")))
        banner_h = 64
        banner = np.full((banner_h, canvas.shape[1], 3), 248, dtype=np.uint8)
        cv2.line(banner, (0, banner_h - 1), (canvas.shape[1], banner_h - 1), (190, 195, 202), 1)
        line1 = f"{cat_id}  |  {sp_lbl}  |  Authority: {tier_txt}"
        line2 = f"Leaves: {len(grp)} (DBSCAN Ramets: {max(1, len(ramet_boxes))})  |  Phenology: {ph_txt}  |  Scale: {cal_txt}"
        fscale = min(0.50, max(0.36, canvas.shape[1] / 1750.0))
        cv2.putText(banner, line1, (14, 25), cv2.FONT_HERSHEY_SIMPLEX, fscale, (20, 24, 28), 1, cv2.LINE_AA)
        cv2.putText(banner, line2, (14, 50), cv2.FONT_HERSHEY_SIMPLEX, fscale * 0.94, (55, 65, 81), 1, cv2.LINE_AA)
        composed = np.vstack([banner, canvas])

        out_overlay = sheet_dir / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', cat_id)}_lm3_sheet_overlay.png"
        tmp_ov = sheet_dir / f".{out_overlay.stem}_{os.getpid()}.png"
        cv2.imwrite(str(tmp_ov), composed)
        tmp_ov.replace(out_overlay)
        written.append(out_overlay)

    if written:
        logger.info(f"Exported {len(written)} whole-sheet LeafMachine3 diagnostic overlays to {sheet_dir}")
    return written


def _repel_medoid_annotations(
    ax: Any,
    points_xy_label_col: List[Tuple[float, float, str, str]],
    xlim: Tuple[float, float],
    ylim: Tuple[float, float],
) -> None:
    """Place direct centroid/medoid labels inside plot bounds with 2D iterative box repulsion."""
    if not points_xy_label_col:
        return
    x_span = max(xlim[1] - xlim[0], 1e-6)
    y_span = max(ylim[1] - ylim[0], 1e-6)
    n_lbl = len(points_xy_label_col)
    pos = np.array([[p[0], p[1] + 0.035 * y_span] for p in points_xy_label_col], dtype=float)
    widths = np.array([min(0.42, max(0.16, 0.0125 * len(re.sub(r"[$\\{}]|it", "", p[2])))) * x_span for p in points_xy_label_col])
    box_h = 0.072 * y_span

    for _ in range(80):
        for i in range(n_lbl):
            pos[i, 0] = np.clip(pos[i, 0], xlim[0] + 0.52 * widths[i], xlim[1] - 0.52 * widths[i])
            pos[i, 1] = np.clip(pos[i, 1], ylim[0] + 0.65 * box_h, ylim[1] - 0.65 * box_h)
        for i in range(n_lbl):
            for j in range(i + 1, n_lbl):
                dx = pos[i, 0] - pos[j, 0]
                dy = pos[i, 1] - pos[j, 1]
                req_x = 0.54 * (widths[i] + widths[j])
                req_y = 1.05 * box_h
                if abs(dx) < req_x and abs(dy) < req_y:
                    push_y = 0.52 * (req_y - abs(dy)) * (1.0 if dy >= 0 else -1.0)
                    pos[i, 1] += push_y
                    pos[j, 1] -= push_y
                    push_x = 0.18 * (req_x - abs(dx)) * (1.0 if dx >= 0 else -1.0)
                    pos[i, 0] += push_x
                    pos[j, 0] -= push_x

    for idx, (x0, y0, lbl, col) in enumerate(points_xy_label_col):
        xk = float(np.clip(pos[idx, 0], xlim[0] + 0.52 * widths[idx], xlim[1] - 0.52 * widths[idx]))
        yk = float(np.clip(pos[idx, 1], ylim[0] + 0.65 * box_h, ylim[1] - 0.65 * box_h))
        ax.annotate(
            lbl,
            xy=(x0, y0),
            xytext=(xk, yk),
            ha="center",
            va="center",
            fontsize=7.2,
            color="#1f2937",
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=col, alpha=0.90, linewidth=0.75),
            arrowprops=dict(arrowstyle="-", color=col, lw=0.65, shrinkA=2, shrinkB=3),
            zorder=8,
        )


def _draw_solidity_vs_ar_zone_panel(
    ax: Any,
    zone_df: pd.DataFrame,
    active_taxa: List[str],
    style_map: Dict[str, Tuple[str, str]],
    panel_title: str,
    show_legend: bool = True,
) -> None:
    """Render hollow-marker Solidity vs. Aspect Ratio panel with reniform L/W < 1.0 dynamic scaling."""
    sol_vals = pd.to_numeric(zone_df["solidity_clean"], errors="coerce")
    ar_vals = pd.to_numeric(zone_df["ar_clean"], errors="coerce")
    ok = sol_vals.notna() & ar_vals.notna() & np.isfinite(sol_vals) & np.isfinite(ar_vals)
    sub = zone_df.loc[ok].copy()
    if len(sub) == 0:
        ax.text(0.5, 0.5, "No leaves in zone", ha="center", va="center", transform=ax.transAxes, fontsize=9, color="#6b7280")
        ax.set_title(panel_title, fontsize=10, fontweight="bold", loc="left")
        return

    x_min = max(0.15, float(sub["solidity_clean"].min()) - 0.05)
    x_max = min(1.02, float(sub["solidity_clean"].max()) + 0.03)
    if x_max - x_min < 0.20:
        x_min = max(0.15, x_max - 0.25)

    min_ar = float(sub["ar_clean"].min())
    max_ar = float(sub["ar_clean"].max())
    # Preserve reniform/orbicular blades (L/W < 1.0) without clipping at 1.0
    y_min = max(0.30, min(0.92, min_ar - 0.12)) if min_ar < 1.0 else max(0.75, min_ar - 0.18)
    y_max = max(y_min + 0.8, min(10.0, max_ar + 0.35))

    # Subtle morphological regime guides (Dissected vs. Entire/Undivided; Reniform L/W = 1.0)
    if x_min < 0.75 < x_max:
        ax.axvspan(x_min, 0.75, color="#f3f4f6", alpha=0.55, zorder=0)
        ax.axvline(0.75, color="#9ca3af", linestyle="--", linewidth=0.8, zorder=1)
        ax.text(
            x_min + 0.01 * (x_max - x_min),
            y_max - 0.05 * (y_max - y_min),
            "Lobed / Pinnatifid",
            fontsize=7.5,
            color="#4b5563",
            va="top",
        )
    if y_min < 1.0 < y_max:
        ax.axhline(1.0, color="#9ca3af", linestyle=":", linewidth=0.8, zorder=1)

    n_pts = len(sub)
    alpha_pt = 0.82 if n_pts <= 150 else (0.68 if n_pts <= 600 else 0.52)
    medoid_labels: List[Tuple[float, float, str, str]] = []

    for tx in active_taxa:
        tx_rows = sub[sub["species_clean"] == tx]
        if len(tx_rows) == 0:
            continue
        col, mk = style_map[tx]
        lbl_math = f"{format_icn_taxon_mathtext(tx, abbreviate_genus=(len(active_taxa) > 6))} (n={len(tx_rows)})"
        ax.scatter(
            tx_rows["solidity_clean"],
            tx_rows["ar_clean"],
            s=36,
            marker=mk,
            facecolors="none",
            edgecolors=col,
            linewidths=1.0,
            alpha=alpha_pt,
            label=lbl_math,
            zorder=4,
        )
        med_x = float(tx_rows["solidity_clean"].median())
        med_y = float(tx_rows["ar_clean"].median())
        medoid_labels.append((med_x, med_y, format_icn_taxon_mathtext(tx, abbreviate_genus=True), col))

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel("Lamina Solidity (Area / Convex Hull Area)", fontsize=9)
    ax.set_ylabel("Anatomical Blade Aspect Ratio (L / W)", fontsize=9)
    ax.set_title(panel_title, fontsize=10, fontweight="bold", loc="left")
    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.55)

    if len(active_taxa) > 10:
        _repel_medoid_annotations(ax, medoid_labels, (x_min, x_max), (y_min, y_max))
    elif show_legend and len(active_taxa) >= 1:
        ncol_leg = 2 if len(active_taxa) > 5 else 1
        ax.legend(
            loc="best",
            fontsize=7.3,
            frameon=True,
            facecolor="white",
            framealpha=0.88,
            edgecolor="#d1d5db",
            ncol=ncol_leg,
        )


def generate_step02_publication_figures(
    extracted_df: pd.DataFrame,
    failed_df: pd.DataFrame,
    repro_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    figures_dir: Union[str, Path],
    vouchers_df: Optional[pd.DataFrame] = None,
) -> Dict[str, Path]:
    """
    Generate the two canonical Step 02 300-DPI PNG plates + native LM3 asset caches:
      1. `fig02a_lm3_organ_and_ect_exemplars.png` (`saved_paths["fig02a_png"]`)
      2. `fig02b_segmentation_qc_and_phenophase.png` (`saved_paths["fig02b_png"]`)
      3. `ect_overlays/` & `lm3_assets/` (`export_taxon_ect_exemplars` across all taxa)
      4. `lm3_sheet_overlays/` (`export_lm3_sheet_diagnostic_overlays` for priority scans)
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

    _step02_mod = sys.modules.get("scripts.pipeline.02_segment_and_extract")
    _save_fig_fn = (
        getattr(_step02_mod, "_atomic_save_matplotlib_fig", _atomic_save_matplotlib_fig)
        if _step02_mod is not None
        else _atomic_save_matplotlib_fig
    )

    out_dir = Path(figures_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    saved_paths: Dict[str, Path] = {}

    if extracted_df is None or len(extracted_df) == 0:
        return saved_paths

    # Clean up any legacy multi-page/duplicate Step 02 PDF/PNG remnants from prior pipeline versions
    for legacy_stem in (
        "fig02a_segmentation_qc_distributions",
        "fig02b_silhouette_montage",
        "fig02b_silhouette_montage_supplement",
        "fig02c_cauline_silhouette_montage",
        "fig02c_cauline_silhouette_montage_supplement",
        "fig02d_basal_ect_radial_montage",
        "fig02e_cauline_ect_radial_montage",
    ):
        for ext in (".pdf", ".png"):
            rem = out_dir / f"{legacy_stem}{ext}"
            if rem.exists():
                rem.unlink(missing_ok=True)

    df = extracted_df.copy()
    if "catalogNumber" not in df.columns:
        df["catalogNumber"] = [f"SPEC_{i:04d}" for i in range(len(df))]
    df["cat_clean"] = df["catalogNumber"].astype(str).str.strip()
    active_cats: Set[str] = set(df["cat_clean"].unique()) - {"", "nan", "None"}

    # Enrich taxon, infraspecific rank, ruler scale, and scan paths from curated_vouchers.csv if available
    if vouchers_df is None:
        v_cand = _discover_downstream_table("curated_vouchers.csv", out_dir, active_cats)
        if v_cand is not None:
            try:
                vouchers_df = pd.read_csv(v_cand, low_memory=False)
            except Exception:
                vouchers_df = None

    if vouchers_df is not None and len(vouchers_df) > 0 and "catalogNumber" in vouchers_df.columns:
        v_sub = vouchers_df.drop_duplicates(subset=["catalogNumber"]).copy()
        v_sub["cat_clean"] = v_sub["catalogNumber"].astype(str).str.strip()
        # Restrict active species strictly to vouchers present in curated_vouchers.csv when matched
        matched_vcats = set(v_sub["cat_clean"]) & active_cats
        if matched_vcats:
            df = df[df["cat_clean"].isin(matched_vcats)].copy()
            active_cats = matched_vcats
        for v_col in ("species_standardized", "scientificName", "infraspecific_taxon", "determiner_tier", "is_type_specimen", "image_path"):
            if v_col in v_sub.columns and v_col not in df.columns:
                v_map = dict(zip(v_sub["cat_clean"], v_sub[v_col]))
                df[v_col] = df["cat_clean"].map(v_map)

    sp_col = next(
        (c for c in ("species_standardized", "scientificName", "species_name", "species_raw", "taxon") if c in df.columns),
        None,
    )
    inf_col = "infraspecific_taxon" if "infraspecific_taxon" in df.columns else None
    if sp_col is not None:
        if inf_col is not None:
            df["species_clean"] = [
                _standardize_packera_species_py(s, inf) for s, inf in zip(df[sp_col], df[inf_col])
            ]
        else:
            df["species_clean"] = df[sp_col].apply(_standardize_packera_species_py)
    else:
        df["species_clean"] = "Unassigned_Packera"

    df["solidity_clean"] = pd.to_numeric(df.get("solidity", df.get("sol_val", np.nan)), errors="coerce")
    df["ar_clean"] = pd.to_numeric(df.get("aspect_ratio", df.get("ar_val", np.nan)), errors="coerce")
    df["zone_clean"] = (
        df["leaf_zone"].astype(str).str.strip().str.lower()
        if "leaf_zone" in df.columns
        else pd.Series("basal", index=df.index)
    )
    df["zone_clean"] = df["zone_clean"].where(df["zone_clean"].isin({"basal", "cauline"}), "basal")

    # Discover downstream Step 03 contour outlier flags and Step 04 misidentification flags so exemplars are 100% clean
    outlier_leaf_keys: Set[Tuple[str, int]] = set()
    outlier_cats: Set[str] = set()
    misid_cats: Set[str] = set()

    ind_efa_p = _discover_downstream_table("leaf_efa_harmonics_individual.csv", out_dir, active_cats)
    if ind_efa_p is not None:
        try:
            ind_df = pd.read_csv(ind_efa_p, low_memory=False)
            if {"catalogNumber", "is_contour_outlier"}.issubset(ind_df.columns):
                for _, r_i in ind_df.iterrows():
                    if _parse_bool_val(r_i.get("is_contour_outlier", False)):
                        c_s = str(r_i["catalogNumber"]).strip()
                        if "leaf_id" in ind_df.columns and pd.notna(r_i.get("leaf_id")):
                            outlier_leaf_keys.add((c_s, int(float(r_i["leaf_id"]))))
                        else:
                            outlier_cats.add(c_s)
            if "leaf_quality_weight" in ind_df.columns and "leaf_quality_weight" not in df.columns and "catalogNumber" in ind_df.columns:
                lqw_map = {}
                for _, r_w in ind_df.iterrows():
                    c_w = str(r_w.get("catalogNumber", "")).strip()
                    lid_w = int(float(r_w.get("leaf_id", 1))) if pd.notna(r_w.get("leaf_id", 1)) else 1
                    w_v = pd.to_numeric(r_w.get("leaf_quality_weight", None), errors="coerce")
                    if c_w and pd.notna(w_v):
                        lqw_map[(c_w, lid_w)] = float(w_v)
                        lqw_map.setdefault((c_w, -1), float(w_v))
                if lqw_map:
                    df["leaf_quality_weight"] = [
                        lqw_map.get(
                            (str(r.get("cat_clean", "")).strip(), int(float(r.get("leaf_id", 1))) if pd.notna(r.get("leaf_id", 1)) else 1),
                            lqw_map.get((str(r.get("cat_clean", "")).strip(), -1), float("nan")),
                        )
                        for _, r in df.iterrows()
                    ]
        except Exception:
            pass

    otu_efa_p = _discover_downstream_table("leaf_efa_harmonics.csv", out_dir, active_cats)
    if otu_efa_p is not None:
        try:
            otu_df = pd.read_csv(otu_efa_p, low_memory=False)
            if {"catalogNumber", "is_contour_outlier"}.issubset(otu_df.columns):
                for _, r_o in otu_df.iterrows():
                    if _parse_bool_val(r_o.get("is_contour_outlier", False)):
                        outlier_cats.add(str(r_o["catalogNumber"]).strip())
        except Exception:
            pass

    mflags_p = _discover_downstream_table("morphometric_flags.csv", out_dir, active_cats)
    if mflags_p is not None:
        try:
            mf_df = pd.read_csv(mflags_p, low_memory=False)
            if "catalogNumber" in mf_df.columns:
                for _, r_m in mf_df.iterrows():
                    if _parse_bool_val(r_m.get("misidentification_flag", False)) or _parse_bool_val(
                        r_m.get("is_confident_label_noise", False)
                    ):
                        misid_cats.add(str(r_m["catalogNumber"]).strip())
        except Exception:
            pass

    def _is_exemplar_eligible(row: pd.Series) -> bool:
        c_s = str(row.get("cat_clean", "")).strip()
        lid = int(float(row.get("leaf_id", 1))) if pd.notna(row.get("leaf_id", 1)) else 1
        if c_s in misid_cats or c_s in outlier_cats or (c_s, lid) in outlier_leaf_keys:
            return False
        if "lm3_gates_pass" in row and pd.notna(row.get("lm3_gates_pass")):
            if not _parse_bool_val(row.get("lm3_gates_pass", True), default=True):
                return False
        return _is_valid_smooth_closed_mask(row.get("mask_path", ""), out_dir)

    df["exemplar_ok"] = [_is_exemplar_eligible(r) for _, r in df.iterrows()]

    sp_counts = df["species_clean"].value_counts()
    valid_taxa = [s for s in sp_counts.index if s not in {"Unassigned_Packera", "Unknown", ""} and pd.notna(s)]
    if not valid_taxa:
        valid_taxa = [str(s) for s in sp_counts.index if str(s).strip()]
    active_taxa = sorted(valid_taxa)
    style_map = get_combinatorial_taxon_style(active_taxa)
    has_cauline_any = bool((df["zone_clean"] == "cauline").any())

    # -------------------------------------------------------------------------
    # 1. Export Native LeafMachine3 Asset Cache (`ect_overlays/` & `lm3_assets/`)
    # -------------------------------------------------------------------------
    clean_for_ect = df[df["exemplar_ok"]].copy() if df["exemplar_ok"].any() else df
    try:
        export_taxon_ect_exemplars(extracted_df=clean_for_ect, figures_dir=out_dir, zone="basal", max_taxa=None)
        if has_cauline_any:
            export_taxon_ect_exemplars(extracted_df=clean_for_ect, figures_dir=out_dir, zone="cauline", max_taxa=None)
    except Exception as exc_ect:
        logger.warning(f"Could not export taxon LM3/ECT exemplars: {exc_ect}")

    # -------------------------------------------------------------------------
    # 2. Plate 02A: `fig02a_lm3_organ_and_ect_exemplars.png`
    #    Dynamic rows = ALL active taxa/ranks (zero [:8] truncation!)
    # -------------------------------------------------------------------------
    n_rows_a = max(1, len(active_taxa))
    n_cols_a = 5 if has_cauline_any else 4
    fig_w_a = 14.6 if has_cauline_any else 12.4
    fig_h_a = max(3.4, min(95.0, 1.72 * n_rows_a + 0.85))

    fig_a = plt.figure(figsize=(fig_w_a, fig_h_a), dpi=300, layout="constrained")
    gs_a = GridSpec(n_rows_a, n_cols_a, figure=fig_a, wspace=0.08, hspace=0.18)

    col_headers = (
        [
            "20th %ile (Compact)",
            "50th %ile (Medoid)",
            "80th %ile (Dissected)",
            "LM3 Radial ECT",
        ]
        if not has_cauline_any
        else [
            "Basal Medoid (LM3 Cutout)",
            "Basal 80th %ile (Dissected)",
            "Basal Radial ECT (LM3)",
            "Cauline Medoid (LM3)",
            "Cauline Radial ECT (LM3)",
        ]
    )

    def _pick_distinct_quantile_rows(cand_df: pd.DataFrame, quantiles: Sequence[float], window: float = 0.10) -> List[Optional[pd.Series]]:
        if len(cand_df) == 0:
            return [None] * len(quantiles)
        ranked = cand_df.assign(
            _ar=cand_df["ar_clean"].fillna(cand_df["ar_clean"].median() if cand_df["ar_clean"].notna().any() else 1.8),
            _sol=cand_df["solidity_clean"].fillna(0.85),
        ).sort_values(by=["_ar", "_sol"], ascending=[True, False])
        n_r = len(ranked)
        q_pos = [p / max(1, n_r - 1) if n_r > 1 else 0.5 for p in range(n_r)]
        meta_scores: List[Tuple[int, int, int, float]] = []
        for p in range(n_r):
            r_p = ranked.iloc[p]
            is_t = 1 if _parse_bool_val(r_p.get("is_type_specimen", False)) else 0
            t_str = str(r_p.get("determiner_tier", "")).strip()
            t_rk = 3 if ("Tier_1" in t_str or "Gold" in t_str) else (
                2 if ("Tier_2" in t_str or "Silver" in t_str) else (1 if ("Tier_3" in t_str or "Bronze" in t_str) else 0)
            )
            is_c = 1 if (
                _parse_bool_val(r_p.get("is_calibrated_mm", False))
                and str(r_p.get("scale_source", "")).strip().lower() not in {"mp_estimated", "predicted_from_megapixels", "uncalibrated"}
            ) else 0
            lqw_v = pd.to_numeric(r_p.get("leaf_quality_weight", r_p.get("archetype_score", 0.5)), errors="coerce")
            lqw_f = round(float(lqw_v), 4) if pd.notna(lqw_v) else 0.5
            meta_scores.append((is_t, t_rk, is_c, lqw_f))

        used_idx: Set[Any] = set()
        picked: List[Optional[pd.Series]] = []
        for q in quantiles:
            avail = [p for p in range(n_r) if ranked.index[p] not in used_idx]
            if not avail:
                picked.append(None)
                continue
            win = [p for p in avail if abs(q_pos[p] - q) <= window]
            if win:
                chosen_pos = max(win, key=lambda p: (*meta_scores[p], -abs(q_pos[p] - q)))
            else:
                chosen_pos = min(
                    avail,
                    key=lambda p: (
                        round(abs(q_pos[p] - q), 6),
                        -meta_scores[p][0],
                        -meta_scores[p][1],
                        -meta_scores[p][2],
                        -meta_scores[p][3],
                    ),
                )
            used_idx.add(ranked.index[chosen_pos])
            picked.append(ranked.iloc[chosen_pos])
        return picked

    for r_idx, sp_name in enumerate(active_taxa):
        sp_all = df[df["species_clean"] == sp_name]
        sp_clean = sp_all[sp_all["exemplar_ok"]]
        if len(sp_clean) == 0:
            # Warn if assets were missing on disk
            n_missing = sum(
                1 for _, r_chk in sp_all.iterrows() if _resolve_asset_path(r_chk.get("mask_path", ""), out_dir) is None
            )
            if n_missing > 0:
                logger.warning(
                    f"fig02a: {n_missing} / {len(sp_all)} {sp_name} extracted leaves have no resolvable mask/contour asset on disk."
                )
        col_hex, _ = style_map[sp_name]

        basal_pool = sp_clean[sp_clean["zone_clean"] == "basal"]
        if len(basal_pool) == 0 and len(sp_all[sp_all["zone_clean"] == "basal"]) > 0:
            # Only fall back if masks were synthetic/empty path (never admit flagged outliers/clipped masks)
            sp_basal_only = sp_all[sp_all["zone_clean"] == "basal"]
            empty_mask_flags = pd.Series(
                [not bool(str(pv or "").strip()) for pv in sp_basal_only.get("mask_path", pd.Series("", index=sp_basal_only.index))],
                index=sp_basal_only.index,
                dtype=bool,
            )
            basal_pool = sp_basal_only[(~sp_basal_only["cat_clean"].isin(misid_cats | outlier_cats)) & empty_mask_flags]

        cauline_pool = sp_clean[sp_clean["zone_clean"] == "cauline"]

        if not has_cauline_any:
            q_rows = _pick_distinct_quantile_rows(basal_pool, [0.24, 0.50, 0.80])
            med_for_ect = q_rows[1] if q_rows[1] is not None else (q_rows[0] if q_rows[0] is not None else q_rows[2])
            tile_specs = [
                ("leaf", q_rows[0], False, "20th %ile (Compact)"),
                ("leaf", q_rows[1], True, "50th %ile (Medoid)"),
                ("leaf", q_rows[2], False, "80th %ile (Dissected)"),
                ("ect", med_for_ect, False, "LM3 Radial ECT"),
            ]
        else:
            b_rows = _pick_distinct_quantile_rows(basal_pool, [0.50, 0.80])
            c_rows = _pick_distinct_quantile_rows(cauline_pool, [0.50])
            b_med = b_rows[0]
            c_med = c_rows[0]
            tile_specs = [
                ("leaf", b_med, True, "Basal Medoid"),
                ("leaf", b_rows[1], False, "Basal 80th %ile"),
                ("ect", b_med, False, "Basal Radial ECT"),
                ("leaf", c_med, True, "Cauline Medoid"),
                ("ect", c_med, False, "Cauline Radial ECT"),
            ]

        for c_idx, (kind, row_sel, prefer_rgb, sub_tag) in enumerate(tile_specs):
            ax = fig_a.add_subplot(gs_a[r_idx, c_idx])
            ax.set_xticks([])
            ax.set_yticks([])
            for sp_edge in ax.spines.values():
                sp_edge.set_edgecolor("#e5e7eb")
                sp_edge.set_linewidth(0.8)

            if c_idx == 0:
                wrapped_sp = format_icn_taxon_mathtext(sp_name, abbreviate_genus=(len(sp_name) > 20))
                ax.set_ylabel(
                    f"{wrapped_sp}\n(n = {len(sp_all)} leaves)",
                    fontsize=8.4,
                    color="#111827",
                    rotation=0,
                    ha="right",
                    va="center",
                    labelpad=12,
                )

            top_tag = col_headers[c_idx] if r_idx == 0 else sub_tag
            if row_sel is None:
                ax.set_facecolor("#f9fafb")
                ax.text(0.5, 0.5, "No clean organ\nin cohort", ha="center", va="center", fontsize=7.5, color="#9ca3af", transform=ax.transAxes)
                ax.set_title(f"{top_tag}\n—", fontsize=7.4, color="#6b7280", pad=3)
                continue

            cat_lbl = str(row_sel.get("cat_clean", row_sel.get("catalogNumber", "")))
            lid_lbl = int(float(row_sel.get("leaf_id", 1))) if pd.notna(row_sel.get("leaf_id", 1)) else 1
            mask_u8, rgb_patch = _load_mask_or_contour_u8(row_sel, out_dir)
            if mask_u8 is None or int(np.count_nonzero(mask_u8)) == 0:
                ax.set_facecolor("#f9fafb")
                ax.text(0.5, 0.5, f"{cat_lbl}\n(mask unavailable)", ha="center", va="center", fontsize=7.5, color="#9ca3af", transform=ax.transAxes)
                ax.set_title(f"{top_tag}\n{cat_lbl} (L{lid_lbl})", fontsize=7.4, color="#374151", pad=3)
                continue

            if kind == "ect":
                ect_card = _render_lm3_ect_rgb(mask_u8, n_dirs=160, size=320)
                if ect_card is not None:
                    ax.imshow(ect_card)
                else:
                    ax.text(0.5, 0.5, "ECT unavailable", ha="center", va="center", fontsize=7.5, transform=ax.transAxes)
            else:
                px_mm = pd.to_numeric(row_sel.get("pixels_per_mm", np.nan), errors="coerce")
                is_cal = _parse_bool_val(row_sel.get("is_calibrated_mm", False)) and str(
                    row_sel.get("scale_source", "")
                ).strip().lower() not in {"mp_estimated", "predicted_from_megapixels", "uncalibrated"}
                leaf_card = _render_lm3_silhouette_card(
                    mask_u8,
                    rgb_patch,
                    accent_hex=col_hex,
                    pixels_per_mm=float(px_mm) if pd.notna(px_mm) else None,
                    is_calibrated=is_cal,
                    card_size=320,
                    prefer_rgb=prefer_rgb,
                )
                ax.imshow(leaf_card)

            ar_v = pd.to_numeric(row_sel.get("ar_clean", np.nan), errors="coerce")
            sol_v = pd.to_numeric(row_sel.get("solidity_clean", np.nan), errors="coerce")
            stat_s = f"L/W={float(ar_v):.2f}, Sol={float(sol_v):.2f}" if (pd.notna(ar_v) and pd.notna(sol_v)) else f"L{lid_lbl}"
            ax.set_title(f"{top_tag}\n{cat_lbl} ({stat_s})", fontsize=7.3, color="#1f2937", pad=3)

    fig_a.suptitle(
        "Step 02 — Native LeafMachine3 Organ Silhouettes, Pose Landmarks & Radial ECT Exemplars",
        fontsize=11.5,
        fontweight="bold",
    )
    fig02a_png = _save_fig_fn(fig_a, out_dir / "fig02a_lm3_organ_and_ect_exemplars", dpi=300)
    plt.close(fig_a)
    saved_paths["fig02a_png"] = fig02a_png[1] if isinstance(fig02a_png, tuple) else Path(fig02a_png)

    # -------------------------------------------------------------------------
    # 3. Plate 02B: `fig02b_segmentation_qc_and_phenophase.png`
    #    Panel A (A1 Basal / A2 Cauline hollow scatter), Panel B (Scale & QC),
    #    Panel C (Involucral AR_inv hollow-jitter boxplots), Panel D (Phenophases)
    # -------------------------------------------------------------------------
    extra_h = max(0.0, 0.24 * max(0, len(active_taxa) - 8))
    if has_cauline_any:
        fig_h_b = 13.2 + 1.5 * extra_h
        fig_b = plt.figure(figsize=(14.4, fig_h_b), dpi=300, layout="constrained")
        gs_b = GridSpec(3, 2, figure=fig_b, height_ratios=[1.0, 1.0 + extra_h / 5.5, 0.95 + extra_h / 5.5])
        ax_b_basal = fig_b.add_subplot(gs_b[0, 0])
        ax_b_caul = fig_b.add_subplot(gs_b[0, 1])
        _draw_solidity_vs_ar_zone_panel(
            ax_b_basal,
            df[df["zone_clean"] == "basal"],
            active_taxa,
            style_map,
            "A1. Basal Rosette Homology Space (Hollow Markers)",
            show_legend=True,
        )
        _draw_solidity_vs_ar_zone_panel(
            ax_b_caul,
            df[df["zone_clean"] == "cauline"],
            active_taxa,
            style_map,
            "A2. Cauline Stem Homology Space (Hollow Markers)",
            show_legend=False,
        )
        ax_qc = fig_b.add_subplot(gs_b[1, 0])
        ax_ar_inv = fig_b.add_subplot(gs_b[1, 1])
        ax_pheno = fig_b.add_subplot(gs_b[2, 0])
        ax_het = fig_b.add_subplot(gs_b[2, 1])
    else:
        fig_h_b = 9.6 + extra_h
        fig_b = plt.figure(figsize=(14.2, fig_h_b), dpi=300, layout="constrained")
        gs_b = GridSpec(2, 2, figure=fig_b, height_ratios=[1.0, 1.0 + extra_h / 5.2])
        ax_b_topL = fig_b.add_subplot(gs_b[0, 0])
        _draw_solidity_vs_ar_zone_panel(
            ax_b_topL,
            df,
            active_taxa,
            style_map,
            "A. Basal Rosette Solidity vs. Blade Aspect Ratio (L / W)",
            show_legend=True,
        )
        ax_qc = fig_b.add_subplot(gs_b[0, 1])
        ax_ar_inv = fig_b.add_subplot(gs_b[1, 0])
        ax_pheno = fig_b.add_subplot(gs_b[1, 1])
        ax_het = None
    y_pos = np.arange(len(active_taxa))
    pct_ruler, pct_mp, pct_uncal, qc_pass_pct, counts_lbl = [], [], [], [], []
    for tx in active_taxa:
        sub_tx = df[df["species_clean"] == tx]
        n_t = max(1, len(sub_tx))
        sc_raw = sub_tx.get("scale_source", pd.Series("", index=sub_tx.index)).astype(str).str.strip().str.lower()
        cal_b = (
            sub_tx["is_calibrated_mm"].apply(_parse_bool_val)
            if "is_calibrated_mm" in sub_tx.columns
            else pd.Series(False, index=sub_tx.index)
        )
        is_mp = sc_raw.isin({"mp_estimated", "predicted_from_megapixels"})
        is_rul = (cal_b & (~is_mp)) | sc_raw.isin({"measured_from_ruler", "measured_from_fieldprism", "ruler_calibrated"})
        is_unc = ~(is_rul | is_mp)
        pct_ruler.append(100.0 * float(is_rul.sum()) / n_t)
        pct_mp.append(100.0 * float(is_mp.sum()) / n_t)
        pct_uncal.append(100.0 * float(is_unc.sum()) / n_t)
        qc_pass_pct.append(100.0 * float(sub_tx["exemplar_ok"].sum()) / n_t)
        counts_lbl.append(f"{format_icn_taxon_mathtext(tx, abbreviate_genus=(len(active_taxa) > 4))} (n={len(sub_tx)})")

    b1 = np.array(pct_ruler)
    b2 = np.array(pct_mp)
    b3 = np.array(pct_uncal)
    ax_qc.barh(y_pos, b1, height=0.62, color="#0072B2", alpha=0.85, label="Ruler / FieldPrism (mm)")
    ax_qc.barh(y_pos, b2, left=b1, height=0.62, color="#E69F00", alpha=0.85, label="Megapixel Prior")
    ax_qc.barh(y_pos, b3, left=b1 + b2, height=0.62, color="#d1d5db", alpha=0.90, label="Uncalibrated (px)")
    ax_qc.scatter(
        qc_pass_pct,
        y_pos,
        s=42,
        marker="D",
        facecolors="none",
        edgecolors="#111827",
        linewidths=1.2,
        zorder=5,
        label="Clean EFA/LM3 QC %",
    )
    ax_qc.set_yticks(y_pos)
    ax_qc.set_yticklabels(counts_lbl, fontsize=7.8)
    ax_qc.invert_yaxis()
    ax_qc.set_xlim(0, 108)
    ax_qc.set_xlabel("Percentage of Extracted Leaves (%)", fontsize=9)
    n_fail_sheets = len(failed_df) if failed_df is not None else 0
    ax_qc.set_title(f"B. Ruler Calibration & Clean LM3 Gate Yield (Failed: {n_fail_sheets})", fontsize=10, fontweight="bold", loc="left", pad=18)
    ax_qc.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=3, fontsize=7.0, frameon=False)
    ax_qc.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)

    # Panels C & D allocated in clean half-width columns above

    # Panel C: Involucral Aspect Ratio (AR_inv) hollow-jitter boxplots
    rep_work = pd.DataFrame()
    if repro_df is not None and len(repro_df) > 0 and "AR_inv" in repro_df.columns:
        rep_work = repro_df.copy()
        if "species_clean" not in rep_work.columns and "catalogNumber" in rep_work.columns:
            cat_to_sp = dict(zip(df["cat_clean"], df["species_clean"]))
            rep_work["species_clean"] = rep_work["catalogNumber"].astype(str).str.strip().map(cat_to_sp)
        rep_work["AR_inv_num"] = pd.to_numeric(rep_work["AR_inv"], errors="coerce")
        rep_work = rep_work[rep_work["AR_inv_num"].notna() & rep_work["species_clean"].isin(active_taxa)]

    if len(rep_work) > 0:
        taxa_with_rep = [tx for tx in active_taxa if (rep_work["species_clean"] == tx).any()]
        box_data = [rep_work.loc[rep_work["species_clean"] == tx, "AR_inv_num"].to_numpy() for tx in taxa_with_rep]
        pos_r = np.arange(len(taxa_with_rep))
        bp = ax_ar_inv.boxplot(
            box_data,
            positions=pos_r,
            orientation="horizontal",
            widths=0.52,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color="#111827", linewidth=1.4),
        )
        for patch_b in bp["boxes"]:
            patch_b.set_facecolor("none")
            patch_b.set_edgecolor("#4b5563")
            patch_b.set_linewidth(0.95)
        rng_j = np.random.default_rng(42)
        for idx_r, tx in enumerate(taxa_with_rep):
            vals_r = box_data[idx_r]
            col_r, mk_r = style_map[tx]
            jit = rng_j.uniform(-0.16, 0.16, size=len(vals_r))
            ax_ar_inv.scatter(
                vals_r,
                np.full(len(vals_r), idx_r, dtype=float) + jit,
                s=28,
                marker=mk_r,
                facecolors="none",
                edgecolors=col_r,
                linewidths=0.95,
                alpha=0.82,
                zorder=4,
            )
        ax_ar_inv.axvline(1.0, color="#9ca3af", linestyle=":", linewidth=0.8)
        ax_ar_inv.set_yticks(pos_r)
        ax_ar_inv.set_yticklabels(
            [f"{format_icn_taxon_mathtext(tx, abbreviate_genus=True)} (n={len(box_data[i])})" for i, tx in enumerate(taxa_with_rep)],
            fontsize=7.6,
        )
        ax_ar_inv.invert_yaxis()
        ax_ar_inv.set_xlabel("Involucral Aspect Ratio (Height / Width)", fontsize=8.8)
        ax_ar_inv.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)
    else:
        ax_ar_inv.text(0.5, 0.5, "No measurable capitula\nin cohort", ha="center", va="center", fontsize=8.5, color="#6b7280", transform=ax_ar_inv.transAxes)
        ax_ar_inv.set_xticks([])
        ax_ar_inv.set_yticks([])
    ax_ar_inv.set_title("C. Capitulum Involucre Aspect Ratio", fontsize=10, fontweight="bold", loc="left")

    # Panel D: Determinate Cyme Phenophase Proportions
    ph_work = pd.DataFrame()
    if pheno_df is not None and len(pheno_df) > 0:
        ph_work = pheno_df.copy()
        st_col = next(
            (c for c in ("voucher_phenological_state", "phenological_state", "pheno_state") if c in ph_work.columns),
            None,
        )
        if st_col is not None:
            if "species_clean" not in ph_work.columns and "catalogNumber" in ph_work.columns:
                cat_to_sp = dict(zip(df["cat_clean"], df["species_clean"]))
                ph_work["species_clean"] = ph_work["catalogNumber"].astype(str).str.strip().map(cat_to_sp)
            ph_work["state_norm"] = ph_work[st_col].astype(str).str.strip().str.lower()
            ph_work = ph_work[ph_work["species_clean"].isin(active_taxa)]
        else:
            ph_work = pd.DataFrame()

    if len(ph_work) > 0:
        taxa_ph = [tx for tx in active_taxa if (ph_work["species_clean"] == tx).any()]
        pos_p = np.arange(len(taxa_ph))
        cats_ph = [
            ("Anthesis", lambda s: s.str.contains("anthesis|flower|bloom", regex=True), "#E69F00"),
            ("Bud", lambda s: s.str.contains("bud", regex=True) & ~s.str.contains("anthesis", regex=True), "#56B4E9"),
            ("Fruiting", lambda s: s.str.contains("fruit|pappus", regex=True) & ~s.str.contains("anthesis", regex=True), "#009E73"),
            ("Sterile / Senescent", lambda s: ~(s.str.contains("anthesis|flower|bloom|bud|fruit|pappus", regex=True)), "#9ca3af"),
        ]
        left_acc = np.zeros(len(taxa_ph), dtype=float)
        for lbl_p, pred_fn, col_p in cats_ph:
            vals_p = []
            for tx in taxa_ph:
                s_tx = ph_work.loc[ph_work["species_clean"] == tx, "state_norm"]
                m_hit = pred_fn(s_tx)
                vals_p.append(100.0 * float(m_hit.sum()) / max(1, len(s_tx)))
            arr_p = np.array(vals_p, dtype=float)
            ax_pheno.barh(pos_p, arr_p, left=left_acc, height=0.60, color=col_p, alpha=0.88, label=lbl_p)
            left_acc += arr_p
        ax_pheno.set_yticks(pos_p)
        ax_pheno.set_yticklabels(
            [format_icn_taxon_mathtext(tx, abbreviate_genus=True) for tx in taxa_ph],
            fontsize=7.6,
        )
        ax_pheno.invert_yaxis()
        ax_pheno.set_xlim(0, 100)
        ax_pheno.set_xlabel("Voucher Proportion (%)", fontsize=8.8)
        ax_pheno.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=2, fontsize=6.9, frameon=False)
        ax_pheno.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)
    else:
        ax_pheno.text(0.5, 0.5, "No capitulum phenology\nrecords in cohort", ha="center", va="center", fontsize=8.5, color="#6b7280", transform=ax_pheno.transAxes)
        ax_pheno.set_xticks([])
        ax_pheno.set_yticks([])
    ax_pheno.set_title("D. Cyme Phenophase Stratification", fontsize=10, fontweight="bold", loc="left", pad=22)

    if ax_het is not None:
        # Panel E: Basal vs. Cauline Organ Counts & Heteroblastic Solidity Shift per Taxon
        b_cnts = [int(((df["species_clean"] == tx) & (df["zone_clean"] == "basal")).sum()) for tx in active_taxa]
        c_cnts = [int(((df["species_clean"] == tx) & (df["zone_clean"] == "cauline")).sum()) for tx in active_taxa]
        pos_h = np.arange(len(active_taxa))
        ax_het.barh(pos_h - 0.17, b_cnts, height=0.32, color="#0072B2", alpha=0.86, label="Basal Rosette Leaves")
        ax_het.barh(pos_h + 0.17, c_cnts, height=0.32, color="#009E73", alpha=0.86, label="Cauline Stem Leaves")
        ax_het.set_yticks(pos_h)
        ax_het.set_yticklabels([format_icn_taxon_mathtext(tx, abbreviate_genus=True) for tx in active_taxa], fontsize=7.6)
        ax_het.invert_yaxis()
        ax_het.set_xlabel("Extracted Organ Count (n)", fontsize=8.8)
        ax_het.set_title("E. Basal Rosette vs. Cauline Stem Yield by Taxon", fontsize=10, fontweight="bold", loc="left")
        ax_het.legend(loc="lower right", fontsize=7.4, frameon=True, facecolor="white", framealpha=0.9)
        ax_het.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.55)

    fig_b.suptitle(
        "Step 02 — Segmentation Quality Control, Ruler Calibration & Reproductive Traits",
        fontsize=11.5,
        fontweight="bold",
    )
    fig02b_png = _save_fig_fn(fig_b, out_dir / "fig02b_segmentation_qc_and_phenophase", dpi=300)
    plt.close(fig_b)
    saved_paths["fig02b_png"] = fig02b_png[1] if isinstance(fig02b_png, tuple) else Path(fig02b_png)

    # -------------------------------------------------------------------------
    # 4. Capped Whole-Sheet LeafMachine3 Diagnostic Overlays (`lm3_sheet_overlays/`)
    # -------------------------------------------------------------------------
    try:
        overlays_written = export_lm3_sheet_diagnostic_overlays(
            extracted_df=df,
            repro_df=repro_df,
            pheno_df=pheno_df,
            figures_dir=out_dir,
            vouchers_df=vouchers_df,
            max_sheets=24,
        )
        if overlays_written:
            saved_paths["sheet_overlays_dir"] = out_dir / "lm3_sheet_overlays"
    except Exception as exc_sh:
        logger.warning(f"Whole-sheet LeafMachine3 diagnostic overlay export note: {exc_sh}")

    logger.info(f"Step 02 publication figures exported (300-DPI PNG) to: {out_dir}")
    return saved_paths


def regenerate_step02_figures_from_tables(
    output_dir: Union[str, Path],
    vouchers_csv: Optional[Union[str, Path]] = None,
) -> Dict[str, Path]:
    """Reload Step 02 + downstream QC tables and regenerate 300-DPI PNG Step 02 plates."""
    out_data_dir = Path(output_dir)
    tables_dir = out_data_dir / "tables" if (out_data_dir / "tables").exists() else out_data_dir
    figures_dir = out_data_dir.parent / "outputs" / "figures"

    def _safe_read_csv(p: Path) -> pd.DataFrame:
        if p.exists() and p.stat().st_size > 0:
            try:
                return pd.read_csv(p, low_memory=False)
            except Exception as e_csv:
                logger.warning(f"Failed to read {p}: {e_csv}")
        return pd.DataFrame()

    extracted_df = _safe_read_csv(tables_dir / "extracted_leaf_manifest.csv")
    failed_df = _safe_read_csv(tables_dir / "failed_qc_vouchers.csv")
    repro_df = _safe_read_csv(tables_dir / "voucher_reproductive_metrics.csv")
    pheno_df = _safe_read_csv(tables_dir / "voucher_phenological_states.csv")
    vouchers_df: Optional[pd.DataFrame] = None

    if vouchers_csv is not None and Path(vouchers_csv).exists() and len(extracted_df) > 0:
        try:
            vouchers_df = pd.read_csv(vouchers_csv, low_memory=False)
            if "catalogNumber" in vouchers_df.columns and "catalogNumber" in extracted_df.columns:
                active_vcats = set(vouchers_df["catalogNumber"].astype(str).str.strip())
                extracted_df = extracted_df[
                    extracted_df["catalogNumber"].astype(str).str.strip().isin(active_vcats)
                ].copy()
        except Exception:
            vouchers_df = None

    return generate_step02_publication_figures(
        extracted_df=extracted_df,
        failed_df=failed_df,
        repro_df=repro_df,
        pheno_df=pheno_df,
        figures_dir=figures_dir,
        vouchers_df=vouchers_df,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Step 02 300-DPI PNG Figures and LeafMachine3 Exemplar Exporter")
    parser.add_argument("--export-ect-exemplars", action="store_true", help="Export authentic LeafMachine3 radial ECT & silhouette overlays")
    parser.add_argument("--output-dir", type=str, default="outputs/figures", help="Output directory for figures/overlays")
    parser.add_argument("--zone", type=str, default="all", choices=["basal", "cauline", "all"], help="Leaf zone (basal, cauline, or all)")
    parser.add_argument("--vouchers-csv", type=str, default="", help="Path to curated_vouchers.csv")
    parser.add_argument("--manifest-csv", type=str, default="", help="Path to extracted_leaf_manifest.csv")
    parser.add_argument("--regenerate-figures-only", action="store_true", help="Regenerate Step 02 publication figures")
    args = parser.parse_args()

    if args.export_ect_exemplars:
        out_fdir = Path(args.output_dir)
        v_csv = args.vouchers_csv or None
        ext_in = pd.read_csv(args.manifest_csv, low_memory=False) if (args.manifest_csv and Path(args.manifest_csv).exists()) else None
        zones = ["basal", "cauline"] if args.zone == "all" else [args.zone]
        for z in zones:
            export_taxon_ect_exemplars(extracted_df=ext_in, figures_dir=out_fdir, zone=z, max_taxa=None, vouchers_csv=v_csv)
    elif args.regenerate_figures_only:
        regenerate_step02_figures_from_tables(output_dir=args.output_dir, vouchers_csv=args.vouchers_csv or None)
