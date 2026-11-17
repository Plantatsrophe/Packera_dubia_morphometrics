#!/usr/bin/env python3
"""
Module: scripts/data_prep/step06_ploidy_hdplot_utils.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  nQuire Gaussian-plus-uniform-noise EM cytotype inference (`2x`/`3x`/`4x`/`6x`/
  `8x` prior support and signed BIC scores), HDplot paralog filtering with
  Compositae-ParaLoss-1272 weights, and taxon-stratified Empirical Bayes polyRAD
  dosage estimation.
"""

from __future__ import annotations

import csv, logging, math, re
from math import comb
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from scripts.data_prep.step06_voucher_linker_utils import _parse_int_ploidy

logger = logging.getLogger("plantmatics.genomics.ingest")

_CYTO_TO_PLOIDY = {"2x": 2, "3x": 3, "4x": 4, "6x": 6, "8x": 8}


def _em_nquire_component(
    ratios: np.ndarray, centers: Sequence[float], *, update_means: bool = False, max_iter: int = 25,
) -> Tuple[float, float, float, np.ndarray, np.ndarray]:
    """Fit Gaussian + U(0.10, 0.90) noise EM mixture for allele balance ratios."""
    n, c_num = int(ratios.size), len(centers)
    mu = np.asarray(centers, dtype=np.float64).copy()
    pi_0 = 0.05
    pi_c = np.full(c_num, (1.0 - pi_0) / max(1, c_num), dtype=np.float64)
    sigma = 0.07
    log_u_dens = math.log(1.25)

    def _eval_log_joint(cur_pi0: float, cur_pic: np.ndarray, cur_mu: np.ndarray, cur_sig: float) -> Tuple[np.ndarray, np.ndarray, float]:
        log_u_col = np.full((n, 1), math.log(max(1e-9, cur_pi0)) + log_u_dens, dtype=np.float64)
        diffs = ratios[:, None] - cur_mu[None, :]
        inv_var = 1.0 / (2.0 * cur_sig * cur_sig)
        log_norm = -0.5 * math.log(2.0 * math.pi * cur_sig * cur_sig) - (diffs * diffs) * inv_var + np.log(np.maximum(1e-9, cur_pic))[None, :]
        log_joint = np.hstack([log_u_col, log_norm])
        m = np.max(log_joint, axis=1)
        lse = m + np.log(np.sum(np.exp(log_joint - m[:, None]), axis=1))
        return log_joint, lse, float(np.sum(lse))

    prev_logl = -math.inf
    for it in range(max_iter):
        log_joint, lse, logl = _eval_log_joint(pi_0, pi_c, mu, sigma)
        if it > 0 and abs(logl - prev_logl) < 1e-6:
            break
        prev_logl = logl
        resp = np.exp(log_joint - lse[:, None])
        r0, rc = resp[:, 0], resp[:, 1:]
        nkc = np.sum(rc, axis=0)
        pi_0 = float(min(0.50, max(0.01, float(np.sum(r0)) / max(1, n))))
        w_raw = np.maximum(1e-4, nkc / max(1e-9, float(np.sum(nkc))))
        pi_c = (1.0 - pi_0) * (w_raw / max(1e-12, float(np.sum(w_raw))))
        if update_means:
            mu = np.clip(np.sum(rc * ratios[:, None], axis=0) / np.maximum(1e-9, nkc), 0.12, 0.88)
        diffs = ratios[:, None] - mu[None, :]
        sigma = float(min(0.16, max(0.03, math.sqrt(max(1e-12, float(np.sum(rc * (diffs * diffs))) / max(1e-9, float(np.sum(rc))))))))

    _, _, final_logl = _eval_log_joint(pi_0, pi_c, mu, sigma)
    return final_logl, sigma, pi_0, pi_c, mu


def _fit_nquire_em_models(ratios: np.ndarray, max_iter: int = 25) -> Dict[str, Any]:
    """Fit Weiß et al. (2018) EM Gaussian-plus-uniform-noise ploidy models (`2x/3x/4x/6x`)."""
    n = max(2, int(ratios.size))
    log_n = math.log(n)
    specs: List[Tuple[str, Sequence[float], int, bool]] = [
        ("2x", [0.50], 2, False),
        ("3x", [1.0 / 3.0, 2.0 / 3.0], 3, False),
        ("4x", [0.25, 0.50, 0.75], 4, False),
        ("6x", [1.0 / 6.0, 2.0 / 6.0, 0.50, 4.0 / 6.0, 5.0 / 6.0], 6, False),
    ]
    logl_map, bic_map, sigma_map, pi0_map = {}, {}, {}, {}
    for label, centers, k_params, upd_means in specs:
        logl, sig, p0, _, _ = _em_nquire_component(ratios, centers, update_means=upd_means, max_iter=max_iter)
        logl_map[label] = float(logl)
        bic_map[label] = float(-2.0 * logl + k_params * log_n)
        sigma_map[label] = float(sig)
        pi0_map[label] = float(p0)

    best_label = min(("2x", "3x", "4x", "6x"), key=lambda lbl: bic_map[lbl])
    delta_bic_vs_2x = float(max(0.0, bic_map["2x"] - bic_map[best_label]))
    bic_gap_2x_vs_poly = float(min(bic_map["3x"], bic_map["4x"], bic_map["6x"]) - bic_map["2x"])
    return {
        "best_label": best_label, "logL_map": logl_map, "bic_map": bic_map,
        "sigma_map": sigma_map, "pi0_map": pi0_map,
        "delta_bic_vs_2x": delta_bic_vs_2x, "bic_gap_2x_vs_poly": bic_gap_2x_vs_poly,
    }


def _infer_sample_cytotype(
    ref_vec: np.ndarray, alt_vec: np.ndarray, dosage_vec: np.ndarray,
    user_ploidy: Optional[int], species_modal_ploidy: Optional[int] = None, is_multi_cytotype: bool = False,
) -> Tuple[str, Optional[float]]:
    """Infer cytotype (`2x`/`3x`/`4x`/`6x`/`8x`/`mixed/uncertain`) and signed nQuire EM BIC score."""
    valid = np.isfinite(ref_vec) & np.isfinite(alt_vec) & np.isfinite(dosage_vec)
    dp = ref_vec + alt_vec
    het_mask = valid & (dp >= 10.0) & (dosage_vec > 0.05) & (dosage_vec < 0.95)
    ratios = (alt_vec[het_mask] / np.maximum(1.0, dp[het_mask])).astype(np.float64)
    ratios = ratios[(ratios >= 0.10) & (ratios <= 0.90)]
    if ratios.size < 15:
        if user_ploidy in (2, 3, 4, 6, 8):
            return f"{int(user_ploidy)}x", 0.0
        if species_modal_ploidy in (2, 3, 4, 6, 8) and not is_multi_cytotype:
            return f"{int(species_modal_ploidy)}x", 0.0
        return "mixed/uncertain", None

    fit = _fit_nquire_em_models(ratios, max_iter=25)
    best_label = str(fit["best_label"])
    delta_bic_vs_2x = float(fit["delta_bic_vs_2x"])
    bic_gap_2x_vs_poly = float(fit["bic_gap_2x_vs_poly"])
    prior_p: Optional[int] = (
        int(user_ploidy) if user_ploidy in (2, 3, 4, 6, 8)
        else (int(species_modal_ploidy) if (species_modal_ploidy in (2, 3, 4, 6, 8) and not is_multi_cytotype) else None)
    )

    if best_label == "2x":
        if prior_p in (3, 4, 6, 8) and abs(bic_gap_2x_vs_poly) < 4.0:
            return f"{prior_p}x", 0.0
        return "2x", round(-bic_gap_2x_vs_poly, 4)

    strong_thresh = 6.0 if ratios.size < 40 else 4.0
    if delta_bic_vs_2x >= strong_thresh:
        return best_label, round(delta_bic_vs_2x, 4)
    if prior_p == 2:
        return "2x", 0.0
    if prior_p in (3, 4, 6, 8):
        return (best_label, round(delta_bic_vs_2x, 4)) if delta_bic_vs_2x >= 4.0 else (f"{prior_p}x", 0.0)
    return "mixed/uncertain", round(delta_bic_vs_2x, 4)


def _load_paraloss1272_catalog(genomics_dir: Optional[Path]) -> Dict[str, Tuple[int, float]]:
    """Load staged `public_packera_paraloss1272_loci.csv` locus weights when available."""
    if genomics_dir is None:
        return {}
    gdir, fname = Path(genomics_dir), "public_packera_paraloss1272_loci.csv"
    cat_path = next((p for p in (gdir / "public_catalogs" / fname, gdir / fname) if p.is_file() and p.stat().st_size > 0), next((p for p in sorted(gdir.rglob(fname)) if p.is_file() and p.stat().st_size > 0), None))
    if cat_path is None:
        return {}
    catalog: Dict[str, Tuple[int, float]] = {}
    try:
        with open(cat_path, encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh):
                if not (lid := str(row.get("locus_id") or "").strip()):
                    continue
                p_cls, raw_w = str(row.get("paralogy_class") or "").strip().lower(), str(row.get("locus_paralogy_weight") or "").strip()
                fb_wt = 1.5 if "single" in p_cls else (0.5 if any(k in p_cls for k in ("multi", "paralog")) else 1.0)
                try:
                    wt = float(raw_w) if raw_w else fb_wt
                except ValueError:
                    wt = fb_wt
                is_sc = 1 if ("single" in p_cls or wt > 1.05) else 0
                catalog[lid] = catalog[lid.lower()] = (is_sc, wt)
    except Exception as exc:
        logger.warning("Could not parse ParaLoss-1272 catalog %s: %s", cat_path, exc)
    return catalog


def _lookup_paraloss1272_locus(locus_id: str, catalog: Mapping[str, Tuple[int, float]]) -> Tuple[bool, int, float]:
    """Match locus ID or stripped SNP suffix against `public_packera_paraloss1272_loci.csv`."""
    if not (raw := str(locus_id or "").strip()) or not catalog:
        return False, 0, 1.0
    cands: List[str] = [raw] + [raw.split(sep, 1)[0] for sep in (":", ".") if sep in raw]
    if (no_pos := re.split(r"(?i)_pos(?:[\d_:-].*)?$", raw, maxsplit=1)[0]) and no_pos != raw:
        cands.append(no_pos)
    if base := re.split(r"(?i)(?::|_pos|\.)", raw, maxsplit=1)[0]:
        cands.append(base)
    for cand in cands:
        if (hit := catalog.get(cand) or catalog.get(cand.lower())) is not None:
            return True, int(hit[0]), float(hit[1])
    return False, 0, 1.0


def _evaluate_hdplot_loci(
    sub_pass: np.ndarray, sub_ref: np.ndarray, ref_d_ref: np.ndarray, alt_d_ref: np.ndarray,
    locus_ids: Sequence[str], h_max: float = 0.60, d_max: float = 4.0,
    paraloss_catalog: Optional[Mapping[str, Tuple[int, float]]] = None, genomics_dir: Optional[Path] = None,
) -> Tuple[np.ndarray, pd.DataFrame]:
    """Apply McKinney et al. (2017) HDplot paralog filter (`H <= 0.60`, `|D| <= 4.0`) without MAF pruning."""
    cat = paraloss_catalog if paraloss_catalog is not None else _load_paraloss1272_catalog(genomics_dir)
    n_loci = len(locus_ids)
    keep_mask = np.ones(n_loci, dtype=bool)
    diag_rows: List[Dict[str, Any]] = []
    has_ref_depths = ref_d_ref.size > 0 and alt_d_ref.size > 0

    for j in range(n_loci):
        col_p, col_r = sub_pass[:, j], sub_ref[:, j]
        fin_p, fin_mask_r = col_p[np.isfinite(col_p)], np.isfinite(col_r)
        n_called_r = int(np.sum(fin_mask_r))
        het_mask_r = fin_mask_r & (col_r > 0.05) & (col_r < 0.95)
        n_hets_r = int(np.sum(het_mask_r))
        h_obs = float(n_hets_r / float(n_called_r)) if n_called_r > 0 else 0.0
        na_sum = nb_sum = d_stat = 0.0
        if has_ref_depths and n_hets_r > 0:
            rr, aa = ref_d_ref[het_mask_r, j], alt_d_ref[het_mask_r, j]
            if np.any(ok_d := np.isfinite(rr) & np.isfinite(aa)):
                na_sum, nb_sum = float(np.sum(rr[ok_d])), float(np.sum(aa[ok_d]))
                if (n_tot := na_sum + nb_sum) >= 10.0:
                    d_stat = float((na_sum - 0.5 * n_tot) / np.sqrt(0.25 * n_tot))
        reason = "PASS"
        if fin_p.size < 2: keep_mask[j] = False; reason = "monomorphic_or_sparse"
        elif n_called_r >= 4 and h_obs > float(h_max): keep_mask[j] = False; reason = f"HDplot_H>{h_max:.2f}"
        elif (na_sum + nb_sum) >= 10.0 and abs(d_stat) > float(d_max): keep_mask[j] = False; reason = f"HDplot_|D|>{d_max:.1f}"
        elif float(np.std(fin_p)) <= 1e-8: keep_mask[j] = False; reason = "monomorphic_or_sparse"
        matched_pl, is_sc, p_wt = _lookup_paraloss1272_locus(locus_ids[j], cat)
        diag_rows.append({
            "locus_id": locus_ids[j], "n_called_ref": n_called_r, "n_hets_ref": n_hets_r,
            "H_obs": round(h_obs, 4) if n_called_r > 0 else np.nan,
            "het_ref_reads": round(na_sum, 1), "het_alt_reads": round(nb_sum, 1),
            "het_total_reads": round(na_sum + nb_sum, 1), "D_read_ratio": round(d_stat, 4),
            "hdplot_keep": 1 if keep_mask[j] else 0, "filter_reason": reason,
            "is_paraloss1272_single_copy": int(is_sc), "locus_paralogy_weight": round(float(p_wt if matched_pl else 1.0), 3),
        })
    return keep_mask, pd.DataFrame(diag_rows)


def _compute_eb_polyrad_dosage(
    dosage_sub: np.ndarray, ref_sub: np.ndarray, alt_sub: np.ndarray, cytotypes: Sequence[str],
    fallback_ploidies: Optional[Sequence[Optional[int]]] = None, eps: float = 0.005,
    taxa_labels: Optional[Sequence[Optional[str]]] = None,
) -> np.ndarray:
    """Estimate taxon-stratified Empirical Bayes discrete polyploid posterior dosage in `[0, 1]`."""
    out = dosage_sub.copy().astype(np.float32)
    if ref_sub.size == 0 or alt_sub.size == 0 or out.size == 0:
        return out

    p_hat = np.nanmean(dosage_sub, axis=0)
    p_global = np.where(np.isfinite(p_hat), np.clip(p_hat, 0.01, 0.99), 0.5).astype(np.float64)
    p_by_taxon: Dict[str, np.ndarray] = {}
    if taxa_labels is not None:
        uniq_tx = {str(lbl).strip() for lbl in taxa_labels if lbl is not None and str(lbl).strip()}
        for tx in uniq_tx:
            tx_rows = [idx for idx, lbl in enumerate(taxa_labels) if str(lbl or "").strip() == tx]
            if tx_rows:
                sub_tx = dosage_sub[tx_rows, :]
                n_tx = np.sum(np.isfinite(sub_tx), axis=0).astype(np.float64)
                sum_tx = np.nansum(sub_tx, axis=0).astype(np.float64)
                p_by_taxon[tx] = np.clip((sum_tx + 2.0 * p_global) / (n_tx + 2.0), 0.01, 0.99)

    for i, cy in enumerate(cytotypes):
        ploidy_p = _CYTO_TO_PLOIDY.get(str(cy or "").strip().lower())
        if ploidy_p is None and fallback_ploidies is not None and i < len(fallback_ploidies):
            if (fb := _parse_int_ploidy(fallback_ploidies[i])) in (2, 3, 4, 6, 8):
                ploidy_p = int(fb)
        if ploidy_p is None:
            ploidy_p = 2
        ks = np.arange(ploidy_p + 1, dtype=np.float64)
        ratios_k = ks / float(ploidy_p)
        ek = np.clip((1.0 - eps) * ratios_k + eps * (1.0 - ratios_k), 1e-6, 1.0 - 1e-6)
        log_comb = np.log(np.asarray([comb(ploidy_p, int(k)) for k in ks], dtype=np.float64))

        nr, na = ref_sub[i, :].astype(np.float64), alt_sub[i, :].astype(np.float64)
        called = np.isfinite(dosage_sub[i, :]) & np.isfinite(nr) & np.isfinite(na) & ((nr + na) > 0.0)
        if not np.any(called):
            continue
        tx_lbl = str(taxa_labels[i] or "").strip() if (taxa_labels is not None and i < len(taxa_labels)) else ""
        p_vec = p_by_taxon.get(tx_lbl, p_global)
        p_c = p_vec[called][:, None]
        prior = log_comb[None, :] + ks[None, :] * np.log(p_c) + (ploidy_p - ks[None, :]) * np.log(1.0 - p_c)
        ll = na[called][:, None] * np.log(ek[None, :]) + nr[called][:, None] * np.log(1.0 - ek[None, :])
        log_post = prior + ll
        log_post -= np.max(log_post, axis=1, keepdims=True)
        w = np.exp(log_post)
        w /= np.sum(w, axis=1, keepdims=True)
        out[i, called] = np.sum(w * ratios_k[None, :], axis=1).astype(np.float32)

    return np.clip(out, 0.0, 1.0)

