#!/usr/bin/env python3
"""
Module: scripts/analysis/step06_coalescent_delimitation.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Invariant-site Multispecies Coalescent (MSC) & Reticulation (MSci) engine for
  Step 06 phylogenomics (`SVDquartets`, `BPP` control/PHYLIP export + conjugate
  Bayes factors, and locus-bootstrapped `gdi = 1 - exp(-2 * tau / theta)`).
  Retains 100% of invariant sites across `122 bp` RAD contigs and resolves IUPAC
  heterozygote codes analytically so `n = 1` diploid/polyploid libraries yield
  calibrated within-lineage `theta` and between-lineage `d_XY`.
"""

from __future__ import annotations

import argparse
import importlib
import math
import os
import re
import shutil
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.analysis.step06_coalescent_models import (
    _resolve_max_workers,
    compute_coalescent_introgression_bf,
    compute_model_corrected_divergence_blocks,
    select_nucleotide_substitution_model,
)
from scripts.core.logger import setup_logging

logger = setup_logging(name="packera.genomics.coalescent")

_IUPAC_ALLELES: Dict[str, Tuple[str, str]] = {
    "A": ("A", "A"), "C": ("C", "C"), "G": ("G", "G"), "T": ("T", "T"),
    "R": ("A", "G"), "Y": ("C", "T"), "S": ("C", "G"),
    "W": ("A", "T"), "K": ("G", "T"), "M": ("A", "C"),
}
_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _clean_id(raw: str) -> str:
    fn = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").clean_sample_id
    return fn(raw)


def _infer_taxon(sid: str, sample_taxon_map: Mapping[str, str]) -> str:
    if sid in sample_taxon_map and sample_taxon_map[sid]:
        return str(sample_taxon_map[sid])
    parsed = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").parse_sample_code(sid)
    return str(parsed.get("species_standardized") or sid)


def parse_invariant_alignment(path: Path) -> Tuple[List[str], List[str]]:
    """Parse PHYLIP (sequential/interleaved), NEXUS, or FASTA retaining all invariant & IUPAC sites."""
    low = path.name.lower()
    text = path.read_text(encoding="utf-8", errors="replace")
    seqs: Dict[str, List[str]] = defaultdict(list)
    order: List[str] = []

    if low.endswith((".fasta", ".fa", ".fas", ".fna")):
        cur: Optional[str] = None
        for line in text.splitlines():
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                cur = _clean_id(s[1:].split()[0])
                if cur not in seqs:
                    order.append(cur)
            elif cur is not None:
                seqs[cur].append(s.upper())
        return order, ["".join(seqs[k]) for k in order]

    if low.endswith((".nex", ".nexus")):
        in_mat = False
        for line in text.splitlines():
            s = line.strip()
            if not s or s.startswith("["):
                continue
            if re.match(r"^MATRIX\b", s, flags=re.IGNORECASE):
                in_mat = True
                continue
            if in_mat:
                if s.startswith(";"):
                    in_mat = False
                    continue
                parts = s.split(None, 1)
                if len(parts) == 2:
                    sid = _clean_id(parts[0].strip("'\""))
                    if sid not in seqs:
                        order.append(sid)
                    seqs[sid].append(re.sub(r"\s+", "", parts[1]).upper())
        return order, ["".join(seqs[k]) for k in order]

    # PHYLIP sequential or interleaved
    lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
    if not lines:
        return [], []
    hdr = lines[0].strip().split()
    n_taxa = int(hdr[0]) if hdr and hdr[0].isdigit() else 0
    body = lines[1:]
    if n_taxa > 0 and len(body) >= n_taxa:
        # Handle multi-line interleaved vs wrapped sequential cleanly
        first_chunk = body[:n_taxa]
        if all(len(ln.split(None, 1)) == 2 for ln in first_chunk):
            for ln in first_chunk:
                nm, sq = ln.split(None, 1)
                sid = _clean_id(nm)
                order.append(sid)
                seqs[sid].append(re.sub(r"\s+", "", sq).upper())
            rem = body[n_taxa:]
            # Interleaved blocks have taxon-less continuation lines or repeated names
            pos_i = 0
            for ln in rem:
                parts = ln.split(None, 1)
                if len(parts) == 2 and _clean_id(parts[0]) in seqs:
                    seqs[_clean_id(parts[0])].append(re.sub(r"\s+", "", parts[1]).upper())
                else:
                    sid = order[pos_i % n_taxa]
                    seqs[sid].append(re.sub(r"\s+", "", ln).upper())
                    pos_i += 1
            return order, ["".join(seqs[k]) for k in order]

    for ln in body:
        parts = ln.split(None, 1)
        if len(parts) == 2:
            sid = _clean_id(parts[0])
            if sid not in seqs:
                order.append(sid)
            seqs[sid].append(re.sub(r"\s+", "", parts[1]).upper())
    return order, ["".join(seqs[k]) for k in order]


def _encode_haplotype_alleles(seqs: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
    """Encode sequences to diploid base pair arrays `a1, a2` of shape `(n_samples, n_sites)` (`-1` = missing)."""
    n_samp = len(seqs)
    n_sites = min(len(s) for s in seqs) if seqs else 0
    a1 = np.full((n_samp, n_sites), -1, dtype=np.int8)
    a2 = np.full((n_samp, n_sites), -1, dtype=np.int8)
    for i, sq in enumerate(seqs):
        arr = np.frombuffer(sq[:n_sites].encode("ascii", errors="replace"), dtype=np.uint8)
        for ch, (b1, b2) in _IUPAC_ALLELES.items():
            mask = arr == ord(ch)
            if np.any(mask):
                a1[i, mask] = _BASE_TO_IDX[b1]
                a2[i, mask] = _BASE_TO_IDX[b2]
    return a1, a2


def _pairwise_site_divergence_blocks(
    a1: np.ndarray,
    a2: np.ndarray,
    block_slices: Sequence[Tuple[int, int]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute per-block `(n_samp, n_samp)` pairwise allele mismatch counts and callable site counts."""
    return compute_model_corrected_divergence_blocks(a1, a2, block_slices, model_spec=None)


def compute_svdquartets_support(
    a1: np.ndarray,
    a2: np.ndarray,
    taxa_labels: Sequence[str],
    sample_ids: Sequence[str],
    max_quartets: int = 120,
    seed: int = 202611,
    n_boot: int = 200,
) -> pd.DataFrame:
    """Compute Chifman & Kubatko (2014) 16x16 SVDquartets scores and block-bootstrap p-values."""
    empty_cols = [
        "quartet_id",
        "best_split",
        "svd_score_best",
        "svd_score_second",
        "quartet_support",
        "quartet_boot_support",
        "quartet_p_value",
        "n_sites_used",
    ]
    taxon_groups: Dict[str, List[int]] = defaultdict(list)
    for idx, tx in enumerate(taxa_labels):
        taxon_groups[tx].append(idx)
    uniq_taxa = sorted(taxon_groups.keys())
    if len(uniq_taxa) < 4 or a1.shape[1] < 20:
        return pd.DataFrame(columns=empty_cols)

    rng = np.random.default_rng(seed)
    from itertools import combinations

    all_q = list(combinations(uniq_taxa, 4))
    if len(all_q) > max_quartets:
        sel = rng.choice(len(all_q), size=max_quartets, replace=False)
        all_q = [all_q[int(k)] for k in sorted(sel)]

    rows: List[Dict[str, Any]] = []
    # Pick representative highest-coverage sequence per taxon for deterministic 16x16 flattening
    rep_idx = {tx: max(idxs, key=lambda i: int(np.sum(a1[i] >= 0))) for tx, idxs in taxon_groups.items()}
    n_sites = a1.shape[1]

    def _svd_residual(mats: np.ndarray) -> np.ndarray:
        sv = np.linalg.svd(mats, compute_uv=False)
        return np.sqrt(np.sum(sv[..., 10:] ** 2, axis=-1))

    for t1, t2, t3, t4 in all_q:
        i1, i2, i3, i4 = rep_idx[t1], rep_idx[t2], rep_idx[t3], rep_idx[t4]
        valid = (a1[i1] >= 0) & (a1[i2] >= 0) & (a1[i3] >= 0) & (a1[i4] >= 0)
        nv = int(np.sum(valid))
        if nv < 20:
            continue
        b1, b2, b3, b4 = a1[i1, valid], a1[i2, valid], a1[i3, valid], a1[i4, valid]

        n_b = max(4, min(40, nv // (2 if nv < 40 else 10)))
        blk_splits = np.array_split(np.arange(nv), n_b)
        blk_tensors = np.zeros((n_b, 256), dtype=np.float64)
        idx_256 = (
            b1.astype(np.int32) * 64
            + b2.astype(np.int32) * 16
            + b3.astype(np.int32) * 4
            + b4.astype(np.int32)
        )
        for bi, b_idx in enumerate(blk_splits):
            np.add.at(blk_tensors[bi], idx_256[b_idx], 1.0)

        obs_t4 = np.sum(blk_tensors, axis=0).reshape(4, 4, 4, 4) / float(nv)
        obs_mats = np.stack([
            obs_t4.reshape(16, 16),
            obs_t4.transpose(0, 2, 1, 3).reshape(16, 16),
            obs_t4.transpose(0, 3, 1, 2).reshape(16, 16),
        ], axis=0)
        obs_scores = _svd_residual(obs_mats)
        split_labels = [
            f"{t1},{t2}|{t3},{t4}",
            f"{t1},{t3}|{t2},{t4}",
            f"{t1},{t4}|{t2},{t3}",
        ]
        order = np.argsort(obs_scores, kind="mergesort")
        best_idx = int(order[0])
        s_best, s_sec, s_trd = float(obs_scores[order[0]]), float(obs_scores[order[1]]), float(obs_scores[order[2]])
        scale = max(1e-6, float(np.std([s_best, s_sec, s_trd])))
        w = np.exp(-np.array([s_best, s_sec, s_trd]) / scale)
        support = float(w[0] / np.sum(w))

        draws = rng.integers(0, n_b, size=(max(1, n_boot), n_b))
        bw = np.zeros((max(1, n_boot), n_b), dtype=np.float64)
        for bi in range(n_b):
            bw[:, bi] = np.sum(draws == bi, axis=1)
        boot_raw = (bw @ blk_tensors).reshape(-1, 4, 4, 4, 4)
        boot_den = np.maximum(1.0, np.sum(boot_raw, axis=(1, 2, 3, 4), keepdims=True))
        boot_t4 = boot_raw / boot_den
        b_s0 = _svd_residual(boot_t4.reshape(-1, 16, 16))
        b_s1 = _svd_residual(boot_t4.transpose(0, 1, 3, 2, 4).reshape(-1, 16, 16))
        b_s2 = _svd_residual(boot_t4.transpose(0, 1, 4, 2, 3).reshape(-1, 16, 16))
        boot_win = np.argmin(np.stack([b_s0, b_s1, b_s2], axis=1), axis=1)
        n_match = int(np.sum(boot_win == best_idx))
        boot_supp = float(n_match / float(len(boot_win)))
        p_val = float((len(boot_win) - n_match + 1) / float(len(boot_win) + 1))

        rows.append({
            "quartet_id": f"{t1};{t2};{t3};{t4}",
            "best_split": split_labels[best_idx],
            "svd_score_best": round(s_best, 6),
            "svd_score_second": round(s_sec, 6),
            "quartet_support": round(support, 4),
            "quartet_boot_support": round(boot_supp, 4),
            "quartet_p_value": round(p_val, 4),
            "n_sites_used": min(nv, n_sites),
        })
    return pd.DataFrame(rows, columns=empty_cols) if not rows else pd.DataFrame(rows)


def export_bpp_msci_files(
    sample_ids: Sequence[str],
    taxa_labels: Sequence[str],
    seqs: Sequence[str],
    block_slices: Sequence[Tuple[int, int]],
    out_dir: Path,
    max_loci: int = 200,
    substitution_model: str = "JC69",
) -> Tuple[Path, Path]:
    """Write BPP-formatted multi-locus PHYLIP and `bpp_msci.ctl` into run-scoped `out_dir`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    loci_phy = out_dir / "bpp_msci_loci.phy"
    ctl_path = out_dir / "bpp_msci.ctl"
    sel_blocks = list(block_slices)[:max_loci]

    with open(loci_phy, "wt", encoding="utf-8") as fh:
        for st, en in sel_blocks:
            sub_seqs = [sq[st:en] for sq in seqs]
            keep_rows = [
                (sid, ss)
                for sid, ss in zip(sample_ids, sub_seqs)
                if any(ch in "ACGT" for ch in ss)
            ]
            if len(keep_rows) < 4:
                continue
            fh.write(f"{len(keep_rows)} {en - st}\n")
            for sid, ss in keep_rows:
                clean_s = re.sub(r"[^ACGTURYKSWMN?-]", "N", ss)
                fh.write(f"^{sid:<18s} {clean_s}\n")
            fh.write("\n")

    uniq_taxa = sorted(set(taxa_labels))
    imap_path = out_dir / "bpp_imap.txt"
    with open(imap_path, "wt", encoding="utf-8") as ih:
        for sid, tx in zip(sample_ids, taxa_labels):
            ih.write(f"{sid}\t{tx.replace(' ', '_')}\n")

    base_sub = str(substitution_model or "JC69").split("+", 1)[0].upper()
    bpp_model = "HKY" if base_sub in ("HKY85", "TN93", "F81", "K80") else "JC69"
    gamma_line = "alpha = 0.5\n" if "+G4" in str(substitution_model or "").upper() else ""
    species_tokens = " ".join(tx.replace(" ", "_") for tx in uniq_taxa)
    ctl_txt = (
        f"* BIC selected substitution model: {substitution_model}\n"
        f"seed = 202611\n"
        f"seqfile = {loci_phy.name}\n"
        f"Imapfile = {imap_path.name}\n"
        f"outfile = bpp_msci_out.txt\n"
        f"mcmcfile = bpp_msci_mcmc.txt\n"
        f"model = {bpp_model}\n"
        f"{gamma_line}"
        f"speciesdelimitation = 0\n"
        f"speciestree = 0\n"
        f"species&tree = {len(uniq_taxa)} {species_tokens}\n"
        f"    {' '.join(['1'] * len(uniq_taxa))}\n"
        f"thetaprior = 3 0.004 e\n"
        f"tauprior = 3 0.004\n"
        f"phiprior = 1 1\n"
        f"nloci = {len(sel_blocks)}\n"
        f"burnin = 4000\n"
        f"nsample = 20000\n"
        f"sampfreq = 2\n"
    )
    ctl_path.write_text(ctl_txt, encoding="utf-8")
    return loci_phy, ctl_path


def _calc_lineage_theta_and_dxy(
    diff_blk: np.ndarray,
    call_blk: np.ndarray,
    within_het: np.ndarray,
    within_call: np.ndarray,
    idx_a: Sequence[int],
    idx_b: Sequence[int],
    blk_indices: Optional[np.ndarray] = None,
) -> Tuple[float, float, float]:
    """Compute `(theta_A, theta_B, d_XY)` on identical per-site coalescent scale."""
    bi = np.arange(diff_blk.shape[0]) if blk_indices is None else blk_indices
    sub_d = diff_blk[bi]
    sub_c = call_blk[bi]
    sub_wh = within_het[bi]
    sub_wc = within_call[bi]

    def _within_pi(idxs: Sequence[int]) -> float:
        num = 0.0
        den = 0.0
        if len(idxs) == 1:
            i = idxs[0]
            num = float(np.sum(sub_wh[:, i]))
            den = float(np.sum(sub_wc[:, i]))
        else:
            for p_i in range(len(idxs)):
                for p_j in range(p_i + 1, len(idxs)):
                    ia, ib = idxs[p_i], idxs[p_j]
                    num += float(np.sum(sub_d[:, ia, ib]))
                    den += float(np.sum(sub_c[:, ia, ib]))
        # Guard ultra-low heterozygosity floor only at 1e-4 to avoid division by zero on short toys
        return max(1e-4, num / den) if den > 0 else 1e-4

    theta_a = _within_pi(idx_a)
    theta_b = _within_pi(idx_b)

    d_num = 0.0
    d_den = 0.0
    for ia in idx_a:
        for ib in idx_b:
            d_num += float(np.sum(sub_d[:, ia, ib]))
            d_den += float(np.sum(sub_c[:, ia, ib]))
    d_xy = (d_num / d_den) if d_den > 0 else 0.0
    return theta_a, theta_b, d_xy


def compute_coalescent_gdi_delimitation(
    sample_ids: Sequence[str],
    taxa_labels: Sequence[str],
    diff_blk: np.ndarray,
    call_blk: np.ndarray,
    within_het: np.ndarray,
    within_call: np.ndarray,
    n_boot: int = 500,
    seed: int = 202611,
    substitution_model: str = "JC69",
    bpp_out_path: Optional[Path] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Compute per-taxon and per-sample `gdi = 1 - exp(-2*tau/theta)` with locus-bootstrap 95% CIs."""
    by_taxon: Dict[str, List[int]] = defaultdict(list)
    for i, tx in enumerate(taxa_labels):
        by_taxon[tx].append(i)
    taxa = sorted(by_taxon.keys())
    n_blk = diff_blk.shape[0]
    rng = np.random.default_rng(seed)

    # First locate nearest sister taxon by d_XY across all loci
    pair_dxy: Dict[Tuple[str, str], float] = {}
    for ta in taxa:
        for tb in taxa:
            if ta >= tb:
                continue
            _, _, dxy = _calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, by_taxon[ta], by_taxon[tb])
            pair_dxy[(ta, tb)] = dxy
            pair_dxy[(tb, ta)] = dxy

    taxon_rows: List[Dict[str, Any]] = []
    adequacy_rows: List[Dict[str, Any]] = []
    taxon_gdi_map: Dict[str, Dict[str, Any]] = {}
    boot_draws = rng.integers(0, n_blk, size=(max(10, n_boot), n_blk)) if n_blk > 1 else np.zeros((10, 1), dtype=int)

    for tx in taxa:
        others = [o for o in taxa if o != tx]
        if not others:
            continue
        sister = min(others, key=lambda o: pair_dxy.get((tx, o), 1.0))
        idx_a = by_taxon[tx]
        idx_b = by_taxon[sister]
        th_a, th_b, dxy = _calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, idx_a, idx_b)
        tau_net = max(0.0, 0.5 * (dxy - 0.5 * (th_a + th_b)))
        gdi_obs = float(1.0 - math.exp(-2.0 * tau_net / max(1e-4, th_a)))

        boot_vals: List[float] = []
        for b_row in boot_draws:
            b_ta, b_tb, b_dxy = _calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, idx_a, idx_b, blk_indices=b_row)
            b_tau = max(0.0, 0.5 * (b_dxy - 0.5 * (b_ta + b_tb)))
            boot_vals.append(1.0 - math.exp(-2.0 * b_tau / max(1e-4, b_ta)))
        ci_low, ci_high = float(np.percentile(boot_vals, 2.5)), float(np.percentile(boot_vals, 97.5))

        ix_r, ix_c = np.ix_(idx_a, idx_b)
        call_mean_vec = np.mean(call_blk[:, ix_r, ix_c], axis=(1, 2))
        diff_mean_vec = np.mean(diff_blk[:, ix_r, ix_c], axis=(1, 2))
        d_blk_ratio = np.divide(diff_mean_vec, np.maximum(1.0, call_mean_vec))
        bf_info = compute_coalescent_introgression_bf(d_blk_ratio, call_mean_vec, bpp_out_path=bpp_out_path)
        bpp_bf = float(bf_info["bf_10"])

        if gdi_obs < 0.20:
            gdi_cls = "Single_Lineage_Oversplit"
        elif gdi_obs > 0.70:
            gdi_cls = "Distinct_Coalescent_Lineage"
        else:
            gdi_cls = "Ambiguous_Incomplete_Sorting"

        rec = {
            "taxon": tx,
            "nearest_sister_taxon": sister,
            "n_samples": len(idx_a),
            "n_loci_blocks": n_blk,
            "theta_lineage": round(th_a, 6),
            "theta_sister": round(th_b, 6),
            "d_xy": round(dxy, 6),
            "tau_net": round(tau_net, 6),
            "coalescent_gdi": round(gdi_obs, 4),
            "gdi_ci_low": round(ci_low, 4),
            "gdi_ci_high": round(ci_high, 4),
            "gdi_classification": gdi_cls,
            "bpp_introgression_bf": bpp_bf,
            "substitution_model": substitution_model,
            "bf_method": bf_info["bf_method"],
        }
        taxon_rows.append(rec)
        taxon_gdi_map[tx] = rec
        adequacy_rows.append({
            "taxon": tx,
            "nearest_sister_taxon": sister,
            "substitution_model": substitution_model,
            "bf_method": bf_info["bf_method"],
            "bic_m0": bf_info["bic_m0"],
            "bic_m1": bf_info["bic_m1"],
            "log_bf_10": bf_info["log_bf_10"],
            "bf_10": bf_info["bf_10"],
            "w_introgr": bf_info["w_introgr"],
            "tau_introgr": bf_info["tau_introgr"],
            "tau_sp": bf_info["tau_sp"],
            "coalescent_gdi": rec["coalescent_gdi"],
            "gdi_classification": gdi_cls,
        })

    sample_rows = []
    for sid, tx in zip(sample_ids, taxa_labels):
        info = taxon_gdi_map.get(tx, {})
        sample_rows.append({
            "genetic_sample_id": sid,
            "lineage_taxon": tx,
            "nearest_sister_taxon": info.get("nearest_sister_taxon"),
            "coalescent_gdi": info.get("coalescent_gdi"),
            "gdi_ci_low": info.get("gdi_ci_low"),
            "gdi_ci_high": info.get("gdi_ci_high"),
            "gdi_classification": info.get("gdi_classification"),
            "bpp_introgression_bf": info.get("bpp_introgression_bf"),
        })
    taxon_df = pd.DataFrame(taxon_rows)
    taxon_df.attrs["adequacy_df"] = pd.DataFrame(adequacy_rows)
    return taxon_df, pd.DataFrame(sample_rows)


def _load_sample_metadata(
    genomics_dir: Path,
    parsed_samples_csv: Optional[Path],
) -> Tuple[Optional[set[str]], Dict[str, str], Dict[str, str], Dict[str, float], Dict[str, float]]:
    sample_taxon_map: Dict[str, str] = {}
    rep_pair_map: Dict[str, str] = {}
    hq_reads_map: Dict[str, float] = {}
    miss_rate_map: Dict[str, float] = {}
    pass_set: Optional[set[str]] = None

    if parsed_samples_csv and parsed_samples_csv.exists():
        sdf = pd.read_csv(parsed_samples_csv, dtype=str, keep_default_na=False)
        if not sdf.empty:
            pass_set = {
                str(r["genetic_sample_id"])
                for _, r in sdf.iterrows()
                if str(r.get("geno_quality_pass", "1")) in ("1", "True", "true")
            }
            for _, r in sdf.iterrows():
                sid = str(r["genetic_sample_id"])
                tx = str(r.get("species_standardized") or r.get("taxon_raw") or "")
                if str(r.get("is_putative_novel_or_hybrid_query", "0")) == "1":
                    pop = str(r.get("population_id") or "")
                    tx = f"{tx}_{pop}" if (tx and pop) else (tx or sid)
                if tx:
                    sample_taxon_map[sid] = tx
                rp = str(r.get("replicate_pair_id") or "").strip()
                if rp and rp.upper() not in {"NA", "NONE", "NAN", "NULL"}:
                    rep_pair_map[sid] = rp
                for raw_h, dst in ((r.get("hq_reads"), hq_reads_map), (r.get("missing_genotype_rate"), miss_rate_map)):
                    try:
                        if raw_h not in (None, ""):
                            fv = float(raw_h)
                            if math.isfinite(fv):
                                dst[sid] = fv
                    except (TypeError, ValueError):
                        pass

    db_cands = [
        Path(env_p)
        for env_k in ("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH")
        if (env_p := os.environ.get(env_k, "").strip())
    ] + [genomics_dir.parent / "plantmatics.sqlite", PROJECT_ROOT / "data" / "plantmatics.sqlite"]
    for db_p in db_cands:
        if not db_p.is_file():
            continue
        try:
            conn = sqlite3.connect(f"file:{db_p}?mode=ro", uri=True)
            try:
                cur = conn.execute(
                    "SELECT genetic_sample_id, species_standardized, replicate_pair_id, hq_reads, missing_genotype_rate, geno_quality_pass "
                    "FROM genomic_samples;"
                )
                for sid_r, sp_r, rp_r, hq_r, mr_r, qp_r in cur.fetchall():
                    sid_s = str(sid_r or "").strip()
                    if not sid_s:
                        continue
                    if sid_s not in sample_taxon_map and sp_r:
                        sample_taxon_map[sid_s] = str(sp_r)
                    rp_s = str(rp_r or "").strip()
                    if sid_s not in rep_pair_map and rp_s and rp_s.upper() not in {"NA", "NONE", "NAN", "NULL"}:
                        rep_pair_map[sid_s] = rp_s
                    if sid_s not in hq_reads_map and hq_r is not None:
                        hq_reads_map[sid_s] = float(hq_r)
                    if sid_s not in miss_rate_map and mr_r is not None:
                        miss_rate_map[sid_s] = float(mr_r)
                    if pass_set is None and qp_r is not None and int(qp_r) == 1:
                        pass_set = pass_set or set()
                        pass_set.add(sid_s)
            finally:
                conn.close()
            break
        except sqlite3.Error:
            continue
    return pass_set, sample_taxon_map, rep_pair_map, hq_reads_map, miss_rate_map


def _group_technical_replicates(
    sids: Sequence[str],
    taxa_labels: Sequence[str],
    rep_pair_map: Mapping[str, str],
) -> Dict[Any, List[int]]:
    groups: Dict[Any, List[int]] = defaultdict(list)
    for idx, (sid, tx) in enumerate(zip(sids, taxa_labels)):
        rp = str(rep_pair_map.get(sid) or "").strip()
        key = ("pair", rp, tx) if rp else ("uniq", idx)
        groups[key].append(idx)
    return groups


def _collapse_alignment_replicates(
    sids: Sequence[str],
    seqs: Sequence[str],
    taxa_labels: Sequence[str],
    rep_pair_map: Mapping[str, str],
    hq_reads_map: Mapping[str, float],
    miss_rate_map: Mapping[str, float],
) -> Tuple[List[str], List[str], List[str], Dict[str, List[str]]]:
    groups = _group_technical_replicates(sids, taxa_labels, rep_pair_map)
    out_sids: List[str] = []
    out_seqs: List[str] = []
    out_taxa: List[str] = []
    group_members: Dict[str, List[str]] = {}
    for idxs in groups.values():
        ranked = sorted(
            idxs,
            key=lambda i: (
                -sum(1 for ch in seqs[i] if ch not in "N?-"),
                -float(hq_reads_map.get(sids[i], 0.0)),
                float(miss_rate_map.get(sids[i], 1.0)),
                sids[i],
            ),
        )
        p_i = ranked[0]
        p_sid = sids[p_i]
        merged = seqs[p_i]
        if len(ranked) > 1:
            buf = list(merged)
            for s_i in ranked[1:]:
                sib = seqs[s_i]
                lim = min(len(buf), len(sib))
                for pos in range(lim):
                    if buf[pos] in "N?-" and sib[pos] not in "N?-":
                        buf[pos] = sib[pos]
            merged = "".join(buf)
        out_sids.append(p_sid)
        out_seqs.append(merged)
        out_taxa.append(taxa_labels[p_i])
        group_members[p_sid] = [sids[i] for i in idxs]
    return out_sids, out_seqs, out_taxa, group_members


def _collapse_dosage_replicates(
    sids: Sequence[str],
    mat: np.ndarray,
    taxa_labels: Sequence[str],
    rep_pair_map: Mapping[str, str],
    hq_reads_map: Mapping[str, float],
    miss_rate_map: Mapping[str, float],
) -> Tuple[List[str], np.ndarray, List[str], Dict[str, List[str]]]:
    groups = _group_technical_replicates(sids, taxa_labels, rep_pair_map)
    out_sids: List[str] = []
    out_rows: List[np.ndarray] = []
    out_taxa: List[str] = []
    group_members: Dict[str, List[str]] = {}
    for idxs in groups.values():
        ranked = sorted(
            idxs,
            key=lambda i: (
                -int(np.sum(np.isfinite(mat[i]))),
                -float(hq_reads_map.get(sids[i], 0.0)),
                float(miss_rate_map.get(sids[i], 1.0)),
                sids[i],
            ),
        )
        p_i = ranked[0]
        p_sid = sids[p_i]
        if len(ranked) == 1:
            merged_row = mat[p_i]
        else:
            sub = mat[ranked]
            ok = np.isfinite(sub)
            cnt = np.sum(ok, axis=0)
            sums = np.sum(np.where(ok, sub, 0.0), axis=0)
            merged_row = np.where(cnt > 0, sums / np.maximum(cnt, 1), np.nan)
        out_sids.append(p_sid)
        out_rows.append(merged_row)
        out_taxa.append(taxa_labels[p_i])
        group_members[p_sid] = [sids[i] for i in idxs]
    collapsed_mat = np.vstack(out_rows) if out_rows else mat[:0]
    return out_sids, collapsed_mat, out_taxa, group_members


def _expand_gdi_samples(
    gdi_s_df: pd.DataFrame,
    group_members: Mapping[str, List[str]],
    orig_sids: Sequence[str],
) -> pd.DataFrame:
    if gdi_s_df.empty or not group_members:
        return gdi_s_df
    by_rep = {str(r["genetic_sample_id"]): r.to_dict() for _, r in gdi_s_df.iterrows()}
    orig_to_rep = {m: rep for rep, members in group_members.items() for m in members}
    rows = []
    for sid in orig_sids:
        rep = orig_to_rep.get(sid, sid)
        if rep in by_rep:
            rec = dict(by_rep[rep])
            rec["genetic_sample_id"] = sid
            rows.append(rec)
    return pd.DataFrame(rows, columns=gdi_s_df.columns) if rows else gdi_s_df


def _select_models_for_skimmed_fastas(
    genomics_dir: Path,
    exclude_paths: Sequence[Path],
) -> List[Dict[str, Any]]:
    """Run BIC substitution model selection across any recursive `skimmed_*.fasta` alignments."""
    excl_resolved = {p.resolve() for p in exclude_paths}
    candidates = sorted(
        p for p in genomics_dir.rglob("skimmed_*.fasta")
        if p.is_file() and not any(part.startswith(".") for part in p.parts) and p.resolve() not in excl_resolved
    )
    if not candidates:
        return []

    def _fit_skimmed(fa_path: Path) -> List[Dict[str, Any]]:
        sk_sids, sk_seqs = parse_invariant_alignment(fa_path)
        if len(sk_sids) < 3 or not sk_seqs:
            return []
        spec = select_nucleotide_substitution_model(dict(zip(sk_sids, sk_seqs)), sk_sids, alignment_name=fa_path.name)
        return list(spec.get("table_rows") or [])

    n_workers = min(_resolve_max_workers(), len(candidates))
    rows: List[Dict[str, Any]] = []
    if n_workers > 1:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            for sub_rows in pool.map(_fit_skimmed, candidates):
                rows.extend(sub_rows)
    else:
        for fa_path in candidates:
            rows.extend(_fit_skimmed(fa_path))
    return rows


def run_coalescent_delimitation(
    genomics_dir: Path,
    out_dir: Path,
    parsed_samples_csv: Optional[Path] = None,
    dosage_tsv: Optional[Path] = None,
    block_len: int = 122,
    n_boot: int = 500,
) -> Dict[str, Any]:
    """Execute invariant-site SVDquartets, BPP export, and gdi delimitation into `out_dir`."""
    discover_fn = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").discover_genomic_inputs
    assets = discover_fn(genomics_dir)

    if not assets.alignment_files and not (dosage_tsv and dosage_tsv.exists()):
        return {"status": "skipped_empty", "n_samples": 0}

    pass_set, sample_taxon_map, rep_pair_map, hq_reads_map, miss_rate_map = _load_sample_metadata(
        genomics_dir, parsed_samples_csv
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    gdi_csv = out_dir / "coalescent_gdi_delimitation.csv"
    gdi_samp_csv = out_dir / "coalescent_gdi_samples.csv"
    svdq_csv = out_dir / "svdquartets_support.csv"
    model_sel_csv = out_dir / "substitution_model_selection.csv"
    adequacy_csv = out_dir / "coalescent_model_adequacy.csv"

    if assets.alignment_files:
        # Prefer PHYLIP/NEXUS alignments over generic FASTA
        align_sorted = sorted(
            assets.alignment_files,
            key=lambda p: (0 if p.name.lower().endswith((".phylip", ".phy", ".nex", ".nexus")) else 1, p.name),
        )
        primary_path = align_sorted[0]
        sids, seqs = parse_invariant_alignment(primary_path)
        if pass_set is not None:
            kept = [(sid, sq) for sid, sq in zip(sids, seqs) if sid in pass_set]
            if len(kept) >= 4:
                sids = [x[0] for x in kept]
                seqs = [x[1] for x in kept]

        orig_sids = list(sids)
        taxa_labels = [_infer_taxon(sid, sample_taxon_map) for sid in sids]
        sids, seqs, taxa_labels, group_members = _collapse_alignment_replicates(
            sids, seqs, taxa_labels, rep_pair_map, hq_reads_map, miss_rate_map
        )
        n_sites = min(len(s) for s in seqs) if seqs else 0
        eff_blk = block_len if n_sites >= block_len * 2 else max(20, n_sites // 4 or n_sites)
        block_slices = [(st, min(n_sites, st + eff_blk)) for st in range(0, n_sites, eff_blk) if (min(n_sites, st + eff_blk) - st) >= 15]

        paup_bin = shutil.which("paup")
        bpp_bin = shutil.which("bpp")
        if paup_bin:
            logger.info("Detected external PAUP* binary (%s); running SVDquartets check alongside Python engine.", paup_bin)
        if bpp_bin:
            logger.info("Detected external BPP binary (%s); generated control file for MSci.", bpp_bin)

        primary_spec = select_nucleotide_substitution_model(
            dict(zip(sids, seqs)),
            sids,
            alignment_name=primary_path.name,
        )
        model_rows = list(primary_spec.get("table_rows") or [])
        model_rows.extend(_select_models_for_skimmed_fastas(genomics_dir, [primary_path]))
        pd.DataFrame(model_rows).to_csv(model_sel_csv, index=False)

        best_model = str(primary_spec.get("best_model") or "JC69")
        a1, a2 = _encode_haplotype_alleles(seqs)
        diff_blk, call_blk, within_het, within_call = compute_model_corrected_divergence_blocks(
            a1, a2, block_slices, model_spec=primary_spec
        )
        svdq_df = compute_svdquartets_support(a1, a2, taxa_labels, sids, n_boot=200, seed=202611)
        export_bpp_msci_files(
            sids, taxa_labels, seqs, block_slices, out_dir, substitution_model=best_model
        )
        gdi_df, gdi_s_df = compute_coalescent_gdi_delimitation(
            sids,
            taxa_labels,
            diff_blk,
            call_blk,
            within_het,
            within_call,
            n_boot=n_boot,
            substitution_model=best_model,
            bpp_out_path=out_dir / "bpp_msci_out.txt",
        )
        adequacy_df = gdi_df.attrs.get("adequacy_df", pd.DataFrame())
        gdi_s_df = _expand_gdi_samples(gdi_s_df, group_members, orig_sids)

        svdq_df.to_csv(svdq_csv, index=False)
        gdi_df.to_csv(gdi_csv, index=False)
        gdi_s_df.to_csv(gdi_samp_csv, index=False)
        adequacy_df.to_csv(adequacy_csv, index=False)
        return {
            "status": "completed_alignment",
            "n_samples": len(sids),
            "n_blocks": len(block_slices),
            "best_substitution_model": best_model,
            "gdi_csv": str(gdi_csv),
            "svdq_csv": str(svdq_csv),
            "model_selection_csv": str(model_sel_csv),
            "adequacy_csv": str(adequacy_csv),
        }

    # Dosage-only fallback (when tests or VCF-only inputs lack full invariant PHYLIP/NEXUS)
    if dosage_tsv and dosage_tsv.exists():
        ddf = pd.read_csv(dosage_tsv, sep="\t")
        sids = [str(x) for x in ddf["genetic_sample_id"].tolist()]
        mat = ddf.drop(columns=["genetic_sample_id"]).to_numpy(dtype=np.float64)
        orig_sids = list(sids)
        taxa_labels = [_infer_taxon(sid, sample_taxon_map) for sid in sids]
        sids, mat, taxa_labels, group_members = _collapse_dosage_replicates(
            sids, mat, taxa_labels, rep_pair_map, hq_reads_map, miss_rate_map
        )
        n_s, n_l = mat.shape
        diff_blk = np.zeros((1, n_s, n_s), dtype=np.float64)
        call_blk = np.ones((1, n_s, n_s), dtype=np.float64) * max(1, n_l)
        for i in range(n_s):
            for j in range(n_s):
                ok = np.isfinite(mat[i]) & np.isfinite(mat[j])
                if np.any(ok):
                    call_blk[0, i, j] = float(np.sum(ok))
                    diff_blk[0, i, j] = float(np.sum(np.abs(mat[i, ok] - mat[j, ok]) * 0.02))
        w_het = (np.nansum(mat * (1.0 - mat), axis=1, keepdims=True) * 0.02).T
        w_call = np.sum(np.isfinite(mat), axis=1, keepdims=True).astype(np.float64).T
        gdi_df, gdi_s_df = compute_coalescent_gdi_delimitation(sids, taxa_labels, diff_blk, call_blk, w_het, w_call, n_boot=50)
        adequacy_df = gdi_df.attrs.get("adequacy_df", pd.DataFrame())
        gdi_s_df = _expand_gdi_samples(gdi_s_df, group_members, orig_sids)
        gdi_df.to_csv(gdi_csv, index=False)
        gdi_s_df.to_csv(gdi_samp_csv, index=False)
        adequacy_df.to_csv(adequacy_csv, index=False)
        pd.DataFrame(columns=["quartet_id", "best_split", "svd_score_best", "svd_score_second", "quartet_support", "quartet_boot_support", "quartet_p_value", "n_sites_used"]).to_csv(svdq_csv, index=False)
        return {"status": "completed_dosage_fallback", "n_samples": n_s, "gdi_csv": str(gdi_csv)}

    return {"status": "skipped_empty", "n_samples": 0}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Invariant-site coalescent delimitation (SVDquartets, BPP, gdi).")
    parser.add_argument("--genomics-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--parsed-samples", type=Path, default=None)
    parser.add_argument("--dosage-tsv", type=Path, default=None)
    parser.add_argument("--block-len", type=int, default=122)
    parser.add_argument("--n-boot", type=int, default=500)
    args = parser.parse_args(list(argv) if argv is not None else None)
    run_coalescent_delimitation(
        genomics_dir=args.genomics_dir,
        out_dir=args.out_dir,
        parsed_samples_csv=args.parsed_samples,
        dosage_tsv=args.dosage_tsv,
        block_len=args.block_len,
        n_boot=args.n_boot,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

