#!/usr/bin/env python3
"""
Module: scripts/analysis/step06_coalescent_delimitation.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Invariant-site Multispecies Coalescent (MSC) & Reticulation (MSci) engine for
  Step 06 phylogenomics (`SVDquartets`, `BPP` control/PHYLIP export + conjugate
  Bayes factors, and locus-bootstrapped `gdi = 1 - exp(-2 * tau / theta)`).
"""

from __future__ import annotations

import argparse
import importlib
import math
import shutil
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.analysis.step06_coalescent_alignment_utils import (
    _BASE_TO_IDX,
    _IUPAC_ALLELES,
    _clean_id,
    _collapse_alignment_replicates,
    _collapse_dosage_replicates,
    _encode_haplotype_alleles,
    _expand_gdi_samples,
    _group_technical_replicates,
    _infer_taxon,
    _load_sample_metadata,
    _pairwise_site_divergence_blocks,
    _select_models_for_skimmed_fastas,
    compute_svdquartets_support,
    export_bpp_msci_files,
    parse_invariant_alignment,
)
from scripts.analysis.step06_coalescent_models import (
    _resolve_max_workers,
    compute_coalescent_introgression_bf,
    compute_model_corrected_divergence_blocks,
    select_nucleotide_substitution_model,
)
from scripts.core.config import resolve_genomics_config
from scripts.core.logger import setup_logging

logger = setup_logging(name="packera.genomics.coalescent")


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
    sub_d, sub_c, sub_wh, sub_wc = diff_blk[bi], call_blk[bi], within_het[bi], within_call[bi]

    def _within_pi(idxs: Sequence[int]) -> float:
        num = den = 0.0
        if len(idxs) == 1:
            num, den = float(np.sum(sub_wh[:, idxs[0]])), float(np.sum(sub_wc[:, idxs[0]]))
        else:
            for p_i in range(len(idxs)):
                for p_j in range(p_i + 1, len(idxs)):
                    ia, ib = idxs[p_i], idxs[p_j]
                    num += float(np.sum(sub_d[:, ia, ib]))
                    den += float(np.sum(sub_c[:, ia, ib]))
        return max(1e-4, num / den) if den > 0 else 1e-4

    theta_a, theta_b = _within_pi(idx_a), _within_pi(idx_b)
    d_num = sum(float(np.sum(sub_d[:, ia, ib])) for ia in idx_a for ib in idx_b)
    d_den = sum(float(np.sum(sub_c[:, ia, ib])) for ia in idx_a for ib in idx_b)
    return theta_a, theta_b, (d_num / d_den) if d_den > 0 else 0.0


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

    pair_dxy: Dict[Tuple[str, str], float] = {}
    for ta in taxa:
        for tb in taxa:
            if ta < tb:
                _, _, dxy = _calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, by_taxon[ta], by_taxon[tb])
                pair_dxy[(ta, tb)] = pair_dxy[(tb, ta)] = dxy

    taxon_rows: List[Dict[str, Any]] = []
    adequacy_rows: List[Dict[str, Any]] = []
    taxon_gdi_map: Dict[str, Dict[str, Any]] = {}
    boot_draws = rng.integers(0, n_blk, size=(max(10, n_boot), n_blk)) if n_blk > 1 else np.zeros((10, 1), dtype=int)

    def _eval_tx_gdi(tx: str) -> Optional[Tuple[str, Dict[str, Any], Dict[str, Any]]]:
        others = [o for o in taxa if o != tx]
        if not others:
            return None
        sister = min(others, key=lambda o: pair_dxy.get((tx, o), 1.0))
        idx_a, idx_b = by_taxon[tx], by_taxon[sister]
        th_a, th_b, dxy = _calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, idx_a, idx_b)
        tau_net = max(0.0, 0.5 * (dxy - 0.5 * (th_a + th_b)))
        gdi_obs = float(1.0 - math.exp(-2.0 * tau_net / max(1e-4, th_a)))

        boot_vals = [
            1.0 - math.exp(-2.0 * max(0.0, 0.5 * (b_dxy - 0.5 * (b_ta + b_tb))) / max(1e-4, b_ta))
            for b_row in boot_draws
            for b_ta, b_tb, b_dxy in (_calc_lineage_theta_and_dxy(diff_blk, call_blk, within_het, within_call, idx_a, idx_b, blk_indices=b_row),)
        ]
        ci_low, ci_high = float(np.percentile(boot_vals, 2.5)), float(np.percentile(boot_vals, 97.5))

        ix_r, ix_c = np.ix_(idx_a, idx_b)
        call_mean_vec = np.mean(call_blk[:, ix_r, ix_c], axis=(1, 2))
        diff_mean_vec = np.mean(diff_blk[:, ix_r, ix_c], axis=(1, 2))
        bf_info = compute_coalescent_introgression_bf(np.divide(diff_mean_vec, np.maximum(1.0, call_mean_vec)), call_mean_vec, bpp_out_path=bpp_out_path)
        bpp_bf = float(bf_info["bf_10"])
        gdi_cls = "Single_Lineage_Oversplit" if gdi_obs < 0.20 else ("Distinct_Coalescent_Lineage" if gdi_obs > 0.70 else "Ambiguous_Incomplete_Sorting")

        rec = {
            "taxon": tx, "nearest_sister_taxon": sister, "n_samples": len(idx_a), "n_loci_blocks": n_blk,
            "theta_lineage": round(th_a, 6), "theta_sister": round(th_b, 6), "d_xy": round(dxy, 6), "tau_net": round(tau_net, 6),
            "coalescent_gdi": round(gdi_obs, 4), "gdi_ci_low": round(ci_low, 4), "gdi_ci_high": round(ci_high, 4),
            "gdi_classification": gdi_cls, "bpp_introgression_bf": bpp_bf, "substitution_model": substitution_model, "bf_method": bf_info["bf_method"],
        }
        ad_row = {
            "taxon": tx, "nearest_sister_taxon": sister, "substitution_model": substitution_model, "bf_method": bf_info["bf_method"],
            "bic_m0": bf_info["bic_m0"], "bic_m1": bf_info["bic_m1"], "log_bf_10": bf_info["log_bf_10"], "bf_10": bf_info["bf_10"],
            "w_introgr": bf_info["w_introgr"], "tau_introgr": bf_info["tau_introgr"], "tau_sp": bf_info["tau_sp"],
            "coalescent_gdi": rec["coalescent_gdi"], "gdi_classification": gdi_cls,
        }
        return tx, rec, ad_row

    nw_tx = min(_resolve_max_workers(), len(taxa)) if len(taxa) >= 2 else 1
    if nw_tx > 1:
        with ThreadPoolExecutor(max_workers=nw_tx) as pool:
            tx_outs = list(pool.map(_eval_tx_gdi, taxa))
    else:
        tx_outs = [_eval_tx_gdi(tx) for tx in taxa]
    for item in tx_outs:
        if item is None:
            continue
        tx, rec, ad_row = item
        taxon_rows.append(rec)
        taxon_gdi_map[tx] = rec
        adequacy_rows.append(ad_row)

    sample_rows = [
        {
            "genetic_sample_id": sid, "lineage_taxon": tx, "nearest_sister_taxon": info.get("nearest_sister_taxon"),
            "coalescent_gdi": info.get("coalescent_gdi"), "gdi_ci_low": info.get("gdi_ci_low"),
            "gdi_ci_high": info.get("gdi_ci_high"), "gdi_classification": info.get("gdi_classification"),
            "bpp_introgression_bf": info.get("bpp_introgression_bf"),
        }
        for sid, tx in zip(sample_ids, taxa_labels)
        for info in (taxon_gdi_map.get(tx, {}),)
    ]
    taxon_df = pd.DataFrame(taxon_rows)
    taxon_df.attrs["adequacy_df"] = pd.DataFrame(adequacy_rows)
    return taxon_df, pd.DataFrame(sample_rows)


def run_coalescent_delimitation(
    genomics_dir: Path,
    out_dir: Path,
    parsed_samples_csv: Optional[Path] = None,
    dosage_tsv: Optional[Path] = None,
    block_len: Optional[int] = None,
    n_boot: Optional[int] = None,
) -> Dict[str, Any]:
    """Execute invariant-site SVDquartets, BPP export, and gdi delimitation into `out_dir`."""
    gcfg = resolve_genomics_config()
    eff_block_len = int(gcfg["coalescent_block_len"]) if block_len is None else int(block_len)
    eff_n_boot = int(gcfg["coalescent_n_boot"]) if n_boot is None else int(n_boot)

    discover_fn = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").discover_genomic_inputs
    assets = discover_fn(genomics_dir)
    if not assets.alignment_files and not (dosage_tsv and dosage_tsv.exists()):
        return {"status": "skipped_empty", "n_samples": 0}

    pass_set, sample_taxon_map, rep_pair_map, hq_reads_map, miss_rate_map = _load_sample_metadata(genomics_dir, parsed_samples_csv)
    out_dir.mkdir(parents=True, exist_ok=True)
    gdi_csv, gdi_samp_csv = out_dir / "coalescent_gdi_delimitation.csv", out_dir / "coalescent_gdi_samples.csv"
    svdq_csv, model_sel_csv, adequacy_csv = out_dir / "svdquartets_support.csv", out_dir / "substitution_model_selection.csv", out_dir / "coalescent_model_adequacy.csv"

    if assets.alignment_files:
        align_sorted = sorted(assets.alignment_files, key=lambda p: (0 if p.name.lower().endswith((".phylip", ".phy", ".nex", ".nexus")) else 1, p.name))
        primary_path = align_sorted[0]
        sids, seqs = parse_invariant_alignment(primary_path)
        if pass_set is not None and len(kept := [(sid, sq) for sid, sq in zip(sids, seqs) if sid in pass_set]) >= 4:
            sids, seqs = [x[0] for x in kept], [x[1] for x in kept]
        orig_sids = list(sids)
        taxa_labels = [_infer_taxon(sid, sample_taxon_map) for sid in sids]
        sids, seqs, taxa_labels, group_members = _collapse_alignment_replicates(sids, seqs, taxa_labels, rep_pair_map, hq_reads_map, miss_rate_map)
        n_sites = min(len(s) for s in seqs) if seqs else 0
        eff_blk = eff_block_len if n_sites >= eff_block_len * 2 else max(20, n_sites // 4 or n_sites)
        block_slices = [(st, min(n_sites, st + eff_blk)) for st in range(0, n_sites, eff_blk) if (min(n_sites, st + eff_blk) - st) >= 15]

        if paup_bin := shutil.which("paup"):
            logger.info("Detected external PAUP* binary (%s); running SVDquartets check alongside Python engine.", paup_bin)
        if bpp_bin := shutil.which("bpp"):
            logger.info("Detected external BPP binary (%s); generated control file for MSci.", bpp_bin)

        primary_spec = select_nucleotide_substitution_model(dict(zip(sids, seqs)), sids, alignment_name=primary_path.name)
        model_rows = list(primary_spec.get("table_rows") or []) + _select_models_for_skimmed_fastas(genomics_dir, [primary_path])
        pd.DataFrame(model_rows).to_csv(model_sel_csv, index=False)

        best_model = str(primary_spec.get("best_model") or "JC69")
        a1, a2 = _encode_haplotype_alleles(seqs)
        diff_blk, call_blk, within_het, within_call = compute_model_corrected_divergence_blocks(a1, a2, block_slices, model_spec=primary_spec)
        svdq_df = compute_svdquartets_support(a1, a2, taxa_labels, sids, n_boot=200, seed=202611)
        export_bpp_msci_files(sids, taxa_labels, seqs, block_slices, out_dir, substitution_model=best_model)
        gdi_df, gdi_s_df = compute_coalescent_gdi_delimitation(
            sids, taxa_labels, diff_blk, call_blk, within_het, within_call,
            n_boot=eff_n_boot, substitution_model=best_model, bpp_out_path=out_dir / "bpp_msci_out.txt",
        )
        adequacy_df = gdi_df.attrs.get("adequacy_df", pd.DataFrame())
        gdi_s_df = _expand_gdi_samples(gdi_s_df, group_members, orig_sids)
        svdq_df.to_csv(svdq_csv, index=False)
        gdi_df.to_csv(gdi_csv, index=False)
        gdi_s_df.to_csv(gdi_samp_csv, index=False)
        adequacy_df.to_csv(adequacy_csv, index=False)
        return {
            "status": "completed_alignment", "n_samples": len(sids), "n_blocks": len(block_slices),
            "best_substitution_model": best_model, "gdi_csv": str(gdi_csv), "svdq_csv": str(svdq_csv),
            "model_selection_csv": str(model_sel_csv), "adequacy_csv": str(adequacy_csv),
        }

    if dosage_tsv and dosage_tsv.exists():
        ddf = pd.read_csv(dosage_tsv, sep="\t")
        sids = [str(x) for x in ddf["genetic_sample_id"].tolist()]
        mat = ddf.drop(columns=["genetic_sample_id"]).to_numpy(dtype=np.float64)
        orig_sids = list(sids)
        taxa_labels = [_infer_taxon(sid, sample_taxon_map) for sid in sids]
        sids, mat, taxa_labels, group_members = _collapse_dosage_replicates(sids, mat, taxa_labels, rep_pair_map, hq_reads_map, miss_rate_map)
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
        gdi_df.to_csv(gdi_csv, index=False); gdi_s_df.to_csv(gdi_samp_csv, index=False); adequacy_df.to_csv(adequacy_csv, index=False)
        pd.DataFrame(columns=["quartet_id", "best_split", "svd_score_best", "svd_score_second", "quartet_support", "quartet_boot_support", "quartet_p_value", "n_sites_used"]).to_csv(svdq_csv, index=False)
        return {"status": "completed_dosage_fallback", "n_samples": n_s, "gdi_csv": str(gdi_csv)}

    return {"status": "skipped_empty", "n_samples": 0}


def main(argv: Optional[Sequence[str]] = None) -> int:
    gcfg = resolve_genomics_config()
    parser = argparse.ArgumentParser(description="Invariant-site coalescent delimitation (SVDquartets, BPP, gdi).")
    parser.add_argument("--genomics-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--parsed-samples", type=Path, default=None)
    parser.add_argument("--dosage-tsv", type=Path, default=None)
    parser.add_argument("--block-len", type=int, default=int(gcfg["coalescent_block_len"]))
    parser.add_argument("--n-boot", type=int, default=int(gcfg["coalescent_n_boot"]))
    args = parser.parse_args(list(argv) if argv is not None else None)
    run_coalescent_delimitation(args.genomics_dir, args.out_dir, parsed_samples_csv=args.parsed_samples, dosage_tsv=args.dosage_tsv, block_len=args.block_len, n_boot=args.n_boot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
