#!/usr/bin/env python3
"""
Module: scripts/analysis/step06_coalescent_alignment_utils.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Alignment parsers (PHYLIP/NEXUS/FASTA), replicate collapsing, SVDquartets
  16x16 tensor bootstrapper, BPP MSci control/PHYLIP exporter, and skimmed-FASTA
  substitution model selection helpers for Step 06 coalescent delimitation.
"""

from __future__ import annotations

import importlib
import math
import os
import re
import sqlite3
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from scripts.analysis.step06_coalescent_models import (
    _resolve_max_workers,
    compute_model_corrected_divergence_blocks,
    select_nucleotide_substitution_model,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]

_IUPAC_ALLELES: Dict[str, Tuple[str, str]] = {
    "A": ("A", "A"), "C": ("C", "C"), "G": ("G", "G"), "T": ("T", "T"),
    "R": ("A", "G"), "Y": ("C", "T"), "S": ("C", "G"),
    "W": ("A", "T"), "K": ("G", "T"), "M": ("A", "C"),
}
_BASE_TO_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _clean_id(raw: str) -> str:
    return importlib.import_module("scripts.data_prep.06_ingest_genetic_data").clean_sample_id(raw)


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
            if not (s := line.strip()):
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
                if len(parts := s.split(None, 1)) == 2:
                    sid = _clean_id(parts[0].strip("'\""))
                    if sid not in seqs:
                        order.append(sid)
                    seqs[sid].append(re.sub(r"\s+", "", parts[1]).upper())
        return order, ["".join(seqs[k]) for k in order]

    lines = [ln.rstrip("\r\n") for ln in text.splitlines() if ln.strip()]
    if not lines:
        return [], []
    hdr = lines[0].strip().split()
    n_taxa = int(hdr[0]) if hdr and hdr[0].isdigit() else 0
    body = lines[1:]
    if n_taxa > 0 and len(body) >= n_taxa and all(len(ln.split(None, 1)) == 2 for ln in body[:n_taxa]):
        for ln in body[:n_taxa]:
            nm, sq = ln.split(None, 1)
            sid = _clean_id(nm)
            order.append(sid)
            seqs[sid].append(re.sub(r"\s+", "", sq).upper())
        pos_i = 0
        for ln in body[n_taxa:]:
            parts = ln.split(None, 1)
            if len(parts) == 2 and _clean_id(parts[0]) in seqs:
                seqs[_clean_id(parts[0])].append(re.sub(r"\s+", "", parts[1]).upper())
            else:
                seqs[order[pos_i % n_taxa]].append(re.sub(r"\s+", "", ln).upper())
                pos_i += 1
        return order, ["".join(seqs[k]) for k in order]

    for ln in body:
        if len(parts := ln.split(None, 1)) == 2:
            sid = _clean_id(parts[0])
            if sid not in seqs:
                order.append(sid)
            seqs[sid].append(re.sub(r"\s+", "", parts[1]).upper())
    return order, ["".join(seqs[k]) for k in order]


def _encode_haplotype_alleles(seqs: Sequence[str]) -> Tuple[np.ndarray, np.ndarray]:
    """Encode sequences to diploid base pair arrays `a1, a2` of shape `(n_samples, n_sites)` (`-1` = missing)."""
    n_samp, n_sites = len(seqs), (min(len(s) for s in seqs) if seqs else 0)
    a1 = np.full((n_samp, n_sites), -1, dtype=np.int8)
    a2 = np.full((n_samp, n_sites), -1, dtype=np.int8)
    for i, sq in enumerate(seqs):
        arr = np.frombuffer(sq[:n_sites].encode("ascii", errors="replace"), dtype=np.uint8)
        for ch, (b1, b2) in _IUPAC_ALLELES.items():
            mask = arr == ord(ch)
            if np.any(mask):
                a1[i, mask], a2[i, mask] = _BASE_TO_IDX[b1], _BASE_TO_IDX[b2]
    return a1, a2


def _pairwise_site_divergence_blocks(
    a1: np.ndarray, a2: np.ndarray, block_slices: Sequence[Tuple[int, int]],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute per-block `(n_samp, n_samp)` pairwise allele mismatch counts and callable site counts."""
    return compute_model_corrected_divergence_blocks(a1, a2, block_slices, model_spec=None)


def compute_svdquartets_support(
    a1: np.ndarray, a2: np.ndarray, taxa_labels: Sequence[str], sample_ids: Sequence[str],
    max_quartets: int = 120, seed: int = 202611, n_boot: int = 200,
) -> pd.DataFrame:
    """Compute Chifman & Kubatko (2014) 16x16 SVDquartets scores and block-bootstrap p-values."""
    _ = (a2, sample_ids)
    empty_cols = ["quartet_id", "best_split", "svd_score_best", "svd_score_second", "quartet_support", "quartet_boot_support", "quartet_p_value", "n_sites_used"]
    taxon_groups: Dict[str, List[int]] = defaultdict(list)
    for idx, tx in enumerate(taxa_labels):
        taxon_groups[tx].append(idx)
    uniq_taxa = sorted(taxon_groups.keys())
    if len(uniq_taxa) < 4 or a1.shape[1] < 20:
        return pd.DataFrame(columns=empty_cols)

    rng = np.random.default_rng(seed)
    all_q = list(combinations(uniq_taxa, 4))
    if len(all_q) > max_quartets:
        sel = rng.choice(len(all_q), size=max_quartets, replace=False)
        all_q = [all_q[int(k)] for k in sorted(sel)]

    rows: List[Dict[str, Any]] = []
    rep_idx = {tx: max(idxs, key=lambda i: int(np.sum(a1[i] >= 0))) for tx, idxs in taxon_groups.items()}
    n_sites = a1.shape[1]

    def _svd_residual(mats: np.ndarray) -> np.ndarray:
        sv = np.linalg.svd(mats, compute_uv=False)
        return np.sqrt(np.sum(sv[..., 10:] ** 2, axis=-1))

    q_jobs: List[Tuple[str, str, str, str, int, int, int, int, np.ndarray, int, np.ndarray]] = []
    for t1, t2, t3, t4 in all_q:
        i1, i2, i3, i4 = rep_idx[t1], rep_idx[t2], rep_idx[t3], rep_idx[t4]
        valid = (a1[i1] >= 0) & (a1[i2] >= 0) & (a1[i3] >= 0) & (a1[i4] >= 0)
        if (nv := int(np.sum(valid))) < 20:
            continue
        n_b = max(4, min(40, nv // (2 if nv < 40 else 10)))
        draws = rng.integers(0, n_b, size=(max(1, n_boot), n_b))
        q_jobs.append((t1, t2, t3, t4, i1, i2, i3, i4, valid, n_b, draws))

    def _eval_quartet_svd(job: Tuple[str, str, str, str, int, int, int, int, np.ndarray, int, np.ndarray]) -> Dict[str, Any]:
        t1, t2, t3, t4, i1, i2, i3, i4, valid, n_b, draws = job
        nv = int(np.sum(valid))
        b1, b2, b3, b4 = a1[i1, valid], a1[i2, valid], a1[i3, valid], a1[i4, valid]
        blk_splits = np.array_split(np.arange(nv), n_b)
        blk_tensors = np.zeros((n_b, 256), dtype=np.float64)
        idx_256 = b1.astype(np.int32) * 64 + b2.astype(np.int32) * 16 + b3.astype(np.int32) * 4 + b4.astype(np.int32)
        for bi, b_idx in enumerate(blk_splits):
            np.add.at(blk_tensors[bi], idx_256[b_idx], 1.0)

        obs_t4 = np.sum(blk_tensors, axis=0).reshape(4, 4, 4, 4) / float(nv)
        obs_mats = np.stack([obs_t4.reshape(16, 16), obs_t4.transpose(0, 2, 1, 3).reshape(16, 16), obs_t4.transpose(0, 3, 1, 2).reshape(16, 16)], axis=0)
        obs_scores = _svd_residual(obs_mats)
        split_labels = [f"{t1},{t2}|{t3},{t4}", f"{t1},{t3}|{t2},{t4}", f"{t1},{t4}|{t2},{t3}"]
        order = np.argsort(obs_scores, kind="mergesort")
        best_idx = int(order[0])
        s_best, s_sec, s_trd = float(obs_scores[order[0]]), float(obs_scores[order[1]]), float(obs_scores[order[2]])
        scale = max(1e-6, float(np.std([s_best, s_sec, s_trd])))
        w = np.exp(-np.array([s_best, s_sec, s_trd]) / scale)
        support = float(w[0] / np.sum(w))

        bw = np.zeros((draws.shape[0], n_b), dtype=np.float64)
        for bi in range(n_b):
            bw[:, bi] = np.sum(draws == bi, axis=1)
        boot_raw = (bw @ blk_tensors).reshape(-1, 4, 4, 4, 4)
        boot_t4 = boot_raw / np.maximum(1.0, np.sum(boot_raw, axis=(1, 2, 3, 4), keepdims=True))
        b_s0 = _svd_residual(boot_t4.reshape(-1, 16, 16))
        b_s1 = _svd_residual(boot_t4.transpose(0, 1, 3, 2, 4).reshape(-1, 16, 16))
        b_s2 = _svd_residual(boot_t4.transpose(0, 1, 4, 2, 3).reshape(-1, 16, 16))
        boot_win = np.argmin(np.stack([b_s0, b_s1, b_s2], axis=1), axis=1)
        n_match = int(np.sum(boot_win == best_idx))
        boot_supp = float(n_match / float(len(boot_win)))
        p_val = float((len(boot_win) - n_match + 1) / float(len(boot_win) + 1))

        return {
            "quartet_id": f"{t1};{t2};{t3};{t4}", "best_split": split_labels[best_idx],
            "svd_score_best": round(s_best, 6), "svd_score_second": round(s_sec, 6),
            "quartet_support": round(support, 4), "quartet_boot_support": round(boot_supp, 4),
            "quartet_p_value": round(p_val, 4), "n_sites_used": min(nv, n_sites),
        }

    nw_q = min(_resolve_max_workers(), len(q_jobs)) if len(q_jobs) >= 4 else 1
    if nw_q > 1:
        with ThreadPoolExecutor(max_workers=nw_q) as pool:
            rows = list(pool.map(_eval_quartet_svd, q_jobs))
    else:
        rows = [_eval_quartet_svd(j) for j in q_jobs]
    return pd.DataFrame(rows, columns=empty_cols) if not rows else pd.DataFrame(rows)


def export_bpp_msci_files(
    sample_ids: Sequence[str], taxa_labels: Sequence[str], seqs: Sequence[str],
    block_slices: Sequence[Tuple[int, int]], out_dir: Path, max_loci: int = 200, substitution_model: str = "JC69",
) -> Tuple[Path, Path]:
    """Write BPP-formatted multi-locus PHYLIP and `bpp_msci.ctl` into run-scoped `out_dir`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    loci_phy, ctl_path, imap_path = out_dir / "bpp_msci_loci.phy", out_dir / "bpp_msci.ctl", out_dir / "bpp_imap.txt"
    sel_blocks = list(block_slices)[:max_loci]

    with open(loci_phy, "wt", encoding="utf-8") as fh:
        for st, en in sel_blocks:
            keep_rows = [(sid, sq[st:en]) for sid, sq in zip(sample_ids, seqs) if any(ch in "ACGT" for ch in sq[st:en])]
            if len(keep_rows) < 4:
                continue
            fh.write(f"{len(keep_rows)} {en - st}\n")
            for sid, ss in keep_rows:
                fh.write(f"^{sid:<18s} {re.sub(r'[^ACGTURYKSWMN?-]', 'N', ss)}\n")
            fh.write("\n")

    uniq_taxa = sorted(set(taxa_labels))
    with open(imap_path, "wt", encoding="utf-8") as ih:
        for sid, tx in zip(sample_ids, taxa_labels):
            ih.write(f"{sid}\t{tx.replace(' ', '_')}\n")

    base_sub = str(substitution_model or "JC69").split("+", 1)[0].upper()
    bpp_model = "HKY" if base_sub in ("HKY85", "TN93", "F81", "K80") else "JC69"
    gamma_line = "alpha = 0.5\n" if "+G4" in str(substitution_model or "").upper() else ""
    species_tokens = " ".join(tx.replace(" ", "_") for tx in uniq_taxa)
    ctl_txt = (
        f"* BIC selected substitution model: {substitution_model}\nseed = 202611\nseqfile = {loci_phy.name}\n"
        f"Imapfile = {imap_path.name}\noutfile = bpp_msci_out.txt\nmcmcfile = bpp_msci_mcmc.txt\nmodel = {bpp_model}\n"
        f"{gamma_line}speciesdelimitation = 0\nspeciestree = 0\nspecies&tree = {len(uniq_taxa)} {species_tokens}\n"
        f"    {' '.join(['1'] * len(uniq_taxa))}\nthetaprior = 3 0.004 e\ntauprior = 3 0.004\nphiprior = 1 1\n"
        f"nloci = {len(sel_blocks)}\nburnin = 4000\nnsample = 20000\nsampfreq = 2\n"
    )
    ctl_path.write_text(ctl_txt, encoding="utf-8")
    return loci_phy, ctl_path


def _load_sample_metadata(
    genomics_dir: Path, parsed_samples_csv: Optional[Path],
) -> Tuple[Optional[set[str]], Dict[str, str], Dict[str, str], Dict[str, float], Dict[str, float]]:
    sample_taxon_map: Dict[str, str] = {}
    rep_pair_map: Dict[str, str] = {}
    hq_reads_map: Dict[str, float] = {}
    miss_rate_map: Dict[str, float] = {}
    pass_set: Optional[set[str]] = None

    if parsed_samples_csv and parsed_samples_csv.exists():
        sdf = pd.read_csv(parsed_samples_csv, dtype=str, keep_default_na=False)
        if not sdf.empty:
            pass_set = {str(r["genetic_sample_id"]) for _, r in sdf.iterrows() if str(r.get("geno_quality_pass", "1")) in ("1", "True", "true")}
            for _, r in sdf.iterrows():
                sid = str(r["genetic_sample_id"])
                tx = str(r.get("species_standardized") or r.get("taxon_raw") or "")
                if str(r.get("is_putative_novel_or_hybrid_query", "0")) == "1":
                    pop = str(r.get("population_id") or "")
                    tx = f"{tx}_{pop}" if (tx and pop) else (tx or sid)
                if tx:
                    sample_taxon_map[sid] = tx
                if (rp := str(r.get("replicate_pair_id") or "").strip()) and rp.upper() not in {"NA", "NONE", "NAN", "NULL"}:
                    rep_pair_map[sid] = rp
                for raw_h, dst in ((r.get("hq_reads"), hq_reads_map), (r.get("missing_genotype_rate"), miss_rate_map)):
                    try:
                        if raw_h not in (None, "") and math.isfinite(fv := float(raw_h)):
                            dst[sid] = fv
                    except (TypeError, ValueError):
                        pass

    db_cands = [Path(env_p) for env_k in ("PLANTMATICS_DB_PATH", "PACKERA_DB_PATH") if (env_p := os.environ.get(env_k, "").strip())] + [genomics_dir.parent / "plantmatics.sqlite", PROJECT_ROOT / "data" / "plantmatics.sqlite"]
    for db_p in db_cands:
        if not db_p.is_file():
            continue
        try:
            conn = sqlite3.connect(f"file:{db_p}?mode=ro", uri=True)
            try:
                for sid_r, sp_r, rp_r, hq_r, mr_r, qp_r in conn.execute("SELECT genetic_sample_id, species_standardized, replicate_pair_id, hq_reads, missing_genotype_rate, geno_quality_pass FROM genomic_samples;").fetchall():
                    if not (sid_s := str(sid_r or "").strip()):
                        continue
                    if sid_s not in sample_taxon_map and sp_r: sample_taxon_map[sid_s] = str(sp_r)
                    rp_s = str(rp_r or "").strip()
                    if sid_s not in rep_pair_map and rp_s and rp_s.upper() not in {"NA", "NONE", "NAN", "NULL"}: rep_pair_map[sid_s] = rp_s
                    if sid_s not in hq_reads_map and hq_r is not None: hq_reads_map[sid_s] = float(hq_r)
                    if sid_s not in miss_rate_map and mr_r is not None: miss_rate_map[sid_s] = float(mr_r)
                    if pass_set is None and qp_r is not None and int(qp_r) == 1:
                        pass_set = pass_set or set()
                        pass_set.add(sid_s)
            finally:
                conn.close()
            break
        except sqlite3.Error:
            continue
    return pass_set, sample_taxon_map, rep_pair_map, hq_reads_map, miss_rate_map


def _group_technical_replicates(sids: Sequence[str], taxa_labels: Sequence[str], rep_pair_map: Mapping[str, str]) -> Dict[Any, List[int]]:
    groups: Dict[Any, List[int]] = defaultdict(list)
    for idx, (sid, tx) in enumerate(zip(sids, taxa_labels)):
        rp = str(rep_pair_map.get(sid) or "").strip()
        groups[("pair", rp, tx) if rp else ("uniq", idx)].append(idx)
    return groups


def _collapse_alignment_replicates(
    sids: Sequence[str], seqs: Sequence[str], taxa_labels: Sequence[str],
    rep_pair_map: Mapping[str, str], hq_reads_map: Mapping[str, float], miss_rate_map: Mapping[str, float],
) -> Tuple[List[str], List[str], List[str], Dict[str, List[str]]]:
    groups = _group_technical_replicates(sids, taxa_labels, rep_pair_map)
    out_sids, out_seqs, out_taxa, group_members = [], [], [], {}
    for idxs in groups.values():
        ranked = sorted(idxs, key=lambda i: (-sum(1 for ch in seqs[i] if ch not in "N?-"), -float(hq_reads_map.get(sids[i], 0.0)), float(miss_rate_map.get(sids[i], 1.0)), sids[i]))
        p_i, p_sid = ranked[0], sids[ranked[0]]
        merged = seqs[p_i]
        if len(ranked) > 1:
            buf = list(merged)
            for s_i in ranked[1:]:
                sib = seqs[s_i]
                for pos in range(min(len(buf), len(sib))):
                    if buf[pos] in "N?-" and sib[pos] not in "N?-":
                        buf[pos] = sib[pos]
            merged = "".join(buf)
        out_sids.append(p_sid); out_seqs.append(merged); out_taxa.append(taxa_labels[p_i])
        group_members[p_sid] = [sids[i] for i in idxs]
    return out_sids, out_seqs, out_taxa, group_members


def _collapse_dosage_replicates(
    sids: Sequence[str], mat: np.ndarray, taxa_labels: Sequence[str],
    rep_pair_map: Mapping[str, str], hq_reads_map: Mapping[str, float], miss_rate_map: Mapping[str, float],
) -> Tuple[List[str], np.ndarray, List[str], Dict[str, List[str]]]:
    groups = _group_technical_replicates(sids, taxa_labels, rep_pair_map)
    out_sids, out_rows, out_taxa, group_members = [], [], [], {}
    for idxs in groups.values():
        ranked = sorted(idxs, key=lambda i: (-int(np.sum(np.isfinite(mat[i]))), -float(hq_reads_map.get(sids[i], 0.0)), float(miss_rate_map.get(sids[i], 1.0)), sids[i]))
        p_i, p_sid = ranked[0], sids[ranked[0]]
        if len(ranked) == 1:
            merged_row = mat[p_i]
        else:
            sub = mat[ranked]
            ok = np.isfinite(sub)
            cnt = np.sum(ok, axis=0)
            merged_row = np.where(cnt > 0, np.sum(np.where(ok, sub, 0.0), axis=0) / np.maximum(cnt, 1), np.nan)
        out_sids.append(p_sid); out_rows.append(merged_row); out_taxa.append(taxa_labels[p_i])
        group_members[p_sid] = [sids[i] for i in idxs]
    return out_sids, (np.vstack(out_rows) if out_rows else mat[:0]), out_taxa, group_members


def _expand_gdi_samples(gdi_s_df: pd.DataFrame, group_members: Mapping[str, List[str]], orig_sids: Sequence[str]) -> pd.DataFrame:
    if gdi_s_df.empty or not group_members:
        return gdi_s_df
    by_rep = {str(r["genetic_sample_id"]): r.to_dict() for _, r in gdi_s_df.iterrows()}
    orig_to_rep = {m: rep for rep, members in group_members.items() for m in members}
    rows = []
    for sid in orig_sids:
        if (rep := orig_to_rep.get(sid, sid)) in by_rep:
            rec = dict(by_rep[rep]); rec["genetic_sample_id"] = sid; rows.append(rec)
    return pd.DataFrame(rows, columns=gdi_s_df.columns) if rows else gdi_s_df


def _select_models_for_skimmed_fastas(genomics_dir: Path, exclude_paths: Sequence[Path]) -> List[Dict[str, Any]]:
    """Run BIC substitution model selection across any recursive `skimmed_*.fasta` alignments."""
    excl_resolved = {p.resolve() for p in exclude_paths}
    candidates = sorted(p for p in genomics_dir.rglob("skimmed_*.fasta") if p.is_file() and not any(part.startswith(".") for part in p.parts) and p.resolve() not in excl_resolved)
    if not candidates:
        return []

    def _fit_skimmed(fa_path: Path) -> List[Dict[str, Any]]:
        sk_sids, sk_seqs = parse_invariant_alignment(fa_path)
        if len(sk_sids) < 3 or not sk_seqs:
            return []
        return list(select_nucleotide_substitution_model(dict(zip(sk_sids, sk_seqs)), sk_sids, alignment_name=fa_path.name).get("table_rows") or [])

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

