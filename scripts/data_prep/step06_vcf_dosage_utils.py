#!/usr/bin/env python3
"""
Module: scripts/data_prep/step06_vcf_dosage_utils.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  VCF genotype/read-depth extraction, precomputed dosage-table parsing,
  companion VCF-to-dosage depth matrix alignment, and organellar FASTA/CSV readers.
"""

from __future__ import annotations

import gzip, re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from scripts.data_prep.step06_sample_metadata_utils import clean_sample_id


def _gt_to_dosage(gt_token: str) -> float:
    """Convert diploid or polyploid GT (`0/0`, `0/1`, `1/1`, `0/0/1/1`, `1|0`) to 0..1 dosage."""
    if not gt_token or "." in gt_token:
        return np.nan
    alleles = re.split(r"[/|]", gt_token)
    if not alleles:
        return np.nan
    try:
        vals = [int(a) for a in alleles]
    except ValueError:
        return np.nan
    return float(sum(1 for a in vals if a > 0)) / float(len(vals)) if vals else np.nan


def _safe_nonneg_int(tok: str) -> Optional[int]:
    if not tok or tok in (".", "./.", ".|."):
        return None
    try:
        v = int(float(tok))
        return v if v >= 0 else None
    except ValueError:
        return None


def _extract_cell_depths(fields: Sequence[str], fmt_idx: Mapping[str, int], gt_dosage: float) -> Tuple[float, float]:
    """Extract `(n_ref, n_alt)` supporting FreeBayes RO/AO, GATK `ref,alt` AD, and BBMap single-int AD."""
    ro_i, ao_i = fmt_idx.get("RO", -1), fmt_idx.get("AO", -1)
    if 0 <= ro_i < len(fields) and 0 <= ao_i < len(fields):
        ro_v = _safe_nonneg_int(fields[ro_i])
        ao_clean = [x for t in fields[ao_i].split(",") if (x := _safe_nonneg_int(t)) is not None]
        if ro_v is not None and ao_clean:
            return float(ro_v), float(sum(ao_clean))

    ad_i, dp_i = fmt_idx.get("AD", -1), fmt_idx.get("DP", -1)
    dp_v = _safe_nonneg_int(fields[dp_i]) if (0 <= dp_i < len(fields)) else None
    if 0 <= ad_i < len(fields) and (ad_tok := fields[ad_i].strip()) and ad_tok != ".":
        if "," in ad_tok:
            ad_parts = [_safe_nonneg_int(t) for t in ad_tok.split(",")]
            if ad_parts[0] is not None and any(x is not None for x in ad_parts[1:]):
                return float(ad_parts[0]), float(sum(x for x in ad_parts[1:] if x is not None))
        elif (n_alt_i := _safe_nonneg_int(ad_tok)) is not None:
            # BBMap v38.69 records a single integer in AD equal to alt allele read count
            return (float(max(0, dp_v - n_alt_i)) if dp_v is not None else 0.0), float(n_alt_i)

    if dp_v is not None and dp_v > 0 and np.isfinite(gt_dosage):
        n_alt_f = float(round(dp_v * float(gt_dosage)))
        return float(max(0.0, float(dp_v) - n_alt_f)), n_alt_f
    return np.nan, np.nan


def _parse_vcf_with_depths(vcf_path: Path) -> Tuple[List[str], List[str], np.ndarray, np.ndarray, np.ndarray]:
    """Stream VCF / `.vcf.gz` into `(sample_ids, locus_ids, dosage_mat, ref_depth_mat, alt_depth_mat)`."""
    opener = gzip.open if vcf_path.name.lower().endswith(".gz") else open
    sample_ids: List[str] = []
    locus_ids: List[str] = []
    rows_loci: List[np.ndarray] = []
    rows_ref: List[np.ndarray] = []
    rows_alt: List[np.ndarray] = []

    with opener(vcf_path, "rt", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_ids = [clean_sample_id(s) for s in line.rstrip("\r\n").split("\t")[9:]]
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 10 or not sample_ids:
                continue
            chrom, pos, vid = parts[0], parts[1], parts[2]
            locus_id = vid if (vid and vid != ".") else f"{chrom}_{pos}"
            if locus_ids and locus_id == locus_ids[-1]:
                locus_id = f"{chrom}_{pos}_{len(locus_ids)}"
            fmt = parts[8].split(":") if len(parts) > 8 else ["GT"]
            fmt_idx = {tag: i for i, tag in enumerate(fmt)}
            gt_idx, dos_idx = fmt_idx.get("GT", 0), fmt_idx.get("DS", -1)

            vec = np.full(len(sample_ids), np.nan, dtype=np.float32)
            rvec = np.full(len(sample_ids), np.nan, dtype=np.float32)
            avec = np.full(len(sample_ids), np.nan, dtype=np.float32)
            for idx, cell in enumerate(parts[9 : 9 + len(sample_ids)]):
                fields = cell.split(":")
                gt_d = _gt_to_dosage(fields[gt_idx] if gt_idx < len(fields) else ".")
                if 0 <= dos_idx < len(fields) and fields[dos_idx] not in ("", "."):
                    try:
                        dv = float(fields[dos_idx])
                        vec[idx] = np.float32(dv if dv <= 1.0 else dv / 2.0)
                    except ValueError:
                        vec[idx] = np.float32(gt_d)
                else:
                    vec[idx] = np.float32(gt_d)
                if np.isfinite(vec[idx]):
                    nr, na = _extract_cell_depths(fields, fmt_idx, float(vec[idx]))
                    rvec[idx], avec[idx] = np.float32(nr), np.float32(na)
            locus_ids.append(locus_id)
            rows_loci.append(vec); rows_ref.append(rvec); rows_alt.append(avec)

    if not sample_ids or not rows_loci:
        empty = np.empty((len(sample_ids), 0), dtype=np.float32)
        return sample_ids, [], empty, empty.copy(), empty.copy()
    return (
        sample_ids, locus_ids,
        np.vstack(rows_loci).T.astype(np.float32),
        np.vstack(rows_ref).T.astype(np.float32),
        np.vstack(rows_alt).T.astype(np.float32),
    )


def _parse_dosage_table(dosage_path: Path) -> Tuple[List[str], List[str], np.ndarray]:
    sep = "\t" if dosage_path.suffix.lower() in (".tsv", ".txt") else ","
    df = pd.read_csv(dosage_path, sep=sep)
    if df.empty or df.shape[1] < 2:
        return [], [], np.empty((0, 0), dtype=np.float32)
    sample_ids = [clean_sample_id(s) for s in df[df.columns[0]].tolist()]
    locus_cols = [str(c) for c in df.columns[1:]]
    mat = df[locus_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    finite = mat[np.isfinite(mat)]
    if finite.size > 0 and float(np.max(finite)) > 1.05:
        max_p = 2.0 if float(np.max(finite)) <= 2.05 else float(np.max(finite))
        mat = mat / np.float32(max_p)
    return sample_ids, locus_cols, mat


def align_vcf_depths_to_dosage_grid(
    d_sids: Sequence[str], d_loci: Sequence[str],
    vcf_sids: Sequence[str], vcf_loci: Sequence[str],
    vcf_ref: np.ndarray, vcf_alt: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Align parsed companion VCF read depths onto precomputed dosage `(d_sids, d_loci)` coordinates."""
    empty = np.empty((0, 0), dtype=np.float32)
    if vcf_ref.size == 0 or vcf_alt.size == 0 or not d_sids or not d_loci:
        return empty, empty.copy()
    s_map = {s: i for i, s in enumerate(vcf_sids)}
    l_map = {l: j for j, l in enumerate(vcf_loci)}
    row_pairs = [(di, s_map[s]) for di, s in enumerate(d_sids) if s in s_map]
    col_pairs = [(dj, l_map[l]) for dj, l in enumerate(d_loci) if l in l_map]
    if not row_pairs or not col_pairs:
        return empty, empty.copy()
    d_ri = np.asarray([p[0] for p in row_pairs], dtype=np.intp)
    v_ri = np.asarray([p[1] for p in row_pairs], dtype=np.intp)
    d_cj = np.asarray([p[0] for p in col_pairs], dtype=np.intp)
    v_cj = np.asarray([p[1] for p in col_pairs], dtype=np.intp)
    out_ref = np.full((len(d_sids), len(d_loci)), np.nan, dtype=np.float32)
    out_alt = np.full((len(d_sids), len(d_loci)), np.nan, dtype=np.float32)
    out_ref[np.ix_(d_ri, d_cj)] = vcf_ref[np.ix_(v_ri, v_cj)]
    out_alt[np.ix_(d_ri, d_cj)] = vcf_alt[np.ix_(v_ri, v_cj)]
    return out_ref, out_alt


def _read_fasta_headers_seqs(path: Path) -> List[Tuple[str, str]]:
    entries: List[Tuple[str, str]] = []
    cur_h: Optional[str] = None
    buf: List[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not (s := line.strip()):
                continue
            if s.startswith(">"):
                if cur_h is not None: entries.append((cur_h, "".join(buf).upper()))
                cur_h, buf = s[1:].split()[0], []
            else:
                buf.append(s)
    if cur_h is not None:
        entries.append((cur_h, "".join(buf).upper()))
    return entries


def _load_organellar_metadata(organellar_files: Sequence[Path]) -> Dict[str, Dict[str, Any]]:
    """Read staged cytonuclear table and skimmed plastid/ITS consensus FASTAs."""
    meta: Dict[str, Dict[str, Any]] = {}
    for p in organellar_files:
        low = p.name.lower()
        if low.startswith("cytonuclear_discordance") and low.endswith((".csv", ".tsv")):
            try:
                df = pd.read_csv(p, sep=("\t" if low.endswith(".tsv") else ","), dtype=str, keep_default_na=False)
            except Exception:
                continue
            for _, row in df.iterrows():
                if not (sid := clean_sample_id(row.get("genetic_sample_id", ""))):
                    continue
                rec = meta.setdefault(sid, {})
                for org_col in ("plastid_clade", "its_ribotype", "public_plastome_match", "public_its_match", "waxy_lineage_match"):
                    if row.get(org_col): rec[org_col] = str(row[org_col]).strip()
                disc = str(row.get("cytonuclear_discordant", "")).strip().lower()
                if disc in ("1", "true", "yes"): rec["cytonuclear_discordant"] = 1
                elif disc in ("0", "false", "no"): rec["cytonuclear_discordant"] = 0
        elif low.endswith((".fasta", ".fa", ".fas")):
            target_col = "plastid_clade" if "plastid" in low else ("its_ribotype" if "its" in low else None)
            if not target_col: continue
            seq_map: Dict[str, str] = {}
            for hdr, seq in _read_fasta_headers_seqs(p):
                sid = clean_sample_id(hdr.split("|")[0])
                clade_tag = next((v.strip() for tok in hdr.split("|")[1:] if "=" in tok and (k := tok.split("=", 1)[0].lower()) in ("clade", "ribotype", target_col) and (v := tok.split("=", 1)[1])), None)
                if sid:
                    rec = meta.setdefault(sid, {})
                    if not rec.get(target_col):
                        if clade_tag:
                            rec[target_col] = clade_tag
                        else:
                            if seq not in seq_map:
                                seq_map[seq] = f"{'cp' if target_col == 'plastid_clade' else 'ITS'}_H{len(seq_map) + 1}"
                            rec[target_col] = seq_map[seq]
    return meta

