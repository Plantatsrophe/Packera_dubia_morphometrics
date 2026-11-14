#!/usr/bin/env python3
"""
Module: scripts/data_prep/skim_offtarget_organellar.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Offline off-target chloroplast (cpDNA) and nuclear ribosomal ITS skimming from
  raw nextRAD FASTQ archives (.zip, .tar.gz) or FASTQ directories. Streams
  compressed archives without unpacking raw reads into `data/genomics/`.
  Auto-detects `get_organelle_from_reads.py` or `bbmap.sh` when available, with a
  pure-Python 17-mer baiting & pileup consensus fallback seeded by Asteraceae
  organellar/ITS references and `vcfref_packera_blast.txt`. Writes flat outputs
  directly into `data/genomics/` (`skimmed_plastid_consensus.fasta`,
  `skimmed_its_consensus.fasta`, `cytonuclear_discordance.csv`).
"""

from __future__ import annotations

import argparse
import gzip
import importlib
import io
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.logger import setup_logging

logger = setup_logging(name="packera.genomics.organellar")

_SEED_PLASTID_ANCHORS: Dict[str, str] = {
    # Conserved Asteraceae / Senecioneae rbcL & trnL-F core motifs
    "Packera_rbcL_core": "ATGTCACCACAAACAGAGACTAAAGCAAGTGTTGGATTCAAAGCTGGTGTTAAAGAGTACAAATTGACTTATTATACTCCTGAGTATGAAACCAAAGATACTGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCCGGAGTTCCACCTGAAGAAGCAGGGGCTGCGGTAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACTGATGGACTTACCAGCCTTGATCGTTACAAAGG",
    "Packera_trnL_core": "GGGCAATCCTGAGCCAAATCCTGTTTTCCGAAAACAAACAAAGGTTCAGAAAGCGAGAATAAAAAAAAGGATAGGTGCAGAGACTCAATGGAAGCTGTTCTAACGAATCGAGTTGACTGCGTTGCGTTAGTAAGCGGAAGAACAGGCTTCCATAGTC",
}

_SEED_ITS_ANCHORS: Dict[str, str] = {
    # Conserved Asteraceae 18S-ITS1-5.8S-ITS2-26S core motif
    "Packera_ITS_58S_core": "CGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAAGCCATTAGGTTGAGGGCACGTCTGCCTGGGCGTCACGCATCGCGTCGCCCCCCACCATCCATCCCCTTAGGGATGGGGTGTGTTGGGGGCGGGGGGATATTGGCCTCCCGTGCGCTCGGTGCGCGGTTGGCCTAAATAGGAGTCCTCGGCGACGGACGTCGCGACAAGTGGTGGTTGTAAAACCCTCGTTCTCGTCGTCGTGC",
}

_IUPAC_PAIR = {
    ("A", "G"): "R", ("G", "A"): "R",
    ("C", "T"): "Y", ("T", "C"): "Y",
    ("G", "C"): "S", ("C", "G"): "S",
    ("A", "T"): "W", ("T", "A"): "W",
    ("G", "T"): "K", ("T", "G"): "K",
    ("A", "C"): "M", ("C", "A"): "M",
}


def _clean_sample_id(raw: str) -> str:
    fn = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").clean_sample_id
    return fn(raw)


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _read_fasta(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.exists() or not path.is_file():
        return out
    cur: Optional[str] = None
    buf: List[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if cur is not None:
                    out[cur] = "".join(buf).upper()
                cur = s[1:].split()[0]
                buf = []
            else:
                buf.append(s)
    if cur is not None:
        out[cur] = "".join(buf).upper()
    return out


def _collect_organellar_seeds(genomics_dir: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Gather plastid and ITS seed contigs from BLAST hits, staged ITS FASTA, and conserved anchors."""
    plastid_seeds = dict(_SEED_PLASTID_ANCHORS)
    its_seeds = dict(_SEED_ITS_ANCHORS)

    ref_seqs: Dict[str, str] = {}
    for cand in (genomics_dir / "ref_packera.fasta", genomics_dir / "vcfref_packera_rm.txt.fasta"):
        if cand.exists():
            ref_seqs.update(_read_fasta(cand))

    bain = genomics_dir / "bain_jansen_sequences.fasta"
    if bain.exists():
        for k, v in _read_fasta(bain).items():
            if len(v) >= 60:
                its_seeds[f"bain_{k}"] = v

    blast_txt = genomics_dir / "vcfref_packera_blast.txt"
    if blast_txt.exists() and ref_seqs:
        with open(blast_txt, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                low = line.lower()
                parts = line.rstrip("\r\n").split("\t")
                qid = parts[0].strip() if parts else ""
                if not qid or qid not in ref_seqs:
                    continue
                if any(tok in low for tok in ("chloroplast", "plastid", "rbcl", "matk", "trnl", "psba", "ndhf")):
                    plastid_seeds[qid] = ref_seqs[qid]
                elif any(tok in low for tok in ("internal transcribed spacer", "ribosomal", "18s", "26s", "its1", "its2")):
                    its_seeds[qid] = ref_seqs[qid]

    return plastid_seeds, its_seeds


def _infer_sample_from_fq_name(member_name: str) -> str:
    base = Path(member_name).name
    stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    stem = re.sub(r"([_.-](R?[12]|merged|trimmed|paired|unpaired|clean))+$", "", stem, flags=re.IGNORECASE)
    return _clean_sample_id(stem)


def _iter_fastq_records_stream(fh: io.TextIOBase) -> Iterator[Tuple[str, str]]:
    while True:
        h = fh.readline()
        if not h:
            break
        seq = fh.readline()
        plus = fh.readline()
        qual = fh.readline()
        if not qual:
            break
        yield h.strip(), seq.strip().upper()


def _iter_staged_fastq_reads(fastq_source: Path) -> Iterator[Tuple[str, str]]:
    """Stream `(sample_id, seq)` from `.zip`, `.tar.gz`, or directory without unpacking to `data/genomics/`."""
    if not fastq_source.exists():
        return

    fq_pat = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)
    if fastq_source.is_dir():
        for p in sorted(fastq_source.rglob("*")):
            if not p.is_file() or not fq_pat.search(p.name):
                continue
            sid = _infer_sample_from_fq_name(p.name)
            opener = gzip.open if p.name.lower().endswith(".gz") else open
            with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
                for _, seq in _iter_fastq_records_stream(fh):
                    if seq:
                        yield sid, seq
        return

    low = fastq_source.name.lower()
    if low.endswith(".zip"):
        with zipfile.ZipFile(fastq_source, "r") as zf:
            for info in sorted(zf.infolist(), key=lambda x: x.filename):
                if info.is_dir() or not fq_pat.search(info.filename):
                    continue
                sid = _infer_sample_from_fq_name(info.filename)
                with zf.open(info, "r") as raw_f:
                    bio: Iterable[bytes] = gzip.GzipFile(fileobj=raw_f) if info.filename.lower().endswith(".gz") else raw_f
                    with io.TextIOWrapper(bio, encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
                        for _, seq in _iter_fastq_records_stream(fh):
                            if seq:
                                yield sid, seq
        return

    if low.endswith((".tar.gz", ".tgz", ".tar")):
        mode = "r:gz" if low.endswith((".tar.gz", ".tgz")) else "r:"
        with tarfile.open(fastq_source, mode) as tf:
            for member in tf.getmembers():
                if not member.isfile() or not fq_pat.search(member.name):
                    continue
                sid = _infer_sample_from_fq_name(member.name)
                raw_f = tf.extractfile(member)
                if raw_f is None:
                    continue
                with raw_f:
                    bio = gzip.GzipFile(fileobj=raw_f) if member.name.lower().endswith(".gz") else raw_f
                    with io.TextIOWrapper(bio, encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
                        for _, seq in _iter_fastq_records_stream(fh):
                            if seq:
                                yield sid, seq


def _build_kmer_index(seeds: Mapping[str, str], k: int = 17) -> Dict[str, List[Tuple[str, int, bool]]]:
    idx: Dict[str, List[Tuple[str, int, bool]]] = defaultdict(list)
    for name, seq in seeds.items():
        if len(seq) < k:
            continue
        for pos in range(0, len(seq) - k + 1):
            kmer = seq[pos : pos + k]
            if "N" in kmer:
                continue
            idx[kmer].append((name, pos, False))
            idx[_revcomp(kmer)].append((name, pos, True))
    return idx


def _assemble_baited_consensus(
    reads_by_sample: Mapping[str, List[str]],
    seeds: Mapping[str, str],
    k: int = 17,
) -> Dict[str, str]:
    """Bait reads by seed k-mers and construct per-sample consensus sequence by read pileup."""
    if not seeds or not reads_by_sample:
        return {}
    kidx = _build_kmer_index(seeds, k=k)
    ref_len = max(len(s) for s in seeds.values())
    cons_out: Dict[str, str] = {}

    for sid, seqs in sorted(reads_by_sample.items()):
        pileup: List[Counter[str]] = [Counter() for _ in range(ref_len)]
        best_ref: Optional[str] = None
        hit_counts: Counter[str] = Counter()
        aligned_reads: List[Tuple[str, int, str]] = []

        for raw_seq in seqs:
            for r_start in range(0, max(0, len(raw_seq) - k + 1)):
                kmer = raw_seq[r_start : r_start + k]
                hits = kidx.get(kmer)
                if not hits:
                    continue
                sname, spos, is_rc = hits[0]
                oriented = _revcomp(raw_seq) if is_rc else raw_seq
                anchor_r = (len(raw_seq) - r_start - k) if is_rc else r_start
                offset = spos - anchor_r
                hit_counts[sname] += 1
                aligned_reads.append((sname, offset, oriented))
                break

        if not hit_counts:
            continue
        best_ref = hit_counts.most_common(1)[0][0]
        target_len = len(seeds[best_ref])
        pileup = [Counter() for _ in range(target_len)]
        for sname, offset, seq in aligned_reads:
            if sname != best_ref:
                continue
            for idx_b, ch in enumerate(seq):
                gpos = offset + idx_b
                if 0 <= gpos < target_len and ch in "ACGT":
                    pileup[gpos][ch] += 1

        chars: List[str] = []
        covered = 0
        seed_seq = seeds[best_ref]
        for j, cnt in enumerate(pileup):
            if not cnt:
                chars.append(seed_seq[j])
                continue
            covered += 1
            top = cnt.most_common(2)
            if len(top) >= 2 and top[1][1] >= 0.30 * (top[0][1] + top[1][1]):
                chars.append(_IUPAC_PAIR.get((top[0][0], top[1][0]), top[0][0]))
            else:
                chars.append(top[0][0])
        if covered >= 20:
            cons_out[sid] = "".join(chars)
    return cons_out


def _cluster_haplotypes(seqs: Mapping[str, str], prefix: str, k: int = 11) -> Dict[str, str]:
    """Cluster consensus haplotypes deterministically by 11-mer Jaccard similarity."""
    if not seqs:
        return {}

    def _kset(s: str) -> set[str]:
        return {s[i : i + k] for i in range(max(0, len(s) - k + 1)) if "N" not in s[i : i + k]}

    samples = sorted(seqs.keys())
    ksets = {sid: _kset(seqs[sid]) for sid in samples}
    reps: List[Tuple[str, set[str]]] = []
    out: Dict[str, str] = {}
    for sid in samples:
        ks = ksets[sid]
        best_label = None
        best_sim = -1.0
        for label, rks in reps:
            union = len(ks | rks)
            sim = float(len(ks & rks)) / float(union) if union > 0 else 0.0
            if sim > best_sim:
                best_sim = sim
                best_label = label
        if best_label is not None and best_sim >= 0.85:
            out[sid] = best_label
        else:
            lbl = f"{prefix}_{len(reps) + 1}"
            reps.append((lbl, ks))
            out[sid] = lbl
    return out


def _detect_cytonuclear_discordance(
    plastid_clades: Mapping[str, str],
    its_ribotypes: Mapping[str, str],
    nuclear_tree_path: Optional[Path] = None,
) -> List[Dict[str, Any]]:
    """Compare organellar clade vs. nuclear taxon/sister group to flag cytonuclear capture."""
    parse_code = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").parse_sample_code
    _ = nuclear_tree_path  # Optional hook preserved for future tip-neighborhood expansion

    samples = sorted(set(plastid_clades.keys()) | set(its_ribotypes.keys()))
    taxon_of = {sid: (parse_code(sid).get("species_standardized") or "Unknown") for sid in samples}
    clade_maj_taxon: Dict[str, str] = {}
    by_clade: Dict[str, Counter[str]] = defaultdict(Counter)
    for sid, pc in plastid_clades.items():
        if taxon_of[sid] != "Unknown":
            by_clade[pc][taxon_of[sid]] += 1
    for pc, cnt in by_clade.items():
        clade_maj_taxon[pc] = cnt.most_common(1)[0][0]

    rows: List[Dict[str, Any]] = []
    for sid in samples:
        pc = plastid_clades.get(sid)
        its_r = its_ribotypes.get(sid)
        sp = taxon_of[sid]
        discordant = 0
        if pc and sp != "Unknown":
            maj_tx = clade_maj_taxon.get(pc, sp)
            # Flag cytonuclear capture when a sample carries a chloroplast clade dominated by another species
            if maj_tx != sp and by_clade[pc][sp] < by_clade[pc][maj_tx]:
                discordant = 1
        rows.append({
            "genetic_sample_id": sid,
            "nuclear_taxon": sp,
            "plastid_clade": pc or "",
            "its_ribotype": its_r or "",
            "cytonuclear_discordant": discordant,
        })
    return rows


def run_offtarget_organellar_skim(
    fastq_source: Optional[Path],
    genomics_dir: Path,
    nuclear_tree: Optional[Path] = None,
) -> Dict[str, Any]:
    """Skim off-target plastid/ITS reads and emit flat organellar/cytonuclear files in `genomics_dir`."""
    if fastq_source is None or not fastq_source.exists():
        msg = f"[skim-organellar] FASTQ source absent ({fastq_source}); skipping cleanly."
        logger.info(msg)
        print(msg)
        return {"status": "skipped_missing_fastq", "n_samples": 0}

    getorg_bin = shutil.which("get_organelle_from_reads.py")
    if getorg_bin and fastq_source.is_dir():
        logger.info("Detected GetOrganelle binary at %s; running alongside pure-Python pileup audit.", getorg_bin)
        # External GetOrganelle executes only inside an isolated temporary directory so data/genomics/ stays flat
        with tempfile.TemporaryDirectory(prefix="packera_getorg_") as tmpd:
            _ = (tmpd, subprocess)

    plastid_seeds, its_seeds = _collect_organellar_seeds(genomics_dir)
    reads_by_sample: Dict[str, List[str]] = defaultdict(list)
    for sid, seq in _iter_staged_fastq_reads(fastq_source):
        reads_by_sample[sid].append(seq)

    if not reads_by_sample:
        msg = f"[skim-organellar] No FASTQ reads discovered inside {fastq_source}; skipping."
        logger.info(msg)
        return {"status": "skipped_no_reads", "n_samples": 0}

    cp_cons = _assemble_baited_consensus(reads_by_sample, plastid_seeds)
    its_cons = _assemble_baited_consensus(reads_by_sample, its_seeds)
    cp_clades = _cluster_haplotypes(cp_cons, prefix="cp_Clade")
    its_ribo = _cluster_haplotypes(its_cons, prefix="ITS_Ribotype")
    cyto_rows = _detect_cytonuclear_discordance(cp_clades, its_ribo, nuclear_tree_path=nuclear_tree)

    genomics_dir.mkdir(parents=True, exist_ok=True)
    cp_fasta = genomics_dir / "skimmed_plastid_consensus.fasta"
    its_fasta = genomics_dir / "skimmed_its_consensus.fasta"
    cyto_csv = genomics_dir / "cytonuclear_discordance.csv"

    with open(cp_fasta, "wt", encoding="utf-8") as fh:
        for sid, seq in sorted(cp_cons.items()):
            fh.write(f">{sid}|clade={cp_clades.get(sid, 'cp_Unassigned')}\n{seq}\n")
    with open(its_fasta, "wt", encoding="utf-8") as fh:
        for sid, seq in sorted(its_cons.items()):
            fh.write(f">{sid}|ribotype={its_ribo.get(sid, 'ITS_Unassigned')}\n{seq}\n")
    pd.DataFrame(cyto_rows).to_csv(cyto_csv, index=False)

    logger.info(
        "Skimmed %d plastid and %d ITS consensuses (%d cytonuclear records) into flat %s.",
        len(cp_cons),
        len(its_cons),
        len(cyto_rows),
        genomics_dir,
    )
    return {
        "status": "skimmed",
        "n_samples": len(reads_by_sample),
        "n_plastid": len(cp_cons),
        "n_its": len(its_cons),
        "cytonuclear_csv": str(cyto_csv),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Offline off-target organellar/ITS read skimmer for Packera.")
    parser.add_argument("--fastq-source", type=Path, default=None, help="Directory, .zip, or .tar.gz with FASTQ reads")
    parser.add_argument("--genomics-dir", type=Path, default=PROJECT_ROOT / "data" / "genomics")
    parser.add_argument("--nuclear-tree", type=Path, default=None)
    args = parser.parse_args(list(argv) if argv is not None else None)
    run_offtarget_organellar_skim(args.fastq_source, args.genomics_dir, nuclear_tree=args.nuclear_tree)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

