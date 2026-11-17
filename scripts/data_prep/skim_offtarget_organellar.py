#!/usr/bin/env python3
"""
Module: scripts/data_prep/skim_offtarget_organellar.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  Offline off-target chloroplast (cpDNA), nuclear ribosomal ITS, and nuclear
  waxy skimming from raw nextRAD FASTQ archives (.zip, .tar.gz) or FASTQ
  directories. Streams compressed archives without unpacking raw reads into
  `data/genomics/` using a pure-Python 17-mer baiting and pileup consensus
  engine seeded by Asteraceae organellar/ITS/waxy references and
  `vcfref_packera_blast.txt`. Writes outputs into
  `data/genomics/public_organellar_and_loci/`.
"""

from __future__ import annotations

import argparse, gzip, importlib, io, multiprocessing as mp, os, re, shutil, subprocess, sys, tarfile, zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.logger import setup_logging
from scripts.data_prep.organellar_cytonuclear_utils import (
    _NON_EPITHET_TOKENS, _SEED_SYNONYM_TAXA, _cluster_haplotypes,
    _detect_cytonuclear_discordance, _extract_species_from_seed_key,
    _match_best_public_lineage, _revcomp, match_best_public_lineage_with_metrics,
)

logger = setup_logging(name="packera.genomics.organellar")

_SEED_PLASTID_ANCHORS: Dict[str, str] = {
    "Packera_rbcL_core": "ATGTCACCACAAACAGAGACTAAAGCAAGTGTTGGATTCAAAGCTGGTGTTAAAGAGTACAAATTGACTTATTATACTCCTGAGTATGAAACCAAAGATACTGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCCGGAGTTCCACCTGAAGAAGCAGGGGCTGCGGTAGCTGCCGAATCTTCTACTGGTACATGGACAACTGTGTGGACTGATGGACTTACCAGCCTTGATCGTTACAAAGG",
    "Packera_trnL_core": "GGGCAATCCTGAGCCAAATCCTGTTTTCCGAAAACAAACAAAGGTTCAGAAAGCGAGAATAAAAAAAAGGATAGGTGCAGAGACTCAATGGAAGCTGTTCTAACGAATCGAGTTGACTGCGTTGCGTTAGTAAGCGGAAGAACAGGCTTCCATAGTC",
}
_SEED_ITS_ANCHORS: Dict[str, str] = {"Packera_ITS_58S_core": "CGTGAACCATCGAGTCTTTGAACGCAAGTTGCGCCCGAAGCCATTAGGTTGAGGGCACGTCTGCCTGGGCGTCACGCATCGCGTCGCCCCCCACCATCCATCCCCTTAGGGATGGGGTGTGTTGGGGGCGGGGGGATATTGGCCTCCCGTGCGCTCGGTGCGCGGTTGGCCTAAATAGGAGTCCTCGGCGACGGACGTCGCGACAAGTGGTGGTTGTAAAACCCTCGTTCTCGTCGTCGTGC"}
_SEED_WAXY_ANCHORS: Dict[str, str] = {"Packera_waxy_GBSSI_core": "ATGGTGGATTGTGTTCCATTTGTTGGTGGTGATGGAAAGGTTTTTGGTACTGCTGATGGACCTTTTGAAGGACTTAAGAAAGAGCTTCTCGTTGATGGTGTTGCTGATCCATCTAAGGTTTTCTTCGTTGCAACTCTTGAGGCTGCTAAGCGTGGTGGTATTGCACCAAAGATTGATCCAGCTCTTATGGGTTCTACTG"}
_MAX_BUCKET_READS = 25000
_FQ_PAT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)
_WORKER_CTX: Dict[str, Any] = {}


def _resolve_threads(threads: Optional[int] = None) -> int:
    if threads is not None and int(threads) >= 1: return int(threads)
    env_w = (os.environ.get("PLANTMATICS_MAX_WORKERS") or os.environ.get("PACKERA_MAX_WORKERS") or "").strip()
    return int(env_w) if env_w.isdigit() and int(env_w) >= 1 else max(1, min(os.cpu_count() or 1, 32))

def _clean_sample_id(raw: str) -> str:
    return importlib.import_module("scripts.data_prep.06_ingest_genetic_data").clean_sample_id(raw)

def _read_fasta(path: Path, keep_label: bool = False) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.is_file(): return out
    cur: Optional[str] = None; buf: List[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not (s := line.strip()): continue
            if s.startswith(">"):
                if cur is not None: out[cur] = "".join(buf).replace("-", "").upper()
                cur = re.sub(r"\s+", "|", s[1:].strip()) if keep_label else s[1:].split()[0]; buf = []
            else: buf.append(s)
    if cur is not None: out[cur] = "".join(buf).replace("-", "").upper()
    return out

def _find_geno_file(genomics_dir: Path, fname: str) -> Optional[Path]:
    direct = genomics_dir / fname
    return direct if direct.is_file() else next((p for p in sorted(genomics_dir.rglob(fname)) if p.is_file()), None)


def _is_ets_only_header(hdr: str) -> bool:
    h = re.sub(r"[|_]+", " ", hdr.lower())
    return ("external transcribed spacer" in h or "ets" in h.split()) and ("internal transcribed spacer" not in h and not any(t in h.split() for t in ("its", "its1", "its2")))


def _collect_organellar_seeds(genomics_dir: Path) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Gather plastid, ITS, and nuclear waxy seed contigs from BLAST hits, staged FASTAs, public Packera archives, and conserved anchors."""
    plastid_seeds, its_seeds, waxy_seeds = dict(_SEED_PLASTID_ANCHORS), dict(_SEED_ITS_ANCHORS), dict(_SEED_WAXY_ANCHORS)
    ref_seqs: Dict[str, str] = {}
    for fname in ("ref_packera.fasta", "vcfref_packera_rm.txt.fasta"):
        if (cand := _find_geno_file(genomics_dir, fname)) is not None:
            ref_seqs.update(_read_fasta(cand))
    if (bain := _find_geno_file(genomics_dir, "bain_jansen_sequences.fasta")) is not None:
        its_seeds.update({f"bain_{k}": v for k, v in _read_fasta(bain, keep_label=True).items() if len(v) >= 60 and not _is_ets_only_header(k)})
    for fname, bucket in (
        ("public_packera_18complete_plastomes.fasta", plastid_seeds), ("public_packera_111mapping_plastomes.fasta", plastid_seeds),
        ("public_packera_plastomes_and_cpdna.fasta", plastid_seeds), ("public_packera_its_ets.fasta", its_seeds), ("public_packera_waxy_nuclear.fasta", waxy_seeds),
    ):
        if (cand := _find_geno_file(genomics_dir, fname)) is not None:
            bucket.update({k: v for k, v in _read_fasta(cand, keep_label=True).items() if len(v) >= 60 and ("111mapping" not in fname or _extract_species_from_seed_key(k)) and (bucket is not its_seeds or not _is_ets_only_header(k))})
    if (bold_fa := _find_geno_file(genomics_dir, "public_packera_bold_barcodes.fasta")) is not None:
        for k, v in _read_fasta(bold_fa, keep_label=True).items():
            if len(v) < 60 or _is_ets_only_header(k): continue
            low_k = k.lower()
            target = its_seeds if re.search(r"(?:^|[|_ \t])(its\d?|5\.8s|ets|ribosomal)", low_k) else (waxy_seeds if re.search(r"(waxy|gbssi)", low_k) else plastid_seeds)
            target[k] = v
    pub_cp, pub_its, pub_wx = len(plastid_seeds) > len(_SEED_PLASTID_ANCHORS), len(its_seeds) > len(_SEED_ITS_ANCHORS), len(waxy_seeds) > len(_SEED_WAXY_ANCHORS)
    if (blast_txt := _find_geno_file(genomics_dir, "vcfref_packera_blast.txt")) is not None and ref_seqs:
        qtexts, cur_q = defaultdict(list), ""
        with open(blast_txt, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                s = line.strip()
                if s.startswith("Query="): cur_q = s[6:].strip().split()[0] if s[6:].strip() else ""
                elif s.startswith("Locus:"):
                    cur_q = ""; parts = s.split("\t")
                    if len(parts) >= 3: qtexts[parts[1].strip()].append(" ".join(parts[2:]))
                elif cur_q: qtexts[cur_q].append(s)
                elif "\t" in line and (p0 := line.rstrip("\r\n").split("\t")[0].strip()) in ref_seqs: qtexts[p0].append(line)
        for qid, chunks in qtexts.items():
            low = " ".join(chunks).lower()
            if qid not in ref_seqs or re.search(r"\b(?:ribosomal\s+protein|proteasome|chloroplastic-like)\b", low): continue
            if not pub_cp and re.search(r"\b(?:chloroplast|plastid|rbcl|matk|trnl|psba|ndhf)\b", low): plastid_seeds[qid] = ref_seqs[qid]
            elif not pub_its and re.search(r"\b(?:internal\s+transcribed\s+spacer|ribosomal|rrna|(?:18|26|28)s|its[12])\b", low): its_seeds[qid] = ref_seqs[qid]
            elif not pub_wx and re.search(r"\b(?:waxy|gbssi|starch\s+synthase)\b", low): waxy_seeds[qid] = ref_seqs[qid]
    return plastid_seeds, its_seeds, waxy_seeds


def _infer_sample_from_fq_name(member_name: str) -> str:
    base = Path(member_name).name
    stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    stem = re.sub(r"([_.-](?:S\d+|L\d{3}|R?[12]|F|R|merged|trimmed|paired|unpaired|clean))+$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"^[A-Z][a-z]+(?=[A-Z]{3}[A-Z0-9]{4,}$)", "", stem)
    return _clean_sample_id(stem)


def _iter_fastq_records_stream(fh: io.TextIOBase) -> Iterator[Tuple[str, str]]:
    while True:
        h, seq, _plus, qual = fh.readline(), fh.readline(), fh.readline(), fh.readline()
        if not (h and qual): break
        yield h.strip(), seq.strip().upper()

def _iter_binary_fastq_stream(raw_f: Any, is_gz: bool) -> Iterator[str]:
    bio: Any = gzip.GzipFile(fileobj=raw_f) if is_gz else raw_f
    with io.TextIOWrapper(bio, encoding="utf-8", errors="replace") as fh:
        yield from (seq for _, seq in _iter_fastq_records_stream(fh) if seq)

def _iter_fastq_file_seqs(p: Path) -> Iterator[str]:
    """Stream uppercase FASTQ sequences from a single `.fq`, `.fastq`, or `.gz` file."""
    is_gz = p.name.lower().endswith(".gz")
    pigz_bin = shutil.which("pigz") if is_gz else None
    gzip_bin = None if (not is_gz or pigz_bin) else shutil.which("gzip")
    if is_gz and (pigz_bin or gzip_bin):
        cmd = [pigz_bin, "-p", "2", "-dc", str(p)] if pigz_bin else [str(gzip_bin), "-dc", str(p)]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, encoding="utf-8", errors="replace", bufsize=1 << 20)
        try:
            assert proc.stdout is not None
            yield from (seq for _, seq in _iter_fastq_records_stream(proc.stdout) if seq)
        finally:
            if proc.stdout is not None: proc.stdout.close()
            if proc.poll() is None: proc.terminate()
            proc.wait()
        return
    with (gzip.open if is_gz else open)(p, "rt", encoding="utf-8", errors="replace") as fh:
        yield from (seq for _, seq in _iter_fastq_records_stream(fh) if seq)

def _iter_staged_fastq_reads(fastq_source: Path) -> Iterator[Tuple[str, str]]:
    """Stream `(sample_id, seq)` from `.zip`, `.tar.gz`, or directory without unpacking to `data/genomics/`."""
    if not fastq_source.exists(): return
    if fastq_source.is_dir():
        for p in sorted(fastq_source.rglob("*")):
            if p.is_file() and _FQ_PAT.search(p.name):
                sid = _infer_sample_from_fq_name(p.name); yield from ((sid, seq) for seq in _iter_fastq_file_seqs(p))
        return
    low = fastq_source.name.lower()
    if low.endswith(".zip"):
        with zipfile.ZipFile(fastq_source, "r") as zf:
            for info in sorted(zf.infolist(), key=lambda x: x.filename):
                if not info.is_dir() and _FQ_PAT.search(info.filename):
                    sid = _infer_sample_from_fq_name(info.filename)
                    with zf.open(info, "r") as raw_f:
                        yield from ((sid, seq) for seq in _iter_binary_fastq_stream(raw_f, info.filename.lower().endswith(".gz")))
        return
    if low.endswith((".tar.gz", ".tgz", ".tar")):
        with tarfile.open(fastq_source, "r:gz" if low.endswith((".tar.gz", ".tgz")) else "r:") as tf:
            for member in tf.getmembers():
                if member.isfile() and _FQ_PAT.search(member.name) and (raw_f := tf.extractfile(member)) is not None:
                    sid = _infer_sample_from_fq_name(member.name)
                    with raw_f: yield from ((sid, seq) for seq in _iter_binary_fastq_stream(raw_f, member.name.lower().endswith(".gz")))


def _build_kmer_index(seeds: Mapping[str, str], k: int = 17, step_override: Optional[int] = None) -> Dict[str, List[Tuple[str, int, bool]]]:
    idx: Dict[str, List[Tuple[str, int, bool]]] = defaultdict(list)
    for name, seq in seeds.items():
        if len(seq) < k: continue
        step = step_override if step_override is not None else (4 if len(seq) > 12000 else 1)
        for pos in range(0, len(seq) - k + 1, step):
            kmer = seq[pos : pos + k]
            if "N" in kmer: continue
            idx[kmer].append((name, pos, False))
            idx[_revcomp(kmer)].append((name, pos, True))
    return idx


def _bait_stream_into_buckets(
    seq_iter: Iterable[str], seed_kmers: Any, cp_kidx: Any, its_kidx: Any, waxy_kidx: Any,
    cp_b: List[str], its_b: List[str], waxy_b: List[str], cap: int = _MAX_BUCKET_READS,
) -> int:
    """Bait reads into plastid, ITS, and waxy buckets using a stride-5 prefilter and early cap stop."""
    n_reads = 0
    for seq in seq_iter:
        if len(cp_b) >= cap and len(its_b) >= cap and len(waxy_b) >= cap: break
        n_reads += 1
        if (n_kmers := len(seq) - 16) <= 0: continue
        stride = 5 if n_kmers >= 16 else 1
        if not any(seq[i : i + 17] in seed_kmers for i in range(0, n_kmers, stride)): continue
        for kidx_b, bucket in ((cp_kidx, cp_b), (its_kidx, its_b), (waxy_kidx, waxy_b)):
            if len(bucket) < cap and any(seq[i : i + 17] in kidx_b for i in range(n_kmers)): bucket.append(seq)
    return n_reads


def _init_skim_worker(ctx: Dict[str, Any]) -> None:
    _WORKER_CTX.clear(); _WORKER_CTX.update(ctx)


def _skim_file_worker(path_str: str) -> Tuple[str, str, int, List[str], List[str], List[str]]:
    p = Path(path_str); cp_b: List[str] = []; its_b: List[str] = []; waxy_b: List[str] = []
    seq_iter = _iter_fastq_file_seqs(p)
    try: n_reads = _bait_stream_into_buckets(seq_iter, _WORKER_CTX["seed_kmers"], _WORKER_CTX["cp_kidx"], _WORKER_CTX["its_kidx"], _WORKER_CTX["waxy_kidx"], cp_b, its_b, waxy_b)
    finally: seq_iter.close()
    return p.name, _infer_sample_from_fq_name(p.name), n_reads, cp_b, its_b, waxy_b


def _assemble_baited_consensus(
    reads_by_sample: Mapping[str, List[str]], seeds: Mapping[str, str], k: int = 17,
    kidx: Optional[Mapping[str, List[Tuple[str, int, bool]]]] = None, min_reads: int = 4, min_bases: int = 60,
    pos_out: Optional[Dict[str, Dict[int, str]]] = None,
) -> Dict[str, str]:
    """Bait reads by seed k-mers and construct per-sample consensus sequence by read pileup."""
    if not seeds or not reads_by_sample: return {}
    kidx = kidx if kidx is not None else _build_kmer_index(seeds, k=k)
    canon_idx = cidx_uniq = canon_len = None
    max_slen = max((len(v) for v in seeds.values()), default=0)
    if pos_out is not None and (max_slen > 20000 or len(reads_by_sample) > 1):
        ckey = "Packera_dubia" if (max_slen > 20000 and "Packera_dubia" in seeds) else max(seeds, key=lambda s: (len(seeds[s]), s))
        canon_len = len(seeds[ckey])
        canon_idx = _build_kmer_index({ckey: seeds[ckey]}, k=k, step_override=1) if (canon_len > 12000 or len(seeds) > 1) else kidx
        if max_slen > 20000: cidx_uniq = {km: hits[0] for km, hits in canon_idx.items() if len(hits) == 1}
    cons_out: Dict[str, str] = {}
    ref_kcache: Dict[str, Mapping[str, List[Tuple[str, int, bool]]]] = {}
    for sid, seqs in sorted(reads_by_sample.items()):
        ref_scores: Dict[str, float] = defaultdict(float)
        for raw_seq in seqs[:400]:
            for r_start in range(0, max(0, len(raw_seq) - k + 1), 16):
                if not (hits := kidx.get(raw_seq[r_start : r_start + k])): continue
                w = 1.0 / float(len(hits))
                for sname, _, _ in hits: ref_scores[sname] += w
        if not ref_scores: continue
        best_ref = max(ref_scores, key=lambda s: (ref_scores[s], len(seeds[s])))
        seed_seq = seeds[best_ref]; target_len = len(seed_seq)
        if target_len > 12000:
            if best_ref not in ref_kcache: ref_kcache[best_ref] = _build_kmer_index({best_ref: seed_seq}, k=k, step_override=1)
            ref_idx = ref_kcache[best_ref]
        else:
            ref_idx = kidx
        sparse, csparse, n_aligned, c_aligned = defaultdict(Counter), defaultdict(Counter), 0, 0
        req_two = bool(cidx_uniq is not None and len(seqs) >= 25)
        for raw_seq in seqs:
            for r_start in range(0, max(0, len(raw_seq) - k + 1)):
                if not (hits := ref_idx.get(raw_seq[r_start : r_start + k])): continue
                if (hit := next((h for h in hits if h[0] == best_ref), None)) is None: continue
                _, spos, is_rc = hit; oriented = _revcomp(raw_seq) if is_rc else raw_seq
                offset = spos - ((len(raw_seq) - r_start - k) if is_rc else r_start); n_aligned += 1
                for idx_b, ch in enumerate(oriented):
                    if 0 <= (gpos := offset + idx_b) < target_len and ch in "ACGT": sparse[gpos][ch] += 1
                break
            if cidx_uniq is not None and canon_len is not None:
                seen_offs: List[Tuple[int, bool]] = []
                picked: Optional[Tuple[int, bool]] = None
                for r_start in range(0, max(0, len(raw_seq) - k + 1), 8):
                    if (ch := cidx_uniq.get(raw_seq[r_start : r_start + k])) is None: continue
                    coff = ch[1] - ((len(raw_seq) - r_start - k) if ch[2] else r_start)
                    if not req_two: picked = (coff, ch[2]); break
                    if any(o[1] == ch[2] and abs(o[0] - coff) <= 2 for o in seen_offs): picked = (coff, ch[2]); break
                    seen_offs.append((coff, ch[2]))
                if picked is not None:
                    cor = _revcomp(raw_seq) if picked[1] else raw_seq; c_aligned += 1
                    for idx_b, ch in enumerate(cor):
                        if 0 <= (cg := picked[0] + idx_b) < canon_len and ch in "ACGT": csparse[cg][ch] += 1
            elif canon_idx is not None and canon_len is not None:
                for r_start in range(0, max(0, len(raw_seq) - k + 1), 8):
                    if not (chits := canon_idx.get(raw_seq[r_start : r_start + k])): continue
                    _, cspos, crc = chits[0]; cor = _revcomp(raw_seq) if crc else raw_seq
                    coff = cspos - ((len(raw_seq) - r_start - k) if crc else r_start)
                    for idx_b, ch in enumerate(cor):
                        if 0 <= (cg := coff + idx_b) < canon_len and ch in "ACGT": csparse[cg][ch] += 1
                    break
        if n_aligned < min_reads or len(sparse) < min_bases: continue
        chars, prev = [], None
        for j in sorted(sparse):
            if prev is not None and (gap := j - prev - 1) > 0: chars.append("N" * min(gap, 17))
            chars.append(sparse[j].most_common(1)[0][0]); prev = j
        if len(cons_s := "".join(chars).strip("N")) >= min_bases and sum(c != "N" for c in cons_s) >= min_bases:
            cons_out[sid] = cons_s
            if pos_out is not None and csparse:
                min_d = (5 if max(n_aligned, c_aligned) >= 500 else (2 if max(n_aligned, c_aligned) >= 30 else 1)) if cidx_uniq is not None else 1
                pos_out[sid] = {p: c.most_common(1)[0][0] for p, c in csparse.items() if (tot := sum(c.values())) >= min_d and c.most_common(1)[0][1] / float(tot) >= (0.80 if min_d > 1 else 0.0)}
    return cons_out


def _discover_empirical_nuclear_taxon_map(genomics_dir: Path, out_dir: Optional[Path]) -> Dict[str, str]:
    """Load high-confidence non-admixed Step 06 nuclear taxon predictions when staged nearby."""
    cands: List[Path] = []
    for base in (out_dir, genomics_dir, genomics_dir.parent / "outputs" / "genomics", genomics_dir.parent / "outputs" / "tables"):
        if base is not None:
            cands.append(Path(base) / "genomic_assignments.csv")
    out: Dict[str, str] = {}
    for p in cands:
        if not p.is_file():
            continue
        try:
            df = pd.read_csv(p, dtype=str, keep_default_na=False)
        except Exception:
            continue
        if "genetic_sample_id" not in df.columns or "geno_predicted_taxon" not in df.columns:
            continue
        for _, row in df.iterrows():
            sid = _clean_sample_id(row.get("genetic_sample_id", ""))
            tx = str(row.get("geno_predicted_taxon", "")).strip()
            if not sid or not tx.startswith("Packera "):
                continue
            try:
                admix = int(float(row.get("is_geno_admixed") or 0))
                prob = float(row.get("geno_posterior_prob") or 0.0)
            except (TypeError, ValueError):
                continue
            if admix == 0 and prob >= 0.70 and sid not in out:
                out[sid] = tx
        if out:
            break
    return out


def run_offtarget_organellar_skim(
    fastq_source: Optional[Path], genomics_dir: Path, nuclear_tree: Optional[Path] = None,
    out_dir: Optional[Path] = None, threads: Optional[int] = None,
) -> Dict[str, Any]:
    """Skim off-target plastid, ITS, and nuclear waxy reads into `public_organellar_and_loci/` (or `--out-dir`)."""
    if fastq_source is None or not fastq_source.exists():
        msg = f"[skim-organellar] FASTQ source absent ({fastq_source}); skipping cleanly."; logger.info(msg); print(msg)
        return {"status": "skipped_missing_fastq", "n_samples": 0}
    plastid_seeds, its_seeds, waxy_seeds = _collect_organellar_seeds(genomics_dir)
    cp_kidx, its_kidx, waxy_kidx = _build_kmer_index(plastid_seeds), _build_kmer_index(its_seeds), _build_kmer_index(waxy_seeds)
    seed_kmers = set(cp_kidx) | set(its_kidx) | set(waxy_kidx); seen_samples: set[str] = set()
    cp_reads, its_reads, waxy_reads = defaultdict(list), defaultdict(list), defaultdict(list)
    eff_threads = _resolve_threads(threads)
    fq_files = sorted((p for p in fastq_source.rglob("*") if p.is_file() and _FQ_PAT.search(p.name)), key=lambda p: (-p.stat().st_size, p.name)) if fastq_source.is_dir() else []
    if len(fq_files) > 1 and eff_threads > 1:
        n_workers = min(eff_threads, len(fq_files))
        logger.info("[skim-organellar] Scanning %d FASTQ file(s) across %d worker(s)...", len(fq_files), n_workers)
        _WORKER_CTX.clear(); _WORKER_CTX.update({"seed_kmers": seed_kmers, "cp_kidx": cp_kidx, "its_kidx": its_kidx, "waxy_kidx": waxy_kidx})
        use_fork = "fork" in mp.get_all_start_methods()
        ctx_mp = mp.get_context("fork") if use_fork else mp.get_context()
        pool_kw: Dict[str, Any] = {"processes": n_workers} if use_fork else {"processes": n_workers, "initializer": _init_skim_worker, "initargs": (dict(_WORKER_CTX),)}
        try:
            with ctx_mp.Pool(**pool_kw) as pool:
                for done_i, (fname, sid, n_reads, cp_b, its_b, waxy_b) in enumerate(pool.imap_unordered(_skim_file_worker, [str(p) for p in fq_files]), start=1):
                    seen_samples.add(sid)
                    for src_b, dst_b in ((cp_b, cp_reads[sid]), (its_b, its_reads[sid]), (waxy_b, waxy_reads[sid])):
                        if src_b and len(dst_b) < _MAX_BUCKET_READS: dst_b.extend(src_b[: _MAX_BUCKET_READS - len(dst_b)])
                    logger.info("[skim-organellar] [%d/%d] %s (%s): scanned %d reads (cp=%d, its=%d, waxy=%d)", done_i, len(fq_files), fname, sid, n_reads, len(cp_b), len(its_b), len(waxy_b))
        finally: _WORKER_CTX.clear()
    else:
        for sid, seq in _iter_staged_fastq_reads(fastq_source):
            seen_samples.add(sid); cp_b, its_b, waxy_b = cp_reads[sid], its_reads[sid], waxy_reads[sid]
            if len(cp_b) < _MAX_BUCKET_READS or len(its_b) < _MAX_BUCKET_READS or len(waxy_b) < _MAX_BUCKET_READS:
                _bait_stream_into_buckets((seq,), seed_kmers, cp_kidx, its_kidx, waxy_kidx, cp_b, its_b, waxy_b)
    if not seen_samples:
        logger.info("[skim-organellar] No FASTQ reads discovered inside %s; skipping.", fastq_source)
        return {"status": "skipped_no_reads", "n_samples": 0}

    cp_pos, its_pos, wx_pos = {}, {}, {}
    cp_cons = _assemble_baited_consensus(cp_reads, plastid_seeds, kidx=cp_kidx, min_reads=4, min_bases=80, pos_out=cp_pos)
    its_cons = _assemble_baited_consensus(its_reads, its_seeds, kidx=its_kidx, min_reads=4, min_bases=60, pos_out=its_pos)
    waxy_cons = _assemble_baited_consensus(waxy_reads, waxy_seeds, kidx=waxy_kidx, min_reads=4, min_bases=60, pos_out=wx_pos)
    cp_clades = _cluster_haplotypes(cp_cons, "cp_Clade", pos_maps=cp_pos, min_id=0.990)
    its_ribo = _cluster_haplotypes(its_cons, "ITS_Ribotype", pos_maps=its_pos, min_id=0.990)
    waxy_clades = _cluster_haplotypes(waxy_cons, "Waxy_Lineage", pos_maps=wx_pos, min_id=0.990)
    cp_cache, its_cache, waxy_cache = {}, {}, {}
    cp_full = {sid: match_best_public_lineage_with_metrics(seq, plastid_seeds, _kcache=cp_cache, q_pos=cp_pos.get(sid)) for sid, seq in cp_cons.items()}
    pub_cp_match = {sid: vals[0] for sid, vals in cp_full.items()}
    cp_metrics = {sid: (vals[1], vals[2]) for sid, vals in cp_full.items()}
    pub_its_match = {sid: _match_best_public_lineage(seq, its_seeds, _kcache=its_cache) for sid, seq in its_cons.items()}
    waxy_match = {sid: _match_best_public_lineage(seq, waxy_seeds, _kcache=waxy_cache) or waxy_clades.get(sid, "") for sid, seq in waxy_cons.items()}
    nuc_taxon_map = _discover_empirical_nuclear_taxon_map(genomics_dir, out_dir)
    cyto_rows = _detect_cytonuclear_discordance(
        cp_clades, its_ribo, nuclear_tree_path=nuclear_tree,
        public_plastome_matches=pub_cp_match, public_its_matches=pub_its_match,
        waxy_lineage_matches=waxy_match, plastome_match_metrics=cp_metrics,
        nuclear_taxon_map=nuc_taxon_map,
    )

    target_dir = out_dir if out_dir is not None else (genomics_dir if genomics_dir.name == "public_organellar_and_loci" else (genomics_dir / "public_organellar_and_loci"))
    target_dir.mkdir(parents=True, exist_ok=True)
    cyto_csv = target_dir / "cytonuclear_discordance.csv"
    for fname_out, cons_map, tag, lbl_map, fallback in (("skimmed_plastid_consensus.fasta", cp_cons, "clade", cp_clades, "cp_Unassigned"), ("skimmed_its_consensus.fasta", its_cons, "ribotype", its_ribo, "ITS_Unassigned"), ("skimmed_waxy_consensus.fasta", waxy_cons, "waxy", waxy_clades, "Waxy_Unassigned")):
        with open(target_dir / fname_out, "wt", encoding="utf-8") as fh:
            for sid, seq in sorted(cons_map.items()): fh.write(f">{sid}|{tag}={lbl_map.get(sid, fallback)}\n{seq}\n")
    cols = ["genetic_sample_id", "nuclear_taxon", "plastid_clade", "its_ribotype", "cytonuclear_discordant", "public_plastome_match", "public_its_match", "waxy_lineage_match"]
    (pd.DataFrame(cyto_rows) if cyto_rows else pd.DataFrame(columns=cols)).to_csv(cyto_csv, index=False)
    logger.info("Skimmed %d plastid, %d ITS, and %d waxy consensuses (%d cytonuclear records) into %s.", len(cp_cons), len(its_cons), len(waxy_cons), len(cyto_rows), target_dir)
    return {"status": "skimmed", "n_samples": len(seen_samples), "n_plastid": len(cp_cons), "n_its": len(its_cons), "n_waxy": len(waxy_cons), "out_dir": str(target_dir), "cytonuclear_csv": str(cyto_csv)}

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Offline off-target organellar/ITS read skimmer for Packera.")
    parser.add_argument("--fastq-source", type=Path, default=None, help="Directory, .zip, or .tar.gz with FASTQ reads")
    parser.add_argument("--genomics-dir", type=Path, default=PROJECT_ROOT / "data" / "genomics")
    parser.add_argument("--nuclear-tree", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=None, help="Output subdirectory (default: <genomics-dir>/public_organellar_and_loci)")
    parser.add_argument("-j", "--threads", "--workers", dest="threads", type=int, default=None, help="Maximum CPU worker processes (default: auto)")
    args = parser.parse_args(list(argv) if argv is not None else None)
    run_offtarget_organellar_skim(args.fastq_source, args.genomics_dir, nuclear_tree=args.nuclear_tree, out_dir=args.out_dir, threads=args.threads)
    return 0

if __name__ == "__main__": raise SystemExit(main())

