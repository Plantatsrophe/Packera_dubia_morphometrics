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
  organellar/ITS/waxy references and `vcfref_packera_blast.txt`. Writes outputs
  into `data/genomics/public_organellar_and_loci/`.
"""

from __future__ import annotations

import argparse, gzip, importlib, io, multiprocessing as mp, os, re, shutil, subprocess, sys, tarfile, tempfile, zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple

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

_SEED_WAXY_ANCHORS: Dict[str, str] = {
    # Conserved Asteraceae granule-bound starch synthase (GBSSI / waxy) core motif
    "Packera_waxy_GBSSI_core": "ATGGTGGATTGTGTTCCATTTGTTGGTGGTGATGGAAAGGTTTTTGGTACTGCTGATGGACCTTTTGAAGGACTTAAGAAAGAGCTTCTCGTTGATGGTGTTGCTGATCCATCTAAGGTTTTCTTCGTTGCAACTCTTGAGGCTGCTAAGCGTGGTGGTATTGCACCAAAGATTGATCCAGCTCTTATGGGTTCTACTG",
}

_SEED_SYNONYM_TAXA: Dict[str, str] = {
    "senecio tomentosus": "Packera dubia", "packera tomentosa": "Packera dubia", "senecio aureus": "Packera aurea",
    "senecio anonymus": "Packera anonyma", "senecio pauperculus": "Packera paupercula", "senecio schweinitzianus": "Packera schweinitziana",
    "senecio millefolium": "Packera millefolium", "packera millefolia": "Packera millefolium", "packera memmingeri": "Packera x memmingeri",
    "packera xmemmingeri": "Packera x memmingeri", "senecio obovatus": "Packera obovata", "senecio plattensis": "Packera plattensis",
    "senecio antennariifolius": "Packera antennariifolia", "senecio glabellus": "Packera glabella", "senecio loratifolius": "Packera loratifolia",
    "senecio pseudaureus": "Packera pseudaurea", "senecio indecorus": "Packera indecora", "senecio canus": "Packera cana",
    "senecio strepthanthifolius": "Packera streptanthifolia", "senecio streptanthifolius": "Packera streptanthifolia",
    "senecio fendleri": "Packera fendleri", "senecio multilobatus": "Packera multilobata", "senecio neomexicanus": "Packera neomexicana",
    "senecio debilis": "Packera debilis", "senecio werneriifolius": "Packera werneriifolia", "senecio hartianus": "Packera hartiana",
}
_NON_EPITHET_TOKENS = frozenset({"rbcl", "trnl", "its", "its1", "its2", "58s", "waxy", "gbssi", "core", "plastid", "chloroplast", "plastome", "ribosomal", "nuclear", "sp", "unk", "unknown", "var", "subsp", "cf", "aff"})

_IUPAC_PAIR = {("A", "G"): "R", ("G", "A"): "R", ("C", "T"): "Y", ("T", "C"): "Y", ("G", "C"): "S", ("C", "G"): "S", ("A", "T"): "W", ("T", "A"): "W", ("G", "T"): "K", ("T", "G"): "K", ("A", "C"): "M", ("C", "A"): "M"}
_MAX_BUCKET_READS = 25000
_FQ_PAT = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)
_WORKER_CTX: Dict[str, Any] = {}


def _resolve_threads(threads: Optional[int] = None) -> int:
    if threads is not None and int(threads) >= 1:
        return int(threads)
    env_w = (os.environ.get("PLANTMATICS_MAX_WORKERS") or os.environ.get("PACKERA_MAX_WORKERS") or "").strip()
    return int(env_w) if env_w.isdigit() and int(env_w) >= 1 else max(1, min(os.cpu_count() or 1, 32))


def _clean_sample_id(raw: str) -> str:
    fn = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").clean_sample_id
    return fn(raw)


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _read_fasta(path: Path, keep_label: bool = False) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.is_file():
        return out
    cur: Optional[str] = None
    buf: List[str] = []
    with open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not (s := line.strip()): continue
            if s.startswith(">"):
                if cur is not None: out[cur] = "".join(buf).replace("-", "").upper()
                cur = re.sub(r"\s+", "|", s[1:].strip()) if keep_label else s[1:].split()[0]
                buf = []
            else:
                buf.append(s)
    if cur is not None: out[cur] = "".join(buf).replace("-", "").upper()
    return out


def _find_geno_file(genomics_dir: Path, fname: str) -> Optional[Path]:
    direct = genomics_dir / fname
    return direct if direct.is_file() else next((p for p in sorted(genomics_dir.rglob(fname)) if p.is_file()), None)


def _collect_organellar_seeds(genomics_dir: Path) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Gather plastid, ITS, and nuclear waxy seed contigs from BLAST hits, staged FASTAs, public Packera archives, and conserved anchors."""
    plastid_seeds, its_seeds, waxy_seeds = dict(_SEED_PLASTID_ANCHORS), dict(_SEED_ITS_ANCHORS), dict(_SEED_WAXY_ANCHORS)
    ref_seqs: Dict[str, str] = {}
    for fname in ("ref_packera.fasta", "vcfref_packera_rm.txt.fasta"):
        if (cand := _find_geno_file(genomics_dir, fname)) is not None:
            ref_seqs.update(_read_fasta(cand))
    if (bain := _find_geno_file(genomics_dir, "bain_jansen_sequences.fasta")) is not None:
        its_seeds.update({f"bain_{k}": v for k, v in _read_fasta(bain, keep_label=True).items() if len(v) >= 60})
    for fname, bucket in (
        ("public_packera_18complete_plastomes.fasta", plastid_seeds), ("public_packera_plastomes_and_cpdna.fasta", plastid_seeds),
        ("public_packera_its_ets.fasta", its_seeds), ("public_packera_waxy_nuclear.fasta", waxy_seeds),
    ):
        if (cand := _find_geno_file(genomics_dir, fname)) is not None:
            bucket.update({k: v for k, v in _read_fasta(cand, keep_label=True).items() if len(v) >= 60})
    if (bold_fa := _find_geno_file(genomics_dir, "public_packera_bold_barcodes.fasta")) is not None:
        for k, v in _read_fasta(bold_fa, keep_label=True).items():
            if len(v) < 60: continue
            low_k = k.lower()
            target = its_seeds if re.search(r"(?:^|[|_ \t])(its\d?|5\.8s|ets|ribosomal)", low_k) else (waxy_seeds if re.search(r"(waxy|gbssi)", low_k) else plastid_seeds)
            target[k] = v
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
            if re.search(r"\b(?:chloroplast|plastid|rbcl|matk|trnl|psba|ndhf)\b", low): plastid_seeds[qid] = ref_seqs[qid]
            elif re.search(r"\b(?:internal\s+transcribed\s+spacer|ribosomal|rrna|(?:18|26|28)s|its[12])\b", low): its_seeds[qid] = ref_seqs[qid]
            elif re.search(r"\b(?:waxy|gbssi|starch\s+synthase)\b", low): waxy_seeds[qid] = ref_seqs[qid]
    return plastid_seeds, its_seeds, waxy_seeds


def _infer_sample_from_fq_name(member_name: str) -> str:
    base = Path(member_name).name
    stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", base, flags=re.IGNORECASE)
    stem = re.sub(r"([_.-](?:S\d+|L\d{3}|R?[12]|F|R|merged|trimmed|paired|unpaired|clean))+$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"^[A-Z][a-z]+(?=[A-Z]{3}[A-Z0-9]{4,}$)", "", stem)
    return _clean_sample_id(stem)


def _iter_fastq_records_stream(fh: io.TextIOBase) -> Iterator[Tuple[str, str]]:
    while True:
        h = fh.readline()
        if not h:
            break
        seq = fh.readline()
        _plus = fh.readline()
        if not fh.readline():
            break
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
    if not fastq_source.exists():
        return
    if fastq_source.is_dir():
        for p in sorted(fastq_source.rglob("*")):
            if p.is_file() and _FQ_PAT.search(p.name):
                sid = _infer_sample_from_fq_name(p.name)
                for seq in _iter_fastq_file_seqs(p):
                    yield sid, seq
        return
    low = fastq_source.name.lower()
    if low.endswith(".zip"):
        with zipfile.ZipFile(fastq_source, "r") as zf:
            for info in sorted(zf.infolist(), key=lambda x: x.filename):
                if not info.is_dir() and _FQ_PAT.search(info.filename):
                    sid = _infer_sample_from_fq_name(info.filename)
                    with zf.open(info, "r") as raw_f:
                        for seq in _iter_binary_fastq_stream(raw_f, info.filename.lower().endswith(".gz")):
                            yield sid, seq
        return
    if low.endswith((".tar.gz", ".tgz", ".tar")):
        with tarfile.open(fastq_source, "r:gz" if low.endswith((".tar.gz", ".tgz")) else "r:") as tf:
            for member in tf.getmembers():
                if member.isfile() and _FQ_PAT.search(member.name) and (raw_f := tf.extractfile(member)) is not None:
                    sid = _infer_sample_from_fq_name(member.name)
                    with raw_f:
                        for seq in _iter_binary_fastq_stream(raw_f, member.name.lower().endswith(".gz")):
                            yield sid, seq


def _build_kmer_index(seeds: Mapping[str, str], k: int = 17) -> Dict[str, List[Tuple[str, int, bool]]]:
    idx: Dict[str, List[Tuple[str, int, bool]]] = defaultdict(list)
    for name, seq in seeds.items():
        if len(seq) < k: continue
        step = 4 if len(seq) > 12000 else 1
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
        if len(cp_b) >= cap and len(its_b) >= cap and len(waxy_b) >= cap:
            break
        n_reads += 1
        n_kmers = len(seq) - 16
        if n_kmers <= 0:
            continue
        stride = 5 if n_kmers >= 16 else 1
        if not any(seq[i : i + 17] in seed_kmers for i in range(0, n_kmers, stride)):
            continue
        for kidx_b, bucket in ((cp_kidx, cp_b), (its_kidx, its_b), (waxy_kidx, waxy_b)):
            if len(bucket) < cap and any(seq[i : i + 17] in kidx_b for i in range(n_kmers)):
                bucket.append(seq)
    return n_reads


def _init_skim_worker(ctx: Dict[str, Any]) -> None:
    _WORKER_CTX.clear()
    _WORKER_CTX.update(ctx)


def _skim_file_worker(path_str: str) -> Tuple[str, str, int, List[str], List[str], List[str]]:
    p = Path(path_str)
    cp_b: List[str] = []; its_b: List[str] = []; waxy_b: List[str] = []
    seq_iter = _iter_fastq_file_seqs(p)
    try:
        n_reads = _bait_stream_into_buckets(
            seq_iter, _WORKER_CTX["seed_kmers"], _WORKER_CTX["cp_kidx"], _WORKER_CTX["its_kidx"], _WORKER_CTX["waxy_kidx"], cp_b, its_b, waxy_b
        )
    finally:
        seq_iter.close()
    return p.name, _infer_sample_from_fq_name(p.name), n_reads, cp_b, its_b, waxy_b


def _assemble_baited_consensus(reads_by_sample: Mapping[str, List[str]], seeds: Mapping[str, str], k: int = 17, kidx: Optional[Mapping[str, List[Tuple[str, int, bool]]]] = None) -> Dict[str, str]:
    """Bait reads by seed k-mers and construct per-sample consensus sequence by read pileup."""
    if not seeds or not reads_by_sample: return {}
    kidx = kidx if kidx is not None else _build_kmer_index(seeds, k=k)
    cons_out: Dict[str, str] = {}
    for sid, seqs in sorted(reads_by_sample.items()):
        hit_counts: Counter[str] = Counter()
        aligned_reads: List[Tuple[str, int, str]] = []
        for raw_seq in seqs:
            for r_start in range(0, max(0, len(raw_seq) - k + 1)):
                if not (hits := kidx.get(raw_seq[r_start : r_start + k])): continue
                sname, spos, is_rc = hits[0]
                oriented = _revcomp(raw_seq) if is_rc else raw_seq
                anchor_r = (len(raw_seq) - r_start - k) if is_rc else r_start
                hit_counts[sname] += 1
                aligned_reads.append((sname, spos - anchor_r, oriented))
                break
        if not hit_counts: continue
        best_ref = hit_counts.most_common(1)[0][0]
        seed_seq = seeds[best_ref]; target_len = len(seed_seq)
        sparse: Dict[int, Counter[str]] = defaultdict(Counter)
        for sname, offset, seq in aligned_reads:
            if sname != best_ref: continue
            for idx_b, ch in enumerate(seq):
                if 0 <= (gpos := offset + idx_b) < target_len and ch in "ACGT": sparse[gpos][ch] += 1
        if len(sparse) < 20: continue
        lo, hi = min(sparse), max(sparse)
        pos_iter = range(target_len) if target_len <= 5000 else (range(lo, hi + 1) if (hi - lo) <= 10000 else sorted(sparse))
        chars = []
        for j in pos_iter:
            if not (cnt := sparse.get(j)):
                chars.append(seed_seq[j])
            else:
                top = cnt.most_common(2)
                chars.append(_IUPAC_PAIR.get((top[0][0], top[1][0]), top[0][0]) if len(top) >= 2 and top[1][1] >= 0.30 * (top[0][1] + top[1][1]) else top[0][0])
        cons_out[sid] = "".join(chars)
    return cons_out


def _cluster_haplotypes(seqs: Mapping[str, str], prefix: str, k: int = 11) -> Dict[str, str]:
    """Cluster consensus haplotypes deterministically by 11-mer Jaccard similarity."""
    if not seqs: return {}
    samples = sorted(seqs.keys())
    ksets = {sid: {seqs[sid][i : i + k] for i in range(max(0, len(seqs[sid]) - k + 1)) if "N" not in seqs[sid][i : i + k]} for sid in samples}
    reps: List[Tuple[str, set[str]]] = []
    out: Dict[str, str] = {}
    for sid in samples:
        ks = ksets[sid]
        best_label, best_sim = max(((lbl, float(len(ks & rks)) / float(u)) for lbl, rks in reps if (u := len(ks | rks)) > 0), key=lambda x: x[1], default=(None, -1.0))
        if best_label is not None and best_sim >= 0.85:
            out[sid] = best_label
        else:
            lbl = f"{prefix}_{len(reps) + 1}"; reps.append((lbl, ks)); out[sid] = lbl
    return out


def _extract_species_from_seed_key(seed_key: str) -> str:
    norm = re.sub(r"[|_/\t]+", " ", seed_key).lower()
    for syn, canon in _SEED_SYNONYM_TAXA.items():
        if syn in norm:
            return canon
    if m := re.search(r"\bpackera\s+[0-9]*([a-z]+(?:-[a-z]+)?)(?=\s|$|[0-9])", norm):
        ep = m.group(1).strip("-")
        if ep and ep not in _NON_EPITHET_TOKENS:
            return "Packera dubia" if ep == "tomentosa" else f"Packera {ep}"
    return ""


def _match_best_public_lineage(seq: str, ref_seeds: Mapping[str, str], k: int = 13, _kcache: Optional[Dict[str, Tuple[str, set[str]]]] = None) -> str:
    """Identify the closest public Packera taxon reference using 13-mer Jaccard similarity."""
    if not seq or len(seq) < k or not ref_seeds: return ""
    q_kmers = {seq[i : i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i : i + k]}
    if not q_kmers: return ""
    best_sp, best_sim = "", 0.0
    for rkey, rseq in ref_seeds.items():
        if _kcache is not None and rkey in _kcache:
            sp, r_kmers = _kcache[rkey]
        else:
            sp = _extract_species_from_seed_key(rkey)
            step = 4 if len(rseq) > 20000 else 1
            r_kmers = {rseq[i : i + k] for i in range(0, len(rseq) - k + 1, step) if "N" not in rseq[i : i + k]} if (sp and len(rseq) >= k) else set()
            if _kcache is not None: _kcache[rkey] = (sp, r_kmers)
        if not sp or not r_kmers or not (inter := len(q_kmers & r_kmers)): continue
        if (union := len(q_kmers | r_kmers)) > 0 and (sim := float(inter) / float(union)) > best_sim:
            best_sim, best_sp = sim, sp
    return best_sp


def _detect_cytonuclear_discordance(
    plastid_clades: Mapping[str, str], its_ribotypes: Mapping[str, str], nuclear_tree_path: Optional[Path] = None,
    public_plastome_matches: Optional[Mapping[str, str]] = None, public_its_matches: Optional[Mapping[str, str]] = None, waxy_lineage_matches: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, Any]]:
    """Compare organellar clade vs. nuclear taxon/sister group to flag cytonuclear capture."""
    parse_code = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").parse_sample_code
    _ = nuclear_tree_path
    samples = sorted(set(plastid_clades.keys()) | set(its_ribotypes.keys()) | set((waxy_lineage_matches or {}).keys()))
    taxon_of = {sid: (parse_code(sid).get("species_standardized") or "Unknown") for sid in samples}
    by_clade: Dict[str, Counter[str]] = defaultdict(Counter)
    for sid, pc in plastid_clades.items():
        if taxon_of[sid].startswith("Packera "): by_clade[pc][taxon_of[sid]] += 1
    clade_maj_taxon: Dict[str, str] = {pc: cnt.most_common(1)[0][0] for pc, cnt in by_clade.items()}

    rows: List[Dict[str, Any]] = []
    for sid in samples:
        pc, sp = plastid_clades.get(sid), taxon_of[sid]
        # Flag cytonuclear capture when a Packera sample carries a chloroplast clade dominated by another species
        discordant = 1 if (pc and sp.startswith("Packera ") and pc in clade_maj_taxon and clade_maj_taxon[pc] != sp and by_clade[pc][sp] < by_clade[pc][clade_maj_taxon[pc]]) else 0
        rows.append({
            "genetic_sample_id": sid, "nuclear_taxon": sp, "plastid_clade": pc or "", "its_ribotype": its_ribotypes.get(sid) or "",
            "cytonuclear_discordant": discordant, "public_plastome_match": (public_plastome_matches or {}).get(sid, ""),
            "public_its_match": (public_its_matches or {}).get(sid, ""), "waxy_lineage_match": (waxy_lineage_matches or {}).get(sid, ""),
        })
    return rows


def run_offtarget_organellar_skim(
    fastq_source: Optional[Path], genomics_dir: Path, nuclear_tree: Optional[Path] = None,
    out_dir: Optional[Path] = None, threads: Optional[int] = None,
) -> Dict[str, Any]:
    """Skim off-target plastid, ITS, and nuclear waxy reads into `public_organellar_and_loci/` (or `--out-dir`)."""
    if fastq_source is None or not fastq_source.exists():
        msg = f"[skim-organellar] FASTQ source absent ({fastq_source}); skipping cleanly."
        logger.info(msg); print(msg)
        return {"status": "skipped_missing_fastq", "n_samples": 0}
    getorg_bin = shutil.which("get_organelle_from_reads.py")
    if getorg_bin and fastq_source.is_dir():
        logger.info("Detected GetOrganelle binary at %s; running alongside pure-Python pileup audit.", getorg_bin)
        with tempfile.TemporaryDirectory(prefix="packera_getorg_") as tmpd:
            _ = (tmpd, subprocess)

    plastid_seeds, its_seeds, waxy_seeds = _collect_organellar_seeds(genomics_dir)
    cp_kidx, its_kidx, waxy_kidx = _build_kmer_index(plastid_seeds), _build_kmer_index(its_seeds), _build_kmer_index(waxy_seeds)
    seed_kmers = set(cp_kidx) | set(its_kidx) | set(waxy_kidx)
    seen_samples: set[str] = set()
    cp_reads, its_reads, waxy_reads = defaultdict(list), defaultdict(list), defaultdict(list)

    eff_threads = _resolve_threads(threads)
    fq_files = sorted((p for p in fastq_source.rglob("*") if p.is_file() and _FQ_PAT.search(p.name)), key=lambda p: (-p.stat().st_size, p.name)) if fastq_source.is_dir() else []
    if len(fq_files) > 1 and eff_threads > 1:
        n_workers = min(eff_threads, len(fq_files))
        logger.info("[skim-organellar] Scanning %d FASTQ file(s) across %d worker(s)...", len(fq_files), n_workers)
        _WORKER_CTX.clear()
        _WORKER_CTX.update({"seed_kmers": seed_kmers, "cp_kidx": cp_kidx, "its_kidx": its_kidx, "waxy_kidx": waxy_kidx})
        use_fork = "fork" in mp.get_all_start_methods()
        ctx_mp = mp.get_context("fork") if use_fork else mp.get_context()
        pool_kw: Dict[str, Any] = {"processes": n_workers} if use_fork else {"processes": n_workers, "initializer": _init_skim_worker, "initargs": (dict(_WORKER_CTX),)}
        try:
            with ctx_mp.Pool(**pool_kw) as pool:
                for done_i, (fname, sid, n_reads, cp_b, its_b, waxy_b) in enumerate(pool.imap_unordered(_skim_file_worker, [str(p) for p in fq_files]), start=1):
                    seen_samples.add(sid)
                    for src_b, dst_b in ((cp_b, cp_reads[sid]), (its_b, its_reads[sid]), (waxy_b, waxy_reads[sid])):
                        if src_b and len(dst_b) < _MAX_BUCKET_READS:
                            dst_b.extend(src_b[: _MAX_BUCKET_READS - len(dst_b)])
                    logger.info("[skim-organellar] [%d/%d] %s (%s): scanned %d reads (cp=%d, its=%d, waxy=%d)", done_i, len(fq_files), fname, sid, n_reads, len(cp_b), len(its_b), len(waxy_b))
        finally:
            _WORKER_CTX.clear()
    else:
        for sid, seq in _iter_staged_fastq_reads(fastq_source):
            seen_samples.add(sid)
            cp_b, its_b, waxy_b = cp_reads[sid], its_reads[sid], waxy_reads[sid]
            if len(cp_b) < _MAX_BUCKET_READS or len(its_b) < _MAX_BUCKET_READS or len(waxy_b) < _MAX_BUCKET_READS:
                _bait_stream_into_buckets((seq,), seed_kmers, cp_kidx, its_kidx, waxy_kidx, cp_b, its_b, waxy_b)

    if not seen_samples:
        logger.info("[skim-organellar] No FASTQ reads discovered inside %s; skipping.", fastq_source)
        return {"status": "skipped_no_reads", "n_samples": 0}

    cp_cons = _assemble_baited_consensus(cp_reads, plastid_seeds, kidx=cp_kidx)
    its_cons = _assemble_baited_consensus(its_reads, its_seeds, kidx=its_kidx)
    waxy_cons = _assemble_baited_consensus(waxy_reads, waxy_seeds, kidx=waxy_kidx)
    cp_clades, its_ribo, waxy_clades = _cluster_haplotypes(cp_cons, prefix="cp_Clade"), _cluster_haplotypes(its_cons, prefix="ITS_Ribotype"), _cluster_haplotypes(waxy_cons, prefix="Waxy_Lineage")
    cp_cache, its_cache, waxy_cache = {}, {}, {}
    pub_cp_match = {sid: _match_best_public_lineage(seq, plastid_seeds, _kcache=cp_cache) for sid, seq in cp_cons.items()}
    pub_its_match = {sid: _match_best_public_lineage(seq, its_seeds, _kcache=its_cache) for sid, seq in its_cons.items()}
    waxy_match = {sid: _match_best_public_lineage(seq, waxy_seeds, _kcache=waxy_cache) or waxy_clades.get(sid, "") for sid, seq in waxy_cons.items()}
    cyto_rows = _detect_cytonuclear_discordance(cp_clades, its_ribo, nuclear_tree_path=nuclear_tree, public_plastome_matches=pub_cp_match, public_its_matches=pub_its_match, waxy_lineage_matches=waxy_match)

    target_dir = out_dir if out_dir is not None else (genomics_dir if genomics_dir.name == "public_organellar_and_loci" else (genomics_dir / "public_organellar_and_loci"))
    target_dir.mkdir(parents=True, exist_ok=True)
    cp_fasta, its_fasta, waxy_fasta = target_dir / "skimmed_plastid_consensus.fasta", target_dir / "skimmed_its_consensus.fasta", target_dir / "skimmed_waxy_consensus.fasta"
    cyto_csv = target_dir / "cytonuclear_discordance.csv"
    for out_fa, cons_map, tag, lbl_map, fallback in (
        (cp_fasta, cp_cons, "clade", cp_clades, "cp_Unassigned"),
        (its_fasta, its_cons, "ribotype", its_ribo, "ITS_Unassigned"),
        (waxy_fasta, waxy_cons, "waxy", waxy_clades, "Waxy_Unassigned"),
    ):
        with open(out_fa, "wt", encoding="utf-8") as fh:
            for sid, seq in sorted(cons_map.items()):
                fh.write(f">{sid}|{tag}={lbl_map.get(sid, fallback)}\n{seq}\n")
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


if __name__ == "__main__":
    raise SystemExit(main())

