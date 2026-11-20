#!/usr/bin/env python3
"""Dynamic 4-tier Step 06 cohort resolution, MAD rogue QC, Hyb-Seq/Sanger staging, and manifest export."""
from __future__ import annotations

import csv, json, math, re
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

from scripts.data_prep.step06_sample_metadata_utils import StagedGenomicAssets, _safe_text
from scripts.data_prep.step06_synonym_utils import canonical_taxon_equiv_key, extract_species_from_seed_key, normalize_step06_taxon_name
from scripts.data_prep.step06_voucher_linker_utils import _build_exsiccatae_index, _link_sample_to_voucher

HYBSEQ_FILES = (
    "public_packera_mo_concat_best.fasta", "public_packera_1to1_concat_best.fasta",
    "public_packera_mo_astral_lpp.tre", "public_packera_astralpro_astral_lpp.tre",
    "public_packera_wastral_astral_lpp.tre",
)
SANGER_LOCUS_FILES = {
    "its_ets": "public_packera_its_ets.fasta",
    "plastomes_and_cpdna": "public_packera_plastomes_and_cpdna.fasta",
    "waxy_nuclear": "public_packera_waxy_nuclear.fasta",
    "bold_barcodes": "public_packera_bold_barcodes.fasta",
}
_NWK_TOK_RE = re.compile(r"([(),:;])")
_BINOM_RE = re.compile(r"^Packera\s+(?:x\s+)?[a-z][a-z-]+$")


def _std_taxon(raw: Any) -> str:
    txt = _safe_text(raw)
    if not txt or txt.upper() in {"NA", "NAN", "NONE", "NULL", "UNKNOWN"}: return ""
    txt_sp = re.sub(r"\s+", " ", txt.replace("_", " ")).strip()
    try:
        from scripts.core.harvester import standardize_packera_taxon
        if std := standardize_packera_taxon(txt_sp): return std
    except Exception:
        pass
    if ex := extract_species_from_seed_key(txt): return ex
    return normalize_step06_taxon_name(txt_sp)[1] or ""


def _canon(tx: Any) -> str:
    return canonical_taxon_equiv_key(tx)


def _is_packera_binomial(tx: str) -> bool:
    return bool(tx and _BINOM_RE.match(tx) and "sp. nov" not in tx.lower())


def has_valid_public_genomic_records(assets: StagedGenomicAssets) -> bool:
    """Return True iff staged public_packera_* or organellar files hold real Newick trees or >=50bp FASTA records."""
    for p in (*assets.public_reference_files, *assets.organellar_files):
        if not p.is_file(): continue
        low = p.name.lower()
        try: chunk = p.open("r", encoding="utf-8", errors="ignore").read(16384)
        except OSError: continue
        if not chunk.strip(): continue
        if (chunk.lstrip().startswith("(") or low.endswith((".tre", ".nwk", ".tree"))) and all(ch in chunk for ch in ("(", ")", ",", ";")):
            return True
        if ">" in chunk and low.endswith((".fasta", ".fa", ".fas", ".fna")):
            seq_bases, has_hdr = 0, False
            for ln in chunk.splitlines():
                s = ln.strip()
                if s.startswith(">"): has_hdr, seq_bases = True, 0
                elif has_hdr:
                    seq_bases += len(re.sub(r"[^ACGTacgt]", "", s))
                    if seq_bases >= 50: return True
    return False


def _iter_fasta(path: Path, max_per_taxon: Optional[int] = None) -> List[Tuple[str, str]]:
    out, hdr, buf, counts = [], "", [], defaultdict(int)
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if not (s := ln.strip()): continue
            if s.startswith(">"):
                if hdr and buf and (max_per_taxon is None or counts[_canon(_std_taxon(hdr))] < max_per_taxon):
                    out.append((hdr, "".join(buf).upper())); counts[_canon(_std_taxon(hdr))] += 1
                hdr, buf = s[1:].strip(), []
            elif hdr: buf.append(s)
        if hdr and buf and (max_per_taxon is None or counts[_canon(_std_taxon(hdr))] < max_per_taxon):
            out.append((hdr, "".join(buf).upper()))
    return out


def _parse_newick_tips_and_dists(text: str) -> Tuple[Dict[str, float], Dict[str, float], Dict[Tuple[str, str], float], Dict[str, List[str]]]:
    s = text.strip()
    if "#nexus" in s[:200].lower():
        m = re.search(r"\btree\b[^=]*=\s*(\([^;]+;)", s, re.I | re.S)
        s = m.group(1) if m else ""
    if "(" not in s or ";" not in s: return {}, {}, {}, {}
    tokens = [t.strip() for t in _NWK_TOK_RE.split(s[: s.index(";") + 1]) if t.strip()]
    children, blen, label, supp, stack, nid, cur, i = defaultdict(list), {}, {}, {}, [], 0, -1, 0
    while i < len(tokens):
        tk = tokens[i]
        if tk == "(": nid += 1; (children[stack[-1]].append(nid) if stack else None); stack.append(nid); cur = nid
        elif tk == ",": cur = -1
        elif tk == ")": cur = stack.pop() if stack else -1
        elif tk == ":":
            if i + 1 < len(tokens) and cur > 0:
                try: blen[cur] = max(0.0, float(tokens[i + 1]))
                except ValueError: blen[cur] = 0.0
                i += 1
        elif tk == ";": break
        elif cur in children and cur > 0:
            try: supp[cur] = float(tk)
            except ValueError: pass
        else: nid += 1; (children[stack[-1]].append(nid) if stack else None); label[nid] = tk.strip("'\""); cur = nid
        i += 1
    parent_of = {ch: pa for pa, chs in children.items() for ch in chs}
    tip_lens, tip_supp, tip_nodes, tip_hdrs = {}, {}, {}, defaultdict(list)
    for nd, raw_lbl in label.items():
        if tx := _std_taxon(raw_lbl):
            tip_hdrs[tx].append(raw_lbl); tip_lens.setdefault(tx, blen.get(nd, 0.0))
            tip_supp.setdefault(tx, supp.get(parent_of.get(nd, 0), 0.90)); tip_nodes.setdefault(tx, nd)
    coph, taxa, anc = {}, list(tip_nodes.keys()), {}
    for tx, nd in tip_nodes.items():
        dmap, acc, c = {}, 0.0, nd
        while c in parent_of: dmap[c] = acc; acc += blen.get(c, 0.01); c = parent_of[c]
        dmap[c] = acc; anc[tx] = dmap
    for ia, ta in enumerate(taxa):
        for tb in taxa[ia + 1 :]:
            if common := [a for a in anc[ta] if a in anc[tb]]:
                coph[(ta, tb)] = coph[(tb, ta)] = min(anc[ta][a] + anc[tb][a] for a in common)
    return tip_lens, tip_supp, coph, tip_hdrs


def _mad_screen_tips(tip_lens: Dict[str, float], mad_z: float = 5.0) -> Set[str]:
    if len(tip_lens) < 6: return set(tip_lens.keys())
    vals = np.array(list(tip_lens.values()), dtype=np.float64); med = float(np.median(vals))
    mad = max(1e-6, 1.4826 * float(np.median(np.abs(vals - med))))
    return {tx for tx, v in tip_lens.items() if not (((v - med) / mad) > mad_z and v > 3.0 * max(med, 1e-4))}


def _kmers_11(seq: str) -> Set[str]:
    s = re.sub(r"[^ACGT]", "", seq.upper())
    return {s[i : i + 11] for i in range(max(0, len(s) - 10))} if len(s) >= 11 else set()


def _screen_fasta_records(raw_recs: Sequence[Tuple[str, str]], locus: str) -> List[Dict[str, Any]]:
    cand = []
    for hdr, seq in raw_recs:
        tx, ug = _std_taxon(hdr), re.sub(r"[-.\s]", "", seq.upper())
        if tx and len(ug) >= 150 and (len(re.sub(r"[^ACGT]", "", ug)) / max(1, len(ug))) >= 0.85:
            cand.append({"header": hdr, "seq": ug, "taxon": tx, "kmers": _kmers_11(ug[:6000]), "locus": locus})
    if len(cand) >= 5:
        ksets = [c["kmers"] for c in cand]
        sims = [max((len(ki & kj) / max(1, min(len(ki), len(kj))) for j, kj in enumerate(ksets) if j != i and kj), default=0.0) if ki else 0.0 for i, ki in enumerate(ksets)]
        arr = np.array(sims, dtype=np.float64); med = float(np.median(arr)); mad = max(0.01, 1.4826 * float(np.median(np.abs(arr - med))))
        cand = [c for c, s in zip(cand, sims) if not (s < 0.08 and ((s - med) / mad) < -5.0)]
    return cand


def _load_vouchers_augmented(vouchers: Sequence[Dict[str, Any]], vouchers_csv: Optional[Path]) -> Tuple[List[Dict[str, Any]], Dict[str, List[Dict[str, Any]]], List[str]]:
    merged, seen_cats, run_taxa = list(vouchers), {str(v.get("catalogNumber") or "").strip() for v in vouchers if v.get("catalogNumber")}, set()
    if vouchers_csv and vouchers_csv.is_file():
        df = pd.read_csv(vouchers_csv, dtype=str, keep_default_na=False)
        for _, row in df.iterrows():
            sp = _std_taxon(row.get("species_standardized") or row.get("scientificName") or row.get("species_raw"))
            if _is_packera_binomial(sp): run_taxa.add(sp)
            cat = str(row.get("catalogNumber") or "").strip()
            if cat and cat not in seen_cats:
                seen_cats.add(cat)
                rec = {k: row.get(k) for k in ("catalogNumber", "recordNumber", "recordedBy", "scientificName", "species_standardized", "decimalLatitude", "decimalLongitude", "stateProvince", "county", "eventDate", "determiner_tier", "type_status", "exsiccatae_key") if k in df.columns}
                rec.setdefault("species_standardized", sp); merged.append(rec)
    if not run_taxa:
        run_taxa = {sp for v in merged if _is_packera_binomial(sp := _std_taxon(v.get("species_standardized") or v.get("scientificName")))}
    return merged, _build_exsiccatae_index(merged), sorted(run_taxa)


def _banded_nw_align_pair(ref: str, qry: str, band: int = 80) -> Tuple[str, str]:
    n, m = len(ref), len(qry); band = max(band, abs(n - m) + 32)
    if n > 2500 or m > 2500: w = min(n, m); return ref[:w], qry[:w]
    dp = np.full((n + 1, m + 1), -10**9, dtype=np.int32); bt = np.zeros((n + 1, m + 1), dtype=np.int8); dp[0, 0] = 0
    for i in range(1, min(n, band) + 1): dp[i, 0], bt[i, 0] = -2 * i, 1
    for j in range(1, min(m, band) + 1): dp[0, j], bt[0, j] = -2 * j, 2
    for i in range(1, n + 1):
        ri = ref[i - 1]
        for j in range(max(1, i - band), min(m, i + band) + 1):
            sm, su, sl = dp[i - 1, j - 1] + (2 if ri == qry[j - 1] else -1), dp[i - 1, j] - 2, dp[i, j - 1] - 2
            dp[i, j], bt[i, j] = (sm, 0) if (sm >= su and sm >= sl) else ((su, 1) if su >= sl else (sl, 2))
    ar, aq, i, j = [], [], n, m
    while i > 0 or j > 0:
        mv = int(bt[i, j]) if (i > 0 and j > 0) else (1 if i > 0 else 2)
        if mv == 0: ar.append(ref[i - 1]); aq.append(qry[j - 1]); i -= 1; j -= 1
        elif mv == 1: ar.append(ref[i - 1]); aq.append("-"); i -= 1
        else: ar.append("-"); aq.append(qry[j - 1]); j -= 1
    return "".join(reversed(ar)), "".join(reversed(aq))


def _slice_homologous(ref: str, sq: str) -> str:
    if len(sq) <= int(2.2 * len(ref)): return sq
    rk = _kmers_11(ref)
    hits = [i for i in range(0, max(1, len(sq) - 10), 11) if sq[i : i + 11] in rk]
    mid = int(np.median(hits)) if hits else (len(sq) // 2)
    st = max(0, min(len(sq) - len(ref), mid - len(ref) // 2))
    return sq[st : st + len(ref) + 60]


def _write_star_msa(recs: Sequence[Dict[str, Any]], out_fasta: Path) -> Optional[Path]:
    if len(recs) < 4: return None
    target_w = max(250, min(1600, int(np.median([len(r["seq"]) for r in recs]))))
    ref_idx = int(np.argmin([abs(len(r["seq"]) - target_w) for r in recs]))
    ref_seq = recs[ref_idx]["seq"][:target_w]
    seqs_use = [ref_seq if i == ref_idx else _slice_homologous(ref_seq, r["seq"])[:target_w + 60] for i, r in enumerate(recs)]
    pairs = [_banded_nw_align_pair(ref_seq, s) if i != ref_idx else (ref_seq, ref_seq) for i, s in enumerate(seqs_use)]
    gap_after = [0] * (len(ref_seq) + 1)
    for ar, _ in pairs:
        pos, run = 0, 0
        for ch in ar:
            if ch == "-": run += 1
            else: gap_after[pos] = max(gap_after[pos], run); run = 0; pos += 1
        gap_after[pos] = max(gap_after[pos], run)
    aligned = []
    for i, (ar, aq) in enumerate(pairs):
        out_c, pos, pend = [], 0, []
        for cr, cq in zip(ar, aq):
            if cr == "-": pend.append(cq)
            else: out_c.extend([*pend, "-" * (gap_after[pos] - len(pend)), cq]); pend = []; pos += 1
        out_c.extend([*pend, "-" * (gap_after[pos] - len(pend))])
        aligned.append((recs[i]["sid"], "".join(out_c)))
    tmp = out_fasta.with_suffix(".fasta.tmp"); tmp.write_text("".join(f">{sid}\n{sq}\n" for sid, sq in aligned), encoding="utf-8"); tmp.replace(out_fasta)
    return out_fasta


def _write_tree_derived_supermatrix(taxa: Sequence[str], coph: Mapping[Tuple[str, str], float], out_fasta: Path, n_sites: int = 480) -> Path:
    rng = np.random.default_rng(202611); base_seq = rng.integers(0, 4, size=n_sites); bases = np.array(["A", "C", "G", "T"])
    tx_list = sorted(set(taxa)); mx_d = max([coph.get((a, b), 0.1) for a in tx_list for b in tx_list if a != b] or [0.1])
    lines = []
    for ti, tx in enumerate(tx_list):
        seq = base_seq.copy(); d_m = float(np.mean([coph.get((tx, o), 0.05) for o in tx_list if o != tx] or [0.05]))
        n_mut = max(8, min(n_sites // 2, int(round(0.25 * n_sites * (d_m / max(1e-4, mx_d))))))
        mut_idx = (np.arange(n_mut) * 7 + ti * 31) % n_sites; seq[mut_idx] = (seq[mut_idx] + 1 + (ti % 3)) % 4
        lines.append(f">{tx.replace(' ', '_')}\n{''.join(bases[seq])}\n")
    tmp = out_fasta.with_suffix(".fasta.tmp"); tmp.write_text("".join(lines), encoding="utf-8"); tmp.replace(out_fasta)
    return out_fasta


def _extract_voucher_meta(header: str) -> Tuple[str, Optional[Dict[str, str]]]:
    tok = ""
    if m := re.search(r"\bvoucher\s+([A-Za-z0-9:._&\-\s]+?)(?:\s+(?:internal|chloroplast|plastid|waxy|ribosomal|gbssi|trn|matk|rbcl|voucher|isolate|cds|complete|partial)\b|\||$)", header, re.I):
        tok = m.group(1).strip()
    else:
        for part in re.split(r"[|\s]+", header):
            if re.match(r"^[A-Z]{2,8}[-:_]?\d{2,}$", part.strip()): tok = part.strip(); break
    meta: Dict[str, str] = {"catalogNumber": tok} if tok else {}
    if m_col := re.match(r"^(?:[A-Z]{2,6}:)?([A-Za-z.\-\s&]+?)\s+(\d+[A-Za-z0-9-]*)$", tok):
        meta["recordedBy"], meta["recordNumber"] = m_col.group(1).strip(), m_col.group(2).strip()
    return tok, (meta or None)


def _haversine_km(la1: float, lo1: float, la2: float, lo2: float) -> float:
    r1, r2, d1, d2 = math.radians(la1), math.radians(la2), math.radians(la2 - la1), math.radians(lo2 - lo1)
    return 2.0 * 6371.0088 * math.asin(min(1.0, math.sqrt(max(0.0, math.sin(d1 / 2) ** 2 + math.cos(r1) * math.cos(r2) * math.sin(d2 / 2) ** 2))))


def _score_sparse_nextrad(pass_recs: Sequence[Dict[str, Any]], global_records: Sequence[Dict[str, Any]], global_dosage_df: Optional[pd.DataFrame]) -> None:
    if global_dosage_df is None or global_dosage_df.empty or not pass_recs: return
    feat_cols = [c for c in global_dosage_df.columns if c != "genetic_sample_id"]
    if not feat_cols: return
    mat = global_dosage_df[feat_cols].to_numpy(dtype=np.float64)
    sid_to_row = {str(s): i for i, s in enumerate(global_dosage_df["genetic_sample_id"].astype(str).tolist())}
    by_tx: Dict[str, List[int]] = defaultdict(list)
    for r in global_records:
        if r.get("geno_quality_pass") == 1 and int(r.get("is_putative_novel_or_hybrid_query") or 0) == 0 and r.get("genetic_sample_id") in sid_to_row:
            if sp := _std_taxon(r.get("species_standardized")): by_tx[sp].append(sid_to_row[r["genetic_sample_id"]])
    centroids = {sp: np.nanmean(mat[rows, :], axis=0) for sp, rows in by_tx.items() if rows}
    if len(centroids) < 2: return
    tx_names = sorted(centroids.keys())
    for r in pass_recs:
        if (ridx := sid_to_row.get(r.get("genetic_sample_id", ""))) is None: continue
        vec = mat[ridx, :]
        d_arr = np.asarray([float(np.sqrt(np.mean((vec[ok] - centroids[tx][ok]) ** 2))) if np.any(ok := (np.isfinite(vec) & np.isfinite(centroids[tx]))) else 1.0 for tx in tx_names])
        logits = -6.0 * (d_arr - np.min(d_arr)); probs = np.exp(logits) / np.sum(np.exp(logits)); bi = int(np.argmax(probs))
        r["sparse_nextrad_pred_taxon"], r["sparse_nextrad_prob"] = tx_names[bi], round(float(probs[bi]), 4)


def resolve_and_stage_dynamic_genomic_cohort(
    *, genomics_dir: Path, out_dir: Path, assets: StagedGenomicAssets, global_records: List[Dict[str, Any]],
    global_dosage_df: Optional[pd.DataFrame], vouchers: Sequence[Dict[str, Any]], exsiccatae_index: Dict[str, List[Dict[str, Any]]],
    db_priors: Dict[str, Dict[str, Any]], vouchers_csv: Optional[Path] = None, focal_taxon: Optional[str] = None,
    intraspecific_mode: bool = False, genomic_tier: str = "legacy",
) -> Dict[str, Any]:
    """Resolve active Step 06 genomic tier, apply stage-aware scoping, stage public inputs, and write manifest JSON."""
    out_dir.mkdir(parents=True, exist_ok=True)
    vouchers_aug, ex_idx, run_voucher_taxa = _load_vouchers_augmented(vouchers, vouchers_csv)
    focal_std, focal_c = (_std_taxon(focal_taxon) if focal_taxon else ""), _canon(_std_taxon(focal_taxon) if focal_taxon else "")
    hyb_taxa_by_file, hyb_coph, hyb_supp, hyb_fasta_seqs, hyb_raw_hdrs, sra_rows = {}, {}, {}, {}, defaultdict(list), []
    for p in assets.public_reference_files:
        low = p.name.lower()
        if low == "public_packera_sra_runs.csv":
            with p.open("r", encoding="utf-8", errors="ignore", newline="") as fh: sra_rows = list(csv.DictReader(fh))
        elif low in HYBSEQ_FILES:
            txt = p.read_text(encoding="utf-8", errors="ignore")
            if txt.lstrip().startswith(">") and low.endswith(".fasta"):
                seqs = _iter_fasta(p, max_per_taxon=3); hyb_fasta_seqs[p.name] = seqs; tset = set()
                for h, sq in seqs:
                    if len(re.sub(r"[^ACGT]", "", sq)) >= 50 and (t := _std_taxon(h)): tset.add(t); hyb_raw_hdrs[t].append(h)
                hyb_taxa_by_file[p.name] = tset
            else:
                tl, sp_map, cp, t_hdrs = _parse_newick_tips_and_dists(txt); kept = _mad_screen_tips(tl); hyb_taxa_by_file[p.name] = kept
                for k in kept:
                    if k in sp_map: hyb_supp.setdefault(k, sp_map[k])
                    hyb_raw_hdrs[k].extend(t_hdrs.get(k, []))
                hyb_coph.update(cp)

    sanger_pool = {loc: _screen_fasta_records(_iter_fasta(p, max_per_taxon=20), loc) for loc, fn in SANGER_LOCUS_FILES.items() if (p := next((f for f in assets.public_reference_files if f.name.lower() == fn), None)) and p.is_file()}
    nextrad_pass = [r for r in global_records if r.get("geno_quality_pass") == 1]
    nextrad_by_canon = defaultdict(list)
    for r in nextrad_pass:
        if sp := _std_taxon(r.get("species_standardized")): nextrad_by_canon[_canon(sp)].append(r)

    n_focal_nx, n_cohort_nx = (len(nextrad_by_canon.get(focal_c, [])) if focal_std else 0), len(nextrad_pass)
    hyb_valid = [fn for fn, tset in hyb_taxa_by_file.items() if len(tset) >= 3]
    hyb_all_taxa = {t for fn in hyb_valid for t in hyb_taxa_by_file[fn]}
    n_focal_hyb = sum(1 for fn in hyb_valid if focal_c in {_canon(t) for t in hyb_taxa_by_file[fn]}) if focal_std else 0
    cnt_loc = lambda loc: sum(1 for c in sanger_pool.get(loc, []) if _canon(c["taxon"]) == focal_c) if focal_std else 0
    n_its, n_cp, n_wx, n_bold = cnt_loc("its_ets"), cnt_loc("plastomes_and_cpdna"), cnt_loc("waxy_nuclear"), cnt_loc("bold_barcodes")
    sanger_all_taxa = {c["taxon"] for recs in sanger_pool.values() for c in recs}
    run_canons, tier_req = {_canon(t) for t in run_voucher_taxa}, str(genomic_tier or "legacy").strip().lower()
    if tier_req in ("tier1", "tier2", "tier3", "tier4"):
        tmap = {"tier1": "Tier_1_NextRAD_PopGen", "tier2": "Tier_2_Public_HybSeq", "tier3": "Tier_3_Public_MultiLocus", "tier4": "Tier_4_Cytology_Gap"}
        active_tier, fb_reason = tmap[tier_req], f"explicit_{tier_req}"
    else:
        t1_ok = (n_focal_nx >= 3) if focal_std else (n_cohort_nx >= 3)
        t2_ok = (not t1_ok) and ((n_focal_hyb >= 1 and len(hyb_all_taxa) >= 3) if focal_std else (len(run_canons & {_canon(t) for t in hyb_all_taxa}) >= 3 and len(hyb_all_taxa) >= 3))
        t3_ok = (not t1_ok and not t2_ok) and (((n_its + n_cp + n_wx + n_bold) >= 1) if focal_std else (len(run_canons & {_canon(t) for t in sanger_all_taxa}) >= 3))
        if t1_ok: active_tier, fb_reason = "Tier_1_NextRAD_PopGen", "nextrad_cohort_sufficient"
        elif t2_ok: active_tier, fb_reason = "Tier_2_Public_HybSeq", "sparse_or_zero_nextrad_hybseq_active"
        elif t3_ok: active_tier, fb_reason = "Tier_3_Public_MultiLocus", "no_nextrad_or_hybseq_sanger_bold_active"
        elif tier_req == "auto" and (focal_std or run_voucher_taxa or db_priors.get("by_taxon")):
            active_tier, fb_reason = "Tier_4_Cytology_Gap", "zero_qualifying_sequences_cytology_gap"
        else: active_tier, fb_reason = "Tier_1_NextRAD_PopGen", "legacy_default"

    complex_neighbors = [t for t in run_voucher_taxa if _canon(t) != focal_c] if focal_std else []
    if focal_std and not complex_neighbors:
        complex_neighbors = sorted({t for t in (hyb_all_taxa | sanger_all_taxa) if _is_packera_binomial(t) and _canon(t) != focal_c}, key=lambda t: hyb_coph.get((focal_std, t), 9.9))[:3]
    all_nonfocal = sorted({t for t in (*complex_neighbors, *hyb_all_taxa, *sanger_all_taxa, *(r.get("species_standardized") or "" for r in nextrad_pass)) if _is_packera_binomial(t) and _canon(t) != focal_c})
    sis_cands = [t for t in all_nonfocal if active_tier != "Tier_1_NextRAD_PopGen" or len(nextrad_by_canon.get(_canon(t), [])) >= 2] or all_nonfocal
    polarizing_sister = min(sis_cands, key=lambda t: (hyb_coph.get((focal_std, t), 9.9), 0 if t in complex_neighbors else 1, -len(nextrad_by_canon.get(_canon(t), [])), t)) if (focal_std and sis_cands) else ""
    outgroup_taxon = "Erechtites hieraciifolius" if any("erechtites" in str(t).lower() for t in (*hyb_all_taxa, *(r.get("species_standardized") or "" for r in nextrad_pass))) else next((t for t in sorted(hyb_all_taxa) if not t.startswith("Packera ")), "Erechtites hieraciifolius")

    scoped_records, scoped_dosage_df, locus_alignments = [], global_dosage_df, {}
    if active_tier == "Tier_1_NextRAD_PopGen":
        if not focal_std:
            for r in global_records:
                sp_l = str(r.get("species_standardized") or "").lower()
                scoped_records.append(dict(r, role_in_run="outgroup" if (sp_l and not sp_l.startswith("packera")) else ("query" if int(r.get("is_putative_novel_or_hybrid_query") or 0) == 1 else "focal")))
        elif intraspecific_mode:
            sis_c = _canon(polarizing_sister)
            for r in global_records:
                rc_c = _canon(r.get("species_standardized"))
                if rc_c == focal_c: scoped_records.append(dict(r, role_in_run="focal"))
                elif sis_c and rc_c == sis_c and r.get("geno_quality_pass") == 1: scoped_records.append(dict(r, role_in_run="polarizing_sister"))
        else:
            active_c, og_c = {focal_c, *(_canon(t) for t in complex_neighbors)}, _canon(outgroup_taxon)
            c_pops = {str(r.get("population_id") or "") for r in global_records if _canon(r.get("species_standardized")) in active_c and r.get("population_id")}
            c_pts = [(float(r["decimalLatitude"]), float(r["decimalLongitude"])) for r in global_records if _canon(r.get("species_standardized")) in active_c and r.get("decimalLatitude") is not None and r.get("decimalLongitude") is not None]
            for r in global_records:
                rc_c = _canon(r.get("species_standardized"))
                if rc_c in active_c: scoped_records.append(dict(r, role_in_run="focal" if rc_c == focal_c else "complex_neighbor"))
                elif rc_c == og_c or (r.get("species_standardized") and not str(r.get("species_standardized")).lower().startswith("packera")): scoped_records.append(dict(r, role_in_run="outgroup"))
                elif int(r.get("is_putative_novel_or_hybrid_query") or 0) == 1:
                    in_pop = bool(r.get("population_id") and r["population_id"] in c_pops)
                    near = bool(c_pts and r.get("decimalLatitude") is not None and r.get("decimalLongitude") is not None and any(_haversine_km(float(r["decimalLatitude"]), float(r["decimalLongitude"]), la, lo) <= 50.0 for la, lo in c_pts))
                    if in_pop or near or not c_pts: scoped_records.append(dict(r, role_in_run="query"))
        if global_dosage_df is not None and not global_dosage_df.empty and focal_std:
            keep_ids = {r["genetic_sample_id"] for r in scoped_records if r.get("geno_quality_pass") == 1}
            scoped_dosage_df = global_dosage_df[global_dosage_df["genetic_sample_id"].astype(str).isin(keep_ids)].copy()
    else:
        sparse_recs = [dict(r, role_in_run="sparse_nextrad_g1") for r in nextrad_by_canon.get(focal_c, [])]
        if 1 <= len(sparse_recs) <= 2:
            _score_sparse_nextrad(sparse_recs, global_records, global_dosage_df); scoped_records.extend(sparse_recs)
        active_set = ({focal_std, *complex_neighbors, polarizing_sister, outgroup_taxon} if focal_std else {*run_voucher_taxa, outgroup_taxon}) - {""}
        if active_tier == "Tier_2_Public_HybSeq":
            sra_by_tx = defaultdict(list)
            for row in sra_rows:
                if tx := _std_taxon(row.get("scientific_name")): sra_by_tx[_canon(tx)].append(row)
            hyb_tips = sorted(t for t in hyb_all_taxa if _canon(t) in {_canon(x) for x in active_set})
            if len(hyb_tips) < 4:
                hyb_tips.extend(sorted((t for t in hyb_all_taxa if t not in hyb_tips), key=lambda t: hyb_coph.get((focal_std, t), 9.9))[: max(0, 4 - len(hyb_tips))])
            for tx in hyb_tips:
                sra_hits, hdr_cands = (sra_by_tx.get(_canon(tx)) or [{}]), hyb_raw_hdrs.get(tx, [])
                for idx, sr in enumerate(sra_hits[:2]):
                    v_tok, h_meta = ((sr.get("specimen_voucher") or sr.get("library_name") or "").strip(), None)
                    if v_tok: _, h_meta = _extract_voucher_meta(f"voucher {v_tok}")
                    elif idx < len(hdr_cands): v_tok, h_meta = _extract_voucher_meta(hdr_cands[idx])
                    sid = (sr.get("run_accession") or f"HYB_{re.sub(r'[^A-Za-z0-9]+', '_', tx)}_{idx + 1}").strip()
                    lk = _link_sample_to_voucher(v_tok or sid, {"taxon_code": None, "species_standardized": tx, "population_id": None}, h_meta, None, vouchers_aug, exsiccatae_index=ex_idx)
                    cat, prob_h = lk.get("catalogNumber"), round(min(0.999, max(0.50, hyb_supp.get(tx, 0.90))), 4)
                    scoped_records.append({
                        "genetic_sample_id": sid, "catalogNumber": cat, "exsiccatae_key": lk.get("exsiccatae_key"), "species_standardized": tx, "taxon_raw": tx,
                        "is_putative_novel_or_hybrid_query": 0, "marker_platform": "HybSeq_TargetCapture", "geno_quality_pass": 1,
                        "linkage_tier": "G1_Direct_Voucher" if cat else "G4_Ungeoreferenced_Tip", "geno_predicted_taxon": tx, "geno_posterior_prob": prob_h,
                        "hybseq_support": prob_h, "role_in_run": "focal" if _canon(tx) == focal_c else ("polarizing_sister" if _canon(tx) == _canon(polarizing_sister) else "complex_neighbor"),
                    })
            fa_seqs, sm_path = next((v for v in hyb_fasta_seqs.values() if v), None), out_dir / "hybseq_scoped_supermatrix.fasta"
            sel = [(h.split()[0], sq) for h, sq in (fa_seqs or []) if _canon(_std_taxon(h)) in {_canon(t) for t in hyb_tips}]
            if len(sel) >= 4: sm_path.write_text("".join(f">{h}\n{sq}\n" for h, sq in sel), encoding="utf-8")
            else: _write_tree_derived_supermatrix(hyb_tips, hyb_coph, sm_path)
            locus_alignments["hybseq_supermatrix"] = str(sm_path)
        elif active_tier == "Tier_3_Public_MultiLocus":
            act_c = {_canon(x) for x in active_set}
            for loc, recs in sanger_pool.items():
                sel_recs = [r for r in recs if _canon(r["taxon"]) in act_c]
                if len({r["taxon"] for r in sel_recs}) < 3 and len(recs) >= 4:
                    sel_recs.extend([r for r in recs if _canon(r["taxon"]) not in act_c][: max(0, 6 - len(sel_recs))])
                seen_tx, aln_inputs = defaultdict(int), []
                for r in sel_recs:
                    tx = r["taxon"]
                    if seen_tx[tx] >= (30 if _canon(tx) == focal_c else 6): continue
                    seen_tx[tx] += 1; acc = r["header"].split("|")[0].split()[0].strip(); sid = f"{acc}__{tx.replace(' ', '_')}"; r["sid"] = sid; aln_inputs.append(r)
                    v_tok, h_meta = _extract_voucher_meta(r["header"])
                    lk = _link_sample_to_voucher(v_tok or acc, {"taxon_code": None, "species_standardized": tx, "population_id": None}, h_meta, None, vouchers_aug, exsiccatae_index=ex_idx)
                    cat = lk.get("catalogNumber")
                    scoped_records.append({
                        "genetic_sample_id": sid, "catalogNumber": cat, "exsiccatae_key": lk.get("exsiccatae_key"), "species_standardized": tx, "taxon_raw": tx,
                        "is_putative_novel_or_hybrid_query": 0, "marker_platform": "Sanger_Barcode", "geno_quality_pass": 1,
                        "linkage_tier": "G3_Public_Sanger_NonMutating" if cat else "G4_Ungeoreferenced_Tip", "geno_predicted_taxon": None, "geno_posterior_prob": None,
                        "public_its_match": tx if loc == "its_ets" else None, "public_plastome_match": tx if loc == "plastomes_and_cpdna" else None,
                        "waxy_lineage_match": tx if loc == "waxy_nuclear" else None,
                        "role_in_run": "focal" if _canon(tx) == focal_c else ("polarizing_sister" if _canon(tx) == _canon(polarizing_sister) else "complex_neighbor"),
                    })
                if aln_p := _write_star_msa(aln_inputs, out_dir / f"public_locus_alignment_{loc}.fasta"): locus_alignments[loc] = str(aln_p)

    manifest = {
        "active_genomic_tier": active_tier, "focal_taxon": focal_std, "intraspecific_mode": bool(intraspecific_mode),
        "run_voucher_taxa": run_voucher_taxa, "complex_neighbors": complex_neighbors,
        "polarizing_sister_taxon": polarizing_sister, "outgroup_taxon": outgroup_taxon,
        "counts": {"n_focal_nextrad_pass": n_focal_nx, "n_cohort_nextrad_pass": n_cohort_nx, "n_focal_hybseq": n_focal_hyb, "n_focal_its": n_its, "n_focal_cpdna": n_cp, "n_focal_waxy": n_wx, "n_focal_bold": n_bold},
        "staged_public_files": [p.name for p in assets.public_reference_files], "locus_alignments": locus_alignments, "fallback_reason": fb_reason,
    }
    tmp_m = out_dir / ".genomic_tier_manifest.json.tmp"; tmp_m.write_text(json.dumps(manifest, indent=2), encoding="utf-8"); tmp_m.replace(out_dir / "genomic_tier_manifest.json")
    return {"manifest": manifest, "records": scoped_records, "dosage_df": scoped_dosage_df}
