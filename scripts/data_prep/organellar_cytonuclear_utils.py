#!/usr/bin/env python3
"""
Module: scripts/data_prep/organellar_cytonuclear_utils.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
  K-mer public reference lineage matching, haplotype clustering, and
  empirical-nuclear cytonuclear discordance detection for skimmed plastid,
  nrITS, and nuclear waxy sequences.
"""

from __future__ import annotations

import importlib
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

from scripts.core.config import resolve_target_genus
from scripts.data_prep.step06_synonym_utils import (
    CANONICAL_SENECIO_PACKERA_SYNONYMS as _SEED_SYNONYM_TAXA,
    _NON_EPITHET_TOKENS,
    canonical_taxon_equiv_key,
    extract_species_from_seed_key as _extract_species_from_seed_key,
)


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _canon_packera_taxon_key(taxon: str) -> str:
    return canonical_taxon_equiv_key(taxon)


def _is_valid_packera_binomial(taxon: str) -> bool:
    t = str(taxon or "").strip()
    genus = (resolve_target_genus() or "Packera").strip().capitalize()
    if not t.startswith(f"{genus} ") or len(t.split()) < 2:
        return False
    g_low = genus.lower()
    return _canon_packera_taxon_key(t) not in (f"{g_low} sp", f"{g_low} sp.", f"{g_low} unknown")


def _cluster_haplotypes(
    seqs: Mapping[str, str], prefix: str, k: int = 13,
    pos_maps: Optional[Mapping[str, Mapping[int, str]]] = None, min_id: float = 0.990,
) -> Dict[str, str]:
    """Cluster consensus haplotypes deterministically by shared-site identity or 13-mer overlap."""
    if not seqs:
        return {}
    samples = sorted(seqs.keys(), key=lambda s: (-sum(c != "N" for c in seqs[s]), s))
    ksets = {sid: {seqs[sid][i : i + k] for i in range(max(0, len(seqs[sid]) - k + 1)) if "N" not in seqs[sid][i : i + k]} for sid in samples}
    reps: List[Tuple[str, set[str], Mapping[int, str]]] = []
    out: Dict[str, str] = {}
    for sid in samples:
        ks, spos = ksets[sid], (pos_maps.get(sid, {}) if pos_maps is not None else {})
        best_label, best_sim = None, -1.0
        for lbl, rks, rpos in reps:
            shared_p = (spos.keys() & rpos.keys()) if pos_maps is not None else set()
            if len(shared_p) >= 150:
                ident = sum(1 for p in shared_p if spos[p] == rpos[p]) / float(len(shared_p))
                if ident >= min_id and ident > best_sim:
                    best_label, best_sim = lbl, ident
            elif min(len(ks), len(rks)) >= 30 and (ov := len(ks & rks) / float(min(len(ks), len(rks)))) >= 0.88 and ov > best_sim:
                best_label, best_sim = lbl, ov
            if best_label is None and prefix != "cp_Clade" and min(len(ks), len(rks)) >= 30 and (ov := len(ks & rks) / float(min(len(ks), len(rks)))) >= 0.94 and ov > best_sim:
                best_label, best_sim = lbl, ov
        if best_label is not None:
            out[sid] = best_label
        else:
            lbl = f"{prefix}_{len(reps) + 1}"
            reps.append((lbl, ks, spos))
            out[sid] = lbl
    return out


def match_best_public_lineage_with_metrics(
    seq: str, ref_seeds: Mapping[str, str], k: int = 13,
    _kcache: Optional[Dict[str, Any]] = None, q_pos: Optional[Mapping[int, str]] = None,
) -> Tuple[str, float, int]:
    """Identify closest public Packera taxon reference returning (taxon, kmer_containment, n_shared)."""
    _ = q_pos
    if not ref_seeds or not seq:
        return "", 0.0, 0
    cache = _kcache if _kcache is not None else {}
    max_ref_len = max((len(s) for s in ref_seeds.values()), default=0)
    use_k = 17 if (len(seq) > 5000 and max_ref_len >= 20000) else k
    if len(seq) < use_k:
        return "", 0.0, 0
    q_kmers = {seq[i : i + use_k] for i in range(len(seq) - use_k + 1) if "N" not in seq[i : i + use_k]}
    min_i = 15 if len(q_kmers) >= 30 else max(3, len(q_kmers) // 3)
    if not q_kmers:
        return "", 0.0, 0
    if "__refs__" not in cache:
        ref_list, sp_kmers = [], defaultdict(set)
        for rkey, rseq in ref_seeds.items():
            if not (sp := _extract_species_from_seed_key(rkey)) or len(rseq) < use_k:
                continue
            rc = _revcomp(rseq)
            rk = {rseq[i : i + use_k] for i in range(len(rseq) - use_k + 1) if "N" not in rseq[i : i + use_k]}
            rk |= {rc[i : i + use_k] for i in range(len(rc) - use_k + 1) if "N" not in rc[i : i + use_k]}
            if rk:
                low_rk = rkey.lower()
                is_pure_its = (180 <= len(rseq) <= 670) and ("large|subunit" not in low_rk) and ("external" not in low_rk)
                ref_list.append((sp, rk, len(rseq), is_pure_its))
                sp_kmers[sp].update(rk)
        cache["__refs__"] = (ref_list, Counter(km for rset in sp_kmers.values() for km in rset))
    ref_list, df_cnt = cache["__refs__"]
    is_plastome_q = len(seq) > 5000 and max_ref_len >= 20000
    complete_refs = [(sp, rk, rl, p_its) for sp, rk, rl, p_its in ref_list if rl >= 125000] if is_plastome_q else []
    pure_its_refs = [(sp, rk, rl, p_its) for sp, rk, rl, p_its in ref_list if p_its] if (not is_plastome_q and len(seq) <= 5000) else []
    pools = ((complete_refs, 0.20), (ref_list, 0.35)) if complete_refs else (((pure_its_refs, 0.35), (ref_list, 0.35)) if pure_its_refs else ((ref_list, 0.35),))
    for pool, min_cov in pools:
        best_sp, best_sc, best_cont, best_n = "", (0.0, 0.0, 0.0), 0.0, 0
        for sp, r_kmers, rl, _ in pool:
            if len(shared := (q_kmers & r_kmers)) < min_i:
                continue
            denom = float(len(q_kmers) if is_plastome_q else min(len(q_kmers), len(r_kmers)))
            if (cov_base := len(shared) / denom) < min_cov:
                continue
            cov_key = round(cov_base, 4) if is_plastome_q else cov_base
            idf_sc = sum(1.0 / float(df_cnt[km]) for km in shared if (not is_plastome_q or df_cnt[km] < 15)) / denom
            if (sc := (cov_key, idf_sc, float(rl >= 500))) > best_sc:
                best_sc, best_sp, best_cont, best_n = sc, sp, round(cov_base, 4), len(shared)
        if best_sp:
            return best_sp, best_cont, best_n
    return "", 0.0, 0


def _match_best_public_lineage(
    seq: str, ref_seeds: Mapping[str, str], k: int = 13,
    _kcache: Optional[Dict[str, Any]] = None, q_pos: Optional[Mapping[int, str]] = None,
) -> str:
    """Identify the closest public Packera taxon reference using strand-aware k-mer recall and IDF specificity."""
    return match_best_public_lineage_with_metrics(seq, ref_seeds, k=k, _kcache=_kcache, q_pos=q_pos)[0]


def _infer_tree_neighbor_taxa(nuclear_tree_path: Optional[Path], raw_sp_of: Mapping[str, str]) -> Dict[str, str]:
    if nuclear_tree_path is None or not Path(nuclear_tree_path).is_file(): return {}
    txt = Path(nuclear_tree_path).read_text(encoding="utf-8", errors="ignore")
    clean_fn = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").clean_sample_id
    tips = [clean_fn(t) for t in re.findall(r"([A-Za-z0-9_.-]+)\s*:[0-9.eE+-]+", txt)]
    out: Dict[str, str] = {}
    for idx in range(1, len(tips) - 1):
        left_sp, right_sp = raw_sp_of.get(tips[idx - 1], ""), raw_sp_of.get(tips[idx + 1], "")
        if _is_valid_packera_binomial(left_sp) and _canon_packera_taxon_key(left_sp) == _canon_packera_taxon_key(right_sp):
            out[tips[idx]] = left_sp
    return out


def _detect_cytonuclear_discordance(
    plastid_clades: Mapping[str, str], its_ribotypes: Mapping[str, str], nuclear_tree_path: Optional[Path] = None,
    public_plastome_matches: Optional[Mapping[str, str]] = None, public_its_matches: Optional[Mapping[str, str]] = None,
    waxy_lineage_matches: Optional[Mapping[str, str]] = None, plastome_match_metrics: Optional[Mapping[str, Tuple[float, int]]] = None,
    nuclear_taxon_map: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, Any]]:
    """Compare organellar clade vs. empirical nuclear taxon and quality-gated reference matches."""
    parse_code = importlib.import_module("scripts.data_prep.06_ingest_genetic_data").parse_sample_code
    samples = sorted(set(plastid_clades) | set(its_ribotypes) | set(waxy_lineage_matches or {}) | set(public_plastome_matches or {}))
    raw_sp_of = {sid: (parse_code(sid).get("species_standardized") or "Unknown") for sid in samples}
    tree_inferred = _infer_tree_neighbor_taxa(nuclear_tree_path, raw_sp_of)
    its_clusters: Dict[str, List[str]] = defaultdict(list)
    for sid, ilbl in its_ribotypes.items():
        if ilbl: its_clusters[ilbl].append(sid)

    genus_prefix = f"{(resolve_target_genus() or 'Packera').strip().capitalize()} "
    taxon_of: Dict[str, str] = {}
    for sid in samples:
        raw_sp = raw_sp_of[sid]
        if not raw_sp.startswith(genus_prefix) and raw_sp != "Unknown":
            taxon_of[sid] = raw_sp; continue
        emp_tx = str((nuclear_taxon_map or {}).get(sid) or "").strip()
        pub_its, pub_wx = (public_its_matches or {}).get(sid, ""), (waxy_lineage_matches or {}).get(sid, "")
        if _is_valid_packera_binomial(emp_tx):
            taxon_of[sid] = emp_tx
        elif _is_valid_packera_binomial(pub_its) and _is_valid_packera_binomial(pub_wx) and _canon_packera_taxon_key(pub_its) == _canon_packera_taxon_key(pub_wx):
            taxon_of[sid] = pub_its
        else:
            others = [raw_sp_of.get(o, "") for o in its_clusters.get(its_ribotypes.get(sid, ""), ()) if o != sid and _is_valid_packera_binomial(raw_sp_of.get(o, ""))]
            if len(its_clusters) >= 2 and len(others) >= 2 and len({_canon_packera_taxon_key(x) for x in others}) == 1 and _canon_packera_taxon_key(others[0]) != _canon_packera_taxon_key(raw_sp) and not (_is_valid_packera_binomial(pub_wx) and _canon_packera_taxon_key(pub_wx) == _canon_packera_taxon_key(raw_sp)):
                taxon_of[sid] = others[0]
            elif sid in tree_inferred and _canon_packera_taxon_key(tree_inferred[sid]) != _canon_packera_taxon_key(raw_sp):
                taxon_of[sid] = tree_inferred[sid]
            else:
                taxon_of[sid] = raw_sp

    by_clade: Dict[str, Counter[str]] = defaultdict(Counter)
    for sid, pc in plastid_clades.items():
        if taxon_of.get(sid, "").startswith(genus_prefix): by_clade[pc][taxon_of[sid]] += 1
    clade_maj_taxon: Dict[str, str] = {pc: cnt.most_common(1)[0][0] for pc, cnt in by_clade.items()}

    rows: List[Dict[str, Any]] = []
    for sid in samples:
        pc, sp = plastid_clades.get(sid), taxon_of[sid]
        pub_cp, pub_its, pub_wx = (public_plastome_matches or {}).get(sid, ""), (public_its_matches or {}).get(sid, ""), (waxy_lineage_matches or {}).get(sid, "")
        discordant = 1 if (pc and sp.startswith(genus_prefix) and pc in clade_maj_taxon and clade_maj_taxon[pc] != sp and by_clade[pc][sp] < by_clade[pc][clade_maj_taxon[pc]]) else 0
        if discordant == 0 and sp.startswith(genus_prefix) and _is_valid_packera_binomial(pub_cp):
            cp_key, sp_key = _canon_packera_taxon_key(pub_cp), _canon_packera_taxon_key(sp)
            metrics = (plastome_match_metrics or {}).get(sid)
            hq_cp = metrics is not None and float(metrics[0]) >= 0.80 and int(metrics[1]) >= 80
            nuc_matches_cp = (_canon_packera_taxon_key(pub_its) == cp_key) or (_canon_packera_taxon_key(pub_wx) == cp_key)
            if cp_key != sp_key and hq_cp and not nuc_matches_cp: discordant = 1
        rows.append({
            "genetic_sample_id": sid, "nuclear_taxon": sp, "plastid_clade": pc or "",
            "its_ribotype": its_ribotypes.get(sid) or "", "cytonuclear_discordant": discordant,
            "public_plastome_match": pub_cp, "public_its_match": pub_its, "waxy_lineage_match": pub_wx,
        })
    return rows
