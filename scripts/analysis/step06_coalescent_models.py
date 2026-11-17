"""BIC substitution model selection, corrected divergence, and coalescent Bayes factors."""

from __future__ import annotations

import importlib
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

_BASE_MAP: Dict[str, int] = {
    "A": 0, "C": 1, "G": 2, "T": 3, "U": 3,
    "R": 0, "Y": 1, "M": 0, "K": 2, "S": 1, "W": 0,
}
_LAG_X, _LAG_W = np.polynomial.laguerre.laggauss(16)
_LAG_LOG_W = np.log(np.maximum(1e-300, _LAG_W))


def _resolve_max_workers() -> int:
    raw = str(os.environ.get("PLANTMATICS_MAX_WORKERS") or os.environ.get("PACKERA_MAX_WORKERS") or "1").strip()
    try:
        cfg = importlib.import_module("scripts.core.config")
        fn = getattr(cfg, "resolve_max_workers", None) or getattr(cfg, "resolve_cpu_threads", None)
        return max(1, int(fn())) if callable(fn) else (max(1, int(raw)) if raw.isdigit() else 1)
    except Exception:
        return max(1, int(raw)) if raw.isdigit() else 1


def _gammap(a: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x < a + 1.0:
        ap, term, total = a, 1.0 / a, 1.0 / a
        for _ in range(200):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) <= abs(total) * 1e-13:
                break
        return min(1.0, total * math.exp(-x + a * math.log(x) - math.lgamma(a)))
    tiny, b = 1e-30, x + 1.0 - a
    c, d = 1.0 / tiny, 1.0 / max(tiny, b)
    h = d
    for i in range(1, 200):
        an = -i * (i - a)
        b += 2.0
        d = 1.0 / (tiny if abs(an * d + b) < tiny else (an * d + b))
        c = tiny if abs(b + an / c) < tiny else (b + an / c)
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-13:
            break
    return max(0.0, 1.0 - math.exp(-x + a * math.log(x) - math.lgamma(a)) * h)


def _yang_gamma4_rates(alpha: float) -> np.ndarray:
    a = min(50.0, max(0.05, float(alpha)))
    qs = [0.0]
    for p in (0.25, 0.50, 0.75):
        lo, hi = 0.0, max(2.0, a + 4.0)
        while _gammap(a, hi) < p and hi < 1e6:
            hi *= 2.0
        for _ in range(36):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if _gammap(a, mid) < p else (lo, mid)
        qs.append(0.5 * (lo + hi))
    cdfs = [0.0] + [_gammap(a + 1.0, q) for q in qs[1:]] + [1.0]
    rates = np.maximum(1e-4, 4.0 * np.diff(np.asarray(cdfs, dtype=np.float64)))
    return rates / max(1e-12, float(np.mean(rates)))


def _encode_aln(seqs: Dict[str, str], sample_ids: Sequence[str]) -> Tuple[List[str], np.ndarray, np.ndarray]:
    present = [sid for sid in sample_ids if sid in seqs and seqs[sid]]
    if len(present) < 3:
        return present, np.empty((len(present), 0), dtype=np.int8), np.full(4, 0.25, dtype=np.float64)
    n_len = min(len(seqs[sid]) for sid in present)
    raw = np.full((len(present), n_len), -1, dtype=np.int8)
    for i, sid in enumerate(present):
        arr = np.frombuffer(seqs[sid][:n_len].upper().encode("ascii", errors="replace"), dtype=np.uint8)
        for ch, idx in _BASE_MAP.items():
            raw[i, arr == ord(ch)] = idx
    mat = raw[:, np.sum(raw >= 0, axis=0) >= 3]
    if mat.shape[1] == 0:
        return present, mat, np.full(4, 0.25, dtype=np.float64)
    counts = np.bincount(mat[mat >= 0], minlength=4).astype(np.float64) + 1.0
    pi = np.maximum(0.01, counts / float(np.sum(counts)))
    return present, mat, pi / float(np.sum(pi))


def _build_nj_tree(mat: np.ndarray) -> List[Tuple[int, int, float, int, float]]:
    n_t = int(mat.shape[0])
    dist = np.zeros((2 * n_t - 1, 2 * n_t - 1), dtype=np.float64)
    ok = mat >= 0
    for i in range(n_t - 1):
        ov = ok[i][None, :] & ok[i + 1 :]
        n_ov = np.sum(ov, axis=1)
        p_diff = np.where(n_ov > 0, np.sum(ov & (mat[i][None, :] != mat[i + 1 :]), axis=1) / np.maximum(1, n_ov), 0.25)
        dist[i, i + 1 : n_t] = dist[i + 1 : n_t, i] = np.maximum(1e-4, -0.75 * np.log(np.maximum(0.05, 1.0 - (4.0 / 3.0) * np.minimum(0.74, p_diff))))
    active, next_id, nodes = list(range(n_t)), n_t, []
    while len(active) > 2:
        r = len(active)
        sub = dist[np.ix_(active, active)]
        rsum = np.sum(sub, axis=1)
        q = (r - 2.0) * sub - rsum[:, None] - rsum[None, :]
        np.fill_diagonal(q, np.inf)
        iu, iv = np.unravel_index(int(np.argmin(q)), q.shape)
        u, v, duv = active[iu], active[iv], float(sub[iu, iv])
        bl_u = max(1e-4, 0.5 * duv + (float(rsum[iu]) - float(rsum[iv])) / (2.0 * (r - 2.0)))
        nodes.append((next_id, u, bl_u, v, max(1e-4, duv - bl_u)))
        for k in active:
            if k != u and k != v:
                dist[next_id, k] = dist[k, next_id] = max(0.0, 0.5 * (float(dist[u, k]) + float(dist[v, k]) - duv))
        active = [k for k in active if k != u and k != v] + [next_id]
        next_id += 1
    bl = max(1e-4, 0.5 * float(dist[active[0], active[1]]))
    nodes.append((next_id, active[0], bl, active[1], bl))
    return nodes


def _eval_tree_logl(
    tip_masks: np.ndarray,
    pat_counts: np.ndarray,
    n_pat: int,
    tree_nodes: Sequence[Tuple[int, int, float, int, float]],
    pi: np.ndarray,
    kappa_r: float,
    kappa_y: float,
    s_tree: float,
    alpha: Optional[float],
    asc_applied: bool,
) -> float:
    q = np.ones((4, 4), dtype=np.float64) * pi[None, :]
    q[0, 2], q[2, 0], q[1, 3], q[3, 1] = kappa_r * pi[2], kappa_r * pi[0], kappa_y * pi[3], kappa_y * pi[1]
    np.fill_diagonal(q, 0.0)
    np.fill_diagonal(q, -np.sum(q, axis=1))
    q /= max(1e-12, -float(np.sum(pi * np.diag(q))))
    sq_pi = np.sqrt(pi)
    inv_sq = 1.0 / sq_pi
    evals, u_vec = np.linalg.eigh((sq_pi[:, None] * q) * inv_sq[None, :])
    rates = np.ones(1, dtype=np.float64) if alpha is None else _yang_gamma4_rates(alpha)
    g_cnt = int(rates.size)
    edge_bls = np.asarray([bl for _, _, bl_l, _, bl_r in tree_nodes for bl in (bl_l, bl_r)], dtype=np.float64)
    x_ge = (rates[:, None] * (s_tree * edge_bls)[None, :]).ravel()
    exp_lam = np.exp(np.clip(x_ge[:, None] * evals[None, :], -700.0, 0.0))
    p_mats = np.maximum(1e-15, np.einsum("ij,mj,jk->mik", inv_sq[:, None] * u_vec, exp_lam, u_vec.T * sq_pi[None, :])).reshape(g_cnt, int(edge_bls.size), 4, 4)

    n_taxa, n_tot, _ = tip_masks.shape
    cl: Dict[int, np.ndarray] = {t: np.broadcast_to(tip_masks[t][None, :, :], (g_cnt, n_tot, 4)) for t in range(n_taxa)}
    log_scale = np.zeros((g_cnt, n_tot), dtype=np.float64)
    for k_idx, (nid, lc, _, rc, _) in enumerate(tree_nodes):
        prod = np.einsum("gba,gpa->gpb", p_mats[:, 2 * k_idx], cl[lc]) * np.einsum("gba,gpa->gpb", p_mats[:, 2 * k_idx + 1], cl[rc])
        mx = np.maximum(1e-300, np.max(prod, axis=-1, keepdims=True))
        cl[nid] = prod / mx
        log_scale += np.log(mx[..., 0])
    logp_g = np.log(np.maximum(1e-300, np.sum(cl[tree_nodes[-1][0]] * pi[None, None, :], axis=-1))) + log_scale
    mx_g = np.max(logp_g, axis=0)
    logp = mx_g + np.log(np.mean(np.exp(logp_g - mx_g[None, :]), axis=0))
    logl_var = float(np.sum(pat_counts * logp[:n_pat]))
    if not asc_applied:
        return logl_var
    return logl_var - float(np.sum(pat_counts)) * math.log(max(1e-9, 1.0 - min(0.99999999, float(np.sum(np.exp(logp[n_pat:]))))))


def _optimize_log_params(obj_fn: Any, init_vals: Sequence[float], bounds: Sequence[Tuple[float, float]]) -> Tuple[List[float], float]:
    cur = [math.log(min(hi, max(lo, v))) for v, (lo, hi) in zip(init_vals, bounds)]
    best_val, gr = float(obj_fn([math.exp(z) for z in cur])), (math.sqrt(5.0) - 1.0) / 2.0
    for _ in range(2):
        for idx, (lo, hi) in enumerate(bounds):
            a, b = math.log(lo), math.log(hi)
            c, d = b - gr * (b - a), a + gr * (b - a)
            def _ev(z: float) -> float:
                vec = list(cur)
                vec[idx] = z
                return float(obj_fn([math.exp(x) for x in vec]))
            fc, fd = _ev(c), _ev(d)
            for _ in range(8):
                if fc < fd:
                    a, c, fc = c, d, fd
                    d, fd = a + gr * (b - a), _ev(a + gr * (b - a))
                else:
                    b, d, fd = d, c, fc
                    c, fc = b - gr * (b - a), _ev(b - gr * (b - a))
            z_opt = 0.5 * (a + b)
            cand = _ev(z_opt)
            if cand > best_val:
                cur[idx], best_val = z_opt, cand
    return [math.exp(z) for z in cur], best_val


def _resolve_iqtree2_bin() -> Optional[str]:
    if hit := shutil.which("iqtree2") or shutil.which("iqtree2.exe"):
        return hit
    fb = Path(os.environ.get("LOCALAPPDATA") or os.environ.get("XDG_DATA_HOME") or (Path.home() / ".local" / "share")) / "Plantmatics"
    try:
        data_dir = Path(importlib.import_module("platformdirs").user_data_dir("Plantmatics"))
    except Exception:
        data_dir = fb
    return next((str(p) for fn in ("iqtree2", "iqtree2.exe") if (p := data_dir / "bin" / fn).is_file()), None)


def select_nucleotide_substitution_model(seqs: Dict[str, str], sample_ids: Sequence[str], alignment_name: str = "packera.phylip") -> Dict[str, Any]:
    """Select best nucleotide substitution model via BIC across 10 CTMC(+G4)(+ASC) models."""
    present, mat, pi_emp = _encode_aln(seqs, sample_ids)
    n_sites = int(mat.shape[1])
    if len(present) < 3 or n_sites < 5:
        row0 = {"alignment_name": alignment_name, "model_name": "JC69", "base_model": "JC69", "has_gamma": False, "asc_applied": False, "k_params": 0, "logL": 0.0, "AICc": 0.0, "BIC": 0.0, "delta_BIC": 0.0, "bic_weight": 1.0, "ti_tv_kappa": 1.0, "gamma_alpha": None, "prop_invariant": 1.0, "n_sites": n_sites, "n_var_sites": 0}
        return {"alignment_name": alignment_name, "best_model": "JC69", "base_model": "JC69", "has_gamma": False, "asc_applied": False, "bic_weight": 1.0, "ti_tv_kappa": 1.0, "kappa_purine": 1.0, "kappa_pyrimidine": 1.0, "gamma_alpha": None, "prop_invariant": 1.0, "n_sites": n_sites, "n_var_sites": 0, "pi_freqs": [0.25] * 4, "iqtree2_model": None, "table_rows": [row0]}

    var_mask = np.max(mat, axis=0) != np.min(np.where(mat < 0, 9, mat), axis=0)
    n_var = int(np.sum(var_mask))
    prop_inv = float((n_sites - n_var) / max(1, n_sites))
    asc_applied = bool(prop_inv < 0.05 and n_sites >= 30)
    active = mat[:, var_mask] if (asc_applied and n_var >= 5) else mat
    n_eff = int(active.shape[1])
    pats, counts = np.unique(active, axis=1, return_counts=True)
    n_pat = int(pats.shape[1])
    pat_tot = np.hstack([pats, np.tile(np.arange(4, dtype=np.int8)[None, :], (len(present), 1))]) if asc_applied else pats
    tip_masks = np.ones((len(present), pat_tot.shape[1], 4), dtype=np.float64)
    for t in range(len(present)):
        called = pat_tot[t] >= 0
        if np.any(called):
            tip_masks[t, called, :] = 0.0
            tip_masks[t, called, pat_tot[t, called]] = 1.0
    tree_nodes, pi_eq = _build_nj_tree(mat), np.full(4, 0.25, dtype=np.float64)
    warm, results = {"s": 1.0, "k": 2.0, "kr": 2.0, "ky": 2.0}, []
    for base, use_k, k_base in [("JC69", False, 0), ("K80", True, 1), ("F81", False, 3), ("HKY85", True, 4), ("TN93", True, 5)]:
        for use_g in (False, True):
            pi_use = pi_eq if base in ("JC69", "K80") else pi_emp
            k_tot = k_base + (1 if use_g else 0)
            bnds = [(0.02, 6.0)] + ([(0.2, 25.0), (0.2, 25.0)] if base == "TN93" else ([(0.2, 25.0)] if use_k else [])) + ([(0.08, 25.0)] if use_g else [])
            inits = [warm["s"]] + ([warm["kr"], warm["ky"]] if base == "TN93" else ([warm["k"]] if use_k else [])) + ([1.0] if use_g else [])
            def _obj(p: Sequence[float], _pi: np.ndarray = pi_use, _b: str = base, _uk: bool = use_k, _ug: bool = use_g) -> float:
                kr, ky = (p[1], p[2]) if _b == "TN93" else ((p[1], p[1]) if _uk else (1.0, 1.0))
                return _eval_tree_logl(tip_masks, counts, n_pat, tree_nodes, _pi, kr, ky, p[0], p[-1] if _ug else None, asc_applied)
            opt, logl = _optimize_log_params(_obj, inits, bnds)
            kr, ky = (opt[1], opt[2]) if base == "TN93" else ((opt[1], opt[1]) if use_k else (1.0, 1.0))
            alp = opt[-1] if use_g else None
            if not use_g:
                warm["s"] = opt[0]
                if base in ("K80", "HKY85"):
                    warm["k"], warm["kr"], warm["ky"] = kr, kr, ky
            bic = -2.0 * logl + k_tot * math.log(max(2, n_eff))
            aicc = -2.0 * logl + 2.0 * k_tot + (2.0 * k_tot * (k_tot + 1.0)) / max(1.0, n_eff - k_tot - 1.0)
            results.append({"alignment_name": alignment_name, "model_name": base + ("+G4" if use_g else "") + ("+ASC" if asc_applied else ""), "base_model": base, "has_gamma": use_g, "asc_applied": asc_applied, "k_params": k_tot, "logL": round(logl, 4), "AICc": round(aicc, 4), "BIC": round(bic, 4), "ti_tv_kappa": round(0.5 * (kr + ky), 4), "_kr": kr, "_ky": ky, "gamma_alpha": round(alp, 4) if alp is not None else None, "prop_invariant": round(prop_inv, 4), "n_sites": n_sites, "n_var_sites": n_var})
    min_bic = min(r["BIC"] for r in results)
    w_raw = [math.exp(max(-50.0, -0.5 * (r["BIC"] - min_bic))) for r in results]
    w_sum = max(1e-12, sum(w_raw))
    for r, w in zip(results, w_raw):
        r["delta_BIC"], r["bic_weight"] = round(r["BIC"] - min_bic, 4), round(w / w_sum, 4)
    results.sort(key=lambda r: (r["BIC"], r["k_params"]))
    best, iq_model = results[0], None
    if (iq_bin := _resolve_iqtree2_bin()) and os.environ.get("PLANTMATICS_SKIP_IQTREE2") != "1" and n_eff >= 10:
        try:
            with tempfile.TemporaryDirectory() as td:
                fa = Path(td) / "aln.fasta"; fa.write_text("".join(f">{sid}\n{seqs[sid]}\n" for sid in present), encoding="utf-8")
                res = subprocess.run([iq_bin, "-s", str(fa), "-m", "MF", "-T", str(min(4, _resolve_max_workers())), "--prefix", str(Path(td) / "iq")], capture_output=True, text=True, timeout=15)
                m = re.search(r"Best-fit model according to BIC:\s*(\S+)", res.stdout); iq_model = m.group(1).strip() if m else None
        except Exception: iq_model = None
    rows = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    return {"alignment_name": alignment_name, "best_model": best["model_name"], "base_model": best["base_model"], "has_gamma": best["has_gamma"], "asc_applied": asc_applied, "bic_weight": best["bic_weight"], "ti_tv_kappa": best["ti_tv_kappa"], "kappa_purine": round(best["_kr"], 4), "kappa_pyrimidine": round(best["_ky"], 4), "gamma_alpha": best["gamma_alpha"], "prop_invariant": round(prop_inv, 4), "n_sites": n_sites, "n_var_sites": n_var, "pi_freqs": [round(float(x), 5) for x in pi_emp], "iqtree2_model": iq_model, "table_rows": rows}


def _apply_evo_distance(p1: np.ndarray, p2: np.ndarray, q: np.ndarray, spec: Optional[Dict[str, Any]]) -> np.ndarray:
    base = str((spec or {}).get("base_model") or "JC69").upper()
    has_g, alpha = bool((spec or {}).get("has_gamma", False)), float((spec or {}).get("gamma_alpha") or 1.0)
    pi = np.maximum(0.02, np.asarray((spec or {}).get("pi_freqs") or [0.25] * 4, dtype=np.float64))
    pi /= float(np.sum(pi))
    def _phi(u: np.ndarray) -> np.ndarray:
        uc = np.clip(u, 0.05, 1.0)
        return (alpha * (np.power(uc, -1.0 / max(0.05, alpha)) - 1.0)) if (has_g and alpha < 40.0) else -np.log(uc)
    if base == "K80":
        return 0.5 * _phi(1.0 - 2.0 * (p1 + p2) - q) + 0.25 * _phi(1.0 - 2.0 * q)
    if base == "F81":
        b = max(0.25, 1.0 - float(np.sum(pi * pi)))
        return b * _phi(1.0 - (p1 + p2 + q) / b)
    if base in ("HKY85", "TN93"):
        pi_r, pi_y = pi[0] + pi[2], pi[1] + pi[3]
        c1, c2 = 2.0 * pi[0] * pi[2] / pi_r, 2.0 * pi[1] * pi[3] / pi_y
        b = max(0.05, 2.0 * pi_r * pi_y - c1 * pi_y - c2 * pi_r)
        return c1 * _phi(1.0 - p1 / c1 - q / (2.0 * pi_r)) + c2 * _phi(1.0 - p2 / c2 - q / (2.0 * pi_y)) + b * _phi(1.0 - q / (2.0 * pi_r * pi_y))
    return 0.75 * _phi(1.0 - (4.0 / 3.0) * (p1 + p2 + q))


def compute_model_corrected_divergence_blocks(
    a1: np.ndarray,
    a2: np.ndarray,
    block_slices: Sequence[Tuple[int, int]],
    model_spec: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute BIC-model-corrected per-block pairwise & within-sample divergence counts."""
    n_samp, n_blk = a1.shape[0], len(block_slices)
    diff_blk, call_blk = np.zeros((n_blk, n_samp, n_samp)), np.zeros((n_blk, n_samp, n_samp))
    w_het, w_call = np.zeros((n_blk, n_samp)), np.zeros((n_blk, n_samp))
    for bi, (st, en) in enumerate(block_slices):
        x1, x2 = a1[:, st:en], a2[:, st:en]
        ok = x1 >= 0
        wc = np.sum(ok, axis=1).astype(np.float64)
        w_call[bi], xd = wc, np.where(ok, x1 ^ x2, 0)
        w_het[bi] = np.where(wc > 0, _apply_evo_distance(np.sum(ok & (xd == 2) & ((x1 & 1) == 0), axis=1) / np.maximum(1.0, wc), np.sum(ok & (xd == 2) & ((x1 & 1) == 1), axis=1) / np.maximum(1.0, wc), np.sum(ok & ((xd & 1) == 1), axis=1) / np.maximum(1.0, wc), model_spec) * wc, 0.0)
        hy = np.stack([x1, x2, x1, x2], axis=1)
        for i in range(n_samp):
            if not np.any(ok[i]):
                continue
            vij = ok & ok[i][None, :]
            cij = np.sum(vij, axis=1).astype(np.float64)
            call_blk[bi, i, :], den = cij, 4.0 * np.maximum(1.0, cij)
            hx = np.stack([x1[i], x1[i], x2[i], x2[i]], axis=0)[None, :, :]
            xor = np.where(vij[:, None, :], hx ^ hy, 0)
            p1 = np.sum(vij[:, None, :] & (xor == 2) & ((hx & 1) == 0), axis=(1, 2))
            p2 = np.sum(vij[:, None, :] & (xor == 2) & ((hx & 1) == 1), axis=(1, 2))
            qs = np.sum(vij[:, None, :] & ((xor & 1) == 1), axis=(1, 2))
            diff_blk[bi, i, :] = np.where(cij > 0, _apply_evo_distance(p1 / den, p2 / den, qs / den, model_spec) * cij, 0.0)
    return diff_blk, call_blk, w_het, w_call


def compute_coalescent_introgression_bf(
    d_blk_vec: np.ndarray, call_blk_vec: np.ndarray, bpp_out_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Compute introgression Bayes factor via BPP stepping-stone logML or nested coalescent mixture BIC."""
    if bpp_out_path is not None and bpp_out_path.exists() and (msc_p := bpp_out_path.with_name("bpp_msc_out.txt")).exists():
        pat = r"log\s*(?:ML|marginal likelihood)\s*[=:]\s*(-?\d+(?:\.\d+)?)"
        m1 = re.search(pat, bpp_out_path.read_text(encoding="utf-8", errors="replace"), re.I)
        m0 = re.search(pat, msc_p.read_text(encoding="utf-8", errors="replace"), re.I)
        if m1 and m0:
            l1, l0 = float(m1.group(1)), float(m0.group(1))
            lbf = max(-6.0, min(10.0, l1 - l0))
            return {"bf_method": "bpp_stepping_stone", "bf_10": round(math.exp(lbf), 4), "log_bf_10": round(lbf, 4), "bic_m0": round(-2.0 * l0, 4), "bic_m1": round(-2.0 * l1, 4), "w_introgr": None, "tau_introgr": None, "tau_sp": None}
    valid = np.isfinite(d_blk_vec) & np.isfinite(call_blk_vec) & (call_blk_vec >= 1.0) & (d_blk_vec >= 0.0)
    if int(np.sum(valid)) < 4:
        return {"bf_method": "bic_coalescent_mixture", "bf_10": 1.0, "log_bf_10": 0.0, "bic_m0": 0.0, "bic_m1": 0.0, "w_introgr": 0.0, "tau_introgr": 0.0, "tau_sp": 0.0}
    l_b, d_b = call_blk_vec[valid].astype(np.float64), d_blk_vec[valid].astype(np.float64)
    k_b = np.rint(np.maximum(0.0, d_b * l_b)).astype(np.int64)
    n_b, max_k = int(k_b.size), int(np.max(k_b))
    j_vec = np.arange(max_k + 1, dtype=np.float64)[None, :] if max_k <= 350 else None
    lgam_j = np.asarray([math.lgamma(j + 1.0) for j in range(max_k + 1)], dtype=np.float64)[None, :] if max_k <= 350 else None
    lgam = None if max_k <= 350 else np.asarray([math.lgamma(float(k) + 1.0) for k in k_b], dtype=np.float64)[:, None]

    def _comp_logp(tau: float, theta: float) -> np.ndarray:
        a, b = np.maximum(1e-12, l_b * max(0.0, float(tau)))[:, None], np.maximum(1e-12, l_b * max(1e-5, float(theta)))[:, None]
        p = b / (1.0 + b)
        if j_vec is not None and lgam_j is not None:
            t = np.where(j_vec <= k_b[:, None], -np.log1p(b) - a + j_vec * np.log(a) - lgam_j + (k_b[:, None] - j_vec) * np.log(p), -np.inf)
        else:
            t = _LAG_LOG_W[None, :] - np.log1p(b) - a + k_b[:, None] * np.log(np.maximum(1e-12, a + p * _LAG_X[None, :])) - lgam
        mx = np.max(t, axis=1)
        return mx + np.log(np.sum(np.exp(t - mx[:, None]), axis=1))

    mean_d, std_d = max(1e-5, float(np.mean(d_b))), max(1e-5, float(np.std(d_b)))
    (tau_0, theta_0), logl_0 = max(
        [((0.0, mean_d), float(np.sum(_comp_logp(0.0, mean_d))))] + [
            ((float(opt0[0]), float(opt0[1])), val0) for f in (0.0, 0.3, 0.6)
            for opt0, val0 in [_optimize_log_params(lambda p: float(np.sum(_comp_logp(p[0], p[1]))), [max(1e-6, f * mean_d), max(1e-4, (1.0 - f) * mean_d + 0.3 * std_d)], [(1e-6, 2.0), (1e-5, 1.0)])]
        ], key=lambda x: x[1],
    )
    q15, q80 = float(np.quantile(d_b, 0.15)), float(np.quantile(d_b, 0.80))
    seeds = [(0.5, tau_0, max(tau_0 * 1.02, tau_0 + 1e-6), theta_0)] + ([(w, max(1e-6, q15), max(1e-5, q80), max(2e-4, 0.35 * std_d)) for w in (0.20, 0.35, 0.50)] if q80 > q15 + 1e-4 else [])
    def _obj1(p: Sequence[float]) -> float:
        w, (t1, t2), th = min(0.98, max(0.02, p[0])), sorted((max(1e-6, p[1]), max(1e-6, p[2]))), max(1e-5, p[3])
        a, b = math.log(w) + _comp_logp(t1, th), math.log(1.0 - w) + _comp_logp(max(t2, t1 + 1e-6), th)
        m = np.maximum(a, b)
        return float(np.sum(m + np.log(np.exp(a - m) + np.exp(b - m))))
    best_m1 = max(
        [(logl_0, 0.5, tau_0, max(tau_0 * 1.02, tau_0 + 1e-6), theta_0)] + [
            (val, min(0.98, max(0.02, opt[0])), *sorted((opt[1], opt[2])), opt[3])
            for seed1 in seeds for opt, val in [_optimize_log_params(_obj1, list(seed1), [(0.02, 0.98), (1e-6, 2.0), (1e-6, 2.0), (1e-5, 1.0)])]
        ], key=lambda x: x[0],
    )
    bic_0, bic_1 = -2.0 * logl_0 + 2.0 * math.log(n_b), -2.0 * max(logl_0, best_m1[0]) + 4.0 * math.log(n_b)
    lbf = max(-6.0, min(10.0, 0.5 * (bic_0 - bic_1)))
    return {"bf_method": "bic_coalescent_mixture", "bf_10": round(math.exp(lbf), 4), "log_bf_10": round(lbf, 4), "bic_m0": round(bic_0, 4), "bic_m1": round(bic_1, 4), "w_introgr": round(best_m1[1], 4), "tau_introgr": round(best_m1[2], 6), "tau_sp": round(best_m1[3], 6)}
