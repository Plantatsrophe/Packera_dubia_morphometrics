#!/usr/bin/env python3
"""Stage public Packera cytological, C-value, and label chromosome records flat into `data/cytology/`."""
from __future__ import annotations

import argparse, csv, hashlib, html, io, json, logging, os, re, sqlite3, sys, tempfile, time, zipfile
import urllib.error, urllib.parse, urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.db import get_db_connection, resolve_db_path, upsert_cytological_records, upsert_taxon_cytology_summary
from scripts.core.harvester import _DYNAMIC_SYNONYM_TO_CANONICAL, standardize_packera_taxon
from scripts.core.io_utils import atomic_write_bytes, atomic_write_csv_rows, ensure_flat_filename
from scripts.core.taxonomy_authority import generate_exsiccatae_key, normalize_collector, normalize_record_number

logger = logging.getLogger("plantmatics.download_public_cytology")
_UA = "Plantmatics-UNC-NCU/2026.11 (+https://herbarium.unc.edu/)"
RECORD_COLUMNS: Tuple[str, ...] = ("cytology_id", "source_db", "source_record_id", "original_taxon_name", "basionym_name", "species_standardized", "catalogNumber", "exsiccatae_key", "cytology_linkage_tier", "raw_count_str", "gametophytic_n", "sporophytic_2n", "base_x", "ploidy_level", "cytotype_label", "b_chromosomes", "is_aneuploid_or_odd_ploid", "is_approximate", "c_value_1c_pg", "c_value_2c_pg", "breeding_system", "genbank_accession", "collector_surname", "record_number_clean", "institution_code", "voucher_citation", "locality_text", "reference_citation")
SUMMARY_COLUMNS: Tuple[str, ...] = ("species_standardized", "n_cytology_records", "n_direct_cytovouchers", "modal_sporophytic_2n", "modal_base_x", "modal_ploidy_level", "reported_cytotypes", "is_multi_cytotype_taxon", "has_dysploid_series", "has_odd_polyploid_or_apomixis", "mean_2c_pg", "breeding_system_summary")
_X22, _X23, _OLD_WORLD_2N = {22 * k: k for k in range(2, 9)}, {23 * k: k for k in range(2, 9)}, {10, 20, 30, 40, 60, 80}
_DB_SYN_CACHE: Dict[str, Dict[str, str]] = {}
_B_CHR_RE = re.compile(r"(?:(?:\+|,|\s)\s*(?:(\d+)\s*[-–]\s*)?(\d+)|(\d+))\s*[Bb]\b")
_APPROX_RE = re.compile(r"\b(?:ca\.?|c\.?|approx\.?)\b|[~?]", re.I)
_MEIOTIC_RE, _MULT_MAP = re.compile(r"(\d+)\s*(IV|III|II|I)s?\b"), {"IV": 4, "III": 3, "II": 2, "I": 1}
_SPLIT_RE = re.compile(r"\s*[;]\s*|\s+\bor\b\s+|\s*,\s*(?![+\s]*(?:\d+\s*[-–]\s*)?\d+\s*[Bb]\b)", re.I)


def _default_db_synonyms() -> Dict[str, str]:
    dbp = resolve_db_path(workspace_root=PROJECT_ROOT)
    ck = f"{os.environ.get('PLANTMATICS_DB_PATH', '')}|{os.environ.get('PACKERA_DB_PATH', '')}|{dbp}"
    if ck not in _DB_SYN_CACHE:
        _DB_SYN_CACHE[ck] = {}
        if Path(dbp).exists():
            try:
                with sqlite3.connect(str(dbp)) as conn:
                    _DB_SYN_CACHE[ck] = {str(k).strip().lower(): str(v).strip() for k, v in conn.execute("SELECT raw_name_lower, canonical_packera FROM taxon_synonyms WHERE canonical_packera LIKE 'Packera %'") if k and v}
            except Exception: pass
    return _DB_SYN_CACHE[ck]
_NON_CNT = r"(?!\s*(?:mm|cm|dm|m|km|in|ft|plants?|specimens?|sheets?|samples?|individuals?|vouchers?|stems?|leaves|heads?|flowers?|collections?|sites?|populations?|plots?|reps?|years?|branched?)\b)"
_RE_2N = re.compile(r"\b2n\s*=\s*(?:(?:ca\.?|c\.?|approx\.?|~)\s*)?\d+\b(?:\s*(?:[-–]\s*\d+\b|\+\s*(?:\d+\s*[-–]\s*)?\d+\s*[Bb]\b|\s*[,;]\s*\d+\b" + _NON_CNT + r"))*" + _NON_CNT)
_RE_N = re.compile(r"\bn\s*=\s*(?:(?:ca\.?|c\.?|approx\.?|~)\s*)?\d+(?:\s*(?:IV|III|II|I)\b(?:\s*\+\s*\d+\s*(?:IV|III|II|I)\b)*|\b(?:\s*(?:\+\s*(?:\d+\s*[-–]\s*)?\d+\s*[Bb]\b|\s*[,;]\s*\d+\b" + _NON_CNT + r"))*)" + _NON_CNT)
_RE_CHR = re.compile(r"\bchromosome\s+(?:count|number|voucher)\b[^.;\n]{0,80}?\b(?:2n\s*=\s*|n\s*=\s*)?(?:ca\.?\s*)?\d+(?:\s*(?:IV|III|II|I)\b|\b)" + _NON_CNT, re.I)
_RE_CVAL = re.compile(r"\b(?:1C|2C)\s*(?:=|:)?\s*\d+(?:\.\d+)?\s*pg\b", re.I)


class _CrossHostAuthStripRedirectHandler(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req: urllib.request.Request, fp: Any, code: int, msg: str, headers: Any, newurl: str):
        nr = super().redirect_request(req, fp, code, msg, headers, newurl)
        if nr and urllib.parse.urlparse(req.full_url).netloc.lower() != urllib.parse.urlparse(newurl).netloc.lower():
            for h in ("Authorization", "Cookie"):
                nr.headers.pop(h, None); nr.unredirected_hdrs.pop(h, None)
        return nr


def _ensure_flat_filename(root: Path, filename: str) -> Path:
    root_res = root.resolve()
    dest = (root_res / ensure_flat_filename(filename)).resolve()
    if dest.parent != root_res:
        raise ValueError(f"Path escapes {root_res}: {dest}")
    return dest


def _http_request(url: str, timeout: int = 25, retries: int = 3) -> bytes:
    opener, last_err = urllib.request.build_opener(_CrossHostAuthStripRedirectHandler), None
    for attempt in range(retries):
        try:
            with opener.open(urllib.request.Request(url, headers={"User-Agent": _UA, "Accept": "application/json,text/html,*/*"}), timeout=timeout) as r: return r.read()
        except Exception as exc:
            last_err = exc
            if (isinstance(exc, urllib.error.HTTPError) and exc.code not in (429, 500, 502, 503, 504)) or attempt == retries - 1: raise
        time.sleep(1.2 * (2 ** attempt))
    raise last_err if last_err else RuntimeError(url)


def _classify_cytotype(n_val: Optional[int], sp_2n: int, b_chr: int, approx: int, raw_str: str, force_odd: bool = False) -> Dict[str, Any]:
    bx, pl = (22, _X22[sp_2n]) if sp_2n in _X22 else ((23, _X23[sp_2n]) if sp_2n in _X23 else (None, None))
    is_odd = int(force_odd or pl is None or pl in (3, 5, 7) or sp_2n % 2 == 1)
    n_eff = None if (force_odd or sp_2n % 2 == 1) else (sp_2n // 2 if (pl is not None and n_val is None) else n_val)
    return {"raw_count_str": raw_str, "gametophytic_n": n_eff, "sporophytic_2n": sp_2n, "base_x": bx, "ploidy_level": pl, "cytotype_label": f"{pl}x" if pl is not None else "aneuploid", "b_chromosomes": b_chr, "is_aneuploid_or_odd_ploid": is_odd, "is_approximate": approx}


def _parse_single_text(text: str, default_mode: Optional[str] = None) -> List[Dict[str, Any]]:
    raw_orig = " ".join(str(text or "").split())
    if not raw_orig: return []
    cleaned = raw_orig
    for _ in range(4): cleaned = re.sub(r"\([^()]*\)|\[[^\[\]]*\]", " ", cleaned)
    if not cleaned.strip(" -?"): return []
    g_approx, out, inherited = int(bool(_APPROX_RE.search(cleaned))), [], default_mode
    for piece in (p.strip() for p in _SPLIT_RE.split(cleaned) if p.strip()):
        approx, bm = int(bool(g_approx or _APPROX_RE.search(piece))), _B_CHR_RE.search(piece)
        b_chr = int(bm.group(2) or bm.group(3) or 0) if bm else 0
        piece = _APPROX_RE.sub("", _B_CHR_RE.sub("", piece) if bm else piece).strip()
        mode_cur = inherited
        if re.match(r"^\s*2n\b", piece, re.I): mode_cur = inherited = "2n"; piece = re.sub(r"^\s*2n\s*=?\s*", "", piece, flags=re.I).strip()
        elif re.match(r"^\s*n\b", piece, re.I): mode_cur = inherited = "n"; piece = re.sub(r"^\s*n\s*=?\s*", "", piece, flags=re.I).strip()
        m_terms = _MEIOTIC_RE.findall(piece)
        if m_terms:
            sp_2n, has_un = sum(int(c) * _MULT_MAP[s] for c, s in m_terms), any(s in ("I", "III") for _, s in m_terms)
            if 8 <= sp_2n <= 250: out.append(_classify_cytotype(None if has_un or sp_2n % 2 else sp_2n // 2, sp_2n, b_chr, approx, raw_orig, has_un))
            continue
        rng = re.match(r"^\s*(\d+)\s*[-–]\s*(\d+)\s*$", piece)
        for val in ([int(rng.group(1)), int(rng.group(2))] if rng else [int(x) for x in re.findall(r"\b\d+\b", piece)]):
            m = mode_cur or ("n" if val < 36 else "2n")
            if m == "n" and 4 <= val <= 125: out.append(_classify_cytotype(val, val * 2, b_chr, approx or int(bool(rng)), raw_orig))
            elif m == "2n" and 8 <= val <= 250: out.append(_classify_cytotype(None if val % 2 else val // 2, val, b_chr, approx or int(bool(rng)), raw_orig))
    return out


def parse_packera_chromosome_record(raw_count: Any, gametophytic_raw: Any = None, sporophytic_raw: Any = None) -> List[Dict[str, Any]]:
    raw_s, g_s, s_s = (" ".join(str(x or "").split()) for x in (raw_count, gametophytic_raw, sporophytic_raw))
    if raw_s:
        parsed = _parse_single_text(raw_s)
        if parsed: return parsed
    g_p, s_p = (_parse_single_text(g_s, "n") if g_s else []), (_parse_single_text(s_s, "2n") if s_s else [])
    raw_fb = raw_s or "; ".join(p for p, v in ((f"n={g_s}", g_p), (f"2n={s_s}", s_p)) if v)
    if len(g_p) == 1 and len(s_p) == 1 and g_p[0]["sporophytic_2n"] == s_p[0]["sporophytic_2n"]:
        return [dict(s_p[0], gametophytic_n=g_p[0]["gametophytic_n"], raw_count_str=raw_fb)]
    seen, uniq = set(), []
    for rec in s_p + g_p:
        k = (rec["sporophytic_2n"], rec["b_chromosomes"])
        if k not in seen: seen.add(k); rec["raw_count_str"] = raw_fb; uniq.append(rec)
    return uniq


def _strip_authors(name: str) -> str:
    s = re.split(r"\s+\bas\b\s+", re.sub(r"\s*\[[^\]]*\]", " ", str(name or "")), maxsplit=1, flags=re.I)[0]
    for _ in range(4):
        s = re.sub(r"\([^()]*\)", " ", s)
    toks = [t for t in s.replace("×", " ").split() if t]
    out, i = toks[:2], 2
    while i < len(toks) - 1 and toks[i].lower() in ("var.", "subsp.", "ssp.", "f."):
        out.extend([toks[i].lower(), toks[i + 1]]); i += 2
    return " ".join(out)


def resolve_packera_cytology_taxon(
    raw_name: Any, synonym_map: Optional[Mapping[str, str]] = None,
    sporophytic_2n: Optional[int] = None, allow_standardize_fallback: bool = True,
) -> Optional[Dict[str, Optional[str]]]:
    orig = " ".join(str(raw_name or "").split())
    clean = _strip_authors(orig)
    words = clean.split()
    if len(words) < 2:
        return None
    smap = dict(_default_db_synonyms())
    for src in (_DYNAMIC_SYNONYM_TO_CANONICAL, synonym_map or {}):
        smap.update({str(k).strip().lower(): v for k, v in src.items() if isinstance(v, str) and v.startswith("Packera ")})
    binom, is_sen = f"{words[0].lower()} {words[1].lower()}", words[0].lower() in ("senecio", "cineraria")
    if re.search(r"\bledeb\.?\b", orig, re.I) or ((is_sen or "tomentos" in orig.lower()) and sporophytic_2n is not None and int(sporophytic_2n) in _OLD_WORLD_2N):
        return None
    if clean.lower().startswith("senecio aureus var. pauperculus"):
        return {"original_taxon_name": orig, "basionym_name": clean, "species_standardized": "Packera paupercula"}
    for cand in (orig.lower(), clean.lower(), binom):
        if cand in smap:
            return {"original_taxon_name": orig, "basionym_name": clean if is_sen else None, "species_standardized": smap[cand]}
    if words[0].lower() == "packera":
        std = standardize_packera_taxon(clean) or f"Packera {words[1].lower()}"
        return {"original_taxon_name": orig, "basionym_name": None, "species_standardized": std} if std.startswith("Packera ") else None
    if is_sen and allow_standardize_fallback and not (sporophytic_2n is not None and int(sporophytic_2n) in _OLD_WORLD_2N):
        std = standardize_packera_taxon(clean)
        if std and std.startswith("Packera "):
            return {"original_taxon_name": orig, "basionym_name": clean, "species_standardized": std}
    return None


def extract_label_cytology_mentions(text: Any) -> List[str]:
    s = " ".join(str(text or "").split())
    spans: List[Tuple[int, int, str]] = []
    for rx in (_RE_2N, _RE_N, _RE_CHR, _RE_CVAL):
        for m in rx.finditer(s):
            a, b = m.span()
            if not any(not (b <= ea or a >= eb) for ea, eb, _ in spans):
                spans.append((a, b, m.group(0).strip()))
    return list(dict.fromkeys(t for _, _, t in sorted(spans, key=lambda x: x[0]))) if s else []


def _load_sqlite_context(db_path: Optional[Path]) -> Tuple[Dict[str, str], Dict[str, Any]]:
    smap, ctx = {}, {"by_cat": {}, "by_exsic": defaultdict(list), "by_cn": defaultdict(list), "gb_to_cat": {}, "vouchers": []}
    if not db_path or not db_path.exists():
        return smap, ctx
    conn = sqlite3.connect(str(db_path)); conn.row_factory = sqlite3.Row
    try:
        tbls = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "taxon_synonyms" in tbls:
            smap.update({str(r[0]).strip().lower(): str(r[1]).strip() for r in conn.execute("SELECT raw_name_lower, canonical_packera FROM taxon_synonyms WHERE canonical_packera LIKE 'Packera %'")})
        if "vouchers" in tbls:
            vcols = {r[1] for r in conn.execute("PRAGMA table_info(vouchers)")}
            sel = [c for c in ("catalogNumber", "species_standardized", "scientificName", "recordedBy", "recordNumber", "institutionCode", "exsiccatae_key", "eventDate", "year", "month", "day", "locality", "occurrenceRemarks", "dynamicProperties", "identificationRemarks") if c in vcols]
            for r in conn.execute(f"SELECT {', '.join(sel)} FROM vouchers"):
                d = {k: r[k] for k in sel}
                cat, sp, sci = str(d.get("catalogNumber") or "").strip(), str(d.get("species_standardized") or "").strip(), str(d.get("scientificName") or "").strip()
                if sci and sp.startswith("Packera "): smap.setdefault(sci.lower(), sp)
                if not cat: continue
                ex = str(d.get("exsiccatae_key") or "").strip() or generate_exsiccatae_key(d.get("recordedBy"), d.get("recordNumber"), d.get("eventDate"), d.get("year"), d.get("month"), d.get("day"))
                d["_ex"] = ex; ctx["vouchers"].append(d); ctx["by_cat"][cat.lower()] = d
                if ex: ctx["by_exsic"][ex].append(d)
                cs, rn = normalize_collector(d.get("recordedBy")), normalize_record_number(d.get("recordNumber"))
                if cs and rn and cs != "unknown": ctx["by_cn"][(cs, rn)].append(d)
        if "genomic_samples" in tbls:
            gcols = {r[1] for r in conn.execute("PRAGMA table_info(genomic_samples)")}
            gsel = [c for c in ("genetic_sample_id", "catalogNumber", "genbank_accession") if c in gcols]
            if "catalogNumber" in gsel:
                for r in conn.execute(f"SELECT {', '.join(gsel)} FROM genomic_samples"):
                    cat = str(r["catalogNumber"] or "").strip()
                    for k in ("genbank_accession", "genetic_sample_id"):
                        if cat and k in gsel and r[k]: ctx["gb_to_cat"][str(r[k]).strip().upper()] = cat
    finally:
        conn.close()
    return smap, ctx


def _opt_num(val: Any, cast: Any = float) -> Any:
    s = str(val or "").strip()
    return cast(s) if re.match(r"^-?\d+(\.\d+)?$" if cast is float else r"^-?\d+$", s) else None


def _make_record(base: Dict[str, Any], parsed: Dict[str, Any], tax: Dict[str, Optional[str]], idx: int = 0) -> Dict[str, Any]:
    row = {**dict.fromkeys(RECORD_COLUMNS, None), **base, **parsed, **tax}
    sig = "|".join(str(row.get(k) or "") for k in ("source_db", "source_record_id", "species_standardized", "sporophytic_2n", "catalogNumber", "voucher_citation", "reference_citation", "c_value_2c_pg")) + f"|{idx}"
    row["cytology_id"] = "CYTO_" + hashlib.sha1(sig.encode("utf-8")).hexdigest()[:16].upper()
    return row


def _row_to_records(low: Dict[str, str], smap: Dict[str, str], default_src: str, rec_id: str, allow_fb: bool) -> List[Dict[str, Any]]:
    pick = lambda *ks: next((low[k] for k in ks if low.get(k) and low[k] != "-"), "")
    sp_raw = pick("original_taxon_name", "species_standardized", "taxon_name", "taxon", "species", "scientific_name", "scientificname", "accepted_name", "binomial")
    if not sp_raw and low.get("genus"):
        gw, sub = re.sub(r"[^A-Za-z]", "", low["genus"].split()[0]).capitalize(), low.get("species / subspecies", "")
        if re.search(r"\s+\bas\b\s+", sub, re.I):
            acc_p, bas_p = re.split(r"\s+\bas\b\s+", sub, maxsplit=1, flags=re.I)
            bas_p = re.sub(r"^s\.\s+", "Senecio ", bas_p.strip(), flags=re.I)
            acc_c, bas_c = _strip_authors(f"{gw} {acc_p}"), _strip_authors(bas_p)
            std_a = standardize_packera_taxon(acc_c) or acc_c if gw == "Packera" else ""
            is_bas = bas_c.lower().startswith(("senecio ", "cineraria "))
            if std_a.startswith("Packera ") and is_bas: smap.setdefault(bas_c.lower(), std_a)
            sp_raw = bas_p if is_bas else f"{gw} {acc_p}".strip()
        else: sp_raw = f"{gw} {sub}".strip()
    var_v = pick("variety")
    vcit = pick("voucher_citation", "voucher") or (re.sub(r"^based on\s+", "", var_v, flags=re.I).strip() if var_v.lower().startswith("based on ") else "")
    if var_v and re.match(r"^(?:var\.|subsp\.|ssp\.|f\.)\b", var_v, re.I) and sp_raw: sp_raw = f"{sp_raw} {var_v}"
    gb_m = re.search(r"\b([A-Z]{2}\d{6})\b", f"{var_v} {low.get('species / subspecies', '')}")
    ref = pick("reference_citation", "reference", "citation", "author name & published year") or " ".join(x for x in (low.get(k, "") for k in ("author", "year", "title", "journal", "vol.(number)", "page")) if x and x != "-")
    counts = parse_packera_chromosome_record(pick("raw_count_str", "raw_count", "count", "chromosome_count", "chromosome_number"), pick("gametophytic_n", "n", "chromosome number n=(no. of plants examined)"), pick("sporophytic_2n", "2n", "chromosome number 2n=(total karyotypic length)(no. of plants examined)"))
    out: List[Dict[str, Any]] = []
    for j, c in enumerate(counts):
        tax = resolve_packera_cytology_taxon(sp_raw, smap, c.get("sporophytic_2n"), allow_standardize_fallback=allow_fb)
        if tax: out.append(_make_record({"source_db": pick("source_db") or default_src, "source_record_id": pick("source_record_id", "id") or f"{rec_id}_{j}", "catalogNumber": pick("catalognumber", "catalog_number") or None, "exsiccatae_key": pick("exsiccatae_key") or None, "c_value_1c_pg": _opt_num(pick("c_value_1c_pg", "1c")), "c_value_2c_pg": _opt_num(pick("c_value_2c_pg", "2c")), "breeding_system": pick("breeding_system", "breeding system") or None, "genbank_accession": pick("genbank_accession", "genbank") or (gb_m.group(1) if gb_m else None), "collector_surname": pick("collector_surname", "collector", "recordedby") or None, "record_number_clean": pick("record_number_clean", "record_number", "recordnumber") or None, "institution_code": pick("institution_code", "institutioncode") or None, "voucher_citation": vcit or None, "locality_text": pick("locality_text", "locality", "collected locality") or None, "reference_citation": ref or None}, c, tax, j))
    return out


def _fetch_ccdb(smap: Dict[str, str]) -> List[Dict[str, Any]]:
    raw_items: List[Dict[str, Any]] = []
    for u in ("http://ccdb.tau.ac.il/services/countsFull/?genus=Packera&format=json", "https://taux.evolseq.net/CCDB_web/services/countsFull/?genus=Packera&format=json"):
        try:
            p = json.loads(_http_request(u, 20, 2).decode("utf-8", errors="replace"))
            raw_items = p if isinstance(p, list) else (p.get("results") or p.get("data") or [])
            if raw_items: break
        except Exception as exc: logger.warning("CCDB endpoint %s unavailable: %s", u, exc)
    out: List[Dict[str, Any]] = []
    for r in (x for x in raw_items if isinstance(x, dict)):
        orig = str(r.get("original_name") or r.get("matched_name") or r.get("resolved_name") or "").strip()
        res_n = str(r.get("resolved_name") or r.get("matched_name") or "").strip()
        if res_n.startswith("Packera ") and orig.lower().startswith(("senecio ", "cineraria ")):
            smap[_strip_authors(orig).lower()] = _strip_authors(res_n)
        cts = parse_packera_chromosome_record("", r.get("gametophytic"), r.get("sporophytic")) or (parse_packera_chromosome_record("", r.get("parsed_n"), "") if r.get("parsed_n") else [])
        for idx, c in enumerate(cts):
            tax = resolve_packera_cytology_taxon(orig or res_n, smap, c["sporophytic_2n"])
            if tax: out.append(_make_record({"source_db": "CCDB", "source_record_id": str(r.get("id") or r.get("source") or ""), "voucher_citation": str(r.get("voucher") or "").strip() or None, "reference_citation": str(r.get("reference") or "").strip() or None}, c, tax, idx))
    return out


def _fetch_icna(smap: Dict[str, str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    try: zbytes = _http_request("https://dbarchive.biosciencedbc.jp/data/asteraceae-chr-index/LATEST/index_chromosome_number_asteraceae.zip", 30, 2)
    except Exception as exc: logger.warning("ICNA archive unavailable: %s", exc); return out
    with tempfile.TemporaryDirectory() as tmp:
        tmp_res = Path(tmp).resolve()
        with zipfile.ZipFile(io.BytesIO(zbytes)) as zf:
            for member in zf.namelist():
                t = (tmp_res / member).resolve()
                if tmp_res not in t.parents and t != tmp_res: raise ValueError(f"Unsafe zip member: {member}")
            zf.extractall(tmp_res)
        for fp in (p for p in tmp_res.rglob("*") if p.is_file() and p.suffix.lower() in (".csv", ".tsv", ".txt")):
            with open(fp, "r", encoding="utf-8", errors="replace", newline="") as fh:
                for i, r in enumerate(csv.DictReader(fh, delimiter="\t" if fp.suffix.lower() in (".tsv", ".txt") else ",")):
                    low = {str(k or "").strip().lower().replace("\uff08", "(").replace("\uff06", "&"): (v or "").strip() for k, v in r.items()}
                    if "genus" in low and re.sub(r"[^a-z]", "", (low["genus"].split()[:1] or [""])[0].lower()) not in ("packera", "senecio", "cineraria"): continue
                    out.extend(_row_to_records(low, smap, "ICNA", f"icna_{i}", False))
    return out


def _fetch_fna(smap: Dict[str, str]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    try: idx_html = _http_request("http://www.efloras.org/florataxon.aspx?flora_id=1&taxon_id=123645", 15, 1).decode("utf-8", errors="replace")
    except Exception as exc: logger.warning("FNA fetch unavailable: %s", exc); return out
    for tid, name in re.findall(r'florataxon\.aspx\?flora_id=1&taxon_id=(\d+)[^>]*>([^<]*Packera[^<]+)<', idx_html, re.I)[:40]:
        try:
            plain = html.unescape(re.sub(r"<[^>]+>", " ", _http_request(f"http://www.efloras.org/florataxon.aspx?flora_id=1&taxon_id={tid}", 12, 1).decode("utf-8", errors="replace")))
            for k, m in enumerate(extract_label_cytology_mentions(plain)):
                for j, c in enumerate(parse_packera_chromosome_record(m)):
                    tax = resolve_packera_cytology_taxon(html.unescape(name), smap, c["sporophytic_2n"])
                    if tax: out.append(_make_record({"source_db": "FNA", "source_record_id": f"fna_{tid}_{k}", "reference_citation": "Flora of North America, Vol. 20 (efloras.org)"}, c, tax, j))
        except Exception: continue
    return out


def _scan_labels(smap: Dict[str, str], ctx: Mapping[str, Any]) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for v in ctx["vouchers"]:
        cat = str(v.get("catalogNumber") or "").strip()
        txt = " ; ".join(str(v.get(c) or "") for c in ("occurrenceRemarks", "dynamicProperties", "identificationRemarks", "locality") if v.get(c))
        for k, mention in enumerate(extract_label_cytology_mentions(txt)):
            c1 = c2 = None; m_c = re.match(r"^(1C|2C)\s*(?:=|:)?\s*(\d+(?:\.\d+)?)\s*pg$", mention, re.I)
            if m_c:
                val = float(m_c.group(2)); c1, c2 = (val, val * 2.0) if m_c.group(1).upper() == "1C" else (val / 2.0, val)
                cnts = [{"raw_count_str": mention, "gametophytic_n": None, "sporophytic_2n": None, "base_x": None, "ploidy_level": None, "cytotype_label": None, "b_chromosomes": 0, "is_aneuploid_or_odd_ploid": 0, "is_approximate": 0}]
            else: cnts = parse_packera_chromosome_record(mention)
            for j, c in enumerate(cnts):
                tax = resolve_packera_cytology_taxon(v.get("species_standardized") or v.get("scientificName") or "", smap, c.get("sporophytic_2n"))
                if tax: out.append(_make_record({"source_db": "Herbarium_Label", "source_record_id": f"{cat}_{k}", "catalogNumber": cat, "exsiccatae_key": v.get("_ex") or None, "c_value_1c_pg": c1, "c_value_2c_pg": c2, "collector_surname": normalize_collector(v.get("recordedBy")), "record_number_clean": normalize_record_number(v.get("recordNumber")), "institution_code": str(v.get("institutionCode") or "").strip() or None, "locality_text": str(v.get("locality") or "").strip() or None}, c, tax, j))
    return out


def _read_local_csv(csv_path: Path, smap: Dict[str, str], *, reload_cache: bool = False) -> List[Dict[str, Any]]:
    out, seen = [], set()
    with open(csv_path, "r", encoding="utf-8", newline="") as fh:
        rows_raw = [{str(k or "").strip().lower(): (str(v).strip() if v is not None else "") for k, v in r.items()} for r in csv.DictReader(fh, delimiter="\t" if csv_path.suffix.lower() in (".tsv", ".txt") else ",")]
    rk_of = lambda lw: tuple(lw.get(x) for x in ("source_db", "source_record_id", "species_standardized", "sporophytic_2n", "raw_count_str", "exsiccatae_key"))
    if reload_cache: seen.update(rk_of(lw) for lw in rows_raw if lw.get("cytology_linkage_tier") == "C1_Direct_Cytovoucher")
    for i, low in enumerate(rows_raw):
        if reload_cache:
            if low.get("source_db") == "Herbarium_Label" or (low.get("cytology_linkage_tier") == "C2_Exsiccatae_Cytoduplicate" and (rk_of(low) in seen or seen.add(rk_of(low)))): continue
            tax = resolve_packera_cytology_taxon(low.get("original_taxon_name") or low.get("species_standardized"), smap, _opt_num(low.get("sporophytic_2n"), int))
            if tax:
                rec = {c: (low.get(c.lower()) or None) for c in RECORD_COLUMNS}
                for ic in ("gametophytic_n", "sporophytic_2n", "base_x", "ploidy_level", "b_chromosomes", "is_aneuploid_or_odd_ploid", "is_approximate"): rec[ic] = _opt_num(low.get(ic), int)
                for fc in ("c_value_1c_pg", "c_value_2c_pg"): rec[fc] = _opt_num(low.get(fc), float)
                if rec["cytology_linkage_tier"] != "C1_Direct_Cytovoucher": rec["catalogNumber"] = None
                rec.update(tax); out.append(rec)
            continue
        out.extend(_row_to_records(low, smap, "Local_CSV", f"local_{i}", True))
    return out


def _cross_link_and_summarize(records: List[Dict[str, Any]], ctx: Mapping[str, Any]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    linked, seen_ids = [], set()
    for rec in records:
        r = dict(rec); base_id = re.sub(r"_[0-9A-F]{6}$", "", str(r.get("cytology_id") or ""))
        cat = str(r.get("catalogNumber") or "").strip() or ctx["gb_to_cat"].get(str(r.get("genbank_accession") or "").strip().upper(), "")
        vc = str(r.get("voucher_citation") or "")
        if not r.get("collector_surname") and vc:
            vm = re.search(r"\b([A-Z][a-zA-Z\-]+)\s*,\s*(?:[A-Z]\.\s*){1,3}(?:et\s+al\.?\s*)?#?\s*(\d+[A-Za-z]?)\b", vc) or re.search(r"\b([A-Z][a-zA-Z\-]+)\s+#?\s*(\d+[A-Za-z]?)\b", vc)
            if vm: r["collector_surname"], r["record_number_clean"] = vm.group(1), vm.group(2)
        if not r.get("institution_code") and vc:
            im = re.search(r"\(([A-Z]{2,8})\)", vc)
            if im: r["institution_code"] = im.group(1)
        cs, rn = normalize_collector(r.get("collector_surname")), normalize_record_number(r.get("record_number_clean"))
        ex = str(r.get("exsiccatae_key") or "").strip() or (generate_exsiccatae_key(cs, rn) if (cs and cs != "unknown" and rn) else "")
        r.update({"collector_surname": cs if cs != "unknown" else r.get("collector_surname"), "record_number_clean": rn or r.get("record_number_clean"), "exsiccatae_key": ex or None})
        emit: List[Dict[str, Any]] = []
        if cat:
            v = ctx["by_cat"].get(cat.lower())
            r.update({"catalogNumber": v["catalogNumber"] if v else cat, "exsiccatae_key": (v.get("_ex") if v else None) or ex or None, "cytology_linkage_tier": "C1_Direct_Cytovoucher", "_bio": True, "cytology_id": base_id})
            emit.append(r)
            for sib in (ctx["by_exsic"].get(r["exsiccatae_key"], []) if r["exsiccatae_key"] else []):
                if sib["catalogNumber"] != r["catalogNumber"]:
                    emit.append(dict(r, catalogNumber=sib["catalogNumber"], institution_code=sib.get("institutionCode") or r.get("institution_code"), cytology_linkage_tier="C2_Exsiccatae_Cytoduplicate", _bio=False, cytology_id=base_id + "_" + hashlib.sha1(sib["catalogNumber"].encode()).hexdigest()[:6].upper()))
        else:
            matches = (ctx["by_exsic"].get(ex, []) if ex else []) or (ctx["by_cn"].get((cs, rn), []) if (cs != "unknown" and rn) else [])
            inst = str(r.get("institution_code") or "").strip().upper()
            emit = [dict(r, catalogNumber=mv["catalogNumber"], exsiccatae_key=mv.get("_ex") or ex or None, cytology_linkage_tier=("C1_Direct_Cytovoucher" if (inst and str(mv.get("institutionCode") or mv.get("institution_code") or "").strip().upper() == inst) else "C2_Exsiccatae_Cytoduplicate"), _bio=(m_idx == 0), cytology_id=base_id + "_" + hashlib.sha1(mv["catalogNumber"].encode()).hexdigest()[:6].upper()) for m_idx, mv in enumerate(matches)] if matches else [dict(r, cytology_linkage_tier="UNLINKED", _bio=True, cytology_id=base_id)]
        for er in emit:
            if er["cytology_id"] not in seen_ids: seen_ids.add(er["cytology_id"]); linked.append(er)
    by_sp: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in linked: by_sp[r["species_standardized"]].append(r)
    summaries: List[Dict[str, Any]] = []
    for sp, grp in sorted(by_sp.items()):
        bio = [g for g in grp if g.get("_bio", True)]
        sp2n = [int(g["sporophytic_2n"]) for g in bio if g.get("sporophytic_2n") is not None]
        modal_2n = min(Counter(sp2n).items(), key=lambda kv: (-kv[1], kv[0]))[0] if sp2n else None
        pl_set, bx_set = sorted({int(g["ploidy_level"]) for g in bio if g.get("ploidy_level") is not None}), sorted({int(g["base_x"]) for g in bio if g.get("base_x") is not None})
        labels = [f"{p}x" for p in pl_set] + (["aneuploid"] if any(g.get("cytotype_label") == "aneuploid" for g in bio) else [])
        is_multi = int(len(pl_set) >= 2 or len(bx_set) >= 2 or len(set(sp2n)) >= 2)
        c2_vals, bs_vals = [float(g["c_value_2c_pg"]) for g in bio if g.get("c_value_2c_pg") is not None], sorted({str(g["breeding_system"]).strip() for g in bio if g.get("breeding_system")})
        for g in grp:
            if g.get("cytology_linkage_tier") == "UNLINKED": g["cytology_linkage_tier"] = "C4_Taxon_MultiCytotype" if is_multi else "C3_Taxon_Monomorphic"
        summaries.append({
            "species_standardized": sp, "n_cytology_records": len(bio), "n_direct_cytovouchers": len({g["catalogNumber"] for g in grp if g.get("cytology_linkage_tier") == "C1_Direct_Cytovoucher" and g.get("catalogNumber")}),
            "modal_sporophytic_2n": modal_2n, "modal_base_x": (22 if modal_2n in _X22 else (23 if modal_2n in _X23 else None)) if modal_2n else None,
            "modal_ploidy_level": (_X22.get(modal_2n) or _X23.get(modal_2n)) if modal_2n else None, "reported_cytotypes": ", ".join(labels) if labels else None,
            "is_multi_cytotype_taxon": is_multi, "has_dysploid_series": int(len(bx_set) >= 2 or any(g.get("cytotype_label") == "aneuploid" for g in bio)),
            "has_odd_polyploid_or_apomixis": int(any(g.get("is_aneuploid_or_odd_ploid") == 1 or "apom" in str(g.get("breeding_system") or "").lower() for g in bio)),
            "mean_2c_pg": round(sum(c2_vals) / len(c2_vals), 4) if c2_vals else None, "breeding_system_summary": "; ".join(bs_vals) if bs_vals else None,
        })
    for r in linked: r.pop("_bio", None)
    return linked, summaries


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    p = argparse.ArgumentParser(description="Download and cross-link public Packera cytological records.")
    p.add_argument("--cytology-dir", default=str(PROJECT_ROOT / "data" / "cytology"))
    p.add_argument("--sources", nargs="+", default=["ccdb", "icna", "fna", "labels"], choices=["ccdb", "icna", "fna", "labels"])
    p.add_argument("--local-cytology-csv", default=None); p.add_argument("--sync-db-only", action="store_true")
    p.add_argument("--db", default=None); p.add_argument("--no-db", action="store_true"); p.add_argument("--force", action="store_true")
    args = p.parse_args(argv); cyto_dir = Path(args.cytology_dir).resolve(); cyto_dir.mkdir(parents=True, exist_ok=True)
    rec_dest, sum_dest, man_dest = (_ensure_flat_filename(cyto_dir, n) for n in ("public_packera_cytology_records.csv", "public_packera_taxon_cytology_summary.csv", "public_packera_cytology_manifest.json"))
    eff_db = None if args.no_db else Path(resolve_db_path(db_path=args.db, workspace_root=PROJECT_ROOT))
    smap, ctx = _load_sqlite_context(eff_db); raw_records: List[Dict[str, Any]] = []
    if args.sync_db_only:
        if rec_dest.is_file(): raw_records.extend(_read_local_csv(rec_dest, smap, reload_cache=True))
        if ctx["vouchers"]: raw_records.extend(_scan_labels(smap, ctx))
    else:
        for src, fn in (("ccdb", _fetch_ccdb), ("icna", _fetch_icna), ("fna", _fetch_fna)):
            if src in args.sources: raw_records.extend(fn(smap))
        if "labels" in args.sources: raw_records.extend(_scan_labels(smap, ctx))
    if args.local_cytology_csv: raw_records.extend(_read_local_csv(Path(args.local_cytology_csv).resolve(), smap))
    linked, summaries = _cross_link_and_summarize(raw_records, ctx)
    if args.force or not rec_dest.exists() or linked or args.sync_db_only:
        atomic_write_csv_rows(rec_dest, RECORD_COLUMNS, linked)
        atomic_write_csv_rows(sum_dest, SUMMARY_COLUMNS, summaries)
        atomic_write_bytes(man_dest, (json.dumps({"generated_at_utc": datetime.now(timezone.utc).isoformat(), "sources": ["sync-db-only"] if args.sync_db_only else list(args.sources), "n_records": len(linked), "n_taxa": len(summaries)}, indent=2) + "\n").encode("utf-8"))
    if not args.no_db and eff_db is not None:
        conn = get_db_connection(db_path=eff_db, init_schema=True, workspace_root=PROJECT_ROOT)
        try:
            if linked and (args.force or args.sync_db_only or set(args.sources) == {"ccdb", "icna", "fna", "labels"}):
                conn.execute("DELETE FROM cytological_records"); conn.execute("DELETE FROM taxon_cytology_summary")
            upsert_cytological_records(conn, linked, demo=False); upsert_taxon_cytology_summary(conn, summaries, demo=False)
        finally: conn.close()
    logger.info("Cytology staging complete: %d records across %d taxa -> %s", len(linked), len(summaries), cyto_dir)
    return 0


if __name__ == "__main__": raise SystemExit(main())
