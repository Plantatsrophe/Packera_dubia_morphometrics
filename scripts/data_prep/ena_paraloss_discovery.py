"""Dynamic ENA taxon-subtree SRA and Compositae-ParaLoss-1272 discovery helpers."""

from __future__ import annotations

import csv
import io
import json
import sys
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

_UA = "Plantmatics/1.0 (https://github.com/plantmatics; mailto:plantmatics@users.noreply.github.com)"
_DEFAULT_PACKERA_TAXID = "189364"
_PARALOSS_DIRECT_URLS = ("https://figshare.com/ndownloader/files/42046956",)
_PARALOSS_FIGSHARE_ARTICLES = ("https://api.figshare.com/v2/articles/23982312", "https://api.figshare.com/v2/articles/24043458")
SRA_MANIFEST_COLS = ("bioproject", "run_accession", "sample_accession", "scientific_name", "library_name", "specimen_voucher", "fastq_ftp", "fastq_bytes")
PARALOSS_COLS = ("locus_id", "gene_family", "paralogy_class", "target_length_bp", "locus_paralogy_weight")


def _default_http_get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": _UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _to_text(payload: Union[bytes, str]) -> str:
    return payload if isinstance(payload, str) else payload.decode("utf-8", errors="replace")


def _to_bytes(payload: Union[bytes, str]) -> bytes:
    return payload.encode("utf-8") if isinstance(payload, str) else payload


def resolve_ena_ncbi_taxid(
    target_genus: str,
    http_get_fn: Optional[Callable[..., Union[bytes, str]]] = None,
) -> Optional[str]:
    """Resolve NCBI/ENA taxon ID for `target_genus`, with offline fallback for Packera."""
    genus = str(target_genus or "").strip()
    if not genus:
        return None
    getter = http_get_fn or _default_http_get
    q_genus = urllib.parse.quote(genus)
    try:
        raw = _to_text(getter(f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{q_genus}", timeout=25))
        data = json.loads(raw)
        items = data if isinstance(data, list) else ([data] if isinstance(data, dict) else [])
        for item in items:
            tid = str(item.get("taxId") or item.get("tax_id") or item.get("id") or "").strip()
            if tid.isdigit():
                return tid
    except Exception:
        pass
    try:
        term = urllib.parse.quote(f"{genus}[Scientific Name]")
        raw = _to_text(getter(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={term}&retmode=json", timeout=25))
        ids = json.loads(raw).get("esearchresult", {}).get("idlist") or []
        for tid in ids:
            if str(tid).strip().isdigit():
                return str(tid).strip()
    except Exception:
        pass
    if genus.lower() == "packera":
        return _DEFAULT_PACKERA_TAXID
    return None


def _parse_ena_tsv_rows(tsv_text: str, default_bioproject: str = "") -> List[Dict[str, str]]:
    out: List[Dict[str, str]] = []
    for row in csv.DictReader(tsv_text.splitlines(), delimiter="\t"):
        run_acc = str(row.get("run_accession") or row.get("Run") or "").strip()
        if not run_acc:
            continue
        bp = str(row.get("study_accession") or row.get("bioproject") or default_bioproject).strip()
        out.append({
            "bioproject": bp,
            "run_accession": run_acc,
            "sample_accession": str(row.get("sample_accession") or "").strip(),
            "scientific_name": str(row.get("scientific_name") or "").strip(),
            "library_name": str(row.get("library_name") or "").strip(),
            "specimen_voucher": str(row.get("specimen_voucher") or "").strip(),
            "fastq_ftp": str(row.get("fastq_ftp") or "").strip(),
            "fastq_bytes": str(row.get("fastq_bytes") or "").strip(),
        })
    return out


def fetch_ena_tax_tree_sra_rows(
    target_genus: str,
    seed_bioprojects: Sequence[str],
    http_get_fn: Optional[Callable[..., Union[bytes, str]]] = None,
    ncbi_get_fn: Optional[Callable[[str, Dict[str, Any]], Union[bytes, str]]] = None,
) -> List[Dict[str, str]]:
    """Query ENA tax_tree read runs and merge curated seed BioProjects cleanly."""
    getter = http_get_fn or _default_http_get
    fields = "study_accession,run_accession,sample_accession,scientific_name,library_name,specimen_voucher,fastq_ftp,fastq_bytes"
    rows: List[Dict[str, str]] = []
    seen: set[Tuple[str, str]] = set()

    def _add(candidates: Sequence[Dict[str, str]]) -> None:
        for rec in candidates:
            key = (rec["bioproject"], rec["run_accession"])
            if key not in seen:
                seen.add(key)
                rows.append(rec)

    tax_id = resolve_ena_ncbi_taxid(target_genus, http_get_fn=getter)
    if tax_id:
        try:
            url = f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_tree({tax_id})&fields={fields}&limit=0&format=tsv"
            _add(_parse_ena_tsv_rows(_to_text(getter(url, timeout=45))))
        except Exception as exc:
            print(f"[INFO] ENA tax_tree({tax_id}) query failed ({exc}); checking seed BioProjects.", file=sys.stderr)

    for prj in seed_bioprojects:
        bp = str(prj or "").strip()
        if not bp:
            continue
        bp_rows: List[Dict[str, str]] = []
        q_bp = urllib.parse.quote(bp)
        for url in (
            f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=study_accession={q_bp}&fields={fields}&limit=0&format=tsv",
            f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={q_bp}&result=read_run&fields={fields}&format=tsv",
        ):
            try:
                bp_rows = _parse_ena_tsv_rows(_to_text(getter(url, timeout=45)), default_bioproject=bp)
                if bp_rows:
                    break
            except Exception as exc:
                print(f"[INFO] ENA BioProject query failed for {bp} ({exc}); trying fallback.", file=sys.stderr)
        if not bp_rows and ncbi_get_fn is not None:
            try:
                es = json.loads(_to_text(ncbi_get_fn("esearch.fcgi", {
                    "db": "sra", "term": f"{bp}[BioProject]", "retmax": 400, "usehistory": "y", "retmode": "json",
                }))).get("esearchresult", {})
                ids = [str(x) for x in (es.get("idlist") or [])]
                ef_p: Dict[str, Any] = {"db": "sra", "rettype": "runinfo", "retmode": "text"}
                if ids:
                    ef_p["id"] = ",".join(ids[:400])
                elif es.get("webenv") and es.get("querykey"):
                    ef_p.update({"WebEnv": es["webenv"], "query_key": es["querykey"], "retmax": 400})
                if "id" in ef_p or "WebEnv" in ef_p:
                    ef_txt = _to_text(ncbi_get_fn("efetch.fcgi", ef_p))
                    for r in csv.DictReader(ef_txt.splitlines()):
                        acc = str(r.get("Run") or r.get("run_accession") or "").strip()
                        if acc:
                            bp_rows.append({
                                "bioproject": str(r.get("BioProject") or bp).strip(),
                                "run_accession": acc,
                                "sample_accession": str(r.get("Sample") or r.get("BioSample") or "").strip(),
                                "scientific_name": str(r.get("ScientificName") or "").strip(),
                                "library_name": str(r.get("LibraryName") or "").strip(),
                                "specimen_voucher": str(r.get("specimen_voucher") or "").strip(),
                                "fastq_ftp": str(r.get("fastq_ftp") or r.get("download_path") or "").strip(),
                                "fastq_bytes": str(r.get("fastq_bytes") or r.get("size_MB") or "").strip(),
                            })
            except Exception as exc:
                print(f"[WARN] NCBI SRA fallback failed for {bp}: {exc}", file=sys.stderr)
        _add(bp_rows)
    return rows


def _classify_paralogy(raw_class: str, default_single: bool = True) -> Tuple[str, str]:
    low = str(raw_class or "").strip().lower()
    if any(tok in low for tok in ("multi", "paralog", "duplic", "ambig", "flag", "non_single", "false")):
        return "multi_copy", "0.5"
    if any(tok in low for tok in ("single", "ortholog", "paraloss", "clean", "pass", "core", "true")) or (not low and default_single):
        return "single_copy", "1.5"
    return "neutral", "1.0"


def _normalize_locus_dict(rec: Dict[str, Any], default_single: bool = True) -> Optional[Dict[str, str]]:
    lowered = {str(k).strip().lower(): str(v or "").strip() for k, v in rec.items() if k is not None}
    locus_id = next((lowered[k] for k in ("locus_id", "locus", "gene_id", "target_id", "id", "name", "marker") if lowered.get(k)), "")
    if not locus_id:
        return None
    gene_family = next((lowered[k] for k in ("gene_family", "family", "orthogroup", "annotation", "function", "gene") if lowered.get(k)), "Compositae-ParaLoss-1272")
    raw_cls = next((lowered[k] for k in ("paralogy_class", "paralogy", "copy_status", "status", "class", "is_paralog", "paralogs") if lowered.get(k)), "")
    p_cls, weight = _classify_paralogy(raw_cls, default_single=default_single)
    raw_len = next((lowered[k] for k in ("target_length_bp", "length_bp", "length", "bp", "len", "sequence_length") if lowered.get(k)), "")
    length_bp = raw_len if raw_len.isdigit() else ""
    return {
        "locus_id": locus_id,
        "gene_family": gene_family,
        "paralogy_class": p_cls,
        "target_length_bp": length_bp,
        "locus_paralogy_weight": weight,
    }


def _parse_paraloss_payload(payload: bytes) -> List[Dict[str, str]]:
    if not payload:
        return []
    if payload[:2] == b"PK":
        rows: List[Dict[str, str]] = []
        with zipfile.ZipFile(io.BytesIO(payload)) as zf:
            for name in sorted(zf.namelist()):
                if name.endswith("/") or "__MACOSX" in name or Path(name).name.startswith("._"):
                    continue
                if any(name.lower().endswith(ext) for ext in (".tsv", ".csv", ".txt", ".fasta", ".fa", ".fna", ".json")):
                    rows.extend(_parse_paraloss_payload(zf.read(name)))
        return rows
    text = _to_text(payload).strip()
    if not text:
        return []
    if text[0] in "{[":
        obj = json.loads(text)
        items = obj if isinstance(obj, list) else (obj.get("loci") or obj.get("records") or [])
        return [n for item in items if isinstance(item, dict) and (n := _normalize_locus_dict(item))]
    if text.startswith(">"):
        out, cur_id, cur_hdr, cur_len = [], "", "", 0
        for line in text.splitlines():
            s = line.strip()
            if s.startswith(">"):
                if cur_id:
                    p_cls, wt = _classify_paralogy(cur_hdr, default_single=True)
                    out.append({"locus_id": cur_id, "gene_family": "Compositae-ParaLoss-1272", "paralogy_class": p_cls, "target_length_bp": str(cur_len), "locus_paralogy_weight": wt})
                hdr = s[1:].strip()
                cur_id = hdr.split()[0].split("|")[0].split(";")[0]
                cur_hdr, cur_len = hdr, 0
            else:
                cur_len += len(s)
        if cur_id:
            p_cls, wt = _classify_paralogy(cur_hdr, default_single=True)
            out.append({"locus_id": cur_id, "gene_family": "Compositae-ParaLoss-1272", "paralogy_class": p_cls, "target_length_bp": str(cur_len), "locus_paralogy_weight": wt})
        return out
    first_line = text.splitlines()[0]
    delim = "\t" if "\t" in first_line else ","
    parsed: List[Dict[str, str]] = []
    for row in csv.DictReader(text.splitlines(), delimiter=delim):
        if norm := _normalize_locus_dict(row):
            parsed.append(norm)
    return parsed


def fetch_compositae_paraloss1272_catalog(
    genomics_dir: Path,
    target_family: str,
    atomic_save_fn: Callable[..., Any],
    http_get_fn: Optional[Callable[..., Union[bytes, str]]] = None,
    force: bool = False,
) -> int:
    """Fetch and stage flat Compositae-ParaLoss-1272 single-copy locus metadata for Asteraceae."""
    if str(target_family or "").strip().lower() != "asteraceae":
        return 0
    genomics_dir = Path(genomics_dir)
    out_path = genomics_dir / "public_packera_paraloss1272_loci.csv"
    if out_path.is_file() and out_path.stat().st_size > 0 and not force:
        try:
            with open(out_path, encoding="utf-8") as fh:
                existing = sum(1 for _ in csv.DictReader(fh))
            if existing > 0:
                return existing
        except Exception:
            pass
    getter = http_get_fn or _default_http_get
    candidate_urls: List[str] = list(_PARALOSS_DIRECT_URLS)
    for api_url in _PARALOSS_FIGSHARE_ARTICLES:
        try:
            meta = json.loads(_to_text(getter(api_url, timeout=25)))
            files = meta.get("files") if isinstance(meta, dict) else None
            if isinstance(files, list):
                ranked = sorted(
                    [f for f in files if isinstance(f, dict) and f.get("download_url")],
                    key=lambda f: (0 if any(k in str(f.get("name", "")).lower() for k in ("1272", "paraloss", "locus", "loci", "target")) else 1, str(f.get("name", ""))),
                )
                candidate_urls.extend(str(f["download_url"]) for f in ranked)
        except Exception:
            continue
    dedup: Dict[str, Dict[str, str]] = {}
    for url in dict.fromkeys(candidate_urls):
        try:
            rows = _parse_paraloss_payload(_to_bytes(getter(url, timeout=45)))
        except Exception:
            continue
        for rec in rows:
            dedup[rec["locus_id"]] = rec
        if dedup:
            break
    if not dedup:
        print("[WARN] Compositae-ParaLoss-1272 metadata download unavailable (continuing offline).", file=sys.stderr)
        return 0
    buf = io.StringIO()
    wr = csv.DictWriter(buf, fieldnames=list(PARALOSS_COLS))
    wr.writeheader()
    wr.writerows(dedup.values())
    atomic_save_fn(genomics_dir, "public_packera_paraloss1272_loci.csv", buf.getvalue().encode("utf-8"), force=force)
    return len(dedup)
