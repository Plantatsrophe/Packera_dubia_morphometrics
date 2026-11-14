#!/usr/bin/env python3
"""
Module: scripts/data_prep/download_public_packera_sequences.py
Project: Plantmatics - Packera Species Delimitation & Morphometrics Pipeline
Description: Retrieves flat public Packera sequences, phylogenomic trees, SRA/ENA
  manifests, and BOLD barcodes into `data/genomics/public_packera_*` via tmpdir staging.
"""

from __future__ import annotations

import argparse, csv, gzip, hashlib, importlib, importlib.util, json, os, re, shutil, sys, tempfile, time, zipfile
import urllib.error, urllib.parse, urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

_UA = "Plantmatics/1.0 (https://github.com/plantmatics; mailto:plantmatics@users.noreply.github.com)"
_DRYAD_DOI = "doi%3A10.5061%2Fdryad.bg79cnppn"
_DRYAD_ZIPS = ("speciestree_files.zip", "PhyTop_results.zip", "RelTime_trees.zip", "alignment_files.zip")
_SRA_BIOPROJECTS = ("PRJNA907383", "PRJNA978568", "PRJNA540287", "PRJNA516161")
_NCBI_TERM = (
    '("Packera"[Organism] OR "Senecio tomentosus"[Organism] OR "Senecio aureus"[Organism] OR '
    '"Senecio anonymus"[Organism] OR "Senecio pauperculus"[Organism] OR "Senecio schweinitzianus"[Organism] OR '
    '"Senecio millefolium"[Organism] OR "Senecio obovatus"[Organism] OR "Senecio plattensis"[Organism] OR '
    '"Senecio antennariifolius"[Organism] OR "Senecio glabellus"[Organism]) AND 120:200000[SLEN]'
)
_SYNONYM_MAP: Dict[str, str] = {
    "senecio tomentosus": "Packera_dubia", "packera tomentosa": "Packera_dubia", "senecio aureus": "Packera_aurea",
    "senecio anonymus": "Packera_anonyma", "senecio pauperculus": "Packera_paupercula",
    "senecio schweinitzianus": "Packera_schweinitziana", "senecio millefolium": "Packera_millefolium",
    "senecio obovatus": "Packera_obovata", "senecio plattensis": "Packera_plattensis",
    "senecio antennariifolius": "Packera_antennariifolia", "senecio glabellus": "Packera_glabella",
}
_DRYAD_EXACT_MAP: Dict[str, str] = {
    "18complete_packera_plastomes.fasta": "public_packera_18complete_plastomes.fasta", "111mapping_packera_plastomes.fasta": "public_packera_111mapping_plastomes.fasta",
    "mo_astral_lpp.tre": "public_packera_MO_Astral_lpp.tre", "astralpro_astral_lpp.tre": "public_packera_ASTRALPro_Astral_lpp.tre",
    "wastral_astral_lpp.tre": "public_packera_wASTRAL_Astral_lpp.tre", "111mapping_packera_plastomes.tre": "public_packera_111mapping_plastomes.tre",
    "mo_concat_best.fasta": "public_packera_MO_concat_best.fasta", "1-to-1_concat_best.fasta": "public_packera_1to1_concat_best.fasta",
}
_LAST_NCBI_CALL = [0.0]


class _CrossHostAuthStripRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Strip Authorization headers on redirects across domains (e.g. Dryad -> S3)."""
    def redirect_request(self, req: Any, fp: Any, code: int, msg: str, headers: Any, newurl: str) -> Any:
        new_req = super().redirect_request(req, fp, code, msg, headers, newurl)
        if new_req is not None and urllib.parse.urlparse(req.full_url).netloc.lower() != urllib.parse.urlparse(newurl).netloc.lower():
            for hdrs in (new_req.headers, new_req.unredirected_hdrs):
                for k in [k for k in hdrs if k.lower() == "authorization"]:
                    hdrs.pop(k, None)
        return new_req


def _ensure_flat_filename(genomics_dir: Path, filename: str) -> Path:
    if not filename or "/" in filename or "\\" in filename or Path(filename).name != filename:
        raise ValueError(f"Refusing non-flat path inside genomics directory: {filename!r}")
    if not filename.startswith("public_packera_"):
        raise ValueError(f"Public genomics filename must start with 'public_packera_': {filename!r}")
    g_res = genomics_dir.resolve()
    dest = g_res / filename
    if dest.parent != g_res or (dest.exists() and dest.is_dir()):
        raise ValueError(f"Refusing subdirectory inside genomics directory: {dest}")
    return genomics_dir / filename


def _atomic_save(genomics_dir: Path, filename: str, source: Union[bytes, Path], force: bool = False) -> Tuple[Path, bool]:
    dest = _ensure_flat_filename(genomics_dir, filename)
    if not force and dest.exists() and dest.is_file() and dest.stat().st_size > 0:
        return dest, False
    genomics_dir.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=".tmp_pub_", suffix=".tmp", dir=str(genomics_dir))
    try:
        with os.fdopen(fd, "wb") as dst_fh:
            if isinstance(source, Path):
                with open(source, "rb") as src_fh:
                    shutil.copyfileobj(src_fh, dst_fh, length=65536)
            else:
                dst_fh.write(source)
            dst_fh.flush()
            os.fsync(dst_fh.fileno())
        os.replace(tmp_name, dest)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise
    return dest, True


def _http_request(url: str, *, headers: Optional[Dict[str, str]] = None, timeout: int = 60,
                  opener: Optional[urllib.request.OpenerDirector] = None, data: Optional[bytes] = None, dest: Optional[Path] = None) -> bytes:
    req = urllib.request.Request(url, data=data, headers={"User-Agent": _UA, **(headers or {})}, method="POST" if data is not None else "GET")
    open_fn = opener.open if opener is not None else urllib.request.urlopen
    with open_fn(req, timeout=timeout) as resp:
        if dest is not None:
            with open(dest, "wb") as out:
                shutil.copyfileobj(resp, out, length=65536)
            return b""
        raw = resp.read()
    return gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw


def _ncbi_get(endpoint: str, params: Dict[str, Any]) -> bytes:
    delay = 0.12 if os.environ.get("NCBI_API_KEY", "").strip() else 0.35
    now = time.monotonic()
    if _LAST_NCBI_CALL[0] > 0.0 and (now - _LAST_NCBI_CALL[0]) < delay:
        time.sleep(delay - (now - _LAST_NCBI_CALL[0]))
    _LAST_NCBI_CALL[0] = time.monotonic()
    q = {**params, "tool": "plantmatics", "email": "plantmatics@users.noreply.github.com"}
    if api_key := os.environ.get("NCBI_API_KEY", "").strip():
        q["api_key"] = api_key
    return _http_request(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}?{urllib.parse.urlencode(q)}", timeout=90)


def _normalize_taxon_binomial(text: str) -> str:
    hay = text.replace("_", " ").lower()
    for syn, canon in _SYNONYM_MAP.items():
        if syn in hay:
            return canon
    if m := re.search(r"\bpackera\s+[0-9]*([a-z-]+)\b", hay):
        return "Packera_dubia" if m.group(1) == "tomentosa" else f"Packera_{m.group(1)}"
    return "Packera_sp"

def _parse_fasta_text(text: str) -> List[Tuple[str, str]]:
    recs, hdr, buf = [], None, []
    for ln in text.splitlines():
        if not (s := ln.strip()):
            continue
        if s.startswith(">"):
            if hdr is not None and buf:
                recs.append((hdr, "".join(buf).upper()))
            hdr, buf = s, []
        else:
            buf.append(s)
    if hdr is not None and buf:
        recs.append((hdr, "".join(buf).upper()))
    return recs

def _classify_ncbi_record(raw_hdr: str, seq: str) -> Tuple[str, str]:
    toks = raw_hdr.lstrip(">").strip().split(None, 1)
    acc, orig_desc = (toks[0] if toks else "UNKNOWN"), (toks[1] if len(toks) > 1 else "")
    out_hdr = f">{acc}|{_normalize_taxon_binomial(f'{acc} {orig_desc}')} {orig_desc}".rstrip()
    hlow = raw_hdr.lower()
    if len(seq) >= 100000 or re.search(r"plastom|chloroplast|plastid|rbcl|matk|psba|trnh|trnl|trnf|ndhf|rpl16|rps16|trnk", hlow):
        cat = "cp"
    elif re.search(r"\bits\b|its1|5\.8s|its2|\bets\b|18s|26s|internal transcribed|external transcribed|ribosomal", hlow):
        cat = "its"
    elif re.search(r"waxy|gbssi|granule-bound starch|starch synthase|nuclear", hlow):
        cat = "waxy"
    elif re.search(r"cpdna|atpb|ycf|rpoc|rpoa|psb|ndh|trn|pet|rpl|rps|accd|ccsa", hlow) or len(seq) >= 10000:
        cat = "cp"
    else:
        cat = "waxy"
    return cat, f"{out_hdr}\n{seq}\n"

def _fetch_ncbi_nuccore(genomics_dir: Path, force: bool) -> Dict[str, int]:
    es = json.loads(_ncbi_get("esearch.fcgi", {"db": "nuccore", "term": _NCBI_TERM, "retmax": 400, "usehistory": "y", "retmode": "json"}).decode("utf-8", errors="replace")).get("esearchresult", {})
    total, webenv, qkey = int(es.get("count") or 0), str(es.get("webenv") or ""), str(es.get("querykey") or "")
    idlist, chunks = [str(x) for x in (es.get("idlist") or [])], []
    if webenv and qkey and total > 0:
        for start in range(0, min(total, 2000), 200):
            chunks.append(_ncbi_get("efetch.fcgi", {"db": "nuccore", "WebEnv": webenv, "query_key": qkey, "retstart": start, "retmax": 200, "rettype": "fasta", "retmode": "text"}).decode("utf-8", errors="replace"))
    elif idlist:
        for start in range(0, len(idlist), 200):
            chunks.append(_ncbi_get("efetch.fcgi", {"db": "nuccore", "id": ",".join(idlist[start : start + 200]), "rettype": "fasta", "retmode": "text"}).decode("utf-8", errors="replace"))
    buckets: Dict[str, List[str]] = {"cp": [], "its": [], "waxy": []}
    for raw_hdr, seq in _parse_fasta_text("\n".join(chunks)):
        cat, entry = _classify_ncbi_record(raw_hdr, seq)
        buckets[cat].append(entry)
    counts: Dict[str, int] = {}
    for key, fname, ckey in (("cp", "public_packera_plastomes_and_cpdna.fasta", "nuccore_plastid"), ("its", "public_packera_its_ets.fasta", "nuccore_its"), ("waxy", "public_packera_waxy_nuclear.fasta", "nuccore_waxy")):
        counts[ckey] = len(buckets[key])
        if buckets[key]:
            _atomic_save(genomics_dir, fname, "".join(buckets[key]).encode("utf-8"), force=force)
    return counts


def _fetch_sra_manifest(genomics_dir: Path, tmpd: Path, force: bool) -> List[Dict[str, str]]:
    cols = ("bioproject", "run_accession", "sample_accession", "scientific_name", "library_name", "fastq_ftp", "fastq_bytes")
    rows, seen = [], set()
    for prj in _SRA_BIOPROJECTS:
        prj_rows: List[Dict[str, str]] = []
        try:
            ena_url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={prj}&result=read_run&fields=run_accession,sample_accession,scientific_name,library_name,fastq_ftp,fastq_bytes&format=tsv"
            tsv = _http_request(ena_url, timeout=45).decode("utf-8", errors="replace")
            prj_rows = [{"bioproject": prj, **{k: str(r.get(k) or "").strip() for k in cols[1:]}}
                        for r in csv.DictReader(tsv.splitlines(), delimiter="\t") if r.get("run_accession")]
        except Exception as exc:
            print(f"[INFO] ENA query failed for {prj} ({exc}); falling back to NCBI SRA.", file=sys.stderr)
        if not prj_rows:
            try:
                es = json.loads(_ncbi_get("esearch.fcgi", {"db": "sra", "term": f"{prj}[BioProject]", "retmax": 400, "usehistory": "y", "retmode": "json"}).decode("utf-8", errors="replace")).get("esearchresult", {})
                ids = [str(x) for x in (es.get("idlist") or [])]
                ef_p: Dict[str, Any] = {"db": "sra", "rettype": "runinfo", "retmode": "text"}
                if ids:
                    ef_p["id"] = ",".join(ids[:400])
                elif es.get("webenv") and es.get("querykey"):
                    ef_p.update({"WebEnv": es["webenv"], "query_key": es["querykey"], "retmax": 400})
                if "id" in ef_p or "WebEnv" in ef_p:
                    for r in csv.DictReader(_ncbi_get("efetch.fcgi", ef_p).decode("utf-8", errors="replace").splitlines()):
                        if acc := str(r.get("Run") or r.get("run_accession") or "").strip():
                            prj_rows.append({"bioproject": str(r.get("BioProject") or prj).strip(), "run_accession": acc,
                                             "sample_accession": str(r.get("Sample") or r.get("BioSample") or "").strip(),
                                             "scientific_name": str(r.get("ScientificName") or "").strip(), "library_name": str(r.get("LibraryName") or "").strip(),
                                             "fastq_ftp": str(r.get("fastq_ftp") or r.get("download_path") or "").strip(), "fastq_bytes": str(r.get("fastq_bytes") or r.get("size_MB") or "").strip()})
            except Exception as exc:
                print(f"[WARN] NCBI SRA fallback failed for {prj}: {exc}", file=sys.stderr)
        for r in prj_rows:
            if (r["bioproject"], r["run_accession"]) not in seen:
                seen.add((r["bioproject"], r["run_accession"]))
                rows.append(r)

    tmp_csv = tmpd / "sra_runs.csv"
    with open(tmp_csv, "w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(cols))
        wr.writeheader()
        wr.writerows(rows)
    _atomic_save(genomics_dir, "public_packera_sra_runs.csv", tmp_csv, force=force)
    return rows


def _skim_sra_accessions(accessions: Sequence[str], sra_rows: Sequence[Dict[str, str]], genomics_dir: Path, tmpd: Path) -> int:
    targets = {a.strip().upper() for a in accessions if a and a.strip()}
    if not targets:
        return 0
    rows = list(sra_rows)
    if not rows and (genomics_dir / "public_packera_sra_runs.csv").is_file():
        with open(genomics_dir / "public_packera_sra_runs.csv", "rt", encoding="utf-8", errors="replace") as fh:
            rows = list(csv.DictReader(fh))
    skim_dir = tmpd / "skim_fastqs"
    skim_dir.mkdir(parents=True, exist_ok=True)
    downloaded = 0
    for r in rows:
        if (acc := str(r.get("run_accession", "")).strip()).upper() not in targets:
            continue
        for idx, raw_u in enumerate(str(r.get("fastq_ftp", "")).split(";")):
            if not (u := raw_u.strip()):
                continue
            dl_url = u if u.startswith(("http://", "https://")) else f"https://{u.lstrip('/')}"
            fname = Path(urllib.parse.urlparse(dl_url).path).name or f"{acc}_{idx}.fastq.gz"
            try:
                _http_request(dl_url, timeout=180, dest=skim_dir / fname)
                downloaded += 1
            except Exception as exc:
                print(f"[WARN] Failed streaming SRA FASTQ {dl_url}: {exc}", file=sys.stderr)
    if downloaded > 0:
        if (pr := str(Path(__file__).resolve().parents[2])) not in sys.path:
            sys.path.insert(0, pr)
        importlib.import_module("scripts.data_prep.skim_offtarget_organellar").run_offtarget_organellar_skim(skim_dir, genomics_dir)
    return downloaded


def _obtain_dryad_token() -> str:
    if tok := os.environ.get("DRYAD_API_TOKEN", "").strip():
        return tok
    cid, csec = os.environ.get("DRYAD_CLIENT_ID", "").strip(), os.environ.get("DRYAD_CLIENT_SECRET", "").strip()
    if cid and csec:
        body = urllib.parse.urlencode({"client_id": cid, "client_secret": csec, "grant_type": "client_credentials"}).encode("utf-8")
        return str(json.loads(_http_request("https://datadryad.org/oauth/token", headers={"Content-Type": "application/x-www-form-urlencoded"}, data=body).decode("utf-8")).get("access_token") or "").strip()
    return ""


def _route_dryad_file(src: Path, zip_name: str) -> Optional[str]:
    low, zlow = src.name.lower(), zip_name.lower()
    if "__macosx" in str(src).lower() or src.name.startswith("._"):
        return None
    if low in _DRYAD_EXACT_MAP:
        return _DRYAD_EXACT_MAP[low]
    if ("phytop" in zlow or "phytop" in low) and low.endswith((".tsv", ".csv", ".tree", ".tre", ".nwk")):
        stem_clean = re.sub(r"_phytop(\.info)?$", "", src.stem, flags=re.IGNORECASE)
        method = (
            "ASTRALPro" if ("astralpro" in low or "astral-pro" in low)
            else "wASTRAL" if "wastral" in low
            else "HybPiper-df" if "hybpiper" in low
            else "MO" if re.search(r"(^|[^a-z0-9])mo([^a-z0-9]|$)", low)
            else (re.sub(r"[^A-Za-z0-9_-]+", "_", stem_clean).strip("_") or "table")
        )
        ext = ".tre" if low.endswith((".tree", ".tre", ".nwk")) else ".tsv"
        return f"public_packera_PhyTop_{method}{ext}"
    if ("reltime" in zlow or "reltime" in low) and low.endswith((".tre", ".nwk", ".tree")):
        return f"public_packera_RelTime_{re.sub(r'[^A-Za-z0-9_.-]+', '_', src.stem).strip('_') or 'tree'}.tre"
    return None


def _fetch_dryad(genomics_dir: Path, tmpd: Path, dryad_zip_dir: Optional[Path], force: bool) -> Tuple[str, int]:
    manifest_files: List[Dict[str, Any]] = []
    manifest_path = genomics_dir / "public_packera_dryad_manifest.json"
    try:
        v_raw = _http_request(f"https://datadryad.org/api/v2/datasets/{_DRYAD_DOI}/versions", timeout=45)
        if v_list := json.loads(v_raw.decode("utf-8")).get("_embedded", {}).get("stash:versions", []):
            v_href = v_list[-1]["_links"]["self"]["href"]
            f_raw = _http_request(urllib.parse.urljoin("https://datadryad.org", v_href.rstrip("/") + "/files?per_page=100"), timeout=45)
            for item in json.loads(f_raw.decode("utf-8")).get("_embedded", {}).get("stash:files", []):
                dl_href = item.get("_links", {}).get("stash:download", {}).get("href") or (f"/api/v2/files/{item['id']}/download" if item.get("id") is not None else "")
                manifest_files.append({k: item.get(k, "") for k in ("id", "path", "size", "digestType", "digest")} | {"download_href": dl_href})
            _atomic_save(genomics_dir, "public_packera_dryad_manifest.json", json.dumps({"doi": "10.5061/dryad.bg79cnppn", "files": manifest_files}, indent=2).encode("utf-8"), force=force)
    except Exception as exc:
        print(f"[WARN] Dryad manifest fetch failed ({exc}).", file=sys.stderr)
        if manifest_path.is_file():
            try:
                manifest_files = json.loads(manifest_path.read_text(encoding="utf-8")).get("files", [])
            except Exception:
                manifest_files = []

    by_name = {Path(str(m.get("path", ""))).name: m for m in manifest_files if m.get("path")}
    staged_zips, status = [], "manifest_only_auth_required"
    if dryad_zip_dir is not None and dryad_zip_dir.is_dir():
        for zname in _DRYAD_ZIPS:
            cand = next((p for p in dryad_zip_dir.iterdir() if p.is_file() and p.name.lower() == zname.lower()), None)
            if cand is not None:
                shutil.copy2(cand, tmpd / zname)
                staged_zips.append(tmpd / zname)
        if staged_zips:
            status = "local_zip_dir"

    if not staged_zips:
        try:
            if token := _obtain_dryad_token():
                opener = urllib.request.build_opener(_CrossHostAuthStripRedirectHandler())
                for zname in _DRYAD_ZIPS:
                    if (meta := by_name.get(zname)) and meta.get("download_href"):
                        _http_request(urllib.parse.urljoin("https://datadryad.org", str(meta["download_href"])), headers={"Authorization": f"Bearer {token}"}, timeout=180, opener=opener, dest=tmpd / zname)
                        staged_zips.append(tmpd / zname)
                if staged_zips:
                    status = "downloaded"
        except urllib.error.HTTPError as exc:
            if exc.code in (401, 403):
                staged_zips.clear()
            else:
                raise

    if not staged_zips:
        pw_avail = importlib.util.find_spec("playwright") is not None
        print(f'[INFO] Dryad archives require authentication (playwright_installed={pw_avail}). Set DRYAD_API_TOKEN and run: export DRYAD_API_TOKEN="<token>" && .venv/bin/python main.py download-public-genomics --genomics-dir data/genomics or pass --dryad-zip-dir <DIR>.')
        return status, 0

    n_extracted = 0
    for zp in staged_zips:
        meta = by_name.get(zp.name, {})
        if str(meta.get("digestType", "")).lower() in ("sha-256", "sha256") and meta.get("digest"):
            h = hashlib.sha256()
            with open(zp, "rb") as fh:
                for chunk in iter(lambda: fh.read(65536), b""):
                    h.update(chunk)
            if h.hexdigest().lower() != str(meta["digest"]).strip().lower():
                raise ValueError(f"SHA-256 mismatch for {zp.name}: expected {meta['digest']}, got {h.hexdigest()}")
        unpack_dir = tmpd / f"unpack_{zp.stem}"
        unpack_dir.mkdir(parents=True, exist_ok=True)
        unpack_res = unpack_dir.resolve()
        with zipfile.ZipFile(zp, "r") as zf:
            for member in zf.infolist():
                target = (unpack_dir / member.filename).resolve()
                if target != unpack_res and unpack_res not in target.parents:
                    raise ValueError(f"Unsafe zip-slip entry in {zp.name}: {member.filename!r}")
            zf.extractall(unpack_dir)
        for fp in sorted(unpack_dir.rglob("*")):
            if fp.is_file() and (out_name := _route_dryad_file(fp, zp.name)):
                _atomic_save(genomics_dir, out_name, fp, force=force)
                n_extracted += 1
    return status, n_extracted


def _format_bold_fasta(raw_text: str) -> Tuple[str, int]:
    if raw_text.lstrip().startswith(">"):
        out = []
        for hdr, seq in _parse_fasta_text(raw_text):
            if seq_c := seq.replace("-", "").strip():
                parts = [p.strip() for p in hdr.lstrip(">").split("|")]
                mk = next((p for p in parts[1:] if re.search(r"^(its|rbcl|matk|psba|trnh)", p, re.I)), "")
                out.append(f">{parts[0].split()[0]}|{_normalize_taxon_binomial(hdr)}{f'|{mk}' if mk else ''}\n{seq_c}\n")
        return "".join(out), len(out)
    try:
        parsed = json.loads(raw_text)
        objs = parsed if isinstance(parsed, list) else next((parsed[k] for k in ("records", "data", "documents") if isinstance(parsed, dict) and isinstance(parsed.get(k), list)), [])
    except json.JSONDecodeError:
        objs = [json.loads(ln) for ln in raw_text.splitlines() if ln.strip().startswith("{")]
    out_lines = []
    for rec in objs:
        if isinstance(rec, dict) and (seq_c := str(rec.get("nuc") or "").replace("-", "").strip().upper()):
            rid = str(rec.get("record_id") or rec.get("processid") or rec.get("id") or "BOLD").strip()
            tax = _normalize_taxon_binomial(str(rec.get("species") or rec.get("identification") or "Packera_sp"))
            mk = re.sub(r"\s+", "_", str(rec.get("marker_code") or rec.get("markercode") or "UNK").strip() or "UNK")
            out_lines.append(f">{rid}|{tax}|{mk}\n{seq_c}\n")
    return "".join(out_lines), len(out_lines)


def _fetch_bold(genomics_dir: Path, force: bool) -> int:
    try:
        q_raw = _http_request("https://portal.boldsystems.org/api/query?query=tax%3Agenus%3APackera&extent=full", timeout=60).decode("utf-8", errors="replace")
        if q_raw.lstrip().startswith(">"):
            fasta_txt, n_rec = _format_bold_fasta(q_raw)
        else:
            if not (qid := str(json.loads(q_raw).get("query_id") or "").strip()):
                raise ValueError("Missing query_id from BOLD v5 response")
            fasta_txt, n_rec = _format_bold_fasta(_http_request(f"https://portal.boldsystems.org/api/documents/{urllib.parse.quote(qid, safe='')}?length=2000", timeout=90).decode("utf-8", errors="replace"))
    except Exception as exc:
        print(f"[WARN] BOLD v5 query failed ({exc}); trying legacy BOLD API.", file=sys.stderr)
        try:
            fasta_txt, n_rec = _format_bold_fasta(_http_request("http://www.boldsystems.org/index.php/API_Public/sequence?taxon=Packera", timeout=60).decode("utf-8", errors="replace"))
        except Exception as exc2:
            print(f"[WARN] BOLD barcode retrieval unavailable: {exc2}", file=sys.stderr)
            return 0
    if fasta_txt:
        _atomic_save(genomics_dir, "public_packera_bold_barcodes.fasta", fasta_txt.encode("utf-8"), force=force)
    return n_rec


def download_public_packera_sequences(
    genomics_dir: Path, sources: Sequence[str] = ("ncbi", "dryad", "bold"),
    dryad_zip_dir: Optional[Path] = None, skim_sra_accessions: Optional[Sequence[str]] = None, force: bool = False,
) -> Dict[str, Any]:
    """Download public Packera reference sequences and phylogenomic trees into flat `genomics_dir`."""
    genomics_dir = Path(genomics_dir)
    genomics_dir.mkdir(parents=True, exist_ok=True)
    src_list = [s.strip().lower() for raw in sources for s in str(raw).split(",") if s.strip()] or ["ncbi", "dryad", "bold"]
    counts = {"nuccore_plastid": 0, "nuccore_its": 0, "nuccore_waxy": 0, "sra_runs": 0, "dryad_files": 0, "bold_records": 0, "sra_skimmed_fastqs": 0}
    dryad_status = "skipped"
    prev_sources: List[str] = []
    manifest_path = genomics_dir / "public_packera_sources_manifest.json"
    if manifest_path.is_file():
        try:
            prev = json.loads(manifest_path.read_text(encoding="utf-8"))
            if isinstance(prev.get("counts"), dict):
                counts.update({k: int(v) for k, v in prev["counts"].items() if k in counts and isinstance(v, int)})
            if "dryad" not in src_list and isinstance(prev.get("dryad_status"), str) and prev["dryad_status"] != "skipped":
                dryad_status = prev["dryad_status"]
            if isinstance(prev.get("sources_run"), list):
                prev_sources = [str(s).strip().lower() for s in prev["sources_run"] if str(s).strip()]
        except Exception:
            pass
    with tempfile.TemporaryDirectory(prefix="plantmatics_pub_") as tmp_str:
        tmpd, sra_rows = Path(tmp_str), []
        if "ncbi" in src_list:
            try:
                counts.update(_fetch_ncbi_nuccore(genomics_dir, force=force))
            except Exception as exc:
                print(f"[WARN] NCBI nuccore fetch failed: {exc}", file=sys.stderr)
            sra_rows = _fetch_sra_manifest(genomics_dir, tmpd, force=force)
            counts["sra_runs"] = len(sra_rows)
        if skim_sra_accessions:
            counts["sra_skimmed_fastqs"] = _skim_sra_accessions([a.strip() for raw in skim_sra_accessions for a in str(raw).split(",") if a.strip()], sra_rows, genomics_dir, tmpd)
        if "dryad" in src_list:
            dryad_status, counts["dryad_files"] = _fetch_dryad(genomics_dir, tmpd, dryad_zip_dir=dryad_zip_dir, force=force)
        if "bold" in src_list:
            counts["bold_records"] = _fetch_bold(genomics_dir, force=force)
    files_meta = {p.name: {"size_bytes": p.stat().st_size} for p in sorted(genomics_dir.iterdir()) if p.is_file() and p.name.startswith("public_packera_") and p.name != "public_packera_sources_manifest.json"}
    summary: Dict[str, Any] = {"timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "genomics_dir": str(genomics_dir), "sources_run": list(dict.fromkeys([*prev_sources, *src_list])), "dryad_status": dryad_status, "counts": counts, "files": files_meta}
    _atomic_save(genomics_dir, "public_packera_sources_manifest.json", json.dumps(summary, indent=2).encode("utf-8"), force=True)
    return summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Download flat public Packera reference sequences into data/genomics/.")
    parser.add_argument("--genomics-dir", type=Path, default=Path(__file__).resolve().parents[2] / "data" / "genomics")
    parser.add_argument("--sources", nargs="+", default=["ncbi", "dryad", "bold"], help="Sources to query (ncbi, dryad, bold)")
    parser.add_argument("--dryad-zip-dir", type=Path, default=None, help="Directory with manually downloaded Dryad .zip archives")
    parser.add_argument("--skim-sra-accessions", nargs="+", default=None, help="SRA/ENA run accessions to stream-skim via skim_offtarget_organellar")
    parser.add_argument("--force", action="store_true", default=False, help="Overwrite existing non-empty public_packera_* files")
    args = parser.parse_args(list(argv) if argv is not None else None)
    summary = download_public_packera_sequences(args.genomics_dir, sources=args.sources, dryad_zip_dir=args.dryad_zip_dir, skim_sra_accessions=args.skim_sra_accessions, force=bool(args.force))
    print(f"[INFO] Public Packera download summary: {json.dumps({'dryad_status': summary.get('dryad_status'), 'counts': summary.get('counts')})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
