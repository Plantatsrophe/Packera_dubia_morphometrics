#!/usr/bin/env python3
"""Hierarchical Eastern North American Packera monograph orchestrator."""

from __future__ import annotations

import argparse, csv, io, json, math, os, re, secrets, shutil, sqlite3, subprocess, sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon
from scripts.core.taxonomy_authority import evaluate_eastern_na_cohort_filter
from scripts.pipeline.resolve_focal_species_complex import (
    _extract_stage_metrics, _fetch_complex_rows_for_run, _read_csv_rows, _stage_run_ids,
    build_four_stage_commands, compile_complex_resolution_summary, slugify_taxon,
)

BINOMIAL_RE = re.compile(r"^Packera [a-z][a-z-]+$")
HYBRID_RE = re.compile(r"(?:\s+[x×]\s+|\s*×\s*)")
NA_TOKENS = {"", "NA", "NAN", "NONE", "NULL"}
MASTER_SYNTHESIS_COLS: Tuple[str, ...] = (
    "parent_run_id", "species_standardized", "species_subdir", "leiden_complex_id", "complex_neighbors",
    "n_raw_otus", "n_misid_purged", "pct_misid_purged", "n_hybrids_stage2", "hybrid_parent_pairs",
    "n_purified_intraspecific_otus", "intraspecific_optimal_k", "max_two_delta_bic", "n_subspecies_ecotype",
    "n_ecophenotypic_plasticity", "n_cohesive_cryptic_otus", "n_confirmed_type_control", "min_d_edaph",
    "max_abs_delta_doy_days", "median_coalescent_gdi", "public_backbone_clade", "n_tier0_abstained",
    "thinned_stage4_otus", "thinned_optimal_k", "primary_verdict",
)

def _write_atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    tmp.replace(path)

def _write_atomic_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, Any]]) -> None:
    buf = io.StringIO(newline="")
    writer = csv.DictWriter(buf, fieldnames=list(fieldnames), extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
    _write_atomic_text(path, buf.getvalue())

def _resolve_ro_db_path(ws: Path, db_path: Optional[Path]) -> Path:
    if db_path is not None:
        return Path(db_path)
    for name in ("plantmatics.sqlite", "packera_pipeline.sqlite"):
        if (cand := ws / "data" / name).is_file():
            return cand
    env_db = os.environ.get("PLANTMATICS_DB_PATH", "").strip() or os.environ.get("PACKERA_DB_PATH", "").strip()
    if env_db and ws.resolve() == PROJECT_ROOT.resolve():
        return Path(env_db)
    return ws / "data" / "plantmatics.sqlite"

def _open_ro_conn(db_file: Path) -> Optional[sqlite3.Connection]:
    if not db_file.is_file():
        return None
    try:
        conn = sqlite3.connect(f"file:{db_file.as_posix()}?mode=ro", uri=True)
        conn.row_factory = sqlite3.Row
        return conn
    except sqlite3.Error:
        return None

def _has_finite_float(val: Any) -> bool:
    try:
        return val is not None and not isinstance(val, bool) and str(val).strip().upper() not in NA_TOKENS and math.isfinite(float(val))
    except (ValueError, TypeError):
        return False

def _extract_uid_hex6(run_id: str) -> str:
    m = re.search(r"_([0-9a-f]{6})$", run_id.strip())
    return m.group(1) if m else ""

def generate_parent_run_id(prefix: str = "eastern_na") -> str:
    """Return a timestamped monograph parent run ID matching '<prefix>_YYYYMMDD_<hex6>'."""
    clean = (prefix or "eastern_na").strip().strip("/") or "eastern_na"
    return f"{clean}_{datetime.now().strftime('%Y%m%d')}_{secrets.token_hex(3)}"

def resolve_parent_run_id(
    workspace_root: Path, prefix: str = "eastern_na", parent_run_id: Optional[str] = None, resume: bool = False
) -> Tuple[str, str]:
    """Resolve `(parent_run_id, unique_id_hex6)` from explicit ID, resume symlink, or fresh generation."""
    if parent_run_id is not None and str(parent_run_id).strip():
        clean = str(parent_run_id).strip().strip("/")
        if clean.startswith("runs/"):
            clean = clean[len("runs/"):].strip("/")
        if not clean:
            raise ValueError("Provided --parent-run-id resolved to an empty string.")
        return clean, _extract_uid_hex6(clean)
    if resume:
        runs_root = (workspace_root / "runs").resolve()
        latest_link = runs_root / "latest_eastern_na"
        if not (latest_link.is_symlink() or latest_link.exists()):
            raise ValueError(f"Cannot --resume without --parent-run-id: {latest_link} is missing.")
        raw_t = Path(os.readlink(latest_link)) if latest_link.is_symlink() else latest_link
        resolved_dir = (latest_link.parent / raw_t).resolve() if not raw_t.is_absolute() else raw_t.resolve()
        if not resolved_dir.is_dir() or runs_root not in resolved_dir.parents:
            raise ValueError(f"Cannot --resume: {latest_link} target {resolved_dir} is invalid.")
        rel_id = resolved_dir.relative_to(runs_root).as_posix()
        return rel_id, _extract_uid_hex6(rel_id)
    fresh = generate_parent_run_id(prefix=prefix)
    return fresh, _extract_uid_hex6(fresh)

def _update_latest_eastern_na_symlink(workspace_root: Path, parent_run_id: str) -> None:
    runs_dir = workspace_root / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    tmp_link = runs_dir / f".latest_eastern_na_{os.getpid()}.tmp"
    if tmp_link.exists() or tmp_link.is_symlink():
        tmp_link.unlink()
    tmp_link.symlink_to(parent_run_id)
    tmp_link.replace(runs_dir / "latest_eastern_na")

def _clean_binomial(raw: str) -> Optional[str]:
    std = standardize_packera_taxon(raw)
    return std if (std and BINOMIAL_RE.match(std)) else None

def _normalize_taxon_inputs(items: Sequence[Any]) -> Tuple[bool, Set[str]]:
    out: Set[str] = set()
    has_nonempty = False
    for item in items:
        for chunk in (item if isinstance(item, (list, tuple)) else ([item] if item is not None else [])):
            for tok in str(chunk).split(","):
                if raw := tok.strip():
                    has_nonempty = True
                    if not HYBRID_RE.search(raw) and (std := _clean_binomial(raw)):
                        out.add(std)
    return has_nonempty, out

def _row_binomial(row: Dict[str, Any]) -> Optional[str]:
    for key in ("species_standardized", "scientificName", "species_raw"):
        raw = "" if row.get(key) is None else str(row[key]).strip()
        if raw.upper() in NA_TOKENS:
            continue
        if HYBRID_RE.search(raw):
            return None
        if std := _clean_binomial(raw):
            return std
    return None

def _norm_pid(val: Any) -> str:
    s = "" if val is None else str(val).strip()
    if s.upper() in NA_TOKENS:
        return "0"
    try:
        return (str(int(round(f))) if f.is_integer() else str(f)) if math.isfinite(f := float(s)) else s
    except ValueError:
        return s

def _count_otus_per_taxon(rows: Sequence[Dict[str, Any]]) -> Dict[str, int]:
    buckets: Dict[str, Set[Tuple[str, str]]] = {}
    for idx, row in enumerate(rows):
        if sp := _row_binomial(row):
            cat = str(row.get("catalogNumber") or "").strip() or f"__row_{idx}"
            buckets.setdefault(sp, set()).add((cat, _norm_pid(row.get("plant_individual_id"))))
    return {sp: len(otus) for sp, otus in buckets.items()}

def _fetch_discovery_run_rows(ws: Path, conn: Optional[sqlite3.Connection], run_id: str) -> List[Dict[str, Any]]:
    if not run_id:
        return []
    if csv_rows := _read_csv_rows(ws / "runs" / run_id / "data" / "tables" / "triage_queue.csv"):
        return csv_rows
    if conn is None:
        return []
    for sql in (
        "SELECT * FROM v_integrative_triage_master WHERE run_id = ?;",
        "SELECT ct.catalogNumber, COALESCE(ct.plant_individual_id, 0) AS plant_individual_id, "
        "v.species_standardized, v.scientificName, v.species_raw "
        "FROM curation_triage ct LEFT JOIN vouchers v ON v.catalogNumber = ct.catalogNumber WHERE ct.run_id = ?;",
    ):
        try:
            if rows := conn.execute(sql, (run_id,)).fetchall():
                return [dict(r) for r in rows]
        except sqlite3.Error:
            continue
    return []

def _query_sqlite_voucher_leaf_rows(conn: sqlite3.Connection) -> List[Dict[str, Any]]:
    try:
        v_cols = {r[1] for r in conn.execute("PRAGMA table_info(vouchers);").fetchall()}
        el_cols = {r[1] for r in conn.execute("PRAGMA table_info(extracted_leaves);").fetchall()}
    except sqlite3.Error:
        return []
    if "catalogNumber" not in v_cols or "catalogNumber" not in el_cols:
        return []
    pid_expr = "COALESCE(el.plant_individual_id, 0) AS plant_individual_id" if "plant_individual_id" in el_cols else "0 AS plant_individual_id"
    opt = [f"v.{c}" for c in ("species_standardized", "scientificName", "species_raw", "stateProvince", "country", "countryCode", "decimalLatitude", "decimalLongitude", "type_status") if c in v_cols]
    try:
        return [dict(r) for r in conn.execute(f"SELECT DISTINCT {', '.join(['v.catalogNumber', pid_expr, *opt])} FROM vouchers v INNER JOIN extracted_leaves el ON el.catalogNumber = v.catalogNumber;").fetchall()]
    except sqlite3.Error:
        return []

def discover_qualifying_eastern_species(
    workspace_root: Path, db_path: Optional[Path] = None, g2_run_id: str = "", g1_run_id: str = "",
    geographic_scope: str = "eastern_na", min_vouchers: int = 15, only_taxa: Sequence[str] = (), exclude_taxa: Sequence[str] = (),
) -> List[Tuple[str, int]]:
    """Discover Packera binomials meeting `min_vouchers` across G2, G1, or SQLite vouchers+leaves."""
    conn = _open_ro_conn(_resolve_ro_db_path(workspace_root, db_path))
    counts: Dict[str, int] = {}
    try:
        for rid in (g2_run_id, g1_run_id):
            if (rows := _fetch_discovery_run_rows(workspace_root, conn, rid.strip())) and (counts := _count_otus_per_taxon(rows)):
                break
        if not counts and conn is not None and (v_rows := _query_sqlite_voucher_leaf_rows(conn)):
            has_geo = any(
                (str(r.get("stateProvince") or "").strip().upper() not in NA_TOKENS or _has_finite_float(r.get("decimalLongitude")))
                for r in v_rows if _row_binomial(r)
            )
            if geographic_scope == "eastern_na" and has_geo:
                import pandas as pd
                keep_mask, _ = evaluate_eastern_na_cohort_filter(pd.DataFrame(v_rows), standardize_packera_taxon)
                v_rows = [r for r, keep in zip(v_rows, keep_mask.tolist()) if keep]
            counts = _count_otus_per_taxon(v_rows)
    finally:
        if conn is not None:
            conn.close()
    has_only, only_set = _normalize_taxon_inputs(only_taxa)
    _, excl_set = _normalize_taxon_inputs(exclude_taxa)
    return sorted((sp, n) for sp, n in counts.items() if BINOMIAL_RE.match(sp) and n >= int(min_vouchers) and (not has_only or sp in only_set) and sp not in excl_set)

def _stage_is_complete(ws: Path, run_id: str) -> bool:
    triage_csv = ws / "runs" / run_id / "data" / "tables" / "triage_queue.csv"
    if not triage_csv.is_file() or triage_csv.stat().st_size == 0:
        return False
    with open(triage_csv, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        return bool(next(reader, None) and next(reader, None))

def _build_genus_commands(run_script: Path, g1_eff: str, g2: str, g3: str, args: argparse.Namespace) -> List[Tuple[str, str, List[str]]]:
    tail = [*(["--db", str(args.db)] if args.db is not None else []), *(["--threads", str(args.threads)] if args.threads and args.threads > 0 else [])]
    gpu = ["--gpu-batch-size", str(args.gpu_batch_size)] if args.gpu_batch_size and args.gpu_batch_size > 0 else []
    cmds: List[Tuple[str, str, List[str]]] = []
    if not args.skip_discovery:
        cmds.append(("Genus_Stage_1_Discovery", g1_eff, ["bash", str(run_script), "--run-id", g1_eff, *(["--skip-harvest"] if args.skip_harvest else []), "--geographic-scope", args.geographic_scope, *tail, *gpu]))
    cmds.append(("Genus_Stage_2_Purified", g2, ["bash", str(run_script), "--run-id", g2, "--skip-harvest", "--geographic-scope", args.geographic_scope, "--exclude-misidentified", "--triage-source", g1_eff, "--gold-anchors", str(args.gold_anchors), *tail, *gpu]))
    cmds.append(("Genus_Stage_3_Gold_Thin", g3, ["bash", str(run_script), "--run-id", g3, "--start-step", "3", "--geographic-scope", args.geographic_scope, "--exclude-misidentified", "--triage-source", g1_eff, "--spatial-thin-km", f"{args.spatial_thin_km:g}", "--gold-only", *tail]))
    return cmds

def _species_commands(run_script: Path, sp: str, sp_prefix: str, g1_eff: str, args: argparse.Namespace) -> List[Tuple[str, str, List[str]]]:
    return build_four_stage_commands(
        run_script=run_script, focal_taxon=sp, prefix=sp_prefix, discovery_run_id=g1_eff, geographic_scope=args.geographic_scope,
        gold_anchors=args.gold_anchors, spatial_thin_km=args.spatial_thin_km, threads=args.threads, gpu_batch_size=args.gpu_batch_size,
        db_path=args.db, skip_discovery=True, skip_harvest=True,
    )

def _compile_genus_summary(ws: Path, genus_dir: Path, g1_eff: str, g2: str, g3: str, conn: Optional[sqlite3.Connection]) -> Tuple[Path, Path]:
    runs_root = ws / "runs"
    rows = [_extract_stage_metrics(stage_name=st, run_id=rid, focal_std="", runs_root=runs_root, conn=conn)[0]
            for st, rid in (("Genus_Stage_1_Discovery", g1_eff), ("Genus_Stage_2_Purified", g2), ("Genus_Stage_3_Gold_Thin", g3))]
    csv_out = genus_dir / "Eastern_NA_Packera_Genus_Three_Stage_Comparison.csv"
    _write_atomic_csv(csv_out, list(rows[0].keys()), rows)
    md_lines = [
        "# Eastern North American *Packera* Genus Baseline Summary (Stages 01 to 03)", "",
        "I ran three genus-wide stages across all qualifying regional sheets before launching per-species sweeps, so misidentified congeners drop out early.", "",
        "| Stage | Run ID | Total OTUs | Taxa Present | Optimal K | Misid. Purged | Hybrids | Min D_edaph |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |",
        *[f"| {r['stage']} | `{r['run_id']}` | {r['n_total_otus']} | {r['n_taxa_present']} | {r['gmm_optimal_k'] or 'N/A'} | "
          f"{r['n_misidentification_reassign']} | {r['n_hybrid_introgressant']} | {r['min_d_edaph'] or 'N/A'} |" for r in rows],
    ]
    md_out = genus_dir / "Eastern_NA_Packera_Genus_Summary.md"
    _write_atomic_text(md_out, "\n".join(md_lines) + "\n")
    return csv_out, md_out

def _lookup_leiden_id(ws: Path, conn: Optional[sqlite3.Connection], sp: str, run_ids: Sequence[str]) -> str:
    runs_root = ws / "runs"
    for rid in filter(None, run_ids):
        for cr in _fetch_complex_rows_for_run(conn, rid, runs_root / rid):
            if standardize_packera_taxon(str(cr.get("species_standardized") or "")) == sp:
                cid = str(cr.get("leiden_complex_id") or cr.get("consensus_complex_id") or "").strip()
                if cid.upper() not in NA_TOKENS:
                    return cid
    return ""

def _build_master_row(
    ws: Path, conn: Optional[sqlite3.Connection], parent_run_id: str, sp: str, n_disc: int, g1_eff: str, g2: str, sweep_ran: bool
) -> Dict[str, Any]:
    slug, runs_root = slugify_taxon(sp), ws / "runs"
    s1, _ = _extract_stage_metrics(stage_name="Stage_1_Discovery", run_id=g1_eff, focal_std=sp, runs_root=runs_root, conn=conn)
    if sweep_ran:
        _, r2, r3, r4 = _stage_run_ids(f"{parent_run_id}/{slug}/", discovery_run_id=g1_eff)
        s2, s3, s4 = (_extract_stage_metrics(stage_name=st, run_id=rid, focal_std=sp, runs_root=runs_root, conn=conn)[0]
                      for st, rid in (("Stage_2_Complex_Triage", r2), ("Stage_3_Intraspecific_Clean", r3), ("Stage_4_Spatially_Thinned", r4)))
    else:
        r2, s2, s3, s4 = g2, dict.fromkeys(s1, ""), dict.fromkeys(s1, ""), {"n_focal_otus": 0, "gmm_optimal_k": ""}
    n_raw, n_misid = int(s1["n_focal_otus"] or 0) or int(n_disc), int(s2.get("n_misidentification_reassign") or 0)
    gdi = next((v for v in (s3.get("median_coalescent_gdi"), s2.get("median_coalescent_gdi"), s1.get("median_coalescent_gdi")) if v is not None), None)
    verdict = ("Genus baseline only (intraspecific sweep skipped)" if not sweep_ran else
               "Candidate cryptic lineage detected" if int(s3.get("n_cohesive_cryptic_otus") or 0) > 0 else
               "Intraspecific ecotypic substructure detected" if int(s3.get("n_subspecies_ecotype") or 0) > 0 else "Single cohesive evolutionary lineage")
    return {
        "parent_run_id": parent_run_id, "species_standardized": sp, "species_subdir": slug,
        "leiden_complex_id": _lookup_leiden_id(ws, conn, sp, (r2, g2, g1_eff)),
        "complex_neighbors": s2.get("complex_neighbors") or s1.get("complex_neighbors") or "",
        "n_raw_otus": n_raw, "n_misid_purged": n_misid, "pct_misid_purged": f"{(100.0 * n_misid / n_raw):.2f}" if n_raw > 0 else "",
        "n_hybrids_stage2": int(s2.get("n_hybrid_introgressant") or 0), "hybrid_parent_pairs": s2.get("hybrid_parent_pairs") or "",
        "n_purified_intraspecific_otus": int(s3.get("n_focal_otus") or 0), "intraspecific_optimal_k": s3.get("gmm_optimal_k", ""),
        "max_two_delta_bic": s3.get("max_two_delta_bic", ""), "n_subspecies_ecotype": int(s3.get("n_subspecies_ecotype") or 0),
        "n_ecophenotypic_plasticity": int(s3.get("n_ecophenotypic_plasticity") or 0),
        "n_cohesive_cryptic_otus": int(s3.get("n_cohesive_cryptic_otus") or 0), "n_confirmed_type_control": int(s3.get("n_confirmed_type_control") or 0),
        "min_d_edaph": s3.get("min_d_edaph", ""), "max_abs_delta_doy_days": s3.get("max_abs_delta_doy_days", ""),
        "median_coalescent_gdi": f"{gdi:.4f}" if gdi is not None else "",
        "public_backbone_clade": s3.get("public_backbone_clade") or s2.get("public_backbone_clade") or s1.get("public_backbone_clade") or "",
        "n_tier0_abstained": int(s3.get("n_tier0_abstained") or 0), "thinned_stage4_otus": int(s4.get("n_focal_otus") or 0),
        "thinned_optimal_k": s4.get("gmm_optimal_k", ""), "primary_verdict": verdict,
    }

def _render_monograph_report(
    parent_run_id: str, args: argparse.Namespace, g1_eff: str, g2: str, g3: str, master_rows: Sequence[Dict[str, Any]]
) -> str:
    lines = [
        "# Eastern North American *Packera* Monograph Synthesis Report", "",
        f"I executed a hierarchical genus-to-species delimitation workflow under parent run `{parent_run_id}` across geographic scope `{args.geographic_scope}`. "
        "Phase A purged interspecific misidentifications across the full genus baseline, whereas Phase B isolated each qualifying taxon inside its own dedicated subdirectory across Stages 02 to 04.", "",
        "## 1. Run Configuration and Genus Stage Identifiers", "",
        f"* **Parent Run ID:** `{parent_run_id}`", f"* **Minimum OTU Threshold:** {args.min_vouchers}",
        f"* **Qualifying Species Evaluated:** {len(master_rows)}", f"* **Genus Stage 1 Discovery:** `{g1_eff}`",
        f"* **Genus Stage 2 Purified:** `{g2}`", f"* **Genus Stage 3 Gold Thin Audit:** `{g3}`", "",
        "## 2. Per-Species Delimitation Matrix", "",
        "| Species | Subdir | Leiden Complex | Raw OTUs | Misid. Purged | Hybrids | Clean OTUs | Opt K | Cohesive Cryptic | Median gdi | Backbone Clade | Thinned OTUs | Verdict |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- | ---: | :--- |",
        *[f"| *{r['species_standardized']}* | `{r['species_subdir']}` | {r['leiden_complex_id'] or 'N/A'} | {r['n_raw_otus']} | "
          f"{r['n_misid_purged']} | {r['n_hybrids_stage2']} | {r['n_purified_intraspecific_otus']} | {r['intraspecific_optimal_k'] or 'N/A'} | "
          f"{r['n_cohesive_cryptic_otus']} | {r['median_coalescent_gdi'] or 'N/A'} | {r['public_backbone_clade'] or 'N/A'} | "
          f"{r['thinned_stage4_otus']} | {r['primary_verdict']} |" for r in master_rows],
        "", "## 3. Lineage Adjudication Notes", "",
    ] + ([f"* *{r['species_standardized']}*: {r['primary_verdict']} (neighbors: {r['complex_neighbors'] or 'None'}; "
          f"confirmed type controls: {r['n_confirmed_type_control']}; Tier 0 abstentions: {r['n_tier0_abstained']})." for r in master_rows]
         if master_rows else ["* No *Packera* binomials exceeded the minimum voucher threshold in this run."])
    return "\n".join(lines) + "\n"

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Hierarchical Eastern North American Packera monograph orchestrator.")
    for flag, dflt in (("--prefix", "eastern_na"), ("--parent-run-id", None), ("--discovery-run-id", None)):
        p.add_argument(flag, default=dflt)
    p.add_argument("--geographic-scope", default="eastern_na", choices=("eastern_na", "north_america", "global"))
    for flag, typ, dflt in (("--min-vouchers", int, 15), ("--gold-anchors", int, 30), ("--spatial-thin-km", float, 15.0), ("--gpu-batch-size", int, None)):
        p.add_argument(flag, type=typ, default=dflt)
    p.add_argument("-j", "--threads", type=int, default=None)
    for flag in ("--db", "--output-dir"):
        p.add_argument(flag, type=Path, default=None)
    p.add_argument("--only-taxa", action="append", nargs="+", default=[])
    p.add_argument("--exclude-taxon", action="append", default=[])
    for flag in ("--resume", "--skip-harvest", "--skip-discovery", "--skip-intraspecific-sweep", "--force", "--summarize-only", "--dry-run"):
        p.add_argument(flag, action="store_true", default=False)
    return p

def _run_stage_batch(ws: Path, prefix_lbl: str, cmds: Sequence[Tuple[str, str, List[str]]], args: argparse.Namespace, statuses: Dict[str, str]) -> int:
    if args.summarize_only:
        statuses.update({rid: "summarize_only" for _, rid, _ in cmds})
        return 0
    for lbl, rid, cmd in cmds:
        full_lbl = f"{prefix_lbl}:{lbl}" if prefix_lbl else lbl
        if not args.force and _stage_is_complete(ws, rid):
            print(f"[Checkpoint Skip] {full_lbl} ({rid}) already has triage_queue.csv.", flush=True)
            statuses[rid] = "skipped_checkpoint"
            continue
        print(f"\n=== Executing {full_lbl} ({rid}) ===", flush=True)
        proc = subprocess.run(cmd, cwd=str(ws), check=False)
        statuses[rid] = "completed" if proc.returncode == 0 else f"failed_{proc.returncode}"
        if proc.returncode != 0:
            print(f"[ERROR] {full_lbl} ({rid}) exited with status {proc.returncode}.", file=sys.stderr)
            return proc.returncode
    return 0

def _git_commit_sha(ws: Path) -> str:
    try:
        res = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(ws), capture_output=True, text=True, timeout=3, check=False)
        return res.stdout.strip() if res.returncode == 0 and res.stdout.strip() else "unknown"
    except Exception:
        return "unknown"

def main(argv: Optional[Sequence[str]] = None, workspace_root: Optional[Path] = None) -> int:
    raw_argv = list(argv) if argv is not None else sys.argv[1:]
    args = build_parser().parse_args(argv)
    ws = (workspace_root or PROJECT_ROOT).resolve()
    parent_run_id, unique_id = resolve_parent_run_id(ws, prefix=args.prefix, parent_run_id=args.parent_run_id, resume=args.resume)
    g1_eff = (args.discovery_run_id or "").strip() or f"{parent_run_id}/00_genus_wide/01_raw_discovery"
    g2, g3 = f"{parent_run_id}/00_genus_wide/02_purified_genus", f"{parent_run_id}/00_genus_wide/03_gold_thin_audit"
    run_script = ws / "run_pipeline.sh"
    genus_cmds = _build_genus_commands(run_script, g1_eff, g2, g3, args)
    disc_kw = dict(db_path=args.db, g2_run_id=g2, g1_run_id=g1_eff, geographic_scope=args.geographic_scope,
                   min_vouchers=args.min_vouchers, only_taxa=args.only_taxa, exclude_taxa=args.exclude_taxon)
    if args.dry_run:
        for lbl, rid, cmd in genus_cmds:
            print(f"[Phase A / {lbl}] ({rid}): {' '.join(cmd)}")
        qualifying = discover_qualifying_eastern_species(ws, **disc_kw)
        if not args.skip_intraspecific_sweep:
            for sp, cnt in qualifying:
                for lbl, rid, cmd in _species_commands(run_script, sp, f"{parent_run_id}/{slugify_taxon(sp)}/", g1_eff, args):
                    print(f"[Phase B / {sp} ({cnt}) / {lbl}] ({rid}): {' '.join(cmd)}")
        print(f"[Dry-Run] Parent={parent_run_id} Qualifying={len(qualifying)} Sweep={not args.skip_intraspecific_sweep}")
        return 0

    parent_dir, genus_dir = ws / "runs" / parent_run_id, ws / "runs" / parent_run_id / "00_genus_wide"
    genus_dir.mkdir(parents=True, exist_ok=True)
    _update_latest_eastern_na_symlink(ws, parent_run_id)
    statuses: Dict[str, str] = {g1_eff: "skipped_discovery_flag"} if args.skip_discovery else {}
    if (rc := _run_stage_batch(ws, "", genus_cmds, args, statuses)) != 0:
        return rc
    conn = _open_ro_conn(_resolve_ro_db_path(ws, args.db))
    try:
        _compile_genus_summary(ws, genus_dir, g1_eff, g2, g3, conn)
        qualifying = discover_qualifying_eastern_species(ws, **disc_kw)
        if not args.skip_intraspecific_sweep:
            for sp, _ in qualifying:
                slug = slugify_taxon(sp)
                sp_prefix = f"{parent_run_id}/{slug}/"
                if (rc := _run_stage_batch(ws, sp, _species_commands(run_script, sp, sp_prefix, g1_eff, args), args, statuses)) != 0:
                    return rc
                compile_complex_resolution_summary(focal_taxon=sp, prefix=sp_prefix, discovery_run_id=g1_eff, workspace_root=ws, db_path=args.db, output_dir=parent_dir / slug)
        master_rows = [_build_master_row(ws, conn, parent_run_id, sp, cnt, g1_eff, g2, not args.skip_intraspecific_sweep) for sp, cnt in qualifying]
    finally:
        if conn is not None:
            conn.close()

    artifacts = (parent_dir / "Eastern_NA_Packera_Master_Synthesis.csv", parent_dir / "Eastern_NA_Packera_Monograph_Report.md", parent_dir / "monograph_run_manifest.json")
    _write_atomic_csv(artifacts[0], MASTER_SYNTHESIS_COLS, master_rows)
    _write_atomic_text(artifacts[1], _render_monograph_report(parent_run_id, args, g1_eff, g2, g3, master_rows))
    _write_atomic_text(artifacts[2], json.dumps({
        "parent_run_id": parent_run_id, "created_date": datetime.now().astimezone().isoformat(), "unique_id": unique_id,
        "git_commit_sha": _git_commit_sha(ws), "cli_arguments": raw_argv,
        "genus_stage_run_ids": {"stage_1_discovery": g1_eff, "stage_2_purified_genus": g2, "stage_3_gold_thin_audit": g3},
        "species_directories": {sp: slugify_taxon(sp) for sp, _ in qualifying}, "stage_completion_statuses": statuses,
    }, indent=2) + "\n")
    if args.output_dir is not None and args.output_dir.resolve() != parent_dir.resolve():
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for src in artifacts:
            shutil.copy2(src, args.output_dir / src.name)
    return 0

if __name__ == "__main__":
    sys.exit(main())
