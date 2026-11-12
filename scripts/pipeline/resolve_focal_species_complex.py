#!/usr/bin/env python3
"""
scripts/pipeline/resolve_focal_species_complex.py
=================================================
Four-stage iterative species-complex discovery and intraspecific delimitation
orchestrator for Packera dubia (and any focal North American Packera taxon).

Executes four isolated, run-scoped passes through ./run_pipeline.sh backed by
shared Tier A SQLite caches (data/packera_pipeline.sqlite), then compiles a
cross-stage quantitative comparison CSV and a structured taxonomic adjudication
report in outputs/reports/:

  Stage 1 (<prefix>_01_eastern_discovery):
      Broad regional harvest, LeafMachine3 organ segmentation, and genus-level
      empirical taxon-complex network discovery.
  Stage 2 (<prefix>_02_complex_triage):
      Targeted focal-taxon boundary audit against auto-discovered high-affinity
      congeneric neighbors to separate Stage 1 Misidentification_Reassign sheets
      from Stage 2 Hybrid_Introgressant and Type_Synonymy_Candidate signals.
  Stage 3 (<prefix>_03_intraspecific_clean):
      Purified intraspecific mode (--intraspecific-mode --exclude-misidentified)
      dedicating 100% of EFA/ECT variance, plasticity ridgeline merging, 0 to 15 cm
      rhizosphere/outcrop contingency, and Hopkins circular phenology to internal
      lineage structure within the focal species sensu stricto.
  Stage 4 (<prefix>_04_gold_thin_audit):
      Fast Step 03+ spatially thinned sensitivity check (--spatial-thin-km) to
      confirm that candidate morphotypes or ecotypes survive collection-density
      de-clumping.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sqlite3
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.core.harvester import standardize_packera_taxon


def slugify_taxon(taxon_name: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9]+", "_", taxon_name.strip()).strip("_")
    return clean or "Packera_focal"


def default_prefix_for_taxon(taxon_name: str) -> str:
    parts = taxon_name.strip().split()
    if len(parts) >= 2 and parts[0].lower() == "packera":
        return re.sub(r"[^a-z0-9]+", "_", parts[1].lower()).strip("_")
    return slugify_taxon(taxon_name).lower()


def build_four_stage_commands(
    *,
    run_script: Path,
    focal_taxon: str,
    prefix: str,
    geographic_scope: str = "eastern_na",
    gold_anchors: int = 30,
    spatial_thin_km: float = 15.0,
    threads: Optional[int] = None,
    gpu_batch_size: Optional[int] = None,
    db_path: Optional[Path] = None,
    skip_discovery: bool = False,
) -> List[Tuple[str, str, List[str]]]:
    """Return ordered (stage_label, run_id, argv) tuples for the four-stage funnel."""
    common_tail: List[str] = []
    if db_path is not None:
        common_tail.extend(["--db", str(db_path)])
    if threads is not None and threads > 0:
        common_tail.extend(["--threads", str(threads)])

    r1 = f"{prefix}_01_eastern_discovery"
    r2 = f"{prefix}_02_complex_triage"
    r3 = f"{prefix}_03_intraspecific_clean"
    r4 = f"{prefix}_04_gold_thin_audit"

    stages: List[Tuple[str, str, List[str]]] = []
    if not skip_discovery:
        cmd1 = [
            "bash",
            str(run_script),
            "--run-id",
            r1,
            "--geographic-scope",
            geographic_scope,
            *common_tail,
        ]
        if gpu_batch_size is not None and gpu_batch_size > 0:
            cmd1.extend(["--gpu-batch-size", str(gpu_batch_size)])
        stages.append(("Stage_1_Discovery", r1, cmd1))

    cmd2 = [
        "bash",
        str(run_script),
        "--run-id",
        r2,
        "--skip-harvest",
        "--focal-taxon",
        focal_taxon,
        "--discover-complexes",
        "--gold-anchors",
        str(gold_anchors),
        *common_tail,
    ]
    stages.append(("Stage_2_Complex_Triage", r2, cmd2))

    cmd3 = [
        "bash",
        str(run_script),
        "--run-id",
        r3,
        "--skip-harvest",
        "--focal-taxon",
        focal_taxon,
        "--intraspecific-mode",
        "--exclude-misidentified",
        "--triage-source",
        r2,
        *common_tail,
    ]
    stages.append(("Stage_3_Intraspecific_Clean", r3, cmd3))

    cmd4 = [
        "bash",
        str(run_script),
        "--run-id",
        r4,
        "--start-step",
        "3",
        "--focal-taxon",
        focal_taxon,
        "--intraspecific-mode",
        "--exclude-misidentified",
        "--triage-source",
        r2,
        "--spatial-thin-km",
        f"{spatial_thin_km:g}",
        *common_tail,
    ]
    stages.append(("Stage_4_Spatially_Thinned", r4, cmd4))
    return stages


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def _safe_float(val: Any) -> Optional[float]:
    if val is None:
        return None
    s = str(val).strip()
    if not s or s.upper() in {"NA", "NAN", "NULL", "NONE"}:
        return None
    try:
        out = float(s)
        return out if out == out else None
    except ValueError:
        return None


def _is_truthy(val: Any) -> bool:
    if val is None:
        return False
    return str(val).strip().upper() in {"1", "TRUE", "T", "YES"}


def _fetch_triage_rows_for_run(
    conn: Optional[sqlite3.Connection],
    run_id: str,
    run_dir: Path,
) -> List[Dict[str, Any]]:
    csv_rows = _read_csv_rows(run_dir / "data" / "tables" / "triage_queue.csv")
    if csv_rows:
        return csv_rows
    if conn is None:
        return []
    try:
        cur = conn.execute(
            "SELECT * FROM v_integrative_triage_master WHERE run_id = ?;",
            (run_id,),
        )
        cols = [d[0] for d in cur.description]
        rows = [dict(zip(cols, r)) for r in cur.fetchall()]
        if rows:
            return rows
    except sqlite3.Error:
        pass
    try:
        cur = conn.execute("SELECT * FROM curation_triage WHERE run_id = ?;", (run_id,))
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]
    except sqlite3.Error:
        return []


def _fetch_complex_rows_for_run(
    conn: Optional[sqlite3.Connection],
    run_id: str,
    run_dir: Path,
) -> List[Dict[str, Any]]:
    csv_rows = _read_csv_rows(run_dir / "outputs" / "reports" / "empirical_taxon_complexes.csv")
    if csv_rows:
        return csv_rows
    if conn is None:
        return []
    try:
        cur = conn.execute("SELECT * FROM taxon_complexes WHERE run_id = ?;", (run_id,))
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, r)) for r in cur.fetchall()]
    except sqlite3.Error:
        return []


def _extract_stage_metrics(
    *,
    stage_name: str,
    run_id: str,
    focal_std: str,
    runs_root: Path,
    conn: Optional[sqlite3.Connection],
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    run_dir = runs_root / run_id
    triage_rows = _fetch_triage_rows_for_run(conn, run_id, run_dir)
    complex_rows = _fetch_complex_rows_for_run(conn, run_id, run_dir)

    gmm_rows = _read_csv_rows(run_dir / "outputs" / "reports" / "gmm_bayes_factors_summary.csv")
    srf_rows = _read_csv_rows(run_dir / "outputs" / "reports" / "multimodal_spatial_rf_summary.csv")
    batch_rows = _read_csv_rows(run_dir / "outputs" / "reports" / "institutional_batch_effect_audit.csv")

    focal_rows = [
        r
        for r in triage_rows
        if (standardize_packera_taxon(str(r.get("species_standardized") or r.get("scientificName") or "")) == focal_std)
    ]
    eval_rows = focal_rows if focal_rows else triage_rows

    def _count_sub(subtype: str) -> int:
        return sum(1 for r in eval_rows if str(r.get("taxonomic_treatment_subtype", "")).strip() == subtype)

    clusters = {
        str(r.get("gmm_cluster", "")).strip()
        for r in eval_rows
        if str(r.get("gmm_cluster", "")).strip() not in {"", "NA", "None"}
    }
    taxa_present = {
        str(r.get("species_standardized", "")).strip()
        for r in triage_rows
        if str(r.get("species_standardized", "")).strip() not in {"", "NA", "None"}
    }

    opt_k: Optional[int] = None
    max_two_delta_bic: Optional[float] = None
    for gr in gmm_rows:
        for k_col in ("optimal_K", "optimal_k", "G_opt", "best_K", "K"):
            kv = _safe_float(gr.get(k_col))
            if kv is not None and opt_k is None:
                opt_k = int(round(kv))
        for b_col in ("two_delta_BIC", "delta_bic_2x", "bayes_factor_2lnB", "delta_BIC"):
            bv = _safe_float(gr.get(b_col))
            if bv is not None:
                max_two_delta_bic = bv if max_two_delta_bic is None else max(max_two_delta_bic, bv)
    if opt_k is None and clusters:
        opt_k = len(clusters)

    d_edaph_vals = [
        v
        for sr in srf_rows
        for v in (_safe_float(sr.get("D_edaph")),)
        if v is not None
    ]
    min_d_edaph = min(d_edaph_vals) if d_edaph_vals else None

    doy_vals = [
        abs(v)
        for r in eval_rows
        for v in (_safe_float(r.get("delta_doy_resid")),)
        if v is not None
    ]
    max_abs_doy = max(doy_vals) if doy_vals else None

    eta_vals = [
        v
        for br in batch_rows
        for v in (_safe_float(br.get("partial_eta_sq")),)
        if v is not None
    ]
    max_eta = max(eta_vals) if eta_vals else None

    cohesive_cryptic = sum(
        1
        for r in eval_rows
        if str(r.get("taxonomic_treatment_subtype", "")).strip() == "Putative_Cryptic_Species"
        and _is_truthy(r.get("cryptic_cluster_cohesion_pass"))
    )

    neighbors: List[str] = []
    focal_cid: Optional[str] = None
    for cr in complex_rows:
        sp = standardize_packera_taxon(str(cr.get("species_standardized") or ""))
        if sp == focal_std:
            focal_cid = str(cr.get("consensus_complex_id") or "").strip()
            nn = standardize_packera_taxon(str(cr.get("nearest_neighbor_taxon") or ""))
            if nn and nn != focal_std:
                neighbors.append(nn)
    if focal_cid:
        for cr in complex_rows:
            sp = standardize_packera_taxon(str(cr.get("species_standardized") or ""))
            cid = str(cr.get("consensus_complex_id") or "").strip()
            if cid == focal_cid and sp and sp != focal_std and sp not in neighbors:
                neighbors.append(sp)

    neigh_txt = run_dir / "outputs" / "reports" / "focal_auto_neighbors.txt"
    if neigh_txt.is_file():
        for ln in neigh_txt.read_text(encoding="utf-8").splitlines():
            sp = standardize_packera_taxon(ln.strip())
            if sp and sp != focal_std and sp not in neighbors:
                neighbors.append(sp)

    metrics = {
        "stage": stage_name,
        "run_id": run_id,
        "n_total_otus": len(triage_rows),
        "n_focal_otus": len(focal_rows),
        "n_taxa_present": len(taxa_present),
        "gmm_optimal_k": opt_k if opt_k is not None else "",
        "max_two_delta_bic": f"{max_two_delta_bic:.2f}" if max_two_delta_bic is not None else "",
        "n_confirmed_species": _count_sub("Confirmed_Species"),
        "n_confirmed_type_control": _count_sub("Confirmed_Type_Control"),
        "n_type_synonymy_candidate": _count_sub("Type_Synonymy_Candidate"),
        "n_misidentification_reassign": _count_sub("Misidentification_Reassign"),
        "n_hybrid_introgressant": _count_sub("Hybrid_Introgressant"),
        "n_subspecies_ecotype": _count_sub("Subspecies_Ecotype"),
        "n_ecophenotypic_plasticity": _count_sub("Ecophenotypic_Plasticity"),
        "n_putative_cryptic_species": _count_sub("Putative_Cryptic_Species"),
        "n_cohesive_cryptic_otus": cohesive_cryptic,
        "min_d_edaph": f"{min_d_edaph:.3f}" if min_d_edaph is not None else "",
        "max_abs_delta_doy_days": f"{max_abs_doy:.1f}" if max_abs_doy is not None else "",
        "max_batch_partial_eta_sq": f"{max_eta:.4f}" if max_eta is not None else "",
        "complex_neighbors": "; ".join(neighbors),
    }
    return metrics, eval_rows


def compile_complex_resolution_summary(
    *,
    focal_taxon: str,
    prefix: str,
    workspace_root: Path = PROJECT_ROOT,
    db_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
) -> Tuple[Path, Path]:
    """Compile cross-run comparison CSV and Markdown synthesis report across Stages 1 to 4."""
    focal_std = standardize_packera_taxon(focal_taxon) or focal_taxon.strip()
    slug = slugify_taxon(focal_std)
    runs_root = workspace_root / "runs"
    out_dir = output_dir or (workspace_root / "outputs" / "reports")
    out_dir.mkdir(parents=True, exist_ok=True)

    resolved_db = db_path or (workspace_root / "data" / "packera_pipeline.sqlite")
    conn: Optional[sqlite3.Connection] = None
    if resolved_db.is_file():
        try:
            conn = sqlite3.connect(f"file:{resolved_db}?mode=ro", uri=True)
        except sqlite3.Error:
            conn = None

    stage_specs = [
        ("Stage_1_Discovery", f"{prefix}_01_eastern_discovery"),
        ("Stage_2_Complex_Triage", f"{prefix}_02_complex_triage"),
        ("Stage_3_Intraspecific_Clean", f"{prefix}_03_intraspecific_clean"),
        ("Stage_4_Spatially_Thinned", f"{prefix}_04_gold_thin_audit"),
    ]

    summary_rows: List[Dict[str, Any]] = []
    stage_triage_map: Dict[str, List[Dict[str, Any]]] = {}
    try:
        for st_name, rid in stage_specs:
            m, eval_rows = _extract_stage_metrics(
                stage_name=st_name,
                run_id=rid,
                focal_std=focal_std,
                runs_root=runs_root,
                conn=conn,
            )
            summary_rows.append(m)
            stage_triage_map[st_name] = eval_rows
    finally:
        if conn is not None:
            conn.close()

    csv_out = out_dir / f"{slug}_Four_Stage_Comparison.csv"
    fieldnames = list(summary_rows[0].keys())
    tmp_csv = csv_out.with_suffix(".csv.tmp")
    with open(tmp_csv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)
        fh.flush()
        os.fsync(fh.fileno())
    tmp_csv.replace(csv_out)

    s1, s2, s3, s4 = summary_rows
    clean_rows = stage_triage_map.get("Stage_3_Intraspecific_Clean", [])

    type_controls = [
        r for r in clean_rows if str(r.get("taxonomic_treatment_subtype", "")).strip() == "Confirmed_Type_Control"
    ]
    priority_pulls = [
        r
        for r in clean_rows
        if str(r.get("taxonomic_treatment_subtype", "")).strip()
        in {"Putative_Cryptic_Species", "Hybrid_Introgressant", "Subspecies_Ecotype", "Type_Synonymy_Candidate"}
        or str(r.get("synthesis_triage_priority", "")).strip() == "CRITICAL"
    ]

    if int(s3["n_cohesive_cryptic_otus"]) > 0:
        verdict = (
            f"Candidate cryptic lineage detected inside *{focal_std}* s.l.: "
            f"{s3['n_cohesive_cryptic_otus']} purified basal-rosette OTUs pass population cohesion "
            f"(optimal K = {s3['gmm_optimal_k'] or 'N/A'}) and remain distinct from the nomenclatural type control."
        )
    elif int(s3["n_subspecies_ecotype"]) > 0:
        verdict = (
            f"Intraspecific ecotypic or geographic substructure detected inside *{focal_std}* "
            f"({s3['n_subspecies_ecotype']} Subspecies_Ecotype OTUs across optimal K = {s3['gmm_optimal_k'] or 1}), "
            f"without cohesive cryptic species barrier closure."
        )
    else:
        verdict = (
            f"*{focal_std}* resolves as a single cohesive evolutionary lineage after purging "
            f"{s2['n_misidentification_reassign']} misidentified congeneric vouchers in Stage 2 "
            f"(purified intraspecific optimal K = {s3['gmm_optimal_k'] or 1}; cohesive cryptic OTUs = 0)."
        )

    thin_note = (
        f"Stage 4 spatial thinning retained {s4['n_focal_otus']} *{focal_std}* OTUs "
        f"with optimal K = {s4['gmm_optimal_k'] or 'N/A'} and {s4['n_cohesive_cryptic_otus']} cohesive cryptic OTUs."
    )

    md_lines = [
        f"# Four-Stage Species Complex Resolution Report: *{focal_std}*",
        "",
        "## 1. Taxonomic Adjudication Verdict",
        "",
        f"* **Primary Delimitation Call:** {verdict}",
        f"* **Complex Neighbor Scope (Stages 1 and 2):** {s2['complex_neighbors'] or s1['complex_neighbors'] or 'None above affinity threshold'}",
        f"* **Interspecific Misidentifications Purged Before Stage 3:** {s2['n_misidentification_reassign']} OTUs",
        f"* **Sympatric Hybrid / Introgressant Candidates (Stage 2):** {s2['n_hybrid_introgressant']} OTUs",
        f"* **Spatial-Thinning Sensitivity Check (Stage 4):** {thin_note}",
        "",
        "## 2. Cross-Stage Quantitative Summary",
        "",
        "| Stage | Run ID | Focal OTUs | Optimal K | 2*Delta BIC | Misid. Purged | Hybrids | Ecotypes | Plasticity | Cohesive Cryptic | Min D_edaph | Max abs(Delta DOY) |",
        "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for r in summary_rows:
        md_lines.append(
            f"| {r['stage']} | `{r['run_id']}` | {r['n_focal_otus']} | {r['gmm_optimal_k'] or 'N/A'} | "
            f"{r['max_two_delta_bic'] or 'N/A'} | {r['n_misidentification_reassign']} | "
            f"{r['n_hybrid_introgressant']} | {r['n_subspecies_ecotype']} | "
            f"{r['n_ecophenotypic_plasticity']} | {r['n_cohesive_cryptic_otus']} | "
            f"{r['min_d_edaph'] or 'N/A'} | {r['max_abs_delta_doy_days'] or 'N/A'} |"
        )

    md_lines.extend([
        "",
        "## 3. Nomenclatural Type Controls (*Sensu Stricto* Anchor)",
        "",
    ])
    if type_controls:
        for tr in type_controls[:10]:
            cat = tr.get("catalogNumber", "Unknown")
            inst = tr.get("institutionCode", "")
            cl = tr.get("gmm_cluster", "N/A")
            md_lines.append(f"* `{cat}` ({inst}): anchored in cluster `{cl}` (`Confirmed_Type_Control`).")
    else:
        md_lines.append("* No explicit nomenclatural type specimen was present in the Stage 3 cohort; verify physical type images before formal nomenclatural assignment.")

    md_lines.extend([
        "",
        "## 4. Priority Herbarium Pull-List for Physical & Genomic Verification",
        "",
    ])
    if priority_pulls:
        for pr in priority_pulls[:25]:
            cat = pr.get("catalogNumber", "Unknown")
            pid = pr.get("plant_individual_id", 0)
            inst = pr.get("institutionCode", "")
            sub = pr.get("taxonomic_treatment_subtype", "")
            rec = pr.get("recommended_determination", "")
            rat = str(pr.get("synthesis_rationale", "")).strip()
            md_lines.append(f"* `{cat}` (individual `{pid}`, {inst}): `{sub}` -> *{rec}*. {rat}")
    else:
        md_lines.append("* Zero CRITICAL or non-confirming vouchers remained in Stage 3 after Stage 2 congener purging.")

    md_out = out_dir / f"{slug}_Complex_Resolution_Report.md"
    tmp_md = md_out.with_suffix(".md.tmp")
    tmp_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    tmp_md.replace(md_out)
    return csv_out, md_out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Execute and synthesize the 4-stage species-complex resolution workflow for Packera dubia or another focal Packera taxon."
    )
    parser.add_argument(
        "--focal-taxon",
        default="Packera dubia",
        help="Focal Packera binomial to resolve (default: 'Packera dubia').",
    )
    parser.add_argument(
        "--prefix",
        default="",
        help="Run ID prefix for the four stages (default: derived from focal epithet, e.g. 'dubia').",
    )
    parser.add_argument(
        "--geographic-scope",
        default="eastern_na",
        choices=("eastern_na", "north_america", "global"),
        help="Geographic scope for Stage 1 baseline discovery (default: eastern_na).",
    )
    parser.add_argument(
        "--gold-anchors",
        type=int,
        default=30,
        help="Target Tier 1 Gold / Silver anchors per neighbor taxon in Stage 2 (default: 30).",
    )
    parser.add_argument(
        "--spatial-thin-km",
        type=float,
        default=15.0,
        help="Spatial thinning distance in km for the Stage 4 sensitivity check (default: 15.0).",
    )
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=None,
        help="CPU worker thread limit forwarded to run_pipeline.sh.",
    )
    parser.add_argument(
        "--gpu-batch-size",
        type=int,
        default=None,
        help="GPU inference batch size for Stage 1 LeafMachine3 segmentation.",
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=None,
        help="Optional custom SQLite WAL database path.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for cross-run comparison CSV and synthesis Markdown report (default: outputs/reports).",
    )
    parser.add_argument(
        "--skip-discovery",
        action="store_true",
        default=False,
        help="Skip Stage 1 when an Eastern NA discovery baseline has already been harvested and segmented.",
    )
    parser.add_argument(
        "--summarize-only",
        action="store_true",
        default=False,
        help="Skip pipeline execution and compile the cross-stage comparison CSV and Markdown report from existing runs.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print the exact four-stage run_pipeline.sh commands without executing them.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    focal_std = standardize_packera_taxon(args.focal_taxon) or args.focal_taxon.strip()
    prefix = args.prefix.strip() or default_prefix_for_taxon(focal_std)
    run_script = PROJECT_ROOT / "run_pipeline.sh"

    stages = build_four_stage_commands(
        run_script=run_script,
        focal_taxon=focal_std,
        prefix=prefix,
        geographic_scope=args.geographic_scope,
        gold_anchors=args.gold_anchors,
        spatial_thin_km=args.spatial_thin_km,
        threads=args.threads,
        gpu_batch_size=args.gpu_batch_size,
        db_path=args.db,
        skip_discovery=args.skip_discovery,
    )

    if args.dry_run:
        for st_label, rid, cmd in stages:
            print(f"[{st_label}] ({rid}): {' '.join(cmd)}")
        return 0

    if not args.summarize_only:
        for st_label, rid, cmd in stages:
            print(f"\n=== Executing {st_label} ({rid}) ===", flush=True)
            proc = subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=False)
            if proc.returncode != 0:
                print(
                    f"[ERROR] {st_label} ({rid}) exited with non-zero status {proc.returncode}.",
                    file=sys.stderr,
                )
                return proc.returncode

    csv_out, md_out = compile_complex_resolution_summary(
        focal_taxon=focal_std,
        prefix=prefix,
        workspace_root=PROJECT_ROOT,
        db_path=args.db,
        output_dir=args.output_dir,
    )
    print(f"[Complex Resolution Summary] Comparison CSV : {csv_out}")
    print(f"[Complex Resolution Summary] Synthesis Report: {md_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
