#!/usr/bin/env python3
"""
===============================================================================
Script: prune_eastern_cohort.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Fast, transactional in-place cohort pruner that restricts an existing
    SQLite database (`data/packera_pipeline.sqlite`), canonical CSV manifests
    (`data/tables/*.csv`), and run tables (`runs/*/data/tables/*.csv`) strictly
    to Cis-Mississippi Eastern North American Packera species while retaining
    trans-Mississippi disjunct populations of predominantly Eastern species
    (e.g., Packera dubia, P. obovata, P. aurea, P. anonyma, P. glabella,
    P. plattensis) up to -100.0W and dropping Cordilleran / Great Plains
    species plus western / Great Plains Packera paupercula.

    Also reconstructs `data/tables/extracted_leaf_manifest.csv`,
    `voucher_reproductive_metrics.csv`, and `voucher_phenological_states.csv`
    from SQLite if Step 02 was interrupted before its final CSV manifest flush.
===============================================================================
"""

from __future__ import annotations

import argparse
import os
import shutil
import sqlite3
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd


def _resolve_project_root(explicit_root: Optional[str] = None) -> Path:
    if explicit_root:
        return Path(explicit_root).resolve()
    here = Path(__file__).resolve()
    for cand in (here.parent, here.parents[1] if len(here.parents) > 1 else here.parent, here.parents[2] if len(here.parents) > 2 else here.parent, Path.cwd().resolve()):
        if (cand / "scripts" / "core" / "config.py").is_file():
            return cand
    return Path.cwd().resolve()


def _write_csv_atomic(df: pd.DataFrame, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(target.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8", newline="") as fh:
        df.to_csv(fh, index=False)
        fh.flush()
        os.fsync(fh.fileno())
    tmp.replace(target)


def _rebuild_step02_csvs_from_sqlite_if_missing(
    conn: sqlite3.Connection,
    tables_dir: Path,
    keep_cats: Set[str],
) -> Dict[str, int]:
    """Reconstruct Step 02 CSV manifests from SQLite if missing or incomplete."""
    rebuilt: Dict[str, int] = {}
    if not keep_cats:
        return rebuilt

    existing_tables = {
        str(r[0]) for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
    }
    leaf_csv = tables_dir / "extracted_leaf_manifest.csv"
    if "extracted_leaves" in existing_tables:
        db_leaf_n = int(conn.execute("SELECT COUNT(*) FROM extracted_leaves;").fetchone()[0])
        csv_leaf_n = 0
        if leaf_csv.is_file() and leaf_csv.stat().st_size > 0:
            try:
                csv_leaf_n = sum(1 for _ in open(leaf_csv, "r", encoding="utf-8")) - 1
            except Exception:
                csv_leaf_n = 0
        if db_leaf_n > 0 and (not leaf_csv.is_file() or csv_leaf_n == 0):
            el_df = pd.read_sql_query("SELECT * FROM extracted_leaves;", conn)
            if "catalogNumber" in el_df.columns:
                el_df = el_df[el_df["catalogNumber"].astype(str).str.strip().isin(keep_cats)].copy()
            if "contour_xy_json" in el_df.columns:
                el_df = el_df.drop(columns=["contour_xy_json"])
            _write_csv_atomic(el_df, leaf_csv)
            rebuilt["extracted_leaf_manifest.csv"] = len(el_df)

    repro_csv = tables_dir / "voucher_reproductive_metrics.csv"
    pheno_csv = tables_dir / "voucher_phenological_states.csv"
    if "voucher_reproductive_phenology" in existing_tables:
        db_rp_n = int(conn.execute("SELECT COUNT(*) FROM voucher_reproductive_phenology;").fetchone()[0])
        if db_rp_n > 0 and (not repro_csv.is_file() or not pheno_csv.is_file()):
            rp_df = pd.read_sql_query("SELECT * FROM voucher_reproductive_phenology;", conn)
            if "catalogNumber" in rp_df.columns:
                rp_df = rp_df[rp_df["catalogNumber"].astype(str).str.strip().isin(keep_cats)].copy()
            if not repro_csv.is_file():
                _write_csv_atomic(rp_df, repro_csv)
                rebuilt["voucher_reproductive_metrics.csv"] = len(rp_df)
            if not pheno_csv.is_file():
                _write_csv_atomic(rp_df, pheno_csv)
                rebuilt["voucher_phenological_states.csv"] = len(rp_df)
    return rebuilt


def prune_workspace_cohort(
    workspace_root: Path,
    db_path: Optional[Path] = None,
    *,
    western_longitude_threshold: float = -100.0,
    transcontinental_longitude_threshold: float = -95.5,
    cis_longitude_threshold: float = -90.5,
    min_cis_records: int = 10,
    min_cis_fraction: float = 0.05,
    strict_whitelist: bool = False,
    delete_dropped_assets: bool = False,
    quarantine_raw: bool = False,
    vacuum: bool = True,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """Prune SQLite database, CSV manifests, and optional disk assets to Cis-Mississippi + Disjunct scope."""
    if str(workspace_root) not in sys.path:
        sys.path.insert(0, str(workspace_root))

    from scripts.core.db import load_taxon_synonyms
    from scripts.core.harvester import _DYNAMIC_SYNONYM_TO_CANONICAL, standardize_packera_taxon
    from scripts.core.taxonomy_authority import evaluate_eastern_na_cohort_filter

    eff_db = (db_path or (workspace_root / "data" / "packera_pipeline.sqlite")).resolve()
    tables_dir = workspace_root / "data" / "tables"
    curated_csv = tables_dir / "curated_vouchers.csv"

    has_db = eff_db.is_file() and eff_db.stat().st_size > 0
    if has_db:
        conn_ro = sqlite3.connect(str(eff_db))
        try:
            try:
                for k_syn, v_can in load_taxon_synonyms(conn_ro).items():
                    _DYNAMIC_SYNONYM_TO_CANONICAL[k_syn] = v_can
            except Exception:
                pass
            v_df = pd.read_sql_query("SELECT * FROM vouchers;", conn_ro)
        finally:
            conn_ro.close()
    elif curated_csv.is_file():
        v_df = pd.read_csv(curated_csv, low_memory=False)
    else:
        raise FileNotFoundError(
            f"Neither SQLite database ({eff_db}) nor curated CSV ({curated_csv}) was found."
        )

    if curated_csv.is_file() and curated_csv.stat().st_size > 0:
        csv_v_df = pd.read_csv(curated_csv, low_memory=False)
        if has_db and not csv_v_df.empty and "catalogNumber" in csv_v_df.columns:
            db_cats = set(v_df["catalogNumber"].astype(str).str.strip()) if not v_df.empty else set()
            extra_csv = csv_v_df[~csv_v_df["catalogNumber"].astype(str).str.strip().isin(db_cats)]
            if not extra_csv.empty:
                v_df = pd.concat([v_df, extra_csv], ignore_index=True)
    else:
        csv_v_df = pd.DataFrame()

    keep_mask, stats = evaluate_eastern_na_cohort_filter(
        v_df,
        standardize_packera_taxon,
        western_longitude_threshold=western_longitude_threshold,
        transcontinental_longitude_threshold=transcontinental_longitude_threshold,
        cis_longitude_threshold=cis_longitude_threshold,
        min_cis_records=min_cis_records,
        min_cis_fraction=min_cis_fraction,
        strict_species_whitelist=strict_whitelist,
    )

    keep_cats: Set[str] = set(
        v_df.loc[keep_mask, "catalogNumber"].dropna().astype(str).str.strip().tolist()
    ) if "catalogNumber" in v_df.columns else set()
    drop_cats: Set[str] = set(
        v_df.loc[~keep_mask, "catalogNumber"].dropna().astype(str).str.strip().tolist()
    ) if "catalogNumber" in v_df.columns else set()

    print("=" * 86)
    print("CIS-MISSISSIPPI EASTERN PACKERA + DISJUNCT COHORT PRUNER")
    print(f"  Workspace Root           : {workspace_root}")
    print(f"  SQLite Database          : {eff_db} ({'FOUND' if has_db else 'MISSING'})")
    print(f"  Disjunct West Ceiling    : {western_longitude_threshold:.2f} deg W (Eastern taxa)")
    print(f"  Transcontinental Ceiling : {transcontinental_longitude_threshold:.2f} deg W (P. paupercula / P. indecora)")
    print(f"  Mode                     : {'DRY-RUN (no files modified)' if dry_run else 'LIVE IN-PLACE PRUNE'}")
    print("-" * 86)
    print(f"{'Taxon':<28} {'Total':>7} {'Cis-MS':>7} {'Disjunct':>9} {'Retained':>9} {'Dropped':>8}  {'Status'}")
    print("-" * 86)
    for sp, info in sorted(stats["species_breakdown"].items()):
        verdict = "KEEP (Eastern)" if info["resident_accepted"] else "DROP (Western/Plains)"
        print(
            f"{sp:<28} {info['total']:>7} {info['cis_mississippi']:>7} "
            f"{info['disjunct_retained']:>9} {info['retained']:>9} {info['dropped']:>8}  {verdict}"
        )
    print("-" * 86)
    print(
        f"SUMMARY: Retained {len(keep_cats):,} / {stats['initial_rows']:,} vouchers "
        f"across {len(stats['retained_species'])} Cis-Mississippi Eastern species "
        f"(purging {len(drop_cats):,} vouchers & {len(stats['dropped_species'])} non-Eastern species)."
    )
    print("=" * 86)

    db_deleted_counts: Dict[str, int] = {}
    csv_updated_counts: Dict[str, Tuple[int, int]] = {}
    removed_files = 0
    quarantined_files = 0
    rebuilt_csvs: Dict[str, int] = {}

    if dry_run:
        return {
            **stats,
            "keep_catalog_count": len(keep_cats),
            "drop_catalog_count": len(drop_cats),
            "db_deleted": db_deleted_counts,
            "csv_updated": csv_updated_counts,
            "removed_files": removed_files,
            "quarantined_files": quarantined_files,
        }

    dropped_asset_paths: Set[Path] = set()
    dropped_raw_paths: Set[Path] = set()
    if (delete_dropped_assets or quarantine_raw) and "image_path" in v_df.columns:
        for rel_p in v_df.loc[~keep_mask, "image_path"].dropna().astype(str):
            rp = Path(rel_p)
            cand = rp if rp.is_absolute() else (workspace_root / rp)
            if cand.is_file():
                dropped_raw_paths.add(cand)

    if has_db and drop_cats:
        conn = sqlite3.connect(str(eff_db))
        try:
            conn.execute("PRAGMA foreign_keys = ON;")
            existing_tables = {
                str(r[0]) for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table';").fetchall()
            }
            if delete_dropped_assets and "extracted_leaves" in existing_tables:
                conn.execute("CREATE TEMP TABLE _tmp_asset_cats (catalogNumber TEXT PRIMARY KEY);")
                conn.executemany(
                    "INSERT OR IGNORE INTO _tmp_asset_cats(catalogNumber) VALUES (?);",
                    [(c,) for c in drop_cats],
                )
                el_cols = {
                    str(r[1]) for r in conn.execute("PRAGMA table_info(extracted_leaves);").fetchall()
                }
                path_cols = [c for c in ("mask_path", "contour_path", "ect_path") if c in el_cols]
                if path_cols:
                    rows_p = conn.execute(
                        f"SELECT {', '.join(path_cols)} FROM extracted_leaves "
                        f"WHERE catalogNumber IN (SELECT catalogNumber FROM _tmp_asset_cats);"
                    ).fetchall()
                    for row_p in rows_p:
                        for p_str in row_p:
                            if p_str and str(p_str).strip():
                                p_obj = Path(str(p_str).strip())
                                full_p = p_obj if p_obj.is_absolute() else (workspace_root / p_obj)
                                if full_p.is_file():
                                    dropped_asset_paths.add(full_p)
                conn.execute("DROP TABLE _tmp_asset_cats;")

            with conn:
                conn.execute("CREATE TEMP TABLE _prune_drop_catalogs (catalogNumber TEXT PRIMARY KEY);")
                conn.executemany(
                    "INSERT OR IGNORE INTO _prune_drop_catalogs(catalogNumber) VALUES (?);",
                    [(c,) for c in drop_cats],
                )
                ordered_cat_tables = (
                    "curation_triage",
                    "spatial_rf_predictions",
                    "otu_morphometrics",
                    "leaf_efa_harmonics",
                    "run_active_cohort",
                    "environmental_cache",
                    "voucher_reproductive_phenology",
                    "extracted_leaves",
                    "plant_individuals",
                    "vouchers",
                )
                for tbl in ordered_cat_tables:
                    if tbl not in existing_tables:
                        continue
                    cur = conn.execute(
                        f"DELETE FROM {tbl} WHERE catalogNumber IN (SELECT catalogNumber FROM _prune_drop_catalogs);"
                    )
                    db_deleted_counts[tbl] = int(cur.rowcount if cur.rowcount is not None and cur.rowcount >= 0 else 0)

                if "taxon_complexes" in existing_tables and stats["retained_species"]:
                    placeholders = ",".join("?" for _ in stats["retained_species"])
                    cur_tc = conn.execute(
                        f"DELETE FROM taxon_complexes WHERE species_standardized NOT IN ({placeholders});",
                        tuple(stats["retained_species"]),
                    )
                    db_deleted_counts["taxon_complexes"] = int(
                        cur_tc.rowcount if cur_tc.rowcount is not None and cur_tc.rowcount >= 0 else 0
                    )
                conn.execute("DROP TABLE _prune_drop_catalogs;")

            rebuilt_csvs = _rebuild_step02_csvs_from_sqlite_if_missing(conn, tables_dir, keep_cats)
            if vacuum:
                conn.execute("VACUUM;")
        finally:
            conn.close()
    elif has_db:
        conn = sqlite3.connect(str(eff_db))
        try:
            rebuilt_csvs = _rebuild_step02_csvs_from_sqlite_if_missing(conn, tables_dir, keep_cats)
        finally:
            conn.close()

    # Filter canonical data/tables/*.csv and any runs/*/data/tables/*.csv
    target_table_dirs: List[Path] = [tables_dir]
    runs_root = workspace_root / "runs"
    if runs_root.is_dir():
        for r_sub in sorted(runs_root.iterdir()):
            if r_sub.is_symlink() or not r_sub.is_dir():
                continue
            rt_dir = r_sub / "data" / "tables"
            if rt_dir.is_dir():
                target_table_dirs.append(rt_dir)

    csv_names = (
        "curated_vouchers.csv",
        "extracted_leaf_manifest.csv",
        "failed_qc_vouchers.csv",
        "voucher_reproductive_metrics.csv",
        "voucher_phenological_states.csv",
        "leaf_efa_harmonics.csv",
        "leaf_efa_harmonics_individual.csv",
        "morphometric_flags.csv",
        "ssurgo_edaphic_validation.csv",
        "phenological_anomalies.csv",
        "multimodal_conflict_flags.csv",
        "triage_queue.csv",
    )
    for t_dir in target_table_dirs:
        for cname in csv_names:
            cpath = t_dir / cname
            if not cpath.is_file() or cpath.stat().st_size == 0:
                continue
            try:
                cdf = pd.read_csv(cpath, low_memory=False)
            except Exception:
                continue
            if "catalogNumber" not in cdf.columns:
                continue
            before_n = len(cdf)
            cdf_kept = cdf[cdf["catalogNumber"].astype(str).str.strip().isin(keep_cats)].copy()
            if len(cdf_kept) != before_n:
                _write_csv_atomic(cdf_kept, cpath)
                rel_key = str(cpath.relative_to(workspace_root))
                csv_updated_counts[rel_key] = (before_n, len(cdf_kept))

    # If canonical curated_vouchers.csv was absent, materialize it from retained SQLite vouchers
    if not curated_csv.is_file() and not v_df.empty:
        kept_v_df = v_df.loc[keep_mask].copy()
        _write_csv_atomic(kept_v_df, curated_csv)
        csv_updated_counts["data/tables/curated_vouchers.csv"] = (len(v_df), len(kept_v_df))

    if quarantine_raw and dropped_raw_paths:
        q_dir = workspace_root / "data" / "raw_vouchers_quarantine"
        q_dir.mkdir(parents=True, exist_ok=True)
        for rp in sorted(dropped_raw_paths):
            if rp.is_file():
                shutil.move(str(rp), str(q_dir / rp.name))
                quarantined_files += 1
    elif delete_dropped_assets:
        for p_del in sorted(dropped_asset_paths | dropped_raw_paths):
            try:
                if p_del.is_file():
                    p_del.unlink()
                    removed_files += 1
            except OSError:
                pass

    if db_deleted_counts:
        print(f"[SQLite Prune] Deleted rows by table: {db_deleted_counts}")
    if rebuilt_csvs:
        print(f"[Manifest Sync] Reconstructed missing Step 02 CSVs from SQLite: {rebuilt_csvs}")
    if csv_updated_counts:
        for k_csv, (b_n, a_n) in csv_updated_counts.items():
            print(f"[CSV Prune] {k_csv}: {b_n:,} -> {a_n:,} rows")
    if quarantined_files:
        print(f"[Quarantine] Moved {quarantined_files:,} non-Eastern raw voucher images to data/raw_vouchers_quarantine/")
    if removed_files:
        print(f"[Disk Cleanup] Removed {removed_files:,} dropped mask/contour/ECT/raw files.")

    return {
        **stats,
        "keep_catalog_count": len(keep_cats),
        "drop_catalog_count": len(drop_cats),
        "db_deleted": db_deleted_counts,
        "csv_updated": csv_updated_counts,
        "rebuilt_csvs": rebuilt_csvs,
        "removed_files": removed_files,
        "quarantined_files": quarantined_files,
    }


def build_cli_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Prune SQLite database and CSV manifests in place to Cis-Mississippi Eastern "
            "Packera species plus trans-Mississippi disjuncts up to -100.0W."
        )
    )
    p.add_argument("--workspace-root", type=str, default=None, help="Repository root directory (default: auto-detected).")
    p.add_argument("--db", type=str, default=None, help="Path to SQLite database (default: data/packera_pipeline.sqlite).")
    p.add_argument(
        "--western-lon",
        type=float,
        default=-100.0,
        help="Western longitude ceiling for trans-Mississippi disjuncts of Eastern species (default: -100.0).",
    )
    p.add_argument(
        "--paupercula-lon",
        type=float,
        default=-95.5,
        help="Western longitude ceiling for transcontinental P. paupercula / P. indecora (default: -95.5).",
    )
    p.add_argument(
        "--cis-lon",
        type=float,
        default=-90.5,
        help="Cis-Mississippi longitude threshold for un-stated eastern coordinates (default: -90.5).",
    )
    p.add_argument(
        "--min-cis-records",
        type=int,
        default=10,
        help="Minimum cis-Mississippi records required for non-whitelisted species discovery (default: 10).",
    )
    p.add_argument(
        "--min-cis-fraction",
        type=float,
        default=0.05,
        help="Minimum fraction of species records east of the Mississippi (default: 0.05).",
    )
    p.add_argument(
        "--strict-whitelist",
        action="store_true",
        default=False,
        help="Restrict strictly to the curated CANONICAL_CIS_MISSISSIPPI_PACKERA_SPECIES checklist.",
    )
    p.add_argument(
        "--delete-dropped-assets",
        action="store_true",
        default=False,
        help="Also delete dropped raw .webp images, masks, contours, and ECT files from disk.",
    )
    p.add_argument(
        "--quarantine-raw",
        action="store_true",
        default=False,
        help="Move dropped raw voucher images to data/raw_vouchers_quarantine/ instead of leaving or deleting them.",
    )
    p.add_argument(
        "--no-vacuum",
        action="store_true",
        default=False,
        help="Skip SQLite VACUUM compaction after row deletion.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        help="Print per-species before/after retention audit without modifying SQLite or CSV files.",
    )
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_cli_parser().parse_args(argv)
    ws = _resolve_project_root(args.workspace_root)
    db_p = Path(args.db).resolve() if args.db else None
    prune_workspace_cohort(
        workspace_root=ws,
        db_path=db_p,
        western_longitude_threshold=float(args.western_lon),
        transcontinental_longitude_threshold=float(args.paupercula_lon),
        cis_longitude_threshold=float(args.cis_lon),
        min_cis_records=int(args.min_cis_records),
        min_cis_fraction=float(args.min_cis_fraction),
        strict_whitelist=bool(args.strict_whitelist),
        delete_dropped_assets=bool(args.delete_dropped_assets),
        quarantine_raw=bool(args.quarantine_raw),
        vacuum=not bool(args.no_vacuum),
        dry_run=bool(args.dry_run),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
