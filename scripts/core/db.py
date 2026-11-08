"""
===============================================================================
Module: db.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Central transactional SQLite (WAL mode) + optional DuckDB relational layer
    (`data/packera_pipeline.sqlite`) implementing the Two-Tier Storage schema:
      - Tier A (Canonical Shared Store):
          taxon_synonyms, vouchers, plant_individuals, extracted_leaves,
          voucher_reproductive_phenology,
          environmental_cache + environmental_rtree (R*Tree spatial index).
      - Tier B (Run-Scoped Analytical Store, partitioned by run_id):
          pipeline_runs, run_active_cohort, leaf_efa_harmonics,
          otu_morphometrics, spatial_rf_predictions,
          curation_triage, and v_integrative_triage_master.

    ABI Invariant:
      Relies strictly on Python's built-in standard library (`sqlite3`, `csv`,
      `json`, `hashlib`, etc.) at import time so `.venv_LM3` requires zero new
      binary package installations. `duckdb` is imported lazily only when
      explicitly requested inside `.venv`.
===============================================================================
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import logging
import math
import os
from pathlib import Path
import re
import sqlite3
import tempfile
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

logger = logging.getLogger("PackeraDB")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SQLITE_REL_PATH = Path("data") / "packera_pipeline.sqlite"

DWC_VOUCHER_COLUMNS: Tuple[str, ...] = (
    "catalogNumber",
    "institutionCode",
    "scientificName",
    "species_raw",
    "recordedBy",
    "recordNumber",
    "identifiedBy",
    "determiner_raw",
    "determiner_tier",
    "type_status",
    "county",
    "stateProvince",
    "decimalLatitude",
    "latitude",
    "decimalLongitude",
    "longitude",
    "coordinateUncertainty",
    "year",
    "month",
    "day",
    "eventDate",
    "regional_group",
    "image_path",
    "exsiccatae_key",
    "is_primary_duplicate",
    "duplicate_count",
)

HARMONIC_COLS: Tuple[str, ...] = tuple(
    f"{coef}{idx}" for idx in range(1, 25) for coef in ("A", "B", "C", "D")
)
PC_COLS: Tuple[str, ...] = tuple(f"PC{idx}" for idx in range(1, 16))
SYM_PC_COLS: Tuple[str, ...] = tuple(f"Sym_PC{idx}" for idx in range(1, 16))

WC_BIO_COLUMNS: Tuple[str, ...] = (
    "bio1_temp_mean",
    "bio2_mean_diurnal_range",
    "bio3_isothermality",
    "bio4_temp_seasonality",
    "bio5_max_temp_warmest_month",
    "bio6_min_temp_coldest_month",
    "bio7_temp_annual_range",
    "bio8_mean_temp_wettest_qtr",
    "bio9_mean_temp_driest_qtr",
    "bio10_mean_temp_warmest_qtr",
    "bio11_mean_temp_coldest_qtr",
    "bio12_precip_annual",
    "bio13_precip_wettest_month",
    "bio14_precip_driest_month",
    "bio15_precip_seasonality",
    "bio16_precip_wettest_qtr",
    "bio17_precip_driest_qtr",
    "bio18_precip_warmest_qtr",
    "bio19_precip_coldest_qtr",
)

ENVIRONMENTAL_CACHE_COLUMNS: Tuple[str, ...] = (
    "lat_lon_key",
    "catalogNumber",
    "latitude",
    "longitude",
    "musym",
    "muname",
    "taxsuborder",
    "taxgreatgroup",
    "is_rock_outcrop",
    "outcrop_class",
    "ssurgo_rock_outcrop_pct",
    "ssurgo_bedrock_depth_cm",
    "ssurgo_sand_pct",
    "ssurgo_clay_pct",
    "ssurgo_ph",
    "ssurgo_om_pct",
    "ssurgo_cec",
    "ssurgo_aws_0_50cm",
    "polaris_ph_30m",
    "polaris_sand_30m",
    "polaris_clay_30m",
    "polaris_om_30m",
    "polaris_bd_30m",
    "polaris_ksat_30m",
    "solus_resdept_cm",
    "solus_fragvol_pct",
    "topo_tpi_30m",
    "topo_roughness_30m",
    "topo_slope_deg",
    "topo_twi_30m",
    "topo_elev_m",
    "fine_soil_ph",
    "fine_soil_sand",
    "fine_soil_clay",
    "fine_bedrock_depth_cm",
    "fine_rock_outcrop_prob",
    "bdricm",
    "resdept_r",
    "cfvo",
    "awc",
    "fine_soil_cec",
    "soil_silt_pct",
    "alr_sand_clay",
    "alr_silt_clay",
    "clr_sand_silt",
    "clr_clay_geom",
    "sgmc_lithology_class",
    "is_azonal_lithology",
    "azonal_contact_dist_m",
    "azonal_proximity_index",
    "buffer_radius_m",
    "high_coord_uncertainty_flag",
    "is_training_eligible_coord",
    "ph_q10_buffer",
    "ph_q90_buffer",
    "min_bedrock_depth_buffer_cm",
    "max_cfvo_buffer_pct",
    "azonal_lithology_presence_buffer",
    *WC_BIO_COLUMNS,
    "bio18_precip_warmest_qt",
    "bio19_precip_coldest_qt",
    "soil_ph",
    "soil_cec",
    "soil_sand",
    "soil_bulk_density",
    "soil_ph_soilgrids",
    "soil_cec_soilgrids",
    "soil_sand_soilgrids",
    "soil_bulk_density_soilgrids",
)

ENVIRONMENTAL_TEXT_COLS: Set[str] = {
    "lat_lon_key",
    "catalogNumber",
    "musym",
    "muname",
    "taxsuborder",
    "taxgreatgroup",
    "outcrop_class",
    "sgmc_lithology_class",
}
ENVIRONMENTAL_BOOL_COLS: Set[str] = {
    "is_rock_outcrop",
    "is_azonal_lithology",
    "high_coord_uncertainty_flag",
    "is_training_eligible_coord",
    "azonal_lithology_presence_buffer",
}


# =============================================================================
# Scalar Coercion & Path Resolution Helpers
# =============================================================================

def _is_missing(val: Any) -> bool:
    if val is None:
        return True
    if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
        return True
    if isinstance(val, str) and val.strip().lower() in {"", "na", "nan", "null", "none"}:
        return True
    return False


def _to_str(val: Any) -> Optional[str]:
    if _is_missing(val):
        return None
    s = str(val).strip()
    return s if s else None


def _to_float(val: Any) -> Optional[float]:
    if _is_missing(val):
        return None
    if isinstance(val, bool):
        return 1.0 if val else 0.0
    try:
        f = float(val)
        return None if (math.isnan(f) or math.isinf(f)) else f
    except (TypeError, ValueError):
        return None


def _to_int(val: Any, default: Optional[int] = None) -> Optional[int]:
    if _is_missing(val):
        return default
    if isinstance(val, bool):
        return 1 if val else 0
    try:
        return int(round(float(val)))
    except (TypeError, ValueError):
        return default


def _to_bool_int(val: Any, default: Optional[int] = None) -> Optional[int]:
    if _is_missing(val):
        return default
    if isinstance(val, bool):
        return 1 if val else 0
    s = str(val).strip().upper()
    if s in {"TRUE", "T", "1", "1.0", "YES", "Y"}:
        return 1
    if s in {"FALSE", "F", "0", "0.0", "NO", "N"}:
        return 0
    return default


def _normalize_tier(raw_tier: Any) -> str:
    s = str(raw_tier or "").strip().lower()
    if "2" in s or "reflect" in s:
        return "Tier_2_Reflected"
    return "Tier_1_Direct"


def format_lat_lon_key(lat: Any, lon: Any) -> Optional[str]:
    f_lat = _to_float(lat)
    f_lon = _to_float(lon)
    if f_lat is None or f_lon is None:
        return None
    return f"{f_lat:.5f}_{f_lon:.5f}"


def compute_file_sha256(path: Union[str, Path]) -> Optional[str]:
    p = Path(path)
    if not p.is_file():
        return None
    h = hashlib.sha256()
    try:
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError:
        return None


def resolve_db_path(
    db_path: Optional[Union[str, Path]] = None,
    workspace_root: Path = PROJECT_ROOT,
    demo: bool = False,
    run_dir: Optional[Union[str, Path]] = None,
) -> Union[str, Path]:
    """
    Resolve the SQLite database path while strictly isolating `--demo` runs from
    the canonical `data/packera_pipeline.sqlite` database.
    """
    canonical = (Path(workspace_root).resolve() / DEFAULT_SQLITE_REL_PATH).resolve()

    if demo:
        if db_path is not None and str(db_path).strip() == ":memory:":
            return ":memory:"
        if db_path is not None and str(db_path).strip():
            cand = Path(str(db_path).strip()).resolve()
            if cand != canonical:
                return cand
        eff_run = str(run_dir or os.environ.get("PACKERA_RUN_DIR", "")).strip()
        if eff_run:
            rdir = Path(eff_run).resolve()
            if rdir.exists():
                return rdir / "demo_ephemeral.sqlite"
        return ":memory:"

    if db_path is not None and str(db_path).strip():
        raw = str(db_path).strip()
        if raw == ":memory:":
            return ":memory:"
        p = Path(raw)
        return p.resolve() if p.is_absolute() else (Path(workspace_root).resolve() / p).resolve()

    env_db = os.environ.get("PACKERA_DB_PATH", "").strip()
    if env_db:
        if env_db == ":memory:":
            return ":memory:"
        p = Path(env_db)
        return p.resolve() if p.is_absolute() else (Path(workspace_root).resolve() / p).resolve()

    return canonical


def is_canonical_or_explicit_db_context(
    target_path: Optional[Union[str, Path]] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> bool:
    """
    Return True when `PACKERA_DB_PATH` is explicitly configured OR when `target_path`
    lives inside `workspace_root` (either canonical `data/` or `runs/<run_id>/`),
    and False when a standalone unit test writes into an isolated `/tmp` directory
    without setting `PACKERA_DB_PATH`.
    """
    if os.environ.get("PACKERA_DB_PATH", "").strip():
        return True
    if target_path is None:
        return True
    try:
        tp = Path(target_path).resolve()
        wr = Path(workspace_root).resolve()
        return tp == wr or wr in tp.parents
    except OSError:
        return False


def _apply_sqlite_pragmas(conn: sqlite3.Connection, read_only: bool = False) -> None:
    """Apply high-throughput workstation SQLite PRAGMAs (256 MiB page cache, 2 GiB mmap, RAM temp store)."""
    conn.execute("PRAGMA busy_timeout = 30000;")
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.execute("PRAGMA cache_size = -262144;")
    conn.execute("PRAGMA temp_store = MEMORY;")
    try:
        conn.execute("PRAGMA mmap_size = 2147483648;")
    except sqlite3.Error:
        pass
    if not read_only:
        conn.execute("PRAGMA journal_mode = WAL;")
        conn.execute("PRAGMA synchronous = NORMAL;")
        try:
            conn.execute("PRAGMA wal_autocheckpoint = 4000;")
        except sqlite3.Error:
            pass


def get_db_connection(
    db_path: Optional[Union[str, Path]] = None,
    read_only: bool = False,
    demo: bool = False,
    init_schema: bool = False,
    workspace_root: Path = PROJECT_ROOT,
    check_same_thread: bool = False,
) -> sqlite3.Connection:
    """
    Open a SQLite connection configured with WAL mode, NORMAL synchronous durability,
    strict foreign-key enforcement, 256 MiB page cache, 2 GiB mmap, and a 30-second busy timeout.
    """
    resolved = resolve_db_path(db_path=db_path, workspace_root=workspace_root, demo=demo)

    if read_only and resolved != ":memory:":
        p = Path(resolved)
        if not p.exists():
            raise FileNotFoundError(f"SQLite database not found for read-only connection: {p}")
        conn = sqlite3.connect(
            f"file:{p.as_posix()}?mode=ro",
            uri=True,
            timeout=30.0,
            check_same_thread=check_same_thread,
        )
        conn.row_factory = sqlite3.Row
        _apply_sqlite_pragmas(conn, read_only=True)
        return conn

    if resolved != ":memory:":
        Path(resolved).parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(resolved), timeout=30.0, check_same_thread=check_same_thread)
    else:
        conn = sqlite3.connect(":memory:", timeout=30.0, check_same_thread=check_same_thread)

    conn.row_factory = sqlite3.Row
    _apply_sqlite_pragmas(conn, read_only=False)

    if init_schema:
        init_db_schema(conn)
    return conn


def get_duckdb_connection(
    db_path: Optional[Union[str, Path]] = None,
    read_only: bool = True,
    workspace_root: Path = PROJECT_ROOT,
) -> Any:
    """
    Lazily attach the SQLite database inside an in-memory DuckDB analytical session
    when `duckdb` is installed in `.venv`.
    """
    import duckdb  # Lazy import so .venv_LM3 never requires duckdb

    resolved = resolve_db_path(db_path=db_path, workspace_root=workspace_root, demo=False)
    dconn = duckdb.connect(database=":memory:")
    n_threads = max(1, min(32, os.cpu_count() or 8))
    try:
        dconn.execute(f"PRAGMA threads={n_threads};")
    except Exception:
        pass
    try:
        dconn.execute("LOAD sqlite;")
    except Exception:
        dconn.execute("INSTALL sqlite; LOAD sqlite;")
    ro_clause = " (TYPE SQLITE, READ_ONLY)" if read_only else " (TYPE SQLITE)"
    escaped_path = str(resolved).replace("'", "''")
    dconn.execute(f"ATTACH '{escaped_path}' AS packera{ro_clause};")
    dconn.execute("USE packera;")
    return dconn


# =============================================================================
# Two-Tier Relational Schema DDL (Tier A + Tier B + R*Tree + Master View)
# =============================================================================

from scripts.core.db_schema import (
    MASTER_VIEW_DDL,
    PACKERA_SCHEMA_USER_VERSION,
    SCHEMA_DDL,
    _DYNAMIC_TABLE_MIGRATION_SPECS,
    _MASTER_VIEW_REQUIRED_TABLES,
    _VOUCHERS_MIGRATION_COLUMNS,
    _build_master_view_ddl,
    _build_schema_ddl,
    _migrate_dynamic_schema_columns,
    _migrate_vouchers_schema_columns,
    init_db_schema,
)


# =============================================================================
# Atomic CSV Export & Generic Upsert Utilities
# =============================================================================

def export_table_or_view_to_csv_atomic(
    conn: sqlite3.Connection,
    query: str,
    params: Sequence[Any],
    target_csv: Union[str, Path],
) -> Path:
    """
    Execute `query` with `params` and atomically write the result set to `target_csv`
    via `.tmp` + `flush()` + `os.fsync()` + `Path.replace()`.
    """
    target_path = Path(target_csv)
    target_path.parent.mkdir(parents=True, exist_ok=True)

    cur = conn.execute(query, tuple(params))
    headers = [desc[0] for desc in (cur.description or [])]
    rows = cur.fetchall()

    tmp_fd, tmp_name = tempfile.mkstemp(
        dir=str(target_path.parent),
        prefix=f".{target_path.stem}_",
        suffix=".csv.tmp",
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(tmp_fd, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            if headers:
                writer.writerow(headers)
            for row in rows:
                writer.writerow(["" if val is None else val for val in row])
            f.flush()
            os.fsync(f.fileno())
        tmp_path.replace(target_path)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)
        raise

    return target_path


def get_table_columns(conn: sqlite3.Connection, table_name: str) -> List[str]:
    cur = conn.execute(f"PRAGMA table_info({table_name});")
    return [str(r[1]) for r in cur.fetchall()]


def upsert_rows_generic(
    conn: sqlite3.Connection,
    table_name: str,
    rows: Sequence[Mapping[str, Any]],
    pk_cols: Sequence[str],
    preserve_existing_on_null: bool = True,
) -> int:
    """
    Parameterized batch UPSERT (`ON CONFLICT(...) DO UPDATE`) for arbitrary rows.
    """
    if not rows:
        return 0
    valid_cols = get_table_columns(conn, table_name)
    valid_set = set(valid_cols)

    present_cols: List[str] = []
    seen: Set[str] = set()
    for pk in pk_cols:
        if pk in valid_set and pk not in seen:
            present_cols.append(pk)
            seen.add(pk)
    for r in rows:
        for k in r.keys():
            if k in valid_set and k not in seen:
                present_cols.append(k)
                seen.add(k)

    non_pk_cols = [c for c in present_cols if c not in pk_cols]
    col_list_sql = ", ".join(present_cols)
    placeholders_sql = ", ".join("?" for _ in present_cols)
    pk_list_sql = ", ".join(pk_cols)

    if non_pk_cols:
        if preserve_existing_on_null:
            updates_sql = ", ".join(
                f"{c} = COALESCE(excluded.{c}, {table_name}.{c})" for c in non_pk_cols
            )
        else:
            updates_sql = ", ".join(f"{c} = excluded.{c}" for c in non_pk_cols)
        sql = (
            f"INSERT INTO {table_name} ({col_list_sql}) VALUES ({placeholders_sql}) "
            f"ON CONFLICT ({pk_list_sql}) DO UPDATE SET {updates_sql};"
        )
    else:
        sql = f"INSERT OR IGNORE INTO {table_name} ({col_list_sql}) VALUES ({placeholders_sql});"

    payload: List[Tuple[Any, ...]] = []
    for r in rows:
        vals: List[Any] = []
        for c in present_cols:
            v = r.get(c)
            if _is_missing(v):
                vals.append(None)
            elif isinstance(v, bool):
                vals.append(1 if v else 0)
            else:
                vals.append(v)
        payload.append(tuple(vals))

    with conn:
        conn.executemany(sql, payload)
    return len(payload)


# =============================================================================
# Step 01: Taxon Synonyms & Vouchers Persistence
# =============================================================================

def upsert_taxon_synonyms(
    conn: sqlite3.Connection,
    synonym_to_canonical: Mapping[str, str],
    author_strings_by_canonical: Optional[Mapping[str, Iterable[str]]] = None,
) -> int:
    """Persist discovered GBIF / configuration taxon synonym mappings into `taxon_synonyms`."""
    if not synonym_to_canonical:
        return 0
    now_iso = datetime.now(timezone.utc).isoformat()
    authors_map: Dict[str, str] = {}
    if author_strings_by_canonical:
        for canon, auths in author_strings_by_canonical.items():
            clean_list = sorted({str(a).strip() for a in (auths or []) if str(a).strip()})
            if clean_list:
                authors_map[str(canon).strip()] = json.dumps(clean_list, ensure_ascii=False)

    payload: List[Tuple[str, str, Optional[str], str]] = []
    for raw_low, canon in synonym_to_canonical.items():
        r_key = re.sub(r"\s+", " ", str(raw_low or "").strip().lower())
        c_val = str(canon or "").strip()
        if not r_key or not c_val:
            continue
        payload.append((r_key, c_val, authors_map.get(c_val), now_iso))

    if not payload:
        return 0

    sql = """
    INSERT INTO taxon_synonyms (raw_name_lower, canonical_packera, protologue_authors_json, updated_at)
    VALUES (?, ?, ?, ?)
    ON CONFLICT(raw_name_lower) DO UPDATE SET
        canonical_packera = excluded.canonical_packera,
        protologue_authors_json = COALESCE(excluded.protologue_authors_json, taxon_synonyms.protologue_authors_json),
        updated_at = excluded.updated_at;
    """
    with conn:
        conn.executemany(sql, payload)
    return len(payload)


def load_taxon_synonyms(conn: sqlite3.Connection) -> Dict[str, str]:
    """Load persisted `raw_name_lower -> canonical_packera` dictionary from `taxon_synonyms`."""
    cur = conn.execute("SELECT raw_name_lower, canonical_packera FROM taxon_synonyms;")
    return {str(r[0]): str(r[1]) for r in cur.fetchall() if r[0] and r[1]}


def _infer_species_standardized_fallback(sci_name: Optional[str], sp_raw: Optional[str], syn_map: Mapping[str, str]) -> Optional[str]:
    for cand in (sci_name, sp_raw):
        if not cand:
            continue
        low = re.sub(r"\s+", " ", cand.strip().lower())
        if low in syn_map:
            return syn_map[low]
        m = re.match(r"^([A-Za-z]+)(?:\s*×\s*|\s+[xX×]\s+|\s+)([a-zA-Z][a-zA-Z-]*)", cand.strip())
        if m:
            genus = m.group(1).capitalize()
            ep = m.group(2).lower().strip("-")
            binom_low = f"{genus.lower()} {ep}"
            if binom_low in syn_map:
                return syn_map[binom_low]
            if genus == "Packera" and len(ep) >= 2:
                return f"Packera {ep}"
    return None


def upsert_vouchers_records(
    conn: sqlite3.Connection,
    records: Sequence[Mapping[str, Any]],
) -> int:
    """
    UPSERT Darwin Core voucher rows into `vouchers`, preserving existing Step 02
    calibration and QC columns when not explicitly provided in `records`.
    """
    if not records:
        return 0
    syn_map = load_taxon_synonyms(conn)
    norm_rows: List[Dict[str, Any]] = []

    for rec in records:
        cat = _to_str(rec.get("catalogNumber"))
        if not cat:
            continue
        sci = _to_str(rec.get("scientificName")) or _to_str(rec.get("species_raw"))
        sp_raw = _to_str(rec.get("species_raw")) or sci
        lat = _to_float(rec.get("latitude"))
        dec_lat = _to_float(rec.get("decimalLatitude"))
        if lat is None:
            lat = dec_lat
        if dec_lat is None:
            dec_lat = lat

        lon = _to_float(rec.get("longitude"))
        dec_lon = _to_float(rec.get("decimalLongitude"))
        if lon is None:
            lon = dec_lon
        if dec_lon is None:
            dec_lon = lon

        ident = _to_str(rec.get("identifiedBy")) or _to_str(rec.get("determiner_raw"))
        det_raw = _to_str(rec.get("determiner_raw")) or ident
        sp_std = _to_str(rec.get("species_standardized")) or _infer_species_standardized_fallback(sci, sp_raw, syn_map)

        row_dict: Dict[str, Any] = {
            "catalogNumber": cat,
            "institutionCode": _to_str(rec.get("institutionCode")),
            "scientificName": sci,
            "species_raw": sp_raw,
            "recordedBy": _to_str(rec.get("recordedBy")),
            "recordNumber": _to_str(rec.get("recordNumber")),
            "identifiedBy": ident,
            "determiner_raw": det_raw,
            "determiner_tier": _to_str(rec.get("determiner_tier")),
            "type_status": _to_str(rec.get("type_status")),
            "county": _to_str(rec.get("county")),
            "stateProvince": _to_str(rec.get("stateProvince")),
            "locality": _to_str(rec.get("locality")),
            "decimalLatitude": dec_lat,
            "latitude": lat,
            "decimalLongitude": dec_lon,
            "longitude": lon,
            "coordinateUncertainty": _to_float(rec.get("coordinateUncertainty")),
            "year": _to_int(rec.get("year")),
            "month": _to_int(rec.get("month")),
            "day": _to_int(rec.get("day")),
            "eventDate": _to_str(rec.get("eventDate")),
            "doy": _to_float(rec.get("doy")),
            "pheno_sin": _to_float(rec.get("pheno_sin")),
            "pheno_cos": _to_float(rec.get("pheno_cos")),
            "regional_group": _to_str(rec.get("regional_group")),
            "image_path": _to_str(rec.get("image_path")),
            "exsiccatae_key": _to_str(rec.get("exsiccatae_key")),
            "species_standardized": sp_std,
            "pixels_per_mm": _to_float(rec.get("pixels_per_mm")),
            "is_calibrated_mm": _to_bool_int(rec.get("is_calibrated_mm")),
        }
        if "is_primary_duplicate" in rec and not _is_missing(rec.get("is_primary_duplicate")):
            row_dict["is_primary_duplicate"] = _to_bool_int(rec.get("is_primary_duplicate"), default=1)
        if "duplicate_count" in rec and not _is_missing(rec.get("duplicate_count")):
            row_dict["duplicate_count"] = _to_int(rec.get("duplicate_count"), default=1)
        if "exsiccatae_taxonomic_conflict" in rec and not _is_missing(rec.get("exsiccatae_taxonomic_conflict")):
            row_dict["exsiccatae_taxonomic_conflict"] = _to_bool_int(rec.get("exsiccatae_taxonomic_conflict"), default=0)
        if "exsiccatae_gold_taxon" in rec and not _is_missing(rec.get("exsiccatae_gold_taxon")):
            row_dict["exsiccatae_gold_taxon"] = _to_str(rec.get("exsiccatae_gold_taxon"))
        qc_val = _to_str(rec.get("qc_status"))
        if qc_val in {"PASS", "FAILED_QC", "PENDING"}:
            row_dict["qc_status"] = qc_val
        if "failure_reason" in rec and not _is_missing(rec.get("failure_reason")):
            row_dict["failure_reason"] = _to_str(rec.get("failure_reason"))
        if "failure_details" in rec and not _is_missing(rec.get("failure_details")):
            row_dict["failure_details"] = _to_str(rec.get("failure_details"))
        norm_rows.append(row_dict)

    return upsert_rows_generic(
        conn,
        table_name="vouchers",
        rows=norm_rows,
        pk_cols=("catalogNumber",),
        preserve_existing_on_null=True,
    )


_CSV_BOOL_VOUCHER_COLUMNS = frozenset({"is_primary_duplicate", "exsiccatae_taxonomic_conflict", "is_calibrated_mm"})


def _format_voucher_csv_value(col: str, val: Any) -> Any:
    """Render a SQLite `vouchers` value in the pandas-style `curated_vouchers.csv` convention."""
    if _is_missing(val):
        return ""
    if col in _CSV_BOOL_VOUCHER_COLUMNS and not isinstance(val, str):
        b = _to_bool_int(val)
        return "" if b is None else ("True" if b else "False")
    return val


def reconcile_curated_vouchers_csv_with_db(
    conn: sqlite3.Connection,
    csv_path: Union[str, Path],
    taxa: Optional[Iterable[str]] = None,
    taxon_resolver: Optional[Any] = None,
) -> int:
    """
    Append canonical SQLite `vouchers` rows that are missing from a
    `curated_vouchers.csv` snapshot (matched on `catalogNumber`).

    Taxon-restricted harvests (`--taxon` / `--focal-taxon`) UPSERT into `vouchers`
    but are intentionally never promoted to the root CSV snapshot, so the snapshot
    can silently lack requested taxa. When `taxa` is given, only vouchers whose
    resolved taxon is in `taxa` are restored; `taxon_resolver(scientificName,
    species_raw, species_standardized)` resolves each row (defaulting to the
    stored `species_standardized`).

    Existing CSV rows are preserved byte-for-byte; appended rows are projected onto
    the existing CSV header (or the DwC voucher schema when the CSV does not exist)
    and written atomically (`.tmp` + `fsync` + `replace`). Returns the number of
    appended rows.
    """
    target = Path(csv_path)
    db_cols = get_table_columns(conn, "vouchers")
    if "catalogNumber" not in db_cols:
        return 0

    taxa_set = {str(t).strip() for t in (taxa or []) if t is not None and str(t).strip()}
    restrict_taxa = taxa is not None
    if restrict_taxa and not taxa_set:
        return 0

    original_text = ""
    header: List[str] = []
    existing_cats: Set[str] = set()
    if target.is_file() and target.stat().st_size > 0:
        with open(target, "r", encoding="utf-8", newline="") as fh:
            original_text = fh.read()
        reader = csv.DictReader(io.StringIO(original_text, newline=""))
        header = list(reader.fieldnames or [])
        if header and "catalogNumber" not in header:
            raise ValueError(f"{target} is missing the required 'catalogNumber' column.")
        for row in reader:
            cat = str(row.get("catalogNumber") or "").strip()
            if cat:
                existing_cats.add(cat)
    if not header:
        original_text = ""
        header = [c for c in DWC_VOUCHER_COLUMNS if c in db_cols]

    taxon_cols = [c for c in ("scientificName", "species_raw", "species_standardized") if c in db_cols]
    select_cols: List[str] = []
    for c in ["catalogNumber"] + [h for h in header if h in db_cols] + taxon_cols:
        if c not in select_cols:
            select_cols.append(c)

    cur = conn.execute(f"SELECT {', '.join(select_cols)} FROM vouchers ORDER BY catalogNumber;")
    new_rows: List[List[Any]] = []
    for rec in cur.fetchall():
        rec_map = dict(zip(select_cols, tuple(rec)))
        cat = str(rec_map.get("catalogNumber") or "").strip()
        if not cat or cat in existing_cats:
            continue
        if restrict_taxa:
            if taxon_resolver is not None:
                resolved = taxon_resolver(
                    rec_map.get("scientificName"),
                    rec_map.get("species_raw"),
                    rec_map.get("species_standardized"),
                )
            else:
                resolved = _to_str(rec_map.get("species_standardized"))
            if not resolved or str(resolved).strip() not in taxa_set:
                continue
        existing_cats.add(cat)
        new_rows.append([_format_voucher_csv_value(c, rec_map.get(c)) for c in header])

    if not new_rows:
        return 0

    target.parent.mkdir(parents=True, exist_ok=True)
    tmp_fd, tmp_name = tempfile.mkstemp(dir=str(target.parent), prefix=f".{target.stem}_", suffix=".csv.tmp")
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(tmp_fd, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, lineterminator="\n")
            if original_text:
                f.write(original_text)
                if not original_text.endswith(("\n", "\r")):
                    f.write("\n")
            else:
                writer.writerow(header)
            writer.writerows(new_rows)
            f.flush()
            os.fsync(f.fileno())
        tmp_path.replace(target)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink(missing_ok=True)
        raise
    return len(new_rows)


# =============================================================================
# Step 02: Atomic Per-Voucher Segmentation Transactions & O(1) Resume Queries
# =============================================================================

def get_completed_qc_vouchers(
    conn: sqlite3.Connection,
    catalog_filter: Optional[Iterable[str]] = None,
) -> Set[str]:
    """Query completed vouchers in O(1) via indexed `qc_status IN ('PASS', 'FAILED_QC')`."""
    cur = conn.execute(
        "SELECT catalogNumber FROM vouchers WHERE qc_status IN ('PASS', 'FAILED_QC');"
    )
    done = {str(r[0]).strip() for r in cur.fetchall() if r[0]}
    if catalog_filter is not None:
        allowed = {str(c).strip() for c in catalog_filter}
        done &= allowed
    return done


def read_contour_csv_xy(contour_path: Union[str, Path]) -> Optional[List[List[float]]]:
    """
    Read a 120-vertex contour CSV from disk and return compact `[[x, y], ...]`
    coordinates (preferring raw pixel `x, y` over `x_norm, y_norm` so Step 03
    metric scale calibration `coo / resolved_ppm` and `area_px` remain exact).
    """
    p = Path(contour_path)
    if not p.is_file() or p.stat().st_size == 0:
        return None
    try:
        with open(p, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            fnames = set(reader.fieldnames or [])
            if {"x", "y"}.issubset(fnames):
                xk, yk = "x", "y"
            elif {"x_norm", "y_norm"}.issubset(fnames):
                xk, yk = "x_norm", "y_norm"
            else:
                return None
            pts: List[List[float]] = []
            for row in reader:
                fx = _to_float(row.get(xk))
                fy = _to_float(row.get(yk))
                if fx is not None and fy is not None:
                    pts.append([round(fx, 4), round(fy, 4)])
            return pts if len(pts) >= 3 else None
    except OSError:
        return None


def _resolve_contour_file(
    raw_contour_path: Optional[str],
    output_dir: Optional[Path] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> Optional[Path]:
    if not raw_contour_path:
        return None
    cp = Path(raw_contour_path)
    if cp.is_file():
        return cp
    fname = cp.name
    candidates: List[Path] = []
    if output_dir is not None:
        candidates.extend([
            output_dir / "contours" / fname,
            output_dir / "data" / "contours" / fname,
            output_dir.parent / raw_contour_path,
        ])
    candidates.extend([
        workspace_root / raw_contour_path,
        workspace_root / "data" / "contours" / fname,
    ])
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def record_voucher_segmentation_transaction(
    conn: sqlite3.Connection,
    catalog_number: str,
    extracted_records: Sequence[Mapping[str, Any]],
    failure_record: Optional[Mapping[str, Any]] = None,
    repro_record: Optional[Mapping[str, Any]] = None,
    pheno_record: Optional[Mapping[str, Any]] = None,
    crowns_by_pid: Optional[Mapping[int, Tuple[float, float]]] = None,
    rosette_centroids_by_pid: Optional[Mapping[int, Tuple[float, float]]] = None,
    by_plant_petiole: Optional[Mapping[int, Mapping[str, Any]]] = None,
    contour_coords_by_leaf: Optional[Mapping[int, Sequence[Sequence[float]]]] = None,
    output_dir: Optional[Path] = None,
    image_path: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
    commit: bool = True,
) -> None:
    """
    Persist a single voucher's Step 02 segmentation results in one atomic SQLite
    transaction across `vouchers`, `plant_individuals`, `extracted_leaves`, and
    `voucher_reproductive_phenology`.
    """
    cat = str(catalog_number).strip()
    if not cat:
        return

    repro = dict(repro_record or {})
    pheno = dict(pheno_record or {})
    crowns = dict(crowns_by_pid or {})
    centroids = dict(rosette_centroids_by_pid or {})
    petiole_by_pid = dict(by_plant_petiole or {})
    coords_map = dict(contour_coords_by_leaf or {})

    # Determine voucher-level calibration & QC status
    ppm_val: Optional[float] = None
    is_cal_val: Optional[int] = None
    for r in extracted_records:
        p = _to_float(r.get("pixels_per_mm"))
        if p is not None and p > 0:
            ppm_val = p
            is_cal_val = 1
            break
        c_flag = _to_bool_int(r.get("is_calibrated_mm"))
        if c_flag is not None and is_cal_val is None:
            is_cal_val = c_flag
    if ppm_val is None and failure_record:
        p_fail = _to_float(failure_record.get("pixels_per_mm"))
        if p_fail is not None and p_fail > 0:
            ppm_val = p_fail
            is_cal_val = 1
        elif is_cal_val is None:
            c_fail = _to_bool_int(failure_record.get("is_calibrated_mm"))
            if c_fail is not None:
                is_cal_val = c_fail

    qc_status = "PASS" if len(extracted_records) > 0 else "FAILED_QC"
    fail_reason = _to_str(failure_record.get("failure_reason")) if failure_record else None
    fail_details = _to_str(failure_record.get("details") or failure_record.get("failure_details")) if failure_record else None
    eff_img_path = _to_str(image_path) or (_to_str(failure_record.get("image_path")) if failure_record else None)

    # Collect all distinct plant_individual_id values (never collapsing multi-plant sheets to 0)
    pids: Set[int] = set()
    for r in extracted_records:
        pids.add(_to_int(r.get("plant_individual_id"), default=0) or 0)
    for k_pid in crowns.keys():
        pids.add(int(k_pid))
    for k_pid in centroids.keys():
        pids.add(int(k_pid))
    if not pids:
        pids.add(0)

    out_base = Path(output_dir) if output_dir is not None else (workspace_root / "data")

    plant_rows: List[Tuple[Any, ...]] = []
    repro_pheno_rows: List[Tuple[Any, ...]] = []
    leaf_rows: List[Tuple[Any, ...]] = []

    p_height_cm = _to_float(repro.get("plant_height_cm"))
    p_height_mm = (p_height_cm * 10.0) if p_height_cm is not None else _to_float(repro.get("plant_height_mm"))
    p_height_px = _to_float(repro.get("plant_height_px"))
    if p_height_px is None and p_height_mm is not None and ppm_val:
        p_height_px = p_height_mm * ppm_val
    scape_ratio = _to_float(repro.get("scape_to_basal_leaf_ratio"))
    is_bent = _to_bool_int(repro.get("is_bent_stem"), default=0)

    # Group leaf-level petiole/blade metrics by pid as fallback
    leaves_by_pid: Dict[int, List[Mapping[str, Any]]] = {pid: [] for pid in pids}
    for r in extracted_records:
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        leaves_by_pid.setdefault(pid, []).append(r)

    for pid in sorted(pids):
        rc_xy = centroids.get(pid)
        cr_xy = crowns.get(pid)
        p_pet_info = dict(petiole_by_pid.get(pid) or {})

        pet_mm = _to_float(p_pet_info.get("petiole_length_mm"))
        bld_mm = _to_float(p_pet_info.get("blade_length_mm"))
        p2b = _to_float(p_pet_info.get("petiole_to_blade_ratio"))
        pet_px = _to_float(p_pet_info.get("petiole_length_px"))
        bld_px = _to_float(p_pet_info.get("blade_length_px"))

        if pet_mm is None and leaves_by_pid.get(pid):
            vals = [_to_float(x.get("petiole_length_mm")) for x in leaves_by_pid[pid]]
            vals_clean = sorted(v for v in vals if v is not None)
            if vals_clean:
                pet_mm = vals_clean[min(len(vals_clean) - 1, int(0.75 * (len(vals_clean) - 1)))]
        if bld_mm is None and leaves_by_pid.get(pid):
            vals = [_to_float(x.get("blade_length_mm")) for x in leaves_by_pid[pid]]
            vals_clean = sorted(v for v in vals if v is not None)
            if vals_clean:
                bld_mm = vals_clean[min(len(vals_clean) - 1, int(0.75 * (len(vals_clean) - 1)))]
        if p2b is None and leaves_by_pid.get(pid):
            vals = [_to_float(x.get("petiole_to_blade_ratio")) for x in leaves_by_pid[pid]]
            vals_clean = sorted(v for v in vals if v is not None)
            if vals_clean:
                p2b = vals_clean[min(len(vals_clean) - 1, int(0.75 * (len(vals_clean) - 1)))]

        r_patch: Optional[str] = None
        r_mask: Optional[str] = None

        plant_rows.append((
            cat,
            pid,
            _to_float(rc_xy[0]) if rc_xy else None,
            _to_float(rc_xy[1]) if rc_xy else None,
            _to_float(cr_xy[0]) if cr_xy else None,
            _to_float(cr_xy[1]) if cr_xy else None,
            p_height_px,
            p_height_mm,
            p_height_cm,
            is_bent,
            pet_px,
            pet_mm,
            bld_px,
            bld_mm,
            p2b,
            scape_ratio,
            r_patch,
            r_mask,
        ))

        if repro or pheno:
            ar_inv = _to_float(repro.get("AR_inv"))
            cap_ar = _to_float(repro.get("capitulum_aspect_ratio"))
            if ar_inv is None:
                ar_inv = cap_ar
            if cap_ar is None:
                cap_ar = ar_inv
            n_fr = _to_int(pheno.get("n_fruit"))
            n_fr_ing = _to_int(pheno.get("n_fruiting"))
            if n_fr is None:
                n_fr = n_fr_ing
            if n_fr_ing is None:
                n_fr_ing = n_fr
            ph_state = _to_str(pheno.get("voucher_phenological_state"))
            is_ver_anth = _to_bool_int(pheno.get("is_verified_anthesis"))
            if is_ver_anth is None and ph_state:
                low_st = ph_state.lower()
                is_ver_anth = 1 if ("anthesis" in low_st and "pre" not in low_st) else 0

            repro_pheno_rows.append((
                cat,
                pid,
                _to_int(repro.get("capitula_count")),
                _to_float(repro.get("involucre_height_px")),
                _to_float(repro.get("involucre_width_px")),
                _to_float(repro.get("involucre_height_mm")),
                _to_float(repro.get("involucre_width_mm")),
                ar_inv,
                cap_ar,
                _to_int(pheno.get("total_capitula"), default=_to_int(repro.get("capitula_count"))),
                _to_int(pheno.get("n_anthesis")),
                _to_int(pheno.get("n_bud")),
                n_fr_ing,
                n_fr,
                _to_int(pheno.get("n_senescent")),
                ph_state,
                is_ver_anth,
                _to_float(pheno.get("reproductive_stage_confidence")),
                _to_float(pheno.get("dominant_pappus_ratio")),
            ))

    for r in extracted_records:
        leaf_id = _to_int(r.get("leaf_id"), default=1) or 1
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        assigned_tier = _normalize_tier(r.get("assigned_tier") or r.get("tier_assignment") or r.get("tier"))
        tier_str = _to_str(r.get("tier")) or ("Tier 2" if assigned_tier == "Tier_2_Reflected" else "Tier 1")
        ucs_val = _to_float(r.get("ucs_score"))
        if ucs_val is None:
            ucs_val = _to_float(r.get("ucs"))

        xy_pts: Optional[List[List[float]]] = None
        if r.get("contour_xy_json"):
            raw_json = str(r["contour_xy_json"])
            try:
                parsed = json.loads(raw_json)
                if isinstance(parsed, list) and len(parsed) >= 3:
                    xy_pts = [[round(float(pt[0]), 4), round(float(pt[1]), 4)] for pt in parsed]
            except (ValueError, TypeError, KeyError, IndexError):
                xy_pts = None
        if xy_pts is None and leaf_id in coords_map:
            seq = coords_map[leaf_id]
            xy_pts = [[round(float(pt[0]), 4), round(float(pt[1]), 4)] for pt in seq]
        if xy_pts is None:
            c_file = _resolve_contour_file(_to_str(r.get("contour_path")), output_dir=out_base, workspace_root=workspace_root)
            if c_file is not None:
                xy_pts = read_contour_csv_xy(c_file)
        if not xy_pts:
            # Skip unbacked leaf rows lacking real contour vertices (Rule 2: Zero Silent Synthetic Data)
            continue

        contour_json = json.dumps(xy_pts, separators=(",", ":"))
        crop_p = _to_str(r.get("crop_path"))

        shape_id = _to_str(r.get("shape_id")) or f"{cat}_p{pid}_leaf{leaf_id}"
        l_ppm = _to_float(r.get("pixels_per_mm"))
        if l_ppm is None:
            l_ppm = ppm_val
        l_cal = _to_bool_int(r.get("is_calibrated_mm"))
        if l_cal is None:
            l_cal = 1 if (l_ppm is not None and l_ppm > 0) else 0

        leaf_rows.append((
            cat,
            pid,
            leaf_id,
            shape_id,
            tier_str,
            assigned_tier,
            _to_str(r.get("tier_assignment")) or assigned_tier,
            _to_float(r.get("ucs")) or ucs_val,
            ucs_val,
            _to_float(r.get("solidity")),
            _to_float(r.get("circularity")),
            _to_float(r.get("extent")),
            _to_float(r.get("aspect_ratio")),
            _to_float(r.get("area_px")),
            _to_float(r.get("blade_area_cm2")),
            _to_bool_int(r.get("is_folded"), default=0),
            _to_bool_int(r.get("is_dissected"), default=0),
            _to_bool_int(r.get("reflection_applied"), default=(1 if assigned_tier == "Tier_2_Reflected" else 0)),
            _to_float(r.get("midrib_angle_deg")),
            _to_float(r.get("petiole_length_px")),
            _to_float(r.get("petiole_length_mm")),
            _to_float(r.get("blade_length_px")),
            _to_float(r.get("blade_length_mm")),
            _to_float(r.get("petiole_to_blade_ratio")),
            l_ppm,
            l_cal,
            _to_str(r.get("mask_path")),
            crop_p,
            _to_str(r.get("contour_path")),
            len(xy_pts),
            contour_json,
        ))

    if extracted_records and not leaf_rows:
        # All extracted records lacked readable contours
        qc_status = "FAILED_QC"
        fail_reason = fail_reason or "MISSING_CONTOUR_COORDINATES"

    def _apply_sql() -> None:
        conn.execute(
            """
            INSERT INTO vouchers (catalogNumber, image_path, pixels_per_mm, is_calibrated_mm, qc_status, failure_reason, failure_details)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(catalogNumber) DO UPDATE SET
                image_path = COALESCE(vouchers.image_path, excluded.image_path),
                pixels_per_mm = COALESCE(excluded.pixels_per_mm, vouchers.pixels_per_mm),
                is_calibrated_mm = COALESCE(excluded.is_calibrated_mm, vouchers.is_calibrated_mm),
                qc_status = excluded.qc_status,
                failure_reason = excluded.failure_reason,
                failure_details = excluded.failure_details;
            """,
            (cat, eff_img_path, ppm_val, is_cal_val, qc_status, fail_reason, fail_details),
        )
        conn.execute("DELETE FROM extracted_leaves WHERE catalogNumber = ?;", (cat,))
        conn.execute("DELETE FROM voucher_reproductive_phenology WHERE catalogNumber = ?;", (cat,))
        conn.execute("DELETE FROM plant_individuals WHERE catalogNumber = ?;", (cat,))

        conn.executemany(
            """
            INSERT INTO plant_individuals (
                catalogNumber, plant_individual_id,
                rosette_centroid_x, rosette_centroid_y, crown_x, crown_y,
                plant_height_px, plant_height_mm, plant_height_cm, is_bent_stem,
                petiole_length_px, petiole_length_mm, blade_length_px, blade_length_mm,
                petiole_to_blade_ratio, scape_to_basal_leaf_ratio,
                rosette_patch_path, rosette_mask_path
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
            """,
            plant_rows,
        )
        if repro_pheno_rows:
            conn.executemany(
                """
                INSERT INTO voucher_reproductive_phenology (
                    catalogNumber, plant_individual_id,
                    capitula_count, involucre_height_px, involucre_width_px,
                    involucre_height_mm, involucre_width_mm, AR_inv, capitulum_aspect_ratio,
                    total_capitula, n_anthesis, n_bud, n_fruiting, n_fruit, n_senescent,
                    voucher_phenological_state, is_verified_anthesis,
                    reproductive_stage_confidence, dominant_pappus_ratio
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                """,
                repro_pheno_rows,
            )
        if leaf_rows:
            conn.executemany(
                """
                INSERT INTO extracted_leaves (
                    catalogNumber, plant_individual_id, leaf_id, shape_id,
                    tier, assigned_tier, tier_assignment,
                    ucs, ucs_score, solidity, circularity, extent, aspect_ratio,
                    area_px, blade_area_cm2, is_folded, is_dissected, reflection_applied,
                    midrib_angle_deg, petiole_length_px, petiole_length_mm,
                    blade_length_px, blade_length_mm, petiole_to_blade_ratio,
                    pixels_per_mm, is_calibrated_mm,
                    mask_path, crop_path, contour_path, contour_n_pts, contour_xy_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
                """,
                leaf_rows,
            )

    if commit:
        with conn:
            _apply_sql()
    else:
        _apply_sql()


def load_step02_records_from_db(
    conn: sqlite3.Connection,
    catalog_numbers: Sequence[str],
) -> Tuple[
    Dict[str, List[Dict[str, Any]]],
    Dict[str, Dict[str, Any]],
    Dict[str, Dict[str, Any]],
    Dict[str, Dict[str, Any]],
]:
    """
    Fetch existing Step 02 leaf, failed-QC, reproductive, and phenology records
    from SQLite for a list of `catalog_numbers` when resuming Step 02.
    """
    cats = [str(c).strip() for c in catalog_numbers if str(c).strip()]
    extracted_by_cat: Dict[str, List[Dict[str, Any]]] = {}
    failed_by_cat: Dict[str, Dict[str, Any]] = {}
    repro_by_cat: Dict[str, Dict[str, Any]] = {}
    pheno_by_cat: Dict[str, Dict[str, Any]] = {}
    if not cats:
        return extracted_by_cat, failed_by_cat, repro_by_cat, pheno_by_cat

    chunk_size = 500
    for i in range(0, len(cats), chunk_size):
        chunk = cats[i : i + chunk_size]
        ph = ", ".join("?" for _ in chunk)

        cur_l = conn.execute(
            f"""
            SELECT
                el.catalogNumber, el.leaf_id, el.tier, el.solidity,
                COALESCE(el.ucs, el.ucs_score) AS ucs,
                el.is_folded, el.is_dissected, el.reflection_applied, el.contour_path,
                el.plant_individual_id, el.assigned_tier,
                COALESCE(el.tier_assignment, el.assigned_tier) AS tier_assignment,
                el.ucs_score, el.midrib_angle_deg, el.pixels_per_mm, el.is_calibrated_mm,
                el.mask_path,
                COALESCE(el.petiole_length_mm, pi.petiole_length_mm) AS petiole_length_mm,
                COALESCE(el.blade_length_mm, pi.blade_length_mm) AS blade_length_mm,
                COALESCE(el.petiole_to_blade_ratio, pi.petiole_to_blade_ratio) AS petiole_to_blade_ratio,
                vrp.capitula_count, vrp.involucre_height_px, vrp.involucre_width_px,
                vrp.involucre_height_mm, vrp.involucre_width_mm,
                COALESCE(vrp.capitulum_aspect_ratio, vrp.AR_inv) AS capitulum_aspect_ratio
            FROM extracted_leaves el
            LEFT JOIN plant_individuals pi
                ON pi.catalogNumber = el.catalogNumber
               AND pi.plant_individual_id = el.plant_individual_id
            LEFT JOIN voucher_reproductive_phenology vrp
                ON vrp.catalogNumber = el.catalogNumber
               AND vrp.plant_individual_id = el.plant_individual_id
            WHERE el.catalogNumber IN ({ph})
            ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;
            """,
            chunk,
        )
        for row in cur_l.fetchall():
            d = dict(row)
            d["is_folded"] = bool(d.get("is_folded"))
            d["is_dissected"] = bool(d.get("is_dissected"))
            d["reflection_applied"] = bool(d.get("reflection_applied"))
            d["is_calibrated_mm"] = bool(d.get("is_calibrated_mm"))
            extracted_by_cat.setdefault(str(d["catalogNumber"]), []).append(d)

        cur_f = conn.execute(
            f"""
            SELECT catalogNumber, image_path, failure_reason, failure_details AS details
            FROM vouchers
            WHERE qc_status = 'FAILED_QC' AND catalogNumber IN ({ph});
            """,
            chunk,
        )
        for row in cur_f.fetchall():
            d = dict(row)
            failed_by_cat[str(d["catalogNumber"])] = d

        cur_rp = conn.execute(
            f"""
            SELECT
                vrp.catalogNumber,
                vrp.capitula_count, vrp.involucre_height_px, vrp.involucre_width_px,
                vrp.involucre_height_mm, vrp.involucre_width_mm,
                COALESCE(vrp.capitulum_aspect_ratio, vrp.AR_inv) AS capitulum_aspect_ratio,
                pi.plant_height_cm, pi.scape_to_basal_leaf_ratio, pi.is_bent_stem,
                vrp.total_capitula, vrp.n_anthesis, COALESCE(vrp.n_fruit, vrp.n_fruiting) AS n_fruit,
                vrp.n_bud, vrp.voucher_phenological_state, vrp.dominant_pappus_ratio
            FROM voucher_reproductive_phenology vrp
            LEFT JOIN plant_individuals pi
                ON pi.catalogNumber = vrp.catalogNumber
               AND pi.plant_individual_id = vrp.plant_individual_id
            WHERE vrp.catalogNumber IN ({ph})
            ORDER BY vrp.plant_individual_id ASC;
            """,
            chunk,
        )
        for row in cur_rp.fetchall():
            d = dict(row)
            c = str(d["catalogNumber"])
            if c not in repro_by_cat:
                repro_by_cat[c] = {
                    "catalogNumber": c,
                    "capitula_count": d.get("capitula_count"),
                    "involucre_height_px": d.get("involucre_height_px"),
                    "involucre_width_px": d.get("involucre_width_px"),
                    "involucre_height_mm": d.get("involucre_height_mm"),
                    "involucre_width_mm": d.get("involucre_width_mm"),
                    "capitulum_aspect_ratio": d.get("capitulum_aspect_ratio"),
                    "plant_height_cm": d.get("plant_height_cm"),
                    "scape_to_basal_leaf_ratio": d.get("scape_to_basal_leaf_ratio"),
                    "is_bent_stem": bool(d.get("is_bent_stem")),
                }
            if c not in pheno_by_cat and d.get("voucher_phenological_state") is not None:
                pheno_by_cat[c] = {
                    "catalogNumber": c,
                    "total_capitula": d.get("total_capitula"),
                    "n_anthesis": d.get("n_anthesis"),
                    "n_fruit": d.get("n_fruit"),
                    "n_bud": d.get("n_bud"),
                    "voucher_phenological_state": d.get("voucher_phenological_state"),
                    "dominant_pappus_ratio": d.get("dominant_pappus_ratio"),
                }

    return extracted_by_cat, failed_by_cat, repro_by_cat, pheno_by_cat


# =============================================================================
# Run Registration & SQL Preflight Cohort Selection
# =============================================================================

def ensure_pipeline_run(
    conn: sqlite3.Connection,
    run_id: str,
    run_dir: Optional[Union[str, Path]] = None,
    start_timestamp: Optional[str] = None,
    end_timestamp: Optional[str] = None,
    git_commit_sha: Optional[str] = None,
    cli_arguments: Optional[Sequence[str]] = None,
    focal_taxon: Optional[str] = None,
    intraspecific_mode: bool = False,
) -> None:
    """Register or update a run record inside `pipeline_runs`."""
    rid = str(run_id).strip() or "canonical"
    rdir_str = str(run_dir) if run_dir is not None else None
    cli_json = json.dumps(list(cli_arguments)) if cli_arguments is not None else None
    now_iso = datetime.now(timezone.utc).isoformat()

    with conn:
        conn.execute(
            """
            INSERT INTO pipeline_runs (
                run_id, run_dir, start_timestamp, end_timestamp,
                git_commit_sha, cli_arguments_json, focal_taxon, intraspecific_mode
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id) DO UPDATE SET
                run_dir = COALESCE(excluded.run_dir, pipeline_runs.run_dir),
                start_timestamp = COALESCE(pipeline_runs.start_timestamp, excluded.start_timestamp),
                end_timestamp = COALESCE(excluded.end_timestamp, pipeline_runs.end_timestamp),
                git_commit_sha = COALESCE(excluded.git_commit_sha, pipeline_runs.git_commit_sha),
                cli_arguments_json = COALESCE(excluded.cli_arguments_json, pipeline_runs.cli_arguments_json),
                focal_taxon = COALESCE(excluded.focal_taxon, pipeline_runs.focal_taxon),
                intraspecific_mode = COALESCE(excluded.intraspecific_mode, pipeline_runs.intraspecific_mode);
            """,
            (
                rid,
                rdir_str,
                start_timestamp or now_iso,
                end_timestamp,
                git_commit_sha,
                cli_json,
                _to_str(focal_taxon),
                1 if intraspecific_mode else 0,
            ),
        )


def populate_run_active_cohort(
    conn: sqlite3.Connection,
    run_id: str,
    voucher_rows: Sequence[Mapping[str, Any]],
    clear_existing: bool = True,
) -> int:
    """
    Populate `run_active_cohort` for `run_id`, expanding every retained voucher
    across all of its `(catalogNumber, plant_individual_id)` rows in `plant_individuals`
    (or `plant_individual_id = 0` for unsegmented/historical sheets).
    """
    rid = str(run_id).strip() or "canonical"
    ensure_pipeline_run(conn, run_id=rid)
    _tier_b_tables = (
        "leaf_efa_harmonics",
        "otu_morphometrics",
        "spatial_rf_predictions",
        "curation_triage",
    )
    if not voucher_rows:
        if clear_existing:
            with conn:
                conn.execute("DELETE FROM run_active_cohort WHERE run_id = ?;", (rid,))
                for _tbl in _tier_b_tables:
                    conn.execute(f"DELETE FROM {_tbl} WHERE run_id = ?;", (rid,))
        return 0

    _ensure_dynamic_table_migrated(conn, "run_active_cohort")
    upsert_vouchers_records(conn, voucher_rows)

    cats = [str(r.get("catalogNumber", "")).strip() for r in voucher_rows if str(r.get("catalogNumber", "")).strip()]
    pids_by_cat: Dict[str, List[int]] = {}
    chunk_size = 500
    for i in range(0, len(cats), chunk_size):
        chunk = cats[i : i + chunk_size]
        ph = ", ".join("?" for _ in chunk)
        cur = conn.execute(
            f"SELECT catalogNumber, plant_individual_id FROM plant_individuals WHERE catalogNumber IN ({ph});",
            chunk,
        )
        for row in cur.fetchall():
            pids_by_cat.setdefault(str(row[0]).strip(), []).append(int(row[1]))

    cohort_payload: List[Tuple[Any, ...]] = []
    seen_keys: Set[Tuple[str, int]] = set()
    for r in voucher_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        sp_std = _to_str(r.get("species_standardized"))
        det_tier = _to_str(r.get("determiner_tier"))
        type_stat = _to_str(r.get("type_status"))
        is_type_sp = _to_bool_int(r.get("is_type_specimen"))
        if is_type_sp is None:
            is_type_sp = (
                1
                if (
                    type_stat
                    and type_stat.upper()
                    not in {
                        "NONE", "NOT A TYPE", "NOTATYPE", "NON-TYPE", "NONTYPE",
                        "UNSPECIFIED", "NULL", "NAN", "NA", "N/A", "FALSE", "0", "UNKNOWN",
                    }
                    and not re.search(r"\b(?:not\s+a?\s*type|non[-\s]?type)\b", type_stat, re.IGNORECASE)
                )
                else 0
            )
        is_prom = _to_bool_int(r.get("is_promoted_bronze") or r.get("is_promoted_anchor"), default=0) or 0
        auth_wt = _to_float(r.get("authority_weight"))
        sp_wt = _to_float(r.get("spatial_balance_weight"))
        if sp_wt is None:
            sp_wt = 1.0

        explicit_pid = _to_int(r.get("plant_individual_id"))
        target_pids = [explicit_pid] if explicit_pid is not None else sorted(set(pids_by_cat.get(cat, [0])))
        for pid in target_pids:
            if (cat, pid) in seen_keys:
                continue
            seen_keys.add((cat, pid))
            cohort_payload.append((rid, cat, pid, sp_std, det_tier, type_stat, is_type_sp, is_prom, auth_wt, sp_wt))

    with conn:
        if clear_existing:
            conn.execute("DELETE FROM run_active_cohort WHERE run_id = ?;", (rid,))
        conn.executemany(
            """
            INSERT INTO run_active_cohort (
                run_id, catalogNumber, plant_individual_id,
                species_standardized, determiner_tier,
                type_status, is_type_specimen,
                is_promoted_bronze, authority_weight, spatial_balance_weight
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id, catalogNumber, plant_individual_id) DO UPDATE SET
                species_standardized = COALESCE(excluded.species_standardized, run_active_cohort.species_standardized),
                determiner_tier = COALESCE(excluded.determiner_tier, run_active_cohort.determiner_tier),
                type_status = COALESCE(excluded.type_status, run_active_cohort.type_status),
                is_type_specimen = COALESCE(excluded.is_type_specimen, run_active_cohort.is_type_specimen),
                is_promoted_bronze = COALESCE(excluded.is_promoted_bronze, run_active_cohort.is_promoted_bronze),
                authority_weight = COALESCE(excluded.authority_weight, run_active_cohort.authority_weight),
                spatial_balance_weight = COALESCE(excluded.spatial_balance_weight, run_active_cohort.spatial_balance_weight);
            """,
            cohort_payload,
        )
        if clear_existing:
            for _tbl in _tier_b_tables:
                conn.execute(
                    f"""
                    DELETE FROM {_tbl}
                    WHERE run_id = ?
                      AND NOT EXISTS (
                          SELECT 1 FROM run_active_cohort rac
                          WHERE rac.run_id = {_tbl}.run_id
                            AND rac.catalogNumber = {_tbl}.catalogNumber
                            AND rac.plant_individual_id = {_tbl}.plant_individual_id
                      );
                    """,
                    (rid,),
                )
            conn.execute(
                """
                UPDATE otu_morphometrics
                SET species_standardized = (
                    SELECT rac.species_standardized
                    FROM run_active_cohort rac
                    WHERE rac.run_id = otu_morphometrics.run_id
                      AND rac.catalogNumber = otu_morphometrics.catalogNumber
                      AND rac.plant_individual_id = otu_morphometrics.plant_individual_id
                )
                WHERE run_id = ?
                  AND EXISTS (
                      SELECT 1
                      FROM run_active_cohort rac
                      WHERE rac.run_id = otu_morphometrics.run_id
                        AND rac.catalogNumber = otu_morphometrics.catalogNumber
                        AND rac.plant_individual_id = otu_morphometrics.plant_individual_id
                        AND rac.species_standardized IS NOT NULL
                        AND TRIM(rac.species_standardized) <> ''
                  );
                """,
                (rid,),
            )
    return len(cohort_payload)


def query_misidentified_vouchers_from_db(
    conn: sqlite3.Connection,
    triage_sources: Optional[Sequence[str]] = None,
    exclude_run_id: Optional[str] = None,
    allowed_subtypes: Optional[Set[str]] = None,
) -> Tuple[Set[str], int]:
    """
    Query `curation_triage` directly in SQLite to identify vouchers flagged as
    `Misidentification_Reassign` (or outside `allowed_subtypes`) for `--exclude-misidentified`.
    Returns `(flagged_catalog_numbers, matched_run_count)`.
    """
    allowed = allowed_subtypes or {
        "Confirmed_Species",
        "Confirmed_Type_Control",
        "Type_Synonymy_Candidate",
        "Hybrid_Introgressant",
        "Subspecies_Ecotype",
        "Ecophenotypic_Plasticity",
        "Putative_Cryptic_Species",
    }
    tokens = [str(t).strip() for t in (triage_sources or ["ALL"]) if str(t).strip()]
    if not tokens:
        tokens = ["ALL"]

    where_clauses: List[str] = []
    params: List[Any] = []
    if exclude_run_id:
        where_clauses.append("run_id != ?")
        params.append(str(exclude_run_id).strip())

    if not any(t.upper() in {"ALL", "*"} for t in tokens):
        run_ids = [Path(t).name if "/" in t else t for t in tokens]
        ph = ", ".join("?" for _ in run_ids)
        where_clauses.append(f"run_id IN ({ph})")
        params.extend(run_ids)

    where_sql = (" WHERE " + " AND ".join(where_clauses)) if where_clauses else ""
    cur_runs = conn.execute(f"SELECT COUNT(DISTINCT run_id) FROM curation_triage{where_sql};", params)
    n_runs = int((cur_runs.fetchone() or [0])[0] or 0)
    if n_runs == 0:
        return set(), 0

    cur = conn.execute(
        f"SELECT catalogNumber, taxonomic_treatment_subtype FROM curation_triage{where_sql};",
        params,
    )
    flagged: Set[str] = set()
    for row in cur.fetchall():
        cat = str(row[0] or "").strip()
        sub = str(row[1] or "").strip()
        if not cat or not sub:
            continue
        if sub == "Misidentification_Reassign" or sub not in allowed:
            flagged.add(cat)
    return flagged, n_runs


# =============================================================================
# Environmental Cache Upsert & R*Tree Spatial Bounding-Box Query
# =============================================================================

def _coalesce_wc_bio_fields(src: Mapping[str, Any], dst: Dict[str, Any]) -> None:
    """
    Bidirectionally coalesce canonical WorldClim `_qtr` (`bio18_precip_warmest_qtr`,
    `bio19_precip_coldest_qtr`) and legacy `_qt` (`bio18_precip_warmest_qt`,
    `bio19_precip_coldest_qt`) fields so 19-variable runs and legacy 6-variable
    offline cache fixtures interoperate seamlessly.
    """
    for qtr_col, qt_col in (
        ("bio18_precip_warmest_qtr", "bio18_precip_warmest_qt"),
        ("bio19_precip_coldest_qtr", "bio19_precip_coldest_qt"),
    ):
        if qtr_col in src or qt_col in src or qtr_col in dst or qt_col in dst:
            v_qtr = _to_float(dst.get(qtr_col) if qtr_col in dst else src.get(qtr_col))
            v_qt = _to_float(dst.get(qt_col) if qt_col in dst else src.get(qt_col))
            resolved_qtr = v_qtr if v_qtr is not None else v_qt
            resolved_qt = v_qt if v_qt is not None else v_qtr
            if resolved_qtr is not None or qtr_col in src or qt_col in src:
                dst[qtr_col] = resolved_qtr
                dst[qt_col] = resolved_qt


def _normalize_semicolon_list_str(val: Any) -> Optional[str]:
    if val is None:
        return None
    if isinstance(val, (list, tuple, set)):
        items = [str(x).strip() for x in val if not _is_missing(x) and str(x).strip()]
        return ";".join(items) if items else None
    return _to_str(val)


def _ensure_dynamic_table_migrated(conn: sqlite3.Connection, table_name: str) -> None:
    specs = _DYNAMIC_TABLE_MIGRATION_SPECS.get(table_name)
    if not specs:
        return
    existing = set(get_table_columns(conn, table_name))
    if existing and any(col not in existing for col, _ in specs):
        _migrate_dynamic_schema_columns(conn)


def upsert_environmental_cache_records(
    conn: sqlite3.Connection,
    records: Sequence[Mapping[str, Any]],
    demo: bool = False,
) -> int:
    """
    UPSERT fine-grain SSURGO / POLARIS / SOLUS100 / 3DEP / SGMC / 19-variable WorldClim
    rows into `environmental_cache` (automatically updating `environmental_rtree` via
    triggers and coalescing `bio18/19` `_qtr <-> _qt`). Strictly no-ops when `demo=True`.
    """
    if demo or not records:
        return 0
    _ensure_dynamic_table_migrated(conn, "environmental_cache")
    norm_rows: List[Dict[str, Any]] = []
    for r in records:
        lat = _to_float(r.get("latitude") if "latitude" in r else r.get("decimalLatitude"))
        lon = _to_float(r.get("longitude") if "longitude" in r else r.get("decimalLongitude"))
        ll_key = _to_str(r.get("lat_lon_key")) or format_lat_lon_key(lat, lon)
        if not ll_key:
            continue
        d: Dict[str, Any] = {
            "lat_lon_key": ll_key,
            "catalogNumber": _to_str(r.get("catalogNumber")),
            "latitude": lat,
            "longitude": lon,
        }
        for col in ENVIRONMENTAL_CACHE_COLUMNS:
            if col in {"lat_lon_key", "catalogNumber", "latitude", "longitude"}:
                continue
            if col not in r:
                continue
            if col in ENVIRONMENTAL_TEXT_COLS:
                d[col] = _to_str(r.get(col))
            elif col in ENVIRONMENTAL_BOOL_COLS:
                d[col] = _to_bool_int(r.get(col))
            else:
                d[col] = _to_float(r.get(col))
        _coalesce_wc_bio_fields(r, d)
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="environmental_cache",
        rows=norm_rows,
        pk_cols=("lat_lon_key",),
        preserve_existing_on_null=True,
    )


def query_environmental_bbox(
    conn: sqlite3.Connection,
    min_lon: float,
    max_lon: float,
    min_lat: float,
    max_lat: float,
) -> List[Dict[str, Any]]:
    """Query `environmental_cache` in O(log N) using the `environmental_rtree` spatial index."""
    cur = conn.execute(
        """
        SELECT ec.*
        FROM environmental_cache ec
        INNER JOIN environmental_rtree rt ON ec.rowid = rt.id
        WHERE rt.min_lon >= ? AND rt.max_lon <= ?
          AND rt.min_lat >= ? AND rt.max_lat <= ?;
        """,
        (float(min_lon), float(max_lon), float(min_lat), float(max_lat)),
    )
    out: List[Dict[str, Any]] = []
    for r in cur.fetchall():
        d = dict(r)
        _coalesce_wc_bio_fields(d, d)
        out.append(d)
    return out


# =============================================================================
# Steps 03, 04 & 06: Dynamic Morphometric & Spatial RF Persistence Helpers
# =============================================================================

_LEAF_EFA_BOOL_COLS: Tuple[str, ...] = (
    "is_dissected",
    "is_calibrated_mm",
    "is_contour_outlier",
)
_LEAF_EFA_INT_COLS: Tuple[str, ...] = (
    "efd_harmonic_opt",
    "efd_k95_pcs",
)
_LEAF_EFA_STR_COLS: Tuple[str, ...] = ("assigned_tier",)
_LEAF_EFA_FLOAT_COLS: Tuple[str, ...] = (
    "solidity",
    "ucs_score",
    "aspect_ratio",
    "area_px",
    "pixels_per_mm",
    "petiole_length_mm",
    "blade_length_mm",
    "petiole_to_blade_ratio",
    *HARMONIC_COLS,
    "sym_var",
    "asym_var",
    "pct_asymmetry",
    "harmonic_roughness",
    "leaf_quality_weight",
    "centroid_size",
    "log_centsize",
    *PC_COLS,
    *SYM_PC_COLS,
)


def persist_leaf_efa_harmonics(
    conn: sqlite3.Connection,
    run_id: str,
    efa_rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
) -> int:
    """
    Persist Step 03 leaf-level Elliptic Fourier harmonics (`A1..D24`), leaf PCs
    (`PC1..PC15`, `Sym_PC1..Sym_PC15`), and adaptive calibration metadata
    (`efd_harmonic_opt`, `efd_k95_pcs`) into `leaf_efa_harmonics`.
    """
    if demo or not efa_rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "leaf_efa_harmonics")
    eff_run_id = str(run_id or os.environ.get("RUN_ID") or "canonical").strip()
    ensure_pipeline_run(conn, run_id=eff_run_id)

    norm_rows: List[Dict[str, Any]] = []
    for idx, r in enumerate(efa_rows, start=1):
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        lid = _to_int(r.get("leaf_id"), default=idx) or idx
        d: Dict[str, Any] = {
            "run_id": eff_run_id,
            "catalogNumber": cat,
            "plant_individual_id": pid,
            "leaf_id": lid,
        }
        for col in _LEAF_EFA_STR_COLS:
            if col in r:
                d[col] = _to_str(r.get(col))
        for col in _LEAF_EFA_BOOL_COLS:
            if col in r:
                d[col] = _to_bool_int(r.get(col))
        for col in _LEAF_EFA_INT_COLS:
            if col in r:
                d[col] = _to_int(r.get(col))
        for col in _LEAF_EFA_FLOAT_COLS:
            if col in r:
                d[col] = _to_float(r.get(col))
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="leaf_efa_harmonics",
        rows=norm_rows,
        pk_cols=("run_id", "catalogNumber", "plant_individual_id", "leaf_id"),
        preserve_existing_on_null=False,
    )


_OTU_MORPH_BOOL_COLS: Tuple[str, ...] = (
    "is_promoted_anchor",
    "is_passive",
    "is_calibrated_mm",
    "is_dissected",
    "is_contour_outlier",
    "morphometric_outlier_flag",
    "misidentification_flag",
    "is_type_specimen",
    "has_empirical_type_anchor",
    "is_typical_morphotype",
)
_OTU_MORPH_INT_COLS: Tuple[str, ...] = (
    "leaf_count",
    "n_tier1_leaves",
    "n_tier2_leaves",
    "efd_harmonic_opt",
    "efd_k95_pcs",
    "basal_leaf_count",
)
_OTU_MORPH_STR_COLS: Tuple[str, ...] = (
    "assigned_tier",
    "species_standardized",
    "determiner_tier",
    "type_status",
    "reproductive_trait_status",
    "morph_tier_source",
    "gmm_cluster",
    "cda_predicted_taxon",
    "cda_second_taxon",
    "ecoregion_cda_predicted",
    "triage_priority",
    "discordance_reason",
    "nearest_type_morph_taxon",
)
_OTU_MORPH_FLOAT_COLS: Tuple[str, ...] = (
    "authority_weight",
    "pixels_per_mm",
    "solidity",
    "dissection_index",
    "q25_solidity",
    "dissected_leaf_fraction",
    "ucs_score",
    "aspect_ratio",
    "area_px",
    "petiole_length_mm",
    "blade_length_mm",
    "petiole_to_blade_ratio",
    "centroid_size",
    "log_centsize",
    "AR_inv",
    "AR_inv_std",
    *HARMONIC_COLS,
    "sym_var",
    "asym_var",
    "pct_asymmetry",
    "harmonic_roughness",
    "leaf_quality_weight",
    "foliar_variance",
    "heteroblasty_index",
    *PC_COLS,
    *SYM_PC_COLS,
    "gmm_uncertainty",
    "cda_posterior_prob",
    "cda_second_posterior_prob",
    "cda_posterior_margin",
    "loocv_accuracy",
    "can1",
    "can2",
    "can3",
    "ecoregion_cda_posterior_prob",
    "ecoregion_can1",
    "ecoregion_can2",
    "type_morph_mahal_d2",
    "type_morph_concordance",
)


def persist_otu_morphometrics(
    conn: sqlite3.Connection,
    run_id: str,
    otu_rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    preserve_existing_on_null: bool = True,
) -> int:
    """
    Persist Step 03/04 specimen-level OTU medians (`A1..D24`, `PC1..PC15`,
    `Sym_PC1..Sym_PC15`), adaptive EFD/feature metadata (`efd_harmonic_opt`,
    `efd_k95_pcs`, `selected_morph_features`), GMM clusters, and CDA scores
    into `otu_morphometrics`.
    """
    if demo or not otu_rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "otu_morphometrics")
    eff_run_id = str(run_id or os.environ.get("RUN_ID") or "canonical").strip()
    ensure_pipeline_run(conn, run_id=eff_run_id)

    norm_rows: List[Dict[str, Any]] = []
    for r in otu_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        d: Dict[str, Any] = {
            "run_id": eff_run_id,
            "catalogNumber": cat,
            "plant_individual_id": pid,
        }
        for col in _OTU_MORPH_STR_COLS:
            if col in r:
                d[col] = _to_str(r.get(col))
        if "selected_morph_features" in r:
            d["selected_morph_features"] = _normalize_semicolon_list_str(
                r.get("selected_morph_features")
            )
        for col in _OTU_MORPH_BOOL_COLS:
            if col in r:
                d[col] = _to_bool_int(r.get(col))
        for col in _OTU_MORPH_INT_COLS:
            if col in r:
                d[col] = _to_int(r.get(col))
        for col in _OTU_MORPH_FLOAT_COLS:
            if col in r:
                d[col] = _to_float(r.get(col))
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="otu_morphometrics",
        rows=norm_rows,
        pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
        preserve_existing_on_null=preserve_existing_on_null,
    )


_SRF_BOOL_COLS: Tuple[str, ...] = (
    "substrate_specialist_flag",
    "piedmont_outcrop_specialist_flag",
    "multimodal_conflict_flag",
    "matches_type_lithology",
)
_SRF_STR_COLS: Tuple[str, ...] = (
    "rf_predicted_taxon",
    "edaphic_best_fit_taxon",
    "nearest_type_niche_taxon",
    "spatial_block_id",
    "pheno_timing_category",
    "anthesis_filter_status",
    "multimodal_concordance",
    "triage_category",
    "conflict_discordance_notes",
)
_SRF_LIST_STR_COLS: Tuple[str, ...] = (
    "primary_env_predictors",
    "primary_clim_predictors",
    "primary_edaph_predictors",
)
_SRF_FLOAT_COLS: Tuple[str, ...] = (
    "rf_predicted_prob",
    "rf_given_taxon_prob",
    "spatial_block_cv_acc",
    "spatial_autocorr_range_km",
    "doy",
    "expected_doy",
    "hopkins_expected",
    "delta_doy",
    "hopkins_adjusted_doy",
    "edaphic_z_score",
    "climatic_z_score",
    "type_edaphic_z_score",
    *WC_BIO_COLUMNS,
    "bio18_precip_warmest_qt",
    "bio19_precip_coldest_qt",
)


def persist_spatial_rf_predictions(
    conn: sqlite3.Connection,
    run_id: str,
    srf_rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
    also_upsert_env_cache: bool = False,
) -> int:
    """
    Persist Step 06 Spatial Random Forest predictions, Hopkins phenology residuals,
    all 19 WorldClim bioclimatic variables (with `_qtr <-> _qt` coalescing),
    `climatic_z_score`, and `primary_*_predictors` into `spatial_rf_predictions`.
    """
    if demo or not srf_rows:
        return 0
    _ensure_dynamic_table_migrated(conn, "spatial_rf_predictions")
    if also_upsert_env_cache:
        upsert_environmental_cache_records(conn, srf_rows, demo=demo)

    eff_run_id = str(run_id or os.environ.get("RUN_ID") or "canonical").strip()
    ensure_pipeline_run(conn, run_id=eff_run_id)

    norm_rows: List[Dict[str, Any]] = []
    for r in srf_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        lat = _to_float(r.get("latitude") if "latitude" in r else r.get("decimalLatitude"))
        lon = _to_float(r.get("longitude") if "longitude" in r else r.get("decimalLongitude"))
        ll_key = _to_str(r.get("lat_lon_key")) or format_lat_lon_key(lat, lon)
        d: Dict[str, Any] = {
            "run_id": eff_run_id,
            "catalogNumber": cat,
            "plant_individual_id": pid,
        }
        if ll_key is not None or "lat_lon_key" in r:
            d["lat_lon_key"] = ll_key
        for col in _SRF_STR_COLS:
            if col in r:
                d[col] = _to_str(r.get(col))
        for col in _SRF_LIST_STR_COLS:
            if col in r:
                d[col] = _normalize_semicolon_list_str(r.get(col))
        for col in _SRF_BOOL_COLS:
            if col in r:
                d[col] = _to_bool_int(r.get(col))
        for col in _SRF_FLOAT_COLS:
            if col in r:
                d[col] = _to_float(r.get(col))
        _coalesce_wc_bio_fields(r, d)
        norm_rows.append(d)

    return upsert_rows_generic(
        conn,
        table_name="spatial_rf_predictions",
        rows=norm_rows,
        pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
        preserve_existing_on_null=False,
    )


# =============================================================================
# Bidirectional Bootstrap Migration (`ingest_existing_csvs_to_db`)
# =============================================================================

def _read_csv_dicts(csv_path: Path) -> List[Dict[str, str]]:
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        return []
    with open(csv_path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f))


def export_table_to_csv(
    conn: sqlite3.Connection,
    table_name: str,
    csv_path: Union[str, Path],
    where_clause: str = "",
    params: Sequence[Any] = (),
) -> int:
    """
    Atomically export rows of a SQLite table or view to a CSV snapshot file.
    Returns the number of exported rows.
    """
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]*$", table_name):
        raise ValueError(f"Invalid SQLite table/view identifier: {table_name!r}")
    dest = Path(csv_path).resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    sql = f'SELECT * FROM "{table_name}"'
    if where_clause.strip():
        sql += f" WHERE {where_clause.strip()}"
    cur = conn.execute(sql, tuple(params))
    col_names = [desc[0] for desc in (cur.description or [])]
    rows = cur.fetchall()
    tmp_path = dest.parent / f".{dest.name}.tmp.{os.getpid()}"
    try:
        with open(tmp_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            if col_names:
                writer.writerow(col_names)
            for row in rows:
                writer.writerow(["" if val is None else val for val in row])
        tmp_path.replace(dest)
    except Exception:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        raise
    return len(rows)


def ingest_existing_run_tables_into_db(
    conn: sqlite3.Connection,
    tables_dir: Union[str, Path],
    run_id: str = "canonical",
    demo: bool = False,
) -> Dict[str, int]:
    """
    Ingest Step 03-07 run-scoped CSV outputs (`leaf_efa_harmonics_individual.csv`,
    `leaf_efa_harmonics.csv`, `morphometric_flags.csv`,
    `multimodal_conflict_flags.csv`, `triage_queue.csv` / `final_curation_triage.csv`)
    into Tier B SQLite tables with full dynamic `A1..D24`, `PC1..PC15`,
    `Sym_PC1..Sym_PC15`, and 19-variable `bio1..bio19` (`_qtr <-> _qt`) support.
    """
    if demo:
        return {}
    t_dir = Path(tables_dir).resolve()
    eff_run = str(run_id).strip() or "canonical"
    ensure_pipeline_run(conn, run_id=eff_run, run_dir=t_dir.parent)
    run_counts: Dict[str, int] = {}

    ind_efa = _read_csv_dicts(t_dir / "leaf_efa_harmonics_individual.csv")
    if ind_efa:
        run_counts["leaf_efa_harmonics"] = persist_leaf_efa_harmonics(
            conn, eff_run, ind_efa, demo=False
        )

    step03_otu = _read_csv_dicts(t_dir / "leaf_efa_harmonics.csv")
    step04_otu = _read_csv_dicts(t_dir / "morphometric_flags.csv")
    n_otu = 0
    if step03_otu:
        n_otu = persist_otu_morphometrics(
            conn, eff_run, step03_otu, demo=False, preserve_existing_on_null=True
        )
    if step04_otu:
        n_otu = persist_otu_morphometrics(
            conn, eff_run, step04_otu, demo=False, preserve_existing_on_null=True
        )
    if step03_otu or step04_otu:
        run_counts["otu_morphometrics"] = n_otu

    srf_rows = _read_csv_dicts(t_dir / "multimodal_conflict_flags.csv")
    if srf_rows:
        run_counts["spatial_rf_predictions"] = persist_spatial_rf_predictions(
            conn, eff_run, srf_rows, demo=False, also_upsert_env_cache=False
        )

    tri_path = t_dir / "triage_queue.csv"
    if not tri_path.is_file():
        tri_path = t_dir.parent / "outputs" / "reports" / "final_curation_triage.csv"
    tri_rows = _read_csv_dicts(tri_path)
    if tri_rows:
        norm_tri: List[Dict[str, Any]] = []
        for r in tri_rows:
            cat = _to_str(r.get("catalogNumber"))
            if not cat:
                continue
            d_tri: Dict[str, Any] = dict(r)
            d_tri["run_id"] = eff_run
            d_tri["catalogNumber"] = cat
            d_tri["plant_individual_id"] = _to_int(r.get("plant_individual_id"), default=0) or 0
            if "is_completed_voucher" in d_tri:
                d_tri["is_completed_voucher"] = _to_bool_int(d_tri.get("is_completed_voucher"))
            norm_tri.append(d_tri)
        run_counts["curation_triage"] = upsert_rows_generic(
            conn, "curation_triage", norm_tri, ("run_id", "catalogNumber", "plant_individual_id")
        )

    return run_counts


def ingest_existing_csvs_to_db(
    workspace_root: Union[str, Path] = PROJECT_ROOT,
    run_id: str = "canonical",
    db_path: Optional[Union[str, Path]] = None,
    tables_dir: Optional[Union[str, Path]] = None,
    contours_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, int]:
    """
    Bootstrap or synchronize `data/packera_pipeline.sqlite` from existing CSV manifests,
    120-vertex contour CSVs (`data/contours/*.csv`), and `fine_grain_edaphic_cache.csv`.
    """
    ws = Path(workspace_root).resolve()
    t_dir = Path(tables_dir).resolve() if tables_dir is not None else (ws / "data" / "tables")
    c_dir = Path(contours_dir).resolve() if contours_dir is not None else (ws / "data" / "contours")
    env_csv = ws / "data" / "environmental" / "fine_grain_edaphic_cache.csv"

    conn = get_db_connection(db_path=db_path, read_only=False, demo=False, init_schema=True, workspace_root=ws)
    counts: Dict[str, int] = {}

    try:
        # 1. Seed taxon_synonyms from harvester if importable
        try:
            from scripts.core.harvester import _DYNAMIC_SYNONYM_TO_CANONICAL, _DYNAMIC_TAXON_AUTHOR_STRINGS

            if _DYNAMIC_SYNONYM_TO_CANONICAL:
                counts["taxon_synonyms"] = upsert_taxon_synonyms(
                    conn, _DYNAMIC_SYNONYM_TO_CANONICAL, _DYNAMIC_TAXON_AUTHOR_STRINGS
                )
        except Exception:
            counts["taxon_synonyms"] = 0

        # 2. Vouchers (`curated_vouchers.csv`)
        vouchers_rows = _read_csv_dicts(t_dir / "curated_vouchers.csv")
        counts["vouchers"] = upsert_vouchers_records(conn, vouchers_rows)

        # 3. Step 02 tables + 120-vertex contours
        failed_rows = _read_csv_dicts(t_dir / "failed_qc_vouchers.csv")
        manifest_path = t_dir / "extracted_leaf_manifest.csv"
        if not manifest_path.is_file():
            manifest_path = t_dir / "extracted_leaves_manifest.csv"
        leaf_manifest_rows = _read_csv_dicts(manifest_path)
        repro_rows = _read_csv_dicts(t_dir / "voucher_reproductive_metrics.csv")
        pheno_rows = _read_csv_dicts(t_dir / "voucher_phenological_states.csv")

        repro_by_cat = {str(r.get("catalogNumber", "")).strip(): r for r in repro_rows if r.get("catalogNumber")}
        pheno_by_cat = {str(r.get("catalogNumber", "")).strip(): r for r in pheno_rows if r.get("catalogNumber")}
        failed_by_cat = {str(r.get("catalogNumber", "")).strip(): r for r in failed_rows if r.get("catalogNumber")}
        leaves_by_cat: Dict[str, List[Dict[str, Any]]] = {}
        for lr in leaf_manifest_rows:
            cat = str(lr.get("catalogNumber", "")).strip()
            if cat:
                leaves_by_cat.setdefault(cat, []).append(lr)

        all_step02_cats = sorted(
            set(leaves_by_cat.keys())
            | set(failed_by_cat.keys())
            | set(repro_by_cat.keys())
            | set(pheno_by_cat.keys())
        )

        # Index contour files by filename for O(1) resolution
        contour_index: Dict[str, Path] = {}
        for scan_cdir in (c_dir, ws / "data" / "contours"):
            if scan_cdir.is_dir():
                with os.scandir(scan_cdir) as it:
                    for entry in it:
                        if entry.is_file() and entry.name.endswith(".csv"):
                            contour_index.setdefault(entry.name, Path(entry.path))

        migrated_leaves = 0
        for cat in all_step02_cats:
            cat_leaves = leaves_by_cat.get(cat, [])
            coords_by_leaf: Dict[int, List[List[float]]] = {}
            valid_cat_leaves: List[Dict[str, Any]] = []
            for lr in cat_leaves:
                lid = _to_int(lr.get("leaf_id"), default=1) or 1
                raw_cp = _to_str(lr.get("contour_path"))
                c_file: Optional[Path] = None
                if raw_cp:
                    cand = Path(raw_cp)
                    if cand.is_file():
                        c_file = cand
                    elif cand.name in contour_index:
                        c_file = contour_index[cand.name]
                if c_file is None:
                    fallback_name = f"{cat}_leaf{lid}.csv"
                    c_file = contour_index.get(fallback_name)
                if c_file is not None:
                    pts = read_contour_csv_xy(c_file)
                    if pts:
                        coords_by_leaf[lid] = pts
                        valid_cat_leaves.append(lr)

            fail_rec = failed_by_cat.get(cat)
            if not valid_cat_leaves and fail_rec is None and not cat_leaves:
                # Voucher only appeared in reproductive/phenology table; ensure stub individual without marking FAILED_QC
                with conn:
                    conn.execute(
                        "INSERT OR IGNORE INTO vouchers (catalogNumber, qc_status) VALUES (?, 'PENDING');",
                        (cat,),
                    )
                    conn.execute(
                        "INSERT OR IGNORE INTO plant_individuals (catalogNumber, plant_individual_id) VALUES (?, 0);",
                        (cat,),
                    )
                if cat in repro_by_cat or cat in pheno_by_cat:
                    rep = repro_by_cat.get(cat, {})
                    phe = pheno_by_cat.get(cat, {})
                    ar_val = _to_float(rep.get("AR_inv")) or _to_float(rep.get("capitulum_aspect_ratio"))
                    upsert_rows_generic(
                        conn,
                        "voucher_reproductive_phenology",
                        [{
                            "catalogNumber": cat,
                            "plant_individual_id": 0,
                            "capitula_count": _to_int(rep.get("capitula_count")),
                            "involucre_height_px": _to_float(rep.get("involucre_height_px")),
                            "involucre_width_px": _to_float(rep.get("involucre_width_px")),
                            "involucre_height_mm": _to_float(rep.get("involucre_height_mm")),
                            "involucre_width_mm": _to_float(rep.get("involucre_width_mm")),
                            "AR_inv": ar_val,
                            "capitulum_aspect_ratio": ar_val,
                            "total_capitula": _to_int(phe.get("total_capitula")),
                            "n_anthesis": _to_int(phe.get("n_anthesis")),
                            "n_bud": _to_int(phe.get("n_bud")),
                            "n_fruiting": _to_int(phe.get("n_fruit")),
                            "n_fruit": _to_int(phe.get("n_fruit")),
                            "voucher_phenological_state": _to_str(phe.get("voucher_phenological_state")),
                            "dominant_pappus_ratio": _to_float(phe.get("dominant_pappus_ratio")),
                        }],
                        pk_cols=("catalogNumber", "plant_individual_id"),
                    )
                continue

            record_voucher_segmentation_transaction(
                conn=conn,
                catalog_number=cat,
                extracted_records=valid_cat_leaves,
                failure_record=fail_rec,
                repro_record=repro_by_cat.get(cat),
                pheno_record=pheno_by_cat.get(cat),
                contour_coords_by_leaf=coords_by_leaf,
                output_dir=t_dir.parent,
                workspace_root=ws,
            )
            migrated_leaves += len(valid_cat_leaves)

        # Ensure every voucher in `vouchers` has at least one row in `plant_individuals`
        with conn:
            conn.execute(
                """
                INSERT OR IGNORE INTO plant_individuals (catalogNumber, plant_individual_id)
                SELECT catalogNumber, 0 FROM vouchers;
                """
            )

        counts["extracted_leaves"] = migrated_leaves
        cur_pi = conn.execute("SELECT COUNT(*) FROM plant_individuals;")
        counts["plant_individuals"] = int((cur_pi.fetchone() or [0])[0])

        # 4. Environmental cache (`data/environmental/fine_grain_edaphic_cache.csv`)
        if env_csv.is_file():
            env_rows = _read_csv_dicts(env_csv)
            counts["environmental_cache"] = upsert_environmental_cache_records(conn, env_rows, demo=False)
        else:
            counts["environmental_cache"] = 0

        # 5. Register active run & cohort (`run_id`)
        eff_run = str(run_id).strip() or "canonical"
        ensure_pipeline_run(conn, run_id=eff_run, run_dir=t_dir.parent)
        if vouchers_rows:
            counts["run_active_cohort"] = populate_run_active_cohort(
                conn, run_id=eff_run, voucher_rows=vouchers_rows, clear_existing=True
            )
        else:
            cur_all_v = conn.execute("SELECT catalogNumber, species_standardized, determiner_tier FROM vouchers;")
            v_stub = [dict(r) for r in cur_all_v.fetchall()]
            counts["run_active_cohort"] = populate_run_active_cohort(
                conn, run_id=eff_run, voucher_rows=v_stub, clear_existing=True
            )

        # 6. Optional downstream run-scoped tables if present in `t_dir`
        counts.update(
            ingest_existing_run_tables_into_db(
                conn=conn, tables_dir=t_dir, run_id=eff_run, demo=False
            )
        )

        logger.info("SQLite migration summary (%s): %s", eff_run, counts)
        return counts
    finally:
        conn.close()


def purge_segments_and_runs(
    workspace_root: Optional[Union[str, Path]] = None,
    db_path: Optional[Union[str, Path]] = None,
    output_dir: Optional[Union[str, Path]] = None,
    purge_runs: bool = True,
    vacuum: bool = True,
) -> Dict[str, int]:
    """
    Purges all Step 02 to 07 segmentation artifacts (masks, contours, ECT matrices,
    cropped patches, non-Step-01 tables, figures, reports, and run directories) and
    clears Step 02 to 07 SQLite tables while preserving Step 01 harvested vouchers,
    raw voucher imagery, taxon synonyms, and the environmental cache.
    """
    import shutil as _shutil

    ws = Path(workspace_root).resolve() if workspace_root is not None else PROJECT_ROOT
    data_root = Path(output_dir).resolve() if output_dir is not None else (ws / "data")

    target_dirs = [
        data_root / "contours",
        data_root / "masks",
        data_root / "ect",
        data_root / "cropped_patches",
    ]
    if data_root == (ws / "data").resolve():
        target_dirs.extend([
            ws / "outputs" / "figures",
            ws / "outputs" / "reports",
        ])
        if purge_runs:
            target_dirs.append(ws / "runs")

    removed_files = 0
    for d in target_dirs:
        if d.exists() or d.is_symlink():
            if d.is_symlink() or d.is_file():
                try:
                    d.unlink(missing_ok=True)
                    removed_files += 1
                except OSError:
                    pass
            else:
                removed_files += sum(1 for p in d.rglob("*") if p.is_file() or p.is_symlink())
                _shutil.rmtree(d, ignore_errors=True)
        d.mkdir(parents=True, exist_ok=True)

    for sub in ("tier1_intact", "tier1_pristine", "tier2_reflected"):
        (data_root / "masks" / sub).mkdir(parents=True, exist_ok=True)
    for sub in ("capitula", "rosettes"):
        (data_root / "cropped_patches" / sub).mkdir(parents=True, exist_ok=True)

    tables_dir = data_root / "tables"
    if tables_dir.exists():
        for f in tables_dir.iterdir():
            if f.is_file() and f.name != "curated_vouchers.csv":
                try:
                    f.unlink(missing_ok=True)
                    removed_files += 1
                except OSError:
                    pass

    resolved_db = resolve_db_path(db_path=db_path, workspace_root=ws)
    purged_rows = 0
    if resolved_db != ":memory:" and Path(resolved_db).exists():
        conn = get_db_connection(db_path=resolved_db, init_schema=True, workspace_root=ws)
        try:
            conn.execute("PRAGMA foreign_keys = OFF;")
            purge_tables = (
                "curation_triage",
                "spatial_rf_predictions",
                "vision_xai_audit",
                "otu_morphometrics",
                "leaf_efa_harmonics",
                "run_active_cohort",
                "pipeline_runs",
                "voucher_reproductive_phenology",
                "extracted_leaves",
                "plant_individuals",
                "dinov2_feature_cache",
            )
            for tbl in purge_tables:
                try:
                    cur = conn.execute(f"DELETE FROM {tbl};")
                    if cur.rowcount and cur.rowcount > 0:
                        purged_rows += int(cur.rowcount)
                except sqlite3.Error:
                    pass
            try:
                conn.execute(
                    "UPDATE vouchers SET qc_status = NULL, failure_reason = NULL, failure_details = NULL;"
                )
            except sqlite3.Error:
                pass
            conn.commit()
            if vacuum:
                try:
                    conn.execute("VACUUM;")
                    conn.commit()
                except sqlite3.Error:
                    pass
        finally:
            conn.close()

    summary = {"removed_files": int(removed_files), "purged_db_rows": int(purged_rows)}
    logger.info("Purged Step 02-07 segments and run artifacts: %s", summary)
    return summary
