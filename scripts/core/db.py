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
          voucher_reproductive_phenology, dinov2_feature_cache,
          environmental_cache + environmental_rtree (R*Tree spatial index).
      - Tier B (Run-Scoped Analytical Store, partitioned by run_id):
          pipeline_runs, run_active_cohort, leaf_efa_harmonics,
          otu_morphometrics, vision_xai_audit, spatial_rf_predictions,
          curation_triage, and v_integrative_triage_master.

    ABI Invariant:
      Relies strictly on Python's built-in standard library (`sqlite3`, `csv`,
      `json`, `hashlib`, etc.) at import time so `.venv_LM2` requires zero new
      binary package installations. `duckdb` is imported lazily only when
      explicitly requested inside `.venv`.
===============================================================================
"""

from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
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
    f"{coef}{idx}" for idx in range(1, 13) for coef in ("A", "B", "C", "D")
)
PC_COLS: Tuple[str, ...] = ("PC1", "PC2", "PC3", "PC4", "PC5")
SYM_PC_COLS: Tuple[str, ...] = ("Sym_PC1", "Sym_PC2", "Sym_PC3", "Sym_PC4", "Sym_PC5")

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
    "bio1_temp_mean",
    "bio4_temp_seasonality",
    "bio12_precip_annual",
    "bio15_precip_seasonality",
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


def get_db_connection(
    db_path: Optional[Union[str, Path]] = None,
    read_only: bool = False,
    demo: bool = False,
    init_schema: bool = False,
    workspace_root: Path = PROJECT_ROOT,
) -> sqlite3.Connection:
    """
    Open a SQLite connection configured with WAL mode, NORMAL synchronous durability,
    strict foreign-key enforcement, and a 30-second busy timeout.
    """
    resolved = resolve_db_path(db_path=db_path, workspace_root=workspace_root, demo=demo)

    if read_only and resolved != ":memory:":
        p = Path(resolved)
        if not p.exists():
            raise FileNotFoundError(f"SQLite database not found for read-only connection: {p}")
        conn = sqlite3.connect(f"file:{p.as_posix()}?mode=ro", uri=True, timeout=30.0)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA busy_timeout = 30000;")
        conn.execute("PRAGMA foreign_keys = ON;")
        return conn

    if resolved != ":memory:":
        Path(resolved).parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(resolved), timeout=30.0)
    else:
        conn = sqlite3.connect(":memory:", timeout=30.0)

    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode = WAL;")
    conn.execute("PRAGMA synchronous = NORMAL;")
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.execute("PRAGMA busy_timeout = 30000;")

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
    import duckdb  # Lazy import so .venv_LM2 never requires duckdb

    resolved = resolve_db_path(db_path=db_path, workspace_root=workspace_root, demo=False)
    dconn = duckdb.connect(database=":memory:")
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

def _build_schema_ddl() -> str:
    harm_ddl = ",\n        ".join(f"{h} REAL" for h in HARMONIC_COLS)
    pc_ddl = ",\n        ".join(f"{p} REAL" for p in PC_COLS)
    sym_pc_ddl = ",\n        ".join(f"{sp} REAL" for sp in SYM_PC_COLS)

    env_col_defs: List[str] = ["lat_lon_key TEXT PRIMARY KEY"]
    for col in ENVIRONMENTAL_CACHE_COLUMNS:
        if col == "lat_lon_key":
            continue
        if col in ENVIRONMENTAL_TEXT_COLS:
            env_col_defs.append(f"{col} TEXT")
        elif col in ENVIRONMENTAL_BOOL_COLS:
            env_col_defs.append(f"{col} INTEGER")
        else:
            env_col_defs.append(f"{col} REAL")
    env_ddl = ",\n        ".join(env_col_defs)

    env_view_selects = ",\n        ".join(
        f"ec.{col} AS {col}"
        for col in ENVIRONMENTAL_CACHE_COLUMNS
        if col not in {"lat_lon_key", "catalogNumber", "latitude", "longitude"}
    )
    harm_view_selects = ",\n        ".join(f"om.{h} AS {h}" for h in HARMONIC_COLS)

    return f"""
    -- =========================================================================
    -- TIER A: CANONICAL SHARED STORE
    -- =========================================================================

    CREATE TABLE IF NOT EXISTS taxon_synonyms (
        raw_name_lower TEXT PRIMARY KEY,
        canonical_packera TEXT NOT NULL,
        protologue_authors_json TEXT,
        updated_at TEXT
    );

    CREATE TABLE IF NOT EXISTS vouchers (
        catalogNumber TEXT PRIMARY KEY,
        institutionCode TEXT,
        scientificName TEXT,
        species_raw TEXT,
        recordedBy TEXT,
        recordNumber TEXT,
        identifiedBy TEXT,
        determiner_raw TEXT,
        determiner_tier TEXT,
        type_status TEXT,
        county TEXT,
        stateProvince TEXT,
        locality TEXT,
        decimalLatitude REAL,
        latitude REAL,
        decimalLongitude REAL,
        longitude REAL,
        coordinateUncertainty REAL,
        year INTEGER,
        month INTEGER,
        day INTEGER,
        eventDate TEXT,
        doy REAL,
        pheno_sin REAL,
        pheno_cos REAL,
        regional_group TEXT,
        image_path TEXT,
        exsiccatae_key TEXT,
        is_primary_duplicate INTEGER,
        duplicate_count INTEGER,
        exsiccatae_taxonomic_conflict INTEGER DEFAULT 0,
        exsiccatae_gold_taxon TEXT,
        species_standardized TEXT,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        qc_status TEXT DEFAULT 'PENDING' CHECK(qc_status IN ('PASS', 'FAILED_QC', 'PENDING')),
        failure_reason TEXT,
        failure_details TEXT
    );

    CREATE INDEX IF NOT EXISTS idx_vouchers_species_std ON vouchers(species_standardized);
    CREATE INDEX IF NOT EXISTS idx_vouchers_qc_status ON vouchers(qc_status);
    CREATE INDEX IF NOT EXISTS idx_vouchers_det_tier ON vouchers(determiner_tier);

    CREATE TABLE IF NOT EXISTS plant_individuals (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        rosette_centroid_x REAL,
        rosette_centroid_y REAL,
        crown_x REAL,
        crown_y REAL,
        plant_height_px REAL,
        plant_height_mm REAL,
        plant_height_cm REAL,
        is_bent_stem INTEGER,
        petiole_length_px REAL,
        petiole_length_mm REAL,
        blade_length_px REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        scape_to_basal_leaf_ratio REAL,
        rosette_patch_path TEXT,
        rosette_mask_path TEXT,
        PRIMARY KEY (catalogNumber, plant_individual_id),
        FOREIGN KEY (catalogNumber) REFERENCES vouchers(catalogNumber) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS extracted_leaves (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        leaf_id INTEGER NOT NULL,
        shape_id TEXT UNIQUE,
        tier TEXT,
        assigned_tier TEXT NOT NULL CHECK(assigned_tier IN ('Tier_1_Direct', 'Tier_2_Reflected')),
        tier_assignment TEXT,
        ucs REAL,
        ucs_score REAL,
        solidity REAL,
        circularity REAL,
        extent REAL,
        aspect_ratio REAL,
        area_px REAL,
        blade_area_cm2 REAL,
        is_folded INTEGER,
        is_dissected INTEGER,
        reflection_applied INTEGER,
        midrib_angle_deg REAL,
        petiole_length_px REAL,
        petiole_length_mm REAL,
        blade_length_px REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        mask_path TEXT,
        crop_path TEXT,
        contour_path TEXT,
        contour_n_pts INTEGER DEFAULT 120,
        contour_xy_json TEXT NOT NULL,
        PRIMARY KEY (catalogNumber, plant_individual_id, leaf_id),
        FOREIGN KEY (catalogNumber, plant_individual_id)
            REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
    );

    CREATE INDEX IF NOT EXISTS idx_extracted_leaves_tier ON extracted_leaves(assigned_tier);
    CREATE INDEX IF NOT EXISTS idx_extracted_leaves_cat ON extracted_leaves(catalogNumber);

    CREATE TABLE IF NOT EXISTS voucher_reproductive_phenology (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        capitula_count INTEGER,
        involucre_height_px REAL,
        involucre_width_px REAL,
        involucre_height_mm REAL,
        involucre_width_mm REAL,
        AR_inv REAL,
        capitulum_aspect_ratio REAL,
        total_capitula INTEGER,
        n_anthesis INTEGER,
        n_bud INTEGER,
        n_fruiting INTEGER,
        n_fruit INTEGER,
        n_senescent INTEGER,
        voucher_phenological_state TEXT,
        is_verified_anthesis INTEGER,
        reproductive_stage_confidence REAL,
        dominant_pappus_ratio REAL,
        PRIMARY KEY (catalogNumber, plant_individual_id),
        FOREIGN KEY (catalogNumber, plant_individual_id)
            REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS dinov2_feature_cache (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        model_name TEXT NOT NULL,
        sigma REAL NOT NULL,
        patch_sha256 TEXT,
        embedding_json TEXT NOT NULL,
        foreground_attention_ratio REAL,
        is_artifact_driven INTEGER,
        surface_orientation TEXT,
        mean_luminance REAL,
        mean_saturation REAL,
        high_brightness_ratio REAL,
        PRIMARY KEY (catalogNumber, plant_individual_id, model_name, sigma)
    );

    CREATE TABLE IF NOT EXISTS environmental_cache (
        {env_ddl}
    );

    CREATE INDEX IF NOT EXISTS idx_environmental_cache_catalog
        ON environmental_cache(catalogNumber);

    CREATE VIRTUAL TABLE IF NOT EXISTS environmental_rtree USING rtree(
        id,
        min_lon, max_lon,
        min_lat, max_lat
    );

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_insert
    AFTER INSERT ON environmental_cache
    WHEN NEW.latitude IS NOT NULL AND NEW.longitude IS NOT NULL
    BEGIN
        INSERT OR REPLACE INTO environmental_rtree(id, min_lon, max_lon, min_lat, max_lat)
        VALUES (NEW.rowid, NEW.longitude, NEW.longitude, NEW.latitude, NEW.latitude);
    END;

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_update
    AFTER UPDATE ON environmental_cache
    BEGIN
        DELETE FROM environmental_rtree WHERE id = OLD.rowid;
        INSERT INTO environmental_rtree(id, min_lon, max_lon, min_lat, max_lat)
        SELECT NEW.rowid, NEW.longitude, NEW.longitude, NEW.latitude, NEW.latitude
        WHERE NEW.latitude IS NOT NULL AND NEW.longitude IS NOT NULL;
    END;

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_delete
    AFTER DELETE ON environmental_cache
    BEGIN
        DELETE FROM environmental_rtree WHERE id = OLD.rowid;
    END;

    -- =========================================================================
    -- TIER B: RUN-SCOPED ANALYTICAL TABLES & MASTER VIEW
    -- =========================================================================

    CREATE TABLE IF NOT EXISTS pipeline_runs (
        run_id TEXT PRIMARY KEY,
        run_dir TEXT,
        start_timestamp TEXT,
        end_timestamp TEXT,
        git_commit_sha TEXT,
        cli_arguments_json TEXT,
        focal_taxon TEXT,
        intraspecific_mode INTEGER DEFAULT 0
    );

    CREATE TABLE IF NOT EXISTS run_active_cohort (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        species_standardized TEXT,
        determiner_tier TEXT,
        is_promoted_bronze INTEGER DEFAULT 0,
        authority_weight REAL,
        spatial_balance_weight REAL DEFAULT 1.0,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE,
        FOREIGN KEY (catalogNumber) REFERENCES vouchers(catalogNumber) ON DELETE CASCADE
    );

    CREATE INDEX IF NOT EXISTS idx_run_active_cohort_run ON run_active_cohort(run_id);

    CREATE TABLE IF NOT EXISTS leaf_efa_harmonics (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        leaf_id INTEGER NOT NULL,
        assigned_tier TEXT,
        solidity REAL,
        ucs_score REAL,
        aspect_ratio REAL,
        area_px REAL,
        is_dissected INTEGER,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        petiole_length_mm REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        {harm_ddl},
        sym_var REAL,
        asym_var REAL,
        pct_asymmetry REAL,
        harmonic_roughness REAL,
        is_contour_outlier INTEGER,
        leaf_quality_weight REAL,
        centroid_size REAL,
        log_centsize REAL,
        {pc_ddl},
        {sym_pc_ddl},
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id, leaf_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS otu_morphometrics (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        assigned_tier TEXT,
        species_standardized TEXT,
        determiner_tier TEXT,
        is_promoted_anchor INTEGER,
        authority_weight REAL,
        is_passive INTEGER,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        solidity REAL,
        ucs_score REAL,
        aspect_ratio REAL,
        area_px REAL,
        is_dissected INTEGER,
        petiole_length_mm REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        centroid_size REAL,
        log_centsize REAL,
        AR_inv REAL,
        AR_inv_std REAL,
        reproductive_trait_status TEXT,
        {harm_ddl},
        sym_var REAL,
        asym_var REAL,
        pct_asymmetry REAL,
        harmonic_roughness REAL,
        is_contour_outlier INTEGER,
        leaf_quality_weight REAL,
        leaf_count INTEGER,
        n_tier1_leaves INTEGER,
        n_tier2_leaves INTEGER,
        morph_tier_source TEXT,
        foliar_variance REAL,
        heteroblasty_index REAL,
        {pc_ddl},
        {sym_pc_ddl},
        gmm_cluster TEXT,
        gmm_uncertainty REAL,
        cda_predicted_taxon TEXT,
        cda_posterior_prob REAL,
        cda_second_taxon TEXT,
        cda_second_posterior_prob REAL,
        cda_posterior_margin REAL,
        loocv_accuracy REAL,
        can1 REAL,
        can2 REAL,
        can3 REAL,
        ecoregion_cda_predicted TEXT,
        ecoregion_cda_posterior_prob REAL,
        ecoregion_can1 REAL,
        ecoregion_can2 REAL,
        morphometric_outlier_flag INTEGER,
        misidentification_flag INTEGER,
        triage_priority TEXT,
        discordance_reason TEXT,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS vision_xai_audit (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        patch_path TEXT,
        given_label TEXT,
        vision_predicted_label TEXT,
        predicted_label TEXT,
        confidence_given_class REAL,
        confidence_predicted_class REAL,
        vision_second_prob REAL,
        vision_margin REAL,
        label_quality_score REAL,
        c_error REAL,
        is_cleanlab_issue INTEGER,
        is_label_corrupted INTEGER,
        vision_ood_score REAL,
        vision_is_ood INTEGER,
        foreground_attention_ratio REAL,
        is_artifact_driven INTEGER,
        is_background_shortcut INTEGER,
        surface_orientation TEXT,
        mean_luminance REAL,
        mean_saturation REAL,
        high_brightness_ratio REAL,
        tomentum_taxon_discordance INTEGER,
        triage_action TEXT,
        discordance_reason TEXT,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS spatial_rf_predictions (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        lat_lon_key TEXT,
        rf_predicted_taxon TEXT,
        rf_predicted_prob REAL,
        rf_given_taxon_prob REAL,
        edaphic_best_fit_taxon TEXT,
        spatial_block_cv_acc REAL,
        spatial_autocorr_range_km REAL,
        spatial_block_id TEXT,
        doy REAL,
        expected_doy REAL,
        hopkins_expected REAL,
        delta_doy REAL,
        hopkins_adjusted_doy REAL,
        pheno_timing_category TEXT,
        anthesis_filter_status TEXT,
        edaphic_z_score REAL,
        substrate_specialist_flag INTEGER,
        piedmont_outcrop_specialist_flag INTEGER,
        multimodal_concordance TEXT,
        multimodal_conflict_flag INTEGER,
        triage_category TEXT,
        conflict_discordance_notes TEXT,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS curation_triage (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        w_morph_kappa REAL,
        w_vision_kappa REAL,
        w_edaphic_kappa REAL,
        w_pheno_kappa REAL,
        n_active_modalities INTEGER,
        n_available_modalities INTEGER,
        probabilistic_consensus_conflict REAL,
        composite_conflict_score REAL,
        probabilistic_consensus_score REAL,
        hopkins_adjusted_doy REAL,
        delta_doy_resid REAL,
        taxonomic_status_call TEXT,
        recommended_determination TEXT,
        taxonomic_treatment TEXT,
        taxonomic_treatment_broad TEXT,
        taxonomic_treatment_subtype TEXT,
        synthesis_triage_priority TEXT,
        synthesis_triage_action TEXT,
        synthesis_rationale TEXT,
        is_completed_voucher INTEGER,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    DROP VIEW IF EXISTS v_integrative_triage_master;
    CREATE VIEW v_integrative_triage_master AS
    SELECT
        rac.run_id AS run_id,
        rac.catalogNumber AS catalogNumber,
        rac.plant_individual_id AS plant_individual_id,
        v.institutionCode AS institutionCode,
        v.scientificName AS scientificName,
        v.species_raw AS species_raw,
        v.recordedBy AS recordedBy,
        v.recordNumber AS recordNumber,
        v.identifiedBy AS identifiedBy,
        v.determiner_raw AS determiner_raw,
        COALESCE(rac.determiner_tier, om.determiner_tier, v.determiner_tier) AS determiner_tier,
        v.type_status AS type_status,
        v.county AS county,
        v.stateProvince AS stateProvince,
        v.locality AS locality,
        COALESCE(v.decimalLatitude, v.latitude) AS decimalLatitude,
        COALESCE(v.latitude, v.decimalLatitude) AS latitude,
        COALESCE(v.decimalLongitude, v.longitude) AS decimalLongitude,
        COALESCE(v.longitude, v.decimalLongitude) AS longitude,
        v.coordinateUncertainty AS coordinateUncertainty,
        v.year AS year,
        v.month AS month,
        v.day AS day,
        v.eventDate AS eventDate,
        v.regional_group AS regional_group,
        v.image_path AS image_path,
        v.exsiccatae_key AS exsiccatae_key,
        v.is_primary_duplicate AS is_primary_duplicate,
        v.duplicate_count AS duplicate_count,
        COALESCE(rac.species_standardized, om.species_standardized, v.species_standardized) AS species_standardized,
        COALESCE(rac.is_promoted_bronze, om.is_promoted_anchor, 0) AS is_promoted_anchor,
        COALESCE(rac.authority_weight, om.authority_weight) AS authority_weight,
        COALESCE(rac.spatial_balance_weight, 1.0) AS spatial_balance_weight,
        COALESCE(om.pixels_per_mm, v.pixels_per_mm) AS pixels_per_mm,
        COALESCE(om.is_calibrated_mm, v.is_calibrated_mm) AS is_calibrated_mm,
        v.qc_status AS qc_status,
        -- Plant individual & linear petiole/scape traits
        pi.rosette_centroid_x AS rosette_centroid_x,
        pi.rosette_centroid_y AS rosette_centroid_y,
        pi.crown_x AS crown_x,
        pi.crown_y AS crown_y,
        pi.plant_height_px AS plant_height_px,
        pi.plant_height_mm AS plant_height_mm,
        pi.plant_height_cm AS plant_height_cm,
        pi.is_bent_stem AS is_bent_stem,
        pi.petiole_length_px AS petiole_length_px,
        COALESCE(om.petiole_length_mm, pi.petiole_length_mm) AS petiole_length_mm,
        pi.blade_length_px AS blade_length_px,
        COALESCE(om.blade_length_mm, pi.blade_length_mm) AS blade_length_mm,
        COALESCE(om.petiole_to_blade_ratio, pi.petiole_to_blade_ratio) AS petiole_to_blade_ratio,
        pi.scape_to_basal_leaf_ratio AS scape_to_basal_leaf_ratio,
        pi.rosette_patch_path AS rosette_patch_path,
        pi.rosette_mask_path AS rosette_mask_path,
        -- Reproductive & phenological traits (with sheet-level pid=0 fallback)
        COALESCE(vrp.capitula_count, vrp0.capitula_count) AS capitula_count,
        COALESCE(vrp.involucre_height_px, vrp0.involucre_height_px) AS involucre_height_px,
        COALESCE(vrp.involucre_width_px, vrp0.involucre_width_px) AS involucre_width_px,
        COALESCE(vrp.involucre_height_mm, vrp0.involucre_height_mm) AS involucre_height_mm,
        COALESCE(vrp.involucre_width_mm, vrp0.involucre_width_mm) AS involucre_width_mm,
        COALESCE(om.AR_inv, vrp.AR_inv, vrp.capitulum_aspect_ratio, vrp0.AR_inv, vrp0.capitulum_aspect_ratio) AS AR_inv,
        COALESCE(vrp.capitulum_aspect_ratio, vrp.AR_inv, vrp0.capitulum_aspect_ratio, vrp0.AR_inv) AS capitulum_aspect_ratio,
        COALESCE(vrp.total_capitula, vrp0.total_capitula, vrp.capitula_count, vrp0.capitula_count) AS total_capitula,
        COALESCE(vrp.n_anthesis, vrp0.n_anthesis) AS n_anthesis,
        COALESCE(vrp.n_bud, vrp0.n_bud) AS n_bud,
        COALESCE(vrp.n_fruiting, vrp.n_fruit, vrp0.n_fruiting, vrp0.n_fruit) AS n_fruiting,
        COALESCE(vrp.n_fruit, vrp.n_fruiting, vrp0.n_fruit, vrp0.n_fruiting) AS n_fruit,
        COALESCE(vrp.n_senescent, vrp0.n_senescent) AS n_senescent,
        COALESCE(vrp.voucher_phenological_state, vrp0.voucher_phenological_state) AS voucher_phenological_state,
        COALESCE(vrp.is_verified_anthesis, vrp0.is_verified_anthesis) AS is_verified_anthesis,
        COALESCE(vrp.reproductive_stage_confidence, vrp0.reproductive_stage_confidence) AS reproductive_stage_confidence,
        COALESCE(vrp.dominant_pappus_ratio, vrp0.dominant_pappus_ratio) AS dominant_pappus_ratio,
        -- Step 03 & Step 04 OTU Morphometrics & CDA
        om.assigned_tier AS assigned_tier,
        om.solidity AS solidity,
        om.ucs_score AS ucs_score,
        om.aspect_ratio AS aspect_ratio,
        om.area_px AS area_px,
        om.is_dissected AS is_dissected,
        om.centroid_size AS centroid_size,
        om.log_centsize AS log_centsize,
        {harm_view_selects},
        om.sym_var AS sym_var,
        om.asym_var AS asym_var,
        om.pct_asymmetry AS pct_asymmetry,
        om.harmonic_roughness AS harmonic_roughness,
        om.is_contour_outlier AS is_contour_outlier,
        om.leaf_quality_weight AS leaf_quality_weight,
        om.leaf_count AS leaf_count,
        om.n_tier1_leaves AS n_tier1_leaves,
        om.n_tier2_leaves AS n_tier2_leaves,
        om.morph_tier_source AS morph_tier_source,
        om.foliar_variance AS foliar_variance,
        om.heteroblasty_index AS heteroblasty_index,
        om.PC1 AS PC1,
        om.PC2 AS PC2,
        om.PC3 AS PC3,
        om.PC4 AS PC4,
        om.PC5 AS PC5,
        om.Sym_PC1 AS Sym_PC1,
        om.Sym_PC2 AS Sym_PC2,
        om.Sym_PC3 AS Sym_PC3,
        om.Sym_PC4 AS Sym_PC4,
        om.Sym_PC5 AS Sym_PC5,
        om.gmm_cluster AS gmm_cluster,
        om.gmm_uncertainty AS gmm_uncertainty,
        om.cda_predicted_taxon AS cda_predicted_taxon,
        om.cda_posterior_prob AS cda_posterior_prob,
        om.cda_second_taxon AS cda_second_taxon,
        om.cda_second_posterior_prob AS cda_second_posterior_prob,
        om.cda_posterior_margin AS cda_posterior_margin,
        om.loocv_accuracy AS loocv_accuracy,
        om.can1 AS can1,
        om.can2 AS can2,
        om.can3 AS can3,
        om.ecoregion_cda_predicted AS ecoregion_cda_predicted,
        om.ecoregion_cda_posterior_prob AS ecoregion_cda_posterior_prob,
        om.morphometric_outlier_flag AS morphometric_outlier_flag,
        COALESCE(om.misidentification_flag, om.morphometric_outlier_flag) AS misidentification_flag,
        om.triage_priority AS triage_priority,
        om.discordance_reason AS morph_discordance_reason,
        -- Step 05 Vision XAI & Cleanlab Audit
        COALESCE(vx.vision_predicted_label, vx.predicted_label) AS vision_predicted_label,
        COALESCE(vx.predicted_label, vx.vision_predicted_label) AS predicted_label,
        vx.confidence_given_class AS confidence_given_class,
        vx.confidence_predicted_class AS confidence_predicted_class,
        vx.vision_second_prob AS vision_second_prob,
        vx.vision_margin AS vision_margin,
        vx.label_quality_score AS label_quality_score,
        vx.c_error AS c_error,
        vx.is_cleanlab_issue AS is_cleanlab_issue,
        vx.is_label_corrupted AS is_label_corrupted,
        vx.vision_ood_score AS vision_ood_score,
        vx.vision_is_ood AS vision_is_ood,
        vx.foreground_attention_ratio AS foreground_attention_ratio,
        vx.is_artifact_driven AS is_artifact_driven,
        vx.is_background_shortcut AS is_background_shortcut,
        vx.surface_orientation AS surface_orientation,
        vx.mean_luminance AS mean_luminance,
        vx.mean_saturation AS mean_saturation,
        vx.high_brightness_ratio AS high_brightness_ratio,
        vx.tomentum_taxon_discordance AS tomentum_taxon_discordance,
        -- Step 06 Spatial RF & Hopkins Phenology
        COALESCE(srf.lat_lon_key, ec.lat_lon_key) AS lat_lon_key,
        srf.rf_predicted_taxon AS rf_predicted_taxon,
        srf.rf_predicted_prob AS rf_predicted_prob,
        srf.rf_given_taxon_prob AS rf_given_taxon_prob,
        srf.edaphic_best_fit_taxon AS edaphic_best_fit_taxon,
        srf.spatial_block_cv_acc AS spatial_block_cv_acc,
        srf.spatial_autocorr_range_km AS spatial_autocorr_range_km,
        srf.spatial_block_id AS spatial_block_id,
        COALESCE(srf.doy, v.doy) AS doy,
        v.pheno_sin AS pheno_sin,
        v.pheno_cos AS pheno_cos,
        srf.expected_doy AS expected_doy,
        srf.hopkins_expected AS hopkins_expected,
        srf.delta_doy AS delta_doy,
        srf.pheno_timing_category AS pheno_timing_category,
        srf.anthesis_filter_status AS anthesis_filter_status,
        srf.edaphic_z_score AS edaphic_z_score,
        srf.substrate_specialist_flag AS substrate_specialist_flag,
        srf.piedmont_outcrop_specialist_flag AS piedmont_outcrop_specialist_flag,
        srf.multimodal_concordance AS multimodal_concordance,
        srf.multimodal_conflict_flag AS multimodal_conflict_flag,
        srf.triage_category AS triage_category,
        srf.conflict_discordance_notes AS conflict_discordance_notes,
        -- Environmental Cache (Fine-grain SSURGO / POLARIS / SOLUS / 3DEP / SGMC / WorldClim)
        {env_view_selects},
        -- Step 07 Curation Triage
        ct.w_morph_kappa AS w_morph_kappa,
        ct.w_vision_kappa AS w_vision_kappa,
        ct.w_edaphic_kappa AS w_edaphic_kappa,
        ct.w_pheno_kappa AS w_pheno_kappa,
        ct.n_active_modalities AS n_active_modalities,
        ct.n_available_modalities AS n_available_modalities,
        ct.probabilistic_consensus_conflict AS probabilistic_consensus_conflict,
        ct.composite_conflict_score AS composite_conflict_score,
        ct.probabilistic_consensus_score AS probabilistic_consensus_score,
        COALESCE(ct.hopkins_adjusted_doy, srf.hopkins_adjusted_doy) AS hopkins_adjusted_doy,
        ct.delta_doy_resid AS delta_doy_resid,
        ct.taxonomic_status_call AS taxonomic_status_call,
        ct.recommended_determination AS recommended_determination,
        ct.taxonomic_treatment AS taxonomic_treatment,
        ct.taxonomic_treatment_broad AS taxonomic_treatment_broad,
        ct.taxonomic_treatment_subtype AS taxonomic_treatment_subtype,
        ct.synthesis_triage_priority AS synthesis_triage_priority,
        ct.synthesis_triage_action AS synthesis_triage_action,
        ct.synthesis_rationale AS synthesis_rationale,
        ct.is_completed_voucher AS is_completed_voucher
    FROM run_active_cohort rac
    LEFT JOIN vouchers v
        ON v.catalogNumber = rac.catalogNumber
    LEFT JOIN plant_individuals pi
        ON pi.catalogNumber = rac.catalogNumber
       AND pi.plant_individual_id = rac.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp
        ON vrp.catalogNumber = rac.catalogNumber
       AND vrp.plant_individual_id = rac.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp0
        ON vrp0.catalogNumber = rac.catalogNumber
       AND vrp0.plant_individual_id = 0
    LEFT JOIN otu_morphometrics om
        ON om.run_id = rac.run_id
       AND om.catalogNumber = rac.catalogNumber
       AND om.plant_individual_id = rac.plant_individual_id
    LEFT JOIN vision_xai_audit vx
        ON vx.run_id = rac.run_id
       AND vx.catalogNumber = rac.catalogNumber
       AND vx.plant_individual_id = rac.plant_individual_id
    LEFT JOIN spatial_rf_predictions srf
        ON srf.run_id = rac.run_id
       AND srf.catalogNumber = rac.catalogNumber
       AND srf.plant_individual_id = rac.plant_individual_id
    LEFT JOIN environmental_cache ec
        ON ec.lat_lon_key = COALESCE(
            srf.lat_lon_key,
            CASE
                WHEN COALESCE(v.latitude, v.decimalLatitude) IS NOT NULL
                 AND COALESCE(v.longitude, v.decimalLongitude) IS NOT NULL
                THEN printf('%.5f_%.5f', COALESCE(v.latitude, v.decimalLatitude), COALESCE(v.longitude, v.decimalLongitude))
                ELSE NULL
            END,
            (SELECT ec2.lat_lon_key FROM environmental_cache ec2 WHERE ec2.catalogNumber = rac.catalogNumber LIMIT 1)
        )
    LEFT JOIN curation_triage ct
        ON ct.run_id = rac.run_id
       AND ct.catalogNumber = rac.catalogNumber
       AND ct.plant_individual_id = rac.plant_individual_id;
    """


SCHEMA_DDL: str = _build_schema_ddl()

_VOUCHERS_MIGRATION_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("exsiccatae_taxonomic_conflict", "INTEGER DEFAULT 0"),
    ("exsiccatae_gold_taxon", "TEXT"),
)


def _migrate_vouchers_schema_columns(conn: sqlite3.Connection) -> None:
    """Idempotently add optional exsiccatae conflict audit columns to existing `vouchers` tables."""
    existing_cols = {
        str(r[1]) for r in conn.execute("PRAGMA table_info(vouchers);").fetchall()
    }
    if not existing_cols:
        return
    for col_name, col_ddl in _VOUCHERS_MIGRATION_COLUMNS:
        if col_name not in existing_cols:
            conn.execute(f"ALTER TABLE vouchers ADD COLUMN {col_name} {col_ddl};")
            existing_cols.add(col_name)


def init_db_schema(conn: sqlite3.Connection) -> None:
    """Create all Tier A and Tier B tables, indexes, triggers, R*Tree, and master view."""
    _migrate_vouchers_schema_columns(conn)
    conn.executescript(SCHEMA_DDL)
    _migrate_vouchers_schema_columns(conn)
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.commit()


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

    # Discover saved rosette patch / companion mask paths per plant individual
    out_base = Path(output_dir) if output_dir is not None else (workspace_root / "data")
    rosettes_dir = (
        out_base / "cropped_patches" / "rosettes_dense"
        if (out_base / "cropped_patches" / "rosettes_dense").exists()
        else (out_base / "data" / "cropped_patches" / "rosettes_dense")
    )
    basal_dir = (
        out_base / "cropped_patches" / "basal_leaves_raw"
        if (out_base / "cropped_patches" / "basal_leaves_raw").exists()
        else (out_base / "data" / "cropped_patches" / "basal_leaves_raw")
    )

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
        if rosettes_dir.exists():
            cand_p = rosettes_dir / f"{cat}_p{pid}_rosette.jpg"
            cand_p0 = rosettes_dir / f"{cat}_rosette.jpg"
            if cand_p.exists():
                r_patch = str(cand_p)
            elif pid == 0 and cand_p0.exists():
                r_patch = str(cand_p0)
            cand_m = rosettes_dir / f"{cat}_p{pid}_rosette_mask.png"
            if cand_m.exists():
                r_mask = str(cand_m)

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
        if crop_p is None and basal_dir.exists():
            cand_crop = basal_dir / f"{cat}_leaf{leaf_id}.jpg"
            if cand_crop.exists():
                crop_p = str(cand_crop)

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

    with conn:
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
# Step 05: DINOv2 Feature Cache & Vision XAI Audit Persistence
# =============================================================================

def fetch_cached_dinov2_features(
    conn: sqlite3.Connection,
    specimen_keys: Sequence[Tuple[str, int, Optional[str]]],
    model_name: str,
    sigma: float,
) -> Dict[Tuple[str, int], Dict[str, Any]]:
    """
    Fetch cached 768-d DINOv2 embeddings and foreground attention ratios from
    `dinov2_feature_cache` for `(catalogNumber, plant_individual_id, model_name, sigma)`.
    If `patch_sha256` is supplied in `specimen_keys` and differs from a non-null cached
    checksum, treats that specimen as a cache miss so modified patches re-embed cleanly.
    """
    if not specimen_keys:
        return {}
    cats = sorted({str(k[0]).strip() for k in specimen_keys if str(k[0]).strip()})
    if not cats:
        return {}

    expected_sha: Dict[Tuple[str, int], Optional[str]] = {
        (str(cat).strip(), int(pid)): sha for cat, pid, sha in specimen_keys
    }
    hits: Dict[Tuple[str, int], Dict[str, Any]] = {}
    sigma_val = round(float(sigma), 4)

    chunk_size = 500
    for i in range(0, len(cats), chunk_size):
        chunk = cats[i : i + chunk_size]
        ph = ", ".join("?" for _ in chunk)
        cur = conn.execute(
            f"""
            SELECT
                catalogNumber, plant_individual_id, patch_sha256, embedding_json,
                foreground_attention_ratio, is_artifact_driven,
                surface_orientation, mean_luminance, mean_saturation, high_brightness_ratio
            FROM dinov2_feature_cache
            WHERE model_name = ? AND ABS(sigma - ?) < 1e-6
              AND catalogNumber IN ({ph});
            """,
            [str(model_name), sigma_val, *chunk],
        )
        for row in cur.fetchall():
            key = (str(row["catalogNumber"]).strip(), int(row["plant_individual_id"]))
            if key not in expected_sha:
                continue
            exp_hash = expected_sha[key]
            cached_hash = row["patch_sha256"]
            if exp_hash and cached_hash and str(exp_hash) != str(cached_hash):
                continue
            try:
                emb_vec = json.loads(row["embedding_json"])
                if not isinstance(emb_vec, list) or len(emb_vec) == 0:
                    continue
            except (TypeError, ValueError):
                continue
            hits[key] = {
                "embedding": emb_vec,
                "patch_sha256": cached_hash,
                "foreground_attention_ratio": _to_float(row["foreground_attention_ratio"]),
                "is_artifact_driven": bool(row["is_artifact_driven"]) if row["is_artifact_driven"] is not None else False,
                "surface_orientation": row["surface_orientation"],
                "mean_luminance": _to_float(row["mean_luminance"]),
                "mean_saturation": _to_float(row["mean_saturation"]),
                "high_brightness_ratio": _to_float(row["high_brightness_ratio"]),
            }
    return hits


def store_dinov2_features_cache(
    conn: sqlite3.Connection,
    cache_rows: Sequence[Mapping[str, Any]],
    model_name: str,
    sigma: float,
    demo: bool = False,
) -> int:
    """
    Batch-insert newly computed 768-d DINOv2 embeddings and foreground attention
    metrics into `dinov2_feature_cache`. Strictly no-ops when `demo=True`.
    """
    if demo or not cache_rows:
        return 0
    sigma_val = round(float(sigma), 4)
    payload: List[Tuple[Any, ...]] = []
    for r in cache_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        emb = r.get("embedding")
        if emb is None and r.get("embedding_json"):
            emb_json = str(r["embedding_json"])
        elif emb is not None:
            emb_list = [round(float(x), 7) for x in emb]
            emb_json = json.dumps(emb_list, separators=(",", ":"))
        else:
            continue
        payload.append((
            cat,
            pid,
            str(model_name),
            sigma_val,
            _to_str(r.get("patch_sha256")),
            emb_json,
            _to_float(r.get("foreground_attention_ratio")),
            _to_bool_int(r.get("is_artifact_driven"), default=0),
            _to_str(r.get("surface_orientation")),
            _to_float(r.get("mean_luminance")),
            _to_float(r.get("mean_saturation")),
            _to_float(r.get("high_brightness_ratio")),
        ))

    if not payload:
        return 0

    sql = """
    INSERT INTO dinov2_feature_cache (
        catalogNumber, plant_individual_id, model_name, sigma,
        patch_sha256, embedding_json, foreground_attention_ratio, is_artifact_driven,
        surface_orientation, mean_luminance, mean_saturation, high_brightness_ratio
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ON CONFLICT(catalogNumber, plant_individual_id, model_name, sigma) DO UPDATE SET
        patch_sha256 = excluded.patch_sha256,
        embedding_json = excluded.embedding_json,
        foreground_attention_ratio = excluded.foreground_attention_ratio,
        is_artifact_driven = excluded.is_artifact_driven,
        surface_orientation = COALESCE(excluded.surface_orientation, dinov2_feature_cache.surface_orientation),
        mean_luminance = COALESCE(excluded.mean_luminance, dinov2_feature_cache.mean_luminance),
        mean_saturation = COALESCE(excluded.mean_saturation, dinov2_feature_cache.mean_saturation),
        high_brightness_ratio = COALESCE(excluded.high_brightness_ratio, dinov2_feature_cache.high_brightness_ratio);
    """
    with conn:
        conn.executemany(sql, payload)
    return len(payload)


def persist_vision_xai_audit(
    conn: sqlite3.Connection,
    run_id: str,
    audit_rows: Sequence[Mapping[str, Any]],
    demo: bool = False,
) -> int:
    """Persist run-scoped Step 05 Cleanlab & DINOv2 XAI results into `vision_xai_audit`."""
    if demo or not audit_rows:
        return 0
    eff_run_id = str(run_id or os.environ.get("RUN_ID") or "canonical").strip()
    ensure_pipeline_run(conn, run_id=eff_run_id)

    norm_rows: List[Dict[str, Any]] = []
    for r in audit_rows:
        cat = _to_str(r.get("catalogNumber"))
        if not cat:
            continue
        pid = _to_int(r.get("plant_individual_id"), default=0) or 0
        pred_lbl = _to_str(r.get("vision_predicted_label")) or _to_str(r.get("predicted_label"))
        norm_rows.append({
            "run_id": eff_run_id,
            "catalogNumber": cat,
            "plant_individual_id": pid,
            "patch_path": _to_str(r.get("patch_path")),
            "given_label": _to_str(r.get("given_label")),
            "vision_predicted_label": pred_lbl,
            "predicted_label": _to_str(r.get("predicted_label")) or pred_lbl,
            "confidence_given_class": _to_float(r.get("confidence_given_class")),
            "confidence_predicted_class": _to_float(r.get("confidence_predicted_class")),
            "vision_second_prob": _to_float(r.get("vision_second_prob")),
            "vision_margin": _to_float(r.get("vision_margin")),
            "label_quality_score": _to_float(r.get("label_quality_score")),
            "c_error": _to_float(r.get("c_error")),
            "is_cleanlab_issue": _to_bool_int(r.get("is_cleanlab_issue"), default=0),
            "is_label_corrupted": _to_bool_int(r.get("is_label_corrupted"), default=0),
            "vision_ood_score": _to_float(r.get("vision_ood_score")),
            "vision_is_ood": _to_bool_int(r.get("vision_is_ood"), default=0),
            "foreground_attention_ratio": _to_float(r.get("foreground_attention_ratio")),
            "is_artifact_driven": _to_bool_int(r.get("is_artifact_driven"), default=0),
            "is_background_shortcut": _to_bool_int(r.get("is_background_shortcut"), default=0),
            "surface_orientation": _to_str(r.get("surface_orientation")),
            "mean_luminance": _to_float(r.get("mean_luminance")),
            "mean_saturation": _to_float(r.get("mean_saturation")),
            "high_brightness_ratio": _to_float(r.get("high_brightness_ratio")),
            "tomentum_taxon_discordance": _to_bool_int(r.get("tomentum_taxon_discordance"), default=0),
            "triage_action": _to_str(r.get("triage_action")),
            "discordance_reason": _to_str(r.get("discordance_reason")),
        })

    return upsert_rows_generic(
        conn,
        table_name="vision_xai_audit",
        rows=norm_rows,
        pk_cols=("run_id", "catalogNumber", "plant_individual_id"),
        preserve_existing_on_null=False,
    )


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
    if not voucher_rows:
        if clear_existing:
            with conn:
                conn.execute("DELETE FROM run_active_cohort WHERE run_id = ?;", (rid,))
        return 0

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
            cohort_payload.append((rid, cat, pid, sp_std, det_tier, is_prom, auth_wt, sp_wt))

    with conn:
        if clear_existing:
            conn.execute("DELETE FROM run_active_cohort WHERE run_id = ?;", (rid,))
        conn.executemany(
            """
            INSERT INTO run_active_cohort (
                run_id, catalogNumber, plant_individual_id,
                species_standardized, determiner_tier,
                is_promoted_bronze, authority_weight, spatial_balance_weight
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(run_id, catalogNumber, plant_individual_id) DO UPDATE SET
                species_standardized = COALESCE(excluded.species_standardized, run_active_cohort.species_standardized),
                determiner_tier = COALESCE(excluded.determiner_tier, run_active_cohort.determiner_tier),
                is_promoted_bronze = COALESCE(excluded.is_promoted_bronze, run_active_cohort.is_promoted_bronze),
                authority_weight = COALESCE(excluded.authority_weight, run_active_cohort.authority_weight),
                spatial_balance_weight = COALESCE(excluded.spatial_balance_weight, run_active_cohort.spatial_balance_weight);
            """,
            cohort_payload,
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

def upsert_environmental_cache_records(
    conn: sqlite3.Connection,
    records: Sequence[Mapping[str, Any]],
    demo: bool = False,
) -> int:
    """
    UPSERT fine-grain SSURGO / POLARIS / SOLUS100 / 3DEP / SGMC / WorldClim rows into
    `environmental_cache` (automatically updating `environmental_rtree` via triggers).
    Strictly no-ops when `demo=True`.
    """
    if demo or not records:
        return 0
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
    return [dict(r) for r in cur.fetchall()]


# =============================================================================
# Bidirectional Bootstrap Migration (`ingest_existing_csvs_to_db`)
# =============================================================================

def _read_csv_dicts(csv_path: Path) -> List[Dict[str, str]]:
    if not csv_path.is_file() or csv_path.stat().st_size == 0:
        return []
    with open(csv_path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f))


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
        ind_efa = _read_csv_dicts(t_dir / "leaf_efa_harmonics_individual.csv")
        if ind_efa:
            for r in ind_efa:
                r["run_id"] = eff_run
                r["plant_individual_id"] = _to_int(r.get("plant_individual_id"), default=0) or 0
                r["leaf_id"] = _to_int(r.get("leaf_id"), default=1) or 1
            counts["leaf_efa_harmonics"] = upsert_rows_generic(
                conn, "leaf_efa_harmonics", ind_efa, ("run_id", "catalogNumber", "plant_individual_id", "leaf_id")
            )

        otu_rows = _read_csv_dicts(t_dir / "morphometric_flags.csv") or _read_csv_dicts(t_dir / "leaf_efa_harmonics.csv")
        if otu_rows:
            for r in otu_rows:
                r["run_id"] = eff_run
                r["plant_individual_id"] = _to_int(r.get("plant_individual_id"), default=0) or 0
                for bcol in ("is_calibrated_mm", "is_dissected", "is_contour_outlier", "morphometric_outlier_flag", "misidentification_flag", "is_passive", "is_promoted_anchor"):
                    if bcol in r:
                        r[bcol] = _to_bool_int(r.get(bcol))
            counts["otu_morphometrics"] = upsert_rows_generic(
                conn, "otu_morphometrics", otu_rows, ("run_id", "catalogNumber", "plant_individual_id")
            )

        vis_rows = _read_csv_dicts(t_dir / "label_noise_audit.csv")
        if vis_rows:
            counts["vision_xai_audit"] = persist_vision_xai_audit(conn, eff_run, vis_rows, demo=False)

        srf_rows = _read_csv_dicts(t_dir / "multimodal_conflict_flags.csv")
        if srf_rows:
            for r in srf_rows:
                r["run_id"] = eff_run
                r["plant_individual_id"] = _to_int(r.get("plant_individual_id"), default=0) or 0
                for bcol in ("multimodal_conflict_flag", "substrate_specialist_flag", "piedmont_outcrop_specialist_flag"):
                    if bcol in r:
                        r[bcol] = _to_bool_int(r.get(bcol))
            counts["spatial_rf_predictions"] = upsert_rows_generic(
                conn, "spatial_rf_predictions", srf_rows, ("run_id", "catalogNumber", "plant_individual_id")
            )

        tri_path = t_dir / "triage_queue.csv"
        if not tri_path.is_file():
            tri_path = t_dir.parent / "outputs" / "reports" / "final_curation_triage.csv"
        tri_rows = _read_csv_dicts(tri_path)
        if tri_rows:
            for r in tri_rows:
                r["run_id"] = eff_run
                r["plant_individual_id"] = _to_int(r.get("plant_individual_id"), default=0) or 0
                if "is_completed_voucher" in r:
                    r["is_completed_voucher"] = _to_bool_int(r.get("is_completed_voucher"))
            counts["curation_triage"] = upsert_rows_generic(
                conn, "curation_triage", tri_rows, ("run_id", "catalogNumber", "plant_individual_id")
            )

        logger.info("SQLite migration summary (%s): %s", eff_run, counts)
        return counts
    finally:
        conn.close()
