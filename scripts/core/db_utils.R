# ==============================================================================
# Module: scripts/core/db_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Cross-language SQLite (WAL mode) DBI interface for R pipeline stages (03-08):
#     - Connection management with WAL, foreign keys, busy_timeout = 30000, and
#       transparent fallback to CSV when RSQLite is absent or standalone unit
#       test fixtures outside the workspace are invoked.
#     - Parameterized transactional upserts (`db_upsert_dataframe`).
#     - Single-query 120-vertex LeafMachine3 contour ingestion for Step 03.
#     - Step 04 OTU + plant_individuals + reproductive SQLite joins & updates.
#     - Step 05 environmental_cache + environmental_rtree spatial bbox queries.
#     - Step 07 & Step 08 `fetch_run_master_view()` (`v_integrative_triage_master`)
#       with seamless CSV fallback merging.
# ==============================================================================

PACKERA_SCHEMA_USER_VERSION <- 20261113L

resolve_packera_workers <- function(max_workers = 32L, n_tasks = NULL) {
  if (!identical(.Platform$OS.type, "unix") || !requireNamespace("parallel", quietly = TRUE)) {
    return(1L)
  }
  env_w <- suppressWarnings(as.integer(Sys.getenv("PACKERA_MAX_WORKERS", unset = "")))
  if (is.na(env_w) || env_w < 1L) {
    det_phys <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
    det_log <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
    base_cores <- if (!is.na(det_phys) && det_phys >= 1L) {
      as.integer(det_phys)
    } else if (!is.na(det_log) && det_log >= 1L) {
      as.integer(det_log)
    } else {
      1L
    }
    env_w <- min(as.integer(max_workers), base_cores)
  } else {
    env_w <- as.integer(env_w)
  }
  if (!is.null(n_tasks) && length(n_tasks) == 1L && is.finite(n_tasks)) {
    env_w <- min(env_w, max(1L, as.integer(n_tasks)))
  }
  max(1L, as.integer(env_w))
}

packera_mclapply <- function(X, FUN, mc.cores = 1L, mc.preschedule = TRUE) {
  n_x <- length(X)
  eff_cores <- max(1L, min(as.integer(mc.cores), max(1L, n_x)))
  if (eff_cores <= 1L || n_x <= 1L || !identical(.Platform$OS.type, "unix") || !requireNamespace("parallel", quietly = TRUE)) {
    return(lapply(X, FUN))
  }
  blas_vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
  prev_vals <- Sys.getenv(blas_vars, unset = NA_character_)
  on.exit({
    for (vn in blas_vars) {
      pv <- prev_vals[[vn]]
      if (is.na(pv)) {
        Sys.unsetenv(vn)
      } else {
        do.call(Sys.setenv, stats::setNames(list(pv), vn))
      }
    }
  }, add = TRUE)
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  res <- parallel::mclapply(X, FUN, mc.cores = eff_cores, mc.preschedule = mc.preschedule)
  bad_idx <- which(vapply(res, function(el) inherits(el, "try-error"), logical(1)))
  if (length(bad_idx) > 0L) {
    res[bad_idx] <- lapply(X[bad_idx], FUN)
  }
  res
}

.get_sqlite_user_version_r <- function(conn) {
  if (is.null(conn)) return(0L)
  uv <- tryCatch(DBI::dbGetQuery(conn, "PRAGMA user_version;"), error = function(e) NULL)
  if (is.null(uv) || nrow(uv) == 0L || ncol(uv) == 0L) return(0L)
  val <- suppressWarnings(as.integer(uv[[1]][1]))
  if (is.na(val)) 0L else val
}

.apply_packera_sqlite_pragmas_r <- function(conn, read_only = FALSE) {
  if (is.null(conn)) return(invisible(NULL))
  tryCatch({
    DBI::dbExecute(conn, "PRAGMA busy_timeout = 30000;")
    DBI::dbExecute(conn, "PRAGMA foreign_keys = ON;")
    DBI::dbExecute(conn, "PRAGMA cache_size = -262144;")
    DBI::dbExecute(conn, "PRAGMA temp_store = MEMORY;")
    tryCatch(DBI::dbExecute(conn, "PRAGMA mmap_size = 2147483648;"), error = function(e) NULL)
    if (!isTRUE(read_only)) {
      DBI::dbExecute(conn, "PRAGMA journal_mode = WAL;")
      DBI::dbExecute(conn, "PRAGMA synchronous = NORMAL;")
      tryCatch(DBI::dbExecute(conn, "PRAGMA wal_autocheckpoint = 4000;"), error = function(e) NULL)
    }
  }, error = function(e) NULL)
  invisible(NULL)
}

PACKERA_DB_LOGICAL_COLS <- c(
  "is_primary_duplicate", "is_calibrated_mm", "is_bent_stem",
  "is_dissected", "lm3_gates_pass", "mask_is_upright", "is_verified_anthesis",
  "is_rock_outcrop", "is_azonal_lithology", "high_coord_uncertainty_flag",
  "is_training_eligible_coord", "is_clim_eligible_coord", "is_edaphic_eligible_coord", "azonal_lithology_presence_buffer",
  "is_promoted_bronze", "is_promoted_anchor", "is_passive",
  "is_contour_outlier", "morphometric_outlier_flag", "misidentification_flag",
  "substrate_specialist_flag", "multimodal_conflict_flag",
  "is_completed_voucher", "is_type_specimen", "has_empirical_type_anchor",
  "is_typical_morphotype", "matches_type_lithology",
  "is_confident_label_noise", "cryptic_cluster_cohesion_pass"
)

PACKERA_DB_HARMONIC_COLS <- as.vector(vapply(
  1:24,
  function(i) paste0(c("A", "B", "C", "D"), i),
  character(4)
))
PACKERA_DB_PC_COLS <- paste0("PC", 1:15)
PACKERA_DB_SYM_PC_COLS <- paste0("Sym_PC", 1:15)
PACKERA_DB_WC_BIO_COLS <- c(
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
  "bio19_precip_coldest_qtr"
)

.normalize_semicolon_col_r <- function(v) {
  if (is.null(v)) return(NULL)
  if (is.list(v)) {
    return(vapply(v, function(el) {
      if (is.null(el) || length(el) == 0L) return(NA_character_)
      toks <- trimws(as.character(unlist(el)))
      toks <- toks[!is.na(toks) & nzchar(toks)]
      if (length(toks) == 0L) NA_character_ else paste(toks, collapse = ";")
    }, character(1)))
  }
  as.character(v)
}

.dynamic_sqlite_table_specs_r <- function(table_name) {
  if (identical(table_name, "extracted_leaves")) {
    return(c(
      ect_path = "TEXT",
      leaf_zone = "TEXT DEFAULT 'basal'",
      leaf_zone_confidence = "REAL",
      dist_from_crown_norm = "REAL",
      cy_rel_plant_bbox = "REAL",
      archetype_score = "REAL",
      bilateral_dice = "REAL",
      si_a = "REAL",
      lm3_gates_pass = "INTEGER",
      n_lamina_components = "INTEGER",
      scale_source = "TEXT",
      petiole_width_px = "REAL",
      petiole_width_mm = "REAL",
      blade_width_px = "REAL",
      blade_width_mm = "REAL",
      mask_is_upright = "INTEGER DEFAULT 1"
    ))
  }
  if (identical(table_name, "run_active_cohort")) {
    return(c(
      type_status = "TEXT",
      is_type_specimen = "INTEGER"
    ))
  }
  if (identical(table_name, "leaf_efa_harmonics")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_HARMONIC_COLS)), PACKERA_DB_HARMONIC_COLS),
      stats::setNames(rep("REAL", length(PACKERA_DB_PC_COLS)), PACKERA_DB_PC_COLS),
      stats::setNames(rep("REAL", length(PACKERA_DB_SYM_PC_COLS)), PACKERA_DB_SYM_PC_COLS),
      c(
        efd_harmonic_opt = "INTEGER",
        efd_k95_pcs = "INTEGER",
        dissection_index = "REAL",
        leaf_zone = "TEXT",
        leaf_zone_confidence = "REAL",
        dist_from_crown_norm = "REAL",
        cy_rel_plant_bbox = "REAL",
        ECT_PC1 = "REAL",
        ECT_PC2 = "REAL",
        ECT_PC3 = "REAL",
        ect_lobe_complexity = "REAL",
        ect_bilateral_asym = "REAL",
        mask_is_upright = "INTEGER",
        area_mm2 = "REAL",
        petiole_width_mm = "REAL",
        blade_width_mm = "REAL",
        type_status = "TEXT",
        is_type_specimen = "INTEGER DEFAULT 0",
        is_promoted_bronze = "INTEGER DEFAULT 0",
        is_promoted_anchor = "INTEGER DEFAULT 0",
        anchor_tier_source = "TEXT",
        authority_weight = "REAL",
        is_passive = "INTEGER DEFAULT 0",
        spatial_balance_weight = "REAL"
      )
    ))
  }
  if (identical(table_name, "otu_morphometrics")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_HARMONIC_COLS)), PACKERA_DB_HARMONIC_COLS),
      stats::setNames(rep("REAL", length(PACKERA_DB_PC_COLS)), PACKERA_DB_PC_COLS),
      stats::setNames(rep("REAL", length(PACKERA_DB_SYM_PC_COLS)), PACKERA_DB_SYM_PC_COLS),
      c(
        efd_harmonic_opt = "INTEGER",
        efd_k95_pcs = "INTEGER",
        selected_morph_features = "TEXT",
        dissection_index = "REAL",
        q25_solidity = "REAL",
        dissected_leaf_fraction = "REAL",
        basal_leaf_count = "INTEGER",
        cauline_leaf_count = "INTEGER",
        cauline_solidity_med = "REAL",
        cauline_dissection_index_med = "REAL",
        cauline_aspect_ratio_med = "REAL",
        cauline_blade_length_mm = "REAL",
        basal_vs_cauline_length_ratio = "REAL",
        basal_vs_cauline_solidity_delta = "REAL",
        heteroblasty_efa_distance = "REAL",
        Cauline_PC1 = "REAL",
        Cauline_PC2 = "REAL",
        Cauline_PC3 = "REAL",
        ECT_PC1 = "REAL",
        ECT_PC2 = "REAL",
        ECT_PC3 = "REAL",
        Cauline_ECT_PC1 = "REAL",
        Cauline_ECT_PC2 = "REAL",
        ect_lobe_complexity = "REAL",
        cauline_ect_lobe_complexity = "REAL",
        heteroblasty_ect_distance = "REAL",
        type_status = "TEXT",
        is_type_specimen = "INTEGER",
        has_empirical_type_anchor = "INTEGER",
        type_morph_mahal_d2 = "REAL",
        nearest_type_morph_taxon = "TEXT",
        type_morph_concordance = "REAL",
        is_typical_morphotype = "INTEGER",
        min_congener_mahal_d2 = "REAL",
        min_congener_taxon = "TEXT",
        open_set_novelty_p = "REAL",
        is_confident_label_noise = "INTEGER",
        is_promoted_bronze = "INTEGER",
        is_promoted_anchor = "INTEGER",
        is_passive = "INTEGER",
        spatial_balance_weight = "REAL",
        area_mm2 = "REAL",
        petiole_width_mm = "REAL",
        blade_width_mm = "REAL"
      )
    ))
  }
  if (identical(table_name, "curation_triage")) {
    return(c(
      cryptic_cluster_cohesion_pass = "INTEGER"
    ))
  }
  if (identical(table_name, "spatial_rf_predictions")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_WC_BIO_COLS)), PACKERA_DB_WC_BIO_COLS),
      c(
        rf_prediction_source = "TEXT",
        n_observed_modalities = "INTEGER",
        climatic_z_score = "REAL",
        primary_env_predictors = "TEXT",
        primary_clim_predictors = "TEXT",
        primary_edaph_predictors = "TEXT",
        type_edaphic_z_score = "REAL",
        nearest_type_niche_taxon = "TEXT",
        matches_type_lithology = "INTEGER"
      )
    ))
  }
  if (identical(table_name, "environmental_cache")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_WC_BIO_COLS)), PACKERA_DB_WC_BIO_COLS),
      c(
        solus_sand_pct = "REAL",
        solus_clay_pct = "REAL",
        solus_silt_pct = "REAL",
        solus_cec = "REAL",
        topo_hli_30m = "REAL",
        elevation_m = "REAL",
        ilr_texture_1 = "REAL",
        ilr_texture_2 = "REAL",
        fine_cfvo_pct = "REAL",
        is_clim_eligible_coord = "INTEGER",
        is_edaphic_eligible_coord = "INTEGER",
        edaphic_coord_weight = "REAL",
        soil_clay = "REAL",
        soil_silt = "REAL",
        soil_cfvo = "REAL",
        soil_clay_soilgrids = "REAL",
        soil_silt_soilgrids = "REAL",
        soil_cfvo_soilgrids = "REAL"
      )
    ))
  }
  character(0)
}

.migrate_dynamic_sqlite_table_r <- function(conn, table_name) {
  if (is.null(conn) || !db_table_exists_safe(conn, table_name)) return(FALSE)
  specs <- .dynamic_sqlite_table_specs_r(table_name)
  if (length(specs) == 0L) return(FALSE)
  existing <- tryCatch(DBI::dbListFields(conn, table_name), error = function(e) character(0))
  if (length(existing) == 0L) return(FALSE)
  missing_cols <- setdiff(names(specs), existing)
  if (length(missing_cols) == 0L) return(FALSE)

  q_tbl <- as.character(DBI::dbQuoteIdentifier(conn, table_name))
  added_any <- FALSE
  for (cn in missing_cols) {
    col_type <- specs[[cn]]
    q_col <- as.character(DBI::dbQuoteIdentifier(conn, cn))
    ok <- tryCatch({
      DBI::dbExecute(conn, sprintf("ALTER TABLE %s ADD COLUMN %s %s;", q_tbl, q_col, col_type))
      TRUE
    }, error = function(e) FALSE)
    if (isTRUE(ok)) added_any <- TRUE
  }

  added_any
}

.migrate_dynamic_sqlite_schema_r <- function(conn) {
  if (is.null(conn)) return(FALSE)
  changed <- FALSE
  for (tbl in c("extracted_leaves", "run_active_cohort", "leaf_efa_harmonics", "otu_morphometrics", "spatial_rf_predictions", "curation_triage", "environmental_cache")) {
    if (isTRUE(.migrate_dynamic_sqlite_table_r(conn, tbl))) {
      changed <- TRUE
    }
  }
  changed
}

resolve_packera_run_id <- function(run_id = NULL) {
  if (!is.null(run_id) && nzchar(trimws(as.character(run_id)[1]))) {
    return(trimws(as.character(run_id)[1]))
  }
  env_rid <- trimws(Sys.getenv("RUN_ID", unset = ""))
  if (nzchar(env_rid)) return(env_rid)
  env_rdir <- trimws(Sys.getenv("PACKERA_RUN_DIR", unset = ""))
  if (nzchar(env_rdir)) return(basename(normalizePath(env_rdir, mustWork = FALSE)))
  "canonical"
}

resolve_packera_db_path <- function(
  db_path = NULL,
  demo = FALSE
) {
  canonical_abs <- normalizePath("data/packera_pipeline.sqlite", mustWork = FALSE)
  if (isTRUE(demo)) {
    if (!is.null(db_path) && nzchar(trimws(db_path))) {
      if (identical(trimws(db_path), ":memory:")) return(":memory:")
      cand_abs <- normalizePath(trimws(db_path), mustWork = FALSE)
      if (!identical(cand_abs, canonical_abs)) return(cand_abs)
    }
    env_rdir <- trimws(Sys.getenv("PACKERA_RUN_DIR", unset = ""))
    if (nzchar(env_rdir) && dir.exists(env_rdir)) {
      return(file.path(normalizePath(env_rdir, mustWork = FALSE), "demo_ephemeral.sqlite"))
    }
    return(":memory:")
  }

  if (!is.null(db_path) && nzchar(trimws(db_path))) {
    raw <- trimws(db_path)
    if (identical(raw, ":memory:")) return(":memory:")
    return(normalizePath(raw, mustWork = FALSE))
  }
  env_db <- trimws(Sys.getenv("PACKERA_DB_PATH", unset = ""))
  if (nzchar(env_db)) {
    if (identical(env_db, ":memory:")) return(":memory:")
    return(normalizePath(env_db, mustWork = FALSE))
  }
  canonical_abs
}

is_sqlite_active_for_inputs <- function(input_paths = character(0), demo = FALSE) {
  if (isTRUE(demo)) return(FALSE)
  env_db <- trimws(Sys.getenv("PACKERA_DB_PATH", unset = ""))
  if (nzchar(env_db)) return(TRUE)

  valid_paths <- as.character(input_paths)
  valid_paths <- valid_paths[!is.na(valid_paths) & nzchar(trimws(valid_paths))]
  if (length(valid_paths) == 0L) return(TRUE)

  ws_root <- normalizePath(getwd(), mustWork = FALSE)
  ws_prefix <- paste0(ws_root, .Platform$file.sep)
  for (p in valid_paths) {
    np <- normalizePath(p, mustWork = FALSE)
    if (!identical(np, ws_root) && !startsWith(np, ws_prefix)) {
      return(FALSE)
    }
  }
  TRUE
}

.run_python_db_init_schema <- function(eff_path) {
  if (identical(eff_path, ":memory:")) return(invisible(FALSE))
  py_cands <- unique(c(
    if (file.exists(".venv/bin/python")) ".venv/bin/python" else character(0),
    Sys.which("python3")
  ))
  py_cands <- py_cands[nzchar(py_cands)]
  if (length(py_cands) == 0L) return(invisible(FALSE))
  py_code <- "import sys; from scripts.core.db import get_db_connection; c = get_db_connection(sys.argv[1], init_schema=True); c.close()"
  for (py_bin in py_cands) {
    rc <- tryCatch(
      system2(
        py_bin,
        args = c("-c", shQuote(py_code), shQuote(eff_path)),
        stdout = FALSE,
        stderr = FALSE
      ),
      error = function(e) 1L
    )
    if (identical(as.integer(rc), 0L)) return(invisible(TRUE))
  }
  invisible(FALSE)
}

connect_packera_db <- function(
  db_path = Sys.getenv("PACKERA_DB_PATH", ""),
  read_only = FALSE,
  demo = FALSE
) {
  if (!requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RSQLite", quietly = TRUE)) {
    return(NULL)
  }
  eff_path <- resolve_packera_db_path(
    db_path = if (nzchar(trimws(db_path))) db_path else NULL,
    demo = demo
  )

  if (isTRUE(read_only) && !identical(eff_path, ":memory:")) {
    if (!file.exists(eff_path)) return(NULL)
    conn <- tryCatch(
      DBI::dbConnect(RSQLite::SQLite(), dbname = eff_path, flags = RSQLite::SQLITE_RO),
      error = function(e) NULL
    )
    if (is.null(conn)) return(NULL)
    .apply_packera_sqlite_pragmas_r(conn, read_only = TRUE)
    return(conn)
  }

  if (!identical(eff_path, ":memory:")) {
    dir.create(dirname(eff_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(eff_path)) {
      # Initialize full schema via Python scripts.core.db if creating a fresh DB file from R
      .run_python_db_init_schema(eff_path)
    }
  }

  conn <- tryCatch(
    DBI::dbConnect(RSQLite::SQLite(), dbname = eff_path),
    error = function(e) NULL
  )
  if (is.null(conn)) return(NULL)
  .apply_packera_sqlite_pragmas_r(conn, read_only = FALSE)

  cur_ver <- .get_sqlite_user_version_r(conn)
  el_fields_chk <- if (db_table_exists_safe(conn, "extracted_leaves")) tryCatch(DBI::dbListFields(conn, "extracted_leaves"), error = function(e) character(0)) else c("leaf_zone", "ect_path")
  om_fields_chk <- if (db_table_exists_safe(conn, "otu_morphometrics")) tryCatch(DBI::dbListFields(conn, "otu_morphometrics"), error = function(e) character(0)) else c("cauline_leaf_count", "ECT_PC1")
  el_has_zone <- all(c("leaf_zone", "ect_path") %in% el_fields_chk)
  om_has_cauline <- all(c("cauline_leaf_count", "ECT_PC1") %in% om_fields_chk)
  if (cur_ver >= PACKERA_SCHEMA_USER_VERSION &&
      el_has_zone && om_has_cauline &&
      (identical(eff_path, ":memory:") || db_table_exists_safe(conn, "v_integrative_triage_master"))) {
    return(conn)
  }

  schema_migrated <- .migrate_dynamic_sqlite_schema_r(conn)
  if (!identical(eff_path, ":memory:")) {
    view_needs_refresh <- schema_migrated
    if (!view_needs_refresh && db_table_exists_safe(conn, "v_integrative_triage_master")) {
      v_cols <- tryCatch(DBI::dbListFields(conn, "v_integrative_triage_master"), error = function(e) character(0))
      req_v_cols <- c(
        "PC15", "Sym_PC15", "selected_morph_features", "climatic_z_score",
        "bio19_precip_coldest_qtr", "dissection_index", "q25_solidity",
        "dissected_leaf_fraction", "basal_leaf_count", "cauline_leaf_count",
        "basal_vs_cauline_length_ratio", "heteroblasty_efa_distance", "is_type_specimen",
        "ECT_PC1", "Cauline_ECT_PC1", "ect_lobe_complexity", "heteroblasty_ect_distance",
        "type_morph_mahal_d2", "nearest_type_morph_taxon", "type_morph_concordance",
        "is_typical_morphotype", "min_congener_mahal_d2", "min_congener_taxon",
        "open_set_novelty_p", "is_confident_label_noise",
        "type_edaphic_z_score", "nearest_type_niche_taxon", "matches_type_lithology",
        "cryptic_cluster_cohesion_pass"
      )
      if (length(v_cols) > 0L && !all(req_v_cols %in% v_cols)) {
        view_needs_refresh <- TRUE
      }
    }
    if (isTRUE(view_needs_refresh)) {
      disconnect_packera_db(conn)
      .run_python_db_init_schema(eff_path)
      conn <- tryCatch(
        DBI::dbConnect(RSQLite::SQLite(), dbname = eff_path),
        error = function(e) NULL
      )
      if (!is.null(conn)) {
        .apply_packera_sqlite_pragmas_r(conn, read_only = FALSE)
      }
    } else if (!is.null(conn) && db_table_exists_safe(conn, "v_integrative_triage_master")) {
      tryCatch(
        DBI::dbExecute(conn, sprintf("PRAGMA user_version = %d;", PACKERA_SCHEMA_USER_VERSION)),
        error = function(e) NULL
      )
    }
  }
  conn
}

disconnect_packera_db <- function(conn) {
  if (!is.null(conn)) {
    is_val <- tryCatch(DBI::dbIsValid(conn), error = function(e) FALSE)
    if (isTRUE(is_val)) {
      tryCatch(suppressWarnings(DBI::dbDisconnect(conn)), error = function(e) NULL)
    }
  }
  invisible(NULL)
}

db_table_exists_safe <- function(conn, table_name) {
  if (is.null(conn)) return(FALSE)
  isTRUE(tryCatch(DBI::dbExistsTable(conn, table_name), error = function(e) FALSE))
}

ensure_pipeline_run_r <- function(conn, run_id = resolve_packera_run_id(), focal_taxon = NULL, intraspecific_mode = FALSE) {
  if (is.null(conn) || !db_table_exists_safe(conn, "pipeline_runs")) return(invisible(FALSE))
  rid <- resolve_packera_run_id(run_id)
  now_iso <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  rdir <- Sys.getenv("PACKERA_RUN_DIR", unset = NA_character_)
  if (!is.na(rdir) && !nzchar(rdir)) rdir <- NA_character_
  ft <- if (!is.null(focal_taxon) && nzchar(trimws(focal_taxon))) trimws(focal_taxon) else NA_character_
  im <- if (isTRUE(intraspecific_mode)) 1L else 0L

  tryCatch(
    DBI::dbExecute(
      conn,
      paste0(
        "INSERT INTO pipeline_runs (run_id, run_dir, start_timestamp, focal_taxon, intraspecific_mode) ",
        "VALUES ($1, $2, $3, $4, $5) ",
        "ON CONFLICT(run_id) DO UPDATE SET ",
        "run_dir = COALESCE(excluded.run_dir, pipeline_runs.run_dir), ",
        "focal_taxon = COALESCE(excluded.focal_taxon, pipeline_runs.focal_taxon), ",
        "intraspecific_mode = COALESCE(excluded.intraspecific_mode, pipeline_runs.intraspecific_mode);"
      ),
      params = list(rid, rdir, now_iso, ft, im)
    ),
    error = function(e) NULL
  )
  invisible(TRUE)
}

coerce_sqlite_logical_cols <- function(df, cols = PACKERA_DB_LOGICAL_COLS) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  hit_cols <- intersect(cols, names(df))
  for (cn in hit_cols) {
    v <- df[[cn]]
    if (is.logical(v)) next
    if (is.numeric(v) || is.integer(v)) {
      df[[cn]] <- ifelse(is.na(v), NA, v != 0)
    } else {
      s <- toupper(trimws(as.character(v)))
      out <- rep(NA, length(s))
      out[s %in% c("TRUE", "T", "1", "1.0", "YES")] <- TRUE
      out[s %in% c("FALSE", "F", "0", "0.0", "NO")] <- FALSE
      df[[cn]] <- out
    }
  }
  df
}

db_upsert_dataframe <- function(
  conn,
  table_name,
  df,
  pk_cols,
  preserve_existing_on_null = TRUE
) {
  if (is.null(conn) || is.null(df) || nrow(df) == 0L) return(0L)
  if (!db_table_exists_safe(conn, table_name)) return(0L)

  .migrate_dynamic_sqlite_table_r(conn, table_name)
  tbl_fields <- DBI::dbListFields(conn, table_name)
  keep_cols <- intersect(names(df), tbl_fields)
  if (!all(pk_cols %in% keep_cols)) {
    stop(sprintf("db_upsert_dataframe: missing primary key columns (%s) for table '%s'",
                 paste(setdiff(pk_cols, keep_cols), collapse = ", "), table_name))
  }

  sub_df <- as.data.frame(df[, keep_cols, drop = FALSE], stringsAsFactors = FALSE)
  for (cn in keep_cols) {
    col_v <- sub_df[[cn]]
    if (is.list(col_v)) {
      col_v <- .normalize_semicolon_col_r(col_v)
      sub_df[[cn]] <- col_v
    }
    if (is.factor(col_v)) {
      sub_df[[cn]] <- as.character(col_v)
    } else if (is.logical(col_v)) {
      sub_df[[cn]] <- ifelse(is.na(col_v), NA_integer_, ifelse(col_v, 1L, 0L))
    } else if (is.numeric(col_v)) {
      bad_num <- is.nan(col_v) | is.infinite(col_v)
      if (any(bad_num)) col_v[bad_num] <- NA_real_
      sub_df[[cn]] <- col_v
    } else if (is.character(col_v)) {
      trimmed <- trimws(col_v)
      trimmed[!nzchar(trimmed) | tolower(trimmed) %in% c("na", "nan", "null")] <- NA_character_
      sub_df[[cn]] <- trimmed
    }
  }

  # Deduplicate on pk_cols keeping the last occurrence so transactional batch binds never collide
  pk_sig <- do.call(paste, c(lapply(pk_cols, function(k) as.character(sub_df[[k]])), sep = "\r"))
  sub_df <- sub_df[!duplicated(pk_sig, fromLast = TRUE), , drop = FALSE]

  non_pk <- setdiff(keep_cols, pk_cols)
  q_table <- DBI::dbQuoteIdentifier(conn, table_name)
  q_cols <- vapply(keep_cols, function(x) as.character(DBI::dbQuoteIdentifier(conn, x)), character(1))
  q_pks <- vapply(pk_cols, function(x) as.character(DBI::dbQuoteIdentifier(conn, x)), character(1))
  placeholders <- paste0("$", seq_along(keep_cols), collapse = ", ")

  if (length(non_pk) > 0L) {
    q_non_pk <- vapply(non_pk, function(x) as.character(DBI::dbQuoteIdentifier(conn, x)), character(1))
    if (isTRUE(preserve_existing_on_null)) {
      upd_clause <- paste(
        sprintf("%s = COALESCE(excluded.%s, %s.%s)", q_non_pk, q_non_pk, q_table, q_non_pk),
        collapse = ", "
      )
    } else {
      upd_clause <- paste(sprintf("%s = excluded.%s", q_non_pk, q_non_pk), collapse = ", ")
    }
    sql <- sprintf(
      "INSERT INTO %s (%s) VALUES (%s) ON CONFLICT (%s) DO UPDATE SET %s;",
      q_table,
      paste(q_cols, collapse = ", "),
      placeholders,
      paste(q_pks, collapse = ", "),
      upd_clause
    )
  } else {
    sql <- sprintf(
      "INSERT OR IGNORE INTO %s (%s) VALUES (%s);",
      q_table,
      paste(q_cols, collapse = ", "),
      placeholders
    )
  }

  n_sub <- nrow(sub_df)
  chunk_sz <- 5000L
  DBI::dbWithTransaction(conn, {
    if (n_sub <= chunk_sz) {
      DBI::dbExecute(conn, sql, params = unname(as.list(sub_df)))
    } else {
      starts <- seq.int(1L, n_sub, by = chunk_sz)
      for (s_idx in starts) {
        e_idx <- min(n_sub, s_idx + chunk_sz - 1L)
        DBI::dbExecute(conn, sql, params = unname(as.list(sub_df[s_idx:e_idx, , drop = FALSE])))
      }
    }
  })
  n_sub
}

# ==============================================================================
# Step 03: Single-Query 120-Vertex Contour & Leaf Metadata Ingestion + Persistence
# ==============================================================================

# Reproducible Step 03 EFA inclusion rule on read-only LeafMachine3 SQLite / manifest records:
#   - Respects explicit `efa_qc_pass` when present, and evaluates LeafMachine3's native `bilateral_dice`
#     & `si_a` with relaxed bounds on genuinely dissected Packera leaves (`dice >= 0.45`, `|si_a| < 0.48`)
#     while excluding fused overlapping clumps (`UNRESOLVED_CLUMP`) and unbridged multi-island pinnae
#     (`n_lamina_components > 1`, which remain in SQLite for multi-component ECT).
evaluate_step03_lm3_efa_qc <- function(df, allow_multi_component_ect = FALSE) {
  if (is.null(df) || nrow(df) == 0L) return(logical(0))
  n <- nrow(df)
  keep <- rep(TRUE, n)
  if ("efa_qc_pass" %in% names(df)) {
    raw_eqc <- tolower(trimws(as.character(df$efa_qc_pass)))
    keep[!is.na(raw_eqc) & raw_eqc %in% c("0", "false", "f", "no")] <- FALSE
  }
  if ("qc_flag_reason" %in% names(df)) {
    raw_qfr <- toupper(trimws(as.character(df$qc_flag_reason)))
    keep[!is.na(raw_qfr) & nzchar(raw_qfr) & !(raw_qfr %in% c("PASS", "OK", "NONE", "NA"))] <- FALSE
  }
  sol <- if ("solidity" %in% names(df)) suppressWarnings(as.numeric(df$solidity)) else rep(NA_real_, n)
  n_cc_vec <- if ("n_lamina_components" %in% names(df)) suppressWarnings(as.integer(df$n_lamina_components)) else rep(1L, n)
  is_multi_cc <- !is.na(n_cc_vec) & n_cc_vec > 1L
  is_diss <- (!is.na(sol) & is.finite(sol) & (sol < 0.72)) | is_multi_cc
  if ("lm3_gates_pass" %in% names(df)) {
    raw_gp <- tolower(trimws(as.character(df$lm3_gates_pass)))
    gp_fail <- !is.na(raw_gp) & raw_gp %in% c("0", "false", "f", "no")
    keep[!is_diss & gp_fail] <- FALSE
  }
  if ("bilateral_dice" %in% names(df)) {
    bd <- suppressWarnings(as.numeric(df$bilateral_dice))
    min_bd <- ifelse(is_diss, 0.45, 0.70)
    keep[!is.na(bd) & is.finite(bd) & bd < min_bd] <- FALSE
  }
  if ("si_a" %in% names(df)) {
    sia <- abs(suppressWarnings(as.numeric(df$si_a)))
    max_sia <- ifelse(is_diss, 0.48, 0.35)
    keep[!is.na(sia) & is.finite(sia) & sia >= max_sia] <- FALSE
  }
  if (!isTRUE(allow_multi_component_ect) && any(is_multi_cc)) {
    keep[is_multi_cc] <- FALSE
  }
  keep
}

fetch_step03_contours_from_db <- function(
  vouchers_path = NULL,
  manifest_path = NULL,
  contours_dir = NULL,
  active_catalogs = NULL,
  align_fn = NULL,
  centsize_fn = NULL,
  run_id = NULL
) {
  if (!is_sqlite_active_for_inputs(c(vouchers_path, manifest_path, contours_dir))) {
    return(NULL)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)

  conn <- connect_packera_db(read_only = TRUE)
  if (is.null(conn)) return(NULL)
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "extracted_leaves")) return(NULL)

  rid <- resolve_packera_run_id(run_id)
  has_cohort <- FALSE
  if (db_table_exists_safe(conn, "run_active_cohort")) {
    cnt_df <- tryCatch(
      DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM run_active_cohort WHERE run_id = $1;", params = list(rid)),
      error = function(e) data.frame(n = 0L)
    )
    has_cohort <- isTRUE(cnt_df$n[1] > 0L)
  }

  el_fields <- tryCatch(DBI::dbListFields(conn, "extracted_leaves"), error = function(e) character(0))
  pi_fields <- tryCatch(DBI::dbListFields(conn, "plant_individuals"), error = function(e) character(0))
  sel_ect <- if ("ect_path" %in% el_fields) "el.ect_path" else "NULL AS ect_path"
  sel_zone <- if ("leaf_zone" %in% el_fields) "el.leaf_zone" else "'basal' AS leaf_zone"
  sel_zone_conf <- if ("leaf_zone_confidence" %in% el_fields) "el.leaf_zone_confidence" else "NULL AS leaf_zone_confidence"
  sel_dist_cr <- if ("dist_from_crown_norm" %in% el_fields) "el.dist_from_crown_norm" else "NULL AS dist_from_crown_norm"
  sel_cy_rel <- if ("cy_rel_plant_bbox" %in% el_fields) "el.cy_rel_plant_bbox" else "NULL AS cy_rel_plant_bbox"
  sel_arch <- if ("archetype_score" %in% el_fields) "el.archetype_score" else "NULL AS archetype_score"
  sel_bdice <- if ("bilateral_dice" %in% el_fields) "el.bilateral_dice" else "NULL AS bilateral_dice"
  sel_sia <- if ("si_a" %in% el_fields) "el.si_a" else "NULL AS si_a"
  sel_lgates <- if ("lm3_gates_pass" %in% el_fields) "el.lm3_gates_pass" else "NULL AS lm3_gates_pass"
  sel_ncc <- if ("n_lamina_components" %in% el_fields) "el.n_lamina_components" else "1 AS n_lamina_components"
  sel_scale_src <- if ("scale_source" %in% el_fields) "el.scale_source" else "NULL AS scale_source"
  sel_upright <- if ("mask_is_upright" %in% el_fields) "el.mask_is_upright" else "NULL AS mask_is_upright"
  sel_amm2 <- if ("area_mm2" %in% el_fields) "el.area_mm2" else "NULL AS area_mm2"
  sel_bw_mm <- if ("blade_width_mm" %in% el_fields) "el.blade_width_mm" else "NULL AS blade_width_mm"
  sel_pw_mm <- if ("petiole_width_mm" %in% el_fields && "petiole_width_mm" %in% pi_fields) {
    "COALESCE(el.petiole_width_mm, pi.petiole_width_mm) AS petiole_width_mm"
  } else if ("petiole_width_mm" %in% el_fields) {
    "el.petiole_width_mm AS petiole_width_mm"
  } else if ("petiole_width_mm" %in% pi_fields) {
    "pi.petiole_width_mm AS petiole_width_mm"
  } else {
    "NULL AS petiole_width_mm"
  }
  sel_eqc <- if ("efa_qc_pass" %in% el_fields) "el.efa_qc_pass" else "NULL AS efa_qc_pass"
  sel_qfr <- if ("qc_flag_reason" %in% el_fields) "el.qc_flag_reason" else "NULL AS qc_flag_reason"

  base_sql <- paste0(
    "SELECT ",
    "  el.catalogNumber, el.plant_individual_id, el.leaf_id, el.shape_id, ",
    "  el.solidity, el.aspect_ratio, el.area_px, el.blade_area_cm2, ",
    "  ", sel_ect, ", ", sel_zone, ", ", sel_zone_conf, ", ", sel_dist_cr, ", ", sel_cy_rel, ", ",
    "  ", sel_arch, ", ", sel_bdice, ", ", sel_sia, ", ", sel_lgates, ", ", sel_ncc, ", ", sel_scale_src, ", ",
    "  ", sel_upright, ", ", sel_amm2, ", ", sel_bw_mm, ", ", sel_pw_mm, ", ", sel_eqc, ", ", sel_qfr, ", ",
    "  COALESCE(el.pixels_per_mm, v.pixels_per_mm) AS pixels_per_mm, ",
    "  COALESCE(el.is_calibrated_mm, v.is_calibrated_mm) AS is_calibrated_mm, ",
    "  COALESCE(el.petiole_length_mm, pi.petiole_length_mm) AS petiole_length_mm, ",
    "  COALESCE(el.blade_length_mm, pi.blade_length_mm) AS blade_length_mm, ",
    "  COALESCE(el.petiole_to_blade_ratio, pi.petiole_to_blade_ratio) AS petiole_to_blade_ratio, ",
    "  COALESCE(vrp.AR_inv, vrp.capitulum_aspect_ratio, vrp0.AR_inv, vrp0.capitulum_aspect_ratio) AS AR_inv, ",
    "  COALESCE(vrp.capitula_count, vrp0.capitula_count) AS capitula_count, ",
    "  COALESCE(vrp.involucre_height_px, vrp0.involucre_height_px) AS involucre_height_px, ",
    "  COALESCE(vrp.involucre_width_px, vrp0.involucre_width_px) AS involucre_width_px, ",
    "  COALESCE(vrp.involucre_height_mm, vrp0.involucre_height_mm) AS involucre_height_mm, ",
    "  COALESCE(vrp.involucre_width_mm, vrp0.involucre_width_mm) AS involucre_width_mm, ",
    "  el.contour_xy_json ",
    "FROM extracted_leaves el ",
    "LEFT JOIN vouchers v ON v.catalogNumber = el.catalogNumber ",
    "LEFT JOIN plant_individuals pi ON pi.catalogNumber = el.catalogNumber AND pi.plant_individual_id = el.plant_individual_id ",
    "LEFT JOIN voucher_reproductive_phenology vrp ON vrp.catalogNumber = el.catalogNumber AND vrp.plant_individual_id = el.plant_individual_id ",
    "LEFT JOIN voucher_reproductive_phenology vrp0 ON vrp0.catalogNumber = el.catalogNumber AND vrp0.plant_individual_id = 0 "
  )

  contour_pred <- "el.contour_xy_json IS NOT NULL AND el.contour_xy_json != ''"

  if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
    cats <- unique(trimws(as.character(active_catalogs)))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) == 0L) return(NULL)
    q_cats <- paste(vapply(cats, function(x) as.character(DBI::dbQuoteString(conn, x)), character(1)), collapse = ", ")
    sql <- paste0(
      base_sql,
      sprintf("WHERE %s AND el.catalogNumber IN (%s) ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;", contour_pred, q_cats)
    )
    rows <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  } else if (has_cohort) {
    sql <- paste0(
      base_sql,
      "INNER JOIN run_active_cohort rac ON rac.catalogNumber = el.catalogNumber AND rac.plant_individual_id = el.plant_individual_id ",
      "WHERE rac.run_id = $1 AND ", contour_pred, " ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;"
    )
    rows <- tryCatch(DBI::dbGetQuery(conn, sql, params = list(rid)), error = function(e) NULL)
  } else {
    sql <- paste0(base_sql, "WHERE ", contour_pred, " ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;")
    rows <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  }

  # Close DBI handle immediately before any parallel fork
  disconnect_packera_db(conn)
  conn <- NULL

  if (is.null(rows) || nrow(rows) == 0L) return(NULL)

  # If explicit manifest_path exists on disk with rows that are NOT yet in SQLite, defer to disk fallback;
  # otherwise build O(1) lookup tables for positional homology & cohort-adaptive blade area columns
  eff_mpath <- if (!is.null(manifest_path) && nzchar(manifest_path) && file.exists(manifest_path)) {
    manifest_path
  } else if (!is.null(vouchers_path) && nzchar(vouchers_path) && file.exists(file.path(dirname(vouchers_path), "extracted_leaf_manifest.csv"))) {
    file.path(dirname(vouchers_path), "extracted_leaf_manifest.csv")
  } else {
    ""
  }
  m_hdr <- NULL; m_idx_3k <- NULL; m_idx_2k <- NULL
  if (nzchar(eff_mpath) && file.exists(eff_mpath)) {
    m_hdr <- if (requireNamespace("data.table", quietly = TRUE)) {
      tryCatch(data.table::fread(eff_mpath, showProgress = FALSE, data.table = FALSE), error = function(e) NULL)
    } else {
      tryCatch(as.data.frame(readr::read_csv(eff_mpath, show_col_types = FALSE, progress = FALSE)), error = function(e) NULL)
    }
    if (!is.null(m_hdr) && nrow(m_hdr) > 0L && "catalogNumber" %in% names(m_hdr)) {
      if (!is.null(manifest_path) && nzchar(manifest_path) && identical(eff_mpath, manifest_path)) {
        disk_cats <- unique(trimws(as.character(m_hdr$catalogNumber)))
        if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
          disk_cats <- intersect(disk_cats, trimws(as.character(active_catalogs)))
        }
        if (length(disk_cats) > 0L && length(intersect(disk_cats, unique(rows$catalogNumber))) == 0L) {
          return(NULL)
        }
      }
      if ("leaf_id" %in% names(m_hdr)) {
        mc_v <- trimws(as.character(m_hdr$catalogNumber))
        ml_v <- suppressWarnings(as.integer(m_hdr$leaf_id))
        mp_v <- if ("plant_individual_id" %in% names(m_hdr)) ifelse(is.na(m_hdr$plant_individual_id), 0L, suppressWarnings(as.integer(m_hdr$plant_individual_id))) else rep(0L, nrow(m_hdr))
        mp_v[is.na(mp_v)] <- 0L
        k3 <- paste(mc_v, mp_v, ml_v, sep = "__")
        k2 <- paste(mc_v, ml_v, sep = "__")
        u3 <- !duplicated(k3, fromLast = TRUE); m_idx_3k <- which(u3); names(m_idx_3k) <- k3[u3]
        u2 <- !duplicated(k2, fromLast = TRUE); m_idx_2k <- which(u2); names(m_idx_2k) <- k2[u2]
      }
    }
  }

  n_rows <- nrow(rows)
  keep_qc <- evaluate_step03_lm3_efa_qc(rows, allow_multi_component_ect = TRUE) & !is.na(rows$contour_xy_json) & nzchar(as.character(rows$contour_xy_json))
  n_excl_qc <- sum(!keep_qc, na.rm = TRUE)
  if (n_excl_qc > 0L) {
    message(sprintf("Step 03 SQLite EFA QC filter: excluded %d clump/multi-island contours (retained in SQLite extracted_leaves for audit/ECT).", n_excl_qc))
  }
  if (!any(keep_qc)) return(NULL)
  if (!all(keep_qc)) {
    rows <- rows[keep_qc, , drop = FALSE]
    n_rows <- nrow(rows)
  }

  ppm_vec <- suppressWarnings(as.numeric(rows$pixels_per_mm))
  raw_cal_col <- if ("is_calibrated_mm" %in% names(rows)) rows$is_calibrated_mm else rep(1L, nrow(rows))
  cal_flag_vec <- ifelse(
    is.na(raw_cal_col),
    TRUE,
    toupper(trimws(as.character(raw_cal_col))) %in% c("1", "TRUE", "T", "YES")
  )
  raw_ss_col <- if ("scale_source" %in% names(rows)) tolower(trimws(as.character(rows$scale_source))) else rep("", nrow(rows))
  raw_ss_col[is.na(raw_ss_col)] <- ""
  is_cal_vec <- !is.na(ppm_vec) & is.finite(ppm_vec) & (ppm_vec > 0) & cal_flag_vec &
    !(raw_ss_col %in% c("mp_estimated", "predicted_from_megapixels", "uncalibrated"))
  db_ar_vec <- if ("aspect_ratio" %in% names(rows)) suppressWarnings(as.numeric(rows$aspect_ratio)) else rep(NA_real_, nrow(rows))
  db_upright_vec <- rep(NA, nrow(rows))
  if ("mask_is_upright" %in% names(rows)) {
    raw_up <- toupper(trimws(as.character(rows$mask_is_upright)))
    db_upright_vec[raw_up %in% c("1", "1.0", "TRUE", "T", "YES")] <- TRUE
    db_upright_vec[raw_up %in% c("0", "0.0", "FALSE", "F", "NO")] <- FALSE
  }
  if (!is.null(m_hdr) && "mask_is_upright" %in% names(m_hdr) && any(is.na(db_upright_vec))) {
    rc_v <- trimws(as.character(rows$catalogNumber))
    rp_v <- ifelse(is.na(rows$plant_individual_id), 0L, suppressWarnings(as.integer(rows$plant_individual_id)))
    rl_v <- suppressWarnings(as.integer(rows$leaf_id))
    mh_up <- if (!is.null(m_idx_3k)) unname(m_idx_3k[paste(rc_v, rp_v, rl_v, sep = "__")]) else rep(NA_integer_, nrow(rows))
    if (any(is.na(mh_up)) && !is.null(m_idx_2k)) {
      miss_u <- is.na(mh_up)
      mh_up[miss_u] <- unname(m_idx_2k[paste(rc_v[miss_u], rl_v[miss_u], sep = "__")])
    }
    fill_u <- which(is.na(db_upright_vec) & !is.na(mh_up))
    if (length(fill_u) > 0L) {
      m_up_raw <- toupper(trimws(as.character(m_hdr$mask_is_upright[mh_up[fill_u]])))
      db_upright_vec[fill_u[m_up_raw %in% c("1", "1.0", "TRUE", "T", "YES")]] <- TRUE
      db_upright_vec[fill_u[m_up_raw %in% c("0", "0.0", "FALSE", "F", "NO")]] <- FALSE
    }
  }
  db_area_vec <- suppressWarnings(as.numeric(rows$area_px))
  db_sol_vec <- suppressWarnings(as.numeric(rows$solidity))
  js_vec <- as.character(rows$contour_xy_json)
  has_momocs <- requireNamespace("Momocs", quietly = TRUE)

  process_contour_chunk <- function(idx_vec) {
    m_len <- length(idx_vec)
    c_list <- vector("list", m_len)
    v_ok <- logical(m_len)
    v_sol <- rep(NA_real_, m_len)
    v_di <- rep(NA_real_, m_len)
    v_cs <- rep(NA_real_, m_len)
    v_ar <- rep(NA_real_, m_len)
    v_area <- rep(NA_real_, m_len)
    v_circ <- rep(NA_real_, m_len)

    for (k in seq_len(m_len)) {
      i <- idx_vec[k]
      coo <- tryCatch(jsonlite::fromJSON(js_vec[i], simplifyMatrix = TRUE), error = function(e) NULL)
      if (is.null(coo) || !is.matrix(coo) || ncol(coo) < 2L || nrow(coo) < 10L) next
      coo <- matrix(as.numeric(coo[, 1:2]), ncol = 2L)
      if (any(!is.finite(coo))) next
      colnames(coo) <- c("x", "y")
      was_tipup_db <- FALSE
      if (exists("align_contour_horizontal_xax", mode = "function")) {
        coo <- align_contour_horizontal_xax(coo, is_upright = db_upright_vec[i])
        was_tipup_db <- isTRUE(attr(coo, "lm3_from_tipup"))
      }
      if (has_momocs) coo <- tryCatch(Momocs::coo_cw(coo), error = function(e) coo)

      x_next <- c(coo[-1L, 1L], coo[1L, 1L])
      y_next <- c(coo[-1L, 2L], coo[1L, 2L])
      area_px <- 0.5 * abs(sum(coo[, 1L] * y_next - x_next * coo[, 2L]))
      perim_px <- sum(sqrt((x_next - coo[, 1L])^2 + (y_next - coo[, 2L])^2))
      circ_val <- if (is.finite(perim_px) && perim_px > 1e-6) round(min(1.0, (4 * pi * area_px) / (perim_px^2)), 4) else NA_real_
      di_val <- if (is.finite(perim_px) && is.finite(area_px) && area_px > 1e-6) round(max(1.0, perim_px / (2.0 * sqrt(pi * area_px))), 4) else NA_real_

      if (is_cal_vec[i]) {
        coo <- coo / ppm_vec[i]
        cs_val <- if (!is.null(centsize_fn)) {
          centsize_fn(coo)
        } else {
          cx <- mean(coo[, 1L]); cy <- mean(coo[, 2L])
          sqrt(sum((coo[, 1L] - cx)^2 + (coo[, 2L] - cy)^2))
        }
      } else {
        cs_val <- NA_real_
      }

      if (!is.null(align_fn)) {
        coo <- align_fn(coo, circularity = circ_val, has_momocs = has_momocs)
      }
      attr(coo, "lm3_xax_aligned") <- TRUE
      attr(coo, "lm3_from_tipup") <- was_tipup_db

      db_ar_i <- db_ar_vec[i]
      if (is.finite(db_ar_i) && db_ar_i > 0.0) {
        ar_val <- round(as.numeric(db_ar_i), 4)
      } else if (was_tipup_db && diff(range(coo[, 2L], na.rm = TRUE)) > 1e-6) {
        ar_val <- round(as.numeric(diff(range(coo[, 1L], na.rm = TRUE)) / diff(range(coo[, 2L], na.rm = TRUE))), 4)
      } else {
        lw <- if (has_momocs) tryCatch(Momocs::coo_lw(coo), error = function(e) NULL) else NULL
        if (is.null(lw) || length(lw) < 2L || is.na(lw[2]) || lw[2] <= 1e-6) {
          coo_centered <- scale(coo, center = TRUE, scale = FALSE)
          svd_rot <- coo_centered %*% svd(coo_centered)$v
          lw <- c(diff(range(svd_rot[, 1], na.rm = TRUE)), diff(range(svd_rot[, 2], na.rm = TRUE)))
        }
        ar_val <- round(as.numeric(lw[1] / lw[2]), 4)
      }

      db_a <- db_area_vec[i]
      eff_area <- if (!is.na(db_a) && is.finite(db_a) && db_a > 0) db_a else round(area_px, 1)
      sol_val <- db_sol_vec[i]
      if (is.na(sol_val) || !is.finite(sol_val)) {
        ch_idx <- tryCatch(grDevices::chull(coo[, 1L], coo[, 2L]), error = function(e) integer(0))
        if (length(ch_idx) >= 3L) {
          hc <- coo[ch_idx, , drop = FALSE]
          h_area <- 0.5 * abs(sum(hc[, 1L] * c(hc[-1L, 2L], hc[1L, 2L]) - c(hc[-1L, 1L], hc[1L, 1L]) * hc[, 2L]))
          if (h_area > 1e-6) sol_val <- round(min(1.0, area_px / h_area), 4)
        }
      }

      v_ok[k] <- TRUE
      c_list[[k]] <- coo
      v_sol[k] <- sol_val
      v_di[k] <- di_val
      v_cs[k] <- as.numeric(cs_val)
      v_ar[k] <- as.numeric(ar_val)
      v_area[k] <- as.numeric(eff_area)
      v_circ[k] <- as.numeric(circ_val)
    }

    list(
      ok = v_ok, coo = c_list, sol = v_sol, di = v_di,
      cs = v_cs, ar = v_ar, area = v_area, circ = v_circ
    )
  }

  n_workers <- if (n_rows >= 160L) resolve_packera_workers(32L, ceiling(n_rows / 80L)) else 1L
  if (n_workers <= 1L) {
    merged <- process_contour_chunk(seq_len(n_rows))
  } else {
    chunk_ids <- split(seq_len(n_rows), cut(seq_len(n_rows), breaks = n_workers, labels = FALSE))
    chunk_out <- packera_mclapply(chunk_ids, process_contour_chunk, mc.cores = n_workers)
    merged <- list(
      ok = unlist(lapply(chunk_out, `[[`, "ok"), use.names = FALSE),
      coo = unlist(lapply(chunk_out, `[[`, "coo"), recursive = FALSE, use.names = FALSE),
      sol = unlist(lapply(chunk_out, `[[`, "sol"), use.names = FALSE),
      di = unlist(lapply(chunk_out, `[[`, "di"), use.names = FALSE),
      cs = unlist(lapply(chunk_out, `[[`, "cs"), use.names = FALSE),
      ar = unlist(lapply(chunk_out, `[[`, "ar"), use.names = FALSE),
      area = unlist(lapply(chunk_out, `[[`, "area"), use.names = FALSE),
      circ = unlist(lapply(chunk_out, `[[`, "circ"), use.names = FALSE)
    )
  }

  ok_idx <- which(merged$ok)
  if (length(ok_idx) == 0L) return(NULL)

  sub_rows <- rows[ok_idx, , drop = FALSE]
  n_kept <- length(ok_idx)
  cat_vec <- as.character(sub_rows$catalogNumber)
  pid_vec <- ifelse(is.na(sub_rows$plant_individual_id), 0L, as.integer(sub_rows$plant_individual_id))
  lid_vec <- ifelse(is.na(sub_rows$leaf_id), ok_idx, as.integer(sub_rows$leaf_id))
  raw_sid_vec <- as.character(sub_rows$shape_id)
  sid_vec <- ifelse(!is.na(raw_sid_vec) & nzchar(raw_sid_vec), raw_sid_vec, sprintf("%s_p%d_leaf%d", cat_vec, pid_vec, lid_vec))
  if (anyDuplicated(sid_vec) > 0L) {
    seen_sids <- new.env(hash = TRUE, parent = emptyenv())
    for (k in seq_len(n_kept)) {
      s_cur <- sid_vec[k]
      if (exists(s_cur, envir = seen_sids, inherits = FALSE)) {
        s_cur <- sprintf("%s_%d", s_cur, ok_idx[k])
        sid_vec[k] <- s_cur
      }
      assign(s_cur, TRUE, envir = seen_sids)
    }
  }

  m_hit_vec <- rep(NA_integer_, n_kept)
  if (!is.null(m_hdr)) {
    t_cat <- trimws(cat_vec)
    if (!is.null(m_idx_3k)) {
      m_hit_vec <- unname(m_idx_3k[paste(t_cat, pid_vec, lid_vec, sep = "__")])
    }
    miss_m <- is.na(m_hit_vec)
    if (any(miss_m) && !is.null(m_idx_2k)) {
      m_hit_vec[miss_m] <- unname(m_idx_2k[paste(t_cat[miss_m], lid_vec[miss_m], sep = "__")])
    }
  }
  pull_manifest_num <- function(cn, fb_vec = rep(NA_real_, n_kept)) {
    out_v <- fb_vec
    if (!is.null(m_hdr) && cn %in% names(m_hdr) && any(!is.na(m_hit_vec))) {
      hit_mask <- !is.na(m_hit_vec)
      cand_v <- suppressWarnings(as.numeric(m_hdr[[cn]][m_hit_vec[hit_mask]]))
      fin_m <- is.finite(cand_v)
      if (any(fin_m)) {
        hit_pos <- which(hit_mask)
        out_v[hit_pos[fin_m]] <- cand_v[fin_m]
      }
    }
    out_v
  }

  db_ba_vec <- if ("blade_area_cm2" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_area_cm2)) else rep(NA_real_, n_kept)
  db_ba_vec[!is.finite(db_ba_vec)] <- NA_real_
  db_cy_rel_vec <- if ("cy_rel_plant_bbox" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$cy_rel_plant_bbox)) else rep(NA_real_, n_kept)
  db_cy_rel_vec[!is.finite(db_cy_rel_vec)] <- NA_real_
  db_zone_conf_vec <- if ("leaf_zone_confidence" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$leaf_zone_confidence)) else rep(NA_real_, n_kept)
  db_zone_conf_vec[!is.finite(db_zone_conf_vec)] <- NA_real_
  db_dist_cr_vec <- if ("dist_from_crown_norm" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$dist_from_crown_norm)) else rep(NA_real_, n_kept)
  db_dist_cr_vec[!is.finite(db_dist_cr_vec)] <- NA_real_
  db_zone_vec <- if ("leaf_zone" %in% names(sub_rows)) trimws(as.character(sub_rows$leaf_zone)) else rep("basal", n_kept)
  db_zone_vec[is.na(db_zone_vec) | !nzchar(db_zone_vec)] <- "basal"
  if (!is.null(m_hdr) && "leaf_zone" %in% names(m_hdr) && any(!is.na(m_hit_vec))) {
    hit_m_z <- !is.na(m_hit_vec)
    cand_z <- trimws(as.character(m_hdr[["leaf_zone"]][m_hit_vec[hit_m_z]]))
    valid_z <- !is.na(cand_z) & nzchar(cand_z)
    if (any(valid_z)) {
      hit_pos_z <- which(hit_m_z)
      db_zone_vec[hit_pos_z[valid_z]] <- cand_z[valid_z]
    }
  }
  db_ect_vec <- if ("ect_path" %in% names(sub_rows)) trimws(as.character(sub_rows$ect_path)) else rep(NA_character_, n_kept)
  db_ect_vec[is.na(db_ect_vec) | !nzchar(db_ect_vec)] <- NA_character_
  if (!is.null(m_hdr) && "ect_path" %in% names(m_hdr) && any(!is.na(m_hit_vec))) {
    hit_m_e <- !is.na(m_hit_vec) & is.na(db_ect_vec)
    if (any(hit_m_e)) {
      cand_e <- trimws(as.character(m_hdr[["ect_path"]][m_hit_vec[hit_m_e]]))
      valid_e <- !is.na(cand_e) & nzchar(cand_e)
      if (any(valid_e)) {
        hit_pos_e <- which(hit_m_e)
        db_ect_vec[hit_pos_e[valid_e]] <- cand_e[valid_e]
      }
    }
  }
  sol_kept <- merged$sol[ok_idx]
  eff_diss_vec <- is.finite(sol_kept) & (sol_kept < 0.72)
  cs_kept <- merged$cs[ok_idx]
  log_cs_kept <- ifelse(!is.na(cs_kept) & cs_kept > 0, round(log10(cs_kept), 6), NA_real_)
  arch_kept <- if ("archetype_score" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$archetype_score)) else rep(NA_real_, n_kept)

  coords_list <- merged$coo[ok_idx]
  names(coords_list) <- sid_vec
  meta_df <- data.frame(
    shape_id = sid_vec,
    catalogNumber = cat_vec,
    plant_individual_id = pid_vec,
    leaf_id = lid_vec,
    n_lamina_components = if ("n_lamina_components" %in% names(sub_rows)) pmax(1L, ifelse(is.na(suppressWarnings(as.integer(sub_rows$n_lamina_components))), 1L, suppressWarnings(as.integer(sub_rows$n_lamina_components)))) else rep(1L, n_kept),
    solidity = sol_kept,
    dissection_index = merged$di[ok_idx],
    archetype_score = arch_kept,
    pixels_per_mm = ppm_vec[ok_idx],
    is_calibrated_mm = is_cal_vec[ok_idx],
    scale_source = ifelse(is_cal_vec[ok_idx], "ruler_calibrated", ifelse(!is.na(ppm_vec[ok_idx]) & ppm_vec[ok_idx] > 0, "mp_estimated", "uncalibrated")),
    mask_is_upright = ifelse(is.na(db_upright_vec[ok_idx]), isTRUE(attr(coords_list[[1L]], "lm3_from_tipup")), as.logical(db_upright_vec[ok_idx])),
    area_mm2 = {
      db_am <- if ("area_mm2" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$area_mm2)) else rep(NA_real_, n_kept)
      db_am[!is.finite(db_am)] <- NA_real_
      am_m <- pull_manifest_num("area_mm2", db_am)
      ifelse(is.finite(am_m), am_m, ifelse(is_cal_vec[ok_idx] & is.finite(ppm_vec[ok_idx]) & ppm_vec[ok_idx] > 0, round(merged$area[ok_idx] / (ppm_vec[ok_idx]^2), 2), NA_real_))
    },
    petiole_length_mm = suppressWarnings(as.numeric(sub_rows$petiole_length_mm)),
    petiole_width_mm = pull_manifest_num("petiole_width_mm", if ("petiole_width_mm" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$petiole_width_mm)) else rep(NA_real_, n_kept)),
    blade_length_mm = suppressWarnings(as.numeric(sub_rows$blade_length_mm)),
    blade_width_mm = pull_manifest_num("blade_width_mm", if ("blade_width_mm" %in% names(sub_rows)) suppressWarnings(as.numeric(sub_rows$blade_width_mm)) else rep(NA_real_, n_kept)),
    petiole_to_blade_ratio = suppressWarnings(as.numeric(sub_rows$petiole_to_blade_ratio)),
    is_dissected = eff_diss_vec,
    ect_path = db_ect_vec,
    leaf_zone = db_zone_vec,
    leaf_zone_confidence = pull_manifest_num("leaf_zone_confidence", db_zone_conf_vec),
    dist_from_crown_norm = pull_manifest_num("dist_from_crown_norm", db_dist_cr_vec),
    cy_rel_plant_bbox = pull_manifest_num("cy_rel_plant_bbox", db_cy_rel_vec),
    cy_rel_sheet = pull_manifest_num("cy_rel_sheet"),
    blade_area_cm2 = pull_manifest_num("blade_area_cm2", db_ba_vec),
    AR_inv = suppressWarnings(as.numeric(sub_rows$AR_inv)),
    capitula_count = suppressWarnings(as.numeric(sub_rows$capitula_count)),
    involucre_height_px = suppressWarnings(as.numeric(sub_rows$involucre_height_px)),
    involucre_width_px = suppressWarnings(as.numeric(sub_rows$involucre_width_px)),
    involucre_height_mm = suppressWarnings(as.numeric(sub_rows$involucre_height_mm)),
    involucre_width_mm = suppressWarnings(as.numeric(sub_rows$involucre_width_mm)),
    centroid_size = cs_kept,
    log_centsize = log_cs_kept,
    aspect_ratio = merged$ar[ok_idx],
    area_px = merged$area[ok_idx],
    circularity = merged$circ[ok_idx],
    mask_source = paste0("sqlite://extracted_leaves/", sid_vec),
    stringsAsFactors = FALSE
  )
  rownames(meta_df) <- sid_vec
  list(
    coords_list = coords_list,
    meta_df = meta_df
  )
}

.normalize_morph_dynamic_cols_r <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  for (ic in intersect(c("plant_individual_id", "leaf_id", "efd_harmonic_opt", "efd_k95_pcs", "basal_leaf_count", "cauline_leaf_count"), names(df))) {
    df[[ic]] <- suppressWarnings(as.integer(df[[ic]]))
  }
  if ("selected_morph_features" %in% names(df)) {
    df[["selected_morph_features"]] <- .normalize_semicolon_col_r(df[["selected_morph_features"]])
  }
  dyn_num_cols <- intersect(
    c(PACKERA_DB_HARMONIC_COLS, PACKERA_DB_PC_COLS, PACKERA_DB_SYM_PC_COLS,
      "dissection_index", "q25_solidity", "dissected_leaf_fraction",
      "leaf_zone_confidence", "dist_from_crown_norm", "cy_rel_plant_bbox",
      "cauline_solidity_med", "cauline_dissection_index_med", "cauline_aspect_ratio_med",
      "cauline_blade_length_mm", "basal_vs_cauline_length_ratio",
      "basal_vs_cauline_solidity_delta", "heteroblasty_efa_distance",
      "Cauline_PC1", "Cauline_PC2", "Cauline_PC3",
      "ECT_PC1", "ECT_PC2", "ECT_PC3",
      "Cauline_ECT_PC1", "Cauline_ECT_PC2",
      "ect_lobe_complexity", "cauline_ect_lobe_complexity",
      "ect_bilateral_asym", "heteroblasty_ect_distance",
      "area_mm2", "petiole_width_mm", "blade_width_mm", "spatial_balance_weight"),
    names(df)
  )
  for (nc in dyn_num_cols) {
    df[[nc]] <- suppressWarnings(as.numeric(df[[nc]]))
  }
  df
}

.prune_stale_run_rows_r <- function(conn, table_name, rid, df) {
  if (is.null(conn) || is.null(df) || nrow(df) == 0L || !("catalogNumber" %in% names(df))) return(invisible(0L))
  cats <- trimws(as.character(df$catalogNumber))
  pids <- if ("plant_individual_id" %in% names(df)) ifelse(is.na(df$plant_individual_id), 0L, suppressWarnings(as.integer(df$plant_individual_id))) else rep(0L, nrow(df))
  pids[is.na(pids)] <- 0L
  keep_keys <- unique(paste(cats, pids, sep = "\r"))
  ex <- tryCatch(DBI::dbGetQuery(conn, sprintf("SELECT rowid AS _rid, catalogNumber, plant_individual_id FROM %s WHERE run_id = $1;", table_name), params = list(rid)), error = function(e) NULL)
  if (is.null(ex) || nrow(ex) == 0L) return(invisible(0L))
  ex_pids <- ifelse(is.na(ex$plant_individual_id), 0L, suppressWarnings(as.integer(ex$plant_individual_id)))
  ex_pids[is.na(ex_pids)] <- 0L
  ex_keys <- paste(trimws(as.character(ex$catalogNumber)), ex_pids, sep = "\r")
  drop_rids <- as.integer(ex[["_rid"]][!(ex_keys %in% keep_keys)])
  if (length(drop_rids) > 0L) {
    for (s in seq.int(1L, length(drop_rids), by = 500L)) {
      chunk <- drop_rids[s:min(length(drop_rids), s + 499L)]
      tryCatch(DBI::dbExecute(conn, sprintf("DELETE FROM %s WHERE rowid IN (%s);", table_name, paste(chunk, collapse = ","))), error = function(e) NULL)
    }
  }
  invisible(length(drop_rids))
}

persist_step03_to_sqlite <- function(
  efa_df,
  specimen_df,
  vouchers_path = NULL,
  run_id = NULL
) {
  if (!is_sqlite_active_for_inputs(vouchers_path)) return(invisible(FALSE))
  conn <- connect_packera_db(read_only = FALSE)
  if (is.null(conn)) return(invisible(FALSE))
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  if (!is.null(efa_df) && nrow(efa_df) > 0L) {
    ind_out <- .normalize_morph_dynamic_cols_r(as.data.frame(efa_df, stringsAsFactors = FALSE))
    ind_out$run_id <- rid
    if (!("plant_individual_id" %in% names(ind_out))) ind_out$plant_individual_id <- 0L
    ind_out$plant_individual_id[is.na(ind_out$plant_individual_id)] <- 0L
    tryCatch(DBI::dbExecute(conn, "DELETE FROM leaf_efa_harmonics WHERE run_id = $1;", params = list(rid)), error = function(e) NULL)
    db_upsert_dataframe(
      conn,
      table_name = "leaf_efa_harmonics",
      df = ind_out,
      pk_cols = c("run_id", "catalogNumber", "plant_individual_id", "leaf_id"),
      preserve_existing_on_null = FALSE
    )
  }

  if (!is.null(specimen_df) && nrow(specimen_df) > 0L) {
    otu_out <- .normalize_morph_dynamic_cols_r(as.data.frame(specimen_df, stringsAsFactors = FALSE))
    otu_out$run_id <- rid
    if (!("plant_individual_id" %in% names(otu_out))) otu_out$plant_individual_id <- 0L
    otu_out$plant_individual_id[is.na(otu_out$plant_individual_id)] <- 0L
    if (("is_promoted_bronze" %in% names(otu_out)) && !("is_promoted_anchor" %in% names(otu_out))) {
      otu_out$is_promoted_anchor <- as.logical(otu_out$is_promoted_bronze)
    } else if (("is_promoted_anchor" %in% names(otu_out)) && !("is_promoted_bronze" %in% names(otu_out))) {
      otu_out$is_promoted_bronze <- as.logical(otu_out$is_promoted_anchor)
    }
    if (!("is_passive" %in% names(otu_out)) && ("authority_weight" %in% names(otu_out))) {
      aw_v <- suppressWarnings(as.numeric(otu_out$authority_weight))
      otu_out$is_passive <- !is.finite(aw_v) | (aw_v <= 0.0)
    }
    tryCatch(DBI::dbExecute(conn, "DELETE FROM otu_morphometrics WHERE run_id = $1;", params = list(rid)), error = function(e) NULL)
    db_upsert_dataframe(
      conn,
      table_name = "otu_morphometrics",
      df = otu_out,
      pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
      preserve_existing_on_null = FALSE
    )
  }
  invisible(TRUE)
}

# ==============================================================================
# Step 04: OTU Morphometrics + Plant Individuals + Reproductive Joins & Updates
# ==============================================================================

fetch_step04_morphometrics_from_db <- function(
  vouchers_path = NULL,
  harmonics_path = NULL,
  active_catalogs = NULL,
  run_id = NULL
) {
  if (!is_sqlite_active_for_inputs(c(vouchers_path, harmonics_path))) return(NULL)
  conn <- connect_packera_db(read_only = TRUE)
  if (is.null(conn)) return(NULL)
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "otu_morphometrics")) return(NULL)

  rid <- resolve_packera_run_id(run_id)
  cnt <- tryCatch(
    DBI::dbGetQuery(conn, "SELECT COUNT(*) AS n FROM otu_morphometrics WHERE run_id = $1;", params = list(rid))$n[1],
    error = function(e) 0L
  )
  if (is.na(cnt) || cnt == 0L) return(NULL)

  sql <- paste0(
    "SELECT ",
    "  om.*, ",
    "  COALESCE(om.petiole_length_mm, pi.petiole_length_mm) AS petiole_length_mm_joined, ",
    "  COALESCE(om.blade_length_mm, pi.blade_length_mm) AS blade_length_mm_joined, ",
    "  COALESCE(om.petiole_to_blade_ratio, pi.petiole_to_blade_ratio) AS petiole_to_blade_ratio_joined, ",
    "  COALESCE(om.AR_inv, vrp.AR_inv, vrp.capitulum_aspect_ratio, vrp0.AR_inv, vrp0.capitulum_aspect_ratio) AS AR_inv_joined, ",
    "  COALESCE(vrp.capitula_count, vrp0.capitula_count) AS capitula_count_joined, ",
    "  COALESCE(vrp.involucre_height_px, vrp0.involucre_height_px) AS involucre_height_px_joined, ",
    "  COALESCE(vrp.involucre_width_px, vrp0.involucre_width_px) AS involucre_width_px_joined, ",
    "  COALESCE(vrp.involucre_height_mm, vrp0.involucre_height_mm) AS involucre_height_mm_joined, ",
    "  COALESCE(vrp.involucre_width_mm, vrp0.involucre_width_mm) AS involucre_width_mm_joined ",
    "FROM otu_morphometrics om ",
    "LEFT JOIN plant_individuals pi ON pi.catalogNumber = om.catalogNumber AND pi.plant_individual_id = om.plant_individual_id ",
    "LEFT JOIN voucher_reproductive_phenology vrp ON vrp.catalogNumber = om.catalogNumber AND vrp.plant_individual_id = om.plant_individual_id ",
    "LEFT JOIN voucher_reproductive_phenology vrp0 ON vrp0.catalogNumber = om.catalogNumber AND vrp0.plant_individual_id = 0 ",
    "WHERE om.run_id = $1 ",
    "ORDER BY om.catalogNumber, om.plant_individual_id;"
  )
  df <- tryCatch(DBI::dbGetQuery(conn, sql, params = list(rid)), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return(NULL)

  if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
    cats <- unique(trimws(as.character(active_catalogs)))
    df <- df[trimws(as.character(df$catalogNumber)) %in% cats, , drop = FALSE]
    if (nrow(df) == 0L) return(NULL)
  }

  for (pair in list(
    c("petiole_length_mm", "petiole_length_mm_joined"),
    c("blade_length_mm", "blade_length_mm_joined"),
    c("petiole_to_blade_ratio", "petiole_to_blade_ratio_joined"),
    c("AR_inv", "AR_inv_joined"),
    c("capitula_count", "capitula_count_joined"),
    c("involucre_height_px", "involucre_height_px_joined"),
    c("involucre_width_px", "involucre_width_px_joined"),
    c("involucre_height_mm", "involucre_height_mm_joined"),
    c("involucre_width_mm", "involucre_width_mm_joined")
  )) {
    df[[pair[1]]] <- df[[pair[2]]]
    df[[pair[2]]] <- NULL
  }
  df$run_id <- NULL
  coerce_sqlite_logical_cols(df)
}

persist_step04_to_sqlite <- function(
  out_df,
  vouchers_path = NULL,
  run_id = NULL,
  focal_taxon = NULL,
  intraspecific_mode = FALSE
) {
  if (!is_sqlite_active_for_inputs(vouchers_path) || is.null(out_df) || nrow(out_df) == 0L) {
    return(invisible(FALSE))
  }
  conn <- connect_packera_db(read_only = FALSE)
  if (is.null(conn)) return(invisible(FALSE))
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid, focal_taxon = focal_taxon, intraspecific_mode = intraspecific_mode)

  df_db <- .normalize_morph_dynamic_cols_r(as.data.frame(out_df, stringsAsFactors = FALSE))
  df_db$run_id <- rid
  if (!("plant_individual_id" %in% names(df_db))) df_db$plant_individual_id <- 0L
  df_db$plant_individual_id[is.na(df_db$plant_individual_id)] <- 0L

  db_upsert_dataframe(
    conn,
    table_name = "otu_morphometrics",
    df = df_db,
    pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
    preserve_existing_on_null = TRUE
  )
  .prune_stale_run_rows_r(conn, "otu_morphometrics", rid, df_db)

  if (db_table_exists_safe(conn, "run_active_cohort") && db_table_exists_safe(conn, "vouchers")) {
    # Update authority_weight and is_promoted_bronze on existing cohort rows without violating vouchers FK
    v_cats <- tryCatch(DBI::dbGetQuery(conn, "SELECT catalogNumber FROM vouchers;")$catalogNumber, error = function(e) character(0))
    coh_sub <- df_db[df_db$catalogNumber %in% v_cats, , drop = FALSE]
    if (nrow(coh_sub) > 0L) {
      if ("is_promoted_anchor" %in% names(coh_sub) && !("is_promoted_bronze" %in% names(coh_sub))) {
        coh_sub$is_promoted_bronze <- coh_sub$is_promoted_anchor
      }
      coh_cols <- intersect(
        c("run_id", "catalogNumber", "plant_individual_id", "species_standardized",
          "determiner_tier", "type_status", "is_type_specimen",
          "is_promoted_bronze", "authority_weight"),
        names(coh_sub)
      )
      db_upsert_dataframe(
        conn,
        table_name = "run_active_cohort",
        df = coh_sub[, coh_cols, drop = FALSE],
        pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
        preserve_existing_on_null = TRUE
      )
    }
  }
  invisible(TRUE)
}

# ==============================================================================
# Step 05: Environmental Cache + R*Tree Spatial Index & Spatial RF Predictions
# ==============================================================================

query_environmental_cache_sqlite <- function(
  catalog_numbers = NULL,
  lat_lon_keys = NULL,
  bbox = NULL,
  cache_path = NULL,
  db_path = NULL,
  demo = FALSE
) {
  has_explicit_db <- !is.null(db_path) && nzchar(trimws(as.character(db_path)))
  if (isTRUE(demo) || (!has_explicit_db && !is_sqlite_active_for_inputs(cache_path, demo = demo))) return(NULL)
  conn <- connect_packera_db(db_path = if (has_explicit_db) as.character(db_path) else "", read_only = TRUE, demo = FALSE)
  if (is.null(conn)) return(NULL)
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "environmental_cache")) return(NULL)

  if (!is.null(bbox) && length(bbox) == 4L && all(is.finite(as.numeric(bbox))) &&
      db_table_exists_safe(conn, "environmental_rtree")) {
    bb <- as.numeric(bbox)
    res <- tryCatch(
      DBI::dbGetQuery(
        conn,
        paste0(
          "SELECT ec.* FROM environmental_cache ec ",
          "INNER JOIN environmental_rtree rt ON ec.rowid = rt.id ",
          "WHERE rt.min_lon >= $1 AND rt.max_lon <= $2 AND rt.min_lat >= $3 AND rt.max_lat <= $4;"
        ),
        params = list(bb[1], bb[2], bb[3], bb[4])
      ),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0L) return(NULL)
    return(coerce_sqlite_logical_cols(res))
  }

  clauses <- character(0)
  if (!is.null(catalog_numbers) && length(catalog_numbers) > 0L) {
    cats <- unique(trimws(as.character(catalog_numbers)))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) > 0L) {
      qc <- paste(vapply(cats, function(x) as.character(DBI::dbQuoteString(conn, x)), character(1)), collapse = ", ")
      clauses <- c(clauses, sprintf("catalogNumber IN (%s)", qc))
    }
  }
  if (!is.null(lat_lon_keys) && length(lat_lon_keys) > 0L) {
    keys <- unique(trimws(as.character(lat_lon_keys)))
    keys <- keys[!is.na(keys) & nzchar(keys)]
    if (length(keys) > 0L) {
      qk <- paste(vapply(keys, function(x) as.character(DBI::dbQuoteString(conn, x)), character(1)), collapse = ", ")
      clauses <- c(clauses, sprintf("lat_lon_key IN (%s)", qk))
    }
  }

  sql <- if (length(clauses) > 0L) {
    sprintf("SELECT * FROM environmental_cache WHERE %s;", paste(clauses, collapse = " OR "))
  } else {
    "SELECT * FROM environmental_cache;"
  }
  res <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  coerce_sqlite_logical_cols(res)
}

upsert_environmental_cache_sqlite <- function(
  df,
  cache_path = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(df) || nrow(df) == 0L) return(invisible(0L))
  if (!is_sqlite_active_for_inputs(cache_path, demo = demo)) return(invisible(0L))

  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  if (is.null(conn)) return(invisible(0L))
  on.exit(disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "environmental_cache")) return(invisible(0L))

  env_df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (bc in intersect(PACKERA_DB_WC_BIO_COLS, names(env_df))) {
    env_df[[bc]] <- suppressWarnings(as.numeric(env_df[[bc]]))
  }
  if (!("latitude" %in% names(env_df)) && "decimalLatitude" %in% names(env_df)) {
    env_df$latitude <- suppressWarnings(as.numeric(env_df$decimalLatitude))
  }
  if (!("longitude" %in% names(env_df)) && "decimalLongitude" %in% names(env_df)) {
    env_df$longitude <- suppressWarnings(as.numeric(env_df$decimalLongitude))
  }
  lat_num <- suppressWarnings(as.numeric(env_df$latitude))
  lon_num <- suppressWarnings(as.numeric(env_df$longitude))
  valid_xy <- !is.na(lat_num) & !is.na(lon_num) & is.finite(lat_num) & is.finite(lon_num)

  if (!("lat_lon_key" %in% names(env_df))) {
    env_df$lat_lon_key <- NA_character_
  }
  missing_key <- is.na(env_df$lat_lon_key) | !nzchar(trimws(as.character(env_df$lat_lon_key)))
  fill_idx <- which(missing_key & valid_xy)
  if (length(fill_idx) > 0L) {
    env_df$lat_lon_key[fill_idx] <- sprintf("%.5f_%.5f", lat_num[fill_idx], lon_num[fill_idx])
  }
  env_df <- env_df[!is.na(env_df$lat_lon_key) & nzchar(trimws(as.character(env_df$lat_lon_key))), , drop = FALSE]
  if (nrow(env_df) == 0L) return(invisible(0L))

  n_up <- db_upsert_dataframe(
    conn,
    table_name = "environmental_cache",
    df = env_df,
    pk_cols = c("lat_lon_key"),
    preserve_existing_on_null = TRUE
  )
  invisible(n_up)
}

persist_step05_to_sqlite <- function(
  vouchers_df,
  vouchers_path = NULL,
  run_id = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(vouchers_df) || nrow(vouchers_df) == 0L) return(invisible(FALSE))
  if (!is_sqlite_active_for_inputs(vouchers_path, demo = demo)) return(invisible(FALSE))

  srf_df <- as.data.frame(vouchers_df, stringsAsFactors = FALSE)
  upsert_environmental_cache_sqlite(srf_df, cache_path = vouchers_path, demo = demo)

  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  if (is.null(conn)) return(invisible(FALSE))
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  srf_df$run_id <- rid
  if (!("plant_individual_id" %in% names(srf_df))) srf_df$plant_individual_id <- 0L
  srf_df$plant_individual_id <- suppressWarnings(as.integer(srf_df$plant_individual_id))
  srf_df$plant_individual_id[is.na(srf_df$plant_individual_id)] <- 0L

  for (pc in intersect(c("primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors"), names(srf_df))) {
    srf_df[[pc]] <- .normalize_semicolon_col_r(srf_df[[pc]])
  }
  for (nc in intersect(c("edaphic_z_score", "climatic_z_score", "type_edaphic_z_score", PACKERA_DB_WC_BIO_COLS), names(srf_df))) {
    srf_df[[nc]] <- suppressWarnings(as.numeric(srf_df[[nc]]))
  }

  if (!("lat_lon_key" %in% names(srf_df))) {
    lat_src <- if ("latitude" %in% names(srf_df)) srf_df$latitude else if ("decimalLatitude" %in% names(srf_df)) srf_df$decimalLatitude else NULL
    lon_src <- if ("longitude" %in% names(srf_df)) srf_df$longitude else if ("decimalLongitude" %in% names(srf_df)) srf_df$decimalLongitude else NULL
    if (!is.null(lat_src) && !is.null(lon_src)) {
      lat_n <- suppressWarnings(as.numeric(lat_src))
      lon_n <- suppressWarnings(as.numeric(lon_src))
      ok_ll <- !is.na(lat_n) & !is.na(lon_n) & is.finite(lat_n) & is.finite(lon_n)
      srf_df$lat_lon_key <- NA_character_
      srf_df$lat_lon_key[ok_ll] <- sprintf("%.5f_%.5f", lat_n[ok_ll], lon_n[ok_ll])
    }
  }

  db_upsert_dataframe(
    conn,
    table_name = "spatial_rf_predictions",
    df = srf_df,
    pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
    preserve_existing_on_null = FALSE
  )
  .prune_stale_run_rows_r(conn, "spatial_rf_predictions", rid, srf_df)
  invisible(TRUE)
}

# ==============================================================================
# Steps 07 & 08: Master View Query (`v_integrative_triage_master`) & Persistence
# ==============================================================================

fetch_run_master_view <- function(conn, run_id = resolve_packera_run_id()) {
  if (is.null(conn)) return(NULL)
  rid <- resolve_packera_run_id(run_id)
  res <- tryCatch(
    DBI::dbGetQuery(
      conn,
      "SELECT * FROM v_integrative_triage_master WHERE run_id = $1 ORDER BY catalogNumber, plant_individual_id;",
      params = list(rid)
    ),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  res$run_id <- NULL
  coerce_sqlite_logical_cols(res)
}

try_fetch_master_view_for_cohort <- function(
  input_paths = character(0),
  active_catalogs = NULL,
  run_id = NULL,
  demo = FALSE,
  vouchers_df = NULL,
  vouchers_path = NULL,
  require_step04 = FALSE,
  require_step05 = FALSE,
  require_step06 = FALSE,
  ...
) {
  eff_paths <- unique(c(input_paths, vouchers_path))
  eff_paths <- eff_paths[!is.na(eff_paths) & nzchar(eff_paths)]
  if (!is_sqlite_active_for_inputs(eff_paths, demo = demo)) return(NULL)
  conn <- connect_packera_db(read_only = TRUE, demo = demo)
  if (is.null(conn)) return(NULL)
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  master_df <- fetch_run_master_view(conn, run_id = rid)
  if (is.null(master_df) || nrow(master_df) == 0L) return(NULL)

  if (is.null(active_catalogs) && !is.null(vouchers_df) && "catalogNumber" %in% names(vouchers_df)) {
    active_catalogs <- vouchers_df$catalogNumber
  }
  if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
    cats <- unique(trimws(as.character(active_catalogs)))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) > 0L) {
      if (!all(cats %in% trimws(as.character(master_df$catalogNumber)))) {
        return(NULL)
      }
      master_df <- master_df[trimws(as.character(master_df$catalogNumber)) %in% cats, , drop = FALSE]
    }
  }
  if (nrow(master_df) == 0L) return(NULL)

  has_m <- ("cda_predicted_taxon" %in% names(master_df) && any(!is.na(master_df$cda_predicted_taxon))) ||
           ("PC1" %in% names(master_df) && any(!is.na(master_df$PC1)))
  has_s <- ("rf_predicted_taxon" %in% names(master_df) && any(!is.na(master_df$rf_predicted_taxon))) ||
           ("edaphic_best_fit_taxon" %in% names(master_df) && any(!is.na(master_df$edaphic_best_fit_taxon)))

  if (isTRUE(require_step04) && !has_m) return(NULL)
  if ((isTRUE(require_step05) || isTRUE(require_step06)) && !has_s) return(NULL)
  master_df
}

persist_step07_to_sqlite <- function(
  triage_df,
  vouchers_path = NULL,
  run_id = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(triage_df) || nrow(triage_df) == 0L) return(invisible(FALSE))
  if (!is_sqlite_active_for_inputs(vouchers_path, demo = demo)) return(invisible(FALSE))

  conn <- connect_packera_db(read_only = FALSE, demo = FALSE)
  if (is.null(conn)) return(invisible(FALSE))
  on.exit(disconnect_packera_db(conn), add = TRUE)

  rid <- resolve_packera_run_id(run_id)
  ensure_pipeline_run_r(conn, run_id = rid)

  ct_df <- as.data.frame(triage_df, stringsAsFactors = FALSE)
  ct_df$run_id <- rid
  if (!("plant_individual_id" %in% names(ct_df))) ct_df$plant_individual_id <- 0L
  ct_df$plant_individual_id[is.na(ct_df$plant_individual_id)] <- 0L

  db_upsert_dataframe(
    conn,
    table_name = "curation_triage",
    df = ct_df,
    pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
    preserve_existing_on_null = FALSE
  )
  .prune_stale_run_rows_r(conn, "curation_triage", rid, ct_df)
  invisible(TRUE)
}

.discover_step08_eco_cols <- function(col_names) {
  if (is.null(col_names) || length(col_names) == 0L) return(character(0))
  eco_pattern <- paste0(
    "^(bio[0-9]+_|fine_|ssurgo_|soil_|topo_|polaris_|solus_|sgmc_|ilr_texture_|elevation_m$|primary_.*_predictors$)",
    "|(^|_)(schoener_d|edaphic_z_score|climatic_z_score|type_edaphic_z_score|nearest_type_niche_taxon|matches_type_lithology|delta_doy|hopkins_|substrate_specialist_flag|rf_prediction_source|n_observed_modalities)"
  )
  setdiff(
    grep(eco_pattern, col_names, value = TRUE, ignore.case = TRUE),
    c("catalogNumber", "plant_individual_id")
  )
}

.discover_step08_dynamic_morph_cols <- function(col_names) {
  if (is.null(col_names) || length(col_names) == 0L) return(character(0))
  morph_pattern <- "^(Sym_PC[0-9]+|PC[0-9]+|Cauline_PC[0-9]+|ECT_PC[0-9]+|Cauline_ECT_PC[0-9]+|[ABCD][0-9]+|efd_harmonic_opt|efd_k95_pcs|selected_morph_features|ect_lobe_complexity|cauline_ect_lobe_complexity|heteroblasty_ect_distance)$"
  grep(morph_pattern, col_names, value = TRUE)
}

.backfill_cols_by_composite_key <- function(target_df, src_df, cols_to_merge) {
  if (is.null(target_df) || nrow(target_df) == 0L || is.null(src_df) || nrow(src_df) == 0L) {
    return(target_df)
  }
  if (!("catalogNumber" %in% names(target_df)) || !("catalogNumber" %in% names(src_df))) {
    return(target_df)
  }
  cols_eff <- intersect(setdiff(cols_to_merge, c("catalogNumber", "plant_individual_id")), names(src_df))
  if (length(cols_eff) == 0L) return(target_df)

  src_df$catalogNumber <- as.character(src_df$catalogNumber)
  has_pid <- ("plant_individual_id" %in% names(target_df)) && ("plant_individual_id" %in% names(src_df))
  idx_pid <- if (has_pid) {
    t_pid <- ifelse(is.na(target_df$plant_individual_id), 0L, as.integer(target_df$plant_individual_id))
    s_pid <- ifelse(is.na(src_df$plant_individual_id), 0L, as.integer(src_df$plant_individual_id))
    match(paste0(target_df$catalogNumber, "__", t_pid), paste0(src_df$catalogNumber, "__", s_pid))
  } else {
    rep(NA_integer_, nrow(target_df))
  }
  src_cat <- src_df[!duplicated(src_df$catalogNumber), , drop = FALSE]
  idx_cat <- match(target_df$catalogNumber, src_cat$catalogNumber)

  for (cl in cols_eff) {
    cur <- if (cl %in% names(target_df)) target_df[[cl]] else rep(NA, nrow(target_df))
    miss_cur <- is.na(cur) | (is.character(cur) & !nzchar(trimws(cur)))
    hp <- miss_cur & !is.na(idx_pid)
    if (any(hp)) cur[hp] <- src_df[[cl]][idx_pid[hp]]
    miss_cur2 <- is.na(cur) | (is.character(cur) & !nzchar(trimws(cur)))
    hc <- miss_cur2 & !is.na(idx_cat)
    if (any(hc)) cur[hc] <- src_cat[[cl]][idx_cat[hc]]
    target_df[[cl]] <- cur
  }
  target_df
}

load_step07_merged_inputs <- function(opts, standardize_fn) {
  in_paths <- c(opts$morphometrics, opts$harmonics, opts$vouchers, opts$multimodal)
  active_cats <- NULL
  if (!is.null(opts$vouchers) && nzchar(opts$vouchers) && file.exists(opts$vouchers)) {
    v_tmp <- tryCatch(utils::read.csv(opts$vouchers, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(v_tmp) && "catalogNumber" %in% names(v_tmp)) {
      active_cats <- unique(trimws(as.character(v_tmp$catalogNumber)))
    }
  }

  master_df <- try_fetch_master_view_for_cohort(
    input_paths = in_paths,
    active_catalogs = active_cats,
    require_step04 = TRUE
  )

  if (!is.null(master_df) && nrow(master_df) > 0L) {
    morph_df <- as.data.frame(master_df, stringsAsFactors = FALSE)
    morph_df$catalogNumber <- as.character(morph_df$catalogNumber)

    # Backfill any dynamic morphometric or ecological columns from disk snapshots if absent/NA in view
    for (supp_m in c(opts$morphometrics, opts$harmonics)) {
      if (!is.null(supp_m) && nzchar(supp_m) && file.exists(supp_m)) {
        sm_df <- tryCatch(utils::read.csv(supp_m, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(sm_df) && nrow(sm_df) > 0L) {
          dyn_m_cols <- union(
            setdiff(names(sm_df), names(morph_df)),
            .discover_step08_dynamic_morph_cols(names(sm_df))
          )
          morph_df <- .backfill_cols_by_composite_key(morph_df, sm_df, dyn_m_cols)
        }
      }
    }
    if (!is.null(opts$multimodal) && nzchar(opts$multimodal) && file.exists(opts$multimodal)) {
      mm_disk <- tryCatch(utils::read.csv(opts$multimodal, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(mm_disk) && nrow(mm_disk) > 0L) {
        eco_cols_disk <- .discover_step08_eco_cols(names(mm_disk))
        morph_df <- .backfill_cols_by_composite_key(morph_df, mm_disk, eco_cols_disk)
      }
    }

    raw_sp <- if ("species_standardized" %in% names(morph_df)) morph_df$species_standardized else if ("scientificName" %in% names(morph_df)) morph_df$scientificName else morph_df$species_raw
    morph_df$species_standardized <- standardize_fn(raw_sp)
    morph_df <- morph_df[!is.na(morph_df$species_standardized) & nzchar(morph_df$species_standardized) & morph_df$species_standardized != "Unknown", , drop = FALSE]
    if ("cda_predicted_taxon" %in% names(morph_df)) morph_df$cda_predicted_taxon <- standardize_fn(morph_df$cda_predicted_taxon)
    return(list(morph_df = morph_df))
  }

  morph_path <- opts$morphometrics
  if (!file.exists(morph_path)) {
    alt_m <- file.path(dirname(morph_path), "morphometrics_misidentification_flags.csv")
    if (file.exists(alt_m)) morph_path <- alt_m
  }
  if (!file.exists(morph_path)) stop("Morphometrics table not found at: ", opts$morphometrics)
  morph_df <- utils::read.csv(morph_path, stringsAsFactors = FALSE)
  morph_df$catalogNumber <- as.character(morph_df$catalogNumber)

  for (supp_path in c(opts$harmonics, opts$vouchers)) {
    if (!is.null(supp_path) && nzchar(supp_path) && file.exists(supp_path)) {
      s_df <- tryCatch(utils::read.csv(supp_path, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(s_df) && nrow(s_df) > 0 && "catalogNumber" %in% names(s_df)) {
        merge_cols <- union(
          setdiff(names(s_df), names(morph_df)),
          .discover_step08_dynamic_morph_cols(names(s_df))
        )
        morph_df <- .backfill_cols_by_composite_key(morph_df, s_df, merge_cols)
      }
    }
  }

  pheno_path <- file.path(dirname(morph_path), "voucher_phenological_states.csv")
  if (file.exists(pheno_path)) {
    ph_df <- tryCatch(utils::read.csv(pheno_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(ph_df) && nrow(ph_df) > 0 && "catalogNumber" %in% names(ph_df)) {
      add_ph <- setdiff(names(ph_df), names(morph_df))
      morph_df <- .backfill_cols_by_composite_key(morph_df, ph_df, add_ph)
    }
  }

  if (!is.null(opts$multimodal) && nzchar(opts$multimodal) && file.exists(opts$multimodal)) {
    mm_df <- tryCatch(utils::read.csv(opts$multimodal, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(mm_df) && nrow(mm_df) > 0 && "catalogNumber" %in% names(mm_df)) {
      eco_cols <- .discover_step08_eco_cols(names(mm_df))
      morph_df <- .backfill_cols_by_composite_key(morph_df, mm_df, eco_cols)
    }
  }

  raw_sp <- if ("species_standardized" %in% names(morph_df)) morph_df$species_standardized else if ("scientificName" %in% names(morph_df)) morph_df$scientificName else if ("species_raw" %in% names(morph_df)) morph_df$species_raw else if ("taxon" %in% names(morph_df)) morph_df$taxon else rep("Unknown", nrow(morph_df))
  morph_df$species_standardized <- standardize_fn(raw_sp)
  morph_df <- morph_df[!is.na(morph_df$species_standardized) & nzchar(morph_df$species_standardized) & morph_df$species_standardized != "Unknown", , drop = FALSE]
  if ("cda_predicted_taxon" %in% names(morph_df)) morph_df$cda_predicted_taxon <- standardize_fn(morph_df$cda_predicted_taxon)

  list(morph_df = morph_df)
}
load_step08_merged_inputs <- load_step07_merged_inputs
