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
#     - Single-query 120-vertex Tier_1_Direct contour ingestion for Step 03.
#     - Step 04 OTU + plant_individuals + reproductive SQLite joins & updates.
#     - Step 06 environmental_cache + environmental_rtree spatial bbox queries.
#     - Step 07 & Step 08 `fetch_run_master_view()` (`v_integrative_triage_master`)
#       with seamless CSV fallback merging.
# ==============================================================================

PACKERA_DB_LOGICAL_COLS <- c(
  "is_primary_duplicate", "is_calibrated_mm", "is_bent_stem", "is_folded",
  "is_dissected", "reflection_applied", "is_verified_anthesis",
  "is_artifact_driven", "is_background_shortcut", "is_rock_outcrop",
  "is_azonal_lithology", "high_coord_uncertainty_flag",
  "is_training_eligible_coord", "azonal_lithology_presence_buffer",
  "is_promoted_bronze", "is_promoted_anchor", "is_passive",
  "is_contour_outlier", "morphometric_outlier_flag", "misidentification_flag",
  "is_cleanlab_issue", "is_label_corrupted", "vision_is_ood",
  "tomentum_taxon_discordance", "substrate_specialist_flag",
  "piedmont_outcrop_specialist_flag", "multimodal_conflict_flag",
  "is_completed_voucher"
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

.coalesce_wc_bio_qt_qtr_df <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  pairs <- list(
    c("bio18_precip_warmest_qtr", "bio18_precip_warmest_qt"),
    c("bio19_precip_coldest_qtr", "bio19_precip_coldest_qt")
  )
  for (pr in pairs) {
    qtr_c <- pr[1]
    qt_c <- pr[2]
    has_qtr <- qtr_c %in% names(df)
    has_qt <- qt_c %in% names(df)
    if (!has_qtr && !has_qt) next
    v_qtr <- if (has_qtr) suppressWarnings(as.numeric(df[[qtr_c]])) else rep(NA_real_, nrow(df))
    v_qt <- if (has_qt) suppressWarnings(as.numeric(df[[qt_c]])) else rep(NA_real_, nrow(df))
    v_qtr[!is.finite(v_qtr)] <- NA_real_
    v_qt[!is.finite(v_qt)] <- NA_real_
    df[[qtr_c]] <- ifelse(!is.na(v_qtr), v_qtr, v_qt)
    df[[qt_c]] <- ifelse(!is.na(v_qt), v_qt, v_qtr)
  }
  df
}

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
  if (identical(table_name, "leaf_efa_harmonics")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_HARMONIC_COLS)), PACKERA_DB_HARMONIC_COLS),
      stats::setNames(rep("REAL", length(PACKERA_DB_PC_COLS)), PACKERA_DB_PC_COLS),
      stats::setNames(rep("REAL", length(PACKERA_DB_SYM_PC_COLS)), PACKERA_DB_SYM_PC_COLS),
      c(efd_harmonic_opt = "INTEGER", efd_k95_pcs = "INTEGER", dissection_index = "REAL")
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
        basal_leaf_count = "INTEGER"
      )
    ))
  }
  if (identical(table_name, "spatial_rf_predictions")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_WC_BIO_COLS)), PACKERA_DB_WC_BIO_COLS),
      c(
        bio18_precip_warmest_qt = "REAL",
        bio19_precip_coldest_qt = "REAL",
        climatic_z_score = "REAL",
        primary_env_predictors = "TEXT",
        primary_clim_predictors = "TEXT",
        primary_edaph_predictors = "TEXT"
      )
    ))
  }
  if (identical(table_name, "environmental_cache")) {
    return(c(
      stats::setNames(rep("REAL", length(PACKERA_DB_WC_BIO_COLS)), PACKERA_DB_WC_BIO_COLS),
      c(
        bio18_precip_warmest_qt = "REAL",
        bio19_precip_coldest_qt = "REAL"
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

  updated_fields <- tryCatch(DBI::dbListFields(conn, table_name), error = function(e) existing)
  if (all(c("bio18_precip_warmest_qtr", "bio18_precip_warmest_qt",
            "bio19_precip_coldest_qtr", "bio19_precip_coldest_qt") %in% updated_fields)) {
    tryCatch(
      DBI::dbExecute(
        conn,
        sprintf(
          paste0(
            "UPDATE %s SET ",
            "bio18_precip_warmest_qtr = COALESCE(bio18_precip_warmest_qtr, bio18_precip_warmest_qt), ",
            "bio18_precip_warmest_qt = COALESCE(bio18_precip_warmest_qt, bio18_precip_warmest_qtr), ",
            "bio19_precip_coldest_qtr = COALESCE(bio19_precip_coldest_qtr, bio19_precip_coldest_qt), ",
            "bio19_precip_coldest_qt = COALESCE(bio19_precip_coldest_qt, bio19_precip_coldest_qtr) ",
            "WHERE (bio18_precip_warmest_qtr IS NULL AND bio18_precip_warmest_qt IS NOT NULL) ",
            "   OR (bio18_precip_warmest_qt IS NULL AND bio18_precip_warmest_qtr IS NOT NULL) ",
            "   OR (bio19_precip_coldest_qtr IS NULL AND bio19_precip_coldest_qt IS NOT NULL) ",
            "   OR (bio19_precip_coldest_qt IS NULL AND bio19_precip_coldest_qtr IS NOT NULL);"
          ),
          q_tbl
        )
      ),
      error = function(e) NULL
    )
  }
  added_any
}

.migrate_dynamic_sqlite_schema_r <- function(conn) {
  if (is.null(conn)) return(FALSE)
  changed <- FALSE
  for (tbl in c("leaf_efa_harmonics", "otu_morphometrics", "spatial_rf_predictions", "environmental_cache")) {
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
  py_bin <- if (file.exists(".venv/bin/python")) ".venv/bin/python" else Sys.which("python3")
  if (!nzchar(py_bin)) return(invisible(FALSE))
  tryCatch({
    system2(
      py_bin,
      args = c(
        "-c",
        "import sys; from scripts.core.db import get_db_connection; c = get_db_connection(sys.argv[1], init_schema=True); c.close()",
        eff_path
      ),
      stdout = FALSE,
      stderr = FALSE
    )
    TRUE
  }, error = function(e) FALSE)
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
    tryCatch({
      DBI::dbExecute(conn, "PRAGMA busy_timeout = 30000;")
      DBI::dbExecute(conn, "PRAGMA foreign_keys = ON;")
    }, error = function(e) NULL)
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
  tryCatch({
    DBI::dbExecute(conn, "PRAGMA journal_mode = WAL;")
    DBI::dbExecute(conn, "PRAGMA synchronous = NORMAL;")
    DBI::dbExecute(conn, "PRAGMA foreign_keys = ON;")
    DBI::dbExecute(conn, "PRAGMA busy_timeout = 30000;")
  }, error = function(e) NULL)

  schema_migrated <- .migrate_dynamic_sqlite_schema_r(conn)
  if (!identical(eff_path, ":memory:")) {
    view_needs_refresh <- schema_migrated
    if (!view_needs_refresh && db_table_exists_safe(conn, "v_integrative_triage_master")) {
      v_cols <- tryCatch(DBI::dbListFields(conn, "v_integrative_triage_master"), error = function(e) character(0))
      if (length(v_cols) > 0L && !all(c("PC15", "Sym_PC15", "selected_morph_features", "climatic_z_score", "bio19_precip_coldest_qtr") %in% v_cols)) {
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
        tryCatch({
          DBI::dbExecute(conn, "PRAGMA journal_mode = WAL;")
          DBI::dbExecute(conn, "PRAGMA synchronous = NORMAL;")
          DBI::dbExecute(conn, "PRAGMA foreign_keys = ON;")
          DBI::dbExecute(conn, "PRAGMA busy_timeout = 30000;")
        }, error = function(e) NULL)
      }
    }
  }
  conn
}

disconnect_packera_db <- function(conn) {
  if (!is.null(conn)) {
    tryCatch(DBI::dbDisconnect(conn), error = function(e) NULL)
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

  DBI::dbWithTransaction(conn, {
    DBI::dbExecute(conn, sql, params = unname(as.list(sub_df)))
  })
  nrow(sub_df)
}

load_taxon_synonyms_sqlite <- function(conn = NULL) {
  own_conn <- is.null(conn)
  if (own_conn) conn <- connect_packera_db(read_only = TRUE)
  if (is.null(conn)) return(character(0))
  on.exit(if (own_conn) disconnect_packera_db(conn), add = TRUE)
  if (!db_table_exists_safe(conn, "taxon_synonyms")) return(character(0))
  res <- tryCatch(
    DBI::dbGetQuery(conn, "SELECT raw_name_lower, canonical_packera FROM taxon_synonyms;"),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0L) return(character(0))
  stats::setNames(as.character(res$canonical_packera), tolower(trimws(as.character(res$raw_name_lower))))
}

# ==============================================================================
# Step 03: Single-Query 120-Vertex Contour & Leaf Metadata Ingestion + Persistence
# ==============================================================================

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

  base_sql <- paste0(
    "SELECT ",
    "  el.catalogNumber, el.plant_individual_id, el.leaf_id, el.shape_id, ",
    "  el.assigned_tier, COALESCE(el.tier_assignment, el.assigned_tier) AS tier_assignment, ",
    "  el.solidity, COALESCE(el.ucs_score, el.ucs) AS ucs_score, ",
    "  el.aspect_ratio, el.area_px, el.blade_area_cm2, el.is_dissected, ",
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

  if (!is.null(active_catalogs) && length(active_catalogs) > 0L) {
    cats <- unique(trimws(as.character(active_catalogs)))
    cats <- cats[!is.na(cats) & nzchar(cats)]
    if (length(cats) == 0L) return(NULL)
    q_cats <- paste(vapply(cats, function(x) as.character(DBI::dbQuoteString(conn, x)), character(1)), collapse = ", ")
    sql <- paste0(base_sql, sprintf("WHERE el.catalogNumber IN (%s) ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;", q_cats))
    rows <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  } else if (has_cohort) {
    sql <- paste0(
      base_sql,
      "INNER JOIN run_active_cohort rac ON rac.catalogNumber = el.catalogNumber AND rac.plant_individual_id = el.plant_individual_id ",
      "WHERE rac.run_id = $1 ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;"
    )
    rows <- tryCatch(DBI::dbGetQuery(conn, sql, params = list(rid)), error = function(e) NULL)
  } else {
    sql <- paste0(base_sql, "ORDER BY el.catalogNumber, el.plant_individual_id, el.leaf_id;")
    rows <- tryCatch(DBI::dbGetQuery(conn, sql), error = function(e) NULL)
  }

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
    m_hdr <- tryCatch(readr::read_csv(eff_mpath, show_col_types = FALSE, progress = FALSE), error = function(e) NULL)
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

  coords_list <- vector("list", nrow(rows))
  coord_names <- character(nrow(rows))
  seen_sids <- new.env(hash = TRUE, parent = emptyenv())
  meta_rows <- vector("list", nrow(rows))
  kept_k <- 0L

  has_momocs <- requireNamespace("Momocs", quietly = TRUE)
  for (i in seq_len(nrow(rows))) {
    raw_tier <- as.character(ifelse(!is.na(rows$tier_assignment[i]) && nzchar(rows$tier_assignment[i]), rows$tier_assignment[i], rows$assigned_tier[i]))
    resolved_tier <- ifelse(grepl("^tier_?1", raw_tier, ignore.case = TRUE), "Tier_1_Direct", ifelse(grepl("^tier_?2", raw_tier, ignore.case = TRUE), "Tier_2_Reflected", raw_tier))
    if (is.na(resolved_tier) || resolved_tier != "Tier_1_Direct") next
    js <- rows$contour_xy_json[i]
    if (is.na(js) || !nzchar(js)) next
    coo <- tryCatch(jsonlite::fromJSON(js, simplifyMatrix = TRUE), error = function(e) NULL)
    if (is.null(coo) || !is.matrix(coo) || ncol(coo) < 2L || nrow(coo) < 10L) next
    coo <- matrix(as.numeric(coo[, 1:2]), ncol = 2L)
    if (any(!is.finite(coo))) next
    colnames(coo) <- c("x", "y")
    if (has_momocs) coo <- tryCatch(Momocs::coo_cw(coo), error = function(e) coo)

    x_next <- c(coo[-1L, 1L], coo[1L, 1L])
    y_next <- c(coo[-1L, 2L], coo[1L, 2L])
    area_px <- 0.5 * abs(sum(coo[, 1L] * y_next - x_next * coo[, 2L]))
    perim_px <- sum(sqrt((x_next - coo[, 1L])^2 + (y_next - coo[, 2L])^2))
    circ_val <- if (is.finite(perim_px) && perim_px > 1e-6) round(min(1.0, (4 * pi * area_px) / (perim_px^2)), 4) else NA_real_
    di_val <- if (is.finite(perim_px) && is.finite(area_px) && area_px > 1e-6) round(max(1.0, perim_px / (2.0 * sqrt(pi * area_px))), 4) else NA_real_

    resolved_ppm <- suppressWarnings(as.numeric(rows$pixels_per_mm[i]))
    is_cal_raw <- rows$is_calibrated_mm[i]
    is_cal <- (!is.na(is_cal_raw) && as.integer(is_cal_raw) != 0L)
    if (!is.na(resolved_ppm) && is.finite(resolved_ppm) && resolved_ppm > 0) {
      is_cal <- TRUE
    } else {
      is_cal <- FALSE
      resolved_ppm <- NA_real_
    }

    if (is_cal) coo <- coo / resolved_ppm
    cs_val <- if (!is.null(centsize_fn)) {
      centsize_fn(coo)
    } else {
      cx <- mean(coo[, 1L]); cy <- mean(coo[, 2L])
      sqrt(sum((coo[, 1L] - cx)^2 + (coo[, 2L] - cy)^2))
    }

    if (!is.null(align_fn)) {
      coo <- align_fn(coo, circularity = circ_val, has_momocs = has_momocs)
    }

    lw <- if (has_momocs) tryCatch(Momocs::coo_lw(coo), error = function(e) NULL) else NULL
    if (is.null(lw) || length(lw) < 2L || is.na(lw[2]) || lw[2] <= 1e-6) {
      coo_centered <- scale(coo, center = TRUE, scale = FALSE)
      svd_rot <- coo_centered %*% svd(coo_centered)$v
      lw <- c(diff(range(svd_rot[, 1], na.rm = TRUE)), diff(range(svd_rot[, 2], na.rm = TRUE)))
    }
    ar_val <- round(as.numeric(lw[1] / lw[2]), 4)

    cat_num <- as.character(rows$catalogNumber[i])
    pid_val <- as.integer(ifelse(is.na(rows$plant_individual_id[i]), 0L, rows$plant_individual_id[i]))
    lid_val <- as.integer(ifelse(is.na(rows$leaf_id[i]), i, rows$leaf_id[i]))
    raw_sid <- as.character(rows$shape_id[i])
    shape_id <- if (!is.na(raw_sid) && nzchar(raw_sid)) raw_sid else sprintf("%s_p%d_leaf%d", cat_num, pid_val, lid_val)
    if (exists(shape_id, envir = seen_sids, inherits = FALSE)) {
      shape_id <- sprintf("%s_%d", shape_id, i)
    }
    assign(shape_id, TRUE, envir = seen_sids)

    db_area <- suppressWarnings(as.numeric(rows$area_px[i]))
    eff_area <- if (!is.na(db_area) && is.finite(db_area) && db_area > 0) db_area else round(area_px, 1)
    diss_raw <- rows$is_dissected[i]
    eff_diss <- if (is.na(diss_raw)) NA else (as.integer(diss_raw) != 0L)
    sol_val <- suppressWarnings(as.numeric(rows$solidity[i]))
    if (is.na(sol_val) || !is.finite(sol_val)) {
      ch_idx <- tryCatch(grDevices::chull(coo[, 1L], coo[, 2L]), error = function(e) integer(0))
      if (length(ch_idx) >= 3L) {
        hc <- coo[ch_idx, , drop = FALSE]
        h_area <- 0.5 * abs(sum(hc[, 1L] * c(hc[-1L, 2L], hc[1L, 2L]) - c(hc[-1L, 1L], hc[1L, 1L]) * hc[, 2L]))
        if (h_area > 1e-6) sol_val <- round(min(1.0, area_px / h_area), 4)
      }
    }

    kept_k <- kept_k + 1L
    coords_list[[kept_k]] <- coo
    coord_names[kept_k] <- shape_id
    m_hit <- NA_integer_
    if (!is.null(m_hdr)) {
      if (!is.null(m_idx_3k)) m_hit <- unname(m_idx_3k[paste(trimws(cat_num), pid_val, lid_val, sep = "__")])
      if (is.na(m_hit) && !is.null(m_idx_2k)) m_hit <- unname(m_idx_2k[paste(trimws(cat_num), lid_val, sep = "__")])
    }
    m_get_num <- function(cn, fb = NA_real_) {
      if (!is.na(m_hit) && !is.null(m_hdr) && cn %in% names(m_hdr)) {
        v <- suppressWarnings(as.numeric(m_hdr[[cn]][m_hit]))
        if (is.finite(v)) return(v)
      }
      fb
    }
    db_ba_cm2 <- if ("blade_area_cm2" %in% names(rows)) suppressWarnings(as.numeric(rows$blade_area_cm2[i])) else NA_real_
    if (!is.finite(db_ba_cm2)) db_ba_cm2 <- NA_real_
    meta_rows[[kept_k]] <- data.frame(
      shape_id = shape_id,
      catalogNumber = cat_num,
      plant_individual_id = pid_val,
      leaf_id = lid_val,
      tier_assignment = resolved_tier,
      assigned_tier = resolved_tier,
      reconstruction_tier = resolved_tier,
      solidity = sol_val,
      dissection_index = as.numeric(di_val),
      ucs_score = suppressWarnings(as.numeric(rows$ucs_score[i])),
      pixels_per_mm = resolved_ppm,
      is_calibrated_mm = is_cal,
      petiole_length_mm = suppressWarnings(as.numeric(rows$petiole_length_mm[i])),
      blade_length_mm = suppressWarnings(as.numeric(rows$blade_length_mm[i])),
      petiole_to_blade_ratio = suppressWarnings(as.numeric(rows$petiole_to_blade_ratio[i])),
      is_dissected = eff_diss,
      cy_rel_plant_bbox = m_get_num("cy_rel_plant_bbox", NA_real_),
      cy_rel_sheet = m_get_num("cy_rel_sheet", NA_real_),
      blade_area_cm2 = m_get_num("blade_area_cm2", db_ba_cm2),
      effective_min_blade_area_cm2 = m_get_num("effective_min_blade_area_cm2", NA_real_),
      AR_inv = suppressWarnings(as.numeric(rows$AR_inv[i])),
      capitula_count = suppressWarnings(as.numeric(rows$capitula_count[i])),
      involucre_height_px = suppressWarnings(as.numeric(rows$involucre_height_px[i])),
      involucre_width_px = suppressWarnings(as.numeric(rows$involucre_width_px[i])),
      involucre_height_mm = suppressWarnings(as.numeric(rows$involucre_height_mm[i])),
      involucre_width_mm = suppressWarnings(as.numeric(rows$involucre_width_mm[i])),
      centroid_size = as.numeric(cs_val),
      log_centsize = if (!is.na(cs_val) && cs_val > 0) log(cs_val) else NA_real_,
      aspect_ratio = as.numeric(ar_val),
      area_px = as.numeric(eff_area),
      circularity = as.numeric(circ_val),
      mask_source = paste0("sqlite://extracted_leaves/", shape_id),
      stringsAsFactors = FALSE
    )
  }

  if (kept_k == 0L) return(NULL)
  coords_list <- coords_list[seq_len(kept_k)]
  names(coords_list) <- coord_names[seq_len(kept_k)]
  meta_df <- as.data.frame(dplyr::bind_rows(meta_rows[seq_len(kept_k)]), stringsAsFactors = FALSE)
  rownames(meta_df) <- coord_names[seq_len(kept_k)]
  list(
    coords_list = coords_list,
    meta_df = meta_df
  )
}

.normalize_morph_dynamic_cols_r <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  for (ic in intersect(c("plant_individual_id", "leaf_id", "efd_harmonic_opt", "efd_k95_pcs", "basal_leaf_count"), names(df))) {
    df[[ic]] <- suppressWarnings(as.integer(df[[ic]]))
  }
  if ("selected_morph_features" %in% names(df)) {
    df[["selected_morph_features"]] <- .normalize_semicolon_col_r(df[["selected_morph_features"]])
  }
  dyn_num_cols <- intersect(
    c(PACKERA_DB_HARMONIC_COLS, PACKERA_DB_PC_COLS, PACKERA_DB_SYM_PC_COLS,
      "dissection_index", "q25_solidity", "dissected_leaf_fraction"),
    names(df)
  )
  for (nc in dyn_num_cols) {
    df[[nc]] <- suppressWarnings(as.numeric(df[[nc]]))
  }
  df
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
    db_upsert_dataframe(
      conn,
      table_name = "otu_morphometrics",
      df = otu_out,
      pk_cols = c("run_id", "catalogNumber", "plant_individual_id"),
      preserve_existing_on_null = TRUE
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
          "determiner_tier", "is_promoted_bronze", "authority_weight"),
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
# Step 06: Environmental Cache + R*Tree Spatial Index & Spatial RF Predictions
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
    return(.coalesce_wc_bio_qt_qtr_df(coerce_sqlite_logical_cols(res)))
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
  .coalesce_wc_bio_qt_qtr_df(coerce_sqlite_logical_cols(res))
}

load_environmental_cache_sqlite <- function(
  cache_path = NULL,
  db_path = NULL,
  catalog_numbers = NULL,
  lat_lon_keys = NULL,
  bbox = NULL,
  demo = FALSE
) {
  query_environmental_cache_sqlite(
    catalog_numbers = catalog_numbers,
    lat_lon_keys = lat_lon_keys,
    bbox = bbox,
    cache_path = cache_path,
    db_path = db_path,
    demo = demo
  )
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

  env_df <- .coalesce_wc_bio_qt_qtr_df(as.data.frame(df, stringsAsFactors = FALSE))
  for (bc in intersect(c(PACKERA_DB_WC_BIO_COLS, "bio18_precip_warmest_qt", "bio19_precip_coldest_qt"), names(env_df))) {
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

persist_step06_to_sqlite <- function(
  vouchers_df,
  vouchers_path = NULL,
  run_id = NULL,
  demo = FALSE
) {
  if (isTRUE(demo) || is.null(vouchers_df) || nrow(vouchers_df) == 0L) return(invisible(FALSE))
  if (!is_sqlite_active_for_inputs(vouchers_path, demo = demo)) return(invisible(FALSE))

  srf_df <- .coalesce_wc_bio_qt_qtr_df(as.data.frame(vouchers_df, stringsAsFactors = FALSE))
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
  for (nc in intersect(c("edaphic_z_score", "climatic_z_score", PACKERA_DB_WC_BIO_COLS, "bio18_precip_warmest_qt", "bio19_precip_coldest_qt"), names(srf_df))) {
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
  .coalesce_wc_bio_qt_qtr_df(coerce_sqlite_logical_cols(res))
}

try_fetch_master_view_for_cohort <- function(
  input_paths = character(0),
  active_catalogs = NULL,
  require_morph_or_vision = FALSE,
  run_id = NULL,
  demo = FALSE,
  vouchers_df = NULL,
  vouchers_path = NULL,
  require_step04 = FALSE,
  require_step05 = FALSE,
  require_step06 = FALSE
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
  has_v <- ("vision_predicted_label" %in% names(master_df) && any(!is.na(master_df$vision_predicted_label)))
  has_s <- ("rf_predicted_taxon" %in% names(master_df) && any(!is.na(master_df$rf_predicted_taxon))) ||
           ("edaphic_best_fit_taxon" %in% names(master_df) && any(!is.na(master_df$edaphic_best_fit_taxon)))

  if (isTRUE(require_step04) && !has_m) return(NULL)
  if (isTRUE(require_step05) && !has_v) return(NULL)
  if (isTRUE(require_step06) && !has_s) return(NULL)
  if (isTRUE(require_morph_or_vision) && !has_m && !has_v) return(NULL)
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
  invisible(TRUE)
}

.discover_step08_eco_cols <- function(col_names) {
  if (is.null(col_names) || length(col_names) == 0L) return(character(0))
  eco_pattern <- paste0(
    "^(bio[0-9]+_|fine_|ssurgo_|soil_|topo_|polaris_|solus_|sgmc_|primary_.*_predictors$)",
    "|(^|_)(schoener_d|edaphic_z_score|climatic_z_score|delta_doy|hopkins_)"
  )
  setdiff(
    grep(eco_pattern, col_names, value = TRUE, ignore.case = TRUE),
    c("catalogNumber", "plant_individual_id")
  )
}

.discover_step08_dynamic_morph_cols <- function(col_names) {
  if (is.null(col_names) || length(col_names) == 0L) return(character(0))
  morph_pattern <- "^(Sym_PC[0-9]+|PC[0-9]+|[ABCD][0-9]+|efd_harmonic_opt|efd_k95_pcs|selected_morph_features)$"
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

load_step08_merged_inputs <- function(opts, standardize_fn) {
  in_paths <- c(opts$morphometrics, opts$harmonics, opts$vision_audit, opts$vouchers, opts$multimodal)
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
    morph_df <- .coalesce_wc_bio_qt_qtr_df(as.data.frame(master_df, stringsAsFactors = FALSE))
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
        mm_disk <- .coalesce_wc_bio_qt_qtr_df(mm_disk)
        eco_cols_disk <- .discover_step08_eco_cols(names(mm_disk))
        morph_df <- .backfill_cols_by_composite_key(morph_df, mm_disk, eco_cols_disk)
        morph_df <- .coalesce_wc_bio_qt_qtr_df(morph_df)
      }
    }

    raw_sp <- if ("species_standardized" %in% names(morph_df)) morph_df$species_standardized else if ("scientificName" %in% names(morph_df)) morph_df$scientificName else morph_df$species_raw
    morph_df$species_standardized <- standardize_fn(raw_sp)
    morph_df <- morph_df[!is.na(morph_df$species_standardized) & nzchar(morph_df$species_standardized) & morph_df$species_standardized != "Unknown", , drop = FALSE]
    if ("cda_predicted_taxon" %in% names(morph_df)) morph_df$cda_predicted_taxon <- standardize_fn(morph_df$cda_predicted_taxon)
    vis_df <- NULL
    if ("vision_predicted_label" %in% names(morph_df) && any(!is.na(morph_df$vision_predicted_label) & nzchar(as.character(morph_df$vision_predicted_label)))) {
      vis_df <- morph_df[!is.na(morph_df$vision_predicted_label) & nzchar(as.character(morph_df$vision_predicted_label)), , drop = FALSE]
      vis_df$vision_predicted_label <- standardize_fn(vis_df$vision_predicted_label)
    } else if (!is.null(opts$vision_audit) && nzchar(opts$vision_audit) && file.exists(opts$vision_audit)) {
      vis_df <- tryCatch(utils::read.csv(opts$vision_audit, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(vis_df) && nrow(vis_df) > 0) {
        vis_sp <- if ("species_standardized" %in% names(vis_df)) vis_df$species_standardized else if ("given_label" %in% names(vis_df)) vis_df$given_label else if ("taxon" %in% names(vis_df)) vis_df$taxon else vis_df$species_raw
        vis_df$species_standardized <- standardize_fn(vis_sp)
        pred_col <- if ("vision_predicted_label" %in% names(vis_df)) vis_df$vision_predicted_label else if ("predicted_label" %in% names(vis_df)) vis_df$predicted_label else vis_df$species_standardized
        vis_df$vision_predicted_label <- standardize_fn(pred_col)
      }
    }
    return(list(morph_df = morph_df, vis_df = vis_df))
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
      mm_df <- .coalesce_wc_bio_qt_qtr_df(mm_df)
      eco_cols <- .discover_step08_eco_cols(names(mm_df))
      morph_df <- .backfill_cols_by_composite_key(morph_df, mm_df, eco_cols)
      morph_df <- .coalesce_wc_bio_qt_qtr_df(morph_df)
    }
  }

  raw_sp <- if ("species_standardized" %in% names(morph_df)) morph_df$species_standardized else if ("scientificName" %in% names(morph_df)) morph_df$scientificName else if ("species_raw" %in% names(morph_df)) morph_df$species_raw else if ("taxon" %in% names(morph_df)) morph_df$taxon else rep("Unknown", nrow(morph_df))
  morph_df$species_standardized <- standardize_fn(raw_sp)
  morph_df <- morph_df[!is.na(morph_df$species_standardized) & nzchar(morph_df$species_standardized) & morph_df$species_standardized != "Unknown", , drop = FALSE]
  if ("cda_predicted_taxon" %in% names(morph_df)) morph_df$cda_predicted_taxon <- standardize_fn(morph_df$cda_predicted_taxon)

  vis_df <- NULL
  if (!is.null(opts$vision_audit) && nzchar(opts$vision_audit) && file.exists(opts$vision_audit)) {
    vis_df <- tryCatch(utils::read.csv(opts$vision_audit, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(vis_df) && nrow(vis_df) > 0) {
      vis_sp <- if ("species_standardized" %in% names(vis_df)) vis_df$species_standardized else if ("given_label" %in% names(vis_df)) vis_df$given_label else if ("taxon" %in% names(vis_df)) vis_df$taxon else vis_df$species_raw
      vis_df$species_standardized <- standardize_fn(vis_sp)
      pred_col <- if ("vision_predicted_label" %in% names(vis_df)) vis_df$vision_predicted_label else if ("predicted_label" %in% names(vis_df)) vis_df$predicted_label else vis_df$species_standardized
      vis_df$vision_predicted_label <- standardize_fn(pred_col)
    }
  }
  list(morph_df = morph_df, vis_df = vis_df)
}
