# ==============================================================================
# Module: scripts/morphometrics/step03_config_scale_utils.R
# Step 03 Config Loading, Calibration Provenance, Photometry Discovery & CLI
# ==============================================================================

STEP03_UNCALIBRATED_SOURCES <- c(
  "mp_estimated", "predicted_from_megapixels", "tick_lattice_unboxed", "uncalibrated"
)

.step03_coerce_bool <- function(v) {
  if (is.null(v) || length(v) == 0L) return(NA)
  val <- v[1L]
  if (is.logical(val)) return(val)
  if (is.numeric(val)) return(if (!is.finite(val)) NA else (val != 0))
  s <- tolower(trimws(as.character(val)))
  if (s %in% c("true", "t", "1", "1.0", "yes", "y")) return(TRUE)
  if (s %in% c("false", "f", "0", "0.0", "no", "n")) return(FALSE)
  NA
}

.step03_coerce_int <- function(v, min_val = 1L) {
  if (is.null(v) || length(v) == 0L || is.logical(v[1L])) return(NULL)
  val <- suppressWarnings(as.numeric(v[1L]))
  if (!is.finite(val)) return(NULL)
  iv <- as.integer(round(val))
  if (!is.na(min_val) && iv < as.integer(min_val)) return(NULL)
  iv
}

.step03_coerce_prob <- function(v) {
  if (is.null(v) || length(v) == 0L || is.logical(v[1L])) return(NULL)
  val <- suppressWarnings(as.numeric(v[1L]))
  if (!is.finite(val) || val <= 0.0 || val > 1.0) return(NULL)
  as.numeric(val)
}

.step03_parse_yaml_scalar <- function(raw_val) {
  s <- trimws(raw_val)
  if (!nzchar(s) || s == "{}") return(NULL)
  if (grepl("^[\"'].*[\"']$", s)) s <- substr(s, 2L, nchar(s) - 1L)
  b <- .step03_coerce_bool(s)
  if (!is.na(b)) return(b)
  num <- suppressWarnings(as.numeric(s))
  if (is.finite(num)) return(num)
  s
}

.step03_parse_morph_yaml_fallback <- function(path) {
  lines <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(e) character(0))
  m_out <- list(); f_out <- list()
  in_morph <- FALSE; morph_ind <- NA_integer_; subsec <- ""; sub_ind <- NA_integer_
  for (ln in lines) {
    if (grepl("^[[:space:]]*(#|$)", ln)) next
    clean_ln <- sub("[[:space:]]+#.*$", "", ln)
    if (!grepl("^[[:space:]]*[A-Za-z0-9_]+:", clean_ln)) next
    ind <- nchar(sub("^([[:space:]]*).*$", "\\1", clean_ln))
    if (ind == 0L) {
      key0 <- trimws(sub(":.*$", "", clean_ln))
      if (key0 == "morphometrics") { in_morph <- TRUE; next }
      if (in_morph) break
      next
    }
    if (!in_morph) next
    if (is.na(morph_ind)) morph_ind <- ind
    key <- trimws(sub(":.*$", "", clean_ln))
    rest <- trimws(sub("^[^:]+:", "", clean_ln))
    if (!nzchar(rest)) {
      subsec <- key; sub_ind <- ind; next
    }
    if (ind <= morph_ind || (!is.na(sub_ind) && ind <= sub_ind)) subsec <- ""
    val <- .step03_parse_yaml_scalar(rest)
    if (is.null(val)) next
    if (!nzchar(subsec)) m_out[[key]] <- val else if (subsec == "fourier") f_out[[key]] <- val
  }
  m_out$fourier <- f_out
  list(morphometrics = m_out)
}

.step03_apply_cfg_map <- function(cfg, src) {
  if (!is.list(src) || length(src) == 0L) return(cfg)
  if (is.null(src$harmonics) && !is.null(src$nb_harmonics)) src$harmonics <- src$nb_harmonics
  for (bk in c("cohort_adaptive_mode", "symmetrize_efa_harmonics")) {
    bv <- .step03_coerce_bool(src[[bk]])
    if (!is.na(bv)) cfg[[bk]] <- bv
  }
  for (ik in c("min_anchor_n_for_empirical", "harmonics", "min_harmonics", "max_harmonics",
               "num_pcs", "min_export_pcs", "max_export_pcs")) {
    iv <- .step03_coerce_int(src[[ik]], min_val = 1L)
    if (!is.null(iv)) cfg[[ik]] <- iv
  }
  rsv <- .step03_coerce_int(src$random_seed, min_val = NA_integer_)
  if (is.null(rsv)) rsv <- .step03_coerce_int(src$seed, min_val = NA_integer_)
  if (!is.null(rsv)) cfg$random_seed <- rsv
  for (pk in c("target_harmonic_power", "target_pca_variance")) {
    pv <- .step03_coerce_prob(src[[pk]])
    if (!is.null(pv)) cfg[[pk]] <- pv
  }
  cfg
}

load_step03_fourier_config <- function(cfg_path = "config/config.yaml", cli_opts = list()) {
  out <- list(
    cohort_adaptive_mode = TRUE, min_anchor_n_for_empirical = 5L,
    harmonics = 12L, min_harmonics = 6L, max_harmonics = 24L,
    target_harmonic_power = 0.99, target_pca_variance = 0.95,
    num_pcs = 5L, min_export_pcs = 5L, max_export_pcs = 15L,
    symmetrize_efa_harmonics = TRUE, random_seed = 42L
  )
  eff_cfg <- if (!is.null(cfg_path) && length(cfg_path) >= 1L && !is.na(cfg_path[1L]) && nzchar(as.character(cfg_path[1L]))) {
    as.character(cfg_path[1L])
  } else {
    e1 <- Sys.getenv("PLANTMATICS_CONFIG_PATH", unset = "")
    e2 <- Sys.getenv("PACKERA_CONFIG_PATH", unset = "")
    if (nzchar(e1)) e1 else if (nzchar(e2)) e2 else "config/config.yaml"
  }
  if (nzchar(eff_cfg) && file.exists(eff_cfg)) {
    raw_yaml <- if (requireNamespace("yaml", quietly = TRUE)) {
      tryCatch(yaml::read_yaml(eff_cfg), error = function(e) NULL)
    } else {
      .step03_parse_morph_yaml_fallback(eff_cfg)
    }
    if (is.list(raw_yaml) && is.list(raw_yaml$morphometrics)) {
      out <- .step03_apply_cfg_map(out, raw_yaml$morphometrics)
      if (is.list(raw_yaml$morphometrics$fourier)) {
        out <- .step03_apply_cfg_map(out, raw_yaml$morphometrics$fourier)
      }
    }
  }
  if (is.list(cli_opts) && length(cli_opts) > 0L) {
    out <- .step03_apply_cfg_map(out, cli_opts)
  }
  out$max_harmonics <- max(out$min_harmonics, out$max_harmonics)
  out$max_export_pcs <- max(out$min_export_pcs, out$max_export_pcs)
  out
}

resolve_step03_calibration_provenance <- function(pixels_per_mm, is_calibrated_mm = NULL, scale_source = NULL) {
  n <- max(length(pixels_per_mm), length(is_calibrated_mm), length(scale_source))
  if (n == 0L) return(list(is_calibrated_mm = logical(0), scale_source = character(0)))
  recycle_v <- function(x, fill) if (is.null(x) || length(x) == 0L) rep(fill, n) else rep_len(x, n)
  ppm <- suppressWarnings(as.numeric(recycle_v(pixels_per_mm, NA_real_)))
  ss_raw <- tolower(trimws(as.character(recycle_v(scale_source, ""))))
  ss_raw[is.na(ss_raw)] <- ""
  if (is.null(is_calibrated_mm) || length(is_calibrated_mm) == 0L) {
    cal_in <- rep(NA, n)
  } else {
    cal_raw <- recycle_v(is_calibrated_mm, NA)
    if (is.logical(cal_raw)) {
      cal_in <- cal_raw
    } else {
      cal_chr <- tolower(trimws(as.character(cal_raw)))
      cal_in <- rep(NA, n)
      cal_in[cal_chr %in% c("1", "1.0", "true", "t", "yes")] <- TRUE
      cal_in[cal_chr %in% c("0", "0.0", "false", "f", "no")] <- FALSE
    }
  }
  uncal_src <- ss_raw %in% STEP03_UNCALIBRATED_SOURCES
  is_cal <- !is.na(ppm) & is.finite(ppm) & (ppm > 0) & !uncal_src & ifelse(is.na(cal_in), !uncal_src, cal_in)
  ss_out <- ifelse(
    is_cal,
    ifelse(nzchar(ss_raw), ss_raw, "ruler_calibrated"),
    ifelse(
      ss_raw == "tick_lattice_unboxed", "tick_lattice_unboxed",
      ifelse(
        ss_raw %in% c("mp_estimated", "predicted_from_megapixels"), "mp_estimated",
        ifelse(!is.na(ppm) & is.finite(ppm) & (ppm > 0) & !is.na(cal_in) & !cal_in & !nzchar(ss_raw), "mp_estimated", "uncalibrated")
      )
    )
  )
  list(is_calibrated_mm = as.logical(is_cal), scale_source = as.character(ss_out))
}

discover_step03_photometry_cols <- function(col_names) {
  cn <- as.character(col_names)
  unique(intersect(
    c("blade_cie_l_norm", "blade_cie_l_sd", "blade_cie_l_sheet_delta",
      "blade_cie_a_delta", "blade_cie_b_delta", grep("^blade_(cie_|.*browning)", cn, value = TRUE)),
    cn
  ))
}

strip_deprecated_step03_cli_args <- function(cli_args) {
  args <- as.character(cli_args)
  if (length(args) == 0L) return(character(0))
  out <- character(0)
  i <- 1L
  val_flags <- c("--permutations", "-r", "--report-out", "-g", "--plot-out")
  while (i <= length(args)) {
    tok <- args[i]
    if (grepl("^(--permutations|--report-out|--plot-out|-r|-g)=", tok)) {
      i <- i + 1L
      next
    }
    if (tok %in% val_flags) {
      i <- if (i < length(args) && !is.na(args[i + 1L]) && !startsWith(args[i + 1L], "-")) i + 2L else i + 1L
      next
    }
    out <- c(out, tok)
    i <- i + 1L
  }
  out
}
