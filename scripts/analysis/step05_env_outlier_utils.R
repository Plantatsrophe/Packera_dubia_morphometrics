# Step 05 Coordinate Quarantine & Azonal-Aware Environmental MCD Screening
.NA_HERBARIUM_ROOFTOPS <- data.frame(
  code = c("NCU", "DUKE", "NCSC", "NY", "GH", "US", "MO", "F", "PH", "BRIT", "TEX", "WIS",
           "MIN", "TENN", "GA", "WILLI", "LSU", "CM", "MU", "CSCN", "DAO", "CAN", "MT", "TRT"),
  lat = c(35.9108, 36.0014, 35.7877, 40.8623, 42.3785, 38.8889, 38.6128, 41.8663,
          39.9571, 32.7397, 30.2850, 43.0742, 44.9845, 35.9544, 33.9480, 37.2707,
          30.4133, 40.4433, 39.5097, 42.8186, 45.3920, 45.3920, 45.5569, 43.6677),
  lon = c(-79.0505, -78.9382, -78.6703, -73.8803, -71.1160, -77.0260, -90.2594, -87.6170,
          -75.1714, -97.3631, -97.7340, -89.4012, -93.1832, -83.9295, -83.3773, -76.7153,
          -91.1800, -79.9496, -84.7345, -103.0005, -75.7160, -75.7160, -73.5575, -79.3948),
  stringsAsFactors = FALSE
)

.coerce_bool_flag_vec <- function(x, n) {
  if (n <= 0L) return(logical(0))
  if (is.null(x) || length(x) == 0L) return(rep(FALSE, n))
  if (length(x) != n) x <- rep_len(x, n)
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & is.finite(x) & (x != 0))
  sx <- tolower(trimws(as.character(x)))
  !is.na(sx) & (sx %in% c("true", "1", "t", "yes"))
}

.resolve_type_flag_vec <- function(df, is_type_vec = NULL) {
  n <- if (is.data.frame(df)) nrow(df) else 0L
  if (n == 0L) return(logical(0))
  if (!is.null(is_type_vec) && length(is_type_vec) == n) return(.coerce_bool_flag_vec(is_type_vec, n))
  if (exists(".identify_type_specimen_vec", mode = "function")) {
    t_res <- tryCatch(.identify_type_specimen_vec(df), error = function(e) NULL)
    if (!is.null(t_res$is_type_specimen) && length(t_res$is_type_specimen) == n) {
      return(.coerce_bool_flag_vec(t_res$is_type_specimen, n))
    }
  }
  if ("is_type_specimen" %in% names(df)) return(.coerce_bool_flag_vec(df$is_type_specimen, n))
  rep(FALSE, n)
}

.coalesce_coord_col <- function(df, cols, n) {
  out <- rep(NA_real_, n)
  for (cn in intersect(cols, names(df))) {
    v <- suppressWarnings(as.numeric(df[[cn]]))
    miss <- !is.finite(out) & is.finite(v)
    if (any(miss)) out[miss] <- v[miss]
  }
  out
}

detect_r_centroid_or_inst_quarantine <- function(df, is_type_vec = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    return(list(is_quarantined = logical(0), reason = character(0)))
  }
  n <- nrow(df)
  is_t <- .resolve_type_flag_vec(df, is_type_vec = is_type_vec)
  is_q <- rep(FALSE, n)
  reason <- rep("", n)

  if ("is_centroid_or_inst_coord" %in% names(df)) {
    is_q <- .coerce_bool_flag_vec(df$is_centroid_or_inst_coord, n)
    if ("coord_quarantine_reason" %in% names(df)) {
      r_raw <- trimws(as.character(df$coord_quarantine_reason)); r_raw[is.na(r_raw)] <- ""
      reason <- r_raw
    }
    reason[!is_q] <- ""
    need_def <- is_q & !nzchar(reason)
    if (any(need_def)) reason[need_def] <- "centroid_or_inst_flagged"
  } else if ("locality" %in% names(df) &&
             any(c("decimalLatitude", "latitude") %in% names(df)) &&
             any(c("decimalLongitude", "longitude") %in% names(df))) {
    lat_v <- .coalesce_coord_col(df, c("decimalLatitude", "latitude"), n)
    lon_v <- .coalesce_coord_col(df, c("decimalLongitude", "longitude"), n)
    loc_v <- as.character(df$locality); loc_v[is.na(loc_v)] <- ""
    hab_v <- if ("habitat" %in% names(df)) as.character(df$habitat) else rep("", n); hab_v[is.na(hab_v)] <- ""
    loc_hab <- trimws(paste(loc_v, hab_v))
    wild_rx <- paste0(
      "(^|[^A-Za-z])([0-9]+(\\.[0-9]+)?[[:space:]]*(km|mi|miles?|meters?|m|ft|yds?)([^A-Za-z]|$)|",
      "(north|south|east|west|ne|nw|se|sw|nne|nnw|ene|ese|sse|ssw|wnw|wsw)[[:space:]]+of|",
      "along|creek|stream|river|branch|brook|seep|fen|bog|swamp|marsh|floodplain|bottomland|wetland|",
      "outcrop|flatrock|glade|bluff|cliff|talus|ridge|slope|ravine|hollow|gorge|cove|forest|woods|",
      "woodland|meadow|prairie|savanna|barren|trail|roadside|hwy|highway|preserve|natural[[:space:]]+area|",
      "granite|gneiss|limestone|shale|serpentine|diabase|sandstone)([^A-Za-z]|$)"
    )
    cand_idx <- which(
      is.finite(lat_v) & is.finite(lon_v) & abs(lat_v) <= 90.0 & abs(lon_v) <= 180.0 &
      !is_t & !grepl(wild_rx, loc_hab, ignore.case = TRUE)
    )
    if (length(cand_idx) > 0L) {
      r_rad <- pi / 180.0
      lat_c <- lat_v[cand_idx] * r_rad; lon_c <- lon_v[cand_idx] * r_rad
      h_lat <- .NA_HERBARIUM_ROOFTOPS$lat * r_rad; h_lon <- .NA_HERBARIUM_ROOFTOPS$lon * r_rad
      a_m <- sin(outer(lat_c, h_lat, "-") * 0.5)^2 +
        outer(cos(lat_c), cos(h_lat), "*") * sin(outer(lon_c, h_lon, "-") * 0.5)^2
      d_m <- 2.0 * 6371000.0 * atan2(sqrt(pmax(0.0, a_m)), sqrt(pmax(0.0, 1.0 - a_m)))
      min_j <- apply(d_m, 1L, which.min)
      min_d <- d_m[cbind(seq_along(cand_idx), min_j)]
      hit <- which(is.finite(min_d) & (min_d <= 250.0))
      if (length(hit) > 0L) {
        orig_i <- cand_idx[hit]
        is_q[orig_i] <- TRUE
        reason[orig_i] <- sprintf("herbarium_rooftop_%s_%.0fm", .NA_HERBARIUM_ROOFTOPS$code[min_j[hit]], min_d[hit])
      }
    }
  }

  is_q[is_t] <- FALSE
  reason[is_t] <- ""
  list(is_quarantined = is_q, reason = reason)
}

has_corroborated_azonal_habitat <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(logical(0))
  n <- nrow(df)
  az <- rep(FALSE, n)
  for (bc in intersect(c("is_azonal_lithology", "azonal_lithology_presence_buffer", "is_rock_outcrop", "is_nwi_wetland"), names(df))) {
    az <- az | .coerce_bool_flag_vec(df[[bc]], n)
  }
  if ("soil_ultramafic_ni_cr_index" %in% names(df)) {
    ni_cr <- suppressWarnings(as.numeric(df$soil_ultramafic_ni_cr_index))
    az <- az | (is.finite(ni_cr) & (ni_cr >= 3.85))
  }
  if ("sgmc_lithology_class" %in% names(df)) {
    sg_v <- as.character(df$sgmc_lithology_class); sg_v[is.na(sg_v)] <- ""
    az <- az | grepl("(Ultramafic|Serpentine|Mafic|Diabase|Gabbro|Carbonate|Calcareous)", sg_v, ignore.case = TRUE)
  }
  if ("outcrop_class" %in% names(df)) {
    oc_v <- trimws(as.character(df$outcrop_class)); oc_v[is.na(oc_v)] <- ""
    az <- az | (nzchar(oc_v) & (tolower(oc_v) != "non-outcrop matrix") & (tolower(oc_v) != "na"))
  }
  az_rx <- paste0(
    "(^|[^A-Za-z])(outcrop|flatrock|granite|gneiss|serpentine|ultramafic|diabase|gabbro|",
    "limestone|dolomite|marble|cedar[[:space:]]+glade|shale[[:space:]]+barren|",
    "bluff|cliff|talus|seep|fen|bog|pocosin|savanna|bald)([^A-Za-z]|$)"
  )
  for (tc in intersect(c("locality", "habitat", "muname"), names(df))) {
    tv <- as.character(df[[tc]]); tv[is.na(tv)] <- ""
    az <- az | grepl(az_rx, tv, ignore.case = TRUE)
  }
  az
}

screen_within_taxon_env_mcd_outliers <- function(
  df, X_env_mat, group_vec, active_anchor_mask,
  is_type_vec = NULL, p_cutoff = 0.005, min_n_group = 12L, seed = 42L
) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(as.integer(seed))

  n <- if (is.data.frame(df)) nrow(df) else if (!is.null(X_env_mat)) nrow(as.matrix(X_env_mat)) else 0L
  env_mahal_d2 <- rep(NA_real_, n); env_mcd_outlier_flag <- rep(FALSE, n)
  env_outlier_rescued_azonal <- rep(FALSE, n); raw_mcd_exceed_flag <- rep(FALSE, n)
  X_mat <- if (!is.null(X_env_mat)) as.matrix(X_env_mat) else matrix(numeric(0), nrow = 0L, ncol = 0L)
  storage.mode(X_mat) <- "double"
  is_t <- .resolve_type_flag_vec(df, is_type_vec = is_type_vec)
  azonal_vec <- has_corroborated_azonal_habitat(df)
  if (n == 0L || nrow(X_mat) != n || ncol(X_mat) == 0L) {
    return(list(env_mahal_d2 = env_mahal_d2, env_mcd_outlier_flag = env_mcd_outlier_flag,
                env_outlier_rescued_azonal = env_outlier_rescued_azonal, raw_mcd_exceed_flag = raw_mcd_exceed_flag,
                is_type_specimen = is_t, is_azonal_corroborated = azonal_vec))
  }

  act_m <- if (is.null(active_anchor_mask)) rep(TRUE, n) else .coerce_bool_flag_vec(active_anchor_mask, n)
  grp_chr <- as.character(group_vec)
  act_rows <- which(act_m & !is.na(grp_chr) & nzchar(trimws(grp_chr)))

  for (g in unique(grp_chr[act_rows])) {
    idx_g <- act_rows[grp_chr[act_rows] == g]; n_g <- length(idx_g)
    if (n_g < as.integer(min_n_group)) next
    Xg <- X_mat[idx_g, , drop = FALSE]
    keep_j <- vapply(seq_len(ncol(Xg)), function(j) {
      vj <- Xg[, j]; fin <- is.finite(vj)
      sum(fin) >= 2L && length(unique(vj[fin])) >= 2L && is.finite(stats::sd(vj[fin])) && stats::sd(vj[fin]) > 1e-12
    }, logical(1))
    if (!any(keep_j)) next
    Xg <- Xg[, keep_j, drop = FALSE]
    Zg <- matrix(0.0, nrow = n_g, ncol = ncol(Xg))
    for (j in seq_len(ncol(Xg))) {
      vj <- Xg[, j]; fin <- is.finite(vj); med_j <- stats::median(vj[fin], na.rm = TRUE); vj[!fin] <- med_j
      Zg[, j] <- (vj - med_j) / max(1.4826 * stats::mad(vj, constant = 1, na.rm = TRUE), 0.15 * stats::sd(vj), 1e-6)
    }

    pr <- tryCatch(stats::prcomp(Zg, center = TRUE, scale. = FALSE), error = function(e) NULL)
    if (is.null(pr) || is.null(pr$x) || ncol(pr$x) == 0L) next
    eig_v <- pr$sdev^2; tot_v <- sum(eig_v[is.finite(eig_v)])
    n_avail <- sum(is.finite(pr$sdev) & (pr$sdev > 1e-6))
    if (!is.finite(tot_v) || tot_v <= 1e-12 || n_avail == 0L) next
    k_cap <- max(2L, min(4L, as.integer(floor((n_g - 4L) / 2L))))
    k_85 <- which((cumsum(eig_v) / tot_v) >= 0.85)[1L]; if (is.na(k_85)) k_85 <- length(eig_v)
    k_pc <- min(k_cap, n_avail, k_85)
    if (ncol(pr$x) >= 2L && n_avail >= 2L && k_pc < 2L) k_pc <- 2L
    Z_pc <- pr$x[, seq_len(k_pc), drop = FALSE]

    mcd_fit <- if (requireNamespace("MASS", quietly = TRUE)) {
      tryCatch({
        fit <- MASS::cov.mcd(Z_pc, quantile.used = as.integer(floor((n_g + k_pc + 1L) / 2L)))
        cov_m <- matrix(as.numeric(fit$cov), nrow = k_pc, ncol = k_pc); ctr_v <- as.numeric(fit$center)
        if (length(ctr_v) != k_pc || any(!is.finite(ctr_v)) || any(!is.finite(cov_m))) stop("non-finite")
        if (is.null(tryCatch(solve(cov_m), error = function(e) NULL))) stop("singular")
        list(center = ctr_v, cov = cov_m)
      }, error = function(e) NULL)
    } else NULL
    if (is.null(mcd_fit)) {
      ctr_fb <- apply(Z_pc, 2L, stats::median, na.rm = TRUE)
      mad_fb <- pmax(1.4826 * apply(Z_pc, 2L, stats::mad, constant = 1, na.rm = TRUE), 1e-3)
      mcd_fit <- list(center = ctr_fb, cov = diag(mad_fb^2, nrow = k_pc, ncol = k_pc))
    }

    d2_env <- tryCatch(stats::mahalanobis(Z_pc, mcd_fit$center, mcd_fit$cov + diag(1e-6, nrow = k_pc, ncol = k_pc)),
                       error = function(e) rep(NA_real_, n_g))
    cand <- is.finite(d2_env) & (d2_env > stats::qchisq(1.0 - p_cutoff, df = k_pc))
    env_mahal_d2[idx_g] <- round(d2_env, 4); raw_mcd_exceed_flag[idx_g] <- cand
    env_outlier_rescued_azonal[idx_g] <- cand & (!is_t[idx_g]) & azonal_vec[idx_g]
    env_mcd_outlier_flag[idx_g] <- cand & (!is_t[idx_g]) & (!azonal_vec[idx_g])
  }
  list(env_mahal_d2 = env_mahal_d2, env_mcd_outlier_flag = env_mcd_outlier_flag,
       env_outlier_rescued_azonal = env_outlier_rescued_azonal, raw_mcd_exceed_flag = raw_mcd_exceed_flag,
       is_type_specimen = is_t, is_azonal_corroborated = azonal_vec)
}

write_env_outlier_quarantine_audit <- function(df, screen_res, out_csv, final_anchor_mask = NULL) {
  if (is.null(df) || !is.data.frame(df)) return(invisible(NULL))
  n <- nrow(df)
  pick_scr <- function(nm, def_v) {
    if (is.list(screen_res) && (nm %in% names(screen_res)) && length(screen_res[[nm]]) == n) screen_res[[nm]]
    else if (nm %in% names(df)) df[[nm]] else def_v
  }
  is_t <- .coerce_bool_flag_vec(pick_scr("is_type_specimen", .resolve_type_flag_vec(df)), n)
  q_res <- detect_r_centroid_or_inst_quarantine(df, is_type_vec = is_t)
  d2_v <- suppressWarnings(as.numeric(pick_scr("env_mahal_d2", rep(NA_real_, n))))
  mcd_f <- .coerce_bool_flag_vec(pick_scr("env_mcd_outlier_flag", rep(FALSE, n)), n)
  res_f <- .coerce_bool_flag_vec(pick_scr("env_outlier_rescued_azonal", rep(FALSE, n)), n)
  final_anc <- if (!is.null(final_anchor_mask) && length(final_anchor_mask) == n) .coerce_bool_flag_vec(final_anchor_mask, n) else !(q_res$is_quarantined | mcd_f) | is_t

  status <- rep("Passive_Excluded", n)
  status[final_anc] <- "Active_Anchor"
  status[(!final_anc) & mcd_f] <- "Passive_Quarantined_Env_MCD"
  status[(!final_anc) & q_res$is_quarantined] <- "Passive_Quarantined_Coord"
  status[res_f & (final_anc | !q_res$is_quarantined)] <- "Azonal_Rescued"
  status[is_t] <- "Type_Immune"

  sp_v <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else rep(NA_character_, n)
  audit_df <- data.frame(
    catalogNumber = if ("catalogNumber" %in% names(df)) as.character(df$catalogNumber) else as.character(seq_len(n)),
    species_standardized = sp_v,
    analysis_group = if ("analysis_group" %in% names(df)) as.character(df$analysis_group) else sp_v,
    is_type_specimen = is_t, is_centroid_or_inst_coord = q_res$is_quarantined, coord_quarantine_reason = q_res$reason,
    env_mahal_d2 = d2_v, env_mcd_outlier_flag = mcd_f, env_outlier_rescued_azonal = res_f,
    final_anchor_status = status, stringsAsFactors = FALSE
  )
  if (!is.null(out_csv) && nzchar(trimws(as.character(out_csv)))) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    if (exists("write_csv_atomic", mode = "function")) write_csv_atomic(audit_df, out_csv, demo = FALSE)
    else { tmp_csv <- paste0(out_csv, ".tmp"); utils::write.csv(audit_df, tmp_csv, row.names = FALSE, na = ""); file.rename(tmp_csv, out_csv) }
  }
  invisible(audit_df)
}
