# ==============================================================================
# Step 05 Module: Config, Haversine Distance, SQLite Phenology, & Circular Tests
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

load_step05_macroecology_config <- function(cfg_path = "config/config.yaml") {
  defaults <- list(
    cor_threshold = 0.70, vif_threshold = 5.0, cumulative_importance_cutoff = 0.85,
    min_node_size = 12L, max_depth = 6L, sample_fraction = 0.70,
    ecotype_min_silhouette = 0.40, ecotype_min_fold_jaccard = 0.60,
    min_flowering_doy = 30.0, max_flowering_doy = 220.0, min_flowering_lat = 24.0, max_flowering_lat = 55.0
  )
  if (is.null(cfg_path) || length(cfg_path) == 0L || is.na(cfg_path[1L])) return(defaults)
  cpath <- trimws(as.character(cfg_path[1L]))
  if (!nzchar(cpath) || !file.exists(cpath)) return(defaults)

  is_valid_cfg <- function(k, v) {
    if (length(v) != 1L || !is.finite(v)) return(FALSE)
    if (k %in% c("cor_threshold", "cumulative_importance_cutoff", "sample_fraction", "ecotype_min_silhouette", "ecotype_min_fold_jaccard")) return(v > 0.0 && v <= 1.0)
    if (identical(k, "vif_threshold")) return(v > 1.0)
    if (k %in% c("min_node_size", "max_depth")) return(v >= 1.0)
    if (k %in% c("min_flowering_doy", "max_flowering_doy")) return(v >= 1.0 && v <= 366.0)
    if (k %in% c("min_flowering_lat", "max_flowering_lat")) return(v >= -90.0 && v <= 90.0)
    TRUE
  }

  out <- defaults; found <- character(0)
  if (requireNamespace("yaml", quietly = TRUE)) {
    cfg <- tryCatch(yaml::read_yaml(cpath), error = function(e) NULL)
    if (is.list(cfg)) {
      sec_map <- c(
        cor_threshold = "spatial_rf", vif_threshold = "spatial_rf", cumulative_importance_cutoff = "spatial_rf",
        min_node_size = "spatial_rf", max_depth = "spatial_rf", sample_fraction = "spatial_rf",
        ecotype_min_silhouette = "spatial_rf", ecotype_min_fold_jaccard = "spatial_rf",
        min_flowering_doy = "phenology", max_flowering_doy = "phenology", min_flowering_lat = "phenology", max_flowering_lat = "phenology"
      )
      for (k in names(sec_map)) {
        sec <- sec_map[[k]]
        cands <- list(
          if (is.list(cfg$macroecology) && is.list(cfg$macroecology[[sec]])) cfg$macroecology[[sec]][[k]] else NULL,
          if (is.list(cfg$macroecology)) cfg$macroecology[[k]] else NULL,
          if (is.list(cfg[[sec]])) cfg[[sec]][[k]] else NULL
        )
        for (cv in cands) {
          nv <- suppressWarnings(as.numeric(cv))
          if (length(nv) == 1L && is_valid_cfg(k, nv)) { out[[k]] <- nv; found <- c(found, k); break }
        }
      }
    }
  }
  miss_keys <- setdiff(names(defaults), found)
  if (length(miss_keys) > 0L) {
    lines <- tryCatch(readLines(cpath, warn = FALSE), error = function(e) character(0))
    for (k in miss_keys) {
      for (hl in grep(paste0("^[[:space:]]*", k, "[[:space:]]*:"), lines, value = TRUE)) {
        nv <- suppressWarnings(as.numeric(gsub("^['\"]|['\"]$", "", trimws(sub("#.*$", "", sub("^[^:]+:[[:space:]]*", "", hl))))))
        if (length(nv) == 1L && is_valid_cfg(k, nv)) { out[[k]] <- nv; break }
      }
    }
  }
  out$min_node_size <- as.integer(round(out$min_node_size))
  out$max_depth <- as.integer(round(out$max_depth))
  out
}

compute_haversine_dist_km_matrix <- function(coords_a, coords_b = NULL) {
  is_self <- is.null(coords_b)
  to_coord_mat <- function(x) {
    if (is.null(x)) return(matrix(numeric(0), nrow = 0L, ncol = 2L))
    if (!is.matrix(x) && !is.data.frame(x) && is.numeric(x) && length(x) == 2L) return(matrix(as.double(x), nrow = 1L, ncol = 2L))
    m <- as.matrix(x); storage.mode(m) <- "double"; m
  }
  ca <- to_coord_mat(coords_a); cb <- if (is_self) ca else to_coord_mat(coords_b)
  if (nrow(ca) == 0L || nrow(cb) == 0L) return(matrix(numeric(0), nrow = nrow(ca), ncol = nrow(cb)))
  if (ncol(ca) < 2L || ncol(cb) < 2L) stop("Coordinate matrices must contain longitude and latitude columns.")
  deg2rad <- pi / 180.0
  lon_a <- ca[, 1L] * deg2rad; lat_a <- ca[, 2L] * deg2rad
  lon_b <- cb[, 1L] * deg2rad; lat_b <- cb[, 2L] * deg2rad
  a <- sin(outer(lat_a, lat_b, "-") / 2.0)^2 + outer(cos(lat_a), cos(lat_b), "*") * sin(outer(lon_a, lon_b, "-") / 2.0)^2
  a <- pmin(pmax(a, 0.0), 1.0)
  d_mat <- 2.0 * 6371.0 * asin(sqrt(a))
  dim(d_mat) <- c(nrow(ca), nrow(cb))
  if (is_self) { d_mat <- 0.5 * (d_mat + t(d_mat)); diag(d_mat) <- 0.0 }
  d_mat
}

is_post_selection_ecotype_groups <- function(group_vec, intraspecific_mode = FALSE, df = NULL) {
  if (!is.null(df) && isTRUE(attr(df, "is_post_selection_ecotype"))) return(TRUE)
  if (!isTRUE(intraspecific_mode) || is.null(group_vec)) return(FALSE)
  vals <- unique(trimws(as.character(group_vec)))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  length(vals) >= 2L && all(grepl("^Ecotype_[0-9]+_", vals))
}

compute_fold_reclustered_ecotype_ref <- function(X_tr_imp, y_tr, X_te_imp, seed = 42L) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({ if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv) }, add = TRUE)
  set.seed(as.integer(seed))

  if (is.null(X_te_imp)) return(character(0))
  X_te_df <- as.data.frame(X_te_imp); n_te <- nrow(X_te_df)
  if (n_te == 0L) return(character(0))
  y_fac <- droplevels(factor(y_tr)); lvls <- levels(y_fac); k_lvls <- length(lvls)
  if (k_lvls < 2L) return(rep(if (k_lvls == 1L) lvls[1L] else as.character(y_tr[1L]), n_te))

  X_tr_m <- as.matrix(X_tr_imp); storage.mode(X_tr_m) <- "double"
  X_te_m <- if (!is.null(colnames(X_tr_m)) && all(colnames(X_tr_m) %in% names(X_te_df))) as.matrix(X_te_df[, colnames(X_tr_m), drop = FALSE]) else as.matrix(X_te_imp)
  storage.mode(X_te_m) <- "double"
  if (ncol(X_tr_m) == 0L || ncol(X_te_m) != ncol(X_tr_m)) return(rep(lvls[1L], n_te))

  col_mu <- apply(X_tr_m, 2L, function(v) { fv <- v[is.finite(v)]; if (length(fv) > 0L) mean(fv) else 0.0 })
  col_sd <- apply(X_tr_m, 2L, function(v) { fv <- v[is.finite(v)]; s <- if (length(fv) >= 2L) stats::sd(fv) else 1.0; if (is.finite(s)) max(s, 1e-6) else 1.0 })
  for (j in seq_len(ncol(X_tr_m))) {
    b_tr <- !is.finite(X_tr_m[, j]); if (any(b_tr)) X_tr_m[b_tr, j] <- col_mu[j]
    b_te <- !is.finite(X_te_m[, j]); if (any(b_te)) X_te_m[b_te, j] <- col_mu[j]
  }
  Z_tr <- sweep(sweep(X_tr_m, 2L, col_mu, "-"), 2L, col_sd, "/"); Z_tr[!is.finite(Z_tr)] <- 0.0
  Z_te <- sweep(sweep(X_te_m, 2L, col_mu, "-"), 2L, col_sd, "/"); Z_te[!is.finite(Z_te)] <- 0.0
  if (nrow(Z_tr) <= k_lvls || ncol(Z_tr) == 0L || nrow(unique(round(Z_tr, 8L))) < k_lvls) return(rep(lvls[1L], n_te))

  nearest_idx <- function(Z, C) {
    d_m <- matrix(vapply(seq_len(nrow(C)), function(u) rowSums(sweep(Z, 2L, C[u, ], "-")^2), numeric(nrow(Z))), nrow = nrow(Z), ncol = nrow(C))
    max.col(-d_m, ties.method = "first")
  }
  sup_c <- matrix(0.0, nrow = k_lvls, ncol = ncol(Z_tr))
  for (s_i in seq_len(k_lvls)) {
    idx_s <- which(!is.na(y_fac) & as.integer(y_fac) == s_i)
    sup_c[s_i, ] <- if (length(idx_s) > 0L) colMeans(Z_tr[idx_s, , drop = FALSE]) else colMeans(Z_tr)
  }
  km <- tryCatch(stats::kmeans(Z_tr, centers = k_lvls, nstart = 20L, iter.max = 50L), error = function(e) NULL)
  if (is.null(km) || is.null(km$centers) || nrow(km$centers) < k_lvls || any(!is.finite(km$centers))) return(as.character(lvls[nearest_idx(Z_te, sup_c)]))

  C_mat <- matrix(0.0, nrow = k_lvls, ncol = k_lvls)
  for (u in seq_len(k_lvls)) for (s_i in seq_len(k_lvls)) C_mat[u, s_i] <- sum((km$centers[u, ] - sup_c[s_i, ])^2)
  if (k_lvls <= 4L) {
    gen_perms <- function(vals) if (length(vals) <= 1L) list(vals) else unlist(lapply(seq_along(vals), function(i) lapply(gen_perms(vals[-i]), function(sp) c(vals[i], sp))), recursive = FALSE)
    perms <- gen_perms(seq_len(k_lvls))
    costs <- vapply(perms, function(p) sum(vapply(seq_len(k_lvls), function(u) C_mat[u, p[u]], numeric(1))), numeric(1))
    perm <- perms[[which.min(costs)[1L]]]
  } else {
    perm <- integer(k_lvls); C_work <- C_mat
    for (step in seq_len(k_lvls)) {
      pos <- which.min(C_work); u_b <- ((pos - 1L) %% k_lvls) + 1L; s_b <- ((pos - 1L) %/% k_lvls) + 1L
      perm[u_b] <- s_b; C_work[u_b, ] <- Inf; C_work[, s_b] <- Inf
    }
  }
  as.character(lvls[perm[nearest_idx(Z_te, km$centers)]])
}

fetch_step02_phenology_from_sqlite <- function(catalog_numbers = NULL, vouchers_path = NULL, demo = FALSE) {
  if (isTRUE(demo)) return(NULL)
  sql_ok <- exists("connect_packera_db", mode = "function") && exists("db_table_exists_safe", mode = "function") &&
    ((exists("is_sqlite_active_for_inputs", mode = "function") && is_sqlite_active_for_inputs(vouchers_path, demo = FALSE)) ||
     nzchar(Sys.getenv("PLANTMATICS_DB_PATH", "")) || nzchar(Sys.getenv("PACKERA_DB_PATH", "")) || nzchar(Sys.getenv("RUN_ID", "")))

  df_p <- NULL
  if (isTRUE(sql_ok)) {
    conn <- tryCatch(connect_packera_db(read_only = TRUE, demo = FALSE), error = function(e) NULL)
    if (!is.null(conn)) {
      on.exit(if (exists("disconnect_packera_db", mode = "function")) tryCatch(disconnect_packera_db(conn), error = function(e) NULL), add = TRUE)
      if (db_table_exists_safe(conn, "voucher_reproductive_phenology")) {
        fields <- tryCatch(DBI::dbListFields(conn, "voucher_reproductive_phenology"), error = function(e) character(0))
        if (length(fields) > 0L) df_p <- tryCatch(DBI::dbGetQuery(conn, "SELECT * FROM voucher_reproductive_phenology;"), error = function(e) NULL)
      }
    }
  }
  if ((is.null(df_p) || !is.data.frame(df_p) || nrow(df_p) == 0L) && !isTRUE(demo) && file.exists("data/tables/voucher_phenological_states.csv")) {
    df_p <- tryCatch(utils::read.csv("data/tables/voucher_phenological_states.csv", stringsAsFactors = FALSE), error = function(e) NULL)
  }
  if (is.null(df_p) || !is.data.frame(df_p) || nrow(df_p) == 0L || !("catalogNumber" %in% names(df_p))) return(NULL)

  df_p$catalogNumber <- as.character(df_p$catalogNumber)
  df_p <- df_p[!is.na(df_p$catalogNumber) & nzchar(trimws(df_p$catalogNumber)), , drop = FALSE]
  if (!is.null(catalog_numbers) && length(catalog_numbers) > 0L) df_p <- df_p[df_p$catalogNumber %in% as.character(catalog_numbers), , drop = FALSE]
  if (nrow(df_p) == 0L) return(NULL)

  cats_u <- unique(df_p$catalogNumber)
  row_groups <- split(seq_len(nrow(df_p)), factor(df_p$catalogNumber, levels = cats_u))
  sum_col <- function(col_nm) {
    if (!(col_nm %in% names(df_p))) return(rep(NA_real_, length(cats_u)))
    v <- suppressWarnings(as.numeric(df_p[[col_nm]]))
    vapply(row_groups, function(idx) { f <- v[idx]; f <- f[is.finite(f)]; if (length(f) > 0L) sum(f) else NA_real_ }, numeric(1))
  }
  max_col_num <- function(col_nm, is_bool = FALSE) {
    if (!(col_nm %in% names(df_p))) return(rep(NA_real_, length(cats_u)))
    raw <- df_p[[col_nm]]
    v <- if (is_bool) {
      if (is.logical(raw)) as.numeric(raw)
      else if (is.numeric(raw) || is.integer(raw)) ifelse(is.na(raw), NA_real_, as.numeric(raw != 0))
      else { u <- toupper(trimws(as.character(raw))); ifelse(u %in% c("TRUE", "T", "1", "YES"), 1.0, ifelse(u %in% c("FALSE", "F", "0", "NO"), 0.0, NA_real_)) }
    } else suppressWarnings(as.numeric(raw))
    vapply(row_groups, function(idx) { f <- v[idx]; f <- f[is.finite(f)]; if (length(f) > 0L) max(f) else NA_real_ }, numeric(1))
  }
  coalesce_num <- function(...) {
    out <- rep(NA_real_, length(cats_u))
    for (vc in list(...)) { miss <- !is.finite(out) & is.finite(vc); if (any(miss)) out[miss] <- vc[miss] }
    out
  }
  pick_state <- function(idx) {
    if (!("voucher_phenological_state" %in% names(df_p))) return(NA_character_)
    s <- trimws(as.character(df_p$voucher_phenological_state[idx]))
    s <- s[!is.na(s) & nzchar(s) & !(tolower(s) %in% c("na", "nan", "null", "none"))]
    if (length(s) == 0L) return(NA_character_)
    rk <- ifelse(grepl("anthesis", s, ignore.case = TRUE), 1L, ifelse(grepl("fruit", s, ignore.case = TRUE), 2L, ifelse(grepl("bud", s, ignore.case = TRUE), 3L, ifelse(grepl("sterile", s, ignore.case = TRUE), 4L, 5L))))
    s[which.min(rk)[1L]]
  }

  r_tot <- sum_col("total_capitula"); r_cap <- sum_col("capitula_count"); r_cor <- sum_col("corymb_head_count_est")
  r_ant <- sum_col("n_anthesis"); r_ahc <- sum_col("anthesis_head_count"); r_fru <- sum_col("n_fruiting"); r_fhc <- sum_col("fruiting_head_count")
  out <- data.frame(
    catalogNumber = cats_u, total_capitula = coalesce_num(r_tot, r_cap, r_cor), corymb_head_count_est = coalesce_num(r_cor, r_tot, r_cap),
    n_anthesis = coalesce_num(r_ant, r_ahc), anthesis_head_count = coalesce_num(r_ahc, r_ant),
    n_fruiting = coalesce_num(r_fru, r_fhc), n_bud = sum_col("n_bud"), n_senescent = sum_col("n_senescent"),
    voucher_phenological_state = vapply(row_groups, pick_state, character(1)),
    reproductive_stage_confidence = max_col_num("reproductive_stage_confidence"),
    dominant_pappus_ratio = max_col_num("dominant_pappus_ratio"), stringsAsFactors = FALSE
  )
  if ("capitula_count" %in% names(df_p)) out$capitula_count <- r_cap
  if ("fruiting_head_count" %in% names(df_p)) out$fruiting_head_count <- r_fhc
  if ("is_verified_anthesis" %in% names(df_p)) out$is_verified_anthesis <- as.integer(max_col_num("is_verified_anthesis", is_bool = TRUE))
  rownames(out) <- NULL; out
}

summarize_weighted_circular_phenology <- function(target_flowering, active_tx, circ_phase_test, is_unif = FALSE, B = 999L, seed = 42L) {
  day_per_rad <- 365.25 / (2.0 * pi)
  to_rad <- if (exists("doy_to_rad", mode = "function")) doy_to_rad else function(d) 2.0 * pi * (as.numeric(d) / 365.25)
  mk_row <- function(sec, comp, n, md = NA_real_, sdd = NA_real_, med = NA_real_, iqr = NA_real_, lo = NA_real_, hi = NA_real_, s1 = NA_real_, s2 = NA_real_, eff = NA_real_, pv = "NA", interp = "") {
    data.frame(Analysis_Section = sec, Taxon_or_Comparison = comp, Sample_Size = as.integer(n), Mean_Delta_DOY = md, SD_Delta_DOY = sdd,
               Median_Delta_DOY = med, IQR_Delta_DOY = iqr, CI_Lower_95 = lo, CI_Upper_95 = hi, Statistic_Primary = s1,
               Statistic_Secondary = s2, Effect_or_R2 = eff, P_Value = pv, Interpretation = interp, stringsAsFactors = FALSE)
  }

  n_tot <- if (!is.null(target_flowering) && is.data.frame(target_flowering)) nrow(target_flowering) else 0L
  wf <- if (!is.null(circ_phase_test$watson_f) && is.finite(circ_phase_test$watson_f)) as.numeric(circ_phase_test$watson_f) else NA_real_
  ng <- if (!is.null(circ_phase_test$n_groups) && is.finite(circ_phase_test$n_groups)) as.integer(circ_phase_test$n_groups) else 0L
  pp <- if (!is.null(circ_phase_test$perm_p) && is.finite(circ_phase_test$perm_p)) as.numeric(circ_phase_test$perm_p) else NA_real_
  df_b <- max(0L, ng - 1L); df_w <- max(1L, n_tot - ng)
  eff_r2 <- if (is.finite(wf) && df_b > 0L && (wf * df_b + df_w) > 0) round((wf * df_b) / (wf * df_b + df_w), 4) else NA_real_
  r1_int <- if (isTRUE(is_unif)) "Gated: uniform circular phase null (aseasonal cohort; divergence calls suppressed)"
    else if (is.finite(pp) && pp < 0.05) "Significant circular phenological phase divergence between taxa"
    else "No significant circular phenological divergence or insufficient sample size"
  rows <- list(mk_row("Circular_Phase_Test", "Across Target Taxa (Verified Anthesis)", n_tot,
                      s1 = if (is.finite(wf)) round(wf, 3) else NA_real_,
                      s2 = if (!is.null(circ_phase_test$n_groups) && is.finite(circ_phase_test$n_groups)) df_b else NA_real_,
                      eff = eff_r2, pv = if (is.finite(pp)) sprintf("%.4e", pp) else "NA", interp = r1_int))

  tx_stats <- list()
  if (n_tot > 0L && length(active_tx) > 0L) {
    tx_col <- as.character(target_flowering$species_standardized)
    for (tx in active_tx) {
      sub_df <- target_flowering[!is.na(tx_col) & tx_col == tx, , drop = FALSE]
      dd_v <- suppressWarnings(as.numeric(sub_df$delta_doy)); ok_dd <- is.finite(dd_v)
      sub_df <- sub_df[ok_dd, , drop = FALSE]; sub_vals <- dd_v[ok_dd]; n_tx <- length(sub_vals)
      if (n_tx == 0L) next
      aw <- if ("authority_weight" %in% names(sub_df)) suppressWarnings(as.numeric(sub_df$authority_weight)) else rep(1.0, n_tx)
      sw <- if ("spatial_balance_weight" %in% names(sub_df)) suppressWarnings(as.numeric(sub_df$spatial_balance_weight)) else rep(1.0, n_tx)
      aw[!is.finite(aw)] <- 1.0; sw[!is.finite(sw)] <- 1.0
      w <- pmax(aw, 0.25) * pmax(sw, 1e-6); w <- w / sum(w)
      th <- 2.0 * pi * sub_vals / 365.25; S_w <- sum(w * sin(th)); C_w <- sum(w * cos(th))
      mu_rad <- atan2(S_w, C_w); mu_doy <- mu_rad * day_per_rad
      sd_doy <- sqrt(max(0.0, -2.0 * log(min(1.0, pmax(sqrt(S_w^2 + C_w^2), 1e-8))))) * day_per_rad
      se_doy <- sd_doy / sqrt(n_tx)
      tx_stats[[tx]] <- list(sub_df = sub_df, n = n_tx, mu_rad = mu_rad, sd_doy = sd_doy)
      tx_int <- if (isTRUE(is_unif)) "Gated: aseasonal / uniform phase null (timing anomaly suppressed)"
        else if (mu_doy < -3.0) "Early blooming relative to latitude" else if (mu_doy > 3.0) "Late blooming relative to latitude" else "Conforms to latitudinal cline"
      rows[[length(rows) + 1L]] <- mk_row("Taxon_Distribution", tx, n_tx, md = round(mu_doy, 3), sdd = round(sd_doy, 3),
                                          med = round(stats::median(sub_vals), 3), iqr = round(stats::IQR(sub_vals), 3),
                                          lo = round(mu_doy - 1.96 * se_doy, 3), hi = round(mu_doy + 1.96 * se_doy, 3), pv = "", interp = tx_int)
    }
  }

  pair_recs <- list()
  if (length(active_tx) >= 2L) {
    eval_pair <- function(pr) {
      tx1 <- pr[1L]; tx2 <- pr[2L]; s1 <- tx_stats[[tx1]]; s2 <- tx_stats[[tx2]]
      if (is.null(s1) || is.null(s2) || s1$n < 3L || s2$n < 3L) return(NULL)
      diff_doy <- atan2(sin(s1$mu_rad - s2$mu_rad), cos(s1$mu_rad - s2$mu_rad)) * day_per_rad
      se_pair <- sqrt((s1$sd_doy^2 / s1$n) + (s2$sd_doy^2 / s2$n))
      pair_df <- rbind(s1$sub_df, s2$sub_df); pair_grp <- c(rep(tx1, s1$n), rep(tx2, s2$n))
      has_hdoy <- ("hopkins_adjusted_doy" %in% names(pair_df)) && all(is.finite(suppressWarnings(as.numeric(pair_df$hopkins_adjusted_doy))))
      pair_th <- if (isTRUE(has_hdoy)) to_rad(pair_df$hopkins_adjusted_doy) else as.numeric(pair_df$delta_doy) * 2.0 * pi / 365.25
      pair_res <- if (exists("run_circular_multigroup_phase_test", mode = "function")) run_circular_multigroup_phase_test(pair_th, pair_grp, B = B, seed = seed) else list(watson_f = NA_real_, perm_p = NA_real_)
      list(
        tx1 = tx1, tx2 = tx2, n_pair = as.integer(s1$n + s2$n), diff_doy = diff_doy, se_pair = se_pair,
        wf = if (!is.null(pair_res$watson_f) && is.finite(pair_res$watson_f)) as.numeric(pair_res$watson_f) else NA_real_,
        raw_p = if (!is.null(pair_res$perm_p) && is.finite(pair_res$perm_p)) as.numeric(pair_res$perm_p) else NA_real_
      )
    }
    raw_pairs <- utils::combn(active_tx, 2L, simplify = FALSE)
    nw_p <- if (length(raw_pairs) >= 2L && exists("packera_mclapply", mode = "function") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, length(raw_pairs)) else 1L
    pair_recs <- Filter(Negate(is.null), if (nw_p > 1L) packera_mclapply(raw_pairs, eval_pair, mc.cores = nw_p) else lapply(raw_pairs, eval_pair))
  }
  if (length(pair_recs) > 0L) {
    adj_ps <- stats::p.adjust(vapply(pair_recs, function(r) r$raw_p, numeric(1)), method = "BH")
    for (idx in seq_along(pair_recs)) {
      pr <- pair_recs[[idx]]; p_bh <- adj_ps[idx]
      p_int <- if (isTRUE(is_unif)) sprintf("Gated under Rayleigh uniform phase null (raw diff %.1f days)", abs(pr$diff_doy))
        else if (is.finite(p_bh) && p_bh < 0.05) sprintf("Statistically significant circular phase separation (%.1f days)", abs(pr$diff_doy))
        else "Phenologically synchronous / overlapping"
      rows[[length(rows) + 1L]] <- mk_row("Circular_Pairwise_PostHoc", paste(pr$tx1, "vs", pr$tx2), pr$n_pair, md = round(pr$diff_doy, 3),
                                          lo = round(pr$diff_doy - 1.96 * pr$se_pair, 3), hi = round(pr$diff_doy + 1.96 * pr$se_pair, 3),
                                          s1 = round(pr$diff_doy, 3), s2 = if (is.finite(pr$wf)) round(pr$wf, 3) else NA_real_,
                                          pv = if (is.finite(p_bh)) sprintf("%.4e", p_bh) else "NA", interp = p_int)
    }
  }
  out_df <- do.call(rbind, rows); rownames(out_df) <- NULL; out_df
}
