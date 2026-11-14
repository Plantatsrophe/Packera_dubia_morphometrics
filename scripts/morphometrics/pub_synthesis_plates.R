#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_synthesis_plates.R
# Project: Packera Genus-Wide Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Step 05 ecological niche & geographic map publication figures, Step 06
#   integrative synthesis and complex discovery plates, pure-R
#   Gaussian Mixture Model EM solver, and individual-level trait mergers.
# ==============================================================================

# ==============================================================================
# Step 06 Figures:
#   - spatial_rf_niche_importance.png (300-DPI companion to .pdf)
#   - fig05_ecological_niche_phenology.{pdf,png} (4-panel ecological niche plate)
#   - fig05d_geographic_map_and_hopkins_phenology.{pdf,png} (US state map + Hopkins)
# ==============================================================================
compute_dynamic_na_map_bounds <- function(lon_v, lat_v) {
  lon_num <- suppressWarnings(as.numeric(lon_v))
  lat_num <- suppressWarnings(as.numeric(lat_v))
  valid_geo <- is.finite(lon_num) & is.finite(lat_num) &
               (lon_num >= -170.0) & (lon_num <= -50.0) &
               (lat_num >= 14.0) & (lat_num <= 72.0)
  if (!any(valid_geo)) {
    return(list(valid_geo = valid_geo, xlim = c(-126.0, -66.0), ylim = c(24.0, 50.0)))
  }
  lon_obs <- range(lon_num[valid_geo], na.rm = TRUE)
  lat_obs <- range(lat_num[valid_geo], na.rm = TRUE)
  lon_pad <- 0.12 * max(diff(lon_obs), 1e-3)
  lat_pad <- 0.12 * max(diff(lat_obs), 1e-3)
  xlim_map <- c(lon_obs[1L] - lon_pad, lon_obs[2L] + lon_pad)
  ylim_map <- c(lat_obs[1L] - lat_pad, lat_obs[2L] + lat_pad)
  if (diff(xlim_map) < 8.0) {
    lon_mid <- mean(lon_obs)
    xlim_map <- c(lon_mid - 4.0, lon_mid + 4.0)
  }
  if (diff(ylim_map) < 6.0) {
    lat_mid <- mean(lat_obs)
    ylim_map <- c(lat_mid - 3.0, lat_mid + 3.0)
  }
  xlim_map <- c(max(-172.0, xlim_map[1L]), min(-48.0, xlim_map[2L]))
  ylim_map <- c(max(12.0, ylim_map[1L]), min(74.0, ylim_map[2L]))
  list(valid_geo = valid_geo, xlim = xlim_map, ylim = ylim_map)
}


# ==============================================================================
# Step 06 Integrative Taxonomic Synthesis & Complex Network 300-DPI PNG Wrappers
#   - Stream 3: "3. Topological ECT & OSR"
#   - Consensus weighting: probabilistic_consensus_score with position_fill(reverse = TRUE)
# ==============================================================================
.ensure_step06_fig_modules <- function() {
  if (!exists("render_step06a_integrative_synthesis_png", mode = "function")) {
    for (.f in c("scripts/figures/theme_packera_pub.R", "scripts/figures/step06_figures.R")) {
      if (file.exists(.f)) source(.f)
    }
  }
}

export_step05d_geographic_and_hopkins <- function(df, out_plot_path = NULL, ...) {
  # Consolidated into canonical fig05b_hopkins_circular_phenology.png & fig05c_sdm_suitability_and_distribution_map.png
  invisible("")
}

generate_step06_pub_figures <- function(triage_df, out_fig_pdf = NULL) {
  # Consolidated into canonical fig06a_integrative_taxonomic_synthesis.png
  invisible("")
}

render_step06_synthesis_plate_pdf <- function(df, out_pdf, cda_summary = NULL, focal_taxon = "") {
  if (identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  .ensure_step06_fig_modules()
  out_png <- if (!is.null(out_pdf) && nzchar(out_pdf)) normalize_png_path(out_pdf) else "outputs/figures/fig06a_integrative_taxonomic_synthesis.png"
  render_step06a_integrative_synthesis_png(
    df = df,
    out_plot_path = out_png,
    cda_summary = cda_summary,
    focal_taxon = focal_taxon
  )
}

render_step06_complex_discovery_pdf <- function(W, part_res, out_pdf, min_affinity = 0.03) {
  if (identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  .ensure_step06_fig_modules()
  out_png <- if (!is.null(out_pdf) && nzchar(out_pdf)) normalize_png_path(out_pdf) else "outputs/figures/fig06b_taxon_complex_delimitation_network.png"
  cdf_in <- if (is.list(part_res) && !is.null(part_res$complexes_df)) part_res$complexes_df else as.data.frame(part_res, stringsAsFactors = FALSE)
  render_step06b_complex_network_png(
    complexes_df = cdf_in,
    affinity_mat = W,
    out_complex_path = out_png,
    min_affinity = min_affinity
  )
}

# ==============================================================================
# Shared Cross-Stage Morphometric Helpers (GMM EM Fallback & Composite-Key Joins)
# ==============================================================================
fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150, X_eval = NULL, seed = 42L) {
  if (is.null(X_eval)) X_eval <- X
  N <- nrow(X); p <- ncol(X); best_bic <- -Inf; best_model <- NULL
  base_seed <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  bic_table <- data.frame(K = seq_len(max_k), LogLik = NA_real_, n_params = NA_integer_, BIC = NA_real_)
  for (k in seq_len(max_k)) {
    set.seed(base_seed + k)
    if (k == 1) {
      mu <- matrix(colMeans(X), nrow = 1); sigma <- list(stats::cov(X) + diag(1e-6, p)); pi_k <- c(1.0)
    } else {
      km <- stats::kmeans(X, centers = k, nstart = 10, iter.max = 30); mu <- km$centers
      pi_k <- as.vector(table(factor(km$cluster, levels = seq_len(k)))) / N
      sigma <- lapply(seq_len(k), function(c) {
        sub_x <- X[km$cluster == c, , drop = FALSE]
        if (nrow(sub_x) > p) stats::cov(sub_x) + diag(1e-5, p) else diag(1e-2, p)
      })
    }
    loglik_old <- -Inf; resp <- matrix(0, nrow = N, ncol = k)
    for (iter in seq_len(max_iter)) {
      log_dens <- matrix(0, nrow = N, ncol = k)
      for (c in seq_len(k)) {
        diff_x <- t(t(X) - mu[c, ]); sig_inv <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]])))
        quad <- rowSums((diff_x %*% sig_inv) * diff_x)
        log_dens[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(max(det(sigma[[c]]), 1e-12)) + quad)
      }
      max_l <- apply(log_dens, 1, max); log_sum_exp <- max_l + log(rowSums(exp(log_dens - max_l)))
      resp <- exp(log_dens - log_sum_exp); loglik <- sum(log_sum_exp)
      if (abs(loglik - loglik_old) < tol) break
      loglik_old <- loglik; N_k <- colSums(resp); pi_k <- N_k / N
      for (c in seq_len(k)) if (N_k[c] > 1e-6) {
        mu[c, ] <- colSums(resp[, c] * X) / N_k[c]; diff_x <- t(t(X) - mu[c, ])
        sigma[[c]] <- (t(diff_x) %*% (resp[, c] * diff_x)) / N_k[c] + diag(1e-5, p)
      }
    }
    n_params <- (k - 1) + k * p + k * (p * (p + 1) / 2); bic_val <- 2 * loglik - n_params * log(N)
    bic_table$LogLik[k] <- round(loglik, 2); bic_table$n_params[k] <- n_params; bic_table$BIC[k] <- round(bic_val, 2)
    if (bic_val > best_bic) {
      best_bic <- bic_val
      ld_ev <- matrix(0, nrow = nrow(X_eval), ncol = k)
      for (c in seq_len(k)) {
        dx_e <- t(t(X_eval) - mu[c, ]); si_e <- tryCatch(solve(sigma[[c]]), error = function(e) diag(1 / diag(sigma[[c]])))
        ld_ev[, c] <- log(max(pi_k[c], 1e-12)) - 0.5 * (p * log(2 * pi) + log(max(det(sigma[[c]]), 1e-12)) + rowSums((dx_e %*% si_e) * dx_e))
      }
      ml_e <- apply(ld_ev, 1, max); resp_ev <- exp(ld_ev - (ml_e + log(rowSums(exp(ld_ev - ml_e)))))
      best_model <- list(k = k, classification = apply(resp_ev, 1, which.max), uncertainty = 1 - apply(resp_ev, 1, max), z = resp_ev, bic = bic_val)
    }
  }
  return(list(best = best_model, bic_table = bic_table))
}

merge_manifest_petiole_traits_by_individual <- function(agg_df, man_cand) {
  if (is.null(agg_df) || nrow(agg_df) == 0L || is.null(man_cand) || !nzchar(man_cand) || !file.exists(man_cand)) return(agg_df)
  target_pcols <- c("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "blade_cie_l_norm", "blade_cie_l_sd")
  need_cols <- target_pcols[vapply(target_pcols, function(cl) !(cl %in% names(agg_df)) || any(is.na(agg_df[[cl]])), logical(1))]
  if (length(need_cols) == 0L) return(agg_df)
  man_df <- tryCatch(utils::read.csv(man_cand, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(man_df) || nrow(man_df) == 0L || !("catalogNumber" %in% names(man_df))) return(agg_df)
  p_cols <- intersect(target_pcols, names(man_df))
  if (length(p_cols) == 0L) return(agg_df)
  man_df$catalogNumber <- as.character(man_df$catalogNumber)
  agg_df$catalogNumber <- as.character(agg_df$catalogNumber)
  q75_fn <- function(x) {
    v <- suppressWarnings(as.numeric(x)); v <- v[is.finite(v)]
    if (length(v) > 0L) round(stats::quantile(v, probs = 0.75, na.rm = TRUE, names = FALSE), 4) else NA_real_
  }
  med_fn <- function(x) {
    v <- suppressWarnings(as.numeric(x)); v <- v[is.finite(v)]
    if (length(v) > 0L) round(stats::median(v), 4) else NA_real_
  }
  if ("leaf_zone" %in% names(man_df)) {
    lz_low <- tolower(trimws(as.character(man_df$leaf_zone)))
    lz_low[is.na(lz_low) | !nzchar(lz_low)] <- "basal"
    if (any(lz_low != "cauline")) {
      grp_k <- if ("plant_individual_id" %in% names(man_df)) {
        mp_tmp <- ifelse(is.na(man_df$plant_individual_id), 0L, suppressWarnings(as.integer(man_df$plant_individual_id)))
        mp_tmp[is.na(mp_tmp)] <- 0L
        paste0(man_df$catalogNumber, "__", mp_tmp)
      } else man_df$catalogNumber
      has_basal_grp <- tapply(lz_low != "cauline", grp_k, any)
      keep_row <- (lz_low != "cauline") | !has_basal_grp[grp_k]
      man_df <- man_df[which(keep_row), , drop = FALSE]
    }
  }
  has_pid <- ("plant_individual_id" %in% names(man_df)) && ("plant_individual_id" %in% names(agg_df))
  q75_cols <- setdiff(p_cols, c("blade_cie_l_norm", "blade_cie_l_sd"))
  med_cols <- intersect(p_cols, c("blade_cie_l_norm", "blade_cie_l_sd"))
  agg_by_fn <- function(by_list) {
    parts <- list()
    if (length(q75_cols) > 0L) parts[[length(parts) + 1L]] <- stats::aggregate(man_df[, q75_cols, drop = FALSE], by = by_list, FUN = q75_fn)
    if (length(med_cols) > 0L) parts[[length(parts) + 1L]] <- stats::aggregate(man_df[, med_cols, drop = FALSE], by = by_list, FUN = med_fn)
    if (length(parts) == 1L) return(parts[[1L]])
    merge(parts[[1L]], parts[[2L]], by = names(by_list), all = TRUE, sort = FALSE)
  }
  agg_cat <- agg_by_fn(list(catalogNumber = man_df$catalogNumber))
  agg_pid <- if (has_pid) {
    mpid <- ifelse(is.na(man_df$plant_individual_id), 0L, suppressWarnings(as.integer(man_df$plant_individual_id)))
    mpid[is.na(mpid)] <- 0L
    agg_by_fn(list(catalogNumber = man_df$catalogNumber, plant_individual_id = mpid))
  } else NULL
  idx_pid <- if (has_pid && !is.null(agg_pid)) {
    apid <- ifelse(is.na(agg_df$plant_individual_id), 0L, suppressWarnings(as.integer(agg_df$plant_individual_id)))
    apid[is.na(apid)] <- 0L
    match(paste0(agg_df$catalogNumber, "__", apid), paste0(agg_pid$catalogNumber, "__", agg_pid$plant_individual_id))
  } else rep(NA_integer_, nrow(agg_df))
  idx_cat <- match(agg_df$catalogNumber, agg_cat$catalogNumber)
  single_cat_in_man <- if (has_pid) {
    mpid_all <- ifelse(is.na(man_df$plant_individual_id), 0L, suppressWarnings(as.integer(man_df$plant_individual_id)))
    mpid_all[is.na(mpid_all)] <- 0L
    n_pids <- tapply(mpid_all, man_df$catalogNumber, function(z) length(unique(z)))
    !(agg_df$catalogNumber %in% names(n_pids[n_pids > 1L]))
  } else rep(TRUE, nrow(agg_df))
  for (cl in intersect(need_cols, p_cols)) {
    val_out <- if (cl %in% names(agg_df)) suppressWarnings(as.numeric(agg_df[[cl]])) else rep(NA_real_, nrow(agg_df))
    if (!is.null(agg_pid)) {
      hit_p <- is.na(val_out) & !is.na(idx_pid)
      if (any(hit_p)) val_out[hit_p] <- agg_pid[[cl]][idx_pid[hit_p]]
    }
    hit_c <- is.na(val_out) & !is.na(idx_cat) & single_cat_in_man
    if (any(hit_c)) val_out[hit_c] <- agg_cat[[cl]][idx_cat[hit_c]]
    agg_df[[cl]] <- val_out
  }
  agg_df
}

merge_reproductive_traits_by_individual <- function(agg_df, repro_path) {
  if (is.null(agg_df) || nrow(agg_df) == 0L || is.null(repro_path) || !nzchar(repro_path) || !file.exists(repro_path)) return(agg_df)
  repro_df <- tryCatch(utils::read.csv(repro_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(repro_df) || nrow(repro_df) == 0L || !("catalogNumber" %in% names(repro_df))) return(agg_df)
  if ("plant_height_mm" %in% names(repro_df) && !("plant_height_cm" %in% names(repro_df))) {
    repro_df$plant_height_cm <- round(suppressWarnings(as.numeric(repro_df$plant_height_mm)) / 10.0, 2)
  }
  if (!("capitula_count" %in% names(repro_df)) && ("capitulum_count" %in% names(repro_df))) {
    repro_df$capitula_count <- repro_df$capitulum_count
  }
  if ("capitula_count" %in% names(repro_df) && !("log1p_capitulum_count" %in% names(repro_df))) {
    cc_num <- suppressWarnings(as.numeric(repro_df$capitula_count))
    repro_df$log1p_capitulum_count <- ifelse(is.finite(cc_num) & cc_num >= 1, round(log1p(cc_num), 5), NA_real_)
  }
  cand_r <- c("AR_inv", "capitula_count", "log1p_capitulum_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm",
              "involucre_width_mm", "plant_height_mm", "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem", "petiole_to_blade_ratio")
  r_val_cols <- intersect(names(repro_df), cand_r)
  r_val_cols <- r_val_cols[vapply(r_val_cols, function(cl) !(cl %in% names(agg_df)) || any(is.na(agg_df[[cl]])), logical(1))]
  if (length(r_val_cols) == 0L) return(agg_df)
  repro_df$catalogNumber <- as.character(repro_df$catalogNumber)
  agg_df$catalogNumber <- as.character(agg_df$catalogNumber)
  has_pid <- ("plant_individual_id" %in% names(repro_df)) && ("plant_individual_id" %in% names(agg_df))
  if (has_pid) {
    rpid <- ifelse(is.na(repro_df$plant_individual_id), 0L, suppressWarnings(as.integer(repro_df$plant_individual_id))); rpid[is.na(rpid)] <- 0L
    apid <- ifelse(is.na(agg_df$plant_individual_id), 0L, suppressWarnings(as.integer(agg_df$plant_individual_id))); apid[is.na(apid)] <- 0L
    r_sub_pid <- repro_df[!duplicated(paste0(repro_df$catalogNumber, "__", rpid)), , drop = FALSE]
    rpid_sub <- ifelse(is.na(r_sub_pid$plant_individual_id), 0L, suppressWarnings(as.integer(r_sub_pid$plant_individual_id))); rpid_sub[is.na(rpid_sub)] <- 0L
    idx_pid <- match(paste0(agg_df$catalogNumber, "__", apid), paste0(r_sub_pid$catalogNumber, "__", rpid_sub))
    n_pids <- tapply(rpid, repro_df$catalogNumber, function(z) length(unique(z)))
    single_cat <- !(agg_df$catalogNumber %in% names(n_pids[n_pids > 1L]))
  } else {
    r_sub_pid <- NULL; idx_pid <- rep(NA_integer_, nrow(agg_df)); single_cat <- rep(TRUE, nrow(agg_df))
  }
  r_sub_cat <- repro_df[!duplicated(repro_df$catalogNumber), , drop = FALSE]
  idx_cat <- match(agg_df$catalogNumber, r_sub_cat$catalogNumber)
  for (cl in r_val_cols) {
    cur <- if (cl %in% names(agg_df)) agg_df[[cl]] else rep(NA, nrow(agg_df))
    if (!is.null(r_sub_pid)) {
      hit_p <- is.na(cur) & !is.na(idx_pid)
      if (any(hit_p)) cur[hit_p] <- r_sub_pid[[cl]][idx_pid[hit_p]]
    }
    hit_c <- is.na(cur) & !is.na(idx_cat) & single_cat
    if (any(hit_c)) cur[hit_c] <- r_sub_cat[[cl]][idx_cat[hit_c]]
    agg_df[[cl]] <- cur
  }
  agg_df
}

render_focal_revision_synthesis_plate_png <- function(df, out_png, full_df = NULL, focal_taxon = "") {
  if (identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) return(invisible(""))
  .ensure_step06_fig_modules()
  fdir <- if (!is.null(out_png) && nzchar(out_png)) dirname(out_png) else "outputs/figures"
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  foc_df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(foc_df) == 0L) return(invisible(""))
  foc_name <- if (nzchar(trimws(as.character(focal_taxon)))) trimws(as.character(focal_taxon)) else as.character(stats::na.omit(foc_df$species_standardized)[1L])
  sub_call <- if ("taxonomic_treatment_subtype" %in% names(foc_df)) as.character(foc_df$taxonomic_treatment_subtype) else rep("Confirmed_Species", nrow(foc_df))
  sub_call[is.na(sub_call) | !nzchar(trimws(sub_call)) | !(sub_call %in% names(PACKERA_TREATMENT_PALETTE))] <- "Confirmed_Species"
  foc_df$treatment_call <- factor(sub_call, levels = intersect(names(PACKERA_TREATMENT_PALETTE), unique(sub_call)))

  # Panel A: Focal vs. Complex Neighbors Basal & Cauline Morphospace
  use_can <- ("can1" %in% names(foc_df)) && sum(is.finite(suppressWarnings(as.numeric(foc_df$can1)))) >= 3L
  mx_col <- if (use_can) "can1" else if ("PC1" %in% names(foc_df)) "PC1" else "ECT_PC1"
  my_col <- if (use_can && "can2" %in% names(foc_df) && isTRUE(stats::sd(suppressWarnings(as.numeric(foc_df$can2)), na.rm = TRUE) > 1e-6)) "can2" else if ("PC2" %in% names(foc_df)) "PC2" else if ("ECT_PC1" %in% names(foc_df)) "ECT_PC1" else mx_col
  foc_df$mx <- suppressWarnings(as.numeric(foc_df[[mx_col]]))
  foc_df$my <- suppressWarnings(as.numeric(foc_df[[my_col]]))
  nb_df <- data.frame()
  if (!is.null(full_df) && is.data.frame(full_df) && nrow(full_df) > nrow(foc_df) && all(c(mx_col, my_col) %in% names(full_df))) {
    nb_raw <- as.data.frame(full_df, stringsAsFactors = FALSE)
    nb_raw$mx <- suppressWarnings(as.numeric(nb_raw[[mx_col]])); nb_raw$my <- suppressWarnings(as.numeric(nb_raw[[my_col]]))
    nb_df <- nb_raw[as.character(nb_raw$species_standardized) != foc_name & is.finite(nb_raw$mx) & is.finite(nb_raw$my), , drop = FALSE]
  }
  foc_m <- foc_df[is.finite(foc_df$mx) & is.finite(foc_df$my), , drop = FALSE]
  p_a <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4)
  if (nrow(nb_df) > 0L) {
    nb_cent <- stats::aggregate(cbind(mx, my) ~ species_standardized, data = nb_df, FUN = stats::median)
    nb_cent$lbl <- abbreviate_icn_taxon(nb_cent$species_standardized, abbreviate_genus = TRUE)
    p_a <- p_a + ggplot2::geom_point(data = nb_df, ggplot2::aes(x = .data$mx, y = .data$my), shape = 1, size = 1.6, stroke = 0.55, color = "#94A3B8", alpha = 0.42)
    p_a <- add_centroid_repel_labels_gg(p_a, utils::head(nb_cent, 8L), "mx", "my", "lbl", size = 2.45)
  }
  p_a <- p_a +
    ggplot2::geom_point(data = foc_m, ggplot2::aes(x = .data$mx, y = .data$my, color = .data$treatment_call, shape = .data$treatment_call), size = 2.4, stroke = 0.9, alpha = 0.88) +
    ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE, labels = function(x) gsub("_", " ", x), name = "Treatment") +
    ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES, labels = function(x) gsub("_", " ", x), name = "Treatment") +
    theme_packera_pub(base_size = 9.4) +
    ggplot2::labs(
      title = sprintf("A. %s vs. Complex Neighbors Morphospace", foc_name),
      subtitle = "Grey circles = congeneric complex; colored hollow symbols = focal treatments (* = Type)",
      x = if (mx_col == "can1") "Canonical Variate 1 (Can1)" else format_morph_trait_labels(mx_col),
      y = if (my_col == "can2") "Canonical Variate 2 (Can2)" else format_morph_trait_labels(my_col)
    )
  if ("is_type_specimen" %in% names(foc_m) && any(!is.na(foc_m$is_type_specimen) & as.logical(foc_m$is_type_specimen))) {
    t_m <- foc_m[!is.na(foc_m$is_type_specimen) & as.logical(foc_m$is_type_specimen), , drop = FALSE]
    p_a <- p_a + ggplot2::geom_point(data = t_m, ggplot2::aes(x = .data$mx, y = .data$my), shape = 8, size = 3.8, stroke = 1.2, color = "#111827")
  }

  # Panel B: EPSG:5070 Geographic Treatment & Type-Control Map
  lon_col <- if ("decimalLongitude" %in% names(foc_df)) "decimalLongitude" else "longitude"
  lat_col <- if ("decimalLatitude" %in% names(foc_df)) "decimalLatitude" else "latitude"
  m_ctx <- build_na_projected_basemap_gg(foc_df[[lon_col]], foc_df[[lat_col]])
  pts_g <- project_lonlat_points_r(foc_df, m_ctx, lon_col = lon_col, lat_col = lat_col)
  p_b <- m_ctx$plot
  if (nrow(pts_g) > 0L) {
    p_b <- p_b +
      ggplot2::geom_point(data = pts_g, ggplot2::aes(x = .data$map_x, y = .data$map_y, color = .data$treatment_call, shape = .data$treatment_call), size = 2.3, stroke = 0.85, alpha = 0.86, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE) +
      ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES)
    if ("is_type_specimen" %in% names(pts_g) && any(!is.na(pts_g$is_type_specimen) & as.logical(pts_g$is_type_specimen))) {
      tg <- pts_g[!is.na(pts_g$is_type_specimen) & as.logical(pts_g$is_type_specimen), , drop = FALSE]
      p_b <- p_b + ggplot2::geom_point(data = tg, ggplot2::aes(x = .data$map_x, y = .data$map_y), shape = 8, size = 3.8, stroke = 1.2, color = "#111827")
    }
  }
  p_b <- p_b + ggplot2::labs(title = "B. Geographic Distribution & Type Localities (EPSG:5070)", subtitle = "Projected CONUS Albers; star = nomenclatural type control")

  # Panel C: Rhizosphere Edaphic & Hopkins Phenology Envelope
  ed_col <- intersect(c("lithic_outcrop_index", "fine_bedrock_depth_cm", "fine_soil_ph", "soil_ph"), names(foc_df))[1L]
  ph_col <- intersect(c("delta_doy_resid", "hopkins_adjusted_doy", "doy"), names(foc_df))[1L]
  p_c <- if (!is.na(ed_col) && !is.na(ph_col)) {
    foc_df$ex_v <- suppressWarnings(as.numeric(foc_df[[ed_col]]))
    foc_df$ph_v <- suppressWarnings(as.numeric(foc_df[[ph_col]]))
    ec <- foc_df[is.finite(foc_df$ex_v) & is.finite(foc_df$ph_v), , drop = FALSE]
    pe <- ggplot2::ggplot()
    if (!is.null(full_df) && is.data.frame(full_df) && all(c(ed_col, ph_col) %in% names(full_df))) {
      nb_e <- as.data.frame(full_df, stringsAsFactors = FALSE)
      nb_e <- nb_e[as.character(nb_e$species_standardized) != foc_name, , drop = FALSE]
      nb_e$ex_v <- suppressWarnings(as.numeric(nb_e[[ed_col]])); nb_e$ph_v <- suppressWarnings(as.numeric(nb_e[[ph_col]]))
      nb_e <- nb_e[is.finite(nb_e$ex_v) & is.finite(nb_e$ph_v), , drop = FALSE]
      if (nrow(nb_e) > 0L) pe <- pe + ggplot2::geom_point(data = nb_e, ggplot2::aes(x = .data$ex_v, y = .data$ph_v), shape = 1, size = 1.6, stroke = 0.5, color = "#94A3B8", alpha = 0.35)
    }
    if (ph_col == "delta_doy_resid") pe <- pe + ggplot2::geom_hline(yintercept = c(-21, 0, 21), linetype = c("dashed", "dotted", "dashed"), color = c("#D55E00", "#9CA3AF", "#D55E00"), linewidth = 0.5)
    pe +
      ggplot2::geom_point(data = ec, ggplot2::aes(x = .data$ex_v, y = .data$ph_v, color = .data$treatment_call, shape = .data$treatment_call), size = 2.3, stroke = 0.85, alpha = 0.85, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE) +
      ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES) +
      theme_packera_pub(base_size = 9.4) +
      ggplot2::labs(
        title = "C. Rhizosphere Edaphic & Hopkins Phenology Envelope",
        subtitle = "Focal vouchers vs. complex background along edaphic & phenological axes",
        x = if (exists("format_step05_predictor_labels", mode = "function")) format_step05_predictor_labels(ed_col) else ed_col,
        y = if (ph_col == "delta_doy_resid") "Hopkins Phenological Anomaly (Delta DOY, d)" else "Hopkins-Adjusted Anthesis DOY"
      )
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "C. Rhizosphere Edaphic & Phenology Envelope")
  }

  # Panel D: Focal Taxon Treatment Breakdown & Consensus Score Summary
  trt_tab <- as.data.frame(table(treatment_call = foc_df$treatment_call), stringsAsFactors = FALSE)
  trt_tab$pct <- 100.0 * trt_tab$Freq / max(1L, sum(trt_tab$Freq))
  trt_tab <- trt_tab[order(trt_tab$Freq), , drop = FALSE]
  trt_tab$lbl <- factor(gsub("_", " ", trt_tab$treatment_call), levels = gsub("_", " ", trt_tab$treatment_call))
  trt_tab$txt <- sprintf("%d (%.1f%%)", trt_tab$Freq, trt_tab$pct)
  pcs <- if ("probabilistic_consensus_score" %in% names(foc_df)) suppressWarnings(as.numeric(foc_df$probabilistic_consensus_score)) else numeric(0)
  pcs_med <- if (any(is.finite(pcs))) sprintf("Median consensus = %.3f", stats::median(pcs, na.rm = TRUE)) else ""
  n_basal_src <- if ("morph_tier_source" %in% names(foc_df)) sum(foc_df$morph_tier_source == "Basal_Rosette", na.rm = TRUE) else nrow(foc_df)
  p_d <- ggplot2::ggplot(trt_tab, ggplot2::aes(x = .data$Freq, y = .data$lbl, fill = .data$treatment_call)) +
    ggplot2::geom_col(width = 0.62, color = "white", linewidth = 0.35, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$txt), hjust = -0.10, size = 2.8, color = "#111827") +
    ggplot2::scale_fill_manual(values = PACKERA_TREATMENT_PALETTE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) +
    theme_packera_pub(base_size = 9.4) +
    ggplot2::labs(
      title = sprintf("D. %s Treatment Breakdown (N = %d)", foc_name, nrow(foc_df)),
      subtitle = paste(c(sprintf("Basal_Rosette = %d / %d", n_basal_src, nrow(foc_df)), if (nzchar(pcs_med)) pcs_med else NULL), collapse = " | "),
      x = "Voucher Count", y = NULL
    )

  plate_f <- ((p_a | p_b) / (p_c | p_d)) +
    patchwork::plot_layout(heights = c(1.12, 1.0)) +
    patchwork::plot_annotation(
      title = sprintf("Integrative Taxonomic Revision Plate: %s", foc_name),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5))
    )
  save_packera_png_300(plate_f, normalize_png_path(out_png), width = 14.4, height = 10.6)
}
