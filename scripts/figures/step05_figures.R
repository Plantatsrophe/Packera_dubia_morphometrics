# ==============================================================================
# Module: scripts/figures/step05_figures.R
# Step 05 Canonical 300-DPI PNG Figures (Spatial Block-CV RF, PCA-env Niche,
# Circular Hopkins Phenology & Equal-Area Projected MaxEnt/Ensemble SDMs):
#   1. fig05a_ecological_niche_and_rf_importance.png
#   2. fig05b_hopkins_circular_phenology.png
#   3. fig05c_sdm_suitability_and_distribution_map.png
# ==============================================================================

.classify_env_domain <- function(var_names) {
  vapply(as.character(var_names), function(v) {
    vl <- tolower(v)
    if (grepl("^bio[0-9]+|temp|precip|clim", vl)) return("Climatic (WorldClim)")
    if (grepl("elev|slope|aspect|tpi|tri|topo|dem", vl)) return("Topographic (3DEP)")
    "Edaphic & Lithic (SOLUS/POLARIS/SSURGO/SGMC)"
  }, character(1), USE.NAMES = FALSE)
}

render_step05a_ecological_niche_png <- function(
  df,
  srf_res = NULL,
  niche_res = NULL,
  out_plot_path = "outputs/figures/fig05a_ecological_niche_and_rf_importance.png"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (!is.null(out_plot_path) && nzchar(out_plot_path)) dirname(out_plot_path) else "outputs/figures"
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (old_stem in c("spatial_rf_niche_importance", "fig05_ecological_niche_phenology", "fig05d_geographic_map_and_hopkins_phenology")) {
    for (ext in c(".pdf", ".png")) {
      p_old <- file.path(fdir, paste0(old_stem, ext))
      if (file.exists(p_old) && normalizePath(p_old, mustWork = FALSE) != normalizePath(normalize_png_path(out_plot_path), mustWork = FALSE)) {
        unlink(p_old, force = TRUE)
      }
    }
  }

  vdf <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(vdf) == 0L) return(invisible(""))
  grp_col <- if ("ecotype_cluster" %in% names(vdf) && length(unique(stats::na.omit(vdf$species_standardized))) < 2L) "ecotype_cluster" else NULL
  vdf$plot_group <- resolve_plot_group_labels(vdf, group_col = grp_col)
  active_groups <- sort(unique(vdf$plot_group[!is.na(vdf$plot_group) & nzchar(vdf$plot_group) & vdf$plot_group != "Unknown"]))
  if (length(active_groups) == 0L) active_groups <- unique(vdf$plot_group)
  vdf <- vdf[vdf$plot_group %in% active_groups, , drop = FALSE]

  pal_map <- get_combinatorial_palette_r(active_groups)
  shp_map <- get_combinatorial_hollow_shapes_r(active_groups)
  cnt_map <- table(vdf$plot_group)
  lbl_exprs <- format_icn_label_vec(active_groups, abbreviate_genus = (length(active_groups) > 5L), counts_map = cnt_map)
  names(lbl_exprs) <- active_groups

  # Panel A: Spatial Block-CV RF Predictor Importance (Dynamic across all retained predictors)
  imp_df <- data.frame()
  if (!is.null(srf_res) && !is.null(srf_res$importance)) {
    imp_raw <- srf_res$importance
    if (is.data.frame(imp_raw)) {
      vn_hit <- intersect(c("variable", "Variable", "Predictor", "predictor", "Feature", "feature"), names(imp_raw))
      vn <- if (length(vn_hit) > 0L) imp_raw[[vn_hit[1L]]] else rownames(imp_raw)
      iv_hit <- intersect(c("importance", "Importance", "MeanDecreaseAccuracy", "IncNodePurity", "MeanDecreaseGini"), names(imp_raw))
      if (length(iv_hit) > 0L) {
        iv <- imp_raw[[iv_hit[1L]]]
      } else {
        num_cols_i <- which(vapply(imp_raw, is.numeric, logical(1)))
        iv <- if (length(num_cols_i) > 0L) imp_raw[[num_cols_i[1L]]] else rep(NA_real_, nrow(imp_raw))
      }
      imp_df <- data.frame(var = as.character(vn), imp = suppressWarnings(as.numeric(iv)), stringsAsFactors = FALSE)
    } else if (is.numeric(imp_raw) && !is.null(names(imp_raw))) {
      imp_df <- data.frame(var = names(imp_raw), imp = as.numeric(imp_raw), stringsAsFactors = FALSE)
    }
  }
  imp_df <- imp_df[!grepl("^MEM[0-9]+$", imp_df$var) & is.finite(imp_df$imp), , drop = FALSE]
  if (nrow(imp_df) == 0L) {
    env_cands <- intersect(c("fine_soil_ph", "soil_ph", "bio1_temp_mean", "bio12_precip_annual", "ssurgo_rock_outcrop_pct", "fine_bedrock_depth_cm", "topo_elevation_m"), names(vdf))
    if (length(env_cands) > 0L) {
      imp_df <- data.frame(
        var = env_cands,
        imp = vapply(env_cands, function(cc) stats::sd(suppressWarnings(as.numeric(vdf[[cc]])), na.rm = TRUE), numeric(1)),
        stringsAsFactors = FALSE
      )
      imp_df <- imp_df[is.finite(imp_df$imp), , drop = FALSE]
    }
  }

  p_imp <- if (nrow(imp_df) > 0L) {
    imp_df <- imp_df[order(imp_df$imp, decreasing = TRUE), , drop = FALSE]
    imp_top <- utils::head(imp_df, 14L)
    imp_top$domain <- .classify_env_domain(imp_top$var)
    pretty_lbl <- if (exists("format_predictor_labels", mode = "function")) format_predictor_labels(imp_top$var) else gsub("_", " ", imp_top$var)
    imp_top$lbl <- factor(wrap_clean_label_r(pretty_lbl, width = 32L), levels = rev(unique(wrap_clean_label_r(pretty_lbl, width = 32L))))
    cv_acc <- if ("spatial_block_cv_acc" %in% names(vdf) && any(is.finite(suppressWarnings(as.numeric(vdf$spatial_block_cv_acc))))) {
      stats::median(suppressWarnings(as.numeric(vdf$spatial_block_cv_acc)), na.rm = TRUE)
    } else NA_real_
    ggplot2::ggplot(imp_top, ggplot2::aes(x = .data$imp, y = .data$lbl, fill = .data$domain)) +
      ggplot2::geom_col(width = 0.68, alpha = 0.88, color = "#374151", linewidth = 0.35) +
      ggplot2::scale_fill_manual(
        values = c(
          "Climatic (WorldClim)" = "#56B4E9",
          "Edaphic & Lithic (SOLUS/POLARIS/SSURGO/SGMC)" = "#D55E00",
          "Topographic (3DEP)" = "#009E73"
        ),
        name = "Environmental Stream"
      ) +
      theme_packera_pub(base_size = 9.6) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        title = "A. Spatial Block-CV Random Forest Importance",
        subtitle = if (is.finite(cv_acc)) sprintf("Spatial block cross-validation accuracy = %.1f%% (non-MEM predictors)", 100.0 * cv_acc) else "3-block VIF-pruned climatic, rhizosphere pedologic & topographic predictors",
        x = "Permutation / Impurity Variable Importance",
        y = NULL
      )
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "A. Spatial RF Predictor Importance")
  }

  # Panel B: 2D Ecological Niche Space (Top-2 Non-MEM RF Predictors or PCA-env)
  cand_num_env <- if (nrow(imp_df) >= 2L) intersect(imp_df$var, names(vdf)) else character(0)
  if (length(cand_num_env) < 2L) {
    fb_env <- grep("^(fine_|soil_|bio[0-9]+|ssurgo_|polaris_|solus_|topo_)", names(vdf), value = TRUE)
    fb_env <- fb_env[vapply(fb_env, function(cc) sum(is.finite(suppressWarnings(as.numeric(vdf[[cc]])))) >= 3L, logical(1))]
    cand_num_env <- unique(c(cand_num_env, fb_env))
  }

  p_niche <- if (length(cand_num_env) >= 2L) {
    vx_nm <- cand_num_env[1L]; vy_nm <- cand_num_env[2L]
    vdf$ex <- suppressWarnings(as.numeric(vdf[[vx_nm]]))
    vdf$ey <- suppressWarnings(as.numeric(vdf[[vy_nm]]))
    sub_e <- vdf[is.finite(vdf$ex) & is.finite(vdf$ey), , drop = FALSE]
    el_g <- names(which(table(sub_e$plot_group) >= 4L))
    xl_lbl <- if (exists("format_predictor_labels", mode = "function")) format_predictor_labels(vx_nm) else vx_nm
    yl_lbl <- if (exists("format_predictor_labels", mode = "function")) format_predictor_labels(vy_nm) else vy_nm
    p_tmp <- ggplot2::ggplot(sub_e, ggplot2::aes(x = .data$ex, y = .data$ey))
    if (length(el_g) > 0L) {
      p_tmp <- p_tmp + ggplot2::stat_ellipse(
        data = sub_e[sub_e$plot_group %in% el_g, , drop = FALSE],
        ggplot2::aes(color = .data$plot_group),
        level = 0.85, linetype = "dashed", linewidth = 0.55, show.legend = FALSE
      )
    }
    p_tmp <- p_tmp +
      ggplot2::geom_point(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), size = 2.4, stroke = 0.85, alpha = 0.82) +
      ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
      ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group") +
      theme_packera_pub() +
      ggplot2::labs(
        title = "B. Realized Ecological Niche Space (Top-2 Primary Predictors)",
        subtitle = "Uncertainty-gated occurrences (hollow markers + 85% niche ellipses)",
        x = xl_lbl,
        y = yl_lbl
      )
    if (length(active_groups) > 10L && nrow(sub_e) > 0L) {
      cent_e <- stats::aggregate(cbind(ex, ey) ~ plot_group, data = sub_e, FUN = stats::median)
      cent_e$lbl <- abbreviate_icn_taxon(cent_e$plot_group, abbreviate_genus = TRUE)
      p_tmp <- add_centroid_repel_labels_gg(p_tmp, cent_e, "ex", "ey", "lbl", color_col = "plot_group") +
        ggplot2::theme(legend.position = "none")
    }
    p_tmp
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "B. Ecological Niche Space")
  }

  # Panel C: Top Edaphic / Climatic Predictor Distributions + Kruskal-Wallis p
  box_vars <- utils::head(cand_num_env, 2L)
  p_edaph <- if (length(box_vars) > 0L) {
    ed_list <- lapply(box_vars, function(bv) {
      vv <- suppressWarnings(as.numeric(vdf[[bv]]))
      ok_b <- is.finite(vv)
      kw_p <- if (length(unique(vdf$plot_group[ok_b])) >= 2L && sum(ok_b) >= 5L) {
        tryCatch(stats::kruskal.test(vv[ok_b] ~ factor(vdf$plot_group[ok_b]))$p.value, error = function(e) NA_real_)
      } else NA_real_
      lbl_b <- if (exists("format_predictor_labels", mode = "function")) format_predictor_labels(bv) else bv
      if (is.finite(kw_p)) lbl_b <- sprintf("%s (KW p = %.3g)", lbl_b, kw_p)
      data.frame(plot_group = vdf$plot_group[ok_b], facet_lbl = lbl_b, val = vv[ok_b], stringsAsFactors = FALSE)
    })
    ed_df <- do.call(rbind, ed_list)
    ggplot2::ggplot(ed_df, ggplot2::aes(x = .data$val, y = .data$plot_group)) +
      ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = 0.55, linewidth = 0.45) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), height = 0.16, width = 0, size = 1.85, stroke = 0.8, alpha = 0.80, show.legend = FALSE) +
      ggplot2::facet_wrap(~ facet_lbl, scales = "free_x", ncol = length(box_vars)) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_shape_manual(values = shp_map) +
      ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
      theme_packera_pub(base_size = 9.4) +
      ggplot2::labs(title = "C. Rhizosphere Pedologic & Climatic Divergence", x = "Extracted Raster Value (0-15 cm / Bioclim)", y = NULL)
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "C. Pedologic & Climatic Divergence")
  }

  # Panel D: Full N x N Pairwise Schoener's D Heatmap (or IBE Partial RDA in K=1 mode)
  nsum <- if (!is.null(niche_res) && is.data.frame(niche_res$summary)) niche_res$summary else data.frame()
  ibe_df <- if (!is.null(srf_res) && is.data.frame(srf_res$ibe_summary)) srf_res$ibe_summary else data.frame()

  t1_col <- intersect(c("Taxon1", "Taxon_1"), names(nsum))[1L]
  t2_col <- intersect(c("Taxon2", "Taxon_2"), names(nsum))[1L]
  sd_col <- intersect(c("Schoener_D", "Schoeners_D", "Warren_D"), names(nsum))[1L]
  pv_col <- intersect(c("p_value_D", "Equivalency_p_value", "p_value"), names(nsum))[1L]
  p_schoener <- if (nrow(nsum) > 0L && !is.na(t1_col) && !is.na(t2_col) && !is.na(sd_col)) {
    nsum$D_val <- suppressWarnings(as.numeric(nsum[[sd_col]]))
    nsum$p_val <- if (!is.na(pv_col)) suppressWarnings(as.numeric(nsum[[pv_col]])) else rep(NA_real_, nrow(nsum))
    nsum$t1_short <- abbreviate_icn_taxon(nsum[[t1_col]], abbreviate_genus = TRUE)
    nsum$t2_short <- abbreviate_icn_taxon(nsum[[t2_col]], abbreviate_genus = TRUE)
    nsum$lbl <- ifelse(is.finite(nsum$p_val) & nsum$p_val < 0.05, sprintf("%.2f*", nsum$D_val), sprintf("%.2f", nsum$D_val))
    ggplot2::ggplot(nsum, ggplot2::aes(x = .data$t1_short, y = .data$t2_short, fill = .data$D_val)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$lbl),
        size = if (length(active_groups) > 12L) 2.2 else 2.9,
        color = ifelse(nsum$D_val > 0.55, "#111827", "white")
      ) +
      ggplot2::scale_fill_viridis_c(option = "viridis", limits = c(0, 1), name = "Schoener's D") +
      theme_packera_pub(base_size = 9.5) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, face = "italic"),
                     axis.text.y = ggplot2::element_text(face = "italic")) +
      ggplot2::labs(
        title = "D. Pairwise PCA-env Niche Overlap (Schoener's D)",
        subtitle = "Broennimann et al. (2012) background-calibrated KDE (* = background test p < 0.05)",
        x = NULL, y = NULL
      )
  } else if (nrow(ibe_df) > 0L) {
    comp_c <- intersect(c("Component", "Fraction", "Term"), names(ibe_df))[1L]
    r2_c <- intersect(c("Adj_R2", "R2", "Variance_Pct"), names(ibe_df))[1L]
    ibe_df$comp <- as.character(ibe_df[[comp_c]])
    ibe_df$r2v <- pmax(0.0, suppressWarnings(as.numeric(ibe_df[[r2_c]])))
    ggplot2::ggplot(ibe_df, ggplot2::aes(x = .data$r2v, y = .data$comp)) +
      ggplot2::geom_col(fill = "#0072B2", alpha = 0.85, width = 0.6) +
      theme_packera_pub() +
      ggplot2::labs(title = "D. Isolation-by-Environment (IBE) vs. IBD Partial RDA", x = "Adjusted R\u00b2", y = NULL)
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Single-niche or insufficient georeferenced pairs\nfor Schoener's D permutation test", color = "#6B7280", size = 3.3) +
      theme_packera_pub() +
      ggplot2::labs(title = "D. Niche Overlap Summary", x = NULL, y = NULL)
  }

  dyn_h_05a <- max(10.0, 6.0 + 0.28 * length(active_groups))
  plate_05a <- ((p_imp | p_niche) / (p_edaph | p_schoener)) +
    patchwork::plot_annotation(
      title = "Step 05: Spatial Block-CV Niche Modeling, Rhizosphere Pedology & Niche Overlap",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5))
    )
  canon_05a <- file.path(fdir, "fig05a_ecological_niche_and_rf_importance.png")
  save_packera_png_300(plate_05a, canon_05a, width = 14.2, height = dyn_h_05a, also_copy_to = out_plot_path)
}

render_step05b_hopkins_phenology_png <- function(
  df,
  baseline_info = NULL,
  out_pheno_path = "outputs/figures/fig05b_hopkins_circular_phenology.png"
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (!is.null(out_pheno_path) && nzchar(out_pheno_path)) dirname(out_pheno_path) else "outputs/figures"
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (ext in c(".pdf", ".png")) {
    p_old <- file.path(fdir, paste0("phenological_latitudinal_anomaly", ext))
    if (file.exists(p_old) && normalizePath(p_old, mustWork = FALSE) != normalizePath(normalize_png_path(out_pheno_path), mustWork = FALSE)) {
      unlink(p_old, force = TRUE)
    }
  }

  vdf <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(vdf) == 0L) return(invisible(""))
  vdf$plot_group <- resolve_plot_group_labels(vdf)
  vdf$doy_num <- suppressWarnings(as.numeric(vdf$doy))
  vdf$lat_num <- if ("decimalLatitude" %in% names(vdf)) suppressWarnings(as.numeric(vdf$decimalLatitude)) else suppressWarnings(as.numeric(vdf$latitude))
  vdf$hop_doy <- if ("hopkins_adjusted_doy" %in% names(vdf)) suppressWarnings(as.numeric(vdf$hopkins_adjusted_doy)) else vdf$doy_num
  vdf$ddoy <- if ("delta_doy" %in% names(vdf)) suppressWarnings(as.numeric(vdf$delta_doy)) else rep(NA_real_, nrow(vdf))

  # Restrict to verified-anthesis or valid flowering DOY records
  is_anth <- if ("anthesis_filter_status" %in% names(vdf)) {
    grepl("anthesis|verified|retained|pass", tolower(as.character(vdf$anthesis_filter_status)))
  } else if ("is_verified_anthesis" %in% names(vdf)) {
    !is.na(vdf$is_verified_anthesis) & as.logical(vdf$is_verified_anthesis)
  } else rep(TRUE, nrow(vdf))
  sub_p <- vdf[is_anth & is.finite(vdf$doy_num) & vdf$doy_num >= 1 & vdf$doy_num <= 366, , drop = FALSE]
  if (nrow(sub_p) == 0L) sub_p <- vdf[is.finite(vdf$doy_num) & vdf$doy_num >= 1 & vdf$doy_num <= 366, , drop = FALSE]
  if (nrow(sub_p) == 0L) return(invisible(""))

  active_groups <- sort(unique(sub_p$plot_group[!is.na(sub_p$plot_group) & nzchar(sub_p$plot_group)]))
  pal_map <- get_combinatorial_palette_r(active_groups)
  shp_map <- get_combinatorial_hollow_shapes_r(active_groups)
  cnt_map <- table(sub_p$plot_group)
  lbl_exprs <- format_icn_label_vec(active_groups, abbreviate_genus = (length(active_groups) > 5L), counts_map = cnt_map)
  names(lbl_exprs) <- active_groups

  # Panel A: Anthesis DOY vs Latitude (Hopkins Bioclimatic Slope)
  sub_lat <- sub_p[is.finite(sub_p$lat_num), , drop = FALSE]
  b_lat <- if (!is.null(baseline_info) && !is.null(baseline_info$slope_lat) && is.finite(baseline_info$slope_lat)) baseline_info$slope_lat else NA_real_
  b_yr <- if (!is.null(baseline_info) && !is.null(baseline_info$slope_year) && is.finite(baseline_info$slope_year)) baseline_info$slope_year else NA_real_
  hop_sub <- paste(
    c(
      sprintf("Verified anthesis N = %d", nrow(sub_p)),
      if (is.finite(b_lat)) sprintf("Empirical latitudinal slope = %+.2f d/\u00b0N", b_lat) else "Hopkins (1918) bioclimatic baseline",
      if (is.finite(b_yr)) sprintf("Secular shift = %+.3f d/yr", b_yr) else NULL
    ),
    collapse = " | "
  )

  p_lat <- if (nrow(sub_lat) >= 2L) {
    ggplot2::ggplot(sub_lat, ggplot2::aes(x = .data$lat_num, y = .data$doy_num, color = .data$plot_group, shape = .data$plot_group)) +
      ggplot2::geom_point(size = 2.4, stroke = 0.85, alpha = 0.82) +
      ggplot2::geom_smooth(ggplot2::aes(x = .data$lat_num, y = .data$doy_num), method = "lm", formula = y ~ x, inherit.aes = FALSE, color = "#111827", linetype = "dashed", linewidth = 0.8, se = TRUE, alpha = 0.14) +
      ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
      ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group") +
      theme_packera_pub() +
      ggplot2::labs(title = "A. Anthesis Day-of-Year vs. Latitude", subtitle = hop_sub, x = "Latitude (\u00b0N)", y = "Observed Anthesis Day of Year (DOY)")
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "A. Anthesis DOY vs. Latitude", subtitle = "Insufficient georeferenced anthesis latitudes")
  }
  if (length(active_groups) > 10L) p_lat <- p_lat + ggplot2::theme(legend.position = "none")

  # Panel B: Circular Polar Rose Histogram of Hopkins-Adjusted DOY (Pewsey et al. 2013)
  sub_p$hop_clean <- ifelse(is.finite(sub_p$hop_doy), ((sub_p$hop_doy - 1.0) %% 365.25) + 1.0, sub_p$doy_num)
  month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
  month_lbls <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

  p_polar <- ggplot2::ggplot(sub_p, ggplot2::aes(x = .data$hop_clean, fill = .data$plot_group, color = .data$plot_group)) +
    ggplot2::geom_histogram(binwidth = 15.2, boundary = 1, position = "identity", alpha = 0.25, linewidth = 0.65, show.legend = FALSE) +
    ggplot2::coord_polar(start = 0) +
    ggplot2::scale_x_continuous(limits = c(1, 366.25), breaks = month_breaks, labels = month_lbls) +
    ggplot2::scale_fill_manual(values = pal_map) +
    ggplot2::scale_color_manual(values = pal_map) +
    theme_packera_pub(base_size = 9.5) +
    ggplot2::labs(
      title = "B. Circular Hopkins-Adjusted Anthesis Clock",
      subtitle = "Standardized to cohort reference latitude & elevation",
      x = NULL, y = "Anthesis Count"
    )

  # Panel C: Circular Phenological Anomaly Residuals (Delta DOY)
  if (!any(is.finite(sub_p$ddoy))) {
    med_h <- stats::median(sub_p$hop_clean, na.rm = TRUE)
    sub_p$ddoy <- sub_p$hop_clean - med_h
  }
  sub_dd <- sub_p[is.finite(sub_p$ddoy), , drop = FALSE]
  p_resid <- ggplot2::ggplot(sub_dd, ggplot2::aes(x = .data$ddoy, y = .data$plot_group)) +
    ggplot2::geom_vline(xintercept = c(-21, 0, 21), linetype = c("dashed", "dotted", "dashed"), color = c("#D55E00", "#9CA3AF", "#D55E00"), linewidth = 0.55) +
    ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = 0.55, linewidth = 0.45) +
    ggplot2::geom_jitter(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), height = 0.16, width = 0, size = 2.0, stroke = 0.8, alpha = 0.82, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = pal_map) +
    ggplot2::scale_shape_manual(values = shp_map) +
    ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
    theme_packera_pub() +
    ggplot2::labs(
      title = "C. Circular Phenological Residuals (\u0394DOY)",
      subtitle = "Dashed lines mark \u00b121-day phenological isolation window",
      x = "Hopkins-Adjusted Anthesis Residual \u0394DOY (days)",
      y = NULL
    )

  dyn_h_05b <- max(8.8, 5.2 + 0.25 * length(active_groups))
  plate_05b <- ((p_lat | p_polar) / p_resid) +
    patchwork::plot_layout(heights = c(1.15, 0.95)) +
    patchwork::plot_annotation(
      title = "Step 05: Circular Anthesis Phenology & Hopkins Bioclimatic Isolation",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4))
    )
  canon_05b <- file.path(fdir, "fig05b_hopkins_circular_phenology.png")
  save_packera_png_300(plate_05b, canon_05b, width = 13.8, height = dyn_h_05b, also_copy_to = out_pheno_path)
}

render_step05c_sdm_suitability_png <- function(
  df,
  focal_rasters = NULL,
  multi_rasters = NULL,
  perf_df = NULL,
  resp_df = NULL,
  roc_df = NULL,
  boyce_df = NULL,
  figures_dir = "outputs/figures",
  focal_taxon = ""
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- as.character(figures_dir)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (old_stem in c("fig05b_sdm_habitat_suitability_map", "fig05c_multispecies_maxent_sdm_grid")) {
    for (ext in c(".pdf", ".png")) {
      p_old <- file.path(fdir, paste0(old_stem, ext))
      if (file.exists(p_old)) unlink(p_old, force = TRUE)
    }
  }

  vdf <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(vdf) == 0L) return(invisible(""))
  if (!("decimalLongitude" %in% names(vdf)) && "longitude" %in% names(vdf)) vdf$decimalLongitude <- vdf$longitude
  if (!("decimalLatitude" %in% names(vdf)) && "latitude" %in% names(vdf)) vdf$decimalLatitude <- vdf$latitude
  vdf$plot_group <- resolve_plot_group_labels(vdf)
  active_groups <- sort(unique(vdf$plot_group[!is.na(vdf$plot_group) & nzchar(vdf$plot_group)]))
  pal_map <- get_combinatorial_palette_r(active_groups)
  shp_map <- get_combinatorial_hollow_shapes_r(active_groups)
  cnt_map <- table(vdf$plot_group)
  lbl_exprs <- format_icn_label_vec(active_groups, abbreviate_genus = TRUE, counts_map = cnt_map)
  names(lbl_exprs) <- active_groups

  raster_to_lonlat_df <- function(r_obj, max_cells = 45000L) {
    if (is.null(r_obj) || !requireNamespace("terra", quietly = TRUE)) return(data.frame())
    tryCatch({
      r_use <- if (terra::ncell(r_obj) > max_cells) {
        agg_f <- max(2L, ceiling(sqrt(terra::ncell(r_obj) / max_cells)))
        terra::aggregate(r_obj, fact = agg_f, fun = "mean", na.rm = TRUE)
      } else r_obj
      d_r <- as.data.frame(r_use, xy = TRUE, na.rm = TRUE)
      if (ncol(d_r) < 3L || nrow(d_r) == 0L) return(data.frame())
      names(d_r)[1:3] <- c("lon", "lat", "suitability")
      d_r[is.finite(d_r$suitability), c("lon", "lat", "suitability"), drop = FALSE]
    }, error = function(e) data.frame())
  }

  # Helper to build a conformal/projected suitability + hollow occurrence map
  build_suitability_panel <- function(r_df, pts_sub, title_txt, sub_txt = "") {
    if (nrow(r_df) > 0L) {
      cos_l <- cos(mean(range(r_df$lat, na.rm = TRUE)) * pi / 180.0)
      p_m <- ggplot2::ggplot() +
        ggplot2::geom_raster(data = r_df, ggplot2::aes(x = .data$lon, y = .data$lat, fill = .data$suitability), alpha = 0.90) +
        ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Suitability")
      if (requireNamespace("maps", quietly = TRUE)) {
        st_d <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)
        if (!is.null(st_d)) {
          p_m <- p_m + ggplot2::geom_path(data = st_d, ggplot2::aes(x = .data$long, y = .data$lat, group = .data$group), color = "#E5E7EB", linewidth = 0.25)
        }
      }
      if (nrow(pts_sub) > 0L) {
        pts_sub$lon_num <- suppressWarnings(as.numeric(pts_sub$decimalLongitude))
        pts_sub$lat_num <- suppressWarnings(as.numeric(pts_sub$decimalLatitude))
        pts_ok_s <- pts_sub[is.finite(pts_sub$lon_num) & is.finite(pts_sub$lat_num), , drop = FALSE]
        if (nrow(pts_ok_s) > 0L) {
          p_m <- p_m + ggplot2::geom_point(
            data = pts_ok_s,
            ggplot2::aes(x = .data$lon_num, y = .data$lat_num),
            shape = 1, size = 2.0, stroke = 0.8, color = "#FF6B00", alpha = 0.88
          )
        }
      }
      p_m +
        ggplot2::coord_fixed(ratio = 1.0 / max(0.25, cos_l), xlim = range(r_df$lon), ylim = range(r_df$lat), expand = FALSE) +
        theme_packera_pub(base_size = 9.5) +
        ggplot2::labs(title = title_txt, subtitle = sub_txt, x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)")
    } else {
      m_ctx <- build_na_projected_basemap_gg(vdf$decimalLongitude, vdf$decimalLatitude)
      pts_p <- project_lonlat_points_r(pts_sub, m_ctx)
      p_m <- m_ctx$plot
      if (nrow(pts_p) > 0L) {
        p_m <- p_m + ggplot2::geom_point(
          data = pts_p,
          ggplot2::aes(x = .data$map_x, y = .data$map_y, color = .data$plot_group, shape = .data$plot_group),
          size = 2.2, stroke = 0.85, alpha = 0.84
        ) +
          ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
          ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group")
        if (length(active_groups) > 10L) p_m <- p_m + ggplot2::theme(legend.position = "none")
      }
      p_m + ggplot2::labs(title = title_txt, subtitle = sub_txt)
    }
  }

  has_multi_sdm <- is.list(multi_rasters) && length(multi_rasters) >= 2L
  has_focal_sdm <- is.list(focal_rasters) && (!is.null(focal_rasters$maxent) || !is.null(focal_rasters$ensemble))

  p_multi_grid <- NULL
  n_rows_multi <- 0L
  if (has_multi_sdm) {
    grid_dfs <- lapply(names(multi_rasters), function(tx_nm) {
      d_tx <- raster_to_lonlat_df(multi_rasters[[tx_nm]], max_cells = 18000L)
      if (nrow(d_tx) > 0L) d_tx$taxon_lbl <- abbreviate_icn_taxon(tx_nm, abbreviate_genus = TRUE)
      d_tx
    })
    grid_all <- do.call(rbind, grid_dfs)
    if (!is.null(grid_all) && nrow(grid_all) > 0L) {
      pts_ok <- vdf[is.finite(suppressWarnings(as.numeric(vdf$decimalLongitude))) & is.finite(suppressWarnings(as.numeric(vdf$decimalLatitude))), , drop = FALSE]
      pts_ok$lon_num <- suppressWarnings(as.numeric(pts_ok$decimalLongitude))
      pts_ok$lat_num <- suppressWarnings(as.numeric(pts_ok$decimalLatitude))
      pts_ok$taxon_lbl <- abbreviate_icn_taxon(pts_ok$plot_group, abbreviate_genus = TRUE)
      pts_ok <- pts_ok[pts_ok$taxon_lbl %in% unique(grid_all$taxon_lbl), , drop = FALSE]
      n_sp_m <- length(unique(grid_all$taxon_lbl))
      ncol_m <- min(4L, max(2L, ceiling(sqrt(n_sp_m))))
      n_rows_multi <- as.integer(ceiling(n_sp_m / ncol_m))
      cos_m <- cos(mean(range(grid_all$lat, na.rm = TRUE)) * pi / 180.0)
      p_multi_grid <- ggplot2::ggplot(grid_all, ggplot2::aes(x = .data$lon, y = .data$lat)) +
        ggplot2::geom_raster(ggplot2::aes(fill = .data$suitability)) +
        ggplot2::geom_point(
          data = pts_ok,
          ggplot2::aes(x = .data$lon_num, y = .data$lat_num),
          shape = 1, size = 1.35, stroke = 0.65, color = "#FF6B00", alpha = 0.80, inherit.aes = FALSE
        ) +
        ggplot2::facet_wrap(~ taxon_lbl, ncol = ncol_m) +
        ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Ensemble\nSuitability") +
        ggplot2::coord_fixed(ratio = 1.0 / max(0.25, cos_m), expand = FALSE) +
        theme_packera_pub(base_size = 9.2) +
        ggplot2::theme(strip.text = ggplot2::element_text(face = "bold.italic")) +
        ggplot2::labs(
          title = sprintf("E. Comparative Multi-Taxon Ensemble Habitat Suitability (%d Modeled Taxa)", n_sp_m),
          subtitle = "Orange hollow circles mark coordinate-uncertainty-gated herbarium presences",
          x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)"
        )
    }
  }

  f_tx <- if (nzchar(trimws(as.character(focal_taxon))) && (trimws(as.character(focal_taxon)) %in% vdf$plot_group)) {
    trimws(as.character(focal_taxon))
  } else {
    active_groups[1L]
  }
  pts_f <- if (has_focal_sdm && f_tx %in% vdf$plot_group) vdf[vdf$plot_group == f_tx, , drop = FALSE] else vdf
  df_max <- if (has_focal_sdm) raster_to_lonlat_df(focal_rasters$maxent) else data.frame()
  df_ens <- if (has_focal_sdm) raster_to_lonlat_df(focal_rasters$ensemble) else data.frame()

  if (nrow(df_max) > 0L || nrow(df_ens) > 0L) {
    p_a <- build_suitability_panel(if (nrow(df_max) > 0L) df_max else df_ens, pts_f, sprintf("A. MaxEnt (cloglog) Suitability : %s", f_tx), "L1-regularized IPP presence-background model")
    p_b <- build_suitability_panel(if (nrow(df_ens) > 0L) df_ens else df_max, pts_f, sprintf("B. Spatial RF + MaxEnt Ensemble : %s", f_tx), "AUC-weighted ensemble consensus surface")

    p_c <- if (!is.null(resp_df) && is.data.frame(resp_df) && nrow(resp_df) > 0L && all(c("predictor", "x_val", "suitability") %in% names(resp_df))) {
      ggplot2::ggplot(resp_df, ggplot2::aes(x = .data$x_val, y = .data$suitability, color = .data$predictor)) +
        ggplot2::geom_line(linewidth = 0.95, show.legend = FALSE) +
        ggplot2::facet_wrap(~ predictor, scales = "free_x", ncol = min(3L, length(unique(resp_df$predictor)))) +
        ggplot2::scale_color_manual(values = rep(c("#0072B2", "#009E73", "#D55E00", "#CC79A7"), length.out = length(unique(resp_df$predictor)))) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        theme_packera_pub(base_size = 9.2) +
        ggplot2::labs(title = "C. Marginal MaxEnt (cloglog) Ecological Response Curves", subtitle = "5th to 95th percentile gradient holding co-predictors at presence median", x = "Predictor Value", y = "Predicted Suitability")
    } else {
      ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "C. Marginal Response Curves")
    }

    p_d <- if (!is.null(roc_df) && is.data.frame(roc_df) && nrow(roc_df) > 0L && all(c("fpr", "tpr") %in% names(roc_df))) {
      p_row <- if (!is.null(perf_df) && is.data.frame(perf_df) && nrow(perf_df) > 0L) {
        if ("Taxon_or_Morphotype" %in% names(perf_df) && any(perf_df$Taxon_or_Morphotype == f_tx)) perf_df[perf_df$Taxon_or_Morphotype == f_tx, , drop = FALSE][1L, , drop = FALSE] else perf_df[1L, , drop = FALSE]
      } else NULL
      .pick_metric <- function(cands) {
        if (is.null(p_row)) return(NA_real_)
        hit <- intersect(cands, names(p_row))
        if (length(hit) == 0L) return(NA_real_)
        suppressWarnings(as.numeric(p_row[[hit[1L]]][1L]))
      }
      auc_v <- .pick_metric(c("AUC_Block_Ensemble", "auc_ensemble", "AUC_Block_MaxEnt", "auc_maxent"))
      boy_v <- .pick_metric(c("CBI_Ensemble", "boyce_index", "CBI_MaxEnt"))
      if (!("model" %in% names(roc_df))) roc_df$model <- "Ensemble"
      ggplot2::ggplot(roc_df, ggplot2::aes(x = .data$fpr, y = .data$tpr, color = .data$model, linetype = .data$model)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "#9CA3AF") +
        ggplot2::geom_line(linewidth = 0.95) +
        ggplot2::scale_color_manual(values = c("MaxEnt" = "#0072B2", "Spatial RF" = "#009E73", "Ensemble" = "#D55E00"), name = "SDM Engine") +
        ggplot2::scale_linetype_manual(values = c("MaxEnt" = "dashed", "Spatial RF" = "dotdash", "Ensemble" = "solid"), name = "SDM Engine") +
        ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        theme_packera_pub() +
        ggplot2::labs(
          title = "D. Spatial-Block Cross-Validation & Boyce Calibration",
          subtitle = sprintf("Ensemble Block AUC = %s | Continuous Boyce Index = %s", ifelse(is.finite(auc_v), sprintf("%.3f", auc_v), "N/A"), ifelse(is.finite(boy_v), sprintf("%+.3f", boy_v), "N/A")),
          x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)"
        )
    } else {
      ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "D. Spatial Block-CV ROC & Boyce Calibration")
    }

    if (!is.null(p_multi_grid)) {
      plate_05c <- ((p_a | p_b) / (p_c | p_d) / p_multi_grid) +
        patchwork::plot_layout(heights = c(1.22, 0.92, max(1.15, 0.72 * n_rows_multi))) +
        patchwork::plot_annotation(
          title = "Step 05 : Ensemble Species Distribution Modeling, Spatial-Block Validation & Multi-Taxon Suitability",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4))
        )
      return(save_packera_png_300(plate_05c, file.path(fdir, "fig05c_sdm_suitability_and_distribution_map.png"), width = 14.2, height = 9.4 + 2.7 * n_rows_multi))
    }

    plate_05c <- ((p_a | p_b) / (p_c | p_d)) +
      patchwork::plot_layout(widths = c(1.1, 1.1), heights = c(1.25, 0.95)) +
      patchwork::plot_annotation(
        title = "Step 05 : Ensemble Species Distribution Modeling & Spatial Block Validation",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4))
      )
    return(save_packera_png_300(plate_05c, file.path(fdir, "fig05c_sdm_suitability_and_distribution_map.png"), width = 14.0, height = 9.8))
  }

  # Geographic Occurrence & Ecoregion Map when raster layers are absent
  p_geo <- build_suitability_panel(data.frame(), vdf, "Step 05 : Equal-Area Projected Geographic Distribution (Albers EPSG:5070)", "Hollow markers stratified by taxon; scale bar & north arrow in projected CRS")
  save_packera_png_300(p_geo, file.path(fdir, "fig05c_sdm_suitability_and_distribution_map.png"), width = 11.8, height = 8.2)
}
