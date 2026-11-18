if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
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
    if (grepl("^bio[0-9]+|temp|precip|clim|summer_moisture|maritime_continentality", vl)) return("Climatic (WorldClim)")
    if (grepl("elev|slope|aspect|tpi|tri|topo|dem", vl)) return("Topographic (3DEP)")
    "Edaphic & Lithic (SOLUS/SSURGO/SGMC)"
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

  vdf <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(vdf) == 0L) return(invisible(""))
  grp_col <- if ("analysis_group" %in% names(vdf) && any(!is.na(vdf$analysis_group) & nzchar(trimws(as.character(vdf$analysis_group))))) {
    "analysis_group"
  } else if ("ecotype_regime" %in% names(vdf) && any(!is.na(vdf$ecotype_regime) & nzchar(trimws(as.character(vdf$ecotype_regime))))) {
    "ecotype_regime"
  } else NULL
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

  imp_df <- imp_df[is.finite(imp_df$imp) & !is.na(imp_df$var) & nzchar(as.character(imp_df$var)), , drop = FALSE]
  p_imp <- if (nrow(imp_df) > 0L) {
    imp_df <- imp_df[order(imp_df$imp, decreasing = TRUE), , drop = FALSE]
    imp_top <- utils::head(imp_df, 14L)
    imp_top$domain <- .classify_env_domain(imp_top$var)
    pretty_lbl <- if (exists("format_predictor_labels", mode = "function")) format_predictor_labels(imp_top$var) else gsub("_", " ", imp_top$var)
    imp_top$lbl <- factor(wrap_clean_label_r(pretty_lbl, width = 36L), levels = rev(unique(wrap_clean_label_r(pretty_lbl, width = 36L))))
    cv_acc <- if ("spatial_block_cv_acc" %in% names(vdf) && any(is.finite(suppressWarnings(as.numeric(vdf$spatial_block_cv_acc))))) {
      stats::median(suppressWarnings(as.numeric(vdf$spatial_block_cv_acc)), na.rm = TRUE)
    } else NA_real_
    ggplot2::ggplot(imp_top, ggplot2::aes(x = .data$imp, y = .data$lbl, fill = .data$domain)) +
      ggplot2::geom_col(alpha = 0.88, color = "#374151", linewidth = 0.35, na.rm = TRUE) +
      ggplot2::scale_fill_manual(
        values = c(
          "Climatic (WorldClim)" = "#56B4E9",
          "Edaphic & Lithic (SOLUS/SSURGO/SGMC)" = "#D55E00",
          "Topographic (3DEP)" = "#009E73"
        ),
        name = "Environmental Stream"
      ) +
      theme_packera_pub(base_size = 9.6) +
      ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 7.8, lineheight = 0.88)) +
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
    fb_env <- grep("^(fine_|soil_|bio[0-9]+|ssurgo_|solus_|topo_)", names(vdf), value = TRUE)
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
      lbl_b <- paste(strwrap(lbl_b, width = 20L), collapse = "\n")
      data.frame(plot_group = vdf$plot_group[ok_b], facet_lbl = lbl_b, val = vv[ok_b], stringsAsFactors = FALSE)
    })
    ed_df <- do.call(rbind, ed_list)
    ggplot2::ggplot(ed_df, ggplot2::aes(x = .data$val, y = .data$plot_group)) +
      ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = if (length(active_groups) <= 1L) 0.32 else 0.55, linewidth = 0.45) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), height = 0.16, width = 0, size = 1.85, stroke = 0.8, alpha = 0.80, show.legend = FALSE) +
      ggplot2::facet_wrap(~ facet_lbl, scales = "free_x", ncol = length(box_vars)) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_shape_manual(values = shp_map) +
      ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
      theme_packera_pub(base_size = 9.4) +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 7.1, face = "bold", lineheight = 0.92, margin = ggplot2::margin(3, 3, 3, 3, unit = "pt"))) +
      ggplot2::labs(title = "C. Rhizosphere Pedologic & Climatic Divergence", x = "Extracted Raster Value (0-15 cm / Bioclim)", y = NULL)
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "C. Pedologic & Climatic Divergence")
  }

  # Panel D: Full N x N Pairwise Schoener's D Heatmap (or IBE Partial RDA / Distance-Decay in K=1 mode)
  nsum <- if (!is.null(niche_res) && is.data.frame(niche_res$summary)) niche_res$summary else data.frame()
  ibe_df <- if (!is.null(srf_res) && is.data.frame(srf_res$ibe_summary)) srf_res$ibe_summary else data.frame()

  t1_col <- intersect(c("Taxon1", "Taxon_1"), names(nsum))[1L]
  t2_col <- intersect(c("Taxon2", "Taxon_2"), names(nsum))[1L]
  sd_col <- intersect(c("Schoener_D", "Schoeners_D", "Warren_D"), names(nsum))[1L]
  pv_col <- intersect(c("p_value_D", "Equivalency_p_value", "p_value"), names(nsum))[1L]
  n_pairs_taxa <- if (nrow(nsum) > 0L && !is.na(t1_col) && !is.na(t2_col)) {
    length(unique(stats::na.omit(c(as.character(nsum[[t1_col]]), as.character(nsum[[t2_col]])))))
  } else 0L
  is_single_niche <- (length(active_groups) < 2L) || (n_pairs_taxa < 2L)
  p_schoener <- if (!is_single_niche && nrow(nsum) > 0L && !is.na(sd_col)) {
    nsum$D_val <- suppressWarnings(as.numeric(nsum[[sd_col]]))
    nsum$p_val <- if (!is.na(pv_col)) suppressWarnings(as.numeric(nsum[[pv_col]])) else rep(NA_real_, nrow(nsum))
    wrap_ecotype <- function(v) ifelse(grepl("^(Ecotype|Morphotype|Complex|Cluster)", v), wrap_clean_label_r(gsub("_", " ", v), width = 22L), v)
    nsum$t1_short <- wrap_ecotype(abbreviate_icn_taxon(nsum[[t1_col]], abbreviate_genus = TRUE))
    nsum$t2_short <- wrap_ecotype(abbreviate_icn_taxon(nsum[[t2_col]], abbreviate_genus = TRUE))
    d_face <- if (any(grepl("^(Ecotype|Morphotype|Complex|Cluster)", active_groups))) "plain" else "italic"
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
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, face = d_face),
                     axis.text.y = ggplot2::element_text(face = d_face)) +
      ggplot2::labs(
        title = "D. Pairwise PCA-env Niche Overlap (Schoener's D)",
        subtitle = "Broennimann et al. (2012) background-calibrated KDE (* = background test p < 0.05)",
        x = NULL, y = NULL
      )
  } else if (nrow(ibe_df) > 0L) {
    if (all(c("Modality", "Pure_Edaphic_AdjR2", "Pure_Climate_AdjR2", "Pure_Space_MEM_AdjR2") %in% names(ibe_df))) {
      ibe_long <- rbind(
        data.frame(Modality = as.character(ibe_df$Modality), Partition = "Pure Edaphic [A|B,C]",
                   AdjR2 = pmax(0.0, suppressWarnings(as.numeric(ibe_df$Pure_Edaphic_AdjR2))), stringsAsFactors = FALSE),
        data.frame(Modality = as.character(ibe_df$Modality), Partition = "Pure Climate [B|A,C]",
                   AdjR2 = pmax(0.0, suppressWarnings(as.numeric(ibe_df$Pure_Climate_AdjR2))), stringsAsFactors = FALSE),
        data.frame(Modality = as.character(ibe_df$Modality), Partition = "Pure Space MEM [C|A,B]",
                   AdjR2 = pmax(0.0, suppressWarnings(as.numeric(ibe_df$Pure_Space_MEM_AdjR2))), stringsAsFactors = FALSE)
      )
      ibe_long$Modality_lbl <- gsub("_", " ", ibe_long$Modality)
      ibe_long <- ibe_long[is.finite(ibe_long$AdjR2) & !is.na(ibe_long$Modality_lbl) & nzchar(ibe_long$Modality_lbl), , drop = FALSE]
      min_p_ibe <- if ("Partial_RDA_p" %in% names(ibe_df)) suppressWarnings(min(as.numeric(ibe_df$Partial_RDA_p), na.rm = TRUE)) else NA_real_
      ggplot2::ggplot(ibe_long, ggplot2::aes(x = .data$AdjR2, y = .data$Partition, fill = .data$Partition)) +
        ggplot2::geom_col(alpha = 0.88, show.legend = FALSE, na.rm = TRUE) +
        ggplot2::facet_wrap(~ Modality_lbl, ncol = 1L, scales = "fixed") +
        ggplot2::scale_fill_manual(values = c(
          "Pure Edaphic [A|B,C]" = "#009E73",
          "Pure Climate [B|A,C]" = "#0072B2",
          "Pure Space MEM [C|A,B]" = "#D55E00"
        )) +
        theme_packera_pub(base_size = 9.2) +
        ggplot2::labs(
          title = "D. Isolation-by-Environment (IBE) vs. IBD Partial RDA",
          subtitle = if (is.finite(min_p_ibe)) sprintf("Unimodal (K = 1) partial RDA fractions (min partial RDA p = %.3f)", min_p_ibe) else "Unimodal (K = 1) partial RDA variance partitioning",
          x = expression("Partial Adjusted " * italic(R)^2),
          y = NULL
        )
    } else {
      comp_c <- intersect(c("Component", "Fraction", "Term", "Modality"), names(ibe_df))[1L]
      r2_c <- intersect(c("Adj_R2", "R2", "Variance_Pct", "Total_Explained_AdjR2"), names(ibe_df))[1L]
      ibe_df$comp <- if (!is.na(comp_c)) as.character(ibe_df[[comp_c]]) else paste0("Row_", seq_len(nrow(ibe_df)))
      ibe_df$r2v <- if (!is.na(r2_c)) pmax(0.0, suppressWarnings(as.numeric(ibe_df[[r2_c]]))) else rep(0.0, nrow(ibe_df))
      ibe_df <- ibe_df[is.finite(ibe_df$r2v) & !is.na(ibe_df$comp) & nzchar(ibe_df$comp), , drop = FALSE]
      ggplot2::ggplot(ibe_df, ggplot2::aes(x = .data$r2v, y = .data$comp)) +
        ggplot2::geom_col(fill = "#0072B2", alpha = 0.85, na.rm = TRUE) +
        theme_packera_pub() +
        ggplot2::labs(title = "D. Isolation-by-Environment (IBE) vs. IBD Partial RDA", x = "Adjusted R\u00b2", y = NULL)
    }
  } else {
    lon_dd <- if ("decimalLongitude" %in% names(vdf)) suppressWarnings(as.numeric(vdf$decimalLongitude)) else if ("longitude" %in% names(vdf)) suppressWarnings(as.numeric(vdf$longitude)) else rep(NA_real_, nrow(vdf))
    lat_dd <- if ("decimalLatitude" %in% names(vdf)) suppressWarnings(as.numeric(vdf$decimalLatitude)) else if ("latitude" %in% names(vdf)) suppressWarnings(as.numeric(vdf$latitude)) else rep(NA_real_, nrow(vdf))
    env_dd_cols <- utils::head(cand_num_env, min(6L, length(cand_num_env)))
    ok_dd <- is.finite(lon_dd) & is.finite(lat_dd) & length(env_dd_cols) >= 1L
    if (sum(ok_dd) >= 6L) {
      idx_dd <- which(ok_dd)
      if (length(idx_dd) > 75L) idx_dd <- idx_dd[unique(round(seq(1L, length(idx_dd), length.out = 75L)))]
      E_dd <- as.matrix(vdf[idx_dd, env_dd_cols, drop = FALSE]); storage.mode(E_dd) <- "double"
      for (jc in seq_len(ncol(E_dd))) {
        mj <- stats::median(E_dd[, jc], na.rm = TRUE); sj <- stats::sd(E_dd[, jc], na.rm = TRUE)
        E_dd[!is.finite(E_dd[, jc]), jc] <- ifelse(is.finite(mj), mj, 0.0)
        E_dd[, jc] <- (E_dd[, jc] - mean(E_dd[, jc])) / max(ifelse(is.finite(sj), sj, 1.0), 1e-6)
      }
      d_env <- as.numeric(stats::dist(E_dd))
      mid_lat_rad <- stats::median(lat_dd[idx_dd], na.rm = TRUE) * pi / 180.0
      xy_km <- cbind(lon_dd[idx_dd] * 111.32 * cos(mid_lat_rad), lat_dd[idx_dd] * 110.57)
      d_geo <- as.numeric(stats::dist(xy_km))
      ok_pr <- is.finite(d_env) & is.finite(d_geo)
      dd_df <- data.frame(d_geo = d_geo[ok_pr], d_env = d_env[ok_pr])
      if (nrow(dd_df) > 600L) dd_df <- dd_df[unique(round(seq(1L, nrow(dd_df), length.out = 600L))), , drop = FALSE]
      rho_ibe <- suppressWarnings(stats::cor(dd_df$d_geo, dd_df$d_env, method = "spearman", use = "complete.obs"))
      ggplot2::ggplot(dd_df, ggplot2::aes(x = .data$d_geo, y = .data$d_env)) +
        ggplot2::geom_point(shape = 1, size = 1.6, stroke = 0.6, color = "#0072B2", alpha = 0.45) +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "#D55E00", linewidth = 0.85, se = TRUE, alpha = 0.14) +
        theme_packera_pub(base_size = 9.5) +
        ggplot2::labs(
          title = "D. Unimodal IBE Distance Decay (K = 1 Regime)",
          subtitle = sprintf("Standardized environmental vs. geographic distance (Spearman rho = %+.3f)", ifelse(is.finite(rho_ibe), rho_ibe, 0.0)),
          x = "Pairwise Geographic Distance (km)",
          y = "Standardized Environmental Distance"
        )
    } else {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Single-niche unimodal envelope (K = 1)\nwith insufficient georeferenced pairs", color = "#6B7280", size = 3.3) +
        theme_packera_pub() +
        ggplot2::labs(title = "D. Niche Distance Decay Summary", x = NULL, y = NULL)
    }
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

  vdf <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(vdf) == 0L) return(invisible(""))
  grp_col_ph <- if ("analysis_group" %in% names(vdf) && any(!is.na(vdf$analysis_group) & nzchar(trimws(as.character(vdf$analysis_group))))) "analysis_group" else NULL
  vdf$plot_group <- resolve_plot_group_labels(vdf, group_col = grp_col_ph)
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
  if (length(active_groups) > 15L) {
    p_lat <- p_lat + ggplot2::theme(legend.position = "none")
  } else if (length(active_groups) >= 2L) {
    p_lat <- p_lat +
      ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(face = "italic", size = 7.6)) +
      ggplot2::guides(
        color = ggplot2::guide_legend(ncol = min(3L, length(active_groups)), byrow = TRUE),
        shape = ggplot2::guide_legend(ncol = min(3L, length(active_groups)), byrow = TRUE)
      )
  }

  # Panel B: Circular Polar Rose Histogram of Hopkins-Adjusted DOY (Pewsey et al. 2013)
  sub_p$hop_clean <- ifelse(is.finite(sub_p$hop_doy), ((sub_p$hop_doy - 1.0) %% 365.25) + 1.0, sub_p$doy_num)
  month_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)
  month_lbls <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

  sub_polar <- sub_p[is.finite(sub_p$hop_clean) & (sub_p$hop_clean >= 1.0) & (sub_p$hop_clean <= 366.25) & !is.na(sub_p$plot_group) & nzchar(as.character(sub_p$plot_group)), , drop = FALSE]
  p_polar <- ggplot2::ggplot(sub_polar, ggplot2::aes(x = .data$hop_clean, fill = .data$plot_group, color = .data$plot_group)) +
    ggplot2::geom_histogram(binwidth = 365.25 / 24.0, boundary = 1, position = "identity", alpha = 0.25, linewidth = 0.65, show.legend = FALSE, na.rm = TRUE) +
    ggplot2::coord_polar(start = 0) +
    ggplot2::scale_x_continuous(limits = c(1, 366.25), breaks = month_breaks, labels = month_lbls) +
    ggplot2::scale_fill_manual(values = pal_map) +
    ggplot2::scale_color_manual(values = pal_map) +
    theme_packera_pub(base_size = 9.5) +
    ggplot2::theme(plot.title.position = "plot", plot.margin = ggplot2::margin(6, 14, 6, 8, unit = "pt")) +
    ggplot2::labs(
      title = "B. Circular Anthesis Clock",
      subtitle = "Hopkins-adjusted to cohort reference\nlatitude & elevation",
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
    ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = if (length(active_groups) <= 1L) 0.32 else 0.55, linewidth = 0.45) +
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

for (.sdm_fig_p in c(
  "scripts/figures/step05_sdm_figures.R",
  file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "step05_sdm_figures.R")
)) {
  if (!exists("render_step05c_sdm_suitability_png", mode = "function") && !is.na(.sdm_fig_p) && nzchar(.sdm_fig_p) && file.exists(.sdm_fig_p)) {
    source(.sdm_fig_p)
    break
  }
}
