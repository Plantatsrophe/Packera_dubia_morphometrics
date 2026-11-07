# ==============================================================================
# Step 06 Module: Ecological Niche Figures, Spatial RF Plots, & MaxEnt/SRF SDMs
# Extracted from 06_multimodal_spatial_rf.R
# ==============================================================================

export_dynamic_fig06_ecological_niche_phenology <- function(df, srf_res = NULL, niche_res = NULL, out_plot_pdf = NULL) {
  fdir <- if (exists("resolve_active_figures_dir", mode = "function")) {
    resolve_active_figures_dir(reference_path = out_plot_pdf)
  } else if (!is.null(out_plot_pdf) && nzchar(out_plot_pdf)) {
    dirname(out_plot_pdf)
  } else {
    "outputs/figures"
  }
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

  tx_clean <- if (exists("resolve_plot_group_column", mode = "function")) {
    resolve_plot_group_column(df)
  } else if ("analysis_group" %in% names(df) && length(unique(stats::na.omit(df$analysis_group))) >= 2L) {
    as.character(df$analysis_group)
  } else {
    as.character(df$species_standardized)
  }
  present_tx <- sort(unique(tx_clean[!is.na(tx_clean) & nzchar(tx_clean) & tx_clean != "Unknown"]))
  if (length(present_tx) == 0L) present_tx <- get_active_packera_taxa(df$species_standardized)
  if (length(present_tx) == 0L) return(invisible(NULL))

  dyn_pal <- resolve_dynamic_group_palette(present_tx)
  hollow_pchs <- c(1L, 2L, 5L, 0L, 6L, 4L, 3L, 8L)
  dyn_pch <- if (exists("get_dynamic_packera_shapes", mode = "function")) {
    sh <- get_dynamic_packera_shapes(present_tx)
    if (any(is.na(sh[present_tx]))) stats::setNames(rep_len(hollow_pchs, length(present_tx)), present_tx) else sh
  } else {
    stats::setNames(rep_len(hollow_pchs, length(present_tx)), present_tx)
  }
  alpha_val <- if (nrow(df) > 1000L) 0.45 else 0.78
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = alpha_val)
  pch_tx  <- dyn_pch[tx_clean]; pch_tx[is.na(pch_tx)] <- 1L

  if (!("fine_soil_sand" %in% names(df)) && "soil_sand" %in% names(df)) df$fine_soil_sand <- df$soil_sand
  if (!("fine_soil_ph" %in% names(df)) && "soil_ph" %in% names(df)) df$fine_soil_ph <- df$soil_ph

  imp_df <- if (!is.null(srf_res) && !is.null(srf_res$importance)) srf_res$importance else NULL
  top_sel <- select_top_rf_predictors(imp_df, df)
  top_env_1 <- top_sel$top_env_1
  top_env_2 <- top_sel$top_env_2
  top_edaph_1 <- top_sel$top_edaph_1
  top_env_lbls <- format_step06_predictor_labels(c(top_env_1, top_env_2))
  top_edaph_lbl <- format_step06_predictor_labels(top_edaph_1)[1]

  render_plate <- function() {
    graphics::par(mfrow = c(2, 2), mar = c(4.5, 11.0, 3.0, 1.2), oma = c(0.5, 0.5, 1.8, 0.5))

    # Panel A: Spatial RF Predictor Importance Ranking
    pred_raw <- c(top_env_2, top_env_1)
    imp_vals <- c(0.15, 0.25)
    if (!is.null(imp_df) && nrow(imp_df) > 0L && all(c("Predictor", "MeanDecreaseAccuracy") %in% names(imp_df))) {
      ord <- utils::tail(order(as.numeric(imp_df$MeanDecreaseAccuracy), decreasing = FALSE), 12)
      pred_raw <- as.character(imp_df$Predictor[ord])
      imp_vals <- pmax(0, as.numeric(imp_df$MeanDecreaseAccuracy[ord]))
    }
    pred_lbls <- format_step06_predictor_labels(pred_raw)
    bar_cols <- ifelse(
      grepl("^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bdricm|resdept_r|cfvo|awc|alr_|sgmc_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)", pred_raw),
      "#009E73",
      ifelse(grepl("^MEM", pred_raw), "#CC79A7", "#0072B2")
    )
    graphics::barplot(
      imp_vals, names.arg = pred_lbls, horiz = TRUE, las = 1,
      col = bar_cols, border = "white", cex.names = 0.70,
      xlim = c(0, max(imp_vals, 0.05, na.rm = TRUE) * 1.18),
      xlab = "Permutation Importance (Mean Decrease Accuracy)",
      main = "A. Spatial RF Predictor Importance (10m-30m Stack)", font.main = 2
    )
    graphics::grid(nx = NULL, ny = NA, col = "grey88")
    graphics::legend(
      "bottomright",
      legend = c("10m-30m Edaphic/Lithology", "WorldClim Macro-Climate", "Spatial MEMs"),
      fill = c("#009E73", "#0072B2", "#CC79A7"), border = NA, cex = 0.68, bty = "n"
    )

    # Panel B: Dynamic Top-2 RF Environmental Predictors (Hollow Markers + 85% Ellipses)
    graphics::par(mar = c(4.5, 4.8, 3.0, 1.2))
    xv <- suppressWarnings(as.numeric(df[[top_env_1]]))
    yv <- suppressWarnings(as.numeric(df[[top_env_2]]))
    xv[!is.finite(xv)] <- stats::median(xv[is.finite(xv)], na.rm = TRUE)
    yv[!is.finite(yv)] <- stats::median(yv[is.finite(yv)], na.rm = TRUE)
    xv[!is.finite(xv)] <- 0.0; yv[!is.finite(yv)] <- 0.0

    graphics::plot(
      xv, yv, type = "n",
      xlab = top_env_lbls[1], ylab = top_env_lbls[2],
      main = "B. Top Spatial RF Niche Predictors (85% Confidence Ellipses)", font.main = 2
    )
    graphics::grid(col = "grey88")
    graphics::points(xv, yv, col = cols_tx, pch = pch_tx, bg = NA, cex = 0.85, lwd = 1.0)
    if ("is_type_specimen" %in% names(df)) {
      m_t_b <- which(as.logical(df$is_type_specimen) == TRUE & is.finite(xv) & is.finite(yv))
      if (length(m_t_b) > 0L) {
        graphics::points(xv[m_t_b], yv[m_t_b], pch = 5L, bg = NA, col = cols_tx[m_t_b], cex = 1.45, lwd = 2.0)
        graphics::points(xv[m_t_b], yv[m_t_b], pch = 1L, bg = NA, col = "#111111", cex = 1.85, lwd = 1.4)
      }
    }
    t_seq <- seq(0, 2 * pi, length.out = 100)
    chi_r <- sqrt(stats::qchisq(0.85, df = 2))
    for (tx in present_tx) {
      m_g <- which(tx_clean == tx & is.finite(xv) & is.finite(yv))
      if (length(m_g) >= 3L) {
        cov_g <- tryCatch(stats::cov(cbind(xv[m_g], yv[m_g])), error = function(e) NULL)
        if (!is.null(cov_g) && all(is.finite(cov_g)) && det(cov_g) > 1e-12) {
          eg <- eigen(cov_g, symmetric = TRUE)
          pts_el <- t(c(mean(xv[m_g]), mean(yv[m_g])) + chi_r * eg$vectors %*% diag(sqrt(pmax(eg$values, 0))) %*% rbind(cos(t_seq), sin(t_seq)))
          graphics::lines(pts_el[, 1], pts_el[, 2], col = dyn_pal[[tx]], lwd = 1.5)
        }
      }
    }
    short_tx <- gsub("^Packera\\s+", "P. ", present_tx)
    graphics::legend(
      "topright", legend = short_tx, col = dyn_pal[present_tx],
      pch = dyn_pch[present_tx], pt.bg = NA, cex = 0.72, bty = "n"
    )

    # Panel C: Top Edaphic/Lithic Predictor Boxplots + Kruskal-Wallis p-value
    ed_v <- suppressWarnings(as.numeric(df[[top_edaph_1]]))
    ed_v[!is.finite(ed_v)] <- stats::median(ed_v[is.finite(ed_v)], na.rm = TRUE)
    ed_v[!is.finite(ed_v)] <- 0.0
    kw_p <- tryCatch(stats::kruskal.test(ed_v ~ factor(tx_clean))$p.value, error = function(e) NA_real_)
    p_str <- if (is.finite(kw_p)) sprintf("Kruskal-Wallis p = %.2e", kw_p) else "Edaphic Divergence"
    tx_fac <- factor(tx_clean, levels = present_tx, labels = short_tx)
    graphics::boxplot(
      ed_v ~ tx_fac,
      col = grDevices::adjustcolor(dyn_pal[present_tx], alpha.f = 0.25),
      border = dyn_pal[present_tx], outpch = NA,
      xlab = "", ylab = top_edaph_lbl,
      main = sprintf("C. %s (%s)", top_edaph_lbl, p_str),
      font.axis = 3, cex.axis = 0.80, font.main = 2
    )
    graphics::grid(nx = NA, ny = NULL, col = "grey88")
    for (s_i in seq_along(present_tx)) {
      m_s <- which(tx_clean == present_tx[s_i])
      if (length(m_s) > 0L) {
        jit_x <- s_i + (((seq_along(m_s) %% 7L) - 3.0) * 0.045)
        graphics::points(
          jit_x, ed_v[m_s],
          pch = dyn_pch[[present_tx[s_i]]], bg = NA,
          col = grDevices::adjustcolor(dyn_pal[[present_tx[s_i]]], alpha.f = alpha_val),
          cex = 0.75, lwd = 0.9
        )
      }
    }

    # Panel D: Pairwise Schoener's D Niche Overlap Matrix (Benjamini-Hochberg FDR Adjusted)
    n_sp <- length(present_tx)
    d_mat <- matrix(1.0, nrow = n_sp, ncol = n_sp, dimnames = list(short_tx, short_tx))
    p_mat <- matrix(NA_real_, nrow = n_sp, ncol = n_sp, dimnames = list(short_tx, short_tx))
    niche_df <- if (!is.null(niche_res) && is.list(niche_res) && !is.data.frame(niche_res) && !is.null(niche_res$summary)) {
      niche_res$summary
    } else if (is.data.frame(niche_res)) {
      niche_res
    } else {
      NULL
    }
    if (!is.null(niche_df) && nrow(niche_df) > 0L) {
      p_col_use <- intersect(c("p_value_D", "Equivalency_p_value"), names(niche_df))[1]
      p_adj_vec <- if (!is.na(p_col_use)) stats::p.adjust(as.numeric(niche_df[[p_col_use]]), method = "BH") else rep(NA_real_, nrow(niche_df))
      t1_col <- intersect(c("Taxon1", "Taxon_1"), names(niche_df))[1]
      t2_col <- intersect(c("Taxon2", "Taxon_2"), names(niche_df))[1]
      d_col <- intersect(c("Schoener_D", "Schoeners_D", "Warren_D"), names(niche_df))[1]
      if (!is.na(t1_col) && !is.na(t2_col) && !is.na(d_col)) {
        for (r_i in seq_len(nrow(niche_df))) {
          t1 <- sub("^Packera\\s+", "P. ", as.character(niche_df[[t1_col]][r_i]))
          t2 <- sub("^Packera\\s+", "P. ", as.character(niche_df[[t2_col]][r_i]))
          dv <- suppressWarnings(as.numeric(niche_df[[d_col]][r_i]))
          pv <- p_adj_vec[r_i]
          if (t1 %in% short_tx && t2 %in% short_tx && is.finite(dv)) {
            d_mat[t1, t2] <- dv; d_mat[t2, t1] <- dv
            p_mat[t1, t2] <- pv; p_mat[t2, t1] <- pv
          }
        }
      }
    }
    graphics::plot(
      c(0.5, n_sp + 0.5), c(0.5, n_sp + 0.5), type = "n", axes = FALSE,
      xlab = "Active Group", ylab = "",
      main = "D. Pairwise Schoener's D Niche Overlap & Equivalency", font.main = 2
    )
    ramp_d <- grDevices::colorRampPalette(c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))(101)
    for (i in seq_len(n_sp)) {
      y_pos <- n_sp - i + 1L
      for (j in seq_len(n_sp)) {
        dv <- d_mat[i, j]
        c_idx <- max(1L, min(101L, round(dv * 100) + 1L))
        graphics::rect(j - 0.48, y_pos - 0.46, j + 0.48, y_pos + 0.46, col = if (i == j) "#EEEEEE" else ramp_d[c_idx], border = "white", lwd = 1.5)
        lbl_cell <- if (i == j) "1.00\n(Self)" else if (is.finite(p_mat[i, j])) sprintf("D = %.2f\n(q = %.3f)", dv, p_mat[i, j]) else sprintf("D = %.2f", dv)
        graphics::text(j, y_pos, labels = lbl_cell, cex = 0.74, font = if (i == j) 1 else 2, col = if (i != j && dv > 0.62) "white" else "#111111")
      }
    }
    graphics::axis(1, at = seq_len(n_sp), labels = short_tx, tick = FALSE, font = 3, cex.axis = 0.80)
    graphics::axis(2, at = rev(seq_len(n_sp)), labels = short_tx, las = 1, tick = FALSE, font = 3, cex.axis = 0.80)
    graphics::box(col = "grey75")
    u_sp_hdr <- get_active_packera_taxa(df$species_standardized)
    hdr_lbl <- if (length(u_sp_hdr) == 1L) u_sp_hdr[1] else sprintf("Packera Multi-Taxon Cohort (%d Active Groups)", length(present_tx))
    graphics::title(sprintf("%s: Multi-Scale Ecological Niche, Top RF Predictors & Equivalency", hdr_lbl), outer = TRUE, cex.main = 1.20)
  }

  pdf_out <- file.path(fdir, "fig06_ecological_niche_phenology.pdf")
  pdf_tmp <- paste0(pdf_out, ".tmp")
  grDevices::pdf(pdf_tmp, width = 13.0, height = 9.5)
  tryCatch(render_plate(), finally = grDevices::dev.off())
  if (file.exists(pdf_tmp)) file.rename(pdf_tmp, pdf_out)

  png_out <- file.path(fdir, "fig06_ecological_niche_phenology.png")
  png_tmp <- paste0(png_out, ".tmp")
  grDevices::png(png_tmp, width = 13.0, height = 9.5, units = "in", res = 300)
  tryCatch(render_plate(), finally = grDevices::dev.off())
  if (file.exists(png_tmp)) file.rename(png_tmp, png_out)
  invisible(list(pdf = pdf_out, png = png_out))
}

export_spatial_rf_figures <- function(df, srf_res, niche_res, out_pdf, pheno_res = NULL) {
  message(sprintf("Exporting publication-quality 5-panel figure to %s...", out_pdf))
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  
  # Theme setup
  thm <- ggplot2::theme_bw(base_size = 9) +
         ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10),
                        plot.subtitle = ggplot2::element_text(size = 8.2, color = "grey25"),
                        panel.grid.minor = ggplot2::element_blank())

  # Panel A: Spatial RF Variable Importance (Human-Readable Ecological Labels & Okabe-Ito Palette)
  imp_df <- srf_res$importance
  imp_df$Label <- if (exists("format_step06_predictor_labels", mode = "function")) {
    format_step06_predictor_labels(as.character(imp_df$Predictor))
  } else if (exists("format_predictor_labels", mode = "function")) {
    format_predictor_labels(as.character(imp_df$Predictor))
  } else {
    as.character(imp_df$Predictor)
  }
  p1 <- ggplot2::ggplot(imp_df, ggplot2::aes(x = stats::reorder(Label, MeanDecreaseAccuracy),
                                             y = MeanDecreaseAccuracy, fill = Type)) +
    ggplot2::geom_col(width = 0.72, color = "#222222", linewidth = 0.25) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(
      "Fine-Grain Edaphic & Outcrop (10-30m)" = "#009E73",
      "Environmental" = "#0072B2",
      "Spatial (MEMs)" = "#CC79A7"
    )) +
    thm + ggplot2::labs(
      title = "A. Spatial RF Variable Importance (10m-30m Edaphic + Climate + MEMs)",
      subtitle = sprintf("OOB Accuracy = %.1f%% | Spatial Block CV Accuracy = %.1f%%",
                         100 * (1 - srf_res$OOB_Error), 100 * srf_res$Spatial_Block_CV_Accuracy),
      x = NULL, y = "Mean Decrease Accuracy (OOB)", fill = "Predictor Class"
    )

  # Panel B: Dynamic Top-2 Non-MEM Spatial RF Environmental Predictors (Hollow Markers + 85% Ellipses)
  df$plot_group <- if (exists("resolve_plot_group_column", mode = "function")) {
    resolve_plot_group_column(df)
  } else {
    sp_v <- as.character(df$species_standardized)
    if (length(unique(stats::na.omit(sp_v))) < 2L ||
        ("analysis_group" %in% names(df) && any(grepl("Morphotype_K|^LatTertile_|Ecotype|Subcluster", df$analysis_group, ignore.case = TRUE))) ||
        ("grouping_col" %in% names(df) && any(grepl("Morphotype_K|^LatTertile_|Ecotype|Subcluster", df$grouping_col, ignore.case = TRUE)))) {
      if ("analysis_group" %in% names(df) && length(unique(stats::na.omit(df$analysis_group))) >= 2L) as.character(df$analysis_group)
      else if ("grouping_col" %in% names(df) && length(unique(stats::na.omit(df$grouping_col))) >= 2L) as.character(df$grouping_col)
      else if ("gmm_cluster" %in% names(df) && any(!is.na(df$gmm_cluster))) paste0("Morphotype_K", ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster)))
      else sp_v
    } else sp_v
  }
  active_tx <- sort(unique(stats::na.omit(df$plot_group)))
  active_tx <- active_tx[nzchar(active_tx) & active_tx != "Unknown"]
  if (length(active_tx) == 0) active_tx <- get_active_packera_taxa(df$species_standardized)
  is_subcluster_mode <- any(grepl("Morphotype_K|\\[K[0-9]+\\]|^LatTertile_|Ecotype|Subcluster", active_tx, ignore.case = TRUE)) ||
    length(unique(stats::na.omit(df$species_standardized))) < 2L
  grp_lbl <- if (is_subcluster_mode) "Ecotype / Morphotype" else "Taxon"
  dyn_pal <- resolve_dynamic_group_palette(active_tx)
  dyn_shapes <- stats::setNames(rep_len(c(21, 24, 22, 23, 1, 2, 5, 0, 6, 4), length(active_tx)), active_tx)
  plot_df <- df[df$plot_group %in% active_tx, , drop = FALSE]
  if (!("fine_soil_sand" %in% names(plot_df)) && "soil_sand" %in% names(plot_df)) plot_df$fine_soil_sand <- plot_df$soil_sand
  if (!("fine_soil_ph" %in% names(plot_df)) && "soil_ph" %in% names(plot_df)) plot_df$fine_soil_ph <- plot_df$soil_ph
  if (!("fine_bedrock_depth_cm" %in% names(plot_df))) plot_df$fine_bedrock_depth_cm <- 150.0
  if (!("ssurgo_rock_outcrop_pct" %in% names(plot_df))) plot_df$ssurgo_rock_outcrop_pct <- 0.0
  if (!("topo_tpi_30m" %in% names(plot_df))) plot_df$topo_tpi_30m <- 0.0

  top_sel <- select_top_rf_predictors(imp_df, plot_df)
  top_env_1 <- top_sel$top_env_1
  top_env_2 <- top_sel$top_env_2
  top_edaph_1 <- top_sel$top_edaph_1

  env_axis_lbls <- format_predictor_labels(c(top_env_1, top_env_2))
  if (exists("format_step06_predictor_labels", mode = "function")) {
    env_axis_lbls <- format_step06_predictor_labels(c(top_env_1, top_env_2))
  }
  edaph_lbl <- format_predictor_labels(top_edaph_1)[1]
  if (exists("format_step06_predictor_labels", mode = "function")) {
    edaph_lbl <- format_step06_predictor_labels(top_edaph_1)[1]
  }

  plot_df$top_env_1_val <- suppressWarnings(as.numeric(plot_df[[top_env_1]]))
  plot_df$top_env_2_val <- suppressWarnings(as.numeric(plot_df[[top_env_2]]))
  plot_df$top_edaph_1_val <- suppressWarnings(as.numeric(plot_df[[top_edaph_1]]))

  pt_alpha <- if (nrow(plot_df) > 1000) 0.45 else 0.78
  p2 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = top_env_1_val, y = top_env_2_val, color = plot_group, shape = plot_group)) +
    ggplot2::stat_ellipse(level = 0.85, linewidth = 1.0) +
    ggplot2::geom_point(fill = NA, stroke = 0.75, alpha = pt_alpha, size = 1.55) +
    ggplot2::scale_shape_manual(values = dyn_shapes) +
    ggplot2::scale_color_manual(values = dyn_pal) +
    thm + ggplot2::labs(
      title = "B. Top Spatial RF Ecological Niche Envelopes (85% Confidence Ellipses)",
      subtitle = sprintf("Top Non-MEM Predictors: %s vs. %s", env_axis_lbls[1], env_axis_lbls[2]),
      x = env_axis_lbls[1], y = env_axis_lbls[2], color = grp_lbl, shape = grp_lbl
    )
  if ("is_type_specimen" %in% names(plot_df) && any(as.logical(plot_df$is_type_specimen) == TRUE, na.rm = TRUE)) {
    p2 <- p2 + ggplot2::geom_point(
      data = plot_df[which(as.logical(plot_df$is_type_specimen) == TRUE), , drop = FALSE],
      ggplot2::aes(x = top_env_1_val, y = top_env_2_val),
      shape = 5, fill = NA, stroke = 1.3, size = 2.8, color = "#111111", inherit.aes = FALSE, show.legend = FALSE
    )
  }

  # Panel C: Top-Importance Edaphic/Lithic Predictor Diagnostic Boxplots + Kruskal-Wallis p-value
  kw_edaph <- tryCatch(stats::kruskal.test(top_edaph_1_val ~ factor(plot_group), data = plot_df)$p.value, error = function(e) NA_real_)
  p3 <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = plot_group, y = top_edaph_1_val, fill = plot_group, color = plot_group)) +
    ggplot2::geom_boxplot(alpha = 0.30, outlier.shape = 1, outlier.size = 1.2, width = 0.58) +
    ggplot2::scale_fill_manual(values = dyn_pal) +
    ggplot2::scale_color_manual(values = dyn_pal) +
    thm + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1, face = if (is_subcluster_mode) "bold" else "italic"), legend.position = "none") +
    ggplot2::labs(
      title = sprintf("C. %s", edaph_lbl),
      subtitle = if (!is.na(kw_edaph) && is.finite(kw_edaph)) sprintf("Kruskal-Wallis p = %.2e across active groups", kw_edaph) else sprintf("Top-ranked edaphic/lithic predictor (%s)", top_edaph_1),
      x = NULL, y = edaph_lbl
    )

  # Panel D: Warren's Niche Identity Permutations & Empirical Schoener's D
  if (!is.null(niche_res$permutations) && is.data.frame(niche_res$permutations) && nrow(niche_res$permutations) > 0L) {
    p4 <- ggplot2::ggplot(niche_res$permutations, ggplot2::aes(x = D_null, fill = Pair)) +
      ggplot2::geom_histogram(bins = 20, alpha = 0.65, position = "identity", color = "white", linewidth = 0.2) +
      ggplot2::geom_vline(data = niche_res$summary, ggplot2::aes(xintercept = Warren_D, color = Pair),
                          linetype = "dashed", linewidth = 0.95, show.legend = FALSE) +
      ggplot2::coord_cartesian(xlim = c(0, 1)) +
      thm + ggplot2::labs(
        title = "D. Pairwise Niche Equivalency (Schoener's D vs. Null Permutations)",
        subtitle = "Dashed vertical lines = Empirical Schoener's D; Histogram = Null niche identity distribution",
        x = "Schoener's D Niche Overlap (0 = Disjoint, 1 = Identical)", y = "Permutation Count", fill = "Taxon Pair"
      )
  } else {
    p4 <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Pairwise niche equivalency requires >= 2 groups with N >= 5") +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      thm + ggplot2::labs(
        title = "D. Pairwise Niche Equivalency (Schoener's D vs. Null Permutations)",
        x = "Schoener's D Niche Overlap", y = "Permutation Count"
      )
  }

  # Panel E: Horizontal Stacked Variance-Partitioning Bar Chart (Replaces Pie Chart)
  vp_df <- srf_res$variance_partition
  vp_df$Label <- ifelse(vp_df$Variance_Pct >= 5, sprintf("%.1f%%", vp_df$Variance_Pct), "")
  vp_df$Component <- factor(vp_df$Component, levels = rev(vp_df$Component))
  vp_cols <- stats::setNames(c("#009E73", "#56B4E9", "#0072B2", "#CC79A7", "#999999")[seq_len(nrow(vp_df))], levels(vp_df$Component))
  p5 <- ggplot2::ggplot(vp_df, ggplot2::aes(x = "Variance\nPartition", y = Variance_Pct, fill = Component)) +
    ggplot2::geom_col(width = 0.52, color = "white", linewidth = 0.6) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(ggplot2::aes(label = Label), position = ggplot2::position_stack(vjust = 0.5),
                       size = 3.1, fontface = "bold", color = "white") +
    ggplot2::scale_fill_manual(values = vp_cols) +
    thm + ggplot2::theme(legend.position = "bottom", legend.direction = "vertical") +
    ggplot2::labs(
      title = "E. Multi-Scale Variance Partitioning",
      subtitle = "Relative contribution of [a] Edaphic/Outcrop, [b] Climate & [c] Spatial MEMs",
      x = NULL, y = "Percentage of Explained Niche Variance (%)", fill = "Variance Component"
    )

  tmp_pdf <- paste0(out_pdf, ".tmp")
  pdf(tmp_pdf, width = 13, height = 8.5)
  tryCatch({
    gridExtra::grid.arrange(p1, p2, p3, p4, p5, layout_matrix = rbind(c(1, 2, 5), c(3, 4, 5)), widths = c(1.15, 1.15, 0.95))
  }, finally = {
    dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, out_pdf)

  # Export auxiliary Step 06 figures via generate_step06_pub_figures if available,
  # then atomically render dynamic top-RF-importance fig06_ecological_niche_phenology.{pdf,png}
  if (exists("generate_step06_pub_figures", mode = "function")) {
    tryCatch({
      format_predictor_labels <<- format_step06_predictor_labels
      generate_step06_pub_figures(df, srf_res, niche_res, pheno_res, out_pdf)
    }, error = function(e) warning("generate_step06_pub_figures warning: ", e$message))
  }
  tryCatch(
    export_dynamic_fig06_ecological_niche_phenology(df, srf_res, niche_res, out_pdf),
    error = function(e) warning("export_dynamic_fig06_ecological_niche_phenology warning: ", e$message)
  )
}

# ------------------------------------------------------------------------------
# 6b. Default MaxEnt (maxnet / dismo / IPP-Cloglog GLM) & Spatial RF Ensemble SDM
#     Projection, Continuous Boyce Index (Hirzel et al. 2006), & Habitat Maps
# ------------------------------------------------------------------------------
compute_continuous_boyce_index <- function(suitability_all, suitability_pres, n_bins = 20L) {
  s_all <- suppressWarnings(as.numeric(suitability_all))
  s_all <- s_all[is.finite(s_all)]
  s_pres <- suppressWarnings(as.numeric(suitability_pres))
  s_pres <- s_pres[is.finite(s_pres)]
  empty_curve <- data.frame(suitability = numeric(0), pe_ratio = numeric(0), stringsAsFactors = FALSE)
  if (length(s_pres) < 2L || length(s_all) < 4L) {
    res <- NA_real_
    attr(res, "boyce_curve") <- empty_curve
    return(res)
  }
  s_min <- min(c(s_all, s_pres), na.rm = TRUE)
  s_max <- max(c(s_all, s_pres), na.rm = TRUE)
  if (!is.finite(s_min) || !is.finite(s_max) || (s_max - s_min) < 1e-8) {
    res <- 0.0
    attr(res, "boyce_curve") <- empty_curve
    return(res)
  }
  n_b <- max(5L, as.integer(n_bins))
  win_w <- max((s_max - s_min) / 10.0, 1e-4)
  centers <- seq(s_min + 0.5 * win_w, s_max - 0.5 * win_w, length.out = n_b)
  pe_ratio <- rep(NA_real_, n_b)
  n_p <- length(s_pres)
  n_a <- length(s_all)
  for (i in seq_len(n_b)) {
    lo <- centers[i] - 0.5 * win_w
    hi <- centers[i] + 0.5 * win_w
    p_i <- sum(s_pres >= lo & s_pres <= hi) / n_p
    e_i <- sum(s_all >= lo & s_all <= hi) / n_a
    if (e_i > 0) pe_ratio[i] <- p_i / e_i
  }
  ok <- is.finite(pe_ratio) & is.finite(centers)
  curve_df <- data.frame(suitability = centers[ok], pe_ratio = pe_ratio[ok], stringsAsFactors = FALSE)
  if (sum(ok) < 3L || stats::sd(pe_ratio[ok]) < 1e-10) {
    res <- 0.0
    attr(res, "boyce_curve") <- curve_df
    return(res)
  }
  cbi <- suppressWarnings(stats::cor(centers[ok], pe_ratio[ok], method = "spearman", use = "complete.obs"))
  if (!is.finite(cbi)) cbi <- 0.0
  res <- round(as.numeric(cbi), 4)
  attr(res, "boyce_curve") <- curve_df
  res
}

.compute_roc_auc_curve <- function(y01, probs, wts = NULL) {
  y_v <- suppressWarnings(as.numeric(y01))
  p_v <- suppressWarnings(as.numeric(probs))
  w_v <- if (is.null(wts) || length(wts) != length(y_v)) rep(1.0, length(y_v)) else suppressWarnings(as.numeric(wts))
  w_v[!is.finite(w_v) | w_v <= 0] <- 1.0
  ok <- is.finite(y_v) & is.finite(p_v) & (y_v %in% c(0, 1))
  if (sum(y_v[ok] == 1) < 1L || sum(y_v[ok] == 0) < 1L) {
    return(list(auc = NA_real_, fpr = c(0, 1), tpr = c(0, 1)))
  }
  y_v <- y_v[ok]; p_v <- p_v[ok]; w_v <- w_v[ok]
  ord <- order(p_v, decreasing = TRUE)
  y_s <- y_v[ord]; p_s <- p_v[ord]; w_s <- w_v[ord]
  pos_w <- sum(w_s[y_s == 1])
  neg_w <- sum(w_s[y_s == 0])
  if (pos_w <= 0 || neg_w <= 0) return(list(auc = 0.5, fpr = c(0, 1), tpr = c(0, 1)))
  tp_cum <- c(0, cumsum(ifelse(y_s == 1, w_s, 0)) / pos_w)
  fp_cum <- c(0, cumsum(ifelse(y_s == 0, w_s, 0)) / neg_w)
  dups <- !duplicated(c(Inf, p_s), fromLast = TRUE)
  tpr <- tp_cum[dups]
  fpr <- fp_cum[dups]
  if ( utils::tail(fpr, 1) < 1.0 || utils::tail(tpr, 1) < 1.0 ) {
    fpr <- c(fpr, 1.0); tpr <- c(tpr, 1.0)
  }
  auc_val <- sum(0.5 * (tpr[-1] + tpr[-length(tpr)]) * diff(fpr))
  list(auc = round(min(1.0, max(0.0, auc_val)), 4), fpr = fpr, tpr = tpr)
}

.fit_single_maxent_model <- function(p_vec, X_df, w_vec = NULL, use_java_maxent = TRUE) {
  p_num <- as.integer(as.numeric(p_vec) > 0.5)
  X_mat <- as.data.frame(X_df)
  pred_cols <- colnames(X_mat)
  n_rows <- nrow(X_mat)
  w_use <- if (is.null(w_vec) || length(w_vec) != n_rows) rep(1.0, n_rows) else as.numeric(w_vec)
  w_use[!is.finite(w_use) | w_use <= 0] <- 1.0
  n_pres <- sum(p_num == 1L)

  # 1. Primary Engine: maxnet (Phillips et al. 2017) with lqh -> lq -> l cascade
  if (requireNamespace("maxnet", quietly = TRUE) && n_pres >= 3L && ncol(X_mat) >= 1L) {
    # Weight high-authority & Nomenclatural Type presences via importance replication
    pres_idx <- which(p_num == 1L)
    bg_idx <- which(p_num == 0L)
    w_p <- w_use[pres_idx]
    w_p_max <- max(w_p, 1e-6)
    rep_p <- pmax(1L, pmin(3L, as.integer(round(1.0 + 1.5 * (w_p / w_p_max)))))
    exp_pres <- rep(pres_idx, times = rep_p)
    idx_mx <- c(exp_pres, bg_idx)
    p_mx <- p_num[idx_mx]
    X_mx <- X_mat[idx_mx, , drop = FALSE]
    cls_cascade <- if (n_pres >= 15L) c("lqh", "lq", "l") else if (n_pres >= 8L) c("lq", "l") else c("l")
    for (cl_opt in cls_cascade) {
      mod_mx <- tryCatch({
        f_mx <- maxnet::maxnet.formula(p_mx, X_mx, classes = cl_opt)
        suppressWarnings(maxnet::maxnet(p = p_mx, data = X_mx, f = f_mx))
      }, error = function(e) NULL)
      if (!is.null(mod_mx)) {
        test_pr <- tryCatch(
          as.numeric(stats::predict(mod_mx, newdata = X_mat[seq_len(min(5L, n_rows)), , drop = FALSE], type = "cloglog")),
          error = function(e) NULL
        )
        if (!is.null(test_pr) && any(is.finite(test_pr))) {
          return(list(
            engine = sprintf("maxnet (%s/cloglog)", cl_opt),
            model = mod_mx,
            predict_cloglog = function(new_df) {
              nd <- as.data.frame(new_df)[, pred_cols, drop = FALSE]
              for (cn in pred_cols) {
                v <- suppressWarnings(as.numeric(nd[[cn]]))
                v[!is.finite(v)] <- stats::median(X_mat[[cn]], na.rm = TRUE)
                nd[[cn]] <- v
              }
              pr <- tryCatch(as.numeric(stats::predict(mod_mx, newdata = nd, type = "cloglog")), error = function(e) rep(0.5, nrow(nd)))
              pr[!is.finite(pr)] <- 0.0
              pmin(1.0, pmax(0.0, pr))
            }
          ))
        }
      }
    }
  }

  # 2. Secondary Engine: dismo::maxent when maxent.jar is available
  if (isTRUE(use_java_maxent) && requireNamespace("dismo", quietly = TRUE) && requireNamespace("rJava", quietly = TRUE)) {
    jar_p <- system.file("java", "maxent.jar", package = "dismo")
    if (nzchar(jar_p) && file.exists(jar_p) && n_pres >= 4L) {
      mx_args <- c("nothreshold", "noproduct", "noplots", "nopictures", "threads=4")
      if (n_pres < 15L) mx_args <- c(mx_args, "nohinge")
      mod_dm <- tryCatch(suppressWarnings(dismo::maxent(x = X_mat, p = p_num, args = mx_args)), error = function(e) NULL)
      if (!is.null(mod_dm)) {
        return(list(
          engine = "dismo::maxent (cloglog)",
          model = mod_dm,
          predict_cloglog = function(new_df) {
            nd <- as.data.frame(new_df)[, pred_cols, drop = FALSE]
            for (cn in pred_cols) {
              v <- suppressWarnings(as.numeric(nd[[cn]]))
              v[!is.finite(v)] <- stats::median(X_mat[[cn]], na.rm = TRUE)
              nd[[cn]] <- v
            }
            pr <- tryCatch(as.numeric(dismo::predict(mod_dm, nd, args = c("outputformat=cloglog"))), error = function(e) rep(0.5, nrow(nd)))
            if (length(pr) != nrow(nd)) pr <- rep(0.5, nrow(nd))
            pr[!is.finite(pr)] <- 0.0
            pmin(1.0, pmax(0.0, pr))
          }
        ))
      }
    }
  }

  # 3. Deterministic Fallback Engine: Authority-Weighted Inhomogeneous Poisson Process (IPP)
  #    Quadratic-Hinge Cloglog Regularized GLM (Phillips et al. 2017 equivalence)
  mu_vec <- vapply(pred_cols, function(cn) mean(as.numeric(X_mat[[cn]]), na.rm = TRUE), numeric(1))
  sd_vec <- vapply(pred_cols, function(cn) max(stats::sd(as.numeric(X_mat[[cn]]), na.rm = TRUE), 1e-4), numeric(1))
  q50_vec <- vapply(pred_cols, function(cn) stats::median(as.numeric(X_mat[[cn]]), na.rm = TRUE), numeric(1))
  use_quad <- n_pres >= 6L
  use_hinge <- n_pres >= 15L

  .build_basis <- function(df_in) {
    Z <- matrix(0.0, nrow = nrow(df_in), ncol = length(pred_cols), dimnames = list(NULL, pred_cols))
    for (j in seq_along(pred_cols)) {
      cn <- pred_cols[j]
      v <- if (cn %in% names(df_in)) suppressWarnings(as.numeric(df_in[[cn]])) else rep(mu_vec[j], nrow(df_in))
      v[!is.finite(v)] <- mu_vec[j]
      Z[, j] <- pmin(4.5, pmax(-4.5, (v - mu_vec[j]) / sd_vec[j]))
    }
    B <- Z
    if (isTRUE(use_quad)) {
      Z2 <- Z^2
      dim(Z2) <- dim(Z)
      colnames(Z2) <- paste0(pred_cols, "__q2")
      B <- cbind(B, Z2)
    }
    if (isTRUE(use_hinge)) {
      Zh <- sweep(Z, 2L, (q50_vec - mu_vec) / sd_vec, "-")
      Zh[Zh < 0.0] <- 0.0
      dim(Zh) <- dim(Z)
      colnames(Zh) <- paste0(pred_cols, "__hg")
      B <- cbind(B, Zh)
    }
    B
  }

  Phi <- .build_basis(X_mat)
  # Balance total background weight against authority-weighted presences
  w_glm <- w_use
  sum_wp <- max(sum(w_glm[p_num == 1L]), 1e-6)
  sum_wb <- max(sum(w_glm[p_num == 0L]), 1e-6)
  w_glm[p_num == 0L] <- w_glm[p_num == 0L] * (2.0 * sum_wp / sum_wb)

  if (requireNamespace("glmnet", quietly = TRUE) && ncol(Phi) >= 2L && n_pres >= 4L && sum(p_num == 0L) >= 4L) {
    fit_gn <- tryCatch(
      suppressWarnings(glmnet::glmnet(x = Phi, y = p_num, family = "binomial", weights = w_glm, alpha = 0.85)),
      error = function(e) NULL
    )
    if (!is.null(fit_gn) && length(fit_gn$lambda) > 0L) {
      s_use <- fit_gn$lambda[max(1L, floor(0.65 * length(fit_gn$lambda)))]
      eta_tr <- as.numeric(stats::predict(fit_gn, newx = Phi, s = s_use, type = "link"))
      # Calibrate entropy offset so median presence cloglog suitability = 1 - exp(-1) ~ 0.632 (Phillips et al. 2017)
      med_eta_p <- stats::median(eta_tr[p_num == 1L], na.rm = TRUE)
      if (!is.finite(med_eta_p)) med_eta_p <- 0.0
      c_ent <- -med_eta_p
      return(list(
        engine = sprintf("glmnet_maxent (%s/cloglog)", if (use_hinge) "lqh" else if (use_quad) "lq" else "l"),
        model = fit_gn,
        predict_cloglog = function(new_df) {
          B_new <- .build_basis(as.data.frame(new_df))
          eta <- tryCatch(as.numeric(stats::predict(fit_gn, newx = B_new, s = s_use, type = "link")) + c_ent, error = function(e) rep(0.0, nrow(B_new)))
          eta <- pmin(10.0, pmax(-15.0, eta))
          pr <- 1.0 - exp(-exp(eta))
          pr[!is.finite(pr)] <- 0.0
          pmin(1.0, pmax(0.0, pr))
        }
      ))
    }
  }

  # Ridge-stabilized weighted cloglog GLM fallback
  X_des <- cbind(Intercept = 1.0, Phi)
  W_diag <- w_glm / max(mean(w_glm), 1e-6)
  P_pen <- diag(c(0.0, rep(0.25, ncol(Phi))), nrow = ncol(X_des), ncol = ncol(X_des))
  y_target <- ifelse(p_num == 1L, 0.85, -1.25)
  beta_hat <- tryCatch(
    solve(crossprod(X_des * sqrt(W_diag)) + P_pen, crossprod(X_des, W_diag * y_target)),
    error = function(e) rep(0.0, ncol(X_des))
  )
  beta_hat <- as.numeric(beta_hat)
  eta_tr <- as.numeric(X_des %*% beta_hat)
  med_eta_p <- stats::median(eta_tr[p_num == 1L], na.rm = TRUE)
  if (!is.finite(med_eta_p)) med_eta_p <- 0.0
  c_ent <- -med_eta_p

  list(
    engine = sprintf("weighted_glm_maxent (%s/cloglog)", if (use_quad) "lq" else "l"),
    model = list(beta = beta_hat, c_ent = c_ent),
    predict_cloglog = function(new_df) {
      B_new <- cbind(Intercept = 1.0, .build_basis(as.data.frame(new_df)))
      eta <- pmin(10.0, pmax(-15.0, as.numeric(B_new %*% beta_hat) + c_ent))
      pr <- 1.0 - exp(-exp(eta))
      pr[!is.finite(pr)] <- 0.0
      pmin(1.0, pmax(0.0, pr))
    }
  )
}

.fit_single_srf_binary_model <- function(p_vec, X_df, w_vec = NULL, n_trees = 250L) {
  p_num <- as.integer(as.numeric(p_vec) > 0.5)
  X_mat <- as.data.frame(X_df)
  pred_cols <- colnames(X_mat)
  n_rows <- nrow(X_mat)
  w_use <- if (is.null(w_vec) || length(w_vec) != n_rows) rep(1.0, n_rows) else as.numeric(w_vec)
  w_use[!is.finite(w_use) | w_use <= 0] <- 1.0
  # Inverse class-frequency weighting (Valavi et al. 2021)
  n_p <- max(sum(p_num == 1L), 1L)
  n_b <- max(sum(p_num == 0L), 1L)
  w_rf <- w_use * ifelse(p_num == 1L, n_rows / (2.0 * n_p), n_rows / (2.0 * n_b))
  y_fac <- factor(ifelse(p_num == 1L, "Presence", "Background"), levels = c("Background", "Presence"))

  if (requireNamespace("ranger", quietly = TRUE) && sum(p_num == 1L) >= 2L && sum(p_num == 0L) >= 2L) {
    rf_threads <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(16L) else 1L
    tr_df <- data.frame(y = y_fac, X_mat, check.names = FALSE)
    fit_rg <- tryCatch(
      ranger::ranger(y ~ ., data = tr_df, num.trees = n_trees, probability = TRUE, case.weights = w_rf, num.threads = rf_threads, seed = 42),
      error = function(e) NULL
    )
    if (!is.null(fit_rg)) {
      return(list(
        model = fit_rg,
        predict_prob = function(new_df) {
          nd <- as.data.frame(new_df)[, pred_cols, drop = FALSE]
          for (cn in pred_cols) {
            v <- suppressWarnings(as.numeric(nd[[cn]]))
            v[!is.finite(v)] <- stats::median(X_mat[[cn]], na.rm = TRUE)
            nd[[cn]] <- v
          }
          pr_m <- tryCatch(stats::predict(fit_rg, data = nd, num.threads = rf_threads)$predictions, error = function(e) NULL)
          if (is.matrix(pr_m) && "Presence" %in% colnames(pr_m)) {
            pr <- as.numeric(pr_m[, "Presence"])
          } else {
            pr <- rep(0.5, nrow(nd))
          }
          pr[!is.finite(pr)] <- 0.0
          pmin(1.0, pmax(0.0, pr))
        }
      ))
    }
  }
  if (requireNamespace("randomForest", quietly = TRUE) && sum(p_num == 1L) >= 2L && sum(p_num == 0L) >= 2L) {
    fit_rf <- tryCatch(randomForest::randomForest(x = X_mat, y = y_fac, ntree = n_trees), error = function(e) NULL)
    if (!is.null(fit_rf)) {
      return(list(
        model = fit_rf,
        predict_prob = function(new_df) {
          nd <- as.data.frame(new_df)[, pred_cols, drop = FALSE]
          pr_m <- tryCatch(stats::predict(fit_rf, newdata = nd, type = "prob"), error = function(e) NULL)
          pr <- if (is.matrix(pr_m) && "Presence" %in% colnames(pr_m)) as.numeric(pr_m[, "Presence"]) else rep(0.5, nrow(nd))
          pr[!is.finite(pr)] <- 0.0
          pmin(1.0, pmax(0.0, pr))
        }
      ))
    }
  }
  # Fallback to MaxEnt cloglog if neither ranger nor randomForest is available
  mx_fb <- .fit_single_maxent_model(p_vec, X_df, w_vec)
  list(model = mx_fb$model, predict_prob = mx_fb$predict_cloglog)
}

.build_regional_sdm_raster_stack <- function(df, env_dir, kept_vars, ncols = 110L, nrows = 90L) {
  if (!requireNamespace("terra", quietly = TRUE)) return(NULL)
  lon_v <- suppressWarnings(as.numeric(df$longitude))
  lat_v <- suppressWarnings(as.numeric(df$latitude))
  ok_ll <- is.finite(lon_v) & is.finite(lat_v)
  if (!any(ok_ll)) return(NULL)
  if (sum(ok_ll) >= 20L) {
    lon_rng <- as.numeric(stats::quantile(lon_v[ok_ll], probs = c(0.005, 0.995), na.rm = TRUE))
    lat_rng <- as.numeric(stats::quantile(lat_v[ok_ll], probs = c(0.005, 0.995), na.rm = TRUE))
  } else {
    lon_rng <- range(lon_v[ok_ll], na.rm = TRUE)
    lat_rng <- range(lat_v[ok_ll], na.rm = TRUE)
  }
  pad_x <- max(0.75, 0.10 * diff(lon_rng))
  pad_y <- max(0.75, 0.10 * diff(lat_rng))
  med_lat <- stats::median(lat_v[ok_ll], na.rm = TRUE)
  min_lat_clamp <- if (is.finite(med_lat) && med_lat >= 28.0) 24.5 else 14.0
  xmin <- max(-130.0, lon_rng[1] - pad_x)
  xmax <- min(-50.0,  lon_rng[2] + pad_x)
  ymin <- max(min_lat_clamp, lat_rng[1] - pad_y)
  ymax <- min(66.0,   lat_rng[2] + pad_y)
  if (xmax - xmin < 1.0) { xmin <- xmin - 0.75; xmax <- xmax + 0.75 }
  if (ymax - ymin < 1.0) { ymin <- ymin - 0.75; ymax <- ymax + 0.75 }

  study_ext <- terra::ext(xmin, xmax, ymin, ymax)
  tmpl_rast <- terra::rast(ext = study_ext, ncols = ncols, nrows = nrows, crs = "EPSG:4326")

  # Discover any local cropped rasters in env_dir
  loaded_layers <- list()
  if (!is.null(env_dir) && nzchar(env_dir) && dir.exists(env_dir)) {
    tif_files <- discover_env_rasters(env_dir, recursive = TRUE, crop_global = TRUE)
    if (length(tif_files) > 0L) {
      for (tf in tif_files) {
        r_t <- tryCatch(terra::rast(tf), error = function(e) NULL)
        if (!is.null(r_t)) {
          nms_t <- names(r_t)
          for (li in seq_along(nms_t)) {
            lnm <- nms_t[li]
            mapped_nm <- lnm
            if (!(mapped_nm %in% kept_vars)) {
              m_bio <- regmatches(lnm, regexec("bio[_.-]?([0-9]{1,2})$", tolower(lnm)))[[1]]
              if (length(m_bio) == 2L && m_bio[2] %in% names(WC_BIO_LOOKUP)) {
                mapped_nm <- WC_BIO_LOOKUP[[m_bio[2]]]
              }
            }
            if (mapped_nm %in% kept_vars && !(mapped_nm %in% names(loaded_layers))) {
              r_sub <- tryCatch({
                rc <- terra::crop(r_t[[li]], study_ext)
                rs <- terra::resample(rc, tmpl_rast, method = "bilinear")
                names(rs) <- mapped_nm
                rs
              }, error = function(e) NULL)
              if (!is.null(r_sub)) loaded_layers[[mapped_nm]] <- r_sub
            }
          }
        }
      }
    }
  }

  land_mask_rast <- if (length(loaded_layers) > 0L) loaded_layers[[1L]] else NULL

  # Interpolate any remaining kept_vars across template_rast via fast vectorized IDW
  miss_vars <- setdiff(kept_vars, names(loaded_layers))
  if (length(miss_vars) > 0L) {
    u_key <- paste(round(lon_v[ok_ll], 4), round(lat_v[ok_ll], 4), sep = "_")
    u_idx <- which(ok_ll)[!duplicated(u_key)]
    if (length(u_idx) > 450L) {
      u_idx <- u_idx[unique(round(seq(1L, length(u_idx), length.out = 450L)))]
    }
    vx <- lon_v[u_idx]
    vy <- lat_v[u_idx]
    V_mat <- matrix(0.0, nrow = length(u_idx), ncol = length(miss_vars), dimnames = list(NULL, miss_vars))
    for (j in seq_along(miss_vars)) {
      vn <- miss_vars[j]
      vv <- if (vn %in% names(df)) suppressWarnings(as.numeric(df[[vn]][u_idx])) else rep(0.0, length(u_idx))
      med_v <- stats::median(vv[is.finite(vv)], na.rm = TRUE)
      if (!is.finite(med_v)) med_v <- 0.0
      vv[!is.finite(vv)] <- med_v
      V_mat[, j] <- vv
    }
    g_xy <- terra::xyFromCell(tmpl_rast, seq_len(terra::ncell(tmpl_rast)))
    D2 <- outer(g_xy[, 1], vx, "-")^2 + outer(g_xy[, 2], vy, "-")^2 + 0.0025
    W <- 1.0 / (D2^1.15)
    W <- W / pmax(rowSums(W), 1e-12)
    G_mat <- W %*% V_mat
    for (j in seq_along(miss_vars)) {
      vn <- miss_vars[j]
      r_idw <- terra::rast(tmpl_rast)
      terra::values(r_idw) <- as.numeric(G_mat[, j])
      if (!is.null(land_mask_rast)) {
        r_idw <- tryCatch(terra::mask(r_idw, land_mask_rast), error = function(e) r_idw)
      }
      names(r_idw) <- vn
      loaded_layers[[vn]] <- r_idw
    }
  }

  ord_layers <- loaded_layers[kept_vars]
  terra::rast(ord_layers)
}

fit_and_project_maxent_sdm_maps <- function(
  df,
  srf_res,
  env_dir,
  out_dir = "outputs",
  kept_vars = NULL,
  primary_clim_vars = NULL,
  primary_edaph_vars = NULL,
  focal_taxon = "Packera dubia",
  intraspecific_mode = FALSE,
  demo = FALSE
) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  if (is.null(out_dir) || !nzchar(out_dir)) out_dir <- "outputs"
  message("=== Fitting MaxEnt (cloglog) & Spatial RF Ensemble SDMs & Projecting Suitability Rasters ===")

  # Resolve non-MEM environmental predictors
  if (is.null(kept_vars) || length(kept_vars) == 0L) {
    if (!is.null(srf_res) && !is.null(srf_res$collinearity_screening$kept_vars)) {
      kept_vars <- srf_res$collinearity_screening$kept_vars
    } else {
      imp_csv <- file.path(out_dir, "reports", "srf_feature_importance.csv")
      if (file.exists(imp_csv)) {
        imp_tbl <- tryCatch(utils::read.csv(imp_csv, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(imp_tbl) && "Predictor" %in% names(imp_tbl)) {
          kept_vars <- as.character(imp_tbl$Predictor)
        }
      }
      if (is.null(kept_vars) || length(kept_vars) == 0L) {
        cand_eco <- discover_candidate_eco_predictors(df, include_pheno = FALSE, min_finite = 2L)
        if (length(cand_eco) > 12L && exists("screen_collinear_predictors_by_block", mode = "function")) {
          kept_vars <- tryCatch(
            screen_collinear_predictors_by_block(df[, cand_eco, drop = FALSE], cor.threshold = 0.70, vif.threshold = 5.0)$kept_vars,
            error = function(e) NULL
          )
        }
        if (is.null(kept_vars) || length(kept_vars) == 0L) {
          clim_c <- cand_eco[grepl("^bio[0-9]+", cand_eco)]
          soil_c <- setdiff(cand_eco, clim_c)
          kept_vars <- c(utils::head(clim_c, 6L), utils::head(soil_c, 6L))
        }
      }
    }
  }
  kept_vars <- intersect(kept_vars[!grepl("^MEM", kept_vars)], names(df))
  if (length(kept_vars) > 14L) {
    kv_clim <- kept_vars[grepl("^bio[0-9]+", kept_vars)]
    kv_edaph <- setdiff(kept_vars, kv_clim)
    kept_vars <- c(utils::head(kv_clim, 7L), utils::head(kv_edaph, 7L))
  }
  if (length(kept_vars) == 0L) {
    kept_vars <- intersect(c("bio1_temp_mean", "bio12_precip_annual", "fine_soil_ph", "soil_ph", "fine_soil_sand", "soil_sand"), names(df))
  }
  if (length(kept_vars) == 0L) return(invisible(NULL))

  # Ensure numeric finite predictor values in df
  for (kv in kept_vars) {
    v <- suppressWarnings(as.numeric(df[[kv]]))
    med_v <- stats::median(v[is.finite(v)], na.rm = TRUE)
    if (!is.finite(med_v)) med_v <- 0.0
    v[!is.finite(v)] <- med_v
    df[[kv]] <- v
  }

  # Identify top climatic and edaphic predictors
  imp_df <- if (!is.null(srf_res) && is.data.frame(srf_res$importance)) srf_res$importance else NULL
  edaph_rx <- "^(fine_|ssurgo_|polaris_|solus_|topo_|soil_|bdricm|resdept_r|cfvo|awc|alr_|sgmc_|is_azonal|azonal_|ph_q|min_bedrock|max_cfvo)"
  clim_cand <- if (!is.null(primary_clim_vars) && length(intersect(primary_clim_vars, kept_vars)) > 0L) {
    intersect(primary_clim_vars, kept_vars)
  } else {
    kept_vars[!grepl(edaph_rx, kept_vars)]
  }
  edaph_cand <- if (!is.null(primary_edaph_vars) && length(intersect(primary_edaph_vars, kept_vars)) > 0L) {
    intersect(primary_edaph_vars, kept_vars)
  } else {
    kept_vars[grepl(edaph_rx, kept_vars)]
  }
  if (!is.null(imp_df) && all(c("Predictor", "MeanDecreaseAccuracy") %in% names(imp_df))) {
    imp_map <- stats::setNames(as.numeric(imp_df$MeanDecreaseAccuracy), as.character(imp_df$Predictor))
    if (length(clim_cand) > 1L) clim_cand <- clim_cand[order(ifelse(clim_cand %in% names(imp_map), imp_map[clim_cand], -Inf), decreasing = TRUE)]
    if (length(edaph_cand) > 1L) edaph_cand <- edaph_cand[order(ifelse(edaph_cand %in% names(imp_map), imp_map[edaph_cand], -Inf), decreasing = TRUE)]
  }
  top_clim_pred <- if (length(clim_cand) > 0L) clim_cand[1L] else kept_vars[1L]
  top_edaph_pred <- if (length(edaph_cand) > 0L) edaph_cand[1L] else kept_vars[min(2L, length(kept_vars))]

  # Resolve authority prior weights (Combined_Prior_Weight / authority_weight, Nomenclatural Types pinned at 1.00 on [0,1])
  raw_auth_w <- if ("Combined_Prior_Weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$Combined_Prior_Weight))
  } else if ("authority_weight" %in% names(df)) {
    suppressWarnings(as.numeric(df$authority_weight))
  } else {
    resolve_step06_authority_weights(df)
  }
  raw_auth_w[!is.finite(raw_auth_w) | raw_auth_w <= 0] <- 0.5
  w_max_curr <- max(raw_auth_w, 1e-6)
  auth_w_01 <- pmin(1.0, pmax(0.20, raw_auth_w / w_max_curr))
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) {
    .identify_type_specimen_vec(df)
  } else {
    is_t <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, nrow(df))
    is_t[is.na(is_t)] <- FALSE
    list(is_type_specimen = is_t, type_rank = ifelse(is_t, 2L, 0L))
  }
  auth_w_01[type_info$is_type_specimen] <- 1.00
  df$sdm_prior_weight_01 <- auth_w_01

  # Determine training-eligible coordinates (coordinateUncertaintyInMeters <= 10000 & finite lat/lon)
  lon_all <- suppressWarnings(as.numeric(df$longitude))
  lat_all <- suppressWarnings(as.numeric(df$latitude))
  unc_col <- if ("coordinateUncertaintyInMeters" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertaintyInMeters))
  } else if ("coordinateUncertainty" %in% names(df)) {
    suppressWarnings(as.numeric(df$coordinateUncertainty))
  } else {
    rep(NA_real_, nrow(df))
  }
  elig_mask <- is.finite(lon_all) & is.finite(lat_all) & (is.na(unc_col) | unc_col <= 10000.0 | type_info$is_type_specimen)
  if (sum(elig_mask) < 4L) elig_mask <- is.finite(lon_all) & is.finite(lat_all)

  # Build regional environmental SpatRaster stack across the study bounding box
  env_stack <- .build_regional_sdm_raster_stack(df[elig_mask, , drop = FALSE], env_dir = env_dir, kept_vars = kept_vars)

  # Generate supplemental random background cells from env_stack so target-group + random background >= 1,000 points
  set.seed(42)
  bg_pool_df <- NULL
  if (!is.null(env_stack) && requireNamespace("terra", quietly = TRUE)) {
    n_cells_tot <- terra::ncell(env_stack)
    s_idx <- sample.int(n_cells_tot, size = min(1200L, n_cells_tot), replace = FALSE)
    s_xy <- terra::xyFromCell(env_stack, s_idx)
    s_vals <- as.data.frame(terra::values(env_stack)[s_idx, , drop = FALSE])
    colnames(s_vals) <- names(env_stack)
    bg_pool_df <- data.frame(
      longitude = as.numeric(s_xy[, 1]),
      latitude = as.numeric(s_xy[, 2]),
      s_vals[, kept_vars, drop = FALSE],
      sdm_prior_weight_01 = 1.0,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    ok_bg <- stats::complete.cases(bg_pool_df[, c("longitude", "latitude", kept_vars), drop = FALSE])
    bg_pool_df <- bg_pool_df[ok_bg, , drop = FALSE]
  }

  # Determine target taxa / morphotypes to model
  is_intra <- isTRUE(intraspecific_mode) || (length(unique(stats::na.omit(df$species_standardized))) < 2L)
  grp_vec <- if (is_intra && "analysis_group" %in% names(df)) {
    as.character(df$analysis_group)
  } else {
    as.character(df$species_standardized)
  }
  grp_vec[is.na(grp_vec) | !nzchar(grp_vec)] <- "Unknown"
  tbl_elig <- table(grp_vec[elig_mask & grp_vec != "Unknown"])
  target_groups <- names(tbl_elig[tbl_elig >= 5L])
  if (length(target_groups) == 0L) target_groups <- names(tbl_elig[tbl_elig >= 2L])
  if (length(target_groups) == 0L && any(grp_vec != "Unknown")) target_groups <- unique(grp_vec[grp_vec != "Unknown"])

  # Also include focal_taxon explicitly if present in species_standardized with >= 3 records and not already in target_groups
  if (!is.null(focal_taxon) && nzchar(focal_taxon) && !(focal_taxon %in% target_groups)) {
    if (sum(df$species_standardized == focal_taxon & elig_mask, na.rm = TRUE) >= 3L) {
      target_groups <- c(focal_taxon, target_groups)
    }
  }
  # Order focal_taxon first
  if (!is.null(focal_taxon) && focal_taxon %in% target_groups) {
    target_groups <- c(focal_taxon, setdiff(target_groups, focal_taxon))
  }

  block_deg_use <- if (!is.null(srf_res) && is.finite(srf_res$spatial_autocorr_range_km)) {
    max(0.25, min(2.5, as.numeric(srf_res$spatial_autocorr_range_km) / 111.0))
  } else {
    0.50
  }

  raster_dir <- file.path(out_dir, "rasters")
  write_prod_rasters <- !(isTRUE(demo) && is_canonical_prod_path(raster_dir))
  if (write_prod_rasters) dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
  gdal_deflate <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6")

  perf_rows <- list()
  focal_artifacts <- NULL
  all_sdm_artifacts <- list()

  for (tg in target_groups) {
    pres_mask <- if (tg %in% df$species_standardized && !(tg %in% grp_vec)) {
      elig_mask & (df$species_standardized == tg)
    } else {
      elig_mask & (grp_vec == tg)
    }
    pres_df <- df[pres_mask, c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
    n_pres <- nrow(pres_df)
    if (n_pres < 2L) next

    # Target-group background points (Phillips et al. 2009): all other Packera localities + supplemental raster cells up to >= 1,000
    tg_bg_mask <- elig_mask & !pres_mask
    tg_bg_df <- df[tg_bg_mask, c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
    if (nrow(tg_bg_df) > 0L) tg_bg_df$sdm_prior_weight_01 <- 1.0
    n_need_supp <- max(0L, 1000L - nrow(tg_bg_df))
    if (n_need_supp > 0L && !is.null(bg_pool_df) && nrow(bg_pool_df) > 0L) {
      take_n <- min(n_need_supp, nrow(bg_pool_df))
      supp_df <- bg_pool_df[seq_len(take_n), c("longitude", "latitude", kept_vars, "sdm_prior_weight_01"), drop = FALSE]
      bg_df <- rbind(tg_bg_df, supp_df)
    } else if (nrow(tg_bg_df) > 0L) {
      bg_df <- tg_bg_df
    } else if (!is.null(bg_pool_df) && nrow(bg_pool_df) > 0L) {
      bg_df <- bg_pool_df
    } else {
      next
    }
    n_bg <- nrow(bg_df)

    comb_df <- rbind(pres_df, bg_df)
    p_vec <- c(rep(1L, n_pres), rep(0L, n_bg))
    w_vec <- as.numeric(comb_df$sdm_prior_weight_01)
    X_comb <- comb_df[, kept_vars, drop = FALSE]
    coords_comb <- as.matrix(comb_df[, c("longitude", "latitude")])

    # Spatial-Block Cross-Validation (k = 5 spatial blocks)
    k_folds <- min(5L, max(2L, n_pres))
    block_ids <- assign_spatial_grid_blocks(coords_comb, block_deg = block_deg_use, k = k_folds)
    # Ensure presence points span at least 2 folds when possible
    if (length(unique(block_ids[p_vec == 1L])) < 2L && n_pres >= 2L) {
      block_ids[p_vec == 1L] <- ((seq_len(n_pres) - 1L) %% min(k_folds, n_pres)) + 1L
    }

    oof_mx <- rep(NA_real_, length(p_vec))
    oof_rf <- rep(NA_real_, length(p_vec))
    u_folds <- sort(unique(block_ids))
    for (fid in u_folds) {
      tr_idx <- which(block_ids != fid)
      te_idx <- which(block_ids == fid)
      if (sum(p_vec[tr_idx] == 1L) >= 2L && sum(p_vec[tr_idx] == 0L) >= 2L && length(te_idx) >= 1L) {
        mx_cv <- .fit_single_maxent_model(p_vec[tr_idx], X_comb[tr_idx, , drop = FALSE], w_vec[tr_idx], use_java_maxent = FALSE)
        rf_cv <- .fit_single_srf_binary_model(p_vec[tr_idx], X_comb[tr_idx, , drop = FALSE], w_vec[tr_idx], n_trees = 200L)
        oof_mx[te_idx] <- mx_cv$predict_cloglog(X_comb[te_idx, , drop = FALSE])
        oof_rf[te_idx] <- rf_cv$predict_prob(X_comb[te_idx, , drop = FALSE])
      }
    }

    # Fit full models on all presence + target-group background data
    mx_full <- .fit_single_maxent_model(p_vec, X_comb, w_vec)
    rf_full <- .fit_single_srf_binary_model(p_vec, X_comb, w_vec, n_trees = 300L)

    # Fill any fold that lacked training presences with full model predictions
    na_mx <- !is.finite(oof_mx)
    if (any(na_mx)) oof_mx[na_mx] <- mx_full$predict_cloglog(X_comb[na_mx, , drop = FALSE])
    na_rf <- !is.finite(oof_rf)
    if (any(na_rf)) oof_rf[na_rf] <- rf_full$predict_prob(X_comb[na_rf, , drop = FALSE])

    roc_mx <- .compute_roc_auc_curve(p_vec, oof_mx, w_vec)
    roc_rf <- .compute_roc_auc_curve(p_vec, oof_rf, w_vec)
    auc_mx <- if (is.finite(roc_mx$auc)) roc_mx$auc else 0.5
    auc_rf <- if (is.finite(roc_rf$auc)) roc_rf$auc else 0.5

    # AUC_block-weighted ensemble (default 0.5 * MaxEnt + 0.5 * SpatialRF weighted by AUC_block)
    w_ens_mx <- if ((auc_mx + auc_rf) > 1e-6) auc_mx / (auc_mx + auc_rf) else 0.5
    w_ens_rf <- 1.0 - w_ens_mx
    oof_ens <- pmin(1.0, pmax(0.0, w_ens_mx * oof_mx + w_ens_rf * oof_rf))
    roc_ens <- .compute_roc_auc_curve(p_vec, oof_ens, w_vec)
    auc_ens <- if (is.finite(roc_ens$auc)) roc_ens$auc else round(0.5 * (auc_mx + auc_rf), 4)

    cbi_mx <- compute_continuous_boyce_index(oof_mx, oof_mx[p_vec == 1L], n_bins = 20L)
    cbi_rf <- compute_continuous_boyce_index(oof_rf, oof_rf[p_vec == 1L], n_bins = 20L)
    cbi_ens <- compute_continuous_boyce_index(oof_ens, oof_ens[p_vec == 1L], n_bins = 20L)

    perf_rows[[length(perf_rows) + 1L]] <- data.frame(
      Taxon_or_Morphotype = tg,
      N_Presences = as.integer(n_pres),
      N_Background = as.integer(n_bg),
      MaxEnt_Engine = as.character(mx_full$engine),
      AUC_Block_MaxEnt = round(auc_mx, 4),
      CBI_MaxEnt = round(as.numeric(cbi_mx), 4),
      AUC_Block_SpatialRF = round(auc_rf, 4),
      CBI_SpatialRF = round(as.numeric(cbi_rf), 4),
      AUC_Block_Ensemble = round(auc_ens, 4),
      CBI_Ensemble = round(as.numeric(cbi_ens), 4),
      Top_Climatic_Predictor = as.character(top_clim_pred),
      Top_Edaphic_Predictor = as.character(top_edaph_pred),
      stringsAsFactors = FALSE
    )

    # Project continuous habitat suitability across env_stack via terra::predict()
    rast_mx <- NULL
    rast_rf <- NULL
    rast_ens <- NULL
    if (!is.null(env_stack) && requireNamespace("terra", quietly = TRUE)) {
      rast_mx <- tryCatch(
        terra::predict(env_stack, mx_full, fun = function(model, data, ...) model$predict_cloglog(data), na.rm = TRUE),
        error = function(e) NULL
      )
      rast_rf <- tryCatch(
        terra::predict(env_stack, rf_full, fun = function(model, data, ...) model$predict_prob(data), na.rm = TRUE),
        error = function(e) NULL
      )
      if (!is.null(rast_mx) && !is.null(rast_rf)) {
        names(rast_mx) <- "maxent_cloglog"
        names(rast_rf) <- "spatial_rf_suitability"
        rast_ens <- w_ens_mx * rast_mx + w_ens_rf * rast_rf
        names(rast_ens) <- "ensemble_suitability"

        safe_name <- gsub("^_+|_+$", "", tolower(gsub("[^A-Za-z0-9]+", "_", tg)))
        mx_tif_out <- file.path(raster_dir, sprintf("sdm_%s_maxent_cloglog.tif", safe_name))
        ens_tif_out <- file.path(raster_dir, sprintf("sdm_%s_ensemble.tif", safe_name))

        .write_deflate_tif <- function(r_obj, path_target) {
          tmp_t <- paste0(path_target, ".tmp.tif")
          dir.create(dirname(path_target), recursive = TRUE, showWarnings = FALSE)
          tryCatch(
            terra::writeRaster(r_obj, filename = tmp_t, overwrite = TRUE, gdal = gdal_deflate),
            error = function(e) terra::writeRaster(r_obj, filename = tmp_t, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "ZLEVEL=6"))
          )
          if (file.exists(tmp_t)) file.rename(tmp_t, path_target)
        }

        if (write_prod_rasters) {
          tryCatch(.write_deflate_tif(rast_mx, mx_tif_out), error = function(e) NULL)
          tryCatch(.write_deflate_tif(rast_ens, ens_tif_out), error = function(e) NULL)
        }
      }
    }

    tg_art <- list(
      taxon = tg,
      pres_mask = pres_mask,
      mx_full = mx_full,
      rf_full = rf_full,
      rast_mx = rast_mx,
      rast_rf = rast_rf,
      rast_ens = rast_ens,
      roc_mx = roc_mx,
      roc_rf = roc_rf,
      roc_ens = roc_ens,
      cbi_mx = cbi_mx,
      cbi_rf = cbi_rf,
      cbi_ens = cbi_ens,
      X_comb = X_comb,
      p_vec = p_vec
    )
    all_sdm_artifacts[[tg]] <- tg_art
    if (is.null(focal_artifacts)) {
      focal_artifacts <- tg_art
    }
  }

  perf_df <- if (length(perf_rows) > 0L) {
    do.call(rbind, perf_rows)
  } else {
    data.frame(
      Taxon_or_Morphotype = character(0), N_Presences = integer(0), N_Background = integer(0),
      MaxEnt_Engine = character(0), AUC_Block_MaxEnt = numeric(0), CBI_MaxEnt = numeric(0),
      AUC_Block_SpatialRF = numeric(0), CBI_SpatialRF = numeric(0), AUC_Block_Ensemble = numeric(0),
      CBI_Ensemble = numeric(0), Top_Climatic_Predictor = character(0), Top_Edaphic_Predictor = character(0),
      stringsAsFactors = FALSE
    )
  }

  perf_csv_out <- file.path(out_dir, "reports", "sdm_model_performance.csv")
  if (write_csv_atomic(perf_df, perf_csv_out, demo = demo)) {
    message(sprintf("Saved SDM performance metrics (%d modeled taxa/morphotypes) to %s.", nrow(perf_df), perf_csv_out))
  }

  # Render 4-panel publication figure fig06b_sdm_habitat_suitability_map.{pdf,png}
  if (!is.null(focal_artifacts) && !is.null(focal_artifacts$rast_mx) && !is.null(focal_artifacts$rast_ens)) {
    fig_dir <- file.path(out_dir, "figures")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_out <- file.path(fig_dir, "fig06b_sdm_habitat_suitability_map.pdf")
    png_out <- file.path(fig_dir, "fig06b_sdm_habitat_suitability_map.png")

    # Select top 3 climatic + edaphic predictors for Panel C marginal response curves
    resp_vars <- unique(c(top_clim_pred, top_edaph_pred, clim_cand, edaph_cand, kept_vars))
    resp_vars <- resp_vars[seq_len(min(3L, length(resp_vars)))]
    resp_lbls <- format_step06_predictor_labels(resp_vars)

    # Optional EPA Level III Ecoregion shapefile overlay
    epa_vect <- NULL
    if (requireNamespace("terra", quietly = TRUE)) {
      epa_shp <- list.files("data/environmental/epa_ecoregions_l3", pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
      if (length(epa_shp) > 0L) {
        epa_vect <- tryCatch({
          v_e <- terra::vect(epa_shp[1L])
          v_e <- terra::project(v_e, "EPSG:4326")
          terra::crop(v_e, terra::ext(focal_artifacts$rast_mx))
        }, error = function(e) NULL)
      }
    }

    render_sdm_plate <- function() {
      graphics::par(mfrow = c(2, 2), mar = c(4.4, 4.6, 3.0, 1.5), oma = c(0.5, 0.5, 1.8, 0.5))
      suit_ramp <- grDevices::colorRampPalette(c("#F7FCF0", "#E0F3DB", "#A8DDB5", "#4EB3D3", "#0868AC", "#084081"))(100)

      .plot_suit_rast <- function(r_layer, panel_title, sub_txt) {
        ex_r <- terra::ext(r_layer)
        xs <- terra::xFromCol(r_layer, seq_len(terra::ncol(r_layer)))
        ys <- rev(terra::yFromRow(r_layer, seq_len(terra::nrow(r_layer))))
        mat_v <- terra::as.matrix(r_layer, wide = TRUE)
        mat_plot <- t(mat_v[rev(seq_len(nrow(mat_v))), , drop = FALSE])
        graphics::image(
          x = xs, y = ys, z = mat_plot, zlim = c(0, 1), col = suit_ramp,
          xlab = "Decimal Longitude (deg E)", ylab = "Decimal Latitude (deg N)",
          main = panel_title, font.main = 2
        )
        graphics::mtext(sub_txt, side = 3, line = 0.25, cex = 0.72, col = "grey25")
        if (!is.null(epa_vect)) {
          tryCatch(terra::lines(epa_vect, col = grDevices::adjustcolor("#555555", alpha.f = 0.35), lwd = 0.55), error = function(e) NULL)
        }
        if (requireNamespace("maps", quietly = TRUE)) {
          tryCatch(
            maps::map("state", xlim = c(terra::xmin(ex_r), terra::xmax(ex_r)),
                      ylim = c(terra::ymin(ex_r), terra::ymax(ex_r)),
                      add = TRUE, col = grDevices::adjustcolor("#222222", alpha.f = 0.50), lwd = 0.75),
            error = function(e) NULL
          )
        }
        tryCatch(
          graphics::contour(x = xs, y = ys, z = mat_plot, levels = c(0.50, 0.75),
                            add = TRUE, col = c("#E69F00", "#D55E00"), lwd = c(0.9, 1.2), drawlabels = TRUE),
          error = function(e) NULL
        )
        graphics::box(col = "grey40")
      }

      # Panel A: MaxEnt (cloglog) Continuous Habitat Suitability Map + Authority Tiers + Type Localities
      p_row_1 <- perf_df[perf_df$Taxon_or_Morphotype == focal_artifacts$taxon, , drop = FALSE]
      sub_a <- sprintf("%s | Block AUC = %.3f | Continuous Boyce Index = %.3f",
                       focal_artifacts$mx_full$engine, p_row_1$AUC_Block_MaxEnt[1], p_row_1$CBI_MaxEnt[1])
      .plot_suit_rast(
        focal_artifacts$rast_mx,
        sprintf("A. MaxEnt (cloglog) Habitat Suitability: %s", focal_artifacts$taxon),
        sub_a
      )
      pm <- which(focal_artifacts$pres_mask)
      det_t <- if ("determiner_tier" %in% names(df)) as.character(df$determiner_tier[pm]) else rep("Tier_1_Gold", length(pm))
      pch_t <- ifelse(grepl("Gold|Tier_1", det_t, ignore.case = TRUE), 1L,
               ifelse(grepl("Silver|Tier_2", det_t, ignore.case = TRUE), 2L, 0L))
      col_t <- ifelse(grepl("Gold|Tier_1", det_t, ignore.case = TRUE), "#D55E00",
               ifelse(grepl("Silver|Tier_2", det_t, ignore.case = TRUE), "#E69F00", "#555555"))
      graphics::points(lon_all[pm], lat_all[pm], pch = pch_t, col = col_t, bg = NA, cex = 0.90, lwd = 1.1)
      t_idx <- pm[which(type_info$is_type_specimen[pm])]
      if (length(t_idx) > 0L) {
        graphics::points(lon_all[t_idx], lat_all[t_idx], pch = 5L, col = "#D55E00", bg = NA, cex = 1.55, lwd = 2.2)
        graphics::points(lon_all[t_idx], lat_all[t_idx], pch = 1L, col = "#111111", bg = NA, cex = 1.95, lwd = 1.5)
      }
      graphics::legend(
        "bottomright",
        legend = c("Tier 1 Gold", "Tier 2 Silver", "Tier 3 Bronze", "Nomenclatural Type", "Suitability >= 0.75"),
        pch = c(1L, 2L, 0L, 5L, NA), lty = c(NA, NA, NA, NA, 1L),
        col = c("#D55E00", "#E69F00", "#555555", "#111111", "#D55E00"),
        pt.lwd = c(1.2, 1.2, 1.0, 2.0, NA), lwd = c(NA, NA, NA, NA, 1.3),
        cex = 0.68, bg = grDevices::adjustcolor("white", alpha.f = 0.85)
      )

      # Panel B: Spatial RF + MaxEnt Ensemble Suitability Map
      sub_b <- sprintf("AUC-Weighted Ensemble (MaxEnt + Spatial RF) | Block AUC = %.3f | CBI = %.3f",
                       p_row_1$AUC_Block_Ensemble[1], p_row_1$CBI_Ensemble[1])
      .plot_suit_rast(
        focal_artifacts$rast_ens,
        sprintf("B. Ensemble Habitat Suitability (MaxEnt + Spatial RF): %s", focal_artifacts$taxon),
        sub_b
      )
      oth_mask <- which(elig_mask & !focal_artifacts$pres_mask)
      if (length(oth_mask) > 0L) {
        graphics::points(lon_all[oth_mask], lat_all[oth_mask], pch = 4L,
                         col = grDevices::adjustcolor("#444444", alpha.f = 0.45), cex = 0.65, lwd = 0.7)
      }
      graphics::points(lon_all[pm], lat_all[pm], pch = 1L, col = "#D55E00", bg = NA, cex = 0.85, lwd = 1.1)
      if (length(t_idx) > 0L) {
        graphics::points(lon_all[t_idx], lat_all[t_idx], pch = 5L, col = "#111111", bg = NA, cex = 1.60, lwd = 2.0)
      }
      graphics::legend(
        "bottomright",
        legend = c(sprintf("%s Presences", gsub("^Packera\\s+", "P. ", focal_artifacts$taxon)),
                   "Target-Group Background", "0.50 / 0.75 Iso-Suitability"),
        pch = c(1L, 4L, NA), lty = c(NA, NA, 1L),
        col = c("#D55E00", "#444444", "#E69F00"),
        cex = 0.68, bg = grDevices::adjustcolor("white", alpha.f = 0.85)
      )

      # Panel C: Marginal MaxEnt (cloglog) Response Curves for Top 3 Climatic & Edaphic Predictors
      graphics::plot(
        c(0, 1), c(0, 1), type = "n",
        xlab = "Normalized Predictor Gradient (5th to 95th Percentile, Others at Median)",
        ylab = "Marginal MaxEnt Suitability (cloglog)",
        main = "C. Marginal MaxEnt (cloglog) Ecological Response Curves", font.main = 2
      )
      graphics::grid(col = "grey88")
      med_row <- as.data.frame(lapply(focal_artifacts$X_comb, function(col_x) stats::median(as.numeric(col_x), na.rm = TRUE)), check.names = FALSE)
      resp_cols <- c("#0072B2", "#009E73", "#D55E00")
      resp_ltys <- c(1L, 2L, 4L)
      leg_c <- character(length(resp_vars))
      u_grid <- seq(0, 1, length.out = 100L)
      for (vi in seq_along(resp_vars)) {
        vn <- resp_vars[vi]
        v_raw <- as.numeric(focal_artifacts$X_comb[[vn]])
        q_lo <- as.numeric(stats::quantile(v_raw, probs = 0.05, na.rm = TRUE))
        q_hi <- as.numeric(stats::quantile(v_raw, probs = 0.95, na.rm = TRUE))
        if (!is.finite(q_lo) || !is.finite(q_hi) || abs(q_hi - q_lo) < 1e-8) {
          q_lo <- min(v_raw, na.rm = TRUE); q_hi <- max(v_raw, na.rm = TRUE)
        }
        eval_df <- med_row[rep(1L, 100L), , drop = FALSE]
        eval_df[[vn]] <- q_lo + u_grid * (q_hi - q_lo)
        y_resp <- focal_artifacts$mx_full$predict_cloglog(eval_df)
        graphics::lines(u_grid, y_resp, col = resp_cols[vi], lty = resp_ltys[vi], lwd = 2.2)
        leg_c[vi] <- sprintf("%s [%.1f to %.1f]", resp_lbls[vi], q_lo, q_hi)
      }
      graphics::legend(
        "topright", legend = leg_c, col = resp_cols[seq_along(resp_vars)],
        lty = resp_ltys[seq_along(resp_vars)], lwd = 2.2, cex = 0.68,
        bg = grDevices::adjustcolor("white", alpha.f = 0.85)
      )

      # Panel D: Spatial-Block Cross-Validated ROC/AUC & Continuous Boyce Index (P/E Ratio)
      graphics::plot(
        c(0, 1), c(0, 1), type = "n",
        xlab = "False Positive Rate (1 - Specificity) / Suitability Bin",
        ylab = "True Positive Rate (Sensitivity) / Scaled P/E",
        main = "D. Spatial-Block CV ROC Curves & Boyce Calibration (P/E)", font.main = 2
      )
      graphics::grid(col = "grey88")
      graphics::abline(0, 1, lty = 3, col = "grey55", lwd = 1.1)
      graphics::lines(focal_artifacts$roc_mx$fpr, focal_artifacts$roc_mx$tpr, col = "#0072B2", lwd = 2.1, lty = 1L)
      graphics::lines(focal_artifacts$roc_rf$fpr, focal_artifacts$roc_rf$tpr, col = "#009E73", lwd = 2.0, lty = 2L)
      graphics::lines(focal_artifacts$roc_ens$fpr, focal_artifacts$roc_ens$tpr, col = "#D55E00", lwd = 2.3, lty = 1L)

      b_curve <- attr(focal_artifacts$cbi_ens, "boyce_curve")
      if (is.data.frame(b_curve) && nrow(b_curve) >= 2L && any(is.finite(b_curve$pe_ratio))) {
        pe_max <- max(b_curve$pe_ratio, 1.0, na.rm = TRUE)
        pe_scaled <- pmin(1.0, pmax(0.0, b_curve$pe_ratio / pe_max))
        graphics::lines(b_curve$suitability, pe_scaled, col = "#CC79A7", lwd = 1.8, lty = 4L, type = "b", pch = 16, cex = 0.65)
      }
      graphics::legend(
        "bottomright",
        legend = c(
          sprintf("MaxEnt (AUC = %.3f, CBI = %.2f)", p_row_1$AUC_Block_MaxEnt[1], p_row_1$CBI_MaxEnt[1]),
          sprintf("Spatial RF (AUC = %.3f, CBI = %.2f)", p_row_1$AUC_Block_SpatialRF[1], p_row_1$CBI_SpatialRF[1]),
          sprintf("Ensemble (AUC = %.3f, CBI = %.2f)", p_row_1$AUC_Block_Ensemble[1], p_row_1$CBI_Ensemble[1]),
          "Boyce P/E Ratio Curve (Scaled)"
        ),
        col = c("#0072B2", "#009E73", "#D55E00", "#CC79A7"),
        lty = c(1L, 2L, 1L, 4L), lwd = c(2.1, 2.0, 2.3, 1.8),
        cex = 0.68, bg = grDevices::adjustcolor("white", alpha.f = 0.85)
      )
      graphics::title(
        sprintf("Step 06 Species Distribution Modeling (MaxEnt & Spatial RF): %s", focal_artifacts$taxon),
        outer = TRUE, cex.main = 1.18
      )
    }

    if (!(isTRUE(demo) && is_canonical_prod_path(fig_dir))) {
      pdf_tmp <- paste0(pdf_out, ".tmp")
      grDevices::pdf(pdf_tmp, width = 13.0, height = 9.5)
      tryCatch(render_sdm_plate(), finally = grDevices::dev.off())
      if (file.exists(pdf_tmp)) file.rename(pdf_tmp, pdf_out)

      png_tmp <- paste0(png_out, ".tmp")
      grDevices::png(png_tmp, width = 13.0, height = 9.5, units = "in", res = 300)
      tryCatch(render_sdm_plate(), finally = grDevices::dev.off())
      if (file.exists(png_tmp)) file.rename(png_tmp, png_out)
      message(sprintf("Exported 4-panel MaxEnt & Ensemble SDM figure to %s (.pdf & .png).", pdf_out))

      # When multiple taxa/morphotypes are modeled, also render a multi-species MaxEnt SDM comparison grid
      valid_multi_arts <- Filter(function(a) !is.null(a$rast_mx), all_sdm_artifacts)
      if (length(valid_multi_arts) > 1L) {
        grid_pdf_out <- file.path(fig_dir, "fig06c_multispecies_maxent_sdm_grid.pdf")
        grid_png_out <- file.path(fig_dir, "fig06c_multispecies_maxent_sdm_grid.png")
        n_sp_m <- length(valid_multi_arts)
        n_cols_m <- if (n_sp_m <= 4L) 2L else if (n_sp_m <= 9L) 3L else 4L
        n_rows_m <- as.integer(ceiling(n_sp_m / n_cols_m))
        fig_w_m <- max(10.0, 3.6 * n_cols_m)
        fig_h_m <- max(7.0, 3.1 * n_rows_m + 0.8)

        render_multispecies_sdm_grid <- function() {
          graphics::par(mfrow = c(n_rows_m, n_cols_m), mar = c(3.2, 3.4, 2.6, 1.0), oma = c(0.5, 0.5, 2.0, 0.5))
          suit_ramp <- grDevices::colorRampPalette(c("#F7FCF0", "#E0F3DB", "#A8DDB5", "#4EB3D3", "#0868AC", "#084081"))(100)
          for (art_i in valid_multi_arts) {
            r_layer <- art_i$rast_mx
            ex_r <- terra::ext(r_layer)
            xs <- terra::xFromCol(r_layer, seq_len(terra::ncol(r_layer)))
            ys <- rev(terra::yFromRow(r_layer, seq_len(terra::nrow(r_layer))))
            mat_v <- terra::as.matrix(r_layer, wide = TRUE)
            mat_plot <- t(mat_v[rev(seq_len(nrow(mat_v))), , drop = FALSE])
            p_row_i <- perf_df[perf_df$Taxon_or_Morphotype == art_i$taxon, , drop = FALSE]
            auc_v <- if (nrow(p_row_i) > 0L) p_row_i$AUC_Block_MaxEnt[1] else NA_real_
            cbi_v <- if (nrow(p_row_i) > 0L) p_row_i$CBI_MaxEnt[1] else NA_real_
            n_pr_v <- sum(art_i$pres_mask, na.rm = TRUE)
            graphics::image(
              x = xs, y = ys, z = mat_plot, zlim = c(0, 1), col = suit_ramp,
              xlab = "Lon (deg E)", ylab = "Lat (deg N)",
              main = sprintf("%s (n = %d)", art_i$taxon, n_pr_v), font.main = 4, cex.main = 0.95
            )
            graphics::mtext(
              sprintf("MaxEnt AUC = %.3f | CBI = %.2f", auc_v, cbi_v),
              side = 3, line = 0.15, cex = 0.65, col = "grey25"
            )
            if (requireNamespace("maps", quietly = TRUE)) {
              tryCatch(
                maps::map("state", xlim = c(terra::xmin(ex_r), terra::xmax(ex_r)),
                          ylim = c(terra::ymin(ex_r), terra::ymax(ex_r)),
                          add = TRUE, col = grDevices::adjustcolor("#222222", alpha.f = 0.45), lwd = 0.6),
                error = function(e) NULL
              )
            }
            tryCatch(
              graphics::contour(x = xs, y = ys, z = mat_plot, levels = c(0.50, 0.75),
                                add = TRUE, col = c("#E69F00", "#D55E00"), lwd = c(0.7, 1.0), drawlabels = FALSE),
              error = function(e) NULL
            )
            pm_i <- which(art_i$pres_mask)
            graphics::points(lon_all[pm_i], lat_all[pm_i], pch = 1L, col = "#D55E00", cex = 0.55, lwd = 0.8)
            graphics::box(col = "grey40")
          }
          graphics::title(
            sprintf("Step 06 Multi-Species MaxEnt (cloglog) Habitat Suitability Maps (%d Modeled Taxa)", n_sp_m),
            outer = TRUE, cex.main = 1.15
          )
        }

        g_pdf_tmp <- paste0(grid_pdf_out, ".tmp")
        grDevices::pdf(g_pdf_tmp, width = fig_w_m, height = fig_h_m)
        tryCatch(render_multispecies_sdm_grid(), finally = grDevices::dev.off())
        if (file.exists(g_pdf_tmp)) file.rename(g_pdf_tmp, grid_pdf_out)

        g_png_tmp <- paste0(grid_png_out, ".tmp")
        grDevices::png(g_png_tmp, width = fig_w_m, height = fig_h_m, units = "in", res = 300)
        tryCatch(render_multispecies_sdm_grid(), finally = grDevices::dev.off())
        if (file.exists(g_png_tmp)) file.rename(g_png_tmp, grid_png_out)
        message(sprintf("Exported multi-species MaxEnt SDM grid (%d taxa) to %s (.pdf & .png).", n_sp_m, grid_pdf_out))
      }
    }
  }

  invisible(list(performance = perf_df, focal = focal_artifacts))
}

