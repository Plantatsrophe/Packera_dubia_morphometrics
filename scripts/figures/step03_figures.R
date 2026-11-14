# ==============================================================================
# Module: scripts/figures/step03_figures.R
# Step 03 Canonical 300-DPI PNG Figures (Momocs EFA + LeafMachine3 Radial ECT):
#   1. fig03a_efa_morphospace_and_harmonics.png
#   2. fig03b_topological_ect_and_heteroblasty.png
# ==============================================================================

render_step03_figures_suite <- function(
  efa_df,
  specimen_df,
  harm_mat = NULL,
  pca_fit = NULL,
  figures_dir = "outputs/figures",
  output_csv_path = NULL,
  sym_pca_fit = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    warning("Step 03 figures skipped: ggplot2 and patchwork are required.")
    return(invisible(character(0)))
  }
  fdir <- if (exists("resolve_active_figures_dir", mode = "function")) {
    resolve_active_figures_dir(figures_dir, output_csv_path)
  } else {
    as.character(figures_dir)
  }
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

  # Remove legacy multi-file Step 03 PDF/PNG clutter
  for (old_stem in c(
    "fig03_efa_morphospace_reconstruction", "fig03a_efa_morphospace_pca",
    "fig03b_harmonic_power_and_asymmetry", "fig03c_cauline_efa_morphospace",
    "fig03d_basal_ect_morphospace", "fig03e_cauline_ect_morphospace"
  )) {
    for (ext in c(".pdf", ".png")) {
      p_old <- file.path(fdir, paste0(old_stem, ext))
      if (file.exists(p_old)) unlink(p_old, force = TRUE)
    }
  }

  spec_df <- as.data.frame(if (!is.null(specimen_df) && nrow(specimen_df) > 0L) specimen_df else efa_df, stringsAsFactors = FALSE)
  if (is.null(spec_df) || nrow(spec_df) == 0L) return(invisible(character(0)))

  # Trigger Step 02 native LM3 asset/ECT export if overlays do not yet exist in fdir
  if (!dir.exists(file.path(fdir, "ect_overlays")) || length(list.files(file.path(fdir, "ect_overlays"), pattern = "\\.png$")) == 0L) {
    py_bin <- Sys.which("python3"); if (!nzchar(py_bin)) py_bin <- Sys.which("python")
    if (nzchar(py_bin) && file.exists("scripts/vision/step02_figures.py")) {
      tryCatch(
        system2(py_bin, c("scripts/vision/step02_figures.py", "--export-ect-exemplars", "--output-dir", fdir, "--zone", "all"), stdout = FALSE, stderr = FALSE),
        error = function(e) NULL
      )
    }
  }

  spec_df$plot_group <- resolve_plot_group_labels(spec_df)
  active_groups <- sort(unique(spec_df$plot_group[!is.na(spec_df$plot_group) & nzchar(spec_df$plot_group) & spec_df$plot_group != "Unknown"]))
  if (length(active_groups) == 0L) active_groups <- unique(spec_df$plot_group)
  spec_df <- spec_df[spec_df$plot_group %in% active_groups, , drop = FALSE]

  pal_map <- get_combinatorial_palette_r(active_groups)
  shp_map <- get_combinatorial_hollow_shapes_r(active_groups)
  cnt_map <- table(spec_df$plot_group)
  lbl_exprs <- format_icn_label_vec(active_groups, abbreviate_genus = (length(active_groups) > 5L), counts_map = cnt_map)
  names(lbl_exprs) <- active_groups

  # ----------------------------------------------------------------------------
  # Plate 03A: fig03a_efa_morphospace_and_harmonics.png
  #   Panel A: Cumulative Fourier Harmonic Power (H_opt cutoff) & Asymmetry
  #   Panel B: 2D Symmetric EFA PCA (PC1 vs PC2) + Inverse-Fourier Background Grid
  #   Panel C: Dynamic PC1..PC{k_95} Variance Scree + Native LM3 Medoid Silhouettes
  # ----------------------------------------------------------------------------
  harm_cols <- grep("^[ABCD][0-9]+$", names(spec_df), value = TRUE)
  max_h <- if (length(harm_cols) > 0L) max(suppressWarnings(as.integer(gsub("^[ABCD]", "", harm_cols))), na.rm = TRUE) else 12L
  max_h <- max(6L, min(24L, ifelse(is.finite(max_h), max_h, 12L)))
  h_opt_val <- if ("efd_harmonic_opt" %in% names(spec_df) && any(is.finite(suppressWarnings(as.numeric(spec_df$efd_harmonic_opt))))) {
    as.integer(stats::median(suppressWarnings(as.numeric(spec_df$efd_harmonic_opt)), na.rm = TRUE))
  } else min(12L, max_h)

  # Compute empirical cumulative harmonic power across harmonics 2..max_h
  pow_vec <- vapply(seq_len(max_h), function(k) {
    ck <- intersect(paste0(c("A", "B", "C", "D"), k), names(spec_df))
    if (length(ck) == 0L || k == 1L) return(0.0)
    m_k <- as.matrix(spec_df[, ck, drop = FALSE]); storage.mode(m_k) <- "double"
    mean(rowSums(m_k^2, na.rm = TRUE), na.rm = TRUE)
  }, numeric(1))
  if (sum(pow_vec) <= 1e-12) pow_vec <- (1.0 / (seq_len(max_h)^2.4))
  cum_pow <- 100.0 * cumsum(pow_vec) / max(sum(pow_vec), 1e-12)
  hp_df <- data.frame(harmonic = seq_len(max_h), cum_power = cum_pow)

  asym_med <- NA_real_
  for (acol in c("bilateral_asymmetry", "asym_power_pct", "pct_asymmetry")) {
    for (src_df in list(efa_df, spec_df)) {
      if (!is.null(src_df) && is.data.frame(src_df) && (acol %in% names(src_df))) {
        av <- suppressWarnings(as.numeric(src_df[[acol]]))
        av <- av[is.finite(av)]
        if (length(av) > 0L) {
          m_av <- stats::median(av, na.rm = TRUE)
          if (m_av <= 1.0 && max(av, na.rm = TRUE) <= 1.0 && !grepl("power_pct", acol)) m_av <- 100.0 * m_av
          asym_med <- m_av
          break
        }
      }
    }
    if (is.finite(asym_med)) break
  }
  asym_sub <- if (is.finite(asym_med)) sprintf("Calibrated H_opt = %d (>=99%% power) | Median bilateral asymmetry = %.1f%%", h_opt_val, asym_med) else sprintf("Calibrated H_opt = %d (>=99%% Fourier power)", h_opt_val)

  p03a_power <- ggplot2::ggplot(hp_df, ggplot2::aes(x = .data$harmonic, y = .data$cum_power)) +
    ggplot2::geom_hline(yintercept = 99.0, linetype = "dashed", color = "#6B7280", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = h_opt_val, linetype = "dotted", color = "#D55E00", linewidth = 0.7) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.9) +
    ggplot2::geom_point(shape = 1, size = 2.4, stroke = 0.9, color = "#0072B2") +
    ggplot2::annotate("label", x = h_opt_val, y = 84, label = sprintf("H_opt = %d", h_opt_val), hjust = if (h_opt_val >= max_h - 1L) 1.05 else 0.5, size = 2.9, color = "#D55E00", fill = "white") +
    ggplot2::scale_x_continuous(limits = c(1, max_h + 1.3), breaks = pretty(c(1, max_h))) +
    ggplot2::scale_y_continuous(limits = c(min(50, min(cum_pow)), 101)) +
    theme_packera_pub() +
    ggplot2::labs(title = "A. Fourier Harmonic Power Calibration", subtitle = asym_sub, x = "Harmonic Order (n)", y = "Cumulative Shape Power (%)")

  # Panel B: 2D Symmetric EFA PCA with Background 5x3 Inverse-Fourier Deformation Grid
  pca_use <- if (!is.null(sym_pca_fit)) sym_pca_fit else pca_fit
  pc_cols <- grep("^PC[0-9]+$", names(spec_df), value = TRUE)
  pc_cols <- pc_cols[order(as.integer(sub("^PC", "", pc_cols)))]
  var_pcts <- if (!is.null(pca_use) && !is.null(pca_use$sdev)) {
    v <- (pca_use$sdev^2) / max(sum(pca_use$sdev^2), 1e-12) * 100.0
    v
  } else if (length(pc_cols) >= 2L) {
    vv <- vapply(pc_cols, function(cl) stats::var(suppressWarnings(as.numeric(spec_df[[cl]])), na.rm = TRUE), numeric(1))
    vv / max(sum(vv, na.rm = TRUE), 1e-12) * 100.0
  } else c(50.0, 25.0)

  spec_df$pc1_num <- if ("PC1" %in% names(spec_df)) suppressWarnings(as.numeric(spec_df$PC1)) else rep(0.0, nrow(spec_df))
  spec_df$pc2_num <- if ("PC2" %in% names(spec_df)) suppressWarnings(as.numeric(spec_df$PC2)) else rep(0.0, nrow(spec_df))
  ok_pc <- is.finite(spec_df$pc1_num) & is.finite(spec_df$pc2_num)
  pc_sub <- spec_df[ok_pc, , drop = FALSE]

  # Build background 5x3 inverse-Fourier outlines when rotation + center exist
  bg_poly_df <- data.frame()
  if (nrow(pc_sub) >= 3L && !is.null(pca_use) && !is.null(pca_use$rotation) && !is.null(pca_use$center) &&
      ncol(pca_use$rotation) >= 2L && exists("reconstruct_efa_outline", mode = "function")) {
    xr <- range(pc_sub$pc1_num, na.rm = TRUE); yr <- range(pc_sub$pc2_num, na.rm = TRUE)
    x_span <- max(diff(xr), 0.2); y_span <- max(diff(yr), 0.2)
    gx_seq <- seq(xr[1] + 0.08 * x_span, xr[2] - 0.08 * x_span, length.out = 5L)
    gy_seq <- seq(yr[1] + 0.12 * y_span, yr[2] - 0.12 * y_span, length.out = 3L)
    cell_scale <- 0.075 * min(x_span, y_span)
    poly_list <- list()
    gid <- 0L
    for (gx in gx_seq) {
      for (gy in gy_seq) {
        gid <- gid + 1L
        h_rec <- pca_use$center + pca_use$rotation[, 1L] * gx + pca_use$rotation[, 2L] * gy
        pts <- tryCatch(reconstruct_efa_outline(h_rec, nb_pts = 80L), error = function(e) NULL)
        if (!is.null(pts) && is.matrix(pts) && nrow(pts) >= 10L) {
          pts_c <- scale(pts, center = TRUE, scale = FALSE)
          max_r <- max(sqrt(rowSums(pts_c^2)), 1e-6)
          poly_list[[gid]] <- data.frame(
            poly_id = gid,
            x = gx + (pts_c[, 1L] / max_r) * cell_scale,
            y = gy + (pts_c[, 2L] / max_r) * cell_scale
          )
        }
      }
    }
    if (length(poly_list) > 0L) bg_poly_df <- do.call(rbind, poly_list)
  }

  p03a_pca <- ggplot2::ggplot(pc_sub, ggplot2::aes(x = .data$pc1_num, y = .data$pc2_num))
  if (nrow(bg_poly_df) > 0L) {
    p03a_pca <- p03a_pca + ggplot2::geom_polygon(
      data = bg_poly_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_id),
      fill = "#F3F4F6", color = "#CBD5E1", linewidth = 0.35, alpha = 0.75, inherit.aes = FALSE
    )
  }
  p03a_pca <- p03a_pca +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4)

  el_groups <- names(which(table(pc_sub$plot_group) >= 4L))
  if (length(el_groups) > 0L) {
    p03a_pca <- p03a_pca + ggplot2::stat_ellipse(
      data = pc_sub[pc_sub$plot_group %in% el_groups, , drop = FALSE],
      ggplot2::aes(color = .data$plot_group),
      level = 0.95, linewidth = 0.55, linetype = "dashed", show.legend = FALSE
    )
  }
  p03a_pca <- p03a_pca +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$plot_group, shape = .data$plot_group),
      size = 2.4, stroke = 0.85, alpha = 0.82
    ) +
    ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
    ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group") +
    theme_packera_pub() +
    ggplot2::labs(
      title = "B. Symmetric Elliptic Fourier Morphospace (Inverse-Fourier Deformation Grid)",
      subtitle = "Authority-anchored PCA (hollow markers; passive Bronze weight = 0.0; 95% confidence ellipses)",
      x = sprintf("Symmetric EFA PC1 (%.1f%%)", ifelse(length(var_pcts) >= 1L, var_pcts[1L], 0.0)),
      y = sprintf("Symmetric EFA PC2 (%.1f%%)", ifelse(length(var_pcts) >= 2L, var_pcts[2L], 0.0))
    )

  if ("is_type_specimen" %in% names(pc_sub)) {
    is_t <- !is.na(pc_sub$is_type_specimen) & as.logical(pc_sub$is_type_specimen)
    if (any(is_t)) {
      p03a_pca <- p03a_pca + ggplot2::geom_point(
        data = pc_sub[is_t, , drop = FALSE],
        ggplot2::aes(x = .data$pc1_num, y = .data$pc2_num, color = .data$plot_group),
        shape = 8, size = 3.6, stroke = 1.1, show.legend = FALSE
      )
    }
  }
  if (length(active_groups) > 10L && nrow(pc_sub) > 0L) {
    cent_pc <- stats::aggregate(cbind(pc1_num, pc2_num) ~ plot_group, data = pc_sub, FUN = stats::median)
    cent_pc$lbl <- abbreviate_icn_taxon(cent_pc$plot_group, abbreviate_genus = TRUE)
    p03a_pca <- add_centroid_repel_labels_gg(p03a_pca, cent_pc, "pc1_num", "pc2_num", "lbl", color_col = "plot_group") +
      ggplot2::theme(legend.position = "none")
  }

  # Panel C: Dynamic Retained PCs Scree (PC1..PC{k_export})
  k_show <- min(length(var_pcts), max(5L, if ("efd_k95_pcs" %in% names(spec_df) && any(is.finite(suppressWarnings(as.numeric(spec_df$efd_k95_pcs))))) as.integer(max(suppressWarnings(as.numeric(spec_df$efd_k95_pcs)), na.rm = TRUE)) else length(pc_cols)))
  k_show <- max(2L, min(15L, k_show))
  scree_df <- data.frame(
    pc_idx = seq_len(k_show),
    pc_lbl = factor(paste0("PC", seq_len(k_show)), levels = paste0("PC", seq_len(k_show))),
    pct_var = var_pcts[seq_len(k_show)],
    cum_var = cumsum(var_pcts[seq_len(k_show)])
  )
  scree_df$bar_vjust <- ifelse(scree_df$pct_var >= 22, 1.45, -0.35)
  p03a_scree <- ggplot2::ggplot(scree_df, ggplot2::aes(x = .data$pc_lbl, y = .data$pct_var)) +
    ggplot2::geom_col(fill = "#E5E7EB", color = "#374151", width = 0.65, linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(x = .data$pc_idx, y = .data$cum_var), color = "#009E73", linewidth = 0.85) +
    ggplot2::geom_point(ggplot2::aes(x = .data$pc_idx, y = .data$cum_var), shape = 2, size = 2.2, stroke = 0.85, color = "#009E73") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", .data$pct_var), vjust = .data$bar_vjust), size = 2.6, color = "#1F2937") +
    ggplot2::scale_y_continuous(limits = c(0, 108), expand = ggplot2::expansion(mult = c(0, 0.03))) +
    theme_packera_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = if (k_show > 8L) 35 else 0, hjust = if (k_show > 8L) 1 else 0.5)) +
    ggplot2::labs(
      title = sprintf("C. Dynamic Subspace Variance Spectrum (PC1..PC%d)", k_show),
      subtitle = sprintf("Green curve = cumulative variance (PC1 = %.1f%% to PC%d = %.1f%%)", scree_df$pct_var[1L], k_show, max(scree_df$cum_var, na.rm = TRUE)),
      x = "Principal Component Axis",
      y = "Variance Explained (%)"
    )

  plate_03a <- (p03a_pca / (p03a_power | p03a_scree)) +
    patchwork::plot_layout(heights = c(1.35, 1.0)) +
    patchwork::plot_annotation(
      title = "Step 03: Elliptic Fourier Harmonic Calibration & Symmetric Leaf Morphospace",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5))
    )
  fig03a_path <- save_packera_png_300(plate_03a, file.path(fdir, "fig03a_efa_morphospace_and_harmonics.png"), width = 13.4, height = 10.2)

  # ----------------------------------------------------------------------------
  # Plate 03B: fig03b_topological_ect_and_heteroblasty.png
  #   Panel A: Basal Topological ECT Ordination (ECT_PC1 vs ECT_PC2)
  #   Panel B: Authentic LeafMachine3 Radial ECT Exemplar Gallery (ect_overlays/)
  #   Panel C: Taxon-Stratified ECT Lobe Complexity / Dissection Index Boxplots
  #   Panel D: Basal-to-Cauline Heteroblasty Shift (or Within-Rosette Dispersion)
  # ----------------------------------------------------------------------------
  spec_df$ect1 <- if ("ECT_PC1" %in% names(spec_df)) suppressWarnings(as.numeric(spec_df$ECT_PC1)) else spec_df$pc1_num
  spec_df$ect2 <- if ("ECT_PC2" %in% names(spec_df)) suppressWarnings(as.numeric(spec_df$ECT_PC2)) else spec_df$pc2_num
  ok_ect <- is.finite(spec_df$ect1) & is.finite(spec_df$ect2)
  ect_sub <- spec_df[ok_ect, , drop = FALSE]

  p03b_ord <- ggplot2::ggplot(ect_sub, ggplot2::aes(x = .data$ect1, y = .data$ect2)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), size = 2.4, stroke = 0.85, alpha = 0.82) +
    ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
    ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group") +
    theme_packera_pub() +
    ggplot2::labs(
      title = "A. Topological Euler Characteristic Transform (ECT) Space",
      subtitle = "32 x 360 polar directional thresholds capturing deep pinnatifid lobing & sinuses",
      x = "Basal Topological ECT_PC1",
      y = "Basal Topological ECT_PC2"
    )
  if (length(active_groups) > 10L && nrow(ect_sub) > 0L) {
    cent_e <- stats::aggregate(cbind(ect1, ect2) ~ plot_group, data = ect_sub, FUN = stats::median)
    cent_e$lbl <- abbreviate_icn_taxon(cent_e$plot_group, abbreviate_genus = TRUE)
    p03b_ord <- add_centroid_repel_labels_gg(p03b_ord, cent_e, "ect1", "ect2", "lbl", color_col = "plot_group") +
      ggplot2::theme(legend.position = "none")
  }

  # Panel B: Native LeafMachine3 Radial ECT Gallery (embedded PNG grobs ordered by ECT_PC1 / lobing)
  lm3_man <- resolve_lm3_exemplar_manifest_r(fdir, zone = "basal")
  gallery_plots <- list()
  if (length(lm3_man) > 0L) {
    man_names <- intersect(names(lm3_man), active_groups)
    if (length(man_names) == 0L) man_names <- names(lm3_man)
    man_names <- utils::head(man_names, 12L)
    for (m_tx in man_names) {
      entry <- lm3_man[[m_tx]]
      grob_e <- read_png_raster_grob_safe(entry$image_path)
      if (!is.null(grob_e)) {
        gallery_plots[[length(gallery_plots) + 1L]] <- ggplot2::ggplot() +
          ggplot2::annotation_custom(grob_e, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
          ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
          ggplot2::theme_void() +
          ggplot2::labs(title = abbreviate_icn_taxon(m_tx, abbreviate_genus = TRUE)) +
          ggplot2::theme(plot.title = ggplot2::element_text(size = 7.6, face = "italic", hjust = 0.5))
      }
    }
  }
  p03b_gal <- if (length(gallery_plots) > 0L) {
    patchwork::wrap_plots(gallery_plots, ncol = min(4L, max(2L, ceiling(sqrt(length(gallery_plots)))))) +
      patchwork::plot_annotation(title = "B. Native LeafMachine3 Radial ECT Medoids")
  } else {
    # Clean fallback: ECT_PC1 vs Dissection Index hollow scatter
    spec_df$dis_num <- if ("dissection_index" %in% names(spec_df)) suppressWarnings(as.numeric(spec_df$dissection_index)) else (1.0 - suppressWarnings(as.numeric(spec_df$solidity)))
    d_sub <- spec_df[is.finite(spec_df$ect1) & is.finite(spec_df$dis_num), , drop = FALSE]
    ggplot2::ggplot(d_sub, ggplot2::aes(x = .data$dis_num, y = .data$ect1, color = .data$plot_group, shape = .data$plot_group)) +
      ggplot2::geom_point(size = 2.3, stroke = 0.85, alpha = 0.82, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_shape_manual(values = shp_map) +
      theme_packera_pub() +
      ggplot2::labs(title = "B. Topological ECT_PC1 vs. Foliar Dissection Index", x = "Perimeter Dissection Ratio (P / [2 * sqrt(pi * Area)])", y = "Basal ECT_PC1")
  }

  # Panel C: ECT Lobe Complexity (or Dissection Index) Boxplots + Hollow Jitter
  lobe_col <- if ("ect_lobe_complexity" %in% names(spec_df) && any(is.finite(suppressWarnings(as.numeric(spec_df$ect_lobe_complexity))))) "ect_lobe_complexity" else if ("dissection_index" %in% names(spec_df)) "dissection_index" else "solidity"
  spec_df$lobe_val <- suppressWarnings(as.numeric(spec_df[[lobe_col]]))
  lobe_sub <- spec_df[is.finite(spec_df$lobe_val), , drop = FALSE]
  ord_lvls <- if (nrow(lobe_sub) > 0L) names(sort(tapply(lobe_sub$lobe_val, lobe_sub$plot_group, stats::median, na.rm = TRUE))) else active_groups
  lobe_sub$grp_ord <- factor(lobe_sub$plot_group, levels = ord_lvls)
  ord_exprs <- format_icn_label_vec(ord_lvls, abbreviate_genus = TRUE, counts_map = cnt_map)

  p03b_lobe <- ggplot2::ggplot(lobe_sub, ggplot2::aes(x = .data$lobe_val, y = .data$grp_ord)) +
    ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = 0.55, linewidth = 0.5) +
    ggplot2::geom_jitter(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), height = 0.16, width = 0, size = 2.0, stroke = 0.8, alpha = 0.80, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = pal_map) +
    ggplot2::scale_shape_manual(values = shp_map) +
    ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE, counts_map = cnt_map)) +
    theme_packera_pub() +
    ggplot2::labs(
      title = "C. Foliar Lobe & Sinus Complexity by Taxon",
      x = if (lobe_col == "ect_lobe_complexity") "Topological ECT Lobe Complexity" else "Lamina Dissection Index",
      y = NULL
    )

  # Panel D: Basal-to-Cauline Heteroblasty Divergence (when cauline track present) or Within-Rosette Dispersion
  has_cauline <- ("cauline_leaf_count" %in% names(spec_df) && any(suppressWarnings(as.numeric(spec_df$cauline_leaf_count)) > 0, na.rm = TRUE)) ||
                 ("heteroblasty_ect_distance" %in% names(spec_df) && any(is.finite(suppressWarnings(as.numeric(spec_df$heteroblasty_ect_distance)))))
  het_col <- if (has_cauline && "heteroblasty_ect_distance" %in% names(spec_df) && any(is.finite(suppressWarnings(as.numeric(spec_df$heteroblasty_ect_distance))))) {
    "heteroblasty_ect_distance"
  } else if (has_cauline && "basal_vs_cauline_solidity_delta" %in% names(spec_df) && any(is.finite(suppressWarnings(as.numeric(spec_df$basal_vs_cauline_solidity_delta))))) {
    "basal_vs_cauline_solidity_delta"
  } else if ("within_rosette_efa_dispersion" %in% names(spec_df) && any(is.finite(suppressWarnings(as.numeric(spec_df$within_rosette_efa_dispersion))))) {
    "within_rosette_efa_dispersion"
  } else "foliar_variance"

  spec_df$het_val <- if (het_col %in% names(spec_df)) suppressWarnings(as.numeric(spec_df[[het_col]])) else rep(NA_real_, nrow(spec_df))
  het_sub <- spec_df[is.finite(spec_df$het_val), , drop = FALSE]
  p03b_het <- if (nrow(het_sub) > 0L) {
    het_sub$grp_ord <- factor(het_sub$plot_group, levels = intersect(ord_lvls, unique(het_sub$plot_group)))
    ggplot2::ggplot(het_sub, ggplot2::aes(x = .data$het_val, y = .data$grp_ord)) +
      ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = 0.55, linewidth = 0.5) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), height = 0.16, width = 0, size = 2.0, stroke = 0.8, alpha = 0.80, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_shape_manual(values = shp_map) +
      ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
      theme_packera_pub() +
      ggplot2::labs(
        title = if (has_cauline) "D. Basal-to-Cauline Heteroblastic Divergence" else "D. Within-Rosette Foliar Shape Dispersion",
        x = if (het_col == "heteroblasty_ect_distance") "Within-Voucher Basal-to-Cauline ECT Divergence" else if (exists("format_morph_trait_labels", mode = "function")) format_morph_trait_labels(het_col) else gsub("_", " ", het_col),
        y = NULL
      )
  } else {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Single-leaf rosette cohort\n(dispersion N/A)", color = "#6B7280", size = 3.4) + theme_packera_pub() + ggplot2::labs(title = "D. Within-Rosette Dispersion", x = NULL, y = NULL)
  }

  dyn_h_03b <- max(9.6, 5.2 + 0.28 * length(active_groups))
  plate_03b <- ((p03b_ord | p03b_gal) / (p03b_lobe | p03b_het)) +
    patchwork::plot_annotation(
      title = "Step 03: Topological Euler Characteristic Transform (ECT) & Heteroblasty",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5))
    )
  fig03b_path <- save_packera_png_300(plate_03b, file.path(fdir, "fig03b_topological_ect_and_heteroblasty.png"), width = 13.8, height = dyn_h_03b)

  invisible(c(fig03a_path, fig03b_path))
}
