# ==============================================================================
# Module: scripts/figures/step04_figures.R
# Step 04 Canonical 300-DPI PNG Figures (RRPP Allometry, GMM, Shrinkage CDA,
# Supervised EPA Level III Ecoregion CDA & Institutional Batch Variance Audit):
#   1. fig04a_foliar_allometry_audit.png
#   2. fig04b_gmm_cda_morphometrics.png
#   3. fig04c_ecoregion_and_batch_audit.png
# ==============================================================================

render_step04_allometry_png <- function(
  sub_df,
  r2_size = 0.0,
  p_size = 1.0,
  r2_int = 0.0,
  output_path = "outputs/figures/fig04a_foliar_allometry_audit.png",
  p_int = 1.0,
  pairwise_df = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (exists("resolve_active_figures_dir", mode = "function")) resolve_active_figures_dir(dirname(output_path), output_path) else dirname(output_path)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (old_f in c("shape_allometry_regression.pdf", "shape_allometry_regression.png", "Allometry_Lamina_Dissection_Trajectory.pdf", "Allometry_Lamina_Dissection_Trajectory.png")) {
    if (file.exists(file.path(fdir, old_f))) unlink(file.path(fdir, old_f), force = TRUE)
  }

  df <- as.data.frame(sub_df, stringsAsFactors = FALSE)
  if (nrow(df) == 0L) return(invisible(""))
  df$plot_group <- resolve_plot_group_labels(df)
  active_groups <- sort(unique(df$plot_group[!is.na(df$plot_group) & nzchar(df$plot_group)]))
  pal_map <- get_combinatorial_palette_r(active_groups)
  shp_map <- get_combinatorial_hollow_shapes_r(active_groups)
  cnt_map <- table(df$plot_group)
  lbl_exprs <- format_icn_label_vec(active_groups, abbreviate_genus = (length(active_groups) > 5L), counts_map = cnt_map)
  names(lbl_exprs) <- active_groups

  df$cs_x <- suppressWarnings(as.numeric(df$log_centsize))
  y_cand <- if ("PC1_raw" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$PC1_raw))))) "PC1_raw" else if ("PC1" %in% names(df)) "PC1" else "solidity"
  df$shape_y <- suppressWarnings(as.numeric(df[[y_cand]]))
  ok_a <- is.finite(df$cs_x) & is.finite(df$shape_y)
  sub_a <- df[ok_a, , drop = FALSE]

  sub_stats <- sprintf(
    "RRPP Size R\u00b2 = %.3f (p = %.4g) | Size x Taxon R\u00b2 = %.3f (p = %.4g; slopes %s)",
    ifelse(is.finite(r2_size), r2_size, 0.0),
    ifelse(is.finite(p_size), p_size, 1.0),
    ifelse(is.finite(r2_int), r2_int, 0.0),
    ifelse(is.finite(p_int), p_int, 1.0),
    ifelse(!is.finite(p_int) || p_int > 0.05, "homogeneous", "heterogeneous")
  )

  p_reg <- ggplot2::ggplot(sub_a, ggplot2::aes(x = .data$cs_x, y = .data$shape_y, color = .data$plot_group, shape = .data$plot_group)) +
    ggplot2::geom_point(size = 2.4, stroke = 0.85, alpha = 0.82) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.75, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
    ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group") +
    theme_packera_pub() +
    ggplot2::labs(
      title = "A. Foliar Allometry Trajectory (Ruler-Calibrated Cohort)",
      subtitle = sub_stats,
      x = expression(log[10]("Centroid Size [mm]")),
      y = sprintf("Primary Foliar Shape Score (%s)", y_cand)
    )
  if (length(active_groups) > 10L && nrow(sub_a) > 0L) {
    cent_a <- stats::aggregate(cbind(cs_x, shape_y) ~ plot_group, data = sub_a, FUN = stats::median)
    cent_a$lbl <- abbreviate_icn_taxon(cent_a$plot_group, abbreviate_genus = TRUE)
    p_reg <- add_centroid_repel_labels_gg(p_reg, cent_a, "cs_x", "shape_y", "lbl", color_col = "plot_group") +
      ggplot2::theme(legend.position = "none")
  }

  # Panel B: Pairwise Slope-Angle Matrix (if available) or Residualized Morphospace (PC1_corrected vs PC2_corrected)
  has_corr <- all(c("PC1_corrected", "PC2_corrected") %in% names(df)) &&
              any(is.finite(suppressWarnings(as.numeric(df$PC1_corrected)))) &&
              any(is.finite(suppressWarnings(as.numeric(df$PC2_corrected))))
  p_rhs <- if (has_corr) {
    df$pc1_c <- suppressWarnings(as.numeric(df$PC1_corrected))
    df$pc2_c <- suppressWarnings(as.numeric(df$PC2_corrected))
    sub_c <- df[is.finite(df$pc1_c) & is.finite(df$pc2_c), , drop = FALSE]
    el_g <- names(which(table(sub_c$plot_group) >= 4L))
    p_tmp <- ggplot2::ggplot(sub_c, ggplot2::aes(x = .data$pc1_c, y = .data$pc2_c, color = .data$plot_group, shape = .data$plot_group)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4)
    if (length(el_g) > 0L) {
      p_tmp <- p_tmp + ggplot2::stat_ellipse(
        data = sub_c[sub_c$plot_group %in% el_g, , drop = FALSE],
        ggplot2::aes(color = .data$plot_group),
        level = 0.95, linetype = "dashed", linewidth = 0.55, show.legend = FALSE
      )
    }
    p_tmp +
      ggplot2::geom_point(size = 2.4, stroke = 0.85, alpha = 0.82, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_shape_manual(values = shp_map) +
      theme_packera_pub() +
      ggplot2::labs(
        title = "B. Size-Free Allometric Residual Morphospace",
        subtitle = "Pooled within-taxon slope centered on authority-weighted reference mean size",
        x = "Allometry-Corrected PC1",
        y = "Allometry-Corrected PC2"
      )
  } else {
    ggplot2::ggplot(sub_a, ggplot2::aes(x = .data$cs_x, y = .data$plot_group, color = .data$plot_group, shape = .data$plot_group)) +
      ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = 0.55, inherit.aes = FALSE, ggplot2::aes(x = .data$cs_x, y = .data$plot_group)) +
      ggplot2::geom_jitter(height = 0.16, width = 0, size = 2.1, stroke = 0.8, alpha = 0.82, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_shape_manual(values = shp_map) +
      ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
      theme_packera_pub() +
      ggplot2::labs(title = "B. Calibrated Lamina Size Distribution by Taxon", x = expression(log[10]("Centroid Size [mm]")), y = NULL)
  }

  plate_04a <- (p_reg | p_rhs) +
    patchwork::plot_annotation(
      title = "Step 04: Multivariate Foliar Allometry & Ontogenetic Scaling Audit",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.2))
    )
  canon_04a <- file.path(fdir, "fig04a_foliar_allometry_audit.png")
  save_packera_png_300(plate_04a, canon_04a, width = 13.2, height = max(5.8, 3.8 + 0.22 * length(active_groups)), also_copy_to = output_path)
}

render_step04_gmm_cda_png <- function(
  audit_df,
  cda_res = NULL,
  gmm_res = NULL,
  output_path = "outputs/figures/fig04b_gmm_cda_morphometrics.png",
  feature_cols = NULL,
  manova_res = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (exists("resolve_active_figures_dir", mode = "function")) resolve_active_figures_dir(dirname(output_path), output_path) else dirname(output_path)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (old_stem in c("cda_passive_projection", "fig04_gmm_cda_clustering", "fig04a_gmm_bic_and_uncertainty", "fig04b_cda_biplot_and_allometry")) {
    for (ext in c(".pdf", ".png")) {
      p_old <- file.path(fdir, paste0(old_stem, ext))
      if (file.exists(p_old) && normalizePath(p_old, mustWork = FALSE) != normalizePath(normalize_png_path(output_path), mustWork = FALSE)) {
        unlink(p_old, force = TRUE)
      }
    }
  }

  df <- as.data.frame(audit_df, stringsAsFactors = FALSE)
  if (nrow(df) == 0L) return(invisible(""))
  df$plot_group <- resolve_plot_group_labels(df, cda_res = cda_res, gmm_res = gmm_res)
  active_groups <- sort(unique(df$plot_group[!is.na(df$plot_group) & nzchar(df$plot_group) & df$plot_group != "Unknown"]))
  if (length(active_groups) == 0L) active_groups <- unique(df$plot_group)
  df <- df[df$plot_group %in% active_groups, , drop = FALSE]

  pal_map <- get_combinatorial_palette_r(active_groups)
  shp_map <- get_combinatorial_hollow_shapes_r(active_groups)
  cnt_map <- table(df$plot_group)
  lbl_exprs <- format_icn_label_vec(active_groups, abbreviate_genus = (length(active_groups) > 5L), counts_map = cnt_map)
  names(lbl_exprs) <- active_groups

  manova_sub <- if (!is.null(manova_res) && exists("format_manova_pillai_subtitle", mode = "function")) {
    format_manova_pillai_subtitle(manova_res, group_pattern = "Species|Morphotype")
  } else ""

  # Panel A: Unsupervised GMM Stepwise 2*Delta*BIC Model Selection Curve
  bic_tb <- if (!is.null(gmm_res) && is.data.frame(gmm_res$bic_table)) gmm_res$bic_table else data.frame(K = 1:3, BIC = c(-10, -8, -12), Delta_BIC_prev = c(0, 2, -4), Two_Delta_BIC = c(0, 2, -2))
  bic_tb$K_num <- suppressWarnings(as.numeric(bic_tb$K))
  bic_tb <- bic_tb[order(bic_tb$K_num), , drop = FALSE]
  bic_tb$y_val <- if ("Delta_BIC_prev" %in% names(bic_tb) && any(is.finite(suppressWarnings(as.numeric(bic_tb$Delta_BIC_prev))))) {
    suppressWarnings(as.numeric(bic_tb$Delta_BIC_prev))
  } else if ("BIC" %in% names(bic_tb) && nrow(bic_tb) >= 2L) {
    c(0.0, diff(suppressWarnings(as.numeric(bic_tb$BIC))))
  } else if ("Two_Delta_BIC" %in% names(bic_tb)) {
    suppressWarnings(as.numeric(bic_tb$Two_Delta_BIC))
  } else rep(0.0, nrow(bic_tb))
  bic_tb$y_val[!is.finite(bic_tb$y_val)] <- 0.0
  best_k <- if (!is.null(gmm_res) && !is.null(gmm_res$best_k) && is.finite(gmm_res$best_k)) as.integer(gmm_res$best_k) else 1L
  is_merged_k1 <- (best_k == 1L) && (any(bic_tb$y_val[bic_tb$K_num >= 2] >= 10, na.rm = TRUE) ||
                   ("Evidence" %in% names(bic_tb) && any(grepl("Merged", as.character(bic_tb$Evidence), ignore.case = TRUE))))

  p_bic <- ggplot2::ggplot(bic_tb, ggplot2::aes(x = .data$K_num, y = .data$y_val)) +
    ggplot2::geom_hline(yintercept = c(0, 6, 10), linetype = c("solid", "dotted", "dashed"), color = c("#9CA3AF", "#6B7280", "#D55E00"), linewidth = 0.45) +
    ggplot2::geom_vline(xintercept = best_k, linetype = "dashed", color = "#0072B2", linewidth = 0.7) +
    ggplot2::geom_line(color = "#1F2937", linewidth = 0.85) +
    ggplot2::geom_point(shape = 1, size = 2.6, stroke = 0.95, color = "#0072B2") +
    ggplot2::scale_x_continuous(breaks = unique(as.integer(bic_tb$K_num))) +
    theme_packera_pub() +
    ggplot2::labs(
      title = sprintf("A. Unsupervised GMM Selection (Selected K* = %d)", best_k),
      subtitle = "Deviance-scale Kass & Raftery threshold (stepwise 2*Delta*BIC >= 6 strong, >= 10 decisive)",
      x = "Number of Mixture Components (K)",
      y = "Stepwise 2 * Delta BIC (vs. K - 1)"
    )
  if (is_merged_k1) {
    mid_k <- mean(range(bic_tb$K_num, na.rm = TRUE))
    mid_y <- 0.55 * max(bic_tb$y_val, 12.0, na.rm = TRUE)
    p_bic <- p_bic + ggplot2::annotate(
      "label", x = mid_k, y = mid_y,
      label = "Hennig (2010) unimodal ridgeline\nmerging -> K* = 1 retained",
      size = 2.75, color = "#0072B2", fill = grDevices::adjustcolor("white", alpha.f = 0.92)
    )
  }

  # Panel B: Shrinkage CDA Ordination (2D Can1 vs Can2 when G >= 3; 1D Can1 density + hollow raincloud when G == 2)
  c1 <- if ("can1" %in% names(df)) suppressWarnings(as.numeric(df$can1)) else if ("PC1" %in% names(df)) suppressWarnings(as.numeric(df$PC1)) else rep(NA_real_, nrow(df))
  c2 <- if ("can2" %in% names(df)) suppressWarnings(as.numeric(df$can2)) else rep(NA_real_, nrow(df))
  has_2d_can <- sum(is.finite(c1) & is.finite(c2)) >= 3L && isTRUE(stats::sd(c2, na.rm = TRUE) > 1e-6) && length(active_groups) >= 3L
  if (!has_2d_can && length(active_groups) != 2L && "PC2" %in% names(df)) {
    c2_fb <- suppressWarnings(as.numeric(df$PC2))
    if (sum(is.finite(c1) & is.finite(c2_fb)) >= 3L) {
      c2 <- c2_fb
      has_2d_can <- TRUE
    }
  }
  df$ord1 <- c1
  df$ord2 <- c2

  loocv_val <- if (!is.null(cda_res) && "loocv_accuracy" %in% names(cda_res)) {
    suppressWarnings(as.numeric(cda_res$loocv_accuracy[1L]))
  } else if ("loocv_accuracy" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$loocv_accuracy))))) {
    stats::median(suppressWarnings(as.numeric(df$loocv_accuracy)), na.rm = TRUE)
  } else NA_real_
  loocv_tag_s04 <- if (is.finite(loocv_val)) sprintf("Deme-blocked CV accuracy = %.1f%%", 100.0 * loocv_val) else "Descriptive Post-GMM Ordination"
  cda_sub_raw <- paste(
    c(
      loocv_tag_s04,
      if (isTRUE(nzchar(manova_sub))) manova_sub else NULL
    ),
    collapse = " | "
  )
  cda_sub_line <- if (nzchar(cda_sub_raw)) paste(strwrap(cda_sub_raw, width = 72L), collapse = "\n") else ""

  p_cda <- if (length(active_groups) == 2L && sum(is.finite(df$ord1)) >= 3L && !has_2d_can) {
    # Standard Koutecky (2015) MorphoTools2 G = 2 protocol: 1D Canonical Axis Density + Hollow Raincloud
    sub_1d <- df[is.finite(df$ord1), , drop = FALSE]
    ggplot2::ggplot(sub_1d, ggplot2::aes(x = .data$ord1, color = .data$plot_group)) +
      ggplot2::geom_density(linewidth = 0.9, fill = NA) +
      ggplot2::geom_jitter(
        ggplot2::aes(y = -0.04, shape = .data$plot_group),
        height = 0.025, width = 0, size = 2.3, stroke = 0.85, alpha = 0.82
      ) +
      ggplot2::geom_hline(yintercept = 0, color = "#9CA3AF", linewidth = 0.4) +
      ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
      ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group") +
      theme_packera_pub() +
      ggplot2::labs(
        title = "B. Shrinkage Canonical Discriminant Distribution (G = 2 Taxa)",
        subtitle = if (isTRUE(nzchar(cda_sub_line))) cda_sub_line else "Single canonical discriminant axis (Can1) density & hollow specimen scores",
        x = "Canonical Variate 1 (Can1)",
        y = "Kernel Density"
      )
  } else {
    sub_2d <- df[is.finite(df$ord1) & is.finite(df$ord2), , drop = FALSE]
    el_g <- names(which(table(sub_2d$plot_group) >= 4L))
    p_tmp <- ggplot2::ggplot(sub_2d, ggplot2::aes(x = .data$ord1, y = .data$ord2)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4)
    if (length(el_g) > 0L) {
      p_tmp <- p_tmp + ggplot2::stat_ellipse(
        data = sub_2d[sub_2d$plot_group %in% el_g, , drop = FALSE],
        ggplot2::aes(color = .data$plot_group),
        level = 0.95, linetype = "dashed", linewidth = 0.55, show.legend = FALSE
      )
    }
    p_tmp <- p_tmp +
      ggplot2::geom_point(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), size = 2.4, stroke = 0.85, alpha = 0.82) +
      ggplot2::scale_color_manual(values = pal_map, labels = lbl_exprs, name = "Taxon / Group") +
      ggplot2::scale_shape_manual(values = shp_map, labels = lbl_exprs, name = "Taxon / Group") +
      theme_packera_pub() +
      ggplot2::labs(
        title = "B. Shrinkage Canonical Discriminant Biplot (Active + Passive Projection)",
        subtitle = if (isTRUE(nzchar(cda_sub_line))) cda_sub_line else "Kish-N_eff Ledoit-Wolf shrinkage CDA (hollow markers; * = Nomenclatural Type)",
        x = if ("can1" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$can1))))) "Canonical Variate 1 (Can1)" else "Morphometric PC1",
        y = if ("can2" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$can2))))) "Canonical Variate 2 (Can2)" else "Morphometric PC2"
      )

    # Overlay non-overlapping canonical structure loading vectors across dynamic selected features
    cand_f <- unique(c(feature_cols, if ("selected_morph_features" %in% names(df)) unlist(strsplit(as.character(stats::na.omit(df$selected_morph_features)), ";")) else character(0)))
    cand_f <- intersect(trimws(cand_f), names(sub_2d))
    cand_f <- cand_f[vapply(cand_f, function(cl) sum(is.finite(suppressWarnings(as.numeric(sub_2d[[cl]])))) >= 3L, logical(1))]
    if (length(cand_f) >= 2L && nrow(sub_2d) >= 5L) {
      num_mat_2d <- do.call(cbind, lapply(cand_f, function(cl) suppressWarnings(as.numeric(sub_2d[[cl]]))))
      colnames(num_mat_2d) <- cand_f
      cor_mat <- tryCatch(suppressWarnings(stats::cor(num_mat_2d, cbind(sub_2d$ord1, sub_2d$ord2), use = "pairwise.complete.obs")), error = function(e) NULL)
      if (is.matrix(cor_mat) && ncol(cor_mat) == 2L) {
        mag <- sqrt(rowSums(cor_mat^2))
        ok_v <- which(is.finite(mag) & mag > 0.18)
        if (length(ok_v) > 0L) {
          top_v <- ok_v[order(mag[ok_v], decreasing = TRUE)[seq_len(min(6L, length(ok_v)))]]
          rad_x <- 0.72 * max(abs(range(sub_2d$ord1, na.rm = TRUE)))
          rad_y <- 0.72 * max(abs(range(sub_2d$ord2, na.rm = TRUE)))
          vec_df <- data.frame(
            var_raw = rownames(cor_mat)[top_v],
            vx = cor_mat[top_v, 1L] * rad_x,
            vy = cor_mat[top_v, 2L] * rad_y,
            stringsAsFactors = FALSE
          )
          vec_df$lbl <- if (exists("format_morph_trait_labels", mode = "function")) {
            format_morph_trait_labels(vec_df$var_raw)
          } else if (exists("format_cda_loading_labels", mode = "function")) {
            format_cda_loading_labels(vec_df$var_raw)
          } else gsub("_std$", "", vec_df$var_raw)
          p_tmp <- p_tmp +
            ggplot2::geom_segment(
              data = vec_df,
              ggplot2::aes(x = 0, y = 0, xend = .data$vx, yend = .data$vy),
              arrow = grid::arrow(length = grid::unit(0.14, "cm")),
              color = "#374151", linewidth = 0.6, inherit.aes = FALSE
            )
          p_tmp <- add_centroid_repel_labels_gg(p_tmp, vec_df, "vx", "vy", "lbl", size = 2.65)
        }
      }
    }
    if ("is_type_specimen" %in% names(sub_2d)) {
      is_t <- !is.na(sub_2d$is_type_specimen) & as.logical(sub_2d$is_type_specimen)
      if (any(is_t)) {
        p_tmp <- p_tmp + ggplot2::geom_point(
          data = sub_2d[is_t, , drop = FALSE],
          ggplot2::aes(x = .data$ord1, y = .data$ord2, color = .data$plot_group),
          shape = 8, size = 3.6, stroke = 1.1, show.legend = FALSE
        )
      }
    }
    if (length(active_groups) > 10L && nrow(sub_2d) > 0L) {
      cent_c <- stats::aggregate(cbind(ord1, ord2) ~ plot_group, data = sub_2d, FUN = stats::median)
      cent_c$lbl <- abbreviate_icn_taxon(cent_c$plot_group, abbreviate_genus = TRUE)
      p_tmp <- add_centroid_repel_labels_gg(p_tmp, cent_c, "ord1", "ord2", "lbl", color_col = "plot_group") +
        ggplot2::theme(legend.position = "none")
    }
    p_tmp
  }

  # Panel C: Dynamic Diagnostic Scale-Invariant Trait Boxplots
  cand_traits <- c("AR_inv", "petiole_to_blade_ratio", "scape_to_basal_leaf_ratio", "q25_solidity", "dissection_index", "ECT_PC1", "aspect_ratio", "PC1")
  avail_traits <- cand_traits[vapply(cand_traits, function(tc) tc %in% names(df) && sum(is.finite(suppressWarnings(as.numeric(df[[tc]])))) >= 3L, logical(1))]
  show_traits <- utils::head(avail_traits, 4L)
  p_box <- if (length(show_traits) > 0L) {
    long_list <- lapply(show_traits, function(tc) {
      vv <- suppressWarnings(as.numeric(df[[tc]]))
      ok_t <- is.finite(vv)
      data.frame(
        plot_group = df$plot_group[ok_t],
        trait = paste(strwrap(if (exists("format_morph_trait_labels", mode = "function")) format_morph_trait_labels(tc) else if (exists("format_cda_loading_labels", mode = "function")) format_cda_loading_labels(tc) else tc, width = 28L), collapse = "\n"),
        val = vv[ok_t],
        stringsAsFactors = FALSE
      )
    })
    long_df <- do.call(rbind, long_list)
    ggplot2::ggplot(long_df, ggplot2::aes(x = .data$val, y = .data$plot_group)) +
      ggplot2::geom_boxplot(fill = NA, color = "#4B5563", outlier.shape = NA, width = 0.55, linewidth = 0.45) +
      ggplot2::geom_jitter(ggplot2::aes(color = .data$plot_group, shape = .data$plot_group), height = 0.16, width = 0, size = 1.75, stroke = 0.75, alpha = 0.78, show.legend = FALSE) +
      ggplot2::facet_wrap(~ trait, scales = "free_x", ncol = min(2L, length(show_traits))) +
      ggplot2::scale_color_manual(values = pal_map) +
      ggplot2::scale_shape_manual(values = shp_map) +
      ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
      theme_packera_pub(base_size = 9.4) +
      ggplot2::labs(title = "C. Diagnostic Scale-Invariant Traits", x = "Trait Value", y = NULL)
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "C. Diagnostic Traits")
  }

  # Panel D: Taxon x Unsupervised GMM Matrix (multi-group/K >= 2) or OSR / Self-Mahalanobis D^2 Tail Audit (K* = 1)
  n_noise <- if ("is_confident_label_noise" %in% names(df)) sum(!is.na(df$is_confident_label_noise) & as.logical(df$is_confident_label_noise)) else 0L
  n_osr <- if ("open_set_novelty_p" %in% names(df)) sum(suppressWarnings(as.numeric(df$open_set_novelty_p)) < 0.01, na.rm = TRUE) else 0L
  d2_col <- intersect(c("self_morph_mahal_d2", "type_morph_d2", "open_set_min_congener_d2", "min_congener_mahal_d2"), names(df))[1L]
  p_cont <- if (length(active_groups) == 1L && best_k <= 1L && !is.na(d2_col) && sum(is.finite(suppressWarnings(as.numeric(df[[d2_col]])))) >= 5L) {
    d2_vals <- suppressWarnings(as.numeric(df[[d2_col]]))
    d2_vals <- d2_vals[is.finite(d2_vals) & d2_vals >= 0]
    df_chi <- if ("self_morph_df" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$self_morph_df))))) {
      stats::median(suppressWarnings(as.numeric(df$self_morph_df)), na.rm = TRUE)
    } else if ("open_set_df" %in% names(df) && any(is.finite(suppressWarnings(as.numeric(df$open_set_df))))) {
      stats::median(suppressWarnings(as.numeric(df$open_set_df)), na.rm = TRUE)
    } else 6.0
    df_chi <- max(1.0, df_chi)
    chi95 <- stats::qchisq(0.95, df = df_chi)
    n_tail <- sum(d2_vals > chi95)
    d2_df <- data.frame(d2 = d2_vals)
    ggplot2::ggplot(d2_df, ggplot2::aes(x = .data$d2)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 24L, fill = "#E2E8F0", color = "#334155", linewidth = 0.4) +
      ggplot2::geom_density(color = "#0072B2", linewidth = 0.9) +
      ggplot2::geom_vline(xintercept = chi95, linetype = "dashed", color = "#D55E00", linewidth = 0.75) +
      ggplot2::annotate("label", x = chi95, y = Inf, vjust = 1.4, hjust = 0.05, label = sprintf("chi^2_0.95 (df=%g) = %.1f", df_chi, chi95), size = 2.65, color = "#D55E00", fill = "white") +
      theme_packera_pub(base_size = 9.5) +
      ggplot2::labs(
        title = "D. Unimodal Self-Mahalanobis D^2 & OSR Tail Audit (K* = 1)",
        subtitle = sprintf("95%% ellipsoid outliers = %d / %d | Confident label-noise = %d | OSR (p < 0.01) = %d", n_tail, length(d2_vals), n_noise, n_osr),
        x = sprintf("Morphometric Squared Mahalanobis Distance (%s)", gsub("_", " ", d2_col)),
        y = "Kernel Density"
      )
  } else {
    cl_vec <- if ("gmm_cluster" %in% names(df)) as.character(df$gmm_cluster) else rep("1", nrow(df))
    cl_vec[is.na(cl_vec) | !nzchar(cl_vec)] <- "Unassigned"
    tab_df <- as.data.frame(table(Taxon = df$plot_group, Cluster = paste0("G", cl_vec)), stringsAsFactors = FALSE)
    row_tots <- stats::ave(tab_df$Freq, tab_df$Taxon, FUN = function(x) max(sum(x), 1L))
    tab_df$Pct <- 100.0 * tab_df$Freq / row_tots
    tab_pos <- tab_df[tab_df$Freq > 0, , drop = FALSE]
    if (nrow(tab_pos) > 0L) tab_pos$lbl <- sprintf("%d\n(%.0f%%)", tab_pos$Freq, tab_pos$Pct)
    pt <- ggplot2::ggplot(tab_df, ggplot2::aes(x = .data$Cluster, y = .data$Taxon, fill = .data$Pct)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5)
    if (nrow(tab_pos) > 0L) {
      pt <- pt + ggplot2::geom_text(
        data = tab_pos,
        ggplot2::aes(x = .data$Cluster, y = .data$Taxon, label = .data$lbl),
        size = if (length(active_groups) > 15L) 2.2 else 2.7,
        color = ifelse(tab_pos$Pct > 55, "white", "#111827"),
        inherit.aes = FALSE
      )
    }
    pt +
      ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 100), name = "Row %") +
      ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
      theme_packera_pub(base_size = 9.6) +
      ggplot2::labs(
        title = "D. Taxon x GMM Cluster Concordance",
        subtitle = sprintf("Confident label-noise flags = %d | Open-set novelty (p < 0.01) = %d", n_noise, n_osr),
        x = "Unsupervised Gaussian Mixture Cluster",
        y = NULL
      )
  }

  dyn_h_04b <- max(10.2, 6.2 + 0.28 * length(active_groups))
  plate_04b <- ((p_bic | p_cda) / (p_box | p_cont)) +
    patchwork::plot_layout(widths = c(1.0, 1.25), heights = c(1.05, 1.15)) +
    patchwork::plot_annotation(
      title = "Step 04: Unsupervised Gaussian Mixture Clustering & Shrinkage Discriminant Audit",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5))
    )
  canon_04b <- file.path(fdir, "fig04b_gmm_cda_morphometrics.png")
  save_packera_png_300(plate_04b, canon_04b, width = 14.2, height = dyn_h_04b, also_copy_to = output_path)
}

render_step04_ecoregion_batch_png <- function(
  audit_df,
  eco_cda_res = NULL,
  gmm_res = NULL,
  figures_dir = "outputs/figures",
  manova_res = NULL,
  batch_csv = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- as.character(figures_dir)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (ext in c(".pdf", ".png")) {
    p_old <- file.path(fdir, paste0("fig04c_supervised_ecoregion_cda", ext))
    if (file.exists(p_old)) unlink(p_old, force = TRUE)
  }

  df <- as.data.frame(audit_df, stringsAsFactors = FALSE)
  eco_vec <- if ("regional_group" %in% names(df)) as.character(df$regional_group) else rep("Unassigned", nrow(df))
  eco_vec[is.na(eco_vec) | !nzchar(trimws(eco_vec)) | eco_vec == "NA"] <- "Unassigned"
  df$eco_group <- gsub("_", " ", eco_vec)
  eco_sub <- df[df$eco_group != "Unassigned", , drop = FALSE]
  if (nrow(eco_sub) == 0L) eco_sub <- df

  eco_lvls <- sort(unique(eco_sub$eco_group))
  pal_eco <- get_combinatorial_palette_r(eco_lvls)
  shp_eco <- get_combinatorial_hollow_shapes_r(eco_lvls)
  eco_manova_sub <- if (!is.null(manova_res) && exists("format_manova_pillai_subtitle", mode = "function")) {
    format_manova_pillai_subtitle(manova_res, group_pattern = "Ecoregion|Region")
  } else ""

  has_ec1 <- "ecoregion_can1" %in% names(eco_sub) && any(is.finite(suppressWarnings(as.numeric(eco_sub$ecoregion_can1))))
  has_ec2 <- "ecoregion_can2" %in% names(eco_sub) && any(is.finite(suppressWarnings(as.numeric(eco_sub$ecoregion_can2))))
  ec1 <- if (has_ec1) suppressWarnings(as.numeric(eco_sub$ecoregion_can1)) else if ("PC1" %in% names(eco_sub)) suppressWarnings(as.numeric(eco_sub$PC1)) else rep(NA_real_, nrow(eco_sub))
  ec2 <- if (has_ec2) suppressWarnings(as.numeric(eco_sub$ecoregion_can2)) else if ("PC2" %in% names(eco_sub)) suppressWarnings(as.numeric(eco_sub$PC2)) else rep(NA_real_, nrow(eco_sub))
  eco_sub$ec1 <- ec1; eco_sub$ec2 <- ec2

  p_eco <- if (length(eco_lvls) == 2L && sum(is.finite(ec1)) >= 3L && (!any(is.finite(ec2)) || isTRUE(stats::sd(ec2, na.rm = TRUE) < 1e-6))) {
    s1 <- eco_sub[is.finite(eco_sub$ec1), , drop = FALSE]
    ggplot2::ggplot(s1, ggplot2::aes(x = .data$ec1, color = .data$eco_group)) +
      ggplot2::geom_density(linewidth = 0.9, fill = NA) +
      ggplot2::geom_jitter(ggplot2::aes(y = -0.03, shape = .data$eco_group), height = 0.02, width = 0, size = 2.2, stroke = 0.85) +
      ggplot2::scale_color_manual(values = pal_eco, name = "EPA Level III Ecoregion") +
      ggplot2::scale_shape_manual(values = shp_eco, name = "EPA Level III Ecoregion") +
      theme_packera_pub() +
      ggplot2::labs(title = "A. Supervised EPA Level III Ecoregion CDA (1D)", subtitle = if (isTRUE(nzchar(eco_manova_sub))) eco_manova_sub else NULL, x = "Ecoregion Can1", y = "Density")
  } else if (sum(is.finite(eco_sub$ec1) & is.finite(eco_sub$ec2)) >= 2L) {
    s2 <- eco_sub[is.finite(eco_sub$ec1) & is.finite(eco_sub$ec2), , drop = FALSE]
    p_t <- ggplot2::ggplot(s2, ggplot2::aes(x = .data$ec1, y = .data$ec2, color = .data$eco_group, shape = .data$eco_group)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
      ggplot2::geom_point(size = 2.4, stroke = 0.85, alpha = 0.82) +
      ggplot2::scale_color_manual(values = pal_eco, name = "EPA Level III Ecoregion") +
      ggplot2::scale_shape_manual(values = shp_eco, name = "EPA Level III Ecoregion") +
      theme_packera_pub() +
      ggplot2::labs(
        title = "A. Supervised EPA Level III Ecoregion Discriminant Space",
        subtitle = if (isTRUE(nzchar(eco_manova_sub))) eco_manova_sub else "Tests geographic/ecotypic structuring across EPA Level III ecoregions",
        x = if (has_ec1) "Ecoregion Canonical Variate 1" else "Morphometric PC1",
        y = if (has_ec2) "Ecoregion Canonical Variate 2" else "Morphometric PC2"
      )
    if (length(eco_lvls) > 10L && nrow(s2) > 0L) {
      tab_er <- sort(table(s2$eco_group), decreasing = TRUE)
      top_er <- names(tab_er[tab_er >= 5L])
      if (length(top_er) == 0L) top_er <- names(tab_er)
      top_er <- utils::head(top_er, 10L)
      cent_er <- stats::aggregate(cbind(ec1, ec2) ~ eco_group, data = s2[s2$eco_group %in% top_er, , drop = FALSE], FUN = stats::median)
      if (requireNamespace("ggrepel", quietly = TRUE) && nrow(cent_er) > 0L) {
        p_t <- p_t + ggrepel::geom_label_repel(
          data = cent_er,
          ggplot2::aes(x = .data$ec1, y = .data$ec2, label = .data$eco_group, color = .data$eco_group),
          fill = grDevices::adjustcolor("white", alpha.f = 0.88),
          size = 2.65, max.overlaps = 14L, box.padding = 0.45, force = 3.0,
          show.legend = FALSE, inherit.aes = FALSE
        ) + ggplot2::theme(legend.position = "none")
      } else {
        p_t <- add_centroid_repel_labels_gg(p_t, cent_er, "ec1", "ec2", "eco_group", color_col = "eco_group") +
          ggplot2::theme(legend.position = "none")
      }
    }
    p_t
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Insufficient multi-ecoregion coordinates\nfor supervised ecoregion CDA", color = "#6B7280", size = 3.3) +
      theme_packera_pub() +
      ggplot2::labs(title = "A. Supervised EPA Level III Ecoregion Space", x = NULL, y = NULL)
  }

  # Panel B: Institutional Herbarium Batch-Effect Variance Audit
  cand_b <- unique(c(
    batch_csv,
    file.path(dirname(fdir), "reports", "institutional_batch_effect_audit.csv"),
    file.path(dirname(fdir), "reports", "batch_variance_audit.csv"),
    file.path(dirname(dirname(fdir)), "data", "tables", "batch_variance_audit.csv"),
    "outputs/reports/institutional_batch_effect_audit.csv",
    "data/tables/batch_variance_audit.csv"
  ))
  cand_b <- cand_b[!is.na(cand_b) & nzchar(cand_b) & file.exists(cand_b)]
  batch_df <- if (length(cand_b) > 0L) tryCatch(utils::read.csv(cand_b[1L], stringsAsFactors = FALSE), error = function(e) NULL) else NULL

  p_batch <- if (!is.null(batch_df) && nrow(batch_df) > 0L) {
    if ("Term" %in% names(batch_df) && any(batch_df$Term == "institutionCode", na.rm = TRUE)) {
      batch_df <- batch_df[!is.na(batch_df$Term) & batch_df$Term == "institutionCode", , drop = FALSE]
    }
    f_col <- intersect(c("Response", "Feature", "Variable", "trait", "Term"), names(batch_df))[1L]
    r2_col <- intersect(c("Partial_R2_Institution", "partial_eta_sq", "R_squared", "R2_Institution", "R2", "eta_sq"), names(batch_df))[1L]
    q_col <- intersect(c("p_value_bh", "p_adj_bh", "q_value_BH", "p_adj", "q_value", "p_value", "P_value"), names(batch_df))[1L]
    if (!is.na(f_col) && !is.na(r2_col)) {
      bdf <- data.frame(
        feat = as.character(batch_df[[f_col]]),
        r2 = suppressWarnings(as.numeric(batch_df[[r2_col]])),
        qval = if (!is.na(q_col)) suppressWarnings(as.numeric(batch_df[[q_col]])) else rep(1.0, nrow(batch_df)),
        stringsAsFactors = FALSE
      )
      bdf <- bdf[is.finite(bdf$r2) & !is.na(bdf$feat) & nzchar(bdf$feat), , drop = FALSE]
      bdf <- bdf[!duplicated(bdf$feat), , drop = FALSE]
    } else bdf <- data.frame()
    if (nrow(bdf) > 0L) {
      bdf$sig <- ifelse(!is.na(bdf$qval) & bdf$qval < 0.05 & bdf$r2 >= 0.05, "Flagged (R2 >= 0.05, q < 0.05)", "Negligible Batch Effect")
      bdf$feat_lbl <- if (exists("format_morph_trait_labels", mode = "function")) format_morph_trait_labels(bdf$feat) else bdf$feat
      bdf$feat_lbl <- factor(bdf$feat_lbl, levels = unique(bdf$feat_lbl[order(bdf$r2)]))
      ggplot2::ggplot(bdf, ggplot2::aes(x = .data$r2, y = .data$feat_lbl, color = .data$sig)) +
        ggplot2::geom_vline(xintercept = 0.05, linetype = "dashed", color = "#D55E00", linewidth = 0.65) +
        ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$r2, y = .data$feat_lbl, yend = .data$feat_lbl), linewidth = 0.75) +
        ggplot2::geom_point(shape = 1, size = 2.8, stroke = 1.0) +
        ggplot2::scale_color_manual(values = c("Negligible Batch Effect" = "#0072B2", "Flagged (R2 >= 0.05, q < 0.05)" = "#D55E00"), name = "Batch Audit") +
        theme_packera_pub() +
        ggplot2::labs(
          title = "B. Institutional Herbarium Batch-Effect Audit",
          subtitle = "Taxon-controlled ANCOVA partial R\u00b2 across herbaria (dashed line = 5% threshold)",
          x = "Institutional Partial R\u00b2",
          y = NULL
        )
    } else {
      ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "B. Institutional Batch-Effect Audit")
    }
  } else {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Single-institution cohort\n(institutional batch variance N/A)", color = "#6B7280", size = 3.4) +
      theme_packera_pub() +
      ggplot2::labs(title = "B. Institutional Batch-Effect Audit", x = NULL, y = NULL)
  }

  plate_04c <- (p_eco | p_batch) +
    patchwork::plot_annotation(
      title = "Step 04: Ecoregional Morphometric Structuring & Herbarium Preservation Audit",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.2))
    )
  save_packera_png_300(plate_04c, file.path(fdir, "fig04c_ecoregion_and_batch_audit.png"), width = 13.6, height = 6.4)
}
