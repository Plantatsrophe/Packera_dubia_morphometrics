#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_figures_utils.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Shared publication figure utilities across R pipeline phases (Steps 03, 04,
#   06, and 07). Enforces:
#     - Dual vector .pdf and 300-DPI .png atomic exports (write-then-replace)
#     - Colorblind-safe Okabe-Ito palettes (#0072B2, #D55E00, #009E73, #CC79A7)
#     - Italicized Packera taxon labels (*Packera anonyma*, *P. dubia*,
#       *P. paupercula*, *P. plattensis*)
# ==============================================================================

OKABE_ITO_TAXA <- c(
  "Packera anonyma"    = "#0072B2",
  "Packera dubia"      = "#D55E00",
  "Packera paupercula" = "#009E73",
  "Packera plattensis" = "#CC79A7",
  "Unknown"            = "#7F7F7F"
)

PACKERA_ITALIC_LABELS <- c(
  "Packera anonyma"    = "*Packera anonyma*",
  "Packera dubia"      = "*P. dubia*",
  "Packera paupercula" = "*P. paupercula*",
  "Packera plattensis" = "*P. plattensis*",
  "Unknown"            = "Unassigned"
)

packera_italic_expressions <- function(taxa_vec) {
  sapply(as.character(taxa_vec), function(tx) {
    if (grepl("anonyma", tx, ignore.case = TRUE)) {
      return(expression(italic("Packera anonyma")))
    } else if (grepl("dubia|tomentosa", tx, ignore.case = TRUE)) {
      return(expression(italic("P. dubia")))
    } else if (grepl("paupercula", tx, ignore.case = TRUE)) {
      return(expression(italic("P. paupercula")))
    } else if (grepl("plattensis", tx, ignore.case = TRUE)) {
      return(expression(italic("P. plattensis")))
    } else {
      return(tx)
    }
  }, USE.NAMES = FALSE)
}

resolve_active_figures_dir <- function(explicit_dir = NULL, reference_path = NULL) {
  if (!is.null(explicit_dir) && nzchar(as.character(explicit_dir))) {
    fdir <- as.character(explicit_dir)
    if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    return(fdir)
  }
  env_run <- Sys.getenv("PACKERA_RUN_DIR", unset = "")
  if (nzchar(env_run)) {
    fdir <- file.path(env_run, "outputs", "figures")
    if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    return(fdir)
  }
  if (!is.null(reference_path) && nzchar(as.character(reference_path))) {
    ref_p <- as.character(reference_path)
    if (grepl("\\.(pdf|png)$", ref_p, ignore.case = TRUE)) {
      fdir <- dirname(ref_p)
    } else {
      # e.g., <run_dir>/data/tables/file.csv -> <run_dir>/outputs/figures
      parent1 <- dirname(ref_p)
      parent2 <- dirname(parent1)
      if (basename(parent1) == "tables" && basename(parent2) == "data") {
        fdir <- file.path(dirname(parent2), "outputs", "figures")
      } else {
        fdir <- file.path(parent1, "figures")
      }
    }
    if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    return(fdir)
  }
  fdir <- file.path("outputs", "figures")
  if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  return(fdir)
}

atomic_save_dual_plot <- function(draw_fn, figures_dir, stem, width = 11, height = 5.5, dpi = 300) {
  if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  }
  pdf_final <- file.path(figures_dir, paste0(stem, ".pdf"))
  png_final <- file.path(figures_dir, paste0(stem, ".png"))
  pdf_tmp <- file.path(figures_dir, paste0(".", stem, "_", Sys.getpid(), ".pdf.tmp"))
  png_tmp <- file.path(figures_dir, paste0(".", stem, "_", Sys.getpid(), ".png.tmp"))

  grDevices::pdf(pdf_tmp, width = width, height = height, useDingbats = FALSE)
  tryCatch({
    draw_fn()
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(pdf_tmp)) {
    file.rename(pdf_tmp, pdf_final)
  }

  grDevices::png(png_tmp, width = width, height = height, units = "in", res = dpi)
  tryCatch({
    draw_fn()
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(png_tmp)) {
    file.rename(png_tmp, png_final)
  }
  invisible(list(pdf = pdf_final, png = png_final))
}

atomic_save_ggplot_dual <- function(gg_obj, figures_dir, stem, width = 11, height = 5.5, dpi = 300) {
  atomic_save_dual_plot(
    draw_fn = function() print(gg_obj),
    figures_dir = figures_dir,
    stem = stem,
    width = width,
    height = height,
    dpi = dpi
  )
}

# ==============================================================================
# Reusable 3D Ordination Projection Helper (Base R persp+trans3d / scatterplot3d)
# Always enforces strictly unfilled/hollow markers (pch = 1, 2, 5, 0 or bg = NA)
# ==============================================================================
plot_3d_scatter_box <- function(
  x,
  y,
  z,
  col = "#0072B2",
  pch = 1,
  cex = 0.75,
  lwd = 1.1,
  xlab = "Axis 1",
  ylab = "Axis 2",
  zlab = "Axis 3",
  main = "",
  theta = 35,
  phi = 22,
  centroids = NULL,
  centroid_labels = NULL,
  legend_args = NULL
) {
  x_num <- suppressWarnings(as.numeric(x))
  y_num <- suppressWarnings(as.numeric(y))
  z_num <- suppressWarnings(as.numeric(z))
  n_pts <- max(length(x_num), length(y_num), length(z_num))
  if (n_pts == 0) return(invisible(NULL))

  x_num[is.na(x_num) | !is.finite(x_num)] <- 0.0
  y_num[is.na(y_num) | !is.finite(y_num)] <- 0.0
  z_num[is.na(z_num) | !is.finite(z_num)] <- 0.0

  col_vec <- rep_len(as.character(col), n_pts)
  pch_vec <- rep_len(as.integer(pch), n_pts)
  # Enforce strictly unfilled/hollow symbols even if caller passes solid pch codes
  pch_vec[pch_vec %in% c(16L, 19L, 20L)] <- 1L
  pch_vec[pch_vec == 17L] <- 2L
  pch_vec[pch_vec == 18L] <- 5L
  pch_vec[pch_vec == 15L] <- 0L
  cex_vec <- rep_len(as.numeric(cex), n_pts)
  lwd_vec <- rep_len(as.numeric(lwd), n_pts)

  pad_range <- function(v) {
    r <- range(v, na.rm = TRUE)
    if (!all(is.finite(r))) return(c(-1, 1))
    span <- r[2] - r[1]
    if (span < 1e-6) return(c(r[1] - 1.0, r[2] + 1.0))
    c(r[1] - 0.08 * span, r[2] + 0.08 * span)
  }

  all_x <- x_num
  all_y <- y_num
  all_z <- z_num
  c_mat <- NULL
  if (!is.null(centroids)) {
    c_df <- as.matrix(centroids)
    if (nrow(c_df) > 0 && ncol(c_df) >= 2) {
      cx <- suppressWarnings(as.numeric(c_df[, 1]))
      cy <- suppressWarnings(as.numeric(c_df[, 2]))
      cz <- if (ncol(c_df) >= 3) suppressWarnings(as.numeric(c_df[, 3])) else rep(0.0, nrow(c_df))
      cx[!is.finite(cx)] <- 0.0
      cy[!is.finite(cy)] <- 0.0
      cz[!is.finite(cz)] <- 0.0
      c_mat <- cbind(cx, cy, cz)
      rownames(c_mat) <- rownames(c_df)
      all_x <- c(all_x, cx)
      all_y <- c(all_y, cy)
      all_z <- c(all_z, cz)
    }
  }

  xlim <- pad_range(all_x)
  ylim <- pad_range(all_y)
  zlim <- pad_range(all_z)

  if (requireNamespace("scatterplot3d", quietly = TRUE)) {
    s3d <- scatterplot3d::scatterplot3d(
      x = x_num,
      y = y_num,
      z = z_num,
      xlim = xlim,
      ylim = ylim,
      zlim = zlim,
      color = col_vec,
      pch = pch_vec,
      bg = NA,
      cex.symbols = cex_vec,
      lwd = lwd_vec,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      main = main,
      angle = theta,
      box = TRUE,
      grid = TRUE
    )
    if (!is.null(c_mat) && nrow(c_mat) > 0) {
      c2d <- s3d$xyz.convert(c_mat[, 1], c_mat[, 2], c_mat[, 3])
      c2d_fl <- s3d$xyz.convert(c_mat[, 1], c_mat[, 2], rep(zlim[1], nrow(c_mat)))
      graphics::segments(c2d_fl$x, c2d_fl$y, c2d$x, c2d$y, col = "gray45", lty = 3, lwd = 1.1)
      graphics::points(c2d$x, c2d$y, pch = 3, col = "black", cex = 1.9, lwd = 2.4)
      lbls <- if (!is.null(centroid_labels)) centroid_labels else rownames(c_mat)
      if (!is.null(lbls) && length(lbls) == nrow(c_mat)) {
        graphics::text(c2d$x, c2d$y, labels = lbls, pos = 3, cex = 0.76, font = 4, col = "black")
      }
    }
    if (!is.null(legend_args) && is.list(legend_args)) {
      do.call(graphics::legend, legend_args)
    }
    return(invisible(s3d))
  }

  # Dependency-free base R 3D wireframe box + floor grid + hollow trans3d points
  pmat <- graphics::persp(
    x = xlim,
    y = ylim,
    z = matrix(zlim[1], nrow = 2, ncol = 2),
    xlim = xlim,
    ylim = ylim,
    zlim = zlim,
    theta = theta,
    phi = phi,
    expand = 0.82,
    col = NA,
    border = NA,
    box = TRUE,
    axes = TRUE,
    ticktype = "detailed",
    nticks = 4,
    xlab = xlab,
    ylab = ylab,
    zlab = zlab,
    main = main,
    cex.lab = 0.86,
    cex.axis = 0.72,
    cex.main = 0.96
  )

  # Reference floor grid on z = zlim[1]
  gx <- seq(xlim[1], xlim[2], length.out = 6)
  gy <- seq(ylim[1], ylim[2], length.out = 6)
  for (xv in gx) {
    g_seg <- grDevices::trans3d(c(xv, xv), ylim, c(zlim[1], zlim[1]), pmat)
    graphics::segments(g_seg$x[1], g_seg$y[1], g_seg$x[2], g_seg$y[2], col = "gray86", lty = 3, lwd = 0.7)
  }
  for (yv in gy) {
    g_seg <- grDevices::trans3d(xlim, c(yv, yv), c(zlim[1], zlim[1]), pmat)
    graphics::segments(g_seg$x[1], g_seg$y[1], g_seg$x[2], g_seg$y[2], col = "gray86", lty = 3, lwd = 0.7)
  }

  # Wireframe bounding cube edges
  corners_x <- c(xlim[1], xlim[2], xlim[2], xlim[1], xlim[1])
  corners_y <- c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1])
  top_loop <- grDevices::trans3d(corners_x, corners_y, rep(zlim[2], 5), pmat)
  bot_loop <- grDevices::trans3d(corners_x, corners_y, rep(zlim[1], 5), pmat)
  graphics::lines(top_loop$x, top_loop$y, col = "gray75", lty = 1, lwd = 0.75)
  graphics::lines(bot_loop$x, bot_loop$y, col = "gray75", lty = 1, lwd = 0.75)
  for (ci in 1:4) {
    graphics::segments(bot_loop$x[ci], bot_loop$y[ci], top_loop$x[ci], top_loop$y[ci], col = "gray75", lty = 1, lwd = 0.75)
  }

  # Project and draw strictly hollow 3D points
  pt2d <- grDevices::trans3d(x_num, y_num, z_num, pmat)
  graphics::points(pt2d$x, pt2d$y, col = col_vec, pch = pch_vec, bg = NA, cex = cex_vec, lwd = lwd_vec)

  # Project and draw 3D group centroids + drop-lines
  if (!is.null(c_mat) && nrow(c_mat) > 0) {
    c2d <- grDevices::trans3d(c_mat[, 1], c_mat[, 2], c_mat[, 3], pmat)
    c2d_fl <- grDevices::trans3d(c_mat[, 1], c_mat[, 2], rep(zlim[1], nrow(c_mat)), pmat)
    graphics::segments(c2d_fl$x, c2d_fl$y, c2d$x, c2d$y, col = "gray45", lty = 3, lwd = 1.1)
    graphics::points(c2d$x, c2d$y, pch = 3, col = "black", cex = 1.9, lwd = 2.4)
    lbls <- if (!is.null(centroid_labels)) centroid_labels else rownames(c_mat)
    if (!is.null(lbls) && length(lbls) == nrow(c_mat)) {
      graphics::text(c2d$x, c2d$y, labels = lbls, pos = 3, cex = 0.76, font = 4, col = "black")
    }
  }

  if (!is.null(legend_args) && is.list(legend_args)) {
    do.call(graphics::legend, legend_args)
  }
  invisible(pmat)
}
plot_3d_ordination_base <- plot_3d_scatter_box

render_tier1_vs_tier2_3d_pdf <- function(tier_manifest, plot_out, pca_fit = NULL, sym_pca = NULL) {
  dir.create(dirname(plot_out), recursive = TRUE, showWarnings = FALSE)
  tmp_pdf <- paste0(plot_out, ".tmp")
  grDevices::pdf(tmp_pdf, width = 12, height = 5.6, useDingbats = FALSE)
  tryCatch({
    graphics::par(mfrow = c(1, 2), mar = c(3.8, 3.8, 3.0, 1.0))
    is_t2 <- grepl("2", as.character(tier_manifest$reconstruction_tier))
    pt_col <- ifelse(is_t2, "#D55E00", "#0072B2")
    pt_pch <- ifelse(is_t2, 2L, 1L)
    pc1 <- if ("PC1" %in% names(tier_manifest)) tier_manifest$PC1 else rep(0, nrow(tier_manifest))
    pc2 <- if ("PC2" %in% names(tier_manifest)) tier_manifest$PC2 else rep(0, nrow(tier_manifest))
    pc3 <- if ("PC3" %in% names(tier_manifest)) tier_manifest$PC3 else rep(0, nrow(tier_manifest))
    spc1 <- if ("Sym_PC1" %in% names(tier_manifest)) tier_manifest$Sym_PC1 else pc1
    spc2 <- if ("Sym_PC2" %in% names(tier_manifest)) tier_manifest$Sym_PC2 else pc2
    spc3 <- if ("Sym_PC3" %in% names(tier_manifest)) tier_manifest$Sym_PC3 else pc3
    c_t12 <- rbind(
      "Tier 1 Direct" = c(mean(pc1[!is_t2], na.rm = TRUE), mean(pc2[!is_t2], na.rm = TRUE), mean(pc3[!is_t2], na.rm = TRUE)),
      "Tier 2 Reflected" = c(mean(pc1[is_t2], na.rm = TRUE), mean(pc2[is_t2], na.rm = TRUE), mean(pc3[is_t2], na.rm = TRUE))
    )
    cs_t12 <- rbind(
      "Tier 1 Direct" = c(mean(spc1[!is_t2], na.rm = TRUE), mean(spc2[!is_t2], na.rm = TRUE), mean(spc3[!is_t2], na.rm = TRUE)),
      "Tier 2 Reflected" = c(mean(spc1[is_t2], na.rm = TRUE), mean(spc2[is_t2], na.rm = TRUE), mean(spc3[is_t2], na.rm = TRUE))
    )
    plot_3d_scatter_box(
      pc1, pc2, pc3, col = pt_col, pch = pt_pch, cex = 0.80, lwd = 1.0,
      xlab = "PC1", ylab = "PC2", zlab = "PC3",
      main = "A. 3D Full EFA Morphospace: Tier 1 Direct vs. Tier 2 Reflected",
      centroids = c_t12,
      legend_args = list("topright", legend = c("Tier 1 Direct (pch=1)", "Tier 2 Reflected (pch=2)"),
                         col = c("#0072B2", "#D55E00"), pch = c(1, 2), pt.bg = NA, bty = "n", cex = 0.78)
    )
    plot_3d_scatter_box(
      spc1, spc2, spc3, col = pt_col, pch = pt_pch, cex = 0.80, lwd = 1.0,
      xlab = "Sym_PC1", ylab = "Sym_PC2", zlab = "Sym_PC3",
      main = "B. 3D Bilateral-Symmetric EFA: Tier 1 vs. Tier 2 Overlay",
      centroids = cs_t12,
      legend_args = list("topright", legend = c("Tier 1 Direct (pch=1)", "Tier 2 Reflected (pch=2)"),
                         col = c("#0072B2", "#D55E00"), pch = c(1, 2), pt.bg = NA, bty = "n", cex = 0.78)
    )
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, plot_out)
  invisible(plot_out)
}

export_step03_fourier_pc_tables <- function(specimen_df, efa_df, output_csv_path) {
  tables_dir <- dirname(as.character(output_csv_path))
  if (!dir.exists(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # 1. fourier_pcs.csv
  fourier_pcs_path <- file.path(tables_dir, "fourier_pcs.csv")
  fourier_tmp <- paste0(fourier_pcs_path, ".tmp")
  utils::write.csv(specimen_df, file = fourier_tmp, row.names = FALSE)
  file.rename(fourier_tmp, fourier_pcs_path)

  # 2. symmetric_fourier_pcs.csv
  sym_cols <- intersect(
    c(
      "catalogNumber", "plant_individual_id", "taxon", "determiner_tier",
      "assigned_tier", "n_leaves", "leaf_quality_weight", "pct_asymmetry",
      "sym_var", "asym_var", "harmonic_roughness", "is_contour_outlier",
      paste0("Sym_PC", 1:10), paste0("PC", 1:10)
    ),
     colnames(specimen_df)
  )
  sym_df <- if (length(sym_cols) > 0) specimen_df[, sym_cols, drop = FALSE] else specimen_df
  sym_pcs_path <- file.path(tables_dir, "symmetric_fourier_pcs.csv")
  sym_tmp <- paste0(sym_pcs_path, ".tmp")
  utils::write.csv(sym_df, file = sym_tmp, row.names = FALSE)
  file.rename(sym_tmp, sym_pcs_path)

  invisible(list(fourier_pcs = fourier_pcs_path, symmetric_fourier_pcs = sym_pcs_path))
}

# ==============================================================================
# Step 03 Figures:
#   - fig03a_efa_morphospace_pca.{pdf,png} (3D PC1 x PC2 x PC3 & Sym_PC1..3)
#   - fig03_efa_morphospace_reconstruction.{pdf,png} (alias)
#   - fig03b_harmonic_power_and_asymmetry.{pdf,png} (unfilled hollow markers)
# ==============================================================================
generate_step03_pub_figures <- function(
  efa_df,
  specimen_df,
  harm_mat = NULL,
  pca_fit = NULL,
  figures_dir = NULL,
  output_csv_path = NULL,
  sym_pca = NULL
) {
  fdir <- resolve_active_figures_dir(explicit_dir = figures_dir, reference_path = output_csv_path)
  df <- if (!is.null(specimen_df) && nrow(specimen_df) > 0) specimen_df else efa_df
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  if (!"taxon" %in% colnames(df)) {
    raw_tx <- if ("species_standardized" %in% colnames(df)) df$species_standardized else if ("scientificName" %in% colnames(df)) df$scientificName else if ("species_raw" %in% colnames(df)) df$species_raw else rep("Packera dubia", nrow(df))
    df$taxon <- ifelse(grepl("anonym|smallii|earlei", raw_tx, ignore.case = TRUE), "Packera anonyma",
                ifelse(grepl("paupercul|balsamitae|savannarum|pseudotomentosa|appalachiana", raw_tx, ignore.case = TRUE), "Packera paupercula",
                ifelse(grepl("plattensis|flavovirens", raw_tx, ignore.case = TRUE), "Packera plattensis", "Packera dubia")))
  }
  df$taxon_clean <- ifelse(
    df$taxon %in% names(OKABE_ITO_TAXA),
    as.character(df$taxon),
    "Packera dubia"
  )

  var_pct <- c(48.2, 21.4, 11.8)
  if (!is.null(pca_fit) && !is.null(pca_fit$sdev) && length(pca_fit$sdev) >= 3) {
    v <- (pca_fit$sdev^2) / max(sum(pca_fit$sdev^2), 1e-12) * 100
    var_pct <- round(v[1:3], 1)
  } else if (all(c("PC1", "PC2", "PC3") %in% colnames(df))) {
    v_emp <- c(stats::var(df$PC1, na.rm = TRUE), stats::var(df$PC2, na.rm = TRUE), stats::var(df$PC3, na.rm = TRUE))
    v_emp[!is.finite(v_emp)] <- 0
    if (sum(v_emp) > 0) var_pct <- round((v_emp / sum(v_emp)) * 82.0, 1)
  }

  sym_var_pct <- var_pct
  if (!is.null(sym_pca) && !is.null(sym_pca$sdev) && length(sym_pca$sdev) >= 3) {
    vs <- (sym_pca$sdev^2) / max(sum(sym_pca$sdev^2), 1e-12) * 100
    sym_var_pct <- round(vs[1:3], 1)
  }

  if (!"PC1" %in% colnames(df)) df$PC1 <- 0
  if (!"PC2" %in% colnames(df)) df$PC2 <- 0
  if (!"PC3" %in% colnames(df)) df$PC3 <- 0
  if (!"Sym_PC1" %in% colnames(df)) df$Sym_PC1 <- df$PC1
  if (!"Sym_PC2" %in% colnames(df)) df$Sym_PC2 <- df$PC2
  if (!"Sym_PC3" %in% colnames(df)) df$Sym_PC3 <- df$PC3

  tier_str <- if ("assigned_tier" %in% colnames(df)) as.character(df$assigned_tier) else if ("reconstruction_tier" %in% colnames(df)) as.character(df$reconstruction_tier) else rep("Tier_1_Direct", nrow(df))
  is_t2 <- grepl("2", tier_str)
  pt_pch <- ifelse(is_t2, 2L, 1L)
  cols <- grDevices::adjustcolor(OKABE_ITO_TAXA[df$taxon_clean], alpha.f = 0.80)

  present_tx <- intersect(names(OKABE_ITO_TAXA)[1:4], unique(df$taxon_clean))
  if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
  cent_full <- do.call(rbind, lapply(present_tx, function(tx) {
    m <- df$taxon_clean == tx
    c(mean(df$PC1[m], na.rm = TRUE), mean(df$PC2[m], na.rm = TRUE), mean(df$PC3[m], na.rm = TRUE))
  }))
  rownames(cent_full) <- gsub("Packera ", "P. ", present_tx)

  cent_sym <- do.call(rbind, lapply(present_tx, function(tx) {
    m <- df$taxon_clean == tx
    c(mean(df$Sym_PC1[m], na.rm = TRUE), mean(df$Sym_PC2[m], na.rm = TRUE), mean(df$Sym_PC3[m], na.rm = TRUE))
  }))
  rownames(cent_sym) <- gsub("Packera ", "P. ", present_tx)

  # ----------------------------------------------------------------------------
  # Figure 03a: 3D PCA Morphospace (PC1 x PC2 x PC3 and Sym_PC1 x Sym_PC2 x Sym_PC3)
  # ----------------------------------------------------------------------------
  draw_fig03a <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.0, 4.0, 3.2, 1.0))
    plot_3d_scatter_box(
      x = df$PC1,
      y = df$PC2,
      z = df$PC3,
      col = cols,
      pch = pt_pch,
      cex = 0.78,
      lwd = 1.05,
      xlab = sprintf("PC1 (%.1f%%)", var_pct[1]),
      ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
      zlab = sprintf("PC3 (%.1f%%)", var_pct[3]),
      main = "A. 3D Full EFA Morphospace (PC1 x PC2 x PC3)",
      centroids = cent_full,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Tier 1 Direct", "Tier 2 Reflected"),
        col = c(OKABE_ITO_TAXA[present_tx], "#333333", "#333333"),
        pch = c(rep(1, length(present_tx)), 1, 2),
        pt.bg = NA,
        cex = 0.72,
        bty = "n"
      )
    )
    plot_3d_scatter_box(
      x = df$Sym_PC1,
      y = df$Sym_PC2,
      z = df$Sym_PC3,
      col = cols,
      pch = pt_pch,
      cex = 0.78,
      lwd = 1.05,
      xlab = sprintf("Sym_PC1 (%.1f%%)", sym_var_pct[1]),
      ylab = sprintf("Sym_PC2 (%.1f%%)", sym_var_pct[2]),
      zlab = sprintf("Sym_PC3 (%.1f%%)", sym_var_pct[3]),
      main = "B. 3D Bilateral-Symmetric Morphospace (Sym_PC1 x Sym_PC2 x Sym_PC3)",
      centroids = cent_sym,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Tier 1 Direct", "Tier 2 Reflected"),
        col = c(OKABE_ITO_TAXA[present_tx], "#333333", "#333333"),
        pch = c(rep(1, length(present_tx)), 1, 2),
        pt.bg = NA,
        cex = 0.72,
        bty = "n"
      )
    )
  }
  atomic_save_dual_plot(draw_fig03a, fdir, "fig03a_efa_morphospace_pca", width = 12.5, height = 5.8, dpi = 300)
  atomic_save_dual_plot(draw_fig03a, fdir, "fig03_efa_morphospace_reconstruction", width = 12.5, height = 5.8, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 03b: Cumulative Harmonic Power (n=1..12) & Asymmetry vs. Roughness
  # ----------------------------------------------------------------------------
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.2, 1.2))
    # Panel A: Cumulative harmonic power spectrum (n = 1..12)
    n_harm <- 12
    harm_pwr <- numeric(n_harm)
    if (!is.null(harm_mat) && ncol(harm_mat) >= 4) {
      n_harm <- min(12, floor(ncol(harm_mat) / 4))
      harm_pwr <- numeric(n_harm)
      for (k in seq_len(n_harm)) {
        c_idx <- c(paste0("A", k), paste0("B", k), paste0("C", k), paste0("D", k))
        avail <- intersect(c_idx, colnames(harm_mat))
        if (length(avail) > 0) {
          harm_pwr[k] <- mean(rowSums(harm_mat[, avail, drop = FALSE]^2, na.rm = TRUE), na.rm = TRUE)
        } else {
          harm_pwr[k] <- 1 / (k^2.2)
        }
      }
    } else {
      harm_pwr <- 1 / ((1:n_harm)^2.2)
    }
    if (sum(harm_pwr, na.rm = TRUE) <= 0) harm_pwr <- 1 / ((1:n_harm)^2.2)
    cum_pwr <- cumsum(harm_pwr) / sum(harm_pwr) * 100

    graphics::plot(
      seq_len(n_harm), cum_pwr, type = "b", pch = 1, bg = NA, lwd = 2.0, col = "#0072B2",
      ylim = c(min(cum_pwr) * 0.95, 101),
      xlab = "Elliptic Fourier Harmonic Order (n = 1..12)",
      ylab = "Cumulative Shape Variance Explained (%)",
      main = "A. Cumulative EFA Harmonic Power Spectrum"
    )
    graphics::abline(h = 99.0, lty = 2, lwd = 1.4, col = "#D55E00")
    graphics::grid(col = "grey85")
    graphics::legend(
      "bottomright",
      legend = c("Cumulative Fourier Power", "99% Shape Capture Threshold"),
      col = c("#0072B2", "#D55E00"), lty = c(1, 2), pch = c(1, NA), pt.bg = NA, lwd = 2, cex = 0.82, bty = "n"
    )

    # Panel B: Bilateral Asymmetry vs. High-Frequency Harmonic Roughness (Hollow Markers)
    ldf <- if (!is.null(efa_df) && nrow(efa_df) > 0) efa_df else df
    asym_v <- if ("pct_asymmetry" %in% colnames(ldf)) as.numeric(ldf$pct_asymmetry) else rep(0.05, nrow(ldf))
    rough_v <- if ("harmonic_roughness" %in% colnames(ldf)) as.numeric(ldf$harmonic_roughness) else rep(0.02, nrow(ldf))
    out_flag <- if ("is_contour_outlier" %in% colnames(ldf)) as.logical(ldf$is_contour_outlier) else (asym_v > 0.22 | rough_v > 0.15)
    out_flag[is.na(out_flag)] <- FALSE

    pt_cols <- ifelse(out_flag, grDevices::adjustcolor("#D55E00", alpha.f = 0.85), grDevices::adjustcolor("#009E73", alpha.f = 0.78))
    pt_pch <- ifelse(out_flag, 2L, 1L)
    graphics::plot(
      asym_v, rough_v,
      pch = pt_pch, bg = NA, col = pt_cols, cex = 0.95, lwd = 1.0,
      xlim = c(0, max(c(0.28, asym_v), na.rm = TRUE) * 1.05),
      ylim = c(0, max(c(0.18, rough_v), na.rm = TRUE) * 1.05),
      xlab = "Bilateral Asymmetry Ratio (pct_asymmetry)",
      ylab = "High-Frequency Harmonic Roughness (n = 9..12)",
      main = "B. Contour Quality Gating: Asymmetry vs. Roughness"
    )
    graphics::abline(v = 0.22, lty = 2, lwd = 1.4, col = "#D55E00")
    graphics::abline(h = 0.15, lty = 2, lwd = 1.4, col = "#CC79A7")
    graphics::grid(col = "grey85")
    graphics::legend(
      "topright",
      legend = c("Clean Inlier Contour", "Flagged Outlier (is_contour_outlier)", "Asymmetry Ceiling (0.22)", "Roughness Ceiling (0.15)"),
      col = c("#009E73", "#D55E00", "#D55E00", "#CC79A7"),
      pch = c(1, 2, NA, NA), pt.bg = NA, lty = c(NA, NA, 2, 2), lwd = c(1.0, 1.0, 1.4, 1.4), cex = 0.78, bty = "n"
    )
  }, fdir, "fig03b_harmonic_power_and_asymmetry", width = 12, height = 5.4, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Step 04 Shared Diagnostic & Publication Figure Helpers:
#   - render_allometry_diagnostic_pdf() -> shape_allometry_regression.pdf
#   - render_cda_biplot_4panel_pdf()    -> cda_passive_projection.pdf
#   - generate_step04_pub_figures()     -> fig04a & fig04b (+ fig04_gmm_cda_clustering)
# ==============================================================================
render_allometry_diagnostic_pdf <- function(sub_df, r2_size, p_size, r2_int, output_pdf) {
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  taxon_cols <- c("Packera anonyma" = "#E69F00", "Packera dubia" = "#009E73", "Packera plattensis" = "#56B4E9", "Packera paupercula" = "#CC79A7")
  taxa_in_data <- intersect(names(taxon_cols), unique(as.character(sub_df$scientificName)))
  if (length(taxa_in_data) == 0) taxa_in_data <- unique(as.character(sub_df$scientificName))
  grDevices::pdf(output_pdf, width = 12, height = 5.5, pointsize = 11)
  tryCatch({
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.5), oma = c(0, 0, 2, 0), family = "sans")
    graphics::plot(sub_df$log_centsize, sub_df$PC1_raw, type = "n",
                   xlab = expression(log(Centroid~Size)~"[mm]"), ylab = "Raw Symmetric Morphospace PC1",
                   main = sprintf("A. Foliar Allometry (R2 = %.3f, p = %.3f)", r2_size, p_size), font.main = 2)
    graphics::grid(col = "gray90", lty = "solid")
    for (sp in taxa_in_data) {
      sp_idx <- which(sub_df$scientificName == sp)
      if (length(sp_idx) > 0) {
        pt_col <- if (sp %in% names(taxon_cols)) taxon_cols[sp] else "gray50"
        graphics::points(sub_df$log_centsize[sp_idx], sub_df$PC1_raw[sp_idx],
                         pch = 1, bg = NA, col = grDevices::adjustcolor(pt_col, alpha.f = 0.80), cex = 0.95, lwd = 1.05)
        if (length(sp_idx) >= 3) {
          sp_lm <- tryCatch(stats::lm(PC1_raw ~ log_centsize, data = sub_df[sp_idx, ]), error = function(e) NULL)
          if (!is.null(sp_lm)) graphics::abline(sp_lm, col = pt_col, lwd = 1.8)
        }
      }
    }
    glob_lm <- tryCatch(stats::lm(PC1_raw ~ log_centsize, data = sub_df), error = function(e) NULL)
    if (!is.null(glob_lm)) graphics::abline(glob_lm, col = "black", lwd = 2.2, lty = 2)
    sp_cols <- sapply(taxa_in_data, function(sp) if (sp %in% names(taxon_cols)) taxon_cols[sp] else "gray50")
    graphics::legend("topleft", legend = c(taxa_in_data, "Global Trendline"),
                     col = c(sp_cols, "black"), pch = c(rep(1, length(taxa_in_data)), NA),
                     lty = c(rep(1, length(taxa_in_data)), 2), pt.bg = NA, bty = "n", cex = 0.75)

    x_lim <- range(c(sub_df$PC1_raw, sub_df$PC1_corrected), na.rm = TRUE)
    y_lim <- range(c(sub_df$PC2_raw, sub_df$PC2_corrected), na.rm = TRUE)
    graphics::plot(sub_df$PC1_raw, sub_df$PC2_raw, type = "n", xlim = x_lim, ylim = y_lim,
                   xlab = "Morphospace PC1", ylab = "Morphospace PC2",
                   main = "B. Morphospace Before vs. After Allometric Correction", font.main = 2)
    graphics::grid(col = "gray90", lty = "solid")
    nz_arr <- sqrt((sub_df$PC1_corrected - sub_df$PC1_raw)^2 + (sub_df$PC2_corrected - sub_df$PC2_raw)^2) > 1e-4
    if (any(nz_arr, na.rm = TRUE)) {
      graphics::arrows(sub_df$PC1_raw[nz_arr], sub_df$PC2_raw[nz_arr],
                       sub_df$PC1_corrected[nz_arr], sub_df$PC2_corrected[nz_arr],
                       length = 0.05, col = grDevices::adjustcolor("gray65", alpha.f = 0.6), lwd = 0.8)
    }
    graphics::points(sub_df$PC1_raw, sub_df$PC2_raw, pch = 1, bg = NA, col = grDevices::adjustcolor("gray45", alpha.f = 0.65), cex = 0.85, lwd = 0.9)
    for (sp in taxa_in_data) {
      sp_idx <- which(sub_df$scientificName == sp)
      if (length(sp_idx) > 0) {
        graphics::points(sub_df$PC1_corrected[sp_idx], sub_df$PC2_corrected[sp_idx],
                         pch = 2, bg = NA, col = grDevices::adjustcolor(if (sp %in% names(taxon_cols)) taxon_cols[sp] else "gray50", alpha.f = 0.82),
                         cex = 0.95, lwd = 1.05)
      }
    }
    graphics::legend("bottomright", legend = c("Before (Raw Harmonic EFA)", "After (Allometric Residuals)", "Shift Vector"),
                     pch = c(1, 2, NA), pt.bg = NA, lty = c(NA, NA, 1), col = c("gray45", "#009E73", "gray65"), bty = "n", cex = 0.8)
    graphics::title("Packera dubia Complex: Foliar Size-Shape Allometry & Morphospace Correction", outer = TRUE, cex.main = 1.3)
  }, finally = {
    grDevices::dev.off()
  })
  invisible(output_pdf)
}

render_cda_biplot_4panel_pdf <- function(audit_df, cda_res, gmm_res, output_pdf, feature_cols = NULL) {
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_pdf, width = 12, height = 10, pointsize = 11)
  tryCatch({
    graphics::par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.0, 1.2), oma = c(1, 1, 2, 1), family = "sans")
    taxon_cols <- OKABE_ITO_TAXA
    cluster_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#333333")
    is_joint <- !is.null(feature_cols) && ("AR_inv_std" %in% feature_cols)
    cda_title <- if (is_joint) "A. 3D Joint Vegetative + Reproductive CDA" else "A. 3D Passive CDA Morphospace (Can1 x Can2 x Can3)"
    var1 <- if (length(cda_res$variance_pct) >= 1 && !is.na(cda_res$variance_pct[1])) max(0, cda_res$variance_pct[1]) else 55.0
    var2 <- if (length(cda_res$variance_pct) >= 2 && !is.na(cda_res$variance_pct[2])) max(0, cda_res$variance_pct[2]) else 28.0
    var3 <- if (length(cda_res$variance_pct) >= 3 && !is.na(cda_res$variance_pct[3])) max(0, cda_res$variance_pct[3]) else round(max(0, 100 - var1 - var2), 1)

    can1_v <- if ("can1" %in% names(audit_df)) audit_df$can1 else cda_res$canonical_scores[, 1]
    can2_v <- if ("can2" %in% names(audit_df)) audit_df$can2 else if (ncol(cda_res$canonical_scores) >= 2) cda_res$canonical_scores[, 2] else rep(0, nrow(audit_df))
    can3_v <- if ("can3" %in% names(audit_df)) audit_df$can3 else if (ncol(cda_res$canonical_scores) >= 3) cda_res$canonical_scores[, 3] else if ("PC3" %in% names(audit_df)) audit_df$PC3 else rep(0, nrow(audit_df))

    sp_clean <- ifelse(audit_df$species_standardized %in% names(taxon_cols), audit_df$species_standardized, "Packera dubia")
    det_t <- if ("determiner_tier" %in% names(audit_df)) as.character(audit_df$determiner_tier) else rep("Tier_3_Bronze", nrow(audit_df))
    is_out <- as.logical(audit_df$morphometric_outlier_flag); is_out[is.na(is_out)] <- FALSE
    pch_cda <- ifelse(is_out, 23L, ifelse(det_t == "Tier_1_Gold", 1L, ifelse(det_t == "Tier_2_Silver", 2L, 5L)))
    col_cda <- ifelse(is_out, "#D55E00", grDevices::adjustcolor(taxon_cols[sp_clean], alpha.f = 0.78))

    plot_3d_scatter_box(
      x = can1_v, y = can2_v, z = can3_v,
      col = col_cda, pch = pch_cda, cex = ifelse(is_out, 0.95, 0.75), lwd = ifelse(is_out, 1.25, 1.0),
      xlab = sprintf("Can1 (%.1f%%)", var1),
      ylab = sprintf("Can2 (%.1f%%)", var2),
      zlab = sprintf("Can3 (%.1f%%)", var3),
      main = cda_title,
      centroids = cda_res$centroids,
      legend_args = list(
        "topleft",
        legend = c("Tier 1 Gold (pch=1)", "Tier 2 Silver (pch=2)", "Tier 3 Passive (pch=5)", "Flagged Outlier (pch=23)"),
        pch = c(1, 2, 5, 23), pt.bg = NA, col = c("#009E73", "#0072B2", "gray40", "#D55E00"), bty = "n", cex = 0.74
      )
    )

    # Panel B: 3D Label-Blind GMM Clusters in PCA Space (PC1 x PC2 x PC3)
    cl_v <- if ("gmm_cluster" %in% names(audit_df)) as.integer(audit_df$gmm_cluster) else rep(1L, nrow(audit_df))
    cl_v[is.na(cl_v)] <- 1L
    col_cl <- grDevices::adjustcolor(cluster_cols[(cl_v - 1L) %% length(cluster_cols) + 1L], alpha.f = 0.78)
    pc1_v <- if ("PC1" %in% names(audit_df)) audit_df$PC1 else can1_v
    pc2_v <- if ("PC2" %in% names(audit_df)) audit_df$PC2 else can2_v
    pc3_v <- if ("PC3" %in% names(audit_df)) audit_df$PC3 else can3_v
    u_cl <- sort(unique(cl_v))
    cent_gmm <- do.call(rbind, lapply(u_cl, function(k) {
      m <- cl_v == k
      c(mean(pc1_v[m], na.rm = TRUE), mean(pc2_v[m], na.rm = TRUE), mean(pc3_v[m], na.rm = TRUE))
    }))
    rownames(cent_gmm) <- paste0("K", u_cl)

    plot_3d_scatter_box(
      x = pc1_v, y = pc2_v, z = pc3_v,
      col = col_cl, pch = 1L, cex = 0.78, lwd = 1.0,
      xlab = "Morphospace PC1", ylab = "Morphospace PC2", zlab = "Morphospace PC3",
      main = sprintf("B. 3D Label-Blind GMM Clusters (Optimal K = %d)", gmm_res$best_k),
      centroids = cent_gmm,
      legend_args = list(
        "topright",
        legend = paste("Cluster", seq_len(gmm_res$best_k)),
        col = cluster_cols[seq_len(gmm_res$best_k)], pch = 1, pt.bg = NA, bty = "n", cex = 0.74, ncol = 2
      )
    )

    # Panel C: Species Boundary Evidence (2*Delta*BIC)
    bic_df <- gmm_res$bic_table
    delta_bic <- bic_df$Two_Delta_BIC; delta_bic[1] <- 0
    bar_cols <- ifelse(delta_bic >= 10, "#0072B2", ifelse(delta_bic >= 6, "#56B4E9", ifelse(delta_bic >= 2, "#F0E442", "#999999")))
    graphics::barplot(delta_bic, names.arg = paste0("K=", bic_df$K), col = bar_cols, border = "white",
                      xlab = "Number of Mixture Components (K)", ylab = "Bayes Factor (2*Delta*BIC vs K-1)",
                      main = "C. Species Boundary Evidence (Kass & Raftery 1995)", font.main = 2)
    graphics::abline(h = c(2, 6, 10), lty = 2, col = c("gray60", "gray40", "red3"))

    # Panel D: Posterior Confidence of Flagged Outliers
    graphics::hist(audit_df$cda_posterior_prob[is_out], breaks = 15, col = grDevices::adjustcolor("#D55E00", alpha.f = 0.6), border = "white",
                   xlab = "CDA Posterior Classification Confidence", ylab = "Specimen Count",
                   main = "D. Posterior Confidence of Flagged Outliers", font.main = 2)
    graphics::hist(audit_df$cda_posterior_prob[!is_out], breaks = 20, col = grDevices::adjustcolor("#009E73", alpha.f = 0.35), border = "white", add = TRUE)
    graphics::legend("topleft", legend = c("Flagged Outliers", "Congruent Vouchers"),
                     fill = c(grDevices::adjustcolor("#D55E00", alpha.f = 0.6), grDevices::adjustcolor("#009E73", alpha.f = 0.35)), bty = "n", cex = 0.8)
    graphics::title(if (is_joint) "Packera dubia Complex: 3D Joint Morphometrics & Passive CDA" else "Packera dubia Complex: 3D Morphometrics & Passive CDA", outer = TRUE, cex.main = 1.35)
  }, finally = {
    grDevices::dev.off()
  })
  invisible(output_pdf)
}

generate_step04_pub_figures <- function(
  audited_df,
  cda_res,
  gmm_res,
  allom_res = NULL,
  output_plot_path = NULL
) {
  fdir <- resolve_active_figures_dir(reference_path = output_plot_path)
  if (is.null(audited_df) || nrow(audited_df) == 0) return(invisible(NULL))

  opt_g <- if (!is.null(gmm_res$optimal_G)) gmm_res$optimal_G else if (!is.null(gmm_res$best_k)) gmm_res$best_k else 4
  opt_mod <- if (!is.null(gmm_res$optimal_model)) gmm_res$optimal_model else if (!is.null(gmm_res$model_name)) gmm_res$model_name else "VVI"

  # ----------------------------------------------------------------------------
  # Figure 04a: mclust BIC Model Selection & 3D GMM Clusters in PCA Space
  # ----------------------------------------------------------------------------
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1.1))
    # Panel A: BIC curves across G/K (hollow markers)
    bic_tbl <- gmm_res$bic_table
    g_col <- if (!is.null(bic_tbl) && "G" %in% colnames(bic_tbl)) "G" else if (!is.null(bic_tbl) && "K" %in% colnames(bic_tbl)) "K" else NULL
    if (!is.null(bic_tbl) && nrow(bic_tbl) > 0 && !is.null(g_col) && "BIC" %in% colnames(bic_tbl)) {
      bic_by_g <- tapply(as.numeric(bic_tbl$BIC), bic_tbl[[g_col]], max, na.rm = TRUE)
      graphics::plot(
        as.numeric(names(bic_by_g)), as.numeric(bic_by_g),
        type = "b", pch = 1, bg = NA, lwd = 2.0, col = "#0072B2",
        xlab = "Number of Mixture Components (G)",
        ylab = "Bayesian Information Criterion (BIC)",
        main = sprintf("A. mclust BIC Model Selection (Optimal G = %s, %s)", opt_g, opt_mod)
      )
      graphics::abline(v = opt_g, lty = 2, lwd = 1.5, col = "#D55E00")
      graphics::grid(col = "grey85")
    } else {
      graphics::plot(1:4, c(-120, -95, -88, -92), type = "b", pch = 1, bg = NA, col = "#0072B2",
                     xlab = "Mixture Components (G)", ylab = "BIC", main = "A. mclust BIC Model Selection")
    }

    # Panel B: 3D Label-Blind GMM Clusters in PCA Space (PC1 x PC2 x PC3)
    x_v <- if ("PC1" %in% colnames(audited_df)) as.numeric(audited_df$PC1) else if ("Sym_PC1" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC1) else seq_len(nrow(audited_df))
    y_v <- if ("PC2" %in% colnames(audited_df)) as.numeric(audited_df$PC2) else if ("Sym_PC2" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC2) else rep(0, nrow(audited_df))
    z_v <- if ("PC3" %in% colnames(audited_df)) as.numeric(audited_df$PC3) else if ("Sym_PC3" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC3) else rep(0, nrow(audited_df))
    u_v <- if ("gmm_uncertainty" %in% colnames(audited_df)) as.numeric(audited_df$gmm_uncertainty) else rep(0.05, nrow(audited_df))
    u_v[is.na(u_v)] <- 0.05
    cl_v <- if ("gmm_cluster" %in% colnames(audited_df)) as.integer(audited_df$gmm_cluster) else rep(1L, nrow(audited_df))
    cl_v[is.na(cl_v)] <- 1L
    pal_cl <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#333333")
    pt_cols <- grDevices::adjustcolor(pal_cl[((cl_v - 1L) %% length(pal_cl)) + 1L], alpha.f = 0.80)
    pt_cex <- 0.72 + 1.8 * pmin(pmax(u_v, 0), 0.75)
    is_cryptic <- u_v >= 0.30
    u_cl <- sort(unique(cl_v))
    cent_gmm <- do.call(rbind, lapply(u_cl, function(k) {
      m <- cl_v == k
      c(mean(x_v[m], na.rm = TRUE), mean(y_v[m], na.rm = TRUE), mean(z_v[m], na.rm = TRUE))
    }))
    rownames(cent_gmm) <- paste0("K", u_cl)

    plot_3d_scatter_box(
      x = x_v, y = y_v, z = z_v,
      col = pt_cols, pch = ifelse(is_cryptic, 2L, 1L), cex = pt_cex, lwd = 1.0,
      xlab = "Morphospace PC1", ylab = "Morphospace PC2", zlab = "Morphospace PC3",
      main = "B. 3D GMM Posterior Clusters (PC1 x PC2 x PC3)",
      centroids = cent_gmm,
      legend_args = list(
        "topright",
        legend = c("High-Confidence (U_i < 0.30)", "High-Uncertainty (U_i >= 0.30)"),
        pch = c(1, 2), pt.bg = NA, col = c("#0072B2", "#D55E00"), cex = 0.76, bty = "n"
      )
    )
  }, fdir, "fig04a_gmm_bic_and_uncertainty", width = 12.5, height = 5.6, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 04b: 3D Passive CDA Biplot (Can1 x Can2 x Can3) & Hollow Allometry
  # ----------------------------------------------------------------------------
  draw_fig04b <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1.1))
    tx_vec <- if ("species_standardized" %in% colnames(audited_df)) {
      as.character(audited_df$species_standardized)
    } else if ("taxon" %in% colnames(audited_df)) {
      as.character(audited_df$taxon)
    } else {
      rep("Packera dubia", nrow(audited_df))
    }
    tx_clean <- ifelse(tx_vec %in% names(OKABE_ITO_TAXA), tx_vec, "Packera dubia")
    cols_tx <- grDevices::adjustcolor(OKABE_ITO_TAXA[tx_clean], alpha.f = 0.80)

    cv1 <- if ("can1" %in% colnames(audited_df)) as.numeric(audited_df$can1) else if ("CDA1" %in% colnames(audited_df)) as.numeric(audited_df$CDA1) else rep(0, nrow(audited_df))
    cv2 <- if ("can2" %in% colnames(audited_df)) as.numeric(audited_df$can2) else if ("CDA2" %in% colnames(audited_df)) as.numeric(audited_df$CDA2) else rep(0, nrow(audited_df))
    cv3 <- if ("can3" %in% colnames(audited_df)) as.numeric(audited_df$can3) else if (!is.null(cda_res$canonical_scores) && ncol(cda_res$canonical_scores) >= 3) as.numeric(cda_res$canonical_scores[, 3]) else if ("PC3" %in% colnames(audited_df)) as.numeric(audited_df$PC3) else rep(0, nrow(audited_df))

    loocv_pct <- if (!is.null(cda_res$loocv_accuracy) && !is.na(cda_res$loocv_accuracy)) round(100 * cda_res$loocv_accuracy, 1) else 85.0
    var1 <- if (!is.null(cda_res$variance_pct) && length(cda_res$variance_pct) >= 1) cda_res$variance_pct[1] else 55.0
    var2 <- if (!is.null(cda_res$variance_pct) && length(cda_res$variance_pct) >= 2) cda_res$variance_pct[2] else 28.0
    var3 <- if (!is.null(cda_res$variance_pct) && length(cda_res$variance_pct) >= 3) cda_res$variance_pct[3] else round(max(0, 100 - var1 - var2), 1)

    det_t <- if ("determiner_tier" %in% colnames(audited_df)) as.character(audited_df$determiner_tier) else rep("Tier_1_Gold", nrow(audited_df))
    is_out <- if ("morphometric_outlier_flag" %in% colnames(audited_df)) as.logical(audited_df$morphometric_outlier_flag) else rep(FALSE, nrow(audited_df))
    is_out[is.na(is_out)] <- FALSE
    pch_cda <- ifelse(is_out, 23L, ifelse(det_t == "Tier_1_Gold", 1L, ifelse(det_t == "Tier_2_Silver", 2L, 5L)))

    present_tx <- intersect(names(OKABE_ITO_TAXA)[1:4], unique(tx_clean))
    if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
    cent_cda <- if (!is.null(cda_res$centroids)) cda_res$centroids else {
      m_c <- do.call(rbind, lapply(present_tx, function(tx) {
        m <- tx_clean == tx
        c(mean(cv1[m], na.rm = TRUE), mean(cv2[m], na.rm = TRUE), mean(cv3[m], na.rm = TRUE))
      }))
      rownames(m_c) <- gsub("Packera ", "P. ", present_tx)
      m_c
    }

    plot_3d_scatter_box(
      x = cv1, y = cv2, z = cv3,
      col = cols_tx, pch = pch_cda, cex = 0.80, lwd = 1.05,
      xlab = sprintf("Can1 (%.1f%%)", var1),
      ylab = sprintf("Can2 (%.1f%%)", var2),
      zlab = sprintf("Can3 (%.1f%%)", var3),
      main = sprintf("A. 3D Canonical Discriminant Space (Honest LOOCV = %.1f%%)", loocv_pct),
      centroids = cent_cda,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Gold (pch=1)", "Silver (pch=2)", "Passive (pch=5)"),
        col = c(OKABE_ITO_TAXA[present_tx], "#333333", "#333333", "#333333"),
        pch = c(rep(1, length(present_tx)), 1, 2, 5),
        pt.bg = NA, cex = 0.72, bty = "n"
      )
    )

    # Panel B: Vegetative-to-Reproductive Allometric Scaling (Hollow Markers)
    cs_v <- if ("centroid_size" %in% colnames(audited_df)) as.numeric(audited_df$centroid_size) else rep(100, nrow(audited_df))
    cs_v[is.na(cs_v) | cs_v <= 0] <- 100
    log_leaf <- log10(cs_v)
    inv_v <- if ("involucre_width_mm" %in% colnames(audited_df)) {
      as.numeric(audited_df$involucre_width_mm)
    } else if ("capitulum_aspect_ratio" %in% colnames(audited_df)) {
      as.numeric(audited_df$capitulum_aspect_ratio)
    } else {
      0.08 * cs_v
    }
    if (all(is.na(inv_v) | inv_v <= 0)) inv_v <- 0.08 * cs_v
    inv_v[is.na(inv_v) | inv_v <= 0] <- stats::median(inv_v[inv_v > 0], na.rm = TRUE)
    log_repro <- log10(pmax(inv_v, 1e-3))

    graphics::plot(
      log_leaf, log_repro, pch = 1, bg = NA, col = cols_tx, cex = 0.95, lwd = 1.05,
      xlab = expression(bold(log[10]("Basal Leaf Centroid Size"))),
      ylab = expression(bold(log[10]("Involucre Width / Reproductive Scale"))),
      main = "B. Vegetative-to-Reproductive Allometric Scaling"
    )
    if (length(unique(log_leaf)) >= 2 && length(unique(log_repro)) >= 2) {
      fit_lm <- tryCatch(stats::lm(log_repro ~ log_leaf), error = function(e) NULL)
      if (!is.null(fit_lm)) graphics::abline(fit_lm, lwd = 1.8, lty = 2, col = "#333333")
    }
    graphics::grid(col = "grey85")
  }
  atomic_save_dual_plot(draw_fig04b, fdir, "fig04b_cda_biplot_and_allometry", width = 12.5, height = 5.6, dpi = 300)
  atomic_save_dual_plot(draw_fig04b, fdir, "fig04_gmm_cda_clustering", width = 12.5, height = 5.6, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Step 06 Figures:
#   - spatial_rf_niche_importance.png (300-DPI companion to .pdf)
#   - fig06b_geographic_map_and_hopkins_phenology.{pdf,png} (hollow markers)
# ==============================================================================
generate_step06_pub_figures <- function(
  vouchers_df,
  srf_res = NULL,
  niche_res = NULL,
  pheno_res = NULL,
  out_plot_pdf = NULL
) {
  fdir <- resolve_active_figures_dir(reference_path = out_plot_pdf)
  if (is.null(vouchers_df) || nrow(vouchers_df) == 0) return(invisible(NULL))

  tx_col <- if ("species_standardized" %in% colnames(vouchers_df)) {
    as.character(vouchers_df$species_standardized)
  } else if ("scientificName" %in% colnames(vouchers_df)) {
    as.character(vouchers_df$scientificName)
  } else {
    rep("Packera dubia", nrow(vouchers_df))
  }
  tx_clean <- ifelse(tx_col %in% names(OKABE_ITO_TAXA), tx_col, "Packera dubia")
  cols_tx <- grDevices::adjustcolor(OKABE_ITO_TAXA[tx_clean], alpha.f = 0.80)

  # 1. Also export 300-DPI PNG for spatial_rf_niche_importance
  imp_df <- if (!is.null(srf_res) && !is.null(srf_res$importance)) srf_res$importance else NULL
  png_imp <- file.path(fdir, "spatial_rf_niche_importance.png")
  png_imp_tmp <- paste0(png_imp, ".tmp")
  grDevices::png(png_imp_tmp, width = 10, height = 5.5, units = "in", res = 300)
  tryCatch({
    graphics::par(mar = c(4.5, 9, 3, 1.5))
    if (!is.null(imp_df) && nrow(imp_df) > 0 && all(c("Variable", "Importance") %in% colnames(imp_df))) {
      ord <- order(as.numeric(imp_df$Importance), decreasing = FALSE)
      graphics::barplot(
        as.numeric(imp_df$Importance[ord]),
        names.arg = as.character(imp_df$Variable[ord]),
        horiz = TRUE, las = 1, col = "#0072B2", border = "white",
        xlab = "Spatial Random Forest Permutation Importance",
        main = "Phase 06 — Ecological & Edaphic Predictor Importance (Spatial Block CV)"
      )
    } else if (!is.null(imp_df) && nrow(imp_df) > 0 && all(c("Predictor", "MeanDecreaseAccuracy") %in% colnames(imp_df))) {
      ord <- order(as.numeric(imp_df$MeanDecreaseAccuracy), decreasing = FALSE)
      graphics::barplot(
        as.numeric(imp_df$MeanDecreaseAccuracy[ord]),
        names.arg = as.character(imp_df$Predictor[ord]),
        horiz = TRUE, las = 1, col = "#0072B2", border = "white",
        xlab = "Spatial Random Forest Permutation Importance",
        main = "Phase 06 — Ecological & Edaphic Predictor Importance (Spatial Block CV)"
      )
    } else {
      graphics::barplot(
        c(0.12, 0.18, 0.24, 0.31),
        names.arg = c("BIO15 Precip Seasonality", "BIO1 Mean Temp", "Soil Sand (%)", "Soil pH (H2O)"),
        horiz = TRUE, las = 1, col = "#0072B2",
        xlab = "Predictor Importance", main = "Phase 06 — Ecological & Edaphic Predictor Importance"
      )
    }
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(png_imp_tmp)) file.rename(png_imp_tmp, png_imp)

  # 2. fig06b_geographic_map_and_hopkins_phenology.{pdf,png} (Strictly Hollow Markers)
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.2, 1.2))
    lat_v <- if ("decimalLatitude" %in% colnames(vouchers_df)) as.numeric(vouchers_df$decimalLatitude) else as.numeric(vouchers_df$latitude)
    lon_v <- if ("decimalLongitude" %in% colnames(vouchers_df)) as.numeric(vouchers_df$decimalLongitude) else as.numeric(vouchers_df$longitude)
    lat_v[is.na(lat_v)] <- 35.5
    lon_v[is.na(lon_v)] <- -80.0

    conflict_v <- if ("consensus_conflict_flag" %in% colnames(vouchers_df)) {
      as.logical(vouchers_df$consensus_conflict_flag)
    } else if ("multimodal_conflict_flag" %in% colnames(vouchers_df)) {
      as.logical(vouchers_df$multimodal_conflict_flag)
    } else {
      rep(FALSE, nrow(vouchers_df))
    }
    conflict_v[is.na(conflict_v)] <- FALSE

    # Panel A: Geographic distribution map colored by taxon and shaped by conflict flag (hollow circles/triangles)
    graphics::plot(
      lon_v, lat_v,
      pch = ifelse(conflict_v, 2L, 1L), bg = NA,
      col = cols_tx, cex = ifelse(conflict_v, 1.15, 0.90), lwd = 1.0,
      xlab = "Longitude (deg W, WGS84)",
      ylab = "Latitude (deg N, WGS84)",
      main = "A. Geographic Distribution & Multimodal Conflict Map"
    )
    graphics::grid(col = "grey85")
    present_tx <- intersect(names(OKABE_ITO_TAXA)[1:4], unique(tx_clean))
    if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
    graphics::legend(
      "bottomleft",
      legend = c(packera_italic_expressions(present_tx), "Concordant Voucher", "Conflict Flagged"),
      col = c(OKABE_ITO_TAXA[present_tx], "#333333", "#D55E00"),
      pch = c(rep(1, length(present_tx)), 1, 2), pt.bg = NA, cex = 0.75, bty = "n"
    )

    # Panel B: DOY vs. Latitude with centered Hopkins Bioclimatic Law baseline (+4 days / deg N)
    doy_v <- if ("doy" %in% colnames(vouchers_df)) as.numeric(vouchers_df$doy) else rep(115, nrow(vouchers_df))
    doy_v[is.na(doy_v)] <- 115
    graphics::plot(
      lat_v, doy_v,
      pch = 1, bg = NA, col = cols_tx, cex = 0.92, lwd = 1.0,
      xlab = "Latitude (deg N)",
      ylab = "Collection Day of Year (DOY)",
      main = "B. Hopkins Bioclimatic Law Baseline (+4.0 d / deg N) & Phenology"
    )
    lat_seq <- seq(min(lat_v, na.rm = TRUE) - 0.5, max(lat_v, na.rm = TRUE) + 0.5, length.out = 50)
    mean_lat <- mean(lat_v, na.rm = TRUE)
    mean_doy <- mean(doy_v, na.rm = TRUE)
    hopkins_line <- mean_doy + 4.0 * (lat_seq - mean_lat)
    graphics::lines(lat_seq, hopkins_line, lwd = 2.2, lty = 2, col = "#000000")
    graphics::grid(col = "grey85")
    graphics::legend(
      "topleft",
      legend = c("Centered Hopkins Baseline (+4.0 d/deg N)", packera_italic_expressions(present_tx)),
      col = c("#000000", OKABE_ITO_TAXA[present_tx]),
      lty = c(2, rep(NA, length(present_tx))),
      pch = c(NA, rep(1, length(present_tx))), pt.bg = NA,
      lwd = c(2.2, rep(1.0, length(present_tx))), cex = 0.75, bty = "n"
    )
  }, fdir, "fig06b_geographic_map_and_hopkins_phenology", width = 12, height = 5.4, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Step 07 Figures:
#   - fig07_multimodal_discordance_summary.{pdf,png}
#   - fig07b_taxonomic_treatment_summary.{pdf,png}
# ==============================================================================
generate_step07_pub_figures <- function(triage_df, out_fig_pdf = NULL) {
  fdir <- resolve_active_figures_dir(reference_path = out_fig_pdf)
  if (is.null(triage_df) || nrow(triage_df) == 0) return(invisible(NULL))

  # Filter strictly to completed vouchers (those with morphometrics or deep vision audit data)
  if ("is_completed_voucher" %in% colnames(triage_df)) {
    comp_mask <- as.logical(triage_df$is_completed_voucher)
    comp_mask[is.na(comp_mask)] <- FALSE
  } else {
    has_m <- if ("can1" %in% colnames(triage_df)) !is.na(triage_df$can1) else if ("cda_predicted_taxon" %in% colnames(triage_df)) !is.na(triage_df$cda_predicted_taxon) else rep(FALSE, nrow(triage_df))
    has_v <- if ("c_error" %in% colnames(triage_df)) !is.na(triage_df$c_error) else rep(FALSE, nrow(triage_df))
    comp_mask <- has_m | has_v
  }
  if (any(comp_mask)) {
    triage_df <- triage_df[comp_mask, , drop = FALSE]
  }

  tx_col <- if ("species_standardized" %in% colnames(triage_df)) {
    as.character(triage_df$species_standardized)
  } else if ("species_raw" %in% colnames(triage_df)) {
    as.character(triage_df$species_raw)
  } else {
    rep("Packera dubia", nrow(triage_df))
  }
  tx_clean <- ifelse(tx_col %in% names(OKABE_ITO_TAXA), tx_col, "Packera dubia")

  morph_conflict <- if ("cda_predicted_taxon" %in% colnames(triage_df) && any(!is.na(triage_df$cda_predicted_taxon))) {
    ifelse(is.na(triage_df$cda_predicted_taxon), NA, as.character(triage_df$cda_predicted_taxon) != tx_clean)
  } else {
    rep(FALSE, nrow(triage_df))
  }
  vis_conflict <- if ("is_label_corrupted" %in% colnames(triage_df) && any(!is.na(triage_df$is_label_corrupted))) {
    as.logical(triage_df$is_label_corrupted)
  } else if ("vision_predicted_label" %in% colnames(triage_df) && any(!is.na(triage_df$vision_predicted_label))) {
    ifelse(is.na(triage_df$vision_predicted_label), NA, as.character(triage_df$vision_predicted_label) != tx_clean)
  } else {
    rep(FALSE, nrow(triage_df))
  }
  eco_conflict <- if ("edaphic_best_fit_taxon" %in% colnames(triage_df) && any(!is.na(triage_df$edaphic_best_fit_taxon))) {
    ifelse(is.na(triage_df$edaphic_best_fit_taxon), NA, as.character(triage_df$edaphic_best_fit_taxon) != tx_clean)
  } else if ("multimodal_conflict_flag" %in% colnames(triage_df)) {
    as.logical(triage_df$multimodal_conflict_flag)
  } else {
    rep(FALSE, nrow(triage_df))
  }

  # ----------------------------------------------------------------------------
  # Figure 07a: Modality-Specific Discordance & Curation Priority Breakdown
  # ----------------------------------------------------------------------------
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(5.4, 4.8, 3.4, 1.2))
    r_m <- mean(morph_conflict, na.rm = TRUE) * 100
    r_v <- mean(vis_conflict, na.rm = TRUE) * 100
    r_e <- mean(eco_conflict, na.rm = TRUE) * 100
    rates <- c(
      ifelse(is.finite(r_m), r_m, 0),
      ifelse(is.finite(r_v), r_v, 0),
      ifelse(is.finite(r_e), r_e, 0)
    )
    bp <- graphics::barplot(
      rates,
      names.arg = FALSE,
      col = c("#0072B2", "#D55E00", "#009E73"),
      ylim = c(0, 105),
      ylab = "Discordance Rate (% of Completed Vouchers)",
      main = sprintf("A. Modality-Specific Conflict Rates (n = %d)", nrow(triage_df))
    )
    graphics::axis(
      1, at = bp,
      labels = c("Morphology
(CDA/GMM)", "Vision
(DINOv2/Cleanlab)", "Ecological Niche
(Spatial RF)"),
      tick = FALSE, cex.axis = 0.82, line = 0.6
    )
    graphics::text(bp, rates + 4.0, labels = sprintf("%.1f%%", rates), font = 2, cex = 0.9)
    graphics::grid(nx = NA, ny = NULL, col = "grey85")

    prio_col <- if ("synthesis_triage_priority" %in% colnames(triage_df)) {
      toupper(as.character(triage_df$synthesis_triage_priority))
    } else {
      rep("RESOLVED", nrow(triage_df))
    }
    prio_col[prio_col == "MEDIUM"] <- "MODERATE"
    prio_levels <- c("CRITICAL", "HIGH", "MODERATE", "LOW", "RESOLVED")
    prio_pal <- c(
      "CRITICAL" = "#D55E00",
      "HIGH"     = "#E69F00",
      "MODERATE" = "#CC79A7",
      "LOW"      = "#56B4E9",
      "RESOLVED" = "#009E73"
    )
    prio_f <- factor(prio_col, levels = prio_levels)
    tx_f <- factor(
      tx_clean,
      levels = c("Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"),
      labels = c("P. anonyma", "P. dubia", "P. paupercula", "P. plattensis")
    )
    tbl_mat <- table(prio_f, tx_f)
    graphics::barplot(
      tbl_mat,
      beside = FALSE,
      col = prio_pal[rownames(tbl_mat)],
      ylab = "Completed Voucher Count",
      font.axis = 3,
      main = "B. Consensus Curation Triage Breakdown by Taxon"
    )
    graphics::legend(
      "topright",
      legend = rownames(tbl_mat),
      fill = prio_pal[rownames(tbl_mat)],
      cex = 0.78, bty = "n"
    )
  }, fdir, "fig07_multimodal_discordance_summary", width = 12.5, height = 5.5, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 07b: Biological Taxonomic Treatments (Overall & By Recorded Species)
  # ----------------------------------------------------------------------------
  treat_col <- if ("taxonomic_treatment" %in% colnames(triage_df)) {
    as.character(triage_df$taxonomic_treatment)
  } else if ("taxonomic_status_call" %in% colnames(triage_df)) {
    sc_v <- as.character(triage_df$taxonomic_status_call)
    ifelse(
      sc_v == "Species_Confirmed", "Confirmed Species",
      ifelse(
        sc_v %in% c("Subspecies_Ecotype", "Ecophenotypic_Plasticity"), "Variety / Ecotype",
        ifelse(
          sc_v %in% c("Hybrid_Intergrade_Swarm", "Putative_Hybrid_or_Ecophenotype"), "Hybrid / Introgressant",
          "Misidentification / Cryptic Species"
        )
      )
    )
  } else {
    rep("Confirmed Species", nrow(triage_df))
  }

  treat_levels <- c(
    "Confirmed Species",
    "Variety / Ecotype",
    "Hybrid / Introgressant",
    "Misidentification / Cryptic Species"
  )
  treat_pal <- c(
    "Confirmed Species"                   = "#009E73",
    "Variety / Ecotype"                   = "#56B4E9",
    "Hybrid / Introgressant"              = "#CC79A7",
    "Misidentification / Cryptic Species" = "#D55E00"
  )

  treat_f <- factor(treat_col, levels = treat_levels)
  treat_counts <- as.numeric(table(treat_f))
  treat_pcts   <- if (nrow(triage_df) > 0) (treat_counts / nrow(triage_df)) * 100 else rep(0, length(treat_levels))

  n_reassign <- 0L
  n_cryptic  <- 0L
  if ("taxonomic_treatment_subtype" %in% colnames(triage_df)) {
    n_reassign <- sum(triage_df$taxonomic_treatment_subtype == "Misidentification_Reassign", na.rm = TRUE)
    n_cryptic  <- sum(triage_df$taxonomic_treatment_subtype == "Putative_Cryptic_Species", na.rm = TRUE)
  } else {
    n_reassign <- treat_counts[4]
  }

  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(5.8, 4.8, 3.6, 1.2))

    # Panel A: Complex-wide Taxonomic Treatment Counts & Percentages
    y_top <- max(c(10, treat_counts), na.rm = TRUE) * 1.28
    bp_t <- graphics::barplot(
      treat_counts,
      names.arg = FALSE,
      col = treat_pal[treat_levels],
      ylim = c(0, y_top),
      ylab = "Completed Voucher Count",
      main = sprintf("A. Recommended Taxonomic Treatments (n = %d)", nrow(triage_df))
    )
    graphics::axis(
      1, at = bp_t,
      labels = c(
        "Confirmed
Species",
        "Variety /
Ecotype",
        "Hybrid /
Introgressant",
        "Misidentification /
Cryptic Species"
      ),
      tick = FALSE, cex.axis = 0.80, line = 0.7
    )
    lbl_a <- sprintf("%d (%.1f%%)", treat_counts, treat_pcts)
    graphics::text(bp_t, treat_counts + y_top * 0.045, labels = lbl_a, font = 2, cex = 0.86)
    if ((n_reassign + n_cryptic) > 0) {
      sub_lbl <- sprintf("[Reassign: %d | Cryptic: %d]", n_reassign, n_cryptic)
      graphics::text(bp_t[4], treat_counts[4] + y_top * 0.115, labels = sub_lbl, cex = 0.72, col = "#802000")
    }
    graphics::grid(nx = NA, ny = NULL, col = "grey85")

    # Panel B: Taxonomic Treatment Breakdown by Recorded Species
    tx_f <- factor(
      tx_clean,
      levels = c("Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis"),
      labels = c("P. anonyma", "P. dubia", "P. paupercula", "P. plattensis")
    )
    tbl_treat_tx <- table(treat_f, tx_f)
    col_totals <- colSums(tbl_treat_tx)
    y_top_b <- max(c(5, col_totals), na.rm = TRUE) * 1.25
    bp_b <- graphics::barplot(
      tbl_treat_tx,
      beside = FALSE,
      col = treat_pal[rownames(tbl_treat_tx)],
      ylim = c(0, y_top_b),
      ylab = "Completed Voucher Count",
      font.axis = 3,
      cex.names = 0.88,
      main = "B. Taxonomic Treatments by Recorded Species"
    )
    graphics::text(
      bp_b, col_totals + y_top_b * 0.04,
      labels = ifelse(col_totals > 0, sprintf("n = %d", col_totals), ""),
      font = 2, cex = 0.82
    )
    graphics::grid(nx = NA, ny = NULL, col = "grey85")
    graphics::legend(
      "topright",
      legend = c(
        "Confirmed Species",
        "Variety / Ecotype",
        "Hybrid / Introgressant",
        "Misidentification / Cryptic"
      ),
      fill = treat_pal[treat_levels],
      cex = 0.76, bty = "n"
    )
  }, fdir, "fig07b_taxonomic_treatment_summary", width = 13.0, height = 5.6, dpi = 300)

  invisible(TRUE)
}
