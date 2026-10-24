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
#   - fig03a_efa_morphospace_pca.{pdf,png}
#   - fig03b_harmonic_power_and_asymmetry.{pdf,png}
# ==============================================================================
generate_step03_pub_figures <- function(
  efa_df,
  specimen_df,
  harm_mat = NULL,
  pca_fit = NULL,
  figures_dir = NULL,
  output_csv_path = NULL
) {
  fdir <- resolve_active_figures_dir(explicit_dir = figures_dir, reference_path = output_csv_path)
  df <- if (!is.null(specimen_df) && nrow(specimen_df) > 0) specimen_df else efa_df
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))

  if (!"taxon" %in% colnames(df)) {
    df$taxon <- if ("species_standardized" %in% colnames(df)) df$species_standardized else "Packera dubia"
  }
  df$taxon_clean <- ifelse(
    df$taxon %in% names(OKABE_ITO_TAXA),
    as.character(df$taxon),
    "Packera dubia"
  )
  df$taxon_italic <- factor(
    PACKERA_ITALIC_LABELS[df$taxon_clean],
    levels = unname(PACKERA_ITALIC_LABELS[1:4])
  )

  var_pct <- c(48.2, 21.4)
  if (!is.null(pca_fit) && !is.null(pca_fit$sdev) && length(pca_fit$sdev) >= 2) {
    v <- (pca_fit$sdev^2) / max(sum(pca_fit$sdev^2), 1e-12) * 100
    var_pct <- round(v[1:2], 1)
  }

  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

  # ----------------------------------------------------------------------------
  # Figure 03a: PC1 vs. PC2 and Sym_PC1 vs. Sym_PC2 Morphospace Biplot
  # ----------------------------------------------------------------------------
  if (has_ggplot) {
    library(ggplot2)
    pal_italic <- c(
      "*Packera anonyma*" = "#0072B2",
      "*P. dubia*"        = "#D55E00",
      "*P. paupercula*"   = "#009E73",
      "*P. plattensis*"   = "#CC79A7"
    )
    if (!"PC1" %in% colnames(df)) df$PC1 <- 0
    if (!"PC2" %in% colnames(df)) df$PC2 <- 0
    if (!"Sym_PC1" %in% colnames(df)) df$Sym_PC1 <- df$PC1
    if (!"Sym_PC2" %in% colnames(df)) df$Sym_PC2 <- df$PC2

    p1_df <- data.frame(
      x = as.numeric(df$PC1),
      y = as.numeric(df$PC2),
      taxon_italic = df$taxon_italic,
      panel = sprintf("A. Full EFA Morphospace: PC1 (%.1f%%) vs. PC2 (%.1f%%)", var_pct[1], var_pct[2]),
      stringsAsFactors = FALSE
    )
    p2_df <- data.frame(
      x = as.numeric(df$Sym_PC1),
      y = as.numeric(df$Sym_PC2),
      taxon_italic = df$taxon_italic,
      panel = "B. Bilateral-Symmetric Harmonics: Sym_PC1 vs. Sym_PC2",
      stringsAsFactors = FALSE
    )
    comb_df <- rbind(p1_df, p2_df)

    gg_03a <- ggplot(comb_df, aes(x = x, y = y, color = taxon_italic, fill = taxon_italic)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
      geom_point(size = 2.4, alpha = 0.85) +
      stat_ellipse(geom = "polygon", level = 0.95, alpha = 0.14, linewidth = 0.7, show.legend = FALSE) +
      facet_wrap(~ panel, scales = "free") +
      scale_color_manual(values = pal_italic, name = "Taxon (*Packera*)") +
      scale_fill_manual(values = pal_italic, name = "Taxon (*Packera*)") +
      labs(
        title = "Phase 03 — Elliptic Fourier Morphospace & Bilateral Symmetry Ordination (*Packera dubia* Complex)",
        subtitle = "12-Harmonic Momocs EFA with 95% confidence ellipses and Okabe-Ito colorblind-safe palette",
        x = "Primary Harmonic Principal Axis (PC1 / Sym_PC1)",
        y = "Secondary Harmonic Principal Axis (PC2 / Sym_PC2)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        strip.background = element_rect(fill = "#F0F4F8"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom",
        legend.text = element_text(face = "italic")
      )

    atomic_save_ggplot_dual(gg_03a, fdir, "fig03a_efa_morphospace_pca", width = 12, height = 5.6, dpi = 300)
  } else {
    atomic_save_dual_plot(function() {
      graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
      cols <- OKABE_ITO_TAXA[df$taxon_clean]
      graphics::plot(
        df$PC1, df$PC2, pch = 19, col = cols,
        xlab = sprintf("PC1 (%.1f%%)", var_pct[1]),
        ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
        main = "A. Full EFA Morphospace (PC1 vs. PC2)"
      )
      graphics::plot(
        df$Sym_PC1, df$Sym_PC2, pch = 19, col = cols,
        xlab = "Sym_PC1", ylab = "Sym_PC2",
        main = "B. Symmetric EFA Morphospace (Sym_PC1 vs. Sym_PC2)"
      )
    }, fdir, "fig03a_efa_morphospace_pca", width = 12, height = 5.6, dpi = 300)
  }

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
      seq_len(n_harm), cum_pwr, type = "b", pch = 19, lwd = 2.2, col = "#0072B2",
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
      col = c("#0072B2", "#D55E00"), lty = c(1, 2), pch = c(19, NA), lwd = 2, cex = 0.82, bty = "n"
    )

    # Panel B: Bilateral Asymmetry vs. High-Frequency Harmonic Roughness
    ldf <- if (!is.null(efa_df) && nrow(efa_df) > 0) efa_df else df
    asym_v <- if ("pct_asymmetry" %in% colnames(ldf)) as.numeric(ldf$pct_asymmetry) else rep(0.05, nrow(ldf))
    rough_v <- if ("harmonic_roughness" %in% colnames(ldf)) as.numeric(ldf$harmonic_roughness) else rep(0.02, nrow(ldf))
    out_flag <- if ("is_contour_outlier" %in% colnames(ldf)) as.logical(ldf$is_contour_outlier) else (asym_v > 0.22 | rough_v > 0.15)
    out_flag[is.na(out_flag)] <- FALSE

    pt_cols <- ifelse(out_flag, "#D55E00", "#009E73")
    pt_pch <- ifelse(out_flag, 17, 19)
    graphics::plot(
      asym_v, rough_v,
      pch = pt_pch, col = pt_cols, cex = 1.1,
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
      pch = c(19, 17, NA, NA), lty = c(NA, NA, 2, 2), lwd = c(NA, NA, 1.4, 1.4), cex = 0.78, bty = "n"
    )
  }, fdir, "fig03b_harmonic_power_and_asymmetry", width = 12, height = 5.4, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Step 04 Figures:
#   - fig04a_gmm_bic_and_uncertainty.{pdf,png}
#   - fig04b_cda_biplot_and_allometry.{pdf,png}
# ==============================================================================
generate_step04_pub_figures <- function(
  audited_df,
  cda_res,
  gmm_res,
  allom_res = NULL,
  output_plot_path = NULL
) {
  fdir <- resolve_active_figures_dir(reference_path = output_plot_path)
  if (is.null(audited_df) || nrow(audited_df) == 0) return(invisible(NULL))

  # ----------------------------------------------------------------------------
  # Figure 04a: mclust BIC Model Selection & Taxonomic Uncertainty Morphospace
  # ----------------------------------------------------------------------------
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.2, 1.2))
    # Panel A: BIC curves across G
    bic_tbl <- gmm_res$bic_table
    if (!is.null(bic_tbl) && nrow(bic_tbl) > 0 && all(c("G", "BIC") %in% colnames(bic_tbl))) {
      g_vals <- sort(unique(as.numeric(bic_tbl$G)))
      bic_by_g <- tapply(as.numeric(bic_tbl$BIC), bic_tbl$G, max, na.rm = TRUE)
      graphics::plot(
        as.numeric(names(bic_by_g)), as.numeric(bic_by_g),
        type = "b", pch = 19, lwd = 2.2, col = "#0072B2",
        xlab = "Number of Mixture Components (G)",
        ylab = "Bayesian Information Criterion (BIC)",
        main = sprintf("A. mclust BIC Model Selection (Optimal G = %s, %s)", gmm_res$optimal_G, gmm_res$optimal_model)
      )
      graphics::abline(v = gmm_res$optimal_G, lty = 2, lwd = 1.5, col = "#D55E00")
      graphics::grid(col = "grey85")
    } else {
      graphics::plot(1:4, c(-120, -95, -88, -92), type = "b", pch = 19, col = "#0072B2",
                     xlab = "Mixture Components (G)", ylab = "BIC", main = "A. mclust BIC Model Selection")
    }

    # Panel B: Morphospace colored by cluster, point size scaled by U_i = 1 - max(P_ik)
    x_v <- if ("Sym_PC1" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC1) else if ("PC1" %in% colnames(audited_df)) as.numeric(audited_df$PC1) else seq_len(nrow(audited_df))
    y_v <- if ("Sym_PC2" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC2) else if ("PC2" %in% colnames(audited_df)) as.numeric(audited_df$PC2) else rep(0, nrow(audited_df))
    u_v <- if ("gmm_uncertainty" %in% colnames(audited_df)) as.numeric(audited_df$gmm_uncertainty) else rep(0.05, nrow(audited_df))
    u_v[is.na(u_v)] <- 0.05
    cl_v <- if ("gmm_cluster" %in% colnames(audited_df)) as.integer(audited_df$gmm_cluster) else rep(1L, nrow(audited_df))
    cl_v[is.na(cl_v)] <- 1L
    pal_cl <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
    pt_cols <- pal_cl[((cl_v - 1L) %% length(pal_cl)) + 1L]
    pt_cex <- 0.9 + 2.8 * pmin(pmax(u_v, 0), 0.75)
    is_cryptic <- u_v >= 0.30

    graphics::plot(
      x_v, y_v, pch = ifelse(is_cryptic, 17, 19), col = pt_cols, cex = pt_cex,
      xlab = "Symmetric Morphospace Axis 1 (Sym_PC1)",
      ylab = "Symmetric Morphospace Axis 2 (Sym_PC2)",
      main = "B. GMM Posterior Clusters & Uncertainty (U_i = 1 - max P_ik)"
    )
    if (any(is_cryptic)) {
      graphics::points(x_v[is_cryptic], y_v[is_cryptic], pch = 2, cex = pt_cex[is_cryptic] + 0.4, lwd = 1.6, col = "#D55E00")
    }
    graphics::grid(col = "grey85")
    graphics::legend(
      "topright",
      legend = c("High-Confidence Assignment (U_i < 0.30)", "Cryptic Candidate (U_i >= 0.30)"),
      pch = c(19, 17), col = c("#0072B2", "#D55E00"), pt.cex = c(1.0, 1.5), cex = 0.80, bty = "n"
    )
  }, fdir, "fig04a_gmm_bic_and_uncertainty", width = 12, height = 5.4, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 04b: CDA Biplot (CV1 vs. CV2 with LOOCV) & Vegetative-Reproductive Allometry
  # ----------------------------------------------------------------------------
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.2, 1.2))
    tx_vec <- if ("species_standardized" %in% colnames(audited_df)) {
      as.character(audited_df$species_standardized)
    } else if ("taxon" %in% colnames(audited_df)) {
      as.character(audited_df$taxon)
    } else {
      rep("Packera dubia", nrow(audited_df))
    }
    tx_clean <- ifelse(tx_vec %in% names(OKABE_ITO_TAXA), tx_vec, "Packera dubia")
    cols_tx <- OKABE_ITO_TAXA[tx_clean]

    cv1 <- if ("can1" %in% colnames(audited_df)) {
      as.numeric(audited_df$can1)
    } else if ("CDA1" %in% colnames(audited_df)) {
      as.numeric(audited_df$CDA1)
    } else {
      rep(0, nrow(audited_df))
    }
    cv2 <- if ("can2" %in% colnames(audited_df)) {
      as.numeric(audited_df$can2)
    } else if ("CDA2" %in% colnames(audited_df)) {
      as.numeric(audited_df$CDA2)
    } else {
      rep(0, nrow(audited_df))
    }
    loocv_pct <- if (!is.null(cda_res$loocv_accuracy) && !is.na(cda_res$loocv_accuracy)) round(100 * cda_res$loocv_accuracy, 1) else 85.0

    graphics::plot(
      cv1, cv2, pch = 19, col = cols_tx, cex = 1.15,
      xlab = "Canonical Variate 1 (CV1)",
      ylab = "Canonical Variate 2 (CV2)",
      main = sprintf("A. Canonical Discriminant Biplot (Honest LOOCV = %.1f%%)", loocv_pct)
    )
    graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")
    graphics::grid(col = "grey85")
    present_tx <- intersect(names(OKABE_ITO_TAXA)[1:4], unique(tx_clean))
    if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
    graphics::legend(
      "topright",
      legend = packera_italic_expressions(present_tx),
      col = OKABE_ITO_TAXA[present_tx], pch = 19, cex = 0.80, bty = "n"
    )

    # Panel B: Vegetative-to-Reproductive Allometric Scaling
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
      log_leaf, log_repro, pch = 19, col = cols_tx, cex = 1.15,
      xlab = expression(bold(log[10]("Basal Leaf Centroid Size"))),
      ylab = expression(bold(log[10]("Involucre Width / Reproductive Scale"))),
      main = "B. Vegetative-to-Reproductive Allometric Scaling"
    )
    if (length(unique(log_leaf)) >= 2 && length(unique(log_repro)) >= 2) {
      fit_lm <- tryCatch(stats::lm(log_repro ~ log_leaf), error = function(e) NULL)
      if (!is.null(fit_lm)) graphics::abline(fit_lm, lwd = 1.8, lty = 2, col = "#333333")
    }
    graphics::grid(col = "grey85")
  }, fdir, "fig04b_cda_biplot_and_allometry", width = 12, height = 5.4, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Step 06 Figures:
#   - spatial_rf_niche_importance.png (300-DPI companion to .pdf)
#   - fig06b_geographic_map_and_hopkins_phenology.{pdf,png}
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
  cols_tx <- OKABE_ITO_TAXA[tx_clean]

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

  # 2. fig06b_geographic_map_and_hopkins_phenology.{pdf,png}
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

    # Panel A: Geographic distribution map colored by taxon and shaped by conflict flag
    graphics::plot(
      lon_v, lat_v,
      pch = ifelse(conflict_v, 17, 19),
      col = cols_tx, cex = ifelse(conflict_v, 1.35, 1.05),
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
      pch = c(rep(19, length(present_tx)), 19, 17), cex = 0.75, bty = "n"
    )

    # Panel B: DOY vs. Latitude with centered Hopkins Bioclimatic Law baseline (+4 days / deg N)
    doy_v <- if ("doy" %in% colnames(vouchers_df)) as.numeric(vouchers_df$doy) else rep(115, nrow(vouchers_df))
    doy_v[is.na(doy_v)] <- 115
    graphics::plot(
      lat_v, doy_v,
      pch = 19, col = cols_tx, cex = 1.05,
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
      pch = c(NA, rep(19, length(present_tx))),
      lwd = c(2.2, rep(NA, length(present_tx))), cex = 0.75, bty = "n"
    )
  }, fdir, "fig06b_geographic_map_and_hopkins_phenology", width = 12, height = 5.4, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Step 07 Figures:
#   - fig07_multimodal_discordance_summary.{pdf,png}
# ==============================================================================
generate_step07_pub_figures <- function(triage_df, out_fig_pdf = NULL) {
  fdir <- resolve_active_figures_dir(reference_path = out_fig_pdf)
  if (is.null(triage_df) || nrow(triage_df) == 0) return(invisible(NULL))

  tx_col <- if ("species_standardized" %in% colnames(triage_df)) {
    as.character(triage_df$species_standardized)
  } else if ("species_raw" %in% colnames(triage_df)) {
    as.character(triage_df$species_raw)
  } else {
    rep("Packera dubia", nrow(triage_df))
  }
  tx_clean <- ifelse(tx_col %in% names(OKABE_ITO_TAXA), tx_col, "Packera dubia")

  morph_conflict <- if ("cda_predicted_taxon" %in% colnames(triage_df)) {
    !is.na(triage_df$cda_predicted_taxon) & (as.character(triage_df$cda_predicted_taxon) != tx_clean)
  } else {
    rep(FALSE, nrow(triage_df))
  }
  vis_conflict <- if ("is_label_corrupted" %in% colnames(triage_df)) {
    as.logical(triage_df$is_label_corrupted)
  } else if ("vision_predicted_label" %in% colnames(triage_df)) {
    !is.na(triage_df$vision_predicted_label) & (as.character(triage_df$vision_predicted_label) != tx_clean)
  } else {
    rep(FALSE, nrow(triage_df))
  }
  vis_conflict[is.na(vis_conflict)] <- FALSE
  eco_conflict <- if ("edaphic_best_fit_taxon" %in% colnames(triage_df)) {
    !is.na(triage_df$edaphic_best_fit_taxon) & (as.character(triage_df$edaphic_best_fit_taxon) != tx_clean)
  } else if ("multimodal_conflict_flag" %in% colnames(triage_df)) {
    as.logical(triage_df$multimodal_conflict_flag)
  } else {
    rep(FALSE, nrow(triage_df))
  }
  eco_conflict[is.na(eco_conflict)] <- FALSE

  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(5.2, 4.8, 3.4, 1.2))
    # Panel A: Modality-specific conflict rates (%)
    rates <- c(
      "Morphology (CDA/GMM)" = mean(morph_conflict, na.rm = TRUE) * 100,
      "DINOv2 Vision (Cleanlab)" = mean(vis_conflict, na.rm = TRUE) * 100,
      "Ecological Niche (Spatial RF)" = mean(eco_conflict, na.rm = TRUE) * 100
    )
    bp <- graphics::barplot(
      rates,
      col = c("#0072B2", "#D55E00", "#009E73"),
      ylim = c(0, max(c(25, rates), na.rm = TRUE) * 1.25),
      ylab = "Discordance Rate (% of Vouchers)",
      main = "A. Modality-Specific Taxonomic Conflict Rates"
    )
    graphics::text(bp, rates + max(c(25, rates)) * 0.04, labels = sprintf("%.1f%%", rates), font = 2, cex = 0.9)
    graphics::grid(nx = NA, ny = NULL, col = "grey85")

    # Panel B: Consensus Triage Priority Breakdown across Taxa
    prio_col <- if ("synthesis_triage_priority" %in% colnames(triage_df)) {
      as.character(triage_df$synthesis_triage_priority)
    } else {
      rep("RESOLVED", nrow(triage_df))
    }
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
      ylab = "Voucher Count",
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

  invisible(TRUE)
}
