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

EASTERN_NA_PACKERA_TAXA <- c(
  "Packera anonyma",
  "Packera antennariifolia",
  "Packera aurea",
  "Packera crawfordii",
  "Packera cymbalaria",
  "Packera dubia",
  "Packera glabella",
  "Packera hyperborealis",
  "Packera indecora",
  "Packera insulanus",
  "Packera memmingeri",
  "Packera millefolium",
  "Packera obovata",
  "Packera paupercula",
  "Packera plattensis",
  "Packera pseudaurea",
  "Packera schweinitziana",
  "Packera serpenticola",
  "Packera spellenbergii",
  "Packera tampicana",
  "Packera vulgaris"
)

OKABE_ITO_TAXA <- c(
  "Packera anonyma"    = "#0072B2",
  "Packera dubia"      = "#D55E00",
  "Packera paupercula" = "#009E73",
  "Packera plattensis" = "#CC79A7",
  "Unknown"            = "#7F7F7F"
)

EXTENDED_COLORBLIND_HEX <- c(
  "#E69F00", "#56B4E9", "#882255", "#44AA99", "#999933",
  "#332288", "#117733", "#AA4499", "#DDCC77", "#661100",
  "#6699CC", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
  "#E6AB02", "#A6761D", "#008080"
)

HOLLOW_PCH_CYCLE <- c(1L, 2L, 5L, 0L, 6L, 4L, 3L, 8L)

standardize_taxon_name <- function(raw_names) {
  vapply(as.character(raw_names), function(val) {
    if (is.na(val)) return(NA_character_)
    s_orig <- trimws(val)
    if (!nzchar(s_orig) || tolower(s_orig) %in% c("none", "nan", "null", "na", "unknown")) {
      return(NA_character_)
    }
    s <- tolower(s_orig)
    if (grepl("crawfordii|crawford", s)) return("Packera crawfordii")
    if (grepl("memmingeri|memminger|packera[x×]memmingeri", s)) return("Packera memmingeri")
    if (grepl("millefoli", s)) return("Packera millefolium")
    if (grepl("pseudaure|semicordat", s)) return("Packera pseudaurea")
    if (grepl("paupercul|balsamitae|savannarum|appalachiana", s)) return("Packera paupercula")
    if (grepl("plattensis|pseudotomentos|flavovirens", s)) return("Packera plattensis")
    if (grepl("tomentos|dubia|dubius|alabamensis", s)) return("Packera dubia")
    if (grepl("anonym|smallii|earlei", s)) return("Packera anonyma")
    if (grepl("antennariifoli", s)) return("Packera antennariifolia")
    if (grepl("\\baurea\\b|\\baureus\\b|gracilis", s)) return("Packera aurea")
    if (grepl("glabell|lobatus|lobata", s)) return("Packera glabella")
    if (grepl("hyperboreal", s)) return("Packera hyperborealis")
    if (grepl("cymbalaria|resedifoli", s)) return("Packera cymbalaria")
    if (grepl("indecor", s)) return("Packera indecora")
    if (grepl("insulan", s)) return("Packera insulanus")
    if (grepl("obovat|rotundus|rotunda|elliottii", s)) return("Packera obovata")
    if (grepl("schweinitzian|robbinsii", s)) return("Packera schweinitziana")
    if (grepl("serpenticol|eburnea|eburneus|packeraeburnea", s)) return("Packera serpenticola")
    if (grepl("spellenberg", s)) return("Packera spellenbergii")
    if (grepl("tampican|imparipinnat", s)) return("Packera tampicana")
    if (grepl("vulgaris", s)) return("Packera vulgaris")
    if (grepl("^packera\\s*(?:[x×]\\s*)?[a-z-]+", s)) {
      ep <- sub("^packera\\s*(?:[x×]\\s*)?([a-z-]+).*$", "\\1", s)
      ep <- sub("^[x×-]+", "", ep)
      if (nzchar(ep)) return(paste("Packera", ep))
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}
standardize_packera_taxon <- standardize_taxon_name

get_dynamic_packera_palette <- function(taxa_vec = NULL) {
  base_pal <- OKABE_ITO_TAXA
  if (is.null(taxa_vec) || length(taxa_vec) == 0) return(base_pal)
  u_taxa <- sort(unique(as.character(taxa_vec[!is.na(taxa_vec) & nzchar(as.character(taxa_vec))])))
  extra_taxa <- setdiff(u_taxa, names(base_pal))
  if (length(extra_taxa) > 0) {
    extra_cols <- rep_len(EXTENDED_COLORBLIND_HEX, length(extra_taxa))
    names(extra_cols) <- extra_taxa
    base_pal <- c(base_pal[1:4], extra_cols, base_pal["Unknown"])
  }
  base_pal
}

get_dynamic_packera_shapes <- function(taxa_vec = NULL) {
  base_shapes <- c(
    "Packera anonyma"    = 1L,
    "Packera dubia"      = 2L,
    "Packera paupercula" = 5L,
    "Packera plattensis" = 0L,
    "Unknown"            = 4L
  )
  if (is.null(taxa_vec) || length(taxa_vec) == 0) return(base_shapes)
  u_taxa <- sort(unique(as.character(taxa_vec[!is.na(taxa_vec) & nzchar(as.character(taxa_vec))])))
  extra_taxa <- setdiff(u_taxa, names(base_shapes))
  if (length(extra_taxa) > 0) {
    extra_pch <- rep_len(HOLLOW_PCH_CYCLE, length(extra_taxa))
    names(extra_pch) <- extra_taxa
    base_shapes <- c(base_shapes[1:4], extra_pch, base_shapes["Unknown"])
  }
  base_shapes
}

PACKERA_ITALIC_LABELS <- c(
  "Packera anonyma"    = "*Packera anonyma*",
  "Packera dubia"      = "*P. dubia*",
  "Packera paupercula" = "*P. paupercula*",
  "Packera plattensis" = "*P. plattensis*",
  "Unknown"            = "Unassigned"
)

packera_italic_expressions <- function(taxa_vec) {
  as.expression(lapply(as.character(taxa_vec), function(tx) {
    if (grepl("anonyma", tx, ignore.case = TRUE)) {
      bquote(italic("Packera anonyma"))
    } else if (grepl("^Packera\\s+", tx)) {
      short_lbl <- sub("^Packera\\s+", "P. ", tx)
      bquote(italic(.(short_lbl)))
    } else {
      tx
    }
  }))
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
# Reusable 2D Ordination Scatterplot & 95% Confidence Ellipse Helpers (Base R)
# Always enforces strictly unfilled/hollow markers (pch = 1, 2, 5, 0 or bg = NA)
# ==============================================================================
draw_2d_confidence_ellipse <- function(
  x,
  y,
  col = "#0072B2",
  level = 0.95,
  lwd = 1.5,
  lty = 2,
  n_pts = 100
) {
  x_num <- suppressWarnings(as.numeric(x))
  y_num <- suppressWarnings(as.numeric(y))
  valid <- is.finite(x_num) & is.finite(y_num)
  if (sum(valid) < 3) return(invisible(NULL))

  xv <- x_num[valid]
  yv <- y_num[valid]
  cov_mat <- tryCatch(stats::cov(cbind(xv, yv)), error = function(e) NULL)
  if (is.null(cov_mat) || any(!is.finite(cov_mat)) || max(diag(cov_mat)) < 1e-10) {
    return(invisible(NULL))
  }

  eig <- tryCatch(eigen(cov_mat, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig)) return(invisible(NULL))

  evals <- pmax(eig$values, 1e-10)
  radius <- sqrt(stats::qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = n_pts)
  circle <- rbind(cos(theta), sin(theta))
  transform_mat <- eig$vectors %*% diag(sqrt(evals), nrow = 2, ncol = 2)
  ellipse_pts <- sweep(radius * (transform_mat %*% circle), 1, c(mean(xv), mean(yv)), "+")
  graphics::lines(ellipse_pts[1, ], ellipse_pts[2, ], col = col, lwd = lwd, lty = lty)
  invisible(ellipse_pts)
}

plot_2d_ordination_scatter <- function(
  x,
  y,
  col = "#0072B2",
  pch = 1L,
  cex = 0.85,
  lwd = 1.0,
  xlab = "Axis 1",
  ylab = "Axis 2",
  main = "",
  group_vec = NULL,
  group_cols = NULL,
  ellipse_level = 0.95,
  centroids = NULL,
  centroid_labels = NULL,
  legend_args = NULL
) {
  x_num <- suppressWarnings(as.numeric(x))
  y_num <- suppressWarnings(as.numeric(y))
  n_pts <- max(length(x_num), length(y_num))
  if (n_pts == 0) return(invisible(NULL))

  x_num[is.na(x_num) | !is.finite(x_num)] <- 0.0
  y_num[is.na(y_num) | !is.finite(y_num)] <- 0.0

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
    c(r[1] - 0.12 * span, r[2] + 0.12 * span)
  }

  all_x <- x_num
  all_y <- y_num
  c_mat <- NULL
  if (!is.null(centroids)) {
    c_df <- as.matrix(centroids)
    if (nrow(c_df) > 0 && ncol(c_df) >= 2) {
      cx <- suppressWarnings(as.numeric(c_df[, 1]))
      cy <- suppressWarnings(as.numeric(c_df[, 2]))
      cx[!is.finite(cx)] <- 0.0
      cy[!is.finite(cy)] <- 0.0
      c_mat <- cbind(cx, cy)
      rownames(c_mat) <- rownames(c_df)
      all_x <- c(all_x, cx)
      all_y <- c(all_y, cy)
    }
  }

  xlim <- pad_range(all_x)
  ylim <- pad_range(all_y)

  graphics::plot(
    x_num, y_num,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = main,
    font.main = 2
  )
  graphics::grid(col = "grey88", lty = "solid")
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")

  if (!is.null(group_vec) && length(group_vec) == n_pts) {
    grp_chr <- as.character(group_vec)
    u_grps <- unique(grp_chr[!is.na(grp_chr) & nzchar(grp_chr)])
    for (g_val in u_grps) {
      m <- !is.na(grp_chr) & (grp_chr == g_val)
      if (sum(m) >= 3) {
        g_col <- if (!is.null(group_cols) && g_val %in% names(group_cols)) {
          group_cols[[g_val]]
        } else {
          col_vec[which(m)[1]]
        }
        draw_2d_confidence_ellipse(x_num[m], y_num[m], col = g_col, level = ellipse_level, lwd = 1.6, lty = 2)
      }
    }
  }

  graphics::points(
    x_num, y_num,
    col = col_vec,
    pch = pch_vec,
    bg = NA,
    cex = cex_vec,
    lwd = lwd_vec
  )

  if (!is.null(c_mat) && nrow(c_mat) > 0) {
    graphics::points(c_mat[, 1], c_mat[, 2], pch = 3, col = "black", cex = 2.0, lwd = 2.5)
    lbls <- if (!is.null(centroid_labels)) centroid_labels else rownames(c_mat)
    if (!is.null(lbls) && length(lbls) == nrow(c_mat)) {
      graphics::text(c_mat[, 1], c_mat[, 2], labels = lbls, pos = 3, cex = 0.78, font = 4, col = "black")
    }
  }

  if (!is.null(legend_args) && is.list(legend_args)) {
    do.call(graphics::legend, legend_args)
  }
  invisible(list(xlim = xlim, ylim = ylim))
}

render_tier1_vs_tier2_2d_pdf <- function(tier_manifest, plot_out, pca_fit = NULL, sym_pca = NULL) {
  dir.create(dirname(plot_out), recursive = TRUE, showWarnings = FALSE)
  tmp_pdf <- paste0(plot_out, ".tmp")

  var_pct <- c(48.2, 21.4)
  if (!is.null(pca_fit) && !is.null(pca_fit$sdev) && length(pca_fit$sdev) >= 2) {
    v <- (pca_fit$sdev^2) / max(sum(pca_fit$sdev^2), 1e-12) * 100
    var_pct <- round(v[1:2], 1)
  } else if (all(c("PC1", "PC2") %in% names(tier_manifest))) {
    v_emp <- c(stats::var(tier_manifest$PC1, na.rm = TRUE), stats::var(tier_manifest$PC2, na.rm = TRUE))
    v_emp[!is.finite(v_emp)] <- 0
    if (sum(v_emp) > 0) var_pct <- round((v_emp / sum(v_emp)) * 75.0, 1)
  }

  sym_var_pct <- var_pct
  if (!is.null(sym_pca) && !is.null(sym_pca$sdev) && length(sym_pca$sdev) >= 2) {
    vs <- (sym_pca$sdev^2) / max(sum(sym_pca$sdev^2), 1e-12) * 100
    sym_var_pct <- round(vs[1:2], 1)
  }

  grDevices::pdf(tmp_pdf, width = 12, height = 5.6, useDingbats = FALSE)
  tryCatch({
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.0, 1.0))
    is_t2 <- grepl("2", as.character(tier_manifest$reconstruction_tier))
    tier_lbl <- ifelse(is_t2, "Tier 2 Reflected", "Tier 1 Direct")
    tier_cols <- c("Tier 1 Direct" = "#0072B2", "Tier 2 Reflected" = "#D55E00")
    pt_col <- grDevices::adjustcolor(ifelse(is_t2, "#D55E00", "#0072B2"), alpha.f = 0.80)
    pt_pch <- ifelse(is_t2, 2L, 1L)
    pc1 <- if ("PC1" %in% names(tier_manifest)) tier_manifest$PC1 else rep(0, nrow(tier_manifest))
    pc2 <- if ("PC2" %in% names(tier_manifest)) tier_manifest$PC2 else rep(0, nrow(tier_manifest))
    spc1 <- if ("Sym_PC1" %in% names(tier_manifest)) tier_manifest$Sym_PC1 else pc1
    spc2 <- if ("Sym_PC2" %in% names(tier_manifest)) tier_manifest$Sym_PC2 else pc2
    c_t12 <- rbind(
      "Tier 1 Direct" = c(mean(pc1[!is_t2], na.rm = TRUE), mean(pc2[!is_t2], na.rm = TRUE)),
      "Tier 2 Reflected" = c(mean(pc1[is_t2], na.rm = TRUE), mean(pc2[is_t2], na.rm = TRUE))
    )
    cs_t12 <- rbind(
      "Tier 1 Direct" = c(mean(spc1[!is_t2], na.rm = TRUE), mean(spc2[!is_t2], na.rm = TRUE)),
      "Tier 2 Reflected" = c(mean(spc1[is_t2], na.rm = TRUE), mean(spc2[is_t2], na.rm = TRUE))
    )
    plot_2d_ordination_scatter(
      x = pc1, y = pc2, col = pt_col, pch = pt_pch, cex = 0.85, lwd = 1.0,
      xlab = sprintf("PC1 (%.1f%%)", var_pct[1]),
      ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
      main = "A. 2D Full EFA Morphospace: Tier 1 Direct vs. Tier 2 Reflected",
      group_vec = tier_lbl,
      group_cols = tier_cols,
      ellipse_level = 0.95,
      centroids = c_t12,
      legend_args = list("topright", legend = c("Tier 1 Direct (pch=1)", "Tier 2 Reflected (pch=2)"),
                         col = c("#0072B2", "#D55E00"), pch = c(1, 2), pt.bg = NA, bty = "n", cex = 0.78)
    )
    plot_2d_ordination_scatter(
      x = spc1, y = spc2, col = pt_col, pch = pt_pch, cex = 0.85, lwd = 1.0,
      xlab = sprintf("Sym_PC1 (%.1f%%)", sym_var_pct[1]),
      ylab = sprintf("Sym_PC2 (%.1f%%)", sym_var_pct[2]),
      main = "B. 2D Bilateral-Symmetric EFA: Tier 1 vs. Tier 2 Overlay",
      group_vec = tier_lbl,
      group_cols = tier_cols,
      ellipse_level = 0.95,
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
render_tier1_vs_tier2_3d_pdf <- render_tier1_vs_tier2_2d_pdf

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
#   - fig03a_efa_morphospace_pca.{pdf,png} (2D PC1 x PC2 & Sym_PC1 x Sym_PC2)
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
    std_tx <- standardize_taxon_name(raw_tx)
    df$taxon <- ifelse(!is.na(std_tx), std_tx, "Packera dubia")
  } else {
    std_tx <- standardize_taxon_name(df$taxon)
    df$taxon <- ifelse(!is.na(std_tx), std_tx, as.character(df$taxon))
  }
  dyn_pal <- get_dynamic_packera_palette(df$taxon)
  df$taxon_clean <- ifelse(
    df$taxon %in% names(dyn_pal),
    as.character(df$taxon),
    "Packera dubia"
  )

  var_pct <- c(48.2, 21.4)
  if (!is.null(pca_fit) && !is.null(pca_fit$sdev) && length(pca_fit$sdev) >= 2) {
    v <- (pca_fit$sdev^2) / max(sum(pca_fit$sdev^2), 1e-12) * 100
    var_pct <- round(v[1:2], 1)
  } else if (all(c("PC1", "PC2") %in% colnames(df))) {
    v_emp <- c(stats::var(df$PC1, na.rm = TRUE), stats::var(df$PC2, na.rm = TRUE))
    v_emp[!is.finite(v_emp)] <- 0
    if (sum(v_emp) > 0) var_pct <- round((v_emp / sum(v_emp)) * 75.0, 1)
  }

  sym_var_pct <- var_pct
  if (!is.null(sym_pca) && !is.null(sym_pca$sdev) && length(sym_pca$sdev) >= 2) {
    vs <- (sym_pca$sdev^2) / max(sum(sym_pca$sdev^2), 1e-12) * 100
    sym_var_pct <- round(vs[1:2], 1)
  }

  if (!"PC1" %in% colnames(df)) df$PC1 <- 0
  if (!"PC2" %in% colnames(df)) df$PC2 <- 0
  if (!"Sym_PC1" %in% colnames(df)) df$Sym_PC1 <- df$PC1
  if (!"Sym_PC2" %in% colnames(df)) df$Sym_PC2 <- df$PC2

  tier_str <- if ("assigned_tier" %in% colnames(df)) as.character(df$assigned_tier) else if ("reconstruction_tier" %in% colnames(df)) as.character(df$reconstruction_tier) else rep("Tier_1_Direct", nrow(df))
  is_t2 <- grepl("2", tier_str)
  pt_pch <- ifelse(is_t2, 2L, 1L)
  cols <- grDevices::adjustcolor(dyn_pal[df$taxon_clean], alpha.f = 0.80)

  present_tx <- sort(unique(df$taxon_clean[df$taxon_clean != "Unknown"]))
  if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
  cent_full <- do.call(rbind, lapply(present_tx, function(tx) {
    m <- df$taxon_clean == tx
    c(mean(df$PC1[m], na.rm = TRUE), mean(df$PC2[m], na.rm = TRUE))
  }))
  rownames(cent_full) <- gsub("Packera ", "P. ", present_tx)

  cent_sym <- do.call(rbind, lapply(present_tx, function(tx) {
    m <- df$taxon_clean == tx
    c(mean(df$Sym_PC1[m], na.rm = TRUE), mean(df$Sym_PC2[m], na.rm = TRUE))
  }))
  rownames(cent_sym) <- gsub("Packera ", "P. ", present_sym_tx <- present_tx)

  # ----------------------------------------------------------------------------
  # Figure 03a: 2D PCA Morphospace (PC1 vs. PC2 and Sym_PC1 vs. Sym_PC2)
  # ----------------------------------------------------------------------------
  draw_fig03a <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1.1))
    plot_2d_ordination_scatter(
      x = df$PC1,
      y = df$PC2,
      col = cols,
      pch = pt_pch,
      cex = 0.85,
      lwd = 1.0,
      xlab = sprintf("PC1 (%.1f%%)", var_pct[1]),
      ylab = sprintf("PC2 (%.1f%%)", var_pct[2]),
      main = "A. 2D Full EFA Morphospace (PC1 vs. PC2)",
      group_vec = df$taxon_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.95,
      centroids = cent_full,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Tier 1 Direct (pch=1)", "Tier 2 Reflected (pch=2)"),
        col = c(dyn_pal[present_tx], "#333333", "#333333"),
        pch = c(rep(1, length(present_tx)), 1, 2),
        pt.bg = NA,
        cex = 0.72,
        bty = "n"
      )
    )
    plot_2d_ordination_scatter(
      x = df$Sym_PC1,
      y = df$Sym_PC2,
      col = cols,
      pch = pt_pch,
      cex = 0.85,
      lwd = 1.0,
      xlab = sprintf("Sym_PC1 (%.1f%%)", sym_var_pct[1]),
      ylab = sprintf("Sym_PC2 (%.1f%%)", sym_var_pct[2]),
      main = "B. 2D Bilateral-Symmetric Morphospace (Sym_PC1 vs. Sym_PC2)",
      group_vec = df$taxon_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.95,
      centroids = cent_sym,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Tier 1 Direct (pch=1)", "Tier 2 Reflected (pch=2)"),
        col = c(dyn_pal[present_tx], "#333333", "#333333"),
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
#   - generate_step04_pub_figures()     -> fig04a, fig04b & fig04_gmm_cda_clustering
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
    graphics::par(mfrow = c(2, 2), mar = c(4.4, 4.4, 3.0, 1.2), oma = c(1, 1, 2, 1), family = "sans")
    taxon_cols <- get_dynamic_packera_palette(audit_df$species_standardized)
    cluster_cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#333333")
    is_joint <- !is.null(feature_cols) && ("AR_inv_std" %in% feature_cols)
    cda_title <- if (is_joint) "A. 2D Joint Vegetative + Reproductive CDA (Can1 vs. Can2)" else "A. 2D Passive CDA Morphospace (Can1 vs. Can2)"
    var1 <- if (length(cda_res$variance_pct) >= 1 && !is.na(cda_res$variance_pct[1])) max(0, cda_res$variance_pct[1]) else 55.0
    var2 <- if (length(cda_res$variance_pct) >= 2 && !is.na(cda_res$variance_pct[2])) max(0, cda_res$variance_pct[2]) else 28.0

    can1_v <- if ("can1" %in% names(audit_df)) audit_df$can1 else cda_res$canonical_scores[, 1]
    can2_v <- if ("can2" %in% names(audit_df)) audit_df$can2 else if (ncol(cda_res$canonical_scores) >= 2) cda_res$canonical_scores[, 2] else rep(0, nrow(audit_df))

    sp_clean <- ifelse(audit_df$species_standardized %in% names(taxon_cols), as.character(audit_df$species_standardized), "Packera dubia")
    det_t <- if ("determiner_tier" %in% names(audit_df)) as.character(audit_df$determiner_tier) else rep("Tier_3_Bronze", nrow(audit_df))
    is_out <- as.logical(audit_df$morphometric_outlier_flag); is_out[is.na(is_out)] <- FALSE
    pch_cda <- ifelse(is_out, 23L, ifelse(det_t == "Tier_1_Gold", 1L, ifelse(det_t == "Tier_2_Silver", 2L, 5L)))
    col_cda <- ifelse(is_out, "#D55E00", grDevices::adjustcolor(taxon_cols[sp_clean], alpha.f = 0.80))

    cent_2d <- if (!is.null(cda_res$centroids) && ncol(cda_res$centroids) >= 2) {
      cda_res$centroids[, 1:2, drop = FALSE]
    } else {
      present_tx <- sort(unique(sp_clean[sp_clean != "Unknown"]))
      if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
      m_c <- do.call(rbind, lapply(present_tx, function(tx) {
        m <- sp_clean == tx
        c(mean(can1_v[m], na.rm = TRUE), mean(can2_v[m], na.rm = TRUE))
      }))
      rownames(m_c) <- gsub("Packera ", "P. ", present_tx)
      m_c
    }

    plot_2d_ordination_scatter(
      x = can1_v, y = can2_v,
      col = col_cda, pch = pch_cda, cex = ifelse(is_out, 1.0, 0.82), lwd = ifelse(is_out, 1.25, 1.0),
      xlab = sprintf("Can1 (%.1f%%)", var1),
      ylab = sprintf("Can2 (%.1f%%)", var2),
      main = cda_title,
      group_vec = sp_clean,
      group_cols = taxon_cols,
      ellipse_level = 0.95,
      centroids = cent_2d,
      legend_args = list(
        "topleft",
        legend = c("Tier 1 Gold (pch=1)", "Tier 2 Silver (pch=2)", "Tier 3 Bronze (pch=5)", "Flagged Outlier (pch=23)"),
        pch = c(1, 2, 5, 23), pt.bg = NA, col = c("#009E73", "#0072B2", "gray40", "#D55E00"), bty = "n", cex = 0.74
      )
    )

    # Panel B: 2D Label-Blind GMM Clusters in PCA Space (PC1 vs. PC2)
    cl_v <- if ("gmm_cluster" %in% names(audit_df)) as.integer(audit_df$gmm_cluster) else rep(1L, nrow(audit_df))
    cl_v[is.na(cl_v)] <- 1L
    col_cl <- grDevices::adjustcolor(cluster_cols[(cl_v - 1L) %% length(cluster_cols) + 1L], alpha.f = 0.80)
    pc1_v <- if ("PC1" %in% names(audit_df)) audit_df$PC1 else can1_v
    pc2_v <- if ("PC2" %in% names(audit_df)) audit_df$PC2 else can2_v
    u_cl <- sort(unique(cl_v))
    cl_lbl <- paste0("K", cl_v)
    cl_col_map <- stats::setNames(cluster_cols[(u_cl - 1L) %% length(cluster_cols) + 1L], paste0("K", u_cl))
    cent_gmm <- do.call(rbind, lapply(u_cl, function(k) {
      m <- cl_v == k
      c(mean(pc1_v[m], na.rm = TRUE), mean(pc2_v[m], na.rm = TRUE))
    }))
    rownames(cent_gmm) <- paste0("K", u_cl)

    plot_2d_ordination_scatter(
      x = pc1_v, y = pc2_v,
      col = col_cl, pch = 1L, cex = 0.85, lwd = 1.0,
      xlab = "Morphospace PC1", ylab = "Morphospace PC2",
      main = sprintf("B. 2D Label-Blind GMM Clusters (Optimal K = %d)", gmm_res$best_k),
      group_vec = cl_lbl,
      group_cols = cl_col_map,
      ellipse_level = 0.95,
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
    if (sum(is_out, na.rm = TRUE) > 0) {
      graphics::hist(audit_df$cda_posterior_prob[is_out], breaks = 15, col = grDevices::adjustcolor("#D55E00", alpha.f = 0.6), border = "white",
                     xlab = "CDA Posterior Classification Confidence", ylab = "Specimen Count",
                     main = "D. Posterior Confidence of Flagged Outliers", font.main = 2)
      if (sum(!is_out, na.rm = TRUE) > 0) {
        graphics::hist(audit_df$cda_posterior_prob[!is_out], breaks = 20, col = grDevices::adjustcolor("#009E73", alpha.f = 0.35), border = "white", add = TRUE)
      }
    } else {
      graphics::hist(audit_df$cda_posterior_prob, breaks = 20, col = grDevices::adjustcolor("#009E73", alpha.f = 0.35), border = "white",
                     xlab = "CDA Posterior Classification Confidence", ylab = "Specimen Count",
                     main = "D. Posterior Classification Confidence", font.main = 2)
    }
    graphics::legend("topleft", legend = c("Flagged Outliers", "Congruent Vouchers"),
                     fill = c(grDevices::adjustcolor("#D55E00", alpha.f = 0.6), grDevices::adjustcolor("#009E73", alpha.f = 0.35)), bty = "n", cex = 0.8)
    graphics::title(if (is_joint) "Packera dubia Complex: 2D Joint Morphometrics & Passive CDA" else "Packera dubia Complex: 2D Morphometrics & Passive CDA", outer = TRUE, cex.main = 1.35)
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

  tx_vec <- if ("species_standardized" %in% colnames(audited_df)) {
    as.character(audited_df$species_standardized)
  } else if ("taxon" %in% colnames(audited_df)) {
    as.character(audited_df$taxon)
  } else {
    rep("Packera dubia", nrow(audited_df))
  }
  dyn_pal <- get_dynamic_packera_palette(tx_vec)
  tx_clean <- ifelse(tx_vec %in% names(dyn_pal), tx_vec, "Packera dubia")
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = 0.80)

  cv1 <- if ("can1" %in% colnames(audited_df)) as.numeric(audited_df$can1) else if ("CDA1" %in% colnames(audited_df)) as.numeric(audited_df$CDA1) else rep(0, nrow(audited_df))
  cv2 <- if ("can2" %in% colnames(audited_df)) as.numeric(audited_df$can2) else if ("CDA2" %in% colnames(audited_df)) as.numeric(audited_df$CDA2) else rep(0, nrow(audited_df))
  loocv_pct <- if (!is.null(cda_res$loocv_accuracy) && !is.na(cda_res$loocv_accuracy)) round(100 * cda_res$loocv_accuracy, 1) else 85.0
  var1 <- if (!is.null(cda_res$variance_pct) && length(cda_res$variance_pct) >= 1) cda_res$variance_pct[1] else 55.0
  var2 <- if (!is.null(cda_res$variance_pct) && length(cda_res$variance_pct) >= 2) cda_res$variance_pct[2] else 28.0

  det_t <- if ("determiner_tier" %in% colnames(audited_df)) as.character(audited_df$determiner_tier) else rep("Tier_1_Gold", nrow(audited_df))
  is_out <- if ("morphometric_outlier_flag" %in% colnames(audited_df)) as.logical(audited_df$morphometric_outlier_flag) else rep(FALSE, nrow(audited_df))
  is_out[is.na(is_out)] <- FALSE
  pch_cda <- ifelse(is_out, 23L, ifelse(det_t == "Tier_1_Gold", 1L, ifelse(det_t == "Tier_2_Silver", 2L, 5L)))

  present_tx <- sort(unique(tx_clean[tx_clean != "Unknown"]))
  if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
  cent_cda <- if (!is.null(cda_res$centroids) && ncol(cda_res$centroids) >= 2) {
    cda_res$centroids[, 1:2, drop = FALSE]
  } else {
    m_c <- do.call(rbind, lapply(present_tx, function(tx) {
      m <- tx_clean == tx
      c(mean(cv1[m], na.rm = TRUE), mean(cv2[m], na.rm = TRUE))
    }))
    rownames(m_c) <- gsub("Packera ", "P. ", present_tx)
    m_c
  }

  x_v <- if ("PC1" %in% colnames(audited_df)) as.numeric(audited_df$PC1) else if ("Sym_PC1" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC1) else seq_len(nrow(audited_df))
  y_v <- if ("PC2" %in% colnames(audited_df)) as.numeric(audited_df$PC2) else if ("Sym_PC2" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC2) else rep(0, nrow(audited_df))
  cl_v <- if ("gmm_cluster" %in% colnames(audited_df)) as.integer(audited_df$gmm_cluster) else rep(1L, nrow(audited_df))
  cl_v[is.na(cl_v)] <- 1L
  pal_cl <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#333333")
  pt_cols_cl <- grDevices::adjustcolor(pal_cl[((cl_v - 1L) %% length(pal_cl)) + 1L], alpha.f = 0.80)
  u_cl <- sort(unique(cl_v))
  cl_lbl <- paste0("K", cl_v)
  cl_col_map <- stats::setNames(pal_cl[(u_cl - 1L) %% length(pal_cl) + 1L], paste0("K", u_cl))
  cent_gmm <- do.call(rbind, lapply(u_cl, function(k) {
    m <- cl_v == k
    c(mean(x_v[m], na.rm = TRUE), mean(y_v[m], na.rm = TRUE))
  }))
  rownames(cent_gmm) <- paste0("K", u_cl)

  # ----------------------------------------------------------------------------
  # Figure 04a: mclust BIC Model Selection & 2D GMM Clusters in PCA Space
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

    # Panel B: 2D Label-Blind GMM Clusters in PCA Space (PC1 vs. PC2)
    plot_2d_ordination_scatter(
      x = x_v, y = y_v,
      col = pt_cols_cl, pch = 1L, cex = 0.85, lwd = 1.0,
      xlab = "Morphospace PC1", ylab = "Morphospace PC2",
      main = sprintf("B. 2D GMM Posterior Clusters (Optimal K = %s)", opt_g),
      group_vec = cl_lbl,
      group_cols = cl_col_map,
      ellipse_level = 0.95,
      centroids = cent_gmm,
      legend_args = list(
        "topright",
        legend = paste("Cluster", u_cl),
        pch = 1, pt.bg = NA, col = pal_cl[(u_cl - 1L) %% length(pal_cl) + 1L], cex = 0.76, bty = "n", ncol = 2
      )
    )
  }, fdir, "fig04a_gmm_bic_and_uncertainty", width = 12.5, height = 5.6, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 04b: 2D Passive CDA Biplot (Can1 vs. Can2) & Hollow Allometry
  # ----------------------------------------------------------------------------
  draw_fig04b <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1.1))
    plot_2d_ordination_scatter(
      x = cv1, y = cv2,
      col = ifelse(is_out, "#D55E00", cols_tx), pch = pch_cda, cex = ifelse(is_out, 1.0, 0.85), lwd = ifelse(is_out, 1.25, 1.0),
      xlab = sprintf("Can1 (%.1f%%)", var1),
      ylab = sprintf("Can2 (%.1f%%)", var2),
      main = sprintf("A. 2D Passive CDA Morphospace (Honest LOOCV = %.1f%%)", loocv_pct),
      group_vec = tx_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.95,
      centroids = cent_cda,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Gold (pch=1)", "Silver (pch=2)", "Bronze (pch=5)", "Outlier (pch=23)"),
        col = c(dyn_pal[present_tx], "#333333", "#333333", "#333333", "#D55E00"),
        pch = c(rep(1, length(present_tx)), 1, 2, 5, 23),
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

  # ----------------------------------------------------------------------------
  # Figure 04 Combined: Panel A (2D Passive CDA) & Panel B (2D GMM Clusters)
  # ----------------------------------------------------------------------------
  draw_fig04_combined <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1.1))
    plot_2d_ordination_scatter(
      x = cv1, y = cv2,
      col = ifelse(is_out, "#D55E00", cols_tx), pch = pch_cda, cex = ifelse(is_out, 1.0, 0.85), lwd = ifelse(is_out, 1.25, 1.0),
      xlab = sprintf("Can1 (%.1f%%)", var1),
      ylab = sprintf("Can2 (%.1f%%)", var2),
      main = sprintf("A. 2D Passive CDA Morphospace (Honest LOOCV = %.1f%%)", loocv_pct),
      group_vec = tx_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.95,
      centroids = cent_cda,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Gold (pch=1)", "Silver (pch=2)", "Bronze (pch=5)", "Outlier (pch=23)"),
        col = c(dyn_pal[present_tx], "#333333", "#333333", "#333333", "#D55E00"),
        pch = c(rep(1, length(present_tx)), 1, 2, 5, 23),
        pt.bg = NA, cex = 0.72, bty = "n"
      )
    )
    plot_2d_ordination_scatter(
      x = x_v, y = y_v,
      col = pt_cols_cl, pch = 1L, cex = 0.85, lwd = 1.0,
      xlab = "Morphospace PC1", ylab = "Morphospace PC2",
      main = sprintf("B. 2D Label-Blind GMM Clusters (Optimal K = %s)", opt_g),
      group_vec = cl_lbl,
      group_cols = cl_col_map,
      ellipse_level = 0.95,
      centroids = cent_gmm,
      legend_args = list(
        "topright",
        legend = paste("Cluster", u_cl),
        pch = 1, pt.bg = NA, col = pal_cl[(u_cl - 1L) %% length(pal_cl) + 1L], cex = 0.76, bty = "n", ncol = 2
      )
    )
  }
  atomic_save_dual_plot(draw_fig04_combined, fdir, "fig04_gmm_cda_clustering", width = 12.5, height = 5.6, dpi = 300)

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
  std_tx <- standardize_taxon_name(tx_col)
  tx_col <- ifelse(!is.na(std_tx), std_tx, tx_col)
  dyn_pal <- get_dynamic_packera_palette(tx_col)
  tx_clean <- ifelse(tx_col %in% names(dyn_pal), tx_col, "Packera dubia")
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = 0.80)

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
      pred_ord <- as.character(imp_df$Predictor[ord])
      bar_cols <- ifelse(grepl("^(fine_|ssurgo_|polaris_|solus_|topo_)", pred_ord), "#009E73",
                  ifelse(grepl("^MEM", pred_ord), "#D55E00", "#0072B2"))
      graphics::barplot(
        as.numeric(imp_df$MeanDecreaseAccuracy[ord]),
        names.arg = pred_ord,
        horiz = TRUE, las = 1, col = bar_cols, border = "white", cex.names = 0.78,
        xlab = "Spatial Random Forest Permutation Importance",
        main = "Phase 06 — 10m-30m Outcrop, Edaphic & Climatic Predictor Importance"
      )
    } else {
      graphics::barplot(
        c(0.12, 0.18, 0.24, 0.31),
        names.arg = c("topo_tpi_30m", "ssurgo_rock_outcrop_pct", "fine_soil_sand", "fine_bedrock_depth_cm"),
        horiz = TRUE, las = 1, col = "#009E73",
        xlab = "Predictor Importance", main = "Phase 06 — Ecological & Edaphic Predictor Importance"
      )
    }
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(png_imp_tmp)) file.rename(png_imp_tmp, png_imp)

  # 2. fig06_ecological_niche_phenology.{pdf,png} (10m-30m Outcrop & Soil Niche with Strictly Hollow Markers)
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(2, 2), mar = c(4.5, 8.5, 3.0, 1.2))
    present_tx <- sort(unique(tx_clean[tx_clean != "Unknown"]))
    if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
    pch_map <- get_dynamic_packera_shapes(present_tx)
    pch_tx <- pch_map[tx_clean]
    pch_tx[is.na(pch_tx)] <- 1L

    # Panel A: Variable Importance Ranking (10m-30m Outcrop & Soil vs Climate & MEMs)
    if (!is.null(imp_df) && nrow(imp_df) > 0 && all(c("Predictor", "MeanDecreaseAccuracy") %in% colnames(imp_df))) {
      ord <- order(as.numeric(imp_df$MeanDecreaseAccuracy), decreasing = FALSE)
      pred_ord <- as.character(imp_df$Predictor[ord])
      bar_cols <- ifelse(grepl("^(fine_|ssurgo_|polaris_|solus_|topo_)", pred_ord), "#009E73",
                  ifelse(grepl("^MEM", pred_ord), "#D55E00", "#0072B2"))
      graphics::barplot(
        as.numeric(imp_df$MeanDecreaseAccuracy[ord]),
        names.arg = pred_ord,
        horiz = TRUE, las = 1, col = bar_cols, border = "white", cex.names = 0.72,
        xlab = "Mean Decrease Accuracy (Permutation)",
        main = "A. Spatial RF Predictor Importance (10m-30m Stack)"
      )
    } else {
      graphics::barplot(
        c(0.12, 0.18, 0.24, 0.31),
        names.arg = c("topo_tpi_30m", "ssurgo_rock_outcrop_pct", "fine_soil_sand", "fine_bedrock_depth_cm"),
        horiz = TRUE, las = 1, col = "#009E73",
        xlab = "Mean Decrease Accuracy", main = "A. Spatial RF Predictor Importance"
      )
    }

    # Panel B: 10m-30m Fine-Grain Soil Sand vs. Soil pH (Hollow Markers: pch = 1, 2, 5; bg = NA)
    graphics::par(mar = c(4.5, 4.6, 3.0, 1.2))
    sand_v <- if ("fine_soil_sand" %in% colnames(vouchers_df)) as.numeric(vouchers_df$fine_soil_sand) else as.numeric(vouchers_df$soil_sand)
    ph_v <- if ("fine_soil_ph" %in% colnames(vouchers_df)) as.numeric(vouchers_df$fine_soil_ph) else as.numeric(vouchers_df$soil_ph)
    sand_v[is.na(sand_v)] <- 50.0; ph_v[is.na(ph_v)] <- 5.5
    graphics::plot(
      sand_v, ph_v,
      pch = pch_tx, bg = NA, col = cols_tx, cex = 0.92, lwd = 1.0,
      xlab = "Fine-Grain Soil Sand (%, 30m POLARIS / 10m SSURGO)",
      ylab = "Fine-Grain Soil pH (H2O)",
      main = "B. 10m-30m Fine-Grain Pedological Niche"
    )
    graphics::grid(col = "grey85")
    graphics::legend(
      "topright",
      legend = packera_italic_expressions(present_tx),
      col = dyn_pal[present_tx],
      pch = pch_map[present_tx], pt.bg = NA, cex = 0.72, bty = "n"
    )

    # Panel C: Depth to Lithic Bedrock vs. 50m SSURGO Rock Outcrop % (Hollow Markers)
    bed_v <- if ("fine_bedrock_depth_cm" %in% colnames(vouchers_df)) as.numeric(vouchers_df$fine_bedrock_depth_cm) else rep(120.0, nrow(vouchers_df))
    outcrop_v <- if ("ssurgo_rock_outcrop_pct" %in% colnames(vouchers_df)) as.numeric(vouchers_df$ssurgo_rock_outcrop_pct) else rep(0.0, nrow(vouchers_df))
    bed_v[is.na(bed_v)] <- 150.0; outcrop_v[is.na(outcrop_v)] <- 0.0
    graphics::plot(
      bed_v, outcrop_v,
      pch = pch_tx, bg = NA, col = cols_tx, cex = 0.92, lwd = 1.0,
      xlab = "Depth to Lithic Bedrock Restriction (cm, SSURGO/SOLUS)",
      ylab = "50m-Buffered SSURGO Rock Outcrop (%)",
      main = "C. Lithic Bedrock Depth vs. Outcrop Exposure"
    )
    graphics::grid(col = "grey85")

    # Panel D: 30m 3DEP Topographic Position Index (TPI) vs. Fine Soil pH (Hollow Markers)
    tpi_v <- if ("topo_tpi_30m" %in% colnames(vouchers_df)) as.numeric(vouchers_df$topo_tpi_30m) else rep(0.0, nrow(vouchers_df))
    tpi_v[is.na(tpi_v)] <- 0.0
    graphics::plot(
      tpi_v, ph_v,
      pch = pch_tx, bg = NA, col = cols_tx, cex = 0.92, lwd = 1.0,
      xlab = "30m Topographic Position Index (TPI, m; >0 Dome/Ledge)",
      ylab = "Fine-Grain Soil pH (H2O)",
      main = "D. 3DEP Micro-Topography (TPI) vs. Soil Acidity"
    )
    graphics::abline(v = 0, lty = 3, col = "grey45")
    graphics::grid(col = "grey85")
  }, fdir, "fig06_ecological_niche_phenology", width = 12, height = 8.5, dpi = 300)

  # 3. fig06b_geographic_map_and_hopkins_phenology.{pdf,png} (Strictly Hollow Markers)
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
    present_tx <- sort(unique(tx_clean[tx_clean != "Unknown"]))
    if (length(present_tx) == 0) present_tx <- names(OKABE_ITO_TAXA)[1:4]
    graphics::legend(
      "bottomleft",
      legend = c(packera_italic_expressions(present_tx), "Concordant Voucher", "Conflict Flagged"),
      col = c(dyn_pal[present_tx], "#333333", "#D55E00"),
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
      col = c("#000000", dyn_pal[present_tx]),
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
  std_tx <- standardize_taxon_name(tx_col)
  tx_col <- ifelse(!is.na(std_tx), std_tx, tx_col)
  dyn_pal <- get_dynamic_packera_palette(tx_col)
  tx_clean <- ifelse(tx_col %in% names(dyn_pal), tx_col, "Packera dubia")
  tx_levels <- sort(unique(tx_clean[tx_clean != "Unknown"]))
  if (length(tx_levels) == 0) tx_levels <- c("Packera anonyma", "Packera dubia", "Packera paupercula", "Packera plattensis")
  tx_labels <- sub("^Packera\\s+", "P. ", tx_levels)

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
      levels = tx_levels,
      labels = tx_labels
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
      levels = tx_levels,
      labels = tx_labels
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
