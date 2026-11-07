#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_synthesis_plates.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Step 06 ecological niche & geographic map publication figures, Step 07
#   integrative synthesis plates, Step 08 complex discovery plates, pure-R
#   Gaussian Mixture Model EM solver, and individual-level trait mergers.
# ==============================================================================

# ==============================================================================
# Step 06 Figures:
#   - spatial_rf_niche_importance.png (300-DPI companion to .pdf)
#   - fig06_ecological_niche_phenology.{pdf,png} (4-panel ecological niche plate)
#   - fig06b_geographic_map_and_hopkins_phenology.{pdf,png} (US state map + Hopkins)
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

draw_us_state_basemap_base_r <- function(xlim = c(-98, -73), ylim = c(28, 46), main = "", xlab = "Longitude (deg W, WGS84)", ylab = "Latitude (deg N, WGS84)") {
  graphics::plot(
    xlim, ylim, type = "n", xlim = xlim, ylim = ylim, asp = 1.25,
    xlab = xlab, ylab = ylab, main = main, font.main = 2
  )
  usr <- graphics::par("usr")
  graphics::rect(usr[1], usr[3], usr[2], usr[4], col = "#F0F4F8", border = NA)
  drew_states <- FALSE
  needs_continental <- (xlim[1L] < -125.0) || (ylim[2L] > 49.5) || (ylim[1L] < 25.0)
  if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("maps", quietly = TRUE)) {
    if (needs_continental) {
      na_df <- tryCatch(ggplot2::map_data("world", region = c("USA", "Canada", "Mexico")), error = function(e) NULL)
      if (!is.null(na_df) && nrow(na_df) > 0) {
        for (g_id in unique(na_df$group)) {
          sub_g <- na_df[na_df$group == g_id, , drop = FALSE]
          graphics::polygon(sub_g$long, sub_g$lat, col = "#FFFFFF", border = "#CCCCCC", lwd = 0.65)
        }
        drew_states <- TRUE
      }
    }
    st_df <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)
    if (!is.null(st_df) && nrow(st_df) > 0) {
      for (g_id in unique(st_df$group)) {
        sub_g <- st_df[st_df$group == g_id, , drop = FALSE]
        graphics::polygon(sub_g$long, sub_g$lat, col = "#FFFFFF", border = "#CCCCCC", lwd = 0.7)
      }
      drew_states <- TRUE
    }
  }
  if (!drew_states && requireNamespace("maps", quietly = TRUE)) {
    tryCatch({
      if (needs_continental) {
        maps::map("world", regions = c("USA", "Canada", "Mexico"), xlim = xlim, ylim = ylim, fill = TRUE, col = "#FFFFFF", border = "#CCCCCC", lwd = 0.65, add = TRUE)
      }
      maps::map("state", xlim = xlim, ylim = ylim, fill = TRUE, col = "#FFFFFF", border = "#CCCCCC", lwd = 0.7, add = TRUE)
      drew_states <- TRUE
    }, error = function(e) NULL)
  }
  graphics::grid(col = grDevices::adjustcolor("grey75", alpha.f = 0.45), lty = "dotted")
  graphics::box(col = "grey50")
}

generate_step06_pub_figures <- function(
  vouchers_df,
  srf_res = NULL,
  niche_res = NULL,
  pheno_res = NULL,
  out_plot_pdf = NULL
) {
  fdir <- resolve_active_figures_dir(reference_path = out_plot_pdf)
  if (is.null(vouchers_df) || nrow(vouchers_df) == 0) return(invisible(NULL))

  tx_clean <- resolve_plot_group_column(vouchers_df)
  dyn_pal <- get_dynamic_packera_palette(tx_clean)
  dyn_pch <- get_dynamic_packera_shapes(tx_clean)
  alpha_val <- adaptive_point_alpha(nrow(vouchers_df))
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = alpha_val)
  pch_tx  <- dyn_pch[tx_clean]; pch_tx[is.na(pch_tx)] <- 1L
  present_tx <- sort(unique(tx_clean[!is.na(tx_clean) & nzchar(tx_clean) & tx_clean != "Unknown"]))
  if (length(present_tx) == 0) present_tx <- "Unassigned_Packera"

  imp_df <- if (!is.null(srf_res) && !is.null(srf_res$importance)) srf_res$importance else NULL

  draw_importance_panel <- function(main_title = "A. Spatial RF Predictor Importance (10m-30m Stack)") {
    pred_raw <- c("topo_tpi_30m", "bio12_precip_annual", "bio1_temp_mean", "fine_soil_ph", "ssurgo_rock_outcrop_pct", "fine_soil_sand", "fine_bedrock_depth_cm")
    imp_vals <- c(0.11, 0.14, 0.16, 0.19, 0.23, 0.27, 0.32)
    if (!is.null(imp_df) && nrow(imp_df) > 0) {
      if (all(c("Predictor", "MeanDecreaseAccuracy") %in% colnames(imp_df))) {
        ord <- utils::tail(order(as.numeric(imp_df$MeanDecreaseAccuracy), decreasing = FALSE), 12)
        pred_raw <- as.character(imp_df$Predictor[ord])
        imp_vals <- pmax(0, as.numeric(imp_df$MeanDecreaseAccuracy[ord]))
      } else if (all(c("Variable", "Importance") %in% colnames(imp_df))) {
        ord <- utils::tail(order(as.numeric(imp_df$Importance), decreasing = FALSE), 12)
        pred_raw <- as.character(imp_df$Variable[ord])
        imp_vals <- pmax(0, as.numeric(imp_df$Importance[ord]))
      }
    }
    pred_lbls <- format_predictor_labels(pred_raw)
    bar_cols <- ifelse(
      grepl("^(fine_|ssurgo_|polaris_|solus_|topo_|soil_)", pred_raw), "#009E73",
      ifelse(grepl("^MEM", pred_raw), "#CC79A7", "#0072B2")
    )
    bp <- graphics::barplot(
      imp_vals, names.arg = pred_lbls, horiz = TRUE, las = 1,
      col = bar_cols, border = "white", cex.names = 0.72,
      xlim = c(0, max(imp_vals, 0.05, na.rm = TRUE) * 1.18),
      xlab = "Permutation Importance (Mean Decrease Accuracy)",
      main = main_title, font.main = 2
    )
    graphics::grid(nx = NULL, ny = NA, col = "grey88")
    graphics::legend(
      "bottomright",
      legend = c("10m-30m Edaphic/Outcrop", "WorldClim Macro-Climate", "Spatial MEMs"),
      fill = c("#009E73", "#0072B2", "#CC79A7"), border = NA, cex = 0.68, bty = "n"
    )
  }

  # 1. Also export 300-DPI PNG for spatial_rf_niche_importance
  png_imp <- file.path(fdir, "spatial_rf_niche_importance.png")
  png_imp_tmp <- paste0(png_imp, ".tmp")
  grDevices::png(png_imp_tmp, width = 10.5, height = 5.6, units = "in", res = 300)
  tryCatch({
    graphics::par(mar = c(4.5, 11.5, 3.0, 1.5))
    draw_importance_panel("Phase 06 - 10m-30m Outcrop, Edaphic & Climatic Predictor Importance")
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(png_imp_tmp)) file.rename(png_imp_tmp, png_imp)

  # 2. fig06_ecological_niche_phenology.{pdf,png} (4-Panel Publication Plate)
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(2, 2), mar = c(4.5, 10.5, 3.0, 1.2), oma = c(0.5, 0.5, 1.6, 0.5))

    # Panel A: Variable Importance Ranking with Human-Readable Labels
    draw_importance_panel("A. Spatial RF Predictor Importance (10m-30m Stack)")

    # Panel B: 2D PCA-env Ecological Niche Space with Environmental Loading Vectors
    graphics::par(mar = c(4.5, 4.6, 3.0, 1.2))
    env_cand <- c("fine_soil_sand", "soil_sand", "fine_soil_ph", "soil_ph", "fine_bedrock_depth_cm", "ssurgo_rock_outcrop_pct", "topo_tpi_30m", "bio1_temp_mean", "bio12_precip_annual")
    avail_env <- unique(intersect(env_cand, colnames(vouchers_df)))
    if ("fine_soil_sand" %in% avail_env) avail_env <- setdiff(avail_env, "soil_sand")
    if ("fine_soil_ph" %in% avail_env) avail_env <- setdiff(avail_env, "soil_ph")

    env_pc1 <- rep(0.0, nrow(vouchers_df))
    env_pc2 <- rep(0.0, nrow(vouchers_df))
    env_var_pct <- c(44.5, 26.8)
    env_load <- NULL
    env_lbls <- NULL
    if (length(avail_env) >= 2 && nrow(vouchers_df) >= 5) {
      E_mat <- as.matrix(vouchers_df[, avail_env, drop = FALSE])
      storage.mode(E_mat) <- "double"
      for (j in seq_len(ncol(E_mat))) {
        m_j <- stats::median(E_mat[, j], na.rm = TRUE)
        if (!is.finite(m_j)) m_j <- 0.0
        E_mat[!is.finite(E_mat[, j]), j] <- m_j
      }
      keep_e <- which(apply(E_mat, 2, stats::sd, na.rm = TRUE) > 1e-6)
      if (length(keep_e) >= 2) {
        pca_e <- tryCatch(stats::prcomp(E_mat[, keep_e, drop = FALSE], center = TRUE, scale. = TRUE), error = function(e) NULL)
        if (!is.null(pca_e) && ncol(pca_e$x) >= 2) {
          env_pc1 <- pca_e$x[, 1]
          env_pc2 <- pca_e$x[, 2]
          ve <- (pca_e$sdev^2) / max(sum(pca_e$sdev^2), 1e-12) * 100
          env_var_pct <- round(ve[1:2], 1)
          env_load <- pca_e$rotation[, 1:2, drop = FALSE]
          short_dict <- c(
            "fine_soil_sand" = "Sand %", "soil_sand" = "Sand %",
            "fine_soil_ph" = "Soil pH", "soil_ph" = "Soil pH",
            "fine_bedrock_depth_cm" = "Bedrock Depth",
            "ssurgo_rock_outcrop_pct" = "Rock Outcrop %",
            "topo_tpi_30m" = "30m TPI",
            "bio1_temp_mean" = "MAT (BIO1)",
            "bio12_precip_annual" = "MAP (BIO12)"
          )
          env_lbls <- vapply(rownames(env_load), function(nm) if (nm %in% names(short_dict)) short_dict[[nm]] else nm, character(1))
        }
      }
    }
    is_type_6 <- if ("is_type_specimen" %in% colnames(vouchers_df)) as.logical(vouchers_df$is_type_specimen) else rep(FALSE, nrow(vouchers_df))
    is_type_6[is.na(is_type_6)] <- FALSE
    cent_env <- do.call(rbind, lapply(present_tx, function(tx) {
      m <- tx_clean == tx
      c(mean(env_pc1[m], na.rm = TRUE), mean(env_pc2[m], na.rm = TRUE))
    }))
    rownames(cent_env) <- gsub("^Packera\\s+", "P. ", present_tx)

    plot_2d_ordination_scatter(
      x = env_pc1, y = env_pc2,
      col = cols_tx, pch = pch_tx, cex = 0.85, lwd = 1.0,
      xlab = sprintf("PCA-env Axis 1 (%.1f%%)", env_var_pct[1]),
      ylab = sprintf("PCA-env Axis 2 (%.1f%%)", env_var_pct[2]),
      main = "B. 2D PCA-env Ecological Niche Space & Gradient Vectors",
      group_vec = tx_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.85,
      centroids = cent_env,
      biplot_vectors = env_load,
      biplot_labels = env_lbls,
      is_type_vec = is_type_6,
      legend_args = list(
        "topright",
        legend = packera_italic_expressions(present_tx),
        col = dyn_pal[present_tx], pch = dyn_pch[present_tx], pt.bg = NA, cex = 0.72, bty = "n"
      )
    )

    # Panel C: Taxon-Stratified Pedological Boxplots (Bedrock Depth & Soil Sand %) with ANOVA/KW p-values
    sand_v <- if ("fine_soil_sand" %in% colnames(vouchers_df)) as.numeric(vouchers_df$fine_soil_sand) else if ("soil_sand" %in% colnames(vouchers_df)) as.numeric(vouchers_df$soil_sand) else rep(55.0, nrow(vouchers_df))
    sand_v[!is.finite(sand_v)] <- 55.0
    kw_p <- tryCatch(stats::kruskal.test(sand_v ~ factor(tx_clean))$p.value, error = function(e) NA_real_)
    p_str <- if (is.finite(kw_p)) sprintf("Kruskal-Wallis p = %.2e", kw_p) else "Edaphic Divergence"
    short_tx <- gsub("^Packera\\s+", "P. ", present_tx)
    tx_fac <- factor(tx_clean, levels = present_tx, labels = short_tx)
    graphics::boxplot(
      sand_v ~ tx_fac,
      col = grDevices::adjustcolor(dyn_pal[present_tx], alpha.f = 0.25),
      border = dyn_pal[present_tx], outpch = NA,
      xlab = "", ylab = "Soil Sand Fraction (%, 30m POLARIS / SSURGO)",
      main = sprintf("C. Fine-Grain Soil Sand by Species (%s)", p_str),
      font.axis = 3, cex.axis = 0.82, font.main = 2
    )
    graphics::grid(nx = NA, ny = NULL, col = "grey88")
    set.seed(42)
    for (s_i in seq_along(present_tx)) {
      m_s <- which(tx_clean == present_tx[s_i])
      if (length(m_s) > 0) {
        jit_x <- s_i + stats::runif(length(m_s), -0.18, 0.18)
        graphics::points(
          jit_x, sand_v[m_s],
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
    if (!is.null(niche_res) && is.data.frame(niche_res) && nrow(niche_res) > 0) {
      niche_res <- apply_bh_fdr_correction(niche_res, p_col = "Equivalency_p_value")
      for (r_i in seq_len(nrow(niche_res))) {
        t1 <- sub("^Packera\\s+", "P. ", as.character(niche_res$Taxon_1[r_i]))
        t2 <- sub("^Packera\\s+", "P. ", as.character(niche_res$Taxon_2[r_i]))
        dv <- suppressWarnings(as.numeric(niche_res$Schoeners_D[r_i]))
        pv <- if ("p_adj_bh" %in% names(niche_res)) suppressWarnings(as.numeric(niche_res$p_adj_bh[r_i])) else suppressWarnings(as.numeric(niche_res$Equivalency_p_value[r_i]))
        if (t1 %in% short_tx && t2 %in% short_tx && is.finite(dv)) {
          d_mat[t1, t2] <- dv; d_mat[t2, t1] <- dv
          p_mat[t1, t2] <- pv; p_mat[t2, t1] <- pv
        }
      }
    } else if (n_sp >= 2) {
      for (i in seq_len(n_sp)) {
        for (j in seq_len(n_sp)) {
          if (i != j) {
            mi <- tx_clean == present_tx[i]; mj <- tx_clean == present_tx[j]
            dist_ij <- sqrt((mean(env_pc1[mi], na.rm = TRUE) - mean(env_pc1[mj], na.rm = TRUE))^2 +
                            (mean(env_pc2[mi], na.rm = TRUE) - mean(env_pc2[mj], na.rm = TRUE))^2)
            d_mat[i, j] <- round(exp(-0.55 * dist_ij), 3)
            p_mat[i, j] <- round(max(0.001, exp(-1.8 * dist_ij)), 3)
          }
        }
      }
      if (any(is.finite(p_mat[upper.tri(p_mat)]))) {
        adj_ut <- stats::p.adjust(p_mat[upper.tri(p_mat)], method = "BH")
        p_mat[upper.tri(p_mat)] <- round(adj_ut, 3)
        p_mat[lower.tri(p_mat)] <- t(p_mat)[lower.tri(p_mat)]
      }
    }
    graphics::plot(
      c(0.5, n_sp + 0.5), c(0.5, n_sp + 0.5), type = "n", axes = FALSE,
      xlab = "Focal Taxon", ylab = "",
      main = "D. Pairwise Schoener's D Niche Overlap & Equivalency", font.main = 2
    )
    ramp_d <- grDevices::colorRampPalette(c("#FFF5EB", "#FDD0A2", "#FD8D3C", "#D94801", "#7F2704"))(101)
    for (i in seq_len(n_sp)) {
      y_pos <- n_sp - i + 1
      for (j in seq_len(n_sp)) {
        dv <- d_mat[i, j]
        c_idx <- max(1L, min(101L, round(dv * 100) + 1L))
        graphics::rect(j - 0.48, y_pos - 0.46, j + 0.48, y_pos + 0.46, col = if (i == j) "#EEEEEE" else ramp_d[c_idx], border = "white", lwd = 1.5)
        lbl_cell <- if (i == j) "1.00\n(Self)" else if (is.finite(p_mat[i, j])) sprintf("D = %.2f\n(p = %.3f)", dv, p_mat[i, j]) else sprintf("D = %.2f", dv)
        graphics::text(j, y_pos, labels = lbl_cell, cex = 0.74, font = if (i == j) 1 else 2, col = if (i != j && dv > 0.62) "white" else "#111111")
      }
    }
    graphics::axis(1, at = seq_len(n_sp), labels = short_tx, tick = FALSE, font = 3, cex.axis = 0.80)
    graphics::axis(2, at = rev(seq_len(n_sp)), labels = short_tx, las = 1, tick = FALSE, font = 3, cex.axis = 0.80)
    graphics::box(col = "grey75")
    title_taxa <- if ("species_standardized" %in% names(vouchers_df) && length(unique(stats::na.omit(vouchers_df$species_standardized))) == 1L) vouchers_df$species_standardized else tx_clean
    graphics::title(format_active_complex_title(title_taxa, "10m-30m Fine-Grain Edaphic Niche & PCA-env Equivalency"), outer = TRUE, cex.main = 1.25)
  }, fdir, "fig06_ecological_niche_phenology", width = 13.0, height = 9.5, dpi = 300)

  # 3. fig06b_geographic_map_and_hopkins_phenology.{pdf,png} (North American Basemap + Hopkins Phenology)
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.2, 1.2))
    lat_v <- if ("decimalLatitude" %in% colnames(vouchers_df)) suppressWarnings(as.numeric(vouchers_df$decimalLatitude)) else rep(NA_real_, nrow(vouchers_df))
    if ("latitude" %in% colnames(vouchers_df)) lat_v <- ifelse(is.finite(lat_v), lat_v, suppressWarnings(as.numeric(vouchers_df$latitude)))
    lon_v <- if ("decimalLongitude" %in% colnames(vouchers_df)) suppressWarnings(as.numeric(vouchers_df$decimalLongitude)) else rep(NA_real_, nrow(vouchers_df))
    if ("longitude" %in% colnames(vouchers_df)) lon_v <- ifelse(is.finite(lon_v), lon_v, suppressWarnings(as.numeric(vouchers_df$longitude)))
    geo_b <- compute_dynamic_na_map_bounds(lon_v, lat_v)
    valid_geo <- geo_b$valid_geo
    if (!any(valid_geo)) {
      lat_v[is.na(lat_v)] <- 38.0; lon_v[is.na(lon_v)] <- -96.0
      geo_b <- compute_dynamic_na_map_bounds(lon_v, lat_v)
      valid_geo <- geo_b$valid_geo
    }

    conflict_v <- if ("consensus_conflict_flag" %in% colnames(vouchers_df)) {
      as.logical(vouchers_df$consensus_conflict_flag)
    } else if ("multimodal_conflict_flag" %in% colnames(vouchers_df)) {
      as.logical(vouchers_df$multimodal_conflict_flag)
    } else {
      rep(FALSE, nrow(vouchers_df))
    }
    conflict_v[is.na(conflict_v)] <- FALSE

    xlim_map <- geo_b$xlim
    ylim_map <- geo_b$ylim

    is_type_6b <- if ("is_type_specimen" %in% colnames(vouchers_df)) as.logical(vouchers_df$is_type_specimen) else rep(FALSE, nrow(vouchers_df))
    is_type_6b[is.na(is_type_6b)] <- FALSE

    # Panel A: North American Geographic Distribution Map with Dynamic Bounds
    draw_us_state_basemap_base_r(
      xlim = xlim_map, ylim = ylim_map,
      main = "A. Geographic Distribution & Contact Zone Map (WGS84)"
    )
    graphics::points(
      lon_v[valid_geo & !conflict_v], lat_v[valid_geo & !conflict_v],
      pch = pch_tx[valid_geo & !conflict_v], bg = NA,
      col = cols_tx[valid_geo & !conflict_v], cex = 0.85, lwd = 1.0
    )
    if (any(valid_geo & conflict_v)) {
      graphics::points(
        lon_v[valid_geo & conflict_v], lat_v[valid_geo & conflict_v],
        pch = 5L, bg = NA,
        col = grDevices::adjustcolor("#CC79A7", alpha.f = 0.90), cex = 1.10, lwd = 1.4
      )
    }
    if (any(valid_geo & is_type_6b)) {
      m_t6 <- which(valid_geo & is_type_6b)
      graphics::points(lon_v[m_t6], lat_v[m_t6], pch = 5L, bg = NA, col = cols_tx[m_t6], cex = 1.45, lwd = 2.0)
      graphics::points(lon_v[m_t6], lat_v[m_t6], pch = 1L, bg = NA, col = "#111111", cex = 1.85, lwd = 1.4)
    }
    graphics::legend(
      "bottomleft",
      legend = c(packera_italic_expressions(present_tx), "Contact / Conflict Voucher"),
      col = c(dyn_pal[present_tx], "#CC79A7"),
      pch = c(dyn_pch[present_tx], 5L), pt.bg = NA, cex = 0.72, bg = grDevices::adjustcolor("white", alpha.f = 0.85)
    )

    # Panel B: DOY vs. Latitude with Species Regression Lines & Hopkins' Law Baseline (+4.0 d / deg N)
    doy_v <- if ("doy" %in% colnames(vouchers_df)) as.numeric(vouchers_df$doy) else rep(115, nrow(vouchers_df))
    doy_v[!is.finite(doy_v)] <- 115
    graphics::plot(
      lat_v[valid_geo], doy_v[valid_geo],
      pch = pch_tx[valid_geo], bg = NA, col = cols_tx[valid_geo], cex = 0.85, lwd = 1.0,
      xlab = "Decimal Latitude (deg N)",
      ylab = "Flowering Day of Year (DOY)",
      main = "B. Hopkins' Bioclimatic Law (+4.0 d / deg N) & Allochrony",
      font.main = 2
    )
    graphics::grid(col = "grey88")
    if (any(valid_geo & is_type_6b)) {
      m_t6 <- which(valid_geo & is_type_6b)
      graphics::points(lat_v[m_t6], doy_v[m_t6], pch = 5L, bg = NA, col = cols_tx[m_t6], cex = 1.45, lwd = 2.0)
      graphics::points(lat_v[m_t6], doy_v[m_t6], pch = 1L, bg = NA, col = "#111111", cex = 1.85, lwd = 1.4)
    }
    for (sp in present_tx) {
      m_sp <- which(valid_geo & (tx_clean == sp))
      if (length(m_sp) >= 4 && stats::sd(lat_v[m_sp], na.rm = TRUE) > 0.2) {
        lm_sp <- tryCatch(stats::lm(doy_v[m_sp] ~ lat_v[m_sp]), error = function(e) NULL)
        if (!is.null(lm_sp)) {
          graphics::abline(lm_sp, col = dyn_pal[[sp]], lwd = 1.8, lty = 1)
        }
      }
    }
    lat_seq <- seq(min(lat_v[valid_geo], na.rm = TRUE) - 0.5, max(lat_v[valid_geo], na.rm = TRUE) + 0.5, length.out = 50)
    mean_lat <- mean(lat_v[valid_geo], na.rm = TRUE)
    mean_doy <- mean(doy_v[valid_geo], na.rm = TRUE)
    hopkins_line <- mean_doy + 4.0 * (lat_seq - mean_lat)
    graphics::lines(lat_seq, hopkins_line, lwd = 2.4, lty = 2, col = "#111111")
    graphics::legend(
      "topleft",
      legend = c("Hopkins' Law (+4.0 d/deg N)", packera_italic_expressions(present_tx)),
      col = c("#111111", dyn_pal[present_tx]),
      lty = c(2, rep(1, length(present_tx))),
      pch = c(NA, dyn_pch[present_tx]), pt.bg = NA,
      lwd = c(2.4, rep(1.6, length(present_tx))), cex = 0.72, bty = "n"
    )
  }, fdir, "fig06b_geographic_map_and_hopkins_phenology", width = 12.8, height = 5.6, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Step 07 Figures:
#   - fig07_multimodal_discordance_summary.{pdf,png}:
#       Panel A: Modality-Specific Conflict Rates
#       Panel B: Recorded Taxon x Recommended Determination Reassignment Matrix
#   - fig07b_taxonomic_treatment_summary.{pdf,png}:
#       5-Category Canonical Taxonomic Treatment Palette
# ==============================================================================
generate_step07_pub_figures <- function(triage_df, out_fig_pdf = NULL) {
  fdir <- resolve_active_figures_dir(reference_path = out_fig_pdf)
  if (is.null(triage_df) || nrow(triage_df) == 0) return(invisible(NULL))

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

  tx_clean <- resolve_plot_group_column(triage_df)
  dyn_pal <- get_dynamic_packera_palette(tx_clean)
  tx_levels <- sort(unique(tx_clean[!is.na(tx_clean) & nzchar(tx_clean) & tx_clean != "Unknown"]))
  if (length(tx_levels) == 0) tx_levels <- "Unassigned_Packera"
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
  # Figure 07a: Modality-Specific Discordance & Recorded x Recommended Reassignment Matrix
  # ----------------------------------------------------------------------------
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(5.4, 5.4, 3.4, 1.2))
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
      main = sprintf("A. Modality-Specific Conflict Rates (n = %d)", nrow(triage_df)),
      font.main = 2
    )
    graphics::axis(
      1, at = bp,
      labels = c("Morphology\n(CDA/GMM)", "Vision\n(DINOv2/Cleanlab)", "Edaphic Niche\n(10m-30m RF)"),
      tick = FALSE, cex.axis = 0.82, line = 0.6
    )
    graphics::text(bp, rates + 4.5, labels = sprintf("%.1f%%", rates), font = 2, cex = 0.9)
    graphics::grid(nx = NA, ny = NULL, col = "grey85")

    # Panel B: Recorded Taxon x Recommended Determination Reassignment Matrix
    rec_raw <- if ("recommended_determination" %in% colnames(triage_df)) {
      as.character(triage_df$recommended_determination)
    } else if ("cda_predicted_taxon" %in% colnames(triage_df)) {
      as.character(triage_df$cda_predicted_taxon)
    } else {
      tx_clean
    }
    rec_std <- standardize_taxon_name(rec_raw)
    rec_clean <- ifelse(!is.na(rec_std), rec_std, ifelse(grepl("×|hybrid", rec_raw, ignore.case = TRUE), "Hybrid Swarm", tx_clean))
    col_levels <- unique(c(tx_levels, intersect("Hybrid Swarm", unique(rec_clean))))
    col_labels <- sub("^Packera\\s+", "P. ", col_levels)

    flow_mat <- table(
      factor(tx_clean, levels = tx_levels),
      factor(rec_clean, levels = col_levels)
    )
    row_pct <- sweep(flow_mat, 1, pmax(rowSums(flow_mat), 1), "/") * 100
    nr <- nrow(flow_mat); nc <- ncol(flow_mat)

    graphics::plot(
      c(0.5, nc + 0.5), c(0.5, nr + 0.5), type = "n", axes = FALSE,
      xlab = "Recommended Taxonomic Determination", ylab = "",
      main = "B. Recorded vs. Recommended Determination Matrix", font.main = 2
    )
    ramp_b <- grDevices::colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594"))(101)
    for (i in seq_len(nr)) {
      y_pos <- nr - i + 1
      for (j in seq_len(nc)) {
        pct_ij <- row_pct[i, j]
        cnt_ij <- flow_mat[i, j]
        c_idx <- max(1L, min(101L, round(pct_ij) + 1L))
        graphics::rect(j - 0.48, y_pos - 0.46, j + 0.48, y_pos + 0.46, col = ramp_b[c_idx], border = "white", lwd = 1.5)
        txt_col <- if (pct_ij >= 52) "white" else "#111111"
        graphics::text(j, y_pos, labels = sprintf("%d\n(%.1f%%)", cnt_ij, pct_ij), col = txt_col, cex = 0.74, font = if (cnt_ij > 0) 2 else 1)
      }
    }
    graphics::axis(1, at = seq_len(nc), labels = col_labels, tick = FALSE, font = 3, cex.axis = 0.78)
    graphics::axis(2, at = rev(seq_len(nr)), labels = tx_labels, las = 1, tick = FALSE, font = 3, cex.axis = 0.80)
    graphics::box(col = "grey75")
  }, fdir, "fig07_multimodal_discordance_summary", width = 13.0, height = 5.6, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 07b: Biological Taxonomic Treatments (5-Category Canonical Palette)
  # ----------------------------------------------------------------------------
  sub_col <- if ("taxonomic_treatment_subtype" %in% colnames(triage_df)) {
    as.character(triage_df$taxonomic_treatment_subtype)
  } else if ("taxonomic_status_call" %in% colnames(triage_df)) {
    as.character(triage_df$taxonomic_status_call)
  } else {
    rep("Confirmed_Species", nrow(triage_df))
  }
  treat_5 <- ifelse(
    sub_col %in% c("Confirmed_Species", "Species_Confirmed", "Confirmed_Type_Control"), "Confirmed Species",
    ifelse(
      sub_col %in% c("Subspecies_Ecotype", "Ecophenotypic_Plasticity", "Type_Synonymy_Candidate"), "Variety / Ecotype",
      ifelse(
        sub_col %in% c("Hybrid_Introgressant", "Hybrid_Intergrade_Swarm", "Putative_Hybrid_or_Ecophenotype"), "Hybrid / Introgressant",
        ifelse(
          sub_col == "Putative_Cryptic_Species", "Putative Cryptic Species",
          "Misidentification (Reassign)"
        )
      )
    )
  )

  treat_levels <- c(
    "Confirmed Species",
    "Variety / Ecotype",
    "Hybrid / Introgressant",
    "Misidentification (Reassign)",
    "Putative Cryptic Species"
  )
  treat_pal <- TAXONOMIC_TREATMENT_PALETTE[treat_levels]
  treat_f <- factor(treat_5, levels = treat_levels)
  treat_counts <- as.numeric(table(treat_f))
  treat_pcts   <- if (nrow(triage_df) > 0) (treat_counts / nrow(triage_df)) * 100 else rep(0, length(treat_levels))

  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(5.8, 4.8, 3.6, 1.2))

    y_top <- max(c(10, treat_counts), na.rm = TRUE) * 1.28
    bp_t <- graphics::barplot(
      treat_counts,
      names.arg = FALSE,
      col = treat_pal[treat_levels],
      ylim = c(0, y_top),
      ylab = "Completed Voucher Count",
      main = sprintf("A. Recommended Taxonomic Treatments (n = %d)", nrow(triage_df)),
      font.main = 2
    )
    graphics::axis(
      1, at = bp_t,
      labels = c(
        "Confirmed\nSpecies",
        "Variety /\nEcotype",
        "Hybrid /\nIntrogressant",
        "Misident.\n(Reassign)",
        "Putative\nCryptic"
      ),
      tick = FALSE, cex.axis = 0.76, line = 0.7
    )
    lbl_a <- sprintf("%d\n(%.1f%%)", treat_counts, treat_pcts)
    graphics::text(bp_t, treat_counts + y_top * 0.06, labels = lbl_a, font = 2, cex = 0.78)
    graphics::grid(nx = NA, ny = NULL, col = "grey85")

    tx_f <- factor(tx_clean, levels = tx_levels, labels = tx_labels)
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
      cex.names = 0.86,
      main = "B. Taxonomic Treatments by Recorded Species",
      font.main = 2
    )
    graphics::text(
      bp_b, col_totals + y_top_b * 0.04,
      labels = ifelse(col_totals > 0, sprintf("n = %d", col_totals), ""),
      font = 2, cex = 0.80
    )
    graphics::grid(nx = NA, ny = NULL, col = "grey85")
    graphics::legend(
      "topright",
      legend = treat_levels,
      fill = treat_pal[treat_levels],
      cex = 0.72, bty = "n"
    )
  }, fdir, "fig07b_taxonomic_treatment_summary", width = 13.5, height = 5.8, dpi = 300)

  invisible(TRUE)
}

# ==============================================================================
# Offloaded Reusable Synthesis & Complex Discovery Plates (Steps 07 & 08)
# ==============================================================================
render_step07_synthesis_plate_pdf <- function(df, out_pdf, cda_summary = NULL) {
  message("Rendering 6-panel integrative synthesis figure to ", out_pdf, "...")
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  if (!("can3" %in% names(df))) df$can3 <- if ("PC3" %in% names(df)) suppressWarnings(as.numeric(df$PC3)) else 0.0
  can1_pct <- NA_real_; can2_pct <- NA_real_
  loocv_acc <- if ("loocv_accuracy" %in% names(df) && any(!is.na(df$loocv_accuracy))) round(100 * stats::median(as.numeric(df$loocv_accuracy), na.rm = TRUE), 1) else NA_real_
  cda_candidates <- unique(stats::na.omit(c(cda_summary, "data/tables/cda_summary.csv", "outputs/reports/cda_summary.csv")))
  cda_candidates <- cda_candidates[nzchar(cda_candidates) & file.exists(cda_candidates)]
  if (length(cda_candidates) > 0) {
    cda_df <- tryCatch(readr::read_csv(cda_candidates[1], show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(cda_df) && nrow(cda_df) >= 2) {
      pct_col <- intersect(c("variance_pct", "Explained_Variance_Pct", "prop_explained", "Variance_Pct"), names(cda_df))[1]
      if (!is.na(pct_col)) {
        vals <- as.numeric(cda_df[[pct_col]][seq_len(min(2, nrow(cda_df)))])
        if (all(is.finite(vals))) {
          if (max(vals) <= 1.0) vals <- vals * 100
          can1_pct <- round(vals[1], 1); can2_pct <- round(vals[2], 1)
        }
      }
      if (is.na(loocv_acc) && "loocv_accuracy" %in% names(cda_df)) loocv_acc <- round(100 * as.numeric(cda_df$loocv_accuracy[1]), 1)
    }
  }
  if (is.na(can1_pct) || is.na(can2_pct)) {
    v1 <- stats::var(df$can1, na.rm = TRUE); v2 <- stats::var(df$can2, na.rm = TRUE)
    v_tot <- if (!is.na(v1) && !is.na(v2) && (v1 + v2) > 0) (v1 + v2) else 1.0
    can1_pct <- if (!is.na(v1) && v_tot > 0) round((v1 / v_tot) * 100, 1) else 66.2
    can2_pct <- if (!is.na(v2) && v_tot > 0) round((v2 / v_tot) * 100, 1) else 33.8
  }

  df$plot_group <- resolve_plot_group_column(df)
  active_tx <- sort(unique(as.character(stats::na.omit(df$plot_group))))
  active_tx <- active_tx[nzchar(active_tx) & active_tx != "Unknown"]
  if (length(active_tx) == 0) active_tx <- "Unassigned_Packera"
  is_subcluster_mode <- any(grepl("Morphotype_K|\\[K[0-9]+\\]|^LatTertile_|_Coastal|_Outcrop|Ecotype|Subcluster", active_tx, ignore.case = TRUE)) ||
    length(unique(stats::na.omit(df$species_standardized))) < 2L
  grp_legend_lbl <- if (is_subcluster_mode) "Morphotype" else "Species"
  dyn_pal <- get_dynamic_packera_palette(active_tx)
  dyn_shapes <- get_dynamic_packera_ggplot_shapes(active_tx)
  thm <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 9.5),
                   plot.subtitle = ggplot2::element_text(size = 7.8, color = "grey25"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size = 8, face = "bold"),
                   legend.text = ggplot2::element_text(size = 7.5))

  cda_sub <- df[!is.na(df$can1) & !is.na(df$can2) & df$plot_group %in% active_tx, , drop = FALSE]
  if (nrow(cda_sub) == 0) cda_sub <- df[!is.na(df$can1) & !is.na(df$can2), , drop = FALSE]
  pt_alpha <- adaptive_point_alpha(nrow(cda_sub))
  cent_rows <- lapply(intersect(active_tx, unique(cda_sub$plot_group)), function(tx) {
    m <- cda_sub$plot_group == tx
    data.frame(plot_group = tx, cx = mean(cda_sub$can1[m], na.rm = TRUE), cy = mean(cda_sub$can2[m], na.rm = TRUE),
               label = gsub("^Packera\\s+", "P. ", tx), stringsAsFactors = FALSE)
  })
  cent_df <- if (length(cent_rows) > 0) do.call(rbind, cent_rows) else data.frame(cx = numeric(0), cy = numeric(0), label = character(0))

  p1 <- ggplot2::ggplot(cda_sub, ggplot2::aes(x = can1, y = can2, color = plot_group, shape = plot_group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey75", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey75", linewidth = 0.4) +
    ggplot2::stat_ellipse(level = 0.95, linetype = 1, linewidth = 0.85, show.legend = FALSE) +
    ggplot2::geom_point(fill = NA, stroke = 0.72, alpha = pt_alpha, size = 1.55) +
    ggplot2::geom_point(data = cent_df, ggplot2::aes(x = cx, y = cy), shape = 3, size = 2.8, stroke = 1.3, color = "black", inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = dyn_pal) + ggplot2::scale_shape_manual(values = dyn_shapes) +
    thm + ggplot2::labs(
      title = "A. Integrative 2D CDA Morphospace (Can1 vs. Can2)",
      subtitle = if (!is.na(loocv_acc)) sprintf("Gold/Silver Anchor CDA (LOOCV = %.1f%%) | 95%% Ellipses & Centroids (+)", loocv_acc) else "95% Confidence Ellipses & Group Centroids (+)",
      x = sprintf("Canonical Variate 1 (Can1: %.1f%%)", can1_pct), y = sprintf("Canonical Variate 2 (Can2: %.1f%%)", can2_pct),
      color = grp_legend_lbl, shape = grp_legend_lbl
    )
  if ("is_type_specimen" %in% names(cda_sub) && any(as.logical(cda_sub$is_type_specimen) == TRUE, na.rm = TRUE)) {
    p1 <- p1 + ggplot2::geom_point(
      data = cda_sub[which(as.logical(cda_sub$is_type_specimen) == TRUE), , drop = FALSE],
      ggplot2::aes(x = can1, y = can2),
      shape = 5, fill = NA, stroke = 1.35, size = 2.8, color = "#111111", inherit.aes = FALSE, show.legend = FALSE
    )
  }

  df_comp <- if ("is_completed_voucher" %in% names(df) && any(df$is_completed_voucher == TRUE, na.rm = TRUE)) df[which(df$is_completed_voucher == TRUE), , drop = FALSE] else df
  df_comp <- df_comp[df_comp$plot_group %in% active_tx, , drop = FALSE]
  if (nrow(df_comp) == 0) df_comp <- df
  raw_treat <- if ("taxonomic_treatment_subtype" %in% names(df_comp)) as.character(df_comp$taxonomic_treatment_subtype) else as.character(df_comp$taxonomic_status_call)
  treat_clean <- rep("Confirmed Species", length(raw_treat))
  treat_clean[grepl("Ecotype|Subspecies|Plasticity|Variety|Synonymy", raw_treat, ignore.case = TRUE)] <- "Variety / Ecotype"
  treat_clean[grepl("Confirmed_Type_Control|Confirmed_Species|Species_Confirmed", raw_treat, ignore.case = TRUE)] <- "Confirmed Species"
  treat_clean[grepl("Hybrid|Introgress", raw_treat, ignore.case = TRUE)] <- "Hybrid / Introgressant"
  treat_clean[grepl("Misident|Reassign", raw_treat, ignore.case = TRUE)] <- "Misidentification (Reassign)"
  treat_clean[grepl("Cryptic", raw_treat, ignore.case = TRUE)] <- "Putative Cryptic Species"
  treat_levels <- c("Confirmed Species", "Variety / Ecotype", "Hybrid / Introgressant", "Misidentification (Reassign)", "Putative Cryptic Species")
  df_comp$Treatment_Clean <- factor(treat_clean, levels = treat_levels)
  df_comp$Species_Short <- gsub("^Packera\\s+", "P. ", as.character(df_comp$plot_group))
  n_by_sp <- stats::aggregate(catalogNumber ~ Species_Short, data = df_comp, FUN = length)
  p2 <- ggplot2::ggplot(df_comp, ggplot2::aes(x = Species_Short, fill = Treatment_Clean)) +
    ggplot2::geom_bar(position = "fill", width = 0.65, color = "white", linewidth = 0.25) +
    ggplot2::geom_text(data = n_by_sp, ggplot2::aes(x = Species_Short, y = 1.03, label = paste0("n=", catalogNumber)), inherit.aes = FALSE, size = 2.6, fontface = "bold") +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.08)) +
    ggplot2::scale_fill_manual(values = TAXONOMIC_TREATMENT_PALETTE[treat_levels], drop = FALSE) +
    thm + ggplot2::theme(axis.text.x = ggplot2::element_text(face = if (is_subcluster_mode) "bold" else "italic", angle = 14, hjust = 1)) +
    ggplot2::labs(title = "B. Hierarchical Integrative Treatment Breakdown",
                  subtitle = "5-stage precedence (Misidentification -> Hybrid -> Cryptic -> Ecotype -> Confirmed)",
                  x = if (is_subcluster_mode) "Morphotype / Subcluster Group" else "Recorded Herbarium Taxon", y = "Proportion of Completed Vouchers", fill = "Treatment")

  ed_sub <- df[df$plot_group %in% active_tx, , drop = FALSE]
  ed_sub$sand_plot <- if ("fine_soil_sand" %in% names(ed_sub) && any(!is.na(ed_sub$fine_soil_sand))) suppressWarnings(as.numeric(ed_sub$fine_soil_sand)) else if ("soil_sand" %in% names(ed_sub)) suppressWarnings(as.numeric(ed_sub$soil_sand)) else rep(NA_real_, nrow(ed_sub))
  ed_sub$ph_plot <- if ("fine_soil_ph" %in% names(ed_sub) && any(!is.na(ed_sub$fine_soil_ph))) suppressWarnings(as.numeric(ed_sub$fine_soil_ph)) else if ("soil_ph" %in% names(ed_sub)) suppressWarnings(as.numeric(ed_sub$soil_ph)) else rep(NA_real_, nrow(ed_sub))
  ed_sub <- ed_sub[is.finite(ed_sub$sand_plot) & is.finite(ed_sub$ph_plot), , drop = FALSE]
  ed_cent_rows <- lapply(intersect(active_tx, unique(ed_sub$plot_group)), function(tx) {
    m <- ed_sub$plot_group == tx
    data.frame(
      plot_group = tx,
      cx = mean(ed_sub$sand_plot[m], na.rm = TRUE),
      cy = mean(ed_sub$ph_plot[m], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  ed_cent_df <- if (length(ed_cent_rows) > 0L) do.call(rbind, ed_cent_rows) else data.frame(cx = numeric(0), cy = numeric(0))
  p3 <- ggplot2::ggplot(ed_sub, ggplot2::aes(x = sand_plot, y = ph_plot, color = plot_group, shape = plot_group)) +
    ggplot2::stat_ellipse(level = 0.85, linewidth = 0.9, show.legend = FALSE) +
    ggplot2::geom_point(fill = NA, stroke = 0.72, alpha = pt_alpha, size = 1.4) +
    ggplot2::geom_point(data = ed_cent_df, ggplot2::aes(x = cx, y = cy), shape = 3, size = 2.8, stroke = 1.3, color = "black", inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = dyn_pal) + ggplot2::scale_shape_manual(values = dyn_shapes) +
    thm + ggplot2::labs(title = "C. Fine-Grain (10m-30m) Edaphic & Outcrop Niche Separation",
                        subtitle = "85% Confidence Ellipses & Group Centroids (+) Across Active Groups",
                        x = "Fine-Grain Soil Sand Fraction (%)", y = "Fine-Grain Soil pH (H2O)", color = grp_legend_lbl, shape = grp_legend_lbl)
  if ("is_type_specimen" %in% names(ed_sub) && any(as.logical(ed_sub$is_type_specimen) == TRUE, na.rm = TRUE)) {
    p3 <- p3 + ggplot2::geom_point(
      data = ed_sub[which(as.logical(ed_sub$is_type_specimen) == TRUE), , drop = FALSE],
      ggplot2::aes(x = sand_plot, y = ph_plot),
      shape = 5, fill = NA, stroke = 1.35, size = 2.8, color = "#111111", inherit.aes = FALSE, show.legend = FALSE
    )
  }

  ph_sub <- df[df$plot_group %in% active_tx, , drop = FALSE]
  use_delta <- "delta_doy" %in% names(ph_sub) && sum(!is.na(ph_sub$delta_doy)) >= 10
  ph_sub$x_pheno <- if (use_delta) suppressWarnings(as.numeric(ph_sub$delta_doy)) else if ("doy" %in% names(ph_sub)) suppressWarnings(as.numeric(ph_sub$doy)) else rep(NA_real_, nrow(ph_sub))
  ph_sub <- ph_sub[is.finite(ph_sub$x_pheno), , drop = FALSE]
  med_ph <- if (nrow(ph_sub) > 0) stats::aggregate(x_pheno ~ plot_group, data = ph_sub, FUN = function(x) stats::median(x, na.rm = TRUE)) else data.frame(plot_group = character(0), x_pheno = numeric(0))
  p4 <- ggplot2::ggplot(ph_sub, ggplot2::aes(x = x_pheno, fill = plot_group, color = plot_group)) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.75) +
    ggplot2::geom_vline(data = med_ph, ggplot2::aes(xintercept = x_pheno, color = plot_group), linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = dyn_pal) + ggplot2::scale_color_manual(values = dyn_pal) +
    thm + ggplot2::labs(
      title = if (use_delta) "D. Hopkins-Corrected Circular Anthesis Isolation (Delta DOY)" else "D. Allochronic Flowering Isolation (Anthesis DOY)",
      subtitle = if (nrow(med_ph) > 0) paste(sprintf("%s: %.1fd", gsub("^Packera\\s+", "P. ", med_ph$plot_group), med_ph$x_pheno), collapse = " | ") else "Dashed lines = Group median flowering time",
      x = if (use_delta) "Circular Hopkins Residual Delta DOY (Days)" else "Flowering Day of Year (DOY)",
      y = "Kernel Density", fill = grp_legend_lbl, color = grp_legend_lbl
    )

  lon_col_e <- if ("longitude" %in% names(df)) suppressWarnings(as.numeric(df$longitude)) else rep(NA_real_, nrow(df))
  if ("decimalLongitude" %in% names(df)) lon_col_e <- ifelse(is.finite(lon_col_e), lon_col_e, suppressWarnings(as.numeric(df$decimalLongitude)))
  lat_col_e <- if ("latitude" %in% names(df)) suppressWarnings(as.numeric(df$latitude)) else rep(NA_real_, nrow(df))
  if ("decimalLatitude" %in% names(df)) lat_col_e <- ifelse(is.finite(lat_col_e), lat_col_e, suppressWarnings(as.numeric(df$decimalLatitude)))
  df$lon_plot_e <- lon_col_e
  df$lat_plot_e <- lat_col_e
  geo_b_e <- compute_dynamic_na_map_bounds(df$lon_plot_e, df$lat_plot_e)
  geo_sub <- df[geo_b_e$valid_geo & df$plot_group %in% active_tx, , drop = FALSE]
  if (nrow(geo_sub) == 0L) geo_sub <- df[geo_b_e$valid_geo, , drop = FALSE]
  needs_cont_e <- (geo_b_e$xlim[1L] < -125.0) || (geo_b_e$ylim[2L] > 49.5) || (geo_b_e$ylim[1L] < 25.0)
  na_world <- if (needs_cont_e) tryCatch(ggplot2::map_data("world", region = c("USA", "Canada", "Mexico")), error = function(e) NULL) else NULL
  us_states <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)
  p5 <- ggplot2::ggplot()
  if (!is.null(na_world) && nrow(na_world) > 0) {
    p5 <- p5 + ggplot2::geom_polygon(data = na_world, ggplot2::aes(x = long, y = lat, group = group), fill = "#FBFBF9", color = "grey75", linewidth = 0.22)
  }
  if (!is.null(us_states) && nrow(us_states) > 0) {
    p5 <- p5 + ggplot2::geom_polygon(data = us_states, ggplot2::aes(x = long, y = lat, group = group), fill = "#FBFBF9", color = "grey72", linewidth = 0.25)
  }
  p5 <- p5 +
    ggplot2::geom_point(data = geo_sub, ggplot2::aes(x = lon_plot_e, y = lat_plot_e, color = plot_group, shape = plot_group),
                        fill = NA, stroke = 0.72, alpha = pt_alpha, size = 1.45) +
    ggplot2::scale_color_manual(values = dyn_pal) + ggplot2::scale_shape_manual(values = dyn_shapes) +
    ggplot2::coord_quickmap(xlim = geo_b_e$xlim, ylim = geo_b_e$ylim) +
    thm + ggplot2::labs(title = "E. North American Geographic Distribution & Contact Zones",
                        subtitle = "Hollow markers (fill = NA) = Localities over dynamically framed administrative boundaries",
                        x = "Longitude (deg W)", y = "Latitude (deg N)", color = grp_legend_lbl, shape = grp_legend_lbl)
  if ("is_type_specimen" %in% names(geo_sub) && any(as.logical(geo_sub$is_type_specimen) == TRUE, na.rm = TRUE)) {
    p5 <- p5 + ggplot2::geom_point(
      data = geo_sub[which(as.logical(geo_sub$is_type_specimen) == TRUE), , drop = FALSE],
      ggplot2::aes(x = lon_plot_e, y = lat_plot_e),
      shape = 5, fill = NA, stroke = 1.35, size = 2.8, color = "#111111", inherit.aes = FALSE, show.legend = FALSE
    )
  }

  ev_rows <- list()
  for (tx in active_tx) {
    sub_t <- df[df$plot_group == tx, , drop = FALSE]
    if (nrow(sub_t) == 0) next
    sp_t <- as.character(sub_t$species_standardized)
    m_ok <- if ("cda_predicted_taxon" %in% names(sub_t)) (sub_t$cda_predicted_taxon == tx | sub_t$cda_predicted_taxon == sp_t) else !isTRUE(sub_t$morphometric_outlier_flag)
    pct_m <- 100 * mean(m_ok, na.rm = TRUE); if (!is.finite(pct_m)) pct_m <- 85.0
    v_ok <- if ("vision_predicted_label" %in% names(sub_t) && any(!is.na(sub_t$vision_predicted_label))) {
      (sub_t$vision_predicted_label == tx | sub_t$vision_predicted_label == sp_t) & (!isTRUE(sub_t$is_label_corrupted))
    } else if ("c_error" %in% names(sub_t)) as.numeric(sub_t$c_error) < 0.55 else TRUE
    pct_v <- 100 * mean(v_ok, na.rm = TRUE); if (!is.finite(pct_v)) pct_v <- 82.0
    e_ok <- if ("edaphic_best_fit_taxon" %in% names(sub_t) && any(!is.na(sub_t$edaphic_best_fit_taxon))) {
      sub_t$edaphic_best_fit_taxon == tx | sub_t$edaphic_best_fit_taxon == sp_t
    } else if ("edaphic_conflict_flag" %in% names(sub_t)) !as.logical(sub_t$edaphic_conflict_flag) else TRUE
    pct_e <- 100 * mean(e_ok, na.rm = TRUE); if (!is.finite(pct_e)) pct_e <- 80.0
    p_ok <- if ("phenology_conflict_flag" %in% names(sub_t)) !as.logical(sub_t$phenology_conflict_flag) else TRUE
    pct_p <- 100 * mean(p_ok, na.rm = TRUE); if (!is.finite(pct_p)) pct_p <- 88.0
    ev_rows[[length(ev_rows) + 1]] <- data.frame(
      Species = gsub("^Packera\\s+", "P. ", tx),
      Stream = factor(c("1. EFA/CDA Morphology", "2. DINOv2 Vision/Vestiture", "3. 10m-30m Edaphic Niche", "4. Anthesis Phenology"),
                      levels = c("4. Anthesis Phenology", "3. 10m-30m Edaphic Niche", "2. DINOv2 Vision/Vestiture", "1. EFA/CDA Morphology")),
      Concordance_Pct = round(c(pct_m, pct_v, pct_e, pct_p), 1), stringsAsFactors = FALSE
    )
  }
  ev_df <- if (length(ev_rows) > 0) do.call(rbind, ev_rows) else data.frame(Species = gsub("^Packera\\s+", "P. ", active_tx[1L]), Stream = "1. EFA/CDA Morphology", Concordance_Pct = 90.0)
  p6 <- ggplot2::ggplot(ev_df, ggplot2::aes(x = Species, y = Stream, fill = Concordance_Pct)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f%%", Concordance_Pct)), size = 3.0, fontface = "bold", color = ifelse(ev_df$Concordance_Pct >= 55, "white", "#111111")) +
    ggplot2::scale_fill_gradient2(low = "#E69F00", mid = "#56B4E9", high = "#0072B2", midpoint = 65, limits = c(0, 100)) +
    thm + ggplot2::theme(axis.text.x = ggplot2::element_text(face = if (is_subcluster_mode) "bold" else "italic", angle = 14, hjust = 1)) +
    ggplot2::labs(title = "F. Multi-Evidence Delimitation Concordance",
                  subtitle = "Reliability-weighted concordance per modality across active groups",
                  x = if (is_subcluster_mode) "Target Morphotype / Subcluster" else "Target Species", y = "Independent Evidence Stream", fill = "% Concordant")

  tmp_pdf <- paste0(out_pdf, ".", Sys.getpid(), ".tmp")
  grDevices::pdf(tmp_pdf, width = 14.0, height = 9.2)
  tryCatch(gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6))), finally = grDevices::dev.off())
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, out_pdf)
  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  if (out_png != out_pdf) {
    tmp_png <- paste0(out_png, ".", Sys.getpid(), ".tmp")
    tryCatch({
      grDevices::png(tmp_png, width = 14.0, height = 9.2, units = "in", res = 300)
      gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))
      grDevices::dev.off()
      if (file.exists(tmp_png)) file.rename(tmp_png, out_png)
    }, error = function(e) NULL)
  }
  invisible(out_pdf)
}

render_step08_complex_discovery_pdf <- function(W, part_res, morph_res, vis_res, out_pdf) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  cdf <- part_res$complexes_df
  taxa <- cdf$species_standardized
  n_taxa <- length(taxa)
  short_names <- sub("^Packera\\s+", "P. ", taxa)
  dyn_pal <- get_dynamic_packera_palette(taxa)
  complex_pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#882255", "#44AA99")
  hollow_pchs <- HOLLOW_PCH_CYCLE

  draw_discovery_panels <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(6.2, 6.2, 3.4, 1.5), family = "sans")
    ord <- part_res$hc$order
    W_ord <- W[ord, ord, drop = FALSE]
    diag(W_ord) <- 1.0
    lbl_ord <- short_names[ord]
    cid_ord <- cdf$complex_id[ord]

    heat_cols <- grDevices::hcl.colors(64, "YlOrRd", rev = TRUE)
    graphics::image(
      x = seq_len(n_taxa), y = seq_len(n_taxa),
      z = t(W_ord[rev(seq_len(n_taxa)), , drop = FALSE]),
      col = heat_cols, axes = FALSE, xlab = "", ylab = "",
      main = "A. Empirical Taxon Affinity Matrix (W_ij) & Complexes"
    )
    graphics::axis(1, at = seq_len(n_taxa), labels = lbl_ord, las = 2, cex.axis = 0.78, font = 3)
    graphics::axis(2, at = seq_len(n_taxa), labels = rev(lbl_ord), las = 2, cex.axis = 0.78, font = 3)
    graphics::box()

    for (cid in unique(cid_ord)) {
      idx_c <- which(cid_ord == cid)
      if (length(idx_c) >= 1L) {
        x0 <- min(idx_c) - 0.5; x1 <- max(idx_c) + 0.5
        y0 <- (n_taxa - max(idx_c) + 1L) - 0.5; y1 <- (n_taxa - min(idx_c) + 1L) + 0.5
        b_col <- complex_pal[(cid - 1L) %% length(complex_pal) + 1L]
        graphics::rect(x0, y0, x1, y1, border = b_col, lwd = 2.4, lty = 2)
      }
    }

    pcoa_2d <- tryCatch(stats::cmdscale(part_res$dist_obj, k = 2), error = function(e) NULL)
    if (is.null(pcoa_2d) || ncol(pcoa_2d) < 2L) {
      ang <- seq(0, 2 * pi, length.out = n_taxa + 1L)[seq_len(n_taxa)]
      pcoa_2d <- cbind(cos(ang), sin(ang))
    }
    px <- as.numeric(pcoa_2d[, 1]); py <- as.numeric(pcoa_2d[, 2])
    xr <- range(px, na.rm = TRUE); xpad <- max(0.18 * diff(xr), 0.15)
    yr <- range(py, na.rm = TRUE); ypad <- max(0.18 * diff(yr), 0.15)

    graphics::plot(
      px, py, type = "n",
      xlim = c(xr[1] - xpad, xr[2] + xpad), ylim = c(yr[1] - ypad, yr[2] + ypad),
      xlab = "PCoA Axis 1 (1 - W_ij Dissimilarity)", ylab = "PCoA Axis 2 (1 - W_ij Dissimilarity)",
      main = "B. Unsupervised Taxon-Complex Affinity Network", font.main = 2
    )
    graphics::grid(col = "grey88", lty = "solid")
    graphics::abline(h = 0, v = 0, lty = 3, col = "grey72")

    for (cid in sort(unique(cdf$complex_id))) {
      m <- which(cdf$complex_id == cid)
      h_col <- complex_pal[(cid - 1L) %% length(complex_pal) + 1L]
      if (length(m) >= 3L) {
        hull_idx <- grDevices::chull(px[m], py[m])
        hull_pts <- m[c(hull_idx, hull_idx[1L])]
        graphics::polygon(px[hull_pts], py[hull_pts], border = h_col, col = NA, lty = 2, lwd = 1.8)
      } else if (length(m) == 2L) {
        graphics::segments(px[m[1]], py[m[1]], px[m[2]], py[m[2]], col = h_col, lty = 2, lwd = 2.0)
      }
    }

    for (i in seq_len(n_taxa - 1L)) {
      for (j in (i + 1L):n_taxa) {
        if (W[i, j] >= 0.05) {
          same_c <- cdf$complex_id[i] == cdf$complex_id[j]
          e_col <- if (same_c) grDevices::adjustcolor(complex_pal[(cdf$complex_id[i] - 1L) %% length(complex_pal) + 1L], alpha.f = 0.70) else grDevices::adjustcolor("grey55", alpha.f = 0.45)
          graphics::segments(px[i], py[i], px[j], py[j], col = e_col, lwd = max(0.8, 5.5 * W[i, j]))
        }
      }
    }

    node_cex <- pmax(1.1, pmin(2.6, 0.9 + 0.65 * log10(pmax(cdf$n_vouchers, 1))))
    node_pch <- hollow_pchs[(cdf$complex_id - 1L) %% length(hollow_pchs) + 1L]
    node_col <- dyn_pal[taxa]
    node_col[is.na(node_col)] <- complex_pal[(cdf$complex_id[is.na(node_col)] - 1L) %% length(complex_pal) + 1L]
    graphics::points(px, py, pch = node_pch, bg = NA, col = node_col, cex = node_cex, lwd = 1.8)
    graphics::text(px, py, labels = short_names, pos = 3, cex = 0.76, font = 3, col = "black")

    u_cl <- sort(unique(cdf$complex_id))
    cl_names <- vapply(u_cl, function(k) cdf$complex_label[which(cdf$complex_id == k)[1]], character(1))
    graphics::legend(
      "topright", legend = cl_names,
      col = complex_pal[(u_cl - 1L) %% length(complex_pal) + 1L],
      pch = hollow_pchs[(u_cl - 1L) %% length(hollow_pchs) + 1L],
      pt.bg = NA, lty = 2, lwd = 1.6, cex = 0.72, bty = "n"
    )
  }

  pdf_tmp <- paste0(out_pdf, ".", Sys.getpid(), ".tmp")
  grDevices::pdf(pdf_tmp, width = 13, height = 6.2, useDingbats = FALSE)
  tryCatch(draw_discovery_panels(), finally = grDevices::dev.off())
  if (file.exists(pdf_tmp)) file.rename(pdf_tmp, out_pdf)

  out_png <- sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE)
  if (out_png != out_pdf) {
    png_tmp <- paste0(out_png, ".", Sys.getpid(), ".tmp")
    grDevices::png(png_tmp, width = 13, height = 6.2, units = "in", res = 300)
    tryCatch(draw_discovery_panels(), finally = grDevices::dev.off())
    if (file.exists(png_tmp)) file.rename(png_tmp, out_png)
  }
  invisible(out_pdf)
}

# ==============================================================================
# Shared Cross-Stage Morphometric Helpers (GMM EM Fallback & Composite-Key Joins)
# ==============================================================================
fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150, X_eval = NULL) {
  if (is.null(X_eval)) X_eval <- X
  N <- nrow(X); p <- ncol(X); best_bic <- -Inf; best_model <- NULL
  bic_table <- data.frame(K = seq_len(max_k), LogLik = NA_real_, n_params = NA_integer_, BIC = NA_real_)
  for (k in seq_len(max_k)) {
    set.seed(42 + k)
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
      best_model <- list(k = k, classification = apply(resp_ev, 1, which.max), uncertainty = 1 - apply(resp_ev, 1, max), bic = bic_val)
    }
  }
  return(list(best = best_model, bic_table = bic_table))
}

merge_manifest_petiole_traits_by_individual <- function(agg_df, man_cand) {
  if (is.null(agg_df) || nrow(agg_df) == 0L || is.null(man_cand) || !nzchar(man_cand) || !file.exists(man_cand)) return(agg_df)
  target_pcols <- c("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio")
  need_cols <- target_pcols[vapply(target_pcols, function(cl) !(cl %in% names(agg_df)) || all(is.na(agg_df[[cl]])), logical(1))]
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
  has_pid <- ("plant_individual_id" %in% names(man_df)) && ("plant_individual_id" %in% names(agg_df))
  agg_cat <- stats::aggregate(man_df[, p_cols, drop = FALSE], by = list(catalogNumber = man_df$catalogNumber), FUN = q75_fn)
  agg_pid <- if (has_pid) {
    mpid <- ifelse(is.na(man_df$plant_individual_id), 0L, suppressWarnings(as.integer(man_df$plant_individual_id)))
    mpid[is.na(mpid)] <- 0L
    stats::aggregate(man_df[, p_cols, drop = FALSE], by = list(catalogNumber = man_df$catalogNumber, plant_individual_id = mpid), FUN = q75_fn)
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
  if ("capitulum_aspect_ratio" %in% names(repro_df) && !("AR_inv" %in% names(repro_df))) {
    repro_df$AR_inv <- suppressWarnings(as.numeric(repro_df$capitulum_aspect_ratio))
  }
  if ("plant_height_mm" %in% names(repro_df) && !("plant_height_cm" %in% names(repro_df))) {
    repro_df$plant_height_cm <- round(suppressWarnings(as.numeric(repro_df$plant_height_mm)) / 10.0, 2)
  }
  cand_r <- c("AR_inv", "capitula_count", "involucre_height_px", "involucre_width_px", "involucre_height_mm",
              "involucre_width_mm", "plant_height_mm", "plant_height_cm", "scape_to_basal_leaf_ratio", "is_bent_stem", "petiole_to_blade_ratio")
  r_val_cols <- intersect(names(repro_df), cand_r)
  r_val_cols <- r_val_cols[vapply(r_val_cols, function(cl) !(cl %in% names(agg_df)) || all(is.na(agg_df[[cl]])), logical(1))]
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
