if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
# ==============================================================================
# Module: scripts/figures/step05_sdm_figures.R
# Step 05 Equal-Area Projected MaxEnt/Ensemble SDM Suitability Figures:
#   fig05c_sdm_suitability_and_distribution_map.png
# ==============================================================================

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

  raster_to_proj_df <- function(r_obj, max_cells = 45000L, target_crs = "EPSG:5070") {
    if (is.null(r_obj) || !requireNamespace("terra", quietly = TRUE)) return(list(df = data.frame(), projected = FALSE))
    tryCatch({
      r_use <- if (terra::ncell(r_obj) > max_cells) {
        agg_f <- max(2L, ceiling(sqrt(terra::ncell(r_obj) / max_cells)))
        terra::aggregate(r_obj, fact = agg_f, fun = "mean", na.rm = TRUE)
      } else r_obj
      is_proj <- FALSE
      r_out <- r_use
      if (requireNamespace("sf", quietly = TRUE)) {
        if (!nzchar(terra::crs(r_use))) terra::crs(r_use) <- "EPSG:4326"
        r_tmp <- tryCatch(terra::project(r_use, target_crs, method = "bilinear"), error = function(e) NULL)
        if (!is.null(r_tmp) && terra::ncell(r_tmp) > 0L) {
          d_chk <- as.data.frame(r_tmp, xy = TRUE, na.rm = TRUE)
          if (ncol(d_chk) >= 3L && nrow(d_chk) > 0L && any(is.finite(d_chk[[3]]))) {
            r_out <- r_tmp; is_proj <- TRUE
          }
        }
      }
      d_r <- as.data.frame(r_out, xy = TRUE, na.rm = TRUE)
      if (ncol(d_r) < 3L || nrow(d_r) == 0L) return(list(df = data.frame(), projected = FALSE))
      names(d_r)[1:3] <- c("x", "y", "suitability")
      d_r <- d_r[is.finite(d_r$suitability), c("x", "y", "suitability"), drop = FALSE]
      list(df = d_r, projected = is_proj)
    }, error = function(e) list(df = data.frame(), projected = FALSE))
  }

  states_5070 <- if (requireNamespace("sf", quietly = TRUE) && requireNamespace("maps", quietly = TRUE)) {
    tryCatch({
      m_usa <- maps::map("state", plot = FALSE, fill = TRUE)
      sf_u <- sf::st_as_sf(m_usa)
      sf::st_crs(sf_u) <- 4326
      sf::st_transform(sf::st_make_valid(sf_u), 5070L)
    }, error = function(e) NULL)
  } else NULL

  # Helper to build a projected Albers Equal Area (EPSG:5070) suitability + hollow occurrence map
  build_suitability_panel <- function(r_info, pts_sub, title_txt, sub_txt = "") {
    r_df <- if (is.list(r_info) && is.data.frame(r_info$df)) r_info$df else data.frame()
    is_proj <- is.list(r_info) && isTRUE(r_info$projected) && !is.null(states_5070)
    if (nrow(r_df) > 0L) {
      pts_ok_s <- data.frame()
      if (nrow(pts_sub) > 0L) {
        pts_sub$lon_num <- suppressWarnings(as.numeric(pts_sub$decimalLongitude))
        pts_sub$lat_num <- suppressWarnings(as.numeric(pts_sub$decimalLatitude))
        pts_ok_s <- pts_sub[is.finite(pts_sub$lon_num) & is.finite(pts_sub$lat_num), , drop = FALSE]
      }
      if (is_proj && requireNamespace("sf", quietly = TRUE)) {
        p_m <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = r_df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$suitability), alpha = 0.92) +
          ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Suitability") +
          ggplot2::geom_sf(data = states_5070, fill = NA, color = "#E5E7EB", linewidth = 0.28, inherit.aes = FALSE)
        if (nrow(pts_ok_s) > 0L) {
          xy_pts <- tryCatch(sf::st_coordinates(sf::st_transform(sf::st_as_sf(pts_ok_s, coords = c("lon_num", "lat_num"), crs = 4326), 5070L)), error = function(e) NULL)
          if (!is.null(xy_pts) && nrow(xy_pts) == nrow(pts_ok_s)) {
            pts_ok_s$px <- xy_pts[, 1L]; pts_ok_s$py <- xy_pts[, 2L]
            p_m <- p_m + ggplot2::geom_point(
              data = pts_ok_s, ggplot2::aes(x = .data$px, y = .data$py),
              shape = 1, size = 2.0, stroke = 0.8, color = "#FF6B00", alpha = 0.88, inherit.aes = FALSE
            )
          }
        }
        return(
          p_m +
            ggplot2::coord_sf(crs = 5070L, xlim = range(r_df$x), ylim = range(r_df$y), expand = FALSE) +
            theme_packera_pub(base_size = 9.5) +
            ggplot2::labs(title = title_txt, subtitle = sub_txt, x = "Longitude (EPSG:5070)", y = "Latitude (EPSG:5070)")
        )
      }
      cos_l <- cos(mean(range(r_df$y, na.rm = TRUE)) * pi / 180.0)
      p_m <- ggplot2::ggplot() +
        ggplot2::geom_raster(data = r_df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$suitability), alpha = 0.90) +
        ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Suitability")
      if (requireNamespace("maps", quietly = TRUE)) {
        st_d <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)
        if (!is.null(st_d)) p_m <- p_m + ggplot2::geom_path(data = st_d, ggplot2::aes(x = .data$long, y = .data$lat, group = .data$group), color = "#E5E7EB", linewidth = 0.25)
      }
      if (nrow(pts_ok_s) > 0L) {
        p_m <- p_m + ggplot2::geom_point(data = pts_ok_s, ggplot2::aes(x = .data$lon_num, y = .data$lat_num), shape = 1, size = 2.0, stroke = 0.8, color = "#FF6B00", alpha = 0.88)
      }
      p_m +
        ggplot2::coord_fixed(ratio = 1.0 / max(0.25, cos_l), xlim = range(r_df$x), ylim = range(r_df$y), expand = FALSE) +
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
    .fmt_sdm_strip <- function(v) wrap_clean_label_r(gsub("_", " ", abbreviate_icn_taxon(v, abbreviate_genus = TRUE)), width = 22L)
    strip_face_m <- if (any(grepl("^(Ecotype|Morphotype|Complex|Cluster)", names(multi_rasters)))) "bold" else "bold.italic"
    grid_infos <- lapply(names(multi_rasters), function(tx_nm) {
      ri <- raster_to_proj_df(multi_rasters[[tx_nm]], max_cells = 18000L)
      if (nrow(ri$df) > 0L) ri$df$taxon_lbl <- .fmt_sdm_strip(tx_nm)
      ri
    })
    proj_multi_ok <- all(vapply(grid_infos, function(ri) nrow(ri$df) == 0L || isTRUE(ri$projected), logical(1))) && !is.null(states_5070)
    grid_all <- do.call(rbind, lapply(grid_infos, function(ri) ri$df))
    if (!is.null(grid_all) && nrow(grid_all) > 0L) {
      pts_ok <- vdf[is.finite(suppressWarnings(as.numeric(vdf$decimalLongitude))) & is.finite(suppressWarnings(as.numeric(vdf$decimalLatitude))), , drop = FALSE]
      pts_ok$lon_num <- suppressWarnings(as.numeric(pts_ok$decimalLongitude))
      pts_ok$lat_num <- suppressWarnings(as.numeric(pts_ok$decimalLatitude))
      pts_ok$taxon_lbl <- .fmt_sdm_strip(pts_ok$plot_group)
      pts_ok <- pts_ok[pts_ok$taxon_lbl %in% unique(grid_all$taxon_lbl), , drop = FALSE]
      n_sp_m <- length(unique(grid_all$taxon_lbl))
      ncol_m <- min(4L, max(2L, ceiling(sqrt(n_sp_m))))
      n_rows_multi <- as.integer(ceiling(n_sp_m / ncol_m))
      if (proj_multi_ok && requireNamespace("sf", quietly = TRUE)) {
        if (nrow(pts_ok) > 0L) {
          xy_p <- tryCatch(sf::st_coordinates(sf::st_transform(sf::st_as_sf(pts_ok, coords = c("lon_num", "lat_num"), crs = 4326), 5070L)), error = function(e) NULL)
          if (!is.null(xy_p) && nrow(xy_p) == nrow(pts_ok)) {
            pts_ok$px <- xy_p[, 1L]; pts_ok$py <- xy_p[, 2L]
          } else {
            pts_ok <- pts_ok[0L, , drop = FALSE]
          }
        }
        p_multi_grid <- ggplot2::ggplot(grid_all, ggplot2::aes(x = .data$x, y = .data$y)) +
          ggplot2::geom_raster(ggplot2::aes(fill = .data$suitability)) +
          ggplot2::geom_sf(data = states_5070, fill = NA, color = "#E5E7EB", linewidth = 0.22, inherit.aes = FALSE) +
          ggplot2::geom_point(
            data = pts_ok,
            ggplot2::aes(x = .data$px, y = .data$py),
            shape = 1, size = 1.35, stroke = 0.65, color = "#FF6B00", alpha = 0.80, inherit.aes = FALSE
          ) +
          ggplot2::facet_wrap(~ taxon_lbl, ncol = ncol_m) +
          ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 1), name = "Ensemble\nSuitability") +
          ggplot2::coord_sf(crs = 5070L, xlim = range(grid_all$x), ylim = range(grid_all$y), expand = FALSE) +
          theme_packera_pub(base_size = 9.2) +
          ggplot2::theme(strip.text = ggplot2::element_text(face = strip_face_m, size = 7.6, lineheight = 0.92)) +
          ggplot2::labs(
            title = sprintf("E. Comparative Multi-Taxon Ensemble Habitat Suitability (%d Modeled Taxa, EPSG:5070)", n_sp_m),
            subtitle = "Orange hollow circles mark coordinate-uncertainty-gated herbarium presences",
            x = "Longitude (EPSG:5070)", y = "Latitude (EPSG:5070)"
          )
      } else {
        cos_m <- cos(mean(range(grid_all$y, na.rm = TRUE)) * pi / 180.0)
        p_multi_grid <- ggplot2::ggplot(grid_all, ggplot2::aes(x = .data$x, y = .data$y)) +
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
          ggplot2::theme(strip.text = ggplot2::element_text(face = strip_face_m, size = 7.6, lineheight = 0.92)) +
          ggplot2::labs(
            title = sprintf("E. Comparative Multi-Taxon Ensemble Habitat Suitability (%d Modeled Taxa)", n_sp_m),
            subtitle = "Orange hollow circles mark coordinate-uncertainty-gated herbarium presences",
            x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)"
          )
      }
    }
  }

  f_raw <- trimws(as.character(focal_taxon))
  f_tx <- if (nzchar(f_raw)) {
    f_raw
  } else {
    active_groups[1L]
  }
  pts_f <- if (has_focal_sdm && f_tx %in% vdf$plot_group) {
    vdf[vdf$plot_group == f_tx, , drop = FALSE]
  } else if (has_focal_sdm && ("species_standardized" %in% names(vdf)) && (f_tx %in% vdf$species_standardized)) {
    vdf[vdf$species_standardized == f_tx, , drop = FALSE]
  } else {
    vdf
  }
  info_max <- if (has_focal_sdm) raster_to_proj_df(focal_rasters$maxent) else list(df = data.frame(), projected = FALSE)
  info_ens <- if (has_focal_sdm) raster_to_proj_df(focal_rasters$ensemble) else list(df = data.frame(), projected = FALSE)

  if (nrow(info_max$df) > 0L || nrow(info_ens$df) > 0L) {
    f_tx_disp <- gsub("_", " ", f_tx)
    p_a <- build_suitability_panel(if (nrow(info_max$df) > 0L) info_max else info_ens, pts_f, sprintf("A. MaxEnt (cloglog) Suitability: %s", f_tx_disp), "L1-regularized IPP presence-background model (EPSG:5070)")
    p_b <- build_suitability_panel(if (nrow(info_ens$df) > 0L) info_ens else info_max, pts_f, sprintf("B. Spatial RF + MaxEnt Ensemble: %s", f_tx_disp), "AUC-weighted ensemble consensus surface (EPSG:5070)")

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
          title = "Step 05: Ensemble Species Distribution Modeling, Spatial-Block Validation & Multi-Taxon Suitability",
          theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4))
        )
      return(save_packera_png_300(plate_05c, file.path(fdir, "fig05c_sdm_suitability_and_distribution_map.png"), width = 14.2, height = 9.4 + 2.7 * n_rows_multi))
    }

    plate_05c <- ((p_a | p_b) / (p_c | p_d)) +
      patchwork::plot_layout(widths = c(1.1, 1.1), heights = c(1.25, 0.95)) +
      patchwork::plot_annotation(
        title = "Step 05: Ensemble Species Distribution Modeling & Spatial Block Validation",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4))
      )
    return(save_packera_png_300(plate_05c, file.path(fdir, "fig05c_sdm_suitability_and_distribution_map.png"), width = 14.0, height = 9.8))
  }

  # Geographic Occurrence & Ecoregion Map when raster layers are absent
  p_geo <- build_suitability_panel(data.frame(), vdf, "Step 05: Equal-Area Projected Geographic Distribution (Albers EPSG:5070)", "Hollow markers stratified by taxon; scale bar & north arrow in projected CRS")
  save_packera_png_300(p_geo, file.path(fdir, "fig05c_sdm_suitability_and_distribution_map.png"), width = 11.8, height = 8.2)
}
