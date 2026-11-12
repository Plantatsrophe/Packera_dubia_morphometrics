# ==============================================================================
# Module: scripts/figures/theme_packera_pub.R
# Shared 300-DPI PNG Publication Theme, Combinatorial 60-Taxon Hollow Scales,
# ICN Rank-Aware Typography, Equal-Area sf Basemaps & Native LM3 Asset Loader
# ==============================================================================

PACKERA_TOL_OKABE_12 <- c(
  "#0072B2", # Blue
  "#D55E00", # Vermilion
  "#009E73", # Bluish green
  "#CC79A7", # Reddish purple
  "#E69F00", # Orange
  "#56B4E9", # Sky blue
  "#332288", # Indigo
  "#882285", # Wine
  "#44AA99", # Teal
  "#999933", # Olive
  "#666666", # Slate
  "#117733"  # Forest green
)

# Strictly unfilled/hollow point geometries: circle, triangle, diamond, square, inverted-triangle
PACKERA_HOLLOW_SHAPES_5 <- c(1L, 2L, 5L, 0L, 6L)

PACKERA_ICN_ROMAN_TOKENS <- c(
  "subsp.", "subsp", "ssp.", "ssp", "var.", "var", "f.", "forma",
  "x", "\u00d7", "aff.", "aff", "cf.", "cf", "sp.", "sp", "s.s.", "s.l.", "ined."
)

normalize_png_path <- function(out_path) {
  if (is.null(out_path) || !nzchar(trimws(as.character(out_path)))) return("")
  p <- trimws(as.character(out_path))
  if (grepl("\\.pdf$", p, ignore.case = TRUE)) {
    p <- sub("\\.pdf$", ".png", p, ignore.case = TRUE)
  } else if (!grepl("\\.png$", p, ignore.case = TRUE)) {
    p <- paste0(p, ".png")
  }
  p
}

save_packera_png_300 <- function(plot_or_fn, out_path, width = 12.5, height = 8.0, dpi = 300L, also_copy_to = NULL) {
  png_final <- normalize_png_path(out_path)
  if (!nzchar(png_final)) return(invisible(""))
  fdir <- dirname(png_final)
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  stem <- tools::file_path_sans_ext(basename(png_final))
  png_tmp <- file.path(fdir, sprintf(".%s_%d.png.tmp", stem, Sys.getpid()))
  stale_pdf <- file.path(fdir, paste0(stem, ".pdf"))

  open_png_device <- function(target_file) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(filename = target_file, width = width, height = height, units = "in", res = dpi, background = "white")
      return(TRUE)
    }
    ok_cairo <- tryCatch({
      grDevices::png(filename = target_file, width = width, height = height, units = "in", res = dpi, bg = "white", type = "cairo")
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(ok_cairo)) {
      grDevices::png(filename = target_file, width = width, height = height, units = "in", res = dpi, bg = "white")
    }
    TRUE
  }

  ok_render <- tryCatch({
    open_png_device(png_tmp)
    on.exit({
      if (grDevices::dev.cur() > 1L) tryCatch(grDevices::dev.off(), error = function(e) NULL)
    }, add = TRUE)
    if (is.function(plot_or_fn)) {
      res_p <- plot_or_fn()
      if (inherits(res_p, "gg") || inherits(res_p, "ggplot") || inherits(res_p, "patchwork")) {
        print(res_p)
      }
    } else {
      print(plot_or_fn)
    }
    grDevices::dev.off()
    on.exit(NULL)
    TRUE
  }, error = function(e) {
    warning("save_packera_png_300 warning (", basename(png_final), "): ", e$message)
    FALSE
  })

  if (isTRUE(ok_render) && file.exists(png_tmp) && file.info(png_tmp)$size > 0) {
    file.rename(png_tmp, png_final)
    if (file.exists(stale_pdf)) unlink(stale_pdf, force = TRUE)
    if (!is.null(also_copy_to)) {
      for (cp_raw in as.character(also_copy_to)) {
        cp_png <- normalize_png_path(cp_raw)
        if (nzchar(cp_png) && normalizePath(cp_png, mustWork = FALSE) != normalizePath(png_final, mustWork = FALSE)) {
          dir.create(dirname(cp_png), recursive = TRUE, showWarnings = FALSE)
          cp_tmp <- paste0(cp_png, ".", Sys.getpid(), ".tmp")
          if (isTRUE(file.copy(png_final, cp_tmp, overwrite = TRUE))) file.rename(cp_tmp, cp_png)
          cp_stale_pdf <- sub("\\.png$", ".pdf", cp_png, ignore.case = TRUE)
          if (cp_stale_pdf != cp_png && file.exists(cp_stale_pdf)) unlink(cp_stale_pdf, force = TRUE)
        }
      }
    }
  } else if (file.exists(png_tmp)) {
    unlink(png_tmp, force = TRUE)
  }
  invisible(png_final)
}

resolve_taxon_or_rank_vector <- function(df, group_col = NULL) {
  n <- nrow(df)
  if (n == 0L) return(character(0))
  base_v <- if (!is.null(group_col) && group_col %in% names(df)) {
    as.character(df[[group_col]])
  } else if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("scientificName" %in% names(df)) {
    as.character(df$scientificName)
  } else {
    rep("Unknown", n)
  }
  base_v[is.na(base_v) | !nzchar(trimws(base_v))] <- "Unknown"
  base_v <- trimws(base_v)

  # Append infraspecific rank (var. / subsp. / f.) when column is populated and distinct from species epithet
  if ("infraspecific_taxon" %in% names(df) && (is.null(group_col) || group_col %in% c("species_standardized", "scientificName"))) {
    inf_v <- trimws(as.character(df$infraspecific_taxon))
    has_inf <- !is.na(inf_v) & nzchar(inf_v) & !(tolower(inf_v) %in% c("na", "none", "null", "unknown"))
    if (any(has_inf)) {
      for (i in which(has_inf)) {
        b_i <- base_v[i]
        if (!grepl("^Packera\\s+", b_i) || grepl("\\b(subsp\\.|var\\.|f\\.)\\s+", b_i, ignore.case = TRUE)) next
        iv <- inf_v[i]
        if (grepl("^(subsp\\.?|ssp\\.?|var\\.?|f\\.?|forma)\\s+", iv, ignore.case = TRUE)) {
          iv_clean <- sub("^(subsp|ssp)\\.?\\s+", "subsp. ", iv, ignore.case = TRUE)
          iv_clean <- sub("^var\\.?\\s+", "var. ", iv_clean, ignore.case = TRUE)
          iv_clean <- sub("^(f|forma)\\.?\\s+", "f. ", iv_clean, ignore.case = TRUE)
          base_v[i] <- paste(b_i, iv_clean)
        } else if (grepl("^[a-z][a-z-]+$", tolower(iv))) {
          ep <- strsplit(b_i, "\\s+")[[1L]]
          if (length(ep) >= 2L && tolower(iv) != tolower(ep[2L])) {
            base_v[i] <- paste0(b_i, " var. ", tolower(iv))
          }
        }
      }
    }
  }
  base_v
}

abbreviate_icn_taxon <- function(label_vec, abbreviate_genus = FALSE) {
  vapply(as.character(label_vec), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("Unknown")
    s_clean <- trimws(s)
    if (!isTRUE(abbreviate_genus)) return(s_clean)
    if (grepl("^(Morphotype|Ecotype|Complex|Cluster)_", s_clean)) return(gsub("_", " ", s_clean))
    toks <- strsplit(s_clean, "\\s+")[[1L]]
    if (length(toks) >= 2L && grepl("^[A-Z][a-z]+$", toks[1L])) {
      toks[1L] <- paste0(substr(toks[1L], 1L, 1L), ".")
      return(paste(toks, collapse = " "))
    }
    s_clean
  }, character(1), USE.NAMES = FALSE)
}

format_icn_label_expr <- function(tx, abbreviate_genus = FALSE, count_n = NULL) {
  s <- abbreviate_icn_taxon(tx, abbreviate_genus = abbreviate_genus)
  suffix_n <- if (!is.null(count_n) && is.finite(count_n)) sprintf(" (n=%d)", as.integer(count_n)) else ""
  if (grepl("^(Morphotype|Ecotype|Complex|Cluster|Unknown|Unassigned)", s, ignore.case = TRUE) ||
      !grepl("^[A-Z][a-z.]*\\s+", s)) {
    lbl_plain <- paste0(gsub("_", " ", s), suffix_n)
    return(bquote(.(lbl_plain)))
  }

  toks <- strsplit(s, "\\s+")[[1L]]
  expr_parts <- character(0)
  for (i in seq_along(toks)) {
    tk <- toks[i]
    low <- tolower(tk)
    esc <- gsub('"', '\\"', tk, fixed = TRUE)
    if (i == 1L) {
      expr_parts <- c(expr_parts, sprintf('italic("%s")', esc))
    } else if (low %in% PACKERA_ICN_ROMAN_TOKENS) {
      rom <- if (low %in% c("x", "\u00d7")) "\u00d7"
             else if (low %in% c("subsp", "ssp", "subsp.", "ssp.")) "subsp."
             else if (low %in% c("var", "var.")) "var."
             else if (low %in% c("f", "f.", "forma")) "f."
             else esc
      expr_parts <- c(expr_parts, sprintf('"%s"', rom))
    } else if (grepl("^[a-z][a-z-]+$", tk)) {
      expr_parts <- c(expr_parts, sprintf('italic("%s")', esc))
    } else {
      expr_parts <- c(expr_parts, sprintf('"%s"', esc))
    }
  }
  if (nzchar(suffix_n)) {
    expr_parts <- c(expr_parts, sprintf('"%s"', trimws(suffix_n)))
  }
  code_str <- paste(expr_parts, collapse = " ~ ")
  tryCatch(parse(text = code_str)[[1L]], error = function(e) paste0(s, suffix_n))
}

resolve_plot_group_labels <- function(df, cda_res = NULL, gmm_res = NULL, group_col = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  if (!is.null(group_col) && is.character(group_col) && length(group_col) == 1L && group_col %in% names(df)) {
    return(resolve_taxon_or_rank_vector(df, group_col = group_col))
  }
  if (exists("resolve_plot_group_column", mode = "function")) {
    gv <- tryCatch(resolve_plot_group_column(df, cda_res = cda_res, gmm_res = gmm_res), error = function(e) NULL)
    if (is.list(gv) && !is.null(gv$col)) {
      return(resolve_taxon_or_rank_vector(df, group_col = gv$col))
    }
    if (is.character(gv) && length(gv) == nrow(df)) {
      gv_clean <- trimws(gv)
      gv_clean[is.na(gv_clean) | !nzchar(gv_clean)] <- "Unknown"
      non_unk <- gv_clean[gv_clean != "Unknown"]
      if (length(non_unk) > 0L && all(grepl("^Packera\\s+", non_unk)) && "infraspecific_taxon" %in% names(df)) {
        tmp_df <- df
        tmp_df$.__tmp_base_sp <- gv_clean
        return(resolve_taxon_or_rank_vector(tmp_df, group_col = ".__tmp_base_sp"))
      }
      return(gv_clean)
    }
  }
  resolve_taxon_or_rank_vector(df)
}

format_icn_label_vec <- function(taxa_vec, abbreviate_genus = FALSE, counts_map = NULL) {
  tx_chr <- as.character(taxa_vec)
  if (length(tx_chr) == 0L) return(expression())
  out_expr <- as.expression(lapply(tx_chr, function(tx) {
    n_v <- if (!is.null(counts_map) && !is.na(tx) && tx %in% names(counts_map)) counts_map[[tx]] else NULL
    format_icn_label_expr(tx, abbreviate_genus = abbreviate_genus, count_n = n_v)
  }))
  names(out_expr) <- tx_chr
  out_expr
}

icn_discrete_labeller <- function(abbreviate_genus = FALSE, counts_map = NULL) {
  force(abbreviate_genus)
  force(counts_map)
  function(x) {
    unname(format_icn_label_vec(x, abbreviate_genus = abbreviate_genus, counts_map = counts_map))
  }
}

get_combinatorial_palette_r <- function(groups_vec, preserve_order = FALSE) {
  lvls <- unique(as.character(stats::na.omit(groups_vec)))
  lvls <- lvls[nzchar(trimws(lvls))]
  if (!isTRUE(preserve_order)) lvls <- sort(lvls)
  if (length(lvls) == 0L) return(character(0))
  n_c <- length(PACKERA_TOL_OKABE_12)
  cols <- PACKERA_TOL_OKABE_12[((seq_along(lvls) - 1L) %% n_c) + 1L]
  stats::setNames(cols, lvls)
}

get_combinatorial_hollow_shapes_r <- function(groups_vec, preserve_order = FALSE) {
  lvls <- unique(as.character(stats::na.omit(groups_vec)))
  lvls <- lvls[nzchar(trimws(lvls))]
  if (!isTRUE(preserve_order)) lvls <- sort(lvls)
  if (length(lvls) == 0L) return(integer(0))
  n_c <- length(PACKERA_TOL_OKABE_12)
  n_s <- length(PACKERA_HOLLOW_SHAPES_5)
  shps <- if (length(lvls) > 8L) {
    PACKERA_HOLLOW_SHAPES_5[(((seq_along(lvls) - 1L) %/% n_c) %% n_s) + 1L]
  } else {
    PACKERA_HOLLOW_SHAPES_5[((seq_along(lvls) - 1L) %% n_s) + 1L]
  }
  stats::setNames(shps, lvls)
}

theme_packera_pub <- function(base_size = 10.2) {
  ggplot2::theme_bw(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 0.8, hjust = 0, color = "#111827"),
      plot.subtitle = ggplot2::element_text(size = base_size - 1.1, hjust = 0, color = "#4B5563"),
      axis.title = ggplot2::element_text(size = base_size - 0.3, color = "#1F2937"),
      axis.text = ggplot2::element_text(size = base_size - 1.1, color = "#1F2937"),
      panel.border = ggplot2::element_rect(color = "#9CA3AF", fill = NA, linewidth = 0.65),
      panel.grid.major = ggplot2::element_line(color = "#E5E7EB", linewidth = 0.35, linetype = "dotted"),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F3F4F6", color = "#D1D5DB", linewidth = 0.5),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 0.8, color = "#111827"),
      legend.background = ggplot2::element_rect(fill = "white", color = "#D1D5DB", linewidth = 0.4),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1.0),
      legend.text = ggplot2::element_text(size = base_size - 1.3),
      plot.margin = ggplot2::margin(6, 8, 6, 8, unit = "pt")
    )
}

wrap_clean_label_r <- function(txt_vec, width = 28L) {
  vapply(as.character(txt_vec), function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return("")
    paste(strwrap(trimws(s), width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

add_centroid_repel_labels_gg <- function(p, cent_df, x_col, y_col, label_col, color_col = NULL, size = 2.85, allow_ggrepel = TRUE) {
  if (is.null(cent_df) || nrow(cent_df) == 0L) return(p)
  xc <- as.character(x_col)[1L]
  yc <- as.character(y_col)[1L]
  lc <- as.character(label_col)[1L]
  if (!all(c(xc, yc, lc) %in% names(cent_df))) return(p)
  xv <- suppressWarnings(as.numeric(cent_df[[xc]]))
  yv <- suppressWarnings(as.numeric(cent_df[[yc]]))
  lv <- as.character(cent_df[[lc]])
  has_col <- !is.null(color_col) && length(color_col) > 0L && nzchar(as.character(color_col)[1L]) && (as.character(color_col)[1L] %in% names(cent_df))
  cv <- if (has_col) as.character(cent_df[[as.character(color_col)[1L]]]) else rep("#111827", nrow(cent_df))
  ok <- is.finite(xv) & is.finite(yv) & !is.na(lv) & nzchar(trimws(lv))
  if (!any(ok)) return(p)
  rep_df <- data.frame(rx = xv[ok], ry = yv[ok], rlbl = lv[ok], rcol = cv[ok], stringsAsFactors = FALSE)
  is_latin <- all(grepl("^[A-Z][a-z.]*\\s+[a-z]", rep_df$rlbl))
  face_use <- if (is_latin) "italic" else "plain"
  sz_num <- as.numeric(size)[1L]
  if (isTRUE(allow_ggrepel) && requireNamespace("ggrepel", quietly = TRUE)) {
    rep_fmls <- names(formals(ggrepel::geom_label_repel))
    base_args <- list(
      data = rep_df,
      fill = grDevices::adjustcolor("white", alpha.f = 0.88),
      size = sz_num,
      fontface = face_use,
      show.legend = FALSE,
      inherit.aes = FALSE
    )
    if ("max.overlaps" %in% rep_fmls) base_args$max.overlaps <- 50
    if ("box.padding" %in% rep_fmls) base_args$box.padding <- 0.28
    if (has_col) {
      base_args$mapping <- ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl, color = .data$rcol)
    } else {
      base_args$mapping <- ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl)
      base_args$color <- "#111827"
    }
    lyr <- tryCatch(do.call(ggrepel::geom_label_repel, base_args), error = function(e) NULL)
    if (!is.null(lyr)) return(p + lyr)
  }
  if (has_col) {
    return(p + ggplot2::geom_label(
      data = rep_df,
      mapping = ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl, color = .data$rcol),
      fill = grDevices::adjustcolor("white", alpha.f = 0.86),
      size = sz_num,
      fontface = face_use,
      show.legend = FALSE,
      inherit.aes = FALSE
    ))
  }
  p + ggplot2::geom_label(
    data = rep_df,
    mapping = ggplot2::aes(x = .data$rx, y = .data$ry, label = .data$rlbl),
    fill = grDevices::adjustcolor("white", alpha.f = 0.86),
    color = "#111827",
    size = sz_num,
    fontface = face_use,
    show.legend = FALSE,
    inherit.aes = FALSE
  )
}

build_na_projected_basemap_gg <- function(lon_vec, lat_vec, target_crs = 5070L) {
  lon <- suppressWarnings(as.numeric(lon_vec))
  lat <- suppressWarnings(as.numeric(lat_vec))
  ok <- is.finite(lon) & is.finite(lat) & lat >= 10 & lat <= 75 & lon >= -170 & lon <= -50
  lon_use <- if (any(ok)) lon[ok] else c(-96.0, -75.0)
  lat_use <- if (any(ok)) lat[ok] else c(30.0, 44.0)

  pad_lon <- max(2.2, 0.14 * diff(range(lon_use)))
  pad_lat <- max(1.8, 0.14 * diff(range(lat_use)))
  xlim_ll <- c(max(-168, min(lon_use) - pad_lon), min(-52, max(lon_use) + pad_lon))
  ylim_ll <- c(max(18, min(lat_use) - pad_lat), min(72, max(lat_use) + pad_lat))

  can_sf <- requireNamespace("sf", quietly = TRUE) && requireNamespace("maps", quietly = TRUE)
  if (isTRUE(can_sf)) {
    bg_sf <- tryCatch({
      m_usa <- maps::map("state", plot = FALSE, fill = TRUE)
      sf_usa <- sf::st_as_sf(m_usa)
      sf::st_crs(sf_usa) <- 4326
      sf_usa <- sf::st_make_valid(sf_usa)
      sf::st_transform(sf_usa, target_crs)
    }, error = function(e) NULL)

    bbox_proj <- tryCatch({
      pts_bb <- sf::st_sfc(
        sf::st_point(c(xlim_ll[1], ylim_ll[1])),
        sf::st_point(c(xlim_ll[2], ylim_ll[1])),
        sf::st_point(c(xlim_ll[2], ylim_ll[2])),
        sf::st_point(c(xlim_ll[1], ylim_ll[2])),
        crs = 4326
      )
      bb <- sf::st_bbox(sf::st_transform(pts_bb, target_crs))
      bb
    }, error = function(e) NULL)

    if (!is.null(bg_sf) && !is.null(bbox_proj)) {
      # Compute scale-bar length in meters rounded to nice km
      span_x_m <- as.numeric(bbox_proj["xmax"] - bbox_proj["xmin"])
      span_y_m <- as.numeric(bbox_proj["ymax"] - bbox_proj["ymin"])
      nice_km_opts <- c(25, 50, 100, 200, 250, 500, 1000)
      bar_km <- nice_km_opts[which.min(abs(nice_km_opts * 1000 - 0.20 * span_x_m))]
      bar_m <- bar_km * 1000.0
      bx1 <- as.numeric(bbox_proj["xmax"]) - 0.05 * span_x_m
      bx0 <- bx1 - bar_m
      by0 <- as.numeric(bbox_proj["ymin"]) + 0.06 * span_y_m
      nx <- as.numeric(bbox_proj["xmax"]) - 0.06 * span_x_m
      ny0 <- as.numeric(bbox_proj["ymax"]) - 0.14 * span_y_m
      ny1 <- as.numeric(bbox_proj["ymax"]) - 0.05 * span_y_m

      p_base <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = bg_sf, fill = "#FAFAFA", color = "#CBD5E1", linewidth = 0.30) +
        ggplot2::annotate("segment", x = bx0, xend = bx1, y = by0, yend = by0, linewidth = 0.85, color = "#1F2937") +
        ggplot2::annotate("text", x = 0.5 * (bx0 + bx1), y = by0 + 0.032 * span_y_m, label = paste0(bar_km, " km"), size = 2.7, color = "#1F2937") +
        ggplot2::annotate("segment", x = nx, xend = nx, y = ny0, yend = ny1, linewidth = 0.75, color = "#1F2937",
                          arrow = grid::arrow(length = grid::unit(0.14, "cm"), type = "open")) +
        ggplot2::annotate("text", x = nx, y = ny1 + 0.028 * span_y_m, label = "N", fontface = "bold", size = 2.9, color = "#1F2937") +
        ggplot2::coord_sf(
          crs = target_crs,
          xlim = c(as.numeric(bbox_proj["xmin"]), as.numeric(bbox_proj["xmax"])),
          ylim = c(as.numeric(bbox_proj["ymin"]), as.numeric(bbox_proj["ymax"])),
          expand = FALSE
        ) +
        theme_packera_pub() +
        ggplot2::labs(x = "Longitude", y = "Latitude")
      return(list(plot = p_base, projected = TRUE, crs = target_crs, xlim = xlim_ll, ylim = ylim_ll))
    }
  }

  # Conformal fallback when sf/maps projection is unavailable
  p_fb <- ggplot2::ggplot()
  if (requireNamespace("maps", quietly = TRUE)) {
    st_df <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)
    if (!is.null(st_df) && nrow(st_df) > 0L) {
      p_fb <- p_fb + ggplot2::geom_polygon(
        data = st_df,
        ggplot2::aes(x = .data$long, y = .data$lat, group = .data$group),
        fill = "#FAFAFA", color = "#CBD5E1", linewidth = 0.30
      )
    }
  }
  cos_lat <- cos(mean(ylim_ll) * pi / 180.0)
  p_fb <- p_fb +
    ggplot2::coord_fixed(ratio = 1.0 / max(0.25, cos_lat), xlim = xlim_ll, ylim = ylim_ll, expand = FALSE) +
    theme_packera_pub() +
    ggplot2::labs(x = "Longitude (\u00b0E)", y = "Latitude (\u00b0N)")
  list(plot = p_fb, projected = FALSE, crs = 4326L, xlim = xlim_ll, ylim = ylim_ll)
}

project_lonlat_points_r <- function(df, map_ctx, lon_col = "decimalLongitude", lat_col = "decimalLatitude") {
  if (is.null(df) || nrow(df) == 0L || !all(c(lon_col, lat_col) %in% names(df))) return(data.frame())
  lon <- suppressWarnings(as.numeric(df[[lon_col]]))
  lat <- suppressWarnings(as.numeric(df[[lat_col]]))
  ok <- is.finite(lon) & is.finite(lat) & abs(lat) <= 90 & abs(lon) <= 180
  if (!any(ok)) return(data.frame())
  sub <- df[ok, , drop = FALSE]
  sub$map_x <- lon[ok]
  sub$map_y <- lat[ok]
  if (isTRUE(map_ctx$projected) && requireNamespace("sf", quietly = TRUE)) {
    xy_p <- tryCatch({
      pts_sf <- sf::st_as_sf(sub, coords = c(lon_col, lat_col), crs = 4326, remove = FALSE)
      pts_tr <- sf::st_transform(pts_sf, map_ctx$crs)
      sf::st_coordinates(pts_tr)
    }, error = function(e) NULL)
    if (!is.null(xy_p) && nrow(xy_p) == nrow(sub)) {
      sub$map_x <- as.numeric(xy_p[, 1L])
      sub$map_y <- as.numeric(xy_p[, 2L])
    }
  }
  sub
}

resolve_lm3_exemplar_manifest_r <- function(figures_dir, zone = "basal") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) return(list())
  z <- tolower(trimws(as.character(zone)))
  cands <- unique(c(
    file.path(figures_dir, "lm3_assets", sprintf("manifest_%s.json", z)),
    file.path(figures_dir, "ect_overlays", sprintf("manifest_%s.json", z)),
    file.path("outputs/figures/lm3_assets", sprintf("manifest_%s.json", z)),
    file.path("outputs/figures/ect_overlays", sprintf("manifest_%s.json", z))
  ))
  for (cp in cands) {
    if (file.exists(cp) && file.info(cp)$size > 2L) {
      res <- tryCatch(jsonlite::fromJSON(cp, simplifyVector = FALSE), error = function(e) NULL)
      if (is.list(res) && length(res) > 0L) return(res)
    }
  }
  list()
}

read_png_raster_grob_safe <- function(png_path) {
  if (is.null(png_path) || !nzchar(as.character(png_path)) || !file.exists(png_path)) return(NULL)
  if (!requireNamespace("png", quietly = TRUE)) return(NULL)
  arr <- tryCatch(png::readPNG(png_path), error = function(e) NULL)
  if (is.null(arr)) return(NULL)
  grid::rasterGrob(arr, interpolate = TRUE)
}
