#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_synthesis_plates.R
# Step 07 synthesis plates; sources step04_allometry_merger_utils.R if needed.
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

if (!exists("fit_gmm_em_pure", mode = "function")) {
  for (.s04_p in c("scripts/morphometrics/step04_allometry_merger_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "step04_allometry_merger_utils.R"))) {
    if (!is.na(.s04_p) && nzchar(.s04_p) && file.exists(.s04_p)) { source(.s04_p); break }
  }
}

.ensure_step07_fig_modules <- function() {
  if (!exists("render_step07a_integrative_synthesis_png", mode = "function")) {
    for (.f in c("scripts/figures/theme_packera_pub.R", "scripts/figures/step07_figures.R")) if (file.exists(.f)) source(.f)
  }
}

generate_step07_pub_figures <- function(triage_df, out_fig_pdf = NULL) invisible("")

render_step07_synthesis_plate_png <- function(df, out_path, cda_summary = NULL, focal_taxon = "") {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") || identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  .ensure_step07_fig_modules()
  out_png <- if (!is.null(out_path) && nzchar(out_path)) normalize_png_path(out_path) else "outputs/figures/fig07a_integrative_taxonomic_synthesis.png"
  render_step07a_integrative_synthesis_png(df = df, out_plot_path = out_png, cda_summary = cda_summary, focal_taxon = focal_taxon)
}

render_step07_complex_discovery_png <- function(W, part_res, out_path, min_affinity = 0.03) {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") || identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  .ensure_step07_fig_modules()
  out_png <- if (!is.null(out_path) && nzchar(out_path)) normalize_png_path(out_path) else "outputs/figures/fig07b_taxon_complex_delimitation_network.png"
  cdf_in <- if (is.list(part_res) && !is.null(part_res$complexes_df)) part_res$complexes_df else as.data.frame(part_res, stringsAsFactors = FALSE)
  render_step07b_complex_network_png(complexes_df = cdf_in, affinity_mat = W, out_complex_path = out_png, min_affinity = min_affinity)
}

render_focal_revision_synthesis_plate_png <- function(df, out_png, full_df = NULL, focal_taxon = "") {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") || identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) return(invisible(""))
  .ensure_step07_fig_modules()
  fdir <- if (!is.null(out_png) && nzchar(out_png)) dirname(out_png) else "outputs/figures"; dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  foc_df <- as.data.frame(df, stringsAsFactors = FALSE); if (nrow(foc_df) == 0L) return(invisible(""))
  foc_name <- if (nzchar(trimws(as.character(focal_taxon)))) trimws(as.character(focal_taxon)) else as.character(stats::na.omit(foc_df$species_standardized)[1L])
  sub_call <- if ("taxonomic_treatment_subtype" %in% names(foc_df)) as.character(foc_df$taxonomic_treatment_subtype) else rep("Confirmed_Species", nrow(foc_df))
  sub_call[is.na(sub_call) | !nzchar(trimws(sub_call)) | !(sub_call %in% names(PACKERA_TREATMENT_PALETTE))] <- "Confirmed_Species"
  foc_df$treatment_call <- factor(sub_call, levels = intersect(names(PACKERA_TREATMENT_PALETTE), unique(sub_call)))
  use_can <- ("can1" %in% names(foc_df)) && sum(is.finite(suppressWarnings(as.numeric(foc_df$can1)))) >= 3L
  mx_col <- if (use_can) "can1" else if ("PC1" %in% names(foc_df)) "PC1" else "ECT_PC1"
  my_col <- if (use_can && "can2" %in% names(foc_df) && isTRUE(stats::sd(suppressWarnings(as.numeric(foc_df$can2)), na.rm = TRUE) > 1e-6)) "can2" else if ("PC2" %in% names(foc_df)) "PC2" else if ("ECT_PC1" %in% names(foc_df)) "ECT_PC1" else mx_col
  foc_df$mx <- suppressWarnings(as.numeric(foc_df[[mx_col]])); foc_df$my <- suppressWarnings(as.numeric(foc_df[[my_col]]))
  nb_df <- data.frame()
  if (!is.null(full_df) && is.data.frame(full_df) && nrow(full_df) > nrow(foc_df) && all(c(mx_col, my_col) %in% names(full_df))) {
    nb_raw <- as.data.frame(full_df, stringsAsFactors = FALSE); nb_raw$mx <- suppressWarnings(as.numeric(nb_raw[[mx_col]])); nb_raw$my <- suppressWarnings(as.numeric(nb_raw[[my_col]]))
    nb_df <- nb_raw[as.character(nb_raw$species_standardized) != foc_name & is.finite(nb_raw$mx) & is.finite(nb_raw$my), , drop = FALSE]
  }
  foc_m <- foc_df[is.finite(foc_df$mx) & is.finite(foc_df$my), , drop = FALSE]
  p_a <- ggplot2::ggplot() + ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) + ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4)
  nb_cent <- NULL
  if (nrow(nb_df) > 0L) {
    nb_cent <- stats::aggregate(cbind(mx, my) ~ species_standardized, data = nb_df, FUN = stats::median); nb_cnt <- table(nb_df$species_standardized)
    nb_cent <- nb_cent[order(-as.integer(nb_cnt[nb_cent$species_standardized])), , drop = FALSE]; nb_cent$lbl <- abbreviate_icn_taxon(nb_cent$species_standardized, abbreviate_genus = TRUE)
    p_a <- p_a + ggplot2::geom_point(data = nb_df, ggplot2::aes(x = .data$mx, y = .data$my), shape = 1, size = 1.6, stroke = 0.55, color = "#94A3B8", alpha = 0.42)
  }
  p_a <- p_a + ggplot2::geom_point(data = foc_m, ggplot2::aes(x = .data$mx, y = .data$my, color = .data$treatment_call, shape = .data$treatment_call), size = 2.4, stroke = 0.9, alpha = 0.88) +
    ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE, labels = function(x) gsub("_", " ", x), name = "Treatment") +
    ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES, labels = function(x) gsub("_", " ", x), name = "Treatment") + theme_packera_pub(base_size = 9.4) +
    ggplot2::labs(title = bquote("A." ~ italic(.(foc_name)) ~ "vs. Complex Neighbors Morphospace"),
                  subtitle = if (nrow(nb_df) > 0L) "Grey circles = congeneric complex; colored hollow symbols = focal treatments (* = Type)" else "Colored hollow symbols = focal treatments (* = Type)",
                  x = if (mx_col == "can1") "Canonical Variate 1 (Can1)" else format_morph_trait_labels(mx_col), y = if (my_col == "can2") "Canonical Variate 2 (Can2)" else format_morph_trait_labels(my_col))
  if ("is_type_specimen" %in% names(foc_m) && any(!is.na(foc_m$is_type_specimen) & as.logical(foc_m$is_type_specimen))) p_a <- p_a + ggplot2::geom_point(data = foc_m[!is.na(foc_m$is_type_specimen) & as.logical(foc_m$is_type_specimen), , drop = FALSE], ggplot2::aes(x = .data$mx, y = .data$my), shape = 8, size = 3.8, stroke = 1.2, color = "#111827")
  if (!is.null(nb_cent) && nrow(nb_cent) > 0L) p_a <- add_centroid_repel_labels_gg(p_a, utils::head(nb_cent, 8L), "mx", "my", "lbl", size = 2.50)
  lon_col <- if ("decimalLongitude" %in% names(foc_df)) "decimalLongitude" else "longitude"; lat_col <- if ("decimalLatitude" %in% names(foc_df)) "decimalLatitude" else "latitude"
  m_ctx <- build_na_projected_basemap_gg(foc_df[[lon_col]], foc_df[[lat_col]]); pts_g <- project_lonlat_points_r(foc_df, m_ctx, lon_col = lon_col, lat_col = lat_col); p_b <- m_ctx$plot
  if (nrow(pts_g) > 0L) {
    p_b <- p_b + ggplot2::geom_point(data = pts_g, ggplot2::aes(x = .data$map_x, y = .data$map_y, color = .data$treatment_call, shape = .data$treatment_call), size = 2.3, stroke = 0.85, alpha = 0.86, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE) + ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES)
    if ("is_type_specimen" %in% names(pts_g) && any(!is.na(pts_g$is_type_specimen) & as.logical(pts_g$is_type_specimen))) p_b <- p_b + ggplot2::geom_point(data = pts_g[!is.na(pts_g$is_type_specimen) & as.logical(pts_g$is_type_specimen), , drop = FALSE], ggplot2::aes(x = .data$map_x, y = .data$map_y), shape = 8, size = 3.8, stroke = 1.2, color = "#111827")
  }
  p_b <- p_b + ggplot2::labs(title = "B. Geographic Distribution & Type Localities (EPSG:5070)", subtitle = "Projected CONUS Albers; star = nomenclatural type control")
  ed_col <- intersect(c("lithic_outcrop_index", "fine_bedrock_depth_cm", "fine_soil_ph", "soil_ph"), names(foc_df))[1L]; ph_col <- intersect(c("delta_doy_resid", "hopkins_adjusted_doy", "doy"), names(foc_df))[1L]
  p_c <- if (!is.na(ed_col) && !is.na(ph_col)) {
    foc_df$ex_v <- suppressWarnings(as.numeric(foc_df[[ed_col]])); foc_df$ph_v <- suppressWarnings(as.numeric(foc_df[[ph_col]]))
    ec <- foc_df[is.finite(foc_df$ex_v) & is.finite(foc_df$ph_v), , drop = FALSE]; pe <- ggplot2::ggplot()
    if (!is.null(full_df) && is.data.frame(full_df) && all(c(ed_col, ph_col) %in% names(full_df))) {
      nb_e <- as.data.frame(full_df, stringsAsFactors = FALSE); nb_e <- nb_e[as.character(nb_e$species_standardized) != foc_name, , drop = FALSE]
      nb_e$ex_v <- suppressWarnings(as.numeric(nb_e[[ed_col]])); nb_e$ph_v <- suppressWarnings(as.numeric(nb_e[[ph_col]])); nb_e <- nb_e[is.finite(nb_e$ex_v) & is.finite(nb_e$ph_v), , drop = FALSE]
      if (nrow(nb_e) > 0L) pe <- pe + ggplot2::geom_point(data = nb_e, ggplot2::aes(x = .data$ex_v, y = .data$ph_v), shape = 1, size = 1.6, stroke = 0.5, color = "#94A3B8", alpha = 0.35)
    }
    if (ph_col == "delta_doy_resid") pe <- pe + ggplot2::geom_hline(yintercept = c(-21, 0, 21), linetype = c("dashed", "dotted", "dashed"), color = c("#D55E00", "#9CA3AF", "#D55E00"), linewidth = 0.5)
    pe + ggplot2::geom_point(data = ec, ggplot2::aes(x = .data$ex_v, y = .data$ph_v, color = .data$treatment_call, shape = .data$treatment_call), size = 2.3, stroke = 0.85, alpha = 0.85, show.legend = FALSE) +
      ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE) + ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES) + theme_packera_pub(base_size = 9.4) +
      ggplot2::labs(title = "C. Rhizosphere Edaphic & Hopkins Phenology Envelope", subtitle = "Focal vouchers vs. complex background along edaphic & phenological axes",
                    x = if (exists("format_predictor_labels", mode = "function")) format_predictor_labels(ed_col) else ed_col, y = if (ph_col == "delta_doy_resid") "Hopkins Phenological Anomaly (Delta DOY, d)" else "Hopkins-Adjusted Anthesis DOY")
  } else ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "C. Rhizosphere Edaphic & Phenology Envelope")
  trt_tab <- as.data.frame(table(treatment_call = foc_df$treatment_call), stringsAsFactors = FALSE); trt_tab$pct <- 100.0 * trt_tab$Freq / max(1L, sum(trt_tab$Freq)); trt_tab <- trt_tab[order(trt_tab$Freq), , drop = FALSE]
  trt_tab$lbl <- factor(gsub("_", " ", trt_tab$treatment_call), levels = gsub("_", " ", trt_tab$treatment_call)); trt_tab$txt <- sprintf("%d (%.1f%%)", trt_tab$Freq, trt_tab$pct)
  pcs <- if ("probabilistic_consensus_score" %in% names(foc_df)) suppressWarnings(as.numeric(foc_df$probabilistic_consensus_score)) else numeric(0)
  pcs_med <- if (any(is.finite(pcs))) sprintf("Median consensus = %.3f", stats::median(pcs, na.rm = TRUE)) else ""
  n_basal_src <- if ("morph_tier_source" %in% names(foc_df)) sum(foc_df$morph_tier_source == "Basal_Rosette", na.rm = TRUE) else nrow(foc_df)
  p_d <- ggplot2::ggplot(trt_tab, ggplot2::aes(x = .data$Freq, y = .data$lbl, fill = .data$treatment_call)) + ggplot2::geom_col(width = 0.62, color = "white", linewidth = 0.35, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = .data$txt), hjust = -0.10, size = 2.8, color = "#111827") + ggplot2::scale_fill_manual(values = PACKERA_TREATMENT_PALETTE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25))) + theme_packera_pub(base_size = 9.4) +
    ggplot2::labs(title = bquote("D." ~ italic(.(foc_name)) ~ .(sprintf("Treatment Breakdown (N = %d)", nrow(foc_df)))),
                  subtitle = paste(c(sprintf("Basal_Rosette = %d / %d", n_basal_src, nrow(foc_df)), if (nzchar(pcs_med)) pcs_med else NULL), collapse = " | "), x = "Voucher Count", y = NULL)
  plate_f <- ((p_a | p_b) / (p_c | p_d)) + patchwork::plot_layout(heights = c(1.12, 1.0)) +
    patchwork::plot_annotation(title = bquote(bold("Integrative Taxonomic Revision Plate:" ~ italic(.(foc_name)))), theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5)))
  save_packera_png_300(plate_f, normalize_png_path(out_png), width = 14.4, height = 10.6)
}
