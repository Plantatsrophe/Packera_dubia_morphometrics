# ==============================================================================
# Module: scripts/figures/step06_tier_figures.R
# Project: Plantmatics | Affiliation: UNC-Chapel Hill Herbarium (NCU)
# Adaptive 300-DPI PNG-only 3-plate figure suite for Step 06 Tiers 2-4.
# ==============================================================================

if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

if (!exists("save_packera_png_300", mode = "function")) {
  for (.p in c("scripts/figures/theme_packera_pub.R", "../figures/theme_packera_pub.R")) {
    if (file.exists(.p)) { source(.p); break }
  }
}

.tier_empty <- function(title_txt, msg_txt, sub_txt = NULL) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg_txt, size = 3.2, color = "#334155") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_packera_pub(base_size = 9.2) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) +
    ggplot2::labs(title = title_txt, subtitle = sub_txt, x = NULL, y = NULL)
}

.tier_tip_expr <- function(raw_lbl, focal_tx = "") {
  cl <- gsub("'", "", trimws(gsub("_", " ", as.character(raw_lbl))))
  is_f <- nzchar(focal_tx) && tolower(cl) == tolower(trimws(gsub("_", " ", focal_tx)))
  tag <- if (is_f) " [FOCAL]" else ""
  if (grepl("^Packera[[:space:]]+", cl)) {
    sh <- sub("^Packera[[:space:]]+", "P. ", cl)
    sprintf("paste(italic('%s'), '%s')", sh, tag)
  } else sprintf("'%s%s'", cl, tag)
}

.tier_phylo_gg <- function(tr, title_txt, sub_txt = NULL, xlab_txt = "Substitutions / Site", focal_tx = "", show_support = TRUE) {
  if (is.null(tr) || !inherits(tr, "phylo") || length(tr$tip.label) < 3L || !requireNamespace("ape", quietly = TRUE)) {
    return(.tier_empty(title_txt, "Fewer than 3 non-rogue tips staged for phylogeny", sub_txt))
  }
  tr_w <- if (exists(".midpoint_root_phylo_safe", mode = "function")) .midpoint_root_phylo_safe(tr) else tr
  tr_u <- ape::ladderize(tr_w, right = FALSE)
  nt <- length(tr_u$tip.label); tot <- nt + tr_u$Nnode; ord <- ape::reorder.phylo(tr_u, order = "cladewise")
  dx <- if (!is.null(tr_u$edge.length)) ape::node.depth.edgelength(tr_u) else {
    d <- numeric(tot); for (ei in seq_len(nrow(ord$edge))) d[ord$edge[ei, 2L]] <- d[ord$edge[ei, 1L]] + 1.0; d
  }
  yp <- numeric(tot); yp[ord$edge[ord$edge[, 2L] <= nt, 2L]] <- seq_len(nt)
  po <- ape::reorder.phylo(tr_u, order = "postorder"); pp <- unique(po$edge[, 1L])
  for (nd in pp) yp[nd] <- mean(yp[po$edge[po$edge[, 1L] == nd, 2L]])
  hs <- data.frame(x = dx[ord$edge[, 1L]], xend = dx[ord$edge[, 2L]], y = yp[ord$edge[, 2L]], yend = yp[ord$edge[, 2L]])
  vs <- do.call(rbind, lapply(pp, function(nd) { ch <- po$edge[po$edge[, 1L] == nd, 2L]; data.frame(x = dx[nd], xend = dx[nd], y = min(yp[ch]), yend = max(yp[ch])) }))
  tx_c <- trimws(gsub("_", " ", tr_u$tip.label))
  tdf <- data.frame(tip = tr_u$tip.label, taxon = tx_c, x = dx[seq_len(nt)], y = yp[seq_len(nt)],
                    expr = vapply(tx_c, .tier_tip_expr, character(1), focal_tx = focal_tx), stringsAsFactors = FALSE)
  mx <- max(c(tdf$x, 1e-4), na.rm = TRUE)
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = hs, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend), color = "#334155", linewidth = 0.45) +
    ggplot2::geom_segment(data = vs, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend), color = "#334155", linewidth = 0.45) +
    ggplot2::geom_point(data = tdf, ggplot2::aes(x = .data$x, y = .data$y, color = .data$taxon, shape = .data$taxon), size = 2.6, fill = NA, stroke = 0.95) +
    ggplot2::geom_text(data = tdf, ggplot2::aes(x = .data$x + 0.025 * mx, y = .data$y, label = .data$expr), parse = TRUE, hjust = 0, size = 2.45, color = "#0F172A") +
    ggplot2::scale_color_manual(values = get_combinatorial_palette_r(tdf$taxon)) +
    ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(tdf$taxon)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.52))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title_txt, subtitle = sub_txt, x = xlab_txt, y = NULL) +
    theme_packera_pub(base_size = 9.0) +
    ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
  if (isTRUE(show_support) && !is.null(tr_u$node.label) && length(tr_u$node.label) > 0L) {
    nd_lbl <- trimws(as.character(tr_u$node.label)); nd_ids <- nt + seq_along(nd_lbl)
    ok_nd <- nzchar(nd_lbl) & nd_lbl != "NA"
    if (any(ok_nd)) {
      ndf <- data.frame(x = dx[nd_ids[ok_nd]], y = yp[nd_ids[ok_nd]], lbl = nd_lbl[ok_nd], stringsAsFactors = FALSE)
      p <- p + ggplot2::geom_text(data = ndf, ggplot2::aes(x = .data$x, y = .data$y, label = .data$lbl), hjust = 1.15, vjust = -0.35, size = 2.2, color = "#475569")
    }
  }
  p
}

.tier_pcoa_gg <- function(pco_df, title_txt, sub_txt = NULL) {
  if (!is.data.frame(pco_df) || nrow(pco_df) < 3L) return(.tier_empty(title_txt, "Insufficient distance pairs for ordination", sub_txt))
  lab_fn <- if (exists("icn_discrete_labeller", mode = "function")) icn_discrete_labeller(abbreviate_genus = TRUE) else ggplot2::waiver()
  e1 <- if ("eig1" %in% names(pco_df)) sprintf("PCo1 (%.1f%%)", pco_df$eig1[1L]) else "PCo1"
  e2 <- if ("eig2" %in% names(pco_df)) sprintf("PCo2 (%.1f%%)", pco_df$eig2[1L]) else "PCo2"
  ggplot2::ggplot(pco_df, ggplot2::aes(x = .data$PCo1, y = .data$PCo2, color = .data$taxon, shape = .data$taxon)) +
    ggplot2::geom_point(size = 3.0, fill = NA, stroke = 1.05) +
    ggplot2::scale_color_manual(values = get_combinatorial_palette_r(pco_df$taxon), labels = lab_fn, name = "Taxon") +
    ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(pco_df$taxon), labels = lab_fn, name = "Taxon") +
    ggplot2::labs(title = title_txt, subtitle = sub_txt, x = e1, y = e2) +
    theme_packera_pub(base_size = 9.0) +
    ggplot2::theme(legend.position = "right")
}

.tier_gdi_gg <- function(gdi_df, title_txt, sub_txt) {
  if (!is.data.frame(gdi_df) || nrow(gdi_df) == 0L || !any(is.finite(gdi_df$coalescent_gdi))) {
    return(.tier_empty(title_txt, "Coalescent gdi below multi-locus sample threshold", sub_txt))
  }
  gd <- gdi_df[is.finite(gdi_df$coalescent_gdi), , drop = FALSE]
  gd$taxon <- factor(gd$taxon, levels = gd$taxon[order(gd$coalescent_gdi)])
  cls <- sort(unique(gd$gdi_classification))
  ggplot2::ggplot(gd, ggplot2::aes(x = .data$coalescent_gdi, y = .data$taxon, color = .data$gdi_classification, shape = .data$gdi_classification)) +
    ggplot2::geom_vline(xintercept = 0.20, linetype = "dashed", color = "#DC2626") +
    ggplot2::geom_vline(xintercept = 0.70, linetype = "dashed", color = "#059669") +
    ggplot2::geom_segment(ggplot2::aes(x = .data$gdi_ci_low, xend = .data$gdi_ci_high, y = .data$taxon, yend = .data$taxon), linewidth = 0.65) +
    ggplot2::geom_point(size = 3.0, stroke = 1.1, fill = NA) +
    ggplot2::scale_color_manual(values = get_combinatorial_palette_r(cls), name = "MSC Gate") +
    ggplot2::scale_shape_manual(values = stats::setNames(rep(c(1L, 2L, 5L, 0L), length.out = length(cls)), cls), name = "MSC Gate") +
    ggplot2::scale_x_continuous(limits = c(0, 1.02)) +
    ggplot2::labs(title = title_txt, subtitle = sub_txt, x = "Genealogical Divergence Index (gdi +/- 95% CI)", y = NULL) +
    theme_packera_pub(base_size = 9.0) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic", size = 8.2))
}

.tier_cytology_gg <- function(cyto_rec, cyto_sum, active_taxa, title_txt) {
  df <- if (is.data.frame(cyto_rec) && nrow(cyto_rec) > 0L) cyto_rec else data.frame()
  if (nrow(df) > 0L && "sporophytic_2n" %in% names(df)) {
    df$val2n <- suppressWarnings(as.numeric(df$sporophytic_2n))
    df <- df[is.finite(df$val2n), , drop = FALSE]
  } else df <- data.frame()
  if (nrow(df) == 0L && is.data.frame(cyto_sum) && nrow(cyto_sum) > 0L) {
    v2 <- suppressWarnings(as.numeric(cyto_sum$modal_sporophytic_2n))
    ok <- is.finite(v2)
    if (any(ok)) df <- data.frame(species_standardized = cyto_sum$species_standardized[ok], val2n = v2[ok], ploidy_level = paste0(cyto_sum$modal_ploidy_level[ok], "x"), stringsAsFactors = FALSE)
  }
  if (nrow(df) == 0L) return(.tier_empty(title_txt, "No CCDB/ICNA chromosome counts staged in SQLite"))
  if (length(active_taxa) > 0L && any(df$species_standardized %in% active_taxa)) {
    df <- df[df$species_standardized %in% active_taxa, , drop = FALSE]
  }
  lab_fn <- if (exists("icn_discrete_labeller", mode = "function")) icn_discrete_labeller(abbreviate_genus = TRUE) else ggplot2::waiver()
  df$pl_lbl <- if ("ploidy_level" %in% names(df)) as.character(df$ploidy_level) else paste0(round(df$val2n / 23), "x")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$species_standardized, y = .data$val2n, color = .data$pl_lbl, shape = .data$pl_lbl)) +
    ggplot2::geom_point(size = 2.8, fill = NA, stroke = 1.0, position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 42L)) +
    ggplot2::scale_color_manual(values = get_combinatorial_palette_r(unique(df$pl_lbl)), name = "Ploidy") +
    ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(unique(df$pl_lbl)), name = "Ploidy") +
    ggplot2::scale_x_discrete(labels = lab_fn) +
    ggplot2::labs(title = title_txt, subtitle = sprintf("CCDB/ICNA SQLite records (n = %d)", nrow(df)), x = NULL, y = "Sporophytic Count (2n)") +
    theme_packera_pub(base_size = 9.0) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 7.8))
}

.tier_geo_map_gg <- function(geo_df, title_txt, sub_txt = NULL, color_col = "taxon") {
  if (!is.data.frame(geo_df) || nrow(geo_df) == 0L || !all(c("decimalLongitude", "decimalLatitude") %in% names(geo_df))) {
    return(.tier_empty(title_txt, "Zero georeferenced coordinates staged", sub_txt))
  }
  geo_df$lon <- suppressWarnings(as.numeric(geo_df$decimalLongitude)); geo_df$lat <- suppressWarnings(as.numeric(geo_df$decimalLatitude))
  geo_df <- geo_df[is.finite(geo_df$lon) & is.finite(geo_df$lat), , drop = FALSE]
  if (nrow(geo_df) == 0L) return(.tier_empty(title_txt, "Zero finite lon/lat points", sub_txt))
  grp <- if (color_col %in% names(geo_df)) as.character(geo_df[[color_col]]) else rep("Voucher", nrow(geo_df))
  grp[is.na(grp) | !nzchar(grp)] <- "Unassigned"; geo_df$grp_col <- grp
  pal <- get_combinatorial_palette_r(unique(grp)); sh <- get_combinatorial_hollow_shapes_r(unique(grp))
  if (exists("build_na_projected_basemap_gg", mode = "function") && exists("project_lonlat_points_r", mode = "function")) {
    ctx <- build_na_projected_basemap_gg(geo_df$lon, geo_df$lat, target_crs = 5070L)
    pts <- project_lonlat_points_r(geo_df, ctx, lon_col = "lon", lat_col = "lat")
    return(ctx$plot + ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$map_x, y = .data$map_y, color = .data$grp_col, shape = .data$grp_col), size = 2.6, fill = NA, stroke = 1.0, inherit.aes = FALSE) +
      ggplot2::scale_color_manual(values = pal, name = color_col) + ggplot2::scale_shape_manual(values = sh, name = color_col) +
      ggplot2::labs(title = title_txt, subtitle = sub_txt) + theme_packera_pub(base_size = 9.0))
  }
  ggplot2::ggplot(geo_df, ggplot2::aes(x = .data$lon, y = .data$lat, color = .data$grp_col, shape = .data$grp_col)) +
    ggplot2::geom_point(size = 2.6, fill = NA, stroke = 1.0) +
    ggplot2::scale_color_manual(values = pal, name = color_col) + ggplot2::scale_shape_manual(values = sh, name = color_col) +
    ggplot2::labs(title = title_txt, subtitle = sub_txt, x = "Longitude", y = "Latitude") + theme_packera_pub(base_size = 9.0)
}

.tier_phytop_gg <- function(conc_df, title_txt) {
  if (!is.data.frame(conc_df) || !all(c("public_phytop_ils_index", "public_phytop_ih_index") %in% names(conc_df))) {
    return(.tier_empty(title_txt, "PhyTop ILS/IH conflict tables not staged"))
  }
  ils <- suppressWarnings(as.numeric(conc_df$public_phytop_ils_index)); ih <- suppressWarnings(as.numeric(conc_df$public_phytop_ih_index))
  ok <- is.finite(ils) & is.finite(ih)
  if (!any(ok)) return(.tier_empty(title_txt, "No finite PhyTop ILS/IH indices matched"))
  df <- data.frame(taxon = conc_df$taxon[ok], ils = ils[ok], ih = ih[ok], stringsAsFactors = FALSE)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$ils, y = .data$ih, color = .data$taxon, shape = .data$taxon)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#94A3B8") +
    ggplot2::geom_point(size = 3.0, fill = NA, stroke = 1.05) +
    ggplot2::scale_color_manual(values = get_combinatorial_palette_r(df$taxon), name = "Taxon") +
    ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(df$taxon), name = "Taxon") +
    ggplot2::labs(title = title_txt, x = "PhyTop ILS Index", y = "PhyTop Introgression Index (IH)") +
    theme_packera_pub(base_size = 9.0)
}

.render_tier2_plates <- function(fig_a, fig_b, fig_c, p, db_path, run_id) {
  foc <- p$focal_taxon %||% ""
  p6a_a <- .tier_pcoa_gg(p$pco_df, "A. Hyb-Seq Supermatrix Ordination (PCoA)")
  p6a_b <- .tier_phylo_gg(p$tree_primary, "B. Pruned ASTRAL Nuclear Backbone (LPP)", focal_tx = foc)
  rf_sub <- if (is.finite(p$rf_norm)) sprintf("Nuclear-Plastome RF_norm = %.3f", p$rf_norm) else "Plastome reference topology"
  p6a_c <- .tier_phylo_gg(p$tree_secondary, "C. Mirrored 111-Taxon Plastome Topology", sub_txt = rf_sub, focal_tx = foc)
  tile_df <- if (is.data.frame(p$concordance_df) && nrow(p$concordance_df) > 0L) {
    cdf <- p$concordance_df
    rbind(
      data.frame(taxon = cdf$taxon, stream = "MO_Sister", call = cdf$mo_nearest_sister, stringsAsFactors = FALSE),
      data.frame(taxon = cdf$taxon, stream = "ASTRALPro_Sister", call = cdf$astralpro_nearest_sister, stringsAsFactors = FALSE)
    )
  } else data.frame()
  p6a_d <- if (nrow(tile_df) > 0L) {
    tile_df$call[is.na(tile_df$call) | !nzchar(tile_df$call)] <- "Unlinked"
    ggplot2::ggplot(tile_df, ggplot2::aes(x = .data$stream, y = .data$taxon, fill = .data$call)) +
      ggplot2::geom_tile(color = "#FFFFFF") +
      ggplot2::scale_fill_manual(values = c(get_combinatorial_palette_r(setdiff(unique(tile_df$call), "Unlinked")), Unlinked = "#E2E8F0"), name = "Sister") +
      ggplot2::labs(title = "D. Multi-Tree Sister & Cross-Modal Tile", x = "Stream", y = NULL) + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("D. Multi-Tree Sister & Cross-Modal Tile", "No concordance rows")
  save_packera_png_300((p6a_a | p6a_b) / (p6a_c | p6a_d) + patchwork::plot_annotation(title = "Step 06 [Tier 2 - Public Hyb-Seq]: Phylogenomic Concordance & Ordination"), fig_a, width = 15.0, height = 11.0, dpi = 300L)

  p6b_a <- .tier_phytop_gg(p$concordance_df, "A. PhyTop Conflict (ILS vs. Introgression)")
  p6b_b <- if (is.data.frame(p$svdq_df) && nrow(p$svdq_df) > 0L) {
    sd <- utils::head(p$svdq_df, 18L)
    ggplot2::ggplot(sd, ggplot2::aes(x = stats::reorder(.data$quartet_id, .data$quartet_support), y = .data$quartet_support)) +
      ggplot2::geom_point(shape = 1L, size = 2.8, stroke = 1.0, fill = NA, color = "#0072B2") +
      ggplot2::coord_flip() + ggplot2::labs(title = "B. SVDquartets Topology Support", x = "Quartet", y = "Quartet Support") + theme_packera_pub(base_size = 8.5)
  } else .tier_empty("B. SVDquartets Topology Support", "Fewer than 4 taxa for quartet scoring")
  p6b_c <- .tier_gdi_gg(p$gdi_df, "C. Hyb-Seq Coalescent Delimitation (gdi)", "Oversplit < 0.20 | Distinct Species > 0.70")
  p6b_d <- .tier_phylo_gg(p$reltime_tree, "D. Pruned RelTime Divergence Chronogram", xlab_txt = "Relative Divergence Time (Ma)", focal_tx = foc, show_support = FALSE)
  save_packera_png_300((p6b_a | p6b_b) / (p6b_c | p6b_d) + patchwork::plot_annotation(title = "Step 06 [Tier 2 - Public Hyb-Seq]: Reticulation, gdi & Divergence Timing"), fig_b, width = 15.0, height = 10.5, dpi = 300L)

  p6c_a <- if (is.data.frame(p$completeness_df) && nrow(p$completeness_df) > 0L) {
    ggplot2::ggplot(p$completeness_df, ggplot2::aes(x = .data$site_completeness, y = .data$informative_density, color = .data$taxon, shape = .data$taxon)) +
      ggplot2::geom_point(size = 2.8, fill = NA, stroke = 1.0) + ggplot2::scale_color_manual(values = get_combinatorial_palette_r(p$completeness_df$taxon)) +
      ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(p$completeness_df$taxon)) +
      ggplot2::labs(title = "A. Supermatrix Completeness vs. Informative Density", x = "Site Completeness", y = "Informative Density") + theme_packera_pub(base_size = 8.8) + ggplot2::theme(legend.position = "none")
  } else .tier_empty("A. Supermatrix Completeness", "Tree-derived Hyb-Seq supermatrix active")
  p6c_b <- if (is.data.frame(p$concordance_df) && nrow(p$concordance_df) > 0L) {
    cdf <- p$concordance_df
    m_df <- rbind(
      data.frame(taxon = cdf$taxon, tree = "MO", present = as.numeric(cdf$in_mo_backbone), stringsAsFactors = FALSE),
      data.frame(taxon = cdf$taxon, tree = "ASTRALPro", present = as.numeric(cdf$in_astralpro_backbone), stringsAsFactors = FALSE)
    )
    ggplot2::ggplot(m_df, ggplot2::aes(x = .data$tree, y = .data$taxon, fill = factor(.data$present))) +
      ggplot2::geom_tile(color = "#FFFFFF") + ggplot2::scale_fill_manual(values = c("0" = "#E2E8F0", "1" = "#009E73"), name = "Present") +
      ggplot2::labs(title = "B. Nuclear-Plastome Cophylogeny Coverage", x = "Backbone", y = NULL) + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("B. Nuclear-Plastome Cophylogeny Coverage", "No concordance table")
  p6c_c <- if (is.data.frame(p$concordance_df) && any(is.finite(p$concordance_df$mo_astralpro_jaccard))) {
    jc <- p$concordance_df[is.finite(p$concordance_df$mo_astralpro_jaccard), , drop = FALSE]
    ggplot2::ggplot(jc, ggplot2::aes(x = stats::reorder(.data$taxon, .data$mo_astralpro_jaccard), y = .data$mo_astralpro_jaccard)) +
      ggplot2::geom_col(fill = "#0072B2") + ggplot2::coord_flip() + ggplot2::labs(title = "C. Multi-ASTRAL Neighbor Jaccard Stability", x = NULL, y = "Jaccard Similarity") + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("C. Multi-ASTRAL Neighbor Jaccard Stability", "Single nuclear backbone tree")
  p6c_d <- .tier_cytology_gg(p$cyto_records_df, p$cyto_summary_df, p$active_taxa, "D. CCDB/ICNA Sporophytic 2n & Ploidy Audit")
  save_packera_png_300((p6c_a | p6c_b) / (p6c_c | p6c_d) + patchwork::plot_annotation(title = "Step 06 [Tier 2 - Public Hyb-Seq]: Completeness, Cophylogeny & Cytotype Audit"), fig_c, width = 15.0, height = 10.5, dpi = 300L)
}

.render_tier3_plates <- function(fig_a, fig_b, fig_c, p, db_path, run_id) {
  foc <- p$focal_taxon %||% ""
  p6a_a <- .tier_pcoa_gg(p$pco_df, "A. Multi-Locus Pairwise Divergence PCoA")
  p6a_b <- .tier_phylo_gg(p$tree_primary, "B. Rooted nrITS/ETS BIONJ Gene Tree (B=200)", focal_tx = foc)
  p6a_c <- .tier_phylo_gg(p$tree_secondary, "C. Plastid / waxy BIONJ Gene Tree (B=200)", focal_tx = foc)
  p6a_d <- if (is.data.frame(p$monophyly_df) && nrow(p$monophyly_df) > 0L) {
    ggplot2::ggplot(p$monophyly_df, ggplot2::aes(x = .data$locus, y = .data$taxon, fill = .data$monophyly_status)) +
      ggplot2::geom_tile(color = "#FFFFFF") + ggplot2::scale_fill_manual(values = get_combinatorial_palette_r(unique(p$monophyly_df$monophyly_status)), name = "Monophyly") +
      ggplot2::labs(title = "D. Locus Coverage & Monophyly Status", x = "Locus", y = NULL) + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("D. Locus Coverage & Monophyly Status", "No locus coverage table")
  save_packera_png_300((p6a_a | p6a_b) / (p6a_c | p6a_d) + patchwork::plot_annotation(title = "Step 06 [Tier 3 - Public Multi-Locus]: Gene Trees & Divergence PCoA"), fig_a, width = 15.0, height = 11.0, dpi = 300L)

  p6b_a <- if (is.data.frame(p$barcode_df) && nrow(p$barcode_df) > 0L) {
    ggplot2::ggplot(p$barcode_df, ggplot2::aes(x = .data$intra_pi, y = .data$min_inter_dxy, color = .data$locus, shape = .data$locus)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#94A3B8") +
      ggplot2::geom_point(size = 3.0, fill = NA, stroke = 1.05) + ggplot2::scale_color_manual(values = get_combinatorial_palette_r(unique(p$barcode_df$locus))) +
      ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(unique(p$barcode_df$locus))) +
      ggplot2::labs(title = "A. Barcode-Gap (Max Intra pi vs. Min Inter d_XY)", x = "Max Intraspecific pi", y = "Min Interspecific d_XY") + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("A. Barcode-Gap Plot", "Insufficient pairwise distances")
  p6b_b <- if (is.data.frame(p$dxy_df) && nrow(p$dxy_df) > 0L) {
    ggplot2::ggplot(p$dxy_df, ggplot2::aes(x = .data$taxon_a, y = .data$taxon_b, fill = .data$mean_dxy)) +
      ggplot2::geom_tile(color = "#FFFFFF") + ggplot2::scale_fill_gradient(low = "#F8FAFC", high = "#B91C1C", name = "d_XY") +
      ggplot2::labs(title = "B. Pairwise Taxon Mean Divergence (d_XY)", x = NULL, y = NULL) + theme_packera_pub(base_size = 8.4) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, face = "italic"), axis.text.y = ggplot2::element_text(face = "italic"))
  } else .tier_empty("B. Pairwise Taxon Mean d_XY", "Insufficient taxon pairs")
  p6b_c <- if (is.data.frame(p$barcode_df) && nrow(p$barcode_df) > 0L) {
    ggplot2::ggplot(p$barcode_df, ggplot2::aes(x = .data$locus, y = .data$exclusivity_rate, fill = .data$taxon)) +
      ggplot2::geom_col(position = "dodge") + ggplot2::scale_fill_manual(values = get_combinatorial_palette_r(unique(p$barcode_df$taxon)), name = "Taxon") +
      ggplot2::labs(title = "C. Per-Locus Sequence Exclusivity", x = "Locus", y = "Exclusivity") + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("C. Per-Locus Exclusivity", "No barcode summary")
  p6b_d <- .tier_geo_map_gg(p$geo_df, "D. Herbarium & Sequenced Cohort Geography")
  save_packera_png_300((p6b_a | p6b_b) / (p6b_c | p6b_d) + patchwork::plot_annotation(title = "Step 06 [Tier 3 - Public Multi-Locus]: Barcode Gap, d_XY & Geographic Scope"), fig_b, width = 15.0, height = 10.5, dpi = 300L)

  p6c_a <- if (is.data.frame(p$placement_df) && nrow(p$placement_df) > 0L) {
    ggplot2::ggplot(p$placement_df, ggplot2::aes(x = .data$n_shared_11mers, y = .data$kmer_containment, color = .data$locus, shape = .data$locus)) +
      ggplot2::geom_point(size = 2.8, fill = NA, stroke = 1.0) + ggplot2::scale_color_manual(values = get_combinatorial_palette_r(unique(p$placement_df$locus))) +
      ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(unique(p$placement_df$locus))) +
      ggplot2::labs(title = "A. 11-mer Containment QC & Rogue Audit", x = "Shared 11-mers", y = "11-mer Containment") + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("A. 11-mer Containment QC", "No placement rows")
  p6c_b <- if (is.data.frame(p$monophyly_df) && nrow(p$monophyly_df) > 0L) {
    ggplot2::ggplot(p$monophyly_df, ggplot2::aes(x = .data$locus, y = .data$taxon, fill = .data$monophyly_status)) +
      ggplot2::geom_tile(color = "#FFFFFF") + ggplot2::scale_fill_manual(values = get_combinatorial_palette_r(unique(p$monophyly_df$monophyly_status)), name = "Clade") +
      ggplot2::labs(title = "B. Multi-Locus Clade Concordance", x = "Locus", y = NULL) + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("B. Multi-Locus Clade Concordance", "Single locus only")
  p6c_c <- if (is.data.frame(p$model_sel_df) && nrow(p$model_sel_df) > 0L && "bic_weight" %in% names(p$model_sel_df)) {
    ms <- p$model_sel_df
    ms$mod_lbl <- if ("model_name" %in% names(ms)) as.character(ms$model_name) else if ("model" %in% names(ms)) as.character(ms$model) else "JC69"
    ms$aln_lbl <- if ("alignment_name" %in% names(ms)) as.character(ms$alignment_name) else if ("alignment" %in% names(ms)) as.character(ms$alignment) else "Locus"
    ggplot2::ggplot(ms, ggplot2::aes(x = .data$mod_lbl, y = as.numeric(.data$bic_weight), fill = .data$aln_lbl)) +
      ggplot2::geom_col(position = "dodge") + ggplot2::scale_fill_manual(values = get_combinatorial_palette_r(unique(ms$aln_lbl)), name = "Locus") +
      ggplot2::labs(title = "C. BIC Nucleotide Model Weights", x = "Model", y = "BIC Weight") + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("C. BIC Nucleotide Model Weights", "Model selection table unavailable")
  p6c_d <- .tier_cytology_gg(p$cyto_records_df, p$cyto_summary_df, p$active_taxa, "D. CCDB/ICNA Sporophytic 2n & Ploidy Audit")
  save_packera_png_300((p6c_a | p6c_b) / (p6c_c | p6c_d) + patchwork::plot_annotation(title = "Step 06 [Tier 3 - Public Multi-Locus]: QC, Model Weights & Cytotype Audit"), fig_c, width = 15.0, height = 10.5, dpi = 300L)
}

.render_tier4_plates <- function(fig_a, fig_b, fig_c, p, db_path, run_id) {
  foc <- p$focal_taxon %||% ""
  p6a_a <- .tier_phylo_gg(p$tree_primary, sprintf("A. Neighbor Nuclear Backbone (Unsequenced: %s)", foc), focal_tx = foc)
  p6a_b <- .tier_phylo_gg(p$tree_secondary, "B. Neighbor Plastome / Reference Topology", focal_tx = foc)
  p6a_c <- if (is.data.frame(p$gap_recs_df) && nrow(p$gap_recs_df) > 0L) {
    top_g <- utils::head(p$gap_recs_df, 12L)
    ggplot2::ggplot(top_g, ggplot2::aes(x = stats::reorder(.data$catalogNumber, .data$priority_score), y = .data$priority_score, fill = .data$determiner_tier)) +
      ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::scale_fill_manual(values = get_combinatorial_palette_r(unique(top_g$determiner_tier)), name = "Tier") +
      ggplot2::labs(title = "C. Top Ranked DNA Extraction Candidates (>=1995)", x = "Voucher", y = "Priority Score") + theme_packera_pub(base_size = 8.4)
  } else .tier_empty("C. Top DNA Extraction Candidates", "No candidate vouchers >= 1995")
  p6a_d <- if (is.data.frame(p$coverage_matrix_df) && nrow(p$coverage_matrix_df) > 0L) {
    ggplot2::ggplot(p$coverage_matrix_df, ggplot2::aes(x = .data$modality, y = .data$taxon, fill = .data$present_frac)) +
      ggplot2::geom_tile(color = "#FFFFFF") + ggplot2::scale_fill_gradient(low = "#F1F5F9", high = "#0072B2", limits = c(0, 1), name = "Coverage") +
      ggplot2::labs(title = "D. Complex Modality Coverage Matrix", x = "Modality", y = NULL) + theme_packera_pub(base_size = 8.4)
  } else .tier_empty("D. Complex Modality Coverage", "No coverage rows")
  save_packera_png_300((p6a_a | p6a_b) / (p6a_c | p6a_d) + patchwork::plot_annotation(title = "Step 06 [Tier 4 - Cytology & Sampling Gap]: Extraction Priorities & Backbone Context"), fig_a, width = 15.0, height = 11.0, dpi = 300L)

  p6b_a <- .tier_phytop_gg(p$concordance_df, "A. Neighbor-Complex PhyTop Reticulation Context")
  p6b_b <- if (is.data.frame(p$morph_eco_dist_df) && nrow(p$morph_eco_dist_df) > 0L) {
    ggplot2::ggplot(p$morph_eco_dist_df, ggplot2::aes(x = .data$morph_dist, y = .data$eco_dist, color = .data$taxon, shape = .data$taxon)) +
      ggplot2::geom_point(size = 3.0, fill = NA, stroke = 1.05) + ggplot2::scale_color_manual(values = get_combinatorial_palette_r(p$morph_eco_dist_df$taxon)) +
      ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(p$morph_eco_dist_df$taxon)) +
      ggplot2::labs(title = "B. Step 04/05 Morphological vs. Ecological Affinity", x = "Morph Distance", y = "Edaphic Distance") + theme_packera_pub(base_size = 8.8)
  } else .tier_empty("B. Morphological vs. Ecological Affinity", "Step 04/05 tables unavailable")
  p6b_c <- .tier_geo_map_gg(p$gap_recs_df, "C. Recommended DNA-Sampling Vouchers Map", color_col = "determiner_tier")
  p6b_d <- if (is.data.frame(p$gap_recs_df) && nrow(p$gap_recs_df) > 0L && "institutionCode" %in% names(p$gap_recs_df)) {
    ggplot2::ggplot(p$gap_recs_df, ggplot2::aes(x = .data$institutionCode, fill = .data$determiner_tier)) +
      ggplot2::geom_bar() + ggplot2::scale_fill_manual(values = get_combinatorial_palette_r(unique(p$gap_recs_df$determiner_tier)), name = "Tier") +
      ggplot2::labs(title = "D. Herbarium Sampling Coverage", x = "Herbarium", y = "Candidate Sheets") + theme_packera_pub(base_size = 8.6)
  } else .tier_empty("D. Herbarium Sampling Coverage", "No candidate vouchers")
  save_packera_png_300((p6b_a | p6b_b) / (p6b_c | p6b_d) + patchwork::plot_annotation(title = "Step 06 [Tier 4 - Cytology & Sampling Gap]: Phenotypic Affinity & Sampling Map"), fig_b, width = 15.0, height = 10.5, dpi = 300L)

  p6c_a <- .tier_cytology_gg(p$cyto_records_df, p$cyto_summary_df, p$active_taxa, "A. Complex Sporophytic 2n Distribution")
  p6c_b <- if (is.data.frame(p$cyto_records_df) && nrow(p$cyto_records_df) > 0L && "cytology_linkage_tier" %in% names(p$cyto_records_df)) {
    ggplot2::ggplot(p$cyto_records_df, ggplot2::aes(x = .data$cytology_linkage_tier, fill = .data$species_standardized)) +
      ggplot2::geom_bar() + ggplot2::scale_fill_manual(values = get_combinatorial_palette_r(unique(p$cyto_records_df$species_standardized)), name = "Taxon") +
      ggplot2::labs(title = "B. C1/C2 Cytovoucher Linkage Summary", x = "Linkage Tier", y = "Records") + theme_packera_pub(base_size = 8.4)
  } else .tier_empty("B. C1/C2 Cytovoucher Linkage Summary", "No cytovoucher linkage records")
  p6c_c <- .tier_geo_map_gg(p$cyto_records_df, "C. Georeferenced Cytotype Map (2x-8x)", color_col = "ploidy_level")
  p6c_d <- if (is.data.frame(p$cyto_summary_df) && nrow(p$cyto_summary_df) >= 2L) {
    cs <- p$cyto_summary_df[p$cyto_summary_df$species_standardized %in% p$active_taxa, , drop = FALSE]
    if (nrow(cs) < 2L) cs <- utils::head(p$cyto_summary_df, 8L)
    mat_df <- expand.grid(tx_a = cs$species_standardized, tx_b = cs$species_standardized, stringsAsFactors = FALSE)
    pl_m <- stats::setNames(suppressWarnings(as.numeric(cs$modal_ploidy_level)), cs$species_standardized)
    mat_df$barrier <- ifelse(is.finite(pl_m[mat_df$tx_a]) & is.finite(pl_m[mat_df$tx_b]) & pl_m[mat_df$tx_a] != pl_m[mat_df$tx_b], "Ploidy_Barrier", "Compatible")
    ggplot2::ggplot(mat_df, ggplot2::aes(x = .data$tx_a, y = .data$tx_b, fill = .data$barrier)) +
      ggplot2::geom_tile(color = "#FFFFFF") + ggplot2::scale_fill_manual(values = c(Ploidy_Barrier = "#D55E00", Compatible = "#009E73"), name = "Cytotype") +
      ggplot2::labs(title = "D. Cytotype Reproductive Barrier Matrix", x = NULL, y = NULL) + theme_packera_pub(base_size = 8.2) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, face = "italic"), axis.text.y = ggplot2::element_text(face = "italic"))
  } else .tier_empty("D. Cytotype Reproductive Barrier Matrix", "Fewer than 2 taxa with modal ploidy")
  save_packera_png_300((p6c_a | p6c_b) / (p6c_c | p6c_d) + patchwork::plot_annotation(title = "Step 06 [Tier 4 - Cytology & Sampling Gap]: Cytotype Barriers & Voucher Audit"), fig_c, width = 15.0, height = 10.5, dpi = 300L)
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a

render_step06_tier_plates_png <- function(active_tier, fig_a, fig_b, fig_c, tier_payload, db_path = "", run_id = "canonical") {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") || identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) return(invisible(""))
  for (fp in c(fig_a, fig_b, fig_c)) if (nzchar(fp)) dir.create(dirname(fp), recursive = TRUE, showWarnings = FALSE)
  if (identical(active_tier, "Tier_2_Public_HybSeq")) .render_tier2_plates(fig_a, fig_b, fig_c, tier_payload, db_path, run_id)
  else if (identical(active_tier, "Tier_3_Public_MultiLocus")) .render_tier3_plates(fig_a, fig_b, fig_c, tier_payload, db_path, run_id)
  else .render_tier4_plates(fig_a, fig_b, fig_c, tier_payload, db_path, run_id)
  invisible(list(fig_a = fig_a, fig_b = fig_b, fig_c = fig_c))
}
