# ==============================================================================
# Module: scripts/figures/step06_figures.R
# Step 06 Canonical 300-DPI PNG Figures (5-Stage Integrative Taxonomic Synthesis
# & Unsupervised Leiden / Metric Ward.D2 Taxon-Complex Delimitation):
#   1. fig06a_integrative_taxonomic_synthesis.png
#   2. fig06b_taxon_complex_delimitation_network.png
# ==============================================================================

PACKERA_TREATMENT_PALETTE <- c(
  "Confirmed_Species"        = "#0072B2",
  "Confirmed_Type_Control"   = "#332288",
  "Subspecies_Ecotype"       = "#009E73",
  "Ecophenotypic_Plasticity" = "#56B4E9",
  "Hybrid_Introgressant"     = "#E69F00",
  "Putative_Cryptic_Species" = "#CC79A7",
  "Type_Synonymy_Candidate"  = "#882255",
  "Misidentification_Reassign" = "#D55E00"
)

PACKERA_TREATMENT_HOLLOW_SHAPES <- c(
  "Confirmed_Species"        = 1L,
  "Confirmed_Type_Control"   = 8L,
  "Subspecies_Ecotype"       = 2L,
  "Ecophenotypic_Plasticity" = 0L,
  "Hybrid_Introgressant"     = 5L,
  "Putative_Cryptic_Species" = 6L,
  "Type_Synonymy_Candidate"  = 4L,
  "Misidentification_Reassign" = 3L
)

render_step06a_integrative_synthesis_png <- function(
  df,
  out_plot_path = "outputs/figures/fig06a_integrative_taxonomic_synthesis.png",
  cda_summary = NULL,
  focal_taxon = ""
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (!is.null(out_plot_path) && nzchar(out_plot_path)) dirname(out_plot_path) else "outputs/figures"
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (old_stem in c("fig06_multimodal_discordance_summary", "fig06b_taxonomic_treatment_summary", "Figure_Integrative_Packera_Revision")) {
    for (ext in c(".pdf", ".png")) {
      p_old <- file.path(fdir, paste0(old_stem, ext))
      if (file.exists(p_old) && normalizePath(p_old, mustWork = FALSE) != normalizePath(normalize_png_path(out_plot_path), mustWork = FALSE)) {
        unlink(p_old, force = TRUE)
      }
    }
  }

  tdf <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(tdf) == 0L) return(invisible(""))
  sp_u_06 <- unique(stats::na.omit(as.character(tdf$species_standardized)))
  sp_u_06 <- sp_u_06[nzchar(trimws(sp_u_06)) & sp_u_06 != "Unknown"]
  grp_col_06 <- if (length(sp_u_06) < 2L && "ecotype_regime" %in% names(tdf) && length(unique(stats::na.omit(tdf$ecotype_regime[grepl("^Ecotype_", as.character(tdf$ecotype_regime))]))) >= 2L) {
    "ecotype_regime"
  } else if (length(sp_u_06) >= 2L) {
    "species_standardized"
  } else NULL
  tdf$plot_group <- resolve_plot_group_labels(tdf, group_col = grp_col_06)
  active_taxa <- sort(unique(tdf$plot_group[!is.na(tdf$plot_group) & nzchar(tdf$plot_group) & tdf$plot_group != "Unknown"]))
  if (length(active_taxa) == 0L) active_taxa <- unique(tdf$plot_group)

  sub_call <- if ("taxonomic_treatment_subtype" %in% names(tdf)) as.character(tdf$taxonomic_treatment_subtype) else if ("taxonomic_status_call" %in% names(tdf)) as.character(tdf$taxonomic_status_call) else rep("Confirmed_Species", nrow(tdf))
  sub_call[is.na(sub_call) | !nzchar(trimws(sub_call))] <- "Confirmed_Species"
  sub_call[!(sub_call %in% names(PACKERA_TREATMENT_PALETTE))] <- "Confirmed_Species"
  tdf$treatment_call <- factor(sub_call, levels = intersect(names(PACKERA_TREATMENT_PALETTE), unique(sub_call)))

  # Panel A: 2D Integrative Morphospace Stratified by 5-Stage Decision Calls
  x1 <- if ("can1" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$can1))))) suppressWarnings(as.numeric(tdf$can1)) else if ("PC1" %in% names(tdf)) suppressWarnings(as.numeric(tdf$PC1)) else rep(0.0, nrow(tdf))
  x2 <- if ("can2" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$can2)))) && isTRUE(stats::sd(suppressWarnings(as.numeric(tdf$can2)), na.rm = TRUE) > 1e-6)) {
    suppressWarnings(as.numeric(tdf$can2))
  } else if ("PC2" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$PC2))))) {
    suppressWarnings(as.numeric(tdf$PC2))
  } else rnorm(nrow(tdf), 0, 0.01)
  tdf$mx <- x1; tdf$my <- x2
  sub_m <- tdf[is.finite(tdf$mx) & is.finite(tdf$my), , drop = FALSE]

  p_morph <- ggplot2::ggplot(sub_m, ggplot2::aes(x = .data$mx, y = .data$my, color = .data$treatment_call, shape = .data$treatment_call)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#9CA3AF", linewidth = 0.4) +
    ggplot2::geom_point(size = 2.5, stroke = 0.9, alpha = 0.84) +
    ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE, labels = function(x) gsub("_", " ", x), name = "5-Stage Call") +
    ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES, labels = function(x) gsub("_", " ", x), name = "5-Stage Call") +
    theme_packera_pub(base_size = 9.6) +
    ggplot2::labs(
      title = "A. Integrative Morphospace & Nomenclatural Type Baseline",
      subtitle = "Hollow symbols encode 5-stage decision calls (* = Type control)",
      x = if ("can1" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$can1))))) "Canonical Variate 1 (Can1)" else "Morphometric PC1",
      y = if ("can2" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$can2))))) "Canonical Variate 2 (Can2)" else "Morphometric PC2"
    )

  # Panel B: Projected North American Geographic Contact Zones & Triage Map
  lon_col <- if ("decimalLongitude" %in% names(tdf)) "decimalLongitude" else "longitude"
  lat_col <- if ("decimalLatitude" %in% names(tdf)) "decimalLatitude" else "latitude"
  p_map <- if (all(c(lon_col, lat_col) %in% names(tdf))) {
    m_ctx <- build_na_projected_basemap_gg(tdf[[lon_col]], tdf[[lat_col]])
    pts_g <- project_lonlat_points_r(tdf, m_ctx, lon_col = lon_col, lat_col = lat_col)
    pm <- m_ctx$plot
    if (nrow(pts_g) > 0L) {
      pm <- pm +
        ggplot2::geom_point(
          data = pts_g,
          ggplot2::aes(x = .data$map_x, y = .data$map_y, color = .data$treatment_call, shape = .data$treatment_call),
          size = 2.3, stroke = 0.85, alpha = 0.85, show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE) +
        ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES)
    }
    pm + ggplot2::labs(
      title = "B. Geographic Contact Zones, Anomalies & Type Localities",
      subtitle = "NAD83 / Conus Albers Equal Area (EPSG:5070)"
    )
  } else {
    ggplot2::ggplot() + theme_packera_pub() + ggplot2::labs(title = "B. Geographic Contact Zones")
  }

  # Panel C: 5-Stage Taxonomic Treatment Proportions by Taxon
  cnt_tx <- table(tdf$plot_group)
  bar_df <- as.data.frame(table(plot_group = tdf$plot_group, treatment_call = tdf$treatment_call), stringsAsFactors = FALSE)
  bar_df <- bar_df[bar_df$Freq > 0, , drop = FALSE]
  bar_df$treatment_call <- factor(bar_df$treatment_call, levels = names(PACKERA_TREATMENT_PALETTE))

  p_treat <- ggplot2::ggplot(bar_df, ggplot2::aes(x = .data$Freq, y = .data$plot_group, fill = .data$treatment_call)) +
    ggplot2::geom_col(position = ggplot2::position_fill(reverse = TRUE), width = 0.65, color = "white", linewidth = 0.25) +
    ggplot2::scale_x_continuous(labels = function(x) paste0(round(100 * x), "%"), expand = c(0, 0)) +
    ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE, counts_map = cnt_tx)) +
    ggplot2::scale_fill_manual(values = PACKERA_TREATMENT_PALETTE, labels = function(x) gsub("_", " ", x), name = "5-Stage Call", drop = TRUE) +
    theme_packera_pub(base_size = 9.5) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      title = "C. Hierarchical 5-Stage Treatment Outcomes",
      subtitle = "Proportion of specimens per taxon across integrative categories",
      x = "Proportion of Evaluated Specimens (%)",
      y = NULL
    )

  # Panel D: Independent Multi-Evidence Wettstein Concordance Matrix
  w_m <- if ("w_morph_kappa" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$w_morph_kappa))))) stats::median(suppressWarnings(as.numeric(tdf$w_morph_kappa)), na.rm = TRUE) else 0.50
  w_e <- if ("w_edaphic_kappa" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$w_edaphic_kappa))))) stats::median(suppressWarnings(as.numeric(tdf$w_edaphic_kappa)), na.rm = TRUE) else 0.30
  w_p <- if ("w_pheno_kappa" %in% names(tdf) && any(is.finite(suppressWarnings(as.numeric(tdf$w_pheno_kappa))))) stats::median(suppressWarnings(as.numeric(tdf$w_pheno_kappa)), na.rm = TRUE) else 0.20

  wett_rows <- lapply(active_taxa, function(tx) {
    sub_tx <- tdf[tdf$plot_group == tx, , drop = FALSE]
    n_tx <- max(1L, nrow(sub_tx))
    is_single_sp_06 <- length(sp_u_06) < 2L
    # Independent stream 1: Shrinkage CDA (multi-taxon) or Self/Type Morphospace Cohesion (single-taxon)
    if (is_single_sp_06) {
      mo_f <- if ("morphometric_outlier_flag" %in% names(sub_tx)) as.logical(sub_tx$morphometric_outlier_flag) else rep(FALSE, n_tx)
      mo_f[is.na(mo_f)] <- FALSE
      sd2 <- if ("self_morph_mahal_d2" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$self_morph_mahal_d2)) else rep(NA_real_, n_tx)
      sdf <- if ("self_morph_df" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$self_morph_df)) else rep(6.0, n_tx)
      sdf[!is.finite(sdf) | sdf < 1.0] <- 6.0
      in_self95 <- ifelse(is.finite(sd2), sd2 <= stats::qchisq(0.95, df = sdf), !mo_f)
      pct_morph <- 100.0 * mean(!mo_f & in_self95, na.rm = TRUE)
    } else {
      c_pred <- if ("cda_predicted_taxon" %in% names(sub_tx)) as.character(sub_tx$cda_predicted_taxon) else rep(NA_character_, n_tx)
      m_ok <- !is.na(c_pred) & nzchar(c_pred)
      pct_morph <- if (any(m_ok)) 100.0 * mean(c_pred[m_ok] == sub_tx$species_standardized[m_ok], na.rm = TRUE) else NA_real_
    }

    # Independent stream 2: Spatial RF / Edaphic niche concordance (against Ecotype_* in single-species multi-regime mode)
    r_pred <- if ("rf_predicted_taxon" %in% names(sub_tx)) as.character(sub_tx$rf_predicted_taxon) else if ("edaphic_best_fit_taxon" %in% names(sub_tx)) as.character(sub_tx$edaphic_best_fit_taxon) else rep(NA_character_, n_tx)
    ez_v <- if ("edaphic_z_score" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$edaphic_z_score)) else rep(NA_real_, n_tx)
    if (is_single_sp_06) {
      er_ref <- if ("ecotype_regime" %in% names(sub_tx)) as.character(sub_tx$ecotype_regime) else if ("analysis_group" %in% names(sub_tx)) as.character(sub_tx$analysis_group) else as.character(sub_tx$species_standardized)
      has_multi_er <- any(grepl("^Ecotype_", er_ref), na.rm = TRUE)
      if (has_multi_er && any(!is.na(r_pred) & nzchar(r_pred))) {
        e_ok <- !is.na(r_pred) & nzchar(r_pred) & !is.na(er_ref) & nzchar(er_ref)
        pct_eco <- if (any(e_ok)) 100.0 * mean(r_pred[e_ok] == er_ref[e_ok], na.rm = TRUE) else NA_real_
      } else if (any(is.finite(ez_v))) {
        pct_eco <- 100.0 * mean(ez_v[is.finite(ez_v)] <= 2.5, na.rm = TRUE)
      } else {
        pct_eco <- NA_real_
      }
    } else {
      e_ok <- !is.na(r_pred) & nzchar(r_pred)
      pct_eco <- if (any(e_ok)) 100.0 * mean(r_pred[e_ok] == sub_tx$species_standardized[e_ok], na.rm = TRUE) else NA_real_
    }

    # Independent stream 3: Topological ECT & Open-Set Recognition (OSR) in-distribution support
    osr_p <- if ("open_set_novelty_p" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$open_set_novelty_p)) else rep(NA_real_, n_tx)
    out_f <- if ("morphometric_outlier_flag" %in% names(sub_tx)) as.logical(sub_tx$morphometric_outlier_flag) else rep(NA, n_tx)
    pct_ect_osr <- if (any(is.finite(osr_p))) 100.0 * mean(osr_p >= 0.05, na.rm = TRUE) else if (any(!is.na(out_f))) 100.0 * mean(!out_f, na.rm = TRUE) else pct_morph

    # Independent stream 4: Circular Hopkins phenological synchrony (|delta_doy_resid| <= 21 d)
    dd_r <- if ("delta_doy_resid" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$delta_doy_resid)) else if ("delta_doy" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$delta_doy)) else rep(NA_real_, n_tx)
    pct_pheno <- if (any(is.finite(dd_r))) 100.0 * mean(abs(dd_r) <= 21.0, na.rm = TRUE) else NA_real_

    # Independent stream 5: Mean probabilistic_consensus_score (scaled 0-100)
    pcs_v <- if ("probabilistic_consensus_score" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$probabilistic_consensus_score)) else rep(NA_real_, n_tx)
    pct_cons <- if (any(is.finite(pcs_v))) 100.0 * mean(pmin(1.0, pmax(0.0, pcs_v)), na.rm = TRUE) else 100.0 * mean(sub_tx$treatment_call %in% c("Confirmed_Species", "Confirmed_Type_Control"))

    data.frame(
      plot_group = tx,
      stream = factor(
        c("1. Leaf Morphometrics", "2. Ecological & Edaphic", "3. Topological ECT & OSR", "4. Circular Phenology", "5. Consensus Score"),
        levels = c("1. Leaf Morphometrics", "2. Ecological & Edaphic", "3. Topological ECT & OSR", "4. Circular Phenology", "5. Consensus Score")
      ),
      score = c(pct_morph, pct_eco, pct_ect_osr, pct_pheno, pct_cons),
      stringsAsFactors = FALSE
    )
  })
  wett_df <- do.call(rbind, wett_rows)
  wett_df$lbl <- ifelse(is.finite(wett_df$score), sprintf("%.0f%%", wett_df$score), "N/A")

  p_wett <- ggplot2::ggplot(wett_df, ggplot2::aes(x = .data$stream, y = .data$plot_group, fill = .data$score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$lbl),
      size = if (length(active_taxa) > 15L) 2.2 else 2.75,
      color = ifelse(!is.finite(wett_df$score) | wett_df$score > 45, "#111827", "white")
    ) +
    ggplot2::scale_fill_viridis_c(option = "cividis", limits = c(0, 100), na.value = "#F3F4F6", name = "Support %") +
    ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
    theme_packera_pub(base_size = 9.4) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)) +
    ggplot2::labs(
      title = "D. Multi-Evidence Wettstein Concordance Grid",
      subtitle = sprintf("OOF \u03ba weights: Morph = %.2f | Edaphic = %.2f | Phenology = %.2f", w_m, w_e, w_p),
      x = NULL, y = NULL
    )

  dyn_h_06a <- max(10.2, 6.2 + 0.28 * length(active_taxa))
  hdr_title <- if (nzchar(trimws(as.character(focal_taxon)))) {
    sprintf("Step 06: Integrative Species Delimitation & Curatorial Synthesis (%s)", trimws(as.character(focal_taxon)))
  } else {
    sprintf("Step 06: Integrative Species Delimitation Across Packera Complexes (N = %d Taxa)", length(active_taxa))
  }
  plate_06a <- ((p_morph | p_map) / (p_treat | p_wett)) +
    patchwork::plot_layout(heights = c(1.12, 1.05)) +
    patchwork::plot_annotation(
      title = hdr_title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5))
    )

  canon_06a <- file.path(fdir, "fig06a_integrative_taxonomic_synthesis.png")
  extra_targets <- character(0)
  if (!is.null(out_plot_path) && nzchar(out_plot_path)) extra_targets <- c(extra_targets, out_plot_path)
  if (nzchar(trimws(as.character(focal_taxon)))) {
    f_slug <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", trimws(as.character(focal_taxon))))
    if (nzchar(f_slug)) {
      focal_old_pdf <- file.path(fdir, sprintf("Figure_Integrative_%s_Revision.pdf", f_slug))
      if (file.exists(focal_old_pdf)) unlink(focal_old_pdf, force = TRUE)
    }
  }
  save_packera_png_300(plate_06a, canon_06a, width = 14.4, height = dyn_h_06a, also_copy_to = extra_targets)
}

render_step06b_complex_network_png <- function(
  complexes_df,
  affinity_mat,
  out_complex_path = "outputs/figures/fig06b_taxon_complex_delimitation_network.png",
  min_affinity = 0.03
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (!is.null(out_complex_path) && nzchar(out_complex_path)) dirname(out_complex_path) else "outputs/figures"
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  for (ext in c(".pdf", ".png")) {
    p_old <- file.path(fdir, paste0("fig06c_taxon_complex_network", ext))
    if (file.exists(p_old) && normalizePath(p_old, mustWork = FALSE) != normalizePath(normalize_png_path(out_complex_path), mustWork = FALSE)) {
      unlink(p_old, force = TRUE)
    }
  }

  cdf <- as.data.frame(complexes_df, stringsAsFactors = FALSE)
  W <- as.matrix(affinity_mat)
  if (nrow(cdf) == 0L || nrow(W) == 0L) return(invisible(""))
  taxa <- rownames(W)
  if (is.null(taxa)) taxa <- as.character(cdf$species_standardized)
  n_t <- length(taxa)

  w_vec <- suppressWarnings(as.numeric(W))
  w_vec[!is.finite(w_vec)] <- 0.0
  w_vec <- pmin(1.0, pmax(0.0, w_vec))
  W_num <- matrix(w_vec, nrow = n_t, ncol = n_t, dimnames = list(taxa, taxa))
  diag(W_num) <- 0.0

  # Order taxa by metric Ward.D2 hierarchical clustering so complexes form contiguous blocks
  ord_idx <- if (n_t >= 3L) {
    D_ward <- matrix(sqrt(pmax(1.0 - as.numeric(W_num), 0.0)), nrow = n_t, ncol = n_t, dimnames = list(taxa, taxa))
    diag(D_ward) <- 0.0
    tryCatch(stats::hclust(stats::as.dist(D_ward), method = "ward.D2")$order, error = function(e) seq_len(n_t))
  } else seq_len(n_t)
  taxa_ord <- taxa[ord_idx]
  W_disp <- W_num
  diag(W_disp) <- 1.0
  W_ord <- W_disp[ord_idx, ord_idx, drop = FALSE]

  # Panel A: Ward.D2 Clustered Multimodal Affinity Heatmap
  hm_grid <- expand.grid(Row = taxa_ord, Col = taxa_ord, stringsAsFactors = FALSE)
  hm_grid$Affinity <- pmin(1.0, pmax(0.0, as.vector(W_ord)))
  hm_grid$Row <- factor(hm_grid$Row, levels = rev(taxa_ord))
  hm_grid$Col <- factor(hm_grid$Col, levels = taxa_ord)

  sil_mean <- if ("silhouette_width" %in% names(cdf) && any(is.finite(suppressWarnings(as.numeric(cdf$silhouette_width))))) {
    mean(suppressWarnings(as.numeric(cdf$silhouette_width)), na.rm = TRUE)
  } else NA_real_
  n_comp <- if ("consensus_complex_id" %in% names(cdf)) length(unique(cdf$consensus_complex_id)) else 1L

  p_heat <- ggplot2::ggplot(hm_grid, ggplot2::aes(x = .data$Col, y = .data$Row, fill = .data$Affinity)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_fill_viridis_c(option = "magma", limits = c(0, 1), name = "Affinity\nW_ij") +
    ggplot2::scale_x_discrete(labels = function(x) abbreviate_icn_taxon(x, abbreviate_genus = TRUE)) +
    ggplot2::scale_y_discrete(labels = icn_discrete_labeller(abbreviate_genus = TRUE)) +
    theme_packera_pub(base_size = 9.4) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, face = "italic")) +
    ggplot2::labs(
      title = sprintf("A. Ward.D2 Multimodal Affinity Matrix (%d Complexes)", n_comp),
      subtitle = if (is.finite(sil_mean)) sprintf("Mean consensus silhouette width = %.3f", sil_mean) else "Morphometric + Ecological KDE + Circular Hopkins affinity",
      x = NULL, y = NULL
    )
  if (n_t <= 12L) {
    hm_sub <- hm_grid[as.character(hm_grid$Row) != as.character(hm_grid$Col) & is.finite(hm_grid$Affinity) & hm_grid$Affinity >= min_affinity, , drop = FALSE]
    if (nrow(hm_sub) > 0L) {
      hm_sub$lbl <- sprintf("%.2f", hm_sub$Affinity)
      p_heat <- p_heat + ggplot2::geom_text(
        data = hm_sub,
        ggplot2::aes(x = .data$Col, y = .data$Row, label = .data$lbl),
        size = 2.4,
        color = ifelse(hm_sub$Affinity > 0.55, "#111827", "white"),
        inherit.aes = FALSE
      )
    }
  }

  # Panel B: Multimodal Affinity Network Graph (Leiden/Consensus Complex Communities)
  idx_m <- match(taxa, as.character(cdf$species_standardized))
  comp_ids <- if ("consensus_complex_id" %in% names(cdf)) as.character(cdf$consensus_complex_id[idx_m]) else rep("1", n_t)
  comp_ids[is.na(comp_ids)] <- "1"
  comp_lbls <- paste0("Complex ", comp_ids)
  n_vouch <- if ("n_vouchers" %in% names(cdf)) suppressWarnings(as.numeric(cdf$n_vouchers[idx_m])) else rep(10.0, n_t)
  n_vouch[!is.finite(n_vouch) | n_vouch <= 0] <- 5.0

  # Deterministic classical MDS 2D embedding from metric chord dissimilarity sqrt(1 - W)
  coords <- matrix(0.0, nrow = n_t, ncol = 2L)
  if (n_t == 1L) {
    coords[1L, ] <- c(0.0, 0.0)
  } else if (n_t == 2L) {
    coords[, 1L] <- c(-0.6, 0.6); coords[, 2L] <- c(0.0, 0.0)
  } else {
    D_net <- matrix(sqrt(pmax(1.0 - as.numeric(W_num), 0.04)), nrow = n_t, ncol = n_t)
    diag(D_net) <- 0.0
    mds_fit <- tryCatch(suppressWarnings(stats::cmdscale(stats::as.dist(D_net), k = 2L)), error = function(e) NULL)
    if (!is.null(mds_fit) && is.matrix(mds_fit) && ncol(mds_fit) == 2L && all(is.finite(mds_fit)) &&
        stats::sd(mds_fit[, 1L]) > 1e-6 && stats::sd(mds_fit[, 2L]) > 1e-6) {
      coords <- mds_fit
    } else {
      ang <- seq(0, 2 * pi, length.out = n_t + 1L)[seq_len(n_t)]
      coords <- cbind(cos(ang), sin(ang))
    }
  }
  log_nv <- log10(n_vouch + 1.0)
  pt_sz <- if (max(log_nv) > min(log_nv) + 1e-6) {
    3.6 + 2.8 * (log_nv - min(log_nv)) / (max(log_nv) - min(log_nv))
  } else {
    rep(4.4, n_t)
  }
  node_df <- data.frame(
    taxon = taxa,
    lbl = abbreviate_icn_taxon(taxa, abbreviate_genus = TRUE),
    complex = comp_lbls,
    n_vouchers = n_vouch,
    pt_sz = pt_sz,
    x = as.numeric(coords[, 1L]),
    y = as.numeric(coords[, 2L]),
    stringsAsFactors = FALSE
  )

  edge_list <- list()
  if (n_t >= 2L) {
    for (i in seq_len(n_t - 1L)) {
      for (j in (i + 1L):n_t) {
        wij <- W_num[i, j]
        if (is.finite(wij) && wij >= min_affinity) {
          edge_list[[length(edge_list) + 1L]] <- data.frame(
            x = node_df$x[i], y = node_df$y[i],
            xend = node_df$x[j], yend = node_df$y[j],
            lw = 0.55 + 2.05 * pmin(1.0, pmax(0.0, wij)),
            edge_kind = ifelse(node_df$complex[i] == node_df$complex[j], "intra", "inter"),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  edge_df <- if (length(edge_list) > 0L) do.call(rbind, edge_list) else data.frame()

  pal_comp <- get_combinatorial_palette_r(sort(unique(node_df$complex)))
  shp_comp <- get_combinatorial_hollow_shapes_r(sort(unique(node_df$complex)))

  build_net_panel <- function(use_repel = TRUE) {
    pn <- ggplot2::ggplot(node_df, ggplot2::aes(x = .data$x, y = .data$y))
    if (nrow(edge_df) > 0L) {
      e_intra <- edge_df[edge_df$edge_kind == "intra", , drop = FALSE]
      e_inter <- edge_df[edge_df$edge_kind == "inter", , drop = FALSE]
      if (nrow(e_inter) > 0L) {
        pn <- pn + ggplot2::geom_segment(
          data = e_inter,
          ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
          linewidth = e_inter$lw, linetype = "dotted", color = "#9CA3AF", alpha = 0.65,
          show.legend = FALSE, inherit.aes = FALSE
        )
      }
      if (nrow(e_intra) > 0L) {
        pn <- pn + ggplot2::geom_segment(
          data = e_intra,
          ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
          linewidth = e_intra$lw, linetype = "solid", color = "#64748B", alpha = 0.78,
          show.legend = FALSE, inherit.aes = FALSE
        )
      }
    }
    pn <- pn +
      ggplot2::geom_point(
        ggplot2::aes(color = .data$complex, shape = .data$complex),
        size = node_df$pt_sz, stroke = 1.25
      ) +
      ggplot2::scale_color_manual(values = pal_comp, name = "Delimited Complex") +
      ggplot2::scale_shape_manual(values = shp_comp, name = "Delimited Complex") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.22)) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.22)) +
      theme_packera_pub(base_size = 9.5) +
      ggplot2::labs(
        title = sprintf("B. Multimodal Taxon-Complex Network (W_ij >= %.2f)", min_affinity),
        subtitle = "Solid edges connect intra-complex members; node size ~ log10(N)",
        x = "Multimodal MDS Axis 1", y = "Multimodal MDS Axis 2"
      )
    add_centroid_repel_labels_gg(pn, node_df, "x", "y", "lbl", color_col = "complex", size = 2.95, allow_ggrepel = use_repel)
  }

  dyn_h_06b <- max(6.6, 4.4 + 0.24 * n_t)
  canon_06b <- file.path(fdir, "fig06b_taxon_complex_delimitation_network.png")
  plate_06b <- (p_heat | build_net_panel(TRUE)) +
    patchwork::plot_layout(widths = c(1.08, 1.0)) +
    patchwork::plot_annotation(
      title = "Step 06: Unsupervised Quantitative Taxon-Complex Delimitation",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4))
    )
  out_saved <- save_packera_png_300(plate_06b, canon_06b, width = 14.2, height = dyn_h_06b, also_copy_to = out_complex_path)
  if (!file.exists(canon_06b) || isTRUE(file.info(canon_06b)$size == 0)) {
    p_heat_fb <- p_heat +
      ggplot2::scale_y_discrete(labels = function(x) abbreviate_icn_taxon(x, abbreviate_genus = TRUE)) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
    plate_06b_fb <- (p_heat_fb | build_net_panel(FALSE)) +
      patchwork::plot_layout(widths = c(1.08, 1.0)) +
      patchwork::plot_annotation(
        title = "Step 06: Unsupervised Quantitative Taxon-Complex Delimitation",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.4))
      )
    out_saved <- save_packera_png_300(plate_06b_fb, canon_06b, width = 14.2, height = dyn_h_06b, also_copy_to = out_complex_path)
  }
  out_saved
}
