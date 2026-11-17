# ==============================================================================
# Module: scripts/figures/step07_figures.R
# Step 07 Canonical 300-DPI PNG Figures (5-Stage Integrative Taxonomic Synthesis
# & Unsupervised Leiden / Metric Ward.D2 Taxon-Complex Delimitation):
#   1. fig07a_integrative_taxonomic_synthesis.png
#   2. fig07b_taxon_complex_delimitation_network.png
# ==============================================================================
if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

for (.nf_p in c("scripts/figures/step07b_network_figures.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "step07b_network_figures.R"))) {
  if (!is.na(.nf_p) && nzchar(.nf_p) && file.exists(.nf_p)) { source(.nf_p); break }
}

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

.build_cytology_fig_subtitle <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return("")
  vp_x <- if ("voucher_ploidy_level" %in% names(df)) suppressWarnings(as.integer(df$voucher_ploidy_level)) else rep(NA_integer_, nrow(df))
  vp_2n <- if ("voucher_sporophytic_2n" %in% names(df)) suppressWarnings(as.integer(df$voucher_sporophytic_2n)) else rep(NA_integer_, nrow(df))
  tx_cy <- if ("taxon_reported_cytotypes" %in% names(df)) trimws(as.character(df$taxon_reported_cytotypes)) else rep("", nrow(df))
  tx_cy[is.na(tx_cy) | tx_cy %in% c("NA", "NULL", "None")] <- ""
  has_x <- is.finite(vp_x) & vp_x > 0L
  has_2n <- is.finite(vp_2n) & vp_2n > 0L
  has_tx <- nzchar(tx_cy)
  if (!any(has_x) && !any(has_tx)) return("")
  n_cyto_v <- sum(has_x | has_2n)
  n_bar <- if ("cytotype_reproductive_barrier_pass" %in% names(df)) {
    sum(!is.na(df$cytotype_reproductive_barrier_pass) & suppressWarnings(as.integer(df$cytotype_reproductive_barrier_pass)) == 1L)
  } else 0L
  u_2n <- sort(unique(vp_2n[has_2n]))
  u_x <- sort(unique(vp_x[has_x]))
  lbl_x <- if (length(u_x) > 0L) paste0(paste(utils::head(u_x, 3L), collapse = "/"), "x") else paste(utils::head(unique(tx_cy[has_tx]), 2L), collapse = "; ")
  lbl_ploidy <- if (length(u_2n) > 0L) sprintf("2n=%s | %s", paste(utils::head(u_2n, 3L), collapse = "/"), lbl_x) else lbl_x
  sprintf("[%s; %d cytovoucher(s), %d barrier(s)]", lbl_ploidy, n_cyto_v, n_bar)
}

render_step07a_integrative_synthesis_png <- function(
  df,
  out_plot_path = "outputs/figures/fig07a_integrative_taxonomic_synthesis.png",
  cda_summary = NULL,
  focal_taxon = ""
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (!is.null(out_plot_path) && nzchar(out_plot_path)) dirname(out_plot_path) else "outputs/figures"
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

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
  } else rep(0.0, nrow(tdf))
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
  tdf_map <- tdf[!.is_step07_centroid_coord(tdf), , drop = FALSE]
  p_map <- if (nrow(tdf_map) > 0L && all(c(lon_col, lat_col) %in% names(tdf_map))) {
    m_ctx <- build_na_projected_basemap_gg(tdf_map[[lon_col]], tdf_map[[lat_col]])
    pts_g <- project_lonlat_points_r(tdf_map, m_ctx, lon_col = lon_col, lat_col = lat_col)
    pm <- m_ctx$plot
    has_halo_anchors <- FALSE
    if (nrow(pts_g) > 0L) {
      g_halo_tiers <- c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Conditional_Deme_Hull")
      c_halo_tiers <- c("C1_Direct_Cytovoucher", "C2_Exsiccatae_Cytoduplicate")
      gl_g <- if ("geno_linkage_tier" %in% names(pts_g)) trimws(as.character(pts_g$geno_linkage_tier)) else if ("linkage_tier" %in% names(pts_g)) trimws(as.character(pts_g$linkage_tier)) else rep("", nrow(pts_g))
      gp_g <- if ("genomic_propagation_tier" %in% names(pts_g)) trimws(as.character(pts_g$genomic_propagation_tier)) else rep("", nrow(pts_g))
      cy_g <- if ("cytology_linkage_tier" %in% names(pts_g)) trimws(as.character(pts_g$cytology_linkage_tier)) else rep("", nrow(pts_g))
      geno_halo_pts <- pts_g[(!is.na(gl_g) & gl_g %in% g_halo_tiers) | (!is.na(gp_g) & gp_g %in% g_halo_tiers), , drop = FALSE]
      cyto_halo_pts <- pts_g[!is.na(cy_g) & cy_g %in% c_halo_tiers, , drop = FALSE]
      has_halo_anchors <- (nrow(geno_halo_pts) > 0L) || (nrow(cyto_halo_pts) > 0L)
      if (nrow(geno_halo_pts) > 0L) {
        pm <- pm + ggplot2::geom_point(
          data = geno_halo_pts, ggplot2::aes(x = .data$map_x, y = .data$map_y),
          inherit.aes = FALSE, shape = 1L, fill = NA, size = 4.0, stroke = 0.7, color = "#1E293B", show.legend = FALSE
        )
      }
      if (nrow(cyto_halo_pts) > 0L) {
        pm <- pm + ggplot2::geom_point(
          data = cyto_halo_pts, ggplot2::aes(x = .data$map_x, y = .data$map_y),
          inherit.aes = FALSE, shape = 5L, fill = NA, size = 4.0, stroke = 0.7, color = "#1E293B", show.legend = FALSE
        )
      }
      pm <- pm +
        ggplot2::geom_point(
          data = pts_g,
          ggplot2::aes(x = .data$map_x, y = .data$map_y, color = .data$treatment_call, shape = .data$treatment_call),
          size = 2.3, stroke = 0.85, alpha = 0.85, show.legend = FALSE
        ) +
        ggplot2::scale_color_manual(values = PACKERA_TREATMENT_PALETTE) +
        ggplot2::scale_shape_manual(values = PACKERA_TREATMENT_HOLLOW_SHAPES)
    }
    map_sub <- if (has_halo_anchors) "NAD83 / Conus Albers Equal Area (EPSG:5070) (O = G1/G2 DNA, <> = C1/C2 2n)" else "NAD83 / Conus Albers Equal Area (EPSG:5070)"
    pm + ggplot2::labs(
      title = "B. Geographic Contact Zones, Anomalies & Type Localities",
      subtitle = map_sub
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
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8.0)) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2L, title = "5-Stage Call")) +
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
  wg_raw <- if ("w_geno_kappa" %in% names(tdf)) suppressWarnings(as.numeric(tdf$w_geno_kappa)) else numeric(0)
  wg_pos <- wg_raw[is.finite(wg_raw) & wg_raw > 0.0]
  w_g <- if (length(wg_pos) > 0L) mean(wg_pos) else if (any(is.finite(wg_raw))) mean(wg_raw[is.finite(wg_raw)]) else NA_real_

  active_geno_tiers <- c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", "G2_Conditional_Deme_Hull")
  .active_geno_mask_07 <- function(d) {
    if (!("geno_predicted_taxon" %in% names(d))) return(rep(FALSE, nrow(d)))
    gp_v <- trimws(as.character(d$geno_predicted_taxon))
    gl_v <- if ("geno_linkage_tier" %in% names(d)) trimws(as.character(d$geno_linkage_tier)) else if ("linkage_tier" %in% names(d)) trimws(as.character(d$linkage_tier)) else rep("", nrow(d))
    pt_v <- if ("genomic_propagation_tier" %in% names(d)) trimws(as.character(d$genomic_propagation_tier)) else rep("", nrow(d))
    !is.na(gp_v) & nzchar(gp_v) & gp_v != "Unknown" & ((!is.na(gl_v) & gl_v %in% active_geno_tiers) | (!is.na(pt_v) & pt_v %in% active_geno_tiers))
  }
  has_geno_stream <- any(.active_geno_mask_07(tdf))

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

    # Independent stream 5/6: Mean probabilistic_consensus_score (scaled 0-100) + optional active Genomic stream
    pcs_v <- if ("probabilistic_consensus_score" %in% names(sub_tx)) suppressWarnings(as.numeric(sub_tx$probabilistic_consensus_score)) else rep(NA_real_, n_tx)
    pct_cons <- if (any(is.finite(pcs_v))) 100.0 * mean(pmin(1.0, pmax(0.0, pcs_v)), na.rm = TRUE) else 100.0 * mean(sub_tx$treatment_call %in% c("Confirmed_Species", "Confirmed_Type_Control"))

    if (has_geno_stream) {
      g_ok <- .active_geno_mask_07(sub_tx)
      pct_geno <- if (any(g_ok)) {
        if (is_single_sp_06) {
          adm_sub <- if ("is_geno_admixed" %in% names(sub_tx)) (!is.na(sub_tx$is_geno_admixed) & toupper(trimws(as.character(sub_tx$is_geno_admixed))) %in% c("TRUE", "1", "T", "YES")) else rep(FALSE, n_tx)
          100.0 * mean(!adm_sub[g_ok], na.rm = TRUE)
        } else {
          100.0 * mean(trimws(as.character(sub_tx$geno_predicted_taxon[g_ok])) == trimws(as.character(sub_tx$species_standardized[g_ok])), na.rm = TRUE)
        }
      } else NA_real_
      stream_lvls <- c("1. Leaf Morphometrics", "2. Ecological & Edaphic", "3. Topological ECT & OSR", "4. Circular Phenology", "5. Genomic & Lineage", "6. Consensus Score")
      score_vec <- c(pct_morph, pct_eco, pct_ect_osr, pct_pheno, pct_geno, pct_cons)
    } else {
      stream_lvls <- c("1. Leaf Morphometrics", "2. Ecological & Edaphic", "3. Topological ECT & OSR", "4. Circular Phenology", "5. Consensus Score")
      score_vec <- c(pct_morph, pct_eco, pct_ect_osr, pct_pheno, pct_cons)
    }

    data.frame(
      plot_group = tx,
      stream = factor(stream_lvls, levels = stream_lvls),
      score = score_vec,
      stringsAsFactors = FALSE
    )
  })
  wett_df <- do.call(rbind, wett_rows)
  wett_df$lbl <- ifelse(is.finite(wett_df$score), sprintf("%.0f%%", wett_df$score), "N/A")
  wett_sub <- if (has_geno_stream) {
    sprintf("OOF \u03ba weights: Morph = %.2f | Edaphic = %.2f | Phenology = %.2f | Genomic = %s", w_m, w_e, w_p, ifelse(is.finite(w_g), sprintf("%.2f", w_g), "N/A"))
  } else {
    sprintf("OOF \u03ba weights: Morph = %.2f | Edaphic = %.2f | Phenology = %.2f", w_m, w_e, w_p)
  }

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
      subtitle = wett_sub,
      x = NULL, y = NULL
    )

  dyn_h_06a <- max(10.2, 6.2 + 0.28 * length(active_taxa))
  hdr_title <- if (length(active_taxa) >= 2L && nzchar(trimws(as.character(focal_taxon)))) {
    sprintf("Step 07: Integrative Species Delimitation Across Packera Complexes (%s Focus, N = %d Taxa)", trimws(as.character(focal_taxon)), length(active_taxa))
  } else if (nzchar(trimws(as.character(focal_taxon)))) {
    sprintf("Step 07: Integrative Species Delimitation & Curatorial Synthesis (%s)", trimws(as.character(focal_taxon)))
  } else if (length(active_taxa) == 1L) {
    sprintf("Step 07: Integrative Species Delimitation & Curatorial Synthesis (%s)", active_taxa[1L])
  } else {
    sprintf("Step 07: Integrative Species Delimitation Across Packera Complexes (N = %d Taxa)", length(active_taxa))
  }
  prop_tier_v <- if ("genomic_propagation_tier" %in% names(df)) trimws(as.character(df$genomic_propagation_tier)) else rep(NA_character_, nrow(df))
  t0_tx_v <- if ("tier0_calibrated_taxon" %in% names(df)) trimws(as.character(df$tier0_calibrated_taxon)) else rep(NA_character_, nrow(df))
  t0_abs_v <- if ("tier0_abstention_flag" %in% names(df)) suppressWarnings(as.integer(df$tier0_abstention_flag)) else rep(NA_integer_, nrow(df))
  cl_gdi_v <- if ("cluster_propagated_gdi" %in% names(df)) suppressWarnings(as.numeric(df$cluster_propagated_gdi)) else rep(NA_real_, nrow(df))
  n_exs_dup <- sum(!is.na(prop_tier_v) & prop_tier_v == "G1_Exsiccatae_Isoduplicate")
  n_deme_hull <- sum(!is.na(prop_tier_v) & prop_tier_v == "G2_Conditional_Deme_Hull")
  n_t0_cal <- sum(!is.na(t0_tx_v) & nzchar(t0_tx_v) & t0_tx_v != "Sensu_Lato_DNA_Required" & !is.na(t0_abs_v) & t0_abs_v == 0L)
  n_cl_gdi <- sum(is.finite(cl_gdi_v))
  n_dna_abst <- sum(!is.na(t0_abs_v) & t0_abs_v == 1L)
  has_prop_fig <- (n_exs_dup + n_deme_hull + n_t0_cal + n_cl_gdi + n_dna_abst) > 0L
  cyto_sub_07a <- .build_cytology_fig_subtitle(df)
  has_sub_07a <- has_prop_fig || nzchar(cyto_sub_07a)

  plate_06a <- if (has_sub_07a) {
    sub_parts <- character(0)
    if (has_prop_fig) {
      sub_parts <- c(sub_parts, sprintf(
        "[Exs-Dup: %d | Deme-Hull: %d | Tier0-Cal: %d | Cluster-gdi: %d | DNA-Abstain: %d]",
        n_exs_dup, n_deme_hull, n_t0_cal, n_cl_gdi, n_dna_abst
      ))
    }
    if (nzchar(cyto_sub_07a)) sub_parts <- c(sub_parts, cyto_sub_07a)
    prop_sub <- paste(sub_parts, collapse = " ")
    ((p_morph | p_map) / (p_treat | p_wett)) +
      patchwork::plot_layout(heights = c(1.12, 1.05)) +
      patchwork::plot_annotation(
        title = hdr_title,
        subtitle = prop_sub,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 12.5),
          plot.subtitle = ggplot2::element_text(size = 9.0, color = "#374151")
        )
      )
  } else {
    ((p_morph | p_map) / (p_treat | p_wett)) +
      patchwork::plot_layout(heights = c(1.12, 1.05)) +
      patchwork::plot_annotation(
        title = hdr_title,
        theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.5))
      )
  }

  canon_07a <- file.path(fdir, "fig07a_integrative_taxonomic_synthesis.png")
  extra_targets <- character(0)
  if (!is.null(out_plot_path) && nzchar(out_plot_path)) extra_targets <- c(extra_targets, out_plot_path)
  save_packera_png_300(plate_06a, canon_07a, width = 14.4, height = dyn_h_06a, also_copy_to = extra_targets)
}
