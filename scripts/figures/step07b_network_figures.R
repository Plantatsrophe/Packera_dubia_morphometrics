if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

if (!exists(".is_step07_centroid_coord", mode = "function")) {
  .is_step07_centroid_coord <- function(df) {
    n <- nrow(df); if (n == 0L) return(logical(0))
    if (!("is_centroid_or_inst_coord" %in% names(df))) return(rep(FALSE, n))
    raw_c <- df$is_centroid_or_inst_coord
    c_flag <- if (is.logical(raw_c)) !is.na(raw_c) & raw_c
    else if (is.numeric(raw_c) || is.integer(raw_c)) !is.na(raw_c) & (raw_c != 0)
    else tolower(trimws(as.character(raw_c))) %in% c("1", "true", "t", "yes", "y")
    is_t <- if ("is_type_specimen" %in% names(df)) (!is.na(df$is_type_specimen) & as.logical(df$is_type_specimen)) else rep(FALSE, n)
    if (exists(".identify_type_specimen_vec", mode = "function")) {
      ti <- .identify_type_specimen_vec(df)
      if (!is.null(ti$is_type_specimen) && length(ti$is_type_specimen) == n) {
        it2 <- as.logical(ti$is_type_specimen); it2[is.na(it2)] <- FALSE; is_t <- is_t | it2
      }
    }
    c_flag & !is_t
  }
}

render_step07b_complex_network_png <- function(
  complexes_df,
  affinity_mat,
  out_complex_path = "outputs/figures/fig07b_taxon_complex_delimitation_network.png",
  min_affinity = 0.03,
  voucher_df = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  fdir <- if (!is.null(out_complex_path) && nzchar(out_complex_path)) dirname(out_complex_path) else "outputs/figures"
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

  cdf <- as.data.frame(complexes_df, stringsAsFactors = FALSE)
  W <- as.matrix(affinity_mat)
  is_deme_net <- FALSE
  if ((nrow(cdf) <= 1L || nrow(W) <= 1L) && !is.null(voucher_df) && is.data.frame(voucher_df) && nrow(voucher_df) >= 6L) {
    is_cent_v <- .is_step07_centroid_coord(voucher_df)
    vdf_d <- voucher_df[!is_cent_v, , drop = FALSE]
    if (nrow(vdf_d) >= 6L) {
      deme_v <- if ("population_deme_id" %in% names(vdf_d) && sum(!is.na(vdf_d$population_deme_id) & nzchar(as.character(vdf_d$population_deme_id))) >= 6L) {
        as.character(vdf_d$population_deme_id)
      } else if ("stateProvince" %in% names(vdf_d) && length(unique(stats::na.omit(as.character(vdf_d$stateProvince)))) >= 2L) {
        as.character(vdf_d$stateProvince)
      } else {
        lat_d <- suppressWarnings(as.numeric(if ("decimalLatitude" %in% names(vdf_d)) vdf_d$decimalLatitude else vdf_d$latitude))
        lon_d <- suppressWarnings(as.numeric(if ("decimalLongitude" %in% names(vdf_d)) vdf_d$decimalLongitude else vdf_d$longitude))
        ifelse(is.finite(lat_d) & is.finite(lon_d), sprintf("%.1f,%.1f", round(lat_d, 1), round(lon_d, 1)), NA_character_)
      }
      vdf_d$.deme_lbl <- trimws(deme_v)
      vdf_d <- vdf_d[!is.na(vdf_d$.deme_lbl) & nzchar(vdf_d$.deme_lbl) & vdf_d$.deme_lbl != "NA", , drop = FALSE]
      d_tab <- sort(table(vdf_d$.deme_lbl), decreasing = TRUE)
      top_demes <- names(d_tab[d_tab >= 2L])
      if (length(top_demes) < 3L) top_demes <- names(d_tab)
      top_demes <- utils::head(top_demes, 12L)
      m_cols <- intersect(c("PC1", "PC2", "PC3", "ECT_PC1", "AR_inv", "solidity", "can1", "can2"), names(vdf_d))
      m_cols <- m_cols[vapply(m_cols, function(cl) sum(is.finite(suppressWarnings(as.numeric(vdf_d[[cl]])))) >= 4L, logical(1))]
      if (length(top_demes) >= 3L && length(m_cols) >= 1L) {
        vdf_d <- vdf_d[vdf_d$.deme_lbl %in% top_demes, , drop = FALSE]
        M_d <- do.call(cbind, lapply(m_cols, function(cl) {
          vv <- suppressWarnings(as.numeric(vdf_d[[cl]])); vv[!is.finite(vv)] <- stats::median(vv, na.rm = TRUE)
          s_v <- stats::sd(vv); if (!is.finite(s_v) || s_v <= 1e-6) s_v <- 1.0
          (vv - mean(vv)) / s_v
        }))
        cent_mat <- do.call(rbind, lapply(top_demes, function(dm) colMeans(M_d[vdf_d$.deme_lbl == dm, , drop = FALSE], na.rm = TRUE)))
        rownames(cent_mat) <- top_demes
        D_m <- as.matrix(stats::dist(cent_mat))
        sig_m <- max(stats::median(D_m[upper.tri(D_m)], na.rm = TRUE), 0.5)
        lon_c <- suppressWarnings(as.numeric(if ("decimalLongitude" %in% names(vdf_d)) vdf_d$decimalLongitude else vdf_d$longitude))
        lat_c <- suppressWarnings(as.numeric(if ("decimalLatitude" %in% names(vdf_d)) vdf_d$decimalLatitude else vdf_d$latitude))
        geo_cent <- do.call(rbind, lapply(top_demes, function(dm) {
          idx_dm <- which(vdf_d$.deme_lbl == dm)
          c(stats::median(lon_c[idx_dm], na.rm = TRUE), stats::median(lat_c[idx_dm], na.rm = TRUE))
        }))
        W_geo <- matrix(1.0, nrow = length(top_demes), ncol = length(top_demes))
        if (all(is.finite(geo_cent))) {
          k_dm <- nrow(geo_cent)
          D_km <- matrix(0.0, nrow = k_dm, ncol = k_dm)
          r_rad <- pi / 180.0
          for (ii in seq_len(k_dm)) {
            a_hv <- sin((geo_cent[, 2L] - geo_cent[ii, 2L]) * r_rad / 2.0)^2 +
              cos(geo_cent[ii, 2L] * r_rad) * cos(geo_cent[, 2L] * r_rad) *
              sin((geo_cent[, 1L] - geo_cent[ii, 1L]) * r_rad / 2.0)^2
            D_km[ii, ] <- 2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, a_hv))))
          }
          W_geo <- exp(-D_km / 600.0)
        }
        W <- exp(-(D_m^2) / (2.0 * sig_m^2)) * (0.5 + 0.5 * W_geo)
        dimnames(W) <- list(top_demes, top_demes)
        cdf <- data.frame(
          species_standardized = top_demes,
          consensus_complex_id = rep(1L, length(top_demes)),
          n_vouchers = as.integer(d_tab[top_demes]),
          silhouette_width = NA_real_,
          stringsAsFactors = FALSE
        )
        is_deme_net <- TRUE
      }
    }
  }
  if (nrow(cdf) == 0L || nrow(W) == 0L) return(invisible(""))
  taxa <- rownames(W)
  if (is.null(taxa)) taxa <- as.character(cdf$species_standardized)
  n_t <- length(taxa)

  w_vec <- suppressWarnings(as.numeric(W))
  w_vec[!is.finite(w_vec)] <- 0.0
  w_vec <- pmin(1.0, pmax(0.0, w_vec))
  W_num <- matrix(w_vec, nrow = n_t, ncol = n_t, dimnames = list(taxa, taxa))
  diag(W_num) <- 0.0

  ord_idx <- if (n_t >= 3L) {
    D_ward <- matrix(sqrt(pmax(1.0 - as.numeric(W_num), 0.0)), nrow = n_t, ncol = n_t, dimnames = list(taxa, taxa))
    diag(D_ward) <- 0.0
    tryCatch(stats::hclust(stats::as.dist(D_ward), method = "ward.D2")$order, error = function(e) seq_len(n_t))
  } else seq_len(n_t)
  taxa_ord <- taxa[ord_idx]
  W_disp <- W_num
  diag(W_disp) <- 1.0
  W_ord <- W_disp[ord_idx, ord_idx, drop = FALSE]

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
    ggplot2::scale_x_discrete(labels = if (is_deme_net) ggplot2::waiver() else icn_discrete_labeller(abbreviate_genus = TRUE)) +
    ggplot2::scale_y_discrete(labels = if (is_deme_net) ggplot2::waiver() else icn_discrete_labeller(abbreviate_genus = TRUE)) +
    theme_packera_pub(base_size = 9.4) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, face = "plain")) +
    ggplot2::labs(
      title = sprintf("A. Ward.D2 Multimodal Affinity Matrix (%d %s)", n_comp, if (n_comp == 1L) "Complex" else "Complexes"),
      subtitle = if (is_deme_net) {
        "Spatial-deme morphometric + geographic affinity (unimodal K* = 1)"
      } else if (is.finite(sil_mean)) {
        sprintf("Mean consensus silhouette width = %.3f", sil_mean)
      } else {
        "Morphometric + Ecological KDE + Circular Hopkins affinity"
      },
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

  idx_m <- match(taxa, as.character(cdf$species_standardized))
  comp_ids <- if ("consensus_complex_id" %in% names(cdf)) as.character(cdf$consensus_complex_id[idx_m]) else rep("1", n_t)
  comp_ids[is.na(comp_ids)] <- "1"
  comp_lbls <- paste0("Complex ", comp_ids)
  n_vouch <- if ("n_vouchers" %in% names(cdf)) suppressWarnings(as.numeric(cdf$n_vouchers[idx_m])) else rep(10.0, n_t)
  n_vouch[!is.finite(n_vouch) | n_vouch <= 0] <- 5.0

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
    lbl = if (is_deme_net) paste0(taxa, " (n=", n_vouch, ")") else abbreviate_icn_taxon(taxa, abbreviate_genus = TRUE),
    complex = if (is_deme_net) rep("Single Species (K* = 1)", n_t) else comp_lbls,
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
      for (es in list(
        list(d = edge_df[edge_df$edge_kind == "inter", , drop = FALSE], lty = "dotted", col = "#9CA3AF"),
        list(d = edge_df[edge_df$edge_kind == "intra", , drop = FALSE], lty = "solid", col = "#64748B")
      )) {
        if (nrow(es$d) > 0L) {
          pn <- pn + ggplot2::geom_segment(
            data = es$d, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
            linewidth = es$d$lw, linetype = es$lty, color = es$col, alpha = 0.35,
            show.legend = FALSE, inherit.aes = FALSE
          )
        }
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
  canon_07b <- file.path(fdir, "fig07b_taxon_complex_delimitation_network.png")
  cyto_sub_07b <- .build_cytology_fig_subtitle(if (!is.null(voucher_df) && is.data.frame(voucher_df) && nrow(voucher_df) > 0L) voucher_df else cdf)
  anno_07b <- patchwork::plot_annotation(
    title = "Step 07: Unsupervised Quantitative Taxon-Complex Delimitation",
    subtitle = if (nzchar(cyto_sub_07b)) cyto_sub_07b else NULL,
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12.4),
      plot.subtitle = ggplot2::element_text(size = 9.0, color = "#374151")
    )
  )
  plate_06b <- (p_heat | build_net_panel(TRUE)) +
    patchwork::plot_layout(widths = c(1.08, 1.0)) +
    anno_07b
  out_saved <- save_packera_png_300(plate_06b, canon_07b, width = 14.2, height = dyn_h_06b, also_copy_to = out_complex_path)
  if (!file.exists(canon_07b) || isTRUE(file.info(canon_07b)$size == 0)) {
    p_heat_fb <- if (is_deme_net) p_heat else (
      p_heat +
        ggplot2::scale_y_discrete(labels = function(x) abbreviate_icn_taxon(x, abbreviate_genus = TRUE)) +
        ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))
    )
    plate_06b_fb <- (p_heat_fb | build_net_panel(FALSE)) +
      patchwork::plot_layout(widths = c(1.08, 1.0)) +
      anno_07b
    out_saved <- save_packera_png_300(plate_06b_fb, canon_07b, width = 14.2, height = dyn_h_06b, also_copy_to = out_complex_path)
  }
  out_saved
}
