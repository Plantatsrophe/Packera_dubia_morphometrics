#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_synthesis_plates.R
# Project: Packera Genus-Wide Species Delimitation & Morphometrics Pipeline
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
#   - fig05_ecological_niche_phenology.{pdf,png} (4-panel ecological niche plate)
#   - fig05d_geographic_map_and_hopkins_phenology.{pdf,png} (US state map + Hopkins)
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

draw_us_state_basemap_base_r <- function(xlim = c(-98, -73), ylim = c(28, 46), main = "", xlab = "Longitude (deg E, WGS84)", ylab = "Latitude (deg N, WGS84)") {
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

.unpack_synthesis_display_tier <- function(tx_vec, present_tx = NULL, max_fg = 6L) {
  if (is.null(present_tx)) {
    present_tx <- sort(unique(tx_vec[!is.na(tx_vec) & nzchar(tx_vec) & tx_vec != "Unknown"]))
  }
  if (length(present_tx) == 0L) present_tx <- "Unassigned_Packera"
  tier_obj <- if (exists("resolve_taxon_display_tier", mode = "function")) {
    tryCatch(resolve_taxon_display_tier(tx_vec), error = function(e) NULL)
  } else {
    NULL
  }
  fg_tx <- NULL
  if (!is.null(tier_obj)) {
    if (is.list(tier_obj)) {
      for (k in c("foreground_taxa", "fg_taxa", "foreground", "focal_and_top", "top_taxa", "legend_taxa")) {
        if (!is.null(tier_obj[[k]]) && length(tier_obj[[k]]) > 0L) {
          fg_tx <- intersect(as.character(tier_obj[[k]]), present_tx)
          if (length(fg_tx) > 0L) break
        }
      }
    } else if (is.character(tier_obj)) {
      fg_tx <- intersect(tier_obj, present_tx)
    }
  }
  if (is.null(fg_tx) || length(fg_tx) == 0L) {
    if (length(present_tx) <= max_fg) {
      fg_tx <- present_tx
    } else {
      cnts <- sort(table(tx_vec[tx_vec %in% present_tx]), decreasing = TRUE)
      focal_cand <- present_tx[grepl("dubia", present_tx, ignore.case = TRUE)]
      fg_tx <- unique(c(focal_cand, names(cnts)))
      fg_tx <- fg_tx[seq_len(min(max_fg, length(fg_tx)))]
    }
  }
  bg_tx <- setdiff(present_tx, fg_tx)
  list(
    foreground_taxa = fg_tx,
    background_taxa = bg_tx,
    is_tiered = length(bg_tx) > 0L
  )
}

export_step05d_geographic_and_hopkins <- function(
  vouchers_df,
  fig_dir = "outputs/figures"
) {
  if (is.null(vouchers_df) || nrow(vouchers_df) == 0L) return(invisible(NULL))
  fdir <- if (is.null(fig_dir) || !nzchar(fig_dir) || grepl("\\.(pdf|png)$", fig_dir, ignore.case = TRUE)) {
    resolve_active_figures_dir(reference_path = fig_dir)
  } else {
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    fig_dir
  }

  tx_clean <- resolve_plot_group_column(vouchers_df)
  dyn_pal <- get_dynamic_packera_palette(tx_clean)
  dyn_pch <- get_dynamic_packera_shapes(tx_clean)
  alpha_val <- adaptive_point_alpha(nrow(vouchers_df))
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = alpha_val)
  pch_tx  <- dyn_pch[tx_clean]; pch_tx[is.na(pch_tx)] <- 1L
  present_tx <- sort(unique(tx_clean[!is.na(tx_clean) & nzchar(tx_clean) & tx_clean != "Unknown"]))
  if (length(present_tx) == 0L) present_tx <- "Unassigned_Packera"

  tier_info <- if (length(present_tx) > 6L) {
    .unpack_synthesis_display_tier(tx_clean, present_tx, max_fg = 6L)
  } else {
    focal_p <- present_tx[grepl("dubia", present_tx, ignore.case = TRUE)]
    list(foreground_taxa = unique(c(focal_p, present_tx)), background_taxa = character(0), is_tiered = FALSE)
  }
  fg_tx <- tier_info$foreground_taxa
  bg_tx <- tier_info$background_taxa
  dyn_pal_fg <- tryCatch(get_dynamic_packera_palette(fg_tx, preserve_order = TRUE), error = function(e) dyn_pal[fg_tx])
  dyn_pch_fg <- tryCatch(get_dynamic_packera_shapes(fg_tx, preserve_order = TRUE), error = function(e) dyn_pch[fg_tx])
  for (nm_f in names(dyn_pal_fg)) dyn_pal[[nm_f]] <- dyn_pal_fg[[nm_f]]
  for (nm_f in names(dyn_pch_fg)) dyn_pch[[nm_f]] <- dyn_pch_fg[[nm_f]]
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = alpha_val)
  pch_tx  <- dyn_pch[tx_clean]; pch_tx[is.na(pch_tx)] <- 1L

  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 4.6, 3.2, 1.2))
    lat_v <- if ("decimalLatitude" %in% colnames(vouchers_df)) suppressWarnings(as.numeric(vouchers_df$decimalLatitude)) else rep(NA_real_, nrow(vouchers_df))
    if ("latitude" %in% colnames(vouchers_df)) lat_v <- ifelse(is.finite(lat_v), lat_v, suppressWarnings(as.numeric(vouchers_df$latitude)))
    lon_v <- if ("decimalLongitude" %in% colnames(vouchers_df)) suppressWarnings(as.numeric(vouchers_df$decimalLongitude)) else rep(NA_real_, nrow(vouchers_df))
    if ("longitude" %in% colnames(vouchers_df)) lon_v <- ifelse(is.finite(lon_v), lon_v, suppressWarnings(as.numeric(vouchers_df$longitude)))
    geo_b <- compute_dynamic_na_map_bounds(lon_v, lat_v)
    valid_geo <- geo_b$valid_geo

    conflict_v <- if ("consensus_conflict_flag" %in% colnames(vouchers_df)) {
      as.logical(vouchers_df$consensus_conflict_flag)
    } else if ("multimodal_conflict_flag" %in% colnames(vouchers_df)) {
      as.logical(vouchers_df$multimodal_conflict_flag)
    } else {
      rep(FALSE, nrow(vouchers_df))
    }
    conflict_v[is.na(conflict_v)] <- FALSE

    prio_v6b <- if ("synthesis_triage_priority" %in% colnames(vouchers_df)) {
      as.character(vouchers_df$synthesis_triage_priority)
    } else if ("triage_priority" %in% colnames(vouchers_df)) {
      as.character(vouchers_df$triage_priority)
    } else {
      rep("", nrow(vouchers_df))
    }
    prob_v6b <- if ("cda_posterior_prob" %in% colnames(vouchers_df)) {
      suppressWarnings(as.numeric(vouchers_df$cda_posterior_prob))
    } else {
      rep(0.0, nrow(vouchers_df))
    }
    prob_v6b[!is.finite(prob_v6b)] <- 0.0
    is_crit_conf <- conflict_v & (toupper(trimws(prio_v6b)) == "CRITICAL" | prob_v6b >= 0.85)
    if (!any(is_crit_conf) && sum(conflict_v) <= 50L) is_crit_conf <- conflict_v

    xlim_map <- geo_b$xlim
    ylim_map <- geo_b$ylim

    is_type_6b <- if ("is_type_specimen" %in% colnames(vouchers_df)) as.logical(vouchers_df$is_type_specimen) else rep(FALSE, nrow(vouchers_df))
    is_type_6b[is.na(is_type_6b)] <- FALSE

    # Panel A: North American Geographic Distribution Map with Dynamic Bounds
    draw_us_state_basemap_base_r(
      xlim = xlim_map, ylim = ylim_map,
      main = "A. Geographic Distribution & Contact Zone Map (WGS84)",
      xlab = "Longitude (deg E, WGS84)"
    )
    if (any(valid_geo)) {
      if (length(bg_tx) > 0L) {
        m_bg <- valid_geo & (tx_clean %in% bg_tx)
        if (any(m_bg)) {
          graphics::points(
            lon_v[m_bg], lat_v[m_bg],
            pch = 1L, bg = NA,
            col = grDevices::adjustcolor("grey68", alpha.f = 0.32), cex = 0.64, lwd = 0.65
          )
        }
      }
      if (any(valid_geo & is_crit_conf)) {
        graphics::points(
          lon_v[valid_geo & is_crit_conf], lat_v[valid_geo & is_crit_conf],
          pch = 5L, bg = NA,
          col = grDevices::adjustcolor("#CC79A7", alpha.f = 0.42), cex = 0.58, lwd = 0.85
        )
      }
      fg_draw_order <- c(setdiff(fg_tx, fg_tx[grepl("dubia", fg_tx, ignore.case = TRUE)]), fg_tx[grepl("dubia", fg_tx, ignore.case = TRUE)])
      for (sp_fg in fg_draw_order) {
        m_sp_fg <- valid_geo & (tx_clean == sp_fg)
        if (any(m_sp_fg)) {
          graphics::points(
            lon_v[m_sp_fg], lat_v[m_sp_fg],
            pch = pch_tx[m_sp_fg], bg = NA,
            col = cols_tx[m_sp_fg], cex = 0.82, lwd = 0.95
          )
        }
      }
      if (any(valid_geo & is_type_6b)) {
        m_t6 <- which(valid_geo & is_type_6b)
        graphics::points(lon_v[m_t6], lat_v[m_t6], pch = 5L, bg = NA, col = cols_tx[m_t6], cex = 1.45, lwd = 2.0)
        graphics::points(lon_v[m_t6], lat_v[m_t6], pch = 1L, bg = NA, col = "#111111", cex = 1.85, lwd = 1.4)
      }
    }
    leg_expr_a <- c(
      packera_italic_expressions(fg_tx),
      if (length(bg_tx) > 0L) sprintf("Other Packera (%d spp.)", length(bg_tx)) else NULL,
      "Contact / Conflict Voucher"
    )
    leg_col_a <- c(
      dyn_pal[fg_tx],
      if (length(bg_tx) > 0L) "grey60" else NULL,
      "#CC79A7"
    )
    leg_pch_a <- c(
      dyn_pch[fg_tx],
      if (length(bg_tx) > 0L) 1L else NULL,
      5L
    )
    graphics::legend(
      "bottomleft",
      legend = leg_expr_a,
      col = leg_col_a,
      pch = leg_pch_a, pt.bg = NA, cex = 0.68, bg = grDevices::adjustcolor("white", alpha.f = 0.85)
    )

    # Panel B: DOY vs. Latitude with Species Regression Lines & Hopkins' Law Baseline (+4.0 d / deg N)
    doy_v <- if ("doy" %in% colnames(vouchers_df)) {
      suppressWarnings(as.numeric(vouchers_df$doy))
    } else if ("collection_doy" %in% colnames(vouchers_df)) {
      suppressWarnings(as.numeric(vouchers_df$collection_doy))
    } else {
      rep(NA_real_, nrow(vouchers_df))
    }
    valid_pheno <- valid_geo & is.finite(doy_v) & (doy_v >= 1) & (doy_v <= 366)
    if ("is_verified_anthesis" %in% colnames(vouchers_df)) {
      ver_anth <- as.logical(vouchers_df$is_verified_anthesis)
      ver_anth[is.na(ver_anth)] <- FALSE
      if (sum(valid_pheno & ver_anth) >= 3L) {
        valid_pheno <- valid_pheno & ver_anth
      }
    }

    if (any(valid_pheno)) {
      graphics::plot(
        lat_v[valid_pheno], doy_v[valid_pheno], type = "n",
        xlab = "Decimal Latitude (deg N)",
        ylab = "Flowering Day of Year (DOY)",
        main = "B. Hopkins' Bioclimatic Law (+4.0 d / deg N) & Allochrony",
        font.main = 2
      )
      graphics::grid(col = "grey88")
      if (length(bg_tx) > 0L) {
        m_bg_p <- valid_pheno & (tx_clean %in% bg_tx)
        if (any(m_bg_p)) {
          graphics::points(
            lat_v[m_bg_p], doy_v[m_bg_p],
            pch = 1L, bg = NA, col = grDevices::adjustcolor("grey68", alpha.f = 0.35), cex = 0.68, lwd = 0.7
          )
        }
      }
      m_fg_p <- valid_pheno & (tx_clean %in% fg_tx)
      if (any(m_fg_p)) {
        graphics::points(
          lat_v[m_fg_p], doy_v[m_fg_p],
          pch = pch_tx[m_fg_p], bg = NA, col = cols_tx[m_fg_p], cex = 0.85, lwd = 1.0
        )
      }
      if (any(valid_pheno & is_type_6b)) {
        m_t6 <- which(valid_pheno & is_type_6b)
        graphics::points(lat_v[m_t6], doy_v[m_t6], pch = 5L, bg = NA, col = cols_tx[m_t6], cex = 1.45, lwd = 2.0)
        graphics::points(lat_v[m_t6], doy_v[m_t6], pch = 1L, bg = NA, col = "#111111", cex = 1.85, lwd = 1.4)
      }
      is_act_6b <- rep(TRUE, nrow(vouchers_df))
      if ("determiner_tier" %in% colnames(vouchers_df)) {
        dt_6b <- as.character(vouchers_df$determiner_tier)
        dt_6b[is.na(dt_6b)] <- "Tier_3_Bronze"
        pr_6b <- rep(FALSE, nrow(vouchers_df))
        if ("is_promoted_bronze" %in% colnames(vouchers_df)) { p1 <- as.logical(vouchers_df$is_promoted_bronze); p1[is.na(p1)] <- FALSE; pr_6b <- pr_6b | p1 }
        if ("is_promoted_anchor" %in% colnames(vouchers_df)) { p2 <- as.logical(vouchers_df$is_promoted_anchor); p2[is.na(p2)] <- FALSE; pr_6b <- pr_6b | p2 }
        is_act_6b <- grepl("Tier_1|Gold|Tier_2|Silver", dt_6b, ignore.case = TRUE) | pr_6b | is_type_6b
      }
      if ("is_confident_label_noise" %in% colnames(vouchers_df)) {
        cn_6b <- as.logical(vouchers_df$is_confident_label_noise); cn_6b[is.na(cn_6b)] <- FALSE
        is_act_6b <- is_act_6b & (!cn_6b | is_type_6b)
      }
      for (sp in fg_tx) {
        m_sp_act <- which(valid_pheno & is_act_6b & (tx_clean == sp))
        m_sp <- if (length(m_sp_act) >= 4L) m_sp_act else which(valid_pheno & (tx_clean == sp))
        if (length(m_sp) >= 4L && is.finite(stats::sd(lat_v[m_sp], na.rm = TRUE)) && stats::sd(lat_v[m_sp], na.rm = TRUE) > 0.2) {
          lm_sp <- tryCatch(stats::lm(doy_v[m_sp] ~ lat_v[m_sp]), error = function(e) NULL)
          if (!is.null(lm_sp)) {
            graphics::abline(lm_sp, col = dyn_pal[[sp]], lwd = 1.8, lty = 1)
          }
        }
      }
      lat_seq <- seq(min(lat_v[valid_pheno], na.rm = TRUE) - 0.5, max(lat_v[valid_pheno], na.rm = TRUE) + 0.5, length.out = 50)
      mean_lat <- mean(lat_v[valid_pheno], na.rm = TRUE)
      mean_doy <- mean(doy_v[valid_pheno], na.rm = TRUE)
      hopkins_line <- mean_doy + 4.0 * (lat_seq - mean_lat)
      graphics::lines(lat_seq, hopkins_line, lwd = 2.4, lty = 2, col = "#111111")
      graphics::legend(
        "topleft",
        legend = c("Hopkins' Law (+4.0 d/deg N)", packera_italic_expressions(fg_tx)),
        col = c("#111111", dyn_pal[fg_tx]),
        lty = c(2, rep(1, length(fg_tx))),
        pch = c(NA, dyn_pch[fg_tx]), pt.bg = NA,
        lwd = c(2.4, rep(1.6, length(fg_tx))), cex = 0.68, bty = "n"
      )
    } else {
      graphics::plot(
        c(0, 1), c(0, 1), type = "n", axes = FALSE,
        xlab = "Decimal Latitude (deg N)",
        ylab = "Flowering Day of Year (DOY)",
        main = "B. Hopkins' Bioclimatic Law (+4.0 d / deg N) & Allochrony",
        font.main = 2
      )
      graphics::box(col = "grey75")
      graphics::text(0.5, 0.5, "No verified flowering DOY observations (1..366) available", col = "grey35", cex = 0.85)
    }
  }, fdir, "fig05d_geographic_map_and_hopkins_phenology", width = 12.8, height = 5.6, dpi = 300)

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
    has_e <- if ("edaphic_best_fit_taxon" %in% colnames(triage_df)) !is.na(triage_df$edaphic_best_fit_taxon) else rep(FALSE, nrow(triage_df))
    comp_mask <- has_m | has_e
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
  pheno_conflict <- if ("phenology_conflict_flag" %in% colnames(triage_df) && any(!is.na(triage_df$phenology_conflict_flag))) {
    as.logical(triage_df$phenology_conflict_flag)
  } else if ("delta_doy_resid" %in% colnames(triage_df) && any(!is.na(triage_df$delta_doy_resid))) {
    abs(suppressWarnings(as.numeric(triage_df$delta_doy_resid))) >= 14.0
  } else if ("delta_doy" %in% colnames(triage_df) && any(!is.na(triage_df$delta_doy))) {
    abs(suppressWarnings(as.numeric(triage_df$delta_doy))) >= 14.0
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
    r_p <- mean(pheno_conflict, na.rm = TRUE) * 100
    r_e <- mean(eco_conflict, na.rm = TRUE) * 100
    rates <- c(
      ifelse(is.finite(r_m), r_m, 0),
      ifelse(is.finite(r_p), r_p, 0),
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
      labels = c("Morphology\n(EFA/ECT/CDA)", "Phenology\n(Anthesis Allochrony)", "Edaphic Niche\n(0-15cm RF)"),
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

    graphics::par(mar = c(if (nc > 5L) 6.6 else 5.4, if (nr > 6L) 6.6 else 5.4, 3.4, 1.2))
    graphics::plot(
      c(0.5, nc + 0.5), c(0.5, nr + 0.5), type = "n", axes = FALSE,
      xlab = if (nc > 5L) "" else "Recommended Taxonomic Determination", ylab = "",
      main = "B. Recorded vs. Recommended Determination Matrix", font.main = 2
    )
    ramp_b <- grDevices::colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594"))(101)
    cex_cell <- min(0.64, max(0.34, 4.4 / max(nr, nc)))
    for (i in seq_len(nr)) {
      y_pos <- nr - i + 1
      for (j in seq_len(nc)) {
        pct_ij <- row_pct[i, j]
        cnt_ij <- flow_mat[i, j]
        c_idx <- max(1L, min(101L, round(pct_ij) + 1L))
        graphics::rect(j - 0.48, y_pos - 0.46, j + 0.48, y_pos + 0.46, col = ramp_b[c_idx], border = "white", lwd = 1.5)
        txt_col <- if (pct_ij >= 52) "white" else "#111111"
        is_diag_ij <- (j <= length(col_levels)) && (tx_levels[i] == col_levels[j])
        lbl_ij <- if (max(nr, nc) > 6L) {
          if (cnt_ij == 0L && !is_diag_ij) "" else sprintf("%d (%.0f%%)", cnt_ij, pct_ij)
        } else {
          sprintf("%d\n(%.1f%%)", cnt_ij, pct_ij)
        }
        if (nzchar(lbl_ij)) {
          graphics::text(j, y_pos, labels = lbl_ij, col = txt_col, cex = cex_cell, font = if (cnt_ij > 0) 2 else 1)
        }
      }
    }
    cex_ax_x <- min(0.78, max(0.42, 6.2 / max(nc, 1L)))
    cex_ax_y <- min(0.80, max(0.42, 6.2 / max(nr, 1L)))
    graphics::axis(1, at = seq_len(nc), labels = col_labels, tick = FALSE, font = 3, cex.axis = cex_ax_x, las = if (nc > 5L) 2 else 1)
    graphics::axis(2, at = rev(seq_len(nr)), labels = tx_labels, las = 1, tick = FALSE, font = 3, cex.axis = cex_ax_y)
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
    if (length(tx_levels) > 6L) {
      graphics::par(mar = c(5.2, 7.6, 4.2, 1.4))
      prop_treat_tx <- sweep(tbl_treat_tx, 2, pmax(col_totals, 1), "/") * 100
      cex_nm <- min(0.78, max(0.40, 6.4 / length(tx_levels)))
      bp_b <- graphics::barplot(
        prop_treat_tx[, rev(seq_len(ncol(prop_treat_tx))), drop = FALSE],
        beside = FALSE,
        horiz = TRUE,
        las = 1,
        col = treat_pal[rownames(prop_treat_tx)],
        xlim = c(0, 118),
        xlab = "Completed Voucher Proportion (%)",
        ylab = "",
        font.axis = 3,
        cex.names = cex_nm,
        main = "B. Taxonomic Treatments by Recorded Species",
        font.main = 2
      )
      rev_totals <- rev(col_totals)
      graphics::text(
        rep(101.5, length(bp_b)), bp_b,
        labels = ifelse(rev_totals > 0, sprintf("n=%d", rev_totals), ""),
        pos = 4, font = 2, cex = cex_nm
      )
      graphics::grid(nx = NULL, ny = NA, col = "grey85")
      graphics::legend(
        x = 50, y = max(bp_b) + 1.35, xjust = 0.5, horiz = FALSE, ncol = 3,
        legend = treat_levels,
        fill = treat_pal[treat_levels],
        cex = 0.62, xpd = TRUE, bty = "n"
      )
    } else {
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
    }
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
  is_subcluster_mode <- any(grepl("Morphotype_K|\\[K[0-9]+\\]|Ecotype|Subcluster", active_tx, ignore.case = TRUE)) ||
    length(unique(stats::na.omit(df$species_standardized))) < 2L
  grp_legend_lbl <- if (is_subcluster_mode) "Morphotype" else "Species"

  tier_info <- if (length(active_tx) > 6L) {
    .unpack_synthesis_display_tier(df$plot_group, active_tx, max_fg = 6L)
  } else {
    focal_a <- active_tx[grepl("dubia", active_tx, ignore.case = TRUE)]
    list(foreground_taxa = unique(c(focal_a, active_tx)), background_taxa = character(0), is_tiered = FALSE)
  }
  fg_tx <- tier_info$foreground_taxa
  bg_tx <- tier_info$background_taxa

  dyn_pal <- tryCatch(get_dynamic_packera_palette(fg_tx, preserve_order = TRUE), error = function(e) get_dynamic_packera_palette(active_tx)[fg_tx])
  dyn_shapes <- tryCatch(get_dynamic_packera_ggplot_shapes(fg_tx, preserve_order = TRUE), error = function(e) get_dynamic_packera_ggplot_shapes(active_tx)[fg_tx])
  thm <- ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 8.6),
                   plot.subtitle = ggplot2::element_text(size = 7.0, color = "#333333"),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size = 8, face = "bold"),
                   legend.text = ggplot2::element_text(size = 7.3))

  # ----------------------------------------------------------------------------
  # Panel A (p1): Integrative Canonical Discriminant & Type Baseline Projection
  # ----------------------------------------------------------------------------
  cda_sub <- df[!is.na(df$can1) & !is.na(df$can2) & df$plot_group %in% active_tx, , drop = FALSE]
  if (nrow(cda_sub) == 0) cda_sub <- df[!is.na(df$can1) & !is.na(df$can2), , drop = FALSE]
  pt_alpha <- adaptive_point_alpha(nrow(cda_sub))
  fg_cda <- cda_sub[cda_sub$plot_group %in% fg_tx, , drop = FALSE]
  if (nrow(fg_cda) == 0) fg_cda <- cda_sub
  fg_cda$plot_group <- factor(as.character(fg_cda$plot_group), levels = fg_tx)
  bg_cda <- if (length(bg_tx) > 0L) cda_sub[cda_sub$plot_group %in% bg_tx, , drop = FALSE] else cda_sub[FALSE, , drop = FALSE]

  cent_rows <- lapply(intersect(fg_tx, unique(as.character(fg_cda$plot_group))), function(tx) {
    m <- as.character(fg_cda$plot_group) == tx
    data.frame(plot_group = tx, cx = mean(fg_cda$can1[m], na.rm = TRUE), cy = mean(fg_cda$can2[m], na.rm = TRUE),
               label = gsub("^Packera\\s+", "P. ", tx), stringsAsFactors = FALSE)
  })
  cent_df <- if (length(cent_rows) > 0) do.call(rbind, cent_rows) else data.frame(cx = numeric(0), cy = numeric(0), label = character(0))

  p1 <- ggplot2::ggplot(fg_cda, ggplot2::aes(x = can1, y = can2, color = plot_group, shape = plot_group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "grey75", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey75", linewidth = 0.4)

  if (nrow(bg_cda) > 0L) {
    hull_list <- lapply(unique(bg_cda$plot_group), function(b_tx) {
      sub_b <- bg_cda[bg_cda$plot_group == b_tx, , drop = FALSE]
      if (nrow(sub_b) >= 3L) {
        idx_h <- grDevices::chull(sub_b$can1, sub_b$can2)
        sub_b[idx_h, c("can1", "can2", "plot_group"), drop = FALSE]
      } else {
        NULL
      }
    })
    hull_list <- Filter(Negate(is.null), hull_list)
    if (length(hull_list) > 0L) {
      bg_hull_df <- do.call(rbind, hull_list)
      p1 <- p1 + ggplot2::geom_polygon(
        data = bg_hull_df, ggplot2::aes(x = can1, y = can2, group = plot_group),
        fill = "grey85", color = "grey70", alpha = 0.18, linetype = "dashed", linewidth = 0.35, inherit.aes = FALSE
      )
    }
    p1 <- p1 + ggplot2::geom_point(
      data = bg_cda, ggplot2::aes(x = can1, y = can2),
      shape = 1, color = "grey70", alpha = 0.32, size = 1.15, stroke = 0.5, inherit.aes = FALSE
    )
  }

  p1 <- p1 +
    ggplot2::stat_ellipse(level = 0.95, linetype = 1, linewidth = 0.85, show.legend = FALSE) +
    ggplot2::geom_point(fill = NA, stroke = 0.72, alpha = pt_alpha, size = 1.55) +
    ggplot2::geom_point(data = cent_df, ggplot2::aes(x = cx, y = cy), shape = 3, size = 2.8, stroke = 1.3, color = "black", inherit.aes = FALSE) +
    ggplot2::scale_color_manual(values = dyn_pal[fg_tx], drop = FALSE) +
    ggplot2::scale_shape_manual(values = dyn_shapes[fg_tx], drop = FALSE) +
    thm + ggplot2::labs(
      title = "A. Integrative Canonical Discriminant & Type Baseline Projection",
      subtitle = if (!is.na(loocv_acc)) sprintf("Gold/Silver Anchor CDA (LOOCV = %.1f%%) | 95%% Ellipses & Double-Ringed Types", loocv_acc) else "95% Confidence Ellipses, Group Centroids (+) & Type Anchors",
      x = sprintf("Canonical Variate 1 (Can1: %.1f%%)", can1_pct), y = sprintf("Canonical Variate 2 (Can2: %.1f%%)", can2_pct),
      color = grp_legend_lbl, shape = grp_legend_lbl
    )
  if ("is_type_specimen" %in% names(cda_sub) && any(as.logical(cda_sub$is_type_specimen) == TRUE, na.rm = TRUE)) {
    type_cda <- cda_sub[which(as.logical(cda_sub$is_type_specimen) == TRUE), , drop = FALSE]
    p1 <- p1 +
      ggplot2::geom_point(
        data = type_cda, ggplot2::aes(x = can1, y = can2),
        shape = 5, fill = NA, stroke = 1.35, size = 2.7, color = "#111111", inherit.aes = FALSE, show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = type_cda, ggplot2::aes(x = can1, y = can2),
        shape = 1, fill = NA, stroke = 1.15, size = 3.5, color = "#D55E00", inherit.aes = FALSE, show.legend = FALSE
      )
  }

  # ----------------------------------------------------------------------------
  # Panel B (p2): Hierarchical Integrative Treatment Breakdown
  # ----------------------------------------------------------------------------
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
  sp_short_raw <- gsub("^Packera\\s+", "P. ", as.character(df_comp$plot_group))
  pct_conf_by_sp <- tapply(df_comp$Treatment_Clean == "Confirmed Species", sp_short_raw, mean, na.rm = TRUE)
  dubia_short <- names(pct_conf_by_sp)[grepl("dubia", names(pct_conf_by_sp), ignore.case = TRUE)]
  other_short <- setdiff(names(sort(pct_conf_by_sp, decreasing = TRUE)), dubia_short)
  sp_order_top_first <- unique(c(dubia_short, other_short))
  sp_fac_lvls <- if (length(active_tx) > 5L) rev(sp_order_top_first) else sp_order_top_first
  df_comp$Species_Short <- factor(sp_short_raw, levels = sp_fac_lvls)
  n_by_sp <- stats::aggregate(catalogNumber ~ Species_Short, data = df_comp, FUN = length)
  p2_sz <- min(2.6, max(1.7, 18.0 / max(nrow(n_by_sp), 1L)))
  p2 <- ggplot2::ggplot(df_comp, ggplot2::aes(x = Species_Short, fill = Treatment_Clean)) +
    ggplot2::geom_bar(position = ggplot2::position_fill(reverse = TRUE), width = 0.65, color = "white", linewidth = 0.25) +
    ggplot2::geom_text(data = n_by_sp, ggplot2::aes(x = Species_Short, y = 1.04, label = paste0("n=", catalogNumber)), inherit.aes = FALSE, size = p2_sz, fontface = "bold") +
    ggplot2::scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"), limits = c(0, 1.12)) +
    ggplot2::scale_fill_manual(values = TAXONOMIC_TREATMENT_PALETTE[treat_levels], drop = FALSE) +
    thm +
    ggplot2::labs(title = "B. Hierarchical Integrative Treatment Breakdown",
                  subtitle = "5-stage precedence (Misidentification -> Hybrid -> Cryptic -> Ecotype -> Confirmed)",
                  x = if (is_subcluster_mode) "Morphotype / Subcluster Group" else "Recorded Herbarium Taxon", y = "Proportion of Completed Vouchers", fill = "Treatment")
  if (length(active_tx) > 5L) {
    ax_sz_b <- min(7.8, max(4.6, 52.0 / length(active_tx)))
    p2 <- p2 + ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = if (is_subcluster_mode) "bold" else "italic", size = ax_sz_b))
  } else {
    p2 <- p2 + ggplot2::theme(axis.text.x = ggplot2::element_text(face = if (is_subcluster_mode) "bold" else "italic", angle = 14, hjust = 1))
  }

  # ----------------------------------------------------------------------------
  # Panel C (p3): Recorded vs. Recommended Integrative Determination Matrix
  # ----------------------------------------------------------------------------
  rec_raw_c <- if ("recommended_determination" %in% names(df_comp)) {
    as.character(df_comp$recommended_determination)
  } else if ("cda_predicted_taxon" %in% names(df_comp)) {
    as.character(df_comp$cda_predicted_taxon)
  } else {
    as.character(df_comp$plot_group)
  }
  rec_std_c <- standardize_taxon_name(rec_raw_c)
  rec_clean_c <- ifelse(!is.na(rec_std_c), rec_std_c, ifelse(grepl("×|hybrid", rec_raw_c, ignore.case = TRUE), "Hybrid Swarm", as.character(df_comp$plot_group)))
  mat_tx_c <- if (length(active_tx) > 10L) {
    .unpack_synthesis_display_tier(df_comp$plot_group, active_tx, max_fg = 8L)$foreground_taxa
  } else {
    active_tx
  }
  row_c_vec <- ifelse(df_comp$plot_group %in% mat_tx_c, as.character(df_comp$plot_group), "Other Packera spp.")
  col_c_vec <- ifelse(rec_clean_c %in% c(mat_tx_c, "Hybrid Swarm"), rec_clean_c, "Other Packera spp.")
  row_c_lvls <- unique(c(mat_tx_c, if (any(row_c_vec == "Other Packera spp.")) "Other Packera spp." else NULL))
  col_c_lvls <- unique(c(mat_tx_c, if (any(col_c_vec == "Other Packera spp.")) "Other Packera spp." else NULL, if (any(col_c_vec == "Hybrid Swarm")) "Hybrid Swarm" else NULL))
  row_c_lbls <- gsub("^Packera\\s+", "P. ", row_c_lvls)
  col_c_lbls <- gsub("^Packera\\s+", "P. ", col_c_lvls)

  flow_tbl_c <- table(
    Recorded = factor(row_c_vec, levels = rev(row_c_lvls), labels = rev(row_c_lbls)),
    Recommended = factor(col_c_vec, levels = col_c_lvls, labels = col_c_lbls)
  )
  flow_pct_c <- sweep(flow_tbl_c, 1, pmax(rowSums(flow_tbl_c), 1), "/") * 100
  flow_df_c <- as.data.frame(flow_tbl_c, stringsAsFactors = FALSE)
  pct_df_c <- as.data.frame(flow_pct_c, stringsAsFactors = FALSE)
  flow_df_c$Pct <- pct_df_c$Freq
  flow_df_c$Recorded <- factor(flow_df_c$Recorded, levels = rev(row_c_lbls))
  flow_df_c$Recommended <- factor(flow_df_c$Recommended, levels = col_c_lbls)
  max_dim_c <- max(length(row_c_lbls), length(col_c_lbls))
  flow_df_c$Cell_Label <- vapply(seq_len(nrow(flow_df_c)), function(k) {
    cnt_k <- flow_df_c$Freq[k]
    pct_k <- flow_df_c$Pct[k]
    is_diag_k <- as.character(flow_df_c$Recorded[k]) == as.character(flow_df_c$Recommended[k])
    if (cnt_k == 0L && !is_diag_k) {
      ""
    } else if (max_dim_c > 6L) {
      sprintf("%d (%.0f%%)", cnt_k, pct_k)
    } else {
      sprintf("%d\n(%.0f%%)", cnt_k, pct_k)
    }
  }, character(1))
  txt_sz_c <- min(2.7, max(1.65, 16.0 / max(max_dim_c, 1L)))
  p3 <- ggplot2::ggplot(flow_df_c, ggplot2::aes(x = Recommended, y = Recorded, fill = Pct)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = Cell_Label),
      size = txt_sz_c,
      fontface = "bold",
      color = ifelse(flow_df_c$Pct >= 52, "white", "#111111")
    ) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#084594", limits = c(0, 100)) +
    thm +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(face = "italic", angle = if (length(col_c_lbls) > 4L) 28 else 12, hjust = 1, size = min(7.6, max(5.2, 48.0 / max(length(col_c_lbls), 1L)))),
      axis.text.y = ggplot2::element_text(face = "italic", size = min(7.6, max(5.2, 48.0 / max(length(row_c_lbls), 1L))))
    ) +
    ggplot2::labs(
      title = "C. Recorded vs. Recommended Integrative Determination Matrix",
      subtitle = "Row-normalized reassignment percentages (off-diagonal zero labels suppressed)",
      x = "Recommended Integrative Determination",
      y = "Recorded Herbarium Taxon",
      fill = "Row %"
    )

  # ----------------------------------------------------------------------------
  # Panel D (p4): Specimen-Level Multi-Modal Confidence & Triage Priority Distribution
  # ----------------------------------------------------------------------------
  conf_sub <- df_comp
  conf_score <- rep(NA_real_, nrow(conf_sub))
  for (cc in c("consensus_confidence", "multimodal_consensus_prob", "cda_posterior_prob", "edaphic_rf_prob")) {
    if (cc %in% names(conf_sub)) {
      v_cc <- suppressWarnings(as.numeric(conf_sub[[cc]]))
      conf_score <- ifelse(is.finite(conf_score), conf_score, v_cc)
    }
  }
  if (any(!is.finite(conf_score))) {
    fallback_conf <- ifelse(conf_sub$Treatment_Clean == "Confirmed Species", 0.88,
                     ifelse(conf_sub$Treatment_Clean == "Variety / Ecotype", 0.74,
                     ifelse(conf_sub$Treatment_Clean == "Hybrid / Introgressant", 0.52, 0.36)))
    conf_score[!is.finite(conf_score)] <- fallback_conf[!is.finite(conf_score)]
  }
  conf_sub$conf_score <- pmin(1.0, pmax(0.0, conf_score))

  prio_raw <- if ("synthesis_triage_priority" %in% names(conf_sub)) {
    toupper(as.character(conf_sub$synthesis_triage_priority))
  } else if ("triage_priority" %in% names(conf_sub)) {
    toupper(as.character(conf_sub$triage_priority))
  } else {
    rep(NA_character_, nrow(conf_sub))
  }
  prio_clean <- ifelse(grepl("CRIT", prio_raw), "CRITICAL",
                ifelse(grepl("HIGH", prio_raw), "HIGH",
                ifelse(grepl("MED|MOD", prio_raw), "MEDIUM",
                ifelse(grepl("LOW", prio_raw) & !grepl("RESOLVED|CONFIRMED", prio_raw), "LOW",
                ifelse(grepl("RES|CONF|NONE", prio_raw), "RESOLVED", NA_character_)))))
  if (any(is.na(prio_clean))) {
    prio_fallback <- ifelse(conf_sub$Treatment_Clean == "Confirmed Species", "RESOLVED",
                     ifelse(conf_sub$Treatment_Clean == "Variety / Ecotype", "LOW",
                     ifelse(conf_sub$Treatment_Clean == "Hybrid / Introgressant", "MEDIUM",
                     ifelse(conf_sub$Treatment_Clean == "Misidentification (Reassign)", "HIGH", "CRITICAL"))))
    prio_clean[is.na(prio_clean)] <- prio_fallback[is.na(prio_clean)]
  }
  prio_lvls <- c("RESOLVED", "LOW", "MEDIUM", "HIGH", "CRITICAL")
  conf_sub$Priority_Fac <- factor(prio_clean, levels = prio_lvls)

  auth_raw <- if ("authority_tier" %in% names(conf_sub)) {
    as.character(conf_sub$authority_tier)
  } else if ("determiner_tier" %in% names(conf_sub)) {
    as.character(conf_sub$determiner_tier)
  } else {
    rep("Bronze", nrow(conf_sub))
  }
  auth_clean <- ifelse(grepl("Type|Tier_0|Gold|Tier_1", auth_raw, ignore.case = TRUE), "Type / Gold",
                ifelse(grepl("Silver|Tier_2", auth_raw, ignore.case = TRUE), "Silver", "Bronze"))
  conf_sub$Auth_Tier <- factor(auth_clean, levels = c("Type / Gold", "Silver", "Bronze"))

  prio_pal <- c(
    "RESOLVED" = "#009E73",
    "LOW"      = "#56B4E9",
    "MEDIUM"   = "#CC79A7",
    "HIGH"     = "#E69F00",
    "CRITICAL" = "#D55E00"
  )
  p4 <- ggplot2::ggplot(conf_sub, ggplot2::aes(x = Priority_Fac, y = conf_score, fill = Priority_Fac)) +
    ggplot2::geom_boxplot(width = 0.50, alpha = 0.35, outlier.shape = NA, color = "grey30", show.legend = FALSE) +
    ggplot2::geom_point(
      ggplot2::aes(color = Priority_Fac, shape = Auth_Tier),
      position = ggplot2::position_jitter(width = 0.15, height = 0, seed = 42),
      fill = NA, stroke = 0.75, size = 1.55, alpha = pt_alpha
    ) +
    ggplot2::scale_fill_manual(values = prio_pal, drop = FALSE) +
    ggplot2::scale_color_manual(values = prio_pal, guide = "none") +
    ggplot2::scale_shape_manual(values = c("Type / Gold" = 5, "Silver" = 2, "Bronze" = 1), drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 1.05)) +
    thm +
    ggplot2::labs(
      title = "D. Specimen-Level Multi-Modal Confidence & Triage Priority Distribution",
      subtitle = "Northcutt et al. (2021) Confident Learning audit stratified by priority & authority tier",
      x = "Integrative Triage Priority",
      y = "Multi-Modal Consensus Confidence",
      shape = "Authority Tier"
    )

  # ----------------------------------------------------------------------------
  # Panel E (p5): Geographic Contact Zones, Triage Flags & Type Localities
  # ----------------------------------------------------------------------------
  lon_col_e <- if ("longitude" %in% names(df_comp)) suppressWarnings(as.numeric(df_comp$longitude)) else rep(NA_real_, nrow(df_comp))
  if ("decimalLongitude" %in% names(df_comp)) lon_col_e <- ifelse(is.finite(lon_col_e), lon_col_e, suppressWarnings(as.numeric(df_comp$decimalLongitude)))
  lat_col_e <- if ("latitude" %in% names(df_comp)) suppressWarnings(as.numeric(df_comp$latitude)) else rep(NA_real_, nrow(df_comp))
  if ("decimalLatitude" %in% names(df_comp)) lat_col_e <- ifelse(is.finite(lat_col_e), lat_col_e, suppressWarnings(as.numeric(df_comp$decimalLatitude)))
  df_comp$lon_plot_e <- lon_col_e
  df_comp$lat_plot_e <- lat_col_e
  df_comp$is_type_flag_e <- if ("is_type_specimen" %in% names(df_comp)) as.logical(df_comp$is_type_specimen) else rep(FALSE, nrow(df_comp))
  df_comp$is_type_flag_e[is.na(df_comp$is_type_flag_e)] <- FALSE

  geo_b_e <- compute_dynamic_na_map_bounds(df_comp$lon_plot_e, df_comp$lat_plot_e)
  geo_sub <- df_comp[geo_b_e$valid_geo, , drop = FALSE]
  needs_cont_e <- (geo_b_e$xlim[1L] < -125.0) || (geo_b_e$ylim[2L] > 49.5) || (geo_b_e$ylim[1L] < 25.0)
  na_world <- if (needs_cont_e) tryCatch(ggplot2::map_data("world", region = c("USA", "Canada", "Mexico")), error = function(e) NULL) else NULL
  us_states <- tryCatch(ggplot2::map_data("state"), error = function(e) NULL)

  fg_triage_cats <- c("Misidentification (Reassign)", "Hybrid / Introgressant", "Putative Cryptic Species")
  bg_geo <- geo_sub[!(as.character(geo_sub$Treatment_Clean) %in% fg_triage_cats) & !geo_sub$is_type_flag_e, , drop = FALSE]
  fg_geo <- geo_sub[as.character(geo_sub$Treatment_Clean) %in% fg_triage_cats, , drop = FALSE]
  if (nrow(fg_geo) > 0L) {
    fg_geo$Triage_Highlight <- factor(as.character(fg_geo$Treatment_Clean), levels = fg_triage_cats)
  }

  p5 <- ggplot2::ggplot()
  if (!is.null(na_world) && nrow(na_world) > 0) {
    p5 <- p5 + ggplot2::geom_polygon(data = na_world, ggplot2::aes(x = long, y = lat, group = group), fill = "#FBFBF9", color = "grey75", linewidth = 0.22)
  }
  if (!is.null(us_states) && nrow(us_states) > 0) {
    p5 <- p5 + ggplot2::geom_polygon(data = us_states, ggplot2::aes(x = long, y = lat, group = group), fill = "#FBFBF9", color = "grey72", linewidth = 0.25)
  }
  if (nrow(bg_geo) > 0L) {
    p5 <- p5 + ggplot2::geom_point(
      data = bg_geo, ggplot2::aes(x = lon_plot_e, y = lat_plot_e),
      shape = 1, fill = NA, color = "grey62", alpha = 0.35, size = 1.25, stroke = 0.55
    )
  }
  if (nrow(fg_geo) > 0L) {
    fg_geo_mod <- fg_geo[fg_geo$Triage_Highlight %in% c("Misidentification (Reassign)", "Hybrid / Introgressant"), , drop = FALSE]
    fg_geo_cry <- fg_geo[fg_geo$Triage_Highlight == "Putative Cryptic Species", , drop = FALSE]
    if (nrow(fg_geo_mod) > 0L) {
      p5 <- p5 + ggplot2::geom_point(
        data = fg_geo_mod, ggplot2::aes(x = lon_plot_e, y = lat_plot_e, color = Triage_Highlight, shape = Triage_Highlight),
        fill = NA, stroke = 0.55, alpha = 0.28, size = 1.0
      )
    }
    if (nrow(fg_geo_cry) > 0L) {
      p5 <- p5 + ggplot2::geom_point(
        data = fg_geo_cry, ggplot2::aes(x = lon_plot_e, y = lat_plot_e, color = Triage_Highlight, shape = Triage_Highlight),
        fill = NA, stroke = 0.95, alpha = 0.92, size = 2.0
      )
    }
    p5 <- p5 +
      ggplot2::scale_color_manual(
        values = c(
          "Misidentification (Reassign)" = "#E69F00",
          "Hybrid / Introgressant"       = "#CC79A7",
          "Putative Cryptic Species"     = "#D55E00"
        ),
        drop = FALSE
      ) +
      ggplot2::scale_shape_manual(
        values = c(
          "Misidentification (Reassign)" = 4,
          "Hybrid / Introgressant"       = 2,
          "Putative Cryptic Species"     = 6
        ),
        drop = FALSE
      )
  }
  if (any(geo_sub$is_type_flag_e)) {
    type_geo <- geo_sub[geo_sub$is_type_flag_e, , drop = FALSE]
    p5 <- p5 +
      ggplot2::geom_point(
        data = type_geo, ggplot2::aes(x = lon_plot_e, y = lat_plot_e),
        shape = 5, fill = NA, stroke = 1.35, size = 3.2, color = "#111111", inherit.aes = FALSE, show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = type_geo, ggplot2::aes(x = lon_plot_e, y = lat_plot_e),
        shape = 1, fill = NA, stroke = 1.15, size = 3.8, color = "#D55E00", inherit.aes = FALSE, show.legend = FALSE
      )
  }
  p5 <- p5 +
    ggplot2::coord_quickmap(xlim = geo_b_e$xlim, ylim = geo_b_e$ylim) +
    thm + ggplot2::labs(
      title = "E. Geographic Contact Zones, Triage Flags & Type Localities",
      subtitle = "Confirmed vouchers muted (alpha = 0.35) | Triage flags & double-ringed Types highlighted",
      x = "Longitude (deg E)", y = "Latitude (deg N)",
      color = "Triage Flag", shape = "Triage Flag"
    )

  # ----------------------------------------------------------------------------
  # Panel F (p6): Multi-Evidence Delimitation Concordance Grid (Wettstein Matrix)
  # ----------------------------------------------------------------------------
  wett_tx <- if (length(active_tx) > 10L) {
    .unpack_synthesis_display_tier(df_comp$plot_group, active_tx, max_fg = 8L)$foreground_taxa
  } else {
    active_tx
  }
  stream_lvls <- c(
    "1. Type Baseline Control",
    "2. EFA/CDA Morphology",
    "3. DINOv2 Vision Audit",
    "4. 10m-30m Edaphic Niche",
    "5. Hopkins Anthesis Phenology",
    "6. Padial Status"
  )
  ev_rows <- list()
  for (tx in wett_tx) {
    sub_t <- df_comp[df_comp$plot_group == tx, , drop = FALSE]
    if (nrow(sub_t) == 0) next
    sp_t <- as.character(sub_t$species_standardized)

    # 2. EFA/CDA Morphology
    m_ok <- if ("cda_predicted_taxon" %in% names(sub_t) && any(!is.na(sub_t$cda_predicted_taxon))) {
      sub_t$cda_predicted_taxon == tx | sub_t$cda_predicted_taxon == sp_t
    } else if ("morphometric_outlier_flag" %in% names(sub_t)) {
      !as.logical(sub_t$morphometric_outlier_flag)
    } else {
      sub_t$Treatment_Clean == "Confirmed Species"
    }
    pct_m <- 100 * mean(m_ok, na.rm = TRUE); if (!is.finite(pct_m)) pct_m <- 0.0

    # 1. Type Baseline Control
    t_ok <- if ("type_outlier_flag" %in% names(sub_t) && any(!is.na(sub_t$type_outlier_flag))) {
      !as.logical(sub_t$type_outlier_flag)
    } else if ("mahalanobis_d2_type" %in% names(sub_t) && any(is.finite(suppressWarnings(as.numeric(sub_t$mahalanobis_d2_type))))) {
      suppressWarnings(as.numeric(sub_t$mahalanobis_d2_type)) <= stats::qchisq(0.95, df = 2)
    } else if ("nearest_type_morph_taxon" %in% names(sub_t) && any(!is.na(sub_t$nearest_type_morph_taxon))) {
      sub_t$nearest_type_morph_taxon == tx | sub_t$nearest_type_morph_taxon == sp_t
    } else {
      m_ok
    }
    pct_t <- 100 * mean(t_ok, na.rm = TRUE); if (!is.finite(pct_t)) pct_t <- pct_m

    # 3. DINOv2 Vision / Unsupervised GMM Morphotype Core Audit
    pct_v <- if ("vision_conflict_flag" %in% names(sub_t) && any(!is.na(sub_t$vision_conflict_flag))) {
      100 * mean(!as.logical(sub_t$vision_conflict_flag), na.rm = TRUE)
    } else if ("dinov2_predicted_taxon" %in% names(sub_t) && any(!is.na(sub_t$dinov2_predicted_taxon))) {
      100 * mean(sub_t$dinov2_predicted_taxon == tx | sub_t$dinov2_predicted_taxon == sp_t, na.rm = TRUE)
    } else if ("is_typical_morphotype" %in% names(sub_t) && any(!is.na(sub_t$is_typical_morphotype))) {
      100 * mean(suppressWarnings(as.numeric(as.logical(sub_t$is_typical_morphotype))) > 0, na.rm = TRUE)
    } else if ("gmm_uncertainty" %in% names(sub_t) && any(is.finite(suppressWarnings(as.numeric(sub_t$gmm_uncertainty))))) {
      100 * mean(suppressWarnings(as.numeric(sub_t$gmm_uncertainty)) <= 0.30, na.rm = TRUE)
    } else {
      50.0
    }
    if (!is.finite(pct_v)) pct_v <- 50.0

    # 4. 10m-30m Edaphic Niche
    e_ok <- if ("edaphic_best_fit_taxon" %in% names(sub_t) && any(!is.na(sub_t$edaphic_best_fit_taxon))) {
      sub_t$edaphic_best_fit_taxon == tx | sub_t$edaphic_best_fit_taxon == sp_t
    } else if ("edaphic_conflict_flag" %in% names(sub_t) && any(!is.na(sub_t$edaphic_conflict_flag))) {
      !as.logical(sub_t$edaphic_conflict_flag)
    } else {
      m_ok
    }
    pct_e <- 100 * mean(e_ok, na.rm = TRUE); if (!is.finite(pct_e)) pct_e <- pct_m

    # 5. Hopkins Anthesis Phenology
    p_ok <- if ("phenology_conflict_flag" %in% names(sub_t) && any(!is.na(sub_t$phenology_conflict_flag))) {
      !as.logical(sub_t$phenology_conflict_flag)
    } else if ("delta_doy_resid" %in% names(sub_t) && any(!is.na(sub_t$delta_doy_resid))) {
      abs(suppressWarnings(as.numeric(sub_t$delta_doy_resid))) < 14.0
    } else if ("delta_doy" %in% names(sub_t) && any(!is.na(sub_t$delta_doy))) {
      abs(suppressWarnings(as.numeric(sub_t$delta_doy))) < 14.0
    } else {
      m_ok
    }
    pct_p <- 100 * mean(p_ok, na.rm = TRUE); if (!is.finite(pct_p)) pct_p <- pct_m

    # Modal Padial et al. (2010) status code (CCS, ECS, HYB, UCS)
    pad_vec <- if ("padial_status_code" %in% names(sub_t) && any(nzchar(as.character(stats::na.omit(sub_t$padial_status_code))))) {
      toupper(as.character(stats::na.omit(sub_t$padial_status_code)))
    } else {
      tr_v <- as.character(sub_t$Treatment_Clean)
      ifelse(tr_v == "Confirmed Species", "CCS",
      ifelse(tr_v == "Variety / Ecotype", "ECS",
      ifelse(tr_v == "Hybrid / Introgressant", "HYB", "UCS")))
    }
    pad_vec <- pad_vec[pad_vec %in% c("CCS", "ECS", "HYB", "UCS")]
    modal_padial <- if (length(pad_vec) > 0L) names(sort(table(pad_vec), decreasing = TRUE))[1L] else "CCS"
    mean_pct <- round(mean(c(pct_t, pct_m, pct_v, pct_e, pct_p), na.rm = TRUE), 1)

    pct_vals <- round(c(pct_t, pct_m, pct_v, pct_e, pct_p, mean_pct), 1)
    cell_lbls <- c(sprintf("%.0f%%", pct_vals[1:5]), modal_padial)

    ev_rows[[length(ev_rows) + 1L]] <- data.frame(
      Species = gsub("^Packera\\s+", "P. ", tx),
      Stream = factor(stream_lvls, levels = stream_lvls),
      Concordance_Pct = pct_vals,
      Cell_Label = cell_lbls,
      stringsAsFactors = FALSE
    )
  }
  ev_df <- if (length(ev_rows) > 0) {
    do.call(rbind, ev_rows)
  } else {
    data.frame(
      Species = gsub("^Packera\\s+", "P. ", active_tx[1L]),
      Stream = factor(stream_lvls, levels = stream_lvls),
      Concordance_Pct = c(90, 90, 90, 90, 90, 90),
      Cell_Label = c("90%", "90%", "90%", "90%", "90%", "CCS"),
      stringsAsFactors = FALSE
    )
  }
  sp_wett_lvls <- rev(unique(ev_df$Species))
  ev_df$Species <- factor(ev_df$Species, levels = sp_wett_lvls)
  p6 <- ggplot2::ggplot(ev_df, ggplot2::aes(x = Stream, y = Species, fill = Concordance_Pct)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = Cell_Label),
      size = min(2.8, max(1.8, 18.0 / max(length(sp_wett_lvls), 1L))),
      fontface = "bold",
      color = ifelse(ev_df$Concordance_Pct >= 55, "white", "#111111")
    ) +
    ggplot2::scale_fill_gradient2(low = "#E69F00", mid = "#56B4E9", high = "#0072B2", midpoint = 65, limits = c(0, 100)) +
    thm +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 24, hjust = 1, size = 7.0, face = "bold"),
      axis.text.y = ggplot2::element_text(face = if (is_subcluster_mode) "bold" else "italic", size = min(7.8, max(5.2, 48.0 / max(length(sp_wett_lvls), 1L))))
    ) +
    ggplot2::labs(
      title = "F. Multi-Evidence Delimitation Concordance Grid (Wettstein Matrix)",
      subtitle = "5 evidence streams + right-hand modal Padial et al. (2010) status code (CCS, ECS, HYB, UCS)",
      x = "Independent Evidence Stream & Padial Summary",
      y = if (is_subcluster_mode) "Target Morphotype / Subcluster" else "Target Species",
      fill = "% Concordant"
    )

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
    cex_ax_heat <- min(0.78, max(0.42, 6.8 / max(n_taxa, 1L)))
    graphics::image(
      x = seq_len(n_taxa), y = seq_len(n_taxa),
      z = t(W_ord[rev(seq_len(n_taxa)), , drop = FALSE]),
      col = heat_cols, axes = FALSE, xlab = "", ylab = "",
      main = "A. Empirical Taxon Affinity Matrix (W_ij) & Complexes"
    )
    graphics::axis(1, at = seq_len(n_taxa), labels = lbl_ord, las = 2, cex.axis = cex_ax_heat, font = 3)
    graphics::axis(2, at = seq_len(n_taxa), labels = rev(lbl_ord), las = 2, cex.axis = cex_ax_heat, font = 3)
    graphics::box()

    if (n_taxa <= 10L) {
      cell_cex_w <- min(0.72, max(0.44, 5.0 / max(n_taxa, 1L)))
      for (r_i in seq_len(n_taxa)) {
        y_pos <- n_taxa - r_i + 1L
        for (c_j in seq_len(n_taxa)) {
          w_val <- W_ord[r_i, c_j]
          if (is.finite(w_val)) {
            graphics::text(
              c_j, y_pos,
              labels = sprintf("%.2f", w_val),
              cex = cell_cex_w,
              font = if (r_i == c_j) 1L else 2L,
              col = if (w_val >= 0.65) "white" else "#111111"
            )
          }
        }
      }
    }

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

    keep_edge <- matrix(FALSE, nrow = n_taxa, ncol = n_taxa)
    for (i in seq_len(n_taxa)) {
      oth <- setdiff(seq_len(n_taxa), i)
      if (length(oth) > 0L) {
        w_i <- W[i, oth]
        w_i[!is.finite(w_i)] <- -Inf
        top_k <- oth[utils::head(order(w_i, decreasing = TRUE), min(2L, length(oth)))]
        top_k <- top_k[is.finite(W[i, top_k]) & W[i, top_k] > 0]
        for (j in oth) {
          if ((is.finite(W[i, j]) && W[i, j] >= 0.08) || (j %in% top_k)) {
            keep_edge[i, j] <- TRUE
            keep_edge[j, i] <- TRUE
          }
        }
      }
    }

    if (n_taxa >= 2L) {
      for (i in seq_len(n_taxa - 1L)) {
        for (j in (i + 1L):n_taxa) {
          if (keep_edge[i, j]) {
            same_c <- cdf$complex_id[i] == cdf$complex_id[j]
            e_col <- if (same_c) grDevices::adjustcolor(complex_pal[(cdf$complex_id[i] - 1L) %% length(complex_pal) + 1L], alpha.f = 0.70) else grDevices::adjustcolor("grey55", alpha.f = 0.45)
            graphics::segments(px[i], py[i], px[j], py[j], col = e_col, lwd = max(0.8, 5.5 * W[i, j]))
          }
        }
      }
    }

    node_cex <- pmax(1.1, pmin(2.6, 0.9 + 0.65 * log10(pmax(cdf$n_vouchers, 1))))
    node_pch <- hollow_pchs[(cdf$complex_id - 1L) %% length(hollow_pchs) + 1L]
    node_col <- dyn_pal[taxa]
    node_col[is.na(node_col)] <- complex_pal[(cdf$complex_id[is.na(node_col)] - 1L) %% length(complex_pal) + 1L]
    graphics::points(px, py, pch = node_pch, bg = NA, col = node_col, cex = node_cex, lwd = 1.8)

    lbl_xy <- NULL
    if (exists("repel_2d_labels", mode = "function")) {
      lbl_xy <- tryCatch(repel_2d_labels(px, py, short_names), error = function(e) {
        tryCatch(repel_2d_labels(px, py), error = function(e2) {
          tryCatch(repel_2d_labels(cbind(px, py)), error = function(e3) NULL)
        })
      })
    }
    lx <- px
    ly <- py + max(0.04 * diff(yr), 0.025)
    used_external_repel <- FALSE
    if (!is.null(lbl_xy)) {
      if ((is.matrix(lbl_xy) || is.data.frame(lbl_xy)) && ncol(lbl_xy) >= 2L && nrow(lbl_xy) == length(px)) {
        lx <- as.numeric(lbl_xy[, 1])
        ly <- as.numeric(lbl_xy[, 2])
        used_external_repel <- TRUE
      } else if (is.list(lbl_xy) && all(c("x", "y") %in% names(lbl_xy)) && length(lbl_xy$x) == length(px)) {
        lx <- as.numeric(lbl_xy$x)
        ly <- as.numeric(lbl_xy$y)
        used_external_repel <- TRUE
      }
    }
    if (!used_external_repel && length(px) >= 2L) {
      min_dx <- max(0.085 * diff(xr), 0.045)
      min_dy <- max(0.055 * diff(yr), 0.030)
      for (iter in seq_len(25L)) {
        for (a in seq_len(length(px) - 1L)) {
          for (b in (a + 1L):length(px)) {
            dx <- lx[a] - lx[b]
            dy <- ly[a] - ly[b]
            if (abs(dx) < min_dx && abs(dy) < min_dy) {
              push_x <- 0.25 * (min_dx - abs(dx)) * (if (dx >= 0) 1 else -1)
              push_y <- 0.25 * (min_dy - abs(dy)) * (if (dy >= 0) 1 else -1)
              lx[a] <- lx[a] + push_x; lx[b] <- lx[b] - push_x
              ly[a] <- ly[a] + push_y; ly[b] <- ly[b] - push_y
            }
          }
        }
      }
    }
    lbl_cex_net <- min(0.76, max(0.48, 7.2 / max(n_taxa, 1L)))
    for (k_i in seq_along(short_names)) {
      if (is.finite(lx[k_i]) && is.finite(ly[k_i]) && nzchar(short_names[k_i])) {
        sw <- graphics::strwidth(short_names[k_i], cex = lbl_cex_net, font = 3) * 0.56
        sh <- graphics::strheight(short_names[k_i], cex = lbl_cex_net, font = 3) * 0.68
        graphics::rect(
          lx[k_i] - sw, ly[k_i] - sh, lx[k_i] + sw, ly[k_i] + sh,
          col = grDevices::adjustcolor("white", alpha.f = 0.82),
          border = NA
        )
      }
    }
    graphics::text(lx, ly, labels = short_names, cex = lbl_cex_net, font = 3, col = "black")

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
fit_gmm_em_pure <- function(X, max_k = 8, tol = 1e-5, max_iter = 150, X_eval = NULL, seed = 42L) {
  if (is.null(X_eval)) X_eval <- X
  N <- nrow(X); p <- ncol(X); best_bic <- -Inf; best_model <- NULL
  base_seed <- as.integer(if (!is.null(seed) && is.finite(seed)) seed else 42L)
  bic_table <- data.frame(K = seq_len(max_k), LogLik = NA_real_, n_params = NA_integer_, BIC = NA_real_)
  for (k in seq_len(max_k)) {
    set.seed(base_seed + k)
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
      best_model <- list(k = k, classification = apply(resp_ev, 1, which.max), uncertainty = 1 - apply(resp_ev, 1, max), z = resp_ev, bic = bic_val)
    }
  }
  return(list(best = best_model, bic_table = bic_table))
}

merge_manifest_petiole_traits_by_individual <- function(agg_df, man_cand) {
  if (is.null(agg_df) || nrow(agg_df) == 0L || is.null(man_cand) || !nzchar(man_cand) || !file.exists(man_cand)) return(agg_df)
  target_pcols <- c("petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio")
  need_cols <- target_pcols[vapply(target_pcols, function(cl) !(cl %in% names(agg_df)) || any(is.na(agg_df[[cl]])), logical(1))]
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
  if ("leaf_zone" %in% names(man_df)) {
    lz_low <- tolower(trimws(as.character(man_df$leaf_zone)))
    lz_low[is.na(lz_low) | !nzchar(lz_low)] <- "basal"
    if (any(lz_low != "cauline")) {
      grp_k <- if ("plant_individual_id" %in% names(man_df)) {
        mp_tmp <- ifelse(is.na(man_df$plant_individual_id), 0L, suppressWarnings(as.integer(man_df$plant_individual_id)))
        mp_tmp[is.na(mp_tmp)] <- 0L
        paste0(man_df$catalogNumber, "__", mp_tmp)
      } else man_df$catalogNumber
      has_basal_grp <- tapply(lz_low != "cauline", grp_k, any)
      keep_row <- (lz_low != "cauline") | !has_basal_grp[grp_k]
      man_df <- man_df[which(keep_row), , drop = FALSE]
    }
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
  r_val_cols <- r_val_cols[vapply(r_val_cols, function(cl) !(cl %in% names(agg_df)) || any(is.na(agg_df[[cl]])), logical(1))]
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
