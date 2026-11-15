# ==============================================================================
# Module: scripts/figures/step06_figures.R
# Project: Plantmatics
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   300-DPI PNG-only publication figure suite for Standalone Step 06 Phylogenomics:
#     - fig06a_genomic_pcoa_tree_concordance.png (PCoA + ML Tree + Admixture + Congruence)
#     - fig06b_reticulation_gdi_landscape.png ((h, Q12) + f-branch + gdi + FEEMS/MMRR)
#   Strictly uses unfilled hollow scatterplot markers (pch = c(1,2,5,0,6,4), fill = NA)
#   and ICN italicized Packera binomials. Zero PDF outputs.
# ==============================================================================

.format_tree_tip_expr <- function(sample_id, taxon_name) {
  tx <- ifelse(is.na(taxon_name) | !nzchar(taxon_name), "", as.character(taxon_name))
  sid <- gsub("'", "", as.character(sample_id))
  if (!nzchar(tx)) return(paste0("'", sid, "'"))
  if (grepl("sp\\.\\s*nov", tx, ignore.case = TRUE)) {
    rest <- trimws(sub("^Packera\\s+", "", tx))
    return(sprintf("paste(italic('P.'), ' %s [', '%s', ']')", gsub("'", "", rest), sid))
  }
  short_tx <- sub("^Packera\\s+", "P. ", tx)
  sprintf("paste(italic('%s'), ' [', '%s', ']')", gsub("'", "", short_tx), sid)
}

.midpoint_root_phylo_safe <- function(tr) {
  if (is.null(tr) || length(tr$tip.label) < 3L) return(tr)
  if (requireNamespace("phangorn", quietly = TRUE)) {
    res <- tryCatch(phangorn::midpoint(tr), error = function(e) NULL)
    if (!is.null(res) && inherits(res, "phylo")) return(res)
  }
  if (isTRUE(ape::is.rooted(tr)) || is.null(tr$edge.length)) return(tr)
  tryCatch({
    n_t <- length(tr$tip.label)
    dmat <- ape::cophenetic.phylo(tr)
    idx_mx <- which.max(dmat)
    ij <- arrayInd(idx_mx, dim(dmat))
    i_tip <- match(rownames(dmat)[ij[1L]], tr$tip.label)
    j_tip <- match(colnames(dmat)[ij[2L]], tr$tip.label)
    if (!is.finite(i_tip) || !is.finite(j_tip) || i_tip == j_tip) return(tr)
    half_d <- 0.5 * dmat[ij[1L], ij[2L]]
    if (!is.finite(half_d) || half_d <= 0.0) return(tr)
    np <- ape::nodepath(tr, from = i_tip, to = j_tip)
    if (length(np) < 2L) return(tr)
    cum <- 0.0
    split_edge <- NA_integer_
    for (k in seq_len(length(np) - 1L)) {
      a <- np[k]; b <- np[k + 1L]
      e_row <- which((tr$edge[, 1L] == a & tr$edge[, 2L] == b) | (tr$edge[, 1L] == b & tr$edge[, 2L] == a))
      if (length(e_row) < 1L) next
      el <- max(0.0, tr$edge.length[e_row[1L]])
      if (cum + el >= half_d) {
        split_edge <- e_row[1L]
        break
      }
      cum <- cum + el
    }
    if (is.na(split_edge)) return(tr)
    ch_node <- tr$edge[split_edge, 2L]
    rt <- if (ch_node <= n_t) {
      ape::root(tr, outgroup = ch_node, resolve.root = TRUE)
    } else {
      ape::root(tr, node = ch_node, resolve.root = TRUE)
    }
    if (!is.null(rt$edge.length) && length(rt$tip.label) >= 3L) {
      n_rt_tips <- length(rt$tip.label)
      root_id <- n_rt_tips + 1L
      r_rows <- which(rt$edge[, 1L] == root_id)
      if (length(r_rows) == 2L) {
        tot_r <- sum(pmax(0.0, rt$edge.length[r_rows]))
        nd_depths <- ape::node.depth.edgelength(rt)
        sub_tips <- function(nd) {
          queue <- nd; res_t <- integer(0)
          while (length(queue) > 0L) {
            cur <- queue[1L]; queue <- queue[-1L]
            if (cur <= n_rt_tips) res_t <- c(res_t, cur)
            else queue <- c(queue, rt$edge[rt$edge[, 1L] == cur, 2L])
          }
          res_t
        }
        d1 <- max(nd_depths[sub_tips(rt$edge[r_rows[1L], 2L])] - nd_depths[rt$edge[r_rows[1L], 2L]])
        d2 <- max(nd_depths[sub_tips(rt$edge[r_rows[2L], 2L])] - nd_depths[rt$edge[r_rows[2L], 2L]])
        new_e1 <- max(0.0, min(tot_r, 0.5 * (tot_r + d2 - d1)))
        rt$edge.length[r_rows[1L]] <- new_e1
        rt$edge.length[r_rows[2L]] <- max(0.0, tot_r - new_e1)
      }
    }
    rt
  }, error = function(e) tr)
}

.build_phylo_panel_gg <- function(tr, plot_df, pal_vec, sh_vec) {
  if (is.null(tr) || length(tr$tip.label) < 3L || !requireNamespace("ape", quietly = TRUE)) {
    return(
      ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$PCo1, y = .data$PCo2, color = .data$display_taxon, shape = .data$display_taxon)) +
        ggplot2::geom_point(size = 3.0, fill = NA, stroke = 1.0) +
        ggplot2::labs(title = "B. Phylogeny Not Supplied (PCoA Mirror)", x = "PCo1", y = "PCo2") +
        theme_packera_pub(base_size = 9.5) +
        ggplot2::theme(legend.position = "none")
    )
  }

  tr_work <- tr
  tr_work$tip.label <- gsub("_sorted$", "", trimws(tr_work$tip.label))
  tr_work <- .midpoint_root_phylo_safe(tr_work)
  tr_u <- ape::ladderize(tr_work, right = FALSE)
  n_tips <- length(tr_u$tip.label)
  tot_n <- n_tips + tr_u$Nnode

  ord <- ape::reorder.phylo(tr_u, order = "cladewise")
  depth_x <- if (!is.null(tr_u$edge.length)) {
    ape::node.depth.edgelength(tr_u)
  } else {
    dx <- numeric(tot_n)
    for (e_i in seq_len(nrow(ord$edge))) {
      dx[ord$edge[e_i, 2L]] <- dx[ord$edge[e_i, 1L]] + 1.0
    }
    dx
  }

  y_pos <- numeric(tot_n)
  tip_ord <- ord$edge[ord$edge[, 2L] <= n_tips, 2L]
  y_pos[tip_ord] <- seq_len(n_tips)
  post_ord <- ape::reorder.phylo(tr_u, order = "postorder")
  post_parents <- unique(post_ord$edge[, 1L])
  for (nd in post_parents) {
    chs <- post_ord$edge[post_ord$edge[, 1L] == nd, 2L]
    y_pos[nd] <- mean(y_pos[chs])
  }

  h_segs <- data.frame(
    x = depth_x[ord$edge[, 1L]], xend = depth_x[ord$edge[, 2L]],
    y = y_pos[ord$edge[, 2L]], yend = y_pos[ord$edge[, 2L]], stringsAsFactors = FALSE
  )
  v_segs <- do.call(rbind, lapply(post_parents, function(nd) {
    chs <- post_ord$edge[post_ord$edge[, 1L] == nd, 2L]
    data.frame(x = depth_x[nd], xend = depth_x[nd], y = min(y_pos[chs]), yend = max(y_pos[chs]), stringsAsFactors = FALSE)
  }))

  sp_lkp <- stats::setNames(plot_df$display_taxon, plot_df$genetic_sample_id)
  tip_df <- data.frame(
    tip = tr_u$tip.label, x = depth_x[seq_len(n_tips)], y = y_pos[seq_len(n_tips)],
    display_taxon = ifelse(tr_u$tip.label %in% names(sp_lkp), sp_lkp[tr_u$tip.label], "Other"),
    stringsAsFactors = FALSE
  )
  tip_df$expr_lbl <- vapply(seq_len(nrow(tip_df)), function(k) .format_tree_tip_expr(tip_df$tip[k], tip_df$display_taxon[k]), character(1))
  max_x <- max(c(tip_df$x, 1e-4), na.rm = TRUE)

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = h_segs, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend), color = "#334155", linewidth = 0.45) +
    ggplot2::geom_segment(data = v_segs, ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend), color = "#334155", linewidth = 0.45) +
    ggplot2::geom_point(data = tip_df, ggplot2::aes(x = .data$x, y = .data$y, color = .data$display_taxon, shape = .data$display_taxon), size = 2.6, fill = NA, stroke = 0.95) +
    ggplot2::geom_text(data = tip_df, ggplot2::aes(x = .data$x + 0.025 * max_x, y = .data$y, label = .data$expr_lbl), parse = TRUE, hjust = 0, size = 2.55, color = "#0f172a") +
    ggplot2::scale_color_manual(values = pal_vec) +
    ggplot2::scale_shape_manual(values = sh_vec) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.56))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "B. Maximum-Likelihood Phylogeny (RAxML / Patristic Tips)", x = "Substitutions per Site", y = NULL) +
    theme_packera_pub(base_size = 9.5) +
    ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
}

.build_concordance_panel_gg <- function(plot_df, concordance_df, db_path, run_id) {
  linked <- plot_df[!is.na(plot_df$catalogNumber) & nzchar(plot_df$catalogNumber), , drop = FALSE]
  rows <- list()
  if (nrow(linked) > 0L) {
    morph_map <- NULL; rf_map <- NULL
    if (!is.null(db_path) && nzchar(db_path) && file.exists(db_path) && exists("db_table_exists_safe", mode = "function")) {
      conn <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), dbname = db_path, flags = RSQLite::SQLITE_RO), error = function(e) NULL)
      if (!is.null(conn)) {
        on.exit(tryCatch(DBI::dbDisconnect(conn), error = function(e) NULL), add = TRUE)
        q_cats <- paste(vapply(unique(linked$catalogNumber), function(s) as.character(DBI::dbQuoteString(conn, s)), character(1)), collapse = ", ")
        if (db_table_exists_safe(conn, "otu_morphometrics")) {
          morph_map <- tryCatch(DBI::dbGetQuery(conn, sprintf("SELECT catalogNumber, cda_predicted_taxon FROM otu_morphometrics WHERE run_id = %s AND catalogNumber IN (%s) GROUP BY catalogNumber;", as.character(DBI::dbQuoteString(conn, run_id)), q_cats)), error = function(e) NULL)
        }
        if (db_table_exists_safe(conn, "spatial_rf_predictions")) {
          rf_map <- tryCatch(DBI::dbGetQuery(conn, sprintf("SELECT catalogNumber, edaphic_best_fit_taxon, rf_predicted_taxon FROM spatial_rf_predictions WHERE run_id = %s AND catalogNumber IN (%s) GROUP BY catalogNumber;", as.character(DBI::dbQuoteString(conn, run_id)), q_cats)), error = function(e) NULL)
        }
      }
    }
    for (k in seq_len(nrow(linked))) {
      cat_k <- linked$catalogNumber[k]; sid_k <- linked$genetic_sample_id[k]
      m_tx <- if (!is.null(morph_map) && cat_k %in% morph_map$catalogNumber) morph_map$cda_predicted_taxon[match(cat_k, morph_map$catalogNumber)] else NA_character_
      e_tx <- if (!is.null(rf_map) && cat_k %in% rf_map$catalogNumber) {
        rr <- rf_map[match(cat_k, rf_map$catalogNumber), , drop = FALSE]
        ifelse(!is.na(rr$edaphic_best_fit_taxon) & nzchar(rr$edaphic_best_fit_taxon), rr$edaphic_best_fit_taxon, rr$rf_predicted_taxon)
      } else NA_character_
      rows[[length(rows) + 1L]] <- data.frame(sample = sid_k, stream = c("Herbarium_ID", "Morph_CDA", "Eco_Niche", "Genomic_DAPC"), taxon = c(linked$species_standardized[k], m_tx, e_tx, linked$geno_predicted_taxon[k]), stringsAsFactors = FALSE)
    }
  } else {
    for (k in seq_len(nrow(plot_df))) {
      rows[[length(rows) + 1L]] <- data.frame(sample = plot_df$genetic_sample_id[k], stream = c("Field_Label", "Genomic_DAPC"), taxon = c(plot_df$species_standardized[k], plot_df$geno_predicted_taxon[k]), stringsAsFactors = FALSE)
    }
  }
  tile_df <- dplyr::bind_rows(rows)
  tile_df$taxon <- ifelse(is.na(tile_df$taxon) | !nzchar(tile_df$taxon), "Unlinked", tile_df$taxon)
  t_levels <- sort(unique(tile_df$taxon))
  t_pal <- c(get_combinatorial_palette_r(setdiff(t_levels, "Unlinked")), Unlinked = "#e2e8f0")
  fin_conc <- concordance_df[is.finite(concordance_df$mantel_r), , drop = FALSE]
  mantel_sub <- if (nrow(fin_conc) > 0L) paste(sprintf("%s: r=%.2f (p=%.3f)", gsub("_", " ", fin_conc$comparison), fin_conc$mantel_r, fin_conc$mantel_p), collapse = " | ") else "Strict upstream independence preserved"
  lab_fn <- if (exists("icn_discrete_labeller", mode = "function")) icn_discrete_labeller() else ggplot2::waiver()

  ggplot2::ggplot(tile_df, ggplot2::aes(x = .data$stream, y = .data$sample, fill = .data$taxon)) +
    ggplot2::geom_tile(color = "#ffffff", linewidth = 0.45) +
    ggplot2::scale_fill_manual(values = t_pal, labels = lab_fn, name = "Stream Call") +
    ggplot2::labs(title = "D. Voucher Cross-Modal Agreement (Post-Hoc)", subtitle = mantel_sub, x = "Independent Stream", y = "Linked Sample") +
    theme_packera_pub(base_size = 9.0) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1), axis.text.y = ggplot2::element_text(size = 7.4))
}

render_fig06a_genomic_pcoa_tree_png <- function(out_png, pass_df, assign_df, q_mat, pco_eig_pct, tr, concordance_df, db_path, run_id) {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") ||
      identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) return(invisible(""))

  plot_df <- dplyr::left_join(pass_df, assign_df[, setdiff(names(assign_df), c("catalogNumber", "population_id")), drop = FALSE], by = "genetic_sample_id")
  plot_df$display_taxon <- ifelse(!is.na(plot_df$species_standardized) & nzchar(plot_df$species_standardized), plot_df$species_standardized, ifelse(!is.na(plot_df$geno_predicted_taxon), plot_df$geno_predicted_taxon, "Unknown"))
  tx_levels <- sort(unique(plot_df$display_taxon))
  pal_vec <- get_combinatorial_palette_r(tx_levels)
  sh_vec <- get_combinatorial_hollow_shapes_r(tx_levels)
  lab_fn <- if (exists("icn_discrete_labeller", mode = "function")) icn_discrete_labeller() else ggplot2::waiver()

  rep_segs <- plot_df %>%
    dplyr::filter(!is.na(.data$replicate_pair_id), nzchar(.data$replicate_pair_id)) %>%
    dplyr::group_by(.data$replicate_pair_id) %>%
    dplyr::filter(dplyr::n() == 2L) %>%
    dplyr::summarise(x = .data$PCo1[1L], y = .data$PCo2[1L], xend = .data$PCo1[2L], yend = .data$PCo2[2L], .groups = "drop")

  pair_key <- ifelse(!is.na(plot_df$replicate_pair_id) & nzchar(plot_df$replicate_pair_id), plot_df$replicate_pair_id, plot_df$genetic_sample_id)
  plot_df$rep_label <- ifelse(!duplicated(pair_key), plot_df$genetic_sample_id, "")

  p_a <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$PCo1, y = .data$PCo2)) +
    ggplot2::geom_segment(data = rep_segs, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend), linetype = "dashed", color = "#64748b", linewidth = 0.45, inherit.aes = FALSE) +
    ggplot2::geom_point(ggplot2::aes(color = .data$display_taxon, shape = .data$display_taxon), size = 3.1, fill = NA, stroke = 1.05) +
    ggplot2::scale_color_manual(values = pal_vec, labels = lab_fn, name = "Taxon") +
    ggplot2::scale_shape_manual(values = sh_vec, labels = lab_fn, name = "Taxon") +
    ggplot2::labs(title = "A. Polyploid-Dosage Genomic PCoA (M/P Pairs Dashed)", x = sprintf("PCo1 (%0.1f%%)", pco_eig_pct[1L]), y = sprintf("PCo2 (%0.1f%%)", pco_eig_pct[min(2L, length(pco_eig_pct))])) +
    theme_packera_pub(base_size = 9.5) +
    ggplot2::theme(legend.position = "right")

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_a <- p_a + ggrepel::geom_text_repel(data = plot_df, ggplot2::aes(label = .data$rep_label, color = .data$display_taxon), size = 2.55, box.padding = 0.42, point.padding = 0.28, min.segment.length = 0.15, force = 5.0, max.overlaps = Inf, seed = 42, show.legend = FALSE)
  }

  p_b <- .build_phylo_panel_gg(tr, plot_df, pal_vec, sh_vec)

  q_df_long <- do.call(rbind, lapply(seq_len(nrow(q_mat)), function(i) {
    data.frame(genetic_sample_id = rownames(q_mat)[i], ancestry_taxon = colnames(q_mat), q_prop = as.numeric(q_mat[i, ]), stringsAsFactors = FALSE)
  }))
  ord_ids <- plot_df$genetic_sample_id[order(plot_df$is_putative_novel_or_hybrid_query, plot_df$display_taxon, plot_df$genetic_sample_id)]
  q_df_long$genetic_sample_id <- factor(q_df_long$genetic_sample_id, levels = ord_ids)
  anc_pal <- get_combinatorial_palette_r(sort(unique(c(tx_levels, colnames(q_mat)))))

  p_c <- ggplot2::ggplot(q_df_long, ggplot2::aes(x = .data$genetic_sample_id, y = .data$q_prop, fill = .data$ancestry_taxon)) +
    ggplot2::geom_col(width = 0.82, color = "#ffffff", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = anc_pal, labels = lab_fn, name = "Ancestry") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.0, 0.03))) +
    ggplot2::labs(title = "C. Genomic Ancestry Simplex & Co-ancestry Structure", x = "QC-Passing Samples", y = "Ancestry Proportion (q)") +
    theme_packera_pub(base_size = 9.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 7.2), legend.position = "right")

  p_d <- .build_concordance_panel_gg(plot_df, concordance_df, db_path, run_id)

  combo <- (p_a | p_b) / (p_c | p_d) +
    patchwork::plot_annotation(
      title = "Step 06: Independent Nuclear & Polyploid-Dosage Phylogenomics",
      subtitle = sprintf("QC-passing genomic cohort (n = %d; strict upstream independence preserved)", nrow(plot_df)),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.0), plot.subtitle = ggplot2::element_text(size = 9.5, color = "#475569"))
    )
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  save_packera_png_300(combo, out_png, width = 15.0, height = 11.0, dpi = 300L)
  invisible(out_png)
}

render_fig06b_reticulation_gdi_png <- function(out_png, tri_df, fbranch_df, gdi_df, feems_edges_df, mmrr_summary_df, organellar_placement_df = NULL) {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") ||
      identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) return(invisible(""))

  tri_bound <- data.frame(h = c(0.0, 0.5, 1.0, 0.0), q12 = c(0.0, 1.0, 0.0, 0.0))
  tri_ok <- tri_df[is.finite(tri_df$hybrid_index_h) & is.finite(tri_df$interspecific_ancestry_q12), , drop = FALSE]
  p_tri <- ggplot2::ggplot() +
    ggplot2::geom_path(data = tri_bound, ggplot2::aes(x = .data$h, y = .data$q12), linetype = "dashed", color = "#4B5563", linewidth = 0.7) +
    ggplot2::labs(title = "A. Interspecific Ancestry Triangle (h vs. Q12)", x = "Hybrid Index (h)", y = "Interspecific Ancestry (Q12)") +
    ggplot2::coord_cartesian(xlim = c(-0.04, 1.04), ylim = c(-0.04, 1.05)) +
    theme_packera_pub(base_size = 9.5) +
    ggplot2::theme(legend.title = ggplot2::element_blank())
  if (nrow(tri_ok) > 0L) {
    cls_lev <- sort(unique(tri_ok$hybrid_generation_class))
    sh_c <- stats::setNames(rep(c(1L, 2L, 5L, 0L, 6L, 4L), length.out = length(cls_lev)), cls_lev)
    p_tri <- p_tri +
      ggplot2::geom_point(data = tri_ok, ggplot2::aes(x = .data$hybrid_index_h, y = .data$interspecific_ancestry_q12, color = .data$hybrid_generation_class, shape = .data$hybrid_generation_class), size = 3.0, stroke = 1.15, fill = NA) +
      ggplot2::scale_shape_manual(values = sh_c)
  }

  pol_vals <- if (!is.null(fbranch_df) && "polarization_source" %in% names(fbranch_df)) {
    v <- as.character(fbranch_df$polarization_source); unique(v[!is.na(v) & nzchar(v)])
  } else character(0)
  fb_sub <- if (length(pol_vals) > 0L) sprintf("Polarization: %s", pol_vals[1L]) else NULL
  p_fb <- if (!is.null(fbranch_df) && nrow(fbranch_df) > 0L) {
    ggplot2::ggplot(fbranch_df, ggplot2::aes(x = .data$donor_taxon, y = .data$recipient_taxon, fill = .data$f_branch_fb)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.25) +
      ggplot2::scale_fill_gradient(low = "#F8FAFC", high = "#B91C1C", limits = c(0, 1), name = "f_b") +
      ggplot2::labs(title = "B. Tree-Guided Introgression (f-branch)", subtitle = fb_sub, x = "Candidate Donor (P3)", y = "Recipient Branch") +
      theme_packera_pub(base_size = 9.5) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 8.0, color = "#475569"), axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, face = "italic", size = 8), axis.text.y = ggplot2::element_text(face = "italic", size = 8))
  } else {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Fewer than 4 lineages for f-branch ABBA/BABA", size = 3.5) + ggplot2::theme_void() + ggplot2::labs(title = "B. Tree-Guided Introgression (f-branch)")
  }

  p_gdi <- if (!is.null(gdi_df) && nrow(gdi_df) > 0L && any(is.finite(gdi_df$coalescent_gdi))) {
    gd <- gdi_df[is.finite(gdi_df$coalescent_gdi), , drop = FALSE]
    gd$taxon <- factor(gd$taxon, levels = gd$taxon[order(gd$coalescent_gdi)])
    cls_gdi <- sort(unique(gd$gdi_classification))
    sh_g <- stats::setNames(rep(c(1L, 2L, 5L, 0L), length.out = length(cls_gdi)), cls_gdi)
    pal_g <- get_combinatorial_palette_r(cls_gdi)
    ggplot2::ggplot(gd, ggplot2::aes(x = .data$coalescent_gdi, y = .data$taxon, color = .data$gdi_classification, shape = .data$gdi_classification)) +
      ggplot2::geom_vline(xintercept = 0.20, linetype = "dashed", color = "#DC2626") +
      ggplot2::geom_vline(xintercept = 0.70, linetype = "dashed", color = "#059669") +
      ggplot2::geom_segment(ggplot2::aes(x = .data$gdi_ci_low, xend = .data$gdi_ci_high, y = .data$taxon, yend = .data$taxon), linewidth = 0.65) +
      ggplot2::geom_point(size = 3.0, stroke = 1.15, fill = NA) +
      ggplot2::scale_color_manual(values = pal_g, name = "MSC Gate") +
      ggplot2::scale_shape_manual(values = sh_g, name = "MSC Gate") +
      ggplot2::scale_x_continuous(limits = c(0, 1.02)) +
      ggplot2::labs(title = "C. Invariant-Site Coalescent Delimitation (gdi)", subtitle = "Oversplit < 0.20 | Distinct Species > 0.70", x = "Genealogical Divergence Index (gdi +/- 95% CI)", y = NULL) +
      theme_packera_pub(base_size = 9.5) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic", size = 8.5))
  } else {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Coalescent gdi metrics unavailable", size = 3.5) + ggplot2::theme_void() + ggplot2::labs(title = "C. Coalescent Delimitation (gdi)")
  }

  mm_sub <- if (!is.null(mmrr_summary_df) && nrow(mmrr_summary_df) > 0L && is.finite(mmrr_summary_df$mmrr_beta_ibd[1L])) {
    sprintf("MMRR: b_IBD=%.2f (p=%.3f), b_IBE=%.2f (p=%.3f)", mmrr_summary_df$mmrr_beta_ibd[1L], mmrr_summary_df$mmrr_p_ibd[1L], mmrr_summary_df$mmrr_beta_ibe[1L], mmrr_summary_df$mmrr_p_ibe[1L])
  } else {
    "Insufficient georeferenced samples for MMRR/FEEMS"
  }
  p_feems <- if (!is.null(feems_edges_df) && nrow(feems_edges_df) > 0L) {
    pts <- unique(rbind(data.frame(lon = feems_edges_df$lon_a, lat = feems_edges_df$lat_a), data.frame(lon = feems_edges_df$lon_b, lat = feems_edges_df$lat_b)))
    feems_plot <- NULL
    if (exists("build_na_projected_basemap_gg", mode = "function") && exists("project_lonlat_points_r", mode = "function")) {
      map_ctx <- build_na_projected_basemap_gg(pts$lon, pts$lat, target_crs = 5070L)
      pts_p <- project_lonlat_points_r(pts, map_ctx, lon_col = "lon", lat_col = "lat")
      ea_p <- project_lonlat_points_r(feems_edges_df, map_ctx, lon_col = "lon_a", lat_col = "lat_a")
      eb_p <- project_lonlat_points_r(feems_edges_df, map_ctx, lon_col = "lon_b", lat_col = "lat_b")
      if (nrow(ea_p) == nrow(feems_edges_df) && nrow(eb_p) == nrow(feems_edges_df)) {
        seg_p <- data.frame(x = ea_p$map_x, y = ea_p$map_y, xend = eb_p$map_x, yend = eb_p$map_y, feems_log_migration = feems_edges_df$feems_log_migration)
        feems_plot <- map_ctx$plot +
          ggplot2::geom_segment(data = seg_p, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, color = .data$feems_log_migration), linewidth = 1.0, inherit.aes = FALSE) +
          ggplot2::geom_point(data = pts_p, ggplot2::aes(x = .data$map_x, y = .data$map_y), shape = 1L, size = 2.8, stroke = 1.15, fill = NA, color = "#111827", inherit.aes = FALSE) +
          ggplot2::scale_color_gradient2(low = "#B91C1C", mid = "#E5E7EB", high = "#1D4ED8", midpoint = 0.0, name = "log(m)") +
          ggplot2::labs(title = "D. Effective Migration Surface (FEEMS) & MMRR", subtitle = mm_sub) +
          theme_packera_pub(base_size = 9.5)
      }
    }
    if (is.null(feems_plot)) {
      feems_plot <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = feems_edges_df, ggplot2::aes(x = .data$lon_a, y = .data$lat_a, xend = .data$lon_b, yend = .data$lat_b, color = .data$feems_log_migration), linewidth = 1.0) +
        ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$lon, y = .data$lat), shape = 1L, size = 2.6, stroke = 1.1, fill = NA, color = "#111827") +
        ggplot2::scale_color_gradient2(low = "#B91C1C", mid = "#E5E7EB", high = "#1D4ED8", midpoint = 0.0, name = "log(m)") +
        ggplot2::labs(title = "D. Effective Migration Surface (FEEMS) & MMRR", subtitle = mm_sub, x = "Longitude", y = "Latitude") +
        theme_packera_pub(base_size = 9.5)
    }
    feems_plot
  } else {
    ggplot2::ggplot() + ggplot2::annotate("text", x = 0.5, y = 0.5, label = mm_sub, size = 3.5) + ggplot2::theme_void() + ggplot2::labs(title = "D. Landscape Genomics (FEEMS / MMRR)")
  }

  ann_sub <- if (!is.null(organellar_placement_df) && nrow(organellar_placement_df) > 0L && ("locus" %in% names(organellar_placement_df))) {
    ltab <- table(as.character(organellar_placement_df$locus))
    sprintf("GenBank/BOLD 11-mer locus placement: n=%d (%s)", nrow(organellar_placement_df), paste(sprintf("%s=%d", names(ltab), as.integer(ltab)), collapse = ", "))
  } else NULL
  combo <- (p_tri | p_fb) / (p_gdi | p_feems) +
    patchwork::plot_annotation(
      title = "Step 06: Reticulation Dynamics, Coalescent gdi, & Landscape Migration",
      subtitle = ann_sub,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12.0), plot.subtitle = ggplot2::element_text(size = 9.0, color = "#475569"))
    )
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  save_packera_png_300(combo, out_png, width = 15.0, height = 10.5, dpi = 300L)
  invisible(out_png)
}
