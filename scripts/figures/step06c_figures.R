# ==============================================================================
# Module: scripts/figures/step06c_figures.R
# Project: Plantmatics | Affiliation: UNC-Chapel Hill Herbarium (NCU)
# 300-DPI PNG-only 4-panel Step 06 audit plate (`fig06c_cytonuclear_paralog_cytotype_audit.png`)
# ==============================================================================

if (!exists("save_packera_png_300", mode = "function")) {
  for (.p in c("scripts/figures/theme_packera_pub.R", "../figures/theme_packera_pub.R")) {
    if (file.exists(.p)) {
      source(.p)
      break
    }
  }
}

.clean_str_06c <- function(v) {
  x <- trimws(as.character(v))
  x[is.na(x) | !nzchar(x) | tolower(x) %in% c("na", "nan", "null", "none", "unknown", "unassigned")] <- NA_character_
  x
}

.empty_panel_06c <- function(title_txt, msg_txt, sub_txt = NULL) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg_txt, size = 3.25, color = "#334155") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_packera_pub(base_size = 9.3) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank()) +
    ggplot2::labs(title = title_txt, subtitle = sub_txt, x = NULL, y = NULL)
}

.build_panel06c_hdplot_gg <- function(parsed_dir) {
  t_a <- "A. HDplot Paralog Screening & ParaLoss-1272 Locus Weights"
  hd_path <- if (!is.null(parsed_dir) && nzchar(parsed_dir)) file.path(parsed_dir, "hdplot_locus_diagnostics.tsv") else ""
  if (!nzchar(hd_path) || !file.exists(hd_path)) {
    return(.empty_panel_06c(t_a, "HDplot diagnostics (hdplot_locus_diagnostics.tsv) not staged", "Run 06_ingest_genetic_data.py from VCF/GDS to populate per-locus H_obs & D"))
  }
  hd_df <- tryCatch(utils::read.delim(hd_path, stringsAsFactors = FALSE), error = function(e) data.frame())
  if (nrow(hd_df) == 0L || !all(c("H_obs", "D_read_ratio") %in% names(hd_df))) {
    return(.empty_panel_06c(t_a, "HDplot diagnostics contain zero H_obs / D_read_ratio rows"))
  }
  lw_path <- file.path(parsed_dir, "locus_weights.tsv")
  if (file.exists(lw_path) && ("locus_id" %in% names(hd_df))) {
    lw_df <- tryCatch(utils::read.delim(lw_path, stringsAsFactors = FALSE), error = function(e) data.frame())
    if (nrow(lw_df) > 0L && ("locus_id" %in% names(lw_df))) {
      m_lw <- match(as.character(hd_df$locus_id), as.character(lw_df$locus_id))
      ok_lw <- !is.na(m_lw)
      if ("locus_paralogy_weight" %in% names(lw_df)) hd_df$locus_paralogy_weight[ok_lw] <- lw_df$locus_paralogy_weight[m_lw[ok_lw]]
      if ("is_paraloss1272_single_copy" %in% names(lw_df)) hd_df$is_paraloss1272_single_copy[ok_lw] <- lw_df$is_paraloss1272_single_copy[m_lw[ok_lw]]
    }
  }
  h_v <- suppressWarnings(as.numeric(hd_df$H_obs))
  d_v <- suppressWarnings(as.numeric(hd_df$D_read_ratio))
  keep_v <- if ("hdplot_keep" %in% names(hd_df)) suppressWarnings(as.integer(hd_df$hdplot_keep)) else rep(1L, nrow(hd_df))
  pl_v <- if ("is_paraloss1272_single_copy" %in% names(hd_df)) suppressWarnings(as.integer(hd_df$is_paraloss1272_single_copy)) else rep(0L, nrow(hd_df))
  wt_v <- if ("locus_paralogy_weight" %in% names(hd_df)) suppressWarnings(as.numeric(hd_df$locus_paralogy_weight)) else rep(1.0, nrow(hd_df))
  fin <- is.finite(h_v) & is.finite(d_v)
  if (!any(fin)) {
    return(.empty_panel_06c(t_a, "No finite H_obs / D_read_ratio locus records"))
  }
  keep_v[is.na(keep_v)] <- 1L
  pl_v[is.na(pl_v)] <- 0L
  is_pl <- (keep_v == 1L) & (pl_v == 1L | (is.finite(wt_v) & wt_v >= 1.25))
  cls_lvls <- c("ParaLoss-1272 Ortholog (w=1.5)", "Retained RAD Locus (w=1.0)", "Filtered Paralog / Sparse")
  cls_v <- ifelse(is_pl, cls_lvls[1L], ifelse(keep_v == 1L, cls_lvls[2L], cls_lvls[3L]))
  n_pl <- sum(fin & is_pl)
  n_ret <- sum(fin & keep_v == 1L)
  n_filt <- sum(fin & keep_v != 1L)

  p_df <- data.frame(H_obs = pmin(1.0, pmax(0.0, h_v[fin])), D_read_ratio = d_v[fin], locus_class = factor(cls_v[fin], levels = cls_lvls), stringsAsFactors = FALSE)
  if (nrow(p_df) > 25000L) {
    set.seed(42L)
    p_df <- p_df[sort(unique(c(which(p_df$locus_class == cls_lvls[1L]), sample(which(p_df$locus_class != cls_lvls[1L]), min(20000L, sum(p_df$locus_class != cls_lvls[1L])))))), , drop = FALSE]
  }
  pal_hd <- stats::setNames(c("#0072B2", "#009E73", "#D55E00"), cls_lvls)
  shp_hd <- stats::setNames(c(1L, 2L, 4L), cls_lvls)
  ggplot2::ggplot(p_df, ggplot2::aes(x = .data$H_obs, y = .data$D_read_ratio, color = .data$locus_class, shape = .data$locus_class)) +
    ggplot2::geom_vline(xintercept = 0.60, linetype = "dashed", color = "#475569", linewidth = 0.45) +
    ggplot2::geom_hline(yintercept = c(-4.0, 4.0), linetype = "dashed", color = "#475569", linewidth = 0.45) +
    ggplot2::geom_point(size = 2.0, stroke = 0.85, fill = NA, alpha = 0.82) +
    ggplot2::scale_color_manual(values = pal_hd, name = "Locus Tier", drop = FALSE) +
    ggplot2::scale_shape_manual(values = shp_hd, name = "Locus Tier", drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = t_a, subtitle = sprintf("Retained = %d/%d loci (ParaLoss-1272 w=1.5: %d; Filtered: %d)", n_ret, sum(fin), n_pl, n_filt), x = "Observed Heterozygosity (H_obs)", y = "Allele Read-Ratio Deviation (D)") +
    theme_packera_pub(base_size = 9.3) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 7.6)) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2L))
}

.find_cytonuclear_csv_06c <- function(parsed_dir, out_dir) {
  cands <- c(file.path(parsed_dir, "cytonuclear_discordance.csv"), file.path(out_dir, "cytonuclear_discordance.csv"))
  mf_path <- if (!is.null(parsed_dir) && nzchar(parsed_dir)) file.path(parsed_dir, "genomic_manifest.json") else ""
  if (nzchar(mf_path) && file.exists(mf_path) && requireNamespace("jsonlite", quietly = TRUE)) {
    mf <- tryCatch(jsonlite::fromJSON(mf_path), error = function(e) list())
    if (!is.null(mf$genomics_dir) && nzchar(as.character(mf$genomics_dir))) {
      cands <- c(cands, list.files(as.character(mf$genomics_dir), pattern = "^cytonuclear_discordance\\.csv$", recursive = TRUE, full.names = TRUE))
    }
  }
  env_g <- Sys.getenv("PLANTMATICS_GENOMICS_DIR", unset = "data/genomics")
  if (nzchar(env_g) && dir.exists(env_g)) cands <- c(cands, list.files(env_g, pattern = "^cytonuclear_discordance\\.csv$", recursive = TRUE, full.names = TRUE))
  for (p in unique(cands)) if (!is.na(p) && nzchar(p) && file.exists(p)) {
    return(p)
  }
  ""
}

.build_panel06c_cytonuclear_gg <- function(parsed_dir, out_dir, full_audit_df, organellar_placement_df) {
  t_b <- "B. Cytonuclear (cpDNA / nrITS) & Nuclear waxy Concordance"
  df <- if (is.data.frame(full_audit_df) && nrow(full_audit_df) > 0L) full_audit_df else data.frame()
  if (nrow(df) == 0L || !("genetic_sample_id" %in% names(df))) {
    return(.empty_panel_06c(t_b, "No passing genomic sample table supplied"))
  }
  sids <- as.character(df$genetic_sample_id)
  pick_col <- function(nm_a, nm_b = NULL) {
    va <- if (nm_a %in% names(df)) .clean_str_06c(df[[nm_a]]) else rep(NA_character_, nrow(df))
    if (!is.null(nm_b) && (nm_b %in% names(df))) {
      vb <- .clean_str_06c(df[[nm_b]])
      va[is.na(va)] <- vb[is.na(va)]
    }
    va
  }
  nuc <- pick_col("geno_predicted_taxon", "species_standardized")
  cp <- pick_col("public_plastome_match", "plastid_clade")
  its <- pick_col("public_its_match", "its_ribotype")
  waxy <- pick_col("waxy_lineage_match")
  disc <- if ("cytonuclear_discordant" %in% names(df)) suppressWarnings(as.integer(df$cytonuclear_discordant)) else rep(0L, nrow(df))
  disc[is.na(disc)] <- 0L

  cy_path <- .find_cytonuclear_csv_06c(parsed_dir, out_dir)
  if (nzchar(cy_path)) {
    cy_df <- tryCatch(utils::read.csv(cy_path, stringsAsFactors = FALSE), error = function(e) data.frame())
    if (nrow(cy_df) > 0L && ("genetic_sample_id" %in% names(cy_df))) {
      mc <- match(sids, as.character(cy_df$genetic_sample_id))
      hit <- !is.na(mc)
      if ("public_plastome_match" %in% names(cy_df)) {
        v <- .clean_str_06c(cy_df$public_plastome_match[mc])
        cp[hit & is.na(cp)] <- v[hit & is.na(cp)]
      }
      if ("plastid_clade" %in% names(cy_df)) {
        v <- .clean_str_06c(cy_df$plastid_clade[mc])
        cp[hit & is.na(cp)] <- v[hit & is.na(cp)]
      }
      if ("public_its_match" %in% names(cy_df)) {
        v <- .clean_str_06c(cy_df$public_its_match[mc])
        its[hit & is.na(its)] <- v[hit & is.na(its)]
      }
      if ("its_ribotype" %in% names(cy_df)) {
        v <- .clean_str_06c(cy_df$its_ribotype[mc])
        its[hit & is.na(its)] <- v[hit & is.na(its)]
      }
      if ("waxy_lineage_match" %in% names(cy_df)) {
        v <- .clean_str_06c(cy_df$waxy_lineage_match[mc])
        waxy[hit & is.na(waxy)] <- v[hit & is.na(waxy)]
      }
      if ("cytonuclear_discordant" %in% names(cy_df)) {
        dv <- suppressWarnings(as.integer(cy_df$cytonuclear_discordant[mc]))
        disc[hit & !is.na(dv) & dv == 1L] <- 1L
      }
    }
  }
  if (is.data.frame(organellar_placement_df) && nrow(organellar_placement_df) > 0L && all(c("genetic_sample_id", "locus", "best_match_taxon") %in% names(organellar_placement_df))) {
    for (spec in list(c("plastid", "cp"), c("its", "its"), c("waxy", "waxy"))) {
      sub_p <- organellar_placement_df[tolower(trimws(organellar_placement_df$locus)) == spec[1L], , drop = FALSE]
      if (nrow(sub_p) > 0L) {
        mp <- match(sids, as.character(sub_p$genetic_sample_id))
        okp <- !is.na(mp)
        vals <- .clean_str_06c(sub_p$best_match_taxon[mp])
        if (spec[2L] == "cp") cp[okp & is.na(cp)] <- vals[okp & is.na(cp)]
        if (spec[2L] == "its") its[okp & is.na(its)] <- vals[okp & is.na(its)]
        if (spec[2L] == "waxy") waxy[okp & is.na(waxy)] <- vals[okp & is.na(waxy)]
      }
    }
  }
  has_loc <- !is.na(cp) | !is.na(its) | !is.na(waxy)
  if (!any(has_loc)) {
    return(.empty_panel_06c(t_b, "No skimmed plastome, nrITS, or nuclear waxy placements detected", "Stage skimmed FASTA references or cytonuclear_discordance.csv in cohort"))
  }
  disc[disc == 0L & !is.na(cp) & !is.na(nuc) & grepl("^Packera ", cp) & cp != nuc] <- 1L

  keep_idx <- which(has_loc)
  if (length(keep_idx) > 26L) keep_idx <- keep_idx[order(-disc[keep_idx], sids[keep_idx])[seq_len(26L)]]
  loc_lvls <- c("Nuclear NextRAD", "Plastome cpDNA", "nrITS", "Nuclear waxy")
  sid_ord <- rev(sids[keep_idx][order(-disc[keep_idx], sids[keep_idx])])
  tile_df <- rbind(
    data.frame(sample = sids[keep_idx], locus = loc_lvls[1L], call = nuc[keep_idx], disc = disc[keep_idx], stringsAsFactors = FALSE),
    data.frame(sample = sids[keep_idx], locus = loc_lvls[2L], call = cp[keep_idx], disc = disc[keep_idx], stringsAsFactors = FALSE),
    data.frame(sample = sids[keep_idx], locus = loc_lvls[3L], call = its[keep_idx], disc = disc[keep_idx], stringsAsFactors = FALSE),
    data.frame(sample = sids[keep_idx], locus = loc_lvls[4L], call = waxy[keep_idx], disc = disc[keep_idx], stringsAsFactors = FALSE)
  )
  tile_df$sample <- factor(tile_df$sample, levels = sid_ord)
  tile_df$locus <- factor(tile_df$locus, levels = loc_lvls)
  c_lvls <- sort(unique(stats::na.omit(tile_df$call)))
  tile_df$call <- factor(tile_df$call, levels = c_lvls)
  pal_cy <- get_combinatorial_palette_r(c_lvls)
  lab_fn <- if (exists("icn_discrete_labeller", mode = "function")) icn_discrete_labeller(abbreviate_genus = TRUE) else ggplot2::waiver()
  disc_df <- tile_df[tile_df$disc == 1L, , drop = FALSE]
  pres_df <- tile_df[!is.na(tile_df$call), , drop = FALSE]

  p_b <- ggplot2::ggplot(tile_df, ggplot2::aes(x = .data$locus, y = .data$sample, fill = .data$call)) +
    ggplot2::geom_tile(color = "#FFFFFF", linewidth = 0.45) +
    ggplot2::geom_point(data = pres_df, shape = 1L, fill = NA, color = "#0F172A", size = 1.9, stroke = 0.75, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal_cy, labels = lab_fn, na.value = "#F1F5F9", name = "Lineage Match", na.translate = FALSE) +
    ggplot2::labs(title = t_b, subtitle = sprintf("Placed n = %d samples | Cytonuclear discordance (red border) = %d", length(keep_idx), sum(disc[keep_idx] == 1L)), x = NULL, y = "Sample ID") +
    theme_packera_pub(base_size = 9.3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, size = 8.0), axis.text.y = ggplot2::element_text(size = 7.2))
  if (nrow(disc_df) > 0L) p_b <- p_b + ggplot2::geom_tile(data = disc_df, fill = NA, color = "#B91C1C", linewidth = 0.95, show.legend = FALSE)
  p_b
}

.build_panel06c_phytop_gg <- function(out_dir, pub_concordance_df) {
  t_c <- "C. Public Backbone PhyTop Conflict (ILS vs. Introgression)"
  cdf <- if (is.data.frame(pub_concordance_df) && nrow(pub_concordance_df) > 0L) {
    pub_concordance_df
  } else {
    cp <- if (!is.null(out_dir) && nzchar(out_dir)) file.path(out_dir, "public_backbone_concordance.csv") else ""
    if (nzchar(cp) && file.exists(cp)) tryCatch(utils::read.csv(cp, stringsAsFactors = FALSE), error = function(e) data.frame()) else data.frame()
  }
  if (nrow(cdf) == 0L || !all(c("public_phytop_ils_index", "public_phytop_ih_index") %in% names(cdf))) {
    return(.empty_panel_06c(t_c, "< 2 taxa matched public PhyTop conflict nodes", "Stage public_packera_* Dryad PhyTop metrics to populate ILS vs. IH coordinates"))
  }
  ils_v <- suppressWarnings(as.numeric(cdf$public_phytop_ils_index))
  ih_v <- suppressWarnings(as.numeric(cdf$public_phytop_ih_index))
  tx_col <- intersect(c("taxon", "normalized_taxon", "species_standardized"), names(cdf))[1L]
  tx_v <- if (!is.na(tx_col)) .clean_str_06c(cdf[[tx_col]]) else rep("Unknown", nrow(cdf))
  fin <- is.finite(ils_v) & is.finite(ih_v) & !is.na(tx_v)
  if (sum(fin) < 2L) {
    return(.empty_panel_06c(t_c, sprintf("Only %d taxon matched public PhyTop conflict nodes (< 2 required)", sum(fin)), "Download public Packera backbone PhyTop tables for multi-taxon conflict overlay"))
  }

  pt_df <- stats::aggregate(cbind(ils, ih) ~ taxon, data = data.frame(taxon = tx_v[fin], ils = ils_v[fin], ih = ih_v[fin], stringsAsFactors = FALSE), FUN = mean)
  pt_df$lbl <- abbreviate_icn_taxon(pt_df$taxon, abbreviate_genus = TRUE)
  p_c <- ggplot2::ggplot(pt_df, ggplot2::aes(x = .data$ils, y = .data$ih, color = .data$taxon, shape = .data$taxon)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#94A3B8", linewidth = 0.45) +
    ggplot2::geom_point(size = 3.1, stroke = 1.05, fill = NA) +
    ggplot2::scale_color_manual(values = get_combinatorial_palette_r(pt_df$taxon)) +
    ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(pt_df$taxon)) +
    ggplot2::labs(title = t_c, subtitle = sprintf("Moore-Pollard et al. (2026) backbone conflict across n = %d taxa", nrow(pt_df)), x = "PhyTop ILS Concordance Index", y = "PhyTop Introgression Index (IH)") +
    theme_packera_pub(base_size = 9.3) +
    ggplot2::theme(legend.position = "none")
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p_c + ggrepel::geom_text_repel(ggplot2::aes(label = .data$lbl), fontface = "italic", size = 2.8, box.padding = 0.35, point.padding = 0.25, max.overlaps = Inf, seed = 42L, show.legend = FALSE)
  } else {
    p_c + ggplot2::geom_text(ggplot2::aes(label = .data$lbl), fontface = "italic", size = 2.8, vjust = -0.7, show.legend = FALSE)
  }
}

.build_panel06c_cytotype_gg <- function(full_audit_df, db_path) {
  t_d <- "D. Cytotype Inference (nQuire Delta-LogL & CCDB/ICNA Priors)"
  df <- if (is.data.frame(full_audit_df) && nrow(full_audit_df) > 0L) full_audit_df else data.frame()
  if (nrow(df) == 0L) {
    return(.empty_panel_06c(t_d, "No sample rows supplied"))
  }
  tx_d <- if ("species_standardized" %in% names(df)) .clean_str_06c(df$species_standardized) else rep(NA_character_, nrow(df))
  tx_g <- if ("geno_predicted_taxon" %in% names(df)) .clean_str_06c(df$geno_predicted_taxon) else rep(NA_character_, nrow(df))
  tx_d[is.na(tx_d)] <- tx_g[is.na(tx_d)]
  tx_d[is.na(tx_d)] <- "Unassigned"
  pl_col <- if ("geno_ploidy_level" %in% names(df)) df$geno_ploidy_level else if ("ploidy_level" %in% names(df)) df$ploidy_level else rep(NA_real_, nrow(df))
  pl_num <- suppressWarnings(as.integer(sub("x$", "", tolower(trimws(as.character(pl_col))))))
  cy_raw <- if ("inferred_cytotype" %in% names(df)) .clean_str_06c(df$inferred_cytotype) else rep(NA_character_, nrow(df))
  cy_raw[is.na(cy_raw) & is.finite(pl_num) & pl_num >= 2L] <- paste0(pl_num[is.na(cy_raw) & is.finite(pl_num) & pl_num >= 2L], "x")
  dlogl <- if ("nquire_delta_logl" %in% names(df)) suppressWarnings(as.numeric(df$nquire_delta_logl)) else rep(NA_real_, nrow(df))
  ohet <- if ("observed_het" %in% names(df)) suppressWarnings(as.numeric(df$observed_het)) else rep(NA_real_, nrow(df))
  use_dlogl <- any(is.finite(dlogl))
  y_v <- if (use_dlogl) dlogl else ohet
  ok <- is.finite(y_v) & !is.na(cy_raw)
  if (!any(ok)) {
    return(.empty_panel_06c(t_d, "No finite nQuire delta-LogL or observed heterozygosity cytotype calls"))
  }

  prior_sub <- "Empirical Bayes allele-balance screen"
  if (!is.null(db_path) && nzchar(db_path) && db_path != ":memory:" && file.exists(db_path) && requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {
    con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), dbname = db_path, flags = RSQLite::SQLITE_RO), error = function(e) NULL)
    if (!is.null(con)) {
      on.exit(tryCatch(DBI::dbDisconnect(con), error = function(e) NULL), add = TRUE)
      if (isTRUE("taxon_cytology_summary" %in% DBI::dbListTables(con))) {
        cs <- tryCatch(DBI::dbGetQuery(con, "SELECT species_standardized, reported_cytotypes, modal_sporophytic_2n FROM taxon_cytology_summary"), error = function(e) data.frame())
        cs_hit <- cs[cs$species_standardized %in% unique(tx_d[ok]), , drop = FALSE]
        if (nrow(cs_hit) > 0L) {
          m2n <- sort(unique(stats::na.omit(suppressWarnings(as.integer(cs_hit$modal_sporophytic_2n)))))
          prior_sub <- sprintf("CCDB/ICNA SQLite priors active (%d taxa%s)", nrow(cs_hit), if (length(m2n) > 0L) sprintf("; modal 2n=%s", paste(m2n, collapse = ",")) else "")
        }
      }
    }
  }
  cy_std <- c("2x", "3x", "4x", "5x", "6x", "8x")
  cy_lvls <- c(intersect(cy_std, unique(cy_raw[ok])), setdiff(sort(unique(cy_raw[ok])), cy_std))
  c_df <- data.frame(taxon = tx_d[ok], cyto_f = factor(cy_raw[ok], levels = cy_lvls), y_val = y_v[ok], stringsAsFactors = FALSE)
  med_df <- stats::aggregate(y_val ~ cyto_f, data = c_df, FUN = stats::median)
  med_df$x_num <- as.numeric(med_df$cyto_f)
  lab_fn <- if (exists("icn_discrete_labeller", mode = "function")) icn_discrete_labeller(abbreviate_genus = TRUE) else ggplot2::waiver()

  p_d <- ggplot2::ggplot(c_df, ggplot2::aes(x = .data$cyto_f, y = .data$y_val, color = .data$taxon, shape = .data$taxon)) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.14, height = 0, seed = 42L), size = 2.6, stroke = 0.95, fill = NA) +
    ggplot2::geom_segment(data = med_df, ggplot2::aes(x = .data$x_num - 0.22, xend = .data$x_num + 0.22, y = .data$y_val, yend = .data$y_val), inherit.aes = FALSE, color = "#1E293B", linewidth = 0.75) +
    ggplot2::scale_color_manual(values = get_combinatorial_palette_r(c_df$taxon), labels = lab_fn, name = "Taxon") +
    ggplot2::scale_shape_manual(values = get_combinatorial_hollow_shapes_r(c_df$taxon), labels = lab_fn, name = "Taxon") +
    ggplot2::labs(title = t_d, subtitle = prior_sub, x = "Inferred Cytotype", y = if (use_dlogl) "nQuire Delta Log-Likelihood" else "Observed Heterozygosity (H_obs)") +
    theme_packera_pub(base_size = 9.3) +
    ggplot2::theme(legend.position = "right")
  if (use_dlogl) p_d <- p_d + ggplot2::geom_hline(yintercept = if (any(c_df$y_val < 0)) c(-4.0, 4.0) else 4.0, linetype = "dashed", color = "#475569", linewidth = 0.45)
  p_d
}

render_fig06c_cytonuclear_paralog_cytotype_png <- function(
  out_png, parsed_dir, out_dir, full_audit_df, pub_concordance_df = NULL, organellar_placement_df = NULL, db_path = NULL
) {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") || identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) {
    return(invisible(""))
  }
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    return(invisible(""))
  }
  combo <- (.build_panel06c_hdplot_gg(parsed_dir) | .build_panel06c_cytonuclear_gg(parsed_dir, out_dir, full_audit_df, organellar_placement_df)) /
    (.build_panel06c_phytop_gg(out_dir, pub_concordance_df) | .build_panel06c_cytotype_gg(full_audit_df, db_path)) +
    patchwork::plot_annotation(
      title = "Step 06: Cytonuclear Discordance, Paralog Filtering & Cytotype Audit",
      subtitle = sprintf("Audit across n = %d genomic samples | HDplot / ParaLoss-1272, cpDNA/ITS/waxy, PhyTop & nQuire/CCDB", if (is.data.frame(full_audit_df)) nrow(full_audit_df) else 0L),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11.8), plot.subtitle = ggplot2::element_text(size = 9.0, color = "#475569"))
    )
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  save_packera_png_300(combo, out_png, width = 15.0, height = 10.5, dpi = 300L)
  invisible(out_png)
}
