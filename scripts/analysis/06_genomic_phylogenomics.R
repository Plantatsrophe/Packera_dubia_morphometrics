#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/06_genomic_phylogenomics.R
# Project: Plantmatics
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Standalone Step 06 Nuclear & Organellar Phylogenomics Stream.
#   1. Executes as a clean zero-overhead no-op when `data/genomics/` holds no payloads.
#   2. Computes pairwise allele-sharing distances, M/P replicate concordance,
#      PCoA, PAM clustering, LORO-CV shrinkage DAPC, convex admixture simplex,
#      122 bp microhaplotype co-ancestry, (h, Q12) hybrid classifier, tree-guided
#      f-branch, SNPs2CF/SNaQ quartets, and invariant-site SVDquartets/BPP/gdi.
#   3. Persists Tier B SQLite `genetic_lineage_assignments` (`20261127`) BEFORE
#      running post-hoc MMRR/FEEMS/Mantel cross-modal congruence.
#   4. Reconciles empirical cytonuclear discordance and renders 300-DPI PNG plates.
# ==============================================================================

suppressPackageStartupMessages({
  if (requireNamespace("ape", quietly = TRUE)) library(ape)
  if (requireNamespace("cluster", quietly = TRUE)) library(cluster)
  if (requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)
  if (requireNamespace("patchwork", quietly = TRUE)) library(patchwork)
  if (requireNamespace("dplyr", quietly = TRUE)) library(dplyr)
  if (requireNamespace("readr", quietly = TRUE)) library(readr)
  if (requireNamespace("optparse", quietly = TRUE)) library(optparse)
})

if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

.src_step06_rel <- function(rel_path) {
  cands <- c(
    rel_path,
    file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), basename(rel_path)),
    file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "..", "..", rel_path)
  )
  for (p in cands) {
    if (!is.na(p) && nzchar(p) && file.exists(p)) {
      source(p)
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

if (!exists("persist_step06_genomics_to_sqlite", mode = "function")) .src_step06_rel("scripts/core/db_utils.R")
if (!exists("save_packera_png_300", mode = "function")) .src_step06_rel("scripts/figures/theme_packera_pub.R")
if (!exists(".compute_allele_sharing_distance", mode = "function")) .src_step06_rel("scripts/analysis/step06_phylogenomics_utils.R")
if (!exists("compute_microhaplotype_coancestry", mode = "function")) .src_step06_rel("scripts/analysis/step06_reticulation_coancestry.R")
if (!exists("run_step06_landscape_genomics", mode = "function")) .src_step06_rel("scripts/analysis/step06_spatial_landscape_genomics.R")
if (!exists("load_public_packera_backbone", mode = "function")) .src_step06_rel("scripts/analysis/step06_public_backbone.R")
if (!exists("select_admixture_k_cross_entropy", mode = "function")) .src_step06_rel("scripts/analysis/step06_model_selection.R")
if (!exists("render_fig06a_genomic_pcoa_tree_png", mode = "function")) .src_step06_rel("scripts/figures/step06_figures.R")
if (!exists("render_fig06c_cytonuclear_paralog_cytotype_png", mode = "function")) .src_step06_rel("scripts/figures/step06c_figures.R")

.resolve_step06_db_env <- function(default = "") {
  for (k in c("PLANTMATICS_DB_PATH", "PACKERA_SQLITE_DB", "PACKERA_DB_PATH")) {
    v <- Sys.getenv(k, unset = "")
    if (nzchar(trimws(v))) return(trimws(v))
  }
  default
}
.emit_s06_prog <- function(...) {
  if (exists("emit_r_progress", mode = "function")) tryCatch(emit_r_progress(6L, ...), error = function(e) NULL)
}

parse_step06_args <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  parser <- optparse::OptionParser()
  parser <- optparse::add_option(parser, c("--genomics-dir"), type = "character", default = "data/genomics")
  parser <- optparse::add_option(parser, c("--raw-genomics-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--parsed-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--genomic-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--out-dir"), type = "character", default = "outputs/genomics")
  parser <- optparse::add_option(parser, c("--tables-dir"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--output-tables"), type = "character", default = "outputs/tables")
  parser <- optparse::add_option(parser, c("--output-figures"), type = "character", default = "outputs/figures")
  parser <- optparse::add_option(parser, c("--figure-a"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--figure-b"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--figure-c"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--run-id"), type = "character", default = Sys.getenv("RUN_ID", unset = "canonical"))
  parser <- optparse::add_option(parser, c("--db"), type = "character", default = "")
  parser <- optparse::add_option(parser, c("--db-path"), type = "character", default = .resolve_step06_db_env(""))
  parser <- optparse::add_option(parser, c("--demo"), action = "store_true", default = FALSE)
  optparse::parse_args(parser, args = cli_args)
}

.resolve_first_nonempty <- function(...) {
  for (v in list(...)) if (!is.null(v) && nzchar(trimws(as.character(v)))) return(trimws(as.character(v)))
  ""
}

run_step06_genomic_phylogenomics <- function(opts) {
  raw_geno_dir <- .resolve_first_nonempty(opts[["raw-genomics-dir"]], opts[["genomics-dir"]], "data/genomics")
  parsed_dir <- .resolve_first_nonempty(opts[["parsed-dir"]], opts[["genomic-dir"]], "outputs/reports")
  out_dir <- .resolve_first_nonempty(opts[["out-dir"]], "outputs/genomics")
  tables_dir <- .resolve_first_nonempty(opts[["tables-dir"]], opts[["output-tables"]], "")
  fig_dir <- .resolve_first_nonempty(opts[["output-figures"]], "outputs/figures")
  eff_db_in <- .resolve_first_nonempty(opts[["db"]], opts[["db-path"]], "")

  samples_csv <- file.path(parsed_dir, "parsed_genomic_samples.csv")
  dosage_tsv <- file.path(parsed_dir, "allele_dosage_matrix.tsv")
  manifest_json <- file.path(parsed_dir, "genomic_manifest.json")
  if (!file.exists(samples_csv) || !file.exists(dosage_tsv)) {
    message("[Step 06/07] No genomic files staged in ", parsed_dir, " - skipping standalone Step 06")
    .emit_s06_prog("genomic_phylogenomics", 0L, 0L, unit = "samples", status = "skipped_cached")
    return(invisible(NULL))
  }

  samples_df <- utils::read.csv(samples_csv, stringsAsFactors = FALSE)
  pass_df <- samples_df[!is.na(samples_df$geno_quality_pass) & as.integer(samples_df$geno_quality_pass) == 1L, , drop = FALSE]
  if (nrow(pass_df) < 3L) {
    message("[Step 06/07] Fewer than 3 QC-passing genomic samples - skipping standalone Step 06")
    .emit_s06_prog("genomic_phylogenomics", 0L, 0L, unit = "samples", status = "skipped_cached")
    return(invisible(NULL))
  }

  dos_df <- readr::read_tsv(dosage_tsv, show_col_types = FALSE, progress = FALSE)
  dos_df$genetic_sample_id <- as.character(dos_df$genetic_sample_id)
  common_ids <- intersect(pass_df$genetic_sample_id, dos_df$genetic_sample_id)
  pass_df <- pass_df[match(common_ids, pass_df$genetic_sample_id), , drop = FALSE]
  dos_mat <- as.matrix(dos_df[match(common_ids, dos_df$genetic_sample_id), setdiff(names(dos_df), "genetic_sample_id"), drop = FALSE])
  storage.mode(dos_mat) <- "numeric"; rownames(dos_mat) <- common_ids

  dist_res <- .compute_allele_sharing_distance(dos_mat); d_mat <- dist_res$dist; d_raw <- dist_res$raw_dist
  rep_conc <- rep(NA_real_, nrow(pass_df)); rep_ids <- as.character(pass_df$replicate_pair_id)
  for (i in seq_len(nrow(pass_df))) {
    if (!is.na(rep_ids[i]) && nzchar(rep_ids[i])) {
      mates <- setdiff(which(rep_ids == rep_ids[i]), i)
      if (length(mates) > 0L && any(is.finite(d_raw[i, mates]))) rep_conc[i] <- max(0.0, min(1.0, 1.0 - mean(d_raw[i, mates], na.rm = TRUE)))
    }
  }

  pco_fit <- stats::cmdscale(stats::as.dist(d_mat), k = min(8L, max(2L, nrow(pass_df) - 1L)), eig = TRUE)
  eig_pos <- pmax(0.0, pco_fit$eig)
  eig_pct <- if (sum(eig_pos) > 0) 100.0 * eig_pos / sum(eig_pos) else rep(100.0 / length(eig_pos), length(eig_pos))
  k_use <- min(6L, ncol(pco_fit$points), max(2L, sum(eig_pos > 1e-9)))
  pco_coords <- pco_fit$points[, seq_len(k_use), drop = FALSE]; colnames(pco_coords) <- paste0("PCo", seq_len(k_use))
  pass_df$PCo1 <- pco_coords[, 1L]; pass_df$PCo2 <- if (k_use >= 2L) pco_coords[, 2L] else 0.0

  ce_res <- select_admixture_k_cross_entropy(dos_mat, d_mat = d_mat, k_max = min(8L, max(1L, nrow(dos_mat) - 1L)), n_folds = 5L, seed = 202611L)
  best_k <- as.integer(ce_res$best_k)
  sel_ce_row <- ce_res$table_df[ce_res$table_df$K == best_k, , drop = FALSE]
  best_sil <- if (nrow(sel_ce_row) > 0L && is.finite(sel_ce_row$pam_silhouette[1L])) as.numeric(sel_ce_row$pam_silhouette[1L]) else NA_real_
  best_clus <- rep(1L, nrow(pass_df))
  if (best_k >= 2L && nrow(pass_df) > best_k) {
    pam_fit <- if (requireNamespace("cluster", quietly = TRUE)) tryCatch(cluster::pam(stats::as.dist(d_mat), k = best_k, diss = TRUE), error = function(e) NULL) else NULL
    best_clus <- if (!is.null(pam_fit)) pam_fit$clustering else max.col(ce_res$Q_opt, ties.method = "first")
  }
  gen_cluster <- sprintf("Geno_Cluster_%d", as.integer(best_clus))
  .emit_s06_prog("pcoa_admixture_k", 1L, 4L, unit = "stages", substage_weight = 0.25, substage_offset = 0.0)

  is_ref_anchor <- (is.na(pass_df$is_putative_novel_or_hybrid_query) | as.integer(pass_df$is_putative_novel_or_hybrid_query) == 0L) &
    !is.na(pass_df$species_standardized) & nzchar(pass_df$species_standardized)
  ref_idx <- which(is_ref_anchor); ref_tbl <- table(pass_df$species_standardized[ref_idx]); keep_taxa <- names(ref_tbl[ref_tbl >= 1L])
  if (length(keep_taxa) >= 2L) ref_idx <- ref_idx[pass_df$species_standardized[ref_idx] %in% keep_taxa] else keep_taxa <- names(ref_tbl)

  m_dapc <- min(5L, ncol(pco_coords), max(2L, length(ref_idx) - length(keep_taxa)))
  X_dapc <- pco_coords[, seq_len(m_dapc), drop = FALSE]; y_ref <- pass_df$species_standardized[ref_idx]
  full_model <- .fit_shrinkage_dapc(X_dapc[ref_idx, , drop = FALSE], y_ref, X_dapc)
  pred_taxon <- full_model$pred; post_prob <- full_model$max_p; margin_vec <- full_model$margin; probs_mat <- full_model$probs
  oof_pred <- rep(NA_character_, length(ref_idx)); oof_logloss <- rep(NA_real_, length(ref_idx))
  for (ri in seq_along(ref_idx)) {
    i_global <- ref_idx[ri]; rep_i <- rep_ids[i_global]
    hold_mask <- if (!is.na(rep_i) && nzchar(rep_i)) (rep_ids[ref_idx] == rep_i) else (seq_along(ref_idx) == ri)
    tr_sub <- which(!hold_mask)
    if (sum(y_ref[tr_sub] == y_ref[i_global], na.rm = TRUE) >= 1L && length(unique(y_ref[tr_sub])) >= 2L) {
      cv_fit <- .fit_shrinkage_dapc(X_dapc[ref_idx[tr_sub], , drop = FALSE], y_ref[tr_sub], X_dapc[i_global, , drop = FALSE])
      oof_pred[ri] <- cv_fit$pred[1L]; pred_taxon[i_global] <- cv_fit$pred[1L]; post_prob[i_global] <- cv_fit$max_p[1L]; margin_vec[i_global] <- cv_fit$margin[1L]
      probs_mat[i_global, ] <- 0.0; common_cols <- intersect(colnames(probs_mat), colnames(cv_fit$probs))
      if (length(common_cols) > 0L) probs_mat[i_global, common_cols] <- cv_fit$probs[1L, common_cols]
      p_true <- if (as.character(y_ref[i_global]) %in% colnames(cv_fit$probs)) as.numeric(cv_fit$probs[1L, as.character(y_ref[i_global])]) else 1e-6
      oof_logloss[ri] <- -log(min(1.0 - 1e-6, max(1e-6, p_true)))
    }
  }
  valid_oof <- !is.na(oof_pred)
  loro_acc <- if (any(valid_oof)) mean(oof_pred[valid_oof] == y_ref[valid_oof], na.rm = TRUE) else NA_real_
  loro_kappa <- if (sum(valid_oof) >= 2L) .compute_kappa_r(y_ref[valid_oof], oof_pred[valid_oof]) else NA_real_
  loro_logloss <- if (any(is.finite(oof_logloss))) mean(oof_logloss, na.rm = TRUE) else NA_real_
  .emit_s06_prog("shrinkage_dapc", 2L, 4L, unit = "stages", substage_weight = 0.25, substage_offset = 0.25)

  tree_files <- character(0); q_files <- character(0)
  if (file.exists(manifest_json)) {
    m_txt <- gsub("[\r\n]+", " ", paste(readLines(manifest_json, warn = FALSE), collapse = " "))
    ext_paths <- function(key) {
      m <- regexec(sprintf('"%s"[[:space:]]*:[[:space:]]*\\[(.*?)\\]', key), m_txt); hit <- regmatches(m_txt, m)[[1L]]
      if (length(hit) < 2L) character(0) else gsub('^"|"$', "", regmatches(hit[2L], gregexpr('"([^"]+)"', hit[2L]))[[1L]])
    }
    tree_files <- ext_paths("tree_files"); q_files <- ext_paths("admixture_q_files")
  }
  if (length(tree_files) > 0L) tree_files <- tree_files[!grepl("^public_packera_", basename(as.character(tree_files)), ignore.case = TRUE)]
  if (length(tree_files) == 0L && dir.exists(raw_geno_dir)) {
    cand_tf <- sort(list.files(raw_geno_dir, full.names = TRUE, recursive = TRUE)); b_tf <- basename(cand_tf)
    tree_files <- cand_tf[grepl("^RAxML_bestTree|\\.(treefile|tre|nwk|tree)$", b_tf, ignore.case = TRUE) & !grepl("^public_packera_", b_tf, ignore.case = TRUE)]
  }

  q_ext <- .read_admixture_q_if_present(q_files, common_ids, keep_taxa, probs_mat)
  q_mat <- if (!is.null(q_ext)) q_ext else {
    q_out <- ce_res$Q_opt
    if (length(unique(y_ref[!is.na(y_ref) & nzchar(y_ref)])) >= 2L && nrow(full_model$centroids) >= 2L) {
      q_sup <- .solve_admixture_simplex(X_dapc, full_model$centroids, rownames(full_model$centroids)); rownames(q_sup) <- common_ids
      ce_sup <- evaluate_q_reconstruction_ce(dos_mat, q_sup, n_folds = 5L, seed = 202611L); ce_opt <- ce_res$selected_k_cross_entropy
      if (ncol(q_out) == 1L || (is.finite(ce_sup) && length(ce_opt) == 1L && is.finite(ce_opt) && ce_sup <= 1.05 * ce_opt)) q_out <- q_sup
    }
    q_out
  }
  if (ncol(q_mat) >= 2L && length(ref_idx) >= 2L && length(keep_taxa) >= 2L && (is.null(colnames(q_mat)) || !all(colnames(q_mat) %in% keep_taxa))) {
    if (is.null(colnames(q_mat))) colnames(q_mat) <- paste0("Ancestry_K", seq_len(ncol(q_mat)))
    M_work <- matrix(0.0, nrow = length(keep_taxa), ncol = ncol(q_mat), dimnames = list(keep_taxa, colnames(q_mat)))
    for (ti in seq_along(keep_taxa)) {
      t_rows <- ref_idx[y_ref == keep_taxa[ti]]
      if (length(t_rows) > 0L) M_work[ti, ] <- colMeans(q_mat[t_rows, , drop = FALSE], na.rm = TRUE)
    }
    M_work[!is.finite(M_work)] <- -1.0
    for (step_m in seq_len(min(nrow(M_work), ncol(M_work)))) {
      mx_v <- max(M_work); if (!is.finite(mx_v) || mx_v < 0.50) break
      rc <- which(M_work == mx_v, arr.ind = TRUE)[1L, ]
      colnames(q_mat)[rc[2L]] <- keep_taxa[rc[1L]]; M_work[rc[1L], ] <- -1.0; M_work[, rc[2L]] <- -1.0
    }
  }
  q_cols <- colnames(q_mat); q_max <- apply(q_mat, 1L, max)
  is_admixed <- integer(nrow(pass_df)); p1 <- character(nrow(pass_df)); p2 <- character(nrow(pass_df)); q_json <- character(nrow(pass_df))
  for (i in seq_len(nrow(pass_df))) {
    wi <- q_mat[i, ]; ord_w <- order(wi, decreasing = TRUE); nz <- which(wi >= 0.02)
    q_json[i] <- paste0("{", paste(vapply(nz[order(wi[nz], decreasing = TRUE)], function(j) .json_escape_kv(q_cols[j], wi[j]), character(1)), collapse = ", "), "}")
    p1[i] <- q_cols[ord_w[1L]]; p2[i] <- if (length(ord_w) >= 2L && wi[ord_w[2L]] >= 0.15) q_cols[ord_w[2L]] else NA_character_
    is_admixed[i] <- if (q_max[i] < 0.80 && !is.na(p2[i])) 1L else 0L
  }
  phylo_res <- .extract_phylo_metrics(tree_files, pass_df)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  phylo_res$clade_summary <- evaluate_ml_tree_model_and_support(raw_geno_dir, phylo_res$tree, pass_df, out_dir)
  pub_bb <- load_public_packera_backbone(raw_geno_dir)
  coanc_res <- compute_microhaplotype_coancestry(common_ids, dos_mat, raw_geno_dir)
  tri_res <- compute_hybrid_index_q12(pass_df, dos_mat, q_mat, p1, p2, d_mat)
  fb_res <- compute_tree_guided_fbranch(pass_df, dos_mat, raw_geno_dir, backbone = pub_bb)
  cf_res <- compute_quartet_concordance_snps2cf(pass_df, dos_mat)
  pub_eval <- evaluate_public_backbone_and_organellar(pass_df, raw_geno_dir, out_dir)
  is_admixed[tri_res$hybrid_class == "Divergent_NonHybrid_Lineage"] <- 0L
  .emit_s06_prog("reticulation_backbone", 3L, 4L, unit = "stages", substage_weight = 0.25, substage_offset = 0.50)

  py_bin <- if (exists("resolve_r_python_bin", mode = "function")) resolve_r_python_bin() else if (file.exists(".venv/bin/python")) ".venv/bin/python" else Sys.which("python3")
  coal_script <- "scripts/analysis/step06_coalescent_delimitation.py"
  if (nzchar(py_bin) && file.exists(coal_script)) {
    c_args <- c(shQuote(coal_script), "--genomics-dir", shQuote(raw_geno_dir), "--out-dir", shQuote(out_dir), "--parsed-samples", shQuote(samples_csv), "--dosage-tsv", shQuote(dosage_tsv), "--n-boot", "250")
    tryCatch(system2(py_bin, args = c_args, stdout = TRUE, stderr = TRUE), error = function(e) warning("[Step 06 Coalescent] ", e$message))
  }
  gdi_samp_csv <- file.path(out_dir, "coalescent_gdi_samples.csv"); gdi_tax_csv <- file.path(out_dir, "coalescent_gdi_delimitation.csv")
  gdi_df <- if (file.exists(gdi_tax_csv)) utils::read.csv(gdi_tax_csv, stringsAsFactors = FALSE) else data.frame()
  gdi_s_df <- if (file.exists(gdi_samp_csv)) utils::read.csv(gdi_samp_csv, stringsAsFactors = FALSE) else data.frame()
  gdi_vec <- rep(NA_real_, nrow(pass_df)); bf_vec <- rep(NA_real_, nrow(pass_df))
  if (nrow(gdi_s_df) > 0L && "genetic_sample_id" %in% names(gdi_s_df)) {
    gm <- match(pass_df$genetic_sample_id, gdi_s_df$genetic_sample_id)
    if ("coalescent_gdi" %in% names(gdi_s_df)) gdi_vec <- suppressWarnings(as.numeric(gdi_s_df$coalescent_gdi[gm]))
    if ("bpp_introgression_bf" %in% names(gdi_s_df)) bf_vec <- suppressWarnings(as.numeric(gdi_s_df$bpp_introgression_bf[gm]))
  }

  assign_df <- data.frame(
    genetic_sample_id = pass_df$genetic_sample_id,
    catalogNumber = ifelse(is.na(pass_df$catalogNumber) | !nzchar(pass_df$catalogNumber), NA_character_, pass_df$catalogNumber),
    population_id = pass_df$population_id, geno_predicted_taxon = pred_taxon, geno_posterior_prob = round(post_prob, 4),
    geno_margin = round(margin_vec, 4), genetic_cluster = gen_cluster, q_max = round(q_max, 4), q_admixture_json = q_json,
    is_geno_admixed = is_admixed, geno_parent_1 = p1, geno_parent_2 = p2,
    nearest_tree_sister_taxon = phylo_res$sample_metrics$nearest_tree_sister_taxon,
    patristic_dist_to_dubia = round(phylo_res$sample_metrics$patristic_dist_to_dubia, 6),
    replicate_pair_concordance = round(rep_conc, 4), hybrid_index_h = tri_res$hybrid_index_h,
    interspecific_ancestry_q12 = tri_res$q12, hybrid_generation_class = tri_res$hybrid_class,
    coancestry_cluster = coanc_res$cluster, f_branch_max = fb_res$f_branch_max, f_branch_donor_taxon = fb_res$f_branch_donor_taxon,
    snaq_gamma_minor = cf_res$sample_gamma, coalescent_gdi = gdi_vec, bpp_introgression_bf = bf_vec,
    mmrr_beta_ibd = NA_real_, mmrr_beta_ibe = NA_real_, mmrr_p_ibe = NA_real_,
    public_backbone_clade = pub_eval$public_backbone_clade, public_phytop_ils_index = pub_eval$public_phytop_ils_index,
    public_phytop_ih_index = pub_eval$public_phytop_ih_index, stringsAsFactors = FALSE
  )

  if (!isTRUE(opts$demo) && exists("persist_step06_genomics_to_sqlite", mode = "function")) {
    tryCatch(persist_step06_genomics_to_sqlite(assign_df, run_id = opts[["run-id"]], db_path = if (nzchar(eff_db_in)) eff_db_in else NULL, demo = FALSE),
             error = function(e) warning("[Step 06] SQLite persistence warning: ", e$message))
  }
  eff_db_ro <- if (nzchar(eff_db_in)) eff_db_in else .resolve_step06_db_env("data/plantmatics.sqlite")
  land_res <- if ("out_dir" %in% names(formals(run_step06_landscape_genomics))) {
    run_step06_landscape_genomics(pass_df, d_mat, db_path = eff_db_ro, run_id = opts[["run-id"]], demo = opts$demo, out_dir = out_dir)
  } else run_step06_landscape_genomics(pass_df, d_mat, db_path = eff_db_ro, run_id = opts[["run-id"]], demo = opts$demo)
  if (nrow(land_res$summary_df) > 0L) {
    assign_df$mmrr_beta_ibd <- land_res$summary_df$mmrr_beta_ibd[1L]; assign_df$mmrr_beta_ibe <- land_res$summary_df$mmrr_beta_ibe[1L]; assign_df$mmrr_p_ibe <- land_res$summary_df$mmrr_p_ibe[1L]
    if (!isTRUE(opts$demo) && exists("persist_step06_genomics_to_sqlite", mode = "function")) {
      tryCatch(persist_step06_genomics_to_sqlite(assign_df, run_id = opts[["run-id"]], db_path = if (nzchar(eff_db_in)) eff_db_in else NULL, demo = FALSE), error = function(e) NULL)
    }
  }
  congr_df <- .compute_posthoc_congruence(pass_df, d_mat, pco_coords, eff_db_ro, opts[["run-id"]], opts$demo)

  full_audit <- dplyr::left_join(samples_df, assign_df[, setdiff(names(assign_df), c("catalogNumber", "population_id")), drop = FALSE], by = "genetic_sample_id")
  m_pass <- match(full_audit$genetic_sample_id, pass_df$genetic_sample_id)
  for (mcol in c("public_plastome_match", "public_its_match", "waxy_lineage_match")) {
    if (!(mcol %in% names(full_audit))) full_audit[[mcol]] <- NA_character_
    need_fill <- !is.na(m_pass) & (is.na(full_audit[[mcol]]) | !nzchar(as.character(full_audit[[mcol]])))
    if (any(need_fill)) full_audit[[mcol]][need_fill] <- pub_eval[[mcol]][m_pass[need_fill]]
  }
  full_audit <- reconcile_step06_cytonuclear_discordance(full_audit, assign_df = assign_df, placement_df = pub_eval$placement_df)
  if (!isTRUE(opts$demo) && nzchar(eff_db_ro) && file.exists(eff_db_ro) && "cytonuclear_discordant" %in% names(full_audit)) {
    conn_up <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), dbname = eff_db_ro), error = function(e) NULL)
    if (!is.null(conn_up)) {
      tryCatch({
        if (DBI::dbExistsTable(conn_up, "genomic_samples")) {
          up_rows <- full_audit[!is.na(full_audit$cytonuclear_discordant), c("genetic_sample_id", "cytonuclear_discordant"), drop = FALSE]
          for (ri in seq_len(nrow(up_rows))) DBI::dbExecute(conn_up, "UPDATE genomic_samples SET cytonuclear_discordant = ? WHERE genetic_sample_id = ?", params = list(as.integer(up_rows$cytonuclear_discordant[ri]), as.character(up_rows$genetic_sample_id[ri])))
        }
      }, error = function(e) NULL)
      tryCatch(DBI::dbDisconnect(conn_up), error = function(e) NULL)
    }
  }

  metrics_df <- data.frame(
    n_samples_total = nrow(samples_df), n_samples_pass = nrow(pass_df), n_reference_anchors = length(ref_idx),
    n_novel_queries = sum(as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L, na.rm = TRUE), n_loci_retained = ncol(dos_mat),
    optimal_k_pam = best_k, pam_silhouette = if (is.finite(best_sil)) round(best_sil, 4) else NA_real_,
    selected_k_cross_entropy = if (length(ce_res$selected_k_cross_entropy) > 0L && is.finite(ce_res$selected_k_cross_entropy)) round(as.numeric(ce_res$selected_k_cross_entropy), 6) else NA_real_,
    selected_shrinkage_lambda = round(as.numeric(full_model$lambda_shrink), 4),
    loro_cv_accuracy = round(loro_acc, 4), loro_cv_kappa = round(loro_kappa, 4), loro_cv_logloss = round(loro_logloss, 4),
    mean_replicate_concordance = round(mean(rep_conc, na.rm = TRUE), 4), stringsAsFactors = FALSE
  )

  readr::write_csv(full_audit, file.path(out_dir, "genomic_assignments.csv"), na = "")
  readr::write_csv(phylo_res$clade_summary, file.path(out_dir, "genomic_clade_cohesion_summary.csv"), na = "")
  readr::write_csv(metrics_df, file.path(out_dir, "genomic_model_metrics.csv"), na = "")
  readr::write_csv(ce_res$table_df, file.path(out_dir, "admixture_k_model_selection.csv"), na = "")
  readr::write_csv(congr_df, file.path(out_dir, "genomic_multimodal_congruence_summary.csv"), na = "")
  readr::write_csv(coanc_res$coanc_df, file.path(out_dir, "coancestry_matrix_finerad.csv"), na = "")
  readr::write_csv(tri_res$triangle_df, file.path(out_dir, "hybrid_index_q12_triangle.csv"), na = "")
  readr::write_csv(fb_res$fbranch_df, file.path(out_dir, "introgressive_fbranch_summary.csv"), na = "")
  readr::write_csv(cf_res$cf_df, file.path(out_dir, "quartet_concordance_factors_snps2cf.csv"), na = "")
  readr::write_csv(land_res$summary_df, file.path(out_dir, "landscape_mmrr_feems_summary.csv"), na = "")
  mmrr_sel_p <- file.path(out_dir, "mmrr_model_selection.csv")
  if (!is.null(land_res$mmrr_models_df) && is.data.frame(land_res$mmrr_models_df) && nrow(land_res$mmrr_models_df) > 0L && !file.exists(mmrr_sel_p)) readr::write_csv(land_res$mmrr_models_df, mmrr_sel_p, na = "")
  if (nrow(land_res$feems_edges_df) > 0L) readr::write_csv(land_res$feems_edges_df, file.path(out_dir, "feems_migration_edges.csv"), na = "")

  if (nzchar(tables_dir) && normalizePath(tables_dir, mustWork = FALSE) != normalizePath(out_dir, mustWork = FALSE)) {
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(full_audit, file.path(tables_dir, "genomic_assignments.csv"), na = "")
    readr::write_csv(phylo_res$clade_summary, file.path(tables_dir, "genomic_clade_cohesion_summary.csv"), na = "")
    readr::write_csv(metrics_df, file.path(tables_dir, "genomic_model_metrics.csv"), na = "")
    readr::write_csv(ce_res$table_df, file.path(tables_dir, "admixture_k_model_selection.csv"), na = "")
    readr::write_csv(congr_df, file.path(tables_dir, "genomic_multimodal_congruence_summary.csv"), na = "")
    readr::write_csv(pub_eval$concordance_df, file.path(tables_dir, "public_backbone_concordance.csv"), na = "")
    readr::write_csv(pub_eval$placement_df, file.path(tables_dir, "public_organellar_placement.csv"), na = "")
  }

  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  fig_a <- .resolve_first_nonempty(opts[["figure-a"]], file.path(fig_dir, "fig06a_genomic_pcoa_tree_concordance.png"))
  fig_b <- .resolve_first_nonempty(opts[["figure-b"]], file.path(fig_dir, "fig06b_reticulation_gdi_landscape.png"))
  fig_c <- .resolve_first_nonempty(opts[["figure-c"]], file.path(fig_dir, "fig06c_cytonuclear_paralog_cytotype_audit.png"))
  pass_df$coancestry_cluster <- coanc_res$cluster
  render_fig06a_genomic_pcoa_tree_png(fig_a, pass_df, assign_df, q_mat, eig_pct, phylo_res$tree, congr_df, eff_db_ro, opts[["run-id"]])
  render_fig06b_reticulation_gdi_png(fig_b, tri_res$triangle_df, fb_res$fbranch_df, gdi_df, land_res$feems_edges_df, land_res$summary_df, organellar_placement_df = pub_eval$placement_df)
  render_fig06c_cytonuclear_paralog_cytotype_png(fig_c, parsed_dir, out_dir, full_audit, pub_concordance_df = pub_eval$concordance_df, organellar_placement_df = pub_eval$placement_df, db_path = eff_db_ro)

  message(sprintf("[Step 06/07] Phylogenomics complete: %d/%d passing | LORO acc=%.3f kappa=%.3f | M/P conc=%.4f | median gdi=%.3f",
                  nrow(pass_df), nrow(samples_df), ifelse(is.finite(loro_acc), loro_acc, NA_real_), ifelse(is.finite(loro_kappa), loro_kappa, NA_real_),
                  ifelse(is.finite(mean(rep_conc, na.rm = TRUE)), mean(rep_conc, na.rm = TRUE), NA_real_), ifelse(any(is.finite(gdi_vec)), stats::median(gdi_vec, na.rm = TRUE), NA_real_)))
  .emit_s06_prog("phylogenomics_complete", 4L, 4L, unit = "stages", substage_weight = 0.25, substage_offset = 0.75, status = "completed")
  invisible(assign_df)
}

main <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  run_step06_genomic_phylogenomics(parse_step06_args(cli_args))
}

if (sys.nframe() == 0L) main()
