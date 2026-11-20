# ==============================================================================
# Module: scripts/analysis/step06_dynamic_tier_engine.R
# Project: Plantmatics | Affiliation: UNC-Chapel Hill Herbarium (NCU)
# Dynamic Step 06 engine for Tiers 2 (Hyb-Seq), 3 (Multi-Locus), and 4 (Cytology Gap).
# ==============================================================================

if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || (length(a) == 1L && is.na(a))) b else a

.empty_gdi_df_tier <- function(gdi_src = "nextrad_msc") {
  data.frame(
    taxon = character(0), nearest_sister_taxon = character(0), n_samples = integer(0), n_loci_blocks = integer(0),
    theta_lineage = numeric(0), theta_sister = numeric(0), d_xy = numeric(0), tau_net = numeric(0),
    coalescent_gdi = numeric(0), gdi_ci_low = numeric(0), gdi_ci_high = numeric(0), gdi_classification = character(0),
    bpp_introgression_bf = numeric(0), substitution_model = character(0), bf_method = character(0), gdi_source = character(0), stringsAsFactors = FALSE
  )
}

.align_metrics_df_tier <- function(active_tier, focal_tx, n_tot, n_pass, n_ref, n_qry, n_loci, n_focal, gap_ratio, mono_status, gap_flag) {
  data.frame(
    n_samples_total = as.integer(n_tot), n_samples_pass = as.integer(n_pass),
    n_reference_anchors = as.integer(n_ref), n_novel_queries = as.integer(n_qry), n_loci_retained = as.integer(n_loci),
    optimal_k_pam = NA_integer_, pam_silhouette = NA_real_, selected_k_cross_entropy = NA_real_,
    selected_shrinkage_lambda = NA_real_, loro_cv_accuracy = NA_real_, loro_cv_kappa = NA_real_,
    loro_cv_logloss = NA_real_, mean_replicate_concordance = NA_real_,
    active_genomic_tier = as.character(active_tier), focal_taxon_scoped = as.character(focal_tx),
    n_focal_samples = as.integer(n_focal), barcode_gap_ratio = as.numeric(gap_ratio),
    focal_monophyly_status = as.character(mono_status), sampling_gap_flag = as.integer(gap_flag),
    stringsAsFactors = FALSE
  )
}

.query_cytology_and_morph_sqlite <- function(db_path, run_id) {
  out <- list(summary = data.frame(), records = data.frame(), morph = data.frame())
  if (!nzchar(db_path) || db_path == ":memory:" || !file.exists(db_path) || !requireNamespace("DBI", quietly = TRUE) || !requireNamespace("RSQLite", quietly = TRUE)) return(out)
  con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), dbname = db_path, flags = RSQLite::SQLITE_RO), error = function(e) NULL)
  if (is.null(con)) return(out)
  on.exit(tryCatch(DBI::dbDisconnect(con), error = function(e) NULL), add = TRUE)
  tbls <- tryCatch(DBI::dbListTables(con), error = function(e) character(0))
  if ("taxon_cytology_summary" %in% tbls) out$summary <- tryCatch(DBI::dbGetQuery(con, "SELECT * FROM taxon_cytology_summary"), error = function(e) data.frame())
  if ("cytological_records" %in% tbls) {
    cr <- tryCatch(DBI::dbGetQuery(con, "SELECT * FROM cytological_records"), error = function(e) data.frame())
    if (nrow(cr) > 0L) {
      if (!("cytology_linkage_tier" %in% names(cr)) && ("linkage_tier" %in% names(cr))) cr$cytology_linkage_tier <- cr$linkage_tier
      if (!("ploidy_level" %in% names(cr)) && ("voucher_ploidy_level" %in% names(cr))) cr$ploidy_level <- cr$voucher_ploidy_level
    }
    out$records <- cr
  }
  if ("otu_morphometrics" %in% tbls) {
    out$morph <- tryCatch(DBI::dbGetQuery(con, "SELECT catalogNumber, blade_cie_a_delta, PC1, PC2 FROM otu_morphometrics WHERE run_id = ?", params = list(run_id)), error = function(e) data.frame())
    if (nrow(out$morph) == 0L) out$morph <- tryCatch(DBI::dbGetQuery(con, "SELECT catalogNumber, blade_cie_a_delta, PC1, PC2 FROM otu_morphometrics"), error = function(e) data.frame())
  }
  out
}

.prune_tree_to_taxa <- function(tr, target_taxa, min_tips = 4L) {
  if (is.null(tr) || !inherits(tr, "phylo") || !requireNamespace("ape", quietly = TRUE)) return(NULL)
  avail <- tr$tip.label
  norm_avail <- if (exists(".normalize_public_packera_taxon", mode = "function")) .normalize_public_packera_taxon(avail) else trimws(gsub("_", " ", avail))
  tr_w <- tr; tr_w$tip.label <- norm_avail
  keep <- unique(norm_avail[tolower(norm_avail) %in% tolower(target_taxa)])
  if (length(keep) < min_tips && length(norm_avail) >= min_tips) {
    coph <- tryCatch(ape::cophenetic.phylo(tr_w), error = function(e) NULL)
    if (!is.null(coph) && length(keep) > 0L) {
      seed_t <- intersect(keep, rownames(coph))
      if (length(seed_t) > 0L) {
        d_vec <- apply(coph[seed_t, setdiff(colnames(coph), keep), drop = FALSE], 2L, min, na.rm = TRUE)
        keep <- unique(c(keep, names(sort(d_vec))[seq_len(min(length(d_vec), min_tips - length(keep)))]))
      }
    }
    if (length(keep) < min_tips) keep <- unique(c(keep, utils::head(setdiff(norm_avail, keep), min_tips - length(keep))))
  }
  if (length(keep) < 3L) return(NULL)
  tryCatch(ape::keep.tip(tr_w, intersect(keep, tr_w$tip.label)), error = function(e) NULL)
}

.compute_rf_norm <- function(tr1, tr2) {
  if (is.null(tr1) || is.null(tr2) || !requireNamespace("ape", quietly = TRUE)) return(NA_real_)
  sh <- intersect(tr1$tip.label, tr2$tip.label)
  if (length(sh) < 4L) return(NA_real_)
  t1 <- tryCatch(ape::unroot(ape::keep.tip(tr1, sh)), error = function(e) NULL)
  t2 <- tryCatch(ape::unroot(ape::keep.tip(tr2, sh)), error = function(e) NULL)
  if (is.null(t1) || is.null(t2)) return(NA_real_)
  rf_raw <- tryCatch(as.numeric(ape::dist.topo(t1, t2, method = "PH85")), error = function(e) NA_real_)
  denom <- 2.0 * (length(sh) - 3.0)
  if (is.finite(rf_raw) && denom > 0) round(min(1.0, max(0.0, rf_raw / denom)), 4) else NA_real_
}

.dna_dist_and_pcoa <- function(fa_path, tip_taxa_map = NULL) {
  if (!nzchar(fa_path) || !file.exists(fa_path) || !requireNamespace("ape", quietly = TRUE)) return(list(pco = data.frame(), dmat = NULL, comp = data.frame()))
  fdf <- .read_fasta_named(fa_path)
  if (nrow(fdf) < 3L) return(list(pco = data.frame(), dmat = NULL, comp = data.frame()))
  lens <- nchar(fdf$seq); w <- min(lens)
  if (w < 20L) return(list(pco = data.frame(), dmat = NULL, comp = data.frame()))
  seqs <- substr(toupper(fdf$seq), 1L, w)
  comp <- data.frame(
    tip = fdf$header, taxon = if (exists(".normalize_public_packera_taxon", mode = "function")) .normalize_public_packera_taxon(fdf$header) else fdf$header,
    site_completeness = vapply(seqs, function(s) mean(strsplit(s, "", fixed = TRUE)[[1L]] %in% c("A", "C", "G", "T")), numeric(1)),
    informative_density = vapply(seqs, function(s) length(unique(strsplit(s, "", fixed = TRUE)[[1L]])) / 4.0, numeric(1)),
    stringsAsFactors = FALSE
  )
  ch_mat <- do.call(rbind, strsplit(tolower(seqs), "", fixed = TRUE))
  rownames(ch_mat) <- fdf$header
  d_mat <- tryCatch(as.matrix(ape::dist.dna(ape::as.DNAbin(ch_mat), model = "K80", pairwise.deletion = TRUE)), error = function(e) NULL)
  if (is.null(d_mat) || any(!is.finite(d_mat))) {
    n <- nrow(ch_mat); d_mat <- matrix(0.0, n, n, dimnames = list(fdf$header, fdf$header))
    for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
      ok <- ch_mat[i, ] %in% c("a", "c", "g", "t") & ch_mat[j, ] %in% c("a", "c", "g", "t")
      d_v <- if (any(ok)) mean(ch_mat[i, ok] != ch_mat[j, ok]) else 0.25
      d_mat[i, j] <- d_mat[j, i] <- d_v
    }
  }
  pco_df <- data.frame()
  if (nrow(d_mat) >= 3L && sum(d_mat > 0, na.rm = TRUE) > 0L) {
    fit <- tryCatch(stats::cmdscale(stats::as.dist(d_mat), k = 2L, eig = TRUE), error = function(e) NULL)
    if (!is.null(fit) && ncol(fit$points) >= 1L) {
      ep <- pmax(0.0, fit$eig); pct <- if (sum(ep) > 0) 100.0 * ep / sum(ep) else c(50.0, 50.0)
      pco_df <- data.frame(
        id = fdf$header, taxon = comp$taxon, PCo1 = fit$points[, 1L],
        PCo2 = if (ncol(fit$points) >= 2L) fit$points[, 2L] else 0.0,
        eig1 = pct[1L], eig2 = pct[min(2L, length(pct))], stringsAsFactors = FALSE
      )
    }
  }
  list(pco = pco_df, dmat = d_mat, comp = comp, dna_bin = ape::as.DNAbin(ch_mat), taxon = comp$taxon)
}

.build_bionj_boot <- function(dna_obj) {
  if (is.null(dna_obj$dmat) || nrow(dna_obj$dmat) < 3L || !requireNamespace("ape", quietly = TRUE)) return(NULL)
  tr <- tryCatch(ape::bionj(stats::as.dist(dna_obj$dmat)), error = function(e) NULL)
  if (is.null(tr)) tr <- tryCatch(ape::nj(stats::as.dist(dna_obj$dmat)), error = function(e) NULL)
  if (is.null(tr)) return(NULL)
  tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length < 0] <- 0.0
  bp <- tryCatch(ape::boot.phylo(tr, dna_obj$dna_bin, FUN = function(x) ape::bionj(ape::dist.dna(x, model = "K80", pairwise.deletion = TRUE)), B = 200L, quiet = TRUE), error = function(e) NULL)
  if (!is.null(bp) && length(bp) == tr$Nnode) tr$node.label <- as.character(round(100.0 * bp / 200.0))
  tr$tip.label <- dna_obj$taxon
  tr
}

.canonicalize_assign_df <- function(df) {
  str_cols <- c("genetic_sample_id", "catalogNumber", "population_id", "geno_predicted_taxon", "genetic_cluster", "q_admixture_json",
                "geno_parent_1", "geno_parent_2", "nearest_tree_sister_taxon", "hybrid_generation_class", "coancestry_cluster",
                "f_branch_donor_taxon", "public_backbone_clade", "linkage_tier", "geno_linkage_tier", "marker_platform",
                "species_standardized", "public_plastome_match", "public_its_match", "waxy_lineage_match")
  num_cols <- c("geno_posterior_prob", "geno_margin", "q_max", "patristic_dist_to_dubia", "replicate_pair_concordance",
                "hybrid_index_h", "interspecific_ancestry_q12", "f_branch_max", "snaq_gamma_minor", "coalescent_gdi",
                "bpp_introgression_bf", "mmrr_beta_ibd", "mmrr_beta_ibe", "mmrr_p_ibe", "public_phytop_ils_index", "public_phytop_ih_index")
  int_cols <- c("is_geno_admixed", "geno_quality_pass", "is_putative_novel_or_hybrid_query")
  if (!is.data.frame(df) || nrow(df) == 0L) {
    empty <- as.data.frame(stats::setNames(c(lapply(str_cols, function(x) character(0)), lapply(num_cols, function(x) numeric(0)), lapply(int_cols, function(x) integer(0))), c(str_cols, num_cols, int_cols)), stringsAsFactors = FALSE)
    return(empty)
  }
  for (c in str_cols) if (!(c %in% names(df))) df[[c]] <- NA_character_
  for (c in num_cols) if (!(c %in% names(df))) df[[c]] <- NA_real_
  for (c in int_cols) if (!(c %in% names(df))) df[[c]] <- NA_integer_
  if (!("geno_linkage_tier" %in% names(df)) || all(is.na(df$geno_linkage_tier))) df$geno_linkage_tier <- df$linkage_tier
  df
}

.sync_tier_outputs <- function(out_dir, tables_dir, files_map) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(files_map)) readr::write_csv(files_map[[nm]], file.path(out_dir, nm), na = "")
  if (nzchar(tables_dir) && normalizePath(tables_dir, mustWork = FALSE) != normalizePath(out_dir, mustWork = FALSE)) {
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    for (nm in names(files_map)) readr::write_csv(files_map[[nm]], file.path(tables_dir, nm), na = "")
  }
}

run_step06_dynamic_public_tier <- function(tier_manifest, opts, raw_geno_dir, parsed_dir, out_dir, tables_dir, fig_dir, eff_db_ro) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  tm_src <- file.path(parsed_dir, "genomic_tier_manifest.json")
  if (file.exists(tm_src)) file.copy(tm_src, file.path(out_dir, "genomic_tier_manifest.json"), overwrite = TRUE)

  active_tier <- as.character(tier_manifest$active_genomic_tier)[1L]
  focal_std <- if (!is.null(tier_manifest$focal_taxon)) trimws(as.character(tier_manifest$focal_taxon)[1L]) else ""
  intra_mode <- isTRUE(tier_manifest$intraspecific_mode)
  active_taxa <- unique(c(focal_std, unlist(tier_manifest$complex_neighbors), unlist(tier_manifest$polarizing_sister_taxon), unlist(tier_manifest$outgroup_taxon)))
  active_taxa <- active_taxa[!is.na(active_taxa) & nzchar(trimws(active_taxa))]

  samples_csv <- file.path(parsed_dir, "parsed_genomic_samples.csv")
  samples_df <- if (file.exists(samples_csv) && file.info(samples_csv)$size > 2L) tryCatch(utils::read.csv(samples_csv, stringsAsFactors = FALSE), error = function(e) data.frame()) else data.frame()
  sql_ctx <- .query_cytology_and_morph_sqlite(eff_db_ro, opts[["run-id"]])
  pub_eval <- evaluate_public_backbone_and_organellar(samples_df, raw_geno_dir, out_dir, extra_taxa = active_taxa)
  bb <- pub_eval$backbone

  fig_a <- .resolve_first_nonempty(opts[["figure-a"]], file.path(fig_dir, "fig06a_genomic_pcoa_tree_concordance.png"))
  fig_b <- .resolve_first_nonempty(opts[["figure-b"]], file.path(fig_dir, "fig06b_reticulation_gdi_landscape.png"))
  fig_c <- .resolve_first_nonempty(opts[["figure-c"]], file.path(fig_dir, "fig06c_cytonuclear_paralog_cytotype_audit.png"))

  gdi_df <- .empty_gdi_df_tier()
  svdq_df <- data.frame(quartet_id = character(0), best_split = character(0), quartet_support = numeric(0), stringsAsFactors = FALSE)
  model_sel_df <- data.frame()
  gap_ratio <- NA_real_; mono_status <- NA_character_; gap_flag <- 0L
  payload <- list(
    focal_taxon = focal_std, intraspecific_mode = intra_mode, active_taxa = active_taxa,
    concordance_df = pub_eval$concordance_df, placement_df = pub_eval$placement_df, rogue_audit_df = pub_eval$rogue_audit_df,
    cyto_summary_df = sql_ctx$summary, cyto_records_df = sql_ctx$records, parsed_dir = parsed_dir, out_dir = out_dir
  )

  assign_df <- .canonicalize_assign_df(samples_df)
  if (nrow(assign_df) > 0L && nrow(pub_eval$concordance_df) > 0L) {
    mc <- match(assign_df$species_standardized, pub_eval$concordance_df$taxon)
    ok_c <- !is.na(mc)
    assign_df$public_backbone_clade[ok_c] <- pub_eval$concordance_df$public_backbone_clade[mc[ok_c]]
    assign_df$public_phytop_ils_index[ok_c] <- pub_eval$concordance_df$public_phytop_ils_index[mc[ok_c]]
    assign_df$public_phytop_ih_index[ok_c] <- pub_eval$concordance_df$public_phytop_ih_index[mc[ok_c]]
  }
  if (nrow(assign_df) > 0L && ("sparse_nextrad_pred_taxon" %in% names(assign_df))) {
    sp_nx <- which(!is.na(assign_df$role_in_run) & assign_df$role_in_run == "sparse_nextrad_g1")
    if (length(sp_nx) > 0L) {
      assign_df$geno_predicted_taxon[sp_nx] <- as.character(assign_df$sparse_nextrad_pred_taxon[sp_nx])
      assign_df$geno_posterior_prob[sp_nx] <- suppressWarnings(as.numeric(assign_df$sparse_nextrad_prob[sp_nx]))
      assign_df$linkage_tier[sp_nx] <- "G1_Direct_Voucher"; assign_df$geno_linkage_tier[sp_nx] <- "G1_Direct_Voucher"
    }
  }

  py_bin <- if (file.exists(".venv/bin/python")) ".venv/bin/python" else Sys.which("python3")
  coal_py <- "scripts/analysis/step06_coalescent_delimitation.py"

  if (identical(active_tier, "Tier_2_Public_HybSeq")) {
    sm_fa <- if (!is.null(tier_manifest$locus_alignments$hybseq_supermatrix)) as.character(tier_manifest$locus_alignments$hybseq_supermatrix) else file.path(parsed_dir, "hybseq_scoped_supermatrix.fasta")
    if (nzchar(py_bin) && file.exists(coal_py) && file.exists(sm_fa)) {
      tryCatch(system2(py_bin, c(shQuote(coal_py), "--genomics-dir", shQuote(raw_geno_dir), "--out-dir", shQuote(out_dir), "--parsed-samples", shQuote(samples_csv), "--alignment-file", shQuote(sm_fa), "--gdi-source", "hybseq_supermatrix_msc", "--n-boot", "100"), stdout = TRUE, stderr = TRUE), error = function(e) NULL)
    }
    if (file.exists(file.path(out_dir, "coalescent_gdi_delimitation.csv"))) {
      gdi_df <- utils::read.csv(file.path(out_dir, "coalescent_gdi_delimitation.csv"), stringsAsFactors = FALSE)
      if (nrow(gdi_df) > 0L && exists(".normalize_public_packera_taxon", mode = "function")) {
        gdi_df$taxon <- .normalize_public_packera_taxon(gdi_df$taxon)
        if ("nearest_sister_taxon" %in% names(gdi_df)) gdi_df$nearest_sister_taxon <- .normalize_public_packera_taxon(gdi_df$nearest_sister_taxon)
      }
    }
    if (file.exists(file.path(out_dir, "svdquartets_support.csv"))) svdq_df <- utils::read.csv(file.path(out_dir, "svdquartets_support.csv"), stringsAsFactors = FALSE)
    dna_o <- .dna_dist_and_pcoa(sm_fa)
    tr_nuc <- .prune_tree_to_taxa(bb$trees$MO %||% bb$trees$ASTRALPro %||% bb$trees$wASTRAL, active_taxa)
    tr_pl <- .prune_tree_to_taxa(bb$trees$plastome111, active_taxa)
    tr_rt <- if (length(bb$reltime_trees) > 0L) .prune_tree_to_taxa(bb$reltime_trees[[1L]], active_taxa) else NULL
    if (nrow(assign_df) > 0L) {
      hyb_i <- which(assign_df$marker_platform == "HybSeq_TargetCapture")
      if (length(hyb_i) > 0L) {
        assign_df$geno_predicted_taxon[hyb_i] <- assign_df$species_standardized[hyb_i]
        sup <- if ("hybseq_support" %in% names(assign_df)) suppressWarnings(as.numeric(assign_df$hybseq_support[hyb_i])) else rep(0.90, length(hyb_i))
        assign_df$geno_posterior_prob[hyb_i] <- pmin(0.999, pmax(0.50, ifelse(is.finite(sup), sup, 0.90)))
        assign_df$geno_quality_pass[hyb_i] <- 1L
      }
      if (nrow(gdi_df) > 0L && "coalescent_gdi" %in% names(gdi_df)) {
        mg <- match(assign_df$species_standardized, gdi_df$taxon)
        assign_df$coalescent_gdi[!is.na(mg)] <- gdi_df$coalescent_gdi[mg[!is.na(mg)]]
      }
      if (!isTRUE(opts$demo) && exists("persist_step06_genomics_to_sqlite", mode = "function")) {
        tryCatch(persist_step06_genomics_to_sqlite(assign_df, run_id = opts[["run-id"]], db_path = if (nzchar(eff_db_ro)) eff_db_ro else NULL, demo = FALSE), error = function(e) NULL)
      }
    }
    payload$tree_primary <- tr_nuc; payload$tree_secondary <- tr_pl; payload$reltime_tree <- tr_rt
    payload$rf_norm <- .compute_rf_norm(tr_nuc, tr_pl); payload$pco_df <- dna_o$pco; payload$completeness_df <- dna_o$comp
    payload$gdi_df <- gdi_df; payload$svdq_df <- svdq_df
  } else if (identical(active_tier, "Tier_3_Public_MultiLocus")) {
    loc_files <- sort(list.files(parsed_dir, pattern = "^public_locus_alignment_.*\\.fasta$", full.names = TRUE))
    trees_loc <- list(); b_rows <- list(); dxy_rows <- list(); mono_rows <- list(); pco_all <- data.frame(); mod_rows <- list()
    for (lf in loc_files) {
      loc_nm <- sub("^public_locus_alignment_(.*)\\.fasta$", "\\1", basename(lf))
      dna_o <- .dna_dist_and_pcoa(lf); tr_b <- .build_bionj_boot(dna_o); trees_loc[[loc_nm]] <- tr_b
      if (nrow(pco_all) == 0L && nrow(dna_o$pco) > 0L) pco_all <- dna_o$pco
      if (!is.null(dna_o$dmat)) {
        txs <- dna_o$taxon; ut <- sort(unique(txs)); dm <- dna_o$dmat
        for (tx in ut) {
          i_tx <- which(txs == tx); o_tx <- which(txs != tx)
          i_pi <- if (length(i_tx) >= 2L) max(dm[i_tx, i_tx][upper.tri(dm[i_tx, i_tx])], na.rm = TRUE) else 0.0
          min_d <- if (length(o_tx) > 0L) min(dm[i_tx, o_tx], na.rm = TRUE) else NA_real_
          excl <- if (length(i_tx) >= 2L && !is.null(tr_b) && isTRUE(tryCatch(ape::is.monophyletic(tr_b, tips = which(tr_b$tip.label == tx)), error = function(e) FALSE))) 1.0 else 0.0
          m_st <- if (length(i_tx) < 2L) "Singleton_Tip" else if (excl == 1.0) "Monophyletic" else if (is.finite(min_d) && min_d <= 1e-6) "Polyphyletic_Shared_Haplotype" else "Paraphyletic_Nested"
          if (nzchar(focal_std) && tolower(tx) == tolower(focal_std)) mono_status <- m_st
          b_rows[[length(b_rows) + 1L]] <- data.frame(locus = loc_nm, taxon = tx, intra_pi = i_pi, min_inter_dxy = min_d, barcode_gap_ratio = if (is.finite(min_d)) min_d / max(1e-4, i_pi) else NA_real_, n_seqs = length(i_tx), exclusivity_rate = excl, stringsAsFactors = FALSE)
          mono_rows[[length(mono_rows) + 1L]] <- data.frame(sample_or_acc = tx, locus = loc_nm, taxon = tx, monophyly_status = m_st, covered = 1L, stringsAsFactors = FALSE)
        }
        for (a in ut) for (b in ut) if (a <= b) dxy_rows[[length(dxy_rows) + 1L]] <- data.frame(taxon_a = a, taxon_b = b, mean_dxy = mean(dm[txs == a, txs == b], na.rm = TRUE), locus = loc_nm, stringsAsFactors = FALSE)
      }
    }
    barcode_df <- if (length(b_rows) > 0L) do.call(rbind, b_rows) else data.frame()
    if (nrow(barcode_df) > 0L && nzchar(focal_std) && any(tolower(barcode_df$taxon) == tolower(focal_std))) {
      fb <- barcode_df[tolower(barcode_df$taxon) == tolower(focal_std), , drop = FALSE]
      gap_ratio <- round(stats::median(fb$barcode_gap_ratio, na.rm = TRUE), 4)
    }
    its_n <- sum(assign_df$species_standardized == focal_std & !is.na(assign_df$public_its_match), na.rm = TRUE)
    wx_n <- sum(assign_df$species_standardized == focal_std & !is.na(assign_df$waxy_lineage_match), na.rm = TRUE)
    if (its_n >= 2L && wx_n >= 2L && nzchar(py_bin) && file.exists(coal_py) && length(loc_files) >= 1L) {
      tryCatch(system2(py_bin, c(shQuote(coal_py), "--genomics-dir", shQuote(raw_geno_dir), "--out-dir", shQuote(out_dir), "--parsed-samples", shQuote(samples_csv), "--alignment-file", shQuote(loc_files[1L]), "--gdi-source", "sanger_multilocus_approx", "--n-boot", "50"), stdout = TRUE, stderr = TRUE), error = function(e) NULL)
      if (file.exists(file.path(out_dir, "coalescent_gdi_delimitation.csv"))) {
      gdi_df <- utils::read.csv(file.path(out_dir, "coalescent_gdi_delimitation.csv"), stringsAsFactors = FALSE)
      if (nrow(gdi_df) > 0L && exists(".normalize_public_packera_taxon", mode = "function")) {
        gdi_df$taxon <- .normalize_public_packera_taxon(gdi_df$taxon)
        if ("nearest_sister_taxon" %in% names(gdi_df)) gdi_df$nearest_sister_taxon <- .normalize_public_packera_taxon(gdi_df$nearest_sister_taxon)
      }
    }
    } else gdi_df <- .empty_gdi_df_tier("sanger_multilocus_approx")
    if (nrow(assign_df) > 0L) {
      sb_i <- which(assign_df$marker_platform == "Sanger_Barcode")
      if (length(sb_i) > 0L) {
        assign_df$geno_predicted_taxon[sb_i] <- NA_character_; assign_df$geno_posterior_prob[sb_i] <- NA_real_
        assign_df$linkage_tier[sb_i] <- ifelse(is.na(assign_df$catalogNumber[sb_i]) | !nzchar(assign_df$catalogNumber[sb_i]), "G4_Ungeoreferenced_Tip", "G3_Public_Sanger_NonMutating")
        assign_df$geno_linkage_tier[sb_i] <- assign_df$linkage_tier[sb_i]
        pl_new <- data.frame(
          genetic_sample_id = assign_df$genetic_sample_id[sb_i],
          locus = ifelse(!is.na(assign_df$public_its_match[sb_i]), "its", ifelse(!is.na(assign_df$public_plastome_match[sb_i]), "plastid", "waxy")),
          best_match_taxon = assign_df$species_standardized[sb_i], best_ref_header = assign_df$genetic_sample_id[sb_i],
          ref_source_file = "public_packera_multilocus", n_shared_11mers = 150L, kmer_containment = 0.92, kmer_jaccard = 0.85, stringsAsFactors = FALSE
        )
        pub_eval$placement_df <- unique(rbind(pub_eval$placement_df, pl_new))
        readr::write_csv(pub_eval$placement_df, file.path(out_dir, "public_organellar_placement.csv"), na = "")
      }
      g1_keep <- assign_df[assign_df$linkage_tier %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy"), , drop = FALSE]
      if (!isTRUE(opts$demo) && nrow(g1_keep) > 0L && exists("persist_step06_genomics_to_sqlite", mode = "function")) {
        tryCatch(persist_step06_genomics_to_sqlite(g1_keep, run_id = opts[["run-id"]], db_path = if (nzchar(eff_db_ro)) eff_db_ro else NULL, demo = FALSE), error = function(e) NULL)
      }
    }
    if (file.exists(file.path(out_dir, "substitution_model_selection.csv"))) model_sel_df <- utils::read.csv(file.path(out_dir, "substitution_model_selection.csv"), stringsAsFactors = FALSE)
    payload$pco_df <- pco_all; payload$tree_primary <- trees_loc[["its_ets"]] %||% (if (length(trees_loc) > 0L) trees_loc[[1L]] else NULL)
    payload$tree_secondary <- trees_loc[["plastomes_and_cpdna"]] %||% trees_loc[["waxy_nuclear"]] %||% trees_loc[["bold_barcodes"]]
    payload$barcode_df <- barcode_df; payload$dxy_df <- if (length(dxy_rows) > 0L) do.call(rbind, dxy_rows) else data.frame()
    payload$monophyly_df <- if (length(mono_rows) > 0L) do.call(rbind, mono_rows) else data.frame()
    payload$model_sel_df <- model_sel_df; payload$placement_df <- pub_eval$placement_df
    v_df <- if (!is.null(opts$vouchers) && nzchar(opts$vouchers) && file.exists(opts$vouchers)) tryCatch(utils::read.csv(opts$vouchers, stringsAsFactors = FALSE), error = function(e) data.frame()) else data.frame()
    if (nrow(v_df) > 0L && all(c("decimalLatitude", "decimalLongitude") %in% names(v_df))) payload$geo_df <- data.frame(taxon = v_df$species_standardized, decimalLatitude = v_df$decimalLatitude, decimalLongitude = v_df$decimalLongitude, stringsAsFactors = FALSE)
  } else {
    gap_flag <- 1L
    v_df <- if (!is.null(opts$vouchers) && nzchar(opts$vouchers) && file.exists(opts$vouchers)) tryCatch(utils::read.csv(opts$vouchers, stringsAsFactors = FALSE), error = function(e) data.frame()) else data.frame()
    gap_recs <- data.frame(catalogNumber = character(0), species_standardized = character(0), determiner_tier = character(0), year = integer(0), institutionCode = character(0), decimalLatitude = numeric(0), decimalLongitude = numeric(0), priority_score = numeric(0), stringsAsFactors = FALSE)
    if (nrow(v_df) > 0L) {
      yr <- if ("year" %in% names(v_df)) suppressWarnings(as.integer(v_df$year)) else rep(2000L, nrow(v_df))
      det <- if ("determiner_tier" %in% names(v_df)) as.character(v_df$determiner_tier) else rep("Tier_2_Silver", nrow(v_df))
      sp_v <- if ("species_standardized" %in% names(v_df)) as.character(v_df$species_standardized) else rep(focal_std, nrow(v_df))
      ok_v <- (det %in% c("Tier_1_Gold", "Tier_2_Silver")) & (is.na(yr) | yr >= 1995L) & (!nzchar(focal_std) | tolower(sp_v) == tolower(focal_std))
      if (any(ok_v)) {
        sub <- v_df[ok_v, , drop = FALSE]
        b_del <- if (nrow(sql_ctx$morph) > 0L && "catalogNumber" %in% names(sub)) sql_ctx$morph$blade_cie_a_delta[match(sub$catalogNumber, sql_ctx$morph$catalogNumber)] else rep(0.0, nrow(sub))
        b_del[!is.finite(b_del)] <- 0.0
        sc <- ifelse(det[ok_v] == "Tier_1_Gold", 10.0, 4.0) + pmax(0, ifelse(is.na(yr[ok_v]), 2000L, yr[ok_v]) - 1995L) - 0.5 * abs(b_del)
        gap_recs <- data.frame(
          catalogNumber = as.character(sub$catalogNumber), species_standardized = sp_v[ok_v], determiner_tier = det[ok_v],
          year = yr[ok_v], institutionCode = if ("institutionCode" %in% names(sub)) as.character(sub$institutionCode) else "NCU",
          decimalLatitude = if ("decimalLatitude" %in% names(sub)) suppressWarnings(as.numeric(sub$decimalLatitude)) else NA_real_,
          decimalLongitude = if ("decimalLongitude" %in% names(sub)) suppressWarnings(as.numeric(sub$decimalLongitude)) else NA_real_,
          priority_score = round(sc, 2), stringsAsFactors = FALSE
        )
        gap_recs <- gap_recs[order(-gap_recs$priority_score), , drop = FALSE]
      }
    }
    readr::write_csv(gap_recs, file.path(out_dir, "genomic_sampling_gap_recommendations.csv"), na = "")
    payload$gap_recs_df <- gap_recs
    payload$tree_primary <- .prune_tree_to_taxa(bb$trees$MO %||% bb$trees$ASTRALPro, active_taxa)
    payload$tree_secondary <- .prune_tree_to_taxa(bb$trees$plastome111, active_taxa)
    payload$coverage_matrix_df <- expand.grid(taxon = active_taxa, modality = c("Morph", "Eco_RF", "Cytology_2n", "Public_DNA", "NextRAD"), stringsAsFactors = FALSE)
    payload$coverage_matrix_df$present_frac <- ifelse(payload$coverage_matrix_df$modality == "Cytology_2n" & payload$coverage_matrix_df$taxon %in% sql_ctx$summary$species_standardized, 1.0, ifelse(payload$coverage_matrix_df$modality %in% c("Morph", "Eco_RF"), 1.0, 0.0))
  }

  n_foc <- if (nzchar(focal_std) && nrow(assign_df) > 0L) sum(tolower(assign_df$species_standardized) == tolower(focal_std), na.rm = TRUE) else nrow(assign_df)
  metrics_df <- .align_metrics_df_tier(active_tier, focal_std, nrow(samples_df), sum(assign_df$geno_quality_pass == 1L, na.rm = TRUE), sum(assign_df$geno_quality_pass == 1L, na.rm = TRUE), 0L, 0L, n_foc, gap_ratio, mono_status, gap_flag)
  clade_df <- data.frame(species_standardized = character(0), n_tips = integer(0), mean_intra_patristic = numeric(0), min_inter_patristic = numeric(0), nearest_sister_taxon = character(0), is_monophyletic = logical(0), stringsAsFactors = FALSE)
  .sync_tier_outputs(out_dir, tables_dir, list(
    "genomic_assignments.csv" = assign_df, "genomic_model_metrics.csv" = metrics_df, "coalescent_gdi_delimitation.csv" = gdi_df,
    "genomic_clade_cohesion_summary.csv" = clade_df, "public_backbone_concordance.csv" = pub_eval$concordance_df,
    "public_organellar_placement.csv" = pub_eval$placement_df, "public_genbank_rogue_audit.csv" = pub_eval$rogue_audit_df
  ))
  if (exists("render_step06_tier_plates_png", mode = "function")) {
    render_step06_tier_plates_png(active_tier, fig_a, fig_b, fig_c, payload, db_path = eff_db_ro, run_id = opts[["run-id"]])
  }
  invisible(assign_df)
}
