#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/analysis/07_discover_taxon_complexes.R
# Project: Packera Species Delimitation & Morphometrics Pipeline (NCU)
# Description: Unsupervised Quantitative Taxon-Complex Discovery & Integrative
#   Delimitation Engine (Phase 07): Cohen's kappa weights (Dawid & Skene 1979),
#   Empirical Bayes circular Hopkins phenology, 5-Stage Decision Matrix (Padial
#   et al. 2010; Schlick-Steiner et al. 2010), and Leiden/Louvain + Ward.D2.
# ==============================================================================
suppressPackageStartupMessages({
  for (.pkg in c("optparse", "igraph", "cluster", "MASS")) if (requireNamespace(.pkg, quietly = TRUE)) library(.pkg, character.only = TRUE)
})
for (.pf_p in c("scripts/morphometrics/pub_figures_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "morphometrics", "pub_figures_utils.R"))) {
  if (!exists("standardize_taxon_name", mode = "function") && !is.na(.pf_p) && nzchar(.pf_p) && file.exists(.pf_p)) { source(.pf_p); break }
}
for (.aw_p in c("scripts/morphometrics/authority_weighting_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "morphometrics", "authority_weighting_utils.R"))) {
  if (!exists(".identify_type_specimen_vec", mode = "function") && !is.na(.aw_p) && nzchar(.aw_p) && file.exists(.aw_p)) { source(.aw_p); break }
}
for (.db_p in c("scripts/core/db_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "..", "core", "db_utils.R"))) {
  if (!is.na(.db_p) && nzchar(.db_p) && file.exists(.db_p)) { source(.db_p); break }
}
# 1. CLI Argument Parser & Atomic Export Helpers
parse_discovery_args <- function() {
  run_dir <- Sys.getenv("PACKERA_RUN_DIR", unset = "")
  rp <- function(rel) if (nzchar(run_dir)) file.path(run_dir, rel) else rel
  defaults <- list(
    morphometrics = rp("data/tables/morphometric_flags.csv"), harmonics = rp("data/tables/leaf_efa_harmonics.csv"),
    vouchers = rp("data/tables/curated_vouchers.csv"),
    multimodal = rp("data/tables/multimodal_conflict_flags.csv"), focal_taxon = "", min_affinity = 0.03,
    output_complexes = rp("outputs/reports/empirical_taxon_complexes.csv"),
    output_affinity = rp("outputs/reports/pairwise_taxon_affinity_matrix.csv"),
    output_neighbors = rp("outputs/reports/focal_auto_neighbors.txt"),
    output_figure = rp("outputs/figures/fig_empirical_taxon_complex_network.pdf")
  )
  raw_args <- commandArgs(trailingOnly = TRUE)
  if (length(raw_args) == 0) return(defaults)
  opts <- defaults; i <- 1L
  flag_map <- c(
    "-m" = "morphometrics", "--morphometrics" = "morphometrics", "-e" = "harmonics", "--harmonics" = "harmonics",
    "-v" = "vouchers", "--vouchers" = "vouchers",
    "-c" = "multimodal", "--multimodal" = "multimodal", "--focal-taxon" = "focal_taxon",
    "--output-complexes" = "output_complexes", "--output-affinity" = "output_affinity",
    "--output-neighbors" = "output_neighbors", "--output-figure" = "output_figure"
  )
  while (i <= length(raw_args)) {
    a <- raw_args[i]
    if (a %in% names(flag_map) && i < length(raw_args)) { opts[[flag_map[[a]]]] <- raw_args[i + 1L]; i <- i + 2L }
    else if (a == "--min-affinity" && i < length(raw_args)) { opts$min_affinity <- as.numeric(raw_args[i + 1L]); i <- i + 2L }
    else i <- i + 1L
  }
  opts
}
atomic_write_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp"); utils::write.csv(df, file = tmp, row.names = FALSE)
  file.rename(tmp, path); invisible(path)
}
atomic_write_lines <- function(lines_vec, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp"); writeLines(as.character(lines_vec), con = tmp)
  file.rename(tmp, path); invisible(path)
}
std_taxon <- function(x) {
  if (exists("standardize_taxon_name", mode = "function")) { res <- standardize_taxon_name(x); return(ifelse(is.na(res), as.character(x), res)) }
  as.character(x)
}
# 2. Step A: Load & Filter Basal_Rosette Gold/Silver Anchor Specimens
load_and_merge_inputs <- function(opts) {
  res <- load_step07_merged_inputs(opts, std_taxon)
  if (!is.null(opts$multimodal) && nzchar(opts$multimodal) && file.exists(opts$multimodal) && !is.null(res$morph_df) && nrow(res$morph_df) > 0L) {
    mm <- tryCatch(utils::read.csv(opts$multimodal, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (!is.null(mm) && nrow(mm) > 0L && "catalogNumber" %in% names(mm)) {
      mm$catalogNumber <- as.character(mm$catalogNumber)
      mm_pid <- if ("plant_individual_id" %in% names(mm)) ifelse(is.na(mm$plant_individual_id), 0L, suppressWarnings(as.integer(mm$plant_individual_id))) else rep(0L, nrow(mm))
      mm_key <- paste0(mm$catalogNumber, "__", ifelse(is.na(mm_pid), 0L, mm_pid))
      m_df <- res$morph_df
      m_pid <- if ("plant_individual_id" %in% names(m_df)) ifelse(is.na(m_df$plant_individual_id), 0L, suppressWarnings(as.integer(m_df$plant_individual_id))) else rep(0L, nrow(m_df))
      m_key <- if ("key" %in% names(m_df)) as.character(m_df$key) else paste0(m_df$catalogNumber, "__", ifelse(is.na(m_pid), 0L, m_pid))
      idx_m <- match(m_key, mm_key); s_cats <- names(which(table(mm$catalogNumber) == 1L))
      idx_fb <- ifelse(m_df$catalogNumber %in% s_cats, match(m_df$catalogNumber, mm$catalogNumber), NA_integer_)
      idx_m <- ifelse(!is.na(idx_m), idx_m, idx_fb)
      extra_cols <- grep("^(primary_(env|clim|edaph)_predictors|edaphic_z_score|climatic_z_score|type_edaphic_z_score|nearest_type_niche_taxon|matches_type_lithology|edaphic_schoener_d|schoener_d|edaphic_niche_overlap_d|bio[0-9]+_|fine_|ilr_texture_|ssurgo_|polaris_|solus_|topo_|soil_)", names(mm), value = TRUE)
      for (ec in extra_cols) {
        if (!(ec %in% names(m_df))) m_df[[ec]] <- mm[[ec]][idx_m]
        else { na_e <- is.na(m_df[[ec]]) & !is.na(idx_m); if (any(na_e)) m_df[[ec]][na_e] <- mm[[ec]][idx_m[na_e]] }
      }
      res$morph_df <- m_df
    }
  }
  # Hydrate canonical Step 06 curation_triage subtypes when Step 06 outputs exist (eliminates Step 06 vs Step 07 drift)
  if (!is.null(res$morph_df) && nrow(res$morph_df) > 0L && !("taxonomic_treatment_subtype" %in% names(res$morph_df))) {
    run_d <- Sys.getenv("PACKERA_RUN_DIR", unset = "")
    tri_cands <- unique(c(
      if (nzchar(run_d)) file.path(run_d, "outputs/reports/final_curation_triage.csv") else character(0),
      if (nzchar(run_d)) file.path(run_d, "data/tables/triage_queue.csv") else character(0),
      file.path(dirname(opts$morphometrics), "triage_queue.csv"),
      "outputs/reports/final_curation_triage.csv",
      "data/tables/triage_queue.csv"
    ))
    tri_hit <- tri_cands[file.exists(tri_cands) & file.info(tri_cands)$size > 0][1L]
    if (!is.na(tri_hit) && nzchar(tri_hit)) {
      tdf <- tryCatch(utils::read.csv(tri_hit, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (!is.null(tdf) && "catalogNumber" %in% names(tdf) && "taxonomic_treatment_subtype" %in% names(tdf)) {
        t_pid <- if ("plant_individual_id" %in% names(tdf)) ifelse(is.na(tdf$plant_individual_id), 0L, suppressWarnings(as.integer(tdf$plant_individual_id))) else rep(0L, nrow(tdf))
        m_pid <- if ("plant_individual_id" %in% names(res$morph_df)) ifelse(is.na(res$morph_df$plant_individual_id), 0L, suppressWarnings(as.integer(res$morph_df$plant_individual_id))) else rep(0L, nrow(res$morph_df))
        idx_t <- match(paste0(res$morph_df$catalogNumber, "__", m_pid), paste0(tdf$catalogNumber, "__", t_pid))
        s_tri <- names(which(table(tdf$catalogNumber) == 1L))
        idx_t_fb <- ifelse(res$morph_df$catalogNumber %in% s_tri, match(res$morph_df$catalogNumber, tdf$catalogNumber), NA_integer_)
        idx_t <- ifelse(!is.na(idx_t), idx_t, idx_t_fb)
        res$morph_df$taxonomic_treatment_subtype <- as.character(tdf$taxonomic_treatment_subtype[idx_t])
      }
    }
  }
  res
}
select_taxon_anchors <- function(morph_df) {
  is_t1_pure <- rep(TRUE, nrow(morph_df))
  if ("morph_tier_source" %in% names(morph_df)) {
    mts <- as.character(morph_df$morph_tier_source); if (any(!is.na(mts) & nzchar(mts))) is_t1_pure <- !is.na(mts) & (mts == "Basal_Rosette")
  }
  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(morph_df) else list(type_rank = rep(0L, nrow(morph_df)), is_type_specimen = rep(FALSE, nrow(morph_df)))
  is_type <- as.logical(type_info$is_type_specimen); is_type[is.na(is_type)] <- FALSE
  det_tier <- if ("determiner_tier" %in% names(morph_df)) ifelse(is.na(morph_df$determiner_tier) | !nzchar(trimws(as.character(morph_df$determiner_tier))), "Tier_3_Bronze", as.character(morph_df$determiner_tier)) else rep("Tier_1_Gold", nrow(morph_df))
  is_gs <- (det_tier %in% c("Tier_1_Gold", "Tier_2_Silver")) | is_type; taxa <- sort(unique(morph_df$species_standardized))
  anchor_mask <- rep(FALSE, nrow(morph_df)); n_gs_vec <- stats::setNames(integer(length(taxa)), taxa); n_type_vec <- stats::setNames(integer(length(taxa)), taxa)
  for (tx in taxa) {
    tx_idx <- which(morph_df$species_standardized == tx); type_idx <- tx_idx[is_type[tx_idx]]
    n_type_vec[tx] <- length(type_idx); gs_pure_idx <- tx_idx[is_t1_pure[tx_idx] & is_gs[tx_idx]]; n_gs_vec[tx] <- length(gs_pure_idx)
    if (length(gs_pure_idx) >= 5L) anchor_mask[gs_pure_idx] <- TRUE
    else if (length(tx_idx[is_t1_pure[tx_idx]]) >= 3L) anchor_mask[tx_idx[is_t1_pure[tx_idx]]] <- TRUE
    else anchor_mask[tx_idx] <- TRUE
    if (length(type_idx) > 0L) anchor_mask[type_idx] <- TRUE
  }
  list(anchor_df = morph_df[anchor_mask, , drop = FALSE], anchor_mask = anchor_mask, n_gold_silver = n_gs_vec, n_type_specimens = n_type_vec)
}
# 3. Step B: Empirical Cohen's Kappa Reliability Weights (Dawid & Skene 1979)
compute_cohens_kappa_stat <- function(y_true, y_pred) {
  valid <- !is.na(y_true) & nzchar(y_true) & !is.na(y_pred) & nzchar(y_pred) & !(y_pred %in% c("Unknown", "Unresolved", "NA"))
  if (sum(valid) < 2L) return(NA_real_)
  yt <- as.character(y_true[valid]); yp <- as.character(y_pred[valid]); lvls <- sort(unique(c(yt, yp)))
  if (length(lvls) < 2L || length(unique(yt)) < 2L) return(max(0.05, min(0.95, mean(yt == yp))))
  tab <- table(factor(yt, levels = lvls), factor(yp, levels = lvls)); n_tot <- sum(tab)
  p_o <- sum(diag(tab)) / n_tot; p_e <- sum((rowSums(tab) / n_tot) * (colSums(tab) / n_tot))
  if (p_e >= 1.0 - 1e-12) return(1.0)
  (p_o - p_e) / (1.0 - p_e)
}
compute_modality_kappa_weights <- function(morph_df, anchor_mask = NULL, ...) {
  if (is.logical(anchor_mask) && length(anchor_mask) == nrow(morph_df)) {
    anch_df <- morph_df[anchor_mask, , drop = FALSE]
  } else {
    dots <- list(...)
    if (length(dots) >= 1L && is.logical(dots[[1L]]) && length(dots[[1L]]) == nrow(morph_df)) {
      anch_df <- morph_df[dots[[1L]], , drop = FALSE]
    } else {
      anch_df <- morph_df
    }
  }
  kappa_morph <- if ("cda_predicted_taxon" %in% names(anch_df)) compute_cohens_kappa_stat(anch_df$species_standardized, anch_df$cda_predicted_taxon) else NA_real_
  w_morph <- if (is.finite(kappa_morph)) max(0.05, min(1.0, kappa_morph)) else 0.70
  ed_col <- intersect(c("edaphic_best_fit_taxon", "rf_predicted_taxon"), names(anch_df))[1L]
  kappa_eco <- if (!is.na(ed_col) && any(!is.na(anch_df[[ed_col]]))) compute_cohens_kappa_stat(anch_df$species_standardized, anch_df[[ed_col]]) else NA_real_
  w_eco <- if (is.finite(kappa_eco)) max(0.05, min(1.0, kappa_eco)) else NA_real_
  present_w <- c(morph = w_morph, eco = w_eco); avail <- !is.na(present_w); norm_w <- present_w
  norm_w[avail] <- present_w[avail] / sum(present_w[avail]); norm_w[!avail] <- 0.0
  list(w_morph_raw = w_morph, w_eco_raw = w_eco, w_morph = unname(norm_w["morph"]), w_eco = unname(norm_w["eco"]))
}
# 4. Step C: Dynamic Morphometric & Ecological Affinity Matrices
.parse_step08_semicolon <- function(df, col_nm) {
  if (!(col_nm %in% names(df))) return(character(0))
  raw <- trimws(as.character(stats::na.omit(df[[col_nm]]))); raw <- raw[nzchar(raw)]
  if (length(raw) == 0L) return(character(0))
  toks <- unique(trimws(unlist(strsplit(raw, ";", fixed = TRUE))))
  intersect(toks[nzchar(toks) & !(tolower(toks) %in% c("na", "none", "null"))], names(df))
}
.sort_step08_pcs <- function(cols, pfx) {
  if (length(cols) <= 1L) return(cols)
  cols[order(suppressWarnings(as.integer(sub(paste0("^", pfx), "", cols))), cols)]
}
resolve_step08_morph_features <- function(morph_df, anchor_df = morph_df, n_taxa = 2L) {
  max_p <- max(2L, min(15L, nrow(anchor_df) - as.integer(n_taxa) - 1L))
  ok_num <- function(cc) {
    if (!(cc %in% names(anchor_df)) || !(cc %in% names(morph_df)) || is.logical(anchor_df[[cc]]) || is.character(anchor_df[[cc]])) return(FALSE)
    v <- suppressWarnings(as.numeric(anchor_df[[cc]])); v <- v[is.finite(v)]
    length(v) >= 3L && is.finite(stats::sd(v)) && stats::sd(v) > 1e-8
  }
  sel <- .parse_step08_semicolon(morph_df, "selected_morph_features"); sel <- sel[vapply(sel, ok_num, logical(1))]
  if (length(sel) >= 2L) {
    is_pc <- grepl("^(Sym_)?PC[0-9]+$", sel)
    ord_pc <- c(.sort_step08_pcs(grep("^Sym_PC[0-9]+$", sel, value = TRUE), "Sym_PC"),
                .sort_step08_pcs(grep("^PC[0-9]+$", sel, value = TRUE), "PC"))
    return(utils::head(unique(c(ord_pc, sel[!is_pc])), max_p))
  }
  sym_cols <- .sort_step08_pcs(grep("^Sym_PC[0-9]+$", names(anchor_df), value = TRUE), "Sym_PC"); sym_cols <- sym_cols[vapply(sym_cols, ok_num, logical(1))]
  pc_cols <- .sort_step08_pcs(grep("^PC[0-9]+$", names(anchor_df), value = TRUE), "PC"); pc_cols <- pc_cols[vapply(pc_cols, ok_num, logical(1))]
  can_cols <- intersect(c("can1", "can2", "can3"), names(anchor_df)); can_cols <- can_cols[vapply(can_cols, ok_num, logical(1))]
  cands <- if (length(sym_cols) >= 2L) sym_cols else if (length(pc_cols) >= 2L) pc_cols else unique(c(can_cols, sym_cols, pc_cols))
  if (length(cands) <= 2L) return(cands)
  X_c <- as.matrix(anchor_df[, cands, drop = FALSE]); storage.mode(X_c) <- "double"
  ev <- apply(X_c, 2L, stats::var, na.rm = TRUE); ev[!is.finite(ev) | ev < 0] <- 0.0
  m <- length(ev); ev_prop <- ev / max(sum(ev), 1e-12); bs_prop <- rev(cumsum(1.0 / (m:1L))) / m
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  n_r <- max(4L, nrow(X_c))
  null_mat <- replicate(20L, sort(eigen(stats::cor(matrix(stats::rnorm(n_r * m), nrow = n_r, ncol = m)), symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE) / m)
  pa_q95 <- apply(null_mat, 1L, stats::quantile, probs = 0.95, na.rm = TRUE, names = FALSE)
  keep_k <- max(2L, min(max_p, sum(cumprod(as.integer(ev_prop >= bs_prop | ev_prop >= pa_q95)) == 1L)))
  cands[seq_len(min(length(cands), keep_k))]
}
resolve_step08_eco_predictors <- function(df) {
  excl <- c("edaphic_z_score", "climatic_z_score", "edaphic_schoener_d", "schoener_d", "edaphic_niche_overlap_d",
            "sgmc_lithology_class", "Dominant_Lithology_Class", "edaphic_best_fit_taxon", "rf_predicted_taxon",
            "primary_env_predictors", "primary_clim_predictors", "primary_edaph_predictors")
  ok_env <- function(cands) {
    cands <- setdiff(unique(cands), excl)
    cands[vapply(cands, function(cc) {
      if (!(cc %in% names(df)) || is.logical(df[[cc]]) || is.character(df[[cc]])) return(FALSE)
      v <- suppressWarnings(as.numeric(df[[cc]])); v <- v[is.finite(v)]
      length(v) >= 4L && is.finite(stats::sd(v)) && stats::sd(v) > 1e-8
    }, logical(1))]
  }
  meta_cols <- ok_env(c(.parse_step08_semicolon(df, "primary_env_predictors"), .parse_step08_semicolon(df, "primary_edaph_predictors"), .parse_step08_semicolon(df, "primary_clim_predictors")))
  if (length(meta_cols) > 0L) return(meta_cols)
  cands <- ok_env(grep("^(bio[0-9]+_|fine_|ssurgo_|polaris_|solus_|topo_)", names(df), value = TRUE))
  if (!any(grepl("^fine_", cands))) cands <- unique(c(cands, ok_env(grep("^soil_", names(df), value = TRUE))))
  if (length(cands) <= 1L) return(cands)
  E <- as.matrix(df[, cands, drop = FALSE]); storage.mode(E) <- "double"
  cands <- cands[order(-colSums(is.finite(E)), cands)]; E <- E[, cands, drop = FALSE]; kept <- cands[1L]
  for (cc in cands[-1L]) {
    cors <- vapply(kept, function(kc) { r <- suppressWarnings(stats::cor(E[, cc], E[, kc], use = "pairwise.complete.obs")); if (is.finite(r)) abs(r) else 0.0 }, numeric(1))
    if (all(cors < 0.70)) kept <- c(kept, cc)
  }
  kept
}
compute_ecological_affinity <- function(morph_df, taxa) {
  n_taxa <- length(taxa); S_env <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  D_env <- matrix(1.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa)); diag(D_env) <- 0.0
  eco_cols <- resolve_step08_eco_predictors(morph_df)
  if (length(eco_cols) == 0L) return(list(S_env = S_env, D_env = D_env, has_env = FALSE, eco_predictors = character(0), specimen_schoener_d = rep(NA_real_, nrow(morph_df))))
  E <- as.matrix(morph_df[, eco_cols, drop = FALSE]); storage.mode(E) <- "double"; has_row_env <- rowSums(is.finite(E)) > 0L
  for (j in seq_len(ncol(E))) {
    vj <- E[is.finite(E[, j]), j]; med_j <- if (length(vj) > 0L) stats::median(vj) else 0.0
    sd_j <- if (length(vj) >= 2L && stats::sd(vj) > 1e-8) stats::sd(vj) else 1.0
    E[!is.finite(E[, j]), j] <- med_j; E[, j] <- (E[, j] - med_j) / sd_j
  }
  Z <- if (ncol(E) >= 2L) stats::prcomp(E, center = FALSE, scale. = FALSE)$x[, 1:2, drop = FALSE] else cbind(E[, 1L], 0.0)
  grid_pts <- as.matrix(expand.grid(seq(min(Z[, 1L]) - 1.0, max(Z[, 1L]) + 1.0, length.out = 20L), seq(min(Z[, 2L]) - 1.0, max(Z[, 2L]) + 1.0, length.out = 20L)))
  dens_by_tx <- stats::setNames(lapply(taxa, function(tx) {
    idx <- which(morph_df$species_standardized == tx & has_row_env)
    if (length(idx) == 0L) return(rep(1.0 / nrow(grid_pts), nrow(grid_pts)))
    Zt <- Z[idx, , drop = FALSE]; bw <- max(0.35, 1.06 * mean(apply(Zt, 2L, function(v) if (length(v) >= 2L) max(stats::sd(v), 0.25) else 0.5)) * (length(idx)^(-0.2)))
    p <- rowMeans(exp(-0.5 * (outer(grid_pts[, 1L], Zt[, 1L], "-")^2 + outer(grid_pts[, 2L], Zt[, 2L], "-")^2) / (bw^2))); p / max(sum(p), 1e-12)
  }), taxa)
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j) {
    S_env[i, j] <- max(0.0, min(1.0, 1.0 - 0.5 * sum(abs(dens_by_tx[[i]] - dens_by_tx[[j]])))); D_env[i, j] <- 1.0 - S_env[i, j]
  }
  spec_d <- rep(NA_real_, nrow(morph_df))
  for (tx in taxa) {
    idx <- which(morph_df$species_standardized == tx & has_row_env)
    if (length(idx) > 0L) spec_d[idx] <- exp(-0.5 * rowMeans(sweep(E[idx, , drop = FALSE], 2L, apply(E[idx, , drop = FALSE], 2L, stats::median), "-")^2))
  }
  list(S_env = S_env, D_env = D_env, has_env = sum(has_row_env) >= 4L, eco_predictors = eco_cols, specimen_schoener_d = spec_d)
}
compute_morphometric_affinity <- function(morph_df, anchor_df, taxa) {
  n_taxa <- length(taxa); feat_cols <- resolve_step08_morph_features(morph_df, anchor_df, n_taxa)
  if (length(feat_cols) == 0L) feat_cols <- intersect(c("can1", "can2", "can3"), names(anchor_df))
  for (rc in setdiff(intersect(c("involucre_height_mm", "involucre_width_mm", "capitulum_aspect_ratio", "AR_inv_std", "W_inv_std"), names(anchor_df)), feat_cols)) {
    vals <- suppressWarnings(as.numeric(anchor_df[[rc]]))
    if (sum(is.finite(vals) & vals > 0) >= 0.5 * nrow(anchor_df) && stats::sd(vals, na.rm = TRUE) > 1e-6) {
      col_nm <- paste0("log_", rc); anchor_df[[col_nm]] <- log10(pmax(suppressWarnings(as.numeric(anchor_df[[rc]])), 1e-3))
      morph_df[[col_nm]] <- log10(pmax(suppressWarnings(as.numeric(morph_df[[rc]])), 1e-3)); feat_cols <- c(feat_cols, col_nm)
    }
  }
  X_anchor <- as.matrix(anchor_df[, feat_cols, drop = FALSE]); storage.mode(X_anchor) <- "double"
  X_all <- as.matrix(morph_df[, feat_cols, drop = FALSE]); storage.mode(X_all) <- "double"
  for (j in seq_len(ncol(X_anchor))) {
    med_j <- stats::median(X_anchor[, j], na.rm = TRUE); if (!is.finite(med_j)) med_j <- 0.0
    X_anchor[!is.finite(X_anchor[, j]), j] <- med_j; X_all[!is.finite(X_all[, j]), j] <- med_j
  }
  p <- ncol(X_anchor); centroids <- matrix(0.0, nrow = n_taxa, ncol = p, dimnames = list(taxa, feat_cols))
  anc_type <- if (exists(".identify_type_specimen_vec", mode = "function")) as.logical(.identify_type_specimen_vec(anchor_df)$is_type_specimen) else rep(FALSE, nrow(anchor_df))
  anc_type[is.na(anc_type)] <- FALSE
  S_pooled <- matrix(0.0, nrow = p, ncol = p); df_total <- 0L
  for (tx in taxa) {
    tx_rows <- which(anchor_df$species_standardized == tx); sub_x <- X_anchor[tx_rows, , drop = FALSE]
    if (nrow(sub_x) == 0L) sub_x <- X_all[morph_df$species_standardized == tx, , drop = FALSE]
    mu_base <- colMeans(sub_x, na.rm = TRUE); type_rows_tx <- tx_rows[anc_type[tx_rows]]
    if (length(type_rows_tx) >= 1L) {
      gamma_k <- if (exists(".compute_type_shrinkage_gamma", mode = "function")) .compute_type_shrinkage_gamma(if (exists(".count_effective_type_gatherings", mode = "function")) .count_effective_type_gatherings(anchor_df, type_rows_tx) else length(type_rows_tx), nrow(sub_x)) else min(0.35, (2.5 * length(type_rows_tx)) / (max(1.0, nrow(sub_x)) + 2.5 * length(type_rows_tx)))
      mu_base <- (1.0 - gamma_k) * mu_base + gamma_k * colMeans(X_anchor[type_rows_tx, , drop = FALSE], na.rm = TRUE)
    }
    centroids[tx, ] <- mu_base
    if (nrow(sub_x) >= 2L) {
      cov_tx <- stats::cov(sub_x); cov_tx[!is.finite(cov_tx)] <- 0.0
      S_pooled <- S_pooled + (nrow(sub_x) - 1L) * cov_tx; df_total <- df_total + (nrow(sub_x) - 1L)
    }
  }
  S_pooled <- if (df_total > 0L) S_pooled / df_total else stats::cov(X_anchor); S_pooled[!is.finite(S_pooled)] <- 0.0
  S_reg <- S_pooled + diag(max(1e-4, 0.02 * max(mean(diag(S_pooled), na.rm = TRUE), 1e-4)), nrow = p, ncol = p)
  S_inv <- tryCatch(solve(S_reg), error = function(e) MASS::ginv(S_reg))

  D2_mahal <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j) {
    diff_ij <- centroids[i, ] - centroids[j, ]; D2_mahal[i, j] <- max(0.0, as.numeric(t(diff_ij) %*% S_inv %*% diff_ij))
  }
  off_diag_d2 <- D2_mahal[upper.tri(D2_mahal)]; med_d2 <- stats::median(off_diag_d2[off_diag_d2 > 0], na.rm = TRUE)
  if (!is.finite(med_d2) || med_d2 <= 1e-8) med_d2 <- 1.0
  K_mahal <- exp(-D2_mahal / (2.0 * med_d2)); diag(K_mahal) <- 0.0

  C_morph <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  has_type_nn <- "nearest_type_morph_taxon" %in% names(morph_df)
  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]; idx_i <- which(morph_df$species_standardized == tx_i); if (length(idx_i) == 0L) next
    X_i <- X_all[idx_i, , drop = FALSE]
    cda_pred_i <- if ("cda_predicted_taxon" %in% names(morph_df)) as.character(morph_df$cda_predicted_taxon[idx_i]) else rep(tx_i, length(idx_i))
    type_nn_i <- if (has_type_nn) as.character(morph_df$nearest_type_morph_taxon[idx_i]) else rep(NA_character_, length(idx_i))
    log_likes <- matrix(0.0, nrow = length(idx_i), ncol = n_taxa)
    for (j in seq_len(n_taxa)) { diff_mat <- sweep(X_i, 2, centroids[j, ], "-"); log_likes[, j] <- -0.5 * pmax(0.0, rowSums((diff_mat %*% S_inv) * diff_mat)) }
    probs_i <- exp(sweep(log_likes, 1, apply(log_likes, 1, max), "-")); probs_i <- sweep(probs_i, 1, pmax(rowSums(probs_i), 1e-12), "/")
    for (j in seq_len(n_taxa)) if (i != j) C_morph[i, j] <- mean((!is.na(cda_pred_i) & cda_pred_i == taxa[j]) | (!is.na(type_nn_i) & type_nn_i == taxa[j] & probs_i[, j] >= 0.10) | (probs_i[, j] >= 0.15), na.rm = TRUE)
  }
  C_morph_sym <- 0.5 * (C_morph + t(C_morph)); diag(C_morph_sym) <- 0.0
  S_morph <- 0.60 * K_mahal + 0.40 * C_morph_sym; diag(S_morph) <- 0.0
  list(S_morph = S_morph, K_mahal = K_mahal, C_morph_sym = C_morph_sym, D2_mahal = D2_mahal, morph_features = feat_cols)
}

# 5. Step D: Hopkins-Corrected Circular Anthesis Phenology (Hopkins 1918; Pewsey et al. 2013)
compute_hopkins_phenological_affinity <- function(morph_df, taxa) {
  n_taxa <- length(taxa); n <- nrow(morph_df)
  S_pheno <- matrix(0.0, nrow = n_taxa, ncol = n_taxa, dimnames = list(taxa, taxa))
  if (!("doy" %in% names(morph_df))) return(list(S_pheno = S_pheno, has_pheno = FALSE, delta_doy_resid = rep(NA_real_, n)))
  doy_vec <- suppressWarnings(as.numeric(morph_df$doy))
  anthesis_mask <- is.finite(doy_vec) & doy_vec >= 1 & doy_vec <= 366
  for (pcol in c("voucher_phenological_state", "phenological_Stage", "phenological_state")) {
    if (pcol %in% names(morph_df)) {
      st <- tolower(trimws(as.character(morph_df[[pcol]]))); has_st <- !is.na(st) & nzchar(st) & st != "na"
      if (any(has_st)) anthesis_mask <- anthesis_mask & (!has_st | st %in% c("anthesis", "flowering", "open_florets", "active_anthesis"))
      break
    }
  }
  if (sum(anthesis_mask) < 4L) return(list(S_pheno = S_pheno, has_pheno = FALSE, delta_doy_resid = rep(NA_real_, n)))

  th_pool <- 2.0 * pi * (doy_vec[anthesis_mask] / 365.25)
  ref_theta <- atan2(mean(sin(th_pool)), mean(cos(th_pool)))
  if (!is.finite(ref_theta)) ref_theta <- 2.0 * pi * (120.0 / 365.25)
  if (ref_theta < 0.0) ref_theta <- ref_theta + 2.0 * pi
  ref_doy <- ref_theta * (365.25 / (2.0 * pi))
  doy_clean <- ifelse(is.finite(doy_vec) & doy_vec >= 1 & doy_vec <= 366, doy_vec, ref_doy)

  lat_col <- intersect(c("decimalLatitude", "latitude"), names(morph_df))
  lat_vec <- if (length(lat_col) > 0L) suppressWarnings(as.numeric(morph_df[[lat_col[1L]]])) else rep(NA_real_, n)
  elev_col <- intersect(c("elevation_m", "elevation", "minimumElevationInMeters"), names(morph_df))
  elev_vec <- if (length(elev_col) > 0L) suppressWarnings(as.numeric(morph_df[[elev_col[1L]]])) else rep(NA_real_, n)
  bio1_col <- intersect(c("bio1_temp_mean", "bio1", "mat_c"), names(morph_df))
  bio1_vec <- if (length(bio1_col) > 0L) suppressWarnings(as.numeric(morph_df[[bio1_col[1L]]])) else rep(NA_real_, n)

  med_emp <- function(v, safety) {
    va <- v[anthesis_mask & is.finite(v)]; if (length(va) == 0L) va <- v[is.finite(v)]
    if (length(va) > 0L) stats::median(va, na.rm = TRUE) else safety
  }
  ref_lat <- med_emp(lat_vec, 35.0); ref_elev <- med_emp(elev_vec, 150.0); ref_bio1 <- med_emp(bio1_vec, 14.5)

  lat_shift <- ifelse(is.finite(lat_vec), 4.0 * (lat_vec - ref_lat), 0.0)
  elev_shift <- ifelse(is.finite(elev_vec), 1.25 * ((elev_vec - ref_elev) / 100.0), 0.0)
  bio1_shift <- ifelse(is.finite(bio1_vec), -1.5 * (bio1_vec - ref_bio1), 0.0)
  has_geo <- is.finite(lat_vec) | is.finite(elev_vec); has_mat <- is.finite(bio1_vec); has_any_cov <- has_geo | has_mat
  theo_shift <- ifelse(has_geo & has_mat, 0.70 * (lat_shift + elev_shift) + 0.30 * bio1_shift,
                       ifelse(has_geo, lat_shift + elev_shift, ifelse(has_mat, bio1_shift, 0.0)))

  theta_obs <- 2.0 * pi * (doy_clean / 365.25); sin_doy <- sin(theta_obs); cos_doy <- cos(theta_obs)
  theta_theo <- ref_theta + 2.0 * pi * (theo_shift / 365.25)

  anc_mask <- select_taxon_anchors(morph_df)$anchor_mask
  mo_m <- if ("morphometric_outlier_flag" %in% names(morph_df)) (!is.na(morph_df$morphometric_outlier_flag) & as.logical(morph_df$morphometric_outlier_flag)) else rep(FALSE, n)
  fit_mask <- anthesis_mask & has_any_cov & anc_mask & !mo_m
  if (sum(fit_mask) < 5L) fit_mask <- anthesis_mask & has_any_cov & anc_mask
  if (sum(fit_mask) < 5L) fit_mask <- anthesis_mask & has_any_cov

  cov_raw <- list(decimalLatitude = lat_vec, elevation = elev_vec, bio1_temp_mean = bio1_vec)
  cov_terms <- names(cov_raw)[vapply(cov_raw, function(v) {
    vf <- v[fit_mask & is.finite(v)]; length(vf) >= 3L && is.finite(stats::sd(vf)) && stats::sd(vf) > 1e-6
  }, logical(1))]
  alpha_eb <- if (length(cov_terms) > 0L) min(1.0, max(0.0, (sum(fit_mask) - 5.0) / 15.0)) else 0.0
  theta_blend <- theta_theo; theta_blend_ref <- ref_theta

  if (alpha_eb > 0.0 && length(cov_terms) > 0L) {
    hop_df <- data.frame(sin_doy = sin_doy, cos_doy = cos_doy,
                         decimalLatitude = ifelse(is.finite(lat_vec), lat_vec, ref_lat),
                         elevation = ifelse(is.finite(elev_vec), elev_vec, ref_elev),
                         bio1_temp_mean = ifelse(is.finite(bio1_vec), bio1_vec, ref_bio1))
    ref_df <- data.frame(decimalLatitude = ref_lat, elevation = ref_elev, bio1_temp_mean = ref_bio1)
    rhs <- paste(cov_terms, collapse = " + ")
    lm_s <- tryCatch(stats::lm(stats::as.formula(paste("sin_doy ~", rhs)), data = hop_df[fit_mask, , drop = FALSE]), error = function(e) NULL)
    lm_c <- tryCatch(stats::lm(stats::as.formula(paste("cos_doy ~", rhs)), data = hop_df[fit_mask, , drop = FALSE]), error = function(e) NULL)
    if (!is.null(lm_s) && !is.null(lm_c)) {
      ps <- tryCatch(as.numeric(stats::predict(lm_s, newdata = hop_df)), error = function(e) rep(NA_real_, n))
      pc <- tryCatch(as.numeric(stats::predict(lm_c, newdata = hop_df)), error = function(e) rep(NA_real_, n))
      rs <- tryCatch(as.numeric(stats::predict(lm_s, newdata = ref_df)), error = function(e) NA_real_)
      rc <- tryCatch(as.numeric(stats::predict(lm_c, newdata = ref_df)), error = function(e) NA_real_)
      if (all(is.finite(ps)) && all(is.finite(pc)) && is.finite(rs) && is.finite(rc)) {
        th_emp <- atan2(ps, pc); th_emp_r <- atan2(rs, rc)
        theta_blend <- atan2(alpha_eb * sin(th_emp) + (1.0 - alpha_eb) * sin(theta_theo), alpha_eb * cos(th_emp) + (1.0 - alpha_eb) * cos(theta_theo))
        theta_blend_ref <- atan2(alpha_eb * sin(th_emp_r) + (1.0 - alpha_eb) * sin(ref_theta), alpha_eb * cos(th_emp_r) + (1.0 - alpha_eb) * cos(ref_theta))
      }
    }
  }
  delta_theta_env <- ifelse(has_any_cov, atan2(sin(theta_blend - theta_blend_ref), cos(theta_blend - theta_blend_ref)), 0.0)
  theta_adj <- theta_obs - delta_theta_env

  mu_theta <- numeric(n_taxa); names(mu_theta) <- taxa
  delta_doy_resid <- rep(NA_real_, n)
  for (tx in taxa) {
    idx_tx <- which(morph_df$species_standardized == tx & anthesis_mask)
    if (length(idx_tx) > 0L) {
      mu_tx <- atan2(mean(sin(theta_adj[idx_tx]), na.rm = TRUE), mean(cos(theta_adj[idx_tx]), na.rm = TRUE))
      mu_theta[tx] <- mu_tx
      diff_rad <- atan2(sin(theta_adj[idx_tx] - mu_tx), cos(theta_adj[idx_tx] - mu_tx))
      delta_doy_resid[idx_tx] <- abs(diff_rad) * (365.25 / (2.0 * pi))
    } else mu_theta[tx] <- ref_theta
  }

  .circular_watson_williams_or_perm_p <- function(th1, th2) {
    th1 <- th1[is.finite(th1)]; th2 <- th2[is.finite(th2)]
    n1 <- length(th1); n2 <- length(th2); N <- n1 + n2
    if (n1 < 2L || n2 < 2L) return(1.0)
    R1 <- sqrt(sum(sin(th1))^2 + sum(cos(th1))^2)
    R2 <- sqrt(sum(sin(th2))^2 + sum(cos(th2))^2)
    th_all <- c(th1, th2)
    R_tot <- sqrt(sum(sin(th_all))^2 + sum(cos(th_all))^2)
    rbar <- (R1 + R2) / N
    # Watson-Williams two-sample circular F-statistic (Mardia & Jupp 2000) with permutation verification
    denom_ww <- max(N - R1 - R2, 1e-8)
    num_ww <- max(0.0, R1 + R2 - R_tot)
    k_corr <- if (rbar < 0.95) (1.0 + 3.0 / max(8.0 * max(rbar, 0.15), 1.0)) else 1.0
    F_obs <- k_corr * (N - 2.0) * num_ww / denom_ww
    p_ww <- stats::pf(F_obs, df1 = 1.0, df2 = max(1.0, N - 2.0), lower.tail = FALSE)
    if (N <= 60L) {
      B_perm <- 199L
      s_all <- sin(th_all); c_all <- cos(th_all)
      ge_cnt <- 0L
      for (b in seq_len(B_perm)) {
        idx1 <- sample.int(N, size = n1, replace = FALSE)
        R1_b <- sqrt(sum(s_all[idx1])^2 + sum(c_all[idx1])^2)
        R2_b <- sqrt(sum(s_all[-idx1])^2 + sum(c_all[-idx1])^2)
        if ((R1_b + R2_b - R_tot) >= (num_ww - 1e-12)) ge_cnt <- ge_cnt + 1L
      }
      return((ge_cnt + 1.0) / (B_perm + 1.0))
    }
    pmin(1.0, pmax(0.0, as.numeric(p_ww)))
  }

  raw_p <- numeric(0)
  for (i in seq_len(n_taxa)) for (j in seq_len(n_taxa)) if (i != j) {
    arc_rad <- abs(atan2(sin(mu_theta[i] - mu_theta[j]), cos(mu_theta[i] - mu_theta[j])))
    delta_circ <- arc_rad * (365.25 / (2.0 * pi))
    S_pheno[i, j] <- exp(-(delta_circ^2) / (2.0 * (21.0^2)))
    if (i < j) {
      th_i <- theta_adj[morph_df$species_standardized == taxa[i] & anthesis_mask]
      th_j <- theta_adj[morph_df$species_standardized == taxa[j] & anthesis_mask]
      raw_p <- c(raw_p, .circular_watson_williams_or_perm_p(th_i, th_j))
    }
  }
  pheno_fdr_p <- if (length(raw_p) > 0L) stats::p.adjust(raw_p, method = "BH") else numeric(0)
  diag(S_pheno) <- 0.0
  list(S_pheno = S_pheno, has_pheno = TRUE, sin_doy = sin_doy, cos_doy = cos_doy,
       delta_doy_resid = delta_doy_resid, pheno_fdr_p = pheno_fdr_p, alpha_eb = alpha_eb)
}

# ------------------------------------------------------------------------------
# 6. Step E: Five-Stage Hierarchical Integrative Decision Matrix
#    (Padial et al. 2010; Schlick-Steiner et al. 2010; Rieseberg 1995; Slenker et al. 2022)
# ------------------------------------------------------------------------------
evaluate_hierarchical_complex_delimitation <- function(morph_df, pheno_res = NULL, eco_res = NULL, ...) {
  # Prefer canonical Step 06 `taxonomic_treatment_subtype` whenever precomputed (zero cross-stage decision drift)
  if ("taxonomic_treatment_subtype" %in% names(morph_df)) {
    st_pre <- trimws(as.character(morph_df$taxonomic_treatment_subtype))
    if (length(st_pre) == nrow(morph_df) && all(!is.na(st_pre) & nzchar(st_pre) & st_pre != "NA")) {
      return(st_pre)
    }
  }
  if (!is.list(pheno_res) || !("delta_doy_resid" %in% names(pheno_res))) {
    dots <- list(...)
    if (length(dots) >= 1L && is.list(dots[[1L]]) && ("delta_doy_resid" %in% names(dots[[1L]]))) {
      pheno_res <- dots[[1L]]
      if (length(dots) >= 2L) eco_res <- dots[[2L]]
    }
  }
  n <- nrow(morph_df)
  sp_given <- as.character(morph_df$species_standardized)
  cda_pred <- if ("cda_predicted_taxon" %in% names(morph_df)) as.character(morph_df$cda_predicted_taxon) else sp_given
  cda_prob <- if ("cda_posterior_prob" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_posterior_prob)) else rep(0.90, n)
  cda_prob[!is.finite(cda_prob)] <- 0.90
  cda_sec <- if ("cda_second_posterior_prob" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_second_posterior_prob)) else pmax(0.0, 1.0 - cda_prob)
  cda_sec[!is.finite(cda_sec)] <- pmax(0.0, 1.0 - cda_prob[!is.finite(cda_sec)])
  cda_margin <- if ("cda_posterior_margin" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$cda_posterior_margin)) else abs(cda_prob - cda_sec)
  cda_margin[!is.finite(cda_margin)] <- abs(cda_prob[!is.finite(cda_margin)] - cda_sec[!is.finite(cda_margin)])

  ed_pred <- if ("edaphic_best_fit_taxon" %in% names(morph_df)) as.character(morph_df$edaphic_best_fit_taxon) else if ("rf_predicted_taxon" %in% names(morph_df)) as.character(morph_df$rf_predicted_taxon) else sp_given
  ed_pred[is.na(ed_pred) | !nzchar(ed_pred)] <- sp_given[is.na(ed_pred) | !nzchar(ed_pred)]

  m_flag <- if ("morphometric_outlier_flag" %in% names(morph_df)) as.logical(morph_df$morphometric_outlier_flag) else rep(FALSE, n)
  m_flag[is.na(m_flag)] <- FALSE
  gmm_cl <- if ("gmm_cluster" %in% names(morph_df)) suppressWarnings(as.integer(morph_df$gmm_cluster)) else rep(1L, n)
  gmm_cl[is.na(gmm_cl)] <- 1L
  n_known_sp <- max(2L, length(unique(sp_given)))
  gmm_bf <- if ("gmm_2_delta_bic" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$gmm_2_delta_bic)) else ifelse(gmm_cl > n_known_sp, 8.0, 0.0)
  gmm_bf[!is.finite(gmm_bf)] <- ifelse(gmm_cl[!is.finite(gmm_bf)] > n_known_sp, 8.0, 0.0)
  d_edaph <- if ("edaphic_schoener_d" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$edaphic_schoener_d)) else if ("schoener_d" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$schoener_d)) else rep(NA_real_, n)
  if (!is.null(eco_res) && !is.null(eco_res$specimen_schoener_d)) d_edaph <- ifelse(is.finite(d_edaph), d_edaph, eco_res$specimen_schoener_d)
  d_edaph[!is.finite(d_edaph)] <- 0.65
  d_doy <- if (!is.null(pheno_res) && !is.null(pheno_res$delta_doy_resid)) pheno_res$delta_doy_resid else rep(NA_real_, n)

  circ_sd_d <- function(v) {
    v <- v[is.finite(v)]; if (length(v) < 2L) return(NA_real_)
    ph <- 2.0 * pi * (v / 365.25); rb <- max(1e-6, min(0.999999, sqrt(mean(sin(ph))^2 + mean(cos(ph))^2)))
    (365.25 / (2.0 * pi)) * sqrt(-2.0 * log(rb))
  }
  coh_ph <- d_doy[!m_flag & is.finite(d_doy)]; if (length(coh_ph) < 5L) coh_ph <- d_doy[is.finite(d_doy)]
  sig_coh <- if (length(coh_ph) >= 5L) circ_sd_d(coh_ph) else NA_real_
  tau_eco_vec <- rep(if (is.finite(sig_coh)) min(28.0, max(10.0, 2.0 * sig_coh)) else 14.0, n)
  tau_cry_vec <- rep(if (is.finite(sig_coh)) min(35.0, max(14.0, 2.5 * sig_coh)) else 18.0, n)
  for (tx in unique(sp_given[!is.na(sp_given)])) {
    itx <- which(sp_given == tx); v_tx <- d_doy[itx[!m_flag[itx] & is.finite(d_doy[itx])]]
    if (length(v_tx) < 5L) v_tx <- d_doy[itx[is.finite(d_doy[itx])]]
    if (length(v_tx) >= 5L && is.finite(stx <- circ_sd_d(v_tx))) {
      tau_eco_vec[itx] <- min(28.0, max(10.0, 2.0 * stx)); tau_cry_vec[itx] <- min(35.0, max(14.0, 2.5 * stx))
    }
  }

  is_t1_pure <- rep(TRUE, n)
  if ("leaf_count" %in% names(morph_df)) {
    nt1 <- suppressWarnings(as.numeric(morph_df$leaf_count))
    if (any(!is.na(nt1))) is_t1_pure <- !is.na(nt1) & (nt1 >= 1)
  }

  type_info <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(morph_df) else list(type_rank = rep(0L, n), is_type_specimen = rep(FALSE, n))
  is_type_spec <- as.logical(type_info$is_type_specimen); is_type_spec[is.na(is_type_spec)] <- FALSE
  type_morph_conc <- if ("type_morph_concordance" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$type_morph_concordance)) else rep(NA_real_, n)
  type_edaph_z <- if ("type_edaphic_z_score" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$type_edaphic_z_score)) else rep(NA_real_, n)
  type_concordant_block <- is.finite(type_morph_conc) & (type_morph_conc >= 0.70) & (!is.finite(type_edaph_z) | type_edaph_z <= 2.0)
  is_typ_morph <- if ("is_typical_morphotype" %in% names(morph_df)) as.logical(morph_df$is_typical_morphotype) else rep(NA, n)
  non_typical_morph <- !is.na(is_typ_morph) & !is_typ_morph

  morph_disagrees <- !is.na(cda_pred) & (cda_pred != sp_given) & (cda_prob >= 0.50)
  edaph_disagrees <- !is.na(ed_pred) & (ed_pred != sp_given) & (d_edaph < 0.45)
  both_alt_agree <- morph_disagrees & edaph_disagrees & (cda_pred == ed_pred)
  inter_modal_mosaic <- morph_disagrees & edaph_disagrees & (cda_pred != ed_pred)
  min_cong_d2 <- if ("min_congener_mahal_d2" %in% names(morph_df)) suppressWarnings(as.numeric(morph_df$min_congener_mahal_d2)) else rep(NA_real_, n)
  min_cong_tx <- if ("min_congener_taxon" %in% names(morph_df)) as.character(morph_df$min_congener_taxon) else rep(NA_character_, n)
  congener_osr_reassign <- is.finite(min_cong_d2) & (min_cong_d2 <= stats::qchisq(0.95, df = 3L)) & !is.na(min_cong_tx) & nzchar(min_cong_tx) & (min_cong_tx != sp_given) & m_flag
  single_strong_reassign <- !inter_modal_mosaic & !type_concordant_block & ((morph_disagrees & cda_prob >= 0.80) | congener_osr_reassign)
  # Stage 1: Misidentification_Reassign (evaluated BEFORE hybridity; Nomenclatural Type Immunity Rule A)
  is_stage1_misid <- !is_type_spec & (both_alt_agree | single_strong_reassign)
  # Stage 2: Hybrid_Introgressant (Rieseberg 1995; Slenker et al. 2022)
  is_stage2_hybrid <- !is_type_spec & !is_stage1_misid & (inter_modal_mosaic | (cda_margin < 0.25 & (morph_disagrees | edaph_disagrees) & d_edaph >= 0.30))
  # Stage 3: Putative_Cryptic_Species via Integration by Congruence (Padial et al. 2010)
  has_evol_barrier <- (!is.na(d_doy) & abs(d_doy) >= tau_cry_vec) | (d_edaph < 0.30) | (is.finite(type_edaph_z) & type_edaph_z >= 3.0)
  coh_pass <- if ("cryptic_cluster_cohesion_pass" %in% names(morph_df)) { cp <- as.logical(morph_df$cryptic_cluster_cohesion_pass); cp[is.na(cp)] <- FALSE; cp } else rep(TRUE, n)
  is_stage3_cryptic <- !is_type_spec & !is_stage1_misid & !is_stage2_hybrid & is_t1_pure & coh_pass & (cda_prob < 0.70) &
    (m_flag | non_typical_morph | (gmm_cl > length(unique(sp_given)) & gmm_bf >= 6.0)) & has_evol_barrier
  # Stage 4: Subspecies_Ecotype vs. Ecophenotypic_Plasticity
  is_stage4_ecotype <- !is_type_spec & !is_stage1_misid & !is_stage2_hybrid & !is_stage3_cryptic &
    ((d_edaph < 0.30) | (is.finite(type_edaph_z) & type_edaph_z >= 2.5) | (!is.na(d_doy) & abs(d_doy) >= tau_eco_vec)) & (m_flag | non_typical_morph | gmm_cl > 1L)
  is_stage4_plastic <- !is_type_spec & !is_stage1_misid & !is_stage2_hybrid & !is_stage3_cryptic & !is_stage4_ecotype & m_flag

  ifelse(is_stage1_misid, "Misidentification_Reassign",
    ifelse(is_stage2_hybrid, "Hybrid_Introgressant",
      ifelse(is_stage3_cryptic, "Putative_Cryptic_Species",
        ifelse(is_stage4_ecotype, "Subspecies_Ecotype",
          ifelse(is_stage4_plastic, "Ecophenotypic_Plasticity", "Confirmed_Species")))))
}

# ------------------------------------------------------------------------------
# 7. Step F: Consensus Affinity Matrix & Unsupervised Graph / Ward Partitioning
# ------------------------------------------------------------------------------
relabel_clusters_sequential <- function(cl_vec, taxa) match(cl_vec, unique(cl_vec[order(taxa)]))

partition_taxon_complexes <- function(W, taxa, n_vouchers_vec, n_gs_vec) {
  n_taxa <- length(taxa); diag(W) <- 0.0
  D_mat <- 1.0 - pmin(pmax(W, 0.0), 1.0); diag(D_mat) <- 0.0
  dist_obj <- stats::as.dist(D_mat)

  off_vals <- W[upper.tri(W)]
  tau <- min(0.03, stats::quantile(off_vals[off_vals > 0], probs = 0.50, na.rm = TRUE))
  if (!is.finite(tau) || tau <= 0) tau <- 0.01
  W_adj <- W; W_adj[W_adj < tau] <- 0.0
  for (i in seq_len(n_taxa)) {
    if (sum(W_adj[i, ] > 0) == 0L && max(W[i, ]) > 0) {
      best_j <- which.max(W[i, ]); W_adj[i, best_j] <- W[i, best_j]; W_adj[best_j, i] <- W[best_j, i]
    }
  }

  leiden_comm <- rep(1L, n_taxa); best_mod <- -Inf
  if (requireNamespace("igraph", quietly = TRUE) && sum(W_adj > 0) > 0) {
    g <- igraph::graph_from_adjacency_matrix(W_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
    for (gm in c(0.6, 0.8, 1.0, 1.2)) {
      cl_try <- tryCatch(suppressWarnings(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight, resolution = gm, n_iterations = 10)),
        error = function(e) tryCatch(suppressWarnings(igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight, resolution_parameter = gm, n_iterations = 10)), error = function(e2) NULL))
      if (!is.null(cl_try)) {
        memb <- as.integer(igraph::membership(cl_try)); k_m <- length(unique(memb))
        q_val <- tryCatch(igraph::modularity(g, memb, weights = igraph::E(g)$weight), error = function(e) -1.0)
        if (k_m >= 2L && k_m <= min(8L, n_taxa - 1L) && is.finite(q_val) && q_val > best_mod) { best_mod <- q_val; leiden_comm <- memb }
      }
    }
    if (length(unique(leiden_comm)) < 2L) {
      lv_try <- tryCatch(igraph::cluster_louvain(g, weights = igraph::E(g)$weight), error = function(e) NULL)
      if (!is.null(lv_try)) leiden_comm <- as.integer(igraph::membership(lv_try))
    }
  }
  leiden_comm <- relabel_clusters_sequential(leiden_comm, taxa)

  hc <- stats::hclust(dist_obj, method = "ward.D2"); k_max <- min(8L, n_taxa - 1L)
  ward_cluster <- rep(1L, n_taxa); best_sil_mean <- -Inf; best_sil_vec <- rep(0.0, n_taxa)
  if (n_taxa >= 3L && k_max >= 2L) {
    k_candidates <- 2L:k_max
    sil_scores <- vapply(k_candidates, function(kv) mean(cluster::silhouette(stats::cutree(hc, k = kv), dist_obj)[, "sil_width"], na.rm = TRUE), numeric(1))
    best_idx <- which.max(sil_scores); best_k <- k_candidates[best_idx]; best_sil_mean <- sil_scores[best_idx]
    ward_cluster <- relabel_clusters_sequential(as.integer(stats::cutree(hc, k = best_k)), taxa)
    best_sil_vec <- as.numeric(cluster::silhouette(ward_cluster, dist_obj)[, "sil_width"])
  }

  leiden_sil_mean <- -Inf; leiden_sil_vec <- best_sil_vec
  if (length(unique(leiden_comm)) >= 2L && length(unique(leiden_comm)) < n_taxa) {
    l_sil <- cluster::silhouette(leiden_comm, dist_obj)
    leiden_sil_mean <- mean(l_sil[, "sil_width"], na.rm = TRUE); leiden_sil_vec <- as.numeric(l_sil[, "sil_width"])
  }
  if (leiden_sil_mean > best_sil_mean + 0.02) { complex_id <- leiden_comm; sil_width_vec <- leiden_sil_vec }
  else { complex_id <- ward_cluster; sil_width_vec <- best_sil_vec }
  complex_id <- relabel_clusters_sequential(complex_id, taxa)

  complex_labels <- character(n_taxa); deg_vec <- rowSums(W)
  for (cid in sort(unique(complex_id))) {
    m_idx <- which(complex_id == cid)
    ord <- m_idx[order(-deg_vec[m_idx], -n_vouchers_vec[m_idx], taxa[m_idx])]
    top_ep <- gsub("\\s*\\[K([0-9]+)\\]", "_K\\1", sub("^(?:Packera|P\\.)\\s+", "", taxa[ord[seq_len(min(2L, length(ord)))]]))
    complex_labels[m_idx] <- paste0("Complex_", cid, "_", paste(top_ep, collapse = "_"))
  }

  top_n1 <- character(n_taxa); top_a1 <- numeric(n_taxa); top_n2 <- character(n_taxa); top_a2 <- numeric(n_taxa); within_mean_aff <- numeric(n_taxa)
  for (i in seq_len(n_taxa)) {
    w_row <- W[i, ]; w_row[i] <- -1.0; ord_n <- order(w_row, decreasing = TRUE)
    top_n1[i] <- if (length(ord_n) >= 1L) taxa[ord_n[1L]] else ""; top_a1[i] <- if (length(ord_n) >= 1L) round(max(0.0, W[i, ord_n[1L]]), 4) else 0.0
    top_n2[i] <- if (length(ord_n) >= 2L) taxa[ord_n[2L]] else ""; top_a2[i] <- if (length(ord_n) >= 2L) round(max(0.0, W[i, ord_n[2L]]), 4) else 0.0
    peers_i <- which(complex_id == complex_id[i] & seq_len(n_taxa) != i)
    within_mean_aff[i] <- if (length(peers_i) > 0L) round(mean(W[i, peers_i], na.rm = TRUE), 4) else 0.0
  }

  complexes_df <- data.frame(
    species_standardized = taxa, taxon = taxa,
    n_vouchers = as.integer(n_vouchers_vec[taxa]), n_gold_silver_anchors = as.integer(n_gs_vec[taxa]), n_Pure_Gold_Silver = as.integer(n_gs_vec[taxa]),
    complex_id = as.integer(complex_id), consensus_complex_id = as.integer(complex_id), complex_label = complex_labels,
    leiden_community = as.integer(leiden_comm), leiden_complex_id = as.integer(leiden_comm),
    ward_silhouette_cluster = as.integer(ward_cluster), ward_complex_id = as.integer(ward_cluster),
    silhouette_width = round(sil_width_vec, 4), within_complex_mean_affinity = within_mean_aff,
    top_neighbor_1 = top_n1, nearest_neighbor_taxon = top_n1, top_neighbor_1_affinity = top_a1, nearest_neighbor_affinity = top_a1,
    top_neighbor_2 = top_n2, top_neighbor_2_affinity = top_a2, stringsAsFactors = FALSE
  )
  list(complexes_df = complexes_df, hc = hc, dist_obj = dist_obj)
}

# 8. Step G: Taxon x Morphotype Subcluster Builder & Publication Figure Wrapper
build_taxon_morphotype_subclusters <- function(morph_df, min_vouchers = 5L) {
  raw_sp <- as.character(morph_df$species_standardized)
  u_sp <- sort(unique(raw_sp[!is.na(raw_sp) & nzchar(raw_sp)])); n_sp <- length(u_sp)
  has_gmm <- "gmm_cluster" %in% names(morph_df) && length(unique(stats::na.omit(morph_df$gmm_cluster))) >= 2L
  within_sp_gmm <- if (has_gmm) any(vapply(u_sp, function(sp) {
    ks <- stats::na.omit(as.integer(morph_df$gmm_cluster[raw_sp == sp]))
    length(unique(ks)) >= 2L && sum(table(ks) >= min_vouchers) >= 2L
  }, logical(1))) else FALSE
  if (n_sp > 3L || (n_sp > 2L && (!within_sp_gmm || nrow(morph_df) > 100L))) {
    return(list(morph_df = morph_df, subcluster_mode = FALSE))
  }

  k_vec <- if (has_gmm) ifelse(is.na(morph_df$gmm_cluster), 1L, as.integer(morph_df$gmm_cluster)) else rep(1L, nrow(morph_df))
  short_sp <- sub("^Packera\\s+", "P. ", raw_sp); cand_lbl <- sprintf("%s [K%d]", short_sp, k_vec)
  if (sum(table(cand_lbl) >= min_vouchers) < 3L && nrow(morph_df) >= 3L * min_vouchers) {
    fcols <- resolve_step08_morph_features(morph_df, morph_df, max(2L, n_sp))
    if (length(fcols) == 0L) fcols <- intersect(c("can1", "can2", "can3"), names(morph_df))
    X_m <- scale(as.matrix(morph_df[, fcols, drop = FALSE])); X_m[!is.finite(X_m)] <- 0.0
    for (sp in u_sp) {
      idx_s <- which(raw_sp == sp)
      k_s <- if (n_sp == 1L) min(3L, max(2L, floor(length(idx_s) / min_vouchers))) else min(2L, max(1L, floor(length(idx_s) / min_vouchers)))
      if (k_s >= 2L && length(unique(k_vec[idx_s])) < k_s) {
        km <- tryCatch(stats::kmeans(X_m[idx_s, , drop = FALSE], centers = k_s, nstart = 10L), error = function(e) NULL)
        if (!is.null(km)) k_vec[idx_s] <- as.integer(km$cluster)
      }
    }
    cand_lbl <- sprintf("%s [K%d]", short_sp, k_vec)
  }
  tab_c <- table(cand_lbl); valid_c <- names(tab_c)[tab_c >= min_vouchers]
  if (length(valid_c) < 2L) valid_c <- names(tab_c)[tab_c >= max(1L, min(tab_c))]
  keep_m <- cand_lbl %in% valid_c; morph_df <- morph_df[keep_m, , drop = FALSE]
  morph_df$base_species <- raw_sp[keep_m]; morph_df$species_standardized <- cand_lbl[keep_m]
  if ("cda_predicted_taxon" %in% names(morph_df)) {
    cp <- std_taxon(morph_df$cda_predicted_taxon); cp_short <- sub("^Packera\\s+", "P. ", cp)
    morph_df$cda_predicted_taxon <- vapply(seq_len(nrow(morph_df)), function(i) {
      if (!is.na(cp[i]) && cp[i] %in% valid_c) return(cp[i])
      if (!is.na(cp[i]) && cp[i] == morph_df$base_species[i]) return(morph_df$species_standardized[i])
      m_alt <- valid_c[startsWith(valid_c, paste0(cp_short[i], " ["))]
      if (length(m_alt) > 0L) m_alt[1L] else morph_df$species_standardized[i]
    }, character(1))
  }
  message(sprintf("[Subcluster Mode] Built %d Taxon x Morphotype clusters (>= %d vouchers): %s", length(valid_c), min_vouchers, paste(valid_c, collapse = ", ")))
  list(morph_df = morph_df, subcluster_mode = TRUE)
}

render_complex_discovery_figure <- function(W, part_res, morph_res, out_pdf) {
  if (exists("render_step08_complex_discovery_pdf", mode = "function")) return(render_step08_complex_discovery_pdf(W = W, part_res = part_res, out_pdf = out_pdf))
  invisible(out_pdf)
}

# 9. Main Execution Entrypoint
main <- function() {
  set.seed(42); opts <- parse_discovery_args()
  message("=== Phase 08: Unsupervised Quantitative Taxon-Complex Discovery Engine ===")
  inputs <- load_and_merge_inputs(opts)
  sub_res <- build_taxon_morphotype_subclusters(inputs$morph_df, min_vouchers = 5L)
  morph_df <- sub_res$morph_df
  taxa <- sort(unique(morph_df$species_standardized)); n_taxa <- length(taxa)
  if (n_taxa < 2L) stop("At least 2 distinct Packera taxa or morphotype clusters are required for discovery (found ", n_taxa, ").")

  n_vouchers_vec <- table(morph_df$species_standardized)[taxa]
  anc <- select_taxon_anchors(morph_df)
  kappa_wts <- compute_modality_kappa_weights(morph_df, anc$anchor_mask)
  message(sprintf("Evaluated %d vouchers across %d Packera groups (%d Basal_Rosette Gold/Silver anchors; kappa weights: morph=%.3f, eco=%.3f).",
                  nrow(morph_df), n_taxa, nrow(anc$anchor_df), kappa_wts$w_morph, kappa_wts$w_eco))

  morph_res <- compute_morphometric_affinity(morph_df, anc$anchor_df, taxa)
  pheno_res <- compute_hopkins_phenological_affinity(morph_df, taxa)
  eco_res <- compute_ecological_affinity(morph_df, taxa)
  specimen_subtypes <- evaluate_hierarchical_complex_delimitation(morph_df, pheno_res, eco_res)

  W <- morph_res$S_morph
  diag(W) <- 0.0; W <- round(W, 6)

  part_res <- partition_taxon_complexes(W = W, taxa = taxa, n_vouchers_vec = n_vouchers_vec, n_gs_vec = anc$n_gold_silver)
  cdf <- part_res$complexes_df
  hop_ov <- numeric(n_taxa); eco_ov <- rep(NA_real_, n_taxa); int_cong <- logical(n_taxa)
  n_misid <- integer(n_taxa); n_hyb <- integer(n_taxa); n_crypt <- integer(n_taxa)
  for (i in seq_len(n_taxa)) {
    tx_i <- taxa[i]; idx_tx <- which(morph_df$species_standardized == tx_i)
    n_misid[i] <- sum(specimen_subtypes[idx_tx] == "Misidentification_Reassign", na.rm = TRUE)
    n_hyb[i] <- sum(specimen_subtypes[idx_tx] == "Hybrid_Introgressant", na.rm = TRUE)
    n_crypt[i] <- sum(specimen_subtypes[idx_tx] == "Putative_Cryptic_Species", na.rm = TRUE)
    nn_tx <- cdf$nearest_neighbor_taxon[i]
    if (nzchar(nn_tx) && nn_tx %in% taxa) {
      hop_ov[i] <- if (isTRUE(pheno_res$has_pheno)) round(pheno_res$S_pheno[tx_i, nn_tx], 4) else NA_real_
      eco_ov[i] <- if (isTRUE(eco_res$has_env)) round(eco_res$S_env[tx_i, nn_tx], 4) else NA_real_
      n_mod <- as.integer(morph_res$S_morph[tx_i, nn_tx] >= opts$min_affinity) +
               as.integer(isTRUE(pheno_res$has_pheno) && pheno_res$S_pheno[tx_i, nn_tx] >= 0.25) +
               as.integer(isTRUE(eco_res$has_env) && eco_res$S_env[tx_i, nn_tx] >= 0.25)
      int_cong[i] <- n_mod >= 2L
    }
  }
  part_res$complexes_df$n_type_specimens <- if (!is.null(anc$n_type_specimens)) as.integer(anc$n_type_specimens[taxa]) else integer(n_taxa)
  part_res$complexes_df$hopkins_pheno_overlap <- hop_ov
  part_res$complexes_df$eco_schoener_overlap <- eco_ov
  part_res$complexes_df$integrative_congruence <- int_cong
  part_res$complexes_df$kappa_morph_weight <- round(kappa_wts$w_morph, 4)
  part_res$complexes_df$kappa_eco_weight <- round(kappa_wts$w_eco, 4)
  part_res$complexes_df$n_misidentification_reassign <- n_misid
  part_res$complexes_df$n_hybrid_introgressant <- n_hyb
  part_res$complexes_df$n_putative_cryptic <- n_crypt

  atomic_write_csv(part_res$complexes_df, opts$output_complexes)
  message("Saved empirical taxon complexes -> ", opts$output_complexes)
  aff_df <- data.frame(taxon = taxa, species_standardized = taxa, as.data.frame(W, check.names = FALSE), stringsAsFactors = FALSE, check.names = FALSE)
  atomic_write_csv(aff_df, opts$output_affinity)
  message("Saved pairwise taxon affinity matrix -> ", opts$output_affinity)

  focal_raw <- trimws(as.character(opts$focal_taxon))
  if (nzchar(focal_raw)) {
    focal_std <- std_taxon(focal_raw); focal_short <- sub("^Packera\\s+", "P. ", focal_std)
    focal_nodes <- taxa[taxa == focal_std | startsWith(taxa, paste0(focal_short, " ["))]
    if (length(focal_nodes) > 0L) {
      cdf <- part_res$complexes_df; focal_cids <- unique(cdf$complex_id[cdf$species_standardized %in% focal_nodes])
      same_complex <- cdf$species_standardized[cdf$complex_id %in% focal_cids]
      aff_mask <- apply(morph_res$C_morph_sym[focal_nodes, , drop = FALSE] >= opts$min_affinity, 2, any)
      selected_neighbors <- sort(setdiff(union(same_complex, taxa[aff_mask]), focal_nodes))
    } else selected_neighbors <- character(0)
    atomic_write_lines(selected_neighbors, opts$output_neighbors)
    message("Saved focal auto-neighbors (", length(selected_neighbors), " taxa) -> ", opts$output_neighbors)
    cat(paste(selected_neighbors, collapse = "\n"), "\n", sep = "")
  } else atomic_write_lines(part_res$complexes_df$species_standardized, opts$output_neighbors)

  render_complex_discovery_figure(W, part_res, morph_res, opts$output_figure)
  message("Saved taxon-complex discovery figure -> ", opts$output_figure)
}

if (sys.nframe() == 0L) main()
