#!/usr/bin/env Rscript
# Script: 03_fourier_extractor.R - Packera Morphometrics Pipeline (NCU)
# Homologized dynamic-harmonic EFA (Momocs::efourier), ECT ordination, OTU-anchored weighted PCA & basal median aggregation.

if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)
suppressPackageStartupMessages({ for (.pkg in c("Momocs", "RRPP", "vegan", "ggplot2", "dplyr", "readr", "tibble", "optparse")) if (requireNamespace(.pkg, quietly = TRUE)) library(.pkg, character.only = TRUE) })
.s03_dir <- dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile))))
for (.sd in list(c("scripts/morphometrics/step03_config_scale_utils.R", "step03_config_scale_utils.R", "load_step03_fourier_config"),
  c("scripts/morphometrics/pub_figures_utils.R", "pub_figures_utils.R", "generate_step03_pub_figures"),
  c("scripts/morphometrics/authority_weighting_utils.R", "authority_weighting_utils.R", "compute_authority_precision_weights"),
  c("scripts/morphometrics/dissected_leaf_utils.R", "dissected_leaf_utils.R", "compute_dissection_index_from_coo"),
  c("scripts/core/db_utils.R", file.path("..", "core", "db_utils.R"), "persist_step03_to_sqlite"),
  c("scripts/morphometrics/fourier_parallel_utils.R", "fourier_parallel_utils.R", "align_contour_horizontal_xax"))) {
  if (!exists(.sd[3L], mode = "function")) for (.p in c(.sd[1L], file.path(.s03_dir, .sd[2L]))) if (!is.na(.p) && nzchar(.p) && file.exists(.p)) { source(.p); break }
}

# 1. CLI Argument Parsing
parse_args_robust <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  cli_args <- if (exists("strip_deprecated_step03_cli_args", mode = "function")) strip_deprecated_step03_cli_args(cli_args) else as.character(cli_args)
  has_cfg <- any(grepl("^--config(=|$)", cli_args))
  e1 <- Sys.getenv("PLANTMATICS_CONFIG_PATH", ""); e2 <- Sys.getenv("PACKERA_CONFIG_PATH", "")
  e_cfg <- if (!has_cfg && nzchar(e1)) e1 else if (!has_cfg && nzchar(e2)) e2 else "config/config.yaml"
  specs <- list(
    list(c("-i", "--input"), "input", "character", "data/contours/"), list(c("-c", "--contours-dir"), "contours_dir", "character", "data/contours/"),
    list(c("-f", "--manifest"), "manifest", "character", "data/tables/extracted_leaf_manifest.csv"), list(c("-m", "--masks-dir"), "masks_dir", "character", "data/masks/"),
    list(c("-v", "--vouchers"), "vouchers", "character", "data/tables/curated_vouchers.csv"), list(c("-o", "--output"), "output", "character", "data/tables/leaf_efa_harmonics.csv"),
    list(c("--output-individual", "--output-ind"), "output_individual", "character", "data/tables/leaf_efa_harmonics_individual.csv"),
    list(c("--figures-dir"), "figures_dir", "character", ""), list(c("--config"), "config", "character", e_cfg),
    list(c("-k", "--harmonics"), "harmonics", "integer", 12L), list(c("-p", "--num-pcs"), "num_pcs", "integer", 5L),
    list(c("--symmetrize-efa-harmonics"), "symmetrize_efa_harmonics", "logical", TRUE), list(c("--seed"), "seed", "integer", 42L)
  )
  if (requireNamespace("optparse", quietly = TRUE)) {
    opts <- optparse::parse_args(optparse::OptionParser(usage = "%prog [options]",
      option_list = lapply(specs, function(s) optparse::make_option(s[[1L]][s[[1L]] != "--output-ind"], type = s[[3L]], default = s[[4L]]))), args = cli_args)
    names(opts) <- gsub("-", "_", names(opts))
  } else {
    opts <- stats::setNames(lapply(specs, `[[`, 4L), vapply(specs, `[[`, character(1L), 2L)); i <- 1L
    while (i <= length(cli_args)) {
      tok <- cli_args[i]; eq <- regexpr("=", tok, fixed = TRUE); arg <- if (eq > 0L) substr(tok, 1L, eq - 1L) else tok
      val <- if (eq > 0L) substr(tok, eq + 1L, nchar(tok)) else if (i < length(cli_args)) cli_args[i + 1L] else ""
      hit <- which(vapply(specs, function(s) arg %in% s[[1L]], logical(1L)))
      if (length(hit) >= 1L) { s <- specs[[hit[1L]]]; tp <- s[[3L]]; opts[[s[[2L]]]] <- if (tp == "integer") as.integer(val) else if (tp == "logical") as.logical(val) else val; i <- i + if (eq > 0L) 1L else 2L } else i <- i + 1L
    }
  }
  if (!is.null(opts$contours_dir) && opts$contours_dir != "data/contours/") opts$input <- opts$contours_dir else if (!is.null(opts$input) && opts$contours_dir == "data/contours/") opts$contours_dir <- opts$input
  ex <- list()
  if (any(grepl("^(--harmonics(=|$)|-k)", cli_args))) ex$harmonics <- opts$harmonics
  if (any(grepl("^(--num-pcs(=|$)|-p)", cli_args))) ex$num_pcs <- opts$num_pcs
  if (any(grepl("^--symmetrize-efa-harmonics(=|$)", cli_args))) ex$symmetrize_efa_harmonics <- opts$symmetrize_efa_harmonics
  if (any(grepl("^--seed(=|$)", cli_args))) ex$seed <- opts$seed
  if (exists("load_step03_fourier_config", mode = "function")) {
    fc <- load_step03_fourier_config(cfg_path = opts$config, cli_opts = ex)
    opts$harmonics <- fc$harmonics; opts$num_pcs <- fc$num_pcs; opts$symmetrize_efa_harmonics <- fc$symmetrize_efa_harmonics; opts$seed <- fc$random_seed
  }
  attr(opts, "explicit_cli") <- ex
  attr(opts, "parsed_snapshot") <- list(harmonics = opts$harmonics, num_pcs = opts$num_pcs, symmetrize_efa_harmonics = opts$symmetrize_efa_harmonics, seed = opts$seed)
  opts
}

# 2. Pure-R Mathematical Engine: Homologized EFA (Momocs / Kuhl & Giardina Fallback)
compute_efourier_coords_pure <- function(coo, nb_h = 12, norm = FALSE, start = FALSE) {
  if (is.null(coo) || nrow(coo) < 5) return(rep(NA_real_, nb_h * 4))
  coo <- sweep(coo, 2, colMeans(coo, na.rm = TRUE), "-"); cs <- sqrt(sum(coo[, 1]^2 + coo[, 2]^2, na.rm = TRUE)); if (is.finite(cs) && cs > 1e-8) coo <- coo / cs
  if (coo[1, 1] != coo[nrow(coo), 1] || coo[1, 2] != coo[nrow(coo), 2]) coo <- rbind(coo, coo[1, ])
  dx <- diff(coo[, 1]); dy <- diff(coo[, 2]); dt <- sqrt(dx^2 + dy^2); v <- dt > 1e-7; if (sum(v) < 5) return(rep(NA_real_, nb_h * 4))
  dx <- dx[v]; dy <- dy[v]; dt <- dt[v]; tv <- c(0, cumsum(dt)); Tp <- tv[length(tv)]; if (Tp <= 0) return(rep(NA_real_, nb_h * 4))
  w0 <- 2 * pi / Tp; A <- B <- C <- D <- numeric(nb_h)
  for (n in seq_len(nb_h)) {
    cf <- Tp / (2 * (n^2) * (pi^2)); cd <- cos(n * w0 * tv[-1]) - cos(n * w0 * tv[-length(tv)]); sd <- sin(n * w0 * tv[-1]) - sin(n * w0 * tv[-length(tv)])
    A[n] <- cf * sum((dx / dt) * cd); B[n] <- cf * sum((dx / dt) * sd); C[n] <- cf * sum((dy / dt) * cd); D[n] <- cf * sum((dy / dt) * sd)
  }
  if (!norm) return(as.vector(rbind(A, B, C, D)))
  th1 <- if (start) 0.5 * atan2(2 * (A[1] * B[1] + C[1] * D[1]), (A[1]^2 + C[1]^2 - B[1]^2 - D[1]^2)) else 0.0
  a1s <- A[1] * cos(th1) + B[1] * sin(th1); c1s <- C[1] * cos(th1) + D[1] * sin(th1); psi1 <- atan2(c1s, a1s); E1 <- sqrt(a1s^2 + c1s^2); if (E1 < 1e-7) return(rep(NA_real_, nb_h * 4))
  Rp <- matrix(c(cos(psi1), sin(psi1), -sin(psi1), cos(psi1)), 2, byrow = TRUE); nA <- nB <- nC <- nD <- numeric(nb_h)
  for (n in seq_len(nb_h)) { nM <- (1 / E1) * (Rp %*% matrix(c(A[n], B[n], C[n], D[n]), 2, byrow = TRUE) %*% matrix(c(cos(n * th1), -sin(n * th1), sin(n * th1), cos(n * th1)), 2, byrow = TRUE)); nA[n] <- nM[1, 1]; nB[n] <- nM[1, 2]; nC[n] <- nM[2, 1]; nD[n] <- nM[2, 2] }
  as.vector(rbind(nA, nB, nC, nD))
}

compute_polygon_solidity <- function(coo, area_px) {
  if (is.null(coo) || nrow(coo) < 3 || is.na(area_px) || area_px <= 0) return(NA_real_)
  h <- tryCatch(grDevices::chull(coo[, 1], coo[, 2]), error = function(e) NULL); if (is.null(h) || length(h) < 3) return(NA_real_)
  hx <- coo[h, 1]; hy <- coo[h, 2]; ha <- abs(0.5 * sum(hx * c(hy[-1], hy[1]) - hy * c(hx[-1], hx[1])))
  if (is.finite(ha) && ha > 1e-6) round(min(1.0, area_px / ha), 4) else NA_real_
}

compute_centsize_robust <- function(coo) {
  if (is.null(coo) || nrow(coo) < 3) return(NA_real_)
  if (requireNamespace("Momocs", quietly = TRUE)) { cs <- tryCatch(Momocs::coo_centsize(coo), error = function(e) NA_real_); if (!is.na(cs) && is.finite(cs) && cs > 0) return(cs) }
  val <- sqrt(sum((coo[, 1] - mean(coo[, 1], na.rm = TRUE))^2 + (coo[, 2] - mean(coo[, 2], na.rm = TRUE))^2, na.rm = TRUE)); if (is.finite(val) && val > 0) val else NA_real_
}

verify_harmonic_power_calibration <- function(out_obj = NULL, harm_mat = NULL, nb_h = 12L, min_h = 6L, max_h = 24L, target_power = 0.99) {
  min_h <- max(1L, as.integer(min_h)); max_h <- max(min_h, as.integer(max_h))
  fb_h <- max(min_h, min(max_h, as.integer(if (is.null(nb_h) || is.na(nb_h)) 12L else nb_h))); h_opt <- fb_h; cum_pct <- NA_real_
  if (exists("calibrate_dissection_aware_harmonics", mode = "function")) {
    fmls <- names(formals(calibrate_dissection_aware_harmonics)); ca <- list(out_obj = out_obj, harm_mat = harm_mat, nb_h = nb_h, min_h = min_h, max_h = max_h)
    if ("target_power" %in% fmls || "..." %in% fmls) ca$target_power <- target_power
    cd <- do.call(calibrate_dissection_aware_harmonics, ca); h_opt <- cd$h_opt; cum_pct <- cd$cum_pct
  } else if (!is.null(out_obj) && requireNamespace("Momocs", quietly = TRUE)) {
    cr <- tryCatch(Momocs::calibrate_harmonicpower_efourier(if (length(out_obj$coo) > 150L) Momocs::slice(out_obj, seq_len(150L)) else out_obj, nb.h = max_h, plot = FALSE), error = function(e) NULL)
    if (!is.null(cr) && length(cr$q) > 0L) h_opt <- fb_h
  }
  if ((is.na(cum_pct) || !is.finite(cum_pct)) && !is.null(harm_mat) && ncol(harm_mat) >= 4L) {
    hp <- vapply(seq_len(max(suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(harm_mat)))), fb_h, na.rm = TRUE)), function(n) { cn <- intersect(paste0(c("A", "B", "C", "D"), n), colnames(harm_mat)); if (length(cn) > 0L) mean(rowSums(harm_mat[, cn, drop = FALSE]^2, na.rm = TRUE) / 2.0, na.rm = TRUE) else 0.0 }, numeric(1))
    te <- sum(hp) + tail(hp, 1) * 0.15; cum_pct <- if (te > 1e-12) min(99.99, (sum(hp[seq_len(min(h_opt, length(hp)))]) / te) * 100.0) else 99.50
  }
  if (is.na(cum_pct) || !is.finite(cum_pct)) cum_pct <- 99.50
  thr <- if (is.finite(target_power) && target_power > 0 && target_power <= 1.0) target_power * 100.0 else 99.0
  message(sprintf("Harmonic Power Calibration (Bonhomme et al. 2014): H_opt = %d captures %.2f%% cumulative harmonic power (threshold >= %.1f%%: %s).", as.integer(h_opt), cum_pct, thr, ifelse(cum_pct >= thr, "PASS", "WARN")))
  as.integer(h_opt)
}

prune_within_taxon_mahalanobis <- function(X, df, ref_mask = NULL, q = 0.95, min_anchors_per_taxon = 5L, taxon_col = "species_standardized") {
  elig <- if (!is.null(ref_mask) && length(ref_mask) == nrow(df)) { m <- as.logical(ref_mask); m[is.na(m)] <- FALSE; m } else rep(TRUE, nrow(df))
  ti <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(df) else list(is_type_specimen = rep(FALSE, nrow(df)))
  if (sum(elig) < 2L || is.null(X) || ncol(X) < 1L) {
    attr(elig, "is_promoted_bronze") <- rep(FALSE, nrow(df)); attr(elig, "is_type_specimen") <- as.logical(ti$is_type_specimen)
    attr(elig, "anchor_tier_source") <- ifelse(ti$is_type_specimen, "Tier_0_Type_Anchor", "Passive_Unpromoted"); return(elig)
  }
  anc <- select_adaptive_reference_anchors(df[elig, , drop = FALSE], feature_matrix = as.matrix(X[elig, , drop = FALSE]),
    taxon_col = if (taxon_col %in% names(df)) taxon_col else "scientificName", min_anchors_per_taxon = min_anchors_per_taxon, mahal_q_gold_silver = q, mahal_q_bronze_promo = 0.80)
  out <- rep(FALSE, nrow(df)); out[elig] <- anc$active_mask
  attr(out, "is_promoted_bronze") <- { p <- rep(FALSE, nrow(df)); p[elig] <- anc$is_promoted_bronze; p }; attr(out, "is_type_specimen") <- as.logical(ti$is_type_specimen)
  attr(out, "anchor_tier_source") <- { a <- rep("Passive_Unpromoted", nrow(df)); a[elig] <- as.character(anc$anchor_tier_source); a }; out
}

.fit_weighted_morphospace_pca <- function(X_mat, fit_idx, w_vec, target_var = 0.95) {
  cv <- apply(X_mat[fit_idx, , drop = FALSE], 2L, stats::var, na.rm = TRUE); cols <- names(cv[!is.na(cv) & cv > 1e-10]); if (length(cols) < 1L) return(NULL)
  w_cl <- w_vec[fit_idx]; if (sum(w_cl, na.rm = TRUE) <= 0) w_cl <- rep(1.0, length(fit_idx))
  wt_cov <- stats::cov.wt(as.matrix(X_mat[fit_idx, cols, drop = FALSE]), wt = w_cl, center = TRUE); eig <- eigen(wt_cov$cov, symmetric = TRUE)
  colnames(eig$vectors) <- paste0("PC", seq_len(ncol(eig$vectors))); rownames(eig$vectors) <- cols
  sd_v <- sqrt(pmax(eig$values, 0.0)); k95 <- which((cumsum(sd_v^2) / max(sum(sd_v^2), 1e-12)) >= target_var)[1L]; if (is.na(k95) || k95 < 1L) k95 <- ncol(eig$vectors)
  list(pca = structure(list(sdev = sd_v, rotation = eig$vectors, center = wt_cov$center, scale = FALSE), class = "prcomp"), cols = cols, k_95 = as.integer(k95))
}

extract_symmetric_fourier_pcs <- function(ef, fac_df = NULL, num_pcs = 5, seed = 42, fit_pca = TRUE, ...) {
  coe <- if (is.list(ef) && !is.null(ef$coe)) as.matrix(ef$coe) else as.matrix(ef)
  if (is.null(fac_df) && is.list(ef) && !is.null(ef$fac)) fac_df <- as.data.frame(ef$fac)
  sr <- grep("^[AD]_?[0-9]+$", colnames(coe), value = TRUE); ar <- grep("^[BC]_?[0-9]+$", colnames(coe), value = TRUE)
  vc <- complete.cases(coe[, if (length(sr) > 0) sr else seq_len(ncol(coe)), drop = FALSE])
  sv <- if (length(sr) > 0) ifelse(vc, rowSums(coe[, sr, drop = FALSE]^2, na.rm = TRUE), NA_real_) else rep(0.0, nrow(coe))
  av <- if (length(ar) > 0) ifelse(vc, rowSums(coe[, ar, drop = FALSE]^2, na.rm = TRUE), NA_real_) else rep(0.0, nrow(coe))
  pa <- ifelse(vc, av / pmax(sv + av, 1e-12), NA_real_)
  out_df <- if (!is.null(fac_df)) as.data.frame(fac_df) else data.frame(row_id = seq_len(nrow(coe)))
  ts <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(out_df)$is_type_specimen else rep(FALSE, nrow(out_df))
  sol <- if ("solidity" %in% names(out_df)) ifelse(is.na(out_df$solidity), 1.0, as.numeric(out_df$solidity)) else rep(1.0, nrow(out_df))
  cc <- coe; cc[!is.finite(cc)] <- 0.0; rc <- compute_dissection_aware_roughness(cc, solidity_vec = sol, pct_asymmetry_vec = ifelse(is.finite(pa), pa, 0.0))
  iso <- ifelse(vc, rc$is_contour_outlier, FALSE); iso[as.logical(ts)] <- FALSE; coe_s <- coe
  if (requireNamespace("Momocs", quietly = TRUE) && inherits(ef, "OutCoe") && any(vc)) {
    if (all(vc)) { es <- tryCatch(Momocs::rm_asym(ef), error = function(e) NULL); if (!is.null(es$coe)) coe_s <- as.matrix(es$coe) }
    else { se <- tryCatch(Momocs::slice(ef, which(vc)), error = function(e) NULL); ss <- if (!is.null(se)) tryCatch(Momocs::rm_asym(se), error = function(e) NULL) else NULL; if (!is.null(ss$coe)) coe_s[which(vc), colnames(ss$coe)] <- as.matrix(ss$coe) }
  }
  if (length(ar) > 0L && any(vc)) coe_s[which(vc), ar] <- 0.0
  scols <- grep("^[AD]_?[0-9]+$", colnames(coe_s), value = TRUE); smat <- coe_s[, if (length(scols) > 0) scols else colnames(coe_s), drop = FALSE]; if (any(!vc)) smat[!vc, ] <- NA_real_
  aw <- if ("archetype_score" %in% names(out_df)) { a <- suppressWarnings(as.numeric(out_df$archetype_score)); ifelse(is.na(a) | !is.finite(a) | a <= 0, 1.0, pmin(pmax(a, 0.20), 1.0)) } else rep(1.0, nrow(out_df))
  out_df$sym_var <- round(sv, 6); out_df$asym_var <- round(av, 6); out_df$pct_asymmetry <- round(pa, 6)
  out_df$harmonic_roughness <- round(ifelse(vc, rc$harmonic_roughness, NA_real_), 6); out_df$is_contour_outlier <- iso
  out_df$leaf_quality_weight <- round(pmax(0.05, aw * ((1.0 - pmin(pa, 0.95))^2) * ifelse(iso, 0.10, 1.0)), 6)
  if (!isTRUE(fit_pca)) return(list(df = out_df, sym_mat = smat, pca = NULL, k_95 = as.integer(num_pcs), permanova = NULL, active_ref_mask = NULL))

  snum <- smat; snum[!is.finite(snum)] <- 0.0; cm <- vc & (if (sum(vc) >= 10 && sum(vc & !iso) >= 5) !iso else rep(TRUE, nrow(snum)))
  lzs <- if ("leaf_zone" %in% names(out_df)) tolower(trimws(as.character(out_df$leaf_zone))) else rep("basal", nrow(out_df))
  fm <- if (sum(cm & (is.na(lzs) | lzs == "basal"), na.rm = TRUE) >= 2) cm & (is.na(lzs) | lzs == "basal") else cm
  arm <- if (sum(fm, na.rm = TRUE) >= 2L) prune_within_taxon_mahalanobis(snum, out_df, fm, q = 0.95) else fm
  out_df$is_promoted_bronze <- if (!is.null(attr(arm, "is_promoted_bronze"))) attr(arm, "is_promoted_bronze") else rep(FALSE, nrow(out_df))
  out_df$is_type_specimen <- if (!is.null(attr(arm, "is_type_specimen"))) attr(arm, "is_type_specimen") else ts
  out_df$anchor_tier_source <- if (!is.null(attr(arm, "anchor_tier_source"))) attr(arm, "anchor_tier_source") else ifelse(out_df$is_type_specimen, "Tier_0_Type_Anchor", "Passive_Unpromoted")
  tw <- compute_authority_precision_weights(out_df, active_mask = arm); if (sum(arm, na.rm = TRUE) < 2L) arm <- fm
  sw <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(out_df, active_mask = arm, recompute = TRUE) else rep(1.0, nrow(out_df))
  out_df$authority_weight <- tw; out_df$spatial_balance_weight <- round(sw, 6); for (pc in paste0("PC", 1:5)) out_df[[pc]] <- ifelse(vc, 0.0, NA_real_)
  pca_s <- NULL; k95 <- 5L; fr <- if (sum(arm, na.rm = TRUE) >= 2L) .fit_weighted_morphospace_pca(snum, which(arm), tw * sw, 0.95) else NULL
  if (!is.null(fr)) {
    pca_s <- fr$pca; k95 <- fr$k_95; sc <- sweep(snum[, fr$cols, drop = FALSE], 2L, pca_s$center, "-") %*% pca_s$rotation; sc[!vc, ] <- NA_real_; pca_s$x <- sc; smat <- smat[, fr$cols, drop = FALSE]
    for (p in seq_len(max(5L, min(max(as.integer(num_pcs), k95), 15L, ncol(sc))))) out_df[[paste0("PC", p)]] <- if (p <= ncol(sc)) ifelse(vc, round(sc[, p], 6), NA_real_) else ifelse(vc, 0.0, NA_real_)
  }
  list(df = out_df, sym_mat = smat, pca = pca_s, k_95 = as.integer(k95), permanova = NULL, active_ref_mask = arm)
}

# 3. Specimen Parsing & Contour Loading
parse_specimen_id <- function(filepath) {
  base <- sub("\\.[^.]+$", "", basename(filepath)); cat_num <- strsplit(base, "(_p|_leaf|_reflected|_curve)")[[1]][1]; p_id <- 0; l_id <- 1
  m_p <- regmatches(base, regexpr("_p([0-9]+)", base)); if (length(m_p) > 0) p_id <- as.integer(sub("_p", "", m_p))
  m_l <- regmatches(base, regexpr("leaf([0-9]+)|leaf_([0-9]+)|_([0-9]+)$", base)); if (length(m_l) > 0) l_id <- as.integer(gsub("[^0-9]", "", m_l))
  list(catalogNumber = cat_num, plant_id = p_id, leaf_id = l_id, shape_id = base)
}

load_contour_file <- function(filepath) {
  df <- tryCatch(read.csv(filepath, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(df) && nrow(df) >= 10 && all(c("x", "y") %in% names(df))) as.matrix(cbind(as.numeric(df$x), as.numeric(df$y))) else NULL
}

.ingest_step03_csv_contours <- function(opts, sec_dup_cats, active_voucher_cats) {
  mdf <- if (!is.null(opts$manifest) && file.exists(opts$manifest)) read.csv(opts$manifest, stringsAsFactors = FALSE) else NULL
  if (!is.null(mdf) && "catalogNumber" %in% names(mdf)) {
    if (length(sec_dup_cats) > 0) mdf <- mdf[!as.character(mdf$catalogNumber) %in% sec_dup_cats, , drop = FALSE]
    if (length(active_voucher_cats) > 0) mdf <- mdf[trimws(as.character(mdf$catalogNumber)) %in% active_voucher_cats, , drop = FALSE]
  }
  cfiles <- if (!is.null(mdf) && nrow(mdf) > 0 && "contour_path" %in% names(mdf)) { mdf <- mdf[!duplicated(mdf$contour_path, fromLast = TRUE), , drop = FALSE]; unique(mdf$contour_path[file.exists(mdf$contour_path)]) } else character(0)
  if (length(cfiles) == 0 && !is.null(opts$contours_dir) && dir.exists(opts$contours_dir)) cfiles <- list.files(opts$contours_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(cfiles) == 0) stop("No contour CSVs found in ", opts$contours_dir, " or manifest.")
  cfiles <- cfiles[!duplicated(sub("\\.[^.]+$", "", basename(cfiles)), fromLast = TRUE)]
  if (length(active_voucher_cats) > 0) cfiles <- cfiles[vapply(cfiles, function(f) trimws(as.character(parse_specimen_id(f)$catalogNumber)), character(1)) %in% active_voucher_cats]
  mb <- mp <- mc <- NULL
  if (!is.null(mdf) && nrow(mdf) > 0) {
    if ("contour_path" %in% names(mdf)) mb <- stats::setNames(seq_len(nrow(mdf)), basename(as.character(mdf$contour_path)))
    if (all(c("catalogNumber", "leaf_id") %in% names(mdf))) {
      pm <- if ("plant_individual_id" %in% names(mdf)) ifelse(is.na(mdf$plant_individual_id), 0L, as.integer(mdf$plant_individual_id)) else rep(0L, nrow(mdf))
      mp <- stats::setNames(seq_len(nrow(mdf)), paste0(mdf$catalogNumber, "__", pm, "__", as.integer(mdf$leaf_id)))
      mc <- stats::setNames(seq_len(nrow(mdf)), paste0(mdf$catalogNumber, "__", as.integer(mdf$leaf_id)))
    }
  }
  mc_all <- if (!is.null(mdf)) names(mdf) else character(0)
  dyn_cols <- unique(c("blade_cie_l_norm", "blade_cie_l_sd", discover_step03_photometry_cols(mc_all), "cy_rel_sheet", intersect(c("cx_px", "cy_px", "leaf_cx_px", "leaf_cy_px", "cx_rel_sheet", "cy_rel_sheet"), mc_all)))
  has_m <- requireNamespace("Momocs", quietly = TRUE); message("Ingesting contour coordinates from ", length(cfiles), " CSV files...")

  ingest_one <- function(f) {
    meta <- parse_specimen_id(f); mm <- NULL
    if (!is.null(mdf) && nrow(mdf) > 0) {
      ih <- if (!is.null(mb)) unname(mb[basename(f)]) else NA_integer_
      if (is.na(ih) && !is.null(mp)) ih <- unname(mp[paste0(meta$catalogNumber, "__", as.integer(meta$plant_id), "__", as.integer(meta$leaf_id))])
      if (is.na(ih) && !is.null(mc)) ih <- unname(mc[paste0(meta$catalogNumber, "__", as.integer(meta$leaf_id))])
      if (!is.na(ih) && ih >= 1L && ih <= nrow(mdf)) mm <- mdf[ih, , drop = FALSE]
    }
    if (!is.null(mm) && exists("evaluate_step03_lm3_efa_qc", mode = "function") && !isTRUE(evaluate_step03_lm3_efa_qc(mm, allow_multi_component_ect = TRUE)[1L])) return(NULL)
    coo <- load_contour_file(f); if (is.null(coo)) return(NULL)
    ms <- function(cn, fb = NA_character_) if (!is.null(mm) && cn %in% names(mm) && !is.na(mm[[cn]][1L])) as.character(mm[[cn]][1L]) else fb
    mn <- function(cn) if (!is.null(mm) && cn %in% names(mm) && !is.na(mm[[cn]][1L])) suppressWarnings(as.numeric(mm[[cn]][1L])) else NA_real_
    ml <- function(cn, fb = NA) if (!is.null(mm) && cn %in% names(mm) && !is.na(mm[[cn]][1L])) as.logical(mm[[cn]][1L]) else fb
    m_up <- ml("mask_is_upright", NA); was_tu <- FALSE
    if (exists("align_contour_horizontal_xax", mode = "function")) { coo <- align_contour_horizontal_xax(coo, is_upright = m_up); was_tu <- isTRUE(attr(coo, "lm3_from_tipup")) }
    if (has_m) coo <- tryCatch(Momocs::coo_cw(coo), error = function(e) coo); attr(coo, "lm3_xax_aligned") <- TRUE; attr(coo, "lm3_from_tipup") <- was_tu

    ppm <- mn("pixels_per_mm"); cr <- resolve_step03_calibration_provenance(ppm, if (!is.null(mm) && "is_calibrated_mm" %in% names(mm)) mm$is_calibrated_mm[1L] else NULL, ms("scale_source", ""))
    is_cal <- cr$is_calibrated_mm[1L]; mmn <- function(cn) if (is_cal) mn(cn) else NA_real_
    area_px <- abs(0.5 * sum(coo[, 1] * c(coo[-1, 2], coo[1, 2]) - coo[, 2] * c(coo[-1, 1], coo[1, 1])))
    perim_px <- sum(sqrt(diff(c(coo[, 1], coo[1, 1]))^2 + diff(c(coo[, 2], coo[1, 2]))^2))
    circ <- if (is.finite(perim_px) && perim_px > 1e-6) round(min(1.0, (4 * pi * area_px) / (perim_px^2)), 4) else NA_real_
    di <- compute_dissection_index_from_coo(area_val = area_px, perim_val = perim_px)
    sol <- if (!is.na(mn("solidity"))) round(mn("solidity"), 4) else compute_polygon_solidity(coo, area_px)
    ncc <- if (!is.na(mn("n_lamina_components"))) max(1L, as.integer(mn("n_lamina_components"))) else 1L
    if (ncc > 1L && !is.na(mn("dissection_index"))) di <- round(mn("dissection_index"), 4)
    cs_px <- compute_centsize_robust(coo); if (is_cal) coo <- coo / ppm; cs_mm <- if (is_cal) compute_centsize_robust(coo) else NA_real_
    coo <- if (has_m) tryCatch(Momocs::coo_center(coo), error = function(e) sweep(coo, 2, colMeans(coo, na.rm = TRUE), "-")) else sweep(coo, 2, colMeans(coo, na.rm = TRUE), "-")
    attr(coo, "lm3_xax_aligned") <- TRUE; attr(coo, "lm3_from_tipup") <- was_tu
    lw <- if (has_m) tryCatch(Momocs::coo_lw(coo), error = function(e) NULL) else NULL
    if (is.null(lw) || length(lw) < 2 || is.na(lw[2]) || lw[2] <= 1e-6) { cc <- scale(coo, TRUE, FALSE); sr <- cc %*% svd(cc)$v; lw <- c(diff(range(sr[, 1], na.rm = TRUE)), diff(range(sr[, 2], na.rm = TRUE))) }
    alw <- if (was_tu && diff(range(coo[, 2], na.rm = TRUE)) > 1e-6) diff(range(coo[, 1], na.rm = TRUE)) / diff(range(coo[, 2], na.rm = TRUE)) else NA_real_
    ar <- round(as.numeric(if (is.finite(mn("aspect_ratio")) && mn("aspect_ratio") > 0) mn("aspect_ratio") else if (is.finite(alw) && alw > 0) alw else (lw[1] / lw[2])), 4)
    rec <- data.frame(shape_id = meta$shape_id, catalogNumber = meta$catalogNumber,
      plant_individual_id = if (!is.na(mn("plant_individual_id"))) as.integer(mn("plant_individual_id")) else meta$plant_id,
      leaf_id = meta$leaf_id, n_lamina_components = ncc, enclosed_void_frac = if (!is.na(mn("enclosed_void_frac"))) mn("enclosed_void_frac") else 0.0,
      qc_flag_reason = ms("qc_flag_reason", "PASS"), mask_is_upright = if (!is.na(m_up)) m_up else was_tu, ect_path = ms("ect_path"), leaf_zone = ms("leaf_zone"),
      leaf_zone_confidence = mn("leaf_zone_confidence"), dist_from_crown_norm = mn("dist_from_crown_norm"),
      area_mm2 = if (!is_cal) NA_real_ else if (!is.na(mn("area_mm2"))) mn("area_mm2") else round(area_px / (ppm^2), 2),
      pixels_per_mm = ppm, is_calibrated_mm = is_cal, scale_source = cr$scale_source[1L], aspect_ratio = ar, circularity = circ, solidity = sol,
      dissection_index = di, archetype_score = mn("archetype_score"), petiole_length_mm = mmn("petiole_length_mm"), petiole_width_mm = mmn("petiole_width_mm"),
      blade_length_mm = mmn("blade_length_mm"), blade_width_mm = mmn("blade_width_mm"), petiole_to_blade_ratio = mn("petiole_to_blade_ratio"),
      is_dissected = ml("is_dissected", (ncc > 1L) || (is.finite(sol) && (sol < 0.72))), cy_rel_plant_bbox = mn("cy_rel_plant_bbox"), area_px = round(area_px, 1),
      centroid_size_px = if (!is.na(cs_px)) round(cs_px, 2) else NA_real_, centroid_size = if (is_cal && !is.na(cs_mm)) round(cs_mm, 4) else NA_real_,
      log_centsize = if (is_cal && !is.na(cs_mm) && cs_mm > 0) round(log10(cs_mm), 6) else NA_real_, mask_source = f, stringsAsFactors = FALSE)
    for (dc in dyn_cols) rec[[dc]] <- mn(dc)
    list(shape_name = meta$shape_id, coo = coo, rec = rec)
  }
  ing_l <- if (length(cfiles) >= 100 && exists("packera_mclapply", mode = "function")) packera_mclapply(cfiles, ingest_one, mc.cores = resolve_packera_workers(32L, length(cfiles))) else lapply(cfiles, ingest_one)
  valid <- Filter(Negate(is.null), ing_l); if (length(valid) == 0) stop("Failed to extract valid QC-passing coordinate contours.")
  list(fac_df = do.call(rbind, lapply(valid, `[[`, "rec")), coo_list = lapply(valid, `[[`, "coo"), shape_names = vapply(valid, `[[`, character(1), "shape_name"))
}

# 4. Standardized Workflow Orchestrator
run_fourier_extraction <- function(opts) {
  if (is.null(opts$contours_dir) || !nzchar(opts$contours_dir)) opts$contours_dir <- opts$input
  ex_cli <- attr(opts, "explicit_cli"); snap <- attr(opts, "parsed_snapshot")
  if (is.null(ex_cli)) ex_cli <- opts[intersect(c("harmonics", "num_pcs", "symmetrize_efa_harmonics", "seed", "random_seed", "min_export_pcs", "max_export_pcs", "min_harmonics", "max_harmonics", "target_harmonic_power", "target_pca_variance", "min_anchor_n_for_empirical", "cohort_adaptive_mode"), names(opts))]
  else if (!is.null(snap)) for (k in c("harmonics", "num_pcs", "symmetrize_efa_harmonics", "seed")) if (!is.null(opts[[k]]) && !identical(opts[[k]], snap[[k]])) ex_cli[[k]] <- opts[[k]]
  f_cfg <- load_step03_fourier_config(cfg_path = if (!is.null(opts$config)) opts$config else "config/config.yaml", cli_opts = ex_cli)
  message(sprintf("Starting Step 03 EFA | Contours: %s | H_init=%d | PCs=%d | Out: %s", opts$contours_dir, f_cfg$harmonics, f_cfg$num_pcs, opts$output))

  vouchers_df <- if (!is.null(opts$vouchers) && file.exists(opts$vouchers)) read.csv(opts$vouchers, stringsAsFactors = FALSE) else NULL
  sec_dup_cats <- character(0); active_voucher_cats <- character(0)
  if (!is.null(vouchers_df)) {
    if ("is_primary_duplicate" %in% names(vouchers_df)) { is_sec <- toupper(trimws(as.character(vouchers_df$is_primary_duplicate))) %in% c("FALSE", "F", "0"); sec_dup_cats <- unique(as.character(vouchers_df$catalogNumber[is_sec])); vouchers_df <- vouchers_df[!is_sec, , drop = FALSE] }
    if ("catalogNumber" %in% names(vouchers_df)) active_voucher_cats <- unique(trimws(as.character(vouchers_df$catalogNumber)))
    message("Loaded ", nrow(vouchers_df), " curated primary voucher records.")
  }

  harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(f_cfg$harmonics), function(x, y) paste0(x, y)))
  db_ing <- if (exists("fetch_step03_contours_from_db", mode = "function")) fetch_step03_contours_from_db(opts$vouchers, opts$manifest, opts$contours_dir, active_catalogs = if (length(active_voucher_cats) > 0L) setdiff(active_voucher_cats, sec_dup_cats) else NULL, align_fn = function(coo, ...) sweep(coo, 2, colMeans(coo, na.rm = TRUE), "-"), centsize_fn = compute_centsize_robust) else NULL
  if (!is.null(db_ing)) {
    coo_list <- db_ing$coords_list; fac_df <- db_ing$meta_df; shape_names_vec <- names(coo_list)
    cal_db <- resolve_step03_calibration_provenance(if ("pixels_per_mm" %in% names(fac_df)) fac_df$pixels_per_mm else NA_real_, if ("is_calibrated_mm" %in% names(fac_df)) fac_df$is_calibrated_mm else NULL, if ("scale_source" %in% names(fac_df)) fac_df$scale_source else NULL)
    fac_df$is_calibrated_mm <- cal_db$is_calibrated_mm; fac_df$scale_source <- cal_db$scale_source
    if (any(!cal_db$is_calibrated_mm)) for (gc in intersect(c("centroid_size", "log_centsize", "area_mm2", "petiole_length_mm", "petiole_width_mm", "blade_length_mm", "blade_width_mm"), names(fac_df))) fac_df[[gc]][!cal_db$is_calibrated_mm] <- NA_real_
    message("Ingested ", length(coo_list), " QC-passing 120-vertex contours directly from SQLite extracted_leaves.")
  } else { csv_ing <- .ingest_step03_csv_contours(opts, sec_dup_cats, active_voucher_cats); fac_df <- csv_ing$fac_df; coo_list <- csv_ing$coo_list; shape_names_vec <- csv_ing$shape_names }

  keep_u <- !duplicated(shape_names_vec, fromLast = TRUE) & !(as.character(fac_df$catalogNumber) %in% sec_dup_cats) & (length(active_voucher_cats) == 0L | trimws(as.character(fac_df$catalogNumber)) %in% active_voucher_cats)
  fac_df <- fac_df[keep_u, , drop = FALSE]; coo_list <- coo_list[keep_u]; names(coo_list) <- shape_names_vec[keep_u]; rownames(fac_df) <- names(coo_list)
  had_explicit_lz <- ("leaf_zone" %in% names(fac_df) && any(!is.na(fac_df$leaf_zone) & nzchar(trimws(as.character(fac_df$leaf_zone))))) || ("dist_from_crown_norm" %in% names(fac_df) && any(is.finite(suppressWarnings(as.numeric(fac_df$dist_from_crown_norm)))))
  if (exists("filter_basal_lamina_contours", mode = "function")) { bg <- filter_basal_lamina_contours(fac_df, coo_list, retain_cauline_tagged = TRUE); fac_df <- bg$fac_df; coo_list <- bg$coo_list }
  if (exists("align_contour_horizontal_xax", mode = "function")) coo_list <- Map(align_contour_horizontal_xax, coo_list, if ("mask_is_upright" %in% names(fac_df)) as.logical(fac_df$mask_is_upright) else rep(NA, length(coo_list)))
  message(sprintf("Assembled %d valid QC-passing contours for Elliptic Fourier Analysis & ECT.", length(coo_list)))
  if (exists("emit_r_progress", mode = "function")) emit_r_progress(3L, "contour_load", length(coo_list), length(coo_list), unit = "contours", substage_weight = 0.25, substage_offset = 0.0, status = "completed")

  # 3. Dynamic Harmonic Power Calibration & Homologized EFA
  eff_min_h <- if (isTRUE(f_cfg$cohort_adaptive_mode)) f_cfg$min_harmonics else f_cfg$harmonics; eff_max_h <- if (isTRUE(f_cfg$cohort_adaptive_mode)) f_cfg$max_harmonics else f_cfg$harmonics
  cal_fn <- function(o = NULL, harm_mat = NULL, nb_h = f_cfg$harmonics) verify_harmonic_power_calibration(o, harm_mat = harm_mat, nb_h = nb_h, min_h = eff_min_h, max_h = eff_max_h, target_power = f_cfg$target_harmonic_power)
  efa_ok <- FALSE; harm_mat <- out_obj <- ef_res <- NULL; h_opt <- max(eff_min_h, min(eff_max_h, as.integer(f_cfg$harmonics)))
  if (requireNamespace("Momocs", quietly = TRUE)) {
    if (length(coo_list) >= 240L && exists("run_parallel_momocs_efa", mode = "function") && exists("resolve_packera_workers", mode = "function") && resolve_packera_workers(32L) > 1L) {
      p_efa <- tryCatch(run_parallel_momocs_efa(coo_list, fac_df, f_cfg$harmonics, cal_fn), error = function(e) NULL)
      if (!is.null(p_efa$harm_mat)) { out_obj <- p_efa$out_obj; fac_df <- p_efa$fac_df; h_opt <- p_efa$h_opt; ef_res <- p_efa$ef_res; harm_mat <- p_efa$harm_mat; efa_ok <- TRUE; harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(h_opt), function(x, y) paste0(x, y))) }
    }
    if (!efa_ok && !is.null(out_obj <- tryCatch(Momocs::Out(coo_list, fac = fac_df), error = function(e) NULL))) {
      cs_v <- tryCatch(Momocs::coo_centsize(out_obj), error = function(e) NULL)
      if (!is.null(cs_v) && length(cs_v) == nrow(fac_df)) {
        cal_ok <- if ("is_calibrated_mm" %in% names(fac_df)) (!is.na(fac_df$is_calibrated_mm) & as.logical(fac_df$is_calibrated_mm)) else rep(TRUE, nrow(fac_df))
        fac_df$centroid_size <- ifelse(cal_ok, round(as.numeric(cs_v), 4), NA_real_)
        fac_df$log_centsize <- ifelse(cal_ok & !is.na(fac_df$centroid_size) & fac_df$centroid_size > 0, round(log10(fac_df$centroid_size), 6), NA_real_)
      }
      out_obj <- Momocs::coo_center(out_obj) %>% Momocs::coo_scale(); h_opt <- cal_fn(out_obj, nb_h = f_cfg$harmonics); harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(h_opt), function(x, y) paste0(x, y)))
      ef_res <- tryCatch(Momocs::efourier(out_obj, nb.h = h_opt, norm = FALSE, start = FALSE), error = function(e) NULL)
      if (!is.null(ef_res$coe)) { harm_mat <- as.matrix(ef_res$coe); efa_ok <- TRUE }
    }
  }
  if (!efa_ok || is.null(harm_mat)) {
    h_opt <- cal_fn(out_obj, nb_h = f_cfg$harmonics); harmonic_names <- as.vector(outer(c("A", "B", "C", "D"), seq_len(h_opt), function(x, y) paste0(x, y)))
    h_l <- if (length(coo_list) >= 240L && exists("packera_mclapply", mode = "function")) packera_mclapply(coo_list, function(c) compute_efourier_coords_pure(c, nb_h = h_opt, norm = FALSE, start = FALSE), mc.cores = resolve_packera_workers(32L, length(coo_list))) else lapply(coo_list, function(c) compute_efourier_coords_pure(c, nb_h = h_opt, norm = FALSE, start = FALSE))
    harm_mat <- do.call(rbind, h_l); colnames(harm_mat) <- harmonic_names
  }
  if (all(harmonic_names %in% colnames(harm_mat))) harm_mat <- harm_mat[, harmonic_names, drop = FALSE] else if (ncol(harm_mat) == length(harmonic_names)) colnames(harm_mat) <- harmonic_names
  if (exists("emit_r_progress", mode = "function")) emit_r_progress(3L, "efa_compute", nrow(harm_mat), nrow(harm_mat), unit = "contours", substage_weight = 0.50, substage_offset = 0.25, status = "completed")

  non_jordan_mask <- rep(FALSE, nrow(fac_df))
  if ("n_lamina_components" %in% names(fac_df)) non_jordan_mask <- non_jordan_mask | (!is.na(fac_df$n_lamina_components) & as.integer(fac_df$n_lamina_components) > 1L)
  if ("enclosed_void_frac" %in% names(fac_df)) { ev <- suppressWarnings(as.numeric(fac_df$enclosed_void_frac)); non_jordan_mask <- non_jordan_mask | (!is.na(ev) & is.finite(ev) & ev > 0.020) }
  if ("qc_flag_reason" %in% names(fac_df)) non_jordan_mask <- non_jordan_mask | (!is.na(fac_df$qc_flag_reason) & grepl("ENCLOSED_INTERLOBE_VOID|MULTI_COMPONENT_PINNAE", toupper(trimws(as.character(fac_df$qc_flag_reason)))))
  if (any(non_jordan_mask)) { harm_mat[which(non_jordan_mask), ] <- NA_real_; if (!is.null(ef_res$coe)) { ef_res$coe <- harm_mat; ef_res$fac <- fac_df } }

  if (!is.null(vouchers_df)) {
    mc <- intersect(colnames(vouchers_df), c("catalogNumber", "taxon", "species_standardized", "scientificName", "species_raw", "institutionCode", "identifiedBy", "determiner_tier", "authority_tier", "type_status", "dawid_skene_error_rate", "label_error_rate", "stateProvince", "county", "locality", "decimalLatitude", "decimalLongitude", "usda_symbol", "usda_county_status", "exsiccatae_key", "exsiccatae_taxonomic_conflict", "exsiccatae_gold_taxon", "coordinateUncertainty", "year", "month", "day", "eventDate", "regional_group", "spatial_balance_weight"))
    fac_df$.row_ord <- seq_len(nrow(fac_df)); fac_df <- merge(fac_df, vouchers_df[!duplicated(vouchers_df$catalogNumber), mc, drop = FALSE], by = "catalogNumber", all.x = TRUE, sort = FALSE)
    fac_df <- fac_df[order(fac_df$.row_ord), , drop = FALSE]; fac_df$.row_ord <- NULL; rownames(fac_df) <- fac_df$shape_id
  }
  if (!"species_raw" %in% names(fac_df)) fac_df$species_raw <- if ("taxon" %in% names(fac_df)) fac_df$taxon else NA_character_
  fac_df$scientificName <- if (!"scientificName" %in% names(fac_df)) fac_df$species_raw else ifelse(is.na(fac_df$scientificName) | fac_df$scientificName == "", fac_df$species_raw, fac_df$scientificName)
  if (exists("standardize_taxon_name", mode = "function")) {
    sp_s <- rep(NA_character_, nrow(fac_df)); for (sc in intersect(c("species_standardized", "scientificName", "species_raw", "taxon"), names(fac_df))) { sv <- standardize_taxon_name(fac_df[[sc]]); h <- is.na(sp_s) & !is.na(sv) & nzchar(sv); sp_s[h] <- sv[h] }
    fac_df$species_standardized <- sp_s; fac_df$scientificName <- ifelse(!is.na(sp_s) & nzchar(sp_s), sp_s, fac_df$scientificName)
  }
  fb_tier <- if ("authority_tier" %in% names(fac_df)) ifelse(is.na(fac_df$authority_tier) | fac_df$authority_tier == "", "Tier_3_Bronze", as.character(fac_df$authority_tier)) else rep("Tier_3_Bronze", nrow(fac_df))
  fac_df$determiner_tier <- if (!"determiner_tier" %in% names(fac_df)) fb_tier else ifelse(is.na(fac_df$determiner_tier) | fac_df$determiner_tier == "", fb_tier, fac_df$determiner_tier)
  fac_df$is_type_specimen <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(fac_df)$is_type_specimen else rep(FALSE, nrow(fac_df))
  if (exists("consolidate_fragmented_sheet_individuals", mode = "function")) fac_df <- consolidate_fragmented_sheet_individuals(fac_df)
  fac_df$authority_weight <- compute_authority_precision_weights(fac_df); fac_df$spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(fac_df) else rep(1.0, nrow(fac_df))

  # 4. Symmetric Fourier Decomposition (Leaf QC) & ECT Extraction
  sym_names <- c(paste0("A", seq_len(h_opt)), paste0("D", seq_len(h_opt)))
  ef_input <- if (!is.null(ef_res)) { ef_res$fac <- fac_df; ef_res } else list(coe = harm_mat, fac = fac_df)
  sym_res <- extract_symmetric_fourier_pcs(ef_input, fac_df = fac_df, num_pcs = f_cfg$num_pcs, seed = f_cfg$random_seed, fit_pca = FALSE)
  if (isTRUE(f_cfg$symmetrize_efa_harmonics)) { harm_mat <- sym_res$sym_mat; if (any(non_jordan_mask)) harm_mat[which(non_jordan_mask), ] <- NA_real_; harmonic_names <- colnames(harm_mat) }
  efa_df <- cbind(fac_df, as.data.frame(harm_mat))
  for (scol in c("sym_var", "asym_var", "pct_asymmetry", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight")) if (scol %in% names(sym_res$df)) efa_df[[scol]] <- sym_res$df[[scol]]
  efa_df$efd_harmonic_opt <- as.integer(h_opt)
  if (exists("compute_leaf_ect_features", mode = "function")) efa_df <- compute_leaf_ect_features(efa_df, coo_list, opts$contours_dir)
  ect_attrs <- attributes(efa_df)[c("ect_sym_mats", "ect_sym_vecs", "ect_coo_list", "ect_basal_var_exp", "ect_sym_dim")]

  # 5. Basal OTU Aggregation & Authoritative OTU-Level Weighted PCA
  valid_closed <- complete.cases(efa_df[, intersect(sym_names, names(efa_df)), drop = FALSE]) & !non_jordan_mask
  ind_keys <- paste(efa_df$catalogNumber, efa_df$plant_individual_id, sep = "__"); agg_keep <- valid_closed | !(ind_keys %in% unique(ind_keys[valid_closed]))
  build_clf <- function(edf) { clf <- edf[agg_keep, , drop = FALSE]; if (!is.null(ect_attrs$ect_sym_vecs)) attr(clf, "ect_sym_vecs") <- ect_attrs$ect_sym_vecs[agg_keep]; attr(clf, "ect_basal_var_exp") <- ect_attrs$ect_basal_var_exp; attr(clf, "ect_sym_dim") <- ect_attrs$ect_sym_dim; clf }
  closed_leaf_df <- build_clf(efa_df); all_harm_cols <- intersect(harmonic_names, names(efa_df))
  specimen_df <- aggregate_step03_otu_medians_fast(closed_leaf_df, all_harm_cols, sym_names)
  if (exists("apply_otu_confidence_routing_tier", mode = "function")) specimen_df <- apply_otu_confidence_routing_tier(specimen_df, closed_leaf_df)

  specimen_df$is_type_specimen <- if (exists(".identify_type_specimen_vec", mode = "function")) .identify_type_specimen_vec(specimen_df)$is_type_specimen else rep(FALSE, nrow(specimen_df))
  sym_cols_spec <- intersect(sym_names, names(specimen_df)); otu_has_efa <- if (length(sym_cols_spec) >= 1L) complete.cases(specimen_df[, sym_cols_spec, drop = FALSE]) else rep(FALSE, nrow(specimen_df))
  otu_base_mask <- otu_has_efa & !(specimen_df$morph_tier_source %in% c("Cauline_Only_Fallback", "Transitional_Passive")) & (!as.logical(specimen_df$is_contour_outlier) | as.logical(specimen_df$is_type_specimen))
  if (sum(otu_base_mask, na.rm = TRUE) < 2L) otu_base_mask <- otu_has_efa
  otu_h_cols <- sym_cols_spec
  if (sum(otu_base_mask, na.rm = TRUE) >= 2L && length(sym_cols_spec) >= 1L) { hv <- apply(as.matrix(specimen_df[otu_base_mask, sym_cols_spec, drop = FALSE]), 2L, stats::var, na.rm = TRUE); nz <- names(hv[!is.na(hv) & hv > 1e-10]); if (length(nz) >= 1L) otu_h_cols <- nz }
  otu_ref_mask <- if (sum(otu_base_mask, na.rm = TRUE) >= 2L && length(otu_h_cols) >= 1L) prune_within_taxon_mahalanobis(as.matrix(specimen_df[, otu_h_cols, drop = FALSE]), specimen_df, otu_base_mask, q = 0.95, min_anchors_per_taxon = f_cfg$min_anchor_n_for_empirical, taxon_col = "species_standardized") else otu_base_mask
  specimen_df$is_promoted_bronze <- if (!is.null(attr(otu_ref_mask, "is_promoted_bronze"))) as.logical(attr(otu_ref_mask, "is_promoted_bronze")) else rep(FALSE, nrow(specimen_df))
  specimen_df$anchor_tier_source <- if (!is.null(attr(otu_ref_mask, "anchor_tier_source"))) as.character(attr(otu_ref_mask, "anchor_tier_source")) else ifelse(specimen_df$is_type_specimen, "Tier_0_Type_Anchor", "Passive_Unpromoted")
  specimen_df$authority_weight <- compute_authority_precision_weights(specimen_df, active_mask = otu_ref_mask)
  nb_p <- !is.na(specimen_df$morph_tier_source) & (specimen_df$morph_tier_source %in% c("Cauline_Only_Fallback", "Transitional_Passive")) & !as.logical(specimen_df$is_type_specimen)
  if (any(nb_p)) specimen_df$authority_weight[nb_p] <- 0.0
  specimen_df$is_passive <- !(is.finite(specimen_df$authority_weight) & specimen_df$authority_weight > 0.0); if (sum(otu_ref_mask, na.rm = TRUE) < 2L) otu_ref_mask <- otu_base_mask
  specimen_df$spatial_balance_weight <- if (exists("compute_spatial_balance_weights", mode = "function")) round(compute_spatial_balance_weights(specimen_df, active_mask = otu_ref_mask, recompute = TRUE), 6) else rep(1.0, nrow(specimen_df))

  spec_k <- paste(specimen_df$catalogNumber, specimen_df$plant_individual_id, sep = "__"); m_otu <- match(ind_keys, spec_k)
  for (pcol in c("is_promoted_bronze", "anchor_tier_source", "authority_weight", "is_passive", "spatial_balance_weight")) {
    vp <- specimen_df[[pcol]][m_otu]; if (is.logical(vp)) vp[is.na(vp)] <- (pcol == "is_passive") else if (is.numeric(vp)) vp[!is.finite(vp)] <- if (pcol == "spatial_balance_weight") 1.0 else 0.0 else vp[is.na(vp)] <- "Passive_Unpromoted"; efa_df[[pcol]] <- vp
  }
  lz_l <- if ("leaf_zone" %in% names(efa_df)) tolower(trimws(as.character(efa_df$leaf_zone))) else rep("basal", nrow(efa_df))
  ol <- if ("is_contour_outlier" %in% names(efa_df)) (!is.na(efa_df$is_contour_outlier) & as.logical(efa_df$is_contour_outlier)) else rep(FALSE, nrow(efa_df))
  tl <- if ("is_type_specimen" %in% names(efa_df)) (!is.na(efa_df$is_type_specimen) & as.logical(efa_df$is_type_specimen)) else rep(FALSE, nrow(efa_df))
  fp <- ((!is.na(lz_l) & lz_l != "basal") | (ol & !tl)); if (any(fp)) { efa_df$authority_weight[fp] <- 0.0; efa_df$is_passive[fp] <- TRUE }

  pca_fit <- NULL; otu_cols <- otu_h_cols; otu_t1 <- which(otu_ref_mask & otu_has_efa); if (length(otu_t1) < 2L) otu_t1 <- which(otu_has_efa)
  k_95 <- as.integer(f_cfg$num_pcs); k_export <- max(as.integer(f_cfg$min_export_pcs), min(as.integer(f_cfg$max_export_pcs), as.integer(f_cfg$num_pcs)))
  fit_otu <- if (length(otu_t1) >= 2L && length(otu_h_cols) >= 1L) .fit_weighted_morphospace_pca(specimen_df[, otu_h_cols, drop = FALSE], otu_t1, specimen_df$authority_weight * specimen_df$spatial_balance_weight, f_cfg$target_pca_variance) else NULL
  if (!is.null(fit_otu)) {
    pca_fit <- fit_otu$pca; otu_cols <- fit_otu$cols; k_95 <- fit_otu$k_95
    k_export <- max(as.integer(f_cfg$min_export_pcs), min(as.integer(max(k_95, f_cfg$num_pcs)), as.integer(f_cfg$max_export_pcs), ncol(pca_fit$rotation)))
    sp_pcs <- sweep(as.matrix(specimen_df[, otu_cols, drop = FALSE]), 2L, pca_fit$center, "-") %*% pca_fit$rotation; sp_pcs[!otu_has_efa, ] <- NA_real_; pca_fit$x <- sp_pcs
    message(sprintf("=== Dynamic %d-Harmonic OTU Morphospace PCA (k_95 = %d, k_export = %d) ===", h_opt, k_95, k_export))
  }
  proj_pcs <- function(d, vm) {
    for (p in seq_len(k_export)) d[[paste0("PC", p)]] <- NA_real_; vi <- which(vm)
    if (length(vi) > 0L) { sm <- if (!is.null(pca_fit)) sweep(as.matrix(d[vi, otu_cols, drop = FALSE]), 2L, pca_fit$center, "-") %*% pca_fit$rotation else matrix(0.0, length(vi), k_export); for (p in seq_len(k_export)) d[vi, paste0("PC", p)] <- if (p <= ncol(sm)) round(sm[, p], 6) else 0.0 }
    d
  }
  c_idx <- which(complete.cases(as.matrix(efa_df[, intersect(otu_cols, names(efa_df)), drop = FALSE])) & !non_jordan_mask)
  specimen_df <- proj_pcs(specimen_df, otu_has_efa); efa_df <- proj_pcs(efa_df, seq_len(nrow(efa_df)) %in% c_idx)

  if (length(c_idx) >= 15L && !is.null(pca_fit)) {
    lpcs <- as.matrix(efa_df[c_idx, paste0("PC", seq_len(min(k_export, ncol(pca_fit$rotation)))), drop = FALSE])
    fpl <- which(efa_df$authority_weight[c_idx] > 0 & (is.na(lz_l[c_idx]) | lz_l[c_idx] == "basal") & !efa_df$is_contour_outlier[c_idx])
    kmcd <- max(2L, min(k_95, ncol(lpcs), max(2L, floor(nrow(lpcs) / 3L))))
    mo <- screen_regime_stratified_mcd_outliers(lpcs, fit_pos = fpl, solidity_vec = efa_df$solidity[c_idx], is_dissected_vec = efa_df$is_dissected[c_idx], k_mcd = kmcd, p_cutoff = 0.001, seed = f_cfg$random_seed)
    if ("is_type_specimen" %in% names(efa_df)) mo[as.logical(efa_df$is_type_specimen[c_idx])] <- FALSE
    new_f <- mo & !as.logical(efa_df$is_contour_outlier[c_idx])
    if (any(new_f) && sum(!mo) >= 5L) {
      efa_df$is_contour_outlier[c_idx] <- as.logical(efa_df$is_contour_outlier[c_idx]) | mo; efa_df$authority_weight[c_idx[new_f]] <- 0.0; efa_df$is_passive[c_idx[new_f]] <- TRUE
      message(sprintf("Regime-stratified MCD robust Mahalanobis outlier filter (top %d PCs, k_95=%d, p < 0.001): flagged %d leaves.", kmcd, k_95, sum(new_f)))
      prev_s <- specimen_df; closed_leaf_df <- build_clf(efa_df); specimen_df <- aggregate_step03_otu_medians_fast(closed_leaf_df, all_harm_cols, sym_names)
      if (exists("apply_otu_confidence_routing_tier", mode = "function")) specimen_df <- apply_otu_confidence_routing_tier(specimen_df, closed_leaf_df)
      mp <- match(paste(specimen_df$catalogNumber, specimen_df$plant_individual_id, sep = "__"), spec_k)
      for (pc in c("is_type_specimen", "is_promoted_bronze", "anchor_tier_source", "authority_weight", "is_passive", "spatial_balance_weight")) if (pc %in% names(prev_s)) specimen_df[[pc]] <- prev_s[[pc]][mp]
      nb2 <- !is.na(specimen_df$morph_tier_source) & (specimen_df$morph_tier_source %in% c("Cauline_Only_Fallback", "Transitional_Passive")) & !as.logical(specimen_df$is_type_specimen)
      if (any(nb2)) { specimen_df$authority_weight[nb2] <- 0.0; specimen_df$is_passive[nb2] <- TRUE }
      specimen_df <- proj_pcs(specimen_df, if (length(sym_cols_spec) >= 1L) complete.cases(specimen_df[, sym_cols_spec, drop = FALSE]) else rep(FALSE, nrow(specimen_df)))
    }
  }
  efa_df$efd_harmonic_opt <- specimen_df$efd_harmonic_opt <- as.integer(h_opt); efa_df$efd_k95_pcs <- specimen_df$efd_k95_pcs <- as.integer(k_95)
  if (exists("compute_separate_cauline_otu_track", mode = "function")) specimen_df <- compute_separate_cauline_otu_track(efa_df, specimen_df, sym_names)

  # 6. Export CSVs, SQLite Persistence & Step 03 Figures
  exp_pcs <- paste0("PC", seq_len(k_export))
  lead_cols <- intersect(c("catalogNumber", "plant_individual_id", "leaf_id", "mask_is_upright", "leaf_zone", "leaf_zone_confidence", "dist_from_crown_norm", "cy_rel_plant_bbox", "pixels_per_mm", "is_calibrated_mm", "scale_source", "species_standardized", "scientificName", "species_raw", "identifiedBy", "determiner_tier", "type_status", "is_type_specimen", "is_promoted_bronze", "anchor_tier_source", "authority_weight", "is_passive", "spatial_balance_weight", "efd_harmonic_opt", "efd_k95_pcs", exp_pcs, "ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity", "ect_bilateral_asym", "aspect_ratio", "circularity", "solidity", "dissection_index", "archetype_score", "petiole_length_mm", "petiole_width_mm", "blade_length_mm", "blade_width_mm", "petiole_to_blade_ratio", discover_step03_photometry_cols(names(efa_df)), "is_dissected", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_mm2", "area_px", "centroid_size_px", "centroid_size", "log_centsize", "mask_source"), names(efa_df))
  efa_df <- efa_df[, c(lead_cols, setdiff(names(efa_df), lead_cols))]; for (anm in intersect(c("ect_sym_mats", "ect_sym_vecs", "ect_coo_list", "ect_basal_var_exp", "ect_sym_dim"), names(ect_attrs))) attr(efa_df, anm) <- ect_attrs[[anm]]
  spec_attrs <- attributes(specimen_df)[c("ect_basal_var_exp", "ect_cauline_var_exp")]
  lead_spec <- intersect(c("catalogNumber", "plant_individual_id", "leaf_count", "basal_leaf_count", "cauline_leaf_count", "morph_tier_source", "foliar_variance", "within_rosette_efa_dispersion", "heteroblasty_efa_distance", "heteroblasty_ect_distance", "pixels_per_mm", "is_calibrated_mm", "scale_source", "mask_is_upright", "species_standardized", "scientificName", "species_raw", "identifiedBy", "determiner_tier", "type_status", "is_type_specimen", "is_promoted_bronze", "anchor_tier_source", "authority_weight", "is_passive", "spatial_balance_weight", "efd_harmonic_opt", "efd_k95_pcs", exp_pcs, "Cauline_PC1", "Cauline_PC2", "Cauline_PC3", "ECT_PC1", "ECT_PC2", "ECT_PC3", "Cauline_ECT_PC1", "Cauline_ECT_PC2", "ect_lobe_complexity", "cauline_ect_lobe_complexity", "aspect_ratio", "circularity", "solidity", "q25_solidity", "dissection_index", "dissected_leaf_fraction", "cauline_solidity_med", "cauline_dissection_index_med", "cauline_aspect_ratio_med", "cauline_blade_length_mm", "basal_vs_cauline_length_ratio", "basal_vs_cauline_solidity_delta", "archetype_score", "petiole_length_mm", "petiole_width_mm", "blade_length_mm", "blade_width_mm", "petiole_to_blade_ratio", discover_step03_photometry_cols(names(specimen_df)), "is_dissected", "harmonic_roughness", "is_contour_outlier", "leaf_quality_weight", "area_mm2", "area_px", "centroid_size", "log_centsize"), names(specimen_df))
  specimen_df <- specimen_df[, c(lead_spec, setdiff(names(specimen_df), lead_spec))]; for (anm in intersect(c("ect_basal_var_exp", "ect_cauline_var_exp"), names(spec_attrs))) attr(specimen_df, anm) <- spec_attrs[[anm]]

  if (!is.null(opts$output_individual) && nzchar(opts$output_individual)) {
    ind_df <- if (isTRUE(had_explicit_lz) || !("leaf_zone" %in% names(efa_df)) || !any(efa_df$leaf_zone == "basal", na.rm = TRUE)) efa_df else efa_df[is.na(efa_df$leaf_zone) | efa_df$leaf_zone == "basal", , drop = FALSE]
    write_step03_csv_fast(ind_df, opts$output_individual); message("Individual leaf EFA harmonics archive exported: ", opts$output_individual, " (Rows: ", nrow(ind_df), ")")
  }
  write_step03_csv_fast(specimen_df, opts$output); message("Specimen-aggregated EFA harmonics table exported: ", opts$output, " (Rows: ", nrow(specimen_df), ")")
  if (exists("persist_step03_to_sqlite", mode = "function")) tryCatch(persist_step03_to_sqlite(efa_df, specimen_df, opts$vouchers), error = function(e) warning("SQLite Step 03 persistence note: ", e$message))
  if (exists("emit_r_progress", mode = "function")) emit_r_progress(3L, "export", nrow(specimen_df), nrow(specimen_df), unit = "otus", substage_weight = 0.25, substage_offset = 0.75, status = "completed")
  tryCatch({ if (!exists("generate_step03_pub_figures", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R"); if (exists("generate_step03_pub_figures", mode = "function")) generate_step03_pub_figures(efa_df, specimen_df, harm_mat, pca_fit, opts$figures_dir, opts$output, pca_fit) }, error = function(e) warning("Step 03 publication figures warning: ", e$message))
  invisible(efa_df)
}

if (sys.nframe() == 0) run_fourier_extraction(parse_args_robust())
