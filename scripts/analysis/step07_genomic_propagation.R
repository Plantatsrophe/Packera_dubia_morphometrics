# Step 07 Genomic Propagation, Tier_0 Calibration, Deme Priors & gdi Governance
.safe_bool_prop <- function(x) if (is.logical(x)) !is.na(x) & x else if (is.numeric(x) || is.integer(x)) !is.na(x) & (x != 0) else tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
.pick_chr_prop <- function(df, cols, default = "") { for (c in cols) if (c %in% names(df)) { v <- as.character(df[[c]]); v[is.na(v)] <- default; return(v) }; rep(default, nrow(df)) }
.pick_num_prop <- function(df, cols) { for (c in cols) if (c %in% names(df)) return(suppressWarnings(as.numeric(df[[c]]))); rep(NA_real_, nrow(df)) }
.pick_bool_prop <- function(df, cols, default = FALSE) { for (c in cols) if (c %in% names(df)) return(.safe_bool_prop(df[[c]])); rep(default, nrow(df)) }

.is_type_prop <- function(df) {
  if (nrow(df) == 0L) return(logical(0))
  ts <- .pick_chr_prop(df, c("typeStatus", "type_status"), "")
  .pick_bool_prop(df, "is_type_specimen", FALSE) | (.pick_chr_prop(df, c("authority_tier", "determiner_tier"), "") == "Tier_0_Type_Anchor") |
    (grepl("type", ts, ignore.case = TRUE) & !grepl("non[- ]?type|not a type", ts, ignore.case = TRUE))
}

.init_propagation_cols <- function(df, reset = FALSE) {
  n <- nrow(df)
  for (c in c("genomic_propagation_tier", "tier0_calibrated_taxon")) df[[c]] <- if (reset || !(c %in% names(df))) rep(NA_character_, n) else as.character(df[[c]])
  for (c in c("deme_mahal_d2", "tier0_calibrated_prob", "cluster_propagated_gdi")) df[[c]] <- if (reset || !(c %in% names(df))) rep(NA_real_, n) else suppressWarnings(as.numeric(df[[c]]))
  df[["tier0_abstention_flag"]] <- if (reset || !("tier0_abstention_flag" %in% names(df))) rep(NA_integer_, n) else suppressWarnings(as.integer(df[["tier0_abstention_flag"]]))
  df
}

.resolve_prop_csv <- function(df, opts, attr_name, opt_names, filename) {
  for (on in opt_names) if (!is.null(opts[[on]]) && is.data.frame(opts[[on]])) return(opts[[on]])
  at <- attr(df, attr_name); if (!is.null(at) && is.data.frame(at)) return(at)
  od <- if (is.character(opts$output_dir) && length(opts$output_dir) == 1L && nzchar(opts$output_dir)) opts$output_dir else ""
  dirs <- unique(Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x), c(opts$genomics_dir, opts$genomic_dir, opts$tables_dir, if (nzchar(od)) file.path(od, c("tables", "reports")) else NULL)))
  for (d in dirs) if (file.exists(p <- file.path(d, filename))) {
    res <- tryCatch(utils::read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(res) && is.data.frame(res)) return(res)
  }
  NULL
}

.extract_collector_surname <- function(s) {
  vapply(as.character(s), function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return("")
    x <- trimws(strsplit(x, "[;|&]|\\bet\\b|\\band\\b", perl = TRUE)[[1L]][1L])
    if (grepl(",", x, fixed = TRUE)) part <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]][1L])
    else { toks <- unlist(strsplit(x, "\\s+")); toks <- toks[!grepl("^[A-Za-z]\\.?$|^[A-Z]{1,3}\\.$", toks)]; part <- if (length(toks) > 0L) toks[length(toks)] else "" }
    tolower(gsub("[^a-zA-Z0-9]+", "", part))
  }, character(1L), USE.NAMES = FALSE)
}

.norm_record_number <- function(x) {
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(""); raw <- tolower(trimws(s)); if (!nzchar(raw) || grepl("^0+$", raw)) return("")
    cln <- gsub("[^a-z0-9]+", "", raw); if (!nzchar(cln) || cln %in% c("sn", "unnumbered", "none", "na", "null", "unknown", "missing")) "" else cln
  }, character(1L), USE.NAMES = FALSE)
}

.haversine_km_vec <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2.0 * 6371.0088 * asin(pmin(1.0, sqrt(pmax(0.0, a))))
}

.is_basal_mask_prop <- function(df) {
  if ("morph_tier_source" %in% names(df) && any(!is.na(df$morph_tier_source) & nzchar(as.character(df$morph_tier_source)))) return(!is.na(df$morph_tier_source) & as.character(df$morph_tier_source) == "Basal_Rosette")
  if ("basal_leaf_count" %in% names(df) && any(!is.na(blc <- suppressWarnings(as.integer(df$basal_leaf_count))))) return(!is.na(blc) & blc >= 1L)
  rep(TRUE, nrow(df))
}

.compute_morph_pool_prop <- function(df, pc_cols = NULL, ridge = 1e-4) {
  if (is.null(pc_cols)) pc_cols <- intersect(paste0("PC", 1:6), names(df)); if (length(pc_cols) < 2L || nrow(df) == 0L) return(NULL)
  X <- matrix(NA_real_, nrow = nrow(df), ncol = length(pc_cols), dimnames = list(NULL, pc_cols))
  for (k in seq_along(pc_cols)) X[, k] <- suppressWarnings(as.numeric(df[[pc_cols[k]]]))
  is_basal <- .is_basal_mask_prop(df); cc <- stats::complete.cases(X); idx <- which(is_basal & cc); if (length(idx) < max(3L, length(pc_cols) + 1L)) idx <- which(cc)
  p <- ncol(X); if (length(idx) < 2L) return(list(X = X, S_inv = diag(1.0, p), pc_cols = pc_cols, is_basal = is_basal))
  S <- stats::cov(X[idx, , drop = FALSE]); S[!is.finite(S)] <- 0.0; S <- 0.5 * (S + t(S))
  S_reg <- 0.90 * S + 0.10 * diag(pmax(diag(S), 1e-3), nrow = p) + diag(ridge, p)
  list(X = X, S_inv = tryCatch(solve(S_reg), error = function(e) diag(1.0 / pmax(diag(S_reg), 1e-3), nrow = p)), pc_cols = pc_cols, is_basal = is_basal)
}

# Step A: Propagate G1_Direct_Voucher across duplicate exsiccatae sheets
propagate_exsiccatae_isoduplicates <- function(df, ga_df = NULL) {
  df <- .init_propagation_cols(df, reset = FALSE); n <- nrow(df); if (n < 2L) return(df)
  tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), ""); qc <- .pick_bool_prop(df, "geno_quality_pass", FALSE)
  adm <- .pick_bool_prop(df, "is_geno_admixed", FALSE); pred <- trimws(.pick_chr_prop(df, "geno_predicted_taxon", "")); prob <- .pick_num_prop(df, "geno_posterior_prob")

  if (!is.null(ga_df) && is.data.frame(ga_df) && nrow(ga_df) > 0L && "catalogNumber" %in% names(df) && "catalogNumber" %in% names(ga_df)) {
    ga_keep <- (.pick_chr_prop(ga_df, c("geno_linkage_tier", "linkage_tier"), "") == "G1_Direct_Voucher") & .pick_bool_prop(ga_df, c("geno_quality_pass", "qc_pass"), FALSE)
    if (any(ga_keep)) {
      ga_sub <- ga_df[ga_keep, , drop = FALSE]; m_idx <- match(as.character(df$catalogNumber), as.character(ga_sub$catalogNumber))
      fill_idx <- which(!nzchar(tier) & !is.na(m_idx))
      if (length(fill_idx) > 0L) for (col in intersect(c("geno_predicted_taxon", "geno_posterior_prob", "geno_quality_pass", "is_geno_admixed", "genetic_sample_id"), names(ga_sub))) if (!(col %in% names(df))) df[[col]] <- rep(NA, n)
      for (r in fill_idx) {
        g <- m_idx[r]; tier[r] <- "G1_Direct_Voucher"; qc[r] <- TRUE
        adm[r] <- .pick_bool_prop(ga_sub[g, , drop = FALSE], "is_geno_admixed", FALSE)[1L]
        pred[r] <- trimws(.pick_chr_prop(ga_sub[g, , drop = FALSE], "geno_predicted_taxon", "")[1L])
        prob[r] <- .pick_num_prop(ga_sub[g, , drop = FALSE], "geno_posterior_prob")[1L]
        for (col in setdiff(intersect(names(ga_sub), names(df)), "catalogNumber")) {
          if (is.na(df[[col]][r]) || (is.character(df[[col]][r]) && !nzchar(trimws(df[[col]][r])))) df[[col]][r] <- ga_sub[[col]][g]
        }
      }
    }
  }

  ex_k <- trimws(.pick_chr_prop(df, "exsiccatae_key", "")); sur_c <- trimws(.pick_chr_prop(df, "collector_surname", ""))
  sur <- .extract_collector_surname(ifelse(nzchar(sur_c), sur_c, .pick_chr_prop(df, "recordedBy", "")))
  recn <- .norm_record_number(.pick_chr_prop(df, "recordNumber", ""))
  evdt <- tolower(substr(trimws(.pick_chr_prop(df, "eventDate", "")), 1L, 10L))
  gkey <- ifelse(nzchar(ex_k), ex_k, ifelse(nzchar(sur) & nzchar(recn), paste(sur, recn, evdt, sep = "|"), ""))
  nonempty_k <- gkey[nzchar(gkey)]; if (length(nonempty_k) < 2L) return(df)
  ktab <- table(nonempty_k); is_donor <- (tier == "G1_Direct_Voucher") & qc & !adm & nzchar(pred) & (pred != "Unknown") & nzchar(gkey)
  if (!any(is_donor)) return(df)

  is_type <- .is_type_prop(df); cat_vec <- if ("catalogNumber" %in% names(df)) as.character(df$catalogNumber) else as.character(seq_len(n))
  n_pl <- if ("n_plant_individuals" %in% names(df)) suppressWarnings(as.integer(df$n_plant_individuals)) else rep(1L, n)
  is_multi <- (as.integer(table(cat_vec)[cat_vec]) > 1L) | (!is.na(n_pl) & n_pl > 1L)
  is_mo <- .pick_bool_prop(df, c("morphometric_outlier_flag", "is_morph_outlier"), FALSE); pool <- .compute_morph_pool_prop(df)
  copy_cols <- intersect(c(
    "genetic_sample_id", "geno_predicted_taxon", "geno_posterior_prob", "geno_second_taxon", "geno_margin", "geno_admixture_qmax", "q_max", "q_admixture_json",
    "genetic_cluster", "geno_sNMF_entropy", "is_geno_admixed", "is_putative_novel_or_hybrid_query", "geno_quality_pass", "geno_parent_1", "geno_parent_2",
    "nearest_tree_sister_taxon", "patristic_dist_to_dubia", "replicate_pair_concordance", "hdplot_retained_loci", "inferred_cytotype", "nquire_delta_logl",
    "plastid_clade", "its_ribotype", "cytonuclear_discordant", "hybrid_index_h", "interclass_het_Q12", "interspecific_ancestry_q12", "hybrid_generation_class",
    "coancestry_cluster_id", "coancestry_cluster", "fbranch_introgression_flag", "f_branch_max", "f_branch_donor_taxon", "snaq_gamma_minor", "coalescent_gdi",
    "bpp_introgression_bf", "mmrr_beta_ibd", "mmrr_beta_ibe", "mmrr_p_ibe", "bisse_net_div_rate", "public_backbone_clade", "public_phytop_ils_index",
    "public_phytop_ih_index", "public_plastome_match", "public_its_match", "waxy_lineage_match"
  ), names(df))

  for (k in unique(gkey[is_donor])) {
    if (is.na(ktab[k]) || ktab[k] < 2L) next
    d_idx <- which(is_donor & gkey == k); if (length(unique(pred[d_idx])) != 1L) next
    d <- d_idx[which.max(ifelse(is.finite(prob[d_idx]), prob[d_idx], -Inf))]
    for (j in setdiff(which(gkey == k), d_idx)) {
      if ((tier[j] == "G1_Direct_Voucher") || is_type[j] || (is_multi[j] && is_mo[j])) next
      if (!is.null(pool) && length(pool$pc_cols) >= 2L && pool$is_basal[d] && pool$is_basal[j] && all(is.finite(pool$X[d, ])) && all(is.finite(pool$X[j, ]))) {
        dv <- pool$X[j, ] - pool$X[d, ]
        if (is.finite(d2 <- as.numeric(crossprod(dv, pool$S_inv %*% dv))) && d2 > stats::qchisq(0.99, df = length(pool$pc_cols))) next
      }
      for (col in copy_cols) df[[col]][j] <- df[[col]][d]
      if (!("geno_linkage_tier" %in% names(df))) df$geno_linkage_tier <- rep(NA_character_, n)
      df$geno_linkage_tier[j] <- "G1_Exsiccatae_Isoduplicate"
      if ("linkage_tier" %in% names(df)) df$linkage_tier[j] <- "G1_Exsiccatae_Isoduplicate"
      df$genomic_propagation_tier[j] <- "G1_Exsiccatae_Isoduplicate"
    }
  }
  df
}

# Step B: Resolve Tier_0 DNA-verified anchors and anchor evictions
resolve_tier0_anchor_overrides <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(list(is_tier0_dna_anchor = logical(0), evict_from_given_anchor = logical(0), rectified_taxon = character(0), tier0_weight = numeric(0)))
  is_type <- .is_type_prop(df); tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), "")
  qc <- .pick_bool_prop(df, "geno_quality_pass", FALSE); adm <- .pick_bool_prop(df, "is_geno_admixed", FALSE)
  nov <- .pick_bool_prop(df, "is_putative_novel_or_hybrid_query", FALSE); pred <- trimws(.pick_chr_prop(df, "geno_predicted_taxon", ""))
  sp <- trimws(.pick_chr_prop(df, "species_standardized", "")); prob <- .pick_num_prop(df, "geno_posterior_prob")
  margin <- if ("geno_margin" %in% names(df)) suppressWarnings(as.numeric(df$geno_margin)) else rep(NA_real_, n)
  q12 <- .pick_num_prop(df, c("interclass_het_Q12", "interspecific_ancestry_q12"))
  is_t0 <- !is_type & (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc & !adm & !nov &
    is.finite(prob) & (prob >= 0.85) & (is.na(margin) | margin >= 0.40) & nzchar(pred) & (pred != "Unknown")
  evict <- !is_type & (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) & qc &
    ((nzchar(pred) & (pred != "Unknown") & (pred != sp) & is.finite(prob) & (prob >= 0.85)) | (is.finite(q12) & (q12 >= 0.35)))
  list(is_tier0_dna_anchor = is_t0, evict_from_given_anchor = evict, rectified_taxon = ifelse(is_t0, pred, sp), tier0_weight = ifelse(is_t0, 4.0, NA_real_))
}

get_tier0_rectified_species_and_evictions <- function(df) resolve_tier0_anchor_overrides(df)

.build_taxon_centroids <- function(df, morph_X, is_basal, t0) {
  is_type <- .is_type_prop(df); det <- .pick_chr_prop(df, c("authority_tier", "determiner_tier"), "Tier_3_Bronze")
  wt <- ifelse(t0$is_tier0_dna_anchor | is_type, 4.0, ifelse(grepl("tier_?1|gold", det, ignore.case = TRUE), 2.0,
               ifelse(grepl("tier_?2|silver", det, ignore.case = TRUE), 0.75, ifelse(.pick_bool_prop(df, "is_promoted_bronze", FALSE), 0.50, 0.25))))
  rect_sp <- trimws(t0$rectified_taxon); cc <- stats::complete.cases(morph_X) & is_basal; pure_ev <- t0$evict_from_given_anchor & !t0$is_tier0_dna_anchor
  ok_base <- cc & !pure_ev & !.pick_bool_prop(df, "is_confident_label_noise", FALSE) &
    (!.pick_bool_prop(df, "is_passive", FALSE) | is_type | t0$is_tier0_dna_anchor) & (!.pick_bool_prop(df, c("morphometric_outlier_flag", "is_morph_outlier"), FALSE) | is_type | t0$is_tier0_dna_anchor)
  mu <- list()
  for (tx in sort(unique(rect_sp[nzchar(rect_sp) & rect_sp != "Unknown"]))) {
    idx <- which(ok_base & rect_sp == tx); if (length(idx) == 0L) idx <- which(cc & rect_sp == tx & !pure_ev)
    if (length(idx) > 0L) mu[[tx]] <- as.numeric(colSums(morph_X[idx, , drop = FALSE] * wt[idx]) / sum(wt[idx]))
  }
  mu
}

# Step C: Conditional spatial deme priors (4-gate transfer)
propagate_conditional_deme_priors <- function(df, ga_df = NULL, morph_X = NULL, morph_S_inv = NULL) {
  df <- .init_propagation_cols(df, reset = FALSE); n <- nrow(df); if (n == 0L) return(df)
  extract_anc <- function(tbl, allow_g3, track_row) {
    tr <- .pick_chr_prop(tbl, c("geno_linkage_tier", "linkage_tier"), ""); qcf <- .pick_bool_prop(tbl, c("geno_quality_pass", "qc_pass"), FALSE)
    la <- .pick_num_prop(tbl, "decimalLatitude"); lo <- .pick_num_prop(tbl, "decimalLongitude"); ad <- .pick_bool_prop(tbl, "is_geno_admixed", FALSE)
    pr <- .pick_num_prop(tbl, "geno_posterior_prob"); hg <- trimws(.pick_chr_prop(tbl, "hybrid_generation_class", ""))
    pd <- trimws(.pick_chr_prop(tbl, "geno_predicted_taxon", "")); spv <- trimws(.pick_chr_prop(tbl, "species_standardized", ""))
    ok <- (tr %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy", if (allow_g3) "G3_Georeferenced_Unvouchered")) &
      qcf & is.finite(la) & is.finite(lo) & !ad & ((is.finite(pr) & pr >= 0.80) | (hg == "Divergent_NonHybrid_Lineage"))
    tx <- ifelse(nzchar(pd) & pd != "Unknown", pd, ifelse(hg == "Divergent_NonHybrid_Lineage", spv, "")); ok <- ok & nzchar(tx) & (tx != "Unknown")
    data.frame(lat = la[ok], lon = lo[ok], taxon = tx[ok], prob = ifelse(is.finite(pr[ok]), pr[ok], 0.85), azonal = .pick_num_prop(tbl, "is_azonal_lithology")[ok],
               cat = .pick_chr_prop(tbl, "catalogNumber", "")[ok], sid = .pick_chr_prop(tbl, "genetic_sample_id", "")[ok], df_row = if (track_row) which(ok) else rep(NA_integer_, sum(ok)), stringsAsFactors = FALSE)
  }
  anc_df <- extract_anc(df, FALSE, TRUE)
  if (!is.null(ga_df) && is.data.frame(ga_df) && nrow(ga_df) > 0L) {
    ga_rows <- extract_anc(ga_df, TRUE, FALSE)
    if (nrow(ga_rows) > 0L) {
      df_keys <- paste(anc_df$cat, anc_df$sid, sep = "::")[nzchar(anc_df$cat) | nzchar(anc_df$sid)]
      ga_keys <- paste(ga_rows$cat, ga_rows$sid, sep = "::")
      anc_df <- rbind(anc_df, ga_rows[!((nzchar(ga_rows$cat) | nzchar(ga_rows$sid)) & (ga_keys %in% df_keys)), , drop = FALSE])
    }
  }
  if (nrow(anc_df) == 0L) return(df)
  if (is.null(morph_X) || is.null(morph_S_inv)) {
    pool <- .compute_morph_pool_prop(df); if (is.null(pool)) return(df); morph_X <- pool$X; morph_S_inv <- pool$S_inv; is_basal <- pool$is_basal
  } else is_basal <- .is_basal_mask_prop(df)
  p <- ncol(morph_X); if (p < 2L) return(df)

  t0 <- resolve_tier0_anchor_overrides(df); mu_taxon <- .build_taxon_centroids(df, morph_X, is_basal, t0)
  for (tx in unique(anc_df$taxon)) if (is.null(mu_taxon[[tx]])) {
    rws <- anc_df$df_row[anc_df$taxon == tx & !is.na(anc_df$df_row)]; rws <- rws[is_basal[rws] & stats::complete.cases(morph_X[rws, , drop = FALSE])]
    if (length(rws) > 0L) mu_taxon[[tx]] <- colMeans(morph_X[rws, , drop = FALSE])
  }
  tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), ""); sp <- trimws(.pick_chr_prop(df, "species_standardized", ""))
  lat <- .pick_num_prop(df, "decimalLatitude"); lon <- .pick_num_prop(df, "decimalLongitude"); az <- .pick_num_prop(df, "is_azonal_lithology")
  ph <- .pick_num_prop(df, "fine_soil_ph"); ed_z <- .pick_num_prop(df, c("type_edaphic_z_score", "edaphic_z_score")); unc <- .pick_num_prop(df, "coordinateUncertaintyInMeters")
  doy <- .pick_num_prop(df, "hopkins_adjusted_doy"); rstat <- tolower(.pick_chr_prop(df, c("reproductive_trait_status", "capitulum_phenology_state"), ""))
  is_anth <- if ("is_active_anthesis" %in% names(df)) .safe_bool_prop(df$is_active_anthesis) else (is.finite(doy) & !grepl("vegetative|bud_only|post_fruiting", rstat))

  pure_ev <- t0$evict_from_given_anchor & !t0$is_tier0_dna_anchor; mu_doy <- list(); ph_ref <- list()
  for (tx in unique(c(names(mu_taxon), anc_df$taxon))) {
    a_idx <- which(t0$rectified_taxon == tx & !pure_ev & is_anth & is.finite(doy))
    if (length(a_idx) >= 2L) { ang <- doy[a_idx] * 2.0 * pi / 365.25; mu_doy[[tx]] <- ((atan2(mean(sin(ang)), mean(cos(ang))) * 365.25 / (2.0 * pi)) %% 365.25) }
    p_idx <- which(t0$rectified_taxon == tx & !pure_ev & is.finite(ph))
    if (length(p_idx) >= 2L) ph_ref[[tx]] <- c(mean(ph[p_idx]), max(stats::sd(ph[p_idx]), 0.4))
  }

  cands <- which(!(tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy")) & !.is_type_prop(df) & is.finite(lat) & is.finite(lon))
  df_geo <- which(is.finite(lat) & is.finite(lon) & nzchar(t0$rectified_taxon) & (t0$rectified_taxon != "Unknown") & !pure_ev)
  crit_d2 <- stats::qchisq(0.95, df = p); qc_num <- is.numeric(df$geno_quality_pass) || is.integer(df$geno_quality_pass); adm_num <- is.numeric(df$is_geno_admixed) || is.integer(df$is_geno_admixed)

  for (i in cands) {
    dists <- .haversine_km_vec(lat[i], lon[i], anc_df$lat, anc_df$lon); in_2km <- which(dists <= 2.0); if (length(in_2km) == 0L) next
    min_d <- min(dists[in_2km]); if (!((is.finite(unc[i]) && unc[i] <= 1000.0) || (is.na(unc[i]) && min_d <= 1.0))) next
    np <- in_2km[dists[in_2km] <= min_d + 1e-6]; near_a <- np[1L]
    if (length(np) > 1L && is_basal[i] && all(is.finite(morph_X[i, ]))) {
      cd2 <- vapply(np, function(a) if (is.null(mu_taxon[[anc_df$taxon[a]]])) Inf else { dv0 <- morph_X[i, ] - mu_taxon[[anc_df$taxon[a]]]; as.numeric(crossprod(dv0, morph_S_inv %*% dv0)) }, numeric(1L))
      if (any(is.finite(cd2))) near_a <- np[which.min(cd2)]
    }
    cand_sp <- anc_df$taxon[near_a]
    if (is.finite(az[i]) && is.finite(anc_df$azonal[near_a]) && as.integer(az[i]) != as.integer(anc_df$azonal[near_a])) next
    if (cand_sp == sp[i] && is.finite(ed_z[i]) && ed_z[i] > 1.75) next
    if (cand_sp != sp[i] && is.finite(ph[i]) && !is.null(ph_ref[[cand_sp]]) && abs(ph[i] - ph_ref[[cand_sp]][1L]) / ph_ref[[cand_sp]][2L] > 2.5) next
    if (!is_basal[i] || !all(is.finite(morph_X[i, ])) || is.null(mu_taxon[[cand_sp]])) next
    dv <- morph_X[i, ] - mu_taxon[[cand_sp]]; d2_cand <- as.numeric(crossprod(dv, morph_S_inv %*% dv)); df$deme_mahal_d2[i] <- d2_cand
    if (!is.finite(d2_cand) || d2_cand > crit_d2) next
    if (is_anth[i] && is.finite(doy[i]) && !is.null(mu_doy[[cand_sp]])) { dd <- abs(doy[i] - mu_doy[[cand_sp]]) %% 365.25; if (pmin(dd, 365.25 - dd) > 18.0) next }
    df_2km <- if (length(df_geo) > 0L) t0$rectified_taxon[df_geo[.haversine_km_vec(lat[i], lon[i], lat[df_geo], lon[df_geo]) <= 2.0]] else character(0)
    alt_taxa <- setdiff(unique(c(anc_df$taxon[dists <= 5.0], df_2km)), c(cand_sp, "", "Unknown"))
    if (length(alt_taxa) > 0L) {
      alt_d2 <- vapply(alt_taxa, function(atx) if (is.null(mu_taxon[[atx]])) NA_real_ else { av <- morph_X[i, ] - mu_taxon[[atx]]; as.numeric(crossprod(av, morph_S_inv %*% av)) }, numeric(1L))
      if (!any(is.finite(alt_d2)) || (min(alt_d2, na.rm = TRUE) - d2_cand) < 4.0) next
    }
    if (!("geno_linkage_tier" %in% names(df))) df$geno_linkage_tier <- rep(NA_character_, n)
    if (!("geno_predicted_taxon" %in% names(df))) df$geno_predicted_taxon <- rep(NA_character_, n)
    if (!("geno_quality_pass" %in% names(df))) df$geno_quality_pass <- rep(NA_integer_, n)
    if (!("geno_posterior_prob" %in% names(df))) df$geno_posterior_prob <- rep(NA_real_, n)
    df$geno_linkage_tier[i] <- "G2_Conditional_Deme_Hull"
    if ("linkage_tier" %in% names(df)) df$linkage_tier[i] <- "G2_Conditional_Deme_Hull"
    df$genomic_propagation_tier[i] <- "G2_Conditional_Deme_Hull"; df$geno_predicted_taxon[i] <- cand_sp
    df$geno_quality_pass[i] <- if (qc_num || !("geno_quality_pass" %in% names(df))) 1L else TRUE
    if ("is_geno_admixed" %in% names(df)) df$is_geno_admixed[i] <- if (adm_num) 0L else FALSE
    df$geno_posterior_prob[i] <- round(pmin(0.90, pmax(0.65, anc_df$prob[near_a] * (1.0 - 0.10 * (min_d / 2.0)))), 4L)
    if ("coalescent_gdi" %in% names(df)) df$coalescent_gdi[i] <- NA_real_
  }
  df
}

.fit_gda_predict <- function(X_tr, y_tr, X_ev, classes, lambda = 0.20) {
  d <- ncol(X_tr); k_cnt <- length(classes)
  for (j in seq_len(d)) { v <- X_tr[is.finite(X_tr[, j]), j]; med <- if (length(v) > 0L) stats::median(v) else 0.0; X_tr[!is.finite(X_tr[, j]), j] <- med; X_ev[!is.finite(X_ev[, j]), j] <- med }
  mu_g <- colMeans(X_tr); sd_g <- apply(X_tr, 2L, stats::sd); sd_g[!is.finite(sd_g) | sd_g < 1e-6] <- 1.0
  Z_tr <- sweep(sweep(X_tr, 2L, mu_g, "-"), 2L, sd_g, "/"); Z_ev <- sweep(sweep(X_ev, 2L, mu_g, "-"), 2L, sd_g, "/")
  mu_k <- matrix(0.0, nrow = k_cnt, ncol = d); n_k <- integer(k_cnt); res_list <- vector("list", k_cnt)
  for (idx in seq_len(k_cnt)) { sub <- Z_tr[y_tr == classes[idx], , drop = FALSE]; n_k[idx] <- nrow(sub); mu_k[idx, ] <- colMeans(sub); res_list[[idx]] <- sweep(sub, 2L, mu_k[idx, ], "-") }
  S_w <- crossprod(do.call(rbind, res_list)) / max(1.0, nrow(Z_tr) - k_cnt); S_w[!is.finite(S_w)] <- 0.0; S_w <- 0.5 * (S_w + t(S_w))
  S_lw <- (1.0 - lambda) * S_w + lambda * diag(pmax(diag(S_w), 1e-4), nrow = d) + diag(1e-5, d)
  W_inv <- tryCatch(solve(S_lw), error = function(e) diag(1.0 / pmax(diag(S_lw), 1e-4), nrow = d))
  priors <- pmax(n_k / sum(n_k), 1.0 / (2.0 * k_cnt)); priors <- priors / sum(priors)
  logits <- matrix(0.0, nrow = nrow(Z_ev), ncol = k_cnt); d2_mat <- matrix(0.0, nrow = nrow(Z_ev), ncol = k_cnt)
  for (idx in seq_len(k_cnt)) {
    a_k <- W_inv %*% mu_k[idx, ]; logits[, idx] <- as.numeric(Z_ev %*% a_k - 0.5 * sum(mu_k[idx, ] * a_k) + log(priors[idx]))
    diff_k <- sweep(Z_ev, 2L, mu_k[idx, ], "-"); d2_mat[, idx] <- rowSums((diff_k %*% W_inv) * diff_k)
  }
  list(logits = logits, d2 = d2_mat)
}

.softmax_T <- function(logits, T_val = 1.0) { sc <- logits / max(T_val, 1e-4); mx <- apply(sc, 1L, max); ex <- exp(sweep(sc, 1L, mx, "-")); sweep(ex, 1L, rowSums(ex), "/") }

# Step D: Tier_0 LOPO-CV Regularized GDA Projection & Abstention Gate
fit_tier0_calibrated_projection <- function(df, opts = list()) {
  df <- .init_propagation_cols(df, reset = FALSE); n <- nrow(df)
  empty_sum <- list(active = FALSE, lopo_cv_kappa = NA_real_, n_calibrated = 0L, n_abstained = 0L, confusion_matrix = NULL)
  t0 <- resolve_tier0_anchor_overrides(df); is_basal <- .is_basal_mask_prop(df)
  cand_feats <- intersect(c(paste0("PC", 1:6), paste0("ECT_PC", 1:3), "ect_lobe_complexity", "circularity", "aspect_ratio", "AR_inv_std", "tooth_density", "sinus_depth_ratio"), names(df))
  if (n == 0L || length(cand_feats) < 2L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }

  X_raw <- matrix(NA_real_, nrow = n, ncol = length(cand_feats))
  for (j in seq_along(cand_feats)) X_raw[, j] <- suppressWarnings(as.numeric(df[[cand_feats[j]]]))
  t0_raw <- which(t0$is_tier0_dna_anchor & is_basal); if (length(t0_raw) < 4L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }
  keep_f <- vapply(seq_len(ncol(X_raw)), function(j) { v <- X_raw[t0_raw, j]; mean(is.finite(v)) >= 0.60 && isTRUE(stats::sd(v, na.rm = TRUE) > 1e-8) }, logical(1L))
  if (sum(keep_f) < 2L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }
  X_use <- X_raw[, keep_f, drop = FALSE]; d <- ncol(X_use); row_fin <- apply(X_use, 1L, function(r) mean(is.finite(r)) >= 0.60); y_all <- trimws(t0$rectified_taxon)
  cls_tab <- table(y_all[t0_raw[row_fin[t0_raw]]]); classes <- sort(names(cls_tab[cls_tab >= 2L]))
  train_idx <- which(t0$is_tier0_dna_anchor & is_basal & row_fin & (y_all %in% classes))
  if (length(classes) < 2L || length(train_idx) < 4L) { attr(df, "tier0_projection_summary") <- empty_sum; return(df) }

  pop_col <- .pick_chr_prop(df, c("population_deme_id", "coancestry_cluster", "coancestry_cluster_id"), "")
  lat <- .pick_num_prop(df, "decimalLatitude"); lon <- .pick_num_prop(df, "decimalLongitude")
  pop_vec <- ifelse(nzchar(pop_col), pop_col, ifelse(is.finite(lat) & is.finite(lon), sprintf("%.1f_%.1f", round(lat, 1), round(lon, 1)), as.character(seq_len(n))))
  oof_logits <- matrix(NA_real_, nrow = length(train_idx), ncol = length(classes))
  for (pop in unique(pop_vec[train_idx])) {
    va_pos <- which(pop_vec[train_idx] == pop); tr_pos <- setdiff(seq_along(train_idx), va_pos)
    if (length(unique(y_all[train_idx[tr_pos]])) == length(classes)) {
      oof_logits[va_pos, ] <- .fit_gda_predict(X_use[train_idx[tr_pos], , drop = FALSE], y_all[train_idx[tr_pos]], X_use[train_idx[va_pos], , drop = FALSE], classes)$logits
    } else for (vp in va_pos) {
      tr1 <- setdiff(seq_along(train_idx), vp)
      if (length(unique(y_all[train_idx[tr1]])) == length(classes)) oof_logits[vp, ] <- .fit_gda_predict(X_use[train_idx[tr1], , drop = FALSE], y_all[train_idx[tr1]], X_use[train_idx[vp], , drop = FALSE], classes)$logits[1L, ]
    }
  }

  oof_ok <- stats::complete.cases(oof_logits); y_tr <- y_all[train_idx]; T_opt <- 1.0; kappa <- NA_real_
  cm <- matrix(0L, nrow = length(classes), ncol = length(classes), dimnames = list(classes, classes))
  if (sum(oof_ok) >= 4L) {
    y_idx <- match(y_tr[oof_ok], classes)
    T_opt <- max(1.0, min(5.0, stats::optimize(function(Tv) { pr <- .softmax_T(oof_logits[oof_ok, , drop = FALSE], Tv); -mean(log(pmax(pr[cbind(seq_len(sum(oof_ok)), y_idx)], 1e-12))) }, c(1.0, 5.0))$minimum))
    oof_pred <- classes[max.col(.softmax_T(oof_logits[oof_ok, , drop = FALSE], T_opt), ties.method = "first")]
    cm <- as.matrix(table(factor(y_tr[oof_ok], levels = classes), factor(oof_pred, levels = classes)))
    n_tot <- sum(cm); pe <- sum(rowSums(cm) * colSums(cm)) / max(1.0, n_tot^2)
    kappa <- if (n_tot == 0L) NA_real_ else if (abs(1.0 - pe) < 1e-12) 1.0 else (sum(diag(cm)) / n_tot - pe) / (1.0 - pe)
  }

  proj_idx <- which(is_basal & row_fin & !t0$is_tier0_dna_anchor)
  if (length(proj_idx) > 0L) {
    full_fit <- .fit_gda_predict(X_use[train_idx, , drop = FALSE], y_tr, X_use[proj_idx, , drop = FALSE], classes)
    probs <- .softmax_T(full_fit$logits, T_opt); crit_osr <- stats::qchisq(0.975, df = d); rs <- rowSums(cm)
    for (r in seq_along(proj_idx)) {
      i <- proj_idx[r]; p_row <- probs[r, ]; ord <- order(p_row, decreasing = TRUE); c1 <- ord[1L]; c2 <- ord[2L]; p1 <- p_row[c1]; p2 <- p_row[c2]
      pair_err <- if ((rs[c1] + rs[c2]) > 0L) (cm[c1, c2] + cm[c2, c1]) / (rs[c1] + rs[c2]) else 0.0
      min_d2 <- min(full_fit$d2[r, ])
      abstain <- (p1 < 0.70) || ((p1 - p2) < 0.15) || (pair_err > 0.20 && p1 < 0.90) || (!is.finite(min_d2) || min_d2 > crit_osr)
      df$tier0_calibrated_prob[i] <- round(p1, 4L); df$tier0_abstention_flag[i] <- if (abstain) 1L else 0L
      df$tier0_calibrated_taxon[i] <- if (abstain) "Sensu_Lato_DNA_Required" else classes[c1]
      if (!abstain && is.na(df$genomic_propagation_tier[i])) df$genomic_propagation_tier[i] <- "Tier0_Projected_Calibrated"
    }
  }
  attr(df, "tier0_projection_summary") <- list(
    active = TRUE, lopo_cv_kappa = kappa, confusion_matrix = cm,
    n_calibrated = sum(!is.na(df$tier0_calibrated_taxon) & df$tier0_abstention_flag == 0L, na.rm = TRUE),
    n_abstained = sum(df$tier0_abstention_flag == 1L, na.rm = TRUE)
  )
  df
}

# Step E: Propagate coalescent gdi to cohesive non-abstained clusters/lineages
propagate_cluster_gdi_governance <- function(df, gdi_df = NULL, morph_X = NULL, morph_S_inv = NULL) {
  df <- .init_propagation_cols(df, reset = FALSE); n <- nrow(df); if (n == 0L) return(df)
  tx_gdi <- list(); cl_gdi <- list(); mcl_med <- list()
  tier <- .pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), ""); qc <- .pick_bool_prop(df, "geno_quality_pass", FALSE)
  adm <- .pick_bool_prop(df, "is_geno_admixed", FALSE); pred <- trimws(.pick_chr_prop(df, "geno_predicted_taxon", ""))
  g_col <- .pick_num_prop(df, "coalescent_gdi"); cl_col <- trimws(.pick_chr_prop(df, c("coancestry_cluster_id", "coancestry_cluster"), ""))
  gmm_v <- .pick_num_prop(df, "gmm_cluster"); eco_v <- trimws(.pick_chr_prop(df, "ecotype_regime", ""))
  sp <- trimws(.pick_chr_prop(df, "species_standardized", ""))
  mcl_col <- ifelse(is.finite(gmm_v) & nzchar(sp) & sp != "Unknown", paste(sp, as.integer(gmm_v), eco_v, sep = "||"), "")
  src_ok <- (tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy")) & qc & !adm & is.finite(g_col) & nzchar(pred) & (pred != "Unknown")
  for (i in which(src_ok)) { tx_gdi[[pred[i]]] <- c(tx_gdi[[pred[i]]], g_col[i]); if (nzchar(cl_col[i])) cl_gdi[[cl_col[i]]] <- c(cl_gdi[[cl_col[i]]], g_col[i]) }
  for (mk in unique(mcl_col[src_ok & nzchar(mcl_col)])) {
    a_idx <- which(src_ok & mcl_col == mk); if (length(a_idx) < 2L) next
    tb <- sort(table(pred[a_idx]), decreasing = TRUE)
    if (length(tb) > 0L && (as.numeric(tb[1L]) / length(a_idx)) >= 0.75) mcl_med[[mk]] <- stats::median(g_col[a_idx[pred[a_idx] == names(tb)[1L]]])
  }
  if (!is.null(gdi_df) && is.data.frame(gdi_df) && nrow(gdi_df) > 0L) {
    gt <- trimws(.pick_chr_prop(gdi_df, c("taxon", "species_standardized", "lineage"), ""))
    gv <- .pick_num_prop(gdi_df, c("coalescent_gdi", "gdi")); gcl <- trimws(.pick_chr_prop(gdi_df, c("coancestry_cluster_id", "coancestry_cluster"), ""))
    for (r in which(is.finite(gv))) { if (nzchar(gt[r])) tx_gdi[[gt[r]]] <- c(tx_gdi[[gt[r]]], gv[r]); if (nzchar(gcl[r])) cl_gdi[[gcl[r]]] <- c(cl_gdi[[gcl[r]]], gv[r]) }
  }
  if (length(tx_gdi) == 0L && length(cl_gdi) == 0L && length(mcl_med) == 0L) return(df)
  tx_med <- lapply(tx_gdi, stats::median); cl_med <- lapply(cl_gdi, stats::median)
  if (is.null(morph_X) || is.null(morph_S_inv)) {
    pool <- .compute_morph_pool_prop(df)
    if (!is.null(pool)) { morph_X <- pool$X; morph_S_inv <- pool$S_inv; is_basal <- pool$is_basal } else is_basal <- rep(FALSE, n)
  } else is_basal <- .is_basal_mask_prop(df)
  t0 <- resolve_tier0_anchor_overrides(df)
  mu_taxon <- if (!is.null(morph_X)) .build_taxon_centroids(df, morph_X, is_basal, t0) else list()
  crit_d2 <- if (!is.null(morph_X)) stats::qchisq(0.95, df = ncol(morph_X)) else NA_real_
  cal_tx <- trimws(.pick_chr_prop(df, "tier0_calibrated_taxon", "")); abst <- df$tier0_abstention_flag
  pure_ev <- t0$evict_from_given_anchor & !t0$is_tier0_dna_anchor
  for (i in seq_len(n)) {
    if (adm[i] || pure_ev[i] || (!is.na(abst[i]) && abst[i] == 1L)) next
    gval <- if (nzchar(cl_col[i]) && !is.null(cl_med[[cl_col[i]]])) cl_med[[cl_col[i]]] else if (nzchar(mcl_col[i]) && !is.null(mcl_med[[mcl_col[i]]])) mcl_med[[mcl_col[i]]] else {
      cand <- ""
      if ((df$genomic_propagation_tier[i] %in% c("G1_Exsiccatae_Isoduplicate", "G2_Conditional_Deme_Hull") ||
           tier[i] %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy")) && nzchar(pred[i]) && pred[i] != "Unknown") {
        cand <- pred[i]
      } else if (!is.na(abst[i]) && abst[i] == 0L && nzchar(cal_tx[i]) && cal_tx[i] != "Sensu_Lato_DNA_Required") {
        cand <- cal_tx[i]
      } else if (nzchar(sp[i]) && !is.null(tx_med[[sp[i]]]) && is_basal[i] && !is.null(morph_X) && all(is.finite(morph_X[i, ])) && !is.null(mu_taxon[[sp[i]]])) {
        dv <- morph_X[i, ] - mu_taxon[[sp[i]]]; if (is.finite(d2 <- as.numeric(crossprod(dv, morph_S_inv %*% dv))) && d2 <= crit_d2) cand <- sp[i]
      }
      if (nzchar(cand) && !is.null(tx_med[[cand]])) tx_med[[cand]] else NA_real_
    }
    if (is.finite(gval)) {
      df$cluster_propagated_gdi[i] <- round(gval, 4L)
      if (is.na(df$genomic_propagation_tier[i]) && !(tier[i] %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy"))) df$genomic_propagation_tier[i] <- "Cluster_GDI_Governed"
    }
  }
  df
}

# Top-level orchestrator
run_step07_genomic_propagation <- function(df, opts = list()) {
  df <- .init_propagation_cols(df, reset = FALSE)
  ga_df <- .resolve_prop_csv(df, opts, "genomic_assignments_df", c("genomic_assignments_df", "ga_df"), "genomic_assignments.csv")
  gdi_df <- .resolve_prop_csv(df, opts, "genomic_gdi_delimitation", c("genomic_gdi_delimitation", "gdi_df"), "coalescent_gdi_delimitation.csv")
  df_act <- if (nrow(df) > 0L) any(.pick_bool_prop(df, "geno_quality_pass", FALSE) & (.pick_chr_prop(df, c("geno_linkage_tier", "linkage_tier"), "") %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy", "G1_Exsiccatae_Isoduplicate", "G2_Conditional_Deme_Hull"))) else FALSE
  ga_act <- (!is.null(ga_df) && is.data.frame(ga_df) && nrow(ga_df) > 0L) && any(.pick_bool_prop(ga_df, c("geno_quality_pass", "qc_pass"), FALSE) & (.pick_chr_prop(ga_df, c("geno_linkage_tier", "linkage_tier"), "") %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy", "G3_Georeferenced_Unvouchered", "G1_Exsiccatae_Isoduplicate", "G2_Conditional_Deme_Hull")))
  if (nrow(df) == 0L || (!df_act && !ga_act)) {
    attr(df, "genomic_propagation_summary") <- list(active = FALSE, n_g1_exsiccatae = 0L, n_g2_deme_hull = 0L, n_tier0_anchors = 0L, n_tier0_evicted = 0L, n_tier0_calibrated = 0L, n_tier0_abstained = 0L, n_cluster_gdi = 0L, lopo_cv_kappa = NA_real_, table = data.frame())
    return(df)
  }

  df <- .init_propagation_cols(df, reset = TRUE)
  df <- propagate_exsiccatae_isoduplicates(df, ga_df = ga_df)
  pool <- .compute_morph_pool_prop(df)
  df <- propagate_conditional_deme_priors(df, ga_df = ga_df, morph_X = pool$X, morph_S_inv = pool$S_inv)
  t0 <- resolve_tier0_anchor_overrides(df)
  df <- fit_tier0_calibrated_projection(df, opts = opts)
  df <- propagate_cluster_gdi_governance(df, gdi_df = gdi_df, morph_X = pool$X, morph_S_inv = pool$S_inv)

  proj_sum <- attr(df, "tier0_projection_summary"); sp <- trimws(.pick_chr_prop(df, "species_standardized", "Unknown")); sp[!nzchar(sp)] <- "Unknown"
  tx_levels <- sort(unique(c(sp, t0$rectified_taxon[nzchar(t0$rectified_taxon)]))); cnt_tx <- function(m) vapply(tx_levels, function(x) sum(m(x), na.rm = TRUE), integer(1L))
  pred_v <- .pick_chr_prop(df, "geno_predicted_taxon", "")
  tbl <- data.frame(
    taxon = tx_levels, n_total = cnt_tx(function(x) sp == x),
    n_tier0_anchors = cnt_tx(function(x) t0$is_tier0_dna_anchor & t0$rectified_taxon == x),
    n_tier0_evicted = cnt_tx(function(x) t0$evict_from_given_anchor & sp == x),
    n_g1_exsiccatae = cnt_tx(function(x) df$genomic_propagation_tier == "G1_Exsiccatae_Isoduplicate" & t0$rectified_taxon == x),
    n_g2_deme_hull = cnt_tx(function(x) df$genomic_propagation_tier == "G2_Conditional_Deme_Hull" & pred_v == x),
    n_tier0_calibrated = cnt_tx(function(x) df$tier0_calibrated_taxon == x & df$tier0_abstention_flag == 0L),
    n_tier0_abstained = cnt_tx(function(x) sp == x & df$tier0_abstention_flag == 1L),
    n_cluster_gdi = cnt_tx(function(x) t0$rectified_taxon == x & is.finite(df$cluster_propagated_gdi)), stringsAsFactors = FALSE
  )
  attr(df, "genomic_propagation_summary") <- list(
    active = TRUE, n_g1_exsiccatae = sum(df$genomic_propagation_tier == "G1_Exsiccatae_Isoduplicate", na.rm = TRUE),
    n_g2_deme_hull = sum(df$genomic_propagation_tier == "G2_Conditional_Deme_Hull", na.rm = TRUE),
    n_tier0_anchors = sum(t0$is_tier0_dna_anchor, na.rm = TRUE), n_tier0_evicted = sum(t0$evict_from_given_anchor, na.rm = TRUE),
    n_tier0_calibrated = sum(!is.na(df$tier0_calibrated_taxon) & df$tier0_abstention_flag == 0L, na.rm = TRUE),
    n_tier0_abstained = sum(df$tier0_abstention_flag == 1L, na.rm = TRUE), n_cluster_gdi = sum(is.finite(df$cluster_propagated_gdi)),
    lopo_cv_kappa = if (!is.null(proj_sum)) proj_sum$lopo_cv_kappa else NA_real_, table = tbl
  )
  df
}
