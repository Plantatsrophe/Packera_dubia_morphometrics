# ==============================================================================
# Module: scripts/analysis/step07_organellar_triage_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Quality-gated organellar reference placement hydration, small-batch
#   reference cytonuclear capture flagging, confirmed Stage 2 hybrid maternal
#   cpDNA polarity annotation, and non-mutating out-of-panel ITS + waxy cues.
# ==============================================================================

.organellar_safe_bool <- function(x) {
  if (exists(".safe_bool_vec_r", mode = "function")) return(.safe_bool_vec_r(x))
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)
  toupper(trimws(as.character(x))) %in% c("1", "TRUE", "T", "YES", "Y")
}

.canon_packera_key_r <- function(x) {
  s <- tolower(trimws(gsub("[[:space:]]+", " ", gsub("_+", " ", as.character(x)))))
  s[is.na(s)] <- ""
  s <- ifelse(s == "packera millefolia", "packera millefolium", s)
  s <- ifelse(s %in% c("packera memmingeri", "packera x memmingeri", "packera xmemmingeri", "packera ×memmingeri", "packera × memmingeri"), "packera x memmingeri", s)
  ifelse(s %in% c("packera tomentosa", "senecio tomentosus", "packera dubia"), "packera dubia", s)
}

.is_valid_packera_name_r <- function(x) {
  s <- trimws(as.character(x))
  k <- .canon_packera_key_r(s)
  !is.na(s) & grepl("^Packera [[:alpha:]x]", s) & !(k %in% c("", "packera sp", "packera sp.", "packera unknown"))
}

.resolve_organellar_placement_df <- function(df, opts = list()) {
  if (!is.null(opts$organellar_placement_df) && is.data.frame(opts$organellar_placement_df)) {
    return(opts$organellar_placement_df)
  }
  at <- attr(df, "public_organellar_placement_df")
  if (!is.null(at) && is.data.frame(at)) return(at)
  base_dirs <- c(
    if (exists(".resolve_step07_genomic_tables_dir", mode = "function")) .resolve_step07_genomic_tables_dir(opts) else character(0),
    if (exists(".resolve_step07_genomic_reports_dir", mode = "function")) .resolve_step07_genomic_reports_dir(opts) else character(0),
    if (!is.null(opts$genomic_dir) && nzchar(opts$genomic_dir)) opts$genomic_dir else character(0),
    if (!is.null(opts$genomics_dir) && nzchar(opts$genomics_dir)) opts$genomics_dir else character(0),
    if (!is.null(opts$out_dir) && nzchar(opts$out_dir)) opts$out_dir else character(0),
    if (!is.null(opts$output_dir) && nzchar(opts$output_dir)) opts$output_dir else character(0),
    if (!is.null(opts$output_queue) && nzchar(opts$output_queue)) dirname(opts$output_queue) else character(0),
    if (!is.null(opts$output_report) && nzchar(opts$output_report)) dirname(opts$output_report) else character(0)
  )
  cands <- character(0)
  for (d in unique(base_dirs[!is.na(base_dirs) & nzchar(base_dirs)])) {
    cands <- c(cands, d, file.path(d, "genomics"), file.path(dirname(d), "genomics"), file.path(dirname(d), "reports"))
  }
  for (d in unique(cands)) {
    p <- file.path(d, "public_organellar_placement.csv")
    if (file.exists(p)) {
      res <- tryCatch(utils::read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(res) && is.data.frame(res)) return(res)
    }
  }
  NULL
}

hydrate_step07_organellar_placement_metrics <- function(df, opts = list()) {
  if (!is.data.frame(df)) return(df)
  n <- nrow(df)
  metric_cols <- list(
    plastid_kmer_containment = NA_real_, plastid_shared_11mers = NA_integer_,
    its_kmer_containment = NA_real_, its_shared_11mers = NA_integer_,
    waxy_kmer_containment = NA_real_, waxy_shared_11mers = NA_integer_
  )
  for (cn in names(metric_cols)) {
    if (!(cn %in% names(df))) {
      df[[cn]] <- rep(metric_cols[[cn]], n)
    } else if (grepl("11mers$", cn)) {
      df[[cn]] <- suppressWarnings(as.integer(df[[cn]]))
    } else {
      df[[cn]] <- suppressWarnings(as.numeric(df[[cn]]))
    }
  }
  if (n == 0L) return(df)

  if ("genetic_sample_id" %in% names(df) && "exsiccatae_key" %in% names(df)) {
    sid_v <- trimws(ifelse(is.na(df$genetic_sample_id), "", as.character(df$genetic_sample_id)))
    ex_v <- trimws(ifelse(is.na(df$exsiccatae_key), "", as.character(df$exsiccatae_key)))
    tier_v <- if ("geno_linkage_tier" %in% names(df)) as.character(df$geno_linkage_tier) else if ("genomic_propagation_tier" %in% names(df)) as.character(df$genomic_propagation_tier) else rep("", n)
    need_sid <- which((!nzchar(sid_v)) & nzchar(ex_v) & (tier_v == "G1_Exsiccatae_Isoduplicate"))
    for (ri in need_sid) {
      don <- which(nzchar(sid_v) & (ex_v == ex_v[ri]))
      if (length(don) > 0L) {
        df$genetic_sample_id[ri] <- sid_v[don[1L]]
        sid_v[ri] <- sid_v[don[1L]]
      }
    }
  }

  pop_df <- .resolve_organellar_placement_df(df, opts)
  if (is.null(pop_df) || nrow(pop_df) == 0L || !("genetic_sample_id" %in% names(df)) ||
      !all(c("genetic_sample_id", "locus") %in% names(pop_df))) {
    return(df)
  }
  df_sid <- trimws(ifelse(is.na(df$genetic_sample_id), "", as.character(df$genetic_sample_id)))
  pop_df$genetic_sample_id <- trimws(as.character(pop_df$genetic_sample_id))
  pop_df$locus <- tolower(trimws(as.character(pop_df$locus)))
  pop_cont <- if ("kmer_containment" %in% names(pop_df)) suppressWarnings(as.numeric(pop_df$kmer_containment)) else rep(NA_real_, nrow(pop_df))
  pop_k11 <- if ("n_shared_11mers" %in% names(pop_df)) suppressWarnings(as.integer(pop_df$n_shared_11mers)) else rep(NA_integer_, nrow(pop_df))
  pop_tx <- if ("best_match_taxon" %in% names(pop_df)) trimws(as.character(pop_df$best_match_taxon)) else rep(NA_character_, nrow(pop_df))

  locus_map <- list(
    plastid = c("plastid_kmer_containment", "plastid_shared_11mers", "public_plastome_match"),
    its = c("its_kmer_containment", "its_shared_11mers", "public_its_match"),
    waxy = c("waxy_kmer_containment", "waxy_shared_11mers", "waxy_lineage_match")
  )
  for (loc in names(locus_map)) {
    l_idx <- which(pop_df$locus == loc & nzchar(pop_df$genetic_sample_id))
    if (length(l_idx) == 0L) next
    ord <- order(-ifelse(is.finite(pop_cont[l_idx]), pop_cont[l_idx], -1), -ifelse(is.finite(pop_k11[l_idx]), pop_k11[l_idx], -1L))
    l_idx <- l_idx[ord]
    l_idx <- l_idx[!duplicated(pop_df$genetic_sample_id[l_idx])]
    m <- match(df_sid, pop_df$genetic_sample_id[l_idx])
    hit <- !is.na(m) & nzchar(df_sid)
    if (!any(hit)) next
    cols <- locus_map[[loc]]
    fill_c <- hit & !is.finite(df[[cols[1L]]])
    if (any(fill_c)) df[[cols[1L]]][fill_c] <- pop_cont[l_idx[m[fill_c]]]
    fill_k <- hit & is.na(df[[cols[2L]]])
    if (any(fill_k)) df[[cols[2L]]][fill_k] <- pop_k11[l_idx[m[fill_k]]]
    if (!(cols[3L] %in% names(df))) df[[cols[3L]]] <- rep(NA_character_, n)
    cur_tx <- ifelse(is.na(df[[cols[3L]]]), "", trimws(as.character(df[[cols[3L]]])))
    fill_t <- hit & !nzchar(cur_tx) & !is.na(pop_tx[l_idx[m]]) & nzchar(pop_tx[l_idx[m]])
    if (any(fill_t)) df[[cols[3L]]][fill_t] <- pop_tx[l_idx[m[fill_t]]]
  }
  df
}

evaluate_step07_organellar_and_oop_signals <- function(df) {
  n <- if (is.data.frame(df)) nrow(df) else 0L
  if (n == 0L) {
    return(list(
      ref_cytonuc_capture = logical(0), ref_cytonuc_note = character(0),
      hq_donor_hint = character(0), oop_rationale_note = character(0),
      format_stage2_maternal_note = function(row_i, p_a, p_b) ""
    ))
  }
  nms <- names(df)
  chr_col <- function(cands) {
    for (cn in cands) if (cn %in% nms) return(trimws(as.character(df[[cn]])))
    rep(NA_character_, n)
  }
  num_col <- function(cands) {
    for (cn in cands) if (cn %in% nms) return(suppressWarnings(as.numeric(df[[cn]])))
    rep(NA_real_, n)
  }
  given_sp <- chr_col(c("species_standardized"))
  geno_pred <- chr_col(c("geno_predicted_taxon"))
  gp2 <- chr_col(c("geno_parent_2"))
  cp_tx <- chr_col(c("public_plastome_match"))
  its_tx <- chr_col(c("public_its_match"))
  wx_tx <- chr_col(c("waxy_lineage_match"))
  ref_taxa <- unique(c(given_sp, geno_pred, chr_col(c("geno_parent_1")), gp2))
  if (exists(".align_packera_orthography", mode = "function")) {
    geno_pred <- .align_packera_orthography(geno_pred, ref_taxa)
    gp2 <- .align_packera_orthography(gp2, ref_taxa)
    cp_tx <- .align_packera_orthography(cp_tx, ref_taxa)
    its_tx <- .align_packera_orthography(its_tx, ref_taxa)
    wx_tx <- .align_packera_orthography(wx_tx, ref_taxa)
  }

  cp_cont <- num_col(c("plastid_kmer_containment"))
  cp_k11 <- num_col(c("plastid_shared_11mers"))
  its_cont <- num_col(c("its_kmer_containment"))
  its_k11 <- num_col(c("its_shared_11mers"))
  wx_cont <- num_col(c("waxy_kmer_containment"))
  wx_k11 <- num_col(c("waxy_shared_11mers"))
  geno_prob <- num_col(c("geno_posterior_prob"))
  q12_v <- num_col(c("interspecific_ancestry_q12", "interclass_het_Q12"))
  is_adm <- if ("is_geno_admixed" %in% nms) .organellar_safe_bool(df$is_geno_admixed) else rep(FALSE, n)
  is_qry <- if ("is_putative_novel_or_hybrid_query" %in% nms) .organellar_safe_bool(df$is_putative_novel_or_hybrid_query) else rep(FALSE, n)
  qc_ok <- if ("geno_quality_pass" %in% nms) .organellar_safe_bool(df$geno_quality_pass) else rep(TRUE, n)
  tier_col <- if ("geno_linkage_tier" %in% nms) trimws(as.character(df$geno_linkage_tier)) else if ("linkage_tier" %in% nms) trimws(as.character(df$linkage_tier)) else rep("G1_Direct_Voucher", n)
  active_tier <- tier_col %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy")

  cp_hq <- .is_valid_packera_name_r(cp_tx) & is.finite(cp_cont) & (cp_cont >= 0.80) & is.finite(cp_k11) & (cp_k11 >= 80)
  its_hq <- .is_valid_packera_name_r(its_tx) & is.finite(its_cont) & (its_cont >= 0.80) & is.finite(its_k11) & (its_k11 >= 80)
  waxy_hq <- .is_valid_packera_name_r(wx_tx) & is.finite(wx_cont) & (wx_cont >= 0.80) & is.finite(wx_k11) & (wx_k11 >= 80)

  cp_k <- .canon_packera_key_r(cp_tx)
  its_k <- .canon_packera_key_r(its_tx)
  wx_k <- .canon_packera_key_r(wx_tx)
  pred_k <- .canon_packera_key_r(geno_pred)
  given_k <- .canon_packera_key_r(given_sp)

  pure_nuc <- (!is.finite(q12_v) & !is_adm) | (is.finite(q12_v) & (q12_v < 0.15))
  nuc_matches_cp <- (its_hq & (its_k == cp_k)) | (waxy_hq & (wx_k == cp_k))
  ref_cytonuc_capture <- active_tier & qc_ok & cp_hq & .is_valid_packera_name_r(geno_pred) &
    is.finite(geno_prob) & (geno_prob >= 0.80) & pure_nuc & (cp_k != pred_k) & !nuc_matches_cp
  ref_cytonuc_note <- ifelse(ref_cytonuc_capture, sprintf("cpDNA=%s (cont=%.2f, 11mers=%d)", cp_tx, cp_cont, as.integer(round(cp_k11))), "")

  proven_adm <- (is.finite(q12_v) & (q12_v >= 0.35)) | is_adm
  gp2_missing <- is.na(gp2) | !nzchar(gp2) | (gp2 == "Unknown") | (.canon_packera_key_r(gp2) == given_k)
  hq_donor_hint <- rep(NA_character_, n)
  for (i in which(proven_adm & gp2_missing)) {
    if (its_hq[i] && its_k[i] != given_k[i]) {
      hq_donor_hint[i] <- its_tx[i]
    } else if (waxy_hq[i] && wx_k[i] != given_k[i]) {
      hq_donor_hint[i] <- wx_tx[i]
    } else if (cp_hq[i] && cp_k[i] != given_k[i]) {
      hq_donor_hint[i] <- cp_tx[i]
    }
  }

  seq_tier_ok <- tier_col %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate", "G2_Population_Voucher_Proxy", "G3_Georeferenced_Unvouchered")
  seq_row <- qc_ok & seq_tier_ok & (.is_valid_packera_name_r(geno_pred) | ("genetic_sample_id" %in% nms & !is.na(df$genetic_sample_id) & nzchar(trimws(as.character(df$genetic_sample_id)))))
  ref_anchor_seq <- seq_row & .is_valid_packera_name_r(geno_pred) & !is_qry & !is_adm
  panel_raw <- c(geno_pred[ref_anchor_seq], given_sp[ref_anchor_seq])
  panel_keys <- unique(.canon_packera_key_r(panel_raw[.is_valid_packera_name_r(panel_raw)]))
  oop_rationale_note <- rep("", n)
  if (length(panel_keys) > 0L) {
    oop_mask <- seq_row & its_hq & is.finite(its_cont) & (its_cont >= 0.85) &
      waxy_hq & is.finite(wx_cont) & (wx_cont >= 0.85) & (its_k == wx_k) & !(its_k %in% panel_keys)
    for (oi in which(oop_mask)) {
      oop_rationale_note[oi] <- sprintf(
        " [Out-of-Panel Nuclear Locus Match: ITS+waxy=%s (ITS_cont=%.2f, waxy_cont=%.2f); taxon absent from NextRAD reference panel]",
        its_tx[oi], its_cont[oi], wx_cont[oi]
      )
    }
  }

  format_stage2_maternal_note <- function(row_i, p_a, p_b) {
    if (is.na(row_i) || row_i < 1L || row_i > n || !isTRUE(cp_hq[row_i])) return("")
    if (is.na(p_a) || !nzchar(trimws(p_a)) || is.na(p_b) || !nzchar(trimws(p_b))) return("")
    pa_k <- .canon_packera_key_r(p_a)
    pb_k <- .canon_packera_key_r(p_b)
    if (cp_k[row_i] == pa_k) {
      mat_sp <- p_a; pat_sp <- p_b; pat_k <- pb_k
    } else if (cp_k[row_i] == pb_k) {
      mat_sp <- p_b; pat_sp <- p_a; pat_k <- pa_k
    } else {
      return("")
    }
    nuc_ok <- (isTRUE(its_hq[row_i]) && its_k[row_i] == pat_k) || (isTRUE(waxy_hq[row_i]) && wx_k[row_i] == pat_k)
    nuc_corrob_str <- if (nuc_ok) ", ITS/waxy corroborated" else ""
    sprintf(
      " [Maternal cpDNA Parent: %s (cont=%.2f, 11mers=%d); Putative Pollen Donor: %s%s]",
      mat_sp, cp_cont[row_i], as.integer(round(cp_k11[row_i])), pat_sp, nuc_corrob_str
    )
  }

  list(
    ref_cytonuc_capture = ref_cytonuc_capture,
    ref_cytonuc_note = ref_cytonuc_note,
    hq_donor_hint = hq_donor_hint,
    oop_rationale_note = oop_rationale_note,
    format_stage2_maternal_note = format_stage2_maternal_note
  )
}

