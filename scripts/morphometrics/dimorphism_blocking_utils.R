#!/usr/bin/env Rscript
# ==============================================================================
# Module: scripts/morphometrics/dimorphism_blocking_utils.R
# Project: Plantmatics Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Within-taxon dioecious sexual dimorphism (staminate vs. pistillate) and
#   seasonal/heteroblastic shoot-phase blocking (chasmogamous vs. cleistogamous,
#   vernal vs. autumnal, sterile_frond vs. fertile_frond) for Step 04.
# ==============================================================================

.canon_sex_state_r <- function(raw_sex, raw_text = NULL) {
  n <- length(raw_sex)
  out <- rep(NA_character_, n)
  sx <- tolower(trimws(as.character(raw_sex)))
  sx[is.na(raw_sex) | !nzchar(sx) | sx %in% c("na", "nan", "null", "unknown")] <- ""
  has_sx <- nzchar(sx)
  out[has_sx & grepl("(^|[^a-z])(pistillate|female|carpellate|gynoecious)([^a-z]|$)|\u2640", sx)] <- "pistillate"
  out[has_sx & is.na(out) & grepl("(^|[^a-z])(staminate|male|androecious)([^a-z]|$)|\u2642", sx)] <- "staminate"
  if (!is.null(raw_text) && length(raw_text) == n) {
    tx <- tolower(trimws(as.character(raw_text)))
    tx[is.na(raw_text)] <- ""
    need <- is.na(out) & nzchar(tx)
    if (any(need)) {
      out[need & grepl("(^|[^a-z])(pistillate|female|carpellate|gynoecious)([^a-z]|$)|\u2640", tx)] <- "pistillate"
      still <- is.na(out) & nzchar(tx)
      out[still & grepl("(^|[^a-z])(staminate|male|androecious)([^a-z]|$)|\u2642", tx)] <- "staminate"
    }
  }
  out
}

.canon_shoot_phase_r <- function(raw_phase, raw_text = NULL) {
  n <- length(raw_phase)
  out <- rep(NA_character_, n)
  ph <- tolower(trimws(as.character(raw_phase)))
  ph[is.na(raw_phase) | !nzchar(ph) | ph %in% c("na", "nan", "null", "unknown")] <- ""
  classify_v <- function(v) {
    res <- rep(NA_character_, length(v))
    ok <- nzchar(v)
    if (!any(ok)) return(res)
    res[ok & grepl("cleistogam", v)] <- "cleistogamous"
    res[ok & is.na(res) & grepl("chasmogam", v)] <- "chasmogamous"
    res[ok & is.na(res) & grepl("fertile[[:space:]_-]*frond|sporophyll|sporangiate|fertile[[:space:]_-]*leaf", v)] <- "fertile_frond"
    res[ok & is.na(res) & grepl("sterile[[:space:]_-]*frond|trophophyll|vegetative[[:space:]_-]*frond", v)] <- "sterile_frond"
    res[ok & is.na(res) & grepl("autumnal|fall[[:space:]_-]*shoot|late[[:space:]_-]*season", v)] <- "autumnal"
    res[ok & is.na(res) & grepl("vernal|spring[[:space:]_-]*shoot|early[[:space:]_-]*season", v)] <- "vernal"
    res
  }
  out <- classify_v(ph)
  if (!is.null(raw_text) && length(raw_text) == n) {
    tx <- tolower(trimws(as.character(raw_text)))
    tx[is.na(raw_text)] <- ""
    need <- is.na(out) & nzchar(tx)
    if (any(need)) out[need] <- classify_v(tx[need])
  }
  out
}

merge_sex_and_shoot_phase_columns_r <- function(df, vouchers_df = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(df)
  n <- nrow(df)
  for (cn in c("specimen_sex", "seasonal_shoot_phase")) {
    if (!(cn %in% names(df))) df[[cn]] <- rep(NA_character_, n) else df[[cn]] <- as.character(df[[cn]])
  }
  if (!is.null(vouchers_df) && is.data.frame(vouchers_df) && nrow(vouchers_df) > 0L &&
      "catalogNumber" %in% names(vouchers_df) && "catalogNumber" %in% names(df)) {
    v_cat <- trimws(as.character(vouchers_df$catalogNumber))
    d_cat <- trimws(as.character(df$catalogNumber))
    has_v_pid <- "plant_individual_id" %in% names(vouchers_df)
    v_pid <- if (has_v_pid) ifelse(is.na(vouchers_df$plant_individual_id), 0L, suppressWarnings(as.integer(vouchers_df$plant_individual_id))) else rep(0L, nrow(vouchers_df))
    d_pid <- if ("plant_individual_id" %in% names(df)) ifelse(is.na(df$plant_individual_id), 0L, suppressWarnings(as.integer(df$plant_individual_id))) else rep(0L, n)
    v_pid[is.na(v_pid)] <- 0L; d_pid[is.na(d_pid)] <- 0L
    k_v <- paste(v_cat, v_pid, sep = "__")
    k_d <- paste(d_cat, d_pid, sep = "__")
    n_pid_v <- tapply(v_pid, v_cat, function(z) length(unique(z)))
    single_v <- !(d_cat %in% names(n_pid_v[n_pid_v > 1L]))
    m_exact <- match(k_d, k_v)
    m_sheet <- match(d_cat, v_cat)
    m_use <- ifelse(!is.na(m_exact), m_exact, ifelse(single_v, m_sheet, NA_integer_))
    for (vc in intersect(c("specimen_sex", "seasonal_shoot_phase", "sex", "reproductiveCondition", "occurrenceRemarks", "fieldNotes", "dynamicProperties"), names(vouchers_df))) {
      v_vals <- as.character(vouchers_df[[vc]][m_use])
      if (!(vc %in% names(df))) {
        df[[vc]] <- v_vals
      } else {
        miss <- is.na(df[[vc]]) | !nzchar(trimws(as.character(df[[vc]])))
        if (any(miss)) df[[vc]][miss] <- v_vals[miss]
      }
    }
  }
  dwc_sex <- if ("sex" %in% names(df)) as.character(df$sex) else rep("", n)
  dwc_rep <- if ("reproductiveCondition" %in% names(df)) as.character(df$reproductiveCondition) else rep("", n)
  dwc_txt <- paste(
    ifelse(is.na(dwc_sex), "", dwc_sex),
    ifelse(is.na(dwc_rep), "", dwc_rep),
    if ("occurrenceRemarks" %in% names(df)) ifelse(is.na(df$occurrenceRemarks), "", as.character(df$occurrenceRemarks)) else "",
    if ("fieldNotes" %in% names(df)) ifelse(is.na(df$fieldNotes), "", as.character(df$fieldNotes)) else "",
    if ("dynamicProperties" %in% names(df)) ifelse(is.na(df$dynamicProperties), "", as.character(df$dynamicProperties)) else "",
    sep = " ; "
  )
  df$specimen_sex <- .canon_sex_state_r(df$specimen_sex, dwc_txt)
  df$seasonal_shoot_phase <- .canon_shoot_phase_r(df$seasonal_shoot_phase, paste(ifelse(is.na(dwc_rep), "", dwc_rep), dwc_txt, sep = " ; "))
  df
}

build_sex_and_shoot_phase_covariates_r <- function(df, ref_idx = NULL) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) return(list())
  df <- merge_sex_and_shoot_phase_columns_r(df)
  n <- nrow(df)
  if (is.null(ref_idx) || length(ref_idx) == 0L) ref_idx <- seq_len(n)
  ref_idx <- ref_idx[ref_idx >= 1L & ref_idx <= n]
  if (length(ref_idx) < 6L) return(list())

  fb_sp <- paste0(if (exists("resolve_r_target_genus", mode = "function")) resolve_r_target_genus() else "Packera", "_sp")
  sp_vec <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else if ("scientificName" %in% names(df)) as.character(df$scientificName) else rep(fb_sp, n)
  sp_vec[is.na(sp_vec) | !nzchar(trimws(sp_vec))] <- fb_sp

  aw <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n)
  aw[!is.finite(aw) | aw <= 0.0] <- 1.0
  sb <- if ("spatial_balance_weight" %in% names(df)) suppressWarnings(as.numeric(df$spatial_balance_weight)) else rep(1.0, n)
  sb[!is.finite(sb) | sb <= 0.0] <- 1.0
  eff_w <- aw * sb

  specs <- list(
    sex_pistillate_ind = c("specimen_sex", "staminate", "pistillate"),
    phase_cleistogamous_ind = c("seasonal_shoot_phase", "chasmogamous", "cleistogamous"),
    phase_autumnal_ind = c("seasonal_shoot_phase", "vernal", "autumnal"),
    phase_fertile_frond_ind = c("seasonal_shoot_phase", "sterile_frond", "fertile_frond")
  )
  out_covs <- list()
  u_taxa <- unique(sp_vec[ref_idx])
  for (cov_nm in names(specs)) {
    sp_spec <- specs[[cov_nm]]
    st <- as.character(df[[sp_spec[1L]]])
    s0 <- sp_spec[2L]; s1 <- sp_spec[3L]
    dim_taxa <- character(0)
    for (tx in u_taxa) {
      r_tx <- ref_idx[sp_vec[ref_idx] == tx]
      if (sum(!is.na(st[r_tx]) & st[r_tx] == s0) >= 3L && sum(!is.na(st[r_tx]) & st[r_tx] == s1) >= 3L) {
        dim_taxa <- c(dim_taxa, tx)
      }
    }
    if (length(dim_taxa) == 0L) next
    vec <- rep(0.0, n)
    for (tx in dim_taxa) {
      r_tx_k <- ref_idx[sp_vec[ref_idx] == tx & !is.na(st[ref_idx]) & (st[ref_idx] %in% c(s0, s1))]
      raw_r <- ifelse(st[r_tx_k] == s1, 1.0, 0.0)
      w_r <- eff_w[r_tx_k]
      if (sum(w_r) <= 1e-12) w_r <- rep(1.0, length(w_r))
      mu_tx <- stats::weighted.mean(raw_r, w = w_r)
      all_tx_k <- which(sp_vec == tx & !is.na(st) & (st %in% c(s0, s1)))
      vec[all_tx_k] <- ifelse(st[all_tx_k] == s1, 1.0, 0.0) - mu_tx
    }
    if (is.finite(stats::sd(vec[ref_idx])) && stats::sd(vec[ref_idx]) > 1e-5) {
      out_covs[[cov_nm]] <- vec
    }
  }
  out_covs
}

has_active_dimorphism_covariates_r <- function(use_covs) {
  is.list(use_covs) && length(use_covs) > 0L && any(
    names(use_covs) %in% c("sex_pistillate_ind", "phase_cleistogamous_ind", "phase_autumnal_ind", "phase_fertile_frond_ind")
  )
}
