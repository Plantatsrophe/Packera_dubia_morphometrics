# Step 07 non-circular cytological integration, heteroploid barrier gate, and gigas audit.
# Cytology never enters Step 03/04/05 feature matrices; it acts solely in Step 07.

.norm_cyto_str <- function(x) {
  out <- trimws(as.character(x)); out[is.na(out)] <- ""; out
}

.norm_exsiccatae_cyto <- function(x) tolower(gsub("\\s+", " ", .norm_cyto_str(x)))

.parse_ploidy_int <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s[is.na(x) | s %in% c("", "na", "nan", "none", "null", "mixed/uncertain")] <- NA_character_
  vals <- suppressWarnings(as.integer(round(as.numeric(sub("x$", "", s)))))
  vals[is.na(vals) | vals < 2L | vals > 16L] <- NA_integer_; vals
}

.parse_2n_int <- function(x) {
  vals <- suppressWarnings(as.integer(round(as.numeric(as.character(x)))))
  vals[is.na(vals) | vals < 4L | vals > 300L] <- NA_integer_; vals
}

.safe_bool_cyto <- function(x) {
  if (is.null(x)) return(logical(0))
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  s <- tolower(trimws(as.character(x))); !is.na(s) & s %in% c("1", "true", "t", "yes")
}

.pick_col_chr <- function(df, candidates, default = "") {
  for (nm in candidates) if (nm %in% names(df)) return(.norm_cyto_str(df[[nm]]))
  rep(default, nrow(df))
}

.pick_col_raw <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  NULL
}

.init_cytology_cols <- function(df) {
  n <- nrow(df)
  if (!"cytology_linkage_tier" %in% names(df)) df$cytology_linkage_tier <- rep(NA_character_, n)
  if (!"voucher_sporophytic_2n" %in% names(df)) df$voucher_sporophytic_2n <- rep(NA_integer_, n)
  if (!"voucher_ploidy_level" %in% names(df)) df$voucher_ploidy_level <- rep(NA_integer_, n)
  if (!"taxon_reported_cytotypes" %in% names(df)) df$taxon_reported_cytotypes <- rep(NA_character_, n)
  if (!"is_multi_cytotype_taxon" %in% names(df)) df$is_multi_cytotype_taxon <- rep(NA_integer_, n)
  if (!"cytotype_reproductive_barrier_pass" %in% names(df)) df$cytotype_reproductive_barrier_pass <- rep(NA_integer_, n)
  df
}

.load_cytology_table <- function(df_opt, conn, db_path, is_demo, tbl_name, csv_path) {
  if (is.data.frame(df_opt)) return(df_opt)
  has_dbi <- requireNamespace("DBI", quietly = TRUE)
  if (!is.null(conn) && has_dbi && isTRUE(tryCatch(DBI::dbIsValid(conn), error = function(e) FALSE))) {
    res <- tryCatch(if (DBI::dbExistsTable(conn, tbl_name)) DBI::dbReadTable(conn, tbl_name) else NULL, error = function(e) NULL)
    if (is.data.frame(res)) return(res)
  }
  if (!isTRUE(is_demo) && is.character(db_path) && length(db_path) == 1L && nzchar(db_path) && file.exists(db_path) && has_dbi) {
    res <- tryCatch({
      c2 <- if (exists("connect_packera_db", mode = "function")) connect_packera_db(db_path = db_path, read_only = TRUE, demo = FALSE)
      else if (requireNamespace("RSQLite", quietly = TRUE)) DBI::dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RO) else NULL
      if (!is.null(c2)) {
        on.exit(tryCatch(DBI::dbDisconnect(c2), error = function(e) NULL), add = TRUE)
        if (DBI::dbExistsTable(c2, tbl_name)) DBI::dbReadTable(c2, tbl_name) else NULL
      } else NULL
    }, error = function(e) NULL)
    if (is.data.frame(res)) return(res)
  }
  if (is.character(csv_path) && length(csv_path) == 1L && nzchar(csv_path) && file.exists(csv_path)) {
    res <- tryCatch(utils::read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.data.frame(res)) return(res)
  }
  data.frame()
}

.unanimous_or_na_int <- function(vals) {
  u <- unique(vals[!is.na(vals)])
  if (length(u) == 1L) as.integer(u[1L]) else NA_integer_
}

.resolve_primary_cluster <- function(sp_rows, cl_vec, df) {
  fin_cl <- cl_vec[sp_rows]; fin_cl <- fin_cl[!is.na(fin_cl)]
  if (length(fin_cl) == 0L) return(NA_integer_)
  t_stat <- tolower(.pick_col_chr(df, c("type_status"), "")[sp_rows])
  a_tier <- tolower(.pick_col_chr(df, c("authority_tier"), "")[sp_rows])
  is_t <- grepl("\\b(holo|iso|lecto|neo|syn)type\\b", t_stat) | (a_tier == "tier_0_type")
  type_cls <- cl_vec[sp_rows][is_t & !is.na(cl_vec[sp_rows])]
  if (length(type_cls) > 0L) {
    tb <- sort(table(type_cls), decreasing = TRUE)
    return(min(suppressWarnings(as.integer(names(tb)[tb == tb[1L]])), na.rm = TRUE))
  }
  if (1L %in% fin_cl) return(1L)
  tb <- sort(table(fin_cl), decreasing = TRUE)
  min(suppressWarnings(as.integer(names(tb)[tb == tb[1L]])), na.rm = TRUE)
}

.summarize_cluster_anchors <- function(idx, is_typed, v_pl, row_bx, d_ex, d_cat) {
  t_idx <- idx[is_typed[idx] & (!is.na(v_pl[idx]) | !is.na(row_bx[idx]))]
  if (length(t_idx) == 0L) return(list(n_anchors = 0L, ploidies = integer(0), base_xs = integer(0)))
  anc_id <- ifelse(nzchar(d_ex[t_idx]), paste0("EX:", d_ex[t_idx]), ifelse(nzchar(d_cat[t_idx]), paste0("CAT:", d_cat[t_idx]), paste0("ROW:", t_idx)))
  list(n_anchors = length(unique(anc_id)), ploidies = sort(unique(v_pl[t_idx][!is.na(v_pl[t_idx])])), base_xs = sort(unique(row_bx[t_idx][!is.na(row_bx[t_idx])])))
}

.compute_gigas_audit <- function(df, v_pl, sp_vec) {
  traits <- intersect(c("area_mm2", "petiole_width_mm", "log1p_capitulum_count"), names(df))
  empty_df <- data.frame(
    species_standardized = character(0), trait = character(0), ploidy_low = integer(0), ploidy_high = integer(0),
    n_low = integer(0), n_high = integer(0), median_low = numeric(0), median_high = numeric(0),
    effect_ratio = numeric(0), wilcox_p = numeric(0), stringsAsFactors = FALSE
  )
  if (length(traits) == 0L || nrow(df) == 0L) return(empty_df)
  rows_out <- list()
  for (sp in unique(sp_vec[nzchar(sp_vec) & sp_vec != "Unknown" & !is.na(v_pl)])) {
    sp_idx <- which(sp_vec == sp & !is.na(v_pl)); pls <- sort(unique(v_pl[sp_idx])); if (length(pls) < 2L) next
    p_lo <- pls[1L]
    for (p_hi in pls[-1L]) {
      i_lo <- sp_idx[v_pl[sp_idx] == p_lo]; i_hi <- sp_idx[v_pl[sp_idx] == p_hi]
      for (tr in traits) {
        x_lo <- suppressWarnings(as.numeric(df[[tr]][i_lo])); x_lo <- x_lo[is.finite(x_lo)]
        x_hi <- suppressWarnings(as.numeric(df[[tr]][i_hi])); x_hi <- x_hi[is.finite(x_hi)]
        if (length(x_lo) == 0L || length(x_hi) == 0L) next
        m_lo <- stats::median(x_lo); m_hi <- stats::median(x_hi)
        eff <- if (tr == "log1p_capitulum_count") (m_hi - m_lo) else if (abs(m_lo) > 1e-9) (m_hi / m_lo) else NA_real_
        wp <- if (length(x_lo) >= 2L && length(x_hi) >= 2L) suppressWarnings(tryCatch(stats::wilcox.test(x_hi, x_lo)$p.value, error = function(e) NA_real_)) else NA_real_
        rows_out[[length(rows_out) + 1L]] <- data.frame(
          species_standardized = sp, trait = tr, ploidy_low = as.integer(p_lo), ploidy_high = as.integer(p_hi),
          n_low = length(x_lo), n_high = length(x_hi), median_low = round(m_lo, 4), median_high = round(m_hi, 4),
          effect_ratio = round(eff, 4), wilcox_p = signif(wp, 4), stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows_out) == 0L) empty_df else do.call(rbind, rows_out)
}

run_step07_cytology_integration <- function(df, conn = NULL, cytology_dir = "data/cytology", opts = list()) {
  if (!is.data.frame(df)) return(df)
  df <- .init_cytology_cols(df); n <- nrow(df)
  if (n == 0L) {
    attr(df, "cytology_markdown_lines") <- character(0); attr(df, "cytotype_barrier_labels") <- character(0)
    return(df)
  }

  rec_df <- .load_cytology_table(opts$cytological_records_df, conn, opts$db_path, opts$demo, "cytological_records", file.path(cytology_dir, "public_packera_cytology_records.csv"))
  tax_df <- .load_cytology_table(opts$taxon_cytology_summary_df, conn, opts$db_path, opts$demo, "taxon_cytology_summary", file.path(cytology_dir, "public_packera_taxon_cytology_summary.csv"))
  has_rec <- is.data.frame(rec_df) && nrow(rec_df) > 0L
  has_tax <- is.data.frame(tax_df) && nrow(tax_df) > 0L

  nq_pl_raw <- .pick_col_raw(df, c("inferred_cytotype"))
  nq_pl <- if (is.null(nq_pl_raw)) rep(NA_integer_, n) else .parse_ploidy_int(nq_pl_raw)
  nq_dl_raw <- .pick_col_raw(df, c("nquire_delta_logl"))
  nq_dl <- if (is.null(nq_dl_raw)) rep(NA_real_, n) else suppressWarnings(as.numeric(as.character(nq_dl_raw)))
  g_tier <- .pick_col_chr(df, c("geno_linkage_tier", "linkage_tier"), "")
  g_qc <- if (any(c("geno_quality_pass", "qc_pass") %in% names(df))) .safe_bool_cyto(.pick_col_raw(df, c("geno_quality_pass", "qc_pass"))) else rep(TRUE, n)
  q12_raw <- .pick_col_raw(df, c("interclass_het_Q12", "interspecific_ancestry_q12"))
  q12_v <- if (is.null(q12_raw)) rep(NA_real_, n) else suppressWarnings(as.numeric(as.character(q12_raw)))
  g_adm_recent <- if ("is_geno_admixed" %in% names(df)) (.safe_bool_cyto(df$is_geno_admixed) & (is.na(q12_v) | q12_v >= 0.35)) else rep(FALSE, n)
  hq_nq_mask <- (g_tier %in% c("G1_Direct_Voucher", "G1_Exsiccatae_Isoduplicate")) &
    g_qc & !g_adm_recent & !is.na(nq_pl) & is.finite(nq_dl) &
    ((nq_pl == 2L & nq_dl <= -4.0) | (nq_pl >= 3L & nq_dl >= 4.0))

  if (!has_rec && !has_tax && !any(hq_nq_mask)) {
    attr(df, "cytology_markdown_lines") <- character(0); attr(df, "cytotype_barrier_labels") <- character(0)
    return(df)
  }

  d_cat <- .pick_col_chr(df, c("catalogNumber"), "")
  d_ex <- .norm_exsiccatae_cyto(.pick_col_chr(df, c("exsiccatae_key"), ""))
  d_gb <- .pick_col_chr(df, c("genbank_accession"), "")
  d_sp <- .pick_col_chr(df, c("species_standardized"), "")
  tier <- rep(NA_character_, n); v_2n <- rep(NA_integer_, n); v_pl <- rep(NA_integer_, n)
  row_bx <- rep(NA_integer_, n); is_typed <- rep(FALSE, n)
  r_pl <- integer(0); r_2n <- integer(0); r_bx <- integer(0); r_sp <- character(0)

  if (has_rec) {
    r_cat <- .pick_col_chr(rec_df, c("catalogNumber"), ""); r_ex <- .norm_exsiccatae_cyto(.pick_col_chr(rec_df, c("exsiccatae_key"), ""))
    r_gb <- .pick_col_chr(rec_df, c("genbank_accession"), ""); r_sp <- .pick_col_chr(rec_df, c("species_standardized"), "")
    r_pl <- .parse_ploidy_int(.pick_col_raw(rec_df, c("ploidy_level", "voucher_ploidy_level")))
    r_2n <- .parse_2n_int(.pick_col_raw(rec_df, c("sporophytic_2n", "chromosome_count_2n", "voucher_sporophytic_2n")))
    r_bx_raw <- suppressWarnings(as.integer(round(as.numeric(.pick_col_raw(rec_df, c("base_x"))))))
    r_bx <- if (length(r_bx_raw) == nrow(rec_df)) r_bx_raw else rep(NA_integer_, nrow(rec_df))
    calc_bx <- is.na(r_bx) & !is.na(r_2n) & !is.na(r_pl) & r_pl >= 2L & (r_2n %% r_pl == 0L)
    r_bx[calc_bx] <- as.integer(r_2n[calc_bx] %/% r_pl[calc_bx])
    r_vok <- !grepl("^C[34]_", toupper(.pick_col_chr(rec_df, c("cytology_linkage_tier", "linkage_tier"), "")))

    for (i in seq_len(n)) {
      c_hits <- if (nzchar(d_cat[i])) which(r_vok & r_cat == d_cat[i]) else integer(0)
      if (length(c_hits) == 0L && nzchar(d_gb[i])) c_hits <- which(r_vok & r_gb == d_gb[i])
      if (length(c_hits) > 0L) {
        tier[i] <- "C1_Direct_Cytovoucher"; v_pl[i] <- .unanimous_or_na_int(r_pl[c_hits])
        v_2n[i] <- .unanimous_or_na_int(r_2n[c_hits]); row_bx[i] <- .unanimous_or_na_int(r_bx[c_hits])
        is_typed[i] <- !is.na(v_pl[i]) || !is.na(row_bx[i])
      } else if (nzchar(d_ex[i])) {
        e_hits <- which(r_vok & r_ex == d_ex[i])
        if (length(e_hits) > 0L) {
          pl_e <- .unanimous_or_na_int(r_pl[e_hits]); bx_e <- .unanimous_or_na_int(r_bx[e_hits])
          if (!is.na(pl_e) || !is.na(bx_e)) {
            tier[i] <- "C2_Exsiccatae_Cytoduplicate"; v_pl[i] <- pl_e
            v_2n[i] <- .unanimous_or_na_int(r_2n[e_hits]); row_bx[i] <- bx_e; is_typed[i] <- TRUE
          }
        }
      }
    }
  }

  is_nq_anchor <- rep(FALSE, n)
  nq_direct <- which(!is_typed & is.na(tier) & hq_nq_mask)
  if (length(nq_direct) > 0L) {
    v_pl[nq_direct] <- nq_pl[nq_direct]; is_typed[nq_direct] <- TRUE; is_nq_anchor[nq_direct] <- TRUE
  }
  if (any(hq_nq_mask)) {
    t_stat_all <- tolower(.pick_col_chr(df, c("type_status"), ""))
    a_tier_all <- tolower(.pick_col_chr(df, c("authority_tier", "determiner_tier"), ""))
    is_type_col <- if ("is_type_specimen" %in% names(df)) .safe_bool_cyto(df$is_type_specimen) else rep(FALSE, n)
    is_nom_type <- is_type_col | grepl("(^|[^a-z])(holo|iso|lecto|neo|syn)type([^a-z]|$)", t_stat_all) | (a_tier_all == "tier_0_type")
    for (j in which(!is_typed & is.na(tier) & nzchar(d_ex) & !is_nom_type)) {
      nq_ex_hits <- which(hq_nq_mask & d_ex == d_ex[j])
      if (length(nq_ex_hits) > 0L) {
        pl_nq_ex <- .unanimous_or_na_int(nq_pl[nq_ex_hits])
        if (!is.na(pl_nq_ex)) {
          v_pl[j] <- pl_nq_ex; is_typed[j] <- TRUE; is_nq_anchor[j] <- TRUE
        }
      }
    }
  }

  tax_rep_map <- list(); tax_multi_map <- list()
  if (has_tax) {
    t_sp <- .pick_col_chr(tax_df, c("species_standardized"), "")
    t_rep <- .pick_col_chr(tax_df, c("taxon_reported_cytotypes", "reported_cytotypes"), "")
    t_mod <- .parse_ploidy_int(.pick_col_raw(tax_df, c("modal_ploidy_level")))
    t_mult <- .safe_bool_cyto(.pick_col_raw(tax_df, c("is_multi_cytotype_taxon")))
    for (j in which(nzchar(t_sp))) {
      sp <- t_sp[j]
      tax_rep_map[[sp]] <- if (nzchar(t_rep[j])) t_rep[j] else if (!is.na(t_mod[j])) paste0(t_mod[j], "x") else ""
      tax_multi_map[[sp]] <- isTRUE(t_mult[j])
    }
  }

  all_taxa <- unique(d_sp[nzchar(d_sp) & d_sp != "Unknown"])
  tx_cyt_str <- rep(NA_character_, n); tx_multi_vec <- rep(NA_integer_, n)
  for (sp in all_taxa) {
    idx_sp <- which(d_sp == sp)
    obs_pl <- unique(c(v_pl[idx_sp][is_typed[idx_sp] & !is.na(v_pl[idx_sp])], if (length(r_sp) > 0L) r_pl[r_sp == sp & !is.na(r_pl)] else integer(0)))
    obs_bx <- unique(c(row_bx[idx_sp][is_typed[idx_sp] & !is.na(row_bx[idx_sp])], if (length(r_sp) > 0L) r_bx[r_sp == sp & !is.na(r_bx)] else integer(0)))
    rep_txt <- tax_rep_map[[sp]]; has_tax_entry <- !is.null(rep_txt)
    if (!has_tax_entry && length(obs_pl) == 0L && length(obs_bx) == 0L) next
    parsed_rep <- if (has_tax_entry && nzchar(rep_txt)) .parse_ploidy_int(unlist(strsplit(rep_txt, "[,;/ ]+"))) else integer(0)
    all_pl <- sort(unique(c(parsed_rep[!is.na(parsed_rep)], obs_pl[!is.na(obs_pl)])))
    rep_final <- if (length(all_pl) > 0L) paste(paste0(all_pl, "x"), collapse = ", ") else if (has_tax_entry && nzchar(rep_txt)) rep_txt else NA_character_
    is_m <- isTRUE(tax_multi_map[[sp]]) || length(all_pl) >= 2L || length(unique(obs_bx[!is.na(obs_bx)])) >= 2L
    tx_cyt_str[idx_sp] <- rep_final; tx_multi_vec[idx_sp] <- if (is_m) 1L else 0L
    need_tier <- idx_sp[is.na(tier[idx_sp])]
    if (length(need_tier) > 0L) tier[need_tier] <- if (is_m) "C4_Taxon_MultiCytotype" else "C3_Taxon_Monomorphic"
  }

  barrier_pass <- rep(0L, n); barrier_labels <- character(0)
  if ("gmm_cluster" %in% names(df)) {
    cl_vec <- suppressWarnings(as.integer(.pick_col_raw(df, c("gmm_cluster"))))
    for (sp in all_taxa) {
      sp_rows <- which(d_sp == sp & !is.na(cl_vec)); if (length(unique(cl_vec[sp_rows])) < 2L) next
      prim_cl <- .resolve_primary_cluster(sp_rows, cl_vec, df); if (is.na(prim_cl)) next
      p_sum <- .summarize_cluster_anchors(sp_rows[cl_vec[sp_rows] == prim_cl], is_typed, v_pl, row_bx, d_ex, d_cat)
      if (p_sum$n_anchors < 1L || (length(p_sum$ploidies) == 0L && length(p_sum$base_xs) == 0L)) next
      for (sec_cl in setdiff(sort(unique(cl_vec[sp_rows])), prim_cl)) {
        sec_idx <- sp_rows[cl_vec[sp_rows] == sec_cl]
        s_sum <- .summarize_cluster_anchors(sec_idx, is_typed, v_pl, row_bx, d_ex, d_cat)
        if (s_sum$n_anchors < 2L || length(s_sum$ploidies) > 1L || length(s_sum$base_xs) > 1L || (length(s_sum$ploidies) == 0L && length(s_sum$base_xs) == 0L)) next
        pl_disjoint <- length(s_sum$ploidies) == 1L && length(p_sum$ploidies) >= 1L && !(s_sum$ploidies[1L] %in% p_sum$ploidies)
        bx_disjoint <- length(s_sum$base_xs) == 1L && length(p_sum$base_xs) >= 1L && !(s_sum$base_xs[1L] %in% p_sum$base_xs)
        if (!pl_disjoint && !bx_disjoint) next
        barrier_pass[sec_idx] <- 1L
        barrier_labels[paste0(sp, "||", sec_cl)] <- if (pl_disjoint) {
          sprintf("%dx vs %dx", s_sum$ploidies[1L], p_sum$ploidies[1L])
        } else if (length(s_sum$ploidies) == 1L && length(p_sum$ploidies) >= 1L) {
          sprintf("%dx (x=%d) vs %dx (x=%d)", s_sum$ploidies[1L], s_sum$base_xs[1L], p_sum$ploidies[1L], p_sum$base_xs[1L])
        } else {
          sprintf("x=%d vs x=%d", s_sum$base_xs[1L], p_sum$base_xs[1L])
        }
      }
    }
  }

  gigas_df <- .compute_gigas_audit(df, ifelse(is_typed, v_pl, NA_integer_), d_sp)
  csv_target <- if (is.character(opts$cytology_summary_path) && nzchar(opts$cytology_summary_path)) opts$cytology_summary_path
  else if (is.character(opts$report_out) && nzchar(opts$report_out)) file.path(dirname(opts$report_out), "cytotype_concordance_summary.csv")
  else if (is.character(opts$reports_dir) && nzchar(opts$reports_dir)) file.path(opts$reports_dir, "cytotype_concordance_summary.csv")
  else if (is.character(opts$output_dir) && nzchar(opts$output_dir)) file.path(opts$output_dir, "cytotype_concordance_summary.csv") else ""
  if (nzchar(csv_target) && nrow(gigas_df) > 0L) {
    tryCatch({ dir.create(dirname(csv_target), recursive = TRUE, showWarnings = FALSE); utils::write.csv(gigas_df, csv_target, row.names = FALSE) }, error = function(e) NULL)
  }

  n_c1 <- sum(!is.na(tier) & tier == "C1_Direct_Cytovoucher"); n_c2 <- sum(!is.na(tier) & tier == "C2_Exsiccatae_Cytoduplicate")
  n_nq <- sum(is_nq_anchor); n_mc <- length(unique(d_sp[!is.na(tx_multi_vec) & tx_multi_vec == 1L]))
  nq_txt <- if (n_nq > 0L) sprintf(", %d Step06 `nQuire` high-confidence anchor(s)", n_nq) else ""
  md_lines <- c(sprintf("- **Step 07 Cytological Linkage:** %d `C1` direct cytovouchers, %d `C2` exsiccatae cytoduplicates%s, %d typed anchors total across %d multi-cytotype taxa.", n_c1, n_c2, nq_txt, sum(is_typed), n_mc))
  if (length(barrier_labels) > 0L) md_lines <- c(md_lines, sprintf("- **Postzygotic Heteroploid Barriers:** %d cohesive cluster(s) passed heteroploid/dysploid isolation (%s).", length(barrier_labels), paste(paste(names(barrier_labels), barrier_labels, sep = ": "), collapse = "; ")))
  if (nrow(gigas_df) > 0L) md_lines <- c(md_lines, sprintf("- **Post-Hoc Gigas Morphometric Audit:** Evaluated %d cytotype-trait comparison(s) across `area_mm2`, `petiole_width_mm`, and capitula without circular trait feedback.", nrow(gigas_df)))

  df$cytology_linkage_tier <- as.character(tier)
  df$voucher_sporophytic_2n <- as.integer(v_2n)
  df$voucher_ploidy_level <- as.integer(v_pl)
  df$taxon_reported_cytotypes <- as.character(tx_cyt_str)
  df$is_multi_cytotype_taxon <- as.integer(tx_multi_vec)
  df$cytotype_reproductive_barrier_pass <- as.integer(barrier_pass)
  attr(df, "cytotype_barrier_labels") <- barrier_labels
  attr(df, "cytology_markdown_lines") <- md_lines
  attr(df, "gigas_concordance_df") <- gigas_df
  df
}
