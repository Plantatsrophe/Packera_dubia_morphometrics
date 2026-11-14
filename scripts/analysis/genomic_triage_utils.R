# ==============================================================================
# Module: scripts/analysis/genomic_triage_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Helper utilities bridging the independent Step 06 Genomic & Phylogenomic
#   Evidence Stream into Step 07 Multi-Evidence Synthesis without altering upstream
#   Step 03-05 modalities or legacy 3-stream scores on unsequenced vouchers.
# ==============================================================================

.safe_bool_vec_r <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x) || is.integer(x)) return(!is.na(x) & x != 0)
  s <- toupper(trimws(as.character(x)))
  s %in% c("1", "TRUE", "T", "YES", "Y")
}

.resolve_step07_genomic_tables_dir <- function(opts) {
  q_dir <- if (!is.null(opts$output_queue) && nzchar(opts$output_queue)) dirname(opts$output_queue) else "outputs/tables"
  cands <- unique(c(
    q_dir,
    if (!is.null(opts$genomic_dir) && nzchar(opts$genomic_dir)) opts$genomic_dir else character(0),
    file.path(dirname(q_dir), "reports"),
    file.path(dirname(dirname(q_dir)), "outputs", "reports")
  ))
  for (cd in cands) {
    if (nzchar(cd) && (file.exists(file.path(cd, "genomic_assignments.csv")) || file.exists(file.path(cd, "genomic_model_metrics.csv")))) {
      return(cd)
    }
  }
  q_dir
}

.resolve_step07_genomic_reports_dir <- function(opts) {
  q_dir <- if (!is.null(opts$output_queue) && nzchar(opts$output_queue)) dirname(opts$output_queue) else "outputs/tables"
  rep_dir <- if (!is.null(opts$output_report) && nzchar(opts$output_report)) dirname(opts$output_report) else character(0)
  cands <- unique(c(
    if (!is.null(opts$genomic_dir) && nzchar(opts$genomic_dir)) opts$genomic_dir else character(0),
    rep_dir,
    file.path(dirname(q_dir), "reports"),
    file.path(dirname(dirname(q_dir)), "outputs", "reports"),
    q_dir
  ))
  for (cd in cands) {
    if (nzchar(cd) && (
      file.exists(file.path(cd, "coalescent_gdi_delimitation.csv")) ||
      file.exists(file.path(cd, "landscape_mmrr_feems_summary.csv")) ||
      file.exists(file.path(cd, "genomic_assignments.csv"))
    )) {
      return(cd)
    }
  }
  .resolve_step07_genomic_tables_dir(opts)
}

has_active_genomic_row_mask <- function(df) {
  n <- nrow(df)
  if (n == 0L) return(logical(0))
  nms <- names(df)
  if (!("geno_predicted_taxon" %in% nms)) return(rep(FALSE, n))
  pred <- as.character(df$geno_predicted_taxon)
  pred_ok <- !is.na(pred) & nzchar(trimws(pred)) & (pred != "Unknown")
  tier_col <- if ("geno_linkage_tier" %in% nms) as.character(df$geno_linkage_tier) else if ("linkage_tier" %in% nms) as.character(df$linkage_tier) else rep("G1_Direct_Voucher", n)
  tier_ok <- tier_col %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy")
  qc_ok <- if ("geno_quality_pass" %in% nms) .safe_bool_vec_r(df$geno_quality_pass) else rep(TRUE, n)
  prob_ok <- if ("geno_posterior_prob" %in% nms) is.finite(suppressWarnings(as.numeric(df$geno_posterior_prob))) else rep(TRUE, n)
  pred_ok & tier_ok & qc_ok & prob_ok
}

load_step07_genomic_metrics <- function(opts, df) {
  t_dir <- .resolve_step07_genomic_tables_dir(opts)
  r_dir <- .resolve_step07_genomic_reports_dir(opts)
  ga_path <- file.path(t_dir, "genomic_assignments.csv")
  if (!file.exists(ga_path)) ga_path <- file.path(r_dir, "genomic_assignments.csv")
  gm_path <- file.path(t_dir, "genomic_model_metrics.csv")
  if (!file.exists(gm_path)) gm_path <- file.path(r_dir, "genomic_model_metrics.csv")

  geno_keep_schema <- c(
    "catalogNumber", "genetic_sample_id", "geno_linkage_tier", "geno_quality_pass",
    "is_putative_novel_or_hybrid_query", "geno_predicted_taxon", "geno_posterior_prob",
    "geno_margin", "genetic_cluster", "q_max", "q_admixture_json", "is_geno_admixed",
    "geno_parent_1", "geno_parent_2", "nearest_tree_sister_taxon",
    "patristic_dist_to_dubia", "replicate_pair_concordance",
    "hdplot_retained_loci", "inferred_cytotype", "nquire_delta_logl",
    "plastid_clade", "its_ribotype", "cytonuclear_discordant",
    "hybrid_index_h", "interspecific_ancestry_q12", "hybrid_generation_class",
    "coancestry_cluster", "f_branch_max", "f_branch_donor_taxon",
    "snaq_gamma_minor", "coalescent_gdi", "bpp_introgression_bf",
    "mmrr_beta_ibd", "mmrr_beta_ibe", "mmrr_p_ibe"
  )

  has_geno_cols <- ("geno_predicted_taxon" %in% names(df)) && any(!is.na(df$geno_predicted_taxon) & nzchar(as.character(df$geno_predicted_taxon)))
  missing_phase2_cols <- setdiff(geno_keep_schema, names(df))
  need_csv_hydrate <- file.exists(ga_path) && ("catalogNumber" %in% names(df)) && (!has_geno_cols || length(missing_phase2_cols) > 0L)

  if (need_csv_hydrate) {
    ga_df <- tryCatch(utils::read.csv(ga_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(ga_df) && nrow(ga_df) > 0L && "catalogNumber" %in% names(ga_df)) {
      ga_df$catalogNumber <- as.character(ga_df$catalogNumber)
      pass_ga <- ga_df[
        !is.na(ga_df$catalogNumber) & nzchar(trimws(ga_df$catalogNumber)) &
          .safe_bool_vec_r(ga_df$geno_quality_pass) &
          as.character(ga_df$linkage_tier) %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy"),
        , drop = FALSE
      ]
      if (nrow(pass_ga) > 0L) {
        tier_ord <- ifelse(pass_ga$linkage_tier == "G1_Direct_Voucher", 1L, 2L)
        hq_ord <- ifelse(is.finite(suppressWarnings(as.numeric(pass_ga$hq_reads))), suppressWarnings(as.numeric(pass_ga$hq_reads)), 0.0)
        pass_ga <- pass_ga[order(tier_ord, -hq_ord, pass_ga$genetic_sample_id), , drop = FALSE]
        pass_ga <- pass_ga[!duplicated(pass_ga$catalogNumber), , drop = FALSE]
        if ("linkage_tier" %in% names(pass_ga) && !("geno_linkage_tier" %in% names(pass_ga))) pass_ga$geno_linkage_tier <- pass_ga$linkage_tier
        keep_cols <- intersect(geno_keep_schema, names(pass_ga))
        m_idx <- match(as.character(df$catalogNumber), pass_ga$catalogNumber)
        hit <- !is.na(m_idx)
        if (any(hit)) {
          for (cn in setdiff(keep_cols, "catalogNumber")) {
            if (!(cn %in% names(df))) {
              df[[cn]] <- NA
              df[[cn]][hit] <- pass_ga[[cn]][m_idx[hit]]
            } else if (!has_geno_cols) {
              df[[cn]][hit] <- pass_ga[[cn]][m_idx[hit]]
            } else {
              fill_idx <- hit & is.na(df[[cn]])
              if (any(fill_idx)) df[[cn]][fill_idx] <- pass_ga[[cn]][m_idx[fill_idx]]
            }
          }
        }
      }
    }
  }

  read_opt_csv <- function(fname) {
    p <- file.path(r_dir, fname)
    if (!file.exists(p)) p <- file.path(t_dir, fname)
    if (file.exists(p)) tryCatch(utils::read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  }

  attr(df, "genomic_model_metrics") <- if (file.exists(gm_path)) tryCatch(utils::read.csv(gm_path, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  attr(df, "genomic_gdi_delimitation") <- read_opt_csv("coalescent_gdi_delimitation.csv")
  attr(df, "genomic_fbranch_summary") <- read_opt_csv("introgressive_fbranch_summary.csv")
  attr(df, "genomic_landscape_summary") <- read_opt_csv("landscape_mmrr_feems_summary.csv")
  df
}

compute_genomic_kappa_weight <- function(df, metrics_df = NULL) {
  active_mask <- has_active_genomic_row_mask(df)
  if (!any(active_mask)) return(NA_real_)
  if (is.null(metrics_df)) metrics_df <- attr(df, "genomic_model_metrics")
  if (!is.null(metrics_df) && nrow(metrics_df) > 0L && "loro_cv_kappa" %in% names(metrics_df)) {
    k_val <- suppressWarnings(as.numeric(metrics_df$loro_cv_kappa[1L]))
    if (is.finite(k_val)) return(max(0.05, min(1.0, k_val)))
  }
  sub_g <- df[active_mask, , drop = FALSE]
  yt <- as.character(sub_g$species_standardized)
  yp <- as.character(sub_g$geno_predicted_taxon)
  ok <- !is.na(yt) & nzchar(yt) & yt != "Unknown" & !is.na(yp) & nzchar(yp)
  if (sum(ok) >= 3L && length(unique(yt[ok])) >= 2L && exists("compute_cohens_kappa", mode = "function")) {
    return(max(0.05, min(1.0, compute_cohens_kappa(yt[ok], yp[ok], default_kappa = 0.85))))
  }
  0.85
}

export_genomic_field_triage_csv <- function(opts, triage_df) {
  t_dir <- .resolve_step07_genomic_tables_dir(opts)
  r_dir <- .resolve_step07_genomic_reports_dir(opts)
  ga_path <- file.path(t_dir, "genomic_assignments.csv")
  if (!file.exists(ga_path)) ga_path <- file.path(r_dir, "genomic_assignments.csv")
  if (!file.exists(ga_path)) return(invisible(NULL))
  ga_df <- tryCatch(utils::read.csv(ga_path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(ga_df) || nrow(ga_df) == 0L) return(invisible(NULL))

  ga_df$catalogNumber <- as.character(ga_df$catalogNumber)
  n_g <- nrow(ga_df)
  ga_df$voucher_recommended_determination <- NA_character_
  ga_df$voucher_treatment_subtype <- NA_character_
  ga_df$voucher_triage_priority <- NA_character_
  ga_df$genomic_field_status <- rep("Genomic_Reference_Sample", n_g)

  if (!is.null(triage_df) && nrow(triage_df) > 0L && "catalogNumber" %in% names(triage_df)) {
    tri_u <- triage_df[!duplicated(triage_df$catalogNumber), , drop = FALSE]
    m <- match(ga_df$catalogNumber, as.character(tri_u$catalogNumber))
    hit <- !is.na(m)
    if (any(hit)) {
      ga_df$voucher_recommended_determination[hit] <- as.character(tri_u$recommended_determination[m[hit]])
      ga_df$voucher_treatment_subtype[hit] <- as.character(tri_u$taxonomic_treatment_subtype[m[hit]])
      ga_df$voucher_triage_priority[hit] <- as.character(tri_u$synthesis_triage_priority[m[hit]])
    }
  }

  qc_ok <- .safe_bool_vec_r(ga_df$geno_quality_pass)
  is_q <- .safe_bool_vec_r(ga_df$is_putative_novel_or_hybrid_query)
  hyb_cls <- if ("hybrid_generation_class" %in% names(ga_df)) as.character(ga_df$hybrid_generation_class) else rep(NA_character_, n_g)
  is_div <- !is.na(hyb_cls) & (hyb_cls == "Divergent_NonHybrid_Lineage")
  is_adm <- if ("is_geno_admixed" %in% names(ga_df)) .safe_bool_vec_r(ga_df$is_geno_admixed) & !is_div else rep(FALSE, n_g)
  sp_std <- as.character(ga_df$species_standardized)
  g_pred <- as.character(ga_df$geno_predicted_taxon)
  g_prob <- suppressWarnings(as.numeric(ga_df$geno_posterior_prob))

  ga_df$genomic_field_status[!qc_ok] <- "Failed_Extraction_QC_Dropout"
  ga_df$genomic_field_status[qc_ok & is_q & is_div] <- "Divergent_NonHybrid_Lineage_Query"
  ga_df$genomic_field_status[qc_ok & is_q & !is_div] <- "Putative_Novel_Or_Hybrid_Query"
  ga_df$genomic_field_status[qc_ok & !is_q & is_adm] <- "Genomic_Admixed_Or_Introgressant"
  disc_mask <- qc_ok & !is_q & !is_adm & !is.na(sp_std) & nzchar(sp_std) & !is.na(g_pred) & nzchar(g_pred) & (sp_std != g_pred) & is.finite(g_prob) & (g_prob >= 0.70)
  ga_df$genomic_field_status[disc_mask] <- "Genomic_Discordant_Reassign_Candidate"
  conc_mask <- qc_ok & !is_q & !is_adm & !is.na(sp_std) & nzchar(sp_std) & !is.na(g_pred) & nzchar(g_pred) & (sp_std == g_pred)
  ga_df$genomic_field_status[conc_mask] <- "Genomic_Concordant_Anchor"

  out_p <- file.path(t_dir, "genomic_field_triage.csv")
  tmp_p <- paste0(out_p, ".tmp")
  readr::write_csv(ga_df, tmp_p, na = "")
  file.rename(tmp_p, out_p)
  invisible(out_p)
}

build_genomic_markdown_summary <- function(opts, df) {
  t_dir <- .resolve_step07_genomic_tables_dir(opts)
  r_dir <- .resolve_step07_genomic_reports_dir(opts)
  resolve_csv <- function(fname) {
    p1 <- file.path(r_dir, fname)
    if (file.exists(p1)) return(p1)
    p2 <- file.path(t_dir, fname)
    if (file.exists(p2)) return(p2)
    ""
  }
  ga_path <- resolve_csv("genomic_assignments.csv")
  gm_path <- resolve_csv("genomic_model_metrics.csv")
  gc_path <- resolve_csv("genomic_multimodal_congruence_summary.csv")
  if (!nzchar(ga_path) && !any(has_active_genomic_row_mask(df))) return(character(0))

  ga_df <- if (nzchar(ga_path)) tryCatch(utils::read.csv(ga_path, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  gm_df <- if (nzchar(gm_path)) tryCatch(utils::read.csv(gm_path, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  gc_df <- if (nzchar(gc_path)) tryCatch(utils::read.csv(gc_path, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  gdi_df <- attr(df, "genomic_gdi_delimitation")
  if (is.null(gdi_df) && nzchar(resolve_csv("coalescent_gdi_delimitation.csv"))) {
    gdi_df <- tryCatch(utils::read.csv(resolve_csv("coalescent_gdi_delimitation.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
  }
  fb_df <- attr(df, "genomic_fbranch_summary")
  if (is.null(fb_df) && nzchar(resolve_csv("introgressive_fbranch_summary.csv"))) {
    fb_df <- tryCatch(utils::read.csv(resolve_csv("introgressive_fbranch_summary.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
  }
  land_df <- attr(df, "genomic_landscape_summary")
  if (is.null(land_df) && nzchar(resolve_csv("landscape_mmrr_feems_summary.csv"))) {
    land_df <- tryCatch(utils::read.csv(resolve_csv("landscape_mmrr_feems_summary.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
  }

  lines <- c("", "## 5. Independent Step 06 Genomic & Phylogenomic Evidence Stream")
  if (!is.null(gm_df) && nrow(gm_df) > 0L) {
    lines <- c(
      lines,
      sprintf(
        "- **Cohort QC, HDplot Paralog Filter & Shrinkage DAPC (LORO-CV):** %d/%d sequenced libraries passed QC (`hq_reads >= 50000`, `missing_genotype_rate <= 0.80`; %d reference anchors, %d novel/hybrid queries, %d HDplot-clean loci retained; LORO-CV accuracy = %.3f, Cohen's kappa = %.3f, mean `M`/`P` replicate concordance = %.4f).",
        as.integer(gm_df$n_samples_pass[1L]),
        as.integer(gm_df$n_samples_total[1L]),
        as.integer(gm_df$n_reference_anchors[1L]),
        as.integer(gm_df$n_novel_queries[1L]),
        as.integer(gm_df$n_loci_retained[1L]),
        suppressWarnings(as.numeric(gm_df$loro_cv_accuracy[1L])),
        suppressWarnings(as.numeric(gm_df$loro_cv_kappa[1L])),
        suppressWarnings(as.numeric(gm_df$mean_replicate_concordance[1L]))
      )
    )
  }
  if (!is.null(ga_df) && nrow(ga_df) > 0L) {
    n_g12 <- sum(as.character(ga_df$linkage_tier) %in% c("G1_Direct_Voucher", "G2_Population_Voucher_Proxy"), na.rm = TRUE)
    n_g3 <- sum(as.character(ga_df$linkage_tier) == "G3_Georeferenced_Unvouchered", na.rm = TRUE)
    n_g4 <- sum(as.character(ga_df$linkage_tier) == "G4_Ungeoreferenced_Tip", na.rm = TRUE)
    lines <- c(
      lines,
      sprintf("- **Voucher Linkage Tiers:** %d voucher-linked (`G1`/`G2`), %d georeferenced field samples (`G3`), %d ungeoreferenced sequence tips (`G4`). Field and unvouchered query dispositions exported to `genomic_field_triage.csv`.", n_g12, n_g3, n_g4)
    )
    if ("inferred_cytotype" %in% names(ga_df)) {
      ct <- as.character(ga_df$inferred_cytotype[!is.na(ga_df$inferred_cytotype) & nzchar(ga_df$inferred_cytotype)])
      if (length(ct) > 0L) {
        ct_tab <- table(ct)
        ct_str <- paste(sprintf("%s=%d", names(ct_tab), as.integer(ct_tab)), collapse = ", ")
        lines <- c(lines, sprintf("- **Allele-Balance Cytotype & Hybrid Architecture:** Inferred ploidies (%s)%s.", ct_str, if ("hybrid_generation_class" %in% names(ga_df)) {
          hc <- as.character(ga_df$hybrid_generation_class[!is.na(ga_df$hybrid_generation_class) & nzchar(ga_df$hybrid_generation_class)])
          if (length(hc) > 0L) sprintf("; `(h, Q12)` classes: %s", paste(sprintf("%s=%d", names(table(hc)), as.integer(table(hc))), collapse = ", ")) else ""
        } else ""))
      }
    }
  }
  if (!is.null(gdi_df) && nrow(gdi_df) > 0L && "coalescent_gdi" %in% names(gdi_df)) {
    gdv <- suppressWarnings(as.numeric(gdi_df$coalescent_gdi))
    gdv <- gdv[is.finite(gdv)]
    if (length(gdv) > 0L) {
      lines <- c(
        lines,
        sprintf(
          "- **Coalescent Delimitation (`SVDquartets` / `BPP` `MSci` / `gdi`):** Across %d evaluated lineages, median `gdi` = %.3f (range %.3f to %.3f; Jackson et al. 2017: `<0.20` conspecific oversplit block, `>0.70` distinct lineage).",
          nrow(gdi_df), stats::median(gdv), min(gdv), max(gdv)
        )
      )
    }
  }
  if (!is.null(fb_df) && nrow(fb_df) > 0L && "f_branch" %in% names(fb_df)) {
    fmax <- suppressWarnings(max(as.numeric(fb_df$f_branch), na.rm = TRUE))
    if (is.finite(fmax)) {
      top_fb <- fb_df[which.max(suppressWarnings(as.numeric(fb_df$f_branch))), , drop = FALSE]
      donor_s <- if ("donor_taxon" %in% names(top_fb)) as.character(top_fb$donor_taxon[1L]) else "NA"
      lines <- c(lines, sprintf("- **Tree-Guided Introgression (`f-branch`):** Maximum `f_b` = %.4f (donor `%s`).", fmax, donor_s))
    }
  }
  if (!is.null(land_df) && nrow(land_df) > 0L && "mmrr_beta_ibd" %in% names(land_df)) {
    b_ibd <- suppressWarnings(as.numeric(land_df$mmrr_beta_ibd[1L]))
    b_ibe <- suppressWarnings(as.numeric(land_df$mmrr_beta_ibe[1L]))
    p_ibe <- suppressWarnings(as.numeric(land_df$mmrr_p_ibe[1L]))
    if (is.finite(b_ibd) || is.finite(b_ibe)) {
      lines <- c(
        lines,
        sprintf(
          "- **Landscape Genomics (`MMRR` / `conStruct` / `FEEMS`):** `beta_IBD` = %s, `beta_IBE` = %s (`p_IBE` = %s).",
          ifelse(is.finite(b_ibd), sprintf("%.3f", b_ibd), "NA"),
          ifelse(is.finite(b_ibe), sprintf("%.3f", b_ibe), "NA"),
          ifelse(is.finite(p_ibe), sprintf("%.4f", p_ibe), "NA")
        )
      )
    }
  }
  if (!is.null(gc_df) && nrow(gc_df) > 0L) {
    c_strs <- vapply(seq_len(nrow(gc_df)), function(k) {
      sprintf("%s (n=%d, Mantel r=%s)", gc_df$comparison[k], as.integer(gc_df$n_samples[k]), ifelse(is.finite(suppressWarnings(as.numeric(gc_df$mantel_r[k]))), sprintf("%.3f", as.numeric(gc_df$mantel_r[k])), "NA"))
    }, character(1))
    lines <- c(lines, paste0("- **Post-Hoc Multi-Modal Congruence:** ", paste(c_strs, collapse = "; "), "."))
  }
  lines
}

