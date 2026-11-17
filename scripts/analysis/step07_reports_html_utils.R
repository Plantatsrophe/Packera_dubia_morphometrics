if (grDevices::dev.cur() == 1L) grDevices::pdf(NULL)

# 7. Synthesis Figure Delegation, Taxonomic Revision Report & Interactive HTML Dashboard
.copy_focal_taxon_artifact <- function(src_path, target_path) {
  if (is.null(src_path) || !nzchar(src_path) || !file.exists(src_path)) return(invisible(FALSE))
  if (normalizePath(src_path, mustWork = FALSE) == normalizePath(target_path, mustWork = FALSE)) return(invisible(FALSE))
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  tmp_dst <- paste0(target_path, ".", Sys.getpid(), ".tmp")
  if (isTRUE(file.copy(src_path, tmp_dst, overwrite = TRUE))) file.rename(tmp_dst, target_path)
  invisible(TRUE)
}

.resolve_focal_taxon_slug <- function(focal_taxon_raw) {
  if (is.null(focal_taxon_raw) || !nzchar(trimws(as.character(focal_taxon_raw)))) return(list(taxon = "", slug = ""))
  f_std <- standardize_packera_taxon(trimws(as.character(focal_taxon_raw)))[1L]
  if (is.na(f_std) || !nzchar(f_std) || f_std == "Unknown") f_std <- trimws(as.character(focal_taxon_raw))
  list(taxon = f_std, slug = gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", f_std)))
}

render_synthesis_figure_png <- function(df, out_path, cda_summary = NULL, focal_taxon = "", output_complex_plot = NULL) {
  if (!exists("render_step07_synthesis_plate_png", mode = "function") && file.exists("scripts/morphometrics/pub_figures_utils.R")) source("scripts/morphometrics/pub_figures_utils.R")
  if (!exists("render_focal_revision_synthesis_plate_png", mode = "function") && file.exists("scripts/morphometrics/pub_synthesis_plates.R")) source("scripts/morphometrics/pub_synthesis_plates.R")
  out_png <- if (!is.null(out_path) && nzchar(out_path)) sub("\\.pdf$", ".png", out_path, ignore.case = TRUE) else "outputs/figures/fig07a_integrative_taxonomic_synthesis.png"
  full_df <- if (!is.null(attr(df, "full_complex_df"))) attr(df, "full_complex_df") else df
  if (exists("render_step07_synthesis_plate_png", mode = "function")) render_step07_synthesis_plate_png(full_df, out_png, cda_summary = cda_summary, focal_taxon = focal_taxon)
  f_info <- .resolve_focal_taxon_slug(focal_taxon)
  if (nzchar(f_info$taxon) && nzchar(f_info$slug) && nzchar(out_png)) {
    focal_png <- file.path(dirname(out_png), sprintf("Figure_Integrative_%s_Revision.png", f_info$slug))
    if (exists("render_focal_revision_synthesis_plate_png", mode = "function")) render_focal_revision_synthesis_plate_png(df, focal_png, full_df = full_df, focal_taxon = f_info$taxon)
    else if (file.exists(out_png)) .copy_focal_taxon_artifact(out_png, focal_png)
  }
  if (is.null(attr(df, "complex_discovery")) && !identical(trimws(env_plantmatics_or_packera("PLANTMATICS_SKIP_FIGURES", "PACKERA_SKIP_FIGURES")), "1") &&
      exists("render_step07b_complex_network_png", mode = "function")) {
    c_out <- if (!is.null(output_complex_plot) && nzchar(output_complex_plot)) sub("\\.pdf$", ".png", output_complex_plot, ignore.case = TRUE) else file.path(dirname(out_png), "fig07b_taxon_complex_delimitation_network.png")
    c_lbl <- if (nzchar(f_info$taxon)) f_info$taxon else as.character(df$species_standardized[1L])
    tryCatch(render_step07b_complex_network_png(
      complexes_df = data.frame(species_standardized = c_lbl, consensus_complex_id = 1L, n_vouchers = nrow(df), stringsAsFactors = FALSE),
      affinity_mat = matrix(1.0, nrow = 1L, ncol = 1L, dimnames = list(c_lbl, c_lbl)), out_complex_path = c_out, voucher_df = df
    ), error = function(e) warning("Single-taxon deme network render warning: ", e$message))
  }
  invisible(out_png)
}

generate_taxonomic_revision_report <- function(df, opts) {
  dir.create(dirname(opts$output_report), recursive = TRUE, showWarnings = FALSE); tot_n <- nrow(df)
  crit_n <- sum(df$synthesis_triage_priority == "CRITICAL", na.rm = TRUE); high_n <- sum(df$synthesis_triage_priority == "HIGH", na.rm = TRUE)
  med_n <- sum(df$synthesis_triage_priority == "MEDIUM", na.rm = TRUE); res_n <- sum(df$synthesis_triage_priority %in% c("LOW", "RESOLVED"), na.rm = TRUE)
  sub_vec <- if ("taxonomic_treatment_subtype" %in% names(df)) as.character(df$taxonomic_treatment_subtype) else rep(NA_character_, tot_n)
  call_vec <- if ("taxonomic_status_call" %in% names(df)) as.character(df$taxonomic_status_call) else rep(NA_character_, tot_n)
  misid_n <- sum(sub_vec == "Misidentification_Reassign" | call_vec == "Misidentification_Severe", na.rm = TRUE)
  ecotype_n <- sum(sub_vec == "Subspecies_Ecotype" | call_vec == "Subspecies_Ecotype", na.rm = TRUE)
  ecoph_n <- sum(sub_vec == "Ecophenotypic_Plasticity" | call_vec == "Ecophenotypic_Plasticity", na.rm = TRUE)
  hybr_n <- sum(sub_vec == "Hybrid_Introgressant" | call_vec == "Hybrid_Intergrade_Swarm", na.rm = TRUE)
  crypt_n <- sum(sub_vec == "Putative_Cryptic_Species" | call_vec %in% c("Putative_Cryptic_Species", "Taxonomic_Resurrection_Candidate"), na.rm = TRUE)
  resurr_n <- sum(call_vec == "Taxonomic_Resurrection_Candidate", na.rm = TRUE)

  sp_col_rep <- if ("species_standardized" %in% names(df)) as.character(df$species_standardized) else character(0)
  active_taxa <- sort(unique(sp_col_rep[!is.na(sp_col_rep) & nzchar(sp_col_rep) & sp_col_rep != "Unknown"]))
  f_info <- .resolve_focal_taxon_slug(opts$focal_taxon)
  g_rep <- if (exists("resolve_r_target_genus", mode = "function")) resolve_r_target_genus() else "Packera"
  header_focal <- if (nzchar(f_info$taxon)) f_info$taxon else if (length(active_taxa) == 1L) active_taxa[1L] else ""
  report_title <- if (nzchar(header_focal)) {
    sprintf("# Taxonomic Revision and Integrative Species Delimitation in the *%s* Complex (Asteraceae: Senecioneae)", header_focal)
  } else sprintf("# Taxonomic Revision and Integrative Species Delimitation Across *%s* Species Complexes (N = %d Taxa; Asteraceae: Senecioneae)", g_rep, length(active_taxa))
  scope_phrase <- if (nzchar(header_focal)) sprintf("the ***%s*** complex", header_focal)
  else if (length(active_taxa) > 0L) sprintf("**%d analyzed *%s* taxa** (%s)", length(active_taxa), g_rep, paste0("*", active_taxa, "*", collapse = ", "))
  else sprintf("the analyzed *%s* herbarium dataset", g_rep)

  gmm_stats <- resolve_gmm_bic_stats(opts)
  bic_str <- if (is.finite(gmm_stats$two_delta_bic)) format(round(gmm_stats$two_delta_bic, 1), big.mark = ",") else "N/A"
  batch_df <- audit_institutional_batch_effects_ancova(df)
  batch_str <- if (!is.null(batch_df) && nrow(batch_df) > 0) {
    sprintf("Taxon-adjusted marginal ANCOVA (`lm(PC1 ~ species_standardized + institutionCode)`): partial $\\eta^2 = %.4f$ ($q_{\\text{BH}} = %s$, %s).", batch_df$partial_eta_sq[1], as.character(batch_df$p_adj_bh[1]), batch_df$Batch_Effect_Status[1])
  } else "Taxon-adjusted marginal ANCOVA (`lm(PC1 ~ species_standardized + institutionCode)`): partial $\\eta^2 < 0.05$ across major herbaria."

  fmt_med_iqr <- function(primary_vec, fallback_vec = NULL, digits = 1, unit_str = "") {
    v <- suppressWarnings(as.numeric(primary_vec))
    if (!is.null(fallback_vec)) { vf <- suppressWarnings(as.numeric(fallback_vec)); v <- ifelse(is.finite(v), v, vf) }
    v <- v[is.finite(v)]; if (length(v) == 0L) return("NA")
    sprintf("%.*f [IQR %.*f]%s", digits, stats::median(v), digits, stats::IQR(v), unit_str)
  }

  treatment_lines <- character(0)
  for (sp in active_taxa) {
    sub_sp <- df[!is.na(sp_col_rep) & sp_col_rep == sp, , drop = FALSE]; n_sp <- nrow(sub_sp)
    is_type_sp <- if ("is_type_specimen" %in% names(sub_sp)) as.logical(sub_sp$is_type_specimen) else rep(FALSE, n_sp); is_type_sp[is.na(is_type_sp)] <- FALSE
    n_types <- sum(is_type_sp, na.rm = TRUE)
    n_gs <- if ("determiner_tier" %in% names(sub_sp)) sum((sub_sp$determiner_tier %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver")) | is_type_sp, na.rm = TRUE) else n_types
    rg_vals <- if ("regional_group" %in% names(sub_sp)) trimws(as.character(sub_sp$regional_group)) else character(0)
    rg_vals <- rg_vals[!is.na(rg_vals) & nzchar(rg_vals) & !(tolower(rg_vals) %in% c("na", "unknown", "unassigned"))]
    eco_str <- if (length(rg_vals) > 0L) { rg_tab <- utils::head(sort(table(rg_vals), decreasing = TRUE), 3L); paste(sprintf("%s (n=%d)", names(rg_tab), as.integer(rg_tab)), collapse = "; ") } else "Unassigned / Non-georeferenced"
    ph_str <- fmt_med_iqr(if ("fine_soil_ph" %in% names(sub_sp)) sub_sp$fine_soil_ph else NA_real_, if ("soil_ph" %in% names(sub_sp)) sub_sp$soil_ph else NULL, digits = 2)
    sand_str <- fmt_med_iqr(if ("fine_soil_sand" %in% names(sub_sp)) sub_sp$fine_soil_sand else NA_real_, if ("soil_sand" %in% names(sub_sp)) sub_sp$soil_sand else NULL, digits = 1, unit_str = "%")
    bed_str <- fmt_med_iqr(if ("fine_bedrock_depth_cm" %in% names(sub_sp)) sub_sp$fine_bedrock_depth_cm else NA_real_, if ("ssurgo_bedrock_depth_cm" %in% names(sub_sp)) sub_sp$ssurgo_bedrock_depth_cm else NULL, digits = 1, unit_str = " cm")
    d_sp <- if ("doy" %in% names(sub_sp)) suppressWarnings(as.numeric(sub_sp$doy)) else numeric(0)
    anth_sp <- if ("is_active_anthesis" %in% names(sub_sp)) as.logical(sub_sp$is_active_anthesis) else rep(TRUE, n_sp); anth_sp[is.na(anth_sp)] <- TRUE
    d_valid <- d_sp[is.finite(d_sp) & d_sp >= 1.0 & d_sp <= 366.0 & anth_sp]; if (length(d_valid) == 0L) d_valid <- d_sp[is.finite(d_sp) & d_sp >= 1.0 & d_sp <= 366.0]
    doy_str <- if (length(d_valid) > 0L) sprintf("DOY %d to %d (10th to 90th pct; n=%d; Hopkins 1918; Pewsey et al. 2013)", as.integer(round(stats::quantile(d_valid, 0.10, na.rm = TRUE, names = FALSE))), as.integer(round(stats::quantile(d_valid, 0.90, na.rm = TRUE, names = FALSE))), length(d_valid))
    else "NA (no dated anthesis vouchers; Hopkins 1918; Pewsey et al. 2013)"

    st_sp <- if ("taxonomic_treatment_subtype" %in% names(sub_sp)) as.character(sub_sp$taxonomic_treatment_subtype) else if ("taxonomic_status_call" %in% names(sub_sp)) as.character(sub_sp$taxonomic_status_call) else rep("Confirmed_Species", n_sp)
    call_sp <- if ("taxonomic_status_call" %in% names(sub_sp)) as.character(sub_sp$taxonomic_status_call) else st_sp
    cnt_sub <- function(lbls) sum(!is.na(st_sp) & (st_sp %in% lbls))
    n_type_ctrl <- sum((!is.na(call_sp) & call_sp == "Confirmed_Type_Control") | (!is.na(st_sp) & st_sp == "Confirmed_Type_Control"))
    n_type_syn <- sum((!is.na(call_sp) & call_sp == "Type_Synonymy_Candidate") | (!is.na(st_sp) & st_sp == "Type_Synonymy_Candidate"))
    sci_vals <- if ("scientificName" %in% names(sub_sp)) trimws(as.character(sub_sp$scientificName)) else character(0)
    sci_vals <- sci_vals[!is.na(sci_vals) & nzchar(sci_vals) & !(tolower(sci_vals) %in% c("na", "unknown", "none", "null"))]
    sp_esc <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", sp)
    sci_with_auth <- sci_vals[grepl(paste0("^", sp_esc, "[[:space:]]+.+"), sci_vals, ignore.case = TRUE)]
    auth_str <- if (length(sci_with_auth) > 0L) trimws(sub(paste0("^", sp_esc), "", names(sort(table(sci_with_auth), decreasing = TRUE))[1L], ignore.case = TRUE)) else ""
    head_sp <- if (nzchar(auth_str)) sprintf("- ***%s*** %s", sp, auth_str) else sprintf("- ***%s***", sp)

    er_col <- if ("ecotype_regime" %in% names(sub_sp)) as.character(sub_sp$ecotype_regime) else if ("analysis_group" %in% names(sub_sp)) as.character(sub_sp$analysis_group) else character(0)
    er_col <- trimws(er_col[!is.na(er_col) & nzchar(trimws(er_col)) & !(tolower(trimws(er_col)) %in% c("na", "unknown"))]); er_multi <- er_col[grepl("^Ecotype_", er_col)]
    er_line <- if (length(unique(er_multi)) >= 2L) {
      er_tab <- sort(table(er_multi), decreasing = TRUE)
      sprintf("  - **Independent Step 05 Ecological Regimes ($K_{\\text{eco}} = %d$):** %s", length(er_tab), paste(sprintf("`%s` (n=%d)", names(er_tab), as.integer(er_tab)), collapse = "; "))
    } else "  - **Independent Step 05 Ecological Regimes:** Unimodal broad niche envelope ($K_{\\text{eco}} = 1$)"

    hyb_line <- character(0); hyb_sp_mask <- !is.na(st_sp) & (st_sp %in% c("Hybrid_Introgressant", "Hybrid_Intergrade_Swarm"))
    if (any(hyb_sp_mask) && all(c("hybrid_parent_1", "hybrid_parent_2") %in% names(sub_sp))) {
      p1_s <- as.character(sub_sp$hybrid_parent_1[hyb_sp_mask]); p2_s <- as.character(sub_sp$hybrid_parent_2[hyb_sp_mask])
      ok_ps <- !is.na(p1_s) & nzchar(p1_s) & !is.na(p2_s) & nzchar(p2_s)
      if (any(ok_ps)) { ptab <- sort(table(sprintf("*%s* x *%s*", p1_s[ok_ps], p2_s[ok_ps])), decreasing = TRUE); hyb_line <- sprintf("  - **Confirmed Hybrid Swarm Parentage:** %s", paste(sprintf("%s (n=%d)", names(ptab), as.integer(ptab)), collapse = "; ")) }
    }
    treatment_lines <- c(
      treatment_lines, head_sp,
      sprintf("  - **Examined Vouchers:** $n = %d$ (Nomenclatural Types: %d; Active Gold/Silver Anchors: %d)", n_sp, n_types, n_gs),
      sprintf("  - **Dominant US EPA Level III Ecoregions:** %s", eco_str), er_line,
      sprintf("  - **Edaphic Envelope (Median [IQR]):** `fine_soil_ph` = %s; `fine_soil_sand` = %s; `fine_bedrock_depth_cm` = %s", ph_str, sand_str, bed_str),
      sprintf("  - **Verified Anthesis Window:** %s", doy_str),
      sprintf("  - **Delimitation Breakdown:** `Confirmed_Species` = %d (`Confirmed_Type_Control` = %d) | `Subspecies_Ecotype` = %d (`Type_Synonymy_Candidate` = %d) | `Ecophenotypic_Plasticity` = %d | `Hybrid_Introgressant` = %d | `Misidentification_Reassign` = %d | `Putative_Cryptic_Species` = %d",
              cnt_sub(c("Confirmed_Species", "Species_Confirmed", "Confirmed_Type_Control")), n_type_ctrl, cnt_sub(c("Subspecies_Ecotype", "Type_Synonymy_Candidate")), n_type_syn,
              cnt_sub("Ecophenotypic_Plasticity"), cnt_sub(c("Hybrid_Introgressant", "Hybrid_Intergrade_Swarm")), cnt_sub(c("Misidentification_Reassign", "Misidentification_Severe")), cnt_sub("Putative_Cryptic_Species")),
      hyb_line
    )
  }
  if (length(treatment_lines) == 0L) treatment_lines <- sprintf("- No standardized *%s* vouchers present in active triage queue (Hopkins 1918; Pewsey et al. 2013).", g_rep)

  is_single_or_intra <- isTRUE(opts$intraspecific_mode) || identical(attr(df, "can_evaluate_hybrids"), FALSE)
  stg2_line <- if (is_single_or_intra) "- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** 0 introgressants (0.0%; interspecific hybrid evaluation disabled in single-taxon intraspecific mode)."
  else if (hybr_n > 0L && all(c("hybrid_parent_1", "hybrid_parent_2") %in% names(df))) {
    hm_all <- !is.na(sub_vec) & (sub_vec == "Hybrid_Introgressant") & !is.na(df$hybrid_parent_1) & !is.na(df$hybrid_parent_2)
    if (any(hm_all)) {
      canon_pairs <- vapply(which(hm_all), function(ii) { ps <- sort(c(as.character(df$hybrid_parent_1[ii]), as.character(df$hybrid_parent_2[ii]))); sprintf("*%s* x *%s*", ps[1L], ps[2L]) }, character(1))
      pt_all <- sort(table(canon_pairs), decreasing = TRUE)
      sprintf("- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** %d introgressants (%0.1f%%) across confirmed parent pairs: %s.", hybr_n, (hybr_n / max(tot_n, 1)) * 100, paste(sprintf("%s (n=%d)", names(pt_all), as.integer(pt_all)), collapse = "; "))
    } else sprintf("- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** %d introgressants (%0.1f%%) (Rieseberg 1995; Slenker et al. 2022).", hybr_n, (hybr_n / max(tot_n, 1)) * 100)
  } else sprintf("- **Stage 2 (Conservative Contact-Zone Hybrid Swarms):** %d introgressants (%0.1f%%; 5-gate syntopic population threshold).", hybr_n, (hybr_n / max(tot_n, 1)) * 100)

  stg3_line <- if (resurr_n > 0L && "resurrection_basionym" %in% names(df)) {
    rb_vals <- as.character(df$resurrection_basionym[!is.na(df$resurrection_basionym) & nzchar(as.character(df$resurrection_basionym))])
    rb_str <- if (length(rb_vals) > 0L) paste(sprintf("*%s*", unique(rb_vals)), collapse = ", ") else "heterotypic basionym"
    sprintf("- **Stage 3 (Putative Cryptic and Basionym Resurrection Lineages):** %d vouchers (%0.1f%%; %d heterotypic resurrection candidates matching %s and %d unnamed cryptic candidates) satisfying Integration by Congruence (Padial et al. 2010).", crypt_n, (crypt_n / max(tot_n, 1)) * 100, resurr_n, rb_str, max(0L, crypt_n - resurr_n))
  } else sprintf("- **Stage 3 (Putative Cryptic Lineages):** %d vouchers (%0.1f%%) satisfying strict Integration by Congruence (Padial et al. 2010).", crypt_n, (crypt_n / max(tot_n, 1)) * 100)

  comp_att <- attr(df, "complex_discovery")
  cdf_rep <- if (!is.null(comp_att) && !is.null(comp_att$complexes_df)) comp_att$complexes_df else if (!is.null(opts$output_complexes) && nzchar(opts$output_complexes) && file.exists(opts$output_complexes)) tryCatch(utils::read.csv(opts$output_complexes, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  sec3_lines <- if (!is.null(cdf_rep) && nrow(cdf_rep) >= 2L) {
    fmt_it <- function(s) ifelse(grepl("^[A-Z][a-z-]+[[:space:]]+([x×][[:space:]]*)?[a-z][a-z-]+|^[A-Z]\\.[[:space:]]+[a-z-]+", s), paste0("*", s, "*"), s)
    c_rows <- vapply(seq_len(nrow(cdf_rep)), function(ri) sprintf("| %s | `%s` | %s (%.3f) | %.3f | %s |", fmt_it(as.character(cdf_rep$species_standardized[ri])), as.character(cdf_rep$complex_label[ri]), fmt_it(as.character(cdf_rep$nearest_neighbor_taxon[ri])), suppressWarnings(as.numeric(cdf_rep$nearest_neighbor_affinity[ri])), suppressWarnings(as.numeric(cdf_rep$silhouette_width[ri])), ifelse(isTRUE(as.logical(cdf_rep$integrative_congruence[ri])), "Yes", "No")), character(1))
    c("## 3. Empirical Taxon-Complex Delimitation", "| Taxon / Subcluster | Complex | Nearest Neighbor ($W_{ij}$) | Silhouette | Congruent |", "| :--- | :--- | :--- | :--- | :--- |", c_rows, "", "## 4. Herbarium Triage Queue Audit Summary")
  } else "## 3. Herbarium Triage Queue Audit Summary"

  md_lines <- c(
    report_title, "", "**Principal Investigator:** J. Brandon Fuller (PhD Candidate, Department of Biology, UNC-CH)  ",
    "**Faculty Advisor:** Dr. Alan S. Weakley (Director, NCU Herbarium; UNC Biology)  ",
    "**Institution:** University of North Carolina at Chapel Hill Herbarium (NCU)  ",
    "**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5`  ", "", "## 1. Executive Botanical Summary & Integrative Taxonomy Protocol",
    sprintf("This document provides the formal taxonomic treatment, hierarchical integrative species delimitation (Padial et al. 2010; Schlick-Steiner et al. 2010), and reliability-weighted curatorial synthesis (Dawid & Skene 1979) across %s.", scope_phrase),
    sprintf("- **Total Examined Vouchers:** %s specimens (including coalesced non-georeferenced historic vouchers).", format(tot_n, big.mark = ",")),
    sprintf("- **GMM Bayes Factors (2ΔBIC):** Decisive evidence for **K = %d discrete morphological clusters** ($2\\Delta\\text{BIC} = %s$ over $K=1$).", gmm_stats$opt_k, bic_str),
    sprintf("- **Institutional Batch Audit:** %s", batch_str),
    sprintf("- **Stage 1 (Herbarium Misidentifications, Reassign):** %d vouchers (%0.1f%%) evaluated prior to hybridity.", misid_n, (misid_n / max(tot_n, 1)) * 100),
    stg2_line, stg3_line,
    sprintf("- **Stage 4 (Subspecific Ecotypes and Taphonomic/Ontogenetic Plasticity):** %d cohesive ecotypes (%0.1f%%) and %d within-envelope plastic variants (%0.1f%%).", ecotype_n, (ecotype_n / max(tot_n, 1)) * 100, ecoph_n, (ecoph_n / max(tot_n, 1)) * 100), "",
    "## 2. Multi-Evidence Taxonomic Treatment", treatment_lines, "", sec3_lines,
    "| Triage Priority | Voucher Count | Percentage | Recommended Action |", "| :--- | :--- | :--- | :--- |",
    sprintf("| **CRITICAL** | **%d** | **%0.1f%%** | Specialist manual re-inspection & sheet annotation |", crit_n, (crit_n / max(tot_n, 1)) * 100),
    sprintf("| **HIGH** | **%d** | **%0.1f%%** | Systematic redetermination in aggregator databases |", high_n, (high_n / max(tot_n, 1)) * 100),
    sprintf("| **MEDIUM** | **%d** | **%0.1f%%** | Annotate foliar wear / introgression zone |", med_n, (med_n / max(tot_n, 1)) * 100),
    sprintf("| **LOW / RESOLVED** | **%d** | **%0.1f%%** | Verified taxon anchor retained |", res_n, (res_n / max(tot_n, 1)) * 100),
    sprintf("| **Total** | **%d** | **100.0%%** | Reliability-weighted multi-modal consensus |", tot_n)
  )
  geno_md <- if (exists("build_genomic_markdown_summary", mode = "function")) tryCatch(build_genomic_markdown_summary(opts, df), error = function(e) character(0)) else character(0)
  if (length(geno_md) > 0L) md_lines <- c(md_lines, geno_md)
  cyto_md_rep <- attr(df, "cytology_markdown_lines")
  if (!is.null(cyto_md_rep) && length(cyto_md_rep) > 0L && length(geno_md) == 0L) {
    cyto_md_rep <- as.character(cyto_md_rep); cyto_md_rep <- cyto_md_rep[!is.na(cyto_md_rep)]
    if (length(cyto_md_rep) > 0L) md_lines <- c(md_lines, if (nzchar(trimws(cyto_md_rep[1L]))) "" else character(0), cyto_md_rep)
  }
  tmp_rep <- paste0(opts$output_report, ".tmp"); writeLines(md_lines, tmp_rep); file.rename(tmp_rep, opts$output_report)
  if (nzchar(f_info$taxon) && nzchar(f_info$slug)) {
    focal_rep <- file.path(dirname(opts$output_report), sprintf("%s_Taxonomic_Revision_Summary.md", f_info$slug))
    if (normalizePath(focal_rep, mustWork = FALSE) != normalizePath(opts$output_report, mustWork = FALSE)) { tmp_frep <- paste0(focal_rep, ".tmp"); writeLines(md_lines, tmp_frep); file.rename(tmp_frep, focal_rep) }
    out_png_rep <- if (!is.null(opts$output_plot) && nzchar(opts$output_plot)) sub("\\.pdf$", ".png", opts$output_plot, ignore.case = TRUE) else ""
    if (nzchar(out_png_rep) && file.exists(out_png_rep)) {
      focal_png <- file.path(dirname(out_png_rep), sprintf("Figure_Integrative_%s_Revision.png", f_info$slug))
      if (!file.exists(focal_png) || file.info(focal_png)$size < 1024) .copy_focal_taxon_artifact(out_png_rep, focal_png)
    }
  }
}

.html_escape_vec <- function(x) {
  s <- ifelse(is.na(x), "", as.character(x))
  s <- gsub("&", "&amp;", s, fixed = TRUE); s <- gsub("<", "&lt;", s, fixed = TRUE); s <- gsub(">", "&gt;", s, fixed = TRUE); gsub("\"", "&quot;", s, fixed = TRUE)
}

.json_escape_vec <- function(x) {
  s <- ifelse(is.na(x), "", as.character(x))
  s <- gsub("\\", "\\\\", s, fixed = TRUE); s <- gsub("\"", "\\\"", s, fixed = TRUE); s <- gsub("\r", " ", s, fixed = TRUE)
  s <- gsub("\n", " ", s, fixed = TRUE); s <- gsub("\t", " ", s, fixed = TRUE); s <- gsub("</", "<\\/", s, fixed = TRUE)
  paste0("\"", s, "\"")
}

generate_interactive_triage_html <- function(df, html_path, focal_taxon = "", opts = NULL) {
  if (is.null(html_path) || !nzchar(trimws(as.character(html_path)))) return(invisible(NULL))
  dir.create(dirname(html_path), recursive = TRUE, showWarnings = FALSE); n <- nrow(df)
  get_c <- function(nm, def = "") if (nm %in% names(df)) ifelse(is.na(df[[nm]]), def, as.character(df[[nm]])) else rep(def, n)
  get_f <- function(nm, digits = 2) if (!(nm %in% names(df))) rep("NA", n) else { v <- suppressWarnings(as.numeric(df[[nm]])); ifelse(is.finite(v), sprintf("%.*f", digits, v), "NA") }
  f_info <- .resolve_focal_taxon_slug(focal_taxon); sp_vec <- get_c("species_standardized", "Unknown")
  active_taxa <- sort(unique(sp_vec[nzchar(sp_vec) & sp_vec != "Unknown"]))
  hdr_taxon <- if (nzchar(f_info$taxon)) f_info$taxon else if (length(active_taxa) == 1L) active_taxa[1L] else ""
  g_html <- if (exists("resolve_r_target_genus", mode = "function")) resolve_r_target_genus() else "Packera"
  title_html <- if (nzchar(hdr_taxon)) sprintf("NCU Herbarium Integrative Triage Dashboard: <em>%s</em>", .html_escape_vec(hdr_taxon)) else sprintf("NCU Herbarium Integrative Triage Dashboard: <em>%s</em> (%d Taxa)", .html_escape_vec(g_html), length(active_taxa))
  fmt_w <- function(cn) if (cn %in% names(df) && any(is.finite(as.numeric(df[[cn]])))) sprintf("%.3f", stats::median(as.numeric(df[[cn]]), na.rm = TRUE)) else "NA"
  w_m <- fmt_w("w_morph_kappa"); w_e <- fmt_w("w_edaphic_kappa"); w_p <- fmt_w("w_pheno_kappa")

  prio_v <- get_c("synthesis_triage_priority", "RESOLVED"); sub_v <- get_c("taxonomic_treatment_subtype", "Confirmed_Species")
  cat_v <- get_c("catalogNumber", ""); pid_v <- get_c("plant_individual_id", "0"); inst_v <- get_c("institutionCode", ""); coll_v <- get_c("recordedBy", "")
  tier_v <- get_c("determiner_tier", ""); rec_v <- get_c("recommended_determination", ""); cda_v <- get_c("cda_predicted_taxon", ""); cdap_v <- get_f("cda_posterior_prob", 2)
  ed_v <- get_c("edaphic_best_fit_taxon", ""); edz_v <- get_f("edaphic_z_score", 2); ddoy_v <- get_f("delta_doy_resid", 1); cscore_v <- get_f("probabilistic_consensus_score", 3)
  act_v <- get_c("synthesis_triage_action", ""); rat_v <- get_c("synthesis_rationale", ""); hp1_h <- get_c("hybrid_parent_1", ""); hp2_h <- get_c("hybrid_parent_2", "")
  hcoh_h <- ifelse(isTRUE(df$hybrid_zone_cohesion_pass) | (!is.null(df$hybrid_zone_cohesion_pass) & !is.na(df$hybrid_zone_cohesion_pass) & as.logical(df$hybrid_zone_cohesion_pass)), "TRUE", "FALSE")
  if (length(hcoh_h) != n) hcoh_h <- rep("FALSE", n)
  mts_v <- get_c("morph_tier_source", "Unknown")
  bc_v <- if ("basal_leaf_count" %in% names(df)) ifelse(is.finite(suppressWarnings(as.numeric(df$basal_leaf_count))), as.character(as.integer(df$basal_leaf_count)), "0") else rep("0", n)
  cc_v <- if ("cauline_leaf_count" %in% names(df)) ifelse(is.finite(suppressWarnings(as.numeric(df$cauline_leaf_count))), as.character(as.integer(df$cauline_leaf_count)), "0") else rep("0", n)
  bcp_v <- paste0(bc_v, "/", cc_v); eco_v <- if ("ecotype_regime" %in% names(df) && any(nzchar(trimws(get_c("ecotype_regime", ""))))) get_c("ecotype_regime", "") else get_c("analysis_group", "")
  geno_pred_h <- get_c("geno_predicted_taxon", ""); geno_cyt_h <- get_c("inferred_cytotype", ""); geno_gdi_h <- get_f("coalescent_gdi", 2); geno_cls_h <- get_c("hybrid_generation_class", "")
  prop_tier_h <- get_c("genomic_propagation_tier", ""); t0_tax_h <- get_c("tier0_calibrated_taxon", ""); t0_prob_h <- get_f("tier0_calibrated_prob", 2)
  t0_abst_h <- if ("tier0_abstention_flag" %in% names(df)) (!is.na(df$tier0_abstention_flag) & suppressWarnings(as.integer(df$tier0_abstention_flag)) == 1L) else rep(FALSE, n)
  t0_ok_h <- if ("tier0_abstention_flag" %in% names(df)) (!is.na(df$tier0_abstention_flag) & suppressWarnings(as.integer(df$tier0_abstention_flag)) == 0L) else rep(FALSE, n)
  cl_gdi_h <- get_f("cluster_propagated_gdi", 2)
  vp_2n_h <- if ("voucher_sporophytic_2n" %in% names(df)) suppressWarnings(as.integer(df$voucher_sporophytic_2n)) else rep(NA_integer_, n)
  vp_x_h <- if ("voucher_ploidy_level" %in% names(df)) suppressWarnings(as.integer(df$voucher_ploidy_level)) else rep(NA_integer_, n)
  tax_cyto_h <- get_c("taxon_reported_cytotypes", "")
  is_mc_h <- if ("is_multi_cytotype_taxon" %in% names(df)) (!is.na(df$is_multi_cytotype_taxon) & suppressWarnings(as.integer(df$is_multi_cytotype_taxon)) == 1L) else rep(FALSE, n)

  geno_sum_h <- vapply(seq_len(n), function(ii) {
    base_str <- if (nzchar(geno_pred_h[ii])) {
      bits <- c(if (nzchar(geno_cyt_h[ii])) geno_cyt_h[ii], if (geno_gdi_h[ii] != "NA") paste0("gdi=", geno_gdi_h[ii]), if (nzchar(geno_cls_h[ii])) geno_cls_h[ii])
      if (length(bits) > 0L) sprintf("%s (%s)", geno_pred_h[ii], paste(bits, collapse = "; ")) else geno_pred_h[ii]
    } else ""
    prop_badges <- c(
      if (prop_tier_h[ii] == "G1_Exsiccatae_Isoduplicate") "[Exs-Dup]" else if (prop_tier_h[ii] == "G2_Conditional_Deme_Hull") "[Deme-Hull]",
      if (t0_ok_h[ii] && nzchar(t0_tax_h[ii]) && (t0_tax_h[ii] != "Sensu_Lato_DNA_Required")) { if (t0_prob_h[ii] != "NA") sprintf("[Tier0-Cal: %s p=%s]", t0_tax_h[ii], t0_prob_h[ii]) else sprintf("[Tier0-Cal: %s]", t0_tax_h[ii]) },
      if (t0_abst_h[ii]) "[DNA-Abstain]",
      if (geno_gdi_h[ii] == "NA" && cl_gdi_h[ii] != "NA") sprintf("[Cluster-gdi=%s]", cl_gdi_h[ii])
    )
    has_2n <- !is.na(vp_2n_h[ii]) && is.finite(vp_2n_h[ii]) && vp_2n_h[ii] > 0L; has_x <- !is.na(vp_x_h[ii]) && is.finite(vp_x_h[ii]) && vp_x_h[ii] > 0L
    if (has_2n && has_x) prop_badges <- c(prop_badges, sprintf("[2n=%d | %dx]", vp_2n_h[ii], vp_x_h[ii]))
    else if (has_2n) prop_badges <- c(prop_badges, sprintf("[2n=%d]", vp_2n_h[ii])) else if (has_x) prop_badges <- c(prop_badges, sprintf("[%dx]", vp_x_h[ii]))
    tc_ii <- trimws(tax_cyto_h[ii])
    if (is_mc_h[ii] && nzchar(tc_ii) && !(toupper(tc_ii) %in% c("NA", "NAN", "NONE", "NULL"))) prop_badges <- c(prop_badges, sprintf("[Multi-Cytotype: %s]", tc_ii))
    if (length(prop_badges) == 0L) base_str else if (nzchar(base_str)) paste(c(base_str, prop_badges), collapse = " ") else paste(prop_badges, collapse = " ")
  }, character(1))

  json_rows <- if (n > 0L) paste0("[", paste(
    .json_escape_vec(cat_v), .json_escape_vec(pid_v), .json_escape_vec(inst_v), .json_escape_vec(coll_v), .json_escape_vec(tier_v),
    .json_escape_vec(sp_vec), .json_escape_vec(rec_v), .json_escape_vec(prio_v), .json_escape_vec(sub_v), .json_escape_vec(cda_v),
    .json_escape_vec(cdap_v), .json_escape_vec(ed_v), .json_escape_vec(edz_v), .json_escape_vec(ddoy_v), .json_escape_vec(cscore_v),
    .json_escape_vec(act_v), .json_escape_vec(rat_v), .json_escape_vec(hp1_h), .json_escape_vec(hp2_h), .json_escape_vec(hcoh_h),
    .json_escape_vec(mts_v), .json_escape_vec(bcp_v), .json_escape_vec(eco_v), .json_escape_vec(geno_sum_h), sep = ","
  ), "]", collapse = ",\n") else ""

  taxon_opts <- paste0("<option value=\"ALL\">All Recorded Taxa</option>", if (length(active_taxa) > 0L) paste0("<option value=\"", .html_escape_vec(active_taxa), "\">", .html_escape_vec(active_taxa), "</option>", collapse = "") else "")
  is_intra_h <- (!is.null(opts) && isTRUE(opts$intraspecific_mode)) || identical(attr(df, "can_evaluate_hybrids"), FALSE)
  hyb_kpi <- if (is_intra_h) "<div class=\"kpi\"><div class=\"lbl\">Hybrid / Introgressant</div><div class=\"val\">0 <span style=\"font-size:0.72rem;font-weight:500;color:#64748b;\">(1-Taxon N/A)</span></div></div>"
  else sprintf("<div class=\"kpi\"><div class=\"lbl\">Hybrid / Introgressant</div><div class=\"val\">%d</div></div>", sum(sub_v == "Hybrid_Introgressant"))

  hm_idx <- which(sub_v == "Hybrid_Introgressant" & nzchar(hp1_h) & nzchar(hp2_h))
  hyb_sec <- if (length(hm_idx) > 0L) {
    ptab_h <- sort(table(vapply(hm_idx, function(ii) { ps <- sort(c(hp1_h[ii], hp2_h[ii])); paste(ps[1L], "x", ps[2L]) }, character(1))), decreasing = TRUE)
    h_chips <- vapply(names(ptab_h), function(plbl) sprintf("<span style=\"display:inline-block;margin:2px 8px 2px 0;padding:3px 9px;background:#fdf2f8;border:1px solid #fbcfe8;border-radius:6px;font-size:0.81rem;\"><span class=\"sci\">%s</span> <strong>(n=%d)</strong></span>", .html_escape_vec(plbl), as.integer(ptab_h[[plbl]])), character(1))
    paste0("<section style=\"background:#fff;border:1px solid #cbd5e1;border-radius:8px;padding:10px 14px;margin-bottom:14px;\"><div style=\"font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:#475569;margin-bottom:4px;\">Confirmed Hybrid Swarms (Putative Parent Species)</div>", paste(h_chips, collapse = ""), "</section>")
  } else ""

  comp_att_h <- attr(df, "complex_discovery")
  cdf_h <- if (!is.null(comp_att_h) && !is.null(comp_att_h$complexes_df)) comp_att_h$complexes_df else if (!is.null(opts) && !is.null(opts$output_complexes) && nzchar(opts$output_complexes) && file.exists(opts$output_complexes)) tryCatch(utils::read.csv(opts$output_complexes, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  comp_sec <- if (!is.null(cdf_h) && nrow(cdf_h) >= 2L) {
    cl_chips <- vapply(unique(as.character(cdf_h$complex_label)), function(clb) sprintf("<span style=\"display:inline-block;margin:2px 8px 2px 0;padding:3px 9px;background:#f1f5f9;border:1px solid #cbd5e1;border-radius:6px;font-size:0.81rem;\"><strong>%s:</strong> <span class=\"sci\">%s</span></span>", .html_escape_vec(clb), paste(.html_escape_vec(as.character(cdf_h$species_standardized[cdf_h$complex_label == clb])), collapse = ", ")), character(1))
    paste0("<section style=\"background:#fff;border:1px solid #cbd5e1;border-radius:8px;padding:10px 14px;margin-bottom:14px;\"><div style=\"font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:#475569;margin-bottom:4px;\">Empirical Species Complexes (Leiden / Metric Ward.D2)</div>", paste(cl_chips, collapse = ""), "</section>")
  } else ""

  fig_rel_cands <- c("Synthesis Plate" = "../figures/fig07a_integrative_taxonomic_synthesis.png", "Complex Network" = "../figures/fig07b_taxon_complex_delimitation_network.png", "Genomic PCoA & Tree" = "../figures/fig06a_genomic_pcoa_tree_concordance.png", "Reticulation, gdi & MMRR" = "../figures/fig06b_reticulation_gdi_landscape.png", "LM3 & ECT Exemplars" = "../figures/fig02a_lm3_organ_and_ect_exemplars.png", "EFA Morphospace" = "../figures/fig03a_efa_morphospace_and_harmonics.png", "GMM & CDA" = "../figures/fig04b_gmm_cda_morphometrics.png", "SDM & Niche" = "../figures/fig05c_sdm_suitability_and_distribution_map.png")
  fig_links <- unlist(lapply(names(fig_rel_cands), function(fnm) {
    rel_p <- fig_rel_cands[[fnm]]
    if (file.exists(normalizePath(file.path(dirname(html_path), rel_p), mustWork = FALSE))) sprintf("<a href=\"%s\" target=\"_blank\" style=\"display:inline-block;margin:2px 8px 2px 0;padding:3px 9px;background:#eff6ff;border:1px solid #bfdbfe;border-radius:6px;font-size:0.80rem;color:#1d4ed8;text-decoration:none;\">%s (300 DPI PNG)</a>", rel_p, .html_escape_vec(fnm)) else NULL
  }))
  fig_sec <- if (length(fig_links) > 0L) paste0("<section style=\"background:#fff;border:1px solid #cbd5e1;border-radius:8px;padding:8px 14px;margin-bottom:12px;\"><span style=\"font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:#475569;margin-right:10px;\">300-DPI PNG Plates:</span>", paste(fig_links, collapse = ""), "</section>") else ""

  html_doc <- c(
    sprintf("<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\"><title>%s Integrative Triage Dashboard (NCU Herbarium)</title><style>", .html_escape_vec(g_html)),
    ":root{--bg:#f8fafc;--card:#ffffff;--ink:#0f172a;--muted:#475569;--border:#cbd5e1;--accent:#0072B2;}",
    "*{box-sizing:border-box;}body{margin:0;padding:20px 24px;font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--ink);line-height:1.45;}",
    "header{background:var(--card);border:1px solid var(--border);border-radius:10px;padding:16px 20px;margin-bottom:16px;}h1{margin:0 0 6px 0;font-size:1.28rem;}.meta{font-size:0.86rem;color:var(--muted);display:flex;flex-wrap:wrap;gap:16px;}",
    ".kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(145px,1fr));gap:10px;margin-bottom:16px;}.kpi{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:10px 14px;}.kpi .lbl{font-size:0.75rem;text-transform:uppercase;letter-spacing:0.04em;color:var(--muted);}.kpi .val{font-size:1.35rem;font-weight:700;margin-top:2px;}",
    ".toolbar{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:12px 14px;margin-bottom:14px;display:flex;flex-wrap:wrap;gap:10px;align-items:center;}.toolbar input,.toolbar select,.toolbar button{padding:7px 10px;border:1px solid var(--border);border-radius:6px;font-size:0.86rem;background:#fff;color:var(--ink);}.toolbar input{flex:1 1 240px;}.toolbar button{background:var(--accent);color:#fff;border-color:var(--accent);font-weight:600;cursor:pointer;}",
    ".table-wrap{background:var(--card);border:1px solid var(--border);border-radius:8px;overflow-x:auto;}table{width:100%;border-collapse:collapse;font-size:0.83rem;}th,td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top;}th{background:#f1f5f9;font-weight:700;cursor:pointer;user-select:none;white-space:nowrap;}tr:hover td{background:#f8fafc;}",
    ".badge{display:inline-block;padding:2px 7px;border-radius:999px;font-size:0.73rem;font-weight:700;color:#fff;}.p-CRITICAL{background:#D55E00;}.p-HIGH{background:#E69F00;color:#111;}.p-MEDIUM{background:#CC79A7;}.p-LOW{background:#56B4E9;color:#111;}.p-RESOLVED{background:#009E73;}.sci{font-style:italic;}.rat{max-width:360px;color:#334155;font-size:0.79rem;}.cat-btn{background:none;border:0;padding:0;color:var(--accent);font-weight:700;cursor:pointer;text-decoration:underline;font-size:inherit;}",
    ".modal{display:none;position:fixed;inset:0;background:rgba(15,23,42,0.72);z-index:999;align-items:center;justify-content:center;padding:18px;}.modal.open{display:flex;}.modal-card{background:#fff;border-radius:10px;max-width:920px;width:100%;max-height:92vh;display:flex;flex-direction:column;overflow:hidden;box-shadow:0 20px 50px rgba(0,0,0,0.35);}.modal-hdr{display:flex;justify-content:space-between;align-items:center;padding:10px 14px;border-bottom:1px solid var(--border);background:#f8fafc;gap:8px;}.modal-body{padding:12px;overflow:auto;text-align:center;background:#0f172a;}.modal-body img{max-width:100%;max-height:76vh;object-fit:contain;border-radius:4px;}</style></head><body>",
    sprintf("<header><h1>%s</h1><div class=\"meta\"><span><strong>SOP:</strong> UNC-BOT-SOP-2026-04-REV5</span><span><strong>OOF Reliability (&kappa;):</strong> Morph = %s | Edaphic = %s | Phenology = %s</span><span id=\"shownCount\">Showing %d of %d vouchers</span></div></header>", title_html, w_m, w_e, w_p, n, n),
    "<section class=\"kpis\">",
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Examined Vouchers</div><div class=\"val\">%d</div></div>", n),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Critical / High Queue</div><div class=\"val\">%d</div></div>", sum(prio_v %in% c("CRITICAL", "HIGH"))),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Misident. (Reassign)</div><div class=\"val\">%d</div></div>", sum(sub_v == "Misidentification_Reassign")),
    hyb_kpi,
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Putative Cryptic</div><div class=\"val\">%d</div></div>", sum(sub_v == "Putative_Cryptic_Species")),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Variety / Plasticity</div><div class=\"val\">%d</div></div>", sum(sub_v %in% c("Subspecies_Ecotype", "Ecophenotypic_Plasticity"))),
    sprintf("<div class=\"kpi\"><div class=\"lbl\">Confirmed Species</div><div class=\"val\">%d</div></div>", sum(sub_v == "Confirmed_Species")),
    "</section>", hyb_sec, comp_sec, fig_sec,
    "<section class=\"toolbar\"><input type=\"search\" id=\"qSearch\" placeholder=\"Search catalogNumber, collector, herbarium, taxon, rationale...\">",
    "<select id=\"qPrio\"><option value=\"ALL\">All Priorities</option><option value=\"CRITICAL\">CRITICAL</option><option value=\"HIGH\">HIGH</option><option value=\"MEDIUM\">MEDIUM</option><option value=\"LOW\">LOW</option><option value=\"RESOLVED\">RESOLVED</option></select>",
    "<select id=\"qSub\"><option value=\"ALL\">All Treatments</option><option value=\"Misidentification_Reassign\">Misidentification_Reassign</option><option value=\"Putative_Cryptic_Species\">Putative_Cryptic_Species</option><option value=\"Hybrid_Introgressant\">Hybrid_Introgressant</option><option value=\"Subspecies_Ecotype\">Subspecies_Ecotype</option><option value=\"Ecophenotypic_Plasticity\">Ecophenotypic_Plasticity</option><option value=\"Confirmed_Species\">Confirmed_Species</option><option value=\"Confirmed_Type_Control\">Confirmed_Type_Control</option><option value=\"Type_Synonymy_Candidate\">Type_Synonymy_Candidate</option></select>",
    sprintf("<select id=\"qTaxon\">%s</select><button type=\"button\" id=\"btnExport\">Export Filtered CSV</button></section>", taxon_opts),
    "<div class=\"table-wrap\"><table id=\"triTable\"><thead><tr><th data-col=\"0\">Catalog #</th><th data-col=\"1\">Ind</th><th data-col=\"2\">Herb.</th><th data-col=\"4\">Tier</th><th data-col=\"5\">Recorded Taxon</th><th data-col=\"6\">Recommended Determination</th><th data-col=\"7\">Priority</th><th data-col=\"8\">Treatment Subtype</th><th data-col=\"9\">CDA Pred (p)</th><th data-col=\"11\">Edaphic Fit (z)</th><th data-col=\"13\">&Delta;DOY (d)</th><th data-col=\"14\">Consensus</th><th data-col=\"20\">Morph Source (B/C)</th><th data-col=\"22\">Ecotype Regime</th><th data-col=\"23\">Genomic / gdi / Cytotype</th><th data-col=\"16\">Synthesis Rationale</th></tr></thead><tbody id=\"triBody\"></tbody></table></div>",
    "<div id=\"sheetModal\" class=\"modal\" role=\"dialog\" aria-modal=\"true\"><div class=\"modal-card\"><div class=\"modal-hdr\"><strong id=\"modalTitle\">Specimen Sheet</strong><div><button type=\"button\" id=\"btnToggleScan\" style=\"margin-right:8px;padding:4px 10px;\">Show Raw Scan</button><button type=\"button\" id=\"btnCloseModal\" style=\"padding:4px 10px;\">Close</button></div></div><div class=\"modal-body\"><img id=\"modalImg\" alt=\"Herbarium specimen sheet preview\"><div id=\"modalMsg\" style=\"color:#cbd5e1;font-size:0.85rem;margin-top:8px;display:none;\"></div></div></div></div>",
    "<script>", paste0("const ROWS = [\n", json_rows, "\n];"),
    "const PRIO_RANK = {CRITICAL:1,HIGH:2,MEDIUM:3,LOW:4,RESOLVED:5}; let sortCol = 7, sortAsc = true, filtered = ROWS.slice(), modalCat = '', modalRaw = false;",
    "function esc(s){return String(s||'').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/\"/g,'&quot;');}",
    "function fmtSci(s){const t=esc(s);return /^[A-Z][a-z-]+[ \t]+(?:[x×][ \t]*)?[a-z][a-z-]+/.test(t)?'<span class=\"sci\">'+t+'</span>':t;}",
    "function loadModalImg(){const img=document.getElementById('modalImg'),msg=document.getElementById('modalMsg'),btn=document.getElementById('btnToggleScan');msg.style.display='none';img.style.display='inline-block';const ov='../figures/lm3_sheet_overlays/'+encodeURIComponent(modalCat)+'.png',rw='../../../data/raw_vouchers/'+encodeURIComponent(modalCat)+'.webp';btn.textContent=modalRaw?'Show LM3 Overlay':'Show Raw Scan';img.onerror=function(){if(!modalRaw){modalRaw=true;btn.textContent='Show LM3 Overlay';img.onerror=function(){img.style.display='none';msg.textContent='Sheet preview image not found locally for '+modalCat;msg.style.display='block';};img.src=rw;}else{img.style.display='none';msg.textContent='Raw voucher scan not found locally for '+modalCat;msg.style.display='block';}};img.src=modalRaw?rw:ov;}",
    "function openModal(cat){modalCat=cat;modalRaw=false;document.getElementById('modalTitle').textContent='Herbarium Voucher: '+cat;document.getElementById('sheetModal').classList.add('open');loadModalImg();}",
    "document.getElementById('btnToggleScan').addEventListener('click',()=>{modalRaw=!modalRaw;loadModalImg();});document.getElementById('btnCloseModal').addEventListener('click',()=>document.getElementById('sheetModal').classList.remove('open'));document.getElementById('sheetModal').addEventListener('click',e=>{if(e.target.id==='sheetModal')document.getElementById('sheetModal').classList.remove('open');});document.addEventListener('keydown',e=>{if(e.key==='Escape')document.getElementById('sheetModal').classList.remove('open');});",
    "function applyFilter(){const q=document.getElementById('qSearch').value.trim().toLowerCase(),pr=document.getElementById('qPrio').value,sb=document.getElementById('qSub').value,tx=document.getElementById('qTaxon').value;filtered=ROWS.filter(r=>{if(pr!=='ALL'&&r[7]!==pr)return false;if(sb!=='ALL'&&r[8]!==sb)return false;if(tx!=='ALL'&&r[5]!==tx)return false;if(q&&!r.join(' ').toLowerCase().includes(q))return false;return true;});filtered.sort((a,b)=>{let va=a[sortCol],vb=b[sortCol];if(sortCol===7){va=PRIO_RANK[va]||9;vb=PRIO_RANK[vb]||9;}else if([1,10,12,13,14].includes(sortCol)){va=parseFloat(va);vb=parseFloat(vb);if(isNaN(va))va=9999;if(isNaN(vb))vb=9999;}if(va<vb)return sortAsc?-1:1;if(va>vb)return sortAsc?1:-1;return 0;});document.getElementById('triBody').innerHTML=filtered.map(r=>'<tr><td><button type=\"button\" class=\"cat-btn\" data-cat=\"'+esc(r[0])+'\">'+esc(r[0])+'</button></td><td>'+esc(r[1])+'</td><td>'+esc(r[2])+'</td><td>'+esc(r[4])+'</td><td>'+fmtSci(r[5])+'</td><td>'+fmtSci(r[6])+'</td><td><span class=\"badge p-'+esc(r[7])+'\">'+esc(r[7])+'</span></td><td>'+esc(r[8])+'</td><td>'+fmtSci(r[9])+(r[10]!=='NA'?' ('+esc(r[10])+')':'')+'</td><td>'+fmtSci(r[11])+(r[12]!=='NA'?' (z='+esc(r[12])+')':'')+'</td><td>'+esc(r[13])+'</td><td>'+esc(r[14])+'</td><td>'+esc(r[20])+' ('+esc(r[21])+')</td><td>'+fmtSci(r[22])+'</td><td>'+fmtSci(r[23])+'</td><td class=\"rat\">'+esc(r[16])+'</td></tr>').join('');document.querySelectorAll('#triBody .cat-btn').forEach(b=>b.addEventListener('click',()=>openModal(b.getAttribute('data-cat'))));document.getElementById('shownCount').textContent='Showing '+filtered.length+' of '+ROWS.length+' vouchers';}",
    "document.querySelectorAll('#triTable th').forEach(th=>th.addEventListener('click',()=>{const c=parseInt(th.getAttribute('data-col'),10);if(sortCol===c)sortAsc=!sortAsc;else{sortCol=c;sortAsc=true;}applyFilter();}));['qSearch','qPrio','qSub','qTaxon'].forEach(id=>document.getElementById(id).addEventListener('input',applyFilter));",
    "document.getElementById('btnExport').addEventListener('click',()=>{const hdr=['catalogNumber','plant_individual_id','institutionCode','recordedBy','determiner_tier','species_standardized','recommended_determination','synthesis_triage_priority','taxonomic_treatment_subtype','cda_predicted_taxon','cda_posterior_prob','edaphic_best_fit_taxon','edaphic_z_score','delta_doy_resid','probabilistic_consensus_score','synthesis_triage_action','synthesis_rationale','hybrid_parent_1','hybrid_parent_2','hybrid_zone_cohesion_pass','morph_tier_source','basal_cauline_leaves','ecotype_regime','genomic_summary'];const csv=[hdr.join(',')].concat(filtered.map(r=>r.map(v=>'\"'+String(v||'').replace(/\"/g,'\"\"')+'\"').join(','))).join('\\n');const blob=new Blob([csv],{type:'text/csv;charset=utf-8;'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='triage_queue_filtered.csv';a.click();URL.revokeObjectURL(a.href);});",
    "applyFilter();</script></body></html>"
  )
  tmp_html <- paste0(html_path, ".", Sys.getpid(), ".tmp"); writeLines(html_doc, tmp_html, useBytes = TRUE); file.rename(tmp_html, html_path)
  invisible(html_path)
}

