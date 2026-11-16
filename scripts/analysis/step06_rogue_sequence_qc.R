# ==============================================================================
# Module: scripts/analysis/step06_rogue_sequence_qc.R
# Project: Plantmatics | Affiliation: UNC-Chapel Hill Herbarium (NCU)
# Purpose: MAD rogue tip quarantine (validated Newick/NEXUS RelTime trees) and
#   sub-marker-partitioned strand-aware FASTA outlier quarantine for Step 06.
# ==============================================================================

.empty_rogue_audit_df <- function() {
  data.frame(
    source_file = character(0), locus_or_tree = character(0), accession_or_tip = character(0),
    taxon = character(0), metric_value = numeric(0), mad_z_score = numeric(0),
    outlier_reason = character(0), quarantined = logical(0), stringsAsFactors = FALSE
  )
}

.kmers_11_qc <- function(s) {
  s_u <- toupper(as.character(s)[1L]); if (is.na(s_u) || !nzchar(s_u)) return(character(0))
  if (exists(".kmers_11", mode = "function")) return(.kmers_11(s_u))
  n <- nchar(s_u); if (n < 11L) return(character(0))
  st <- seq_len(n - 10L); k <- substring(s_u, st, st + 10L); unique(k[!grepl("[^ACGT]", k)])
}

.revcomp_qc <- function(s) {
  s_u <- toupper(as.character(s)[1L]); if (is.na(s_u) || !nzchar(s_u)) return("")
  if (exists(".revcomp_seq", mode = "function")) return(.revcomp_seq(s_u))
  chartr("ACGTU", "TGCAA", paste(rev(strsplit(s_u, "", fixed = TRUE)[[1L]]), collapse = ""))
}

.write_audit_csv_safe <- function(df, path) {
  if (is.null(path) || !nzchar(path)) return(invisible(NULL))
  if (exists(".write_csv_safe", mode = "function")) return(.write_csv_safe(df, path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(df, path) else utils::write.csv(df, path, row.names = FALSE)
}

.norm_taxon_qc <- function(raw) {
  if (exists(".normalize_public_packera_taxon", mode = "function")) return(.normalize_public_packera_taxon(raw))
  if (is.null(raw) || length(raw) == 0L) return(character(0))
  sen_map <- c(tomentosus = "dubia", tomentosa = "dubia", aureus = "aurea", anonymus = "anonyma", obovatus = "obovata",
    pauperculus = "paupercula", schweinitzianus = "schweinitziana", millefolium = "millefolium", millefolia = "millefolium",
    antennariifolius = "antennariifolia", glabellus = "glabella", plattensis = "plattensis", pseudaureus = "pseudaurea",
    indecorus = "indecora", canus = "cana", loratifolius = "loratifolia", strepthanthifolius = "streptanthifolia",
    streptanthifolius = "streptanthifolia", memmingeri = "xmemmingeri", xmemmingeri = "xmemmingeri")
  vapply(as.character(raw), function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return("")
    s <- trimws(gsub("[\"']", "", sub("^>", "", trimws(x))))
    if (grepl("\\|", s)) {
      parts <- Filter(nzchar, trimws(strsplit(s, "\\|", fixed = FALSE)[[1L]]))
      hit <- parts[grepl("^(Packera|Senecio|P\\.)([_[:space:]]|$)", parts, ignore.case = TRUE)]
      s <- if (length(hit) > 0L) hit[1L] else if (length(parts) >= 2L) parts[2L] else parts[1L]
    }
    s <- trimws(gsub("[[:space:]]+", " ", gsub("_+", " ", s)))
    s <- sub("^([A-Za-z]+)[[:space:]]+[0-9]+([A-Za-z][A-Za-z-]*)$", "\\1 \\2", sub("^(Packera|P\\.)[[:space:]]+sp\\.?[[:space:]]+(Packera|Senecio|P\\.)", "\\2", s, ignore.case = TRUE))
    reg <- regmatches(s, regexec("(^|[^A-Za-z])(Packera|Senecio|P\\.)[[:space:]]+[0-9]*([A-Za-z][A-Za-z.-]*)", s, ignore.case = TRUE))[[1L]]
    if (length(reg) < 4L) return(s)
    ep <- tolower(sub("[-.].*$", "", reg[4L])); if (ep %in% c("sp", "cf", "aff", "unk", "unknown")) return("")
    if (tolower(reg[3L]) == "senecio") return(if (ep %in% names(sen_map)) paste("Packera", sen_map[[ep]]) else paste("Senecio", ep))
    paste("Packera", if (ep %in% names(sen_map)) sen_map[[ep]] else ep)
  }, character(1L), USE.NAMES = FALSE)
}

.is_designated_outgroup_taxon <- function(norm_taxon) {
  if (is.null(norm_taxon) || length(norm_taxon) == 0L) return(logical(0))
  x <- trimws(as.character(norm_taxon)); ok <- !is.na(x) & nzchar(x); res <- rep(FALSE, length(x))
  pat <- "(^|[^A-Za-z])(Jacobaea|Senecio|Tephroseris|Pericallis|Doronicum|Erechtites|Emilia|Crassocephalum|Gynura|Arrhenechthites|Roldana|Telanthophora|Digitacalia|Pippenalia|Blennosperma|Barkleyanthus|Petasites|Pittocaulon|Psacalium|Robinsonecio|Tussilago|Sinacalia|Cacaliopsis|Brachyglottis|Cotula|Artemisia|Achillea|Anthemis|Elekmania|Rugelia|Luina)[[:space:]_]"
  if (any(ok)) res[ok] <- grepl(pat, x[ok], ignore.case = TRUE); res
}

.is_valid_phylo_qc <- function(tr) {
  if (is.null(tr) || !inherits(tr, "phylo") || is.null(tr$edge) || !is.matrix(tr$edge) || ncol(tr$edge) != 2L) return(FALSE)
  nt <- length(tr$tip.label); nn <- suppressWarnings(as.integer(tr$Nnode)[1L]); if (nt < 3L || is.na(nn) || nn < 1L || anyNA(tr$edge)) return(FALSE)
  p <- suppressWarnings(as.integer(tr$edge[, 1L])); ch <- suppressWarnings(as.integer(tr$edge[, 2L])); max_id <- nt + nn
  if (anyNA(p) || anyNA(ch) || any(p <= nt | p > max_id | ch < 1L | ch > max_id | p == ch)) return(FALSE)
  if (!identical(sort(ch[ch <= nt]), seq_len(nt)) || anyDuplicated(ch) > 0L) return(FALSE)
  roots <- setdiff(unique(p), unique(ch)); length(roots) == 1L && roots[1L] == (nt + 1L)
}

.read_nexus_phylo_qc <- function(path) {
  if (!nzchar(path) || !file.exists(path) || !requireNamespace("ape", quietly = TRUE)) return(NULL)
  lns <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (length(lns) == 0L || !any(grepl("^#NEXUS", trimws(head(lns, 10L)), ignore.case = TRUE))) return(NULL)
  tr_idx <- grep("^[[:space:]]*tree[[:space:]]+", lns, ignore.case = TRUE); if (length(tr_idx) == 0L) return(NULL)
  nwk <- gsub("\\[[^]]*\\]", "", sub("^[[:space:]]*tree[[:space:]]+[^=]+=[[:space:]]*", "", lns[tr_idx[1L]], ignore.case = TRUE))
  if (!grepl(";[[:space:]]*$", nwk)) nwk <- paste0(nwk, ";")
  tr <- tryCatch(ape::read.tree(text = nwk), error = function(e) NULL); if (inherits(tr, "multiPhylo") && length(tr) > 0L) tr <- tr[[1L]]
  if (!.is_valid_phylo_qc(tr)) return(NULL)
  tl_s <- grep("^[[:space:]]*translate([[:space:]]|;|$)", lns, ignore.case = TRUE)
  if (length(tl_s) > 0L && tl_s[1L] < tr_idx[1L]) {
    s_i <- tl_s[1L] + 1L; semi <- grep(";", lns[s_i:tr_idx[1L]]); e_i <- if (length(semi) > 0L) (s_i + semi[1L] - 1L) else (tr_idx[1L] - 1L)
    if (e_i >= s_i) {
      tmap <- character(0)
      for (ln in lns[s_i:e_i]) { cl <- trimws(sub("[,;][[:space:]]*$", "", trimws(ln))); m <- regmatches(cl, regexec("^([0-9]+)[[:space:]]+(.+)$", cl))[[1L]]; if (length(m) >= 3L) tmap[[m[2L]]] <- gsub("^['\"]|['\"]$", "", trimws(m[3L])) }
      hit <- tr$tip.label %in% names(tmap); if (any(hit)) tr$tip.label[hit] <- unname(tmap[tr$tip.label[hit]])
    }
  }
  storage.mode(tr$edge) <- "integer"; if (is.null(tr$edge.length)) tr$edge.length <- rep(1.0, nrow(tr$edge)); tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length < 0] <- 0.0
  nt <- .norm_taxon_qc(tr$tip.label); tr$tip.label <- ifelse(nzchar(nt), nt, tr$tip.label); tr
}

screen_public_tree_rogue_tips <- function(tr, tree_label = "public_tree", mad_z_thresh = 5.0, min_ingroup_tips = 6L, source_file = tree_label) {
  if (is.null(tr) || !inherits(tr, "phylo")) return(list(tree = tr, audit_df = .empty_rogue_audit_df()))
  if (!requireNamespace("ape", quietly = TRUE) || !.is_valid_phylo_qc(tr)) return(list(tree = NULL, audit_df = .empty_rogue_audit_df()))
  storage.mode(tr$edge) <- "integer"
  if (length(tr$tip.label) < min_ingroup_tips || is.null(tr$edge.length)) return(list(tree = tr, audit_df = .empty_rogue_audit_df()))
  el_sd <- suppressWarnings(stats::sd(tr$edge.length)); if (!is.finite(el_sd) || el_sd < 1e-8) return(list(tree = tr, audit_df = .empty_rogue_audit_df()))
  orig_tips <- tr$tip.label; norm_tips <- .norm_taxon_qc(orig_tips)
  is_pk <- grepl("^(Packera|P\\.)([[:space:]_]|$)", norm_tips, ignore.case = TRUE)
  is_out <- .is_designated_outgroup_taxon(norm_tips) | (!is_pk & .is_designated_outgroup_taxon(orig_tips))
  ing_idx <- which(!is_out & is_pk); if (length(ing_idx) < min_ingroup_tips) ing_idx <- which(!is_out)
  if (length(ing_idx) < min_ingroup_tips) return(list(tree = tr, audit_df = .empty_rogue_audit_df()))

  L_term <- tr$edge.length[match(ing_idx, tr$edge[, 2L])]; coph <- ape::cophenetic.phylo(tr)
  D_nn <- vapply(ing_idx, function(ii) { dv <- coph[ii, setdiff(ing_idx, ii)]; dv <- dv[is.finite(dv)]; if (length(dv) == 0L) NA_real_ else min(dv) }, numeric(1L))
  RTT <- ape::node.depth.edgelength(tr)[ing_idx]; med_L <- stats::median(L_term, na.rm = TRUE)
  z_term <- (L_term - med_L) / max(c(1.4826 * stats::mad(L_term, constant = 1, na.rm = TRUE), 1e-6), na.rm = TRUE)
  med_RTT <- stats::median(RTT, na.rm = TRUE); z_rtt <- abs(RTT - med_RTT) / max(c(1.4826 * stats::mad(RTT, constant = 1, na.rm = TRUE), 1e-6), na.rm = TRUE)
  med_Dnn <- stats::median(D_nn, na.rm = TRUE)
  cand_pos <- which(is.finite(z_term) & (z_term > mad_z_thresh) & (L_term > 3.0 * max(med_L, 1e-4)) & is.finite(D_nn) & is.finite(med_Dnn) & (D_nn > 2.5 * med_Dnn))
  if (length(cand_pos) == 0L) return(list(tree = tr, audit_df = .empty_rogue_audit_df()))
  cand_pos <- cand_pos[order(z_term[cand_pos], decreasing = TRUE)]; q_flag <- seq_along(cand_pos) <= max(0L, length(ing_idx) - 4L); drop_idx <- ing_idx[cand_pos[q_flag]]
  if (length(drop_idx) > 0L) {
    tmp_lbls <- paste0("__ROGUE_TIP_", seq_along(orig_tips), "__"); tr_work <- tr; tr_work$tip.label <- tmp_lbls
    tr_drop <- tryCatch(ape::drop.tip(tr_work, tmp_lbls[drop_idx]), error = function(e) NULL)
    if (!is.null(tr_drop) && inherits(tr_drop, "phylo") && length(tr_drop$tip.label) >= 4L) { tr_drop$tip.label <- orig_tips[match(tr_drop$tip.label, tmp_lbls)]; tr <- tr_drop } else q_flag[] <- FALSE
  }
  audit_df <- data.frame(
    source_file = rep(as.character(source_file)[1L], length(cand_pos)), locus_or_tree = rep(as.character(tree_label)[1L], length(cand_pos)),
    accession_or_tip = orig_tips[ing_idx[cand_pos]], taxon = norm_tips[ing_idx[cand_pos]], metric_value = round(L_term[cand_pos], 6),
    mad_z_score = round(z_term[cand_pos], 4), outlier_reason = rep("rogue_long_pendant_tip", length(cand_pos)), quarantined = q_flag, stringsAsFactors = FALSE
  )
  list(tree = tr, audit_df = audit_df)
}

.submarker_qc <- function(hdr) {
  h <- as.character(hdr); out <- rep("main", length(h))
  has_i1 <- grepl("(ITS1|spacer[[:space:]_]*1)([^0-9]|$)", h, ignore.case = TRUE); has_i2 <- grepl("(ITS2|spacer[[:space:]_]*2)([^0-9]|$)", h, ignore.case = TRUE)
  has_its <- grepl("internal[[:space:]]+transcribed|(^|[^A-Za-z])ITS", h, ignore.case = TRUE); has_ets <- grepl("external[[:space:]]+transcribed|(^|[^A-Za-z])ETS([^A-Za-z]|$)", h, ignore.case = TRUE)
  out[grepl("rbcL", h, ignore.case = TRUE)] <- "rbcL"; out[grepl("matK", h, ignore.case = TRUE)] <- "matK"
  out[grepl("rpoC", h, ignore.case = TRUE)] <- "rpoC"; out[grepl("trnL|trnH|psbA|ndhF|rps16", h, ignore.case = TRUE)] <- "cp_spacer"
  out[has_i1 & !has_i2 & out == "main"] <- "ITS1_only"; out[has_i2 & !has_i1 & out == "main"] <- "ITS2_only"; out[has_ets & !has_its & out == "main"] <- "ETS_only"; out
}

.screen_single_fasta_pool <- function(refs_df, loc_str, fwd_km, acc_vec, norm_tx, is_pk, is_out, mad_z_thresh = 5.0) {
  n <- nrow(refs_df); seq_u <- toupper(gsub("[-._[:space:]]", "", as.character(refs_df$seq))); seq_u[is.na(seq_u)] <- ""
  acgt_n <- nchar(gsub("[^ACGT]", "", seq_u)); ambig_frac <- ifelse(nchar(seq_u) > 0L, 1.0 - acgt_n / nchar(seq_u), 1.0); med_len <- stats::median(acgt_n, na.rm = TRUE)
  drop_g1 <- (acgt_n < 120L) | (n >= 5L & acgt_n < min(2000.0, 0.30 * med_len)) | (ambig_frac > 0.15)
  mad_len <- max(c(1.4826 * stats::mad(acgt_n, constant = 1, na.rm = TRUE), 1.0), na.rm = TRUE); mad_amb <- max(c(1.4826 * stats::mad(ambig_frac, constant = 1, na.rm = TRUE), 0.01), na.rm = TRUE)
  g1_metric <- ifelse(ambig_frac > 0.15, round(ambig_frac, 4), as.numeric(acgt_n))
  g1_z <- ifelse(ambig_frac > 0.15, (ambig_frac - stats::median(ambig_frac, na.rm = TRUE)) / mad_amb, (acgt_n - med_len) / mad_len)

  non_g1 <- which(!drop_g1); sk_fwd <- fwd_km; sk_rev <- NULL
  if (length(non_g1) > 0L && stats::median(vapply(fwd_km[non_g1], length, integer(1L))) > 6000L) {
    for (suf in list(c("AAA", "TTT"), c("AA", "TT"))) {
      cand_f <- lapply(fwd_km, function(k) k[endsWith(k, suf[1L])])
      if (all(vapply(cand_f[non_g1], length, integer(1L)) >= 40L)) {
        sk_fwd <- cand_f; sk_rev <- lapply(fwd_km, function(k) { kt <- k[startsWith(k, suf[2L])]; if (length(kt) > 0L) vapply(kt, .revcomp_qc, character(1L), USE.NAMES = FALSE) else character(0) }); break
      }
    }
  }
  if (is.null(sk_rev)) sk_rev <- lapply(seq_len(n), function(i) if (drop_g1[i]) character(0) else if (nzchar(seq_u[i])) .kmers_11_qc(.revcomp_qc(seq_u[i])) else if (length(sk_fwd[[i]]) > 0L) unique(vapply(sk_fwd[[i]], .revcomp_qc, character(1L), USE.NAMES = FALSE)) else character(0))

  drop_g2 <- rep(FALSE, n); g2_metric <- rep(NA_real_, n); g2_z <- rep(NA_real_, n)
  eval_idx <- which(!drop_g1 & !is_out); pool_idx <- which(!drop_g1 & is_pk); if (length(pool_idx) < 3L) pool_idx <- eval_idx
  if (n >= 5L && length(eval_idx) >= 4L && length(pool_idx) >= 3L) {
    max_sim <- vapply(eval_idx, function(ii) {
      li <- length(sk_fwd[[ii]]); others <- setdiff(pool_idx, ii); if (li == 0L || length(others) == 0L) return(0.0)
      max(vapply(others, function(jj) { lj <- length(sk_fwd[[jj]]); if (lj == 0L) 0.0 else { sh <- max(sum(sk_fwd[[ii]] %in% sk_fwd[[jj]]), sum(sk_rev[[ii]] %in% sk_fwd[[jj]])); max(sh / min(li, lj), sh / max(1L, li + lj - sh)) } }, numeric(1L)))
    }, numeric(1L))
    med_s <- stats::median(max_sim, na.rm = TRUE); mad_s <- max(c(1.4826 * stats::mad(max_sim, constant = 1, na.rm = TRUE), 0.01), na.rm = TRUE)
    z_sim <- (max_sim - med_s) / mad_s; drop_g2[eval_idx[is.finite(z_sim) & (max_sim < 0.08) & (z_sim < -mad_z_thresh)]] <- TRUE; g2_metric[eval_idx] <- max_sim; g2_z[eval_idx] <- z_sim
  }

  drop_g3 <- rep(FALSE, n); g3_metric <- rep(NA_real_, n); g3_z <- rep(NA_real_, n); surv_idx <- which(!drop_g1 & !drop_g2 & is_pk)
  if (length(surv_idx) >= 5L) {
    m <- length(surv_idx); D <- matrix(0.0, nrow = m, ncol = m); same_str <- matrix(TRUE, nrow = m, ncol = m)
    for (a in seq_len(m - 1L)) for (b in (a + 1L):m) {
      ia <- surv_idx[a]; ib <- surv_idx[b]; la <- length(sk_fwd[[ia]]); lb <- length(sk_fwd[[ib]])
      if (la == 0L || lb == 0L) { D[a, b] <- D[b, a] <- 1.0; next }
      sf <- sum(sk_fwd[[ia]] %in% sk_fwd[[ib]]); sr <- sum(sk_rev[[ia]] %in% sk_fwd[[ib]]); if (sr > sf) same_str[a, b] <- same_str[b, a] <- FALSE
      sh <- max(sf, sr); D[a, b] <- D[b, a] <- 1.0 - (sh / (if (max(la, lb) > 1.5 * min(la, lb)) max(1L, min(la, lb)) else max(1L, la + lb - sh)))
    }
    if (requireNamespace("ape", quietly = TRUE)) {
      lens_s <- nchar(seq_u[surv_idx])
      for (ul in unique(lens_s[lens_s >= 100L & lens_s <= 8000L])) {
        sub_p <- which(lens_s == ul); if (length(sub_p) < 5L) next
        ch_mat <- do.call(rbind, strsplit(tolower(seq_u[surv_idx[sub_p]]), "", fixed = TRUE))
        k80 <- suppressWarnings(tryCatch(as.matrix(ape::dist.dna(ape::as.DNAbin(ch_mat), model = "K80", pairwise.deletion = TRUE)), error = function(e) NULL))
        if (!is.null(k80) && nrow(k80) == length(sub_p)) for (u in seq_len(length(sub_p) - 1L)) for (v in (u + 1L):length(sub_p)) {
          a <- sub_p[u]; b <- sub_p[v]; if (is.finite(k80[u, v]) && same_str[a, b]) D[a, b] <- D[b, a] <- 0.5 * D[a, b] + 0.5 * min(1.0, k80[u, v])
        }
      }
    }
    sp_surv <- norm_tx[surv_idx]; sp_tab <- table(sp_surv); other_all <- setdiff(names(sp_tab[sp_tab >= 2L]), "")
    for (sp in setdiff(names(sp_tab[sp_tab >= 3L]), "")) {
      sp_pos <- which(sp_surv == sp); d_own <- vapply(seq_along(sp_pos), function(k) stats::median(D[sp_pos[k], sp_pos[-k]], na.rm = TRUE), numeric(1L))
      c_med <- stats::median(d_own, na.rm = TRUE); s_mad <- max(c(1.4826 * stats::mad(d_own, constant = 1, na.rm = TRUE), 0.015), na.rm = TRUE); thresh <- c_med + mad_z_thresh * s_mad
      for (k in seq_along(sp_pos)) {
        if (!is.finite(d_own[k]) || d_own[k] <= thresh) next
        thr_o <- 0.40 * d_own[k]
        if (any(vapply(setdiff(other_all, sp), function(spo) { d_o <- D[sp_pos[k], which(sp_surv == spo)]; d_o <- d_o[is.finite(d_o)]; length(d_o) >= 2L && min(d_o) < thr_o && (stats::median(d_o) < thr_o || sum(d_o < thr_o) >= 2L) }, logical(1L)))) {
          gi <- surv_idx[sp_pos[k]]; drop_g3[gi] <- TRUE; g3_metric[gi] <- d_own[k]; g3_z[gi] <- (d_own[k] - c_med) / s_mad
        }
      }
    }
  }
  drop_all <- drop_g1 | drop_g2 | drop_g3; if (!any(drop_all)) return(list(drop_mask = drop_all, audit_df = .empty_rogue_audit_df()))
  f_idx <- which(drop_all)
  audit_df <- data.frame(
    source_file = as.character(refs_df$ref_source[f_idx]), locus_or_tree = rep(loc_str, length(f_idx)), accession_or_tip = acc_vec[f_idx], taxon = norm_tx[f_idx],
    metric_value = round(ifelse(drop_g1[f_idx], g1_metric[f_idx], ifelse(drop_g2[f_idx], g2_metric[f_idx], g3_metric[f_idx])), 6),
    mad_z_score = round(ifelse(drop_g1[f_idx], g1_z[f_idx], ifelse(drop_g2[f_idx], g2_z[f_idx], g3_z[f_idx])), 4),
    outlier_reason = ifelse(drop_g1[f_idx], "excess_ambiguity_or_fragment", ifelse(drop_g2[f_idx], "non_homologous_or_contaminant_accession", "mislabeled_genbank_conspecific_outlier")),
    quarantined = rep(TRUE, length(f_idx)), stringsAsFactors = FALSE
  )
  list(drop_mask = drop_all, audit_df = audit_df)
}

screen_public_fasta_outliers <- function(refs_df, locus_name, kmers_list = NULL, mad_z_thresh = 5.0) {
  if (is.null(refs_df) || !is.data.frame(refs_df) || nrow(refs_df) == 0L || !("seq" %in% names(refs_df))) {
    empty_r <- if (is.data.frame(refs_df)) refs_df else data.frame(header = character(0), seq = character(0), stringsAsFactors = FALSE)
    return(list(clean_refs = empty_r, clean_kmers = list(), audit_df = .empty_rogue_audit_df()))
  }
  n <- nrow(refs_df); loc_str <- as.character(locus_name)[1L]
  if (!("header" %in% names(refs_df))) refs_df$header <- if ("accession" %in% names(refs_df)) as.character(refs_df$accession) else paste0("seq_", seq_len(n))
  if (!("taxon" %in% names(refs_df))) refs_df$taxon <- .norm_taxon_qc(refs_df$header)
  miss_tx <- is.na(refs_df$taxon) | !nzchar(trimws(as.character(refs_df$taxon))); if (any(miss_tx)) refs_df$taxon[miss_tx] <- .norm_taxon_qc(refs_df$header[miss_tx])
  if (!("ref_source" %in% names(refs_df))) refs_df$ref_source <- loc_str; refs_df$ref_source[is.na(refs_df$ref_source) | !nzchar(trimws(as.character(refs_df$ref_source)))] <- loc_str
  acc_vec <- sub("\\|.*$", "", sub("[[:space:]].*$", "", sub("^>", "", trimws(as.character(refs_df$header))))); acc_vec[!nzchar(acc_vec)] <- as.character(refs_df$header[!nzchar(acc_vec)])
  norm_tx <- .norm_taxon_qc(refs_df$taxon); is_pk <- grepl("^(Packera|P\\.)([[:space:]_]|$)", norm_tx, ignore.case = TRUE)
  is_out <- .is_designated_outgroup_taxon(norm_tx) | (!is_pk & .is_designated_outgroup_taxon(refs_df$header)); if (!any(is_pk)) is_pk <- !is_out
  fwd_km <- if (!is.null(kmers_list) && length(kmers_list) == n) kmers_list else lapply(toupper(gsub("[-._[:space:]]", "", as.character(refs_df$seq))), .kmers_11_qc)
  grp <- paste(refs_df$ref_source, .submarker_qc(refs_df$header), sep = "::"); drop_mask <- rep(FALSE, n); audits <- list()
  for (idx in split(seq_len(n), grp)) {
    sub_res <- .screen_single_fasta_pool(refs_df[idx, , drop = FALSE], loc_str, fwd_km[idx], acc_vec[idx], norm_tx[idx], is_pk[idx], is_out[idx], mad_z_thresh)
    drop_mask[idx] <- sub_res$drop_mask; if (nrow(sub_res$audit_df) > 0L) audits[[length(audits) + 1L]] <- sub_res$audit_df
  }
  audit_df <- if (length(audits) > 0L) do.call(rbind, audits) else .empty_rogue_audit_df()
  rownames(audit_df) <- NULL; all_drop <- all(drop_mask); if (all_drop && nrow(audit_df) > 0L) audit_df$quarantined <- FALSE
  keep <- if (all_drop) rep(TRUE, n) else !drop_mask
  list(clean_refs = refs_df[keep, , drop = FALSE], clean_kmers = fwd_km[keep], audit_df = audit_df)
}

.screen_public_backbone_trees_bundle <- function(trees, found_files, phytop_tables, all_gf, read_tree_fn) {
  audits <- list(); found_files <- if (is.null(found_files)) character(0) else found_files; all_gf <- if (is.null(all_gf)) character(0) else all_gf
  for (nm in names(trees)) {
    if (is.null(trees[[nm]])) next
    src <- if (!is.null(names(found_files)) && (nm %in% names(found_files))) basename(found_files[[nm]]) else { hit <- found_files[grepl(nm, basename(found_files), ignore.case = TRUE)]; if (length(hit) > 0L) basename(hit[1L]) else nm }
    qc <- screen_public_tree_rogue_tips(trees[[nm]], tree_label = nm, source_file = src)
    trees[[nm]] <- qc$tree; if (nrow(qc$audit_df) > 0L) audits[[length(audits) + 1L]] <- qc$audit_df
  }
  for (nm in names(phytop_tables)) {
    pt_tr <- attr(phytop_tables[[nm]], "phytop_tree"); if (is.null(pt_tr)) next
    qc <- screen_public_tree_rogue_tips(pt_tr, tree_label = paste0("phytop:", nm), source_file = sub("\\.(tsv|csv)$", ".tre", nm, ignore.case = TRUE))
    attr(phytop_tables[[nm]], "phytop_tree") <- qc$tree; if (nrow(qc$audit_df) > 0L) audits[[length(audits) + 1L]] <- qc$audit_df
  }
  rt_files <- setdiff(all_gf[grepl("reltime", all_gf, ignore.case = TRUE) & grepl("\\.(tre|nwk|tree)$", basename(all_gf), ignore.case = TRUE)], unname(found_files))
  reltime_trees <- list()
  for (rf in rt_files) {
    r_tr <- .read_nexus_phylo_qc(rf); if (is.null(r_tr) && is.function(read_tree_fn)) r_tr <- read_tree_fn(rf)
    if (is.null(r_tr) || !.is_valid_phylo_qc(r_tr)) next
    qc <- screen_public_tree_rogue_tips(r_tr, tree_label = paste0("reltime:", basename(rf)), source_file = basename(rf))
    if (!is.null(qc$tree)) reltime_trees[[basename(rf)]] <- qc$tree; if (nrow(qc$audit_df) > 0L) audits[[length(audits) + 1L]] <- qc$audit_df
  }
  list(trees = trees, reltime_trees = reltime_trees, phytop_tables = phytop_tables, audit_df = if (length(audits) > 0L) do.call(rbind, audits) else .empty_rogue_audit_df())
}

.write_public_genbank_rogue_audit <- function(tree_audit_df, fasta_audit_df, out_csv) {
  parts <- Filter(function(d) is.data.frame(d) && nrow(d) > 0L, list(tree_audit_df, fasta_audit_df))
  comb <- if (length(parts) > 0L) do.call(rbind, parts) else .empty_rogue_audit_df()
  req_cols <- names(.empty_rogue_audit_df()); if (!all(req_cols %in% names(comb))) comb <- .empty_rogue_audit_df()
  comb <- comb[, req_cols, drop = FALSE]; rownames(comb) <- NULL; .write_audit_csv_safe(comb, out_csv); comb
}
