# ==============================================================================
# Module: scripts/analysis/step06_public_backbone.R
# Project: Plantmatics
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
# Purpose: Public Packera nuclear species-tree backbone (MO, ASTRAL-Pro, wASTRAL),
#   PhyTop ILS/IH reticulation tables, and 11-mer organellar/ITS/waxy placement.
# ==============================================================================

if (!exists("screen_public_tree_rogue_tips", mode = "function") || !exists(".screen_public_backbone_trees_bundle", mode = "function")) {
  .qc_srcs <- c("scripts/analysis/step06_rogue_sequence_qc.R", "step06_rogue_sequence_qc.R",
                file.path(dirname(as.character(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
                                                 unlist(lapply(sys.frames(), function(fr) fr$ofile))))), "step06_rogue_sequence_qc.R"))
  for (.qc_p in .qc_srcs) if (!is.na(.qc_p) && nzchar(.qc_p) && file.exists(.qc_p)) { source(.qc_p); break }
}

.write_csv_safe <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(df, path)
  else utils::write.csv(df, path, row.names = FALSE)
}

.read_delim_safe <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  is_tsv <- grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)
  tryCatch({
    if (requireNamespace("readr", quietly = TRUE)) {
      if (is_tsv) readr::read_tsv(path, show_col_types = FALSE) else readr::read_csv(path, show_col_types = FALSE)
    } else {
      utils::read.delim(path, sep = if (is_tsv) "\t" else ",", stringsAsFactors = FALSE, check.names = FALSE)
    }
  }, error = function(e) NULL)
}

.normalize_public_packera_taxon <- function(raw) {
  if (is.null(raw) || length(raw) == 0L) return(character(0))
  vapply(as.character(raw), function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return("")
    s <- trimws(gsub("[\"']", "", sub("^>", "", trimws(x))))
    if (grepl("\\|", s)) {
      parts <- trimws(strsplit(s, "\\|", fixed = FALSE)[[1L]])
      parts <- parts[nzchar(parts)]
      hit <- parts[grepl("^(Packera|Senecio|P\\.)([_[:space:]]|$)", parts, ignore.case = TRUE)]
      s <- if (length(hit) > 0L) hit[1L] else if (length(parts) >= 2L) parts[2L] else parts[1L]
    }
    s <- trimws(gsub("[[:space:]]+", " ", gsub("_+", " ", s)))
    s <- gsub("(?i)\\b(Packera|Senecio|P\\.)\\s+(?:x\\s+|×\\s*)memmingeri\\b", "\\1 xmemmingeri", s, perl = TRUE)
    s <- sub("^(?i)(?:Packera|P\\.)\\s+sp\\.?(?=\\s+(?:Packera|Senecio|P\\.)\\s+[A-Za-z])\\s+", "", s, perl = TRUE)
    s <- sub("^([A-Za-z]+)\\s+[0-9]+([A-Za-z][A-Za-z-]*)$", "\\1 \\2", s, perl = TRUE)
    if (grepl("(?i)\\b(Packera|P\\.)\\s+sp\\.?\\s*nov\\b", s, perl = TRUE)) return("Packera sp. nov.")
    if (grepl("(?i)^(?:sp\\.?|cf\\.?|aff\\.?|unk|unknown)$", s, perl = TRUE)) return("")
    m <- regexec("(?i)\\b(Packera|Senecio|P\\.)\\s+[0-9]*([A-Za-z][A-Za-z0-9.-]*)", s, perl = TRUE)
    reg <- regmatches(s, m)[[1L]]
    if (length(reg) >= 3L) {
      gen_low <- tolower(reg[2L])
      ep_raw <- sub("[-.].*$", "", reg[3L])
      ep_low <- tolower(ep_raw)
      if (ep_low %in% c("sp", "cf", "aff", "unk", "unknown")) return("")
      sen_map <- c(
        tomentosus = "dubia", tomentosa = "dubia", aureus = "aurea",
        anonymus = "anonyma", obovatus = "obovata", pauperculus = "paupercula",
        schweinitzianus = "schweinitziana", millefolium = "millefolium", millefolia = "millefolium",
        antennariifolius = "antennariifolia", glabellus = "glabella",
        plattensis = "plattensis", pseudaureus = "pseudaurea",
        indecorus = "indecora", canus = "cana", loratifolius = "loratifolia",
        strepthanthifolius = "streptanthifolia", streptanthifolius = "streptanthifolia",
        fendleri = "fendleri", multilobatus = "multilobata", neomexicanus = "neomexicana",
        debilis = "debilis", werneriifolius = "werneriifolia", hartianus = "hartiana",
        sanguisorboides = "sanguisorboides", dimorphophyllus = "dimorphophylla",
        crocanthes = "crocanthes", hyperborealis = "hyperborealis", porteri = "porteri",
        bellidifolius = "bellidifolia", ganderi = "ganderi", bernardinus = "bernardina",
        hesperius = "hesperia", macounii = "macounii", cardamine = "cardamine",
        cynthioides = "cynthioides", cymbalaria = "cymbalaria", flettii = "flettii",
        layneae = "layneae", cleavelandii = "clevelandii", eurycephalus = "eurycephala",
        bolanderi = "bolanderi", breweri = "breweri", memmingeri = "xmemmingeri", xmemmingeri = "xmemmingeri"
      )
      if (gen_low == "senecio") {
        if (ep_low %in% names(sen_map)) return(paste("Packera", sen_map[[ep_low]]))
        return(paste("Senecio", ep_low))
      }
      if (ep_low %in% c("tomentosa", "tomentosus")) return("Packera dubia")
      if (ep_low == "millefolia") return("Packera millefolium")
      if (ep_low %in% c("memmingeri", "xmemmingeri")) return("Packera xmemmingeri")
      if (ep_low == "strepthanthifolius") return("Packera streptanthifolia")
      return(paste("Packera", ep_low))
    }
    if (tolower(s) %in% c("tomentosa", "tomentosus")) return("Packera dubia")
    s
  }, character(1L), USE.NAMES = FALSE)
}

.read_clean_public_tree <- function(tre_path) {
  if (!nzchar(tre_path) || !file.exists(tre_path) || !requireNamespace("ape", quietly = TRUE)) return(NULL)
  if (exists(".read_nexus_phylo_qc", mode = "function")) {
    nx_tr <- .read_nexus_phylo_qc(tre_path)
    if (!is.null(nx_tr)) return(nx_tr)
  }
  tr <- tryCatch(ape::read.tree(tre_path), error = function(e) NULL)
  if (is.null(tr) || inherits(tr, "multiPhylo") || length(tr$tip.label) < 3L) return(NULL)
  norm_lbls <- .normalize_public_packera_taxon(tr$tip.label)
  keep_idx <- which(nzchar(norm_lbls) & !duplicated(norm_lbls))
  if (length(keep_idx) < 3L) return(NULL)
  if (length(keep_idx) < length(tr$tip.label)) {
    tmp_lbls <- paste0("__PK_TIP_", seq_along(tr$tip.label), "__")
    tr$tip.label <- tmp_lbls
    tr <- ape::drop.tip(tr, tmp_lbls[-keep_idx])
    tr$tip.label <- norm_lbls[keep_idx][match(tr$tip.label, tmp_lbls[keep_idx])]
  } else {
    tr$tip.label <- norm_lbls
  }
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1.0, nrow(tr$edge))
  tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length < 0] <- 0.0
  storage.mode(tr$edge) <- "integer"
  if (exists(".is_valid_phylo_qc", mode = "function") && !.is_valid_phylo_qc(tr)) return(NULL)
  tr
}

.find_geno_file <- function(dirs, fname) {
  for (d in dirs) {
    if (is.null(d) || !nzchar(d) || !dir.exists(d)) next
    p <- file.path(d, fname)
    if (file.exists(p)) return(p)
    cand <- sort(list.files(d, full.names = TRUE, recursive = TRUE))
    hits <- cand[tolower(basename(cand)) == tolower(fname)]
    if (length(hits) > 0L) return(hits[1L])
  }
  ""
}

load_public_packera_backbone <- function(raw_geno_dir) {
  empty_trees <- list(MO = NULL, ASTRALPro = NULL, wASTRAL = NULL, plastome111 = NULL)
  if (is.null(raw_geno_dir) || !nzchar(raw_geno_dir) || !dir.exists(raw_geno_dir)) {
    return(list(has_backbone = FALSE, trees = empty_trees, reltime_trees = list(),
                phytop_tables = list(), tree_files_found = character(0), tree_rogue_audit_df = .empty_rogue_audit_df()))
  }
  tree_map <- c(
    MO = "public_packera_MO_Astral_lpp.tre",
    ASTRALPro = "public_packera_ASTRALPro_Astral_lpp.tre",
    wASTRAL = "public_packera_wASTRAL_Astral_lpp.tre",
    plastome111 = "public_packera_111mapping_plastomes.tre"
  )
  trees <- empty_trees
  found_files <- character(0)
  for (nm in names(tree_map)) {
    p <- .find_geno_file(raw_geno_dir, tree_map[[nm]])
    if (nzchar(p) && file.exists(p)) {
      found_files <- c(found_files, stats::setNames(p, nm))
      trees[[nm]] <- .read_clean_public_tree(p)
    }
  }
  all_gf <- sort(list.files(raw_geno_dir, full.names = TRUE, recursive = TRUE))
  pt_files <- all_gf[grepl("^public_packera_PhyTop_.*\\.(tsv|csv)$", basename(all_gf), ignore.case = TRUE)]
  phytop_tables <- list()
  for (pf in pt_files) {
    df_pt <- .read_delim_safe(pf)
    if (!is.null(df_pt) && nrow(df_pt) > 0L) {
      df_res <- as.data.frame(df_pt, stringsAsFactors = FALSE)
      tre_p <- sub("\\.(tsv|csv)$", ".tre", pf, ignore.case = TRUE)
      if (file.exists(tre_p) && requireNamespace("ape", quietly = TRUE)) {
        pt_tr <- tryCatch(ape::read.tree(tre_p), error = function(e) NULL)
        if (inherits(pt_tr, "multiPhylo") && length(pt_tr) > 0L) pt_tr <- pt_tr[[1L]]
        if (!is.null(pt_tr)) attr(df_res, "phytop_tree") <- pt_tr
      }
      phytop_tables[[basename(pf)]] <- df_res
    }
  }
  qc <- .screen_public_backbone_trees_bundle(trees, found_files, phytop_tables, all_gf, .read_clean_public_tree)
  has_bb <- any(vapply(qc$trees, function(x) !is.null(x), logical(1L)))
  list(has_backbone = isTRUE(has_bb), trees = qc$trees, reltime_trees = qc$reltime_trees,
       phytop_tables = qc$phytop_tables, tree_files_found = unname(found_files), tree_rogue_audit_df = qc$audit_df)
}

select_backbone_guided_outgroup <- function(uniq_taxa, p_dist, backbone = NULL) {
  if (is.null(backbone) || is.null(backbone$trees) || length(uniq_taxa) < 4L || !requireNamespace("ape", quietly = TRUE)) {
    return(NULL)
  }
  ref_tr <- backbone$trees$MO
  if (is.null(ref_tr)) ref_tr <- backbone$trees$ASTRALPro
  if (is.null(ref_tr)) ref_tr <- backbone$trees$wASTRAL
  if (is.null(ref_tr)) return(NULL)
  norm_u <- .normalize_public_packera_taxon(uniq_taxa)
  norm_to_orig <- character(0)
  for (idx in seq_along(uniq_taxa)) {
    nk <- norm_u[idx]
    if (nzchar(nk) && (!nk %in% names(norm_to_orig) || identical(uniq_taxa[idx], nk))) {
      norm_to_orig[[nk]] <- uniq_taxa[idx]
    }
  }
  ov <- intersect(names(norm_to_orig), ref_tr$tip.label)
  if (length(ov) < 3L) return(NULL)
  pruned <- tryCatch(ape::keep.tip(ref_tr, ov), error = function(e) NULL)
  if (is.null(pruned) || length(pruned$tip.label) < 3L) return(NULL)
  if (is.null(pruned$edge.length)) pruned$edge.length <- rep(1.0, nrow(pruned$edge))
  b_coph <- ape::cophenetic.phylo(pruned)
  m_dist <- rowMeans(b_coph, na.rm = TRUE)
  max_d <- max(m_dist, na.rm = TRUE)
  cands <- names(m_dist)[abs(m_dist - max_d) <= 1e-9]
  if (length(cands) > 1L && any(cands != "Packera dubia")) cands <- cands[cands != "Packera dubia"]
  outgroup <- cands[1L]
  lad <- ape::reorder.phylo(ape::ladderize(pruned), order = "cladewise")
  tip_ord <- lad$tip.label[lad$edge[lad$edge[, 2L] <= length(lad$tip.label), 2L]]
  base_ord <- setdiff(tip_ord, outgroup)
  outgroup <- unname(norm_to_orig[outgroup])
  base_ord <- unname(norm_to_orig[base_ord])
  extra <- setdiff(setdiff(uniq_taxa, outgroup), base_ord)
  for (ex in extra) {
    if (length(base_ord) == 0L) { base_ord <- c(base_ord, ex); next }
    d_vec <- if (ex %in% rownames(p_dist) && all(base_ord %in% colnames(p_dist))) p_dist[ex, base_ord] else rep(1, length(base_ord))
    base_ord <- append(base_ord, ex, after = which.min(d_vec))
  }
  list(outgroup = outgroup, tree_order = base_ord, source = "public_nuclear_backbone")
}

.nearest_and_neighbors <- function(tr, norm_tx, k = 3L) {
  if (is.null(tr) || !(norm_tx %in% tr$tip.label)) return(list(present = FALSE, sister = NA_character_, neigh = character(0)))
  coph <- ape::cophenetic.phylo(tr)
  others <- setdiff(colnames(coph), norm_tx)
  if (length(others) == 0L) return(list(present = TRUE, sister = NA_character_, neigh = character(0)))
  d <- sort(coph[norm_tx, others], na.last = TRUE)
  list(present = TRUE, sister = names(d)[1L], neigh = head(names(d), min(k, length(d))))
}

.match_phytop_taxon_metrics <- function(norm_tx, phytop_tables) {
  if (length(phytop_tables) == 0L || !nzchar(norm_tx)) return(list(ils = NA_real_, ih = NA_real_))
  ils_vals <- numeric(0); ih_vals <- numeric(0)
  for (df in phytop_tables) {
    cn <- names(df); lc <- tolower(cn)
    t_col <- cn[lc %in% c("taxon", "node", "branch", "clade", "species", "lineage", "tip")][1L]
    i_col <- cn[lc %in% c("ils", "ils_index", "ils_score", "phytop_ils")][1L]
    h_col <- cn[lc %in% c("ih", "ih_index", "introgression", "introgression_index", "phytop_ih")][1L]
    if (is.na(t_col) || (is.na(i_col) && is.na(h_col))) next
    keys <- as.character(df[[t_col]])
    hits <- which(vapply(keys, function(kv) {
      if (is.na(kv) || !nzchar(trimws(kv))) return(FALSE)
      toks <- trimws(unlist(strsplit(kv, "[;,+/|]")))
      any(.normalize_public_packera_taxon(toks) == norm_tx)
    }, logical(1L)))
    if (length(hits) == 0L) {
      pt_tr <- attr(df, "phytop_tree")
      if (!is.null(pt_tr) && length(pt_tr$tip.label) >= 3L && length(pt_tr$node.label) > 0L) {
        n_tips <- length(pt_tr$tip.label)
        norm_tips <- .normalize_public_packera_taxon(pt_tr$tip.label)
        tip_idx <- which(norm_tips == norm_tx)
        if (length(tip_idx) > 0L) {
          exact_mask <- grepl(sprintf("^%s$", gsub(" ", "_", norm_tx)),
                              sub("_[0-9]+", "_", pt_tr$tip.label[tip_idx]), ignore.case = TRUE)
          if (any(exact_mask)) tip_idx <- tip_idx[exact_mask]
          matched_labels <- character(0)
          for (ti in tip_idx) {
            cur <- ti
            for (step in seq_len(length(pt_tr$node.label))) {
              par <- pt_tr$edge[pt_tr$edge[, 2L] == cur, 1L]
              if (length(par) == 0L) break
              par <- par[1L]
              lbl <- if (par > n_tips) as.character(pt_tr$node.label[par - n_tips]) else ""
              if (nzchar(lbl) && lbl %in% keys) {
                matched_labels <- c(matched_labels, lbl)
                break
              }
              cur <- par
            }
          }
          if (length(matched_labels) > 0L) hits <- which(keys %in% unique(matched_labels))
        }
      }
    }
    if (length(hits) > 0L) {
      if (!is.na(i_col)) { v <- suppressWarnings(as.numeric(df[[i_col]][hits])); ils_vals <- c(ils_vals, v[is.finite(v)]) }
      if (!is.na(h_col)) { v <- suppressWarnings(as.numeric(df[[h_col]][hits])); ih_vals <- c(ih_vals, v[is.finite(v)]) }
    }
  }
  list(
    ils = if (length(ils_vals) > 0L) round(mean(ils_vals), 5) else NA_real_,
    ih = if (length(ih_vals) > 0L) round(mean(ih_vals), 5) else NA_real_
  )
}

.read_fasta_named <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(data.frame(header = character(0), seq = character(0), stringsAsFactors = FALSE))
  lns <- readLines(path, warn = FALSE)
  hdr_idx <- grep("^>", lns)
  if (length(hdr_idx) == 0L) return(data.frame(header = character(0), seq = character(0), stringsAsFactors = FALSE))
  ends <- c(hdr_idx[-1L] - 1L, length(lns))
  hdrs <- sub("^>", "", trimws(lns[hdr_idx]))
  seqs <- vapply(seq_along(hdr_idx), function(i) {
    if (hdr_idx[i] + 1L > ends[i]) return("")
    gsub("[^A-Z]", "", toupper(paste(lns[(hdr_idx[i] + 1L):ends[i]], collapse = "")))
  }, character(1L))
  keep <- nzchar(hdrs) & nchar(seqs) >= 11L
  data.frame(header = hdrs[keep], seq = seqs[keep], stringsAsFactors = FALSE)
}

.kmers_11 <- function(s) {
  n <- nchar(s)
  if (is.na(n) || n < 11L) return(character(0))
  st <- seq_len(n - 10L)
  k <- substring(s, st, st + 10L)
  unique(k[!grepl("[^ACGT]", k)])
}

.revcomp_seq <- function(s) {
  ch <- rev(strsplit(s, "", fixed = TRUE)[[1L]])
  chartr("ACGTU", "TGCAA", paste(ch, collapse = ""))
}

.place_organellar_loci <- function(raw_geno_dir, out_dir, sample_ids) {
  search_dirs <- unique(Filter(function(d) nzchar(d) && dir.exists(d),
                               c(raw_geno_dir, out_dir, dirname(out_dir), file.path(dirname(out_dir), "genomics"))))
  locus_specs <- list(
    list(locus = "plastid", qfile = "skimmed_plastid_consensus.fasta",
         refs = c("public_packera_111mapping_plastomes.fasta", "public_packera_18complete_plastomes.fasta", "public_packera_plastomes_and_cpdna.fasta")),
    list(locus = "its", qfile = "skimmed_its_consensus.fasta",
         refs = c("public_packera_its_ets.fasta", "public_packera_bold_barcodes.fasta", "bain_jansen_sequences.fasta")),
    list(locus = "waxy", qfile = "skimmed_waxy_consensus.fasta",
         refs = c("public_packera_waxy_nuclear.fasta"))
  )
  rows <- list()
  audit_list <- list()
  for (sp in locus_specs) {
    ref_list <- list()
    for (rf in sp$refs) {
      rp <- .find_geno_file(search_dirs, rf)
      rdf <- .read_fasta_named(rp)
      if (nrow(rdf) > 0L) {
        rdf$taxon <- .normalize_public_packera_taxon(rdf$header)
        rdf <- rdf[!is.na(rdf$taxon) & nzchar(rdf$taxon), , drop = FALSE]
        if (nrow(rdf) > 0L) {
          rdf$ref_source <- rf
          ref_list[[length(ref_list) + 1L]] <- rdf
        }
      }
    }
    if (length(ref_list) == 0L) next
    refs <- do.call(rbind, ref_list)
    r_km <- lapply(refs$seq, .kmers_11)
    qc_fa <- screen_public_fasta_outliers(refs, sp$locus, kmers_list = r_km)
    if (!is.null(qc_fa$audit_df) && nrow(qc_fa$audit_df) > 0L) audit_list[[length(audit_list) + 1L]] <- qc_fa$audit_df
    refs <- qc_fa$clean_refs; r_km <- qc_fa$clean_kmers
    qp <- .find_geno_file(search_dirs, sp$qfile)
    qdf <- .read_fasta_named(qp)
    if (nrow(qdf) == 0L || nrow(refs) == 0L) next
    r_len <- vapply(r_km, length, integer(1L))
    valid_r <- which(r_len > 0L)
    if (length(valid_r) == 0L) next
    r_sk <- lapply(seq_along(r_km), function(ri) if (r_len[ri] > 6000L) r_km[[ri]][endsWith(r_km[[ri]], "AA")] else r_km[[ri]])
    r_sk_len <- vapply(r_sk, length, integer(1L))
    for (qi in seq_len(nrow(qdf))) {
      sid <- sub("\\|.*$", "", sub("\\s.*$", "", qdf$header[qi]))
      if (length(sample_ids) > 0L && !(sid %in% sample_ids)) next
      qf <- .kmers_11(qdf$seq[qi]); qr <- .kmers_11(.revcomp_seq(qdf$seq[qi]))
      q_len <- length(qf)
      if (q_len < 3L) next
      use_q_sk <- q_len > 6000L
      qf_u <- if (use_q_sk) qf[endsWith(qf, "AA")] else qf
      qr_u <- if (use_q_sk) qr[endsWith(qr, "AA")] else qr
      q_eff <- length(qf_u); if (q_eff < 3L) { qf_u <- qf; qr_u <- qr; q_eff <- q_len; use_q_sk <- FALSE }
      best <- list(score = -1, sh = 0L, cont = 0, jac = 0, ri = NA_integer_)
      for (ri in valid_r) {
        rk <- if (use_q_sk) r_sk[[ri]] else r_km[[ri]]; rl <- if (use_q_sk) r_sk_len[ri] else r_len[ri]
        if (rl < 3L) next
        sh <- max(sum(qf_u %in% rk), sum(qr_u %in% rk))
        if (sh < 3L) next
        cont <- sh / q_eff
        if (cont < 0.05) next
        jac <- sh / (q_eff + rl - sh)
        sc <- cont + 1e-4 * jac
        if (sc > best$score) best <- list(score = sc, sh = sh, cont = cont, jac = jac, ri = ri)
      }
      if (!is.na(best$ri)) {
        rows[[length(rows) + 1L]] <- data.frame(
          genetic_sample_id = sid, locus = sp$locus,
          best_match_taxon = refs$taxon[best$ri], best_ref_header = refs$header[best$ri],
          ref_source_file = refs$ref_source[best$ri], n_shared_11mers = as.integer(best$sh),
          kmer_containment = round(best$cont, 4), kmer_jaccard = round(best$jac, 4), stringsAsFactors = FALSE
        )
      }
    }
  }
  res_df <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame(
    genetic_sample_id = character(0), locus = character(0), best_match_taxon = character(0),
    best_ref_header = character(0), ref_source_file = character(0), n_shared_11mers = integer(0),
    kmer_containment = numeric(0), kmer_jaccard = numeric(0), stringsAsFactors = FALSE
  )
  attr(res_df, "fasta_rogue_audit") <- if (length(audit_list) > 0L) do.call(rbind, audit_list) else .empty_rogue_audit_df()
  res_df
}

evaluate_public_backbone_and_organellar <- function(pass_df, raw_geno_dir, out_dir, extra_taxa = character(0)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  bb <- load_public_packera_backbone(raw_geno_dir)
  n <- nrow(pass_df)
  sp_col <- intersect(c("species_standardized", "taxon", "taxon_name"), names(pass_df))[1L]
  sp_raw <- if (n > 0L && !is.na(sp_col)) as.character(pass_df[[sp_col]]) else rep(NA_character_, n)
  ex_clean <- as.character(extra_taxa[!is.na(extra_taxa) & nzchar(trimws(as.character(extra_taxa)))])
  uniq_tx <- unique(c(sp_raw[!is.na(sp_raw) & nzchar(trimws(sp_raw))], ex_clean))

  mo_tr <- bb$trees$MO; ap_tr <- bb$trees$ASTRALPro
  shared_tips <- if (!is.null(mo_tr) && !is.null(ap_tr)) intersect(mo_tr$tip.label, ap_tr$tip.label) else character(0)
  mo_sh <- if (length(shared_tips) >= 4L) ape::keep.tip(mo_tr, shared_tips) else NULL
  ap_sh <- if (length(shared_tips) >= 4L) ape::keep.tip(ap_tr, shared_tips) else NULL

  conc_rows <- lapply(uniq_tx, function(tx) {
    ntx <- .normalize_public_packera_taxon(tx)
    m_info <- .nearest_and_neighbors(mo_tr, ntx, k = 3L)
    a_info <- .nearest_and_neighbors(ap_tr, ntx, k = 3L)
    jacc <- NA_real_; conc_flag <- NA
    if (!is.null(mo_sh) && !is.null(ap_sh) && (ntx %in% shared_tips)) {
      nm <- .nearest_and_neighbors(mo_sh, ntx, k = 3L)$neigh
      na <- .nearest_and_neighbors(ap_sh, ntx, k = 3L)$neigh
      if (length(union(nm, na)) > 0L) {
        jacc <- round(length(intersect(nm, na)) / length(union(nm, na)), 4)
        conc_flag <- isTRUE(jacc >= 0.5)
      }
    }
    pt <- .match_phytop_taxon_metrics(ntx, bb$phytop_tables)
    clade_str <- if (isTRUE(m_info$present) && !is.na(m_info$sister)) {
      sprintf("MO:%s+%s", ntx, m_info$sister)
    } else if (isTRUE(a_info$present) && !is.na(a_info$sister)) {
      sprintf("ASTRALPro:%s+%s", ntx, a_info$sister)
    } else NA_character_
    data.frame(
      taxon = tx, normalized_taxon = ntx,
      in_mo_backbone = isTRUE(m_info$present), in_astralpro_backbone = isTRUE(a_info$present),
      public_backbone_clade = clade_str, mo_nearest_sister = m_info$sister,
      astralpro_nearest_sister = a_info$sister, mo_astralpro_jaccard = jacc,
      mo_astralpro_concordant = conc_flag,
      public_phytop_ils_index = pt$ils, public_phytop_ih_index = pt$ih,
      stringsAsFactors = FALSE
    )
  })
  conc_df <- if (length(conc_rows) > 0L) do.call(rbind, conc_rows) else data.frame(
    taxon = character(0), normalized_taxon = character(0), in_mo_backbone = logical(0),
    in_astralpro_backbone = logical(0), public_backbone_clade = character(0),
    mo_nearest_sister = character(0), astralpro_nearest_sister = character(0),
    mo_astralpro_jaccard = numeric(0), mo_astralpro_concordant = logical(0),
    public_phytop_ils_index = numeric(0), public_phytop_ih_index = numeric(0),
    stringsAsFactors = FALSE
  )
  .write_csv_safe(conc_df, file.path(out_dir, "public_backbone_concordance.csv"))

  sids <- if (n > 0L && "genetic_sample_id" %in% names(pass_df)) as.character(pass_df$genetic_sample_id) else character(0)
  place_df <- .place_organellar_loci(raw_geno_dir, out_dir, sids)
  .write_csv_safe(place_df, file.path(out_dir, "public_organellar_placement.csv"))
  rogue_audit_df <- .write_public_genbank_rogue_audit(
    bb$tree_rogue_audit_df, attr(place_df, "fasta_rogue_audit"),
    file.path(out_dir, "public_genbank_rogue_audit.csv")
  )

  idx_tx <- match(sp_raw, conc_df$taxon)
  locus_vec <- function(loc) {
    out <- rep(NA_character_, n)
    sub <- place_df[place_df$locus == loc, , drop = FALSE]
    if (nrow(sub) > 0L && length(sids) > 0L) {
      m <- match(sids, sub$genetic_sample_id)
      ok <- !is.na(m)
      out[ok] <- sub$best_match_taxon[m[ok]]
    }
    out
  }
  list(
    public_backbone_clade = if (n > 0L) conc_df$public_backbone_clade[idx_tx] else character(0),
    public_phytop_ils_index = if (n > 0L) conc_df$public_phytop_ils_index[idx_tx] else numeric(0),
    public_phytop_ih_index = if (n > 0L) conc_df$public_phytop_ih_index[idx_tx] else numeric(0),
    public_plastome_match = locus_vec("plastid"),
    public_its_match = locus_vec("its"),
    waxy_lineage_match = locus_vec("waxy"),
    backbone = bb, concordance_df = conc_df, placement_df = place_df, rogue_audit_df = rogue_audit_df
  )
}
