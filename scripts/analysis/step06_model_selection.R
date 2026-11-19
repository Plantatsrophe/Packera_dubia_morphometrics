# Step 06 Evolutionary Model Selection & Hypothesis Testing Utilities
# Sources into 06_genomic_phylogenomics.R; preserves schema 20261127 (0 DDL changes).

.restore_rng_state <- function(had_seed, old_seed) {
  if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
}

.write_csv_safe <- function(df, path) {
  if (is.null(path) || !nzchar(trimws(as.character(path)))) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (requireNamespace("readr", quietly = TRUE)) readr::write_csv(df, path, na = "") else utils::write.csv(df, path, row.names = FALSE)
}

.clean_dosage_matrix <- function(dos_mat) {
  X <- as.matrix(dos_mat); storage.mode(X) <- "numeric"
  if (nrow(X) == 0L || ncol(X) == 0L) return(X)
  fin <- is.finite(X)
  if (any(fin)) {
    mx <- max(X[fin])
    if (mx > 1.0 + 1e-6) X[fin] <- if (mx <= 2.0 + 1e-6) X[fin] / 2.0 else pmin(1.0, pmax(0.0, X[fin]))
    X[fin] <- pmin(1.0, pmax(0.0, X[fin]))
  }
  keep <- apply(X, 2L, function(v) { fv <- v[is.finite(v)]; length(fv) >= 2L && stats::var(fv) > 1e-12 })
  X[, keep, drop = FALSE]
}

.impute_col_means <- function(X_na) {
  col_mu <- colMeans(X_na, na.rm = TRUE); glob_mu <- mean(X_na, na.rm = TRUE)
  col_mu[!is.finite(col_mu)] <- if (is.finite(glob_mu)) glob_mu else 0.5
  col_mu <- pmin(1.0 - 1e-4, pmax(1e-4, col_mu))
  X_imp <- X_na; na_idx <- which(!is.finite(X_imp), arr.ind = TRUE)
  if (nrow(na_idx) > 0L) X_imp[na_idx] <- col_mu[na_idx[, 2L]]
  list(X_imp = X_imp, col_mu = col_mu, obs_mask = is.finite(X_na))
}

.simplex_rows <- function(Q_raw) {
  n <- nrow(Q_raw); k <- ncol(Q_raw)
  if (k <= 1L) return(matrix(1.0, nrow = n, ncol = 1L, dimnames = list(rownames(Q_raw), colnames(Q_raw))))
  out <- matrix(0.0, nrow = n, ncol = k, dimnames = dimnames(Q_raw)); ind <- seq_len(k)
  for (i in seq_len(n)) {
    v <- Q_raw[i, ]; v[!is.finite(v)] <- 0.0; u <- sort(v, decreasing = TRUE); cssv <- cumsum(u)
    cond <- (u - (cssv - 1.0) / ind) > 0.0
    rho <- if (any(cond)) max(ind[cond]) else 1L
    w <- pmin(1.0 - 1e-4, pmax(1e-4, pmax(v - (cssv[rho] - 1.0) / rho, 0.0)))
    out[i, ] <- w / sum(w)
  }
  out
}

.build_speckled_folds <- function(X, n_folds = 5L) {
  n <- nrow(X); nf <- max(1L, as.integer(n_folds)); n_buckets <- max(nf, 10L)
  folds <- vector("list", nf)
  for (i in seq_len(n)) {
    cols <- which(is.finite(X[i, ])); mi <- length(cols); if (mi < 2L) next
    perm <- cols[sample.int(mi)]; b_id <- ((seq_len(mi) - 1L) %% n_buckets) + 1L
    for (f in seq_len(nf)) {
      sel <- perm[b_id == f]
      if (length(sel) > 0L) folds[[f]] <- c(folds[[f]], i + (sel - 1L) * n)
    }
  }
  folds
}

.fit_projected_nmf <- function(X_tr_na, k_val, n_iter = 35L) {
  n <- nrow(X_tr_na); m <- ncol(X_tr_na); imp <- .impute_col_means(X_tr_na)
  if (k_val <= 1L) {
    return(list(Q = matrix(1.0, nrow = n, ncol = 1L, dimnames = list(rownames(X_tr_na), "Ancestry_K1")),
                F = matrix(imp$col_mu, nrow = 1L, ncol = m)))
  }
  X_imp <- imp$X_imp; obs_mask <- imp$obs_mask
  km <- if (nrow(unique(X_imp)) >= k_val) suppressWarnings(tryCatch(stats::kmeans(X_imp, centers = k_val, nstart = 5L, iter.max = 25L), error = function(e) NULL)) else NULL
  Q0 <- if (!is.null(km)) {
    qi <- matrix(0.0, n, k_val); qi[cbind(seq_len(n), km$cluster)] <- 1.0
    qi + matrix(stats::runif(n * k_val, 0.0, 0.05), n, k_val)
  } else matrix(stats::runif(n * k_val, 0.05, 1.0), n, k_val)
  Q <- .simplex_rows(Q0); reg_k <- 1e-4 * diag(k_val)
  F_mat <- matrix(imp$col_mu, nrow = k_val, ncol = m, byrow = TRUE)
  for (it in seq_len(n_iter)) {
    F_raw <- tryCatch(solve(crossprod(Q) + reg_k, crossprod(Q, X_imp)), error = function(e) F_mat)
    F_mat[] <- pmin(pmax(as.numeric(F_raw), 1e-4), 1.0 - 1e-4)
    Q_new <- .simplex_rows(tryCatch(X_imp %*% t(F_mat) %*% solve(tcrossprod(F_mat) + reg_k), error = function(e) Q))
    if (any(!obs_mask)) X_imp[!obs_mask] <- pmin(pmax((Q_new %*% F_mat)[!obs_mask], 0.0), 1.0)
    if (it >= 8L && max(abs(Q_new - Q)) < 1e-5) { Q <- Q_new; break }
    Q <- Q_new
  }
  colnames(Q) <- paste0("Ancestry_K", seq_len(k_val)); rownames(Q) <- rownames(X_tr_na)
  list(Q = Q, F = F_mat)
}

.compute_pam_silhouette <- function(d_mat, k_val) {
  if (k_val < 2L || is.null(d_mat) || k_val >= nrow(d_mat)) return(NA_real_)
  dm <- as.matrix(d_mat); storage.mode(dm) <- "numeric"; if (any(!is.finite(dm))) return(NA_real_)
  if (requireNamespace("cluster", quietly = TRUE)) {
    return(as.numeric(tryCatch(cluster::pam(stats::as.dist(dm), k = k_val, diss = TRUE)$silinfo$avg.width, error = function(e) NA_real_)))
  }
  cl <- stats::cutree(stats::hclust(stats::as.dist(dm), method = "average"), k = k_val); n <- nrow(dm); sil <- numeric(n)
  for (i in seq_len(n)) {
    in_c <- which(cl == cl[i] & seq_len(n) != i); a_i <- if (length(in_c) > 0L) mean(dm[i, in_c]) else 0.0
    oks <- setdiff(unique(cl), cl[i]); b_i <- if (length(oks) > 0L) min(vapply(oks, function(ok) mean(dm[i, cl == ok]), numeric(1L))) else 0.0
    denom <- max(a_i, b_i); sil[i] <- if (length(in_c) == 0L || denom <= 1e-12) 0.0 else (b_i - a_i) / denom
  }
  mean(sil)
}

.try_lea_snmf_ce <- function(X, k_grid, d_mat, seed) {
  if (!requireNamespace("LEA", quietly = TRUE) || nrow(X) < 12L || ncol(X) < 50L ||
      identical(Sys.getenv("PLANTMATICS_FORCE_NMF_CE"), "1")) return(NULL)
  lea_dir <- tempfile("plantmatics_snmf_"); dir.create(lea_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(lea_dir, recursive = TRUE, force = TRUE), add = TRUE)
  geno_p <- file.path(lea_dir, "dos.geno"); G_int <- matrix(9L, nrow = ncol(X), ncol = nrow(X))
  Xt <- t(X); fin <- is.finite(Xt); G_int[fin] <- as.integer(pmin(2L, pmax(0L, round(2.0 * Xt[fin]))))
  writeLines(apply(G_int, 1L, paste0, collapse = ""), geno_p)
  n_cpu <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers() else 1L
  tryCatch({
    set.seed(as.integer(seed))
    proj <- LEA::snmf(geno_p, K = k_grid, entropy = TRUE, repetitions = 3L, project = "new", CPU = n_cpu)
    ce_means <- vapply(k_grid, function(kv) mean(LEA::cross.entropy(proj, K = kv), na.rm = TRUE), numeric(1L))
    ce_ses <- vapply(k_grid, function(kv) { v <- LEA::cross.entropy(proj, K = kv); if (length(v) >= 2L) stats::sd(v) / sqrt(length(v)) else 0.0 }, numeric(1L))
    if (any(!is.finite(ce_means))) return(NULL)
    idx_min <- which.min(ce_means); best_k <- min(k_grid[ce_means <= (ce_means[idx_min] + ce_ses[idx_min])])
    best_run <- which.min(LEA::cross.entropy(proj, K = best_k)); if (length(best_run) == 0L) best_run <- 1L
    Q_opt <- as.matrix(LEA::Q(proj, K = best_k, run = best_run))
    rownames(Q_opt) <- rownames(X); colnames(Q_opt) <- paste0("Ancestry_K", seq_len(best_k))
    list(best_k = as.integer(best_k), Q_opt = Q_opt, ce_means = ce_means,
         rmse_means = rep(NA_real_, length(k_grid)), sils = vapply(k_grid, function(kv) .compute_pam_silhouette(d_mat, kv), numeric(1L)))
  }, error = function(e) NULL)
}

select_admixture_k_cross_entropy <- function(dos_mat, d_mat = NULL, k_max = 8L, n_folds = 5L, seed = 202611L, out_csv = NULL) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(.restore_rng_state(had_seed, old_seed), add = TRUE); set.seed(as.integer(seed))

  X <- .clean_dosage_matrix(dos_mat); n <- nrow(X); m <- ncol(X)
  if (n < 3L || m < 2L) {
    Q_opt <- matrix(1.0, nrow = n, ncol = 1L, dimnames = list(rownames(X), "Ancestry_K1"))
    tbl <- data.frame(K = 1L, masked_cross_entropy = NA_real_, masked_rmse = NA_real_,
                      pam_silhouette = NA_real_, delta_ce_from_min = 0.0, is_selected_k = TRUE, stringsAsFactors = FALSE)
    .write_csv_safe(tbl, out_csv)
    return(list(best_k = 1L, Q_opt = Q_opt, table_df = tbl, selected_k_cross_entropy = NA_real_))
  }
  k_grid <- seq_len(max(1L, min(as.integer(k_max), n - 1L)))
  lea_res <- .try_lea_snmf_ce(X, k_grid, d_mat, seed)
  if (!is.null(lea_res)) {
    best_k <- lea_res$best_k; Q_opt <- lea_res$Q_opt; ce_means <- lea_res$ce_means; rmse_means <- lea_res$rmse_means; sils <- lea_res$sils
  } else {
    folds <- .build_speckled_folds(X, n_folds = n_folds)
    eval_k <- function(kv) {
      set.seed(as.integer(seed) + as.integer(kv))
      ce_f <- rep(NA_real_, length(folds)); rmse_f <- rep(NA_real_, length(folds))
      for (f in seq_along(folds)) {
        m_idx <- folds[[f]]; if (length(m_idx) == 0L) next
        X_tr <- X; X_tr[m_idx] <- NA_real_; fit <- .fit_projected_nmf(X_tr, kv, n_iter = 35L)
        P_hat <- pmin(1.0 - 1e-5, pmax(1e-5, (fit$Q %*% fit$F)[m_idx])); x_obs <- X[m_idx]
        ce_f[f] <- -mean(x_obs * log(P_hat) + (1.0 - x_obs) * log(1.0 - P_hat)); rmse_f[f] <- sqrt(mean((x_obs - P_hat)^2))
      }
      vc <- is.finite(ce_f)
      list(ce_mean = if (any(vc)) mean(ce_f[vc]) else NA_real_, ce_se = if (sum(vc) >= 2L) stats::sd(ce_f[vc]) / sqrt(sum(vc)) else 0.0,
           rmse_mean = if (any(is.finite(rmse_f))) mean(rmse_f, na.rm = TRUE) else NA_real_, sil = .compute_pam_silhouette(d_mat, kv))
    }
    nw <- if ((n * m) > 1000L && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, length(k_grid)) else 1L
    k_res <- if (nw > 1L && exists("packera_mclapply", mode = "function")) packera_mclapply(k_grid, eval_k, mc.cores = nw) else lapply(k_grid, eval_k)
    ce_means <- vapply(k_res, `[[`, numeric(1L), "ce_mean"); ce_ses <- vapply(k_res, `[[`, numeric(1L), "ce_se")
    rmse_means <- vapply(k_res, `[[`, numeric(1L), "rmse_mean"); sils <- vapply(k_res, `[[`, numeric(1L), "sil")
    best_k <- if (all(!is.finite(ce_means))) 1L else {
      min_i <- which.min(ce_means)
      min(k_grid[is.finite(ce_means) & ce_means <= (ce_means[min_i] + ifelse(is.finite(ce_ses[min_i]), ce_ses[min_i], 0.0))])
    }
    set.seed(as.integer(seed) + as.integer(best_k) + 999L)
    Q_opt <- .fit_projected_nmf(X, best_k, n_iter = 35L)$Q
  }
  min_ce <- if (any(is.finite(ce_means))) min(ce_means, na.rm = TRUE) else NA_real_
  tbl <- data.frame(K = as.integer(k_grid), masked_cross_entropy = round(ce_means, 6), masked_rmse = round(rmse_means, 6),
                    pam_silhouette = round(sils, 4), delta_ce_from_min = round(ce_means - min_ce, 6),
                    is_selected_k = (k_grid == as.integer(best_k)), stringsAsFactors = FALSE)
  .write_csv_safe(tbl, out_csv)
  list(best_k = as.integer(best_k), Q_opt = Q_opt, table_df = tbl,
       selected_k_cross_entropy = as.numeric(ce_means[which(k_grid == as.integer(best_k))[1L]]))
}

evaluate_q_reconstruction_ce <- function(dos_mat, q_mat, n_folds = 5L, seed = 202611L) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(.restore_rng_state(had_seed, old_seed), add = TRUE); set.seed(as.integer(seed))
  X <- .clean_dosage_matrix(dos_mat); Q <- as.matrix(q_mat); storage.mode(Q) <- "numeric"
  if (!is.null(rownames(X)) && !is.null(rownames(Q)) && all(rownames(X) %in% rownames(Q))) Q <- Q[rownames(X), , drop = FALSE]
  if (nrow(X) < 3L || ncol(X) < 2L || nrow(Q) != nrow(X) || ncol(Q) < 1L) return(NA_real_)
  Q[!is.finite(Q)] <- 0.0; Q <- .simplex_rows(Q)
  folds <- .build_speckled_folds(X, n_folds = n_folds); reg_k <- 1e-4 * diag(ncol(Q)); ce_f <- rep(NA_real_, length(folds))
  for (f in seq_along(folds)) {
    m_idx <- folds[[f]]; if (length(m_idx) == 0L) next
    X_tr <- X; X_tr[m_idx] <- NA_real_; X_imp <- .impute_col_means(X_tr)$X_imp
    F_fit <- tryCatch(solve(crossprod(Q) + reg_k, crossprod(Q, X_imp)), error = function(e) NULL)
    if (is.null(F_fit)) next
    F_fit[] <- pmin(pmax(as.numeric(F_fit), 1e-4), 1.0 - 1e-4)
    P_hat <- pmin(pmax((Q %*% F_fit)[m_idx], 1e-5), 1.0 - 1e-5); x_obs <- X[m_idx]
    ce_f[f] <- -mean(x_obs * log(P_hat) + (1.0 - x_obs) * log(1.0 - P_hat))
  }
  if (any(is.finite(ce_f))) mean(ce_f, na.rm = TRUE) else NA_real_
}

estimate_ledoit_wolf_lambda <- function(X_centered) {
  X <- as.matrix(X_centered); storage.mode(X) <- "numeric"; n <- nrow(X); p <- ncol(X)
  if (is.null(n) || is.null(p) || n < 2L || p < 1L) return(0.95)
  if (p == 1L) return(0.05)
  X <- scale(X, center = TRUE, scale = FALSE); X[!is.finite(X)] <- 0.0
  G <- tcrossprod(X); tr_S <- sum(diag(G)) / n; mu <- tr_S / p; tr_S2 <- sum(G * G) / (n * n)
  d2 <- (tr_S2 / p) - (mu * mu); if (!is.finite(d2) || d2 <= 1e-12) return(0.95)
  b2_raw <- (sum(diag(G)^2) - n * tr_S2) / (p * (n^2))
  as.numeric(min(0.95, max(0.01, min(d2, max(0.0, b2_raw)) / d2)))
}

.parse_alignment_dnabin <- function(path, target_tips) {
  lines <- readLines(path, warn = FALSE); lines <- lines[nzchar(trimws(lines))]; if (length(lines) < 4L) return(NULL)
  seqs <- list()
  if (startsWith(trimws(lines[1L]), ">")) {
    cur <- ""
    for (ln in lines) {
      tln <- trimws(ln)
      if (startsWith(tln, ">")) { cur <- sub("[[:space:]].*$", "", sub("^>", "", tln)); seqs[[cur]] <- character(0) }
      else if (nzchar(cur)) seqs[[cur]] <- c(seqs[[cur]], gsub("[[:space:]]+", "", tln))
    }
    seqs <- lapply(seqs, paste0, collapse = "")
  } else {
    hdr <- strsplit(trimws(lines[1L]), "[[:space:]]+")[[1L]]; if (length(hdr) < 2L || any(!grepl("^[0-9]+$", hdr[1:2]))) return(NULL)
    for (ln in lines[-1L]) {
      toks <- strsplit(trimws(ln), "[[:space:]]+")[[1L]]
      if (length(toks) >= 2L && (toks[1L] %in% target_tips || grepl("^[A-Za-z0-9_.-]+$", toks[1L]))) {
        seqs[[toks[1L]]] <- paste0(c(seqs[[toks[1L]]], toks[-1L]), collapse = "")
      }
    }
  }
  keep <- intersect(names(seqs), target_tips); if (length(keep) < 4L) return(NULL)
  lens <- vapply(seqs[keep], nchar, integer(1L)); if (length(unique(lens)) != 1L || lens[1L] < 4L) return(NULL)
  cmat <- do.call(rbind, lapply(seqs[keep], function(s) strsplit(tolower(s), "")[[1L]]))
  rownames(cmat) <- keep; ape::as.DNAbin(cmat)
}

.compute_base_clade_cohesion <- function(tree_obj, pass_df) {
  empty <- data.frame(species = character(0), n_samples = integer(0), is_monophyletic = logical(0),
                      mean_cophenetic_dist = numeric(0), median_branch_length = numeric(0),
                      tree_status = character(0), stringsAsFactors = FALSE)
  if (is.null(tree_obj) || !requireNamespace("ape", quietly = TRUE)) return(empty)
  d_coph <- ape::cophenetic.phylo(tree_obj)
  rows <- lapply(sort(unique(as.character(pass_df$species_standardized))), function(tx) {
    ids <- intersect(pass_df$genetic_sample_id[pass_df$species_standardized == tx], tree_obj$tip.label)
    if (length(ids) < 2L) return(NULL)
    sub_d <- d_coph[ids, ids, drop = FALSE]; tip_edges <- which(tree_obj$edge[, 2L] %in% match(ids, tree_obj$tip.label))
    data.frame(species = tx, n_samples = length(ids),
               is_monophyletic = tryCatch(ape::is.monophyletic(tree_obj, ids), error = function(e) NA),
               mean_cophenetic_dist = round(mean(sub_d[upper.tri(sub_d)], na.rm = TRUE), 6),
               median_branch_length = round(if (!is.null(tree_obj$edge.length) && length(tip_edges) > 0L) stats::median(tree_obj$edge.length[tip_edges], na.rm = TRUE) else NA_real_, 6),
               tree_status = "observed_tree", stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows); if (length(rows) == 0L) empty else do.call(rbind, rows)
}

.fit_dna_model_and_boot <- function(dna_mat, tr_sub) {
  tr_sub$node.label <- NULL
  tr_sub$edge.length <- if (is.null(tr_sub$edge.length)) rep(0.01, nrow(tr_sub$edge)) else pmax(1e-4, ifelse(is.finite(tr_sub$edge.length), tr_sub$edge.length, 1e-4))
  if (requireNamespace("phangorn", quietly = TRUE)) {
    pr <- tryCatch({
      pd <- phangorn::phyDat(dna_mat, type = "DNA"); f0 <- phangorn::pml(tr_sub, data = pd, k = 4L); ctrl <- phangorn::pml.control(trace = 0)
      fits <- list(JC = phangorn::optim.pml(f0, model = "JC", rearrangement = "none", control = ctrl),
                   K80 = phangorn::optim.pml(f0, model = "K80", rearrangement = "none", control = ctrl),
                   HKY = phangorn::optim.pml(f0, model = "HKY", rearrangement = "none", control = ctrl),
                   "GTR+G" = phangorn::optim.pml(f0, model = "GTR", optGamma = TRUE, rearrangement = "none", control = ctrl))
      bics <- vapply(fits, function(f) { ll <- stats::logLik(f); dk <- attr(ll, "df"); if (is.null(dk)) dk <- f$df; -2 * as.numeric(ll) + as.numeric(dk) * log(max(1, ncol(dna_mat))) }, numeric(1L))
      bm <- names(which.min(bics))
      nw_bs <- if (identical(.Platform$OS.type, "unix") && exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(64L, 100L) else 1L
      bs <- tryCatch(phangorn::bootstrap.pml(fits[[bm]], bs = 100L, optNni = FALSE, multicore = (nw_bs > 1L), mc.cores = nw_bs, control = ctrl),
                     error = function(e) phangorn::bootstrap.pml(fits[[bm]], bs = 100L, optNni = FALSE, multicore = FALSE, control = ctrl))
      list(model = bm, boot = round(mean(as.numeric(ape::prop.clades(tr_sub, lapply(bs, ape::unroot))), na.rm = TRUE), 2))
    }, error = function(e) NULL)
    if (!is.null(pr) && is.finite(pr$boot)) return(pr)
  }
  mods <- c(JC69 = 0L, K80 = 1L, F84 = 4L, TN93 = 5L)
  bics <- vapply(names(mods), function(m_nm) {
    dm <- suppressWarnings(tryCatch(as.matrix(ape::dist.dna(dna_mat, model = m_nm, pairwise.deletion = TRUE)), error = function(e) NULL))
    if (is.null(dm) || any(!is.finite(dm))) return(Inf)
    nj_m <- tryCatch(ape::nj(stats::as.dist(dm)), error = function(e) NULL); if (is.null(nj_m)) return(Inf)
    dp <- ape::cophenetic.phylo(nj_m)[rownames(dm), colnames(dm)]; lt <- lower.tri(dm); np <- sum(lt)
    np * log(max(1e-12, sum((dm[lt] - dp[lt])^2)) / np) + (nrow(nj_m$edge) + mods[[m_nm]]) * log(np)
  }, numeric(1L))
  bm <- if (any(is.finite(bics))) names(which.min(bics)) else "JC69"
  bfn <- function(xb) {
    dm <- suppressWarnings(as.matrix(ape::dist.dna(xb, model = bm, pairwise.deletion = TRUE)))
    if (any(!is.finite(dm))) { fv <- dm[is.finite(dm) & dm > 0]; dm[!is.finite(dm)] <- if (length(fv) > 0L) max(fv) else 0.05; diag(dm) <- 0 }
    ape::nj(dm)
  }
  bp <- tryCatch(ape::boot.phylo(tr_sub, dna_mat, FUN = bfn, B = 100L, quiet = TRUE), error = function(e) numeric(0))
  list(model = bm, boot = if (any(is.finite(bp))) round(mean(as.numeric(bp), na.rm = TRUE), 2) else NA_real_)
}

.compute_backbone_rf <- function(all_files, tree_obj, pass_df) {
  bb_cand <- all_files[grepl("^public_packera_MO_Astral_lpp\\.tre$", basename(all_files), ignore.case = TRUE)]
  if (length(bb_cand) == 0L) bb_cand <- all_files[grepl("public_nuclear_backbone/.*\\.tre$", all_files, ignore.case = TRUE)]
  if (length(bb_cand) == 0L) return(NA_real_)
  pub_tr <- tryCatch(ape::read.tree(bb_cand[1L]), error = function(e) NULL)
  if (is.null(pub_tr) || length(pub_tr$tip.label) < 4L) return(NA_real_)
  norm_fn <- if (exists(".normalize_public_packera_taxon", mode = "function")) .normalize_public_packera_taxon else function(x) trimws(gsub("_+", " ", x))
  p_lbl <- norm_fn(pub_tr$tip.label); pk <- which(nzchar(p_lbl) & !duplicated(p_lbl)); if (length(pk) < 4L) return(NA_real_)
  pub_tr <- ape::keep.tip(pub_tr, pub_tr$tip.label[pk]); pub_tr$tip.label <- p_lbl[pk]
  sp_map <- stats::setNames(as.character(pass_df$species_standardized), as.character(pass_df$genetic_sample_id))
  qry_v <- if ("is_putative_novel_or_hybrid_query" %in% names(pass_df)) as.integer(pass_df$is_putative_novel_or_hybrid_query) == 1L else rep(FALSE, nrow(pass_df))
  qry_map <- stats::setNames(qry_v, as.character(pass_df$genetic_sample_id))
  keep_s <- tree_obj$tip.label[tree_obj$tip.label %in% names(sp_map) & !is.na(sp_map[tree_obj$tip.label]) & nzchar(sp_map[tree_obj$tip.label]) & (is.na(qry_map[tree_obj$tip.label]) | !qry_map[tree_obj$tip.label])]
  sp_rep <- norm_fn(sp_map[keep_s]); sk <- which(nzchar(sp_rep) & !duplicated(sp_rep)); if (length(sk) < 4L) return(NA_real_)
  sp_tr <- ape::keep.tip(tree_obj, keep_s[sk]); sp_tr$tip.label <- sp_rep[sk]
  sh_sp <- intersect(sp_tr$tip.label, pub_tr$tip.label); if (length(sh_sp) < 4L) return(NA_real_)
  t1 <- ape::unroot(ape::keep.tip(sp_tr, sh_sp)); t2 <- ape::unroot(ape::keep.tip(pub_tr, sh_sp))
  round(min(1.0, max(0.0, as.numeric(ape::dist.topo(t1, t2, method = "PH85")) / max(1.0, as.numeric((t1$Nnode - 1L) + (t2$Nnode - 1L))))), 4)
}

evaluate_ml_tree_model_and_support <- function(raw_geno_dir, tree_obj, pass_df, out_dir = NULL) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit(.restore_rng_state(had_seed, old_seed), add = TRUE); set.seed(202611L)

  sum_df <- if (exists(".evaluate_tree_clade_cohesion", mode = "function")) .evaluate_tree_clade_cohesion(tree_obj, pass_df) else .compute_base_clade_cohesion(tree_obj, pass_df)
  best_model <- NA_character_; mean_boot <- NA_real_; rf_dist <- NA_real_
  if (!is.null(tree_obj) && length(tree_obj$tip.label) >= 4L && requireNamespace("ape", quietly = TRUE) &&
      !is.null(raw_geno_dir) && nzchar(raw_geno_dir) && dir.exists(raw_geno_dir)) {
    all_f <- sort(list.files(raw_geno_dir, full.names = TRUE, recursive = TRUE))
    aln_f <- all_f[grepl("\\.(phy|phylip|fa|fas|fasta|aln|fna)$", basename(all_f), ignore.case = TRUE) & !grepl("^public_packera_", basename(all_f), ignore.case = TRUE)]
    dna_mat <- NULL
    for (af in aln_f) { dna_mat <- tryCatch(.parse_alignment_dnabin(af, tree_obj$tip.label), error = function(e) NULL); if (!is.null(dna_mat)) break }
    if (!is.null(dna_mat) && nrow(dna_mat) >= 4L) {
      res <- .fit_dna_model_and_boot(dna_mat, ape::unroot(ape::keep.tip(tree_obj, rownames(dna_mat))))
      best_model <- res$model; mean_boot <- res$boot
    }
    rf_dist <- .compute_backbone_rf(all_f, tree_obj, pass_df)
  }
  sum_df$best_tree_model <- rep(best_model, nrow(sum_df))
  sum_df$mean_clade_bootstrap <- rep(mean_boot, nrow(sum_df))
  sum_df$rf_dist_public_backbone <- rep(rf_dist, nrow(sum_df))
  if (!is.null(out_dir) && nzchar(trimws(as.character(out_dir)))) .write_csv_safe(sum_df, file.path(out_dir, "genomic_clade_cohesion_summary.csv"))
  sum_df
}
