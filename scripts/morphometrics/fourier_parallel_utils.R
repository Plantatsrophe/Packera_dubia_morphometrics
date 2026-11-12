# ==============================================================================
# Module: scripts/morphometrics/fourier_parallel_utils.R
# Project: Packera Species Delimitation & Morphometrics Pipeline (NCU)
#
# Description:
#   Multi-threaded Step 03 Elliptic Fourier & OTU median aggregation utilities:
#     - Fast atomic CSV export via data.table::fwrite (fallback utils::write.csv).
#     - Dynamic multi-core chunked Momocs centering, scaling, & homologized
#       efourier(norm = FALSE, start = FALSE) decomposition (bit-identical to single-core).
#     - Unweighted basal-rosette OTU coordinate-wise median harmonic, ECT, and
#       dissection-metric aggregation across (catalogNumber, plant_individual_id)
#       with full Tier_0 nomenclatural type immunity from outlier dropping.
# ==============================================================================

write_step03_csv_fast <- function(df, target_path) {
  dir.create(dirname(target_path), recursive = TRUE, showWarnings = FALSE)
  tmp_p <- paste0(target_path, ".", Sys.getpid(), ".tmp")
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(df, file = tmp_p, na = "", quote = "auto", logical01 = FALSE, showProgress = FALSE)
  } else {
    utils::write.csv(df, file = tmp_p, row.names = FALSE, na = "")
  }
  file.rename(tmp_p, target_path)
}

align_contour_horizontal_xax <- function(coo, is_upright = NA) {
  if (is.null(coo) || !is.matrix(coo) || nrow(coo) < 5L || ncol(coo) < 2L) return(coo)
  if (isTRUE(attr(coo, "lm3_xax_aligned"))) return(coo)
  cx <- mean(coo[, 1L], na.rm = TRUE)
  cy <- mean(coo[, 2L], na.rm = TRUE)
  dx <- coo[, 1L] - cx
  dy <- coo[, 2L] - cy
  span_x <- diff(range(dx, na.rm = TRUE))
  span_y <- diff(range(dy, na.rm = TRUE))
  max_abs_x <- max(abs(dx), 1e-8)
  max_abs_y <- max(abs(dy), 1e-8)

  # 1. Step 02 canonical image-space vertical tip-up contour:
  #    When Step 02 LM3 momocs_prep explicitly certifies mask_is_upright == TRUE,
  #    or when unflagged contours carry the Step 02 tip-up signature (including deeply
  #    cordate/reniform blades whose basal auricles hang up to 55% below the sinus),
  #    apply the exact chart-preserving rigid orthogonal map (X = -dy + cx, Y = -dx + cy)
  #    so lamina_base sits at (-X, 0), lamina_tip at (+X, 0), anatomical left at (+Y),
  #    clockwise in Momocs Cartesian space with 0.000-deg midrib drift.
  is_step02_tipup <- isTRUE(is_upright) || (
    !identical(is_upright, FALSE) && is.finite(dy[1L]) && (dy[1L] > 0.0) && (min(dy, na.rm = TRUE) < 0.0) &&
    ((((max(dy, na.rm = TRUE) - dy[1L]) <= 0.30 * max(span_y, 1e-8)) && (abs(dx[1L]) <= 0.28 * max(span_x, 1e-8))) ||
     ((abs(dx[1L]) <= 0.12 * max(span_x, 1e-8)) && ((max(dy, na.rm = TRUE) - dy[1L]) <= 0.58 * max(span_y, 1e-8))))
  )
  if (isTRUE(is_step02_tipup)) {
    out_coo <- cbind(-dy + cx, -dx + cy)
    colnames(out_coo) <- colnames(coo)
    attr(out_coo, "lm3_xax_aligned") <- TRUE
    attr(out_coo, "lm3_from_tipup") <- TRUE
    return(out_coo)
  }

  # 2. Already horizontally aligned along X-axis (Point 1 at basal/apical X pole):
  #    Preserve exact midrib orientation without vertex-61 perimeter-asymmetry drift.
  if (abs(dx[1L]) >= 0.42 * max_abs_x && abs(dy[1L]) <= 0.22 * max_abs_y) {
    attr(coo, "lm3_xax_aligned") <- TRUE
    return(coo)
  }

  # 3. Fallback for arbitrary oblique legacy contours lacking Step 02 tip-up frame:
  #    Align Point 1 -> apical cap median onto +X without relying on 50% perimeter index.
  d_from_p1 <- sqrt((coo[, 1L] - coo[1L, 1L])^2 + (coo[, 2L] - coo[1L, 2L])^2)
  max_d <- max(d_from_p1, na.rm = TRUE)
  if (!is.finite(max_d) || max_d < 1e-8) {
    attr(coo, "lm3_xax_aligned") <- TRUE
    return(coo)
  }
  mid_tip_idx <- if (is.finite(dy[1L]) && dy[1L] > 0.0 && min(dy, na.rm = TRUE) < 0.0 && abs(dx[1L]) <= 0.18 * max(span_x, 1e-8)) {
    which(abs(dx) <= 0.18 * max(span_x, 1e-8) & dy <= (min(dy, na.rm = TRUE) + 0.08 * max(span_y, 1e-8)))
  } else integer(0L)
  cap_idx <- if (length(mid_tip_idx) >= 1L) mid_tip_idx else which(d_from_p1 >= 0.985 * max_d)
  ap_x <- stats::median(coo[cap_idx, 1L], na.rm = TRUE)
  ap_y <- stats::median(coo[cap_idx, 2L], na.rm = TRUE)
  v_apical <- c(ap_x - coo[1L, 1L], ap_y - coo[1L, 2L])
  v_norm <- sqrt(sum(v_apical^2))
  if (!is.finite(v_norm) || v_norm < 1e-8) {
    attr(coo, "lm3_xax_aligned") <- TRUE
    return(coo)
  }
  theta <- -atan2(v_apical[2L], v_apical[1L])
  cs_t <- cos(theta)
  sn_t <- sin(theta)
  out_coo <- cbind(dx * cs_t - dy * sn_t + cx, dx * sn_t + dy * cs_t + cy)
  colnames(out_coo) <- colnames(coo)
  attr(out_coo, "lm3_xax_aligned") <- TRUE
  out_coo
}

run_parallel_momocs_efa <- function(coo_list, fac_df, nb_h, calibrate_fn) {
  coo_list <- lapply(coo_list, align_contour_horizontal_xax)
  n_out <- length(coo_list)
  n_w <- if (exists("resolve_packera_workers", mode = "function")) {
    resolve_packera_workers(32L, ceiling(n_out / 80L))
  } else {
    1L
  }
  chunk_ids <- split(seq_len(n_out), cut(seq_len(n_out), breaks = max(2L, n_w), labels = FALSE))
  prep_fn <- function(idx_vec) {
    m <- length(idx_vec)
    cs_v <- numeric(m)
    sc_list <- vector("list", m)
    for (k in seq_len(m)) {
      c_i <- coo_list[[idx_vec[k]]]
      cs_v[k] <- tryCatch(Momocs::coo_centsize(c_i), error = function(e) NA_real_)
      sc_list[[k]] <- Momocs::coo_scale(Momocs::coo_center(c_i))
    }
    list(cs = cs_v, coo = sc_list)
  }
  prep_chunks <- if (exists("packera_mclapply", mode = "function")) {
    packera_mclapply(chunk_ids, prep_fn, mc.cores = n_w)
  } else {
    lapply(chunk_ids, prep_fn)
  }
  cs_vals <- unlist(lapply(prep_chunks, `[[`, "cs"), use.names = FALSE)
  scaled_coos <- unlist(lapply(prep_chunks, `[[`, "coo"), recursive = FALSE, use.names = FALSE)
  names(scaled_coos) <- names(coo_list)
  if (length(cs_vals) == nrow(fac_df)) {
    cal_ok <- if ("is_calibrated_mm" %in% names(fac_df)) (!is.na(fac_df$is_calibrated_mm) & as.logical(fac_df$is_calibrated_mm)) else rep(TRUE, nrow(fac_df))
    fac_df$centroid_size <- ifelse(cal_ok, round(as.numeric(cs_vals), 4), NA_real_)
    fac_df$log_centsize <- ifelse(cal_ok & !is.na(fac_df$centroid_size) & fac_df$centroid_size > 0, round(log10(fac_df$centroid_size), 6), NA_real_)
  }
  out_obj <- Momocs::Out(scaled_coos, fac = fac_df)
  h_opt <- calibrate_fn(out_obj, nb_h = nb_h)

  efa_chunk_fn <- function(idx_vec) {
    sub_out <- Momocs::Out(scaled_coos[idx_vec])
    as.matrix(Momocs::efourier(sub_out, nb.h = h_opt, norm = FALSE, start = FALSE)$coe)
  }
  coe_chunks <- if (exists("packera_mclapply", mode = "function")) {
    packera_mclapply(chunk_ids, efa_chunk_fn, mc.cores = n_w)
  } else {
    lapply(chunk_ids, efa_chunk_fn)
  }
  harm_mat <- do.call(rbind, coe_chunks)
  rownames(harm_mat) <- names(scaled_coos)
  ef_res <- tryCatch({
    tpl <- Momocs::efourier(Momocs::slice(out_obj, 1L), nb.h = h_opt, norm = FALSE, start = FALSE)
    tpl$coe <- harm_mat
    tpl$fac <- fac_df
    tpl
  }, error = function(e) list(coe = harm_mat, fac = fac_df))
  list(out_obj = out_obj, fac_df = fac_df, h_opt = as.integer(h_opt), ef_res = ef_res, harm_mat = harm_mat)
}

aggregate_step03_otu_medians_fast <- function(closed_leaf_df, all_harm_cols, sym_names) {
  n_rows <- nrow(closed_leaf_df)
  if (n_rows == 0L) return(data.frame())

  cat_col <- as.character(closed_leaf_df$catalogNumber)
  pid_col <- suppressWarnings(as.integer(closed_leaf_df$plant_individual_id))
  pid_col[is.na(pid_col)] <- 0L
  grp_key <- paste(cat_col, pid_col, sep = "\r")
  u_keys <- unique(grp_key)
  raw_groups <- split(seq_len(n_rows), factor(grp_key, levels = u_keys))

  has_out_col <- "is_contour_outlier" %in% names(closed_leaf_df)
  out_vec <- if (has_out_col) as.logical(closed_leaf_df$is_contour_outlier) else rep(FALSE, n_rows)
  out_vec[is.na(out_vec)] <- FALSE
  type_vec <- if ("is_type_specimen" %in% names(closed_leaf_df)) (!is.na(closed_leaf_df$is_type_specimen) & as.logical(closed_leaf_df$is_type_specimen)) else rep(FALSE, n_rows)
  out_vec[type_vec] <- FALSE
  lz_vec <- if ("leaf_zone" %in% names(closed_leaf_df)) tolower(trimws(as.character(closed_leaf_df$leaf_zone))) else rep("basal", n_rows)
  lz_vec[is.na(lz_vec) | !(lz_vec %in% c("basal", "cauline"))] <- "basal"

  grp_list <- vector("list", length(raw_groups))
  kept_g <- 0L
  first_idx <- integer(length(raw_groups))
  is_caul_fb_vec <- logical(length(raw_groups))
  for (g in seq_along(raw_groups)) {
    idx <- raw_groups[[g]]
    if (length(idx) == 0L) next
    if (has_out_col && any(!out_vec[idx], na.rm = TRUE)) {
      idx <- idx[!out_vec[idx]]
    }
    if (length(idx) == 0L) next
    bas_idx <- idx[lz_vec[idx] == "basal"]
    kept_g <- kept_g + 1L
    if (length(bas_idx) >= 1L) {
      grp_list[[kept_g]] <- bas_idx
      first_idx[kept_g] <- bas_idx[1L]
      is_caul_fb_vec[kept_g] <- FALSE
    } else {
      grp_list[[kept_g]] <- idx
      first_idx[kept_g] <- idx[1L]
      is_caul_fb_vec[kept_g] <- TRUE
    }
  }
  if (kept_g == 0L) return(data.frame())
  grp_list <- grp_list[seq_len(kept_g)]
  first_idx <- first_idx[seq_len(kept_g)]
  is_caul_fb_vec <- is_caul_fb_vec[seq_len(kept_g)]
  ord_g <- order(pid_col[first_idx], cat_col[first_idx])
  grp_list <- grp_list[ord_g]
  first_idx <- first_idx[ord_g]
  is_caul_fb_vec <- is_caul_fb_vec[ord_g]
  n_spec <- length(grp_list)

  harm_mat <- as.matrix(closed_leaf_df[, all_harm_cols, drop = FALSE])
  storage.mode(harm_mat) <- "double"
  sym_col_pos <- match(sym_names, all_harm_cols)
  sym_col_pos <- sym_col_pos[!is.na(sym_col_pos)]
  n_harm <- ncol(harm_mat)

  pull_num <- function(cn) {
    if (cn %in% names(closed_leaf_df)) suppressWarnings(as.numeric(closed_leaf_df[[cn]])) else rep(NA_real_, n_rows)
  }
  cs_v <- pull_num("centroid_size")
  sol_v <- pull_num("solidity")
  di_v <- pull_num("dissection_index")
  ar_v <- pull_num("aspect_ratio")
  circ_v <- pull_num("circularity")
  arch_v <- pull_num("archetype_score")
  area_v <- pull_num("area_px")
  has_amm2 <- "area_mm2" %in% names(closed_leaf_df)
  amm2_v <- if (has_amm2) pull_num("area_mm2") else NULL
  svar_v <- pull_num("sym_var")
  avar_v <- pull_num("asym_var")
  pasym_v <- pull_num("pct_asymmetry")
  hrough_v <- pull_num("harmonic_roughness")
  lqw_v <- pull_num("leaf_quality_weight")
  has_pet <- "petiole_length_mm" %in% names(closed_leaf_df)
  has_petw <- "petiole_width_mm" %in% names(closed_leaf_df)
  has_bld <- "blade_length_mm" %in% names(closed_leaf_df)
  has_bldw <- "blade_width_mm" %in% names(closed_leaf_df)
  pet_v <- if (has_pet) pull_num("petiole_length_mm") else NULL
  petw_v <- if (has_petw) pull_num("petiole_width_mm") else NULL
  bld_v <- if (has_bld) pull_num("blade_length_mm") else NULL
  bldw_v <- if (has_bldw) pull_num("blade_width_mm") else NULL
  pbr_raw <- if ("petiole_to_blade_ratio" %in% names(closed_leaf_df)) pull_num("petiole_to_blade_ratio") else rep(NA_real_, n_rows)
  if (has_pet && has_bld) {
    miss_pbr <- !is.finite(pbr_raw) & is.finite(pet_v) & is.finite(bld_v) & (bld_v > 1e-6)
    if (any(miss_pbr)) pbr_raw[miss_pbr] <- round(pet_v[miss_pbr] / bld_v[miss_pbr], 4)
  }
  has_pbr <- ("petiole_to_blade_ratio" %in% names(closed_leaf_df)) || any(is.finite(pbr_raw))
  pbr_v <- if (has_pbr) pbr_raw else NULL
  has_ect1 <- "ECT_PC1" %in% names(closed_leaf_df)
  has_ect2 <- "ECT_PC2" %in% names(closed_leaf_df)
  has_ect3 <- "ECT_PC3" %in% names(closed_leaf_df)
  has_ectlc <- "ect_lobe_complexity" %in% names(closed_leaf_df)
  has_blln <- "blade_cie_l_norm" %in% names(closed_leaf_df)
  has_bllsd <- "blade_cie_l_sd" %in% names(closed_leaf_df)
  blln_v <- if (has_blln) pull_num("blade_cie_l_norm") else NULL
  bllsd_v <- if (has_bllsd) pull_num("blade_cie_l_sd") else NULL
  ect1_v <- if (has_ect1) pull_num("ECT_PC1") else NULL
  ect2_v <- if (has_ect2) pull_num("ECT_PC2") else NULL
  ect3_v <- if (has_ect3) pull_num("ECT_PC3") else NULL
  ectlc_v <- if (has_ectlc) pull_num("ect_lobe_complexity") else NULL
  has_dis <- "is_dissected" %in% names(closed_leaf_df)
  dis_v <- if (has_dis) as.logical(closed_leaf_df$is_dissected) else NULL

  med_fin <- function(v, digits) {
    vf <- v[is.finite(v)]
    if (length(vf) == 0L) NA_real_ else round(stats::median(vf), digits)
  }
  q75_fin <- function(v) {
    vf <- v[is.finite(v)]
    if (length(vf) == 0L) NA_real_ else round(as.numeric(stats::quantile(vf, probs = 0.75, na.rm = TRUE, names = FALSE)), 4)
  }

  agg_chunk_fn <- function(g_positions) {
    m <- length(g_positions)
    h_out <- matrix(NA_real_, nrow = m, ncol = n_harm)
    n_lv <- integer(m)
    f_var <- numeric(m)
    h_idx <- numeric(m)
    cs_m <- rep(NA_real_, m)
    sol_m <- rep(NA_real_, m)
    q25_m <- rep(NA_real_, m)
    di_m <- rep(NA_real_, m)
    dfrac_m <- rep(NA_real_, m)
    ar_m <- rep(NA_real_, m)
    circ_m <- rep(NA_real_, m)
    arch_m <- rep(NA_real_, m)
    area_m <- rep(NA_real_, m)
    amm2_m <- rep(NA_real_, m)
    sv_m <- rep(NA_real_, m)
    av_m <- rep(NA_real_, m)
    pa_m <- rep(NA_real_, m)
    hr_m <- rep(NA_real_, m)
    lq_m <- rep(NA_real_, m)
    pet_m <- rep(NA_real_, m)
    petw_m <- rep(NA_real_, m)
    bld_m <- rep(NA_real_, m)
    bldw_m <- rep(NA_real_, m)
    pbr_m <- rep(NA_real_, m)
    ect1_m <- rep(NA_real_, m)
    ect2_m <- rep(NA_real_, m)
    ect3_m <- rep(NA_real_, m)
    ectlc_m <- rep(NA_real_, m)
    blln_m <- rep(NA_real_, m)
    bllsd_m <- rep(NA_real_, m)
    dis_m <- rep(NA, m)
    out_m <- logical(m)

    for (k in seq_len(m)) {
      idx <- grp_list[[g_positions[k]]]
      nl <- length(idx)
      n_lv[k] <- nl
      if (nl == 1L) {
        hm <- harm_mat[idx, ]
        h_out[k, ] <- hm
        f_var[k] <- 0.0
        h_idx[k] <- 0.0
      } else {
        sub_h <- harm_mat[idx, , drop = FALSE]
        hm <- apply(sub_h, 2L, stats::median, na.rm = TRUE)
        h_out[k, ] <- hm
        sub_sym <- sub_h[, sym_col_pos, drop = FALSE]
        hm_sym <- hm[sym_col_pos]
        f_var[k] <- round(mean(apply(sub_sym, 2L, stats::var, na.rm = TRUE), na.rm = TRUE), 8)
        h_idx[k] <- round(mean(sqrt(rowSums(sweep(sub_sym, 2L, hm_sym, "-")^2, na.rm = TRUE)), na.rm = TRUE), 6)
      }
      mc <- stats::median(cs_v[idx], na.rm = TRUE)
      cs_m[k] <- if (is.finite(mc) && mc > 0) round(mc, 4) else NA_real_

      sv <- sol_v[idx]
      sv_f <- sv[is.finite(sv)]
      sol_m[k] <- if (length(sv_f) == 0L) NA_real_ else round(stats::median(sv_f), 4)
      q25_m[k] <- if (length(sv_f) == 0L) NA_real_ else round(as.numeric(stats::quantile(sv_f, probs = 0.25, na.rm = TRUE, names = FALSE)), 4)
      dr <- if (has_dis) dis_v[idx] else rep(NA, nl)
      dl <- ifelse(!is.na(dr), dr, FALSE) | (is.finite(sv) & sv < 0.72)
      dfrac_m[k] <- if (any(!is.na(dr)) || length(sv_f) > 0L) round(mean(dl, na.rm = TRUE), 4) else NA_real_

      di_m[k] <- med_fin(di_v[idx], 4L)
      ar_m[k] <- med_fin(ar_v[idx], 4L)
      circ_m[k] <- med_fin(circ_v[idx], 4L)
      arch_m[k] <- med_fin(arch_v[idx], 4L)
      area_m[k] <- med_fin(area_v[idx], 1L)
      if (has_amm2) amm2_m[k] <- med_fin(amm2_v[idx], 2L)
      sv_m[k] <- med_fin(svar_v[idx], 6L)
      av_m[k] <- med_fin(avar_v[idx], 6L)
      pa_m[k] <- med_fin(pasym_v[idx], 6L)
      hr_m[k] <- med_fin(hrough_v[idx], 6L)
      lq_m[k] <- med_fin(lqw_v[idx], 6L)
      if (has_pet) pet_m[k] <- q75_fin(pet_v[idx])
      if (has_petw) petw_m[k] <- q75_fin(petw_v[idx])
      if (has_bld) bld_m[k] <- q75_fin(bld_v[idx])
      if (has_bldw) bldw_m[k] <- q75_fin(bldw_v[idx])
      if (has_pbr) pbr_m[k] <- q75_fin(pbr_v[idx])
      if (has_ect1) ect1_m[k] <- med_fin(ect1_v[idx], 6L)
      if (has_ect2) ect2_m[k] <- med_fin(ect2_v[idx], 6L)
      if (has_ect3) ect3_m[k] <- med_fin(ect3_v[idx], 6L)
      if (has_ectlc) ectlc_m[k] <- med_fin(ectlc_v[idx], 4L)
      if (has_blln) blln_m[k] <- med_fin(blln_v[idx], 4L)
      if (has_bllsd) bllsd_m[k] <- med_fin(bllsd_v[idx], 4L)
      if (has_dis) dis_m[k] <- if (any(!is.na(dr))) any(dr, na.rm = TRUE) else NA
      if (has_out_col) out_m[k] <- all(out_vec[idx], na.rm = TRUE)
    }

    list(
      h = h_out, nl = n_lv, fv = f_var, hi = h_idx, cs = cs_m,
      sol = sol_m, q25 = q25_m, di = di_m, df = dfrac_m,
      ar = ar_m, circ = circ_m, arch = arch_m, area = area_m, amm2 = amm2_m,
      sv = sv_m, av = av_m, pa = pa_m, hr = hr_m, lq = lq_m,
      pet = pet_m, petw = petw_m, bld = bld_m, bldw = bldw_m, pbr = pbr_m,
      ect1 = ect1_m, ect2 = ect2_m, ect3 = ect3_m, ectlc = ectlc_m,
      blln = blln_m, bllsd = bllsd_m,
      dis = dis_m, out = out_m
    )
  }

  n_w <- if (n_spec >= 300L && exists("resolve_packera_workers", mode = "function")) {
    resolve_packera_workers(32L, ceiling(n_spec / 100L))
  } else {
    1L
  }
  if (n_w <= 1L) {
    res <- agg_chunk_fn(seq_len(n_spec))
  } else {
    g_chunks <- split(seq_len(n_spec), cut(seq_len(n_spec), breaks = n_w, labels = FALSE))
    c_list <- packera_mclapply(g_chunks, agg_chunk_fn, mc.cores = n_w)
    pull_c <- function(nm) unlist(lapply(c_list, `[[`, nm), use.names = FALSE)
    res <- list(
      h = do.call(rbind, lapply(c_list, `[[`, "h")),
      nl = pull_c("nl"), fv = pull_c("fv"), hi = pull_c("hi"), cs = pull_c("cs"),
      sol = pull_c("sol"), q25 = pull_c("q25"), di = pull_c("di"), df = pull_c("df"),
      ar = pull_c("ar"), circ = pull_c("circ"), arch = pull_c("arch"), area = pull_c("area"), amm2 = pull_c("amm2"),
      sv = pull_c("sv"), av = pull_c("av"), pa = pull_c("pa"), hr = pull_c("hr"), lq = pull_c("lq"),
      pet = pull_c("pet"), petw = pull_c("petw"), bld = pull_c("bld"), bldw = pull_c("bldw"), pbr = pull_c("pbr"),
      ect1 = pull_c("ect1"), ect2 = pull_c("ect2"), ect3 = pull_c("ect3"), ectlc = pull_c("ectlc"),
      blln = pull_c("blln"), bllsd = pull_c("bllsd"),
      dis = pull_c("dis"), out = pull_c("out")
    )
  }
  colnames(res$h) <- all_harm_cols

  det_vec <- if ("determiner_tier" %in% names(closed_leaf_df)) as.character(closed_leaf_df$determiner_tier[first_idx]) else rep("Tier_3_Bronze", n_spec)
  det_vec[is.na(det_vec) | !nzchar(det_vec)] <- "Tier_3_Bronze"
  auth_init <- if ("authority_weight" %in% names(closed_leaf_df)) suppressWarnings(as.numeric(closed_leaf_df$authority_weight[first_idx])) else rep(0.0, n_spec)
  auth_init[!is.finite(auth_init)] <- 0.0
  log_cs <- ifelse(!is.na(res$cs) & res$cs > 0, round(log10(res$cs), 6), NA_real_)

  spec_df <- data.frame(
    catalogNumber = cat_col[first_idx],
    plant_individual_id = pid_col[first_idx],
    leaf_count = res$nl,
    morph_tier_source = ifelse(is_caul_fb_vec, "Cauline_Only_Fallback", "Basal_Rosette"),
    foliar_variance = res$fv,
    within_rosette_efa_dispersion = res$hi,
    scientificName = closed_leaf_df$scientificName[first_idx],
    species_raw = closed_leaf_df$species_raw[first_idx],
    determiner_tier = det_vec,
    authority_weight = auth_init,
    centroid_size = res$cs,
    log_centsize = log_cs,
    stringsAsFactors = FALSE
  )

  other_meta <- intersect(
    c("pixels_per_mm", "is_calibrated_mm", "scale_source", "mask_is_upright", "species_standardized", "taxon", "institutionCode", "identifiedBy", "authority_tier", "type_status", "is_type_specimen", "is_promoted_bronze", "anchor_tier_source",
      "stateProvince", "county", "locality", "decimalLatitude", "decimalLongitude", "coordinateUncertainty",
      "year", "month", "day", "eventDate", "regional_group", "spatial_balance_weight"),
    names(closed_leaf_df)
  )
  for (m_col in other_meta) {
    col_all <- closed_leaf_df[[m_col]]
    v_first <- col_all[first_idx]
    na_g <- which(is.na(v_first))
    if (length(na_g) > 0L) {
      for (gi in na_g) {
        vals <- col_all[grp_list[[gi]]]
        vals <- vals[!is.na(vals)]
        if (length(vals) > 0L) v_first[gi] <- vals[1L]
      }
    }
    spec_df[[m_col]] <- v_first
  }

  spec_df$basal_leaf_count <- ifelse(is_caul_fb_vec, 0L, as.integer(res$nl))
  spec_df$solidity <- res$sol
  spec_df$q25_solidity <- res$q25
  spec_df$dissection_index <- res$di
  spec_df$dissected_leaf_fraction <- res$df
  spec_df$aspect_ratio <- res$ar
  spec_df$circularity <- res$circ
  spec_df$archetype_score <- res$arch
  spec_df$area_px <- res$area
  if (has_amm2) spec_df$area_mm2 <- res$amm2
  spec_df$sym_var <- res$sv
  spec_df$asym_var <- res$av
  spec_df$pct_asymmetry <- res$pa
  spec_df$harmonic_roughness <- res$hr
  spec_df$leaf_quality_weight <- res$lq
  if (has_pet) spec_df$petiole_length_mm <- res$pet
  if (has_petw) spec_df$petiole_width_mm <- res$petw
  if (has_bld) spec_df$blade_length_mm <- res$bld
  if (has_bldw) spec_df$blade_width_mm <- res$bldw
  if (has_pbr) spec_df$petiole_to_blade_ratio <- res$pbr
  if (has_ect1) spec_df$ECT_PC1 <- res$ect1
  if (has_ect2) spec_df$ECT_PC2 <- res$ect2
  if (has_ect3) spec_df$ECT_PC3 <- res$ect3
  if (has_ectlc) spec_df$ect_lobe_complexity <- res$ectlc
  if (has_blln) spec_df$blade_cie_l_norm <- res$blln
  if (has_bllsd) spec_df$blade_cie_l_sd <- res$bllsd
  if (has_dis) spec_df$is_dissected <- res$dis
  if (has_out_col) spec_df$is_contour_outlier <- res$out

  out_spec_df <- cbind(spec_df, as.data.frame(res$h))
  attr(out_spec_df, "step03_fast_aggregated") <- TRUE
  out_spec_df
}
