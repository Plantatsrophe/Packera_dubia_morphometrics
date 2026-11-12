#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_figures_utils.R
# Project: Packera Genus-Wide Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Shared publication figure utilities across R pipeline phases (Steps 03, 04,
#   06, and 07). Enforces:
#     - Dual vector .pdf and 300-DPI .png atomic exports (write-then-replace)
#     - Colorblind-safe Okabe-Ito palettes (#0072B2, #D55E00, #009E73, #CC79A7)
#     - Dynamically italicized Packera taxon labels
# ==============================================================================

EXTENDED_COLORBLIND_HEX <- c(
  "#E69F00", "#56B4E9", "#F0E442", "#333333", "#882255",
  "#44AA99", "#999933", "#332288", "#117733", "#AA4499",
  "#DDCC77", "#661100", "#6699CC", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#008080"
)

HOLLOW_PCH_CYCLE <- c(1L, 2L, 5L, 0L, 6L, 4L)
if (!exists("fast_euclidean_permanova_trace", mode = "function")) {
  for (.dl_p in c("scripts/morphometrics/dissected_leaf_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "dissected_leaf_utils.R"))) {
    if (!is.na(.dl_p) && nzchar(.dl_p) && file.exists(.dl_p)) { source(.dl_p); break }
  }
}

apply_bh_fdr_correction <- function(df_or_p, p_col = NULL) {
  if (is.null(df_or_p)) return(NULL)
  if (is.numeric(df_or_p)) {
    dims <- dim(df_or_p)
    dnames <- dimnames(df_or_p)
    adj <- stats::p.adjust(as.numeric(df_or_p), method = "BH")
    if (!is.null(dims)) {
      dim(adj) <- dims
      dimnames(adj) <- dnames
    }
    return(adj)
  }
  if (is.data.frame(df_or_p)) {
    cand_cols <- if (!is.null(p_col)) p_col else intersect(
      c("p_value", "P_value", "P_Value", "Equivalency_p_value", "p_value_D", "p_value_I", "Pr(>F)"),
      colnames(df_or_p)
    )
    for (col in cand_cols) {
      if (col %in% colnames(df_or_p)) {
        raw_p <- suppressWarnings(as.numeric(df_or_p[[col]]))
        adj_p <- stats::p.adjust(raw_p, method = "BH")
        df_or_p[[paste0(col, "_bh")]] <- signif(adj_p, 4)
        if (!("p_adj_bh" %in% colnames(df_or_p))) {
          df_or_p$p_adj_bh <- signif(adj_p, 4)
        }
      }
    }
    return(df_or_p)
  }
  df_or_p
}

TAXONOMIC_TREATMENT_PALETTE <- c(
  "Confirmed Species"                   = "#009E73",
  "Species_Confirmed"                   = "#009E73",
  "Variety / Ecotype"                   = "#56B4E9",
  "Subspecies_Ecotype"                  = "#56B4E9",
  "Ecophenotypic_Plasticity"            = "#56B4E9",
  "Hybrid / Introgressant"              = "#CC79A7",
  "Hybrid_Intergrade_Swarm"             = "#CC79A7",
  "Putative_Hybrid_or_Ecophenotype"     = "#CC79A7",
  "Misidentification (Reassign)"        = "#E69F00",
  "Misidentification_Reassign"          = "#E69F00",
  "Putative Cryptic Species"            = "#D55E00",
  "Putative_Cryptic_Species"            = "#D55E00",
  "Misidentification / Cryptic Species" = "#D55E00"
)

compute_spatial_balance_weights <- function(df, bandwidth_km = 25.0, regional_bandwidth_km = 75.0, active_mask = NULL, recompute = FALSE) {
  if (is.null(df) || nrow(df) == 0L) return(numeric(0))
  n_rows <- nrow(df)
  if (!isTRUE(recompute) && is.null(active_mask) && ("spatial_balance_weight" %in% names(df))) {
    w_in <- suppressWarnings(as.numeric(df$spatial_balance_weight))
    if (any(is.finite(w_in) & w_in > 0)) {
      w_in[!is.finite(w_in) | w_in <= 0] <- 1.0
      return(pmin(3.0, pmax(0.20, w_in)))
    }
  }
  sp_col <- if ("species_standardized" %in% names(df)) "species_standardized" else if ("scientificName" %in% names(df)) "scientificName" else if ("species_raw" %in% names(df)) "species_raw" else NULL
  taxa <- if (!is.null(sp_col)) as.character(df[[sp_col]]) else rep("Unknown", n_rows)
  taxa[is.na(taxa) | !nzchar(trimws(taxa))] <- "Unknown"

  if (!("decimalLatitude" %in% names(df)) || !("decimalLongitude" %in% names(df))) return(rep(1.0, n_rows))
  lats <- suppressWarnings(as.numeric(df$decimalLatitude))
  lons <- suppressWarnings(as.numeric(df$decimalLongitude))
  ok_geo <- is.finite(lats) & is.finite(lons) & (abs(lats) <= 90.0) & (abs(lons) <= 180.0)

  eval_cohort_kernel <- function(cohort_idx) {
    w_c <- rep(1.0, length(cohort_idx))
    if (length(cohort_idx) == 0L) return(w_c)
    tx_c <- taxa[cohort_idx]
    for (sp in unique(tx_c)) {
      pos_sp <- which(tx_c == sp)
      rows_sp <- cohort_idx[pos_sp]
      geo_pos <- pos_sp[ok_geo[rows_sp]]
      if (length(geo_pos) < 2L) next
      geo_rows <- cohort_idx[geo_pos]
      rlat <- (lats[geo_rows] * (pi / 180.0)) + 0.0
      rlon <- (lons[geo_rows] * (pi / 180.0)) + 0.0
      coord_key <- paste(rlat, rlon, sep = "_")
      u_keys <- unique(coord_key)
      inv_idx <- match(coord_key, u_keys)
      first_idx <- match(u_keys, coord_key)
      counts_u <- as.numeric(tabulate(inv_idx, nbins = length(u_keys)))
      u_rlat <- rlat[first_idx]
      u_rlon <- rlon[first_idx]
      n_u <- length(u_keys)
      dlat <- outer(u_rlat, u_rlat, "-")
      dlon <- outer(u_rlon, u_rlon, "-")
      a_mat <- sin(dlat / 2.0)^2 + (cos(u_rlat) %o% cos(u_rlat)) * (sin(dlon / 2.0)^2)
      dist_km <- 2.0 * 6371.0088 * asin(pmin(1.0, pmax(0.0, sqrt(pmax(0.0, a_mat)))))
      kern_fine <- exp(-(dist_km^2) / (2.0 * (bandwidth_km^2)))
      kern_reg <- exp(-(dist_km^2) / (2.0 * (regional_bandwidth_km^2)))
      kern_u <- matrix(kern_fine + 0.35 * kern_reg, nrow = n_u, ncol = n_u)
      diag(kern_u) <- 1.35
      # Exact identity to 1.0 + sum_{j != i} (kern_fine[i,j] + 0.35*kern_reg[i,j]):
      d_local <- (as.numeric(kern_u %*% counts_u) - 0.35)[inv_idx]
      w_raw <- 1.0 / sqrt(pmax(d_local, 1e-8))
      w_sp <- pmin(3.0, pmax(0.20, w_raw))
      for (iter in 1:3) {
        m_cur <- mean(w_sp, na.rm = TRUE)
        if (is.finite(m_cur) && m_cur > 0) w_sp <- pmin(3.0, pmax(0.20, w_sp / m_cur))
      }
      w_c[geo_pos] <- round(w_sp, 6)
    }
    w_c
  }

  w_out <- rep(1.0, n_rows)
  use_split <- !is.null(active_mask) && (length(active_mask) == n_rows) && any(as.logical(active_mask), na.rm = TRUE)
  if (!use_split && ("authority_weight" %in% names(df))) {
    aw_v <- suppressWarnings(as.numeric(df$authority_weight))
    if (any(is.finite(aw_v) & aw_v > 0.0) && any(!is.finite(aw_v) | aw_v <= 0.0)) {
      active_mask <- is.finite(aw_v) & (aw_v > 0.0)
      use_split <- TRUE
    }
  }
  if (use_split) {
    act_m <- !is.na(active_mask) & as.logical(active_mask)
    act_idx <- which(act_m)
    pas_idx <- which(!act_m)
    w_out[act_idx] <- eval_cohort_kernel(act_idx)
    if (length(pas_idx) > 0L) w_out[pas_idx] <- eval_cohort_kernel(pas_idx)
  } else {
    w_out <- eval_cohort_kernel(seq_len(n_rows))
  }
  w_out
}

format_predictor_labels <- function(var_vec) {
  dict <- c(
    "fine_bedrock_depth_cm"   = "Depth to Lithic Bedrock (cm)",
    "ssurgo_bedrock_depth_cm" = "Depth to Lithic Bedrock (cm)",
    "ssurgo_rock_outcrop_pct" = "SSURGO 50m Rock Outcrop (%)",
    "fine_soil_sand"          = "Soil Sand Fraction (%)",
    "soil_sand"               = "Soil Sand Fraction (%)",
    "fine_soil_ph"            = "Soil pH (H2O)",
    "soil_ph"                 = "Soil pH (H2O)",
    "fine_soil_clay"          = "Soil Clay Fraction (%)",
    "soil_clay"               = "Soil Clay Fraction (%)",
    "ssurgo_aws_0_50cm"       = "Available Water Storage 0-50cm (cm)",
    "topo_tpi_30m"            = "30m Topographic Position Index (TPI, m)",
    "topo_twi_30m"            = "30m Topographic Wetness Index (TWI)",
    "topo_hli_30m"            = "30m Heat Load Index (HLI)",
    "ilr_texture_1"           = "Soil Texture ILR-1 (Sand vs Silt-Clay)",
    "ilr_texture_2"           = "Soil Texture ILR-2 (Silt vs Clay)",
    "fine_cfvo_pct"           = "Coarse Fragment Volume (%, 0-15cm)",
    "topo_slope_deg"          = "30m Slope Gradient (deg)",
    "soil_cec"                = "Cation Exchange Capacity (cmol/kg)",
    "ssurgo_cec"              = "Cation Exchange Capacity (cmol/kg)",
    "soil_bulk_density"       = "Soil Bulk Density (g/cm3)",
    "bio1_temp_mean"          = "Mean Annual Temp. (BIO1, C)",
    "bio4_temp_seasonality"   = "Temp. Seasonality (BIO4)",
    "bio12_precip_annual"     = "Annual Precipitation (BIO12, mm)",
    "bio15_precip_seasonality"= "Precip. Seasonality (BIO15)",
    "delta_doy"               = "Phenological Anomaly (Delta DOY, d)",
    "doy"                     = "Flowering Day of Year (DOY)"
  )
  vapply(as.character(var_vec), function(v) {
    if (is.na(v) || !nzchar(v)) return("Unknown Predictor")
    if (v %in% names(dict)) return(unname(dict[v]))
    if (grepl("^MEM[0-9]+$", v)) return(sprintf("Spatial Eigenfunction (%s)", v))
    v
  }, character(1), USE.NAMES = FALSE)
}

adaptive_point_alpha <- function(n_pts) {
  if (is.na(n_pts) || n_pts <= 1000L) 0.75 else 0.45
}

.pub_figures_script_dir <- tryCatch({
  ofiles <- unlist(lapply(sys.frames(), function(fr) fr$ofile))
  cand <- if (length(ofiles) > 0L) tail(ofiles, 1L) else sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1L]
  if (!is.na(cand) && nzchar(cand)) dirname(normalizePath(cand, mustWork = FALSE)) else ""
}, error = function(e) "")

.packera_taxon_syn_env <- new.env(parent = emptyenv())
.packera_taxon_syn_env$db_sig <- ""
.packera_taxon_syn_env$syn_map <- character(0)
.packera_taxon_syn_env$py_memo <- character(0)

.expand_synonym_binomial_keys <- function(syn_vec) {
  if (length(syn_vec) == 0L) return(character(0))
  nms <- tolower(gsub("\\s+", " ", trimws(names(syn_vec))))
  vals <- gsub("^Packera\\s+(?:x|×)\\s+", "Packera ", trimws(as.character(syn_vec)), ignore.case = TRUE)
  keep <- !is.na(nms) & nzchar(nms) & !is.na(vals) & nzchar(vals)
  if (!any(keep)) return(character(0))
  out <- stats::setNames(vals[keep], nms[keep])
  bin_pat <- "^([a-z]+\\s+(?:(?:x|×)\\s*)?[a-z][a-z-]*)\\b.*$"
  is_bin <- grepl(bin_pat, names(out))
  if (any(is_bin)) {
    b_keys <- sub(bin_pat, "\\1", names(out)[is_bin])
    b_vals <- unname(out[is_bin])
    new_m <- !(b_keys %in% names(out)) & !duplicated(b_keys)
    if (any(new_m)) out <- c(out, stats::setNames(b_vals[new_m], b_keys[new_m]))
  }
  out
}

.get_packera_synonym_map <- function() {
  cand_dbs <- character(0)
  env_db <- trimws(Sys.getenv("PACKERA_DB_PATH", unset = ""))
  if (nzchar(env_db) && env_db != ":memory:") cand_dbs <- c(cand_dbs, env_db)
  if (exists("resolve_packera_db_path", mode = "function")) {
    rdb <- tryCatch(as.character(resolve_packera_db_path()), error = function(e) "")
    if (!is.na(rdb) && nzchar(rdb) && rdb != ":memory:") cand_dbs <- c(cand_dbs, rdb)
  }
  cand_dbs <- c(cand_dbs, "data/packera_pipeline.sqlite")
  if (nzchar(.pub_figures_script_dir)) {
    cand_dbs <- c(cand_dbs, file.path(.pub_figures_script_dir, "..", "..", "data", "packera_pipeline.sqlite"))
  }
  cand_dbs <- unique(normalizePath(cand_dbs[file.exists(cand_dbs)], mustWork = FALSE))
  sig <- if (length(cand_dbs) == 0L) "none" else paste(cand_dbs, as.numeric(file.info(cand_dbs)$mtime), sep = "@", collapse = "|")
  if (identical(sig, .packera_taxon_syn_env$db_sig)) return(.packera_taxon_syn_env$syn_map)

  merged <- character(0)
  if (length(cand_dbs) > 0L && requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {
    for (dbp in cand_dbs) {
      res <- tryCatch({
        con <- DBI::dbConnect(RSQLite::SQLite(), dbname = dbp, flags = RSQLite::SQLITE_RO)
        on.exit(tryCatch(DBI::dbDisconnect(con), error = function(e) NULL), add = TRUE)
        if ("taxon_synonyms" %in% DBI::dbListTables(con)) {
          DBI::dbGetQuery(con, "SELECT raw_name_lower, canonical_packera FROM taxon_synonyms;")
        } else NULL
      }, error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0L) {
        cur <- .expand_synonym_binomial_keys(stats::setNames(as.character(res$canonical_packera), as.character(res$raw_name_lower)))
        add_nms <- setdiff(names(cur), names(merged))
        if (length(add_nms) > 0L) merged <- c(merged, cur[add_nms])
      }
    }
  }
  .packera_taxon_syn_env$db_sig <- sig
  .packera_taxon_syn_env$syn_map <- merged
  merged
}

.resolve_taxa_via_harvester_py <- function(query_vec) {
  if (length(query_vec) == 0L) return(character(0))
  u_q <- unique(as.character(query_vec))
  memo <- .packera_taxon_syn_env$py_memo
  miss <- setdiff(u_q, names(memo))
  if (length(miss) > 0L) {
    ws_root <- if (nzchar(.pub_figures_script_dir)) normalizePath(file.path(.pub_figures_script_dir, "..", ".."), mustWork = FALSE) else normalizePath(".", mustWork = FALSE)
    py_cands <- c(file.path(ws_root, ".venv", "bin", "python"), file.path(".venv", "bin", "python"), Sys.which("python3"))
    py_bin <- py_cands[nzchar(py_cands) & file.exists(py_cands)][1L]
    resolved_miss <- stats::setNames(rep(NA_character_, length(miss)), miss)
    if (!is.na(py_bin) && nzchar(py_bin) && file.exists(file.path(ws_root, "scripts", "core", "harvester.py"))) {
      tf_in <- tempfile("packera_tax_in_"); on.exit(unlink(tf_in), add = TRUE)
      writeLines(miss, con = tf_in, useBytes = TRUE)
      py_code <- sprintf(
        "import sys; sys.path.insert(0, %s); from scripts.core.harvester import standardize_packera_taxon as f\nwith open(%s, 'r', encoding='utf-8') as fh:\n for ln in fh:\n  r = f(ln.rstrip('\\r\\n'))\n  print(r if r else '')",
        shQuote(ws_root), shQuote(tf_in)
      )
      py_out <- tryCatch(suppressWarnings(system2(py_bin, args = c("-c", shQuote(py_code)), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
      if (length(py_out) == length(miss)) {
        clean_out <- trimws(as.character(py_out))
        resolved_miss <- stats::setNames(ifelse(nzchar(clean_out), clean_out, NA_character_), miss)
      }
    }
    memo <- c(memo, resolved_miss)
    .packera_taxon_syn_env$py_memo <- memo
  }
  unname(memo[as.character(query_vec)])
}

standardize_taxon_name <- function(raw_vec) {
  if (is.null(raw_vec) || length(raw_vec) == 0L) return(character(0))
  chr_vec <- as.character(raw_vec)
  u_vals <- unique(chr_vec)
  syn_map <- .get_packera_synonym_map()
  has_syn <- length(syn_map) > 0L
  need_py_idx <- integer(0)
  fb_packera <- rep(NA_character_, length(u_vals))
  open_ep <- c("sp", "spp", "cf", "aff", "indet", "unknown", "na")
  syn_gen_pat <- "^(packera|p\\.|senecio|s\\.|cineraria|jacobaea|aureolaria)\\s*(?:(?:cf\\.?|aff\\.?)\\s+)?(x\\s+|×\\s*)?([a-z][a-z-]*).*$"

  std_u <- character(length(u_vals))
  for (i in seq_along(u_vals)) {
    val <- u_vals[i]
    if (is.na(val)) { std_u[i] <- NA_character_; next }
    s_orig <- trimws(val)
    s <- tolower(gsub("\\s+", " ", s_orig))
    if (!nzchar(s_orig) || s %in% c("none", "nan", "null", "na", "unknown", "unassigned", "unassigned_packera", "packera")) {
      std_u[i] <- NA_character_; next
    }
    if (grepl("^Morphotype_K[0-9]+$|^Ecotype_|\\[K[0-9]+\\]$", s_orig, ignore.case = TRUE)) {
      std_u[i] <- s_orig; next
    }
    if (has_syn && (s %in% names(syn_map))) {
      std_u[i] <- unname(syn_map[[s]]); next
    }
    if (has_syn && grepl("^[a-z]+\\s+[a-z-]+\\s+(?:var\\.?|subsp\\.?|ssp\\.?|f\\.?)\\s+[a-z-]+", s)) {
      v_raw <- sub("^([a-z]+\\s+[a-z-]+\\s+(?:var\\.?|subsp\\.?|ssp\\.?|f\\.?)\\s+[a-z-]+).*$", "\\1", s)
      v_dot <- gsub("\\b(var|subsp|ssp|f)\\.?\\s+", "\\1. ", v_raw)
      if (v_raw %in% names(syn_map)) { std_u[i] <- unname(syn_map[[v_raw]]); next }
      if (v_dot %in% names(syn_map)) { std_u[i] <- unname(syn_map[[v_dot]]); next }
    }
    if (grepl(syn_gen_pat, s)) {
      gen_tok <- sub(syn_gen_pat, "\\1", s)
      hyb_tok <- trimws(sub(syn_gen_pat, "\\2", s))
      ep <- sub("^[×-]+|[-]+$", "", sub(syn_gen_pat, "\\3", s))
      if (!nzchar(ep) || nchar(ep) < 2L || (ep %in% open_ep)) { std_u[i] <- NA_character_; next }
      gen_full <- if (gen_tok == "p.") "packera" else if (gen_tok == "s.") "senecio" else gen_tok
      is_hyb <- nzchar(hyb_tok)
      b_plain <- paste(gen_full, ep)
      b_hyb <- if (is_hyb) paste(gen_full, "×", ep) else b_plain
      if (has_syn && (b_hyb %in% names(syn_map))) { std_u[i] <- unname(syn_map[[b_hyb]]); next }
      if (has_syn && (b_plain %in% names(syn_map))) { std_u[i] <- unname(syn_map[[b_plain]]); next }
      if (gen_full == "packera") {
        p_can <- paste0("Packera ", ep)
        if (has_syn) { std_u[i] <- p_can; next }
        fb_packera[i] <- p_can
      }
      need_py_idx <- c(need_py_idx, i)
      std_u[i] <- NA_character_
      next
    }
    std_u[i] <- NA_character_
  }

  if (length(need_py_idx) > 0L) {
    py_res <- .resolve_taxa_via_harvester_py(u_vals[need_py_idx])
    std_u[need_py_idx] <- ifelse(!is.na(py_res) & nzchar(py_res), py_res, fb_packera[need_py_idx])
  }
  std_u[match(chr_vec, u_vals)]
}
standardize_packera_taxon <- standardize_taxon_name

format_active_complex_title <- function(taxa_vec, suffix) {
  clean_vec <- as.character(taxa_vec)
  clean_vec <- trimws(clean_vec[!is.na(clean_vec) & nzchar(trimws(clean_vec))])
  clean_vec <- clean_vec[!(clean_vec %in% c("Unknown", "Unassigned_Packera", "NA"))]
  unique_taxa <- unique(clean_vec)
  env_focal <- trimws(Sys.getenv("PACKERA_FOCAL_TAXON", unset = ""))
  if (nzchar(env_focal)) {
    env_std <- standardize_taxon_name(env_focal)[1L]
    if (is.na(env_std) || !nzchar(env_std)) env_std <- env_focal
    return(sprintf("%s Complex: %s", env_std, suffix))
  }
  if (length(unique_taxa) == 0L) {
    return(sprintf("Packera Taxon Complex: %s", suffix))
  }
  if (length(unique_taxa) == 1L) {
    return(sprintf("%s: %s", unique_taxa[1L], suffix))
  }
  sprintf("Packera Taxon Complex (%d Taxa): %s", length(unique_taxa), suffix)
}

.order_unique_packera_taxa <- function(taxa_vec, preserve_order = FALSE) {
  clean_v <- trimws(as.character(taxa_vec[!is.na(taxa_vec) & nzchar(trimws(as.character(taxa_vec)))]))
  if (length(clean_v) == 0L) return(character(0))
  u_raw <- unique(clean_v)
  u_non_unk <- setdiff(u_raw, "Unknown")
  if (length(u_non_unk) == 0L) return(character(0))
  if (isTRUE(preserve_order)) {
    return(u_non_unk)
  }
  sort(u_non_unk)
}

get_dynamic_packera_palette <- function(taxa_vec = NULL, preserve_order = FALSE) {
  if (is.null(taxa_vec) || length(taxa_vec) == 0L) return(c("Unknown" = "#7F7F7F"))
  u_non_unk <- .order_unique_packera_taxa(taxa_vec, preserve_order = preserve_order)
  if (length(u_non_unk) == 0L) return(c("Unknown" = "#7F7F7F"))
  pal_pool <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", EXTENDED_COLORBLIND_HEX)
  if (length(u_non_unk) > length(pal_pool)) {
    more_n <- length(u_non_unk) - length(pal_pool)
    hcl_ext <- toupper(grDevices::hcl.colors(more_n + 8L, palette = "Dark 3"))
    hcl_ext <- setdiff(hcl_ext, c(toupper(pal_pool), "#7F7F7F"))
    pal_pool <- c(pal_pool, hcl_ext)
  }
  cols <- rep_len(pal_pool, length(u_non_unk))
  names(cols) <- u_non_unk
  c(cols, "Unknown" = "#7F7F7F")
}

get_dynamic_packera_shapes <- function(taxa_vec = NULL, preserve_order = FALSE) {
  if (is.null(taxa_vec) || length(taxa_vec) == 0L) return(c("Unknown" = 0L))
  u_non_unk <- .order_unique_packera_taxa(taxa_vec, preserve_order = preserve_order)
  if (length(u_non_unk) == 0L) return(c("Unknown" = 0L))
  shps <- rep_len(c(1L, 2L, 5L, 0L, 6L, 4L), length(u_non_unk))
  names(shps) <- u_non_unk
  c(shps, "Unknown" = 0L)
}

get_dynamic_packera_ggplot_shapes <- function(taxa_vec = NULL, preserve_order = FALSE) {
  if (is.null(taxa_vec) || length(taxa_vec) == 0L) return(c("Unknown" = 25))
  u_non_unk <- .order_unique_packera_taxa(taxa_vec, preserve_order = preserve_order)
  if (length(u_non_unk) == 0L) return(c("Unknown" = 0L))
  shps <- rep_len(c(21, 24, 22, 23, 25), length(u_non_unk))
  names(shps) <- u_non_unk
  c(shps, "Unknown" = 25)
}

resolve_taxon_display_tier <- function(
  df,
  group_vec = NULL,
  focal_taxon = Sys.getenv("PACKERA_FOCAL_TAXON", "Packera dubia"),
  max_fg = 6L,
  feat_cols = NULL
) {
  df_obj <- if (is.data.frame(df)) df else NULL
  if (is.null(df_obj)) {
    if (is.null(group_vec) && !is.null(df)) {
      group_vec <- as.character(df)
    }
  } else if (is.null(group_vec)) {
    cand_cols <- c("taxon_clean", "species_standardized", "species_clean", "taxon", "scientificName", "given_label", "recorded_taxon")
    hit_cols <- intersect(cand_cols, colnames(df_obj))
    if (length(hit_cols) > 0L) {
      group_vec <- as.character(df_obj[[hit_cols[1L]]])
    } else if (exists("resolve_plot_group_column", mode = "function")) {
      group_vec <- resolve_plot_group_column(df_obj)
    } else {
      group_vec <- rep("Unassigned_Packera", nrow(df_obj))
    }
  }
  grp_chr <- trimws(as.character(group_vec))
  valid_m <- !is.na(grp_chr) & nzchar(grp_chr) & !(tolower(grp_chr) %in% c("unknown", "na", "none"))
  u_grps <- unique(grp_chr[valid_m])
  if (length(u_grps) == 0L) {
    u_grps <- "Unassigned_Packera"
  }

  n_g <- vapply(u_grps, function(g) sum(valid_m & grp_chr == g, na.rm = TRUE), integer(1))
  names(n_g) <- u_grps

  focal_chr <- trimws(as.character(focal_taxon)[1L])
  if (is.na(focal_chr) || !nzchar(focal_chr)) focal_chr <- "Packera dubia"

  exact_hit <- u_grps[tolower(u_grps) == tolower(focal_chr)]
  dubia_hit <- u_grps[grepl("^Packera\\s+dubia|^P\\.\\s*dubia", u_grps, ignore.case = TRUE)]
  focal_candidates <- unique(c(exact_hit, dubia_hit))

  if (length(focal_candidates) > 0L) {
    focal_candidates <- focal_candidates[order(-n_g[focal_candidates], focal_candidates)]
    slot1 <- focal_candidates[1L]
  } else {
    ord_n <- u_grps[order(-n_g[u_grps], u_grps)]
    slot1 <- ord_n[1L]
  }

  rem_grps <- setdiff(u_grps, slot1)
  if (length(rem_grps) > 0L) {
    d2_vec <- rep(Inf, length(rem_grps))
    names(d2_vec) <- rem_grps

    if (!is.null(df_obj) && nrow(df_obj) == length(grp_chr) && nrow(df_obj) >= 3L) {
      cand_f <- if (!is.null(feat_cols)) {
        intersect(as.character(feat_cols), colnames(df_obj))
      } else {
        grep("^(Sym_)?PC[1-5]$|^can[1-3]$|^PC[1-5]_corrected$", colnames(df_obj), value = TRUE)
      }
      if (length(cand_f) > 0L) {
        num_f <- cand_f[vapply(cand_f, function(cn) {
          v <- suppressWarnings(as.numeric(df_obj[[cn]]))
          sum(is.finite(v)) >= 3L && isTRUE(stats::sd(v, na.rm = TRUE) > 1e-7)
        }, logical(1))]
        if (length(num_f) > 0L) {
          num_f <- utils::head(num_f, 5L)
          X_mat <- do.call(cbind, lapply(num_f, function(cn) {
            v <- suppressWarnings(as.numeric(df_obj[[cn]]))
            v[!is.finite(v)] <- mean(v[is.finite(v)], na.rm = TRUE)
            v[!is.finite(v)] <- 0.0
            v
          }))
          m1 <- which(valid_m & grp_chr == slot1)
          if (length(m1) >= 1L) {
            mu1 <- colMeans(X_mat[m1, , drop = FALSE], na.rm = TRUE)
            W_cov <- tryCatch(stats::cov(X_mat[valid_m, , drop = FALSE]), error = function(e) NULL)
            W_inv <- NULL
            if (!is.null(W_cov) && all(is.finite(W_cov))) {
              W_reg <- W_cov + diag(max(1e-4, 0.05 * mean(diag(W_cov), na.rm = TRUE)), ncol(X_mat))
              W_inv <- tryCatch(solve(W_reg), error = function(e) NULL)
            }
            for (rg in rem_grps) {
              mg <- which(valid_m & grp_chr == rg)
              if (length(mg) >= 1L) {
                mug <- colMeans(X_mat[mg, , drop = FALSE], na.rm = TRUE)
                dv <- mug - mu1
                if (all(is.finite(dv))) {
                  if (!is.null(W_inv)) {
                    val_d2 <- as.numeric(t(dv) %*% W_inv %*% dv)
                    d2_vec[rg] <- if (is.finite(val_d2) && val_d2 >= 0) val_d2 else sum(dv^2)
                  } else {
                    d2_vec[rg] <- sum(dv^2)
                  }
                }
              }
            }
          }
        }
      }
    }
    has_d2 <- any(is.finite(d2_vec))
    if (has_d2) {
      rem_grps <- rem_grps[order(d2_vec[rem_grps], -n_g[rem_grps], rem_grps)]
    } else {
      rem_grps <- rem_grps[order(-n_g[rem_grps], rem_grps)]
    }
  }

  ordered_groups <- c(slot1, rem_grps)
  max_fg_int <- max(1L, as.integer(max_fg))
  fg_groups <- utils::head(ordered_groups, max_fg_int)
  bg_groups <- if (length(ordered_groups) > max_fg_int) utils::tail(ordered_groups, -max_fg_int) else character(0)

  list(
    ordered_groups = ordered_groups,
    fg_groups = fg_groups,
    bg_groups = bg_groups,
    is_multi_crowded = length(ordered_groups) > max_fg_int,
    palette = get_dynamic_packera_palette(ordered_groups, preserve_order = TRUE),
    shapes = get_dynamic_packera_shapes(ordered_groups, preserve_order = TRUE)
  )
}

format_italic_taxon_expr <- function(tx) {
  tx_chr <- trimws(as.character(tx)[1L])
  if (is.na(tx_chr) || !nzchar(tx_chr) || tx_chr %in% c("Unknown", "Unassigned", "Unassigned_Packera")) {
    return("Unassigned")
  }
  if (grepl("^([A-Z][a-z]*\\.?\\s+[a-z-]+)\\s*(\\[K[0-9]+\\])$", tx_chr)) {
    sp_part <- sub("^([A-Z][a-z]*\\.?\\s+[a-z-]+)\\s*(\\[K[0-9]+\\])$", "\\1", tx_chr)
    k_part  <- sub("^([A-Z][a-z]*\\.?\\s+[a-z-]+)\\s*(\\[K[0-9]+\\])$", "\\2", tx_chr)
    return(bquote(italic(.(sp_part)) ~ .(k_part)))
  }
  if (grepl("^[A-Z][a-z]*\\.?\\s+(?:[x×]\\s*)?[a-z-]+$", tx_chr)) {
    return(bquote(italic(.(tx_chr))))
  }
  tx_chr
}

packera_italic_expressions <- function(taxa_vec) {
  as.expression(lapply(as.character(taxa_vec), format_italic_taxon_expr))
}

resolve_plot_group_column <- function(df, cda_res = NULL, gmm_res = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  n <- nrow(df)
  raw_sp <- if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("taxon" %in% names(df)) {
    as.character(df$taxon)
  } else if ("scientificName" %in% names(df)) {
    as.character(df$scientificName)
  } else if ("species_raw" %in% names(df)) {
    as.character(df$species_raw)
  } else {
    rep(NA_character_, n)
  }
  std_sp <- standardize_taxon_name(raw_sp)
  cand_sp <- ifelse(!is.na(std_sp), std_sp, ifelse(!is.na(raw_sp) & nzchar(trimws(raw_sp)), trimws(raw_sp), NA_character_))
  valid_sp <- cand_sp[!is.na(cand_sp) & nzchar(cand_sp) & !(tolower(cand_sp) %in% c("unknown", "na", "none"))]
  fallback_sp <- if (length(valid_sp) > 0L) valid_sp[1L] else "Unassigned_Packera"
  sp_vec <- ifelse(!is.na(cand_sp) & nzchar(cand_sp), cand_sp, fallback_sp)
  sp_valid <- unique(sp_vec[!is.na(sp_vec) & nzchar(sp_vec) & sp_vec != "Unknown"])

  clean_cand <- function(v) {
    if (is.null(v) || length(v) != n) return(NULL)
    vc <- trimws(as.character(v))
    vc[is.na(vc) | !nzchar(vc)] <- NA_character_
    vc
  }
  has_valid_multi <- function(vc) {
    !is.null(vc) && length(unique(vc[!is.na(vc) & vc != "Unknown"])) >= 2L
  }
  is_subcluster_vec <- function(vc) {
    if (!has_valid_multi(vc)) return(FALSE)
    vv <- unique(vc[!is.na(vc) & vc != "Unknown"])
    any(grepl("Morphotype_K|\\[K[0-9]+\\]|Ecotype|Subcluster|Region_", vv, ignore.case = TRUE)) ||
      !all(grepl("^Packera\\s+", vv))
  }

  cand_cda <- if (!is.null(cda_res) && !is.null(cda_res$grouping_col)) clean_cand(cda_res$grouping_col) else NULL
  cand_grp <- if ("grouping_col" %in% names(df)) clean_cand(df$grouping_col) else NULL
  cand_ana <- if ("analysis_group" %in% names(df)) clean_cand(df$analysis_group) else NULL
  cand_gmm <- NULL
  if ("gmm_cluster" %in% names(df) && length(unique(stats::na.omit(df$gmm_cluster))) >= 2L) {
    cand_gmm <- paste0("Morphotype_K", ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster)))
  } else if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n &&
             length(unique(stats::na.omit(gmm_res$classification))) >= 2L) {
    cand_gmm <- paste0("Morphotype_K", ifelse(is.na(gmm_res$classification), 1L, as.integer(gmm_res$classification)))
  }

  use_subclusters <- (length(sp_valid) < 2L) ||
    is_subcluster_vec(cand_cda) ||
    is_subcluster_vec(cand_grp) ||
    is_subcluster_vec(cand_ana)

  if (use_subclusters) {
    for (cand in list(cand_cda, cand_grp, cand_ana, cand_gmm)) {
      if (has_valid_multi(cand)) {
        if (length(sp_valid) >= 2L && !is_subcluster_vec(cand)) next
        def_lbl <- stats::na.omit(cand)[1]
        return(ifelse(!is.na(cand), cand, def_lbl))
      }
    }
    for (col_fb in c("cda_predicted_taxon", "edaphic_best_fit_taxon", "recommended_determination")) {
      if (col_fb %in% names(df)) {
        cf <- clean_cand(df[[col_fb]])
        if (is_subcluster_vec(cf)) {
          def_lbl <- stats::na.omit(cf)[1]
          return(ifelse(!is.na(cf), cf, def_lbl))
        }
      }
    }
    if ("gmm_cluster" %in% names(df) && length(unique(stats::na.omit(df$gmm_cluster))) >= 2L) {
      return(paste0("Morphotype_K", ifelse(is.na(df$gmm_cluster), 1L, as.integer(df$gmm_cluster))))
    }
  }
  sp_vec
}

resolve_active_figures_dir <- function(explicit_dir = NULL, reference_path = NULL) {
  if (!is.null(explicit_dir) && nzchar(as.character(explicit_dir))) {
    fdir <- as.character(explicit_dir)
    if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    return(fdir)
  }
  env_run <- Sys.getenv("PACKERA_RUN_DIR", unset = "")
  if (nzchar(env_run)) {
    fdir <- file.path(env_run, "outputs", "figures")
    if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    return(fdir)
  }
  if (!is.null(reference_path) && nzchar(as.character(reference_path))) {
    ref_p <- as.character(reference_path)
    if (grepl("\\.(pdf|png)$", ref_p, ignore.case = TRUE)) {
      fdir <- dirname(ref_p)
    } else {
      # e.g., <run_dir>/data/tables/file.csv -> <run_dir>/outputs/figures
      parent1 <- dirname(ref_p)
      parent2 <- dirname(parent1)
      if (basename(parent1) == "tables" && basename(parent2) == "data") {
        fdir <- file.path(dirname(parent2), "outputs", "figures")
      } else {
        fdir <- file.path(parent1, "figures")
      }
    }
    if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    return(fdir)
  }
  fdir <- file.path("outputs", "figures")
  if (!dir.exists(fdir)) dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
  return(fdir)
}

atomic_save_dual_plot <- function(draw_fn, figures_dir, stem, width = 11, height = 5.5, dpi = 300) {
  if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  }
  pdf_legacy <- file.path(figures_dir, paste0(stem, ".pdf"))
  if (file.exists(pdf_legacy)) unlink(pdf_legacy, force = TRUE)
  png_final <- file.path(figures_dir, paste0(stem, ".png"))
  png_tmp <- file.path(figures_dir, paste0(".", stem, "_", Sys.getpid(), ".png.tmp"))

  opened_ragg <- FALSE
  if (requireNamespace("ragg", quietly = TRUE)) {
    opened_ragg <- tryCatch({
      ragg::agg_png(png_tmp, width = width, height = height, units = "in", res = dpi)
      TRUE
    }, error = function(e) FALSE)
  }
  if (!opened_ragg) {
    grDevices::png(png_tmp, width = width, height = height, units = "in", res = dpi)
  }
  tryCatch({
    draw_fn()
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(png_tmp)) {
    file.rename(png_tmp, png_final)
  }
  invisible(list(pdf = NULL, png = png_final))
}

# ==============================================================================
# Inverse Elliptic Fourier Synthesis & Reusable 2D Ordination Helpers (Base R)
# Always enforces strictly unfilled/hollow markers (pch = 1, 2, 0, 5 or bg = NA)
# ==============================================================================
reconstruct_efa_outline <- function(harm_vec, nb_pts = 120L) {
  if (is.null(harm_vec) || length(harm_vec) == 0) {
    theta <- seq(0, 2 * pi, length.out = nb_pts)
    return(cbind(x = 0.42 * cos(theta), y = 1.0 * sin(theta)))
  }
  nms <- names(harm_vec)
  v_num <- suppressWarnings(as.numeric(harm_vec))
  v_num[!is.finite(v_num)] <- 0.0
  names(v_num) <- nms

  if (!is.null(nms) && any(grepl("^[ABCD]_?[0-9]+$", nms, ignore.case = TRUE))) {
    clean_nms <- toupper(gsub("_", "", nms))
    names(v_num) <- clean_nms
    h_idx <- suppressWarnings(as.integer(gsub("[^0-9]", "", clean_nms)))
    K <- max(c(1L, h_idx[is.finite(h_idx)]), na.rm = TRUE)
    A <- numeric(K); B <- numeric(K); C <- numeric(K); D <- numeric(K)
    for (k in seq_len(K)) {
      if (paste0("A", k) %in% clean_nms) A[k] <- v_num[paste0("A", k)][1]
      if (paste0("B", k) %in% clean_nms) B[k] <- v_num[paste0("B", k)][1]
      if (paste0("C", k) %in% clean_nms) C[k] <- v_num[paste0("C", k)][1]
      if (paste0("D", k) %in% clean_nms) D[k] <- v_num[paste0("D", k)][1]
    }
  } else {
    K <- max(1L, floor(length(v_num) / 4L))
    A <- v_num[seq(1, by = 4, length.out = K)]
    B <- v_num[seq(2, by = 4, length.out = K)]
    C <- v_num[seq(3, by = 4, length.out = K)]
    D <- v_num[seq(4, by = 4, length.out = K)]
  }

  if (sum(A^2 + B^2 + C^2 + D^2) < 1e-12) {
    A[1] <- 0.0; B[1] <- -0.42; C[1] <- 0.95; D[1] <- 0.0
  }
  t_seq <- seq(0, 2 * pi, length.out = nb_pts)
  x <- numeric(nb_pts)
  y <- numeric(nb_pts)
  for (n in seq_len(K)) {
    nt <- n * t_seq
    x <- x + A[n] * cos(nt) + B[n] * sin(nt)
    y <- y + C[n] * cos(nt) + D[n] * sin(nt)
  }
  # Orient major blade axis vertically (apex up, petiole base down)
  span_x <- diff(range(x, na.rm = TRUE))
  span_y <- diff(range(y, na.rm = TRUE))
  if (span_x > 1.05 * span_y) {
    tmp <- x; x <- y; y <- -tmp
  }
  # Ensure narrower petiole end is at bottom (y < 0) and broader lamina is at top
  y_centered <- y - mean(y, na.rm = TRUE)
  x_centered <- x - mean(x, na.rm = TRUE)
  top_w <- stats::sd(x_centered[y_centered > 0.25 * max(y_centered, na.rm = TRUE)], na.rm = TRUE)
  bot_w <- stats::sd(x_centered[y_centered < 0.25 * min(y_centered, na.rm = TRUE)], na.rm = TRUE)
  if (is.finite(top_w) && is.finite(bot_w) && bot_w > 1.25 * top_w) {
    y_centered <- -y_centered
  }
  cbind(x = x_centered, y = y_centered)
}

extract_standardized_harmonic_matrix <- function(df, max_k = NULL) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (is.null(max_k) || is.na(max_k)) {
    h_cols <- grep("^[ABCD][0-9]+$", colnames(df), value = TRUE)
    h_nums <- suppressWarnings(as.integer(gsub("[^0-9]", "", h_cols)))
    max_k <- if (any(is.finite(h_nums))) max(12L, min(24L, max(h_nums, na.rm = TRUE))) else 12L
  } else {
    max_k <- max(6L, min(24L, as.integer(max_k)))
  }
  interleaved_nms <- as.vector(outer(c("A", "B", "C", "D"), seq_len(max_k), paste0))
  avail <- intersect(interleaved_nms, colnames(df))
  if (length(avail) < 4L) return(NULL)
  mat <- as.matrix(df[, avail, drop = FALSE])
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- 0.0
  mat
}

draw_outline_in_box <- function(
  xy,
  cx,
  cy,
  w_box,
  h_box,
  col_border = "grey55",
  col_fill = grDevices::adjustcolor("grey85", alpha.f = 0.35),
  lwd = 1.0,
  lty = 1,
  shared_scale = NULL,
  shared_center = NULL
) {
  if (is.null(xy) || nrow(xy) < 3) return(invisible(NULL))
  if (!is.null(shared_scale) && length(shared_scale) == 1L && is.finite(shared_scale) &&
      !is.null(shared_center) && length(shared_center) >= 2L && all(is.finite(shared_center[1:2]))) {
    px <- cx + (xy[, 1] - shared_center[1]) * shared_scale
    py <- cy + (xy[, 2] - shared_center[2]) * shared_scale
  } else {
    rx <- range(xy[, 1], na.rm = TRUE)
    ry <- range(xy[, 2], na.rm = TRUE)
    sx <- max(diff(rx), 1e-6)
    sy <- max(diff(ry), 1e-6)
    scale_f <- min(w_box / sx, h_box / sy) * 0.86
    px <- cx + (xy[, 1] - mean(rx)) * scale_f
    py <- cy + (xy[, 2] - mean(ry)) * scale_f
  }
  if (!is.null(col_fill) && !is.na(col_fill)) {
    graphics::polygon(px, py, col = col_fill, border = col_border, lwd = lwd, lty = lty)
  } else {
    graphics::lines(c(px, px[1]), c(py, py[1]), col = col_border, lwd = lwd, lty = lty)
  }
  invisible(cbind(px, py))
}

.compute_outline_grid_layout <- function(n_disp, draw_dividers = TRUE) {
  n_d <- max(1L, min(9L, as.integer(n_disp)))
  if (n_d == 1L) {
    nr <- 1L; nc <- 1L
    cx_vec <- 1.0; cy_vec <- 1.0
  } else if (n_d == 2L) {
    nr <- 1L; nc <- 2L
    if (draw_dividers) graphics::abline(v = 1.0, col = "grey88", lty = 2)
    cx_vec <- c(0.5, 1.5); cy_vec <- c(1.0, 1.0)
  } else if (n_d == 3L) {
    nr <- 2L; nc <- 2L
    if (draw_dividers) {
      graphics::abline(h = 1.0, col = "grey88", lty = 2)
      graphics::segments(1.0, 1.0, 1.0, 2.0, col = "grey88", lty = 2)
    }
    cx_vec <- c(0.5, 1.5, 1.0); cy_vec <- c(1.48, 1.48, 0.48)
  } else if (n_d == 4L) {
    nr <- 2L; nc <- 2L
    if (draw_dividers) graphics::abline(h = 1.0, v = 1.0, col = "grey88", lty = 2)
    cx_vec <- c(0.5, 1.5, 0.5, 1.5); cy_vec <- c(1.48, 1.48, 0.48, 0.48)
  } else if (n_d <= 6L) {
    nr <- 2L; nc <- 3L
    if (draw_dividers) {
      graphics::abline(h = 1.0, col = "grey88", lty = 2)
      graphics::abline(v = c(2.0 / 3.0, 4.0 / 3.0), col = "grey88", lty = 2)
    }
    cx_all <- rep(c(1.0 / 3.0, 1.0, 5.0 / 3.0), times = 2L)
    cy_all <- rep(c(1.48, 0.48), each = 3L)
    cx_vec <- cx_all[seq_len(n_d)]; cy_vec <- cy_all[seq_len(n_d)]
  } else {
    nr <- 3L; nc <- 3L
    if (draw_dividers) {
      graphics::abline(h = c(2.0 / 3.0, 4.0 / 3.0), col = "grey88", lty = 2)
      graphics::abline(v = c(2.0 / 3.0, 4.0 / 3.0), col = "grey88", lty = 2)
    }
    cx_all <- rep(c(1.0 / 3.0, 1.0, 5.0 / 3.0), times = 3L)
    cy_all <- rep(c(5.0 / 3.0, 1.0, 1.0 / 3.0), each = 3L)
    cx_vec <- cx_all[seq_len(n_d)]; cy_vec <- cy_all[seq_len(n_d)]
  }
  cell_scale <- if (n_d <= 4L) 1.0 else if (n_d <= 6L) 0.72 else 0.58
  list(
    nr = nr, nc = nc,
    cx = cx_vec, cy = cy_vec,
    cell_scale = cell_scale,
    title_dy = if (nr == 3L) 0.26 else 0.41,
    sub_dy   = if (nr == 3L) 0.27 else 0.43,
    body_dy  = if (nr == 3L) 0.02 else 0.05,
    cex_title = if (n_d <= 4L) 0.78 else if (n_d <= 6L) 0.68 else 0.60,
    cex_sub   = if (n_d <= 4L) 0.65 else if (n_d <= 6L) 0.56 else 0.50
  )
}

repel_2d_labels <- function(
  x0,
  y0,
  labels,
  xlim,
  ylim,
  cex = 0.72,
  n_iter = 25L,
  draw_segments = TRUE,
  seg_col = "grey45",
  seg_lwd = 0.75,
  anchor_x = NULL,
  anchor_y = NULL
) {
  x_init <- suppressWarnings(as.numeric(x0))
  y_init <- suppressWarnings(as.numeric(y0))
  n <- min(length(x_init), length(y_init))
  if (n == 0L) return(list(x = numeric(0), y = numeric(0)))
  x_lbl <- x_init[seq_len(n)]
  y_lbl <- y_init[seq_len(n)]
  ax <- if (!is.null(anchor_x) && length(anchor_x) >= n) suppressWarnings(as.numeric(anchor_x[seq_len(n)])) else x_lbl
  ay <- if (!is.null(anchor_y) && length(anchor_y) >= n) suppressWarnings(as.numeric(anchor_y[seq_len(n)])) else y_lbl
  lbl_chr <- as.character(rep_len(labels, n))

  span_x <- max(diff(xlim), 1e-6)
  span_y <- max(diff(ylim), 1e-6)
  x_min <- xlim[1] + 0.04 * span_x
  x_max <- xlim[2] - 0.04 * span_x
  y_min <- ylim[1] + 0.04 * span_y
  y_max <- ylim[2] - 0.04 * span_y

  sw <- tryCatch(graphics::strwidth(lbl_chr, cex = cex), error = function(e) rep(0.0, n))
  sh <- tryCatch(graphics::strheight(lbl_chr, cex = cex), error = function(e) rep(0.0, n))
  sw[!is.finite(sw) | sw <= 0] <- 0.012 * span_x * pmax(nchar(lbl_chr), 3L)
  sh[!is.finite(sh) | sh <= 0] <- 0.040 * span_y
  hw <- 0.56 * sw + 0.012 * span_x
  hh <- 0.64 * sh + 0.010 * span_y

  x_lbl <- pmin(x_max, pmax(x_min, x_lbl))
  y_lbl <- pmin(y_max, pmax(y_min, y_lbl))

  if (n >= 2L && n_iter >= 1L) {
    for (iter in seq_len(as.integer(n_iter))) {
      for (i in 1:(n - 1L)) {
        for (j in (i + 1L):n) {
          if (!is.finite(x_lbl[i]) || !is.finite(y_lbl[i]) || !is.finite(x_lbl[j]) || !is.finite(y_lbl[j])) next
          dx <- x_lbl[i] - x_lbl[j]
          dy <- y_lbl[i] - y_lbl[j]
          ov_x <- (hw[i] + hw[j]) - abs(dx)
          ov_y <- (hh[i] + hh[j]) - abs(dy)
          if (ov_x > 0 && ov_y > 0) {
            if (abs(dx) < 1e-9 && abs(dy) < 1e-9) {
              ang <- (i - j) * 0.75
              dx <- cos(ang) * 1e-4 * span_x
              dy <- sin(ang) * 1e-4 * span_y
            }
            sx <- if (dx >= 0) 1.0 else -1.0
            sy <- if (dy >= 0) 1.0 else -1.0
            px_step <- 0.34 * ov_x * sx
            py_step <- 0.44 * ov_y * sy
            x_lbl[i] <- x_lbl[i] + px_step
            x_lbl[j] <- x_lbl[j] - px_step
            y_lbl[i] <- y_lbl[i] + py_step
            y_lbl[j] <- y_lbl[j] - py_step
          }
        }
      }
      x_lbl <- pmin(x_max, pmax(x_min, x_lbl))
      y_lbl <- pmin(y_max, pmax(y_min, y_lbl))
    }
  }

  if (isTRUE(draw_segments)) {
    move_d <- sqrt((x_lbl - ax)^2 + (y_lbl - ay)^2)
    thresh <- 0.04 * min(span_x, span_y)
    moved_idx <- which(is.finite(move_d) & (move_d > thresh))
    if (length(moved_idx) > 0L) {
      graphics::segments(ax[moved_idx], ay[moved_idx], x_lbl[moved_idx], y_lbl[moved_idx], col = seg_col, lwd = seg_lwd)
    }
  }

  list(x = x_lbl, y = y_lbl)
}

draw_2d_confidence_ellipse <- function(
  x,
  y,
  col = "#0072B2",
  level = 0.95,
  lwd = 2.0,
  lty = 2,
  n_pts = 100,
  weights = NULL
) {
  x_num <- suppressWarnings(as.numeric(x))
  y_num <- suppressWarnings(as.numeric(y))
  valid <- is.finite(x_num) & is.finite(y_num)
  if (sum(valid) < 3) return(invisible(NULL))

  xv <- x_num[valid]
  yv <- y_num[valid]
  cx_el <- mean(xv)
  cy_el <- mean(yv)
  cov_mat <- NULL
  if (!is.null(weights) && length(weights) == length(x_num)) {
    w_sub <- suppressWarnings(as.numeric(weights[valid]))
    w_ok <- is.finite(w_sub) & (w_sub > 0)
    if (sum(w_ok) >= 3L) {
      w_pos <- w_sub[w_ok]
      n_eff <- (sum(w_pos)^2) / max(sum(w_pos^2), 1e-12)
      if (is.finite(n_eff) && n_eff >= 2.2) {
        cwt <- tryCatch(stats::cov.wt(cbind(xv[w_ok], yv[w_ok]), wt = w_pos, method = "unbiased"), error = function(e) NULL)
        if (!is.null(cwt) && !is.null(cwt$cov) && all(is.finite(cwt$cov))) {
          cov_mat <- cwt$cov
          cx_el <- as.numeric(cwt$center[1])
          cy_el <- as.numeric(cwt$center[2])
        }
      }
    }
  }
  if (is.null(cov_mat)) {
    cov_mat <- tryCatch(stats::cov(cbind(xv, yv)), error = function(e) NULL)
  }
  if (is.null(cov_mat) || any(!is.finite(cov_mat)) || max(diag(cov_mat)) < 1e-10) {
    return(invisible(NULL))
  }

  eig <- tryCatch(eigen(cov_mat, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig)) return(invisible(NULL))

  evals <- pmax(eig$values, 1e-10)
  radius <- sqrt(stats::qchisq(level, df = 2))
  theta <- seq(0, 2 * pi, length.out = n_pts)
  circle <- rbind(cos(theta), sin(theta))
  transform_mat <- eig$vectors %*% diag(sqrt(evals), nrow = 2, ncol = 2)
  ellipse_pts <- sweep(radius * (transform_mat %*% circle), 1, c(cx_el, cy_el), "+")
  graphics::lines(ellipse_pts[1, ], ellipse_pts[2, ], col = col, lwd = lwd, lty = lty)
  invisible(ellipse_pts)
}

plot_2d_ordination_scatter <- function(
  x,
  y,
  col = "#0072B2",
  pch = 1L,
  cex = 0.85,
  lwd = 1.0,
  xlab = "Axis 1",
  ylab = "Axis 2",
  main = "",
  subtitle = NULL,
  group_vec = NULL,
  group_cols = NULL,
  ellipse_level = 0.95,
  ellipse_weights = NULL,
  centroids = NULL,
  centroid_labels = NULL,
  legend_args = NULL,
  bg_draw_fn = NULL,
  biplot_vectors = NULL,
  biplot_labels = NULL,
  is_type_vec = NULL
) {
  x_num <- suppressWarnings(as.numeric(x))
  y_num <- suppressWarnings(as.numeric(y))
  n_pts <- max(length(x_num), length(y_num))
  if (n_pts == 0) return(invisible(NULL))

  x_num[is.na(x_num) | !is.finite(x_num)] <- 0.0
  y_num[is.na(y_num) | !is.finite(y_num)] <- 0.0

  col_vec <- rep_len(as.character(col), n_pts)
  pch_vec <- rep_len(as.integer(pch), n_pts)
  # Enforce strictly unfilled/hollow symbols even if caller passes solid pch codes
  pch_vec[pch_vec %in% c(16L, 19L, 20L)] <- 1L
  pch_vec[pch_vec == 17L] <- 2L
  pch_vec[pch_vec == 18L] <- 5L
  pch_vec[pch_vec == 15L] <- 0L
  cex_vec <- rep_len(as.numeric(cex), n_pts)
  lwd_vec <- rep_len(as.numeric(lwd), n_pts)

  pad_range <- function(v) {
    r <- range(v, na.rm = TRUE)
    if (!all(is.finite(r))) return(c(-1, 1))
    span <- r[2] - r[1]
    if (span < 1e-6) return(c(r[1] - 1.0, r[2] + 1.0))
    c(r[1] - 0.14 * span, r[2] + 0.14 * span)
  }

  all_x <- x_num
  all_y <- y_num
  c_mat <- NULL
  if (!is.null(centroids)) {
    c_df <- as.matrix(centroids)
    if (nrow(c_df) > 0 && ncol(c_df) >= 2) {
      cx <- suppressWarnings(as.numeric(c_df[, 1]))
      cy <- suppressWarnings(as.numeric(c_df[, 2]))
      cx[!is.finite(cx)] <- 0.0
      cy[!is.finite(cy)] <- 0.0
      c_mat <- cbind(cx, cy)
      rownames(c_mat) <- rownames(c_df)
      all_x <- c(all_x, cx)
      all_y <- c(all_y, cy)
    }
  }

  xlim <- pad_range(all_x)
  ylim <- pad_range(all_y)
  has_sub <- !is.null(subtitle) && length(subtitle) > 0L && !is.na(subtitle[1]) && nzchar(trimws(as.character(subtitle[1])))

  grp_chr <- if (!is.null(group_vec) && length(group_vec) == n_pts) trimws(as.character(group_vec)) else NULL
  u_grps <- if (!is.null(grp_chr)) unique(grp_chr[!is.na(grp_chr) & nzchar(grp_chr) & grp_chr != "Unknown"]) else character(0)

  # G = 2 Mirrored/Overlaid Canonical Histogram Mode (MorphoTools2::plotHist convention)
  is_can2_degen <- (length(u_grps) == 2L) &&
    grepl("^Can2", as.character(ylab), ignore.case = TRUE) &&
    (isTRUE(stats::sd(y_num, na.rm = TRUE) < 1e-5) || grepl("\\(0\\.0%\\)", as.character(ylab)))

  if (is_can2_degen) {
    brks <- seq(xlim[1], xlim[2], length.out = 22L)
    h_list <- list()
    d_list <- list()
    max_dens <- 0.25
    for (g_val in u_grps) {
      xg <- x_num[!is.na(grp_chr) & grp_chr == g_val]
      hg <- tryCatch(graphics::hist(xg, breaks = brks, plot = FALSE), error = function(e) NULL)
      dg <- if (length(xg) >= 2L && isTRUE(stats::sd(xg, na.rm = TRUE) > 1e-6)) {
        tryCatch(stats::density(xg, from = xlim[1], to = xlim[2]), error = function(e) NULL)
      } else NULL
      h_list[[g_val]] <- hg
      d_list[[g_val]] <- dg
      if (!is.null(hg) && length(hg$density) > 0L) max_dens <- max(max_dens, max(hg$density, na.rm = TRUE))
      if (!is.null(dg) && length(dg$y) > 0L) max_dens <- max(max_dens, max(dg$y, na.rm = TRUE))
    }
    ylim_hist <- c(0.0, max_dens * 1.28)
    graphics::plot(
      xlim, ylim_hist,
      type = "n",
      xlim = xlim,
      ylim = ylim_hist,
      xlab = xlab,
      ylab = "Canonical Score Density (Can1)",
      main = if (has_sub) "" else main,
      font.main = 2
    )
    if (has_sub) {
      if (!is.null(main) && nzchar(trimws(as.character(main)))) {
        graphics::title(main = main, line = 1.45, font.main = 2)
      }
      graphics::mtext(trimws(as.character(subtitle[1])), side = 3, line = 0.28, cex = 0.73, col = "grey20", font = 1)
    }
    graphics::grid(col = "grey88", lty = "solid")
    for (g_val in u_grps) {
      m <- !is.na(grp_chr) & (grp_chr == g_val)
      g_col <- if (!is.null(group_cols) && g_val %in% names(group_cols)) group_cols[[g_val]] else col_vec[which(m)[1L]]
      hg <- h_list[[g_val]]
      dg <- d_list[[g_val]]
      if (!is.null(hg) && length(hg$density) > 0L) {
        graphics::rect(
          hg$breaks[-length(hg$breaks)], 0.0, hg$breaks[-1L], hg$density,
          col = grDevices::adjustcolor(g_col, alpha.f = 0.28),
          border = grDevices::adjustcolor(g_col, alpha.f = 0.75),
          lwd = 1.1
        )
      }
      if (!is.null(dg)) {
        graphics::lines(dg$x, dg$y, col = g_col, lwd = 2.2)
      }
      mu_g <- mean(x_num[m], na.rm = TRUE)
      if (is.finite(mu_g)) {
        graphics::abline(v = mu_g, col = g_col, lty = 2, lwd = 2.0)
        graphics::text(mu_g, max_dens * 1.15, labels = sub("^Packera\\s+", "P. ", g_val), col = g_col, cex = 0.76, font = 4)
      }
    }
    if (!is.null(is_type_vec) && length(is_type_vec) > 0L) {
      t_bool <- as.logical(rep_len(is_type_vec, n_pts))
      t_bool[is.na(t_bool)] <- FALSE
      type_idx <- which(t_bool)
      if (length(type_idx) > 0L) {
        graphics::arrows(
          x_num[type_idx], max_dens * 0.20, x_num[type_idx], 0.0,
          length = 0.07, lwd = 2.0, col = col_vec[type_idx]
        )
        graphics::points(x_num[type_idx], rep(max_dens * 0.20, length(type_idx)), pch = 5L, col = "#111111", cex = 1.25, lwd = 1.6)
      }
    }
    if (!is.null(legend_args) && is.list(legend_args)) {
      do.call(graphics::legend, legend_args)
    }
    return(invisible(list(xlim = xlim, ylim = ylim_hist)))
  }

  graphics::plot(
    x_num, y_num,
    type = "n",
    xlim = xlim,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    main = if (has_sub) "" else main,
    font.main = 2
  )
  if (has_sub) {
    if (!is.null(main) && nzchar(trimws(as.character(main)))) {
      graphics::title(main = main, line = 1.45, font.main = 2)
    }
    graphics::mtext(trimws(as.character(subtitle[1])), side = 3, line = 0.28, cex = 0.73, col = "grey20", font = 1)
  }
  graphics::grid(col = "grey88", lty = "solid")
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")

  if (is.function(bg_draw_fn)) {
    tryCatch(bg_draw_fn(xlim, ylim), error = function(e) NULL)
  }

  # G > 6 Focal-Plus-Top-k Context Layering
  tier_ord <- NULL
  fg_grps <- u_grps
  bg_grps <- character(0)
  if (length(u_grps) > 6L) {
    tier_ord <- resolve_taxon_display_tier(
      data.frame(PC1 = x_num, PC2 = y_num, stringsAsFactors = FALSE),
      group_vec = grp_chr,
      max_fg = 6L,
      feat_cols = c("PC1", "PC2")
    )
    fg_grps <- tier_ord$fg_groups
    bg_grps <- tier_ord$bg_groups

    for (bg_val in bg_grps) {
      m_bg <- which(!is.na(grp_chr) & (grp_chr == bg_val))
      if (length(m_bg) >= 3L) {
        ch_idx <- tryCatch(grDevices::chull(x_num[m_bg], y_num[m_bg]), error = function(e) NULL)
        if (!is.null(ch_idx) && length(ch_idx) >= 3L) {
          ch_pts <- m_bg[c(ch_idx, ch_idx[1L])]
          graphics::polygon(
            x_num[ch_pts], y_num[ch_pts],
            col = NA,
            border = grDevices::adjustcolor("#888888", alpha.f = 0.45),
            lty = 2, lwd = 1.0
          )
        }
      }
    }
    bg_mask <- !is.na(grp_chr) & (grp_chr %in% bg_grps)
    if (any(bg_mask)) {
      graphics::points(
        x_num[bg_mask], y_num[bg_mask],
        col = grDevices::adjustcolor("#999999", alpha.f = 0.32),
        pch = 1L, bg = NA, cex = 0.58, lwd = 0.85
      )
    }
    fg_mask <- !bg_mask
    if (any(fg_mask)) {
      fg_alpha <- adaptive_point_alpha(sum(fg_mask))
      fg_col_vec <- col_vec[fg_mask]
      fg_pch_vec <- pch_vec[fg_mask]
      for (fg_val in fg_grps) {
        mf <- which(grp_chr[fg_mask] == fg_val)
        if (length(mf) > 0L) {
          c_use <- if (!is.null(group_cols) && fg_val %in% names(group_cols)) group_cols[[fg_val]] else tier_ord$palette[[fg_val]]
          fg_col_vec[mf] <- grDevices::adjustcolor(c_use, alpha.f = fg_alpha)
          fg_pch_vec[mf] <- tier_ord$shapes[[fg_val]]
        }
      }
      graphics::points(
        x_num[fg_mask], y_num[fg_mask],
        col = fg_col_vec,
        pch = fg_pch_vec,
        bg = NA,
        cex = cex_vec[fg_mask],
        lwd = lwd_vec[fg_mask]
      )
    }
  } else {
    graphics::points(
      x_num, y_num,
      col = col_vec,
      pch = pch_vec,
      bg = NA,
      cex = cex_vec,
      lwd = lwd_vec
    )
  }

  if (!is.null(is_type_vec) && length(is_type_vec) > 0L) {
    t_bool <- as.logical(rep_len(is_type_vec, n_pts))
    t_bool[is.na(t_bool)] <- FALSE
    type_idx <- which(t_bool)
    if (length(type_idx) > 0L) {
      graphics::points(
        x_num[type_idx], y_num[type_idx],
        col = col_vec[type_idx], pch = 5L, bg = NA, cex = 1.45, lwd = 2.0
      )
      graphics::points(
        x_num[type_idx], y_num[type_idx],
        col = "#111111", pch = 1L, bg = NA, cex = 1.85, lwd = 1.4
      )
    }
  }

  if (!is.null(grp_chr) && length(fg_grps) > 0L) {
    for (g_val in fg_grps) {
      m <- !is.na(grp_chr) & (grp_chr == g_val)
      if (sum(m) >= 3) {
        g_col <- if (!is.null(group_cols) && g_val %in% names(group_cols)) {
          group_cols[[g_val]]
        } else if (!is.null(tier_ord) && g_val %in% names(tier_ord$palette)) {
          tier_ord$palette[[g_val]]
        } else {
          col_vec[which(m)[1]]
        }
        w_m <- if (!is.null(ellipse_weights) && length(ellipse_weights) == n_pts) ellipse_weights[m] else NULL
        draw_2d_confidence_ellipse(x_num[m], y_num[m], col = g_col, level = ellipse_level, lwd = 2.0, lty = 1, weights = w_m)
      }
    }
  }

  if (!is.null(biplot_vectors)) {
    b_mat <- tryCatch(as.matrix(biplot_vectors), error = function(e) NULL)
    if (!is.null(b_mat) && nrow(b_mat) > 0 && ncol(b_mat) >= 2) {
      bx <- suppressWarnings(as.numeric(b_mat[, 1]))
      by <- suppressWarnings(as.numeric(b_mat[, 2]))
      bx[!is.finite(bx)] <- 0.0; by[!is.finite(by)] <- 0.0
      vec_len <- sqrt(bx^2 + by^2)
      max_len <- max(vec_len, 1e-6)
      rad_scale <- 0.36 * min(diff(xlim), diff(ylim)) / max_len
      vx <- bx * rad_scale
      vy <- by * rad_scale
      top_idx <- utils::head(order(vec_len, decreasing = TRUE), 4L)
      keep_v <- seq_along(vec_len) %in% top_idx & (sqrt(vx^2 + vy^2) > 0.06 * min(diff(xlim), diff(ylim)))
      if (any(keep_v)) {
        graphics::arrows(
          0, 0, vx[keep_v], vy[keep_v],
          length = 0.07, lwd = 1.8, col = "#222222"
        )
        v_lbls <- if (!is.null(biplot_labels) && length(biplot_labels) == nrow(b_mat)) {
          biplot_labels[keep_v]
        } else if (!is.null(rownames(b_mat))) {
          rownames(b_mat)[keep_v]
        } else {
          paste0("V", which(keep_v))
        }
        tx_x0 <- vx[keep_v] * 1.14
        tx_y0 <- vy[keep_v] * 1.14
        rep_v <- repel_2d_labels(
          x0 = tx_x0, y0 = tx_y0, labels = v_lbls,
          xlim = xlim, ylim = ylim, cex = 0.72, n_iter = 25L,
          draw_segments = TRUE, seg_col = "grey45", seg_lwd = 0.75,
          anchor_x = vx[keep_v], anchor_y = vy[keep_v]
        )
        graphics::text(rep_v$x, rep_v$y, labels = v_lbls, cex = 0.72, font = 2, col = "#111111")
      }
    }
  }

  if (!is.null(c_mat) && nrow(c_mat) > 0) {
    lbls <- if (!is.null(centroid_labels)) as.character(centroid_labels) else rownames(c_mat)
    c_use_mat <- c_mat
    c_use_lbls <- lbls
    if (!is.null(c_use_lbls) && length(c_use_lbls) == nrow(c_use_mat)) {
      rownames(c_use_mat) <- c_use_lbls
    }
    norm_lbl_key <- function(v) gsub("[_[:space:]]+", " ", tolower(trimws(as.character(v))))
    if (length(u_grps) > 6L && length(fg_grps) > 0L && !is.null(c_use_lbls) && length(c_use_lbls) == nrow(c_use_mat)) {
      short_fg <- sub("^Packera\\s+", "P. ", fg_grps)
      keep_c <- (norm_lbl_key(c_use_lbls) %in% norm_lbl_key(fg_grps)) | (norm_lbl_key(c_use_lbls) %in% norm_lbl_key(short_fg))
      if (!any(keep_c)) keep_c <- seq_len(min(6L, nrow(c_use_mat)))
      c_use_mat <- c_use_mat[keep_c, , drop = FALSE]
      if (nrow(c_use_mat) > 6L) c_use_mat <- c_use_mat[seq_len(6L), , drop = FALSE]
      c_use_lbls <- rownames(c_use_mat)
    }
    graphics::points(c_use_mat[, 1], c_use_mat[, 2], pch = 3, col = "black", cex = 1.8, lwd = 2.4)
    if (!is.null(c_use_lbls) && length(c_use_lbls) == nrow(c_use_mat)) {
      cx0 <- c_use_mat[, 1]
      cy0 <- c_use_mat[, 2]
      rep_c <- repel_2d_labels(
        x0 = cx0,
        y0 = cy0 + 0.032 * diff(ylim),
        labels = c_use_lbls,
        xlim = xlim,
        ylim = ylim,
        cex = 0.76,
        n_iter = 25L,
        draw_segments = TRUE,
        seg_col = "grey45",
        seg_lwd = 0.75,
        anchor_x = cx0,
        anchor_y = cy0
      )
      for (k_c in seq_along(c_use_lbls)) {
        lbl_txt <- c_use_lbls[k_c]
        rx <- rep_c$x[k_c]
        ry <- rep_c$y[k_c]
        sw_lbl <- graphics::strwidth(lbl_txt, cex = 0.72, font = 2)
        sh_lbl <- graphics::strheight(lbl_txt, cex = 0.72, font = 2)
        graphics::rect(
          rx - sw_lbl * 0.56, ry - sh_lbl * 0.65, rx + sw_lbl * 0.56, ry + sh_lbl * 0.65,
          col = grDevices::adjustcolor("white", alpha.f = 0.84), border = NA
        )
        graphics::text(rx, ry, labels = lbl_txt, cex = 0.76, font = 4, col = "black")
      }
    }
  }

  if (!is.null(legend_args) && is.list(legend_args)) {
    if (!is.null(tier_ord) && tier_ord$is_multi_crowded) {
      fg_cols_leg <- vapply(fg_grps, function(g) {
        if (!is.null(group_cols) && g %in% names(group_cols)) group_cols[[g]] else tier_ord$palette[[g]]
      }, character(1))
      fg_pch_leg <- vapply(fg_grps, function(g) tier_ord$shapes[[g]], integer(1))
      pos_leg <- if (length(legend_args) >= 1L && is.character(legend_args[[1L]])) legend_args[[1L]] else "topright"
      graphics::legend(
        pos_leg,
        legend = c(packera_italic_expressions(fg_grps), "Background Congeners (Chull)"),
        col = c(fg_cols_leg, "#888888"),
        pch = c(fg_pch_leg, 1L),
        lty = c(rep(NA, length(fg_grps)), 2L),
        pt.bg = NA,
        cex = if (!is.null(legend_args$cex)) legend_args$cex else 0.70,
        bty = "n"
      )
    } else {
      do.call(graphics::legend, legend_args)
    }
  }
  invisible(list(xlim = xlim, ylim = ylim))
}

export_step03_fourier_pc_tables <- function(specimen_df, efa_df, output_csv_path) {
  tables_dir <- dirname(as.character(output_csv_path))
  if (!dir.exists(tables_dir)) {
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  }

  .write_pub_csv <- function(d_obj, tmp_p) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      data.table::fwrite(d_obj, file = tmp_p, na = "", quote = "auto", logical01 = FALSE)
    } else {
      utils::write.csv(d_obj, file = tmp_p, row.names = FALSE)
    }
  }

  # 1. fourier_pcs.csv
  fourier_pcs_path <- file.path(tables_dir, "fourier_pcs.csv")
  fourier_tmp <- paste0(fourier_pcs_path, ".", Sys.getpid(), ".tmp")
  .write_pub_csv(specimen_df, fourier_tmp)
  file.rename(fourier_tmp, fourier_pcs_path)

  # 2. symmetric_fourier_pcs.csv
  sym_cols <- intersect(
    c(
      "catalogNumber", "plant_individual_id", "taxon", "species_standardized", "scientificName",
      "determiner_tier", "is_type_specimen", "is_promoted_bronze", "authority_weight", "is_passive",
      "leaf_count", "basal_leaf_count", "cauline_leaf_count", "morph_tier_source",
      "leaf_quality_weight", "pct_asymmetry", "sym_var", "asym_var",
      "harmonic_roughness", "is_contour_outlier", "area_mm2",
      "petiole_length_mm", "petiole_width_mm", "blade_length_mm", "blade_width_mm",
      "petiole_to_blade_ratio", "is_dissected",
      "ECT_PC1", "ECT_PC2", "ECT_PC3", "ect_lobe_complexity",
      paste0("PC", 1:15)
    ),
    colnames(specimen_df)
  )
  sym_df <- if (length(sym_cols) > 0) specimen_df[, sym_cols, drop = FALSE] else specimen_df
  sym_pcs_path <- file.path(tables_dir, "symmetric_fourier_pcs.csv")
  sym_tmp <- paste0(sym_pcs_path, ".", Sys.getpid(), ".tmp")
  .write_pub_csv(sym_df, sym_tmp)
  file.rename(sym_tmp, sym_pcs_path)

  invisible(list(fourier_pcs = fourier_pcs_path, symmetric_fourier_pcs = sym_pcs_path))
}

format_cda_loading_labels <- function(raw_vars) {
  dict <- c(
    "PC1"                              = "PC1 [Basal Elongation]",
    "PC2"                              = "PC2 [Basal Base Shape]",
    "PC3"                              = "PC3 [Basal Dissection]",
    "PC4"                              = "PC4 [Basal Apex Taper]",
    "PC5"                              = "PC5 [Basal Lobe Depth]",
    "AR_inv"                           = "AR_inv [Involucre Ratio]",
    "AR_inv_std"                       = "AR_inv [Involucre Ratio]",
    "log_centsize"                     = "CS [Basal Leaf Size]",
    "log_centsize_std"                 = "CS [Basal Leaf Size]",
    "aspect_ratio"                     = "AR [Basal Blade L/W]",
    "aspect_ratio_std"                 = "AR [Basal Blade L/W]",

    "petiole_to_blade_ratio"           = "Petiole/Blade [Basal]",
    "petiole_to_blade_ratio_std"       = "Petiole/Blade [Basal]",
    "solidity"                         = "Solidity [Basal Dissection]",
    "solidity_std"                     = "Solidity [Basal Dissection]",
    "q25_solidity"                     = "Q25 Solidity [Basal Lobing]",
    "q25_solidity_std"                 = "Q25 Solidity [Basal Lobing]",
    "dissection_index"                 = "DI [Basal Margin Complexity]",
    "dissection_index_std"             = "DI [Basal Margin Complexity]",
    "within_rosette_efa_dispersion"    = "Rosette Disp. [Basal EFA]",
    "within_rosette_efa_dispersion_std"= "Rosette Disp. [Basal EFA]",
    "basal_vs_cauline_length_ratio"    = "Basal/Cauline Length",
    "basal_vs_cauline_length_ratio_std"= "Basal/Cauline Length",
    "basal_vs_cauline_solidity_delta"  = "Basal-Cauline Solidity",
    "basal_vs_cauline_solidity_delta_std" = "Basal-Cauline Solidity",
    "heteroblasty_efa_distance"        = "Heteroblasty EFA Dist.",
    "heteroblasty_efa_distance_std"    = "Heteroblasty EFA Dist.",
    "cauline_solidity_med"             = "Cauline Solidity [Stem]",
    "cauline_solidity_med_std"         = "Cauline Solidity [Stem]",
    "cauline_dissection_index_med"     = "Cauline DI [Stem]",
    "cauline_dissection_index_med_std" = "Cauline DI [Stem]",
    "Cauline_PC1"                      = "Cauline PC1 [Stem]",
    "Cauline_PC1_std"                  = "Cauline PC1 [Stem]",
    "Cauline_PC2"                      = "Cauline PC2 [Stem]",
    "Cauline_PC2_std"                  = "Cauline PC2 [Stem]",
    "Cauline_PC3"                      = "Cauline PC3 [Stem]",
    "Cauline_PC3_std"                  = "Cauline PC3 [Stem]",
    "ECT_PC1"                          = "ECT PC1 [Basal Topology]",
    "ECT_PC1_std"                      = "ECT PC1 [Basal Topology]",
    "ECT_PC2"                          = "ECT PC2 [Basal Topology]",
    "ECT_PC2_std"                      = "ECT PC2 [Basal Topology]",
    "ECT_PC3"                          = "ECT PC3 [Basal Topology]",
    "ECT_PC3_std"                      = "ECT PC3 [Basal Topology]",
    "ect_lobe_complexity"              = "ECT Lobe Complexity [Basal]",
    "ect_lobe_complexity_std"          = "ECT Lobe Complexity [Basal]",
    "Cauline_ECT_PC1"                  = "Cauline ECT PC1 [Stem]",
    "Cauline_ECT_PC1_std"              = "Cauline ECT PC1 [Stem]",
    "Cauline_ECT_PC2"                  = "Cauline ECT PC2 [Stem]",
    "Cauline_ECT_PC2_std"              = "Cauline ECT PC2 [Stem]",
    "cauline_ect_lobe_complexity"      = "Cauline ECT Lobing [Stem]",
    "cauline_ect_lobe_complexity_std"  = "Cauline ECT Lobing [Stem]",
    "heteroblasty_ect_distance"        = "Heteroblasty ECT Dist.",
    "heteroblasty_ect_distance_std"    = "Heteroblasty ECT Dist."
  )
  vapply(as.character(raw_vars), function(v) {
    if (v %in% names(dict)) unname(dict[v]) else v
  }, character(1), USE.NAMES = FALSE)
}

compute_ari_pure <- function(labels_true, labels_pred) {
  tab <- table(as.character(labels_true), as.character(labels_pred))
  n <- sum(tab)
  if (n <= 1) return(0.0)
  choose2 <- function(v) sum(as.numeric(v) * (as.numeric(v) - 1.0) / 2.0)
  sum_ij <- choose2(tab)
  sum_i  <- choose2(rowSums(tab))
  sum_j  <- choose2(colSums(tab))
  n_pairs <- n * (n - 1.0) / 2.0
  expected <- (sum_i * sum_j) / max(n_pairs, 1.0)
  max_idx  <- 0.5 * (sum_i + sum_j)
  denom    <- max_idx - expected
  if (abs(denom) < 1e-12) return(0.0)
  round((sum_ij - expected) / denom, 3)
}

# ==============================================================================
# Modular 300-DPI PNG Publication Figure Suite Delegation (Steps 03-06)
# Strictly 300-DPI PNG outputs (zero .pdf files); ggplot2 + patchwork + sf
# ==============================================================================
.resolve_scripts_figures_dir <- function() {
  if (dir.exists("scripts/figures")) return("scripts/figures")
  cand <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", "figures"), mustWork = FALSE),
    error = function(e) ""
  )
  if (nzchar(cand) && dir.exists(cand)) return(cand)
  "scripts/figures"
}

.fig_mod_dir <- .resolve_scripts_figures_dir()
for (.fmod in c("theme_packera_pub.R", "step03_figures.R", "step04_figures.R", "step05_figures.R", "step06_figures.R")) {
  .fp <- file.path(.fig_mod_dir, .fmod)
  if (file.exists(.fp)) source(.fp)
}

generate_step03_pub_figures <- function(
  efa_df,
  specimen_df,
  harm_mat = NULL,
  pca_fit = NULL,
  figures_dir = NULL,
  output_csv_path = NULL,
  sym_pca_fit = NULL
) {
  fdir <- resolve_active_figures_dir(figures_dir, output_csv_path)
  render_step03_figures_suite(
    efa_df = efa_df,
    specimen_df = specimen_df,
    harm_mat = harm_mat,
    pca_fit = pca_fit,
    figures_dir = fdir,
    output_csv_path = output_csv_path,
    sym_pca_fit = sym_pca_fit
  )
}

render_allometry_diagnostic_pdf <- function(sub_df, r2_size, p_size, r2_int, output_pdf, p_int = 1.0) {
  out_png <- if (!is.null(output_pdf) && nzchar(output_pdf)) normalize_png_path(output_pdf) else "outputs/figures/fig04a_foliar_allometry_audit.png"
  render_step04_allometry_png(
    sub_df = sub_df,
    r2_size = r2_size,
    p_size = p_size,
    r2_int = r2_int,
    output_path = out_png,
    p_int = p_int
  )
}

render_cda_biplot_4panel_pdf <- function(audit_df, cda_res, gmm_res, output_pdf, feature_cols = NULL, manova_res = NULL) {
  out_png <- if (!is.null(output_pdf) && nzchar(output_pdf)) normalize_png_path(output_pdf) else "outputs/figures/fig04b_gmm_cda_morphometrics.png"
  render_step04_gmm_cda_png(
    audit_df = audit_df,
    cda_res = cda_res,
    gmm_res = gmm_res,
    output_path = out_png,
    feature_cols = feature_cols,
    manova_res = manova_res
  )
}

generate_step04_pub_figures <- function(
  audit_df,
  cda_res,
  gmm_res,
  figures_dir = NULL,
  reference_path = NULL,
  eco_cda_res = NULL,
  manova_res = NULL,
  batch_csv = NULL
) {
  fdir <- resolve_active_figures_dir(figures_dir, reference_path)
  canon_04b <- file.path(fdir, "fig04b_gmm_cda_morphometrics.png")
  if (!file.exists(canon_04b)) {
    render_step04_gmm_cda_png(
      audit_df = audit_df,
      cda_res = cda_res,
      gmm_res = gmm_res,
      output_path = canon_04b,
      feature_cols = if (!is.null(cda_res)) cda_res$selected_features else NULL,
      manova_res = manova_res
    )
  }
  if (exists("render_fig04c_supervised_ecoregion_cda", mode = "function")) {
    render_fig04c_supervised_ecoregion_cda(
      audit_df = audit_df,
      eco_cda_res = eco_cda_res,
      gmm_res = gmm_res,
      figures_dir = fdir,
      output_plot_path = reference_path,
      manova_res = manova_res,
      batch_csv = batch_csv
    )
  }
  invisible(TRUE)
}

.pub_figures_dir <- if (dir.exists("scripts/morphometrics")) {
  "scripts/morphometrics"
} else {
  tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
    error = function(e) "scripts/morphometrics"
  )
}
source(file.path(.pub_figures_dir, "pub_ecoregion_manova_utils.R"))
source(file.path(.pub_figures_dir, "pub_synthesis_plates.R"))
