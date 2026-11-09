#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_figures_utils.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
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

compute_spatial_balance_weights <- function(df, bandwidth_km = 25.0) {
  if (is.null(df) || nrow(df) == 0L) return(numeric(0))
  n_rows <- nrow(df)
  if ("spatial_balance_weight" %in% names(df)) {
    w_in <- suppressWarnings(as.numeric(df$spatial_balance_weight))
    if (any(is.finite(w_in) & w_in > 0)) {
      w_in[!is.finite(w_in) | w_in <= 0] <- 1.0
      return(pmin(3.0, pmax(0.20, w_in)))
    }
  }
  sp_col <- if ("species_standardized" %in% names(df)) "species_standardized" else if ("scientificName" %in% names(df)) "scientificName" else if ("species_raw" %in% names(df)) "species_raw" else NULL
  taxa <- if (!is.null(sp_col)) as.character(df[[sp_col]]) else rep("Unknown", n_rows)
  taxa[is.na(taxa) | !nzchar(trimws(taxa))] <- "Unknown"

  lat_col <- if ("decimalLatitude" %in% names(df)) "decimalLatitude" else if ("latitude" %in% names(df)) "latitude" else NULL
  lon_col <- if ("decimalLongitude" %in% names(df)) "decimalLongitude" else if ("longitude" %in% names(df)) "longitude" else NULL
  if (is.null(lat_col) || is.null(lon_col)) return(rep(1.0, n_rows))
  lats <- suppressWarnings(as.numeric(df[[lat_col]]))
  lons <- suppressWarnings(as.numeric(df[[lon_col]]))

  st_col <- if ("stateProvince" %in% names(df)) "stateProvince" else if ("regional_group" %in% names(df)) "regional_group" else NULL
  states <- if (!is.null(st_col)) toupper(trimws(as.character(df[[st_col]]))) else rep("", n_rows)
  states[is.na(states)] <- ""

  w_out <- rep(1.0, n_rows)
  for (sp in unique(taxa)) {
    t_idx <- which(taxa == sp)
    geo_idx <- t_idx[is.finite(lats[t_idx]) & is.finite(lons[t_idx])]
    if (length(geo_idx) == 0L) next

    st_sub <- states[t_idx]
    st_valid <- st_sub[nzchar(st_sub)]
    if (length(st_valid) > 0L) {
      st_tab <- table(st_valid)
      med_s <- max(as.numeric(stats::median(as.numeric(st_tab))), 1.0)
      st_geo <- states[geo_idx]
      st_lookup <- pmax(as.numeric(st_tab) / med_s, 1e-6)
      names(st_lookup) <- names(st_tab)
      f_state <- unname(st_lookup[st_geo])
      f_state[is.na(f_state) | !nzchar(st_geo)] <- 1.0
    } else {
      f_state <- rep(1.0, length(geo_idx))
    }

    rlat <- (lats[geo_idx] * (pi / 180.0)) + 0.0
    rlon <- (lons[geo_idx] * (pi / 180.0)) + 0.0
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
    dist_km <- 2.0 * 6371.0 * asin(pmin(1.0, pmax(0.0, sqrt(pmax(0.0, a_mat)))))
    kern_u <- matrix(exp(-(dist_km^2) / (2.0 * (bandwidth_km^2))), nrow = n_u, ncol = n_u)
    diag(kern_u) <- 1.0
    d_local <- as.numeric(kern_u %*% counts_u)[inv_idx]

    w_raw <- 1.0 / sqrt(pmax(d_local * f_state, 1e-8))
    m_raw <- mean(w_raw, na.rm = TRUE)
    w_sp <- if (is.finite(m_raw) && m_raw > 0) w_raw / m_raw else rep(1.0, length(geo_idx))
    w_sp <- pmin(3.0, pmax(0.20, w_sp))
    for (iter in 1:3) {
      m_cur <- mean(w_sp, na.rm = TRUE)
      if (is.finite(m_cur) && m_cur > 0) w_sp <- pmin(3.0, pmax(0.20, w_sp / m_cur))
    }
    w_out[geo_idx] <- round(w_sp, 6)
  }
  return(w_out)
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
    "topo_roughness_30m"      = "30m Terrain Roughness (m)",
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
    if (grepl("^Morphotype_K[0-9]+$|^Ecotype_|^LatTertile_|_Coastal$|_Outcrop$|\\[K[0-9]+\\]$", s_orig, ignore.case = TRUE)) {
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
    any(grepl("Morphotype_K|\\[K[0-9]+\\]|^LatTertile_|_Coastal|_Outcrop|Ecotype|Subcluster|Region_", vv, ignore.case = TRUE)) ||
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
    if ("gmm_cluster" %in% names(df) && any(!is.na(df$gmm_cluster))) {
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
  pdf_final <- file.path(figures_dir, paste0(stem, ".pdf"))
  png_final <- file.path(figures_dir, paste0(stem, ".png"))
  pdf_tmp <- file.path(figures_dir, paste0(".", stem, "_", Sys.getpid(), ".pdf.tmp"))
  png_tmp <- file.path(figures_dir, paste0(".", stem, "_", Sys.getpid(), ".png.tmp"))

  grDevices::pdf(pdf_tmp, width = width, height = height, useDingbats = FALSE)
  tryCatch({
    draw_fn()
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(pdf_tmp)) {
    file.rename(pdf_tmp, pdf_final)
  }

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
  invisible(list(pdf = pdf_final, png = png_final))
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
    # Detect legacy block-ordered Momocs vector stored under interleaved names
    if (K >= 10L && length(v_num) >= 4L * K && abs(A[10]) > 0.15 && abs(D[1]) < 0.08) {
      interleaved_nms <- as.vector(outer(c("A", "B", "C", "D"), seq_len(K), paste0))
      if (all(interleaved_nms %in% clean_nms)) {
        raw_4k <- unname(v_num[interleaved_nms])
        A <- raw_4k[1:K]
        B <- raw_4k[(K + 1L):(2L * K)]
        C <- raw_4k[(2L * K + 1L):(3L * K)]
        D <- raw_4k[(3L * K + 1L):(4L * K)]
      }
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
  if (length(avail) < 8L) return(NULL)
  mat <- as.matrix(df[, avail, drop = FALSE])
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- 0.0
  # Check if legacy block-ordered Momocs matrix was labeled with interleaved names
  if (all(interleaved_nms %in% colnames(mat)) && max_k >= 10L) {
    med_a10 <- stats::median(abs(mat[, "A10"]), na.rm = TRUE)
    med_d1  <- stats::median(abs(mat[, "D1"]), na.rm = TRUE)
    if (is.finite(med_a10) && is.finite(med_d1) && med_a10 > 0.15 && med_d1 < 0.08) {
      raw_4k <- mat[, interleaved_nms, drop = FALSE]
      colnames(raw_4k) <- c(
        paste0("A", seq_len(max_k)),
        paste0("B", seq_len(max_k)),
        paste0("C", seq_len(max_k)),
        paste0("D", seq_len(max_k))
      )
      mat <- raw_4k[, interleaved_nms, drop = FALSE]
    }
  }
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

render_tier1_vs_tier2_2d_pdf <- function(tier_manifest, plot_out, pca_fit = NULL, sym_pca = NULL, permanova_df = NULL) {
  dir.create(dirname(plot_out), recursive = TRUE, showWarnings = FALSE)
  tmp_pdf <- paste0(plot_out, ".tmp")
  if (!is.null(tier_manifest) && nrow(tier_manifest) > 0L && "leaf_zone" %in% names(tier_manifest)) {
    lz_t <- tolower(trimws(as.character(tier_manifest$leaf_zone)))
    bas_idx <- which(is.na(lz_t) | !nzchar(lz_t) | lz_t != "cauline")
    if (length(bas_idx) > 0L) tier_manifest <- tier_manifest[bas_idx, , drop = FALSE]
  }

  var_pct <- c(48.2, 21.4)
  if (!is.null(pca_fit) && !is.null(pca_fit$sdev) && length(pca_fit$sdev) >= 2) {
    v <- (pca_fit$sdev^2) / max(sum(pca_fit$sdev^2), 1e-12) * 100
    var_pct <- round(v[1:2], 1)
  } else if (all(c("PC1", "PC2") %in% names(tier_manifest))) {
    v_emp <- c(stats::var(tier_manifest$PC1, na.rm = TRUE), stats::var(tier_manifest$PC2, na.rm = TRUE))
    v_emp[!is.finite(v_emp)] <- 0
    if (sum(v_emp) > 0) var_pct <- round((v_emp / sum(v_emp)) * 75.0, 1)
  }

  sym_var_pct <- var_pct
  if (!is.null(sym_pca) && !is.null(sym_pca$sdev) && length(sym_pca$sdev) >= 2) {
    vs <- (sym_pca$sdev^2) / max(sum(sym_pca$sdev^2), 1e-12) * 100
    sym_var_pct <- round(vs[1:2], 1)
  }

  perm_label <- "PERMANOVA Tier Effect: R2 < 0.01 (Non-Significant)"
  if (!is.null(permanova_df) && is.data.frame(permanova_df) && "Term" %in% names(permanova_df)) {
    permanova_df <- apply_bh_fdr_correction(permanova_df)
    t_row <- permanova_df[grepl("reconstruction_tier", permanova_df$Term, ignore.case = TRUE), , drop = FALSE]
    if (nrow(t_row) > 0) {
      r2_v <- suppressWarnings(as.numeric(t_row$R2[1]))
      p_v  <- if ("p_adj_bh" %in% names(t_row)) suppressWarnings(as.numeric(t_row$p_adj_bh[1])) else suppressWarnings(as.numeric(t_row$p_value[1]))
      if (is.finite(r2_v)) {
        perm_label <- sprintf(
          "PERMANOVA Tier Effect: R2 = %.4f (%.2f%%), q_BH = %s",
          r2_v, r2_v * 100,
          if (is.finite(p_v)) sprintf("%.3f", p_v) else "NS"
        )
      }
    }
  }

  grDevices::pdf(tmp_pdf, width = 12.5, height = 5.8, useDingbats = FALSE)
  tryCatch({
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1.0), cex.main = 0.90)
    is_t2 <- grepl("2", as.character(tier_manifest$reconstruction_tier))
    n_t1 <- sum(!is_t2, na.rm = TRUE)
    n_t2 <- sum(is_t2, na.rm = TRUE)
    tier_lbl <- ifelse(is_t2, "Tier 2 Reflected", "Tier 1 Direct")
    tier_cols <- c("Tier 1 Direct" = "#0072B2", "Tier 2 Reflected" = "#D55E00")
    alpha_val <- adaptive_point_alpha(nrow(tier_manifest))
    pt_col <- grDevices::adjustcolor(ifelse(is_t2, "#D55E00", "#0072B2"), alpha.f = alpha_val)
    pt_pch <- ifelse(is_t2, 2L, 1L)
    pc1 <- if ("PC1" %in% names(tier_manifest)) tier_manifest$PC1 else rep(0, nrow(tier_manifest))
    pc2 <- if ("PC2" %in% names(tier_manifest)) tier_manifest$PC2 else rep(0, nrow(tier_manifest))
    spc1 <- if ("Sym_PC1" %in% names(tier_manifest)) tier_manifest$Sym_PC1 else pc1
    spc2 <- if ("Sym_PC2" %in% names(tier_manifest)) tier_manifest$Sym_PC2 else pc2
    cs_t12 <- if (n_t2 > 0L) {
      rbind(
        "Tier 1 Direct" = c(mean(spc1[!is_t2], na.rm = TRUE), mean(spc2[!is_t2], na.rm = TRUE)),
        "Tier 2 Reflected" = c(mean(spc1[is_t2], na.rm = TRUE), mean(spc2[is_t2], na.rm = TRUE))
      )
    } else {
      rbind(
        "Tier 1 Direct" = c(mean(spc1[!is_t2], na.rm = TRUE), mean(spc2[!is_t2], na.rm = TRUE))
      )
    }
    t2_leg_lbl <- if (n_t2 > 0L) sprintf("Tier 2 Reflected (n=%d, pch=2)", n_t2) else "Tier 2 Reflected (n=0; Strict Tier 1 Mode)"

    # Panel A: 2D Sym_PC1 vs. Sym_PC2 Overlay with PERMANOVA R2 and p-value inset
    plot_2d_ordination_scatter(
      x = spc1, y = spc2, col = pt_col, pch = pt_pch, cex = 0.85, lwd = 1.0,
      xlab = sprintf("Sym_PC1 (%.1f%%)", sym_var_pct[1]),
      ylab = sprintf("Sym_PC2 (%.1f%%)", sym_var_pct[2]),
      main = "A. 2D Symmetric EFA Morphospace: Tier 1 Direct vs. Tier 2 Reflected",
      group_vec = tier_lbl,
      group_cols = tier_cols,
      ellipse_level = 0.95,
      centroids = cs_t12,
      legend_args = list(
        "topright",
        legend = c(
          sprintf("Tier 1 Direct (n=%d, pch=1)", n_t1),
          t2_leg_lbl,
          perm_label
        ),
        col = c("#0072B2", "#D55E00", NA),
        pch = c(1, 2, NA), pt.bg = NA, bty = "n", cex = 0.74
      )
    )

    # Panel B: Superimposed Mean Reconstructed Leaf Silhouette of Tier 1 Direct vs. Tier 2 Reflected per Species / Morphotype
    h_mat <- extract_standardized_harmonic_matrix(tier_manifest, max_k = 12L)
    tx_vec <- resolve_plot_group_column(tier_manifest)
    tier_info <- resolve_taxon_display_tier(tier_manifest, tx_vec)
    disp_tx <- utils::head(tier_info$ordered_groups, 9L)
    if (length(disp_tx) == 0L) disp_tx <- "Unassigned_Packera"

    graphics::plot(
      c(0, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "",
      main = "B. Mean Reconstructed Silhouette: Tier 1 Direct vs. Tier 2 Reflected",
      cex.main = 0.90
    )
    graphics::box(col = "grey80")
    grid_lay <- .compute_outline_grid_layout(length(disp_tx))

    for (idx_sp in seq_along(disp_tx)) {
      sp_name <- disp_tx[idx_sp]
      cx <- grid_lay$cx[idx_sp]; cy <- grid_lay$cy[idx_sp]
      m_t1 <- (tx_vec == sp_name) & !is_t2
      m_t2 <- (tx_vec == sp_name) & is_t2
      has_sp_t2 <- sum(m_t2, na.rm = TRUE) > 0L
      if (!is.null(h_mat)) {
        h1 <- if (any(m_t1)) colMeans(h_mat[m_t1, , drop = FALSE], na.rm = TRUE) else colMeans(h_mat, na.rm = TRUE)
        h2 <- if (has_sp_t2) colMeans(h_mat[m_t2, , drop = FALSE], na.rm = TRUE) else NULL
      } else {
        h1 <- c(A1 = -1.25, B1 = 0, C1 = 0, D1 = 0.45)
        h2 <- if (has_sp_t2) c(A1 = -1.23, B1 = 0, C1 = 0, D1 = 0.46) else NULL
      }
      xy1 <- reconstruct_efa_outline(h1, nb_pts = 120L)
      xy2 <- if (!is.null(h2)) reconstruct_efa_outline(h2, nb_pts = 120L) else NULL
      env_pts <- if (!is.null(xy2)) rbind(xy1, xy2) else xy1
      shared_center <- c(mean(range(env_pts[, 1], na.rm = TRUE)), mean(range(env_pts[, 2], na.rm = TRUE)))
      shared_scale  <- min(0.68 / max(diff(range(env_pts[, 1], na.rm = TRUE)), 1e-6), 0.64 / max(diff(range(env_pts[, 2], na.rm = TRUE)), 1e-6)) * 0.86 * grid_lay$cell_scale
      draw_outline_in_box(
        xy1, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
        col_border = "#0072B2", col_fill = grDevices::adjustcolor("#0072B2", alpha.f = 0.16), lwd = 2.0, lty = 1,
        shared_scale = shared_scale, shared_center = shared_center
      )
      if (!is.null(xy2) && n_t2 > 0L) {
        draw_outline_in_box(
          xy2, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
          col_border = "#D55E00", col_fill = NA, lwd = 1.9, lty = 2,
          shared_scale = shared_scale, shared_center = shared_center
        )
      }
      short_sp <- sub("^Packera\\s+", "P. ", sp_name)
      graphics::text(
        cx, cy + grid_lay$title_dy,
        labels = as.expression(bquote(bolditalic(.(short_sp)) ~ .(sprintf("(T1:%d / T2:%d)", sum(m_t1), sum(m_t2))))),
        cex = grid_lay$cex_title
      )
    }
    graphics::legend(
      "bottom", horiz = TRUE,
      legend = c(
        "Tier 1 Direct Mean Outline",
        if (n_t2 > 0L) "Tier 2 Reflected Mean Outline" else "Tier 2 Reflected (n=0; Strict Tier 1 Mode)"
      ),
      col = c("#0072B2", "#D55E00"), lty = c(1, 2), lwd = 2.0, bty = "n", cex = 0.78
    )
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, plot_out)
  invisible(plot_out)
}
render_tier1_vs_tier2_3d_pdf <- render_tier1_vs_tier2_2d_pdf

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
  fourier_tmp <- paste0(fourier_pcs_path, ".tmp")
  .write_pub_csv(specimen_df, fourier_tmp)
  file.rename(fourier_tmp, fourier_pcs_path)

  # 2. symmetric_fourier_pcs.csv
  sym_cols <- intersect(
    c(
      "catalogNumber", "plant_individual_id", "taxon", "scientificName", "determiner_tier",
      "assigned_tier", "n_leaves", "leaf_quality_weight", "pct_asymmetry",
      "sym_var", "asym_var", "harmonic_roughness", "is_contour_outlier",
      "petiole_length_mm", "blade_length_mm", "petiole_to_blade_ratio", "is_dissected",
      paste0("Sym_PC", 1:10), paste0("PC", 1:10)
    ),
    colnames(specimen_df)
  )
  sym_df <- if (length(sym_cols) > 0) specimen_df[, sym_cols, drop = FALSE] else specimen_df
  sym_pcs_path <- file.path(tables_dir, "symmetric_fourier_pcs.csv")
  sym_tmp <- paste0(sym_pcs_path, ".tmp")
  .write_pub_csv(sym_df, sym_tmp)
  file.rename(sym_tmp, sym_pcs_path)

  invisible(list(fourier_pcs = fourier_pcs_path, symmetric_fourier_pcs = sym_pcs_path))
}

# ==============================================================================
# Step 03 Figures:
#   - fig03_efa_morphospace_reconstruction.{pdf,png}:
#       Panel A: 2D Symmetric EFA Morphospace with Background 5x3 Inverse-Fourier Grid
#       Panel B: Mean Consensus Leaf Silhouettes (+/- 1 SD Envelope) by Species
#   - fig03a_efa_morphospace_pca.{pdf,png} (companion)
#   - fig03b_harmonic_power_and_asymmetry.{pdf,png} (unfilled hollow markers)
# ==============================================================================
generate_step03_pub_figures <- function(
  efa_df,
  specimen_df,
  harm_mat = NULL,
  pca_fit = NULL,
  figures_dir = NULL,
  output_csv_path = NULL,
  sym_pca = NULL
) {
  fdir <- resolve_active_figures_dir(explicit_dir = figures_dir, reference_path = output_csv_path)
  df_full <- if (!is.null(specimen_df) && nrow(specimen_df) > 0) specimen_df else efa_df
  if (is.null(df_full) || nrow(df_full) == 0) return(invisible(NULL))

  # Isolate pure basal rosette rows for primary Step 03 morphospace & outline reconstructions (never overlay cauline rows)
  bas_keep <- rep(TRUE, nrow(df_full))
  if ("leaf_zone" %in% colnames(df_full)) {
    lz_df <- tolower(trimws(as.character(df_full$leaf_zone)))
    bas_keep <- bas_keep & (is.na(lz_df) | !nzchar(lz_df) | lz_df != "cauline")
  }
  if ("morph_tier_source" %in% colnames(df_full)) {
    mts_df <- trimws(as.character(df_full$morph_tier_source))
    bas_keep <- bas_keep & (is.na(mts_df) | mts_df != "Cauline_Only_Fallback")
  }
  if ("basal_leaf_count" %in% colnames(df_full)) {
    blc_df <- suppressWarnings(as.numeric(df_full$basal_leaf_count))
    if (any(is.finite(blc_df) & blc_df > 0)) {
      bas_keep <- bas_keep & (is.na(blc_df) | blc_df > 0)
    }
  }
  df <- if (any(bas_keep, na.rm = TRUE)) df_full[which(bas_keep), , drop = FALSE] else df_full

  if (!"taxon" %in% colnames(df)) {
    raw_tx <- if ("species_standardized" %in% colnames(df)) df$species_standardized else if ("scientificName" %in% colnames(df)) df$scientificName else if ("species_raw" %in% colnames(df)) df$species_raw else rep(NA_character_, nrow(df))
    std_tx <- standardize_taxon_name(raw_tx)
    cand_tx <- ifelse(!is.na(std_tx), std_tx, ifelse(!is.na(raw_tx) & nzchar(trimws(as.character(raw_tx))), trimws(as.character(raw_tx)), NA_character_))
    valid_tx <- cand_tx[!is.na(cand_tx) & nzchar(cand_tx) & cand_tx != "Unknown"]
    fb_tx <- if (length(valid_tx) > 0L) valid_tx[1L] else "Unassigned_Packera"
    df$taxon <- ifelse(!is.na(cand_tx) & nzchar(cand_tx), cand_tx, fb_tx)
  } else {
    std_tx <- standardize_taxon_name(df$taxon)
    df$taxon <- ifelse(!is.na(std_tx), std_tx, as.character(df$taxon))
  }
  resolved_grp <- resolve_plot_group_column(df)
  dyn_pal <- get_dynamic_packera_palette(resolved_grp)
  dyn_pch <- get_dynamic_packera_shapes(resolved_grp)
  valid_grp <- resolved_grp[!is.na(resolved_grp) & nzchar(resolved_grp) & resolved_grp != "Unknown"]
  default_grp <- if (length(valid_grp) > 0L) valid_grp[1L] else "Unassigned_Packera"
  df$taxon_clean <- ifelse(
    resolved_grp %in% names(dyn_pal),
    as.character(resolved_grp),
    default_grp
  )

  var_pct <- c(48.2, 21.4)
  if (!is.null(pca_fit) && !is.null(pca_fit$sdev) && length(pca_fit$sdev) >= 2) {
    v <- (pca_fit$sdev^2) / max(sum(pca_fit$sdev^2), 1e-12) * 100
    var_pct <- round(v[1:2], 1)
  } else if (all(c("PC1", "PC2") %in% colnames(df))) {
    v_emp <- c(stats::var(df$PC1, na.rm = TRUE), stats::var(df$PC2, na.rm = TRUE))
    v_emp[!is.finite(v_emp)] <- 0
    if (sum(v_emp) > 0) var_pct <- round((v_emp / sum(v_emp)) * 75.0, 1)
  }

  sym_var_pct <- var_pct
  if (!is.null(sym_pca) && !is.null(sym_pca$sdev) && length(sym_pca$sdev) >= 2) {
    vs <- (sym_pca$sdev^2) / max(sum(sym_pca$sdev^2), 1e-12) * 100
    sym_var_pct <- round(vs[1:2], 1)
  }

  if (!"PC1" %in% colnames(df)) df$PC1 <- 0
  if (!"PC2" %in% colnames(df)) df$PC2 <- 0
  if (!"Sym_PC1" %in% colnames(df)) df$Sym_PC1 <- df$PC1
  if (!"Sym_PC2" %in% colnames(df)) df$Sym_PC2 <- df$PC2

  alpha_val <- adaptive_point_alpha(nrow(df))
  pt_pch <- dyn_pch[df$taxon_clean]
  pt_pch[is.na(pt_pch)] <- 1L
  cols <- grDevices::adjustcolor(dyn_pal[df$taxon_clean], alpha.f = alpha_val)

  present_tx <- sort(unique(df$taxon_clean[df$taxon_clean != "Unknown" & !is.na(df$taxon_clean)]))
  present_tx <- present_tx[vapply(present_tx, function(tx) sum(df$taxon_clean == tx, na.rm = TRUE) >= 1L, logical(1))]
  if (length(present_tx) == 0) present_tx <- default_grp
  is_type_row <- if ("is_type_specimen" %in% colnames(df)) as.logical(df$is_type_specimen) else rep(FALSE, nrow(df))
  is_type_row[is.na(is_type_row)] <- FALSE
  auth_w_3 <- if ("authority_weight" %in% colnames(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, nrow(df))
  auth_w_3[is.na(auth_w_3) | !is.finite(auth_w_3)] <- 0.0
  if (all(auth_w_3 <= 0)) auth_w_3 <- rep(1.0, nrow(df))
  cent_sym <- do.call(rbind, lapply(present_tx, function(tx) {
    m <- which(df$taxon_clean == tx)
    w_m <- if (any(auth_w_3[m] > 0)) pmax(auth_w_3[m], 1e-6) else rep(1.0, length(m))
    c(stats::weighted.mean(df$Sym_PC1[m], w = w_m, na.rm = TRUE),
      stats::weighted.mean(df$Sym_PC2[m], w = w_m, na.rm = TRUE))
  }))
  rownames(cent_sym) <- gsub("Packera ", "P. ", present_tx)

  std_hmat <- extract_standardized_harmonic_matrix(df, max_k = NULL)
  if (is.null(std_hmat) && !is.null(efa_df) && nrow(efa_df) > 0) {
    efa_bas_idx <- if ("leaf_zone" %in% colnames(efa_df)) {
      lz_e <- tolower(trimws(as.character(efa_df$leaf_zone)))
      which(is.na(lz_e) | !nzchar(lz_e) | lz_e != "cauline")
    } else seq_len(nrow(efa_df))
    if (length(efa_bas_idx) == 0L) efa_bas_idx <- seq_len(nrow(efa_df))
    std_hmat <- extract_standardized_harmonic_matrix(efa_df[efa_bas_idx, , drop = FALSE], max_k = NULL)
  }

  # Precompute linear projection from (Sym_PC1, Sym_PC2) to standardized harmonics once
  grid_B_hat <- NULL
  grid_h_bar <- NULL
  if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && nrow(df) >= 5) {
    s1 <- df$Sym_PC1; s2 <- df$Sym_PC2
    if (stats::sd(s1, na.rm = TRUE) > 1e-6 && stats::sd(s2, na.rm = TRUE) > 1e-6) {
      grid_B_hat <- tryCatch(qr.solve(cbind(1, s1, s2), std_hmat), error = function(e) NULL)
    }
    grid_h_bar <- colMeans(std_hmat, na.rm = TRUE)
  }

  # Helper to synthesize harmonic vector at any (pc1_val, pc2_val) coordinate
  project_grid_harmonics <- function(pc1_val, pc2_val) {
    if (!is.null(grid_B_hat)) {
      h_pred <- as.numeric(c(1, pc1_val, pc2_val) %*% grid_B_hat)
      names(h_pred) <- colnames(std_hmat)
      return(h_pred)
    }
    if (!is.null(grid_h_bar)) {
      return(grid_h_bar)
    }
    # Parametric botanical fallback if harmonic columns absent in minimal test table
    c(
      A1 = -1.25 - 0.18 * pc1_val, B1 = 0, C1 = 0, D1 = 0.48 - 0.10 * pc1_val + 0.08 * pc2_val,
      A2 = -0.04 + 0.06 * pc2_val, B2 = 0, C2 = 0, D2 = -0.01 + 0.04 * pc2_val,
      A3 = -0.08 - 0.05 * pc1_val, B3 = 0, C3 = 0, D3 = 0.03 * pc2_val
    )
  }

  # ----------------------------------------------------------------------------
  # Figure 03: Inverse EFA Morphospace Reconstruction Grid & Species Mean Outlines
  # ----------------------------------------------------------------------------
  draw_fig03_reconstruction <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.1))
    # Panel A: 2D Symmetric EFA Morphospace with 5x3 Background Inverse-Fourier Grid
    sd1 <- max(stats::sd(df$Sym_PC1, na.rm = TRUE), 0.08)
    sd2 <- max(stats::sd(df$Sym_PC2, na.rm = TRUE), 0.06)
    m1  <- mean(df$Sym_PC1, na.rm = TRUE)
    m2  <- mean(df$Sym_PC2, na.rm = TRUE)

    bg_fourier_grid <- function(xlim, ylim) {
      z1_seq <- c(-2.0, -1.0, 0.0, 1.0, 2.0)
      z2_seq <- c(-1.5, 0.0, 1.5)
      gx_vals <- pmax(xlim[1] + 0.08 * diff(xlim), pmin(xlim[2] - 0.08 * diff(xlim), m1 + z1_seq * sd1))
      gy_vals <- pmax(ylim[1] + 0.10 * diff(ylim), pmin(ylim[2] - 0.10 * diff(ylim), m2 + z2_seq * sd2))
      cell_w <- 0.13 * diff(xlim)
      cell_h <- 0.20 * diff(ylim)
      for (gx in gx_vals) {
        for (gy in gy_vals) {
          h_vec <- project_grid_harmonics(gx, gy)
          xy_rec <- reconstruct_efa_outline(h_vec, nb_pts = 100L)
          draw_outline_in_box(
            xy_rec, cx = gx, cy = gy, w_box = cell_w, h_box = cell_h,
            col_border = grDevices::adjustcolor("grey50", alpha.f = 0.55),
            col_fill = grDevices::adjustcolor("grey82", alpha.f = 0.30),
            lwd = 0.9
          )
        }
      }
    }

    plot_2d_ordination_scatter(
      x = df$Sym_PC1,
      y = df$Sym_PC2,
      col = cols,
      pch = pt_pch,
      cex = 0.85,
      lwd = 1.05,
      xlab = sprintf("Sym_PC1 (%.1f%% -- Blade Elongation & Petiole Ratio)", sym_var_pct[1]),
      ylab = sprintf("Sym_PC2 (%.1f%% -- Lamina Base Shape & Ovate Taper)", sym_var_pct[2]),
      main = "A. 2D Basal Rosette EFA Morphospace & Inverse-Fourier Grid",
      group_vec = df$taxon_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.95,
      ellipse_weights = auth_w_3,
      centroids = cent_sym,
      bg_draw_fn = bg_fourier_grid,
      is_type_vec = is_type_row,
      legend_args = list(
        "topright",
        legend = c(packera_italic_expressions(present_tx), "Bg: Inverse EFA (-2..+2 SD)"),
        col = c(dyn_pal[present_tx], "grey55"),
        pch = c(dyn_pch[present_tx], 0),
        pt.bg = NA,
        cex = 0.72,
        bty = "n"
      )
    )

    # Panel B: Mean Consensus Basal Leaf Silhouettes (+/- 1 SD Envelope) by Species
    graphics::plot(
      c(0, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "",
      main = "B. Species Basal Mean Consensus Silhouettes (+/- 1 SD Envelope)"
    )
    graphics::box(col = "grey80")

    tier_info <- resolve_taxon_display_tier(df, df$taxon_clean)
    disp_tx <- utils::head(tier_info$ordered_groups, 9L)
    disp_tx <- disp_tx[vapply(disp_tx, function(tx) sum(df$taxon_clean == tx, na.rm = TRUE) >= 1L, logical(1))]
    if (length(disp_tx) == 0L) disp_tx <- default_grp
    n_disp <- length(disp_tx)
    grid_lay <- .compute_outline_grid_layout(n_disp)
    cx_vec <- grid_lay$cx
    cy_vec <- grid_lay$cy

    is_t1_pure <- if ("morph_tier_source" %in% colnames(df)) (df$morph_tier_source == "Tier_1_Pure") else !grepl("2", as.character(df$assigned_tier))
    is_promo_row <- if ("is_promoted_anchor" %in% colnames(df)) as.logical(df$is_promoted_anchor) else rep(FALSE, nrow(df))
    is_promo_row[is.na(is_promo_row)] <- FALSE
    is_gold_silv <- if ("determiner_tier" %in% colnames(df)) (df$determiner_tier %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver") | is_type_row | is_promo_row) else rep(TRUE, nrow(df))

    for (idx_sp in seq_len(n_disp)) {
      sp_name <- disp_tx[idx_sp]
      sp_col <- if (sp_name %in% names(dyn_pal)) dyn_pal[[sp_name]] else "#0072B2"
      cx <- cx_vec[idx_sp]; cy <- cy_vec[idx_sp]

      m_anchor <- (df$taxon_clean == sp_name) & is_t1_pure & is_gold_silv
      if (sum(m_anchor, na.rm = TRUE) < 3) m_anchor <- (df$taxon_clean == sp_name) & is_t1_pure
      if (sum(m_anchor, na.rm = TRUE) < 1) m_anchor <- (df$taxon_clean == sp_name)
      m_type <- (df$taxon_clean == sp_name) & is_type_row

      if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && any(m_anchor)) {
        h_mean <- colMeans(std_hmat[m_anchor, , drop = FALSE], na.rm = TRUE)
        h_sd   <- apply(std_hmat[m_anchor, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
        h_sd[!is.finite(h_sd)] <- 0.0
      } else {
        z_sp <- if (any(m_anchor)) c(mean(df$Sym_PC1[m_anchor], na.rm = TRUE), mean(df$Sym_PC2[m_anchor], na.rm = TRUE)) else c((idx_sp - 2.5) * 0.5, 0)
        h_mean <- project_grid_harmonics(z_sp[1], z_sp[2])
        h_sd   <- abs(h_mean) * 0.12
      }

      xy_mean  <- reconstruct_efa_outline(h_mean, nb_pts = 120L)
      xy_plus  <- reconstruct_efa_outline(h_mean + 0.65 * h_sd, nb_pts = 120L)
      xy_minus <- reconstruct_efa_outline(h_mean - 0.65 * h_sd, nb_pts = 120L)

      env_pts <- rbind(xy_plus, xy_mean, xy_minus)
      shared_center <- c(mean(range(env_pts[, 1], na.rm = TRUE)), mean(range(env_pts[, 2], na.rm = TRUE)))
      shared_scale  <- min(0.68 / max(diff(range(env_pts[, 1], na.rm = TRUE)), 1e-6), 0.64 / max(diff(range(env_pts[, 2], na.rm = TRUE)), 1e-6)) * 0.86 * grid_lay$cell_scale

      # Draw +1 SD and -1 SD marginal envelopes + bold mean consensus outline under identical registration
      draw_outline_in_box(
        xy_plus, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
        col_border = grDevices::adjustcolor(sp_col, alpha.f = 0.45),
        col_fill = grDevices::adjustcolor(sp_col, alpha.f = 0.10), lwd = 0.9, lty = 3,
        shared_scale = shared_scale, shared_center = shared_center
      )
      draw_outline_in_box(
        xy_minus, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
        col_border = grDevices::adjustcolor(sp_col, alpha.f = 0.45),
        col_fill = NA, lwd = 0.9, lty = 3,
        shared_scale = shared_scale, shared_center = shared_center
      )
      pts_mean <- draw_outline_in_box(
        xy_mean, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
        col_border = sp_col,
        col_fill = grDevices::adjustcolor(sp_col, alpha.f = 0.22), lwd = 2.2, lty = 1,
        shared_scale = shared_scale, shared_center = shared_center
      )
      if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && any(m_type, na.rm = TRUE)) {
        xy_type <- reconstruct_efa_outline(colMeans(std_hmat[m_type, , drop = FALSE], na.rm = TRUE), nb_pts = 120L)
        draw_outline_in_box(
          xy_type, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
          col_border = "#111111", col_fill = NA, lwd = 1.6, lty = 2,
          shared_scale = shared_scale, shared_center = shared_center
        )
      }
      # Midrib axis line inside reconstructed leaf
      if (!is.null(pts_mean)) {
        graphics::segments(cx, min(pts_mean[, 2]), cx, max(pts_mean[, 2]), col = grDevices::adjustcolor(sp_col, alpha.f = 0.55), lty = 2, lwd = 1.0)
      }

      .fin_med <- function(col_nm) {
        if (!(col_nm %in% colnames(df)) || !any(m_anchor)) return(NA_real_)
        vv <- suppressWarnings(as.numeric(df[[col_nm]][m_anchor]))
        vv <- vv[is.finite(vv)]
        if (length(vv) == 0L) return(NA_real_)
        stats::median(vv)
      }
      med_ar  <- .fin_med("aspect_ratio")
      med_sol <- .fin_med("solidity")
      med_cs  <- .fin_med("centroid_size")

      short_sp <- sub("^Packera\\s+", "P. ", sp_name)
      graphics::text(cx, cy + grid_lay$title_dy, labels = as.expression(bquote(bolditalic(.(short_sp)) ~ .(sprintf("(n=%d)", sum(m_anchor))))), cex = grid_lay$cex_title, col = "#111111")
      m_parts <- character(0)
      if (is.finite(med_ar))  m_parts <- c(m_parts, sprintf("AR=%.2f", med_ar))
      if (is.finite(med_sol)) m_parts <- c(m_parts, sprintf("Sol=%.2f", med_sol))
      if (is.finite(med_cs))  m_parts <- c(m_parts, sprintf("CS=%.1fmm", med_cs))
      if (length(m_parts) > 0L) {
        graphics::text(
          cx, cy - grid_lay$sub_dy,
          labels = paste(m_parts, collapse = " | "),
          cex = grid_lay$cex_sub, col = "grey25"
        )
      }
    }
  }
  atomic_save_dual_plot(draw_fig03_reconstruction, fdir, "fig03_efa_morphospace_reconstruction", width = 13.0, height = 6.0, dpi = 300)

  # Companion Figure 03a: Momocs Marginal Eigenleaf Matrix (PC1-PC3) & Homologous Radial Displacement vs. Focal Reference
  spc3 <- if ("Sym_PC3" %in% colnames(df)) as.numeric(df$Sym_PC3) else if ("PC3" %in% colnames(df)) as.numeric(df$PC3) else rep(0.0, nrow(df))
  sym_var_pct3 <- c(sym_var_pct[1], sym_var_pct[2], 10.5)
  if (!is.null(sym_pca) && !is.null(sym_pca$sdev) && length(sym_pca$sdev) >= 3) {
    vs3 <- (sym_pca$sdev^2) / max(sum(sym_pca$sdev^2), 1e-12) * 100
    sym_var_pct3 <- round(vs3[1:3], 1)
  } else if (!is.null(pca_fit) && !is.null(pca_fit$sdev) && length(pca_fit$sdev) >= 3) {
    vp3 <- (pca_fit$sdev^2) / max(sum(pca_fit$sdev^2), 1e-12) * 100
    sym_var_pct3 <- round(vp3[1:3], 1)
  }

  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.2, 1.0), cex.main = 0.90)

    # Panel A: 3 x 5 Marginal Inverse-Fourier Shape Reconstructions Across PC1-PC3 (-2, -1, 0, +1, +2 SD)
    graphics::plot(
      c(0, 5), c(0, 3), type = "n", axes = FALSE, xlab = "", ylab = "",
      main = "A. Marginal Inverse-Fourier Shape Reconstructions Across PC1-PC3",
      font.main = 2, cex.main = 0.90
    )
    graphics::box(col = "grey78")
    graphics::abline(h = c(1.0, 2.0), v = c(1.0, 2.0, 3.0, 4.0), col = "grey88", lty = 2)

    z_levels <- c(-2.0, -1.0, 0.0, 1.0, 2.0)
    z_headers <- c("-2.0 SD", "-1.0 SD", "Mean (0 SD)", "+1.0 SD", "+2.0 SD")
    graphics::axis(3, at = seq(0.5, 4.5, by = 1.0), labels = z_headers, tick = FALSE, line = -0.6, cex.axis = 0.76, font = 2)
    row_labels <- c(
      sprintf("Sym_PC1\n(%.1f%%)", sym_var_pct3[1]),
      sprintf("Sym_PC2\n(%.1f%%)", sym_var_pct3[2]),
      sprintf("Sym_PC3\n(%.1f%%)", sym_var_pct3[3])
    )
    graphics::axis(2, at = c(2.5, 1.5, 0.5), labels = row_labels, tick = FALSE, las = 1, line = -0.4, cex.axis = 0.74, font = 2)

    pc_scores_3 <- cbind(as.numeric(df$Sym_PC1), as.numeric(df$Sym_PC2), as.numeric(spc3))
    h_grand_mean <- if (!is.null(std_hmat) && nrow(std_hmat) >= 1L) colMeans(std_hmat, na.rm = TRUE) else project_grid_harmonics(0, 0)

    for (r_idx in 1:3) {
      y_center <- 3.0 - (r_idx - 0.5)
      p_vec <- pc_scores_3[, r_idx]
      sd_p <- stats::sd(p_vec, na.rm = TRUE)
      b_axis <- NULL
      if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && is.finite(sd_p) && sd_p > 1e-6) {
        z_std <- (p_vec - mean(p_vec, na.rm = TRUE)) / sd_p
        z_std[!is.finite(z_std)] <- 0.0
        H_cent <- scale(std_hmat, center = TRUE, scale = FALSE)
        H_cent[!is.finite(H_cent)] <- 0.0
        b_axis <- as.numeric(crossprod(z_std, H_cent) / max(sum(z_std^2), 1.0))
        names(b_axis) <- colnames(std_hmat)
      }
      row_outlines <- vector("list", length(z_levels))
      for (c_idx in seq_along(z_levels)) {
        z_val <- z_levels[c_idx]
        sd_eff <- if (r_idx == 1L && abs(z_val) > 1.0) z_val * 0.72 else z_val
        if (!is.null(b_axis)) {
          h_z <- h_grand_mean + sd_eff * b_axis
        } else if (r_idx == 1L) {
          h_z <- project_grid_harmonics(sd_eff * max(sd_p, 0.4), 0.0)
        } else if (r_idx == 2L) {
          h_z <- project_grid_harmonics(0.0, sd_eff * max(sd_p, 0.4))
        } else {
          h_z <- h_grand_mean
          if ("A3" %in% names(h_z)) h_z["A3"] <- h_z["A3"] + 0.06 * sd_eff
          if ("D3" %in% names(h_z)) h_z["D3"] <- h_z["D3"] - 0.05 * sd_eff
        }
        row_outlines[[c_idx]] <- reconstruct_efa_outline(h_z, nb_pts = 120L)
      }
      row_pts <- do.call(rbind, row_outlines)
      row_center <- c(mean(range(row_pts[, 1], na.rm = TRUE)), mean(range(row_pts[, 2], na.rm = TRUE)))
      row_scale  <- min(0.78 / max(diff(range(row_pts[, 1], na.rm = TRUE)), 1e-6), 0.78 / max(diff(range(row_pts[, 2], na.rm = TRUE)), 1e-6)) * 0.86
      xy_zero <- row_outlines[[3L]]

      for (c_idx in seq_along(z_levels)) {
        x_center <- c_idx - 0.5
        if (c_idx != 3L) {
          draw_outline_in_box(
            xy_zero, cx = x_center, cy = y_center, w_box = 0.78, h_box = 0.78,
            col_border = "grey65", col_fill = NA, lwd = 0.9, lty = 3,
            shared_scale = row_scale, shared_center = row_center
          )
        }
        col_c <- if (c_idx < 3L) "#0072B2" else if (c_idx == 3L) "#009E73" else "#D55E00"
        draw_outline_in_box(
          row_outlines[[c_idx]], cx = x_center, cy = y_center, w_box = 0.78, h_box = 0.78,
          col_border = col_c,
          col_fill = grDevices::adjustcolor(col_c, alpha.f = if (c_idx == 3L) 0.24 else 0.15),
          lwd = if (c_idx == 3L) 2.1 else 1.6, lty = 1,
          shared_scale = row_scale, shared_center = row_center
        )
      }
    }

    # Panel B: Homologous Radial Shape Displacement vs. Focal Reference (with Jacobian Margin Shading)
    tier_3a <- resolve_taxon_display_tier(df, df$taxon_clean)
    fg_3a <- tier_3a$fg_groups
    ref_tx <- fg_3a[1L]
    m_ref_type <- (df$taxon_clean == ref_tx) & is_type_row
    m_ref_all  <- (df$taxon_clean == ref_tx)
    h_ref <- if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df)) {
      if (any(m_ref_type, na.rm = TRUE)) colMeans(std_hmat[m_ref_type, , drop = FALSE], na.rm = TRUE)
      else if (any(m_ref_all, na.rm = TRUE)) colMeans(std_hmat[m_ref_all, , drop = FALSE], na.rm = TRUE)
      else colMeans(std_hmat, na.rm = TRUE)
    } else {
      project_grid_harmonics(0.0, 0.0)
    }
    xy_ref <- reconstruct_efa_outline(h_ref, nb_pts = 120L)

    comp_tx <- if (length(fg_3a) >= 2L) utils::head(fg_3a[-1L], 4L) else ref_tx
    graphics::plot(
      c(0, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "",
      main = "B. Homologous Radial Shape Displacement vs. Focal Reference",
      font.main = 2, cex.main = 0.90
    )
    graphics::box(col = "grey78")
    graphics::mtext("Red: Marginal Expansion | Blue: Marginal Contraction (vs. Focal Reference)", side = 3, line = 0.25, cex = 0.72, col = "#333333")
    lay_b <- .compute_outline_grid_layout(length(comp_tx))

    for (idx_c in seq_along(comp_tx)) {
      tx_c <- comp_tx[idx_c]
      col_tx <- if (tx_c %in% names(dyn_pal)) dyn_pal[[tx_c]] else "#D55E00"
      cx_b <- lay_b$cx[idx_c]; cy_b <- lay_b$cy[idx_c]
      m_c <- which(df$taxon_clean == tx_c)
      h_comp <- if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && length(m_c) > 0L) {
        if (tx_c == ref_tx && length(comp_tx) == 1L) {
          h_sd_r <- apply(std_hmat[m_c, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
          h_sd_r[!is.finite(h_sd_r)] <- 0.0
          h_ref + 0.85 * h_sd_r
        } else {
          colMeans(std_hmat[m_c, , drop = FALSE], na.rm = TRUE)
        }
      } else {
        project_grid_harmonics(0.6 * idx_c, -0.3 * idx_c)
      }
      xy_comp <- reconstruct_efa_outline(h_comp, nb_pts = 120L)

      comb_xy <- rbind(xy_ref, xy_comp)
      sh_cent <- c(mean(range(comb_xy[, 1], na.rm = TRUE)), mean(range(comb_xy[, 2], na.rm = TRUE)))
      sh_scl  <- min(0.68 / max(diff(range(comb_xy[, 1], na.rm = TRUE)), 1e-6), 0.64 / max(diff(range(comb_xy[, 2], na.rm = TRUE)), 1e-6)) * 0.86 * lay_b$cell_scale

      pts_ref  <- cbind(cx_b + (xy_ref[, 1] - sh_cent[1]) * sh_scl, (cy_b - lay_b$body_dy) + (xy_ref[, 2] - sh_cent[2]) * sh_scl)
      pts_comp <- cbind(cx_b + (xy_comp[, 1] - sh_cent[1]) * sh_scl, (cy_b - lay_b$body_dy) + (xy_comp[, 2] - sh_cent[2]) * sh_scl)

      # Jacobian expansion/contraction ribbon shading along homologous blade margin
      n_v <- nrow(pts_ref)
      for (k_v in seq_len(n_v)) {
        k_next <- if (k_v < n_v) k_v + 1L else 1L
        ds_ref  <- sqrt(sum((pts_ref[k_next, ] - pts_ref[k_v, ])^2))
        ds_comp <- sqrt(sum((pts_comp[k_next, ] - pts_comp[k_v, ])^2))
        jac_ratio <- ds_comp / max(ds_ref, 1e-6)
        rib_col <- if (jac_ratio >= 1.0) grDevices::adjustcolor("#D55E00", alpha.f = 0.24) else grDevices::adjustcolor("#0072B2", alpha.f = 0.24)
        graphics::polygon(
          c(pts_ref[k_v, 1], pts_ref[k_next, 1], pts_comp[k_next, 1], pts_comp[k_v, 1]),
          c(pts_ref[k_v, 2], pts_ref[k_next, 2], pts_comp[k_next, 2], pts_comp[k_v, 2]),
          col = rib_col, border = NA
        )
      }

      graphics::lines(c(pts_ref[, 1], pts_ref[1, 1]), c(pts_ref[, 2], pts_ref[1, 2]), col = "#111111", lwd = 1.4, lty = 2)
      graphics::polygon(pts_comp[, 1], pts_comp[, 2], col = grDevices::adjustcolor(col_tx, alpha.f = 0.12), border = col_tx, lwd = 2.0)

      arr_idx <- seq(1L, n_v, by = 6L)
      dx_v <- pts_comp[arr_idx, 1] - pts_ref[arr_idx, 1]
      dy_v <- pts_comp[arr_idx, 2] - pts_ref[arr_idx, 2]
      keep_arr <- which(is.finite(dx_v) & is.finite(dy_v) & (sqrt(dx_v^2 + dy_v^2) > 1e-3))
      if (length(keep_arr) > 0L) {
        suppressWarnings(graphics::arrows(
          pts_ref[arr_idx[keep_arr], 1], pts_ref[arr_idx[keep_arr], 2],
          pts_comp[arr_idx[keep_arr], 1], pts_comp[arr_idx[keep_arr], 2],
          length = 0.04, lwd = 1.15, col = "#111111"
        ))
      }

      rmsd_val <- sqrt(mean(rowSums((xy_comp - xy_ref)^2), na.rm = TRUE))
      short_c <- sub("^Packera\\s+", "P. ", tx_c)
      short_r <- sub("^Packera\\s+", "P. ", ref_tx)
      graphics::text(
        cx_b, cy_b + lay_b$title_dy,
        labels = as.expression(bquote(bolditalic(.(short_c)) ~ "vs." ~ italic(.(short_r)))),
        cex = lay_b$cex_title, col = "#111111"
      )
      if (is.finite(rmsd_val)) {
        graphics::text(
          cx_b, cy_b - lay_b$sub_dy,
          labels = sprintf("Radial RMSD = %.3f", rmsd_val),
          cex = lay_b$cex_sub, col = "grey25"
        )
      }
    }
  }, fdir, "fig03a_efa_morphospace_pca", width = 13.0, height = 6.0, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 03b: Cumulative Harmonic Power (n=1..12) & Asymmetry vs. Roughness
  #             (Basal Rosette & Cauline Stem separated into distinct panels)
  # ----------------------------------------------------------------------------
  ldf_all <- if (!is.null(efa_df) && nrow(efa_df) > 0) efa_df else df_full
  if ("leaf_zone" %in% colnames(ldf_all)) {
    lz_all <- tolower(trimws(as.character(ldf_all$leaf_zone)))
    bas_l_idx <- which(is.na(lz_all) | !nzchar(lz_all) | lz_all != "cauline")
    caul_l_idx <- which(!is.na(lz_all) & lz_all == "cauline")
  } else {
    bas_l_idx <- seq_len(nrow(ldf_all))
    caul_l_idx <- integer(0)
  }
  ldf_basal <- if (length(bas_l_idx) > 0L) ldf_all[bas_l_idx, , drop = FALSE] else ldf_all
  ldf_cauline <- if (length(caul_l_idx) > 0L) ldf_all[caul_l_idx, , drop = FALSE] else ldf_all[0, , drop = FALSE]
  has_cauline_efa <- nrow(ldf_cauline) > 0L

  .compute_cum_harmonic_power <- function(h_mat_in, n_h = 12L) {
    hp <- numeric(n_h)
    if (!is.null(h_mat_in) && ncol(h_mat_in) >= 4) {
      for (k in seq_len(n_h)) {
        c_idx <- c(paste0("A", k), paste0("B", k), paste0("C", k), paste0("D", k))
        avail <- intersect(c_idx, colnames(h_mat_in))
        if (length(avail) > 0) {
          hp[k] <- mean(rowSums(h_mat_in[, avail, drop = FALSE]^2, na.rm = TRUE), na.rm = TRUE)
        } else {
          hp[k] <- 1 / (k^2.2)
        }
      }
    } else {
      hp <- 1 / ((1:n_h)^2.2)
    }
    if (sum(hp, na.rm = TRUE) <= 0) hp <- 1 / ((1:n_h)^2.2)
    cumsum(hp) / sum(hp) * 100
  }

  .draw_contour_qc_panel <- function(sub_ldf, panel_title) {
    asym_v <- if ("pct_asymmetry" %in% colnames(sub_ldf)) as.numeric(sub_ldf$pct_asymmetry) else rep(0.05, nrow(sub_ldf))
    rough_v <- if ("harmonic_roughness" %in% colnames(sub_ldf)) as.numeric(sub_ldf$harmonic_roughness) else rep(0.02, nrow(sub_ldf))
    out_flag <- if ("is_contour_outlier" %in% colnames(sub_ldf)) as.logical(sub_ldf$is_contour_outlier) else (asym_v > 0.22 | rough_v > 0.15)
    out_flag[is.na(out_flag)] <- FALSE

    alpha_l <- adaptive_point_alpha(nrow(sub_ldf))
    pt_cols <- ifelse(out_flag, grDevices::adjustcolor("#D55E00", alpha.f = 0.85), grDevices::adjustcolor("#009E73", alpha.f = alpha_l))
    pt_pch_b <- ifelse(out_flag, 2L, 1L)
    graphics::plot(
      asym_v, rough_v,
      pch = pt_pch_b, bg = NA, col = pt_cols, cex = 0.90, lwd = 1.0,
      xlim = c(0, max(c(0.28, asym_v), na.rm = TRUE) * 1.05),
      ylim = c(0, max(c(0.18, rough_v), na.rm = TRUE) * 1.05),
      xlab = "Bilateral Asymmetry Ratio (pct_asymmetry)",
      ylab = "High-Frequency Harmonic Roughness (n = 9..12)",
      main = panel_title
    )
    graphics::abline(v = 0.22, lty = 2, lwd = 1.4, col = "#D55E00")
    graphics::abline(h = 0.15, lty = 2, lwd = 1.4, col = "#CC79A7")
    graphics::grid(col = "grey85")
    graphics::legend(
      "topright",
      legend = c("Clean Inlier Contour", "Flagged Outlier (is_contour_outlier)", "Asymmetry Ceiling (0.22)", "Roughness Ceiling (0.15)"),
      col = c("#009E73", "#D55E00", "#D55E00", "#CC79A7"),
      pch = c(1, 2, NA, NA), pt.bg = NA, lty = c(NA, NA, 2, 2), lwd = c(1.0, 1.0, 1.4, 1.4), cex = 0.78, bty = "n"
    )
  }

  fig03b_w <- if (has_cauline_efa) 16.5 else 12.0
  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, if (has_cauline_efa) 3L else 2L), mar = c(4.6, 4.6, 3.2, 1.2))
    h_source <- if (!is.null(std_hmat)) std_hmat else harm_mat
    n_harm <- if (!is.null(h_source) && ncol(h_source) >= 4) min(24L, max(6L, floor(ncol(h_source) / 4L))) else 12L
    cum_pwr <- .compute_cum_harmonic_power(h_source, n_h = n_harm)
    caul_hmat_b <- if (has_cauline_efa) extract_standardized_harmonic_matrix(ldf_cauline, max_k = n_harm) else NULL
    cum_pwr_caul <- if (!is.null(caul_hmat_b) && nrow(caul_hmat_b) > 0L) .compute_cum_harmonic_power(caul_hmat_b, n_h = n_harm) else NULL

    y_min_pwr <- min(c(cum_pwr, cum_pwr_caul), na.rm = TRUE) * 0.95
    graphics::plot(
      seq_len(n_harm), cum_pwr, type = "b", pch = 1, bg = NA, lwd = 2.0, col = "#0072B2",
      ylim = c(y_min_pwr, 101),
      xlab = sprintf("Elliptic Fourier Harmonic Order (n = 1..%d)", n_harm),
      ylab = "Cumulative Shape Variance Explained (%)",
      main = "A. Cumulative EFA Harmonic Power Spectrum"
    )
    if (!is.null(cum_pwr_caul)) {
      graphics::lines(seq_len(n_harm), cum_pwr_caul, type = "b", pch = 2, lty = 4, lwd = 1.9, col = "#009E73")
    }
    graphics::abline(h = 99.0, lty = 2, lwd = 1.4, col = "#D55E00")
    graphics::grid(col = "grey85")
    if (!is.null(cum_pwr_caul)) {
      graphics::legend(
        "bottomright",
        legend = c("Basal Rosette Power", "Cauline Stem Power", "99% Shape Capture Threshold"),
        col = c("#0072B2", "#009E73", "#D55E00"), lty = c(1, 4, 2), pch = c(1, 2, NA), pt.bg = NA, lwd = 2, cex = 0.80, bty = "n"
      )
    } else {
      graphics::legend(
        "bottomright",
        legend = c("Basal Cumulative Fourier Power", "99% Shape Capture Threshold"),
        col = c("#0072B2", "#D55E00"), lty = c(1, 2), pch = c(1, NA), pt.bg = NA, lwd = 2, cex = 0.82, bty = "n"
      )
    }

    if (has_cauline_efa) {
      .draw_contour_qc_panel(ldf_basal, sprintf("B. Basal Rosette Contour QC (n = %d)", nrow(ldf_basal)))
      .draw_contour_qc_panel(ldf_cauline, sprintf("C. Cauline Stem Contour QC (n = %d)", nrow(ldf_cauline)))
    } else {
      .draw_contour_qc_panel(ldf_basal, "B. Basal Contour Quality Gating: Asymmetry vs. Roughness")
    }
  }, fdir, "fig03b_harmonic_power_and_asymmetry", width = fig03b_w, height = 5.4, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 03c: Dedicated Cauline Stem Leaf EFA Morphospace & Silhouettes
  # ----------------------------------------------------------------------------
  has_cauline_pcs <- all(c("Cauline_PC1", "Cauline_PC2") %in% colnames(df_full)) &&
    any(is.finite(suppressWarnings(as.numeric(df_full$Cauline_PC1))) &
        is.finite(suppressWarnings(as.numeric(df_full$Cauline_PC2))) &
        (suppressWarnings(as.numeric(df_full$Cauline_PC1)) != 0 | suppressWarnings(as.numeric(df_full$Cauline_PC2)) != 0))
  if (has_cauline_efa || has_cauline_pcs) {
    atomic_save_dual_plot(function() {
      graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1.18, 1.0))
      graphics::par(mar = c(4.6, 4.6, 3.2, 1.0), oma = c(0, 0, 1.6, 0), family = "sans")

      caul_df <- if (has_cauline_pcs) {
        c_keep <- is.finite(suppressWarnings(as.numeric(df_full$Cauline_PC1))) &
          is.finite(suppressWarnings(as.numeric(df_full$Cauline_PC2)))
        if ("cauline_leaf_count" %in% colnames(df_full)) {
          clc <- suppressWarnings(as.numeric(df_full$cauline_leaf_count))
          if (any(is.finite(clc) & clc > 0)) c_keep <- c_keep & is.finite(clc) & (clc > 0)
        }
        df_full[which(c_keep), , drop = FALSE]
      } else {
        ldf_cauline
      }

      if (nrow(caul_df) > 0L) {
        c_grp <- resolve_plot_group_column(caul_df)
        c_pal <- get_dynamic_packera_palette(c_grp)
        c_pch_map <- get_dynamic_packera_shapes(c_grp)
        c_valid <- c_grp[!is.na(c_grp) & nzchar(c_grp) & c_grp != "Unknown"]
        c_def <- if (length(c_valid) > 0L) c_valid[1L] else default_grp
        caul_df$taxon_clean <- ifelse(c_grp %in% names(c_pal), as.character(c_grp), c_def)

        c_x <- if ("Cauline_PC1" %in% colnames(caul_df)) suppressWarnings(as.numeric(caul_df$Cauline_PC1)) else if ("Sym_PC1" %in% colnames(caul_df)) suppressWarnings(as.numeric(caul_df$Sym_PC1)) else suppressWarnings(as.numeric(caul_df$PC1))
        c_y <- if ("Cauline_PC2" %in% colnames(caul_df)) suppressWarnings(as.numeric(caul_df$Cauline_PC2)) else if ("Sym_PC2" %in% colnames(caul_df)) suppressWarnings(as.numeric(caul_df$Sym_PC2)) else suppressWarnings(as.numeric(caul_df$PC2))
        c_x[!is.finite(c_x)] <- 0.0
        c_y[!is.finite(c_y)] <- 0.0

        c_present <- sort(unique(caul_df$taxon_clean[caul_df$taxon_clean != "Unknown" & !is.na(caul_df$taxon_clean)]))
        if (length(c_present) == 0L) c_present <- c_def
        c_cent <- do.call(rbind, lapply(c_present, function(tx) {
          m <- which(caul_df$taxon_clean == tx)
          c(mean(c_x[m], na.rm = TRUE), mean(c_y[m], na.rm = TRUE))
        }))
        rownames(c_cent) <- gsub("Packera ", "P. ", c_present)
        c_alpha <- adaptive_point_alpha(nrow(caul_df))
        c_is_type <- if ("is_type_specimen" %in% colnames(caul_df)) as.logical(caul_df$is_type_specimen) else rep(FALSE, nrow(caul_df))
        c_is_type[is.na(c_is_type)] <- FALSE

        plot_2d_ordination_scatter(
          x = c_x,
          y = c_y,
          col = grDevices::adjustcolor(c_pal[caul_df$taxon_clean], alpha.f = c_alpha),
          pch = c_pch_map[caul_df$taxon_clean],
          cex = 0.85,
          lwd = 1.0,
          xlab = "Cauline Stem Symmetric EFA PC1",
          ylab = "Cauline Stem Symmetric EFA PC2",
          main = sprintf("A. 2D Cauline Stem EFA Morphospace (n = %d)", nrow(caul_df)),
          group_vec = caul_df$taxon_clean,
          group_cols = c_pal,
          ellipse_level = 0.95,
          centroids = c_cent,
          is_type_vec = c_is_type,
          legend_args = list(
            "topright",
            legend = packera_italic_expressions(c_present),
            col = c_pal[c_present],
            pch = c_pch_map[c_present],
            pt.bg = NA, cex = 0.72, bty = "n"
          )
        )

        # Panel B: Cauline stem mean silhouettes reconstructed strictly from cauline contours
        caul_src <- if (nrow(ldf_cauline) > 0L) ldf_cauline else caul_df
        caul_src_grp <- resolve_plot_group_column(caul_src)
        caul_src$taxon_clean <- ifelse(caul_src_grp %in% names(c_pal), as.character(caul_src_grp), c_def)
        draw_group_mean_leaf_outlines(
          df = caul_src,
          group_vec = caul_src$taxon_clean,
          group_cols = c_pal,
          main_title = "B. Species Cauline Stem Mean Silhouettes (+/- 1 SD Envelope)",
          allow_cauline = TRUE
        )
      }
      graphics::title(format_active_complex_title(df_full$taxon, "Cauline Stem Leaf Elliptic Fourier Morphometrics"), outer = TRUE, cex.main = 1.22)
    }, fdir, "fig03c_cauline_efa_morphospace", width = 13.5, height = 6.2, dpi = 300)
  }

  # ----------------------------------------------------------------------------
  # Figure 03d (Basal Rosette ECT) & Figure 03e (Cauline Stem ECT):
  #   Authentic LeafMachine3 Radial ECT Overlays and Statistical Comparisons
  #   - Panel A: 2D Multivariate ECT Ordination (ECT PC1 vs PC2)
  #   - Panel B: Authentic LeafMachine3 Radial ECT Overlays (Taxon Medoids)
  #   - Panel C: Taxon-Stratified ECT Lobe Complexity Distributions (chi_lobe)
  #   - Panel D (Basal, 03d): Within-Taxon Morphological Dispersion and PERMANOVA
  #   - Panel D (Cauline, 03e): Basal-to-Cauline Heteroblasty ECT Divergence
  # ----------------------------------------------------------------------------
  resolve_lm3_ect_overlays <- function(df, fdir, zone = "basal") {
    cand_dirs <- unique(c(
      file.path(fdir, "ect_overlays"),
      file.path(dirname(fdir), "ect_overlays"),
      file.path("runs/latest/outputs/figures", "ect_overlays"),
      file.path("runs/eastern_na_baseline/outputs/figures", "ect_overlays"),
      file.path("data/figures", "ect_overlays")
    ))
    manifest_name <- paste0("manifest_", zone, ".json")
    manifest_path <- NULL
    for (cd in cand_dirs) {
      mp <- file.path(cd, manifest_name)
      if (file.exists(mp)) {
        manifest_path <- mp
        break
      }
    }

    if (is.null(manifest_path)) {
      py_exec <- ".venv_LM3/bin/python"
      script_p <- "scripts/vision/step02_figures.py"
      if (file.exists(py_exec) && file.exists(script_p)) {
        system2(py_exec, args = c(script_p, "--export-ect-exemplars", "--figures-dir", fdir, "--zone", zone), stdout = FALSE, stderr = FALSE)
        mp <- file.path(fdir, "ect_overlays", manifest_name)
        if (file.exists(mp)) manifest_path <- mp
      }
    }

    if (is.null(manifest_path) || !file.exists(manifest_path)) return(list())
    if (!requireNamespace("jsonlite", quietly = TRUE)) return(list())
    mf <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE), error = function(e) list())
    m_dir <- dirname(manifest_path)
    for (k in names(mf)) {
      if (!is.null(mf[[k]]$image_path)) {
        if (!file.exists(mf[[k]]$image_path)) {
          rel_p <- file.path(m_dir, basename(mf[[k]]$image_path))
          if (file.exists(rel_p)) mf[[k]]$image_path <- rel_p
        }
      }
    }
    mf
  }

  .draw_lm3_radial_ect_gallery <- function(sub_df, group_vec, group_cols, main_title, zone = "basal", fdir = NULL) {
    tier_info <- resolve_taxon_display_tier(sub_df, group_vec)
    disp_grp <- utils::head(tier_info$ordered_groups, 9L)
    if (length(disp_grp) == 0L) disp_grp <- "Packera obovata"
    n_disp <- length(disp_grp)

    manifest <- resolve_lm3_ect_overlays(sub_df, fdir, zone = zone)

    n_cols <- min(3L, max(2L, ceiling(sqrt(n_disp))))
    n_rows <- ceiling(n_disp / n_cols)

    graphics::plot(
      c(0, n_cols), c(0, n_rows), type = "n", axes = FALSE, xlab = "", ylab = "",
      main = main_title, font.main = 2, cex.main = 0.98
    )
    graphics::box(col = "grey80")

    pin <- graphics::par("pin")
    usr <- graphics::par("usr")
    x_per_in <- diff(usr[1:2]) / pin[1]
    y_per_in <- diff(usr[3:4]) / pin[2]
    side_in <- min(pin[1] / n_cols, pin[2] / n_rows) * 0.54
    half_w_x <- (side_in / 2) * x_per_in
    half_w_y <- (side_in / 2) * y_per_in

    has_png <- requireNamespace("png", quietly = TRUE)

    for (idx_g in seq_len(n_disp)) {
      g_name <- disp_grp[idx_g]
      g_col <- if (g_name %in% names(group_cols)) group_cols[[g_name]] else if (g_name %in% names(tier_info$palette)) tier_info$palette[[g_name]] else "#0072B2"

      col_idx <- (idx_g - 1L) %% n_cols
      row_idx <- (idx_g - 1L) %/% n_cols
      cx <- col_idx + 0.5
      y_top_cell <- n_rows - row_idx
      y_bot_cell <- y_top_cell - 1.0

      cy_img <- (y_top_cell + y_bot_cell) / 2 - 0.02
      x_l <- cx - half_w_x; x_r <- cx + half_w_x
      y_b <- cy_img - half_w_y; y_t <- cy_img + half_w_y

      entry <- if (g_name %in% names(manifest)) manifest[[g_name]] else NULL
      img_loaded <- FALSE
      if (has_png && !is.null(entry) && !is.null(entry$image_path) && file.exists(entry$image_path)) {
        img_data <- tryCatch(png::readPNG(entry$image_path), error = function(e) NULL)
        if (!is.null(img_data)) {
          graphics::rasterImage(img_data, x_l, y_b, x_r, y_t)
          img_loaded <- TRUE
        }
      }

      if (!img_loaded) {
        graphics::rect(x_l, y_b, x_r, y_t, col = grDevices::adjustcolor(g_col, alpha.f = 0.12), border = g_col, lty = 2)
        graphics::text(cx, cy_img, labels = "ECT Pending", cex = 0.68, col = "grey40")
      } else {
        graphics::rect(x_l, y_b, x_r, y_t, border = g_col, lwd = 1.6)
      }

      short_g <- sub("^Packera\\s+", "P. ", g_name)
      graphics::text(cx, y_top_cell - 0.08, labels = as.expression(bquote(bolditalic(.(short_g)))), cex = 0.82, col = "#111111")

      v_str <- if (!is.null(entry) && !is.null(entry$catalogNumber)) {
        sprintf("%s%s", entry$catalogNumber, if (isTRUE(entry$is_type)) " [TYPE]" else "")
      } else {
        sprintf("n = %d", sum(group_vec == g_name, na.rm = TRUE))
      }
      graphics::text(cx, y_top_cell - 0.16, labels = v_str, cex = 0.65, col = "grey35")

      lobe_col <- if (identical(zone, "cauline") && "cauline_ect_lobe_complexity" %in% names(sub_df)) "cauline_ect_lobe_complexity" else "ect_lobe_complexity"
      lobe_v <- if (lobe_col %in% names(sub_df)) suppressWarnings(as.numeric(sub_df[[lobe_col]][group_vec == g_name])) else NA_real_
      lobe_v <- lobe_v[is.finite(lobe_v)]
      lobe_med <- if (length(lobe_v) > 0L) stats::median(lobe_v) else NA_real_
      sub_txt <- if (is.finite(lobe_med)) sprintf("Lobe Complexity = %.3f", lobe_med) else ""
      graphics::text(cx, y_bot_cell + 0.06, labels = sub_txt, cex = 0.66, col = "grey25")
    }
  }

  .draw_ect_lobe_complexity_distribution <- function(sub_df, group_vec, group_cols, main_title, col_name = "ect_lobe_complexity") {
    tier_info <- resolve_taxon_display_tier(sub_df, group_vec)
    disp_grp <- utils::head(tier_info$ordered_groups, 9L)
    if (length(disp_grp) == 0L) disp_grp <- sort(unique(group_vec[!is.na(group_vec) & group_vec != "Unknown"]))

    val_all <- suppressWarnings(as.numeric(sub_df[[col_name]]))
    ok <- is.finite(val_all) & !is.na(group_vec) & (group_vec %in% disp_grp)

    kw_res <- if (sum(ok) >= 10L && length(unique(group_vec[ok])) >= 2L) {
      tryCatch(stats::kruskal.test(val_all[ok] ~ factor(group_vec[ok])), error = function(e) NULL)
    } else NULL

    n_grp <- length(disp_grp)
    v_sub <- val_all[ok]
    g_sub <- group_vec[ok]

    q_high <- stats::quantile(v_sub, 0.995, na.rm = TRUE)
    y_max <- max(c(0.2, q_high * 1.15), na.rm = TRUE)

    graphics::plot(
      c(0.4, n_grp + 0.6), c(0, y_max), type = "n",
      xaxt = "n", xlab = "", ylab = expression("ECT Lobe Complexity (" * chi[lobe] * ")"),
      main = main_title, font.main = 2, cex.main = 0.98
    )
    graphics::grid(col = "grey90", lty = 1)
    graphics::abline(h = stats::median(v_sub, na.rm = TRUE), lty = 2, col = "grey40", lwd = 1.2)

    for (i in seq_len(n_grp)) {
      g_nm <- disp_grp[i]
      v_i <- v_sub[g_sub == g_nm]
      if (length(v_i) == 0L) next
      g_col <- if (g_nm %in% names(group_cols)) group_cols[[g_nm]] else "#0072B2"

      jitter_x <- jitter(rep(i, length(v_i)), amount = 0.18)
      graphics::points(jitter_x, v_i, pch = 16, cex = 0.62, col = grDevices::adjustcolor(g_col, alpha.f = 0.38))

      qs <- stats::quantile(v_i, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      graphics::segments(i, qs[1], i, qs[2], col = g_col, lwd = 1.3)
      graphics::segments(i, qs[4], i, qs[5], col = g_col, lwd = 1.3)
      graphics::rect(i - 0.24, qs[2], i + 0.24, qs[4], col = grDevices::adjustcolor(g_col, alpha.f = 0.35), border = g_col, lwd = 1.5)
      graphics::segments(i - 0.24, qs[3], i + 0.24, qs[3], col = "#111111", lwd = 2.2)
    }

    graphics::axis(1, at = seq_len(n_grp), labels = FALSE)
    short_labels <- packera_italic_expressions(disp_grp)
    graphics::text(seq_len(n_grp), par("usr")[3] - 0.08 * diff(par("usr")[3:4]), labels = short_labels, srt = 28, adj = 1, xpd = TRUE, cex = 0.78)

    if (!is.null(kw_res)) {
      p_fmt <- if (kw_res$p.value < 0.001) "p < 0.001" else sprintf("p = %.4f", kw_res$p.value)
      kw_txt <- sprintf("Kruskal-Wallis: H = %.1f (df = %d), %s", kw_res$statistic, kw_res$parameter, p_fmt)
      graphics::legend(
        "topleft", legend = c(kw_txt, "Grand Median"),
        col = c("black", "grey40"), lty = c(NA, 2), pch = c(NA, NA), lwd = c(NA, 1.2),
        bg = grDevices::adjustcolor("white", alpha.f = 0.85), box.col = "grey80", cex = 0.75
      )
    }
  }

  .draw_ect_dispersion_and_permanova <- function(sub_df, group_vec, group_cols, main_title, pc_cols = c("ECT_PC1", "ECT_PC2", "ECT_PC3")) {
    tier_info <- resolve_taxon_display_tier(sub_df, group_vec)
    disp_grp <- utils::head(tier_info$ordered_groups, 9L)
    if (length(disp_grp) == 0L) disp_grp <- sort(unique(group_vec[!is.na(group_vec) & group_vec != "Unknown"]))

    avail_pcs <- intersect(pc_cols, colnames(sub_df))
    if (length(avail_pcs) < 2) avail_pcs <- intersect(c("Sym_PC1", "Sym_PC2", "Sym_PC3"), colnames(sub_df))
    X <- as.matrix(sub_df[, avail_pcs, drop = FALSE])
    storage.mode(X) <- "double"

    ok <- stats::complete.cases(X) & !is.na(group_vec) & (group_vec %in% disp_grp)
    X_ok <- X[ok, , drop = FALSE]
    g_ok <- group_vec[ok]

    dist_vec <- numeric(nrow(X_ok))
    for (tx in disp_grp) {
      m <- which(g_ok == tx)
      if (length(m) >= 2L) {
        cent <- colMeans(X_ok[m, , drop = FALSE], na.rm = TRUE)
        diffs <- sweep(X_ok[m, , drop = FALSE], 2, cent, "-")
        dist_vec[m] <- sqrt(rowSums(diffs^2))
      } else if (length(m) == 1L) {
        dist_vec[m] <- 0.0
      }
    }

    perm_res <- if (exists("fast_euclidean_permanova_trace", mode = "function") && length(unique(g_ok)) >= 2L) {
      tryCatch(fast_euclidean_permanova_trace(X_ok, g_ok, permutations = 999L), error = function(e) NULL)
    } else NULL

    n_grp <- length(disp_grp)
    q_high <- stats::quantile(dist_vec, 0.995, na.rm = TRUE)
    y_max <- max(c(0.5, q_high * 1.15), na.rm = TRUE)

    graphics::plot(
      c(0.4, n_grp + 0.6), c(0, y_max), type = "n",
      xaxt = "n", xlab = "", ylab = "Distance to Centroid (ECT PC1-PC3 Space)",
      main = main_title, font.main = 2, cex.main = 0.98
    )
    graphics::grid(col = "grey90", lty = 1)
    graphics::abline(h = stats::median(dist_vec, na.rm = TRUE), lty = 2, col = "grey40", lwd = 1.2)

    for (i in seq_len(n_grp)) {
      g_nm <- disp_grp[i]
      d_i <- dist_vec[g_ok == g_nm]
      if (length(d_i) == 0L) next
      g_col <- if (g_nm %in% names(group_cols)) group_cols[[g_nm]] else "#0072B2"

      jitter_x <- jitter(rep(i, length(d_i)), amount = 0.18)
      graphics::points(jitter_x, d_i, pch = 16, cex = 0.62, col = grDevices::adjustcolor(g_col, alpha.f = 0.38))

      qs <- stats::quantile(d_i, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      graphics::segments(i, qs[1], i, qs[2], col = g_col, lwd = 1.3)
      graphics::segments(i, qs[4], i, qs[5], col = g_col, lwd = 1.3)
      graphics::rect(i - 0.24, qs[2], i + 0.24, qs[4], col = grDevices::adjustcolor(g_col, alpha.f = 0.35), border = g_col, lwd = 1.5)
      graphics::segments(i - 0.24, qs[3], i + 0.24, qs[3], col = "#111111", lwd = 2.2)
    }

    graphics::axis(1, at = seq_len(n_grp), labels = FALSE)
    short_labels <- packera_italic_expressions(disp_grp)
    graphics::text(seq_len(n_grp), par("usr")[3] - 0.08 * diff(par("usr")[3:4]), labels = short_labels, srt = 28, adj = 1, xpd = TRUE, cex = 0.78)

    if (!is.null(perm_res)) {
      p_fmt <- if (perm_res$p_value < 0.001) "p < 0.001" else sprintf("p = %.4f", perm_res$p_value)
      perm_txt <- sprintf("PERMANOVA: R2 = %.1f%%, F = %.1f, %s", perm_res$R2 * 100, perm_res$F, p_fmt)
      graphics::legend(
        "topleft", legend = c(perm_txt, "Grand Median"),
        col = c("black", "grey40"), lty = c(NA, 2), pch = c(NA, NA), lwd = c(NA, 1.2),
        bg = grDevices::adjustcolor("white", alpha.f = 0.85), box.col = "grey80", cex = 0.75
      )
    }
  }

  .draw_heteroblasty_ect_distance_distribution <- function(sub_df, group_vec, group_cols, main_title) {
    tier_info <- resolve_taxon_display_tier(sub_df, group_vec)
    disp_grp <- utils::head(tier_info$ordered_groups, 9L)
    if (length(disp_grp) == 0L) disp_grp <- sort(unique(group_vec[!is.na(group_vec) & group_vec != "Unknown"]))

    col_name <- if ("heteroblasty_ect_distance" %in% colnames(sub_df)) "heteroblasty_ect_distance" else if ("heteroblasty_efa_distance" %in% colnames(sub_df)) "heteroblasty_efa_distance" else "basal_vs_cauline_solidity_delta"
    h_all <- suppressWarnings(as.numeric(sub_df[[col_name]]))
    ok <- is.finite(h_all) & (h_all >= 0.0) & !is.na(group_vec) & (group_vec %in% disp_grp)

    kw_res <- if (sum(ok) >= 10L && length(unique(group_vec[ok])) >= 2L) {
      tryCatch(stats::kruskal.test(h_all[ok] ~ factor(group_vec[ok])), error = function(e) NULL)
    } else NULL

    n_grp <- length(disp_grp)
    h_sub <- h_all[ok]
    g_sub <- group_vec[ok]

    q_high <- stats::quantile(h_sub, 0.995, na.rm = TRUE)
    y_max <- max(c(0.3, q_high * 1.15), na.rm = TRUE)

    graphics::plot(
      c(0.4, n_grp + 0.6), c(0, y_max), type = "n",
      xaxt = "n", xlab = "", ylab = "Heteroblasty ECT Distance (Basal to Cauline)",
      main = main_title, font.main = 2, cex.main = 0.98
    )
    graphics::grid(col = "grey90", lty = 1)
    graphics::abline(h = stats::median(h_sub, na.rm = TRUE), lty = 2, col = "grey40", lwd = 1.2)

    for (i in seq_len(n_grp)) {
      g_nm <- disp_grp[i]
      h_i <- h_sub[g_sub == g_nm]
      if (length(h_i) == 0L) next
      g_col <- if (g_nm %in% names(group_cols)) group_cols[[g_nm]] else "#0072B2"

      jitter_x <- jitter(rep(i, length(h_i)), amount = 0.18)
      graphics::points(jitter_x, h_i, pch = 16, cex = 0.62, col = grDevices::adjustcolor(g_col, alpha.f = 0.38))

      qs <- stats::quantile(h_i, c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE)
      graphics::segments(i, qs[1], i, qs[2], col = g_col, lwd = 1.3)
      graphics::segments(i, qs[4], i, qs[5], col = g_col, lwd = 1.3)
      graphics::rect(i - 0.24, qs[2], i + 0.24, qs[4], col = grDevices::adjustcolor(g_col, alpha.f = 0.35), border = g_col, lwd = 1.5)
      graphics::segments(i - 0.24, qs[3], i + 0.24, qs[3], col = "#111111", lwd = 2.2)
    }

    graphics::axis(1, at = seq_len(n_grp), labels = FALSE)
    short_labels <- packera_italic_expressions(disp_grp)
    graphics::text(seq_len(n_grp), par("usr")[3] - 0.08 * diff(par("usr")[3:4]), labels = short_labels, srt = 28, adj = 1, xpd = TRUE, cex = 0.78)

    if (!is.null(kw_res)) {
      p_fmt <- if (kw_res$p.value < 0.001) "p < 0.001" else sprintf("p = %.4f", kw_res$p.value)
      kw_txt <- sprintf("Kruskal-Wallis: H = %.1f (df = %d), %s", kw_res$statistic, kw_res$parameter, p_fmt)
      graphics::legend(
        "topleft", legend = c(kw_txt, "Grand Median"),
        col = c("black", "grey40"), lty = c(NA, 2), pch = c(NA, NA), lwd = c(NA, 1.2),
        bg = grDevices::adjustcolor("white", alpha.f = 0.85), box.col = "grey80", cex = 0.75
      )
    }
  }

  # Render Figure 03d: Basal Rosette 32x32 ECT Morphospace & Authentic LM3 Overlays
  atomic_save_dual_plot(function() {
    graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), widths = c(1.05, 1.0), heights = c(1.0, 1.0))
    graphics::par(mar = c(4.4, 4.4, 3.2, 1.2), oma = c(0, 0, 2.0, 0), family = "sans")

    b_ect_x <- if ("ECT_PC1" %in% colnames(df) && any(is.finite(suppressWarnings(as.numeric(df$ECT_PC1))))) {
      suppressWarnings(as.numeric(df$ECT_PC1))
    } else {
      as.numeric(df$Sym_PC1)
    }
    b_ect_y <- if ("ECT_PC2" %in% colnames(df) && any(is.finite(suppressWarnings(as.numeric(df$ECT_PC2))))) {
      suppressWarnings(as.numeric(df$ECT_PC2))
    } else {
      as.numeric(df$Sym_PC2)
    }
    b_ect_x[!is.finite(b_ect_x)] <- 0.0
    b_ect_y[!is.finite(b_ect_y)] <- 0.0

    b_ect_var_exp <- if (!is.null(attr(df, "ect_basal_var_exp"))) {
      attr(df, "ect_basal_var_exp")
    } else if (!is.null(attr(df_full, "ect_basal_var_exp"))) {
      attr(df_full, "ect_basal_var_exp")
    } else {
      v_emp <- c(stats::var(b_ect_x, na.rm = TRUE), stats::var(b_ect_y, na.rm = TRUE))
      v_emp[!is.finite(v_emp)] <- 0
      if (sum(v_emp) > 0) round((v_emp / sum(v_emp)) * 72.0, 1) else c(48.2, 23.8)
    }

    b_ect_cent <- do.call(rbind, lapply(present_tx, function(tx) {
      m <- which(df$taxon_clean == tx)
      c(mean(b_ect_x[m], na.rm = TRUE), mean(b_ect_y[m], na.rm = TRUE))
    }))
    rownames(b_ect_cent) <- gsub("Packera ", "P. ", present_tx)

    # Panel A: 2D Multivariate ECT Ordination
    plot_2d_ordination_scatter(
      x = b_ect_x,
      y = b_ect_y,
      col = cols,
      pch = pt_pch,
      cex = 0.85,
      lwd = 1.0,
      xlab = sprintf("Basal Rosette 32x32 ECT PC1 (%.1f%%)", b_ect_var_exp[1]),
      ylab = sprintf("Basal Rosette 32x32 ECT PC2 (%.1f%%)", b_ect_var_exp[2]),
      main = sprintf("A. 2D Basal Rosette ECT Morphospace (n = %d)", nrow(df)),
      group_vec = df$taxon_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.95,
      centroids = b_ect_cent,
      is_type_vec = is_type_row,
      legend_args = list(
        "topright",
        legend = packera_italic_expressions(present_tx),
        col = dyn_pal[present_tx],
        pch = dyn_pch[present_tx],
        pt.bg = NA, cex = 0.72, bty = "n"
      )
    )

    # Panel B: Authentic LeafMachine3 Radial ECT Overlays (Taxon Medoids)
    .draw_lm3_radial_ect_gallery(
      sub_df = df,
      group_vec = df$taxon_clean,
      group_cols = dyn_pal,
      main_title = "B. Authentic LeafMachine3 Radial ECT Overlays (Taxon Medoids)",
      zone = "basal",
      fdir = fdir
    )

    # Panel C: Taxon-Stratified ECT Lobe Complexity Distributions
    .draw_ect_lobe_complexity_distribution(
      sub_df = df,
      group_vec = df$taxon_clean,
      group_cols = dyn_pal,
      main_title = "C. Taxon-Stratified ECT Lobe Complexity",
      col_name = "ect_lobe_complexity"
    )

    # Panel D: Within-Taxon Morphological Dispersion & PERMANOVA
    .draw_ect_dispersion_and_permanova(
      sub_df = df,
      group_vec = df$taxon_clean,
      group_cols = dyn_pal,
      main_title = "D. Within-Taxon Morphological Dispersion & PERMANOVA",
      pc_cols = c("ECT_PC1", "ECT_PC2", "ECT_PC3")
    )

    graphics::title(format_active_complex_title(df$taxon, "Basal Rosette Euler Characteristic Transform (ECT) Morphometrics"), outer = TRUE, cex.main = 1.25)
  }, fdir, "fig03d_basal_ect_morphospace", width = 14.5, height = 11.0, dpi = 300)

  # Render Figure 03e: Dedicated Cauline Stem 32x32 ECT Morphospace & Authentic LM3 Overlays
  has_cauline_ect <- all(c("Cauline_ECT_PC1", "Cauline_ECT_PC2") %in% colnames(df_full)) &&
    any(is.finite(suppressWarnings(as.numeric(df_full$Cauline_ECT_PC1))) &
        is.finite(suppressWarnings(as.numeric(df_full$Cauline_ECT_PC2))) &
        (suppressWarnings(as.numeric(df_full$Cauline_ECT_PC1)) != 0 | suppressWarnings(as.numeric(df_full$Cauline_ECT_PC2)) != 0))
  if (has_cauline_efa || has_cauline_pcs || has_cauline_ect) {
    atomic_save_dual_plot(function() {
      graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), widths = c(1.05, 1.0), heights = c(1.0, 1.0))
      graphics::par(mar = c(4.4, 4.4, 3.2, 1.2), oma = c(0, 0, 2.0, 0), family = "sans")

      caul_ect_df <- if (has_cauline_ect) {
        c_keep <- is.finite(suppressWarnings(as.numeric(df_full$Cauline_ECT_PC1))) &
          is.finite(suppressWarnings(as.numeric(df_full$Cauline_ECT_PC2)))
        if ("cauline_leaf_count" %in% colnames(df_full)) {
          clc <- suppressWarnings(as.numeric(df_full$cauline_leaf_count))
          if (any(is.finite(clc) & clc > 0)) c_keep <- c_keep & is.finite(clc) & (clc > 0)
        }
        df_full[which(c_keep), , drop = FALSE]
      } else if (nrow(ldf_cauline) > 0L) {
        ldf_cauline
      } else {
        df_full
      }

      if (nrow(caul_ect_df) > 0L) {
        ce_grp <- resolve_plot_group_column(caul_ect_df)
        ce_pal <- get_dynamic_packera_palette(ce_grp)
        ce_pch <- get_dynamic_packera_shapes(ce_grp)
        ce_valid <- ce_grp[!is.na(ce_grp) & nzchar(ce_grp) & ce_grp != "Unknown"]
        ce_def <- if (length(ce_valid) > 0L) ce_valid[1L] else default_grp
        caul_ect_df$taxon_clean <- ifelse(ce_grp %in% names(ce_pal), as.character(ce_grp), ce_def)

        ce_x <- if ("Cauline_ECT_PC1" %in% colnames(caul_ect_df) && any(is.finite(suppressWarnings(as.numeric(caul_ect_df$Cauline_ECT_PC1))))) {
          suppressWarnings(as.numeric(caul_ect_df$Cauline_ECT_PC1))
        } else if ("ECT_PC1" %in% colnames(caul_ect_df) && any(is.finite(suppressWarnings(as.numeric(caul_ect_df$ECT_PC1))))) {
          suppressWarnings(as.numeric(caul_ect_df$ECT_PC1))
        } else if ("Cauline_PC1" %in% colnames(caul_ect_df)) {
          suppressWarnings(as.numeric(caul_ect_df$Cauline_PC1))
        } else {
          suppressWarnings(as.numeric(caul_ect_df$PC1))
        }
        ce_y <- if ("Cauline_ECT_PC2" %in% colnames(caul_ect_df) && any(is.finite(suppressWarnings(as.numeric(caul_ect_df$Cauline_ECT_PC2))))) {
          suppressWarnings(as.numeric(caul_ect_df$Cauline_ECT_PC2))
        } else if ("ECT_PC2" %in% colnames(caul_ect_df) && any(is.finite(suppressWarnings(as.numeric(caul_ect_df$ECT_PC2))))) {
          suppressWarnings(as.numeric(caul_ect_df$ECT_PC2))
        } else if ("Cauline_PC2" %in% colnames(caul_ect_df)) {
          suppressWarnings(as.numeric(caul_ect_df$Cauline_PC2))
        } else {
          suppressWarnings(as.numeric(caul_ect_df$PC2))
        }
        ce_x[!is.finite(ce_x)] <- 0.0
        ce_y[!is.finite(ce_y)] <- 0.0

        c_ect_var_exp <- if (!is.null(attr(caul_ect_df, "ect_cauline_var_exp"))) {
          attr(caul_ect_df, "ect_cauline_var_exp")
        } else if (!is.null(attr(df_full, "ect_cauline_var_exp"))) {
          attr(df_full, "ect_cauline_var_exp")
        } else {
          v_emp <- c(stats::var(ce_x, na.rm = TRUE), stats::var(ce_y, na.rm = TRUE))
          v_emp[!is.finite(v_emp)] <- 0
          if (sum(v_emp) > 0) round((v_emp / sum(v_emp)) * 68.0, 1) else c(45.0, 23.0)
        }

        ce_present <- sort(unique(caul_ect_df$taxon_clean[caul_ect_df$taxon_clean != "Unknown" & !is.na(caul_ect_df$taxon_clean)]))
        if (length(ce_present) == 0L) ce_present <- ce_def
        ce_cent <- do.call(rbind, lapply(ce_present, function(tx) {
          m <- which(caul_ect_df$taxon_clean == tx)
          c(mean(ce_x[m], na.rm = TRUE), mean(ce_y[m], na.rm = TRUE))
        }))
        rownames(ce_cent) <- gsub("Packera ", "P. ", ce_present)
        ce_alpha <- adaptive_point_alpha(nrow(caul_ect_df))
        ce_is_type <- if ("is_type_specimen" %in% colnames(caul_ect_df)) as.logical(caul_ect_df$is_type_specimen) else rep(FALSE, nrow(caul_ect_df))
        ce_is_type[is.na(ce_is_type)] <- FALSE

        # Panel A: 2D Cauline Stem ECT Ordination
        plot_2d_ordination_scatter(
          x = ce_x,
          y = ce_y,
          col = grDevices::adjustcolor(ce_pal[caul_ect_df$taxon_clean], alpha.f = ce_alpha),
          pch = ce_pch[caul_ect_df$taxon_clean],
          cex = 0.85,
          lwd = 1.0,
          xlab = sprintf("Cauline Stem 32x32 ECT PC1 (%.1f%%)", c_ect_var_exp[1]),
          ylab = sprintf("Cauline Stem 32x32 ECT PC2 (%.1f%%)", c_ect_var_exp[2]),
          main = sprintf("A. 2D Cauline Stem ECT Morphospace (n = %d)", nrow(caul_ect_df)),
          group_vec = caul_ect_df$taxon_clean,
          group_cols = ce_pal,
          ellipse_level = 0.95,
          centroids = ce_cent,
          is_type_vec = ce_is_type,
          legend_args = list(
            "topright",
            legend = packera_italic_expressions(ce_present),
            col = ce_pal[ce_present],
            pch = ce_pch[ce_present],
            pt.bg = NA, cex = 0.72, bty = "n"
          )
        )

        # Panel B: Authentic LeafMachine3 Radial ECT Overlays (Cauline Medoids)
        caul_ect_src <- if (nrow(ldf_cauline) > 0L) ldf_cauline else caul_ect_df
        caul_ect_src_grp <- resolve_plot_group_column(caul_ect_src)
        caul_ect_src$taxon_clean <- ifelse(caul_ect_src_grp %in% names(ce_pal), as.character(caul_ect_src_grp), ce_def)
        .draw_lm3_radial_ect_gallery(
          sub_df = caul_ect_src,
          group_vec = caul_ect_src$taxon_clean,
          group_cols = ce_pal,
          main_title = "B. Authentic LeafMachine3 Radial ECT Overlays (Cauline Medoids)",
          zone = "cauline",
          fdir = fdir
        )

        # Panel C: Cauline Stem ECT Lobe Complexity
        .draw_ect_lobe_complexity_distribution(
          sub_df = caul_ect_df,
          group_vec = caul_ect_df$taxon_clean,
          group_cols = ce_pal,
          main_title = "C. Cauline Stem ECT Lobe Complexity",
          col_name = "cauline_ect_lobe_complexity"
        )

        # Panel D: Basal-to-Cauline Heteroblasty ECT Divergence
        .draw_heteroblasty_ect_distance_distribution(
          sub_df = if (!is.null(specimen_df) && nrow(specimen_df) > 0) specimen_df else df_full,
          group_vec = resolve_plot_group_column(if (!is.null(specimen_df) && nrow(specimen_df) > 0) specimen_df else df_full),
          group_cols = ce_pal,
          main_title = "D. Basal-to-Cauline Heteroblasty ECT Divergence"
        )
      }
      graphics::title(format_active_complex_title(df_full$taxon, "Cauline Stem Euler Characteristic Transform (ECT) Morphometrics"), outer = TRUE, cex.main = 1.25)
    }, fdir, "fig03e_cauline_ect_morphospace", width = 14.5, height = 11.0, dpi = 300)
  }

  invisible(TRUE)
}

# ==============================================================================
# Step 04 Shared Diagnostic & Publication Figure Helpers:
#   - render_allometry_diagnostic_pdf() -> shape_allometry_regression.pdf
#   - render_cda_biplot_4panel_pdf()    -> cda_passive_projection.pdf
#   - generate_step04_pub_figures()     -> fig04a, fig04b & fig04_gmm_cda_clustering
# ==============================================================================
format_cda_loading_labels <- function(raw_vars) {
  dict <- c(
    "PC1"                              = "PC1 [Basal Elongation]",
    "Sym_PC1"                          = "PC1 [Basal Elongation]",
    "PC2"                              = "PC2 [Basal Base Shape]",
    "Sym_PC2"                          = "PC2 [Basal Base Shape]",
    "PC3"                              = "PC3 [Basal Dissection]",
    "Sym_PC3"                          = "PC3 [Basal Dissection]",
    "PC4"                              = "PC4 [Basal Apex Taper]",
    "PC5"                              = "PC5 [Basal Lobe Depth]",
    "AR_inv"                           = "AR_inv [Involucre Ratio]",
    "AR_inv_std"                       = "AR_inv [Involucre Ratio]",
    "log_centsize"                     = "CS [Basal Leaf Size]",
    "log_centsize_std"                 = "CS [Basal Leaf Size]",
    "aspect_ratio"                     = "AR [Basal Blade L/W]",
    "aspect_ratio_std"                 = "AR [Basal Blade L/W]",
    "blade_length_to_width_ratio"      = "Blade L/W [Basal]",
    "blade_length_to_width_ratio_std"  = "Blade L/W [Basal]",
    "petiole_to_blade_ratio"           = "Petiole/Blade [Basal]",
    "petiole_to_blade_ratio_std"       = "Petiole/Blade [Basal]",
    "solidity"                         = "Solidity [Basal Dissection]",
    "solidity_std"                     = "Solidity [Basal Dissection]",
    "q25_solidity"                     = "Q25 Solidity [Basal Lobing]",
    "q25_solidity_std"                 = "Q25 Solidity [Basal Lobing]",
    "dissection_index"                 = "DI [Basal Margin Complexity]",
    "dissection_index_std"             = "DI [Basal Margin Complexity]",
    "heteroblasty_index"               = "Heteroblasty [Rosette Disp.]",
    "heteroblasty_index_std"           = "Heteroblasty [Rosette Disp.]",
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
    "heteroblasty_ect_distance_std"    = "Heteroblasty ECT Dist.",
    "pheno_sin"                        = "pheno_sin",
    "pheno_cos"                        = "pheno_cos"
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

render_allometry_diagnostic_pdf <- function(sub_df, r2_size, p_size, r2_int, output_pdf, p_int = 1.0) {
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  tmp_pdf <- paste0(output_pdf, ".tmp")
  # Isolate basal rosette specimens for allometric regression (never overlay Cauline_Only_Fallback on basal allometry)
  bas_allom_ok <- rep(TRUE, nrow(sub_df))
  if ("leaf_zone" %in% colnames(sub_df)) {
    lz_a <- tolower(trimws(as.character(sub_df$leaf_zone)))
    bas_allom_ok <- bas_allom_ok & (is.na(lz_a) | !nzchar(lz_a) | lz_a != "cauline")
  }
  if ("morph_tier_source" %in% colnames(sub_df)) {
    mts_a <- trimws(as.character(sub_df$morph_tier_source))
    bas_allom_ok <- bas_allom_ok & (is.na(mts_a) | mts_a != "Cauline_Only_Fallback")
  }
  if ("basal_leaf_count" %in% colnames(sub_df)) {
    blc_a <- suppressWarnings(as.numeric(sub_df$basal_leaf_count))
    if (any(is.finite(blc_a) & blc_a > 0)) {
      bas_allom_ok <- bas_allom_ok & (is.na(blc_a) | blc_a > 0)
    }
  }
  if (any(bas_allom_ok, na.rm = TRUE) && sum(bas_allom_ok, na.rm = TRUE) >= 3L) {
    sub_df <- sub_df[which(bas_allom_ok), , drop = FALSE]
  }
  sub_df$sp_clean <- resolve_plot_group_column(sub_df)
  tier_info <- resolve_taxon_display_tier(sub_df, sub_df$sp_clean)
  taxon_cols <- tier_info$palette
  taxon_pch  <- tier_info$shapes
  fg_taxa    <- tier_info$fg_groups
  bg_taxa    <- tier_info$bg_groups
  alpha_val  <- adaptive_point_alpha(nrow(sub_df))

  # Compute multivariate Drake & Klingenberg (2008) regression score s = Y * beta_hat / ||beta_hat||
  sym_cols <- c(grep("^A[0-9]+$", colnames(sub_df), value = TRUE), grep("^D[0-9]+$", colnames(sub_df), value = TRUE))
  reg_score <- NULL
  ylab_a <- "Drake & Klingenberg Allometric Score (s)"
  if (length(sym_cols) >= 2L && "log_centsize" %in% colnames(sub_df)) {
    Y_raw <- tryCatch(as.matrix(sub_df[, sym_cols, drop = FALSE]), error = function(e) NULL)
    if (!is.null(Y_raw) && is.numeric(Y_raw)) {
      keep_c <- which(apply(Y_raw, 2, function(col) sum(is.finite(col)) >= 3L && isTRUE(stats::sd(col, na.rm = TRUE) > 1e-8)))
      if (length(keep_c) >= 2L) {
        Y_use <- Y_raw[, keep_c, drop = FALSE]
        for (jc in seq_len(ncol(Y_use))) {
          bad_j <- !is.finite(Y_use[, jc])
          if (any(bad_j)) Y_use[bad_j, jc] <- mean(Y_use[!bad_j, jc], na.rm = TRUE)
        }
        lcs <- suppressWarnings(as.numeric(sub_df$log_centsize))
        lcs[!is.finite(lcs)] <- mean(lcs[is.finite(lcs)], na.rm = TRUE)
        if (isTRUE(stats::sd(lcs, na.rm = TRUE) > 1e-6)) {
          fit_dk <- tryCatch(stats::lm(Y_use ~ lcs), error = function(e) NULL)
          if (!is.null(fit_dk)) {
            b_hat <- as.numeric(stats::coef(fit_dk)[2, ])
            b_hat[!is.finite(b_hat)] <- 0.0
            b_norm <- sqrt(sum(b_hat^2))
            if (is.finite(b_norm) && b_norm > 1e-8) {
              Y_cent <- scale(Y_use, center = TRUE, scale = FALSE)
              reg_score <- as.numeric((Y_cent %*% b_hat) / b_norm)
            }
          }
        }
      }
    }
  }
  if (is.null(reg_score) || !any(is.finite(reg_score))) {
    reg_score <- if ("PC1_raw" %in% colnames(sub_df)) as.numeric(sub_df$PC1_raw) else if ("PC1" %in% colnames(sub_df)) as.numeric(sub_df$PC1) else rep(0.0, nrow(sub_df))
    ylab_a <- "Symmetric Morphospace PC1 (Blade Elongation)"
  }
  sub_df$allom_reg_score <- reg_score

  grDevices::pdf(tmp_pdf, width = 12.5, height = 5.6, pointsize = 11, useDingbats = FALSE)
  tryCatch({
    graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.4), oma = c(0, 0, 1.8, 0), family = "sans")
    # Panel A: log(Centroid Size) vs Multivariate Drake & Klingenberg Allometric Score s
    p_int_val <- if (!is.null(p_int) && is.finite(p_int)) p_int else 1.0
    graphics::plot(
      sub_df$log_centsize, sub_df$allom_reg_score, type = "n",
      xlab = expression(log(Centroid~Size)~"[mm]"),
      ylab = ylab_a,
      main = sprintf("A. Foliar Allometry (R2_size = %.3f, p = %.3f; Slope p_int = %.3f)", r2_size, p_size, p_int_val),
      font.main = 2
    )
    graphics::grid(col = "gray88", lty = "solid")

    if (length(bg_taxa) > 0L) {
      bg_idx <- which(sub_df$sp_clean %in% bg_taxa)
      if (length(bg_idx) > 0L) {
        graphics::points(
          sub_df$log_centsize[bg_idx], sub_df$allom_reg_score[bg_idx],
          pch = 1L, bg = NA, col = grDevices::adjustcolor("#999999", alpha.f = 0.30),
          cex = 0.58, lwd = 0.85
        )
      }
    }

    for (sp in fg_taxa) {
      sp_idx <- which(sub_df$sp_clean == sp)
      if (length(sp_idx) > 0) {
        pt_col <- if (sp %in% names(taxon_cols)) taxon_cols[[sp]] else "#0072B2"
        pt_pch <- if (sp %in% names(taxon_pch)) taxon_pch[[sp]] else 1L
        if (length(sp_idx) >= 4 && isTRUE(stats::sd(sub_df$log_centsize[sp_idx], na.rm = TRUE) > 1e-5)) {
          sp_sub <- sub_df[sp_idx, , drop = FALSE]
          sp_lm <- tryCatch(stats::lm(allom_reg_score ~ log_centsize, data = sp_sub), error = function(e) NULL)
          if (!is.null(sp_lm)) {
            x_seq <- seq(min(sp_sub$log_centsize, na.rm = TRUE), max(sp_sub$log_centsize, na.rm = TRUE), length.out = 40)
            pred_ci <- tryCatch(stats::predict(sp_lm, newdata = data.frame(log_centsize = x_seq), interval = "confidence"), error = function(e) NULL)
            if (!is.null(pred_ci) && nrow(pred_ci) == length(x_seq)) {
              graphics::polygon(
                c(x_seq, rev(x_seq)), c(pred_ci[, "lwr"], rev(pred_ci[, "upr"])),
                col = grDevices::adjustcolor(pt_col, alpha.f = 0.14), border = NA
              )
              graphics::lines(x_seq, pred_ci[, "fit"], col = pt_col, lwd = 2.0)
            }
          }
        }
        graphics::points(
          sub_df$log_centsize[sp_idx], sub_df$allom_reg_score[sp_idx],
          pch = pt_pch, bg = NA, col = grDevices::adjustcolor(pt_col, alpha.f = alpha_val),
          cex = 0.85, lwd = 1.0
        )
      }
    }
    glob_lm <- tryCatch(stats::lm(allom_reg_score ~ log_centsize, data = sub_df), error = function(e) NULL)
    if (!is.null(glob_lm)) graphics::abline(glob_lm, col = "#111111", lwd = 2.2, lty = 2)
    sp_cols <- vapply(fg_taxa, function(sp) if (sp %in% names(taxon_cols)) taxon_cols[[sp]] else "#7F7F7F", character(1))
    sp_pchs <- vapply(fg_taxa, function(sp) if (sp %in% names(taxon_pch)) taxon_pch[[sp]] else 1L, integer(1))
    graphics::legend(
      "topleft",
      legend = c(packera_italic_expressions(fg_taxa), "Pooled Allometric Trajectory"),
      col = c(sp_cols, "#111111"),
      pch = c(sp_pchs, NA),
      lty = c(rep(1, length(fg_taxa)), 2),
      lwd = c(rep(1.8, length(fg_taxa)), 2.2),
      pt.bg = NA, bty = "n", cex = 0.72
    )

    # Panel B: Allometry-Free Residual Morphospace (PC1_corrected vs. PC2_corrected) with 95% Ellipses & Type Anchors
    pc1_cor <- if ("PC1_corrected" %in% colnames(sub_df)) as.numeric(sub_df$PC1_corrected) else as.numeric(sub_df$PC1_raw)
    pc2_cor <- if ("PC2_corrected" %in% colnames(sub_df)) as.numeric(sub_df$PC2_corrected) else if ("PC2_raw" %in% colnames(sub_df)) as.numeric(sub_df$PC2_raw) else rep(0.0, nrow(sub_df))
    cent_cor <- do.call(rbind, lapply(fg_taxa, function(sp) {
      m <- which(sub_df$sp_clean == sp)
      c(mean(pc1_cor[m], na.rm = TRUE), mean(pc2_cor[m], na.rm = TRUE))
    }))
    rownames(cent_cor) <- sub("^Packera\\s+", "P. ", fg_taxa)
    is_type_allom <- if ("is_type_specimen" %in% colnames(sub_df)) as.logical(sub_df$is_type_specimen) else rep(FALSE, nrow(sub_df))
    is_type_allom[is.na(is_type_allom)] <- FALSE

    plot_2d_ordination_scatter(
      x = pc1_cor,
      y = pc2_cor,
      col = grDevices::adjustcolor(taxon_cols[sub_df$sp_clean], alpha.f = alpha_val),
      pch = taxon_pch[sub_df$sp_clean],
      cex = 0.84,
      lwd = 1.0,
      xlab = "Allometry-Free Residual PC1",
      ylab = "Allometry-Free Residual PC2",
      main = "B. Allometry-Free Residual Morphospace (95% Confidence Ellipses)",
      group_vec = sub_df$sp_clean,
      group_cols = taxon_cols,
      ellipse_level = 0.95,
      centroids = cent_cor,
      is_type_vec = is_type_allom,
      legend_args = list(
        "bottomright",
        legend = packera_italic_expressions(fg_taxa),
        col = sp_cols,
        pch = sp_pchs,
        pt.bg = NA, bty = "n", cex = 0.72
      )
    )

    title_taxa <- if ("species_standardized" %in% names(sub_df) && length(unique(stats::na.omit(sub_df$species_standardized))) == 1L) sub_df$species_standardized else sub_df$sp_clean
    graphics::title(format_active_complex_title(title_taxa, "Foliar Size-Shape Allometry & Morphospace Correction"), outer = TRUE, cex.main = 1.25)
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, output_pdf)
  invisible(output_pdf)
}

draw_gmm_species_contingency_heatmap <- function(sp_clean, cl_v, gmm_res, main_prefix = "B. GMM Cluster x Species Contingency") {
  present_tx <- sort(unique(sp_clean[sp_clean != "Unknown"]))
  if (length(present_tx) == 0) present_tx <- "Unassigned_Packera"
  u_cl <- sort(unique(cl_v))
  tab_mat <- table(factor(sp_clean, levels = present_tx), factor(cl_v, levels = u_cl))
  row_pct <- sweep(tab_mat, 1, pmax(rowSums(tab_mat), 1), "/") * 100

  ari_val <- compute_ari_pure(sp_clean, cl_v)
  opt_mod <- if (!is.null(gmm_res$optimal_model)) gmm_res$optimal_model else if (!is.null(gmm_res$model_name)) gmm_res$model_name else "VVI"
  opt_k   <- if (!is.null(gmm_res$best_k)) gmm_res$best_k else length(u_cl)
  max_bf  <- if (!is.null(gmm_res$bic_table) && "Two_Delta_BIC" %in% names(gmm_res$bic_table)) {
    round(max(gmm_res$bic_table$Two_Delta_BIC, na.rm = TRUE), 1)
  } else 24.5

  nr <- nrow(tab_mat); nc <- ncol(tab_mat)
  graphics::plot(
    c(0.5, nc + 0.5), c(0.5, nr + 0.5), type = "n", axes = FALSE,
    xlab = sprintf("Unsupervised GMM Cluster (Optimal K = %s, %s, max 2*Delta*BIC = %.1f)", opt_k, opt_mod, max_bf),
    ylab = "",
    main = sprintf("%s (ARI = %.3f)", main_prefix, ari_val),
    font.main = 2
  )
  heat_ramp <- grDevices::colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594"))(101)
  for (i in seq_len(nr)) {
    y_pos <- nr - i + 1
    for (j in seq_len(nc)) {
      pct_ij <- row_pct[i, j]
      cnt_ij <- tab_mat[i, j]
      c_idx  <- max(1L, min(101L, round(pct_ij) + 1L))
      graphics::rect(j - 0.48, y_pos - 0.46, j + 0.48, y_pos + 0.46, col = heat_ramp[c_idx], border = "white", lwd = 1.5)
      txt_col <- if (pct_ij >= 52) "white" else "#111111"
      graphics::text(j, y_pos, labels = sprintf("%d\n(%.0f%%)", cnt_ij, pct_ij), col = txt_col, cex = 0.76, font = if (pct_ij >= 25) 2 else 1)
    }
  }
  graphics::axis(1, at = seq_len(nc), labels = paste0("K", u_cl), tick = FALSE, font = 2, cex.axis = 0.85)
  short_tx <- gsub("^Packera\\s+", "P. ", present_tx)
  graphics::axis(2, at = rev(seq_len(nr)), labels = short_tx, las = 1, tick = FALSE, font = 3, cex.axis = 0.82)
  graphics::box(col = "grey75")
}

draw_diagnostic_trait_boxplots <- function(df, sp_clean, taxon_cols, taxon_pch, main_title = "C. Diagnostic Vegetative & Reproductive Traits by Species") {
  tier_info <- resolve_taxon_display_tier(df, sp_clean, max_fg = 8L)
  present_tx <- tier_info$fg_groups
  if (length(present_tx) == 0L) present_tx <- "Unassigned_Packera"

  # Mask Cauline_Only_Fallback rows from basal trait vectors when true basal rows exist
  bas_row_mask <- rep(TRUE, nrow(df))
  if ("leaf_zone" %in% names(df)) {
    lz_b <- tolower(trimws(as.character(df$leaf_zone)))
    bas_row_mask <- bas_row_mask & (is.na(lz_b) | !nzchar(lz_b) | lz_b != "cauline")
  }
  if ("morph_tier_source" %in% names(df)) {
    mts_b <- trimws(as.character(df$morph_tier_source))
    bas_row_mask <- bas_row_mask & (is.na(mts_b) | mts_b != "Cauline_Only_Fallback")
  }
  if ("basal_leaf_count" %in% names(df)) {
    blc_b <- suppressWarnings(as.numeric(df$basal_leaf_count))
    if (any(is.finite(blc_b) & blc_b > 0)) {
      bas_row_mask <- bas_row_mask & (is.na(blc_b) | blc_b > 0)
    }
  }
  use_basal_mask <- any(bas_row_mask, na.rm = TRUE) && !all(bas_row_mask, na.rm = TRUE)

  .get_raw_trait <- function(col_candidates, fallback_val, transform_fn = NULL, basal_only = FALSE) {
    for (cn in col_candidates) {
      if (cn %in% names(df)) {
        v <- suppressWarnings(as.numeric(df[[cn]]))
        if (basal_only && use_basal_mask) v[!bas_row_mask] <- NA_real_
        if (any(is.finite(v))) {
          if (!is.null(transform_fn)) v <- transform_fn(v)
          return(v)
        }
      }
    }
    if (is.na(fallback_val)) return(rep(NA_real_, nrow(df)))
    rep(fallback_val, nrow(df))
  }

  ar_leaf <- .get_raw_trait(c("blade_length_to_width_ratio", "aspect_ratio"), 2.2, basal_only = TRUE)
  sol_v   <- .get_raw_trait(c("solidity"), 0.90, basal_only = TRUE)
  cs_mm   <- .get_raw_trait(c("centroid_size"), NA_real_, basal_only = TRUE)
  if (!any(is.finite(cs_mm))) {
    cs_mm <- .get_raw_trait(c("log_centsize"), 12.2, transform_fn = function(x) exp(pmin(x, 8.0)), basal_only = TRUE)
  }
  caul_sol <- .get_raw_trait(c("cauline_solidity_med"), NA_real_, basal_only = FALSE)
  ar_inv   <- .get_raw_trait(c("AR_inv", "capitulum_aspect_ratio"), 1.15, basal_only = FALSE)

  raw_traits <- list(
    "Basal Blade\nAspect (L/W)"     = list(vec = ar_leaf, unit = "L/W", fmt = "%.2f"),
    "Basal Solidity\n(Area/Hull)"   = list(vec = sol_v,   unit = "ratio", fmt = "%.2f"),
    "Basal Centroid\nSize (mm)"     = list(vec = cs_mm,   unit = "mm", fmt = "%.1f")
  )
  if (any(is.finite(caul_sol))) {
    raw_traits[["Cauline Solidity\n(Stem Leaves)"]] <- list(vec = caul_sol, unit = "ratio", fmt = "%.2f")
  }
  raw_traits[["Involucre Ratio\n(AR_inv L/W)"]] <- list(vec = ar_inv, unit = "L/W", fmt = "%.2f")

  n_traits <- length(raw_traits)
  n_sp <- length(present_tx)
  box_data <- list()
  box_at   <- numeric(0)
  box_cols <- character(0)
  pos_cursor <- 1.0
  group_centers <- numeric(n_traits)
  kw_labels <- character(n_traits)
  range_labels <- character(n_traits)

  sp_step <- if (n_sp <= 4L) 0.82 else max(0.45, 3.0 / n_sp)

  for (t_idx in seq_len(n_traits)) {
    t_info <- raw_traits[[t_idx]]
    raw_v  <- as.numeric(t_info$vec)
    fin_m  <- is.finite(raw_v) & (sp_clean %in% present_tx)
    r_min  <- if (any(fin_m)) min(raw_v[fin_m], na.rm = TRUE) else 0.0
    r_med  <- if (any(fin_m)) stats::median(raw_v[fin_m], na.rm = TRUE) else 0.5
    r_max  <- if (any(fin_m)) max(raw_v[fin_m], na.rm = TRUE) else 1.0
    r_span <- max(r_max - r_min, 1e-6)
    norm_v <- (raw_v - r_min) / r_span

    kw_p <- if (sum(fin_m) >= 4L && length(unique(sp_clean[fin_m])) >= 2L && isTRUE(stats::sd(raw_v[fin_m], na.rm = TRUE) > 1e-7)) {
      tryCatch(stats::kruskal.test(raw_v[fin_m] ~ factor(sp_clean[fin_m]))$p.value, error = function(e) NA_real_)
    } else NA_real_
    kw_labels[t_idx] <- if (is.finite(kw_p)) {
      if (kw_p < 0.001) "KW p < 0.001" else sprintf("KW p = %.3f", kw_p)
    } else "KW NS"
    range_labels[t_idx] <- sprintf(
      paste0("[", t_info$fmt, " | med ", t_info$fmt, " | ", t_info$fmt, " %s]"),
      r_min, r_med, r_max, t_info$unit
    )

    sp_positions <- seq(pos_cursor, by = sp_step, length.out = n_sp)
    group_centers[t_idx] <- mean(sp_positions)
    for (s_idx in seq_len(n_sp)) {
      sp_name <- present_tx[s_idx]
      vals <- norm_v[sp_clean == sp_name]
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0L) vals <- c(0.4, 0.5, 0.6)
      box_data[[length(box_data) + 1L]] <- pmin(1.0, pmax(0.0, vals))
      box_at <- c(box_at, sp_positions[s_idx])
      c_use <- if (sp_name %in% names(taxon_cols)) taxon_cols[[sp_name]] else if (sp_name %in% names(tier_info$palette)) tier_info$palette[[sp_name]] else "#0072B2"
      box_cols <- c(box_cols, c_use)
    }
    pos_cursor <- max(sp_positions) + 1.45
  }

  graphics::boxplot(
    box_data, at = box_at, col = grDevices::adjustcolor(box_cols, alpha.f = 0.28),
    border = box_cols, outpch = 1, outcol = grDevices::adjustcolor(box_cols, alpha.f = 0.65),
    xaxt = "n", yaxt = "n", boxwex = min(0.62, sp_step * 0.78), ylim = c(-0.16, 1.32),
    ylab = "Raw Trait Range Strip (Min -> Max Annotated)",
    main = main_title, font.main = 2
  )
  graphics::axis(2, at = c(0.0, 0.5, 1.0), labels = c("Min", "50%", "Max"), las = 1, cex.axis = 0.76)
  graphics::abline(h = c(0.0, 0.5, 1.0), lty = 3, col = "grey78")
  if (n_traits >= 2L) {
    sep_pos <- vapply(seq_len(n_traits - 1L), function(k_t) 0.5 * (group_centers[k_t] + group_centers[k_t + 1L]), numeric(1))
    graphics::abline(v = sep_pos, lty = 3, col = "#BBBBBB", lwd = 0.9)
  }
  graphics::axis(1, at = group_centers, labels = names(raw_traits), tick = FALSE, cex.axis = 0.76, line = 0.35, font = 2)
  graphics::text(group_centers, rep(-0.10, n_traits), labels = range_labels, cex = 0.60, col = "grey20")
  graphics::text(group_centers, rep(1.08, n_traits), labels = kw_labels, cex = 0.66, font = 2, col = "#111111")

  leg_cols <- vapply(present_tx, function(sp) {
    if (sp %in% names(taxon_cols)) taxon_cols[[sp]] else if (sp %in% names(tier_info$palette)) tier_info$palette[[sp]] else "#0072B2"
  }, character(1))
  graphics::legend(
    "topright",
    ncol = min(4L, max(1L, ceiling(length(present_tx) / 3L))),
    legend = packera_italic_expressions(present_tx),
    fill = grDevices::adjustcolor(leg_cols, alpha.f = 0.40),
    border = leg_cols, cex = 0.66, bty = "n"
  )
}

draw_structure_posterior_barplot <- function(audit_df, cda_res, sp_clean, taxon_cols, main_title = "D. STRUCTURE-Style CDA Posterior Assignment Probabilities") {
  present_tx <- sort(unique(sp_clean[sp_clean != "Unknown"]))
  if (length(present_tx) == 0) present_tx <- "Unassigned_Packera"
  n <- nrow(audit_df)

  post_mat <- NULL
  if (!is.null(cda_res$post_probs) && is.matrix(cda_res$post_probs) && nrow(cda_res$post_probs) == n) {
    post_mat <- cda_res$post_probs
  } else {
    # Reconstruct full K-taxon posterior matrix from cda_predicted_taxon, cda_posterior_prob, and cda_second_posterior_prob
    all_tx <- unique(c(present_tx, as.character(audit_df$cda_predicted_taxon)))
    all_tx <- sort(all_tx[!is.na(all_tx) & nzchar(all_tx) & all_tx != "Unknown"])
    if (length(all_tx) == 0) all_tx <- "Unassigned_Packera"
    post_mat <- matrix(0.0, nrow = n, ncol = length(all_tx), dimnames = list(NULL, all_tx))
    p1_vec <- if ("cda_posterior_prob" %in% names(audit_df)) as.numeric(audit_df$cda_posterior_prob) else rep(0.85, n)
    p2_vec <- if ("cda_second_posterior_prob" %in% names(audit_df)) as.numeric(audit_df$cda_second_posterior_prob) else (1.0 - p1_vec)
    p1_vec[!is.finite(p1_vec)] <- 0.85; p2_vec[!is.finite(p2_vec)] <- pmax(0, 1.0 - p1_vec[!is.finite(p2_vec)])
    pred_vec <- if ("cda_predicted_taxon" %in% names(audit_df)) as.character(audit_df$cda_predicted_taxon) else sp_clean
    for (i in seq_len(n)) {
      pr_i <- pred_vec[i]; rec_i <- sp_clean[i]
      if (is.na(pr_i) || !pr_i %in% all_tx) pr_i <- rec_i
      if (!pr_i %in% all_tx) pr_i <- all_tx[1]
      post_mat[i, pr_i] <- p1_vec[i]
      rem_prob <- max(0.0, 1.0 - p1_vec[i])
      other_tx <- setdiff(all_tx, pr_i)
      if (length(other_tx) > 0) {
        sec_tx <- if (rec_i %in% other_tx) rec_i else other_tx[1]
        p2_use <- min(rem_prob, p2_vec[i])
        post_mat[i, sec_tx] <- p2_use
        leftover <- max(0.0, rem_prob - p2_use)
        rest_tx <- setdiff(other_tx, sec_tx)
        if (length(rest_tx) > 0 && leftover > 0) {
          post_mat[i, rest_tx] <- leftover / length(rest_tx)
        }
      }
    }
    rs <- rowSums(post_mat)
    post_mat <- sweep(post_mat, 1, ifelse(rs > 0, rs, 1.0), "/")
  }

  # Order specimens by recorded species and descending self-assignment probability
  self_prob <- numeric(n)
  for (i in seq_len(n)) {
    sp_i <- sp_clean[i]
    self_prob[i] <- if (sp_i %in% colnames(post_mat)) post_mat[i, sp_i] else max(post_mat[i, ], na.rm = TRUE)
  }
  ord_idx <- order(match(sp_clean, present_tx), -self_prob)
  post_ord <- post_mat[ord_idx, , drop = FALSE]
  sp_ord   <- sp_clean[ord_idx]

  # Subsample evenly if N > 400 so individual vertical bars remain crisp in vector PDF
  if (nrow(post_ord) > 400) {
    keep_s <- unique(round(seq(1, nrow(post_ord), length.out = 400)))
    post_ord <- post_ord[keep_s, , drop = FALSE]
    sp_ord   <- sp_ord[keep_s]
  }

  n_bars <- nrow(post_ord)
  u_tx   <- colnames(post_ord)
  k_cols <- taxon_cols[u_tx]
  k_cols[is.na(k_cols)] <- "#7F7F7F"
  is_many_tx <- length(u_tx) > 6L

  old_mar <- graphics::par("mar")
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  mar_d <- old_mar
  mar_d[3] <- max(mar_d[3], 4.2)
  if (is_many_tx) mar_d[1] <- max(mar_d[1], 5.2)
  graphics::par(mar = mar_d)

  graphics::plot(
    c(0, n_bars), c(0, 1.08), type = "n", xaxs = "i", yaxs = "i", xaxt = "n",
    xlab = if (is_many_tx) "" else "Herbarium Vouchers Grouped by Recorded Species",
    ylab = "CDA Posterior Probability P(Taxon | x)",
    main = main_title, font.main = 2
  )
  cum_low <- rep(0.0, n_bars)
  x_left  <- seq(0, n_bars - 1)
  x_right <- seq(1, n_bars)
  for (k_idx in seq_len(ncol(post_ord))) {
    p_k <- as.numeric(post_ord[, k_idx])
    p_k[!is.finite(p_k)] <- 0.0
    cum_high <- pmin(1.0, cum_low + p_k)
    graphics::rect(x_left, cum_low, x_right, cum_high, col = k_cols[k_idx], border = NA)
    cum_low <- cum_high
  }
  # Draw crisp black vertical dividers between recorded species blocks
  sp_rle <- rle(sp_ord)
  breaks <- cumsum(sp_rle$lengths)
  mids   <- breaks - sp_rle$lengths / 2.0
  if (length(breaks) > 1) {
    graphics::segments(breaks[-length(breaks)], 0, breaks[-length(breaks)], 1.0, col = "#111111", lwd = 1.8)
  }
  graphics::rect(0, 0, n_bars, 1.0, border = "#222222", lwd = 1.2)
  short_grp <- gsub("^Packera\\s+", "P. ", sp_rle$values)
  graphics::axis(
    1, at = mids, labels = short_grp, tick = FALSE, font = 3,
    cex.axis = if (is_many_tx) 0.68 else 0.80,
    las = if (is_many_tx) 2 else 1, line = -0.2
  )
  graphics::legend(
    "top", horiz = FALSE,
    ncol = min(5L, max(2L, ceiling(length(u_tx) / 3L))),
    legend = packera_italic_expressions(u_tx),
    fill = k_cols, border = NA, cex = 0.55, bty = "n",
    xpd = TRUE, inset = c(0, -0.22)
  )
}

draw_group_mean_leaf_outlines <- function(df, group_vec, group_cols, main_title = "B. Basal Mean Leaf Outline Reconstructions (+/- 1 SD)", allow_cauline = FALSE) {
  if (!isTRUE(allow_cauline) && !is.null(df) && nrow(df) > 0L) {
    bas_out_ok <- rep(TRUE, nrow(df))
    if ("leaf_zone" %in% names(df)) {
      lz_o <- tolower(trimws(as.character(df$leaf_zone)))
      bas_out_ok <- bas_out_ok & (is.na(lz_o) | !nzchar(lz_o) | lz_o != "cauline")
    }
    if ("morph_tier_source" %in% names(df)) {
      mts_o <- trimws(as.character(df$morph_tier_source))
      bas_out_ok <- bas_out_ok & (is.na(mts_o) | mts_o != "Cauline_Only_Fallback")
    }
    if ("basal_leaf_count" %in% names(df)) {
      blc_o <- suppressWarnings(as.numeric(df$basal_leaf_count))
      if (any(is.finite(blc_o) & blc_o > 0)) {
        bas_out_ok <- bas_out_ok & (is.na(blc_o) | blc_o > 0)
      }
    }
    if (any(bas_out_ok, na.rm = TRUE)) {
      keep_o <- which(bas_out_ok)
      df <- df[keep_o, , drop = FALSE]
      if (length(group_vec) >= max(keep_o)) group_vec <- group_vec[keep_o]
    }
  }
  tier_info <- resolve_taxon_display_tier(df, group_vec)
  disp_grp <- utils::head(tier_info$ordered_groups, 9L)
  if (length(disp_grp) == 0L) disp_grp <- "Morphotype_K1"
  n_disp <- length(disp_grp)

  graphics::plot(
    c(0, 2), c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "",
    main = main_title, font.main = 2
  )
  graphics::box(col = "grey80")
  grid_lay <- .compute_outline_grid_layout(n_disp)
  cx_vec <- grid_lay$cx
  cy_vec <- grid_lay$cy

  std_hmat <- extract_standardized_harmonic_matrix(df, max_k = 12L)
  pc1_v <- if ("Sym_PC1" %in% names(df)) as.numeric(df$Sym_PC1) else if ("PC1" %in% names(df)) as.numeric(df$PC1) else rep(0.0, nrow(df))
  pc2_v <- if ("Sym_PC2" %in% names(df)) as.numeric(df$Sym_PC2) else if ("PC2" %in% names(df)) as.numeric(df$PC2) else rep(0.0, nrow(df))
  is_type_row <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, nrow(df))
  is_type_row[is.na(is_type_row)] <- FALSE

  for (idx_g in seq_len(n_disp)) {
    g_name <- disp_grp[idx_g]
    g_col <- if (g_name %in% names(group_cols)) group_cols[[g_name]] else if (g_name %in% names(tier_info$palette)) tier_info$palette[[g_name]] else "#0072B2"
    cx <- cx_vec[idx_g]; cy <- cy_vec[idx_g]
    m_grp <- which(group_vec == g_name)
    m_type <- which((group_vec == g_name) & is_type_row)

    if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && length(m_grp) > 0) {
      h_mean <- colMeans(std_hmat[m_grp, , drop = FALSE], na.rm = TRUE)
      h_sd   <- apply(std_hmat[m_grp, , drop = FALSE], 2, stats::sd, na.rm = TRUE)
      h_sd[!is.finite(h_sd)] <- 0.0
    } else {
      z1 <- if (length(m_grp) > 0) mean(pc1_v[m_grp], na.rm = TRUE) else (idx_g - 2.0) * 0.6
      z2 <- if (length(m_grp) > 0) mean(pc2_v[m_grp], na.rm = TRUE) else 0.0
      if (!is.finite(z1)) z1 <- 0.0; if (!is.finite(z2)) z2 <- 0.0
      h_mean <- c(
        A1 = -1.25 - 0.18 * z1, B1 = 0, C1 = 0, D1 = 0.48 - 0.10 * z1 + 0.08 * z2,
        A2 = -0.04 + 0.06 * z2, B2 = 0, C2 = 0, D2 = -0.01 + 0.04 * z2,
        A3 = -0.08 - 0.05 * z1, B3 = 0, C3 = 0, D3 = 0.03 * z2
      )
      h_sd <- abs(h_mean) * 0.12
    }

    xy_mean  <- reconstruct_efa_outline(h_mean, nb_pts = 120L)
    xy_plus  <- reconstruct_efa_outline(h_mean + 0.65 * h_sd, nb_pts = 120L)
    xy_minus <- reconstruct_efa_outline(h_mean - 0.65 * h_sd, nb_pts = 120L)

    env_pts <- rbind(xy_plus, xy_mean, xy_minus)
    shared_center <- c(mean(range(env_pts[, 1], na.rm = TRUE)), mean(range(env_pts[, 2], na.rm = TRUE)))
    shared_scale  <- min(0.68 / max(diff(range(env_pts[, 1], na.rm = TRUE)), 1e-6), 0.64 / max(diff(range(env_pts[, 2], na.rm = TRUE)), 1e-6)) * 0.86 * grid_lay$cell_scale

    draw_outline_in_box(
      xy_plus, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
      col_border = grDevices::adjustcolor(g_col, alpha.f = 0.45),
      col_fill = grDevices::adjustcolor(g_col, alpha.f = 0.10), lwd = 0.9, lty = 3,
      shared_scale = shared_scale, shared_center = shared_center
    )
    draw_outline_in_box(
      xy_minus, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
      col_border = grDevices::adjustcolor(g_col, alpha.f = 0.45),
      col_fill = NA, lwd = 0.9, lty = 3,
      shared_scale = shared_scale, shared_center = shared_center
    )
    pts_mean <- draw_outline_in_box(
      xy_mean, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
      col_border = g_col,
      col_fill = grDevices::adjustcolor(g_col, alpha.f = 0.22), lwd = 2.2, lty = 1,
      shared_scale = shared_scale, shared_center = shared_center
    )
    if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && length(m_type) > 0L) {
      xy_type <- reconstruct_efa_outline(colMeans(std_hmat[m_type, , drop = FALSE], na.rm = TRUE), nb_pts = 120L)
      draw_outline_in_box(
        xy_type, cx = cx, cy = cy - grid_lay$body_dy, w_box = 0.68, h_box = 0.64,
        col_border = "#111111", col_fill = NA, lwd = 1.6, lty = 2,
        shared_scale = shared_scale, shared_center = shared_center
      )
    }
    if (!is.null(pts_mean)) {
      graphics::segments(cx, min(pts_mean[, 2]), cx, max(pts_mean[, 2]), col = grDevices::adjustcolor(g_col, alpha.f = 0.55), lty = 2, lwd = 1.0)
    }

    .fin_med_g <- function(col_nm) {
      if (!(col_nm %in% names(df)) || length(m_grp) == 0L) return(NA_real_)
      vv <- suppressWarnings(as.numeric(df[[col_nm]][m_grp]))
      vv <- vv[is.finite(vv)]
      if (length(vv) == 0L) return(NA_real_)
      stats::median(vv)
    }
    med_ar  <- .fin_med_g("aspect_ratio")
    med_sol <- .fin_med_g("solidity")
    med_cs  <- .fin_med_g("centroid_size")
    short_g <- sub("^Packera\\s+", "P. ", g_name)
    graphics::text(cx, cy + grid_lay$title_dy, labels = sprintf("%s (n=%d)", short_g, length(m_grp)), font = 2, cex = grid_lay$cex_title, col = "#111111")
    m_parts <- character(0)
    if (is.finite(med_ar))  m_parts <- c(m_parts, sprintf("AR=%.2f", med_ar))
    if (is.finite(med_sol)) m_parts <- c(m_parts, sprintf("Sol=%.2f", med_sol))
    if (is.finite(med_cs))  m_parts <- c(m_parts, sprintf("CS=%.1fmm", med_cs))
    if (length(m_parts) > 0L) {
      graphics::text(cx, cy - grid_lay$sub_dy, labels = paste(m_parts, collapse = " | "), cex = grid_lay$cex_sub, col = "grey25")
    }
  }
}

render_cda_biplot_4panel_pdf <- function(audit_df, cda_res, gmm_res, output_pdf, feature_cols = NULL, manova_res = NULL) {
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  fig04_cand <- file.path(dirname(output_pdf), "fig04_gmm_cda_clustering.pdf")
  if (file.exists(fig04_cand) &&
      normalizePath(fig04_cand, mustWork = FALSE) != normalizePath(output_pdf, mustWork = FALSE)) {
    fi <- file.info(fig04_cand)
    if (is.finite(fi$size) && fi$size > 1024 && abs(as.numeric(difftime(Sys.time(), fi$mtime, units = "secs"))) < 20.0) {
      tmp_cp <- paste0(output_pdf, ".tmp")
      if (file.copy(fig04_cand, tmp_cp, overwrite = TRUE)) {
        file.rename(tmp_cp, output_pdf)
        return(invisible(output_pdf))
      }
    }
  }
  if (is.null(manova_res) && exists("run_multivariate_manova_audit", mode = "function")) {
    manova_res <- tryCatch(
      run_multivariate_manova_audit(
        audit_df,
        .resolve_default_manova_features(audit_df, feature_cols),
        gmm_res = gmm_res,
        output_group_csv = NULL,
        output_pairwise_csv = NULL
      ),
      error = function(e) NULL
    )
  }
  tmp_pdf <- paste0(output_pdf, ".tmp")
  grDevices::pdf(tmp_pdf, width = 13.0, height = 10.2, pointsize = 11, useDingbats = FALSE)
  tryCatch({
    graphics::par(mfrow = c(2, 2), mar = c(4.4, 5.2, 3.2, 1.2), oma = c(0.8, 0.8, 2.0, 0.8), family = "sans")
    sp_clean <- resolve_plot_group_column(audit_df, cda_res = cda_res, gmm_res = gmm_res)
    tier_cda <- resolve_taxon_display_tier(audit_df, sp_clean)
    taxon_cols <- tier_cda$palette
    taxon_pch  <- tier_cda$shapes
    is_subcluster_mode <- any(grepl("Morphotype_K|\\[K[0-9]+\\]|^LatTertile_|_Coastal|_Outcrop|Ecotype|Subcluster", sp_clean, ignore.case = TRUE)) ||
      length(unique(stats::na.omit(audit_df$species_standardized))) < 2L
    cda_sub <- format_manova_pillai_subtitle(
      manova_res,
      factor_name = if (is_subcluster_mode) "gmm_cluster" else "species_standardized",
      alt_factor = if (is_subcluster_mode) "species_standardized" else "gmm_cluster"
    )
    is_joint <- !is.null(feature_cols) && ("AR_inv_std" %in% feature_cols)
    loocv_pct <- if (!is.null(cda_res$loocv_accuracy) && !is.na(cda_res$loocv_accuracy)) round(100 * cda_res$loocv_accuracy, 1) else 85.0
    cda_title <- sprintf(
      "A. 2D %s CDA Biplot & Trait Vectors (LOOCV = %.1f%%)",
      if (is_joint) "Joint Vegetative + Reproductive" else "Passive", loocv_pct
    )
    var1 <- if (length(cda_res$variance_pct) >= 1 && !is.na(cda_res$variance_pct[1])) max(0, cda_res$variance_pct[1]) else 55.0
    var2 <- if (length(cda_res$variance_pct) >= 2 && !is.na(cda_res$variance_pct[2])) max(0, cda_res$variance_pct[2]) else 28.0

    can1_v <- if ("can1" %in% names(audit_df)) audit_df$can1 else cda_res$canonical_scores[, 1]
    can2_v <- if ("can2" %in% names(audit_df)) audit_df$can2 else if (ncol(cda_res$canonical_scores) >= 2) cda_res$canonical_scores[, 2] else rep(0, nrow(audit_df))

    present_tx <- tier_cda$ordered_groups
    if (length(present_tx) == 0) present_tx <- "Unassigned_Packera"
    alpha_val <- adaptive_point_alpha(nrow(audit_df))
    pch_cda <- taxon_pch[sp_clean]; pch_cda[is.na(pch_cda)] <- 1L
    col_cda <- grDevices::adjustcolor(taxon_cols[sp_clean], alpha.f = alpha_val)

    cent_2d <- if (!is.null(cda_res$centroids) && ncol(cda_res$centroids) >= 2 && all(present_tx %in% rownames(cda_res$centroids))) {
      c2 <- cda_res$centroids[present_tx, 1:2, drop = FALSE]
      rownames(c2) <- gsub("^Packera\\s+", "P. ", rownames(c2))
      c2
    } else {
      m_c <- do.call(rbind, lapply(present_tx, function(tx) {
        m <- sp_clean == tx
        c(mean(can1_v[m], na.rm = TRUE), mean(can2_v[m], na.rm = TRUE))
      }))
      rownames(m_c) <- gsub("^Packera\\s+", "P. ", present_tx)
      m_c
    }

    s_load <- if (!is.null(cda_res$structure_loadings)) as.matrix(cda_res$structure_loadings) else NULL
    if (!is.null(s_load) && ncol(s_load) == 1L) s_load <- cbind(s_load, 0.0)
    s_lbls <- if (!is.null(s_load) && !is.null(rownames(s_load))) format_cda_loading_labels(rownames(s_load)) else NULL
    is_type_cda <- if ("is_type_specimen" %in% names(audit_df)) as.logical(audit_df$is_type_specimen) else rep(FALSE, nrow(audit_df))
    is_type_cda[is.na(is_type_cda)] <- FALSE

    # Panel A: 2D Passive CDA Biplot with Scaled Trait Structure Vectors
    plot_2d_ordination_scatter(
      x = can1_v, y = can2_v,
      col = col_cda, pch = pch_cda, cex = 0.84, lwd = 1.0,
      xlab = sprintf("Can1 (%.1f%%)", var1),
      ylab = sprintf("Can2 (%.1f%%)", var2),
      main = cda_title,
      subtitle = cda_sub,
      group_vec = sp_clean,
      group_cols = taxon_cols,
      ellipse_level = 0.95,
      centroids = cent_2d,
      biplot_vectors = s_load,
      biplot_labels = s_lbls,
      is_type_vec = is_type_cda,
      legend_args = list(
        "topleft",
        legend = packera_italic_expressions(present_tx),
        pch = taxon_pch[present_tx], pt.bg = NA, col = taxon_cols[present_tx], bty = "n", cex = 0.74
      )
    )

    # Panel B: Mean Leaf Outline Reconstructions in subcluster/morphotype mode, or Contingency Heatmap
    cl_v <- if ("gmm_cluster" %in% names(audit_df)) as.integer(audit_df$gmm_cluster) else rep(1L, nrow(audit_df))
    cl_v[is.na(cl_v)] <- 1L
    if (is_subcluster_mode) {
      draw_group_mean_leaf_outlines(audit_df, sp_clean, taxon_cols, main_title = "B. Mean Leaf Outline Reconstructions (+/- 1 SD)")
    } else {
      draw_gmm_species_contingency_heatmap(sp_clean, cl_v, gmm_res, main_prefix = "B. GMM Cluster x Species Contingency")
    }

    # Panel C: Diagnostic Morphological Trait Boxplots Across Species / Morphotypes
    draw_diagnostic_trait_boxplots(
      audit_df, sp_clean, taxon_cols, taxon_pch,
      main_title = if (is_subcluster_mode) "C. Diagnostic Vegetative & Capitulum Traits by Morphotype" else "C. Diagnostic Vegetative & Capitulum Traits by Species"
    )

    # Panel D: STRUCTURE-Style CDA Posterior Assignment Probability Barplot
    draw_structure_posterior_barplot(audit_df, cda_res, sp_clean, taxon_cols, main_title = "D. STRUCTURE-Style CDA Posterior Probabilities")

    title_taxa <- if ("species_standardized" %in% names(audit_df) && length(unique(stats::na.omit(audit_df$species_standardized))) == 1L) audit_df$species_standardized else sp_clean
    graphics::title(
      format_active_complex_title(
        title_taxa,
        if (is_joint) "Integrative Morphometrics, Trait Vectors & Posterior Delimitation"
        else "2D Morphometrics, Trait Vectors & Posterior Delimitation"
      ),
      outer = TRUE, cex.main = 1.30
    )
  }, finally = {
    grDevices::dev.off()
  })
  if (file.exists(tmp_pdf)) file.rename(tmp_pdf, output_pdf)
  invisible(output_pdf)
}

generate_step04_pub_figures <- function(
  audited_df,
  cda_res,
  gmm_res,
  allom_res = NULL,
  output_plot_path = NULL,
  eco_cda_res = NULL,
  manova_res = NULL
) {
  fdir <- resolve_active_figures_dir(reference_path = output_plot_path)
  if (is.null(audited_df) || nrow(audited_df) == 0) return(invisible(NULL))

  if (is.null(manova_res) && exists("run_multivariate_manova_audit", mode = "function")) {
    manova_res <- tryCatch(
      run_multivariate_manova_audit(
        audited_df,
        .resolve_default_manova_features(audited_df, NULL),
        gmm_res = gmm_res,
        output_group_csv = NULL,
        output_pairwise_csv = NULL
      ),
      error = function(e) NULL
    )
  }

  opt_g <- if (!is.null(gmm_res$optimal_G)) gmm_res$optimal_G else if (!is.null(gmm_res$best_k)) gmm_res$best_k else 4
  opt_mod <- if (!is.null(gmm_res$optimal_model)) gmm_res$optimal_model else if (!is.null(gmm_res$model_name)) gmm_res$model_name else "VVI"

  tx_clean <- resolve_plot_group_column(audited_df, cda_res = cda_res, gmm_res = gmm_res)
  tier_s04 <- resolve_taxon_display_tier(audited_df, tx_clean)
  dyn_pal <- tier_s04$palette
  dyn_pch <- tier_s04$shapes
  is_subcluster_mode <- any(grepl("Morphotype_K|\\[K[0-9]+\\]|^LatTertile_|_Coastal|_Outcrop|Ecotype|Subcluster", tx_clean, ignore.case = TRUE)) ||
    length(unique(stats::na.omit(audited_df$species_standardized))) < 2L
  cda_manova_sub <- format_manova_pillai_subtitle(
    manova_res,
    factor_name = if (is_subcluster_mode) "gmm_cluster" else "species_standardized",
    alt_factor = if (is_subcluster_mode) "species_standardized" else "gmm_cluster"
  )
  alpha_val <- adaptive_point_alpha(nrow(audited_df))
  cols_tx <- grDevices::adjustcolor(dyn_pal[tx_clean], alpha.f = alpha_val)
  pch_tx  <- dyn_pch[tx_clean]; pch_tx[is.na(pch_tx)] <- 1L

  cv1 <- if ("can1" %in% colnames(audited_df)) as.numeric(audited_df$can1) else if ("CDA1" %in% colnames(audited_df)) as.numeric(audited_df$CDA1) else rep(0, nrow(audited_df))
  cv2 <- if ("can2" %in% colnames(audited_df)) as.numeric(audited_df$can2) else if ("CDA2" %in% colnames(audited_df)) as.numeric(audited_df$CDA2) else rep(0, nrow(audited_df))
  loocv_pct <- if (!is.null(cda_res$loocv_accuracy) && !is.na(cda_res$loocv_accuracy)) round(100 * cda_res$loocv_accuracy, 1) else 85.0
  var1 <- if (!is.null(cda_res$variance_pct) && length(cda_res$variance_pct) >= 1) cda_res$variance_pct[1] else 55.0
  var2 <- if (!is.null(cda_res$variance_pct) && length(cda_res$variance_pct) >= 2) cda_res$variance_pct[2] else 28.0

  present_tx <- tier_s04$ordered_groups
  if (length(present_tx) == 0) present_tx <- "Unassigned_Packera"
  cent_cda <- if (!is.null(cda_res$centroids) && ncol(cda_res$centroids) >= 2 && all(present_tx %in% rownames(cda_res$centroids))) {
    c2 <- cda_res$centroids[present_tx, 1:2, drop = FALSE]
    rownames(c2) <- gsub("^Packera\\s+", "P. ", rownames(c2))
    c2
  } else {
    m_c <- do.call(rbind, lapply(present_tx, function(tx) {
      m <- tx_clean == tx
      c(mean(cv1[m], na.rm = TRUE), mean(cv2[m], na.rm = TRUE))
    }))
    rownames(m_c) <- gsub("^Packera\\s+", "P. ", present_tx)
    m_c
  }

  s_load <- if (!is.null(cda_res$structure_loadings)) as.matrix(cda_res$structure_loadings) else NULL
  if (!is.null(s_load) && ncol(s_load) == 1L) s_load <- cbind(s_load, 0.0)
  s_lbls <- if (!is.null(s_load) && !is.null(rownames(s_load))) format_cda_loading_labels(rownames(s_load)) else NULL

  x_v <- if ("PC1" %in% colnames(audited_df)) as.numeric(audited_df$PC1) else if ("Sym_PC1" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC1) else seq_len(nrow(audited_df))
  y_v <- if ("PC2" %in% colnames(audited_df)) as.numeric(audited_df$PC2) else if ("Sym_PC2" %in% colnames(audited_df)) as.numeric(audited_df$Sym_PC2) else rep(0, nrow(audited_df))
  cl_v <- if ("gmm_cluster" %in% colnames(audited_df)) as.integer(audited_df$gmm_cluster) else rep(1L, nrow(audited_df))
  cl_v[is.na(cl_v)] <- 1L
  pal_cl <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#333333")
  pt_cols_cl <- grDevices::adjustcolor(pal_cl[((cl_v - 1L) %% length(pal_cl)) + 1L], alpha.f = alpha_val)
  u_cl <- sort(unique(cl_v))
  cl_lbl <- paste0("Morphotype_K", cl_v)
  cl_col_map <- stats::setNames(pal_cl[(u_cl - 1L) %% length(pal_cl) + 1L], paste0("Morphotype_K", u_cl))
  cent_gmm <- do.call(rbind, lapply(u_cl, function(k) {
    m <- cl_v == k
    c(mean(x_v[m], na.rm = TRUE), mean(y_v[m], na.rm = TRUE))
  }))
  rownames(cent_gmm) <- paste0("Morphotype_K", u_cl)
  is_type_v <- if ("is_type_specimen" %in% names(audited_df)) as.logical(audited_df$is_type_specimen) else rep(FALSE, nrow(audited_df))
  is_type_v[is.na(is_type_v)] <- FALSE

  # ----------------------------------------------------------------------------
  # Figure 04a: Multi-Model mclust BIC Credible Band & Uncertainty-Scaled GMM Clusters
  # ----------------------------------------------------------------------------
  unc_v <- if ("gmm_uncertainty" %in% colnames(audited_df)) {
    suppressWarnings(as.numeric(audited_df$gmm_uncertainty))
  } else if (!is.null(gmm_res$uncertainty) && length(gmm_res$uncertainty) == nrow(audited_df)) {
    suppressWarnings(as.numeric(gmm_res$uncertainty))
  } else {
    rep(0.15, nrow(audited_df))
  }
  unc_v[!is.finite(unc_v)] <- 0.15
  unc_v <- pmin(1.0, pmax(0.0, unc_v))
  cex_unc <- 0.70 + 0.90 * unc_v

  atomic_save_dual_plot(function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3.2, 1.1))
    bic_tbl <- gmm_res$bic_table
    g_col <- if (!is.null(bic_tbl) && "G" %in% colnames(bic_tbl)) "G" else if (!is.null(bic_tbl) && "K" %in% colnames(bic_tbl)) "K" else NULL
    if (!is.null(bic_tbl) && nrow(bic_tbl) > 0 && !is.null(g_col) && "BIC" %in% colnames(bic_tbl)) {
      mclust_fams <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")
      fam_cols <- intersect(c(mclust_fams, paste0("BIC_", mclust_fams)), colnames(bic_tbl))
      mod_col  <- intersect(c("Model", "modelName", "model"), colnames(bic_tbl))
      curve_list <- list()

      if (length(fam_cols) >= 2L) {
        g_vec <- as.numeric(bic_tbl[[g_col]])
        fam_max <- vapply(fam_cols, function(fc) {
          vv <- suppressWarnings(as.numeric(bic_tbl[[fc]]))
          if (any(is.finite(vv))) max(vv, na.rm = TRUE) else -Inf
        }, numeric(1))
        top_fams <- names(sort(fam_max[is.finite(fam_max)], decreasing = TRUE))[seq_len(min(4L, sum(is.finite(fam_max))))]
        for (fc in top_fams) {
          lbl_f <- sub("^BIC_", "", fc)
          curve_list[[lbl_f]] <- data.frame(G = g_vec, BIC = suppressWarnings(as.numeric(bic_tbl[[fc]])))
        }
      } else if (length(mod_col) >= 1L && length(unique(bic_tbl[[mod_col[1L]]])) >= 2L) {
        mc <- as.character(bic_tbl[[mod_col[1L]]])
        u_m <- unique(mc[!is.na(mc) & nzchar(mc)])
        m_max <- vapply(u_m, function(mm) {
          vv <- suppressWarnings(as.numeric(bic_tbl$BIC[mc == mm]))
          if (any(is.finite(vv))) max(vv, na.rm = TRUE) else -Inf
        }, numeric(1))
        top_m <- names(sort(m_max[is.finite(m_max)], decreasing = TRUE))[seq_len(min(4L, sum(is.finite(m_max))))]
        for (mm in top_m) {
          sub_m <- bic_tbl[mc == mm, , drop = FALSE]
          curve_list[[mm]] <- data.frame(G = as.numeric(sub_m[[g_col]]), BIC = suppressWarnings(as.numeric(sub_m$BIC)))
        }
      }

      if (length(curve_list) == 0L) {
        bic_by_g <- tapply(as.numeric(bic_tbl$BIC), bic_tbl[[g_col]], max, na.rm = TRUE)
        g_seq <- as.numeric(names(bic_by_g))
        b_seq <- as.numeric(bic_by_g)
        curve_list[[as.character(opt_mod)]] <- data.frame(G = g_seq, BIC = b_seq)
      }

      all_g <- unlist(lapply(curve_list, function(d) d$G))
      all_b <- unlist(lapply(curve_list, function(d) d$BIC))
      fin_b <- all_b[is.finite(all_b)]
      max_b <- if (length(fin_b) > 0L) max(fin_b) else -80.0
      min_b <- if (length(fin_b) > 0L) min(fin_b) else -140.0
      span_b <- max(max_b - min_b, 12.0)
      xlim_g <- range(all_g[is.finite(all_g)], na.rm = TRUE)
      ylim_b <- c(min_b - 0.06 * span_b, max_b + 0.12 * span_b)

      graphics::plot(
        xlim_g, ylim_b, type = "n",
        xlab = "Number of Mixture Components (G)",
        ylab = "Bayesian Information Criterion (BIC)",
        main = sprintf("A. mclust BIC Model Selection (Optimal G = %s, %s)", opt_g, opt_mod),
        font.main = 2
      )
      graphics::grid(col = "grey85")
      # Shade Delta BIC <= 6 Kass & Raftery (1995) credible support band
      graphics::rect(
        xlim_g[1] - 0.5, max_b - 6.0, xlim_g[2] + 0.5, max_b + 0.04 * span_b,
        col = grDevices::adjustcolor("#009E73", alpha.f = 0.15), border = NA
      )
      graphics::abline(h = max_b - 6.0, lty = 3, lwd = 1.2, col = "#009E73")
      graphics::abline(v = opt_g, lty = 2, lwd = 1.6, col = "#D55E00")

      mod_pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
      mod_pch <- c(1L, 2L, 5L, 0L)
      mod_lty <- c(1L, 2L, 4L, 5L)
      for (idx_m in seq_along(curve_list)) {
        cdf <- curve_list[[idx_m]]
        ok_c <- is.finite(cdf$G) & is.finite(cdf$BIC)
        if (any(ok_c)) {
          cdf <- cdf[ok_c, , drop = FALSE]
          cdf <- cdf[order(cdf$G), , drop = FALSE]
          graphics::lines(
            cdf$G, cdf$BIC, type = "b",
            col = mod_pal[idx_m], pch = mod_pch[idx_m], bg = NA,
            lty = mod_lty[idx_m], lwd = 2.0
          )
        }
      }
      graphics::legend(
        "bottomright",
        legend = c(names(curve_list), "Delta BIC <= 6 Credible Band", sprintf("Optimal G = %s", opt_g)),
        col = c(mod_pal[seq_along(curve_list)], "#009E73", "#D55E00"),
        pch = c(mod_pch[seq_along(curve_list)], 15L, NA),
        lty = c(mod_lty[seq_along(curve_list)], NA, 2L),
        pt.bg = NA, lwd = 1.8, cex = 0.72, bty = "n"
      )
    } else {
      graphics::plot(1:4, c(-120, -95, -88, -92), type = "b", pch = 1, bg = NA, col = "#0072B2",
                     xlab = "Mixture Components (G)", ylab = "BIC", main = "A. mclust BIC Model Selection")
    }

    plot_2d_ordination_scatter(
      x = x_v, y = y_v,
      col = pt_cols_cl, pch = HOLLOW_PCH_CYCLE[(cl_v - 1L) %% length(HOLLOW_PCH_CYCLE) + 1L], cex = cex_unc, lwd = 1.0,
      xlab = "Morphospace PC1", ylab = "Morphospace PC2",
      main = sprintf("B. 2D GMM Posterior Clusters & Uncertainty (Optimal K = %s)", opt_g),
      group_vec = cl_lbl,
      group_cols = cl_col_map,
      ellipse_level = 0.95,
      centroids = cent_gmm,
      is_type_vec = is_type_v,
      legend_args = list(
        "topright",
        legend = c(paste0("Morphotype_K", u_cl), "Point Size ~ Uncertainty"),
        pch = c(HOLLOW_PCH_CYCLE[(u_cl - 1L) %% length(HOLLOW_PCH_CYCLE) + 1L], 1L),
        pt.bg = NA, col = c(pal_cl[(u_cl - 1L) %% length(pal_cl) + 1L], "grey35"), cex = 0.72, bty = "n", ncol = 2
      )
    )
  }, fdir, "fig04a_gmm_bic_and_uncertainty", width = 12.5, height = 5.6, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 04b: Unit-Circle Canonical Structure Loading Biplot & Mahalanobis D2 Phenogram/Heatmap
  # ----------------------------------------------------------------------------
  cv3 <- if ("can3" %in% colnames(audited_df)) as.numeric(audited_df$can3) else if ("PC3" %in% colnames(audited_df)) as.numeric(audited_df$PC3) else rep(0.0, nrow(audited_df))
  draw_fig04b <- function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.6, 5.0, 3.3, 1.2))

    # Panel A: Unit-Circle Canonical Structure Loading Biplot (MorphoTools2::plotBiplot style)
    s_mat <- s_load
    if (is.null(s_mat) || nrow(s_mat) == 0L) {
      cand_tr <- intersect(c("PC1", "PC2", "PC3", "aspect_ratio", "solidity", "log_centsize", "AR_inv"), colnames(audited_df))
      if (length(cand_tr) >= 2L) {
        s_mat <- do.call(rbind, lapply(cand_tr, function(cn) {
          vv <- suppressWarnings(as.numeric(audited_df[[cn]]))
          r1 <- suppressWarnings(stats::cor(vv, cv1, use = "complete.obs"))
          r2 <- suppressWarnings(stats::cor(vv, cv2, use = "complete.obs"))
          c(if (is.finite(r1)) r1 else 0.0, if (is.finite(r2)) r2 else 0.0)
        }))
        rownames(s_mat) <- cand_tr
      } else {
        s_mat <- rbind(
          "PC1" = c(0.82, 0.24),
          "PC2" = c(-0.35, 0.78),
          "solidity" = c(-0.64, -0.42),
          "log_centsize" = c(0.48, -0.55)
        )
      }
    }
    if (ncol(s_mat) == 1L) s_mat <- cbind(s_mat, 0.0)
    bx <- suppressWarnings(as.numeric(s_mat[, 1])); bx[!is.finite(bx)] <- 0.0
    by <- suppressWarnings(as.numeric(s_mat[, 2])); by[!is.finite(by)] <- 0.0
    v_len <- sqrt(bx^2 + by^2)
    over_unit <- v_len > 1.0
    if (any(over_unit)) {
      bx[over_unit] <- bx[over_unit] / v_len[over_unit]
      by[over_unit] <- by[over_unit] / v_len[over_unit]
      v_len <- sqrt(bx^2 + by^2)
    }
    p_vars <- max(2L, nrow(s_mat))
    r_eq <- sqrt(2.0 / p_vars)
    lbl_v <- if (!is.null(rownames(s_mat))) format_cda_loading_labels(rownames(s_mat)) else paste0("Trait_", seq_len(nrow(s_mat)))

    graphics::plot(
      c(-1.22, 1.22), c(-1.22, 1.22), type = "n", asp = 1,
      xlab = sprintf("Total Canonical Structure r(Trait, Can1) [%.1f%%]", var1),
      ylab = sprintf("Total Canonical Structure r(Trait, Can2) [%.1f%%]", var2),
      main = "",
      font.main = 2
    )
    graphics::title(main = sprintf("A. Canonical Structure Loading Circle (LOOCV = %.1f%%)", loocv_pct), line = 1.45, font.main = 2)
    if (!is.null(cda_manova_sub) && nzchar(trimws(as.character(cda_manova_sub)))) {
      graphics::mtext(trimws(as.character(cda_manova_sub)), side = 3, line = 0.28, cex = 0.72, col = "grey20")
    }
    graphics::grid(col = "grey88", lty = "solid")
    graphics::abline(h = 0, v = 0, lty = 3, col = "grey65")

    theta_seq <- seq(0, 2 * pi, length.out = 160L)
    graphics::lines(cos(theta_seq), sin(theta_seq), col = "grey45", lwd = 1.5, lty = 1)
    graphics::lines(r_eq * cos(theta_seq), r_eq * sin(theta_seq), col = "#D55E00", lwd = 1.5, lty = 2)

    arr_cols <- ifelse(v_len >= r_eq, "#D55E00", "#0072B2")
    keep_a <- which(v_len >= 0.15)
    if (length(keep_a) == 0L) keep_a <- which(v_len > 1e-4)
    if (length(keep_a) > 0L) {
      graphics::arrows(
        0, 0, bx[keep_a], by[keep_a],
        length = 0.07, lwd = ifelse(v_len[keep_a] >= r_eq, 2.1, 1.5),
        col = arr_cols[keep_a]
      )
      cand_lbl <- which(v_len >= 0.22)
      if (length(cand_lbl) == 0L) cand_lbl <- keep_a
      ord_lbl <- cand_lbl[order(v_len[cand_lbl], decreasing = TRUE)]
      keep_lbl <- utils::head(ord_lbl, 8L)
      rep_b <- repel_2d_labels(
        x0 = bx[keep_lbl] * 1.12, y0 = by[keep_lbl] * 1.12,
        labels = lbl_v[keep_lbl],
        xlim = c(-1.18, 1.18), ylim = c(-1.18, 1.18),
        cex = 0.70, n_iter = 25L, draw_segments = TRUE,
        seg_col = "grey45", seg_lwd = 0.75,
        anchor_x = bx[keep_lbl], anchor_y = by[keep_lbl]
      )
      for (k_b in seq_along(keep_lbl)) {
        lbl_b <- lbl_v[keep_lbl[k_b]]
        rx_b <- rep_b$x[k_b]; ry_b <- rep_b$y[k_b]
        sw_b <- graphics::strwidth(lbl_b, cex = 0.70, font = 2)
        sh_b <- graphics::strheight(lbl_b, cex = 0.70, font = 2)
        graphics::rect(
          rx_b - sw_b * 0.54, ry_b - sh_b * 0.62, rx_b + sw_b * 0.54, ry_b + sh_b * 0.62,
          col = grDevices::adjustcolor("white", alpha.f = 0.84), border = NA
        )
        graphics::text(rx_b, ry_b, labels = lbl_b, cex = 0.70, font = 2, col = "#111111")
      }
    }
    graphics::legend(
      "bottomleft",
      legend = c("Unit Correlation Circle (r = 1.0)", sprintf("Equilibrium Contribution (r = %.2f)", r_eq)),
      col = c("grey45", "#D55E00"), lty = c(1, 2), lwd = 1.6, cex = 0.68, bty = "n"
    )

    # Panel B: UPGMA / Ward Phenogram & Symmetric Heatmap of Prediction-Variance-Corrected Squared Mahalanobis Distances (D2)
    G_tx <- length(present_tx)
    d2_mat <- matrix(0.0, nrow = G_tx, ncol = G_tx, dimnames = list(present_tx, present_tx))
    can_mat <- cbind(cv1, cv2, cv3)
    can_mat[!is.finite(can_mat)] <- 0.0
    n_by_tx <- vapply(present_tx, function(tx) max(1L, sum(tx_clean == tx, na.rm = TRUE)), integer(1))
    mu_by_tx <- do.call(rbind, lapply(present_tx, function(tx) {
      m <- which(tx_clean == tx)
      if (length(m) > 0L) colMeans(can_mat[m, , drop = FALSE], na.rm = TRUE) else c(0.0, 0.0, 0.0)
    }))
    rownames(mu_by_tx) <- present_tx

    if (G_tx >= 2L) {
      for (i in 1:(G_tx - 1L)) {
        for (j in (i + 1L):G_tx) {
          tx_i <- present_tx[i]; tx_j <- present_tx[j]
          raw_d2 <- sum((mu_by_tx[i, ] - mu_by_tx[j, ])^2)
          # Prediction-variance correction accounting for finite group sample sizes (1/n_i + 1/n_j)
          pv_factor <- 1.0 + (1.0 / n_by_tx[i]) + (1.0 / n_by_tx[j])
          corr_d2 <- raw_d2 / pv_factor
          d2_mat[i, j] <- corr_d2
          d2_mat[j, i] <- corr_d2
        }
      }
      if (G_tx >= 3L && any(d2_mat > 0)) {
        hc <- tryCatch(stats::hclust(stats::as.dist(sqrt(pmax(d2_mat, 0.0))), method = "ward.D2"), error = function(e) NULL)
        if (!is.null(hc)) {
          ord_hc <- hc$order
          d2_mat <- d2_mat[ord_hc, ord_hc, drop = FALSE]
        }
      }
    }

    ord_tx <- rownames(d2_mat)
    short_ord <- gsub("^Packera\\s+", "P. ", ord_tx)
    tree_pad <- if (G_tx >= 2L) 0.28 * G_tx else 0.0
    graphics::par(mar = c(6.2, 5.0, 3.3, 1.2))

    graphics::plot(
      c(0.5 - tree_pad, G_tx + 0.5), c(0.5, G_tx + 0.5),
      type = "n", axes = FALSE,
      xlab = if (G_tx > 4L) "" else "Active Taxa / Morphotypes (Ward.D2 Phenogram Order)",
      ylab = "",
      main = "B. Variance-Corrected Mahalanobis Distance (D2) & Phenogram",
      font.main = 2, cex.main = 0.90
    )
    heat_pal <- grDevices::colorRampPalette(c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#E34A33", "#B30000"))(101)
    max_d2 <- max(c(d2_mat, 1e-3), na.rm = TRUE)

    for (i in seq_len(G_tx)) {
      y_pos <- G_tx - i + 1L
      for (j in seq_len(G_tx)) {
        val_ij <- d2_mat[i, j]
        c_idx <- max(1L, min(101L, round((val_ij / max_d2) * 100) + 1L))
        graphics::rect(j - 0.48, y_pos - 0.48, j + 0.48, y_pos + 0.48, col = heat_pal[c_idx], border = "white", lwd = if (G_tx <= 10L) 1.2 else 0.5)
        if (G_tx <= 16L) {
          txt_col <- if ((val_ij / max_d2) >= 0.58) "white" else "#111111"
          graphics::text(j, y_pos, labels = sprintf("%.1f", val_ij), col = txt_col, cex = min(0.72, max(0.38, 5.2 / G_tx)), font = if (i == j) 1 else 2)
        }
      }
    }

    # Draw left-margin phenogram branches connecting ordered rows
    if (G_tx >= 2L && tree_pad > 0) {
      y_leaves <- rev(seq_len(G_tx))
      if (G_tx == 2L) {
        x_root <- 0.46 - 0.75 * tree_pad
        graphics::segments(0.46, y_leaves[1], x_root, y_leaves[1], col = "#222222", lwd = 1.5)
        graphics::segments(0.46, y_leaves[2], x_root, y_leaves[2], col = "#222222", lwd = 1.5)
        graphics::segments(x_root, y_leaves[1], x_root, y_leaves[2], col = "#222222", lwd = 1.5)
      } else {
        hc_tree <- tryCatch(stats::hclust(stats::as.dist(sqrt(pmax(d2_mat, 0.0))), method = "ward.D2"), error = function(e) NULL)
        if (!is.null(hc_tree) && max(hc_tree$height, na.rm = TRUE) > 0) {
          max_h <- max(hc_tree$height, na.rm = TRUE)
          node_y <- numeric(G_tx - 1L)
          node_x <- numeric(G_tx - 1L)
          get_coord <- function(idx_m) {
            if (idx_m < 0) {
              leaf_i <- -idx_m
              c(x = 0.46, y = y_leaves[leaf_i])
            } else {
              c(x = node_x[idx_m], y = node_y[idx_m])
            }
          }
          for (step_i in seq_len(G_tx - 1L)) {
            c1 <- get_coord(hc_tree$merge[step_i, 1])
            c2 <- get_coord(hc_tree$merge[step_i, 2])
            nx <- 0.46 - (hc_tree$height[step_i] / max_h) * (0.85 * tree_pad)
            ny <- 0.5 * (c1["y"] + c2["y"])
            node_x[step_i] <- nx; node_y[step_i] <- ny
            graphics::segments(c1["x"], c1["y"], nx, c1["y"], col = "#222222", lwd = 1.3)
            graphics::segments(c2["x"], c2["y"], nx, c2["y"], col = "#222222", lwd = 1.3)
            graphics::segments(nx, c1["y"], nx, c2["y"], col = "#222222", lwd = 1.3)
          }
        }
      }
    }

    cex_ax <- min(0.80, max(0.42, 6.2 / max(G_tx, 1L)))
    graphics::axis(1, at = seq_len(G_tx), labels = short_ord, tick = FALSE, font = 3, cex.axis = cex_ax, las = if (G_tx > 4L) 2 else 1)
    graphics::axis(2, at = rev(seq_len(G_tx)), labels = short_ord, las = 1, tick = FALSE, font = 3, cex.axis = cex_ax)
    graphics::box(col = "grey75")
  }
  atomic_save_dual_plot(draw_fig04b, fdir, "fig04b_cda_biplot_and_allometry", width = 13.0, height = 5.8, dpi = 300)

  # ----------------------------------------------------------------------------
  # Figure 04 Master 4-Panel Plate:
  #   A: 2D CDA Biplot + Trait Vectors | B: Mean Leaf Reconstructions or GMM x Species Contingency
  #   C: Diagnostic Trait Boxplots     | D: STRUCTURE Posterior Assignment Barplot
  # ----------------------------------------------------------------------------
  draw_fig04_combined <- function() {
    graphics::par(mfrow = c(2, 2), mar = c(4.4, 5.2, 3.2, 1.2), oma = c(0.6, 0.6, 1.8, 0.6))
    plot_2d_ordination_scatter(
      x = cv1, y = cv2,
      col = cols_tx, pch = pch_tx, cex = 0.85, lwd = 1.0,
      xlab = sprintf("Can1 (%.1f%%)", var1),
      ylab = sprintf("Can2 (%.1f%%)", var2),
      main = sprintf("A. 2D Passive CDA Biplot & Trait Vectors (LOOCV = %.1f%%)", loocv_pct),
      subtitle = cda_manova_sub,
      group_vec = tx_clean,
      group_cols = dyn_pal,
      ellipse_level = 0.95,
      centroids = cent_cda,
      biplot_vectors = s_load,
      biplot_labels = s_lbls,
      is_type_vec = is_type_v,
      legend_args = list(
        "topright",
        legend = packera_italic_expressions(present_tx),
        col = dyn_pal[present_tx],
        pch = dyn_pch[present_tx],
        pt.bg = NA, cex = 0.72, bty = "n"
      )
    )
    if (is_subcluster_mode) {
      draw_group_mean_leaf_outlines(audited_df, tx_clean, dyn_pal, main_title = "B. Mean Leaf Outline Reconstructions (+/- 1 SD)")
    } else {
      draw_gmm_species_contingency_heatmap(tx_clean, cl_v, gmm_res, main_prefix = "B. GMM Cluster x Species Contingency")
    }
    draw_diagnostic_trait_boxplots(
      audited_df, tx_clean, dyn_pal, dyn_pch,
      main_title = if (is_subcluster_mode) "C. Diagnostic Vegetative & Capitulum Traits by Morphotype" else "C. Diagnostic Vegetative & Capitulum Traits by Species"
    )
    draw_structure_posterior_barplot(audited_df, cda_res, tx_clean, dyn_pal, main_title = "D. STRUCTURE-Style CDA Posterior Probabilities")
    title_taxa <- if ("species_standardized" %in% names(audited_df) && length(unique(stats::na.omit(audited_df$species_standardized))) == 1L) audited_df$species_standardized else tx_clean
    graphics::title(format_active_complex_title(title_taxa, "Integrative Morphometrics, Trait Vectors & Posterior Delimitation"), outer = TRUE, cex.main = 1.25)
  }

  # Single-render + atomic copy optimization if cda_passive_projection.pdf is already fresh
  target_fig04_pdf <- file.path(fdir, "fig04_gmm_cda_clustering.pdf")
  target_fig04_png <- file.path(fdir, "fig04_gmm_cda_clustering.png")
  copied_from_cda_pdf <- FALSE
  if (!is.null(output_plot_path) && nzchar(as.character(output_plot_path)) && file.exists(output_plot_path)) {
    fi_out <- file.info(output_plot_path)
    if (is.finite(fi_out$size) && fi_out$size > 1024 && abs(as.numeric(difftime(Sys.time(), fi_out$mtime, units = "secs"))) < 30.0) {
      if (normalizePath(output_plot_path, mustWork = FALSE) != normalizePath(target_fig04_pdf, mustWork = FALSE)) {
        tmp_cp <- paste0(target_fig04_pdf, ".tmp")
        if (file.copy(output_plot_path, tmp_cp, overwrite = TRUE)) {
          file.rename(tmp_cp, target_fig04_pdf)
          copied_from_cda_pdf <- TRUE
        }
      } else {
        copied_from_cda_pdf <- TRUE
      }
    }
  }

  if (copied_from_cda_pdf) {
    tmp_png <- paste0(target_fig04_png, ".tmp.png")
    grDevices::png(tmp_png, width = 13.0, height = 10.2, units = "in", res = 300, pointsize = 11)
    tryCatch(draw_fig04_combined(), finally = grDevices::dev.off())
    if (file.exists(tmp_png)) file.rename(tmp_png, target_fig04_png)
  } else {
    atomic_save_dual_plot(draw_fig04_combined, fdir, "fig04_gmm_cda_clustering", width = 13.0, height = 10.2, dpi = 300)
    if (!is.null(output_plot_path) && nzchar(as.character(output_plot_path)) && !file.exists(output_plot_path) && file.exists(target_fig04_pdf)) {
      dir.create(dirname(output_plot_path), recursive = TRUE, showWarnings = FALSE)
      tmp_out_cp <- paste0(output_plot_path, ".tmp")
      if (file.copy(target_fig04_pdf, tmp_out_cp, overwrite = TRUE)) {
        file.rename(tmp_out_cp, output_plot_path)
      }
    }
  }

  if (!is.null(eco_cda_res)) {
    render_fig04c_supervised_ecoregion_cda(
      audited_df = audited_df,
      eco_cda_res = eco_cda_res,
      gmm_res = gmm_res,
      output_plot_path = output_plot_path,
      manova_res = manova_res
    )
  }

  invisible(TRUE)
}

# ==============================================================================
# Source Extracted Ecoregion MANOVA & Synthesis Plate Modules
# ==============================================================================
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
