#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_figures_utils.R
# Project: Packera Genus-Wide Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   Shared publication figure utilities across R pipeline phases (Steps 03, 04,
#   06, and 07). Enforces:
#     - 300-DPI .png atomic exports (write-then-replace)
#     - Colorblind-safe Okabe-Ito palettes (#0072B2, #D55E00, #009E73, #CC79A7)
#     - Dynamically italicized Plantmatics/Packera taxon labels
# ==============================================================================

if (grDevices::dev.cur() == 1L) {
  grDevices::pdf(NULL)
}

EXTENDED_COLORBLIND_HEX <- c(
  "#E69F00", "#56B4E9", "#F0E442", "#333333", "#882255",
  "#44AA99", "#999933", "#332288", "#117733", "#AA4499",
  "#DDCC77", "#661100", "#6699CC", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#008080"
)

if (!exists("fast_euclidean_permanova_trace", mode = "function")) {
  for (.dl_p in c("scripts/morphometrics/dissected_leaf_utils.R", file.path(dirname(c(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)), unlist(lapply(sys.frames(), function(fr) fr$ofile)))), "dissected_leaf_utils.R"))) {
    if (!is.na(.dl_p) && nzchar(.dl_p) && file.exists(.dl_p)) { source(.dl_p); break }
  }
}

compute_spatial_balance_weights <- function(df, bandwidth_km = 25.0, regional_bandwidth_km = 75.0, active_mask = NULL, recompute = FALSE) {
  if (is.null(df) || nrow(df) == 0L) return(numeric(0))
  n_rows <- nrow(df)
  if (!isTRUE(recompute) && is.null(active_mask) && ("spatial_balance_weight" %in% names(df))) {
    w_in <- suppressWarnings(as.numeric(df$spatial_balance_weight))
    pos_w <- w_in[is.finite(w_in) & w_in > 0]
    if (length(pos_w) > 0L && any(abs(pos_w - 1.0) > 1e-6)) {
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
  if ("is_centroid_or_inst_coord" %in% names(df)) {
    c_raw <- df$is_centroid_or_inst_coord
    is_cent <- if (is.logical(c_raw)) !is.na(c_raw) & c_raw else tolower(trimws(as.character(c_raw))) %in% c("true", "1", "t", "yes")
    ok_geo <- ok_geo & !is_cent
  }
  if ("authority_weight" %in% names(df)) {
    aw_num <- suppressWarnings(as.numeric(df$authority_weight))
    if (any(is.finite(aw_num) & aw_num > 0.0)) ok_geo <- ok_geo & is.finite(aw_num) & (aw_num > 0.0)
  } else if (is.null(active_mask) && ("determiner_tier" %in% names(df))) {
    dt_chr <- trimws(as.character(df$determiner_tier))
    is_tp <- if ("is_type_specimen" %in% names(df)) {
      tr <- df$is_type_specimen
      if (is.logical(tr)) !is.na(tr) & tr else tolower(trimws(as.character(tr))) %in% c("true", "1", "t", "yes")
    } else rep(FALSE, n_rows)
    is_pr <- if ("is_promoted_bronze" %in% names(df)) {
      pr <- df$is_promoted_bronze
      if (is.logical(pr)) !is.na(pr) & pr else tolower(trimws(as.character(pr))) %in% c("true", "1", "t", "yes")
    } else rep(FALSE, n_rows)
    if (any((dt_chr != "Tier_3_Bronze") | is_tp | is_pr, na.rm = TRUE)) {
      ok_geo <- ok_geo & ((dt_chr != "Tier_3_Bronze") | is_tp | is_pr)
    }
  }

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
    "fine_soil_ph"                     = "Rhizosphere Soil pH in H2O (0-15cm)",
    "ilr_texture_1"                    = "Soil Texture ILR-1: Sand vs Silt-Clay (0-15cm)",
    "ilr_texture_2"                    = "Soil Texture ILR-2: Silt vs Clay (0-15cm)",
    "fine_soil_sand"                   = "Soil Sand Fraction (%, 0-15cm)",
    "fine_soil_clay"                   = "Soil Clay Fraction (%, 0-15cm)",
    "fine_soil_cec"                    = "Rhizosphere Soil CEC (cmol+/kg, 0-15cm)",
    "fine_cfvo_pct"                    = "Coarse Fragment Volume (%, 0-15cm)",
    "fine_bedrock_depth_cm"            = "Depth to Lithic Bedrock (cm)",
    "is_azonal_lithology"              = "USGS SGMC Azonal Lithology Indicator",
    "ssurgo_rock_outcrop_pct"          = "SSURGO Rock Outcrop Component (%)",
    "ssurgo_bedrock_depth_cm"          = "SSURGO Lithic Bedrock Depth (cm)",
    "ssurgo_aws_0_50cm"                = "Available Water Storage 0-50cm (cm)",
    "topo_slope_deg"                   = "30m Terrain Slope (deg)",
    "topo_tpi_30m"                     = "30m Topographic Position Index (TPI, m)",
    "topo_hli_30m"                     = "30m Heat Load Index (HLI, McCune & Keon 2002)",
    "topo_twi_30m"                     = "30m D8 Topographic Wetness Index (TWI)",
    "soil_ph"                          = "SoilGrids 250m Soil pH (H2O)",
    "soil_sand"                        = "SoilGrids 250m Sand Fraction (%)",
    "soil_clay"                        = "SoilGrids 250m Clay Fraction (%)",
    "soil_silt"                        = "SoilGrids 250m Silt Fraction (%)",
    "soil_cec"                         = "SoilGrids 250m Cation Exchange (cmol+/kg)",
    "soil_cfvo"                        = "SoilGrids 250m Coarse Fragments (%)",
    "soil_bulk_density"                = "SoilGrids 250m Bulk Density (g/cm3)",
    "bio1_temp_mean"                   = "Mean Annual Temp (BIO1, deg C)",
    "bio2_mean_diurnal_range"          = "Mean Diurnal Temp Range (BIO2, deg C)",
    "bio3_isothermality"               = "Isothermality (BIO3, BIO2/BIO7 x 100)",
    "bio4_temp_seasonality"            = "Temp Seasonality (BIO4, SD x 100)",
    "bio5_max_temp_warmest_month"      = "Max Temp Warmest Month (BIO5, deg C)",
    "bio6_min_temp_coldest_month"      = "Min Temp Coldest Month (BIO6, deg C)",
    "bio7_temp_annual_range"           = "Temp Annual Range (BIO7, deg C)",
    "bio8_mean_temp_wettest_qtr"       = "Mean Temp Wettest Quarter (BIO8, deg C)",
    "bio9_mean_temp_driest_qtr"        = "Mean Temp Driest Quarter (BIO9, deg C)",
    "bio10_mean_temp_warmest_qtr"      = "Mean Temp Warmest Quarter (BIO10, deg C)",
    "bio11_mean_temp_coldest_qtr"      = "Mean Temp Coldest Quarter (BIO11, deg C)",
    "bio12_precip_annual"              = "Annual Precipitation (BIO12, mm)",
    "bio13_precip_wettest_month"       = "Wettest-Month Precip (BIO13, mm)",
    "bio14_precip_driest_month"        = "Driest-Month Precip (BIO14, mm)",
    "bio15_precip_seasonality"         = "Precipitation Seasonality (BIO15, CV)",
    "bio16_precip_wettest_qtr"         = "Wettest-Quarter Precip (BIO16, mm)",
    "bio17_precip_driest_qtr"          = "Driest-Quarter Precip (BIO17, mm)",
    "bio18_precip_warmest_qtr"         = "Warmest-Quarter Precip (BIO18, mm)",
    "bio19_precip_coldest_qtr"         = "Coldest-Quarter Precip (BIO19, mm)",
    "lithic_outcrop_index"             = "Lithic Outcrop Exposure Index (0-100)",
    "summer_moisture_deficit_idx"      = "Summer Moisture Deficit Index",
    "maritime_continentality_idx"      = "Maritime Continentality Index",
    "min_bedrock_depth_buffer_cm"      = "Min Buffered Bedrock Depth (cm)",
    "ph_q10_buffer"                    = "Buffered Soil pH 10th Pct",
    "ph_q90_buffer"                    = "Buffered Soil pH 90th Pct",
    "max_cfvo_buffer_pct"              = "Max Buffered Coarse Fragments (%)",
    "azonal_lithology_presence_buffer" = "Buffered Azonal Lithology Presence",
    "fine_soil_soc"                    = "Soil Organic Carbon (g/kg, 0-15cm)",
    "topo_elev_m"                      = "30m Elevation (m a.s.l.)",
    "delta_doy"                        = "Hopkins Phenological Anomaly (Delta DOY)",
    "doy"                              = "Flowering Day of Year (DOY)",
    "soil_ca_mg_molar_ratio"           = "Soil Ca:Mg Molar Ratio",
    "soil_ultramafic_ni_cr_index"      = "Soil Ultramafic Ni-Cr Index",
    "soil_carbonate_c_pct"             = "Soil Carbonate Carbon (%)",
    "nwi_wetland_dist_m"               = "NWI Wetland Distance (m)",
    "ultramafic_calcareous_gradient"   = "Ultramafic-Calcareous Geochem Gradient",
    "wetland_seepage_index"            = "Wetland & Hydrologic Seepage Index"
  )
  vapply(as.character(var_vec), function(v) {
    if (is.na(v) || !nzchar(v)) return("Unknown Predictor")
    if (v %in% names(dict)) return(unname(dict[[v]]))
    if (grepl("^MEM[0-9]+$", v)) return(sub("^MEM", "Spatial Eigenvector MEM-", v))
    gsub("_", " ", v)
  }, character(1), USE.NAMES = FALSE)
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
  nms <- tolower(gsub("[[:space:]]+", " ", trimws(names(syn_vec))))
  vals <- gsub("^Packera[[:space:]]+(?:x|×)[[:space:]]+", "Packera ", trimws(as.character(syn_vec)), ignore.case = TRUE)
  keep <- !is.na(nms) & nzchar(nms) & !is.na(vals) & nzchar(vals)
  if (!any(keep)) return(character(0))
  out <- stats::setNames(vals[keep], nms[keep])
  bin_pat <- "^([a-z]+[[:space:]]+(?:(?:x|×)[[:space:]]*)?[a-z][a-z-]*)\\b.*$"
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
  env_db <- trimws(Sys.getenv("PLANTMATICS_DB_PATH", unset = ""))
  if (!nzchar(env_db)) env_db <- trimws(Sys.getenv("PACKERA_DB_PATH", unset = ""))
  if (nzchar(env_db) && env_db != ":memory:") cand_dbs <- c(cand_dbs, env_db)
  if (exists("resolve_plantmatics_db_path", mode = "function")) {
    rdb <- tryCatch(as.character(resolve_plantmatics_db_path()), error = function(e) "")
    if (!is.na(rdb) && nzchar(rdb) && rdb != ":memory:") cand_dbs <- c(cand_dbs, rdb)
  } else if (exists("resolve_packera_db_path", mode = "function")) {
    rdb <- tryCatch(as.character(resolve_packera_db_path()), error = function(e) "")
    if (!is.na(rdb) && nzchar(rdb) && rdb != ":memory:") cand_dbs <- c(cand_dbs, rdb)
  }
  cand_dbs <- c(cand_dbs, "data/plantmatics.sqlite")
  if (!file.exists("data/plantmatics.sqlite") && file.exists("data/packera_pipeline.sqlite")) {
    cand_dbs <- c(cand_dbs, "data/packera_pipeline.sqlite")
  }
  if (nzchar(.pub_figures_script_dir)) {
    rel_new <- file.path(.pub_figures_script_dir, "..", "..", "data", "plantmatics.sqlite")
    rel_old <- file.path(.pub_figures_script_dir, "..", "..", "data", "packera_pipeline.sqlite")
    cand_dbs <- c(cand_dbs, rel_new)
    if (!file.exists(rel_new) && file.exists(rel_old)) {
      cand_dbs <- c(cand_dbs, rel_old)
    }
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
  syn_gen_pat <- "^(packera|p\\.|senecio|s\\.|cineraria|jacobaea|aureolaria)[[:space:]]*(?:(?:cf\\.?|aff\\.?)[[:space:]]+)?(x[[:space:]]+|×[[:space:]]*)?([a-z][a-z-]*).*$"

  std_u <- character(length(u_vals))
  for (i in seq_along(u_vals)) {
    val <- u_vals[i]
    if (is.na(val)) { std_u[i] <- NA_character_; next }
    s_orig <- trimws(val)
    s <- tolower(gsub("[[:space:]]+", " ", s_orig))
    if (!nzchar(s_orig) || s %in% c("none", "nan", "null", "na", "unknown", "unassigned", "unassigned_packera", "packera")) {
      std_u[i] <- NA_character_; next
    }
    if (grepl("^Morphotype_K[0-9]+$|^Ecotype_|^LatTertile_|_Coastal$|_Outcrop$|\\[K[0-9]+\\]$", s_orig, ignore.case = TRUE)) {
      std_u[i] <- s_orig; next
    }
    if (has_syn && (s %in% names(syn_map))) {
      std_u[i] <- unname(syn_map[[s]]); next
    }
    if (has_syn && grepl("^[a-z]+[[:space:]]+[a-z-]+[[:space:]]+(?:var\\.?|subsp\\.?|ssp\\.?|f\\.?)[[:space:]]+[a-z-]+", s)) {
      v_raw <- sub("^([a-z]+[[:space:]]+[a-z-]+[[:space:]]+(?:var\\.?|subsp\\.?|ssp\\.?|f\\.?)[[:space:]]+[a-z-]+).*$", "\\1", s)
      v_dot <- gsub("\\b(var|subsp|ssp|f)\\.?[[:space:]]+", "\\1. ", v_raw)
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
standardize_packera_taxon <- function(species_vec, na_value = "Unknown") {
  res <- standardize_taxon_name(species_vec)
  if (!is.null(na_value)) {
    res[is.na(res)] <- na_value
  }
  res
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
      !all(grepl("^Packera[[:space:]]+", vv))
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
  env_run <- trimws(Sys.getenv("PLANTMATICS_RUN_DIR", unset = ""))
  if (!nzchar(env_run)) env_run <- trimws(Sys.getenv("PACKERA_RUN_DIR", unset = ""))
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

# ==============================================================================
# EFA Outline Reconstruction Helper (300-DPI PNG Suites)
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
# Modular 300-DPI PNG Publication Figure Suite Delegation (Steps 03-07)
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
for (.fmod in c("theme_packera_pub.R", "step03_figures.R", "step04_figures.R", "step05_figures.R", "step06_figures.R", "step07_figures.R")) {
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
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") ||
      identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
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

render_allometry_diagnostic_png <- function(sub_df, r2_size, p_size, r2_int, output_png, p_int = 1.0) {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") ||
      identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  out_png <- if (!is.null(output_png) && nzchar(output_png)) normalize_png_path(output_png) else "outputs/figures/fig04a_foliar_allometry_audit.png"
  render_step04_allometry_png(
    sub_df = sub_df,
    r2_size = r2_size,
    p_size = p_size,
    r2_int = r2_int,
    output_path = out_png,
    p_int = p_int
  )
}

render_cda_biplot_4panel_png <- function(audit_df, cda_res, gmm_res, output_png, feature_cols = NULL, manova_res = NULL) {
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") ||
      identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(""))
  out_png <- if (!is.null(output_png) && nzchar(output_png)) normalize_png_path(output_png) else "outputs/figures/fig04b_gmm_cda_morphometrics.png"
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
  if (identical(Sys.getenv("PLANTMATICS_SKIP_FIGURES", unset = ""), "1") ||
      identical(Sys.getenv("PACKERA_SKIP_FIGURES", unset = ""), "1")) return(invisible(TRUE))
  fdir <- resolve_active_figures_dir(figures_dir, reference_path)
  canon_04b <- file.path(fdir, "fig04b_gmm_cda_morphometrics.png")
  ref_png <- if (!is.null(reference_path) && nzchar(reference_path)) normalize_png_path(reference_path) else ""
  if (!file.exists(canon_04b)) {
    if (nzchar(ref_png) && file.exists(ref_png)) {
      dir.create(dirname(canon_04b), recursive = TRUE, showWarnings = FALSE)
      file.copy(ref_png, canon_04b, overwrite = TRUE)
    } else {
      render_step04_gmm_cda_png(
        audit_df = audit_df,
        cda_res = cda_res,
        gmm_res = gmm_res,
        output_path = canon_04b,
        feature_cols = if (!is.null(cda_res)) cda_res$selected_features else NULL,
        manova_res = manova_res
      )
    }
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
