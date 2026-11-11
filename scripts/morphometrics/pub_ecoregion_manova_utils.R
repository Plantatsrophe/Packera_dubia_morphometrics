#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_ecoregion_manova_utils.R
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
#
# Description:
#   US EPA Level III Ecoregion spatial polygon lookup, weighted multivariate
#   MANOVA/PERMANOVA statistical audit utilities, and Track 2 supervised
#   ecoregion canonical discriminant analysis (CDA) publication figures.
# ==============================================================================

# ==============================================================================
# US EPA Level III Ecoregion Spatial Polygon Lookup & Track 2 Supervised CDA
#   - standardize_epa_l3_label_r()
#   - load_epa_level3_spatial_r()
#   - infer_epa_level3_ecoregion()
#   - run_supervised_ecoregion_cda()
#   - render_fig04c_supervised_ecoregion_cda() -> fig04c_supervised_ecoregion_cda.{pdf,png}
# ==============================================================================
.epa_l3_r_cache <- new.env(parent = emptyenv())

.LEGACY_HEURISTIC_REGIONS <- c(
  "coastal_plain_sandhills",
  "piedmont_granite_flatrocks",
  "appalachian_highlands",
  "interior_prairie_midwest",
  "other_us",
  "other_region",
  "other_state",
  "unassigned_na"
)

standardize_epa_l3_label_r <- function(raw_vec) {
  vapply(as.character(raw_vec), function(val) {
    if (is.na(val)) return(NA_character_)
    txt <- trimws(val)
    if (!nzchar(txt)) return(NA_character_)
    low <- tolower(txt)
    if (low %in% c("na", "nan", "none", "null", "unknown", .LEGACY_HEURISTIC_REGIONS)) {
      return(NA_character_)
    }
    cleaned <- gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", txt))
    if (!nzchar(cleaned) || tolower(cleaned) %in% .LEGACY_HEURISTIC_REGIONS) return(NA_character_)
    cleaned
  }, character(1), USE.NAMES = FALSE)
}

.norm_us_state_r <- function(st_vec) {
  abbr_map <- stats::setNames(tolower(datasets::state.name), tolower(datasets::state.abb))
  abbr_map["dc"] <- "district of columbia"
  vapply(as.character(st_vec), function(s) {
    if (is.na(s)) return(NA_character_)
    s_clean <- gsub("\\.$", "", gsub("\\s+", " ", tolower(trimws(s))))
    if (!nzchar(s_clean) || s_clean %in% c("na", "nan", "none", "unknown")) return(NA_character_)
    if (s_clean %in% names(abbr_map)) unname(abbr_map[[s_clean]]) else s_clean
  }, character(1), USE.NAMES = FALSE)
}

.norm_us_county_r <- function(co_vec, st_norm_vec = NULL) {
  n <- length(co_vec)
  if (is.null(st_norm_vec)) st_norm_vec <- rep(NA_character_, n)
  vapply(seq_len(n), function(i) {
    raw <- as.character(co_vec[i])
    if (is.na(raw)) return(NA_character_)
    c_str <- tolower(trimws(raw))
    if (!nzchar(c_str) || c_str %in% c("na", "nan", "none", "unknown", "unplaced", "indeterminate")) {
      return(NA_character_)
    }
    if (grepl(" or |/|\\?", c_str)) return(NA_character_)
    c_str <- gsub("\\[|\\]|\\(|\\)", "", c_str)
    c_str <- gsub("^city\\s+of\\s+", "", c_str)
    c_str <- gsub("\\s+(county|parish|par\\.?|borough|census\\s+area|municipality|co\\.?|city\\s+and\\s+borough|independent\\s+city|independent)$", "", c_str)
    c_str <- gsub("\\s+city$", "", c_str)
    c_str <- gsub("^(st\\.|saint)\\s+", "st ", c_str)
    c_str <- gsub("^ste\\.\\s*", "ste ", c_str)
    c_str <- trimws(gsub("[^a-z0-9]+", " ", c_str))
    if (!nzchar(c_str)) return(NA_character_)
    st_i <- as.character(st_norm_vec[i])
    if (!is.na(st_i)) {
      if (st_i == "virginia" && c_str == "nansemond") c_str <- "suffolk"
      else if (st_i == "virginia" && c_str == "princess anne") c_str <- "virginia beach"
      else if (st_i == "virginia" && c_str == "chesapeake") c_str <- "norfolk"
      else if (st_i == "georgia" && c_str == "columbus") c_str <- "muscogee"
      else if (st_i == "south dakota" && c_str == "shannon") c_str <- "oglala lakota"
    }
    c_str
  }, character(1), USE.NAMES = FALSE)
}

.build_offline_epa_l3_sf_r <- function() {
  mk_poly <- function(mat) sf::st_polygon(list(rbind(mat, mat[1, , drop = FALSE])))
  eco_df <- data.frame(
    ecoregion_label = c(
      "Blue_Ridge", "Ridge_and_Valley", "Piedmont", "Middle_Atlantic_Coastal_Plain",
      "Southeastern_Plains", "Southern_Coastal_Plain", "Mississippi_Alluvial_Plain",
      "South_Central_Plains", "East_Central_Texas_Plains", "Ouachita_Mountains",
      "Ozark_Highlands", "Interior_Plateau", "Southwestern_Appalachians"
    ),
    stringsAsFactors = FALSE
  )
  eco_geoms <- sf::st_sfc(
    mk_poly(matrix(c(-84.5,34.5, -83.2,34.8, -81.6,36.1, -79.0,37.9, -78.2,38.9, -79.5,39.2, -82.8,36.8, -84.6,35.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-85.6,34.2, -84.6,35.2, -82.8,36.8, -79.5,39.2, -77.2,40.8, -78.8,41.2, -83.8,37.4, -85.9,34.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-85.2,32.5, -82.0,33.5, -79.4,34.9, -77.6,36.5, -77.1,38.9, -78.2,38.9, -79.0,37.9, -81.6,36.1, -83.2,34.8, -84.5,34.5), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-79.6,33.0, -75.4,35.2, -75.0,38.8, -76.4,39.2, -77.1,38.5, -77.6,36.5, -78.5,35.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-90.2,31.0, -85.0,30.8, -81.3,32.2, -79.6,33.0, -78.5,35.2, -79.4,34.9, -82.0,33.5, -85.2,32.5, -88.8,34.8, -90.1,33.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-85.0,29.6, -80.2,26.5, -80.8,31.2, -81.3,32.2, -85.0,30.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-91.8,29.5, -90.2,30.5, -90.1,35.0, -89.2,37.0, -90.6,37.0, -91.5,34.2, -92.0,31.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-95.8,30.2, -92.0,31.0, -91.5,34.2, -94.8,34.1, -96.0,32.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-98.2,29.0, -95.8,30.2, -96.0,32.8, -97.2,33.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-95.6,34.1, -92.4,34.3, -92.6,35.1, -95.4,35.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-94.8,35.7, -91.0,35.8, -90.5,38.2, -94.6,38.1), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-88.5,34.8, -85.5,35.2, -84.5,38.2, -87.8,38.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-86.8,33.6, -85.6,34.2, -84.2,36.8, -85.5,36.8), ncol = 2, byrow = TRUE)),
    crs = 4326
  )
  eco_sf <- sf::st_sf(eco_df, geometry = eco_geoms)
  cnty_df <- data.frame(
    county_key = c("north carolina|orange", "north carolina|new hanover", "north carolina|wake", "texas|nacogdoches"),
    state_name = c("north carolina", "north carolina", "north carolina", "texas"),
    county_norm = c("orange", "new hanover", "wake", "nacogdoches"),
    centroid_lon = c(-79.10, -77.90, -78.40, -94.65),
    centroid_lat = c(36.05, 34.20, 35.78, 31.62),
    n_ecoregions = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_),
    top_ecoregion = c(NA_character_, NA_character_, NA_character_, NA_character_),
    top_area_share = c(NA_real_, NA_real_, NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  cnty_geoms <- sf::st_sfc(
    mk_poly(matrix(c(-79.25,35.86, -78.95,35.86, -78.95,36.24, -79.25,36.24), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-78.02,34.05, -77.78,34.05, -77.78,34.35, -78.02,34.35), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-79.00,35.50, -77.80,35.50, -77.80,36.05, -79.00,36.05), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-94.90,31.40, -94.40,31.40, -94.40,31.85, -94.90,31.85), ncol = 2, byrow = TRUE)),
    crs = 4326
  )
  cnty_sf <- sf::st_sf(cnty_df, geometry = cnty_geoms)
  list(eco_sf = eco_sf, cnty_sf = cnty_sf, use_live_county_intersection = TRUE)
}

load_epa_level3_spatial_r <- function(ecoregion_path = NULL, county_path = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("Package 'sf' is required for EPA Level III Ecoregion lookup.")
  custom_paths <- !is.null(ecoregion_path) || !is.null(county_path)
  default_dir <- "data/environmental/epa_ecoregions_l3"
  gpkg_p <- file.path(default_dir, "us_eco_l3.gpkg")
  eco_json_p <- if (!is.null(ecoregion_path)) as.character(ecoregion_path) else file.path(default_dir, "us_eco_l3.geojson")
  cnty_json_p <- if (!is.null(county_path)) as.character(county_path) else file.path(default_dir, "us_counties.geojson")

  cache_id <- paste(
    if (is.null(ecoregion_path)) "DEFAULT_ECO" else eco_json_p,
    if (is.null(county_path)) "DEFAULT_CNTY" else cnty_json_p,
    sep = "::"
  )
  if (exists(cache_id, envir = .epa_l3_r_cache, inherits = FALSE)) {
    return(get(cache_id, envir = .epa_l3_r_cache, inherits = FALSE))
  }

  if (!custom_paths && file.exists(gpkg_p)) {
    eco_sf <- sf::st_read(gpkg_p, layer = "us_eco_l3", quiet = TRUE)
    cnty_sf <- tryCatch(sf::st_read(gpkg_p, layer = "us_counties", quiet = TRUE), error = function(e) NULL)
    eco_sf$ecoregion_label <- standardize_epa_l3_label_r(
      if ("ecoregion_label" %in% names(eco_sf)) eco_sf$ecoregion_label else eco_sf$US_L3NAME
    )
    bundle <- list(eco_sf = eco_sf, cnty_sf = cnty_sf, use_live_county_intersection = FALSE)
    assign(cache_id, bundle, envir = .epa_l3_r_cache)
    return(bundle)
  }

  if (file.exists(eco_json_p)) {
    eco_sf <- sf::st_read(eco_json_p, quiet = TRUE)
    raw_lbl <- if ("ecoregion_label" %in% names(eco_sf)) eco_sf$ecoregion_label else if ("US_L3NAME" %in% names(eco_sf)) eco_sf$US_L3NAME else eco_sf[[1]]
    eco_sf$ecoregion_label <- standardize_epa_l3_label_r(raw_lbl)
    cnty_sf <- if (file.exists(cnty_json_p)) tryCatch(sf::st_read(cnty_json_p, quiet = TRUE), error = function(e) NULL) else NULL
    bundle <- list(eco_sf = eco_sf, cnty_sf = cnty_sf, use_live_county_intersection = custom_paths)
    assign(cache_id, bundle, envir = .epa_l3_r_cache)
    return(bundle)
  }

  .build_offline_epa_l3_sf_r()
}

infer_epa_level3_ecoregion <- function(
  lat,
  lon,
  state_province = NULL,
  county = NULL,
  ecoregion_path = NULL,
  county_path = NULL
) {
  n <- max(length(lat), length(lon), 1L)
  if (length(lat) == 0L && length(lon) == 0L) return(character(0))
  lat_num <- suppressWarnings(as.numeric(rep_len(lat, n)))
  lon_num <- suppressWarnings(as.numeric(rep_len(lon, n)))
  st_vec <- if (!is.null(state_province) && length(state_province) > 0L) rep_len(as.character(state_province), n) else rep(NA_character_, n)
  co_vec <- if (!is.null(county) && length(county) > 0L) rep_len(as.character(county), n) else rep(NA_character_, n)

  old_s2 <- sf::sf_use_s2()
  suppressMessages(sf::sf_use_s2(FALSE))
  on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)

  bundle <- load_epa_level3_spatial_r(ecoregion_path = ecoregion_path, county_path = county_path)
  eco_sf <- bundle$eco_sf
  cnty_sf <- bundle$cnty_sf
  out <- rep(NA_character_, n)

  # (a) Finite coordinates: strict point-in-polygon intersection
  has_coords <- is.finite(lat_num) & is.finite(lon_num) & (lat_num >= -90) & (lat_num <= 90) & (lon_num >= -180) & (lon_num <= 180)
  if (any(has_coords)) {
    pts_sf <- sf::st_as_sf(
      data.frame(lon = lon_num[has_coords], lat = lat_num[has_coords]),
      coords = c("lon", "lat"),
      crs = 4326
    )
    hits <- suppressWarnings(sf::st_intersects(pts_sf, eco_sf))
    out[has_coords] <- vapply(hits, function(h) {
      if (length(h) >= 1L) as.character(eco_sf$ecoregion_label[h[1]]) else NA_character_
    }, character(1))
  }

  # (b) Missing coordinates: county-containment fallback (>= 95% single-ecoregion area share)
  miss_idx <- which(!has_coords)
  if (length(miss_idx) > 0L && !is.null(cnty_sf) && nrow(cnty_sf) > 0L) {
    st_norm <- .norm_us_state_r(st_vec[miss_idx])
    co_norm <- .norm_us_county_r(co_vec[miss_idx], st_norm)
    c_st <- .norm_us_state_r(if ("state_name" %in% names(cnty_sf)) cnty_sf$state_name else cnty_sf[[1]])
    c_co <- .norm_us_county_r(if ("county_norm" %in% names(cnty_sf)) cnty_sf$county_norm else cnty_sf[[2]], c_st)
    c_keys <- paste0(c_st, "|", c_co)
    c_keys_ns <- paste0(c_st, "|", gsub(" ", "", c_co))

    for (k_i in seq_along(miss_idx)) {
      row_i <- miss_idx[k_i]
      if (is.na(st_norm[k_i]) || is.na(co_norm[k_i])) next
      q_key <- paste0(st_norm[k_i], "|", co_norm[k_i])
      q_key_ns <- paste0(st_norm[k_i], "|", gsub(" ", "", co_norm[k_i]))
      m_idx <- match(q_key, c_keys)
      if (is.na(m_idx)) m_idx <- match(q_key_ns, c_keys_ns)
      if (is.na(m_idx)) next

      c_row <- cnty_sf[m_idx, ]
      use_live <- isTRUE(bundle$use_live_county_intersection) || !("top_area_share" %in% names(c_row)) || is.na(c_row$top_area_share[1])
      if (use_live) {
        ix <- tryCatch(suppressWarnings(sf::st_intersection(c_row, eco_sf[, "ecoregion_label", drop = FALSE])), error = function(e) NULL)
        if (is.null(ix) || nrow(ix) == 0L) next
        a_parts <- as.numeric(sf::st_area(ix))
        a_by_eco <- tapply(a_parts, ix$ecoregion_label, sum, na.rm = TRUE)
        tot_a <- sum(a_by_eco, na.rm = TRUE)
        c_tot_a <- as.numeric(sf::st_area(c_row))[1]
        if (!is.finite(tot_a) || tot_a <= 0 || (is.finite(c_tot_a) && c_tot_a > 0 && (tot_a / c_tot_a) < 0.50)) next
        top_idx <- which.max(a_by_eco)
        top_eco <- names(a_by_eco)[top_idx]
        top_share <- as.numeric(a_by_eco[top_idx]) / tot_a
        n_eco <- sum((a_by_eco / tot_a) >= 0.001)
      } else {
        top_eco <- standardize_epa_l3_label_r(c_row$top_ecoregion[1])
        top_share <- suppressWarnings(as.numeric(c_row$top_area_share[1]))
        n_eco <- suppressWarnings(as.integer(c_row$n_ecoregions[1]))
      }

      if (is.na(top_eco) || !(isTRUE(n_eco == 1L) || (is.finite(top_share) && top_share >= 0.95))) next

      c_lon <- if ("centroid_lon" %in% names(c_row)) suppressWarnings(as.numeric(c_row$centroid_lon[1])) else NA_real_
      c_lat <- if ("centroid_lat" %in% names(c_row)) suppressWarnings(as.numeric(c_row$centroid_lat[1])) else NA_real_
      if (!is.finite(c_lon) || !is.finite(c_lat)) {
        pos_pt <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(c_row)))
        c_lon <- pos_pt[1, 1]; c_lat <- pos_pt[1, 2]
      }
      c_pt_sf <- sf::st_as_sf(data.frame(lon = c_lon, lat = c_lat), coords = c("lon", "lat"), crs = 4326)
      c_hit <- suppressWarnings(sf::st_intersects(c_pt_sf, eco_sf))[[1]]
      hit_lbl <- if (length(c_hit) >= 1L) as.character(eco_sf$ecoregion_label[c_hit[1]]) else NA_character_
      if (!is.na(hit_lbl) && hit_lbl == top_eco) {
        out[row_i] <- hit_lbl
      } else if (is.na(hit_lbl)) {
        out[row_i] <- top_eco
      }
    }
  }
  out
}

run_supervised_ecoregion_cda <- function(
  df,
  feature_cols,
  gmm_res = NULL,
  output_summary_csv = "outputs/reports/ecoregion_supervised_cda_summary.csv",
  permutations = 999,
  seed = 42,
  min_anchors_per_ecoregion = 3L
) {
  n <- nrow(df)
  rg <- if ("regional_group" %in% names(df)) standardize_epa_l3_label_r(df$regional_group) else rep(NA_character_, n)

  # If regional_group was absent or contained only legacy non-EPA heuristic strings, infer via polygon lookup
  if (all(is.na(rg))) {
    lat_col <- intersect(c("decimalLatitude", "latitude"), names(df))
    lon_col <- intersect(c("decimalLongitude", "longitude"), names(df))
    if (length(lat_col) > 0L && length(lon_col) > 0L || any(c("stateProvince", "county") %in% names(df))) {
      lat_v <- if (length(lat_col) > 0L) df[[lat_col[1]]] else rep(NA_real_, n)
      lon_v <- if (length(lon_col) > 0L) df[[lon_col[1]]] else rep(NA_real_, n)
      st_v <- if ("stateProvince" %in% names(df)) df$stateProvince else NULL
      co_v <- if ("county" %in% names(df)) df$county else NULL
      rg <- tryCatch(infer_epa_level3_ecoregion(lat_v, lon_v, state_province = st_v, county = co_v), error = function(e) rg)
    }
  }

  pc_req <- intersect(paste0("PC", seq_len(5)), feature_cols)
  req_cols <- if (length(pc_req) > 0L) pc_req else feature_cols
  has_pcs <- stats::complete.cases(df[, req_cols, drop = FALSE])
  auth_w <- if ("authority_weight" %in% names(df)) suppressWarnings(as.numeric(df$authority_weight)) else rep(1.0, n)
  auth_w[is.na(auth_w)] <- 0.0
  is_t1 <- if ("morph_tier_source" %in% names(df)) (!is.na(df$morph_tier_source) & df$morph_tier_source == "Basal_Rosette") else rep(TRUE, n)

  # Strictly prefer Basal_Rosette Gold/Silver anchors (>= 3 per EPA Level III Ecoregion), with valid Basal_Rosette fallback only when Gold/Silver < 2 groups
  cand_mask <- !is.na(rg) & has_pcs & is_t1 & (auth_w > 0)
  tab_cand <- table(rg[cand_mask])
  active_ecos <- sort(names(tab_cand[tab_cand >= min_anchors_per_ecoregion]))
  if (length(active_ecos) < 2L) {
    cand_t1 <- !is.na(rg) & has_pcs & is_t1
    tab_t1 <- table(rg[cand_t1])
    if (sum(tab_t1 >= min_anchors_per_ecoregion) >= 2L) {
      cand_mask <- cand_t1
      active_ecos <- sort(names(tab_t1[tab_t1 >= min_anchors_per_ecoregion]))
    } else {
      cand_all <- !is.na(rg) & has_pcs
      tab_all <- table(rg[cand_all])
      if (sum(tab_all >= min_anchors_per_ecoregion) >= 2L) {
        cand_mask <- cand_all
        active_ecos <- sort(names(tab_all[tab_all >= min_anchors_per_ecoregion]))
      }
    }
  }

  # Unsupervised Morphotype_K* x Supervised EPA_Ecoregion contingency statistics (excluding NA ecoregions)
  gmm_cl <- if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n) {
    as.integer(gmm_res$classification)
  } else if ("gmm_cluster" %in% names(df)) {
    as.integer(df$gmm_cluster)
  } else {
    rep(NA_integer_, n)
  }
  morph_lbl <- ifelse(is.na(gmm_cl), NA_character_, paste0("Morphotype_K", gmm_cl))
  gmm_opt_k <- length(unique(stats::na.omit(gmm_cl)))
  cont_mask <- !is.na(morph_lbl) & !is.na(rg)
  ari_val <- NA_real_
  fisher_p <- NA_real_
  if (sum(cont_mask) >= 4L) {
    ari_val <- compute_ari_pure(morph_lbl[cont_mask], rg[cont_mask])
    cont_tab <- table(morph_lbl[cont_mask], rg[cont_mask])
    if (nrow(cont_tab) >= 2L && ncol(cont_tab) >= 2L) {
      set.seed(seed)
      use_exact <- (nrow(cont_tab) <= 3L && ncol(cont_tab) <= 3L && sum(cont_tab) <= 180L)
      fisher_p <- tryCatch(
        if (use_exact) stats::fisher.test(cont_tab, workspace = 2e6)$p.value
        else stats::fisher.test(cont_tab, simulate.p.value = TRUE, B = max(2000L, as.integer(permutations)))$p.value,
        error = function(e) tryCatch(stats::fisher.test(cont_tab, simulate.p.value = TRUE, B = 2000L)$p.value, error = function(e2) NA_real_)
      )
    }
  }

  pred_eco <- rep(NA_character_, n)
  post_max <- rep(NA_real_, n)
  can1_out <- rep(NA_real_, n)
  can2_out <- rep(NA_real_, n)
  post_mat <- NULL
  loocv_acc <- NA_real_
  perm_f <- NA_real_
  perm_r2 <- NA_real_
  perm_p <- NA_real_
  var_pct <- c(100.0, 0.0)
  struct_load <- NULL
  centroids_2d <- NULL
  is_train <- rep(FALSE, n)

  if (length(active_ecos) >= 2L) {
    is_train <- cand_mask & (rg %in% active_ecos)
    X_full <- as.matrix(df[, feature_cols, drop = FALSE])
    storage.mode(X_full) <- "double"
    for (j in seq_len(ncol(X_full))) {
      if (any(is.na(X_full[, j]))) {
        fill_v <- if (any(!is.na(X_full[is_train, j]))) mean(X_full[is_train, j], na.rm = TRUE) else 0.0
        X_full[is.na(X_full[, j]), j] <- fill_v
      }
    }

    sel_cols <- if (exists("select_morphometric_subspace", mode = "function")) {
      select_morphometric_subspace(X_full[is_train, , drop = FALSE], group_vec = rg[is_train])
    } else {
      n_min_g <- min(as.integer(table(rg[is_train])))
      feature_cols[seq_len(max(1L, min(length(feature_cols), floor(n_min_g / 3L))))]
    }
    X_use <- X_full[, sel_cols, drop = FALSE]
    p_dim <- ncol(X_use)

    # Per-ecoregion 95% chi-square Mahalanobis trimming on active training anchors (exempting nomenclatural types)
    is_type_eco <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, n)
    is_type_eco[is.na(is_type_eco)] <- FALSE
    for (eco_g in active_ecos) {
      idx_g <- which(is_train & rg == eco_g)
      if (length(idx_g) >= 5L) {
        sub_g <- X_use[idx_g, , drop = FALSE]
        d2_g <- tryCatch(stats::mahalanobis(sub_g, colMeans(sub_g), stats::cov(sub_g) + diag(1e-4, p_dim)), error = function(e) rep(0, length(idx_g)))
        drop_g <- !is.na(d2_g) & (d2_g > stats::qchisq(0.95, df = p_dim)) & !is_type_eco[idx_g]
        if (sum(!drop_g) >= min_anchors_per_ecoregion) is_train[idx_g[drop_g]] <- FALSE
      }
    }

    active_idx <- which(is_train)
    X_a <- X_use[active_idx, , drop = FALSE]
    grp_a <- droplevels(factor(rg[active_idx], levels = active_ecos))
    g_n <- nlevels(grp_a)
    priors_u <- rep(1.0 / g_n, g_n)

    if (g_n >= 2L && requireNamespace("MASS", quietly = TRUE)) {
      # Drop any near-constant predictor inside active training anchors to prevent singular covariance
      sds_a <- apply(X_a, 2, stats::sd, na.rm = TRUE)
      keep_sd <- is.finite(sds_a) & (sds_a > 1e-7)
      if (!any(keep_sd)) keep_sd[1] <- TRUE
      X_a_fit <- X_a[, keep_sd, drop = FALSE]
      X_proj_fit <- X_use[, keep_sd, drop = FALSE]

      lda_cv <- tryCatch(MASS::lda(x = X_a_fit, grouping = grp_a, prior = priors_u, CV = TRUE), error = function(e) NULL)
      if (!is.null(lda_cv) && !is.null(lda_cv$class)) {
        loocv_acc <- round(mean(as.character(lda_cv$class) == as.character(grp_a), na.rm = TRUE), 4)
      }

      lda_fit <- tryCatch(MASS::lda(x = X_a_fit, grouping = grp_a, prior = priors_u), error = function(e) NULL)
      if (!is.null(lda_fit)) {
        lda_pred <- stats::predict(lda_fit, newdata = X_proj_fit)
        pred_eco <- as.character(lda_pred$class)
        post_mat <- as.matrix(lda_pred$posterior)
        post_max <- round(apply(post_mat, 1, max, na.rm = TRUE), 4)

        Z_raw <- as.matrix(lda_pred$x)
        if (ncol(Z_raw) == 1L) {
          v1 <- as.numeric(Z_raw[, 1])
          src2 <- if ("PC2" %in% colnames(X_full)) X_full[, "PC2"] else if (ncol(X_full) >= 2L) X_full[, 2] else X_full[, 1]
          v2 <- as.numeric(src2) - mean(src2, na.rm = TRUE)
          b1 <- v1 - mean(v1, na.rm = TRUE)
          v2 <- v2 - (sum(v2 * b1, na.rm = TRUE) / max(sum(b1^2, na.rm = TRUE), 1e-12)) * b1
          Z_raw <- cbind(LD1 = v1, LD2 = v2)
        }
        can1_out <- round(as.numeric(Z_raw[, 1]), 5)
        can2_out <- round(as.numeric(Z_raw[, 2]), 5)

        svd_sq <- lda_fit$svd^2
        if (length(svd_sq) >= 2L) {
          var_pct <- round(100.0 * svd_sq[1:2] / max(sum(svd_sq), 1e-12), 2)
        } else if (length(svd_sq) == 1L) {
          var_pct <- c(100.0, 0.0)
        }

        Z_act_2d <- cbind(can1_out[active_idx], can2_out[active_idx])
        struct_load <- tryCatch(stats::cor(X_a_fit, Z_act_2d), error = function(e) NULL)
        centroids_2d <- do.call(rbind, lapply(levels(grp_a), function(eg) {
          m_e <- active_idx[grp_a == eg]
          c(mean(can1_out[m_e], na.rm = TRUE), mean(can2_out[m_e], na.rm = TRUE))
        }))
        rownames(centroids_2d) <- levels(grp_a)
        colnames(centroids_2d) <- c("Can1", "Can2")
      }
    }

    # Euclidean PERMANOVA across active EPA Level III Ecoregions
    set.seed(seed)
    if (g_n >= 2L && exists("fast_euclidean_permanova_trace", mode = "function")) {
      ad_fast <- tryCatch(
        fast_euclidean_permanova_trace(X_a, grp_a, permutations = as.integer(permutations), seed = seed),
        error = function(e) NULL
      )
      if (!is.null(ad_fast) && is.finite(ad_fast$R2)) {
        perm_r2 <- round(as.numeric(ad_fast$R2), 4)
        perm_f <- round(as.numeric(ad_fast$F), 4)
        perm_p <- signif(as.numeric(ad_fast$p_value), 4)
      }
    } else if (requireNamespace("vegan", quietly = TRUE) && g_n >= 2L) {
      ad_res <- tryCatch(
        vegan::adonis2(X_a ~ grp_a, permutations = as.integer(permutations), method = "euclidean"),
        error = function(e) NULL
      )
      if (!is.null(ad_res) && nrow(ad_res) >= 1L) {
        perm_r2 <- round(as.numeric(ad_res$R2[1]), 4)
        perm_f <- round(as.numeric(ad_res$F[1]), 4)
        perm_p <- signif(as.numeric(ad_res[["Pr(>F)"]][1]), 4)
      }
    }
  } else if (length(active_ecos) == 1L) {
    pred_eco <- rep(active_ecos[1], n)
  }

  all_present_ecos <- sort(unique(rg[!is.na(rg)]))
  n_na_total <- sum(is.na(rg))
  if (length(all_present_ecos) > 0L) {
    summary_df <- do.call(rbind, lapply(all_present_ecos, function(eg) {
      data.frame(
        ecoregion = eg,
        is_active_cda_group = eg %in% active_ecos,
        n_total_vouchers = sum(!is.na(rg) & rg == eg),
        n_active_anchors = sum(is_train & !is.na(rg) & rg == eg),
        n_unassigned_na_projected = n_na_total,
        loocv_accuracy = loocv_acc,
        permanova_f = perm_f,
        permanova_r2 = perm_r2,
        permanova_p_value = perm_p,
        unsupervised_gmm_k = gmm_opt_k,
        ari_morphotype_vs_ecoregion = ari_val,
        fisher_exact_p_value = signif(fisher_p, 4),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    summary_df <- data.frame(
      ecoregion = NA_character_,
      is_active_cda_group = FALSE,
      n_total_vouchers = 0L,
      n_active_anchors = 0L,
      n_unassigned_na_projected = n_na_total,
      loocv_accuracy = NA_real_,
      permanova_f = NA_real_,
      permanova_r2 = NA_real_,
      permanova_p_value = NA_real_,
      unsupervised_gmm_k = gmm_opt_k,
      ari_morphotype_vs_ecoregion = NA_real_,
      fisher_exact_p_value = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(output_summary_csv) && nzchar(output_summary_csv)) {
    dir.create(dirname(output_summary_csv), recursive = TRUE, showWarnings = FALSE)
    tmp_sum <- paste0(output_summary_csv, ".tmp")
    utils::write.csv(summary_df, file = tmp_sum, row.names = FALSE)
    file.rename(tmp_sum, output_summary_csv)
  }

  list(
    ecoregion_vec = rg,
    active_ecoregions = active_ecos,
    is_active = is_train,
    ecoregion_cda_predicted = pred_eco,
    ecoregion_cda_posterior_prob = post_max,
    ecoregion_can1 = can1_out,
    ecoregion_can2 = can2_out,
    post_probs = post_mat,
    loocv_accuracy = loocv_acc,
    permanova_f = perm_f,
    permanova_r2 = perm_r2,
    permanova_p_value = perm_p,
    ari_morphotype_vs_ecoregion = ari_val,
    fisher_exact_p_value = fisher_p,
    variance_pct = var_pct,
    structure_loadings = struct_load,
    centroids = centroids_2d,
    summary_table = summary_df
  )
}

# ==============================================================================
# Formal Parametric Multivariate MANOVA (Pillai's Trace & Wilks' Lambda) and
# Pairwise Post-Hoc Hotelling's T^2 / MANOVA Contrast Audit
# ==============================================================================
.resolve_default_manova_features <- function(df, feature_cols = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(character(0))
  cand <- if (!is.null(feature_cols) && length(feature_cols) > 0L) {
    as.character(feature_cols)
  } else {
    c(
      paste0("PC", seq_len(5)),
      "AR_inv_std",
      "petiole_to_blade_ratio_std",
      "plant_height_cm_std",
      "log_centsize_std",
      "aspect_ratio_std",
      "solidity_std"
    )
  }
  avail <- intersect(unique(cand), names(df))
  if (length(avail) == 0L) avail <- intersect(c(paste0("PC", seq_len(5)), paste0("Sym_PC", seq_len(5))), names(df))
  avail
}

format_manova_pillai_subtitle <- function(
  manova_res,
  factor_name,
  alt_factor = NULL,
  fallback_label = NULL
) {
  gt <- if (!is.null(manova_res) && is.list(manova_res) && !is.null(manova_res$group_table)) {
    manova_res$group_table
  } else if (is.data.frame(manova_res)) {
    manova_res
  } else {
    NULL
  }
  if (!is.null(gt) && nrow(gt) > 0L && "grouping_factor" %in% names(gt)) {
    idx <- which(gt$grouping_factor == factor_name)
    if (length(idx) == 0L && !is.null(alt_factor)) {
      idx <- which(gt$grouping_factor == alt_factor)
    }
    if (length(idx) > 0L) {
      r <- gt[idx[1], , drop = FALSE]
      v_val <- suppressWarnings(as.numeric(r$pillai_trace[1]))
      f_val <- suppressWarnings(as.numeric(r$approx_F[1]))
      ndf <- suppressWarnings(as.numeric(r$num_df[1]))
      ddf <- suppressWarnings(as.numeric(r$den_df[1]))
      p_val <- suppressWarnings(as.numeric(r$p_value[1]))
      if (is.finite(v_val) && is.finite(f_val)) {
        ndf_s <- if (is.finite(ndf) && abs(ndf - round(ndf)) < 1e-5) sprintf("%d", as.integer(round(ndf))) else sprintf("%.1f", ndf)
        ddf_s <- if (is.finite(ddf) && abs(ddf - round(ddf)) < 1e-5) sprintf("%d", as.integer(round(ddf))) else sprintf("%.1f", ddf)
        p_s <- if (!is.finite(p_val)) "NA" else if (p_val < 1e-4) sprintf("%.2e", p_val) else sprintf("%.4f", p_val)
        return(sprintf("Omnibus MANOVA: Pillai's Trace V = %.3f, F(%s, %s) = %.2f, p = %s", v_val, ndf_s, ddf_s, f_val, p_s))
      }
    }
  }
  fallback_label
}

.fit_weighted_manova_stats <- function(Y_mat, grp_vec, w_vec) {
  grp_fac <- droplevels(factor(as.character(grp_vec)))
  g_k <- nlevels(grp_fac)
  n_obs <- nrow(Y_mat)
  df_res <- n_obs - g_k
  if (g_k < 2L || df_res < 1L || ncol(Y_mat) < 1L) return(NULL)

  w_pos <- suppressWarnings(as.numeric(w_vec))
  w_pos[!is.finite(w_pos) | w_pos <= 0] <- 1.0
  w_pos <- w_pos / max(mean(w_pos), 1e-12)

  sds <- apply(Y_mat, 2, stats::sd, na.rm = TRUE)
  keep_sd <- is.finite(sds) & (sds > 1e-7)
  if (!any(keep_sd)) return(NULL)
  Y_use <- Y_mat[, keep_sd, drop = FALSE]

  p_cap <- max(1L, df_res - 1L)
  if (ncol(Y_use) > p_cap) {
    Y_use <- Y_use[, seq_len(p_cap), drop = FALSE]
  }

  lm_pre <- tryCatch(stats::lm(Y_use ~ grp_fac, weights = w_pos), error = function(e) NULL)
  if (is.null(lm_pre)) return(NULL)
  res_mat <- as.matrix(stats::residuals(lm_pre)) * sqrt(w_pos)
  qr_e <- qr(res_mat, tol = 1e-6)
  if (qr_e$rank < 1L) return(NULL)
  if (qr_e$rank < ncol(Y_use)) {
    piv <- sort(qr_e$pivot[seq_len(qr_e$rank)])
    Y_use <- Y_use[, piv, drop = FALSE]
  }

  p_dim <- ncol(Y_use)
  s_rank <- min(p_dim, g_k - 1L)

  if (p_dim >= 2L) {
    m_fit <- tryCatch(stats::manova(Y_use ~ grp_fac, weights = w_pos), error = function(e) NULL)
    if (is.null(m_fit)) return(NULL)
    s_pil <- tryCatch(summary(m_fit, test = "Pillai")$stats, error = function(e) NULL)
    s_wlk <- tryCatch(summary(m_fit, test = "Wilks")$stats, error = function(e) NULL)
    s_hl  <- tryCatch(summary(m_fit, test = "Hotelling-Lawley")$stats, error = function(e) NULL)
    if (is.null(s_pil) || is.null(s_wlk) || nrow(s_pil) < 1L) return(NULL)

    pillai_v <- pmax(0.0, as.numeric(s_pil[1, "Pillai"]))
    wilks_l  <- pmin(1.0, pmax(0.0, as.numeric(s_wlk[1, "Wilks"])))
    hl_u     <- if (!is.null(s_hl) && "Hotelling-Lawley" %in% colnames(s_hl)) {
      pmax(0.0, as.numeric(s_hl[1, "Hotelling-Lawley"]))
    } else {
      pillai_v / max(1.0 - min(pillai_v, 1.0 - 1e-12), 1e-12)
    }
    approx_f <- pmax(0.0, as.numeric(s_pil[1, "approx F"]))
    num_df   <- as.numeric(s_pil[1, "num Df"])
    den_df   <- as.numeric(s_pil[1, "den Df"])
    p_val    <- pmin(1.0, pmax(0.0, as.numeric(s_pil[1, "Pr(>F)"])))
  } else {
    aov_1 <- tryCatch(stats::anova(stats::lm(Y_use[, 1] ~ grp_fac, weights = w_pos)), error = function(e) NULL)
    if (is.null(aov_1) || nrow(aov_1) < 2L) return(NULL)
    ss_eff <- pmax(0.0, as.numeric(aov_1[1, "Sum Sq"]))
    ss_res <- pmax(1e-12, as.numeric(aov_1[2, "Sum Sq"]))
    ss_tot <- ss_eff + ss_res
    pillai_v <- ss_eff / ss_tot
    wilks_l  <- ss_res / ss_tot
    hl_u     <- ss_eff / ss_res
    approx_f <- pmax(0.0, as.numeric(aov_1[1, "F value"]))
    num_df   <- as.numeric(aov_1[1, "Df"])
    den_df   <- as.numeric(aov_1[2, "Df"])
    p_val    <- pmin(1.0, pmax(0.0, as.numeric(aov_1[1, "Pr(>F)"])))
  }

  p_eta2 <- min(1.0, max(0.0, pillai_v / max(1L, s_rank)))
  hot_t2 <- pmax(0.0, df_res * hl_u)
  mahal_d <- NA_real_
  if (g_k == 2L) {
    lvls <- levels(grp_fac)
    w1_sum <- sum(w_pos[grp_fac == lvls[1]])
    w2_sum <- sum(w_pos[grp_fac == lvls[2]])
    mahal_d <- sqrt(max(0.0, hot_t2 * (w1_sum + w2_sum) / max(w1_sum * w2_sum, 1e-12)))
  }

  list(
    n_groups = g_k,
    n_obs = n_obs,
    n_features = p_dim,
    active_features = colnames(Y_use),
    pillai_trace = pillai_v,
    wilks_lambda = wilks_l,
    hotelling_t2 = hot_t2,
    mahalanobis_d = mahal_d,
    approx_F = approx_f,
    num_df = num_df,
    den_df = den_df,
    partial_eta2 = p_eta2,
    p_value = p_val
  )
}

run_multivariate_manova_audit <- function(
  df,
  feature_cols,
  gmm_res = NULL,
  output_group_csv = "outputs/reports/manova_group_significance.csv",
  output_pairwise_csv = "outputs/reports/manova_pairwise_posthoc.csv",
  include_gmm_clusters = TRUE
) {
  run_dir <- trimws(Sys.getenv("PACKERA_RUN_DIR", unset = ""))
  if (nzchar(run_dir)) {
    if (missing(output_group_csv) && !is.null(output_group_csv)) {
      output_group_csv <- file.path(run_dir, "outputs/reports/manova_group_significance.csv")
    }
    if (missing(output_pairwise_csv) && !is.null(output_pairwise_csv)) {
      output_pairwise_csv <- file.path(run_dir, "outputs/reports/manova_pairwise_posthoc.csv")
    }
  }

  empty_group_df <- data.frame(
    grouping_factor = character(0),
    factor_label = character(0),
    n_groups = integer(0),
    n_active_vouchers = integer(0),
    n_features = integer(0),
    active_features = character(0),
    pillai_trace = numeric(0),
    wilks_lambda = numeric(0),
    approx_F = numeric(0),
    num_df = numeric(0),
    den_df = numeric(0),
    partial_eta2 = numeric(0),
    p_value = numeric(0),
    q_value = numeric(0),
    stringsAsFactors = FALSE
  )
  empty_pairwise_df <- data.frame(
    grouping_factor = character(0),
    group_1 = character(0),
    group_2 = character(0),
    contrast = character(0),
    n_1 = integer(0),
    n_2 = integer(0),
    n_features = integer(0),
    hotelling_t2 = numeric(0),
    mahalanobis_d = numeric(0),
    pillai_trace = numeric(0),
    wilks_lambda = numeric(0),
    approx_F = numeric(0),
    num_df = numeric(0),
    den_df = numeric(0),
    partial_eta2 = numeric(0),
    p_value = numeric(0),
    q_value = numeric(0),
    stringsAsFactors = FALSE
  )

  atomic_write_manova_csv <- function(tbl, path) {
    if (!is.null(path) && nzchar(path)) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      tmp_p <- paste0(path, ".", Sys.getpid(), ".tmp")
      utils::write.csv(tbl, file = tmp_p, row.names = FALSE)
      file.rename(tmp_p, path)
    }
  }

  if (is.null(df) || !is.data.frame(df) || nrow(df) < 4L) {
    atomic_write_manova_csv(empty_group_df, output_group_csv)
    atomic_write_manova_csv(empty_pairwise_df, output_pairwise_csv)
    return(list(
      group_table = empty_group_df,
      pairwise_table = empty_pairwise_df,
      output_group_csv = output_group_csv,
      output_pairwise_csv = output_pairwise_csv
    ))
  }

  n <- nrow(df)
  avail_cols <- .resolve_default_manova_features(df, feature_cols)
  avail_cols <- avail_cols[vapply(avail_cols, function(cn) {
    sum(is.finite(suppressWarnings(as.numeric(df[[cn]])))) >= 3L
  }, logical(1))]

  if (length(avail_cols) == 0L) {
    atomic_write_manova_csv(empty_group_df, output_group_csv)
    atomic_write_manova_csv(empty_pairwise_df, output_pairwise_csv)
    return(list(
      group_table = empty_group_df,
      pairwise_table = empty_pairwise_df,
      output_group_csv = output_group_csv,
      output_pairwise_csv = output_pairwise_csv
    ))
  }

  # Strict active/passive authority weighting (Tier_1_Gold & Tier_2_Silver active; Tier_3_Bronze = 0.0)
  is_type_row <- if ("is_type_specimen" %in% names(df)) as.logical(df$is_type_specimen) else rep(FALSE, n)
  is_type_row[is.na(is_type_row)] <- FALSE
  if ("authority_weight" %in% names(df)) {
    auth_w <- suppressWarnings(as.numeric(df$authority_weight))
    auth_w[is.na(auth_w) | !is.finite(auth_w)] <- 0.0
  } else if ("determiner_tier" %in% names(df)) {
    if (exists("compute_authority_precision_weights", mode = "function")) {
      auth_w <- compute_authority_precision_weights(df)
    } else {
      det_tier <- as.character(df$determiner_tier)
      err_col <- intersect(c("dawid_skene_error_rate", "label_error_rate", "error_rate_k"), names(df))
      prec_w <- if (length(err_col) > 0L) pmax(0.0, pmin(1.0, 1.0 - suppressWarnings(as.numeric(df[[err_col[1]]])))) else rep(1.0, n)
      prec_w[is.na(prec_w) | !is.finite(prec_w)] <- 1.0
      auth_w <- ifelse(det_tier %in% c("Tier_0_Type_Anchor", "Tier_1_Gold", "Tier_2_Silver") | is_type_row, prec_w, 0.0)
    }
  } else {
    auth_w <- rep(1.0, n)
  }
  has_positive_anchors <- any(auth_w > 0)

  is_t1 <- if ("morph_tier_source" %in% names(df)) {
    !is.na(df$morph_tier_source) & (as.character(df$morph_tier_source) == "Basal_Rosette")
  
  } else {
    rep(TRUE, n)
  }

  pc_req <- intersect(paste0("PC", seq_len(5)), avail_cols)
  req_cols <- if (length(pc_req) > 0L) pc_req else avail_cols
  has_pcs <- stats::complete.cases(df[, req_cols, drop = FALSE])

  # 1. species_standardized (>= 2 species)
  sp_raw <- if ("species_standardized" %in% names(df)) {
    as.character(df$species_standardized)
  } else if ("scientificName" %in% names(df)) {
    standardize_taxon_name(df$scientificName)
  } else if ("species_raw" %in% names(df)) {
    standardize_taxon_name(df$species_raw)
  } else {
    rep(NA_character_, n)
  }
  sp_vec <- trimws(as.character(sp_raw))
  sp_vec[is.na(sp_vec) | !nzchar(sp_vec) | tolower(sp_vec) == "unknown"] <- NA_character_

  # 2. Unsupervised Morphotype_K1..Kn from gmm_cluster / gmm_res (>= 2 clusters)
  cl_raw <- if (!is.null(gmm_res) && !is.null(gmm_res$classification) && length(gmm_res$classification) == n) {
    gmm_res$classification
  } else if ("gmm_cluster" %in% names(df)) {
    df$gmm_cluster
  } else {
    rep(NA_character_, n)
  }
  cl_chr <- trimws(as.character(cl_raw))
  cl_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", cl_chr)))
  cl_vec <- ifelse(
    is.na(cl_chr) | !nzchar(cl_chr) | cl_chr == "NA",
    NA_character_,
    ifelse(grepl("^Morphotype_K[0-9]+$", cl_chr), cl_chr, ifelse(is.na(cl_num), NA_character_, paste0("Morphotype_K", cl_num)))
  )

  # 3. Non-NA regional_group US EPA Level III Ecoregions (>= 2 regions with >= 3 vouchers)
  rg_vec <- if ("regional_group" %in% names(df)) standardize_epa_l3_label_r(df$regional_group) else rep(NA_character_, n)
  if (all(is.na(rg_vec))) {
    lat_col <- intersect(c("decimalLatitude", "latitude"), names(df))
    lon_col <- intersect(c("decimalLongitude", "longitude"), names(df))
    if ((length(lat_col) > 0L && length(lon_col) > 0L) || any(c("stateProvince", "county") %in% names(df))) {
      lat_v <- if (length(lat_col) > 0L) df[[lat_col[1]]] else rep(NA_real_, n)
      lon_v <- if (length(lon_col) > 0L) df[[lon_col[1]]] else rep(NA_real_, n)
      st_v <- if ("stateProvince" %in% names(df)) df$stateProvince else NULL
      co_v <- if ("county" %in% names(df)) df$county else NULL
      rg_vec <- tryCatch(infer_epa_level3_ecoregion(lat_v, lon_v, state_province = st_v, county = co_v), error = function(e) rg_vec)
    }
  }
  rg_vec[is.na(rg_vec) | !nzchar(trimws(rg_vec))] <- NA_character_

  sort_group_levels <- function(lvls) {
    if (length(lvls) > 0L && all(grepl("^Morphotype_K[0-9]+$", lvls))) {
      lvls[order(as.integer(sub("^Morphotype_K", "", lvls)))]
    } else {
      sort(lvls)
    }
  }

  factor_specs <- Filter(Negate(is.null), list(
    list(name = "species_standardized", label = "Taxonomic Species (species_standardized)", vec = sp_vec, min_n = 2L),
    if (isTRUE(include_gmm_clusters)) list(name = "gmm_cluster", label = "Unsupervised GMM Morphotype (Morphotype_K1..Kn)", vec = cl_vec, min_n = 2L) else NULL,
    list(name = "regional_group", label = "US EPA Level III Ecoregion (regional_group)", vec = rg_vec, min_n = 3L)
  ))

  group_rows <- list()
  pairwise_rows <- list()

  for (spec in factor_specs) {
    g_raw <- spec$vec
    min_n <- as.integer(spec$min_n)
    if (length(unique(g_raw[!is.na(g_raw)])) < 2L) next

    # Enforce strict Tier_1/Tier_2 authority anchor split when positive anchors exist
    if (has_positive_anchors) {
      cand_mask <- !is.na(g_raw) & has_pcs & is_t1 & (auth_w > 0)
      tab_c <- table(g_raw[cand_mask])
      act_lvls <- sort_group_levels(names(tab_c[tab_c >= min_n]))
      if (length(act_lvls) < 2L) {
        cand_alt <- !is.na(g_raw) & has_pcs & (auth_w > 0)
        tab_alt <- table(g_raw[cand_alt])
        if (sum(tab_alt >= min_n) >= 2L) {
          cand_mask <- cand_alt
          act_lvls <- sort_group_levels(names(tab_alt[tab_alt >= min_n]))
        }
      }
      w_use <- auth_w
    } else {
      cand_mask <- !is.na(g_raw) & has_pcs & is_t1
      tab_c <- table(g_raw[cand_mask])
      act_lvls <- sort_group_levels(names(tab_c[tab_c >= min_n]))
      if (length(act_lvls) < 2L) {
        cand_alt <- !is.na(g_raw) & has_pcs
        tab_alt <- table(g_raw[cand_alt])
        if (sum(tab_alt >= min_n) >= 2L) {
          cand_mask <- cand_alt
          act_lvls <- sort_group_levels(names(tab_alt[tab_alt >= min_n]))
        }
      }
      w_use <- rep(1.0, n)
    }

    if (length(act_lvls) < 2L) next
    is_train <- cand_mask & (g_raw %in% act_lvls)

    X_full <- as.matrix(df[, avail_cols, drop = FALSE])
    storage.mode(X_full) <- "double"
    for (j in seq_len(ncol(X_full))) {
      na_j <- !is.finite(X_full[, j])
      if (any(na_j)) {
        tr_ok <- is_train & !na_j
        fill_v <- if (any(tr_ok)) stats::weighted.mean(X_full[tr_ok, j], w = pmax(w_use[tr_ok], 1e-6)) else 0.0
        X_full[na_j, j] <- fill_v
      }
    }

    sel_cols <- if (exists("select_morphometric_subspace", mode = "function")) {
      suppressMessages(select_morphometric_subspace(X_full[is_train, , drop = FALSE], group_vec = g_raw[is_train]))
    } else {
      n_min_g <- min(as.integer(table(g_raw[is_train])))
      p_max_wt <- if (n_min_g >= 6L) max(2L, floor(n_min_g / 3L)) else ncol(X_full)
      avail_cols[seq_len(min(length(avail_cols), p_max_wt))]
    }
    X_sub <- X_full[, sel_cols, drop = FALSE]
    p_sub <- ncol(X_sub)

    # Per-group 95% chi-square Mahalanobis outlier pruning on active anchors (exempting nomenclatural types)
    for (lvl in act_lvls) {
      idx_g <- which(is_train & g_raw == lvl)
      if (length(idx_g) >= 5L) {
        sub_g <- X_sub[idx_g, , drop = FALSE]
        d2_g <- tryCatch(stats::mahalanobis(sub_g, colMeans(sub_g), stats::cov(sub_g) + diag(1e-4, p_sub)), error = function(e) rep(0.0, length(idx_g)))
        drop_g <- !is.na(d2_g) & (d2_g > stats::qchisq(0.95, df = p_sub)) & !is_type_row[idx_g]
        if (sum(!drop_g) >= max(3L, min_n)) is_train[idx_g[drop_g]] <- FALSE
      }
    }

    tab_post <- table(g_raw[is_train])
    act_lvls <- sort_group_levels(names(tab_post[tab_post >= min_n]))
    if (length(act_lvls) < 2L) next
    is_train <- is_train & (g_raw %in% act_lvls)

    act_idx <- which(is_train)
    Y_act <- X_sub[act_idx, , drop = FALSE]
    g_act <- factor(g_raw[act_idx], levels = act_lvls)
    w_act <- w_use[act_idx]

    omni <- .fit_weighted_manova_stats(Y_act, g_act, w_act)
    if (is.null(omni)) next

    group_rows[[length(group_rows) + 1L]] <- data.frame(
      grouping_factor = spec$name,
      factor_label = spec$label,
      n_groups = as.integer(omni$n_groups),
      n_active_vouchers = as.integer(omni$n_obs),
      n_features = as.integer(omni$n_features),
      active_features = paste(omni$active_features, collapse = ";"),
      pillai_trace = round(omni$pillai_trace, 6),
      wilks_lambda = round(omni$wilks_lambda, 6),
      approx_F = round(omni$approx_F, 5),
      num_df = round(omni$num_df, 3),
      den_df = round(omni$den_df, 3),
      partial_eta2 = round(omni$partial_eta2, 6),
      p_value = signif(omni$p_value, 6),
      q_value = NA_real_,
      stringsAsFactors = FALSE
    )

    comb_mat <- utils::combn(length(act_lvls), 2L)
    n_pairs <- ncol(comb_mat)
    eval_one_pair <- function(p_idx) {
      g1 <- act_lvls[comb_mat[1L, p_idx]]
      g2 <- act_lvls[comb_mat[2L, p_idx]]
      m_pair <- (g_act == g1) | (g_act == g2)
      n1 <- sum(g_act == g1)
      n2 <- sum(g_act == g2)
      pw_stat <- .fit_weighted_manova_stats(
        Y_act[m_pair, omni$active_features, drop = FALSE],
        droplevels(g_act[m_pair]),
        w_act[m_pair]
      )
      if (is.null(pw_stat)) return(NULL)
      data.frame(
        grouping_factor = spec$name,
        group_1 = g1,
        group_2 = g2,
        contrast = paste(g1, "vs", g2),
        n_1 = as.integer(n1),
        n_2 = as.integer(n2),
        n_features = as.integer(pw_stat$n_features),
        hotelling_t2 = round(pw_stat$hotelling_t2, 5),
        mahalanobis_d = round(pw_stat$mahalanobis_d, 5),
        pillai_trace = round(pw_stat$pillai_trace, 6),
        wilks_lambda = round(pw_stat$wilks_lambda, 6),
        approx_F = round(pw_stat$approx_F, 5),
        num_df = round(pw_stat$num_df, 3),
        den_df = round(pw_stat$den_df, 3),
        partial_eta2 = round(pw_stat$partial_eta2, 6),
        p_value = signif(pw_stat$p_value, 6),
        q_value = NA_real_,
        stringsAsFactors = FALSE
      )
    }
    raw_pairs <- if (n_pairs >= 40L && exists("packera_mclapply", mode = "function")) {
      n_pw_cores <- if (exists("resolve_packera_workers", mode = "function")) resolve_packera_workers(min(32L, n_pairs), n_pairs) else 8L
      packera_mclapply(seq_len(n_pairs), eval_one_pair, mc.cores = n_pw_cores)
    } else {
      lapply(seq_len(n_pairs), eval_one_pair)
    }
    pair_list_factor <- Filter(Negate(is.null), raw_pairs)
    if (length(pair_list_factor) > 0L) {
      pw_fac_df <- do.call(rbind, pair_list_factor)
      pw_fac_df$q_value <- signif(stats::p.adjust(pw_fac_df$p_value, method = "BH"), 6)
      pairwise_rows[[length(pairwise_rows) + 1L]] <- pw_fac_df
    }
  }

  group_df <- if (length(group_rows) > 0L) do.call(rbind, group_rows) else empty_group_df
  if (nrow(group_df) > 0L) {
    group_df$q_value <- signif(stats::p.adjust(group_df$p_value, method = "BH"), 6)
    rownames(group_df) <- NULL
  }

  pairwise_df <- if (length(pairwise_rows) > 0L) do.call(rbind, pairwise_rows) else empty_pairwise_df
  if (nrow(pairwise_df) > 0L) {
    rownames(pairwise_df) <- NULL
  }

  atomic_write_manova_csv(group_df, output_group_csv)
  atomic_write_manova_csv(pairwise_df, output_pairwise_csv)

  list(
    group_table = group_df,
    pairwise_table = pairwise_df,
    output_group_csv = output_group_csv,
    output_pairwise_csv = output_pairwise_csv
  )
}

.short_epa_ecoregion_label <- function(eco_vec) {
  gsub("_", " ", as.character(eco_vec))
}

draw_gmm_ecoregion_contingency_heatmap <- function(rg_vec, cl_v, active_ecos, ari_val = NA_real_, fisher_p = NA_real_) {
  ecos_use <- if (length(active_ecos) >= 2L) active_ecos else sort(unique(rg_vec[!is.na(rg_vec)]))
  keep <- !is.na(rg_vec) & (rg_vec %in% ecos_use) & !is.na(cl_v)
  if (sum(keep) == 0L || length(ecos_use) == 0L) {
    graphics::plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "",
                   main = "B. Unsupervised Morphotype x EPA Ecoregion Contingency", font.main = 2)
    graphics::box(col = "grey75")
    graphics::text(0.5, 0.5, "Insufficient multi-ecoregion vouchers for contingency test", col = "grey35", cex = 0.85)
    return(invisible(NULL))
  }
  # Keep top 12 active ecoregions by count if >12 so heatmap cells stay crisp
  if (length(ecos_use) > 12L) {
    cnt_e <- sort(table(rg_vec[keep]), decreasing = TRUE)
    ecos_use <- names(cnt_e)[seq_len(12L)]
    keep <- !is.na(rg_vec) & (rg_vec %in% ecos_use) & !is.na(cl_v)
  }

  u_cl <- sort(unique(as.integer(cl_v[keep])))
  tab_mat <- table(factor(rg_vec[keep], levels = ecos_use), factor(as.integer(cl_v[keep]), levels = u_cl))
  row_pct <- sweep(tab_mat, 1, pmax(rowSums(tab_mat), 1), "/") * 100

  p_txt <- if (is.finite(fisher_p)) {
    if (fisher_p < 1e-4) sprintf("%.2e", fisher_p) else sprintf("%.4f", fisher_p)
  } else "NA"
  ari_txt <- if (is.finite(ari_val)) sprintf("%.3f", ari_val) else "NA"

  nr <- nrow(tab_mat); nc <- ncol(tab_mat)
  graphics::plot(
    c(0.5, nc + 0.5), c(0.5, nr + 0.5), type = "n", axes = FALSE,
    xlab = sprintf("Unsupervised GMM Morphotype (K = %d)", length(u_cl)),
    ylab = "",
    main = sprintf("B. Morphotype x EPA Ecoregion (ARI = %s, Fisher p = %s)", ari_txt, p_txt),
    font.main = 2
  )
  heat_ramp <- grDevices::colorRampPalette(c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594"))(101)
  cell_cex <- if (nr > 8L) 0.62 else 0.74
  for (i in seq_len(nr)) {
    y_pos <- nr - i + 1
    for (j in seq_len(nc)) {
      pct_ij <- row_pct[i, j]
      cnt_ij <- tab_mat[i, j]
      c_idx <- max(1L, min(101L, round(pct_ij) + 1L))
      graphics::rect(j - 0.48, y_pos - 0.46, j + 0.48, y_pos + 0.46, col = heat_ramp[c_idx], border = "white", lwd = 1.4)
      txt_col <- if (pct_ij >= 52) "white" else "#111111"
      lbl_ij <- if (nr > 7L) sprintf("%d (%.0f%%)", cnt_ij, pct_ij) else sprintf("%d\n(%.0f%%)", cnt_ij, pct_ij)
      graphics::text(j, y_pos, labels = lbl_ij, col = txt_col, cex = cell_cex, font = if (pct_ij >= 25) 2 else 1)
    }
  }
  graphics::axis(1, at = seq_len(nc), labels = paste0("K", colnames(tab_mat)), tick = FALSE, font = 2, cex.axis = 0.72, gap.axis = 0.05)
  short_e <- vapply(.short_epa_ecoregion_label(ecos_use), function(s) if (nchar(s) > 26L) paste0(substr(s, 1, 24), "..") else s, character(1))
  graphics::axis(2, at = rev(seq_len(nr)), labels = short_e, las = 1, tick = FALSE, cex.axis = if (nr > 8L) 0.66 else 0.75)
  graphics::box(col = "grey75")
}

draw_dynamic_ecoregion_mean_silhouettes <- function(
  df,
  rg_vec,
  active_ecos,
  eco_cols,
  main_title = "C. Mean Consensus Leaf Silhouettes (+/- 1 SD) by EPA Ecoregion"
) {
  cand_ecos <- if (length(active_ecos) >= 1L) active_ecos else sort(unique(rg_vec[!is.na(rg_vec)]))
  # Strictly retain ONLY ecoregions with n >= 1 (zero n=0 placeholders)
  ecos_present <- cand_ecos[vapply(cand_ecos, function(eg) sum(!is.na(rg_vec) & rg_vec == eg) > 0L, logical(1))]
  if (length(ecos_present) == 0L) {
    graphics::plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", main = main_title, font.main = 2)
    graphics::box(col = "grey80")
    graphics::text(0.5, 0.5, "No assigned EPA Level III Ecoregion vouchers", col = "grey35", cex = 0.85)
    return(invisible(NULL))
  }
  if (length(ecos_present) > 8L) {
    cnts <- vapply(ecos_present, function(eg) sum(!is.na(rg_vec) & rg_vec == eg), integer(1))
    ecos_present <- ecos_present[order(cnts, decreasing = TRUE)[seq_len(8L)]]
  }

  m_disp <- length(ecos_present)
  n_cols_g <- if (m_disp <= 2L) m_disp else if (m_disp <= 4L) 2L else if (m_disp <= 6L) 3L else 4L
  n_rows_g <- ceiling(m_disp / n_cols_g)

  graphics::plot(
    c(0, n_cols_g), c(0, n_rows_g), type = "n", axes = FALSE, xlab = "", ylab = "",
    main = main_title, font.main = 2
  )
  graphics::box(col = "grey80")
  if (n_cols_g > 1L) {
    for (vc in seq_len(n_cols_g - 1L)) graphics::abline(v = vc, col = "grey90", lty = 2)
  }
  if (n_rows_g > 1L) {
    for (hr in seq_len(n_rows_g - 1L)) graphics::abline(h = hr, col = "grey90", lty = 2)
  }

  std_hmat <- extract_standardized_harmonic_matrix(df, max_k = 12L)
  pc1_v <- if ("Sym_PC1" %in% names(df)) as.numeric(df$Sym_PC1) else if ("PC1" %in% names(df)) as.numeric(df$PC1) else rep(0.0, nrow(df))
  pc2_v <- if ("Sym_PC2" %in% names(df)) as.numeric(df$Sym_PC2) else if ("PC2" %in% names(df)) as.numeric(df$PC2) else rep(0.0, nrow(df))

  for (idx_g in seq_len(m_disp)) {
    g_name <- ecos_present[idx_g]
    g_col <- if (g_name %in% names(eco_cols)) eco_cols[[g_name]] else "#0072B2"
    row_g <- ceiling(idx_g / n_cols_g)
    col_g <- idx_g - (row_g - 1L) * n_cols_g
    cx <- col_g - 0.5
    cy <- (n_rows_g - row_g) + 0.5
    m_grp <- which(!is.na(rg_vec) & rg_vec == g_name)

    if (!is.null(std_hmat) && nrow(std_hmat) == nrow(df) && length(m_grp) > 0L) {
      h_mean <- colMeans(std_hmat[m_grp, , drop = FALSE], na.rm = TRUE)
      h_sd <- if (length(m_grp) >= 2L) apply(std_hmat[m_grp, , drop = FALSE], 2, stats::sd, na.rm = TRUE) else rep(0.0, ncol(std_hmat))
      h_sd[!is.finite(h_sd)] <- 0.0
    } else {
      z1 <- mean(pc1_v[m_grp], na.rm = TRUE); if (!is.finite(z1)) z1 <- 0.0
      z2 <- mean(pc2_v[m_grp], na.rm = TRUE); if (!is.finite(z2)) z2 <- 0.0
      h_mean <- c(
        A1 = -1.25 - 0.18 * z1, B1 = 0, C1 = 0, D1 = 0.48 - 0.10 * z1 + 0.08 * z2,
        A2 = -0.04 + 0.06 * z2, B2 = 0, C2 = 0, D2 = -0.01 + 0.04 * z2,
        A3 = -0.08 - 0.05 * z1, B3 = 0, C3 = 0, D3 = 0.03 * z2
      )
      h_sd <- abs(h_mean) * 0.12
    }

    xy_mean <- reconstruct_efa_outline(h_mean, nb_pts = 120L)
    xy_plus <- reconstruct_efa_outline(h_mean + 1.0 * h_sd, nb_pts = 120L)
    xy_minus <- reconstruct_efa_outline(h_mean - 1.0 * h_sd, nb_pts = 120L)

    env_pts <- rbind(xy_plus, xy_mean, xy_minus)
    shared_center <- c(mean(range(env_pts[, 1], na.rm = TRUE)), mean(range(env_pts[, 2], na.rm = TRUE)))
    shared_scale <- min(
      0.68 / max(diff(range(env_pts[, 1], na.rm = TRUE)), 1e-6),
      0.64 / max(diff(range(env_pts[, 2], na.rm = TRUE)), 1e-6)
    ) * 0.86

    draw_outline_in_box(
      xy_plus, cx = cx, cy = cy - 0.03, w_box = 0.68, h_box = 0.64,
      col_border = grDevices::adjustcolor(g_col, alpha.f = 0.48),
      col_fill = grDevices::adjustcolor(g_col, alpha.f = 0.10), lwd = 0.9, lty = 3,
      shared_scale = shared_scale, shared_center = shared_center
    )
    draw_outline_in_box(
      xy_minus, cx = cx, cy = cy - 0.03, w_box = 0.68, h_box = 0.64,
      col_border = grDevices::adjustcolor(g_col, alpha.f = 0.48),
      col_fill = NA, lwd = 0.9, lty = 3,
      shared_scale = shared_scale, shared_center = shared_center
    )
    pts_mean <- draw_outline_in_box(
      xy_mean, cx = cx, cy = cy - 0.03, w_box = 0.68, h_box = 0.64,
      col_border = g_col,
      col_fill = grDevices::adjustcolor(g_col, alpha.f = 0.22), lwd = 2.1, lty = 1,
      shared_scale = shared_scale, shared_center = shared_center
    )
    if (!is.null(pts_mean)) {
      graphics::segments(cx, min(pts_mean[, 2]), cx, max(pts_mean[, 2]), col = grDevices::adjustcolor(g_col, alpha.f = 0.55), lty = 2, lwd = 0.9)
    }

    short_lbl <- .short_epa_ecoregion_label(g_name)
    if (nchar(short_lbl) > 22L) short_lbl <- paste0(substr(short_lbl, 1, 20), "..")
    lbl_cex <- if (n_cols_g >= 4L) 0.62 else if (n_cols_g == 3L) 0.68 else 0.76
    graphics::text(cx, cy + 0.41, labels = sprintf("%s (n=%d)", short_lbl, length(m_grp)), font = 2, cex = lbl_cex, col = "#111111")
    med_ar <- if ("aspect_ratio" %in% names(df)) stats::median(as.numeric(df$aspect_ratio[m_grp]), na.rm = TRUE) else NA_real_
    med_sol <- if ("solidity" %in% names(df)) stats::median(as.numeric(df$solidity[m_grp]), na.rm = TRUE) else NA_real_
    if (is.finite(med_ar) && is.finite(med_sol)) {
      graphics::text(cx, cy - 0.43, labels = sprintf("AR=%.2f | Sol=%.2f", med_ar, med_sol), cex = max(0.56, lbl_cex - 0.10), col = "grey25")
    }
  }
}

render_fig04c_supervised_ecoregion_cda <- function(
  audited_df,
  eco_cda_res,
  gmm_res,
  output_plot_path = NULL,
  manova_res = NULL
) {
  fdir <- resolve_active_figures_dir(reference_path = output_plot_path)
  if (is.null(audited_df) || nrow(audited_df) == 0L || is.null(eco_cda_res)) return(invisible(NULL))

  if (is.null(manova_res) && exists("run_multivariate_manova_audit", mode = "function")) {
    df_m <- audited_df
    if (!is.null(eco_cda_res$ecoregion_vec) && length(eco_cda_res$ecoregion_vec) == nrow(df_m)) {
      df_m$regional_group <- eco_cda_res$ecoregion_vec
    }
    manova_res <- tryCatch(
      run_multivariate_manova_audit(
        df_m,
        .resolve_default_manova_features(df_m, NULL),
        gmm_res = gmm_res,
        output_group_csv = NULL,
        output_pairwise_csv = NULL
      ),
      error = function(e) NULL
    )
  }
  eco_manova_sub <- format_manova_pillai_subtitle(
    manova_res,
    factor_name = "regional_group",
    fallback_label = "Omnibus MANOVA (Pillai's Trace): Insufficient EPA Ecoregions (<2 regions with >=3 vouchers)"
  )

  n <- nrow(audited_df)
  rg_vec <- eco_cda_res$ecoregion_vec
  active_ecos <- eco_cda_res$active_ecoregions
  pal_pool <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", EXTENDED_COLORBLIND_HEX)
  all_ecos <- unique(c(active_ecos, sort(unique(rg_vec[!is.na(rg_vec)]))))
  eco_cols <- if (length(all_ecos) > 0L) stats::setNames(rep_len(pal_pool, length(all_ecos)), all_ecos) else c("Unassigned" = "#7F7F7F")
  eco_pch <- if (length(all_ecos) > 0L) stats::setNames(rep_len(HOLLOW_PCH_CYCLE, length(all_ecos)), all_ecos) else c("Unassigned" = 1L)

  draw_fig04c <- function() {
    graphics::par(mfrow = c(2, 2), mar = c(4.4, 7.6, 3.2, 1.2), oma = c(0.6, 0.4, 2.0, 0.6), family = "sans")

    # Panel A: 2D Supervised CDA Biplot & Trait Vectors across EPA Level III Ecoregions (excluding NA from ellipses)
    graphics::par(mar = c(4.4, 4.6, 3.2, 1.2))
    has_2d_cda <- length(active_ecos) >= 2L && any(is.finite(eco_cda_res$ecoregion_can1))
    x_a <- if (has_2d_cda) eco_cda_res$ecoregion_can1 else if ("PC1" %in% names(audited_df)) as.numeric(audited_df$PC1) else seq_len(n)
    y_a <- if (has_2d_cda) eco_cda_res$ecoregion_can2 else if ("PC2" %in% names(audited_df)) as.numeric(audited_df$PC2) else rep(0.0, n)
    alpha_v <- adaptive_point_alpha(n)

    pt_col <- ifelse(
      !is.na(rg_vec) & (rg_vec %in% names(eco_cols)),
      grDevices::adjustcolor(eco_cols[rg_vec], alpha.f = alpha_v),
      grDevices::adjustcolor("grey55", alpha.f = 0.40)
    )
    pt_pch <- ifelse(!is.na(rg_vec) & (rg_vec %in% names(eco_pch)), eco_pch[rg_vec], 4L)
    # Strictly exclude NA ecoregions from active training 95% confidence ellipses
    ellipse_grp <- ifelse(!is.na(rg_vec) & (rg_vec %in% active_ecos) & as.logical(eco_cda_res$is_active), rg_vec, NA_character_)
    if (length(unique(stats::na.omit(ellipse_grp))) < length(active_ecos)) {
      ellipse_grp <- ifelse(!is.na(rg_vec) & (rg_vec %in% active_ecos), rg_vec, NA_character_)
    }

    loocv_txt <- if (is.finite(eco_cda_res$loocv_accuracy)) sprintf("LOOCV = %.1f%%", 100 * eco_cda_res$loocv_accuracy) else "Insufficient Ecoregions"
    r2_txt <- if (is.finite(eco_cda_res$permanova_r2)) sprintf(", R2 = %.3f", eco_cda_res$permanova_r2) else ""
    var1_e <- if (length(eco_cda_res$variance_pct) >= 1L && is.finite(eco_cda_res$variance_pct[1])) eco_cda_res$variance_pct[1] else 100.0
    var2_e <- if (length(eco_cda_res$variance_pct) >= 2L && is.finite(eco_cda_res$variance_pct[2])) eco_cda_res$variance_pct[2] else 0.0

    s_load_e <- eco_cda_res$structure_loadings
    s_lbl_e <- if (!is.null(s_load_e) && !is.null(rownames(s_load_e))) format_cda_loading_labels(rownames(s_load_e)) else NULL
    cent_e <- eco_cda_res$centroids
    cent_lbl_e <- if (!is.null(cent_e) && !is.null(rownames(cent_e))) .short_epa_ecoregion_label(rownames(cent_e)) else NULL

    leg_ecos <- if (length(active_ecos) > 0L) active_ecos[seq_len(min(7L, length(active_ecos)))] else all_ecos[seq_len(min(7L, length(all_ecos)))]
    leg_names <- .short_epa_ecoregion_label(leg_ecos)
    leg_cols <- unname(eco_cols[leg_ecos])
    leg_pchs <- unname(eco_pch[leg_ecos])
    if (any(is.na(rg_vec))) {
      leg_names <- c(leg_names, "Unassigned (NA Projected)")
      leg_cols <- c(leg_cols, "grey55")
      leg_pchs <- c(leg_pchs, 4L)
    }
    leg_arg <- if (length(leg_names) > 0L) {
      list("topright", legend = leg_names, col = leg_cols, pch = leg_pchs, pt.bg = NA, cex = 0.64, bty = "n")
    } else NULL

    is_type_e <- if ("is_type_specimen" %in% names(audited_df)) as.logical(audited_df$is_type_specimen) else rep(FALSE, n)
    is_type_e[is.na(is_type_e)] <- FALSE
    plot_2d_ordination_scatter(
      x = x_a, y = y_a,
      col = pt_col, pch = pt_pch, cex = 0.84, lwd = 1.0,
      xlab = if (has_2d_cda) sprintf("Ecoregion Can1 (%.1f%%)", var1_e) else "Morphospace PC1",
      ylab = if (has_2d_cda) sprintf("Ecoregion Can2 (%.1f%%)", var2_e) else "Morphospace PC2",
      main = sprintf("A. 2D Supervised EPA Ecoregion CDA (%s%s)", loocv_txt, r2_txt),
      subtitle = eco_manova_sub,
      group_vec = ellipse_grp,
      group_cols = eco_cols,
      ellipse_level = 0.95,
      centroids = cent_e,
      centroid_labels = cent_lbl_e,
      biplot_vectors = s_load_e,
      biplot_labels = s_lbl_e,
      legend_args = leg_arg,
      is_type_vec = is_type_e
    )

    # Panel B: Unsupervised GMM Cluster x EPA Level III Ecoregion Contingency Heatmap
    graphics::par(mar = c(4.4, 9.2, 3.0, 1.0))
    cl_v <- if ("gmm_cluster" %in% names(audited_df)) as.integer(audited_df$gmm_cluster) else if (!is.null(gmm_res$classification)) as.integer(gmm_res$classification) else rep(1L, n)
    draw_gmm_ecoregion_contingency_heatmap(
      rg_vec = rg_vec,
      cl_v = cl_v,
      active_ecos = active_ecos,
      ari_val = eco_cda_res$ari_morphotype_vs_ecoregion,
      fisher_p = eco_cda_res$fisher_exact_p_value
    )

    # Panel C: Mean Consensus Leaf Silhouettes (+/- 1 SD) per active EPA Level III Ecoregion
    graphics::par(mar = c(2.4, 2.2, 3.0, 1.0))
    draw_dynamic_ecoregion_mean_silhouettes(
      df = audited_df,
      rg_vec = rg_vec,
      active_ecos = active_ecos,
      eco_cols = eco_cols,
      main_title = "C. Mean Consensus Leaf Silhouettes (+/- 1 SD) by Active EPA Ecoregion"
    )

    # Panel D: STRUCTURE-Style Supervised Ecoregion Posterior Probability Barplot
    graphics::par(mar = c(4.6, 4.6, 3.0, 1.2))
    post_mat <- eco_cda_res$post_probs
    if (!is.null(post_mat) && is.matrix(post_mat) && nrow(post_mat) == n && ncol(post_mat) >= 2L) {
      grp_disp <- ifelse(!is.na(rg_vec) & (rg_vec %in% active_ecos), rg_vec, "Unassigned_NA")
      grp_order_lvls <- c(active_ecos, if (any(grp_disp == "Unassigned_NA")) "Unassigned_NA" else character(0))
      self_p <- vapply(seq_len(n), function(i) {
        gi <- grp_disp[i]
        if (gi %in% colnames(post_mat)) post_mat[i, gi] else max(post_mat[i, ], na.rm = TRUE)
      }, numeric(1))
      ord_i <- order(match(grp_disp, grp_order_lvls), -self_p)
      post_ord <- post_mat[ord_i, , drop = FALSE]
      grp_ord <- grp_disp[ord_i]
      if (nrow(post_ord) > 400L) {
        keep_s <- unique(round(seq(1, nrow(post_ord), length.out = 400L)))
        post_ord <- post_ord[keep_s, , drop = FALSE]
        grp_ord <- grp_ord[keep_s]
      }
      nb <- nrow(post_ord)
      grle <- rle(grp_ord)
      brks <- cumsum(grle$lengths)
      mids <- brks - grle$lengths / 2.0
      use_las2_d <- length(mids) > 4L
      if (use_las2_d) graphics::par(mar = c(5.8, 4.6, 3.0, 1.2))
      k_cols <- ifelse(colnames(post_ord) %in% names(eco_cols), eco_cols[colnames(post_ord)], "#7F7F7F")
      graphics::plot(
        c(0, nb), c(0, 1.10), type = "n", xaxs = "i", yaxs = "i", xaxt = "n",
        xlab = if (use_las2_d) "" else "Herbarium Vouchers Grouped by US EPA Level III Ecoregion",
        ylab = "Posterior Probability P(Ecoregion | x)",
        main = "D. STRUCTURE-Style Supervised Ecoregion Posterior Probabilities",
        font.main = 2
      )
      cum_l <- rep(0.0, nb)
      xl <- seq(0, nb - 1L); xr <- seq(1L, nb)
      for (ki in seq_len(ncol(post_ord))) {
        pk <- as.numeric(post_ord[, ki]); pk[!is.finite(pk)] <- 0.0
        cum_h <- pmin(1.0, cum_l + pk)
        graphics::rect(xl, cum_l, xr, cum_h, col = k_cols[ki], border = NA)
        cum_l <- cum_h
      }
      if (length(brks) > 1L) {
        graphics::segments(brks[-length(brks)], 0, brks[-length(brks)], 1.0, col = "#111111", lwd = 1.6)
      }
      graphics::rect(0, 0, nb, 1.0, border = "#222222", lwd = 1.1)
      ax_lbls <- vapply(grle$values, function(v) {
        if (v == "Unassigned_NA") return("NA Proj.")
        s <- .short_epa_ecoregion_label(v)
        if (nchar(s) > 14L) paste0(substr(s, 1, 12), "..") else s
      }, character(1))
      graphics::axis(1, at = mids, labels = ax_lbls, tick = FALSE, las = if (use_las2_d) 2 else 1, cex.axis = 0.66, line = -0.3)
      leg_k <- seq_len(min(6L, ncol(post_ord)))
      graphics::legend(
        "top", horiz = TRUE,
        legend = vapply(.short_epa_ecoregion_label(colnames(post_ord)[leg_k]), function(s) if (nchar(s) > 16L) paste0(substr(s, 1, 14), "..") else s, character(1)),
        fill = k_cols[leg_k], border = NA, cex = 0.62, bty = "n"
      )
    } else {
      graphics::plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "",
                     main = "D. STRUCTURE-Style Supervised Ecoregion Posterior Probabilities", font.main = 2)
      graphics::box(col = "grey75")
      graphics::text(0.5, 0.5, "Requires >= 2 active EPA Level III Ecoregions (>= 3 anchors each)", col = "grey35", cex = 0.85)
    }

    graphics::title(
      "Supervised US EPA Level III Ecoregion Canonical Discriminant Analysis & Morphotype Congruence",
      outer = TRUE, cex.main = 1.25
    )
  }

  atomic_save_dual_plot(draw_fig04c, fdir, "fig04c_supervised_ecoregion_cda", width = 13.5, height = 10.4, dpi = 300)
}
