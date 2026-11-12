#!/usr/bin/env Rscript
# ==============================================================================
# Script: scripts/morphometrics/pub_ecoregion_manova_utils.R
# Project: Packera Genus-Wide Species Delimitation & Morphometrics Pipeline
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
  abbr_map <- c(
    stats::setNames(tolower(datasets::state.name), tolower(datasets::state.abb)),
    "dc" = "district of columbia", "d.c." = "district of columbia",
    "on" = "ontario", "ont" = "ontario",
    "qc" = "quebec", "pq" = "quebec", "que" = "quebec",
    "nb" = "new brunswick", "ns" = "nova scotia",
    "pe" = "prince edward island", "pei" = "prince edward island",
    "nl" = "newfoundland and labrador", "nf" = "newfoundland and labrador", "nfld" = "newfoundland and labrador",
    "mb" = "manitoba", "man" = "manitoba",
    "sk" = "saskatchewan", "sask" = "saskatchewan",
    "ab" = "alberta", "alta" = "alberta",
    "bc" = "british columbia", "yt" = "yukon", "nt" = "northwest territories", "nu" = "nunavut"
  )
  vapply(as.character(st_vec), function(s) {
    if (is.na(s)) return(NA_character_)
    s_ascii <- iconv(s, from = "", to = "ASCII//TRANSLIT")
    if (is.na(s_ascii)) s_ascii <- s
    s_clean <- gsub("\\.$", "", gsub("\\s+", " ", tolower(trimws(s_ascii))))
    if (!nzchar(s_clean) || s_clean %in% c("na", "nan", "none", "unknown")) return(NA_character_)
    if (s_clean %in% names(abbr_map)) unname(abbr_map[[s_clean]]) else s_clean
  }, character(1), USE.NAMES = FALSE)
}

.HOMONYMOUS_INDEP_CITIES_R <- list(
  "virginia" = c("fairfax", "franklin", "richmond", "roanoke"),
  "maryland" = c("baltimore"),
  "missouri" = c("st louis")
)

.norm_us_county_r <- function(co_vec, st_norm_vec = NULL) {
  n <- length(co_vec)
  if (is.null(st_norm_vec)) st_norm_vec <- rep(NA_character_, n)
  vapply(seq_len(n), function(i) {
    raw <- as.character(co_vec[i])
    if (is.na(raw)) return(NA_character_)
    c_ascii <- iconv(raw, from = "", to = "ASCII//TRANSLIT")
    if (is.na(c_ascii)) c_ascii <- raw
    c_str <- tolower(trimws(c_ascii))
    if (!nzchar(c_str) || c_str %in% c("na", "nan", "none", "unknown", "unplaced", "indeterminate")) {
      return(NA_character_)
    }
    if (grepl(" or |/|\\?", c_str)) return(NA_character_)
    c_str <- gsub("\\[|\\]|\\(|\\)", "", c_str)
    was_explicit_city <- grepl("^city\\s+of\\s+|\\b(independent\\s+city|city)$", c_str) && !grepl("\\bcounty\\b", c_str)
    c_str <- gsub("^(regional\\s+municipality\\s+of|united\\s+counties\\s+of|county\\s+of|city\\s+of|mrc\\s+de|mrc\\s+d'|municipalite\\s+regionale\\s+de\\s+comte\\s+de)\\s+", "", c_str)
    c_str <- gsub("\\s+(regional\\s+municipality|united\\s+counties|county\\s+municipality|census\\s+division|census\\s+area|city\\s+and\\s+borough|independent\\s+city|independent|municipality|county|parish|par\\.?|borough|district|mrc|co\\.?)$", "", c_str)
    c_str <- gsub("\\s+city$", "", c_str)
    c_str <- gsub("^(st\\.|saint)\\s+", "st ", c_str)
    c_str <- gsub("^ste\\.\\s*", "ste ", c_str)
    c_str <- trimws(gsub("[^a-z0-9]+", " ", c_str))
    if (!nzchar(c_str)) return(NA_character_)
    st_i <- as.character(st_norm_vec[i])
    if (!is.na(st_i)) {
      if (st_i == "virginia" && c_str == "nansemond") c_str <- "suffolk"
      else if (st_i == "virginia" && c_str == "princess anne") c_str <- "virginia beach"
      else if (st_i == "virginia" && c_str == "bedford") was_explicit_city <- FALSE
      else if (st_i == "georgia" && c_str == "columbus") c_str <- "muscogee"
      else if (st_i == "south dakota" && c_str == "shannon") c_str <- "oglala lakota"
      if (isTRUE(was_explicit_city) && (st_i %in% names(.HOMONYMOUS_INDEP_CITIES_R)) && (c_str %in% .HOMONYMOUS_INDEP_CITIES_R[[st_i]])) {
        c_str <- paste0(c_str, " city")
      }
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
      "Ozark_Highlands", "Interior_Plateau", "Southwestern_Appalachians",
      "Mixed_Wood_Plains", "Northern_Lakes_and_Forests", "Atlantic_Highlands",
      "Lake_Manitoba_and_Lake_Winnipeg_Plain", "Aspen_Parkland_Northern_Glaciated_Plains"
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
    mk_poly(matrix(c(-83.5,41.8, -71.0,45.0, -71.0,47.2, -83.5,44.8), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-95.2,46.0, -76.0,45.2, -76.0,50.5, -95.2,50.5), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-71.0,45.0, -64.0,45.0, -64.0,49.2, -71.0,49.2), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-99.5,49.0, -96.0,49.0, -96.0,53.0, -99.5,53.0), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-104.05,49.0, -99.5,49.0, -99.5,53.5, -104.05,53.5), ncol = 2, byrow = TRUE)),
    crs = 4326
  )
  eco_sf <- sf::st_sf(eco_df, geometry = eco_geoms)
  cnty_df <- data.frame(
    county_key = c("north carolina|orange", "north carolina|new hanover", "north carolina|wake", "texas|nacogdoches", "ontario|essex", "quebec|gatineau"),
    state_name = c("north carolina", "north carolina", "north carolina", "texas", "ontario", "quebec"),
    county_norm = c("orange", "new hanover", "wake", "nacogdoches", "essex", "gatineau"),
    centroid_lon = c(-79.10, -77.90, -78.40, -94.65, -82.80, -75.95),
    centroid_lat = c(36.05, 34.20, 35.78, 31.62, 42.15, 46.10),
    n_ecoregions = rep(NA_integer_, 6L),
    top_ecoregion = rep(NA_character_, 6L),
    top_area_share = rep(NA_real_, 6L),
    stringsAsFactors = FALSE
  )
  cnty_geoms <- sf::st_sfc(
    mk_poly(matrix(c(-79.25,35.86, -78.95,35.86, -78.95,36.24, -79.25,36.24), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-78.02,34.05, -77.78,34.05, -77.78,34.35, -78.02,34.35), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-79.00,35.50, -77.80,35.50, -77.80,36.05, -79.00,36.05), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-94.90,31.40, -94.40,31.40, -94.40,31.85, -94.90,31.85), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-83.10,41.95, -82.50,41.95, -82.50,42.32, -83.10,42.32), ncol = 2, byrow = TRUE)),
    mk_poly(matrix(c(-76.30,45.75, -75.65,45.75, -75.65,46.45, -76.30,46.45), ncol = 2, byrow = TRUE)),
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
  na_eco_p <- file.path(default_dir, "na_cec_eco_l3.geojson")
  ca_cnty_p <- file.path(default_dir, "ca_census_divisions.geojson")
  eco_json_p <- if (!is.null(ecoregion_path)) as.character(ecoregion_path) else if (file.exists(na_eco_p)) na_eco_p else file.path(default_dir, "us_eco_l3.geojson")
  cnty_json_p <- if (!is.null(county_path)) as.character(county_path) else file.path(default_dir, "us_counties.geojson")

  cache_id <- paste(
    if (is.null(ecoregion_path)) "DEFAULT_ECO" else eco_json_p,
    if (is.null(county_path)) "DEFAULT_CNTY" else cnty_json_p,
    sep = "::"
  )
  if (exists(cache_id, envir = .epa_l3_r_cache, inherits = FALSE)) {
    return(get(cache_id, envir = .epa_l3_r_cache, inherits = FALSE))
  }

  .merge_ca_cec_r <- function(eco_sf_in, cnty_sf_in) {
    if (!custom_paths && file.exists(na_eco_p) && !identical(eco_json_p, na_eco_p)) {
      na_sf <- tryCatch(sf::st_read(na_eco_p, quiet = TRUE), error = function(e) NULL)
      if (!is.null(na_sf) && nrow(na_sf) > 0L) {
        na_lbl <- if ("ecoregion_label" %in% names(na_sf)) na_sf$ecoregion_label else if ("NA_L3NAME" %in% names(na_sf)) na_sf$NA_L3NAME else if ("US_L3NAME" %in% names(na_sf)) na_sf$US_L3NAME else na_sf[[1]]
        na_sf$ecoregion_label <- standardize_epa_l3_label_r(na_lbl)
        na_sf <- sf::st_transform(na_sf[, "ecoregion_label", drop = FALSE], sf::st_crs(eco_sf_in))
        gcol_e <- attr(eco_sf_in, "sf_column")
        if (!is.null(gcol_e) && !identical(attr(na_sf, "sf_column"), gcol_e)) sf::st_geometry(na_sf) <- gcol_e
        eco_sf_in <- rbind(eco_sf_in[, "ecoregion_label", drop = FALSE], na_sf[, "ecoregion_label", drop = FALSE])
      }
    }
    if (!is.null(eco_sf_in) && !is.na(sf::st_crs(eco_sf_in)) && !isTRUE(sf::st_is_longlat(eco_sf_in))) {
      eco_sf_in <- sf::st_transform(eco_sf_in, 4326)
    }
    if (!is.null(cnty_sf_in) && !is.na(sf::st_crs(cnty_sf_in)) && !isTRUE(sf::st_is_longlat(cnty_sf_in))) {
      cnty_sf_in <- sf::st_transform(cnty_sf_in, 4326)
    }
    if (!custom_paths && file.exists(ca_cnty_p) && !is.null(cnty_sf_in)) {
      ca_sf <- tryCatch(sf::st_read(ca_cnty_p, quiet = TRUE), error = function(e) NULL)
      if (!is.null(ca_sf) && nrow(ca_sf) > 0L) {
        ca_sf <- sf::st_transform(ca_sf, sf::st_crs(cnty_sf_in))
        gcol_c <- attr(cnty_sf_in, "sf_column")
        if (!is.null(gcol_c) && !identical(attr(ca_sf, "sf_column"), gcol_c)) sf::st_geometry(ca_sf) <- gcol_c
        ca_st_c <- intersect(c("state_name", "PRNAME", "STATE_NAME", "state"), names(ca_sf))
        ca_co_c <- intersect(c("county_norm", "CDNAME", "NAME", "county"), names(ca_sf))
        if (length(ca_st_c) > 0L && length(ca_co_c) > 0L) {
          ca_st_v <- .norm_us_state_r(ca_sf[[ca_st_c[1L]]])
          ca_co_v <- .norm_us_county_r(ca_sf[[ca_co_c[1L]]], ca_st_v)
          ca_sf$state_name <- ca_st_v
          ca_sf$county_norm <- ca_co_v
          ca_sf$county_key <- paste0(ca_st_v, "|", ca_co_v)
        }
        non_geom_c <- setdiff(names(cnty_sf_in), gcol_c)
        for (miss_c in setdiff(non_geom_c, names(ca_sf))) ca_sf[[miss_c]] <- NA
        cnty_sf_in <- rbind(cnty_sf_in[, c(non_geom_c, gcol_c), drop = FALSE], ca_sf[, c(non_geom_c, gcol_c), drop = FALSE])
      }
    }
    list(eco_sf = eco_sf_in, cnty_sf = cnty_sf_in)
  }

  if (!custom_paths && file.exists(gpkg_p)) {
    eco_sf <- sf::st_read(gpkg_p, layer = "us_eco_l3", quiet = TRUE)
    cnty_sf <- tryCatch(sf::st_read(gpkg_p, layer = "us_counties", quiet = TRUE), error = function(e) NULL)
    eco_sf$ecoregion_label <- standardize_epa_l3_label_r(
      if ("ecoregion_label" %in% names(eco_sf)) eco_sf$ecoregion_label else if ("NA_L3NAME" %in% names(eco_sf)) eco_sf$NA_L3NAME else eco_sf$US_L3NAME
    )
    m_ca <- .merge_ca_cec_r(eco_sf, cnty_sf)
    bundle <- list(eco_sf = m_ca$eco_sf, cnty_sf = m_ca$cnty_sf, use_live_county_intersection = FALSE)
    assign(cache_id, bundle, envir = .epa_l3_r_cache)
    return(bundle)
  }

  if (file.exists(eco_json_p)) {
    eco_sf <- sf::st_read(eco_json_p, quiet = TRUE)
    raw_lbl <- if ("ecoregion_label" %in% names(eco_sf)) eco_sf$ecoregion_label else if ("NA_L3NAME" %in% names(eco_sf)) eco_sf$NA_L3NAME else if ("US_L3NAME" %in% names(eco_sf)) eco_sf$US_L3NAME else eco_sf[[1]]
    eco_sf$ecoregion_label <- standardize_epa_l3_label_r(raw_lbl)
    cnty_sf <- if (file.exists(cnty_json_p)) tryCatch(sf::st_read(cnty_json_p, quiet = TRUE), error = function(e) NULL) else NULL
    m_ca <- .merge_ca_cec_r(eco_sf, cnty_sf)
    bundle <- list(eco_sf = m_ca$eco_sf, cnty_sf = m_ca$cnty_sf, use_live_county_intersection = custom_paths || identical(eco_json_p, na_eco_p))
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
    c_st_col <- intersect(c("state_name", "STATE_NAME", "state", "PRNAME"), names(cnty_sf))
    c_st <- .norm_us_state_r(if (length(c_st_col) > 0L) cnty_sf[[c_st_col[1L]]] else if ("county_key" %in% names(cnty_sf)) sub("\\|.*$", "", as.character(cnty_sf$county_key)) else cnty_sf[[1]])
    c_co_col <- intersect(c("county_norm", "NAME", "CDNAME", "county"), names(cnty_sf))
    c_co <- .norm_us_county_r(if (length(c_co_col) > 0L) cnty_sf[[c_co_col[1L]]] else if ("county_key" %in% names(cnty_sf)) sub("^.*\\|", "", as.character(cnty_sf$county_key)) else cnty_sf[[min(2L, ncol(cnty_sf))]], c_st)
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
    lat_col <- intersect(c("decimalLatitude"), names(df))
    lon_col <- intersect(c("decimalLongitude"), names(df))
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
      select_morphometric_subspace(X_full[is_train, , drop = FALSE], group_vec = rg[is_train], allow_supervised_f = FALSE)
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

    if (g_n >= 2L) {
      sds_a <- apply(X_a, 2, stats::sd, na.rm = TRUE)
      keep_sd <- is.finite(sds_a) & (sds_a > 1e-7)
      if (!any(keep_sd)) keep_sd[1] <- TRUE
      X_a_fit <- X_a[, keep_sd, drop = FALSE]
      X_proj_fit <- X_use[, keep_sd, drop = FALSE]
      p_ef <- ncol(X_a_fit)

      sp_w_eco <- if (exists("compute_spatial_balance_weights", mode = "function")) compute_spatial_balance_weights(df, active_mask = is_train, recompute = TRUE) else rep(1.0, n)
      w_eco_a <- auth_w[active_idx] * sp_w_eco[active_idx]
      if (!any(is.finite(w_eco_a) & w_eco_a > 0.0)) w_eco_a <- rep(1.0, length(active_idx))

      if (exists("compute_weighted_shrinkage_cda_loocv", mode = "function")) {
        loocv_acc <- tryCatch(
          compute_weighted_shrinkage_cda_loocv(
            X_a_fit, grp_a, w_a = w_eco_a,
            coords_df = df[active_idx, , drop = FALSE],
            X_cand = X_full[active_idx, , drop = FALSE],
            deme_radius_km = 10.0
          ),
          error = function(e) NA_real_
        )
      }

      eco_lvls <- levels(grp_a)
      gm_e <- colSums(X_a_fit * w_eco_a) / sum(w_eco_a)
      B_e <- matrix(0.0, p_ef, p_ef); W_e <- matrix(0.0, p_ef, p_ef)
      Xc_e <- matrix(0.0, nrow(X_a_fit), p_ef)
      mu_e_mat <- matrix(0.0, g_n, p_ef, dimnames = list(eco_lvls, colnames(X_a_fit)))
      for (ki in seq_len(g_n)) {
        km <- which(grp_a == eco_lvls[ki])
        wk <- w_eco_a[km]; if (sum(wk) <= 0.0) wk <- rep(1.0, length(km))
        mk <- colSums(X_a_fit[km, , drop = FALSE] * wk) / sum(wk)
        mu_e_mat[ki, ] <- mk
        dm <- mk - gm_e; B_e <- B_e + sum(wk) * tcrossprod(dm)
        dx <- sweep(X_a_fit[km, , drop = FALSE], 2L, mk, "-"); Xc_e[km, ] <- dx
        W_e <- W_e + crossprod(dx, dx * wk)
      }
      denom_e <- max(sum(w_eco_a) - sum(vapply(eco_lvls, function(eg) { wk <- w_eco_a[grp_a == eg]; if (sum(wk) > 0) sum(wk^2) / sum(wk) else 1.0 }, numeric(1))), 1.0)
      S_w_e <- W_e / denom_e
      S_reg_e <- if (exists("compute_ledoit_wolf_cov", mode = "function")) {
        compute_ledoit_wolf_cov(S_w_e, X_centered = Xc_e, n_obs = nrow(X_a_fit), weights = w_eco_a)
      } else {
        S_w_e + diag(1e-4, p_ef)
      }
      k_ef <- max(1L, min(p_ef, g_n - 1L))
      eig_e <- solve_shrinkage_cda_eigen(S_reg_e, B_e, max_axes = k_ef)
      Z_full_e <- sweep(X_proj_fit, 2L, gm_e, "-") %*% eig_e$vectors
      cent_full_e <- sweep(mu_e_mat, 2L, gm_e, "-") %*% eig_e$vectors

      lp_mat <- matrix(0.0, nrow = n, ncol = g_n, dimnames = list(NULL, eco_lvls))
      for (ki in seq_len(g_n)) {
        dz <- sweep(Z_full_e, 2L, cent_full_e[ki, ], "-")
        lp_mat[, ki] <- log(priors_u[ki]) - 0.5 * rowSums(dz^2)
      }
      post_mat <- exp(lp_mat - apply(lp_mat, 1L, max, na.rm = TRUE))
      post_mat <- post_mat / pmax(rowSums(post_mat, na.rm = TRUE), 1e-12)
      pred_eco <- eco_lvls[max.col(post_mat, ties.method = "first")]
      post_max <- round(apply(post_mat, 1L, max, na.rm = TRUE), 4)

      Z_raw <- Z_full_e
      if (ncol(Z_raw) == 1L) {
        v1 <- as.numeric(Z_raw[, 1])
        src2 <- if ("PC2" %in% colnames(X_full)) X_full[, "PC2"] else if (ncol(X_full) >= 2L) X_full[, 2] else X_full[, 1]
        v2 <- as.numeric(src2) - mean(src2, na.rm = TRUE)
        b1 <- v1 - mean(v1, na.rm = TRUE)
        v2 <- v2 - (sum(v2 * b1, na.rm = TRUE) / max(sum(b1^2, na.rm = TRUE), 1e-12)) * b1
        Z_raw <- cbind(Can1 = v1, Can2 = v2)
      }
      can1_out <- round(as.numeric(Z_raw[, 1]), 5)
      can2_out <- round(as.numeric(Z_raw[, 2]), 5)

      ev_e <- pmax(eig_e$values, 0.0)
      if (length(ev_e) >= 2L) {
        var_pct <- round(100.0 * ev_e[1:2] / max(sum(ev_e), 1e-12), 2)
      } else {
        var_pct <- c(100.0, 0.0)
      }
      Z_act_2d <- cbind(can1_out[active_idx], can2_out[active_idx])
      struct_load <- tryCatch(stats::cor(X_a_fit, Z_act_2d), error = function(e) NULL)
      centroids_2d <- do.call(rbind, lapply(eco_lvls, function(eg) {
        m_e <- active_idx[grp_a == eg]
        w_m <- w_eco_a[grp_a == eg]
        c(stats::weighted.mean(can1_out[m_e], w = w_m, na.rm = TRUE), stats::weighted.mean(can2_out[m_e], w = w_m, na.rm = TRUE))
      }))
      rownames(centroids_2d) <- eco_lvls
      colnames(centroids_2d) <- c("Can1", "Can2")
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
    tmp_sum <- paste0(output_summary_csv, ".", Sys.getpid(), ".tmp")
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
      "scape_to_basal_leaf_ratio_std",
      "aspect_ratio_std",
      "solidity_std"
    )
  }
  avail <- intersect(unique(cand), names(df))
  if (length(avail) == 0L) avail <- intersect(paste0("PC", seq_len(5)), names(df))
  avail
}

format_manova_pillai_subtitle <- function(
  manova_res,
  factor_name = NULL,
  alt_factor = NULL,
  fallback_label = "",
  group_pattern = NULL
) {
  fb_out <- if (is.null(fallback_label)) "" else as.character(fallback_label)[1L]
  gt <- if (!is.null(manova_res) && is.list(manova_res) && !is.null(manova_res$group_table)) {
    manova_res$group_table
  } else if (is.data.frame(manova_res)) {
    manova_res
  } else {
    NULL
  }
  if (!is.null(gt) && nrow(gt) > 0L && "grouping_factor" %in% names(gt)) {
    idx <- integer(0)
    if (!is.null(factor_name) && nzchar(as.character(factor_name)[1L])) {
      idx <- which(gt$grouping_factor == as.character(factor_name)[1L])
    }
    if (length(idx) == 0L && !is.null(alt_factor) && nzchar(as.character(alt_factor)[1L])) {
      idx <- which(gt$grouping_factor == as.character(alt_factor)[1L])
    }
    if (length(idx) == 0L && !is.null(group_pattern) && nzchar(as.character(group_pattern)[1L])) {
      pat_g <- as.character(group_pattern)[1L]
      m_gf <- grepl(pat_g, as.character(gt$grouping_factor), ignore.case = TRUE)
      m_fl <- if ("factor_label" %in% names(gt)) grepl(pat_g, as.character(gt$factor_label), ignore.case = TRUE) else rep(FALSE, nrow(gt))
      idx <- which(m_gf | m_fl)
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
      } else if (is.finite(v_val) && !is.finite(p_val)) {
        return(sprintf("Descriptive Multivariate Separation: Pillai's Trace V = %.3f (post-clustering p suppressed)", v_val))
      }
    }
  }
  fb_out
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
  include_gmm_clusters = FALSE
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
    lat_col <- intersect(c("decimalLatitude"), names(df))
    lon_col <- intersect(c("decimalLongitude"), names(df))
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
    if (has_positive_anchors && exists("compute_spatial_balance_weights", mode = "function")) w_use <- auth_w * compute_spatial_balance_weights(df, active_mask = is_train, recompute = TRUE)

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
      suppressMessages(select_morphometric_subspace(X_full[is_train, , drop = FALSE], group_vec = g_raw[is_train], allow_supervised_f = FALSE))
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

    is_post_sel_gmm <- identical(spec$name, "gmm_cluster")
    group_rows[[length(group_rows) + 1L]] <- data.frame(
      grouping_factor = spec$name,
      factor_label = spec$label,
      n_groups = as.integer(omni$n_groups),
      n_active_vouchers = as.integer(omni$n_obs),
      n_features = as.integer(omni$n_features),
      active_features = paste(omni$active_features, collapse = ";"),
      pillai_trace = round(omni$pillai_trace, 6),
      wilks_lambda = round(omni$wilks_lambda, 6),
      approx_F = if (is_post_sel_gmm) NA_real_ else round(omni$approx_F, 5),
      num_df = if (is_post_sel_gmm) NA_real_ else round(omni$num_df, 3),
      den_df = if (is_post_sel_gmm) NA_real_ else round(omni$den_df, 3),
      partial_eta2 = round(omni$partial_eta2, 6),
      p_value = if (is_post_sel_gmm) NA_real_ else signif(omni$p_value, 6),
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
        approx_F = if (is_post_sel_gmm) NA_real_ else round(pw_stat$approx_F, 5),
        num_df = if (is_post_sel_gmm) NA_real_ else round(pw_stat$num_df, 3),
        den_df = if (is_post_sel_gmm) NA_real_ else round(pw_stat$den_df, 3),
        partial_eta2 = round(pw_stat$partial_eta2, 6),
        p_value = if (is_post_sel_gmm) NA_real_ else signif(pw_stat$p_value, 6),
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
      fin_pw <- is.finite(pw_fac_df$p_value)
      if (any(fin_pw)) pw_fac_df$q_value[fin_pw] <- signif(stats::p.adjust(pw_fac_df$p_value[fin_pw], method = "BH"), 6)
      pairwise_rows[[length(pairwise_rows) + 1L]] <- pw_fac_df
    }
  }

  group_df <- if (length(group_rows) > 0L) do.call(rbind, group_rows) else empty_group_df
  if (nrow(group_df) > 0L) {
    fin_grp <- is.finite(group_df$p_value)
    if (any(fin_grp)) group_df$q_value[fin_grp] <- signif(stats::p.adjust(group_df$p_value[fin_grp], method = "BH"), 6)
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

render_fig04c_supervised_ecoregion_cda <- function(
  audit_df,
  eco_cda_res = NULL,
  gmm_res = NULL,
  figures_dir = NULL,
  output_plot_path = NULL,
  manova_res = NULL,
  batch_csv = NULL
) {
  fdir <- if (!is.null(figures_dir) && nzchar(figures_dir)) {
    as.character(figures_dir)
  } else if (exists("resolve_active_figures_dir", mode = "function")) {
    resolve_active_figures_dir(NULL, output_plot_path)
  } else {
    "outputs/figures"
  }
  if (!exists("render_step04_ecoregion_batch_png", mode = "function") && file.exists("scripts/figures/step04_figures.R")) {
    if (file.exists("scripts/figures/theme_packera_pub.R")) source("scripts/figures/theme_packera_pub.R")
    source("scripts/figures/step04_figures.R")
  }
  render_step04_ecoregion_batch_png(
    audit_df = audit_df,
    eco_cda_res = eco_cda_res,
    gmm_res = gmm_res,
    figures_dir = fdir,
    manova_res = manova_res,
    batch_csv = batch_csv
  )
}
