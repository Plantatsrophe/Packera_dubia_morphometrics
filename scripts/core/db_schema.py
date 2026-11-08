"""
===============================================================================
Module: db_schema.py
Project: Packera Species Delimitation & Morphometrics Pipeline
Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)

Description:
    Two-Tier Relational Schema DDL (Tier A + Tier B + R*Tree + Master View)
    and idempotent dynamic column migration utilities for SQLite.
===============================================================================
"""

from __future__ import annotations

import sqlite3
from typing import Dict, List, Set, Tuple

HARMONIC_COLS: Tuple[str, ...] = tuple(
    f"{coef}{idx}" for idx in range(1, 25) for coef in ("A", "B", "C", "D")
)
PC_COLS: Tuple[str, ...] = tuple(f"PC{idx}" for idx in range(1, 16))
SYM_PC_COLS: Tuple[str, ...] = tuple(f"Sym_PC{idx}" for idx in range(1, 16))

WC_BIO_COLUMNS: Tuple[str, ...] = (
    "bio1_temp_mean",
    "bio2_mean_diurnal_range",
    "bio3_isothermality",
    "bio4_temp_seasonality",
    "bio5_max_temp_warmest_month",
    "bio6_min_temp_coldest_month",
    "bio7_temp_annual_range",
    "bio8_mean_temp_wettest_qtr",
    "bio9_mean_temp_driest_qtr",
    "bio10_mean_temp_warmest_qtr",
    "bio11_mean_temp_coldest_qtr",
    "bio12_precip_annual",
    "bio13_precip_wettest_month",
    "bio14_precip_driest_month",
    "bio15_precip_seasonality",
    "bio16_precip_wettest_qtr",
    "bio17_precip_driest_qtr",
    "bio18_precip_warmest_qtr",
    "bio19_precip_coldest_qtr",
)

ENVIRONMENTAL_CACHE_COLUMNS: Tuple[str, ...] = (
    "lat_lon_key",
    "catalogNumber",
    "latitude",
    "longitude",
    "musym",
    "muname",
    "taxsuborder",
    "taxgreatgroup",
    "is_rock_outcrop",
    "outcrop_class",
    "ssurgo_rock_outcrop_pct",
    "ssurgo_bedrock_depth_cm",
    "ssurgo_sand_pct",
    "ssurgo_clay_pct",
    "ssurgo_ph",
    "ssurgo_om_pct",
    "ssurgo_cec",
    "ssurgo_aws_0_50cm",
    "polaris_ph_30m",
    "polaris_sand_30m",
    "polaris_clay_30m",
    "polaris_om_30m",
    "polaris_bd_30m",
    "polaris_ksat_30m",
    "solus_resdept_cm",
    "solus_fragvol_pct",
    "topo_tpi_30m",
    "topo_roughness_30m",
    "topo_slope_deg",
    "topo_twi_30m",
    "topo_elev_m",
    "fine_soil_ph",
    "fine_soil_sand",
    "fine_soil_clay",
    "fine_bedrock_depth_cm",
    "fine_rock_outcrop_prob",
    "bdricm",
    "resdept_r",
    "cfvo",
    "awc",
    "fine_soil_cec",
    "soil_silt_pct",
    "alr_sand_clay",
    "alr_silt_clay",
    "clr_sand_silt",
    "clr_clay_geom",
    "sgmc_lithology_class",
    "is_azonal_lithology",
    "azonal_contact_dist_m",
    "azonal_proximity_index",
    "buffer_radius_m",
    "high_coord_uncertainty_flag",
    "is_training_eligible_coord",
    "ph_q10_buffer",
    "ph_q90_buffer",
    "min_bedrock_depth_buffer_cm",
    "max_cfvo_buffer_pct",
    "azonal_lithology_presence_buffer",
    *WC_BIO_COLUMNS,
    "bio18_precip_warmest_qt",
    "bio19_precip_coldest_qt",
    "soil_ph",
    "soil_cec",
    "soil_sand",
    "soil_bulk_density",
    "soil_ph_soilgrids",
    "soil_cec_soilgrids",
    "soil_sand_soilgrids",
    "soil_bulk_density_soilgrids",
)

ENVIRONMENTAL_TEXT_COLS: Set[str] = {
    "lat_lon_key",
    "catalogNumber",
    "musym",
    "muname",
    "taxsuborder",
    "taxgreatgroup",
    "outcrop_class",
    "sgmc_lithology_class",
}
ENVIRONMENTAL_BOOL_COLS: Set[str] = {
    "is_rock_outcrop",
    "is_azonal_lithology",
    "high_coord_uncertainty_flag",
    "is_training_eligible_coord",
    "azonal_lithology_presence_buffer",
}


def _build_master_view_ddl() -> str:
    wc_excluded_from_env: Set[str] = {
        "lat_lon_key",
        "catalogNumber",
        "latitude",
        "longitude",
        "bio18_precip_warmest_qt",
        "bio19_precip_coldest_qt",
        *WC_BIO_COLUMNS,
    }
    env_view_selects = ",\n        ".join(
        f"ec.{col} AS {col}"
        for col in ENVIRONMENTAL_CACHE_COLUMNS
        if col not in wc_excluded_from_env
    )

    wc_select_parts: List[str] = []
    for col in WC_BIO_COLUMNS:
        if col == "bio18_precip_warmest_qtr":
            wc_select_parts.append(
                "COALESCE(srf.bio18_precip_warmest_qtr, srf.bio18_precip_warmest_qt, "
                "ec.bio18_precip_warmest_qtr, ec.bio18_precip_warmest_qt) AS bio18_precip_warmest_qtr"
            )
        elif col == "bio19_precip_coldest_qtr":
            wc_select_parts.append(
                "COALESCE(srf.bio19_precip_coldest_qtr, srf.bio19_precip_coldest_qt, "
                "ec.bio19_precip_coldest_qtr, ec.bio19_precip_coldest_qt) AS bio19_precip_coldest_qtr"
            )
        else:
            wc_select_parts.append(f"COALESCE(srf.{col}, ec.{col}) AS {col}")
    wc_select_parts.append(
        "COALESCE(srf.bio18_precip_warmest_qt, ec.bio18_precip_warmest_qt, "
        "srf.bio18_precip_warmest_qtr, ec.bio18_precip_warmest_qtr) AS bio18_precip_warmest_qt"
    )
    wc_select_parts.append(
        "COALESCE(srf.bio19_precip_coldest_qt, ec.bio19_precip_coldest_qt, "
        "srf.bio19_precip_coldest_qtr, ec.bio19_precip_coldest_qtr) AS bio19_precip_coldest_qt"
    )
    wc_view_selects = ",\n        ".join(wc_select_parts)

    harm_view_selects = ",\n        ".join(f"om.{h} AS {h}" for h in HARMONIC_COLS)
    pc_view_selects = ",\n        ".join(f"om.{p} AS {p}" for p in PC_COLS)
    sym_pc_view_selects = ",\n        ".join(f"om.{sp} AS {sp}" for sp in SYM_PC_COLS)

    return f"""
    DROP VIEW IF EXISTS v_integrative_triage_master;
    CREATE VIEW v_integrative_triage_master AS
    SELECT
        rac.run_id AS run_id,
        rac.catalogNumber AS catalogNumber,
        rac.plant_individual_id AS plant_individual_id,
        v.institutionCode AS institutionCode,
        v.scientificName AS scientificName,
        v.species_raw AS species_raw,
        v.recordedBy AS recordedBy,
        v.recordNumber AS recordNumber,
        v.identifiedBy AS identifiedBy,
        v.determiner_raw AS determiner_raw,
        COALESCE(rac.determiner_tier, om.determiner_tier, v.determiner_tier) AS determiner_tier,
        COALESCE(rac.type_status, v.type_status, om.type_status) AS type_status,
        COALESCE(
            rac.is_type_specimen,
            om.is_type_specimen,
            CASE
                WHEN COALESCE(rac.type_status, v.type_status, '') != ''
                 AND UPPER(TRIM(COALESCE(rac.type_status, v.type_status, ''))) NOT IN ('NONE', 'NOT A TYPE', 'NOTATYPE', 'NON-TYPE', 'NONTYPE', 'UNSPECIFIED', 'NULL', 'NAN', 'NA', 'N/A', 'FALSE', '0', 'UNKNOWN')
                THEN 1
                ELSE 0
            END
        ) AS is_type_specimen,
        v.county AS county,
        v.stateProvince AS stateProvince,
        v.locality AS locality,
        COALESCE(v.decimalLatitude, v.latitude) AS decimalLatitude,
        COALESCE(v.latitude, v.decimalLatitude) AS latitude,
        COALESCE(v.decimalLongitude, v.longitude) AS decimalLongitude,
        COALESCE(v.longitude, v.decimalLongitude) AS longitude,
        v.coordinateUncertainty AS coordinateUncertainty,
        v.year AS year,
        v.month AS month,
        v.day AS day,
        v.eventDate AS eventDate,
        v.regional_group AS regional_group,
        v.image_path AS image_path,
        v.exsiccatae_key AS exsiccatae_key,
        v.is_primary_duplicate AS is_primary_duplicate,
        v.duplicate_count AS duplicate_count,
        COALESCE(rac.species_standardized, om.species_standardized, v.species_standardized) AS species_standardized,
        COALESCE(rac.is_promoted_bronze, om.is_promoted_anchor, 0) AS is_promoted_anchor,
        COALESCE(rac.authority_weight, om.authority_weight) AS authority_weight,
        COALESCE(rac.spatial_balance_weight, 1.0) AS spatial_balance_weight,
        COALESCE(om.pixels_per_mm, v.pixels_per_mm) AS pixels_per_mm,
        COALESCE(om.is_calibrated_mm, v.is_calibrated_mm) AS is_calibrated_mm,
        v.qc_status AS qc_status,
        -- Plant individual & linear petiole/scape traits
        pi.rosette_centroid_x AS rosette_centroid_x,
        pi.rosette_centroid_y AS rosette_centroid_y,
        pi.crown_x AS crown_x,
        pi.crown_y AS crown_y,
        pi.plant_height_px AS plant_height_px,
        pi.plant_height_mm AS plant_height_mm,
        pi.plant_height_cm AS plant_height_cm,
        pi.is_bent_stem AS is_bent_stem,
        pi.petiole_length_px AS petiole_length_px,
        COALESCE(om.petiole_length_mm, pi.petiole_length_mm) AS petiole_length_mm,
        pi.blade_length_px AS blade_length_px,
        COALESCE(om.blade_length_mm, pi.blade_length_mm) AS blade_length_mm,
        COALESCE(om.petiole_to_blade_ratio, pi.petiole_to_blade_ratio) AS petiole_to_blade_ratio,
        pi.scape_to_basal_leaf_ratio AS scape_to_basal_leaf_ratio,
        pi.rosette_patch_path AS rosette_patch_path,
        pi.rosette_mask_path AS rosette_mask_path,
        -- Reproductive & phenological traits (with sheet-level pid=0 fallback)
        COALESCE(vrp.capitula_count, vrp0.capitula_count) AS capitula_count,
        COALESCE(vrp.involucre_height_px, vrp0.involucre_height_px) AS involucre_height_px,
        COALESCE(vrp.involucre_width_px, vrp0.involucre_width_px) AS involucre_width_px,
        COALESCE(vrp.involucre_height_mm, vrp0.involucre_height_mm) AS involucre_height_mm,
        COALESCE(vrp.involucre_width_mm, vrp0.involucre_width_mm) AS involucre_width_mm,
        COALESCE(om.AR_inv, vrp.AR_inv, vrp.capitulum_aspect_ratio, vrp0.AR_inv, vrp0.capitulum_aspect_ratio) AS AR_inv,
        COALESCE(vrp.capitulum_aspect_ratio, vrp.AR_inv, vrp0.capitulum_aspect_ratio, vrp0.AR_inv) AS capitulum_aspect_ratio,
        COALESCE(vrp.total_capitula, vrp0.total_capitula, vrp.capitula_count, vrp0.capitula_count) AS total_capitula,
        COALESCE(vrp.n_anthesis, vrp0.n_anthesis) AS n_anthesis,
        COALESCE(vrp.n_bud, vrp0.n_bud) AS n_bud,
        COALESCE(vrp.n_fruiting, vrp.n_fruit, vrp0.n_fruiting, vrp0.n_fruit) AS n_fruiting,
        COALESCE(vrp.n_fruit, vrp.n_fruiting, vrp0.n_fruit, vrp0.n_fruiting) AS n_fruit,
        COALESCE(vrp.n_senescent, vrp0.n_senescent) AS n_senescent,
        COALESCE(vrp.voucher_phenological_state, vrp0.voucher_phenological_state) AS voucher_phenological_state,
        COALESCE(vrp.is_verified_anthesis, vrp0.is_verified_anthesis) AS is_verified_anthesis,
        COALESCE(vrp.reproductive_stage_confidence, vrp0.reproductive_stage_confidence) AS reproductive_stage_confidence,
        COALESCE(vrp.dominant_pappus_ratio, vrp0.dominant_pappus_ratio) AS dominant_pappus_ratio,
        -- Step 03 & Step 04 OTU Morphometrics & CDA
        om.assigned_tier AS assigned_tier,
        om.solidity AS solidity,
        om.dissection_index AS dissection_index,
        om.q25_solidity AS q25_solidity,
        om.dissected_leaf_fraction AS dissected_leaf_fraction,
        om.basal_leaf_count AS basal_leaf_count,
        om.ucs_score AS ucs_score,
        om.aspect_ratio AS aspect_ratio,
        om.area_px AS area_px,
        om.is_dissected AS is_dissected,
        om.centroid_size AS centroid_size,
        om.log_centsize AS log_centsize,
        {harm_view_selects},
        om.sym_var AS sym_var,
        om.asym_var AS asym_var,
        om.pct_asymmetry AS pct_asymmetry,
        om.harmonic_roughness AS harmonic_roughness,
        om.is_contour_outlier AS is_contour_outlier,
        om.leaf_quality_weight AS leaf_quality_weight,
        om.leaf_count AS leaf_count,
        om.n_tier1_leaves AS n_tier1_leaves,
        om.n_tier2_leaves AS n_tier2_leaves,
        om.morph_tier_source AS morph_tier_source,
        om.foliar_variance AS foliar_variance,
        om.heteroblasty_index AS heteroblasty_index,
        om.efd_harmonic_opt AS efd_harmonic_opt,
        om.efd_k95_pcs AS efd_k95_pcs,
        om.selected_morph_features AS selected_morph_features,
        {pc_view_selects},
        {sym_pc_view_selects},
        om.gmm_cluster AS gmm_cluster,
        om.gmm_uncertainty AS gmm_uncertainty,
        om.cda_predicted_taxon AS cda_predicted_taxon,
        om.cda_posterior_prob AS cda_posterior_prob,
        om.cda_second_taxon AS cda_second_taxon,
        om.cda_second_posterior_prob AS cda_second_posterior_prob,
        om.cda_posterior_margin AS cda_posterior_margin,
        om.loocv_accuracy AS loocv_accuracy,
        om.can1 AS can1,
        om.can2 AS can2,
        om.can3 AS can3,
        om.ecoregion_cda_predicted AS ecoregion_cda_predicted,
        om.ecoregion_cda_posterior_prob AS ecoregion_cda_posterior_prob,
        om.morphometric_outlier_flag AS morphometric_outlier_flag,
        COALESCE(om.misidentification_flag, om.morphometric_outlier_flag) AS misidentification_flag,
        om.triage_priority AS triage_priority,
        om.discordance_reason AS morph_discordance_reason,
        om.has_empirical_type_anchor AS has_empirical_type_anchor,
        om.type_morph_mahal_d2 AS type_morph_mahal_d2,
        om.nearest_type_morph_taxon AS nearest_type_morph_taxon,
        om.type_morph_concordance AS type_morph_concordance,
        om.is_typical_morphotype AS is_typical_morphotype,
        om.min_congener_mahal_d2 AS min_congener_mahal_d2,
        om.min_congener_taxon AS min_congener_taxon,
        om.open_set_novelty_p AS open_set_novelty_p,
        om.is_confident_label_noise AS is_confident_label_noise,
        -- Step 06 Spatial RF & Hopkins Phenology
        COALESCE(srf.lat_lon_key, ec.lat_lon_key) AS lat_lon_key,
        srf.rf_predicted_taxon AS rf_predicted_taxon,
        srf.rf_predicted_prob AS rf_predicted_prob,
        srf.rf_given_taxon_prob AS rf_given_taxon_prob,
        srf.edaphic_best_fit_taxon AS edaphic_best_fit_taxon,
        srf.spatial_block_cv_acc AS spatial_block_cv_acc,
        srf.spatial_autocorr_range_km AS spatial_autocorr_range_km,
        srf.spatial_block_id AS spatial_block_id,
        COALESCE(srf.doy, v.doy) AS doy,
        v.pheno_sin AS pheno_sin,
        v.pheno_cos AS pheno_cos,
        srf.expected_doy AS expected_doy,
        srf.hopkins_expected AS hopkins_expected,
        srf.delta_doy AS delta_doy,
        srf.pheno_timing_category AS pheno_timing_category,
        srf.anthesis_filter_status AS anthesis_filter_status,
        srf.edaphic_z_score AS edaphic_z_score,
        srf.climatic_z_score AS climatic_z_score,
        srf.type_edaphic_z_score AS type_edaphic_z_score,
        srf.nearest_type_niche_taxon AS nearest_type_niche_taxon,
        srf.matches_type_lithology AS matches_type_lithology,
        srf.primary_env_predictors AS primary_env_predictors,
        srf.primary_clim_predictors AS primary_clim_predictors,
        srf.primary_edaph_predictors AS primary_edaph_predictors,
        srf.substrate_specialist_flag AS substrate_specialist_flag,
        srf.piedmont_outcrop_specialist_flag AS piedmont_outcrop_specialist_flag,
        srf.multimodal_concordance AS multimodal_concordance,
        srf.multimodal_conflict_flag AS multimodal_conflict_flag,
        srf.triage_category AS triage_category,
        srf.conflict_discordance_notes AS conflict_discordance_notes,
        -- Environmental Cache (Fine-grain SSURGO / POLARIS / SOLUS / 3DEP / SGMC & 19 WorldClim Rasters)
        {env_view_selects},
        {wc_view_selects},
        -- Step 07 Curation Triage
        ct.w_morph_kappa AS w_morph_kappa,
        ct.w_vision_kappa AS w_vision_kappa,
        ct.w_edaphic_kappa AS w_edaphic_kappa,
        ct.w_pheno_kappa AS w_pheno_kappa,
        ct.n_active_modalities AS n_active_modalities,
        ct.n_available_modalities AS n_available_modalities,
        ct.probabilistic_consensus_conflict AS probabilistic_consensus_conflict,
        ct.composite_conflict_score AS composite_conflict_score,
        ct.probabilistic_consensus_score AS probabilistic_consensus_score,
        COALESCE(ct.hopkins_adjusted_doy, srf.hopkins_adjusted_doy) AS hopkins_adjusted_doy,
        ct.delta_doy_resid AS delta_doy_resid,
        ct.taxonomic_status_call AS taxonomic_status_call,
        ct.recommended_determination AS recommended_determination,
        ct.taxonomic_treatment AS taxonomic_treatment,
        ct.taxonomic_treatment_broad AS taxonomic_treatment_broad,
        ct.taxonomic_treatment_subtype AS taxonomic_treatment_subtype,
        ct.synthesis_triage_priority AS synthesis_triage_priority,
        ct.synthesis_triage_action AS synthesis_triage_action,
        ct.synthesis_rationale AS synthesis_rationale,
        ct.cryptic_cluster_cohesion_pass AS cryptic_cluster_cohesion_pass,
        ct.is_completed_voucher AS is_completed_voucher
    FROM run_active_cohort rac
    LEFT JOIN vouchers v
        ON v.catalogNumber = rac.catalogNumber
    LEFT JOIN plant_individuals pi
        ON pi.catalogNumber = rac.catalogNumber
       AND pi.plant_individual_id = rac.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp
        ON vrp.catalogNumber = rac.catalogNumber
       AND vrp.plant_individual_id = rac.plant_individual_id
    LEFT JOIN voucher_reproductive_phenology vrp0
        ON vrp0.catalogNumber = rac.catalogNumber
       AND vrp0.plant_individual_id = 0
    LEFT JOIN otu_morphometrics om
        ON om.run_id = rac.run_id
       AND om.catalogNumber = rac.catalogNumber
       AND om.plant_individual_id = rac.plant_individual_id
    LEFT JOIN spatial_rf_predictions srf
        ON srf.run_id = rac.run_id
       AND srf.catalogNumber = rac.catalogNumber
       AND srf.plant_individual_id = rac.plant_individual_id
    LEFT JOIN environmental_cache ec
        ON ec.lat_lon_key = COALESCE(
            srf.lat_lon_key,
            CASE
                WHEN COALESCE(v.latitude, v.decimalLatitude) IS NOT NULL
                 AND COALESCE(v.longitude, v.decimalLongitude) IS NOT NULL
                THEN printf('%.5f_%.5f', COALESCE(v.latitude, v.decimalLatitude), COALESCE(v.longitude, v.decimalLongitude))
                ELSE NULL
            END,
            (SELECT ec2.lat_lon_key FROM environmental_cache ec2 WHERE ec2.catalogNumber = rac.catalogNumber LIMIT 1)
        )
    LEFT JOIN curation_triage ct
        ON ct.run_id = rac.run_id
       AND ct.catalogNumber = rac.catalogNumber
       AND ct.plant_individual_id = rac.plant_individual_id;
    """


def _build_schema_ddl() -> str:
    harm_ddl = ",\n        ".join(f"{h} REAL" for h in HARMONIC_COLS)
    pc_ddl = ",\n        ".join(f"{p} REAL" for p in PC_COLS)
    sym_pc_ddl = ",\n        ".join(f"{sp} REAL" for sp in SYM_PC_COLS)
    wc_srf_ddl = ",\n        ".join(
        f"{col} REAL"
        for col in (*WC_BIO_COLUMNS, "bio18_precip_warmest_qt", "bio19_precip_coldest_qt")
    )

    env_col_defs: List[str] = ["lat_lon_key TEXT PRIMARY KEY"]
    for col in ENVIRONMENTAL_CACHE_COLUMNS:
        if col == "lat_lon_key":
            continue
        if col in ENVIRONMENTAL_TEXT_COLS:
            env_col_defs.append(f"{col} TEXT")
        elif col in ENVIRONMENTAL_BOOL_COLS:
            env_col_defs.append(f"{col} INTEGER")
        else:
            env_col_defs.append(f"{col} REAL")
    env_ddl = ",\n        ".join(env_col_defs)

    master_view_ddl = _build_master_view_ddl()

    return f"""
    -- =========================================================================
    -- TIER A: CANONICAL SHARED STORE
    -- =========================================================================

    CREATE TABLE IF NOT EXISTS taxon_synonyms (
        raw_name_lower TEXT PRIMARY KEY,
        canonical_packera TEXT NOT NULL,
        protologue_authors_json TEXT,
        updated_at TEXT
    );

    CREATE TABLE IF NOT EXISTS vouchers (
        catalogNumber TEXT PRIMARY KEY,
        institutionCode TEXT,
        scientificName TEXT,
        species_raw TEXT,
        recordedBy TEXT,
        recordNumber TEXT,
        identifiedBy TEXT,
        determiner_raw TEXT,
        determiner_tier TEXT,
        type_status TEXT,
        county TEXT,
        stateProvince TEXT,
        locality TEXT,
        decimalLatitude REAL,
        latitude REAL,
        decimalLongitude REAL,
        longitude REAL,
        coordinateUncertainty REAL,
        year INTEGER,
        month INTEGER,
        day INTEGER,
        eventDate TEXT,
        doy REAL,
        pheno_sin REAL,
        pheno_cos REAL,
        regional_group TEXT,
        image_path TEXT,
        exsiccatae_key TEXT,
        is_primary_duplicate INTEGER,
        duplicate_count INTEGER,
        exsiccatae_taxonomic_conflict INTEGER DEFAULT 0,
        exsiccatae_gold_taxon TEXT,
        species_standardized TEXT,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        qc_status TEXT DEFAULT 'PENDING' CHECK(qc_status IN ('PASS', 'FAILED_QC', 'PENDING')),
        failure_reason TEXT,
        failure_details TEXT
    );

    CREATE INDEX IF NOT EXISTS idx_vouchers_species_std ON vouchers(species_standardized);
    CREATE INDEX IF NOT EXISTS idx_vouchers_qc_status ON vouchers(qc_status);
    CREATE INDEX IF NOT EXISTS idx_vouchers_det_tier ON vouchers(determiner_tier);

    CREATE TABLE IF NOT EXISTS plant_individuals (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        rosette_centroid_x REAL,
        rosette_centroid_y REAL,
        crown_x REAL,
        crown_y REAL,
        plant_height_px REAL,
        plant_height_mm REAL,
        plant_height_cm REAL,
        is_bent_stem INTEGER,
        petiole_length_px REAL,
        petiole_length_mm REAL,
        blade_length_px REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        scape_to_basal_leaf_ratio REAL,
        rosette_patch_path TEXT,
        rosette_mask_path TEXT,
        PRIMARY KEY (catalogNumber, plant_individual_id),
        FOREIGN KEY (catalogNumber) REFERENCES vouchers(catalogNumber) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS extracted_leaves (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        leaf_id INTEGER NOT NULL,
        shape_id TEXT UNIQUE,
        tier TEXT,
        assigned_tier TEXT NOT NULL CHECK(assigned_tier IN ('Tier_1_Direct', 'Tier_2_Reflected')),
        tier_assignment TEXT,
        ucs REAL,
        ucs_score REAL,
        solidity REAL,
        circularity REAL,
        extent REAL,
        aspect_ratio REAL,
        area_px REAL,
        blade_area_cm2 REAL,
        is_folded INTEGER,
        is_dissected INTEGER,
        reflection_applied INTEGER,
        midrib_angle_deg REAL,
        petiole_length_px REAL,
        petiole_length_mm REAL,
        blade_length_px REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        mask_path TEXT,
        crop_path TEXT,
        contour_path TEXT,
        contour_n_pts INTEGER DEFAULT 120,
        contour_xy_json TEXT NOT NULL,
        PRIMARY KEY (catalogNumber, plant_individual_id, leaf_id),
        FOREIGN KEY (catalogNumber, plant_individual_id)
            REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
    );

    CREATE INDEX IF NOT EXISTS idx_extracted_leaves_tier ON extracted_leaves(assigned_tier);
    CREATE INDEX IF NOT EXISTS idx_extracted_leaves_cat ON extracted_leaves(catalogNumber);
    CREATE INDEX IF NOT EXISTS idx_extracted_leaves_cat_pid_tier
        ON extracted_leaves(catalogNumber, plant_individual_id, assigned_tier);

    CREATE TABLE IF NOT EXISTS voucher_reproductive_phenology (
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        capitula_count INTEGER,
        involucre_height_px REAL,
        involucre_width_px REAL,
        involucre_height_mm REAL,
        involucre_width_mm REAL,
        AR_inv REAL,
        capitulum_aspect_ratio REAL,
        total_capitula INTEGER,
        n_anthesis INTEGER,
        n_bud INTEGER,
        n_fruiting INTEGER,
        n_fruit INTEGER,
        n_senescent INTEGER,
        voucher_phenological_state TEXT,
        is_verified_anthesis INTEGER,
        reproductive_stage_confidence REAL,
        dominant_pappus_ratio REAL,
        PRIMARY KEY (catalogNumber, plant_individual_id),
        FOREIGN KEY (catalogNumber, plant_individual_id)
            REFERENCES plant_individuals(catalogNumber, plant_individual_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS environmental_cache (
        {env_ddl}
    );

    CREATE INDEX IF NOT EXISTS idx_environmental_cache_catalog
        ON environmental_cache(catalogNumber);

    CREATE VIRTUAL TABLE IF NOT EXISTS environmental_rtree USING rtree(
        id,
        min_lon, max_lon,
        min_lat, max_lat
    );

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_insert
    AFTER INSERT ON environmental_cache
    WHEN NEW.latitude IS NOT NULL AND NEW.longitude IS NOT NULL
    BEGIN
        INSERT OR REPLACE INTO environmental_rtree(id, min_lon, max_lon, min_lat, max_lat)
        VALUES (NEW.rowid, NEW.longitude, NEW.longitude, NEW.latitude, NEW.latitude);
    END;

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_update
    AFTER UPDATE ON environmental_cache
    BEGIN
        DELETE FROM environmental_rtree WHERE id = OLD.rowid;
        INSERT INTO environmental_rtree(id, min_lon, max_lon, min_lat, max_lat)
        SELECT NEW.rowid, NEW.longitude, NEW.longitude, NEW.latitude, NEW.latitude
        WHERE NEW.latitude IS NOT NULL AND NEW.longitude IS NOT NULL;
    END;

    CREATE TRIGGER IF NOT EXISTS trg_env_cache_rtree_delete
    AFTER DELETE ON environmental_cache
    BEGIN
        DELETE FROM environmental_rtree WHERE id = OLD.rowid;
    END;

    -- =========================================================================
    -- TIER B: RUN-SCOPED ANALYTICAL TABLES & MASTER VIEW
    -- =========================================================================

    CREATE TABLE IF NOT EXISTS pipeline_runs (
        run_id TEXT PRIMARY KEY,
        run_dir TEXT,
        start_timestamp TEXT,
        end_timestamp TEXT,
        git_commit_sha TEXT,
        cli_arguments_json TEXT,
        focal_taxon TEXT,
        intraspecific_mode INTEGER DEFAULT 0
    );

    CREATE TABLE IF NOT EXISTS run_active_cohort (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        species_standardized TEXT,
        determiner_tier TEXT,
        type_status TEXT,
        is_type_specimen INTEGER DEFAULT 0,
        is_promoted_bronze INTEGER DEFAULT 0,
        authority_weight REAL,
        spatial_balance_weight REAL DEFAULT 1.0,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE,
        FOREIGN KEY (catalogNumber) REFERENCES vouchers(catalogNumber) ON DELETE CASCADE
    );

    CREATE INDEX IF NOT EXISTS idx_run_active_cohort_run ON run_active_cohort(run_id);

    CREATE TABLE IF NOT EXISTS leaf_efa_harmonics (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        leaf_id INTEGER NOT NULL,
        assigned_tier TEXT,
        solidity REAL,
        dissection_index REAL,
        ucs_score REAL,
        aspect_ratio REAL,
        area_px REAL,
        is_dissected INTEGER,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        petiole_length_mm REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        {harm_ddl},
        sym_var REAL,
        asym_var REAL,
        pct_asymmetry REAL,
        harmonic_roughness REAL,
        is_contour_outlier INTEGER,
        leaf_quality_weight REAL,
        centroid_size REAL,
        log_centsize REAL,
        efd_harmonic_opt INTEGER,
        efd_k95_pcs INTEGER,
        {pc_ddl},
        {sym_pc_ddl},
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id, leaf_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS otu_morphometrics (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        assigned_tier TEXT,
        species_standardized TEXT,
        determiner_tier TEXT,
        is_promoted_anchor INTEGER,
        authority_weight REAL,
        is_passive INTEGER,
        pixels_per_mm REAL,
        is_calibrated_mm INTEGER,
        solidity REAL,
        dissection_index REAL,
        q25_solidity REAL,
        dissected_leaf_fraction REAL,
        basal_leaf_count INTEGER,
        ucs_score REAL,
        aspect_ratio REAL,
        area_px REAL,
        is_dissected INTEGER,
        petiole_length_mm REAL,
        blade_length_mm REAL,
        petiole_to_blade_ratio REAL,
        centroid_size REAL,
        log_centsize REAL,
        AR_inv REAL,
        AR_inv_std REAL,
        reproductive_trait_status TEXT,
        {harm_ddl},
        sym_var REAL,
        asym_var REAL,
        pct_asymmetry REAL,
        harmonic_roughness REAL,
        is_contour_outlier INTEGER,
        leaf_quality_weight REAL,
        leaf_count INTEGER,
        n_tier1_leaves INTEGER,
        n_tier2_leaves INTEGER,
        morph_tier_source TEXT,
        foliar_variance REAL,
        heteroblasty_index REAL,
        efd_harmonic_opt INTEGER,
        efd_k95_pcs INTEGER,
        selected_morph_features TEXT,
        {pc_ddl},
        {sym_pc_ddl},
        gmm_cluster TEXT,
        gmm_uncertainty REAL,
        cda_predicted_taxon TEXT,
        cda_posterior_prob REAL,
        cda_second_taxon TEXT,
        cda_second_posterior_prob REAL,
        cda_posterior_margin REAL,
        loocv_accuracy REAL,
        can1 REAL,
        can2 REAL,
        can3 REAL,
        ecoregion_cda_predicted TEXT,
        ecoregion_cda_posterior_prob REAL,
        ecoregion_can1 REAL,
        ecoregion_can2 REAL,
        morphometric_outlier_flag INTEGER,
        misidentification_flag INTEGER,
        triage_priority TEXT,
        discordance_reason TEXT,
        type_status TEXT,
        is_type_specimen INTEGER DEFAULT 0,
        has_empirical_type_anchor INTEGER,
        type_morph_mahal_d2 REAL,
        nearest_type_morph_taxon TEXT,
        type_morph_concordance REAL,
        is_typical_morphotype INTEGER,
        min_congener_mahal_d2 REAL,
        min_congener_taxon TEXT,
        open_set_novelty_p REAL,
        is_confident_label_noise INTEGER,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS spatial_rf_predictions (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        lat_lon_key TEXT,
        rf_predicted_taxon TEXT,
        rf_predicted_prob REAL,
        rf_given_taxon_prob REAL,
        edaphic_best_fit_taxon TEXT,
        spatial_block_cv_acc REAL,
        spatial_autocorr_range_km REAL,
        spatial_block_id TEXT,
        doy REAL,
        expected_doy REAL,
        hopkins_expected REAL,
        delta_doy REAL,
        hopkins_adjusted_doy REAL,
        pheno_timing_category TEXT,
        anthesis_filter_status TEXT,
        edaphic_z_score REAL,
        climatic_z_score REAL,
        type_edaphic_z_score REAL,
        nearest_type_niche_taxon TEXT,
        matches_type_lithology INTEGER,
        primary_env_predictors TEXT,
        primary_clim_predictors TEXT,
        primary_edaph_predictors TEXT,
        {wc_srf_ddl},
        substrate_specialist_flag INTEGER,
        piedmont_outcrop_specialist_flag INTEGER,
        multimodal_concordance TEXT,
        multimodal_conflict_flag INTEGER,
        triage_category TEXT,
        conflict_discordance_notes TEXT,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    CREATE TABLE IF NOT EXISTS curation_triage (
        run_id TEXT NOT NULL,
        catalogNumber TEXT NOT NULL,
        plant_individual_id INTEGER NOT NULL DEFAULT 0,
        w_morph_kappa REAL,
        w_vision_kappa REAL,
        w_edaphic_kappa REAL,
        w_pheno_kappa REAL,
        n_active_modalities INTEGER,
        n_available_modalities INTEGER,
        probabilistic_consensus_conflict REAL,
        composite_conflict_score REAL,
        probabilistic_consensus_score REAL,
        hopkins_adjusted_doy REAL,
        delta_doy_resid REAL,
        taxonomic_status_call TEXT,
        recommended_determination TEXT,
        taxonomic_treatment TEXT,
        taxonomic_treatment_broad TEXT,
        taxonomic_treatment_subtype TEXT,
        synthesis_triage_priority TEXT,
        synthesis_triage_action TEXT,
        synthesis_rationale TEXT,
        cryptic_cluster_cohesion_pass INTEGER,
        is_completed_voucher INTEGER,
        PRIMARY KEY (run_id, catalogNumber, plant_individual_id),
        FOREIGN KEY (run_id) REFERENCES pipeline_runs(run_id) ON DELETE CASCADE
    );

    {master_view_ddl}
    """


MASTER_VIEW_DDL: str = _build_master_view_ddl()
SCHEMA_DDL: str = _build_schema_ddl()

_VOUCHERS_MIGRATION_COLUMNS: Tuple[Tuple[str, str], ...] = (
    ("exsiccatae_taxonomic_conflict", "INTEGER DEFAULT 0"),
    ("exsiccatae_gold_taxon", "TEXT"),
)

_DYNAMIC_TABLE_MIGRATION_SPECS: Dict[str, Tuple[Tuple[str, str], ...]] = {
    "run_active_cohort": (
        ("type_status", "TEXT"),
        ("is_type_specimen", "INTEGER DEFAULT 0"),
    ),
    "leaf_efa_harmonics": tuple(
        [
            *((h, "REAL") for h in HARMONIC_COLS),
            *((p, "REAL") for p in PC_COLS),
            *((sp, "REAL") for sp in SYM_PC_COLS),
            ("efd_harmonic_opt", "INTEGER"),
            ("efd_k95_pcs", "INTEGER"),
            ("dissection_index", "REAL"),
        ]
    ),
    "otu_morphometrics": tuple(
        [
            *((h, "REAL") for h in HARMONIC_COLS),
            *((p, "REAL") for p in PC_COLS),
            *((sp, "REAL") for sp in SYM_PC_COLS),
            ("efd_harmonic_opt", "INTEGER"),
            ("efd_k95_pcs", "INTEGER"),
            ("selected_morph_features", "TEXT"),
            ("dissection_index", "REAL"),
            ("q25_solidity", "REAL"),
            ("dissected_leaf_fraction", "REAL"),
            ("basal_leaf_count", "INTEGER"),
            ("type_status", "TEXT"),
            ("is_type_specimen", "INTEGER DEFAULT 0"),
            ("has_empirical_type_anchor", "INTEGER"),
            ("type_morph_mahal_d2", "REAL"),
            ("nearest_type_morph_taxon", "TEXT"),
            ("type_morph_concordance", "REAL"),
            ("is_typical_morphotype", "INTEGER"),
            ("min_congener_mahal_d2", "REAL"),
            ("min_congener_taxon", "TEXT"),
            ("open_set_novelty_p", "REAL"),
            ("is_confident_label_noise", "INTEGER"),
        ]
    ),
    "spatial_rf_predictions": tuple(
        [
            *((col, "REAL") for col in WC_BIO_COLUMNS),
            ("bio18_precip_warmest_qt", "REAL"),
            ("bio19_precip_coldest_qt", "REAL"),
            ("climatic_z_score", "REAL"),
            ("type_edaphic_z_score", "REAL"),
            ("nearest_type_niche_taxon", "TEXT"),
            ("matches_type_lithology", "INTEGER"),
            ("primary_env_predictors", "TEXT"),
            ("primary_clim_predictors", "TEXT"),
            ("primary_edaph_predictors", "TEXT"),
        ]
    ),
    "curation_triage": (
        ("cryptic_cluster_cohesion_pass", "INTEGER"),
    ),
    "environmental_cache": tuple(
        (
            col,
            "TEXT"
            if col in ENVIRONMENTAL_TEXT_COLS
            else ("INTEGER" if col in ENVIRONMENTAL_BOOL_COLS else "REAL"),
        )
        for col in ENVIRONMENTAL_CACHE_COLUMNS
        if col != "lat_lon_key"
    ),
}

_MASTER_VIEW_REQUIRED_TABLES: Tuple[str, ...] = (
    "run_active_cohort",
    "vouchers",
    "plant_individuals",
    "voucher_reproductive_phenology",
    "otu_morphometrics",
    "spatial_rf_predictions",
    "environmental_cache",
    "curation_triage",
)


def _migrate_vouchers_schema_columns(conn: sqlite3.Connection) -> None:
    """Idempotently add optional exsiccatae conflict audit columns to existing `vouchers` tables."""
    existing_cols = {
        str(r[1]) for r in conn.execute("PRAGMA table_info(vouchers);").fetchall()
    }
    if not existing_cols:
        return
    for col_name, col_ddl in _VOUCHERS_MIGRATION_COLUMNS:
        if col_name not in existing_cols:
            conn.execute(f"ALTER TABLE vouchers ADD COLUMN {col_name} {col_ddl};")
            existing_cols.add(col_name)


PACKERA_SCHEMA_USER_VERSION: int = 20260407


def _migrate_dynamic_schema_columns(
    conn: sqlite3.Connection,
    recreate_view: bool = True,
) -> None:
    """
    Idempotently upgrade existing SQLite databases in-place:
      1. Inspects `PRAGMA table_info(...)` on `leaf_efa_harmonics`, `otu_morphometrics`,
         `spatial_rf_predictions`, and `environmental_cache`.
      2. Issues `ALTER TABLE ... ADD COLUMN ...` for any missing dynamic harmonic (`A1..D24`),
         PCA (`PC1..PC15`, `Sym_PC1..Sym_PC15`), adaptive EFD metadata (`efd_harmonic_opt`,
         `efd_k95_pcs`, `selected_morph_features`), dissection scalars, 19-variable WorldClim
         (`bio1..bio19`), `climatic_z_score`, or `primary_*_predictors` columns.
      3. Backfills bidirectional `bio18_precip_warmest_qtr <-> bio18_precip_warmest_qt` and
         `bio19_precip_coldest_qtr <-> bio19_precip_coldest_qt` values on existing tables.
      4. Optionally drops and recreates `v_integrative_triage_master` whenever all base tables exist.
    """
    existing_tables: Set[str] = {
        str(r[0])
        for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table';"
        ).fetchall()
    }
    if not existing_tables:
        return

    conn.execute("DROP VIEW IF EXISTS v_integrative_triage_master;")

    for table_name, col_specs in _DYNAMIC_TABLE_MIGRATION_SPECS.items():
        if table_name not in existing_tables:
            continue
        existing_cols = {
            str(r[1]) for r in conn.execute(f"PRAGMA table_info({table_name});").fetchall()
        }
        if not existing_cols:
            continue
        for col_name, col_ddl in col_specs:
            if col_name not in existing_cols:
                conn.execute(f"ALTER TABLE {table_name} ADD COLUMN {col_name} {col_ddl};")
                existing_cols.add(col_name)

        if {
            "bio18_precip_warmest_qtr",
            "bio18_precip_warmest_qt",
            "bio19_precip_coldest_qtr",
            "bio19_precip_coldest_qt",
        }.issubset(existing_cols):
            conn.execute(
                f"""
                UPDATE {table_name}
                SET bio18_precip_warmest_qtr = COALESCE(bio18_precip_warmest_qtr, bio18_precip_warmest_qt),
                    bio18_precip_warmest_qt = COALESCE(bio18_precip_warmest_qt, bio18_precip_warmest_qtr),
                    bio19_precip_coldest_qtr = COALESCE(bio19_precip_coldest_qtr, bio19_precip_coldest_qt),
                    bio19_precip_coldest_qt = COALESCE(bio19_precip_coldest_qt, bio19_precip_coldest_qtr)
                WHERE (bio18_precip_warmest_qtr IS NULL AND bio18_precip_warmest_qt IS NOT NULL)
                   OR (bio18_precip_warmest_qt IS NULL AND bio18_precip_warmest_qtr IS NOT NULL)
                   OR (bio19_precip_coldest_qtr IS NULL AND bio19_precip_coldest_qt IS NOT NULL)
                   OR (bio19_precip_coldest_qt IS NULL AND bio19_precip_coldest_qtr IS NOT NULL);
                """
            )

    if recreate_view and set(_MASTER_VIEW_REQUIRED_TABLES).issubset(existing_tables):
        conn.executescript(MASTER_VIEW_DDL)


def init_db_schema(conn: sqlite3.Connection) -> None:
    """Create or idempotently migrate all Tier A and Tier B tables, indexes, triggers, R*Tree, and master view."""
    cur_ver_row = conn.execute("PRAGMA user_version;").fetchone()
    cur_ver = int(cur_ver_row[0]) if (cur_ver_row and cur_ver_row[0] is not None) else 0
    if cur_ver == PACKERA_SCHEMA_USER_VERSION:
        existing_objs = {
            str(r[0])
            for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type IN ('table', 'view');"
            ).fetchall()
        }
        required_objs = set(_MASTER_VIEW_REQUIRED_TABLES) | {
            "taxon_synonyms",
            "extracted_leaves",
            "pipeline_runs",
            "leaf_efa_harmonics",
            "v_integrative_triage_master",
        }
        if required_objs.issubset(existing_objs):
            conn.execute("PRAGMA foreign_keys = ON;")
            return

    _migrate_vouchers_schema_columns(conn)
    _migrate_dynamic_schema_columns(conn, recreate_view=False)
    conn.executescript(SCHEMA_DDL)
    conn.execute(f"PRAGMA user_version = {PACKERA_SCHEMA_USER_VERSION};")
    conn.execute("PRAGMA foreign_keys = ON;")
    conn.commit()
