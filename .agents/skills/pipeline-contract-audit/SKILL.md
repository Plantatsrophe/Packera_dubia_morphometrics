---
name: pipeline-contract-audit
description: Audits cross-stage Python and R CSV schema contracts, checks for forbidden _archive or synthetic fallback leaks, and runs the test suite. Use before committing code changes or when debugging pipeline stage handshakes.
---

# Pipeline Contract & Schema Audit (`UNC-BOT-SOP-2026-04-REV5`)

Use this skill before committing code changes, after modifying any pipeline stage (`01` through `08`), or when diagnosing missing-column / join-key errors across Python and R scripts.

---

## Step 1: Scan Active `.py` and `.R` Files for Forbidden `_archive` & Synthetic Leaks

Active code must never import from, read from, or fall back to `scripts/_archive/` or `data/_archive/`, and must never silently synthesize random embeddings, Gaussian heatmaps, or trigonometric environmental variables outside an explicit `--demo` flag.

Run the following verification script from the repository root:

```bash
python3 - <<'EOF'
import re
import sys
from pathlib import Path

root = Path(".").resolve()
active_files = [
    p for p in list(root.rglob("*.py")) + list(root.rglob("*.R"))
    if "_archive" not in p.parts
    and ".venv" not in p.parts
    and ".venv_LM2" not in p.parts
    and "LeafMachine2" not in p.parts
]

archive_violations = []
for f in sorted(active_files):
    rel = f.relative_to(root)
    # Skip test files that explicitly assert _archive isolation
    if rel.parts[:2] == ("scripts", "tests"):
        continue
    for lineno, line in enumerate(f.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if "_archive" in line:
            archive_violations.append(f"{rel}:{lineno}: {stripped}")

if archive_violations:
    print("[FAIL] Forbidden _archive references found in active code:")
    for v in archive_violations:
        print("  -", v)
    sys.exit(1)
else:
    print(f"[PASS] Scanned {len(active_files)} active .py/.R files: 0 _archive leaks detected.")
EOF
```

---

## Step 2: Audit Cross-Stage Python & R CSV Schema Contracts

Verify that all production stages (`02` through `08`) agree on primary keys, column names, join semantics, and dynamic reporting:

### Core Orchestration & Grouping Invariants

* **Composite Specimen Key Invariant `(catalogNumber, plant_individual_id)`:** Every leaf-level and specimen-level join and aggregation across Steps 02, 03, 04, 05, 06, 07, and 08 must group by `(catalogNumber, plant_individual_id)`. Multi-plant herbarium sheets (`plant_individual_id = 0, 1, ...`) must never be collapsed to `plant_individual_id = 0` or grouped by `catalogNumber` alone.
* **Orchestrator CLI Handoffs (`main.py` & `run_pipeline.sh`):**
  * `main.py` passes `--reproductive` (plus allometry/batch-audit outputs, `--focal-taxon`, and `--intraspecific-mode`) to Step 04, `--manifest-csv` and `--morphometrics-csv` to Step 05, and `--pheno-states` (plus SSURGO, lithology contingency, and phenology outputs) to Step 06.
  * `run_pipeline.sh` passes `--manifest-csv "${RUN_DIR}/data/tables/extracted_leaf_manifest.csv"` to Step 05 and executes Step 08 (`08_discover_taxon_complexes.R`) **after** Step 06 (`06_multimodal_spatial_rf.R`) with `--multimodal "${RUN_DIR}/data/tables/multimodal_conflict_flags.csv"`.

### Canonical Stage-to-Stage Handshake Matrix

| Stage Script | Relational SQLite Tables / Views | Synchronized CSV Outputs | Required Keys & Shared Columns | Critical Handshake Invariants |
| :--- | :--- | :--- | :--- | :--- |
| **`scripts/core/db.py` & `scripts/core/db_utils.R`** | `data/packera_pipeline.sqlite` (All Tier A & Tier B tables, triggers, indexes, and `v_integrative_triage_master`) | Exported via `export_table_or_view_to_csv_atomic()` | `(catalogNumber, plant_individual_id)` | Core database interface. Configures SQLite WAL mode (`PRAGMA journal_mode = WAL; PRAGMA synchronous = NORMAL; PRAGMA foreign_keys = ON;`), establishes transactional integrity across Python (`sqlite3`) and R (`DBI`/`RSQLite`), and manages the write-then-replace atomic CSV snapshot bridge (`.tmp` + `os.fsync()` + `Path.replace()`). |
| **`scripts/data_prep/01_voucher_harvester.py`** | Writes `taxon_synonyms`, `vouchers` (PK: `catalogNumber`, atomic UPSERT) | `data/tables/curated_vouchers.csv` | `catalogNumber`, `institutionCode`, `scientificName`, `species_raw`, `recordedBy`, `recordNumber`, `identifiedBy`, `determiner_raw`, `determiner_tier`, `type_status`, `county`, `stateProvince`, `decimalLatitude`, `latitude`, `decimalLongitude`, `longitude`, `coordinateUncertainty`, `year`, `month`, `day`, `eventDate`, `regional_group`, `image_path`, `exsiccatae_key`, `is_primary_duplicate`, `duplicate_count`, `exsiccatae_taxonomic_conflict`, `exsiccatae_gold_taxon` | 26-column Darwin Core schema. Dynamic taxon standardizer handles all 24+ Eastern North American *Packera* taxa, basionyms, and continent-wide species via `expand_taxon_synonyms_gbif()`. Running with zero taxon flags dynamically queries and harvests all *Packera* species across the entire genus (`discover_genus_packera_species_gbif()`, `target_species: []`), whereas passing `--taxon "<sp1>" "<sp2>" ...` (or `--focal-taxon "<sp>"`) restricts analysis strictly to user-specified species with zero hardcoded defaults or injected outgroups. Configurable geographic filtering via `--geographic-scope {eastern_na,north_america,global}`, `--include-western`, and `--countries`. Evaluates points against official US EPA Level III Ecoregions (`data/environmental/epa_ecoregions_l3/`, Omernik & Griffith 2014); for ungeoreferenced sheets, assigns ecoregion via county polygon containment ($\ge 95\%$ in single ecoregion) or marks `NA`. 3-layer dynamic authority engine (`SPECIALIST_PATTERNS`, `extract_protologue_authors()`, `discover_corpus_authorities()`) stratifies into active reference anchors (`Tier_1_Gold` & `Tier_2_Silver`, `authority_weight = 1.0` or Dawid-Skene precision) vs. passively projected (`Tier_3_Bronze`, `authority_weight = 0.0`). `parse_determiner_tier()` restricts `SPECIALIST_PATTERNS` matching strictly to `identifiedBy` and `identificationHistory` (only falling back to `recordedBy` when `identifiedBy` is unrecorded or indicates the collector via `_should_fallback_to_collector()`), preventing collector-as-determiner false Gold promotions. `stratify_duplicates()` performs cross-duplicate exsiccatae conflict auditing (`exsiccatae_taxonomic_conflict`, `exsiccatae_gold_taxon`), flagging conflicting determinations across multi-herbarium sets. Layer 3 `discover_corpus_authorities()` enforces a cross-duplicate concordance guard, rejecting determiners with $> 25\%$ conflicting duplicate annotations (`det_conflicted_groups / det_dup_groups > 0.25`). |
| **`scripts/pipeline/02_segment_and_extract.py`** | Reads `vouchers`; writes `plant_individuals` (PK: `(catalogNumber, plant_individual_id)`), `extracted_leaves` (PK: `(catalogNumber, plant_individual_id, leaf_id)`, with `contour_xy_json`), `voucher_reproductive_phenology` (PK: `(catalogNumber, plant_individual_id)`) | `data/tables/extracted_leaf_manifest.csv`, `data/tables/failed_qc_vouchers.csv`, `data/tables/voucher_reproductive_metrics.csv`, `data/tables/voucher_phenological_states.csv`, `data/contours/*.csv`, `data/masks/{tier1_intact,tier2_reflected}/*.png`, `data/cropped_patches/{basal_leaves_raw,capitula,rosettes_dense}/` | `catalogNumber`, `plant_individual_id`, `leaf_id`, `assigned_tier`, `pixels_per_mm`, `is_calibrated_mm`, `petiole_length_px`, `petiole_length_mm`, `plant_height_px`, `plant_height_mm`, `is_bent_stem`, `contour_path`, `mask_path` | Stock LeafMachine2 v2.2+ (`PLANT_LeafPriority`, `minimum_confidence_threshold: 0.65`, `treat_leaflet_as_leaf_whole: False`, `use_efds_for_png_masks: False`, `do_apply_landmarks: False`); Stage 1 YOLOv5 processes `Class 0` and `Class 1` boxes with high-DPI upper cauline/tiny fragment pre-filtering; Otsu boundary-continuity lobe rejection (`is_truncated_lobe_or_leaflet`, `max_tissue_border_ratio: 0.15`); ruler-calibrated physical area floor (`min_blade_area_cm2: 0.75`); classical dimensionless descriptors (`Solidity >= 0.72`, `Circularity >= 0.25`, `Extent >= 0.35`) with 100% label-blind pinnatifid ($0.48\text{–}0.72$) and bipinnatisect ($0.18\text{–}0.48$) dissection gates via `is_botanical_dissection()`; resolution-adaptive `w_span` and petiole width thresholds; spatial `rosette_centroids_by_pid` wired via `rosette_centroid`; estimates crown point $p_{\text{crown}}$, per-leaf petiole length $L_{\text{petiole}} = \max(d_{\text{geodesic}}(p_{\text{base}}, p_{\text{crown}}), \|p_{\text{base}} - p_{\text{crown}}\|_2)$ aggregated across outer-rosette leaves at $P_{75}$, straight vs. bent scape height $H_{\text{plant}}$, and regenerates post-QC montage via `--regenerate-figures-only`. Per-voucher atomic SQLite transactions commit to `plant_individuals`, `extracted_leaves`, and `voucher_reproductive_phenology` with $O(1)$ SQL resumption checks (`qc_status IN ('PASS', 'FAILED_QC')`). `crop_dense_rosette_tier4()` executes strictly after cauline/stem rejection (`is_cauline_or_stem_contaminant()`) and the `Tier_1_Direct` gatekeeper, replacing non-basal-leaf pixels outside the 5-pixel dilated union of verified basal leaf masks with neutral mounting paper RGB `[242, 239, 230]`, saving companion binary mask `<catalogNumber>_p<pid>_rosette_mask.png` in `data/cropped_patches/rosettes_dense/`. |
| **`scripts/morphometrics/03_fourier_extractor.R`** | Reads `extracted_leaves` (single-query batch SQL contour ingestion), `vouchers`, `plant_individuals`; writes `leaf_efa_harmonics` | `data/tables/fourier_pcs.csv`, `data/tables/leaf_efa_harmonics.csv`, `data/tables/leaf_efa_harmonics_individual.csv` | `(catalogNumber, plant_individual_id)`, `assigned_tier`, `reconstruction_tier`, `pixels_per_mm`, `is_calibrated_mm`, `centroid_size`, `log_centsize`, `petiole_length_mm`, `blade_length_mm`, `petiole_to_blade_ratio`, `is_dissected`, `institutionCode`, `stateProvince`, `county`, `locality`, `decimalLatitude`, `decimalLongitude`, `coordinateUncertainty`, `year`, `month`, `day`, `eventDate`, `authority_tier`, `type_status`, `PC1`–`PC5`, `Sym_PC1`–`Sym_PC5`, `sym_var`, `asym_var`, `pct_asymmetry`, `harmonic_roughness`, `is_contour_outlier`, `leaf_quality_weight`, `n_tier1_leaves`, `n_tier2_leaves`, `morph_tier_source`, `A1`–`D12` | Strict `Tier_1_Direct` closed-contour restriction (never pooling `Tier_2_Reflected` mirrored hemi-blades or Tier 3 open curves into EFA); batch SQL contour ingestion from `extracted_leaves.contour_xy_json` eliminating 14,000+ CSV reads (with directory fallback for unit tests); inertial axis pre-alignment (`Momocs::coo_alignxax` with circularity $> 0.75$ suborbicular guard); $\ge 99\%$ cumulative harmonic power calibration at $H = 12$ (`calibrate_harmonicpower_efourier`); uniform `Momocs::rm_asym(ef)` bilateral symmetrization (Claude 2008); MCD Mahalanobis outlier filtering (`MASS::cov.rob`, $p < 0.001$) combined with `pct_asymmetry > 0.22 | harmonic_roughness > 0.15`; unweighted specimen-level (OTU) median harmonic and PC aggregation to `(catalogNumber, plant_individual_id)` (`morph_tier_source = "Tier_1_Pure"` when `n_tier1_leaves >= 1`), carrying forward Step 02 linear petiole/blade traits and full Step 01 DwC voucher metadata. Sourced with shared `scripts/morphometrics/authority_weighting_utils.R`. |
| **`scripts/morphometrics/04_gmm_morphotools.R`** | Reads `leaf_efa_harmonics`, `vouchers`, `plant_individuals`; writes `otu_morphometrics` | `data/tables/morphometric_flags.csv`, `outputs/reports/gmm_bayes_factors_summary.csv`, `outputs/reports/allometry_homogeneity_test.csv`, `outputs/reports/institutional_batch_effect_audit.csv`, `outputs/reports/ecoregion_supervised_cda_summary.csv`, `outputs/reports/manova_group_significance.csv`, `outputs/reports/manova_pairwise_posthoc.csv` | `(catalogNumber, plant_individual_id)`, `determiner_tier`, `authority_weight`, `species_standardized`, `morph_tier_source`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_taxon`, `cda_second_posterior_prob`, `cda_posterior_margin`, `can1`, `can2`, `can3`, `gmm_cluster`, `gmm_uncertainty`, `morphometric_outlier_flag`, `loocv_accuracy` | Sources `scripts/morphometrics/authority_weighting_utils.R`: adaptive reference anchor selection via `select_adaptive_reference_anchors()` (pruning $\chi^2_{0.95}$ MCD Mahalanobis outliers once to avoid tail truncation; promoting central $\chi^2_{0.80}$ `Tier_3_Bronze` vouchers when a taxon has $< 8$ clean Gold+Silver anchors); 2-step Dawid-Skene / Beta-Binomial Empirical Bayes shrinkage via `compute_empirical_bayes_authority_weights()` ($\alpha_0 = 5.0$, rescaled to mean $1.0$ for active anchors, $0.0$ for unpromoted Bronze); taxon-agnostic discordance strings (`<raw_sp>_Discordant_to_<pred_sp>`) and `"Unassigned_Packera"` fallbacks; Klingenberg (2016) allometric slope homogeneity testing (`PC ~ log_centsize * species_standardized`, applying common-slope correction only when $p \ge 0.01$ on calibrated `is_calibrated_mm == TRUE` anchors); Horn's Parallel Analysis / Broken-Stick and Williams & Titus (1988) $P \le \lfloor N_{\min} / 3 \rfloor$ subspace selection; full-rank $\min(P, g - 1)$ canonical axes for CDA posteriors (`cda_posterior_prob`, `cda_second_posterior_prob`, `cda_posterior_margin`) with `can1`, `can2`, `can3` export; Kass & Raftery (1995) $2\Delta\text{BIC}$ GMM cluster evidence ($> 6$ strong, $> 10$ very strong); Dual Clustering Architecture: Track 1 Unsupervised GMM Morphotype CDA (`fig04_gmm_cda_clustering.{pdf,png}`) and Track 2 Supervised EPA Level III Ecoregion CDA (`fig04c_supervised_ecoregion_cda.{pdf,png}`); parametric omnibus and post-hoc pairwise MANOVA with BH FDR; taxon-adjusted `audit_institutional_batch_effects` (`lm(PC1 ~ species_standardized + institutionCode)`); exports `cda_second_taxon` (secondary morphological parent) alongside `cda_second_posterior_prob` and `cda_posterior_margin`, and aggregates manifest petiole traits by composite key `(catalogNumber, plant_individual_id)`. |
| **`scripts/analysis/05_cleanlab_vision_xai.py`** | Reads `vouchers`, `plant_individuals`, and `dinov2_feature_cache`; writes `dinov2_feature_cache` & `vision_xai_audit` | `data/tables/label_noise_audit.csv`, `outputs/figures/GradCAM_audit_panel.{png,pdf}`, `outputs/figures/fig05b_cleanlab_quality_and_tomentum.{pdf,png}` | `catalogNumber`, `plant_individual_id`, `predicted_label` (`vision_predicted_label`), `confidence_predicted_class`, `c_error`, `is_label_corrupted`, `surface_orientation`, `mean_luminance`, `mean_saturation`, `vision_second_prob`, `vision_margin`, `tomentum_taxon_discordance`, `vision_ood_score`, `vision_is_ood`, `foreground_attention_ratio`, `is_artifact_driven`, `is_background_shortcut` | `load_and_link_rosette_patches()` links each specimen `(catalogNumber, plant_individual_id)` to verified `Tier_1_Direct` basal leaf crops and binary masks via `--manifest-csv data/tables/extracted_leaf_manifest.csv` (falling back to masked `rosettes_dense` crops with companion `<cat>_p<pid>_rosette_mask.png` masks); foliar surface orientation (`classify_foliar_surface_orientation()`), CIE $L^*$ luminance, and HSV saturation are computed strictly inside verified basal leaf lamina mask; indumentum-blind Gaussian low-pass filter ($\sigma = 1.0\text{ px}$ at $224 \times 224$) preventing taphonomic hair-loss bias; 5-px dilated neutral paper masking `[242, 239, 230]`; `dinov2_vitb14_reg` with foreground-masked $16 \times 16$ spatial patch token average pooling (`compute_foreground_patch_mask`); caches 768-d embeddings and `foreground_attention_ratio` in `dinov2_feature_cache` (PK: `(catalogNumber, plant_individual_id, model_name, sigma)`); dynamic `active_taxa` derivation and Okabe-Ito palettes (`_get_dynamic_okabe_ito_map`); anchor parity incorporating `is_promoted_bronze` anchors for Gold/Silver-poor taxa alongside Gold/Silver anchors; `institutionCode`-grouped `StratifiedGroupKFold`, fold-internal `StandardScaler` + `PCA(0.95)` + `LogisticRegression(C=0.1, class_weight="balanced")`, Platt/temperature calibration, and passive projection of unpromoted `Tier_3_Bronze`; class-conditional Cleanlab `prune_by_noise_rate`; conditional abaxial artifact gate in `filter_abaxial_artifacts` suppressing false label noise on arachnoid tomentum unless high-confidence structural misidentification ($\text{conf} \ge 0.85$, $c_{\text{error}} > 0.88$, $\text{FG} \ge 0.65$); feature-space Mahalanobis + cosine $k$-NN OOD score (`vision_ood_score`, `vision_is_ood = TRUE` when `vision_ood_score >= 0.65`); empirical `foreground_attention_ratio` computed for every specimen, flagging `foreground_attention_ratio < 0.65` as `is_artifact_driven = TRUE` and `is_background_shortcut = TRUE`; 3-panel `GradCAM_audit_panel.pdf` masks background to neutral paper `[242, 239, 230]`, restricts attention overlay in Panel `[2]` strictly to `fg_mask`, and plots calibrated posteriors in Panel `[3]`; single-species intraspecific mode evaluates visual discriminability across unsupervised morphometric clusters (`Morphotype_K1..Kn`). |
| **`scripts/analysis/06_multimodal_spatial_rf.R`** | Reads `vouchers`, `plant_individuals`, `otu_morphometrics`, `vision_xai_audit`; queries & writes `environmental_cache` + `environmental_rtree`; writes `spatial_rf_predictions` | `data/tables/multimodal_conflict_flags.csv`, `outputs/reports/multimodal_spatial_rf_summary.csv`, `outputs/reports/edaphic_niche_equivalency.csv`, `outputs/reports/micro_edaphic_outcrop_contingency.csv`, `data/tables/ssurgo_edaphic_validation.csv` | `catalogNumber`, `plant_individual_id`, `misidentification_flag`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_taxon`, `cda_second_posterior_prob`, `cda_posterior_margin`, `Sym_PC1`, `Sym_PC2`, `ecoregion_cda_predicted`, `ecoregion_cda_posterior_prob`, `vision_predicted_label`, `c_error`, `vision_ood_score`, `vision_is_ood`, `foreground_attention_ratio`, `is_artifact_driven`, `surface_orientation`, `tomentum_taxon_discordance`, `edaphic_best_fit_taxon`, `edaphic_z_score`, `multimodal_concordance`, `triage_category`, `spatial_block_cv_acc`, `ssurgo_rock_outcrop_pct`, `ssurgo_bedrock_depth_cm`, `fine_soil_ph`, `fine_soil_sand`, `fine_bedrock_depth_cm`, `topo_tpi_30m`, `topo_twi_30m`, `sgmc_lithology_class`, `substrate_specialist_flag`, `petiole_length_mm`, `blade_length_mm`, `petiole_to_blade_ratio`, `plant_height_cm`, `scape_to_basal_leaf_ratio`, `is_bent_stem`, `capitula_count`, `involucre_height_mm`, `involucre_width_mm`, `AR_inv`, `voucher_phenological_state`, `is_verified_anthesis`, `n_anthesis` | $0\text{–}15\text{ cm}$ trapezoidal rhizosphere depth-weighting (`(5 * v_0_5 + 10 * v_5_15) / 15`) with strict non-blended coalesce hierarchy (100% SSURGO / 30m POLARIS / 100m SOLUS100 in CONUS; 250m SoilGrids strictly non-CONUS fallback when both NA; raw `soil_*` excluded from active `env_vars`), depth-to-bedrock/CFVO/AWC/clr-texture (`clr_sand_silt`, `clr_clay_geom`), 6-class USGS SGMC bedrock lithology (`sgmc_lithology_class`), uncertainty-buffered focal extraction (`100–1000 m`), $|r| < 0.70$ and $\text{VIF} < 5.0$ collinearity pruning (`prune_collinear_predictors`), exhaustive multi-class SGMC and SSURGO outcrop contingency testing with BH FDR (`micro_edaphic_outcrop_contingency.csv`), organic `PCA-env` intraspecific ecotypes (`Ecotype_1_<TopVar>_Low`, `Ecotype_2_<TopVar>_High`), dynamic top-importance RF predictor plotting with 85% confidence ellipses, partitioned Climatic ($D_{\text{clim}}, I_{\text{clim}}$) vs. Edaphic-Geologic ($D_{\text{edaph}}, I_{\text{edaph}}$) Niche Equivalency & Background Similarity tests (`edaphic_niche_equivalency.csv`), active anchor training with Empirical Bayes weights under $0.5^\circ \times 0.5^\circ$ Spatial Block CV (`spatial_block_cv_acc`), `edaphic_best_fit_taxon` pooling Spatial Block-CV out-of-fold RF posteriors with `PCA-env` Mahalanobis likelihoods, `edaphic_z_score` screening physiological envelope violations ($z > 3.0$), and **`dplyr::left_join`** preserving all Step 02 linear/petiole/reproductive/phenology columns, Step 04 morphometric columns, and Step 05 vision XAI/OOD columns. |
| **`scripts/analysis/07_triage_dashboard_synthesis.R`** | Queries master view `v_integrative_triage_master` (left-joining all modalities onto `plant_individuals`); writes `curation_triage` | `data/tables/triage_queue.csv`, `outputs/reports/final_curation_triage.csv`, `outputs/figures/Figure_Integrative_Packera_Revision.pdf` (or `Figure_Integrative_<Focal_Taxon>_Revision.pdf`), `outputs/reports/Packera_Taxonomic_Revision_Summary.md` (or `<Focal_Taxon>_Taxonomic_Revision_Summary.md`) | `(catalogNumber, plant_individual_id)`, `morphometric_outlier_flag` / `misidentification_flag`, `morph_tier_source`, `can1`, `can2`, `can3`, `cda_predicted_taxon`, `cda_posterior_prob`, `cda_second_taxon`, `cda_second_posterior_prob`, `cda_posterior_margin`, `vision_predicted_label`, `taxonomic_status_call`, `recommended_determination`, `taxonomic_treatment_broad`, `taxonomic_treatment_subtype`, `is_completed_voucher`, `synthesis_triage_priority`, `fine_soil_ph`, `fine_soil_sand`, `fine_bedrock_depth_cm`, `ssurgo_rock_outcrop_pct`, `edaphic_z_score`, `rf_predicted_taxon`, `rf_predicted_prob`, `petiole_length_mm`, `blade_length_mm`, `petiole_to_blade_ratio`, `plant_height_cm`, `scape_to_basal_leaf_ratio`, `is_bent_stem`, `capitula_count`, `involucre_height_mm`, `involucre_width_mm`, `AR_inv`, `voucher_phenological_state`, `is_verified_anthesis`, `n_anthesis`, `vision_ood_score`, `vision_is_ood`, `foreground_attention_ratio`, `is_artifact_driven`, `surface_orientation`, `tomentum_taxon_discordance` | Enforces `validate_upstream_contracts()`; coalesces non-georeferenced vouchers achieving 98.4% completion (`5,117 / 5,201` vouchers); out-of-fold Cohen's $\kappa$ and class-conditional precision ($C^{(m)}_{k,j}$) reliability-weighted consensus scoring dynamically renormalized over observed modalities; Hopkins-adjusted circular anthesis phenology (`atan2` residuals controlling for latitude, elevation, and MAT gated on active anthesis); taxon-adjusted `institutionCode` ANCOVA with BH FDR; retains all Step 02, 04, 05, and 06 fine-scale columns in `final_curation_triage.csv`; taxon-agnostic 5-stage hierarchical decision matrix: Stage 1 (`Misidentification_Reassign`) evaluated **BEFORE** `Hybrid_Introgressant` with cross-modal non-contradiction in `morph_or_vis_corrupt_reassign` (vision-supported reassignments strictly gated on `!is_artifact_driven` with `foreground_attention_ratio >= 0.65`); Stage 2 (`Candidate_Taxonomic_Boundary_Zone`) resolves secondary morphological parent via `cda_second_taxon` and nearest-neighbor Mahalanobis distance (`p_b`); Stage 3 (`Putative_Cryptic_Species` via Integration by Congruence) incorporates `vision_is_ood` / `vision_ood_score >= 0.65` as corroborating deep visual divergence; Stage 4 evaluates `is_cohesive_ecotype` for `Subspecies_Ecotype` and `is_taphonomic_plasticity` for `Ecophenotypic_Plasticity`; dynamically generates per-species empirical treatments in `Packera_Taxonomic_Revision_Summary.md` and `Figure_Integrative_Packera_Revision.pdf` (or `<Focal_Taxon>_Taxonomic_Revision_Summary.md` and `Figure_Integrative_<Focal_Taxon>_Revision.pdf` when `--focal-taxon` is specified) with dynamic North American map bounding boxes (`xlim_map`, `ylim_map` padded 12%). |
| **`scripts/analysis/08_discover_taxon_complexes.R`** | Reads `v_integrative_triage_master` / `multimodal_conflict_flags.csv`, `leaf_efa_harmonics`, `vouchers` | `outputs/reports/discovered_taxon_complexes.csv`, `outputs/reports/species_pairwise_multimodal_distances.csv` | `species_standardized`, `complex_id`, `complex_label`, `n_vouchers`, `silhouette_width`, `d_morph`, `d_env`, `d_pheno`, `d_vis`, `d_total` | Accepts `--multimodal data/tables/multimodal_conflict_flags.csv` to join Step 06 fine-scale SSURGO/SoilGrids edaphic and WorldClim variables (`edaphic_schoener_d`, `schoener_d`); unsupervised graph modularity (Leiden/Louvain on $k$-NN similarity network) and Ward.D2 hierarchical clustering across North American *Packera* species to discover objective morpho-ecological species complexes, resolving Taxon x GMM Morphotype subclusters when $N_{\text{taxa}} \le 2$. |

### Automated Contract Verification Script

Run this Python audit to verify column-name agreement and join-key rules across `02`–`07` and the SQLite database interface:

```bash
python3 - <<'EOF'
import sys
from pathlib import Path

scripts_to_check = [
    "scripts/core/db.py",
    "scripts/core/db_utils.R",
    "scripts/core/harvester.py",
    "scripts/pipeline/02_segment_and_extract.py",
    "scripts/morphometrics/03_fourier_extractor.R",
    "scripts/morphometrics/04_gmm_morphotools.R",
    "scripts/morphometrics/authority_weighting_utils.R",
    "scripts/analysis/05_cleanlab_vision_xai.py",
    "scripts/analysis/06_multimodal_spatial_rf.R",
    "scripts/analysis/07_triage_dashboard_synthesis.R",
]

missing = [s for s in scripts_to_check if not Path(s).exists()]
if missing:
    print(f"[FAIL] Missing pipeline scripts: {missing}")
    sys.exit(1)

s_db = Path("scripts/core/db.py").read_text()
s_db_r = Path("scripts/core/db_utils.R").read_text()
s_harv = Path("scripts/core/harvester.py").read_text()
s02 = Path("scripts/pipeline/02_segment_and_extract.py").read_text()
s03 = Path("scripts/morphometrics/03_fourier_extractor.R").read_text()
s04 = Path("scripts/morphometrics/04_gmm_morphotools.R").read_text()
s_auth = Path("scripts/morphometrics/authority_weighting_utils.R").read_text()
s05 = Path("scripts/analysis/05_cleanlab_vision_xai.py").read_text()
s06 = Path("scripts/analysis/06_multimodal_spatial_rf.R").read_text()
s07 = Path("scripts/analysis/07_triage_dashboard_synthesis.R").read_text()

checks = [
    ("db.py defines schema, WAL mode, tables, and v_integrative_triage_master", "v_integrative_triage_master" in s_db and "PRAGMA journal_mode = WAL" in s_db and "dinov2_feature_cache" in s_db),
    ("db_utils.R exposes R connection and master view helpers", "v_integrative_triage_master" in s_db_r and "RSQLite::SQLite()" in s_db_r),
    ("01 parses determiner tier with collector isolation & duplicate conflict guard", "_should_fallback_to_collector" in s_harv and "exsiccatae_taxonomic_conflict" in s_harv and "det_conflicted_groups" in s_harv),
    ("02 exports assigned_tier, pixels_per_mm & rosette_centroids_by_pid", "assigned_tier" in s02 and "pixels_per_mm" in s02 and "rosette_centroids_by_pid" in s02),
    ("02 supports Class 0 and Class 1 candidate boxes & is_cauline_or_stem_contaminant", "s1_cls not in (0, 1)" in s02 and "LM2_LEAF_CLASS_IDS" in s02 and "is_cauline_or_stem_contaminant" in s02),
    ("03 groups by (catalogNumber, plant_individual_id), uses cov.wt & exports morph_tier_source", "plant_individual_id" in s03 and "catalogNumber" in s03 and "cov.wt" in s03 and "morph_tier_source" in s03),
    ("03 matches tier via case-insensitive regex", "tier_?1" in s03.lower() or "tier_?[12]" in s03.lower()),
    ("03 exports harmonic_roughness, is_contour_outlier & leaf_quality_weight", "harmonic_roughness" in s03 and "is_contour_outlier" in s03 and "leaf_quality_weight" in s03),
    ("04 enforces 2-step Dawid-Skene / EB authority update & full-rank canonical axes", "compute_empirical_bayes_authority_weights" in s_auth and "alpha_0" in s_auth and "can3" in s04),
    ("04 exports morphometric_outlier_flag, loocv_accuracy, cda_posterior_margin & can3", "morphometric_outlier_flag" in s04 and "loocv_accuracy" in s04 and "cda_posterior_margin" in s04 and "can3" in s04),
    ("05 uses promoted bronze anchor parity and conditional abaxial artifact gate", "is_promoted_bronze" in s05 and "filter_abaxial_artifacts" in s05),
    ("05 uses foreground patch pooling, OOD scoring & Pointing Game attribution", ("foreground_masked_patch_pooling" in s05 or "compute_foreground_tissue_mask" in s05) and "vision_ood_score" in s05 and "foreground_attention_ratio" in s05),
    ("06 uses 0-15cm depth weighting, VIF pruning, spatial_block_cv_acc & fine_grain_edaphic_cache", "left_join" in s06 and "prune_collinear_predictors" in s06 and "spatial_block_cv_acc" in s06 and "fine_grain_edaphic_cache" in s06),
    ("06 pools Spatial Block-CV RF posteriors and computes edaphic_z_score", "spatial_block_cv_acc" in s06 and "edaphic_z_score" in s06 and "edaphic_best_fit_taxon" in s06),
    ("07 validates contracts, coalesces non-georeferenced vouchers & enforces Tier_1_Pure cryptic gate", "validate_upstream_contracts" in s07 and "taxonomic_treatment_broad" in s07 and "taxonomic_treatment_subtype" in s07 and "is_completed_voucher" in s07 and "Tier_1_Pure" in s07),
    ("07 enforces cross-modal non-contradiction and class-conditional precision", "morph_or_vis_corrupt_reassign" in s07 and "pred_prec" in s07),
]

failed = False
for label, ok in checks:
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {label}")
    if not ok:
        failed = True

if failed:
    sys.exit(1)
EOF
```

---

## Step 2B: Verify Canonical Two-Tier Storage Architecture (`data/` vs. `runs/<RUN_ID>/`) & Figure Suite Contracts

When executing via `./run_pipeline.sh` or `python main.py --run-id <name>` (or `PACKERA_RUN_DIR`), verify adherence to the **Canonical Two-Tier Storage Architecture**:

* **Tier A — Canonical Shared Store (`data/` at repository root):** Houses all invariant Step 01 and Step 02 specimen, environmental, and organ extraction assets, backed relationally by `data/packera_pipeline.sqlite`:
  * **Central Relational DB (`data/packera_pipeline.sqlite`):** Holds Tier A shared tables: `taxon_synonyms`, `vouchers` (PK: `catalogNumber`), `plant_individuals` (PK: `(catalogNumber, plant_individual_id)`), `extracted_leaves` (PK: `(catalogNumber, plant_individual_id, leaf_id)`, storing the 120-vertex `contour_xy_json`), `voucher_reproductive_phenology` (PK: `(catalogNumber, plant_individual_id)`), `dinov2_feature_cache` (PK: `(catalogNumber, plant_individual_id, model_name, sigma)`, caching 768-d embeddings and `foreground_attention_ratio`), and `environmental_cache` + `environmental_rtree` (PK: `lat_lon_key`, spatial R*Tree index).
  * `data/raw_vouchers/` & `data/environmental/` (including `data/environmental/epa_ecoregions_l3/`, `fine_grain_edaphic_cache.csv`, `polaris_30m/`, `solus_100m/`, `dem_3dep_30m/`, `ssurgo_cache/`)
  * `data/masks/tier1_intact/` & `data/masks/tier2_reflected/`
  * `data/contours/` (120-vertex CSV contour fallback files)
  * `data/cropped_patches/{rosettes_dense,basal_leaves_raw,capitula}/`
  * The 5 Step 01–02 CSV manifests in `data/tables/`: `curated_vouchers.csv`, `extracted_leaf_manifest.csv`, `failed_qc_vouchers.csv`, `voucher_reproductive_metrics.csv`, and `voucher_phenological_states.csv`.
* **Tier B — Run-Scoped Analysis Directories (`runs/<RUN_ID>/` with `runs/latest` symlink):**
  * **Run-Scoped SQLite Tables (`data/packera_pipeline.sqlite`):** Partitioned by `run_id`: `pipeline_runs`, `run_active_cohort`, `leaf_efa_harmonics`, `otu_morphometrics`, `vision_xai_audit`, `spatial_rf_predictions`, `curation_triage`, and the pre-joined master view `v_integrative_triage_master`.
  * Automatically seeds Step 01–02 CSV tables (applying `--exclude-taxon`, `--focal-taxon`, or `--limit` cohort filters via SQL queries on `run_active_cohort` and `curation_triage`) and symlinks `data/masks`, `data/contours`, and `data/cropped_patches/` (via directory-level symlinks) directly from canonical root `data/` without referencing previous `runs/` directories.
  * Isolates all Step 03–08 statistical CSV tables (`leaf_efa_harmonics.csv`, `fourier_pcs.csv`, `symmetric_fourier_pcs.csv`, `morphometric_flags.csv`, `label_noise_audit.csv`, `multimodal_conflict_flags.csv`), reports (`final_curation_triage.csv`, `multimodal_spatial_rf_summary.csv`, `edaphic_niche_equivalency.csv`, `micro_edaphic_outcrop_contingency.csv`, `gmm_bayes_factors_summary.csv`, `ecoregion_supervised_cda_summary.csv`, `manova_group_significance.csv`, `manova_pairwise_posthoc.csv`, `discovered_taxon_complexes.csv`, `Packera_Taxonomic_Revision_Summary.md` or `<Focal_Taxon>_Taxonomic_Revision_Summary.md`, `run_metadata.json`), logs (`logs/pipeline.log`), and figures (`outputs/figures/`) inside `runs/<RUN_ID>/`.
* **Publication Figure Suite, 3D Ordinations & Universal Unfilled Marker Standard (`${RUN_DIR}/outputs/figures/`):** Dual vector `.pdf` and 300-DPI `.png` figures with dynamic Okabe-Ito palettes (`get_dynamic_packera_palette`), dynamic hollow marker shapes (`get_dynamic_packera_shapes`), dynamic complex titles (`format_active_complex_title`), italicized *Packera* taxon labels, 3D PCA/CDA ordinations (`PC1 x PC2 x PC3`, `Sym_PC1 x Sym_PC2 x Sym_PC3`, `Can1 x Can2 x Can3` via `plot_3d_scatter_box` in `pub_figures_utils.R`), dynamic North American map bounding boxes, and unfilled (hollow) scatterplot markers (`pch = 1, 2, 5, 0` / `bg = NA` in base R; `shape = 21..25, fill = NA` in `ggplot2`; `facecolors = 'none'` in `matplotlib`):
  * Step 02: `fig02a_segmentation_qc_distributions.{pdf,png}`, `fig02b_silhouette_montage.{pdf,png}`
  * Step 03: `fig03_efa_morphospace_reconstruction.{pdf,png}`, `fig03a_efa_morphospace_pca.{pdf,png}`, `fig03b_harmonic_power_and_asymmetry.{pdf,png}`
  * Step 04: `fig04_gmm_cda_clustering.{pdf,png}`, `fig04a_gmm_bic_and_uncertainty.{pdf,png}`, `fig04b_cda_biplot_and_allometry.{pdf,png}`, `fig04c_supervised_ecoregion_cda.{pdf,png}`, `shape_allometry_regression.{pdf,png}`
  * Step 05: `GradCAM_audit_panel.{png,pdf}`, `fig05b_cleanlab_quality_and_tomentum.{pdf,png}`
  * Step 06: `spatial_rf_niche_importance.{pdf,png}`, `fig06_ecological_niche_phenology.{pdf,png}`, `fig06b_geographic_map_and_hopkins_phenology.{pdf,png}`
  * Step 07: `fig07_multimodal_discordance_summary.{pdf,png}`, `Figure_Integrative_Packera_Revision.{pdf,png}` (or `Figure_Integrative_<Focal_Taxon>_Revision.{pdf,png}`)
  * Step 08: `fig_empirical_taxon_complex_network.{pdf,png}`

---

## Step 3: Execute the Central Test Suite (`run_audit.py`)

Always run the centralized QA test suite to verify 100% pass rate across all unit and integration tests:

```bash
python3 -m unittest discover -s scripts/tests -p "test_*.py" -v
```

* If any test fails or errors, fix the underlying regression before committing.
