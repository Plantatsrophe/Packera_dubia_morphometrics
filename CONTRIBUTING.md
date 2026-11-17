# Contributing and Developer Guide (`Plantmatics`)

**Standard Operating Procedure:** [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md) (`UNC-BOT-SOP-2026-04-REV5`)
**Architectural Contract:** [`AGENTS.md`](AGENTS.md) | **Schema Reference:** [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md)

I welcome bug fixes, botanical authority-table additions, regional environmental adapters, and statistical improvements that preserve the mathematical and curatorial rigor of Plantmatics. Because Python harvesting and ONNX computer vision stages hand off directly to R morphometric and spatial niche models through a shared SQLite Write-Ahead Logging (`WAL`) store (`data/plantmatics.sqlite`), every pull request must satisfy the architectural boundaries and automated audit checks below before merge.

## 1. Cross-Platform `pixi.toml` Workspace, Dual Virtual Environments, and `LeafMachine3` Boundary

To keep compiled geospatial C++ libraries (`GDAL`, `GEOS`, `PROJ`) from conflicting with pinned ONNX Runtime and `CUDA` wheels while supporting native Linux, macOS, and Windows execution, I provide a `conda-forge`-only [`pixi.toml`](pixi.toml) workspace (`linux-64`, `osx-arm64`, `osx-64`, and `win-64` providing `python-isal`, `Rcpp`, `MSCquartets`, and target features for `onnxruntime-directml` and `onnxruntime-gpu`) alongside the two separated Python environments:

* **`.venv` (Primary Orchestration & Geospatial Python):** Created from `requirements.txt` (or `pixi.toml` default environment), resolved across POSIX `.venv/bin/python` and Windows `.venv/Scripts/python.exe` via `resolve_python_executable()` and `resolve_r_python_bin()`. Runs `main.py`, `run_pipeline.sh` preflight helpers, Step 01 GBIF/iDigBio harvesting, EPA ecoregion and USDA county range lookups, `download_environmental_layers.py`, SQLite/DuckDB migrations, and the Python test runner.
* **`LeafMachine3/.venv` (`.venv_LM3` Vision Isolation):** Provisioned automatically via `uv` by `bash setup_leafmachine3.sh`. Executes Step 02 (`02_segment_and_extract.py`) ONNX inference, native `LeafMachine3` visual rendering, and $32 \times 360$ ECT calculation, writing to SQLite strictly through Python's standard-library `sqlite3` module.
* **Windows Parallelism & Run Pointer Fallbacks:** On non-Unix hosts, `packera_mclapply()` (`plantmatics_par_lapply`) inside `scripts/core/db_utils.R` spawns a cached `parallel::makePSOCKcluster()` worker pool instead of dropping to `1` thread, and `initialize_run_directory()` writes an atomic `runs/latest_run.json` pointer whenever filesystem symlink creation raises `OSError`.

Never modify, patch, or commit changes inside `LeafMachine3/`. All taxon-agnostic computer vision adaptations belong in our read-only adapter layer under `scripts/vision/` and `scripts/pipeline/02_segment_and_extract.py`, and all working-frame coordinates (`coord_work`) must divide by `work_scale` before entering SQLite or CSV outputs.

## 2. Module Ownership and Script Line Budgets

Keep functions cohesive within their designated directories (`scripts/core/`, `scripts/data_prep/`, `scripts/vision/`, `scripts/pipeline/`, `scripts/morphometrics/`, `scripts/figures/`, and `scripts/analysis/`), and never import or `source()` anything from `_archive` folders.

Core R analytical entrypoints and public genomics modules enforce hard line-count ceilings verified automatically by the test suite so statistical and ingestion logic stays modular:

| Script Path | Maximum Line Budget | Companion Helper / Figure Modules |
| :--- | :---: | :--- |
| `scripts/core/progress_telemetry.py` | `< 350` lines | `scripts/core/logger.py`, `scripts/core/db_utils.R`, `scripts/tests/test_v15_taxon_and_telemetry.py` (`< 380`) |
| `scripts/morphometrics/03_fourier_extractor.R` | `< 500` lines | `step03_config_scale_utils.R` (`< 240`), `dissected_leaf_utils.R`, `fourier_parallel_utils.R` (`< 500`), `scripts/figures/step03_figures.R` (`< 385`) |
| `scripts/morphometrics/04_gmm_morphotools.R` | `< 500` lines | `step04_allometry_merger_utils.R` (`< 380`), `authority_weighting_utils.R`, `plasticity_mitigation_utils.R` (`< 470`), `pub_ecoregion_manova_utils.R`, `scripts/figures/step04_figures.R` (`< 500`) |
| `scripts/analysis/07_triage_dashboard_synthesis.R` | `< 500` lines | `step07_ancova_merge_utils.R` (`< 490`), `step07_hopkins_consensus_utils.R` (`< 480`), `step07_biogeo_hybrid_utils.R` (`< 490`), `step07_cluster_cohesion_utils.R` (`< 460`), `step07_decision_stages_utils.R` (`< 490`), `step07_postprocess_guards_utils.R` (`< 460`), `step07_reports_html_utils.R` (`< 490`), `step07_cytology_integration.R` (`< 380`), `step07_genomic_propagation.R` (`< 400`), `step07_deme_propagation_utils.R` (`< 260`), `step07_organellar_triage_utils.R` (`< 260`), `genomic_triage_utils.R` (`< 500`) |
| `scripts/analysis/taxon_complex_discovery_utils.R` | `< 480` lines | `step07_complex_affinity_utils.R` (`< 450`), `scripts/figures/step07_figures.R` (`< 420`), `scripts/figures/step07b_network_figures.R` (`< 320`) |
| `scripts/data_prep/download_public_packera_sequences.py` | `< 500` lines | `06_ingest_genetic_data.py`, `skim_offtarget_organellar.py` |
| `scripts/data_prep/skim_offtarget_organellar.py` | `< 500` lines | `06_ingest_genetic_data.py`, `download_public_packera_sequences.py` |
| `scripts/analysis/step06_public_backbone.R` | `< 500` lines | `06_genomic_phylogenomics.R`, `step06_reticulation_coancestry.R` |

When adding new morphometric diagnostics or plots, place shared plotting logic in `scripts/figures/theme_packera_pub.R` or `scripts/figures/step0X_figures.R` rather than expanding the core modeling entrypoints past their line budgets.

## 3. SQLite Schema (`20261128`), Cross-Language `TaxonContext`, Composite Keys, and Telemetry

1. **Synchronized DDL Across Python and R:** Whenever you add or rename a table column or modify `v_integrative_triage_master`, update `scripts/core/db_schema.py`, `scripts/core/db_utils.R`, and [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) in the same commit, and increment `PLANTMATICS_SCHEMA_USER_VERSION` (preserving `PACKERA_SCHEMA_USER_VERSION` alias compatibility) if breaking migrations are required.
2. **Composite Individual Keys:** Group and join all specimen-level data on `(catalogNumber, plant_individual_id)` (plus `run_id` in Tier B tables). Never deduplicate by `catalogNumber` alone in Steps 02 through 07, as multi-plant herbarium sheets carry distinct DBSCAN individuals (`plant_individual_id = 1, 2, ...`).
3. **Atomic Writes and Honest Missingness:** Write CSV, HTML, Markdown, and GeoTIFF deliverables to process-unique `.tmp` or `.part` files first and `rename()` atomically upon completion. Outside explicit `--demo` runs (`PLANTMATICS_DEMO_MODE=1`, with `PACKERA_DEMO_MODE=1` compatibility, which isolates writes inside `demo_ephemeral.sqlite`), never synthesize mock coordinates, random embeddings, or fabricated raster values when external layers are absent; return `NA` / `NULL` cleanly.
4. **RNG State Hygiene in R:** Any stochastic call inside helper routines (`igraph::cluster_leiden`, `kmeans`, bootstrap/permutation loops) must preserve and restore the caller's `.Random.seed` via `on.exit` (see `.with_preserved_rng_seed()` in `taxon_complex_discovery_utils.R`).
5. **Cross-Language `TaxonContext` & Cache-Subtracted ETA Telemetry:** Resolve the active genus via `PLANTMATICS_TARGET_GENUS` (`resolve_target_genus()` in Python and `resolve_r_target_genus()` in R) rather than hardcoding *Packera* regexes or fallback epithets, and persist `taxon_synonyms.canonical_taxon` alongside `usda_county_ranges.canonical_taxon` and `UNASSIGNED_TAXON_SENTINEL = "Unassigned_Taxon"`. Both `emit_progress()` and `emit_r_progress()` append single-line JSONL records so `StepETATracker` subtracts Tier A cache hits before computing EWMA and SQLite `step_telemetry_benchmarks` ETAs.

## 4. 300-DPI PNG-Only Figure Protocol

All pipeline visualizations export strictly as **300-DPI PNG** files (`0` `.pdf` outputs) through `save_packera_png()` in `scripts/figures/theme_packera_pub.R` (or Matplotlib `dpi=300` in `scripts/vision/step02_figures.py`):

* **Dynamic Taxon and Variable Scaling:** Figures must handle $N = 1$ intraspecific runs (falling back to unsupervised GMM morphotype clusters), $G = 2$ pairwise comparisons (1D `Can1` density + raincloud rather than degenerate 2D axes), and $60+$ species or infraspecific ranks (`subsp.`, `var.`, `f.`) using `build_packera_taxon_aesthetic_map()` (12 colorblind-safe hues $\times$ 5 unfilled hollow shapes) and dynamic legend canvas heights.
* **Botanical & Cartographic Standards:** Format all taxon labels through `format_packera_taxon_markdown()` / ICN expression labellers so genus and species/infraspecific epithets appear in italics while rank markers (`subsp.`, `var.`, `f.`, `s.s.`, `s.l.`, `sp. nov.`) remain Roman. Project all geographic maps to equal-area `EPSG:5070` (`NAD83 / Conus Albers`).

## 5. Documentation and Prose Style (`Writing DNA v5.2`)

All Markdown documentation (`README.md`, `AGENTS.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, and `docs/*.md`) follows the prose rules in [`.agents/skills/write-in-my-voice/SKILL.md`](.agents/skills/write-in-my-voice/SKILL.md):

* Write in clear, active first-person or direct imperative sentences with varied paragraph lengths and explicit causal or contrastive connectors (`, as ...`, `, whereas ...`, `, so ...`).
* Enforce **zero em dashes**, **zero en dashes**, and **zero prose double-hyphens** (write numerical ranges as `"0 to 15 cm"` or `"Steps 01 to 07"`). CLI flags inside code spans (`--focal-taxon`) are fine.
* Cite academic literature using clean `doi: 10.xxxx/...` identifiers rather than raw publisher URLs, and describe active pipeline architecture in positive present-state terms.

## 6. Running the Pre-Commit Verification Suite

Before opening a pull request or generating a patch, run the full verification gate from the repository root:

```bash
# 1a. Run the complete Python, SQLite handshake, warm-R statistical, figure, and hygiene suite
.venv/bin/python scripts/tests/run_audit.py

# 1b. Optional rapid Python-only pre-check during active iteration (~25s)
.venv/bin/python scripts/tests/run_audit.py --fast

# 2. Verify zero em/en dashes across all Markdown and citation files
grep -rnP "[\x{2013}\x{2014}]" README.md AGENTS.md CONTRIBUTING.md CHANGELOG.md CITATION.cff docs/*.md && exit 1 || echo "PASS: Zero em/en dashes"

# 3. Verify R script line-count budgets
wc -l scripts/morphometrics/03_fourier_extractor.R scripts/morphometrics/04_gmm_morphotools.R scripts/analysis/07_triage_dashboard_synthesis.R scripts/analysis/taxon_complex_discovery_utils.R

# 4. Verify git whitespace cleanliness
git diff --check
```
