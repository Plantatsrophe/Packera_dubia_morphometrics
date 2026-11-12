# Contributing and Developer Guide (`Packera_dubia_morphometrics`)

**Standard Operating Procedure:** [`docs/TAXONOMIC_SOP.md`](docs/TAXONOMIC_SOP.md) (`UNC-BOT-SOP-2026-04-REV5`)
**Architectural Contract:** [`AGENTS.md`](AGENTS.md) | **Schema Reference:** [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md)

I welcome bug fixes, botanical authority-table additions, regional environmental adapters, and statistical improvements that preserve the mathematical and curatorial rigor of this pipeline. Because Python harvesting and ONNX computer vision stages hand off directly to R morphometric and spatial niche models through a shared SQLite Write-Ahead Logging (`WAL`) store, every pull request must satisfy the architectural boundaries and automated audit checks below before merge.

## 1. Dual Virtual Environment and `LeafMachine3` Submodule Boundary

To keep compiled geospatial C++ libraries (`GDAL`, `GEOS`, `PROJ`) from conflicting with pinned ONNX Runtime and `CUDA` wheels, the repository maintains two strictly separated Python environments:

* **`.venv` (Primary Orchestration & Geospatial Python):** Created from `requirements.txt`. Runs `main.py`, `run_pipeline.sh` preflight helpers, Step 01 GBIF/iDigBio harvesting, EPA ecoregion and USDA county range lookups, `download_environmental_layers.py`, SQLite/DuckDB migrations, and the Python test runner.
* **`LeafMachine3/.venv` (`.venv_LM3` Vision Isolation):** Provisioned automatically via `uv` by `bash setup_leafmachine3.sh`. Executes Step 02 (`02_segment_and_extract.py`) ONNX inference, native `LeafMachine3` visual rendering, and $32 \times 360$ ECT calculation, writing to SQLite strictly through Python's standard-library `sqlite3` module.

Never modify, patch, or commit changes inside `LeafMachine3/`. All `Packera`-specific computer vision adaptations belong in our read-only adapter layer under `scripts/vision/` and `scripts/pipeline/02_segment_and_extract.py`, and all working-frame coordinates (`coord_work`) must divide by `work_scale` before entering SQLite or CSV outputs.

## 2. Module Ownership and R Script Line Budgets

Keep functions cohesive within their designated directories (`scripts/core/`, `scripts/data_prep/`, `scripts/vision/`, `scripts/pipeline/`, `scripts/morphometrics/`, `scripts/figures/`, and `scripts/analysis/`), and never import or `source()` anything from `_archive` folders.

Three core R analytical scripts enforce hard line-count ceilings verified automatically by the test suite so statistical logic stays modular and readable:

| Script Path | Maximum Line Budget | Companion Helper / Figure Modules |
| :--- | :---: | :--- |
| `scripts/morphometrics/03_fourier_extractor.R` | `<= 800` lines | `dissected_leaf_utils.R`, `fourier_parallel_utils.R`, `scripts/figures/step03_figures.R` |
| `scripts/morphometrics/04_gmm_morphotools.R` | `<= 650` lines | `authority_weighting_utils.R`, `plasticity_mitigation_utils.R`, `pub_ecoregion_manova_utils.R`, `scripts/figures/step04_figures.R` |
| `scripts/analysis/taxon_complex_discovery_utils.R` | `<= 800` lines | `06_triage_dashboard_synthesis.R`, `scripts/figures/step06_figures.R` |

When adding new morphometric diagnostics or plots, place shared plotting logic in `scripts/figures/theme_packera_pub.R` or `scripts/figures/step0X_figures.R` rather than expanding the core modeling entrypoints past their line budgets.

## 3. SQLite Schema (`20261121`), Composite Keys, and Demo Isolation

1. **Synchronized DDL Across Python and R:** Whenever you add or rename a table column or modify `v_integrative_triage_master`, update `scripts/core/db_schema.py`, `scripts/core/db_utils.R`, and [`docs/DATA_DICTIONARY.md`](docs/DATA_DICTIONARY.md) in the same commit, and increment `PACKERA_SCHEMA_USER_VERSION` if breaking migrations are required.
2. **Composite Individual Keys:** Group and join all specimen-level data on `(catalogNumber, plant_individual_id)` (plus `run_id` in Tier B tables). Never deduplicate by `catalogNumber` alone in Steps 02 through 06, as multi-plant herbarium sheets carry distinct DBSCAN individuals (`plant_individual_id = 1, 2, ...`).
3. **Atomic Writes and Honest Missingness:** Write CSV, HTML, Markdown, and GeoTIFF deliverables to process-unique `.tmp` or `.part` files first and `rename()` atomically upon completion. Outside explicit `--demo` runs (`PACKERA_DEMO_MODE=1`, which isolates writes inside `demo_ephemeral.sqlite`), never synthesize mock coordinates, random embeddings, or fabricated raster values when external layers are absent; return `NA` / `NULL` cleanly.
4. **RNG State Hygiene in R:** Any stochastic call inside helper routines (`igraph::cluster_leiden`, `kmeans`, bootstrap/permutation loops) must preserve and restore the caller's `.Random.seed` via `on.exit` (see `.with_preserved_rng_seed()` in `taxon_complex_discovery_utils.R`).

## 4. 300-DPI PNG-Only Figure Protocol

All pipeline visualizations export strictly as **300-DPI PNG** files (`0` `.pdf` outputs) through `save_packera_png()` in `scripts/figures/theme_packera_pub.R` (or Matplotlib `dpi=300` in `scripts/vision/step02_figures.py`):

* **Dynamic Taxon and Variable Scaling:** Figures must handle $N = 1$ intraspecific runs (falling back to unsupervised GMM morphotype clusters), $G = 2$ pairwise comparisons (1D `Can1` density + raincloud rather than degenerate 2D axes), and $60+$ species or infraspecific ranks (`subsp.`, `var.`, `f.`) using `build_packera_taxon_aesthetic_map()` (12 colorblind-safe hues $\times$ 5 unfilled hollow shapes) and dynamic legend canvas heights.
* **Botanical & Cartographic Standards:** Format all taxon labels through `format_packera_taxon_markdown()` / ICN expression labellers so genus and species/infraspecific epithets appear in italics while rank markers (`subsp.`, `var.`, `f.`, `s.s.`, `s.l.`, `sp. nov.`) remain Roman. Project all geographic maps to equal-area `EPSG:5070` (`NAD83 / Conus Albers`).

## 5. Documentation and Prose Style (`Writing DNA v5.2`)

All Markdown documentation (`README.md`, `AGENTS.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, and `docs/*.md`) follows the prose rules in [`.agents/skills/write-in-my-voice/SKILL.md`](.agents/skills/write-in-my-voice/SKILL.md):

* Write in clear, active first-person or direct imperative sentences with varied paragraph lengths and explicit causal or contrastive connectors (`, as ...`, `, whereas ...`, `, so ...`).
* Enforce **zero em dashes**, **zero en dashes**, and **zero prose double-hyphens** (write numerical ranges as `"0 to 15 cm"` or `"Steps 01 to 06"`). CLI flags inside code spans (`--focal-taxon`) are fine.
* Cite academic literature using clean `doi: 10.xxxx/...` identifiers rather than raw publisher URLs, and describe active pipeline architecture in positive present-state terms.

## 6. Running the Pre-Commit Verification Suite (`222` Tests)

Before opening a pull request or generating a patch, run the full verification gate from the repository root:

```bash
# 1a. Run the complete 222-test Python, SQLite handshake, warm-R statistical, figure, and hygiene suite
.venv/bin/python scripts/tests/run_audit.py

# 1b. Optional rapid Python-only pre-check during active iteration (~25s)
.venv/bin/python scripts/tests/run_audit.py --fast

# 2. Verify zero em/en dashes across all Markdown and citation files
grep -rnP "[\x{2013}\x{2014}]" README.md AGENTS.md CONTRIBUTING.md CHANGELOG.md CITATION.cff docs/*.md && exit 1 || echo "PASS: Zero em/en dashes"

# 3. Verify R script line-count budgets
wc -l scripts/morphometrics/03_fourier_extractor.R scripts/morphometrics/04_gmm_morphotools.R scripts/analysis/taxon_complex_discovery_utils.R

# 4. Verify git whitespace cleanliness
git diff --check
```
