---
name: pipeline-contract-audit
description: Audits cross-stage Python and R CSV schema contracts, checks for forbidden _archive or synthetic fallback leaks, and runs the test suite. Use before committing code changes or when debugging pipeline stage handshakes.
---

# Pipeline Contract & Schema Audit (`UNC-BOT-SOP-2026-04-REV4`)

Use this skill before committing code changes, after modifying any pipeline stage (`01` through `07`), or when diagnosing missing-column / join-key errors across Python and R scripts.

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

Verify that all six core production stages (`02` through `07`) agree on primary keys, column names, join semantics, and dynamic reporting:

### Canonical Stage-to-Stage Handshake Matrix

| Stage Script | Primary Inputs | Primary Outputs | Required Keys & Shared Columns | Critical Handshake Invariants |
| :--- | :--- | :--- | :--- | :--- |
| **`scripts/pipeline/02_segment_and_extract.py`** | `data/tables/curated_vouchers.csv` | `data/tables/extracted_leaf_manifest.csv`, `data/tables/voucher_reproductive_metrics.csv`, `data/contours/*.csv` | `catalogNumber`, `plant_individual_id`, `leaf_id`, `assigned_tier`, `pixels_per_mm`, `is_calibrated_mm`, `contour_path`, `mask_path` | Must preserve multi-plant `plant_individual_id` clusters via DBSCAN; `assigned_tier` exports `"Tier_1_Direct"` / `"Tier_2_Reflected"`. |
| **`scripts/morphometrics/03_fourier_extractor.R`** | `data/tables/extracted_leaf_manifest.csv`, `data/contours/*.csv`, `data/tables/curated_vouchers.csv` | `data/tables/leaf_efa_harmonics.csv`, `data/tables/leaf_efa_harmonics_individual.csv` | `(catalogNumber, plant_individual_id)`, `assigned_tier`, `reconstruction_tier`, `pixels_per_mm`, `is_calibrated_mm`, `centroid_size`, `log_centsize`, `PC1`–`PC5`, `A1`–`D12` | Matches tier case-insensitively via `grepl("^tier_?[12]", assigned_tier, ignore.case = TRUE)`; aggregates specimen medians by `(catalogNumber, plant_individual_id)` (never collapsing to `plant_individual_id = 0`). |
| **`scripts/morphometrics/04_gmm_morphotools.R`** | `data/tables/leaf_efa_harmonics.csv`, `data/tables/voucher_reproductive_metrics.csv` | `data/tables/morphometric_flags.csv` (and `morphometrics_misidentification_flags.csv`), `outputs/reports/gmm_bayes_factors_summary.csv` | `(catalogNumber, plant_individual_id)`, `determiner_tier`, `species_standardized`, `cda_predicted_taxon`, `cda_posterior_prob`, `can1`, `can2`, `gmm_cluster`, `gmm_uncertainty`, `morphometric_outlier_flag` | Aggregates by `(catalogNumber, plant_individual_id)`; exports `morphometric_outlier_flag` (aliased to `misidentification_flag` for downstream consumers). |
| **`scripts/analysis/05_cleanlab_vision_xai.py`** | `data/tables/curated_vouchers.csv`, `data/cropped_patches/` | `data/tables/label_noise_audit.csv` | `catalogNumber`, `plant_individual_id`, `predicted_label` (`vision_predicted_label`), `confidence_predicted_class`, `c_error`, `is_label_corrupted`, `surface_orientation` | Must neutralize mounting paper (`RGB(128,128,128)`) and control for `surface_orientation` (`"abaxial"` vs. `"adaxial"`). Never synthesize random embeddings unless `--demo` is passed. |
| **`scripts/analysis/06_multimodal_spatial_rf.R`** | `data/tables/curated_vouchers.csv`, `data/tables/morphometric_flags.csv`, `data/tables/label_noise_audit.csv`, `data/tables/voucher_phenological_states.csv` (optional) | `data/tables/multimodal_conflict_flags.csv`, `outputs/reports/multimodal_spatial_rf_summary.csv`, `data/tables/ssurgo_edaphic_validation.csv` | `catalogNumber`, `plant_individual_id`, `misidentification_flag` (resolved from `morphometric_outlier_flag`), `c_error`, `edaphic_best_fit_taxon`, `multimodal_concordance`, `triage_category` | Must use **`dplyr::left_join`** (never `inner_join`) when merging optional phenology tables (`voucher_phenological_states.csv`) so vegetative-only vouchers are never dropped. Never cache synthetic SSURGO/SoilGrids values to disk. |
| **`scripts/analysis/07_triage_dashboard_synthesis.R`** | `data/tables/curated_vouchers.csv`, `data/tables/morphometric_flags.csv`, `data/tables/label_noise_audit.csv`, `data/tables/multimodal_conflict_flags.csv`, `outputs/reports/gmm_bayes_factors_summary.csv` | `data/tables/triage_queue.csv`, `outputs/figures/Figure_Integrative_Packera_dubia_Revision.pdf`, `outputs/reports/Packera_dubia_Taxonomic_Revision_Summary.md` | `(catalogNumber, plant_individual_id)`, `morphometric_outlier_flag` / `misidentification_flag`, `taxonomic_status_call`, `recommended_determination`, `synthesis_triage_priority` | Must resolve `morphometric_outlier_flag` <-> `misidentification_flag` seamlessly, dynamically compute CDA axis variance labels and GMM $2\Delta\text{BIC}$ values from input tables/summaries (never hardcode static percentages or $2\Delta\text{BIC}$ constants), and write outputs atomically. |

### Automated Contract Verification Script

Run this Python audit to verify column-name agreement and join-key rules across `02`–`07`:

```bash
python3 - <<'EOF'
import sys
from pathlib import Path

scripts_to_check = [
    "scripts/pipeline/02_segment_and_extract.py",
    "scripts/morphometrics/03_fourier_extractor.R",
    "scripts/morphometrics/04_gmm_morphotools.R",
    "scripts/analysis/05_cleanlab_vision_xai.py",
    "scripts/analysis/06_multimodal_spatial_rf.R",
    "scripts/analysis/07_triage_dashboard_synthesis.R",
]

missing = [s for s in scripts_to_check if not Path(s).exists()]
if missing:
    print(f"[FAIL] Missing pipeline scripts: {missing}")
    sys.exit(1)

s02 = Path("scripts/pipeline/02_segment_and_extract.py").read_text()
s03 = Path("scripts/morphometrics/03_fourier_extractor.R").read_text()
s04 = Path("scripts/morphometrics/04_gmm_morphotools.R").read_text()
s05 = Path("scripts/analysis/05_cleanlab_vision_xai.py").read_text()
s06 = Path("scripts/analysis/06_multimodal_spatial_rf.R").read_text()
s07 = Path("scripts/analysis/07_triage_dashboard_synthesis.R").read_text()

checks = [
    ("02 exports assigned_tier & pixels_per_mm", "assigned_tier" in s02 and "pixels_per_mm" in s02),
    ("03 groups by (catalogNumber, plant_individual_id)", "plant_individual_id" in s03 and "catalogNumber" in s03),
    ("03 matches tier via case-insensitive regex", "tier_?1" in s03.lower() or "tier_?[12]" in s03.lower()),
    ("04 exports morphometric_outlier_flag", "morphometric_outlier_flag" in s04),
    ("05 exports c_error & is_label_corrupted", "c_error" in s05 and "is_label_corrupted" in s05),
    ("06 uses left_join", "left_join" in s06),
    ("07 ingests CDA/GMM & Vision flags", "cda_predicted_taxon" in s07 and "c_error" in s07),
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

## Step 3: Execute the Central Test Suite (`run_audit.py`)

Always run the centralized QA test suite to verify 100% pass rate across all unit and integration tests:

```bash
python scripts/tests/run_audit.py
```

* Inspect `/tmp/lm2_audit_results.txt` and stdout.
* If any test fails or errors, fix the underlying regression before committing.
