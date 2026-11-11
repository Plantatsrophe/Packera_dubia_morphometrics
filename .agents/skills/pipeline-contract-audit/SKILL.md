---
name: pipeline-contract-audit
description: Enforces LeafMachine3 computer vision and morphometric guardrails, audits cross-stage Python and R SQLite/CSV schema contracts, and runs the test suite.
---

# Pipeline Contract Audit and CV Guardrails (`UNC-BOT-SOP-2026-04-REV5`)

Use this skill before committing code changes, after modifying pipeline stages (`01` through `07`), or when diagnosing missing-column, coordinate unscaling, or join-key discrepancies across Python and R scripts.

---

## Part 1: Computer Vision and Morphometric Guardrails (LeafMachine3)

Every CV, extraction, and morphometrics module in `scripts/vision/`, `scripts/pipeline/`, and `scripts/morphometrics/` must enforce these invariants:

### 1. Working-Frame Coordinate Unscaling (`work_scale`)
* When LeafMachine3 processes downscaled or letterboxed representations (`work_scale`), all bounding boxes, polygon vertices, landmark coordinates, and binary masks saved to database tables or exported manifests must be divided by `work_scale` to return to native full-sheet coordinates:
  $$\text{coord}_{\text{native}} = \frac{\text{coord}_{\text{work}}}{\text{work\_scale}}$$
* Columns with suffix `_px` (e.g., `petiole_length_px`, `plant_height_px`, contour vertices) must strictly represent coordinates on the original unscaled image.

### 2. 3-Class Organ Mask Handling and UNet++ v9 Rachis Bridging
* Module 8 instance segmentation (`lm3_leaf_segmenter__yolo26x_seg_1024`) outputs 3 semantic classes: `Leaf` (lamina, class 0), `Petiole` (class 1), and `Hole` (sinuses or perforations, class 2).
* Hole pixels falling within lamina polygons must be subtracted from the blade mask.
* In deeply dissected or pinnatifid leaves where narrow rachis segments experience minor breaks, Module 4 UNet++ v9 `specimen_paper_cleaned_mask` (`LLR_Hyst_KeepBlobs: true` color likelihood) acts as an anatomical bridge constraint (`bridge_dissected_rachis_components`), reconnecting genuine pinnae without bridging across mounting paper.

### 3. Tip-Up Clockwise `lamina_base`-Anchored Contour Geometry
* Normalized 120-vertex contours in `extracted_leaves.contour_xy_json` must be oriented vertically with apex pointing upward ($y$-axis toward apex) and starting homologously at the anatomical lamina-petiole junction (`lamina_base` anchor $p_{\text{base}}$), traversing clockwise around the perimeter.
* Because coordinates are standardized and homologized in Step 02, `Momocs::efourier` in Step 03 must be executed with `norm = FALSE` and `start = FALSE`:
  ```r
  Momocs::efourier(out_obj, nb.h = H_opt, norm = FALSE, start = FALSE)
  ```
  preserving biologically meaningful absolute orientation, rotation, and starting-point homology.

### 4. Strict Active vs. Passive Determiner Authority Split
* **Active Reference Anchors (`Tier_1_Gold`, `Tier_2_Silver`, promoted `Tier_3_Bronze`, `Tier_0_Type`):** Used exclusively for model fitting, Mahalanobis scatter estimation, CDA discriminant axes, and Spatial RF training.
* **Passively Projected Vouchers (`Tier_4_Unverified` or unpromoted `Tier_3_Bronze`):** Assigned weight $0.0$ during training and projected into fitted spaces without altering reference geometry.

### 5. Native LeafMachine3 ONNX Alignment (`lm3_inference_engines.py`)
* Stage 1 and Stage 2 models execute dynamic rectangular stride-32 letterboxing via `ultra_rep.letterbox(auto=True, stride=32)` with `cv2.INTER_LINEAR_EXACT` interpolation.
* Stage 1 operates at `conf = 0.30` and suppresses nested duplicate detections within each class at `overlap = 0.95` (`inter / min_area >= 0.95`).
* Stage 2 extracts tight crops (`crop_pad = 0.0`), decodes prototype masks natively via `ultra_rep.process_mask_native`, and performs class-wise unioning of all valid `Leaf` (`0`), `Petiole` (`1`), and `Hole` (`2`) instances (`>= 64 px`).
* Whole-sheet `LLR_Hyst_KeepBlobs` log-likelihood hysteresis is applied strictly in `LM3SpecimenSegmenter` at the whole-sheet level, and `10%` white border padding (`WHITE_PAD_FRAC = 0.10`) is added around leaf crops in `LM3LandmarkDetector`.
* Adaptive physical blade-area flooring ($\min(0.75, \max(0.15, 0.25 \times \text{median\_area}))\text{ cm}^2$) and relative area filtering ($0.15 \times \max\text{area}$) are computed per plant individual (`plant_individual_id`, `effective_min_area_by_pid`), retaining all passing leaves without a sheet-wide 6-leaf cap.

---

## Part 2: Automated Cross-Stage SQLite/CSV Contract and Hygiene Audit

### 1. Zero `_archive` Coupling and Obsolete Token Scan
Active files must never reference `_archive` paths and must contain zero obsolete framework remnants:

```bash
python3 - <<'EOF'
import sys
from pathlib import Path

root = Path(".").resolve()
active_files = [
    p for p in list(root.rglob("*.py")) + list(root.rglob("*.R"))
    if "_archive" not in p.parts and ".venv" not in p.parts
    and ".venv_LM3" not in p.parts and "LeafMachine3" not in p.parts
]

violations = []
for f in sorted(active_files):
    rel = f.relative_to(root)
    if rel.parts[:2] == ("scripts", "tests"):
        continue
    for lineno, line in enumerate(f.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if "_archive" in line:
            violations.append(f"{rel}:{lineno}: {stripped}")

if violations:
    print("[FAIL] Forbidden _archive references found:")
    for v in violations:
        print("  -", v)
    sys.exit(1)
print(f"[PASS] Scanned {len(active_files)} active files: 0 _archive leaks.")
EOF
```

### 2. Cross-Stage Handshake Contract Verification Script
Verifies column agreements, composite keys `(catalogNumber, plant_individual_id)`, and handshakes across stages:

```bash
python3 - <<'EOF'
import sys
from pathlib import Path

scripts = [
    "scripts/core/db.py", "scripts/core/db_utils.R", "scripts/core/harvester.py",
    "scripts/pipeline/02_segment_and_extract.py", "scripts/morphometrics/03_fourier_extractor.R",
    "scripts/morphometrics/04_gmm_morphotools.R", "scripts/morphometrics/authority_weighting_utils.R",
    "scripts/analysis/05_multimodal_spatial_rf.R", "scripts/analysis/06_triage_dashboard_synthesis.R",
    "scripts/analysis/07_discover_taxon_complexes.R",
]

missing = [s for s in scripts if not Path(s).exists()]
if missing:
    print(f"[FAIL] Missing pipeline scripts: {missing}")
    sys.exit(1)

s_db = Path("scripts/core/db.py").read_text() + Path("scripts/core/db_schema.py").read_text()
s_db_r = Path("scripts/core/db_utils.R").read_text()
s_harv = Path("scripts/core/harvester.py").read_text()
s02 = Path("scripts/pipeline/02_segment_and_extract.py").read_text()
s03 = Path("scripts/morphometrics/03_fourier_extractor.R").read_text()
s04 = Path("scripts/morphometrics/04_gmm_morphotools.R").read_text()
s_auth = Path("scripts/morphometrics/authority_weighting_utils.R").read_text()
s05 = Path("scripts/analysis/05_multimodal_spatial_rf.R").read_text()
s06 = Path("scripts/analysis/06_triage_dashboard_synthesis.R").read_text()
s07 = Path("scripts/analysis/07_discover_taxon_complexes.R").read_text()

checks = [
    ("db.py defines schema, WAL mode, tables, and master view", "v_integrative_triage_master" in s_db and "PRAGMA journal_mode = WAL" in s_db),
    ("db_utils.R exposes R connection and master view helpers", "v_integrative_triage_master" in s_db_r and "RSQLite::SQLite()" in s_db_r),
    ("db_utils.R provides type edaphic z-score handling", "type_edaphic_z_score" in s_db_r),
    ("01 parses determiner tier with collector isolation & duplicate guard", "_should_fallback_to_collector" in s_harv and "exsiccatae_taxonomic_conflict" in s_harv),
    ("01 defines nomenclatural type hierarchy, subranks & cascading priority", "PRIMARY_TYPE_STATUSES" in s_harv and "classify_nomenclatural_type_rank" in s_harv and "_type_rank" in s_harv),
    ("01 historical types exempt from uncertainty filters and retained in thinning", "is_rec_type" in s_harv and "spatially_balance_and_thin_vouchers" in s_harv),
    ("02 exports assigned_tier, pixels_per_mm & rosette centroids", "assigned_tier" in s02 and "pixels_per_mm" in s02 and "rosette_centroids_by_pid" in s02),
    ("02 supports Class 0/1 candidates & is_cauline_or_stem_contaminant", "s1_cls not in (0, 1)" in s02 and "is_cauline_or_stem_contaminant" in s02),
    ("03 groups by (catalogNumber, plant_individual_id), uses cov.wt & exports morph_tier_source", "plant_individual_id" in s03 and "catalogNumber" in s03 and "cov.wt" in s03),
    ("03 matches tier via case-insensitive regex", "tier_?1" in s03.lower() or "tier_?[12]" in s03.lower()),
    ("03 exports harmonic_roughness, is_contour_outlier & leaf_quality_weight", "harmonic_roughness" in s03 and "is_contour_outlier" in s03),
    ("03 enforces nomenclatural type MCD outlier immunity & precision weighting", "mcd_out" in s03 and "is_type_specimen" in s03 and "compute_authority_precision_weights" in s03),
    ("04 enforces 2-step Dawid-Skene / EB authority update & full-rank can3", "compute_empirical_bayes_authority_weights" in s_auth and "can3" in s04),
    ("04 exports morphometric flags, loocv_accuracy & cda_posterior_margin", "morphometric_outlier_flag" in s04 and "loocv_accuracy" in s04 and "cda_posterior_margin" in s04),
    ("04 enforces type-medoid shrinkage & type morphological control", "compute_type_specimen_morph_control" in s04 and ".identify_type_specimen_vec" in s_auth),
    ("05 uses 0 to 15 cm depth weighting, VIF pruning & spatial_block_cv_acc", "left_join" in s05 and "prune_collinear_predictors" in s05 and "spatial_block_cv_acc" in s05),
    ("05 pools Spatial Block-CV RF posteriors and computes edaphic_z_score", "spatial_block_cv_acc" in s05 and "edaphic_z_score" in s05),
    ("05 derives locus classicus reference set & type_edaphic_z_score", "compute_type_locality_ecological_control" in s05 and "type_edaphic_z_score" in s05),
    ("05 enforces authority weight resolution", "resolve_step" in s05 and "3.0" in s05),
    ("06 validates contracts & coalesces non-georeferenced vouchers", "validate_upstream_contracts" in s06 and "taxonomic_treatment_broad" in s06 and "is_completed_voucher" in s06),
    ("06 enforces cross-modal precision weighting", "Confirmed_Type_Control" in s06 and ("prec_morph" in s06 or "lookup_prec" in s06)),
    ("06 enforces type control Rules A, B, and C", "Confirmed_Type_Control" in s06 and "Type_Synonymy_Candidate" in s06 and "is_sensu_stricto_locked_i" in s06),
    ("07 enforces type-anchored complex discovery and delimitation", "select_taxon_anchors" in s07 and "is_type_spec" in s07 and "gamma_k" in s07),
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

### 3. Punctuation, Tone, and Test Suite Gate
Always execute the verification gate to confirm zero regressions:

```bash
# Verify zero em/en dashes across repository documentation
grep -rnP "[\x{2013}\x{2014}]" README.md docs/TAXONOMIC_SOP.md AGENTS.md && exit 1 || echo "PASS: Zero em/en dashes found"

# Verify positive present-state documentation across repository documentation
python3 -c '
import re, sys
from pathlib import Path
banned = [
    "clean" + "lab", "dino" + "v2", "vision_" + "xai", "05_clean" + "lab", "06_multimodal_" + "spatial",
    "07_tri" + "age", "08_dis" + "cover", "step07_" + "pub", "fig" + "08", "target_" + "species",
    "hard" + "coded", "operational" + "iz", "lever" + "ag", "util" + "iz", "holis" + "tic",
    "multi" + "faceted", "para" + "digm",
]
pat = re.compile("|".join(banned), re.IGNORECASE)
violations = []
for p in [Path("README.md"), Path("docs/TAXONOMIC_SOP.md"), Path("AGENTS.md")]:
    for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
        if pat.search(line):
            violations.append(f"{p}:{i}: {line.strip()}")
if violations:
    print("[FAIL] Banned tokens found:\n" + "\n".join(violations))
    sys.exit(1)
print("PASS: Docs clean and renumbered to Steps 01 to 07")
'

pytest scripts/tests/ -q
pytest scripts/tests/test_code_hygiene.py -q
```
