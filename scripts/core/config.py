from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
from datetime import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"
RUNS_ROOT = PROJECT_ROOT / "runs"

# Path fields in PathsConfig that are scoped to the active run directory
# (shared read-only caches like raw_vouchers_dir, quarantine_dir, models_dir, and model_weights remain at PROJECT_ROOT)
RUN_SCOPED_PATH_FIELDS: Set[str] = {
    "curated_vouchers_csv",
    "failed_qc_vouchers_csv",
    "extraction_manifest_csv",
    "contours_dir",
    "masks_dir",
    "leaf_efa_harmonics_csv",
    "morphometric_flags_csv",
    "cda_biplot_pdf",
    "revision_report_md",
    "gmm_report_csv",
    "summary_log",
    "tables_dir",
    "figures_dir",
}


def generate_default_run_id() -> str:
    """Generate a timestamped run identifier: run_YYYYMMDD_HHMMSS."""
    return f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"


def get_active_run_dir(
    run_id: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> Optional[Path]:
    """
    Return the resolved run directory if `run_id` is supplied or `PACKERA_RUN_DIR`
    is set in the environment; otherwise return None (default legacy/test mode).
    """
    if run_id is not None and str(run_id).strip():
        clean_id = str(run_id).strip()
        p = Path(clean_id)
        if p.is_absolute():
            return p.resolve()
        if p.parts and p.parts[0] == "runs":
            return (workspace_root / p).resolve()
        return (workspace_root / "runs" / clean_id).resolve()

    env_run_dir = os.environ.get("PACKERA_RUN_DIR", "").strip()
    if env_run_dir:
        p = Path(env_run_dir)
        return p.resolve() if p.is_absolute() else (workspace_root / p).resolve()

    return None


def initialize_run_directory(
    run_id: Optional[str] = None,
    workspace_root: Path = PROJECT_ROOT,
    update_latest_symlink: bool = True,
    export_env: bool = True,
) -> Path:
    """
    Create the standardized self-contained run directory structure under `runs/<run_id>/`:
      - data/masks
      - data/contours
      - data/cropped_patches/capitula
      - data/tables
      - outputs/{figures,reports}
      - logs
    Also maintains the `runs/latest` symlink pointing to the active run directory.
    """
    run_dir = get_active_run_dir(run_id=run_id, workspace_root=workspace_root)
    if run_dir is None:
        eff_id = os.environ.get("RUN_ID", "").strip() or generate_default_run_id()
        run_dir = (workspace_root / "runs" / eff_id).resolve()

    subdirs = [
        run_dir / "data" / "masks",
        run_dir / "data" / "contours",
        run_dir / "data" / "cropped_patches" / "capitula",
        run_dir / "data" / "tables",
        run_dir / "outputs" / "figures",
        run_dir / "outputs" / "reports",
        run_dir / "logs",
    ]
    for d in subdirs:
        d.mkdir(parents=True, exist_ok=True)

    if update_latest_symlink:
        runs_parent = workspace_root / "runs"
        runs_parent.mkdir(parents=True, exist_ok=True)
        latest_link = runs_parent / "latest"
        try:
            tmp_link = runs_parent / f".latest_{os.getpid()}.tmp"
            if tmp_link.exists() or tmp_link.is_symlink():
                tmp_link.unlink()
            target_rel = run_dir.name if run_dir.parent == runs_parent.resolve() else str(run_dir)
            tmp_link.symlink_to(target_rel)
            tmp_link.replace(latest_link)
        except OSError:
            pass

    if export_env:
        os.environ["PACKERA_RUN_DIR"] = str(run_dir)
        os.environ["RUN_ID"] = run_dir.name

    return run_dir


def ensure_run_curated_vouchers_snapshot(
    run_dir: Path,
    workspace_root: Path = PROJECT_ROOT,
) -> Path:
    """
    Ensure `${RUN_DIR}/data/tables/curated_vouchers.csv` exists before Step 02.
    If absent in the run directory, copy the shared `data/tables/curated_vouchers.csv`
    into the run directory atomically so the run folder is a self-contained snapshot.
    """
    run_vouchers = run_dir / "data" / "tables" / "curated_vouchers.csv"
    if run_vouchers.exists() and run_vouchers.stat().st_size > 0:
        return run_vouchers

    shared_vouchers = workspace_root / "data" / "tables" / "curated_vouchers.csv"
    if shared_vouchers.exists() and shared_vouchers.stat().st_size > 0:
        run_vouchers.parent.mkdir(parents=True, exist_ok=True)
        tmp_vouchers = run_vouchers.with_suffix(".csv.tmp")
        shutil.copy2(shared_vouchers, tmp_vouchers)
        tmp_vouchers.replace(run_vouchers)

    return run_vouchers


def _count_csv_rows(csv_path: Path) -> int:
    """Count data rows (excluding header) in a CSV file."""
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return 0
    try:
        with open(csv_path, "r", encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f)
            total = sum(1 for _ in reader)
            return max(0, total - 1)
    except Exception:
        return 0


def write_run_metadata(
    run_dir: Path,
    run_id: Optional[str] = None,
    start_time: Optional[str] = None,
    end_time: Optional[str] = None,
    cli_args: Optional[List[str]] = None,
    workspace_root: Path = PROJECT_ROOT,
) -> Path:
    """
    Atomically write `${RUN_DIR}/outputs/reports/run_metadata.json` recording:
      - run_id
      - run_dir
      - start_timestamp / end_timestamp
      - git_commit_sha
      - cli_arguments
      - table_row_counts (for all generated CSV tables in data/tables and outputs/reports)
    """
    reports_dir = run_dir / "outputs" / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = reports_dir / "run_metadata.json"

    existing: Dict[str, Any] = {}
    if metadata_path.exists():
        try:
            with open(metadata_path, "r", encoding="utf-8") as f:
                loaded = json.load(f)
                if isinstance(loaded, dict):
                    existing = loaded
        except Exception:
            existing = {}

    git_sha = existing.get("git_commit_sha", "unknown")
    if git_sha == "unknown":
        try:
            res = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                cwd=str(workspace_root),
                capture_output=True,
                text=True,
                timeout=3,
                check=False,
            )
            if res.returncode == 0 and res.stdout.strip():
                git_sha = res.stdout.strip()
        except Exception:
            pass

    now_iso = datetime.now().astimezone().isoformat()
    eff_start = start_time or existing.get("start_timestamp") or now_iso
    eff_end = end_time or now_iso
    eff_run_id = run_id or existing.get("run_id") or run_dir.name
    eff_cli = cli_args if cli_args is not None else existing.get("cli_arguments", list(sys.argv))

    table_counts: Dict[str, int] = {}
    for search_dir in (run_dir / "data" / "tables", run_dir / "outputs" / "reports"):
        if search_dir.exists():
            for csv_file in sorted(search_dir.glob("*.csv")):
                rel_key = str(csv_file.relative_to(run_dir))
                table_counts[rel_key] = _count_csv_rows(csv_file)
                table_counts[csv_file.name] = table_counts[rel_key]

    payload = {
        "run_id": eff_run_id,
        "run_dir": str(run_dir),
        "start_timestamp": eff_start,
        "end_timestamp": eff_end,
        "git_commit_sha": git_sha,
        "cli_arguments": eff_cli,
        "table_row_counts": table_counts,
    }

    tmp_path = metadata_path.with_suffix(".json.tmp")
    with open(tmp_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=False)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    tmp_path.replace(metadata_path)
    return metadata_path


@dataclass(frozen=True)
class PathsConfig:
    """Strongly-typed, project-root anchored filesystem paths."""
    workspace_root: Path = PROJECT_ROOT
    db_path: Path = PROJECT_ROOT / "data" / "packera_pipeline.sqlite"
    raw_vouchers_dir: Path = PROJECT_ROOT / "data" / "raw_vouchers"
    curated_vouchers_csv: Path = PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv"
    failed_qc_vouchers_csv: Path = PROJECT_ROOT / "data" / "tables" / "failed_qc_vouchers.csv"
    extraction_manifest_csv: Path = PROJECT_ROOT / "data" / "tables" / "dataset_manifest.csv"
    contours_dir: Path = PROJECT_ROOT / "data" / "contours"
    masks_dir: Path = PROJECT_ROOT / "data" / "masks"
    leaf_efa_harmonics_csv: Path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics.csv"
    morphometric_flags_csv: Path = PROJECT_ROOT / "data" / "tables" / "morphometric_flags.csv"
    cda_biplot_pdf: Path = PROJECT_ROOT / "outputs" / "figures" / "Figure_Integrative_Packera_Revision.pdf"
    revision_report_md: Path = PROJECT_ROOT / "outputs" / "reports" / "Packera_Taxonomic_Revision_Summary.md"
    gmm_report_csv: Path = PROJECT_ROOT / "outputs" / "reports" / "gmm_bayes_factors_summary.csv"
    summary_log: Path = PROJECT_ROOT / "outputs" / "reports" / "voucher_ingestion_summary.log"
    lm3_root: Path = PROJECT_ROOT / "LeafMachine3"
    models_dir: Path = PROJECT_ROOT / "models" / "lm3"
    model_weights: Path = (
        PROJECT_ROOT / "models" / "lm3" / "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
    )
    quarantine_dir: Path = PROJECT_ROOT / "data" / "raw_vouchers_quarantine"
    tables_dir: Path = PROJECT_ROOT / "data" / "tables"
    figures_dir: Path = PROJECT_ROOT / "outputs" / "figures"


@dataclass(frozen=True)
class TaxaConfig:
    """Target focal taxa, synonymy mappings, and outgroup query sets."""
    target_species: List[str] = field(default_factory=list)
    synonyms: Dict[str, List[str]] = field(default_factory=dict)
    outgroups: List[str] = field(default_factory=list)

    @property
    def target_taxa(self) -> List[str]:
        return self.target_species


@dataclass(frozen=True)
class CapitulumPhenologyConfig:
    """Thresholds for capitulum phenology biomarker classification."""
    min_pappus_luminance: float = 75.0
    max_pappus_saturation: float = 0.18
    min_pappus_area_ratio: float = 0.18
    min_floret_yellow_saturation: float = 0.32
    min_floret_yellow_ratio: float = 0.12
    min_anthesis_heads_for_flowering: int = 1


VALID_GEOGRAPHIC_SCOPES: Set[str] = {"eastern_na", "north_america", "global"}


def normalize_countries_list(raw_countries: Any, default: Optional[List[str]] = None) -> List[str]:
    """Normalize country configuration values (list or comma/whitespace-separated string) to uppercase ISO-2 codes."""
    if raw_countries is None:
        return list(default) if default is not None else ["US", "CA"]
    items: List[str] = []
    if isinstance(raw_countries, str):
        raw_items = [tok for part in raw_countries.split(",") for tok in part.strip().split()]
    elif isinstance(raw_countries, (list, tuple, set)):
        raw_items = [
            tok
            for elem in raw_countries
            for part in str(elem).split(",")
            for tok in part.strip().split()
        ]
    else:
        raw_items = [str(raw_countries).strip()]
    for tok in raw_items:
        code = tok.strip().upper()
        if code and code not in items:
            items.append(code)
    return items


@dataclass(frozen=True)
class ThresholdsConfig:
    """Botanical optical quality and harvesting thresholds."""
    pcd_conf: float = 0.25
    minimum_confidence_threshold: float = 0.65
    min_megapixels: float = 8.0
    min_file_size_kb: float = 500.0
    min_sharpness_laplacian: float = 80.0
    max_uncertainty_meters: float = 5000.0
    western_longitude_threshold: float = -98.0
    geographic_scope: str = "eastern_na"
    countries: List[str] = field(default_factory=lambda: ["US", "CA"])
    capitulum_phenology: CapitulumPhenologyConfig = field(default_factory=CapitulumPhenologyConfig)


@dataclass(frozen=True)
class NormalizationConfig:
    """Fourier elliptic contour normalization invariants."""
    align_major_axis: bool = True
    scale_invariant: bool = True
    rotation_invariant: bool = True
    start_point_invariant: bool = True


@dataclass(frozen=True)
class MorphometricsConfig:
    """Elliptic Fourier Analysis and downstream clustering parameters."""
    nb_harmonics: int = 12
    harmonics: int = 12
    num_pcs: int = 5
    max_k: int = 8
    random_seed: int = 42
    normalization: NormalizationConfig = field(default_factory=NormalizationConfig)
    correct_allometry: str = "auto"
    allometry_r2_threshold: float = 0.10
    symmetrize_efa_harmonics: bool = True


@dataclass(frozen=True)
class PhenologyConfig:
    """Phenological modeling, latitudinal spring progression, and anomaly detection."""
    min_flowering_doy: int = 30
    max_flowering_doy: int = 220
    robust_regression: bool = True


@dataclass(frozen=True)
class MacroecologyConfig:
    """Phase 6 Macroecology, spatial random forest, and phenological parameters."""
    phenology: PhenologyConfig = field(default_factory=PhenologyConfig)



@dataclass(frozen=True)
class HarvestingConfig:
    """Phase 1 GBIF query, download, and filtration hyperparameters."""
    max_records_per_taxon: int = 5000
    max_uncertainty_meters: float = 5000.0
    geographic_scope: str = "eastern_na"
    countries: List[str] = field(default_factory=lambda: ["US", "CA"])
    exclude_western: bool = True
    harvest_full_continental_range: bool = True
    min_eastern_residency_records: int = 3
    western_longitude_threshold: float = -98.0
    home_herbarium: str = "NCU"
    min_megapixels: float = 8.0
    min_file_size_kb: float = 500.0
    min_sharpness_laplacian: float = 80.0
    check_sharpness: bool = False
    webp_quality: int = 92
    download_concurrency: int = 15
    usda_county_filter: bool = False
    usda_county_buffer_km: float = 50.0
    usda_exempt_gold: bool = True


@dataclass(frozen=True)
class SegmentationConfig:
    """Phase 2 LeafMachine3 ONNX segmentation and morphological extraction parameters."""
    cohort_adaptive_mode: bool = True
    min_anchor_n_for_empirical: int = 5
    lm3_root: Path = PROJECT_ROOT / "LeafMachine3"
    models_dir: Path = PROJECT_ROOT / "models" / "lm3"
    model_weights: Path = (
        PROJECT_ROOT / "models" / "lm3" / "lm3_leaf_segmenter__yolo26x_seg_1024.onnx"
    )
    plant_detector_model: str = "lm3_plant_detector__yolo26x_det_1280"
    archival_detector_model: str = "lm3_archival_detector__yolo26x_det_1280"
    specimen_segmenter_model: str = "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"
    leaf_segmenter_model: str = "lm3_leaf_segmenter__yolo26x_seg_1024"
    landmark_detector_model: str = "lm3_landmark_detector__yolo26x_pose_640"
    mp_conversion_model: str = "lm3_mp_conversion_factor__sqrt_fit"
    detector_version: str = "lm3_plant_detector__yolo26x_det_1280"
    max_working_dim: int = 3200
    specimen_max_working_dim: int = 3200
    paperclean: bool = True
    llr_hyst_keep_blobs: bool = True
    paperclean_llr_hyst_keepblobs: bool = True
    subtract_holes: bool = True
    use_cf_predicted_by_mp: bool = False
    ect_num_dirs: int = 360
    minimum_confidence_threshold: float = 0.50
    treat_leaflet_as_leaf_whole: bool = False
    use_efds_for_png_masks: bool = False
    do_apply_landmarks: bool = True
    segment_partial_leaves: bool = False
    rosette_plant_bbox_lower_fraction: float = 0.60
    score_thresh: float = 0.25
    device: str = "cuda"


@dataclass(frozen=True)
class ModelsConfig:
    """LeafMachine3 ONNX model identifiers, checkpoint paths, remote URLs, and integrity hashes."""
    lm3_root: Path = PROJECT_ROOT / "LeafMachine3"
    models_dir: Path = PROJECT_ROOT / "models" / "lm3"
    plant_detector_model: str = "lm3_plant_detector__yolo26x_det_1280"
    archival_detector_model: str = "lm3_archival_detector__yolo26x_det_1280"
    specimen_segmenter_model: str = "lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain"
    leaf_segmenter_model: str = "lm3_leaf_segmenter__yolo26x_seg_1024"
    landmark_detector_model: str = "lm3_landmark_detector__yolo26x_pose_640"
    mp_conversion_model: str = "lm3_mp_conversion_factor__sqrt_fit"
    pcd_weights_path: Path = (
        PROJECT_ROOT / "models" / "lm3" / "lm3_plant_detector__yolo26x_det_1280.onnx"
    )
    pcd_weights_url: str = ""
    pcd_weights_sha256: str = ""


@dataclass(frozen=True)
class ComputeConfig:
    """CPU multithreading and GPU inference batch size configuration."""
    threads: Optional[int] = None
    gpu_batch_size: Optional[int] = None


@dataclass
class PipelineConfig:
    """Centralized, type-hinted configuration for Packera morphometrics pipeline."""
    paths: PathsConfig
    taxa: TaxaConfig
    thresholds: ThresholdsConfig
    morphometrics: MorphometricsConfig
    harvesting: HarvestingConfig
    segmentation: SegmentationConfig
    models: ModelsConfig = field(default_factory=ModelsConfig)
    macroecology: MacroecologyConfig = field(default_factory=MacroecologyConfig)
    compute: ComputeConfig = field(default_factory=ComputeConfig)
    _raw_dict: Dict[str, Any] = field(default_factory=dict, repr=False)

    def __getitem__(self, key: str) -> Any:
        return self._raw_dict[key]

    def __contains__(self, key: str) -> bool:
        return key in self._raw_dict

    def get(self, key: str, default: Any = None) -> Any:
        return self._raw_dict.get(key, default)

    def to_dict(self) -> Dict[str, Any]:
        return dict(self._raw_dict)

    def with_run_dir(self, run_dir: Union[str, Path]) -> PipelineConfig:
        """
        Return a new PipelineConfig whose run-scoped output/table/figure paths
        point inside `run_dir`, while shared read-only caches (raw_vouchers_dir,
        quarantine_dir, models_dir, model_weights) remain anchored at PROJECT_ROOT.
        """
        rdir = Path(run_dir).resolve()
        new_paths = replace(
            self.paths,
            curated_vouchers_csv=rdir / "data" / "tables" / "curated_vouchers.csv",
            failed_qc_vouchers_csv=rdir / "data" / "tables" / "failed_qc_vouchers.csv",
            extraction_manifest_csv=rdir / "data" / "tables" / "extracted_leaf_manifest.csv",
            contours_dir=rdir / "data" / "contours",
            masks_dir=rdir / "data" / "masks",
            leaf_efa_harmonics_csv=rdir / "data" / "tables" / "leaf_efa_harmonics.csv",
            morphometric_flags_csv=rdir / "data" / "tables" / "morphometric_flags.csv",
            cda_biplot_pdf=rdir / "outputs" / "figures" / "Figure_Integrative_Packera_Revision.pdf",
            revision_report_md=rdir / "outputs" / "reports" / "Packera_Taxonomic_Revision_Summary.md",
            gmm_report_csv=rdir / "outputs" / "reports" / "gmm_bayes_factors_summary.csv",
            summary_log=rdir / "logs" / "pipeline.log",
            tables_dir=rdir / "data" / "tables",
            figures_dir=rdir / "outputs" / "figures",
        )
        new_raw = dict(self._raw_dict)
        raw_paths = dict(new_raw.get("paths", {}))
        for k in RUN_SCOPED_PATH_FIELDS:
            raw_paths[k] = str(getattr(new_paths, k))
        new_raw["paths"] = raw_paths
        new_raw["resolved_paths"] = {
            k: getattr(new_paths, k) for k in new_paths.__dataclass_fields__
        }
        return replace(self, paths=new_paths, _raw_dict=new_raw)

    @classmethod
    def from_yaml(
        cls,
        yaml_path: Optional[Union[str, Path]] = None,
        run_dir: Optional[Union[str, Path]] = None,
    ) -> PipelineConfig:
        """Dynamically parse YAML and resolve paths relative to project root (or active run_dir if supplied)."""
        target = Path(yaml_path) if yaml_path else DEFAULT_CONFIG_PATH
        target = target if target.is_absolute() else (PROJECT_ROOT / target).resolve()

        if not target.exists():
            raise FileNotFoundError(f"Configuration file not found at: {target}")

        try:
            with open(target, "r", encoding="utf-8") as f:
                raw_cfg = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            raise ValueError(f"Invalid YAML syntax in configuration file '{target}': {exc}") from exc

        if not isinstance(raw_cfg, dict) or not raw_cfg:
            raise ValueError(f"Configuration file at '{target}' is empty or not a valid dictionary mapping.")

        def _resolve(rel_path: Union[str, Path]) -> Path:
            p = Path(rel_path)
            return p if p.is_absolute() else (PROJECT_ROOT / p).resolve()

        p_kwargs = {
            k: _resolve(v) for k, v in raw_cfg.get("paths", {}).items()
            if k in PathsConfig.__dataclass_fields__
        }
        paths = PathsConfig(**p_kwargs)

        taxa_dict = raw_cfg.get("taxa", {}) or {}
        legacy_taxa = (raw_cfg.get("taxonomy") or {}).get("target_taxa", [])
        taxa = TaxaConfig(
            target_species=list(taxa_dict.get("target_species", legacy_taxa) or []),
            synonyms=dict(taxa_dict.get("synonyms", {}) or {}),
            outgroups=list(taxa_dict.get("outgroups", []) or []),
        )

        thresh_dict = raw_cfg.get("thresholds", {})
        pheno_dict = thresh_dict.get("capitulum_phenology")
        if isinstance(pheno_dict, dict):
            pheno_cfg = CapitulumPhenologyConfig(
                min_pappus_luminance=float(pheno_dict.get("min_pappus_luminance", 75.0)),
                max_pappus_saturation=float(pheno_dict.get("max_pappus_saturation", 0.18)),
                min_pappus_area_ratio=float(pheno_dict.get("min_pappus_area_ratio", 0.18)),
                min_floret_yellow_saturation=float(pheno_dict.get("min_floret_yellow_saturation", 0.32)),
                min_floret_yellow_ratio=float(pheno_dict.get("min_floret_yellow_ratio", 0.12)),
                min_anthesis_heads_for_flowering=int(pheno_dict.get("min_anthesis_heads_for_flowering", 1)),
            )
        else:
            pheno_cfg = CapitulumPhenologyConfig()

        harv_dict = raw_cfg.get("harvesting", {})
        non_float_thresh_fields = {"capitulum_phenology", "geographic_scope", "countries"}
        thresh_kwargs: Dict[str, Any] = {
            k: float(v)
            for k, v in thresh_dict.items()
            if k in ThresholdsConfig.__dataclass_fields__
            and k not in non_float_thresh_fields
            and not isinstance(v, (dict, list))
        }
        raw_t_scope = str(thresh_dict.get("geographic_scope", harv_dict.get("geographic_scope", "eastern_na"))).strip().lower()
        thresh_kwargs["geographic_scope"] = raw_t_scope if raw_t_scope in VALID_GEOGRAPHIC_SCOPES else "eastern_na"
        thresh_kwargs["countries"] = normalize_countries_list(
            thresh_dict.get("countries", harv_dict.get("countries", ["US", "CA"])),
            default=["US", "CA"],
        )
        thresh_kwargs["capitulum_phenology"] = pheno_cfg

        thresholds = ThresholdsConfig(**thresh_kwargs)

        norm_dict = raw_cfg.get("morphometrics", {}).get("normalization", {})
        norm = NormalizationConfig(**{
            k: bool(v) for k, v in norm_dict.items() if k in NormalizationConfig.__dataclass_fields__
        })

        m_dict = raw_cfg.get("morphometrics", {})
        nb_h = int(m_dict.get("nb_harmonics", m_dict.get("harmonics", 12)))
        morph = MorphometricsConfig(
            nb_harmonics=nb_h, harmonics=int(m_dict.get("harmonics", nb_h)),
            num_pcs=int(m_dict.get("num_pcs", 5)), max_k=int(m_dict.get("max_k", 8)),
            random_seed=int(m_dict.get("random_seed", 42)), normalization=norm,
            correct_allometry=str(m_dict.get("correct_allometry", "auto")),
            allometry_r2_threshold=float(m_dict.get("allometry_r2_threshold", 0.10)),
            symmetrize_efa_harmonics=bool(m_dict.get("symmetrize_efa_harmonics", True)),
        )

        harv_kwargs: Dict[str, Any] = {
            k: v for k, v in harv_dict.items() if k in HarvestingConfig.__dataclass_fields__
        }
        raw_h_scope = str(harv_kwargs.get("geographic_scope", thresh_kwargs["geographic_scope"])).strip().lower()
        harv_kwargs["geographic_scope"] = raw_h_scope if raw_h_scope in VALID_GEOGRAPHIC_SCOPES else "eastern_na"
        if "countries" in harv_kwargs:
            harv_kwargs["countries"] = normalize_countries_list(harv_kwargs["countries"], default=["US", "CA"])
        elif harv_kwargs["geographic_scope"] == "north_america":
            harv_kwargs["countries"] = ["US", "CA", "MX"]
        elif harv_kwargs["geographic_scope"] == "global":
            harv_kwargs["countries"] = []
        if "home_herbarium" in harv_kwargs:
            harv_kwargs["home_herbarium"] = str(harv_kwargs["home_herbarium"]).strip().upper() or "NCU"
        if "usda_county_filter" in harv_kwargs:
            harv_kwargs["usda_county_filter"] = bool(harv_kwargs["usda_county_filter"])
        if "usda_county_buffer_km" in harv_kwargs and harv_kwargs["usda_county_buffer_km"] is not None:
            harv_kwargs["usda_county_buffer_km"] = float(harv_kwargs["usda_county_buffer_km"])
        if "usda_exempt_gold" in harv_kwargs:
            harv_kwargs["usda_exempt_gold"] = bool(harv_kwargs["usda_exempt_gold"])
        harvesting = HarvestingConfig(**harv_kwargs)

        leaf_seg_fallback = raw_cfg.get("leaf_segmentation", {}) or {}
        seg_dict = {**leaf_seg_fallback, **(raw_cfg.get("segmentation", {}) or {})}
        seg_kwargs: Dict[str, Any] = {
            k: v for k, v in seg_dict.items() if k in SegmentationConfig.__dataclass_fields__
        }
        for path_k in ("model_weights", "lm3_root", "models_dir"):
            if path_k in seg_kwargs and seg_kwargs[path_k] is not None:
                seg_kwargs[path_k] = _resolve(seg_kwargs[path_k])
        for bool_k in (
            "cohort_adaptive_mode",
            "paperclean",
            "llr_hyst_keep_blobs",
            "paperclean_llr_hyst_keepblobs",
            "subtract_holes",
            "use_cf_predicted_by_mp",
            "treat_leaflet_as_leaf_whole",
            "use_efds_for_png_masks",
            "do_apply_landmarks",
            "segment_partial_leaves",
        ):
            if bool_k in seg_kwargs and seg_kwargs[bool_k] is not None:
                seg_kwargs[bool_k] = bool(seg_kwargs[bool_k])
        for int_k in ("min_anchor_n_for_empirical", "max_working_dim", "specimen_max_working_dim", "ect_num_dirs"):
            if int_k in seg_kwargs and seg_kwargs[int_k] is not None:
                seg_kwargs[int_k] = int(seg_kwargs[int_k])
        for float_k in ("minimum_confidence_threshold", "rosette_plant_bbox_lower_fraction", "score_thresh"):
            if float_k in seg_kwargs and seg_kwargs[float_k] is not None:
                seg_kwargs[float_k] = float(seg_kwargs[float_k])
        segmentation = SegmentationConfig(**seg_kwargs)

        models_dict = raw_cfg.get("models", {})
        models_kwargs: Dict[str, Any] = {
            k: v for k, v in models_dict.items() if k in ModelsConfig.__dataclass_fields__
        }
        for path_k in ("pcd_weights_path", "lm3_root", "models_dir"):
            if path_k in models_kwargs and models_kwargs[path_k] is not None:
                models_kwargs[path_k] = _resolve(models_kwargs[path_k])
        for str_k in (
            "plant_detector_model",
            "archival_detector_model",
            "specimen_segmenter_model",
            "leaf_segmenter_model",
            "landmark_detector_model",
            "mp_conversion_model",
            "pcd_weights_url",
            "pcd_weights_sha256",
        ):
            if str_k in models_kwargs and models_kwargs[str_k] is not None:
                models_kwargs[str_k] = str(models_kwargs[str_k])
        models = ModelsConfig(**models_kwargs)

        macro_dict = raw_cfg.get("macroecology", {})
        pheno_dict = macro_dict.get("phenology", {})
        phenology_cfg = PhenologyConfig(
            min_flowering_doy=int(pheno_dict.get("min_flowering_doy", 30)),
            max_flowering_doy=int(pheno_dict.get("max_flowering_doy", 220)),
            robust_regression=bool(pheno_dict.get("robust_regression", True)),
        )
        macroecology = MacroecologyConfig(phenology=phenology_cfg)

        comp_dict = raw_cfg.get("compute", {}) or {}
        raw_threads = comp_dict.get("threads")
        raw_gpu_bs = comp_dict.get("gpu_batch_size")
        compute_cfg = ComputeConfig(
            threads=int(raw_threads) if raw_threads is not None and int(raw_threads) >= 1 else None,
            gpu_batch_size=int(raw_gpu_bs) if raw_gpu_bs is not None and int(raw_gpu_bs) >= 1 else None,
        )

        resolved_raw = dict(raw_cfg)
        resolved_raw["resolved_paths"] = {
            k: getattr(paths, k) for k in paths.__dataclass_fields__
        }

        instance = cls(
            paths=paths, taxa=taxa, thresholds=thresholds,
            morphometrics=morph, harvesting=harvesting,
            segmentation=segmentation, models=models,
            macroecology=macroecology, compute=compute_cfg,
            _raw_dict=resolved_raw,
        )
        if run_dir is not None:
            return instance.with_run_dir(run_dir)
        return instance


# Botanical ontology mappings for Packera phenotyping
CLASS_NAMES: List[str] = [
    "basal_leaf_blade", "leaf_petiole", "cauline_leaf",
    "cauline_stem", "root_rhizome", "basal_rosette_clump", "capitulum",
]
CLASS_MAP: Dict[str, int] = {name: idx for idx, name in enumerate(CLASS_NAMES)}
CLASS_COLORS_BGR: Dict[int, Tuple[int, int, int]] = {
    0: (0, 200, 0), 1: (50, 255, 150), 2: (0, 165, 255),
    3: (0, 215, 255), 4: (50, 100, 200), 5: (0, 140, 70), 6: (230, 30, 230),
}

# Layer 1: Genus-wide North American Packera / Senecioneae monographers & initial-qualified flora authors
# (Common surnames like Freeman, Turner, Fuller, Weber, Golden, Packer, Estes, Welsh, Dorn, Heller, Britton
#  are strictly initial/given-name-qualified; checklist compilers like Kartesz are excluded.)
SPECIALIST_PATTERNS: List[str] = [
    r"\bBarkley\b", r"\bT\.?\s*M\.?\s*Barkley\b", r"\bTheodore\s+M\.?\s+Barkley\b",
    r"\bTrock\b", r"\bD\.?\s*K\.?\s*Trock\b", r"\bDebra\s+K\.?\s+Trock\b",
    r"\bKowal\b", r"\bR\.?\s*R\.?\s*Kowal\b", r"\bRobert\s+R\.?\s+Kowal\b",
    r"\bWeakley\b", r"\bA\.?\s*S\.?\s*Weakley\b", r"\bAlan\s+S\.?\s+Weakley\b",
    r"\bBain\b", r"\bJ\.?\s*F\.?\s*Bain\b", r"\bJohn\s+F\.?\s+Bain\b",
    r"\bMahoney\b", r"\bA\.?\s*M\.?\s*Mahoney\b", r"\bAlison\s+M\.?\s+Mahoney\b",
    r"\bSoderholm\b", r"\bJ\.?\s*Soderholm\b",
    r"\bJ\.?\s*G\.?\s*Packer\b|\bJohn\s+G\.?\s+Packer\b|\bPacker\s*,\s*(?:J\.?\s*G\.?|John)\b",
    r"\bC\.?\s*Jeffrey\b|\bCharles\s+Jeffrey\b|\bJeffrey\s*,\s*(?:C\.?|Charles)\b",
    r"\bMabberley\b", r"\bD\.?\s*J\.?\s*Mabberley\b",
    r"\b(?:J\.?\s*B\.?|(?:J\.?\s+)?Brandon)\s+Fuller\b|\bFuller\s*,\s*(?:J\.?\s*B\.?|J\.?\s*Brandon|Brandon)\b",
    r"\bMandel\b", r"\bJ\.?\s*R\.?\s*Mandel\b", r"\bJennifer\s+R\.?\s+Mandel\b",
    r"\bMoore[-\s]?Pollard\b",
    r"\b(?:E\.?\s*R\.?|Erika\s+R\.?)\s+Moore\b|\bMoore\s*,\s*(?:E\.?\s*R\.?|Erika)\b",
    r"\bBoufford\b", r"\bD\.?\s*E\.?\s*Boufford\b",
    r"\bGramling\b", r"\bA\.?\s*Gramling\b",
    r"\bC\.?\s*C\.?\s*Freeman\b|\bCraig\s+(?:C\.?\s+)?Freeman\b|\bFreeman\s*,\s*(?:C\.?\s*C\.?|Craig)\b",
    r"\bCronquist\b", r"\bA\.?\s*Cronquist\b", r"\bArthur\s+(?:J\.?\s+)?Cronquist\b",
    r"\bGreenman\b", r"\bJ\.?\s*M\.?\s*Greenman\b", r"\bJesse\s+M(?:ore)?\.?\s+Greenman\b",
    r"\bW\.?\s*A\.?\s*Weber\b|\bWilliam\s+A\.?\s+Weber\b|\bWeber\s*,\s*(?:W\.?\s*A\.?|William)\b",
    r"(?:^|[\s|,;(])[ÁA]\.?\s*L[öo]ve\b", r"\bD\.?\s*L[öo]ve\b",
    r"\b[ÁA]skell\s+L[öo]ve\b", r"\bDoris\s+L[öo]ve\b",
    r"\bSpellenberg\b", r"\bR\.?\s*(?:W\.?\s*)?Spellenberg\b", r"\bRichard\s+Spellenberg\b",
    r"\bPippen\b", r"\bR\.?\s*W\.?\s*Pippen\b", r"\bRichard\s+W\.?\s+Pippen\b",
    r"\bUttal\b", r"\bL\.?\s*J\.?\s*Uttal\b", r"\bLeonard\s+J\.?\s+Uttal\b",
    r"\bJudziewicz\b", r"\bE\.?\s*J\.?\s*Judziewicz\b",
    r"\bO['’]?Kennon\b", r"\bR\.?\s*J\.?\s*O['’]?Kennon\b",
    r"\bS\.?\s*L\.?\s*Welsh\b|\bStanley\s+L\.?\s+Welsh\b|\bWelsh\s*,\s*(?:S\.?\s*L\.?|Stanley)\b",
    r"\bN\.?\s*D\.?\s*Atwood\b|\bAtwood\s*,\s*(?:N\.?\s*D\.?|Duane)\b",
    r"\bR\.?\s*D\.?\s*Dorn\b|\bRobert\s+D\.?\s+Dorn\b|\bDorn\s*,\s*(?:R\.?\s*D\.?|Robert)\b",
    r"\bO['’]?Kane\b", r"\bS\.?\s*L\.?\s*O['’]?Kane\b",
    r"\bYeatts\b", r"\bL\.?\s*Yeatts\b",
    r"\bRuyle[-\s]?Douglas\b",
    r"\b(?:G\.?\s*W\.?|George\s+W\.?)\s+Douglas\b|\bDouglas\s*,\s*(?:G\.?\s*W\.?|George)\b",
    r"\bBoivin\b", r"\bB\.?\s*Boivin\b",
    r"\bFernald\b", r"\bM\.?\s*L\.?\s*Fernald\b", r"\bMerritt\s+L(?:yndon)?\.?\s+Fernald\b",
    r"\bRydberg\b", r"\bP\.?\s*A\.?\s*Rydberg\b", r"\bPer\s+Axel\s+Rydberg\b",
    r"\bN\.?\s*L\.?\s*Britton\b|\bNathaniel\s+L(?:ord)?\.?\s+Britton\b|\bBritton\s*,\s*(?:N\.?\s*L\.?|Nathaniel)\b",
    r"\b(?:E\.?\s*L\.?|Edward\s+Lee)\s+Greene\b|\bGreene\s*,\s*(?:E\.?\s*L\.?|Edward)\b",
    r"\b(?:A\.?|Asa)\s+Gray\b|\bGray\s*,\s*(?:A\.?|Asa)\b",
    r"\b(?:J\.?\s*K\.?|John\s+Kunkel)\s+Small\b|\bSmall\s*,\s*(?:J\.?\s*K\.?|John)\b",
    r"\b(?:Alph\.?|Alphonso)\s+Wood\b|\bWood\s*,\s*(?:Alph\.?|Alphonso)\b",
    r"\b(?:S\.?\s*F\.?|Sidney\s+Fay)\s+Blake\b|\bBlake\s*,\s*(?:S\.?\s*F\.?|sidney)\b",
    r"\b(?:S\.?|Sereno)\s+Watson\b|\bWatson\s*,\s*(?:S\.?|Sereno)\b",
    r"\bA\.?\s*A\.?\s*Heller\b|\bAmos\s+A(?:rthur)?\.?\s+Heller\b|\bHeller\s*,\s*(?:A\.?\s*A\.?|Amos)\b",
    r"\bKral\b", r"\bR\.?\s*Kral\b", r"\bRobert\s+Kral\b",
    r"\bGodfrey\b", r"\bR\.?\s*K\.?\s*Godfrey\b", r"\bRobert\s+K\.?\s+Godfrey\b",
    r"\bSorrie\b", r"\bB\.?\s*A\.?\s*Sorrie\b", r"\bBruce\s+A\.?\s+Sorrie\b",
    r"\bNesom\b", r"\bG\.?\s*L\.?\s*Nesom\b", r"\bGuy\s+L\.?\s+Nesom\b",
    r"\bStrother\b", r"\bJ\.?\s*L\.?\s*Strother\b", r"\bJohn\s+L\.?\s+Strother\b",
    r"\bB\.?\s*L\.?\s*Turner\b|\bBillie\s+L(?:ee)?\.?\s+Turner\b|\bTurner\s*,\s*(?:B\.?\s*L\.?|Billie)\b",
    r"\bJ\.?\s*L\.?\s*Golden\b|\bJudith\s+L\.?\s+Golden\b|\bGolden\s*,\s*(?:J\.?\s*L\.?|Judith)\b",
    r"\bBrouillet\b", r"\bL\.?\s*Brouillet\b",
    r"\bSemple\b", r"\bJ\.?\s*C\.?\s*Semple\b",
    r"\bSchilling\b", r"\bE\.?\s*E\.?\s*Schilling\b",
    r"\bFloden\b", r"\bA\.?\s*J\.?\s*Floden\b",
    r"\bL\.?\s*D\.?\s*Estes\b|\bDwayne\s+Estes\b|\bEstes\s*,\s*(?:L\.?\s*D\.?|Dwayne)\b",
    r"\bPoindexter\b", r"\bD\.?\s*B\.?\s*Poindexter\b",
    r"\bWieboldt\b", r"\bT\.?\s*F\.?\s*Wieboldt\b",
    r"\bMcMillan\b", r"\bP\.?\s*D\.?\s*McMillan\b",
    r"\bBeauchamp\b", r"\bR\.?\s*M\.?\s*Beauchamp\b",
]
MAJOR_HERBARIA_CODES: Set[str] = {
    "NCU", "DUKE", "NCSC", "GA", "US", "NY", "BRIT", "MO", "WIS", "VDB", "FLAS", "TEX", "LL", "TENN", "F",
    "GH", "PH", "KANU", "MSC", "COLO", "RM", "UC", "JEPS", "CAS", "DAO", "CAN", "MEXU",
    "USCH", "CLEMS", "VPI", "WILLI", "WVA", "CM", "MICH", "MIN", "OS", "ILL", "ILLS", "IND",
    "NEBC", "VT", "TRT", "MT", "QFA", "ACAD", "UNB", "UARK", "ARK", "LSU", "MISS", "UNA",
    "APSC", "EKY", "KY", "MU", "PAC", "PENN", "CONN", "NYS", "BH",
}
PRIMARY_TYPE_STATUSES: Set[str] = {
    "HOLOTYPE", "LECTOTYPE", "NEOTYPE", "EPITYPE", "ISOTYPE", "ISOLECTOTYPE", "ISONEOTYPE"
}
SECONDARY_TYPE_STATUSES: Set[str] = {
    "SYNTYPE", "ISOSYNTYPE", "PARATYPE", "ISOPARATYPE", "TYPE", "COTYPE"
}
VALID_TYPE_STATUSES: Set[str] = PRIMARY_TYPE_STATUSES | SECONDARY_TYPE_STATUSES
EXCLUDED_WESTERN_STATES: Set[str] = {
    "CO", "COLORADO", "NM", "NEW MEXICO", "WY", "WYOMING", "MT", "MONTANA",
    "UT", "UTAH", "AZ", "ARIZONA", "NV", "NEVADA", "ID", "IDAHO",
    "WA", "WASHINGTON", "OR", "OREGON", "CA", "CALIFORNIA", "AK", "ALASKA", "HI", "HAWAII",
    "BC", "BRITISH COLUMBIA", "AB", "ALBERTA", "SK", "SASKATCHEWAN",
    "YT", "YUKON", "YUKON TERRITORY", "NT", "NORTHWEST TERRITORIES", "NU", "NUNAVUT",
}

# Module-level defaults instantiated dynamically from configuration YAML
_DEFAULT_CONFIG = PipelineConfig.from_yaml()
DEFAULT_WORKSPACE = _DEFAULT_CONFIG.paths.workspace_root
DEFAULT_DB_PATH = _DEFAULT_CONFIG.paths.db_path
DEFAULT_RAW_DIR = _DEFAULT_CONFIG.paths.raw_vouchers_dir
DEFAULT_CURATED_CSV = _DEFAULT_CONFIG.paths.curated_vouchers_csv
DEFAULT_OUTPUT_CSV = _DEFAULT_CONFIG.paths.curated_vouchers_csv
DEFAULT_SUMMARY_LOG = _DEFAULT_CONFIG.paths.summary_log
DEFAULT_QUARANTINE_DIR = _DEFAULT_CONFIG.paths.quarantine_dir
DEFAULT_MIN_MEGAPIXELS = _DEFAULT_CONFIG.thresholds.min_megapixels
DEFAULT_MIN_FILE_SIZE_KB = _DEFAULT_CONFIG.thresholds.min_file_size_kb
DEFAULT_MIN_SHARPNESS_LAPLACIAN = _DEFAULT_CONFIG.thresholds.min_sharpness_laplacian
DEFAULT_TARGET_TAXA: List[str] = []
WESTERN_LONGITUDE_THRESHOLD = _DEFAULT_CONFIG.thresholds.western_longitude_threshold
DEFAULT_GEOGRAPHIC_SCOPE = _DEFAULT_CONFIG.harvesting.geographic_scope
DEFAULT_COUNTRIES = _DEFAULT_CONFIG.harvesting.countries
DEFAULT_HOME_HERBARIUM = _DEFAULT_CONFIG.harvesting.home_herbarium
DEFAULT_USDA_COUNTY_FILTER = _DEFAULT_CONFIG.harvesting.usda_county_filter
DEFAULT_USDA_COUNTY_BUFFER_KM = _DEFAULT_CONFIG.harvesting.usda_county_buffer_km
DEFAULT_USDA_EXEMPT_GOLD = _DEFAULT_CONFIG.harvesting.usda_exempt_gold


def resolve_cpu_threads(
    cli_val: Optional[int] = None,
    cfg: Optional[PipelineConfig] = None,
    default: Optional[int] = None,
) -> Optional[int]:
    """
    Resolve CPU worker thread count via priority chain:
    CLI (--threads / -j) -> ENV (PACKERA_MAX_WORKERS) -> config.yaml (compute.threads) -> default.
    """
    if cli_val is not None:
        try:
            v = int(cli_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    env_val = os.environ.get("PACKERA_MAX_WORKERS", "").strip()
    if env_val:
        try:
            v = int(env_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    active_cfg = cfg if cfg is not None else _DEFAULT_CONFIG
    if active_cfg.compute.threads is not None and int(active_cfg.compute.threads) >= 1:
        return int(active_cfg.compute.threads)
    return default


def resolve_gpu_batch_size(
    cli_val: Optional[int] = None,
    cfg: Optional[PipelineConfig] = None,
    default_cuda: int = 64,
    default_cpu: int = 32,
    use_cuda: bool = True,
    default: Optional[int] = None,
) -> int:
    """
    Resolve GPU/inference batch size via priority chain:
    CLI (--gpu-batch-size / --batch-size) -> ENV (PACKERA_GPU_BATCH_SIZE) ->
    config.yaml (compute.gpu_batch_size) -> stage-optimal default.
    """
    if cli_val is not None:
        try:
            v = int(cli_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    env_val = os.environ.get("PACKERA_GPU_BATCH_SIZE", "").strip()
    if env_val:
        try:
            v = int(env_val)
            if v >= 1:
                return v
        except (TypeError, ValueError):
            pass
    active_cfg = cfg if cfg is not None else _DEFAULT_CONFIG
    if active_cfg.compute.gpu_batch_size is not None and int(active_cfg.compute.gpu_batch_size) >= 1:
        return int(active_cfg.compute.gpu_batch_size)
    if default is not None:
        return max(1, int(default))
    return max(1, int(default_cuda) if use_cuda else int(default_cpu))
