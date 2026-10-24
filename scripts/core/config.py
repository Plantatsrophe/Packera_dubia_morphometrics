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
      - data/masks/{tier1_intact,tier2_reflected}
      - data/contours
      - data/cropped_patches/{basal_leaves_raw,capitula,rosettes_dense}
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
        run_dir / "data" / "masks" / "tier1_intact",
        run_dir / "data" / "masks" / "tier2_reflected",
        run_dir / "data" / "contours",
        run_dir / "data" / "cropped_patches" / "basal_leaves_raw",
        run_dir / "data" / "cropped_patches" / "capitula",
        run_dir / "data" / "cropped_patches" / "rosettes_dense",
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
    raw_vouchers_dir: Path = PROJECT_ROOT / "data" / "raw_vouchers"
    curated_vouchers_csv: Path = PROJECT_ROOT / "data" / "tables" / "curated_vouchers.csv"
    failed_qc_vouchers_csv: Path = PROJECT_ROOT / "data" / "tables" / "failed_qc_vouchers.csv"
    extraction_manifest_csv: Path = PROJECT_ROOT / "data" / "tables" / "dataset_manifest.csv"
    contours_dir: Path = PROJECT_ROOT / "data" / "contours"
    masks_dir: Path = PROJECT_ROOT / "data" / "masks"
    leaf_efa_harmonics_csv: Path = PROJECT_ROOT / "data" / "tables" / "leaf_efa_harmonics.csv"
    morphometric_flags_csv: Path = PROJECT_ROOT / "data" / "tables" / "morphometric_flags.csv"
    cda_biplot_pdf: Path = PROJECT_ROOT / "outputs" / "figures" / "cda_passive_projection.pdf"
    gmm_report_csv: Path = PROJECT_ROOT / "outputs" / "reports" / "gmm_bayes_factors_summary.csv"
    summary_log: Path = PROJECT_ROOT / "outputs" / "reports" / "voucher_ingestion_summary.log"
    model_weights: Path = PROJECT_ROOT / "models" / "lm2_packera_pcd_finetuned.pth"
    quarantine_dir: Path = PROJECT_ROOT / "data" / "raw_vouchers_quarantine"
    tables_dir: Path = PROJECT_ROOT / "data" / "tables"
    figures_dir: Path = PROJECT_ROOT / "outputs" / "figures"
    models_dir: Path = PROJECT_ROOT / "models"


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
class FoldDetectionConfig:
    """Empirical straight-chord fold detection parameters (calibrated)."""
    enabled: bool = True
    max_folded_aspect_ratio: float = 0.42
    max_chord_deviation_ratio: float = 0.035
    min_chord_length_px: int = 150


@dataclass(frozen=True)
class DissectionConfig:
    """Empirical dissection & lyrate sinus parameters (calibrated)."""
    enabled: bool = True
    min_solidity_dissected: float = 0.50
    sinus_defect_min_depth_ratio: float = 0.08
    min_bilateral_sinus_count: int = 3
    taxa_with_lyrate_tendency: List[str] = field(default_factory=lambda: [
        "Packera paupercula",
        "Packera plattensis",
        "Packera paupercula var. paupercula",
        "Packera paupercula var. savannarum",
    ])

    @property
    def default(self) -> float:
        """Backward-compatible alias for standard baseline min_solidity."""
        return 0.72

    @property
    def min_dissected(self) -> float:
        """Backward-compatible alias for min_solidity_dissected."""
        return self.min_solidity_dissected


# Backward-compatibility alias for legacy callers importing SolidityThresholdsConfig
SolidityThresholdsConfig = DissectionConfig


@dataclass(frozen=True)
class CapitulumPhenologyConfig:
    """Thresholds for capitulum phenology biomarker classification."""
    min_pappus_luminance: float = 75.0
    max_pappus_saturation: float = 0.18
    min_pappus_area_ratio: float = 0.18
    min_floret_yellow_saturation: float = 0.32
    min_floret_yellow_ratio: float = 0.12
    min_anthesis_heads_for_flowering: int = 1


@dataclass(frozen=True)
class ThresholdsConfig:
    """Botanical optical quality, geometry, and detection thresholds."""
    pcd_conf: float = 0.72
    min_solidity: float = 0.72
    min_ucs: float = 0.85
    cleanlab_cutoff: float = 0.85
    min_megapixels: float = 8.0
    min_file_size_kb: float = 500.0
    min_sharpness_laplacian: float = 80.0
    max_uncertainty_meters: float = 5000.0
    western_longitude_threshold: float = -106.65
    fold_detection: FoldDetectionConfig = field(default_factory=FoldDetectionConfig)
    dissection: DissectionConfig = field(default_factory=DissectionConfig)
    capitulum_phenology: CapitulumPhenologyConfig = field(default_factory=CapitulumPhenologyConfig)

    @property
    def solidity(self) -> DissectionConfig:
        """Backward-compatible alias for dissection thresholds."""
        return self.dissection


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


@dataclass(frozen=True)
class PhenologyConfig:
    """Phenological modeling, latitudinal spring progression, and anomaly detection."""
    min_flowering_doy: int = 60
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
    exclude_western: bool = True
    western_longitude_threshold: float = -106.65
    min_megapixels: float = 8.0
    min_file_size_kb: float = 500.0
    min_sharpness_laplacian: float = 80.0
    download_concurrency: int = 15


@dataclass(frozen=True)
class SegmentationConfig:
    """Phase 2 PointRend deep segmentation parameters."""
    model_weights: Path = (
        PROJECT_ROOT
        / "LeafMachine2"
        / "leafmachine2"
        / "segmentation"
        / "models"
        / "Group3_Dataset_100000_Iter_1176PTS_512Batch_smooth_l1_LR00025_BGR"
        / "model_final.pth"
    )
    score_thresh: float = 0.35
    min_ucs: float = 0.85
    min_solidity: float = 0.72
    device: str = "cuda"


@dataclass(frozen=True)
class ModelsConfig:
    """Model weights checkpoint paths, remote URLs, and integrity hashes."""
    pcd_weights_path: Path = PROJECT_ROOT / "models" / "lm2_packera_pcd_finetuned.pth"
    pcd_weights_url: str = (
        "https://github.com/Plantatsrophe/Packera_dubia_morphometrics/releases/download/v1.0-weights/lm2_packera_pcd_finetuned.pth"
    )
    pcd_weights_sha256: str = ""


@dataclass
class PipelineConfig:
    """Centralized, type-hinted configuration for Packera dubia pipeline."""
    paths: PathsConfig
    taxa: TaxaConfig
    thresholds: ThresholdsConfig
    morphometrics: MorphometricsConfig
    harvesting: HarvestingConfig
    segmentation: SegmentationConfig
    models: ModelsConfig = field(default_factory=ModelsConfig)
    macroecology: MacroecologyConfig = field(default_factory=MacroecologyConfig)
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
            cda_biplot_pdf=rdir / "outputs" / "figures" / "cda_passive_projection.pdf",
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

        taxa_dict = raw_cfg.get("taxa", {})
        legacy_taxa = raw_cfg.get("taxonomy", {}).get("target_taxa", [])
        taxa = TaxaConfig(
            target_species=taxa_dict.get("target_species", legacy_taxa),
            synonyms=taxa_dict.get("synonyms", {}),
            outgroups=taxa_dict.get("outgroups", []),
        )

        thresh_dict = raw_cfg.get("thresholds", {})

        # Parse dissection parameters with safe fallbacks
        diss_dict = thresh_dict.get("dissection") or thresh_dict.get("solidity")
        if isinstance(diss_dict, dict):
            dissection_cfg = DissectionConfig(
                enabled=bool(diss_dict.get("enabled", True)),
                min_solidity_dissected=float(
                    diss_dict.get("min_solidity_dissected", diss_dict.get("min_dissected", 0.50))
                ),
                sinus_defect_min_depth_ratio=float(
                    diss_dict.get("sinus_defect_min_depth_ratio", 0.08)
                ),
                min_bilateral_sinus_count=int(
                    diss_dict.get("min_bilateral_sinus_count", 3)
                ),
                taxa_with_lyrate_tendency=list(
                    diss_dict.get("taxa_with_lyrate_tendency", [
                        "Packera paupercula",
                        "Packera plattensis",
                        "Packera paupercula var. paupercula",
                        "Packera paupercula var. savannarum",
                    ])
                ),
            )
        else:
            dissection_cfg = DissectionConfig()

        # Parse fold detection parameters with safe fallbacks
        fold_dict = thresh_dict.get("fold_detection")
        if isinstance(fold_dict, dict):
            fold_cfg = FoldDetectionConfig(
                enabled=bool(fold_dict.get("enabled", True)),
                max_folded_aspect_ratio=float(
                    fold_dict.get("max_folded_aspect_ratio", 0.42)
                ),
                max_chord_deviation_ratio=float(
                    fold_dict.get("max_chord_deviation_ratio", 0.035)
                ),
                min_chord_length_px=int(
                    fold_dict.get("min_chord_length_px", 150)
                ),
            )
        else:
            fold_cfg = FoldDetectionConfig()

        # Parse capitulum phenology parameters with safe fallbacks
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

        thresh_kwargs: Dict[str, Any] = {
            k: float(v) for k, v in thresh_dict.items()
            if k in ThresholdsConfig.__dataclass_fields__ and not isinstance(v, dict)
        }
        thresh_kwargs["fold_detection"] = fold_cfg
        thresh_kwargs["dissection"] = dissection_cfg
        thresh_kwargs["capitulum_phenology"] = pheno_cfg
        if "min_solidity" not in thresh_kwargs:
            thresh_kwargs["min_solidity"] = dissection_cfg.default

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
        )

        harv_dict = raw_cfg.get("harvesting", {})
        harvesting = HarvestingConfig(**{
            k: v for k, v in harv_dict.items() if k in HarvestingConfig.__dataclass_fields__
        })

        seg_dict = raw_cfg.get("segmentation", {})
        seg_kwargs = {k: v for k, v in seg_dict.items() if k in SegmentationConfig.__dataclass_fields__}
        if "model_weights" in seg_kwargs:
            seg_kwargs["model_weights"] = _resolve(seg_kwargs["model_weights"])
        segmentation = SegmentationConfig(**seg_kwargs)

        models_dict = raw_cfg.get("models", {})
        models_kwargs: Dict[str, Any] = {}
        if "pcd_weights_path" in models_dict:
            models_kwargs["pcd_weights_path"] = _resolve(models_dict["pcd_weights_path"])
        if "pcd_weights_url" in models_dict:
            models_kwargs["pcd_weights_url"] = str(models_dict["pcd_weights_url"])
        if "pcd_weights_sha256" in models_dict:
            models_kwargs["pcd_weights_sha256"] = str(models_dict["pcd_weights_sha256"])
        models = ModelsConfig(**models_kwargs)

        macro_dict = raw_cfg.get("macroecology", {})
        pheno_dict = macro_dict.get("phenology", {})
        phenology_cfg = PhenologyConfig(
            min_flowering_doy=int(pheno_dict.get("min_flowering_doy", 60)),
            max_flowering_doy=int(pheno_dict.get("max_flowering_doy", 220)),
            robust_regression=bool(pheno_dict.get("robust_regression", True)),
        )
        macroecology = MacroecologyConfig(phenology=phenology_cfg)

        resolved_raw = dict(raw_cfg)
        resolved_raw["resolved_paths"] = {
            k: getattr(paths, k) for k in paths.__dataclass_fields__
        }

        instance = cls(
            paths=paths, taxa=taxa, thresholds=thresholds,
            morphometrics=morph, harvesting=harvesting,
            segmentation=segmentation, models=models,
            macroecology=macroecology, _raw_dict=resolved_raw,
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

# Monographic specialists, herbaria, and type statuses
SPECIALIST_PATTERNS: List[str] = [
    r"\bBarkley\b", r"\bT\.?\s*M\.?\s*Barkley\b", r"\bTheodore\s+M\.?\s+Barkley\b",
    r"\bTrock\b", r"\bD\.?\s*K\.?\s*Trock\b", r"\bDebra\s+K\.?\s+Trock\b",
    r"\bKowal\b", r"\bR\.?\s*R\.?\s*Kowal\b", r"\bRobert\s+R\.?\s+Kowal\b",
    r"\bWeakley\b", r"\bA\.?\s*S\.?\s*Weakley\b", r"\bAlan\s+S\.?\s+Weakley\b",
    r"\bBain\b", r"\bJ\.?\s*F\.?\s*Bain\b", r"\bJohn\s+F\.?\s+Bain\b",
    r"\bMahoney\b", r"\bA\.?\s*M\.?\s*Mahoney\b", r"\bAlison\s+M\.?\s+Mahoney\b",
    r"\bFuller\b", r"\bJ\.?\s*B\.?\s*Fuller\b", r"\bBrandon\s+Fuller\b",
]
MAJOR_HERBARIA_CODES: Set[str] = {
    "NCU", "GA", "US", "NY", "BRIT", "MO", "WIS", "VDB", "FLAS", "TEX", "LL", "TENN", "F"
}
VALID_TYPE_STATUSES: Set[str] = {
    "HOLOTYPE", "ISOTYPE", "LECTOTYPE", "ISOLECTOTYPE", "SYNTYPE", "ISOSYNTYPE",
    "NEOTYPE", "ISONEOTYPE", "PARATYPE", "ISOPARATYPE", "EPITYPE", "TYPE", "COTYPE"
}
EXCLUDED_WESTERN_STATES: Set[str] = {
    "CO", "COLORADO", "NM", "NEW MEXICO", "WY", "WYOMING", "MT", "MONTANA",
    "UT", "UTAH", "AZ", "ARIZONA", "NV", "NEVADA", "ID", "IDAHO",
    "WA", "WASHINGTON", "OR", "OREGON", "CA", "CALIFORNIA", "AK", "ALASKA", "HI", "HAWAII",
}

# Module-level defaults instantiated dynamically from configuration YAML
_DEFAULT_CONFIG = PipelineConfig.from_yaml()
DEFAULT_WORKSPACE = _DEFAULT_CONFIG.paths.workspace_root
DEFAULT_RAW_DIR = _DEFAULT_CONFIG.paths.raw_vouchers_dir
DEFAULT_CURATED_CSV = _DEFAULT_CONFIG.paths.curated_vouchers_csv
DEFAULT_OUTPUT_CSV = _DEFAULT_CONFIG.paths.curated_vouchers_csv
DEFAULT_SUMMARY_LOG = _DEFAULT_CONFIG.paths.summary_log
DEFAULT_QUARANTINE_DIR = _DEFAULT_CONFIG.paths.quarantine_dir
DEFAULT_MIN_MEGAPIXELS = _DEFAULT_CONFIG.thresholds.min_megapixels
DEFAULT_MIN_FILE_SIZE_KB = _DEFAULT_CONFIG.thresholds.min_file_size_kb
DEFAULT_MIN_SHARPNESS_LAPLACIAN = _DEFAULT_CONFIG.thresholds.min_sharpness_laplacian
DEFAULT_TARGET_TAXA = _DEFAULT_CONFIG.taxa.target_species
WESTERN_LONGITUDE_THRESHOLD = _DEFAULT_CONFIG.thresholds.western_longitude_threshold
