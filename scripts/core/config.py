from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


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
    model_weights: Path = PROJECT_ROOT / "models" / "lm2_packera_pcd_finetuned.pth"
    score_thresh: float = 0.40
    min_ucs: float = 0.85
    min_solidity: float = 0.72
    device: str = "cuda"


@dataclass
class PipelineConfig:
    """Centralized, type-hinted configuration for Packera dubia pipeline."""
    paths: PathsConfig
    taxa: TaxaConfig
    thresholds: ThresholdsConfig
    morphometrics: MorphometricsConfig
    harvesting: HarvestingConfig
    segmentation: SegmentationConfig
    _raw_dict: Dict[str, Any] = field(default_factory=dict, repr=False)

    def __getitem__(self, key: str) -> Any:
        return self._raw_dict[key]

    def __contains__(self, key: str) -> bool:
        return key in self._raw_dict

    def get(self, key: str, default: Any = None) -> Any:
        return self._raw_dict.get(key, default)

    def to_dict(self) -> Dict[str, Any]:
        return dict(self._raw_dict)

    @classmethod
    def from_yaml(cls, yaml_path: Optional[Union[str, Path]] = None) -> PipelineConfig:
        """Dynamically parse YAML and resolve paths relative to project root."""
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
        thresholds = ThresholdsConfig(**{
            k: float(v) for k, v in thresh_dict.items() if k in ThresholdsConfig.__dataclass_fields__
        })

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

        resolved_raw = dict(raw_cfg)
        resolved_raw["resolved_paths"] = {
            k: getattr(paths, k) for k in paths.__dataclass_fields__
        }

        return cls(
            paths=paths, taxa=taxa, thresholds=thresholds,
            morphometrics=morph, harvesting=harvesting,
            segmentation=segmentation, _raw_dict=resolved_raw,
        )


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
