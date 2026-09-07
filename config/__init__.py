"""
Centralized Configuration Loader for Packera dubia Morphometrics Pipeline.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional, Union

import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"


def load_config(config_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """
    Load pipeline configuration from a YAML file, resolving relative paths
    against the project root directory.

    Args:
        config_path: Optional path to YAML configuration file. Defaults to config/config.yaml.

    Returns:
        Dict[str, Any]: Parsed configuration dictionary with resolved paths.
    """
    target_path = Path(config_path) if config_path else DEFAULT_CONFIG_PATH

    if not target_path.exists():
        raise FileNotFoundError(f"Configuration file not found at: {target_path}")

    with open(target_path, "r", encoding="utf-8") as f:
        cfg: Dict[str, Any] = yaml.safe_load(f) or {}

    # Resolve paths relative to PROJECT_ROOT
    if "paths" in cfg and isinstance(cfg["paths"], dict):
        resolved_paths: Dict[str, Path] = {}
        for key, val in cfg["paths"].items():
            p = Path(val)
            if not p.is_absolute():
                p = (PROJECT_ROOT / p).resolve()
            resolved_paths[key] = p
        cfg["resolved_paths"] = resolved_paths

    return cfg


__all__ = ["load_config", "PROJECT_ROOT", "DEFAULT_CONFIG_PATH"]
