"""
Centralized Configuration Loader for Packera dubia Morphometrics Pipeline.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Union

from scripts.core.config import DEFAULT_CONFIG_PATH, PROJECT_ROOT, ModelsConfig, PipelineConfig


def load_config(config_path: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """
    Load pipeline configuration from a YAML file, resolving relative paths
    against the project root directory.

    Args:
        config_path: Optional path to YAML configuration file. Defaults to config/config.yaml.

    Returns:
        Dict[str, Any]: Parsed configuration dictionary with resolved paths.
    """
    cfg_obj = PipelineConfig.from_yaml(config_path)
    return cfg_obj.to_dict()


__all__ = ["PipelineConfig", "ModelsConfig", "load_config", "PROJECT_ROOT", "DEFAULT_CONFIG_PATH"]
