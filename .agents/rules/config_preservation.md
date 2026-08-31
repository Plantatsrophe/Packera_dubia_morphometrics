# Configuration Parameters Preservation Rule

## Policy: Custom Configuration Integrity
All configuration parameters in this workspace (including `LeafMachine2.yaml`, configuration files under `LM2_Project/`, YAML/JSON/INI config files, and testing hyperparameters) are custom values established and calibrated from model testing.

## Mandatory Rule
- **NEVER** change, overwrite, tune, or adjust configuration parameters unless given explicit and direct instructions from the user to do so.
- Do not adjust thresholds, model selection paths, detection/segmentation confidence values, batch sizes, or enabled/disabled component flags in configuration files during bug fixes or optimizations unless explicitly requested.
