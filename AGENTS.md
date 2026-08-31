# Workspace Rules: Packera Dubia Morphometrics

## Configuration Parameters Preservation
All configuration parameters in this workspace (including `LeafMachine2.yaml`, configuration files in `LM2_Project/`, YAML/JSON/INI config files, and testing hyperparameters) are custom values calibrated from model testing.

### Mandatory Directive
- **NEVER** change, overwrite, tune, or adjust configuration parameters unless given explicit and direct instructions from the user to do so.
- Do not modify thresholds, model selection paths, confidence values, batch sizes, or enabled/disabled flags in configuration files during bug fixes or optimizations unless explicitly requested.
