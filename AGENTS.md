# Workspace Rules: Packera Dubia Morphometrics

## Configuration Parameters Preservation
All configuration parameters in this workspace (including `LeafMachine2.yaml`, configuration files in `LM2_Project/`, YAML/JSON/INI config files, and testing hyperparameters) are custom values calibrated from model testing.

### Mandatory Directive
- **NEVER** change, overwrite, tune, or adjust configuration parameters unless given explicit and direct instructions from the user to do so.
- Do not modify thresholds, model selection paths, confidence values, batch sizes, or enabled/disabled flags in configuration files during bug fixes or optimizations unless explicitly requested.

## Code Modularity & Script Length Guidelines
### Practical Rules of Thumb for Scientific Pipelines
- **Keep tightly coupled logic together:**
  - If a module's helper functions are only ever called by that module, keep them in the same file (either as internal functions or private methods prefixed with `_`).
  - A 600–800 line file containing the complete logic for a phase (e.g., harvesting or mask extraction) is standard and manageable.
- **Extract only when logic is reusable or distinct:**
  - **Extract:** Generic mathematical routines (e.g., coordinate rotation, unit conversions) that are used across multiple unrelated phases.
  - **Extract:** Standalone clients (e.g., an independent API wrapper for GBIF if used across multiple separate research projects).
  - **Keep Together:** Data-cleaning glue, pipeline orchestration, and format translation specific to this dataset.
- **Use classes and functions, not files, to manage length:**
  - Inside a single script, readability is better achieved through clear functions, type annotations, and docstrings rather than breaking the script into multiple files.
- **Target size guidelines:**
  - **Entry-point scripts / CLI runners:** ~100–300 lines (focused on parsing arguments and calling the core logic).
  - **Core pipeline modules:** ~400–800 lines (cohesive units like `harvester.py` or `segmentation.py`).
  - **Pure utilities:** ~200–400 lines (e.g., `logger.py`, general geometric transformations).

