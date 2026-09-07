# Practical Rules of Thumb for Scientific Pipelines

## 1. Keep Tightly Coupled Logic Together
- If a module's helper functions are only ever called by that module, keep them in the same file (either as internal functions or private methods prefixed with `_`).
- A 600–800 line file containing the complete logic for a phase (e.g., harvesting or mask extraction) is standard and manageable.

## 2. Extract Only When Logic is Reusable or Distinct
- **Extract:** Generic mathematical routines (e.g., coordinate rotation, unit conversions) that are used across multiple unrelated phases.
- **Extract:** Standalone clients (e.g., an independent API wrapper for GBIF if used across multiple separate research projects).
- **Keep Together:** Data-cleaning glue, pipeline orchestration, and format translation specific to this dataset.

## 3. Use Classes and Functions, Not Files, to Manage Length
- Inside a single script, readability is better achieved through clear functions, type annotations, and docstrings rather than breaking the script into multiple files.

## 4. Target Size Guidelines
- **Entry-point scripts / CLI runners:** ~100–300 lines (focused on parsing arguments and calling the core logic).
- **Core pipeline modules:** ~400–800 lines (cohesive units like `harvester.py` or `segmentation.py`).
- **Pure utilities:** ~200–400 lines (e.g., `logger.py`, general geometric transformations).
