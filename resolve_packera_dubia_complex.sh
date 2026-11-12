#!/usr/bin/env bash
# ==============================================================================
# Script: resolve_packera_dubia_complex.sh
# Description:
#   Convenience launcher for the four-stage Packera dubia species-complex
#   discovery, congeneric misidentification purge, purified intraspecific
#   delimitation, and 15 km spatial-thinning sensitivity workflow.
#   Forwards any additional CLI flags (--threads, --gpu-batch-size, --dry-run,
#   --skip-discovery, --summarize-only, --focal-taxon) to
#   scripts/pipeline/resolve_focal_species_complex.py.
# ==============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${PROJECT_ROOT}${PYTHONPATH:+:${PYTHONPATH}}"
cd "${PROJECT_ROOT}"

if [[ -x "${PROJECT_ROOT}/.venv/bin/python" ]]; then
  PYTHON_BIN="${PROJECT_ROOT}/.venv/bin/python"
else
  PYTHON_BIN="python3"
fi

exec "${PYTHON_BIN}" "${PROJECT_ROOT}/scripts/pipeline/resolve_focal_species_complex.py" "$@"
