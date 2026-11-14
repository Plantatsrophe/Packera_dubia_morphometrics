#!/usr/bin/env bash
# ==============================================================================
# Script: resolve_eastern_na_packera.sh
# Description:
#   Convenience launcher for the hierarchical Eastern North American Packera
#   monograph orchestrator (Phase A genus baseline, Phase B per-species complex
#   resolution sweep, and Phase C master synthesis artifacts).
#   Forwards all CLI flags (--prefix, --parent-run-id, --resume,
#   --discovery-run-id, --geographic-scope, --min-vouchers, --only-taxa,
#   --exclude-taxon, --threads, --gpu-batch-size, --dry-run, --summarize-only)
#   to scripts/pipeline/resolve_eastern_na_packera.py.
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

exec "${PYTHON_BIN}" "${PROJECT_ROOT}/scripts/pipeline/resolve_eastern_na_packera.py" "$@"
