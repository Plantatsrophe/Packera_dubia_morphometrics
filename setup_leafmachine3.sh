#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# LeafMachine3 Environment & ONNX Model Setup Script
# Project: Packera dubia Species Delimitation & Morphometrics Pipeline
# Affiliation: University of North Carolina at Chapel Hill Herbarium (NCU)
# ==============================================================================

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LM3_DIR="${PROJECT_ROOT}/LeafMachine3"

echo "=== 1. Checking / Installing uv (0.12.23) ==="
export PATH="${HOME}/.local/bin:${HOME}/.cargo/bin:${PATH}"
if ! command -v uv >/dev/null 2>&1; then
    echo "Installing uv 0.12.23..."
    curl -LsSf https://astral.sh/uv/0.12.23/install.sh | sh
    export PATH="${HOME}/.local/bin:${HOME}/.cargo/bin:${PATH}"
elif ! uv --version | grep -q "0.12.23"; then
    echo "Updating uv ($(uv --version)) to required version 0.12.23..."
    uv self update 0.12.23 || curl -LsSf https://astral.sh/uv/0.12.23/install.sh | sh
    export PATH="${HOME}/.local/bin:${HOME}/.cargo/bin:${PATH}"
else
    echo "uv is already installed: $(uv --version)"
fi

echo "=== 2. Cloning LeafMachine3 repository ==="
if [ ! -d "${LM3_DIR}" ]; then
    git clone https://github.com/Gene-Weaver/LeafMachine3.git "${LM3_DIR}"
else
    echo "LeafMachine3 directory already exists at ${LM3_DIR}. Skipping clone."
fi

echo "=== 3. Detecting hardware acceleration variant ==="
OS_NAME="$(uname -s || echo Linux)"
ARCH_NAME="$(uname -m || echo x86_64)"
LM3_EXTRA="cpu"

if [ "${OS_NAME}" = "Darwin" ] && [ "${ARCH_NAME}" = "arm64" ]; then
    LM3_EXTRA="macos"
elif command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
    LM3_EXTRA="gpu"
fi
echo "Selected LeafMachine3 hardware extra: --extra ${LM3_EXTRA}"

echo "=== 4. Syncing LeafMachine3 environment via uv ==="
pushd "${LM3_DIR}" >/dev/null
uv sync --frozen --extra "${LM3_EXTRA}"

echo "=== 5. Installing pinned Hugging Face ONNX models (~1.7 GB) ==="
uv run --frozen --no-sync lm3 models install --yes

echo "=== 6. Running LeafMachine3 doctor check ==="
DOCTOR_OUT="$(uv run --frozen --no-sync lm3 doctor 2>&1 | tee /dev/stderr)"
if ! printf '%s\n' "${DOCTOR_OUT}" | grep -q "Result: READY"; then
    echo "ERROR: 'lm3 doctor' did not conclude with 'Result: READY'." >&2
    exit 1
fi
popd >/dev/null

echo "=== 6b. Linking canonical LeafMachine3 ONNX models into models/lm3/ and syncing .venv_LM3 ==="
mkdir -p "${PROJECT_ROOT}/models/lm3"
declare -A LM3_STAGE_LINKS=(
    ["lm3_archival_detector__yolo26x_det_1280.onnx"]="${LM3_DIR}/models/archival_detector/model.onnx"
    ["lm3_plant_detector__yolo26x_det_1280.onnx"]="${LM3_DIR}/models/plant_detector/model.onnx"
    ["lm3_landmark_detector__yolo26x_pose_640_canonical.onnx"]="${LM3_DIR}/models/landmark_detector/model.onnx"
    ["lm3_landmark_detector__yolo26x_pose_640.onnx"]="${LM3_DIR}/models/landmark_detector/model.onnx"
    ["lm3_leaf_segmenter__yolo26x_seg_1024.onnx"]="${LM3_DIR}/models/leaf_segmenter/model.onnx"
    ["lm3_specimen_segmenter__unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain.onnx"]="${LM3_DIR}/models/specimen_segmenter/unetpp_effb7_1024_dinov3_vith16plus_2048_selftrain/model.onnx"
    ["lm3_ruler_classifier__dinov3_vits16_tile_four_896_mlp.onnx"]="${LM3_DIR}/models/ruler_classifier/model.onnx"
)
for link_name in "${!LM3_STAGE_LINKS[@]}"; do
    target_file="${LM3_STAGE_LINKS[$link_name]}"
    if [ -f "${target_file}" ]; then
        ln -sfn "${target_file}" "${PROJECT_ROOT}/models/lm3/${link_name}"
    fi
done

if [ -d "${LM3_DIR}/.venv" ]; then
    if [ -L "${PROJECT_ROOT}/.venv_LM3" ]; then
        ln -sfn "${LM3_DIR}/.venv" "${PROJECT_ROOT}/.venv_LM3"
    elif [ -d "${PROJECT_ROOT}/.venv_LM3" ]; then
        if ! "${PROJECT_ROOT}/.venv_LM3/bin/python" -c "import onnxruntime" >/dev/null 2>&1; then
            rm -rf "${PROJECT_ROOT}/.venv_LM3"
            ln -sfn "${LM3_DIR}/.venv" "${PROJECT_ROOT}/.venv_LM3"
        fi
    else
        ln -sfn "${LM3_DIR}/.venv" "${PROJECT_ROOT}/.venv_LM3"
    fi
fi

echo "=== 7. Ensuring main project .venv has ONNX Runtime & core dependencies ==="
if [ ! -d "${PROJECT_ROOT}/.venv" ]; then
    uv venv "${PROJECT_ROOT}/.venv"
fi
if [ "${LM3_EXTRA}" = "gpu" ]; then
    uv pip install --python "${PROJECT_ROOT}/.venv/bin/python" "onnxruntime-gpu>=1.17.0" || \
        uv pip install --python "${PROJECT_ROOT}/.venv/bin/python" "onnxruntime>=1.17.0"
else
    uv pip install --python "${PROJECT_ROOT}/.venv/bin/python" "onnxruntime>=1.17.0"
fi

echo "=== LeafMachine3 setup completed successfully! ==="
