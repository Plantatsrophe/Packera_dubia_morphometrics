#!/bin/bash
set -e

# 1. Install system dependencies
echo "=== Installing system dependencies ==="
sudo apt-get update && sudo apt-get install -y libgl1 libglib2.0-0 python3-venv

# Ensure current user owns project files if previously run with sudo
sudo chown -R "$(id -u):$(id -g)" LeafMachine2 .venv_LM2 2>/dev/null || true

# 2. Set up dedicated virtual environment named .venv_LM2
echo "=== Setting up virtual environment .venv_LM2 ==="
if [ ! -f ".venv_LM2/bin/activate" ]; then
    if [ -d ".venv_LM2" ]; then
        echo "Removing broken/incomplete .venv_LM2 directory..."
        sudo rm -rf .venv_LM2 2>/dev/null || rm -rf .venv_LM2
    fi
    python3 -m venv .venv_LM2
fi

# 3. Clone LeafMachine2 repository
echo "=== Cloning LeafMachine2 repository ==="
if [ ! -d "LeafMachine2" ]; then
    git clone https://github.com/Gene-Weaver/LeafMachine2.git
else
    echo "LeafMachine2 directory already exists. Skipping clone."
fi

# Ensure user ownership of cloned repo
sudo chown -R "$(id -u):$(id -g)" LeafMachine2 2>/dev/null || true

# 4. Activate virtual environment
echo "=== Activating .venv_LM2 virtual environment ==="
source .venv_LM2/bin/activate

# 5. Install wheel, setuptools (<70 for pkg_resources support in legacy setup.py), Cython and upgrade pip
echo "=== Upgrading pip, setuptools (<70), and wheel ==="
pip install --upgrade pip
pip install "setuptools>=65.5.0,<70.0.0" wheel Cython

# Pre-install packages that require legacy build environment / no build isolation
echo "=== Pre-installing packages requiring legacy setup ==="
pip install --no-build-isolation "visdom>=0.2.4"

# 6. Install dependencies from LeafMachine2/requirements.txt
echo "=== Installing dependencies from LeafMachine2 requirements.txt ==="
if [ -f "LeafMachine2/requirements.txt" ]; then
    # Patch scikit-learn pin for Python 3.12 compatibility (1.1.3 lacks cp312 wheels and relies on deprecated numpy.distutils)
    sed -i 's/scikit-learn==1.1.3/scikit-learn>=1.3.0/g' LeafMachine2/requirements.txt
    pip install -r LeafMachine2/requirements.txt
else
    echo "Error: LeafMachine2/requirements.txt not found!" >&2
    exit 1
fi

# 7. Explicitly install required package versions with NumPy 1.x pinned
echo "=== Installing pycocotools, opencv-contrib-python (<=4.10.0.84), and vit-pytorch (with numpy<2.0.0) ==="
pip install "numpy<2.0.0,>=1.26.0" "pycocotools>=2.0.5" "opencv-contrib-python<=4.10.0.84,>=4.7.0.68" "vit-pytorch==0.37.1"

echo "=== LeafMachine2 setup completed successfully! ==="
