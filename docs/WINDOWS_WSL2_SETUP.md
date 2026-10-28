# Windows Subsystem for Linux (WSL 2) Setup & Operational Guide

### Running the *Packera dubia* Morphometrics & Species Delimitation Pipeline on Windows
**Author:** J. Brandon Fuller (PhD Candidate, Department of Biology, University of North Carolina at Chapel Hill)  
**Affiliation:** University of North Carolina Herbarium (NCU)  
**Standard Operating Procedure:** `UNC-BOT-SOP-2026-04-REV5`  
**Target Environment:** Windows 10 (Build 19044+) or Windows 11 with WSL 2 (Ubuntu 22.04 / 24.04 LTS)

---

## 📋 Overview

The computational pipeline couples automated deep learning organ segmentation (stock LeafMachine2 v2.2+ `PLANT_LeafPriority` + 38-class Ruler Detector + `leaf_seg_rcnn_v2_1` PointRend), classical dimensionless geometric gatekeeping (`Tier_1_Direct`), R statistical morphometrics (`Momocs`, `MorphoTools2`, `RRPP`, `mclust`, `clustvarsel`, `igraph`), taphonomic indumentum-blind `dinov2_vitb14_reg` foreground patch-token pooling with Cleanlab confident learning, and $0\text{–}15\text{ cm}$ rhizosphere / USGS SGMC bedrock spatial macroecology (`geodata`, `soilDB`, `ranger`, `adespatial`, `blockCV`, `ecospat`).

While developed natively on Linux, **Windows Subsystem for Linux (WSL 2)** provides first-class support for the entire pipeline on Windows machines:
- **NVIDIA GPU Acceleration:** Direct pass-through of CUDA cores and Tensor cores to PyTorch, DINOv2, and LeafMachine2 via the Windows host display driver.
- **Native GUI Acceleration (WSLg):** X11/Wayland graphics render directly in native Windows desktop windows without requiring external X-servers.
- **Near-Native Filesystem Performance:** When operating inside native ext4 virtual disks (`~/`), batch processing on thousands of high-resolution herbarium scans achieves near bare-metal Linux speeds.

---

## ⚡ 10-Minute Setup Checklist

| Step | Action | Platform | Command / Reference |
| :--- | :--- | :--- | :--- |
| **1** | Install WSL 2 & Ubuntu | Windows PowerShell | `wsl --install -d Ubuntu` |
| **2** | Install Host GPU Driver | Windows Host | NVIDIA Windows Display Driver |
| **3** | Verify GPU in Linux | Ubuntu (WSL) | `nvidia-smi` |
| **4** | Clone to Native ext4 | Ubuntu (WSL) | `cd ~ && git clone ...` (Avoid `/mnt/c/`) |
| **5** | Install System Packages | Ubuntu (WSL) | `sudo apt-get install -y libgl1 libglib2.0-0 r-base ...` |
| **6** | Configure Python & LM2 | Ubuntu (WSL) | `source .venv/bin/activate && bash setup_leafmachine2.sh` |
| **7** | Configure R Packages | Ubuntu (WSL) | `Rscript -e 'install.packages(...); remotes::install_github("MomX/Momocs")'` |
| **8** | Launch Pipeline | Ubuntu (WSL) | `./run_pipeline.sh` |

---

## Step 1: Install WSL 2 & Ubuntu

1. Open **PowerShell** or **Windows Terminal** as **Administrator** (Right-click $\rightarrow$ *Run as Administrator*).
2. Execute the official installation command:
   ```powershell
   wsl --install -d Ubuntu
   ```
   > [!NOTE]
   > If WSL is already installed, ensure your WSL kernel and subsystem packages are fully up-to-date:
   > ```powershell
   > wsl --update
   > wsl --set-default-version 2
   > ```
3. **Reboot your PC** if prompted by Windows.
4. Upon reboot, Ubuntu will launch automatically to finalize initialization.

---

## Step 2: NVIDIA CUDA Driver Setup

WSL 2 uses virtual GPU pass-through (WSL-DirectX) provided directly by the Windows NVIDIA display driver.

> [!CAUTION]
> ### 🛑 DO NOT Install Linux NVIDIA Display Drivers Inside WSL!
> - **Never** install Linux display drivers (`apt-get install nvidia-driver-...`) inside WSL 2. Doing so overwrites WSL's synthetic Direct3D GPU linkage.
> - The Linux environment automatically accesses GPU hardware via `/usr/lib/wsl/drivers` mapped from the Windows host.

1. **Windows Host:** Download and install the latest Game Ready or Studio Driver for your GPU from the [official NVIDIA Driver Download page](https://www.nvidia.com/Download/index.aspx).
2. **Ubuntu Terminal (WSL):** Verify that CUDA pass-through is active by running:
   ```bash
   nvidia-smi
   ```

---

## Step 3: Critical Storage Guidance (ext4 vs. `/mnt/c/`)

> [!WARNING]
> ### ⚠️ Always Clone into `~/` (Native Linux ext4) — Avoid `/mnt/c/`!
> Accessing Windows filesystem drives via `/mnt/c/` passes through the WSL 9P network translation protocol, causing a $10\times$ to $50\times$ reduction in disk I/O performance. Always clone and execute the repository inside your Ubuntu home directory (`~` or `/home/<username>/`).

```bash
cd ~
git clone https://github.com/Plantatsrophe/Packera_dubia_morphometrics.git
cd Packera_dubia_morphometrics
```

---

## Step 4: Install Required System apt Packages

```bash
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    curl \
    git \
    libgl1 \
    libglib2.0-0 \
    python3-dev \
    python3-venv \
    python3-pip \
    r-base \
    r-base-dev \
    gdal-bin \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev
```

---

## Step 5: Python Environments & Stock LeafMachine2 v2.2+ Setup

The pipeline uses two virtual environments to decouple modern PyTorch / `dinov2_vitb14_reg` / Cleanlab dependencies (`.venv`) from stock LeafMachine2 v2.2+ and Detectron2 PointRend ABI pins (`.venv_LM2`):

### A. Create and Populate Primary Pipeline Environment (`.venv`):
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
```

### B. Configure LeafMachine2 Virtual Environment (`.venv_LM2`):
```bash
bash setup_leafmachine2.sh
```

### C. Run Preflight Diagnostic:
```bash
source .venv/bin/activate
python main.py check-env
```

---

## Step 6: R Statistical Computing Environment

Install the required R ($\ge 4.3$) packages for Elliptic Fourier Analysis (`Momocs`), Klingenberg allometric testing, Gaussian Mixture Models (`mclust`, `clustvarsel`), Canonical Discriminant Analysis (`MorphoTools2`), spatial macroecology (`terra`, `sf`, `geodata`, `soilDB`, `spatialRF`, `ranger`, `adespatial`, `blockCV`, `ecospat`, `ENMTools`), and unsupervised taxon-complex discovery (`igraph`):

```bash
Rscript -e 'install.packages(c("remotes", "tidyverse", "optparse", "mclust", "clustvarsel", "MorphoTools2", "RRPP", "vegan", "terra", "sf", "geodata", "soilDB", "spatialRF", "ranger", "adespatial", "blockCV", "ecospat", "ENMTools", "igraph", "patchwork", "viridis"), repos="https://cloud.r-project.org"); if (!requireNamespace("Momocs", quietly=TRUE)) remotes::install_github("MomX/Momocs")'
```

---

## Step 7: Running the Automated Pipeline (`run_pipeline.sh`)

The production pipeline runs completely headlessly using stock LeafMachine2 v2.2+ weights (`PLANT_LeafPriority` and `leaf_seg_rcnn_v2_1`) without requiring manual interactive annotation:

```bash
source .venv/bin/activate

# Full complex-wide run
./run_pipeline.sh --run-id packera_wsl_run --skip-harvest

# Run the complete unit test suite
python3 -m unittest discover -s scripts/tests -p "test_*.py" -v
```

> [!NOTE]
> **Optional Standalone SAM 2 Utility (`annotate_with_sam2.py`):**  
> While not used by the automated production pipeline, if you ever wish to launch the optional standalone SAM 2 GUI (`scripts/annotation_and_training/annotate_with_sam2.py`), WSLg renders its OpenCV window natively on the Windows desktop without requiring VcXsrv or Xming.

---

## 📚 Related Documentation
- Main Operational Pipeline: [`README.md`](../README.md)
- Taxonomic Standard Operating Procedure: [`docs/TAXONOMIC_SOP.md`](TAXONOMIC_SOP.md)
- Complete Pipeline Guide: [`docs/WORKFLOW_GUIDE.md`](WORKFLOW_GUIDE.md)
