# Windows Subsystem for Linux (WSL 2) Setup & Operational Guide

### Running the *Packera dubia* Morphometrics & Species Delimitation Pipeline on Windows
**Author:** J. Brandon Fuller (PhD Candidate, Department of Biology, University of North Carolina at Chapel Hill)  
**Affiliation:** University of North Carolina Herbarium (NCU)  
**Target Environment:** Windows 10 (Build 19044+) or Windows 11 with WSL 2 (Ubuntu 22.04 / 24.04 LTS)

---

## 📋 Overview

The computational pipeline couples automated deep learning morphometrics (LeafMachine2 / PointRend), interactive Segment Anything Model 2 (SAM 2) annotation, high-throughput OpenCV geometric symmetry reconstruction, and R statistical morphometrics (`Momocs`, `MorphoTools2`).

While developed natively on Linux, **Windows Subsystem for Linux (WSL 2)** provides first-class support for the entire pipeline on Windows machines. With modern WSL 2:
- **NVIDIA GPU Acceleration:** Direct pass-through of CUDA cores and Tensor cores to PyTorch and LeafMachine2 via the Windows host display driver.
- **Native GUI Acceleration (WSLg):** Interactive OpenCV/X11 annotation tools like [`annotate_with_sam2.py`](file:///home/brandon/Packera_dubia_morphometrics/scripts/annotation_and_training/annotate_with_sam2.py) render directly in native Windows desktop windows without requiring external X-servers (such as VcXsrv or Xming).
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
| **7** | Configure R Packages | Ubuntu (WSL) | `sudo Rscript -e 'install.packages(...)'` |
| **8** | Launch in VS Code | Ubuntu (WSL) | `code .` (WSL extension) |

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
4. Upon reboot, Ubuntu will launch automatically to finalize initialization. Follow the prompt to create your Linux username and password.

---

## Step 2: NVIDIA CUDA Driver Setup

WSL 2 uses virtual GPU pass-through (WSL-DirectX) provided directly by the Windows NVIDIA display driver.

> [!CAUTION]
> ### 🛑 DO NOT Install Linux NVIDIA Display Drivers Inside WSL!
> A common and destructive error is running `apt-get install nvidia-driver-...` or executing an NVIDIA `.run` driver installer inside the Ubuntu terminal. 
> - **Never** install Linux display drivers inside WSL 2. Doing so overwrites WSL's synthetic Direct3D GPU linkage and **breaks both CUDA pass-through and WSLg desktop graphics**.
> - The Linux environment automatically accesses GPU hardware via `/usr/lib/wsl/drivers` mapped from the Windows host.

### Instructions:
1. **Windows Host:** Download and install the latest Game Ready or Studio Driver for your GPU from the [official NVIDIA Driver Download page](https://www.nvidia.com/Download/index.aspx).
2. **Ubuntu Terminal (WSL):** Verify that CUDA pass-through is active by running:
   ```bash
   nvidia-smi
   ```
   **Expected Output:** A standard NVIDIA System Management Interface table displaying your GPU model (e.g., RTX 3080/4080/A6000), driver version, and CUDA version (e.g., 12.x).

---

## Step 3: Critical Storage Guidance (ext4 vs. `/mnt/c/`)

> [!WARNING]
> ### ⚠️ Always Clone into `~/` (Native Linux ext4) — Avoid `/mnt/c/`!
> The *Packera dubia* pipeline processes **>6,600 high-resolution herbarium vouchers (>25 GB)** and generates tens of thousands of intermediate image patches, binary masks, and contour coordinate files.
> 
> - **The Problem:** Accessing Windows filesystem drives via `/mnt/c/` or mounted cloud drives (e.g., `/mnt/g/My Drive/...`) passes through the WSL 9P network translation protocol. This causes a **10$\times$ to 50$\times$ reduction in disk I/O performance**, turning a 5-minute batch segmentation into hours.
> - **The Solution:** Always clone and execute the repository inside your Ubuntu home directory (`~` or `/home/<username>/`), which resides on a native high-speed ext4 virtual disk (`ext4.vhdx`).

### Clone the Repository:
```bash
# Navigate to your native Linux home directory
cd ~

# Clone the repository
git clone https://github.com/Plantatsrophe/Packera_dubia_morphometrics.git
cd Packera_dubia_morphometrics
```

### Accessing Project Files from Windows File Explorer:
You can seamlessly view or copy files between Windows and WSL at any time:
- In your Ubuntu terminal, run:
  ```bash
  explorer.exe .
  ```
- Or open Windows File Explorer and enter the network path:
  ```text
  \\wsl$\Ubuntu\home\<your_username>\Packera_dubia_morphometrics
  ```

---

## Step 4: Install Required System apt Packages

Update your package repositories and install essential compilers, graphic libraries for OpenCV GUI rendering, R, and geospatial dependencies for environmental rasters (`terra` and `sf`):

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

## Step 5: Python Environments & LeafMachine2 Setup

The pipeline uses two virtual environments to decouple modern PyTorch / SAM 2 dependencies from legacy LeafMachine2 requirements:

### A. Create and Populate Primary Pipeline Environment (`.venv`):
```bash
# From the project root: ~/Packera_dubia_morphometrics
python3 -m venv .venv
source .venv/bin/activate

# Upgrade packaging tools
pip install --upgrade pip setuptools wheel

# Install core pipeline dependencies
pip install -r requirements.txt
```

### B. Configure LeafMachine2 Virtual Environment (`.venv_LM2`):
Execute the automated setup script to provision `.venv_LM2` and configure LeafMachine2:
```bash
bash setup_leafmachine2.sh
```

### C. Run Preflight Diagnostic:
Verify that all system libraries, Python modules, CUDA availability, and R binaries are correctly configured:
```bash
source .venv/bin/activate
python main.py check-env
```

---

## Step 6: R Statistical Computing Environment

Install the required R packages for Elliptic Fourier Analysis (`Momocs`), Gaussian Mixture Models (`mclust`), and Canonical Discriminant Analysis (`MorphoTools2`):

```bash
# Install CRAN dependencies
sudo Rscript -e 'install.packages(c("Momocs", "mclust", "terra", "tidyverse", "sf", "ggplot2", "patchwork", "gridExtra", "optparse", "remotes"), repos="https://cloud.r-project.org/")'

# Install specialized morphometrics and macroecological packages from GitHub
sudo Rscript -e 'remotes::install_github(c("V-Z/MorphoTools2", "danlwarren/ENMTools", "blasbenito/spatialRF"))'
```

---

## Step 7: Running `annotate_with_sam2.py` via Native WSLg Graphics

Windows 11 and updated Windows 10 include **WSLg (Windows Subsystem for Linux GUI)** built-in. WSLg integrates an automated Wayland/X11 display server with GPU acceleration directly into the Windows desktop compositor.

### No Third-Party X-Server Needed!
You do **not** need to install or configure VcXsrv, Xming, or manual `export DISPLAY` IP forwarding. WSLg sets up `$DISPLAY` (typically `:0`) and `$WAYLAND_DISPLAY` automatically.

### Running the Interactive SAM 2 Annotator:
```bash
source .venv/bin/activate

python scripts/annotation_and_training/annotate_with_sam2.py \
    --images-dir data/raw_vouchers/ \
    --output-coco data/annotations/packera_train_coco.json
```

### What to Expect:
- An OpenCV HighGUI window labeled **"Packera dubia - SAM2 Botanical Annotator"** will open as a native Windows application window.
- Smooth mouse wheel zooming ($1\times$ to $16\times$) and middle-click panning operate with hardware-accelerated responsiveness.
- Real-time SAM 2 mask inference will utilize your host NVIDIA GPU via CUDA without latency.
- Full hotkey reference (`0`–`6` class assignment, `k` knife cut, `o` contour toggle, `Tab` multimask cycling) is documented in the main [`README.md`](file:///home/brandon/Packera_dubia_morphometrics/README.md#sam-2-botanical-annotator-hotkey-cheat-sheet).

---

## Step 8: Launching Visual Studio Code with the WSL Extension

Developing inside WSL 2 is completely seamless using Microsoft Visual Studio Code:

1. **Install VS Code on Windows** from [code.visualstudio.com](https://code.visualstudio.com/).
2. In VS Code, open the Extensions view (`Ctrl+Shift+X`) and install the **WSL** extension (Identifier: `ms-vscode-remote.remote-wsl`).
3. Inside your Ubuntu WSL terminal, navigate to the repository directory and type:
   ```bash
   cd ~/Packera_dubia_morphometrics
   code .
   ```
4. VS Code will launch on your Windows desktop, connected directly to your Ubuntu WSL environment:
   - The bottom-left status bar will display `WSL: Ubuntu`.
   - The integrated terminal (`Ctrl+\``) defaults to bash inside your ext4 repository.
   - You can select `.venv/bin/python` as your active Python interpreter for full linting, autocompletion, and debugging.

---

## 🛠️ Troubleshooting & FAQ

### 1. `nvidia-smi` returns `command not found` or fails inside WSL
- **Solution:** Verify that your Windows host has the latest NVIDIA drivers installed. Ensure your WSL instance is running version 2 by opening Windows PowerShell and checking:
  ```powershell
  wsl -l -v
  ```
  If the VERSION column displays `1`, convert it to WSL 2:
  ```powershell
  wsl --set-version Ubuntu 2
  ```

### 2. GUI Window does not appear when running `annotate_with_sam2.py`
- Check that `$DISPLAY` is set:
  ```bash
  echo $DISPLAY
  ```
  It should print `:0`.
- Verify WSLg health by testing a lightweight X11 app:
  ```bash
  sudo apt-get install -y x11-apps
  xclock
  ```
  If `xclock` does not appear, update WSL from PowerShell (`wsl --update`) and restart WSL (`wsl --shutdown`).

### 3. Extremely slow processing or high disk activity
- Ensure that `pwd` shows `/home/...` and **not** `/mnt/c/...`. Always execute from the ext4 filesystem.

---

## 📚 Related Documentation
- Main Operational Pipeline: [`README.md`](file:///home/brandon/Packera_dubia_morphometrics/README.md)
- Complete Pipeline Guide: [`docs/WORKFLOW_GUIDE.md`](file:///home/brandon/Packera_dubia_morphometrics/docs/WORKFLOW_GUIDE.md)
- SAM 2 Botanical Annotation Protocol: [`docs/SAM2_Precision_Botanical_Annotation_Guide.txt`](file:///home/brandon/Packera_dubia_morphometrics/docs/SAM2_Precision_Botanical_Annotation_Guide.txt)
