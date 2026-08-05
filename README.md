# 3d_napari_pipeline

A napari-based image-analysis pipeline for 3D confocal z-stacks of hiPSC-VnTs Cardaic Fibroblasts cultured on a vessel/lumen surface, imaged with FRET tension (mClover donor / mRuby2 acceptor). It turns raw multi-channel stacks into per-adhesion measurements — tension (Clover/mRuby2 ratio and force), size, and geometry (distance to the vessel wall, to nuclei, and local vessel curvature) — and exports them to CSV.

## What it does

For each image it can:

- **Compute FRET efficiency and force** per voxel from the mClover/mRuby2
  channels, using a calibration lookup table (`Fret_LUT/`).
- **Segment nuclei** (DAPI) with Otsu + connected components, and split touching
  nuclei by watershed.
- **Segment the VE-cadherin** junctional network and build a 3D **vessel mesh**
  (alpha-shape) from it.
- **Measure vessel curvature** two ways: local *surface* curvature (quadric fit
  on the mesh) and *in-plane (X-Y) path* curvature of the vessel centerline
  (which distinguishes straight vs. branched vessels).
- **Segment focal adhesions** with a trained random-forest pixel classifier,
  then split/clean them into instances.
- **Quantify each adhesion:** Clover/mRuby2 ratio (on the native-Z grid), FRET
  force, volume, axis lengths, distance to the vessel wall, distance to the
  nearest nucleus, and the local vessel curvature — all exported to a
  timestamped CSV with a settings JSON for full reproducibility.

## Installation

### 1. Install Miniforge (conda)

If you don't already have `conda`, install
[Miniforge](https://github.com/conda-forge/miniforge) — it provides the `conda`
command on the conda-forge channel.

### 2. Create the `image_analysis` environment

The simplest path is the pinned spec, which installs everything below in one go:

```bash
conda env create -f environment.yml      # macOS: use environment-mac.yml
conda activate image_analysis
```

### What gets installed

- **Core scientific stack** — Python 3.11, numpy, scipy, scikit-image,
  tifffile, pandas, matplotlib, zarr, dask.
- **napari GUI** — `napari=0.6.0` + `vispy=0.14.3` (**pinned** — newer vispy
  breaks 3D rendering, so keep these versions) + PyQt.
- **Segmentation / IO (pip)** — scikit-learn + joblib (the random-forest
  adhesion / VE-cad junction classifiers) and `nd2` (reading raw `.nd2` files).
- **GPU acceleration (pip, optional)** — `cupy-cuda12x` plus
  `nvidia-cuda-runtime-cu12` / `nvidia-cuda-nvrtc-cu12` (the CUDA 12 runtime and
  NVRTC libraries cupy needs to JIT-compile its kernels, especially on Windows).
  If no GPU is present these go unused and the code falls back to CPU
  automatically.

### Building it manually (alternative)

If you'd rather assemble the environment by hand instead of from the YAML:

```bash
conda create -n image_analysis -c conda-forge python=3.11 numpy scipy \
    scikit-image tifffile zarr dask pandas matplotlib napari=0.6.0 \
    vispy=0.14.3 pyqt
conda activate image_analysis
pip install nd2 scikit-learn joblib cupy-cuda12x \
    nvidia-cuda-runtime-cu12 nvidia-cuda-nvrtc-cu12
```

### Optional: Cellpose (for `Scripts/fibroblast_cellpose.py`)

```bash
pip install cellpose
```

For **GPU** Cellpose, make sure PyTorch is a CUDA build (the default pip install
can be CPU-only):

```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121
python -c "import torch; print(torch.cuda.is_available())"   # should print True
```

### Verify the GPU path (optional)

```bash
python Scripts/gputest.py     # prints ... OK if the GPU JIT path works
```

## Data layout

```
3d_images_folder/
  params.json                # channel name -> role mapping (FRET, clover, ruby,
                             #   substrate/VE-cadherin, nuclei)
  images/
    WT/
      <image>.tif            # OME-TIFF z-stacks (ZCYX), pixel size + channel
                             #   names in the OME metadata
      <image>_saved/         # per-image intermediates + results (auto-created)
Fret_LUT/                    # efficiency -> force calibration table
```

Raw microscope `.nd2` files can be converted to the expected OME-TIFF format
with:

```bash
python Scripts/nd2_to_ometiff.py
```

## Running it

### Interactive (napari GUI)

```bash
python -m pipeline.view
```

This opens napari with a tab per stage — **Load, Display, Nuclei, Cadherin,
Cell mesh, Vessel mesh, Clover, Adhesion RF, VE-cad RF, ROI, Measure**. Typical
flow: Load an image → segment nuclei → VE-cadherin mask → vessel mesh +
curvature → (load a trained RF model) → Adhesion RF predict → post-process →
Measure → export CSV. Progress prints appear in the terminal **and** in the
in-app **Log** dock.

Every stage has Save/Load buttons (results land in the image's `_saved/`
folder), plus **Load all saved features** and a one-click **Run all** that
chains the whole pipeline.

### Headless batch (whole folder, no window)

```bash
python -m pipeline.batch --dir "3d_images_folder/images/WT" --model adhesion_rf.joblib
```

Runs the full pipeline on every OME-TIFF in a folder with the window hidden —
faster, unattended, freeing GPU memory between images. Each image writes its
`_saved/` folder with all intermediates, a timestamped CSV, and a settings JSON.
A failing image is reported and skipped so one bad file can't stall the batch.

## Helper scripts (`Scripts/`)

- `nd2_to_ometiff.py` — convert raw `.nd2` to OME-TIFF, dropping junk channels
  and writing pixel size + channel names into the OME metadata.
- `gputest.py` — quick check that the cupy/CUDA GPU path works.
- `vecad_mip.py` — batch max-intensity-projection PNGs of the VE-cadherin channel.
- `fibroblast_cellpose.py` — count endothelial cells from the VE-cadherin channel
  with Cellpose (needs `pip install cellpose`), for deriving fibroblast counts.

## Repository layout

```
pipeline/            core package (I/O, FRET, RF, meshing, curvature, GUI tabs,
                     batch runners)
Scripts/             standalone utility scripts (see above)
3d_images_folder/    data root (images are git-ignored; only config is tracked)
Fret_LUT/            FRET efficiency -> force calibration
environment.yml      conda spec (image_analysis)
CHANGES.md           detailed change log
CLAUDE.md            project notes / conventions
```

## Notes

- Image data, per-image `_saved/` outputs, and derived `.npy`/`.csv`/figures are
  git-ignored; only code and small config are tracked.
- Hardware note: developed/run on an RTX 4070 SUPER (12 GB VRAM). Large-volume
  GPU work can approach that limit; the code degrades to CPU per-call on an
  out-of-memory error rather than crashing.
