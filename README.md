# 3d_napari_pipeline

A napari-based image-analysis pipeline for 3D confocal z-stacks of a co-culture
grown on a vessel/lumen surface: **hiPSC cardiac fibroblasts** carrying the VnTs
tension sensor (mClover donor / mRuby2 acceptor) together with **endothelial
cells** marked by VE-cadherin. It turns raw multi-channel stacks into
per-adhesion measurements — a tension readout, size, and geometry (distance to
the vessel wall, to nuclei, and local vessel curvature) — and exports them to
CSV.

> **On the tension number.** The pipeline supports FRET efficiency → force via a
> calibration LUT, but **that path only runs when you configure it** (see
> [Configuring `params.json`](#configuring-paramsjson)). The images shipped with
> this repo have **no valid FRET calibration**, so FRET is skipped and the
> uncalibrated **mClover/mRuby ratio** (`ratio_clover_ruby`) is the tension
> proxy. See [What not to trust yet](#what-not-to-trust-yet).

## What it does

For each image it can:

- **Segment nuclei** (DAPI) with Otsu + connected components, and split touching
  nuclei by watershed.
- **Segment the VE-cadherin** junctional network and build a 3D **vessel mesh**
  (alpha-shape) from it.
- **Measure vessel curvature** two ways: local *surface* curvature (quadric fit
  on the mesh) and *in-plane (X-Y) path* curvature of the vessel centerline
  (which distinguishes straight vs. branched vessels).
- **Segment focal adhesions** with a trained random-forest pixel classifier,
  then split/clean them into instances.
- **Compute FRET efficiency and force** per voxel from the mClover/mRuby2
  channels using the calibration LUT in `Fret_LUT/` — *only when a `fret` block
  is present in `params.json`*.
- **Quantify each adhesion:** the mClover/mRuby ratio (on the native-Z grid),
  FRET force if configured, volume, axis lengths, centroid position, distance to
  the vessel wall, distance to the nearest nucleus, and the local vessel
  curvature — all exported to a timestamped CSV with a settings JSON.

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

## Configuring `params.json`

Image data is not shipped with this repo, so **you bring your own** — and
`params.json` is how you describe it. It lives at
`3d_images_folder/params.json` and has one required key.

```json
{
  "channels": [
    { "name": "FRET",        "role": "fret"      },
    { "name": "mClover",     "role": "clover"    },
    { "name": "mRuby",       "role": "ruby"      },
    { "name": "VE-Cadherin", "role": "substrate" },
    { "name": "Nuclei",      "role": "nuclei"    }
  ]
}
```

### `channels` (required)

A list of `{"name", "role"}` pairs. **`name` must match the channel name inside
the OME-TIFF exactly** — the loader reads names out of the OME XML and looks
each one up; a mismatch raises rather than guessing.
(`Scripts/nd2_to_ometiff.py` writes those names when converting from `.nd2`.)
List only the channels you have; extra channels in the file are ignored.

### The five roles

| Role | Channel | What needs it |
|---|---|---|
| `nuclei` | DAPI | Nuclei tab; `dist_nuc_um` / `nuc_id`; the fibroblast count |
| `substrate` | VE-cadherin | Cadherin tab; vessel mesh; all curvature; `dist_mesh_um` |
| `clover` | mClover (donor) | Adhesion RF segmentation; numerator of the tension ratio |
| `ruby` | mRuby2 (acceptor) | Denominator of the tension ratio |
| `fret` | FRET (sensitised emission) | FRET efficiency / force only |

**Omitting a role** disables what depends on it rather than failing: the tab
prints why it can't run and the affected CSV columns come out blank. `clover`
and `substrate` are the two you effectively can't do without — `clover` drives
adhesion segmentation and `substrate` drives every geometry column.

### Voxel size is *not* configured here

It is read from each image's own OME metadata (`PhysicalSizeZ/Y/X`, parsed by
`io.py::_parse_ome_voxel_zyx`), so files acquired at different z-spacing each
resample onto their own isotropic grid automatically. There is no voxel-size key
in `params.json`; adding one would only create a way for it to disagree with the
image.

### The optional `fret` block

Add this block **only if you have calibrated these constants on your own
instrument with your own sensor:**

```json
"fret": {
  "abt": 0.05338,
  "dbt": 0.3769,
  "G": 1.725,
  "venus_thres": 100,
  "leave_neg": false
}
```

- `abt` / `dbt` — acceptor and donor bleedthrough coefficients.
- `G` — the donor-recovery / sensitised-emission gain factor.
- `venus_thres` — acceptor intensity below which a voxel is excluded.
- `leave_neg` — keep negative efficiencies instead of clipping to `[0, 1]`.

**These values are instrument- and sensor-specific. Do not copy them from this
file or from another lab's.** They encode the spectral bleedthrough of one
particular filter set and detector; borrowed values produce force numbers that
look plausible and mean nothing.

Behaviour is driven by presence, not by a flag: FRET runs only when the `fret`
block exists **and** the `clover`, `ruby`, and `fret` roles are all mapped.
Otherwise it is skipped, the efficiency/force layers are not added, and
`force_*_pN` exports blank. Either way the loader logs which path it took
(`io.py::fret_configured`). Skipping is also substantially faster and uses far
less memory — the FRET volumes are the largest arrays the loader builds.

## Data layout

```
3d_images_folder/
  params.json                # channel name -> role mapping (see above)
  raw_nd2/                   # original .nd2 files (git-ignored)
    <image>.nd2
  images/
    WT/                      # one folder per condition
      <image>.tif            # OME-TIFF z-stacks (ZCYX); channel names and
                             #   pixel size live in the OME metadata
      <image>_saved/         # per-image intermediates + results (auto-created)
  models/                    # trained RF models
Fret_LUT/                    # efficiency -> force calibration table
```

**`images/*/*.tif` is exactly one folder deep.** The loader globs that shape
(`imgs_dir.glob("*/*.tif")`), so a stack sitting directly in `images/`, or
nested two levels down, will not appear in the image picker. Keep the condition
folder.

Convert raw `.nd2` to the expected OME-TIFF with:

```bash
python Scripts/nd2_to_ometiff.py
```

### Filename convention

`V3b_CorrectSettings_40xWater_04step_001.tif` reads as
`<sample>_<settings>_<objective>_<z-step>_<series>`.

⚠️ **`09step` was acquired at 0.9 µm z-spacing**, versus 0.4 µm for the other
five. Nothing breaks — voxel size comes from each file's own metadata and the
stack is resampled accordingly — but it is ~2× coarser in z, so treat its
z-resolved results (adhesion volumes, surface curvature) as less reliable and
don't pool it with the 0.4 µm images without checking.

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
chains the whole pipeline through Measure and the CSV export.

### Headless batch (whole folder, no window)

```bash
python -m pipeline.batch --dir "3d_images_folder/images/WT" --model adhesion_rf.joblib
```

Runs on every OME-TIFF in a folder with the window hidden — faster, unattended,
freeing GPU memory between images. A failing image is reported and skipped so
one bad file can't stall the batch.

⚠️ **The batch stops before Measure and writes no CSV.** Its step list
(`batch.py::_RUN_STEPS`) deliberately omits the Measure + export step, which
crashes headless on a GL/Qt draw. It produces every *intermediate* in each
`<stem>_saved/` folder (proba, instances, mesh, curvature, DAPI labels, cadherin
mask, fibroblast count). To get CSVs, do one of:

```bash
# Headless Measure + export over a folder of already-batched images:
python -m pipeline.batch_measure --dir "3d_images_folder/images/WT"
```

or, per image in the GUI: **Load** → **Load all saved features** → **Measure** →
**Export measures (CSV)**.

## Output columns

One row per adhesion, in each `<stem>_saved/<stem>_adhesion_measures_<ts>.csv`.

| Column | Meaning |
|---|---|
| `image` | Source filename, so concatenated CSVs stay identifiable |
| `label` | Adhesion instance ID within that image |
| `ratio_clover_ruby` | **The tension readout.** Σ(clover) ÷ Σ(ruby) over the adhesion, each background-subtracted per z-plane, on the native-Z grid. Uncalibrated. |
| `dist_mesh_um` | Distance from the adhesion to the nearest vessel-mesh vertex (µm) |
| `centroid_z_um`, `centroid_y_um`, `centroid_x_um` | Adhesion centroid in µm — needed for any spatial check, including the depth confound below |
| `path_curv_xy` | In-plane (X-Y) centerline curvature κ at the nearest centerline point (1/µm) |
| `dist_branch_um` | X-Y distance to the nearest vessel branch node (µm) |
| `curv_mean` | Mean curvature H = (κ1+κ2)/2 of the vessel surface at the nearest surface vertex (1/µm) |
| `curv_k1`, `curv_k2` | Principal curvatures κmax, κmin (1/µm) |
| `curv_gaussian` | Gaussian curvature K = κ1·κ2 (1/µm²) |
| `dist_nuc_um` | Distance to the nearest nucleus (µm) |
| `force_mean_pN`, `force_std_pN` | Mean / SD of the FRET force over the adhesion. **Blank unless FRET is configured.** |
| `force_total_pN` | `force_mean_pN` × volume |
| `volume_um3` | Adhesion volume (µm³) |
| `axis_major_um`, `axis_minor_um` | Fitted ellipsoid axis lengths (µm) |
| `eccentricity` | `axis_major_um / axis_minor_um` |
| `nuc_id` | Label of the nearest nucleus (crude cell assignment) |
| `n_fibroblasts`, `n_nuclei_total` | Image-level counts, repeated on every row |

## What not to trust yet

Several outputs are produced but are **not** yet suitable for interpretation.
`TODO.md` has the full detail and what fixing each involves.

- **`ratio_clover_ruby`** — the tension readout itself carries two caveats,
  both measured on `V3_..._04step_002` (2092 adhesions):
  - **mRuby sits only ~6% above its own background** (30 counts over ~490,
    versus mClover's 664 over ~530). Because `background_subtract` quantises
    the background to ±5 counts, that is a **~16% systematic error on every
    ratio**, and the highest ratios are the adhesions with the least ruby.
    Treat the ratio as a rank, not a calibrated quantity.
  - **It drifts ~23% with imaging depth** (shallowest quartile 25.6, deepest
    19.7). mClover attenuates 24% with depth while mRuby does not, which fully
    accounts for the drift — it is optics, not biology.
- **`force_mean_pN` / `force_std_pN` / `force_total_pN`** — blank on this
  dataset by design (no FRET calibration). If you *do* configure a `fret` block,
  the numbers are only as good as your calibration; the LUT shipped here is for
  a Clover–mRuby2 GGSGGS7 sensor on one specific setup.
- **`n_fibroblasts`** — the fibroblast/endothelial classifier is unvalidated. It
  samples a ~2 µm shell around each nucleus, which is where *both* markers are
  dimmest, and the reported fraction swings by ~2× across one slider step.
- **`curv_mean` / `curv_k1` / `curv_k2`** — the curvature sign inverts on
  branched vessels: normals are oriented "away from the mesh centroid", which a
  Y-junction breaks. `curv_gaussian` and `path_curv_xy` are unaffected (the
  first is a product of two curvatures, so the signs cancel; the second is
  computed in-plane without surface normals).

**Before believing any tension-vs-curvature relationship,** run:

```bash
python Scripts/depth_check.py --dir 3d_images_folder/images/WT --channels
```

The tube lies in xy and is imaged through it, so depth is effectively "near wall
vs far wall" — and curvature covaries with depth. On the one image measured so
far the depth drift is large enough that a curvature relationship would be
substantially an attenuation artefact. Also note that adhesions are not
independent observations — many per cell, few images — so pooled correlation
p-values mostly reflect adhesion count. Report per-image first.

## Helper scripts (`Scripts/`)

- `nd2_to_ometiff.py` — convert raw `.nd2` to OME-TIFF, dropping junk channels
  and writing pixel size + channel names into the OME metadata.
- `gputest.py` — quick check that the cupy/CUDA GPU path works.
- `depth_check.py` — tests whether `ratio_clover_ruby` drifts with imaging depth
  (it does; see above). `--channels` splits the ratio into mClover and mRuby to
  show which channel drives the drift and how much signal each has over
  background.
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
TODO.md              known problems + deferred ideas — read before trusting output
CLAUDE.md            project notes / conventions
```

## Notes

- Image data, per-image `_saved/` outputs, and derived `.npy`/`.csv`/figures are
  git-ignored; only code and small config are tracked.
- Hardware: developed on a 32-core / RTX 3090 (24 GB) workstation; also run on an
  RTX 4070 SUPER (12 GB VRAM). Large-volume GPU work can approach the smaller
  card's limit. On a genuine out-of-memory error the code falls back to CPU for
  that call and prints a note; any other CUDA error is raised rather than
  silently degrading (`pipeline/gpu.py`).
