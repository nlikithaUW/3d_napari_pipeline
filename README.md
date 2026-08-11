# 3d_napari_pipeline

A napari-based image-analysis pipeline for 3D confocal z-stacks of cells grown on
a **vessel / lumen** surface. It turns raw multi-channel stacks into
**per-adhesion measurements** — a tension readout, size, and local geometry
(distance to the vessel wall, distance to nuclei, vessel curvature) — and exports
them to CSV.

The biological question it is built for: **do the forces that fibroblasts exert
through their focal adhesions depend on the local geometry of the vessel they sit
on?**

> ### ⚠️ Read this before using the numbers
>
> This pipeline is under active development and **several exported columns are
> not yet trustworthy** — including the surface-curvature columns and the tension
> readout itself. The defects are known, measured, and documented; they are not
> hidden in the code.
>
> **[`docs/ISSUES.md`](docs/ISSUES.md)** is the authoritative list. A summary
> table is in [What not to trust yet](#what-not-to-trust-yet) below.
>
> The code is shared because the methods and the failure modes are both useful.
> It is not a validated instrument.

---

## The imaging setup

Cells are cultured on a curved 3D vessel/lumen surface and imaged as z-stacks.
**Two cell types share every field, and both are DAPI-positive:**

| Cell type | Marker | Role |
|---|---|---|
| **hiPSC cardiac fibroblasts** | VnTs tension sensor (mClover donor / mRuby2 acceptor) | The cells whose adhesion tension is measured |
| **Endothelial cells** | VE-cadherin at cell–cell junctions | Form the vessel surface |

Because a nucleus alone does not say which cell it belongs to, **separating the
two cell types is an open problem** in this pipeline
([P05](docs/ISSUES.md#p05)).

There is also **no dedicated vessel stain.** The VE-cadherin junction rings are
embedded in the curved surface, so the union of rings approximates the vessel
surface — that is how all geometry is inferred.

### The tension readout

The pipeline supports FRET efficiency → force via a calibration LUT, but **that
path only runs when you configure it** (see
[Configuring `params.json`](#configuring-paramsjson)). Without a calibrated
`fret` block, FRET is skipped, `force_*_pN` exports blank, and the **uncalibrated
mClover/mRuby ratio** (`ratio_clover_ruby`) is the tension proxy. That proxy has
substantial known problems — see [P02](docs/ISSUES.md#p02).

## What it does

For each image it can:

- **Segment nuclei** (DAPI) with Otsu + connected components, splitting touching
  nuclei by watershed.
- **Segment the VE-cadherin** junctional network and build a 3D **vessel mesh**
  (alpha shape or marching cubes) from it.
- **Measure vessel curvature** two ways: local *surface* curvature (quadric fit on
  the mesh) and *in-plane (X-Y) path* curvature of the vessel centerline, which
  distinguishes straight from branched vessels.
- **Segment focal adhesions** with a trained random-forest pixel classifier on the
  mClover channel, then split and clean them into instances.
- **Compute FRET efficiency and force** per voxel from the mClover/mRuby2 channels
  using the calibration LUT in `Fret_LUT/` — *only when a `fret` block is present
  in `params.json`*.
- **Quantify each adhesion** — the mClover/mRuby ratio, FRET force if configured,
  volume, axis lengths, centroid, distance to the vessel wall, distance to the
  nearest nucleus, and local vessel curvature — exported to a timestamped CSV
  with a companion settings JSON.

## Installation

### 1. Install Miniforge (conda)

If you do not already have `conda`, install
[Miniforge](https://github.com/conda-forge/miniforge).

### 2. Create the environment

```bash
conda env create -f environment.yml      # macOS: use environment-mac.yml
conda activate image_analysis
```

### What gets installed

- **Core scientific stack** — Python 3.11, numpy, scipy, scikit-image, tifffile,
  pandas, matplotlib, zarr, dask.
- **napari GUI** — `napari=0.6.0` + `vispy=0.14.3` (**pinned** — newer vispy
  breaks 3D rendering with a `Cannot SIZE object N` GLIR crash) + PyQt.
- **Segmentation / IO (pip)** — scikit-learn + joblib (the random-forest
  classifiers) and `nd2` (reading raw `.nd2` files).
- **GPU acceleration (pip, optional)** — `cupy-cuda12x` plus
  `nvidia-cuda-runtime-cu12` / `nvidia-cuda-nvrtc-cu12`. With no GPU present these
  go unused and the code falls back to CPU automatically.

### Building it manually (alternative)

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

Image data is not shipped with this repository, so **you bring your own** —
and `params.json` is how you describe it. It lives at
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
the OME-TIFF exactly** — the loader reads names out of the OME XML and looks each
one up; a mismatch raises rather than guessing.
(`Scripts/nd2_to_ometiff.py` writes those names when converting from `.nd2`.)
List only the channels you have; extra channels in the file are ignored.

### The five roles

| Role | Channel | What needs it |
|---|---|---|
| `nuclei` | DAPI | Nuclei tab; `dist_nuc_um` / `nuc_id`; the fibroblast count |
| `substrate` | VE-cadherin | Cadherin tab; vessel mesh; all curvature; `dist_mesh_um` |
| `clover` | mClover (donor) | Adhesion segmentation; numerator of the tension ratio |
| `ruby` | mRuby2 (acceptor) | Denominator of the tension ratio |
| `fret` | FRET (sensitised emission) | FRET efficiency / force only |

**Omitting a role** disables what depends on it rather than failing: the tab
prints why it cannot run and the affected CSV columns come out blank. `clover` and
`substrate` are the two you effectively cannot do without — `clover` drives
adhesion segmentation and `substrate` drives every geometry column.

### Voxel size is *not* configured here

It is read from each image's own OME metadata (`PhysicalSizeZ/Y/X`, parsed by
`io.py::_parse_ome_voxel_zyx`), so files acquired at different z-spacing each
resample onto their own isotropic grid automatically. There is deliberately no
voxel-size key in `params.json`; adding one would only create a way for it to
disagree with the image.

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

> **These values are instrument- and sensor-specific. Do not copy them from this
> file or from another lab's.** They encode the spectral bleedthrough of one
> particular filter set and detector; borrowed values produce force numbers that
> look plausible and mean nothing. This has already happened once in this
> project's history — see [P18](docs/ISSUES.md#p18).

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
(`imgs_dir.glob("*/*.tif")`), so a stack sitting directly in `images/`, or nested
two levels down, will not appear in the image picker. Keep the condition folder.

Convert raw `.nd2` to the expected OME-TIFF with:

```bash
python Scripts/nd2_to_ometiff.py
```

### About `adhesion_rf.joblib`

A trained random-forest adhesion classifier is included at the repository root.
It was **trained on one specific dataset** (one microscope, one sensor, one set of
acquisition settings) and is a starting point, not a general-purpose model.
Expect to retrain it on your own data via the Adhesion RF tab.

## Running it

### Interactive (napari GUI)

```bash
python -m pipeline.view
```

Opens napari with a tab per stage — **Load, Display, Nuclei, Cadherin, Cell mesh,
Vessel mesh, Clover, Adhesion RF, VE-cad RF, ROI, Measure**. Typical flow: load an
image → segment nuclei → VE-cadherin mask → vessel mesh + curvature → load a
trained RF model → Adhesion RF predict → post-process → Measure → export CSV.
Progress prints appear in the terminal **and** in the in-app **Log** dock.

Every stage has Save/Load buttons (results land in the image's `_saved/` folder),
plus **Load all saved features** and a one-click **Run all** that chains the whole
pipeline through Measure and the CSV export.

### Headless batch (whole folder, no window)

```bash
python -m pipeline.batch --dir "3d_images_folder/images/WT" --model adhesion_rf.joblib
```

Runs on every OME-TIFF in a folder with the window hidden — faster, unattended,
freeing GPU memory between images. A failing image is reported and skipped so one
bad file cannot stall the batch.

> **The batch stops before Measure and writes no CSV.** Its step list
> (`batch.py::_RUN_STEPS`) deliberately omits the Measure + export step, which
> crashes headless on a GL/Qt draw. It produces every *intermediate* in each
> `<stem>_saved/` folder (proba, instances, mesh, curvature, DAPI labels, cadherin
> mask, fibroblast count).

To get CSVs, either:

```bash
# Headless Measure + export over a folder of already-batched images:
python -m pipeline.batch_measure --dir "3d_images_folder/images/WT"
```

or, per image in the GUI: **Load** → **Load all saved features** → **Measure** →
**Export measures (CSV)**.

Note that `batch_measure` recomputes adhesion instances from the saved probability
map at standard settings and **rewrites `adhesion_instances.npz`** as a side
effect ([P12](docs/ISSUES.md#p12)).

## Output columns

One row per adhesion, in each `<stem>_saved/<stem>_adhesion_measures_<ts>.csv`.

| Column | Meaning |
|---|---|
| `image` | Source filename, so concatenated CSVs stay identifiable |
| `label` | Adhesion instance ID within that image |
| `ratio_clover_ruby` | **The tension readout.** Σ(clover) ÷ Σ(ruby) over the adhesion, each background-subtracted per z-plane, on the native-Z grid. Uncalibrated |
| `dist_mesh_um` | Distance from the adhesion to the nearest vessel-mesh vertex (µm) |
| `centroid_z_um`, `centroid_y_um`, `centroid_x_um` | Adhesion centroid in µm — needed for any spatial check |
| `path_curv_xy` | In-plane (X-Y) centerline curvature κ at the nearest centerline point (1/µm). Unsigned |
| `dist_branch_um` | X-Y distance to the nearest vessel branch node (µm) |
| `curv_mean` | Mean curvature H = (κ1+κ2)/2 of the vessel surface at the nearest surface vertex (1/µm) |
| `curv_k1`, `curv_k2` | Principal curvatures κmax, κmin (1/µm) |
| `curv_gaussian` | Gaussian curvature K = κ1·κ2 (1/µm²) |
| `dist_nuc_um` | Distance to the nearest nucleus (µm) |
| `force_mean_pN`, `force_std_pN` | Mean / SD of the FRET force over the adhesion. **Blank unless FRET is configured** |
| `force_total_pN` | `force_mean_pN` × volume. Units are pN·µm³, not pN ([P10](docs/ISSUES.md#p10)) |
| `volume_um3` | Adhesion volume (µm³) |
| `axis_major_um`, `axis_minor_um` | Fitted ellipsoid axis lengths (µm) |
| `eccentricity` | `axis_major_um / axis_minor_um` |
| `nuc_id` | Label of the nearest nucleus (crude cell assignment) |
| `n_fibroblasts`, `n_nuclei_total` | Image-level counts, repeated on every row |

## What not to trust yet

Summary only. **[`docs/ISSUES.md`](docs/ISSUES.md) is the authoritative list**,
with measurements, mechanisms, and what fixing each involves.
[`docs/audit.md`](docs/audit.md) is the frozen evidence behind them.

**These describe the CSVs as exported.** Two of the defects below are now fixed
in the code, and one image has been re-exported against the fixes
(`V3_..._04step_002`, `..._adhesion_measures_20260811_105032.csv`). **Every other
CSV predates the fixes.** Those rows say which is which; a fix in the code does
not retroactively fix a file already on disk.

| Column(s) | Problem | ID |
|---|---|---|
| `curv_mean`, `curv_k1`, `curv_k2`, `curv_gaussian` | **Do not use any CSV exported before 2026-08-11.** A destroyed triangle winding, plus an away-from-centroid normal flip that overrode the winding on 44% of surface vertices, made the curvature fitting frame near-random; correlation between the old column and the corrected one is ρ = +0.069 (54.8% sign agreement — a coin flip). **Fixed 2026-08-10**: the columns are now self-consistent (κ2 99.5% consistent sign on the re-export, face-normal coherence 0.21 → 0.97). Re-export any other image before using these columns. **But see the row below — self-consistent is not the same as measuring the vessel** | [P01](docs/ISSUES.md#p01) |
| `curv_mean`, `curv_k1`, `curv_k2`, `curv_gaussian` | **The value depends on the scale you fit it at, with no plateau.** Sweeping the curvature fit radius on the shipped mesh moves the implied vessel radius 14.2 → 134.7 µm as `radius_um` goes 5 → 40, tracking the fit window at ~3×. A real tube gives the same κ at any window well below its radius. Independently, the vessel measures ~86 µm deep and ~140 µm wide, so its radius is of order 40 µm and the shipped κ2 (R = 14.2 µm) is 3–5× too curved to be its circumferential curvature. These columns currently describe the reconstruction's roughness at the chosen scale, not the vessel. Not fixable by tuning the slider | [P24](docs/ISSUES.md#p24) |
| `ratio_clover_ruby` | Driven by the ruby channel's noise floor, not by tension — ρ(ratio, ruby) = −0.773 vs ρ(ratio, clover) = +0.068. Spread exceeds what the sensor can physically produce. Treat as a rank at best | [P02](docs/ISSUES.md#p02) |
| *every column* | Adhesion detection is depth-truncated: the far vessel wall yields 37% of the near wall's detection density per unit surface area | [P03](docs/ISSUES.md#p03) |
| `nuc_id`, `n_fibroblasts`, `n_nuclei_total` | Two different nucleus segmentations; `nuc_id` cannot be joined to the fibroblast classification | [P04](docs/ISSUES.md#p04) |
| `n_fibroblasts` | The fibroblast/endothelial classifier is unvalidated; the reported fraction swings ~2× across one slider step | [P05](docs/ISSUES.md#p05) |
| `ratio_clover_ruby` | Drifts with imaging depth (ρ = −0.187); on this geometry depth is confounded with position along the vessel | [P06](docs/ISSUES.md#p06) |
| `n_nuclei_total`, `n_fibroblasts`, adhesion instance boundaries | `peak_local_max` discarded watershed seeds near the volume edge, so nuclei touching the stack boundary were dropped and elongated adhesions failed to split. **Fixed 2026-08-10**; on the re-exported image the counts move 273 → 276 nuclei and 2,092 → 2,147 adhesions (`n_fibroblasts` happened not to change — all three recovered nuclei were non-fibroblast, which is one image, not a general result). Re-export any other image to pick this up | [P07](docs/ISSUES.md#p07) |
| `path_curv_xy` | **Fitted over a 20 µm window on a vessel ~140 µm wide**, so it largely measures wobble in the extracted centerline rather than bending of the vessel — the skeleton is only ~4.5 vessel-widths long in a field ~2.4 widths across. On top of that: contaminated within `smooth_um` of a branch (filter on `dist_branch_um`, but guard for NaN or non-branched images lose every row); ~10× inflated within 5 µm of the xy field boundary, where 42% of the highest-κ points sit; and unsigned, so inside and outside of a bend are indistinguishable | [P08](docs/ISSUES.md#p08) |
| `dist_mesh_um`, `curv_*` | Measured to a ~2 µm-spaced point cloud including geometry the mesh rejected; curvature is assigned with no distance limit. The point cloud's resolution is set by `subsample_um` (5.0 in the shipped run — a 169:1 reduction), not by the image | [P09](docs/ISSUES.md#p09), [P20](docs/ISSUES.md#p20) |
| `force_*_pN` | Blank by design without a FRET calibration. If you configure one, the numbers are only as good as your calibration | [P18](docs/ISSUES.md#p18) |

### Before drawing any conclusion from this CSV

**As of 2026-08-11 there is no trustworthy geometry variable to relate tension
to.** Both curvature measures are now known to be set by the scale they were
fitted at rather than by the vessel — the surface curvature has no plateau across
fit radii ([P24](docs/ISSUES.md#p24)) and the path curvature uses a window a
fraction of the vessel's own width ([P08](docs/ISSUES.md#p08)). That is the
independent variable in "does tension depend on vessel geometry", so the question
the pipeline is built to answer is currently blocked on the vessel
reconstruction, not on the tension side alone.

Adhesions are also **not independent observations** — many per cell, few images,
and curvature varies slowly in space, so the effective sample size is closer to
the number of cells (or images) than the number of adhesions. A pooled correlation
p-value across all adhesions mostly reports adhesion count.

`docs/ISSUES.md` § **Analysis notes** sets out the design that follows: per-image
first, a spatial null rather than a permutation null, a falsification battery, and
pre-specified exclusions.

```bash
# Test whether the tension readout drifts with imaging depth on your own data:
python Scripts/depth_check.py --dir 3d_images_folder/images/WT --channels
```

## Helper scripts (`Scripts/`)

- `nd2_to_ometiff.py` — convert raw `.nd2` to OME-TIFF, dropping junk channels and
  writing pixel size + channel names into the OME metadata.
- `gputest.py` — quick check that the cupy/CUDA GPU path works.
- `depth_check.py` — tests whether `ratio_clover_ruby` drifts with imaging depth.
  `--channels` splits the ratio into mClover and mRuby to show which channel
  drives the drift and how much signal each has over background.
- `vecad_mip.py` — batch max-intensity-projection PNGs of the VE-cadherin channel.
- `fibroblast_cellpose.py` — count endothelial cells from the VE-cadherin channel
  with Cellpose (needs `pip install cellpose`).

## Repository layout

```
pipeline/            core package (I/O, FRET, RF, meshing, curvature, GUI tabs,
                     batch runners)
pipeline/gui/        one module per napari tab
Scripts/             standalone utility scripts (see above)
docs/
  ISSUES.md          known problems — the living development log. Read this
  audit.md           frozen pipeline audit, 2026-08-10, commit 5077918
3d_images_folder/    data root (images are git-ignored; only config is tracked)
Fret_LUT/            FRET efficiency -> force calibration
environment.yml      conda spec (image_analysis)
CLAUDE.md            project conventions and working notes
```

### Filename convention (this project's data)

`V3b_CorrectSettings_40xWater_04step_001.tif` reads as
`<sample>_<settings>_<objective>_<z-step>_<series>`. The `09step` file was
acquired at 0.9 µm z-spacing versus 0.4 µm for the others — nothing breaks, since
voxel size comes from each file's own metadata, but it is ~2× coarser in z, so
treat its z-resolved results (adhesion volumes, surface curvature) as less
reliable and do not pool it with the 0.4 µm images without checking.

## Notes

- Image data, per-image `_saved/` outputs, and derived `.npy`/`.csv`/figures are
  git-ignored; only code and small config are tracked.
- Developed on a 32-core / RTX 3090 (24 GB) workstation; also run on an RTX 4070
  SUPER (12 GB VRAM). Large-volume GPU work can approach the smaller card's
  limit. On a genuine out-of-memory error the code falls back to CPU for that call
  and prints a note; any other CUDA error is raised rather than silently degrading
  (`pipeline/gpu.py`).

## License

**Not yet licensed** — see [H01](docs/ISSUES.md#h01). Until a `LICENSE` file is
added, default copyright applies and the code is not licensed for reuse. If you
would like to use it, please open an issue.
