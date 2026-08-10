# 3d_napari_pipeline — read this first

Refer to this file at the start of every session in this repo.

## What this project is

A clean-slate rewrite of an older `napari_pipeline` (sibling repo at
`/mnt/d/OneDrive - UW/Documents/GitHub/napari_pipeline`). The old one got
tangled handling multiple data types; this one starts from scratch and
stays narrow.

## Biology / imaging

- Cells cultured on a 3D **vessel / lumen** surface, imaged as z-stacks.
- Voxel size **0.325 × 0.325 × 0.4 µm** (xy × z) for five of the six images;
  the `09step` file is **0.9 µm in z**. Every 3D kernel must scale its z
  component for this anisotropy. Voxel size is read per-file from the OME
  metadata, never hard-coded.
- Condition: `WT` (folder under `3d_images_folder/images/`).

### Two cell types — both DAPI-positive

This is the fact that matters most.

| Cell type | Marker | Notes |
|---|---|---|
| **hiPSC cardiac fibroblasts** | VnTs tension sensor (mClover / mRuby2) | The cells whose adhesion tension we want |
| **Endothelial cells** | VE-cadherin at cell–cell junctions | Form the vessel surface |

**Both stain with DAPI, so a nucleus alone does not tell you which cell it
belongs to.** Separating them is **unsolved** — the current classifier is
unvalidated and its output swings ~2× on one slider step
(`docs/ISSUES.md` P05). Never assume a nucleus count, or anything derived from
one, is fibroblast-only.

- Channels (see `3d_images_folder/params.json`):
  - `mClover` — donor, VnTs sensor (fibroblasts)
  - `mRuby` — acceptor, VnTs sensor (fibroblasts)
  - `FRET` — sensitised-emission channel
  - `VE-Cadherin` — endothelial junctions, **also the only proxy for the vessel
    surface** (no dedicated vessel stain)
  - `Nuclei` — DAPI, both cell types
- The VE-cadherin **rings** are embedded in the curved vessel surface, so the
  **union of rings approximates the vessel surface** — that's how the geometry
  is inferred.

## Force readout — read before touching anything tension-related

- **FRET is supported but NOT calibrated for this dataset.** The `fret` block has
  been removed from `params.json`, and FRET now runs only when that block is
  present *and* the clover/ruby/fret roles are mapped
  (`io.py::fret_configured`). Otherwise it's skipped and `force_*_pN` exports
  blank.
- **The uncalibrated `mClover/mRuby` ratio is the current tension proxy**
  (`ratio_clover_ruby`): Σclover ÷ Σruby over each adhesion, both
  background-subtracted per z-plane. Sum-then-divide, *not* a mean of per-voxel
  ratios — one dim voxel otherwise dominates.
- **A future dataset will have valid FRET, so keep the path.** Do not delete
  FRET code because it's currently unused.
- **`pipeline/fret.py` is FROZEN.** Its math is verbatim from the old pipeline.
  Call it; do not modify it without explicitly flagging and asking.

## End goal (long term)

Relate **per-adhesion tension** to **local vessel geometry** (surface curvature,
in-plane path curvature, distance to the wall and to branch nodes). Reaching
that credibly needs the depth confound ruled out first — depth and curvature
covary on a tube, so optical attenuation masquerades as a curvature effect.
`centroid_z_um` is exported for exactly this check.

## What exists now

~20 modules. The shape of it:

**Core (`pipeline/`, no Qt):**
- `io.py` — OME-TIFF loading, role-addressable channels (`Stack.by_role`),
  isotropic resampling, `fret_configured` gate.
- `fret.py` — **FROZEN** FRET math (background subtract, bleedthrough, force LUT).
- `gpu.py` — optional cupy shim for windowed ndimage ops. Falls back to CPU on
  genuine OOM only; any other CUDA error raises.
- `denoise.py` — Otsu / Niblack / Sauvola / Bernsen thresholding, size filtering.
- `rf.py` — random-forest voxel classifier + instance labeling.
- `cell_mesh.py` — per-cell surfaces, vertex normals, quadric surface curvature.
- `vessel_mesh.py` — vessel surface (alpha shape) from the cadherin mask.
- `centerline.py` — in-plane centerline + path curvature + branch nodes.
- `measure_core.py` — headless per-adhesion measurement.
- `batch.py`, `batch_measure.py`, `batch_postprocess.py` — headless runners.
  **`batch.py` writes no CSV** (stops before Measure); `batch_measure.py` does.
- `view.py` — the napari app.

**GUI (`pipeline/gui/`, one module per tab):** `loader`, `display`, `nuclei`,
`cadherin`, `cadherin_rf`, `cell_mesh`, `vessel_mesh`, `clover`, `adhesion_rf`,
`roi`, `measure`, plus `widgets.py` (shared helpers) and `logpanel.py`.

**Scripts (`Scripts/`):** `nd2_to_ometiff.py`, `gputest.py`, `vecad_mip.py`,
`fibroblast_cellpose.py`.

## Environment

Two machines, two separate environments that happen to share a name:

- **Development:** WSL2, conda env `image_analysis` (see `environment.yml`).
- **Running the GUI:** Windows Anaconda terminal, a *separate* env also called
  `image_analysis`. Installing a package in one does not install it in the other.
- **GUI stack is pinned** to `napari=0.6.0` + `vispy=0.14.3`. vispy 0.15
  breaks 3D rendering (`Cannot SIZE object N` GLIR crash). Keep dev and run
  envs on the same pinned versions.
- **Hardware:** 32/64 c/t CPU, RTX 3090. **The GPU is in use** — `pipeline/gpu.py`
  accelerates the windowed filters and resampling via cupy.
- **A new dependency has to be installed by hand in the Windows env.** Adding it
  to `environment.yml` is not enough; say so explicitly whenever one is added,
  or the GUI will fail there and not in development.

## Working agreement (durable preferences)

- **Do not code until the user explicitly says so.** When given a new
  task, propose the plan, ask clarifying questions, wait for approval.
- Tunable parameters → **magicgui sliders in napari** for live tuning,
  not hard-coded defaults. The user wants to see and adjust.
- Keep the pipeline narrow. Resist adding abstractions for hypothetical
  future data types — that's exactly what sank the previous version.
- FRET math is load-bearing and copied from the old pipeline; treat it
  as frozen unless the user asks otherwise.
- **Don't report a number as meaningful when it isn't.** If a column can't be
  trusted, it belongs in `docs/ISSUES.md` and in the README's "What not to trust
  yet" table, not quietly in the CSV.

## Where the known problems live

**`docs/ISSUES.md` is the single development log** — every known defect, its
status, and the dates it was opened, last verified, and resolved. Read it before
trusting or extending any output. Do not keep a second copy of the priority list
here; it will drift.

- **`docs/ISSUES.md`** — living state. One entry per problem, permanent IDs
  (`P01`, `P02`, …). Closed entries stay, including ones closed as `invalid`.
- **`docs/audit.md`** — frozen evidence: a dated snapshot of the full audit
  (2026-08-10, commit `5077918`). Never edited to reflect a fix. `ISSUES.md`
  cites it as `audit F1`, `audit F2`, …
- **`README.md`** — public. Carries only a short limitations table linking here.

When a problem is found, opened, or fixed, update `ISSUES.md` — add a `History`
line with the date rather than rewriting the entry in place, so the record of what
was believed and when survives.

As of the 2026-08-10 audit the top of the queue is **P01 → P04 → P07 → P03 →
P02**; `ISSUES.md` carries the current ordering and the reasoning. Segmentation
quality and the vessel-surface reconstruction are in reasonable shape; the gap is
validation, not features.
