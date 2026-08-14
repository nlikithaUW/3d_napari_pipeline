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
unvalidated, its output swings ~2× on one slider step, and across the six images
the fibroblast fraction spans 5.0%–14.9% (`docs/ISSUES.md` P05). Never assume a
nucleus count, or anything derived from one, is fibroblast-only.

Since 2026-08-11 that classification is joinable per adhesion via the
`nuc_is_fibroblast` column (`docs/ISSUES.md` P04). **Joinable is not validated** —
the column inherits every P05 problem above, now attached to every row rather
than to two image-level scalars. Do not stratify an analysis on it yet.

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
- `vessel_mesh.py` — vessel surface from the cadherin mask (polar / 2D alpha /
  3D alpha / poisson — see `vessel_poisson.py`). **`build_mesh` is the single
  dispatch point** and both the napari tab and `batch_mesh` call it, so the
  interactive and headless paths cannot drift — the same doctrine as
  `saved_check.py`. `LIVE_PARAMS` says which knobs each method actually reads;
  the rest are inert and must be recorded as such. It returns a 4-tuple and must
  keep doing so (`saved_check`, `batch_mesh` and the GUI all destructure it);
  anything extra a method produces goes into the optional `extras=` out-dict.
  **Adding a method means editing four places in one change** — `METHODS` +
  `LIVE_PARAMS` + `ALL_PARAMS` here, and `MESH_DEFAULTS` + the `explicit` dict +
  an argparse flag in `batch_mesh` — because `resolve_params` indexes both by
  `ALL_PARAMS` and raises `KeyError` otherwise.
- `centerline.py` — in-plane centerline + path curvature + branch nodes, and
  `lumen_center_z` (the lumen-centre Z under each point, from a 3D
  distance-from-wall field, plus the **inscribed lumen radius** — the only direct
  caliber measurement in the pipeline). The footprint is replicate-padded before
  `skeletonize` so the field-of-view border does not fork the medial axis.
  **The points it returns are a point set in raster order, not a traced path**
  (`docs/ISSUES.md` P27) — there is no arc length, and κ survives only because it
  is fitted over KD-tree neighbourhoods.
- `vessel_poisson.py` — the fourth reconstruction method, and the only one that
  **approximates** rather than interpolates: sheet normals from the structure
  tensor, oriented into the lumen *volume* (no centerline, no fill-holes, no
  skeletonize), then a Neumann DCT Poisson solve for a smooth indicator function
  whose level set is the surface. Pure functions, and deliberately **no `gpu.*`**
  — cupy and scipy differ in the last bits and that would break `batch_mesh`'s
  byte-identical rebuild. **The surface it returns is open at the x/y domain
  faces by design** (that is the Neumann boundary; the vessel really does leave
  through them), so it carries its own per-vertex normals and
  `cell_mesh.mesh_curvature` must be handed them — the enclosed-volume sign rule
  cannot work on an open surface. It is the first reconstruction in the project
  to pass `docs/ISSUES.md` P24's plateau test (`slope_loglog` +0.009 against
  alpha 3D's −1.428 on the same mask); it does **not** beat alpha 3D on coverage,
  and I06 records the whole four-vector rather than a verdict.
- `mesh_score.py` — how good is a vessel reconstruction: scale stability,
  topology, coverage, support. Pure functions, same doctrine as `saved_check.py`.
  **Four numbers, reported as a vector and never combined into one** — a
  weighting is a hidden constant, and the selection rule belongs to whoever runs
  a sweep. Scoring is not choosing; nothing here picks parameters.
  `coverage` (is the surface everywhere the data is?) and `support` (is it
  anywhere the data is *not*?) are duals and must be read together: each alone
  is gamed by the other's failure mode. Validated against known surfaces, not
  just against this repo's own recorded numbers — see `docs/ISSUES.md` I03.
- `measure_core.py` — headless per-adhesion measurement. **Never prints**; it
  returns its degradation warnings in `result["notes"]` for the caller to print.
- `saved_check.py` — the consistency checks over `_saved` intermediates, as pure
  functions (numpy in, `list[str]` out; no Qt, no napari, no `state`). Both the
  GUI and the headless path call these, so a warning is one string, not two that
  drift apart. Also `describe_saved`, behind `batch_measure --check-only` and
  `batch_mesh --check-only`.
- `batch.py`, `batch_measure.py`, `batch_mesh.py`, `batch_postprocess.py` —
  headless runners. **`batch.py` writes no CSV** (stops before Measure);
  `batch_measure.py` does. `batch_mesh.py` rebuilds only the vessel geometry
  (mesh / surface curvature / centerline) from the saved cadherin mask plus the
  OME *header* — no pixels decoded, no napari — and writes a `vessel_recon.json`
  recording what produced each mesh. Its fourth step, `--steps score`, builds
  nothing and only records `mesh_score`'s four numbers into that JSON; it is
  opt-in because it costs more than the rebuild does. All print a per-image
  summary and set an exit code.
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
- **The GUI and the headless path must both work and must not drift.** When both
  need the same computation, extract it into the core package and have both call
  it (`saved_check.py`, `vessel_mesh.build_mesh`). Never write the second copy.
  Headless equivalents of a slider are argparse flags with **no upper bound** —
  a slider ceiling is a widget limit, not a physical one, and hand-tuned values
  sitting exactly at one may have been truncated rather than converged.
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

**`ISSUES.md` carries the current ordering and the reasoning — read it there.**
This file deliberately does not restate the queue. It used to, and the copy went
stale within a day of the fixes landing, which is the drift this section warns
about.

Two standing facts that are not in the queue and do not change:

- Segmentation quality and the vessel-surface reconstruction are in reasonable
  shape; the gap is **validation, not features**.
- **Line numbers in `ISSUES.md` go stale on every edit.** They are re-verified
  against the code when entries are touched; if one looks wrong, trust the quoted
  symbol and grep for it.
