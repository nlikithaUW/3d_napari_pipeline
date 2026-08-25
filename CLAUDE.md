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

**There is now a second, independent handle on it — `pipeline/nuc_side.py`.**
Fibroblasts are abluminal, so which *side of the wall* a nucleus sits on is a
physical statement about its cell type that shares no machinery with a marker
intensity. The signed nucleus-to-surface distance agrees with the marker
classifier on **5 of 6 images** (AUC 0.64–0.77), while the unsigned version is
chance — see `docs/DIRECTION.md` §8C.4 and P05's 2026-08-22 History.

Two cautions that matter more than the result:

- **Agreement between two unvalidated measures is not validation.** It is the
  first external corroboration P05 has ever had and nothing more.
- **The obvious way to validate it does not work, and the reason is
  structural** (§8C.5): calibrating against nuclei that touch the VE-cadherin
  network fails because fibroblasts sit on the outside of that very wall, so the
  reference population is not independent of the axis being measured. What is
  needed is an **independent endothelial label**; `cadherin_junction_proba`
  (P14, currently consumed by nothing) is the nearest candidate in the repo.

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

**`docs/DIRECTION.md` is the program document — read it first, before this
file's module list and before `ISSUES.md`.** It states what is being measured,
why, what is established, and what to do next. This section is the one-paragraph
version and defers to it on every detail.

Relate **per-adhesion tension** to the **local mechanical environment of the
vessel** — specifically **wall shear stress**. The working hypothesis is that
fibroblasts attach preferentially, or adhere with different force, where the
lumen experiences high shear.

Two things follow that are easy to get wrong:

- **Fibroblasts are abluminal.** They were seeded onto an already-formed,
  already-perfused vessel, so they feel no shear directly; the endothelium is
  the sensor and the signal crosses the wall. The shear covariate for an
  adhesion is therefore the shear at the **nearest luminal surface point**, not
  at the adhesion's own position, and transmural distance is a covariate in its
  own right.
- **Geometry is a proxy, not the target.** A geometric descriptor is worth only
  what it predicts about shear. This reordered the plan in 2026-08-22:
  per-vertex surface curvature left the critical path, because Gaussian
  curvature is K = −r″/(r(1+r′²)²) and is therefore **blind to a taper** and
  flags crotch, flare and waist identically. What matters instead is the
  **cross-section** — and width and height enter the shear formula with
  different exponents, so they must stay separate columns and never be
  collapsed into one "caliber".

Reaching any of this credibly needs the depth confound ruled out — depth and
geometry covary on a tube, so optical attenuation masquerades as a geometric
effect. `centroid_z_um` is exported for exactly this check, and the flat
cross-section gives a real handle on it: predicted shear is **U-shaped in z**
while attenuation is monotonic, and those forms are separable.

**Force is not in this dataset.** FRET is uncalibrated and `ratio_clover_ruby`
tracks the ruby noise floor (`docs/ISSUES.md` P02), so the mechanical half
cannot be tested here at all. What is testable now is the **spatial** half —
where the fibroblasts are.

## What exists now

~20 modules. The shape of it:

**Core (`pipeline/`, no Qt):**
- `io.py` — OME-TIFF loading, role-addressable channels (`Stack.by_role`),
  isotropic resampling, `fret_configured` gate.
- `fret.py` — **FROZEN** FRET math (background subtract, bleedthrough, force LUT).
- `gpu.py` — optional cupy shim for windowed ndimage ops. Falls back to CPU on
  genuine OOM only; any other CUDA error raises.
- `denoise.py` — Otsu / Niblack / Sauvola / Bernsen thresholding, size filtering,
  and `build_mask` — the **single dispatch** for threshold + label +
  size-suggest, shared by the Clover tab and `adhesion_rf`'s gate so the mask
  you tune and the mask the RF gates on cannot differ (they each had their own
  copy). **`sauvola_mask` is inert on this data and the reason is structural**
  (`docs/ISSUES.md` P36): it normalises by half the *global* intensity range, so
  the adaptive term contributes ~0.4% of `k` and the threshold falls below
  background — `sauvola AND floor` is bit-identical to `floor` alone. Every
  cadherin mask built before 2026-08-24 was therefore a pure global cut.
  **Niblack is not inert** (`vol > m + k·s`, no global term) and is now the
  default everywhere. Two sign conventions live side by side and run opposite
  ways: Niblack's higher `k` is **stricter**, Sauvola's is **more permissive**.
  Also `local_contrast_score` — the Niblack statistic as a **continuous rectified
  field** rather than a mask (`docs/ISSUES.md` P28). It and `niblack_mask` both
  call `_local_zscore`, so the threshold and the field are the same quantity by
  construction, not two copies that drift; the identity
  `niblack_mask(v, w, k) == (local_contrast_score(v, w, vx, 0) > k)` is pinned by
  `docs/scratch/p28/score_check.py`. **The field was measured on 2026-08-21 to
  make every reconstruction worse** — it is kept so the experiment stays
  re-runnable, not because it is a better input. Do not reach for it as an
  improvement without reading P28's History first.
- `rf.py` — random-forest voxel classifier + instance labeling.
- `cell_mesh.py` — per-cell surfaces, vertex normals, quadric surface curvature.
- `vessel_mesh.py` — vessel surface from the cadherin mask (polar / 2D alpha /
  3D alpha / poisson / alpha-3D-solid→poisson — see `vessel_poisson.py`).
  `rasterize_alpha_solid` turns the alpha complex's *interior* into a volume so
  the last of those can feed alpha's gap-bridging into Poisson's surface fit. **`build_mesh` is the single
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
- `lumen_graph.py` — **the vessel centerline, as a 3D graph on the Poisson
  lumen.** Replaced `centerline.py`, which is deleted: that module skeletonized
  a *max projection* of the cadherin mask, and making a projection solid needs
  `binary_fill_holes`, which is the single root of seven recorded defects
  (`docs/DIRECTION.md` §8B). Input is χ from `vessel_poisson` via
  `build_mesh(..., extras=)`, so **the centerline and the surface are level sets
  of the same field and cannot disagree about where the vessel is.** It gives
  arc length by construction (`docs/ISSUES.md` P27 closes), caliber everywhere
  from a 3D distance transform, spur pruning **as a ratio of the local
  caliber** rather than an absolute length (P23), endpoints classified as
  field-of-view **crop** vs real terminus, and — the payoff — **cross-sections
  cut perpendicular to the axis**, which is the in-plane caliber §8A.2 could not
  measure by any route and which the shear estimate needs. Pure functions, no
  Qt. **It chooses nothing**: `censored`, `section_clipped` and the z-face
  warning are flags for the caller, not filters it applies.
- `nuc_side.py` — which side of the wall a nucleus is on: **signed** distance,
  exact nearest-**triangle** (two passes; the second is what makes it exact),
  no subsample. Only meaningful on a Poisson-family mesh, whose ∇χ normals make
  positive mean abluminal — `mesh_sign_is_luminal` refuses on anything else,
  because on an alpha solid "outward" points *into* the lumen on the inner sheet
  (P21) and the numbers would look identical. Agrees with the P05 marker
  classifier on 5 of 6 images (AUC 0.64–0.77); the unsigned version is chance,
  which is the whole point. **Two unvalidated measures agreeing is not
  validation** — see P05's 2026-08-22 History before quoting it.
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
  `occupancy_grid` block-means the input **values**, so `occ` is mean wall score
  for a continuous field and the occupancy *fraction* exactly for a 0/1 mask —
  bit-identical, which is the contract behind `batch_mesh`'s byte-identical
  rebuild and the reason this module stays off `gpu.*`.
- `mesh_score.py` — how good is a vessel reconstruction: scale stability,
  topology, coverage, support. Pure functions, same doctrine as `saved_check.py`.
  **Four numbers, reported as a vector and never combined into one** — a
  weighting is a hidden constant, and the selection rule belongs to whoever runs
  a sweep. Scoring is not choosing; nothing here picks parameters.
  `coverage` (is the surface everywhere the data is?) and `support` (is it
  anywhere the data is *not*?) are duals and must be read together: each alone
  is gamed by the other's failure mode. Validated against known surfaces, not
  just against this repo's own recorded numbers — see `docs/ISSUES.md` I03.
- `measure_core.py` — headless per-adhesion measurement. `graph_columns` is
  **shared with the GUI worker** rather than copied — P11 records that P04 and
  P12 each had to be written twice, and that is where it stopped. **Never prints**; it
  returns its degradation warnings in `result["notes"]` for the caller to print.
- `saved_check.py` — the consistency checks over `_saved` intermediates, as pure
  functions (numpy in, `list[str]` out; no Qt, no napari, no `state`). Both the
  GUI and the headless path call these, so a warning is one string, not two that
  drift apart. Also `describe_saved`, behind `batch_measure --check-only` and
  `batch_mesh --check-only`.
- `batch.py`, `batch_measure.py`, `batch_mesh.py`, `batch_postprocess.py` —
  headless runners. **`batch.py --steps`** re-runs part of the pipeline —
  `all`, a comma-separated list, or a `first:last` range over
  `dapi, fibroblasts, cadherin, mesh, curvature, graph, clover, rf,
  postprocess`; order is always the pipeline's own, never the order typed. A
  step re-run alone still reads the tab defaults, so it reproduces an earlier
  run only if those have not moved (P26, and P37 for how that bites).
  **`batch.py` writes no CSV** (stops before Measure);
  `batch_measure.py` does. `batch_mesh.py` rebuilds only the vessel geometry
  (mesh / surface curvature / centerline **graph**; `--steps graph` needs a
  Poisson-family method, since the graph is built on the solved lumen) from the saved cadherin mask plus the
  OME *header* — no pixels decoded, no napari — and writes a `vessel_recon.json`
  recording what produced each mesh. **`--cadherin-field score` is the one
  exception to "no pixels decoded"**: it loads the substrate channel to build
  P28's continuous field. It defaults to `mask`, so every existing invocation
  still decodes nothing, and it is **headless-only — the napari tab has no
  equivalent**. That is a deliberate exception to the no-drift rule below, not an
  oversight: the field is measured to degrade the surface, and a slider would
  invite tuning something that should not be tuned yet. If it ever becomes worth
  seeing interactively, the honest form is a read-only score *layer* in the
  Cadherin tab, not a mesh knob. Its fourth step, `--steps score`, builds
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
  it (`saved_check.py`, `vessel_mesh.build_mesh`, `denoise.build_mask`,
  `vessel_mesh.check_alpha_vs_spacing`). Never write the second copy.
  **This rule has been broken three times and each break cost real work**, so
  it is worth stating what the breaks looked like: a headless runner *forcing a
  different method* than the tab (`_run_cadherin_step` forced Sauvola for both
  `batch.py` and "Run all", so the method being tuned was never the method
  shipped — P36); a check that existed **only** in the headless path, so a rule
  was enforced there and silent in the GUI (`alpha3d_um` vs `subsample_um`); and
  a loader **overwriting a tab default on every image load**, so the recorded
  default never reached a build (P37). None of these is a duplicated function —
  they are all one path quietly deciding for the other.
- **A default is a decision, because `batch.py` runs at widget defaults and has
  no CLI for them.** Changing one changes every headless run and invalidates the
  artifacts it produced; nothing may silently overwrite one.
  Headless equivalents of a slider are argparse flags with **no upper bound** —
  a slider ceiling is a widget limit, not a physical one, and hand-tuned values
  sitting exactly at one may have been truncated rather than converged.
- Keep the pipeline narrow. Resist adding abstractions for hypothetical
  future data types — that's exactly what sank the previous version.
- **napari is the only viewer.** A geometry library may be added for its
  algorithms (Open3D's screened Poisson is the live candidate, `docs/ISSUES.md`
  I06), but never for its rendering: `open3d.visualization` and equivalents are
  not to be imported. Anything a library builds comes back as verts/faces and is
  drawn through the same `add_surface` path as every other method. The end user
  gets one window.
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

- **`docs/DIRECTION.md`** — the **program**: goal, hypothesis, physics, dataset
  and the work queue. Cites issue IDs and quotes only the single number that
  carries an argument; **the `Pnn` entry is always authoritative** for the
  measurement and the date. If the two disagree, `ISSUES.md` is right and
  `DIRECTION.md` is stale.
- **`docs/ISSUES.md`** — living state. One entry per problem, permanent IDs
  (`P01`, `P02`, …). Closed entries stay, including ones closed as `invalid`.
- **`docs/audit.md`** — frozen evidence: a dated snapshot of the full audit
  (2026-08-10, commit `5077918`). Never edited to reflect a fix. `ISSUES.md`
  cites it as `audit F1`, `audit F2`, …
- **`README.md`** — public. Carries only a short limitations table linking here.

When a problem is found, opened, or fixed, update `ISSUES.md` — add a `History`
line with the date rather than rewriting the entry in place, so the record of what
was believed and when survives.

**The work queue lives in `docs/DIRECTION.md` §8, and `ISSUES.md`'s
status-at-a-glance section carries the defect-level ordering and its reasoning —
read them there.**
This file deliberately does not restate the queue. It used to, and the copy went
stale within a day of the fixes landing, which is the drift this section warns
about.

Two standing facts that are not in the queue and do not change:

- Segmentation quality and the vessel-surface reconstruction are in reasonable
  shape; the gap is **validation, not features**.
- **Line numbers in `ISSUES.md` go stale on every edit.** They are re-verified
  against the code when entries are touched; if one looks wrong, trust the quoted
  symbol and grep for it.
