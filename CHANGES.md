# Pipeline changes — Windows PC setup + fixes (2026-07-05)

Notes for the original author. This covers the changes made to get the pipeline
running with GPU acceleration on the Windows workstation and to work correctly
through the adhesion segmentation step. It spans one environment fix, two code
changes needed just to load an image on the GPU (§2, §3), and one post-process
behavior fix found while tuning adhesion instances (§4). Each entry gives the
symptom, the root cause, and what was changed. All code edits are minimal and
reversible; the FRET math itself was **not** touched.

Read alongside `CLAUDE.md`.

---

## 0. Machine / environment

The workstation is:

- **GPU:** NVIDIA GeForce RTX 4070 SUPER — **12 GB VRAM**, Ada Lovelace (CUDA 12
  capable, so `cupy-cuda12x` is correct).
- **CPU:** AMD Ryzen 5 5600X (6c/12t) — not the 32/64 c/t in `CLAUDE.md`.
- **RAM:** 32 GB DDR4-3600.
- **OS:** Windows, native Anaconda-style workflow (Miniforge), **not** WSL.

The 12 GB VRAM is the one spec that matters downstream: it makes concurrent GPU
work and large-volume predicts more likely to hit memory limits. Both code
changes below are, in effect, about staying within that budget gracefully.

conda was not previously installed on this PC. Installed **Miniforge**
(`Miniforge3-Windows-x86_64.exe`, silent install) and built the env from the
original `environment.yml` (name `image_analysis`). The pinned GUI stack
(`napari=0.6.0` + `vispy=0.14.3`) was kept unchanged.

---

## 1. `environment.yml` — CUDA runtime + NVRTC wheels (setup, not a bug)

### Symptom
After `conda env create -f environment.yml`, cupy imported and could allocate
device arrays, but the first kernel launch failed. `gputest.py` produced, in
sequence across two attempts:

1. `RuntimeError: Failed to find CUDA headers. Please install CUDA toolkit
   headers (e.g., pip install cupy-cuda12x[ctk]) or specify CUDA_PATH ...`
2. After adding the runtime wheel: `DynamicLibNotFoundError: Failure finding
   "nvrtc*.dll"`.

### Cause
The `cupy-cuda12x` **pip wheel** bundles much of the CUDA runtime but relies on
NVRTC (runtime kernel compilation) for some kernels — including the
`cupyx.scipy.ndimage` filters this pipeline uses. NVRTC needs CUDA headers and
the `nvrtc` shared library, which the wheel does not ship. On a system with no
system-wide CUDA Toolkit (`CUDA_PATH` unset), cupy locates these from the
`nvidia-*-cu12` wheels via `cuda.pathfinder`. They were missing.

This is a known cupy-pip-wheel-on-Windows issue; it does **not** occur when
cupy is installed from conda-forge (which pulls `cuda-cudart-dev` as a
dependency). See CuPy install docs and cupy/cupy#8734.

### Change
Added two packages under the existing `pip:` block, next to `cupy-cuda12x`:

```yaml
      - cupy-cuda12x
      # CUDA runtime + NVRTC libs the cupy-cuda12x wheel needs to JIT-compile
      # kernels on Windows (without these: "Failed to find CUDA headers" then
      # missing nvrtc*.dll). Harmless/unused on a non-GPU machine.
      - nvidia-cuda-runtime-cu12
      - nvidia-cuda-nvrtc-cu12
```

After this, `gputest.py` prints `uniform / median / zoom / OK`. The residual
`UserWarning: CUDA path could not be detected` is cosmetic — cupy loads its
libraries from the wheels, not from a system CUDA install, and the JIT path
(which `gputest.py` exercises) works.

**Action for the repo:** these two lines should be committed to the canonical
`environment.yml` so the env builds cleanly on any GPU box. (Already applied in
this working copy.) They are inert on a non-GPU machine.

---

## 2. `pipeline/fret.py` → `resample_fret_outputs` — the load crash

### Symptom
Loading any image raised a napari worker exception. The user-visible
notification was `cudaErrorAlreadyMapped: resource already mapped`; the full
traceback showed the underlying cause chained beneath it:

```
CUDARuntimeError: cudaErrorMemoryAllocation: out of memory
  ... cupy.cuda.pinned_memory.PinnedMemoryPool.malloc ...
  ... cupy_backends.cuda.api.runtime.hostAlloc ...
During handling of the above exception, another exception occurred:
CUDARuntimeError: cudaErrorAlreadyMapped: resource already mapped
  loader.py:54  -> resample_fret_outputs(...)
  fret.py:138   -> dict(ex.map(_do, keys))          # ThreadPoolExecutor(3)
  io.py:78      -> gpu.zoom(vol, zf, order=1, mode="reflect")
  gpu.py:46     -> cp.asnumpy(gpu_op(cp.asarray(arr), ...))   # H2D transfer
```

Concrete sizes for this file: raw stack `(263, 5, 1024, 1024)` uint16, voxel
`(0.4, 0.325, 0.325)`; isotropic target `0.325` → z-factor `1.23` →
`(324, ...)`. The three FRET output volumes (`fc_vol`, `eff_vol`, `force_vol`)
are each `263×1024×1024` float32 ≈ **1.05 GB**.

### Cause
`resample_fret_outputs` resamples those three volumes through a
`ThreadPoolExecutor(max_workers=3)`. On the GPU path each `_do` calls
`gpu.zoom`, whose first act is `cp.asarray(arr)` — a host→device copy that
allocates a **pinned (page-locked) host transfer buffer**. Three of these run
concurrently, so three ~1 GB pinned buffers are requested at once from cupy's
shared pinned-memory pool. That pool is exhausted (`hostAlloc` →
`cudaErrorMemoryAllocation`), and the threads racing on the shared pool /
context surface it as `cudaErrorAlreadyMapped`.

Notably, `io.to_isotropic` **already anticipates exactly this** and goes
sequential when the GPU is enabled (its comment: "The GPU is already internally
parallel; concurrent host threads would just contend for the device and risk
OOM"). `resample_fret_outputs` simply never received the same guard — it is the
one remaining place that dispatches concurrent GPU work.

### Change
Mirror the `to_isotropic` pattern — run sequentially when the GPU is on, keep
the thread pool only for the CPU path (where it helps, since scipy releases the
GIL):

```python
    from .io import resample_volume
    from . import gpu
    keys = ("fc_vol", "eff_vol", "force_vol")

    def _do(k: str):
        return k, resample_volume(fret_dict[k], src_zyx_um, target_um, order=1)

    if gpu.gpu_enabled():
        # The GPU is already internally parallel; three concurrent host threads
        # would contend for the device and exhaust the pinned-transfer pool
        # (cudaErrorAlreadyMapped / host-alloc OOM). Mirror io.to_isotropic and
        # go sequential when the GPU is on.
        return dict(_do(k) for k in keys)

    with ThreadPoolExecutor(max_workers=3) as ex:
        return dict(ex.map(_do, keys))
```

This is a **threading/orchestration change only**; `compute_fret_volume` and
the FRET/LUT math are untouched, per the `CLAUDE.md` rule that FRET is frozen.

---

## 3. `pipeline/gpu.py` → `_oom` / `_run` — why the CPU fallback didn't catch it

### Symptom
The module docstring promises that "an out-of-memory error on a single call
degrades to CPU for that call rather than crashing." That did not happen — the
error above propagated all the way out and crashed the load worker instead of
silently falling back to scipy.

### Cause
The fallback guard `except _oom()` only caught
`cupy.cuda.memory.OutOfMemoryError`, which is the **memory-pool** OOM raised by
cupy's device allocator. The failure here was a
`cupy.cuda.runtime.CUDARuntimeError` (`cudaErrorMemoryAllocation` from the
pinned **host** allocation, and `cudaErrorAlreadyMapped`) — a completely
different exception class raised at the CUDA-runtime level. It was never in the
`except` tuple, so it escaped the try/except entirely.

### Change
Two edits, both in the fallback path only:

1. Broaden `_oom()` to also include the runtime error class:

```python
def _oom():
    """cupy exception types that mean "the device couldn't do this call" — a
    signal to retry the call on the CPU. Covers both the memory-pool OOM
    (``OutOfMemoryError``) and runtime-level failures such as a pinned host
    alloc OOM or ``cudaErrorAlreadyMapped`` (``CUDARuntimeError``), which are a
    *different* exception class and would otherwise escape the fallback. Empty
    (never-raised) on a CPU-only env."""
    if HAS_GPU:
        return (cp.cuda.memory.OutOfMemoryError,
                cp.cuda.runtime.CUDARuntimeError)
    return ()
```

2. Also drain the pinned-memory pool on fallback (the pool implicated in this
   crash), alongside the existing device-pool free, inside `_run`:

```python
        except _oom():
            if cp is not None:
                cp.get_default_memory_pool().free_all_blocks()
                cp.get_default_pinned_memory_pool().free_all_blocks()
    return cpu_op(arr, **kw)
```

Because every wrapper (`_run`, `uniform_filter`, `local_mean_std`, …) shares the
`except _oom()` clause, broadening `_oom()` makes all of them degrade to CPU on
a runtime-level device failure rather than crashing.

Note: with fix #2 in place, this fallback should no longer be *needed* during a
normal load — but it is the correct safety net for a 12 GB card on large
volumes, which is exactly where a single-call OOM is now most plausible.

---

## 4. `pipeline/rf.py` → `instance_labels_3d` — size filter re-applied after split

### Symptom
During adhesion post-processing, touching adhesions labeled as one giant
connected component (a ~14,700-voxel blob). Enabling `split_merged` correctly
broke it apart (max instance size dropped from 14661 → ~1982, count rose), but
the smallest instance dropped from 101 → **2** voxels: the watershed carved
sub-threshold slivers that `min_size_vox=100` did not remove.

### Cause
In `instance_labels_3d` the order of operations was: connected-components →
remove components smaller than `min_size_vox` → **then** watershed-split large
components. Because the size filter ran *before* the split, any fragment the
watershed created afterward was never size-checked. So with `split_merged` on,
`min_size_vox` silently stopped applying to exactly the pieces most likely to be
noise (the freshly cut slivers).

### Change
Re-apply the same `min_size_vox` filter a second time, *after* the watershed
block and before the final relabel. Guarded by `split_merged` so the non-split
path is byte-for-byte unchanged:

```python
    # Watershed splitting can carve fragments smaller than min_size_vox. The
    # size filter above ran BEFORE the split, so re-apply it here to drop those
    # post-split slivers too — otherwise sub-threshold pieces leak through and
    # min_size_vox silently stops meaning what it says once split_merged is on.
    if split_merged and min_size_vox > 0 and lbl.max() > 0:
        sizes = np.bincount(lbl.ravel())
        keep = (sizes >= min_size_vox)
        keep[0] = True  # background
        remap = np.where(keep, np.arange(len(keep), dtype=np.int32), 0)
        lbl = remap[lbl]
```

Net effect: with `split_merged` on, `min_size_vox` now means what it says — the
split still separates fused adhesions, but the tiny fragments it produces are
dropped by the same size threshold as everything else.

---

## 5. `pipeline/gui/adhesion_rf.py` — save/load the probability map (new)

### Motivation
RF predict is the slow step, and its output ("Adhesion proba") is held only in
`state["adhesion_proba"]`, in memory. Any viewer restart — e.g. to pick up a
code change — discards it and forces a full re-predict just to keep tuning
post-process parameters.

### Change
Added a `proba_path` field plus **Save proba** / **Load proba** buttons in the
Adhesion RF tab, between predict and post-process. Save writes the uint8 proba
volume to a `.npy`; Load reads it back into `state["adhesion_proba"]` and
recreates the "Adhesion proba" layer, so you can restart the viewer and resume
directly at Apply post-process without re-running predict. Load warns if the
saved volume's shape doesn't match the currently loaded stack. Predict and the
RF model itself are unchanged.

---

## 6. `pipeline/gui/adhesion_rf.py` + `loader.py` — clover/ruby ratio in adhesions (new)

### Motivation
A derived quantity was wanted: the clover / ruby intensity ratio, restricted to
the voxels inside segmented adhesions.

### Change
- **`loader.py`** now also keeps the ruby channel in state
  (`state["ruby_raw"]`, mirroring `clover_raw` / `dapi_raw` / `cadh_raw`) and
  resets `state["clover_ruby_ratio"]` on load. Ruby was not previously retained,
  so **an image must be reloaded** for this feature to have data.
- **`adhesion_rf.py`** adds a **"Clover/Ruby ratio (in adhesions)"** button
  (between Apply post-process and model save/load). In a worker it produces:
  - a per-voxel image layer **"Clover/Ruby (adhesions)"** = clover / ruby at
    each voxel inside the "Adhesion instances" mask; 0 outside, and 0 where
    ruby == 0 (the chosen zero policy). Raw intensities, no background
    subtraction. Display contrast auto-set to the 99th percentile of nonzero
    ratios.
  - a **per-adhesion mean** ratio, averaged over each instance's ruby > 0
    voxels, stored in `state["clover_ruby_ratio"]` (`label`, `ratio_mean`,
    `n_vox`) and summarized to the terminal.

  Runs on CPU/numpy (no GPU needed). Guards on missing channels/instances and on
  a shape mismatch between the instances layer and the clover volume. A
  `ratio_path` field with **Save ratio layer** / **Load ratio layer** buttons
  writes the ratio volume to a `.npy` and reloads it as the same layer, so the
  spatial map can be restored without re-running post-process or the ratio.
- **`measure.py`** (Measure tab) computes the same per-adhesion clover/ruby
  ratio internally — aligned to its own adhesion labels, independent of the
  button above — and plots it. The figure grew from 2×2 to 2×3, adding a
  "Clover/Ruby ratio" distribution and a "Clover/Ruby vs mesh distance" hexbin,
  plus a Clover/Ruby line in the printed summary. NaN where an adhesion has no
  ruby > 0 voxel; the panels are simply empty if ruby isn't loaded.

---

## 7. `pipeline/gui/adhesion_rf.py` — ImageJ-style ratio (comparison, new)

### Motivation
The adhesions were previously analysed in ImageJ (`RatioStacksAveIntensity.ijm`):
whole-volume `(mClover / mRuby) × 1000`, kept where the ratio ≥ ~2.0
(`setThreshold(1999.08, ∞)` + `NaN Background`), then an Average-Intensity
Z-projection. That defines its ROI by ratio *magnitude* (no segmentation) and
collapses to 2D — fundamentally different from the adhesion-segmented 3D ratio
in §6. This feature reproduces the ImageJ path so the two can be compared on the
same image.

### Change
Added an **"ImageJ-style ratio + Z-avg"** button plus an `ij_ratio_threshold`
slider (default 2.0) in the Adhesion RF tab. In a worker it computes, over the
*whole* volume (no adhesion mask):
- a 3D layer **"Clover/Ruby ImageJ (3D)"** = `(clover/ruby) × 1000`, NaN wherever
  `ruby == 0` or the ratio is below the threshold (matches ImageJ's ×1000 scale
  and NaN Background);
- a 2D layer **"Clover/Ruby ImageJ (Z-avg)"** = NaN-excluded mean over Z
  (`np.nanmean`), equivalent to ImageJ's Average-Intensity projection.

Displayed with the `inferno` colormap (fire-like) and contrast starting at the
threshold. CPU/numpy; needs only clover + ruby (reload the image so ruby is in
state). This is intentionally separate from §6 so the ratio-threshold ROI and
the adhesion-segmented ROI can be viewed side by side.

---

## 8. Ratio quantification moved to native Z (`loader.py`, `adhesion_rf.py`, `measure.py`)

### Motivation
Verified empirically: recreating the ImageJ ratio on the **native-Z** ND2
reproduces the ImageJ output pixel-for-pixel (correlation 0.9999, identical
means). napari's ratio differed only because it ran on the isotropic-resampled
stack (Z upsampled 0.9 → 0.325 µm, 129 → ~357 linearly-interpolated planes).
Because clover/ruby is nonlinear, interpolated voxels carry ratio values that
were never measured, so the ratio is more faithful on the native grid.

### Change
- **`loader.py`** now keeps the native-Z (pre-resample) clover and ruby channels
  in state (`clover_raw_native`, `ruby_raw_native`) plus `raw_voxel_zyx_um`.
  Segmentation is unchanged — it still runs on the isotropic stack.
- **Adhesion ratio (§6)** now computes on native Z: the isotropic adhesion labels
  are mapped down to the native planes by nearest-neighbour in Z (native plane
  `n` → isotropic index `round(n · raw_z / iso_z)`; XY is unchanged by
  resampling), then the per-voxel layer and per-adhesion means are taken from
  native intensities. The "Clover/Ruby (adhesions)" layer is displayed at the
  native voxel size so it stays physically aligned with the isotropic scene, and
  its save/load reload at native scale.
- **ImageJ-style ratio (§7)** now reads the native channels directly (whole
  volume, no label mapping), so it reproduces ImageJ exactly.
- **Measure tab** per-adhesion ratio uses the same native-Z mapping; shape,
  force, and distance metrics stay on the isotropic grid.

Trade-off: the two native channels add ~`2 × (Z_native · Y · X)` uint16 to memory.

---

## 9. Save / load the DAPI mask and the vessel mesh (`nuclei.py`, `vessel_mesh.py`)

### Nuclei tab
Added a `labels_path` field with **Save DAPI labels** / **Load DAPI labels**
buttons. Saves the full result of "Compute DAPI Otsu + components" — the
component `labels` and `sizes` — to a compressed `.npz`. Load restores
`state["dapi_state"]`, rebuilds the displayed mask, and re-seeds the `min_size`
slider + histogram, so the Otsu+label step never has to be re-run. The labels
are what Measure and the DAPI gate consume.

### Vessel mesh tab
Added a `mesh_path` field with **Save vessel mesh** / **Load vessel mesh**
buttons. Saves the built mesh (verts + faces) to a small `.npz`; Load restores
`state["vessel_mesh_cache"]` and re-adds the "Vessel mesh" surface (the colormap
values are recomputed from the vertices), so the slow mesh build can be skipped.
Also satisfies Measure's mesh-distance / `gate_by_mesh` prerequisite.

---

## 10. VE-Cadherin + instances save/load, and "Load all" (`cadherin.py`, `adhesion_rf.py`, `loader.py`)

### VE-Cadherin mask (Cadherin tab)
Added `cadh_path` + **Save VE-Cadherin mask** / **Load VE-Cadherin mask**. Saves
the cadherin components (`labels`, `sizes`, `centroids`, `normals`, `label_ids`)
to a compressed `.npz`. Load restores `state["cadh_state"]`, re-seeds the
`min_size` slider, and rebuilds the "VE-Cadherin mask" layer — so the vessel mesh
can be **re-tuned** after a reopen, not just reloaded.

### Adhesion instances (Adhesion RF tab)
Added `inst_path` + **Save instances** / **Load instances**. Saves the "Adhesion
instances" label layer to a compressed `.npz` and reloads it, so Apply
post-process can be skipped entirely.

### "Load all saved features" (Load tab)
New button that, after an image is loaded, triggers every tab's loader in one
click — DAPI labels, VE-Cadherin mask, vessel mesh, adhesion proba, adhesion
instances, clover/ruby ratio — each using its own path field. Every tab
registers its loader in shared `state` (`load_dapi_labels`, `load_cadh`,
`load_vessel_mesh`, `load_proba`, `load_instances`, `load_ratio`); the Load tab
calls whichever are present. Never-saved features (missing file) are reported and
skipped.

**Save/load coverage now:** RF model, adhesion proba, adhesion instances,
clover/ruby ratio layer (Adhesion RF); DAPI labels (Nuclei); VE-Cadherin mask
(Cadherin); vessel mesh (Vessel mesh) — with a one-click "Load all" in the Load
tab.

---

## 11. Per-image save/load paths (`loader.py` + tabs)

Every save/load path field now defaults into a per-image folder created on load —
`images/<cond>/<stem>_saved/` next to the image — so features for different
images no longer collide. Each tab registers its path widget in
`state["save_path_widgets"]`; on load, `loader.py` creates the folder, stores it
as `state["save_dir"]`, and retargets each field to `<save_dir>/<default name>`.
The RF model path is left shared (a trained model is usually reused across
images). Net workflow: pick an image → Load → save features (they land in that
image's folder); reopen → Load → **Load all saved features** pulls them back.

---

## 12. In-napari log panel (`pipeline/gui/logpanel.py`, `view.py`)

Added a **"Log"** dock at the bottom of the viewer that mirrors `stdout`/`stderr`,
so all the pipeline's progress `print(...)` messages appear inside napari, not
only in the launching terminal. A thin stream tees each write to the original
stream and emits it as a Qt signal; a GUI slot inserts the text, so prints from
`@thread_worker` background threads are delivered to the GUI thread safely
(queued cross-thread). Includes a Clear button and a 10k-line cap. `view.py`
installs it right after the viewer is created.

---

## 13. Export per-adhesion measurements to CSV (`measure.py`)

Added an `export_path` field + **Export measures (CSV)** button to the Measure
tab. After "Measure adhesions" runs, the full per-adhesion table is written to
CSV via pandas, with the requested columns first — `ratio_clover_ruby` and
`dist_mesh_um` — then `dist_nuc_um`, force mean/std/total, volume, major/minor
axes, eccentricity, and `nuc_id`. The measurement result is cached in
`state["measure_result"]`, and the export path defaults into the per-image
`_saved` folder (§11). CSV opens directly in Excel.

---

## 14. Vessel curvature + Clover/Ruby-vs-curvature (`cell_mesh.py`, `vessel_mesh.py`, `measure.py`)

### `cell_mesh.mesh_curvature(verts, faces, radius_um)`
Per-vertex curvature by a local quadric fit at a *physical* radius: neighbours
within `radius_um` are expressed in a frame whose z-axis is the outward vertex
normal (consistent sign), a quadric `z = a·x² + b·xy + c·y² + d·x + e·y` is
least-squares fit, and the principal curvatures are the eigenvalues of the shape
operator (second fundamental form vs first). Fitting over a physical radius
smooths the triangle-scale segmentation noise. Returns per-vertex `mean` (H),
`k1` (κmax), `k2` (κmin), `gaussian` (K), in 1/µm.

### Vessel mesh tab
**Compute curvature** button + `curvature_radius_um` slider + a measure combobox.
Computes all four curvatures (stored in the mesh cache) and adds a **"Vessel
curvature"** surface overlay recolored by the selected measure (turbo, 2–98th
percentile contrast); switching the combobox recolors without recomputing.

### Measure tab
For each adhesion, the nearest vessel vertex's curvature (from the KD-tree query
already used for mesh distance) is averaged over the adhesion's voxels and added
as `curv_mean`, `curv_k1`, `curv_k2`, `curv_gaussian` — in the CSV and as two new
plot panels (curvature distribution + **Clover/Ruby vs mean curvature**). NaN if
curvature hasn't been computed. For a tube: κ1 ≈ 1/radius (caliber), κ2 ≈ axial
bending, Gaussian ≈ branch/saddle regions — so the CSV lets you test which aspect
the ratio tracks. On each Measure run it also prints **Spearman ρ** of the ratio
vs each curvature measure and vs mesh distance (ρ, p-value, n) to the log, for an
immediate read without opening the CSV.

---

## 15. "Load all" made serial + heavy layers hidden (`loader.py` + tab loaders)

Clicking "Load all saved features" previously fired every loader at once, so
several ~1 GB volumes loaded in parallel and their large 3D layers were added in
a burst — a RAM + GPU-VRAM spike that could hang / near-crash the viewer
(especially on a 12 GB card). Fixed:
- Each tab loader now accepts `after=` (a completion callback) and `visible=`.
  Worker-based loaders call `after()` from their done/error handler; the
  synchronous vessel-mesh loader calls it at the end.
- "Load all" now runs the loaders **one at a time**, chaining via `after` so the
  next starts only after the previous finishes (peak = one load, not six), and
  passes `visible=False` so loaded layers are added hidden — toggle them on in
  the layer list as needed. Individual Load buttons are unchanged (visible, no
  callback).

---

## 16. Fix: hidden-load could break the canvas (`adhesion_rf.py`, `vessel_mesh.py`)

Loading a layer with `visible=False` (from §15's "Load all") — especially a
multiscale **Labels** layer in 3D — could leave it without a GPU visual, so a
later visibility toggle crashed napari's canvas reorder (`KeyError` in
`_reorder_layers` on the "Adhesion instances" layer). Fix: loaders now add each
layer visible (so vispy builds the visual) and then set `.visible`, and the
"Adhesion instances" Labels layer is added **single-scale** (napari's 3D
multiscale-Labels support is unreliable) in both Apply and Load. Recovery for an
already-broken session: delete the "Adhesion instances" layer, or restart.

---

## 17. Fibroblast count — VnTs+ nuclei (`nuclei.py`)

Added a **"Count fibroblasts (VnTs+ nuclei)"** button plus `nucleus_min_dist_um`
and `vnts_threshold` sliders to the Nuclei tab. After "Compute DAPI Otsu +
components", it (1) splits touching nuclei with a distance-transform watershed
(`nucleus_min_dist_um` sets the peak separation), (2) measures each nucleus's
perinuclear mClover in a ~2 µm shell, and (3) marks it a fibroblast if that
signal exceeds a threshold — auto via Otsu on the per-nucleus means, or a manual
`vnts_threshold` (0 = auto). It prints total nuclei / VnTs+ fibroblasts / other
(EC), adds a single-scale "Nuclei class" layer (1 = fibroblast, 2 = other), and
stores the counts in `state["fibroblast_count"]`. The Measure CSV export now also
carries `n_fibroblasts` and `n_nuclei_total` as image-level columns (broadcast to
every adhesion row; NaN if the count wasn't run for that image).

---

## 18. Fix: stale clover mask crashed RF predict (`adhesion_rf.py`, `rf.py`)

With `gate_by_clover` on, if the "Clover filter mask" layer was left over from a
previously-loaded, differently-sized image, its voxel coordinates were out of
bounds for the current volume and `predict_proba` crashed with an `IndexError`
on the output write (`index N out of bounds for axis 1`). Fix: the gate builder
now checks the clover mask shape matches the clover volume and, if not, warns and
predicts the **full volume** instead of gating; `predict_proba` also drops any
out-of-bounds gate indices as a safety net. Recovery for a stale session: reload
the image (or recompute the clover mask) before predicting.

---

## 19. Fix: curvature was fit to interior point-cloud vertices (`cell_mesh.py`, `vessel_mesh.py`, `measure.py`)

The alpha-3D vessel mesh is a *point cloud*: on a real image only ~18% of its
"vertices" lie on the boundary surface (are referenced by a face); the other
~82% are interior cadherin points with no incident face. `mesh_curvature` was
building its KD-tree over **all** vertices, so every local quadric fit was
contaminated by interior points — giving a physically impossible caliber
(median κ1 → 78–145 µm "radius" for ~10 µm vessels) and curvature that couldn't
distinguish branched from straight vessels.

Fix: `mesh_curvature` now restricts the fit to **surface vertices only** (those
in `faces`), reindexing a compact sub-mesh; interior/rank-deficient vertices are
returned as NaN so callers exclude them. Added optional Taubin (λ|μ) smoothing
(`smooth_iter`, GUI slider `curvature_smooth_iter`, default 10) to low-pass the
bumpy alpha-shape without shrinkage. `vessel_mesh` fills NaN before the GL
overlay and prints the implied caliber; `measure` now maps each adhesion to its
nearest **surface** vertex (a KD-tree over finite-curvature vertices) instead of
the nearest of all vertices. After the fix the recovered caliber is a consistent
~9 µm across all CorrectSettings meshes (radius 5 µm, smooth 10) — physically
correct. **Re-run "Compute curvature" and re-export** to refresh the CSV columns.

---

## 20. New: in-plane (X-Y) vessel path curvature (`centerline.py`, `vessel_mesh.py`, `measure.py`)

The biological question ("does adhesion tension differ on straight vs curved vs
branched vessel"), is about how the vessel *path* bends in the top-down (X-Y)
view — not local wall curvature and not anything in Z (the low-resolution,
PSF-elongated axis). New `centerline.vessel_path_curvature_xy` max-projects a 3D
vessel mask to X-Y, solidifies it, skeletonizes to the medial axis, spur-prunes,
and measures the in-plane path curvature κ = 1/R (1/µm) at a physical scale, plus
branch-node locations.

Crucially it is built from the **cadherin mask**, not the vessel mesh: the
alpha-shape mesh is a rough point cloud whose X-Y silhouette skeletonizes into
spurious branches (a straight vessel showed as many "branches" as a branched one),
whereas the cadherin footprint is smooth. Verified on CorrectSettings data:
branched CS_04_002 → 895 centerline points, median κ 0.0059, 6 branch nodes;
straight CS_04_003 → 394 points, median κ 0.0000 (dead straight), 1 node — a
clean separation the surface/mesh curvature never achieved.

GUI: Vessel-mesh tab button "Compute in-plane path curvature (X-Y)" (slider
`path_smooth_um`, default 20) overlays the centerline coloured by κ and stores
`state["vessel_centerline"]`. Measure adds per-adhesion columns `path_curv_xy`
(κ at the nearest centerline point) and `dist_branch_um` (X-Y distance to the
nearest branch node), both in the CSV and the Spearman-vs-ratio printout.

---

## 21. Fix: harden `finite_surface` against GL access-violation crashes (`widgets.py`)

A surface layer draw crashed with `OSError: exception: access violation` inside
vispy's `glDrawArrays` (seen building the pipeline on V3b_04step_003). These
faults happen when vispy reads past a GL buffer, which has three triggers:
a face index outside the vertex array, a non-finite vertex coordinate, or a
per-vertex `values` array that is non-finite or the wrong length. `finite_surface`
previously guarded only non-finite vertex coords. It now also drops faces with
out-of-range indices and makes `values` finite and length-consistent, so no
surface layer (vessel mesh or curvature overlay) can crash the GL draw.
Recovery for a live crash: restart napari to load the fix, then rebuild.

---

## 22. Fix: fibroblast count auto-threshold was Otsu (unstable) → background-relative (`nuclei.py`)

The VnTs+ fibroblast count classified nuclei by Otsu-thresholding their
perinuclear mClover — which gave nonsense (e.g. 98% fibroblast on a straight
channel). Diagnosis on V3_006: perinuclear mClover is only ~1.4× the global
mClover background (median 1009 vs 727) and essentially unimodal, so Otsu (which
assumes two separated modes) is unstable — it landed at 2% on the unsplit nuclei
and ~98% after the watershed split. Tiny population changes flip the result.

Fix: the auto threshold is now **background-relative** — `mClover_bg × factor`
(new slider `vnts_bg_factor`, default 1.5; on V3_006 that gives ~38% at 1.5×,
19% at 2×). A manual absolute threshold still overrides. Added an optional
`exclude VE-cadherin+ (endothelial)` checkbox: when on, a nucleus must also have
low perinuclear VE-cadherin (< cadherin_bg × 1.3) to count as a fibroblast,
which drops endothelial cells (on V3_006: ~14%). The log now reports the
threshold kind, the mClover background, and the resulting fraction.

---

## 23. New: save/load for curvature and fibroblast count (`vessel_mesh.py`, `nuclei.py`, `loader.py`)

Curvature and the fibroblast count were the only per-image results that had to be
recomputed on every reload (mesh save only stores geometry). Added save/load for
all three, matching the existing pattern (per-image default paths, registered in
`save_path_widgets` and the "Load all saved features" chain):

- **Surface curvature** (Vessel-mesh tab) → `vessel_curvature.npz` (the four
  per-vertex arrays); load re-overlays it if the mesh is present.
- **In-plane path curvature** (Vessel-mesh tab) → `vessel_path_curvature.npz`
  (centerline points, κ, branch nodes); load restores `vessel_centerline` and
  re-draws the overlay. The overlay was refactored into a shared `_show_path`.
- **Fibroblast count** (Nuclei tab) → `fibroblast_count.npz` (total / fibroblast
  / other / threshold), so the CSV's `n_fibroblasts` column is restored without
  re-running.

All three are added to the "Load all saved features" sequence in dependency
order (curvature after the mesh; fibroblast count last).

---

## 24. New: "Run all (full pipeline)" one-click orchestration (`loader.py` + every tab)

A "Run all" button (Loader tab, next to "Load all") runs the whole analysis in
order, each tab's compute **and** save, waiting for each threaded step to finish
before starting the next — the same serial `after`-callback chain as "Load all".

Order: DAPI + components → save; fibroblast count → save; VE-cadherin (forces
**Sauvola**) → save; vessel mesh → save; surface curvature → save; in-plane path
curvature → save; clover Otsu; RF predict → save proba; RF post-process (forces
**split_merged**) → save instances; measure → export CSV.

Implementation: every threaded compute (`_run_dapi`, `_count`, `_run_cadherin`,
`_run_build`, `_compute_curv`, `_compute_path`, clover `_compute`, RF `_predict`
/`_apply`, `_measure`) gained an optional `after` callback, called on **both**
success and error so the chain never stalls; the two magicgui computes (DAPI,
cadherin, vessel build) were refactored into plain `_run_*` helpers that read
their widget values. Each tab registers a bundled `state["run_*"]` step
(compute + save). Preconditions/failures in a step (e.g. no trained RF model
loaded) are reported and the chain continues. Preload the RF model and pre-tune
thresholds first; "Run all" uses the current widget values (only method=Sauvola
and split_merged=on are forced).

---

## 25. New: X-Y region-of-interest (ROI) tab (`roi.py`, `measure.py`, `view.py`)

A new "ROI" tab (between Adhesion RF and Measure) lets you draw a rectangle in
the imaging plane and restrict Measure to adhesions inside it — e.g. to exclude
curved vessel segments near the field edges of a "straight" image. "Add ROI box"
drops a default central rectangle (a napari Shapes layer) you can move/resize (or
draw your own); "Apply ROI" stores its X-Y bounds (µm, all Z) in
`state["roi_yx_um"]`; "Clear ROI" restores whole-field measurement. Measure keeps
an adhesion only if its centroid falls inside the box (applied before all
per-adhesion metrics, so ratio/force/curvature/counts all respect the ROI), and
logs how many of N adhesions were kept.

---

## 26. Fix: "Run all" timing bugs — debounced mask layers not ready for the next step (`nuclei.py`, `cadherin.py`, `clover.py`, `measure.py`)

The "Run all" chain runs each step's compute immediately followed by the next,
but several steps updated their **display mask layer via a 150 ms debounce**, so
the downstream step read a stale/empty layer:
- Fibroblast count read the DAPI mask layer → 0 nuclei / 0 fibroblasts. Fix: it
  now builds the mask from the DAPI labels + min_size directly.
- Vessel mesh read the cadherin mask layer (empty) → no mesh built → Measure then
  had no vessel mesh → no result → export skipped. Fix: cadherin now refreshes its
  mask layer **synchronously** after compute (in addition to the debounced call).
- RF predict's clover gate read the clover mask layer. Fix: clover refreshes its
  mask layer synchronously too.
- Measure/export: the chain now exports only if the measurement produced a fresh
  result (clears `measure_result` first), with a clear message otherwise, instead
  of the misleading "run measure first".
- **Critical:** the vessel-mesh refactor read parameters as `controls.orientation`,
  but `orientation` collides with magicgui's `Container.orientation` layout
  property, so it got the string `"vertical"` instead of the widget →
  `AttributeError: 'str' object has no attribute 'value'` on *every* mesh build
  (manual and Run-all). Fixed by reading all params via item access
  (`controls["orientation"].value`); same applied to the cadherin refactor.

---

## 27. Default post-process settings standardized (`adhesion_rf.py`)

Set the adhesion post-process defaults to the validated standard: `prob_threshold`
0.5, `min_size_vox` 100, `split_merged` **ON** (was OFF), with the
`split_min_size_vox` slider shown by default. On 006 this gives adhesion axis
ratios centred in the expected 2-3 range for the cell type (split ON ~2.8 vs
split OFF ~3.0) and a reproducible size distribution. Matches what "Run all"
already forces, so manual and Run-all runs are now consistent. Use the same
settings across all images for comparability.

---

## 28. New: settings logged alongside every CSV export (`measure.py` + all tabs)

Each compute step now records the parameters it used into `state["run_settings"]`
(RF predict gate, adhesion post-process threshold/min_size/split/gates, cadherin
method+params, DAPI/fibroblast settings, vessel-mesh params, surface- and
path-curvature params). On export, a companion `<stem>_<timestamp>_settings.json`
is written next to the timestamped CSV, containing the image name, export
timestamp, voxel scale, whether an ROI was active, the fibroblast count, the
adhesion count, and the full settings dict — so any export can be fully
backtracked to the exact parameters that produced it.

---

## 29. New: headless batch runner (`batch.py`, `loader.py`)

`python -m pipeline.batch --dir <folder> --model <rf.joblib>` runs the full
pipeline on every OME-TIFF in a folder with the napari window hidden
(`Viewer(show=False)`), so a whole batch processes unattended and faster (no
rendering). It reuses the exact GUI code: builds all tabs once (registering the
`run_*` chain), loads the RF model once, then for each image loads it and runs
the same serial `after`-callback chain as "Run all" (DAPI→fibroblast→cadherin→
mesh→curvatures→clover→RF predict→post-process→measure→export), freeing GPU
memory between images. Each image writes its `<stem>_saved` folder with all
intermediates, the timestamped CSV, and the settings JSON — using the standard
defaults (threshold 0.5, min_size 100, split ON, Sauvola, etc.). `_load_image`
gained an `after` hook + error handler so loads can be sequenced. A failed load
or step is reported and skipped so one bad image can't stall the batch.

---

## Verification

- `python -m py_compile pipeline/gui/loader.py pipeline/fret.py pipeline/gpu.py
  pipeline/rf.py pipeline/gui/adhesion_rf.py pipeline/gui/measure.py` → clean.
- `gputest.py` → `OK` (GPU JIT path works).
- Viewer loads an image on the GPU without the crash (previously failed at
  `resample_fret_outputs`).
- Adhesion RF predict + post-process run on the GPU; `split_merged` separates
  fused adhesions and the post-split size filter (§4) removes the resulting
  slivers.

## Open questions / recommendations for the author

1. **FRET outputs are placeholders here.** The `params.json` FRET constants are
   currently all 1.0 and the Force/efficiency layers are unused. The crash
   occurred while resampling those unused volumes. Worth considering:
   skip `compute_fret_volume` / `resample_fret_outputs` entirely when the FRET
   constants aren't set, instead of computing and resampling meaningless data
   on every load. Fix #2 makes it safe; it doesn't make it useful.
2. **12 GB VRAM headroom.** Concurrency was the trigger this time, but a single
   isotropic resample or the adhesion predict on a large stack could still
   approach the limit. The fallback (fix #3) now covers it, but if predict is
   frequently falling back to CPU it may be worth tiling the volume on the GPU
   path instead.
3. The two `nvidia-*-cu12` pip deps (§1) should land in the canonical
   `environment.yml`.
