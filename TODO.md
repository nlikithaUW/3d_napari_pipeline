# TODO — known problems and deferred ideas

Two sections, deliberately separated: **Known problems** are defects in numbers
the pipeline already exports, so they can mislead. **Deferred ideas** are things
that don't exist yet, which cannot.

---

# Known problems

## 1. Curvature sign inverts on branched vessels

**Where:** `pipeline/cell_mesh.py::vertex_outward_normals`, consumed by
`mesh_curvature`.

**What's wrong:** vertex normals are oriented by pointing them *away from the
mesh centroid*. That is a correct heuristic for a single convex-ish blob. These
vessels are Y-shaped: for a branch that sits on the far side of the junction,
"away from the centroid" points **into** the lumen, so the normal flips and the
signed curvature flips with it.

**What you see:** `curv_mean`, `curv_k1`, `curv_k2` change sign across a
junction on what is geometrically the same kind of surface. Pooling them
averages `+κ` against `−κ` toward zero and destroys any real trend.

**Not affected:** `curv_gaussian` (K = κ1·κ2 — both signs flip, the product
doesn't) and `path_curv_xy` (computed in-plane from the centerline, no surface
normal involved).

**Fixing it** needs a real orientation scheme, not a tweak: propagate a
consistent orientation across the face graph (breadth-first sign flooding), or
orient against the vessel centerline (the lumen axis is the natural "inward"
reference and is already computed in `pipeline/centerline.py`), or take the sign
from a signed distance field of the filled lumen. Deliberately deferred — it is
a design decision about what "outward" means on a branched tube.

## 2. Fibroblast / endothelial classifier is unvalidated

**Where:** `pipeline/gui/nuclei.py`, the "Count fibroblasts" step →
`n_fibroblasts`.

**What's wrong:** several compounding issues.

- It samples a **~2 µm shell around each nucleus** and asks whether mClover is
  above a threshold there. That perinuclear shell is precisely where *both*
  markers are dimmest — the sensor is out at the adhesions, VE-cadherin is out at
  the junctions — so the classifier reads the two cell types where they look
  most alike.
- **The result isn't stable.** The reported fibroblast fraction swings from
  **38% to 19%** across a single step of the threshold slider. A number that
  halves on one slider click is not measuring the biology.
- **The shell isn't masked against neighbouring nuclei.** In a dense field the
  shell around nucleus A overlaps nucleus B and B's cytoplasm, so a fibroblast
  next to an endothelial cell contaminates both readings.
- **Nucleus splitting is likely over-segmenting.** Splitting uses a ~3.9 µm
  minimum separation, but these nuclei are oblong and up to ~20 µm long, so one
  nucleus can be broken into several — inflating `n_nuclei_total` and the
  denominator of every fraction.

**What you see:** `n_fibroblasts` and `n_nuclei_total` in the CSV, currently not
trustworthy in either absolute value or ratio.

**Fixing it** starts with the diagnostic below (Deferred idea 2), not with more
threshold tuning — the current design has no validation signal at all, so
there's no way to tell a fix from a regression. Two cell types that are both
DAPI-positive and distinguished only by which marker they carry is the core
unsolved problem in this project.

## 3. `measure_core.py` duplicates the Measure worker and has diverged

**Where:** `pipeline/measure_core.py::compute_measurements` vs the `_run` worker
inside `pipeline/gui/measure.py`.

**What's wrong:** two near-identical implementations of the same measurement.
They have **already diverged**: `measure_core` has no X-Y ROI filter, so the
headless path measures every adhesion in the volume where the GUI measures only
those inside a drawn region. Any future fix has to be made twice — the
sum-then-divide ratio fix had to be applied in three places for exactly this
reason.

**What you see:** a headless CSV and a GUI CSV for the same image disagree
whenever an ROI is active. Both docstrings now say so, which is a patch on the
symptom.

**Fixing it** means extracting one implementation both call — a restructure that
should not share a commit with behaviour changes. Deliberately deferred.

## 4. Curvature is assigned with no distance limit

**Where:** `measure_core.py` / `gui/measure.py`, the curvature block.

**What's wrong:** every adhesion is assigned the curvature of its **nearest
surface vertex**, with no maximum distance. An adhesion 30 µm from the vessel
wall still gets a curvature value, reported with the same confidence as one
sitting on it. "Local curvature" stops being local.

**What you see:** `curv_*` columns are populated for every row, including rows
whose `dist_mesh_um` is large enough that the association is meaningless.

**Fixing it** is a one-line gate — NaN the curvature where
`dist_mesh_um > threshold` — but the threshold is a scientific choice about how
far an adhesion can sit from the wall and still be "on" it. Needs a number from
the user. **Interim workaround: filter on `dist_mesh_um` yourself when
analysing the CSV.**

## 5. Path curvature is undefined at branch nodes but still reported

**Where:** `pipeline/centerline.py` → `path_curv_xy`.

**What's wrong:** κ is computed by fitting along a 1D centerline. At a Y-junction
the centerline is not a 1D curve — three branches meet — so κ there is an
artefact of whichever branch the skeletonisation happened to trace through.

**What you see:** `path_curv_xy` values near junctions are arbitrary. The
`dist_branch_um` column exists and is the tool for this: **filter out small
`dist_branch_um` when using `path_curv_xy`.** Nothing does that automatically.

## 6. `cadherin_junction_proba` has no consumer

**Where:** `pipeline/gui/cadherin_rf.py`.

**What's wrong:** the VE-cad RF tab trains a classifier and produces a junction
probability map into `state["cadherin_junction_proba"]`, and then nothing reads
it. No measurement, no mesh step, no export consumes it. Its module docstring
also points at a headless `cadherin_rf_batch` script **that does not exist**.

**What you see:** a tab that appears to be part of the pipeline and is currently
a dead end. It is a reasonable path to endothelial-cell segmentation (a clean
junction network is far more segmentable than thresholded VE-cad), but that
consumer has not been written.

## 7. Corrections to the record

Statements in the deleted `CHANGES.md` that were wrong, recorded here so they
don't get rediscovered:

- **"The FRET constants are unset."** They were **set** in `params.json` — and
  not calibrated for these images. Set-but-uncalibrated is the dangerous state,
  because everything downstream runs and produces plausible force numbers. The
  block has now been removed from `params.json` and FRET is gated on its
  presence (`io.py::fret_configured`).
- **"The batch runner exports CSVs."** `pipeline.batch` does **not**. Its
  `_RUN_STEPS` list stops before Measure (which crashes headless), so it writes
  intermediates only. Use `pipeline.batch_measure` for headless CSVs.
- **Undocumented in PR #1:** the force-proxy tab was **removed** (superseded by
  the clover/ruby ratio; `pipeline/gui/force_proxy.py` has now been deleted), and
  the **VE-cad RF tab was added** (see problem 6).

## 8. Statistical caution for the force-vs-curvature goal

Not a code defect — a analysis trap this dataset is shaped to fall into.

**Adhesions are not independent observations.** There are many adhesions per
cell, several cells per image, and only a handful of images. Curvature also
varies slowly in space, so two adhesions near each other share both their cell
and their curvature value. Pooling every adhesion into one correlation treats
n = (number of adhesions) as the sample size when the effective n is closer to
(number of cells), or even (number of images) for anything curvature-related.

**Consequence:** a pooled Spearman p-value across all adhesions is driven almost
entirely by adhesion count. It will be tiny for effects that are not real. The
`Spearman(ratio, ...)` lines the Measure tab prints are a quick read for one
image, not evidence.

**Do instead:** compute the relationship **per image** first and look at whether
the direction is consistent across images. Then, if pooling, use a method that
respects the nesting (mixed-effects model with image and cell as random effects,
or per-image summary statistics combined across images). And check the depth
confound (`ratio_clover_ruby` vs `centroid_z_um`) before any of it — depth and
curvature covary on a tube, so optical attenuation masquerades as curvature.

---

# Deferred ideas

## 1. Progress / completion bars for long-running operations

**Want:** when a tab's compute button is clicked (e.g. Nuclei → "Compute DAPI
Otsu + components", Adhesion RF → "Run RF predict"), show how far along it is
and when it's done, instead of only terminal prints.

**Key nuance — not all steps can show a true percentage:**

- **Determinate (real %) is possible** for steps that loop over items:
  - `rf.predict_proba` already runs in N z-chunks and prints `chunk i/N` — a
    natural progress source.
  - Load-time resampling loops over channels (`io.to_isotropic`) and FRET
    outputs (`fret.resample_fret_outputs`).
  - FRET compute loops per-Z.

  These would need the blocking worker converted to a generator that `yield`s
  progress (and `rf.predict_proba` would need a progress callback threaded out
  of it).
- **Indeterminate only** for atomic single library calls, which have no
  sub-progress hook: `otsu_mask` (`threshold_otsu`), `label_and_sizes`
  (`ndi.label`), `niblack_mask`, watershed splitting. The honest option here is a
  "working…" bar that animates and clears on completion.

**Implementation options:**

1. **Inline `magicgui.ProgressBar` per tab** — bar appears in the tab, right
   where the button is. Matches the request most closely. More work: wire each
   worker's `yielded`/`started`/`finished` signals to a per-tab bar, and make the
   determinate workers yield progress.
2. **napari's built-in activity dock** (`napari.utils.progress`) — the hourglass
   panel, bottom-right. Much less code and idiomatic, but off to the side rather
   than in the tab.

**Open decisions:** inline vs activity dock; scope (all tabs vs just the heavy
steps vs a Nuclei-tab prototype first); whether atomic steps get a
busy/indeterminate bar or none at all.

**Note:** `pipeline/gui/logpanel.py` now implements `isatty`/`fileno`/`encoding`
on its stream, so a tqdm-style bar routed through stdout will at least behave
sanely in the Log dock rather than emitting one line per update.

## 2. Fibroblast / endothelial scatter diagnostic

**Want:** a validation signal for problem 2 above, before touching the
classifier itself.

**Idea:** one point per cell, plotted as **per-cell sensor signal (mClover)** vs
**per-cell-edge VE-cadherin signal** — measuring each marker where it actually
lives (the sensor through the cell body / adhesions, cadherin at the cell
boundary) rather than both in the same dim perinuclear shell. Score each cell
against the other cells **in the same image**, so illumination and expression
differences between images don't move the decision boundary.

**Why it comes first:** if the two populations separate into two clouds, the
classification problem is real and the current features are simply wrong — and
the scatter shows what the boundary should be. If they don't separate, no amount
of threshold tuning will work and the approach has to change (e.g. a
segmentation-based assignment rather than a per-nucleus shell). Either outcome is
more informative than the current situation, where there is nothing to check a
change against.

---

# Explicitly out of scope for the current commit

Recorded so the omissions are deliberate, not forgotten:

- **Curvature normal-sign fix** (problem 1) — needs a design decision.
- **`measure_core` deduplication** (problem 3) — a restructure that shouldn't
  share a commit with behaviour changes.
- **Fibroblast classifier** (problem 2) — needs the diagnostic first.
- **`dist_mesh_um` gate on curvature** (problem 4) — needs a threshold from the
  user.
