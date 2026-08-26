# Pipeline audit — 3d_napari_pipeline

> **This file is frozen.** It is a dated snapshot of what was measured on
> 2026-08-10 against commit `0faf2e5`, and it is **not** edited to reflect
> fixes — its value is that it records what was true before them. Only factual
> corrections to the measurements themselves belong here.
>
> **For current status, see [`ISSUES.md`](ISSUES.md)** — the living development
> log. Its entries cite the findings below by ID (`audit F1`, `audit F2`, …).
> A future audit becomes a new dated file alongside this one rather than an
> edit to it.
>
> **The commit hash changed on 2026-08-26, and nothing else did.** The branch
> `consolidation-fixes-and-docs` was rewritten to strip a co-author trailer from
> its commit messages, which changes every hash on it. The audited commit was
> `5077918` and is now `0faf2e5`; the tree, the author, the dates and the
> measurements below are byte-identical. This is the one kind of edit this file
> accepts — a pointer that would otherwise dangle — and it is recorded here
> rather than made silently.
>
> **Navigating the references below.** This audit was written against
> `TODO.md`, which was retired on 2026-08-10 and migrated into `ISSUES.md`.
> References to "`TODO.md` problem N" are historically accurate; the mapping is:
>
> | audit finding | old `TODO.md` | now |
> |---|---|---|
> | F1 | problem 3 | [P01](ISSUES.md#p01) |
> | F2 | problem 1 | [P02](ISSUES.md#p02) |
> | F3 | — | [P03](ISSUES.md#p03) |
> | F4 | — | [P04](ISSUES.md#p04) |
> | F5 | — | [P07](ISSUES.md#p07) |
> | F6 | problem 7 | [P08](ISSUES.md#p08) |
> | F7 | problem 6 | [P09](ISSUES.md#p09) |
> | F8 | — | [P10](ISSUES.md#p10) |
> | F9 | — | [P12](ISSUES.md#p12) |
> | F10 | — | [P13](ISSUES.md#p13) |
> | F11 | — | [P15](ISSUES.md#p15) |
> | F12 | — | [P16](ISSUES.md#p16), [P17](ISSUES.md#p17) |
> | — | problem 2 | [P06](ISSUES.md#p06) |
> | — | problem 4 | [P05](ISSUES.md#p05) |
> | — | problem 5 | [P11](ISSUES.md#p11) |
> | — | problem 8 | [P14](ISSUES.md#p14) |
> | — | problem 9 | [P18](ISSUES.md#p18) |
> | §4 | problem 10 | `ISSUES.md` § Analysis notes |

## 1. What was audited

| | |
|---|---|
| **Audited at commit** | `0faf2e5` "Record the ruby-contrast and depth-drift findings in TODO.md", branch `consolidation-fixes-and-docs` — **was `5077918` until 2026-08-26**, see the note below |
| **Date** | 2026-08-10 |
| **Scope** | Every module in `pipeline/` and `pipeline/gui/`, the three batch runners, `Scripts/depth_check.py`, and the exported CSV. ~8,600 lines |
| **Evidence image** | `V3_CorrectSettings_40xWater_04step_002` — the only image with a full set of saved intermediates and a shipped CSV (2,092 adhesions) |
| **Code changes made** | **None.** This is a read-only audit |

**Method.** Findings marked *measured* were reproduced against the shipped
artefacts in `V3_..._04step_002_saved/` — the mesh, the curvature arrays, the
adhesion instances, the DAPI labels, and the raw OME-TIFF. Where a defect was
suspected, the corrected version was computed and compared against the shipped
one. Findings marked *inspection* are read from the code without a measurement.

**Severity.** **S1** = a number reported as a result is wrong or meaningless ·
**S2** = a derived metric is wrong · **S3** = correctness/robustness defect with
no effect on the headline numbers · **S4** = documentation.

---

## 2. Headline

Three findings change what the pipeline can currently claim.

1. **All four surface-curvature columns are noise, not sign-flipped.**
   `alpha_shape_faces` sorts each triangle's vertex indices, which destroys the
   winding, so the "vertex normal" that defines the curvature fitting frame is a
   near-random direction. Correlation between the exported `curv_mean` and the
   same computation with winding preserved is **ρ = +0.048**. This supersedes
   `TODO.md` problem 3, and it means `curv_gaussian` is **not** protected as the
   README states.

2. **`ratio_clover_ruby` is driven by its denominator's noise floor.**
   ρ(ratio, ruby per voxel) = **−0.773**; ρ(ratio, clover sum) = **+0.068**.
   10.6% of ruby voxels inside adhesions are clipped to exactly zero by
   background subtraction. The observed p1–p99 spread is **4.8×**, against a
   ceiling of **1.77×** that the shipped LUT allows across 1–6 pN.

3. **Adhesion detection is depth-truncated.** Per unit of vessel surface, the
   far wall yields **37%** of the near wall's adhesion density. This is
   undocumented and it biases every depth-stratified comparison.

Taken together: the tension-vs-curvature question cannot be addressed with the
current columns — not because the effect is small, but because neither variable
is currently measuring what its name says.

Everything else found is smaller. The FRET gating, the ratio's direction and
sum-then-divide construction, the per-plane background subtraction, the
native-Z ratio grid, and the per-file voxel-size handling are all correct — see
§5.

---

## 3. Findings

### F1 — `alpha_shape_faces` destroys triangle winding; curvature is noise · S1 · `cell_mesh.py:156`

`alpha_shape_faces` has to sort each face's vertex indices to hash it and find
the faces that bound exactly one tetrahedron. It then **returns the sorted
faces**:

```python
faces = np.sort(faces, axis=1).astype(np.int64, copy=False)   # cell_mesh.py:156
...
return boundary.view(fc.dtype).reshape(-1, 3)
```

Sorting is correct as a *key*; returning it as the *geometry* is not. Triangle
orientation is carried entirely by vertex order, so after the sort every face
normal `cross(p1-p0, p2-p0)` points in a direction determined by index order
rather than by the surface. `vertex_outward_normals` then averages those
near-randomly-signed normals, and `mesh_curvature` builds its local tangent
frame from the result.

The away-from-centroid flip in `vertex_outward_normals` (which `TODO.md`
problem 3 identifies as the bug) is applied *on top of* an already-meaningless
vector. It is not the root cause.

**Measured** on `V3_..._04step_002` (8,966 verts, 15,786 faces, 7,282 surface
verts):

| | shipped | winding preserved |
|---|---|---|
| face-normal coherence at a vertex (median) | **0.209** | **0.945** |
| vertices with coherence < 0.5 | 93.3% | 19.7% |

Coherence is `‖Σ face normals‖ / (incident face count)`; 1.0 means every
incident face agrees. 0.209 is what randomly-signed unit vectors give.

Recomputing curvature both ways (radius 5 µm, smooth 10, identical code path
otherwise):

| column | shipped | winding preserved | ρ(shipped, correct) |
|---|---|---|---|
| `curv_mean` | median −0.0060, 43.6% positive | median −0.0383, 13.7% positive | **+0.048** |
| `curv_k1` | median +0.0170, 69.0% positive | median −0.0001, 49.9% positive | +0.090 |
| `curv_k2` | median −0.0317, 21.4% positive | median −0.0743, **3.1%** positive | +0.117 |
| `curv_gaussian` | median +0.0000, 52.4% positive | median +0.0001, 53.2% positive | **+0.187** |

`curv_mean` sign agreement between the two: **54.6%** (50% is a coin flip).

Two things confirm which one is right. The shipped values are a coin flip on
sign (43.6% / 52.4% positive) — no geometry does that. And with the winding
preserved the numbers become a textbook tube: **κ1 ≈ 0** (the axial direction),
**κ2 ≈ −0.074 1/µm and 96.9% consistently signed** (circumferential, radius
≈ 13.5 µm), **K ≈ 0**. That is exactly what a cylinder must produce.

**Consequences beyond the columns.** The Vessel-mesh tab's "implied caliber ~
1/median κ1" readout (`gui/vessel_mesh.py:368`) is computed from the broken κ1
and is meaningless; with the winding preserved the caliber comes from |κ2| and
is ~13.5 µm radius. The `Spearman(ratio, …)` lines printed by the Measure tab
for the four curvature columns are correlations against noise.

**Not affected:** `prune_by_normal` compares `|dot|` and so is sign-agnostic —
a face normal with an arbitrary sign is still a valid line direction there.
`mesh_from_volume` (marching cubes, used by the polar and alpha-2D methods)
returns consistently-wound faces, so **only the default method — `alpha 3D
(points)` — is affected.** That is the default in the GUI
(`METHOD_ALPHA3D`, `gui/vessel_mesh.py:96`) and the one that produced every
shipped mesh.

**The fix is small and local:** keep the unsorted face alongside the sorted key,
and orient each surviving boundary face away from the fourth vertex of its
owning tetrahedron. That is what the verification above did. It is a genuine
fix, not a re-derivation of "outward" — the orientation comes from the
tetrahedralization, so a Y-branch and the two-sided wall shell both come out
consistent. `TODO.md` problem 3 currently defers this as "a design decision
about what outward means"; with the winding intact, most of that decision is
already made by the alpha shape.

### F2 — `ratio_clover_ruby` tracks the ruby noise floor, not tension · S1 · data + `fret.py:20`

**Measured**, recomputing the per-adhesion channel sums exactly as the pipeline
does (native Z, per-plane background subtraction, sum-then-divide; reproduces
the shipped column at ρ = 1.0000, max relative difference 0.0000):

| | ρ with `ratio_clover_ruby` |
|---|---|
| clover sum | **+0.068** |
| clover per voxel | +0.337 |
| ruby sum | −0.329 |
| **ruby per voxel** | **−0.773** |
| adhesion voxel count | −0.036 |

A ratio is coupled to its denominator by construction, and `TODO.md` problem 1
correctly notes this — but then uses it to set the observation aside. The
*relative* magnitudes are the point: if both channels carried comparable real
signal, both would correlate substantially. **The numerator contributes
essentially nothing (+0.068).** The exported tension proxy is, to a good
approximation, an inverse map of the ruby channel.

**Why ruby carries no signal.** Measured inside adhesions, after per-plane
background subtraction:

| channel | signal/voxel | plane background | contrast | voxels clipped to exactly 0 |
|---|---|---|---|---|
| mClover | 716.9 | ~553 | 129.6% | 0.0% |
| mRuby | 34.1 | ~499 | 6.8% | **10.6%** |

`background_subtract` ends in `np.clip(img - mode_val, 0, None)`. On a channel
at 6.8% contrast this discards a tenth of the measurements outright and
rectifies the noise on the rest — an asymmetric, signal-dependent bias on
exactly the quantity in the denominator.

**The LUT bounds how much the ratio may vary.** The mRuby channel is direct
acceptor excitation (the file lists FRET, mClover and mRuby as three separate
cubes), so `I_D/I_A ∝ (1 − E)` and the direction is correct: higher ratio =
higher tension. But that also means the shipped calibration
(`Fret_LUT/Clover_mRuby2_GGSGGS7_force_new.txt`) sets a ceiling on the spread,
and the ceiling is scale-free — the gain constants cancel in a ratio of ratios:

| force range | E spans | maximum `clover/ruby` spread |
|---|---|---|
| 2–5 pN | 0.400 → 0.218 | 1.30× |
| 1–6 pN | 0.538 → 0.181 | **1.77×** |
| 0.5–8 pN | 0.677 → 0.124 | 2.71× |

**Observed: p1 = 10.15, p99 = 49.05 → 4.8×** (full range 6.7–151.5 → 22.6×).
Even on the robust p1–p99 figure, and even against the most generous 0.5–8 pN
window, the spread exceeds what the sensor can physically produce. Most of the
variance in the column is not tension.

This is largely an acquisition problem, as `TODO.md` says. The audit's addition
is that it is now *quantified against the sensor's own calibration*, and that
the numerator is contributing almost nothing — which is a stronger statement
than "the ratio is noisy."

### F3 — Adhesion detection is depth-truncated · S1 · `gui/adhesion_rf.py:138,163`

The adhesion RF runs on the clover channel (`vol_f32 = clover_raw.astype(...)`)
and prediction is gated by default on a **global Otsu mask of that same
channel** (`gate_by_clover` defaults `True`, 1 µm dilation). Clover attenuates
with depth. A global threshold on an attenuating channel therefore accepts less
of the volume the deeper you go — a selection effect upstream of every column.

**Measured** — adhesion count normalised by the vessel surface area in each
depth band (triangle areas from the shipped mesh):

| depth band | vessel surface (µm²) | adhesions | per 1000 µm² |
|---|---|---|---|
| 0–20 µm | 134,908 | 1402 | **10.39** |
| 20–40 µm | 46,667 | 98 | 2.10 |
| 40–60 µm | 28,928 | 21 | 0.73 |
| 60–80 µm | 36,480 | 106 | 2.91 |
| 80–105 µm | 128,890 | 465 | 3.61 |

Aggregated: near wall 8.26 per 1000 µm², far wall 3.05 — the **far wall yields
37% of the near wall's detection density across comparable surface area**
(181,575 vs 194,298 µm²).

Median adhesion volume falls only 11% across the same span (14.2 → 12.6 µm³),
so this is not adhesions shrinking with depth. It is adhesions going undetected.

**Consequence.** Any depth-stratified comparison contrasts a near-complete
shallow sample against roughly the brightest third of the deep one. That
truncation is itself a plausible generator of a "depth effect" in the tension
proxy — the surviving deep adhesions are a biased draw. `TODO.md` problem 2
notes that adhesion *volume* does not vary with depth (ρ −0.007) and reads that
as evidence against a segmentation artefact; volume is the wrong statistic,
because truncation removes objects rather than shrinking them.

### F4 — Two different nucleus label maps in one CSV · S1 · `gui/nuclei.py:211` vs `gui/measure.py:57`

The Nuclei tab computes two separate segmentations and only one of them reaches
each consumer:

- **Otsu + connected components**, stored in `state["dapi_state"]["labels"]`.
  Not split, not size-filtered beyond ≥10 voxels.
- **Watershed-split nuclei** (`nuc`, `gui/nuclei.py:211`), built from the
  size-filtered binary. **Never written back to `dapi_state`.**

Measure reads `state["dapi_state"]["labels"]` (`gui/measure.py:57`), so
`dist_nuc_um` and `nuc_id` come from the unsplit map. `n_fibroblasts` /
`n_nuclei_total` come from the split map.

**Measured:** the unsplit map has **362** components; `n_nuclei_total` is
**273**. `nuc_id` values in the CSV run to 355 and index the 362-label map.
**The two cannot be joined** — there is no way to ask whether a given adhesion's
nucleus was classified as a fibroblast, which is precisely the join the project
needs.

Two further consequences:

- **Measure ignores the `min_size` slider.** The user tunes it in the Nuclei
  tab; it drives the display, the fibroblast count, and nothing else. The
  unsplit map's median component is **81.7 µm³** — far below a real nucleus
  (an oblong 20 × 5 µm nucleus is ~260 µm³), so more than half of what
  `dist_nuc_um` measures distance *to* is sub-nuclear debris.
- `TODO.md` problem 4 states that over-splitting "inflates `n_nuclei_total` and
  the denominator of every fraction." Measured, splitting **reduces** the count
  (362 → 273), because the size filter and the border loss in F5 remove more
  than the splitting adds. The over-splitting concern may still be real
  per-nucleus; the stated *direction of the effect on the count* is not.

### F5 — `peak_local_max` drops border seeds, in two places · S2 · `gui/nuclei.py:207`, `rf.py:612`

`peak_local_max` defaults to `exclude_border=True`, which excludes peaks within
`min_distance` of the array edge (verified on the installed scikit-image 0.26.0;
with `min_distance=5`, a peak at index 3 is dropped and only reappears with
`exclude_border=False`).

**(a) Nucleus splitting** — `gui/nuclei.py:207` calls
`peak_local_max(dist, min_distance=md, labels=binary)` with `md ≈ 12` voxels
(3.9 µm at 0.325 µm/voxel). Nuclei whose distance peak falls within 12 voxels of
any face of the volume get **no watershed seed**. A connected mask component
with no marker is left at label 0 by `watershed(..., mask=binary)` — so those
nuclei are dropped entirely, or absorbed into a neighbour that does have a seed.
On a 324-plane stack the excluded z shell is 7.4% of the depth, and it sits
exactly where the near-wall and far-wall adhesion populations are.

**(b) Adhesion splitting** — `rf.py:612` crops each candidate component to its
bounding box with `pad = 2`, then calls
`peak_local_max(dist, min_distance=5, labels=crop_mask)`. **The border exclusion
zone (5 voxels) is larger than the padding (2 voxels)**, so seeds within 3
voxels of the component's own bounding box are suppressed. For elongated
adhesions the end seeds are the ones lost, so merged elongated objects
systematically fail to split. `min_distance=5` is also hard-coded (1.6 µm) and
not exposed.

### F6 — `path_curv_xy` is contaminated over `smooth_um`, and is unsigned · S2 · `centerline.py:128`

Two separate problems in the one column that currently correlates with the
tension proxy.

**(a) Contamination radius.** κ is fitted over
`tree.query_ball_point(p, smooth_um)` — a **Euclidean** ball of 20 µm (default)
over the centerline points. Near a junction that ball contains points from two
or three different branches, so the parabola is fitted to a Y-shaped point set.
The same happens wherever the vessel doubles back within 20 µm of itself.

`TODO.md` problem 7 correctly says κ is undefined *at* branch nodes and points
at `dist_branch_um` as the filter. What it does not say is that **the
contamination radius is `smooth_um`, so the filter threshold must be at least
20 µm** (and should track the slider). Measured on this image, filtering
`dist_branch_um > 20 µm` costs 2.6% of rows — cheap, and currently nobody does
it.

**(b) κ is unsigned.** `centerline.py:139` returns `abs(2·c₂)/(1+c₁²)^1.5`.
Combined with the fact that κ is looked up from the **centerline**, every
adhesion at a given axial position receives the same κ regardless of which side
of the vessel it sits on. So `path_curv_xy` cannot distinguish the inside of a
bend from the outside — arguably the most mechanically interesting contrast
available on a curved tube. This is a design limitation, not a bug, but it
should be stated before the column is used as *the* curvature variable.

Also worth noting: the shipped centerline has median κ = 0.0040 1/µm (radius
250 µm) on a 332 µm field, i.e. **most of the vessel is nearly straight**, with
the signal concentrated in a few high-κ points near the junction — the same
points contaminated by (a).

### F7 — `dist_mesh_um` is measured to vertices the mesh rejected · S2 · `measure_core.py:103`, `gui/measure.py:187`

`vessel_tree = KDTree(vessel_verts)` uses **every** vertex; the curvature lookup
correctly restricts to surface vertices (`gui/measure.py:238`). For the default
alpha-3D method, `vessel_verts` is the whole subsampled cadherin point cloud,
so the two queries run against different sets. Worse, `filter_faces_by_extent`
(`min_component_um = 50` by default) deliberately drops the faces of small
debris islands — and those islands' **vertices stay in the `dist_mesh_um`
tree**. Distance is measured to geometry the mesh step explicitly rejected.

**Measured:** 1,684 of 8,966 vertices (**18.8%**) belong to no face; 15.1% of
adhesions have one of them as their nearest vertex. The numerical effect is
small — median understatement 0.00 µm, p90 0.37 µm — so this is a semantics
defect rather than a large error. The related point that *does* matter: surface
vertex spacing is **2.25 µm** (median nearest-neighbour), so `dist_mesh_um` is
a distance-to-a-point-cloud with ~2 µm granularity, not a distance to a
surface. Any threshold chosen for the `TODO.md` problem 6 gate has to be read
against that number, and it means different things for the alpha-3D method
(point cloud) than for the marching-cubes methods (vertex spacing ≈ 1 voxel).

**On problem 6 itself:** measured, 98.6% of adhesions are within 10 µm of the
mesh and only 2 rows exceed 20 µm (max 54.2 µm). The defect is real but affects
a handful of rows, not the bulk of the table — worth recording so the gate isn't
prioritised above F1–F4.

### F8 — `force_total_pN` is not a force · S2 · `measure_core.py:181`, `gui/measure.py:280`

`force_total_pN = force_mean_pN × props["area"]`, and `props["area"]` is a
volume in µm³ (`regionprops_table` was given `spacing=scale`). The column's
units are pN·µm³. Currently blank on this dataset, so it misleads no one yet —
but it will the moment a `fret` block is configured.

### F9 — Headless measure fails silently where the GUI warns · S3 · `batch_measure.py:143`

Two asymmetries between the paths:

- `gui/vessel_mesh.py:501` prints an explicit warning when a loaded curvature
  file's length does not match the mesh, because Measure will silently NaN every
  curvature column. `batch_measure.py` loads the same two files with no such
  check, so a stale pairing produces a CSV with four blank columns and no
  message.
- `_standard_instances` **overwrites `adhesion_instances.npz`** as a side effect
  of measuring (`batch_measure.py:123`). It is documented in the docstring and
  the reasoning is sound, but a command named "measure" silently rewriting a
  saved segmentation is a surprising default; it deserves a flag.

### F10 — Tunables that are not sliders · S3 · `gui/nuclei.py:241`, `rf.py:612`

The project's stated working agreement is that tunable parameters become
magicgui sliders. Two thresholds that materially change results are hard-coded:

- `is_fib &= cad_means < cbg * 1.3` — the VE-cadherin exclusion multiplier, on a
  channel that is dim in the perinuclear shell for both cell types.
- `min_distance=5` in the adhesion watershed split.

### F11 — Smaller defects · S3 · inspection

| | finding |
|---|---|
| `centerline.py:41` | `_prune_spurs` mutates `skel` while iterating the endpoint list, but `S`/`nbrs` were built before the mutation. Paths can be traced through already-deleted pixels within one pass. Converges, but the pruned result depends on iteration order |
| `measure_core.py:80` | Adhesions thinner than one native z-step can fall entirely between sampled planes, yielding `r_sum = 0` → NaN ratio. Correctly NaN'd rather than wrong, but it is a size-dependent (and therefore depth-dependent) dropout |
| `vessel_mesh.py:172` | `alpha_fill_slice` rasterises kept triangles in a Python loop over `skimage.draw.polygon`. Fine for the default method, slow for alpha-2D on a large mask |
| `vessel_mesh.py:233` | `close_volume` uses a cubic (box) structuring element, so the closing reaches √3× further along diagonals than the requested `close_um` |
| `denoise.py:107` | `np.convolve(..., mode="same")` on an even-length kernel introduces a half-bin phase shift in the polar radius smoothing |

### F12 — Documentation corrections · S4

| where | current text | measured |
|---|---|---|
| `README.md` "What not to trust yet" | "`curv_gaussian` and `path_curv_xy` are unaffected (the first is a product of two curvatures, so the signs cancel)" | ρ(shipped, correct) = **+0.187** for `curv_gaussian`. The reasoning holds for a global sign flip; the actual defect is a corrupted *frame*, which the product does not cancel. Separately, K ≈ 0 on a tube by construction, so even when correct the column has almost no dynamic range — median +0.0000 |
| `TODO.md` #3 | "vertex normals are oriented by pointing them away from the mesh centroid… for a branch on the far side of the junction the normal flips" | The centroid heuristic is real but secondary. Root cause is F1; the columns are noise, not sign-flipped, and 49.9% of vertices are flipped — not a branch-local minority |
| `TODO.md` #4 | "Nucleus splitting is likely over-segmenting… inflating `n_nuclei_total`" | Splitting *reduces* the count on this image (362 unsplit components → 273). Direction of the stated effect is wrong; see F4/F5 for the actual mechanism |
| `TODO.md` #6 | "An adhesion 30 µm from the wall still gets a curvature value" | True but rare: 98.6% within 10 µm, 2 rows beyond 20 µm |
| `TODO.md` #1 | "the ratio is inherently coupled to its denominator, so a strong ratio–ruby correlation is expected; it is the contrast figures that make it a problem" | Correct as far as it goes, but it sets aside the decisive comparison. ρ(ratio, clover sum) = +0.068 says the numerator contributes essentially nothing — see F2 |
| `TODO.md` #2 | "shallow 25% 25.55 → deep 25% 19.71, −23%" | Reproduced, but the quartile series is **non-monotonic**: 25.55, 18.01, 23.45, 19.71. Overall ρ(ratio, z) = −0.187. The effect is real and modest, not a clean gradient — and the quartiles straddle a bimodal depth distribution, so they mix the two walls |
| `gui/cadherin_rf.py` docstring | points at a headless `cadherin_rf_batch` script | Confirmed: does not exist (`TODO.md` #8 already records this) |

---

## 4. Statistics — the design, and a correction to an earlier suggestion

### The confound is worse than "near wall vs far wall"

The working assumption has been that depth on a tube is essentially which wall
you are on, so stratifying by wall would give two internal replicates spanning
the same range of axial curvature. **Measured, that does not hold on this
image.**

The two depth clusters occupy largely *different* regions of the field:

| | n | median y (µm) | median x (µm) |
|---|---|---|---|
| near (z < 40 µm) | 1500 | 164.7 | 172.2 |
| far (z ≥ 40 µm) | 592 | 144.4 | 117.6 |

Gridding the field at ~14 µm, the two populations share only **26%** of the
cells they occupy (Jaccard 0.26). Depth also correlates with in-plane position
directly: ρ(z, y) = −0.293, ρ(z, x) = +0.202, and ρ(z, `dist_branch_um`) =
**−0.461**. The vessel dives through the stack, so **depth is substantially
"where along the vessel you are"**, not just "which wall."

That has a hard consequence: depth cannot be separated from axial position on
this image, and axial position is what determines curvature. A wall-stratified
design has little support here. It should still be checked per image — the
geometry may differ across the six — but it is not the escape route.

### The correlation structure is one spatial gradient

Everything that is a spatial coordinate correlates with the tension proxy, and
nothing else does:

| variable | ρ with `ratio_clover_ruby` | pooled p |
|---|---|---|
| `dist_branch_um` | +0.202 | 9e-21 |
| `path_curv_xy` | +0.191 | 2e-18 |
| `centroid_z_um` | −0.187 | 6e-18 |
| `dist_mesh_um` | +0.043 | 0.05 |
| `curv_k2` | +0.039 | 0.07 |
| `curv_gaussian` | −0.022 | 0.31 |
| `curv_k1` | −0.022 | 0.30 |
| `volume_um3` | −0.036 | 0.10 |
| `curv_mean` | +0.014 | 0.52 |
| `dist_nuc_um` | −0.020 | 0.37 |

The three that "work" are the three that encode position. The four curvature
columns sit at |ρ| ≤ 0.04 — exactly what F1 predicts. Splitting by wall does not
resolve it; it amplifies it (far wall: ρ(ratio, `dist_branch_um`) = +0.432,
ρ(ratio, z) = −0.330, ρ(ratio, `path_curv_xy`) = +0.292 — everything correlates
with everything, the signature of one underlying gradient).

### What follows for the analysis

1. **The unit of analysis is not the adhesion.** Curvature is a property of the
   vessel; two adhesions 5 µm apart share it. The independent-observation count
   for a curvature question is roughly (vessel length ÷ curvature correlation
   length) — tens per image, not thousands. A pooled Spearman over 2,092
   adhesions reports adhesion count.

2. **The cell level is currently unavailable.** `ratio ~ curvature + (1|image/cell)`
   needs cell identity, and F4 shows `nuc_id` indexes a different segmentation
   from the one that was classified — and an unfiltered one at that. Until F4 is
   fixed you have `(1|image)` with n = 6, or a spatial clustering of adhesions
   as a rough cell proxy.

3. **Use a spatial null, not a permutation null.** Shuffling tension across
   adhesions destroys spatial autocorrelation and produces a null far too narrow
   — it will manufacture significance. Circularly shift the tension field
   against the curvature field along centerline arclength within an image and
   rebuild the statistic over many shifts. That preserves both fields'
   autocorrelation and destroys only their alignment. Given the correlation
   structure above, this is not optional.

4. **Run the falsification battery identically**, decided in advance: clover
   alone vs curvature; ruby alone vs curvature; `volume_um3` vs curvature; and
   adhesions from image A against the curvature field of image B (must be null).

5. **Pre-specify the exclusions** before looking at tension: the `dist_mesh_um`
   gate, `dist_branch_um > smooth_um` (F6), the depth band, and whether `09step`
   is in. Each is otherwise a fork.

6. **State the ceiling.** With 6 images, if all agree in direction a sign test
   floors at p = 0.031. "In all six images, tension rank increases with axial
   curvature" is defensible. A slope in pN per unit curvature is not, and will
   not be until FRET is calibrated.

---

## 5. What was checked and is correct

Recorded so it isn't re-audited.

- **`fret.py` math.** `fret_correct` implements the standard three-cube
  correction, and the efficiency and LUT interpolation are consistent with the
  frozen original. Not modified.
- **FRET gating.** `io.fret_configured` keys on the *presence* of the `fret`
  block rather than a boolean flag, so there is no default to get wrong.
  Verified: `force_mean_pN`, `force_std_pN`, `force_total_pN` are all-NaN in the
  shipped CSV and export blank. This is the right design and the right default.
- **Ratio direction and construction.** `c_sum / r_sum` in all three
  implementations (`measure_core.py:91`, `gui/measure.py:173`,
  `gui/adhesion_rf.py:634`) — no divergence. Direction is physically correct
  (donor over direct-acceptor ∝ 1 − E, so higher = more tension). Sum-then-divide
  is right, and so is doing it on the native-Z grid.
- **Per-plane background subtraction.** `_bs_per_z` is correct in principle —
  the modal background does shift with depth, so a per-volume estimate would be
  wrong. Its failure on ruby (F2) is a signal-level problem, not a design error.
- **The quadric curvature math.** Given a correct normal, `mesh_curvature`'s
  Monge-patch second fundamental form (`E = 1 + d²`, `L = 2a/√(d²+e²+1)`, shape
  operator `I⁻¹II`) is correct, as is restricting the fit to face-referenced
  vertices and NaN-ing rank-deficient fits. The defect is entirely in the frame
  it is handed.
- **Voxel-size handling.** Read per-file from OME metadata, no `params.json`
  key to disagree with the image, each file resampled onto its own isotropic
  grid. Correct, and the right call.
- **`gpu.py` OOM discrimination.** Falls back to CPU only on
  `cudaErrorMemoryAllocation` / `cudaErrorAlreadyMapped` and re-raises every
  other CUDA failure. This is the correct and unusually careful choice.
- **`Scripts/depth_check.py`.** Well constructed — it already computes the
  ratio-vs-clover / ratio-vs-ruby comparison that produced F2, sizes its binned
  markers by n, and warns about the bimodal depth distribution. The finding was
  available from the existing tool; it was the interpretation in `TODO.md` that
  set it aside.
- **`TODO.md` #5 (duplicate measurement path) and #10 (pseudoreplication)** are
  accurate as written and are confirmed by this audit.

---

## 6. Suggested order

Dependency-ordered, not severity-ordered.

1. **F1 — preserve the winding in `alpha_shape_faces`.** Smallest fix with the
   largest effect; it is the precondition for any curvature result, and it makes
   `TODO.md` problem 3 mostly moot. Acceptance test: κ1 ≈ 0, κ2 consistently
   signed at ≈ −1/R, K ≈ 0, and ρ(`curv_mean`, `centroid_z_um`) → 0.
2. **F4 — one nucleus map.** Write the split labels back to `dapi_state`, or
   have Measure use the split map. Unlocks the cell level for the statistics and
   makes `n_fibroblasts` joinable to `nuc_id`.
3. **F5 — `exclude_border=False`** in both `peak_local_max` calls, and raise the
   adhesion crop padding above `min_distance`.
4. **F3 — quantify the detection truncation** on the other five images before
   any depth-corrected analysis. If it replicates, the deep sample is a biased
   draw and depth correction alone will not fix it.
5. **F2 — decide what to do about ruby.** Mostly acquisition. On the code side,
   `background_subtract`'s `bin_size=10` and its clip-at-zero are both wrong for
   a 6.8%-contrast channel, but the function is frozen and shared with the FRET
   path, so this needs an explicit decision rather than a patch.
6. **F6/F7/F8** — column semantics; cheap, and they make the CSV honest.
7. Then, and only then, the tension-vs-curvature analysis under §4.
