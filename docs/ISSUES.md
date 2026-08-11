# Known problems — the development log

**This is the development record, not the README.** It is the single place where
defects in this pipeline are tracked from the day they are found to the day they
are fixed. If a number the pipeline exports cannot be trusted, the reason is
here.

Read this before trusting or extending any output.

- **`README.md`** is the public description of what the project is and how to
  run it. It carries a short limitations table that links back here.
- **`docs/audit.md`** is a frozen, dated snapshot of a full pipeline audit
  (2026-08-10, commit `5077918`). Entries below cite it as `audit F1` and so on.
  It is evidence, not state: it is never edited to reflect a fix.
- **This file** is state. It changes every time a problem is opened, verified,
  or closed.

---

## How to use this file

**IDs are permanent.** `P01` always means the same problem. A retired ID is
never reused. Entries stay in the file after they are closed — a claim that
turned out to be wrong is worth as much as one that was right, because it stops
the same wrong claim being rediscovered.

**Three dates, because they answer different questions.**

| field | question it answers |
|---|---|
| **Opened** | when the problem was first recorded |
| **Last verified** | whether the numbers quoted in the entry are still current |
| **Resolved** | when it was fixed, and in which commit |

Evidence goes stale. A measurement quoted under **Last verified: 2026-08-10** was
true of the code and data on that date and may not be true now.

**Statuses**

| status | meaning |
|---|---|
| `open` | confirmed, not fixed |
| `investigating` | being worked on now |
| `fixed` | resolved; carries a date and commit |
| `mitigated` | worked around, root cause still present |
| `wontfix` | deliberate; carries a reason |
| `invalid` | investigated and found not to be a problem |
| `superseded by Pnn` | folded into another entry |

**Severity**

| | |
|---|---|
| **S1** | a number reported as a result is wrong or meaningless |
| **S2** | a derived metric is wrong |
| **S3** | correctness or robustness defect, no effect on the headline numbers |
| **S4** | documentation |

---

## Status at a glance

| ID | title | sev | status |
|---|---|---|---|
| [P01](#p01) | Curvature columns are noise — `alpha_shape_faces` destroys winding | S1 | **fixed** |
| [P02](#p02) | `ratio_clover_ruby` tracks the ruby noise floor, not tension | S1 | open |
| [P03](#p03) | Adhesion detection is depth-truncated | S1 | open |
| [P04](#p04) | Two nucleus label maps; `nuc_id` is not joinable to cell type | S1 | open |
| [P05](#p05) | Fibroblast / endothelial classifier is unvalidated | S1 | open |
| [P06](#p06) | Tension readout drifts with imaging depth | S2 | open |
| [P07](#p07) | `peak_local_max` drops border seeds, in two places | S2 | **fixed** |
| [P08](#p08) | `path_curv_xy` — 20 µm contamination radius, and κ is unsigned | S2 | open |
| [P09](#p09) | `dist_mesh_um` measures to rejected geometry; no gate on curvature | S2 | open |
| [P10](#p10) | `force_total_pN` is not a force | S2 | open |
| [P11](#p11) | `measure_core` duplicates the Measure worker and has diverged | S3 | open |
| [P12](#p12) | Headless measure fails silently and overwrites saved instances | S3 | open |
| [P13](#p13) | Tunable parameters that are not sliders | S3 | open |
| [P14](#p14) | `cadherin_junction_proba` has no consumer | S3 | open |
| [P15](#p15) | Minor defects (five items) | S3 | open |
| [P16](#p16) | Nucleus splitting inflates `n_nuclei_total` | — | **invalid** |
| [P17](#p17) | `curv_gaussian` is unaffected by the curvature bug | — | **invalid** |
| [P18](#p18) | FRET constants set but not calibrated | — | **fixed** |
| [P19](#p19) | `auto_min_size_otsu` is defeated by a hard-coded floor | S3 | open |
| [P20](#p20) | `subsample_um` silently sets the vessel mesh's resolution | S2 | open |
| [P21](#p21) | One global alpha cannot be correct on a tapering vessel | S2 | open |
| [P22](#p22) | Centerline overlay is drawn at a single wrong z | S3 | **fixed** |
| [P23](#p23) | Centerline preprocessing amputates tapering vessels | S2 | open |
| [I01](#i01) | Progress / completion bars for long-running operations | idea | open |
| [I02](#i02) | Fibroblast / endothelial scatter diagnostic | idea | open |
| [I03](#i03) | Objective-scored auto-calibration of the vessel mesh | idea | open |
| [H01](#h01) | Repository has no LICENSE | — | open |

**Suggested order of work** (dependency-ordered, not severity-ordered). The
original ordering was **P01 → P04 → P07 → P03 → P02 → P08/P09/P10**, then the
analysis in [Analysis notes](#analysis-notes) — P01 first because it is the
smallest fix with the largest effect and every curvature result depends on it.

**As of 2026-08-10 P01, P07 and P22 are fixed in the code**, so the queue is now
**P04 → P03 → P02 → P08/P09/P10**. One thing carries over: no CSV has been
regenerated yet, so every curvature column on disk is still the pre-P01 noise and
the [What not to trust yet](../README.md#what-not-to-trust-yet) rows stay as they
are until a new export exists. Fixing the code and fixing the data are separate
steps.

The 2026-08-10 vessel-mesh review (P19–P23, [I03](#i03)) does not displace that
ordering. The border-flare half of [P08](#p08) is independent of everything above
and can land at any time. [P21](#p21) and [I03](#i03) are one piece of work —
reconstruction method plus its calibration — and are now unblocked: corrected
curvature can be part of the calibration objective rather than something checked
afterwards.

---

# Open defects

<a id="p01"></a>
## P01 — Curvature columns are noise: `alpha_shape_faces` destroys triangle winding

| | |
|---|---|
| **Status** | fixed |
| **Severity** | S1 |
| **Opened** | 2026-08-10 (audit F1) |
| **Last verified** | 2026-08-10 |
| **Resolved** | 2026-08-10 (commit `5d0a9c6`) |
| **Where** | `pipeline/cell_mesh.py:137` (`alpha_shape_faces`), `pipeline/cell_mesh.py:246` (`vertex_outward_normals`) |
| **Affects** | `curv_mean`, `curv_k1`, `curv_k2`, `curv_gaussian`; the Vessel-mesh tab's "implied caliber" readout; the Measure tab's `Spearman(ratio, curv_*)` lines |
| **Evidence** | `docs/audit.md` §3 F1 |

`alpha_shape_faces` sorts each face's vertex indices in order to hash it and find
the faces bounding exactly one tetrahedron — then **returns the sorted faces**.
Sorting is correct as a *key* and wrong as *geometry*: triangle orientation is
carried entirely by vertex order, so afterwards every face normal points in a
direction set by index order rather than by the surface.

`vertex_outward_normals` then averages those near-randomly-signed normals, and
`mesh_curvature` builds its local tangent frame from the result. The
away-from-centroid sign flip in `vertex_outward_normals` is applied *on top of* a
vector that is already meaningless — it is not the root cause.

> **Correction, 2026-08-10 (left in place above so the original reading is on the
> record).** The last sentence is wrong, and it mattered. The away-from-centroid
> flip is a *second, independent* root cause, not a harmless passenger. Once the
> winding is preserved it is applied on top of a vector that is now meaningful,
> and it **overrides the true winding on 44.4% of the 7,282 surface vertices** —
> because on a vessel that branches and dives through the stack "away from the
> global centroid" is simply the wrong reference direction. Fixing only the
> winding leaves the curvature broken and fails every acceptance criterion below
> (see the 2026-08-10 History entry). **Both causes had to be fixed together.**

Measured on `V3_..._04step_002` (8,966 verts, 7,282 surface verts):

| | shipped | winding preserved |
|---|---|---|
| face-normal coherence at a vertex (median) | **0.209** | **0.945** |
| ρ(shipped `curv_mean`, correct `curv_mean`) | **+0.048** | — |
| `curv_mean` sign agreement | **54.6%** (50% = coin flip) | — |

With the winding preserved the numbers become a textbook tube: κ1 ≈ 0 (axial),
κ2 ≈ −0.074 1/µm and **96.9% consistently signed** (circumferential, radius
≈ 13.5 µm), K ≈ 0. That is what a cylinder must produce, and it is the strongest
evidence that the corrected version is the right one.

**Only the default mesh method is affected.** `mesh_from_volume` (marching cubes,
used by the polar and alpha-2D methods) returns consistently-wound faces. The
GUI default is `alpha 3D (points)` (`pipeline/gui/vessel_mesh.py:96`), which is
what produced every mesh saved so far.

**Not affected:** `prune_by_normal` compares `|dot|` and is sign-agnostic.

**The fix, part A — preserve the winding (`alpha_shape_faces`).** Keep the
unsorted face alongside the sorted key, and orient each surviving boundary face
away from the fourth vertex of its owning tetrahedron. The orientation then comes
from the tetrahedralization, so a Y-branch and the two-sided wall shell both come
out consistent — which removes most of the "what does *outward* mean on a
branched tube" design question that [P17](#p17)'s predecessor deferred.

**The fix, part B — orient from the mesh, not the centroid
(`vertex_outward_normals`).** Preserving the winding fixes the orientation only
up to one global sign, and the away-from-centroid heuristic that used to supply
that sign gets it wrong on 44.4% of vertices (see the Correction above). Replace
it with a single global decision read off the mesh: test closure (the fraction of
undirected edges shared by exactly two faces — **98.45%** on the shipped mesh, so
gate at 0.90), and on a closed surface flip the whole normal field by the sign of
the enclosed volume `Σ p0·(p1×p2)/6`. One decision for the whole mesh cannot
disagree with itself the way a per-vertex heuristic does. Open surfaces
(marching-cubes patches; a cell mesh after `prune_by_normal` deletes faces)
enclose no volume, so they keep the old centroid rule as a fallback — the new
test cannot apply there and behaviour is unchanged.

**Acceptance test:** ρ(`curv_mean`, `centroid_z_um`) → 0; κ2 consistently signed
at ≈ −1/R; K ≈ 0; face-normal coherence > 0.9. **Passed** — see below.

**History**
- 2026-08-06 — recorded in `TODO.md` as "curvature sign inverts on branched
  vessels", attributing the defect to the away-from-centroid heuristic.
- 2026-08-10 — full audit found the true root cause upstream. The columns are
  **noise, not sign-flipped**, and 49.9% of vertices are flipped, not a
  branch-local minority. Old entry superseded; [P17](#p17) closed as invalid as a
  consequence.
- 2026-08-10 — **fixed, in two parts.** The audit's root cause was real but
  incomplete: the winding fix *alone* moves almost nothing, and the entry's own
  acceptance test fails on it. Measured in memory on the shipped
  `V3_..._04step_002` mesh (8,966 verts / 7,282 surface verts, `alpha3d_um=15.0`,
  `subsample_um=5.0`, `min_component_um=253.03`, curvature `radius_um=5.0`,
  `smooth_iter=10`):

  | | shipped | winding only (part A) | **both (A+B)** | target |
  |---|---|---|---|---|
  | κ1 median (axial, want ≈ 0) | +0.0170 | +0.0183 | **−0.0027** | ≈ 0 |
  | κ2 median (1/µm) | −0.0317 | −0.0388 | **−0.0703** | ≈ −1/R |
  | implied R | 31.6 µm | 25.8 µm | **14.2 µm** | ≈ 13.5 µm |
  | κ2 sign consistency | 78.6% | 77.1% | **98.3%** | ≥ 96.9% |
  | K median (1/µm²) | +0.00002 | +0.00015 | **+0.00015** | ≈ 0 |
  | ρ(`curv_mean`, `centroid_z_um`) | −0.059 | −0.011 | **+0.042** | → 0 |

  Winding checks on the rebuilt mesh: directed edges with a reversed twin
  **53.7% → 100.0%**; face-normal coherence at a vertex (median) **0.209 →
  0.972**; closed-edge fraction unchanged at 98.45%.

  **Regression guard passed:** the face *set* is byte-identical before and after
  part A — 15,790 raw, 15,786 after `filter_faces_by_extent`, and the same set as
  the shipped `vessel_mesh.npz`. Only the winding of each row moved, so
  connectivity and every downstream row-subset are untouched.

  Per-vertex ρ(shipped `curv_mean`, corrected `curv_mean`) measures **+0.069**
  with **54.8%** sign agreement, reproducing `audit.md` F1's +0.048 / 54.6% and
  confirming the shipped column really was noise. Note that the audit's ρ is
  **per-vertex**; the acceptance test's ρ is **per-adhesion** (assigned exactly as
  `measure_core` does — nearest surface vertex per voxel, averaged over the
  adhesion) and the two are not the same number.

  Independent check of part B on a synthetic watertight mesh: a marching-cubes
  sphere (r = 20) is closed-fraction 1.0 and comes out **inward**-wound
  (signed volume −33,360 against a true 33,510), and the volume test flips it so
  100% of its normals point outward. The vessel alpha shape's own winding is
  already outward (6V = +2.82e6), so for that mesh part B amounts to *not*
  applying the centroid flip — which is the whole 44.4%.

  Blast radius traced and accepted: `gui/cell_mesh.py:162,:441` (the "Nuclei
  normals" arrows may reverse; `offset_centroids` only moves if
  `inward_offset_um > 0`, GUI default `0.0`) and `gui/clover.py:52-55`, whose
  one-sided `sd > 0` signed-distance gate can change membership. The other four
  signed-distance consumers are symmetric or `abs()` and are unaffected.

---

<a id="p02"></a>
## P02 — `ratio_clover_ruby` tracks the ruby noise floor, not tension

| | |
|---|---|
| **Status** | open |
| **Severity** | S1 |
| **Opened** | 2026-08-07 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | the data itself, surfacing in `ratio_clover_ruby`; `pipeline/fret.py:20` (`background_subtract`) |
| **Affects** | `ratio_clover_ruby` — the primary tension output |
| **Evidence** | `docs/audit.md` §3 F2. Reproduce with `python Scripts/depth_check.py --dir 3d_images_folder/images/WT --channels` |

Recomputing the per-adhesion channel sums exactly as the pipeline does
(reproduces the shipped column at ρ = 1.0000):

| | ρ with `ratio_clover_ruby` |
|---|---|
| clover sum | **+0.068** |
| clover per voxel | +0.337 |
| ruby sum | −0.329 |
| **ruby per voxel** | **−0.773** |

A ratio is coupled to its denominator by construction, so *some* ruby
correlation is expected. The decisive figure is the other one: **the numerator
contributes essentially nothing (+0.068).** The exported tension proxy is, to a
good approximation, an inverse map of the ruby channel.

Why ruby carries no signal — measured inside adhesions after per-plane
background subtraction:

| channel | signal/voxel | plane background | contrast | voxels clipped to exactly 0 |
|---|---|---|---|---|
| mClover | 716.9 | ~553 | 129.6% | 0.0% |
| mRuby | 34.1 | ~499 | **6.8%** | **10.6%** |

`background_subtract` ends in `np.clip(img - mode_val, 0, None)`. On a 6.8%
contrast channel that discards a tenth of the measurements outright and
rectifies the noise on the rest — an asymmetric, signal-dependent bias on
exactly the quantity in the denominator. Its `bin_size=10` also quantises the
background estimate to ±5 counts, which is ~16% of a 34-count signal.

**The LUT bounds how much the ratio may legitimately vary.** The mRuby channel is
direct acceptor excitation, so `I_D/I_A ∝ (1 − E)` and the shipped calibration
(`Fret_LUT/Clover_mRuby2_GGSGGS7_force_new.txt`) sets a scale-free ceiling:

| force range | E spans | max `clover/ruby` spread |
|---|---|---|
| 2–5 pN | 0.400 → 0.218 | 1.30× |
| 1–6 pN | 0.538 → 0.181 | **1.77×** |
| 0.5–8 pN | 0.677 → 0.124 | 2.71× |

**Observed p1–p99 spread: 4.8×** (full range 22.6×). Even the robust figure, against
the most generous window, exceeds what the sensor can physically produce. Most of
the variance in the column is not tension.

**Fixing it** is mostly acquisition — mRuby2 needs more signal (longer exposure,
higher gain, better acceptor expression). On the code side `background_subtract`
is wrong for a channel this dim in two ways, but `pipeline/fret.py` is **frozen**
and the function is shared with the FRET path, so changing it needs an explicit
decision.

**Until then:** treat `ratio_clover_ruby` as a rank, not a calibrated quantity,
and be wary of the high tail — the highest ratios are the adhesions with the
least ruby.

**History**
- 2026-08-07 — recorded in `TODO.md` from `Scripts/depth_check.py --channels`;
  the ratio-vs-ruby correlation was noted and set aside as expected-by-construction.
- 2026-08-10 — audit established that the *numerator* correlation is ~0, which
  the earlier reading had not tested, and bounded the legitimate spread against
  the shipped LUT. Escalated from "the ratio is noisy" to "the ratio is
  substantially not measuring tension".

---

<a id="p03"></a>
## P03 — Adhesion detection is depth-truncated

| | |
|---|---|
| **Status** | open |
| **Severity** | S1 |
| **Opened** | 2026-08-10 (audit F3) |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/gui/adhesion_rf.py:138` (`gate_by_clover`, default on), `:163` |
| **Affects** | every column — the sample itself is biased |
| **Evidence** | `docs/audit.md` §3 F3 |

The adhesion RF runs on the clover channel, and prediction is gated by default on
a **global Otsu mask of that same channel**. Clover attenuates with depth, so a
global threshold accepts less of the volume the deeper you go — a selection
effect upstream of everything else.

Adhesion count normalised by vessel surface area per depth band (triangle areas
from the saved mesh):

| depth band | vessel surface (µm²) | adhesions | per 1000 µm² |
|---|---|---|---|
| 0–20 µm | 134,908 | 1402 | **10.39** |
| 20–40 µm | 46,667 | 98 | 2.10 |
| 40–60 µm | 28,928 | 21 | 0.73 |
| 60–80 µm | 36,480 | 106 | 2.91 |
| 80–105 µm | 128,890 | 465 | 3.61 |

Aggregated: near wall 8.26 per 1000 µm², far wall 3.05 — **the far wall yields
37% of the near wall's detection density** across comparable surface area
(181,575 vs 194,298 µm²).

Median adhesion volume falls only 11% over the same span (14.2 → 12.6 µm³), so
this is not adhesions shrinking with depth. It is adhesions going undetected.

**Why it matters:** any depth-stratified comparison contrasts a near-complete
shallow sample against roughly the brightest third of the deep one. That
truncation is itself a plausible generator of the drift in [P06](#p06) — the
surviving deep adhesions are a biased draw, not a fair one.

**Note on the earlier reasoning:** [P06](#p06) cited flat adhesion volume vs depth
(ρ −0.007) as evidence against a segmentation artefact. Volume is the wrong
statistic when the mechanism is truncation — truncation removes objects rather
than shrinking them.

**Next:** replicate across the other five images before any depth-corrected
analysis. If it holds, correcting the intensities will not be enough on its own.

---

<a id="p04"></a>
## P04 — Two nucleus label maps; `nuc_id` is not joinable to cell type

| | |
|---|---|
| **Status** | open |
| **Severity** | S1 |
| **Opened** | 2026-08-10 (audit F4) |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/gui/nuclei.py:211` vs `pipeline/gui/measure.py:57` |
| **Affects** | `nuc_id`, `dist_nuc_um`, `n_fibroblasts`, `n_nuclei_total` |
| **Evidence** | `docs/audit.md` §3 F4 |

The Nuclei tab builds two separate segmentations and only one reaches each
consumer:

- **Otsu + connected components** → `state["dapi_state"]["labels"]`. Not split,
  not size-filtered beyond ≥10 voxels.
- **Watershed-split nuclei** (`nuc`, `gui/nuclei.py:211`), built from the
  size-filtered binary. **Never written back to `dapi_state`.**

Measure reads `dapi_state["labels"]`, so `dist_nuc_um` and `nuc_id` come from the
unsplit map. `n_fibroblasts` / `n_nuclei_total` come from the split map.

Measured: the unsplit map has **362** components; `n_nuclei_total` is **273**;
CSV `nuc_id` values run to 355 and index the 362-label map. **The two cannot be
joined** — there is no way to ask whether a given adhesion's nucleus was
classified as a fibroblast, which is exactly the join this project needs.

Two consequences:

- **Measure ignores the `min_size` slider.** The user tunes it in the Nuclei tab;
  it drives the display and the fibroblast count and nothing else. The unsplit
  map's median component is **81.7 µm³**, far below a real nucleus (an oblong
  20 × 5 µm nucleus is ~260 µm³) — so more than half of what `dist_nuc_um`
  measures distance *to* is sub-nuclear debris.
- It blocks the statistics. `ratio ~ curvature + (1|image/cell)` needs cell
  identity; until this is fixed you are limited to `(1|image)` with n = 6. See
  [Analysis notes](#analysis-notes).

**The fix:** write the split labels back to `dapi_state`, or have Measure use the
split map. Either way, one nucleus definition, size-filtered once.

---

<a id="p05"></a>
## P05 — Fibroblast / endothelial classifier is unvalidated

| | |
|---|---|
| **Status** | open |
| **Severity** | S1 |
| **Opened** | 2026-08-06 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/gui/nuclei.py`, the "Count fibroblasts" step |
| **Affects** | `n_fibroblasts` |

Two cell types that are both DAPI-positive and distinguished only by which marker
they carry is the core unsolved problem in this project. The current classifier
has several compounding issues:

- It samples a **~2 µm shell around each nucleus** and asks whether mClover is
  above a threshold there. That perinuclear shell is precisely where *both*
  markers are dimmest — the sensor is out at the adhesions, VE-cadherin is out at
  the junctions — so the classifier reads the two cell types where they look most
  alike.
- **The result isn't stable.** The reported fibroblast fraction swings from 38%
  to 19% across a single step of the threshold slider. A number that halves on
  one slider click is not measuring the biology.
- **The shell isn't masked against neighbouring nuclei.** In a dense field the
  shell around nucleus A overlaps nucleus B and B's cytoplasm, so a fibroblast
  next to an endothelial cell contaminates both readings.
- The dilation is `ndi.binary_dilation(core, iterations=r)` with the default
  6-connectivity structure, so the "2 µm shell" is a diamond reaching ~1.95 µm
  along the axes and ~1.1 µm along the diagonals — not a ball.

**Fixing it** starts with [I02](#i02), not with more threshold tuning: the current
design has no validation signal at all, so there is no way to tell a fix from a
regression.

**Related:** [P04](#p04) (the count is not joinable to `nuc_id`), [P07](#p07)
(border nuclei lose their seed), [P13](#p13) (the cadherin exclusion multiplier
is hard-coded), [P16](#p16) (a claim from this entry that was retired).

**History**
- 2026-08-06 — recorded in `TODO.md`.
- 2026-08-10 — the "over-splitting inflates `n_nuclei_total`" bullet split out
  and closed as [P16](#p16) (invalid); the shell-shape and border-seed mechanisms
  added from the audit.

---

<a id="p06"></a>
## P06 — Tension readout drifts with imaging depth

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-07 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | an optical confound in `ratio_clover_ruby`, not a code defect |
| **Affects** | `ratio_clover_ruby`, and any relationship drawn against a spatial variable |
| **Evidence** | `docs/audit.md` §3 F2, §4. Reproduce with `Scripts/depth_check.py --channels` |

mClover (~515 nm) attenuates with depth; mRuby2 (~600 nm) does not measurably.
Depth quartiles on `V3_..._04step_002`:

| quartile | median ratio | clover/voxel | ruby/voxel |
|---|---|---|---|
| Q1 (shallowest) | 25.55 | 748.3 | 29.2 |
| Q2 | 18.01 | 653.9 | 34.4 |
| Q3 | 23.45 | 680.7 | 29.8 |
| Q4 (deepest) | 19.71 | 571.9 | 29.6 |

Overall ρ(ratio, `centroid_z_um`) = **−0.187**. Note the series is
**non-monotonic** — the effect is real and modest, not a clean gradient, and the
quartiles straddle a bimodal depth distribution so they mix the two walls.

The ruby channel's flatness has two readings and both are bad: either red
genuinely penetrates better, or ruby is too close to its noise floor
([P02](#p02)) to register attenuation at all.

**Why it matters, and how the picture changed.** The working assumption was that
depth on a tube is essentially *which wall you are on*, so stratifying by wall
would give two internal replicates spanning the same axial curvature. Measured,
that does not hold on this image: the two depth clusters share only **26%** of
the xy cells they occupy (Jaccard 0.26), and ρ(z, `dist_branch_um`) = **−0.461**.
The vessel dives through the stack, so **depth is substantially *where along the
vessel you are*.** Depth cannot be separated from axial position here, and axial
position determines curvature.

See [Analysis notes](#analysis-notes) for what follows, and [P03](#p03) for the
detection truncation that compounds it.

**Fixing it** needs, in order: (a) confirm it replicates across all six images;
(b) fix [P01](#p01) so curvature is trustworthy; (c) then correct for depth per
channel, restrict to a depth band, or carry depth as a covariate — but note (b)
and [P03](#p03) both have to land first for any of those to mean anything.

**History**
- 2026-08-07 — recorded in `TODO.md` as a −23% shallow-to-deep quartile drop.
- 2026-08-10 — reproduced, but the quartile series is non-monotonic and ρ is
  −0.187; the "split by wall" mitigation was tested against the data and does not
  hold on this image.

---

<a id="p07"></a>
## P07 — `peak_local_max` drops border seeds, in two places

| | |
|---|---|
| **Status** | fixed |
| **Severity** | S2 |
| **Opened** | 2026-08-10 (audit F5) |
| **Last verified** | 2026-08-10 (scikit-image 0.26.0) |
| **Resolved** | 2026-08-10 (commit `5d0a9c6`) |
| **Where** | `pipeline/gui/nuclei.py:207`, `pipeline/rf.py:606,613` |
| **Affects** | `n_nuclei_total`, `n_fibroblasts`, adhesion instance boundaries (**not** `nuc_id` / `dist_nuc_um` — see History) |
| **Evidence** | `docs/audit.md` §3 F5 |

`peak_local_max` defaults to `exclude_border=True`, which excludes peaks within
`min_distance` of the array edge. Verified on the installed scikit-image 0.26.0:
with `min_distance=5`, a peak at index 3 is dropped and reappears only with
`exclude_border=False`.

**(a) Nucleus splitting** — `gui/nuclei.py:207` uses `min_distance ≈ 12` voxels
(3.9 µm at 0.325 µm/voxel). Nuclei whose distance peak falls within 12 voxels of
any face of the volume get **no watershed seed**, and a mask component with no
marker is left at label 0 by `watershed(..., mask=binary)` — so those nuclei are
dropped entirely, or absorbed into a neighbour that does have a seed. On a
324-plane stack the excluded z shell is 7.4% of the depth, and it sits exactly
where the near-wall and far-wall adhesion populations are.

**(b) Adhesion splitting** — `rf.py:612` crops each candidate to its bounding box
with `pad = 2`, then calls `peak_local_max(dist, min_distance=5, ...)`. **The
border exclusion zone (5 voxels) is larger than the padding (2 voxels)**, so
seeds within 3 voxels of the component's own bounding box are suppressed. For
elongated adhesions the end seeds are the ones lost, so merged elongated objects
systematically fail to split.

**The fix:** pass `exclude_border=False` in both calls, and raise the adhesion
crop padding above `min_distance`.

**History**
- 2026-08-10 — opened from `audit.md` F5.
- 2026-08-10 — **fixed.** `exclude_border=False` added at `gui/nuclei.py:207` and
  `rf.py:613`, and `rf.py:606` `pad` raised `2 → 6` (above `min_distance=5`).
  Both parts of (b) are needed: the pad is clipped at the volume bounds, so
  padding alone cannot rescue components touching the stack edge. Measured on the
  shipped `V3_..._04step_002` intermediates:

  | | before | after |
  |---|---|---|
  | (a) nucleus seeds → `n_nuclei_total` | 273 → **273** | 276 → **276** |
  | (b) adhesion instances | **2,092** | **2,147** (+2.6%) |

  Directional only — there is no ground truth to check these against, so the new
  numbers are recorded, not validated. The post-split size re-filter at
  `rf.py:633-638` already absorbs the extra sub-threshold fragments the additional
  seeds produce.
- 2026-08-10 — **"Affects" corrected.** The original list named `nuc_id` and
  `dist_nuc_um`. Those come from the *unsplit* DAPI label map: the watershed
  result `nuc` is never written to `state` (`gui/nuclei.py` writes only
  `fibroblast_count`, `:248`), so the seeds that (a) restores never reach them.
  P07(a) affects `n_nuclei_total` and `n_fibroblasts` only; the `nuc_id` half of
  the original claim belongs to [P04](#p04).

---

<a id="p08"></a>
## P08 — `path_curv_xy`: 20 µm contamination radius, and κ is unsigned

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-06 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/centerline.py:128` (radius), `:139` (unsigned) |
| **Affects** | `path_curv_xy` |
| **Evidence** | `docs/audit.md` §3 F6 |

**(a) The contamination radius is `smooth_um`, not "at the junction".** κ is
fitted over `tree.query_ball_point(p, smooth_um)` — a **Euclidean** ball of 20 µm
by default over the centerline points. Near a junction that ball contains points
from two or three branches, so the parabola is fitted to a Y-shaped point set.
The same happens wherever the vessel doubles back within 20 µm of itself.

The `dist_branch_um` column is the tool for this, but **the threshold has to be at
least `smooth_um`, and should track the slider.** Measured, filtering
`dist_branch_um > 20 µm` costs 2.6% of rows. Nothing does it automatically.

**(b) κ is unsigned.** `centerline.py:139` returns `abs(2·c₂)/(1+c₁²)^1.5`, and κ
is looked up from the **centerline** — so every adhesion at a given axial
position receives the same κ regardless of which side of the vessel it is on.
`path_curv_xy` therefore cannot distinguish the inside of a bend from the
outside, arguably the most mechanically interesting contrast on a curved tube.
This is a design limitation rather than a bug, but it must be stated before the
column is used as *the* curvature variable.

Also worth knowing: the shipped centerline has median κ = 0.0040 1/µm (radius
250 µm) on a 332 µm field — **most of the vessel is nearly straight**, with the
signal concentrated in a few high-κ points. The audit attributed that
concentration to the junction; see (c), which is the larger contributor.

**(c) Border flare at the vessel exits dominates the κ signal.** The medial axis
of a shape with a flat cut end forks toward the corners of that cut. The vessel
exits through the xy walls of the volume, so the max-projected footprint is
truncated by a straight edge and `skeletonize` produces a fork at each exit.
`_prune_spurs` does not remove them — those forks run longer than
`prune_um = 30 µm`.

Measured on `V3_..._04step_002` (643 centerline points, 2 branch nodes):

| centerline points within… | count | median κ there | median κ interior |
|---|---|---|---|
| 5 µm of the volume edge | 17 (2.6%) | **0.0388** | 0.0038 |
| 10 µm | 28 (4.4%) | 0.0218 | 0.0038 |
| 20 µm | 65 (10.1%) | 0.0087 | 0.0036 |

**10× the interior curvature in the outer 5 µm, and 42% of all top-decile κ
points lie within 30 µm of an edge.** Both detected branch nodes are well
interior (51.7 and 143.6 µm from any edge), so this is not junction
contamination — it is an artefact of where the field of view was cropped.
`path_curv_xy` is the column with the strongest apparent correlation to the
tension proxy (ρ = +0.191), so a large share of its dynamic range is a
consequence of framing.

**The fix:** drop centerline points within `smooth_um` of the xy boundary, or
discard skeleton branches that touch the border before fitting κ. Costs ~10% of
centerline points.

**(d) The interim workaround silently empties a non-branched image.** Some
vessels in this dataset are a single tunnel or a cone rather than a Y. With no
branch node, `dist_branch_um` is left all-NaN — correctly and without crashing,
both paths guard with `if len(branch):` (`gui/measure.py:265`,
`measure_core.py:167`). But `dist_branch_um > 20` is False for NaN, so applying
the workaround below to such an image **drops every row** rather than none.

**Interim workaround:** filter `dist_branch_um > 20` yourself when analysing the
CSV — but guard it, e.g. `dist_branch_um.isna() | (dist_branch_um > 20)`, or the
non-branched images vanish. See (d).

**History**
- 2026-08-06 — recorded in `TODO.md` as "path curvature is undefined at branch
  nodes but still reported".
- 2026-08-10 — audit established the contamination radius (20 µm, not
  point-local) and added the unsigned-κ limitation.
- 2026-08-10 (later) — border flare (c) measured and found to be a larger
  contributor to the high-κ tail than the junction, revising the audit's reading.
  The NaN trap (d) added after the user noted that not every vessel is Y-shaped.

---

<a id="p09"></a>
## P09 — `dist_mesh_um` measures to rejected geometry; curvature has no distance gate

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-06 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/measure_core.py:103`, `pipeline/gui/measure.py:187` |
| **Affects** | `dist_mesh_um`, and the interpretation of every `curv_*` value |
| **Evidence** | `docs/audit.md` §3 F7 |

**(a) The two queries use different vertex sets.** `vessel_tree = KDTree(vessel_verts)`
uses **every** vertex; the curvature lookup correctly restricts to surface
vertices. For the default alpha-3D method `vessel_verts` is the whole subsampled
cadherin point cloud, so `filter_faces_by_extent` (`min_component_um = 50`) drops
debris islands' *faces* while their **vertices stay in the `dist_mesh_um` tree**.
Distance is measured to geometry the mesh step explicitly rejected.

Measured: 1,684 of 8,966 vertices (**18.8%**) belong to no face; 15.1% of
adhesions have one of them as their nearest vertex. The numerical effect is small
— median understatement 0.00 µm, p90 0.37 µm — so this is a semantics defect, not
a large error.

**(b) The resolution is coarser than the column suggests.** Surface vertex spacing
is **2.25 µm** (median nearest-neighbour). `dist_mesh_um` is a distance to a
point cloud with ~2 µm granularity, not to a surface.

> **Correction to `audit.md` F7.** The audit attributes that spacing to
> `subsample_um = 1.0`. The shipped settings JSON records **`subsample_um = 5.0`**
> (`..._settings.json` → `settings.vessel_mesh`), and rebuilding the point cloud
> at 5.0 reproduces the shipped mesh's vertex count exactly (8,966). The measured
> 2.25 µm spacing stands; the parameter attributed to it does not. `audit.md` is
> frozen and is not edited — see [P20](#p20) for what `subsample_um` actually
> controls, and re-read this entry's threshold discussion against a 5 µm grid
> rather than a 1 µm one.

**(c) No distance gate on curvature.** Every adhesion is assigned the curvature of
its nearest surface vertex with no maximum distance, so `curv_*` is populated for
every row regardless of how far the adhesion sits from the wall. Measured, this
is rarer than it sounds: **98.6% of adhesions are within 10 µm and only 2 rows
exceed 20 µm** (max 54.2 µm).

**Fixing (c)** is a one-line gate — NaN the curvature where
`dist_mesh_um > threshold` — but the threshold is a scientific choice, and it has
to be read against (b): it means something different for the alpha-3D method
(point cloud, ~2 µm spacing) than for the marching-cubes methods (spacing ≈ 1
voxel). **Needs a number from the user.**

**Interim workaround:** filter on `dist_mesh_um` yourself when analysing the CSV.

**History**
- 2026-08-06 — recorded in `TODO.md` as "curvature is assigned with no distance
  limit", with an example of an adhesion 30 µm from the wall.
- 2026-08-10 — measured: that case is 2 rows, not the bulk of the table. The
  rejected-geometry and vertex-spacing problems added.
- 2026-08-10 (later) — `audit.md` F7's `subsample_um = 1.0` attribution corrected
  to 5.0 from the shipped settings JSON; the (c) threshold must be read against
  that. Cross-linked to [P20](#p20).

---

<a id="p10"></a>
## P10 — `force_total_pN` is not a force

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-10 (audit F8) |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/measure_core.py:181`, `pipeline/gui/measure.py:280` |
| **Affects** | `force_total_pN` |

`force_total_pN = force_mean_pN × props["area"]`, and `props["area"]` is a volume
in µm³ (`regionprops_table` was given `spacing=scale`). The column's units are
pN·µm³.

Currently blank on this dataset because FRET is not configured ([P18](#p18)), so
it misleads no one yet — but it will the moment a `fret` block is added. Fix the
name, the units, or both before that happens.

---

<a id="p11"></a>
## P11 — `measure_core` duplicates the Measure worker and has diverged

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 |
| **Opened** | 2026-08-06 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/measure_core.py::compute_measurements` vs the `_run` worker in `pipeline/gui/measure.py` |

Two near-identical implementations of the same measurement. They have **already
diverged**: `measure_core` has no X-Y ROI filter, so the headless path measures
every adhesion in the volume where the GUI measures only those inside a drawn
region. Any future fix has to be made twice — the sum-then-divide ratio fix had
to be applied in three places for exactly this reason.

A headless CSV and a GUI CSV for the same image disagree whenever an ROI is
active. Both docstrings now say so, which is a patch on the symptom.

**Fixing it** means extracting one implementation both call — a restructure that
should not share a commit with behaviour changes. **Deliberately deferred.**

---

<a id="p12"></a>
## P12 — Headless measure fails silently and overwrites saved instances

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 |
| **Opened** | 2026-08-10 (audit F9) |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/batch_measure.py:143`, `:123` |

Two asymmetries between the GUI and headless paths:

- `gui/vessel_mesh.py:501` prints an explicit warning when a loaded curvature
  file's length does not match the mesh, because Measure will silently NaN every
  curvature column. `batch_measure.py` loads the same two files with **no such
  check**, so a stale pairing produces a CSV with four blank columns and no
  message.
- `_standard_instances` **overwrites `adhesion_instances.npz`** as a side effect
  of measuring. It is documented in the docstring and the reasoning is sound
  (many folders hold a stale, non-standard labelling), but a command named
  "measure" silently rewriting a saved segmentation deserves a flag rather than
  being the default.

---

<a id="p13"></a>
## P13 — Tunable parameters that are not sliders

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 |
| **Opened** | 2026-08-10 (audit F10) |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/gui/nuclei.py:241`, `pipeline/rf.py:612` |

The project convention is that tunable parameters become magicgui sliders so they
can be tuned live and are recorded in the settings JSON. Two thresholds that
materially change results are hard-coded:

- `is_fib &= cad_means < cbg * 1.3` — the VE-cadherin exclusion multiplier, on a
  channel that is dim in the perinuclear shell for **both** cell types.
- `min_distance=5` in the adhesion watershed split (1.6 µm at 0.325 µm/voxel).
  See also [P07](#p07)(b).

---

<a id="p14"></a>
## P14 — `cadherin_junction_proba` has no consumer

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 |
| **Opened** | 2026-08-06 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/gui/cadherin_rf.py` |

The VE-cad RF tab trains a classifier and produces a junction probability map into
`state["cadherin_junction_proba"]`, and then nothing reads it. No measurement, no
mesh step, no export consumes it. Its module docstring also points at a headless
`cadherin_rf_batch` script that **does not exist**.

The tab appears to be part of the pipeline and is currently a dead end. It is a
reasonable path to endothelial-cell segmentation — a clean junction network is
far more segmentable than thresholded VE-cad, and that bears directly on
[P05](#p05) — but the consumer has not been written.

---

<a id="p15"></a>
## P15 — Minor defects (five items)

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 |
| **Opened** | 2026-08-10 (audit F11) |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Evidence** | `docs/audit.md` §3 F11 |

Grouped because each is small and none affects the headline numbers.

- [ ] `centerline.py:41` — `_prune_spurs` mutates `skel` while iterating the
      endpoint list, but the neighbour set `S` was built before the mutation, so
      paths can be traced through already-deleted pixels within one pass. It
      converges, but the pruned result depends on iteration order.
- [ ] `measure_core.py:80` — adhesions thinner than one native z-step can fall
      entirely between sampled planes, giving `r_sum = 0` → NaN ratio. Correctly
      NaN'd rather than wrong, but it is a size-dependent (therefore
      depth-dependent) dropout that compounds [P03](#p03).
- [ ] `vessel_mesh.py:172` — `alpha_fill_slice` rasterises kept triangles in a
      Python loop over `skimage.draw.polygon`. Slow for alpha-2D on a large mask.
- [ ] `vessel_mesh.py:233` — `close_volume` uses a cubic structuring element, so
      the closing reaches √3× further along diagonals than the requested
      `close_um`.
- [ ] `denoise.py:107` — `np.convolve(..., mode="same")` with an even-length
      kernel introduces a half-bin phase shift in the polar radius smoothing.

---

<a id="p19"></a>
## P19 — `auto_min_size_otsu` is defeated by a hard-coded floor

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/denoise.py:61`, `pipeline/gui/cadherin.py:42,171,191` |
| **Affects** | the cadherin mask, and therefore every vessel-geometry column |

The Cadherin tab does auto-select the despeckle size cut — `auto_min_size_otsu(sizes)`
— but the returned value is then overridden by a module constant:

```python
MIN_SIZE_FLOOR = 12                                    # gui/cadherin.py:42
filter_controls.min_size.value = max(suggested, MIN_SIZE_FLOOR)   # :191
```

Measured on `V3_..._04step_002`, `auto_min_size_otsu` returns **3**, so the value
actually in use is always the constant. The comment at `gui/cadherin.py:38-42`
states the auto value is "~5"; it is 3.

So the operative parameter is a hand-tuned constant, and one that is **worse
placed than a slider**: it is invisible in the GUI and it is not written to
`run_settings`, so the exported settings JSON records the threshold parameters
but not the size cut that shaped the mask.

**Why Otsu fails here.** It splits `log10(sizes)` into two variance-minimising
classes. The size distribution is not two lumps in log space — it is a noise
power-law plus a broad real distribution:

| min component size | components kept | % of components | % of voxels |
|---|---|---|---|
| ≥1 vox | 226,093 | 100% | 100% |
| ≥5 vox (0.17 µm³) | 9,345 | 4.1% | 84.8% |
| ≥12 vox (**in use**) | 3,764 | 1.7% | 82.8% |
| ≥50 vox (1.72 µm³) | 1,616 | 0.7% | 80.1% |

95%+ of components are 1–4 voxels, and that spike dominates the histogram so
completely that Otsu's split lands *inside* the noise. The knee is sharp and
objectively locatable — the issue is the criterion, not the absence of one.

**Note:** whether `MIN_SIZE_FLOOR = 12` is right for the other five images has
never been checked. A single constant across images with different expression and
exposure is the assumption being made silently.

**Related:** [P13](#p13) (other hard-coded tunables), [I03](#i03) (scoring the
mesh rather than the size histogram is the better criterion).

---

<a id="p20"></a>
## P20 — `subsample_um` silently sets the vessel mesh's resolution

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-10 |
| **Resolved** | — |
| **Where** | `pipeline/vessel_mesh.py:47` (`voxel_subsample`), `pipeline/gui/vessel_mesh.py:90` |
| **Affects** | `dist_mesh_um`, all four `curv_*` columns, and the alpha value that is meaningful |

The shipped mesh used `subsample_um = 5.0`. That step turned **1,517,396 mask
voxels into 8,966 points — a 169:1 reduction** (rebuilding at 5.0 reproduces the
shipped vertex count exactly). Every downstream geometric property follows from
it, and none of this is documented:

- `dist_mesh_um` cannot resolve finer than the point lattice ([P09](#p09)(b)).
- `mesh_curvature` at `radius_um = 5.0` fits a **five**-parameter quadric to
  whatever lies within 5 µm — at 5 µm point spacing that is roughly 6–12
  neighbours, right at the rank-deficiency boundary the function guards against
  (`cell_mesh.py:383`).
- `alpha3d_um` must be ≳3× the point spacing simply to connect neighbours, so it
  is not principally tuning "how large a biological gap to bridge" — it is
  compensating for how sparse the subsample made the cloud. The two sliders are
  coupled and are presented as independent.

`voxel_subsample` also **discards the evidence it collects**. It quantises to a
grid and keeps one representative per occupied cell
(`np.unique(keys, axis=0, return_index=True)`), so a cell holding 2,244 wall
voxels and a cell holding 1 stray voxel become the same single point. Measured
occupancy on the 5 µm grid actually in use — median 69 voxels per cell, p10 = 4:

| votes required | cells kept | signal retained |
|---|---|---|
| ≥1 (**current**) | 8,966 (100%) | 100% |
| ≥5 | 7,909 (88.2%) | 99.8% |
| ≥10 | 7,318 (81.6%) | **99.6%** |
| ≥25 | 6,162 (68.7%) | 98.3% |

**Requiring ≥10 votes discards 18.4% of the points for 0.4% of the signal** — and
those are precisely the isolated points the alpha shape bridges to. Keeping the
count would give a noise gate whose failure mode is monotone in one threshold,
unlike the alpha shape's, which is combinatorial (any 4 points within a
circumradius create surface). It would not help against *correlated* error: a
debris clump or a neighbouring membrane votes as strongly as real wall.

**Related:** [P09](#p09)(b), [P21](#p21), [I03](#i03).

---

<a id="p21"></a>
## P21 — One global alpha cannot be correct on a tapering vessel

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-10 (measured on one image; the taper case is predicted, not yet measured) |
| **Resolved** | — |
| **Where** | `pipeline/cell_mesh.py:137` (`alpha_shape_faces`), `pipeline/gui/vessel_mesh.py:91` |
| **Affects** | what surface `dist_mesh_um` and `curv_*` refer to |

`alpha_shape_faces` returns the boundary of a **solid** — the faces bounding
exactly one surviving tetrahedron. Verified on the shipped mesh: **zero** edges
belong to exactly one face (98.4% belong to exactly two), one face-connected
component, and 1,684 of 8,966 vertices (18.8%) belong to no face because they are
strictly *inside* the solid. It is watertight, not a sheet draped over the wall.

Alpha must therefore satisfy two conditions at once: **above** the wall-gap scale
so gaps bridge, and **below** the local lumen radius so the lumen is not filled
with tets. The vessels in this dataset are not all uniform tubes — some are a
single tunnel, some a cone tapering to a point — so the lumen radius varies along
the axis and necessarily sweeps through any chosen alpha. **The narrow end fills
solid while the wide end stays hollow, for every possible alpha.** This is a
structural property of the method, not a tuning difficulty.

Two consequences:

- Where the lumen is hollow the wall is a slab whose boundary has **two sheets**,
  inner and outer. "Outward from the solid" points away from the axis on the outer
  sheet and *toward* it on the inner, so the same tube yields κ2 = −1/R on one and
  **+1/R** on the other. An adhesion receives whichever sheet its nearest vertex
  sits on. Where the lumen is filled there is one sheet and no such ambiguity — so
  a single image can contain both regimes.
- `dist_mesh_um` then measures to the outer envelope in one part of an image and
  to a filled rod in another.

**Status of the evidence.** Measured on `V3_..._04step_002` at the shipped
`alpha3d_um = 15.0`: the lumen is **not** filled (0% of centerline points at
z = 30–80 µm fall inside a surviving tet), so the two-sheet regime is live there.
Against that, `audit.md` F1 measured κ2 at **96.9% consistently signed** once
winding was preserved, which implies one sheet dominates. A partly-filled vessel
is consistent with both observations, but that is a hypothesis and has not been
verified. **Check it when [P01](#p01) lands** — split surface vertices by sheet
and inspect the κ2 sign histogram.

> **Answer, 2026-08-10 ([P01](#p01) landed).** With both parts of the P01 fix in,
> κ2 comes out **98.3% consistently signed** (median −0.0703 1/µm, R ≈ 14.2 µm)
> across all 7,282 surface vertices — up from 78.6% shipped and 77.1% with the
> winding fix alone. **One sheet dominates**, so the two-sheet ambiguity is not a
> material contributor to the curvature columns on this image. That does not make
> the structural argument above wrong: a tapering vessel must still sweep through
> any fixed alpha, so this is evidence about *this* image, not a general result.
>
> The inner/outer split itself **could not be measured directly.** P01 part B
> applies a single *global* orientation (flip the whole normal field by the sign
> of the enclosed volume), which is what makes the sign consistent — and which
> erases exactly the inner-vs-outer distinction this check wanted to see. A direct
> test would need per-sheet labelling (e.g. ray casting from the centerline)
> rather than reading it off the normals. Left open; the 98.3% figure bounds how
> much it can matter.

**Note on the acceptance test.** [P01](#p01)'s test reads "κ1 ≈ 0, κ2 ≈ −1/R
consistently signed, K ≈ 0". A cone is a *developable* surface, so K = 0 exactly,
as for a cylinder, and κ1 ≈ 0 still holds — but κ2 is **−1/R(s), varying along the
axis**. The test should check sign consistency and smooth variation, not agreement
with a single radius, or it will read as a failure on a taper image. K ≈ 0 applies
away from junctions only; a Y-junction is a saddle, K < 0.

The `implied caliber ~ 1/median κ1` readout (`gui/vessel_mesh.py:368`) is
meaningless on a tapering vessel regardless of [P01](#p01) — it averages a radius
that varies by construction.

**Related:** [P01](#p01), [P20](#p20), [I03](#i03).

---

<a id="p22"></a>
## P22 — Centerline overlay is drawn at a single wrong z

| | |
|---|---|
| **Status** | fixed |
| **Severity** | S3 (display only — no exported column changes) |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-10 |
| **Resolved** | 2026-08-10 (commit `5d0a9c6`) |
| **Where** | `pipeline/gui/vessel_mesh.py:397` (`_path_z_um`, `_show_path`) |
| **Affects** | the `Vessel path curvature (X-Y)` layer only |

`_show_path` places every centerline point at one z — the mean of **all** mesh
vertices — so the overlay is a flat sheet at an arbitrary depth rather than
anything tracking the vessel.

Measured on `V3_..._04step_002` (display z = 46.7 µm), taking the mean z of the
cadherin column beneath each centerline point:

| | |
|---|---|
| true per-column wall z along the path | **3.1 → 99.8 µm** (the stack is 105 µm) |
| median \|true z − displayed z\| | **32.8 µm** |
| p90 / max | 48.7 / 53.2 µm |

The vessel traverses essentially the whole depth ([P06](#p06): it dives through
the stack), so the overlay sits on average a third of the stack away from it.
Independently confirmed by the [P21](#p21) measurement: 0% of centerline points at
the mean z fall inside the mesh solid.

**The fix** is small and changes no exported number: give each centerline point the
mean z of the cadherin column beneath it (defined for 574/643 points here; a small
neighbourhood average covers the rest). Worth doing before any 3D-centerline work,
because it shows whether a 3D axis estimate is well-behaved before one is built.

**Not a defect in `path_curv_xy` itself** — κ is computed in xy and z never enters
it. This is the overlay only.

**History**
- 2026-08-10 — opened.
- 2026-08-10 — **fixed.** New helper `_path_z_um` gives each centerline point the
  intensity-weighted mean z of the cadherin voxels in a ±3 µm column beneath it,
  sampled from `_cadh_mask_array(state)` — the same array `_compute_path` fed the
  centerline, so the two agree by construction. Only the ~643 needed columns are
  touched; `_show_path` runs on the GUI thread and must not materialise the
  volume. Measured on `V3_..._04step_002`:

  | | before | after |
  |---|---|---|
  | displayed z | one value, **46.7 µm** | per point, **6.1 → 96.1 µm** (median 38.9, σ 25.2) |
  | \|new − old\| | — | median **23.5 µm**, max **49.4 µm** |
  | columns with cadherin signal | — | **643 / 643** |
  | cost on the GUI thread | — | **0.15 s** |

  The measured span is narrower than the 3.1 → 99.8 µm quoted above because the
  ±3 µm window averages rather than taking a single column — intended, and it is
  what makes all 643 points resolvable (the entry's single-column figure covered
  only 574).

  **Also fixed while here** (not in the original entry): `zc` fell back to `0.0`
  when no vessel mesh was loaded, silently drawing the whole path on the bottom
  slice. Points with an empty column now fall back to the mask's own mean z over
  the points that did hit; only a completely absent cadherin mask can still reach
  the mesh-mean → `0.0` chain.

---

<a id="p23"></a>
## P23 — Centerline preprocessing amputates tapering vessels

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-10 (predicted from the code and the stated geometry; not yet measured on a taper image) |
| **Resolved** | — |
| **Where** | `pipeline/centerline.py:109` (`open_um`), `:118` (`prune_um`) |
| **Affects** | `path_curv_xy`, `dist_branch_um` on non-uniform vessels |

Some vessels in this dataset taper to a fine point rather than running as a
uniform tube. `vessel_path_curvature_xy` removes such a tip twice over:

- `binary_opening(foot, disk(open_um / ppm))` with `open_um = 4.0 µm` deletes any
  part of the footprint narrower than the structuring element, so the taper is cut
  off once its radius falls below ~4 µm.
- `_prune_spurs(..., prune_um = 30.0 µm)` then removes terminal branches shorter
  than 30 µm — which is what a surviving tip looks like to it.

Both are **function defaults with no slider** and neither scales with the vessel.
The tip is where the geometry is most extreme (tightest caliber, highest
curvature), so this removes the most informative part of exactly the images where
it exists.

**Not yet measured** — no taper image has saved intermediates. Quantify it on one
before choosing values, and consider making both thresholds scale-aware rather
than absolute.

**Related:** [P08](#p08)(c) (border flare, same function), [P13](#p13).

---

# Closed

<a id="p16"></a>
## P16 — Nucleus splitting inflates `n_nuclei_total`

| | |
|---|---|
| **Status** | **invalid** |
| **Opened** | 2026-08-06 (as a bullet within [P05](#p05)) |
| **Closed** | 2026-08-10 |
| **Evidence** | `docs/audit.md` §3 F4 |

**The claim:** nucleus splitting uses a ~3.9 µm minimum separation while these
nuclei are oblong and up to ~20 µm long, so one nucleus is broken into several,
inflating `n_nuclei_total` and the denominator of every fraction.

**Why it is closed:** measured on `V3_..._04step_002`, splitting **reduces** the
count — 362 unsplit connected components against `n_nuclei_total` = 273. The size
filter and the border loss in [P07](#p07) remove more than the splitting adds, so
the stated *direction of the effect on the count* is wrong.

The per-nucleus over-splitting concern may still be real, but it is not visible in
the count and there is currently no evidence for it. Recorded so it is not
rediscovered as fact.

---

<a id="p17"></a>
## P17 — `curv_gaussian` is unaffected by the curvature bug

| | |
|---|---|
| **Status** | **invalid** |
| **Opened** | 2026-08-06 (stated in `TODO.md` and `README.md`) |
| **Closed** | 2026-08-10 |
| **Evidence** | `docs/audit.md` §3 F1 |

**The claim:** `curv_gaussian` (K = κ1·κ2) survives the curvature defect because
both signs flip and the product does not, so it can be used where `curv_mean`
cannot.

**Why it is closed:** the reasoning holds for a *global sign flip*, but the actual
defect ([P01](#p01)) is a corrupted fitting **frame**, not a sign — and a product
does not cancel that. Measured ρ(shipped, correctly-wound) for `curv_gaussian` is
**+0.187**, essentially as destroyed as `curv_mean` (+0.048).

A second reason not to reach for it even once [P01](#p01) is fixed: on a tube
κ2 ≈ 0, so K ≈ 0 everywhere by construction. Measured median is +0.0000. The
sign-safe curvature is the one with almost no dynamic range on this geometry.

---

<a id="p18"></a>
## P18 — FRET constants set but not calibrated for this dataset

| | |
|---|---|
| **Status** | **fixed** |
| **Opened** | before 2026-08-06 |
| **Resolved** | 2026-08-06, commit `f8106a5` |
| **Verified** | 2026-08-10 — `force_mean_pN`, `force_std_pN`, `force_total_pN` are all-NaN in the shipped CSV and export blank |
| **Where** | `pipeline/io.py::fret_configured` |

The `fret` block in `params.json` carried values that were **set but not
calibrated for these images** — the dangerous state, because everything
downstream runs and produces plausible force numbers.

**Fix:** the block was removed from `params.json`, and FRET now runs only when
that block is present **and** the `clover`/`ruby`/`fret` roles are all mapped.
Behaviour is driven by presence, not by a boolean flag that could default wrongly,
and the loader logs which path it took either way.

A related claim from the same era — *"the FRET constants are unset"* — was simply
wrong; they were set. Recorded here so the distinction survives.

**Note:** [P10](#p10) is latent behind this gate and will surface the moment a
`fret` block is configured.

---

# Ideas

Not defects. These do not exist yet, so they cannot mislead.

<a id="i01"></a>
## I01 — Progress / completion bars for long-running operations

**Opened** 2026-08-06 · **Status** open

When a tab's compute button is clicked (Nuclei → "Compute DAPI Otsu +
components", Adhesion RF → "Run RF predict"), show how far along it is and when
it is done, instead of only terminal prints.

**Not all steps can show a true percentage:**

- **Determinate (real %)** for steps that loop over items: `rf.predict_proba`
  already runs in N z-chunks and prints `chunk i/N`; load-time resampling loops
  over channels (`io.to_isotropic`) and FRET outputs
  (`fret.resample_fret_outputs`); FRET compute loops per-Z. These need the
  blocking worker converted to a generator that `yield`s progress, and
  `rf.predict_proba` would need a progress callback threaded out of it.
- **Indeterminate only** for atomic single library calls with no sub-progress
  hook: `otsu_mask`, `label_and_sizes`, `niblack_mask`, watershed splitting. The
  honest option is a "working…" bar that animates and clears on completion.

**Options:** inline `magicgui.ProgressBar` per tab (matches the request most
closely, more wiring), or napari's built-in activity dock
(`napari.utils.progress` — much less code and idiomatic, but off to the side).

**Open decisions:** inline vs activity dock; scope (all tabs, just the heavy
steps, or a Nuclei-tab prototype first); whether atomic steps get a busy bar or
none.

**Note:** `pipeline/gui/logpanel.py` implements `isatty`/`fileno`/`encoding` on
its stream, so a tqdm-style bar routed through stdout will behave sanely in the
Log dock rather than emitting one line per update.

---

<a id="i02"></a>
## I02 — Fibroblast / endothelial scatter diagnostic

**Opened** 2026-08-06 · **Status** open · **Blocks** [P05](#p05)

A validation signal for [P05](#p05), before touching the classifier itself.

One point per cell, plotted as **per-cell sensor signal (mClover)** against
**per-cell-edge VE-cadherin signal** — measuring each marker where it actually
lives (the sensor through the cell body and adhesions, cadherin at the cell
boundary) rather than both in the same dim perinuclear shell. Score each cell
against the other cells **in the same image**, so illumination and expression
differences between images do not move the decision boundary.

**Why it comes first:** if the two populations separate into two clouds, the
classification problem is real, the current features are simply wrong, and the
scatter shows where the boundary should be. If they do not separate, no amount of
threshold tuning will work and the approach has to change (e.g. a
segmentation-based assignment rather than a per-nucleus shell). Either outcome is
more informative than the current situation, where there is nothing to check a
change against.

---

<a id="i03"></a>
## I03 — Objective-scored auto-calibration of the vessel mesh

**Opened** 2026-08-10 · **Status** open · **Relates to** [P19](#p19),
[P20](#p20), [P21](#p21)

The vessel surface currently depends on four coupled knobs — `min_size`
(a hidden constant, [P19](#p19)), `subsample_um` (dominant and undocumented,
[P20](#p20)), `alpha3d_um` (structurally unsatisfiable on a taper,
[P21](#p21)), and `min_component_um` — tuned by eye per image. The goal is to
replace that with a measurement.

**Score the mesh, not the mask.** [P19](#p19) shows why: a component-size
histogram does not know what the mask is *for*, which is how Otsu on it comes to
need a magic constant to rescue it. Candidate objectives, none needing a human
eye:

- **Topology.** The strongest, because the answer is known a priori. The shipped
  mesh is one component of genus ≈ 125 (χ = −248); a real vessel wall shell is
  single-digit genus. **State it as minimise-genus subject to a coverage floor,
  not as a target value** — expected genus differs by image (roughly 0 for a cone
  closing to a point, 1 for a single tunnel, 2 for a Y), so a fixed target would
  reintroduce exactly the hidden constant this is meant to remove. Genus alone is
  gameable by collapsing to a blob; the coverage term is what prevents that.
- **Coverage.** Fraction of retained cadherin voxels within *d* of the surface.
- **Parsimony.** Surface area per unit covered signal — penalises grabbing noise.
- **Stability.** Rebuild from a random half of the components; a good setting
  reproduces, one tuned to noise does not.

**Method selection can be automatic too.** Branch count is already computed by
`vessel_path_curvature_xy`. With zero branches, `polar_fill_slice`'s
one-star-shaped-lumen-per-slice assumption holds and it has no gap-bridging knob
at all (it interpolates empty angular bins by construction) — so on the
non-branched images the cheapest and most noise-robust method already in the repo
applies. The Y-focused module docstrings obscure this.

**Cost.** An alpha sweep ran 10–25 s per point on `V3_..._04step_002`, so a
coarse grid over three parameters is minutes. Acceptable for a "calibrate once per
image" button; not for live slider feedback.

**Sequencing.** After [P01](#p01), so that corrected curvature can be part of the
objective rather than something validated afterwards. This and [P21](#p21) are one
piece of work — a reconstruction method and its calibration — and should not share
a commit with [P01](#p01).

---

# Housekeeping

<a id="h01"></a>
## H01 — Repository has no LICENSE

**Opened** 2026-08-10 · **Status** open

The repository is intended to be public and carries no `LICENSE` file. By default
that means all rights reserved: nobody may legally use, modify, or redistribute
the code, which is almost certainly not the intent for a lab tool.

Deliberately left open pending a choice. MIT is the usual default for academic
lab tooling; BSD 3-Clause adds a no-endorsement clause that some universities
prefer; Apache 2.0 adds an express patent grant, which matters if the
tension-measurement method might be patented or if industry collaborators are
involved.

Also noted at the same time: `3d_images_folder/Channel Info.rtf` was removed from
the working tree on 2026-08-10 (lab acquisition notes wrapped in an RTF carrying
embedded authoring metadata; its useful content — the channel index list — is in
`params.json` and `CLAUDE.md`). **It remains in git history.** Scrubbing history
would need `git filter-repo` and a force-push, which is a separate, deliberate
decision.

---

<a id="analysis-notes"></a>
# Analysis notes

Not defects — standing cautions that apply whenever the exported CSV is analysed.
The full reasoning and the measurements behind these are in `docs/audit.md` §4.

## The unit of analysis is not the adhesion

There are many adhesions per cell, several cells per image, and **six images, one
condition**. Curvature is a property of the *vessel*, not of the adhesion: two
adhesions 5 µm apart share it. The number of independent curvature observations
per image is roughly (vessel length ÷ curvature correlation length) — tens, not
thousands.

A pooled Spearman over 2,092 adhesions therefore reports adhesion count, not
evidence. It will be tiny for effects that are not real. **The
`Spearman(ratio, …)` lines the Measure tab prints are a quick read for one image,
not a result.**

The cell level is currently unavailable — see [P04](#p04). Until that is fixed you
have `(1|image)` with n = 6, or a spatial clustering of adhesions as a rough cell
proxy.

## Everything spatial correlates, because it is all one gradient

Measured on `V3_..._04step_002`:

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
columns sit at |ρ| ≤ 0.04 — exactly what [P01](#p01) predicts. Stratifying by wall
does not resolve this; it amplifies it (far wall: `dist_branch_um` +0.432,
`centroid_z_um` −0.330, `path_curv_xy` +0.292 — everything correlating with
everything, the signature of one underlying gradient).

## Use a spatial null, not a permutation null

Shuffling tension values across adhesions destroys spatial autocorrelation and
produces a null that is far too narrow — it will manufacture significance.
Instead, circularly shift the tension field against the curvature field along
centerline arclength within an image, and rebuild the statistic over many shifts.
That preserves both fields' autocorrelation and destroys only their alignment.
Given the correlation structure above, this is not optional.

## Run a falsification battery, decided in advance

- **clover alone vs curvature** — if it reproduces the ratio's effect, the ratio
  contributed nothing.
- **ruby alone vs curvature** — if ruby drives it, the effect is denominator
  noise ([P02](#p02)).
- **`volume_um3` vs curvature** — catches a segmentation-density artefact.
- **cross-image scramble** — adhesions from image A against image B's curvature
  field. Must be null. If it is not, the analysis is broken.

## Pre-specify the exclusions

Before looking at tension, fix and write down: the `dist_mesh_um` gate
([P09](#p09)), `dist_branch_um > smooth_um` ([P08](#p08)), a margin from the xy
field boundary ([P08](#p08)(c) — 42% of top-decile κ is border flare), the depth
band, and whether `09step` (0.9 µm z-spacing) is in. Each is otherwise a fork in
the garden.

**Guard the branch filter.** Not every vessel here is Y-shaped — some are a single
tunnel or a cone. Those images have no branch node, so `dist_branch_um` is all-NaN
and a bare `dist_branch_um > 20` drops **every** row rather than none. Write it as
`dist_branch_um.isna() | (dist_branch_um > 20)`. See [P08](#p08)(d).

## State the ceiling

With six images, if all agree in direction a sign test floors at p = 0.031. That
is the honest best case.

- **Defensible, eventually:** "In all six images, adhesion tension rank increases
  with axial vessel curvature."
- **Not defensible, and will not be until FRET is calibrated:** any slope in pN
  per unit curvature.
