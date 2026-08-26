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
| [P01](#p01) | Curvature columns are noise — `alpha_shape_faces` destroys winding | S1 | **fixed** — part B has a principled replacement for [I06](#i06) |
| [P02](#p02) | `ratio_clover_ruby` tracks the ruby noise floor, not tension | S1 | open |
| [P03](#p03) | Adhesion detection is depth-truncated | S1 | open — **confirmed 6/6**; 2026-08-25 the clover gate is ruled out as the cause ([P43](#p43)); live hypothesis is annotation sampling |
| [P04](#p04) | Two nucleus label maps; `nuc_id` is not joinable to cell type | S1 | **fixed** — but joinable is not correct: see [P42](#p42) |
| [P05](#p05) | Fibroblast / endothelial classifier is unvalidated | S1 | open — **2026-08-25: likely over-calls ~10×**; all three routes to the count now failed ([P41](#p41), [P42](#p42), §8E) |
| [P06](#p06) | Tension readout drifts with imaging depth | S2 | open — **sign flips 6 images** |
| [P07](#p07) | `peak_local_max` drops border seeds, in two places | S2 | **fixed** |
| [P08](#p08) | `path_curv_xy` — 20 µm fit window on a ~140 µm-wide vessel; contamination; unsigned | S1 | **superseded** — the 2D centerline is deleted (2026-08-22) |
| [P09](#p09) | `dist_mesh_um` measures to rejected geometry; no gate on curvature | S2 | **(a) fixed**, (b)(c) open |
| [P10](#p10) | `force_total_pN` is not a force | S2 | open |
| [P11](#p11) | `measure_core` duplicates the Measure worker and has diverged | S3 | open |
| [P12](#p12) | Headless measure fails silently and overwrites saved instances | S3 | **fixed** |
| [P13](#p13) | Tunable parameters that are not sliders | S3 | open |
| [P14](#p14) | `cadherin_junction_proba` has no consumer | S3 | open |
| [P15](#p15) | Minor defects (five items) | S3 | open |
| [P16](#p16) | Nucleus splitting inflates `n_nuclei_total` | — | **invalid** |
| [P17](#p17) | `curv_gaussian` is unaffected by the curvature bug | — | **invalid** |
| [P18](#p18) | FRET constants set but not calibrated | — | **fixed** |
| [P19](#p19) | `auto_min_size_otsu` is defeated by a hard-coded floor | **S2** | open — **raised 2026-08-21: `min_size` deletes 59–85% of [P28](#p28)'s dim wall** |
| [P20](#p20) | `subsample_um` silently sets the vessel mesh's resolution | S2 | open |
| [P21](#p21) | One global alpha cannot be correct on a tapering vessel | S2 | open — **does not apply to [I06](#i06)'s method** |
| [P22](#p22) | Centerline overlay is drawn at a single wrong z | S3 | **fixed** — Z rule replaced 2026-08-13 |
| [P23](#p23) | Centerline preprocessing amputates tapering vessels | S2 | **superseded** — `open_um`/`prune_um` deleted; the graph prunes by caliber ratio |
| [P24](#p24) | `curv_*` has no scale-independent value — `radius_um` sets the answer | S1 | open — **plateau reached 2026-08-13**, CSVs not rebuilt |
| [P25](#p25) | `fibroblast_count` never cleared — a failed image inherits the previous one's counts | S1 | **fixed** |
| [P26](#p26) | Batch builds vessel meshes at widget defaults; the five new images do not match the reference | S2 | open |
| [P27](#p27) | The centerline is a point set in raster order, not a traced path | S3 | **fixed** — the graph carries arc length by construction |
| [P28](#p28) | The cadherin threshold discards real wall — the gaps have signal in them | S2 | open — **downgraded 2026-08-24: nice-to-have for reconstruction, see History** — replicated 5/6 2026-08-21; **the local threshold is innocent, the floor + `min_size` are not**; the continuous-field fix is measured and **does not work** |
| [P29](#p29) | A field-of-view notch fabricates a branch node; `dist_branch_um` is wrong for 18% of rows | S2 | **superseded** — no projection, no `fill_holes`; image 1 now reports its one junction |
| [P30](#p30) | `build_mesh` silently treats an unknown method string as alpha 3D | S3 | open |
| [P31](#p31) | `mesh_curvature` has no boundary gate — curvature biased near an open edge | S2 | open — **7.1% of vertices at `radius_um` 5.0** |
| [I01](#i01) | Progress / completion bars for long-running operations | idea | open |
| [I02](#i02) | Fibroblast / endothelial scatter diagnostic | idea | open |
| [I03](#i03) | Objective-scored auto-calibration of the vessel mesh | idea | open — **scorer landed** (`mesh_score`); choosing still open; **2026-08-21: `slope_loglog` alone was gamed by a worthless mesh** |
| [I04](#i04) | `batch.py` has no way to re-run one step | idea | **done** 2026-08-24 — `batch.py --steps` |
| [I05](#i05) | Is the alpha shape the right reconstruction method at all? | idea | open |
| [I06](#i06) | Oriented-point implicit reconstruction — the plan of record | idea | **in progress — `alpha 3D solid -> poisson` is the direction of record as of 2026-08-21**, visually confirmed; three refinements open (bumps, runtime, Open3D) |
| [I07](#i07) | Local thresholding on the clover channel | idea | **closed 2026-08-25 — premise refuted.** The gate was never binding; it is [P43](#p43)-removed, and [P03](#p03) is unaffected |
| [I08](#i08) | "Refine with cadherin" is a dead end and should be removed | idea | open — **superseded by [I06](#i06); sole consumer of [P34](#p34)'s centroids** |
| [P32](#p32) | The fitted lumen leaks through the z faces on images 4 and 5, wrecking the graph topology | S2 | open — **diagnosed and warned, not fixed** |
| [P33](#p33) | `sauvola_mask` cast to float32 on the host, doubling the PCIe transfer | S3 | **fixed** — now at parity with niblack, bit-identical |
| [P34](#p34) | Cadherin tab computed per-component centroids eagerly for an off-path consumer | S3 | **fixed** — 34.45 s -> 0 at `floor_k` 3.0 |
| [P35](#p35) | `alpha 3D` hangs forever on a permissive mask, with no warning | S3 | **mitigated** — guard landed 2026-08-24; it refuses instead of hanging |
| [P36](#p36) | The cadherin "local threshold" is inert: the mask is a global threshold | **S2** | open — Sauvola inert; **the forcing override is fixed 2026-08-24**, Niblack is now the default everywhere |
| [P37](#p37) | The loader overwrote `floor_k` on every load, so the tab default never ran | S2 | **fixed** 2026-08-24 |
| [P38](#p38) | A failed mesh build left the previous image's lumen field; the graph used it | **S1** | **fixed** 2026-08-24 — one saved graph was from the wrong image |
| [P39](#p39) | 7 days may be short of a confluent monolayer, so some cadherin gaps may be real | S3 | open — hypothesis, untested |
| [P40](#p40) | `fibroblast_cellpose.py` thresholded noise as tissue; wall split degenerated silently | S2 | **fixed** 2026-08-25 — sanity-gate half still open |
| [P41](#p41) | The en-face VE-cadherin wall has no polygonal cells; endothelial segmentation returns nothing | S2 | open — **blocks the subtraction route to [P05](#p05)** |
| [P42](#p42) | `nuc_id` is a nearest-neighbour join presented as ownership | **S1** | open — **a cross-tab was retracted on this**; qualifies [P04](#p04) |
| [P43](#p43) | The clover gate cost bias and saved ~10%; removed | S3 | **fixed** 2026-08-25 — refutes [I07](#i07)'s premise |
| [I09](#i09) | RF inference is dominated by `_eigh3x3` on the CPU (~52%) | idea | open — **~2× available, exact**; profiled, three proposals rejected on measurement |
| [P44](#p44) | The RF's training annotations were never saved; the model cannot be retrained | **S2** | open — **blocks [P03](#p03)'s hypothesis and its fix**; escalates [H02](#h02) |
| [P45](#p45) | `batch_mesh` crashed on every invocation — two names used, never imported | S3 | **fixed** 2026-08-26 — every mesh artifact since `48b4f9d` came from the GUI |
| [P46](#p46) | A digitized wall is longer than the wall it stands for, so absolute τ reads ~5% low | S3 | open — **quantified, common-mode**; cancels in the within-segment profile, not in Pa |
| [H01](#h01) | Repository has no LICENSE | — | open |
| [H02](#h02) | `adhesion_rf.joblib` is a pickle against an unpinned scikit-learn | — | open |

**Suggested order of work** (dependency-ordered, not severity-ordered). The
original ordering was **P01 → P04 → P07 → P03 → P02 → P08/P09/P10**, then the
analysis in [Analysis notes](#analysis-notes) — P01 first because it is the
smallest fix with the largest effect and every curvature result depends on it.

**As of 2026-08-10 P01, P07 and P22 are fixed in the code**, so the queue is now
**P04 → P03 → P02 → P08/P09/P10**. One caveat carries over: on 2026-08-11 only
`V3_..._04step_002` was regenerated (`..._20260811_105032.csv`). It is the one
image that had `_saved` intermediates; the other five have none, so **any older
CSV, and any CSV for another image, still carries the pre-P01 curvature noise**.
Fixing the code and fixing the data are separate steps, and the data half is one
sixth done.

The 2026-08-10 vessel-mesh review (P19–P23, [I03](#i03)) does not displace that
ordering. The border-flare half of [P08](#p08) is independent of everything above
and can land at any time. [P21](#p21) and [I03](#i03) are one piece of work —
reconstruction method plus its calibration — and are now unblocked: corrected
curvature can be part of the calibration objective rather than something checked
afterwards.

**Revision, 2026-08-11 — the geometry work moves up.** [P24](#p24) and
[P08](#p08)(e) establish that **both** curvature columns are set by their own fit
scale rather than by the vessel: `curv_*`'s implied radius tracks `radius_um` with
no plateau, and `path_curv_xy` fits a 20 µm window to a structure ~140 µm wide.
Neither is a defect that a small fix reaches — both follow from what the
reconstruction is. That makes [P21](#p21) + [I03](#i03) the critical path for the
project's end goal rather than a follow-up to it, because there is currently no
trustworthy geometry variable to relate tension to.

The revised queue is **P04 → P03 → P21/I03 (+ P24) → P02 → P08/P09/P10**, with the
five un-batched images ([P03](#p03), [P06](#p06), [P19](#p19), [P23](#p23) all
need them) running unattended alongside. [P04](#p04) stays first: it is small,
it is the only route to a cell-level unit of analysis, and it is independent of
everything the geometry work will change.

**Update, 2026-08-11 — P04 and P12 are fixed; the five-image batch is now the
critical path.** Both were done *before* the batch deliberately: `batch.py` has no
`--steps` flag and only instances and measurement can be redone headlessly
(`batch_measure`, `batch_postprocess`), so anything that changes a **saved
artifact** had to land first. [P04](#p04) adds one (`nuclei_split.npz`);
[P12](#p12) makes an unattended run report what it degraded instead of writing a
plausible-looking CSV in silence; [P25](#p25) — found while tracing P04 — stops
one image's nucleus counts being written into another's file.

The queue is now **the 5-image batch → P03 → P21/I03 (+ P24) → P02 →
P08/P09/P10**. Producing those five `_saved` folders is the single thing blocking
[P03](#p03), [P06](#p06), [P19](#p19), [P21](#p21), [P23](#p23) and the whole
six-image design in [Analysis notes](#analysis-notes).

> **Two operational cautions for that run.**
>
> 1. **Use `--skip-existing`.** The WT folder holds six tifs and one `_saved`.
>    Without the flag the batch reprocesses `V3_..._04step_002` and overwrites
>    the only post-[P01](#p01) reference folder and its CSV — the baseline every
>    before/after in this file is measured against.
> 2. **That leaves image 6 without `nuclei_split.npz`.** Regenerate it afterwards
>    in the GUI (Load → Load all saved → Count fibroblasts → Save), which is
>    cheap and skips the RF predict. Until then that image measures through the
>    documented unsplit fallback and says so.
>
> Run `python -m pipeline.batch_measure --dir <dir> --check-only` after the batch
> and before measuring: it reports which columns would come out blank in about a
> second per folder.

**Update, 2026-08-11 (later) — the batch has run. All six images now have
`_saved` intermediates.** 5 of 5 processed, 0 skipped, no step failures;
`--check-only` reports **0 warnings across all six folders**. `nuclei_split.npz`
was regenerated for the reference image afterwards, as planned, so the P04 join
is available for every image.

Two things came out of it, and they pull in opposite directions:

- **[P05](#p05) now has cross-image evidence for the first time** — the
  fibroblast fraction spans 5.0%–14.9% (3.0×) across six images of one
  condition, the same order as the within-image slider swing. The two cannot
  currently be separated.
- **[P26](#p26) — new.** The batch built every mesh at *widget defaults*, but the
  reference image's mesh was hand-tuned (`subsample_um` 5.0 vs the default 1.0,
  and three more). The five new meshes carry 23–30× more vertices, so
  `dist_mesh_um` and `curv_*` are **not comparable between the reference and the
  rest**. Deliberately left open: there is no known-good parameter set to
  standardise on, which is [I03](#i03)'s entire premise.

**No CSVs have been exported yet, on purpose.** Measuring now would bake the
default-mesh geometry into six exports that would then need redoing, and this
project already has three CSVs in one folder causing exactly that confusion
(see [P01](#p01)'s note). Run `batch_measure` once [P26](#p26) is decided. The
non-geometry columns are unaffected either way, so an interim measure is
defensible if the tension and nucleus columns are wanted sooner — just expect to
re-export.

> **Overtaken by events, 2026-08-12.** All six *have* now been exported
> (`20260812_011*`), taking the "interim measure" option in that last sentence.
> The trigger was [P09](#p09)(a), which needed a re-measure to land and which
> made the tension and nucleus columns available for all six as a side effect.
> The caution still holds for what it was about: **the geometry columns in those
> six exports will need redoing** once [I03](#i03) settles the reconstruction.

So the queue is now **P26 (+ the headless mesh rebuild it needs, [I04](#i04)) →
P03 → P21/I03 (+ P24) → P02 → P08/P09/P10**. [P03](#p03), [P06](#p06),
[P19](#p19) and [P23](#p23) are unblocked *today* on everything that does not
depend on the mesh — the five images they were waiting for now exist.

**Update, 2026-08-12 — [P09](#p09)(a) jumps the queue; the reference mesh is
demoted.** Two things landed on the same day and they reorder the front of the
list.

1. **The reference mesh is not a target.** The user's correction, recorded in
   full in [P26](#p26): it is *"somewhat hand tuned, but not extensively, since
   tuning revealed huge parameter dependence."* Every "bring the five into line
   with the reference" reading in this file is therefore wrong — that would buy
   comparability across six untrustworthy reconstructions, not correctness. It
   also means the parameter dependence [P24](#p24)(a) measured was already known
   at the bench, which is independent support for [I03](#i03).
2. **[P09](#p09)(a) is a real error on five of six images, and fixing it needs no
   mesh rebuild.** Interior (face-less) vertices are 73–85% of the batched
   meshes against 18.8% on the reference, so 39–75% of adhesions measure
   `dist_mesh_um` to a point *inside* the wall (p90 understatement up to
   5.5 µm). Building the KD-tree from surface vertices fixes it; only a
   re-measure is needed, not a re-batch. It is a small signature change to
   `compute_measurements` plus both call sites, not an edit in place — see the
   entry.

**Revised queue: [P09](#p09)(a) → [P03](#p03) → the mesh rebuild tool
([I04](#i04)-shaped) → [P21](#p21)/[I03](#i03) (+ [P24](#p24)) → [P02](#p02) →
[P08](#p08)/[P10](#p10).** [P26](#p26) is no longer a work item in its own right
— with the reference demoted there is nothing to standardise *onto*, so it
collapses into [I03](#i03): pick parameters by score, then rebuild all six. P09(a)
goes first purely on cost-to-value; it is the only open correctness defect whose
fix requires neither a decision nor a recompute of anything expensive.

**Update, 2026-08-12 (later) — [P09](#p09)(a) is fixed and all six images are
exported.** `batch_measure --redo-json` over the folder: **6 measured, 0 failed**,
the only degradation being the expected FRET-not-configured warning on each.
`dist_mesh_um` medians rise 1.4–2.3× on the five batched images; the reference's
re-export leaves 21 of 24 columns byte-identical. Full numbers in
[P09](#p09)'s History.

> ### The current authoritative CSVs are the six `20260812_011*`
>
> One per image, in each `<stem>_saved/`. **These are the first exports in the
> project's history to carry every landed fix at once** — [P01](#p01)
> (curvature winding + orientation), [P07](#p07) (border seeds), [P04](#p04)
> (split nuclei, `nuc_is_fibroblast`), [P25](#p25) (count clearing) and
> [P09](#p09)(a) (surface-only distance) — and the first for five of the six
> images at all. Every earlier CSV in the reference folder is superseded; see
> the box in [P01](#p01) for what each one was.
>
> They are **not** trustworthy for the geometry question. Nothing here touches
> [P24](#p24), [P08](#p08)(e) or [P26](#p26): `curv_*` still has no
> scale-independent value, `path_curv_xy` still fits a 20 µm window to a ~140 µm
> vessel, and the six meshes are still mutually incomparable. What improved is
> `dist_mesh_um`, and the fact that all six images now have an export at all.

**The queue is now [P03](#p03) → the mesh rebuild tool → [P21](#p21)/[I03](#i03)
(+ [P24](#p24)) → [P02](#p02) → [P08](#p08)/[P10](#p10).** [P03](#p03) is
unblocked and needs no new tooling — six CSVs and six meshes now exist, which is
everything its "replicate across the other five images" step was waiting for.

**Update, 2026-08-12 (later still) — [P03](#p03) replicated on all six; it
holds.** Far-wall adhesion detection density is **19–77% of near-wall, in the
same direction in every image** (6/6, both split definitions), and median
adhesion volume does *not* fall consistently with depth (−19% to +16%, 3 up 3
down) — so this is objects going undetected, not objects getting smaller. Full
table in [P03](#p03)'s History. Two consequences:

- **[P06](#p06)(a) is discharged, with the opposite answer.** Its fix list starts
  "confirm it replicates across all six images". Done — and **the tension-vs-depth
  drift does *not* replicate**: ρ runs −0.241 to +0.307 and is negative in only
  2 of 6, while the reference reproduces its recorded −0.187 exactly. A sign that
  flips rules out optical attenuation, which must act the same way in every
  image, and supports the geometric reading [P06](#p06) already contains (the
  vessel dives through the stack, so depth is partly axial position, whose
  orientation is arbitrary per field). **[P03](#p03) and [P06](#p06) were bundled
  into one "depth" problem and they have now come apart:** the detection
  truncation is real and unanimous, the ratio drift is per-image and
  inconsistent.
- **The magnitude is not transferable.** 19%–77% is a 4× spread, so "the far
  wall yields 37%" cannot be used as a correction coefficient — any depth
  correction has to be fitted per image.

**Queue: the mesh rebuild tool → [P21](#p21)/[I03](#i03) (+ [P24](#p24)) →
[P02](#p02) → [P08](#p08)/[P10](#p10).** Both S1s at the front of the list are
now measured rather than open questions: [P03](#p03) is confirmed and
[P09](#p09)(a) is fixed. What remains ahead of the geometry work is nothing —
the reconstruction is next.

**Update, 2026-08-12 (later still) — the mesh rebuild tool and the scorer both
exist; [I03](#i03) is now down to its last third.** `pipeline/batch_mesh.py`
rebuilds mesh / curvature / centerline headlessly ([I04](#i04)), and
`pipeline/mesh_score.py` scores the result on four axes — `scale_stability`,
`topology`, `coverage`, `support` — validated against a sphere and a torus
rather than only against this file's own recorded numbers ([I03](#i03) History).
Three consequences for entries that were previously arguments rather than
measurements:

- **[P24](#p24)'s plateau test is now a calibrated number.** `slope_loglog` = 0
  means κ2 is a property of the surface, −1 that it is a readout of the fit
  window. A sphere returns +0.02; the reference vessel mesh returns **−1.19**,
  and *no* reconstruction tried so far reaches a plateau.
- **[P21](#p21)'s alpha knob is now mechanistically understood.** Max fabricated
  spur ≈ 0.95 × `alpha3d_um`, max edge ≈ 2 × `alpha3d_um`, independent of
  `subsample_um`. `alpha3d_um` **is** "the longest spur out of the wall I
  permit".
- **[P19](#p19) is weaker than it looked and [P20](#p20) is not a knife-edge.**
  The `min_size` floor sits mid-plateau over an 8× range, and `subsample_um`
  1 → 5 moves `slope_loglog` by only 0.185.

**[I05](#i05) is opened by this work, not closed by it.** With comparison finally
possible, the alpha shape was compared — against the other two in-repo methods
(neither completes in an hour) and against a prototype density-field surface
(always manifold, spur-free, but it cannot yet be smooth and covering at the same
time). Verdict: alpha 3D stays the default on the evidence, and is not settled.

**Queue: [I03](#i03)'s selection rule → the six-image rebuild it implies (+
[P24](#p24), [P21](#p21), [I05](#i05)) → [P02](#p02) → [P08](#p08)/[P10](#p10).**
What is left in [I03](#i03) is **choosing**: nothing weights the four scores, and
deliberately so — on the reference/defaults pair they rank the two meshes
*oppositely*, so any weighting would be the thing deciding the answer. Until a
rule exists there is no way to say which of two candidate surfaces is better,
which is the trap this whole line of work was opened to escape.

**Update, 2026-08-13 — centerline work landed off-queue.** Three changes to
`pipeline/centerline.py` and the vessel-mesh tab, none of them from the queue
above, all recorded in their entries today: the skeleton is now padded against
the field-of-view border ([P08](#p08)(c)), the overlay's Z comes from a 3D
distance-from-wall field rather than a column mean ([P22](#p22)), and the same
computation yields an **inscribed lumen radius** — the pipeline's only direct
caliber measurement, which independently confirms [P24](#p24)(b). One new defect
came out of it: [P27](#p27), the centerline is a point *set*, not a traced path.

**Update, 2026-08-13 (later) — the direction changes, and the reason is
[P28](#p28).** A session spent asking what the centerline is actually *for*
ended somewhere else. Three findings, in the order they landed:

1. **[P29](#p29) — the centerline forks at a field-of-view corner**, because a
   notch in the mask that opens onto the image border cannot be closed by
   `binary_fill_holes` and the medial axis runs around both sides of it. The
   border-padding fix does **not** address this, and the false branch node is in
   the saved file, so `dist_branch_um` is wrong for **18.4%** of rows.
2. **The cross-section is not circular and not clipped.** Measured directly
   from the enclosed lumen on the left arm: 151–180 µm wide in plane, 83–92 µm
   deep, **aspect 1.64–2.15**, with only 0.01% of the mask in the first z slice
   and none in the last. Recorded in [P24](#p24)(b), which had been reading two
   measurements of the *minor* axis and one of the *major* axis as three
   estimates of one radius.
3. **[P28](#p28) — the cadherin threshold is throwing away real wall.** The
   gaps the reconstruction has to bridge mostly are not holes: ray-casting the
   raw channel finds signal at **the same radius** as the detected wall, ~40% as
   bright, at 5–6× the lumen noise. The raw ring is continuous; the mask is a
   dashed copy of it.

**What that does to the plan.** The reconstruction work had been heading toward
a centerline-frame sweep, which needs the centerline to be ordered
([P27](#p27)), accurate, branch-aware and clean at the edges — i.e. it needs
every open centerline defect fixed first. It is superseded by [I06](#i06),
which needs the centerline only to say **which side of the wall the lumen is
on**, a job robust to every one of those defects. The corner fork, the notch and
the ordering all stop being blockers to the surface. They remain real defects;
they are no longer on the critical path.

**Scope, decided 2026-08-13:** the deliverable is the **geometry pipeline**, built
in anticipation of a FRET calibration that does not exist yet. The force columns
stay blank and are not a blocker — which is workable precisely because
[I06](#i06)'s acceptance test ([P24](#p24)'s plateau) needs no force at all. One
consequence: anything that would be *silently wrong the moment force is
populated* should be fixed while it is still blank. That is [P10](#p10), and
only [P10](#p10).

**Queue: [I06](#i06) (+ [P28](#p28) as its input) → [P29](#p29) → re-export the
six → [P02](#p02)/[P03](#p03) when the force data exists.** [I03](#i03)'s
selection rule is still unwritten, but [I06](#i06) can proceed ahead of it: a
weighting is only needed to break a trade, and if the new surface dominates
alpha 3D on all four scores no weighting is required to prefer it.

**Update, 2026-08-13 (later still) — [I06](#i06) is built and the plateau is
reached; it does *not* dominate, so the selection rule is now load-bearing.**
`pipeline/vessel_poisson.py` returns `slope_loglog` = **+0.009** on the reference
image against alpha 3D's **−1.428** on the same mask, with 0 non-manifold faces
against 18.7% and χ = −4 against +2,870. **This is the first reconstruction in
the project's history to plateau**, and [P24](#p24)'s pass/fail is met.

But the hoped-for clean sweep did not happen: **alpha 3D still wins both
coverage columns** (77.0% vs 61.6% at 2 µm, 99.9% vs 87.9% at 5 µm) **and
support max reach** (9.87 µm vs 26.70 µm). Full table in [I06](#i06)'s History.
So the paragraph above has been tested and its escape clause did not apply — the
four scores rank the two methods oppositely again, exactly as they did on the
reference/defaults pair, and **[I03](#i03)'s missing selection rule is now what
stands between "the plateau is reachable" and "this is the reconstruction we
use".**

That said, the trade is not symmetric and the next step is not the rule. The
coverage gap is a statement about *how many mask voxels lie near the surface*,
and [P28](#p28) has already measured that a large share of what the mask calls
"gap" has real signal in it while a share of what it calls "wall" is noise. Which
way the coverage comparison should come out therefore depends on a question
[P28](#p28) answers and this comparison assumes. **Run [P28](#p28) first, with
the method held fixed** — that is the controlled experiment [I06](#i06) was
sequenced around, and it is cheaper than writing a weighting nobody can justify.

Two things also came out of the build and are recorded in their entries:
[P30](#p30) (`build_mesh` silently treats an unknown method string as alpha 3D)
and, in [P01](#p01)'s History, a sign check that compared a median against a p99
and thereby inverted a correct curvature field on the real image — caught only
because the plateau came out at −0.447 instead of +0.009.

**Update, 2026-08-21 — [P28](#p28) ran as the controlled experiment, and it
came back negative on the fix while confirming the defect.** Three results, in
the order they matter:

1. **The defect replicates.** Five of six images (the sixth sealed no
   cross-section and is unmeasurable, not contrary): gap-direction peaks
   42–50 counts above lumen background against 9–24 for no-wall controls, sitting
   **0.0–1.0 µm from the wall radius the mask itself reports**. The dim wall is
   real, and the "same radius, so it is the same wall" claim is now a much
   tighter number than the original measurement gave.

2. **But [P28](#p28)'s own proposed fix aims at the wrong gate.** The mask is
   built by three gates and the attribution is unambiguous: the Sauvola **local
   threshold passes 100% of gap-wall peaks in every image**. What deletes them is
   the **global intensity floor** (15–41%) and the **`min_size` component
   filter** (59–85%). So "stop binarising this early" targets the one gate that
   is innocent, and **[P19](#p19) is now the larger half of [P28](#p28)** — it is
   raised to S2 accordingly.

3. **The continuous field was built, measured, and does not work.** The whole
   plumbing landed (`denoise.local_contrast_score`, the field carried through
   `build_mesh` and `occupancy_grid`, `batch_mesh --cadherin-field`), with alpha
   3D verified **byte-identical** and the recorded Poisson JSON reproduced to
   every digit, so the comparison is genuinely controlled. Both methods then
   **degrade** on the field — Poisson's plateau is lost (+0.009 → −0.399) and
   coverage collapses (61.6% → 7.8%); alpha 3D is worse still (77.0% → 4.2%,
   support reach 26.7 → 101.7 µm). The cause is in [P28](#p28)'s own data: the
   local z-score does not separate dim wall from lumen noise, so no `k_lo`
   exists that keeps one and drops the other. Full table in [I06](#i06).

**What that does to the plan.** The escape clause hoped for above is now closed
in both directions: [P28](#p28) has been run with the method held fixed, and the
coverage comparison it was supposed to reinterpret **stands unchanged**. Alpha 3D
still wins coverage, Poisson still owns the plateau, and **[I03](#i03)'s
selection rule is now unavoidable** — there is no longer a cheaper experiment
sitting in front of it. [P28](#p28) also handed [I03](#i03) a concrete warning
about how *not* to write that rule: a mesh with 4.2% coverage and a foam of 85
pockets scored `slope_loglog` = **+0.024**, better than the accepted
reconstruction's +0.009, because a single log-log fit through a non-monotonic
curve averages to zero.

One incidental positive: the cadherin masks are **exactly reproducible** from
the recorded recipe (Jaccard 1.000, all six images), so the [P26](#p26)-shaped
worry about unrecorded threshold parameters does not stand and no entry was
opened for it.

**Queue: [I03](#i03)'s selection rule → [P19](#p19) (now the larger half of
[P28](#p28)) → [P29](#p29) → re-export the six → [P02](#p02)/[P03](#p03) when the
force data exists.** [P28](#p28) is no longer in front of [I03](#i03): it has
been run, and it did not move the comparison the rule has to resolve.

**Update, 2026-08-22 — the 2D centerline is deleted and Tier 2 of
`docs/DIRECTION.md` is built.** `pipeline/centerline.py` is gone, replaced by
`pipeline/lumen_graph.py`, a 3D graph on the Poisson lumen. Four entries change
status and one is new:

- **[P08](#p08), [P23](#p23), [P29](#p29) — superseded.** The code they describe
  no longer exists. They are left open-as-record rather than withdrawn: what
  each measured was true of the code at the time.
- **[P27](#p27) — fixed.** The replacement is ordered along each edge and
  carries arc length by construction, so everything that entry lists as blocked
  is now available.
- **[P32](#p32) — new.** The fitted lumen leaks through the z faces on images 4
  and 5, so the graph invents junctions there and clips half its cross-sections.
  Diagnosed and warned about in code; not fixed.

Two results that change what the queue is *for*, both in `docs/DIRECTION.md`
§8C:

- **In-plane caliber is measurable at last** (§8C.3) — area, both extents and
  aspect, from cuts perpendicular to the axis, agreeing with
  [P24](#p24)(b)'s direct measurement by a route sharing no machinery. This is
  what the shear estimate needed and it was the thing no projection could
  supply.
- **The signed nucleus-to-surface distance agrees with [P05](#p05)'s classifier
  on 5 of 6 images** (§8C.4), while the unsigned version is chance. But the
  attempt to validate it *without* a marker failed structurally (§8C.5), so
  P05's blocker is now "find an independent endothelial label", which is
  [I02](#i02)/[P14](#p14) territory.

**Revised queue: [I03](#i03)'s selection rule (now also gating whether
alpha→poisson becomes the shipped reconstruction, and therefore whether any
geometry column exports at all) → [P32](#p32) → [P19](#p19) → re-export the six
→ [P02](#p02)/[P03](#p03) when the force data exists.** [P29](#p29) leaves the
queue by being superseded rather than fixed.

> **One state to be aware of.** The canonical `vessel_mesh.npz` is still alpha
> 3D and there is no untagged `vessel_graph.npz`, because the graph needs a
> Poisson lumen and adopting alpha→poisson as the default is [I03](#i03)'s
> decision. So **every geometry column currently exports blank** rather than
> mixing two reconstructions in one CSV. That is deliberate; `--check-only`
> says so per folder.

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
| **Resolved** | 2026-08-10 (commit `af8e6c0`) |
| **Where** | `pipeline/cell_mesh.py:137` (`alpha_shape_faces`), `pipeline/cell_mesh.py:305` (`vertex_outward_normals`) |
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

> **Correction, 2026-08-11 (left in place above so the original reading is on the
> record).** The last sentence overstates it, and the parenthetical "radius
> ≈ 13.5 µm" is not a measurement of the vessel — it is κ2 inverted. The vessel's
> real radius is of order **40 µm** ([P24](#p24)(b)), so κ2 = −0.074 is 3–5× too
> curved to be its circumferential curvature, and the "textbook tube" agreement is
> self-confirming. The *strongest* evidence for the fix is the winding and
> orientation evidence — coherence 0.209 → 0.972, directed-edge twins 53.7% →
> 100%, the synthetic-sphere check — none of which depends on caliber. See
> [P24](#p24) for what κ2 is actually tracking.

**Only the default mesh method is affected.** `mesh_from_volume` (marching cubes,
used by the polar and alpha-2D methods) returns consistently-wound faces. The
GUI default is `alpha 3D (points)` (`pipeline/gui/vessel_mesh.py:80`), which is
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

> **Caveat, 2026-08-11.** The `κ2 ≈ −1/R` clause of this test is **circular** and
> should not be read as evidence: the R it is scored against was itself read off
> κ2. The vessel's true caliber is ~5× larger. See the 2026-08-11 History entry
> and [P24](#p24). The other three clauses stand, and the fix itself stands.

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

- 2026-08-11 — **verified end to end on a regenerated CSV.**
  `vessel_mesh.npz` (new winding) and `vessel_curvature.npz` were rewritten for
  `V3_..._04step_002` and `batch_measure` re-exported it as
  **`..._adhesion_measures_20260811_105032.csv`** (2,147 rows) — the authoritative
  post-fix export for this image. The acceptance test passes on the actual
  exported columns, not just in memory:

  | column-level check | old CSV (20260806) | **new CSV (20260811)** | target |
  |---|---|---|---|
  | `curv_k1` median | +0.0167 | **−0.0034** | ≈ 0 |
  | `curv_k2` median (implied R) | −0.0422 (23.7 µm) | **−0.0701 (14.3 µm)** | ≈ −1/R |
  | `curv_k2` sign consistency | 90.1% | **99.5%** | ≥ 96.9% |
  | `curv_gaussian` median | +0.00005 | **+0.00025** | ≈ 0 |
  | ρ(`curv_mean`, `centroid_z_um`) | −0.059 | **+0.040** | → 0 |

  Sign consistency reads higher per-adhesion (90.1 / 99.5%) than per-vertex
  (78.6 / 98.3%) because adhesions do not sample surface vertices uniformly —
  they concentrate where the wall is well reconstructed. Quote the per-vertex
  figure for the mesh and the per-adhesion figure for the CSV; they are different
  measurements of the same fix.

  **Three CSVs now sit in that folder — use `20260811_105032`.** `20260806_145720`
  is the pre-fix baseline, kept deliberately so the before/after above can be
  re-derived. `20260811_104417` is a **superseded intermediate**: it has the
  corrected curvature but was exported before `fibroblast_count.npz` was
  recomputed, so it still carries the pre-[P07](#p07) `n_nuclei_total = 273`.
  Delete it or ignore it; do not analyse it.

  > **Superseded, 2026-08-12 — there are now four, and the one to use is
  > `20260812_011205`.** The paragraph above is left as written because the
  > `20260811_105032` re-export is what the curvature evidence in this entry was
  > measured on, and that provenance should not be rewritten. But it predates
  > the [P04](#p04) nucleus fixes and [P09](#p09)(a), so it is no longer the
  > file to analyse. Current standing of the four:
  >
  > | CSV | status |
  > |---|---|
  > | `20260806_145720` | pre-[P01](#p01) baseline — keep, do not analyse |
  > | `20260811_104417` | superseded intermediate — ignore |
  > | `20260811_105032` | the curvature evidence in this entry; superseded for analysis |
  > | **`20260812_011205`** | **authoritative** — every landed fix |
  >
  > The four `curv_*` columns are **byte-identical** between `20260811_105032`
  > and `20260812_011205`, so nothing in this entry's verification is affected
  > by the move.

- 2026-08-11 — **the acceptance test's radius clause is circular; the fix is
  not affected.** Found while investigating whether the X-Y projection in
  `path_curv_xy` costs anything ([P08](#p08)(e)).

  The entry above reads κ2 = −0.074 → "circumferential, radius ≈ 13.5 µm", then
  adopts 13.5 µm as the target and scores the fix's 14.2 µm against it. **The
  target was derived from κ2 itself**, so the agreement is self-confirming and
  the "textbook tube" reading rests on nothing external. Two independent
  measurements of the vessel's actual caliber, neither using the mesh:

  | measurement | result |
  |---|---|
  | z extent of wall in a 5 µm slab ⟂ to the local axis | **~86 µm**, saturating (87.1 → 85.8 as the in-plane search radius goes 30 → 110 µm) |
  | xy inscribed half-width of the footprint at centerline points | median **70.6 µm** (p10 21.4, p90 87.8) — an upper bound, `binary_fill_holes` inflates it at the Y |

  So the vessel radius is of order **40 µm, not 13.5 µm**, and κ2 = −0.0703 is
  ~3–5× too curved to be this tube's circumferential curvature. What it is
  measuring instead is in [P24](#p24) — the implied radius tracks `radius_um`
  with no plateau, so it is a property of the reconstruction scale.

  **What still stands, and why the fix is not in question.** Every measurement
  that established the fix is an *internal-consistency* check and none of them
  depends on caliber: face-set identity before/after part A, directed-edge twin
  fraction 53.7% → 100.0%, face-normal coherence 0.209 → 0.972, κ2 sign
  consistency 78.6% → 98.3%, ρ(shipped, corrected) = +0.069 at 54.8% sign
  agreement, and the synthetic-sphere check on part B. The columns really were
  noise before and really are self-consistent now. What is withdrawn is only the
  claim that they are *validated against the vessel* — they are not, and
  [P24](#p24) is why.

- 2026-08-13 — **part B now has a principled replacement, for one method.**
  `vessel_poisson` reads the winding of every face from the sign of
  `face_normal · ∇χ` at the face's own centroid, where χ is the fitted indicator
  function. That needs neither a closed surface (the Poisson surface is open at
  the x/y domain faces by construction, so the enclosed-volume branch is
  unavailable to it) nor the star-shaped-about-the-centroid assumption this
  entry measured wrong on **44.4%** of the reference mesh's vertices. On the
  synthetic elliptic cylinder **100.0%** of the resulting vertex normals agree
  in sign with the analytic outward normal, and κ2 sign consistency is 100.0%
  (81.6–94.5% on the real image, where the surface is not an ideal tube).

  This is a per-method fix, not a fix to `vertex_outward_normals`: the
  centroid fallback is still what the three alpha/polar methods get, and it is
  still wrong on 44.4% of vertices there. What changed is that
  `cell_mesh.mesh_curvature` now takes a `normals=` argument, so a method that
  knows its own orientation can supply it instead of having one inferred.

  **A near miss worth recording.** The first version of `vessel_poisson`'s
  global-sign guard compared the *median* χ over the lumen core against the
  *99th percentile* over the exterior. The exterior is ~9× the core's volume and
  includes every cell just outside the wall, so its p99 sits inside the wall's
  own transition: measured on the reference image, core median +0.0762,
  exterior median −0.0037, exterior p99 +0.0985. The guard therefore negated a
  perfectly correct field and inverted every curvature sign downstream — the
  exact failure it existed to catch, and invisible in every score
  (`slope_loglog` read −0.447 rather than +0.009, which looks like a mediocre
  result rather than a sign error). It now compares median to median. The
  lesson generalises: **a sign check must compare like with like.**


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
| **Status** | open — **confirmed on all six images 2026-08-12** |
| **Severity** | S1 |
| **Opened** | 2026-08-10 (audit F3) |
| **Last verified** | 2026-08-12 (all six; far wall 19–77% of near, unanimous) |
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

**History**
- 2026-08-10 — opened from `audit.md` F3, measured on one image.
- 2026-08-12 — **replicated on all six. It holds, in every image.** Same method
  as the entry above: adhesion count per depth band ÷ vessel surface area in
  that band, triangle areas from each saved mesh, adhesions from the
  `20260812_011*` exports. The reference image's five band areas reproduce this
  entry's table **exactly** (134,908 / 46,667 / 28,928 / 36,480 / 128,890 µm²),
  which is what makes the other five comparable to it; its counts read slightly
  higher (2,147 against 2,092) because [P07](#p07) landed in between.

  Far-wall detection density as a percentage of near-wall, split two ways —
  at an absolute 40 µm as the original did, and at each stack's own midpoint,
  because the six stacks are 80–117 µm deep and a fixed cut is not the same
  fraction of each:

  | image | stack | far % of near (40 µm cut) | far % of near (midpoint) |
  |---|---|---|---|
  | `V3_..._04step_002` (reference) | 105 µm | **36%** | 43% |
  | `V3_..._04step_003` | 98 µm | 61% | 68% |
  | `V3_..._09step_001` | 115 µm | 62% | **77%** |
  | `V3b_..._04step_001` | 80 µm | **19%** | **19%** |
  | `V3b_..._04step_002` | 117 µm | 45% | 30% |
  | `V3b_..._04step_003` | 94 µm | 32% | 36% |

  **Unanimous in direction, 6 of 6, on both splits.** By this file's own
  ceiling ([Analysis notes](#analysis-notes) — "if all agree in direction a sign
  test floors at p = 0.031"), that is as strong as a six-image result gets. The
  original 37% is representative rather than extreme; the median across images
  is ~40% on both splits.

  **But the magnitude is not a constant — it spans 4× (19%–77%).** "The far wall
  yields 37%" is a fact about one image, not a coefficient to correct with. Any
  depth correction has to be fitted per image.

  **The truncation reading is strengthened, not just replicated.** This entry
  argued that flat adhesion volume vs depth is the wrong statistic because
  truncation removes objects rather than shrinking them. Across six images the
  median volume shallow → deep moves **−19% to +16%, and is not even consistent
  in sign** (3 down, 3 up) while detection density falls in all six. A mechanism
  that shrinks adhesions cannot produce that; one that fails to detect them can.
  The single image's −11% was within the noise of zero.

  **Two things that make this conservative rather than overstated.** The
  cadherin channel attenuates with depth as well, so the *denominator* — far-wall
  surface area recovered by the mesh — is itself likely underestimated. That
  inflates far-wall density and makes the truncation look milder than it is.
  And the split is by absolute z, while [P06](#p06) established that on this
  geometry the vessel dives through the stack, so "far wall" is partly *further
  along the vessel*; that confound is unresolved and cuts both ways.

  Densities are **not** compared between images here, only the within-image
  ratio — the six meshes are built at different parameters ([P26](#p26)) so the
  absolute per-1000 µm² figures are not commensurable. The alpha-shape mesh is
  also a closed solid whose face area counts both wall sheets, so the absolute
  area is roughly 2× a true surface; constant within an image, so the ratio
  survives.

  `09step_001` (0.9 µm z-spacing) shows the **least** truncation of the six on
  both splits. Not interpreted — one image at that spacing, and it is the
  direction a coarser z would push things if fewer, thicker planes recover more
  deep signal. Recorded because it should not be pooled without checking, the
  same caution [P05](#p05) reached from the other direction.

  **What this settles.** The "confirm it replicates across all six images" step
  that [P06](#p06)(a) is blocked on is **done**. Depth-stratified comparison on
  this dataset contrasts a near-complete shallow sample against a systematically
  incomplete deep one, in every image. Correcting intensities for attenuation
  cannot fix that, because the missing adhesions were never measured.

- 2026-08-24 — **the optical reading is confirmed against its main rival, and
  this entry's argument from flat adhesion volume is corrected.** Scripts:
  `docs/scratch/tier2/{orientation_vs_depth,circumferential}.py`. Scratch is
  gitignored, so the numbers here are the tracked record. Three parts.

  **(a) The flat volume distribution does not mean what this entry says it
  means.** The entry argues "truncation removes objects rather than shrinking
  them" and cites median volume falling only 11%. Checked against the whole
  distribution rather than the median — volume percentiles by depth quartile,
  reference image:

  | band µm | n | p10 | p25 | p50 | p75 | p90 |
  |---|---|---|---|---|---|---|
  | 1–10 | 537 | 4.77 | 6.97 | 13.04 | 20.60 | 29.31 |
  | 10–14 | 536 | 5.34 | 8.29 | 13.49 | 21.70 | 31.12 |
  | 14–74 | 537 | 5.42 | 8.38 | 13.90 | 22.11 | 31.79 |
  | 74–103 | 537 | 5.22 | 7.62 | 12.50 | 20.87 | 28.73 |

  Flat at **every** percentile, in all six images, while the count falls 63%.
  So the observation is stronger than the entry claims — but the *reasoning* is
  wrong. Under a threshold applied directly to a soft-edged object, dimming
  shrinks the object before it vanishes, so shrinkage is what should be seen.
  The reason it is not is **architectural and specific to this pipeline**: the
  global Otsu is not the segmenter, it is a **gate** (`gate_by_clover`, default
  on, dilated 1 µm) deciding which voxels the RF may look at, and the RF
  segments inside it. A gate acts all-or-nothing per object — clear the
  threshold and you are segmented at full extent, miss it and you are invisible.
  Count falls, surviving sizes unaffected. "Truncation removes objects rather
  than shrinking them" is true *here because of the gate*, not as a general
  principle, and the entry should not be quoted as if it were one.

  There is a hard floor at exactly 3.433 µm³ in every image, but p10 sits at
  ~5 and does not drift toward it with depth, so a size cutoff is not what
  produces the flatness.

  **(b) A rival explanation this entry never considered: gravity.** The
  fibroblasts were seeded onto an already-formed vessel
  (`docs/DIRECTION.md` §1). If they settled under gravity the up-facing wall
  gets more cells — real biology, large, and pointing the same way as the
  optical artifact. A falling count with an unchanged size distribution is what
  **both** produce, so nothing previously recorded here separated them.

  **The discriminator, and it comes out optical.** Attenuation depends on
  absolute depth in the stack; gravity depends on which way the wall faces. On
  a tube those are normally locked together, but the vessel axis dives 6–36 µm
  (`docs/scratch/tier2/depth_confound.py`; p90 |dz/ds| 0.33–0.71), so the
  up-facing wall exists across a range of
  absolute depths. Holding orientation fixed at the nearly-flat top wall
  (per-vertex ∇χ normals, `n_z < −0.9`), shallow half vs deep half:

  | image | z span µm | shallow | deep | ratio |
  |---|---|---|---|---|
  | `V3_..._04step_002` | 38 | 21.17 | 15.09 | **1.40** |
  | `V3_..._04step_003` | 18 | 15.37 | 20.85 | 0.74 |
  | `V3_..._09step_001` | 24 | 26.77 | 10.23 | **2.62** |
  | `V3b_..._04step_001` | 67 | 22.95 | 24.67 | 0.93 |
  | `V3b_..._04step_002` | 27 | 37.04 | 17.48 | **2.12** |
  | `V3b_..._04step_003` | 34 | 36.72 | 12.23 | **3.00** |

  **Four of six fall with depth at constant orientation** (up to 4.7× at
  `n_z < −0.95`). Gravity cannot produce that — it is one wall facing one way.
  The two exceptions are the two with least power: `04step_003` has 18 µm of
  span (the flattest vessel, its axis barely dives) and `V3b_..._04step_001` is
  [P32](#p32)'s z-leak image.

  **This does not exclude a gravity component**, and the tempting comparison
  does not work: the up/down gap (1.4–15.8×) exceeds the within-orientation
  gradient, but the down-facing wall is imaged through the whole lumen and
  everything above it, so the two populations do not share an attenuation law
  and the extrapolation is meaningless. Gravity is **not required**; it is not
  ruled out.

  A confound had to be fixed mid-run and is recorded so it is not reintroduced:
  at a loose orientation cut (`n_z < −0.5`) the "up-facing" set still includes
  surface wrapping toward the sides, and that wrap *itself* increases z, so
  binning by depth inside it silently re-measures circumferential position —
  the exact confound the test exists to remove. Stable from 0.8 up.

  **(c) A null worth recording so it is not re-proposed.** This entry's five
  depth bands are not monotonic (10.39 / 2.10 / 0.73 / 2.91 / 3.61), and the
  band areas show why: the extreme bands are the top and bottom walls, the
  middle band is the sides. That invited the reading that part of it is real
  circumferential shear — for an elliptical section τ(minor-axis end) /
  τ(major-axis end) = a/b, so ~2× at the measured aspects, high at **both** z
  extremes where attenuation can only ever run one way.

  **Measured, and it is absent.** Binning adhesion density by angle around the
  perpendicular cross-section, minor-axis ends vs major-axis ends: **1.04 /
  1.25 / 0.95 / 0.87** on the four images with reliable geometry (2.07 and 1.09
  on [P32](#p32)'s two). Predicted ~2, observed ~1. The apparent structure in
  the z bands is depth, not angle.

  Caveat on that null, and on a discarded first version of it: excluding
  `censored` centerline points removed 50–76% of the data non-randomly by angle
  and manufactured a spurious anti-correlation (0.55/1.02/0.81/0.68). The
  radial gate is inert by comparison (2.5× / 5× / none agree to ±0.05). The
  null also has little power at the deep end, where attenuation suppresses
  detection ~2–3× — the same size as the predicted effect — but the shallow end
  is *brighter* than the sides and should have been biased toward the predicted
  result, and still reads 1.0.

  **What this settles.** The root cause named in **Where** is confirmed as
  optical rather than biological, for the first time with evidence that
  separates the two. That promotes the fix: see [I07](#i07).

- 2026-08-25 — **the clover gate is ruled out as the cause.** [I07](#i07)
  proposed local thresholding on the gate as this entry's fix. The gate was
  widened on image 6 (the worst image here, shallow/deep 3.00) from Otsu's
  0.969% of the volume to 93.4%, and the `proba >= 0.5` set came back
  **bit-identical, 0 voxels different** ([P43](#p43)). The gate was never
  binding, so no change to it can fix this. The truncation is downstream — in
  the RF itself or in the post-process `min_size_vox=100`.
  **New live hypothesis: annotation sampling.** `rf.py::train` fits on
  hand-painted voxels over the *raw* volume (`np.argwhere(ann_vol > 0)`), so the
  gate never touched training — but if the painted examples were drawn mostly on
  shallow, bright adhesions, the model learned that brightness profile and will
  not fire on deep dim ones. That would make the fix *retraining with deep
  examples*, not any threshold change. Untested; the discriminating measurement
  is the depth distribution of `proba >= 0.5` voxels against that of the surviving
  instances, which separates an RF cause from a `min_size` cause.

---

<a id="p04"></a>
## P04 — Two nucleus label maps; `nuc_id` is not joinable to cell type

| | |
|---|---|
| **Status** | fixed |
| **Severity** | S1 |
| **Opened** | 2026-08-10 (audit F4) |
| **Last verified** | 2026-08-11 |
| **Resolved** | 2026-08-11 (commit `04262b8`) |
| **Where** | `pipeline/gui/nuclei.py:218` vs `pipeline/gui/measure.py:62` |
| **Affects** | `nuc_id`, `dist_nuc_um`, `n_fibroblasts`, `n_nuclei_total`, and the new `nuc_is_fibroblast` |
| **Evidence** | `docs/audit.md` §3 F4 |

The Nuclei tab builds two separate segmentations and only one reaches each
consumer:

- **Otsu + connected components** → `state["dapi_state"]["labels"]`. Not split,
  not size-filtered beyond ≥10 voxels.
- **Watershed-split nuclei** (`nuc`, `gui/nuclei.py:218`), built from the
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

**The fix (as originally stated — incomplete, see below):** write the split
labels back to `dapi_state`, or have Measure use the split map. Either way, one
nucleus definition, size-filtered once.

> **Correction, 2026-08-11 (the original left in place above).** Neither of those
> delivers the join this entry was opened for. `is_fib` — the per-nucleus
> fibroblast classification — is computed at `gui/nuclei.py:243-246` and then
> **thrown away** at `:249-271`: only four scalars (`total`, `fibroblast`,
> `other`, `threshold`) survive into `state`, and the class map is rendered
> straight into a napari layer. Unifying the label maps would make `nuc_id` point
> at nuclei, which is necessary, but there would still be nothing to join it
> *to*. The split labels, `ids` **and** `is_fib` all have to be co-persisted.

**The fix as implemented, 2026-08-11.** A new per-image artifact,
`nuclei_split.npz`, holding `labels`, `sizes`, `ids`, `is_fib`, `min_size_vox`
and `threshold` — and a new CSV column `nuc_is_fibroblast` (1.0 VnTs+, 0.0 other,
blank unknown).

**Why a second artifact rather than replacing `dapi_state["labels"]`.** These are
two genuinely different objects, and naming them separately is the fix rather
than a compromise:

| | `dapi_labels.npz` — **DAPI components** | `nuclei_split.npz` — **nuclei** |
|---|---|---|
| what it is | Otsu connected components, ≥10 vox | watershed-split, size-filtered |
| count here | 362 | 276 |
| consumers | `min_size` slider, size histogram, displayed mask, the watershed's own input, the RF DAPI gate | `nuc_id`, `dist_nuc_um`, `nuc_is_fibroblast` |

Replacing in place would have made Count non-re-entrant (`_count` reads
`dapi_state` to build its own watershed input), silently changed what the
`min_size` slider and histogram refer to, and — decisively — left a failed count
indistinguishable on disk from a successful one. No step reordering was needed:
`run_fibroblasts` already runs immediately after `run_dapi`.

`labels` and `is_fib` are written to **one file** deliberately: that is what makes
them impossible to stale relative to each other, which is the exact failure class
[P12](#p12) is about.

**Verified not to touch adhesion segmentation.** The concern with changing the
nucleus map is that it feeds the DAPI gate. It does not, headlessly: there is no
`run_cell_mesh` step, and `gate_by_dapi` / `gate_by_mesh` both default `False`
(`gui/adhesion_rf.py:363,368`), so the DAPI cell mesh is never built in the batch
and no nucleus-map change can reach the adhesion instances.

> ### `nuc_id` is not comparable across this fix
>
> Old CSVs' `nuc_id` values run to **355** and index the **362**-component
> unsplit map. New ones run to ~**276** and index nuclei. **The same integer
> means a different object.** Do not join, pool or diff `nuc_id` across the
> boundary; the `nucleus_map` key in each export's settings JSON records which
> map that CSV used. `dist_nuc_um` shifts too, and should — it is no longer
> measuring the distance to the nearest speck of sub-nuclear debris.

**History**
- 2026-08-10 — opened from `audit.md` F4.
- 2026-08-11 — **fixed.** Measured on `V3_..._04step_002`, comparing a
  re-measure of the same `_saved` folder with and without `nuclei_split.npz`
  (2,147 adhesions both times, every other column byte-identical):

  | | unsplit map (before) | split map (after) |
  |---|---|---|
  | `nuc_id` max | 355 | **275** |
  | distinct `nuc_id` | 176 | **170** |
  | `dist_nuc_um` median | 4.24 µm | **4.88 µm** (mean 5.26 → 6.20) |
  | rows where `dist_nuc_um` moved | — | 310 / 2,147 |
  | `nuc_is_fibroblast` | absent | **2,147 / 2,147 populated** |

  The join reconciles against the count: 424 adhesions sit on a VnTs+ nucleus
  and 1,723 on an "other" nucleus; those 424 belong to **25 distinct fibroblast
  nuclei out of the 27** in `fibroblast_count.npz` (the other two carry no
  detected adhesion), and 170 of 276 nuclei carry at least one adhesion.

  `dist_nuc_um` rising is the expected direction and is the clearest single
  confirmation that the entry's "half of what `dist_nuc_um` measures to is
  sub-nuclear debris" reading was right — with the debris gone, the nearest
  *nucleus* is further away.

  **Backward compatibility measured, not assumed.** A `_saved` folder with no
  `nuclei_split.npz` falls back to the unsplit map, and that path re-exported
  `V3_..._04step_002` byte-identically across all 24 pre-existing columns, added
  `nuc_is_fibroblast` blank, and said so with an explicit `WARN`. The fallback is
  loud precisely because it is the old, wrong-unit behaviour.

  Not addressed here, and unchanged: the classifier producing `is_fib` is still
  unvalidated ([P05](#p05)). **This makes the classification joinable, not
  correct** — `nuc_is_fibroblast` inherits every problem in P05, including the
  ~2× swing across one slider step.

- 2026-08-25 — **qualified, not reopened.** `nuc_id` is joinable, which is what
  this entry claimed and fixed. It is not *correct*: it is a nearest-neighbour
  join, and at 19–25 µm nucleus spacing against a 50–100 µm cell footprint it
  misassigns a spread cell's peripheral adhesions. Opened as [P42](#p42); no
  analysis should aggregate adhesions by nucleus until that is resolved.

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

**Related:** [P04](#p04) (the count is now joinable to `nuc_id` — but see the
History note below), [P07](#p07) (border nuclei lose their seed), [P13](#p13)
(the cadherin exclusion multiplier is hard-coded), [P16](#p16) (a claim from this
entry that was retired).

**History**
- 2026-08-06 — recorded in `TODO.md`.
- 2026-08-10 — the "over-splitting inflates `n_nuclei_total`" bullet split out
  and closed as [P16](#p16) (invalid); the shell-shape and border-seed mechanisms
  added from the audit.
- 2026-08-11 — **first cross-image measurement.** Until the five-image batch ran,
  the classifier had only ever been applied to one image, so every statement in
  this entry rested on within-image slider behaviour. All six now:

  | image | VnTs+ / nuclei | fraction |
  |---|---|---|
  | `V3_..._04step_002` | 27 / 276 | 9.8% |
  | `V3_..._04step_003` | 13 / 259 | **5.0%** |
  | `V3_..._09step_001` | 47 / 315 | **14.9%** |
  | `V3b_..._04step_001` | 18 / 230 | 7.8% |
  | `V3b_..._04step_002` | 29 / 290 | 10.0% |
  | `V3b_..._04step_003` | 17 / 174 | 9.8% |

  **A 3.0× spread across six images of one condition** (median 9.8%). This is
  *not* by itself evidence that the classifier is wrong — the real fibroblast
  fraction may genuinely vary between fields, and the threshold is
  background-relative so it adapts per image. It is recorded because it sets the
  scale of the problem: the between-image spread is the same order as the ~2×
  within-image swing across one slider step, so **the two cannot currently be
  told apart**. Any per-image analysis that conditions on cell type inherits
  that.

  The highest fraction is `09step_001`, the 0.9 µm z-spacing file. With one
  image at that spacing it is not separable from field-to-field variation, but
  it is the direction a coarser z would push a perinuclear-shell measurement,
  and it should not be pooled without checking.

- 2026-08-22 — **the first independent corroboration this entry has ever had,
  and it replicates on five of six images.** A signed nucleus-to-vessel-surface
  distance (`pipeline/nuc_side.py`; exact nearest-triangle, on the
  alpha->poisson mesh whose ∇χ normals make positive mean abluminal) agrees with
  this classifier at AUC **0.644–0.773** on five images, p from 6.9e-07 to 0.04,
  with image 3 (`09step`, the 0.9 µm z file) the lone null at 0.506. Full table
  in `docs/DIRECTION.md` §8C.4.

  **This is not validation and must not be quoted as such.** Two unvalidated
  measurements agreeing is exactly that — but they share no machinery and no
  failure mode: one reads mClover in a perinuclear shell, the other reads
  geometry against a Poisson surface. Fibroblasts were seeded onto an existing
  vessel and are abluminal, so the geometric measure has a physical reason to
  work that the marker measure does not supply.

  **The unsigned version is chance** (0.438–0.619, below 0.5 on two images),
  which confirms on six images what `docs/DIRECTION.md` §7.2 argued from one:
  the earlier "the two cell types do not separate by distance" result was the
  right observation read against the wrong test.

  **The attempt to validate it without a marker failed, structurally.** §7.2
  proposed calibrating on nuclei touching the VE-cadherin junction network,
  which are endothelial by definition. Built and measured (§8C.5): contact with
  the *mask* calls 93% of nuclei endothelial (232 of 276 overlap a mask voxel
  outright — [P19](#p19)'s speckle), and restricting to the junction network
  makes it selective but finds **no displacement on any image**, with the sign
  consistently backwards. The reason is not a threshold: fibroblasts sit on the
  outside of the very wall the reference is defined by, so the reference and the
  test populations are not independent of the axis they are compared on. An
  independent endothelial label is what is needed — [P14](#p14)'s unused
  junction-probability map is the nearest candidate in the repo, and
  [I02](#i02) is still the entry for it.

- 2026-08-11 — **this entry got wider, not narrower, when [P04](#p04) was
  fixed.** The classification now reaches the CSV per adhesion as
  `nuc_is_fibroblast`, so an unvalidated call that previously affected only two
  image-level scalars now labels every row. **P04 made the classification
  joinable; it did nothing to make it correct.** Everything above still holds —
  the perinuclear shell reads both cell types where they look most alike, and
  the reported fraction still swings ~2× across one slider step. Do not stratify
  an analysis on `nuc_is_fibroblast` until [I02](#i02) gives it a validation
  signal.

---

<a id="p06"></a>
- 2026-08-25 — **the seeding premise in the entry above is wrong, and correcting
  it points at this classifier over-calling by roughly an order of magnitude.**
  Not measured; queued. Full argument in `docs/DIRECTION.md` §8D.

  The experimenter states the cardiac fibroblasts were **cast into the collagen
  before gelation** (0.5 × 10⁶ cells/mL, 6 mg/mL collagen, fixed at 7 days), not
  seeded onto a finished vessel. Three consequences for this entry:

  - **Casting density predicts ~3–4 fibroblasts per field** (one cell per
    2 × 10⁶ µm³, mean spacing ~126 µm, minus the ~40% of the field the vessel
    occupies) against the 14–41 tabulated above. **This is not evidence on its
    own** — 7 days at a 47–82 h doubling time reproduces the classifier's range,
    and that is an ordinary rate in dense 3D collagen. Recorded as inconclusive.
  - **The depth distribution is the part that discriminates.**
    `docs/DIRECTION.md` §8A.1 found the signed nucleus-to-surface distance
    unimodal, spanning −2.4 to +2.2 µm, with only 2.5% of nuclei more than 5 µm
    abluminal. Cells scattered through a gel at 126 µm spacing should appear at
    every distance out to the edge of the field; there is no such population.
    A single peak at zero of monolayer width is what one wall-bound cell type
    looks like. And **71.4% of the 7 most-abluminal nuclei are called
    fibroblast** — the hit rate climbs steeply with distance.
  - **"Essentially unimodal" is no longer a criticism of the marker.** The code
    comment justifying the background-relative threshold notes perinuclear
    mClover is only ~1.4× background and unimodal. If the true fibroblast share
    is 2–4%, **no marker however good would produce a second mode**, and cutting
    a unimodal distribution at background × 1.5 returns whatever fraction of the
    tail lies above the cut. That is a property of the noise, and it accounts
    for both the ~2× slider swing and the 5.0–14.9% between-image spread without
    requiring the marker itself to be bad.

  **This also weakens the 2026-08-22 corroboration rather than strengthening
  it.** That entry's stated reason for trusting the signed distance —
  "fibroblasts were seeded onto an existing vessel and are abluminal, so the
  geometric measure has a physical reason to work" — rests on the wrong premise.
  The corrected one cuts the other way: cells cast through the gel should be
  separated by distance-from-wall *decisively*, so AUC 0.644–0.773 is a weak
  result on a test that ought to be easy.

  **Status unchanged pending one cheap measurement.** The fibroblasts carry the
  VnTs sensor, so they are directly visible in the clover/ruby channels as whole
  spread cells, isolated at 126 µm spacing and needing no threshold. ~5–10 per
  field confirms the above; ~27 refutes it. Queued as `DIRECTION.md` Tier 3A
  along with a `nuc_is_fibroblast` × signed-distance cross-tab (bounds the
  false-positive rate from below) and a check of whether the fibroblast fraction
  tracks per-image clover noise instead of biology.

  **Not a defect in the shell radius.** Suspected and checked: the dilation
  radius is `round(2.0 / scale[0])` with `scale[0]` the z spacing, which would
  give the 0.9 µm `09step` file a 2.5× thinner in-plane shell. It does not — the
  loader calls `to_isotropic` first, so `scale` is 0.325 isotropic on all six and
  the shell is ~1.95 µm everywhere. Recorded so it is not re-investigated.
- 2026-08-25b — **all three routes to the count have now failed**, so this entry
  is no longer one cheap check from an answer. The eye count assumed the
  fibroblasts would be unmistakable in the sensor channels; the experimenter
  reports they are not distinguishable by eye. Grouping adhesions by nucleus was
  tried as a substitute and **retracted** — a Voronoi artifact that also scores
  the clover photons twice ([P42](#p42)). Endothelial subtraction by Cellpose
  finds zero cells, because the en-face wall has no polygonal tiles
  ([P41](#p41)). See `docs/DIRECTION.md` §8E. What is left is CD31 (an
  independent per-cell endothelial label) and two bounding checks that validate
  nothing — §8D.5 items 2 and 3.

## P06 — Tension readout drifts with imaging depth

| | |
|---|---|
| **Status** | open — **does not replicate; sign flips across images (2026-08-12)** |
| **Severity** | S2 |
| **Opened** | 2026-08-07 |
| **Last verified** | 2026-08-12 (all six; ρ −0.241 to +0.307) |
| **Resolved** | — |
| **Where** | a confound in `ratio_clover_ruby`, not a code defect — **geometric, not optical**; see the 2026-08-12 History entry |
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

> **(a) is discharged, 2026-08-12 — and it came back negative for this entry.**
> Two different things were bundled into one step here. [P03](#p03)'s detection
> truncation replicates unanimously (6/6). **This entry's depth drift does not.**
> See the History entry below: the sign flips across images.

**History**
- 2026-08-07 — recorded in `TODO.md` as a −23% shallow-to-deep quartile drop.
- 2026-08-10 — reproduced, but the quartile series is non-monotonic and ρ is
  −0.187; the "split by wall" mitigation was tested against the data and does not
  hold on this image.
- 2026-08-12 — **measured on all six for the first time. The drift does not
  replicate: the sign flips.** ρ(`ratio_clover_ruby`, `centroid_z_um`) per image,
  on the `20260812_011*` exports:

  | image | n | ρ |
  |---|---|---|
  | `V3_..._04step_002` (reference) | 2,147 | **−0.188** |
  | `V3_..._04step_003` | 2,212 | +0.234 |
  | `V3_..._09step_001` | 3,053 | +0.045 |
  | `V3b_..._04step_001` | 2,019 | **−0.241** |
  | `V3b_..._04step_002` | 4,049 | +0.297 |
  | `V3b_..._04step_003` | 1,872 | +0.307 |

  Range **−0.241 to +0.307**, median +0.140, **negative in only 2 of 6**. The
  reference reproduces this entry's recorded −0.187 to three decimals, so the
  method matches; the other five simply disagree with it.

  **This refutes the optical-attenuation reading of this entry.** mClover
  attenuating with depth while mRuby2 does not is a property of the *optics*, so
  it must push the ratio the same way in every image. It does not. Whatever
  drives the per-image correlation is not primarily attenuation.

  **It confirms the reinterpretation this entry already contains.** The "Why it
  matters" paragraph established that the vessel dives through the stack, so
  depth is substantially *where along the vessel you are* (ρ(z,
  `dist_branch_um`) = −0.461 on the reference). Axial position relative to the
  field of view is arbitrary per image — the vessel can descend or ascend
  through the stack — so a confound of that shape predicts exactly the sign
  inconsistency observed. The six-image result is therefore evidence *for* the
  geometric reading and *against* the optical one.

  **What this changes.** The entry's title overstates the case and is left in
  place for the record; read it as "the tension readout correlates with depth,
  per image, with inconsistent sign". Consequences:

  - **Do not apply a depth correction with a fixed sign or magnitude.** There is
    no common coefficient to fit; the reference image's −0.187 is not the
    dataset's behaviour.
  - **Pooling the six images across depth would cancel much of this to nothing**
    — which would look reassuring and would mean nothing.
  - The per-image pooled p-values are all tiny (1e-2 to 1e-83) and should be
    ignored: [Analysis notes](#analysis-notes) already states that a pooled
    p-value over thousands of adhesions reports adhesion count. Six ρ values
    with inconsistent sign is the informative summary, not any of the p's.
  - **[P03](#p03) is now the better-evidenced depth problem of the two.** It
    replicates 6/6 where this does not, and it is a plausible generator of a
    per-image ratio-vs-depth correlation of arbitrary sign, since which
    adhesions survive detection varies with where the bright wall sits in each
    field.

---

<a id="p07"></a>
## P07 — `peak_local_max` drops border seeds, in two places

| | |
|---|---|
| **Status** | fixed |
| **Severity** | S2 |
| **Opened** | 2026-08-10 (audit F5) |
| **Last verified** | 2026-08-10 (scikit-image 0.26.0) |
| **Resolved** | 2026-08-10 (commit `af8e6c0`) |
| **Where** | `pipeline/gui/nuclei.py:213`, `pipeline/rf.py:606,619` |
| **Affects** | `n_nuclei_total`, `n_fibroblasts`, adhesion instance boundaries (**not** `nuc_id` / `dist_nuc_um` — see History) |
| **Evidence** | `docs/audit.md` §3 F5 |

`peak_local_max` defaults to `exclude_border=True`, which excludes peaks within
`min_distance` of the array edge. Verified on the installed scikit-image 0.26.0:
with `min_distance=5`, a peak at index 3 is dropped and reappears only with
`exclude_border=False`.

**(a) Nucleus splitting** — `gui/nuclei.py:213` uses `min_distance ≈ 12` voxels
(3.9 µm at 0.325 µm/voxel). Nuclei whose distance peak falls within 12 voxels of
any face of the volume get **no watershed seed**, and a mask component with no
marker is left at label 0 by `watershed(..., mask=binary)` — so those nuclei are
dropped entirely, or absorbed into a neighbour that does have a seed. On a
324-plane stack the excluded z shell is 7.4% of the depth, and it sits exactly
where the near-wall and far-wall adhesion populations are.

**(b) Adhesion splitting** — `rf.py:611` crops each candidate to its bounding box
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
- 2026-08-11 — **regenerated for `V3_..._04step_002`.** `fibroblast_count.npz`
  was recomputed by driving the real `state["run_fibroblasts"]` headlessly (the
  GUI nuclei tab is the only implementation of the count, so re-implementing it
  headless would have re-created [P11](#p11)); `auto_min_size_otsu` returned 725,
  reproducing the shipped `dapi_min_size_vox` exactly, and the mClover threshold
  came out at the same 832.5.

  | | before | after |
  |---|---|---|
  | `n_nuclei_total` | 273 | **276** |
  | `n_fibroblasts` | 27 | **27** (unchanged) |
  | "other" (non-VnTs) | 246 | **249** |

  **All three recovered border nuclei classified as non-fibroblast**, so the
  headline `n_fibroblasts` did not move on this image. That is one image and says
  nothing about the general case — the border shell is 7.4% of the depth and
  there is no reason a recovered nucleus should preferentially be endothelial.
  Do not read it as "P07(a) does not matter".

  Both halves are in the re-export
  `..._adhesion_measures_20260811_105032.csv`: `n_nuclei_total = 276` and 2,147
  adhesion rows. **Only this one image has been regenerated** — the other five
  have no `_saved` intermediates at all, so they carry no post-fix export.
- 2026-08-10 — **"Affects" corrected.** The original list named `nuc_id` and
  `dist_nuc_um`. Those come from the *unsplit* DAPI label map: the watershed
  result `nuc` is never written to `state` (`gui/nuclei.py` writes only
  `fibroblast_count`, `:248`), so the seeds that (a) restores never reach them.
  P07(a) affects `n_nuclei_total` and `n_fibroblasts` only; the `nuc_id` half of
  the original claim belongs to [P04](#p04).
- 2026-08-11 — **that correction expires with [P04](#p04).** The watershed result
  *is* now written to `state` and saved as `nuclei_split.npz`, and Measure reads
  it. So from this commit (`04262b8`) onward P07(a) **does** reach `nuc_id` and
  `dist_nuc_um`, and the original "Affects" list is right again for any CSV
  exported after it. The correction remains accurate for every CSV exported
  before it. This is the reason `nuc_id` is not comparable across the P04
  boundary — see the box in [P04](#p04).

---

<a id="p08"></a>
## P08 — `path_curv_xy`: a 20 µm fit window on a ~140 µm-wide vessel, plus contamination and unsigned κ

| | |
|---|---|
| **Status** | open |
| **Severity** | S1 (raised from S2 on 2026-08-11 — see (e)) |
| **Opened** | 2026-08-06 |
| **Last verified** | 2026-08-13 ((c)'s fix landed; the κ table below is **pre-fix** and has not been re-measured) |
| **Resolved** | — ((c) has a landed fix, unverified) |
| **Where** | `pipeline/centerline.py:166` (radius), `:177` (unsigned), `:146-153` (the border pad) — **the file was deleted 2026-08-22**; see `pipeline/lumen_graph.py` for what replaced it |
| **Affects** | `path_curv_xy` |
| **Evidence** | `docs/audit.md` §3 F6 |

> ### Superseded 2026-08-22 — the code this entry describes is deleted
>
> `pipeline/centerline.py` is gone. The in-plane centerline was replaced by
> `pipeline/lumen_graph.py`, a 3D graph built on the Poisson lumen rather than
> on a max projection of the cadherin mask — see `docs/DIRECTION.md` §8B for the
> decision and §8C for the measurements. There is no `binary_fill_holes`, no
> `binary_closing`, no `skeletonize` of a projection and no border pad anywhere
> in the replacement, so the mechanism this entry describes cannot occur.
>
> **Left open as a record, not as work.** The entry is not withdrawn — what it
> measured was true of the code that existed — and the replacement has its own
> defects, recorded as [P32](#p32).

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

> **Superseded by what landed, 2026-08-13 — a third option, and it discards
> nothing.** Neither of the two above was implemented. Instead the footprint is
> **replicate-padded** before `skeletonize` and the pad cropped away afterwards
> (`centerline.py:146-153`), so every arm extends straight past the edge and the
> medial axis carries on instead of forking toward the corners of the cut. The
> pad depth is read off the footprint's own distance transform (p95 of the
> half-width, 81.7 µm on the reference image) rather than chosen, because the
> region over which the border can bend the axis *is* the footprint half-width.
> **The κ table above is pre-fix and has not been re-measured** — see the
> History entry.

**(d) The interim workaround silently empties a non-branched image.** Some
vessels in this dataset are a single tunnel or a cone rather than a Y. With no
branch node, `dist_branch_um` is left all-NaN — correctly and without crashing,
both paths guard with `if len(branch):` (`gui/measure.py:265`,
`measure_core.py:167`). But `dist_branch_um > 20` is False for NaN, so applying
the workaround below to such an image **drops every row** rather than none.

**(e) The fit window is a fraction of the vessel's own width, so κ measures
medial-axis wobble.** This is more fundamental than (a)–(d) and largely subsumes
them. Measured on `V3_..._04step_002`, independently of the mesh:

| | |
|---|---|
| wall z extent in a 5 µm slab ⟂ to the local axis | **~86 µm** (saturates; see [P24](#p24)) |
| xy inscribed half-width at centerline points | median **70.6 µm** → footprint ~140 µm wide |
| total skeleton length | ~640 µm (643 points at ~1 µm spacing) |
| `smooth_um` default | **20 µm** |

The vessel is therefore about **4.5 of its own widths long**, in a field only
~2.4 widths across, and κ is fitted over a window under **a third of the vessel's
half-width**. A medial axis extracted from a blob that wide is dominated by
irregularities in its boundary — including where the field of view cut it, which
is exactly the mechanism behind (c). Fitting κ at 20 µm on it measures how much
the medial axis wobbles, not how much the vessel bends.

**Going to 3D does not fix this** (checked 2026-08-11, see History): an object of
these proportions has no well-defined 1D axis at the 20 µm scale, so the problem
is the definition of the variable, not the projection that produced it. Any real
fix has to come from [P21](#p21)/[I03](#i03) — a reconstruction that yields a
defensible axis — or from abandoning a single-scale κ in favour of one reported
with its scale attached.

**Interim workaround:** filter `dist_branch_um > 20` yourself when analysing the
CSV — but guard it, e.g. `dist_branch_um.isna() | (dist_branch_um > 20)`, or the
non-branched images vanish. See (d). Note this does nothing about (e).

**History**
- 2026-08-06 — recorded in `TODO.md` as "path curvature is undefined at branch
  nodes but still reported".
- 2026-08-10 — audit established the contamination radius (20 µm, not
  point-local) and added the unsigned-κ limitation.
- 2026-08-10 (later) — border flare (c) measured and found to be a larger
  contributor to the high-κ tail than the junction, revising the audit's reading.
  The NaN trap (d) added after the user noted that not every vessel is Y-shaped.
- 2026-08-11 — **(e) added and severity raised S2 → S1**, after the user asked
  why the centerline is built from a 2-D projection at all. The investigation set
  out to test whether the projection costs anything and **could not answer that
  question**, which is itself the finding: there is no trustworthy 3D vessel axis
  to compare against.

  Two attempts at one, both on `V3_..._04step_002`, judged by residual jitter
  about a 20 µm local mean (the κ fit window is 20 µm, so jitter at that scale is
  indistinguishable from bending):

  | Z estimator | \|dz\| to nearest xy neighbour | residual about 20 µm mean |
  |---|---|---|
  | column mean — [P22](#p22)'s `_path_z_um` rule | 2.43 µm | **24.4 µm** |
  | centroid of a 5 µm slab ⟂ to the local xy tangent (median 6,892 voxels/slab) | 1.75 µm | **9.2 µm** |

  The column mean fails structurally: a column through the projection of a
  *hollow* vessel contains **both walls** — median **2 separated z-lobes spanning
  a median 85 µm** of a 105 µm stack — so the count-weighted mean sits in the
  lumen and swings with whichever wall dominates. Weighted median, wall midpoint
  and a robust p5/p95 midpoint were also tried; all leave 20–31 µm of residual.

  The slab centroid is much better and still not enough. κ fitted on it comes out
  **flat across every distance-to-edge band** (0.030 / 0.028 / 0.026 / 0.023 1/µm
  going inward) against the 2D column's 10× variation (0.039 / 0.009 / 0.007 /
  0.004) — a jitter floor, not a measurement. Spearman(κ_2D, κ_3D) = +0.25, which
  is **not** reportable as "the projection destroys the column", because the 3D
  number has no signal in it.

  Chasing why the slab looked wrong produced the caliber measurements in (e), and
  through them [P24](#p24). Ordering note: (c) border flare is still worth fixing
  on its own — it is cheap and it removes a known artefact — but with (e) on the
  record it should no longer be described as making the column trustworthy.

- 2026-08-13 — **(c) has a landed fix, by border padding rather than by
  discarding points.** `vessel_path_curvature_xy` gained a `border_pad_um`
  argument (default: auto, from the footprint's own distance transform) which
  replicate-pads the footprint before `skeletonize` and crops the pad after. The
  effect recorded at the bench, on the reference image: the lower arm had run
  straight down and then kinked hard toward the bottom-right corner, and the
  upper-left arm had forked into two short branches at the top edge; both are
  gone, and 80 µm of pad is where the lower arm measured straight.

  **This is not yet verified as a fix to the defect this entry states.** What was
  observed is the skeleton's *shape*; what (c) measures is the **κ distribution
  by distance-to-edge** (10× the interior curvature in the outer 5 µm, 42% of
  top-decile κ within 30 µm of an edge). That table has not been recomputed, and
  it is the acceptance test. Until it is, treat (c) as *fix landed, unverified* —
  re-running it is one `batch_mesh --steps path` per image plus the same
  banding as the table above.

  **(a), (b), (d) and (e) are untouched.** In particular (e) still stands and
  still dominates: a border-clean medial axis fitted over a 20 µm window on a
  ~140 µm-wide vessel is still measuring medial-axis wobble. Removing a known
  artefact does not make the column trustworthy, exactly as the ordering note
  above says.

  Also landed alongside, and **not** a fix to anything in this entry: the
  centerline overlay now carries a Z ([P22](#p22)) and an inscribed lumen radius
  ([P24](#p24)(b)). Neither enters κ, which is still computed purely in xy.
  A new limitation was recorded while doing it — [P27](#p27).

  > **Correction, 2026-08-13 (same day) — the padding fix is directional, and it
  > fails at corners.** The entry above says the fix removes the border forks.
  > Measured: it removes the one on the *lower* arm, which meets the bottom edge
  > close to square. It does **not** remove the fork where an arm exits near a
  > **corner** — there the footprint touches two edges, and replicate-padding
  > extends each one *perpendicular to itself*, producing one arm in +x and one
  > in +y instead of a single arm running into the corner. Fork present with and
  > without the pad, at identical coordinates. Full measurement and the real
  > cause (a border-connected notch in the mask) in [P29](#p29).

---

<a id="p09"></a>
## P09 — `dist_mesh_um` measures to rejected geometry; curvature has no distance gate

| | |
|---|---|
| **Status** | open — **(a) fixed 2026-08-12**; (b) and (c) open |
| **Severity** | S2 |
| **Opened** | 2026-08-06 |
| **Last verified** | 2026-08-12 (all six images) |
| **Resolved** | (a) 2026-08-12; (b), (c) — |
| **Where** | `pipeline/saved_check.py::surface_verts_for_distance` (the fix), `pipeline/measure_core.py:132`, `pipeline/gui/measure.py:223` |
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
- 2026-08-12 — **(a) is much larger on the other five images; "small numerical
  effect" was a property of the reference mesh, not of the defect.** All six
  measured the same way (adhesion centroids from `adhesion_instances.npz`
  against a KD-tree of all verts vs one of surface verts only). The reference
  row reproduces this entry's original figures — 15.1% and p90 0.36 µm against
  the recorded 0.37 — so the rows below are comparable to it:

  | image | interior verts | adhesions whose nearest vertex is interior | `dist_mesh_um` understated: median / p90 / max |
  |---|---|---|---|
  | `V3_..._04step_002` (reference) | 18.8% | **15.1%** | 0.00 / **0.36** / 5.26 µm |
  | `V3_..._04step_003` | 82.2% | 60.3% | 0.18 / 1.61 / 4.07 µm |
  | `V3_..._09step_001` | 84.9% | 61.3% | 0.25 / 1.98 / 5.34 µm |
  | `V3b_..._04step_001` | 85.3% | **75.4%** | 0.57 / 1.88 / 23.25 µm |
  | `V3b_..._04step_002` | 72.9% | 39.1% | 0.00 / **5.48** / 52.87 µm |
  | `V3b_..._04step_003` | 79.8% | 39.3% | 0.00 / 2.01 / **127.87** µm |

  The mechanism is [P26](#p26)/[P20](#p20): at `subsample_um = 1.0` the point
  cloud is dense enough that the alpha shape buries **73–85%** of its vertices
  inside the solid, against 18.8% on the reference's 5 µm cloud. `dist_mesh_um`
  queries all of them, so on five of six images the majority of adhesions
  measure their distance to a point **inside** the wall rather than on it.

  This changes the entry's severity in practice while leaving the mechanism as
  described. It is no longer fair to call (a) "a semantics defect, not a large
  error" — that sentence is true of the reference image only, and is left in
  place above with this note against it.

  **The fix needs no re-batch, but it is not one line.** Build the distance tree
  from `verts[np.unique(faces)]` instead of `verts`, at `measure_core.py:132`
  and `gui/measure.py:223` — twice, per [P11](#p11). The obstacle is that
  `compute_measurements` **is never given the faces** (signature at
  `measure_core.py:45`: it takes `vessel_verts`, not the mesh), so this is a
  signature change to the shared core plus both call sites, not an edit in
  place.

  **Do not reuse the existing `surf_mask`** at `measure_core.py:189` /
  `gui/measure.py:292`, which is tempting because it is right there. It is
  defined as `np.isfinite(k1_all)` — vertices where *curvature* was computable —
  so it is gated on `surf_ok` and silently absent whenever the curvature file is
  missing or stale, which is one of the degradations [P12](#p12) exists to
  report. Faces come with the mesh and are always present, and
  `np.unique(faces)` is the same definition `mesh_curvature` itself uses
  (`cell_mesh.py:404`). Keep the two sets derived from the same source.

  Backward compatibility: default the new parameter to `None` and fall back to
  all vertices with a note, so a caller that has verts but no faces degrades
  loudly rather than changing answers silently. Expect `dist_mesh_um` to *rise*
  on most rows, by the medians in the table.

- 2026-08-12 — **(a) fixed, and all six images re-measured.** New pure function
  `saved_check.py::surface_verts_for_distance(verts, faces)` returns the
  face-referenced vertex subset plus notes; `compute_measurements` takes a new
  `vessel_faces=None` and builds the distance tree from it. It lives in
  `saved_check.py` rather than in either worker for the reason that module
  exists ([P12](#p12)): the GUI and the batch emit the same string or they
  drift. Both call sites updated — `measure_core.py` and `gui/measure.py`
  ([P11](#p11)'s two-edit tax, paid again) — plus `batch_measure.py`, which
  loads `faces` from the mesh npz guarded by `"faces" in _vm`.

  **Safe because `mesh_nn_idx` is dead.** Both workers unpack
  `mesh_dist_all, mesh_nn_idx = f_mesh.result()` and never use the index, so
  restricting the tree changes distances only. Nothing is looked up by vertex
  index, which is what would otherwise have made this a silent reindexing bug.

  **Measured on every image.** The recomputation below reproduces each exported
  `dist_mesh_um` to **3.6e-15 µm**, so the "before" column is the exact
  counterfactual, not an estimate:

  | image | interior verts | adhesions moved | `dist_mesh_um` median before → after | Δ p90 | Δ max |
  |---|---|---|---|---|---|
  | `V3_..._04step_002` (reference) | 18.8% | 14.7% | 2.106 → **2.298** | 0.30 | 5.01 |
  | `V3_..._04step_003` | 82.2% | 50.5% | 1.028 → **1.489** | 0.98 | 2.92 |
  | `V3_..._09step_001` | 84.9% | 52.2% | 0.727 → **1.379** | 1.22 | 4.71 |
  | `V3b_..._04step_001` | 85.3% | **62.5%** | 0.563 → **1.300** | 1.17 | 21.85 |
  | `V3b_..._04step_002` | 72.9% | 34.5% | 4.275 → **4.550** | 5.71 | 51.92 |
  | `V3b_..._04step_003` | 79.8% | 34.7% | 2.620 → **2.925** | 1.22 | 127.96 |

  On the five batched images the median **rises 1.4–2.3×** — the old column was
  not slightly short, it was reporting the distance to whatever interior point
  of the solid happened to be nearest, and on a dense point cloud that is
  always close. The reference moves least (14.7%) because its 5 µm subsample
  leaves the fewest interior points, which is exactly why the defect looked
  minor when it was only ever measured there.

  Worth noting for [P26](#p26): the fix also *narrows* the between-image spread
  it describes. Pre-fix medians ran 0.56–4.28 µm (7.6×); post-fix they run
  1.30–4.55 µm (3.5×). Part of what looked like a reconstruction mismatch was
  this defect tracking point density. The rest is real and [P26](#p26) stands.

  **Verified against the one controlled before/after available.** The reference
  is the only image with a pre-fix post-[P01](#p01) CSV. Re-measuring it
  (`20260811_105032` → `20260812_011205`, 2,147 rows both) leaves **21 of 24
  columns byte-identical**, including all four `curv_*` — the check that the
  curvature path, which keeps its own `surf_mask`, was not disturbed. The other
  two movers are not this fix: `dist_nuc_um` (4.237 → 4.875, 310 rows) and a
  renumbered `nuc_id` reproduce [P04](#p04)'s recorded figures exactly, because
  that CSV predated the nucleus fixes.

- 2026-08-13 — **(a) cannot arise on the [I06](#i06) Poisson method.** It is
  marching-cubes output, so every vertex is a corner of at least one triangle:
  measured on the reference image, **1.2%** interior vertices against alpha 3D's
  **79.9%** at defaults and 18.8% on the hand-tuned reference. (The residual
  1.2% is vertices whose faces `min_component_um` cut, not buried points.) The
  fix recorded above is still what makes `dist_mesh_um` right for the alpha
  methods; this is a note that one reconstruction does not have the defect to
  begin with.

  (b) is *improved but not gone*: surface-vertex spacing is **0.92 µm** on the
  Poisson mesh against 2.25 µm on the reference, so `dist_mesh_um` is a distance
  to a finer point lattice — still a lattice.

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
| **Where** | `pipeline/measure_core.py:228`, `pipeline/gui/measure.py:335` |
| **Affects** | `force_total_pN` |

`force_total_pN = force_mean_pN × props["area"]`, and `props["area"]` is a volume
in µm³ (`regionprops_table` was given `spacing=scale`). The column's units are
pN·µm³.

Currently blank on this dataset because FRET is not configured ([P18](#p18)), so
it misleads no one yet — but it will the moment a `fret` block is added. Fix the
name, the units, or both before that happens.

> **Priority note, 2026-08-13.** The 2026-08-13 scope decision — build the
> geometry pipeline now, in anticipation of a FRET calibration that does not yet
> exist — makes this **the one thing on the force side worth fixing today**. The
> rule it follows from: fix anything that would be *silently wrong the moment
> force is populated* while it is still blank. This is the only entry that
> qualifies. ([P02](#p02)'s `background_subtract` clipping will also apply to the
> FRET channels when they arrive, but that is frozen code and needs a decision at
> calibration time, not a fix now.)

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

**History**
- 2026-08-06 — opened.
- 2026-08-11 — **the duplication got measurably worse, and the cost is now
  quantified.** [P04](#p04) and [P12](#p12) each had to be written twice: the
  nucleus-map selection, the `nuc_is_fibroblast` join, and the notes plumbing all
  exist once in `measure_core.py` and again in `gui/measure.py::_run`. That is
  three more paired edits on top of the ROI divergence and the sum-then-divide
  ratio fix already on the record.

- 2026-08-22 — **the trend was stopped rather than continued.** Replacing the
  centerline with the 3D graph added eleven new per-adhesion columns, which
  under the old pattern would have been a fourth paired edit. Instead the whole
  block was extracted as `measure_core.graph_columns`, and both
  `measure_core.compute_measurements` and `gui/measure.py::_run` call it. The
  ROI divergence and everything else in this entry are untouched, so it stays
  open — but the newest code is shared rather than duplicated.

  Partially mitigated rather than fixed: the *checks* were extracted into
  `pipeline/saved_check.py`, which both workers call, so the class of divergence
  P12 was about ("the GUI warns, the batch does not") cannot recur silently. The
  *measurement* is still duplicated. Extracting it remains deferred, and the
  argument for doing it is now stronger than it was — each fix to this pipeline
  costs two edits and risks a third divergence.

---

<a id="p12"></a>
## P12 — Headless measure fails silently and overwrites saved instances

| | |
|---|---|
| **Status** | fixed |
| **Severity** | S3 |
| **Opened** | 2026-08-10 (audit F9) |
| **Last verified** | 2026-08-11 |
| **Resolved** | 2026-08-11 (commit `04262b8`) |
| **Where** | `pipeline/batch_measure.py:181`, `:144`; `pipeline/gui/vessel_mesh.py:515` (`_load_curv`'s length check) |

Two asymmetries between the GUI and headless paths:

- `gui/vessel_mesh.py:515` prints an explicit warning when a loaded curvature
  file's length does not match the mesh, because Measure will silently NaN every
  curvature column. `batch_measure.py` loads the same two files with **no such
  check**, so a stale pairing produces a CSV with four blank columns and no
  message.

  > **Reference correction, 2026-08-11.** This entry cited `gui/vessel_mesh.py:501`
  > for that warning. The line number was already stale when the entry was
  > written; the check is at `:547-552`.

- `_standard_instances` **overwrites `adhesion_instances.npz`** as a side effect
  of measuring. It is documented in the docstring and the reasoning is sound
  (a folder may hold a stale, non-standard labelling), but a command named
  "measure" silently rewriting a saved segmentation deserves a flag rather than
  being the default.

**The underlying problem is structural, not two missing checks.** The GUI had
grown nine of these checks inline and the headless path had none, because each
was written at the point of use. Adding them inline a second time would have
recreated the same divergence one release later — the same mechanism that
produced [P11](#p11).

**The fix, 2026-08-11 — one module both paths call.** New `pipeline/saved_check.py`:
pure functions (numpy in, `list[str]` out; no Qt, no napari, no `state`, no
printing) that answer "which columns come out blank, and why?". `measure_core`
accumulates them into `result["notes"]` and **still prints nothing** — it does not
know its caller's indent and is the one module a non-CLI caller could reuse — and
each caller prints at its own indent. `gui/vessel_mesh.py` now calls
`curvature_notes` instead of its inline warning, so both paths emit *the same
string*, verified character-for-character.

Notes were added for every degradation that was previously silent:

| column(s) | silent cause now reported |
|---|---|
| `curv_*` (all four) | no curvature file / length mismatch / all-NaN `k1` |
| `curv_<one>` | the per-key branch, which can blank one column while the others survive |
| `path_curv_xy`, `dist_branch_um` | no centerline — *renamed 2026-08-22; see below* |
| `dist_branch_um` | centerline with no branch nodes (+ the [P08](#p08)(d) NaN-filter trap) |
| `ratio_clover_ruby` | **native clover/ruby vs labels shape mismatch** — the project's primary output, previously blanked with no message at all |
| `ratio_clover_ruby` (per row) | the [P15](#p15) `r_sum = 0` dropout, now reported as a count, so a depth-dependent effect is visible |
| `force_*_pN` | no force volume (previously the *only* announced one; moving it here makes the set uniform) |
| `nuc_is_fibroblast` | no/inconsistent `nuclei_split.npz` ([P04](#p04)) |

> **2026-08-22 — those two rows now read `vessel_graph.npz`.** The columns were
> renamed with the centerline rewrite, and `saved_check.centerline_notes` became
> `graph_notes`. It warns for a missing graph (every geometry column blank), for
> a graph with no junction (`dist_junction_um` blank — the same NaN-filter trap,
> which is a property of NaN and outlived the code), for endpoints that stop
> inside the volume where no real terminus is expected, and for a run where more
> than a quarter of cross-sections were clipped, which makes the caliber a lower
> bound rather than a measurement.

Also in this fix:

- **The instances rewrite is now opt-in** (`--rewrite-instances`).
  Recompute-from-proba stays the default for correctness on legacy folders, but
  without the flag the recomputed labels are *compared* against the file and a
  difference is reported, not applied. Verified to be a **no-op for
  freshly-batched images**: `run_rf_postprocess` forces `split_merged=True`
  (`gui/adhesion_rf.py:550`) and the widget defaults (0.5 / 100 / 1000) are
  byte-identical to `STD_POSTPROCESS`, so the previous behaviour was rewriting a
  file it had just correctly written. Measured on `V3_..._04step_002`: 2,147
  adhesions recomputed, 2,147 on disk, npz left untouched.
- **`--check-only`**, the pre-flight for an unattended run: `describe_saved` over
  every `*_saved` folder, printing what is present and which columns would come
  out blank — **before** the compute rather than after. Cheap by construction
  (an npz is a zip, so only small members are ever decompressed; array shapes are
  read from `.npy` headers without inflating the payload). Measured at **1.0 s**
  over a folder containing a 340 MB probability map.
- **Explicit messages on the two previously silent skips** (no instances;
  already measured).
- **Exit codes.** `batch_measure` returns 1 if any image FAILED and, with
  `--strict`, 2 if any CSV was written with a silently blank column. A
  degradation warns but does not fail by default — a legitimately FRET-less
  dataset (this one) must not fail the run. `batch.py` tracks per-image outcome
  in `_advance` (already the single idempotent funnel), prints
  `N processed, N skipped` with reasons, and exits non-zero if any image was
  skipped.
- **`fibroblast_count.npz` can no longer be a foot-gun.** `int(fc[k])` raised
  `TypeError` on any non-scalar member, turning a slightly-different file into a
  hard failure of the whole image. Non-scalar members are now reported and
  skipped.

**History**
- 2026-08-10 — opened from `audit.md` F9.
- 2026-08-11 — **fixed.** Verified end to end on a copy of
  `V3_..._04step_002_saved` with `vessel_curvature.npz` truncated to 500 of 8,966
  values: `--check-only` reports it in 1.0 s without measuring; the GUI's
  `Load curvature` and the headless measure emit **the identical string**; the
  run reports `1 degraded`; `--strict` exits 2 while the plain run exits 0.
  Scope note: the uppercase `WARN:` / `FAILED:` convention was applied only in
  the files this change touches, not swept across the repository.

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
| **Where** | `pipeline/gui/nuclei.py:258`, `pipeline/rf.py:619` |

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

- [x] ~~`centerline.py:41` — `_prune_spurs` mutates `skel` while iterating the
      endpoint list, but the neighbour set `S` was built before the mutation, so
      paths can be traced through already-deleted pixels within one pass. It
      converges, but the pruned result depends on iteration order.~~
      **Moot 2026-08-22: the function and its module are deleted.**
      `lumen_graph.prune_spurs` removes a whole leaf edge at a time from the
      graph and rebuilds the degrees after each removal, so there is no
      partially-mutated skeleton to trace through and no order dependence.
- [ ] `measure_core.py:114` — adhesions thinner than one native z-step can fall
      entirely between sampled planes, giving `r_sum = 0` → NaN ratio. Correctly
      NaN'd rather than wrong, but it is a size-dependent (therefore
      depth-dependent) dropout that compounds [P03](#p03).
      **Now instrumented and measured (2026-08-11):** [P12](#p12) reports the
      dropout as a count, and on `V3_..._04step_002` it is **0 of 2,147
      adhesions** — the mechanism is real but is not biting on this image at the
      standard post-process (`min_size_vox = 100`, so surviving adhesions are
      large enough to intersect a sampled plane). Worth re-reading on the five
      un-batched images, and on `09step`, whose 0.9 µm z-spacing is ~2× coarser
      and is exactly the condition that would trigger it.
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
| **Last verified** | 2026-08-12 (the counterfactual sweep — the floor sits mid-plateau) |
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

**History**
- 2026-08-10 — opened.
- 2026-08-12 — **the cut is now recorded and overridable, and the "never
  checked" note above is half answered.** `pipeline.batch_mesh` takes
  `--min-size <int>|auto|auto-raw` (`auto` = `max(otsu, 12)`, exactly what
  `gui/cadherin.py` applies; `auto-raw` = the Otsu cut with no floor — the
  experiment this entry asks for) and writes the resolved value, the raw Otsu
  value and how it was chosen into `vessel_recon.json`. So the size cut is no
  longer absent from the record.

  Running it across all six saved masks, `auto_min_size_otsu` returns:

  | image | Otsu cut | components ≥ 12 | mask voxels at 12 |
  |---|---|---|---|
  | `V3_..._04step_002` | **3** | 3,764 / 226,093 | 1,517,396 |
  | `V3_..._04step_003` | **2** | 3,330 / 329,511 | 2,056,144 |
  | `V3_..._09step_001` | **2** | 5,841 / 265,173 | 2,715,578 |
  | `V3b_..._04step_001` | **2** | 3,473 / 356,494 | 1,746,715 |
  | `V3b_..._04step_002` | **4** | 2,396 / 69,482 | 423,197 |
  | `V3b_..._04step_003` | **2** | 2,512 / 251,045 | 435,300 |

  The auto value is **2–4 on every image**, so the floor of 12 is 3–6× it
  everywhere and the operative parameter is the constant in all six cases, not
  just the one measured in 2026-08-10. The entry's claim generalises.

  **Still open**, and nothing here changes it: the criterion is still Otsu on a
  distribution that is not bimodal. Recording a number is not choosing a better
  one. `auto-raw` now makes the counterfactual runnable — nobody has run it yet.
- 2026-08-12 — **the counterfactual has now been run, and the magic constant is
  doing real work.** `min_size` swept with the mesh knobs pinned at the reference
  set, every mesh scored by `mesh_score` against **one fixed mask**
  (`min_size = 12`) so a filter that discards real wall cannot score as an
  improvement by being graded on the subset it kept:

  | `min_size` | components kept | mask voxels | surface verts | χ | coverage ≤2 µm | support max reach |
  |---|---|---|---|---|---|---|
  | **0** (none) | 226,093 | 1,833,005 | 3,478 | **−3** | **2.5%** | **78.27 µm** |
  | **3** (≈ the Otsu cut) | 21,808 | 1,595,609 | 8,075 | −289 | 40.3% | 23.43 µm |
  | **12** (the floor, shipped) | 3,764 | 1,517,396 | 7,282 | −248 | 42.9% | 13.72 µm |
  | 30 | 2,153 | 1,488,270 | 6,907 | −223 | 43.3% | 12.82 µm |
  | 100 | 1,069 | 1,429,001 | 6,183 | −209 | 42.4% | 12.61 µm |
  | 500 | 321 | 1,269,637 | 4,781 | −137 | 36.7% | 10.66 µm |

  Three things follow, none of them previously measurable:

  - **Otsu alone is materially worse than the floor.** At `min_size = 3` the
    reconstruction carries **23.4 µm** spurs against 13.7 µm at 12, and covers
    less of the wall (40.3% vs 42.9%). So "the floor defeats Otsu" is true as
    stated, but the thing it defeats is worse than it is. Running `auto-raw`
    would degrade this image, not improve it.
  - **12 is not a knife-edge.** From 12 to 100 — an 8× range — coverage stays
    42.4–43.3% and max reach 12.6–13.7 µm. The constant is unjustified but it
    sits in the middle of a broad plateau, which is a much weaker defect than
    this entry assumed.
  - **`min_size = 0` is a cliff, not a slope.** With no filtering at all the
    alpha shape bridges the 226,093 speckle components into one smooth envelope:
    χ = −3, coverage **2.5%**, spurs to 78 µm. Note that this mesh also has the
    *best* `slope_loglog` of any built in this project (**+0.07**, i.e. formally
    "scale-free") — because a featureless blob genuinely is. That is the
    blob-gaming failure [I03](#i03) predicts, observed rather than hypothesised,
    and it is caught only by `coverage`. It is the concrete argument for scoring
    with a vector rather than any single number.

  **Still open**: the criterion. This measures that the shipped value is
  defensible and that the plateau around it is wide; it does not supply a
  principled way to pick one, and the distribution is still not bimodal.
- 2026-08-13 — **a larger problem sits upstream of this one: [P28](#p28).** This
  entry is about which *components* to discard after binarisation. Measured on
  the raw channel, the binarisation itself is discarding real wall — the gaps in
  the mask carry signal at the same radius as the detected wall, ~40% as bright,
  at 5–6× the lumen noise. So tuning the size floor is optimising downstream of
  a lossy step. Nothing here is retracted; the plateau measurement above stands.
  But "the mask is roughly right and the floor is an unjustified constant on it"
  is too generous a reading of the mask.

- 2026-08-21 — **raised to S2: this entry, not the binarisation, is where most of
  [P28](#p28)'s lost wall goes.** [P28](#p28)'s attribution measured which of the
  three gates deletes the dim wall it found, per gap-wall peak, on five images.
  The local threshold fails **0.000** of them, everywhere. The global floor takes
  15–41%. **`min_size` takes 59–85%** — i.e. the majority of the deleted wall
  passes both thresholds, lands in `labels`, and is then discarded as a small
  component by the constant this entry is about.

  That inverts the 2026-08-13 line above. It is *not* true that "tuning the size
  floor is optimising downstream of a lossy step"; the size floor **is** the
  lossy step, and the binarisation upstream of it is keeping the signal. The
  plateau measurement still stands, and so does the complaint that the value is
  an invisible unrecorded constant — but its stakes are now much higher than
  "which speckle to discard".

  **The counterfactual is harder than it looks, and that is also new.** Lowering
  `min_size` cannot be justified on this evidence alone: the components it drops
  have a mean local z-score of 2.42–3.91 against 3.02–4.97 for the components it
  keeps — **overlapping distributions**, with the gap between them smaller than
  the spread within either. On the last two images the dropped voxels are a
  quarter to three-quarters of what is kept (`V3b_..._04step_003`: 327,269
  dropped against 435,300 kept). So there is no brightness cut that recovers the
  dim wall without also recovering the speckle, which is the same wall this
  entry's plateau measurement has been running into from the other side.
  Evidence: `docs/scratch/p28/rays.py`, `out/summary.json`.

---

<a id="p20"></a>
## P20 — `subsample_um` silently sets the vessel mesh's resolution

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-12 (scored across `subsample_um` 1 → 5; not a knife-edge) |
| **Resolved** | — |
| **Where** | `pipeline/vessel_mesh.py:52` (`voxel_subsample`), `pipeline/gui/vessel_mesh.py:85` (the slider default) |
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

**History**
- 2026-08-10 — opened.
- 2026-08-12 — **`subsample_um` is no longer silent, and its effect is now
  measured in isolation.** Every mesh built by `pipeline.batch_mesh` writes a
  `vessel_recon.json` recording `subsample_um` next to the vertex counts it
  produced, and `describe_saved` warns for any mesh that has no such file. On one
  image with the mask held fixed, moving it 5.0 → 1.0 µm takes the mesh from
  7,282 to 33,537 **surface** vertices and κ2's median from −0.0701 to −0.1632
  1/µm — full table in [P26](#p26)'s History.

  The coupling this entry describes is also enforced now rather than only
  documented: `batch_mesh` warns when `alpha3d_um < 3 × subsample_um`, which is
  the commonest route to "N verts, 0 faces". Verified — at `subsample_um = 1.0`
  with `alpha3d_um = 0.4` it warns, and the build then fails with exactly that.

  **Unchanged:** the vote-count idea in the table above is still not implemented.
  `voxel_subsample` still keeps one representative per occupied cell and discards
  the occupancy, and [P09](#p09)(b)'s resolution limit is still set by this knob
  rather than by the image.
- 2026-08-12 (later) — **this knob is now scored, and it is not a knife-edge.**
  With `alpha3d_um` and `min_component_um` pinned, `subsample_um` 1 → 5 moves
  `mesh_score`'s `slope_loglog` by only 0.185 and its coverage median smoothly
  and monotonically (1.62 → 2.16 µm) — so the reconstruction is not balanced on
  a hair-trigger setting of the dominant knob, which was an open worry rather
  than a measured fact. Full table in [I03](#i03)'s History.

  The same sweep shows the knob is **not** what produces the spurs reaching out
  of the vessel: max fabricated reach stays at 14.26 → 13.72 µm and max edge at
  ~29 µm across the whole range. That is set by `alpha3d_um` ([P21](#p21)), not
  by this. The `subsample_um = 5.0` rebuild is byte-identical to the shipped
  `vessel_mesh.npz`, which independently confirms the reference parameter set a
  second time.

  **Still unchanged:** the vote count. And note that the isolated points this
  entry says the alpha shape bridges to are exactly what `mesh_score.support`
  now measures the consequence of — the evidence for the vote-count idea is
  therefore now collectable, though it has still not been collected.

---

<a id="p21"></a>
## P21 — One global alpha cannot be correct on a tapering vessel

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-12 (the alpha↔spur relation measured on one image; the taper case is still predicted, not measured) |
| **Resolved** | — |
| **Where** | `pipeline/cell_mesh.py:137` (`alpha_shape_faces`), `pipeline/gui/vessel_mesh.py:86` (the slider default) |
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

> **Addendum, 2026-08-11 — the caliber this entry reasons about is now measured.**
> On `V3_..._04step_002` the wall spans **~86 µm in z** in a slab cut
> perpendicular to the local axis, saturating against the search radius, and the
> footprint is ~140 µm wide ([P24](#p24)(b)). A vessel radius of order 40 µm sits
> **far above** the shipped `alpha3d_um = 15.0`, which is consistent with this
> entry's measurement that the lumen is not filled here — alpha is nowhere near
> the local lumen radius, so the two-sheet regime is live, as stated. It also
> means the *other* half of the alpha condition is the binding one on this image:
> alpha is being set by the wall-gap scale (~8–10 µm) with a wide margin below the
> lumen, not squeezed between the two. On a taper the narrow end still closes the
> gap, so the structural argument is unchanged.

> **Addendum, 2026-08-12 — what alpha costs is now measured, and it is what the
> eye sees.** `mesh_score.support` measures how far the reconstructed surface
> sits from any retained cadherin voxel — sampled over triangle *interiors*,
> since the alpha-3D vertices are cadherin voxels and so are at distance zero by
> construction. Across `alpha3d_um` = 5 → 30 on `V3_..._04step_002`, the maximum
> fabricated reach tracks **≈ 0.95 × `alpha3d_um`** and the longest edge
> **≈ 2 × `alpha3d_um`** (the latter exactly as a circumradius bound requires),
> while `subsample_um` moves neither. So `alpha3d_um` is precisely "the longest
> spur out of the wall that I permit", and at the shipped 15.0 that is a
> **13.7 µm** spur with **19.2%** of the surface area sitting more than 5 µm from
> any voxel.
>
> This is the same trade-off this entry describes, now with both sides on one
> axis: alpha must be large enough to bridge the ~8–10 µm wall gaps (below ~7.5
> coverage collapses) and every micron above that is bought as fabricated
> surface. It does not resolve the taper argument — a varying lumen radius still
> sweeps through any fixed alpha — but it does mean the cost of a too-large alpha
> is no longer invisible. Table in [I03](#i03)'s History.

**Related:** [P01](#p01), [P20](#p20), [P24](#p24), [I03](#i03).

**History**
- 2026-08-13 — **`alpha3d_um` as "the longest spur I permit" does not apply to
  the [I06](#i06) Poisson method at all**, because that method has no alpha and
  bridges gaps with a harmonic solution rather than a triangle. Its longest
  triangle edge on the reference image is **3.36 µm** against alpha 3D's
  **19.90 µm** — the fabricated-surface mechanism this entry describes is simply
  absent.

  That does not make its bridging free, and the numbers go the other way on
  reach: `support` max reach is **26.70 µm** for Poisson against **9.87 µm** for
  alpha 3D at defaults (and 13.72 µm for the hand-tuned reference), with
  unsupported area >5 µm at 19.6% against 16.5%. So the Poisson surface reaches
  *further* from the data than the alpha shape does, but as smooth sheet rather
  than as spurs. The taper argument this entry makes is untouched — it is an
  argument about a single global length scale, and `poisson_smooth_um` is also a
  single global length scale, so a tapering vessel sweeps through it too. What
  changes is only that the cost is no longer paid in spikes.

---

<a id="p22"></a>
## P22 — Centerline overlay is drawn at a single wrong z

| | |
|---|---|
| **Status** | fixed — the 2026-08-10 fix was itself **replaced** on 2026-08-13 |
| **Severity** | S3 (display only — no exported column changes) |
| **Opened** | 2026-08-10 |
| **Last verified** | 2026-08-13 |
| **Resolved** | 2026-08-10 (commit `af8e6c0`); superseded 2026-08-13 |
| **Where** | `pipeline/centerline.py:195` (`lumen_center_z`), `pipeline/gui/vessel_mesh.py::_show_path`. **`_path_z_um` no longer exists** — the pre-2026-08-13 fix. **Both the file and `_show_path` were deleted 2026-08-22** |
| **Affects** | the `Vessel path curvature (X-Y)` layer only — **which no longer exists**; see the 2026-08-22 note |

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

> **2026-08-22 — the layer this entry is about is gone.** `lumen_center_z` and
> the whole of `pipeline/centerline.py` were deleted with the 2D centerline
> (`docs/DIRECTION.md` §8B). The replacement has no Z problem to solve: the
> graph's points are 3D positions in the lumen by construction, so there is no
> separate Z rule to get wrong and nothing to draw at a "single wrong z". Left
> `fixed` rather than reopened — the defect was real and the fix worked; the
> code simply no longer exists.

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

- 2026-08-11 — **the entry's own question answered: this Z is a good overlay and
  an unusable axis.** The fix above notes it is "worth doing before any
  3D-centerline work, because it shows whether a 3D axis estimate is
  well-behaved before one is built." It is not. A ±3 µm column through the
  projection of a hollow vessel holds **both walls** — median **2 separated
  z-lobes, median 85 µm apart** — so the count-weighted mean lands in the lumen
  and moves with whichever wall has more voxels in that column. Residual about a
  20 µm local mean is **24.4 µm**; a cross-sectional slab centroid gets that to
  9.2 µm and is still too noisy to fit curvature on. Full numbers in
  [P08](#p08)'s 2026-08-11 History entry.

  **No change to this fix or to `_path_z_um`.** The overlay only has to sit on
  the vessel, and it does — median displacement from the old single-z overlay was
  23.5 µm, all 643 columns resolve, 0.15 s on the GUI thread. This note exists so
  the helper is not later mistaken for an axis estimator, which is a use it was
  never measured for.

- 2026-08-13 — **`_path_z_um` is gone; the Z now comes from a 3D
  distance-from-wall field.** The 2026-08-11 note above concluded the column mean
  was "a good overlay and an unusable axis" and left it in place. It turns out it
  was not a good overlay either. Two measurements taken at the bench when the
  replacement was written:

  | | column mean (`_path_z_um`) | `lumen_center_z` |
  |---|---|---|
  | occupied columns detecting **only one** of the two walls | **86%** | n/a |
  | spatial roughness of the resulting Z | **2.95 µm** | **0.00 µm** |

  The structural reason is the one [P08](#p08)'s History already gives — a column
  through a hollow vessel crosses the near **and** far wall — but with a sparse,
  depth-attenuated signal the two are rarely both detected, so the mean does not
  land in the lumen. It lands on whichever wall was detected, and flickers
  between the bottom and top of the stack from point to point.

  **The replacement** (`centerline.lumen_center_z`) block-reduces the mask to a
  local occupancy density, thresholds it, takes a 3D distance-from-wall field,
  and puts each centerline point where that distance is greatest inside the
  footprint. Three properties are what make it work, each recorded because the
  obvious alternatives were tried and failed:

  - it is **3D**, so a column detecting neither wall still inherits a sensible
    distance from its neighbours;
  - it needs **no watertight seal**. Flood-filling the lumen was tried and
    returns nothing — the wall sheet has genus ≈ 125, so every fill leaks. A hole
    perturbs a distance field smoothly instead of breaking it;
  - the stack has only ~10 µm of headroom above the upper wall, so the distance
    maximum lands in the lumen rather than in the empty space above the vessel.

  **Two things it adds beyond the overlay.** The inscribed **lumen radius** at
  each point — see [P24](#p24)(b), where it is an independent caliber check — and
  an explicit `valid` flag marking points where the lumen is unbounded above or
  below, i.e. **where the vessel leaves the field of view**. Those points are
  drawn (Z extrapolated from the nearest bounded neighbour) and ringed in a
  second layer rather than dropped, because gaps along the path read as a broken
  centerline rather than as a caveat. The Z, the radius and the flag travel with
  the saved `path_curvature.npz` as extra keys, so an older reader is unaffected
  and a file predating them says so on load.

  **Still display-only.** κ is computed in xy and Z never enters it, so no
  exported column moves — the same scope this entry has always had.

  > **Caveat on "roughness 0.00 µm", added the same day.** Z is the `argmax` of a
  > distance field on a **block-reduced grid**, so it is quantized to roughly
  > `4 × voxel` ≈ 1.3 µm. Read 0.00 as "below the quantum", not as exact. Two
  > further limits worth stating before this is used as anything more than an
  > overlay: on an elongated cross-section ([P24](#p24)(b), aspect 1.6–2.2) the
  > distance field is nearly flat along the long axis, so the centre is weakly
  > determined in that direction; and `valid` tests only that the maximum is
  > bounded in **z** and inside the footprint — it does **not** flag proximity to
  > the x or y image edge, which is the relevant boundary at a branch end. For
  > [I06](#i06)'s purposes (which side of the wall is the lumen) none of this
  > matters; for anything finer it does.

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
| **Where** | `pipeline/centerline.py:122` (`open_um`), `:156` (`prune_um`) — **the file was deleted 2026-08-22**; `lumen_graph.prune_spurs` replaces both |
| **Affects** | `path_curv_xy`, `dist_branch_um` on non-uniform vessels |

> ### Superseded 2026-08-22 — the code this entry describes is deleted
>
> `pipeline/centerline.py` is gone. The in-plane centerline was replaced by
> `pipeline/lumen_graph.py`, a 3D graph built on the Poisson lumen rather than
> on a max projection of the cadherin mask — see `docs/DIRECTION.md` §8B for the
> decision and §8C for the measurements. There is no `binary_fill_holes`, no
> `binary_closing`, no `skeletonize` of a projection and no border pad anywhere
> in the replacement, so the mechanism this entry describes cannot occur.
>
> **Left open as a record, not as work.** The entry is not withdrawn — what it
> measured was true of the code that existed — and the replacement has its own
> defects, recorded as [P32](#p32).

> The replacement prunes a leaf branch when it is shorter than `spur_ratio`
> times **the median inscribed radius along that branch**, so the threshold
> scales with the vessel instead of being 4 µm and 30 µm everywhere. That is
> what this entry asks for in its last line ("consider making both thresholds
> scale-aware rather than absolute"), and it means image 6's waist — the
> tightest caliber in the dataset and the highest-shear region in it — is no
> longer amputated by a constant tuned on a wide vessel.

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

**History**
- 2026-08-10 — opened.
- 2026-08-12 — **all four parameters are now reachable, so this is measurable
  for the first time.** `pipeline.batch_mesh` exposes `--path-open-um`,
  `--path-prune-um`, `--path-close-um` and `--path-target-px-um` alongside the
  one knob that had a slider (`--path-smooth-um`), records all five in
  `vessel_recon.json`, and `--steps path` recomputes only the centerline. The
  premise of the entry — that these are "function defaults with no slider" —
  no longer holds headlessly; it still holds in the GUI, where only `smooth_um`
  is exposed.

  All six images now have saved intermediates, so the "no taper image has saved
  intermediates" blocker is gone too. **The measurement itself has not been
  made** — sweeping `open_um` and `prune_um` on a taper image and watching where
  the centerline ends is the next step, and it is now one command per setting.

- 2026-08-13 — **the measurement got a better dependent variable, and one
  confound was removed.** `centerline.lumen_center_z` ([P22](#p22)) returns the
  inscribed lumen radius at every centerline point, so "where does the taper get
  cut off?" can be read as *the caliber at the last surviving point* rather than
  only as the path's endpoint — which is the quantity the entry is actually
  about, since both thresholds are caliber thresholds (`open_um` = 4 µm against a
  radius that is falling). Separately, the border pad ([P08](#p08)(c)) means an
  arm truncated by the field of view no longer forks there, so an amputated taper
  and a clipped one are now distinguishable. **Still not measured.**

---

<a id="p24"></a>
## P24 — `curv_*` has no scale-independent value; `radius_um` sets the answer

| | |
|---|---|
| **Status** | open |
| **Severity** | S1 |
| **Opened** | 2026-08-11 |
| **Last verified** | 2026-08-13 (caliber re-confirmed by a third, independent method) |
| **Resolved** | — |
| **Where** | `pipeline/cell_mesh.py:379` (`mesh_curvature`), the `radius_um` slider in `pipeline/gui/vessel_mesh.py` |
| **Affects** | `curv_mean`, `curv_k1`, `curv_k2`, `curv_gaussian` |

[P01](#p01) fixed the curvature *frame*, and the columns are now internally
consistent — 98.3% consistent κ2 sign, face-normal coherence 0.972. This entry is
the separate question of whether they measure **the vessel**. On the evidence
below they do not: they measure the scale at which they were fitted.

**(a) The implied radius tracks the fit window.** Sweeping `radius_um` on the
shipped `V3_..._04step_002` mesh — same 8,966 verts / 7,282 surface verts, same
`smooth_iter = 10`, geometry untouched:

| `radius_um` | κ2 median (1/µm) | implied R | κ2 sign consistency |
|---|---|---|---|
| **5.0 (shipped)** | **−0.0703** | **14.2 µm** | 98.3% |
| 10 | −0.0606 | 16.5 µm | 98.3% |
| 15 | −0.0277 | 36.2 µm | 98.3% |
| 20 | −0.0171 | 58.5 µm | 98.0% |
| 30 | −0.0097 | 102.9 µm | 97.0% |
| 40 | −0.0074 | 134.7 µm | 96.6% |

For a surface that is genuinely a tube of radius R, **κ2 ≈ −1/R is independent of
the fit scale** while the window stays well below R. Here the implied radius rises
monotonically with the window at roughly 3× and **never plateaus**. That is the
signature of fitting a quadric to a surface rough at every scale: the estimate is
a property of the window, not of the object.

Note that sign consistency stays ≥96.6% throughout. **A consistently-signed
number can still be scale-dependent** — [P01](#p01)'s sign-consistency evidence
does not bear on this, and neither refutes the other.

**(b) The vessel is ~5× larger than the shipped setting implies.** Two
measurements, neither using the mesh or its curvature:

| measurement | result |
|---|---|
| wall z extent in a 5 µm slab ⟂ to the local xy tangent | **~86 µm**, and it **saturates** — 87.1 / 87.4 / 86.8 / 86.4 / 85.8 µm as the in-plane search radius goes 30 / 45 / 60 / 80 / 110 µm |
| xy inscribed half-width of the solidified footprint at centerline points | median **70.6 µm** (p10 21.4, p90 87.8) |

Read the z figure as the robust one: it involves no `binary_fill_holes`, no mesh,
and it stops growing with the search window. The xy figure is an upper bound
because `fill_holes` closes the region enclosed by the Y.

> **Superseded the same day — see the direct measurement below. The three
> figures do not measure the same axis.**
>
> **A third measurement, 2026-08-13 — and it agrees.** `centerline.lumen_center_z`
> ([P22](#p22)) returns the **inscribed radius of the lumen** at each centerline
> point: the distance from the chosen centre to the nearest wall in a 3D
> distance-from-wall field. On the reference image its median is **38.1 µm**,
> against the ~40 µm the z-extent figure above gives. This is the pipeline's only
> direct caliber measurement and it is independent of both rows in the table —
> different input (a density-thresholded wall, not a slab), different operation
> (a distance transform, not an extent), no mesh, no `fill_holes`. **The
> ~40 µm figure is now confirmed by three methods that share no machinery**, so
> the 3–5× gap between it and the shipped κ2 = −0.0703 (R = 14.2 µm) is not an
> artefact of how caliber was measured.

> ### Corrected, 2026-08-13 — the cross-section is an ellipse, and "~40 µm" is its **minor** axis
>
> The paragraph above is wrong about *what* agrees. Sorted by which direction
> each figure is sensitive to, the three measurements are **the short axis twice
> and the long axis once**:
>
> | measurement | value | sensitive to |
> |---|---|---|
> | wall z extent ⟂ to local axis | ~86 µm → **43 µm** half | the **z** direction |
> | inscribed lumen radius (`lumen_center_z`) | **38.1 µm** | the **shortest** direction — an inscribed sphere is capped by the minor axis |
> | xy inscribed half-width of the footprint | **70.6 µm** | the **in-plane** direction |
>
> Measured directly, cutting the enclosed lumen out of the cadherin shell on the
> left arm (`docs/scratch/centerline_corner/`):
>
> | x (µm) | lumen area | y extent | z extent | **aspect** | solidity |
> |---|---|---|---|---|---|
> | 25 | 10,679 µm² | 151.5 | 92.3 | **1.64** | 0.968 |
> | 55 | 10,975 | 152.1 | 91.0 | **1.67** | 0.979 |
> | 70 | 11,067 | 154.1 | 87.8 | **1.76** | 0.979 |
> | 90 | 11,296 | 160.2 | 89.0 | **1.80** | 0.964 |
> | 115 | 12,637 | 179.7 | 83.5 | **2.15** | 0.985 |
>
> So the lumen is ~150–180 µm across in plane and ~85–92 µm deep — a flattened,
> slightly D-shaped ellipse. **Not an imaging artefact:** 0.01% of the mask lies
> in the first z slice and **zero** in the last, 3.1% within 5 µm of either
> limit, so the vessel fits inside the stack and the flattening is real.
>
> **Two consequences for this entry.** (a) is unaffected — scale-independence is
> a property of any smooth surface whatever its cross-section, so the sweep, the
> absent plateau and `slope_loglog` all stand. **(b) and (c) need restating.**
> On an ellipse with semi-axes ~70 and ~40, curvature genuinely varies around
> the circumference by ~5× — κ ≈ a/b² ≈ 0.044 (R ≈ 23 µm) at the sharp ends
> against κ ≈ b/a² ≈ 0.008 (R ≈ 123 µm) on the flat sides — and it goes
> **negative** wherever the wall folds inward. So "κ2 = −0.0703 is 3–5× too
> curved" compares against a circle the vessel is not, and "there is no good
> `radius_um` to choose" cannot be argued from a single true caliber. The
> defect this entry describes is real; the arithmetic used to size it was
> circular geometry applied to an elliptical vessel.
>
> **A concavity is also real.** The user reported a bean-shaped cross-section at
> a branch end, and it is visible in the **raw** channel as an arc cutting
> inward while staying connected to the ring at both ends. It is not visible in
> the mask, which is [P28](#p28). Note the solidity column above cannot resolve
> it: sealing the ring needs a 12–16 µm closing, which fills any indentation
> shallower than that.

So the vessel radius is of order **40 µm**. The shipped κ2 = −0.0703 (R = 14.2 µm)
is 3–5× too curved to be its circumferential curvature. Consistent with
[P03](#p03)'s bimodal surface-area-vs-depth table (large area at 0–20 µm and
80–105 µm, little between — the top and bottom walls of a vessel that fills the
stack) and with [P21](#p21)'s finding that the lumen is open at z = 30–80 µm.

**(c) Why there is no good `radius_um` to choose.** Past ~15 µm the estimate
overshoots the true caliber and keeps falling (R = 103 µm at `radius_um` = 30),
because at those scales each wall is a broad, gently-curved sheet rather than a
tube — locally almost flat. So the sweep does not contain a correct answer that a
better default would find. Below the vessel scale it reports reconstruction
roughness; above it, sheet flatness. **This cannot be closed by tuning the
slider.**

**What this does not overturn.** [P01](#p01)'s fix. Every measurement behind it is
an internal-consistency check independent of caliber (face-set identity, directed
edge twins 53.7% → 100%, coherence 0.209 → 0.972, ρ(shipped, corrected) = +0.069).
The pre-fix columns were noise; the post-fix columns are self-consistent. Only the
claim that they are *validated against the vessel* is withdrawn — see the caveat
on [P01](#p01)'s acceptance test.

**Where a fix comes from.** Not from `mesh_curvature`, which is doing what it is
asked. It needs a reconstruction whose surface is smooth at the scale curvature is
read at — [P21](#p21) + [I03](#i03) — with the roughness–scale relationship above
usable directly as an objective: a good reconstruction should show a **plateau**
in implied R across `radius_um`, and the shipped one has none. That is a sharper,
cheaper test than the topology/coverage scores [I03](#i03) proposes, and it needs
no ground truth.

**Interim:** do not read `curv_*` as a radius or compare its magnitude across
images or settings. Its *rank within one mesh at one `radius_um`* is the most that
is defensible, and even that is a rank of surface roughness of unknown relation to
the vessel. Any exported `curv_*` should be read together with the `radius_um`
and `subsample_um` that produced it ([P20](#p20)) — neither is currently recorded
next to the column.

**Related:** [P01](#p01) (the fix, and its circular acceptance clause),
[P20](#p20) (`subsample_um` sets the point spacing this is fitted over),
[P21](#p21), [I03](#i03), [P08](#p08)(e) (the same scale-mismatch failure in the
other curvature column).

**History**
- 2026-08-11 — opened. Found while investigating whether `path_curv_xy`'s X-Y
  projection costs anything ([P08](#p08)(e)); the caliber measurements taken to
  check a suspect cross-section are what exposed [P01](#p01)'s circular radius
  target, and the `radius_um` sweep followed from that.
- 2026-08-12 — **the sweep in (a) is now one command, and it reproduces
  exactly.** Curvature needs only the mesh file, so
  `pipeline.batch_mesh --steps curvature` reads neither the image nor the mask:

  ```
  python -m pipeline.batch_mesh --dir 3d_images_folder/images/WT \
      --steps curvature --curv-radius-um 5,10,15,20,30,40 --tag r{radius}
  ```

  Re-run on `V3_..._04step_002` it returns κ2 medians −0.0703, −0.0606, −0.0277,
  −0.0171, −0.0097, −0.0074 — **every value in the table above, to four
  decimals**, from a completely separate code path. The 2026-08-11 measurement is
  reproducible and the entry's evidence stands unchanged.

  Each run writes `vessel_curvature__r<radius>.npz` beside the canonical file
  plus a `vessel_recon__r<radius>.json` recording the radius, the smoothing, the
  mesh it was computed on and that mesh's own parameters — so the plateau test
  this entry names as "a sharper, cheaper test than the topology and coverage
  scores [I03](#i03) proposes" is now routine across all six images rather than a
  bespoke script.

  **This does not advance the entry.** The tool measures the plateau; it does not
  produce a reconstruction that has one. What it removes is the cost of asking.
- 2026-08-12 (later) — **the plateau test this entry asks for is now a number,
  and it is calibrated.** `mesh_score.scale_stability` reports `slope_loglog`,
  the least-squares slope of log|median κ2| against log `radius_um`: **0 means
  the value is a property of the surface, −1 that it is a pure readout of the
  fit window.** It reproduces this entry's table exactly (all six medians to
  four decimals, and the sign-consistency column with them) from an independent
  code path.

  What makes the scale meaningful rather than asserted: on a **sphere of radius
  R**, where both principal curvatures genuinely are 1/R at every scale, it
  recovers κ2 to within 0.2% at every fit radius and returns
  `slope_loglog = +0.02`. The reference vessel mesh returns **−1.191**. So the
  claim in (a) — that the estimate is a property of the window rather than of
  the object — is now quantified against a surface that has the property the
  vessel lacks.

  Two refinements this entry could not previously state. The **per-pair local
  slopes** are −0.21 / −1.94 / −1.67 / −1.39 / −0.94 across 5→10→15→20→30→40 µm,
  so the near-flat first segment is not a plateau in the physical sense: at
  2.25 µm point spacing a 5 µm window holds 6–12 neighbours, exactly the
  rank-deficiency boundary [P20](#p20) names, and the widget-defaults mesh —
  whose 0.97 µm spacing makes that window well-conditioned — shows no flattening
  at all (−1.33). And **no reconstruction tried so far has a plateau**:
  `slope_loglog` is −1.19 to −1.43 across every mesh scored, moving smoothly
  with `subsample_um` and never approaching 0.

  **The entry is unchanged.** `curv_*` still measures the scale it was fitted
  at. What is new is that the defect is now measured per mesh rather than argued
  from one table, so a reconstruction that fixed it could be recognised.

- 2026-08-13 (later the same day) — **the plateau is reached, by the
  [I06](#i06) Poisson reconstruction, on the reference image.** Same mask
  (`min_size=12`, 1,517,396 voxels), same `mesh_score.scale_stability` sweep,
  same six radii:

  | fit radius µm | 5 | 10 | 15 | 20 | 30 | 40 |
  |---|---|---|---|---|---|---|
  | κ2 median (poisson) | −0.0151 | −0.0149 | −0.0149 | −0.0148 | −0.0150 | −0.0155 |
  | κ2 median (alpha 3D, defaults) | −0.1463 | −0.0581 | −0.0298 | −0.0191 | −0.0107 | −0.0082 |

  `slope_loglog` = **+0.009** against alpha 3D's **−1.428** on the identical
  mask, with per-pair local slopes −0.012 / −0.003 / −0.034 / +0.046 / +0.101 —
  flat across an eightfold change of the ruler, where alpha 3D runs −1.33 to
  −0.96 and never flattens. The implied radius sits at **66 µm** and stops
  moving; the previous table had it climbing 14 → 135 µm.

  **What this changes and what it does not.** It changes the *first* clause of
  this entry: on a Poisson mesh `curv_*` is no longer a readout of `radius_um`.
  It does **not** yet make `curv_*` in any exported CSV trustworthy — every CSV
  on disk was measured against an alpha-3D mesh, and this entry stays `open`
  until the six images are rebuilt and re-joined ([I06](#i06) step 6). Nor does
  a plateau make the *value* right: 66 µm implied radius against a 38.1 µm
  measured inscribed radius means κ2 is reading something flatter than the local
  caliber, which is a separate question this test was never able to ask before.
  Stability is a necessary condition, not a sufficient one.

  Also worth recording as a caution: the synthetic elliptic cylinder plateaus at
  `slope_loglog` = **+0.011** and its κ2 recovers the analytic curvature to
  within 10% per angular bin — but under 45% dropout the plateau *survives*
  while the per-bin accuracy degrades to 37%. **A plateau does not imply an
  accurate curvature**, and P24's test cannot tell the two apart. See
  [I06](#i06)'s acceptance table.

- 2026-08-21 — **a second, sharper caution on the same statistic: a plateau does
  not even imply a *plateau*.** The 2026-08-13 line above records that a flat
  `slope_loglog` can accompany an inaccurate curvature. This is worse. Alpha 3D
  built on [P28](#p28)'s local-contrast field returns `slope_loglog` = **+0.024**
  on the reference image — flatter than the Poisson reconstruction this test
  accepted at +0.009 — from a curve whose per-radius local slopes are −1.062,
  −1.496, −1.358, **+1.913, +5.518**. There is no plateau anywhere in it; the
  single log-log fit averages a steep negative half against a steeper positive
  one.

  So the pass/fail this entry defines can be met by a mesh that fails it. Reading
  `slope_loglog` **without** the per-radius table underneath it is not safe, and
  the number should not be quoted alone in this file or anywhere else. The
  per-radius values are already printed by `--steps score`; what is missing is
  any programmatic check that they are monotone or that the fit is good — see
  [I03](#i03).

---

<a id="p25"></a>
## P25 — `fibroblast_count` is never cleared, so a failed image inherits the previous image's counts

| | |
|---|---|
| **Status** | fixed |
| **Severity** | S1 |
| **Opened** | 2026-08-11 |
| **Last verified** | 2026-08-11 |
| **Resolved** | 2026-08-11 (commit `da9f15e`) |
| **Where** | `pipeline/gui/loader.py:186` (the fix), `pipeline/batch.py:227` |
| **Affects** | `n_fibroblasts`, `n_nuclei_total`, and the saved `fibroblast_count.npz` |

Found while tracing [P04](#p04), not by an audit.

The per-image reset in `loader.py` clears `dapi_state`, `cadh_state`,
`clover_state`, `vessel_mesh_cache`, `nuc_cache`, `adhesion_proba` and
`clover_ruby_ratio` — **but not `state["fibroblast_count"]`**, which lives at the
top level of `state` rather than inside one of those dicts.
`batch.py::_free_image_memory` did not clear it either.

The leak has a path to disk. `state["run_fibroblasts"]` chains `_save_fib`
through `call_all`, which fires **even when `_count` bails out early** — and
`_count` bails whenever there is no clover channel (`gui/nuclei.py:175-177`) or
no DAPI labels (`:179-182`), the latter being exactly what happens when
`run_dapi` fails. `batch.py:280-285` catches a failed step and continues to the
next one. So:

> if `run_dapi` fails on image N, image N−1's counts are written into image N's
> `fibroblast_count.npz`, and from there into image N's CSV.

Silent cross-image contamination, with nothing in the log to distinguish it from
a real count — in precisely the unattended five-image run this work was
preparing for.

**Reproduced** headlessly on 2026-08-11 by driving the real
`state["run_fibroblasts"]` with a bailing `_count`: with a stale count in state
the bailed-out step still wrote `27/276` for the next image; with the count
cleared it wrote nothing and the chain continued normally in both cases.

**The fix** is two lines: `state["fibroblast_count"] = None` in the loader's
per-image reset, and `"fibroblast_count"` added to `_free_image_memory`'s key
list. The loader line fixes the GUI and the batch at once, since the batch drives
the same loader; the `batch.py` line covers the free-between-images path.

**Why the fix is placed there rather than in `_save_fib`.** Guarding the save
("don't save if the count is from another image") would need the count to carry
its image identity, which is a larger change for a narrower fix. Clearing state
per image is what every other artifact already does — `fibroblast_count` was
simply missed. This is also why [P04](#p04)'s split nuclei were deliberately put
*inside* `dapi_state` rather than at the top level: that dict is replaced
wholesale by all three of `loader.py`, `batch.py` and `view.py`, so it cannot
develop this failure mode.

**Related:** [P04](#p04) (found while tracing it), [P07](#p07) (the `27/276`
figures above), [P12](#p12) (the same "unattended run, silent wrong number"
class).

---

<a id="p26"></a>
## P26 — The batch builds vessel meshes at widget defaults, so the five new images do not match the reference

| | |
|---|---|
| **Status** | open (decision deferred deliberately — see below) |
| **Severity** | S2 |
| **Opened** | 2026-08-11 |
| **Last verified** | 2026-08-11 |
| **Resolved** | — |
| **Where** | `pipeline/gui/vessel_mesh.py:80-89` (the `controls` defaults), `pipeline/batch.py::_RUN_STEPS` |
| **Affects** | `dist_mesh_um` and all four `curv_*`, **across images** |

Found when the five-image batch produced its first mesh.

`pipeline.batch` drives the GUI's `run_*` steps with **whatever the widgets
default to** — it never sets a mesh parameter. The shipped reference mesh for
`V3_..._04step_002` was hand-tuned, and its settings JSON records values that are
mostly *not* the defaults:

| parameter | reference (`20260806` settings JSON) | widget default = what the batch used | read by `alpha 3D (points)`? |
|---|---|---|---|
| `subsample_um` | **5.0** | **1.0** | yes |
| `alpha3d_um` | **15.0** | **10.0** | yes |
| `min_component_um` | **253.03** | **50.0** | yes |
| `close_um` | 5.0 | 1.0 | **no — inert** |
| `r_smooth_bins` | 12 | 3 | **no — inert** |

> **The last column was added 2026-08-12 and it removes two rows from the
> argument.** Both meshes were built with `alpha 3D (points)`, whose dispatch
> reads only `subsample_um`, `alpha3d_um` and `min_component_um`
> (`vessel_mesh.py::LIVE_PARAMS`, derived from `build_mesh`). `close_um` and
> `r_smooth_bins` belong to the two per-slice methods and did nothing here. The
> settings JSON records every widget regardless, so it *looks* like five
> hand-tuned values; **three of them were.** This is the same class of error the
> entry is about — a recorded number that shaped nothing — one level up.

The effect is exactly what [P20](#p20) predicts, since `subsample_um` sets the
point lattice:

| image | mesh verts | faces | settings |
|---|---|---|---|
| `V3_..._04step_002` (reference) | **8,966** | 15,786 | hand-tuned |
| `V3_..._04step_003` | 205,961 | 77,968 | defaults |
| `V3_..._09step_001` | 270,599 | 90,866 | defaults |
| `V3b_..._04step_001` | 179,468 | 54,922 | defaults |
| `V3b_..._04step_002` | 53,220 | 30,112 | defaults |
| `V3b_..._04step_003` | 57,379 | 24,142 | defaults |

**23–30× more vertices.** [P20](#p20) and [P24](#p24) both establish that
`dist_mesh_um` and `curv_*` are set by the point spacing and the fit radius
rather than by the vessel — so **the geometry columns of the five new images are
not comparable with the reference's**, which is the six-image comparison the
batch exists to enable. `curv_*` is the worse of the two: [P24](#p24)(a) measured
the implied radius tracking `radius_um` at ~3× with no plateau, and `radius_um`
is fixed at 5.0 while the point spacing underneath it moved 5× — so the same
number now means something different per image.

**A second finding, independent of the mismatch:** vertex count varies **5×
between images at identical settings** (53,220 → 270,599), tracking the cadherin
mask size (3.1 → 8.8 MB). A fixed parameter set does not produce a fixed
reconstruction quality across these images, which is direct evidence for
[I03](#i03)'s premise that these knobs need to be *scored* per image rather than
set once. [P19](#p19)'s note that `MIN_SIZE_FLOOR = 12` had never been checked on
the other five images is the same point one stage upstream.

**Not a defect in the batch so much as an unstated assumption.** Nothing records
that the reference mesh was hand-tuned; the values survive only in one settings
JSON, and the two later re-exports of that image
(`20260811_104417`, `20260811_105032`) record **only** `adhesion_postprocess`,
because `batch_measure` does not know the mesh parameters. So the provenance of
the best-available mesh in the project is a single 2026-08-06 file.

> **Correction, 2026-08-12 (user).** This entry, and several others, lean on the
> reference mesh as though it were a known-good target the other five should be
> brought into line with. **It is not.** In the user's words: it is *"somewhat
> hand tuned, but not extensively, since tuning revealed huge parameter
> dependence — keep it since it's better than default widget slider values, but
> don't put too much credence to it."*
>
> Two consequences, and they point in opposite directions from the reading above:
>
> - **"Standardise the five onto the reference" is not the fix.** It was never a
>   validated setting, only a better-than-default one. Re-running the five at the
>   reference's parameters would make the six *comparable* without making any of
>   them *right* — and comparability across six equally-untrustworthy
>   reconstructions is worth less than this entry implied.
> - **The parameter dependence was observed directly at the bench, not just
>   inferred from [P24](#p24)(a)'s sweep.** That is independent, first-hand
>   support for [I03](#i03)'s premise, and it predates the measurements in this
>   file. It is the strongest argument on record that these knobs need to be
>   *scored*, not chosen.
>
> Keep the file — it is still the best reconstruction available and re-running
> the batch is expensive. Stop treating it as the baseline.

**Deliberately not fixed yet (user decision, 2026-08-11).** There is no known-good
parameter set to standardise on: `min_component_um = 253.030303` is visibly tuned
to one image, and [P21](#p21) argues a single global alpha cannot suit vessels
that taper differently. Choosing one blindly could be worse than the defaults.

**What makes this cheap to fix later:** the mesh is built from the cadherin mask,
and `cadherin_labels.npz` is saved for every image. Rebuilding meshes and
curvature therefore needs an image load plus that file — **not** the RF predict,
which is the expensive step. There is no headless script for it today (see
[I04](#i04)); that is the piece of work this entry implies.

> **2026-08-12 — that script exists now: `pipeline/batch_mesh.py`.** It needs no
> image load either, only the OME *header* (see [I04](#i04)'s History for what it
> does and does not settle). It does **not** close this entry: it records
> parameters, it does not choose them, and the five batched meshes are still
> unattributable until they are rebuilt through it. What it changes is that an
> unattributable mesh is now visible — `describe_saved` emits a `WARN` for any
> `vessel_mesh.npz` with no `vessel_recon.json` beside it, so this defect cannot
> recur silently.

**Until then:** do not compare `dist_mesh_um` or any `curv_*` *between* images,
and treat the reference image as being on a different reconstruction from the
other five. Within-image rank remains what [P24](#p24) already allows, and no
more.

**Related:** [P19](#p19), [P20](#p20), [P21](#p21), [P24](#p24), [I03](#i03),
[I04](#i04).

**History**
- 2026-08-11 — opened when the five-image batch produced its first mesh.
- 2026-08-12 — **the reference mesh is demoted from baseline to
  best-available** on the user's correction (see the box above). The vertex-count
  table stands; what it means changes — it is a mismatch between six untuned
  reconstructions, not a drift away from a good one.
- 2026-08-12 — **the vertex-count table above overstates the mismatch, and hides
  a second one.** It quotes *total* vertices. Curvature is fitted on **surface**
  vertices only (`cell_mesh.py:404`, `np.unique(faces)`), and `dist_mesh_um`
  queries a KD-tree over **all** of them ([P09](#p09)(a)):

  | image | total verts | surface verts | **interior** |
  |---|---|---|---|
  | `V3_..._04step_002` (reference) | 8,966 | 7,282 | **18.8%** |
  | `V3_..._04step_003` | 205,961 | 36,642 | **82.2%** |
  | `V3_..._09step_001` | 270,599 | 40,972 | **84.9%** |
  | `V3b_..._04step_001` | 179,468 | 26,298 | **85.4%** |
  | `V3b_..._04step_002` | 53,220 | 14,400 | **72.9%** |
  | `V3b_..._04step_003` | 57,379 | 11,588 | **79.8%** |

  So the curvature mismatch is **5×** (7,282 → 40,972 surface verts), not the
  23–30× this entry claimed — still large, but a third of the headline. The
  claim is corrected downward.

  The `dist_mesh_um` half moves the other way and is the more serious of the
  two. [P09](#p09)(a) measured the rejected-geometry effect as "small — median
  understatement 0.00 µm, p90 0.37 µm", but that was measured on the **one mesh
  where interior points are 18.8%**. On the five batched meshes they are
  **73–85%** of the tree `dist_mesh_um` queries. Measured the same day: the
  share of adhesions whose nearest vertex is an interior point rises 15.1% →
  **39–75%**, and the p90 understatement 0.36 µm → **1.6–5.5 µm**. P09(a)'s
  reassuring numbers describe the reference mesh only. Full table and the fix in
  [P09](#p09)'s 2026-08-12 History entry — that fix needs no mesh rebuild and is
  therefore independent of this entry.
- 2026-08-12 (later) — **[P09](#p09)(a) landed, and it absorbed part of what
  this entry was measuring.** With `dist_mesh_um` now taken to the mesh surface,
  the between-image spread of its median falls from **7.6×** (0.56–4.28 µm) to
  **3.5×** (1.30–4.55 µm). Some of the "not comparable between images" effect
  was the interior-vertex defect tracking point density rather than the
  reconstruction differing. The remaining 3.5× is this entry, and it stands —
  as does everything said here about `curv_*`, which P09(a) does not touch.
- 2026-08-12 (later still) — **the parameter effect is now measured *within one
  image*, which the entry could not do before.** `pipeline.batch_mesh` rebuilt
  `V3_..._04step_002` from its own saved `cadherin_labels.npz` twice — identical
  mask (1,517,396 voxels at `min_size = 12`), identical image, only the mesh
  knobs moved:

  | | reference set (`subsample_um` 5.0, `alpha3d_um` 15.0, `min_component_um` 253.03) | widget defaults (1.0 / 10.0 / 50.0) |
  |---|---|---|
  | total verts | 8,966 | **166,721** |
  | surface verts | 7,282 | **33,537** |
  | faces | 15,786 | 70,176 |
  | `dist_mesh_um` median | 2.298 µm | **1.625 µm** |
  | `curv_k2` median | −0.0701 1/µm | **−0.1632 1/µm** |
  | `path_curv_xy` median | 0.0041 1/µm | 0.0041 1/µm (unchanged) |

  Both CSVs were produced by `batch_measure` over the same adhesions (n = 2,147),
  so this isolates the reconstruction from every other source of variation —
  where the six-image table above confounds the parameter change with a 5×
  spread in mask size. **The surface-vertex ratio is 4.6× and κ2 more than
  doubles**, which is the entry's claim reproduced under control.

  `path_curv_xy` is identical to the last digit, and that is the useful negative:
  the centerline is derived from the cadherin mask, not from the mesh, so it is
  the one geometry column a mesh rebuild cannot move. Its problems are
  [P08](#p08)'s and [P23](#p23)'s alone.

  The reference-set rebuild is **byte-identical** to the `vessel_mesh.npz` the
  GUI itself wrote on 2026-08-11, so none of this is an artifact of the new
  headless path.

- 2026-08-21 — **the same worry does *not* extend to the cadherin threshold, and
  that is now measured.** `cadherin_labels.npz` records no threshold parameters,
  and only one of the six images has a settings JSON that does, so the cadherin
  mask looked like a second instance of this entry: artifacts whose parameters
  are unknown. It is not. The masks reproduce **exactly** — `sauvola_mask(vol, 50,
  0.2, scale) & (vol > median + floor_k·MAD)` with
  `floor_k = min(suggest_floor_k(vol), 10.0)`, which is the loader's own
  auto-seed, gives **Jaccard 1.000 to the voxel on all six images**. `floor_k`
  varies 3.69–5.75 between images precisely because it is measured from each
  volume rather than defaulted, which is why it looked unrecorded.

  So this entry stays about **mesh** parameters. No new entry was opened for the
  cadherin ones; the finding is that there is no defect there. Evidence:
  `docs/scratch/p28/rays.py`, `out/summary.json`.

---

<a id="p27"></a>
## P27 — The centerline is a point set in raster order, not a traced path

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 (no exported column is currently wrong because of it) |
| **Opened** | 2026-08-13 |
| **Last verified** | 2026-08-13 |
| **Resolved** | — |
| **Where** | `pipeline/centerline.py:157` (`np.nonzero(skel)`) — **the file was deleted 2026-08-22**; `lumen_graph.skeleton_graph` returns ordered edges |
| **Affects** | nothing exported today; blocks anything needing arc length or an ordered path |

> ### Fixed 2026-08-22 — the replacement is ordered by construction
>
> `pipeline/lumen_graph.py` returns points **in path order**, edge by edge, with
> `s_um` (arc length from that edge's start), `edge_id`, and a node/edge graph
> whose edges are ordered polylines. Everything this entry lists as blocked is
> now available: arc length exists, there is a parameterisation to resample
> along, and the centerline is drawn in napari as a `Shapes` **path** rather
> than as loose points — the scribble this entry predicts does not happen
> because consecutive points really are consecutive.
>
> The three consumers it worried about were written against the graph directly
> (`measure_core.graph_columns`), so there was never a flat array to migrate.

`vessel_path_curvature_xy` returns `pts_yx_um` straight from `np.nonzero(skel)`,
so the points come out in **raster order** — row by row down the image — not in
the order they occur along the vessel. It is a point *set* that happens to lie on
a curve, and the array order carries no geometric meaning at all.

Measured on the reference image: **40% of consecutive points jump more than
5 µm**, and summing the stored order end to end gives **31,250 µm** across a
structure roughly 300 µm long — two orders of magnitude out, which is what
raster order over a curve that doubles back through many rows produces.

**Nothing exported is wrong because of this.** κ is fitted over KD-tree
neighbourhoods (`query_ball_point`), which are order-independent by construction,
and `dist_branch_um` is a nearest-neighbour distance. Both are geometric, not
sequential. The consumers that would break are the ones that do not exist yet.

**What it blocks**, and why it is worth an entry rather than a comment:

- **There is no arc length**, so nothing can be expressed *per unit length along
  the vessel* — including the axial position variable that [P06](#p06) needs to
  separate depth from position-along-the-vessel, which is currently the main
  unresolved confound in the whole dataset.
- **There is no parameterisation to smooth or resample along**, which is what any
  fix to [P08](#p08)(e) that reports κ with its scale attached would need.
- **Drawing the path as a line is not possible** — it is drawn as points for
  exactly this reason. A `Shapes` path layer would connect raster neighbours and
  produce a scribble.

**A fix is not free on a Y.** Ordering a branched skeleton is not a sort; it is a
graph traversal producing one ordered sequence *per branch* between junction and
endpoint, with the junction shared. That is a different return type from the flat
array every caller reads today, so it is a small piece of design rather than a
one-liner — which is why this is recorded rather than fixed in passing.

**Related:** [P08](#p08) (κ is fitted over neighbourhoods, so it is unaffected —
and that is also why the defect stayed invisible), [P22](#p22) (the Z is per
point and likewise order-free), [P06](#p06) (axial position is the reason to want
arc length).

**History**
- 2026-08-13 — opened, from the module docstring written while adding
  `lumen_center_z`. Found by inspection rather than by a failure: nothing has
  yet asked the centerline for an order.
- 2026-08-13 (later) — **demoted off the critical path.** This was a
  prerequisite for the centerline-frame sweep, which needed arc length and a
  tangent frame. [I06](#i06) replaced that approach and needs neither. Still
  open, still the only route to an axial-position covariate for [P06](#p06).

---

<a id="p28"></a>
## P28 — The cadherin threshold discards real wall; the "gaps" have signal in them

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-13 |
| **Last verified** | 2026-08-21 (five images, both arms, 12 cross-sections; the sixth is unmeasurable) |
| **Resolved** | — |
| **Where** | `pipeline/gui/cadherin.py:164-167` — **the global floor and the `min_size` filter**, not the local threshold (see the 2026-08-21 attribution). `pipeline/denoise.py`'s `sauvola_mask` is exonerated |
| **Affects** | the cadherin mask, and therefore every vessel-geometry column and every reconstruction built on it |
| **Evidence** | `docs/scratch/p28/` (`rays.py`, `klo_sweep.py`, `score_check.py`, `out/`); originally `docs/scratch/centerline_corner/`. Both are gitignored per the project's scratch convention — **the numbers themselves live in the History below**, which is the tracked record |

The whole reconstruction problem has been framed as "the wall has 8–10 µm gaps in
it, so the method must bridge them." Measured on the raw channel, **most of those
gaps are not holes**. They are wall the threshold did not keep.

Ray-cast outward from the lumen centre on the left arm of
`V3_..._04step_002`, comparing directions where the mask found wall against
directions where it did not:

| | mask present | **mask gap** |
|---|---|---|
| peak raw above lumen background | 53–63 counts | **21–24 counts** |
| SNR against lumen noise (~4 counts) | ~14 | **~5.6–6.5** |
| radius at which the peak sits | 59.5–62.5 µm | **59.6–62.4 µm** |

**Same radius, so it is the same wall** — about 40% as bright. The raw ring is
continuous around the entire circumference; the mask is a dashed copy of it, and
the gaps are finely interleaved (mostly a few degrees wide) rather than a few
large arcs.

**Not optical bleed from the neighbouring bright wall.** The signal does not
decay with angular distance from a detected ray (24.2 counts at 1–5° away, 21–23
at 5–10°), and the peak sits at the wall radius rather than smeared inward.

**Why it matters more than it looks.** Every reconstruction method in this
project is judged on how well it bridges gaps, and `alpha3d_um` exists to do
nothing else — [P21](#p21) shows it buys each bridged micron as up to 0.95 µm of
fabricated spur. If the gaps are mostly recoverable signal, the method is being
asked to invent geometry that the data already contains, and the scoring
([I03](#i03)) is grading it on a task it should not have.

~~**The fix is to stop binarising this early.**~~ **Superseded 2026-08-21 —
this paragraph was wrong twice over, and both corrections are measured.** It read:
"`denoise.py` already computes a local mean and standard deviation on the way to
each threshold; the continuous score … keeps dim-but-real wall alive with a small
positive value instead of deleting it … That is the input [I06](#i06) wants.
Raising the threshold's `k` instead would recover the same wall *and* more
speckle, which is the trade a continuous score avoids having to make."

Both halves failed on measurement. **(a)** The local threshold deletes none of
this wall — the floor and `min_size` do — so replacing it changes nothing about
the gaps. **(b)** The continuous score does *not* avoid the speckle trade: the
local z-score does not separate dim wall from lumen noise at all, so the field
buys the same wall *and* the same speckle, which is exactly the trade the
paragraph claimed it escaped. Built, measured, and it degrades every
reconstruction. See the History.

> **Reworded by [P36](#p36), 2026-08-24.** "The local threshold deletes none of
> this wall" is numerically right and the *reading* is wrong: it deletes none of
> it because it is **inert** — `sauvola AND floor` is bit-identical to `floor`
> alone, and `sauvola_k` is a dead knob. "Correctly permissive" and "does
> nothing" are different claims. The attribution to the floor and `min_size`
> below is unaffected; it is now the whole story rather than most of it.

**Related:** [P19](#p19) — this is a stronger charge than that entry makes. P19
argues the size floor is an unjustified constant sitting on a wide plateau; this
says the binarisation upstream of it is losing real wall regardless of where the
floor sits. [P21](#p21), [I03](#i03), [I06](#i06).

~~**Not yet established:**~~ **Both open questions were answered on 2026-08-21;
see the History.** Whether the effect holds on the other images, other arms and
away from the field-of-view edge: **yes**, on five of six (the sixth seals no
cross-section and is unmeasurable). What the false-positive cost of a lower
effective threshold would be: **prohibitive** — the dropped speckle's brightness
distribution overlaps the kept wall's, so no cut separates them.

**Still open:** what *would* separate dim wall from noise, given that neither
raw intensity, the local z-score, nor component size does. The gap wall's one
robust signature so far is **positional**: it sits 0–1 µm from the wall radius
its neighbours report. Nothing in the pipeline currently uses that, and a
neighbourhood-aware recovery is untried.

**A second untried candidate, recorded 2026-08-24: sheetness.** Note what the
three failed discriminators have in common — **all three are intensity or
size**. The wall's remaining distinguishing property is its *shape*: it is a
**sheet**, and lumen speckle is an isotropic **blob**. That is local geometry,
and it is neither of the axes already ruled out.

The repo already computes it. `vessel_poisson.sheet_normals` builds the
structure tensor ``J_ab = G_σ * (∂_a d)(∂_b d)`` and returns
``sheetness = (λ3 − λ2)/λ3`` in [0, 1] — large when the density varies in one
direction only (a sheet), near zero when it varies equally in all (a speck).
Same family as the Frangi/Sato "plate-like" measures; standard technique, not a
new invention.

**It needs no surface, and that is the point** — it is gradients and a
smoothing kernel, nothing else. It *looks* like it needs one only because it
currently lives inside `reconstruct`, on the coarse occupancy grid, after
masking; that is where it happened to be needed (it supplies the normals the
Poisson solve orients, with sheetness as their confidence weight), not a
dependency. **The place to compute it for this purpose is upstream of the
threshold, on the raw intensity**, which also removes any dependence on the
mask being right.

**Two reasons it may fail, both empirical:**

- **Scale.** The structure tensor needs a neighbourhood to see planarity. The
  existing call uses σ = `poisson_smooth_um` = 5 µm on a 2 µm grid; a 1–2 voxel
  speck at 0.325 µm is far below that and may simply vanish into the smoothing
  rather than be classified isotropic. A σ small enough to resolve specks may
  be too small to read the wall as a plane. **There may be no σ where both
  hold.**
- Computed on the *mask* it inherits the mask's errors, which is the argument
  for computing it on the raw channel instead.

**The test is cheap and needs no reconstruction**, and it is the same test this
entry already ran for brightness and for component size: take the gap-wall
peaks recorded here as real wall and the speckle `min_size` drops, compute
sheetness on the raw channel, and ask whether those two populations separate.
Brightness overlapped; size overlapped; this is one number on existing data.

**Priority: low, and deliberately so.** Under the precision-over-recall
asymmetry recorded in this entry's 2026-08-24 History, recovering the dim wall
buys a better-*evidenced* surface rather than a different one. This is the
candidate to reach for **if** the far wall's evidence ever becomes load-bearing
— e.g. for [P32](#p32), where missing wall stops being a confidence question
and lets the level set escape.

> **Attempted 2026-08-24, twice, and both attempts are invalid. The test is
> harder to construct than it looks and the reason is structural — record this
> before re-running it.** `docs/scratch/p28/sheetness.py`.
>
> Populations were built in 3D from this entry's own "the gap peak sits at the
> wall radius": **wall** = voxels a permissive threshold newly admits within
> 2 µm of the strict mask; **speckle** = newly admitted beyond 15 µm. Sheetness
> on the raw channel, σ ∈ {1, 2, 4} µm, with **λ3 (gradient magnitude) as a
> control** — an intensity-like quantity, so if it separates as well as
> sheetness then shape added nothing.
>
> | σ µm | sheetness AUC | sheetness AUC, intensity-matched | λ3 control |
> |---|---|---|---|
> | 1.0 | 0.902 | 0.898 | **1.000** |
> | 2.0 | 0.938 | 0.938 | **1.000** |
> | 4.0 | 0.945 | 0.945 | **1.000** |
>
> **The control reads 1.000, so neither run measured shape.** The populations
> are trivially separable on a gradient-like quantity, and sheetness scores
> *below* the trivial control — it is a noisier proxy for the same thing.
>
> Matching the two populations exactly on raw intensity (40 brightness strata,
> median 558.0 in both) **moved sheetness by 0.004 and left the control at
> 1.000.** That is the informative part: λ3 is a **gradient** measure, so a
> voxel on the steep flank of the bright detected wall carries a large gradient
> *whatever its own brightness*. Brightness was never the confound.
>
> **The confound is the population definition itself.** "Distance from the
> strict mask" is precisely what sets the local gradient, so *every* local
> structure measure separates these two sets by construction. This entry's own
> comparison avoids it — its gap peaks lie in directions where the mask found
> **no** wall, so they are not beside detected wall at all.
>
> **The design that would work**, untried: draw **A** = newly admitted, near the
> `alpha→poisson` *surface* but >5 µm from any strict-mask voxel (in the gaps,
> not the halo), and **B** = newly admitted, far from that surface. Both are
> then dim *and* away from bright wall, so neither sits on a flank. Mild
> circularity — the surface was fit to the strict mask — but it interpolates
> across gaps, so "near surface, far from mask" is this entry's gap wall.
> **Keep λ3, and treat a control near 0.5 as the precondition for believing any
> sheetness number at all.** At 1.000 the test is reporting on something else.

**History**
- 2026-08-13 — opened. Found while checking a claim made earlier the same day
  that the bean-shaped cross-section at a branch end was "missing wall". It is
  missing from the *mask*; the raw channel has it. The first measurement of this
  returned the opposite answer (contrast +1 count in the gaps) because it
  sampled along a moment-fitted ellipse that sits inside the true ring — a
  reminder that a shape assumption in the *measurement* can manufacture the
  result. Ray-casting without assuming a shape reversed it.

- 2026-08-21 — **replicated on five images; the sixth could not be measured.**
  `docs/scratch/p28/rays.py`, generalising the `centerline_corner/` scripts:
  cross-sections at constant x *and* constant y (±4 µm slabs, max-projected),
  the lumen found by closing the mask and labelling the enclosed background (no
  shape assumed anywhere — see the 2026-08-13 line for why that matters), 720
  rays per section, and the wall radius for each direction read from **the
  mask's own longest run of hits in that direction** rather than from a fitted
  curve. The peak is then sought in a ±12 µm window at that radius. 12 sections
  in total across the five measurable images, on both arms and away from the
  field-of-view edge.

  | image | present | **gap** | lumen ctrl | outside ctrl | gap SNR | gap Δr from wall |
  |---|---|---|---|---|---|---|
  | `V3_..._04step_002` | 140.7 | **50.3** | 23.5 | 9.4 | 4.85 | +0.25 µm |
  | `V3_..._04step_003` | 149.4 | **47.8** | 22.7 | 9.9 | 4.60 | +1.00 µm |
  | `V3_..._09step_001` | 133.0 | **45.1** | 23.1 | 17.5 | 4.35 | +1.00 µm |
  | `V3b_..._04step_001` | 144.8 | **46.5** | 19.0 | 13.3 | 4.48 | +0.25 µm |
  | `V3b_..._04step_002` | — | — | — | — | — | *no cross-section sealed* |
  | `V3b_..._04step_003` | 84.4 | **42.0** | — | 14.2 | 4.47 | +0.00 µm |

  Counts are raw above the lumen background, medians per image. **Direction is
  the same in 5/5 measurable images**; magnitudes are not transferable
  ([P03](#p03)'s precedent) and are not offered as such. The gap peak sits
  **0.0–1.0 µm from the wall radius the mask itself reports** — a much tighter
  restatement of "same radius, so it is the same wall" than the original
  measurement made — and ≤2.5% of peaks land on the search window's edge, so
  the window is not clipping the answer.

  **The plan's angular control arm came out empty: 0 rays in every section.**
  The mask ring is closed at every cross-section that seals, so no direction
  qualifies as "no wall anywhere this way". Two *radial* controls were measured
  instead — the same ray one window further in (lumen) and one further out
  (exterior), each required to hold no mask at all. Those are the better null
  for a *peak* statistic anyway: a maximum over ~50 noise samples is biased
  upward, and only a peak taken the same way over the same window says by how
  much. `V3b_..._04step_002` produced **no usable cross-section at all** — its
  mask is the sparsest of the six (524,795 voxels before size filtering) and the
  wall never sealed a lumen at closing radii up to 20 µm. That is a
  measurement that could not be made, not a result pointing the other way.

- 2026-08-21 — **the attribution, and it reverses this entry's proposed fix.**
  Of the gap-wall peaks, the fraction failing each gate:

  | image | n | local threshold only | **global floor only** | both | **`min_size` filter** |
  |---|---|---|---|---|---|
  | `V3_..._04step_002` | 567 | **0.000** | 0.411 | 0.000 | 0.589 |
  | `V3_..._04step_003` | 186 | **0.000** | 0.199 | 0.000 | 0.801 |
  | `V3_..._09step_001` | 79 | **0.000** | 0.405 | 0.000 | 0.595 |
  | `V3b_..._04step_001` | 112 | **0.000** | 0.179 | 0.000 | 0.821 |
  | `V3b_..._04step_003` | 206 | **0.000** | 0.150 | 0.000 | 0.850 |

  **The Sauvola local threshold passes every single gap-wall peak, in every
  image.** It deletes none of this. The mask is built by three gates, not one
  (`gui/cadherin.py:164-167` plus the size slider), and the two that do the
  deleting are the **global intensity floor** (15–41%) and the **`min_size`
  component filter** (59–85%) — the latter being [P19](#p19)'s territory, not
  the binarisation's. So this entry's own recommendation — "the fix is to stop
  binarising this early", replacing the *local threshold* with a continuous
  score — aims at the one gate that is innocent. It stays worth doing for other
  reasons (a score also keeps the *value*, which the Poisson method can use),
  but it is not the fix for the gaps, and [P19](#p19) is now the larger part of
  this defect.

- 2026-08-21 — **the `k_lo` measurement failed, and that is the finding.** The
  plan was to read `k_lo` off this experiment: below the gap wall's local
  z-score, above the no-wall controls'. There is no such value.

  | image | present z | gap z | control z | gap p10 | control p95 |
  |---|---|---|---|---|---|
  | `V3_..._04step_002` | 5.64 | 3.02 | **3.35** | 2.17 | 4.22 |
  | `V3_..._04step_003` | 5.53 | 2.96 | **3.36** | 2.13 | 4.17 |
  | `V3_..._09step_001` | 5.51 | 3.07 | **3.64** | 2.06 | 4.47 |
  | `V3b_..._04step_001` | 5.51 | 2.74 | **3.23** | 1.78 | 4.02 |
  | `V3b_..._04step_003` | 4.75 | 3.23 | **3.30** | 1.74 | 4.08 |

  **The controls score *higher* than the gap wall on the very statistic meant to
  separate them**, in all five images. The reason is that the z-score is a
  contrast measure and the peak is a maximum: in the lumen the local standard
  deviation is small, so a noise peak stands as many local σ above its
  neighbourhood as dim wall does. Raw contrast against the lumen background
  separates cleanly (42–50 vs 9–24 counts); the local z-score does not separate
  at all at `window_xy = 50`.

  Consequence: a rectified z-score field cannot be tuned to keep the dim wall
  without keeping as much non-wall. `docs/scratch/p28/klo_sweep.py` measures the
  cost directly on the reference image, method held at the Poisson defaults:

  | field | voxels | ×mask | occupied cells | lumen core | surface |
  |---|---|---|---|---|---|
  | binary mask | 1,517,396 | 1.00 | 47,102 | 31 | 58,844 v / 114,499 f |
  | `k_lo` = 2.0 (gap p10) | 9,443,248 | 6.22 | 1,366,921 | **0** | **none — refuses to orient** |
  | `k_lo` = 2.5 | 4,230,961 | 2.79 | 1,109,064 | 4 | 65,716 v / 126,770 f |
  | `k_lo` = 3.0 | 2,034,588 | 1.34 | 675,895 | 10 | 59,039 v / 111,480 f |
  | `k_lo` = 4.0 | 639,401 | 0.42 | 145,831 | 18 | 66,269 v / 128,913 f |
  | `k_lo` = 5.0 | 275,645 | 0.18 | 39,503 | 52 | 64,063 v / 124,149 f |

  At the value that would keep 90% of the measured gap wall the field covers
  6.2× the mask, the "wall" swallows the lumen and `lumen_core` refuses to guess
  an orientation — the refusal working exactly as designed. **2.5 is the lowest
  `k_lo` that reconstructs at all**, and it is `batch_mesh`'s default for that
  reason and no other; 3.0 is already back at roughly the gap wall's *median*,
  i.e. throwing half of it away again.

- 2026-08-21 — **the speckle cost.** What `min_size` drops, and how bright it is
  (mean local z-score per component):

  | image | components dropped | voxels dropped | voxels kept | dropped z p50 | kept z p50 |
  |---|---|---|---|---|---|
  | `V3_..._04step_002` | 222,329 | 315,609 | 1,517,396 | 2.93 | 3.91 |
  | `V3_..._04step_003` | 326,181 | 453,315 | 2,056,144 | 2.42 | 3.09 |
  | `V3_..._09step_001` | 259,332 | 483,188 | 2,715,578 | 2.63 | 3.02 |
  | `V3b_..._04step_001` | 353,021 | 488,149 | 1,746,715 | 2.59 | 3.21 |
  | `V3b_..._04step_002` | 67,086 | 101,598 | 423,197 | 3.91 | 4.97 |
  | `V3b_..._04step_003` | 248,533 | 327,269 | 435,300 | 3.28 | 3.52 |

  The distributions **overlap heavily** — dropped speckle sits about 0.3–1.0 z
  below kept wall, less than the spread within either. This is the same message
  as the `k_lo` table from the component side: there is no brightness cut that
  keeps the dim wall and drops the speckle. Note the last two images, where the
  dropped voxels are a quarter to three-quarters of what is kept rather than a
  fifth — `V3b_..._04step_003` drops 327,269 voxels to keep 435,300.

- 2026-08-21 — **the cadherin mask is exactly reproducible, so the [P26](#p26)-
  shaped worry about it does not stand.** `cadherin_labels.npz` stores no
  threshold parameters, and only one of the six images has a settings JSON that
  records them, so the concern was that the masks could not be regenerated. They
  can, exactly: `sauvola_mask(vol, 50, 0.2, scale) & (vol > median + floor_k·MAD)`
  with `floor_k = min(suggest_floor_k(vol), 10.0)` — the loader's auto-seed —
  reproduces the saved `labels > 0` at **Jaccard 1.000, to the voxel, on all
  six images** (1,833,005 / 2,509,459 / 3,198,766 / 2,234,864 / 524,795 /
  762,569 voxels). `floor_k` is per-image (3.69–5.75), which is why it looked
  unrecorded: it is not a default, it is a measurement of each volume. No new
  entry is opened for this — the finding is that there is no defect here.

- 2026-08-24 — **priority downgraded: for reconstruction this is a
  nice-to-have, not a defect to fix.** Recorded on the user's judgement, and it
  reframes what this entry is for.

  Reconstruction is **asymmetric in precision vs recall**, and the two errors
  are not comparable:

  - a **false positive** (lumen speckle) *actively moves the surface* — it
    pulls the level set toward itself, adds bumps, and can break orientation
    outright by destroying the distance-transform contrast `lumen_core` rests
    on;
  - a **false negative** (dropped dim wall) leaves a patch locally
    unconstrained, and Poisson **approximates** rather than interpolates, so it
    bridges smoothly. The surface is still there — it is *inferred* at that
    spot rather than *measured*.

  So getting rid of signal far from the surface matters more than keeping
  signal near it. Recovering the dim wall would not produce a different
  surface, it would produce a **better-evidenced** one. That is worth having
  and is not a prerequisite for anything.

  This is why `min_size` deleting 59–85% of the dim wall is defensible rather
  than scandalous, and it is consistent with every measurement on this entry:
  the more permissive input collapsed Poisson coverage 61.6% → 7.8% and alpha
  3D's 77.0% → 4.2%. **A user-side confirmation, 2026-08-24:** re-running with
  a deliberately permissive cadherin filter (`sauvola_k` 0.43, `floor_k` 2.4)
  produces a *worse* Poisson surface, exactly as this entry predicts — and
  hangs alpha 3D outright, because 94% of the volume's 5 µm bins end up
  occupied ([P35](#p35)).

  **Two caveats keep it from being purely optional.**

  1. **The missing wall is not random — it is depth-systematic.** Attenuation
     means what is dropped is preferentially the *deep* wall
     ([P03](#p03)). "Poisson bridges a missing patch" is fair for scattered
     gaps; "Poisson bridges the entire far side" is a weaker claim. The far
     wall is the least-evidenced part of every reconstruction and it will not
     average out across images.
  2. **Past a certain hole size the level set escapes**, and then it is not a
     confidence question at all — it is [P32](#p32): on images 4 and 5 the
     lumen leaks through the z faces, the graph invents junctions and 44–56% of
     cross-sections clip. That is the one place recovering the wall would be
     *necessary* rather than nice.

- 2026-08-24 — **the floor's 15–41% and [P19](#p19)'s 59–85% are not
  independent: the first partly *causes* the second.** Ordering, in
  `gui/cadherin.py::_run`: the floor is ANDed in **before** `label_and_sizes`,
  while `min_size` is applied **after**, in `_compute_mask`. So the floor does
  not merely delete voxels — it **fragments** a dim wall segment into several
  small components, each of which then falls below `min_size` and is deleted
  whole.

  Two consequences for anyone tuning these:

  - **The two numbers must not be added.** Part of `min_size`'s 59–85% is
    material the floor broke up rather than material that was independently too
    small.
  - **Testing them one at a time understates both.** Lowering `floor_k` should
    return more than its own 15–41%, because it also stops feeding fragments to
    the size filter. The paired experiment is the informative one.

  Mechanism only — not yet measured. It follows from the ordering in the code,
  and the measurement to do is the `floor_k` x `min_size` grid rather than
  either axis alone.

---

<a id="p29"></a>
## P29 — A field-of-view notch fabricates a branch node; `dist_branch_um` is wrong for 18% of rows

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 |
| **Opened** | 2026-08-13 |
| **Last verified** | 2026-08-13 |
| **Resolved** | — |
| **Where** | `pipeline/centerline.py:120` (`binary_fill_holes`), `:152` (the pad), `:181-190` (branch nodes) — **the file was deleted 2026-08-22**; there is no projection and no fill in the replacement |
| **Affects** | `dist_branch_um`, the branch count, and the [P08](#p08)(a) workaround built on both |
| **Evidence** | `docs/scratch/centerline_corner/` |

> ### Superseded 2026-08-22 — the code this entry describes is deleted
>
> `pipeline/centerline.py` is gone. The in-plane centerline was replaced by
> `pipeline/lumen_graph.py`, a 3D graph built on the Poisson lumen rather than
> on a max projection of the cadherin mask — see `docs/DIRECTION.md` §8B for the
> decision and §8C for the measurements. There is no `binary_fill_holes`, no
> `binary_closing`, no `skeletonize` of a projection and no border pad anywhere
> in the replacement, so the mechanism this entry describes cannot occur.
>
> **Left open as a record, not as work.** The entry is not withdrawn — what it
> measured was true of the code that existed — and the replacement has its own
> defects, recorded as [P32](#p32).

> Measured on the replacement, same image: **1 junction**, which is the number
> the vessel has. The old code reported 2 in the saved file and 6 at widget
> defaults (`docs/DIRECTION.md` §8A.3). Across the six the graph gets the
> topology right on four; the two it fails are a different defect with a
> different cause ([P32](#p32)).

`V3_..._04step_002` reports **two** branch nodes. The vessel has one. The saved
`vessel_path_curvature.npz` holds:

```
[[ 51.7 264.7]   <- fabricated, near the top-right corner
 [189.2 182.1]]  <- the real Y junction
```

**Cause: a notch in the footprint that opens onto the image border.**
`binary_fill_holes` fills only *enclosed* background, and a concavity touching
the border is not enclosed — verified, zero enclosed holes survive and both
background components in that window touch the border. So the notch persists and
the medial axis runs around **both** sides of it, producing one arm in +x and one
in +y where a single arm should run into the corner. The notch is ~8–10 µm, the
characteristic cadherin wall-gap scale, so it is dropout ([P28](#p28)) rather
than vessel structure.

**The border padding does not fix it**, and sharpens it. Replicate-padding
extends the footprint *perpendicular to each edge*; near a corner the footprint
touches two edges, so it grows one arm along each. Fork present with and without
the pad, at identical coordinates.

**It is `close_um`, and the default sits just below the threshold that works:**

| `close_um` | branch nodes |
|---|---|
| **8.0 (default)** | 2 — the real Y **plus the fabricated fork** |
| 10.0 | 1 — real Y only; one clean diagonal arm into the corner |
| 12.0 | 1 — same |
| 15.0 | 3 — new spurious nodes appear elsewhere |
| 20.0 | 3 — worse |

**Exported consequence.** `dist_branch_um` is nearest-of-all-branch-nodes, so the
fabricated node captures a whole quadrant: **396 of 2,147 adhesions (18.4%)**
measure to it, understated by a **median 115.5 µm** (53.3 → 160.4 µm corrected).
**43 rows (2.0%)** would be wrongly discarded by [P08](#p08)(a)'s own recommended
`dist_branch_um > 20` filter — they are actually ~160 µm from any real junction.

> Caveat on that reproduction: recomputing the column from the saved nodes
> matches the export to within 4.9 µm but not exactly (median column value
> 107 µm), so measure's lookup differs slightly — probably a downsample-origin
> convention. The which-node-is-nearer classification the 18.4% rests on is
> robust to a few µm.

**The fix is not a bigger `close_um`** — that is a constant tuned on one field,
the [P19](#p19) pattern, and the sweep above shows it breaks elsewhere by 15.
Two structural options:

- **Prune branches by caliber.** The fabricated branch hugs the notch and so runs
  through a shallow part of the distance transform, while the real arm runs along
  its ridge. Dropping child branches whose median inscribed radius is far below
  the parent's kills notch-induced forks generally, and the distance transform is
  already computed a few lines earlier for the pad.
- **Flag branch nodes inside the border-influence zone.** The code already
  computes that distance (p95 of the footprint half-width, 81.7 µm here). The
  fabricated node is 52 µm from the top edge — inside it. The real node is
  143 µm — outside. Wouldn't have to be the fix, but would have caught this.

**This is the third symptom at the field-of-view boundary**, after
[P08](#p08)(c)'s 10× κ inflation in the outer 5 µm and the skeleton forking the
padding addressed. One boundary policy applied consistently across the
centerline, the curvature and the reconstruction would be better than three
patches.

**History**
- 2026-08-13 — opened, reported by the user from the napari view at
  approximately (z 59, y 52, x 266) µm and reproduced exactly at (51.7, 264.7).

---

<a id="p32"></a>
## P32 — The fitted lumen leaks through the z faces on images 4 and 5, and the centerline graph then invents junctions

| | |
|---|---|
| **Status** | open — diagnosed and warned about, not fixed |
| **Severity** | S2 |
| **Opened** | 2026-08-22 |
| **Last verified** | 2026-08-22 (all six images, `alpha 3D solid -> poisson` at one recorded parameter set) |
| **Resolved** | — |
| **Where** | `pipeline/vessel_poisson.py` (the super-level set of χ), surfacing in `pipeline/lumen_graph.py::face_openings` |
| **Affects** | the centerline graph's topology, and every caliber column, on `V3b_..._04step_001` and `V3b_..._04step_002` |
| **Evidence** | `docs/DIRECTION.md` §8C.2 |

The lumen `lumen_graph` builds on is `chi > iso_level`. Where the wall is
missing from the mask, χ's super-level set bleeds through the hole. Some of that
is expected and harmless — cadherin attenuates with depth ([P03](#p03)) so the
deep wall is under-detected everywhere — but on two images the lumen leaves
**mainly** through the z faces rather than through the x/y ones:

| image | lumen voxels on x/y faces | on z faces | junctions found | correct |
|---|---|---|---|---|
| `V3_..._04step_002` | 12,974 | 3,476 | **1** | 1 |
| `V3_..._04step_003` | — | — | **0** | 0 |
| `V3_..._09step_001` | — | — | **1** | 1 |
| `V3b_..._04step_001` | 7,688 | **18,328** | 4 | 1 |
| `V3b_..._04step_002` | 7,955 | **9,748** | 2 | 0 |
| `V3b_..._04step_003` | — | — | **0** | 0 |

**The two failures are exactly the two leaks.** A lumen that has escaped the
stack is not a tube, so its skeleton branches. Both images also clip 44–56% of
their cross-sections, so their caliber is a lower bound as well.

Neither is a surprise on its own terms: [P24](#p24)(b) records image 4 with
1.224% of its mask in the last z slice (150× the reference), which is the
z-containment breach `lumen_core` and the whole orientation step rest on; and
image 5 is the sparsest mask in the set, the one [P28](#p28) could not seal a
single cross-section on.

**`lumen_graph.face_openings` warns whenever z-face lumen exceeds x/y-face
lumen**, and z faces are excluded from the crop test, so the failure is visible
rather than silent. That is a mitigation, not a fix.

**Intersecting the lumen with the z-enclosed core does not work, and was
measured rather than assumed.** On image 4 the core is 2,584 voxels against a
498,469-voxel lumen — 0.5%, a sliver. On image 5 the intersection makes the
topology *worse*: 3 junctions and 13 endpoints against 2 and 4.
`vessel_poisson.lumen_core` is a **seed for orientation**, not a lumen estimate,
and it should not be pressed into service as one.

**Where a fix might come from**, none of it tried:

- A per-image `poisson_smooth_um`. `docs/DIRECTION.md` §8B already flags σ = 5 µm
  as marginal at image 6's caliber, and a σ that lets χ spread through a hole is
  the same lever. **Per-image parameters must be recorded, not chosen by eye**
  ([P26](#p26)) — so this needs the [I03](#i03) selection rule, not a session at
  the sliders.
- A Neumann condition on the z faces as well as x/y, which would state
  "the vessel does not leave through the top or bottom" as a boundary condition
  rather than hoping the wall holds it in.
- [P19](#p19)/[P28](#p28): both images are mask-poor, and the leak is downstream
  of the missing wall.

**Related:** [P24](#p24)(b), [P19](#p19), [P28](#p28), [P26](#p26), [I03](#i03).

---

# Closed

<a id="p30"></a>
## P30 — `build_mesh` silently treats an unknown method string as alpha 3D

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 |
| **Opened** | 2026-08-13 |
| **Last verified** | 2026-08-13 |
| **Resolved** | — |
| **Where** | `pipeline/vessel_mesh.py` — `build_mesh`, the fall-through after the `METHOD_POISSON` branch |
| **Affects** | any caller that passes a `method` string not in `METHODS` |

`build_mesh` dispatches with a chain of `if method == ...` and then **falls
through to alpha 3D** for anything that matched none of them. A typo, a renamed
constant, or a `--from-settings` JSON written by an older version therefore does
not raise — it silently builds a different reconstruction, and
`vessel_recon.json` records the *requested* method beside geometry that method
did not produce. That is the same class of failure as [P26](#p26), where knobs
were recorded as though they had shaped a mesh that never read them.

It is currently unreachable from either shipped caller: `batch_mesh`'s
`--method` is an argparse `choices=METHODS`, and the GUI's is a `ComboBox` over
the same list. The exposure is `--from-settings`, which reads a method string
out of a JSON with no validation, and any future caller.

Found while adding [I06](#i06)'s method, whose branch had to be placed
*before* the fall-through for exactly this reason. Not fixed in that change on
purpose: raising on an unknown method changes what an existing caller gets, and
that belongs in its own change with its own test rather than riding along with a
new reconstruction method.

**Fix:** replace the fall-through with an explicit `if method == METHOD_ALPHA3D`
and raise `ValueError` on anything else, listing `METHODS`.

**Related:** [P26](#p26), [I06](#i06).

---

<a id="p31"></a>
## P31 — `mesh_curvature` has no boundary gate, so curvature is biased near an open edge

**Opened** 2026-08-21 · **Status** open · **Severity** S2 · **Relates to**
[P24](#p24), [I06](#i06)

`cell_mesh.mesh_curvature` fits a quadric to every surface vertex's neighbours
within `radius_um`. It excludes interior and rank-deficient vertices (left NaN),
but it has **no notion of a surface boundary** — a vertex on an open edge gets a
one-sided, roughly half-disc neighbourhood, and the quadric fitted to it is
biased. Those values are returned as ordinary numbers with nothing marking them.

This was harmless while every reconstruction was a closed solid. It stops being
harmless with [I06](#i06)'s Poisson family, which is **open by construction** at
the x/y domain faces — that is the Neumann boundary, and it is correct, because
the vessel genuinely leaves through them.

Measured on `vessel_mesh__alphapoisson.npz` (62,092 verts, 122,799 faces,
1,285 boundary edges = 2.07% of vertices):

| `radius_um` | surface vertices within that distance of an open boundary |
|---|---|
| **5.0** (the shipped default) | **4,387 — 7.07%** |
| 10.0 | 7,776 — 12.53% |
| 15.0 | 11,719 — 18.89% |

So at the shipped radius about **one vertex in fourteen** carries a curvature
fitted to a half-neighbourhood, rising to one in five at the radii
`mesh_score.scale_stability` sweeps. It is not a large fraction and it does not
threaten the plateau result, but these are exactly the vertices at the ends of
the vessel, and nothing downstream can currently tell them apart.

**Fix:** compute distance-to-boundary once per mesh (boundary vertices are those
on edges used by a single face — `mesh_score._edge_use_counts` already finds
them) and either set curvature to NaN within `radius_um` of one, matching the
existing treatment of rank-deficient vertices, or export a
`curv_near_boundary` flag so the join can decide. NaN is the more consistent
choice: the module's stated contract is that untrustworthy vertices are NaN.

**Not** a defect of the Poisson reconstruction — it is a gap in the curvature
code that only an honest open surface exposes.

---

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

**Opened** 2026-08-10 · **Status** open — **the scorer exists and is validated
(`pipeline/mesh_score.py`, 2026-08-12); choosing is what is left** · **Relates
to** [P19](#p19), [P20](#p20), [P21](#p21), [P24](#p24), [I05](#i05)

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

**History**
- 2026-08-10 — opened.
- 2026-08-12 — **the instrument is built; the calibration is not.**
  `pipeline/batch_mesh.py` makes all four knobs this entry names explicit,
  recorded and cheap to vary — `--min-size` (the hidden constant,
  [P19](#p19)), `--subsample-um`, `--alpha3d-um`, `--min-component-um` — with no
  upper bounds, since two of the reference image's values (`subsample_um` 5.0,
  `close_um` 5.0) sit exactly at their slider ceilings and may have been
  truncated by the widget rather than converged. `--suggest-alpha` gives the
  per-image alpha seed headlessly for the first time, and `--dry-run` prints what
  actually resolved.

  Cost, measured on `V3_..._04step_002`: **9.5 s** for mesh + curvature +
  centerline end to end (0.4 s header, 5.0 s mask load, 2.1 s mesh, 1.2 s
  curvature, 0.6 s centerline) at the reference parameters. This entry budgeted
  "10–25 s per point" for an alpha sweep and called minutes acceptable; a coarse
  three-parameter grid is now well inside that.

  **What is still entirely missing is the objective.** Nothing computes genus,
  coverage, parsimony or stability, and nothing scores anything. The tool makes a
  sweep cheap; choosing from it is this entry.
- 2026-08-12 (later) — **the objective now exists and is validated;
  choosing still does not.** `pipeline/mesh_score.py` computes four scores,
  reported as a vector and never combined (`batch_mesh --steps score`, recorded
  into the `scores` block of `vessel_recon.json`):

  | score | question | gameable alone by | caught by |
  |---|---|---|---|
  | `scale_stability` | does κ2 stop depending on the fit radius? | — | — |
  | `topology` | is it shaped like a vessel wall? | collapsing the lumen to a blob | `coverage` |
  | `coverage` | is the surface everywhere the data is? | draping surface over everything | `support` |
  | `support` | is the surface anywhere the data is **not**? | shrinking the surface to nothing | `coverage` |

  **Validated against ground truth, not only against this repo's own history.**
  Reproducing recorded numbers only shows two code paths agree, so the scorer
  was also run on surfaces whose answers are known from theory: genus 0 for a
  sphere, 1 for a torus and unchanged under a 2.8× finer triangulation; and on a
  sphere of radius R it recovers κ2 = 1/R to **0.2% at every fit radius**, giving
  `slope_loglog = +0.02`. That is what calibrates the slope — **0 means the
  number is a property of the surface, −1 means it is a readout of the ruler** —
  and the reference vessel mesh reads **−1.19**. It also reproduces every
  recorded figure: χ = −248, genus 125, the six κ2 medians of [P24](#p24)(a) to
  four decimals with their sign-consistency column, 2.25 µm spacing and 18.8%
  interior vertices ([P09](#p09)).

  **`support` exists because the eye caught what the other three could not.**
  Asked to compare the hand-tuned reference against the widget defaults on the
  renders rather than the numbers, the user reported: *"the reference is more
  solid but includes noise in the mesh, as spiky outliers perpendicular from the
  vessel; the defaults have fewer of those spiky inclusions, but is more
  hole-y."* Both halves were then measurable, and one of them was invisible to
  the scorer as it stood:

  - *hole-y* — `topology` had this and was hiding it. χ's **sign** is the
    distinction: χ ≪ 0 is one shell with tunnels bored through it, χ ≫ 0 is a
    foam of disconnected closed pockets. Cutting the non-manifold edges that
    glue them on, the defaults mesh falls into **1,826 pieces** (1,825 of them
    tiny, median 2 faces) against the reference's **42**. The entry was
    reporting "genus 125" and "genus undefined" and throwing that away.
  - *spiky* — nothing measured it. `coverage` asks "is every wall voxel near the
    mesh?"; a spur is near no wall voxel at all, which is the wrong direction.
    Confirmed by construction: bolting 50 spikes onto a meshed sphere leaves
    `coverage` at **100.0%, median 0.000 — completely unmoved** — while
    `support`'s reach goes 1.12 → 6.04 µm.

  Two dead ends worth recording so they are not re-tried. Surface **vertex**
  distance to data is identically zero on every mesh this project builds (the
  alpha-3D vertices *are* cadherin voxels), so spurs are fabricated *triangles*
  between real points, not stray vertices — `support` therefore samples triangle
  interiors, area-weighted. And **aspect ratio does not separate spikes from
  legitimate bridging**; it ranks the two meshes the opposite way (19.6% sliver
  area on the defaults against 8.5% on the reference). Absolute reach does.

  **What `support` measures turns out to be `alpha3d_um`, near-deterministically**
  — which makes that knob interpretable for the first time. Holding
  `alpha3d_um = 15` and moving `subsample_um` 1 → 5 (a 5× change in the knob
  [P20](#p20) calls dominant) leaves max reach at 14.26 → 13.72 µm and max edge
  at ~29 µm, flat. Across alpha instead, max reach ≈ 0.95 × `alpha3d_um` and max
  edge ≈ 2 × `alpha3d_um`, the latter exactly as a circumradius bound requires.
  So `alpha3d_um` **is** "the longest spur I permit", and [P21](#p21)'s structural
  argument now has a number attached to it.

  **Coverage and support are a genuine trade, not a free win.** Sweeping
  `alpha3d_um` at `subsample_um = 5.0` on `V3_..._04step_002`:

  | `alpha3d_um` | coverage ≤2 µm | support max reach | unsupported area >5 µm | χ |
  |---|---|---|---|---|
  | 5.0 | 26.2% | 4.73 µm | 0.0% | +52 |
  | **7.5** | **47.6%** | 7.06 µm | 4.6% | −263 |
  | 10.0 | 46.7% | 9.27 µm | 13.6% | −358 |
  | **15.0 (reference)** | 42.9% | 13.72 µm | 19.2% | −248 |
  | 20.0 | 39.6% | 19.61 µm | 23.6% | −165 |
  | 30.0 | 30.2% | 29.38 µm | 38.4% | −39 |

  Coverage is **not monotone** — it peaks near 7.5 and falls in both directions —
  so minimising `support` cannot simply drive alpha to zero without coverage
  collapsing. That is the pair working as counterweights. Note also that the
  reference's `alpha3d_um = 15.0` is beaten on *both* of those scores by 7.5
  (coverage 47.6% vs 42.9%, unsupported area 4.6% vs 19.2%), while `topology`
  marginally prefers 15.0 (χ −248 vs −263). **This is one image at one
  `subsample_um` and is not a recommendation** — it is evidence that the scores
  now discriminate between settings that were previously indistinguishable.

  **Knob sensitivity, which [P20](#p20) makes a first-class question.**
  `subsample_um` swept 1 → 5 with the other two knobs pinned: `slope_loglog`
  moves only 0.185, and coverage moves smoothly and monotonically (1.62 → 2.16 µm
  median). Not a knife-edge. The exception is topology — χ runs +2,036 → −248
  without a jump, but **genus is only *defined* over part of that range**, so the
  headline number blinks out between `subsample_um` 2 and 3. Of the four,
  topology is the least suitable for automated selection; χ and the non-manifold
  fraction are its continuous parts, and the non-manifold fraction is itself
  non-monotone (3.7 → 0.4 → 1.6%).

  Cost, measured: `coverage` and `topology` are free (<0.3 s); `support` is ~1 s;
  **`scale_stability` is the entire budget** (6.7 s on the reference's 7,282
  surface vertices, 39.9 s on the defaults' 33,537), because `mesh_curvature`
  loops in Python per vertex with a ball query that grows with the fit radius.
  A five-point build-and-score sweep is ~2 minutes.

  **What is still missing is choosing.** Nothing here picks parameters, and
  nothing weights the four scores — deliberately, since on the reference/defaults
  pair they rank the two meshes *oppositely*, so any weighting would be the thing
  deciding the answer. The selection rule, and the sweep it is applied to,
  remain this entry. Also still absent: the resample-stability score this entry
  lists ("rebuild from a random half of the components"), and any of this run on
  the other five images.

- 2026-08-21 — **the selection rule is now unavoidable, and [P28](#p28) supplied
  a concrete warning about how not to write it.** Two things changed.

  First, the cheaper experiment in front of this entry is spent. [I06](#i06) was
  sequenced so that [P28](#p28) would run before any weighting was written, on
  the argument that [P28](#p28) might reinterpret the coverage column and dissolve
  the trade. It ran, with the method held fixed and the incumbent verified
  byte-identical, and **the comparison came out unchanged**. Nothing cheaper than
  the rule is left.

  Second, and more useful: **a worthless mesh scored best on the acceptance
  criterion.** Alpha 3D on [P28](#p28)'s field returned `slope_loglog` = **+0.024**
  — numerically better than the accepted Poisson reconstruction's **+0.009** —
  while covering **4.2%** of mask voxels at 2 µm, reaching **101.70 µm** from any
  data, and presenting χ = +171, which this scorer correctly describes as "a foam
  of roughly 85 disconnected closed pockets around the main shell, not one wall".

  The mechanism matters for the rule's design. `slope_loglog` is a single
  least-squares fit through log κ against log r, and that mesh's *local* slopes
  are −1.062, −1.496, −1.358, **+1.913, +5.518**, with `IQR/|median|` reaching
  31.3 at r = 20 µm. The fit averages the negative and positive halves to near
  zero. A flat curve and a violently non-monotonic one produce the same summary
  number, and **the summary cannot tell them apart**.

  Two consequences for this entry:

  - **`scale_stability` should expose monotonicity or fit quality alongside the
    slope**, so a non-plateau cannot masquerade as a plateau. The local slopes
    are already computed and printed; nothing reads them programmatically.
  - **The four-vector doctrine is now empirically load-bearing, not just
    principled.** The only reason this was caught is that coverage, support and
    topology all disagreed with the slope. Any weighting that could let one
    column outvote the other three would have accepted this mesh.

---

<a id="i04"></a>
## I04 — `batch.py` has no way to re-run one step

**Opened** 2026-08-11 · **Status** open · **Relates to** [P04](#p04), [P12](#p12)

`pipeline.batch` runs `_RUN_STEPS` all or nothing: there is no `--steps` flag.
Only two things can be redone headlessly afterwards — adhesion instances
(`batch_postprocess`) and measurement (`batch_measure`). Everything else (DAPI,
the split nuclei, the cadherin mask, the vessel mesh, both curvatures, the RF
probability map) can only be regenerated by re-running the whole image, which is
the expensive part.

**Why it is worth recording rather than just fixing.** It is the reason
[P04](#p04), [P12](#p12) and [P25](#p25) had to land *before* the five-image
batch instead of alongside it: a run is effectively one-shot, so any change to a
**saved artifact** must be in place first or the compute is wasted. That
sequencing constraint is invisible in the code and cost real ordering work.

A `--steps` flag is not hard — the step list is already a list of
`(label, state_key)` pairs, and each `run_*` is independently chainable through
`call_all`. What needs thought is dependencies: `run_fibroblasts` needs
`run_dapi`'s output in `state`, so re-running one step alone means loading the
prerequisite artifacts first, which is `load_*`-shaped work rather than a bare
filter over the list. The GUI's "Load all saved features" already does exactly
that and could be reused.

**Until then:** `--skip-existing` is the only safety mechanism, and it keys on
`adhesion_instances.npz` — the *last* artifact written — so it is all-or-nothing
per image too.

**History**
- 2026-08-11 — opened.
- 2026-08-12 — **the vessel-geometry third of it is done, by a separate runner
  rather than a `--steps` flag on `batch.py`.** `pipeline/batch_mesh.py` re-runs
  the mesh, the surface curvature and the in-plane centerline for any image, in
  any combination, from the saved `cadherin_labels.npz`.

  It sidesteps the dependency problem this entry flags ("re-running one step
  alone means loading the prerequisite artifacts first, which is `load_*`-shaped
  work") by not going through `state` at all. Three things made that possible,
  and all three were checked rather than assumed:

  - the isotropic **scale** is `min(voxel_zyx_um)` from the OME header
    (`io.to_isotropic`), so no pixels are decoded — ~0.3 s per image;
  - `round(Z · pz / target)` reproduces the saved mask's Z **exactly on all six
    images**, including the 0.9 µm `09step` file, so the same header read
    doubles as a check that the mask belongs to the image;
  - `cell_mesh.mesh_curvature` needs only the mesh, so a curvature-only run
    touches neither the header nor the mask.

  Drift between the GUI and the headless path is prevented structurally rather
  than by discipline: the three closing strategies were moved into
  `vessel_mesh.build_mesh`, and `gui/vessel_mesh.py`'s worker is now a single
  call to it — the same arrangement `saved_check.py` was created under. Verified
  end to end: rebuilding the reference image at its recorded parameters produces
  an npz **byte-identical** to the one the GUI wrote, and a headless smoke test
  of the real napari tab confirms its mesh equals a direct `build_mesh` call.

  **Still open:** DAPI, the split nuclei, the cadherin mask itself and the RF
  probability map remain all-or-nothing under `pipeline.batch`. Those are the
  steps that need the image pixels, which is exactly why they are the ones left.


> **Done 2026-08-24.** `batch.py --steps` takes `all`, a comma-separated list,
> or a `first:last` range over the stable short names `dapi, fibroblasts,
> cadherin, mesh, curvature, graph, clover, rf, postprocess`. Order is always
> the pipeline's own, never the order typed, so a later step consuming an
> earlier one's output cannot be run first by typo; unknown names and empty
> ranges are argparse errors, not silent no-ops.
>
> **The caveat is in the help text and matters:** a step re-run alone still
> reads the tab defaults, so it reproduces an earlier run only if those have not
> moved ([P26](#p26)) — which [P37](#p37) shows is not something to assume.
---

<a id="i05"></a>
## I05 — Is the alpha shape the right reconstruction method at all?

**Opened** 2026-08-12 · **Status** open · **Relates to** [P19](#p19),
[P20](#p20), [P21](#p21), [P24](#p24), [I03](#i03)

The requirement, stated plainly: **a smooth, continuous 3D surface from which to
extract geometry measures, built from a noisy cadherin channel.** Alpha 3D is the
default because it "bridges gaps in 3D and handles branch + curvature natively"
(`vessel_mesh.py` docstring) — it has never been compared against anything. This
entry opens that question and records what has actually been measured, now that
[I03](#i03)'s scorer makes comparison possible. **Nothing here changes the
default.**

### The structural objection

An alpha shape **interpolates**. Its vertices *are* the cadherin voxels (a
`voxel_subsample` of them), so it is constrained to pass through its input and
therefore **cannot smooth**. Every surviving piece of noise becomes a vertex, and
`alpha3d_um` decides only how far the surface may reach out to meet it. This is
not a tuning difficulty; it is what the method is.

Measured consequences on `V3_..._04step_002` at the reference parameters:
8.3% of faces sit on non-manifold edges, the surface falls into 42 pieces when
those are cut, genus is 125, and **19.2% of its area is more than 5 µm from any
cadherin voxel**.

### What was measured

**1. The alpha knob is now mechanistically understood** ([I03](#i03),
[P21](#p21)). Max spur ≈ 0.95 × `alpha3d_um`, max edge ≈ 2 × `alpha3d_um`
(exactly, as a circumradius bound requires), independent of `subsample_um`. So
the knob is "the longest spur out of the wall I permit". Coverage is
**non-monotone** in it, peaking near 7.5 µm.

**2. Sensitivity to the mask is a cliff, not a slope — and it is smaller than
the sensitivity to alpha.** This partly answers the worry that prompted this
entry. Across `min_size` 3 → 500, a 167× range, coverage at 2 µm moves only
40.3% → 36.7% and max spur 23.4 → 10.7 µm ([P19](#p19) History for the table).
Across `alpha3d_um` 5 → 30 coverage moves 26.2% → 47.6% → 30.2%. **The mesh knob
is the more sensitive of the two**, not the mask. The one exception is
`min_size = 0`, where alpha 3D bridges 226,093 speckle components into a single
smooth envelope — χ = −3, coverage **2.5%**, spurs to 78 µm.

**3. The two other in-repo methods are not practical to compare at full
resolution.** `polar (per-slice)` was run on the same mask with the same scores
and **did not complete in over an hour**, against **2.1 s** for alpha 3D. The
per-slice fills run over ~1,024 slices in two orientations, then `close_volume`
and marching cubes operate on the full 340 M-voxel volume. Whatever their merits,
they cannot participate in the kind of sweep [I03](#i03) needs without being
reworked to run on a decimated volume. **Not yet measured**, therefore: whether
either produces a better surface.

**4. A density-field alternative was prototyped, and it does not yet win.** The
approach that answers the structural objection is to **approximate** rather than
interpolate: block-reduce the mask to a local occupancy density on a coarser
grid, Gaussian-smooth it, take an isosurface (`docs/scratch/mesh_score/`,
`m3`–`m6`). Isolated speckle contributes little density and falls below the
level; real wall is dense over a neighbourhood and survives.

Three of its four predicted advantages hold, measured against the reference
alpha-3D mesh on the same mask and the same fixed scoring mask:

| | alpha 3D (reference) | density (1.3 µm grid, σ 1.5, level p25) |
|---|---|---|
| closed edges / non-manifold | 98.4% / 8.3% of faces | **100.0% / 0.0%** |
| pieces after cutting non-manifold edges | 42 | 1–17 |
| genus | 125 | **25** |
| support max reach | 13.72 µm | **3.00 µm** |
| unsupported area > 5 µm | 19.2% | **0.0%** |
| coverage ≤ 2 µm | 42.9% | 40.4% |
| coverage ≤ 5 µm | **99.1%** | 73.0% |
| build cost | 2.1 s | ~2 s |

- **It is manifold in every configuration tried** — 100.0% closed edges and zero
  non-manifold edges at every sigma, level and `min_component_um`. Alpha 3D is
  never manifold.
- **It has essentially no spurs** — 0.0% of area beyond 5 µm, against 19.2%.
- **Its genus/coverage trade is monotone in one knob** (the iso-level): genus
  4 → 25 → 56 → 141 → 203 as coverage ≤5 µm rises 54% → 73% → 88% → 93% → 95%.
  Alpha 3D's is **not** monotone (χ = +52, −263, −358, −248, −165, −39 across
  alpha 5 → 30), which matters a great deal for automated calibration.
- **It degrades rather than collapses without any component filtering.** At
  `min_size = 0` it goes from genus 10 to 94 with coverage roughly unchanged,
  where alpha 3D collapses to a blob at 2.5% coverage.

**5. Smoothness is controllable by the smoothing scale — and trades against
coverage.** This is the measurement that bears directly on the requirement, and
it took three attempts because the first configuration tried was the worst one.

| surface | `slope_loglog` (5–40 µm) | 5→10 µm local slope | coverage ≤2 µm |
|---|---|---|---|
| alpha 3D (reference) | −1.191 | −0.214 | 42.9% |
| density, σ = 1.5 µm (1.3 µm grid) | **−2.252** | −1.570 | 40.4% |
| density, σ = 3 µm (2.6 µm grid) | −1.622 | −0.253 | 7.8% |
| density, σ = 5 µm (2.6 µm grid) | **−1.099** | +0.151 | 3.6% |
| density, σ = 8 µm (2.6 µm grid) | not measured | **+0.226** | not measured |

σ = 1.5 µm on a 1.3 µm grid is barely one cell of smoothing, and it reads −2.25
— far *worse* than alpha 3D. Taken alone that looked like a refutation. It is
not: smoothness improves **monotonically** with σ (−2.25 → −1.62 → −1.10), and by
σ = 5 µm the density surface matches alpha 3D (−1.10 against −1.19) while also
being manifold and spur-free. At σ = 5 and 8 the 5→10 µm local slope goes
*positive*, meaning |κ2| stops falling as the ruler grows — the beginnings of the
plateau [P24](#p24) says no reconstruction here has ever shown.

**But coverage collapses exactly where smoothness arrives**: 42.9% → 7.8% → 3.6%
at 2 µm as σ goes 3 → 5. So **no configuration tested is simultaneously smooth
and covering.** (Caveat on the table: the σ = 1.5 row is on the fine grid at
level p25, the others on the coarse grid at p50, so grid and level are not held
fixed across rows. The monotone trend in σ is clear; the exact values are not a
controlled comparison.)

**And it still loses the coverage tail.** At the level where its genus is 5×
better, 27% of wall voxels are more than 5 µm away, against 0.9% for alpha 3D. It
reconstructs part of the wall tightly and misses the rest.

### One confound to keep in mind

`coverage` structurally favours interpolating methods: alpha 3D's vertices *are*
data points, so it is being graded partly on a tautology. `support` penalises the
same property. That the density surface nearly matches on coverage while beating
support by 4.6× is therefore more meaningful than the raw numbers suggest — and
is also why neither score should be read alone here.

### Verdict, on the evidence so far

**Alpha 3D remains the best method actually available, but the case against it is
now empirical as well as structural, and the density field is the more promising
direction.** Alpha 3D is the only one of the three in-repo methods cheap enough
to calibrate, its dominant knob is now interpretable, and it is less
mask-sensitive than feared. Against that, the density field is better on
manifoldness (always, categorically), spurs (0.0% against 19.2%), monotonicity of
its knobs, robustness to masking, and — at σ ≥ 5 µm — smoothness. What it has not
done is deliver smoothness and coverage **at the same time**.

The honest summary is that the two methods fail in different places and neither
is adequate: alpha 3D covers the wall but is a non-manifold sponge with 13.7 µm
spurs; the density field is a clean manifold surface that either misses much of
the wall or is no smoother. **Nothing here justifies switching the default, and
nothing here justifies treating alpha 3D as settled.**

### What would settle it, in order

1. **Find whether a smooth-and-covering density configuration exists.** This is
   now the whole question. The two properties were measured on different grids
   and levels; the trade has never been mapped as a surface over (grid, σ, level,
   `min_component_um`). Since the smoothness/σ relation is monotone and each
   point costs seconds to build, this is a small sweep — but it needs
   [I03](#i03)'s selection rule first, or there is no way to say which point on
   the trade is best.
2. **Localise the coverage misses.** 27% of wall voxels more than 5 µm away looks
   like whole regions dropped by the level or by `min_component_um`, not uniform
   error. If it is a few dropped regions, it is fixable and the verdict changes.
3. **Try a genuine surface *fit* rather than an isosurface.** Both methods tested
   are still thresholding operations. Poisson reconstruction or a moving
   least-squares fit would approximate in the sense the requirement wants;
   neither is in the dependency set today, and adding one must be installed by
   hand in the Windows env.
4. **Make the per-slice methods measurable** by running them on a decimated
   volume, so all three in-repo options can be compared at all.

### A practical obstacle this exposed

`scale_stability` is too slow to be the inner loop of a method search on
marching-cubes meshes. `mesh_curvature` loops per surface vertex and its ball
query grows with the fit radius, so a 50,000-vertex surface at `radius_um = 40`
takes tens of minutes — several attempts here had to be abandoned. The 5→10 µm
**local** slope is comparable across meshes, costs seconds, and separated the
configurations just as well (alpha 3D −0.214, density σ 1.5 −1.570, σ 5 +0.151).
Any method sweep should use the cheap local slope and reserve the full sweep for
finalists.

None of this should start before [I03](#i03)'s selection rule exists — otherwise
there is again no way to say which of two candidate surfaces is better, which is
the trap this whole line of work was opened to escape.

**Update, 2026-08-21 — post-hoc filtering of the alpha complex is closed off.**
The structural objection above says an alpha shape cannot smooth because its
vertices *are* its input. The obvious escape is to delete the bad tets instead:
`alpha3d_um` is one scalar doing two jobs (max edge = 2α is the widest gap it
bridges, max spur ≈ 0.95α is the longest thing it lets stick out), and since the
wall is a thin **sheet**, a gap-bridging tet should lie *in* the sheet plane while
a spur should run *along* the local sheet normal. Filtering on "extent along the
sheet normal" would split the one knob into two, letting α rise for coverage while
spurs were forbidden separately — which is exactly the bind [P21](#p21) describes.

Measured (`docs/scratch/mesh_score/m7_tet_perp.py`, reference image, reference
parameters, sheet plane from `vessel_poisson.sheet_normals` at the Poisson
defaults). The harness reproduces this entry's reference alpha-3D row exactly
(local slope −0.214, χ −248, genus 125, 8.3% non-manifold, 42 pieces, coverage
42.9%/99.1%, support 13.72 µm), so the result is about the idea, not the
measurement. **It fails three ways:**

- **Perpendicular reach is not bimodal** — p25 1.82, p50 3.01, p90 7.85 µm, one
  smooth mode with a monotone tail. No valley, so no principled threshold.
- **The "narrow within" half runs backwards.** In-plane extent *rises* with
  perpendicular reach (median 6.91 → 10.82 µm). High-reach tets are fat in every
  direction, not thin spikes; there is no long-and-narrow population to select.
- **The sweep degrades the mesh monotonically.** Genus 125 → 194, non-manifold
  8.3% → 16.3%, pieces 42 → 165 as the threshold tightens ∞ → 6 µm, with
  `slope_loglog` moving −0.214 → −0.373 *away* from the plateau. Support max
  barely moves (13.72 → 12.61, then flat) — the filter does not even remove the
  spurs it was aimed at. The best sheetness gate is the one that filters least.

Face count goes *up* as tets are removed (15,786 → 20,728), so the tets being
deleted are **interior** to the wall slab: removing one opens a passage rather
than shaving a spur. **Any scheme that prunes the alpha complex after the fact
will hit this** — the tets it wants gone are load bearing, not peripheral. Alpha
3D's noise is the tail of one continuum, not a separable geometric class.

**Correction, same day.** The first version of this update explained that as
*enclosed voids* — sealed bubbles whose walls became new boundary faces. That was
asserted without measurement and it is **wrong**.
`docs/scratch/mesh_score/m8_fill.py` rasterized the alpha solid exactly (via
`Delaunay.find_simplex`) and ran `binary_fill_holes`: **0 cells added, in 0
enclosed components.** There are no voids. The genus is **tunnels through the
slab**, which connect exterior to exterior and which no fill operation can close.
The refutation above survives — it rests on the reach distribution, not on the
mechanism — but the mechanism was mis-stated.

**And measuring that turned up something larger, which is now the actionable
result of the session.** Arms 1 and 2 of that script are *the same alpha solid*,
differing only in how the surface is read off it — `alpha_shape_faces` extracting
tet boundary faces, versus rasterizing at 1 µm and running marching cubes:

| | non-manifold faces | closed edges | coverage ≤2 µm | support max |
|---|---|---|---|---|
| `alpha_shape_faces` (today) | 8.3% | 98.4% | 42.9% | 13.72 µm |
| same solid, rasterized | **0.5%** | **99.9%** | **87.1%** | 14.03 µm |

Non-manifoldness essentially gone, coverage **doubled**, support unchanged, for
~1 s — so a real part of what this file records as "alpha 3D is a non-manifold
sponge" belongs to the face-extraction step rather than to the alpha shape, which
bears on [P01](#p01), whose entire subject is `alpha_shape_faces` mangling the
surface.

**"Free win" was too strong, and is withdrawn.** Re-measured undecimated at a
grid matched to Poisson's (`RESULTS.md` §12), the rasterized surface reads
`slope_loglog` **−1.167** against the baseline's −0.214 — far worse, and worse
than a pure ruler readout. Marching cubes on a binary volume is a voxel
staircase, which is exactly the roughness a 5 µm quadric fit cannot see past. The
topology and coverage gains stand; they are **bought with smoothness**, which is
the one axis [I06](#i06)'s acceptance test is about.

**A second pruning criterion was tried and failed the same way.** The wall is a
sheet of measured thickness ~2.5 µm, so a tet belonging to it must be thin in at
least one direction, and that direction should be the sheet normal. The physical
picture is **confirmed** — median |thin axis · sheet normal| decays monotonically
0.897 → 0.477 as minimum width grows — but deleting the thick tets gives 42.0%
non-manifold faces and 999 pieces at a 2.5 µm cut. **Two independent criteria,
one with sound reasoning and one without, produce the same failure.** The alpha
complex is not prunable through `alpha_shape_faces` at all: removing any tet
turns its former interior into new boundary, and which tets you choose is not the
issue. Considered closed.

**Update, 2026-08-13 — a fourth candidate, and it is now the plan of record.**
See [I06](#i06). Two structural findings from this session close off part of the
option space above:

- **A slice can contain two lumens.** Downstream of the Y a plane cuts both
  daughter branches, so any explicit per-station contour method has to assign
  contours to branches, track merges, and decide what the surface does through
  the junction. That is bookkeeping an implicit representation never performs —
  each branch contributes a field, the fields are unioned, and marching cubes
  runs once.
- **The cross-section is neither circular nor reliably convex** ([P24](#p24)(b)):
  aspect 1.64–2.15, with a real inward fold visible in the raw channel at a
  branch end. Any radial parameterization — `r(s, θ)`, per-slice polar, the
  straightened-frame sweep this entry's item 1 proposes — assumes a single-valued
  radius per direction and would report the vessel as an ellipse everywhere,
  erasing the concavity **with no warning**. That rules out the whole family, not
  just tunes it.

Item 3 above ("try a genuine surface *fit*") is what [I06](#i06) is, and its
dependency objection turns out not to bind: the grid form is a divergence and a
Laplace solve, both already in scipy.

---

<a id="i06"></a>
## I06 — Oriented-point implicit reconstruction (the plan of record)

**Opened** 2026-08-13 · **Status** in progress · **Supersedes** the
centerline-frame sweep sketched in [I05](#i05) · **Relates to** [P01](#p01),
[P21](#p21), [P24](#p24), [P28](#p28), [I03](#i03)

### What the goal actually requires

Working backwards from "curvature at any point on the lumen surface, related to
adhesion force", the surface must be: **smooth at the scale curvature is read at**
(else [P24](#p24)), **correctly oriented** (else the sign of curvature has no
meaning), **covering without fabricating** (the [I03](#i03) pair), and
**joinable** to each adhesion. The centerline is not on that list. It is a means,
and the question is only which means is cheapest.

### The method

1. **Dense oriented points.** The wall is a thin sheet, so locally the cadherin
   lies in a plane and the perpendicular to it is the surface normal. Computed
   from the gradient of a lightly smoothed field. The saved per-component normals
   in `cadherin_labels.npz` are unit-length and non-degenerate but there are only
   **3,764** of them over ~190,000 µm² of wall — one per ~7 µm, too sparse for a
   wall with 8–10 µm gaps and for curvature read at 5–15 µm.

2. **Orient them from the centerline.** Local data gives the *line* a normal lies
   along, never which end is the lumen. Pointing each normal at the nearest
   centerline point resolves it — and **this is the correct fix to
   [P01](#p01)'s part B**, which used away-from-global-centroid and was wrong on
   44.4% of vertices precisely because a global centroid is the wrong reference
   on a branched tube. Measured on the saved normals: agreement with
   away-from-centroid is **42.1%**, i.e. a coin flip. Orienting *inward* is also
   what selects the **luminal** face, rather than the alpha shape's closed solid
   whose area counts both wall sheets ([P03](#p03)).

3. **Solve for the surface.** Take the indicator that is 1 inside the lumen and 0
   outside; the normals are measurements of the direction in which it changes
   fastest. Recovering it is one Poisson solve over the volume, then isosurface
   at the half level. It **approximates instead of interpolating**, which is
   [I05](#i05)'s structural objection to alpha 3D; noise averages out rather than
   becoming vertices, gaps are spanned by a smooth solution rather than by a
   fabricated triangle, and topology, branching and concavity are handled without
   being told about them.

### Why this de-risks every open centerline defect

Step 2 is the centerline's only job here, and it is forgiving: the tube is
85–180 µm across, so the axis need only be *somewhere inside it*. A normal flips
only if the axis lands on the wrong side of the wall. **[P29](#p29)'s fabricated
branch node, [P27](#p27)'s raster ordering, [P08](#p08)(c)'s border flare and the
~1.3 µm quantization of `lumen_center_z`'s Z are all far below that threshold.**
They remain defects; they stop blocking the surface.

That is the whole argument for this over the sweep, which needed the centerline
ordered, accurate, branch-aware and edge-clean — i.e. needed every one of those
fixed first.

### Acceptance

**[P24](#p24)(a)'s plateau test, which is already implemented and already
calibrated.** `mesh_score.scale_stability`'s `slope_loglog` is 0 when curvature
is a property of the surface and −1 when it is a readout of the fit window; a
sphere returns +0.02, the reference vessel mesh **−1.19**, and no reconstruction
in this project has ever plateaued. Pass/fail for the whole geometry effort, and
it needs **no force data** — which is what makes this runnable now under the
2026-08-13 scope decision.

Scored alongside on the full [I03](#i03) four-vector, against alpha 3D and the
density prototype, **on the same mask** so the comparison is honest.

### Order of work

1. **Method first, on the existing binary mask** — deliberately, despite
   [P28](#p28) being the bigger defect. Two reasons: `mesh_score`'s coverage and
   support are defined over *retained mask voxels*, so building on the same mask
   is the only way the comparison against alpha 3D and every number already
   recorded in this file is honest; and the acceptance test is the **plateau**,
   which is about the smoothness of the surface rather than about gaps, so it is
   answerable without [P28](#p28). Changing the input and the method together
   would leave any improvement unattributable.
2. Dense normals, oriented from the centerline; check what fraction of
   orientations are unstable under perturbation of the axis.
3. Grid Poisson solve → isosurface. Knobs: grid resolution, regularisation
   weight, iso-level.
4. Score, and accept or reject on the plateau.
4a. **Then [P28](#p28), as a controlled second experiment** with the method held
   fixed: swap the binary mask for the continuous local-contrast field. If it
   helps, it helps by shrinking the gaps, and it should show up as coverage
   rising and the surface having less to invent. Note this moves the yardstick —
   see Open questions.
5. Sign convention for `curv_*` relative to the lumen — only meaningful once
   normals are oriented — and export `curv_*` with the `radius_um` and
   reconstruction parameters that produced it, which [P24](#p24) notes is
   currently missing.
6. Add the method to `vessel_mesh.build_mesh` so the GUI and `batch_mesh` share
   one path (the no-drift doctrine), record it in `vessel_recon.json`, re-join
   and re-export all six.

**Added alongside alpha 3D, not replacing it** — nothing existing breaks, and the
comparison stays runnable.

### Deliberately deferred

[P27](#p27)'s ordering and arc length (only needed for an axial-position
covariate), the remaining centerline edge work beyond [P29](#p29), the
bean question, and everything in [P02](#p02)/[P05](#p05).

### Open questions

Whether one normal per *component* would do after all (cheap to try, expected too
sparse); whether the Poisson solve wants screening/regularisation or plain; and
whether [P28](#p28)'s continuous input changes what `mesh_score.coverage` is
measured against, since coverage is currently defined over *retained mask
voxels* and that yardstick moves if the mask does.

**History**
- 2026-08-13 — opened as the plan of record, after the session recorded in the
  Status-at-a-glance queue notes.

- 2026-08-21 — **the plateau reproduced independently, and feeding Poisson a
  gap-bridged input was tried and does not pay for itself.**
  `docs/scratch/mesh_score/m9_minwidth_poisson.py`, five arms on the reference
  image, undecimated, raster grid matched to the Poisson grid:

  | arm | `slope_loglog` | χ | closed / non-mf | pieces | cov ≤2 / ≤5 µm | sup max |
  |---|---|---|---|---|---|---|
  | alpha 3D | −0.214 | −248 | 98.4% / 8.3% | 42 | 42.9% / 99.1% | **13.72** |
  | rasterized alpha | −1.167 | −2,658 | 99.5% / 1.5% | 98 | **80.8%** / 98.4% | 14.16 |
  | **Poisson on binary** | **−0.018** | −5 | 99.0% / **0.0%** | **1** | 61.4% / 87.5% | 26.70 |
  | Poisson on rasterized alpha | **−0.024** | −3 | 99.3% / **0.0%** | **1** | 61.6% / **91.3%** | 46.47 |

  The Poisson row reproduces this entry's recorded numbers (χ −5 vs −4, coverage
  61.4%/87.5% vs 61.6%/87.9%) from a separate harness, which is the first
  independent confirmation of the plateau result.

  Feeding it the alpha-bridged solid **works** — coverage ≤5 µm 87.9% → 91.3%
  and ≤10 µm 93.3% → **98.2%**, and `lumen_core` reported **5** z-enclosed
  fragments instead of **31**, i.e. alpha really did close wall gaps.

  **Corrected 2026-08-21, later the same day.** This entry first read the
  `support` **max** (26.70 → 46.47 µm) as a near-doubling of fabricated surface
  and called the trade not worth taking. **That was wrong.** A max is one worst
  point. Area-weighted, the two surfaces are **19.6% vs 19.9%** of area beyond
  5 µm, with medians 2.27 vs 2.43 µm and p90 6.90 vs 7.01 — indistinguishable
  until the top ~5% tail. `support` is an area-weighted score and must be read
  that way; quoting its max alone mischaracterises a surface.

  Two further corrections in the same direction. `support` is measured against
  the **retained** mask, which [P19](#p19)/[P28](#p28) measure as missing 59–85%
  of the dim wall — so bridged surface is scored unsupported wherever the mask is
  wrong, which is precisely where this method operates. And the surface was
  inspected in napari and read as visually clean, corroborating that the max is
  an outlier.

  **Corrected verdict: it dominates Poisson-on-mask** — same plateau
  (`slope_loglog` −0.052 on the full sweep), same topology (1 piece, 0.0%
  non-manifold), materially better coverage, statistically indistinguishable
  fabrication outside a small tail. Still not made the default: that is a
  selection decision [I03](#i03) leaves to whoever runs the sweep, and it has
  been measured on **one image**.

- 2026-08-21 (late) — **`alpha 3D solid -> poisson` is the direction of record.**
  The surface was inspected in napari and read as clean; the numbers agree
  (plateau `slope_loglog` −0.052 full sweep, 1 piece, 0.0% non-manifold,
  coverage ≤5 µm 91.3% / ≤10 µm 98.2%, fabrication indistinguishable from
  Poisson-on-mask by area). The method is **built and shipping as a choice**;
  alpha 3D is still the default and changing that is a separate decision.

  **Curvature off this surface is well defined** — the question was asked
  directly and checked against the code and the mesh:

  - **0 non-manifold edges**, so the normal and hence the sign of κ exists
    everywhere. Alpha 3D never managed this (8.3% of faces).
  - **The normals are the real ones.** `build_mesh` returns ∇χ normals in
    `extras["normals"]` and they are passed explicitly into `mesh_curvature`.
    This is load-bearing: on an open surface `Σ p0·(p1×p2)/6` is not
    origin-independent, so the enclosed-volume branch would put a coin-flip
    global sign on κ2 — the very quantity the plateau test reads.
  - **Scale-stable.** k2 median across `radius_um` 5→40 µm: −0.0153, −0.0151,
    −0.0148, −0.0143, −0.0138, −0.0140; implied |R| ≈ 65–72 µm throughout.
    62,039 / 62,092 finite, 88.0% negative rising to 95.3% at r = 40 — a
    consistent sign, concave toward the lumen as a tube should be.
  - **One gap, now [P31](#p31):** no boundary gate, so 7.07% of vertices at
    `radius_um` 5.0 carry a quadric fitted to a half-neighbourhood.

  ### The three refinements this still needs

  **(a) Small bumps from bad tets.** Visible in the viewer: a spur tet makes a
  lump of solid → a lump of occupancy → the level set follows it. The lever
  already exists and is `alpha_raster_um`. **Beware the interaction**, which cost
  a wasted run: `occupancy_grid` uses `k = round(grid_um / raster)` and
  `h = k · raster`, so the raster grid *changes the realised solve grid*:

  | `alpha_raster_um` | k | realised `h` | effect |
  |---|---|---|---|
  | 2.0 | 1 | 2.00 µm | binary solid — **no spur filtering at all** |
  | 0.8 | 2 | **1.60 µm** | finer solve, ~2x the vertices, much slower |
  | **0.65** | 3 | **1.95 µm** | 27-level fractional occupancy at ~the same cost |

  **0.65 is the setting to use**: it is exactly 2 × the 0.325 µm voxel, and
  `h = 1.95 µm` matches what Poisson-on-mask already runs at, so runtime is
  unchanged while a thin spur fills only a fraction of its cell and drops below
  the iso level. Not yet swept or scored.

  **(b) Runtime — and it is not the reconstruction.** Measured on the
  `alphapoisson` batch run: mesh build **3.2 s**, curvature at one radius
  **9.4 s**, `--steps score` over six radii **96 s**. The Poisson solve is ~3% of
  the total. The rest is `cell_mesh.mesh_curvature`, which loops **in Python over
  every vertex** doing a KD-tree ball query plus a least-squares quadric fit, and
  which `scale_stability` calls six times — the same cost [I05](#i05)'s "a
  practical obstacle" section records. Two real fixes, in order of value:
  **vectorise `mesh_curvature`** (batch the ball queries and the lstsq — same
  maths, no per-vertex Python, and it speeds up *every* method), or **read
  curvature off the implicit function** (the surface is a level set of χ, so the
  principal curvatures come from the Hessian of χ projected onto the tangent
  plane — pure array ops, no tree, no fit). The second is a fast path for the
  shipped `radius_um` only: it gives curvature at the smoothing scale, and
  [P24](#p24)'s plateau test needs the multi-scale probe. Meanwhile `--steps
  mesh` alone is ~3–5 s and is all a visual check needs.

  **(c) Open3D — worth having, for two things, and it cannot replace the
  module.** `open3d.geometry.TriangleMesh.create_from_point_cloud_poisson` is
  Kazhdan's *screened* Poisson and is genuinely more refined than
  `vessel_poisson`'s uniform-grid finite-difference Neumann solve: an **adaptive
  octree** (ours is uniform 2 µm, which leaves the ~2.5 µm wall about *one cell*
  thick), a **B-spline basis** rather than finite differences, a **screening
  term** that pulls the level set through the samples — the lever for the
  coverage shortfall — and a **per-vertex density** output, which is the standard
  tool for trimming exactly the low-confidence fabricated tail measured above.

  **But it cannot supply the orientation.** Kazhdan's method *requires* oriented
  normals as input, and Open3D's own
  `orient_normals_consistent_tangent_plane` is generic propagation that has no
  concept of a lumen — this entry already measured a naive orientation rule at
  **42.1%**, a coin flip, which would silently flip the sign of κ2.
  `lumen_core` + `orient_inward` encode [P24](#p24)'s z-enclosure measurement,
  and `lumen_core` *refuses* rather than guessing when the vessel is not
  z-contained. **The seam is therefore: keep our normal estimation and
  orientation, swap the solver**, handing Open3D an already-oriented point cloud.

  It is **not** the fix for (b) — the solve is not the bottleneck.

  **Library only, never a second viewer.** napari stays the single viewer; any
  Open3D mesh returns as verts/faces through the same `add_surface` path as every
  other method. `open3d.visualization` is not to be imported.

  ### Order of work from here

  1. Sweep `alpha_raster_um` at **0.65** against the 2.0 run — cheap, no new
     dependency, targets (a) directly.
  2. [P31](#p31)'s boundary gate.
  3. Vectorise `mesh_curvature` — (b), and it pays back across every method.
  4. Open3D screened Poisson as a sixth method, fed our oriented normals,
     install-gated.
  5. Re-measure on the other five images before any of this becomes the default;
     **everything above is one image** ([P26](#p26)).

  **Unexplained and blocking any quotation of the number:** `support` max for
  Poisson-on-binary reads 26.70 µm here where this entry records **9.87 µm** for
  the same arm — and 26.70 is exactly the value recorded for *alpha 3D* in the
  same table. One of the two is a transcription error and it must be resolved
  before either is cited.

- 2026-08-21 — **`alpha 3D solid -> poisson` added as a fifth method**, on
  request, with the other four untouched. CLAUDE.md's four-place rule was
  followed in one change: `METHODS` + `LIVE_PARAMS` + `ALL_PARAMS` in
  `vessel_mesh`, `MESH_DEFAULTS` + the `explicit` dict + `--alpha-raster-um` in
  `batch_mesh`, plus the GUI slider. New public function
  `vessel_mesh.rasterize_alpha_solid` returns the alpha complex's *interior* as
  a volume (exact, via `Delaunay.find_simplex`); the one new knob
  `alpha_raster_um` defaults to the Poisson grid, where `occupancy_grid` reduces
  to the identity and the input Poisson sees is the binary solid.

  Verified through the shipping path with `batch_mesh --tag alphapoisson`: it
  reproduces `RESULTS.md` §12's scratch numbers exactly (62,092 verts, 122,799
  faces, coverage 61.6% / 91.3%, support max 46.47 µm), and alpha 3D, Poisson
  and [P30](#p30)'s unknown-method fall-through all still build bit-identical
  meshes. Full-sweep `slope_loglog` is **−0.052** (the 5→10 µm local slope
  −0.024 matches §12).

  **Added, not adopted.** Alpha 3D remains the default. The coverage-for-support
  trade this method makes is not adjudicated in code — `mesh_score` returns a
  vector and never a total, and [I03](#i03) puts the selection rule with whoever
  runs the sweep.

- 2026-08-21 — **[P28](#p28)'s controlled experiment ran, and the continuous
  field made the reconstruction worse on every column.** Method held fixed at
  the Poisson defaults (grid 2.0 / smooth 5.0 / iso 50.0), reference image,
  scored against the pinned **binary** point set in every arm so the numbers stay
  comparable to everything above.

  | arm | `slope_loglog` | χ | coverage ≤2 µm / ≤5 µm | support max |
  |---|---|---|---|---|
  | alpha 3D on binary | −1.428 | +2,870 | 77.0% / 99.9% | 26.70 µm |
  | Poisson on binary | **+0.009** | −4 | 61.6% / 87.9% | 9.87 µm |
  | **Poisson on the score field** (`k_lo` 2.5) | **−0.399** | −11 | **7.8% / 25.3%** | 42.10 µm |
  | alpha 3D on `score > k_lo` | +0.024 † | +171 | **4.2% / 17.1%** | 101.70 µm |

  **The fourth arm degrades the incumbent too, and harder.** Alpha 3D
  interpolates, so it takes `score > k_lo` as a looser mask. That mask holds 2.8x
  the voxels — but after `subsample_um = 1.0` it holds **2,352,268 points against
  the binary mask's 166,721, i.e. 14.1x**. The two ratios disagreeing is
  diagnostic on its own: recovered *wall* would land on the 1 µm lattice cells
  the mask already occupies and would barely move the point count, so a 14x
  increase says the added voxels are **diffusely scattered through the volume**.
  That is the same conclusion [P28](#p28)'s z-score tables reach from the other
  direction, arrived at independently.

  Everything downstream follows. Coverage falls from 77.0%/99.9% to 4.2%/17.1%,
  support reach grows from 26.70 µm to **101.70 µm**, χ goes to +171 — which the
  scorer correctly refuses to call a genus, reporting instead "a foam of roughly
  85 disconnected closed pockets around the main shell, not one wall". Only
  46,018 of the 2.35M vertices end up on the surface (98.0% interior, against
  [P09](#p09)(a)'s already-bad 80%). And it cost **600 s in `build_mesh` against
  ~12 s for the binary arm** — Qhull tetrahedralising 14x the points,
  single-threaded.

  † **`slope_loglog` = +0.024 here is not a plateau and must not be read as
  one.** The per-radius local slopes are −1.062, −1.496, −1.358, **+1.913,
  +5.518** — violently non-monotonic, with `IQR/|median|` reaching 31.3 at
  r = 20 µm. A single log-log fit through that averages the negative and positive
  halves to something near zero. The binary Poisson arm's +0.009 comes from a
  curve that is actually flat; this +0.024 comes from a curve that is not a
  curve. See [I03](#i03) — this is the clearest example yet of a scorer column
  being gamed by a summary statistic, and it was caught only because the other
  three columns disagreed with it.

  **The coverage gap does not close; it collapses.** And it is not an artifact of
  the pinned yardstick: measured against the *field's own* 4.2M nonzero voxels
  the same surface covers 5.5% at ≤2 µm and 15.9% at ≤5 µm, i.e. worse still.
  The plateau — the one thing this method was accepted on — is **lost**
  (+0.009 → −0.399), χ moves the wrong way, and the surface reaches 42 µm from
  any data where the binary arm reached 9.9.

  **Why, and it is [P28](#p28)'s own finding.** The local z-score does not
  separate dim wall from lumen noise (see that entry's tables), so at the lowest
  `k_lo` the method tolerates the field carries 2.8× the mask's voxels, most of
  the addition being noise rather than wall. `sheet_normals` then fits sheets to
  noise, `lumen_core` fragments, and the level set follows the fabricated
  density. This is the failure mode the plumbing was built to expose, and it
  exposed it.

  **So the question this entry deferred is answered in the negative.** The open
  item above — "whether [P28](#p28)'s continuous input changes what
  `mesh_score.coverage` is measured against" — is moot: the input is not usable
  as it stands, so the yardstick does not have to move. Alpha 3D still wins both
  coverage columns and Poisson still owns the plateau, exactly as on 2026-08-13,
  and **the head-to-head is unchanged by this experiment** — which is a result,
  not a null. What would change it is a discriminator that actually separates
  dim wall from noise; the z-score at `window_xy = 50` is not one.

- 2026-08-13 (later the same day) — **built, measured, and accepted on its own
  acceptance test.** `pipeline/vessel_poisson.py`, wired into
  `vessel_mesh.build_mesh` as `METHOD_POISSON` (`"poisson (oriented points)"`),
  into `batch_mesh` (three flags, `MESH_DEFAULTS`, `LIVE_PARAMS`) and into the
  napari tab (three sliders). Alpha 3D is untouched and still the default;
  rebuilding the reference mesh at its recorded parameters is still
  byte-identical to `vessel_mesh.npz`.

  **The headline.** On the reference image, same mask, `slope_loglog` =
  **+0.009** against alpha 3D's **−1.428** — see [P24](#p24)'s History for the
  full table. The plateau test passes. That is the pass/fail this entry set in
  advance, and nothing in this project had ever plateaued before.

  **The four-vector, same mask, same scorer** (alpha 3D column is at *widget
  defaults*, rebuilt for this comparison rather than quoted; the hand-tuned
  reference mesh's figures are in [I03](#i03)):

  | | poisson | alpha 3D (defaults) |
  |---|---|---|
  | `slope_loglog` | **+0.009** | −1.428 |
  | local slopes 5→40 µm | −0.03 … +0.10 | −1.33 … −0.96 |
  | χ | **−4** (≈3 tunnels) | +2,870 (a foam of ~1,435 pockets) |
  | non-manifold faces | **0 (0.0%)** | 13,122 (18.7%) |
  | pieces after the non-manifold cut | **2** (largest 99.3%) | 1,826 (largest 77.1%) |
  | interior vertices | **1.2%** | 79.9% |
  | coverage ≤2 µm | 61.6% | **77.0%** |
  | coverage ≤5 µm | 87.9% | **99.9%** |
  | support max reach | 26.70 µm | **9.87 µm** |
  | unsupported area >5 µm | 19.6% | **16.5%** |
  | longest triangle edge | **3.36 µm** | 19.90 µm |
  | mesh build | **3.1 s** | 13.4 s |
  | rebuild determinism | byte-identical | byte-identical |

  **Read as a vector, which is the whole point of [I03](#i03)'s refusal to
  total it.** Poisson wins stability and topology outright and is 4× cheaper;
  alpha 3D wins both coverage columns and support reach. That is not a
  contradiction, it is the approximate/interpolate trade stated numerically: a
  surface forced through every voxel is by construction near every voxel
  (coverage), and a surface that averages instead leaves the noisiest voxels
  behind. **Which one is right depends on whether those voxels are wall or
  noise, and this entry does not settle that** — [P28](#p28) is the experiment
  that would, with the method now held fixed.

  **Synthetic ground truth** (`docs/scratch/vessel_poisson/synthetic_check.py`,
  full log beside it): an elliptic cylinder with the measured 152×88 µm
  cross-section, 45% of the shell deleted in correlated patches with the arc gap
  capped at 12 µm, 50,000 speckle blobs, plus a Y variant and an analytic
  reference mesh at matched vertex spacing.

  | | criterion | result |
  |---|---|---|
  | (a) | plateau, \|slope_loglog\| < 0.15 | **PASS** +0.011 |
  | (b) | median\|κ1\| < 0.2 median\|κ2\|; κ2 sign ≥ 95%; K<0 at the Y crotch | **PASS** (ratio 0.013; 100%; 50.8% K<0) |
  | (c) | per-bin κ2 within 15%; end-to-side ratio within 25%; estimator drift < 0.10 | **PASS** (worst 9.9%; 12.1% off; 0.068) |
  | (d) | ≥99% of normals agree with the analytic outward normal | **PASS** 100.0% |
  | (e) | p95 distance < 0.5 µm both directions, max < 3 µm | **FAIL** 0.752 / 0.825 µm |
  | (e2) | \|signed half-cell bias\| < 0.15 µm | **PASS** −0.003 µm |
  | (f) | 0 non-manifold faces; 1 piece; 2 boundary loops, both on domain faces | **PASS** (Y variant: 3 loops) |
  | (g) | (c) and (e) survive dropout; drift < 1 µm; no interior loops | **FAIL** — see below (but (a) survives: +0.0099) |
  | (h) | two runs byte-identical | **PASS** |
  | (i) | alpha 3D on the same phantom fails (a), (c), (f) | **PASS** (−1.196; 860%; 93.6% non-manifold) |

  **(e) is a resolution statement, not a defect.** The residual tracks
  `poisson_grid_um` and not `poisson_smooth_um`: at a fixed h = 1.95 µm, σ of
  3 / 5 / 8 µm all give p95 0.72–1.03 µm, while dropping the grid moves it
  directly — p95 **0.752 µm at grid 2.0**, **0.507 at 1.0**, **0.400 at 0.65**.
  The criterion as written never mentioned the grid; at the default it is asking
  for sub-⅓-cell accuracy. The half-cell coordinate correction it was really
  aimed at is (e2), which passes at −0.003 µm where dropping the `−h/2` term
  would have put it at +0.98 µm.

  **(g) is a real limitation and is recorded as one.** Under 45% dropout the
  plateau *survives* — `slope_loglog` = **+0.0099**, against +0.0106 on the
  clean shell — but the surface wobbles ±1.5–3 µm at the
  gaps: p95 mesh→analytic **3.19 µm**, worst per-bin κ error **37%**, drift from
  the clean mesh p95 **3.12 µm**, and **2 interior boundary loops** where the
  harmonic solution did not span a hole. Finer grids do not help (tested at 1.0
  and 0.65 µm), which is what identifies it as genuine surface error rather than
  discretisation. **So a plateau does not certify an accurate curvature** — the
  two criteria come apart exactly here, and [P24](#p24) carries that caution.

  **Deviations from the plan of record, all deliberate:**

  1. `build_mesh` gained an optional `extras=` out-dict rather than the planned
     "call `diagnose()` separately". The per-vertex normals are needed for
     *correctness* — `mesh_curvature` cannot recover the global sign on a
     surface this open — and re-deriving them would mean building the whole
     reconstruction twice, which a sweep cannot afford. The return is still a
     4-tuple, so `saved_check`, `batch_mesh` and the GUI are unchanged.
     `diagnose()` still exists.
  2. Every percentile in the module is **occupancy-weighted**. The unweighted
     form counts a cell holding one speckle voxel as loudly as a cell packed
     with wall; on the phantom, speckle cells outnumber wall cells ~6:1 and the
     unweighted median put the wall threshold at a speckle density, after which
     `lumen_core` found nothing at all.
  3. The wall level in `lumen_core` is the occupancy-weighted **p10** ("the band
     holding 90% of the wall mass"), not the planned p50. At p50 the band is
     only the top half of the density ridge and the lumen leaks into the
     exterior through it — measured on the *clean* phantom, which has no gaps,
     so the leak was the threshold and not the data.
  4. `lumen_core` requires wall **above and below each cell in its own z
     column**, in addition to the planned connectivity + z-face rule. Pure
     connectivity leaks: thresholding widens a hole by ~σ on each side, so a
     [P28](#p28)-scale 12 µm gap becomes a ~22 µm channel and the lumen merges
     with the exterior through it. Per-column enclosure is the direct statement
     of what [P24](#p24) measured (0.01% of the mask in the first z slice, zero
     in the last) and a lateral leak cannot defeat it. Cost: the lumen comes
     back as fragments — 31 on the reference image — whose union is used.
  5. The global-sign guard compares median to median, not median to p99. See
     [P01](#p01)'s History; the p99 form silently inverted the real image.

  **Measured and rejected, recorded so it is not rediscovered.** Normalising the
  splatted field by the local sample density (Kazhdan's density compensation)
  was tried, because χ's jump across the wall is proportional to local wall
  density and a *global* iso-level therefore sits off-surface where the wall is
  sparse. It helps the dropout geometry (p95 2.44 → 1.82 µm, max 8.04 → 3.75 µm,
  bias +0.29 → +0.09 µm) and costs a little on the clean shell (p95 0.752 →
  0.868 µm) — but it does **not** rescue (c) under dropout (61% vs 66%), and it
  needs a clamp constant on the boost. Not adopted: a hidden constant that does
  not buy the criterion it was aimed at is exactly what this project has paid
  for twice.

  **Cost**, reference image: 3.1 s to build (against alpha 3D's 13.4 s at
  defaults, 2.3 s at the reference's `subsample_um=5`), so it is usable as the
  inner loop of an [I03](#i03) sweep, which was the other thing this had to be.

  **Status stays `in progress`,** because steps 4a (the [P28](#p28) experiment)
  and 6 (re-join and re-export all six images) are not done, and because the
  four-vector does not yet say which method to adopt. What is settled is that
  the plateau is reachable.

### Still open, in order

1. **[P28](#p28) as the controlled second experiment** — the method is now the
   fixed half. Step 2 of the algorithm needs no change for it: `occ` is already
   a continuous field.
2. **Boundary-aware genus in `mesh_score`** — `genus = (2 − χ − b)/2` for `b`
   boundary loops. Without it every Poisson mesh reports `genus: None` (it is an
   open sheet), so the head-to-head against alpha 3D's genus 125 cannot be made
   on that axis. A loop counter exists in `synthetic_check.py` and should move
   into `mesh_score` rather than be written twice.
3. Whether the 66 µm implied radius is right, against the 38.1 µm measured
   inscribed radius. New question, only askable now that κ2 stopped moving.
4. Re-join and re-export all six images; sign convention for `curv_*`; record
   `radius_um` and the reconstruction parameters beside every `curv_*`.

---

<a id="p33"></a>
## P33 — `sauvola_mask` cast to float32 on the host, doubling the PCIe transfer

| | |
|---|---|
| **Status** | **fixed** 2026-08-24 |
| **Severity** | S3 (performance only — output is bit-identical) |
| **Opened** | 2026-08-24 |
| **Last verified** | 2026-08-24 (reference image, 324x1024x1024 uint16) |
| **Resolved** | 2026-08-24, `pipeline/denoise.py::sauvola_mask` |
| **Where** | `pipeline/denoise.py` — `sauvola_mask`'s `v = vol.astype(np.float32)` |
| **Affects** | every cadherin mask build; `_run_cadherin_step` forces Sauvola, so this was the interactive *and* batch path |

`gpu.local_mean_std` is deliberately written to transfer the volume in its
**native dtype** and cast to float32 **on the device** — its docstring: "the
volume crosses the PCIe bus only once — transferred as uint16 (half the bytes
of float32) ... host<->device copies are not [cheap]".

`niblack_mask` honours that (`_local_zscore` passes `vol` through untouched).
**`sauvola_mask` cast on the host first**, so it shipped float32 and defeated
the one optimisation that function exists for. Both call `_mean_std` exactly
once, so this was the entire difference between them.

It compounded rather than merely doubling: on the reference stack native is
**679 MB** against **1359 MB** as float32, and the larger figure exceeds the
pinned-memory pool, so the copy silently degraded to a **synchronous** transfer
(and to a CPU fallback beyond that). The cupy warning names the exact byte
count — `1358954496 bytes` — which is the float32 size to the byte, and it is
what identified the cause.

| | time | foreground |
|---|---|---|
| before (host float32 cast) | **8.79 s** | 339,738,527 |
| after (native dtype), same process | 5.40 s | 339,738,527 |
| after (native dtype), **clean warm run** | **2.89 s** | — |

**Bit-identical** — verified by array comparison, not inferred. Nothing
downstream changes and no artifact is invalidated.

**Read the third row, not the second.** The first two were measured
back-to-back in one process, so the "after" run started cold *and* in the
pinned-memory pool the "before" run had just exhausted — it is an overestimate
of the fixed path. Measured cleanly the fix is far larger than the 1.63x that
pairing suggests.

**The decisive check is parity with Niblack**, which is what the diagnosis
predicts: both call `_mean_std` exactly once and differ in nothing else costly,
so if the transfer was the whole story they must land together. Warm, same
volume, same window:

| | time |
|---|---|
| `niblack_mask` | 2.92 s |
| `sauvola_mask` | **2.89 s** |

They do. (A first call in a fresh process costs ~4.7 s regardless of method —
cupy context and memory-pool init, not a property of either threshold. Benchmark
either one cold and the comparison is meaningless.)

`v` was only used for `.max()`/`.min()` and the `_mean_std` call — the final
comparison already used `vol` — so removing the cast was free.

**`bernsen_mask` has the same host-side cast (`denoise.py:227`) and is NOT
fixed**, deliberately: there the cast does real work, because
`midgray = 0.5 * (hi + lo)` would overflow in uint16. The correct fix there is
to filter in native dtype and cast afterwards. Bernsen is also inherently far
costlier than the other two at the same window — min/max are rank filters,
where `uniform_filter` is an O(1) running sum — which is expected behaviour
rather than a defect.

**Related:** [P28](#p28) (the Sauvola threshold is the innocent gate),
[I07](#i07) (the same host/device transfer discipline on the clover path).

---

<a id="p34"></a>
## P34 — The Cadherin tab computed per-component centroids eagerly, for a consumer that is not on the main path

| | |
|---|---|
| **Status** | **fixed** 2026-08-24 |
| **Severity** | S3 (performance only — no output changes) |
| **Opened** | 2026-08-24 |
| **Last verified** | 2026-08-24 (reference image) |
| **Resolved** | 2026-08-24 — `gui/cadherin.py::_ensure_centroids`, `gui/cell_mesh.py` |
| **Where** | `gui/cadherin.py::_run` called `cell_mesh.label_centroids_normals(labels, scale, min_voxels=1)` on every run |
| **Affects** | every Cadherin run, interactive and batch; the reason a permissive `floor_k` was unusable |

`label_centroids_normals` is a **Python loop over every component**, doing
`labels[sl]`, `np.argwhere(sub == lid)`, a mean, a covariance and a 3x3
`np.linalg.eigh` each. The tab called it with **`min_voxels=1`**, so it ran on
every 1-voxel speck — and the function's own docstring says PCA on clouds that
small is unstable, so most of the output was unusable by construction.

Cost against a **constant ~3.7 s** for `label_and_sizes` on the same volume:

| `floor_k` | components | `label_and_sizes` | `label_centroids_normals` |
|---|---|---|---|
| 5.0 | 170,236 | 3.70 s | 4.98 s |
| 4.0 | 375,467 | 3.65 s | 9.92 s |
| 3.0 | 1,335,402 | 3.74 s | **34.45 s** |

**Labelling itself is not what scales — this loop is.** Component count grows
superlinearly as the floor drops (7.8x for 3.4x the voxels), because what a
lower floor admits is speckle. At `floor_k` 2.0 the process reached **2654% CPU
and 13 GB RSS** and was killed; the memory is `ndi.find_objects` plus millions
of small `argwhere` allocations.

**Nothing on the main path reads the result.** Traced every consumer:

- the **vessel mesh consumes the binary mask**, not the labels
  (`gui/vessel_mesh.py::_cadh_mask_array` reads the mask *layer*);
- `sizes` alone drives the size filter, the histogram and
  `auto_min_size_otsu`;
- `centroids`/`normals` have **exactly one consumer in the codebase** —
  `gui/cell_mesh.py`'s optional "Refine with cadherin" button, which is **not in
  `batch.py`'s `_RUN_STEPS`**.

So every headless run computed this and discarded it.

**Fix:** `_ensure_centroids()` memoises on first use; `_cadh_vertices()` is the
single entry point for outside consumers, materialising centroids, normals and
`keep_mask` together so callers need not know anything is lazy. `_compute_mask`
materialises only when the band filter will genuinely use it (`band_um > 0` and
a DAPI mesh exists); at the default `band_um = 0` it stays lazy.
`cell_mesh.py`'s refine fetches **inside its worker**, so a first-use
computation cannot freeze the GUI.

**No behaviour change** — verified: after a run the fields are `None`; the
accessor returns (549, 3) centroids, (549, 3) normals and a consistent
`keep_mask`; the second call is cached at 0.0000 s.

If it ever needs computing often, the `eigh` calls should be **batched** —
measured 1.44 s looped against 0.33 s batched for 200k 3x3 matrices, a free
4.4x on that portion.

**Related:** [I08](#i08) (the sole consumer is itself a candidate for removal),
[P19](#p19), [P28](#p28).

---

<a id="p35"></a>
## P35 — `alpha 3D` hangs forever on a permissive mask, with no warning

| | |
|---|---|
| **Status** | open |
| **Severity** | S3 (no wrong number — it never returns one) |
| **Opened** | 2026-08-24 |
| **Last verified** | 2026-08-24 (reference image) |
| **Where** | `vessel_mesh.py::build_mesh` — `Delaunay(voxel_subsample(mask_to_points_um(...)))` |
| **Affects** | `alpha 3D (points)` **and** `alpha 3D solid -> poisson` — `rasterize_alpha_solid` Delaunays the same cloud |

`alpha 3D` does not see components; it sees **points**, one per occupied
`subsample_um` bin. `voxel_subsample`'s docstring states the design range:
"caps the count to ~1e4–5e4". A permissive cadherin mask blows straight through
it, because speckle is spread over the whole volume rather than concentrated on
a sheet — so the occupied-bin count saturates toward the *total* bin count
instead of scaling with wall area.

Reference image, extent 105.3 x 332.8 x 332.8 µm (~98,700 bins at
`subsample_um` 5):

| cadherin setting | voxels | components | points | % of all bins |
|---|---|---|---|---|
| `sauvola_k` 0.2, `floor_k` 5.0 (working) | 1,508,959 | 170,236 | **19,426** | 19.7% |
| `sauvola_k` 0.2, `floor_k` 3.0 | 5,157,789 | 1,335,402 | 85,500 | 86.6% |
| `sauvola_k` 0.43, `floor_k` 2.4 | 9,693,371 | 3,364,112 | **92,894** | **94.1%** |

At 94% occupancy the input is **not a surface** — it is a solid block filling
the field of view — so alpha shape has nothing to extract even in principle,
and any `alpha3d_um` above the bin spacing bridges everything to everything.
The points also sit on the 0.325 µm voxel lattice (`voxel_subsample` keeps an
original voxel, not a bin centre), and a regular lattice is near worst-case
input for Qhull: massively cospherical, so degeneracy handling degrades badly.
Together: it hangs, indefinitely, with no output.

**`min_size` is the lever, and it is dramatic**, because the offending material
is 3.35M specks that contribute little mass but occupy most of the bins:

| `min_size` | voxels | components | points | % of bins |
|---|---|---|---|---|
| 0 | 9,693,371 | 3,364,112 | 92,894 | 94.1% |
| **12** | 5,147,101 | 10,262 | **14,906** | **15.1%** |
| 200 | 4,870,230 | 458 | 13,353 | 13.5% |
| 2000 | 4,654,348 | 59 | 12,148 | 12.3% |

Wide plateau from 12 upward — [P19](#p19)'s point, visible directly.
`gui/cadherin.py`'s `MIN_SIZE_FLOOR = 12` exists precisely for this and its
comment already says the auto value "wrecks the downstream vessel alpha shape".

**A second defect surfaced by the same run:** `auto_min_size_otsu` returns **2**
at these settings, against the ~5 its comment records at stricter ones. **The
auto-seed falls as the mask gets more permissive** — it degrades exactly when
it is needed most, because speckle comes to dominate the log-size distribution
it is cutting. It cannot be relied on to protect this path.

**The fix wanted here is a guard, not a parameter.** `build_mesh` should refuse,
or warn loudly, when the subsampled point count exceeds the documented range,
instead of hanging silently. Turning an unbounded hang into a message naming
`min_size` is cheap and is the whole of the defect.

> **Guard landed 2026-08-24.** `vessel_mesh.check_point_cloud` (pure, testable
> against a count with no volume built) plus `_guard_points` at **both**
> `Delaunay` entry points — `alpha 3D` and `alpha 3D solid -> poisson`, which
> share the cloud. It raises `PointCloudTooDense`.
>
> **The refusal is on bin occupancy, not on the count**, because occupancy says
> *what* went wrong rather than merely how big it got: a surface touches a small
> fraction of the volume's bins, and something touching most of them is a solid
> block from which no `alpha3d_um` can extract a surface. Threshold 50%, against
> 19.7% measured working and 94.1% measured hanging. The count keeps a separate
> soft **warning** at 50,000, `voxel_subsample`'s documented ceiling.
>
> Verified both ways on synthetic volumes: a hollow spherical shell passes at
> 8.9% occupancy and meshes normally; a solid block is refused at 100%. The
> module stays **print-free** — notes go into `extras` for the caller, the
> `measure_core`/`mesh_score` contract — and the refusal carries its reason in
> the exception, so a caller passing no `extras` still learns why.

Note this does **not** argue for a permissive mask plus aggressive despeckling:
[P28](#p28)'s 2026-08-24 History records that the permissive input produces a
worse *surface* regardless, for reasons `min_size` does not touch.

**Related:** [P19](#p19), [P28](#p28), [P20](#p20), [I06](#i06).

---

<a id="p36"></a>
## P36 — The cadherin "local threshold" is inert: the mask is a global threshold

| | |
|---|---|
| **Status** | open |
| **Severity** | S2 (S1 if the depth bias below is confirmed) |
| **Opened** | 2026-08-24 |
| **Last verified** | 2026-08-24 (reference image, bit-identical over 6 parameter combinations) |
| **Where** | `pipeline/denoise.py::sauvola_mask` — `r = max(0.5 * (vol.max() - vol.min()), 1e-6)` |
| **Affects** | every cadherin mask, and therefore every vessel surface, centerline graph and geometry column derived from one |

**`sauvola AND floor` is bit-identical to `floor` alone**, at every setting
tested. The Sauvola local threshold contributes nothing:

| `sauvola_k` | `floor_k` | floor only | Sauvola AND floor | identical |
|---|---|---|---|---|
| 0.20 | 5.0 | 1,508,959 | 1,508,959 | **True** |
| 0.43 | 5.0 | 1,508,959 | 1,508,959 | **True** |
| 0.20 | 3.0 | 5,157,789 | 5,157,789 | **True** |
| 0.43 | 3.0 | 5,157,789 | 5,157,789 | **True** |
| 0.20 | 2.4 | 9,693,371 | 9,693,371 | **True** |
| 0.43 | 2.4 | 9,693,371 | 9,693,371 | **True** |

`sauvola_k` is a **dead knob**: 0.2 and 0.43 give the same mask to the voxel.

**The cause is the normaliser.** Classic Sauvola divides by `R`, the dynamic
range of grey levels, so that `s/R` lands in [0, 1] and the adaptive term
actually modulates. This uses **half the global intensity range**, set by the
single brightest voxel in the stack:

| | reference image |
|---|---|
| volume range | 170 – 7,725 → **R = 3,777.5** |
| median local mean `m` | 493.4 |
| median local std `s` | **16.1** |
| median `s/R` | **0.0043** |

The adaptive term therefore contributes 0.4% of `k`, and the formula collapses
to ``thr ≈ m·(1 − k)`` — measured at 1.000 / 0.900 / 0.801 / 0.701 / 0.572 for
k = 0 / 0.1 / 0.2 / 0.3 / 0.43. At the default k = 0.2 that is **0.8 × the local
mean, which is below the background level itself**: background over the 51-voxel
window is ~490 ± 16, so a voxel must fall ~6σ below its own neighbourhood mean
to be rejected. **97 voxels out of 339,738,624 manage it.** At k ≥ 0.3 it passes
all but one voxel in the volume.

The brightest voxel is ~15x the typical local mean and alone sets `R`, so this
is also not reproducible across images with different dynamic range — the same
`k` means a different threshold depending on each volume's hottest pixel.

**Two consequences.**

1. **[P28](#p28)'s exoneration of the local threshold needs rewording, not
   retracting.** It records that the local threshold deletes **0.000** of the
   gap-wall peaks and concludes it is "innocent". The number is right; the
   reading is not. It deletes 0.000 because it deletes *essentially nothing at
   all*. "Correctly permissive" and "inert" are different claims and only the
   second is true. P28's attribution of the loss to the floor and `min_size`
   stands — indeed it is now the whole story.
2. **The cadherin mask therefore has no local adaptivity, and inherits the
   defect [P03](#p03) describes.** One scalar over a stack with a systematic
   depth gradient necessarily over-accepts shallow and under-accepts deep. This
   is exactly the position the clover gate was in before [I07](#i07), where the
   measured acceptance ratio was **3.06** shallow-third to deep-third. **If it
   holds here, a substantial part of P28's missing 59–85% of dim wall is not
   "dim wall" at all but the deep wall, lost to a global cut.**

**The test, not yet run**, is the one [I07](#i07) used: accepted fraction of the
floor, shallow third vs deep third, on the substrate channel. It needs no
reconstruction and takes minutes.

**The fix** is to give `sauvola_mask` a normaliser on the scale of the local
variation rather than half the global range — a robust range (e.g. 1st–99th
percentile) at minimum, since a single hot voxel currently sets it. **Do not
change it silently**: every cadherin mask, every vessel surface and every
geometry column on disk was built through the current behaviour, so a change
invalidates them, and `sauvola_k` becoming live means the recorded 0.2 stops
meaning what it meant.

**Niblack is NOT inert — only Sauvola is**, and the batch path uses the broken
one. `niblack_mask` is ``vol > m + k·s`` with **no global normaliser at all**,
so nothing degenerates:

| method | k | passes alone | % of volume |
|---|---|---|---|
| sauvola | 0.20 | 339,738,527 | **100.0000%** |
| sauvola | 0.43 | 339,738,623 | **100.0000%** |
| niblack | 1.00 | 46,503,236 | 13.69% |
| niblack | 2.00 | 9,443,248 | 2.78% |
| niblack | 3.00 | 2,034,588 | 0.60% |
| niblack | 4.00 | 639,401 | 0.19% |

ANDed with the floor, Niblack does substantial work where Sauvola does none:
at `floor_k` 5.0 it takes 1,508,959 → 980,683, and at `floor_k` 3.0
**5,157,789 → 1,605,014 (69% removed)**.

**The drift makes this worse than a dead knob.** `gui/cadherin.py:200` sets the
tab's interactive default to **Niblack**, but `:379` — `_run_cadherin_step`,
which is what `state["run_cadherin"]` and therefore `batch.py`'s pipeline step
call — **forces Sauvola**. So the method you tune with works and the method that
builds every saved artifact does not. `denoise.py`'s own docstring records the
consequence without drawing it: "every cadherin mask in this project was built
with Sauvola". This was noticed earlier as a tidy-up; it is not one.

**A lead worth testing, and the first in this area not already measured to
fail.** Two compositions reach almost the same voxel count by opposite routes:

| | voxels |
|---|---|
| shipped: `floor_k` 5.0 + inert Sauvola | 1,508,959 |
| alternative: `floor_k` 3.0 + Niblack k=3.0 | 1,605,014 |

The shipped mask is a strict global cut with no local suppression. The
alternative admits far more — including the dim wall [P28](#p28) is about — and
strips speckle **locally**. That is structurally the same fix [I07](#i07)
applied to the clover gate. **Every previous attempt at a more permissive
cadherin mask was made with the local threshold switched off**, which is why
they all degraded the reconstruction.

**This is not a claim that it is better.** Equal counts say nothing about
*which* voxels, and by [P28](#p28)'s precision-over-recall asymmetry only the
surface can judge. The test is `mesh_score`'s four numbers on a reconstruction
from each composition, plus the shallow/deep acceptance ratio.

**Related:** [P28](#p28), [P19](#p19), [P03](#p03), [I07](#i07), [P26](#p26).

---

<a id="p37"></a>
## P37 — The loader overwrote `floor_k` on every image load, so the tab default never ran

| | |
|---|---|
| **Status** | **fixed** 2026-08-24 |
| **Severity** | S2 |
| **Opened** | 2026-08-24 |
| **Resolved** | 2026-08-24 — `gui/loader.py` no longer calls `seed_cadh_floor`; it is a button on the Cadherin tab |
| **Where** | `gui/loader.py` (the load path) calling `gui/cadherin.py::_seed_floor` |

Caught mid-rebuild: the Cadherin tab's `floor_k` default said **1.0** and the
run logged **`floor_k=4.595517`**. `_seed_floor` ran from the loader on every
image load and replaced whatever was there with `suggest_floor_k(raw)` — so the
shipped default never reached a build, and a value typed for the previous image
was silently overwritten when the next one loaded.

Same shape as [P26](#p26): **the recorded default is not what runs.** It is
worse than a stale default, because nothing in the settings JSON distinguishes
"the user chose this" from "the loader imposed it".

Auto-seeding was defensible while the floor was the **primary** gate — which it
was for exactly as long as the local threshold was Sauvola and therefore inert
([P36](#p36)). With Niblack live the floor is a weak backstop and its value is a
deliberate choice, so it is no longer imposed. `suggest_floor_k` remains, as an
explicit **"Suggest floor_k"** button that prints what it set.

**Related:** [P26](#p26), [P36](#p36), [P28](#p28).

---

<a id="p38"></a>
## P38 — A failed vessel-mesh build left the previous image's lumen field in `state`, and the centerline graph used it

| | |
|---|---|
| **Status** | **fixed** 2026-08-24 |
| **Severity** | **S1** — a saved artifact was computed from a different image |
| **Opened** | 2026-08-24 |
| **Where** | `gui/vessel_mesh.py::_done` (the failure early-returns) and `batch.py::_free_image_memory` |
| **Affects** | `vessel_graph.npz` for any image whose mesh build failed inside a multi-image batch |

`state["vessel_lumen_field"]` is set in `_done` **after** the empty/no-faces
early returns, and it was **not** in `batch.py`'s inter-image clear list. So a
failed build left the *previous* build's field in `state`, and the centerline
graph — which reads exactly that field — built on it and saved under the failed
image's name.

**Observed, not theorised.** In the 2026-08-24 six-image rebuild,
`V3_..._09step_001`'s mesh refused:

> `WARN: no lumen core found — 1 candidate region(s), none of which both avoided
> the z faces and reached 20.0 µm from a wall. ... Refusing to guess an
> orientation.` → `No points/vertices to mesh.`

and the very next step reported `Centerline graph: 418 points on 1 edge(s), 0
junction(s), ... 409 µm of vessel` and saved it. Those numbers are image 2's
vessel. **Nothing in the output looks wrong** — that is the danger: a plausible
graph, on the wrong image, under the right filename.

`vessel_poisson.lumen_core`'s refusal did its job. The bug is that a refusal
downstream of it was not propagated.

**Fixed in both places**: `_fail` now clears `vessel_lumen_field`, and
`batch.py` clears it (and `vessel_mesh_cache`) between images.

The contaminated artifact is quarantined as
`vessel_graph.npz.CONTAMINATED_P38` rather than deleted, so the failure stays
inspectable.

**Two things this exposes, both still open.**

1. **The new cadherin defaults break `lumen_core` on image 3.**
   `docs/DIRECTION.md` §8A.3 recorded that it "refuses on none of the six" —
   that was with the old, stricter mask. With Niblack k=2.0 / floor_k=1.0 the
   lumen no longer reaches 4 x `poisson_smooth_um` on `09step_001`. So the
   rebuild produced 5 usable meshes, not 6. Either that image needs its own
   `poisson_smooth_um`, or the defaults are wrong for it — and per
   [P26](#p26) a per-image parameter has to be recorded, not eyeballed.
2. **`batch.py`'s mesh step writes no `vessel_recon.json`.** Only
   `batch_mesh.py` does. So a mesh built through `batch.py` is unprovenanced —
   and worse, it can *inherit a stale JSON from an earlier `batch_mesh` run*,
   which is what happened here: a 20:09 mesh sat beside a 12:38 JSON naming a
   different method. `saved_check` warns when the JSON is **missing**, so this
   case passes the check while lying. The stale file was deleted; the writer is
   not yet fixed.

**Related:** [P26](#p26), [P32](#p32), [I04](#i04).

---

<a id="i07"></a>
## I07 — Local thresholding on the clover channel

**Opened** 2026-08-24 · **Status** open · **Fixes** [P03](#p03)

`denoise.otsu_mask` is `threshold_otsu(vol)` — **one scalar for the entire 3D
stack** — and it is what builds the clover mask (`gui/clover.py:80`) and the
adhesion gate (`gui/adhesion_rf.py:261`). On a channel with a systematic depth
gradient, a single global number necessarily over-accepts shallow and
under-accepts deep. [P03](#p03)'s 2026-08-24 History confirms that gradient is
optical rather than biological, so this is the root cause of an S1 defect and
not a robustness nicety.

**The tools already exist and are already correct for this data.** `denoise.py`
has `niblack_mask`, `sauvola_mask` and `bernsen_mask`, all local, all taking
`voxel_zyx_um` so they scale the z window for the 0.325/0.4 (and 0.9)
anisotropy. They are used on the cadherin channel. **Nothing on the clover path
calls them.**

**This has an unusually concrete success criterion**, which is rare in this
repo: [P03](#p03)'s shallow/deep ratio at fixed wall orientation (1.40 / 0.74 /
2.62 / 0.93 / 2.12 / 3.00) should move toward **1.0** if the fix works. Rerun
`docs/scratch/tier2/orientation_vs_depth.py` against the new instances. The
far-wall-vs-near-wall figure in [P03](#p03) is the coarser version of the same
check.

**Measured 2026-08-24, reference image** (`docs/scratch/tier2/{threshold_depth,
clover_gate_sweep}.py`, isotropic grid — the RF, the gate and
`adhesion_instances.npz` all live there, 324 planes against 263 native).

**(a) The gate's depth bias is real and the local threshold removes it.**
Accepted voxel fraction, shallow third vs deep third:

| gate | shallow | deep | sh/dp |
|---|---|---|---|
| **global Otsu (current)** | 2.141% | 0.698% | **3.06** |
| niblack k=+0.0 | 40.333% | 42.238% | 0.95 |
| niblack k=+0.5 | 22.123% | 24.247% | **0.91** |
| niblack k=+1.0 | 11.876% | 13.173% | 0.90 |
| niblack k=+3.0 | 0.962% | 0.937% | 1.03 |

The Otsu ratio **3.06** independently corroborates this entry's parent
[P03](#p03), whose far-wall-vs-near-wall figure is 37% (1/0.37 = 2.7) by a
different route. **Depth-uniformity is a property of the method, not of k** —
every k lands at 0.90–1.03 — so k can be chosen purely on recall vs cost. At
k ≤ 1.5 the gate is slightly *more* permissive deep, partially compensating for
attenuation rather than merely not amplifying it.

**(b) k must be permissive, and the cadherin default is catastrophic here.**
The clover mask is **a gate on the RF, not a segmentation**: its job is to
discard only what is certainly not an adhesion. Recall of the 992,011 voxels in
`adhesion_instances.npz`:

| gate | volume retained | recall |
|---|---|---|
| global Otsu (current) | 0.969% | 94.82% |
| niblack k=+0.0 | 42.275% | 99.98% |
| niblack k=+0.5 | 24.254% | **98.85%** |
| niblack k=+1.0 | 13.181% | 91.63% |
| **niblack k=+3.0** (cadherin default) | 0.921% | **36.57%** |

**k=3.0 is the worst available choice** — the same volume as Otsu for a third
of the recall. It was briefly the shipped default on the strength of its 1.03
depth ratio alone; uniformity at that k is uniform *discarding*, and
uniformity is not a sufficient criterion on its own.

**Recommended: `niblack_k = 0.5`, `window_xy = 50`** — 98.85% recall, 0.91
depth ratio, 24% volume. The cost is ~25× more voxels through
`rf.predict_proba`, and that runtime is the only real argument for a stricter k.

**Two limits on (b), both structural.** The reference instances were produced
*through* the Otsu gate, so they cannot contain anything it excluded: recall
against them is a necessary condition, never sufficient, and it is blind to
exactly the deep population [P03](#p03) says is missing. Otsu's own recall is
94.82% rather than 100% only because the gate is dilated 1 µm before use. And
Otsu's recall is **depth-flat (95.41 / 93.51)**, which looks reassuring and
means nothing — adhesions it excluded were never detected, so they are not in
the set to be missed. **The comparison is rigged in Otsu's favour** and must not
be read as "Otsu is better"; niblack keeps different voxels, losing Otsu-found
ones while gaining ones this metric cannot see. Only re-detection settles it.

**Cautions.**

- **[P28](#p28) is not evidence against this.** P28 measured a continuous
  local-contrast *field* degrading the cadherin surface fit. That is a
  different use — a field feeding a reconstruction, not a mask replacing a
  global threshold — and the two should not be conflated. Someone will raise
  it; the answer is on the record here.
- **[P19](#p19) applies directly.** `auto_min_size_otsu` + `filter_by_size` run
  immediately after the threshold on this same path, and P19 records that
  `min_size` deletes 59–85% of the dim wall on the cadherin channel. A local
  threshold that recovers dim deep adhesions will hand them straight to a size
  filter tuned against the global-threshold population. **Both have to be
  looked at together or the fix will be masked.**
- Per the working agreement the knob is a slider in the Clover tab with an
  argparse equivalent of **no upper bound**, and the GUI and headless paths
  call one shared function in `pipeline/` — the `build_mesh` doctrine. The gate
  in `adhesion_rf` and the mask in `clover` must not acquire two
  implementations.
- Changing this **invalidates every existing adhesion CSV**. The
  `20260812_011*` exports are the current authoritative six; a re-measure is
  part of the work, not a follow-up.

**History**

- 2026-08-24 — opened, measured, and **the Clover tab default changed from
  global Otsu to Niblack k = 0.5** (`window_xy` 50) on the user's decision.
  `denoise.build_mask` is now the single dispatch for threshold + label +
  size-suggest, and both `gui/clover.py` and `gui/adhesion_rf.py`'s gate
  fallback call it — they each had their own copy of
  `otsu_mask -> label_and_sizes -> auto_min_size_otsu`, so the mask the user
  tuned and the mask the RF gated on could differ silently. The RF gate path
  itself is unchanged: `adhesion_rf.py` still reads `clover_mask_layer.data`,
  dilates it by `clover_dilation_um`, and passes `gate_indices` to
  `rf.predict_proba`.

  **`batch.py` runs every tab at widget defaults and has no CLI for them, so
  this default is now the headless setting**, and the `20260812_011*` CSVs are
  superseded in principle even though nothing has been re-run yet. **The fix is
  not yet demonstrated** — everything measured above is upstream of the RF.
  Confirming it needs re-detection (RF predict -> post-process -> re-measure)
  and then rerunning `docs/scratch/tier2/orientation_vs_depth.py` to see
  whether P03's shallow/deep ratios move toward 1.0. Equal *acceptance* is not
  proof the right voxels were accepted.

  **The [P19](#p19) interaction was checked and is benign at this setting.** The
  gate reads the *size-filtered* mask, so `min_size` is a second strictness knob
  on the same gate, and a permissive threshold changes the component structure
  completely — 12,247,014 components at k = 0.5 against Otsu's 89,303. It could
  therefore have clawed back the permissiveness. It does not:
  `auto_min_size_otsu` seeds **3**, not something large, taking volume 24.254%
  -> 20.955% for a recall cost of **0.01 points** (98.85% -> 98.84%) with the
  depth split unchanged. That is a size filter removing 1-2 voxel speckle, which
  is what it is for. No override of the auto-seed is needed on this tab.
  **This says nothing about P19 on the cadherin channel**, where the filter is
  measured to delete 59-85% of the dim wall; the two are the same code and
  different regimes.

- 2026-08-25 — **closed: the premise is refuted, though the sub-measurements
  stand.** This entry's claim was that the global Otsu gate is "the root cause of
  an S1 defect", [P03](#p03). The clover tab did move to niblack (commit
  `9aec319`), and this entry's depth-ratio table (Otsu 3.06 vs niblack ~0.91) is
  still correct as far as it goes. But it measured the *gate*, never the
  *output*: re-running detection on image 6 with the gate at 93.4% of the volume
  produced a bit-identical `proba >= 0.5` set ([P43](#p43)). The gate was never
  binding, so its depth bias never reached the result.
  **Two of this entry's own cautions turn out to have been the important part.**
  It warned that recall was measured against instances produced *through* the
  Otsu gate and so was "rigged in Otsu's favour", and that "only re-detection
  settles it". Re-detection settled it the other way. And its rejection of
  `k=3.0` as "uniform *discarding*" has a mirror it did not consider: `k=0.5` is
  uniform *passing* — a fixed noise quantile is depth-uniform by construction
  whether or not any signal is present, which is what the 0.91 ratio was
  measuring. See [P43](#p43) for the percolation argument.
  The gate is now removed entirely, so this idea has no remaining consumer.

---

<a id="i08"></a>
## I08 — "Refine with cadherin" is a dead end and should be removed

**Opened** 2026-08-24 · **Status** open · **Type** cleanup, not a defect

`gui/cell_mesh.py`'s "Refine with cadherin" button adds cadherin component
centroids as extra vertices to the **DAPI cell mesh**, where nuclei are too
sparse to close the surface: keep cadherin points within `refine_band_um` of
the existing mesh, drop any within `min_sep_um` of a nucleus, rebuild the alpha
shape.

**It exists to get a better vessel surface out of an alpha shape over cell
vertices. That is no longer how the vessel surface is built.**
`alpha 3D solid -> poisson` ([I06](#i06)) fits the wall directly and is the
first reconstruction to pass [P24](#p24)'s plateau test; the centerline graph
([P27](#p27)) and `nuc_side` are built on the same solved field. The nucleus-
derived alpha shape is a proxy that has been superseded.

**Why it is worth actively removing rather than leaving:**

- It is the **only** consumer of the Cadherin tab's centroids/normals, which
  [P34](#p34) has just had to make lazy to keep that tab usable. Remove it and
  `label_centroids_normals` can leave the cadherin path entirely.
- The mesh it refines (`mesh_cache`) is also what gates adhesions by distance
  (`gui/clover.py`'s `gate_by_mesh`/`d_max_um`, and the RF post-process gate).
  **That gate should measure against the reconstructed vessel surface**, for
  which there is now an exact signed nearest-*triangle* distance in
  `pipeline/nuc_side.py`, rather than against an alpha shape over nucleus
  centroids. Doing so would share machinery with `docs/DIRECTION.md` §8C.4
  instead of maintaining a second, weaker surface.
- It is not in `batch.py`'s `_RUN_STEPS`, so nothing headless exercises it and
  it is not covered by any check.

**Recorded on the user's judgement, 2026-08-24:** the project should be
clearing dead ends rather than adding to them. This is the first one named.
Removal is not done here — it touches the cell-mesh tab, the clover gate and
`mesh_cache`, and the replacement gate is a design decision in its own right.

**Related:** [P34](#p34), [I06](#i06), [P09](#p09), `docs/DIRECTION.md` §8C.

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

---

<a id="h02"></a>
## H02 — `adhesion_rf.joblib` is a pickle against an unpinned scikit-learn

**Opened** 2026-08-11 · **Status** open (measured harmless *today*)

`environment.yml` pins `napari=0.6.0` and `vispy=0.14.3` but lists
`scikit-learn` with no version. The shipped `adhesion_rf.joblib` was trained
under **scikit-learn 1.9.0**; the WSL development env currently resolves to
**1.8.0**, so loading it prints:

> `InconsistentVersionWarning: Trying to unpickle estimator RandomForestClassifier
> from version 1.9.0 when using version 1.8.0. This might lead to breaking code or
> invalid results.`

A joblib pickle of a fitted estimator is not a portable model format — it is a
snapshot of sklearn's internal objects. There is nothing in the repo that would
catch a future version actually changing them, and the failure would be
**silent**: a forest that loads but traverses differently produces a plausible
probability map, not an exception.

**Measured on 2026-08-11 — no effect on the current pair.** Features were
recomputed at 30,000 gated voxels of `V3_..._04step_002` and predicted with the
1.8.0-loaded forest, then compared against the stored `adhesion_proba.npy`:

| | |
|---|---|
| byte-identical probabilities | 29,935 / 30,000 (**99.78%**) |
| max \|difference\| | 7 / 255 |
| **voxels flipping the `>= 128` segmentation threshold** | **0** |

And the residual 0.22% is **not** an sklearn effect. `predict_proba` computes
features on 20-plane z-chunks padded by `ceil(4σ)+1`; the check recomputed them
over the whole volume, and a 4σ truncation differs slightly from an untruncated
filter near a chunk edge. The differing voxels are **6.7× enriched within one
plane of a chunk boundary** (0.92% at the edge plane, decaying monotonically to
0.00% by six planes in) — chunk geometry, which a version effect could not
respect. A random forest's `predict_proba` is a mean over leaf class fractions;
it has no reason to know where a z-chunk starts.

**So the five-image batch is safe to run in either environment as things stand.**

**What to do:** pin `scikit-learn` in `environment.yml` to the version the shipped
model was trained under, or record that version next to the model and re-verify
after any upgrade. Better still, treat the check above as the regression test —
it is cheap (one image, ~4 minutes) and it is the only thing standing between a
silently-degraded classifier and five images of output.

**Related:** [P12](#p12) (same failure class: an unattended run that produces
plausible numbers with no message).

---

Also noted under [H01](#h01): `3d_images_folder/Channel Info.rtf` was removed from
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

> **2026-08-11 — [P04](#p04) is fixed, so `(1|image/nucleus)` is now available**
> on any CSV exported after commit `04262b8`: `nuc_id` indexes nuclei rather than
> DAPI components, and `nuc_is_fibroblast` says which class each belongs to. On
> `V3_..._04step_002` that is 2,147 adhesions across **170 nuclei**, median ~8
> adhesions per nucleus.
>
> Three cautions before using it as the grouping factor:
>
> - **A nucleus is not a cell.** `nuc_id` is still nearest-nucleus assignment,
>   with no cell segmentation and no maximum distance. It is a better proxy than
>   before, not a measurement of cell identity.
> - **Do not stratify on `nuc_is_fibroblast`** — the classifier behind it is
>   unvalidated and swings ~2× on one slider step ([P05](#p05)). As a grouping
>   variable it is fine; as a *biological* split it is not yet.
> - **The independent variable is still the problem.** Adding a cell level
>   improves the error structure; it does nothing for the fact that neither
>   curvature column currently measures the vessel ([P24](#p24),
>   [P08](#p08)(e)).

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

> **2026-08-22 — the column names in this table are pre-rewrite.**
> `path_curv_xy` and `dist_branch_um` were deleted with the 2D centerline; the
> measurements below are left exactly as taken, because they are evidence about
> what was measured on that date. Nothing here has been re-run against the graph
> columns, and it should be before any of it is quoted again.

> **2026-08-11 — there is currently no trustworthy geometry variable in this
> table.** `path_curv_xy`, the strongest apparent correlate at +0.191, is fitted
> over a 20 µm window on a vessel ~140 µm wide, so it largely measures
> medial-axis wobble ([P08](#p08)(e)); the four `curv_*` columns have no
> scale-independent value ([P24](#p24)). Both are now known to be set by their
> own fit scale rather than by the vessel. Until [P21](#p21)/[I03](#i03) lands,
> **the geometry side of "tension vs geometry" has no defensible independent
> variable**, and the correlations above should be read as a description of one
> spatial gradient, not as candidate results.

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
([P09](#p09)), a distance-to-junction gate, a field-of-view margin, the depth
band, and whether `09step` (0.9 µm z-spacing) is in. Each is otherwise a fork in
the garden.

> **Updated 2026-08-22 for the column rename.** `dist_branch_um` and
> `path_curv_xy` no longer exist; the 2D centerline was replaced by the 3D graph
> (`docs/DIRECTION.md` §8B/§8C). The equivalents are **`dist_junction_um`**
> (distance *along the vessel*, not Euclidean in projection) and
> **`path_curv_1_um`** / `path_curv_xy_1_um`.
>
> The field-of-view margin is now a **column rather than a hand-rolled filter**:
> `dist_crop_um` is the distance along the vessel to where the image cuts it,
> and `censored` is already 1 where the lumen is not fully imaged. Prefer
> `censored == 0` to a margin you pick yourself — it scales with the local
> caliber, which a fixed µm margin does not.

**Guard the junction filter.** Not every vessel here is Y-shaped — three of the
six are a single tunnel. Those images have no junction, so `dist_junction_um` is
all-NaN and a bare `dist_junction_um > 20` drops **every** row rather than none.
Write it as `dist_junction_um.isna() | (dist_junction_um > 20)`. The mechanism is
[P08](#p08)(d)'s and it survived the rewrite unchanged, because it is a property
of NaN and not of the old code.

## State the ceiling

With six images, if all agree in direction a sign test floors at p = 0.031. That
is the honest best case.

- **Defensible, eventually:** "In all six images, adhesion tension rank increases
  with axial vessel curvature."
- **Not defensible, and will not be until FRET is calibrated:** any slope in pN
  per unit curvature.

<a id="p39"></a>
## P39 — 7 days may be short of a confluent monolayer, so some cadherin gaps may be real

| | |
|---|---|
| **Status** | open — hypothesis, not yet tested |
| **Severity** | S3 |
| **Opened** | 2026-08-25 |
| **Last verified** | — |
| **Resolved** | — |
| **Where** | the dataset, not the code; bears on `pipeline/denoise.py` mask tuning |
| **Affects** | the interpretation of [P28](#p28) / [P19](#p19), and vessel-surface coverage |

The construct was fixed after **7 days** (`docs/DIRECTION.md` §3). The reference
protocol this method descends from cultures endothelialised microvessels for
**1–2 weeks** to reach "a confluent and continuous endothelium"
(`docs/papers/ncomms8858.pdf`, Methods). Seven days sits at or below the short
end of that.

[P28](#p28) and [P19](#p19) established that the cadherin wall's apparent gaps
are mostly *not* holes — mask-gap directions show the same peak radius as
mask-present ones at ~40% the brightness, and the global floor plus `min_size`
delete 59–85% of the dim wall. That evidence stands and is not in question here.

What is in question is whether it accounts for **all** of the gaps. If the
monolayer was still filling in at fixation, some fraction of the missing wall is
genuinely absent rather than dim, and no threshold change will recover it.
Chasing those with a more permissive mask would add speckle for nothing.

**Why it matters even though it is S3.** It is an alternative hypothesis that
nobody is currently holding, and it bounds how much [P28](#p28) can be expected
to buy. It also changes what a coverage score means: `mesh_score`'s `coverage`
asks "is the surface everywhere the data is", which is the right question only
if the data is everywhere the *vessel* is.

**How to test it.** Cheapest first:

- Ask the experimenter whether the endothelium was confluent at 7 days, and
  whether that was checked (`DIRECTION.md` §9). One sentence closes this.
- Look at whether the dim-wall directions in [P28](#p28)'s analysis have a
  bimodal residual — dim-but-present versus genuinely empty — rather than a
  single dim population.
- A real hole should be **contiguous across z** and bounded by bright wall on
  all sides; thresholding loss should be scattered and depth-correlated
  ([P03](#p03)).

**Related:** [P28](#p28) (the continuous field, and the 59–85% dim-wall loss),
[P19](#p19), [P03](#p03) (depth bias, the competing explanation),
[I03](#i03) (what coverage is validated against).

**History**
- 2026-08-25 — opened while recording `docs/papers/ncomms8858.pdf` against this
  project's protocol. Hypothesis only; no measurement.

---

<a id="p40"></a>
## P40 — `fibroblast_cellpose.py` thresholded noise as tissue, and the wall split degenerated silently

| | |
|---|---|
| **Status** | **fixed** 2026-08-25 (both defects); the no-sanity-gate part is **open** |
| **Severity** | S2 |
| **Opened** | 2026-08-25 |
| **Last verified** | 2026-08-25 |
| **Resolved** | 2026-08-25 — `_saved_wall_mask`, `_voxel_zyx`, `_split_z` guard |
| **Where** | `Scripts/fibroblast_cellpose.py` — `_tissue_mask`, `_split_z` |
| **Affects** | every endothelial count the script could produce; nothing else in the pipeline |

The script had **never been run** before 2026-08-25 and was cited in neither
`DIRECTION.md` nor this file. Running it surfaced two defects.

**(a) `_tissue_mask` thresholded the noise.** Otsu on a smoothed copy, then the
largest connected component, marked **24,484,433 of 61,603,840 voxels — 40% of
the volume** — as tissue on image 6. The per-z profile was consequently flat at
~40% on every plane. This is the same global-threshold failure [P36](#p36) and
[P28](#p28) record for this channel: the wall is a thin, dim minority of a
speckly volume, and a global cut lands inside the noise. The pipeline's own
Niblack mask on the same image is **2.7–3.0%** nonzero.

**(b) `_split_z` degenerated without saying so.** The valley search returned the
**endpoint** of the tissue range (`zs = 234` on a 235-plane stack — a legal
index, which is why it did not crash). At lines 235–236 that sends the whole
stack to one "wall" and a single plane to the other, so a **whole-volume MIP was
reported and plotted as an en-face wall projection**. A flat profile from (a)
guarantees this, but it can happen whenever the two walls are not resolved as two
humps in depth.

**The fix.** `_saved_wall_mask` reuses the pipeline's saved `cadherin_labels.npz`
rather than re-thresholding — this is CLAUDE.md's "never write the second copy",
and it also means the endothelial count refers to the *same* wall as the mesh and
the graph. Nearest-neighbour resampling from the isotropic 0.325 µm grid to the
script's grid (boolean masks must not be interpolated); voxel size read from the
OME header via a new `_voxel_zyx`, never hard-coded, since `09step` is 0.9 µm in
z. Otsu remains as a labelled fallback for an image with no saved mask.
`_split_z` now warns and falls back to the tissue-range midpoint when the minimum
lands on an endpoint.

**Still open: no sanity gate on the output.** The script wrote
`0 endothelial / 174 nuclei / 100% fibroblast` to `endothelial_counts.csv`
without complaint. A count of zero, or a fibroblast fraction at 0% or 100%,
should refuse to write — same class as [P25](#p25), a nonsense number reaching a
file. The garbage CSV produced during this session was deleted.

**Also cosmetic but misleading:** `_qc_overlay`'s title prints `int(m.max())`,
the highest **label id**, not the object count. It read "68 cells on this slice"
where there was one.

**Related:** [P41](#p41) (what the fixed script then showed), [P36](#p36),
[P28](#p28), [P25](#p25).

**History**
- 2026-08-25 — opened and both defects fixed in the same session, while trying to
  answer `DIRECTION.md` §8D.5's fibroblast count. See §8E.3.

---

<a id="p41"></a>
## P41 — the en-face VE-cadherin wall has no polygonal cells, so endothelial segmentation returns nothing

| | |
|---|---|
| **Status** | open — measured, not explained |
| **Severity** | S2 |
| **Opened** | 2026-08-25 |
| **Last verified** | 2026-08-25 (image 6 only) |
| **Resolved** | — |
| **Where** | the data; `Scripts/fibroblast_cellpose.py` is the tool that surfaced it |
| **Affects** | the endothelial-subtraction route to [P05](#p05); `DIRECTION.md` Tier 2 item 5a |

With [P40](#p40)'s defects fixed — tissue mask at 3.0% of volume, a genuine
interior valley at Z=205, no warning — **Cellpose returns zero objects on both
wall projections** of image 6. Still zero at `--min-vox 0 --cellprob -2`, so it is
not a size or probability filter deleting them.

Two separate observations from the QC projection:

- **The depth split is wrong for this geometry.** Image 6's vessel curves
  diagonally across the field, so it has no top and bottom wall at separable
  depths; the vessel landed almost entirely in one half. Doing this properly
  means **unwrapping the surface using the Poisson mesh the project already
  has**, not splitting by z.
- **Where the wall *is* seen en-face, the texture is fibrous filaments and
  speckle, with no closed polygonal tiles.** There is nothing cell-shaped to
  segment. This is independent of the split problem and is the reason the route
  fails.

**Two candidate explanations, neither tested.** [P39](#p39) — the monolayer was
not confluent at 7 days, so there are no complete junctional polygons to find; or
[P28](#p28) — the junctions are present but too dim and too dashed to close a
contour. They are distinguishable, and the same evidence P39 lists would separate
them.

**Do not tune this script further.** The blocker is image content, not
parameters, and eye-tuning a parameter that then has to be recorded is exactly
what [P26](#p26) warns against.

**One caveat not chased:** cellpose 4.x is the SAM rewrite and ignores
`model_type=cyto3`; 3.x might behave differently. That is a dependency downgrade
on speculation.

**Related:** [P39](#p39), [P28](#p28), [P05](#p05), [P40](#p40),
[I02](#i02).

**History**
- 2026-08-25 — opened. Image 6 only; the other five were not run, because the
  failure is structural rather than image-specific.

---

<a id="p42"></a>
## P42 — `nuc_id` is a nearest-neighbour join presented as ownership

| | |
|---|---|
| **Status** | open |
| **Severity** | **S1** |
| **Opened** | 2026-08-25 |
| **Last verified** | 2026-08-25 (all six images) |
| **Resolved** | — |
| **Where** | `pipeline/measure_core.py` — `nuc_assignment`, `dist_nuc_um` |
| **Affects** | `nuc_id`, `nuc_is_fibroblast` per adhesion; any per-cell aggregation |

[P04](#p04) is closed as fixed because cell type became **joinable** per adhesion
via `nuc_id`. It is joinable. It is not **correct**, and the difference matters
for any analysis that groups adhesions by cell.

`nuc_id` is the nearest nucleus. Measured on all six images:

- **Nucleus spacing is 19–25 µm** on the wall (√(area/N); 174–315 nuclei over a
  ~332 × 332 µm field), ~34 µm volumetric.
- A spread fibroblast footprint is **50–100 µm**, so its peripheral adhesions sit
  25–50 µm from their own nucleus while two or three other nuclei are nearer.
- The artifact is visible directly: per-nucleus adhesion groups have **median
  radius 16–17 µm in every image** — half the nucleus spacing, i.e. the radius of
  a Voronoi cell. The groups are a tessellation of the field, not cells.

**What this invalidated.** A cross-tab showing classifier-called fibroblast
nuclei owning 1.5–6× more adhesions than endothelial ones, consistent in
direction across all six images, was withdrawn on these grounds. It is confounded
a second way as well: the [P05](#p05) classifier thresholds perinuclear mClover
while the adhesions are detected in that same channel, so it scores the same
photons twice — §8C.5's structural trap in a new place.

**What survives and is worth measuring.** `dist_nuc_um` has median **4.5–6.2 µm**
across the six, which appears **below** what chance would give at this nucleus
density. That statistic does not care which nucleus is the correct one, so it
survives the defect. If it holds against a proper CSR null it says the detected
adhesions are **perinuclear** — backwards for focal adhesions, and consistent
with the detector fragmenting perinuclear cytoplasmic signal, which would reach
[P03](#p03) and [P02](#p02). **The null has not been computed**; the expectation
above was estimated by hand and must not be quoted as a result.

**What a fix would need.** Ownership requires a segmented cell body — sensor
signal grouped into cells, then adhesion-to-cell by containment rather than
proximity. Whether that is possible on this data is open: the fibroblasts are not
distinguishable by eye in the sensor channels (§8E.1). Until then, `nuc_id`
should be documented as a proximity join and **no analysis should aggregate
adhesions by nucleus**.

**Related:** [P04](#p04) (which this qualifies), [P05](#p05), [P03](#p03),
[P11](#p11), `DIRECTION.md` §8E.2.

**History**
- 2026-08-25 — opened after the adhesions-per-nucleus cross-tab was retracted.

---

<a id="p43"></a>
## P43 — the clover gate cost bias and saved ~10%; removed

| | |
|---|---|
| **Status** | **fixed** 2026-08-25 — gate default off, `chunk_z` auto-sized |
| **Severity** | S3 (compute), but it was the assumed cause of an S1 — see [P03](#p03) |
| **Opened** | 2026-08-25 |
| **Last verified** | 2026-08-25 (image 6) |
| **Resolved** | 2026-08-25 |
| **Where** | `pipeline/gui/adhesion_rf.py` (`gate_by_clover`), `pipeline/rf.py` (`predict_proba`) |
| **Affects** | RF inference cost; the premise of [I07](#i07) |

`gate_by_clover` restricted RF inference to voxels passing a threshold on the
clover channel. It was **only ever a compute shortcut** — but it decided *a
priori* which voxels were allowed to be positive, which is exactly the bias
[P03](#p03) is about.

**It did not save what it appeared to.** `predict_proba` chunks over **all** z
regardless of the gate, and `_extract_scale_*` filters the whole padded slab
densely before gathering at the gate's indices. Since adhesions span all z, a
sparse gate still touches every chunk and pays the same filtering. Per chunk:

| | time |
|---|---|
| dense feature filtering | ~80–100 s |
| RF `predict_proba` | ~9–12 s |

The gate only shrinks the second. **It saved ~10% of runtime.**

**And it was never binding.** Widening the gate on image 6 from Otsu's 0.969% of
the volume to niblack `k=0.5`'s 93.4% (post-dilation) left the `proba >= 0.5`
set **bit-identical — 0 voxels differed**:

| gate | nonzero proba | proba >= 0.5 | components | >=100 vox |
|---|---|---|---|---|
| Otsu, dil 1 µm | 4,577,402 | **843,386** | 2,026 | 829 |
| niblack k=0.5, dil 1 µm | 4,882,386 | **843,386** | 2,026 | 829 |

305,000 extra voxels got a nonzero probability and **not one reached 0.5**. The
resulting `adhesion_instances.npz` was byte-identical to the previous one, which
is the correct consequence rather than a stale-artifact bug (the proba files do
differ).

**Confirmed at full volume.** A third run with the gate off entirely (100% of
the volume, `chunk_z` auto-sized to 24) gives:

| | proba >= 0.5 | instances |
|---|---|---|
| niblack gate 93.4%, `chunk_z` 20 | 843,386 | 1872 (mean 438, min 100, max 2161) |
| **no gate, `chunk_z` 24** | **843,372** | **1872 (mean 438, min 100, max 2161)** |

**130 voxels differ — 0.015%.** They are scattered over 78 z-planes rather than
concentrated at chunk boundaries, and this run changed two things at once (gate
*and* `chunk_z`), so they are consistent with float-level chunking effects
rather than with the gate having withheld anything. `adhesion_instances.npz`
differs in 5,197 voxels, which is larger only because a handful of mask voxels
can merge, split or **renumber** whole components downstream; the instance count
and every summary statistic are unchanged.

**Runtime, as predicted:** gated 15 × ~89 s ≈ 22 min, ungated 13 × ~112 s ≈ 24
min — about 10% more for removing the bias.

**Two mechanisms had to be understood to get here, both worth keeping:**

- **The 1 µm dilation destroys a speckly gate.** niblack `k=0.5` selects 19.51%
  of the volume as **11.1 M components**; dilating by a 9×9×9 element expands it
  to 93.4%. `min_size` cannot rescue it (see below), so the gate became vacuous.
- **niblack at low `k` is a noise quantile, not a detector.** Foreground is
  `vol > mean + k·std`, so on pure noise it passes `P(Z>k)` by construction:
  30.9% at `k=0.5` against [I07](#i07)'s measured 24.25%. 3D site percolation on
  a simple cubic lattice (which is `ndi.label`'s default 6-connectivity) has
  **p_c ≈ 0.3116** — so at `k=0.5` the background condenses into one spanning
  cluster. That is why it "flags basically everything" in the viewer, why
  `min_size` is useless there, and why the dilation closes the rest.
  **Design rule: keep `P(Z>k)` well below 0.31, i.e. `k >= 1.5`.**

**The fix.** `gate_by_clover` now defaults **off**. Removing it makes the
feature matrix `X` the binding memory constraint — it is
`(chunk voxels × 60) float32`, i.e. 5.0 GB for a 20-plane chunk of a 1024²
volume — so `predict_proba` grew `_auto_chunk_z`, which sizes the chunk to fit
both a host budget (`HOST_BUDGET_GB`, default 8.0) and a share of free VRAM
(`VRAM_FRACTION`, default 0.8), and **logs the choice**. `X` cannot be streamed
in batches: the dense filtering is shared across every voxel in the chunk, so
streaming would recompute it per batch.

**Note for 12 GB cards** (the target this was sized for): the largest scale needs
a `2·max_pad = 126`-plane slab whatever `chunk_z` is, so there is a **~6.5 GB
device-memory floor** for a 1024² volume at 5 µm/0.325 µm that shrinking the
chunk cannot get under. `_auto_chunk_z` warns when it does not fit rather than
letting it OOM silently.

**Related:** [P03](#p03), [I07](#i07) (whose premise this refutes),
[I09](#i09) (where the runtime actually goes), [P36](#p36), [P19](#p19).

**History**
- 2026-08-25 — opened and fixed. Found while re-running detection to test
  [I07](#i07); see `docs/DIRECTION.md` §8F.

---

<a id="i09"></a>
## I09 — RF inference is dominated by `_eigh3x3` on the CPU

**Opened** 2026-08-25 · **Status** open · **Relates to** [P43](#p43)

Feature extraction is ~90% of RF inference cost. Per-scale timing on image 6,
one chunk (24 z-planes of a 1024² volume, GPU):

| scale | σ (voxels) | taps/axis | time |
|---|---|---|---|
| 0.3 µm | 0.9 | 9 | 16.4 s |
| 0.7 µm | 2.2 | 19 | 14.0 s |
| 1.0 µm | 3.1 | 27 | 16.5 s |
| 1.6 µm | 4.9 | 41 | 15.6 s |
| 3.5 µm | 10.8 | 89 | 16.9 s |
| 5.0 µm | 15.4 | 125 | 17.9 s |

**Taps vary 14×; time varies 1.3×.** The cost is per-filter-operation memory
traffic over the slab, not kernel arithmetic. Everything below follows from
that one table.

**Two optimisations were proposed and are rejected on this measurement** —
recorded so they are not re-proposed:

- **`truncate` 4.0 → 3.0** gives ~24% fewer taps at every scale. On something
  that is not tap-bound that is worth almost nothing, and it changes every
  feature value while `adhesion_rf.joblib` was fitted at 4.0 — a
  train/inference mismatch for no gain.
- **Decimating the large scales** (blur, subsample, filter small, upsample) was
  premised on large σ being expensive. It is not: 5.0 µm is 17.9 s of ~97 s, so
  making it *free* caps the saving at 18%. It also has a real failure mode —
  decimating without prefiltering aliases broadband content, and this channel is
  noise-dominated, which is the worst case for exactly that.

**Profiled 2026-08-25, and the answer is none of the above.** One chunk
(24 z-planes of 1024², GPU), each step timed with an explicit sync:

| step | 0.3 µm | 5.0 µm |
|---|---|---|
| one `gaussian_filter` | 4.4 ms | 206.7 ms |
| gather `g[iz,iy,ix]` | 2.6 ms | 2.7 ms |
| `.get()` 100 MB to host | 18.5 ms | 20.2 ms |
| `X[:, col] = ...` (stride 240 B) | 140.4 ms | 143.2 ms |

Per-chunk budget:

| | time | share |
|---|---|---|
| **`_eigh3x3` + (N,3,3) assembly** | **~52 s** | **~52%** |
| strided `X[:, col]` writes | ~8.4 s | 8% |
| all 19×6 `gaussian_filter` calls | ~8.1 s | 8% |
| gathers + `.get()` | ~1.3 s | 1% |
| unaccounted (allocation, page faults, GPU products) | ~30 s | 30% |

**`_eigh3x3` costs 3.77 s per call and is called twice per scale** (Hessian and
structure tensor), on a 0.91 GB `(N,3,3)` float32 array assembled by 9 strided
host writes, then computed in float64 on the CPU. **It does not depend on sigma
at all, which is exactly why the per-scale timings are flat** — that flatness was
the clue, and it is what makes taps, slab size and padding all red herrings.

**The pass-sharing idea this entry originally recommended is therefore demoted.**
Sharing the separable 1D passes attacks the `gaussian_filter` line — 8% of
runtime — so even a perfect 2× there is ~4% overall. Recording that plainly
because proposing it was the same error the clover gate made: optimising the
part that is not the cost. It remains correct, just nearly worthless.

**What to do instead, in order:**

1. **Move `_eigh3x3` to the GPU, and drop the `(N,3,3)` assembly.** The six
   components are already on the device when they are transferred off; computing
   the eigenvalues there and returning only the 3 results removes the host math,
   6 transfers per call, and the 0.91 GB temporary. Elementwise and
   embarrassingly parallel. **~2×, exact, no retraining.** Note it is
   memory-bound rather than FLOP-bound (~4 GB of traffic per call), so float64
   costs little even on a 3090's 1/64-rate FP64 units — prefer keeping float64
   over risking cancellation in the `q`/`p`/`r` terms in float32.
2. **Write X in blocks, not columns.** 10 single-column writes cost 9.2 s/chunk;
   one 10-column block costs 2.9 s. ~6 s/chunk, trivial change.
3. **Allocate X once and reuse it** across chunks rather than `np.empty` per
   chunk — 6.04 GB of first-touch page faults per chunk, repeated 13 times.
4. **Gather `gzz`/`gyy`/`gxx` once.** They are gathered for the LoG sum and then
   again individually for the Hessian, and the LoG builds a full slab-sized
   `(gzz+gyy+gxx)` temporary (629 MB at 5 µm) that is only ever read at the
   gathered indices.

**Dropping the 3.5/5.0 µm scales is re-scored by this profile too.** Its speed
argument was the filters, which are only 8% — so it is ~1.3×, not 2–3×. Its real
remaining value is that `max_pad` 63 → 21 collapses the ~6.5 GB device-memory
floor, which is what constrains 12 GB cards ([P43](#p43)). That is a memory
argument, not a speed one, and it still needs a retrain.

**Related:** [P43](#p43), [P03](#p03).

**History**
- 2026-08-25 — opened from the per-scale timing above, replacing two rejected
  proposals with a measured one.

---

<a id="p44"></a>
## P44 — the RF's training annotations were never saved; the model cannot be retrained

| | |
|---|---|
| **Status** | open |
| **Severity** | **S2** — no exported number is wrong, but the detector behind every adhesion number is unreproducible |
| **Opened** | 2026-08-25 |
| **Last verified** | 2026-08-25 |
| **Resolved** | — |
| **Where** | `pipeline/gui/adhesion_rf.py` — the "Adhesion labels" layer, `_train` |
| **Affects** | [P03](#p03) (blocks its live hypothesis *and* its fix), [H02](#h02), every adhesion column |

`_train` reads `ann_vol` from a live napari layer named **"Adhesion labels"**
(`adhesion_rf.py:101`) and fits on it. **Nothing ever writes that layer to
disk.** The tab persists `adhesion_proba.npy`, `clover_ruby_ratio.npy`,
`adhesion_instances.npz` and `adhesion_rf.joblib`; the annotations are not among
them. The `*_labels.npz` files in `_saved/` are cadherin and DAPI segmentation
output, unrelated.

So the hand-painted training set exists **only** as weights inside
`adhesion_rf.joblib`, and only for as long as that file keeps loading.

**Three consequences, in increasing order of how much they hurt:**

- **The model cannot be retrained or reproduced.** Any change that needs a
  refit — a different feature set, a different scale list, more examples —
  requires re-painting from scratch.
- **[H02](#h02) escalates.** That entry records an unpinned-scikit-learn pickle
  risk. With no training data it is not a *reproducibility* risk but a
  *permanent loss* risk: the file already warns on load (fitted 1.9.0, loaded
  1.8.0), and if a future version refuses it the detector is gone. Note the
  warning appears to be benign so far — a run under 1.8.0 reproduced the 1.9.0
  `proba >= 0.5` set to the voxel ([P43](#p43)) — but that is luck, not a
  guarantee.
- **It blocks [P03](#p03) from both ends.** P03's live hypothesis as of
  2026-08-25 is that the annotations were painted mostly on shallow, bright
  adhesions, so the RF never learned deep dim ones. **Testing that needs the
  annotations** — their depth distribution is the whole measurement. And if the
  hypothesis is right, the fix is retraining with deep examples, which needs
  them too. An S1 is currently gated on a file that was never written.

**The fix is small and should land before any retraining is attempted.** Persist
the annotation volume whenever "Train RF" runs — sparse (only nonzero voxels:
`np.argwhere(ann_vol > 0)` plus the class value), which for a few thousand
painted voxels is a trivial file. Write it per image into `_saved/` and
alongside the model, and load it back into the layer so a session can be
resumed. This is the same doctrine as `adhesion_proba.npy` being persisted
because "RF predict is the slow step" — hand-painting is far slower than RF
predict and is the one input that cannot be recomputed at all.

**Until then, treat `adhesion_rf.joblib` as irreplaceable**: back it up outside
the repo, and do not upgrade scikit-learn in either environment without
verifying it still loads and still reproduces a known `proba >= 0.5` set.

**Related:** [P03](#p03), [H02](#h02), [P43](#p43), [I09](#i09).

**History**
- 2026-08-25 — opened after asking whether the RF could simply be retrained to
  sidestep a feature-set change. It cannot.

---

## P45 — `batch_mesh` crashed on every invocation: two names used but never imported

| | |
|---|---|
| **Status** | **fixed** 2026-08-26 |
| **Severity** | S3 |
| **Opened** | 2026-08-26 |
| **Last verified** | 2026-08-26 |
| **Resolved** | 2026-08-26 |
| **Where** | `pipeline/batch_mesh.py` — the `from .vessel_mesh import (...)` block; both names are used inside `parameter_warnings` |
| **Affects** | every headless mesh rebuild |

`parameter_warnings` tests `method in (METHOD_ALPHA3D, METHOD_ALPHA_POISSON)`,
and the second name was never added to the import list. `parameter_warnings`
runs on **every** invocation, before any work, so `batch_mesh` raised

    NameError: name 'METHOD_ALPHA_POISSON' is not defined

unconditionally — `--dry-run`, `--check-only` and every real rebuild alike. The
same name is also referenced in the `--steps graph` guard's message, which would
have failed the same way one line later.

**And it was two names, not one.** With `METHOD_ALPHA_POISSON` imported the next
line failed the same way on `check_alpha_vs_spacing` — the very function CLAUDE.md
cites as an example of *not* letting the GUI and the headless path drift. It was
extracted into `vessel_mesh` for both callers to share and then not imported by
the second one, so the check ran in the tab and crashed the runner.

Introduced in `48b4f9d` ("Stop the cadherin GUI/headless drift; scope sliders to
the active method"), which added both uses. The fix is the two missing import
lines.

**Why nothing caught it.** There is no smoke test that runs `batch_mesh --dry-run`
— which is exactly what a `NameError` in the argument-resolution path would have
been caught by, in under a second, on any of the six images.

**What this says about the last week of work.** Every `vessel_graph.npz` and
`vessel_mesh.npz` written since `48b4f9d` came from the **napari tab**, not from
`batch_mesh` — the headless path could not have produced one. That is worth
holding next to [P26](#p26) and [P37](#p37): the artifacts in
`WT/..._saved/` carry no `vessel_recon.json` for the untagged build precisely
because nothing headless wrote them, so their parameters are unrecorded and a
rebuild does not reproduce them. Image 6's saved untagged graph is one edge of
455 µm; the same image rebuilt headless at defaults gives three edges and 433 µm.
**Neither is known to be right**, and the difference is not explained here.

**Related:** [P26](#p26), [P37](#p37), [P30](#p30).

**History**
- 2026-08-26 — found while running the first `--shear` build for
  `DIRECTION.md` §8 item 8a; fixed in the same change. `--dry-run` over all six
  images now reports `6 ok, 0 degraded, 0 failed`.

---

## P46 — a digitized wall is longer than the wall it stands for, so absolute τ reads low

| | |
|---|---|
| **Status** | open — **quantified, common-mode, not corrected** |
| **Severity** | S3 |
| **Opened** | 2026-08-26 |
| **Last verified** | 2026-08-26 (analytic disc and rectangle, `docs/scratch/shear/solver_check.py`) |
| **Resolved** | — |
| **Where** | `pipeline/shear.py` — `section_contour`, `wall_profile`; inherited by `lumen_graph.cross_section`'s `perim_um` |
| **Affects** | `perim_um`, `hydraulic_diam_um`, `tau_hat_*`, and any absolute τ in Pa. **Not** `tau_rel_1_um3` read as a within-segment profile |

The shear solve runs on a **pixelated** cross-section, and a pixelated boundary
is longer than the smooth boundary it approximates. For a digitized disc the
marching-squares perimeter is **5.6–7.9% high**, and — this is the part that
matters — **it does not converge**: 5.63% at 120 pixels per radius, essentially
the same as at 60. It is a digitization bias, not a discretization error, and
finer sampling of the same voxel grid cannot remove it.

The solve is faithful to the duct it is handed, so the extra wetted perimeter
comes out as τ̂ reading about **5% low** on a curved wall. On the analytic
rectangle, whose walls are axis-aligned and therefore not digitized at all,
`Q̂` lands within 0.08% and the wall profile within 1% — which is how we know
this is the boundary representation and not the solver.

**Why it is left uncorrected.** It is common-mode: it multiplies every section by
roughly the same factor, so it cancels in the within-segment τ ratio, and by
[DIRECTION](DIRECTION.md) §4.1 the within-segment profile is where the entire
spatial hypothesis lives. It does **not** cancel in an absolute Pa number, and
that is the number a perfusion rate would produce. A Crofton perimeter fixes the
disc (0.15% at fine resolution) but is 5.4% *low* on an axis-aligned rectangle,
so it trades one bias for another on the D-shaped sections this project actually
has; picking between them is not a decision to make by eye
([P26](#p26)).

**The real fix, if it is ever worth it,** is to stop handing the solve a binary
mask: `vessel_poisson` already produces a smooth indicator χ, and cutting the
section from χ with linear interpolation would give a sub-voxel wall with no
digitization bias at all. That changes `cross_section`'s input contract, and it
would have to change it for the area and the extents too, or the section the
shear sees and the section the geometry columns describe would be different
shapes — which is the drift this project's standing rule exists to prevent.

**The free check that catches everything else.** ∮τ dP = A·(dp/ds) is an exact
force balance, so with the normalisation used here it must equal 1. It is
recorded per section as `force_balance` and needs no ground truth. On image 6 it
runs 0.993–1.077 over 60 usable sections — i.e. the residual is this bias and
nothing larger is hiding under it.

**Related:** [P24](#p24), [P26](#p26).

**History**
- 2026-08-26 — opened when `solver_check.py` measured it. Two τ estimators were
  discarded on the way (contour-normal sampling, 36% high and converging the
  wrong way; first-layer closed form, ±28% scatter on a digitized curve); the
  near-wall quadratic fit that shipped is what left this as the only bias
  standing.
