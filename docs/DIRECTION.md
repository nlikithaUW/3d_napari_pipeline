# Direction — the scientific program and the state of the evidence

**Opened 2026-08-22.** This is the document to read first when picking the
project back up. It states what we are trying to measure, why, what is known,
what is not, and what to do next.

## How this relates to the other docs

| file | what it is | who owns detail |
|---|---|---|
| **this file** | the **program**: goal, hypothesis, physics, dataset, and the work queue | the direction and the ordering |
| `docs/ISSUES.md` | the **defect log**: one entry per problem, permanent IDs (`P01`…), full history | the detail and the dates for every `Pnn`/`Inn` |
| `docs/audit.md` | frozen evidence, 2026-08-10, commit `0faf2e5` | never edited |
| `README.md` | public description + short limitations table | links to `ISSUES.md` |

**The rule that keeps these from drifting:** this file cites issue IDs and quotes
only the single number that carries the argument. **The `Pnn` entry is always
authoritative** for the measurement, the date and the history. If a number here
disagrees with `ISSUES.md`, `ISSUES.md` is right and this file is stale. Do not
copy an entry's body here, and do not maintain a second status table — the "at a
glance" table lives in `ISSUES.md` and nowhere else.

New defects still get opened in `ISSUES.md` with a new ID. What lands here is
what changes the *program*.

---

# 1. The goal

Relate **per-adhesion tension in hiPSC cardiac fibroblasts** to **the local
mechanical environment of the vessel they sit on**.

The working hypothesis, as of 2026-08-22:

> Fibroblasts attach preferentially, and/or adhere with different force, at
> parts of the lumen that experience **high wall shear stress**.

**Fibroblasts are abluminal** — outside the lumen, **distributed through the
bulk collagen** rather than seeded onto the vessel (§3, corrected 2026-08-25).
They therefore feel **no shear directly.**
The causal chain has to be:

> flow → wall shear stress → the **endothelium** senses it → signals across the
> wall → the **fibroblast** responds (in where it attaches, and/or in how hard
> it pulls)

Two consequences that shape every measurement below:

1. The shear covariate for a fibroblast adhesion is the shear at the **nearest
   luminal surface point**, not at the adhesion's own position.
2. The **transmural distance** — how far the adhesion sits outside the wall — is
   a covariate in its own right, since any signal should attenuate with it.

**Geometry is a proxy, not the target.** Every geometric descriptor in this
project is only worth what it predicts about shear. That single sentence
reorders most of the previous plan; see §7.

# 2. Force is not in this dataset

FRET is **not calibrated** for these images — the `fret` block is absent from
`params.json` and `force_*_pN` exports blank. The uncalibrated
`ratio_clover_ruby` is the current tension proxy and it **does not work**
([P02](ISSUES.md#p02)): the numerator contributes ρ = +0.068 while
ruby-per-voxel contributes ρ = −0.773, because ruby sits at 6.8% contrast with
10.6% of its voxels clipped to exactly zero by `background_subtract`. The
exported "tension" is close to an inverse map of the ruby noise floor.

**So the mechanical half of the hypothesis cannot be tested on this dataset at
all.** What *can* be tested now is the **spatial** half — where fibroblasts are —
and that is the whole point of the current exploratory phase. See §8.

# 3. What we know about the experiment, and what we don't

| | |
|---|---|
| lumen geometry | **patterned/fabricated**, then endothelialised |
| perfusion | **yes, steady** |
| flow rate Q | **unknown now**, obtainable for future datasets |
| flow direction | **unknown now**, obtainable — and possibly recoverable from the existing data (§8.5) |
| fibroblast placement | **cast into the collagen before it gelled** — 0.5 × 10⁶ cells/mL in 6 mg/mL collagen. Distributed through the bulk, **abluminal**, mean spacing ~126 µm. **Corrected 2026-08-25**; previously recorded here as "seeded onto an already-formed, already-perfused vessel", which was wrong |
| culture duration | **7 days**, then fixed |
| flow during/after seeding | **question retired** — there was no separate seeding event; the fibroblasts predate the endothelium. See §9 for what replaces it |
| pattern design file | **unknown whether it exists — ask**; and it is a prior, not the lumen — §3.1 |

**Why "patterned" matters more than it sounds.**

- **Geometry is exogenous.** A fabricated mold fixes the geometry, so shear
  cannot have remodelled it and the cells cannot have remodelled it. In vivo,
  geometry and shear co-determine each other and the causality is circular. Here
  it runs one way. That is a genuinely strong design position and belongs in the
  eventual methods section.
- **There may be a ground-truth file.** A mask/CAD/mold spec would give intended
  caliber, branch angles, junction positions, inlet vs outlet — and therefore
  **flow direction**, **bifurcation vs confluence**, and the **flow split** at
  each Y. It would also validate the reconstruction against something external,
  which is what [I03](ISSUES.md#i03) has wanted all along and has so far only
  approximated with synthetic phantoms.
- **It probably fixes the cross-section shape.** See §4.

**What the corrected fibroblast placement changes.** "Abluminal" survives intact
and so does the whole causal chain in §1 — the endothelium is still the sensor
and the signal still crosses the wall. Three things do change:

- **Transmural distance is a wide covariate, not an offset.** A layer seeded on
  the outer wall would sit at a roughly constant small distance. Cells cast
  through the gel sit at *every* distance out to the edge of the field. That is
  the difference between a nuisance parameter and a primary axis — see §7.5,
  where it also supplies an in-image control.
- **Their adhesions are mostly to collagen, not to the endothelium.** The
  tension proxy is therefore reading fibroblast-on-matrix grip for most cells.
- **They compact the matrix.** Contracting collagen is what fibroblasts do, and
  6 mg/mL is softer than the 7.5 mg/mL of the reference protocol. §3.1's
  design-vs-actual gap is correspondingly larger here than in a
  fibroblast-free system.

## 3.1 The design file would be a prior, not the lumen

*Added 2026-08-25, from `docs/papers/ncomms8858.pdf` (Zheng, Chen & Lopez, Nat
Commun 6:7858, 2015) — the lab's own description of this microvessel system.*

§3 treats the pattern design file as an external validation target. It is worth
one, but it is **not** the lumen. The mould defines the *collagen* channel; the
lumen is the space inside a confluent endothelium, and the two differ
systematically rather than randomly:

- **The monolayer is subtracted from both walls.** Endothelium a few µm thick,
  nuclei bulging further, so the true lumen is **always smaller than the
  design** by roughly a constant absolute amount — single-digit percent at our
  ~90 µm heights, worse on image 6. The paper's only acknowledgement of any
  mismatch is exactly this scale-dependence: strand deviations attributed to
  "surface irregularities on the endothelial monolayer", "more prominent in
  vessels <200 µm in diameter". It never compares an imaged lumen to a design.
- **Corners round off.** A planar stamp gives a rectangular groove capped by a
  flat second collagen layer, and cells do not sit in sharp corners. So the
  low-shear corners a rigid rectangular-duct solve predicts are the first thing
  that is not real — which matters, because "low shear, therefore no signal" is
  read off precisely there.
- **There is a seam** at the mid-plane where the two collagen layers are sealed.
  No design geometry contains it.
- **Cells compact the gel**, and **the gel is compliant** — the paper drives
  10–1,000 Pa across collagen and solves on rigid walls. **This is worse in our
  system than in the reference one**: ours contains a contractile fibroblast
  population cast through the bulk (§3) and uses 6 mg/mL collagen against their
  7.5, so there is more pulling on a softer matrix. Our 7 days is shorter than
  their 1–2 weeks, which cuts the other way.

**Our own six images argue both ways.** The perpendicular cross-sections (§8C.3)
give aspect **2.02 / 1.91 / 1.89 / 1.80** on images 1, 2, 3 and 5 — four images
within ±6%, which is more reproducibility than remodelling would leave, and an
argument that the geometry *is* design-dominated. Against that, §8A.3's heights
span **49–97 µm** and a lithographic process should not vary in height. One of
those two readings is wrong and we do not know which, which is the sharpest
current reason to ask for the file (§9 item 2) — not to treat it as the answer.

**Which is the weaker link splits by image.** On images 1–3 (6–20% of sections
clipped, image 1 agreeing with [P24](ISSUES.md#p24)(b)'s independent direct cut)
the reconstruction is probably tighter than the design-vs-actual gap: solve on
the measured lumen. On images 4 and 5 (44–56% clipped, widths are lower bounds)
and image 6 (inscribed radius moved 60% between methods, §8C.3) the
reconstruction is the larger error, and a design geometry would be worth more
than a measurement.

# 4. The physics, and the one number that decides the formula

For a **circular** tube: τ ≈ 4μQ / πR³ — shear goes as caliber **cubed**.

For a **wide flat duct** of fixed height h and width w ≫ h: τ ≈ 6μQ / (w·h²) —
shear goes as **1/w**, and only as 1/h².

Lithographically patterned channels are typically fixed-height, varying-width,
which would make the second formula the right one. The measured cross-sections
are consistent with that: ~150–180 µm in plane by ~85–92 µm deep, aspect
1.64–2.15 on the reference image and ≈3 on `V3b_..._04step_002`
([P24](ISSUES.md#p24)'s 2026-08-13 correction).

**This is not cosmetic.** On image 6 (wide → thin waist → wide), R⁻³ predicts a
~27× shear rise at the waist; 1/w predicts ~3×. Both are usable ranges; they are
an order of magnitude apart.

**The check that settles it** (Tier 1, §8): measure lumen **height and width
separately** across all six. If height is flat at ~85–90 µm while width varies
freely, it is a lithographic channel and the flat-duct formula applies.

Regardless of which: **keep width and height as separate quantities.** They enter
the shear formula with different exponents, and collapsing them into one
"caliber" destroys the distinction.

## 4.1 Relative shear needs no flow rate

Along an unbranched segment, Q is constant by conservation. So the **shape of the
shear profile along a vessel** is fully determined by the caliber profile, with
no knowledge of Q at all. Absolute τ and cross-image comparison need Q; the
within-image gradient does not. **Every result available today is a within-image
relative-shear result**, and that is enough for a spatial hypothesis.

## 4.2 At fixed Δp, wall shear is linear in hydraulic diameter

§4's two formulas are both **fixed-Q** statements, and fixed Q is the regime
where a caliber error is punishing. This system is not run at fixed Q. It is run
at fixed pressure — a gravity head between reservoirs, and the matching COMSOL
boundary conditions in `docs/papers/ncomms8858.pdf` are a constant-pressure inlet
(10–1,000 Pa) and a zero-pressure outlet.

For fully developed flow in a uniform duct, a force balance on the control volume
gives τ_w·P·L = Δp·A, hence

    τ_w = (Δp / L) · D_h / 4,    D_h = 4A / P

**Wall shear is then linear in hydraulic diameter, and independent of Q and of
viscosity entirely.** A 7% shrink of the lumen from the monolayer costs 7% in
shear, not 23%. §7.3's warning that R⁻³ amplifies a caliber error is a fixed-Q
warning and does not transfer to this regime.

**The junction is where it does transfer.** In a network the per-segment Δp/L is
not given — it is set by the resistance of every other segment, and
rectangular-duct resistance goes as 1/(w·h³). So for a Y:

- **absolute shear along one segment** is forgiving of caliber error (linear),
  but
- **the flow split between the two daughter branches** is not: it rides on h³.

§4.1's "relative shear needs no flow rate" holds along an *unbranched* segment by
conservation, and does **not** rescue the split at the crotch. The crotch is the
structure this program cares about most, so that is where lumen accuracy has to
be earned — see §3.1 for how far the design geometry can be trusted to supply it.

Both results assume fully developed flow in a locally uniform duct. Turns, tapers
and the junction itself violate that, which is why the paper solves numerically
rather than algebraically.

## 4.3 How to actually compute τ_w on this reconstruction

*Added 2026-08-25.* §4 and §4.2 give the algebra for idealised ducts. This
section says what to compute from the geometry we actually have.

**The exact local statement.** For fully developed flow down a duct of arbitrary
cross-section, the axial velocity `u` on the section solves

    µ ∇²u = dp/ds,   u = 0 on the wall

and the wall shear at each boundary point is `τ_w = µ|∇u|`. Everything in §4 and
§4.2 is a special case of that one equation.

**Solve it; do not look up a formula.** Two reasons, both material here:

- **No shape correlation fits this data.** Every closed form in §4 assumes a
  circle, a rectangle or parallel plates. `cross_section`'s own docstring
  records that the measured sections are **D-shaped, not elliptical**
  ([P24](ISSUES.md#p24)(b)) — which is also why it reports extents rather than
  fitted semi-axes.
- **It gives the circumferential distribution, and no closed form does.** A
  single τ per section is an average around the perimeter, but an adhesion sits
  at *one* circumferential position, and on an aspect-2 section the flat sides
  and the ends differ severalfold. §7.4 already requires circumferential
  position be preserved; collapsing to a per-section scalar throws away exactly
  what the per-adhesion join needs.

It is a 2D Poisson solve on a binary mask of a few thousand pixels, once per
section. Cheap.

**The normalisation — and why the result needs no unknowns.** Solve once per
section with `dp/ds = -1`, `µ = 1`, giving `û` and `Q̂ = ∫û dA`. Along an
**unbranched** segment Q is conserved (§4.1), so the true gradient at station
`s` is `(dp/ds)(s) = -Q/Q̂(s)`, and

    τ_w(s, θ) = |∇û(s, θ)| · Q / Q̂(s)

**Q is a single unknown multiplying the entire profile.** So the *shape* of the
shear map along a segment — which is the whole of the spatial hypothesis — is
fixed by geometry alone, with no flow rate, no viscosity and no pressure.
Normalised within a segment it is a measurement, not an assumption. This is
§4.1's point made local and per-boundary-point.

**What this needs that does not exist yet.** `cross_section` returns
`area_um2`, `major_um`, `minor_um`, `aspect`, `clipped`, `touches_domain`, and
builds the binary section internally and then discards it. Two additions:

- **return the section mask** (or its contour) — that is what the solve runs on;
- **perimeter**, which `D_h = 4A/P` needs before *any* of §4.2's algebra can be
  evaluated numerically. It is not currently computed anywhere.

> **Built, 2026-08-26.** `pipeline/shear.py`, plus both additions above:
> `cross_section` now reports `perim_um` and fills an optional `extras=`
> out-dict with the section mask and its frame. `build_graph(compute_shear=True)`
> — `batch_mesh --shear`, and a checkbox on the Vessel Mesh tab so the two
> cannot drift — adds `hydraulic_diam_um`, `q_hat_um4`, `tau_hat_mean/min/max_um`,
> `tau_rel_1_um3`, `force_balance` and `wall_ok_frac` per station. Off by
> default, because `batch.py` runs at defaults and a default is a decision.
>
> **Validated against closed form, not against itself** ([I03](ISSUES.md#i03)'s
> doctrine): `docs/scratch/shear/solver_check.py` solves the disc analytically
> and the rectangle by its double sine series. `Q̂` lands within 0.08% and the
> wall profile within 1% on the rectangle. The disc is 5% low for a reason that
> is not the solver — [P46](ISSUES.md#p46), a digitized wall is longer than the
> wall it stands for, and it is common-mode.
>
> **The circumferential variation is not a detail.** On an aspect-1.8 rectangle
> the solved τ varies **9× around the perimeter**, and on image 6's real
> sections the median is 6.7×. §7.4 asked that circumferential position be
> preserved; this is the size of what a per-section scalar would have thrown
> away.
>
> **Visible on the surface, 2026-08-26.** `build_graph(..., extras=)` keeps the
> full per-boundary-point profile — every wall point with its 3D position and
> its τ — as `extras["wall"]`, written to `vessel_shear.npz` by
> `batch_mesh --shear` and drawn by the Vessel Mesh tab as a per-vertex surface
> layer beside `Vessel curvature`. On image 6 that is 27,682 wall points
> covering 88% of the mesh vertices (70% after dropping clipped and censored
> sections), and the waist reads red against dark arms — the arc-length plot,
> in place.
>
> **The view is a view.** Vertices take the value of the **nearest measured**
> wall point, never an interpolation, and wall with no measurement within a
> section spacing is **not drawn** rather than filled — `curvature` fills its
> NaNs with the median, which is right for a field defined on every vertex and
> would be a fabricated measurement here. Where the axial sampling shows, the
> fix is a smaller `section_every_um`, which is a real measurement and costs
> seconds.
>
> **One caution that is about statistics, not about the picture.** Cutting
> sections closer together for a better view also multiplies the rows. Adjacent
> sections on a 60 µm-wide vessel are nearly the same cut, so a finer
> `section_every_um` does **not** buy independent measurements — the
> correlation length is the caliber, not the spacing. Anything that counts
> sections, or fits with an n, has to use the recorded `section_every_um` rather
> than the row count.
>
> **And it is the same array the adhesion join will use.** The wall point cloud
> is what an adhesion joins to — nearest luminal wall point, read its τ — so
> what is on screen is literally what gets joined. That join must carry the
> **distance to that wall point** as a column: the fibroblasts are cast through
> the bulk gel at a range of distances (§8D), so some adhesions are nowhere near
> the lumen and must not silently receive a τ.

Take the perimeter from a marching-squares contour, not a count of boundary
pixels: the staircase estimate overestimates a smooth boundary by up to 4/π.

**Where it stops being valid — flagged per section, not argued away:**

- **`clipped` sections are unusable for this, more so than for area.** A clipped
  boundary contributes a straight edge that is not a wall, and the solve would
  impose no-slip on it. Images 4 and 5 clip 44–56% of their sections (§8C.3), so
  they are out until [P32](ISSUES.md#p32) is fixed.
- **Junctions are excluded, not approximated.** Q is conserved only along an
  unbranched segment; the split at a Y rides on h³ (§4.2) and is not recoverable
  from geometry alone.
- **Fully developed flow fails at turns, tapers and the crotch** — §4.2's
  closing caveat, and the reason the reference paper solves numerically. Treating
  each section as locally uniform is a lubrication approximation whose error
  grows with `dA/ds`, which is largest exactly at image 6's waist, i.e. at the
  place the first result wants to read.

**Recommended first result.** Image 6, one unbranched segment, τ along arc
length, normalised within the segment. It is z-contained, clips only 5% of
sections, and carries the largest caliber gradient in the set (§5).

# 5. The dataset

Six images, `3d_images_folder/images/WT/`. Voxel 0.325 × 0.325 × 0.4 µm except
`09step` which is **0.9 µm in z**.

| # | file | structure | role in the program |
|---|---|---|---|
| 1 | `V3_..._04step_002` | Y branch | the calibrated reference — every number in `ISSUES.md` is from this image |
| 2 | `V3_..._04step_003` | tube in view, **junctions off both ends, flaring** | flare without a visible crotch; junction distance is **censored**, not missing |
| 3 | `V3_..._09step_001` | Y, **asymmetric** — large caliber one end, wide crotch branch | topology generalisation — **and the 0.9 µm z-step file** |
| 4 | `V3b_..._04step_001` | Y, near-equilateral, equal-caliber branches | symmetric crotch |
| 5 | `V3b_..._04step_002` | straight tube, constant ellipse, **aspect ≈ 3** | **negative control** — K ≈ 0 everywhere is the known right answer |
| 6 | `V3b_..._04step_003` | curved; wide → **very thin waist** → widens, still small | **best shear image** — large gradient inside one field of view |

## 5.1 What the inventory implies

- **Three crotches, and no two alike.** Good for generalisation, bad for
  statistics: adhesions within one crotch share that crotch's geometry,
  illumination, depth and segmentation quality. **Effective n for "crotch" is 3,
  clustered** — a future force analysis must cluster by image or it will report a
  confidently wrong p-value off what looks like thousands of samples.
- **Under a shear framing the predictor varies within all six images**, including
  the three non-branching ones. That is a large statistical improvement over the
  crotch-vs-tube framing and is the main argument for the reframe.
- **Image 3 confounds crotch topology with z resolution.** It is one of only
  three crotches and it is the only 0.9 µm-in-z image; after `to_isotropic` its
  z detail is interpolated, not measured, and every quantity here is z-sensitive.
  Three images cannot design this away. **Measure it and state it.**
- **Image 5 is a free hard negative control.** A constant cross-section is a
  generalised cylinder, so Gaussian curvature K ≈ 0 everywhere is correct. Any
  structure a descriptor shows on image 5 is reconstruction roughness.
  It is unforgiving in a second way: at aspect 3, κ around the circumference
  varies by a³/b³ = **27×** between the flat sides and the sharp ends, on a
  vessel that is structurally uniform end to end. **Nothing should be built on
  |κ|.**
- **No image appears to contain a real vessel terminus.** The tubes leave through
  the x/y field-of-view walls and image 6's thin section is a waist, not an end.
  So "tip" is **not** a class to detect; what the pipeline must do instead is
  **recognise a crop** (surface meeting an x/y domain face) and mark that region
  censored. `vessel_poisson`'s Neumann boundary is already doing the right thing
  physically — it lets the vessel leave without fabricating a cap.

# 6. What is carried forward from `ISSUES.md`

Only the live, load-bearing findings. IDs link to the authoritative entry.

## 6.1 Numbers that cannot be trusted

| column / claim | why | entry |
|---|---|---|
| `ratio_clover_ruby` | tracks the ruby noise floor, not tension (§2) | [P02](ISSUES.md#p02) |
| `force_total_pN` | units are pN·µm³; blank now, wrong the moment FRET lands | [P10](ISSUES.md#p10) |
| `n_fibroblasts`, `nuc_is_fibroblast` | classifier unvalidated; 38%→19% on one slider step; 5.0–14.9% across the six | [P05](ISSUES.md#p05) |
| every column, as a *sample* | adhesion detection is depth-truncated — far wall yields 37% of near-wall density, confirmed 6/6 | [P03](ISSUES.md#p03) |
| `curv_*` | no scale-independent value on the alpha mesh: implied R climbs 14.2 → 134.7 µm as `radius_um` goes 5 → 40 | [P24](ISSUES.md#p24) |
| ~~`path_curv_xy`~~, ~~`dist_branch_um`~~ | **gone** — both columns were deleted with the 2D centerline and redefined on the graph (§8C). See the two rows below | [P08](ISSUES.md#p08), [P29](ISSUES.md#p29) |
| `path_curv_*`, `dist_junction_um`, `caliber_um`, `lumen_*` | the replacements, and they are only as good as the lumen: on images 4 and 5 it leaks through the z faces, so the junction count is wrong and 44–56% of cross-sections are clipped. Right on the other four | [P32](ISSUES.md#p32) |
| every geometry column, today | no untagged `vessel_graph.npz` exists, because adopting alpha→poisson as the shipped reconstruction is still an open decision (§8C.7) — so they export **blank** rather than mixing two reconstructions | [I03](ISSUES.md#i03) |
| depth relationships | ρ(ratio, z) ranges −0.241 to +0.307 across the six — **the sign flips** | [P06](ISSUES.md#p06) |

## 6.2 The reconstruction, as it stands

`alpha 3D solid -> poisson` ([I06](ISSUES.md#i06), `RESULTS.md` §12) is the
direction of record. Alpha bridges the wall gaps in 3D; its *solid* interior is
rasterised and handed to a Neumann DCT Poisson solve, which fits a smooth
indicator function whose level set is the surface.

| | alpha 3D (old default) | poisson on mask | **alpha → raster → poisson** |
|---|---|---|---|
| `slope_loglog` (0 = plateau) | −0.214 | −0.018 | **−0.024** |
| non-manifold faces | 8.3% | 0.0% | **0.0%** |
| pieces | 42 | 1 | **1** |
| coverage ≤ 5 µm | 99.1% | 87.5% | **91.3%** |
| coverage ≤ 10 µm | — | 93.3% | **98.2%** |
| area > 5 µm from data | 19.2% | 19.6% | **19.9%** |

**This is the first reconstruction in the project to reach [P24](ISSUES.md#p24)'s
plateau.** Read `support` by **area, not by max** — the max (46.47 vs 26.70 µm)
looks like a doubling and isn't; the two are indistinguishable through p90.

Alpha 3D remains the shipped default. `mesh_score` returns four numbers and never
combines them; **choosing** is what is still open in [I03](ISSUES.md#i03).

## 6.3 Geometry facts, established and reusable

- **Caliber ~150–180 µm in plane × ~85–92 µm deep**, aspect 1.64–2.15 on the
  reference, ≈3 on image 5. Confirmed by three methods sharing no machinery.
  ([P24](ISSUES.md#p24))
- **The vessel is contained in z** *on the reference*: 0.01% of the mask in the
  first z slice, zero in the last, ~7 µm of headroom. This is what
  `vessel_poisson.lumen_core`'s z-enclosure discriminator rests on.
  ([P24](ISSUES.md#p24)) **It does not hold on all six** — image 4 puts 1.224%
  of its mask in the last slice, and on images 4 and 5 the fitted lumen leaves
  mainly through the z faces ([P32](ISSUES.md#p32)).
- **Inscribed lumen radius, median per image:** 40.6 / 42.0 / 46.0 / 28.0 / 35.2
  / 30.1 µm, from `lumen_graph`'s 3D distance transform of the lumen (§8C.3).
  **It is the *minor* semi-axis** (an inscribed sphere is capped by the shortest
  direction), so it is not the caliber shear depends on — see §7.3. The earlier
  38.1 µm figure came from `lumen_center_z`, which is deleted; the two agree on
  the reference to within 7%.
- **Cross-section, measured perpendicular to the axis** (§8C.3, the first time
  this was measurable at all): ~13,600 µm², 169–186 µm by 90–100 µm, aspect
  ~1.9–2.0 on images 1–3. Agrees with [P24](ISSUES.md#p24)(b)'s direct cut by a
  route sharing no machinery. Images 4 and 5 clip 44–56% of their sections and
  are lower bounds.
- **The wall's "gaps" mostly are not holes.** Mask-gap directions show the same
  peak radius as mask-present ones at ~40% the brightness — same wall, not kept.
  The global floor plus `min_size` delete **59–85% of the dim wall**.
  ([P28](ISSUES.md#p28), [P19](ISSUES.md#p19))
- **A cross-section concavity is real** and visible in the raw channel but not in
  the mask. ([P24](ISSUES.md#p24), [P28](ISSUES.md#p28))

## 6.4 Structural obstacles this program will hit

**The first three rows are gone as of §8C** — they were symptoms of the 2D
projection, and the replacement removed them rather than fixing them. They are
struck through rather than deleted so the shape of the old problem stays legible
next to what replaced it.

| obstacle | effect on this program | entry |
|---|---|---|
| ~~centerline is a **point set in raster order** — no arc length~~ | **resolved.** The graph carries arc length by construction; dr/ds is unblocked | [P27](ISSUES.md#p27) *fixed* |
| ~~branch nodes come from the **xy silhouette**~~ | **resolved.** Junctions come from a 3D skeleton; image 1 now reports its one real Y | [P29](ISSUES.md#p29) *superseded* |
| ~~`open_um`/`prune_um` **amputate tapering vessels**~~ | **resolved.** Pruning is now a ratio of the local caliber, so a taper is not cut off by a constant | [P23](ISSUES.md#p23) *superseded* |
| **the fitted lumen leaks through the z faces** on images 4 and 5 | the graph invents junctions there and clips half its cross-sections; 2 of 6 images unusable for caliber | [P32](ISSUES.md#p32) |
| `mesh_curvature` has **no boundary gate** | 7.07% of vertices at `radius_um` 5.0 are fitted to half-neighbourhoods, and they are the vessel ends | [P31](ISSUES.md#p31) |
| `subsample_um` sets mesh resolution silently — 169:1, 2.25 µm spacing | distances and curvature are quantised at ~2 µm | [P20](ISSUES.md#p20), [P09](ISSUES.md#p09)(b) |
| ~~`lumen_core` requires the lumen to reach 4·σ~~ | **checked, §8A.3: it refuses on none of the six.** Image 6 clears it by the least margin | `vessel_poisson.lumen_core` |
| ~~batch built the five non-reference meshes at **widget defaults**~~ | **resolved §8C.1** — all six rebuilt at one recorded parameter set. Comparable; still not validated | [P26](ISSUES.md#p26) |

---

# 7. New findings, 2026-08-22

These came out of the design discussion and are not in `ISSUES.md`. The first two
should be opened as entries.

## 7.1 Gaussian curvature is a caliber-acceleration detector, not a branch detector

For any tube-like surface with cross-sectional radius r(s) along its axis:

> **K = −r″ / ( r · (1 + r′²)² )**

K depends on **r″ and nothing else.** Therefore:

| structure | r″ | K |
|---|---|---|
| straight tube | 0 | **0** |
| **pure taper / cone** | 0 | **0 — indistinguishable from a straight tube** |
| flare | > 0 | **< 0** |
| waist | > 0 | **< 0** |
| crotch | — | **< 0** |

**So K flags the crotch, the flare and the waist identically, and is blind to a
taper.** It cannot deliver the crotch-vs-flare-vs-tube label the program wants,
and it is blind to the taper distinction that matters most for shear.

Two consequences:

- **Crotch vs flare is topological, not differential** — at a crotch the lumen
  cross-section *splits*, at a flare it does not. No local surface quantity sees
  that. It requires a real 3D centerline **graph**.
- The per-vertex quadric curvature machinery
  (`cell_mesh.mesh_curvature`) is **not on the critical path** for this program.
  That is a reversal of where this discussion started.

Also worth recording: you cannot reach a coarse curvature simply by raising
`radius_um`. The quadric assumes the patch is a **height field over the tangent
plane**, which fails once the window approaches the caliber; and the
neighbourhood is a **Euclidean** ball, so at a Y crotch it scoops in the other
branch's wall across the gap. The one location of interest is where the
neighbourhood definition is most wrong.

## 7.2 `measure_core`'s `dist_mesh_um` is unsigned — and that likely explains a failed experiment

`pipeline/measure_core.py:147` builds a `KDTree` over surface vertices and calls
`.query()`, which returns an **unsigned** distance.
`cell_mesh.signed_distance_to_mesh_vertices` exists and is used in `rf.py`,
`gui/cadherin.py`, `gui/cell_mesh.py`, `gui/clover.py` and `gui/adhesion_rf.py`
— **but not in `measure_core`.**

An attempt was made to separate the two cell types by nucleus distance to the
reconstructed mesh, expecting bimodality (endothelial nuclei in/near the wall,
fibroblast nuclei outside). It came out unimodal. **An unsigned distance folds
the inside and outside populations onto each other by construction** — a nucleus
6 µm inside and one 6 µm outside both read `6.0`. If the test used this path, it
could not have succeeded.

**Even with the sign fixed, absence of bimodality would be weak evidence.**
Expected separation is perhaps 5–10 µm; the error budget against it:

- nearest-**vertex** rather than nearest-triangle — ~2.25 µm spacing, and the
  estimator is **biased outward**, not merely noisy ([P09](ISSUES.md#p09)(b));
- a uniform translation of the whole surface is **invisible to every check that
  exists** — `vessel_poisson.reconstruct` carries an explicit
  `(k−1)·voxel/2 − h/2` correction precisely because omitting it shifts the
  surface ~1 µm and nothing in `mesh_score` catches it;
- ~19.6% of surface area sits >5 µm from any data, and the mask is missing
  59–85% of the dim wall ([P28](ISSUES.md#p28)), so the surface is locally
  displaced wherever the stain was weak;
- `measure_core` subsamples the nucleus volume `::3` before the query.

Two Gaussians 8 µm apart with 5 µm spread each are one broad unimodal peak.

**The right test is not a bimodality hunt.** VE-cadherin rings are definitionally
endothelial: use ring-adjacent nuclei to measure the endothelial distribution's
location and spread **empirically**, then ask whether the remaining nuclei are
displaced outward relative to that calibrated reference. That is a one-sided test
against a measured control and it works even under heavy overlap.

**If side-of-wall separates the two cell types, it is a physical solution to
[P05](ISSUES.md#p05)** — independent of marker intensity, and therefore
cross-checkable against the current classifier rather than sharing its failure
modes. That would be a bigger result than the geometry work.

## 7.3 The inscribed radius is not caliber, and R⁻³ amplifies the error

`lumen_center_z`'s `radius_um` returns the **minor** semi-axis. On image 5 at
aspect 3 that understates the vessel threefold — and under R⁻³ that is a **27×**
error in estimated shear. Caliber for shear purposes must be the cross-sectional
**area** (or width and height separately, §4) on a plane perpendicular to the
local axis. Not optional.

> **Acted on, 2026-08-22.** `lumen_center_z` is deleted with the rest of the 2D
> centerline. `lumen_graph` still returns an inscribed radius — `caliber_um`,
> from a 3D distance transform — and it is **still the minor semi-axis**, for
> the same geometric reason (an inscribed sphere is capped by the shortest
> direction). What changed is that it is no longer the only caliber available:
> `lumen_area_um2`, `lumen_major_um` and `lumen_minor_um` come from a real
> perpendicular cut (§8C.3). `caliber_um` is kept because it is the right thing
> to prune spurs against, and it is documented as the wrong thing to compute
> shear from.

> **Qualified, 2026-08-25.** The **27×** above is a *fixed-Q* number. This system
> is pressure-driven, where τ_w is linear in hydraulic diameter and the same
> threefold error costs 3×, not 27× (§4.2). The conclusion is unchanged — the
> minor semi-axis is still the wrong caliber and area or width×height is still
> required — but the amplification quoted here is an upper bound, and the regime
> where it actually bites is the **flow split at the crotch**, not shear along a
> segment.

## 7.4 The depth confound becomes testable

In a non-circular duct wall shear is **not uniform around the circumference** —
for an ellipse it is highest at the ends of the **minor** axis. Here the minor
axis is z. So the model predicts **high shear on the top and bottom walls, low
shear on the sides**:

- predicted shear is **U-shaped in z** — high at both extremes;
- optical attenuation is **monotonic in z**.

**Those functional forms are separable.** This is the strongest handle the
project has had on the depth confound that
[P03](ISSUES.md#p03)/[P06](ISSUES.md#p06) describe, and it comes for free from
the vessel being flat rather than round. **Design the export to preserve where on
the circumference an adhesion sits**, not merely its depth.

## 7.5 Shear and geometry are collinear — only a bench control breaks it

Shear is *computed from* geometry. So "cells accumulate where shear is high" and
"cells accumulate where the vessel is narrow" are **the same prediction**.
Nothing in a perfused-only dataset separates a shear-signalling effect from cells
settling differently in a narrow groove (matrix thickness, outer-surface
curvature, gravity, seeding hydrodynamics).

**The control is identical geometry, no flow**: same patterned chip,
endothelialised, never perfused, fibroblasts seeded the same way. If the spatial
preference disappears, it is shear. If it persists, it is geometry.

This is a bench request with a long lead time and it is the single most valuable
thing this exploratory phase can hand to the next experiment.

> **A partial control may already be in the data — added 2026-08-25.** The
> fibroblasts are cast through the bulk gel (§3), not layered onto the vessel, so
> the ones sitting far from the wall cannot be responding to anything the vessel
> is doing. That is a comparison group in the *same image*, with the same optics,
> staining and threshold — no new experiment. If a tension or density measure
> varies with distance from the vessel and then flattens beyond some range, the
> flat part is the local null.
>
> **It does not replace the no-flow control**, because distance-from-vessel is
> itself confounded with depth ([P03](ISSUES.md#p03)) and with the collagen
> compaction gradient a contractile population creates around a stiff inclusion.
> It is cheaper and available now, which the bench control is not.
>
> **It is also contingent on there being cells out there to compare against.**
> §8D finds almost no nuclei more than 5 µm abluminal, which either means the
> fibroblasts have migrated onto the vessel or means there are very few of them.
> Settle §8D's count question first — this control is only real if the far
> population exists.

## 7.6 `lumen_center_z`'s inscribed radius is gated on the filled footprint

`centerline.lumen_center_z` rebuilds the same
`binary_fill_holes(binary_closing(...))` footprint the 2D centerline uses and
sets ``valid = foot[cy, cx]``. That fill closes the region enclosed by a Y, so a
centerline point sitting **between the arms** is marked valid and its
``radius_um`` — the max of the distance-from-wall field down that column — is a
clearance measured from *outside* the vessel, not a lumen caliber.

Same root cause as §8A.2's width failure. Audited in §8A.2A; the effect on the
reported medians turns out to be small on five of six images, but the mechanism
is real and it is one more reason the 2D centerline is being replaced (§8B).

> **Resolved 2026-08-22.** The function is deleted. `lumen_graph`'s caliber is
> the distance transform of the **lumen volume itself**, with no footprint, no
> fill and therefore no gate — a point between the arms of a Y is simply not in
> the lumen, so the failure mode cannot arise. Image 6, which §8A.2A left
> explicitly unsettled, now reads 30.1 µm against the old 18.8; §8C.3 flags that
> 60% change as the first thing to check by eye.

---

# 8. The queue

## Tier 1 — **done 2026-08-22**. Results in §8A; scripts in `docs/scratch/tier1/`.

1. ~~Signed vs unsigned in the failed nuclei test~~ — **answered, and it was the
   problem.** §8A.1.
2. ~~Lumen height vs width across all six~~ — height measured; **width could not
   be measured reliably by any available route.** §8A.2.
3. ~~Inscribed-radius profile + Poisson go/no-go~~ — Poisson runs on all six;
   two new defects found. §8A.3.
4. ~~`dist_mesh_um` distribution~~ — measured; uninformative for the abluminal
   question because it is unsigned. §8A.4.

## Tier 2 — **items 5, 6a and 7 done 2026-08-22**. Results in §8C.

5. ~~Redo the cell-type separation properly~~ — **done, and it replicates on 5 of
   6 images** (AUC 0.644–0.773, best p 6.9e-07); unsigned is chance, confirming
   the diagnosis. §8C.4. **But §7.2's endothelial calibration failed
   structurally** (§8C.5), so this is still an agreement between two unvalidated
   measures, not a validated one. See 5a.
5a. **Find an independent endothelial label.** New, and it is now the blocker on
   [P05](ISSUES.md#p05): the geometric measure cannot be validated against a
   reference defined by position relative to the same wall it measures against.
   [P14](ISSUES.md#p14)'s unused `cadherin_junction_proba` is the nearest
   candidate in the repo; [I02](ISSUES.md#i02) is the entry.
6. **VE-cadherin as an internal shear reporter.** Endothelial elongation and
   junction alignment under sustained flow give (a) independent validation of the
   geometric shear map and (b) **flow direction, now, without knowing Q**. Uses a
   channel already imaged and already segmented. *Not started.*
6a. ~~`batch_mesh` run: alpha→poisson meshes for the other five images.~~
   **Done**, all six at one recorded parameter set, image 1 byte-identical.
   §8C.1.
7. ~~The centerline replacement (§8B)~~ — **done.** `pipeline/lumen_graph.py`;
   `pipeline/centerline.py` deleted. In-plane caliber is measurable at last, so
   the shear estimate is unblocked on four of six images. §8C.2, §8C.3.
7a. **Fix or bound [P32](ISSUES.md#p32)** — the lumen leaks through the z faces
   on images 4 and 5, which is what breaks their topology and clips their
   cross-sections. Until then those two images carry no usable caliber. The
   candidate levers are in P32; none is a slider session, because a per-image
   parameter has to be recorded rather than chosen by eye
   ([P26](ISSUES.md#p26)).
8. **Fibroblast density vs relative shear** — the actual exploratory result
   available today, and it is now reachable: image 6 has a measured
   cross-section profile and is unbranched, clean and z-contained. This is the
   next substantive thing to do.
   **Promoted above Tier 3A on 2026-08-25** (§8E.4): every route to the
   fibroblast count is blocked, while the shear half is pure geometry and needs
   no cell-type call. **Split it and do the geometry first** —
   8a. ~~**Build the shear map itself** (§4.3)~~ — **done 2026-08-26.**
   `pipeline/shear.py`; `docs/scratch/shear/image6_profile.py` draws it.
   **Image 6's shear varies 75× within the one field of view**, with no flow
   rate assumed, over 60 usable sections of 411 µm — rising into the waist at
   ~250 µm of arc length and falling again, exactly the form §5 predicted from
   the caliber profile. An independent check on the solve: the area ratio is
   17.9× and 17.9^1.5 = 75.5 against the solved 75.1, i.e. these near-circular
   sections reproduce the R⁻³ scaling without any formula having been applied.
   **The absolute scale is the one thing still missing** and it is a single
   multiplier µQ; the perfusion rate is being asked for.
   Two cautions on the number: the vessel is unbranched only after a 21 µm stub
   ending in a §5.1 `terminus` is read as a reconstruction artifact rather than
   a branch (a dead end carries no flow, so Q is still conserved); and 23 of the
   83 solved sections are censored near a field-of-view crop and are excluded
   above.
   8b. Join it to fibroblast density — which does need [P05](ISSUES.md#p05), so
   it waits.

## Tier 3A — the fibroblast count. **Blocked, moved below item 8 — 2026-08-25**

Opened 2026-08-25 ahead of Tier 3 because it was an eye check that cost nothing.
**All three routes to it have now failed** (§8E), so it no longer gates anything
cheaply and no longer sits ahead of the shear work, which needs no cell-type call
at all.

0. ~~Count sensor-positive cell bodies on one image~~ — **cannot be done.** The
   fibroblasts are not distinguishable by eye in the sensor channels; the premise
   that they would be was wrong. §8E.1.
0a. ~~Group adhesions by nucleus as a substitute~~ — **retracted.** A Voronoi
   artifact, and it double-counts the clover photons. §8E.2.
0c. ~~Endothelial count by Cellpose, fibroblasts by subtraction~~ — **fails.**
   Two script defects fixed ([P40](ISSUES.md#p40)); with the method working it
   finds zero cells, because the en-face wall shows no polygonal tiles
   ([P41](ISSUES.md#p41)). §8E.3.
0b. **Still live, and now the only cheap items left:** cross-tab
   `nuc_is_fibroblast` × `nuc_side` signed distance, and fibroblast fraction vs
   per-image clover noise statistics — §8D.5 items 2 and 3. Both run on data
   already exported. Neither validates the classifier; both bound how wrong it
   is.
0d. **Compute the CSR null for `dist_nuc_um`** (§8E.2). Cheap, robust to the
   misassignment that killed 0a, and it audits the *adhesion detector* rather
   than the classifier — which is upstream of more of the program.
0e. **Ask for CD31** — §9 item 4, now the main bench ask rather than a
   nice-to-have.

## Tier 3 — validation that needs no phantom

9. **Image 5 as a hard negative control**: K ≈ 0 everywhere is correct.
10. **Images 1 / 3 / 4 test whether caliber-normalisation normalises**: three
    crotches at three calibers should score the same.

## Deferred

The per-vertex curvature descriptor (§7.1). It is not on the critical path and
should not be worked on until the shear pipeline exists and shows it needs one.

---

# 8A. Tier 1 results — 2026-08-22

Scripts: `docs/scratch/tier1/{geom_profile,width_fix,nuc_signed,nuc_effect}.py`,
outputs in `docs/scratch/tier1/out/`. Scratch is gitignored per project
convention, so **the numbers below are the tracked record.** Masks and the OME
header only; no pixels decoded.

## 8A.1 The failed nuclei test used an unsigned distance — and the signed version works

**There is no nucleus-to-vessel-mesh distance anywhere in the codebase.** The
only vessel-surface distance is `measure_core`'s `dist_mesh_um`, which is
adhesion-to-mesh and unsigned (`measure_core.py:147`).
`signed_distance_to_mesh_vertices` is used in five modules but always against the
**DAPI cell mesh**, never against the vessel.

Computed both, on the one mesh where the sign is physically meaningful — image
1's `vessel_mesh__alphapoisson.npz`, whose stored per-vertex normals come from
∇χ, so positive = away from the lumen = **abluminal**. 276 nuclei, 27 called
fibroblast by the [P05](ISSUES.md#p05) classifier (9.78%):

| measure | AUC | p (permutation, 20k) | median fib | median other | Δ |
|---|---|---|---|---|---|
| **signed** | **0.680** | **0.0009** | +0.95 µm | −1.25 µm | **2.20 µm** |
| unsigned | 0.602 | 0.040 | +2.87 µm | +1.66 µm | 1.21 µm |

Signing the distance nearly doubles the separation above chance (0.180 vs 0.102)
and moves the p-value ~45×. And the relationship is **monotone** — the further
abluminal a nucleus sits, the more likely the marker-based classifier calls it a
fibroblast, against a 9.8% base rate:

| signed distance | nuclei | of which classifier-fibroblast |
|---|---|---|
| > 0 µm | 92 (33.3%) | 17.4% |
| > 1 µm | 59 (21.4%) | 22.0% |
| > 2 µm | 29 (10.5%) | 34.5% |
| > 3 µm | 15 (5.4%) | 46.7% |
| > 5 µm | 7 (2.5%) | **71.4%** |

**Two unrelated measurements agreeing.** One reads marker intensity in a
perinuclear shell; the other reads geometry against a Poisson surface. Neither is
validated, and agreement between two unvalidated measures is not validation — but
they share no machinery and no failure mode, which is the first external
corroboration [P05](ISSUES.md#p05) has ever had. Note also that 10.5% of nuclei
sit >2 µm abluminal against the classifier's 9.78% fibroblast fraction; that
closeness is suggestive and nothing more, on n = 27.

**The distribution is unimodal even when signed** (single peak spanning −2.4 to
+2.2 µm with a long negative tail), which is exactly what §7.2 predicted. The
original null result was the right observation read against the wrong test.

**Blocked for the other five images.** On a plain alpha mesh the sign is *not*
luminal/abluminal: [P21](ISSUES.md#p21) established it is a watertight solid
around the wall slab, so "outward from the solid" points **into** the lumen on
the inner sheet. Testing the other five needs Poisson meshes built for them,
which is a `batch_mesh` run, not new code.

## 8A.2 Height is measurable; in-plane width is not

Lumen height per (y,x) column — the wall present both above and below, so the
column is inside the lumen. The vessel axis lies in xy, so a z column is already
perpendicular to it and this needs no centerline.

| # | image | z extent | **height p90** | height max | inscribed R p50 | **2R / height** |
|---|---|---|---|---|---|---|
| 1 | V3_04step_002 | 105.3 | **82.9** | 100.4 | 36.4 | 0.88 |
| 2 | V3_04step_003 | 98.5 | **85.8** | 97.5 | 34.1 | 0.80 |
| 3 | V3_09step_001 | 115.0 | **96.8** | 114.1 | 42.3 | 0.87 |
| 4 | V3b_04step_001 | 80.3 | **49.4** | 79.6 | 24.7 | 1.00 |
| 5 | V3b_04step_002 | 117.3 | 23.7 † | 109.8 | 36.2 | 3.05 † |
| 6 | V3b_04step_003 | 93.9 | **64.3** | 93.3 | ~18.8 ‡ | — |

† image 5's mask is the sparsest in the set, so few columns satisfy the
wall-above-and-below test and the height is under-read. ‡ see §8A.3.

**The last column is the real finding.** On the four images where both quantities
are trustworthy, **twice the inscribed radius equals the lumen height** (0.80,
0.87, 0.88, 1.00). That is what you get if and only if z is the **minor** axis of
the cross-section — which independently confirms [P24](ISSUES.md#p24)'s reading
that the inscribed radius is the minor semi-axis, by a route sharing no machinery
with it. It is consistent with a flat duct. It does **not** prove a lithographic
fixed height: the measured heights span 49–97 µm, so either the fabricated height
varies between chips or the sparse-mask images under-read. **§4's shear-exponent
question is therefore still open**, and the pattern design file remains the way to
settle it.

**In-plane width could not be measured, by either available route**, and the two
failures point in opposite directions:

- Filling the **raw** max-projection under-reads badly on a sparse mask — median
  0.9 µm on image 5, because the dashed cadherin ring ([P28](ISSUES.md#p28)) never
  closes into a solid.
- Using the **solidified** footprint `centerline.vessel_path_curvature_xy` builds
  over-reads badly — medians 145–254 µm, max 431 µm — because its `binary_fill_holes`
  closes the region enclosed by the Y. This is exactly the caveat
  [P24](ISSUES.md#p24) attached to its own 70.6 µm figure, now confirmed
  quantitatively and on all six images.

**In-plane caliber needs a perpendicular cross-section, which needs the 3D
centerline graph.** That promotes §8 item 7 from "substantial build" to
**prerequisite**: no shear estimate is possible without it, because τ depends on
the cross-section and the max-projected footprint cannot supply one.

## 8A.2A Audit: is the inscribed radius contaminated by the filled Y interior?

Raised by §7.6 and tested with a discriminator sharing no machinery with the
footprint: a column is over lumen only if the **raw** mask has wall both above
and below it. Restricting the radius to those columns:

| # | image | n valid (shipped) | R p50 shipped | n over-lumen | **R p50 restricted** | shift |
|---|---|---|---|---|---|---|
| 1 | V3_04step_002 | 575 | 36.4 | 168 | **35.2** | −3% |
| 2 | V3_04step_003 | 179 | 34.1 | 57 | **36.5** | +7% |
| 3 | V3_09step_001 | 331 | 42.3 | 57 | **41.8** | −1% |
| 4 | V3b_04step_001 | 241 | 24.7 | 69 | **23.4** | −5% |
| 5 | V3b_04step_002 | 438 | 36.2 | 34 | **34.6** | −4% |
| 6 | V3b_04step_003 | 56 | 18.8 | 15 | **12.3** | **−35%** |

**68–92% of "valid" points fail the over-lumen test**, but the medians barely
move on five of six. Read that carefully: the discriminator is itself degraded by
[P28](ISSUES.md#p28) — wherever the dashed mask misses the wall above or below,
a column inside the *true* lumen is flagged not-lumen — so the large flagged
fraction is mostly mask sparsity, not Y-interior contamination. **The audit
cannot cleanly separate the two**, and that is its honest limit.

What it does establish:

- **The radius medians are robust** to the restriction on images 1–5 (≤7%).
- **Image 6 is not** (18.8 → 12.3 µm, −35%), which is the smallest-caliber image
  and the one that also needed a forced border pad. Treat every image-6 caliber
  number as unsettled until it is recomputed from the χ lumen (§8B).
- **§8A.2's headline result survives and tightens.** Recomputing 2R/height with
  the restricted radius on the four images where both are trustworthy gives
  **0.85 / 0.85 / 0.86 / 0.95** against the unrestricted 0.88 / 0.80 / 0.87 /
  1.00. One caveat on that: the restricted version shares the z-enclosure gate
  with the height measurement, so it is the *tighter* comparison but the less
  independent one. The unrestricted figures are the ones to quote as
  independent.

## 8A.3 Poisson runs on all six — and two new defects

`vessel_poisson.diagnose` at shipped defaults, on the saved cadherin mask:

| # | image | occupied cells | sheetness med | **core fragments** | candidates | max D | verts | faces | sign flipped |
|---|---|---|---|---|---|---|---|---|---|
| 1 | V3_04step_002 | 111,701 | 0.412 | 33 | 40 | 42.6 | 58,999 | 115,365 | no |
| 2 | V3_04step_003 | 151,090 | 0.489 | 19 | 28 | 41.4 | 46,560 | 88,680 | no |
| 3 | V3_09step_001 | — | — | 27 | — | — | 59,878 | — | no |
| 4 | V3b_04step_001 | — | — | 7 | — | — | 37,769 | — | no |
| 5 | V3b_04step_002 | — | — | 21 | — | — | 36,117 | — | no |
| 6 | V3b_04step_003 | 138,211 | 0.447 | **3** | 8 | **28.4** | 25,837 | 50,185 | no |

**No refusals — the 4σ rule was cleared on every image**, which was the go/no-go
this task existed to answer. Image 6 clears it by the least margin (max D 28.4 µm
against the 20 µm requirement) and yields only 3 z-enclosed core fragments.
Checkerboard ratios are ~1e-4 throughout, so the solve is clean everywhere.

Inscribed lumen radius (median, from `lumen_center_z`): 36.4 / 34.1 / 42.3 /
24.7 / 36.2 / **~18.8** µm. **Image 6 is roughly half the caliber of the others
and its minimum reaches 5.4 µm** — the waist is real and measurable, which is
what makes it §5's best shear image.

### New defect A — the auto border pad can exceed the field of view, silently

**`centerline.vessel_path_curvature_xy` returns zero points on image 6 at default
parameters.** `border_pad_um` defaults to the footprint distance transform's p95,
which on image 6 is **348.8 µm against a 332.8 µm field of view**. The
edge-replicated pad then dominates the padded image and the medial axis leaves
the cropped window entirely. Nothing warns; the empty result surfaces only as a
downstream `lumen_center_z` note.

Forcing any smaller pad recovers it — and the pad changes the answer everywhere,
not just at the edges:

| border_pad_um | centerline points |
|---|---|
| auto (348.8) | **0** |
| 0 | 203 |
| 20 | 266 |
| 40 | 326 |
| 80 | 343 |
| 160 | 466 |

The auto rule is also fed by the same `fill_holes` inflation as §8A.2, so it is
one root cause with two symptoms. Needs its own `ISSUES.md` entry: **a pad
ceiling relative to the field of view, and a warning when the centerline comes
back empty.**

### New defect B — the branch-node count is worse than [P29](ISSUES.md#p29) recorded

Branch nodes at default parameters: **6 / 0 / 1 / 1 / 4 / 0**. Image 1 has one
real Y junction and reports **six**; [P29](ISSUES.md#p29) recorded two for the
saved run. Image 5, a straight tube, reports four. The count is strongly
parameter-sensitive and is not usable as-is. `dist_branch_um` inherits this.

### Also worth recording

**Image 4 breaks the z-containment premise.** 1.224% of its mask lies in the last
z slice, against 0.008% on the reference — 150×. The vessel is leaving through
the bottom of the stack, which is the assumption `lumen_core`'s z-enclosure
discriminator and [P24](ISSUES.md#p24)(b) both rest on. Poisson still ran, but
image 4's orientation is the least well-founded of the six.

**Centerline validity is low across the board.** The fraction of centerline
points with a bounded lumen above and below: 59% / 28% / 67% / 48% / 57% / n/a.
On image 2 — the flare image, whose junctions are off-frame — only 28% of the
centerline has a defined Z.

## 8A.4 `dist_mesh_um` cannot answer the abluminal question

| # | image | n | p10 | median | p90 | p99 | <1 µm | <5 µm |
|---|---|---|---|---|---|---|---|---|
| 1 | V3_04step_002 | 2,147 | 0.33 | 2.30 | 5.68 | 11.34 | 23% | 86% |
| 2 | V3_04step_003 | 2,212 | 0.00 | 1.49 | 4.87 | 8.14 | 39% | 91% |
| 3 | V3_09step_001 | 3,053 | 0.00 | 1.38 | 4.89 | 15.89 | 40% | 91% |
| 4 | V3b_04step_001 | 2,019 | 0.00 | 1.30 | 4.21 | 17.64 | 42% | 93% |
| 5 | V3b_04step_002 | 4,049 | 0.46 | 4.55 | 31.90 | 45.58 | 17% | 53% |
| 6 | V3b_04step_003 | 1,872 | 0.00 | 2.92 | 33.73 | 112.40 | 26% | 67% |

Medians are 1.3–4.6 µm with 23–42% of adhesions inside 1 µm — the mode is *on*
the surface, not at a consistent abluminal offset. **But the column is unsigned,
so it could not show an offset even if there were one**, and that is the whole
point of §7.2. This table is a baseline, not evidence either way.

Images 5 and 6 are visibly different — p90 of 31.9 and 33.7 µm — but their meshes
were built at widget defaults ([P26](ISSUES.md#p26)), so this is at least as
likely to be a bad mesh as a real biological difference. Not interpretable until
those two are rebuilt.

## 8A.5 What Tier 1 changed

- **The signed-distance route is worth pursuing** (§8A.1) — it is the only lead
  in the project with two independent measurements agreeing, and it points at
  [P05](ISSUES.md#p05).
- **The 3D centerline graph is now a prerequisite, not a Tier 2 item** (§8A.2):
  no in-plane caliber, therefore no shear, without it.
- **The shear exponent is still undetermined** (§8A.2). The design file is now the
  cheapest way to settle it.
- **Poisson is viable on all six** (§8A.3) — the reconstruction is not the
  blocker.
- **Two new centerline defects** (§8A.3), which together with §7.6 triggered the
  decision in §8B: **the 2D centerline is replaced, not repaired.**

---

# 8B. The centerline is replaced, not repaired — decided 2026-08-22


## Nothing upstream depends on it

`build_mesh`'s alpha→poisson branch is
`mask → mask_to_points_um → voxel_subsample → rasterize_alpha_solid →
vessel_poisson.reconstruct`. **No centerline.** `vessel_poisson` never imports
it; the only mentions in that module are docstring prose explaining the
deliberate avoidance ("Orientation from the lumen volume, not from a
centerline … uses no `binary_fill_holes`, no `binary_closing` and no
`skeletonize`").

Every reference to `centerline` in the repo is a **consumer**:

| caller | what it does with it |
|---|---|
| `measure_core.py:222,228` | `path_curv_xy`, `dist_branch_um` |
| `gui/measure.py:322,329` | the same two columns |
| `gui/vessel_mesh.py:425` | `lumen_center_z`, to draw the overlay at a sensible Z |
| `saved_check.py:79,82` | asserts the saved dict is non-empty |
| `batch_mesh.py:83` | *produces* it as a separate output step |

Both columns it feeds are already flagged untrusted
([P08](ISSUES.md#p08), [P29](ISSUES.md#p29)), so redefining them costs nothing.

## One root cause, seven symptoms

The 2D centerline exists because when it was written **there was no trustworthy
3D surface** — the alpha shape was 42 pieces at genus 125 — so the max-projected
silhouette was the only thing solid enough to skeletonize. Making a dashed
projection solid requires `binary_fill_holes`, and that fill is the root of:

[P08](ISSUES.md#p08)(a) Euclidean junction contamination · [P08](ISSUES.md#p08)(c)
border flare · **new defect A** (the auto pad, derived from the inflated
footprint EDT) · [P29](ISSUES.md#p29) the fabricated branch node ·
[P23](ISSUES.md#p23) amputated tapers · [P27](ISSUES.md#p27) raster ordering ·
§8A.2's unmeasurable width · §7.6's radius gate.

**That reason has expired.** There is now a smooth, gap-bridged, single-component
lumen. Patching seven symptoms is worse value than replacing the input.

## The replacement, built on χ

1. **Lumen volume = χ above the iso level.** Smooth by construction, gaps spanned
   by the harmonic solve, one component, no fill, no closing, no projection. It
   already exists inside `reconstruct` and is simply not returned — the
   `extras=` out-dict is the mechanism for exposing it without touching the
   4-tuple contract `saved_check`/`batch_mesh`/the GUI all destructure.
2. **Caliber = the 3D distance transform inside that lumen** — the inscribed
   radius everywhere, not sampled at skeleton points and not gated on a
   footprint. Removes §7.6 outright.
3. **3D curve skeleton of the lumen**, spur-pruned by a **scale-aware** rule: a
   branch shorter than the local caliber is a spur. Replaces
   [P23](ISSUES.md#p23)'s absolute 4 µm / 30 µm constants with something that
   adapts to a taper, which is the thing P23 says absolute thresholds cannot do.
4. **A graph, not a point set.** Nodes = junctions and endpoints; edges =
   ordered polylines with **arc length by construction**
   ([P27](ISSUES.md#p27) closes). Distance-to-junction becomes distance *along
   the vessel* rather than Euclidean ([P08](ISSUES.md#p08)(a) closes).
5. **Endpoints get classified.** An endpoint reaching an x/y domain face is a
   **crop** → censored; one terminating inside is a real terminus. This is
   §5.1's "recognise a crop" requirement, free once the graph exists.
6. **The payoff — perpendicular cross-sections.** An ordered path has a tangent,
   so the lumen can be cut by a plane ⟂ to it, giving **area, both principal
   axes, and aspect**. That is exactly what §8A.2 could not measure by any route
   and exactly what §4's shear formula needs. **This one replacement unblocks the
   shear estimate.**

Then dr/ds, junction topology, relative shear per segment, and per-adhesion
attachment by nearest point on the graph all follow.

**Five open entries retire as superseded** rather than being fixed one at a time:
[P08](ISSUES.md#p08), [P23](ISSUES.md#p23), [P27](ISSUES.md#p27),
[P29](ISSUES.md#p29), and new defect A.

## Risks, stated up front

- **3D voxel skeletons still sprout spurious branches.** The scale-aware prune
  should beat a 30 µm constant, but it needs validating — image 5, a straight
  tube currently reporting **4** branch nodes, is the test case, and the correct
  answer there is 0.
- **Caliber inherits the Poisson bridging.** Where [P28](ISSUES.md#p28)'s missing
  59–85% of dim wall lies, the lumen boundary is the solver's inference rather
  than a measurement. Flag it per point; do not average it in silently.
- **σ = `poisson_smooth_um` = 5 µm smooths the lumen.** Fine at 20–40 µm
  calibers, marginal at image 6's minimum, which may need a smaller σ recorded
  per image — and per-image parameters have to be recorded, not chosen by eye
  ([P26](ISSUES.md#p26)).
- **Sequencing:** only image 1 has an alpha→poisson mesh today. **The other five
  need a `batch_mesh` run before any of this can be measured across the set.**
  That is the first step, not a detail — and it also unblocks replicating
  §8A.1's signed-distance result beyond n = 1.

## What gets deleted

`vessel_path_curvature_xy`, `_prune_spurs`, `lumen_center_z`. `path_curv_xy` and
`dist_branch_um` get redefined against the graph.

> **Done 2026-08-22.** All three are gone — the whole of `pipeline/centerline.py`
> was deleted, since those three functions were the module. The two columns were
> **renamed** rather than redefined in place (`path_curv_1_um`,
> `dist_junction_um`, plus `caliber_um` and the `lumen_*` set): the same header
> meaning a different thing across a boundary is [P04](ISSUES.md#p04)'s `nuc_id`
> trap, and it is not worth repeating to save a rename. Results in §8C. Per the project's working
agreement the new knobs are sliders in the Vessel-mesh tab with argparse
equivalents in `batch_mesh` that have **no upper bound**, and the GUI and
headless paths call one shared function — the `build_mesh` doctrine.

---

# 8C. Tier 2 built and measured — 2026-08-22

All three Tier 2 prerequisites are done: the six-image `batch_mesh` run (§8 item
6a), the cell-type separation redone properly (item 5), and the centerline
replacement (item 7, §8B). Code: `pipeline/lumen_graph.py`,
`pipeline/nuc_side.py`, `measure_core.graph_columns`, the Vessel-mesh tab's
graph section, `batch_mesh --steps graph`. **`pipeline/centerline.py` is
deleted.**

## 8C.1 The six-image rebuild

`alpha 3D solid -> poisson` on all six at **one recorded parameter set** —
image 1's (`subsample_um` 5.0, `alpha3d_um` 15.0, `min_component_um` 253.03,
grid 2.0, σ 5.0, iso 50, raster 2.0), min_size 12. Image 1 rebuilt
**byte-identical** to the 2026-08-21 run, so the path is deterministic and the
other five are comparable to it.

| # | image | verts | faces | interior |
|---|---|---|---|---|
| 1 | V3_04step_002 | 62,092 | 122,799 | 0.1% |
| 2 | V3_04step_003 | 54,471 | 107,744 | 0.0% |
| 3 | V3_09step_001 | 67,868 | 134,648 | 0.0% |
| 4 | V3b_04step_001 | 43,452 | 85,413 | 0.2% |
| 5 | V3b_04step_002 | 45,663 | 83,706 | 7.1% |
| 6 | V3b_04step_003 | 31,842 | 63,142 | 0.0% |

**Images 5 and 6 carry a coverage warning** — only 36.8% and 31.4% of retained
cadherin voxels lie within 2 µm of the surface, against ~62% on the reference.
They are the two sparsest masks. Their geometry is the least well-supported of
the six and every number below inherits that.

Choosing one parameter set for all six buys comparability, not correctness —
[P26](ISSUES.md#p26)'s standing caveat. It is recorded rather than tuned, which
is the only defensible option while [I03](ISSUES.md#i03)'s selection rule is
still unwritten.

## 8C.2 The 3D graph gets the topology right on four of six — and says which

Correct answers from §5: images 1, 3, 4 are Y branches (1 junction), images 2, 5,
6 are unbranched (0).

| # | image | junctions | expected | ends | crop | terminus | length µm | inscribed R µm |
|---|---|---|---|---|---|---|---|---|
| 1 | V3_04step_002 | **1** | 1 | 3 | 3 | 0 | 541 | 40.6 |
| 2 | V3_04step_003 | **0** | 0 | 2 | 2 | 0 | 410 | 42.0 |
| 3 | V3_09step_001 | **1** | 1 | 3 | 3 | 0 | 647 | 46.0 |
| 4 | V3b_04step_001 | 4 ✗ | 1 | 6 | 3 | 3 | 1010 | 28.0 |
| 5 | V3b_04step_002 | 2 ✗ | 0 | 4 | 2 | 2 | 699 | 35.2 |
| 6 | V3b_04step_003 | **0** | 0 | 2 | 1 | 1 | 377 | 30.1 |

Against the old 2D centerline this is a large improvement: image 1 reported
**2** branch nodes in the saved file and **6** at defaults (§8A.3, [P29](ISSUES.md#p29));
it now reports the one junction the vessel has. Image 5 went 4 → 2 and is still
wrong.

**The two failures are the two images that trip the z-leak warning, and that is
a diagnosis rather than a coincidence.** On images 4 and 5 the fitted lumen
leaves through the **z** faces more than through the x/y ones (18,328 vs 7,688
voxels; 9,748 vs 7,955). Once the lumen bleeds out of the stack it is not a tube
any more, so its skeleton sprouts junctions. Image 4 is the known
z-containment breach ([P24](ISSUES.md#p24)(b): 1.224% of its mask in the last z
slice, 150× the reference); image 5 is the sparsest mask, the one [P28](ISSUES.md#p28)
could not seal a cross-section on at all. The module warns and does not quietly
report the number.

**Restricting the lumen to the z-enclosed core does not fix it, and that is
measured, not assumed.** On image 4 the core is 2,584 voxels of a 498,469-voxel
lumen — 0.5%, a sliver, not a vessel. On image 5 the intersection makes the
topology *worse* (3 junctions, 13 endpoints, against 2 and 4). `lumen_core` is a
seed for orientation, not a lumen. Recorded so it is not retried.

## 8C.3 In-plane caliber is measurable now — §8A.2's blocker is gone

Perpendicular cross-sections, median per image:

| # | image | area µm² | major µm | minor µm | aspect | clipped |
|---|---|---|---|---|---|---|
| 1 | V3_04step_002 | 13,836 | 185.8 | 90.1 | 2.02 | 20% |
| 2 | V3_04step_003 | 13,538 | 175.7 | 92.3 | 1.91 | 6% |
| 3 | V3_09step_001 | 13,610 | 169.1 | 100.1 | 1.89 | 17% |
| 4 | V3b_04step_001 | 7,208 | 156.8 | 58.7 | 2.48 | 56% |
| 5 | V3b_04step_002 | 12,836 | 170.7 | 93.9 | 1.80 | 44% |
| 6 | V3b_04step_003 | 4,090 | 82.4 | 68.3 | 1.24 | 5% |

**Image 1 agrees with the only direct measurement on record.**
[P24](ISSUES.md#p24)(b) cut the enclosed lumen out of the cadherin shell on the
left arm and got 151–180 µm in plane × 83–92 µm deep, aspect 1.64–2.15, area
10,679–12,637 µm². This gives 185.8 × 90.1, aspect 2.02, area 13,836 — the minor
axis lands inside the measured range and the major and the area sit just above
it, which is what including the junction region should do. **Two routes sharing
no machinery**: one cut constant-x slabs out of the mask, the other cuts planes
⟂ to a fitted axis through the Poisson lumen.

Three cautions, all flagged per point rather than argued away:

- **Images 4 and 5 clip 44–56% of their sections**, so their width and area are
  lower bounds. Both are the z-leak images.
- **Image 5's aspect reads 1.80 where §5 records ≈3.** Its lumen is the least
  supported of the set; do not read this as a correction to §5.
- **Image 6's inscribed radius is 30.1 µm against §8A.2's 18.8 µm** from
  `lumen_center_z`. §8A.2A explicitly left image 6 unsettled pending exactly
  this recomputation, so the new number supersedes it — but image 6 is the
  smallest-caliber image and this is a 60% change, so it is the one to check by
  eye first.

**§4's shear-exponent question is now answerable and is not yet answered** —
though as of 2026-08-26 it is also less pressing than it was, because §4.3's
solve does not use an exponent: it reads the shape of each section and solves.
The exponent question survives as a question about the *chip*, not about the
calculation.
Height and width are separate columns on all six; what is missing is the
along-vessel profile rather than the median, and image 6's waist is where to
read it.

## 8C.4 The cell-type separation: the signed distance replicates 5 of 6

`nuc_side.measure` — signed distance, **exact nearest-triangle** (two-pass, not
nearest-vertex), no subsample, on the alpha→poisson mesh whose ∇χ normals make
"positive" mean abluminal. AUC is against the [P05](ISSUES.md#p05) marker
classifier; 0.5 is chance.

| # | image | nuclei | marker-fib | **AUC signed** | p | AUC vertex | AUC unsigned |
|---|---|---|---|---|---|---|---|
| 1 | V3_04step_002 | 276 | 27 | **0.695** | 4.5e-04 | 0.680 | 0.581 |
| 2 | V3_04step_003 | 259 | 13 | **0.644** | 4.0e-02 | 0.639 | 0.562 |
| 3 | V3_09step_001 | 315 | 47 | 0.506 | 0.45 | 0.504 | 0.619 |
| 4 | V3b_04step_001 | 230 | 18 | **0.675** | 6.8e-03 | 0.673 | 0.438 |
| 5 | V3b_04step_002 | 290 | 29 | **0.773** | 6.9e-07 | 0.771 | 0.456 |
| 6 | V3b_04step_003 | 174 | 17 | **0.716** | 1.8e-03 | 0.713 | 0.601 |

**Five of six, p from 6.9e-07 to 0.04.** §8A.1's n = 1 result was not a fluke.
The vertex column reproduces §8A.1's recorded 0.680 exactly, which validates the
implementation against the earlier one.

Three readings:

- **The sign is the whole effect.** Unsigned runs 0.438–0.619 with no consistent
  direction and is *below* chance on two images. §7.2's diagnosis of the original
  failed experiment is confirmed on six images, not one.
- **Nearest-triangle is a correctness fix, not the source of the signal**
  (0.695 vs 0.680; 0.773 vs 0.771). Worth having — the estimator was biased
  outward and the whole effect is 5–10 µm — but it should not be credited with
  the result.
- **Image 3 is the null**, and it is the 0.9 µm z-step file. It also has the
  highest marker-fibroblast fraction (14.9%, [P05](ISSUES.md#p05)) and the least
  depth truncation ([P03](ISSUES.md#p03)). One image at that spacing; not
  separable from field-to-field variation, and it should not be pooled without
  saying so.

## 8C.5 The endothelial calibration in §7.2 does not work — and the reason is structural

§7.2 proposed calibrating against cadherin-ring-adjacent nuclei, which are
endothelial by definition of the marker, then asking whether the rest are
displaced outward. Built and measured. **It fails at both steps.**

**Contact with the cadherin *mask* is not a discriminator.** On image 1, **232
of 276 nuclei overlap a mask voxel outright** (contact distance exactly 0.0 µm)
and 255 lie within 1 µm — 93% of nuclei "endothelial". The mask is 3,764
components and mostly speckle ([P19](ISSUES.md#p19)).

Restricting to the **junction network** — components above 150 µm³, which on the
reference is 52 of 3,764 holding 59% of the voxels — makes it selective (7–53% of
nuclei per image). The displacement test on that split:

| # | image | n reference | AUC | p (one-sided) | Δ median µm |
|---|---|---|---|---|---|
| 1 | V3_04step_002 | 114 | 0.407 | 1.00 | −0.55 |
| 2 | V3_04step_003 | 134 | 0.376 | 1.00 | −0.96 |
| 3 | V3_09step_001 | 156 | 0.501 | 0.49 | −0.01 |
| 4 | V3b_04step_001 | 123 | 0.379 | 1.00 | −4.63 |
| 5 | V3b_04step_002 | 20 | 0.348 | 0.99 | −1.13 |
| 6 | V3b_04step_003 | 29 | 0.476 | 0.66 | −0.12 |

**No displacement anywhere, and the sign is consistently the wrong way**: the
non-contacting nuclei sit slightly *more luminal* than the contacting ones.

The reason looks structural rather than fixable by a better threshold.
**Fibroblasts sit directly on the outside of the same wall the endothelium
forms**, so being near the junction network does not distinguish them — it is
where both cell types are. Touching cadherin is a statement about position
relative to the wall, which is the very thing being measured, so the "reference"
and the "test" populations are not independent of the axis they are compared on.

**This does not weaken §8C.4.** The marker comparison stands on its own and got
stronger. What is retired is §7.2's specific proposal for validating it without
a marker. Finding an independent endothelial label is now the open problem —
[P14](ISSUES.md#p14)'s unused junction-probability map is the nearest thing in
the repo, and [I02](ISSUES.md#i02) is still the entry for it.

## 8C.6 What Tier 2 changed

- **The signed-distance lead is real across the dataset** (§8C.4), and it is the
  strongest thing in this project. It is still an agreement between two
  unvalidated measures.
- **The geometric route to validating it is closed** (§8C.5). Validation needs an
  independent endothelial label, not a smarter distance.
- **In-plane caliber exists** (§8C.3), so the shear estimate is unblocked on four
  of six images.
- **Two images have an unreliable lumen** (§8C.2), diagnosed and flagged, not
  silently averaged in.
- **Five entries retire** — [P08](ISSUES.md#p08), [P23](ISSUES.md#p23),
  [P27](ISSUES.md#p27), [P29](ISSUES.md#p29) and §8A.3's defect A — because the
  code they describe is deleted.

## 8C.7 The one decision left open, deliberately — **taken 2026-08-24**

> **Resolved: `alpha 3D solid -> poisson` is adopted as the shipped
> reconstruction**, on the user's explicit decision of 2026-08-24. The Vessel
> mesh tab's default method is now `alpha 3D solid -> poisson`, so the untagged
> `vessel_mesh.npz` holds it and the `__alphapoisson` tag stops being the only
> place it lives.
>
> **Adopted by decision, not by a scored sweep.** [I03](ISSUES.md#i03)'s
> selection rule is still unwritten, and `mesh_score`'s four numbers were not
> used to choose. What the decision rests on is §6.2's table — it is the first
> reconstruction to reach [P24](ISSUES.md#p24)'s plateau (`slope_loglog`
> -0.024 against alpha 3D's -0.214), it is manifold and single-piece where
> alpha 3D is 8.3% non-manifold across 42 pieces, and it does not lose on
> coverage once read at 10 µm (98.2% vs 93.3%) — plus a visual confirmation on
> the reference image. **It does not beat alpha 3D on coverage at 5 µm**
> (91.3% vs 99.1%), and that remains true after adoption.
>
> Consequence: `batch_measure` will now find an untagged `vessel_graph.npz`, so
> **the geometry columns stop exporting blank**. Every geometry column in a CSV
> exported from here on is alpha->poisson, and every one in the
> `20260812_011*` set is not — the two must not be pooled.

## The original text, kept because the reasoning still applies to the next such decision

The canonical `vessel_mesh.npz` is **still alpha 3D**. The graph is built on the
Poisson lumen, so the new artifacts are tagged `__alphapoisson` and the
untagged files are untouched. Adopting alpha→poisson as the shipped
reconstruction is [I03](ISSUES.md#i03)'s selection decision and it belongs to
whoever runs the sweep, not to the code — so it has not been taken here.

Until it is, `batch_measure` finds no untagged `vessel_graph.npz` and the
geometry columns export blank. That is the honest state: the CSVs are not
silently carrying a mixture of two reconstructions.

---

# 8D. The fibroblast population is smaller than the classifier says — 2026-08-25

Opened by a protocol correction, not by new code. The experimenter supplied
three numbers this project did not have: the cardiac fibroblasts were **cast
into the collagen before gelation** at **0.5 × 10⁶ cells/mL** in **6 mg/mL**
collagen, and the construct was fixed after **7 days**. §3 previously recorded
the opposite arrangement — cells seeded onto a finished, perfused vessel — and
every downstream expectation about where fibroblast nuclei should sit was built
on it.

**Nothing here is measured from new code, and nothing here is settled.** It is
arithmetic over existing numbers plus a re-reading of §8A.1, assembled because
it points at [P05](ISSUES.md#p05) hard enough to be worth writing down.

## 8D.1 What the casting density predicts

0.5 × 10⁶ cells/mL is one cell per 2 × 10⁶ µm³ — **mean spacing ~126 µm**,
comparable to the vessel caliber itself and to a third of the field of view. The
imaged volume is 1024 × 1024 px at 0.325 µm and 128–293 z-slices:

| # | image | depth µm | FOV volume µm³ | expected CF, whole FOV | expected CF, gel only (~60%) |
|---|---|---|---|---|---|
| 1 | V3_04step_002 | 105.2 | 1.17e7 | 5.8 | 3.5 |
| 2 | V3_04step_003 | 98.4 | 1.09e7 | 5.4 | 3.3 |
| 3 | V3_09step_001 | 115.2 | 1.28e7 | 6.4 | 3.8 |
| 4 | V3b_04step_001 | 80.4 | 8.9e6 | 4.5 | 2.7 |
| 5 | V3b_04step_002 | 117.2 | 1.30e7 | 6.5 | 3.9 |
| 6 | V3b_04step_003 | 94.0 | 1.04e7 | 5.2 | 3.1 |

The gel-only column removes the ~40% of the field the vessel itself occupies, on
image 1's 13,836 µm² cross-section (§8C.3). So **~3–4 fibroblasts per field at
casting density**, against [P05](ISSUES.md#p05)'s classifier calling **14–41**.

## 8D.2 The count check does not survive 7 days — state that plainly

Starting from 3.4 cells and dividing for 7 days:

| doubling time | fold | CF per field |
|---|---|---|
| quiescent | 1× | 3.4 |
| 96 h | 3.4× | 11 |
| 72 h | 5× | 17 |
| 60 h | 7× | 24 |
| 48 h | 11× | 39 |

The classifier's 14–41 corresponds to doubling times of **47–82 h**, which is an
unremarkable rate for fibroblasts in dense 3D collagen. **The count discrepancy
is therefore not evidence of anything on its own** — 7 days is exactly long
enough to make "the classifier is right and the cells expanded ~8×"
indistinguishable from "the classifier over-calls and they barely divided."

An endothelial budget does not break the tie either. Image 1's perimeter (447 µm
for a 185.8 × 90.1 µm section) over a Y-branch's worth of vessel in a 333 µm
field gives ~2.7e5 µm² of wall, which at 1,000–2,000 µm² per spread cell is
**134–267 endothelial cells** against 276 nuclei observed. Suggestive that the
endothelium nearly accounts for every nucleus, but the per-cell area is a
literature guess and the range is too wide to conclude from.

## 8D.3 The depth distribution is what actually discriminates

This is the part that is not arithmetic. §8A.1 measured signed
nucleus-to-surface distance on image 1 and recorded:

- the distribution is **unimodal, one peak spanning −2.4 to +2.2 µm**
- **29 of 276 (10.5%)** sit more than 2 µm abluminal
- **7 of 276 (2.5%)** sit more than 5 µm abluminal

Under the old seeding picture that was unremarkable. Under casting, it is the
finding: cells distributed through a gel at 126 µm spacing should appear at
*every* distance, many of them tens of µm out into the collagen and well inside
a 333 µm field. **There is no such population.** A unimodal peak at zero with
±2.4 µm width is what a monolayer looks like — one cell type tiling a surface,
spread over its own thickness.

Note also that **71.4% of the 7 most-abluminal nuclei are classifier-called
fibroblasts** — the classifier's hit rate climbs steeply with distance, which is
what you would expect if a small true population sits out in the gel and the
rest of its calls are noise on wall-bound nuclei.

## 8D.4 The reading that makes the pieces consistent

Assembled, not measured. **There may be only ~5–10 fibroblasts per field**, they
are the ones sitting well abluminal, and most of the classifier's remaining ~20
calls per image are noise on endothelial nuclei. That number is close to the
*no-proliferation* casting estimate, which would mean little division in a dense
6 mg/mL gel — plausible, but read backwards from the conclusion and therefore
weak.

**It explains [P05](ISSUES.md#p05)'s instability without appealing to a bad
marker.** `pipeline/gui/nuclei.py` records that perinuclear mClover is "only
~1.4× background and essentially unimodal", and that was read as the measurement
being too poor to resolve two groups. But **a 2–4% minority cannot produce a
second mode no matter how good the marker is** — unimodality is the expected
appearance, not a defect. Cutting a single smooth distribution at
background × 1.5 then returns whatever share of the tail lies above the cut,
which is a property of the noise rather than of the cells. A 9.8% median with a
5.0–14.9% spread is consistent with that and with nothing biological.

**It also weakens §8C.4 rather than strengthening it.** [P05](ISSUES.md#p05)'s
2026-08-22 History justified the signed-distance corroboration by saying
fibroblasts "were seeded onto an existing vessel and are abluminal, so the
geometric measure has a physical reason to work." The premise is wrong, and the
corrected one cuts the other way: cells scattered through the gel should be
separated by distance-from-wall **decisively**, not at AUC 0.64–0.77. Weak
agreement on a test that ought to be easy is evidence against the measure being
agreed with.

## 8D.5 Three tests, all cheap, none affected by the proliferation confound

1. **Count sensor-positive cell bodies by eye.** The fibroblasts carry the VnTs
   sensor, so they are directly visible in the clover/ruby channels as whole
   spread cells — and at 126 µm spacing they are isolated and unmistakable. This
   needs no threshold and no code. **~5–10 per field confirms §8D.4; ~27 refutes
   it.** Do this one first; it settles the rest.
2. **Cross-tabulate `nuc_is_fibroblast` against `nuc_side`'s signed distance.**
   Both exist and are joinable per nucleus. A cell cast in bulk gel has no
   business sitting *inside* the endothelial monolayer, so classifier-called
   fibroblasts at ≈0 distance bound the false-positive rate from below, without
   validating anything else. Caveat: fibroblasts are motile and 7 days is long
   enough for some to reach and wrap the vessel, so "at the wall" is not
   automatically wrong — "inside the monolayer" still is.
3. **Check whether the fibroblast fraction tracks image noise rather than
   biology.** If 9.8% is a tail fraction, it should correlate with per-image
   brightness and noise statistics of the clover channel. Six images, both
   quantities already in hand.

## 8D.6 Also from this session

- **The 2 µm perinuclear shell is *not* inconsistent across images.** Suspected
  and checked: `nuclei.py` computes the dilation radius as
  `round(2.0 / scale[0])`, and `scale[0]` is the z spacing, which would have
  given the 0.9 µm `09step` file a 2.5× thinner in-plane shell than the rest.
  It does not, because the loader calls `to_isotropic` first, so `scale` is
  (0.325, 0.325, 0.325) on all six and the shell is ~1.95 µm everywhere.
  **Recorded so it is not re-investigated.**
- **7 days is short for the endothelium**, against the 1–2 weeks the reference
  protocol uses to reach confluence. This opens
  [P39](ISSUES.md#p39): some cadherin-mask gaps may be real holes in an
  immature monolayer rather than dim wall lost to thresholding, which is an
  alternative to [P28](ISSUES.md#p28)'s reading and is not currently on the
  table.
- **CD31 is a candidate independent endothelial label** for §8C.5's blocker —
  the reference protocol runs it routinely, so it is a cheap ask. Caveat: CD31
  is also junctional, so a "does this nucleus touch the CD31 network" test
  inherits §8C.5's structural problem. It only helps if CD31 marks enough of the
  cell body to ask whether a nucleus is *enclosed by* endothelial signal, which
  is a per-cell rather than a per-junction question. Checkable on one test image
  before committing to a restain.

---

# 8E. Three routes to the fibroblast count, all failed — 2026-08-25

§8D queued one cheap test to settle whether [P05](ISSUES.md#p05) over-calls by an
order of magnitude. It cannot be run, and the two substitutes tried in its place
also fail. **All three failures are about attribution: nothing in this pipeline
can currently say which cell an object belongs to.**

## 8E.1 By eye — the premise was wrong

§8D.5 item 1 assumed the fibroblasts would be "isolated and unmistakable" in the
clover/ruby channels, being the only cells carrying the VnTs sensor. Reported by
the experimenter, 2026-08-25: **they are not distinguishable by eye.**

That is not a neutral null. It is consistent with `pipeline/gui/nuclei.py`'s
record that perinuclear mClover sits at ~1.4× background and is essentially
unimodal, and it cuts in §8D.4's favour: a marker that cannot resolve a cell body
by eye is not going to resolve cell *type* by threshold.

## 8E.2 By nucleus proximity — a Voronoi artifact, retracted

`nuc_id` and `nuc_is_fibroblast` are exported per adhesion, so grouping
adhesions by nucleus looked like a free substitute. It produced a clean result —
classifier-called fibroblast nuclei owning **1.5–6× more adhesions** than
endothelial ones, same direction in all six images — and **that result is
withdrawn.** Two independent reasons:

- **`nuc_id` is a nearest-neighbour join, not ownership.** Nucleus spacing is
  **19–25 µm** on the wall (√(area/N), 174–315 nuclei over ~332×332 µm) against a
  spread fibroblast footprint of 50–100 µm, so a cell's peripheral adhesions are
  routinely nearer some other nucleus. The data shows it directly: per-nucleus
  adhesion groups have **median radius 16–17 µm in every image**, which is half
  the nucleus spacing — i.e. the groups are Voronoi cells, and their size is set
  by nucleus density rather than by cells.
- **It scores the same photons twice.** The classifier thresholds perinuclear
  mClover; the adhesions are detected in that same channel. A nucleus with more
  clover signal is called fibroblast *and* accumulates more clover-channel
  detections. This is §8C.5's structural trap in a new place, and it would
  produce the observed effect under the noise hypothesis too.

**This qualifies [P04](ISSUES.md#p04).** That entry is closed as fixed because
`nuc_id` made cell type joinable per adhesion. It is joinable; it is not
*correct*. Opened as [P42](ISSUES.md#p42).

One quantity survives and is worth pursuing: `dist_nuc_um` has median
**4.5–6.2 µm** across the six, which looks *below* what chance would give at this
nucleus density. It is robust to the misassignment above, since it does not
depend on which nucleus is the right one. If it holds against a proper CSR null
it says the detected objects are **perinuclear** — backwards for focal adhesions,
and consistent with the detector fragmenting perinuclear cytoplasmic signal. The
null has not been computed; the expectation was estimated by hand.

## 8E.3 By endothelial subtraction — the monolayer has no cells to segment

`Scripts/fibroblast_cellpose.py` counts endothelial cells from VE-cadherin and
gets fibroblasts as `total nuclei − endothelial cells`. It had **never been
run** and was cited in neither document. Run on image 6 it fails, and two real
defects were found and fixed on the way ([P40](ISSUES.md#p40)):

- default 2D+stitch segments **axial** slices, where the endothelium is edge-on
  and there are no polygons to find — it returned 128 "cells" that were stitched
  fragments of the wall line;
- `--enface` degenerated because `_tissue_mask`'s Otsu marked **40% of the
  volume** as tissue (it thresholds the noise — the same global-threshold failure
  [P36](ISSUES.md#p36)/[P28](ISSUES.md#p28) record for this channel), flattening
  the depth profile so the top/bottom split collapsed into a whole-stack MIP.

Both are fixed — the tissue mask now reuses the pipeline's own
`cadherin_labels.npz` (3.0% of volume), and the split guards against an endpoint
minimum. **With the method working as designed, Cellpose returns zero objects**,
and still zero at `--min-vox 0 --cellprob -2`.

Two readings of the QC projection, and the second is the one that matters:

- The depth split is **wrong for this geometry** — image 6's vessel curves
  diagonally, so it has no top and bottom wall at separable depths. Doing this
  properly means unwrapping the surface using the Poisson mesh the project
  already has, not splitting by z.
- **Where the wall *is* seen en-face, the texture is fibrous filaments and
  speckle, with no closed polygonal tiles.** There is nothing cell-shaped to
  segment. This is independent of the split problem, and it is consistent with
  [P39](ISSUES.md#p39) (not confluent at 7 days) or with the junctions being too
  dim and dashed ([P28](ISSUES.md#p28)).

Recorded as [P41](ISSUES.md#p41). **Do not spend more time tuning this
script** — the blocker is image content, not parameters.

## 8E.4 What this changes

- **§8D is unresolved and its test is gone.** Whether P05 over-calls ~10× is
  still the gating question and now has no cheap answer in hand.
- **Tier 2 item 5a hardens.** An independent endothelial label was the blocker;
  the cheapest in-repo candidate has now been tried and failed. CD31 (§9 item 4)
  moves from "worth asking" to **the main ask**, with §8D.6's caveat that it only
  helps if it marks cell body rather than junction.
- **The spatial half of the program does not depend on any of this.** §4.3's
  shear map is pure geometry and needs no cell-type call. It should proceed
  ahead of the count, not behind it — that reorders Tier 3A below Tier 2 item 8.

---

# 8F. P03's assumed fix is dead, and the adhesion map is unchanged — 2026-08-25

Short, because the detail belongs in [P43](ISSUES.md#p43), [I07](ISSUES.md#i07)
and [I09](ISSUES.md#i09). What matters at program level is three things.

**1. The adhesion artifacts on disk were two weeks stale, and re-running them
changed nothing.** Every `adhesion_instances.npz` dated 2026-08-11, against the
clover-gate fix in commit `46f99cc` (2026-08-24). Re-detection on image 6 with a
30× wider gate returned a **bit-identical** `proba >= 0.5` set — 0 voxels
different. So the CSVs' adhesion positions were never affected by the gate, and
everything measured off them in §8E stands.

**2. [P03](ISSUES.md#p03) is still open and its planned fix is refuted.**
[I07](ISSUES.md#i07) held that the global Otsu gate was the root cause of the
depth truncation. It was not the cause; it was never even binding. The
truncation is downstream — the RF itself, or the post-process `min_size`. The
live hypothesis is now **annotation sampling**: training fits hand-painted
voxels on the raw volume, so if the examples were painted on shallow bright
adhesions the model simply never learned the deep dim ones. If that is right the
fix is *retraining*, not a threshold — a materially bigger job than the queue
assumed.

**3. This does not block the shear work, and that is the point.** §4.3's map is
pure geometry. P03 bites only when adhesions are *joined* to it (queue 8b), and
it bites as a depth-dependent detection efficiency — which is exactly the
confound §7.4 already requires be modelled rather than assumed away. **8a is
still clear to run.**

Also from this session: the clover gate is removed
([P43](ISSUES.md#p43)) because it saved ~10% of runtime while deciding a priori
which voxels could be positive; `chunk_z` is now sized automatically against
host and VRAM budgets, with 12 GB cards as the target; and two proposed
speed-ups were measured and rejected in favour of one that changes no feature
value ([I09](ISSUES.md#i09)).

---

# 9. Open questions for the bench

Send as one message, not three:

1. ~~**Was flow running during and after fibroblast seeding?**~~ **Retired
   2026-08-25** — the premise was wrong. The fibroblasts were cast into the
   collagen before gelation (§3), so they predate the endothelium and there was
   no separate seeding event for flow to be running during. What replaces it:
   **(a) how much flow ran, and for how long, over the 7 days?** — still decides
   whether shear is a live force or a fossil imprint. **(b) Were the fibroblasts
   counted or viability-checked at any point?** — §8D says the expected number
   per field is ~3–4 before proliferation and the classifier reports 14–41, and
   any independent count collapses that question.
2. **Does the pattern design file exist, and can we get it?** Ground-truth
   caliber, inlet/outlet (→ flow direction), bifurcation vs confluence, flow
   split. Would also give [I03](ISSUES.md#i03) an external validation target.
3. **Can a no-flow control be run** — identical chip and seeding, never perfused?
   §7.5. Long lead time, so ask now.

4. **Is there a VE-cadherin-independent endothelial label available** — CD31
   would do, and the reference protocol runs it routinely? §8C.5's blocker is
   exactly this, and it is a far cheaper ask than item 3.

Also worth confirming: flow rate Q and flow direction for future datasets, both
of which the export should accept as optional inputs rather than assume. And
whether the endothelium was confluent at 7 days ([P39](ISSUES.md#p39)).

---

# 10. Design principles for what gets built

- **Q and flow direction are optional inputs, not assumptions.** With neither:
  within-image relative shear (well-defined by conservation, §4.1). With
  direction: converging vs diverging becomes signed, carina vs lateral wall
  becomes assignable. With Q: absolute τ, comparable across images. Current data
  produces real results; the future dataset drops into the same code path with
  two more numbers filled in.
- **Keep width and height separate** (§4) — different exponents.
- **Preserve circumferential position**, not just depth (§7.4).
- **Report the vector, never a total.** Same doctrine as `mesh_score`: a
  weighting is a hidden constant, and the selection rule belongs to whoever runs
  the analysis.
- **Cluster by image.** Effective n for crotch is 3 (§5.1).
- Everything else in `CLAUDE.md`'s working agreement still holds — sliders in
  napari, no GUI/headless drift, napari as the only viewer, keep it narrow.

---

# History

- **2026-08-25c** (latest) — **the adhesion map was stale, re-running it changed
  nothing, and [P03](ISSUES.md#p03)'s planned fix is refuted.** New **§8F**.
  Every adhesion artifact predated the 2026-08-24 clover-gate fix; re-detection
  with a 30× wider gate gave a bit-identical `proba >= 0.5` set, so the gate was
  never binding and [I07](ISSUES.md#i07) closes with its premise wrong. The
  depth truncation is downstream, and the live hypothesis is annotation
  sampling — which would make the fix retraining rather than a threshold.
  The gate is removed ([P43](ISSUES.md#p43)); `chunk_z` is auto-sized for 12 GB
  cards; and two proposed speed-ups were measured and rejected for one that
  changes no feature value ([I09](ISSUES.md#i09)). **None of this blocks queue
  item 8a**, which is pure geometry.

- **2026-08-25b** — **every route to the fibroblast count failed, and
  the shear half was promoted past it.** §8D queued an eye count to settle
  [P05](ISSUES.md#p05); the experimenter reports the fibroblasts are not
  distinguishable by eye, so the premise was wrong. Two substitutes were tried
  and both fail: grouping adhesions by `nuc_id` is a Voronoi artifact that also
  double-counts the clover photons (**retracted**, and it qualifies
  [P04](ISSUES.md#p04) — new [P42](ISSUES.md#p42)), and Cellpose endothelial
  subtraction finds zero cells because the en-face wall has no polygonal tiles
  ([P41](ISSUES.md#p41)), after two real script defects were fixed
  ([P40](ISSUES.md#p40)). New **§8E**.

  **The reordering is the point.** §4.3 was added answering how τ_w is actually
  computed on this reconstruction — a per-section Poisson solve, `µ∇²u = dp/ds`,
  `τ_w = µ|∇u|`, normalised by conservation along an unbranched segment so **Q
  drops out of the profile shape**. That is pure geometry and needs no cell-type
  call, so queue item 8 was split: **8a (the shear map) now runs ahead of Tier
  3A**, which is blocked. Two things §4.3 needs do not exist yet:
  `cross_section` discards the section mask and never computes perimeter.

  One quantity worth chasing survives §8E.2: `dist_nuc_um` median 4.5–6.2 µm
  looks below chance for this nucleus density, which would mean the detected
  adhesions are perinuclear. The null has not been computed.

- **2026-08-25** — **a protocol correction that moves
  [P05](ISSUES.md#p05), plus a paper read.** Two separate things, recorded
  together because the second set up the first.

  **Read `docs/papers/ncomms8858.pdf`** (Zheng, Chen & Lopez, Nat Commun 6:7858,
  2015) — a *different lab*, whose in-vitro microvessel system this project's
  method descends from but whose question is not ours: it has endothelium only,
  no fibroblasts and no tension readout. New sections **§3.1** (a pattern design
  file would be a prior, not the lumen — the mould defines the collagen channel
  and the endothelium is subtracted from inside it) and **§4.2** (this system is
  pressure-driven, so τ_w = (Δp/L)·D_h/4 and wall shear is *linear* in hydraulic
  diameter; §7.3's R⁻³ amplification is a fixed-Q result, and the cubic
  sensitivity reappears only at the crotch, where the flow split between
  daughters goes as h³). §7.3 carries a qualifier block so it does not read as
  contradicting §4.2.

  **The correction, and it is a real one: the fibroblasts were cast into the
  bulk collagen before gelation** — 0.5 × 10⁶ cells/mL, 6 mg/mL collagen, fixed
  at 7 days. §3, §1 and `CLAUDE.md` all recorded them as *seeded onto an
  already-formed, already-perfused vessel*, which was wrong. "Abluminal" and the
  whole causal chain survive; what changes is that transmural distance is a wide
  covariate rather than a small offset, their adhesions are largely to collagen
  rather than to the endothelium, and a contractile population cast through a
  6 mg/mL gel compacts the matrix more than §3.1's reference system would.

  **§8D is new and follows from it.** Casting density predicts ~3–4 fibroblasts
  per field; the classifier reports 14–41. That gap is **not** evidence by
  itself — 7 days at a 47–82 h doubling time reproduces the classifier's numbers
  exactly, and the check is stated as inconclusive. What does discriminate is
  §8A.1's depth distribution: unimodal, ±2.4 µm, with only 2.5% of nuclei more
  than 5 µm abluminal. Cells scattered through a gel at 126 µm spacing should
  appear at every distance and do not. The assembled reading — ~5–10 real
  fibroblasts per field, the rest of the calls noise — also explains why
  perinuclear mClover is unimodal (a 2–4% minority cannot make a second mode)
  and why the fraction swings 2× on a slider step (it is a tail fraction). It
  **weakens** §8C.4 rather than strengthening it: distance-from-wall ought to
  separate cast-in-gel cells decisively, so AUC 0.64–0.77 is a poor result on an
  easy test.

  Three cheap tests are queued as Tier 3A, the first of which is an eye count of
  sensor-positive cell bodies and settles the rest. **Nothing in §8D is
  measured** and P05's status is unchanged pending that count.

  Also: §7.5 gains a partial in-image control (far-from-vessel fibroblasts, if
  they exist), §9 item 1 is retired as unanswerable and replaced, CD31 is added
  as a candidate independent endothelial label for §8C.5, and
  [P39](ISSUES.md#p39) opens on 7 days being short for a confluent monolayer.

- **2026-08-25** (latest) — **read `docs/papers/ncomms8858.pdf`** (Zheng, Chen &
  Lopez, Nat Commun 6:7858, 2015), the lab's own account of this microvessel
  system, and recorded two things it changes. **§3.1**: the pattern design file
  §3 has been asking for is a *prior*, not the lumen — the mould defines the
  collagen channel and the endothelium is subtracted from inside it, so the true
  lumen is always smaller, and our aspect-ratio reproducibility and our 49–97 µm
  height spread argue in opposite directions about how much the design would
  settle. **§4.2**: this system is pressure-driven, not flow-driven, so
  τ_w = (Δp/L)·D_h/4 and wall shear is *linear* in hydraulic diameter —
  §7.3's R⁻³ amplification is a fixed-Q result and does not apply. The cubic
  sensitivity reappears at the crotch, where the flow split between daughters
  goes as h³. Nothing measured changed; what changed is how much a caliber error
  is believed to cost, and where.

  Also worth noting for the eventual methods section: the paper's grid
  experiment is the template for the analysis §8 is heading toward — solve the
  flow once, reduce it to one wall-shear number per vessel segment (312 of
  them), and correlate against a per-segment imaging measurement. It has no
  fibroblasts and no tension, so its thresholds are not transferable as
  covariate cutoffs.

- **2026-08-24** — **a session of defect-finding in the mask pipeline,
  and one adoption.** Nothing here changes the scientific program; it changes
  what the program has been running on.

  **The finding that matters: the cadherin "local threshold" was inert**
  ([P36](ISSUES.md#p36)). `sauvola AND floor` is bit-identical to `floor` alone
  at every setting tested, because Sauvola's normaliser is half the *global*
  intensity range — 3,777 against a median local σ of 16 — so the adaptive term
  contributes 0.4% of `k` and the threshold falls below background. **Every
  cadherin mask in this project was a pure global intensity cut**, and
  therefore inherited [P03](ISSUES.md#p03)'s depth bias without anyone
  accounting for it. Niblack is *not* inert and was the tab's interactive
  default, but `_run_cadherin_step` forced Sauvola for both `batch.py` and the
  GUI's "Run all" — so the method being tuned and the method being shipped were
  different, and the shipped one did nothing.

  This reframes [P28](ISSUES.md#p28) rather than overturning it. Its
  attribution of the loss to `floor_k` and `min_size` stands and is now the
  whole story; its conclusion that "the local threshold is innocent" was
  numerically right and read wrong. It also means **every previous attempt at a
  more permissive cadherin mask was made with the local threshold switched
  off**, which is why they all degraded the reconstruction.

  **[P28](ISSUES.md#p28) is downgraded to a nice-to-have for reconstruction**,
  on the user's argument: a false positive (lumen speckle) *moves* the surface
  and can break orientation outright, while a false negative (dropped dim wall)
  only leaves a patch unconstrained, and Poisson approximates rather than
  interpolates, so it bridges smoothly. Recovering the dim wall buys a
  better-*evidenced* surface, not a different one. Two caveats keep it from
  being optional: the missing wall is depth-systematic rather than random, and
  past a certain hole size the level set escapes — which is
  [P32](ISSUES.md#p32), not a confidence question.

  **The clover gate's depth bias is fixed and the cause is confirmed optical**
  ([I07](ISSUES.md#i07), [P03](ISSUES.md#p03)). Holding wall orientation fixed,
  adhesion density still falls with absolute depth on 4 of 6 images — gravity,
  the only rival explanation, cannot do that. The gate is now a local threshold:
  acceptance shallow/deep goes **3.06 → 0.91** at 98.85% recall. **Not yet
  demonstrated end-to-end** — that needs re-detection.

  **§8C.7 is taken**: alpha->poisson is adopted as the shipped reconstruction.

  Also: [P33](ISSUES.md#p33)–[P37](ISSUES.md#p37) opened (a Sauvola PCIe
  transfer bug, an eager centroid loop costing 34 s a run, an unbounded alpha
  3D hang, the inert threshold, and a loader that overwrote `floor_k` on every
  load), [I04](ISSUES.md#i04) closed with `batch.py --steps`, and
  [I08](ISSUES.md#i08) opened to remove "Refine with cadherin". Sliders in the
  Cadherin and Vessel-mesh tabs are now scoped to the active dropdown.

- **2026-08-22** — **Tier 2 built and measured (§8C).** All six images
  rebuilt with `alpha 3D solid -> poisson` at one recorded parameter set; the 2D
  centerline **deleted** and replaced by `lumen_graph`, a 3D graph on the
  Poisson lumen with arc length, scale-aware pruning, crop classification and
  perpendicular cross-sections; the cell-type separation redone with a signed
  exact nearest-triangle distance. Three results: the signed-distance agreement
  with the marker classifier **replicates on 5 of 6 images** (AUC 0.64–0.77,
  best p 6.9e-07) and the unsigned version is chance, confirming §7.2's
  diagnosis on the whole dataset; **in-plane caliber is measurable at last**,
  agreeing with [P24](ISSUES.md#p24)(b)'s direct cut on image 1 by a route
  sharing no machinery; and the graph recovers the right topology on 4 of 6,
  failing on exactly the two images whose lumen leaks through the z faces, which
  it warns about rather than reporting silently. One negative worth as much as
  the positives: **§7.2's endothelial calibration does not work** — fibroblasts
  sit on the outside of the very wall the calibration uses as a reference, so
  the two populations are not independent of the axis they are compared on.
- **2026-08-22** — **decision: the 2D centerline is replaced, not
  repaired (§8B).** Confirmed nothing upstream of the mesh uses it — every
  reference is a consumer, and `vessel_poisson` orients from the lumen volume
  with no centerline at all. Seven known defects trace to one decision (working
  on a max projection of a dashed mask), and the reason for that decision — no
  trustworthy 3D surface — no longer holds. Also audited §7.6's radius-gate
  concern (§8A.2A): the medians are robust on five of six images and the
  2R ≈ height result tightens rather than weakening, but image 6 moves −35% and
  stays unsettled.
- **2026-08-22** (later) — **Tier 1 run, all four items; results in §8A.** The
  signed-distance lead is real (AUC 0.680, p 0.0009, monotone dose-response
  against the marker classifier); in-plane width proved unmeasurable from the
  max-projected footprint, which promotes the 3D centerline graph to a
  prerequisite; Poisson clears the 4σ rule on all six; two new centerline defects
  found (an auto border pad that exceeds the field of view and silently returns
  nothing, and a branch-node count worse than [P29](ISSUES.md#p29) recorded).
- **2026-08-22** — opened. Consolidates the live findings from `ISSUES.md` with
  the 2026-08-22 design discussion that reframed the program from *curvature* to
  *shear*. The reframe chain: curvature is fitted at 5 µm and measures roughness
  → cannot be coarsened by raising `radius_um` → Gaussian curvature is the
  orientation-free coarse invariant → **but K = −r″/(r(1+r′²)²), so it is blind
  to a taper and flags crotch/flare/waist alike** (§7.1) → crotch vs flare is
  topological → **and the real target is shear, whose dominant term is dr/ds,
  which K cannot see** → patterned fabrication changes the shear exponent and
  offers ground truth → fibroblasts are abluminal, so the endothelium is the
  sensor → the failed nuclei separation has a mechanical explanation (§7.2) and
  fixing it may close [P05](ISSUES.md#p05). Net effect: the curvature descriptor
  leaves the critical path; cell-type classification and the caliber profile
  join it.
