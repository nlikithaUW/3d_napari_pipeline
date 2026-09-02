# The two Poisson solves: the vessel surface, and the wall shear

Two places in this pipeline solve a Poisson equation:

- `pipeline/vessel_poisson.py` fits the **vessel surface** to the cadherin data.
- `pipeline/shear.py` computes the **wall shear stress** on a cross-section of
  the lumen.

They share one mathematical operator and nothing else — different unknowns,
different sources, different regions, different boundary conditions, different
solvers, and no shared code. This document explains both from the ground up and
then states precisely how they relate.

Every technical term is defined at first use. Nothing here is an analogy; every
statement is meant literally.

---

## 1. Vocabulary

These are used by both solves. Read this section once; the rest refers back to
it.

### 1.1 Scalar field, vector field

A **scalar field** is a single number assigned to every point of space. On our
data, space is a grid of cells, so a scalar field is a 3D array of numbers. The
cadherin occupancy is a scalar field; so is χ, defined below.

A **vector field** assigns a triple of numbers (one per axis: z, y, x) to every
point. It has both a direction and a length at each point.

### 1.2 Gradient (∇f)

For a scalar field f, the **gradient** ∇f is the vector field whose components
are the rate of change of f along each axis:

    ∇f = (∂f/∂z, ∂f/∂y, ∂f/∂x)

On a grid, ∂f/∂z at a cell is computed as a difference between neighboring cell
values divided by the cell spacing. The gradient vector points in the direction
along which f increases fastest, and its length is how fast f increases in that
direction. Where f is flat, the gradient is zero.

### 1.3 Divergence (∇·V)

For a vector field V = (V_z, V_y, V_x), the **divergence** ∇·V is a scalar field:

    ∇·V = ∂V_z/∂z + ∂V_y/∂y + ∂V_x/∂x

It is a single number per cell, measuring the net outward change of V there.

### 1.4 Laplacian (∇²f)

The **Laplacian** of a scalar field f is the divergence of its gradient:

    ∇²f = ∇·(∇f) = ∂²f/∂z² + ∂²f/∂y² + ∂²f/∂x²

On a grid with spacing h, the standard discrete form in 2D is

    ∇²f at pixel i = ( f(left) + f(right) + f(up) + f(down) − 4·f(i) ) / h²

and in 3D the same with six neighbors and −6·f(i). So the Laplacian at a cell is
(the sum of the neighbor values) minus (the number of neighbors times the value
at the cell), divided by h². It is zero exactly when the value at a cell equals
the average of its neighbors.

The 2D version is called the **five-point stencil** (four neighbors plus the
center); the 3D version is the **seven-point stencil**.

### 1.5 The Poisson equation

A **Poisson equation** is

    ∇²f = g

where g is a known scalar field ("the source" or "the right-hand side") and f is
the unknown field to be found. Solving it means finding the f whose Laplacian,
everywhere, equals the prescribed g.

This is one equation per grid cell, and each equation couples a cell to its
neighbors. So a Poisson solve is a large system of simultaneous linear
equations — one unknown and one equation per cell — and it must be solved
globally: the value at one cell depends on the source everywhere.

### 1.6 Boundary conditions

The equation ∇²f = g alone does not pick out a unique f. Something must also be
said about f at the edge of the region being solved on. Two kinds appear here.

- **Dirichlet condition**: the *value* of f is prescribed on the boundary. Both
  solves in this document use the Dirichlet value zero where they use it at all.
- **Neumann condition**: the *derivative of f across the boundary* is
  prescribed — here always zero, meaning f is required to have no slope
  perpendicular to the boundary face. This is also called a "zero-flux"
  condition.

A purely Neumann problem does not determine f uniquely either: if f is a
solution, so is f plus any constant, because adding a constant changes neither
the Laplacian nor the boundary slope. That leftover constant is called the
**null space** of the problem, and something outside the equation has to fix it.
Section 2.6 says how.

### 1.7 Level set

Given a scalar field f defined throughout a 3D volume and a chosen number c, the
**level set of f at c** is the surface consisting of all the points where
f = c. It separates the region where f > c from the region where f < c.

**Marching cubes** is the standard algorithm for extracting a level set from a
3D grid as a triangle mesh: it walks over each 2×2×2 group of cells, finds where
the value c falls between cell values, interpolates the crossing position along
each grid edge, and emits triangles connecting those crossings. The 2D version,
which finds a curve rather than a surface, is **marching squares**. Both give
positions between cell centers rather than snapping to them, which is what
"sub-pixel" means below.

### 1.8 Eigenvector and eigenvalue

A 3×3 symmetric matrix M acts on a vector x to give another vector Mx. An
**eigenvector** of M is a nonzero vector x whose direction M does not change:

    M x = λ x

for some number λ, called the corresponding **eigenvalue**. A symmetric 3×3
matrix has three mutually perpendicular eigenvectors with three real
eigenvalues. Sorting the eigenvalues from smallest to largest, λ1 ≤ λ2 ≤ λ3,
the eigenvector belonging to λ3 is the direction in which whatever M measures is
largest, and λ3 says how large.

`numpy.linalg.eigh` computes these for a symmetric matrix, and it is what
`sheet_normals` uses (`vessel_poisson.py:210`).

---

## 2. The surface solve — `pipeline/vessel_poisson.py`

### 2.1 What the input is and what the output is

The input is the **cadherin mask**: a 3D binary array marking the voxels where
VE-cadherin signal was segmented. VE-cadherin sits at the junctions between
endothelial cells, and those cells tile the vessel wall, so the mask is a thin,
broken, noisy shell of marked voxels lying on the vessel surface. It is broken
because the junction network is not a solid sheet and because segmentation drops
weak signal: there are gaps of roughly 8–10 µm.

The desired output is a triangle mesh of the vessel wall — a closed, smooth
surface that passes through the cadherin data, spans its gaps, and does not
follow every speck of noise.

### 2.2 χ: what it is

**χ (chi) is a scalar field defined on the whole image volume, whose value
indicates whether a point is inside the lumen or outside it.** It is intended to
be near a high value inside the lumen, near a low value outside, and to pass
smoothly through intermediate values across the wall. A function like this — one
that indicates membership in a region by its value — is called an **indicator
function**. Because it is smooth rather than a hard 0/1 step, it is sometimes
called a *smooth* or *soft* indicator.

χ is not a physical quantity. It is not fluorescence, density, or probability.
It is a constructed field whose only purpose is that **one of its level sets is
the vessel surface**. Nothing about χ's absolute values matters; only where it
crosses the chosen level matters.

Two things in the pipeline are read off the same χ:

- The **surface**: the level set of χ, extracted by marching cubes.
- The **lumen volume**: the set of cells where χ exceeds that same level. This
  is what `lumen_graph.py` builds the centerline on, and it is what the shear
  solve in §3 ultimately cuts its cross-sections out of.

Because both come from one field, the surface and the centerline cannot disagree
about where the vessel is. That is the reason `reconstruct` hands χ back to the
caller through its `extras` out-dict (`vessel_poisson.py:522`) rather than
discarding it.

### 2.3 Step 1 — the coarse grid (occupancy)

The full image is hundreds of millions of voxels, which is too many to solve on
and finer than the wall is smooth. `occupancy_grid` (`vessel_poisson.py:166`)
therefore averages the mask over blocks of k×k×k voxels, where k is chosen so
the resulting cell spacing h is close to `grid_um` (default 2.0 µm).

The resulting field is called **occ**. Because the average is taken over the
mask *values*, and the mask is 0 or 1, occ at a coarse cell is exactly the
fraction of that cell's voxels that were marked as cadherin — a number between 0
and 1. Everything downstream uses occ, not the original mask.

### 2.4 Step 2 — sheet normals from the structure tensor

**A normal is a unit vector perpendicular to a surface at a point.** The solve
in §2.6 needs, at each occupied cell, the direction perpendicular to the local
patch of vessel wall. That direction is the **sheet normal**; "sheet" because
the cadherin occupancy forms a thin two-dimensional layer embedded in three
dimensions.

The obvious candidate, ∇occ, does not work, and the reason is structural. The
cadherin density is a **ridge**: it is low in the lumen, rises to a maximum
inside the wall, and falls again outside. So ∇occ points one way on the luminal
face of the wall and the opposite way on the abluminal face, and it is zero at
the ridge crest, which is the middle of the wall. It describes the two faces of
the wall separately, and it is not the normal of the sheet.

What is needed is the **axis** along which occ varies, without regard to which
way along that axis it happens to be increasing. The **structure tensor**
supplies exactly that. At each cell it is the 3×3 symmetric matrix

    J_ab = G_σ * ( (∂_a d)(∂_b d) )

read as follows:

- d is occ pre-smoothed by one cell, to suppress single-cell noise.
- (∂_a d)(∂_b d) is the product of two components of the gradient of d. For
  a, b each ranging over z, y, x, these nine products (six distinct, since the
  matrix is symmetric) form a 3×3 matrix at every cell. This construction is
  called the **outer product** of the gradient with itself.
- G_σ * (…) means each of those six product fields is smoothed with a Gaussian
  filter of width σ = `smooth_um` (default 5.0 µm). This is the essential step:
  it averages the gradient's *directional information* over a neighborhood
  instead of averaging the gradient itself. Averaging gradients across the wall
  would cancel the two opposite-facing sides against each other; averaging their
  outer products does not, because the outer product of a vector with itself is
  unchanged when the vector is negated.

The result at each cell is a matrix summarizing which directions the local
density varies along, over a neighborhood of size σ. Take its eigenvectors and
eigenvalues (§1.8):

- The eigenvector for the **largest** eigenvalue λ3 is the direction of greatest
  variation. For a thin sheet, that is the direction perpendicular to it. **That
  eigenvector is the sheet normal.**
- Its sign is arbitrary — an eigenvector negated is still an eigenvector — so
  the normals come out pointing to one side or the other at random. Step 3 fixes
  that.

The eigenvalues also give a free quality measure. `sheet_normals` returns

    sheetness = (λ3 − λ2) / λ3,     a number in [0, 1]

- For a genuine sheet, occ varies strongly across it and weakly along it, so
  λ3 is much larger than λ2 and sheetness is near 1.
- For an isolated speck of noise, occ varies about equally in all three
  directions, so λ1 ≈ λ2 ≈ λ3 and sheetness is near 0.

Sheetness later multiplies each normal's weight, so specks contribute almost
nothing to the solve. It requires no threshold and no tunable parameter, which
is why it is the pipeline's speckle rejector.

### 2.5 Step 3 — orienting the normals inward

Each normal's sign is arbitrary, but the solve needs them all pointing the same
way relative to the vessel: **into the lumen**. To do that, the code first has
to identify the lumen as a volume.

`lumen_core` (`vessel_poisson.py:261`) does this without a centerline, without
`binary_fill_holes`, and without skeletonization — all of which have caused
recorded defects elsewhere (see `docs/DIRECTION.md` §8B):

1. Threshold the smoothed occupancy to get a **wall** region.
2. Compute D, the **distance transform**: at every cell, the distance to the
   nearest wall cell.
3. Candidate interiors are the connected regions where D > 2σ — cells well away
   from any wall.

That last test alone selects the wrong thing: the corners of the image are more
than 100 µm from any wall, further than the lumen's own inscribed radius of
38.1 µm, so the *exterior* passes it more strongly than the lumen does. The
discriminator is **z-enclosure**. The vessel does not touch the top or bottom of
the imaged stack — 0.01% of the mask lies in the first z slice and none in the
last (`docs/ISSUES.md` P24), so there is roughly 7 µm of clearance above and
below. Therefore, for a cell in the lumen there is wall both above and below it
in its own z column, while any exterior region runs the full z extent. A
candidate region is kept only if it never touches a z face and reaches D ≥ 4σ
somewhere.

If nothing qualifies, `lumen_core` refuses and returns nothing. There is
deliberately **no fallback**: a guessed orientation would put a silent global
sign error on everything downstream.

With the lumen core in hand, `orient_inward` (`vessel_poisson.py:374`) computes
Dc, the distance from each cell to the core, and flips any normal whose
direction agrees with ∇Dc — that is, any normal pointing away from the core.
Distance to the core *volume* is used rather than to a one-dimensional
centerline, because at the Y junction the nearest point of a 1D axis can lie
across the wall on the other branch, whereas the nearest cell of the lumen
volume cannot.

### 2.6 Step 4 — splat, and the Poisson equation being solved

Now the connection between the normals and χ.

Suppose χ were exactly 1 inside the lumen and 0 outside, with a smooth
transition across the wall. Then:

- ∇χ is zero everywhere except at the wall, since χ is flat elsewhere.
- At the wall, ∇χ points from low χ to high χ, that is, **inward**.

So ∇χ is precisely an inward-pointing vector field supported on the wall — which
is what the previous two steps have just measured. Call the measured field V:
each occupied cell contributes its inward unit normal, weighted by occ (how much
wall is in this cell) and by sheetness (whether it is wall at all), and the
result is smoothed by a Gaussian of width σ. This deposition of per-cell vectors
onto a grid is called **splatting**.

We want the χ whose gradient is V. Written directly, ∇χ = V is three equations
per cell for one unknown per cell — overdetermined, and V is measured with noise
so no exact solution exists. Taking the divergence of both sides collapses it to
one equation per cell:

    ∇·(∇χ) = ∇·V,   that is,   ∇²χ = ∇·V

**This is the Poisson equation the surface solve solves.** The unknown is χ; the
source is the divergence of the measured normal field. The solution is the χ
whose gradient matches V as closely as possible in a least-squares sense, which
is exactly the "approximating rather than interpolating" property that motivates
this module: noise in V averages out, and a gap in the data is bridged by the
smoothest field consistent with the normals on either side rather than by a
fabricated triangle.

**Boundary condition: Neumann, zero flux, on all six faces of the image box.**
Not Dirichlet. The vessel genuinely runs out through the x and y faces of the
field of view, so requiring χ = 0 there would force the surface to close over
and fabricate a cap exactly where the data is thinnest. The consequence, which
must be carried downstream, is that **the reconstructed surface is open at the
x/y faces** — it encloses no volume, so the usual trick of orienting normals by
"outward from the enclosed volume" is unavailable, and `reconstruct` returns its
own per-vertex normals derived from ∇χ for `cell_mesh.mesh_curvature` to use.

**The solver: a discrete cosine transform** (`poisson_dct`,
`vessel_poisson.py:414`). A **discrete cosine transform (DCT)** re-expresses a
grid of numbers as a sum of cosine waves of increasing frequency, and it is
invertible. Its relevance is this: on a rectangular box with zero-flux
boundaries, each cosine wave is an eigenvector of the discrete Laplacian
(§1.8) — applying ∇² to one of these waves returns the same wave multiplied by a
known number λ, here

    λ = −(4/h²) · Σ_axes sin²(π k / 2N)

for wave index k along an axis of N cells. So the whole coupled system becomes
independent one-number divisions:

1. Take the DCT of the source ∇·V.
2. Divide each transformed coefficient by its λ.
3. Take the inverse DCT.

That is a **direct, exact** solve of the discrete system in three fast
operations — not an iterative approximation, and there is nothing to converge or
fall back to.

Two numerical details are load-bearing:

- **The difference operators are staggered**: a forward difference for the
  gradient, a backward difference for the divergence. Their composition is
  exactly the seven-point Neumann Laplacian. Using centered differences for both
  gives an operator whose λ is zero not only for the constant field but also for
  the highest-frequency alternating pattern (the **Nyquist mode**, a field that
  flips sign at every cell). An unconstrained mode with λ = 0 is invisible to
  the equation and rides along on the solution as a checkerboard artifact.
  `_checkerboard_ratio` (`vessel_poisson.py:467`) exists to detect that.
- **The λ = 0 constant mode is left at zero.** That is the Neumann null space
  from §1.6: the solve cannot determine χ's overall offset. Leaving it at zero
  fixes the offset arbitrarily, which is harmless because the next step chooses
  the level from the data, absorbing whatever offset was chosen.

### 2.7 Step 5 — the level set, and reading the winding off ∇χ

The surface is the level set of χ (§1.7) at a level chosen from the data: an
occupancy-weighted percentile of χ over the occupied cells (default 50, i.e. the
weighted median). In words, the surface is drawn so that it passes through the
median *wall sample*. Weighting by occ means a cell containing one stray voxel
does not count as loudly as a cell packed with wall. Because the level is
derived from the same field it slices, it remains correct when σ, the grid
spacing, or the input field changes — and it absorbs the undetermined constant
from §2.6 exactly.

Marching cubes then extracts the triangle mesh. Finally, each triangle's
**winding** — the order its three vertices are listed in, which is what defines
which side of the triangle counts as the front — is set from the sign of
(face normal · ∇χ). This uses only the solved field, so it needs neither a
closed surface nor any assumption about the vessel's shape.

A global sign check follows: the median of χ inside the lumen core must exceed
the median outside. If it does not, the field is negated and a warning is
emitted, because an inverted χ would silently flip the sign of every curvature
computed downstream.

---

## 3. The shear solve — `pipeline/shear.py`

### 3.1 What is being computed and why

The quantity wanted is τ_w, the **force per unit area that the moving fluid
exerts on the wall it moves along**, directed along the wall rather than into
it. Its units are force divided by area, i.e. pascals. It is the mechanical
quantity this project's hypothesis is about: whether fibroblasts adhere
differently where the vessel wall experiences a large one
(`docs/DIRECTION.md` §4.3).

It cannot be measured from the images — nothing in this dataset records a flow.
It is computed from the shape of the vessel. §3.3 sets up the equation, and §3.5
explains the part that is easy to disbelieve: why the *relative* answer needs no
flow rate, no pressure and no material constant.

Three assumptions hold throughout §3, and everything after them follows:

1. The flow is **steady** — nothing depends on time.
2. The fluid moves **only along the vessel axis**. There is no sideways motion.
3. The cross-section **does not change** along the length being considered, so
   the speed depends on where you are *across* the section but not on how far
   along you are.

Assumption 3 is the one the vessel violates, since it tapers and branches. §3.9
and §3.8 say where that matters and what is done about it.

### 3.2 Where the cross-section comes from

`lumen_graph.py` builds the vessel centerline on the lumen volume χ > iso
(§2.2), then at stations along it cuts a plane **perpendicular to the local
centerline direction** (`lumen_graph.cross_section`). The result is the open
channel at that station, with a known pixel spacing in µm.

**What is sampled on that plane is χ itself, not the lumen mask** (since
2026-08-28). The plane is filled by linear interpolation of χ, and the section
is where the interpolated χ exceeds the same iso level the surface uses. So each
section carries two things: a 0/1 region, which says *which* channel is being
measured, and the continuous field on the same grid, whose level curve is the
**exact** position of the wall — the same curve the displayed mesh is built
from, because marching cubes interpolates the same χ at the same level.

The earlier version sampled the already-binarised lumen with nearest-neighbour
lookup. That is the same surface, but every wall in the section then landed on
a pixel edge, and a staircase is longer than the curve it approximates: 4.4% on
a test disc at 2 µm spacing and **5.2% at 0.25 µm** — *worse* with resolution,
because it is a property of the representation rather than of the sampling.
On the real sections the perimeter ran about 11% long. See §3.7 and
`docs/ISSUES.md` P46.

Perpendicular cutting matters. An oblique cut through a tube gives a section
wider than the true one, and the true in-plane cross-section is what the shear
depends on. The measured sections are **D-shaped** — not circles, not
ellipses — which is why the code reports extents rather than fitted axes, and
which is why the shear must be solved rather than looked up: every closed-form
shear formula in the literature assumes a circular or rectangular duct.

### 3.3 The equation, derived rather than quoted

Let u be the speed of the fluid along the axis. By assumptions 2 and 3 of §3.1,
u is a scalar field on the flat cross-section: one number per point of a 2D
region, u(y, x). The goal of this subsection is to find the equation u obeys.

**One measured fact about fluids is needed, and this is it.** If two parallel
sheets of fluid slide past each other, each exerts a force on the other along
the direction of sliding, and that force per unit area is proportional to how
fast the speed changes with distance across the sliding direction:

    force per unit area = µ · (du/dn)

where n measures distance perpendicular to the sliding and µ is a constant of
the fluid, in pascal-seconds. Larger µ means a given difference in speed
produces a larger force. This is the definition of µ, and it is the only
physical input; everything below is calculus.

The second fact is not about fluids at all: **pressure**, p, is a force per unit
area pushing on any surface you care to draw, and a difference in pressure
between the two ends of a region produces a net force on the fluid in it.

**Now balance the forces on a small box.** Take a box of side lengths dy, dx and
length ds along the axis, sitting inside the fluid. By assumption 1 the flow is
steady, and by assumption 2 the fluid is not accelerating along the axis, so the
forces along the axis must sum to zero.

Along the axis, the box is pushed by:

- **Pressure on its two end faces.** The pressure is p at one end and p + (dp/ds)·ds
  at the other, each acting on an area dy·dx. The net force along the axis is
  −(dp/ds)·ds·dy·dx.
- **The four side faces**, each dragged by the fluid outside it. Take the two
  faces perpendicular to y. On one, the force per unit area is µ·(∂u/∂y)
  evaluated there; on the other, µ·(∂u/∂y) evaluated a distance dy away, with
  the opposite sign convention. Their sum is the difference of the two, which to
  first order is µ·(∂²u/∂y²)·dy, acting on area ds·dx. The same applies to the
  two faces perpendicular to x, giving µ·(∂²u/∂x²)·dx·ds·dy.

Setting the total to zero and dividing through by the volume ds·dy·dx:

    µ·(∂²u/∂y² + ∂²u/∂x²) = dp/ds

The bracket is the Laplacian of u (§1.4). So

    µ ∇²u = dp/ds

**This is a Poisson equation** (§1.5): the unknown field is u, and the source is
dp/ds ÷ µ, which by assumption 3 is the same number everywhere on the section.
The Laplacian appears not by analogy with anything but because the force on a
sheet of fluid is a *first* derivative of speed, and the net force on a box is
the *difference* of that across the box — a second derivative.

**The boundary condition.** Fluid in contact with a solid does not move relative
to it, so u = 0 everywhere on the wall. That is a Dirichlet condition (§1.6),
and it is what makes the solution unique: the equation says how curved u must
be, and the boundary says where u must vanish.

**Reading the answer.** Once u is known, apply the same measured fact to the
wall itself. At a wall point, let n measure distance straight into the fluid.
Since u = 0 on the wall and u > 0 inside, the direction of steepest increase of
u at the wall *is* that inward direction, so du/dn = |∇u| there. The force per
unit area on the wall is therefore

    τ_w = µ |∇u|    evaluated at the wall.

Three things follow from this that are worth stating plainly:

- τ_w is a **local** quantity. It has a value at every point around the
  perimeter of the section, not one value per section.
- It depends on u only through *how steeply u rises at the wall*, not on how
  fast the fluid moves in the middle.
- A wide, slow flow and a narrow, fast one can produce the same τ_w. The shape
  of the section is what sets the relationship between them, which is why this
  is computable from geometry at all.

### 3.4 How it is solved on the grid

`solve_shape_function` (`shear.py`) solves it with dp/ds = −1 and µ = 1, so
the equation becomes ∇²û = −1 with û = 0 on the wall. The solution û is called
the **shape function** and has units of µm².

The region is the irregular D-shaped mask, not a rectangle, so the cosine-
transform method of §2.6 does not apply — that trick requires a box. Instead the
code builds the linear system explicitly:

- One unknown per pixel inside the mask.
- One equation per pixel: the five-point stencil (§1.4), i.e. the four neighbors
  minus four times the center, set equal to −h².
- These coefficients are stored as a **sparse matrix** — a matrix in which only
  the nonzero entries are stored, appropriate here because each row has at most
  five nonzeros out of thousands of columns — and solved directly with
  `scipy.sparse.linalg.spsolve`.

The one subtlety is **where the wall is**, and it has been got wrong twice.

A pixel at the edge of the mask has neighbors that lie outside it. Setting those
to zero would place u = 0 at the *centers* of the outside pixels, a full pixel
beyond the last interior pixel — that version was implemented first and inflates
Q̂ by 5.8% on a 160 × 55 µm duct at 1 µm spacing, because the flow rate scales
as the cube of the height.

The second version set each missing neighbor to −uᵢ, so the zero lands on the
**face** between them, half a pixel out. That is self-consistent — it is where
marching squares puts a binary contour and where the pixel-count area implies —
but half a pixel is only right when the wall really is at the pixel face, and on
a curved wall it is not: the wall is somewhere between the last inside pixel and
the first outside one, and where exactly varies all the way around the section.

The current version uses the field. Given χ on the plane, the fraction θ of the
way to the outside neighbor at which χ crosses the level is known per pixel and
per direction, so the stencil becomes the **Shortley–Weller** one for unequal
spacing: with s₊ and s₋ the distances to the two neighbors along an axis (each
either a full pixel, or θ·h to the wall),

    u_xx ≈ 2/(s₊+s₋) · [u₊/s₊ + u₋/s₋ − uᵢ(1/s₊ + 1/s₋)]

and u = 0 is imposed on whichever side is the wall. On the analytic disc this
takes Q̂ from 0.29–3.95% error to **0.00%** at every spacing tested.

**This had to move together with the contour.** Measuring the true wall length
while still solving inside the staircase makes τ **25% low** on the disc — worse
than either error alone — because the perimeter and the solve are then
describing different ducts. The half-pixel rule is still available
(`section_from_field=False`) and reproduces a pre-2026-08-28 build exactly.
`docs/scratch/shear/solver_check.py` is what caught all of this.

### 3.5 Why no flow rate, viscosity, or pressure is needed

This is the step that makes the whole thing possible, and it is three lines of
algebra.

The equation of §3.3 is **linear** in u: the unknown appears only to the first
power, never squared or inside another function. So if û solves ∇²û = −1 with
û = 0 on the wall, then for any constant k, the field k·û solves ∇²(k·û) = −k
with the same boundary condition. Every solution on this section is a multiple
of the one already computed. The shape of u is fixed by the section's geometry;
only its scale is open.

Fix that scale by the total flow. The **flow rate** Q is the volume of fluid
crossing the section per unit time, which is the integral of the speed over the
section, Q = ∫u dA. Write

    Q̂ = ∫û dA          (units µm⁴, computed by summing û over the pixels)

If the true solution is u = k·û, then Q = k·Q̂, so k = Q/Q̂. Matching the sources
in ∇²u = (1/µ)(dp/ds) and ∇²(kû) = −k also gives k = −(dp/ds)/µ. Both
expressions for k are the same number, which is the statement that a given
pressure gradient produces a definite flow rate through this particular shape.

Now substitute u = (Q/Q̂)·û into τ_w = µ|∇u|:

    τ_w(s, θ) = µ · Q · |∇û(s, θ)| / Q̂(s)

where s is position along the vessel and θ is position around the section's
perimeter. The pressure gradient has disappeared entirely, and what is left of
the fluid and the flow is the single product µQ.

Along an **unbranched** segment, Q is conserved — the same fluid passes every
section — so **µQ is a single unknown number multiplying the entire profile**.
Every relative statement is therefore fixed by geometry alone: how shear varies
around one section, and how it varies from station to station along a segment.
Only the overall scale is missing. The column `tau_rel_1_um3` carries that
scale-free quantity; `tau_wall_pa` converts to pascals once
someone supplies Q and µ, and neither argument has a default, because in this
project a default is a decision.

**Flow direction is not needed either.** The duct axis comes from the centerline
tangent; the only remaining unknown is which way along it the fluid moves.
Reversing the flow negates u everywhere, and τ_w = µ|∇u| depends only on the
magnitude, so it is unchanged. Direction does matter at a junction, where the
split between branches depends on the branch geometry, and for any
upstream/downstream asymmetry — neither of which arises on a single unbranched
segment.

### 3.6 Getting τ at the wall, and the value of keeping θ

`wall_profile` (`shear.py`) evaluates the wall gradient. The wall is found by
marching squares (§1.7) on the interpolated χ of §3.2, at the same level the
surface uses — so it is a curve passing between pixel centres, not a path along
pixel edges. The distinction is worth the arithmetic: a path that follows pixel
edges around a circle has length 4/π ≈ 1.27 times the circle's, and the
**perimeter P** enters directly into the quantity D_h = 4A/P used elsewhere in
the project.

At each contour point, τ̂ is estimated by fitting û against distance from the
wall. Writing n for the distance inward, a solution of ∇²û = −1 behaves near the
wall as

    û(n) = τ̂·n + ½·û_nn·n²

so a least-squares fit of û against n and n² over a small near-wall
neighborhood, taking the linear coefficient, gives τ̂ directly. The quadratic
term is fitted rather than assumed, because on a curved wall û_nn is not −1 but
−1 + κτ̂ (κ is the wall curvature), and assuming −1 reads about 3% high one pixel
in and worse further in. The neighborhood radius (3 px) was chosen by a sweep
against the analytic disc and rectangle, not by eye: below about 3 px the fit is
starved, and above about 5 px it smooths over a rectangle's corner and collapses
the true variation around the perimeter from 9× to 4.5×.

**n is the distance to the wall *curve*, not to the nearest contour vertex.**
That distinction is worth a paragraph because it was wrong until 2026-08-28 and
it is not a small effect. The contour is a polyline, and its vertices sit
wherever the level crossing happens to fall on a pixel edge — up to one pixel
apart. A pixel sitting 0.1 px from the wall but half a vertex-spacing *along* it
records its distance to the nearest vertex, which can be three to five times its
true distance to the wall. Those near-wall pixels are precisely the ones that
set the linear coefficient of the fit. Handed an exactly correct û for the disc,
the vertex version returned τ̂ **10.8% high**. The distance is now computed
against the contour's line segments, which is straightforward point-to-segment
algebra with nothing to tune.

**Every wall point keeps its own value.** τ varies 6–9× around a single section,
and an adhesion sits at one specific circumferential position, so a single
per-section average would discard exactly the information the analysis needs.
`solve_section` returns the full per-point profile; the per-section scalars are a
summary of it, never a replacement. This is also the project's best handle on the
depth confound: the predicted shear is U-shaped in z while optical attenuation
is monotonic in z, and those two shapes are separable
(`docs/DIRECTION.md` §7.4).

### 3.7 Checks that need no ground truth

- **Force balance.** For the whole cross-section, the total wall friction must
  balance the pressure force driving the fluid: ∮τ dP = A·G exactly, where the
  integral runs around the perimeter, A is the section area and G the pressure
  gradient. With G = 1 this ratio must be 1. `solve_section` reports it for every
  section as `force_balance`, and it requires no reference solution — it would
  catch a case where the solve and the contour disagreed about where the wall
  is.
- **Validation against closed forms.** `docs/scratch/shear/solver_check.py`
  compares the solver against the analytic solution for a disc and against the
  rectangle's double sine series. Two earlier τ estimators were rejected by that
  check before the current one shipped.
- **The bias this document used to end on is fixed** (2026-08-28,
  `docs/ISSUES.md` P46). A pixelated wall is about 6% longer than the smooth
  wall it represents *at every resolution* — a digitization bias, not a
  discretization error — and the solve, being faithful to the duct it is given,
  read τ̂ about 5% low because of it. It is gone now that the section is cut from
  χ (§3.2), the solve stops where χ says (§3.4), and the fit measures distance
  to the wall curve (§3.6). On the analytic disc:

  | pixel spacing | perimeter error | τ̂ error | Q̂ error | force balance |
  |---|---|---|---|---|
  | 2.00 µm | 0.03% | 0.79% | 0.07% | 1.0219 |
  | 1.00 µm | 0.01% | 0.42% | 0.01% | 1.0107 |
  | 0.50 µm | 0.00% | 0.20% | 0.00% | 1.0049 |
  | 0.25 µm | 0.00% | 0.10% | 0.00% | 1.0016 |

  It now **converges** — the τ̂ error halves when the spacing does — where the
  old version got worse as the grid was refined. The part that matters most is
  not in the table: the spurious variation of τ̂ around a shape that has none
  falls from **14–24% to 0.2–0.9%**, and it is the variation around a section
  that §7.4 of `docs/DIRECTION.md` is trying to read.

  On the real sections the force balance above — the check that needs no
  reference answer — moves from a median of **1.041 to 0.9958**. Perimeter fell
  11.2%, hydraulic diameter rose 11.3%, and τ/µQ rose 10.6%. **Every geometric
  number produced before that date is on the old wall** and should not be
  compared with one produced after; `section_from_field=False` reproduces the
  old behaviour exactly if an old artifact has to be rebuilt.

### 3.8 What the flags mean

Neither `lumen_graph` nor `shear` filters anything. A section flagged `clipped`
(the cut ran off the sampling disc or off the edge of the volume) is solved
anyway and the number is returned with the flag beside it. Read the flag: **on a
clipped section, the missing side is treated as a wall**, so û is confined by a
boundary that does not exist and τ̂ there is fiction, not merely uncertain.

There are three flags, they mean different things, and the napari tab gives each
its own control (they shared one until 2026-08-28, which hid how different they
are — on one image 4 of 88 sections are clipped against 23 censored, and no
section is both):

- **`section_clipped`** — the cut ran off the sampling disc or the volume, so
  the section's extent is a lower bound. Not a measurement. Dropped by default.
  In practice these cluster at junctions, where the inscribed radius collapses
  while the cut plane still spans both arms.
- **`censored`** — the section sits within `crop_ratio` local calibers of a
  point where the vessel leaves the field of view. This is a statement about the
  *premise*, not the solve: the force balance there is as good as anywhere, and
  the section may simply not be the whole vessel. Since Q̂ scales roughly as the
  square of the area, absolute τ/µQ out there is untrustworthy while the
  variation *around* the section is not. Shown by default, with a count printed
  every time.
- **near a junction** — hidden within `hide_junction_ratio` local vessel radii
  (default 2.0), because the equation itself does not apply there. Assumptions 2
  and 3 of §3.1 both fail at a branch: the flow divides, which is motion across
  the axis, and the cross-section changes rapidly. See §3.9.

The layer also asks how the branches are put on one scale, which is not a
display detail: by §3.5 each branch carries its **own unknown µQ**, so a single
colour bar across a junction is comparing numbers multiplied by different
constants. §3.9 computes those constants; that is the default. The alternatives
are dividing each branch by its own median, which claims nothing about how they
compare, and the raw scale, which prints a warning saying so.

### 3.9 Along the vessel: resistance, and the split at a branch

§3.5 leaves one unknown number, µQ, per stretch of vessel over which Q is
constant. That is fine on an unbranched vessel and awkward at a branch, where
the flow arriving divides between two outgoing branches and each carries its own
Q. Two branches then have two unknown scales, and no colour, number or plot may
be compared between them.

The division can be computed, and it needs nothing new. Rearrange the two
expressions for k in §3.5:

    dp/ds = −µQ / Q̂(s)

The pressure drop over a stretch of vessel is the integral of that along it. On
a stretch where Q is constant, µ and Q come outside the integral:

    Δp = µQ · ∫ ds / Q̂(s)

The integral depends only on the geometry, since Q̂ was computed at every station
in §3.5. Call it

    R̂ = ∫ ds / Q̂(s)        (units µm⁻³)

so that Δp = µQ·R̂. This says the pressure drop across a branch is proportional
to the flow through it, with a constant fixed by the branch's shape alone.

Two consequences finish the argument.

- **The flow is determined by the pressures at the ends.** Q = Δp/(µR̂).
- **Flow is conserved at a branch point**: what arrives equals what leaves, since
  fluid is not created or destroyed there.

Those two statements are one linear equation per branch and one per junction —
the same structure as any network of linear resistances — so with the pressure
prescribed at each opening of the vessel the whole system is a small set of
simultaneous linear equations, solved exactly. For the simplest case, one branch
dividing into two whose far ends sit at the same pressure, the answer reduces to

    Q₁ / Q₂ = R̂₂ / R̂₁

the flows dividing in inverse proportion to the two resistances.

The outcome is that **one unknown µQ covers the whole vessel** rather than one
per branch. `shear.edge_resistance` computes R̂ by integrating over the stations
that carry a section, and `shear.flow_split` solves the network.

**What is assumed, and it is not nothing.** The pressures at the openings are
not in the image. Three things must be supplied from outside: that the outlets
sit at a common pressure; which opening is the inlet; and whether a branch that
appears to stop inside the imaged volume is an opening or a dead end. The code
takes all three as arguments and the headless runner records the answer under
*every* combination rather than choosing one.

**What this does not give is the shear at the branch point itself.** Assumption
3 of §3.1 — an unchanging cross-section — fails badly there: the flow divides,
which is sideways motion, and that is assumption 2 as well. The equation of §3.3
is simply not the right equation within roughly one vessel width of a junction,
and the sections there are excluded from the display rather than corrected.
Putting the branches on one scale and modelling the junction are different
problems, and only the first is solved.

---

## 4. Are they the same equation?

They share the operator ∇² and the form ∇²f = g. Nothing else.

| | surface solve (§2) | shear solve (§3) |
|---|---|---|
| module | `vessel_poisson.py` | `shear.py` |
| unknown | χ, a constructed indicator field | u, the fluid speed |
| meaning of the unknown | no physical units; only its level set is used | a physical velocity |
| source (right-hand side) | ∇·V, computed from measured wall normals | a constant, the pressure gradient |
| region | the full 3D image box | one 2D cross-section of the lumen |
| dimensionality | 3D, seven-point stencil | 2D, five-point stencil |
| boundary condition | Neumann (zero flux) on the box faces | Dirichlet (u = 0) on the vessel wall |
| solver | discrete cosine transform, direct and exact | sparse matrix, direct via `spsolve` |
| null space | yes — the constant, fixed by the data-derived level | no — Dirichlet pins the solution |
| what is read from the answer | its level set: the wall's position | its gradient at the wall: τ |

Two further points:

- **No shared code.** `shear.py` does not import `vessel_poisson.py`, and
  `vessel_poisson.py` does not import `shear.py`. Each writes its own Laplacian,
  because the two situations need different discretizations and different
  boundary handling — a shared "Poisson helper" would have to be parameterized
  by everything in the table above, at which point it would share nothing.
- **They are chained, not duplicated.** The surface solve produces χ; the lumen
  is χ > iso; `lumen_graph` builds the centerline on that lumen and cuts
  perpendicular sections **out of χ itself**; the shear solve treats the level
  curve of χ on that plane as its Dirichlet wall. The surface displayed in
  napari and the wall the shear is solved against are the same object *and now
  the same representation of it*, so they cannot disagree about where the vessel
  is. Until 2026-08-28 the second solve read the binarised lumen instead — the
  same object, written down more coarsely — and §3.2, §3.4 and §3.6 are what
  that cost. That is deliberate, and it is the same doctrine as
  `denoise.build_mask` and `vessel_mesh.build_mesh`: when two paths need the
  same quantity, they call one implementation rather than keeping two.

---

## 5. Where to look in the code

| Concept | Location |
|---|---|
| coarse occupancy grid | `vessel_poisson.py:166` `occupancy_grid` |
| structure tensor, sheet normals, sheetness | `vessel_poisson.py:210` `sheet_normals` |
| lumen interior by z-enclosure | `vessel_poisson.py:261` `lumen_core` |
| flipping normals inward | `vessel_poisson.py:374` `orient_inward` |
| the 3D Poisson solve (DCT) | `vessel_poisson.py:414` `poisson_dct` |
| checkerboard (Nyquist) detector | `vessel_poisson.py:467` `_checkerboard_ratio` |
| splat, solve, level set, winding, sign check | `vessel_poisson.py:522` `reconstruct` |
| perpendicular cross-section of the lumen | `lumen_graph.py` `cross_section` |
| the 2D Poisson solve (sparse, Shortley-Weller) | `shear.py` `solve_shape_function` |
| sub-pixel wall contour and perimeter | `shear.py` `section_contour` |
| distance from a pixel to the wall curve | `shear.py` `_dist_to_segments` |
| τ̂ at the wall by near-wall fit | `shear.py` `wall_profile` |
| per-section assembly, Q̂, force balance | `shear.py` `solve_section` |
| conversion to pascals | `shear.py` `tau_wall_pa` |
| branch resistance R̂ = ∫ds/Q̂, and the flow split | `shear.py` `edge_resistance`, `flow_split` |
| carrying τ onto the mesh for display | `shear.py` `vertex_field`, `covered_faces` |
| the napari layer | `gui/vessel_mesh.py` `_show_shear` |

Line numbers go stale on every edit; trust the symbol name and grep for it.

---

*Written 2026-08-27; revised 2026-08-28 for `docs/ISSUES.md` P46 (the section
is cut from χ, not from the binarised lumen), P49 (the three section flags, and
why the branches need putting on one scale) and I10 (the flow split from the
branch resistances). Descriptive, not normative: `docs/ISSUES.md` remains
authoritative on what is trusted, and `docs/DIRECTION.md` on what is being
measured and why.*
