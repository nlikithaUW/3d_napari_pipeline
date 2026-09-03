"""Vessel surface as an implicit function fitted to oriented points (Poisson).

**Why this exists.** The default reconstruction is a 3D alpha shape, which
*interpolates*: the surface is forced through every cadherin voxel, so noise
becomes surface and gaps are bridged with flat triangles that stick out into
space. Measured on the reference image (`docs/ISSUES.md` [I06], [P24]): 8.3% of
faces non-manifold, genus 125, 19.2% of the area more than 5 µm from any data,
spurs reaching 13.7 µm, and — the one that matters — **no scale plateau** in the
curvature, `slope_loglog` = −1.19 where a sphere gives +0.02. `curv_*` therefore
measures the fit window rather than the vessel.

This module *approximates* instead. It fits a smooth indicator function χ to the
cadherin sheet and takes a level set of it, so noise averages out, gaps are
spanned by a harmonic solution rather than a fabricated triangle, and branching
plus the measured non-convex cross-section are handled without being told about
them. It is added **alongside** alpha 3D, never in place of it — the head-to-head
has to stay runnable.

The five steps, and the one non-obvious thing about each:

1. **Coarse grid.** Block-mean the mask by ``k`` voxels, one output z-slice at a
   time (the pattern `centerline.lumen_center_z` uses, for the same reason:
   never materialise a float copy of 340M voxels). The block mean is taken over
   the **values**, so ``occ`` is mean wall score in general and the occupancy
   fraction exactly when the input is 0/1 — nothing downstream changes when
   [P28]'s local-contrast field replaces the binary mask, and every number
   recorded on the binary mask still reproduces bit-for-bit. That is what makes
   P28 a controlled experiment.

2. **Sheet normals from the structure tensor, not from the gradient.** The
   cadherin density is a **ridge** — low in the lumen *and* low outside — so
   ``grad(occ)`` flips sign across the wall and is not a normal. The principal
   eigenvector of the smoothed outer product of the gradient is, and it is
   sign-free by construction. It also hands us ``sheetness = (λ3 − λ2)/λ3`` for
   free: an isolated speckle has λ1≈λ2≈λ3 and so contributes ≈ nothing. **That
   is the speckle rejector, and it needs no knob.**

3. **Orientation from the lumen volume, not from a centerline.** Thresholding
   the distance-from-wall field alone selects the *exterior* — the field corners
   sit >100 µm from any wall, further than the lumen's 38.1 µm inscribed radius.
   The discriminator is **z-face contact**: [P24] measured 0.01% of the mask in
   the first z slice and zero in the last, so the vessel is walled above and
   below and every exterior region runs the full z extent, while the lumen does
   not. The vessel exiting through x/y is fine; those faces are allowed.

   This uses no ``binary_fill_holes``, no ``binary_closing`` and no
   ``skeletonize``, so [P29]'s fabricated branch node, [P27]'s raster ordering,
   [P08](c)'s border flare and `lumen_center_z`'s 1.3 µm Z quantization are
   bypassed rather than merely tolerated. It also beats a 1D axis at the Y
   crotch, where the nearest *axis point* can sit across the wall on the other
   branch while the nearest *core voxel* cannot.

4. **Splat and solve, with Neumann boundaries.** χ = 1 inside the lumen and 0
   outside, so ∇χ points inward and the splatted field ``V`` **is** the inward
   normal field with ∇·V = ∇²χ. The solve is DCT-II (Neumann), not Dirichlet:
   the vessel exits through the x/y faces and χ = 0 there would fabricate a cap
   exactly where the pipeline is already weakest. For the same reason the field
   is **not padded** before marching cubes — that pad *is* a Dirichlet cap.

   ``close_faces`` makes that cap available **per face, and only when asked**.
   Its use is the z faces: [P03]'s deep wall is under-detected, so the fitted
   lumen bleeds out through a wall that physically exists, and the resulting
   opening is an artefact of detection rather than a place the vessel leaves.
   Capping asserts wall that was never imaged, which is why it is opt-in,
   recorded in ``info["close_faces"]``, and never the default.

   Two numerics that have to be right. Central differences give a discrete
   Laplacian whose eigenvalue vanishes at the Nyquist mode as well as at DC, so
   an unconstrained checkerboard rides along on χ; the staggered pair (forward
   difference for the gradient, backward for the divergence) does not. And the
   λ = 0 DC mode is the Neumann null space — it is left at zero, which fixes χ's
   mean and is absorbed exactly by the data-derived iso-level.

   Plain Poisson, not screened: the only null mode is a global constant,
   regularisation is already σ_V's job, and a screening term would add a second
   length scale confounded with σ_V — precisely the confound [I05] complains its
   own prototype never controlled.

5. **Surface, then winding read off ∇χ.** The iso-level is an occupancy-weighted
   percentile of χ over the occupied cells — literally "pass through the median
   wall *sample*" — so it stays calibrated when σ, the grid or the input field
   changes, and a cell holding one speckle voxel does not vote as loudly as a
   cell packed with wall (see `_weighted_percentile`). Winding is then set
   per face from the sign of ``face_normal · ∇χ``, which is the real fix to
   [P01] part B: it needs neither a closed surface nor the star-shaped-about-
   centroid assumption that was measured wrong on 44.4% of vertices.

**Purity.** Arrays in, arrays and dicts out. No Qt, no napari, no printing, and
deliberately **no `gpu.*`**: `gpu.gaussian_filter` dispatches to cupy where a GPU
exists, cupy and scipy differ in the last bits, and that would break
`batch_mesh`'s byte-identical rebuild contract. All nine Gaussians run on the
coarse grid and cost ~0.24 s on CPU, so the GPU would buy nothing. Callers print;
this module never does — `diagnose` returns the numbers for them to format.

**The surface is open at the x/y domain faces by design.** That is the Neumann
boundary showing through, not a defect, and it is why `mesh_curvature` must be
handed these normals rather than left to infer a sign from an enclosed volume
(see `cell_mesh.mesh_curvature`'s ``normals=`` argument).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from .cell_mesh import filter_faces_by_extent

# Defaults, mirrored by `vessel_mesh.build_mesh` and `batch_mesh.MESH_DEFAULTS`.
# Four knobs; everything else is derived. They are deliberately *not* merged
# into one: [I05]'s own caveat is that its table never held the grid and the
# level fixed while varying σ, so un-confounding cost from smoothness is part of
# what this method is for.
DEFAULT_GRID_UM = 2.0
DEFAULT_SMOOTH_UM = 5.0
DEFAULT_ISO_PCT = 50.0
# The **lumen detection scale**, split out of `smooth_um` on 2026-09-02
# ([P51](docs/ISSUES.md#p51)). 0 means "follow `smooth_um`", which is what every
# mesh on disk was built under, so the default changes nothing.
#
# Why it had to be split: `smooth_um` was doing three unrelated jobs at once —
# the structure-tensor width in `sheet_normals`, the smoothing of the splatted
# normal field before the solve (the *surface regularisation*), and, in
# `lumen_core`, the scale that decides what counts as lumen at all
# (``deep = D > 2σ``, kept only if ``D >= 4σ``). The first two want a large
# value on a noisy wall; the third wants a small one on a vessel with a narrow
# waist, because 2σ is a **floor on the inscribed radius the method can see**.
# Measured on V3b_04step_003, whose lumen radius runs 1-31 µm: at σ = 5 µm,
# 46.4% of the vessel's length sits below 2σ and is invisible to the trace, the
# surface is severed at x ≈ 40-60 µm and the centerline starts ~90 px inside
# the cadherin; at σ = 2.99 µm it reaches the end but the surface is
# measurably rougher (mean dihedral 3.02° against 2.47°). No single value gave
# both, because the two jobs move in opposite directions across the whole
# usable range.
DEFAULT_CORE_UM = 0.0

# The occupancy-weighted percentile of the smoothed density that defines "wall"
# in `lumen_core`. **Not** the median, which is where this started and which is
# a knife edge: at p50 the band is only the top half of the density ridge, thin
# enough that the lumen leaks into the exterior through it and the whole
# orientation step then refuses (measured on the clean synthetic cylinder, which
# has no gaps at all — the leak was the threshold, not the data). p10 reads
# "the band holding 90% of the wall mass", which is a statement about the data
# rather than a tuned constant, and it errs toward a *thicker* barrier — that
# only shrinks the core, and the 4σ rule below already guards a core that has
# shrunk to nothing.
_WALL_MASS_PCT = 10.0

_EMPTY = (np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64),
          np.zeros((0,)), np.zeros((0, 3)))


# The outermost grid-cell slab of each domain face. Named the same way
# `lumen_graph`'s `face_openings` names them, so an operator reading "the lumen
# leaves mainly through the z faces" can pass exactly that string back in.
_FACE_SLICES = {
    "z_lo": (0, slice(None), slice(None)),
    "z_hi": (-1, slice(None), slice(None)),
    "y_lo": (slice(None), 0, slice(None)),
    "y_hi": (slice(None), -1, slice(None)),
    "x_lo": (slice(None), slice(None), 0),
    "x_hi": (slice(None), slice(None), -1),
}
FACE_NAMES = tuple(_FACE_SLICES)


def _weighted_percentile(values: np.ndarray, weights: np.ndarray,
                         pct: float) -> float:
    """The `pct` percentile of `values`, each weighted by `weights`.

    **Why every percentile in this module is occupancy-weighted.** An unweighted
    percentile over the occupied cells counts a cell holding one speckle voxel
    the same as a cell packed with wall. The synthetic phantom (50,000 isolated
    1–3 voxel blobs, `docs/scratch/vessel_poisson/synthetic_check.py`) makes that
    fatal: speckle cells outnumber wall cells 6:1 there, so the unweighted median
    lands at a speckle density, the "wall" swallows the lumen, and `lumen_core`
    finds nothing at all. Weighting by ``occ`` states what was always meant —
    *the median wall **sample***, not the median occupied *cell* — and a cell
    holding 1 voxel of 216 then carries 1/216 of the vote it used to.
    """
    values = np.asarray(values, dtype=np.float64).ravel()
    weights = np.asarray(weights, dtype=np.float64).ravel()
    if values.size == 0:
        return float("nan")
    order = np.argsort(values, kind="stable")
    v, w = values[order], weights[order]
    total = w.sum()
    if not np.isfinite(total) or total <= 0:
        return float(np.percentile(values, pct))
    c = np.cumsum(w) / total
    return float(v[np.searchsorted(c, float(pct) / 100.0, side="left")
                   .clip(0, len(v) - 1)])


def _voxel_um(scale) -> float:
    """The isotropic µm/voxel. The loader resamples to isotropic, so the three
    agree; the mean is used so a caller passing a marginally-anisotropic scale
    gets a defensible number rather than silently the z one.
    """
    return float(np.mean(np.asarray(scale, dtype=np.float64)))


# --------------------------------------------------------------------------- #
# (1) the coarse occupancy grid
# --------------------------------------------------------------------------- #

def occupancy_grid(mask, scale, grid_um: float) -> tuple[np.ndarray, float, int]:
    """Block-mean `mask` onto a ~`grid_um` grid → ``(occ, h_um, k)``.

    `mask` is either the 0/1 cadherin mask or a continuous wall field
    (`denoise.local_contrast_score`, `docs/ISSUES.md` [P28]). The block mean is
    taken over the **values**, so ``occ`` is mean wall score in general and the
    familiar occupancy *fraction* exactly when the input is 0/1 — the same
    arithmetic on the same 0.0/1.0 summands in the same order, which is what
    keeps every recorded Poisson number reproducible and `batch_mesh`'s
    byte-identical rebuild intact. (A caller handing this a *label* image would
    get block-mean label ids; nothing in the repo does, and `build_mesh` passes
    either a 0/1 mask or a rectified score.)

    ``k = max(1, round(grid_um / voxel))`` and ``h = k · voxel``, so the realised
    spacing is a whole number of voxels and the offset back to µm is exact (see
    `reconstruct` for that correction, which nothing in `mesh_score` can catch).

    Built one output z-slice at a time: the full volume is hundreds of millions
    of voxels and a whole-array mean would materialise a float copy of all of it.
    A trailing partial block on any axis is dropped rather than padded — up to
    ``k − 1`` voxels, ≲2 µm at the default — so every coarse cell is the mean of
    exactly ``k³`` voxels and no edge cell is diluted by invented zeros.
    """
    mask = np.asarray(mask)
    voxel = _voxel_um(scale)
    k = max(1, int(round(float(grid_um) / voxel)))
    h = k * voxel
    nz0, ny0, nx0 = mask.shape
    nz, ny, nx = nz0 // k, ny0 // k, nx0 // k
    if min(nz, ny, nx) < 3:
        return np.zeros((0, 0, 0), np.float32), h, k
    occ = np.empty((nz, ny, nx), np.float32)
    for i in range(nz):
        blk = mask[i * k:(i + 1) * k, :ny * k, :nx * k].astype(np.float32,
                                                               copy=False)
        occ[i] = blk.reshape(k, ny, k, nx, k).mean(axis=(0, 2, 4),
                                                   dtype=np.float32)
    return occ, h, k


# --------------------------------------------------------------------------- #
# (2) sheet normals via the structure tensor
# --------------------------------------------------------------------------- #

def sheet_normals(occ: np.ndarray, h_um: float, smooth_um: float):
    """Unit sheet normal and sheetness at every occupied cell.

    Returns ``(idx, n, sheetness)``: ``idx`` is the ``(M, 3)`` integer index of
    the cells with ``occ > 0``, ``n`` their ``(M, 3)`` unit normals (sign
    arbitrary — `orient_inward` fixes it), and ``sheetness`` the ``(M,)``
    weight ``(λ3 − λ2)/λ3`` in [0, 1].

    The tensor is ``J_ab = G_σ * (∂_a d)(∂_b d)`` with ``d`` the occupancy
    pre-smoothed by one cell. Its **largest** eigenvector is the direction the
    density varies fastest in, i.e. the sheet normal; ``λ3 − λ2`` large relative
    to λ3 means "varies in one direction only", which a sheet does and a blob
    does not. One batched `numpy.linalg.eigh` over the occupied cells does the
    whole field — `eigh`, not a closed form, because it is deterministic and
    LAPACK-stable where a discriminant-based formula is neither.
    """
    idx = np.argwhere(occ > 0)
    if len(idx) == 0:
        return idx, np.zeros((0, 3)), np.zeros((0,))

    # `mode="nearest"` everywhere in this module: the vessel genuinely leaves
    # through the x/y faces, and "constant" would fabricate a density cliff
    # there — the same reason the solve below is Neumann rather than Dirichlet.
    d = ndi.gaussian_filter(occ.astype(np.float64), sigma=1.0, mode="nearest")
    g = [np.gradient(d, float(h_um), axis=a) for a in range(3)]

    sig = float(smooth_um) / float(h_um)
    zi, yi, xi = idx[:, 0], idx[:, 1], idx[:, 2]
    J = np.empty((len(idx), 3, 3), dtype=np.float64)
    for a in range(3):
        for b in range(a, 3):
            Jab = ndi.gaussian_filter(g[a] * g[b], sigma=sig, mode="nearest")
            v = Jab[zi, yi, xi]
            J[:, a, b] = v
            J[:, b, a] = v

    evals, evecs = np.linalg.eigh(J)          # ascending: λ1 ≤ λ2 ≤ λ3
    n = evecs[:, :, 2]
    lam3 = evals[:, 2]
    lam2 = evals[:, 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        sheet = np.where(lam3 > 0, (lam3 - lam2) / lam3, 0.0)
    sheet = np.clip(np.nan_to_num(sheet, nan=0.0), 0.0, 1.0)
    n = n / (np.linalg.norm(n, axis=1, keepdims=True) + 1e-12)
    return idx, n, sheet


# --------------------------------------------------------------------------- #
# (3) the lumen core, and orienting the normals into it
# --------------------------------------------------------------------------- #

def lumen_core(occ: np.ndarray, h_um: float, core_um: float):
    """The lumen interior as a boolean field → ``(core, info)``.

    Wall = where the smoothed occupancy exceeds the level holding
    ``100 − _WALL_MASS_PCT``% of the wall mass (see that constant for why it is
    not the median); ``D`` = distance from wall; candidate cores are the
    components of ``D > 2σ`` **that have wall both above and below them in their
    own z column**. A candidate is kept only if it touches **neither** z face
    and reaches ``D ≥ 4σ`` somewhere.

    **σ here is `core_um`, the lumen detection scale, not the surface smoothing
    scale.** They were one parameter until 2026-09-02 and pull in opposite
    directions on a vessel that both tapers and has a noisy wall — see
    `DEFAULT_CORE_UM` and [P51](docs/ISSUES.md#p51). `reconstruct` passes
    `smooth_um` through when `core_um` is 0 or None, so the split changes
    nothing by default.

    The ``4σ`` rule is a real constraint on the method, not a formality: a
    vessel whose inscribed lumen radius is under ``4 · core_um`` cannot be
    oriented by this route at that σ, and the refusal below says so — with the
    inscribed radius it actually measured, so the reader can see how far off it
    was rather than only that it failed. The reference vessel's measured
    38.1 µm clears the default 20 µm comfortably; V3b_04step_003's 31.2 µm does
    not clear the 60 µm that `poisson_smooth_um = 15` demands, which is what
    [P51] was opened on.

    The z-face rule is the whole discriminator and it is empirical, not
    decorative: [P24] measured 0.01% of the reference mask in the first z slice
    and **zero** in the last, i.e. ~7 µm of headroom above and below the vessel,
    so an exterior region necessarily runs the full z extent while the lumen
    cannot. Two disconnected lumens past the Y junction need no bookkeeping —
    every qualifying component is kept and step (4) takes the distance to their
    union.

    ``info`` carries the per-component volume and max-D that the "spurious core
    pocket" detector reads, plus a note when nothing qualifies. **There is no
    fallback to a centroid rule**: if the vessel is not contained in z the
    orientation is not recoverable from this data, and inventing one would put a
    silent global sign on every curvature downstream.
    """
    info: dict = {"notes": [], "components": [], "core_um": float(core_um)}
    sig = float(core_um) / float(h_um)
    sm = ndi.gaussian_filter(occ.astype(np.float64), sigma=sig, mode="nearest")
    occupied = occ > 0
    empty = np.zeros(occ.shape, bool)
    if not occupied.any():
        info["notes"].append("WARN: the cadherin mask is empty — no lumen core.")
        info["n_core_components"] = 0
        return empty, info

    wall = sm > _weighted_percentile(sm[occupied], occ[occupied],
                                     _WALL_MASS_PCT)
    if not wall.any():
        info["notes"].append("WARN: no wall survived the density threshold — "
                             "no lumen core.")
        info["n_core_components"] = 0
        return empty, info

    D = ndi.distance_transform_edt(~wall, sampling=float(h_um))
    deep = D > 2.0 * float(core_um)

    # **z-enclosure, applied per column rather than by connectivity.** Pure
    # connectivity was tried first and leaks: a hole in the wall band widens by
    # ~σ on each side when the density is thresholded, so a 12 µm gap becomes a
    # ~22 µm channel and the "lumen" merges with the exterior through it — which
    # made the method refuse outright on a phantom with nothing worse than
    # [P28]-scale gaps. Requiring wall both above *and* below a cell in its own
    # z column is the direct statement of what [P24] measured (0.01% of the mask
    # in the first z slice, zero in the last, ~7 µm of headroom), and a lateral
    # leak cannot defeat it because it is not a path property.
    above = np.flip(np.cumsum(np.flip(wall, 0), axis=0), 0) > 0
    below = np.cumsum(wall, axis=0) > 0
    enclosed = above & below

    lab, nlab = ndi.label(deep & enclosed)
    core = np.zeros(occ.shape, bool)
    cell_vol = float(h_um) ** 3
    for lid in range(1, nlab + 1):
        sel = lab == lid
        touches_z = bool(sel[0].any() or sel[-1].any())
        dmax = float(D[sel].max())
        vol = float(sel.sum()) * cell_vol
        keep = (not touches_z) and dmax >= 4.0 * float(core_um)
        info["components"].append({
            "label": int(lid), "volume_um3": vol, "max_D_um": dmax,
            "touches_z_face": touches_z, "kept": bool(keep)})
        if keep:
            core |= sel

    # The exterior, for `reconstruct`'s global sign check: the deep regions that
    # run the full z extent. Labelled *without* the enclosure test, because that
    # is exactly the property the exterior fails.
    olab, onl = ndi.label(deep)
    outside = np.zeros(occ.shape, bool)
    for lid in range(1, onl + 1):
        sel = olab == lid
        if sel[0].any() or sel[-1].any():
            outside |= sel
    info["outside"] = outside

    kept = [c for c in info["components"] if c["kept"]]
    info["n_core_components"] = len(kept)
    info["core_volume_um3"] = float(sum(c["volume_um3"] for c in kept))
    # The largest inscribed radius anywhere z-enclosed, measured whether or not
    # anything qualified. Without it the refusal below says only that the test
    # failed; with it, it says by how much — which is the difference between
    # "this vessel is not contained in z" and "no value of core_um could ever
    # have passed on this image" ([P51] was two hours of the wrong hypothesis
    # for exactly this reason).
    d_best = float(D[enclosed].max()) if enclosed.any() else 0.0
    info["max_enclosed_D_um"] = d_best
    if not kept:
        need = 4.0 * float(core_um)
        info["notes"].append(
            f"WARN: no lumen core found — {nlab} candidate region(s), none of "
            f"which both avoided the z faces and reached {need:.1f} µm from a "
            f"wall. The largest inscribed radius anywhere z-enclosed in this "
            f"mask is {d_best:.1f} µm"
            + (f", which is under the {need:.1f} µm required, so **no z-enclosed "
               f"region could have passed**: lower poisson_core_um to at most "
               f"{d_best / 4.0:.1f} µm (it defaults to poisson_smooth_um)."
               if d_best < need else
               ", which clears the requirement — so this is the z-enclosure "
               "test failing, not the size test: the vessel is probably not "
               "contained in z.")
            + " Refusing to guess an orientation.")
    elif len(kept) > 2:
        # Not alarming on its own: the enclosure test cuts the lumen wherever
        # the wall drops out of a z column, so this counts *fragments*, and is
        # an upper bound on the number of distinct lumens rather than a count of
        # them. Their union is what gets used, which is the same answer either
        # way.
        info["notes"].append(
            f"NOTE: the lumen **core** came out as {len(kept)} fragments, cut "
            f"wherever the wall drops out of a z column. This is a count of "
            f"seed fragments, not of lumens: the core only picks which "
            f"super-level components are lumen and which way the normals "
            f"face. Their union is used, and the lumen itself is labelled "
            f"with plain 3D connectivity and is not fragmented by this — see "
            f"`lumen_graph.lumen_volume`'s `n_components_kept` for the number "
            f"of lumen components, which is the one that matters.")
    return core, info


def orient_inward(n: np.ndarray, idx: np.ndarray, core: np.ndarray,
                  h_um: float):
    """Flip each sheet normal to point toward the lumen core → ``(n, info)``.

    ``Dc`` is the distance to the core, so ``−∇Dc`` points at it; a normal is
    flipped where ``n · ∇Dc > 0``. Distance to the core *volume* rather than to
    a 1D axis is what makes this survive the Y crotch, where the nearest axis
    point can sit across the wall on the other branch and the nearest core cell
    cannot.
    """
    info: dict = {"notes": []}
    if len(idx) == 0 or not core.any():
        info["flipped_fraction"] = None
        return n, info
    Dc = ndi.distance_transform_edt(~core, sampling=float(h_um))
    gz, gy, gx = (np.gradient(Dc, float(h_um), axis=a) for a in range(3))
    zi, yi, xi = idx[:, 0], idx[:, 1], idx[:, 2]
    grad = np.stack([gz[zi, yi, xi], gy[zi, yi, xi], gx[zi, yi, xi]], axis=1)
    dot = np.einsum("ij,ij->i", n, grad)
    flip = dot > 0
    n = n.copy()
    n[flip] = -n[flip]
    info["flipped_fraction"] = float(flip.mean())
    # A cell whose |∇Dc| is ~0 sits at a ridge of the distance field and its
    # orientation is a coin flip. Report the share rather than pretend.
    ambiguous = float((np.linalg.norm(grad, axis=1) < 1e-9).mean())
    info["ambiguous_fraction"] = ambiguous
    if ambiguous > 0.02:
        info["notes"].append(
            f"NOTE: {100 * ambiguous:.1f}% of wall cells sit on a ridge of the "
            f"distance-to-core field, where the inward direction is not "
            f"determined. Their normals are arbitrary but are down-weighted by "
            f"sheetness only if they are also blob-like.")
    return n, info


# --------------------------------------------------------------------------- #
# (4) the Poisson solve
# --------------------------------------------------------------------------- #

def poisson_dct(V: np.ndarray, h_um: float) -> np.ndarray:
    """Solve ``∇²χ = ∇·V`` on a Neumann box, exactly, by DCT-II.

    ``V`` is ``(3, nz, ny, nx)``: a cell-centred vector field read as living on
    the ``i+½`` faces. The operator pair is staggered — **forward** difference
    for the gradient (zero flux at the outer faces), **backward** for the
    divergence — so ``D⁻D⁺`` is exactly the 7-point Neumann Laplacian, whose
    DCT-II eigenvalues are ``λ = −(4/h²) Σ_axes sin²(πk/2N)``. Central
    differences would instead give an operator that also annihilates the Nyquist
    mode, leaving an unconstrained checkerboard riding on χ.

    The λ = 0 mode is the Neumann null space and is left at zero, which fixes
    χ's mean; the data-derived iso-level in `reconstruct` absorbs that exactly.
    ``V`` is zeroed on the last cell of each axis first, which is the discrete
    statement of the zero-flux boundary and makes the system exactly solvable.

    This is a **direct exact** solve of the discrete system, not an iterative
    approximation — there is nothing to fall back to. ``workers=-1`` is bitwise
    identical to ``workers=1`` (pocketfft splits over independent lines), so
    determinism holds.
    """
    from scipy.fft import dctn, idctn

    V = np.asarray(V, dtype=np.float64)
    h = float(h_um)
    shape = V.shape[1:]

    # ∇· by backward difference, with V[-1] = 0 (inflow face) and the outer
    # face flux set to zero (Neumann).
    f = np.zeros(shape, dtype=np.float64)
    for a in range(3):
        v = V[a].copy()
        sl_last = [slice(None)] * 3
        sl_last[a] = -1
        v[tuple(sl_last)] = 0.0
        f += np.diff(v, axis=a, prepend=0.0) / h

    F = dctn(f, type=2, workers=-1)
    lam = np.zeros(shape, dtype=np.float64)
    for a, N in enumerate(shape):
        kk = np.arange(N, dtype=np.float64)
        la = -(4.0 / (h * h)) * np.sin(np.pi * kk / (2.0 * N)) ** 2
        lam = lam + la.reshape([-1 if i == a else 1 for i in range(3)])
    nz = lam != 0
    X = np.zeros(shape, dtype=np.float64)
    X[nz] = F[nz] / lam[nz]
    return idctn(X, type=2, workers=-1)


# --------------------------------------------------------------------------- #
# (5) surface, winding, and the whole thing
# --------------------------------------------------------------------------- #

def _checkerboard_ratio(chi: np.ndarray) -> float:
    """``median |χ − 6-neighbour mean| / std(χ)``. ≪ 1 unless the difference
    operators are mismatched and a Nyquist mode has ridden along.
    """
    s = float(chi.std())
    if not np.isfinite(s) or s <= 0:
        return float("nan")
    kern = np.zeros((3, 3, 3))
    for a in range(3):
        for d in (-1, 1):
            p = [1, 1, 1]
            p[a] += d
            kern[tuple(p)] = 1.0 / 6.0
    nb = ndi.correlate(chi, kern, mode="nearest")
    return float(np.median(np.abs(chi - nb)) / s)


def _wind_faces_by_gradient(verts_idx, faces, chi, h_um):
    """Swap two indices per face so face normals point **down** ∇χ.

    χ is 1 inside the lumen, so "down ∇χ" is away from the lumen — the same
    convention the alpha shape carries (κ2 ≈ −1/R on a tube). Read per face from
    a trilinear sample of ∇χ at its own centroid, so it needs neither a closed
    surface nor [P01](B)'s star-shaped-about-the-centroid assumption, which was
    measured wrong on 44.4% of this vessel's vertices.
    """
    if len(faces) == 0:
        return faces
    grads = [np.gradient(chi, float(h_um), axis=a) for a in range(3)]
    cen = verts_idx[faces].mean(axis=1)          # (M, 3) in grid-index units
    coords = cen.T
    g = np.stack([ndi.map_coordinates(gg, coords, order=1, mode="nearest")
                  for gg in grads], axis=1)
    p = verts_idx[faces]
    fn = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    flip = np.einsum("ij,ij->i", fn, g) > 0
    faces = faces.copy()
    faces[flip] = faces[flip][:, [0, 2, 1]]
    return faces


def _vertex_normals(verts, faces):
    """Per-vertex unit normal = normalised sum of incident face normals. The
    winding is already correct, so no sign heuristic is applied or needed.
    """
    vn = np.zeros_like(verts)
    if len(faces) == 0:
        return vn
    p = verts[faces]
    fn = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    for k in range(3):
        np.add.at(vn, faces[:, k], fn)
    return vn / (np.linalg.norm(vn, axis=1, keepdims=True) + 1e-12)


def reconstruct(mask, scale, *, grid_um: float = DEFAULT_GRID_UM,
                smooth_um: float = DEFAULT_SMOOTH_UM,
                iso_pct: float = DEFAULT_ISO_PCT,
                min_component_um: float = 50.0,
                core_um: float | None = DEFAULT_CORE_UM,
                close_faces: tuple[str, ...] = (),
                extras: dict | None = None):
    """The vessel surface as a level set of the fitted χ.

    ``extras``, when a dict is passed, is filled with the **solved field
    itself** — ``chi`` (the grid), ``iso_level``, ``h_um``, ``origin_um`` (the
    µm coordinate of grid cell ``(0,0,0)``, carrying the same half-cell
    correction the vertices get) and ``core`` (the z-enclosed lumen fragments
    `lumen_core` found). It is filled on the empty-return paths too, with
    whatever had been computed, so a caller can see *why* nothing came back.

    This exists because the lumen is already solved for here and throwing it
    away means recovering it later by a worse route. `lumen_graph` builds the
    3D centerline on ``chi > iso_level``: smooth by construction, gaps spanned
    by the harmonic solve, one component, and needing none of the
    `binary_fill_holes` / `binary_closing` / max-projection machinery whose
    absence is the whole point of this module (`docs/ISSUES.md` §8B). Nothing
    is recomputed for it and no parameter is added — it is the same field the
    surface is a level set of, so the surface and the centerline cannot
    disagree about where the vessel is.

    ``chi`` is a full float64 grid (~54×168×168 at the default 2 µm), so it is
    handed over only when asked for.

    Returns ``(verts_um, faces, values, normals, info)``. ``values`` is the
    distance from the mesh centroid, matching `vessel_mesh.mesh_from_volume` so
    the GUI colormap behaves identically; ``normals`` are the per-vertex outward
    normals implied by ∇χ, which `cell_mesh.mesh_curvature` must be handed (this
    surface is open at the x/y faces, so it encloses no volume and the sign
    cannot be read back off the mesh); ``info`` is the diagnostic record,
    including ``notes`` for the caller to print.

    ``core_um`` is the **lumen detection scale** `lumen_core` runs at; 0 or
    ``None`` means "follow ``smooth_um``", which is what every mesh built before
    2026-09-02 used and is why the split is a no-op by default. Raise
    ``smooth_um`` for a smoother surface, lower ``core_um`` to see a narrower
    lumen — see `DEFAULT_CORE_UM` for the measurement that separated them.

    A deterministic pure function of ``(mask, scale, params)`` — required so
    `batch_mesh`'s byte-identical rebuild check keeps working.

    Because this is marching-cubes output, **every vertex is referenced by a
    face**, so [P09](a)'s interior-vertex defect cannot arise here at all.
    """
    from skimage.measure import marching_cubes

    # 0 / None both mean "follow smooth_um"; resolved once, here, so the
    # recorded value is the one that actually ran rather than the sentinel.
    core_res = float(core_um) if core_um else float(smooth_um)
    info: dict = {"notes": [], "method": "poisson",
                  "grid_um": float(grid_um), "smooth_um": float(smooth_um),
                  "iso_pct": float(iso_pct),
                  "min_component_um": float(min_component_um),
                  "core_um": core_res,
                  "core_um_follows_smooth": not bool(core_um)}
    def _ex(**kw):
        """Record into the caller's out-dict, if there is one."""
        if extras is not None:
            extras.update(kw)

    mask = np.asarray(mask)
    if mask.ndim != 3 or not (mask > 0).any():
        info["notes"].append("WARN: the cadherin mask is empty — no surface.")
        return (*_EMPTY, info)

    voxel = _voxel_um(scale)
    occ, h, k = occupancy_grid(mask, scale, grid_um)
    # Grid cell (i,j,k) sits at i·h + origin. The half-cell term is the same
    # correction the vertices carry below; without it the lumen volume and the
    # surface would sit ~1 µm apart, which is exactly the silent offset the
    # comment on `verts` warns about.
    _ex(h_um=float(h), k_voxels=int(k),
        origin_um=np.full(3, (k - 1) * voxel / 2.0 - h / 2.0))
    info.update({"h_um": float(h), "k_voxels": int(k),
                 "grid_shape": tuple(int(s) for s in occ.shape)})
    if occ.size == 0 or not occ.any():
        info["notes"].append(
            f"WARN: the volume {tuple(mask.shape)} is too small to block-reduce "
            f"by {k}, or is empty on the coarse grid — no surface.")
        return (*_EMPTY, info)
    if smooth_um < h:
        info["notes"].append(
            f"WARN: poisson_smooth_um={smooth_um:g} is below the realised grid "
            f"spacing h={h:.2f} µm, so the smoothing is sub-cell and does "
            f"nothing. Lower poisson_grid_um or raise poisson_smooth_um.")
    if core_res < h:
        info["notes"].append(
            f"WARN: the lumen detection scale {core_res:g} µm is below the "
            f"realised grid spacing h={h:.2f} µm, so `lumen_core` smooths "
            f"sub-cell and 2x/4x it are both under one cell — the detection "
            f"test is no longer resolved. Lower poisson_grid_um to match.")
    info["n_occupied_cells"] = int((occ > 0).sum())

    idx, n, sheet = sheet_normals(occ, h, smooth_um)
    info["sheetness_median"] = float(np.median(sheet)) if len(sheet) else None

    core, core_info = lumen_core(occ, h, core_res)
    # The exterior mask is an array; keep it out of the JSON-serialisable record.
    outside = core_info.pop("outside", np.zeros(occ.shape, bool))
    info["core"] = core_info
    info["notes"].extend(core_info["notes"])
    _ex(core=core)
    if core_info["n_core_components"] == 0:
        # No fallback, by design — see `lumen_core`.
        return (*_EMPTY, info)

    n, orient_info = orient_inward(n, idx, core, h)
    info["orient"] = orient_info
    info["notes"].extend(orient_info["notes"])

    # Splat: weight by occupancy (how much wall is in this cell) and by
    # sheetness (is it wall at all, or a speckle).
    w = occ[idx[:, 0], idx[:, 1], idx[:, 2]].astype(np.float64) * sheet
    V = np.zeros((3, *occ.shape), dtype=np.float64)
    for a in range(3):
        V[a][idx[:, 0], idx[:, 1], idx[:, 2]] = w * n[:, a]
        V[a] = ndi.gaussian_filter(V[a], sigma=float(smooth_um) / h,
                                   mode="nearest")

    chi = poisson_dct(V, h)

    # Global sign check. χ must be *higher* inside the lumen than in the
    # exterior; if the orientation came out inverted the whole curvature field
    # would silently change sign, so this is corrected and warned about, never
    # silently fixed and never silently left.
    # **Median against median, not median against p99.** The p99 form was tried
    # first and is a biased test: the exterior is ~9x the volume of the core and
    # includes every cell just outside the wall, so its 99th percentile sits up
    # inside the wall's own transition and exceeds the core's *median* even when
    # the sign is perfectly right. Measured on the reference image — core median
    # +0.0762, exterior median −0.0037, exterior p99 +0.0985 — so the p99 test
    # negated a correct field and inverted every curvature sign downstream,
    # which is precisely the failure this check exists to catch.
    chi_core = float(np.median(chi[core]))
    chi_out = (float(np.median(chi[outside])) if outside.any()
               else float(np.median(chi)))
    info["chi_core_median"] = chi_core
    info["chi_exterior_median"] = chi_out
    info["chi_exterior_p99"] = (float(np.percentile(chi[outside], 99.0))
                                if outside.any() else None)
    if chi_core < chi_out:
        chi = -chi
        info["sign_flipped"] = True
        info["notes"].append(
            f"WARN: χ came out inverted (median {chi_core:.4g} in the lumen "
            f"core against {chi_out:.4g} outside). The field was negated so "
            f"the surface and its normals are still correct, but the normal "
            f"orientation step disagreed with the lumen volume — check "
            f"`orient.flipped_fraction` before trusting the curvature sign.")
    else:
        info["sign_flipped"] = False
    info["checkerboard_ratio"] = _checkerboard_ratio(chi)
    if np.isfinite(info["checkerboard_ratio"]) \
            and info["checkerboard_ratio"] > 0.5:
        info["notes"].append(
            f"WARN: χ carries a high-frequency component "
            f"(checkerboard ratio {info['checkerboard_ratio']:.2f}, expected "
            f"≪ 1) — the difference operators are not matched.")

    occupied = occ > 0
    level = _weighted_percentile(chi[occupied], occ[occupied], float(iso_pct))
    info["iso_level"] = level

    # ---- optional Dirichlet cap on named domain faces ----------------------
    # The pad this module's docstring calls "a Dirichlet cap" and declines to
    # apply by default, made available per face and **off unless asked**.
    #
    # It is applied HERE — after the sign check and after `level`, before
    # `_ex()` and before marching cubes — for the reason the whole design
    # rests on: the surface, the centerline and every cross-section are level
    # sets of this one array. Capping the mesh afterwards would close the
    # surface while `lumen_graph` went on believing the vessel leaks out of
    # that face.
    #
    # `level` is computed on the UNcapped field, so adding cap cells cannot
    # move the iso-level and a capped run stays comparable to an uncapped one.
    #
    # **This asserts wall that was never imaged.** On the z faces that is the
    # intended reading — [P03]'s deep wall is under-detected, so the fitted
    # lumen bleeds through a wall that is physically there — but it is an
    # assumption and not a measurement, and any shear computed on a section
    # that reaches a capped face is solved against asserted wall. The faces
    # are recorded in `info` so a consumer can tell.
    faces_closed = tuple(str(f) for f in (close_faces or ()))
    if faces_closed:
        bad = [f for f in faces_closed if f not in _FACE_SLICES]
        if bad:
            raise ValueError(
                f"close_faces: unknown face(s) {bad}; expected a subset of "
                f"{tuple(_FACE_SLICES)}")
        # The most-exterior value in the field, so the cap reads as
        # unambiguously outside and the iso-crossing lands within one grid cell
        # of the face. Deterministic, which `batch_mesh`'s byte-identical
        # rebuild depends on.
        cap_value = float(chi.min())
        for f in faces_closed:
            chi[_FACE_SLICES[f]] = cap_value
        info["notes"].append(
            f"NOTE: Dirichlet cap applied on {', '.join(sorted(faces_closed))}"
            f" — χ forced to {cap_value:.4g} on the outermost grid cell of "
            f"each, so the surface closes there. That is asserted wall, not "
            f"measured wall: the vessel is treated as not leaving the volume "
            f"through those faces.")
    info["close_faces"] = list(faces_closed)

    # After the sign check and after any cap, so `chi > iso_level` is the lumen
    # side in the caller's hands too and not a coin flip.
    _ex(chi=chi, iso_level=float(level))
    info["chi_range"] = (float(chi.min()), float(chi.max()))
    if not (chi.min() < level < chi.max()):
        info["notes"].append(
            f"WARN: the iso-level {level:.4g} is outside χ's range "
            f"[{chi.min():.4g}, {chi.max():.4g}] — no surface. Move "
            f"poisson_iso_pct away from its extremes.")
        return (*_EMPTY, info)

    try:
        mc_verts, faces, _nrm, _vals = marching_cubes(
            chi.astype(np.float64), level=level, spacing=(h, h, h))
    except (ValueError, RuntimeError) as exc:
        info["notes"].append(f"WARN: marching cubes found no surface ({exc}).")
        return (*_EMPTY, info)

    faces = np.asarray(faces, dtype=np.int64)
    # Winding is read in grid-index space, where ∇χ is defined.
    verts_idx = mc_verts / h
    faces = _wind_faces_by_gradient(verts_idx, faces, chi, 1.0)

    # The coordinate correction, easy to miss and invisible to every score:
    # block-reduce puts coarse cell i at `i·h + (k−1)·voxel/2`, and treating a
    # cell-centred V as a face value at i+½ shifts the result by h/2. Dropping
    # the −h/2 biases the whole surface ~1 µm outward; nothing in `mesh_score`
    # catches a uniform translation.
    verts = mc_verts + (k - 1) * voxel / 2.0 - h / 2.0

    normals = _vertex_normals(verts, faces)
    faces = filter_faces_by_extent(verts, faces, float(min_component_um))
    values = (np.linalg.norm(verts - verts.mean(0), axis=1) if len(verts)
              else np.zeros((0,)))

    info["n_verts"] = int(len(verts))
    info["n_faces"] = int(len(faces))
    if len(faces) == 0:
        info["notes"].append(
            f"WARN: {len(verts):,} vertices but 0 faces survived "
            f"min_component_um={min_component_um:g} — every component is "
            f"smaller than that.")
    return verts, faces, values, normals, info


def diagnose(mask, scale, *, grid_um: float = DEFAULT_GRID_UM,
             smooth_um: float = DEFAULT_SMOOTH_UM,
             iso_pct: float = DEFAULT_ISO_PCT,
             min_component_um: float = 50.0,
             core_um: float | None = DEFAULT_CORE_UM) -> dict:
    """`reconstruct`'s diagnostic record, plus the mesh counts.

    Separate from `build_mesh` on purpose: `build_mesh` must stay a 4-tuple
    (`saved_check`, `batch_mesh` and the GUI all destructure it), so purity is
    kept and callers that want the numbers ask for them here. This module never
    prints; the caller formats ``notes`` at its own indent.
    """
    verts, faces, _values, normals, info = reconstruct(
        mask, scale, grid_um=grid_um, smooth_um=smooth_um, iso_pct=iso_pct,
        min_component_um=min_component_um, core_um=core_um)
    out = dict(info)
    out["n_verts"] = int(len(verts))
    out["n_faces"] = int(len(faces))
    out["n_normals"] = int(len(normals))
    if len(verts):
        out["extent_um"] = [float(v) for v in
                            (verts.max(0) - verts.min(0))]
    return out
