"""Wall shear stress on a lumen cross-section, from geometry alone.

This is `docs/DIRECTION.md` §4.3 made into code, and item **8a** of §8's queue.

**What it computes and why it needs no flow rate.** For fully developed flow
down a duct of arbitrary cross-section the velocity is purely axial and solves

    µ ∇²u = dp/ds,    u = 0 on the wall,    τ_w = µ |∇u| on the wall

Every closed form in §4 — `4µQ/πR³` for a circle, `6µQ/(w h²)` for a flat duct —
is that one equation solved by hand for a shape we do not have: the measured
sections are **D-shaped, not elliptical** ([P24](b)), which is why
`lumen_graph.cross_section` reports extents rather than fitted semi-axes. So
solve it instead of looking a formula up. It is a 2D Poisson solve on a few
thousand pixels; on image 6 that is milliseconds per section.

Solve once per section with ``dp/ds = -1`` and ``µ = 1``, giving a shape function
``û`` (µm²) and ``Q̂ = ∫û dA`` (µm⁴). Along an **unbranched** segment Q is
conserved (§4.1), so the true gradient at station `s` is ``(dp/ds)(s) = -µQ/Q̂``
and

    τ_w(s, θ) = µ · Q · |∇û(s, θ)| / Q̂(s)

**µQ is a single unknown multiplying the entire profile.** The *shape* of the
shear map along a segment — the whole of the spatial hypothesis — is therefore
fixed by geometry, with no flow rate, no viscosity and no pressure. That shape
is `tau_rel_1_um3` below; `tau_wall_pa` puts absolute units on it once someone
supplies Q and µ.

**Flow direction is not needed.** The duct axis comes from the centerline
tangent; the only unknown is the *sign* along it, and `τ_w = µ|∇u|` is even in
that sign — reverse the flow and every velocity flips while the magnitude of the
wall gradient is unchanged. Direction matters at a **junction** (which branch is
parent, and a flow split that rides on h³, §4.2) and for any upstream/downstream
asymmetry, neither of which arises on a single unbranched segment.

**θ is preserved, deliberately.** A single τ per section is an average around
the perimeter, but an adhesion sits at *one* circumferential position and on an
aspect-2 section the flat sides and the ends differ severalfold. §7.4 turns that
into the best handle the project has on the depth confound: predicted shear is
**U-shaped in z** while optical attenuation is monotonic, and those forms are
separable. `solve_section` therefore returns the full per-boundary-point profile
and `contour_to_3d` maps it back into the image frame; the per-station scalars
are a summary of it, never a replacement.

**Nothing here chooses**, the same doctrine as `lumen_graph` and `mesh_score`.
In particular a section flagged `clipped` — the cut ran off the sampling disc or
off the volume — is solved anyway and the number is returned, because the flag
travels beside it. Read it: on a clipped section the missing side is treated as
a **wall**, so `û` is confined by a boundary that is not real and τ̂ there is
fiction, not merely uncertain.

Pure functions, numpy in and numpy/dict out. No Qt, no napari, no `gpu.*` — the
same reason `vessel_poisson` stays off it: cupy and scipy differ in the last
bits and that would break `batch_mesh`'s byte-identical rebuild.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy import sparse
from scipy.sparse import linalg as spla
from scipy.spatial import cKDTree
from skimage import measure

# --------------------------------------------------------------------------- #
# (1) the solve
# --------------------------------------------------------------------------- #

def solve_shape_function(mask: np.ndarray, step_um: float) -> np.ndarray:
    """``û`` with ``∇²û = -1`` inside ``mask`` and ``û = 0`` on the wall.

    Five-point stencil on the mask pixels. **The wall is taken at the cell
    faces**, half a pixel outside the last interior pixel centre, by mirroring
    each missing neighbour to ``-u_i`` — which is what makes ``û`` zero exactly
    at the face. That is not a detail: it is the same place marching squares
    puts the 0.5 contour and the same place ``area = n·step²`` implies, so the
    solve, the perimeter and the area agree about where the wall is by
    construction rather than by coincidence.

    The obvious alternative — Dirichlet zero at the *centres* of the outside
    pixels — puts the wall a full pixel out instead, and it is wrong by an
    amount that does not vanish: on a 160 × 55 µm duct at 1 µm it inflates the
    effective height by one pixel and `Q̂` by 5.8%, since Q̂ ∝ H³. That version
    was written first and `solver_check.py` is what caught it.

    Returns an array the shape of ``mask``, in µm², zero outside.
    """
    mask = np.asarray(mask, bool)
    h = float(step_um)
    n = int(mask.sum())
    u = np.zeros(mask.shape, float)
    if n == 0:
        return u

    idx = -np.ones(mask.shape, np.int64)
    idx[mask] = np.arange(n)

    rows = [np.arange(n)]
    cols = [np.arange(n)]
    diag = np.full(n, -4.0)
    for axis, shift in ((0, 1), (0, -1), (1, 1), (1, -1)):
        nb = np.roll(idx, shift, axis=axis)
        # A roll wraps; the wrapped edge is outside the duct either way, so
        # blank it rather than letting the far side of the image be a neighbour.
        sl = [slice(None), slice(None)]
        sl[axis] = 0 if shift == 1 else -1
        nb[tuple(sl)] = -1
        j = nb[mask]
        ok = j >= 0
        rows.append(np.arange(n)[ok])
        cols.append(j[ok])
        # Mirror the missing neighbours: u_ghost = -u_i puts zero at the face.
        diag -= ~ok
    vals = [diag] + [np.ones(len(r)) for r in rows[1:]]

    A = sparse.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n)).tocsc()
    x = spla.spsolve(A, np.full(n, -h * h))
    u[mask] = x
    return u


def section_contour(mask: np.ndarray, step_um: float) -> dict:
    """The sub-pixel wall of a section: contour points, arc weights, perimeter.

    Split out so that **the perimeter `lumen_graph.cross_section` reports and
    the wall the shear solve integrates over are the same curve**, not two
    estimates that drift. `D_h = 4A/P` is §4.2's whole input and P is not
    computed anywhere else in this project.

    Marching squares at 0.5, so it is sub-pixel rather than a staircase; each
    vertex owns half the edge on either side of it, so `seg_len_um` sums to the
    perimeter exactly. Multiple contours (a section with a hole, or a mask that
    the plane cut into pieces) are all returned and `n_contours` says so —
    deciding whether that is a defect belongs to the caller.
    """
    mask = np.asarray(mask, bool)
    h = float(step_um)
    out = {"contour_rc": np.zeros((0, 2)), "seg_len_um": np.zeros(0),
           "perim_um": np.nan, "n_contours": 0}
    if not mask.any():
        return out
    contours = measure.find_contours(np.pad(mask.astype(float), 1), 0.5)
    rc, seg = [], []
    for c in contours:
        c = c - 1.0                       # undo the pad
        closed = np.allclose(c[0], c[-1])
        pts = c[:-1] if closed and len(c) > 1 else c
        if len(pts) < 3:
            continue
        nxt = np.roll(pts, -1, axis=0) if closed else np.vstack([pts[1:],
                                                                 pts[-1]])
        prv = np.roll(pts, 1, axis=0) if closed else np.vstack([pts[0],
                                                                pts[:-1]])
        seg.append(0.5 * (np.linalg.norm(nxt - pts, axis=1)
                          + np.linalg.norm(pts - prv, axis=1)) * h)
        rc.append(pts)
    if not rc:
        return out
    out["contour_rc"] = np.concatenate(rc)
    out["seg_len_um"] = np.concatenate(seg)
    out["perim_um"] = float(out["seg_len_um"].sum())
    out["n_contours"] = len(rc)
    return out


# Radius, in pixels, of the near-wall neighbourhood the τ̂ fit runs over, and
# how far in from the wall pixels are taken. **Chosen by a sweep against the
# analytic rectangle and disc, not by eye** — `docs/scratch/shear/solver_check.py`
# reproduces it. Below ~3 px the fit is starved and the disc reads 17% low;
# above ~5 px it smooths the corner of a rectangle and the true max/min around
# the perimeter collapses from 9x to 4.5x, which is exactly the circumferential
# structure §7.4 is trying to see.
DEFAULT_FIT_RADIUS_PX = 3.0
DEFAULT_FIT_DEPTH_PX = 3.0


def wall_profile(u_hat: np.ndarray, mask: np.ndarray, step_um: float,
                 fit_radius_px: float = DEFAULT_FIT_RADIUS_PX,
                 fit_depth_px: float = DEFAULT_FIT_DEPTH_PX) -> dict:
    """``|∇û|`` at sub-pixel points around the wall, plus the perimeter.

    The wall is the marching-squares 0.5 contour of ``mask`` — sub-pixel, and
    the only perimeter estimate here that is not a staircase. (A staircase
    perimeter overestimates a circle by 4/π ≈ 27%, and `D_h = 4A/P` inherits
    every bit of that.)

    **How τ̂ is estimated.** With `n` the distance inward from the wall, a
    solution of ``∇²û = -1`` expands as ``û(n) = τ̂·n + ½û_nn·n²``, so fitting
    ``û`` against ``n`` and ``n²`` over a small near-wall neighbourhood and
    taking the linear coefficient gives τ̂ directly. The quadratic term is
    carried rather than assumed: on a curved wall ``û_nn = -1 + κτ̂``, which for
    a disc of radius 30 at τ̂ = 15 is -0.5, not -1. Fixing it at -1 — i.e. using
    the closed form ``τ̂ = û/n + n/2`` — reads 3% high at one pixel in and drifts
    further with depth. Letting the fit find it removes that entirely.

    Two earlier estimators were tried and are recorded here because
    `solver_check.py` is what killed them, not taste. Sampling û one δ inside
    along the contour normal and dividing reads through bilinear interpolation
    against the zeros stored *outside* the mask, which is not the wall value:
    36% high on the disc, and converging the wrong way in δ. Using only the
    first pixel layer with the closed form is unbiased on a straight wall but
    scatters ±28% around a digitized curve, because that layer's true distance
    to the wall is not the half-pixel the formula assumes.

    **The one bias that does not go away with resolution.** A pixelated wall is
    longer than the smooth wall it stands for — for a digitized disc, ~6% longer
    at every resolution tried, since this is a digitization bias and not a
    discretization error. The solve is faithful to the duct it is *given*, so
    the extra wetted perimeter shows up as τ̂ reading ~5% low on a curved wall.
    It is common-mode: it cancels in every within-segment comparison, which by
    §4.1 is where the science is. It does not cancel in an absolute Pa number.
    """
    mask = np.asarray(mask, bool)
    h = float(step_um)
    out = dict(section_contour(mask, h))
    out["tau_hat_um"] = np.zeros(0)
    pts = out["contour_rc"]
    if len(pts) == 0:
        return out

    # Interior pixels within `fit_depth_px` of the wall, with their **true**
    # distance to the contour rather than a nominal half pixel.
    ctree = cKDTree(pts)
    src = np.argwhere(mask).astype(float)
    d_px, _ = ctree.query(src)
    keep = d_px <= float(fit_depth_px)
    src, d_px = src[keep], d_px[keep]
    u_src = u_hat[mask][keep]
    if len(src) < 4:
        out["tau_hat_um"] = np.full(len(pts), np.nan)
        return out

    tau = np.full(len(pts), np.nan)
    stree = cKDTree(src)
    for i, nb in enumerate(stree.query_ball_point(pts, float(fit_radius_px))):
        if len(nb) < 4:
            continue
        n = d_px[nb] * h
        A = np.column_stack([n, n * n])
        tau[i] = np.linalg.lstsq(A, u_src[nb], rcond=None)[0][0]
    # û > 0 strictly inside and 0 on the wall, so ∂û/∂n at the wall is
    # positive by the maximum principle. A fit that returns ≤ 0 has failed —
    # it is not a small shear — so it is dropped rather than reported. The
    # count travels with it as `n_wall_failed`.
    bad = np.isfinite(tau) & (tau <= 0)
    tau[bad] = np.nan
    out["tau_hat_um"] = tau
    out["n_wall_failed"] = int(bad.sum() + (~np.isfinite(tau)).sum() - bad.sum())
    return out


def solve_section(mask: np.ndarray, step_um: float) -> dict:
    """The whole per-section solve: ``û``, ``Q̂``, and τ̂ around the wall.

    Every returned quantity is the ``dp/ds = -1, µ = 1`` shape function, so
    `tau_hat_um` is in µm and `q_hat_um4` in µm⁴. `tau_rel_1_um3` is the
    combination that survives normalisation along a segment — see the module
    docstring — and is the column an unbranched vessel can be read from with no
    flow rate at all.
    """
    mask = np.asarray(mask, bool)
    h = float(step_um)
    u = solve_shape_function(mask, h)
    prof = wall_profile(u, mask, h)
    area = float(mask.sum()) * h * h
    q_hat = float(u.sum()) * h * h
    w = prof["seg_len_um"]
    tau = prof["tau_hat_um"]
    out = {
        "u_hat_um2": u,
        "area_um2": area,
        "perim_um": prof["perim_um"],
        "hydraulic_diam_um": (4.0 * area / prof["perim_um"]
                              if prof["perim_um"] > 0 else np.nan),
        "q_hat_um4": q_hat,
        "contour_rc": prof["contour_rc"],
        "tau_hat_um": tau,
        "seg_len_um": w,
        "n_contours": prof["n_contours"],
    }
    ok = np.isfinite(tau)
    out["n_wall_points"] = int(ok.sum())
    out["wall_ok_frac"] = (float(ok.sum()) / len(tau)) if len(tau) else np.nan
    if ok.any() and float(w[ok].sum()) > 0:
        # Perimeter-weighted, because the contour vertices are not evenly
        # spaced: marching squares puts more of them on diagonal wall.
        out["tau_hat_mean_um"] = float((tau[ok] * w[ok]).sum() / w[ok].sum())
        out["tau_hat_min_um"] = float(tau[ok].min())
        out["tau_hat_max_um"] = float(tau[ok].max())
        # ∮τ dP = A·G is an exact force balance on the control volume, so with
        # G = 1 this ratio is 1 for a correct solve on a well-resolved wall.
        # It is a **free consistency check** on every section — no ground truth
        # needed — and it is what would catch a mask whose wall the solve and
        # the contour disagree about. On a digitized curve it sits ~1% off for
        # the reason `wall_profile` documents.
        out["force_balance"] = float(
            (tau[ok] * w[ok]).sum() / area) if area > 0 else np.nan
    else:
        out["tau_hat_mean_um"] = np.nan
        out["tau_hat_min_um"] = np.nan
        out["tau_hat_max_um"] = np.nan
        out["force_balance"] = np.nan
    out["tau_rel_1_um3"] = (out["tau_hat_mean_um"] / q_hat
                            if q_hat > 0 else np.nan)
    return out


# --------------------------------------------------------------------------- #
# (2) putting units on it
# --------------------------------------------------------------------------- #

def tau_wall_pa(tau_hat_um, q_hat_um4, q_um3_per_s: float,
                mu_pa_s: float) -> np.ndarray:
    """τ_w in **Pa** from the shape function and a measured flow rate.

    ``τ_w = µ Q |∇û| / Q̂``. Both arguments are required and neither has a
    default: a viscosity default would be a decision, and this project's rule is
    that a default is a decision. Culture medium at 37 °C is ~0.78e-3 Pa·s if
    someone needs a starting number, but it belongs in the caller.

    Q is the flow **through this section**, not the pump's total into the
    device. They are the same number only on a single unbranched channel; past a
    junction the split rides on h³ (§4.2) and needs the branch geometry.
    """
    tau_hat_um = np.asarray(tau_hat_um, float)
    q_hat = np.asarray(q_hat_um4, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(mu_pa_s) * float(q_um3_per_s) * tau_hat_um / q_hat


def tau_wall_pa_from_pressure(hydraulic_diam_um, dp_ds_pa_per_um: float
                              ) -> np.ndarray:
    """§4.2's pressure-driven form: ``τ_w = (Δp/L)·D_h/4``.

    A force balance on the control volume, so it is exact for a uniform duct and
    **independent of Q and of viscosity entirely**. It is a per-section average
    around the perimeter and carries no circumferential distribution; use
    `tau_wall_pa` when θ matters. This is the form to reach for when the rig is
    a gravity head rather than a syringe pump.
    """
    return np.asarray(hydraulic_diam_um, float) * float(dp_ds_pa_per_um) / 4.0


# --------------------------------------------------------------------------- #
# (3) back into the image frame
# --------------------------------------------------------------------------- #

def vertex_field(verts_um: np.ndarray, wall_points_um: np.ndarray,
                 values: np.ndarray, max_dist_um: float) -> np.ndarray:
    """Carry a wall-point measurement onto mesh vertices, by nearest neighbour.

    **Nearest neighbour and nothing else — no interpolation, no smoothing.**
    Every vertex is handed a value that was *measured* at a wall point, so the
    surface a viewer sees is the same set of numbers an adhesion would join to.
    Where the picture looks coarse the fix is to cut sections closer together
    (`section_every_um`), which is a real measurement and costs seconds, not to
    invent values between them — `lumen_graph` leaves NaN between sections for
    the same reason.

    Vertices further than ``max_dist_um`` from any wall point come back **NaN**:
    the reconstructed surface extends past the region sections were cut over,
    and a vertex out there has no shear measurement. Filling it with a median or
    a zero would paint un-measured wall in a colour that reads as a value.

    ``max_dist_um`` has no default on purpose. The right scale is the section
    spacing — a vertex mid-way between two cuts is half a spacing from either —
    and that is the caller's parameter, not this function's to guess.
    """
    verts = np.asarray(verts_um, float)
    pts = np.asarray(wall_points_um, float)
    vals = np.asarray(values, float)
    out = np.full(len(verts), np.nan)
    if len(verts) == 0 or len(pts) == 0:
        return out
    ok = np.isfinite(vals) & np.isfinite(pts).all(1)
    if not ok.any():
        return out
    d, j = cKDTree(pts[ok]).query(verts)
    hit = d <= float(max_dist_um)
    out[hit] = vals[ok][j[hit]]
    return out


def covered_faces(faces: np.ndarray, values: np.ndarray) -> np.ndarray:
    """The faces all of whose vertices carry a measurement.

    For drawing: a surface built from these shows the wall the shear was
    measured on and nothing else, so un-measured wall is absent rather than
    coloured. Curvature fills its NaN vertices with the median, which is fine
    for a field defined everywhere on the mesh and misleading for one that is
    not.
    """
    faces = np.asarray(faces)
    ok = np.isfinite(np.asarray(values, float))
    if len(faces) == 0:
        return faces
    return faces[ok[faces].all(1)]


def contour_to_3d(contour_rc: np.ndarray, extras: dict) -> np.ndarray:
    """Map section contour pixels to zyx µm, using `cross_section`'s basis.

    ``extras`` is the out-dict `lumen_graph.cross_section` fills. This is what
    §7.4 needs: with the wall points in the image frame, an adhesion joins to
    the nearest one, and its **z** — not merely its depth below the coverslip —
    is what separates the U-shaped shear prediction from monotonic attenuation.
    """
    rc = np.asarray(contour_rc, float)
    if rc.size == 0:
        return np.zeros((0, 3))
    a = (rc[:, 0] - extras["center_rc"][0]) * extras["step_um"]
    b = (rc[:, 1] - extras["center_rc"][1]) * extras["step_um"]
    return (np.asarray(extras["point_um"], float)[None, :]
            + a[:, None] * np.asarray(extras["u"], float)[None, :]
            + b[:, None] * np.asarray(extras["v"], float)[None, :])
