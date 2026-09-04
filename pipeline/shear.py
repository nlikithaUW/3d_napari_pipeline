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
from scipy import sparse
from scipy.sparse import linalg as spla
from scipy.spatial import cKDTree
from skimage import measure

# --------------------------------------------------------------------------- #
# (1) the solve
# --------------------------------------------------------------------------- #

def solve_shape_function(mask: np.ndarray, step_um: float,
                         level_field: np.ndarray | None = None,
                         level: float = 0.5) -> np.ndarray:
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

    **With ``level_field``, the wall goes where the field says instead.** The
    half-pixel rule above is right only when the wall really is at the pixel
    face, and on a curved wall it is not: it is somewhere between the last
    inside pixel and the first outside one. Given the continuous field the mask
    was thresholded from, the crossing fraction θ is known per direction, and
    the stencil becomes the Shortley-Weller one for unequal spacing —
    ``u = 0`` at θ·step rather than at 0.5·step.

    This has to move together with `section_contour`'s ``level_field``. The two
    are one change: measuring the true wall length while solving inside the
    staircase makes τ **worse**, not better (25% low on the disc), because the
    solve and the perimeter would then be describing different ducts. Moving
    both puts every number on the real wall — see [P46](docs/ISSUES.md#p46).

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

    def _neighbour(axis, shift):
        """Index of the neighbour in one direction, -1 where there is none."""
        nb = np.roll(idx, shift, axis=axis)
        # A roll wraps; the wrapped edge is outside the duct either way, so
        # blank it rather than letting the far side of the image be a neighbour.
        sl = [slice(None), slice(None)]
        sl[axis] = 0 if shift == 1 else -1
        nb[tuple(sl)] = -1
        return nb[mask]

    if level_field is None:
        rows = [np.arange(n)]
        cols = [np.arange(n)]
        diag = np.full(n, -4.0)
        for axis, shift in ((0, 1), (0, -1), (1, 1), (1, -1)):
            j = _neighbour(axis, shift)
            ok = j >= 0
            rows.append(np.arange(n)[ok])
            cols.append(j[ok])
            # Mirror the missing neighbours: u_ghost = -u_i puts zero at the
            # face.
            diag -= ~ok
        vals = [diag] + [np.ones(len(r)) for r in rows[1:]]
        A = sparse.coo_matrix(
            (np.concatenate(vals),
             (np.concatenate(rows), np.concatenate(cols))),
            shape=(n, n)).tocsc()
        x = spla.spsolve(A, np.full(n, -h * h))
        u[mask] = x
        return u

    # Shortley-Weller. Per axis, the two spacings are `h` toward an interior
    # neighbour and `θ·h` toward the wall, and
    #     u_xx ≈ 2/(s₊+s₋) · [u₊/s₊ + u₋/s₋ − u_i(1/s₊ + 1/s₋)]
    # with u = 0 on the wall side. θ comes from linear interpolation of the
    # field to `level`, which is the same interpolation marching squares uses
    # for the contour — so the wall the solve stops at and the wall
    # `section_contour` measures are the same curve, not two.
    f = np.asarray(level_field, float) - float(level)
    fi = f[mask]
    rows = [np.arange(n)]
    cols = [np.arange(n)]
    diag = np.zeros(n)
    off_r, off_c, off_v = [], [], []
    # θ below this is treated as this: a pixel centre sitting a thousandth of a
    # cell from the wall would otherwise put 10⁶ on the diagonal and wreck the
    # conditioning for a cell that carries no flow anyway.
    THETA_MIN = 0.02
    for axis in (0, 1):
        jp, jm = _neighbour(axis, -1), _neighbour(axis, 1)
        fp = np.where(jp >= 0, 0.0, _shift_field(f, axis, -1)[mask])
        fm = np.where(jm >= 0, 0.0, _shift_field(f, axis, 1)[mask])
        sp = np.where(jp >= 0, 1.0,
                      np.clip(fi / np.maximum(fi - fp, 1e-12),
                              THETA_MIN, 1.0)) * h
        sm = np.where(jm >= 0, 1.0,
                      np.clip(fi / np.maximum(fi - fm, 1e-12),
                              THETA_MIN, 1.0)) * h
        c = 2.0 / (sp + sm)
        diag -= c * (1.0 / sp + 1.0 / sm)
        for j, sj in ((jp, sp), (jm, sm)):
            ok = j >= 0
            off_r.append(np.arange(n)[ok])
            off_c.append(j[ok])
            off_v.append((c / sj)[ok])
    A = sparse.coo_matrix(
        (np.concatenate([diag] + off_v),
         (np.concatenate(rows + off_r), np.concatenate(cols + off_c))),
        shape=(n, n)).tocsc()
    u[mask] = spla.spsolve(A, np.full(n, -1.0))
    return u


def _shift_field(f: np.ndarray, axis: int, shift: int) -> np.ndarray:
    """`np.roll` with the wrapped edge blanked to a large negative value.

    The wrapped edge is outside the duct by construction, so it must read as
    "well outside" rather than as whatever happens to sit on the far side.
    """
    out = np.roll(f, shift, axis=axis)
    sl = [slice(None), slice(None)]
    sl[axis] = 0 if shift == 1 else -1
    out[tuple(sl)] = -np.abs(f).max() - 1.0
    return out


def section_contour(mask: np.ndarray, step_um: float,
                    level_field: np.ndarray | None = None,
                    level: float = 0.5) -> dict:
    """The sub-pixel wall of a section: contour points, arc weights, perimeter.

    Split out so that **the perimeter `lumen_graph.cross_section` reports and
    the wall the shear solve integrates over are the same curve**, not two
    estimates that drift. `D_h = 4A/P` is §4.2's whole input and P is not
    computed anywhere else in this project.

    Marching squares, so it is sub-pixel rather than a staircase; each vertex
    owns half the edge on either side of it, so `seg_len_um` sums to the
    perimeter exactly. Multiple contours (a section with a hole, or a mask that
    the plane cut into pieces) are all returned and `n_contours` says so —
    deciding whether that is a defect belongs to the caller.

    **``level_field`` is the fix for [P46](docs/ISSUES.md#p46) and it is worth
    understanding why it works.** Run on a 0/1 ``mask``, marching squares can
    only place a vertex on a pixel edge, so the contour it returns is a
    staircase and its length is biased high — 4.4% on a digitized disc at 2 µm,
    **5.2% at 0.25 µm**. It gets *worse* with resolution, because it is a
    digitization bias and not a discretization error, and no amount of finer
    sampling of a binary array removes it.

    Hand it the **continuous** field the mask was thresholded from and the
    contour becomes the real level curve, interpolated between pixels: the same
    disc reads **0.03% at 2 µm and 0.00% at 0.5 µm**. The field costs nothing
    to keep — `lumen_graph.cross_section` samples χ on the cut plane and
    thresholds it, so it has both, and passing the one it already had is the
    whole change.

    ``mask`` still defines *which* region is being measured (the area, the
    solve domain, and which component the caller selected); ``level_field``
    only says where its wall is.
    """
    mask = np.asarray(mask, bool)
    h = float(step_um)
    out = {"contour_rc": np.zeros((0, 2)), "seg_len_um": np.zeros(0),
           "perim_um": np.nan, "n_contours": 0,
           "segments_rc": np.zeros((0, 2, 2))}
    if not mask.any():
        return out
    if level_field is None:
        contours = measure.find_contours(np.pad(mask.astype(float), 1), 0.5)
    else:
        # Below-level padding, so a section touching the array edge closes
        # against the pad rather than running off it — the same role the pad
        # plays in the binary branch.
        f = np.asarray(level_field, float)
        lo = float(min(np.nanmin(f), level)) - 1.0
        f = np.where(np.isfinite(f), f, lo)
        contours = measure.find_contours(
            np.pad(f, 1, constant_values=lo), float(level))
    rc, seg, segments = [], [], []
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
        # The wall as **segments**, not only as vertices. `wall_profile` needs
        # the distance from a pixel to the wall *curve*; the distance to the
        # nearest vertex is a different quantity and it is badly wrong exactly
        # where it matters most — a pixel 0.1 px from the wall but half a
        # vertex-spacing along it reads 3-5x too far, and those near-wall
        # pixels are what set the linear coefficient of the τ fit.
        segments.append(np.stack([pts, nxt], axis=1))
    if not rc:
        return out
    out["contour_rc"] = np.concatenate(rc)
    out["seg_len_um"] = np.concatenate(seg)
    out["perim_um"] = float(out["seg_len_um"].sum())
    out["n_contours"] = len(rc)
    out["segments_rc"] = np.concatenate(segments)
    return out


def _dist_to_segments(pts_rc: np.ndarray, segments_rc: np.ndarray,
                      chunk: int = 4096) -> np.ndarray:
    """Exact distance from each point to the nearest segment, in pixels.

    Straight point-to-segment algebra, chunked over the points so the N x M
    intermediate stays bounded. There is no approximation in here and no
    tolerance to tune: the contour is a polyline and this is its distance
    function.
    """
    P = np.asarray(pts_rc, float)
    S = np.asarray(segments_rc, float)
    if len(P) == 0 or len(S) == 0:
        return np.full(len(P), np.inf)
    A, B = S[:, 0, :], S[:, 1, :]
    AB = B - A
    L2 = np.einsum("ij,ij->i", AB, AB)
    L2 = np.where(L2 > 0, L2, 1.0)
    out = np.empty(len(P))
    for k in range(0, len(P), chunk):
        Q = P[k:k + chunk]
        AP = Q[:, None, :] - A[None, :, :]
        t = np.clip(np.einsum("ijk,jk->ij", AP, AB) / L2, 0.0, 1.0)
        d = AP - t[..., None] * AB[None, :, :]
        out[k:k + chunk] = np.sqrt(np.einsum("ijk,ijk->ij", d, d)).min(1)
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
                 fit_depth_px: float = DEFAULT_FIT_DEPTH_PX,
                 level_field: np.ndarray | None = None,
                 level: float = 0.5) -> dict:
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
    out = dict(section_contour(mask, h, level_field, level))
    out["tau_hat_um"] = np.zeros(0)
    pts = out["contour_rc"]
    if len(pts) == 0:
        return out

    # Interior pixels within `fit_depth_px` of the wall, with their **true**
    # distance to the contour rather than a nominal half pixel. The distance is
    # to the contour *polyline*; using the nearest vertex instead reads 10%
    # high on the analytic disc even when handed an exact û, because the
    # closest pixels — the ones the linear coefficient rests on — are the ones
    # a vertex spacing biases most.
    src = np.argwhere(mask).astype(float)
    ctree = cKDTree(pts)
    d_near, _ = ctree.query(src)                 # cheap prefilter, never used
    near = d_near <= float(fit_depth_px) + 2.0   # as the fitted distance
    src = src[near]
    d_px = _dist_to_segments(src, out["segments_rc"])
    keep = d_px <= float(fit_depth_px)
    src, d_px = src[keep], d_px[keep]
    u_src = u_hat[mask][near][keep]
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


def solve_section(mask: np.ndarray, step_um: float,
                  level_field: np.ndarray | None = None,
                  level: float = 0.5) -> dict:
    """The whole per-section solve: ``û``, ``Q̂``, and τ̂ around the wall.

    Every returned quantity is the ``dp/ds = -1, µ = 1`` shape function, so
    `tau_hat_um` is in µm and `q_hat_um4` in µm⁴. `tau_rel_1_um3` is the
    combination that survives normalisation along a segment — see the module
    docstring — and is the column an unbranched vessel can be read from with no
    flow rate at all.
    """
    mask = np.asarray(mask, bool)
    h = float(step_um)
    u = solve_shape_function(mask, h, level_field, level)
    prof = wall_profile(u, mask, h, level_field=level_field, level=level)
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


# --------------------------------------------------------------------------- #
# (4) along the vessel: resistance, and the flow split at a junction
# --------------------------------------------------------------------------- #
#
# §4.2 says the flow split between two daughters "rides on h³ and is not
# recoverable from geometry". That is true of the *closed form* — you cannot
# read the split off a caliber — but it is not true of the vessel: the per-
# section solve already gives the exact conductance of every station, so the
# resistance of a branch is an integral this project can evaluate, and the
# split follows from Kirchhoff. What is genuinely unrecoverable is the single
# **global** µQ, not the ratios between branches.
#
# This is what retires [P49](docs/ISSUES.md#p49)'s per-edge unknowns: with the
# split in hand, one scale covers the whole vessel and τ/µQ becomes comparable
# across a junction, which per-edge normalisation only ever hid.

def edge_resistance(edge_id, s_um, q_hat_um4, usable=None) -> dict:
    """Hydraulic resistance of each edge, per unit viscosity → ``R̂`` in µm⁻³.

    For fully developed flow ``dp/ds = -µQ/Q̂(s)``, so integrating along an
    edge gives ``ΔP = µQ · ∫ds/Q̂`` and the bracket is a resistance that needs
    no flow rate, no viscosity and no pressure — exactly the same trick that
    makes `tau_rel_1_um3` readable, applied along the axis instead of around
    the wall.

    Trapezoid over the stations that actually carry a section, since the
    sections are cut every `section_every_um` and NaN in between. ``usable``
    (default: every finite station) drops stations the caller does not trust;
    the flags are `lumen_graph`'s to supply, not this function's to invent.

    **Read `span_frac` before using R̂.** The integrand is 1/Q̂, so a *narrow*
    station contributes most, and a narrow station is also the one most likely
    to be dropped as clipped. An edge whose sections span 60% of its length has
    a resistance that is a lower bound, and by an unknown amount — which is why
    the fraction is returned beside the number rather than folded into it.
    """
    edge_id = np.asarray(edge_id, int)
    s = np.asarray(s_um, float)
    q = np.asarray(q_hat_um4, float)
    ok = np.isfinite(q) & (q > 0) & np.isfinite(s)
    if usable is not None:
        ok &= np.asarray(usable, bool)
    out: dict[int, dict] = {}
    for e in np.unique(edge_id):
        on = edge_id == e
        full = float(s[on].max() - s[on].min()) if on.any() else 0.0
        sel = on & ok
        n = int(sel.sum())
        if n < 2:
            out[int(e)] = {"r_hat_1_um3": np.nan, "n_sections": n,
                           "span_um": 0.0, "length_um": full,
                           "span_frac": 0.0}
            continue
        o = np.argsort(s[sel])
        ss, qq = s[sel][o], q[sel][o]
        span = float(ss[-1] - ss[0])
        out[int(e)] = {
            "r_hat_1_um3": float(np.trapezoid(1.0 / qq, ss)),
            "n_sections": n, "span_um": span, "length_um": full,
            "span_frac": (span / full) if full > 0 else np.nan}
    return out


def flow_split(edges, r_hat, inlet_node: int, open_nodes) -> dict:
    """Relative flow in every edge, from the resistances and the topology.

    A Kirchhoff solve: unknown pressures at the interior nodes, ``p = 1`` at
    ``inlet_node`` and ``p = 0`` at every other node in ``open_nodes``,
    conductance ``1/R̂`` on each edge. The absolute pressure drop is arbitrary,
    so the answer is reported as ``q_rel`` **normalised to the inlet edge**:
    the ratios are the content, the scale is not.

    ``open_nodes`` is the caller's decision and it matters. An endpoint that is
    a real vessel terminus is a dead end and carries no flow; one that is a
    field-of-view crop is an opening. `lumen_graph` classifies them
    (`node_endpoint_kind`), but [DIRECTION](docs/DIRECTION.md) §5.1 expects no
    real terminus in this dataset, so an endpoint classified as a terminus is
    usually a lumen reconstructed short — treat it as open unless there is a
    reason not to, and say which was done.

    **Every outlet is held at one pressure.** That is a statement about the rig
    — the daughters drain to a common reservoir — and it is the assumption the
    whole split rests on. If the outlets sit at different pressures the ratios
    change and nothing here would notice.
    """
    edges = [(int(a), int(b)) for a, b in edges]
    r = {int(k): float(v) for k, v in dict(r_hat).items()}
    openset = {int(n) for n in open_nodes}
    inlet = int(inlet_node)
    out = {"q_rel": {}, "p": {}, "notes": [], "inlet_node": inlet,
           "open_nodes": sorted(openset)}
    if inlet not in openset:
        out["notes"].append(
            f"WARN: inlet node {inlet} is not one of the open nodes "
            f"{sorted(openset)} — nothing can flow in.")
        return out

    live = [i for i, (a, b) in enumerate(edges)
            if a >= 0 and b >= 0 and np.isfinite(r.get(i, np.nan))
            and r.get(i, 0.0) > 0]
    if not live:
        out["notes"].append("WARN: no edge has a usable resistance.")
        return out
    nodes = sorted({n for i in live for n in edges[i]})
    free = [n for n in nodes if n not in openset]
    pos = {n: k for k, n in enumerate(free)}

    # Σ_edges (p_a - p_b)/R = 0 at every free node; fixed nodes move to the RHS.
    A = np.zeros((len(free), len(free)))
    b = np.zeros(len(free))
    for i in live:
        na, nb = edges[i]
        g = 1.0 / r[i]
        for u, v in ((na, nb), (nb, na)):
            if u in pos:
                A[pos[u], pos[u]] += g
                if v in pos:
                    A[pos[u], pos[v]] -= g
                else:
                    b[pos[u]] += g * (1.0 if v == inlet else 0.0)
    p = {n: (1.0 if n == inlet else 0.0) for n in openset}
    if free:
        try:
            x = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            out["notes"].append(
                "WARN: the pressure system is singular — the graph has a "
                "component with no opening on it. No split.")
            return out
        p.update({n: float(x[pos[n]]) for n in free})

    q = {}
    for i in live:
        na, nb = edges[i]
        q[i] = (p.get(na, 0.0) - p.get(nb, 0.0)) / r[i]
    inlet_edges = [i for i in live if inlet in edges[i]]
    scale = max((abs(q[i]) for i in inlet_edges), default=0.0)
    if scale <= 0:
        out["notes"].append("WARN: no flow leaves the inlet node.")
        return out
    out["p"] = p
    out["q_rel"] = {i: float(q[i] / scale) for i in live}
    out["q_signed"] = {i: float(q[i] / scale) for i in live}
    # Conservation is exact by construction at every free node; report the
    # worst residual anyway, because it is free and it catches a topology the
    # caller assembled wrong.
    resid = 0.0
    for n in free:
        resid = max(resid, abs(sum(
            (q[i] if edges[i][0] == n else -q[i])
            for i in live if n in edges[i])) / scale)
    out["max_node_residual"] = float(resid)
    return out
