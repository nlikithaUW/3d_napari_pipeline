"""Which way do the adhesions point, relative to the vessel?

Pure functions — numpy and scipy only, no Qt, no napari, no `gpu.*`, and no
`_saved` I/O (that belongs to `saved_io.py`). The napari tab
`pipeline/gui/adhesion_orient.py` is the only caller today; anything headless
must call these same functions rather than re-deriving them, which is
CLAUDE.md's standing rule and the reason this file exists at all instead of the
maths living in the tab.

Grown out of `docs/scratch/orientation/`, and it carries that scratch's
conclusions as constraints rather than as prose:

**The stack's z axis is not a neutral direction.** The PSF is elongated in z,
so nearly every adhesion's measured axis leans that way — the median angle
between an adhesion axis and z is about 15° on `04step_003`. Two consequences
are wired into the API:

  * `into_deg` (how far the axis points out of the wall) is **not usable as a
    response variable**. It is ~0.85 rank-correlated with the axis's own angle
    to z, and a further ~0.40 with how z-facing the local wall is, so it mostly
    reports where on the vessel you are standing. It is returned because it is
    needed to know whether an in-wall angle means anything at all, and it is
    named in `LAB_FRAME_CONTAMINATED` so a caller cannot pick it up unaware.
  * anything measured **within the xy plane** is comparatively safe: the blur
    runs along z and imposes no preferred compass direction in xy. That is why
    `local_coherence` and `bin_glyphs` both take a ``plane`` argument, and why
    the xy setting is the one whose result survived the scratch controls.
  * **the wall's own tangent plane is NOT the safe one**, which is the opposite
    of what it looks like. Projecting into it subtracts the normal component,
    and for a z-leaning axis that subtraction removes almost the whole vector,
    leaving a residual whose direction is set by *which way the wall tilts* and
    not by the adhesion. On a tube the wall tilts circumferentially, so this
    manufactures a circumferential preference out of nothing. `lab_axis_control`
    measures exactly how much, and on `04step_003` it manufactures **more than
    the observed effect** (R = 0.61 against the real data's 0.17, both at
    pref ~77 deg). Never quote an `inwall_deg` statistic without it.

**Directions here are axes, not vectors.** An adhesion's major eigenvector has
no sign — "into" and "out of" the wall are the same measurement. So every
average is an eigen-decomposition of an orientation tensor (`mean_axis`), never
a mean of vectors, which would cancel to nothing; and every in-plane angle is
mod 180°, handled by the doubled-angle resultant in `axial_resultant`.
`docs/scratch/orientation/orient.py` got that second point wrong — it doubled a
quantity built from absolute projections, which spans only a quarter turn — and
`inplane.py` records the correction. Do not reintroduce a [0, 90] angle here.

**Nothing in this module gates, filters or decides.** `elongation` and distance
to the wall are returned as columns; which cut to make is the caller's, and the
tab exposes it as a slider for exactly that reason.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .lumen_graph import path_curvature

#: Angles whose value is dominated by the stack's z axis rather than by the
#: vessel. Named so a caller — a colour-by menu, a future CSV column — can mark
#: them rather than presenting them alongside the safe ones.
LAB_FRAME_CONTAMINATED = ("into_deg", "cos_normal", "e1_z_deg", "inwall_deg",
                          "wrap_deg", "cos_axial", "cos_circ")

#: Below this the major axis of an adhesion is not a direction, it is noise.
#: A default, not a rule: every entry point takes it as an argument.
ELONG_MIN_DEFAULT = 1.3


# --------------------------------------------------------------------------- #
# adhesion shape axes
# --------------------------------------------------------------------------- #

def adhesion_axes(labels: np.ndarray, scale_um_zyx) -> dict:
    """Major axis, elongation and centroid of every labelled adhesion, in µm.

    The axis is the major eigenvector of the label's inertia tensor — a second
    *moment*, an integral over the adhesion's voxels, not a derivative of a
    fitted surface. That is why it does not inherit the fine-scale noise that
    makes `curv_*` and `kappa_1_um` untrustworthy ([P24](docs/ISSUES.md#p24)).

    Accumulated with `np.bincount` and diagonalised in one batched `eigh`
    rather than through `skimage.regionprops`: regionprops returns the tensor
    but not its eigen*vectors*, and getting them back means a Python loop over
    every region.

    Returns ``label`` (sorted unique ids, the join key every other artifact in
    the project uses), ``centroid_zyx_um``, ``e1``, ``sigma_um`` (ascending
    semi-axis scales), ``elong`` = sigma3/sigma2, ``flatness`` = sigma1/sigma3
    and ``n_vox``.
    """
    labels = np.asarray(labels)
    flat = labels.ravel()
    idx = np.flatnonzero(flat)
    if idx.size == 0:
        return {"label": np.zeros(0, np.int32),
                "centroid_zyx_um": np.zeros((0, 3)), "e1": np.zeros((0, 3)),
                "sigma_um": np.zeros((0, 3)), "elong": np.zeros(0),
                "flatness": np.zeros(0), "n_vox": np.zeros(0, np.int32)}
    lab = flat[idx].astype(np.int64)
    zz, yy, xx = np.unravel_index(idx, labels.shape)
    P = np.column_stack([zz, yy, xx]).astype(np.float64) * np.asarray(
        scale_um_zyx, float)

    nmax = int(lab.max())
    n = np.bincount(lab, minlength=nmax + 1).astype(np.float64)
    mom1 = np.stack([np.bincount(lab, w, minlength=nmax + 1) for w in P.T],
                    axis=1)
    ids = np.flatnonzero(n > 0)
    cen = mom1[ids] / n[ids, None]

    # Raw second moments then de-meaned, rather than subtracting each label's
    # centroid from its voxels: that would broadcast the centroid back over
    # every voxel, which is the same arithmetic done the slow way.
    C = np.empty((len(ids), 3, 3))
    for i in range(3):
        for j in range(i, 3):
            raw = np.bincount(lab, P[:, i] * P[:, j], minlength=nmax + 1)[ids]
            C[:, i, j] = C[:, j, i] = raw / n[ids] - cen[:, i] * cen[:, j]

    evals, evecs = np.linalg.eigh(C)              # ascending
    sig = np.sqrt(np.maximum(evals, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        elong = np.where(sig[:, 1] > 0, sig[:, 2] / sig[:, 1], np.inf)
        flatness = np.where(sig[:, 2] > 0, sig[:, 0] / sig[:, 2], np.nan)
    return {"label": ids.astype(np.int32), "centroid_zyx_um": cen,
            "e1": evecs[:, :, 2], "sigma_um": sig, "elong": elong,
            "flatness": flatness, "n_vox": n[ids].astype(np.int32)}


# --------------------------------------------------------------------------- #
# the wall-attached frame
# --------------------------------------------------------------------------- #

def centerline_tangents(points_zyx_um, s_um, edge_id, window_um: float):
    """Unit tangent at every centerline point, per edge, via `path_curvature`.

    Calls the shipped `lumen_graph.path_curvature` rather than differencing the
    points here, so the tangent this module resolves angles against and the
    tangent the graph itself is built on are the same quantity. Per edge, because
    a window that reaches across a junction fits into a branch it is not on
    ([P08](docs/ISSUES.md#p08)).
    """
    pts = np.asarray(points_zyx_um, float)
    s = np.asarray(s_um, float)
    eid = np.asarray(edge_id)
    t = np.zeros_like(pts)
    for e in np.unique(eid):
        m = np.flatnonzero(eid == e)
        if len(m) < 2:
            t[m] = np.array([1.0, 0.0, 0.0])
            continue
        o = m[np.argsort(s[m])]
        t[o] = path_curvature(pts[o], window_um)[2]
    return t


def wall_frames(P_um, verts_um, normals, centerline_um, tangents,
                *, degenerate_below: float = 0.2) -> dict:
    """Wall-attached frame ``(n, a, c)`` at each point of ``P_um``.

    ``n`` is the mesh's stored per-vertex normal at the nearest vertex — on a
    Poisson-family mesh that is ∇χ, a first derivative of a smooth field.
    ``a`` is the centerline tangent projected perpendicular to ``n`` (axial),
    and ``c = n × a`` (circumferential).

    Where the tangent is nearly parallel to the wall normal there is no axial
    direction to speak of; those points are **flagged** in ``degenerate`` and
    left in place. Dropping them here would be this module deciding something.
    """
    P = np.asarray(P_um, float)
    V = np.asarray(verts_um, float)
    N = np.asarray(normals, float)
    dmesh, vi = cKDTree(V).query(P)
    n = N[vi]
    n = n / np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)

    C = np.asarray(centerline_um, float)
    dcen, ci = cKDTree(C).query(P)
    t = np.asarray(tangents, float)[ci]

    a = t - np.sum(t * n, axis=1, keepdims=True) * n
    an = np.linalg.norm(a, axis=1, keepdims=True)
    degenerate = an[:, 0] < degenerate_below
    a = a / np.maximum(an, 1e-12)
    return {"n": n, "a": a, "c": np.cross(n, a), "dist_mesh_um": dmesh,
            "dist_cen_um": dcen, "cen_idx": ci, "vert_idx": vi,
            "degenerate": degenerate}


def resolve(e1, n, a, c) -> dict:
    """Direction cosines, the two angles, and the signed in-wall angle.

    ``e1`` is an axis, so the three cosines are absolute values and their
    squares sum to 1.

      * ``wrap_deg``  in [0, 90]: 0 = along the tube, 90 = around it.
      * ``into_deg``  in [0, 90]: 0 = lies in the wall, 90 = normal to it.
        **Lab-frame contaminated** — see the module docstring.
      * ``inwall_deg`` in [0, 180): the same in-wall direction as ``wrap_deg``
        but **signed**, which is what an axial average needs. ``wrap_deg`` folds
        the four quadrants onto one and cannot be averaged as an angle.
      * ``inwall_mag`` in [0, 1]: how much of the axis lies in the wall at all,
        i.e. cos(into_deg). An in-wall angle from a near-zero magnitude is
        noise; it is returned rather than silently NaN'd so the caller can gate.
    """
    e1 = np.asarray(e1, float)
    ca = np.abs(np.sum(e1 * a, axis=1))
    cc = np.abs(np.sum(e1 * c, axis=1))
    cn = np.abs(np.sum(e1 * n, axis=1))
    u = e1 - np.sum(e1 * n, axis=1, keepdims=True) * n
    mag = np.linalg.norm(u, axis=1)
    ok = mag > 1e-9
    inwall = np.full(len(e1), np.nan)
    inwall[ok] = np.degrees(np.arctan2(np.sum(u[ok] * np.asarray(c)[ok], axis=1),
                                       np.sum(u[ok] * np.asarray(a)[ok], axis=1))
                            ) % 180.0
    return {"cos_axial": ca, "cos_circ": cc, "cos_normal": cn,
            "wrap_deg": np.degrees(np.arctan2(cc, np.maximum(ca, 1e-12))),
            "into_deg": np.degrees(np.arcsin(np.clip(cn, 0, 1))),
            "inwall_deg": inwall, "inwall_mag": mag}


def angle_to_z(e1) -> np.ndarray:
    """Angle in degrees between each axis and the stack's z axis, in [0, 90].

    The lab-frame component, carried explicitly. Every result in this module
    should be read next to it — that is what `docs/scratch/orientation` found
    the hard way.
    """
    e1 = np.asarray(e1, float)
    return np.degrees(np.arccos(np.clip(np.abs(e1[:, 0]), 0, 1)))


# --------------------------------------------------------------------------- #
# averaging axes
# --------------------------------------------------------------------------- #

def mean_axis(e1, weights=None):
    """Mean axis and coherence of a set of axes, via the orientation tensor.

    A mean of the vectors is meaningless for axes — two anti-parallel ones
    cancel — so this diagonalises T = Σ w·e1 e1ᵀ / Σ w and takes its major
    eigenvector.

    Coherence is ``(3λ₁ − 1) / 2`` on the normalised eigenvalues: 0 when the
    axes are isotropic in 3D, 1 when they are identical. It is reported
    alongside every mean axis and is not optional — an average direction over a
    bin with no preferred direction is a confident-looking arrow pointing at
    nothing, which is the display failure this whole module exists to avoid.
    """
    e1 = np.asarray(e1, float)
    if len(e1) == 0:
        return np.array([np.nan] * 3), float("nan")
    w = (np.ones(len(e1)) if weights is None
         else np.asarray(weights, float))
    tot = w.sum()
    if not np.isfinite(tot) or tot <= 0:
        return np.array([np.nan] * 3), float("nan")
    T = (e1 * w[:, None]).T @ e1 / tot
    evals, evecs = np.linalg.eigh(T)
    lam1 = float(evals[-1])
    return evecs[:, -1], (3.0 * lam1 - 1.0) / 2.0


def axial_resultant(phi_deg, weights=None):
    """Doubled-angle resultant of a **mod-180** in-plane axial angle.

    ``R`` in [0, 1] is 0 for no preferred axis and 1 for perfect alignment;
    ``pref`` is the preferred angle back in [0, 180).

    The input must span a half turn. Feeding it a [0, 90] quantity — such as
    `resolve`'s ``wrap_deg`` — doubles into half a circle rather than a whole
    one, which inflates ``R`` and makes ``pref`` meaningless;
    `docs/scratch/orientation/orient.py` did exactly that and
    `inplane.py` records the correction. Use ``inwall_deg``.
    """
    phi = np.asarray(phi_deg, float)
    ok = np.isfinite(phi)
    if not ok.any():
        return float("nan"), float("nan")
    z = np.exp(2j * np.radians(phi[ok]))
    w = np.ones(ok.sum()) if weights is None else np.asarray(
        weights, float)[ok]
    if w.sum() <= 0:
        return float("nan"), float("nan")
    m = (z * w).sum() / w.sum()
    return float(np.abs(m)), float(np.degrees(np.angle(m)) / 2.0 % 180.0)


# --------------------------------------------------------------------------- #
# the binned glyph field
# --------------------------------------------------------------------------- #

def _project(e1, plane: str):
    """Axes restricted to ``plane``: "3d" as given, "xy" with z dropped.

    "xy" is the setting whose result survived the scratch controls: the PSF
    elongates along z and imposes no preferred direction within xy, so an
    alignment found there cannot have been put there by the microscope. Axes
    with almost nothing left after the projection are returned as NaN and
    excluded by the caller rather than renormalised up from noise.
    """
    e1 = np.asarray(e1, float)
    if plane == "3d":
        return e1, np.ones(len(e1))
    if plane != "xy":
        raise ValueError(f"plane must be '3d' or 'xy', got {plane!r}")
    v = e1.copy()
    v[:, 0] = 0.0
    mag = np.linalg.norm(v, axis=1)
    out = np.full_like(v, np.nan)
    ok = mag > 1e-9
    out[ok] = v[ok] / mag[ok, None]
    return out, mag


def bin_glyphs(P_um, e1, bin_um: float, *, min_count: int = 5,
               plane: str = "3d", weights=None, axial=None,
               min_plane_mag: float = 0.0) -> dict:
    """One mean axis per cubic spatial bin — the binned-glyph field.

    This is the display that answers the question the per-adhesion arrows
    cannot: *do neighbouring adhesions agree on a direction, and is that
    direction the vessel's?* Neighbour alignment says they agree; the alignment
    is local; and the alignment is not with the vessel axis. Drawing one arrow
    per patch shows all three at once, and shows where each holds.

    ``bin_um`` is the patch size. It is a real choice and it is a slider in the
    tab: the alignment decays over tens of µm, so a bin much larger than that
    averages the structure away and a bin much smaller has too few adhesions to
    average at all. ``min_count`` drops bins with too little to say — those bins
    are **not drawn**, the same rule the wall-shear layer follows for
    un-measured wall ([P49](docs/ISSUES.md#p49)).

    ``axial`` is the local vessel axis at each adhesion (the ``a`` of
    `wall_frames`). When given, each bin also reports ``axial_cos``, the
    absolute cosine between the bin's mean axis and the bin's mean vessel axis
    — the direct test of "do they run along the vessel", per patch.

    Returns arrays over the kept bins: ``center_zyx_um``, ``axis``,
    ``coherence``, ``count``, ``axial_cos``, ``vessel_axis``, plus the
    ``bin_um`` used and the number of bins dropped for being too sparse.

    ``vessel_axis`` is returned and not just consumed because ``axial_cos``
    **cannot be read against a constant**: see `shuffled_axial_cos`.
    """
    P = np.asarray(P_um, float)
    v, pmag = _project(e1, plane)
    keep = np.isfinite(v).all(axis=1) & np.isfinite(P).all(axis=1)
    if min_plane_mag > 0:
        keep &= pmag >= min_plane_mag
    if bin_um <= 0:
        raise ValueError("bin_um must be positive")

    idx = np.floor(P / float(bin_um)).astype(np.int64)
    kept = np.flatnonzero(keep)
    if kept.size == 0:
        return {"center_zyx_um": np.zeros((0, 3)), "axis": np.zeros((0, 3)),
                "coherence": np.zeros(0), "count": np.zeros(0, int),
                "axial_cos": np.zeros(0), "vessel_axis": np.zeros((0, 3)),
                "bin_um": float(bin_um), "n_bins_dropped": 0, "n_used": 0}

    _, inv, cnt = np.unique(idx[kept], axis=0, return_inverse=True,
                            return_counts=True)
    w = (np.ones(len(P)) if weights is None else np.asarray(weights, float))
    ax = None if axial is None else np.asarray(axial, float)

    centers, axes, coh, counts, acos, vaxes = [], [], [], [], [], []
    dropped = 0
    for b in range(len(cnt)):
        sel = kept[inv == b]
        if len(sel) < min_count:
            dropped += 1
            continue
        m, k = mean_axis(v[sel], w[sel])
        centers.append(P[sel].mean(axis=0))
        axes.append(m)
        coh.append(k)
        counts.append(len(sel))
        if ax is None:
            acos.append(np.nan)
            vaxes.append(np.full(3, np.nan))
        else:
            # The bin's vessel direction is itself an axis average, for the same
            # reason the adhesions' is: the tangent's sign flips with the
            # arbitrary direction of travel along an edge.
            va_all, _ = _project(ax[sel], plane)
            va_ok = va_all[np.isfinite(va_all).all(axis=1)]
            va, _ = mean_axis(va_ok)
            vaxes.append(va)
            acos.append(float(np.abs(np.dot(m, va)))
                        if np.isfinite(va).all() else np.nan)
    return {"center_zyx_um": np.array(centers).reshape(-1, 3),
            "axis": np.array(axes).reshape(-1, 3),
            "coherence": np.array(coh, float),
            "count": np.array(counts, int),
            "axial_cos": np.array(acos, float),
            "vessel_axis": np.array(vaxes, float).reshape(-1, 3),
            "bin_um": float(bin_um), "n_bins_dropped": int(dropped),
            "n_used": int(kept.size)}


def lab_axis_control(n, a, c, *, mag_min: float = 0.0, weights=True):
    """What an in-wall statistic reads when the axes carry NO information.

    Substitutes a constant **z axis** for every adhesion and runs the in-wall
    calculation on the real wall frames. A constant axis says nothing about any
    adhesion, so whatever comes back is produced by the wall's geometry alone,
    and any real measurement has to beat it — not beat zero, and not beat a
    shuffle, which cannot see this at all because it preserves each frame's
    relationship to the z axis.

    This is the control that killed the wall-frame in-wall angle as a usable
    quantity on `04step_003`: at ``mag_min`` 0.5 it returns R = 0.607 at
    pref = 77.9 deg where the measured axes give R = 0.170 at 76.0 deg. Same
    preferred direction, larger resultant, from data containing no adhesions.

    Returns ``R``, ``pref_deg`` and ``n`` — the count surviving ``mag_min``,
    which is itself informative: the constant axis lies in the wall only where
    the wall is steeply tilted, so a high ``mag_min`` selects a small and very
    non-random patch of vessel.
    """
    z = np.zeros((len(np.asarray(n)), 3))
    z[:, 0] = 1.0                                  # zyx: index 0 is z
    r = resolve(z, n, a, c)
    sel = np.isfinite(r["inwall_deg"]) & (r["inwall_mag"] >= mag_min)
    if sel.sum() < 3:
        return {"R": float("nan"), "pref_deg": float("nan"),
                "n": int(sel.sum())}
    R, pref = axial_resultant(r["inwall_deg"][sel],
                              r["inwall_mag"][sel] if weights else None)
    return {"R": R, "pref_deg": pref, "n": int(sel.sum())}


def shuffled_axial_cos(axis, vessel_axis, *, reps: int = 400, seed: int = 0):
    """Baseline for ``axial_cos``: the same axes paired with the WRONG patches.

    The obvious comparison — mean |cos| between two independent uniformly random
    axes, 1/2 in 3D and 2/pi in a plane — is **wrong here, and wrong in the
    direction that manufactures a result**. It assumes the two directions being
    compared are each uniform. The vessel's direction is nothing like uniform:
    over one field of view the tube runs one general way, so every patch's
    vessel axis is similar to every other's. Any global lean in the adhesion
    axes then produces a high ``axial_cos`` with no patch-by-patch relationship
    whatsoever.

    Measured on `04step_003`: against the isotropic constant the median
    ``axial_cos`` reads well above it at every bin size, which would say the
    patches follow the vessel — while the per-adhesion permutation test in
    `docs/scratch/orientation/inplane.py` says they do not. This function is the
    difference between those two answers.

    So: keep both sets of axes exactly as measured and shuffle only *which*
    patch's axis meets which patch's vessel direction. Returns the observed
    median and the 2.5/97.5 percentiles of the shuffled median. Inside the band
    means the patches are not aligned with the vessel, however tidy the arrows
    look.
    """
    a = np.asarray(axis, float)
    v = np.asarray(vessel_axis, float)
    ok = np.isfinite(a).all(axis=1) & np.isfinite(v).all(axis=1)
    n = int(ok.sum())
    if n < 3:
        return {"observed": float("nan"), "null": (float("nan"),
                                                   float("nan")), "n": n}
    a, v = a[ok], v[ok]
    obs = float(np.median(np.abs(np.sum(a * v, axis=1))))
    rng = np.random.default_rng(seed)
    meds = [float(np.median(np.abs(np.sum(a * v[rng.permutation(n)], axis=1))))
            for _ in range(reps)]
    lo, hi = np.percentile(meds, [2.5, 97.5])
    return {"observed": obs, "null": (float(lo), float(hi)), "n": n}


def local_coherence(P_um, e1, radius_um: float, *, plane: str = "3d",
                    min_neighbours: int = 2) -> dict:
    """Per-adhesion agreement with its neighbours within ``radius_um``.

    Mean ``|e1_i · e1_j|`` over the neighbours of each adhesion. Unlike
    `bin_glyphs` this needs no bin edges, so it shows the same organisation
    without a binning choice — the two are meant to be read together, and a
    feature that appears in one and not the other is a binning artifact.

    ``NaN`` where an adhesion has fewer than ``min_neighbours`` neighbours;
    those are not drawn rather than drawn as zero.

    The absolute value on the *isotropic* expectation is 0.5 in 3D, but the
    number to compare against is **not** 0.5: with every axis leaning toward z
    the no-structure level is far higher (about 0.64 on `04step_003`). In the
    ``xy`` plane that inflation is gone, which is why xy is the setting the
    scratch result rests on.
    """
    P = np.asarray(P_um, float)
    v, _ = _project(e1, plane)
    ok = np.isfinite(v).all(axis=1) & np.isfinite(P).all(axis=1)
    out = np.full(len(P), np.nan)
    n_nb = np.zeros(len(P), int)
    if ok.sum() < 2:
        return {"coherence": out, "n_neighbours": n_nb}
    sub = np.flatnonzero(ok)
    tree = cKDTree(P[sub])
    for i, nb in enumerate(tree.query_ball_point(P[sub], radius_um)):
        nb = [j for j in nb if j != i]
        n_nb[sub[i]] = len(nb)
        if len(nb) < min_neighbours:
            continue
        out[sub[i]] = float(np.mean(np.abs(v[sub[nb]] @ v[sub[i]])))
    return {"coherence": out, "n_neighbours": n_nb}
