"""In-plane (X-Y) vessel centerline curvature.

The "straight vs curved vs branched" geometry the user cares about is a property
of the vessel *path* as seen from the top-down (X-Y) view — not the local wall
curvature and not anything in Z (the low-resolution, PSF-elongated axis). This
module extracts the vessel centerline from the VE-cadherin footprint and measures
how sharply it bends within the imaging plane.

Why the cadherin mask and not the vessel mesh: the alpha-shape vessel mesh is a
rough point cloud whose X-Y silhouette skeletonizes into spurious branches, so it
cannot tell branched from straight. The cadherin signal gives a smooth, solid
silhouette whose medial axis is a clean centerline (empirically: a branched
vessel yields ~6 junctions and non-zero median curvature, a straight one yields a
single path with median curvature 0).

`lumen_center_z` adds the Z that `vessel_path_curvature_xy` cannot carry, by a
different route — a 3D distance-from-wall field rather than anything projected —
so the path can be drawn where the vessel actually is. It also returns the
inscribed lumen radius, which is the only caliber measurement in the pipeline.

**Two things this module does not do, and callers must not assume.** The points
it returns are in raster order, not path order: they are a point *set*, not a
traced curve (measured, 40% of consecutive points jump more than 5 µm and the
stored order sums to 31,250 µm across a ~300 µm structure). Curvature is fitted
over KD-tree neighbourhoods and so is unaffected, but there is no arc length and
no ordering to parameterise by. And the skeleton is not protected against the
image boundary, so a vessel clipped by the field of view forks there.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

# 8-connected neighbour offsets for skeleton tracing.
_OFF = [(-1, -1), (-1, 0), (-1, 1), (0, -1),
        (0, 1), (1, -1), (1, 0), (1, 1)]


def _prune_spurs(skel: np.ndarray, min_px: float) -> np.ndarray:
    """Iteratively drop skeleton branches shorter than ``min_px`` that end in a
    free endpoint — removes skeletonization hairs while keeping real junctions.
    """
    skel = skel.copy()
    while True:
        ys, xs = np.nonzero(skel)
        S = set(zip(ys.tolist(), xs.tolist()))

        def nbrs(y, x):
            return [(y + dy, x + dx) for dy, dx in _OFF if (y + dy, x + dx) in S]

        endpoints = [(y, x) for (y, x) in S if len(nbrs(y, x)) == 1]
        removed = False
        for e in endpoints:
            path = [e]
            cur, prev = e, None
            while True:
                nxt = [n for n in nbrs(*cur) if n != prev]
                if len(nxt) != 1:
                    break
                prev, cur = cur, nxt[0]
                path.append(cur)
                if len(nbrs(*cur)) != 2:
                    break
            if len(path) < min_px:
                for q in path[:-1]:      # keep the junction pixel itself
                    skel[q] = False
                removed = True
        if not removed:
            return skel


def vessel_path_curvature_xy(
    mask3d: np.ndarray,
    scale_zyx,
    target_px_um: float = 1.0,
    close_um: float = 8.0,
    open_um: float = 4.0,
    prune_um: float = 30.0,
    smooth_um: float = 20.0,
    border_pad_um: float | None = None,
) -> dict:
    """In-plane centerline curvature of a tubular vessel from a 3D mask.

    ``mask3d`` is any 3D labelling/mask of the vessel wall (e.g. the cadherin
    components); non-zero is treated as vessel. ``scale_zyx`` is the voxel size
    in µm. The mask is max-projected to X-Y, solidified, skeletonized to the
    medial axis, spur-pruned, and each centerline point gets the local path
    curvature ``κ = |x'y'' − y'x''| / (x'²+y'²)^{3/2}`` (1/µm) fit over a
    ``smooth_um`` window. Caliber and Z are ignored by construction.

    Returns ``{"pts_yx_um": (M,2), "kappa": (M,), "branch_yx_um": (B,2)}``
    (positions in µm, Y then X to match image axes).
    """
    from scipy.spatial import cKDTree
    from skimage.measure import block_reduce
    from skimage.morphology import (binary_closing, binary_opening, disk,
                                    skeletonize)

    empty = {"pts_yx_um": np.zeros((0, 2), np.float32),
             "kappa": np.zeros(0, np.float32),
             "branch_yx_um": np.zeros((0, 2), np.float32)}
    mask3d = np.asarray(mask3d)
    if mask3d.ndim != 3 or mask3d.size == 0:
        return empty
    py = float(scale_zyx[1])
    px = float(scale_zyx[2])

    foot = np.asarray(mask3d > 0).any(axis=0)
    if not foot.any():
        return empty

    # Downsample to ~target_px_um so the skeleton is smooth and tracing is fast.
    ds = max(1, int(round(target_px_um / max(py, 1e-6))))
    if ds > 1:
        foot = block_reduce(foot, (ds, ds), np.max)
    ppy, ppx = py * ds, px * ds
    ppm = 0.5 * (ppy + ppx)               # ~isotropic px size for radii in px

    foot = ndi.binary_fill_holes(
        binary_closing(foot, disk(max(1, int(round(close_um / ppm))))))
    if open_um > 0:
        foot = binary_opening(foot, disk(max(1, int(round(open_um / ppm)))))
    # Keep the largest connected footprint (drop detached debris).
    lab, n = ndi.label(foot)
    if n > 1:
        biggest = 1 + int(np.argmax(np.bincount(lab.ravel())[1:]))
        foot = lab == biggest
    if not foot.any():
        return empty

    # The vessel is clipped by the field of view, so the footprint's outline
    # contains an artificial straight border edge — and the medial axis of that
    # edge is a branch running along it. Measured on the reference image: the
    # lower arm ran straight down and then kinked hard toward the bottom-right
    # corner, and the upper-left arm forked into two short branches at the top
    # edge. Neither is vessel geometry; both are the border.
    #
    # Replicating the border outward extends every arm straight past the edge,
    # so the medial axis carries on instead of turning, and the pad is cropped
    # away afterwards. The pad has to be about as deep as the region over which
    # the border can influence the axis, which is the footprint's own
    # half-width — so it is read off the footprint's distance transform rather
    # than chosen. (p95 = 81.7 µm on the reference image, and 80 µm is where the
    # lower arm measured straight.)
    if border_pad_um is None:
        d_foot = ndi.distance_transform_edt(foot, sampling=(ppm, ppm))
        border_pad_um = float(np.percentile(d_foot[foot], 95)) if foot.any() \
            else 0.0
    pad = int(round(float(border_pad_um) / ppm))
    if pad > 0:
        skel_full = skeletonize(np.pad(foot, pad, mode="edge"))
        skel = skel_full[pad:-pad, pad:-pad]
    else:
        skel = skeletonize(foot)
    skel = _prune_spurs(skel, prune_um / ppm)
    ys, xs = np.nonzero(skel)
    if len(ys) < 5:
        return empty
    pts = np.column_stack([ys * ppy, xs * ppx]).astype(np.float64)

    # Local in-plane curvature at scale smooth_um (parabola fit in a PCA frame).
    tree = cKDTree(pts)
    kappa = np.zeros(len(pts), np.float32)
    for i, p in enumerate(pts):
        idx = tree.query_ball_point(p, smooth_um)
        if len(idx) < 5:
            _, idx = tree.query(p, k=min(7, len(pts)))
            idx = np.atleast_1d(idx)
        Q = pts[idx] - p
        _, _, vt = np.linalg.svd(Q - Q.mean(0), full_matrices=False)
        a = Q @ vt[0]
        b = Q @ vt[1]
        A = np.column_stack([np.ones_like(a), a, a * a])
        try:
            c, *_ = np.linalg.lstsq(A, b, rcond=None)
            kappa[i] = abs(2.0 * c[2]) / (1.0 + c[1] ** 2) ** 1.5
        except Exception:
            kappa[i] = 0.0

    # Branch nodes: skeleton pixels with >=3 neighbours, collapsed to centroids.
    nb = ndi.convolve(skel.astype(np.int32), np.ones((3, 3), np.int32),
                      mode="constant") - skel.astype(np.int32)
    bmask = skel & (nb >= 3)
    branch = np.zeros((0, 2), np.float32)
    lb, nbn = ndi.label(bmask, structure=np.ones((3, 3)))
    if nbn > 0:
        cen = ndi.center_of_mass(bmask, lb, range(1, nbn + 1))
        branch = (np.asarray(cen) * np.array([ppy, ppx])).astype(np.float32)

    return {"pts_yx_um": pts.astype(np.float32), "kappa": kappa,
            "branch_yx_um": branch}


def lumen_center_z(
    mask3d: np.ndarray,
    scale_zyx,
    pts_yx_um: np.ndarray,
    block: int = 4,
    sigma_um: float = 4.0,
    level_pct: float = 40.0,
    close_um: float = 10.0,
    open_um: float = 4.0,
) -> dict:
    """Z of the **lumen centre** under each xy centerline point (µm).

    `vessel_path_curvature_xy` works on the X-Y projection and so carries no Z.
    Drawing the path at one Z puts most of it off the vessel, so a Z has to come
    from somewhere. The obvious choice — the mean Z of the cadherin under each
    point — does not work here, and the reason is structural: the vessel axis
    lies in X-Y, so a Z column crosses the **near wall and the far wall**, and a
    mean over that bimodal profile only lands in the lumen when both are equally
    detected. They are not: the signal is sparse (measured 2026-08-12, **86% of
    occupied columns detect only one of the two walls**) and attenuates with
    depth. The result flickers between the bottom and top of the stack — spatial
    roughness 2.95 µm, against 0.00 µm for the method below.

    So instead: take a 3D **distance-from-wall** field and, inside the vessel
    footprint, put the centre where that distance is greatest. Three properties
    make this work where the column mean and every enclosure-based method fail:

    - it is **3D**, so a column detecting neither wall still gets a sensible
      distance from walls in neighbouring columns;
    - it needs **no watertight seal**. Flood-filling the lumen was tried and
      returns nothing: the wall sheet is perforated (the alpha-shape surface has
      genus ≈ 125, i.e. ~125 handles), so any fill leaks through. A hole
      perturbs a distance field smoothly instead of breaking it;
    - the imaging volume is barely larger than the vessel — only ~10 µm of
      headroom above the upper wall — so outside the tube the distance stays
      small and the maximum lands in the lumen rather than above it.

    Returns ``{"z_um", "radius_um", "valid", "notes"}``, each array of length
    ``len(pts_yx_um)``. ``radius_um`` is the distance to the wall at the chosen
    point — the **inscribed radius of the lumen**, i.e. the local vessel
    caliber, which nothing else in this pipeline measures. On the reference
    image it comes to a median of **38.1 µm** against the ~40 µm `docs/ISSUES.md`
    [P24](b) derives from the wall's Z extent, which is an independent check.

    ``valid`` is False where the distance maximum runs into the first or last Z
    slice, or the point falls outside the footprint. That is **the vessel
    leaving the field of view**, not a numerical failure — the same clipping
    that forks the 2D skeleton at the image edge — so it is flagged rather than
    hidden or silently clamped.
    """
    from skimage.morphology import binary_closing, binary_opening, disk

    pts = np.asarray(pts_yx_um, dtype=np.float64)
    n = len(pts)
    out = {"z_um": np.zeros(n, np.float32),
           "radius_um": np.zeros(n, np.float32),
           "valid": np.zeros(n, bool), "notes": []}
    mask3d = np.asarray(mask3d)
    if mask3d.ndim != 3 or n == 0:
        out["notes"].append("WARN: no mask or no centerline points — the path "
                            "has no Z and will be drawn at Z = 0.")
        return out

    sz, sy, sx = (float(s) for s in scale_zyx)
    k = max(1, int(block))
    nz0, ny0, nx0 = mask3d.shape
    nz, ny, nx = nz0 // k, ny0 // k, nx0 // k
    if min(nz, ny, nx) < 3:
        out["notes"].append(f"WARN: volume {mask3d.shape} is too small to "
                            f"block-reduce by {k} — Z not computed.")
        return out

    # Local occupancy on a coarse grid, one output slice at a time: the full
    # volume is hundreds of millions of voxels and a whole-array mean over it
    # would materialise a float copy of all of it.
    occ = np.empty((nz, ny, nx), np.float32)
    for i in range(nz):
        blk = mask3d[i * k:(i + 1) * k, :ny * k, :nx * k] > 0
        occ[i] = blk.reshape(k, ny, k, nx, k).mean(axis=(0, 2, 4),
                                                   dtype=np.float32)
    vz, vy, vx = sz * k, sy * k, sx * k
    if not occ.any():
        out["notes"].append("WARN: the cadherin mask is empty — no Z.")
        return out

    sm = ndi.gaussian_filter(occ, sigma=(sigma_um / vz, sigma_um / vy,
                                         sigma_um / vx), mode="constant")
    wall = sm > float(np.percentile(sm[occ > 0], float(level_pct)))
    if not wall.any():
        out["notes"].append("WARN: no wall survived the density threshold — "
                            "no Z.")
        return out
    edt = ndi.distance_transform_edt(~wall, sampling=(vz, vy, vx))

    # The vessel silhouette, solidified — the same operations the X-Y centerline
    # is built from, so the two agree about where the vessel is.
    vm = 0.5 * (vy + vx)
    foot = wall.any(axis=0)
    foot = ndi.binary_fill_holes(
        binary_closing(foot, disk(max(1, int(round(close_um / vm))))))
    if open_um > 0:
        foot = binary_opening(foot, disk(max(1, int(round(open_um / vm)))))
    lab, nlab = ndi.label(foot)
    if nlab > 1:
        foot = lab == 1 + int(np.argmax(np.bincount(lab.ravel())[1:]))

    zarg = np.argmax(edt, axis=0)
    rmax = edt.max(axis=0)
    cy = np.clip((pts[:, 0] / vy).astype(np.int64), 0, ny - 1)
    cx = np.clip((pts[:, 1] / vx).astype(np.int64), 0, nx - 1)
    zi = zarg[cy, cx]

    z = (zi * vz).astype(np.float64)
    valid = foot[cy, cx] & (zi > 0) & (zi < nz - 1)

    # Fill the unbounded points from their nearest bounded neighbour in X-Y
    # rather than leaving them at whatever the argmax escaped to. Two reasons:
    # the raw value there is the top or bottom slice, which is not a lumen
    # centre by any reading; and a caller that hides them leaves **gaps along
    # the centerline**, which is worse than a flagged estimate — the path is
    # continuous in the data, so it should be continuous on screen. `valid`
    # still says which are extrapolated, so nothing is hidden.
    if valid.any() and not valid.all():
        from scipy.spatial import cKDTree
        _, near = cKDTree(pts[valid]).query(pts[~valid], k=1)
        z[~valid] = z[valid][near]

    out["z_um"] = z.astype(np.float32)
    out["radius_um"] = rmax[cy, cx].astype(np.float32)
    out["valid"] = valid

    n_bad = int((~valid).sum())
    if n_bad:
        out["notes"].append(
            f"NOTE: {n_bad}/{n} centerline points ({100.0 * n_bad / n:.0f}%) "
            f"have no bounded lumen above and below — the vessel is leaving the "
            f"field of view there. Their Z is interpolated from the nearest "
            f"bounded point and they are flagged, not dropped.")
    elif not valid.any():
        out["notes"].append(
            "WARN: no centerline point has a bounded lumen above and below — "
            "every Z is unusable. The vessel may not be contained in the stack.")
    return out
