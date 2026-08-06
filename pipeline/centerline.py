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

    skel = _prune_spurs(skeletonize(foot), prune_um / ppm)
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
