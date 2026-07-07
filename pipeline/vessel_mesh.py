"""Vessel surface mesh built directly from the cadherin signal.

The vessel is a curved Y-tube whose wall is coated with VE-cadherin, but the
signal is noisy and non-contiguous (broken arcs, speckle). Framed per
cross-section perpendicular to the tube axis (the axis lies in xy, z is depth,
so a cross-section is an xz or yz slice), the wall reads as a noisy ring;
"closing" that ring means reconstructing a closed curve from scattered boundary
samples, bridging gaps up to a tunable length scale. Stack the closed slices and
mesh them.

Three closing strategies (compared in the GUI tab):
  1. Per-slice polar fill — cheap, noise-robust, assumes one star-shaped lumen
     per slice (breaks at the Y-branch).
  2. Per-slice 2D alpha fill — handles the branch (two lumens → two loops); one
     `alpha_um` knob = max gap bridged.
  3. 3D alpha on the cadherin point cloud — reuses `cell_mesh.alpha_shape_faces`;
     bridges gaps in 3D and handles branch + curvature natively.

For the two per-slice methods, the GUI runs the fill on xz and yz slicing
independently and combines (AND = conservative, OR = complete).

Everything here is isotropic: the loader already runs `to_isotropic`, so the
incoming `scale` is uniform µm/voxel, marching-cubes spacing is uniform, and
mesh verts come out in µm aligned with the image layers (which use `scale=scale`).

Heavy array ops are kept localized (per-slice loops, mask) so a later `cupy`
swap stays narrow; the CPU path uses only the scipy/skimage already in-repo.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import Delaunay

from . import gpu


# --- point cloud -------------------------------------------------------------

def mask_to_points_um(mask, voxel_zyx_um) -> np.ndarray:
    """Foreground voxels of `mask` as (N, 3) float coordinates in µm."""
    coords = np.argwhere(mask).astype(np.float64)
    return coords * np.asarray(voxel_zyx_um, dtype=np.float64)


def voxel_subsample(points_um, bin_um) -> np.ndarray:
    """Quantize points to a `bin_um` grid and keep one representative per cell.

    3D Delaunay on millions of wall voxels is infeasible; this caps the count
    to ~1e4–5e4 before the alpha shape. `bin_um <= 0` is a no-op.
    """
    if bin_um <= 0 or len(points_um) == 0:
        return points_um
    keys = np.floor(np.asarray(points_um) / float(bin_um)).astype(np.int64)
    _, idx = np.unique(keys, axis=0, return_index=True)
    return points_um[np.sort(idx)]


# --- polar fill (per-slice) --------------------------------------------------

def polar_fill_slice(mask2d, scale2d, n_bins, r_smooth_bins) -> np.ndarray:
    """Close one cross-section as a single star-shaped lumen.

    Boundary pixels → µm; centroid; per-bin robust (median) radius r(θ);
    circularly interpolate empty θ bins so the curve always closes; optional
    circular moving-average smoothing of r(θ); rasterize the filled disk by
    testing each pixel's (θ, r) ≤ r_curve(θ). Vectorized over the slice.
    """
    H, W = mask2d.shape
    out = np.zeros((H, W), dtype=bool)
    coords = np.argwhere(mask2d)
    if len(coords) < 3:
        return out
    scale2d = np.asarray(scale2d, dtype=np.float64)

    pts = coords * scale2d
    centroid = pts.mean(0)
    rel = pts - centroid
    theta = np.arctan2(rel[:, 1], rel[:, 0])
    r = np.hypot(rel[:, 0], rel[:, 1])

    n_bins = int(n_bins)
    edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_idx = np.clip(np.digitize(theta, edges) - 1, 0, n_bins - 1)

    r_curve = np.full(n_bins, np.nan)
    for b in range(n_bins):
        sel = bin_idx == b
        if sel.any():
            r_curve[b] = np.median(r[sel])

    # Circularly interpolate empty bins so the closed curve has no gaps.
    valid = ~np.isnan(r_curve)
    if valid.sum() < 2:
        return out
    vc, vr = centers[valid], r_curve[valid]
    ext_c = np.concatenate([vc - 2 * np.pi, vc, vc + 2 * np.pi])
    ext_r = np.concatenate([vr, vr, vr])
    r_curve = np.interp(centers, ext_c, ext_r)

    # Optional circular moving-average smoothing of r(θ).
    s = int(r_smooth_bins)
    if s > 1:
        kernel = np.ones(s) / s
        r_ext = np.concatenate([r_curve[-s:], r_curve, r_curve[:s]])
        r_curve = np.convolve(r_ext, kernel, mode="same")[s:s + n_bins]

    # Rasterize: every pixel inside its angular bin's radius is filled.
    ys, xs = np.mgrid[0:H, 0:W]
    py = ys * scale2d[0] - centroid[0]
    px = xs * scale2d[1] - centroid[1]
    p_theta = np.arctan2(px, py)
    p_r = np.hypot(py, px)
    p_bin = np.clip(np.digitize(p_theta.ravel(), edges) - 1, 0, n_bins - 1)
    return p_r <= r_curve[p_bin].reshape(H, W)


# --- 2D alpha fill (per-slice) -----------------------------------------------

def _tri_circumradii(pts: np.ndarray, tris: np.ndarray) -> np.ndarray:
    """Circumradius of every 2D triangle (µm). 2D analog of
    `cell_mesh._tet_circumradii`. Degenerate triangles → +inf.
    """
    p = pts[tris]                          # (n, 3, 2)
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    A = np.stack([a, b], axis=1)           # (n, 2, 2)
    rhs = 0.5 * np.stack([
        np.einsum("ij,ij->i", a, a),
        np.einsum("ij,ij->i", b, b),
    ], axis=1)                             # (n, 2)
    det = np.linalg.det(A)
    ok = np.abs(det) > 1e-12
    R = np.full(len(tris), np.inf)
    if ok.any():
        center_rel = np.linalg.solve(A[ok], rhs[ok][..., None])[..., 0]
        R[ok] = np.linalg.norm(center_rel, axis=1)
    return R


def alpha_fill_slice(mask2d, scale2d, alpha_um, fill_lumen=True) -> np.ndarray:
    """Close one cross-section as the 2D alpha shape of its boundary pixels.

    2D Delaunay on the pixels (µm); keep triangles whose circumradius < alpha
    (so `alpha_um` = max gap bridged); rasterize the kept triangles. Their
    union fills each lumen, and disconnected clusters yield multiple loops —
    so the Y-branch produces two lumens naturally.

    With `fill_lumen` (default), `binary_fill_holes` then turns each closed wall
    ring into a solid disk. Without it a small `alpha_um` (< lumen radius) only
    fills the thin wall arcs, leaving a hollow shell that marching cubes wraps
    into stringy ribbons; filling gives solid disks the mesh wraps cleanly.
    """
    H, W = mask2d.shape
    out = np.zeros((H, W), dtype=bool)
    coords = np.argwhere(mask2d)
    if len(coords) < 3:
        return out
    pts_um = coords * np.asarray(scale2d, dtype=np.float64)
    try:
        tri = Delaunay(pts_um)
    except Exception:
        return out
    R = _tri_circumradii(pts_um, tri.simplices)
    keep = tri.simplices[R < float(alpha_um)]
    if len(keep) == 0:
        return out

    from skimage.draw import polygon
    for t in keep:
        rr, cc = polygon(coords[t, 0], coords[t, 1], shape=(H, W))
        out[rr, cc] = True
    if fill_lumen:
        out = ndi.binary_fill_holes(out)
    return out


# --- per-slice volume driver -------------------------------------------------

def _per_slice_volume(mask3d, axis, scale, slice_fn) -> np.ndarray:
    """Run `slice_fn(slice2d, scale2d)` over every cross-section along `axis`.

    xz = fix Y → `mask[:, y, :]` (axes Z, X); yz = fix X → `mask[:, :, x]`
    (axes Z, Y). Empty slices are skipped.
    """
    out = np.zeros_like(mask3d, dtype=bool)
    if axis == "xz":
        scale2d = (scale[0], scale[2])
        for y in range(mask3d.shape[1]):
            sl = mask3d[:, y, :]
            if sl.any():
                out[:, y, :] = slice_fn(sl, scale2d)
    elif axis == "yz":
        scale2d = (scale[0], scale[1])
        for x in range(mask3d.shape[2]):
            sl = mask3d[:, :, x]
            if sl.any():
                out[:, :, x] = slice_fn(sl, scale2d)
    else:
        raise ValueError(f"axis must be 'xz' or 'yz', got {axis!r}")
    return out


def polar_fill_volume(mask3d, axis, scale, n_bins, r_smooth_bins) -> np.ndarray:
    """Stack of `polar_fill_slice` over cross-sections along `axis`."""
    return _per_slice_volume(
        mask3d, axis, scale,
        lambda sl, s2: polar_fill_slice(sl, s2, n_bins, r_smooth_bins),
    )


def alpha_fill_volume(mask3d, axis, scale, alpha_um, fill_lumen=True) -> np.ndarray:
    """Stack of `alpha_fill_slice` over cross-sections along `axis`."""
    return _per_slice_volume(
        mask3d, axis, scale,
        lambda sl, s2: alpha_fill_slice(sl, s2, alpha_um, fill_lumen=fill_lumen),
    )


def close_volume(filled, voxel_um, close_um) -> np.ndarray:
    """Isotropic binary closing (dilate then erode) of a filled volume.

    Bridges small inter-slice gaps left between independently-filled
    cross-sections so the stacked mask meshes into one connected surface.
    `close_um <= 0` is a no-op. The volume is isotropic post-resample, so a
    single box size applies on all axes; uses the `gpu` max/min filters (CPU
    fallback) mirroring `denoise.bernsen_mask`.
    """
    if close_um <= 0:
        return filled
    size = 2 * int(round(float(close_um) / float(voxel_um))) + 1
    if size <= 1:
        return filled
    v = filled.astype(np.float32, copy=False)
    v = gpu.maximum_filter(v, size=size)
    v = gpu.minimum_filter(v, size=size)
    return v > 0.5


# --- combine + mesh ----------------------------------------------------------

def combine(vol_a, vol_b, mode) -> np.ndarray:
    """Combine two filled volumes: AND (conservative) or OR (complete)."""
    if mode == "and":
        return vol_a & vol_b
    if mode == "or":
        return vol_a | vol_b
    raise ValueError(f"mode must be 'and' or 'or', got {mode!r}")


def mesh_from_volume(filled3d, scale):
    """Marching cubes on a filled binary volume → (verts_um, faces, values).

    `values` = distance from the mesh centroid (for colormap), as in
    `cell_mesh`. Empty / saturated volumes (no 0.5 isosurface) return empties.
    """
    empty = (np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64), np.zeros((0,)))
    if not filled3d.any() or filled3d.all():
        return empty
    from skimage.measure import marching_cubes
    spacing = tuple(float(s) for s in scale)
    try:
        verts, faces, _normals, _vals = marching_cubes(
            filled3d.astype(np.float32), level=0.5, spacing=spacing,
        )
    except (ValueError, RuntimeError):
        return empty
    values = np.linalg.norm(verts - verts.mean(0), axis=1)
    return verts, faces, values
