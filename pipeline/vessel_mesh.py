"""Vessel surface mesh built directly from the cadherin signal.

The vessel is a curved Y-tube whose wall is coated with VE-cadherin, but the
signal is noisy and non-contiguous (broken arcs, speckle). Framed per
cross-section perpendicular to the tube axis (the axis lies in xy, z is depth,
so a cross-section is an xz or yz slice), the wall reads as a noisy ring;
"closing" that ring means reconstructing a closed curve from scattered boundary
samples, bridging gaps up to a tunable length scale. Stack the closed slices and
mesh them.

Five closing strategies (compared in the GUI tab). The first three
**interpolate** — the surface is forced through the data, so noise becomes
surface; the last two **approximate**, fitting a smooth implicit function to the
same points (`vessel_poisson`, `docs/ISSUES.md` I06):
  1. Per-slice polar fill — cheap, noise-robust, assumes one star-shaped lumen
     per slice (breaks at the Y-branch).
  2. Per-slice 2D alpha fill — handles the branch (two lumens → two loops); one
     `alpha_um` knob = max gap bridged.
  3. 3D alpha on the cadherin point cloud — reuses `cell_mesh.alpha_shape_faces`;
     bridges gaps in 3D and handles branch + curvature natively.
  4. Poisson on oriented points (`vessel_poisson`) — sheet normals from the
     structure tensor, oriented into the lumen volume, then a Neumann DCT solve
     for a smooth indicator function whose level set is the surface. Gaps are
     spanned by a harmonic solution rather than a fabricated triangle.

For the two per-slice methods, the fill runs on xz and yz slicing independently
and combines (AND = conservative, OR = complete).

`build_mesh` at the bottom is the single entry point that dispatches them all:
the napari tab and the headless `batch_mesh` runner both call it, so the two
paths cannot drift.

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

from . import gpu, vessel_poisson
from .cell_mesh import (_tet_circumradii, alpha_shape_faces,
                        filter_faces_by_extent)


# --- point cloud -------------------------------------------------------------

def mask_to_points_um(mask, voxel_zyx_um) -> np.ndarray:
    """Foreground voxels of `mask` as (N, 3) float coordinates in µm."""
    coords = np.argwhere(mask).astype(np.float64)
    return coords * np.asarray(voxel_zyx_um, dtype=np.float64)


class PointCloudTooDense(ValueError):
    """The subsampled cloud is not a surface, so the alpha shape cannot run.

    Raised rather than warned because the failure it replaces is an **unbounded
    hang**, not a wrong number — a warning printed alongside a process that
    never returns is no help (`docs/ISSUES.md` P35).
    """


# `voxel_subsample`'s docstring states the design range: "caps the count to
# ~1e4-5e4". Above that Delaunay slows sharply; the points also sit on the
# 0.325 µm voxel lattice (an original voxel is kept, not a bin centre), and a
# regular lattice is near worst-case for Qhull — massively cospherical.
_POINTS_SOFT = 50_000

# The refusal is on **bin occupancy**, not on the count, because that is the
# quantity that says *what went wrong* rather than merely how big it got. A
# surface should touch a small fraction of the volume's bins; if it touches
# most of them the input is a solid block and no `alpha3d_um` can extract a
# surface from it. Measured on the reference image: 19.7% works, 94.1% hangs
# forever (P35). Half is well clear of any real sheet and well below the
# observed failure.
_OCCUPANCY_REFUSE = 0.5


def check_point_cloud(n_points: int, mask_shape, scale, bin_um: float):
    """Is this cloud usable as an alpha-shape input? Pure; returns notes.

    Returns ``(ok, occupancy, notes)`` — `ok` False means refuse. Split out
    from `build_mesh` so it can be tested against a count without building a
    volume, the same doctrine as `saved_check.py`.
    """
    notes: list[str] = []
    bins = 1
    for n, s in zip(mask_shape, scale):
        bins *= max(1, int(np.ceil(n * float(s) / max(float(bin_um), 1e-9))))
    occ = n_points / bins if bins else 0.0
    if occ > _OCCUPANCY_REFUSE:
        notes.append(
            f"the subsampled cloud fills {occ:.0%} of the volume's "
            f"{bin_um:g} µm bins ({n_points:,} of {bins:,}). That is a solid "
            f"block, not a surface, so the alpha shape has nothing to extract "
            f"and Delaunay will not finish. Raise the Cadherin tab's min_size "
            f"(12 is the documented floor and is usually enough), or raise "
            f"subsample_um. See docs/ISSUES.md P35.")
        return False, occ, notes
    if n_points > _POINTS_SOFT:
        notes.append(
            f"WARN: {n_points:,} points is above voxel_subsample's design "
            f"range (~1e4-5e4); Delaunay may be very slow. Raise min_size or "
            f"subsample_um if it does not return.")
    return True, occ, notes


def _guard_points(pts, mask_shape, scale, bin_um, extras=None):
    """Refuse on `pts`, and record the check in `extras`.

    **This module never prints** — notes go into `extras` for the caller to
    format, the same contract as `measure_core` and `mesh_score`. The refusal
    is the exception, and it carries its reason in the message, so a caller
    that passes no `extras` still sees why it stopped.
    """
    ok, occ, notes = check_point_cloud(len(pts), mask_shape, scale, bin_um)
    if extras is not None:
        extras.setdefault("notes", []).extend(notes)
        extras["point_cloud_occupancy"] = occ
        extras["n_points"] = int(len(pts))
    if not ok:
        raise PointCloudTooDense(notes[0])
    return pts


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


# --- alpha solid as a volume -------------------------------------------------

def rasterize_alpha_solid(points_um, alpha_um, shape_zyx, scale,
                          grid_um) -> np.ndarray:
    """The 3D alpha shape as a filled boolean volume rather than a surface.

    `alpha_shape_faces` returns the *boundary* of the alpha complex; this
    returns its *interior*, on a `grid_um` grid spanning the same extent as the
    source image. A cell is solid iff the Delaunay tetrahedron containing it
    survived the circumradius cut — exact, via `Delaunay.find_simplex`, with no
    scan-conversion and no tolerance.

    This exists because the alpha complex is good at exactly one thing —
    bridging the wall's 8-10 um gaps in 3D, through the branch — and bad at
    being a surface (`docs/ISSUES.md` I05: 8.3% non-manifold faces, 42 pieces,
    genus 125). Handing the *solid* to a method that builds surfaces well keeps
    the first and discards the second. It is deliberately **not** exposed as a
    surface of its own: marching cubes straight off this volume was measured at
    `slope_loglog` -1.167 against alpha 3D's -0.214 (`RESULTS.md` §12), because
    a binary volume meshes to a voxel staircase and that is exactly the
    roughness a 5 um quadric fit cannot see past.

    Built one z-slice at a time; the full coordinate array for a 100 M-cell grid
    would be a multi-GB float64 copy.

    **`grid_um` interacts with the Poisson grid and is not independent of it.**
    `vessel_poisson.occupancy_grid` block-means this volume with
    ``k = round(poisson_grid_um / grid_um)`` and solves at ``h = k * grid_um``,
    so a finer raster does two things at once: it *grades* the occupancy — a thin
    spur fills only a fraction of its cell and drops below the iso level, which
    is the reason to use one — and it can also move the realised solve grid. At
    ``poisson_grid_um = 2.0``:

    ==================  ===  =========  ===========================================
    ``grid_um``          k   ``h``      effect
    ==================  ===  =========  ===========================================
    2.0                  1   2.00 um    binary solid; no spur filtering at all
    0.8                  2   1.60 um    finer solve, ~2x the vertices, much slower
    **0.65**             3   1.95 um    27-level occupancy at the solve grid already used
    ==================  ===  =========  ===========================================

    0.65 is the one to reach for: it is exactly 2 x the 0.325 um voxel and leaves
    the solve where Poisson-on-mask already runs it. See `docs/ISSUES.md` [I06].
    """
    empty = np.zeros((0, 0, 0), dtype=bool)
    if len(points_um) < 4 or grid_um <= 0:
        return empty
    h = float(grid_um)
    shape_um = [float(n) * float(sc) for n, sc in zip(shape_zyx, scale)]
    nz, ny, nx = (int(np.ceil(u / h)) + 1 for u in shape_um)
    if min(nz, ny, nx) < 3:
        return empty
    tri = Delaunay(points_um)
    keep = _tet_circumradii(points_um, tri.simplices) < float(alpha_um)
    out = np.zeros((nz, ny, nx), dtype=bool)
    if not keep.any():
        return out
    yy, xx = np.meshgrid(np.arange(ny) * h, np.arange(nx) * h, indexing="ij")
    flat = np.empty((ny * nx, 3), dtype=np.float64)
    flat[:, 1], flat[:, 2] = yy.ravel(), xx.ravel()
    for iz in range(nz):
        flat[:, 0] = iz * h
        simp = tri.find_simplex(flat)
        inside = simp >= 0
        if inside.any():
            inside[inside] = keep[simp[inside]]
        out[iz] = inside.reshape(ny, nx)
    return out


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


# --- the one build entry point -----------------------------------------------
#
# `build_mesh` below is the *only* place the strategies are dispatched.
# The GUI tab (`gui/vessel_mesh.py`) and the headless runner
# (`batch_mesh.py`) both call it, so the interactive and unattended paths
# cannot drift into building different meshes from the same parameters — the
# same doctrine `saved_check.py` was written under, one implementation and two
# callers. The GUI re-exports the constants below for its ComboBox choices.

METHOD_POLAR = "polar (per-slice)"
METHOD_ALPHA2D = "alpha 2D (per-slice)"
METHOD_ALPHA3D = "alpha 3D (points)"
METHOD_POISSON = "poisson (oriented points)"
METHOD_ALPHA_POISSON = "alpha 3D solid -> poisson"
METHODS = [METHOD_POLAR, METHOD_ALPHA2D, METHOD_ALPHA3D, METHOD_POISSON,
           METHOD_ALPHA_POISSON]
ORIENTATIONS = ["xz", "yz", "both (AND)", "both (OR)"]

# Which parameters each method actually reads. Derived from the dispatch in
# `build_mesh` and kept beside it so the two cannot disagree: a caller that
# records every knob it was handed would otherwise record values that did
# nothing (`docs/ISSUES.md` P26 — `close_um` and `r_smooth_bins` were both
# recorded as "hand-tuned" for a mesh built with alpha 3D, which reads neither).
LIVE_PARAMS = {
    METHOD_POLAR: ("orientation", "n_bins", "r_smooth_bins", "close_um",
                   "min_component_um"),
    METHOD_ALPHA2D: ("orientation", "alpha2d_um", "fill_lumen", "close_um",
                     "min_component_um"),
    METHOD_ALPHA3D: ("subsample_um", "alpha3d_um", "min_component_um"),
    METHOD_POISSON: ("poisson_grid_um", "poisson_smooth_um", "poisson_iso_pct",
                     "min_component_um"),
    METHOD_ALPHA_POISSON: ("subsample_um", "alpha3d_um", "alpha_raster_um",
                           "poisson_grid_um", "poisson_smooth_um",
                           "poisson_iso_pct", "min_component_um"),
}

ALL_PARAMS = ("method", "orientation", "n_bins", "r_smooth_bins", "alpha2d_um",
              "subsample_um", "alpha3d_um", "fill_lumen", "close_um",
              "min_component_um", "poisson_grid_um", "poisson_smooth_um",
              "poisson_iso_pct", "alpha_raster_um")


def live_params(method: str) -> tuple[str, ...]:
    """The parameters `method` reads. Unknown method → every parameter."""
    return LIVE_PARAMS.get(method, tuple(p for p in ALL_PARAMS if p != "method"))


def _per_slice(fill_fn, mask, orientation) -> np.ndarray:
    """Run a per-slice fill over the requested orientation(s) and combine."""
    if orientation in ("xz", "yz"):
        return fill_fn(mask, orientation)
    a = fill_fn(mask, "xz")
    b = fill_fn(mask, "yz")
    mode = "and" if "AND" in orientation else "or"
    return combine(a, b, mode)


def build_mesh(mask, scale, *, method, orientation, n_bins, r_smooth_bins,
               alpha2d_um, subsample_um, alpha3d_um, fill_lumen, close_um,
               min_component_um,
               poisson_grid_um=vessel_poisson.DEFAULT_GRID_UM,
               poisson_smooth_um=vessel_poisson.DEFAULT_SMOOTH_UM,
               poisson_iso_pct=vessel_poisson.DEFAULT_ISO_PCT,
               alpha_raster_um=vessel_poisson.DEFAULT_GRID_UM,
               extras=None):
    """Build the vessel surface from a cadherin mask. Pure: no Qt, no viewer.

    `mask` is any array whose non-zero voxels are wall (the filtered cadherin
    mask); `scale` is the isotropic µm/voxel triple. Returns
    ``(kind, verts, faces, values)`` where `kind` is:

      - ``"empty"``  — the mask has no foreground; the other three are None.
      - ``"mesh"``   — marching cubes on a filled volume (polar / alpha 2D).
      - ``"points"`` — the alpha shape over the subsampled point cloud, so
        `verts` are the cadherin points themselves, most of which end up
        interior to the solid (`docs/ISSUES.md` P09(a), P21).

    The Poisson method also returns ``"mesh"``: it is marching cubes on the
    fitted implicit function, so every vertex is on the surface. **That surface
    is open at the x/y domain faces by design** — the solve is Neumann because
    the vessel genuinely leaves through them, and a Dirichlet box would cap it
    there.

    ``extras``, when a dict is passed, is *filled in* with whatever the chosen
    method produces beyond the four-tuple — for both Poisson methods,
    ``"normals"`` (the per-vertex outward normals read from ∇χ), ``"info"``
    (the diagnostic record, including ``notes``) and ``"lumen"`` (the solved χ
    grid, its iso level and its origin, which `lumen_graph` builds the 3D
    centerline on). The return stays a 4-tuple because
    `saved_check`, `batch_mesh` and the GUI all destructure it. An out-dict
    rather than `vessel_poisson.diagnose`: the normals are needed for
    *correctness*, not decoration — `cell_mesh.mesh_curvature` cannot recover
    the global sign on a surface this open — and asking for them separately
    would mean building the whole reconstruction twice, which a parameter sweep
    cannot afford. `diagnose` remains for callers that want the record without
    a mesh.

    Callers must still check for 0 vertices and 0 faces: both are ordinary
    outcomes of parameters that are too tight, not errors.

    Not every argument is read — see `LIVE_PARAMS`. Alpha 3D ignores
    `orientation`, `n_bins`, `r_smooth_bins`, `alpha2d_um`, `fill_lumen` and
    `close_um` entirely; Poisson ignores everything but its own three knobs and
    `min_component_um`.
    """
    # **Not binarised here.** `mask` may be a continuous wall field
    # (`denoise.local_contrast_score`, `docs/ISSUES.md` [P28]); only Poisson can
    # read one, so the three interpolating branches take `field > 0` themselves
    # and see exactly the boolean array they saw before this line changed. That
    # is what keeps the P28 experiment controlled: adding a *field* is not
    # adding a *method*, and alpha 3D must stay bit-for-bit identical.
    field = np.asarray(mask)
    if not field.any():
        return ("empty", None, None, None)
    if method == METHOD_POLAR:
        filled = _per_slice(
            lambda m, ax: polar_fill_volume(
                m, ax, scale, int(n_bins), int(r_smooth_bins)),
            field > 0, orientation,
        )
        filled = close_volume(filled, scale[0], float(close_um))
        verts, faces, vals = mesh_from_volume(filled, scale)
        faces = filter_faces_by_extent(verts, faces, float(min_component_um))
        return ("mesh", verts, faces, vals)
    if method == METHOD_ALPHA2D:
        filled = _per_slice(
            lambda m, ax: alpha_fill_volume(
                m, ax, scale, float(alpha2d_um), fill_lumen=fill_lumen),
            field > 0, orientation,
        )
        filled = close_volume(filled, scale[0], float(close_um))
        verts, faces, vals = mesh_from_volume(filled, scale)
        faces = filter_faces_by_extent(verts, faces, float(min_component_um))
        return ("mesh", verts, faces, vals)
    if method == METHOD_POISSON:
        lumen: dict = {}
        verts, faces, vals, normals, info = vessel_poisson.reconstruct(
            field, scale, grid_um=float(poisson_grid_um),
            smooth_um=float(poisson_smooth_um),
            iso_pct=float(poisson_iso_pct),
            min_component_um=float(min_component_um),
            extras=lumen)
        if extras is not None:
            extras["normals"] = normals
            extras["info"] = info
            extras["lumen"] = lumen
        return ("mesh", verts, faces, vals)
    if method == METHOD_ALPHA_POISSON:
        # Alpha for what it is good at (bridging gaps in 3D, through the
        # branch), Poisson for what it is good at (a smooth, oriented, single
        # luminal sheet). Measured on the reference image against
        # poisson-on-mask (`RESULTS.md` §12): it holds the plateau and the
        # topology (1 piece, 0.0% non-manifold) and covers materially better --
        # <=5 um 87.9% -> 91.3%, <=10 um 93.3% -> **98.2%**.
        #
        # Its `support` **max** is 46.47 um against 26.70, which read alone
        # looks like a doubling of fabricated surface. It is not: area-weighted,
        # the two are 19.6% vs **19.9%** of area beyond 5 um and their medians
        # are 2.27 vs 2.43 um -- indistinguishable through p90, diverging only
        # in the top ~5% tail. Read `support` by area, not by max.
        #
        # And `support` is scored against the *retained* mask, which [P19]/[P28]
        # measure as missing 59-85% of the dim wall -- so bridged surface is
        # marked unsupported wherever the mask is wrong, which is precisely
        # where this method works. The yardstick is biased against it.
        #
        # Still not adjudicated here: `mesh_score` returns a vector and never a
        # total, and the selection rule belongs to whoever runs the sweep
        # (`docs/ISSUES.md` I03). Alpha 3D remains the default.
        pts = _guard_points(
            voxel_subsample(mask_to_points_um(field > 0, scale),
                            float(subsample_um)),
            field.shape, scale, float(subsample_um), extras)
        solid = rasterize_alpha_solid(pts, float(alpha3d_um), field.shape,
                                      scale, float(alpha_raster_um))
        if not solid.any():
            return ("empty", None, None, None)
        lumen: dict = {}
        verts, faces, vals, normals, info = vessel_poisson.reconstruct(
            solid, (float(alpha_raster_um),) * 3,
            grid_um=float(poisson_grid_um),
            smooth_um=float(poisson_smooth_um),
            iso_pct=float(poisson_iso_pct),
            min_component_um=float(min_component_um),
            extras=lumen)
        if extras is not None:
            extras["normals"] = normals
            extras["info"] = info
            extras["lumen"] = lumen
        return ("mesh", verts, faces, vals)
    # 3D alpha on the subsampled cadherin point cloud. **This is also the
    # fall-through for an unknown method string** — a latent defect recorded as
    # its own entry rather than papered over here, since fixing it would change
    # what an existing caller silently gets.
    # `mask_to_points_um` is `np.argwhere`, which on a float field would take
    # every nonzero voxel — the `> 0` is what stops that.
    pts = _guard_points(
        voxel_subsample(mask_to_points_um(field > 0, scale),
                        float(subsample_um)),
        field.shape, scale, float(subsample_um), extras)
    faces = alpha_shape_faces(pts, float(alpha3d_um))
    faces = filter_faces_by_extent(pts, faces, float(min_component_um))
    vals = (np.linalg.norm(pts - pts.mean(0), axis=1)
            if len(pts) else np.zeros((0,)))
    return ("points", pts, faces, vals)
