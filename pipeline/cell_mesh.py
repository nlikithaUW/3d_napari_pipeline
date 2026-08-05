"""Lumen surface mesh from DAPI nuclei as vertices.

Each cell contributes one vertex (its nucleus centroid in µm), so the mesh
encodes the cell adjacency graph embedded in 3D — useful downstream for
per-cell force vs local geometry.

Pipeline:
  1. Per-nucleus centroid + PCA on the DAPI label map. Smallest eigenvector
     ≈ surface normal (nuclei are oblong, lying tangent to the vessel wall:
     ~20 µm long axes vs ~3 µm short axis → well-conditioned).
  2. 3D alpha shape on the centroids: Delaunay tetrahedralization, drop tets
     whose circumradius exceeds `alpha_um`, keep faces that bound exactly
     one surviving tet. Small alpha → follow lumen curvature and branches;
     large alpha → approach the convex hull.
  3. Optional triangle pruning: drop triangles whose face normal disagrees
     with the average vertex (nuclear) normal — kills spurious caps across
     vessel openings.
  4. Optional projection of centroids outward along their nuclear normal
     onto the cadherin signal peak, to recover the true vessel wall (nuclei
     sit a few µm inside).
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import Delaunay


def label_centroids_normals(
    labels: np.ndarray,
    voxel_zyx_um,
    min_voxels: int = 10,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-label centroid (µm) and short-axis normal (unit, in µm space).

    Returns (centroids_zyx, normals_zyx, label_ids). Labels with fewer than
    `min_voxels` voxels are skipped. PCA on tiny clouds is numerically
    unstable; callers using a small `min_voxels` should treat the returned
    normals as best-effort.
    """
    scale = np.array(voxel_zyx_um, dtype=np.float64)
    slices = ndi.find_objects(labels)

    centroids = []
    normals = []
    kept = []
    for idx, sl in enumerate(slices):
        if sl is None:
            continue
        lid = idx + 1
        sub = labels[sl]
        offset = np.array([s.start for s in sl], dtype=np.float64)
        coords = np.argwhere(sub == lid).astype(np.float64)
        if len(coords) < min_voxels:
            continue
        coords = (coords + offset) * scale  # µm
        c = coords.mean(0)
        X = coords - c
        if len(coords) >= 3:
            cov = (X.T @ X) / len(X)
            _, V = np.linalg.eigh(cov)
            n = V[:, 0]
            n = n / (np.linalg.norm(n) + 1e-12)
        else:
            n = np.array([1.0, 0.0, 0.0])
        centroids.append(c)
        normals.append(n)
        kept.append(lid)

    if not centroids:
        return (np.zeros((0, 3)), np.zeros((0, 3)), np.zeros((0,), dtype=int))
    return (np.array(centroids), np.array(normals), np.array(kept, dtype=int))


def nuclei_vertices(labels: np.ndarray, voxel_zyx_um) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-nucleus centroid (µm) and short-axis normal (unit, in µm space).

    Returns (centroids_zyx, normals_zyx, label_ids). Labels with fewer than
    10 voxels are skipped.
    """
    return label_centroids_normals(labels, voxel_zyx_um, min_voxels=10)


def _tet_circumradii(pts: np.ndarray, tets: np.ndarray) -> np.ndarray:
    """Circumradius of every tetrahedron (µm). Degenerate tets → +inf."""
    p = pts[tets]                         # (n, 4, 3)
    a = p[:, 1] - p[:, 0]
    b = p[:, 2] - p[:, 0]
    c = p[:, 3] - p[:, 0]
    A = np.stack([a, b, c], axis=1)       # (n, 3, 3)
    rhs = 0.5 * np.stack([
        np.einsum("ij,ij->i", a, a),
        np.einsum("ij,ij->i", b, b),
        np.einsum("ij,ij->i", c, c),
    ], axis=1)                            # (n, 3)
    det = np.linalg.det(A)
    ok = np.abs(det) > 1e-9
    R = np.full(len(tets), np.inf)
    if ok.any():
        center_rel = np.linalg.solve(A[ok], rhs[ok][..., None])[..., 0]
        R[ok] = np.linalg.norm(center_rel, axis=1)
    return R


def tet_circumradii(points_um: np.ndarray) -> np.ndarray:
    """Circumradii (µm) of every Delaunay tetrahedron on `points_um`."""
    if len(points_um) < 4:
        return np.array([], dtype=np.float64)
    tri = Delaunay(points_um)
    return _tet_circumradii(points_um, tri.simplices)


def suggest_alpha_um(circumradii: np.ndarray) -> float:
    """Pick the first valley on log10(circumradii) — separates surface tets
    from cross-lumen / hull tets.

    The circumradius distribution is typically trimodal: small surface tets,
    medium cross-lumen tets, and near-degenerate hull tets with huge R.
    Single Otsu lands in the larger gap and overshoots the surface scale,
    so we use multi-Otsu (3 classes) and return the LOWER of the two
    thresholds. Falls back to single Otsu if there aren't enough samples.
    """
    R = circumradii[np.isfinite(circumradii) & (circumradii > 0)]
    if R.size < 8:
        return 0.0
    from skimage.filters import threshold_multiotsu, threshold_otsu
    log_R = np.log10(R)
    try:
        thrs = threshold_multiotsu(log_R, classes=3)
        log_thr = float(thrs[0])
    except ValueError:
        log_thr = float(threshold_otsu(log_R))
    return float(10 ** log_thr)


def alpha_shape_faces(points_um: np.ndarray, alpha_um: float) -> np.ndarray:
    """Boundary triangles of the 3D alpha shape. Returns int array (M, 3)
    indexing into `points_um`. Empty (0, 3) if no tets survive.
    """
    if len(points_um) < 4:
        return np.zeros((0, 3), dtype=np.int64)
    tri = Delaunay(points_um)
    tets = tri.simplices
    R = _tet_circumradii(points_um, tets)
    tets = tets[R < alpha_um]
    if len(tets) == 0:
        return np.zeros((0, 3), dtype=np.int64)

    faces = np.vstack([
        tets[:, [0, 1, 2]],
        tets[:, [0, 1, 3]],
        tets[:, [0, 2, 3]],
        tets[:, [1, 2, 3]],
    ])
    faces = np.sort(faces, axis=1).astype(np.int64, copy=False)
    # Hash each (i, j, k) row as one void value — ~5x faster than np.unique(axis=0).
    fc = np.ascontiguousarray(faces)
    view = fc.view([("", fc.dtype)] * 3).reshape(-1)
    uniq_v, counts = np.unique(view, return_counts=True)
    boundary = uniq_v[counts == 1]
    return boundary.view(fc.dtype).reshape(-1, 3)


def filter_faces_by_extent(
    verts: np.ndarray,
    faces: np.ndarray,
    min_extent_um: float,
) -> np.ndarray:
    """Keep only faces in connected mesh components of sufficient spatial size.

    Components are the connected components of the face-adjacency graph (two
    faces adjacent iff they share an edge). A component's size is the diagonal
    of its vertex bounding box (µm). Faces in components whose diagonal is below
    `min_extent_um` are dropped — this removes the small disconnected islands
    that both the 3D-alpha point mesh and the 2D-alpha marching-cubes mesh leave
    floating around the main vessel wall.

    `min_extent_um <= 0` or empty `faces` → returns `faces` unchanged.
    """
    if min_extent_um <= 0 or len(faces) == 0:
        return faces
    faces = np.asarray(faces)
    nf = len(faces)

    # The 3 undirected edges of every face, endpoints sorted so shared edges
    # compare equal regardless of winding (same edge-hashing style as
    # `alpha_shape_faces`).
    edges = np.concatenate([
        faces[:, [0, 1]],
        faces[:, [1, 2]],
        faces[:, [0, 2]],
    ], axis=0)
    edges = np.sort(edges, axis=1)
    face_of_edge = np.tile(np.arange(nf), 3)

    # Sort edges together so identical edges land adjacent; consecutive-equal
    # rows then link the two faces that share that edge.
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    e_sorted = edges[order]
    f_sorted = face_of_edge[order]
    same = np.all(e_sorted[:-1] == e_sorted[1:], axis=1)
    rows = f_sorted[:-1][same]
    cols = f_sorted[1:][same]

    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    graph = coo_matrix(
        (np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(nf, nf),
    )
    n_comp, face_comp = connected_components(graph, directed=False)

    # Per-component vertex bbox via scatter-min/max over each face's 3 vertices.
    comp_min = np.full((n_comp, 3), np.inf)
    comp_max = np.full((n_comp, 3), -np.inf)
    comp_per_vert = np.repeat(face_comp, 3)
    pts = verts[faces.reshape(-1)]
    np.minimum.at(comp_min, comp_per_vert, pts)
    np.maximum.at(comp_max, comp_per_vert, pts)
    diag = np.linalg.norm(comp_max - comp_min, axis=1)

    keep_comp = diag >= float(min_extent_um)
    return faces[keep_comp[face_comp]]


def prune_by_normal(
    points_um: np.ndarray,
    normals: np.ndarray,
    faces: np.ndarray,
    min_abs_dot: float,
) -> np.ndarray:
    """Drop triangles whose face normal disagrees with the mean vertex
    normal. Vertex normals are sign-ambiguous (PCA), so we compare |dot|.
    """
    if len(faces) == 0 or min_abs_dot <= 0.0:
        return faces
    p = points_um[faces]                                  # (M, 3, 3)
    fn = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12
    vn = normals[faces].mean(axis=1)                      # (M, 3)
    vn /= np.linalg.norm(vn, axis=1, keepdims=True) + 1e-12
    keep = np.abs(np.einsum("ij,ij->i", fn, vn)) >= min_abs_dot
    return faces[keep]


def vertex_outward_normals(
    points_um: np.ndarray,
    faces: np.ndarray,
) -> np.ndarray:
    """Per-vertex outward unit normal from the mesh itself.

    Averages incident face normals and sign-flips each vertex normal so it
    points away from the mesh centroid (works while the lumen is roughly
    star-shaped about its own axis — fine for short vessel segments).
    """
    if len(points_um) == 0 or len(faces) == 0:
        return np.zeros_like(points_um)
    p = points_um[faces]
    fn = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    fn /= np.linalg.norm(fn, axis=1, keepdims=True) + 1e-12

    vn = np.zeros_like(points_um)
    for k in range(3):
        np.add.at(vn, faces[:, k], fn)
    vn /= np.linalg.norm(vn, axis=1, keepdims=True) + 1e-12

    centroid = points_um.mean(0)
    radial = points_um - centroid
    radial /= np.linalg.norm(radial, axis=1, keepdims=True) + 1e-12
    flip = np.einsum("ij,ij->i", vn, radial) < 0
    vn[flip] = -vn[flip]
    return vn


def _taubin_smooth(verts: np.ndarray, faces: np.ndarray, iterations: int,
                   lam: float = 0.5, mu: float = -0.53) -> np.ndarray:
    """Taubin (λ|μ) mesh smoothing — low-pass without the shrinkage that plain
    Laplacian smoothing causes (which would bias the estimated caliber). Returns
    new vertex positions; connectivity (faces) is unchanged.
    """
    from scipy.sparse import coo_matrix

    verts = np.asarray(verts, dtype=np.float64)
    faces = np.asarray(faces)
    n = len(verts)
    if iterations <= 0 or n == 0 or len(faces) == 0:
        return verts
    e = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    e = np.vstack([e, e[:, ::-1]])
    adj = (coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])),
                      shape=(n, n)).tocsr() > 0).astype(np.float64)
    deg = np.asarray(adj.sum(1)).ravel()
    deg[deg == 0] = 1.0
    v = verts.copy()
    for _ in range(int(iterations)):
        for fac in (lam, mu):
            v = v + fac * (adj.dot(v) / deg[:, None] - v)
    return v


def mesh_curvature(verts: np.ndarray, faces: np.ndarray, radius_um: float,
                   smooth_iter: int = 0) -> dict:
    """Per-vertex curvature by fitting a local quadric at scale ``radius_um``.

    For each vertex, the neighbours within ``radius_um`` are expressed in a local
    frame whose z-axis is the outward vertex normal (so the sign is consistent:
    convex-outward is positive), then a quadric
    ``z = a·x² + b·xy + c·y² + d·x + e·y`` is least-squares fit. The principal
    curvatures are the eigenvalues of the shape operator (second fundamental form
    relative to the first). Fitting over a *physical* radius smooths the
    triangle-scale segmentation noise, so ``radius_um`` sets the length scale at
    which "curvature" is measured.

    IMPORTANT: the alpha-shape "point cloud" vessel meshes carry many interior
    vertices that are not part of the boundary surface (they have no incident
    face). Fitting a quadric to those interior points is meaningless and badly
    contaminates the estimate. So the fit is restricted to *surface* vertices
    (those referenced by ``faces``); interior/rank-deficient vertices are left
    NaN so callers can exclude them. ``smooth_iter`` optionally Taubin-smooths
    the surface first.

    Returns a dict of length-``len(verts)`` float32 arrays: ``mean`` = H =
    (κ1+κ2)/2, ``k1`` = κmax, ``k2`` = κmin, ``gaussian`` = K = κ1·κ2.
    """
    from scipy.spatial import cKDTree

    verts = np.asarray(verts, dtype=np.float64)
    faces = np.asarray(faces)
    n = len(verts)
    out = {k: np.full(n, np.nan, dtype=np.float32)
           for k in ("mean", "k1", "k2", "gaussian")}
    if n < 6 or len(faces) == 0:
        return out

    # Surface vertices only (referenced by at least one face) — reindex a
    # compact sub-mesh so the KD-tree / normals never see interior points.
    surf = np.unique(faces)
    if len(surf) < 6:
        return out
    remap = np.full(n, -1, dtype=np.int64)
    remap[surf] = np.arange(len(surf))
    sv = verts[surf]
    sf = remap[faces]
    if smooth_iter and smooth_iter > 0:
        sv = _taubin_smooth(sv, sf, int(smooth_iter))

    normals = vertex_outward_normals(sv, sf)
    tree = cKDTree(sv)
    r = float(radius_um)
    m = len(sv)
    K1 = np.full(m, np.nan, dtype=np.float32)
    K2 = np.full(m, np.nan, dtype=np.float32)
    H = np.full(m, np.nan, dtype=np.float32)
    K = np.full(m, np.nan, dtype=np.float32)

    for i in range(m):
        nz = normals[i]
        if not np.any(nz):
            continue
        idx = tree.query_ball_point(sv[i], r)
        if len(idx) < 6:
            _, idx = tree.query(sv[i], k=min(12, m))
            idx = np.atleast_1d(idx)
        P = sv[idx] - sv[i]
        # Local tangent frame: z = outward normal, t1/t2 span the tangent plane.
        seed = np.array([0.0, 1.0, 0.0]) if abs(nz[0]) > 0.9 else np.array([1.0, 0.0, 0.0])
        t1 = np.cross(nz, seed)
        t1 /= np.linalg.norm(t1) + 1e-12
        t2 = np.cross(nz, t1)
        x = P @ t1
        y = P @ t2
        z = P @ nz
        A = np.column_stack([x * x, x * y, y * y, x, y])
        try:
            coef, *_ = np.linalg.lstsq(A, z, rcond=None)
        except Exception:
            continue
        a, b, c, d, e = coef
        E = 1.0 + d * d
        F = d * e
        G = 1.0 + e * e
        den = np.sqrt(d * d + e * e + 1.0)
        L = 2.0 * a / den
        M = b / den
        N = 2.0 * c / den
        try:
            S = np.linalg.solve(np.array([[E, F], [F, G]]),
                                np.array([[L, M], [M, N]]))
            ev = np.real(np.linalg.eigvals(S))
        except Exception:
            continue
        kmax = float(ev.max())
        kmin = float(ev.min())
        K1[i] = kmax
        K2[i] = kmin
        H[i] = 0.5 * (kmax + kmin)
        K[i] = kmax * kmin

    out["k1"][surf] = K1
    out["k2"][surf] = K2
    out["mean"][surf] = H
    out["gaussian"][surf] = K
    return out


def signed_distance_to_mesh_vertices(
    query_zyx_um: np.ndarray,
    vert_zyx_um: np.ndarray,
    vert_normals_zyx: np.ndarray,
    tree=None,
) -> np.ndarray:
    """Signed distance from each query point to the nearest mesh vertex.

    Sign comes from `dot(query - vertex, vertex_outward_normal)`: positive
    if the query is on the outward side of the surface, negative inside.
    Uses nearest-vertex (not nearest-face) — an approximation that's fine
    when vertex spacing is small compared to the offsets we care about.

    `tree` may be a prebuilt `cKDTree` over `vert_zyx_um` to skip rebuilding
    on every call (e.g. when sliding a slider).
    """
    if len(vert_zyx_um) == 0 or len(query_zyx_um) == 0:
        return np.zeros(len(query_zyx_um), dtype=np.float32)
    if tree is None:
        from scipy.spatial import cKDTree
        tree = cKDTree(vert_zyx_um)
    d, idx = tree.query(query_zyx_um, k=1)
    diff = query_zyx_um - vert_zyx_um[idx]
    sign = np.sign(np.einsum("ij,ij->i", diff, vert_normals_zyx[idx]))
    sign[sign == 0] = 1.0
    return (sign * d).astype(np.float32)
