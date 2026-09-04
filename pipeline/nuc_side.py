"""Which side of the vessel wall is a nucleus on?

Both cell types in this dataset are DAPI-positive, and separating them is the
core unsolved problem ([P05](../docs/ISSUES.md)). The current classifier reads
mClover in a ~2 µm perinuclear shell — precisely where both cell types look most
alike — and swings ~2x across one slider step.

This module measures something else entirely: **geometry**. Fibroblasts were
seeded onto an already-formed, already-perfused vessel, so they are
**abluminal** — outside the wall. Endothelial cells *are* the wall. So the side
of the surface a nucleus sits on is a physical statement about which cell type
it belongs to, sharing no machinery and no failure mode with a marker intensity.

**An earlier attempt at this failed, and the reason is now known.** It looked
for two lumps in the nucleus-to-surface distance and found one. The distance it
used was *unsigned*, which folds the inside and outside populations onto each
other by construction — a nucleus 6 µm inside and one 6 µm outside both read
6.0. Signed, on image 1, the same data separates: AUC 0.680 against the marker
classifier at p = 0.0009, with a monotone dose-response (DIRECTION §8A.1).

Three things this does that the failed version did not, all from §7.2:

1. **Signs the distance**, from the surface normal. On the alpha->poisson mesh
   the stored per-vertex normals come from ∇χ, so positive is away from the
   lumen — genuinely abluminal. **On any other mesh the sign is not
   luminal/abluminal** and this module refuses rather than reporting a number
   that looks fine: [P21] established the alpha shape is a watertight solid
   around the wall *slab*, so "outward" points into the lumen on the inner
   sheet.
2. **Measures to the nearest triangle, not the nearest vertex.** At ~1-2 µm
   vertex spacing the nearest-vertex estimator is biased *outward*, not merely
   noisy ([P09](b)) — and the whole expected effect is 5-10 µm.
3. **Calibrates against nuclei that are definitionally endothelial** — the ones
   in contact with the VE-cadherin junction network — instead of hunting for
   bimodality. That turns a hopeless two-population deconvolution into a
   one-sided comparison against a measured reference, which works even when the
   two distributions overlap heavily. And they do overlap: §8A.1 measured the
   signed distribution as unimodal even when correctly signed.

**This is not validated either.** It agrees with the marker classifier on image
1, and two unvalidated measures agreeing is not validation — but they share no
machinery, which is the first external corroboration [P05] has ever had.
Nothing here writes a cell-type label; it reports a signed distance and the
calibration, and leaves the call to the analysis.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from . import saved_io

# A nucleus whose voxels come within this of the VE-cadherin **junction
# network** is taken as endothelial for calibration. VE-cadherin marks
# endothelial cell-cell junctions, so contact with the network is a positive
# statement about the cell, not an inference from a missing marker.
DEFAULT_CONTACT_UM = 1.0
# ...but "the mask" is not "the network". Measured on the reference image: the
# size-filtered cadherin mask is 3,764 components, and **232 of 276 nuclei
# overlap one outright** (contact distance exactly 0.0 µm), so contact with the
# mask calls 93% of nuclei endothelial and discriminates nothing. Keeping only
# components above this volume leaves **52** of them — holding **59%** of the
# mask's voxels, i.e. the rings rather than the speckle — and the contacting set
# drops to 41%. Stated in µm³ so it does not silently change meaning with voxel
# size.
DEFAULT_NETWORK_MIN_UM3 = 150.0
# How many nearest triangles to test exactly before the correctness pass. Only
# a speed knob: the second pass is what makes the answer exact.
_KNN = 32


# --------------------------------------------------------------------------- #
# exact point-to-triangle distance
# --------------------------------------------------------------------------- #

def _point_triangle_closest(p: np.ndarray, a: np.ndarray, b: np.ndarray,
                            c: np.ndarray) -> np.ndarray:
    """Closest point on each triangle ``(a,b,c)`` to the matching point ``p``.

    Ericson's region-based solution, vectorised: the six Voronoi regions of a
    triangle (three vertices, three edges, the interior) each have a closed
    form, and the barycentric tests select between them. Exact — no iteration.
    """
    ab, ac, ap = b - a, c - a, p - a
    d1 = np.einsum("ij,ij->i", ab, ap)
    d2 = np.einsum("ij,ij->i", ac, ap)
    bp = p - b
    d3 = np.einsum("ij,ij->i", ab, bp)
    d4 = np.einsum("ij,ij->i", ac, bp)
    cp = p - c
    d5 = np.einsum("ij,ij->i", ab, cp)
    d6 = np.einsum("ij,ij->i", ac, cp)

    vc = d1 * d4 - d3 * d2
    vb = d5 * d2 - d1 * d6
    va = d3 * d6 - d5 * d4
    denom = va + vb + vc
    with np.errstate(divide="ignore", invalid="ignore"):
        v_i = np.where(np.abs(denom) > 1e-20, vb / denom, 0.0)
        w_i = np.where(np.abs(denom) > 1e-20, vc / denom, 0.0)
        v_ab = np.where(np.abs(d1 - d3) > 1e-20, d1 / (d1 - d3), 0.0)
        w_ac = np.where(np.abs(d2 - d6) > 1e-20, d2 / (d2 - d6), 0.0)
        den_bc = (d4 - d3) + (d5 - d6)
        w_bc = np.where(np.abs(den_bc) > 1e-20, (d4 - d3) / den_bc, 0.0)

    out = a + v_i[:, None] * ab + w_i[:, None] * ac          # interior
    m_a = (d1 <= 0) & (d2 <= 0)
    m_b = (d3 >= 0) & (d4 <= d3)
    m_c = (d6 >= 0) & (d5 <= d6)
    m_ab = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
    m_ac = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
    m_bc = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)

    out = np.where(m_bc[:, None], b + w_bc[:, None] * (c - b), out)
    out = np.where(m_ac[:, None], a + w_ac[:, None] * ac, out)
    out = np.where(m_ab[:, None], a + v_ab[:, None] * ab, out)
    out = np.where(m_c[:, None], c, out)
    out = np.where(m_b[:, None], b, out)
    out = np.where(m_a[:, None], a, out)
    return out


def signed_distance_to_triangles(points: np.ndarray, verts: np.ndarray,
                                 faces: np.ndarray, normals: np.ndarray):
    """Signed distance from each point to the **surface**, not to its vertices.

    Positive is on the side the normals point to. Returns
    ``(signed, unsigned, face_idx)``.

    Two passes, because one is not exact. The first takes the ``_KNN`` nearest
    triangle *centroids* and measures exactly to those. The second is the one
    that makes it correct: a triangle can only beat the best distance found so
    far if its centroid lies within ``d_best + R_max`` of the point, where
    ``R_max`` is the largest triangle circumradius, so re-testing everything in
    that ball cannot miss the true nearest. Skipping it would leave a
    silently-approximate distance, which is the failure mode this whole module
    exists to correct.

    The sign is read from the **interpolated vertex normal** at the closest
    point rather than from the face normal: the ∇χ normals are a per-vertex
    field, and interpolating keeps the sign stable for a point sitting almost
    exactly on the surface, where adjacent face normals disagree most.
    """
    points = np.asarray(points, float)
    verts = np.asarray(verts, float)
    faces = np.asarray(faces, int)
    normals = np.asarray(normals, float)
    if len(points) == 0 or len(faces) == 0:
        n = len(points)
        return (np.full(n, np.nan), np.full(n, np.nan),
                np.full(n, -1, int))

    tri = verts[faces]                                  # (F,3,3)
    cen = tri.mean(axis=1)
    rmax = float(np.linalg.norm(tri - cen[:, None, :], axis=2).max())
    tree = cKDTree(cen)

    k = min(_KNN, len(cen))
    _d, idx = tree.query(points, k=k)
    idx = np.atleast_2d(idx.T).T if k > 1 else idx.reshape(-1, 1)
    P = np.repeat(points, idx.shape[1], axis=0)
    F = idx.ravel()
    cp = _point_triangle_closest(P, tri[F, 0], tri[F, 1], tri[F, 2])
    d = np.linalg.norm(P - cp, axis=1).reshape(len(points), -1)
    best = np.argmin(d, axis=1)
    d_best = d[np.arange(len(points)), best]
    f_best = idx[np.arange(len(points)), best]
    cp_best = cp.reshape(len(points), -1, 3)[np.arange(len(points)), best]

    # Exactness pass.
    for i, (p, db) in enumerate(zip(points, d_best)):
        cand = tree.query_ball_point(p, db + rmax)
        if len(cand) <= idx.shape[1]:
            continue
        cand = np.asarray(cand, int)
        q = np.repeat(p[None, :], len(cand), axis=0)
        c2 = _point_triangle_closest(q, tri[cand, 0], tri[cand, 1],
                                     tri[cand, 2])
        dd = np.linalg.norm(q - c2, axis=1)
        j = int(np.argmin(dd))
        if dd[j] < d_best[i]:
            d_best[i] = dd[j]
            f_best[i] = cand[j]
            cp_best[i] = c2[j]

    # Sign from the interpolated vertex normal at the closest point.
    tf = tri[f_best]
    w = np.stack([1.0 / (np.linalg.norm(cp_best - tf[:, j], axis=1) + 1e-9)
                  for j in range(3)], axis=1)
    w /= w.sum(axis=1, keepdims=True)
    nrm = np.einsum("ij,ijk->ik", w, normals[faces[f_best]])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12
    sgn = np.sign(np.einsum("ij,ij->i", points - cp_best, nrm))
    sgn[sgn == 0] = 1.0
    return d_best * sgn, d_best, f_best


# --------------------------------------------------------------------------- #
# the measurement
# --------------------------------------------------------------------------- #

def nucleus_centroids(labels: np.ndarray, ids, scale_zyx):
    """Centroid of every nucleus, in µm."""
    ids = np.asarray(ids)
    if len(ids) == 0:
        return np.zeros((0, 3)), np.zeros(0, int)
    cen = np.asarray(ndi.center_of_mass(labels > 0, labels, ids))
    return cen * np.asarray(scale_zyx, float), ids


def junction_network(cadherin_mask: np.ndarray, scale_zyx,
                     min_component_um3: float = DEFAULT_NETWORK_MIN_UM3):
    """The VE-cadherin mask reduced to components big enough to be junctions.

    The mask at the shipped despeckle cut is overwhelmingly small components
    ([P19]: 95%+ of components are 1-4 voxels), and a nucleus sitting anywhere
    near the wall touches one. Dropping to the large components is what makes
    "touches the junction network" a statement about the cell rather than about
    the speckle density.
    """
    vox = float(np.prod(np.asarray(scale_zyx, float)))
    lab, n = ndi.label(np.asarray(cadherin_mask) > 0)
    if n == 0:
        return np.zeros(np.shape(cadherin_mask), bool), {"n_components": 0}
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    keep = sizes * vox >= float(min_component_um3)
    info = {"n_components": int(n), "n_kept": int(keep.sum()),
            "voxels_kept": int(sizes[keep].sum()),
            "voxels_total": int(sizes.sum()),
            "min_component_um3": float(min_component_um3)}
    return keep[lab], info


def cadherin_contact_um(labels: np.ndarray, ids, cadherin_mask: np.ndarray,
                        scale_zyx, max_voxels_per_nucleus: int = 4000):
    """Distance from each nucleus to the nearest VE-cadherin voxel, in µm.

    A KD-tree over the cadherin voxels rather than a distance transform of the
    whole volume: the stacks are ~340 M voxels and the mask is ~1.5 M, so the
    tree is two orders of magnitude smaller than the field it replaces.

    Nuclei larger than ``max_voxels_per_nucleus`` are sampled down — this is a
    *minimum* over the nucleus, which a subsample can only over-estimate, and
    the bias is bounded by the voxel spacing. That is a different situation from
    `measure_core`'s ``::3``, which subsamples the thing being measured rather
    than a set being minimised over.
    """
    ids = np.asarray(ids)
    out = np.full(len(ids), np.nan)
    cad = np.argwhere(np.asarray(cadherin_mask) > 0)
    if len(cad) == 0 or len(ids) == 0:
        return out
    sc = np.asarray(scale_zyx, float)
    tree = cKDTree(cad * sc)
    objs = ndi.find_objects(labels)
    rng = np.random.default_rng(0)
    for i, lid in enumerate(ids):
        sl = objs[int(lid) - 1] if int(lid) - 1 < len(objs) else None
        if sl is None:
            continue
        sub = labels[sl] == lid
        vox = np.argwhere(sub) + np.array([s.start for s in sl])
        if len(vox) == 0:
            continue
        if len(vox) > max_voxels_per_nucleus:
            vox = vox[rng.choice(len(vox), max_voxels_per_nucleus,
                                 replace=False)]
        d, _ = tree.query(vox * sc, k=1)
        out[i] = float(d.min())
    return out


def measure(labels: np.ndarray, ids, is_fib, verts, faces, normals,
            scale_zyx, cadherin_mask=None,
            contact_um: float = DEFAULT_CONTACT_UM,
            network_min_um3: float = DEFAULT_NETWORK_MIN_UM3) -> dict:
    """Signed nucleus-to-vessel distance for one image, plus the calibration.

    ``is_fib`` is the [P05] marker classifier's call, carried through only so
    the two independent measurements can be compared. **It is not used to
    compute anything here.**

    Returns per-nucleus arrays and a summary. The ``endothelial_reference``
    block is the calibration §7.2 asks for: the location and spread of the
    signed distance among nuclei touching the VE-cadherin network, which are
    endothelial by definition of the marker.
    """
    out: dict = {"notes": []}
    cents, ids = nucleus_centroids(labels, ids, scale_zyx)
    out["n_nuclei"] = int(len(cents))
    if len(cents) == 0 or len(faces) == 0:
        out["notes"].append("WARN: no nuclei or no surface — nothing measured.")
        return out

    signed, unsigned, _f = signed_distance_to_triangles(
        cents, verts, faces, normals)
    out["centroid_zyx_um"] = cents.astype(np.float32)
    out["nuc_id"] = np.asarray(ids, np.int32)
    out["signed_um"] = signed.astype(np.float32)
    out["unsigned_um"] = unsigned.astype(np.float32)
    out["is_fib_marker"] = (np.asarray(is_fib).astype(bool)
                            if is_fib is not None
                            else np.zeros(len(cents), bool))

    if cadherin_mask is not None:
        net, ninfo = junction_network(cadherin_mask, scale_zyx,
                                      network_min_um3)
        out["junction_network"] = ninfo
        cd = cadherin_contact_um(labels, ids, net, scale_zyx)
        out["cadherin_contact_um"] = cd.astype(np.float32)
        endo = np.isfinite(cd) & (cd <= contact_um)
        out["is_endothelial_by_contact"] = endo
        if endo.sum() >= 5:
            ref = signed[endo]
            out["endothelial_reference"] = {
                "n": int(endo.sum()), "contact_um": float(contact_um),
                "median": float(np.median(ref)),
                "p10": float(np.percentile(ref, 10)),
                "p90": float(np.percentile(ref, 90)),
                "mad": float(np.median(np.abs(ref - np.median(ref)))),
            }
            rest = signed[~endo]
            if len(rest) >= 5:
                from scipy.stats import mannwhitneyu
                u, p = mannwhitneyu(rest, ref, alternative="greater")
                out["displacement_test"] = {
                    "n_other": int(len(rest)), "n_reference": int(endo.sum()),
                    "auc": float(u / (len(rest) * len(ref))),
                    "p_one_sided": float(p),
                    "median_other": float(np.median(rest)),
                    "median_reference": float(np.median(ref)),
                    "delta_um": float(np.median(rest) - np.median(ref)),
                }
        else:
            out["notes"].append(
                f"WARN: only {int(endo.sum())} nuclei contact the cadherin "
                f"network within {contact_um:g} µm — too few to calibrate an "
                f"endothelial reference distribution.")
    return out


def mesh_sign_is_luminal(mesh_path) -> tuple[bool, str]:
    """Is this mesh's normal field a luminal/abluminal statement?

    Only the Poisson family carries ∇χ normals, which point away from the lumen
    volume. On an alpha-3D mesh the surface is a watertight solid around the
    wall slab ([P21]), so "outward" points **into** the lumen on the inner
    sheet and the sign of a distance to it means nothing. Refusing is the point:
    the same array of numbers comes back either way, and only the provenance
    distinguishes a measurement from noise.
    """
    # `saved_io.has_normals` reads the archive directory only. It must stay
    # that and not become `load_mesh`: inflating a multi-million-vertex mesh to
    # answer a yes/no is the cost this function exists to avoid.
    if not saved_io.has_normals(mesh_path):
        return False, (f"{mesh_path.name} carries no stored normals, so it "
                       f"is not a Poisson reconstruction. On an alpha-3D "
                       f"solid the outward direction is not abluminal "
                       f"([P21]) and the sign would be meaningless.")
    return True, ""
