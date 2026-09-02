"""The vessel centerline as a **3D graph on the Poisson lumen**.

This replaces `centerline.py`'s in-plane centerline (`docs/ISSUES.md` §8B of
`docs/DIRECTION.md`). The old module skeletonized the **max projection** of the
cadherin mask, which had to be made solid first with `binary_closing` +
`binary_fill_holes`; that one decision is the root of seven recorded defects —
[P08](a) Euclidean junction contamination, [P08](c) border flare, [P23]
amputated tapers, [P27] raster ordering, [P29] the fabricated branch node, the
auto border pad that can exceed the field of view, and the unmeasurable in-plane
width. It was the right call when it was written, because there was no
trustworthy 3D surface to work from. There is now.

**The input is χ, not the mask.** `vessel_poisson.reconstruct` already solves a
smooth indicator function whose level set is the vessel surface, and hands it
over through ``extras["lumen"]``. Above the iso level is the lumen: smooth by
construction, its wall gaps spanned by the harmonic solve, a single component,
and produced with no fill, no closing and no projection. Because the surface is
a level set of the *same* field, the centerline and the surface cannot disagree
about where the vessel is — which two separately-derived geometries always
eventually do.

What this buys, in the order the program needs it:

1. **Caliber everywhere** — the 3D distance transform inside the lumen, not
   sampled at skeleton points and not gated on a footprint, which is what made
   the old inscribed radius read a clearance from *outside* the vessel between
   the arms of a Y (§7.6).
2. **A graph, not a point set.** Nodes are junctions and endpoints, edges are
   ordered polylines. **Arc length exists by construction**, which is what
   [P27] blocks and what every "per unit length along the vessel" quantity
   needs.
3. **Scale-aware pruning.** A branch shorter than the local caliber is a spur.
   [P23]'s absolute 4 µm / 30 µm constants cannot adapt to a taper; a rule
   stated as a multiple of the caliber can.
4. **Endpoints are classified.** One that reaches an x/y domain face is a
   **crop** and everything near it is censored, not measured. No image in this
   dataset contains a real vessel terminus, so this is the common case.
5. **Perpendicular cross-sections** — the payoff. An ordered path has a tangent,
   so the lumen can be cut ⟂ to it for **area, both extents and aspect**. That
   is what §8A.2 could not measure by any route from the projected footprint,
   and it is what the wall-shear estimate needs: τ depends on the cross-section,
   and width and height enter it with different exponents.

**Pure functions, numpy in and numpy/dict out.** No Qt, no napari, no `state` —
the same doctrine as `saved_check.py` and `mesh_score.py`, so the GUI tab and
`batch_mesh` call one implementation and cannot drift.

**Nothing here chooses.** Every quantity is returned with the flags needed to
decide whether to trust it (`censored`, `section_clipped`, `bridged`); which
points to drop belongs to whoever runs the analysis, exactly as `mesh_score`
returns four numbers and never a total.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from . import shear as _shear

# --------------------------------------------------------------------------- #
# Defaults. Mirrored by `batch_mesh.GRAPH_DEFAULTS` and the napari sliders.
# --------------------------------------------------------------------------- #

# A leaf branch shorter than this many times the local inscribed radius is a
# spur. **Stated as a ratio on purpose** ([P23]): the vessels here taper, and an
# absolute length that is right at a 40 µm caliber is wrong at image 6's 5 µm
# waist. 1.0 reads "a stub that does not reach as far as the vessel is wide".
DEFAULT_SPUR_RATIO = 1.0
# How close to a field-of-view **opening** an endpoint must be to count as a
# crop rather than a real terminus, as a multiple of the local inscribed radius
# plus an absolute floor.
#
# **Why a ratio and not just a margin.** The medial axis of a tube cut off by a
# flat plane stops about one radius short of the cut — the largest inscribed
# ball near the cut is centred inward — so on these ~40 µm-radius vessels a
# crop endpoint sits 30-40 µm from the face it exits through. A small absolute
# margin therefore calls every crop a terminus, which is what the first version
# of this module did: 4 termini and 0 crops on an image whose lumen meets all
# four x/y faces with thousands of voxels. 2.0 leaves room for that offset
# without reaching across the vessel.
DEFAULT_CROP_RATIO = 2.0
DEFAULT_CROP_MARGIN_UM = 6.0
# Arc-length window for the 3D path curvature fit. κ is a second derivative and
# is meaningless below the caliber; 40 µm is ~one diameter on these vessels.
DEFAULT_KAPPA_WINDOW_UM = 40.0
# Cross-sections are the expensive step, so they are cut every this many µm and
# left NaN in between rather than silently interpolated.
DEFAULT_SECTION_EVERY_UM = 5.0
# The cut plane is sampled out to this many times the local inscribed radius.
# The inscribed radius is the *minor* semi-axis (§7.3), and these ducts run to
# aspect 3, so a disc of 1x caliber would clip every section it measured.
DEFAULT_SECTION_SPAN = 3.5

_NB27 = np.ones((3, 3, 3), np.uint8)


# --------------------------------------------------------------------------- #
# (1) the lumen volume
# --------------------------------------------------------------------------- #

def lumen_volume(chi: np.ndarray, iso_level: float,
                 core: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    """The lumen as ``chi > iso_level``, reduced to the component(s) that are
    actually lumen.

    ``chi`` is higher inside the lumen than outside — `reconstruct` checks that
    and negates the field if it came out inverted — so the super-level set is
    the lumen side. It can still contain pockets outside the vessel where the
    solve ran high, so the result is restricted to components meeting ``core``,
    the z-enclosed fragments `vessel_poisson.lumen_core` found. That is the same
    discriminator the orientation step trusts, reused rather than re-derived.

    With no ``core`` (or no overlap with it) the largest component is kept and
    the fallback is recorded in ``notes`` — never taken silently, because
    "largest component" is a guess and z-enclosure is a measurement.
    """
    info: dict = {"notes": []}
    chi = np.asarray(chi)
    raw = chi > float(iso_level)
    info["n_voxels_super_level"] = int(raw.sum())
    if not raw.any():
        info["notes"].append(
            "WARN: nothing lies above the iso level — no lumen.")
        return np.zeros(chi.shape, bool), info

    lab, n = ndi.label(raw)
    info["n_components_super_level"] = int(n)
    keep = None
    if core is not None and np.asarray(core).any():
        hit = np.unique(lab[np.asarray(core) & raw])
        hit = hit[hit > 0]
        if hit.size:
            keep = hit
            info["selected_by"] = "z-enclosed core overlap"
    if keep is None:
        sizes = np.bincount(lab.ravel())
        sizes[0] = 0
        keep = np.array([int(np.argmax(sizes))])
        info["selected_by"] = "largest component (fallback)"
        info["notes"].append(
            "WARN: no super-level component meets the z-enclosed lumen core, "
            "so the largest was taken instead. The core is a measurement and "
            "this is a guess — check the orientation before trusting anything "
            "downstream of it.")
    lumen = np.isin(lab, keep)
    info["n_components_kept"] = int(len(keep))
    info["n_voxels_lumen"] = int(lumen.sum())
    return lumen, info


def caliber_field(lumen: np.ndarray, h_um: float) -> np.ndarray:
    """Inscribed radius in µm at **every** lumen voxel.

    The distance to the nearest non-lumen voxel, i.e. the radius of the largest
    sphere centred there that fits inside the lumen. Being a sphere it is capped
    by the shortest direction, so on these flattened ducts it is the **minor**
    semi-axis and not the caliber that τ depends on (§7.3, and §8A.2 measured
    2R ≈ the lumen height on four images). It is the right thing to prune spurs
    against and the wrong thing to compute shear from — use the cross-sections
    for that.
    """
    return ndi.distance_transform_edt(lumen, sampling=(h_um,) * 3)


# --------------------------------------------------------------------------- #
# (2) the skeleton and its graph
# --------------------------------------------------------------------------- #

def curve_skeleton(lumen: np.ndarray) -> np.ndarray:
    """A one-voxel-wide 3D curve skeleton of the lumen (Lee's algorithm)."""
    from skimage.morphology import skeletonize
    if not lumen.any():
        return np.zeros(lumen.shape, bool)
    return np.asarray(skeletonize(lumen, method="lee") > 0)


def _degree(skel: np.ndarray) -> np.ndarray:
    """26-connected neighbour count at every skeleton voxel (0 elsewhere)."""
    cnt = ndi.convolve(skel.astype(np.uint8), _NB27, mode="constant")
    return np.where(skel, cnt - skel.astype(np.uint8), 0).astype(np.int16)


def _trace_edges(skel: np.ndarray, node_lab: np.ndarray):
    """Walk every degree-2 chain between node clusters.

    Returns ``(edges, paths)`` where ``edges[i] = (node_a, node_b)`` as 0-based
    node ids and ``paths[i]`` is the ordered ``(L,3)`` array of voxel indices
    from ``node_a`` to ``node_b``, node voxels included at both ends.

    A closed loop carrying no node at all is emitted with ``node_a == node_b ==
    -1`` and its first voxel repeated at the end, so it is still an ordered
    path with a defined arc length. The graph is a description of what the
    skeleton is, not an assertion that it is a tree.
    """
    shape = skel.shape
    off = np.array([(dz, dy, dx)
                    for dz in (-1, 0, 1) for dy in (-1, 0, 1)
                    for dx in (-1, 0, 1) if (dz, dy, dx) != (0, 0, 0)])

    def nbrs(p):
        q = np.asarray(p) + off
        ok = np.all((q >= 0) & (q < np.asarray(shape)), axis=1)
        q = q[ok]
        return [tuple(int(v) for v in r) for r in q if skel[tuple(r)]]

    edges, paths = [], []
    visited = set()          # traversed non-node voxels

    node_vox: dict[int, list] = {}
    for p in zip(*np.nonzero(node_lab)):
        node_vox.setdefault(int(node_lab[p]) - 1, []).append(p)

    for nid, vox in node_vox.items():
        for p in vox:
            for q in nbrs(p):
                if node_lab[q]:
                    # Node-to-node contact: only emit once, and only between
                    # two *different* clusters (inside one cluster it is not an
                    # edge, it is the cluster).
                    other = int(node_lab[q]) - 1
                    if other == nid:
                        continue
                    key = ("nn", min(nid, other), max(nid, other), p, q)
                    rkey = ("nn", min(nid, other), max(nid, other), q, p)
                    if key in visited or rkey in visited:
                        continue
                    visited.add(key)
                    edges.append((nid, other))
                    paths.append(np.array([p, q], np.int64))
                    continue
                if q in visited:
                    continue
                # Walk the degree-2 chain out of this node.
                chain = [p, q]
                visited.add(q)
                cur, prev = q, p
                while True:
                    nxt = [r for r in nbrs(cur)
                           if r != prev and not (r in visited and
                                                 not node_lab[r])]
                    nxt = [r for r in nxt if r != prev]
                    step = None
                    for r in nxt:
                        if node_lab[r]:
                            step = r
                            break
                    if step is None:
                        cand = [r for r in nxt if r not in visited]
                        step = cand[0] if cand else None
                    if step is None:
                        edges.append((nid, -1))      # dangling; should not
                        paths.append(np.array(chain, np.int64))
                        break
                    chain.append(step)
                    if node_lab[step]:
                        edges.append((nid, int(node_lab[step]) - 1))
                        paths.append(np.array(chain, np.int64))
                        break
                    visited.add(step)
                    prev, cur = cur, step

    # Closed loops: skeleton voxels never touched and carrying no node.
    if node_lab.any():
        seen = np.zeros(shape, bool)
        for pth in paths:
            seen[pth[:, 0], pth[:, 1], pth[:, 2]] = True
    else:
        seen = np.zeros(shape, bool)
    leftover = skel & ~seen & (node_lab == 0)
    if leftover.any():
        lab, n = ndi.label(leftover, structure=_NB27)
        for k in range(1, n + 1):
            vox = list(zip(*np.nonzero(lab == k)))
            if len(vox) < 3:
                continue
            chain = [vox[0]]
            used = {vox[0]}
            while True:
                nx = [r for r in nbrs(chain[-1])
                      if lab[r] == k and r not in used]
                if not nx:
                    break
                chain.append(nx[0])
                used.add(nx[0])
            chain.append(chain[0])
            edges.append((-1, -1))
            paths.append(np.array(chain, np.int64))

    return edges, paths


def skeleton_graph(skel: np.ndarray):
    """Nodes and ordered edges of a 3D curve skeleton.

    Node voxels are those of degree != 2 (endpoints at degree 1, junctions at
    degree >= 3). Adjacent node voxels are collapsed into **one** node — a
    26-connected skeleton renders a single Y as a small clump rather than a
    lone voxel, and counting each voxel as a junction is exactly how the old
    2D code came to report six branch nodes on a vessel with one (§8A.3).
    """
    deg = _degree(skel)
    nodes = skel & (deg != 2)
    node_lab, n_nodes = ndi.label(nodes, structure=_NB27)
    edges, paths = _trace_edges(skel, node_lab)
    centres = (np.array(ndi.center_of_mass(nodes, node_lab,
                                           range(1, n_nodes + 1)))
               if n_nodes else np.zeros((0, 3)))
    # A node's degree is how many edges leave it, not how many voxels it has.
    ndeg = np.zeros(n_nodes, int)
    for a, b in edges:
        if a >= 0:
            ndeg[a] += 1
        if b >= 0:
            ndeg[b] += 1
    return {"node_idx": centres, "node_degree": ndeg,
            "edges": edges, "paths": paths}


# --------------------------------------------------------------------------- #
# (3) pruning, in units of the vessel's own caliber
# --------------------------------------------------------------------------- #

def prune_spurs(graph: dict, caliber: np.ndarray, h_um: float,
                spur_ratio: float = DEFAULT_SPUR_RATIO) -> tuple[dict, dict]:
    """Drop leaf edges shorter than ``spur_ratio`` x the local caliber.

    The comparison is against the **median inscribed radius along the branch
    itself**, so a stub off a 40 µm-wide arm and a stub off a 5 µm waist are
    judged by their own surroundings. That is the thing [P23] says an absolute
    threshold cannot do, and the reason the waist of image 6 — the highest-shear
    region in the dataset and the one `open_um`/`prune_um` amputate — survives
    this rule.

    Iterated, because removing a spur can turn a degree-3 junction into a
    degree-2 pass-through whose two edges should then be one.
    """
    info = {"notes": [], "n_removed": 0, "removed_um": []}
    for _ in range(32):
        edges, paths = graph["edges"], graph["paths"]
        deg = graph["node_degree"]
        drop = None
        for i, ((a, b), pth) in enumerate(zip(edges, paths)):
            leaf = [k for k in (a, b) if k >= 0 and deg[k] == 1]
            if not leaf or len(edges) == 1:
                continue
            length = _polyline_length(pth * h_um)
            rad = float(np.median(caliber[pth[:, 0], pth[:, 1], pth[:, 2]]))
            if length < spur_ratio * max(rad, h_um):
                drop = (i, length, rad)
                break
        if drop is None:
            break
        i, length, rad = drop
        info["removed_um"].append(round(float(length), 2))
        info["n_removed"] += 1
        graph = _drop_edge(graph, i)
    else:
        info["notes"].append(
            "WARN: spur pruning hit its 32-pass limit; the graph may still "
            "carry stubs.")
    if info["n_removed"]:
        info["notes"].append(
            f"NOTE: pruned {info['n_removed']} spur(s) shorter than "
            f"{spur_ratio:g}x the local inscribed radius.")
    return graph, info


def simplify(graph: dict) -> tuple[dict, dict]:
    """Collapse degree-2 nodes and drop nodes no edge touches.

    Two things the raw trace leaves behind, both of which corrupt arc length if
    they survive:

    - **Pass-through nodes.** A voxel with three 26-neighbours is not
      necessarily a junction — a diagonal bend makes one — and pruning a spur
      turns a real junction into a pass-through. Either way the node splits one
      vessel into two edges, so ``s_um`` restarts in the middle of a straight
      run and "arc length along the vessel" stops meaning that. The two edges
      are joined back into one.
    - **Orphans.** A node whose every edge was pruned. Harmless to the geometry
      and actively misleading in a count: the first run of this module reported
      an endpoint at degree 0.

    Nodes surviving here are junctions (degree >= 3) and endpoints (degree 1)
    and nothing else, which is what §8B asks a vessel graph to be.
    """
    info = {"n_merged": 0, "n_orphans_dropped": 0, "notes": []}
    for _ in range(4096):
        deg = graph["node_degree"]
        edges, paths = graph["edges"], graph["paths"]
        cand = None
        for n in range(len(deg)):
            if deg[n] != 2:
                continue
            touch = [i for i, (a, b) in enumerate(edges) if n in (a, b)]
            # A node both of whose edges are the same edge is a loop pinned at
            # one point; joining it would produce a path with no ends.
            if len(touch) == 2 and touch[0] != touch[1]:
                cand = (n, touch[0], touch[1])
                break
        if cand is None:
            break
        n, i, j = cand
        pi = paths[i] if edges[i][1] == n else paths[i][::-1]
        oi = edges[i][0] if edges[i][1] == n else edges[i][1]
        pj = paths[j] if edges[j][0] == n else paths[j][::-1]
        oj = edges[j][1] if edges[j][0] == n else edges[j][0]
        merged = np.concatenate([pi, pj[1:]], axis=0)
        keep = [k for k in range(len(edges)) if k not in (i, j)]
        new_edges = [edges[k] for k in keep] + [(oi, oj)]
        new_paths = [paths[k] for k in keep] + [merged]
        graph = _rebuild(graph["node_idx"], new_edges, new_paths)
        info["n_merged"] += 1
    else:
        info["notes"].append(
            "WARN: node simplification hit its pass limit; the graph may still "
            "carry degree-2 nodes.")

    deg = graph["node_degree"]
    live = deg > 0
    info["n_orphans_dropped"] = int((~live).sum())
    if not live.all():
        remap = np.full(len(deg), -1, int)
        remap[live] = np.arange(int(live.sum()))
        edges = [(int(remap[a]) if a >= 0 else -1,
                  int(remap[b]) if b >= 0 else -1)
                 for a, b in graph["edges"]]
        graph = _rebuild(np.asarray(graph["node_idx"])[live], edges,
                         graph["paths"])
    if info["n_merged"] or info["n_orphans_dropped"]:
        info["notes"].append(
            f"NOTE: simplified the graph — merged {info['n_merged']} "
            f"degree-2 node(s), dropped {info['n_orphans_dropped']} orphan(s).")
    return graph, info


def _rebuild(node_idx, edges, paths) -> dict:
    n = len(node_idx)
    deg = np.zeros(n, int)
    for a, b in edges:
        if a >= 0:
            deg[a] += 1
        if b >= 0:
            deg[b] += 1
    return {"node_idx": np.asarray(node_idx), "node_degree": deg,
            "edges": list(edges), "paths": list(paths)}


def face_openings(lumen: np.ndarray) -> tuple[np.ndarray, dict]:
    """Where the lumen leaves the volume through an x/y face.

    Returns a boolean array marking lumen voxels lying on the x or y domain
    faces, and a record of how much lumen sits on each of the six faces.

    **x/y only, deliberately.** The vessel runs in the imaging plane and leaves
    through the side walls; `vessel_poisson`'s Neumann boundary is exactly the
    statement that it may. The z faces are different — [P24](b) measured the
    vessel to be *contained* in z (0.01% of the mask in the first slice, none in
    the last), and that containment is what `lumen_core`'s z-enclosure
    discriminator and the whole normal orientation rest on. So lumen on a z face
    is not an opening to be censored around, it is a **warning that the premise
    failed**, and it is reported as one. Image 4 is the known case (1.224% of
    its mask in the last z slice, 150x the reference).
    """
    m = np.zeros(lumen.shape, bool)
    m[:, 0, :] |= lumen[:, 0, :]
    m[:, -1, :] |= lumen[:, -1, :]
    m[:, :, 0] |= lumen[:, :, 0]
    m[:, :, -1] |= lumen[:, :, -1]
    info = {"y_lo": int(lumen[:, 0, :].sum()), "y_hi": int(lumen[:, -1, :].sum()),
            "x_lo": int(lumen[:, :, 0].sum()), "x_hi": int(lumen[:, :, -1].sum()),
            "z_lo": int(lumen[0].sum()), "z_hi": int(lumen[-1].sum()),
            "notes": []}
    xy = info["y_lo"] + info["y_hi"] + info["x_lo"] + info["x_hi"]
    z = info["z_lo"] + info["z_hi"]
    info["xy_face_voxels"] = int(xy)
    info["z_face_voxels"] = int(z)
    # Some lumen on a z face is expected and is **not** the vessel leaving:
    # cadherin attenuates with depth ([P03]), so the deep wall is missing from
    # the mask, and χ's super-level set bleeds through the hole it leaves. The
    # premise that actually matters — that the vessel is contained in z, which
    # `lumen_core` and every normal orientation rest on ([P24](b)) — fails only
    # when the lumen leaves mainly through z rather than through the sides.
    if z and z <= xy:
        info["notes"].append(
            f"NOTE: {z:,} lumen voxels sit on a z face against {xy:,} on the "
            f"x/y faces. Expected — the deep wall is under-detected ([P03]) so "
            f"the fitted lumen bleeds through it — and not counted as an "
            f"opening. z faces are excluded from the crop test either way.")
    elif z > xy:
        info["notes"].append(
            f"WARN: the lumen leaves mainly through the **z** faces "
            f"({z:,} voxels against {xy:,} on x/y). [P24](b) measures this "
            f"vessel family as contained in z, and `lumen_core`'s orientation "
            f"depends on that — treat this image's normals, and therefore the "
            f"sign of its curvature, as the least well-founded of the set.")
    return m, info


def _drop_edge(graph: dict, i: int) -> dict:
    edges = [e for k, e in enumerate(graph["edges"]) if k != i]
    paths = [p for k, p in enumerate(graph["paths"]) if k != i]
    n = len(graph["node_idx"])
    deg = np.zeros(n, int)
    for a, b in edges:
        if a >= 0:
            deg[a] += 1
        if b >= 0:
            deg[b] += 1
    return {"node_idx": graph["node_idx"], "node_degree": deg,
            "edges": edges, "paths": paths}


def _polyline_length(pts_um: np.ndarray) -> float:
    if len(pts_um) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(pts_um, axis=0), axis=1).sum())


# --------------------------------------------------------------------------- #
# (4) geometry along an ordered path
# --------------------------------------------------------------------------- #

def _resample(pts_um: np.ndarray, step_um: float) -> np.ndarray:
    """Uniform arc-length resampling of an ordered polyline.

    The skeleton's voxel spacing alternates between h, h*sqrt2 and h*sqrt3
    depending on whether a step is face-, edge- or corner-diagonal, so a
    finite-difference derivative on the raw voxel sequence carries that
    staircase straight into κ. Resampling first is what makes the derivative a
    derivative of the *curve*.
    """
    if len(pts_um) < 2:
        return pts_um.astype(np.float64)
    seg = np.linalg.norm(np.diff(pts_um, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(s[-1])
    if total <= 0:
        return pts_um[:1].astype(np.float64)
    m = max(2, int(round(total / step_um)) + 1)
    su = np.linspace(0.0, total, m)
    return np.column_stack([np.interp(su, s, pts_um[:, k]) for k in range(3)])


def path_curvature(pts_um: np.ndarray, window_um: float):
    """3D path curvature κ = |r' x r''| / |r'|³ along an ordered polyline.

    Returned in 1/µm together with the arc length of every point and the unit
    tangent. **This is the quantity `path_curv_xy` was trying to be** and could
    not: the old column fitted a parabola to a *Euclidean* KD-tree ball over an
    unordered point set, so near a junction it fitted across the gap into
    another branch ([P08](a)), and there was no arc length to fit against at all
    ([P27]). Here the neighbourhood is a window along the path, which cannot
    reach into a branch it is not on.

    Smoothing is a Gaussian along arc length with sigma = ``window_um`` / 4, so
    the scale is explicit and travels with the number — [P08](e)'s complaint is
    that a 20 µm window on a 140 µm vessel measures medial-axis wobble, and the
    honest response is to state the scale, not to pick a better constant.
    """
    p = np.asarray(pts_um, float)
    if len(p) < 3:
        return (np.zeros(len(p)), np.zeros(len(p)),
                np.zeros((len(p), 3)))
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    # `np.gradient` divides by the coordinate spacing, so a zero-length segment
    # is a divide-by-zero. Uniformly-resampled 3D paths have none, but the
    # in-plane variant zeroes the z coordinate, and a stretch of vessel running
    # straight down through the stack projects to a repeated xy point. Nudging
    # to strictly increasing costs nothing where the segments are real and
    # turns a NaN into a (correctly) near-zero curvature where they are not.
    s = np.maximum.accumulate(s) + np.arange(len(s)) * 1e-9
    ds = float(np.median(seg)) if len(seg) else 1.0
    sigma = max(window_um / 4.0 / max(ds, 1e-9), 0.5)
    sm = np.column_stack([ndi.gaussian_filter1d(p[:, k], sigma,
                                                mode="nearest")
                          for k in range(3)])
    d1 = np.gradient(sm, s, axis=0)
    d2 = np.gradient(d1, s, axis=0)
    num = np.linalg.norm(np.cross(d1, d2), axis=1)
    den = np.linalg.norm(d1, axis=1) ** 3
    kappa = np.where(den > 1e-12, num / np.maximum(den, 1e-12), 0.0)
    tang = d1 / (np.linalg.norm(d1, axis=1, keepdims=True) + 1e-12)
    return kappa, s, tang


# --------------------------------------------------------------------------- #
# (5) the payoff: a cut perpendicular to the axis
# --------------------------------------------------------------------------- #

def cross_section(lumen: np.ndarray, origin_um: np.ndarray, h_um: float,
                  point_um: np.ndarray, tangent: np.ndarray,
                  radius_um: float, span: float = DEFAULT_SECTION_SPAN,
                  step_um: float | None = None,
                  extras: dict | None = None,
                  field: np.ndarray | None = None,
                  level: float | None = None) -> dict:
    """Cut the lumen with the plane through ``point_um`` normal to ``tangent``.

    Returns area, the two principal extents, aspect, and whether the section ran
    into the edge of the sampling disc or the edge of the volume. **This is the
    measurement §8A.2 could not make by any route from the max-projected
    footprint**, and the one §4's shear formula needs: for a flat duct
    τ ≈ 6 µQ / (w h²), so width and height enter with different exponents and
    must not be collapsed into one caliber.

    ``major``/``minor`` are the full **extents** along the section's principal
    axes, not ellipse semi-axes fitted to it — the measured cross-sections are
    D-shaped rather than elliptical ([P24](b)), and an extent says what it says
    on a shape that is not an ellipse.

    ``perim_um`` is the sub-pixel wall length, from `shear.section_contour`, and
    is what `D_h = 4A/P` needs — §4.2's algebra could not be evaluated
    numerically anywhere in this project before it existed.

    **``field`` is χ and ``level`` its iso-level — hand them over.** Without
    them the plane is cut out of the *binarised* lumen with nearest-neighbour
    sampling, and every wall in the section lands on a pixel edge. That
    staircase is not a rounding error: it makes the wall 4-5% longer than the
    surface it stands for at **every** resolution, and the shear solved inside
    it inherits that ([P46](docs/ISSUES.md#p46)). With them, the plane is cut
    out of χ by linear interpolation and the mask is `plane > level`, so the
    section carries both the region *and* the smooth field whose level curve is
    its true wall — which is the same curve the displayed mesh is made of, since
    marching cubes interpolates the same χ at the same level. `shear` then puts
    the wall where the field says rather than at the pixel edge, and the disc
    test goes from 4.4% long and 4% low in τ to **0.03% and 0.8%**, converging
    with step where the old one did not converge at all.

    ``extras``, if given, is filled with the section **mask** and the frame it
    lives in (``u``, ``v``, ``step_um``, ``center_rc``, ``point_um``,
    ``tangent``), plus ``field``/``level`` when they were supplied. That is what `shear.solve_section` runs on and what
    `shear.contour_to_3d` needs to put a wall point back in the image frame.
    It is an optional out-dict for the same reason `vessel_mesh.build_mesh`
    uses one: the return value has consumers that destructure it, and a mask is
    orders of magnitude larger than the scalars beside it.
    """
    step = float(step_um if step_um else h_um / 2.0)
    R = float(span) * float(radius_um)
    out = {"area_um2": np.nan, "major_um": np.nan, "minor_um": np.nan,
           "aspect": np.nan, "perim_um": np.nan,
           "clipped": True, "touches_domain": False}
    if not np.isfinite(R) or R <= 0:
        return out

    t = np.asarray(tangent, float)
    nt = np.linalg.norm(t)
    if nt < 1e-9:
        return out
    t = t / nt
    # Any vector not parallel to t; the least-aligned axis is the safe pick.
    u = np.zeros(3)
    u[int(np.argmin(np.abs(t)))] = 1.0
    u = np.cross(t, u)
    u /= np.linalg.norm(u)
    v = np.cross(t, u)

    m = int(np.ceil(R / step))
    g = np.arange(-m, m + 1) * step
    A, B = np.meshgrid(g, g, indexing="ij")
    pos = (np.asarray(point_um, float)[None, None, :]
           + A[..., None] * u[None, None, :]
           + B[..., None] * v[None, None, :])
    idx = (pos - np.asarray(origin_um, float)[None, None, :]) / h_um
    inside_vol = np.all((idx >= 0) &
                        (idx <= np.asarray(lumen.shape) - 1), axis=-1)
    plane = None
    if field is not None and level is not None:
        # Linear interpolation, and a cval well below the level so anything
        # sampled outside the volume reads as solidly outside rather than as
        # whatever `cval=0` happens to mean for this field.
        lo = float(np.nanmin(field)) - 1.0
        plane = ndi.map_coordinates(
            np.asarray(field, float),
            np.moveaxis(idx, -1, 0).reshape(3, -1),
            order=1, mode="constant", cval=lo).reshape(A.shape)
        sec = plane > float(level)
    else:
        samp = ndi.map_coordinates(lumen.astype(np.uint8),
                                   np.moveaxis(idx, -1, 0).reshape(3, -1),
                                   order=0, mode="constant", cval=0)
        sec = samp.reshape(A.shape).astype(bool)
    if not sec[m, m]:
        # The centre is not in the lumen: the point sits on a skeleton voxel
        # that the plane cut obliquely out of. Take the nearest component.
        lab, n = ndi.label(sec)
        if n == 0:
            return out
        yy, xx = np.nonzero(sec)
        d = (yy - m) ** 2 + (xx - m) ** 2
        comp = lab[yy[np.argmin(d)], xx[np.argmin(d)]]
    else:
        lab, n = ndi.label(sec)
        comp = lab[m, m]
    sel = lab == comp
    yy, xx = np.nonzero(sel)
    if len(yy) < 3:
        return out

    plane_sel = None
    if plane is not None:
        # Only the selected component's wall. The field is blanked outside a
        # one-pixel dilation of `sel`, so the level curve just outside the
        # component's boundary pixels survives while a neighbouring component's
        # does not. Dilating by exactly one is the smallest window that keeps
        # the crossing and the largest that cannot reach a component a pixel
        # away.
        lo = float(np.nanmin(plane)) - 1.0
        plane_sel = np.where(ndi.binary_dilation(sel), plane, lo)

    if extras is not None:
        extras.update({"mask": sel, "step_um": step, "u": u, "v": v,
                       "center_rc": (float(m), float(m)),
                       "point_um": np.asarray(point_um, float).copy(),
                       "tangent": t.copy(),
                       "field": plane_sel,
                       "level": (None if plane_sel is None else float(level))})

    out["area_um2"] = float(len(yy) * step * step)
    out["perim_um"] = _shear.section_contour(
        sel, step, plane_sel,
        0.5 if plane_sel is None else float(level))["perim_um"]
    a = (yy - m) * step
    b = (xx - m) * step
    P = np.column_stack([a, b])
    P = P - P.mean(0)
    _, _, vt = np.linalg.svd(P, full_matrices=False)
    proj = P @ vt.T
    ext = proj.max(0) - proj.min(0)
    major, minor = float(max(ext)), float(min(ext))
    out["major_um"] = major
    out["minor_um"] = minor
    out["aspect"] = float(major / minor) if minor > 0 else np.nan
    # Did the section reach the edge of the disc we sampled, or the edge of the
    # volume? Either way the extent is a lower bound and must not be read as a
    # measurement. Both are flagged rather than silently dropped.
    edge = np.zeros_like(sel)
    edge[0, :] = edge[-1, :] = edge[:, 0] = edge[:, -1] = True
    r2 = (a ** 2 + b ** 2)
    out["clipped"] = bool((sel & edge).any()
                          or (r2.max() >= (R - step) ** 2))
    out["touches_domain"] = bool((sel & ~inside_vol).any())
    return out


# --------------------------------------------------------------------------- #
# (6) the whole thing
# --------------------------------------------------------------------------- #

def build_graph(chi: np.ndarray, iso_level: float, h_um: float,
                origin_um, *, core: np.ndarray | None = None,
                spur_ratio: float = DEFAULT_SPUR_RATIO,
                crop_ratio: float = DEFAULT_CROP_RATIO,
                crop_margin_um: float = DEFAULT_CROP_MARGIN_UM,
                kappa_window_um: float = DEFAULT_KAPPA_WINDOW_UM,
                section_from_field: bool = True,
                section_every_um: float = DEFAULT_SECTION_EVERY_UM,
                section_span: float = DEFAULT_SECTION_SPAN,
                compute_shear: bool = False,
                domain_shape_um=None, extras: dict | None = None) -> dict:
    """χ in, the whole centerline graph out.

    Every returned per-point array is parallel and in path order, concatenated
    edge by edge — ``edge_id`` says which edge a point belongs to and ``s_um``
    its arc length from that edge's start. **That ordering is the thing [P27]
    is about**, and every consumer may now rely on it.

    ``domain_shape_um`` is the extent of the *original image* in µm, used to
    decide whether an endpoint is a crop. It defaults to the χ grid's own
    extent, which is very slightly smaller.

    ``extras``, if given, receives ``extras["wall"]`` — the **per-boundary-point**
    shear profile, when ``compute_shear`` is on. That is a ragged quantity (a
    different number of wall points per section), so it cannot live in the
    parallel per-point columns and it goes in an optional out-dict instead, the
    same doctrine `vessel_mesh.build_mesh` uses. It is **the** shear
    measurement: the per-station columns are a summary of it. Two consumers want
    it and they want the same numbers — colouring the vessel surface, and
    joining an adhesion to the shear at its nearest luminal wall point.

    ``compute_shear`` adds §4.3's per-section wall-shear solve (`shear.py`).

    ``section_from_field`` cuts the cross-sections out of **χ** rather than out
    of the binarised lumen — see `cross_section`. On by default: it is the same
    surface either way, but the binary cut throws away the sub-voxel precision
    the Poisson solve already produced, and [P46](docs/ISSUES.md#p46) is the
    size of what it throws away. Set it False to reproduce a pre-2026-08-28
    build.
    **Default off, and that is a decision**: `batch.py` and the tab run at
    defaults, so turning it on here would change every headless run and
    invalidate the artifacts already produced. It is opt-in for the same reason
    `batch_mesh --steps score` is — it costs more than the thing it annotates.
    """
    origin_um = np.broadcast_to(np.asarray(origin_um, float), (3,)).copy()
    out: dict = {"notes": [], "params": {
        "iso_level": float(iso_level), "h_um": float(h_um),
        "origin_um": origin_um.tolist(), "spur_ratio": float(spur_ratio),
        "crop_ratio": float(crop_ratio),
        "crop_margin_um": float(crop_margin_um),
        "kappa_window_um": float(kappa_window_um),
        "section_every_um": float(section_every_um),
        "section_from_field": bool(section_from_field),
        "section_span": float(section_span),
        "compute_shear": bool(compute_shear)}}

    lumen, linfo = lumen_volume(chi, iso_level, core)
    out["lumen_info"] = linfo
    out["notes"].extend(linfo["notes"])
    out["lumen_volume_um3"] = float(lumen.sum() * h_um ** 3)
    if not lumen.any():
        return _empty_graph(out)

    caliber = caliber_field(lumen, h_um)
    skel = curve_skeleton(lumen)
    if skel.sum() < 3:
        out["notes"].append(
            "WARN: the lumen skeletonized to fewer than 3 voxels — no "
            "centerline. The lumen is probably a blob rather than a tube.")
        return _empty_graph(out)

    g = skeleton_graph(skel)
    g, pinfo = prune_spurs(g, caliber, h_um, spur_ratio)
    out["notes"].extend(pinfo["notes"])
    out["prune"] = pinfo
    g, sinfo = simplify(g)
    out["notes"].extend(sinfo["notes"])
    out["simplify"] = sinfo

    openings, oinfo = face_openings(lumen)
    out["face_openings"] = oinfo
    out["notes"].extend(oinfo["notes"])
    # Distance from every voxel to the nearest lumen-on-an-x/y-face voxel.
    d_open = (ndi.distance_transform_edt(~openings, sampling=(h_um,) * 3)
              if openings.any() else
              np.full(lumen.shape, np.inf, np.float32))

    if domain_shape_um is None:
        domain_shape_um = (np.asarray(chi.shape) - 1) * h_um + origin_um

    # ---- nodes ----
    node_um = (np.asarray(g["node_idx"]) * h_um + origin_um
               if len(g["node_idx"]) else np.zeros((0, 3)))
    deg = g["node_degree"]
    kind = np.array(["junction" if d >= 3 else
                     ("endpoint" if d <= 1 else "pass") for d in deg],
                    dtype="<U9")
    # An endpoint near a field-of-view **opening** is where the image cut the
    # vessel, not where the vessel ends. §5.1: no image in this dataset appears
    # to contain a real terminus, so "crop" should be the common answer and a
    # "terminus" is the result that needs explaining.
    #
    # Measured against the opening itself and scaled by the local caliber, not
    # against the bounding plane by a fixed margin — see DEFAULT_CROP_RATIO for
    # why the fixed margin called every crop a terminus.
    gi_n = (np.clip(np.rint((node_um - origin_um) / h_um).astype(int), 0,
                    np.asarray(lumen.shape) - 1)
            if len(node_um) else np.zeros((0, 3), int))
    node_cal = (caliber[gi_n[:, 0], gi_n[:, 1], gi_n[:, 2]]
                if len(node_um) else np.zeros(0))
    node_dopen = (d_open[gi_n[:, 0], gi_n[:, 1], gi_n[:, 2]]
                  if len(node_um) else np.zeros(0))
    thresh = np.maximum(crop_margin_um, crop_ratio * node_cal)
    at_face = node_dopen <= thresh
    end_kind = np.where(kind != "endpoint", "",
                        np.where(at_face, "crop", "terminus")).astype("<U9")
    out["node_dist_opening_um"] = node_dopen.astype(np.float32)
    out["node_caliber_um"] = node_cal.astype(np.float32)

    # ---- per-edge geometry ----
    P, E, S, K, KXY, CAL, SEC = [], [], [], [], [], [], []
    WALL: list = []
    edge_len, edge_nodes, edge_censored = [], [], []
    for ei, ((a, b), pth) in enumerate(zip(g["edges"], g["paths"])):
        raw_um = pth.astype(float) * h_um + origin_um
        pts = _resample(raw_um, max(h_um / 2.0, 1e-6))
        if len(pts) < 2:
            continue
        kappa, s, tang = path_curvature(pts, kappa_window_um)
        # The in-plane component, so the new column is readable against the old
        # `path_curv_xy` even though it is not the same quantity.
        kxy, _s2, _t2 = path_curvature(
            np.column_stack([np.zeros(len(pts)), pts[:, 1], pts[:, 2]]),
            kappa_window_um)
        gi = np.clip(np.rint((pts - origin_um) / h_um).astype(int), 0,
                     np.asarray(lumen.shape) - 1)
        cal = caliber[gi[:, 0], gi[:, 1], gi[:, 2]]

        cens = bool((a >= 0 and end_kind[a] == "crop")
                    or (b >= 0 and end_kind[b] == "crop"))
        wall_rows: list | None = [] if (compute_shear
                                        and extras is not None) else None
        sec_rows = _sections_along(lumen, origin_um, h_um, pts, tang, cal, s,
                                   section_every_um, section_span,
                                   compute_shear, wall_rows,
                                   field=(chi if section_from_field else None),
                                   level=(float(iso_level) if section_from_field
                                          else None))
        if wall_rows is not None:
            WALL.append((sum(len(x) for x in P), wall_rows))
        P.append(pts)
        E.append(np.full(len(pts), ei, np.int32))
        S.append(s)
        K.append(kappa)
        KXY.append(kxy)
        CAL.append(cal)
        SEC.append(sec_rows)
        edge_len.append(float(s[-1]))
        edge_nodes.append((int(a), int(b)))
        edge_censored.append(cens)

    if not P:
        return _empty_graph(out)

    out["points_zyx_um"] = np.concatenate(P).astype(np.float32)
    out["edge_id"] = np.concatenate(E)
    out["s_um"] = np.concatenate(S).astype(np.float32)
    out["kappa_1_um"] = np.concatenate(K).astype(np.float32)
    out["kappa_xy_1_um"] = np.concatenate(KXY).astype(np.float32)
    out["caliber_um"] = np.concatenate(CAL).astype(np.float32)
    sec = np.concatenate(SEC)
    for name in _SEC_FLOATS + _SHEAR_FLOATS:
        out[name] = sec[name].astype(np.float32)
    out["section_clipped"] = sec["clipped"]
    out["section_touches_domain"] = sec["touches_domain"]

    out["nodes_zyx_um"] = node_um.astype(np.float32)
    out["node_degree"] = np.asarray(deg, np.int32)
    out["node_kind"] = kind
    out["node_endpoint_kind"] = end_kind
    out["edges"] = np.asarray(edge_nodes, np.int32).reshape(-1, 2)
    out["edge_len_um"] = np.asarray(edge_len, np.float32)
    out["edge_censored"] = np.asarray(edge_censored, bool)

    out["dist_junction_um"] = _dist_to_nodes(
        out["edge_id"], out["s_um"], out["edges"], out["edge_len_um"],
        np.nonzero(kind == "junction")[0], len(node_um))
    out["dist_crop_um"] = _dist_to_nodes(
        out["edge_id"], out["s_um"], out["edges"], out["edge_len_um"],
        np.nonzero(end_kind == "crop")[0], len(node_um))
    # **Censoring is a distance, not a label on the whole edge.** Marking every
    # point of an edge that eventually reaches a crop censors 100% of a Y whose
    # three arms all leave the frame — true, useless as a filter, and it would
    # throw away the junction, which is the one place the geometry is most
    # informative. A point is censored when it sits within `crop_ratio` calibers
    # of the cut, i.e. where the cross-section is being measured on a lumen the
    # image no longer fully contains.
    with np.errstate(invalid="ignore"):
        out["censored"] = np.asarray(
            out["dist_crop_um"] <= crop_ratio * out["caliber_um"], bool)
    out["censored"] &= np.isfinite(out["dist_crop_um"])

    # ---- the per-boundary-point wall profile ---------------------------------
    # Ragged, so it goes in `extras` rather than in the parallel columns. Each
    # wall point carries the **station index** it came from, which is how it
    # inherits that station's `section_clipped` and `censored` flags — the same
    # flags, not a second copy computed some other way. Nothing is filtered
    # here; `lumen_graph` chooses nothing.
    if extras is not None and compute_shear:
        PTS, TAU, SEG, QH, IDX = [], [], [], [], []
        for base, rows in WALL:
            for r in rows:
                n = len(r["tau_hat_um"])
                PTS.append(r["points_zyx_um"])
                TAU.append(r["tau_hat_um"])
                SEG.append(r["seg_len_um"])
                QH.append(np.full(n, r["q_hat_um4"]))
                IDX.append(np.full(n, base + r["i"], np.int64))
        if PTS:
            idx = np.concatenate(IDX)
            tau = np.concatenate(TAU)
            qh = np.concatenate(QH)
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = np.where(qh > 0, tau / qh, np.nan)
            wall = {
                "points_zyx_um": np.concatenate(PTS).astype(np.float32),
                "tau_hat_um": tau.astype(np.float32),
                "tau_rel_1_um3": rel.astype(np.float32),
                "q_hat_um4": qh.astype(np.float32),
                "seg_len_um": np.concatenate(SEG).astype(np.float32),
                "point_index": idx,
                "s_um": out["s_um"][idx],
                "edge_id": out["edge_id"][idx],
                "section_clipped": out["section_clipped"][idx],
                "censored": out["censored"][idx],
                # Carried per wall point so the file stands on its own. A
                # consumer of `tau_rel_1_um3` **must** be able to tell which
                # edge a point is on and how near a junction it is: τ/µQ is
                # normalised by conservation of Q, which holds along an
                # unbranched segment and fails at a Y, where the split rides
                # on h³ and is not recoverable (§4.2, [P49]). Reading those
                # out of `vessel_graph.npz` instead would work right up until
                # someone loads a shear file without its graph beside it.
                "dist_junction_um": out["dist_junction_um"][idx],
                "caliber_um": out["caliber_um"][idx],
            }
        else:
            wall = {"points_zyx_um": np.zeros((0, 3), np.float32),
                    "point_index": np.zeros(0, np.int64)}
            for k in ("tau_hat_um", "tau_rel_1_um3", "q_hat_um4", "seg_len_um",
                      "s_um", "dist_junction_um", "caliber_um"):
                wall[k] = np.zeros(0, np.float32)
            wall["edge_id"] = np.zeros(0, np.int32)
            wall["section_clipped"] = np.zeros(0, bool)
            wall["censored"] = np.zeros(0, bool)
        wall["max_dist_um"] = float(max(section_every_um, 2.0 * h_um))
        extras["wall"] = wall
        out["notes"].append(
            f"NOTE: {len(wall['points_zyx_um']):,} wall points carry a shear "
            f"measurement, over {int(np.unique(wall['point_index']).size)} "
            f"section(s). Sections are cut every {section_every_um:g} um, so "
            f"the surface field is that coarse along the axis — cut them "
            f"closer rather than interpolating if it matters.")

    n_junc = int((kind == "junction").sum())
    n_end = int((kind == "endpoint").sum())
    out["n_junctions"] = n_junc
    out["n_endpoints"] = n_end
    out["n_crop_endpoints"] = int((end_kind == "crop").sum())
    out["n_terminus_endpoints"] = int((end_kind == "terminus").sum())
    out["total_length_um"] = float(np.sum(edge_len))
    out["notes"].append(
        f"NOTE: {len(edge_len)} edge(s), {n_junc} junction(s), {n_end} "
        f"endpoint(s) of which {out['n_crop_endpoints']} are field-of-view "
        f"crops; {out['total_length_um']:.0f} µm of centerline.")
    if out["n_terminus_endpoints"]:
        out["notes"].append(
            f"NOTE: {out['n_terminus_endpoints']} endpoint(s) terminate inside "
            f"the volume rather than at an x/y face. No image in this dataset "
            f"is expected to contain a real vessel terminus (DIRECTION §5.1), "
            f"so check these against the raw channel before believing them.")
    return out


_SEC_FLOATS = ("area_um2", "major_um", "minor_um", "aspect", "perim_um")
# The shear columns, filled only when `build_graph(compute_shear=True)`. Left
# NaN otherwise rather than absent, so a consumer sees the same columns either
# way and a missing solve reads as "not measured" instead of as a KeyError.
_SHEAR_FLOATS = ("hydraulic_diam_um", "q_hat_um4", "tau_hat_mean_um",
                 "tau_hat_min_um", "tau_hat_max_um", "tau_rel_1_um3",
                 "force_balance", "wall_ok_frac")
_SEC_DTYPE = np.dtype([(n, "f8") for n in _SEC_FLOATS + _SHEAR_FLOATS]
                      + [("clipped", "?"), ("touches_domain", "?")])


def _sections_along(lumen, origin_um, h_um, pts, tang, cal, s,
                    every_um, span, compute_shear=False, wall_out=None,
                    field=None, level=None):
    """Cross-sections every ``every_um`` of arc length; NaN in between.

    Left NaN rather than interpolated on purpose: an interpolated cross-section
    is not a measurement, and this project has already been bitten by a number
    that looked measured and was not ([P26]).

    With ``compute_shear``, each section is also handed to `shear.solve_section`
    and the per-station summary of the wall-shear profile is recorded. The
    summary goes in the parallel columns; the **full per-boundary-point
    profile** is ragged, so it goes into ``wall_out`` — one dict per station —
    for `build_graph` to concatenate.
    """
    rows = np.zeros(len(pts), _SEC_DTYPE)
    for name in _SEC_FLOATS + _SHEAR_FLOATS:
        rows[name] = np.nan
    rows["clipped"] = False
    if len(pts) == 0:
        return rows
    want = np.arange(0.0, float(s[-1]) + 1e-9, max(every_um, h_um))
    take = np.unique(np.searchsorted(s, want).clip(0, len(s) - 1))
    for i in take:
        ex: dict = {} if compute_shear else None
        r = cross_section(lumen, origin_um, h_um, pts[i], tang[i],
                          float(cal[i]), span, extras=ex,
                          field=field, level=level)
        for name in _SEC_FLOATS:
            rows[name][i] = r[name]
        rows["clipped"][i] = r["clipped"]
        rows["touches_domain"][i] = r["touches_domain"]
        if compute_shear and ex:
            sh = _shear.solve_section(
                ex["mask"], ex["step_um"], ex.get("field"),
                0.5 if ex.get("level") is None else float(ex["level"]))
            for name in _SHEAR_FLOATS:
                rows[name][i] = sh[name]
            if wall_out is not None and len(sh["contour_rc"]):
                wall_out.append({
                    "i": int(i),
                    "points_zyx_um": _shear.contour_to_3d(sh["contour_rc"], ex),
                    "tau_hat_um": sh["tau_hat_um"],
                    "seg_len_um": sh["seg_len_um"],
                    "q_hat_um4": sh["q_hat_um4"]})
    return rows


def _dist_to_nodes(edge_id, s_um, edges, edge_len, targets, n_nodes):
    """Distance **along the vessel** to the nearest of a set of nodes.

    Used for distance-to-junction (the replacement for `dist_branch_um`, which
    was a Euclidean nearest-node distance in the xy projection — [P08](a),
    [P29]) and for distance-to-crop. Geodesic on the graph: Dijkstra from every
    target over node-to-node edge lengths, then each point takes the better of
    going forward or back along its own edge.

    All-NaN when there is no target node — three of these six images are not
    branched, and a straight tube genuinely has no distance-to-junction.
    Callers must guard for it; a bare ``> 20`` filter drops every row on such an
    image ([P08](d)).
    """
    import heapq
    out = np.full(len(edge_id), np.nan, np.float32)
    junc = np.asarray(targets, int)
    if len(junc) == 0 or n_nodes == 0:
        return out
    adj: dict[int, list] = {}
    for ei, (a, b) in enumerate(edges):
        if a < 0 or b < 0:
            continue
        adj.setdefault(int(a), []).append((int(b), float(edge_len[ei])))
        adj.setdefault(int(b), []).append((int(a), float(edge_len[ei])))
    D = np.full(n_nodes, np.inf)
    pq = [(0.0, int(j)) for j in junc]
    for j in junc:
        D[j] = 0.0
    heapq.heapify(pq)
    while pq:
        d, u = heapq.heappop(pq)
        if d > D[u]:
            continue
        for v, w in adj.get(u, ()):
            if d + w < D[v]:
                D[v] = d + w
                heapq.heappush(pq, (D[v], v))
    for ei, (a, b) in enumerate(edges):
        m = edge_id == ei
        if not m.any():
            continue
        s = s_um[m]
        L = float(edge_len[ei])
        da = (s + D[a]) if a >= 0 and np.isfinite(D[a]) else np.inf
        db = ((L - s) + D[b]) if b >= 0 and np.isfinite(D[b]) else np.inf
        out[m] = np.minimum(da, db)
    out[~np.isfinite(out)] = np.nan
    return out


def _empty_graph(out: dict) -> dict:
    out.update({n: np.zeros(0, np.float32)
                for n in _SEC_FLOATS + _SHEAR_FLOATS})
    out.update({
        "points_zyx_um": np.zeros((0, 3), np.float32),
        "edge_id": np.zeros(0, np.int32),
        "s_um": np.zeros(0, np.float32),
        "kappa_1_um": np.zeros(0, np.float32),
        "kappa_xy_1_um": np.zeros(0, np.float32),
        "caliber_um": np.zeros(0, np.float32),
        "section_clipped": np.zeros(0, bool),
        "section_touches_domain": np.zeros(0, bool),
        "censored": np.zeros(0, bool),
        "dist_junction_um": np.zeros(0, np.float32),
        "dist_crop_um": np.zeros(0, np.float32),
        "nodes_zyx_um": np.zeros((0, 3), np.float32),
        "node_degree": np.zeros(0, np.int32),
        "node_kind": np.zeros(0, "<U9"),
        "node_endpoint_kind": np.zeros(0, "<U9"),
        "node_dist_opening_um": np.zeros(0, np.float32),
        "node_caliber_um": np.zeros(0, np.float32),
        "edges": np.zeros((0, 2), np.int32),
        "edge_len_um": np.zeros(0, np.float32),
        "edge_censored": np.zeros(0, bool),
        "n_junctions": 0, "n_endpoints": 0, "n_crop_endpoints": 0,
        "n_terminus_endpoints": 0, "total_length_um": 0.0,
    })
    return out
