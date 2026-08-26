"""Quality scores for a vessel-surface reconstruction — four numbers, never one.

**Why this exists.** `pipeline/batch_mesh.py` made the reconstruction knobs
explicit, recorded and cheap to vary, but nothing in the package could tell a
good mesh from a bad one, so the knobs stayed un-automatable
(`docs/ISSUES.md` [I03], "what is still entirely missing is the objective").
This module is that missing half. It scores; it does not choose.

**Four scores, reported as a vector.** They are deliberately *not* combined:

  - `scale_stability` — does the curvature stop depending on the ruler?
  - `topology`        — is it shaped like a vessel wall at all?
  - `coverage`        — is the surface everywhere the data is?
  - `support`         — is the surface anywhere the data is *not*?

Collapsing them into one number needs weights, and an invented weight is exactly
the hidden constant this project has already paid for twice — `MIN_SIZE_FLOOR =
12` ([P19]) and [P01]'s circular ``κ2 ≈ −1/R`` acceptance clause. A selection
rule is a separate, stated decision made when a sweep is run, not a default
buried here.

**Read them against each other, not alone.** Each is gameable by itself and the
others are what stop it:

- `topology` alone rewards collapsing the lumen into a solid blob, which has
  almost no handles. `coverage` measures distance to the mesh **surface**, not
  into the solid ([P09](a)), so the blob's wall voxels come out buried far from
  its boundary and it fails there.
- `coverage` alone rewards draping surface over everything, since more surface
  can only put more voxels near it. `support` is the dual — it measures the
  surface's distance to the *data* — so fabric hung across empty space fails
  there.

**`support` exists because the eye found what the other three could not** (the
2026-08-12 by-eye check on the reference/defaults pair). The reference mesh
carries spurs reaching 28.6 µm out of the vessel; `coverage` cannot see them
because they are near no wall voxel at all, and that is the wrong direction for
the question it asks. Both directions have to be asked.

**No target radius appears anywhere in this module, on purpose.** The vessel
radius varies between images and the cross-section can be oval or irregular, so
a good reconstruction is *not* one whose κ2 converges to a single number — on an
oval it genuinely differs around the circumference. `scale_stability` therefore
tests "the curvature distribution stops moving as the ruler changes", never
"it agrees with 1/R".

Pure functions: arrays in, plain dicts out. No Qt, no napari, no ``state``, no
printing — the doctrine `saved_check.py` and `vessel_mesh.build_mesh` were
written under, so the GUI and the headless path can later share one
implementation rather than growing a second copy. Callers format and print at
their own indent, which is why nothing here knows what indent it is at.
"""

from __future__ import annotations

import time

import numpy as np

from .cell_mesh import closed_edge_fraction, face_components, mesh_curvature
from .saved_check import surface_verts_for_distance

# The radii P24's plateau table was measured at, so a stability run is directly
# comparable to the evidence already on record.
DEFAULT_RADII_UM = (5.0, 10.0, 15.0, 20.0, 30.0, 40.0)

# Absolute coverage thresholds (µm). Fixed, so the number is comparable between
# two meshes; see `coverage` for why the spacing-relative ones are not.
DEFAULT_COVER_UM = (1.0, 2.0, 3.0, 5.0, 10.0)

# Coverage thresholds as multiples of the mesh's own surface-vertex spacing.
DEFAULT_COVER_REL = (0.5, 1.0, 2.0, 4.0)


# --------------------------------------------------------------------------- #
# Small shared helpers
# --------------------------------------------------------------------------- #

def _edge_use_counts(faces: np.ndarray) -> np.ndarray:
    """How many faces use each undirected edge. Same hashing as
    `cell_mesh.closed_edge_fraction`: 1 = boundary, 2 = manifold, >2 = a
    non-manifold edge where three or more triangles meet.
    """
    faces = np.asarray(faces)
    if len(faces) == 0:
        return np.zeros(0, dtype=np.int64)
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]],
                        faces[:, [2, 0]]], axis=0)
    e = np.ascontiguousarray(np.sort(e, axis=1).astype(np.int64, copy=False))
    view = e.view([("", e.dtype)] * 2).reshape(-1)
    _uniq, counts = np.unique(view, return_counts=True)
    return counts


def _nonmanifold_face_mask(faces: np.ndarray) -> np.ndarray:
    """Per-face: does this face use an edge shared by three or more faces?"""
    faces = np.asarray(faces)
    if len(faces) == 0:
        return np.zeros(0, dtype=bool)
    e = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]],
                        faces[:, [2, 0]]], axis=0)
    e = np.ascontiguousarray(np.sort(e, axis=1).astype(np.int64, copy=False))
    view = e.view([("", e.dtype)] * 2).reshape(-1)
    _uniq, inv, counts = np.unique(view, return_inverse=True,
                                   return_counts=True)
    return (counts[inv] > 2).reshape(3, len(faces)).any(axis=0)


def _f(x) -> float | None:
    """A finite float, or None — so the value survives ``json.dumps``."""
    if x is None:
        return None
    x = float(x)
    return x if np.isfinite(x) else None


def surface_vertex_spacing_um(verts: np.ndarray, faces: np.ndarray) -> float:
    """Median nearest-neighbour distance between *surface* vertices (µm).

    The mesh's own resolution limit. `docs/ISSUES.md` [P09](b) records 2.25 µm
    for the reference mesh and points out that `dist_mesh_um` is therefore a
    distance to a point cloud of that granularity, not to a surface. It is what
    makes an absolute coverage threshold interpretable: 1 µm means something
    quite different on a 5 µm lattice than on a 1 µm one.
    """
    from scipy.spatial import cKDTree

    verts = np.asarray(verts, dtype=np.float64)
    faces = np.asarray(faces)
    if len(faces) == 0 or len(verts) < 2:
        return float("nan")
    sv = verts[np.unique(faces)]
    if len(sv) < 2:
        return float("nan")
    d, _idx = cKDTree(sv).query(sv, k=2, workers=-1)
    return float(np.median(d[:, 1]))


def mesh_size(verts: np.ndarray, faces: np.ndarray) -> dict:
    """The counts every entry in `docs/ISSUES.md` keeps quoting, computed once.

    `interior_fraction` is the share of vertices belonging to no face. It is the
    [P09](a) figure: 18.8% on the hand-tuned reference mesh but 73–85% on the
    batched ones, where those buried points capture 39–75% of all nearest-vertex
    lookups. Reported free here because it costs one `np.unique` and every
    reading of `coverage` needs it.
    """
    verts = np.asarray(verts)
    faces = np.asarray(faces)
    n_v = int(len(verts))
    n_f = int(len(faces))
    n_surf = int(np.unique(faces).size) if n_f else 0
    return {
        "n_verts": n_v,
        "n_surface_verts": n_surf,
        "n_interior_verts": n_v - n_surf,
        "interior_fraction": _f((n_v - n_surf) / n_v) if n_v else None,
        "n_faces": n_f,
        "surface_vertex_spacing_um": _f(surface_vertex_spacing_um(verts, faces)),
    }


# --------------------------------------------------------------------------- #
# (a) scale stability — does the answer stop depending on the ruler?
# --------------------------------------------------------------------------- #

def _decimate_surface(verts, faces, max_verts, seed):
    """Thin a mesh to at most `max_verts` surface vertices. ``(verts, faces,
    notes)``; `verts` is returned unchanged and only `faces` is cut.

    **This is a cost control, not a free approximation, and it is off by
    default.** Keeping a random fraction *p* of surface vertices keeps only the
    faces whose three vertices all survive — about *p³* of them — so the surface
    fragments fast. Two consequences a caller must weigh:

      - The effective point spacing grows like *p*^(−1/2). At `radius_um = 5` on
        a mesh whose spacing is already 2.25 µm the quadric fit has 6–12
        neighbours ([P20]), i.e. it sits on `mesh_curvature`'s rank-deficiency
        boundary; decimating pushes it over and the small-radius end of the
        sweep degrades first.
      - `closed_edge_fraction` collapses, so `vertex_outward_normals` falls back
        to the "point away from the centroid" rule, which [P01] measured
        overriding the true winding on 44.4% of this vessel's surface vertices.

    So a decimated sweep is comparable to another decimated sweep at the same
    `max_verts` and seed, and is **not** interchangeable with a full one. The
    returned note says so, and it is carried into the result dict.
    """
    faces = np.asarray(faces)
    verts = np.asarray(verts)
    surf = np.unique(faces) if len(faces) else np.zeros(0, dtype=np.int64)
    if max_verts is None or len(surf) <= int(max_verts):
        return verts, faces, []
    rng = np.random.default_rng(int(seed))
    keep = rng.choice(surf, size=int(max_verts), replace=False)
    sel = np.zeros(len(verts), dtype=bool)
    sel[keep] = True
    kf = faces[sel[faces].all(axis=1)]
    n_kept_surf = int(np.unique(kf).size) if len(kf) else 0
    if n_kept_surf < 6 or len(kf) == 0:
        return verts, faces, [
            f"WARN: max_verts={int(max_verts)} left only {n_kept_surf} "
            f"connected surface vertices of {len(surf):,} — too few to fit. "
            f"Scored the full mesh instead."]
    frac = len(kf) / len(faces)
    return verts, kf, [
        f"NOTE: decimated to {n_kept_surf:,} of {len(surf):,} surface vertices "
        f"(seed {int(seed)}), which kept {len(kf):,} of {len(faces):,} faces "
        f"({100 * frac:.1f}%). Effective spacing grows ~"
        f"{(len(surf) / max(n_kept_surf, 1)) ** 0.5:.1f}x, so small radii "
        f"degrade first and the normals lose their closed-surface orientation. "
        f"Comparable only to another run at the same max_verts and seed."]


def scale_stability(verts, faces, radii_um=DEFAULT_RADII_UM, smooth_iter: int = 10,
                    max_verts: int | None = None, seed: int = 0,
                    normals=None) -> dict:
    """Track how the κ2 **distribution** moves as the fit window changes.

    For a surface that really is a tube of radius R, ``κ2 ≈ −1/R`` is
    independent of the fit scale while the window stays well below R. [P24]
    measured the opposite on the reference mesh: the implied radius climbs
    14 → 135 µm as `radius_um` goes 5 → 40 and never plateaus, which is the
    signature of fitting a quadric to a surface that is rough at every scale.

    Summary statistics, all computed from the per-radius medians:

    ``slope_loglog``
        Least-squares slope of ``log|median κ2|`` against ``log radius_um``.
        **0 = the number is a property of the surface; −1 = it is a pure readout
        of the ruler.** The reference mesh's recorded sweep gives ≈ −1.19 — worse
        than a pure ruler readout, because the surface is rougher than the window
        at every scale tested.

    ``slope_local``
        The same slope per *adjacent pair* of radii. A real tube produces a
        plateau (slope ≈ 0) at windows below its caliber and then a rise toward
        −1 as the window approaches it; a single overall slope averages those two
        regimes into one number and hides the shape. This is what tells a
        genuinely good reconstruction from a uniformly bad one.

    ``iqr_ratio``
        Per radius: the κ2 interquartile range over |median|. The spread is the
        part that stays meaningful on an oval or irregular cross-section, where
        the median alone is not expected to mean much — a good reconstruction of
        an oval has a *stable* spread, not a vanishing one.

    `max_verts` (default: every vertex) caps the cost — see `_decimate_surface`
    for the caveats, which are attached to the result rather than left implicit.
    `mesh_curvature` loops in Python over every surface vertex and each iteration
    does a ball query whose size grows with `radius_um`, so cost rises in both;
    the measured seconds per radius are returned so a sweep can be budgeted from
    data rather than guessed at.

    `normals` is a pass-through to `mesh_curvature`: a surface that is open by
    construction (the Poisson reconstruction, open at the x/y domain faces)
    cannot have its global orientation read back off itself, and getting that
    sign wrong flips κ2 — the very quantity this function's slope is computed
    from. `_decimate_surface` only cuts faces and leaves `verts` alone, so the
    indices stay valid through decimation.
    """
    verts = np.asarray(verts, dtype=np.float64)
    faces = np.asarray(faces)
    notes: list[str] = []
    radii = [float(r) for r in radii_um]

    sv_verts, sv_faces, dec_notes = _decimate_surface(verts, faces, max_verts,
                                                      seed)
    notes.extend(dec_notes)

    per_radius = []
    for r in radii:
        t0 = time.perf_counter()
        curv = mesh_curvature(sv_verts, sv_faces, r,
                              smooth_iter=int(smooth_iter), normals=normals)
        secs = time.perf_counter() - t0
        k1 = np.asarray(curv["k1"])
        k2 = np.asarray(curv["k2"])
        fin = np.isfinite(k2)
        n_fin = int(fin.sum())
        row = {"radius_um": r, "n_finite_k2": n_fin,
               "n_finite_k1": int(np.isfinite(k1).sum()), "seconds": round(secs, 2)}
        if n_fin:
            k = k2[fin]
            q1, med, q3 = (float(v) for v in np.percentile(k, [25, 50, 75]))
            neg = float((k < 0).mean())
            row.update({
                "k2_median": _f(med),
                "k2_q1": _f(q1), "k2_q3": _f(q3),
                "k2_iqr": _f(q3 - q1),
                "iqr_ratio": _f((q3 - q1) / abs(med)) if med != 0 else None,
                "k2_sign_consistency": _f(max(neg, 1.0 - neg)),
                "k1_median": _f(np.median(k1[np.isfinite(k1)])
                                if np.isfinite(k1).any() else np.nan),
                # Reported because [P24]'s table is in these units and the
                # comparison must be checkable — NOT because any radius is the
                # target. There is none; see the module docstring.
                "implied_R_um": _f(1.0 / abs(med)) if med != 0 else None,
            })
        else:
            row.update({"k2_median": None, "k2_q1": None, "k2_q3": None,
                        "k2_iqr": None, "iqr_ratio": None,
                        "k2_sign_consistency": None, "k1_median": None,
                        "implied_R_um": None})
            notes.append(f"WARN: no finite κ2 at radius_um={r} — the fit failed "
                         f"at every surface vertex.")
        per_radius.append(row)

    usable = [(row["radius_um"], row["k2_median"]) for row in per_radius
              if row["k2_median"] not in (None, 0.0)]
    slope = r2 = None
    if len(usable) >= 2:
        x = np.log(np.array([r for r, _ in usable], dtype=np.float64))
        y = np.log(np.abs(np.array([m for _, m in usable], dtype=np.float64)))
        slope, intercept = np.polyfit(x, y, 1)
        resid = y - (slope * x + intercept)
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = _f(1.0 - float((resid ** 2).sum()) / ss_tot) if ss_tot > 0 else None
        slope = _f(slope)
    else:
        notes.append("WARN: fewer than two radii produced a usable median κ2 — "
                     "no slope could be fitted.")

    local = []
    for (ra, ma), (rb, mb) in zip(usable[:-1], usable[1:]):
        local.append({"from_um": ra, "to_um": rb,
                      "slope": _f(np.log(abs(mb) / abs(ma))
                                  / np.log(rb / ra))})

    return {
        "radii_um": radii,
        "smooth_iter": int(smooth_iter),
        "max_verts": None if max_verts is None else int(max_verts),
        "seed": int(seed),
        "per_radius": per_radius,
        "slope_loglog": slope,
        "slope_loglog_r2": r2,
        "slope_local": local,
        "seconds": round(sum(row["seconds"] for row in per_radius), 2),
        "notes": notes,
    }


# --------------------------------------------------------------------------- #
# (b) topology — is it shaped like a vessel wall at all?
# --------------------------------------------------------------------------- #

def _component_topology(faces: np.ndarray, genus_nonmanifold_tol: float) -> dict:
    """V − E + F and the genus (when defined) for one face-connected component."""
    n_f = int(len(faces))
    v_ids = np.unique(faces)
    n_v = int(v_ids.size)
    counts = _edge_use_counts(faces)
    n_e = int(counts.size)
    n_boundary = int((counts == 1).sum())
    n_nonmanifold = int((counts > 2).sum())
    chi = n_v - n_e + n_f
    nm_frac = (n_nonmanifold / n_e) if n_e else 0.0

    out = {"n_faces": n_f, "n_verts": n_v, "n_edges": n_e,
           "n_boundary_edges": n_boundary,
           "n_nonmanifold_edges": n_nonmanifold,
           "nonmanifold_edge_fraction": _f(nm_frac),
           "closed_edge_fraction": _f(closed_edge_fraction(faces)),
           "euler_characteristic": int(chi)}

    # Genus is only defined on a closed orientable surface, where χ = 2 − 2g.
    # A component with boundary edges is a sheet, not a shell, and "genus" of a
    # mesh with debris islands or non-manifold junctions is a nonsense number —
    # so it is withheld with the reason attached rather than reported anyway.
    if n_boundary > 0:
        out["genus"] = None
        out["genus_note"] = (f"undefined: {n_boundary:,} boundary edge(s) — this "
                             f"component is an open sheet, not a closed shell")
        return out
    if nm_frac > float(genus_nonmanifold_tol):
        out["genus"] = None
        out["genus_note"] = (f"undefined: {100 * nm_frac:.1f}% of edges are "
                             f"non-manifold (>{100 * genus_nonmanifold_tol:.0f}% "
                             f"tolerance) — not a surface")
        return out
    if chi > 2:
        # Every closed orientable surface has χ ≤ 2 (the sphere), so χ = 2 − 2g
        # cannot be inverted here — it would return a *negative* genus, which is
        # not a thing. Measured: the widget-defaults mesh gives χ = +2870, i.e.
        # genus −1434. This guard is a theorem, not a tuned threshold: it fires
        # exactly when the complex is provably not a closed surface, however
        # small the non-manifold fraction happens to look.
        out["genus"] = None
        out["genus_note"] = (f"undefined: χ = {chi:+,} exceeds 2, which no "
                             f"closed orientable surface does — χ = 2 − 2g "
                             f"would give genus {(2 - chi) / 2:,.0f}. The "
                             f"complex is closed but is not a surface; read "
                             f"the non-manifold edge fraction "
                             f"({100 * nm_frac:.1f}%) instead")
        return out
    out["genus"] = _f((2 - chi) / 2.0)
    out["genus_exact"] = bool(n_nonmanifold == 0)
    if n_nonmanifold:
        out["genus_note"] = (f"nominal: closed, but {n_nonmanifold:,} of "
                             f"{n_e:,} edges ({100 * nm_frac:.1f}%) are shared "
                             f"by 3+ faces, so χ is that of a non-manifold "
                             f"complex and the genus is indicative only")
    return out


def _chi_shape(chi: int) -> str:
    """What the **sign** of χ says about how a reconstruction failed.

    This is the part `topology` used to throw away. For a closed complex,
    ``χ = Σ (2 − 2g)`` over its closed pieces, so:

      - ``χ ≪ 0`` — one shell with many tunnels bored through it. Solid, but
        perforated.
      - ``χ ≫ 0`` — many *disconnected closed pockets*: a foam. Roughly χ/2 of
        them.

    Those are two different defects, and the user's by-eye reading of the
    reference/defaults pair (2026-08-12) named exactly this distinction —
    "the reference is more solid, the defaults are more hole-y" — while the
    scorer was reporting only "genus 125" and "genus undefined". Reporting the
    sign is free and is what makes the number mean something.
    """
    if chi > 2:
        return (f"χ = {chi:+,} ≫ 0: many disconnected closed pockets — a foam "
                f"of roughly {chi // 2:,} of them around the main shell, not "
                f"one wall")
    if chi < 0:
        return (f"χ = {chi:+,} < 0: one shell with roughly {(2 - chi) // 2:,} "
                f"tunnels bored through it — solid but perforated")
    return f"χ = {chi:+,}: sphere-like (χ = 2) or torus-like (χ = 0)"


def nonmanifold_cut(faces: np.ndarray) -> dict:
    """Split the surface at its non-manifold edges and describe the pieces.

    `face_components` treats two faces as adjacent whenever they share an edge,
    including edges where three or more faces meet. That glues genuinely
    separate pockets onto the main shell through single edges, so the component
    count reads **1** for a mesh the eye plainly sees as a foam — measured on
    the widget-defaults mesh, one "component" that falls into 1,826 pieces the
    moment those edges are cut, 1,825 of them tiny.

    So this drops every face touching a non-manifold edge and re-runs the
    labelling on what is left. Reported **without a size threshold**: the piece
    count and the face share of the largest piece say it, and a "pockets are
    pieces under N faces" constant is precisely the kind of hidden number this
    project has already paid for twice.
    """
    faces = np.asarray(faces)
    if len(faces) == 0:
        return {"n_nonmanifold_faces": 0, "nonmanifold_face_fraction": None,
                "n_pieces": 0, "largest_piece_faces": 0,
                "largest_piece_face_fraction": None,
                "median_piece_faces": None}
    nm = _nonmanifold_face_mask(faces)
    keep = faces[~nm]
    if len(keep) == 0:
        return {"n_nonmanifold_faces": int(nm.sum()),
                "nonmanifold_face_fraction": _f(nm.mean()),
                "n_pieces": 0, "largest_piece_faces": 0,
                "largest_piece_face_fraction": None,
                "median_piece_faces": None}
    n_pieces, comp = face_components(keep)
    sizes = np.bincount(comp, minlength=n_pieces)
    biggest = int(sizes.max())
    return {
        "n_nonmanifold_faces": int(nm.sum()),
        "nonmanifold_face_fraction": _f(nm.mean()),
        "n_pieces": int(n_pieces),
        "largest_piece_faces": biggest,
        "largest_piece_face_fraction": _f(biggest / len(faces)),
        "median_piece_faces": int(np.median(sizes)),
    }


def topology(verts, faces, genus_nonmanifold_tol: float = 0.05,
             max_components_reported: int = 10) -> dict:
    """Component structure, Euler characteristic and genus of the surface.

    The strongest of the three scores, because the answer is known a priori
    ([I03]): a real vessel wall shell is single-digit genus — roughly 0 for a
    cone closing to a point, 1 for a single tunnel, 2 for a Y. The reference
    surface is on record at **χ = −248, genus ≈ 125**: a sponge, not a tube, and
    the kind of defect that is also visible by eye.

    Stated as *minimise genus subject to a coverage floor*, never as a target
    value — a fixed target would reintroduce exactly the hidden constant this is
    meant to remove, since the expected genus differs by image. And genus alone
    is gameable by collapsing the lumen to a blob, which is what `coverage` is
    there to catch.

    Genus is reported **only** for components closed enough for it to be
    defined, with the reason given otherwise: any boundary edge disqualifies a
    component outright, and non-manifold edges are tolerated up to
    `genus_nonmanifold_tol` of the component's edges with `genus_exact = False`
    and a note. The reference mesh needs that tolerance — [P21] measured zero
    boundary edges but only 98.4% of edges shared by exactly two faces.
    """
    verts = np.asarray(verts)
    faces = np.asarray(faces)
    size = mesh_size(verts, faces)
    notes: list[str] = []
    if len(faces) == 0:
        return {**size, "n_components": 0, "components": [], "largest": None,
                "euler_characteristic": None, "genus": None,
                "shape": "no faces", "nonmanifold_cut": nonmanifold_cut(faces),
                "closed_edge_fraction": 0.0,
                "notes": ["WARN: the mesh has no faces — there is no surface to "
                          "describe. Every geometry column read off it is "
                          "undefined."]}

    n_comp, face_comp = face_components(faces)
    order = np.argsort(np.bincount(face_comp, minlength=n_comp))[::-1]
    comps = []
    for rank, c in enumerate(order[:int(max_components_reported)]):
        comps.append({"rank": rank, **_component_topology(
            faces[face_comp == c], genus_nonmanifold_tol)})

    all_counts = _edge_use_counts(faces)
    chi_all = int(np.unique(faces).size - all_counts.size + len(faces))
    largest = comps[0] if comps else None
    cut = nonmanifold_cut(faces)
    n_small = n_comp - len(comps)
    if n_small > 0:
        notes.append(f"NOTE: {n_comp:,} face-connected components; only the "
                     f"{len(comps)} largest are described individually.")
    if n_comp > 1:
        big = int(largest["n_faces"]) if largest else 0
        notes.append(f"NOTE: the surface is in {n_comp:,} pieces; the largest "
                     f"holds {big:,} of {len(faces):,} faces "
                     f"({100 * big / len(faces):.1f}%). A vessel wall is one "
                     f"shell — the rest are islands min_component_um did not "
                     f"cut.")
    if largest and largest.get("genus") is not None and largest["genus"] >= 10:
        notes.append(f"WARN: the largest component has genus "
                     f"{largest['genus']:.0f} (χ = "
                     f"{largest['euler_characteristic']}). A vessel wall shell "
                     f"is single-digit genus; this is a sponge "
                     f"(docs/ISSUES.md I03).")
    elif largest and largest.get("genus") is None:
        # Withholding the number must not make the score go quiet: "genus is
        # undefined because this is not a surface" is a *worse* finding than a
        # large genus, and it would otherwise be the one mesh that raised
        # nothing. Lead with what the shape actually is, not with what could not
        # be computed.
        notes.append(f"WARN: {_chi_shape(chi_all)}. Genus is undefined here "
                     f"({largest.get('genus_note', 'reason not recorded')}), "
                     f"which is a stronger failure than a large genus, not a "
                     f"missing measurement (docs/ISSUES.md I03).")
    if cut["n_pieces"] > 1:
        notes.append(
            f"NOTE: cutting the {cut['n_nonmanifold_faces']:,} faces on "
            f"non-manifold edges "
            f"({100 * (cut['nonmanifold_face_fraction'] or 0):.1f}%) splits the "
            f"surface into {cut['n_pieces']:,} pieces; the largest holds "
            f"{100 * (cut['largest_piece_face_fraction'] or 0):.1f}% of the "
            f"faces and the median piece has "
            f"{cut['median_piece_faces']:,}. Face-connectivity reports one "
            f"component because single edges glue those pieces on.")

    return {
        **size,
        "n_components": int(n_comp),
        "components": comps,
        "largest": largest,
        # Over the whole face set, islands included — quoted for continuity with
        # the recorded evidence, but the per-component figures above are what
        # should be read.
        "euler_characteristic": chi_all,
        # The sign of χ is the difference between "solid but perforated" and
        # "a foam of pockets"; reporting only genus threw it away.
        "shape": _chi_shape(chi_all),
        "genus": largest["genus"] if largest else None,
        "closed_edge_fraction": _f(closed_edge_fraction(faces)),
        "n_boundary_edges": int((all_counts == 1).sum()),
        "n_nonmanifold_edges": int((all_counts > 2).sum()),
        "nonmanifold_cut": cut,
        "genus_nonmanifold_tol": float(genus_nonmanifold_tol),
        "notes": notes,
    }


# --------------------------------------------------------------------------- #
# (c) coverage — does the surface actually follow the data?
# --------------------------------------------------------------------------- #

def coverage(verts, faces, mask_points_um, cover_um=DEFAULT_COVER_UM,
             cover_rel=DEFAULT_COVER_REL, max_query: int | None = None,
             seed: int = 0) -> dict:
    """Distance from every retained cadherin voxel to the mesh **surface**.

    This is the counterweight that stops `topology` being gamed. A
    reconstruction that collapses the lumen into a solid has few handles and
    scores wonderfully on genus — but its wall voxels are then buried deep
    inside it, far from its boundary, and that shows up here. Measuring to the
    surface rather than into the solid is the whole point, so the vertex set
    comes from `saved_check.surface_verts_for_distance` — the same
    face-referenced subset `dist_mesh_um` uses since [P09](a) — including its
    warning notes, rather than a second copy of that rule.

    Two families of threshold, because they answer different questions:

    ``covered_fraction_um``
        At fixed absolute distances. **This is the one that is comparable
        between meshes**, since the ruler does not change with the thing being
        measured.

    ``covered_fraction_rel``
        At multiples of *this mesh's own* median surface-vertex spacing. A 1 µm
        criterion means something quite different on a 5 µm point lattice than
        on a 1 µm one, so this judges each reconstruction at its own resolution
        — and is correspondingly **not** comparable between two meshes whose
        spacing differs. Both are reported; a caller must say which it used.

    `mask_points_um` is (N, 3) in µm — `vessel_mesh.mask_to_points_um(mask,
    scale)`. `max_query` subsamples the queries with a fixed seed when N is
    large; the reported quantiles are robust statistics and converge long before
    the full mask is needed.
    """
    from scipy.spatial import cKDTree

    verts = np.asarray(verts, dtype=np.float64)
    faces = np.asarray(faces)
    pts = np.asarray(mask_points_um, dtype=np.float64)
    size = mesh_size(verts, faces)
    notes: list[str] = []

    sv, sv_notes = surface_verts_for_distance(verts, faces)
    notes.extend(sv_notes)
    sv = np.asarray(sv, dtype=np.float64)
    if len(sv) == 0 or len(pts) == 0:
        return {**size, "n_mask_points": int(len(pts)), "n_queried": 0,
                "distance_um": {}, "covered_fraction_um": {},
                "covered_fraction_rel": {},
                "notes": notes + ["WARN: nothing to measure — the mesh has no "
                                  "surface vertices or the mask is empty."]}

    n_total = int(len(pts))
    if max_query is not None and n_total > int(max_query):
        rng = np.random.default_rng(int(seed))
        pts = pts[rng.choice(n_total, size=int(max_query), replace=False)]
        notes.append(f"NOTE: queried a random {len(pts):,} of {n_total:,} mask "
                     f"voxels (seed {int(seed)}).")

    t0 = time.perf_counter()
    d, _idx = cKDTree(sv).query(pts, k=1, workers=-1)
    secs = time.perf_counter() - t0

    q = np.percentile(d, [50, 75, 90, 95])
    spacing = size["surface_vertex_spacing_um"]
    cov_um = {f"{float(t):g}": _f((d <= float(t)).mean()) for t in cover_um}
    cov_rel = {}
    if spacing:
        for k in cover_rel:
            cov_rel[f"{float(k):g}x"] = {
                "d_um": _f(float(k) * spacing),
                "fraction": _f((d <= float(k) * spacing).mean())}
    else:
        notes.append("WARN: surface-vertex spacing is undefined, so the "
                     "spacing-relative coverage thresholds were not computed.")

    worst = _f(cov_um.get(f"{float(DEFAULT_COVER_UM[1]):g}"))
    if worst is not None and worst < 0.5:
        notes.append(f"WARN: only {100 * worst:.1f}% of retained cadherin "
                     f"voxels lie within {DEFAULT_COVER_UM[1]:g} µm of the "
                     f"surface — most of the signal the mesh is built from is "
                     f"not near it.")

    return {
        **size,
        "n_mask_points": n_total,
        "n_queried": int(len(pts)),
        "distance_um": {"median": _f(q[0]), "p75": _f(q[1]), "p90": _f(q[2]),
                        "p95": _f(q[3]), "max": _f(d.max()),
                        "mean": _f(d.mean())},
        "covered_fraction_um": cov_um,
        "covered_fraction_rel": cov_rel,
        "seconds": round(secs, 2),
        "notes": notes,
    }


# --------------------------------------------------------------------------- #
# (d) support — is the surface anywhere the data is not?
# --------------------------------------------------------------------------- #

def support(verts, faces, mask_points_um, reach_um=DEFAULT_COVER_UM,
            reach_rel=DEFAULT_COVER_REL, seed: int = 0) -> dict:
    """Distance from the mesh's **triangle interiors** to the nearest cadherin
    voxel — the dual of `coverage`, and the axis the eye caught that the other
    three scores are blind to.

    **Why this is a fourth score and not part of `coverage`.** `coverage` asks
    "is every wall voxel near the mesh?" A spike — surface reaching out into
    empty space — is near no wall voxel at all, so it is invisible in that
    direction. Both questions have to be asked, and asking only the first is
    what let a reconstruction with 28.6 µm spurs score well.

    **Why the interiors and not the vertices.** For the default alpha-3D method
    the mesh vertices *are* cadherin voxels (a subsample of them), so
    vertex-to-data distance is identically zero on every mesh this project
    builds — a vacuous test. The fabricated part of an alpha shape is the
    *surface between* real points: a triangle bridging from the wall out to an
    isolated debris voxel has three perfectly good corners and a middle that is
    nowhere near anything. So each face is sampled at its centroid and its three
    edge midpoints, and every statistic is **area-weighted**, because a spike's
    significance is the area it adds, not the number of triangles.

    Reported at absolute distances (comparable between meshes) and at multiples
    of the mesh's own surface-vertex spacing (fair to each mesh's resolution) —
    the same dual as `coverage`, and for the same reason. The absolute figure is
    the one that matters physically: a spur sticking out of a vessel does not
    care what the mesh's point spacing is.

    **No threshold warning, deliberately.** Some surface spanning empty space is
    the method working: the cadherin signal is genuinely broken into arcs and
    bridging them is the point of `alpha3d_um`. Separating "legitimate bridge"
    from "spike" needs a distance constant nobody has justified, so this reports
    the distribution and the maximum reach as facts and leaves the cut to
    whoever states a selection rule. Measured 2026-08-12: aspect ratio does
    *not* separate them either (it ranks the two meshes the opposite way);
    absolute reach does.
    """
    from scipy.spatial import cKDTree

    verts = np.asarray(verts, dtype=np.float64)
    faces = np.asarray(faces)
    pts = np.asarray(mask_points_um, dtype=np.float64)
    size = mesh_size(verts, faces)
    notes: list[str] = []

    if len(faces) == 0 or len(pts) == 0:
        return {**size, "area_um2": 0.0, "distance_um": {},
                "unsupported_area_um": {}, "unsupported_area_rel": {},
                "notes": notes + ["WARN: nothing to measure — the mesh has no "
                                  "faces or the mask is empty."]}

    p = verts[faces]                                    # (F, 3, 3)
    area = 0.5 * np.linalg.norm(
        np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    total = float(area.sum())

    # Centroid + the three edge midpoints. A fixed, deterministic pattern: the
    # score has to be reproducible run to run, and random barycentric sampling
    # would make it wobble by more than the differences being measured.
    samples = np.concatenate([
        p.mean(axis=1),
        0.5 * (p[:, 0] + p[:, 1]),
        0.5 * (p[:, 1] + p[:, 2]),
        0.5 * (p[:, 2] + p[:, 0]),
    ], axis=0)
    w = np.tile(area / 4.0, 4)

    t0 = time.perf_counter()
    d, _idx = cKDTree(pts).query(samples, k=1, workers=-1)
    secs = time.perf_counter() - t0

    order = np.argsort(d)
    d_s, w_s = d[order], w[order]
    cw = np.cumsum(w_s)
    cw = cw / cw[-1] if cw[-1] > 0 else cw

    def wq(q):
        """Area-weighted quantile of the sample distances."""
        return _f(float(d_s[np.searchsorted(cw, q)])) if cw[-1] > 0 else None

    # Edge lengths: the sharpest single discriminator measured on the
    # reference/defaults pair (max 28.6 µm against 19.9 µm), and free.
    e = np.concatenate([
        np.linalg.norm(p[:, 1] - p[:, 0], axis=1),
        np.linalg.norm(p[:, 2] - p[:, 1], axis=1),
        np.linalg.norm(p[:, 0] - p[:, 2], axis=1)])

    spacing = size["surface_vertex_spacing_um"]
    # Area-weighted over all four samples rather than the centroid alone, so a
    # triangle that is fine at its middle but reaches far at one edge counts.
    unsup_um = {}
    for t in reach_um:
        far = d > float(t)
        unsup_um[f"{float(t):g}"] = _f(w[far].sum() / total) if total > 0 else None
    unsup_rel = {}
    if spacing:
        for k in reach_rel:
            far = d > float(k) * spacing
            unsup_rel[f"{float(k):g}x"] = {
                "d_um": _f(float(k) * spacing),
                "fraction": _f(w[far].sum() / total) if total > 0 else None}
    else:
        notes.append("WARN: surface-vertex spacing is undefined, so the "
                     "spacing-relative reach thresholds were not computed.")

    notes.append(
        f"NOTE: the surface reaches at most {d.max():.2f} µm from any retained "
        f"cadherin voxel, and its longest edge is {e.max():.2f} µm. Some of "
        f"that is intended gap-bridging, so no threshold is applied here — read "
        f"it against alpha3d_um, which bounds it by construction.")

    return {
        **size,
        "area_um2": _f(total),
        "n_mask_points": int(len(pts)),
        "n_samples": int(len(samples)),
        "distance_um": {"median": wq(0.50), "p75": wq(0.75), "p90": wq(0.90),
                        "p95": wq(0.95), "p99": wq(0.99), "max": _f(d.max())},
        "unsupported_area_um": unsup_um,
        "unsupported_area_rel": unsup_rel,
        "edge_length_um": {"median": _f(np.median(e)),
                           "p99": _f(np.percentile(e, 99)), "max": _f(e.max())},
        "seconds": round(secs, 2),
        "notes": notes,
    }


# --------------------------------------------------------------------------- #
# The four together
# --------------------------------------------------------------------------- #

def score_mesh(verts, faces, mask_points_um=None, radii_um=DEFAULT_RADII_UM,
               smooth_iter: int = 10, max_verts: int | None = None,
               max_query: int | None = None, seed: int = 0,
               normals=None) -> dict:
    """All three scores as one dict — **a vector, not a total.**

    There is deliberately no combined figure of merit here. Weighting stability
    against topology against coverage requires weights nobody has justified, and
    the selection rule belongs to whoever runs a sweep, stated explicitly at that
    point. What this returns is the evidence to state it over.

    `mask_points_um` may be None, in which case coverage is skipped and said to
    have been — scoring a mesh without its mask is a legitimate half-answer
    (curvature and topology need no pixels), not a silent one.

    `normals` is a pass-through to `scale_stability`; see there for why an
    open-by-construction surface must supply its own.
    """
    verts = np.asarray(verts)
    faces = np.asarray(faces)
    t0 = time.perf_counter()
    out = {
        "scale_stability": scale_stability(verts, faces, radii_um=radii_um,
                                           smooth_iter=smooth_iter,
                                           max_verts=max_verts, seed=seed,
                                           normals=normals),
        "topology": topology(verts, faces),
    }
    if mask_points_um is None:
        out["coverage"] = {**mesh_size(verts, faces), "notes": [
            "WARN: no cadherin mask supplied, so coverage was not computed — "
            "topology alone is gameable by collapsing the lumen to a blob."]}
        out["support"] = {**mesh_size(verts, faces), "notes": [
            "WARN: no cadherin mask supplied, so support was not computed — "
            "surface reaching into empty space goes unmeasured."]}
    else:
        out["coverage"] = coverage(verts, faces, mask_points_um,
                                   max_query=max_query, seed=seed)
        out["support"] = support(verts, faces, mask_points_um, seed=seed)
    out["seconds"] = round(time.perf_counter() - t0, 2)
    out["notes"] = [n for part in ("scale_stability", "topology", "coverage",
                                   "support")
                    for n in out[part].get("notes", [])]
    return out


def format_scores(scores: dict) -> list[str]:
    """The three scores as report lines. Pure: dict in, ``list[str]`` out.

    Callers print at their own indent (`saved_check`'s convention), which is why
    nothing here knows what indent it is at.
    """
    lines: list[str] = []
    st = scores.get("scale_stability") or {}
    tp = scores.get("topology") or {}
    cv = scores.get("coverage") or {}

    lines.append(f"mesh: {tp.get('n_verts', 0):,} verts "
                 f"({tp.get('n_surface_verts', 0):,} surface, "
                 f"{100 * (tp.get('interior_fraction') or 0):.1f}% interior), "
                 f"{tp.get('n_faces', 0):,} faces, spacing "
                 f"{tp.get('surface_vertex_spacing_um') or float('nan'):.2f} um")

    slope = st.get("slope_loglog")
    lines.append(f"scale_stability: slope_loglog = "
                 f"{'n/a' if slope is None else f'{slope:+.3f}'}"
                 f" (0 = a property of the surface, -1 = a pure readout of the "
                 f"ruler)")
    for row in st.get("per_radius", []):
        head = f"    r={row['radius_um']:>5g} um  "
        med = row.get("k2_median")
        if med is None:
            lines.append(f"{head}k2 all non-finite  "
                         f"[{row.get('seconds', 0):.1f} s]")
            continue
        lines.append(
            f"{head}k2 med {med:+.4f}  IQR/|med| "
            f"{row.get('iqr_ratio') or float('nan'):5.2f}  sign "
            f"{100 * (row.get('k2_sign_consistency') or 0):5.1f}%  "
            f"n={row.get('n_finite_k2', 0):,}  [{row.get('seconds', 0):.1f} s]")
    for step in st.get("slope_local", []):
        head = f"    local slope {step['from_um']:g}->{step['to_um']:g} um: "
        s = step.get("slope")
        lines.append(head + ("n/a" if s is None else f"{s:+.3f}"))

    g = tp.get("genus")
    cut = tp.get("nonmanifold_cut") or {}
    lines.append(f"topology: {tp.get('shape', '?')}; largest genus "
                 f"{'undefined' if g is None else f'{g:.0f}'}, closed edges "
                 f"{100 * (tp.get('closed_edge_fraction') or 0):.1f}%")
    if cut.get("n_pieces"):
        lines.append(f"    cut at the {cut['n_nonmanifold_faces']:,} "
                     f"non-manifold faces -> {cut['n_pieces']:,} piece(s), "
                     f"largest holds "
                     f"{100 * (cut.get('largest_piece_face_fraction') or 0):.1f}%"
                     f" of the faces, median piece "
                     f"{cut.get('median_piece_faces')} face(s)")

    dist = cv.get("distance_um") or {}
    if dist:
        lines.append(f"coverage: {cv.get('n_queried', 0):,} mask voxels, "
                     f"distance to surface median "
                     f"{dist.get('median') or float('nan'):.2f} / p90 "
                     f"{dist.get('p90') or float('nan'):.2f} / max "
                     f"{dist.get('max') or float('nan'):.2f} um")
        frac = "  ".join(
            f"<={t} um {100 * (v or 0):.1f}%"
            for t, v in (cv.get("covered_fraction_um") or {}).items())
        lines.append(f"    absolute (comparable between meshes): {frac}")
        rel = "  ".join(
            f"{k} (={v['d_um']:.2f} um) {100 * (v['fraction'] or 0):.1f}%"
            for k, v in (cv.get("covered_fraction_rel") or {}).items())
        if rel:
            lines.append(f"    per-mesh spacing (not comparable): {rel}")

    su = scores.get("support") or {}
    sd = su.get("distance_um") or {}
    if sd:
        el = su.get("edge_length_um") or {}
        lines.append(f"support: {su.get('area_um2') or float('nan'):,.0f} um2 of "
                     f"surface, distance to data median "
                     f"{sd.get('median') or float('nan'):.2f} / p95 "
                     f"{sd.get('p95') or float('nan'):.2f} / max "
                     f"{sd.get('max') or float('nan'):.2f} um  (area-weighted)")
        frac = "  ".join(
            f">{t} um {100 * (v or 0):.1f}%"
            for t, v in (su.get("unsupported_area_um") or {}).items())
        lines.append(f"    unsupported area, absolute: {frac}")
        lines.append(f"    edge length: median "
                     f"{el.get('median') or float('nan'):.2f} / p99 "
                     f"{el.get('p99') or float('nan'):.2f} / max "
                     f"{el.get('max') or float('nan'):.2f} um")
    return lines
