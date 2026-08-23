"""Standalone per-adhesion measurement (no napari / no Qt).

Produces the same result dict — and therefore the same CSV columns — as the
Measure tab's worker, given the adhesion label volume + the saved intermediates
(force volume, vessel vertices, DAPI labels, native clover/ruby channels,
surface curvature, centerline).

NOT a byte-for-byte mirror of the GUI worker: this copy has no X-Y ROI filter,
so it measures every adhesion in the volume where the GUI would measure only
those inside a drawn region. The two implementations are duplicated and have
already diverged once — see docs/ISSUES.md P11 ("measure_core duplicates the
Measure worker and has diverged").
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.spatial import cKDTree as KDTree

from .fret import background_subtract
from .saved_check import (
    graph_notes, curvature_notes, force_notes, mesh_notes, native_notes,
    nuclei_notes, ratio_dropout_note, surface_verts_for_distance,
    volume_shape_notes,
)


def _bs_per_z(vol: np.ndarray) -> np.ndarray:
    """Per-Z-plane background subtraction, via the frozen ``fret`` routine.

    Each plane's modal intensity is its own background (the mode shifts with
    depth as the illumination attenuates), so this must be applied plane by
    plane rather than once over the volume.
    """
    out = np.empty(vol.shape, dtype=np.float32)
    for z in range(vol.shape[0]):
        # background_subtract returns float64; write straight into the float32
        # buffer so only one plane's worth of float64 exists at a time (np.stack
        # of every plane would peak at 2x the volume in float64).
        out[z] = background_subtract(vol[z])
    return out


GRAPH_COLUMNS = ("path_curv_1_um", "path_curv_xy_1_um", "dist_junction_um",
                 "dist_crop_um", "censored", "s_um", "caliber_um",
                 "lumen_area_um2", "lumen_major_um", "lumen_minor_um",
                 "lumen_aspect")


def graph_columns(graph, adhesion_um, order, first, counts):
    """Per-adhesion columns read off the 3D centerline graph.

    Extracted so the GUI worker and `measure_core` share one implementation
    rather than a fourth copy — [P11] records that [P04] and [P12] each had to
    be written twice, and this is the change where that stops getting worse.
    Returns ``(columns, notes)``.

    Each adhesion attaches to its **nearest point on the graph in 3D**, and
    inherits that point's whole record. The old code did this in the xy
    projection against a raster-order point set, so an adhesion on the far wall
    and one on the near wall at the same (y,x) received the same curvature, and
    "distance to branch" was Euclidean rather than along the vessel ([P08](a),
    [P27], [P29]).

    ``caliber_um`` and the three ``lumen_*`` columns are new, and they are the
    ones the shear estimate needs: for a flat duct tau ~ 6 uQ/(w h^2), so width
    and height enter with different exponents and must travel as separate
    columns rather than be collapsed into one caliber (DIRECTION §4).
    """
    n_ad = len(first)
    out = {k: np.full(n_ad, np.nan, dtype=np.float32) for k in GRAPH_COLUMNS}
    out["censored"] = np.zeros(n_ad, dtype=np.float32)
    notes = list(graph_notes(graph))
    gpts = (np.asarray(graph.get("points_zyx_um", []))
            if graph is not None else np.zeros((0, 3)))
    if not len(gpts):
        return out, notes

    _, cnn = KDTree(gpts).query(adhesion_um)

    def reduce_to_adhesion(arr, how="mean"):
        a = np.asarray(arr)
        if a.shape[0] != len(gpts):
            return np.full(n_ad, np.nan, dtype=np.float32)
        v = a[cnn][order].astype(np.float64)
        if how == "min":
            return np.minimum.reduceat(v, first).astype(np.float32)
        if how == "max":
            return np.maximum.reduceat(v, first).astype(np.float32)
        return (np.add.reduceat(v, first) / counts).astype(np.float32)

    out["path_curv_1_um"] = reduce_to_adhesion(graph["kappa_1_um"])
    out["path_curv_xy_1_um"] = reduce_to_adhesion(graph["kappa_xy_1_um"])
    out["caliber_um"] = reduce_to_adhesion(graph["caliber_um"])
    out["s_um"] = reduce_to_adhesion(graph["s_um"])
    if "dist_junction_um" in graph:
        out["dist_junction_um"] = reduce_to_adhesion(
            graph["dist_junction_um"], "min")
    if "dist_crop_um" in graph:
        out["dist_crop_um"] = reduce_to_adhesion(graph["dist_crop_um"], "min")
    if "censored" in graph:
        out["censored"] = reduce_to_adhesion(
            np.asarray(graph["censored"], np.float64), "max")

    # The cross-sections are cut every few µm and are NaN in between, and that
    # gap is deliberate — an interpolated section is not a measurement ([P26]).
    # So each adhesion takes the nearest *measured* section rather than
    # averaging into NaN.
    sec_ok = np.isfinite(np.asarray(graph.get("area_um2", [])))
    if sec_ok.any():
        sidx = np.nonzero(sec_ok)[0]
        _, snn = KDTree(gpts[sidx]).query(adhesion_um)
        snn = sidx[snn]
        for col, key in (("lumen_area_um2", "area_um2"),
                         ("lumen_major_um", "major_um"),
                         ("lumen_minor_um", "minor_um"),
                         ("lumen_aspect", "aspect")):
            v = np.asarray(graph[key])[snn][order].astype(np.float64)
            out[col] = (np.add.reduceat(v, first) / counts).astype(np.float32)
    return out, notes


def compute_measurements(labels, force_vol, vessel_verts, nuc_labels_full, scale,
                         clover_native, ruby_native, raw_voxel, curvature,
                         graph, nuc_ids=None, nuc_is_fib=None,
                         vessel_faces=None):
    """Return the per-adhesion measurement dict (or None if no adhesions).

    ``force_vol`` may be None when FRET is not configured; the force columns are
    then all-NaN and export blank.

    ``vessel_faces`` are the mesh faces. Given them, ``dist_mesh_um`` is
    measured to the mesh *surface* (vertices carrying at least one face);
    without them it falls back to every vertex, which on an alpha-shape solid
    includes points buried inside the wall and reads short (`docs/ISSUES.md`
    P09(a)). The fallback is the pre-2026-08-12 behaviour and warns, so an old
    caller degrades loudly rather than changing answers in silence.

    ``nuc_ids`` / ``nuc_is_fib`` are the per-nucleus label ids and their
    fibroblast classification (``nuclei_split.npz``). Given both, ``nuc_id`` can
    be joined to cell type and the ``nuc_is_fibroblast`` column is populated;
    given neither it comes out blank. The raw vectors are passed rather than a
    prebuilt lookup table so both callers get the same join and the same
    validation, instead of duplicating it (`docs/ISSUES.md` P11).

    **This function never prints.** It does not know its caller's indent level
    and is the one module a non-CLI caller could reuse, so every degradation it
    detects is accumulated into ``result["notes"]`` for the caller to report.
    ``notes`` is not a metric column — both writers select columns by name, so
    it never reaches the CSV.
    """
    from skimage.measure import regionprops_table
    import scipy.ndimage as sndi

    notes: list[str] = []

    labels = np.asarray(labels).astype(np.int32)
    unique_ids = np.unique(labels[labels > 0])
    if unique_ids.size == 0:
        return None

    props = regionprops_table(
        labels, spacing=scale,
        properties=["label", "area", "axis_major_length",
                    "axis_minor_length", "centroid"],
    )

    notes += force_notes(force_vol)
    if force_vol is None:
        force_mean = np.full(unique_ids.shape, np.nan, dtype=np.float64)
        force_std = np.full(unique_ids.shape, np.nan, dtype=np.float64)
    else:
        force_mean = np.array(sndi.mean(force_vol, labels, index=unique_ids))
        force_std = np.array(sndi.standard_deviation(force_vol, labels,
                                                     index=unique_ids))

    # Per-adhesion clover/ruby ratio on native Z: total clover / total ruby over
    # the adhesion, both background-subtracted per Z plane first. Summing before
    # dividing keeps one dim voxel from dominating (a mean of per-voxel ratios
    # is driven by its smallest denominators), and the background subtraction
    # stops dim/deep adhesions from drifting toward a fixed offset ratio.
    ratio_bad = native_notes(clover_native, ruby_native, labels.shape)
    notes += ratio_bad
    if not ratio_bad:
        Zn = clover_native.shape[0]
        zf = float(raw_voxel[0]) / float(scale[0])
        iso_idx = np.clip(np.rint(np.arange(Zn) * zf).astype(np.int64),
                          0, labels.shape[0] - 1)
        labels_n = labels[iso_idx]
        clover_bs = _bs_per_z(clover_native)
        ruby_bs = _bs_per_z(ruby_native)
        inside = labels_n > 0
        lv = labels_n[inside]
        nmax = int(labels.max())
        c_sum = np.bincount(lv, weights=clover_bs[inside], minlength=nmax + 1)
        r_sum = np.bincount(lv, weights=ruby_bs[inside], minlength=nmax + 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio_by_label = c_sum / r_sum
        # An adhesion whose ruby is entirely at/below background has no ratio.
        ratio_by_label[r_sum <= 0] = np.nan
        ratio_mean = ratio_by_label[unique_ids].astype(np.float32)
        notes += ratio_dropout_note(int(np.isnan(ratio_mean).sum()),
                                    int(ratio_mean.size))
    else:
        ratio_mean = np.full(unique_ids.shape, np.nan, dtype=np.float32)

    adhesion_vox = np.argwhere(labels > 0)
    adhesion_ids = labels[adhesion_vox[:, 0], adhesion_vox[:, 1],
                          adhesion_vox[:, 2]]
    adhesion_um = adhesion_vox * np.array(scale)

    notes += mesh_notes(vessel_verts)
    # A nucleus volume from a different image would not crash — it would just
    # produce wrong distances — so this has to be said rather than caught.
    notes += volume_shape_notes("the nucleus label volume",
                                np.shape(nuc_labels_full),
                                "the adhesion label volume", labels.shape)
    # dist_mesh_um measures to the mesh surface, not to interior points of the
    # alpha-shape solid (P09a). `mesh_nn_idx` below is unused, so restricting
    # the tree changes distances only — nothing is looked up by vertex index.
    dist_verts, dist_notes = surface_verts_for_distance(vessel_verts,
                                                        vessel_faces)
    notes += dist_notes
    vessel_tree = KDTree(dist_verts)
    nuc_sub = nuc_labels_full[::3, ::3, ::3]
    nuc_vox = np.argwhere(nuc_sub > 0)
    nuc_ids_sub = nuc_sub[nuc_vox[:, 0], nuc_vox[:, 1], nuc_vox[:, 2]]
    nuc_um = nuc_vox * 3 * np.array(scale)
    nuc_tree = KDTree(nuc_um)

    with ThreadPoolExecutor(max_workers=2) as ex:
        f_mesh = ex.submit(vessel_tree.query, adhesion_um)
        f_nuc = ex.submit(nuc_tree.query, adhesion_um)
        mesh_dist_all, mesh_nn_idx = f_mesh.result()
        nuc_dist_all, nuc_nn_idx = f_nuc.result()
    nuc_ids_at_nn = nuc_ids_sub[nuc_nn_idx]

    order = np.argsort(adhesion_ids, kind="stable")
    s_ids = adhesion_ids[order]
    s_mesh = mesh_dist_all[order]
    s_nuc = nuc_dist_all[order]
    s_nuc_id = nuc_ids_at_nn[order]

    _, first = np.unique(s_ids, return_index=True)
    counts = np.diff(np.concatenate([first, [len(s_ids)]]))

    mesh_min = np.minimum.reduceat(s_mesh, first)
    nuc_min = np.minimum.reduceat(s_nuc, first)

    group_argmin = np.array([
        first[i] + np.argmin(s_nuc[first[i]: first[i] + counts[i]])
        for i in range(len(first))
    ])
    nuc_assignment = s_nuc_id[group_argmin]

    # Join each adhesion's assigned nucleus to that nucleus's fibroblast class
    # (P04). float32 + NaN rather than an int sentinel: both writers blank NaN
    # and only NaN, so a -1 would export literally and be read as data.
    nuc_is_fibroblast = np.full(len(first), np.nan, dtype=np.float32)
    if nuc_ids is not None or nuc_is_fib is not None:
        nuc_max = int(np.asarray(nuc_labels_full).max())
        nn = nuclei_notes(nuc_ids, nuc_is_fib, nuc_max)
        notes += nn
        if not nn:
            lut = np.full(nuc_max + 1, np.nan, dtype=np.float32)
            lut[np.asarray(nuc_ids).astype(np.int64)] = np.asarray(
                nuc_is_fib, bool).astype(np.float32)
            nuc_is_fibroblast = lut[np.asarray(nuc_assignment).astype(np.int64)]
    else:
        notes += nuclei_notes(None, None, 0)

    # Surface curvature at nearest surface vertex.
    nverts = len(vessel_verts)
    notes += curvature_notes(curvature, nverts)
    curv_out = {}
    k1_all = None if curvature is None else curvature.get("k1")
    k1_all = None if k1_all is None else np.asarray(k1_all)
    surf_ok = (k1_all is not None and k1_all.shape[0] == nverts
               and bool(np.isfinite(k1_all).any()))
    if surf_ok:
        surf_mask = np.isfinite(k1_all)
        surf_tree = KDTree(vessel_verts[surf_mask])
        _, surf_nn = surf_tree.query(adhesion_um)
        s_surf_nn = surf_nn[order]
    for key in ("mean", "k1", "k2", "gaussian"):
        cv = None if curvature is None else curvature.get(key)
        cv = None if cv is None else np.asarray(cv)
        if not surf_ok or cv is None or cv.shape[0] != nverts:
            curv_out[key] = np.full(len(first), np.nan, dtype=np.float32)
            continue
        s_cv = cv[surf_mask][s_surf_nn].astype(np.float64)
        curv_out[key] = (np.add.reduceat(s_cv, first) / counts).astype(np.float32)

    gcols, gnotes = graph_columns(graph, adhesion_um, order, first, counts)
    notes += gnotes

    return {
        "label": unique_ids,
        "volume_um3": props["area"],
        "axis_major_um": props["axis_major_length"],
        "axis_minor_um": props["axis_minor_length"],
        "eccentricity": (props["axis_major_length"]
                         / np.maximum(props["axis_minor_length"], 1e-9)),
        "force_mean_pN": force_mean,
        "force_std_pN": force_std,
        "force_total_pN": force_mean * props["area"],
        "dist_mesh_um": mesh_min,
        "dist_nuc_um": nuc_min,
        "nuc_id": nuc_assignment,
        # Unconditional, so the CSV header is stable whether or not the join
        # was available. 1.0 = VnTs+ fibroblast, 0.0 = other, blank = unknown.
        "nuc_is_fibroblast": nuc_is_fibroblast,
        "ratio_clover_ruby": ratio_mean,
        # regionprops_table was given `spacing=scale`, so these are µm.
        "centroid_z_um": props["centroid-0"],
        "centroid_y_um": props["centroid-1"],
        "centroid_x_um": props["centroid-2"],
        "curv_mean": curv_out["mean"],
        "curv_k1": curv_out["k1"],
        "curv_k2": curv_out["k2"],
        "curv_gaussian": curv_out["gaussian"],
        # **Renamed, because they are different quantities.** The old
        # `path_curv_xy` was fitted to a Euclidean neighbourhood of an
        # unordered xy point set; these are fitted along arc length on an
        # ordered 3D path. Reusing the name would repeat [P04]'s `nuc_id`
        # trap — the same header meaning a different thing across a boundary.
        **gcols,
        # Not a metric column — callers print it, `metric_cols` filters it out.
        "notes": notes,
    }
