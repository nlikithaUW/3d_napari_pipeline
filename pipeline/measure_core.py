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
    centerline_notes, curvature_notes, force_notes, mesh_notes, native_notes,
    nuclei_notes, ratio_dropout_note, volume_shape_notes,
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


def compute_measurements(labels, force_vol, vessel_verts, nuc_labels_full, scale,
                         clover_native, ruby_native, raw_voxel, curvature,
                         centerline, nuc_ids=None, nuc_is_fib=None):
    """Return the per-adhesion measurement dict (or None if no adhesions).

    ``force_vol`` may be None when FRET is not configured; the force columns are
    then all-NaN and export blank.

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
    vessel_tree = KDTree(vessel_verts)
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

    # In-plane path curvature + distance to branch.
    n_ad = len(first)
    path_curv = np.full(n_ad, np.nan, dtype=np.float32)
    dist_branch = np.full(n_ad, np.nan, dtype=np.float32)
    notes += centerline_notes(centerline)
    if centerline is not None and len(centerline.get("pts_yx_um", [])):
        adh_yx = adhesion_um[:, 1:3]
        cl_tree = KDTree(np.asarray(centerline["pts_yx_um"]))
        _, cnn = cl_tree.query(adh_yx)
        s_kap = np.asarray(centerline["kappa"])[cnn][order].astype(np.float64)
        path_curv = (np.add.reduceat(s_kap, first) / counts).astype(np.float32)
        branch = np.asarray(centerline.get("branch_yx_um", []))
        if len(branch):
            bd, _ = KDTree(branch).query(adh_yx)
            dist_branch = np.minimum.reduceat(bd[order],
                                              first).astype(np.float32)

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
        "path_curv_xy": path_curv,
        "dist_branch_um": dist_branch,
        # Not a metric column — callers print it, `metric_cols` filters it out.
        "notes": notes,
    }
