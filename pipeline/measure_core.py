"""Standalone per-adhesion measurement (no napari / no Qt).

Mirrors the Measure tab's worker exactly so the headless batch produces the same
CSV columns. Given the adhesion label volume + the saved intermediates (force
volume, vessel vertices, DAPI labels, native clover/ruby channels, surface
curvature, centerline), it returns the same result dict the GUI assembles.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.spatial import cKDTree as KDTree


def compute_measurements(labels, force_vol, vessel_verts, nuc_labels_full, scale,
                         clover_native, ruby_native, raw_voxel, curvature,
                         centerline):
    """Return the per-adhesion measurement dict (or None if no adhesions)."""
    from skimage.measure import regionprops_table
    import scipy.ndimage as sndi

    labels = np.asarray(labels).astype(np.int32)
    unique_ids = np.unique(labels[labels > 0])
    if unique_ids.size == 0:
        return None

    props = regionprops_table(
        labels, spacing=scale,
        properties=["label", "area", "axis_major_length",
                    "axis_minor_length", "centroid"],
    )

    force_mean = np.array(sndi.mean(force_vol, labels, index=unique_ids))
    force_std = np.array(sndi.standard_deviation(force_vol, labels,
                                                 index=unique_ids))

    # Per-adhesion clover/ruby ratio on native Z.
    if (clover_native is not None and ruby_native is not None
            and clover_native.shape[1:] == labels.shape[1:]):
        Zn = clover_native.shape[0]
        zf = float(raw_voxel[0]) / float(scale[0])
        iso_idx = np.clip(np.rint(np.arange(Zn) * zf).astype(np.int64),
                          0, labels.shape[0] - 1)
        labels_n = labels[iso_idx]
        ruby_f = ruby_native.astype(np.float32)
        inside = (labels_n > 0) & (ruby_f > 0)
        rv = clover_native.astype(np.float32)[inside] / ruby_f[inside]
        lv = labels_n[inside]
        nmax = int(labels.max())
        sums = np.bincount(lv, weights=rv, minlength=nmax + 1)
        cnts = np.bincount(lv, minlength=nmax + 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio_by_label = sums / cnts
        ratio_mean = ratio_by_label[unique_ids].astype(np.float32)
        ratio_mean[cnts[unique_ids] == 0] = np.nan
    else:
        ratio_mean = np.full(unique_ids.shape, np.nan, dtype=np.float32)

    adhesion_vox = np.argwhere(labels > 0)
    adhesion_ids = labels[adhesion_vox[:, 0], adhesion_vox[:, 1],
                          adhesion_vox[:, 2]]
    adhesion_um = adhesion_vox * np.array(scale)

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

    # Surface curvature at nearest surface vertex.
    nverts = len(vessel_verts)
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
        "ratio_clover_ruby": ratio_mean,
        "curv_mean": curv_out["mean"],
        "curv_k1": curv_out["k1"],
        "curv_k2": curv_out["k2"],
        "curv_gaussian": curv_out["gaussian"],
        "path_curv_xy": path_curv,
        "dist_branch_um": dist_branch,
    }
