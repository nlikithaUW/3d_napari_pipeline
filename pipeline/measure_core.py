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
from . import saved_io as sio
from .saved_check import (
    graph_notes, curvature_notes, force_notes, mesh_notes, native_notes,
    nuclei_notes, ratio_dropout_note, surface_verts_for_distance,
    volume_shape_notes,
)


def bs_per_z(vol: np.ndarray) -> np.ndarray:
    """Per-Z-plane background subtraction, via the frozen ``fret`` routine.

    Shared with the Measure and Adhesion-RF tabs, which each had a
    byte-identical private copy: the ratio the three produce is the tension
    proxy, so they must not be able to drift.

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
                 "lumen_aspect",
                 # Section-level shear, from the nearest *measured* station.
                 # Context for the per-wall-point columns below, not a
                 # substitute for them: tau varies 6-9x around one section, so
                 # the section mean is exactly the circumferential variation
                 # DIRECTION section 7.4 needs thrown away. `sec_force_balance`
                 # and `sec_wall_ok_frac` travel with them because they say
                 # whether that station's solve is trustworthy at all.
                 "sec_tau_rel_1_um3", "sec_tau_hat_mean_um",
                 "sec_force_balance", "sec_wall_ok_frac")

# Per-adhesion wall shear, read from the ragged `vessel_shear.npz` — one row
# per wall point, not per centerline point. CLAUDE.md fixes this as the source:
# "it is also what a per-adhesion join will read, so the surface you look at
# and
# the number an adhesion gets are the same array."
#
# `tau_rel_1_um3` = tau / (mu*Q) and is the column an unbranched vessel can be
# read from with no flow rate at all. The three that follow it are not
# decoration:
#   * `dist_tau_um`     — how far the nearest measured wall point actually was.
#     Sections are cut every few um, so most of the wall carries no
#     measurement;
#     a value from 30 um away is not about this adhesion.
#   * `tau_censored`    — the premise fails (field-of-view crop), not the
#     solve.
#   * `tau_sec_clipped` — the section was clipped: not a measurement at all.
# They mean different things and are kept as separate columns for the reason
# [P49] gives for the display layer's separate controls: one flag for all of
# them hid a third of the wall.
WALL_COLUMNS = ("tau_rel_1_um3", "tau_hat_um", "dist_tau_um",
                "tau_censored", "tau_sec_clipped")


def wall_columns(wall, adhesion_um, order, first, counts):
    """Per-adhesion wall shear from the ragged wall-point cloud.

    Each adhesion voxel takes its nearest wall point; voxels whose nearest wall
    point is further than ``wall["max_dist_um"]`` contribute **nothing**, and
    an
    adhesion with no voxel inside that radius gets NaN rather than a number
    borrowed from the far side of the vessel. That is the same rule
    `shear.vertex_field` / `covered_faces` apply when drawing the surface —
    nearest measured point, never interpolated, un-measured wall not drawn — so
    the number an adhesion gets and the surface the operator looks at agree by
    construction.

    Returns ``(columns, notes)``.
    """
    n_ad = len(first)
    out = {k: np.full(n_ad, np.nan, dtype=np.float32) for k in WALL_COLUMNS}
    if wall is None:
        return out, ["NOTE: no wall-shear cloud — tau_* columns are blank. "
                     "Build the vessel graph with shear on "
                     "(`batch_mesh --shear`, or the Vessel Mesh tab's "
                     "'Vessel wall shear' button)."]
    pts = np.asarray(wall.get("points_zyx_um", np.zeros((0, 3))), float)
    if not len(pts):
        return out, ["WARN: the wall-shear cloud is empty — tau_* columns are "
                     "blank. No cross-section produced a solvable wall."]

    max_d = float(wall.get("max_dist_um", np.inf))
    d, j = KDTree(pts).query(adhesion_um)
    near = d <= max_d
    s_near = near[order]
    s_d = d[order]

    def _reduce(vals, how="mean"):
        """Reduce a per-voxel array to per-adhesion, over in-range voxels only.

        `np.add.reduceat` has no way to skip elements, so out-of-range voxels
        are zeroed and the divisor is the in-range count. An adhesion with no
        in-range voxel divides by zero and is left NaN, which is the intended
        answer: no wall measurement reaches it.
        """
        v = np.asarray(vals, float)[order]
        v = np.where(s_near, v, 0.0)
        tot = np.add.reduceat(v, first)
        n_in = np.add.reduceat(s_near.astype(np.float64), first)
        with np.errstate(invalid="ignore", divide="ignore"):
            if how == "max":
                # A flag is set for the adhesion if any in-range voxel has it.
                return np.where(n_in > 0,
                                np.maximum.reduceat(v, first), np.nan)
            return np.where(n_in > 0, tot / n_in, np.nan)

    for col, key in (("tau_rel_1_um3", "tau_rel_1_um3"),
                     ("tau_hat_um", "tau_hat_um")):
        if key in wall:
            out[col] = _reduce(np.asarray(wall[key], float)[j]).astype(
                np.float32)
    # Distance is the *minimum* over the adhesion's voxels and is reported for
    # every adhesion, in range or not — it is how the caller judges the rest.
    out["dist_tau_um"] = np.minimum.reduceat(s_d, first).astype(np.float32)
    for col, key in (("tau_censored", "censored"),
                     ("tau_sec_clipped", "section_clipped")):
        if key in wall:
            out[col] = _reduce(
                np.asarray(wall[key], float)[j], "max").astype(np.float32)

    n_hit = int(np.isfinite(out["tau_rel_1_um3"]).sum())
    notes = [f"NOTE: wall shear reached {n_hit}/{n_ad} adhesions "
             f"({100.0 * n_hit / max(n_ad, 1):.0f}%) within {max_d:g} um of a "
             f"measured wall point; the rest are blank, not zero."]
    if n_hit:
        cens = np.nansum(out["tau_censored"] > 0)
        if cens:
            notes.append(
                f"NOTE: {int(cens)} of those {n_hit} sit on a censored "
                f"station — tau/muQ is normalised by conservation of Q, which "
                f"fails within a crop. Filter on `tau_censored`, do not "
                f"assume "
                f"it was filtered.")
    return out, notes


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
                         ("lumen_aspect", "aspect"),
                         # All-NaN unless the graph was built with
                         # compute_shear=True. Left as NaN rather than absent
                         # so the CSV header does not depend on the solve —
                         # the same contract `_SHEAR_FLOATS` states upstream.
                         ("sec_tau_rel_1_um3", "tau_rel_1_um3"),
                         ("sec_tau_hat_mean_um", "tau_hat_mean_um"),
                         ("sec_force_balance", "force_balance"),
                         ("sec_wall_ok_frac", "wall_ok_frac")):
            if key not in graph:
                continue
            v = np.asarray(graph[key])[snn][order].astype(np.float64)
            out[col] = (np.add.reduceat(v, first) / counts).astype(np.float32)
    return out, notes


def adhesion_fret(labels_n, clover_native, fret_native, ruby_native,
                  fret_params, lut_eff, lut_force, ids, clover_bs=None,
                  ruby_bs=None):
    """Per-adhesion corrected FRET efficiency and force, **sum-then-correct**.

    For each adhesion the background-subtracted donor (clover), sensitized
    emission (FRET) and acceptor (ruby) are averaged over its voxels on the
    native Z grid, and those three means go **once** through the frozen
    `fret.fret_correct` with the experiment's G and bleedthroughs, then the
    force LUT. Calling the frozen function — not a copy of its formula — keeps
    the math the reviewed one. Efficiency is a ratio, so means and sums give
    the same answer; means keep the frozen function's per-voxel checks (the
    65535 saturation test, `venus_thres` on the acceptor) meaningful.

    Why it exists beside ``force_mean_pN`` (the mean of per-voxel force): that
    averages thousands of noisy single-voxel efficiencies through a steep LUT,
    and counts every voxel the frozen math excluded as 0 pN (14% of a typical
    V04 adhesion). Here the signal is pooled first, so an adhesion gets one
    well-measured efficiency — the same "sum before dividing" reasoning as
    ``ratio_clover_ruby``.

    Voxels saturated (≥ 65535) in any of the three raw channels are left out of
    the means. Returns ``(eff, force_pN, notes)`` aligned to `ids`. An adhesion
    is NaN, not 0, when it has no usable voxels, its mean acceptor is below
    ``venus_thres`` (no measurable sensor), or its efficiency is ≤ 0 — beyond
    the sensor's range, where the LUT has no defined force.
    """
    from .fret import fret_correct

    notes: list[str] = []
    ids = np.asarray(ids)
    nmax = int(labels_n.max()) if labels_n.size else 0
    sat = ((clover_native >= 65535) | (fret_native >= 65535)
           | (ruby_native >= 65535))
    use = (labels_n > 0) & ~sat
    lv = labels_n[use]
    cbs = bs_per_z(clover_native) if clover_bs is None else clover_bs
    rbs = bs_per_z(ruby_native) if ruby_bs is None else ruby_bs
    fbs = bs_per_z(fret_native)
    n = np.bincount(lv, minlength=nmax + 1).astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        idd = np.bincount(lv, weights=cbs[use], minlength=nmax + 1) / n
        dxam = np.bincount(lv, weights=fbs[use], minlength=nmax + 1) / n
        iaa = np.bincount(lv, weights=rbs[use], minlength=nmax + 1) / n
    idd, dxam, iaa = (np.nan_to_num(a[ids]) for a in (idd, dxam, iaa))
    r = fret_correct(idd, dxam, iaa, lut_eff=lut_eff, lut_force=lut_force,
                     **fret_params)
    eff = np.asarray(r["eff"], dtype=np.float64)
    force = np.asarray(r["force"], dtype=np.float64)
    no_vox = n[ids] == 0
    no_acc = iaa < float(fret_params["venus_thres"])
    out_of_range = ~no_vox & ~no_acc & (eff <= 0)
    bad = no_vox | no_acc | out_of_range
    eff[no_vox | no_acc] = np.nan
    force[bad] = np.nan
    n_sat = int(((labels_n > 0) & sat).sum())
    if n_sat:
        notes.append(f"NOTE: {n_sat:,} saturated adhesion voxel(s) left out of "
                     f"fret_eff_adh / force_adh_pN.")
    if int(bad.sum()):
        notes.append(
            f"NOTE: fret_eff_adh / force_adh_pN blank for {int(bad.sum())}/"
            f"{ids.size} adhesions — {int(no_acc.sum())} with mean acceptor "
            f"below venus_thres, {int(out_of_range.sum())} with corrected "
            f"efficiency <= 0 (beyond the sensor's range), {int(no_vox.sum())} "
            f"with no unsaturated voxels.")
    return eff, force, notes


def compute_measurements(labels, force_vol, vessel_verts, nuc_labels_full, scale,
                         clover_native, ruby_native, raw_voxel, curvature,
                         graph, nuc_ids=None, nuc_is_fib=None,
                         vessel_faces=None, wall=None, fret_native=None,
                         fret_params=None, fret_lut=None):
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
        clover_bs = bs_per_z(clover_native)
        ruby_bs = bs_per_z(ruby_native)
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

    # Per-adhesion corrected FRET (sum-then-correct) — only when FRET is
    # configured and on for this image (fret_params given) and the native
    # channels line up with the labels.
    eff_adh = np.full(unique_ids.shape, np.nan, dtype=np.float32)
    force_adh = np.full(unique_ids.shape, np.nan, dtype=np.float32)
    if (fret_params and fret_lut is not None and fret_native is not None
            and not ratio_bad):
        e, f, fn = adhesion_fret(labels_n, clover_native, fret_native,
                                 ruby_native, fret_params, fret_lut[0],
                                 fret_lut[1], unique_ids, clover_bs, ruby_bs)
        eff_adh, force_adh = e.astype(np.float32), f.astype(np.float32)
        notes += fn

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

    wcols, wnotes = wall_columns(wall, adhesion_um, order, first, counts)
    notes += wnotes

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
        "fret_eff_adh": eff_adh,
        "force_adh_pN": force_adh,
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
        **wcols,
        # Not a metric column — callers print it, `metric_cols` filters it out.
        "notes": notes,
    }


# --------------------------------------------------------------------------- #
# The export contract
# --------------------------------------------------------------------------- #
# Shared by the Measure tab's CSV export and `batch_measure`, which used to hold
# two copies of both the column order and the settings-JSON schema. The two
# writers are still separate on purpose — the GUI writes through pandas, the
# headless path through the stdlib `csv` module, because pandas pulls in pyarrow
# whose native DLL access-violates in the Windows env. What must not differ is
# *which* columns are written, in *what* order, and what the companion JSON
# records; that is what lives here.

CSV_COLUMNS = [
    "label", "ratio_clover_ruby", "dist_mesh_um",
    "centroid_z_um", "centroid_y_um", "centroid_x_um",
    "path_curv_1_um", "path_curv_xy_1_um", "dist_junction_um",
    "dist_crop_um", "censored", "s_um", "caliber_um", "lumen_area_um2",
    "lumen_major_um", "lumen_minor_um", "lumen_aspect",
    "tau_rel_1_um3", "tau_hat_um", "dist_tau_um", "tau_censored",
    "tau_sec_clipped", "sec_tau_rel_1_um3", "sec_tau_hat_mean_um",
    "sec_force_balance", "sec_wall_ok_frac",
    "curv_mean", "curv_k1", "curv_k2", "curv_gaussian",
    "dist_nuc_um", "force_mean_pN", "force_std_pN", "force_total_pN",
    "fret_eff_adh", "force_adh_pN",
    "volume_um3", "axis_major_um", "axis_minor_um", "eccentricity", "nuc_id",
    "nuc_is_fibroblast",
]

# The two strings `nucleus_map` may take. They are NOT comparable: the unsplit
# map's ids index DAPI connected components, the split map's index nuclei (P04).
NUCLEUS_MAP_SPLIT = sio.filename(sio.NUCLEI_SPLIT)
NUCLEUS_MAP_UNSPLIT = f"{sio.filename(sio.DAPI_LABELS)} (unsplit fallback)"


def export_meta(*, image, exported, csv_name, n_adhesions, scale_um_zyx,
                roi_active, fibroblast_count, nucleus_map, notes, settings,
                headless=False) -> dict:
    """The companion `*_settings.json` beside an exported CSV.

    One schema for both export paths, so an export can be backtracked the same
    way whichever produced it. `headless` adds the flag the batch path records
    and is otherwise absent, exactly as before.
    """
    meta = {
        "image": image,
        "exported": exported,
        "csv": csv_name,
        "n_adhesions": int(n_adhesions),
        "scale_um_zyx": list(scale_um_zyx),
        "roi_active": bool(roi_active),
    }
    if headless:
        meta["headless_measure"] = True
    meta.update({
        "fibroblast_count": fibroblast_count,
        "nucleus_map": nucleus_map,
        "notes": list(notes),
        "settings": settings,
    })
    return meta
