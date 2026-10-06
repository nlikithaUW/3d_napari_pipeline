"""Measure tab — per-adhesion force, shape, and distance metrics."""

from __future__ import annotations

import os

import numpy as np
from concurrent.futures import ThreadPoolExecutor
from magicgui.widgets import Container, FileEdit, PushButton, SpinBox
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from napari.qt.threading import thread_worker
from scipy.spatial import KDTree

from ..measure_core import (
    CSV_COLUMNS, NUCLEUS_MAP_SPLIT, NUCLEUS_MAP_UNSPLIT,
    bs_per_z, export_meta, graph_columns, wall_columns,
)
from ..saved_check import (
    curvature_notes, force_notes, mesh_notes, native_notes,
    nuclei_notes, ratio_dropout_note, surface_verts_for_distance,
    volume_shape_notes,
)
from .widgets import call_all, with_after

N_WORKERS = max(1, (os.cpu_count() or 8) - 4)

# The 2x4 panel grid at the size it is meant to be read at. The on-screen
# canvas is resized to the dock by Qt, so this is only honoured on export;
# it is what keeps a saved figure from inheriting a squashed dock.
FIG_SIZE_IN = (17.0, 7.0)


def build(state, viewer):
    fig = Figure(figsize=(17, 7), tight_layout=True)
    axes = fig.subplots(2, 4)
    canvas = FigureCanvasQTAgg(fig)

    measure_btn = PushButton(text="Measure adhesions")

    def _measure(*_, after=None):
        # ---- prerequisite checks ----
        # NB: force_vol is legitimately None when FRET isn't configured (see
        # loader.fret_configured), so "is an image loaded?" is asked of `scale`,
        # not of the force volume. The force columns just come out NaN.
        if state.get("scale") is None:
            print("Measure: load an image first.", flush=True)
            call_all(after)
            return
        if (state.get("vessel_mesh_cache") or {}).get("verts") is None:
            print("Measure: build the vessel mesh first (Vessel mesh tab).", flush=True)
            call_all(after)
            return
        ds = state.get("dapi_state") or {}
        if ds.get("nuclei") is None and ds.get("labels") is None:
            print("Measure: compute nuclei first (Nuclei tab → Compute).", flush=True)
            call_all(after)
            return
        layer_names = [lay.name for lay in viewer.layers]
        if "Adhesion instances" not in layer_names:
            print("Measure: run RF predict + Apply first "
                  "(Adhesion RF tab → Apply post-process).", flush=True)
            call_all(after)
            return

        force_vol = state.get("force_vol")
        vessel_verts = state["vessel_mesh_cache"]["verts"]
        # Faces restrict dist_mesh_um to surface vertices (P09a). Captured with
        # the verts so the two always come from the same mesh.
        vessel_faces = (state.get("vessel_mesh_cache") or {}).get("faces")
        # The watershed-split nuclei when available, else the unsplit DAPI
        # components (P04). Captured into locals here, before the worker, so
        # they are closed over rather than re-read off `state` on another
        # thread. The fallback is loud: the two maps mean different things and
        # only the split one can be joined to the fibroblast classification.
        nuc_ids = np.asarray(ds["nuclei_ids"]) if ds.get("nuclei_ids") is not None else None
        nuc_is_fib = np.asarray(ds["nuclei_is_fib"]) if ds.get("nuclei_is_fib") is not None else None
        if ds.get("nuclei") is not None:
            nuc_labels_full = ds["nuclei"]
        else:
            nuc_labels_full = ds["labels"]
            nuc_ids = nuc_is_fib = None
            print("Measure: WARN: no split nuclei — using the UNSPLIT DAPI "
                  "components for nuc_id and dist_nuc_um, and "
                  "nuc_is_fibroblast will be blank. Those are DAPI connected "
                  "components, not nuclei (most are sub-nuclear debris), and "
                  "nuc_id is then not joinable to n_fibroblasts. Run 'Count "
                  "fibroblasts' on the Nuclei tab first.", flush=True)
        scale = state.get("scale", (1.0, 1.0, 1.0))
        clover_native = state.get("clover_raw_native")
        ruby_native = state.get("ruby_raw_native")
        raw_voxel = tuple(state.get("raw_voxel_zyx_um", scale))
        curvature = (state.get("vessel_mesh_cache") or {}).get("curvature")
        graph = state.get("vessel_graph")
        # Captured here, before the worker, for the same reason the nuclei are:
        # so it is closed over rather than re-read off `state` on another
        # thread. None until the Vessel Mesh tab's wall-shear button has run.
        wall = state.get("vessel_shear")

        lyr = viewer.layers["Adhesion instances"]
        adh_data = lyr.data
        labels_raw = adh_data[0] if getattr(lyr, "multiscale", False) else adh_data

        @thread_worker
        def _run():
            from skimage.measure import regionprops_table
            import scipy.ndimage as sndi

            # Degradations are accumulated and printed by `_done`, so the same
            # string reaches the Log dock here and the terminal in the headless
            # path — see pipeline/saved_check.py.
            notes: list[str] = []

            labels = np.asarray(labels_raw).astype(np.int32)

            # Optional X-Y region-of-interest (ROI tab): keep only adhesions
            # whose centroid is inside the drawn region (any polygon outline;
            # a rectangle is just a 4-vertex polygon). All Z.
            polys = state.get("roi_polygons_yx_um")
            roi = state.get("roi_yx_um")
            if polys or roi is not None:
                ids0 = np.unique(labels[labels > 0])
                if ids0.size:
                    cen = np.asarray(sndi.center_of_mass(labels > 0, labels,
                                                         list(ids0)))
                    yc = cen[:, 1] * float(scale[1])
                    xc = cen[:, 2] * float(scale[2])
                    if polys:
                        from matplotlib.path import Path as _MplPath
                        pts = np.column_stack([yc, xc])
                        inside = np.zeros(len(ids0), dtype=bool)
                        for p in polys:
                            inside |= _MplPath(np.asarray(p)).contains_points(pts)
                        keep = ids0[inside]
                    else:
                        y0, y1, x0, x1 = roi
                        keep = ids0[(yc >= y0) & (yc <= y1)
                                    & (xc >= x0) & (xc <= x1)]
                    km = np.zeros(int(labels.max()) + 1, dtype=bool)
                    km[keep] = True
                    labels = np.where(km[labels], labels, 0).astype(np.int32)
                    print(f"ROI: kept {keep.size}/{ids0.size} adhesions inside "
                          f"the region.", flush=True)

            unique_ids = np.unique(labels[labels > 0])
            if unique_ids.size == 0:
                return None

            # Step 2 — shape properties
            props = regionprops_table(
                labels,
                spacing=scale,
                properties=[
                    "label", "area",
                    "axis_major_length", "axis_minor_length",
                    "centroid",
                ],
            )

            # Step 3 — force statistics. NaN when FRET isn't configured for this
            # dataset (no force volume was computed); the columns export blank.
            notes += force_notes(force_vol)
            if force_vol is None:
                force_mean = np.full(unique_ids.shape, np.nan, dtype=np.float64)
                force_std = np.full(unique_ids.shape, np.nan, dtype=np.float64)
            else:
                force_mean = np.array(sndi.mean(force_vol, labels, index=unique_ids))
                force_std  = np.array(sndi.standard_deviation(force_vol, labels,
                                                                index=unique_ids))

            # Step 3b — per-adhesion clover/ruby ratio on NATIVE Z (matches the
            # Adhesion RF tab). The isotropic adhesion labels are mapped down to
            # the native Z planes (nearest neighbour) so the ratio uses measured,
            # non-interpolated voxels.
            #
            # Total clover / total ruby over the adhesion, both background-
            # subtracted per Z plane first. Summing before dividing keeps one dim
            # voxel from dominating (a mean of per-voxel ratios is driven by its
            # smallest denominators), and the background subtraction stops dim
            # (deep) adhesions from drifting toward a fixed offset ratio — which
            # would masquerade as a curvature effect, since depth and curvature
            # covary on a tube. Aligned to unique_ids; NaN where an adhesion's
            # ruby is entirely at/below background.
            ratio_bad = native_notes(clover_native, ruby_native, labels.shape)
            notes += ratio_bad
            if not ratio_bad:
                Zn = clover_native.shape[0]
                zf = float(raw_voxel[0]) / float(scale[0])
                iso_idx = np.clip(np.rint(np.arange(Zn) * zf).astype(np.int64),
                                  0, labels.shape[0] - 1)
                labels_n = labels[iso_idx]                    # (Zn, Y, X)
                clover_bs = bs_per_z(clover_native)
                ruby_bs = bs_per_z(ruby_native)
                inside = labels_n > 0
                lv = labels_n[inside]
                nmax = int(labels.max())
                c_sum = np.bincount(lv, weights=clover_bs[inside],
                                    minlength=nmax + 1)
                r_sum = np.bincount(lv, weights=ruby_bs[inside],
                                    minlength=nmax + 1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    ratio_by_label = c_sum / r_sum
                ratio_by_label[r_sum <= 0] = np.nan
                ratio_mean = ratio_by_label[unique_ids].astype(np.float32)
                notes += ratio_dropout_note(int(np.isnan(ratio_mean).sum()),
                                            int(ratio_mean.size))
            else:
                ratio_mean = np.full(unique_ids.shape, np.nan, dtype=np.float32)

            # Step 4 — all adhesion voxel coordinates
            adhesion_vox = np.argwhere(labels > 0)
            adhesion_ids = labels[
                adhesion_vox[:, 0], adhesion_vox[:, 1], adhesion_vox[:, 2]
            ]
            adhesion_um = adhesion_vox * np.array(scale)

            # Step 5 — vessel mesh KD-tree (built here; no cached tree)
            notes += mesh_notes(vessel_verts)
            notes += volume_shape_notes("the nucleus label volume",
                                        np.shape(nuc_labels_full),
                                        "the adhesion label volume",
                                        labels.shape)
            # Measure to the mesh surface, not to interior points of the
            # alpha-shape solid (P09a). `mesh_nn_idx` below is unused, so this
            # changes distances only — nothing is looked up by vertex index.
            dist_verts, dist_notes = surface_verts_for_distance(vessel_verts,
                                                                vessel_faces)
            notes += dist_notes
            vessel_tree = KDTree(dist_verts)

            # Step 6 — nucleus KD-tree (subsampled every 3rd voxel)
            nuc_sub = nuc_labels_full[::3, ::3, ::3]
            nuc_vox = np.argwhere(nuc_sub > 0)
            nuc_ids_sub = nuc_sub[nuc_vox[:, 0], nuc_vox[:, 1], nuc_vox[:, 2]]
            nuc_um = nuc_vox * 3 * np.array(scale)
            nuc_tree = KDTree(nuc_um)

            # Concurrent queries (both release GIL in scipy C code)
            with ThreadPoolExecutor(max_workers=2) as ex:
                f_mesh = ex.submit(vessel_tree.query, adhesion_um)
                f_nuc  = ex.submit(nuc_tree.query, adhesion_um)
                mesh_dist_all, mesh_nn_idx = f_mesh.result()
                nuc_dist_all, nuc_nn_idx = f_nuc.result()

            nuc_ids_at_nn = nuc_ids_sub[nuc_nn_idx]

            # Step 7 — per-label reduction (vectorised)
            order = np.argsort(adhesion_ids, kind="stable")
            s_ids    = adhesion_ids[order]
            s_mesh   = mesh_dist_all[order]
            s_nuc    = nuc_dist_all[order]
            s_nuc_id = nuc_ids_at_nn[order]

            _, first = np.unique(s_ids, return_index=True)
            counts = np.diff(np.concatenate([first, [len(s_ids)]]))

            mesh_min = np.minimum.reduceat(s_mesh, first)
            nuc_min  = np.minimum.reduceat(s_nuc,  first)

            group_argmin = np.array([
                first[i] + np.argmin(s_nuc[first[i]: first[i] + counts[i]])
                for i in range(len(first))
            ])
            nuc_assignment = s_nuc_id[group_argmin]

            # Step 7a — join each adhesion's nucleus to that nucleus's
            # fibroblast class (P04). float32 + NaN, not an int sentinel: the
            # exporter blanks NaN and only NaN, so a -1 would export literally.
            nuc_is_fibroblast = np.full(len(first), np.nan, dtype=np.float32)
            if nuc_ids is not None or nuc_is_fib is not None:
                nuc_max = int(np.asarray(nuc_labels_full).max())
                nn = nuclei_notes(nuc_ids, nuc_is_fib, nuc_max)
                notes += nn
                if not nn:
                    lut = np.full(nuc_max + 1, np.nan, dtype=np.float32)
                    lut[np.asarray(nuc_ids).astype(np.int64)] = np.asarray(
                        nuc_is_fib, bool).astype(np.float32)
                    nuc_is_fibroblast = lut[
                        np.asarray(nuc_assignment).astype(np.int64)]
            else:
                notes += nuclei_notes(None, None, 0)

            # Step 7b — per-adhesion vessel curvature at the nearest *surface*
            # vertex, averaged over the adhesion's voxels. The alpha-shape vessel
            # mesh is a point cloud whose interior vertices carry no valid
            # curvature (NaN); mapping to the nearest of ALL vertices would land
            # on those. So restrict to surface vertices (finite curvature) and
            # query a tree over just those. NaN if curvature not computed/stale.
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
                curv_out[key] = (np.add.reduceat(s_cv, first)
                                 / counts).astype(np.float32)

            # Step 7c — everything read off the 3D centerline graph. **The
            # shared implementation**, not a second copy: [P11] records that
            # P04 and P12 each had to be written twice because this worker and
            # `measure_core` had diverged, and the no-drift rule in CLAUDE.md
            # says extract rather than duplicate.
            gcols, gnotes = graph_columns(graph, adhesion_um, order, first,
                                          counts)
            notes += gnotes

            # Step 7d — per-adhesion wall shear off the ragged cloud. Shared
            # with `measure_core` for the same reason as 7c.
            wcols, wnotes = wall_columns(wall, adhesion_um, order, first,
                                         counts)
            notes += wnotes

            # Step 8 — assemble result
            return {
                "label":          unique_ids,
                "volume_um3":     props["area"],
                "axis_major_um":  props["axis_major_length"],
                "axis_minor_um":  props["axis_minor_length"],
                "eccentricity":   (props["axis_major_length"]
                                   / np.maximum(props["axis_minor_length"], 1e-9)),
                "force_mean_pN":  force_mean,
                "force_std_pN":   force_std,
                "force_total_pN": force_mean * props["area"],
                "dist_mesh_um":   mesh_min,
                "dist_nuc_um":    nuc_min,
                "nuc_id":         nuc_assignment,
                # 1.0 = VnTs+ fibroblast, 0.0 = other, blank = unknown. Always
                # present so the CSV header does not depend on the join.
                "nuc_is_fibroblast": nuc_is_fibroblast,
                "ratio_clover_ruby": ratio_mean,
                # regionprops_table was given `spacing=scale`, so these are µm.
                "centroid_z_um":  props["centroid-0"],
                "centroid_y_um":  props["centroid-1"],
                "centroid_x_um":  props["centroid-2"],
                "curv_mean":      curv_out["mean"],
                "curv_k1":        curv_out["k1"],
                "curv_k2":        curv_out["k2"],
                "curv_gaussian":  curv_out["gaussian"],
                **gcols,
                **wcols,
                # Not a metric column — `_export` selects by name.
                "notes":          notes,
            }

        def _done(result):
            if result is None:
                print("Measure: no adhesion instances found.", flush=True)
                return
            state["measure_result"] = result

            # Report the silent degradations first — a blank column is easy to
            # miss in a CSV and impossible to miss here.
            for n in result.get("notes", []):
                print(f"  {n}", flush=True)

            n_adh = len(result["label"])
            n_nuc = len(np.unique(result["nuc_id"][result["nuc_id"] > 0]))
            fm = result["force_mean_pN"]
            dm = result["dist_mesh_um"]
            dn = result["dist_nuc_um"]
            vol = result["volume_um3"]
            cr = result["ratio_clover_ruby"]
            cr_f = cr[np.isfinite(cr)]
            cr_line = (f"  Clover/Ruby:     {cr_f.mean():.3f} ± {cr_f.std():.3f}"
                       if cr_f.size else "  Clover/Ruby:     n/a (no ruby data)")
            fm_f = fm[np.isfinite(fm)]
            fm_line = (f"  Force:           {fm_f.mean():.2f} ± "
                       f"{fm_f.std():.2f} pN" if fm_f.size else
                       "  Force:           n/a (FRET not configured)")

            print(
                f"Measured {n_adh} adhesions across {n_nuc} nuclei\n"
                f"{fm_line}\n"
                f"  Dist to mesh:    {dm.mean():.2f} ± {dm.std():.2f} µm\n"
                f"  Dist to nucleus: {dn.mean():.2f} ± {dn.std():.2f} µm\n"
                f"  Volume:          {vol.mean():.2f} ± {vol.std():.2f} µm³\n"
                f"{cr_line}",
                flush=True,
            )

            # Spearman correlation of the clover/ruby ratio vs curvature and
            # mesh distance (finite pairs only) — a quick read without the CSV.
            from scipy.stats import spearmanr
            crr = np.asarray(result["ratio_clover_ruby"])

            def _corr(name, x):
                x = np.asarray(x)
                m = np.isfinite(crr) & np.isfinite(x)
                if int(m.sum()) < 3 or np.ptp(x[m]) == 0:
                    print(f"  Spearman(ratio, {name}): n/a", flush=True)
                    return
                rho, p = spearmanr(crr[m], x[m])
                print(f"  Spearman(ratio, {name}) = {rho:+.3f}  "
                      f"(p={p:.1e}, n={int(m.sum())})", flush=True)

            for _nm, _key in (("mean curv", "curv_mean"), ("κ1", "curv_k1"),
                              ("κ2", "curv_k2"), ("gaussian", "curv_gaussian"),
                              ("path_curv_3d", "path_curv_1_um"),
                              ("caliber", "caliber_um"),
                              ("lumen_aspect", "lumen_aspect"),
                              ("dist_junction", "dist_junction_um"),
                              ("dist_mesh", "dist_mesh_um")):
                _corr(_nm, result[_key])

            # Plotting is display-only; skip it in headless batch mode and never
            # let a plot failure block the CSV export (which runs in `after`).
            if not state.get("headless"):
                try:
                    _update_plots(result)
                except Exception as exc:
                    print(f"Measure: plot update skipped ({exc})", flush=True)

        def _err(exc):
            import traceback
            print(f"Measure error: {exc}\n{traceback.format_exc()}", flush=True)

        print("Measuring adhesions ...", flush=True)
        w = _run()
        w.returned.connect(with_after(_done, after, "Measure"))
        w.errored.connect(with_after(_err, after, "Measure"))
        w.start()

    measure_btn.changed.connect(_measure)

    def _kde_overlay(ax, data, color, ylabel):
        """Histogram with KDE overlay (density-normalised).

        Drops non-finite values first: force is all-NaN when FRET isn't
        configured, and hist/linspace on NaN would blow up the whole figure.

        `density=True` (needed so the `gaussian_kde` PDF overlays on the same
        scale), so the y-axis is a probability density and its unit is the
        *reciprocal* of the x-axis unit -- 1/pN against force, and µm against
        curvature, which is already 1/µm. `ylabel` is passed per panel rather
        than derived because nothing here knows the x unit. The label is set
        even when the panel is empty (< 2 finite points), so an all-NaN force
        panel is still labelled rather than blank and ambiguous.
        """
        ax.set_ylabel(ylabel)
        data = np.asarray(data)
        data = data[np.isfinite(data)]
        if data.size < 2:
            return
        from scipy.stats import gaussian_kde
        ax.hist(data, bins=40, density=True, color=color, alpha=0.6)
        xs = np.linspace(data.min(), data.max(), 300)
        try:
            kde = gaussian_kde(data)
            ax.plot(xs, kde(xs), color=color, lw=2)
        except Exception:
            pass

    def _draw_panels(target_axes, r):
        """Draw the eight panels into `target_axes` (a 2x4 array).

        Shared by the on-screen canvas and the figure export so the file and
        the dock cannot show different plots -- the same doctrine as
        `saved_check.py`: never write the second copy.
        """
        axes = target_axes
        for ax in axes.flat:
            ax.clear()

        fm = r["force_mean_pN"]
        dm = r["dist_mesh_um"]
        dn = r["dist_nuc_um"]
        cr = r["ratio_clover_ruby"]

        _kde_overlay(axes[0, 0], fm, "crimson", "Density (1/pN)")
        axes[0, 0].set_xlabel("Force (pN)")
        axes[0, 0].set_title("Force distribution")

        _kde_overlay(axes[0, 1], dm, "steelblue", "Density (1/µm)")
        axes[0, 1].set_xlabel("Distance to vessel mesh (µm)")
        axes[0, 1].set_title("Mesh distance")

        _kde_overlay(axes[0, 2], dn, "darkorange", "Density (1/µm)")
        axes[0, 2].set_xlabel("Distance to nearest nucleus (µm)")
        axes[0, 2].set_title("Nucleus distance")

        _kde_overlay(axes[1, 0], cr[np.isfinite(cr)], "seagreen",
                     "Density (per unit ratio)")
        axes[1, 0].set_xlabel("Clover / Ruby ratio")
        axes[1, 0].set_title("Clover/Ruby ratio")

        finite = np.isfinite(fm) & np.isfinite(dm)
        if finite.sum() > 1:
            axes[1, 1].hexbin(dm[finite], fm[finite],
                              gridsize=30, cmap="YlOrRd", mincnt=1)
        axes[1, 1].set_xlabel("Distance to vessel mesh (µm)")
        axes[1, 1].set_ylabel("Force (pN)")
        axes[1, 1].set_title("Force vs mesh distance")

        finite_cr = np.isfinite(cr) & np.isfinite(dm)
        if finite_cr.sum() > 1:
            axes[1, 2].hexbin(dm[finite_cr], cr[finite_cr],
                              gridsize=30, cmap="YlGn", mincnt=1)
        axes[1, 2].set_xlabel("Distance to vessel mesh (µm)")
        axes[1, 2].set_ylabel("Clover / Ruby ratio")
        axes[1, 2].set_title("Clover/Ruby vs mesh distance")

        cm = r["curv_mean"]
        # x is 1/µm, so the density is in µm.
        _kde_overlay(axes[0, 3], cm[np.isfinite(cm)], "purple",
                     "Density (µm)")
        axes[0, 3].set_xlabel("Mean curvature (1/µm)")
        axes[0, 3].set_title("Vessel curvature (mean)")

        finite_cm = np.isfinite(cm) & np.isfinite(cr)
        if finite_cm.sum() > 1:
            axes[1, 3].hexbin(cm[finite_cm], cr[finite_cm],
                              gridsize=30, cmap="PuRd", mincnt=1)
        axes[1, 3].set_xlabel("Mean curvature (1/µm)")
        axes[1, 3].set_ylabel("Clover / Ruby ratio")
        axes[1, 3].set_title("Clover/Ruby vs curvature")

    def _update_plots(r):
        _draw_panels(axes, r)
        canvas.draw()

    # ---- Export per-adhesion measurements to CSV ----
    export_path_w = FileEdit(value="adhesion_measures.csv", label="export_path",
                             mode="w")
    export_btn = PushButton(text="Export measures (CSV)")

    def _export(*_):
        result = state.get("measure_result")
        if result is None:
            print("No measurements to export — click 'Measure adhesions' first.",
                  flush=True)
            return
        # Append an export-time timestamp so each export is a unique file and
        # never overwrites a previous one (e.g. ..._adhesion_measures_20260728_143052.csv).
        from datetime import datetime
        from pathlib import Path
        _base = Path(str(export_path_w.value))
        _ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        path = str(_base.with_name(f"{_base.stem}_{_ts}{_base.suffix or '.csv'}"))
        # Column order is `measure_core.CSV_COLUMNS` — shared with
        # `batch_measure` so the two exports cannot drift.
        try:
            import pandas as pd
            df = pd.DataFrame({c: np.asarray(result[c])
                               for c in CSV_COLUMNS if c in result})
            # Source image filename as the first column, so concatenated CSVs
            # from multiple images stay identifiable.
            df.insert(0, "image", state.get("image_name") or "")
            # Image-level fibroblast count (from the Nuclei tab), broadcast to
            # every row so concatenated CSVs keep each image's count. NaN if the
            # count wasn't run.
            fc = state.get("fibroblast_count") or {}
            df["n_fibroblasts"] = fc.get("fibroblast", np.nan)
            df["n_nuclei_total"] = fc.get("total", np.nan)
            df.to_csv(path, index=False)
            print(f"Exported {len(df)} adhesions -> {path}", flush=True)

            # Companion settings file: every step records the parameters it used
            # into state["run_settings"]; dump them next to the CSV (same
            # timestamped stem) so any export can be fully backtracked.
            import json
            meta = export_meta(
                image=state.get("image_name"),
                exported=_ts,
                csv_name=Path(path).name,
                n_adhesions=len(df),
                scale_um_zyx=state.get("scale", ()),
                roi_active=bool(state.get("roi_polygons_yx_um")
                                or state.get("roi_yx_um")),
                fibroblast_count=state.get("fibroblast_count"),
                nucleus_map=(
                    NUCLEUS_MAP_SPLIT
                    if (state.get("dapi_state") or {}).get("nuclei") is not None
                    else NUCLEUS_MAP_UNSPLIT),
                notes=result.get("notes", []),
                settings=state.get("run_settings") or {},
            )
            spath = str(_base.with_name(f"{_base.stem}_{_ts}_settings.json"))
            with open(spath, "w") as fh:
                json.dump(meta, fh, indent=2, default=str)
            print(f"Settings -> {spath}", flush=True)
        except Exception as exc:
            print(f"Export error: {exc}", flush=True)

    export_btn.changed.connect(_export)
    state.setdefault("save_path_widgets", []).append(
        (export_path_w, "{stem}_adhesion_measures.csv"))

    # Chainable "Run all" step: measure, then export the CSV only if the
    # measurement actually produced a result (otherwise the export would just
    # complain "run measure first" — the real problem is usually that RF
    # predict + post-process didn't create adhesion instances).
    def _run_measure_step(after=None):
        state["measure_result"] = None   # so export fires only on a fresh result
        def _export_if_measured():
            if state.get("measure_result") is not None:
                _export()
            else:
                print("Run all: measure produced no result — skipping export. "
                      "Check that RF predict + post-process created adhesion "
                      "instances (look for 'Adhesion instances: N').", flush=True)
        _measure(after=lambda: call_all(_export_if_measured, after,
                                        label="run_measure"))
    state["run_measure"] = _run_measure_step

    # ---- Export the plot panels as PNG + SVG ----
    # Deliberately NOT `fig.savefig(...)`: Qt owns the on-screen figure's size,
    # so the live canvas is whatever the dock happens to be and saving it bakes
    # the squashed dock geometry into the file (overlapping ticks, clipped
    # titles) -- and it would look different every time the window is resized.
    # Instead redraw the same panels into a fresh off-screen figure at
    # FIG_SIZE_IN, which also works when the dock is collapsed or tiny.
    fig_path_w = FileEdit(value="adhesion_figures.png", label="figure_path",
                          mode="w")
    fig_dpi_w = SpinBox(value=300, min=50, max=1200, step=50, label="PNG dpi")
    # The outer container hides labels, so a bare spinbox would be an unlabelled
    # number; nest it in a labelled row of its own.
    fig_dpi_row = Container(widgets=[fig_dpi_w], labels=True)
    fig_btn = PushButton(text="Export figures (PNG+SVG)")

    def _export_figures(*_):
        result = state.get("measure_result")
        if result is None:
            print("No measurements to plot — click 'Measure adhesions' first.",
                  flush=True)
            return
        from datetime import datetime
        from pathlib import Path
        import matplotlib
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        _base = Path(str(fig_path_w.value))
        _ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        png = _base.with_name(f"{_base.stem}_{_ts}.png")
        svg = _base.with_name(f"{_base.stem}_{_ts}.svg")
        try:
            # constrained_layout, not tight_layout: it solves the panel spacing
            # for this exact figure size, so nothing overlaps at 17x7.
            out = Figure(figsize=FIG_SIZE_IN, constrained_layout=True)
            FigureCanvasAgg(out)          # savefig needs a canvas attached
            _draw_panels(out.subplots(2, 4), result)
            title = state.get("image_name") or ""
            if title:
                out.suptitle(title, fontsize=11)
            # svg.fonttype "none" keeps labels as real text rather than
            # outlines, so the SVG stays editable in Illustrator/Inkscape.
            with matplotlib.rc_context({"svg.fonttype": "none"}):
                out.savefig(png, dpi=int(fig_dpi_w.value))
                out.savefig(svg)
            print(f"Figures -> {png}\n           {svg}", flush=True)
        except Exception as exc:
            print(f"Figure export error: {exc}", flush=True)

    fig_btn.changed.connect(_export_figures)
    state.setdefault("save_path_widgets", []).append(
        (fig_path_w, "{stem}_adhesion_figures.png"))

    controls = Container(widgets=[measure_btn, export_path_w, export_btn,
                                  fig_path_w, fig_dpi_row, fig_btn],
                        labels=False)
    return controls, canvas
