"""Measure tab — per-adhesion force, shape, and distance metrics."""

from __future__ import annotations

import os

import numpy as np
from concurrent.futures import ThreadPoolExecutor
from magicgui.widgets import Container, FileEdit, PushButton
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from napari.qt.threading import thread_worker
from scipy.spatial import KDTree

N_WORKERS = max(1, (os.cpu_count() or 8) - 4)


def build(state, viewer):
    fig = Figure(figsize=(17, 7), tight_layout=True)
    axes = fig.subplots(2, 4)
    canvas = FigureCanvasQTAgg(fig)

    measure_btn = PushButton(text="Measure adhesions")

    def _measure(*_, after=None):
        # ---- prerequisite checks ----
        if state.get("force_vol") is None:
            print("Measure: load an image first (force_vol not computed).", flush=True)
            if after:
                after()
            return
        if (state.get("vessel_mesh_cache") or {}).get("verts") is None:
            print("Measure: build the vessel mesh first (Vessel mesh tab).", flush=True)
            if after:
                after()
            return
        if (state.get("dapi_state") or {}).get("labels") is None:
            print("Measure: compute nuclei first (Nuclei tab → Compute).", flush=True)
            if after:
                after()
            return
        layer_names = [lay.name for lay in viewer.layers]
        if "Adhesion instances" not in layer_names:
            print("Measure: run RF predict + Apply first "
                  "(Adhesion RF tab → Apply post-process).", flush=True)
            if after:
                after()
            return

        force_vol = state["force_vol"]
        vessel_verts = state["vessel_mesh_cache"]["verts"]
        nuc_labels_full = state["dapi_state"]["labels"]
        scale = state.get("scale", (1.0, 1.0, 1.0))
        clover_native = state.get("clover_raw_native")
        ruby_native = state.get("ruby_raw_native")
        raw_voxel = tuple(state.get("raw_voxel_zyx_um", scale))
        curvature = (state.get("vessel_mesh_cache") or {}).get("curvature")
        centerline = state.get("vessel_centerline")

        lyr = viewer.layers["Adhesion instances"]
        adh_data = lyr.data
        labels_raw = adh_data[0] if getattr(lyr, "multiscale", False) else adh_data

        @thread_worker
        def _run():
            from skimage.measure import regionprops_table
            import scipy.ndimage as sndi

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

            # Step 3 — force statistics
            force_mean = np.array(sndi.mean(force_vol, labels, index=unique_ids))
            force_std  = np.array(sndi.standard_deviation(force_vol, labels,
                                                            index=unique_ids))

            # Step 3b — per-adhesion clover/ruby ratio on NATIVE Z (matches the
            # Adhesion RF tab). The isotropic adhesion labels are mapped down to
            # the native Z planes (nearest neighbour) so the ratio averages
            # measured, non-interpolated voxels. Aligned to unique_ids; NaN where
            # an adhesion has no ruby > 0 native voxel.
            if (clover_native is not None and ruby_native is not None
                    and clover_native.shape[1:] == labels.shape[1:]):
                Zn = clover_native.shape[0]
                zf = float(raw_voxel[0]) / float(scale[0])
                iso_idx = np.clip(np.rint(np.arange(Zn) * zf).astype(np.int64),
                                  0, labels.shape[0] - 1)
                labels_n = labels[iso_idx]                    # (Zn, Y, X)
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

            # Step 4 — all adhesion voxel coordinates
            adhesion_vox = np.argwhere(labels > 0)
            adhesion_ids = labels[
                adhesion_vox[:, 0], adhesion_vox[:, 1], adhesion_vox[:, 2]
            ]
            adhesion_um = adhesion_vox * np.array(scale)

            # Step 5 — vessel mesh KD-tree (built here; no cached tree)
            vessel_tree = KDTree(vessel_verts)

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

            # Step 7b — per-adhesion vessel curvature at the nearest *surface*
            # vertex, averaged over the adhesion's voxels. The alpha-shape vessel
            # mesh is a point cloud whose interior vertices carry no valid
            # curvature (NaN); mapping to the nearest of ALL vertices would land
            # on those. So restrict to surface vertices (finite curvature) and
            # query a tree over just those. NaN if curvature not computed/stale.
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
                curv_out[key] = (np.add.reduceat(s_cv, first)
                                 / counts).astype(np.float32)

            # Step 7c — in-plane (X-Y) vessel path curvature per adhesion: κ at
            # the nearest centerline point (averaged over the adhesion's voxels)
            # and the X-Y distance to the nearest vessel branch node.
            n_ad = len(first)
            path_curv = np.full(n_ad, np.nan, dtype=np.float32)
            dist_branch = np.full(n_ad, np.nan, dtype=np.float32)
            if centerline is not None and len(centerline.get("pts_yx_um", [])):
                adh_yx = adhesion_um[:, 1:3]
                cl_tree = KDTree(np.asarray(centerline["pts_yx_um"]))
                _, cnn = cl_tree.query(adh_yx)
                s_kap = np.asarray(centerline["kappa"])[cnn][order].astype(np.float64)
                path_curv = (np.add.reduceat(s_kap, first)
                             / counts).astype(np.float32)
                branch = np.asarray(centerline.get("branch_yx_um", []))
                if len(branch):
                    bd, _ = KDTree(branch).query(adh_yx)
                    dist_branch = np.minimum.reduceat(bd[order],
                                                      first).astype(np.float32)

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
                "ratio_clover_ruby": ratio_mean,
                "curv_mean":      curv_out["mean"],
                "curv_k1":        curv_out["k1"],
                "curv_k2":        curv_out["k2"],
                "curv_gaussian":  curv_out["gaussian"],
                "path_curv_xy":   path_curv,
                "dist_branch_um": dist_branch,
            }

        def _done(result):
            if result is None:
                print("Measure: no adhesion instances found.", flush=True)
                if after:
                    after()
                return
            state["measure_result"] = result

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

            print(
                f"Measured {n_adh} adhesions across {n_nuc} nuclei\n"
                f"  Force:           {fm.mean():.2f} ± {fm.std():.2f} pN\n"
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
                              ("path_curv_xy", "path_curv_xy"),
                              ("dist_branch", "dist_branch_um"),
                              ("dist_mesh", "dist_mesh_um")):
                _corr(_nm, result[_key])

            # Plotting is display-only; skip it in headless batch mode and never
            # let a plot failure block the CSV export (which runs in `after`).
            if not state.get("headless"):
                try:
                    _update_plots(result)
                except Exception as exc:
                    print(f"Measure: plot update skipped ({exc})", flush=True)
            if after:
                after()

        def _err(exc):
            import traceback
            print(f"Measure error: {exc}\n{traceback.format_exc()}", flush=True)
            if after:
                after()

        print("Measuring adhesions ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    measure_btn.changed.connect(_measure)

    def _kde_overlay(ax, data, color):
        """Histogram with KDE overlay (density-normalised)."""
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

    def _update_plots(r):
        for ax in axes.flat:
            ax.clear()

        fm = r["force_mean_pN"]
        dm = r["dist_mesh_um"]
        dn = r["dist_nuc_um"]
        cr = r["ratio_clover_ruby"]

        _kde_overlay(axes[0, 0], fm, "crimson")
        axes[0, 0].set_xlabel("Force (pN)")
        axes[0, 0].set_title("Force distribution")

        _kde_overlay(axes[0, 1], dm, "steelblue")
        axes[0, 1].set_xlabel("Distance to vessel mesh (µm)")
        axes[0, 1].set_title("Mesh distance")

        _kde_overlay(axes[0, 2], dn, "darkorange")
        axes[0, 2].set_xlabel("Distance to nearest nucleus (µm)")
        axes[0, 2].set_title("Nucleus distance")

        _kde_overlay(axes[1, 0], cr[np.isfinite(cr)], "seagreen")
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
        _kde_overlay(axes[0, 3], cm[np.isfinite(cm)], "purple")
        axes[0, 3].set_xlabel("Mean curvature (1/µm)")
        axes[0, 3].set_title("Vessel curvature (mean)")

        finite_cm = np.isfinite(cm) & np.isfinite(cr)
        if finite_cm.sum() > 1:
            axes[1, 3].hexbin(cm[finite_cm], cr[finite_cm],
                              gridsize=30, cmap="PuRd", mincnt=1)
        axes[1, 3].set_xlabel("Mean curvature (1/µm)")
        axes[1, 3].set_ylabel("Clover / Ruby ratio")
        axes[1, 3].set_title("Clover/Ruby vs curvature")

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
        # Requested columns first (clover/ruby ratio, distance to vessel mesh),
        # then the rest of the per-adhesion metrics.
        cols = ["label", "ratio_clover_ruby", "dist_mesh_um",
                "path_curv_xy", "dist_branch_um",
                "curv_mean", "curv_k1", "curv_k2", "curv_gaussian",
                "dist_nuc_um", "force_mean_pN", "force_std_pN", "force_total_pN",
                "volume_um3", "axis_major_um", "axis_minor_um",
                "eccentricity", "nuc_id"]
        try:
            import pandas as pd
            df = pd.DataFrame({c: np.asarray(result[c])
                               for c in cols if c in result})
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
            meta = {
                "image": state.get("image_name"),
                "exported": _ts,
                "csv": Path(path).name,
                "n_adhesions": int(len(df)),
                "scale_um_zyx": list(state.get("scale", ())),
                "roi_active": bool(state.get("roi_polygons_yx_um")
                                   or state.get("roi_yx_um")),
                "fibroblast_count": state.get("fibroblast_count"),
                "settings": state.get("run_settings") or {},
            }
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
        def _after():
            if state.get("measure_result") is not None:
                _export()
            else:
                print("Run all: measure produced no result — skipping export. "
                      "Check that RF predict + post-process created adhesion "
                      "instances (look for 'Adhesion instances: N').", flush=True)
            if after:
                after()
        _measure(after=_after)
    state["run_measure"] = _run_measure_step

    controls = Container(widgets=[measure_btn, export_path_w, export_btn],
                        labels=False)
    return controls, canvas
