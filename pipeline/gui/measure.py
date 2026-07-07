"""Measure tab — per-adhesion force, shape, and distance metrics."""

from __future__ import annotations

import os

import numpy as np
from concurrent.futures import ThreadPoolExecutor
from magicgui.widgets import Container, PushButton
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from napari.qt.threading import thread_worker
from scipy.spatial import KDTree

N_WORKERS = max(1, (os.cpu_count() or 8) - 4)


def build(state, viewer):
    fig = Figure(figsize=(9, 7), tight_layout=True)
    axes = fig.subplots(2, 2)
    canvas = FigureCanvasQTAgg(fig)

    measure_btn = PushButton(text="Measure adhesions")

    def _measure(*_):
        # ---- prerequisite checks ----
        if state.get("force_vol") is None:
            print("Measure: load an image first (force_vol not computed).", flush=True)
            return
        if (state.get("vessel_mesh_cache") or {}).get("verts") is None:
            print("Measure: build the vessel mesh first (Vessel mesh tab).", flush=True)
            return
        if (state.get("dapi_state") or {}).get("labels") is None:
            print("Measure: compute nuclei first (Nuclei tab → Compute).", flush=True)
            return
        layer_names = [lay.name for lay in viewer.layers]
        if "Adhesion instances" not in layer_names:
            print("Measure: run RF predict + Apply first "
                  "(Adhesion RF tab → Apply post-process).", flush=True)
            return

        force_vol = state["force_vol"]
        vessel_verts = state["vessel_mesh_cache"]["verts"]
        nuc_labels_full = state["dapi_state"]["labels"]
        scale = state.get("scale", (1.0, 1.0, 1.0))

        lyr = viewer.layers["Adhesion instances"]
        adh_data = lyr.data
        labels_raw = adh_data[0] if getattr(lyr, "multiscale", False) else adh_data

        @thread_worker
        def _run():
            from skimage.measure import regionprops_table
            import scipy.ndimage as sndi

            labels = np.asarray(labels_raw).astype(np.int32)
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
                mesh_dist_all, _ = f_mesh.result()
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
            }

        def _done(result):
            if result is None:
                print("Measure: no adhesion instances found.", flush=True)
                return

            n_adh = len(result["label"])
            n_nuc = len(np.unique(result["nuc_id"][result["nuc_id"] > 0]))
            fm = result["force_mean_pN"]
            dm = result["dist_mesh_um"]
            dn = result["dist_nuc_um"]
            vol = result["volume_um3"]

            print(
                f"Measured {n_adh} adhesions across {n_nuc} nuclei\n"
                f"  Force:           {fm.mean():.2f} ± {fm.std():.2f} pN\n"
                f"  Dist to mesh:    {dm.mean():.2f} ± {dm.std():.2f} µm\n"
                f"  Dist to nucleus: {dn.mean():.2f} ± {dn.std():.2f} µm\n"
                f"  Volume:          {vol.mean():.2f} ± {vol.std():.2f} µm³",
                flush=True,
            )

            _update_plots(result)

        def _err(exc):
            import traceback
            print(f"Measure error: {exc}\n{traceback.format_exc()}", flush=True)

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

        _kde_overlay(axes[0, 0], fm, "crimson")
        axes[0, 0].set_xlabel("Force (pN)")
        axes[0, 0].set_title("Force distribution")

        _kde_overlay(axes[0, 1], dm, "steelblue")
        axes[0, 1].set_xlabel("Distance to vessel mesh (µm)")
        axes[0, 1].set_title("Mesh distance")

        _kde_overlay(axes[1, 0], dn, "darkorange")
        axes[1, 0].set_xlabel("Distance to nearest nucleus (µm)")
        axes[1, 0].set_title("Nucleus distance")

        finite = np.isfinite(fm) & np.isfinite(dm)
        if finite.sum() > 1:
            axes[1, 1].hexbin(dm[finite], fm[finite],
                              gridsize=30, cmap="YlOrRd", mincnt=1)
        axes[1, 1].set_xlabel("Distance to vessel mesh (µm)")
        axes[1, 1].set_ylabel("Force (pN)")
        axes[1, 1].set_title("Force vs mesh distance")

        canvas.draw()

    controls = Container(widgets=[measure_btn], labels=False)
    return controls, canvas
