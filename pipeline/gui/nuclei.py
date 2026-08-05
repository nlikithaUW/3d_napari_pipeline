"""DAPI Otsu + size-filter tab."""

from __future__ import annotations

import numpy as np
from magicgui import magicgui
from magicgui.widgets import (
    CheckBox, Container, FileEdit, FloatSlider, PushButton,
)
from napari.qt.threading import thread_worker

from ..denoise import auto_min_size_otsu, filter_by_size, label_and_sizes, otsu_mask
from .widgets import Histogram, debounce, pyramid_views


def build(state, viewer) -> Container:
    hist = Histogram("Component size histogram", "log10(size in voxels)", "steelblue")
    state.setdefault("hists", []).append(hist)

    def _redraw(min_size: int):
        sizes = state["dapi_state"]["sizes"]
        hist.update(sizes, threshold=min_size if min_size > 0 else None,
                    threshold_label=f"min_size={min_size}")

    @magicgui(
        auto_call=True,
        min_size={"widget_type": "Slider", "min": 0, "max": 5000},
    )
    def size_controls(min_size: int = 0):
        _apply(min_size)

    @debounce(150)
    def _apply(min_size: int):
        labels = state["dapi_state"]["labels"]
        sizes = state["dapi_state"]["sizes"]
        layer = state.get("dapi_mask_layer")
        if labels is None or layer is None:
            return
        layer.data = pyramid_views(filter_by_size(labels, sizes, min_size),
                                   state.get("levels", 1))
        layer.visible = True
        _redraw(min_size)

    def _run_dapi(after=None):
        dapi_raw = state.get("dapi_raw")
        if dapi_raw is None:
            print("Load a stack with a 'nuclei' channel first.", flush=True)
            if after:
                after()
            return

        @thread_worker
        def _run():
            mask = otsu_mask(dapi_raw)
            labels, sizes = label_and_sizes(mask)
            keep = sizes >= 10
            remap = np.zeros(sizes.size + 1, dtype=np.int32)
            remap[1:][keep] = np.arange(1, int(keep.sum()) + 1)
            labels = remap[labels]
            sizes = sizes[keep]
            suggested = auto_min_size_otsu(sizes)
            return labels, sizes, suggested

        def _done(result):
            labels, sizes, suggested = result
            state["dapi_state"]["labels"] = labels
            state["dapi_state"]["sizes"] = sizes
            print(f"DAPI: {sizes.size} components. Suggested min_size = {suggested}",
                  flush=True)
            slider_max = int(sizes.max()) if sizes.size else 5000
            size_controls.min_size.max = slider_max
            size_controls.min_size.value = suggested
            size_controls()
            if after:
                after()

        def _err(exc):
            print(f"DAPI Otsu error: {exc}", flush=True)
            if after:
                after()

        print("DAPI Otsu + labeling running ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    @magicgui(call_button="Compute DAPI Otsu + components")
    def otsu_controls():
        _run_dapi()

    # ---- Save / load the DAPI components (labels + sizes) ----
    # Persists the full result of "Compute DAPI Otsu + components" so it can be
    # restored on reopen — the labels feed Measure and the DAPI gate, and the
    # sizes drive the min_size slider + histogram and the displayed mask.
    labels_path_w = FileEdit(value="dapi_labels.npz", label="labels_path",
                             mode="w")
    save_lbl_w = PushButton(text="Save DAPI labels")
    load_lbl_w = PushButton(text="Load DAPI labels")

    def _save_labels(*_):
        ds = state.get("dapi_state") or {}
        labels = ds.get("labels")
        sizes = ds.get("sizes")
        if labels is None or sizes is None:
            print("No DAPI components to save — run 'Compute DAPI Otsu + "
                  "components' first.", flush=True)
            return
        path = str(labels_path_w.value)
        np.savez_compressed(path, labels=np.asarray(labels).astype(np.int32),
                            sizes=np.asarray(sizes))
        print(f"Saved DAPI components -> {path} "
              f"({np.asarray(sizes).size} components).", flush=True)

    def _load_labels(*_, after=None, visible=True):
        path = str(labels_path_w.value)

        @thread_worker
        def _run():
            d = np.load(path)
            return d["labels"].astype(np.int32), d["sizes"]

        def _done(result):
            labels, sizes = result
            dapi_raw = state.get("dapi_raw")
            if dapi_raw is not None and labels.shape != dapi_raw.shape:
                print(f"Warning: loaded labels {labels.shape} != current DAPI "
                      f"{dapi_raw.shape}; load the matching image first.",
                      flush=True)
            state["dapi_state"]["labels"] = labels
            state["dapi_state"]["sizes"] = sizes
            suggested = auto_min_size_otsu(sizes)
            size_controls.min_size.max = int(sizes.max()) if sizes.size else 5000
            size_controls.min_size.value = suggested
            size_controls()  # refresh the mask layer + histogram
            if not visible and state.get("dapi_mask_layer") is not None:
                state["dapi_mask_layer"].visible = False
            print(f"Loaded DAPI components <- {path}: {sizes.size} components, "
                  f"suggested min_size = {suggested}.", flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Load DAPI labels error: {exc}", flush=True)
            if after:
                after()

        print(f"Loading DAPI components from {path} ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    save_lbl_w.changed.connect(_save_labels)
    load_lbl_w.changed.connect(_load_labels)
    state["load_dapi_labels"] = _load_labels
    state.setdefault("save_path_widgets", []).append(
        (labels_path_w, "dapi_labels.npz"))

    # ---- Count fibroblasts: VnTs (mClover)-positive nuclei ----
    # Splits touching nuclei (watershed on the DAPI distance transform), then
    # marks a nucleus as a fibroblast if its perinuclear mClover (VnTs) signal is
    # high; endothelial cells (VE-Cadherin+, VnTs-negative) fall below.
    #
    # NOTE: the auto threshold is *background-relative* (mClover_bg × factor), not
    # Otsu. Otsu assumes two well-separated intensity modes, but perinuclear
    # mClover here is only ~1.4× background and essentially unimodal, so Otsu's
    # threshold is unstable — it swung to 2% on unsplit nuclei and ~98% after the
    # watershed split. A background multiple is stable and interpretable.
    split_dist_w = FloatSlider(value=4.0, min=1.0, max=15.0, step=0.5,
                               label="nucleus_min_dist_um")
    vnts_factor_w = FloatSlider(value=1.5, min=1.0, max=4.0, step=0.1,
                                label="vnts_bg_factor (auto)")
    vnts_thresh_w = FloatSlider(value=0.0, min=0.0, max=10000.0, step=25.0,
                                label="vnts_threshold (0=auto→bg×factor)")
    exclude_cad_w = CheckBox(value=True,
                             label="exclude VE-cadherin+ (endothelial)")
    count_btn = PushButton(text="Count fibroblasts (VnTs+ nuclei)")

    def _count(*_, after=None):
        clover = state.get("clover_raw")
        dapi_layer = state.get("dapi_mask_layer")
        if clover is None:
            print("Count: load a stack with a 'clover' channel first.", flush=True)
            if after:
                after()
            return
        ds = state.get("dapi_state") or {}
        if ds.get("labels") is None:
            print("Count: run 'Compute DAPI Otsu + components' first.", flush=True)
            if after:
                after()
            return
        # Build the binary mask from the labels + min_size directly (not the
        # displayed mask layer, whose update is debounced 150 ms and so may not
        # be ready yet inside the fast "Run all" chain → would give zero nuclei).
        binary = filter_by_size(ds["labels"], ds["sizes"],
                                int(size_controls.min_size.value)) > 0
        scale = state.get("scale", (1.0, 1.0, 1.0))
        min_dist_um = float(split_dist_w.value)
        thr_manual = float(vnts_thresh_w.value)
        factor = float(vnts_factor_w.value)
        exclude_cad = bool(exclude_cad_w.value)
        cadh = state.get("cadh_raw") if exclude_cad else None
        state.setdefault("run_settings", {})["fibroblast"] = {
            "nucleus_min_dist_um": min_dist_um, "vnts_bg_factor": factor,
            "vnts_threshold_manual": thr_manual, "exclude_cadherin": exclude_cad,
            "dapi_min_size_vox": int(size_controls.min_size.value)}

        @thread_worker
        def _run():
            from scipy import ndimage as ndi
            from skimage.feature import peak_local_max
            from skimage.segmentation import watershed
            sampling = tuple(float(s) for s in scale)
            dist = ndi.distance_transform_edt(binary, sampling=sampling)
            md = max(1, int(round(min_dist_um / float(scale[0]))))
            coords = peak_local_max(dist, min_distance=md, labels=binary)
            markers = np.zeros(binary.shape, dtype=np.int32)
            if len(coords):
                markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)
            nuc = watershed(-dist, markers, mask=binary)
            ids = np.unique(nuc[nuc > 0])
            if ids.size == 0:
                return nuc, ids, np.zeros(0, bool), 0.0, 0.0
            # Perinuclear signal per nucleus (~2 µm shell around each nucleus):
            # mClover (VnTs) always, VE-cadherin too when excluding endothelial.
            cf = clover.astype(np.float32)
            cad = np.asarray(cadh).astype(np.float32) if cadh is not None else None
            r = max(1, int(round(2.0 / float(scale[0]))))
            objs = ndi.find_objects(nuc)
            means = np.zeros(ids.size, dtype=np.float32)
            cad_means = np.zeros(ids.size, dtype=np.float32)
            for k, i in enumerate(ids):
                sl = objs[i - 1]
                if sl is None:
                    continue
                sl2 = tuple(slice(max(0, s.start - r), s.stop + r) for s in sl)
                core = nuc[sl2] == i
                shell = ndi.binary_dilation(core, iterations=r) & ~core
                vals = cf[sl2][shell]
                means[k] = float(vals.mean()) if vals.size else 0.0
                if cad is not None:
                    cvals = cad[sl2][shell]
                    cad_means[k] = float(cvals.mean()) if cvals.size else 0.0
            # Background-relative auto threshold (stable, unlike Otsu here).
            mbg = float(np.median(cf))
            T = thr_manual if thr_manual > 0 else mbg * factor
            is_fib = means > T
            if cad is not None:
                cbg = float(np.median(cad))
                is_fib &= cad_means < cbg * 1.3   # drop VE-cadherin+ endothelial
            return nuc, ids, is_fib, T, mbg

        def _done(result):
            nuc, ids, is_fib, T, mbg = result
            n_tot = int(ids.size)
            n_fib = int(is_fib.sum())
            state["fibroblast_count"] = {"total": n_tot, "fibroblast": n_fib,
                                         "other": n_tot - n_fib, "threshold": T}
            # NB: do NOT write the computed auto threshold back into the manual
            # slider — that would make it "sticky" and be reused (as a manual
            # value, with the wrong background) on the next image.
            # class map: 1 = fibroblast (VnTs+), 2 = other; single-scale Labels.
            remap = np.zeros(int(nuc.max()) + 1, dtype=np.uint8)
            remap[ids] = np.where(is_fib, 1, 2).astype(np.uint8)
            name = "Nuclei class (1=VnTs+ fibroblast, 2=other)"
            if name in [l.name for l in viewer.layers]:
                viewer.layers.remove(name)
            viewer.add_labels(remap[nuc], name=name, scale=scale)
            thr_kind = ("manual" if thr_manual > 0
                        else f"auto = {mbg:.0f}×{float(vnts_factor_w.value):.1f}")
            cad_note = " + VE-cad exclusion" if exclude_cad else ""
            print(f"Fibroblasts: {n_fib} VnTs+ of {n_tot} nuclei "
                  f"({100.0 * n_fib / max(n_tot, 1):.0f}%); {n_tot - n_fib} other "
                  f"(mClover threshold = {T:.0f} [{thr_kind}]{cad_note}; "
                  f"mClover bg = {mbg:.0f}).", flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Count error: {exc}", flush=True)
            if after:
                after()

        print("Counting fibroblasts (splitting nuclei + classifying) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    count_btn.changed.connect(_count)

    # ---- Save / load the fibroblast count so the CSV's n_fibroblasts column can
    # be restored without re-running (the count itself is not otherwise saved).
    fib_path_w = FileEdit(value="fibroblast_count.npz", label="fib_count_path",
                          mode="w")
    save_fib_w = PushButton(text="Save fibroblast count")
    load_fib_w = PushButton(text="Load fibroblast count")

    def _save_fib(*_):
        fc = state.get("fibroblast_count")
        if not fc:
            print("No fibroblast count to save — click 'Count fibroblasts' "
                  "first.", flush=True)
            return
        path = str(fib_path_w.value)
        np.savez(path, total=int(fc.get("total", 0)),
                 fibroblast=int(fc.get("fibroblast", 0)),
                 other=int(fc.get("other", 0)),
                 threshold=float(fc.get("threshold", 0.0)))
        print(f"Saved fibroblast count -> {path} "
              f"({fc.get('fibroblast', 0)}/{fc.get('total', 0)}).", flush=True)

    def _load_fib(*_, after=None, visible=True):
        path = str(fib_path_w.value)
        try:
            d = np.load(path)
            fc = {"total": int(d["total"]), "fibroblast": int(d["fibroblast"]),
                  "other": int(d["other"]), "threshold": float(d["threshold"])}
        except Exception as exc:
            print(f"Load fibroblast count error: {exc}", flush=True)
            if after:
                after()
            return
        state["fibroblast_count"] = fc
        print(f"Loaded fibroblast count <- {path}: {fc['fibroblast']} VnTs+ of "
              f"{fc['total']} nuclei.", flush=True)
        if after:
            after()

    save_fib_w.changed.connect(_save_fib)
    load_fib_w.changed.connect(_load_fib)
    state["load_fibroblast_count"] = _load_fib

    # Chainable "Run all" steps: compute then save, then call `after`.
    state["run_dapi"] = lambda after=None: _run_dapi(
        after=(lambda: (_save_labels(), after() if after else None)))
    state["run_fibroblasts"] = lambda after=None: _count(
        after=(lambda: (_save_fib(), after() if after else None)))
    state.setdefault("save_path_widgets", []).append(
        (fib_path_w, "fibroblast_count.npz"))

    state["dapi_size_controls"] = size_controls
    # Group the fibroblast-count parameters in a labelled sub-container so each
    # slider / checkbox shows its name (the outer panel hides labels).
    fib_params = Container(widgets=[split_dist_w, vnts_factor_w, vnts_thresh_w,
                                    exclude_cad_w], labels=True)
    return (Container(widgets=[otsu_controls, size_controls,
                              labels_path_w, save_lbl_w, load_lbl_w,
                              fib_params, count_btn,
                              fib_path_w, save_fib_w, load_fib_w],
                     labels=False),
            hist.canvas)
