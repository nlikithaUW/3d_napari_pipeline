"""Unified Cadherin tab.

Per run:
  local threshold (Niblack / Sauvola / Bernsen) ANDed with a global intensity
  floor (median + k·MAD, auto-seeded from threshold_triangle on load)
  → labeled components → per-component centroid + PCA normal.

Sliders then AND a size filter with a band filter (|signed distance to the
current DAPI cell mesh| ≤ band_um); the band filter is skipped if no DAPI mesh
exists yet. Centroids/normals feed the Cell mesh "Refine with cadherin" step;
the binary mask feeds the Vessel mesh tab.
"""

from __future__ import annotations

import time

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container, FileEdit, PushButton
from napari.qt.threading import thread_worker

from ..cell_mesh import label_centroids_normals, signed_distance_to_mesh_vertices
from ..denoise import (
    auto_min_size_otsu,
    bernsen_mask,
    label_and_sizes,
    niblack_mask,
    robust_floor_value,
    sauvola_mask,
    suggest_floor_k,
)
from .widgets import Histogram, call_all, debounce, pyramid_views, with_after

METHOD_NIBLACK = "niblack"
METHOD_SAUVOLA = "sauvola"
METHOD_BERNSEN = "bernsen"
METHODS = [METHOD_NIBLACK, METHOD_SAUVOLA, METHOD_BERNSEN]

# Floor for the auto-seeded despeckle size filter. `auto_min_size_otsu` returns
# ~5 on the WT stack, which leaves ~99% of (speckle) components in the mask and
# wrecks the downstream vessel alpha shape; ~12 drops the speckle. Fully tunable
# on the slider afterwards.
MIN_SIZE_FLOOR = 12


def build(state, viewer) -> Container:
    hist = Histogram("Cadherin component sizes",
                     "log10(size in voxels)", "deeppink")
    state.setdefault("hists", []).append(hist)

    no_mesh_notified = {"done": False}

    def _redraw(min_size: int):
        sizes = state["cadh_state"]["sizes"]
        hist.update(sizes, threshold=min_size if min_size > 0 else None,
                    threshold_label=f"min_size={min_size}")

    def _ensure_centroids():
        """Per-component centroids/normals, computed on first use and cached.

        This is the expensive half of the Cadherin tab and it is **not** on the
        main path. `label_centroids_normals` is a Python loop over every
        component doing `argwhere` + covariance + `eigh` each, and the tab calls
        it with `min_voxels=1`, so it runs on every 1-voxel speck. Measured on
        the reference image, against a constant ~3.7 s for `label_and_sizes`:

            floor_k 5.0 -> 170,236 components ->  4.98 s
            floor_k 4.0 -> 375,467 components ->  9.92 s
            floor_k 3.0 -> 1,335,402 components -> 34.45 s

        Component count explodes as the floor drops (7.8x for 3.4x the voxels),
        so this dominated every run and made a permissive `floor_k` unusable.

        Nothing on the main path reads the result: the vessel mesh consumes the
        binary *mask*, not the labels. The only consumer in the codebase is
        `gui/cell_mesh.py`'s optional "Refine with cadherin" button, which is
        not in `batch.py`'s `_RUN_STEPS` — so on every headless run this was
        computed and thrown away.

        Kept memoised rather than deleted because the band filter needs it too.
        """
        cs = state["cadh_state"]
        if cs.get("centroids") is not None:
            return cs["centroids"], cs["normals"], cs["label_ids"]
        labels = cs.get("labels")
        if labels is None:
            return None, None, None
        n = cs["sizes"].size
        print(f"Cadherin: computing centroids/normals for {n:,} components "
              f"(first use) ...", flush=True)
        t0 = time.perf_counter()
        c, nm, lids = label_centroids_normals(labels, state["scale"],
                                              min_voxels=1)
        cs["centroids"], cs["normals"], cs["label_ids"] = c, nm, lids
        print(f"Cadherin: centroids DONE in {time.perf_counter() - t0:.1f} s "
              f"({len(c):,} of {n:,} components).", flush=True)
        return c, nm, lids

    def _cadh_vertices():
        """(centroids, normals, keep_mask), materialising all three on demand.

        The single entry point for consumers outside this tab, so they never
        have to know the fields are lazy or how `keep_mask` is derived.
        """
        c, nm, _ = _ensure_centroids()
        if c is None:
            return None, None, None
        cs = state["cadh_state"]
        if cs.get("keep_mask") is None:
            _compute_mask(int(filter_controls.min_size.value),
                          float(filter_controls.band_um.value))
        return cs["centroids"], cs["normals"], cs["keep_mask"]

    state["cadh_vertices"] = _cadh_vertices
    state["cadh_ensure_centroids"] = _ensure_centroids

    def _compute_mask(min_size: int, band_um: float):
        cs = state["cadh_state"]
        labels = cs["labels"]
        sizes = cs["sizes"]
        layer = state.get("cadh_mask_layer")
        if labels is None or layer is None:
            return None
        n_labels = sizes.size
        size_pass = sizes >= max(int(min_size), 0)
        keep_label = np.zeros(n_labels + 1, dtype=bool)
        keep_label[1:] = size_pass

        mesh = state.get("mesh_cache") or {}
        kdtree = mesh.get("kdtree")
        offset_centroids = mesh.get("offset_centroids")
        outward = mesh.get("outward")
        # Only materialise centroids if the band filter will actually use them.
        # With band_um == 0 (the default) this stays lazy.
        want_band = (band_um > 0 and kdtree is not None
                     and offset_centroids is not None and outward is not None
                     and len(offset_centroids) > 0)
        centroids, lids = cs.get("centroids"), cs.get("label_ids")
        if want_band:
            centroids, _n, lids = _ensure_centroids()
        use_band = want_band and centroids is not None and len(centroids) > 0

        if use_band:
            sd = signed_distance_to_mesh_vertices(
                centroids, offset_centroids, outward, tree=kdtree
            )
            band_pass = np.abs(sd) <= float(band_um)
            surviving = keep_label[lids] & band_pass
            keep_label = np.zeros(n_labels + 1, dtype=bool)
            keep_label[lids] = surviving
            keep_mask = surviving
        else:
            if band_um > 0 and not no_mesh_notified["done"]:
                print("Cadherin: no DAPI mesh yet — band filter skipped "
                      "(build the cell mesh first).", flush=True)
                no_mesh_notified["done"] = True
            # `lids` is None while centroids are unmaterialised. Leave
            # `keep_mask` as None rather than forcing the expensive computation
            # here — `_cadh_vertices` fills it in when a consumer actually asks.
            if lids is None:
                keep_mask = None
            else:
                keep_mask = (size_pass[lids - 1] if lids.size
                             else np.zeros(0, dtype=bool))

        cs["keep_mask"] = keep_mask
        return keep_label[labels].astype(np.uint8)

    @magicgui(
        auto_call=True,
        min_size={"widget_type": "Slider", "min": 0, "max": 5000},
        band_um={"widget_type": "FloatSlider", "min": 0.0, "max": 50.0, "step": 0.5},
    )
    def filter_controls(min_size: int = 0, band_um: float = 0.0):
        _apply(min_size, band_um)

    @debounce(150)
    def _apply(min_size: int, band_um: float):
        layer = state.get("cadh_mask_layer")
        if layer is None or state["cadh_state"]["labels"] is None:
            return
        mask = _compute_mask(min_size, band_um)
        if mask is None:
            return
        layer.data = pyramid_views(mask, state.get("levels", 1))
        layer.visible = True
        _redraw(min_size)

    @magicgui(
        call_button="Apply threshold + label",
        method={"widget_type": "ComboBox", "choices": METHODS},
        window_xy={"widget_type": "Slider", "min": 5, "max": 101, "step": 1},
        niblack_k={"widget_type": "FloatSlider", "min": -1.0, "max": 4.0, "step": 0.05},
        sauvola_k={"widget_type": "FloatSlider", "min": 0.0, "max": 1.0, "step": 0.01},
        bernsen_contrast={"widget_type": "FloatSlider", "min": 0.0, "max": 1.0, "step": 0.01},
        floor_k={"widget_type": "FloatSlider", "min": 0.0, "max": 10.0, "step": 0.1},
    )
    def threshold_controls(
        method: str = METHOD_NIBLACK,
        window_xy: int = 50,
        niblack_k: float = 3.0,
        sauvola_k: float = 0.2,
        bernsen_contrast: float = 0.15,
        floor_k: float = 3.0,
    ):
        _run_cadherin()

    def _run_cadherin(after=None):
        cadh_raw = state.get("cadh_raw")
        if cadh_raw is None:
            print("Load a stack with a 'substrate' channel first.", flush=True)
            call_all(after)
            return
        scale = state["scale"]
        method = threshold_controls["method"].value
        window_xy = int(threshold_controls["window_xy"].value)
        niblack_k = float(threshold_controls["niblack_k"].value)
        sauvola_k = float(threshold_controls["sauvola_k"].value)
        bernsen_contrast = float(threshold_controls["bernsen_contrast"].value)
        floor_k = float(threshold_controls["floor_k"].value)
        state.setdefault("run_settings", {})["cadherin"] = {
            "method": method, "window_xy": window_xy, "niblack_k": niblack_k,
            "sauvola_k": sauvola_k, "bernsen_contrast": bernsen_contrast,
            "floor_k": floor_k}

        @thread_worker
        def _run():
            vol = cadh_raw
            if method == METHOD_NIBLACK:
                mask = niblack_mask(vol, window_xy, niblack_k, scale)
            elif method == METHOD_SAUVOLA:
                mask = sauvola_mask(vol, window_xy, sauvola_k, scale)
            else:
                mask = bernsen_mask(vol, window_xy, bernsen_contrast, scale)
            floor = robust_floor_value(vol, floor_k)
            if floor is not None:
                mask = (mask.astype(bool) & (vol > floor)).astype(np.uint8)
            labels, sizes = label_and_sizes(mask)
            # Centroids/normals are NOT computed here — see `_ensure_centroids`.
            suggested = auto_min_size_otsu(sizes)
            return labels, sizes, suggested

        def _done(result):
            labels, sizes, suggested = result
            cs = state["cadh_state"]
            cs["labels"] = labels
            cs["sizes"] = sizes
            cs["centroids"] = None      # lazy, see `_ensure_centroids`
            cs["normals"] = None
            cs["label_ids"] = None
            cs["keep_mask"] = None
            no_mesh_notified["done"] = False
            slider_max = int(sizes.max()) if sizes.size else 5000
            filter_controls.min_size.max = slider_max
            filter_controls.min_size.value = max(suggested, MIN_SIZE_FLOOR)
            # Refresh the mask layer synchronously — the auto-called filter is
            # debounced 150 ms, but the "Run all" chain's next step (vessel mesh)
            # reads this layer immediately, so it must be current *now*.
            _layer = state.get("cadh_mask_layer")
            _m = _compute_mask(int(filter_controls.min_size.value),
                               float(filter_controls.band_um.value))
            if _layer is not None and _m is not None:
                _layer.data = pyramid_views(_m, state.get("levels", 1))
                _layer.visible = True
            _redraw(int(filter_controls.min_size.value))
            filter_controls()
            # Printed HERE, not when the worker returned. The threshold itself
            # is only part of the wait: `_compute_mask` + `pyramid_views` +
            # `_redraw` run after it on the GUI thread and take seconds on a
            # full stack, so a "done" printed on return arrives while the
            # viewer is still updating and cannot be trusted to mean the image
            # is final. Carries the voxel count and the elapsed time because
            # "is it still going?" is the actual question being asked.
            vox = int(_m.sum()) if _m is not None else 0
            print(
                f"Cadherin: {method} DONE in {time.perf_counter() - _t0:.1f} s "
                f"— {sizes.size:,} components, "
                f"{vox:,} voxels in the displayed mask. "
                f"Suggested min_size = {suggested} "
                f"(applied {int(filter_controls.min_size.value)}).",
                flush=True,
            )

        def _err(exc):
            print(f"Cadherin error: {exc}", flush=True)

        floor = robust_floor_value(cadh_raw, floor_k)
        floor_str = f"{floor:.1f}" if floor is not None else "off"
        print(f"Cadherin: {method} running (window_xy={window_xy}, "
              f"floor_k={floor_k} → floor={floor_str}) ...", flush=True)
        _t0 = time.perf_counter()
        w = _run()
        w.returned.connect(with_after(_done, after, "Cadherin"))
        w.errored.connect(with_after(_err, after, "Cadherin"))
        w.start()

    def _seed_floor():
        raw = state.get("cadh_raw")
        if raw is None:
            return
        threshold_controls.floor_k.value = min(
            suggest_floor_k(raw), threshold_controls.floor_k.max
        )
    state["seed_cadh_floor"] = _seed_floor

    # ---- Save / load the cadherin components (mask + labels for re-tuning) ----
    cadh_path_w = FileEdit(value="cadherin_labels.npz", label="cadh_path",
                           mode="w")
    save_cadh_w = PushButton(text="Save VE-Cadherin mask")
    load_cadh_w = PushButton(text="Load VE-Cadherin mask")

    def _save_cadh(*_):
        cs = state.get("cadh_state") or {}
        if cs.get("labels") is None:
            print("No cadherin components to save — run 'Apply threshold + "
                  "label' first.", flush=True)
            return
        path = str(cadh_path_w.value)
        np.savez_compressed(
            path,
            labels=np.asarray(cs["labels"]).astype(np.int32),
            sizes=np.asarray(cs["sizes"]),
            centroids=np.asarray(cs["centroids"]),
            normals=np.asarray(cs["normals"]),
            label_ids=np.asarray(cs["label_ids"]),
        )
        print(f"Saved VE-Cadherin components -> {path} "
              f"({np.asarray(cs['sizes']).size} components).", flush=True)

    def _load_cadh(*_, after=None, visible=True):
        path = str(cadh_path_w.value)

        @thread_worker
        def _run():
            d = np.load(path)
            return (d["labels"].astype(np.int32), d["sizes"], d["centroids"],
                    d["normals"], d["label_ids"])

        def _done(result):
            labels, sizes, centroids, normals, lids = result
            cadh_raw = state.get("cadh_raw")
            if cadh_raw is not None and labels.shape != cadh_raw.shape:
                print(f"Warning: loaded cadherin {labels.shape} != current "
                      f"{cadh_raw.shape}; load the matching image first.",
                      flush=True)
            cs = state["cadh_state"]
            cs["labels"] = labels
            cs["sizes"] = sizes
            cs["centroids"] = centroids
            cs["normals"] = normals
            cs["label_ids"] = lids
            cs["keep_mask"] = None
            no_mesh_notified["done"] = False
            suggested = auto_min_size_otsu(sizes)
            filter_controls.min_size.max = int(sizes.max()) if sizes.size else 5000
            filter_controls.min_size.value = max(suggested, MIN_SIZE_FLOOR)
            filter_controls()  # rebuild the "VE-Cadherin mask" layer
            if not visible and state.get("cadh_mask_layer") is not None:
                state["cadh_mask_layer"].visible = False
            print(f"Loaded VE-Cadherin components <- {path}: {sizes.size} "
                  f"components.", flush=True)

        def _err(exc):
            print(f"Load cadherin error: {exc}", flush=True)

        print(f"Loading VE-Cadherin components from {path} ...", flush=True)
        w = _run()
        w.returned.connect(with_after(_done, after, "Load cadherin"))
        w.errored.connect(with_after(_err, after, "Load cadherin"))
        w.start()

    save_cadh_w.changed.connect(_save_cadh)
    load_cadh_w.changed.connect(_load_cadh)
    state["load_cadh"] = _load_cadh

    # Chainable "Run all" step: force Sauvola, compute, save, then `after`.
    def _run_cadherin_step(after=None):
        threshold_controls["method"].value = METHOD_SAUVOLA
        _run_cadherin(after=lambda: call_all(_save_cadh, after,
                                             label="run_cadherin"))
    state["run_cadherin"] = _run_cadherin_step
    state.setdefault("save_path_widgets", []).append(
        (cadh_path_w, "cadherin_labels.npz"))

    return (Container(widgets=[threshold_controls, filter_controls,
                              cadh_path_w, save_cadh_w, load_cadh_w],
                     labels=False),
            hist.canvas)
