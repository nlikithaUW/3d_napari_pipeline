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
from .widgets import Histogram, debounce, pyramid_views

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

    def _compute_mask(min_size: int, band_um: float):
        cs = state["cadh_state"]
        labels = cs["labels"]
        sizes = cs["sizes"]
        centroids = cs["centroids"]
        lids = cs["label_ids"]
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
        use_band = (band_um > 0 and kdtree is not None
                    and offset_centroids is not None and outward is not None
                    and len(offset_centroids) > 0 and len(centroids) > 0)

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
            keep_mask = size_pass[lids - 1] if lids.size else np.zeros(0, dtype=bool)

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
            if after:
                after()
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
            centroids, normals, lids = label_centroids_normals(
                labels, scale, min_voxels=1,
            )
            suggested = auto_min_size_otsu(sizes)
            return labels, sizes, centroids, normals, lids, suggested

        def _done(result):
            labels, sizes, centroids, normals, lids, suggested = result
            cs = state["cadh_state"]
            cs["labels"] = labels
            cs["sizes"] = sizes
            cs["centroids"] = centroids
            cs["normals"] = normals
            cs["label_ids"] = lids
            cs["keep_mask"] = None
            no_mesh_notified["done"] = False
            print(
                f"Cadherin: {method} done. {sizes.size} components, "
                f"{len(centroids)} with centroids. Suggested min_size = {suggested}.",
                flush=True,
            )
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
            if after:
                after()

        def _err(exc):
            print(f"Cadherin error: {exc}", flush=True)
            if after:
                after()

        floor = robust_floor_value(cadh_raw, floor_k)
        floor_str = f"{floor:.1f}" if floor is not None else "off"
        print(f"Cadherin: {method} running (window_xy={window_xy}, "
              f"floor_k={floor_k} → floor={floor_str}) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
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
            if after:
                after()

        def _err(exc):
            print(f"Load cadherin error: {exc}", flush=True)
            if after:
                after()

        print(f"Loading VE-Cadherin components from {path} ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    save_cadh_w.changed.connect(_save_cadh)
    load_cadh_w.changed.connect(_load_cadh)
    state["load_cadh"] = _load_cadh

    # Chainable "Run all" step: force Sauvola, compute, save, then `after`.
    def _run_cadherin_step(after=None):
        threshold_controls["method"].value = METHOD_SAUVOLA
        _run_cadherin(after=(lambda: (_save_cadh(), after() if after else None)))
    state["run_cadherin"] = _run_cadherin_step
    state.setdefault("save_path_widgets", []).append(
        (cadh_path_w, "cadherin_labels.npz"))

    return (Container(widgets=[threshold_controls, filter_controls,
                              cadh_path_w, save_cadh_w, load_cadh_w],
                     labels=False),
            hist.canvas)
