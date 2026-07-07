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
from magicgui.widgets import Container
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
        cadh_raw = state.get("cadh_raw")
        if cadh_raw is None:
            print("Load a stack with a 'substrate' channel first.", flush=True)
            return
        scale = state["scale"]

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
            filter_controls()

        floor = robust_floor_value(cadh_raw, floor_k)
        floor_str = f"{floor:.1f}" if floor is not None else "off"
        print(f"Cadherin: {method} running (window_xy={window_xy}, "
              f"floor_k={floor_k} → floor={floor_str}) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.start()

    def _seed_floor():
        raw = state.get("cadh_raw")
        if raw is None:
            return
        threshold_controls.floor_k.value = min(
            suggest_floor_k(raw), threshold_controls.floor_k.max
        )
    state["seed_cadh_floor"] = _seed_floor

    return (Container(widgets=[threshold_controls, filter_controls], labels=False),
            hist.canvas)
