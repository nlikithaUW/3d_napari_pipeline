"""DAPI Otsu + size-filter tab."""

from __future__ import annotations

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container
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

    @magicgui(call_button="Compute DAPI Otsu + components")
    def otsu_controls():
        dapi_raw = state.get("dapi_raw")
        if dapi_raw is None:
            print("Load a stack with a 'nuclei' channel first.", flush=True)
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

        print("DAPI Otsu + labeling running ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.start()

    state["dapi_size_controls"] = size_controls
    return (Container(widgets=[otsu_controls, size_controls], labels=False),
            hist.canvas)
