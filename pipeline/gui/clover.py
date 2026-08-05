"""Clover Otsu + mesh-gated size-filter tab."""

from __future__ import annotations

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container, PushButton
from napari.qt.threading import thread_worker

from ..cell_mesh import signed_distance_to_mesh_vertices
from ..denoise import auto_min_size_otsu, filter_by_size, label_and_sizes, otsu_mask
from .widgets import debounce, pyramid_views


def build(state, viewer) -> Container:
    compute_w = PushButton(text="Compute Otsu")

    @magicgui(
        auto_call=True,
        min_size={"widget_type": "Slider", "min": 0, "max": 5000},
        gate_by_mesh={"widget_type": "CheckBox"},
        d_max_um={"widget_type": "FloatSlider", "min": 0.25, "max": 30.0, "step": 0.25},
    )
    def size_controls(
        min_size: int = 0,
        gate_by_mesh: bool = False,
        d_max_um: float = 5.0,
    ):
        _apply(min_size, gate_by_mesh, d_max_um)

    @debounce(150)
    def _apply(min_size, gate_by_mesh, d_max_um):
        labels = state["clover_state"]["labels"]
        sizes = state["clover_state"]["sizes"]
        layer = state.get("clover_mask_layer")
        if labels is None or layer is None:
            return
        method = state["clover_state"]["method"] or ""
        mask = filter_by_size(labels, sizes, min_size)

        if gate_by_mesh:
            mc = state["mesh_cache"]
            if mc["offset_centroids"] is None or len(mc["offset_centroids"]) == 0:
                print("gate_by_mesh requested but no cell mesh built yet — "
                      "skipping gate.", flush=True)
            else:
                scale = np.asarray(state["scale"], dtype=np.float64)
                idx = np.argwhere(mask.astype(bool))  # (N, 3) zyx voxel indices
                if len(idx):
                    query_um = idx.astype(np.float64) * scale[None, :]
                    sd = signed_distance_to_mesh_vertices(
                        query_um, mc["offset_centroids"], mc["outward"],
                        tree=mc.get("kdtree"),
                    )
                    keep = (sd > 0.0) & (sd <= float(d_max_um))
                    drop_idx = idx[~keep]
                    if len(drop_idx):
                        mask[drop_idx[:, 0], drop_idx[:, 1], drop_idx[:, 2]] = 0
                    print(
                        f"Clover gate: kept {int(keep.sum())}/{len(idx)} "
                        f"voxels in 0 < d ≤ {d_max_um:.2f} µm "
                        f"(inward_offset_um={mc['inward_offset_um']:.2f}).",
                        flush=True,
                    )

        layer.name = f"Clover filter ({method})"
        layer.data = pyramid_views(mask, state.get("levels", 1))
        layer.visible = True

    def _compute(*_, after=None):
        clover_raw = state.get("clover_raw")
        if clover_raw is None:
            print("Load a stack with a 'clover' channel first.", flush=True)
            if after:
                after()
            return
        method = "Otsu"

        @thread_worker
        def _run():
            mask = otsu_mask(clover_raw)
            labels, sizes = label_and_sizes(mask)
            suggested = auto_min_size_otsu(sizes)
            return labels, sizes, suggested

        def _done(result):
            labels, sizes, suggested = result
            state["clover_state"]["labels"] = labels
            state["clover_state"]["sizes"] = sizes
            state["clover_state"]["method"] = method
            print(
                f"Clover {method}: {sizes.size} components. "
                f"Suggested min_size = {suggested}.", flush=True,
            )
            slider_max = int(sizes.max()) if sizes.size else 5000
            size_controls.min_size.max = slider_max
            size_controls.min_size.value = suggested
            # Refresh the mask layer synchronously — size_controls' _apply is
            # debounced 150 ms, but RF predict's clover gate reads this layer
            # right after in "Run all", so it must be current now.
            _layer = state.get("clover_mask_layer")
            if _layer is not None:
                _mask = filter_by_size(labels, sizes,
                                       int(size_controls.min_size.value))
                _layer.name = f"Clover filter ({method})"
                _layer.data = pyramid_views(_mask, state.get("levels", 1))
                _layer.visible = True
            size_controls()
            if after:
                after()

        def _err(exc):
            print(f"Clover Otsu error: {exc}", flush=True)
            if after:
                after()

        print(f"Clover {method} running ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    compute_w.changed.connect(_compute)
    state["run_clover"] = _compute
    return Container(widgets=[compute_w, size_controls], labels=True)
