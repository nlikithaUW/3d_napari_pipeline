"""Clover threshold + mesh-gated size-filter tab.

The threshold was a single global Otsu over the whole stack. Clover attenuates
with depth, so one scalar necessarily over-accepts shallow and under-accepts
deep — the root cause of `docs/ISSUES.md` P03, confirmed optical rather than
biological on 2026-08-24. The local methods are selectable here (I07).

**The default is Niblack at k = 0.5, changed from global Otsu on 2026-08-24.**
`batch.py` runs every tab at its widget defaults and has no CLI for them, so
this default *is* the headless setting, and changing it invalidates the
`20260812_011*` adhesion CSVs. It was changed deliberately, with the
measurement in hand and on the user's decision — not drifted into. The numbers
are in I07; the two that carry it:

- **acceptance by depth** (shallow third vs deep third) goes from Otsu's
  **3.06** to **0.91**, so the gate stops discarding the deep stack;
- **recall** of the known adhesion voxels is **98.85%**, against 36.57% at the
  cadherin default of k = 3.0.

**This mask is a gate on the RF, not a segmentation.** Its job is to discard
only what is certainly not an adhesion — the RF decides the rest — so it is
tuned for recall and the retained volume is merely the price (~24% here,
against Otsu's 0.97%, i.e. ~25x more voxels through `rf.predict_proba`). That
is the opposite of how a segmentation threshold is tuned, and it is why the
Cadherin tab's k is the wrong number to copy.
"""

from __future__ import annotations

import time

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container, PushButton
from napari.qt.threading import thread_worker

from ..cell_mesh import signed_distance_to_mesh_vertices
from ..denoise import (
    MASK_LIVE_PARAMS,
    MASK_METHODS,
    METHOD_NIBLACK,
    build_mask,
    filter_by_size,
)
from .widgets import call_all, debounce, pyramid_views, with_after


def build(state, viewer) -> Container:
    compute_w = PushButton(text="Compute threshold")

    @magicgui(
        auto_call=True,
        method={"widget_type": "ComboBox", "choices": MASK_METHODS},
        window_xy={"widget_type": "Slider", "min": 5, "max": 101, "step": 1},
        niblack_k={"widget_type": "FloatSlider", "min": -1.0, "max": 4.0,
                   "step": 0.05},
        sauvola_k={"widget_type": "FloatSlider", "min": 0.0, "max": 1.0,
                   "step": 0.01},
        bernsen_contrast={"widget_type": "FloatSlider", "min": 0.0, "max": 1.0,
                          "step": 0.01},
    )
    def threshold_controls(
        method: str = METHOD_NIBLACK,
        window_xy: int = 50,
        niblack_k: float = 0.5,
        sauvola_k: float = 0.2,
        bernsen_contrast: float = 0.15,
    ):
        _sync_visibility(method)

    def _sync_visibility(method):
        # Hide the knobs this method does not read, so an inert slider cannot
        # look like it is doing something. MASK_LIVE_PARAMS is the authority.
        live = MASK_LIVE_PARAMS.get(method, ())
        for name in ("window_xy", "niblack_k", "sauvola_k", "bernsen_contrast"):
            threshold_controls[name].visible = name in live

    def threshold_params():
        # The tab's current settings, as build_mask kwargs. adhesion_rf's gate
        # fallback reads this too, so the mask the user tuned here and the mask
        # the RF gates on cannot be built from different knobs.
        return {
            "method": threshold_controls["method"].value,
            "window_xy": int(threshold_controls["window_xy"].value),
            "niblack_k": float(threshold_controls["niblack_k"].value),
            "sauvola_k": float(threshold_controls["sauvola_k"].value),
            "bernsen_contrast": float(
                threshold_controls["bernsen_contrast"].value),
        }

    state["clover_threshold_params"] = threshold_params
    _sync_visibility(METHOD_NIBLACK)

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
            call_all(after)
            return
        params = threshold_params()
        method = params["method"]
        scale = state["scale"]
        # Record the inert knobs as inert rather than logging them as if they
        # had an effect — same contract as vessel_mesh's LIVE_PARAMS.
        live = MASK_LIVE_PARAMS.get(method, ())
        state.setdefault("run_settings", {})["clover"] = {
            **{k: v for k, v in params.items() if k == "method" or k in live},
            "inert": [k for k in params if k != "method" and k not in live],
        }

        @thread_worker
        def _run():
            labels, sizes, suggested, _ = build_mask(
                clover_raw, scale, **params)
            return labels, sizes, suggested

        def _done(result):
            labels, sizes, suggested = result
            state["clover_state"]["labels"] = labels
            state["clover_state"]["sizes"] = sizes
            state["clover_state"]["method"] = method
            slider_max = int(sizes.max()) if sizes.size else 5000
            size_controls.min_size.max = slider_max
            size_controls.min_size.value = suggested
            # Refresh the mask layer synchronously — size_controls' _apply is
            # debounced 150 ms, but RF predict's clover gate reads this layer
            # right after in "Run all", so it must be current now.
            _layer = state.get("clover_mask_layer")
            _mask = None
            if _layer is not None:
                _mask = filter_by_size(labels, sizes,
                                       int(size_controls.min_size.value))
                _layer.name = f"Clover filter ({method})"
                _layer.data = pyramid_views(_mask, state.get("levels", 1))
                _layer.visible = True
            size_controls()
            # Printed HERE, not when the worker returned. The threshold is only
            # part of the wait: the size filter and `pyramid_views` run after it
            # on the GUI thread, so a "done" printed on return arrives while the
            # viewer is still updating and cannot be trusted to mean the image
            # is final. Carries elapsed time and the displayed voxel count
            # because "is it still going?" is the actual question, and because
            # this mask GATES THE RF — how much it retains is the number that
            # matters, not the component count (I07: a permissive gate
            # fragments heavily, 12.2M components at the default k).
            vox = int(_mask.sum()) if _mask is not None else 0
            pct = (100.0 * vox / _mask.size) if _mask is not None else 0.0
            print(
                f"Clover {method}: DONE in {time.perf_counter() - _t0:.1f} s "
                f"— {sizes.size:,} components, {vox:,} voxels "
                f"({pct:.2f}% of the volume) in the gate. "
                f"Suggested min_size = {suggested} "
                f"(applied {int(size_controls.min_size.value)}).",
                flush=True,
            )

        def _err(exc):
            print(f"Clover threshold error: {exc}", flush=True)

        print(f"Clover {method} running ...", flush=True)
        _t0 = time.perf_counter()
        w = _run()
        w.returned.connect(with_after(_done, after, "Clover threshold"))
        w.errored.connect(with_after(_err, after, "Clover threshold"))
        w.start()

    compute_w.changed.connect(_compute)
    state["run_clover"] = _compute
    return Container(widgets=[compute_w, threshold_controls, size_controls],
                     labels=True)
