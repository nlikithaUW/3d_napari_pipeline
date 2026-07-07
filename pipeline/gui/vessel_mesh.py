"""Vessel mesh tab — compare three cadherin-closing strategies.

The wall mask is NOT filtered here: it consumes the live filtered mask from the
Cadherin tab (the "VE-Cadherin mask" layer), so the filter lives in one place
and it is explicit which cadherin is being meshed.

`method` picks the closing strategy and `orientation` picks the per-slice
slicing (ignored by 3D alpha):
  - polar (per-slice): cheap, star-shaped lumen; collapses at the Y-branch.
  - alpha 2D (per-slice): one `alpha2d_um` knob; branch → two loops.
  - alpha 3D (points): reuses the tested `cell_mesh` alpha shape on the
    subsampled cadherin point cloud, with a Suggest button + histogram.
"""

from __future__ import annotations

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container, PushButton
from napari.qt.threading import thread_worker

from ..cell_mesh import (
    alpha_shape_faces,
    filter_faces_by_extent,
    suggest_alpha_um,
    tet_circumradii,
)
from ..vessel_mesh import (
    alpha_fill_volume,
    close_volume,
    combine,
    mask_to_points_um,
    mesh_from_volume,
    polar_fill_volume,
    voxel_subsample,
)
from .widgets import Histogram, finite_surface

METHOD_POLAR = "polar (per-slice)"
METHOD_ALPHA2D = "alpha 2D (per-slice)"
METHOD_ALPHA3D = "alpha 3D (points)"
METHODS = [METHOD_POLAR, METHOD_ALPHA2D, METHOD_ALPHA3D]
ORIENTATIONS = ["xz", "yz", "both (AND)", "both (OR)"]

LAYER_NAME = "Vessel mesh"


def _cadh_mask_array(state):
    """Full-res cadherin mask array from the Cadherin tab's layer, or None.

    Returns the level-0 array of the (possibly multiscale) mask layer without
    copying — the `> 0` bool conversion is deferred to the worker thread.
    """
    layer = state.get("cadh_mask_layer")
    if layer is None:
        return None
    # `layer.data` is a `MultiScaleData` (a Sequence, NOT a list/tuple) when
    # the layer is multiscale; `np.asarray()` on it silently collapses to the
    # coarsest level. Use the public `.multiscale` flag to grab full-res.
    data = layer.data
    return data[0] if layer.multiscale else data


def _per_slice(fill_fn, mask, orientation) -> np.ndarray:
    """Run a per-slice fill over the requested orientation(s) and combine."""
    if orientation in ("xz", "yz"):
        return fill_fn(mask, orientation)
    a = fill_fn(mask, "xz")
    b = fill_fn(mask, "yz")
    mode = "and" if "AND" in orientation else "or"
    return combine(a, b, mode)


def build(state, viewer) -> Container:
    hist = Histogram("Cadherin tet circumradii (alpha_um)",
                     "log10(circumradius µm)", "darkorange")
    state.setdefault("hists", []).append(hist)

    @magicgui(
        call_button="Build vessel mesh",
        method={"widget_type": "ComboBox", "choices": METHODS},
        orientation={"widget_type": "ComboBox", "choices": ORIENTATIONS},
        n_bins={"widget_type": "Slider", "min": 16, "max": 720, "step": 1},
        r_smooth_bins={"widget_type": "Slider", "min": 0, "max": 50, "step": 1},
        alpha2d_um={"widget_type": "FloatSlider", "min": 0.5, "max": 30.0, "step": 0.5},
        subsample_um={"widget_type": "FloatSlider", "min": 0.0, "max": 5.0, "step": 0.25},
        alpha3d_um={"widget_type": "FloatSlider", "min": 1.0, "max": 80.0, "step": 0.5},
        close_um={"widget_type": "FloatSlider", "min": 0.0, "max": 5.0, "step": 0.25},
        min_component_um={"widget_type": "FloatSlider", "min": 0.0, "max": 500.0, "step": 5.0},
    )
    def controls(
        method: str = METHOD_ALPHA3D,
        orientation: str = "both (OR)",
        n_bins: int = 180,
        r_smooth_bins: int = 3,
        alpha2d_um: float = 8.0,
        subsample_um: float = 1.0,
        alpha3d_um: float = 10.0,
        fill_lumen: bool = True,
        close_um: float = 1.0,
        min_component_um: float = 50.0,
    ):
        arr = _cadh_mask_array(state)
        if arr is None:
            print("Load a stack with a 'substrate' channel first.", flush=True)
            return
        scale = state["scale"]

        @thread_worker
        def _run():
            mask = np.asarray(arr) > 0
            if not mask.any():
                return ("empty", None, None, None)
            if method == METHOD_POLAR:
                filled = _per_slice(
                    lambda m, ax: polar_fill_volume(
                        m, ax, scale, int(n_bins), int(r_smooth_bins)),
                    mask, orientation,
                )
                filled = close_volume(filled, scale[0], float(close_um))
                verts, faces, vals = mesh_from_volume(filled, scale)
                faces = filter_faces_by_extent(verts, faces, float(min_component_um))
                return ("mesh", verts, faces, vals)
            if method == METHOD_ALPHA2D:
                filled = _per_slice(
                    lambda m, ax: alpha_fill_volume(
                        m, ax, scale, float(alpha2d_um), fill_lumen=fill_lumen),
                    mask, orientation,
                )
                filled = close_volume(filled, scale[0], float(close_um))
                verts, faces, vals = mesh_from_volume(filled, scale)
                faces = filter_faces_by_extent(verts, faces, float(min_component_um))
                return ("mesh", verts, faces, vals)
            # 3D alpha on the subsampled cadherin point cloud.
            pts = voxel_subsample(mask_to_points_um(mask, scale), float(subsample_um))
            faces = alpha_shape_faces(pts, float(alpha3d_um))
            faces = filter_faces_by_extent(pts, faces, float(min_component_um))
            vals = (np.linalg.norm(pts - pts.mean(0), axis=1)
                    if len(pts) else np.zeros((0,)))
            return ("points", pts, faces, vals)

        def _done(result):
            kind, verts, faces, vals = result
            if LAYER_NAME in viewer.layers:
                viewer.layers.remove(LAYER_NAME)
            if kind == "empty":
                print("The cadherin mask is empty — run the Cadherin tab and "
                      "tune its filter first.", flush=True)
                state["vessel_mesh_cache"] = {"verts": None, "faces": None}
                return
            if verts is None or len(verts) == 0:
                print("No points/vertices to mesh.", flush=True)
                state["vessel_mesh_cache"] = {"verts": None, "faces": None}
                return
            if faces is None or len(faces) == 0:
                hint = ("increase alpha3d_um" if method == METHOD_ALPHA3D
                        else "loosen the fill / orientation")
                print(f"{len(verts)} vertices but 0 faces — {hint}.", flush=True)
                state["vessel_mesh_cache"] = {"verts": verts, "faces": None}
                return
            viewer.add_surface(
                finite_surface(verts, faces, vals), name=LAYER_NAME,
                colormap="viridis", blending="translucent", opacity=0.6,
            )
            state["vessel_mesh_cache"] = {"verts": verts, "faces": faces}
            print(f"Vessel mesh ({method}, orientation={orientation}): "
                  f"{len(verts)} verts, {len(faces)} faces.", flush=True)

        print(f"Building vessel mesh ({method}, "
              f"orientation={orientation}) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.start()

    suggest_w = PushButton(text="Suggest alpha_um (3D)")

    def _suggest(*_):
        arr = _cadh_mask_array(state)
        if arr is None:
            print("Load a stack with a 'substrate' channel first.", flush=True)
            return
        scale = state["scale"]
        subsample_um = float(controls.subsample_um.value)

        @thread_worker
        def _run():
            mask = np.asarray(arr) > 0
            pts = voxel_subsample(mask_to_points_um(mask, scale), subsample_um)
            R = tet_circumradii(pts)
            return pts, R, suggest_alpha_um(R)

        def _done(result):
            pts, R, suggested = result
            if R.size == 0:
                print(f"Only {len(pts)} points — need ≥4 to tetrahedralize "
                      "(raise subsample_um or loosen the filter).", flush=True)
                hist.clear()
                return
            finite = R[np.isfinite(R)]
            print(f"alpha histogram: {R.size} tets, R range "
                  f"{finite.min():.2f}–{finite.max():.2f} µm, "
                  f"suggested alpha_um = {suggested:.2f}.", flush=True)
            R_plot = R[np.isfinite(R) & (R > 0)]
            if R_plot.size:
                cap = np.percentile(R_plot, 99.0)
                R_plot = R_plot[R_plot <= cap]
            hist.update(R_plot, threshold=suggested if suggested > 0 else None,
                        threshold_label=f"alpha_um={suggested:.2f}")
            if suggested > 0:
                if suggested > controls.alpha3d_um.max:
                    controls.alpha3d_um.max = float(suggested) * 2.0
                controls.alpha3d_um.value = suggested

        print("Computing cadherin point circumradii ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.start()

    suggest_w.changed.connect(_suggest)

    return (Container(widgets=[controls, suggest_w], labels=False), hist.canvas)
