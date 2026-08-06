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
from magicgui.widgets import (
    ComboBox, Container, FileEdit, FloatSlider, PushButton, Slider,
)
from napari.qt.threading import thread_worker

from ..cell_mesh import (
    alpha_shape_faces,
    filter_faces_by_extent,
    mesh_curvature,
    suggest_alpha_um,
    tet_circumradii,
)
from ..centerline import vessel_path_curvature_xy
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
        _run_build()

    def _run_build(after=None):
        arr = _cadh_mask_array(state)
        if arr is None:
            print("Load a stack with a 'substrate' channel first.", flush=True)
            if after:
                after()
            return
        scale = state["scale"]
        # Read the parameter widgets by name via item access. Attribute access
        # (controls.orientation) collides with magicgui's Container.orientation
        # layout property and returns the layout string instead of the widget.
        method = controls["method"].value
        orientation = controls["orientation"].value
        n_bins = int(controls["n_bins"].value)
        r_smooth_bins = int(controls["r_smooth_bins"].value)
        alpha2d_um = float(controls["alpha2d_um"].value)
        subsample_um = float(controls["subsample_um"].value)
        alpha3d_um = float(controls["alpha3d_um"].value)
        fill_lumen = bool(controls["fill_lumen"].value)
        close_um = float(controls["close_um"].value)
        min_component_um = float(controls["min_component_um"].value)
        state.setdefault("run_settings", {})["vessel_mesh"] = {
            "method": method, "orientation": orientation, "alpha3d_um": alpha3d_um,
            "alpha2d_um": alpha2d_um, "subsample_um": subsample_um,
            "close_um": close_um, "min_component_um": min_component_um,
            "fill_lumen": fill_lumen, "n_bins": n_bins,
            "r_smooth_bins": r_smooth_bins}

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
                if after:
                    after()
                return
            if verts is None or len(verts) == 0:
                print("No points/vertices to mesh.", flush=True)
                state["vessel_mesh_cache"] = {"verts": None, "faces": None}
                if after:
                    after()
                return
            if faces is None or len(faces) == 0:
                hint = ("increase alpha3d_um" if method == METHOD_ALPHA3D
                        else "loosen the fill / orientation")
                print(f"{len(verts)} vertices but 0 faces — {hint}.", flush=True)
                state["vessel_mesh_cache"] = {"verts": verts, "faces": None}
                if after:
                    after()
                return
            viewer.add_surface(
                finite_surface(verts, faces, vals), name=LAYER_NAME,
                colormap="viridis", blending="translucent", opacity=0.6,
            )
            state["vessel_mesh_cache"] = {"verts": verts, "faces": faces}
            print(f"Vessel mesh ({method}, orientation={orientation}): "
                  f"{len(verts)} verts, {len(faces)} faces.", flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Vessel mesh error: {exc}", flush=True)
            if after:
                after()

        print(f"Building vessel mesh ({method}, "
              f"orientation={orientation}) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
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

    # ---- Save / load the built vessel mesh (verts + faces) ----
    mesh_path_w = FileEdit(value="vessel_mesh.npz", label="mesh_path", mode="w")
    save_mesh_w = PushButton(text="Save vessel mesh")
    load_mesh_w = PushButton(text="Load vessel mesh")

    def _save_mesh(*_):
        cache = state.get("vessel_mesh_cache") or {}
        verts = cache.get("verts")
        faces = cache.get("faces")
        if verts is None or faces is None or len(verts) == 0 or len(faces) == 0:
            print("No vessel mesh to save — build it first.", flush=True)
            return
        path = str(mesh_path_w.value)
        np.savez(path, verts=np.asarray(verts), faces=np.asarray(faces))
        print(f"Saved vessel mesh -> {path} "
              f"({len(verts)} verts, {len(faces)} faces).", flush=True)

    def _load_mesh(*_, after=None, visible=True):
        path = str(mesh_path_w.value)
        try:
            d = np.load(path)
            verts = d["verts"]
            faces = d["faces"]
        except Exception as exc:
            print(f"Load vessel mesh error: {exc}", flush=True)
            if after:
                after()
            return
        if len(verts) == 0 or len(faces) == 0:
            print("Loaded vessel mesh is empty.", flush=True)
            if after:
                after()
            return
        # `values` (centroid distance, for the colormap) is derived from verts.
        vals = np.linalg.norm(verts - verts.mean(0), axis=1)
        if LAYER_NAME in viewer.layers:
            viewer.layers.remove(LAYER_NAME)
        lyr = viewer.add_surface(
            finite_surface(verts, faces, vals), name=LAYER_NAME,
            colormap="viridis", blending="translucent", opacity=0.6,
        )
        lyr.visible = visible
        state["vessel_mesh_cache"] = {"verts": verts, "faces": faces}
        print(f"Loaded vessel mesh <- {path} "
              f"({len(verts)} verts, {len(faces)} faces).", flush=True)
        if after:
            after()

    save_mesh_w.changed.connect(_save_mesh)
    load_mesh_w.changed.connect(_load_mesh)
    state["load_vessel_mesh"] = _load_mesh
    state.setdefault("save_path_widgets", []).append(
        (mesh_path_w, "vessel_mesh.npz"))

    # ---- Per-vertex curvature (local quadric fit at a physical radius) ----
    CURV_CHOICES = {"mean (H)": "mean", "κ1 max": "k1",
                    "κ2 min": "k2", "gaussian (K)": "gaussian"}
    curv_measure_w = ComboBox(choices=list(CURV_CHOICES.keys()),
                              value="mean (H)", label="curvature")
    curv_radius_w = FloatSlider(value=5.0, min=1.0, max=50.0, step=0.5,
                                label="curvature_radius_um")
    curv_smooth_w = Slider(value=10, min=0, max=50,
                           label="curvature_smooth_iter")
    curv_btn = PushButton(text="Compute curvature")
    CURV_LAYER = "Vessel curvature"

    def _show_curvature():
        cache = state.get("vessel_mesh_cache") or {}
        curv = cache.get("curvature")
        verts = cache.get("verts")
        faces = cache.get("faces")
        if curv is None or verts is None or faces is None:
            return
        vals = np.asarray(curv[CURV_CHOICES[curv_measure_w.value]], dtype=float)
        finite = vals[np.isfinite(vals)]
        lo, hi = (np.percentile(finite, [2, 98]) if finite.size else (0.0, 1.0))
        # Interior/rank-deficient vertices are NaN and belong to no face, but
        # vispy still buffers every vertex value — replace NaN with the median.
        fill = float(np.median(finite)) if finite.size else 0.0
        vals = np.where(np.isfinite(vals), vals, fill)
        if CURV_LAYER in viewer.layers:
            viewer.layers.remove(CURV_LAYER)
        lyr = viewer.add_surface(
            finite_surface(verts, faces, vals), name=CURV_LAYER,
            colormap="turbo", blending="translucent", opacity=0.9,
        )
        try:
            lyr.contrast_limits = (float(lo), float(hi))
        except Exception:
            pass

    def _compute_curv(*_, after=None):
        cache = state.get("vessel_mesh_cache") or {}
        verts = cache.get("verts")
        faces = cache.get("faces")
        if verts is None or faces is None or len(verts) == 0 or len(faces) == 0:
            print("Compute curvature: build or load the vessel mesh first.",
                  flush=True)
            if after:
                after()
            return
        radius = float(curv_radius_w.value)
        smooth_iter = int(curv_smooth_w.value)
        state.setdefault("run_settings", {})["surface_curvature"] = {
            "radius_um": radius, "smooth_iter": smooth_iter}

        @thread_worker
        def _run():
            return mesh_curvature(verts, faces, radius, smooth_iter=smooth_iter)

        def _done(curv):
            state["vessel_mesh_cache"]["curvature"] = curv
            k1 = np.asarray(curv["k1"]); k1 = k1[np.isfinite(k1)]
            n_surf = int(k1.size)
            for k in ("mean", "k1", "k2", "gaussian"):
                v = np.asarray(curv[k])
                v = v[np.isfinite(v)]
                if v.size:
                    print(f"  {k}: median={np.median(v):.4f}  range "
                          f"[{v.min():.4f}, {v.max():.4f}] 1/µm", flush=True)
            pos = k1[k1 > 0.01]
            if pos.size:
                print(f"  implied caliber ~ {1.0/np.median(pos):.1f} µm "
                      f"(from median κ1)", flush=True)
            _show_curvature()
            print(f"Curvature computed (radius={radius} µm, smooth={smooth_iter}) "
                  f"on {n_surf} surface vertices.", flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Curvature error: {exc}", flush=True)
            if after:
                after()

        print(f"Computing vessel curvature (radius={radius} µm, "
              f"smooth={smooth_iter}) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    curv_btn.changed.connect(_compute_curv)
    curv_measure_w.changed.connect(lambda *_: _show_curvature())

    # ---- In-plane (X-Y) vessel path curvature (centerline of the cadherin
    # footprint). Measures how sharply the vessel bends within the imaging
    # plane — the "straight vs curved vs branched" geometry — ignoring caliber
    # and the noisy Z axis. Built from the cadherin mask (the vessel mesh point
    # cloud is too rough to skeletonize cleanly). ----
    path_smooth_w = FloatSlider(value=20.0, min=5.0, max=60.0, step=1.0,
                                label="path_smooth_um")
    path_btn = PushButton(text="Compute in-plane path curvature (X-Y)")
    PATH_LAYER = "Vessel path curvature (X-Y)"

    def _show_path(res, visible=True):
        """Overlay the centerline (points at the vessel's mean Z, coloured by κ)."""
        pts = np.asarray(res["pts_yx_um"])
        kap = np.asarray(res["kappa"])
        if len(pts) == 0:
            return
        zc = 0.0
        cache = state.get("vessel_mesh_cache") or {}
        vv = cache.get("verts")
        if vv is not None and len(vv):
            zc = float(np.mean(np.asarray(vv)[:, 0]))
        pts3 = np.column_stack([np.full(len(pts), zc), pts[:, 0], pts[:, 1]])
        if PATH_LAYER in viewer.layers:
            viewer.layers.remove(PATH_LAYER)
        lyr = viewer.add_points(
            pts3, name=PATH_LAYER, size=3.0, face_color=kap,
            face_colormap="turbo", border_width=0.0)
        hi = float(np.percentile(kap, 98)) if kap.size else 1.0
        try:
            lyr.face_contrast_limits = (0.0, max(hi, 1e-4))
        except Exception:
            pass
        lyr.visible = visible

    def _compute_path(*_, after=None):
        mask = _cadh_mask_array(state)
        if mask is None:
            print("Path curvature: compute or load the cadherin mask first "
                  "(Cadherin tab).", flush=True)
            if after:
                after()
            return
        scale = state.get("scale", (1.0, 1.0, 1.0))
        smooth_um = float(path_smooth_w.value)
        state.setdefault("run_settings", {})["path_curvature"] = {
            "smooth_um": smooth_um}

        @thread_worker
        def _run():
            m = np.asarray(mask)
            return vessel_path_curvature_xy(m > 0, scale, smooth_um=smooth_um)

        def _done(res):
            state["vessel_centerline"] = res
            pts = res["pts_yx_um"]
            kap = res["kappa"]
            nbr = len(res["branch_yx_um"])
            if len(pts) == 0:
                print("Path curvature: no centerline found.", flush=True)
                if after:
                    after()
                return
            _show_path(res)
            med = float(np.median(kap))
            print(f"In-plane path curvature: {len(pts)} centerline points, "
                  f"median κ = {med:.4f} 1/µm "
                  f"(radius ~ {1.0/med:.0f} µm)" if med > 0 else
                  f"In-plane path curvature: {len(pts)} centerline points, "
                  f"median κ = 0 (straight); {nbr} branch node(s).", flush=True)
            if med > 0:
                print(f"  {nbr} branch node(s) detected.", flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Path curvature error: {exc}", flush=True)
            if after:
                after()

        print("Computing in-plane vessel path curvature (X-Y) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    path_btn.changed.connect(_compute_path)

    # ---- Save / load the computed curvature so it need not be recomputed each
    # session (only the mesh geometry is saved with the mesh; curvature and the
    # centerline are separate). Surface curvature -> the 4 per-vertex arrays;
    # path curvature -> the centerline points, κ and branch nodes. ----
    curv_path_w = FileEdit(value="vessel_curvature.npz", label="curv_path",
                           mode="w")
    save_curv_w = PushButton(text="Save curvature")
    load_curv_w = PushButton(text="Load curvature")

    def _save_curv(*_):
        curv = (state.get("vessel_mesh_cache") or {}).get("curvature")
        if not curv:
            print("No curvature to save — click 'Compute curvature' first.",
                  flush=True)
            return
        path = str(curv_path_w.value)
        np.savez_compressed(path, **{k: np.asarray(curv[k])
                                     for k in ("mean", "k1", "k2", "gaussian")})
        print(f"Saved surface curvature -> {path}.", flush=True)

    def _load_curv(*_, after=None, visible=True):
        path = str(curv_path_w.value)
        try:
            d = np.load(path)
            curv = {k: d[k] for k in ("mean", "k1", "k2", "gaussian")}
        except Exception as exc:
            print(f"Load curvature error: {exc}", flush=True)
            if after:
                after()
            return
        cache = state.get("vessel_mesh_cache") or {"verts": None, "faces": None}
        cache["curvature"] = curv
        state["vessel_mesh_cache"] = cache
        if cache.get("verts") is not None and cache.get("faces") is not None:
            _show_curvature()
            if not visible and CURV_LAYER in viewer.layers:
                viewer.layers[CURV_LAYER].visible = False
        print(f"Loaded surface curvature <- {path}.", flush=True)
        if after:
            after()

    path_curv_path_w = FileEdit(value="vessel_path_curvature.npz",
                                label="path_curv_path", mode="w")
    save_path_curv_w = PushButton(text="Save path curvature")
    load_path_curv_w = PushButton(text="Load path curvature")

    def _save_path_curv(*_):
        cl = state.get("vessel_centerline")
        if not cl or not len(cl.get("pts_yx_um", [])):
            print("No path curvature to save — click 'Compute in-plane path "
                  "curvature' first.", flush=True)
            return
        path = str(path_curv_path_w.value)
        np.savez_compressed(path, pts_yx_um=np.asarray(cl["pts_yx_um"]),
                            kappa=np.asarray(cl["kappa"]),
                            branch_yx_um=np.asarray(cl["branch_yx_um"]))
        print(f"Saved path curvature -> {path}.", flush=True)

    def _load_path_curv(*_, after=None, visible=True):
        path = str(path_curv_path_w.value)
        try:
            d = np.load(path)
            cl = {"pts_yx_um": d["pts_yx_um"], "kappa": d["kappa"],
                  "branch_yx_um": d["branch_yx_um"]}
        except Exception as exc:
            print(f"Load path curvature error: {exc}", flush=True)
            if after:
                after()
            return
        state["vessel_centerline"] = cl
        _show_path(cl, visible=visible)
        print(f"Loaded path curvature <- {path}: {len(cl['pts_yx_um'])} "
              f"centerline points, {len(cl['branch_yx_um'])} branch node(s).",
              flush=True)
        if after:
            after()

    save_curv_w.changed.connect(_save_curv)
    load_curv_w.changed.connect(_load_curv)
    save_path_curv_w.changed.connect(_save_path_curv)
    load_path_curv_w.changed.connect(_load_path_curv)
    state["load_vessel_curvature"] = _load_curv
    state["load_vessel_path_curvature"] = _load_path_curv
    state.setdefault("save_path_widgets", []).extend([
        (curv_path_w, "vessel_curvature.npz"),
        (path_curv_path_w, "vessel_path_curvature.npz"),
    ])

    # Chainable "Run all" steps: compute then save, then call `after`.
    state["run_vessel_mesh"] = lambda after=None: _run_build(
        after=(lambda: (_save_mesh(), after() if after else None)))
    state["run_vessel_curvature"] = lambda after=None: _compute_curv(
        after=(lambda: (_save_curv(), after() if after else None)))
    state["run_vessel_path_curvature"] = lambda after=None: _compute_path(
        after=(lambda: (_save_path_curv(), after() if after else None)))

    return (Container(widgets=[controls, suggest_w,
                              mesh_path_w, save_mesh_w, load_mesh_w,
                              curv_measure_w, curv_radius_w, curv_smooth_w,
                              curv_btn, curv_path_w, save_curv_w, load_curv_w,
                              path_smooth_w, path_btn,
                              path_curv_path_w, save_path_curv_w,
                              load_path_curv_w],
                     labels=False),
            hist.canvas)
