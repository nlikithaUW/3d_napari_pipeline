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
  - poisson (oriented points): fits a smooth implicit function to the cadherin
    sheet instead of interpolating it, so noise averages out rather than
    becoming surface (`pipeline/vessel_poisson.py`, `docs/ISSUES.md` I06). Three
    knobs, and the surface it returns is **open at the x/y edges of the volume
    by design** — that is the Neumann boundary, not a hole.
"""

from __future__ import annotations

import numpy as np
from magicgui import magicgui
from magicgui.widgets import (
    ComboBox, Container, FileEdit, FloatSlider, PushButton, Slider,
)
from napari.qt.threading import thread_worker

from ..cell_mesh import mesh_curvature, suggest_alpha_um, tet_circumradii
from .. import lumen_graph
from ..saved_check import curvature_notes
# The method constants and the build itself live in the core module so the
# headless `pipeline.batch_mesh` runner builds the *same* mesh from the same
# parameters. Re-exported here (the ComboBox choices below read them, and other
# modules import them from this one).
from ..vessel_mesh import (  # noqa: F401  (re-export)
    LIVE_PARAMS,
    check_alpha_vs_spacing,
    METHOD_ALPHA2D,
    METHOD_ALPHA3D,
    METHOD_ALPHA_POISSON,
    METHOD_POISSON,
    METHOD_POLAR,
    METHODS,
    ORIENTATIONS,
    build_mesh,
    mask_to_points_um,
    voxel_subsample,
)
from .widgets import Histogram, call_all, finite_surface, with_after

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
        poisson_grid_um={"widget_type": "FloatSlider", "min": 0.65, "max": 6.0, "step": 0.325},
        poisson_smooth_um={"widget_type": "FloatSlider", "min": 1.0, "max": 15.0, "step": 0.5},
        poisson_iso_pct={"widget_type": "FloatSlider", "min": 5.0, "max": 95.0, "step": 1.0},
        alpha_raster_um={"widget_type": "FloatSlider", "min": 0.325, "max": 6.0, "step": 0.325},
    )
    def controls(
        method: str = METHOD_ALPHA_POISSON,
        orientation: str = "both (OR)",
        n_bins: int = 180,
        r_smooth_bins: int = 3,
        alpha2d_um: float = 8.0,
        subsample_um: float = 5.0,
        alpha3d_um: float = 15.0,
        fill_lumen: bool = True,
        close_um: float = 1.0,
        min_component_um: float = 50.0,
        poisson_grid_um: float = 2.0,
        poisson_smooth_um: float = 5.0,
        poisson_iso_pct: float = 50.0,
        alpha_raster_um: float = 2.0,
    ):
        _run_build()

    # Every knob the build reads, in widget order. `LIVE_PARAMS` in the core is
    # the authority for which of them a given method actually uses; anything
    # else is inert for that method and is hidden rather than left on screen
    # looking operative. Same contract as the Clover tab.
    _BUILD_KNOBS = ("orientation", "n_bins", "r_smooth_bins", "alpha2d_um",
                    "subsample_um", "alpha3d_um", "fill_lumen", "close_um",
                    "min_component_um", "poisson_grid_um", "poisson_smooth_um",
                    "poisson_iso_pct", "alpha_raster_um")

    def _sync_method_visibility(*_):
        live = LIVE_PARAMS.get(controls["method"].value, _BUILD_KNOBS)
        for name in _BUILD_KNOBS:
            controls[name].visible = name in live
        # "Suggest alpha_um" only means anything where alpha3d_um is live.
        suggest_w.visible = "alpha3d_um" in live

    controls["method"].changed.connect(_sync_method_visibility)

    def _run_build(after=None):
        arr = _cadh_mask_array(state)
        if arr is None:
            print("Load a stack with a 'substrate' channel first.", flush=True)
            call_all(after)
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
        poisson_grid_um = float(controls["poisson_grid_um"].value)
        poisson_smooth_um = float(controls["poisson_smooth_um"].value)
        poisson_iso_pct = float(controls["poisson_iso_pct"].value)
        alpha_raster_um = float(controls["alpha_raster_um"].value)
        state.setdefault("run_settings", {})["vessel_mesh"] = {
            "method": method, "orientation": orientation, "alpha3d_um": alpha3d_um,
            "alpha2d_um": alpha2d_um, "subsample_um": subsample_um,
            "close_um": close_um, "min_component_um": min_component_um,
            "fill_lumen": fill_lumen, "n_bins": n_bins,
            "r_smooth_bins": r_smooth_bins,
            "poisson_grid_um": poisson_grid_um,
            "poisson_smooth_um": poisson_smooth_um,
            "poisson_iso_pct": poisson_iso_pct,
            "alpha_raster_um": alpha_raster_um}

        # The same alpha-vs-spacing rule `batch_mesh` applies. It lived only
        # in the headless path, so the rule was enforced there and silent here
        # — the drift this project keeps paying for.
        for n in check_alpha_vs_spacing(alpha3d_um, subsample_um):
            print(f"  {n}", flush=True)

        @thread_worker
        def _run():
            # One call into the core module — the identical build the headless
            # `pipeline.batch_mesh` runs, so the two paths cannot diverge.
            extras = {}
            out = build_mesh(
                arr, scale, method=method, orientation=orientation,
                n_bins=n_bins, r_smooth_bins=r_smooth_bins,
                alpha2d_um=alpha2d_um, subsample_um=subsample_um,
                alpha3d_um=alpha3d_um, fill_lumen=fill_lumen,
                close_um=close_um, min_component_um=min_component_um,
                poisson_grid_um=poisson_grid_um,
                poisson_smooth_um=poisson_smooth_um,
                poisson_iso_pct=poisson_iso_pct,
                alpha_raster_um=alpha_raster_um, extras=extras)
            return out, extras

        def _done(result):
            (kind, verts, faces, vals), extras = result
            normals = extras.get("normals")
            # `extras["notes"]` is build_mesh's own channel (the point-cloud
            # guard, P35); `extras["info"]["notes"]` is vessel_poisson's. Both
            # must be printed — the guard's soft warning went into a dict
            # nothing read when it was first added, which is exactly the silent
            # failure it exists to prevent.
            for n in extras.get("notes", []):
                print(f"  {n}", flush=True)
            if extras.get("n_points") is not None:
                print(f"  point cloud: {extras['n_points']:,} points, "
                      f"{extras.get('point_cloud_occupancy', 0):.1%} of the "
                      f"volume's subsample bins.", flush=True)
            for n in (extras.get("info") or {}).get("notes", []):
                print(f"  {n}", flush=True)
            # **The old layer is NOT removed until we know we can replace
            # it.** It used to be dropped before these checks, so a failed
            # build blanked the previous mesh and put nothing back — you lost
            # the result you had and got no picture of why.
            def _fail(msg, keep_verts=None):
                print(msg, flush=True)
                state["vessel_mesh_cache"] = {"verts": keep_verts,
                                              "faces": None}
                # **Drop the lumen field too.** It is set below, after these
                # early returns, so a failed build used to leave the PREVIOUS
                # build's field in `state` — and the centerline graph reads it.
                # In a batch that means one image's graph is built on another
                # image's lumen and saved under the failed image's name, which
                # is what happened on 2026-08-24 (P38).
                state["vessel_lumen_field"] = None

            if kind == "empty":
                # Do not blame the cadherin mask unconditionally: on the
                # alpha->poisson path `rasterize_alpha_solid` also returns an
                # empty solid when `alpha3d_um` is too small for the point
                # spacing, and the old message sent you to the wrong tab.
                if method in (METHOD_ALPHA3D, METHOD_ALPHA_POISSON):
                    n_pts = extras.get("n_points")
                    _fail(f"Nothing to mesh: the alpha complex came out empty "
                          f"at alpha3d_um={alpha3d_um:g}"
                          + (f" over {n_pts:,} points" if n_pts else "")
                          + f". Raise alpha3d_um — it must exceed the point "
                            f"spacing, and ~3x subsample_um "
                            f"(={3 * subsample_um:g}) is the working rule — or "
                            f"check the Cadherin tab's mask.")
                else:
                    _fail("The cadherin mask is empty — run the Cadherin tab "
                          "and tune its filter first.")
                return
            if verts is None or len(verts) == 0:
                _fail("No points/vertices to mesh.")
                return
            if faces is None or len(faces) == 0:
                if method in (METHOD_ALPHA3D, METHOD_ALPHA_POISSON):
                    hint = (f"increase alpha3d_um (now {alpha3d_um:g}; it must "
                            f"exceed the point spacing, ~3x subsample_um "
                            f"= {3 * subsample_um:g})")
                elif method == METHOD_POISSON:
                    hint = "lower poisson_iso_pct or raise poisson_smooth_um"
                else:
                    hint = "loosen the fill / orientation"
                _fail(f"{len(verts)} vertices but 0 faces — {hint}.", verts)
                return
            if LAYER_NAME in viewer.layers:
                viewer.layers.remove(LAYER_NAME)
            viewer.add_surface(
                finite_surface(verts, faces, vals), name=LAYER_NAME,
                colormap="viridis", blending="translucent", opacity=0.6,
            )
            state["vessel_mesh_cache"] = {"verts": verts, "faces": faces,
                                          "normals": normals}
            # The solved lumen field, for the 3D centerline graph below. Only
            # the Poisson family produces one; the interpolating methods leave
            # it absent and the graph button then says so rather than falling
            # back to something that looks like a centerline and is not.
            state["vessel_lumen_field"] = extras.get("lumen") or None
            print(f"Vessel mesh ({method}, orientation={orientation}): "
                  f"{len(verts)} verts, {len(faces)} faces.", flush=True)

        def _err(exc):
            print(f"Vessel mesh error: {exc}", flush=True)

        print(f"Building vessel mesh ({method}, "
              f"orientation={orientation}) ...", flush=True)
        w = _run()
        w.returned.connect(with_after(_done, after, "Vessel mesh"))
        w.errored.connect(with_after(_err, after, "Vessel mesh"))
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
        arrays = {"verts": np.asarray(verts), "faces": np.asarray(faces)}
        # An extra key only, the same additive pattern `z_um` uses in
        # `vessel_path_curvature.npz`: an older reader that names verts/faces is
        # unaffected, and a Poisson mesh keeps the normals it cannot re-derive
        # (it is open at the volume edges, so `vertex_outward_normals` has no
        # enclosed volume to read the global sign from).
        if cache.get("normals") is not None:
            arrays["normals"] = np.asarray(cache["normals"])
        np.savez(path, **arrays)
        print(f"Saved vessel mesh -> {path} "
              f"({len(verts)} verts, {len(faces)} faces).", flush=True)

    def _load_mesh(*_, after=None, visible=True):
        path = str(mesh_path_w.value)
        try:
            d = np.load(path)
            verts = d["verts"]
            faces = d["faces"]
            normals = d["normals"] if "normals" in d.files else None
            if len(verts) == 0 or len(faces) == 0:
                print("Loaded vessel mesh is empty.", flush=True)
                return
            # `values` (centroid distance, for the colormap) derives from verts.
            vals = np.linalg.norm(verts - verts.mean(0), axis=1)
            if LAYER_NAME in viewer.layers:
                viewer.layers.remove(LAYER_NAME)
            lyr = viewer.add_surface(
                finite_surface(verts, faces, vals), name=LAYER_NAME,
                colormap="viridis", blending="translucent", opacity=0.6,
            )
            lyr.visible = visible
            state["vessel_mesh_cache"] = {"verts": verts, "faces": faces,
                                          "normals": normals}
            print(f"Loaded vessel mesh <- {path} "
                  f"({len(verts)} verts, {len(faces)} faces).", flush=True)
        except Exception as exc:
            print(f"Load vessel mesh error: {exc}", flush=True)
        finally:
            call_all(after, label="Load vessel mesh")

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
            call_all(after)
            return
        # Carried through to `mesh_curvature`: on the Poisson surface the sign
        # of kappa2 is not recoverable from the mesh itself, so these are not
        # optional there (P01(B)).
        normals = cache.get("normals")
        radius = float(curv_radius_w.value)
        smooth_iter = int(curv_smooth_w.value)
        state.setdefault("run_settings", {})["surface_curvature"] = {
            "radius_um": radius, "smooth_iter": smooth_iter}

        @thread_worker
        def _run():
            return mesh_curvature(verts, faces, radius,
                                  smooth_iter=smooth_iter, normals=normals)

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

        def _err(exc):
            print(f"Curvature error: {exc}", flush=True)

        print(f"Computing vessel curvature (radius={radius} µm, "
              f"smooth={smooth_iter}) ...", flush=True)
        w = _run()
        w.returned.connect(with_after(_done, after, "Surface curvature"))
        w.errored.connect(with_after(_err, after, "Surface curvature"))
        w.start()

    curv_btn.changed.connect(_compute_curv)
    curv_measure_w.changed.connect(lambda *_: _show_curvature())

    # ---- In-plane (X-Y) vessel path curvature (centerline of the cadherin
    # footprint). Measures how sharply the vessel bends within the imaging
    # plane — the "straight vs curved vs branched" geometry — ignoring caliber
    # and the noisy Z axis. Built from the cadherin mask (the vessel mesh point
    # cloud is too rough to skeletonize cleanly). ----
    # ---- the 3D centerline graph (docs/DIRECTION.md §8B) --------------------
    # Replaces the in-plane centerline entirely. The old one skeletonized a max
    # projection of the cadherin mask; this one is built on the Poisson lumen,
    # so it has arc length, a real graph, per-point caliber and — the thing the
    # projection could never give — cross-sections cut perpendicular to the
    # axis.
    #
    # Every knob is a **ratio or a window**, not an absolute length. The old
    # `open_um` / `prune_um` were absolute and [P23] records them amputating
    # exactly the tapering, high-shear regions the program is about.
    spur_ratio_w = FloatSlider(value=lumen_graph.DEFAULT_SPUR_RATIO,
                               min=0.0, max=4.0, step=0.1, label="spur_ratio")
    crop_ratio_w = FloatSlider(value=lumen_graph.DEFAULT_CROP_RATIO,
                               min=0.0, max=6.0, step=0.25, label="crop_ratio")
    kappa_win_w = FloatSlider(value=lumen_graph.DEFAULT_KAPPA_WINDOW_UM,
                              min=5.0, max=120.0, step=5.0,
                              label="kappa_window_um")
    sec_every_w = FloatSlider(value=lumen_graph.DEFAULT_SECTION_EVERY_UM,
                              min=1.0, max=40.0, step=1.0,
                              label="section_every_um")
    sec_span_w = FloatSlider(value=lumen_graph.DEFAULT_SECTION_SPAN,
                             min=1.0, max=8.0, step=0.25, label="section_span")
    graph_color_w = ComboBox(
        choices=["caliber_um", "kappa_1_um", "kappa_xy_1_um", "lumen_area_um2",
                 "lumen_major_um", "lumen_minor_um", "lumen_aspect",
                 "dist_junction_um", "dist_crop_um", "censored"],
        value="caliber_um", label="colour by")
    graph_btn = PushButton(text="Compute 3D centerline graph (lumen)")

    GRAPH_LAYER = "Vessel centerline (3D graph)"
    GRAPH_PATH_LAYER = "Vessel centerline path"
    GRAPH_NODE_LAYER = "Centerline nodes"

    _COLOR_KEY = {"lumen_area_um2": "area_um2", "lumen_major_um": "major_um",
                  "lumen_minor_um": "minor_um", "lumen_aspect": "aspect"}

    def _show_graph(g, visible=True):
        """Draw the graph: the path as lines, the points coloured, nodes marked.

        **Drawing it as a path is new and is the point of [P27] closing.** The
        old centerline was a raster-order point set, so connecting consecutive
        points produced a scribble and it had to be drawn as loose dots. These
        points are ordered along each edge, so a `Shapes` polyline is exactly
        the vessel.
        """
        for name in (GRAPH_LAYER, GRAPH_PATH_LAYER, GRAPH_NODE_LAYER):
            if name in viewer.layers:
                viewer.layers.remove(name)
        pts = np.asarray(g.get("points_zyx_um", []))
        if not len(pts):
            return

        key = _COLOR_KEY.get(graph_color_w.value, graph_color_w.value)
        vals = np.asarray(g.get(key, np.zeros(len(pts))), dtype=np.float64)
        if vals.shape[0] != len(pts):
            vals = np.zeros(len(pts))
        finite = np.isfinite(vals)
        # NaN is a real answer here — sections are cut every few µm and the
        # gaps are deliberately not interpolated — so it is drawn at the low
        # end rather than dropped, which would make the path look broken.
        fill = float(np.nanmin(vals)) if finite.any() else 0.0
        shown = np.where(finite, vals, fill)

        eid = np.asarray(g.get("edge_id", np.zeros(len(pts), int)))
        paths = [pts[eid == e] for e in np.unique(eid)]
        paths = [q for q in paths if len(q) >= 2]
        if paths:
            lyr = viewer.add_shapes(paths, shape_type="path", name=GRAPH_PATH_LAYER,
                                    edge_color="white", edge_width=1.2,
                                    opacity=0.7)
            lyr.visible = visible

        pl = viewer.add_points(pts, name=GRAPH_LAYER, size=3.0,
                               face_color=shown, face_colormap="turbo",
                               border_width=0.0)
        if finite.any():
            lo = float(np.percentile(shown[finite], 2))
            hi = float(np.percentile(shown[finite], 98))
            try:
                pl.face_contrast_limits = (lo, max(hi, lo + 1e-6))
            except Exception:
                pass
        pl.visible = visible

        # Junctions and crop endpoints, which are what you check by eye: a
        # junction that is not a real Y, or a "terminus" that is really the
        # edge of the frame, is visible immediately and invisible in a table.
        nodes = np.asarray(g.get("nodes_zyx_um", []))
        if len(nodes):
            kind = np.asarray(g.get("node_kind", []))
            ek = np.asarray(g.get("node_endpoint_kind", []))
            col = np.where(kind == "junction", "red",
                           np.where(ek == "crop", "deepskyblue", "yellow"))
            nl = viewer.add_points(nodes, name=GRAPH_NODE_LAYER, size=9.0,
                                   face_color=list(col), border_width=0.0)
            nl.opacity = 0.85
            nl.visible = visible

    def _graph_params():
        return {"spur_ratio": float(spur_ratio_w.value),
                "crop_ratio": float(crop_ratio_w.value),
                "kappa_window_um": float(kappa_win_w.value),
                "section_every_um": float(sec_every_w.value),
                "section_span": float(sec_span_w.value)}

    def _compute_graph(*_, after=None):
        lum = state.get("vessel_lumen_field")
        if not lum or lum.get("chi") is None:
            print("Centerline graph: build the mesh first with a Poisson "
                  "method ('poisson (oriented points)' or 'alpha 3D solid -> "
                  "poisson'). The graph is built on the solved lumen, and the "
                  "interpolating methods never solve for one.", flush=True)
            call_all(after)
            return
        mask = _cadh_mask_array(state)
        scale = state.get("scale", (1.0, 1.0, 1.0))
        gp = _graph_params()
        state.setdefault("run_settings", {})["centerline_graph"] = dict(gp)
        dom = (np.asarray(np.asarray(mask).shape) * float(scale[0])
               if mask is not None else None)

        @thread_worker
        def _run():
            return lumen_graph.build_graph(
                lum["chi"], lum["iso_level"], lum["h_um"], lum["origin_um"],
                core=lum.get("core"), domain_shape_um=dom, **gp)

        def _done(g):
            state["vessel_graph"] = g
            for n in g.get("notes", []):
                print(f"  {n}", flush=True)
            if not len(g.get("points_zyx_um", [])):
                print("Centerline graph: nothing found.", flush=True)
                return
            _show_graph(g)
            def med(k):
                a = np.asarray(g.get(k, []), float)
                a = a[np.isfinite(a)]
                return float(np.median(a)) if a.size else float("nan")
            print(f"Centerline graph: {len(g['points_zyx_um'])} points on "
                  f"{len(g['edge_len_um'])} edge(s), {g['n_junctions']} "
                  f"junction(s), {g['n_endpoints']} endpoint(s) "
                  f"({g['n_crop_endpoints']} cropped, "
                  f"{g['n_terminus_endpoints']} terminating inside), "
                  f"{g['total_length_um']:.0f} µm of vessel.", flush=True)
            print(f"  inscribed radius median {med('caliber_um'):.1f} µm; "
                  f"cross-section {med('area_um2'):.0f} µm², "
                  f"{med('major_um'):.0f} x {med('minor_um'):.0f} µm, "
                  f"aspect {med('aspect'):.2f}.", flush=True)

        def _err(exc):
            print(f"Centerline graph error: {exc}", flush=True)

        print("Building the 3D centerline graph on the lumen ...", flush=True)
        w = _run()
        w.returned.connect(with_after(_done, after, "Centerline graph"))
        w.errored.connect(with_after(_err, after, "Centerline graph"))
        w.start()

    graph_color_w.changed.connect(
        lambda *_: (state.get("vessel_graph")
                    and _show_graph(state["vessel_graph"])))
    graph_btn.changed.connect(_compute_graph)

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
            cache = state.get("vessel_mesh_cache") or {"verts": None,
                                                       "faces": None,
                                                       "normals": None}
            cache["curvature"] = curv
            state["vessel_mesh_cache"] = cache
            verts = cache.get("verts")
            # Curvature is per-vertex, so it is only meaningful against the mesh
            # it was computed on. A length mismatch means this file belongs to a
            # different (or re-built) mesh: Measure silently drops every
            # curvature column to NaN in that case, so say so here instead.
            # Delegated to saved_check so this warning and the headless one are
            # the same string rather than two that drift apart (P12).
            if verts is not None:
                for n in curvature_notes(curv, len(verts)):
                    print(f"  {n}", flush=True)
            if verts is not None and cache.get("faces") is not None:
                _show_curvature()
                if not visible and CURV_LAYER in viewer.layers:
                    viewer.layers[CURV_LAYER].visible = False
            print(f"Loaded surface curvature <- {path}.", flush=True)
        except Exception as exc:
            print(f"Load curvature error: {exc}", flush=True)
        finally:
            call_all(after, label="Load curvature")

    path_curv_path_w = FileEdit(value="vessel_graph.npz",
                                label="graph_path", mode="w")
    save_path_curv_w = PushButton(text="Save centerline graph")
    load_path_curv_w = PushButton(text="Load centerline graph")

    def _save_path_curv(*_):
        g = state.get("vessel_graph")
        if not g or not len(g.get("points_zyx_um", [])):
            print("No centerline graph to save — click 'Compute 3D centerline "
                  "graph' first.", flush=True)
            return
        path = str(path_curv_path_w.value)
        # Arrays only. The notes and the nested info dicts are provenance and
        # belong in vessel_recon.json, which is where provenance already lives.
        np.savez_compressed(path, **{k: v for k, v in g.items()
                                     if isinstance(v, np.ndarray)})
        print(f"Saved centerline graph -> {path}: "
              f"{len(g['points_zyx_um'])} points, "
              f"{len(g['edge_len_um'])} edge(s).", flush=True)

    def _load_path_curv(*_, after=None, visible=True):
        path = str(path_curv_path_w.value)
        try:
            d = np.load(path, allow_pickle=False)
            g = {k: d[k] for k in d.files}
            if "points_zyx_um" not in g:
                print(f"Load centerline graph: {path} is an old in-plane "
                      f"`vessel_path_curvature.npz`, which this build no "
                      f"longer reads — the 2D centerline was replaced by the "
                      f"3D lumen graph (docs/DIRECTION.md §8B). Rebuild the "
                      f"mesh with a Poisson method and click 'Compute 3D "
                      f"centerline graph'.", flush=True)
                return
            for k in ("n_junctions", "n_endpoints", "n_crop_endpoints",
                      "n_terminus_endpoints", "total_length_um"):
                g.setdefault(k, 0)
            state["vessel_graph"] = g
            _show_graph(g, visible=visible)
            print(f"Loaded centerline graph <- {path}: "
                  f"{len(g['points_zyx_um'])} points, "
                  f"{len(g.get('edge_len_um', []))} edge(s).", flush=True)
        except Exception as exc:
            print(f"Load centerline graph error: {exc}", flush=True)
        finally:
            call_all(after, label="Load centerline graph")

    save_curv_w.changed.connect(_save_curv)
    load_curv_w.changed.connect(_load_curv)
    save_path_curv_w.changed.connect(_save_path_curv)
    load_path_curv_w.changed.connect(_load_path_curv)
    state["load_vessel_curvature"] = _load_curv
    state["load_vessel_path_curvature"] = _load_path_curv
    state.setdefault("save_path_widgets", []).extend([
        (curv_path_w, "vessel_curvature.npz"),
        (path_curv_path_w, "vessel_graph.npz"),
    ])

    # Chainable "Run all" steps: compute then save, then call `after`.
    state["run_vessel_mesh"] = lambda after=None: _run_build(
        after=lambda: call_all(_save_mesh, after, label="run_vessel_mesh"))
    state["run_vessel_curvature"] = lambda after=None: _compute_curv(
        after=lambda: call_all(_save_curv, after,
                               label="run_vessel_curvature"))
    state["run_vessel_path_curvature"] = lambda after=None: _compute_graph(
        after=lambda: call_all(_save_path_curv, after,
                               label="run_vessel_path_curvature"))

    # Apply the visibility rule once at build time (`suggest_w` exists by now).
    _sync_method_visibility()

    return (Container(widgets=[controls, suggest_w,
                              mesh_path_w, save_mesh_w, load_mesh_w,
                              curv_measure_w, curv_radius_w, curv_smooth_w,
                              curv_btn, curv_path_w, save_curv_w, load_curv_w,
                              spur_ratio_w, crop_ratio_w, kappa_win_w,
                              sec_every_w, sec_span_w, graph_color_w,
                              graph_btn,
                              path_curv_path_w, save_path_curv_w,
                              load_path_curv_w],
                     labels=False),
            hist.canvas)
