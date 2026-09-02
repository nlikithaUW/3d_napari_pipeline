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
    CheckBox, ComboBox, Container, FileEdit, FloatSlider, PushButton, Slider,
)
from napari.qt.threading import thread_worker

from ..cell_mesh import mesh_curvature, suggest_alpha_um, tet_circumradii
from .. import lumen_graph
from .. import shear
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

    # Every knob carries a tooltip saying what it does **and** what it costs
    # when it is wrong. A slider with only a name on it is a slider nobody can
    # set on purpose, and these are not self-explanatory: half of them are inert
    # for any given method, and several have a documented failure mode rather
    # than a smooth trade-off.
    @magicgui(
        call_button="Build vessel mesh",
        method={"widget_type": "ComboBox", "choices": METHODS,
                "tooltip": "How the gaps in the cadherin wall are closed into "
                           "a surface. The three interpolating methods pass "
                           "through the mask; 'poisson' fits a smooth field to "
                           "it and is the only family that also solves for a "
                           "lumen volume, which the centerline graph and every "
                           "shear number below are built on."},
        orientation={"widget_type": "ComboBox", "choices": ORIENTATIONS,
                     "tooltip": "Which plane the per-slice methods work in. "
                                "'both (OR)' unions the two. Ignored entirely "
                                "by the 3D methods."},
        n_bins={"widget_type": "Slider", "min": 16, "max": 720, "step": 1,
                "tooltip": "polar only: how many angular bins the wall radius "
                           "is sampled into per slice. More bins follow the "
                           "wall more closely and are noisier."},
        r_smooth_bins={"widget_type": "Slider", "min": 0, "max": 50, "step": 1,
                       "tooltip": "polar only: how many neighbouring angular "
                                  "bins the radius is smoothed over. 0 is no "
                                  "smoothing."},
        alpha2d_um={"widget_type": "FloatSlider", "min": 0.5, "max": 30.0,
                    "step": 0.5,
                    "tooltip": "alpha 2D only: the alpha radius, in µm, within "
                               "one slice. A gap wider than this is left open."},
        subsample_um={"widget_type": "FloatSlider", "min": 0.0, "max": 5.0,
                      "step": 0.25,
                      "tooltip": "alpha 3D: keep one wall point per cube of "
                                 "this size before building the alpha complex, "
                                 "which is what makes it affordable. 0 keeps "
                                 "every voxel. Must stay well below alpha3d_um "
                                 "or the complex is built on points further "
                                 "apart than the radius meant to connect them."},
        alpha3d_um={"widget_type": "FloatSlider", "min": 1.0, "max": 80.0,
                    "step": 0.5,
                    "tooltip": "alpha 3D: the alpha radius in µm. This is the "
                               "widest wall gap that gets bridged; the gaps in "
                               "this data run 8-10 µm."},
        close_um={"widget_type": "FloatSlider", "min": 0.0, "max": 5.0,
                  "step": 0.25,
                  "tooltip": "Radius of a morphological closing applied to the "
                             "filled slice, in µm. 0 is off."},
        min_component_um={"widget_type": "FloatSlider", "min": 0.0,
                          "max": 500.0, "step": 5.0,
                          "tooltip": "Discard reconstructed pieces smaller "
                                     "than this across, in µm. Removes "
                                     "detached fragments; set it above the "
                                     "vessel and it removes the vessel."},
        poisson_grid_um={"widget_type": "FloatSlider", "min": 0.65, "max": 6.0,
                         "step": 0.325,
                         "tooltip": "poisson: the cell size of the coarse grid "
                                    "the field is solved on, in µm. This is "
                                    "the resolution of everything downstream — "
                                    "the surface, the lumen, the "
                                    "cross-sections and the shear. Cost grows "
                                    "as the cube as it shrinks."},
        poisson_smooth_um={"widget_type": "FloatSlider", "min": 1.0,
                           "max": 15.0, "step": 0.5,
                           "tooltip": "poisson: the width in µm over which the "
                                      "wall direction is estimated. It also "
                                      "sets the lumen test: a vessel whose "
                                      "inscribed radius is under 4x this "
                                      "cannot be oriented at all and the build "
                                      "refuses rather than guessing."},
        poisson_iso_pct={"widget_type": "FloatSlider", "min": 5.0, "max": 95.0,
                         "step": 1.0,
                         "tooltip": "poisson: which level of the solved field "
                                    "becomes the surface, as a percentile of "
                                    "the field sampled where the cadherin "
                                    "actually is. 50 puts the surface at the "
                                    "middle of the wall. Higher shrinks the "
                                    "lumen, lower grows it — and the lumen "
                                    "size drives the shear."},
        alpha_raster_um={"widget_type": "FloatSlider", "min": 0.325, "max": 6.0,
                         "step": 0.325,
                         "tooltip": "alpha 3D solid -> poisson only: the voxel "
                                    "size the alpha solid is filled in at "
                                    "before it is handed to the Poisson fit."},
        fill_lumen={"tooltip": "Per-slice methods: fill the interior of each "
                               "slice's wall before meshing, so the result is "
                               "a solid cross-section rather than a ring."},
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
                _drop_shear()
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
            _drop_shear()
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

    suggest_w = PushButton(
        text="Suggest alpha_um (3D)",
        tooltip="Build one Delaunay triangulation of the subsampled wall points and print a suggested alpha3d_um, with a histogram. Builds no mesh and changes no setting.")

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
    mesh_path_w = FileEdit(
        tooltip="Where 'Save vessel mesh' writes and 'Load vessel mesh' reads. Relative names resolve beside the current image.",value="vessel_mesh.npz", label="mesh_path", mode="w")
    save_mesh_w = PushButton(
        text="Save vessel mesh",
        tooltip="Write the built mesh to the path above, beside this image.")
    load_mesh_w = PushButton(
        text="Load vessel mesh",
        tooltip="Read a mesh back from the path above and draw it. It does not rebuild, so the sliders may not describe what you are looking at.")

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
            _drop_shear()
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

    # A new mesh invalidates the wall-shear field: the shear was measured on
    # the lumen of the *previous* build, and `vertex_field` would happily map it
    # onto whatever vertices it is now handed. That is [P38](docs/ISSUES.md#p38)
    # exactly — one build's geometry surviving into another's — so the cloud and
    # its layer are dropped here rather than left to look plausible.
    SHEAR_LAYER = "Vessel wall shear"

    def _drop_shear():
        state["vessel_shear"] = None
        if SHEAR_LAYER in viewer.layers:
            viewer.layers.remove(SHEAR_LAYER)

    # ---- Per-vertex curvature (local quadric fit at a physical radius) ----
    CURV_CHOICES = {"mean (H)": "mean", "κ1 max": "k1",
                    "κ2 min": "k2", "gaussian (K)": "gaussian"}
    curv_measure_w = ComboBox(
        tooltip="Which curvature of the fitted surface is drawn.",choices=list(CURV_CHOICES.keys()),
                              value="mean (H)", label="curvature")
    curv_radius_w = FloatSlider(value=5.0, min=1.0, max=50.0, step=0.5,
                                label="curvature_radius_um",
        tooltip='Radius in µm of the neighbourhood the surface curvature is fitted over. There is no right value on this data: the answer tracks this slider instead of settling, which is what docs/ISSUES.md P24 is about. Read a curvature from here only if it holds steady as you move this.')
    curv_smooth_w = Slider(value=10, min=0, max=50,
                           label="curvature_smooth_iter",
        tooltip='How many smoothing passes the mesh gets before curvature is fitted. 0 is none.')
    curv_btn = PushButton(
        text="Compute curvature",
        tooltip="Fit curvature on the mesh currently built or loaded, and draw "
                "it as its own layer. Read curvature_radius_um's tooltip "
                "first — the value depends on it.")
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
    spur_ratio_w = FloatSlider(
        value=lumen_graph.DEFAULT_SPUR_RATIO, min=0.0, max=4.0, step=0.1,
        label="spur_ratio",
        tooltip="A dead-end branch shorter than this many local vessel radii "
                "is treated as a bump on the reconstruction and removed. A "
                "ratio and not a length, so it means the same thing at a 40 µm "
                "caliber and at a 5 µm waist. 1.0 reads 'a stub that does not "
                "reach as far as the vessel is wide'. 0 prunes nothing.")
    crop_ratio_w = FloatSlider(
        value=lumen_graph.DEFAULT_CROP_RATIO, min=0.0, max=6.0, step=0.25,
        label="crop_ratio",
        tooltip="Two jobs. A centerline end within this many local radii of "
                "the edge of the image is where the vessel leaves the field of "
                "view, rather than a real vessel ending. And every section "
                "within that distance is flagged `censored`, because the lumen "
                "there is no longer fully imaged. 2.0 is not arbitrary: the "
                "centerline of a tube cut off by a flat plane stops about one "
                "radius short of the cut, so a smaller value calls every "
                "cropped end a real terminus.")
    kappa_win_w = FloatSlider(
        value=lumen_graph.DEFAULT_KAPPA_WINDOW_UM, min=5.0, max=120.0,
        step=5.0, label="kappa_window_um",
        tooltip="Arc length, in µm, over which the bend of the centerline is "
                "fitted. It is a second derivative, so it is meaningless below "
                "the vessel width and is smoothed away far above it; 40 µm is "
                "about one diameter here.")
    sec_every_w = FloatSlider(
        value=lumen_graph.DEFAULT_SECTION_EVERY_UM, min=1.0, max=40.0,
        step=1.0, label="section_every_um",
        tooltip="Cut a cross-section every this many µm along the vessel. "
                "Between cuts nothing is measured and nothing is interpolated, "
                "so this is literally how coarse the shear map is along the "
                "axis. It is also the distance a mesh vertex may be from a "
                "measured wall point before it is left uncoloured. Cutting "
                "closer is a real measurement and costs seconds.")
    sec_span_w = FloatSlider(
        value=lumen_graph.DEFAULT_SECTION_SPAN, min=1.0, max=8.0, step=0.25,
        label="section_span",
        tooltip="How far out the cutting plane is sampled, as a multiple of "
                "the local inscribed radius. That radius is the *short* "
                "half-width of the section, so this must exceed the section's "
                "width-to-height ratio or the cut runs off the sampled disc "
                "and the section is flagged `clipped` — its size then being a "
                "lower bound rather than a measurement.")
    # §4.3's wall-shear solve, matching `batch_mesh --shear`. Present here for
    # the standing rule that the GUI and the headless path must not drift: a
    # check that exists in only one of them is how P36 and the alpha3d_um case
    # happened. Off by default in both.
    shear_w = CheckBox(
        value=False, text="wall shear (sec 4.3)",
        tooltip="Solve the flow in every cross-section and get the wall shear "
                "around each one. Off by default because it is extra work per "
                "section. Nothing below in the shear group does anything until "
                "this has been run once.")
    # Present here only so the tab and `batch_mesh --no-section-from-field`
    # cannot disagree about what a section is. On in both; unticking reproduces
    # a pre-2026-08-28 build, which is the only reason it is a control at all.
    sec_field_w = CheckBox(
        value=True, text="cut sections from chi (sub-voxel wall)",
        tooltip="On: the cut plane is sampled from the solved field chi by "
                "interpolation and the wall is its level curve, so the "
                "section's wall is the same curve the displayed mesh is made "
                "of. Off: the plane is cut from the binarised lumen and every "
                "wall lands on a pixel edge — perimeter ~11% long, force "
                "balance ~1.04 instead of ~1.00 (docs/ISSUES.md P46).")
    graph_color_w = ComboBox(
        tooltip="Which per-station quantity colours the centerline points. "
                "A view of the graph only; nothing here feeds a calculation.",
        choices=["caliber_um", "kappa_1_um", "kappa_xy_1_um", "tau_rel_1_um3",
                 "tau_hat_mean_um", "hydraulic_diam_um", "lumen_area_um2",
                 "lumen_major_um", "lumen_minor_um", "lumen_aspect",
                 "dist_junction_um", "dist_crop_um", "censored"],
        value="caliber_um", label="colour by")
    graph_btn = PushButton(
        text="Compute 3D centerline graph (lumen)",
        tooltip="Trace the vessel centerline through the solved lumen, cut a cross-section every section_every_um along it, and — with 'wall shear' ticked — solve the flow in each. Needs a Poisson-family method, since the interpolating ones never solve for a lumen volume.")

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
                "section_span": float(sec_span_w.value),
                "section_from_field": bool(sec_field_w.value),
                "compute_shear": bool(shear_w.value)}

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
            ex: dict = {}
            g = lumen_graph.build_graph(
                lum["chi"], lum["iso_level"], lum["h_um"], lum["origin_um"],
                core=lum.get("core"), domain_shape_um=dom, extras=ex, **gp)
            return g, ex

        def _done(result):
            g, ex = result
            state["vessel_graph"] = g
            # The wall point cloud, when the shear box is ticked. Cleared
            # otherwise, so an untick cannot leave the previous build's shear
            # painted on a new mesh — that is [P38](docs/ISSUES.md#p38)'s
            # failure mode, one image's field surviving into another's build.
            state["vessel_shear"] = ex.get("wall")
            _show_shear()
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

    # ---- Wall shear on the vessel surface -----------------------------------
    # The same read as `curvature` above: a per-vertex scalar on the mesh, drawn
    # with `add_surface`. Two deliberate differences, both because shear is
    # measured on **part** of the surface rather than defined on all of it:
    #
    #  - **nearest measured wall point, never interpolated.** The colour under
    #    the cursor is a number that was solved for, not a blend of two. Where
    #    the banding along the axis is visible, the honest fix is a smaller
    #    `section_every_um` — a real measurement that costs seconds — not a
    #    smoother that would make a 5 µm sampling look like a continuous field.
    #  - **un-measured wall is not drawn.** `curvature` fills its NaN vertices
    #    with the median, which is fine for a field defined on every vertex and
    #    actively misleading here: it would paint wall that was never sectioned
    #    in a mid-range colour that reads as a measurement.
    #
    # This is a view of `extras["wall"]` and nothing else — the same array
    # `batch_mesh` writes to `vessel_shear.npz` and the same one an adhesion
    # will join to. Nothing here feeds a calculation.
    SHEAR_CHOICES = {"τ / µQ (relative)": "tau_rel_1_um3",
                     "τ̂ (shape function, µm)": "tau_hat_um"}
    shear_measure_w = ComboBox(
        choices=list(SHEAR_CHOICES.keys()), value="τ / µQ (relative)",
        label="wall shear",
        tooltip="Which quantity is painted on the surface. 'τ / µQ' holds the "
                "flow rate fixed and is the one to compare along a vessel: it "
                "rises where the vessel narrows, as it should. 'τ̂' holds the "
                "pressure gradient fixed instead, which inverts that — a "
                "narrow section passes less flow and reads LOW. τ̂ is also "
                "exactly the section's area divided by its perimeter, by force "
                "balance, so it carries no information the geometry columns "
                "do not already have.")
    # Two flags, two boxes. They were one box until 2026-08-28 and that hid the
    # difference between them: on 04step_003, 4 of 88 sections are `clipped`
    # and 23 are `censored`, and no section is both, so a single control that
    # dropped 34% of the wall points read as "the solver failed here".
    #
    #  - **clipped** — the section contour ran off the edge of the sampling
    #    disc (`section_span` x caliber) or the volume, so its extent is a
    #    lower bound. Not a measurement; dropped by default. On 04step_003 all
    #    four are the stub at the Y junction, where the inscribed radius
    #    collapses to ~10 um while the cut plane spans 65-75 um across both
    #    arms.
    #  - **censored** — the section sits within `crop_ratio` calibers of a
    #    field-of-view crop. That is a statement about the *premise*, not the
    #    solve: `force_balance` there is 1.033 against 1.044 in the interior
    #    and `wall_ok_frac` is 1.000, so the solve is healthy and the section
    #    may simply not be the whole vessel. Shown by default, because hiding
    #    a 70 um cuff of wall that has real signal on it is the more
    #    misleading of the two errors — but Q ~ A^2, so the absolute
    #    tau / uQ out there is not to be trusted even though the
    #    circumferential *shape* probably is. The summary line says how much
    #    of what is drawn is censored, every time, so that caveat cannot go
    #    quiet.
    shear_drop_clipped_w = CheckBox(
        value=True, text="drop clipped sections (contour ran off the disc)",
        tooltip="The section contour reached the edge of the sampling disc "
                "(section_span x caliber) or of the volume, so its extent is "
                "a lower bound rather than a measurement. Usually the stub at "
                "a junction, where the inscribed radius collapses but the cut "
                "plane still spans both arms.")
    shear_drop_censored_w = CheckBox(
        value=False, text="drop censored sections (near a field-of-view crop)",
        tooltip="The section sits within crop_ratio calibers of a "
                "field-of-view crop, so the lumen leaves the image here and "
                "the section may not be the whole vessel. The solve itself is "
                "healthy (force_balance ~1.03, wall_ok_frac 1.000); it is the "
                "premise that is flagged. Q ~ A^2, so absolute tau / uQ out "
                "there is untrustworthy even where the circumferential shape "
                "is not.")

    # τ/µQ is normalised by **conservation of Q**, which holds along an
    # unbranched segment and fails at a Y: the split between daughters rides on
    # h³ and is not recoverable from geometry (docs/DIRECTION.md §4.2). So each
    # edge of the graph carries its own unknown µQ, and painting all of them on
    # one colour scale invites exactly the comparison the quantity cannot
    # support — [P49](docs/ISSUES.md#p49). Two controls, both display-only:
    #
    #  - **normalise per edge** divides each edge by its own median, so the
    #    colour means "relative to this edge" and no cross-edge reading is
    #    implied. The per-edge medians are printed, so the absolute scale is
    #    not lost, only kept off the colour bar.
    #  - **hide near junctions** greys the zone where the per-section duct
    #    model does not apply at all — no flow split, no inertia, no entrance
    #    length — the same way censoring greys the mouths. As a ratio of the
    #    local caliber, for [P23]'s reason.
    # How the edges are put on a common scale. Three answers, in decreasing
    # order of how much they claim:
    #
    #  - **flow split** solves the resistance network (`shear.edge_resistance`
    #    + `shear.flow_split`) and multiplies each edge by its own relative
    #    flow. §4.2 says the split "rides on h³ and is not recoverable from
    #    geometry" — true of the closed form, not of the vessel: the per-
    #    section solve gives the conductance of every station, so R̂ = ∫ds/Q̂
    #    is an integral we can evaluate and Kirchhoff does the rest. This
    #    leaves **one** unknown µQ for the whole vessel instead of one per
    #    edge, which is what [P49](docs/ISSUES.md#p49) was working around.
    #  - **per edge median** claims nothing and shows only within-edge
    #    structure. The safe fallback, and what is used when there is no graph
    #    beside the shear file.
    #  - **raw** is the old picture, and prints why it is misleading.
    SHEAR_SCALES = ["flow split (Kirchhoff)", "per edge median", "raw"]
    shear_scale_w = ComboBox(
        choices=SHEAR_SCALES, value=SHEAR_SCALES[0], label="edge scaling",
        tooltip="How the branches are put on one colour scale. Each branch "
                "carries its own unknown flow rate, so the raw numbers are "
                "comparable within a branch and not across one. 'flow split' "
                "computes each branch's resistance and solves for the flows, "
                "which makes them comparable; 'per edge median' divides each "
                "branch by its own middle value and claims nothing about how "
                "they compare; 'raw' is neither and prints a warning.")
    shear_inlet_w = ComboBox(
        choices=[-1], value=-1, label="inlet node",
        tooltip="Which opening the flow enters through. This is a fact about "
                "the perfusion rig and is not in the image, so it is a control "
                "and not an inference; the choices are the graph's own node "
                "numbers, shown as red/blue/yellow points by the graph layer. "
                "It only affects 'flow split'. Defaults to the widest opening.")
    shear_terminus_open_w = CheckBox(
        value=True, text="treat a terminus endpoint as an opening",
        tooltip="An endpoint that stops inside the volume is a dead end and "
                "carries no flow. docs/DIRECTION.md §5.1 expects no real "
                "vessel terminus in this dataset, so one usually means the "
                "lumen was reconstructed short — open by default. Untick to "
                "see how much the split depends on it; on 04step_003 it is "
                "the difference between the stub taking 69% of the flow and "
                "taking none.")
    shear_junction_ratio_w = FloatSlider(
        value=2.0, min=0.0, max=6.0, step=0.25, label="hide_junction_ratio",
        tooltip="Hide wall points within this many local calibers of a "
                "junction. The solve is a fully developed per-section duct "
                "flow; at a Y none of its assumptions hold and Q is not even "
                "conserved. 0 shows them anyway.")

    def _open_nodes(g):
        """Endpoint nodes that flow may enter or leave through.

        A crop endpoint is where the vessel leaves the field of view, so it is
        an opening. A terminus stops inside the volume; that is a dead end
        unless it is a lumen reconstructed short, which docs/DIRECTION.md §5.1
        says is the likely case here — hence the checkbox rather than a rule.
        """
        deg = np.asarray(g["node_degree"], int)
        kind = np.asarray(g["node_endpoint_kind"]).astype(str)
        ends = np.nonzero(deg == 1)[0]
        if shear_terminus_open_w.value:
            return ends
        return ends[kind[ends] == "crop"]

    def _sync_inlet_choices(g):
        """Offer the openings as inlet choices; default the widest.

        The inlet is a fact about the rig, not about the image, so it is a
        control and not an inference. The default is the largest-caliber
        opening — the trunk of a Y usually is — and the choice is printed with
        the split every time, because the split depends on it.
        """
        opens = [int(n) for n in _open_nodes(g)]
        if not opens:
            shear_inlet_w.choices = [-1]
            return
        cal = np.asarray(g.get("node_caliber_um", []), float)
        cur = int(shear_inlet_w.value) if shear_inlet_w.value is not None else -1
        shear_inlet_w.choices = opens
        if cur not in opens:
            widest = max(opens, key=lambda n: (cal[n] if n < len(cal)
                                               else -np.inf))
            shear_inlet_w.value = widest

    def _flow_split():
        """Relative flow per edge, or None when it cannot be solved.

        Resistances come from **all** finite sections, not from the ones the
        display filters kept: dropping a clipped section punches a hole in
        ∫ds/Q̂, and since the integrand is 1/Q̂ the hole is exactly where the
        integral is largest. The flags travel with the answer instead — the
        printed `span_frac` says how much of each edge the sections cover.
        """
        g = state.get("vessel_graph")
        if not g or "edges" not in g or "q_hat_um4" not in g:
            return None
        edges = np.asarray(g["edges"], int)
        if not len(edges):
            return None
        _sync_inlet_choices(g)
        opens = _open_nodes(g)
        inlet = int(shear_inlet_w.value)
        if inlet not in set(int(n) for n in opens):
            return None
        R = shear.edge_resistance(g["edge_id"], g["s_um"], g["q_hat_um4"])
        fs = shear.flow_split(edges, {e: v["r_hat_1_um3"]
                                      for e, v in R.items()}, inlet, opens)
        if not fs.get("q_rel"):
            return None
        fs["r_table"] = [(e, R[e]["r_hat_1_um3"], R[e]) for e in sorted(R)]
        return fs

    def _show_shear():
        wall = state.get("vessel_shear")
        cache = state.get("vessel_mesh_cache") or {}
        verts, faces = cache.get("verts"), cache.get("faces")
        if SHEAR_LAYER in viewer.layers:
            viewer.layers.remove(SHEAR_LAYER)
        if not wall or verts is None or faces is None:
            return
        pts = np.asarray(wall["points_zyx_um"], float)
        vals = np.asarray(wall[SHEAR_CHOICES[shear_measure_w.value]], float)
        clipped = np.asarray(wall["section_clipped"], bool)
        censored = np.asarray(wall["censored"], bool)
        keep = np.ones(len(pts), bool)
        if shear_drop_clipped_w.value:
            keep &= ~clipped
        if shear_drop_censored_w.value:
            keep &= ~censored
        # The junction zone. `dist_junction_um` is all-NaN on an unbranched
        # image (three of six are), and absent from a shear file written
        # before [P49] — neither is an error, both mean "no zone to hide".
        dj = np.asarray(wall.get("dist_junction_um", []), float)
        cal = np.asarray(wall.get("caliber_um", []), float)
        near_junction = np.zeros(len(pts), bool)
        jr = float(shear_junction_ratio_w.value)
        if jr > 0 and len(dj) == len(pts) and len(cal) == len(pts):
            with np.errstate(invalid="ignore"):
                near_junction = np.isfinite(dj) & (dj <= jr * cal)
            keep &= ~near_junction
        elif jr > 0 and len(pts) and len(dj) != len(pts):
            print("Wall shear: this shear file predates the per-point "
                  "junction distance, so 'hide_junction_ratio' is inactive. "
                  "Re-run the graph to record it.", flush=True)

        n_all = len(pts)
        edge = np.asarray(wall.get("edge_id", np.zeros(n_all)), int)
        pts, vals, cens_kept = pts[keep], vals[keep], censored[keep]
        edge_kept = edge[keep] if len(edge) == n_all else np.zeros(len(pts), int)
        if not len(pts):
            print(f"Wall shear: no wall points left — of {n_all:,}, "
                  f"{int(clipped.sum()):,} are clipped, "
                  f"{int(censored.sum()):,} censored and "
                  f"{int(near_junction.sum()):,} within {jr:g} calibers of a "
                  "junction. Untick a filter or set the ratio to 0 to see "
                  "them anyway; each control's tooltip says what it means.",
                  flush=True)
            return

        # Putting the edges on a common scale. Done on the **values**, before
        # they are carried to vertices, so a vertex takes the scaled number of
        # the point it is nearest to rather than a blend of two edges'.
        scale_mode = str(shear_scale_w.value)
        edge_meds: list[tuple[int, float, int]] = []
        split_note = ""
        if scale_mode == "flow split (Kirchhoff)":
            fs = _flow_split()
            if fs is None:
                scale_mode = "per edge median"
                split_note = ("no centerline graph beside this shear file, so "
                              "the flow split cannot be solved — fell back to "
                              "the per-edge median")
            else:
                for i in fs["q_rel"]:
                    m = edge_kept == i
                    if m.any():
                        edge_meds.append((int(i), abs(fs["q_rel"][i]),
                                          int(m.sum())))
                        vals = np.where(m, vals * abs(fs["q_rel"][i]), vals)
                missing = sorted({int(e) for e in np.unique(edge_kept)}
                                 - set(fs["q_rel"]))
                # An edge with no solved flow is **hidden, not drawn at zero**:
                # zero is a value and it would read as "no shear here" when
                # what happened is that the network solve had nothing to say.
                if missing:
                    keep_e = ~np.isin(edge_kept, missing)
                    pts, vals = pts[keep_e], vals[keep_e]
                    cens_kept, edge_kept = cens_kept[keep_e], edge_kept[keep_e]
                    split_note = (f"edge(s) {missing} carry no solved flow and "
                                  f"are not drawn")
                if not len(pts):
                    print("Wall shear: the flow split left no edge with a "
                          "usable scale. Switch 'edge scaling' to per-edge "
                          "median.", flush=True)
                    return
                fs_info = fs
        if scale_mode == "per edge median":
            for e in np.unique(edge_kept):
                m = edge_kept == e
                med = float(np.nanmedian(vals[m])) if np.isfinite(
                    vals[m]).any() else np.nan
                edge_meds.append((int(e), med, int(m.sum())))
                if np.isfinite(med) and med > 0:
                    vals = np.where(m, vals / med, vals)
        per_vertex = shear.vertex_field(verts, pts, vals,
                                        float(wall["max_dist_um"]))
        keep_faces = shear.covered_faces(faces, per_vertex)
        if not len(keep_faces):
            print("Wall shear: no mesh face has a measurement on all three "
                  "vertices.", flush=True)
            return
        finite = per_vertex[np.isfinite(per_vertex)]
        lo, hi = np.percentile(finite, [2, 98])
        fill = float(finite.min())
        vv = np.where(np.isfinite(per_vertex), per_vertex, fill)
        lyr = viewer.add_surface((verts, keep_faces, vv), name=SHEAR_LAYER,
                                 colormap="turbo", blending="translucent",
                                 opacity=0.95)
        try:
            lyr.contrast_limits = (float(lo), float(hi))
        except Exception:
            pass
        cov = float(np.isfinite(per_vertex).mean())
        print(f"Wall shear: {len(pts):,} wall points -> "
              f"{int(np.isfinite(per_vertex).sum()):,} of {len(verts):,} "
              f"vertices ({100 * cov:.0f}%), {len(keep_faces):,} of "
              f"{len(faces):,} faces drawn. Nearest measured point within "
              f"{wall['max_dist_um']:g} µm; the rest is not drawn because it "
              f"was not measured.", flush=True)
        if split_note:
            print(f"  NOTE: {split_note}.", flush=True)
        if scale_mode == "flow split (Kirchhoff)" and edge_meds:
            body = ", ".join(f"edge {e}: x{q:.3f} ({n:,} pts)"
                             for e, q, n in edge_meds)
            print(f"  scaled by the solved flow split, inlet node "
                  f"{fs_info['inlet_node']}, "
                  f"{'termini open' if shear_terminus_open_w.value else 'termini closed'}"
                  f" — one µQ now covers the whole vessel instead of one per "
                  f"edge ([P49], [I10]). Relative flow: {body}. Resistances "
                  f"R̂ = ∫ds/Q̂ in µm^-3: "
                  + ", ".join(f"edge {e}: {r:.3g} (sections span "
                              f"{100 * v['span_frac']:.0f}% of it)"
                              for e, r, v in fs_info["r_table"]) + ".",
                  flush=True)
            for n in fs_info.get("notes", []):
                print(f"  {n}", flush=True)
        elif scale_mode == "per edge median" and edge_meds:
            body = ", ".join(f"edge {e}: {med:.3g} ({n:,} pts)"
                             for e, med, n in edge_meds)
            print(f"  normalised per edge — colour is x the edge's own "
                  f"median. This claims nothing about how the branches "
                  f"compare; the flow-split option does. Medians in the raw "
                  f"unit: {body}.", flush=True)
        elif len(np.unique(edge_kept)) > 1:
            print(f"  {len(np.unique(edge_kept))} edges on one colour scale, "
                  f"**raw**. τ/µQ carries a separate unknown µQ per edge, so "
                  f"read each branch on its own and no colour across the "
                  f"junction ([P49]). 'flow split' solves for the missing "
                  f"factors and makes them comparable ([I10]).", flush=True)
        if near_junction.any():
            print(f"  {int(near_junction.sum()):,} point(s) hidden within "
                  f"{jr:g} calibers of a junction, where the fully developed "
                  f"per-section duct model does not apply.", flush=True)
        # Always, not only when it is large: what is on screen is a mix of two
        # populations and the colour does not distinguish them.
        if cens_kept.any():
            print(f"  {int(cens_kept.sum()):,} of the {len(pts):,} points "
                  f"drawn ({100 * cens_kept.mean():.0f}%) are **censored** — "
                  f"within `crop_ratio` calibers of a field-of-view crop, "
                  f"where the section may not be the whole vessel. Q ~ A², so "
                  f"read their circumferential shape, not their absolute "
                  f"τ / µQ.", flush=True)

    shear_measure_w.changed.connect(lambda *_: _show_shear())
    shear_drop_clipped_w.changed.connect(lambda *_: _show_shear())
    shear_drop_censored_w.changed.connect(lambda *_: _show_shear())
    shear_scale_w.changed.connect(lambda *_: _show_shear())
    shear_inlet_w.changed.connect(lambda *_: _show_shear())
    shear_terminus_open_w.changed.connect(lambda *_: _show_shear())
    shear_junction_ratio_w.changed.connect(lambda *_: _show_shear())

    shear_path_w = FileEdit(
        tooltip="Where the wall-shear point cloud is written and read. Kept separate from the graph file because it has one row per wall point rather than per centerline point.",value="vessel_shear.npz", label="shear_path",
                            mode="w")
    save_shear_w = PushButton(
        text="Save wall shear",
        tooltip="Write the wall-shear point cloud: one row per wall point, with its position, its shear and its flags.")
    load_shear_w = PushButton(
        text="Load wall shear",
        tooltip="Read a wall-shear point cloud back and paint it on the current mesh.")

    def _save_shear(*_):
        wall = state.get("vessel_shear")
        if not wall:
            print("No wall shear to save — tick 'wall shear (sec 4.3)' and "
                  "rebuild the centerline graph.", flush=True)
            return
        path = str(shear_path_w.value)
        np.savez_compressed(
            path, **{k: v for k, v in wall.items()
                     if isinstance(v, np.ndarray)},
            max_dist_um=np.float64(wall["max_dist_um"]))
        print(f"Saved wall shear -> {path}.", flush=True)

    def _load_shear(*_, after=None, visible=True):
        path = str(shear_path_w.value)
        try:
            with np.load(path, allow_pickle=False) as d:
                wall = {k: d[k] for k in d.files}
            wall["max_dist_um"] = float(wall["max_dist_um"])
            state["vessel_shear"] = wall
            if visible:
                _show_shear()
            print(f"Loaded wall shear <- {path} "
                  f"({len(wall['points_zyx_um']):,} wall points).", flush=True)
        except Exception as exc:
            print(f"Load wall shear error: {exc}", flush=True)
        finally:
            call_all(after, label="Load wall shear")

    def _save_shear_if_any(*_):
        """Save only if the shear box was ticked — silent when it was not.

        `_save_shear` prints a complaint when there is nothing to save, which is
        right for the button and wrong for the "Run all" chain, where not
        computing shear is the default and not an error.
        """
        if state.get("vessel_shear"):
            _save_shear()

    save_shear_w.changed.connect(_save_shear)
    load_shear_w.changed.connect(_load_shear)
    state["load_vessel_shear"] = _load_shear

    # ---- Save / load the computed curvature so it need not be recomputed each
    # session (only the mesh geometry is saved with the mesh; curvature and the
    # centerline are separate). Surface curvature -> the 4 per-vertex arrays;
    # path curvature -> the centerline points, κ and branch nodes. ----
    curv_path_w = FileEdit(
        tooltip="Where the per-vertex surface curvature is written and read.",value="vessel_curvature.npz", label="curv_path",
                           mode="w")
    save_curv_w = PushButton(
        text="Save curvature",
        tooltip="Write the per-vertex curvature to the path above.")
    load_curv_w = PushButton(
        text="Load curvature",
        tooltip="Read per-vertex curvature back. It must have been computed on the mesh now loaded, or its length will not match.")

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

    path_curv_path_w = FileEdit(
        tooltip="Where the centerline graph is written and read.",value="vessel_graph.npz",
                                label="graph_path", mode="w")
    save_path_curv_w = PushButton(
        text="Save centerline graph",
        tooltip="Write the centerline graph: points, edges, per-station geometry and the cross-section columns.")
    load_path_curv_w = PushButton(
        text="Load centerline graph",
        tooltip="Read a centerline graph back and draw it.")

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
        (shear_path_w, "vessel_shear.npz"),
    ])

    # Chainable "Run all" steps: compute then save, then call `after`.
    state["run_vessel_mesh"] = lambda after=None: _run_build(
        after=lambda: call_all(_save_mesh, after, label="run_vessel_mesh"))
    state["run_vessel_curvature"] = lambda after=None: _compute_curv(
        after=lambda: call_all(_save_curv, after,
                               label="run_vessel_curvature"))
    state["run_vessel_path_curvature"] = lambda after=None: _compute_graph(
        after=lambda: call_all(_save_path_curv, _save_shear_if_any, after,
                               label="run_vessel_path_curvature"))

    # Apply the visibility rule once at build time (`suggest_w` exists by now).
    _sync_method_visibility()
    # The call button is made by `@magicgui` and takes no options there, so its
    # tooltip is set on the widget afterwards.
    controls.call_button.tooltip = (
        "Build the vessel surface from the Cadherin tab's current mask, with "
        "the method and knobs above. The knobs a method does not read are "
        "hidden rather than left on screen looking operative.")

    return (Container(widgets=[controls, suggest_w,
                              mesh_path_w, save_mesh_w, load_mesh_w,
                              curv_measure_w, curv_radius_w, curv_smooth_w,
                              curv_btn, curv_path_w, save_curv_w, load_curv_w,
                              spur_ratio_w, crop_ratio_w, kappa_win_w,
                              sec_every_w, sec_span_w, sec_field_w, shear_w,
                              graph_color_w,
                              graph_btn,
                              shear_measure_w, shear_drop_clipped_w,
                              shear_drop_censored_w, shear_scale_w,
                              shear_inlet_w, shear_terminus_open_w,
                              shear_junction_ratio_w,
                              shear_path_w, save_shear_w, load_shear_w,
                              path_curv_path_w, save_path_curv_w,
                              load_path_curv_w],
                     # **Labels on.** With `labels=False` a magicgui Container
                     # draws no name beside a slider or a combo, so every knob
                     # in this tab that is not a checkbox (whose text is drawn
                     # inside it) appeared as a bare unlabelled control. The
                     # `label=` values were all there and none of them reached
                     # the screen.
                     labels=True),
            hist.canvas)
