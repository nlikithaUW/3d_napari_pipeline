"""Stack loader tab: project_dir + image picker + Load button."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from magicgui.widgets import ComboBox, Container, FileEdit, PushButton
from napari.qt.threading import thread_worker

from ..fret import compute_fret_volume, load_force_lut, resample_fret_outputs
from ..io import load_params, load_stack, to_isotropic
from .widgets import percentile_limits, pyramid_views

FRET_REQUIRED = {"clover", "ruby", "fret"}

# Percentile contrast defaults applied to raw channels on load (napari's
# min/max auto-contrast magnifies the noise floor — see percentile_limits).
LO_PCT = 1.0
HI_PCT = 99.9

# (role, colormap, contrast_limits_or_None, visible)
DISPLAY = [
    ("nuclei",    "blue",    None,  True),
    ("substrate", "magenta", None,  True),
    ("clover",    "green",   None,  True),
    ("ruby",      "red",     None,  False),
    ("fret",      "gray",    None,  False),
]


def _load_image(state, viewer, tif_path, params_path, lut_path, levels,
                after=None):
    print(f"Loading {tif_path.name} ...", flush=True)

    @thread_worker
    def _run():
        params = load_params(params_path)
        stack_raw = load_stack(tif_path, params)
        role_to_c = {role: stack_raw.channels.index(name)
                     for role, name in stack_raw.roles.items()}
        missing = FRET_REQUIRED - set(stack_raw.roles)
        if missing:
            raise ValueError(
                f"{tif_path.name}: missing roles {sorted(missing)} required "
                f"for FRET/force computation"
            )
        lut_eff, lut_force = load_force_lut(str(lut_path))
        fret = compute_fret_volume(stack_raw.data, role_to_c, params["fret"],
                                   lut_eff, lut_force)
        raw_shape = stack_raw.data.shape
        raw_nbytes = stack_raw.data.nbytes
        raw_voxel = stack_raw.voxel_zyx_um
        # Keep native-Z (pre-resample) clover/ruby for ratio quantification on
        # the original acquisition grid (.copy() so deleting stack_raw frees the
        # rest of the raw stack). Segmentation still runs on the isotropic stack.
        clover_native = (stack_raw.data[:, role_to_c["clover"]].copy()
                         if "clover" in stack_raw.roles else None)
        ruby_native = (stack_raw.data[:, role_to_c["ruby"]].copy()
                       if "ruby" in stack_raw.roles else None)
        stack = to_isotropic(stack_raw)
        fret = resample_fret_outputs(fret, stack_raw.voxel_zyx_um,
                                     stack.voxel_zyx_um[0])
        # Drop the raw uint16 stack — closure would otherwise retain it.
        del stack_raw
        return (raw_shape, raw_nbytes, raw_voxel, stack, role_to_c, fret,
                clover_native, ruby_native)

    def _done(result):
        (raw_shape, raw_nbytes, raw_voxel, stack, role_to_c, fret,
         clover_native, ruby_native) = result
        print(f"  raw   shape={raw_shape} voxel(z,y,x)={raw_voxel} µm "
              f"({raw_nbytes/1e6:.1f} MB)", flush=True)
        print(f"  iso   shape={stack.data.shape} dtype={stack.data.dtype} "
              f"voxel(z,y,x)={stack.voxel_zyx_um} µm "
              f"({stack.data.nbytes/1e6:.1f} MB)", flush=True)

        # Clear the cell-mesh cache BEFORE touching layers. Adding/removing
        # layers changes `viewer.dims` and fires `current_step`, which runs the
        # persistent `_update_section` hook. If the previous image's mesh is
        # still cached, that hook re-creates the "Cell mesh (slice)" layer from
        # geometry sliced at the new image's z-plane — the source of the
        # NaN-segment GL crash when loading a different-depth image.
        state["mesh_cache"] = {"offset_centroids": None, "outward": None,
                               "inward_offset_um": 0.0, "kdtree": None,
                               "verts": None, "faces": None}

        for layer in list(viewer.layers):
            viewer.layers.remove(layer)

        scale = stack.voxel_zyx_um
        ml = levels > 1

        raw_layers = []
        for role, cmap, clim, visible in DISPLAY:
            if role not in stack.roles:
                continue
            c = role_to_c[role]
            data = pyramid_views(stack.data[:, c], levels)
            kwargs = dict(name=stack.roles[role], colormap=cmap, scale=scale,
                          blending="additive", rendering="mip", visible=visible,
                          multiscale=ml)
            if clim is not None:
                kwargs["contrast_limits"] = clim
            lyr = viewer.add_image(data, **kwargs)
            if clim is None:
                # Percentile window instead of napari's noise-floor auto.
                lyr.contrast_limits_range = (0.0, 65535.0)
                lyr.contrast_limits = percentile_limits(stack.data[:, c],
                                                        LO_PCT, HI_PCT)
            raw_layers.append(lyr)

        state["force_vol"] = fret["force_vol"] if fret is not None else None

        if fret is not None:
            print(f"  force range: {fret['force_vol'].min():.2f} – "
                  f"{fret['force_vol'].max():.2f} pN", flush=True)
            print(f"  eff   range: {fret['eff_vol'].min():.3f} – "
                  f"{fret['eff_vol'].max():.3f}", flush=True)
            raw_layers.append(viewer.add_image(
                pyramid_views(fret["force_vol"], levels),
                name="Force (pN)", colormap="turbo",
                scale=scale, contrast_limits=(0.0, 10.0), rendering="mip",
                blending="additive", visible=True, multiscale=ml,
            ))

        state["raw_layers"] = raw_layers

        dapi_raw = stack.data[:, role_to_c["nuclei"]] if "nuclei" in stack.roles else None
        cadh_raw = stack.data[:, role_to_c["substrate"]] if "substrate" in stack.roles else None
        clover_raw = stack.data[:, role_to_c["clover"]] if "clover" in stack.roles else None
        ruby_raw = stack.data[:, role_to_c["ruby"]] if "ruby" in stack.roles else None
        state["dapi_raw"] = dapi_raw
        state["cadh_raw"] = cadh_raw
        seed = state.get("seed_cadh_floor")
        if seed is not None:
            seed()
        state["clover_raw"] = clover_raw
        state["ruby_raw"] = ruby_raw
        state["clover_raw_native"] = clover_native
        state["ruby_raw_native"] = ruby_native
        state["raw_voxel_zyx_um"] = raw_voxel
        state["scale"] = scale
        state["levels"] = levels

        def _mask_layer(raw, name, cmap):
            if raw is None:
                return None
            return viewer.add_image(
                pyramid_views(np.zeros(raw.shape, dtype=np.uint8), levels),
                name=name, colormap=cmap, scale=scale,
                blending="additive", rendering="mip", visible=False,
                contrast_limits=(0, 1), multiscale=ml,
            )

        state["dapi_mask_layer"] = _mask_layer(dapi_raw, "DAPI Otsu mask", "cyan")
        state["cadh_mask_layer"] = _mask_layer(cadh_raw, "VE-Cadherin mask", "magenta")
        state["clover_mask_layer"] = _mask_layer(clover_raw, "Clover filter mask", "green")

        state["dapi_state"] = {"labels": None, "sizes": None}
        state["cadh_state"] = {
            "labels": None, "sizes": None,
            "centroids": None, "normals": None, "label_ids": None,
            "keep_mask": None,
        }
        state["clover_state"] = {"labels": None, "sizes": None, "method": None}
        # mesh_cache is cleared earlier, before the layer churn (see above).
        state["vessel_mesh_cache"] = {"verts": None, "faces": None}
        state["nuc_cache"] = {"labels_id": None, "centroids": None,
                              "normals": None, "lids": None}
        state["adhesion_proba"] = None
        state["clover_ruby_ratio"] = None
        for hist in state.get("hists", []):
            hist.clear()

        # Point every tab's save/load path field into a per-image folder next to
        # the image (e.g. images/WT/<stem>_saved/), so features for different
        # images don't collide. Widgets register themselves in
        # state["save_path_widgets"] as (widget, default_filename).
        save_dir = tif_path.parent / f"{tif_path.stem}_saved"
        try:
            save_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            print(f"  (could not create save folder {save_dir}: {exc})", flush=True)
        state["save_dir"] = save_dir
        state["image_name"] = tif_path.name
        for w, fn in state.get("save_path_widgets", []):
            try:
                # A "{stem}" placeholder in the default filename is filled with
                # the image name (e.g. the CSV -> <stem>_adhesion_measures.csv).
                w.value = str(save_dir / fn.replace("{stem}", tif_path.stem))
            except Exception:
                pass

        print("Load complete.", flush=True)
        if after:
            after()

    def _err(exc):
        print(f"Load error ({tif_path.name}): {exc}", flush=True)
        if after:
            after()

    w = _run()
    w.returned.connect(_done)
    w.errored.connect(_err)
    w.start()


def build(state, viewer, project_dir, params_path, lut_path) -> Container:
    project_dir_w = FileEdit(value=str(project_dir), mode="d", label="project_dir")
    image_w = ComboBox(choices=[], label="image")
    refresh_w = PushButton(text="Refresh images")
    pyramid_w = ComboBox(value=2, choices=[1, 2, 3], label="pyramid_levels")
    load_w = PushButton(text="Load")

    def _scan(*_):
        pd = Path(str(project_dir_w.value))
        imgs_dir = pd / "images"
        if not imgs_dir.exists():
            image_w.choices = []
            return
        tifs = sorted(imgs_dir.glob("*/*.tif"))
        image_w.choices = [str(p.relative_to(imgs_dir).as_posix()) for p in tifs]

    def _do_load(*_):
        pd = Path(str(project_dir_w.value))
        rel = image_w.value
        if not rel:
            print("No image selected.", flush=True)
            return
        tif_path = pd / "images" / rel
        _load_image(state, viewer, tif_path, params_path, lut_path,
                    int(pyramid_w.value))

    # "Load all saved" — after loading an image, trigger every tab's loader
    # (registered in state) so a session can be restored in one click. Each
    # loader uses its own path field's current value and reports/handles a
    # missing file itself, so features that were never saved are just skipped.
    load_all_w = PushButton(text="Load all saved features")

    # (label, state key) in dependency-friendly order.
    _LOADERS = [
        ("DAPI labels", "load_dapi_labels"),
        ("VE-Cadherin mask", "load_cadh"),
        ("Vessel mesh", "load_vessel_mesh"),
        ("Vessel curvature", "load_vessel_curvature"),
        ("Vessel path curvature", "load_vessel_path_curvature"),
        ("Adhesion proba", "load_proba"),
        ("Adhesion instances", "load_instances"),
        ("Clover/Ruby ratio", "load_ratio"),
        ("Fibroblast count", "load_fibroblast_count"),
    ]

    def _load_all(*_):
        if not state.get("raw_layers"):
            print("Load all: load an image first, then click again.", flush=True)
            return
        queue = [(label, key) for label, key in _LOADERS if state.get(key)]
        if not queue:
            print("Load all: no saved-feature loaders available.", flush=True)
            return

        # Load serially — each loader calls `after` when it finishes, which
        # starts the next. This caps peak RAM at one load at a time instead of
        # all at once, and `visible=False` keeps the heavy layers off the GPU
        # until you toggle them on (a burst of large 3D layers can overwhelm the
        # viewer). Missing files are reported by each loader and skipped.
        def _next(i):
            if i >= len(queue):
                print("Load all: finished. Loaded layers are hidden — toggle "
                      "visibility in the layer list as needed.", flush=True)
                return
            label, key = queue[i]
            fn = state.get(key)
            print(f"Load all: {i + 1}/{len(queue)} — {label} ...", flush=True)
            try:
                fn(after=lambda: _next(i + 1), visible=False)
            except Exception as exc:
                print(f"Load all: {label} failed: {exc}", flush=True)
                _next(i + 1)

        _next(0)

    load_all_w.changed.connect(_load_all)

    # "Run all" — run the whole analysis in order, each tab's compute + save,
    # waiting for each threaded step to finish before starting the next (same
    # serial `after`-callback chain as "Load all"). Each tab registers a
    # `run_*` step in state. Preconditions/failures in a step are reported and
    # the chain continues, so a missing model or channel won't stall it.
    run_all_w = PushButton(text="Run all (full pipeline)")
    _RUN_STEPS = [
        ("DAPI + components", "run_dapi"),
        ("Fibroblast count", "run_fibroblasts"),
        ("VE-Cadherin (Sauvola)", "run_cadherin"),
        ("Vessel mesh", "run_vessel_mesh"),
        ("Surface curvature", "run_vessel_curvature"),
        ("In-plane path curvature", "run_vessel_path_curvature"),
        ("Clover Otsu", "run_clover"),
        ("RF predict", "run_rf_predict"),
        ("RF post-process (split-merged)", "run_rf_postprocess"),
        ("Measure + export CSV", "run_measure"),
    ]

    def _run_all(*_):
        if not state.get("raw_layers"):
            print("Run all: load an image first, then click again.", flush=True)
            return
        queue = [(label, key) for label, key in _RUN_STEPS if state.get(key)]
        if not queue:
            print("Run all: no run steps available.", flush=True)
            return

        def _next(i):
            if i >= len(queue):
                print("Run all: finished — full pipeline complete.", flush=True)
                return
            label, key = queue[i]
            fn = state.get(key)
            print(f"Run all: {i + 1}/{len(queue)} — {label} ...", flush=True)
            try:
                fn(after=lambda: _next(i + 1))
            except Exception as exc:
                print(f"Run all: {label} failed: {exc}", flush=True)
                _next(i + 1)

        _next(0)

    run_all_w.changed.connect(_run_all)

    project_dir_w.changed.connect(_scan)
    refresh_w.changed.connect(_scan)
    load_w.changed.connect(_do_load)
    _scan()

    return Container(
        widgets=[project_dir_w, image_w, refresh_w, pyramid_w, load_w,
                 load_all_w, run_all_w],
        labels=True,
    )
