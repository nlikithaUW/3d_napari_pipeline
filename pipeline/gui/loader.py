"""Stack loader tab: project_dir + image picker + Load button."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from magicgui.widgets import CheckBox, ComboBox, Container, FileEdit, Label, PushButton
from napari.qt.threading import thread_worker

from .. import gpu
from ..fret import compute_fret_volume, load_force_lut, resample_fret_outputs
from ..io import fret_configured, load_stack, resolve_params, to_isotropic
from .state import RUN_STEPS
from .widgets import percentile_limits, pyramid_views, with_after


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
                after=None, fret_enabled: bool = True):
    print(f"Loading {tif_path.name} ...", flush=True)

    @thread_worker
    def _run():
        params, params_src, params_note = resolve_params(tif_path, params_path)
        print(f"  {params_note}", flush=True)
        stack_raw = load_stack(tif_path, params)
        role_to_c = {role: stack_raw.channels.index(name)
                     for role, name in stack_raw.roles.items()}
        fret_ok, why = fret_configured(params, stack_raw.roles, fret_enabled)
        if fret_ok:
            print(f"  FRET: computing ({why}).", flush=True)
            lut_eff, lut_force = load_force_lut(str(lut_path))
            fret = compute_fret_volume(stack_raw.data, role_to_c,
                                       params["fret"], lut_eff, lut_force)
        else:
            print(f"  FRET: skipped — {why}. No efficiency/force layers; the "
                  f"force_*_pN columns will export blank. The uncalibrated "
                  f"mClover/mRuby ratio is the tension proxy.", flush=True)
            fret = None
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
        if fret is not None:
            fret = resample_fret_outputs(fret, stack_raw.voxel_zyx_um,
                                         stack.voxel_zyx_um[0])
        # Drop the raw uint16 stack — closure would otherwise retain it.
        del stack_raw
        params_used = {"source": str(params_src), "fret": params.get("fret"),
                       "fret_enabled": bool(fret_enabled)}
        return (raw_shape, raw_nbytes, raw_voxel, stack, role_to_c, fret,
                clover_native, ruby_native, params_used)

    def _done(result):
        (raw_shape, raw_nbytes, raw_voxel, stack, role_to_c, fret,
         clover_native, ruby_native, params_used) = result
        # Into the settings JSON beside every export: which params file this
        # image was analysed with, and the FRET constants actually applied.
        state.setdefault("run_settings", {})["params"] = params_used

        # Give the device its room back *before* the layer churn below, not
        # after. `with_after` drains at the end of every step, but this handler
        # is itself the step: it removes one image's worth of napari layers and
        # allocates the next image's, and those GL textures come out of the same
        # VRAM that cupy's pool is sitting on from the previous image's filters.
        # A GL allocation that fails raises nothing here — it surfaces later as
        # an `access violation` inside `glDrawArrays` ([P48]). The free/total
        # line is printed so the next such fault has a number beside it.
        gpu.free_pools()
        mem = gpu.device_memory_mb()
        if mem is not None:
            print(f"  vram  {mem[0]:.0f} / {mem[1]:.0f} MB free "
                  f"(after draining the cupy pools)", flush=True)
        print(f"  raw   shape={raw_shape} voxel(z,y,x)={raw_voxel} µm "
              f"({raw_nbytes/1e6:.1f} MB)", flush=True)
        print(f"  iso   shape={stack.data.shape} dtype={stack.data.dtype} "
              f"voxel(z,y,x)={stack.voxel_zyx_um} µm "
              f"({stack.data.nbytes/1e6:.1f} MB)", flush=True)

        # Clear the cell-mesh cache on load: it holds the *previous* image's
        # geometry in µm, and Refine/Measure would otherwise read it against the
        # new image. (This also used to stop a `current_step` hook re-creating
        # the "Cell mesh (slice)" overlay at the new image's z; that overlay is
        # gone — see [P47] — but the cache still has to be cleared.)
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
        # `seed_cadh_floor` is NOT called here any more. It overwrote the
        # Cadherin tab's `floor_k` on every load, so the tab default never
        # reached the build and a headless run silently used a different value
        # than the one recorded as the default (P26). It is a button on the
        # Cadherin tab now — press it if you want the suggestion.
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

        # "labels"/"sizes" are the unsplit DAPI connected components; "nuclei"
        # and friends are the watershed-split nuclei + their per-nucleus
        # fibroblast class (P04). Two different objects, one dict, replaced
        # wholesale per image so neither can leak across images.
        state["dapi_state"] = {"labels": None, "sizes": None,
                               "nuclei": None, "nuclei_sizes": None,
                               "nuclei_ids": None, "nuclei_is_fib": None,
                               "nuclei_min_size_vox": None}
        state["cadh_state"] = {
            "labels": None, "sizes": None,
            "centroids": None, "normals": None, "label_ids": None,
            "keep_mask": None,
        }
        state["clover_state"] = {"labels": None, "sizes": None, "method": None}
        # mesh_cache is cleared earlier, before the layer churn (see above).
        state["vessel_mesh_cache"] = {"verts": None, "faces": None}
        # **The graph, the shear and the solved lumen go with it.** Only the
        # mesh cache was cleared here, so loading a second image left image 1's
        # solved lumen field in state and "Compute 3D centerline graph" then
        # traced image 1's centerline and reported it as image 2's. That is
        # [P38](../../docs/ISSUES.md#p38)'s failure mode — one build's geometry
        # surviving into another's — on the path P38 did not cover: it fixed the
        # failed-*build* path, not the image-*load* path.
        state["vessel_lumen_field"] = None
        state["vessel_graph"] = None
        state["vessel_shear"] = None
        state["vessel_build_id"] = 0
        state["nuc_cache"] = {"labels_id": None, "centroids": None,
                              "normals": None, "lids": None}
        state["adhesion_proba"] = None
        state["clover_ruby_ratio"] = None
        # Must be cleared per image (P25). `run_fibroblasts` chains its save
        # through `call_all`, which fires even when the count bails out early —
        # so a stale count from the previous image would be written into this
        # image's fibroblast_count.npz and, from there, into its CSV.
        state["fibroblast_count"] = None
        for hist in state.get("hists", []):
            hist.clear()
        # The Vessel mesh tab's four status lines describe caches this handler
        # just emptied.
        #
        # **Guarded, and the reason is severe.** These callbacks only repaint
        # text, but they run *before* the save-path block below, so an exception
        # in one aborts `_done` and every tab keeps the previous image's path —
        # or, on the first image, its bare relative default, which resolves to
        # the current working directory. Every artifact of the run then lands in
        # the CWD under one name, each image overwriting the last, with the only
        # symptom a traceback that the batch scrolls past. Measured 2026-09-02
        # on a 6-image `--steps cadherin` run: two images written to
        # `./cadherin_labels.npz` before it was caught.
        #
        # The trigger was `batch.py` discarding what `build()` returns, which
        # let Qt collect the labels these closures write to ("wrapped C/C++
        # object of type QLabel has been deleted"). That is fixed there too, but
        # the ordering hazard is the real defect: nothing cosmetic may stand
        # between a load and the save paths that load implies.
        for fn in state.get("vessel_status_refresh", []):
            try:
                fn()
            except Exception as exc:
                print(f"  (status refresh failed, continuing: "
                      f"{type(exc).__name__}: {exc})", flush=True)

        # Point every tab's save/load path field into a per-image folder next to
        # the image (e.g. images/WT/<stem>_saved/), so features for different
        # images don't collide. Widgets register themselves in
        # state["save_path_widgets"] as (widget, default_filename).
        save_dir = tif_path.parent / f"{tif_path.stem}_saved"
        try:
            save_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            print(f"  (could not create save folder {save_dir}: {exc})", flush=True)
        state["image_name"] = tif_path.name
        for w, fn in state.get("save_path_widgets", []):
            try:
                # A "{stem}" placeholder in the default filename is filled with
                # the image name (e.g. the CSV -> <stem>_adhesion_measures.csv).
                w.value = str(save_dir / fn.replace("{stem}", tif_path.stem))
            except Exception:
                pass

        print("Load complete.", flush=True)

    def _err(exc):
        print(f"Load error ({tif_path.name}): {exc}", flush=True)

    w = _run()
    w.returned.connect(with_after(_done, after, "Load"))
    w.errored.connect(with_after(_err, after, "Load"))
    w.start()


def build(state, viewer, project_dir, params_path, lut_path) -> Container:
    project_dir_w = FileEdit(value=str(project_dir), mode="d", label="project_dir")
    ALL = "(all)"

    def _experiments(*_):
        # A callable so every reset_choices() re-reads the disk (a folder the
        # Convert tab just created appears without restarting).
        imgs = Path(str(project_dir_w.value)) / "images"
        dirs = sorted(p.name for p in imgs.iterdir() if p.is_dir()) if imgs.exists() else []
        return [ALL] + dirs

    experiment_w = ComboBox(choices=_experiments, value=ALL, label="experiment")
    # Backed by a callable for the same reason: a list assigned to `.choices`
    # is discarded by the next reset_choices(), emptying the picker silently.
    image_choices: list = []
    image_w = ComboBox(choices=lambda *_: image_choices, label="image")
    refresh_w = PushButton(text="Refresh images")
    pyramid_w = ComboBox(value=2, choices=[1, 2, 3], label="pyramid_levels")
    # Off switch only: FRET still needs a fret block in the image's experiment
    # params.json (io.resolve_params). Applies at the next Load. Headless: --no-fret.
    fret_w = CheckBox(value=True, text="compute FRET")
    fret_status_w = Label(value="")
    load_w = PushButton(text="Load")

    def _fret_status(*_):
        rel = image_w.value
        if not rel:
            fret_status_w.value = ""
            return
        tif = Path(str(project_dir_w.value)) / "images" / rel
        try:
            params, src, _ = resolve_params(tif, params_path)
        except Exception as exc:  # noqa: BLE001
            fret_status_w.value = f"params unreadable: {exc}"
            return
        if not params.get("fret"):
            fret_status_w.value = (f"FRET unavailable — no calibration in "
                                   f"{tif.parent.name}/params.json")
        elif fret_w.value:
            fret_status_w.value = f"FRET will be computed ({src.parent.name}/{src.name})"
        else:
            fret_status_w.value = "FRET available but switched off"

    def _scan(*_):
        pd = Path(str(project_dir_w.value))
        imgs_dir = pd / "images"
        experiment_w.reset_choices()
        exp = experiment_w.value
        tifs = (sorted(imgs_dir.glob("*/*.tif" if exp in (None, ALL) else f"{exp}/*.tif"))
                if imgs_dir.exists() else [])
        # The value stays "<experiment>/<file>.tif" (what `_do_load` joins onto
        # images/); only the label drops the folder once one is picked.
        image_choices[:] = [
            (p.name if exp not in (None, ALL) else p.relative_to(imgs_dir).as_posix(),
             p.relative_to(imgs_dir).as_posix())
            for p in tifs]
        image_w.reset_choices()
        _fret_status()

    def _do_load(*_):
        pd = Path(str(project_dir_w.value))
        rel = image_w.value
        if not rel:
            print("No image selected.", flush=True)
            return
        tif_path = pd / "images" / rel
        _load_image(state, viewer, tif_path, params_path, lut_path,
                    int(pyramid_w.value), fret_enabled=bool(fret_w.value))

    # "Load all saved" and "Run all" are the same machine: walk an ordered
    # (label, state key) list, calling each registered callable and letting its
    # `after` hand-off start the next. Serial, so peak RAM is one step's worth
    # rather than all of them, and a step that fails is reported and skipped
    # instead of stalling the chain.
    def _chain(what: str, steps: list, empty_msg: str, done_msg: str, **kwargs):
        if not state.get("raw_layers"):
            print(f"{what}: load an image first, then click again.", flush=True)
            return
        queue = [(label, key) for label, key in steps if state.get(key)]
        if not queue:
            print(f"{what}: {empty_msg}", flush=True)
            return

        def _next(i):
            if i >= len(queue):
                print(f"{what}: {done_msg}", flush=True)
                return
            label, key = queue[i]
            fn = state.get(key)
            print(f"{what}: {i + 1}/{len(queue)} — {label} ...", flush=True)
            try:
                fn(after=lambda: _next(i + 1), **kwargs)
            except Exception as exc:
                print(f"{what}: {label} failed: {exc}", flush=True)
                _next(i + 1)

        _next(0)

    # ---- "Load all saved": every tab's loader, so a session restores in one
    # click. Each loader uses its own path field's current value and opens with
    # `widgets.skip_if_missing`, so a feature that was never saved prints one
    # skip line and hands off. That guard is shared rather than re-written per
    # tab: it existed in exactly one loader of twelve, and the other eleven
    # turned a missing optional artifact into a full traceback and a napari
    # error popup -- the chain survived it, but it does not read that way.
    # (label, state key) in dependency-friendly order.
    _LOADERS = [
        ("DAPI labels", "load_dapi_labels"),
        ("Split nuclei + classes", "load_nuclei_split"),
        ("VE-Cadherin mask", "load_cadh"),
        ("Vessel mesh", "load_vessel_mesh"),
        ("Vessel curvature", "load_vessel_curvature"),
        ("Vessel centerline graph", "load_vessel_path_curvature"),
        ("Vessel wall shear", "load_vessel_shear"),
        ("Adhesion proba", "load_proba"),
        ("Adhesion instances", "load_instances"),
        ("Clover/Ruby ratio", "load_ratio"),
        ("Fibroblast count", "load_fibroblast_count"),
    ]
    load_all_w = PushButton(text="Load all saved features")
    load_all_w.changed.connect(lambda *_: _chain(
        "Load all", _LOADERS,
        "no saved-feature loaders available.",
        "finished. Loaded layers are hidden — toggle visibility in the layer "
        "list as needed.",
        # `visible=False` keeps the heavy layers off the GPU until you toggle
        # them on; a burst of large 3D layers can overwhelm the viewer.
        visible=False))

    # ---- "Run all": the whole analysis in order. The step table is
    # `gui/state.RUN_STEPS`, shared with the headless `pipeline.batch` runner so
    # a step cannot be added to one and forgotten in the other.
    run_all_w = PushButton(text="Run all (full pipeline)")
    run_all_w.changed.connect(lambda *_: _chain(
        "Run all", [(label, key) for _n, label, key in RUN_STEPS],
        "no run steps available.",
        "finished — full pipeline complete."))

    project_dir_w.changed.connect(_scan)
    experiment_w.changed.connect(_scan)
    image_w.changed.connect(_fret_status)
    fret_w.changed.connect(_fret_status)
    refresh_w.changed.connect(_scan)
    state["refresh_image_list"] = _scan   # the Convert tab calls this after writing
    load_w.changed.connect(_do_load)
    _scan()

    return Container(
        widgets=[project_dir_w, experiment_w, image_w, refresh_w, pyramid_w,
                 fret_w, fret_status_w, load_w,
                 load_all_w, run_all_w],
        labels=True,
    )
