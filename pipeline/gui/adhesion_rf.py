"""Adhesion RF tab.

Workflow:
  1. Create label layer → paint adhesions (1) and background (2).
  2. Train RF.
  3. Run RF predict  — slow; stores probability map in state["adhesion_proba"].
  4. Apply post-process — fast; threshold + DAPI/mesh gate + instance labeling.
     Re-run freely without re-running RF.
  5. Save / load model.
"""

from __future__ import annotations

import numpy as np
from magicgui.widgets import (
    CheckBox,
    Container,
    FileEdit,
    FloatSlider,
    Label,
    PushButton,
    Slider,
    SpinBox,
)
from napari.qt.threading import thread_worker

from ..denoise import auto_min_size_otsu, filter_by_size, label_and_sizes, otsu_mask
from ..rf import AdhesionRF, instance_labels_3d
from .widgets import debounce, pyramid_views


def build(state, viewer) -> Container:

    # ------------------------------------------------------------------ #
    # Section 1 — Create annotation layer                                  #
    # ------------------------------------------------------------------ #

    hint_w = Label(value="Paint: 1 = adhesion,  2 = background")
    create_label_w = PushButton(text="Create label layer")

    def _create_label(*_):
        clover_raw = state.get("clover_raw")
        if clover_raw is None:
            print("Load a stack with a 'clover' channel first.", flush=True)
            return
        name = "Adhesion labels"
        if name in [lay.name for lay in viewer.layers]:
            viewer.layers.selection.active = viewer.layers[name]
            print(f"'{name}' layer already exists — made active.", flush=True)
            return
        ann = np.zeros(clover_raw.shape, dtype=np.uint8)
        scale = state.get("scale", (1.0, 1.0, 1.0))
        viewer.add_labels(ann, name=name, scale=scale)
        print("'Adhesion labels' created. Select it and paint with the brush tool.",
              flush=True)

    create_label_w.changed.connect(_create_label)

    # ------------------------------------------------------------------ #
    # Section 2 — Train                                                    #
    # ------------------------------------------------------------------ #

    n_estimators_w = SpinBox(value=100, min=50, max=500, label="n_estimators")
    train_w = PushButton(text="Train RF")

    def _train(*_):
        clover_raw = state.get("clover_raw")
        if clover_raw is None:
            print("Load a stack first.", flush=True)
            return
        name = "Adhesion labels"
        if name not in [lay.name for lay in viewer.layers]:
            print(f"No '{name}' layer — create it and paint first.", flush=True)
            return

        ann_layer = viewer.layers[name]
        if getattr(ann_layer, "multiscale", False):
            ann_vol = np.asarray(ann_layer.data[0]).astype(np.uint8)
        else:
            ann_vol = np.asarray(ann_layer.data).astype(np.uint8)

        if not np.any(ann_vol > 0):
            print("No annotations found. Paint first.", flush=True)
            return

        n_est = int(n_estimators_w.value)
        voxel_um = state.get("scale", (1.0, 1.0, 1.0))
        vol_f32 = clover_raw.astype(np.float32)

        @thread_worker
        def _run():
            rf = AdhesionRF(n_estimators=n_est)
            rf.train(vol_f32, ann_vol, voxel_um)
            return rf

        def _done(rf):
            state["adhesion_rf"] = rf
            print("Training complete.", flush=True)

        def _err(exc):
            print(f"Training error: {exc}", flush=True)

        print(f"Training RF (n_estimators={n_est}) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    train_w.changed.connect(_train)

    # ------------------------------------------------------------------ #
    # Section 3 — RF predict (slow; run once, stores proba map)           #
    # ------------------------------------------------------------------ #

    gate_by_clover_w = CheckBox(value=True, label="gate_by_clover")
    clover_dilation_um_w = FloatSlider(value=1.0, min=0.0, max=10.0,
                                       label="clover_dilation_um")
    clover_dilation_um_w.visible = True

    def _on_gate_clover(val):
        clover_dilation_um_w.visible = bool(val)

    gate_by_clover_w.changed.connect(_on_gate_clover)

    predict_w = PushButton(text="Run RF predict")

    def _predict(*_, after=None):
        rf: AdhesionRF | None = state.get("adhesion_rf")
        if rf is None or not rf.is_trained:
            print("No trained RF. Train or load a model first.", flush=True)
            if after:
                after()
            return
        clover_raw = state.get("clover_raw")
        if clover_raw is None:
            print("Load a stack first.", flush=True)
            if after:
                after()
            return

        voxel_um = state.get("scale", (1.0, 1.0, 1.0))
        vol_f32 = clover_raw.astype(np.float32)
        levels = state.get("levels", 1)
        scale = voxel_um
        gate_clover = bool(gate_by_clover_w.value)
        dilation_um = float(clover_dilation_um_w.value)
        state.setdefault("run_settings", {})["rf_predict"] = {
            "gate_by_clover": gate_clover, "clover_dilation_um": dilation_um}

        def _start_rf():
            @thread_worker
            def _run():
                gate_indices = None
                if gate_clover:
                    clover_mask_layer = state.get("clover_mask_layer")
                    if clover_mask_layer is not None:
                        ldata = clover_mask_layer.data
                        if getattr(clover_mask_layer, "multiscale", False):
                            mask = np.asarray(ldata[0]).astype(bool)
                        else:
                            mask = np.asarray(ldata).astype(bool)

                        if mask.shape != vol_f32.shape:
                            # Stale clover mask (e.g. left over from a previously
                            # loaded image of a different size) — its coordinates
                            # would be out of bounds for this volume.
                            print(f"gate_by_clover: clover mask {mask.shape} != "
                                  f"volume {vol_f32.shape} (stale mask) — predicting "
                                  "the full volume. Recompute the clover mask (or "
                                  "reload the image) to gate.", flush=True)
                        else:
                            if dilation_um > 0 and mask.any() and not mask.all():
                                vox = tuple(float(v) for v in voxel_um) \
                                      if hasattr(voxel_um, "__len__") \
                                      else (float(voxel_um),) * 3
                                rz = max(1, int(np.ceil(dilation_um / vox[0])))
                                ry = max(1, int(np.ceil(dilation_um / vox[1])))
                                rx = max(1, int(np.ceil(dilation_um / vox[2])))
                                Z, Y, X = np.ogrid[-rz:rz+1, -ry:ry+1, -rx:rx+1]
                                struct = (Z/rz)**2 + (Y/ry)**2 + (X/rx)**2 <= 1.0
                                print(f"Clover dilation: {dilation_um} µm "
                                      f"(SE {2*rz+1}×{2*ry+1}×{2*rx+1}) ...",
                                      flush=True)
                                from .. import gpu as _gpu
                                if _gpu.HAS_GPU:
                                    import cupy as cp
                                    import cupyx.scipy.ndimage as _cndi
                                    mask = _cndi.binary_dilation(
                                        cp.asarray(mask),
                                        structure=cp.asarray(struct),
                                    ).get()
                                else:
                                    from scipy.ndimage import binary_dilation
                                    mask = binary_dilation(mask, structure=struct)

                            gate_indices = np.argwhere(mask).astype(np.int32)
                            pct = 100.0 * len(gate_indices) / mask.size
                            print(f"Clover gate: {len(gate_indices):,} / "
                                  f"{mask.size:,} voxels ({pct:.1f}%)", flush=True)
                    else:
                        print("gate_by_clover: no clover mask found — "
                              "predicting full volume.", flush=True)

                proba = rf.predict_proba(vol_f32, voxel_um, gate_indices=gate_indices)
                return proba

            def _done(proba):
                state["adhesion_proba"] = proba
                name = "Adhesion proba"
                lyr_data = pyramid_views(proba, levels)
                ml = levels > 1
                if name in [lay.name for lay in viewer.layers]:
                    viewer.layers.remove(name)
                viewer.add_image(
                    lyr_data, name=name, scale=scale,
                    colormap="magma", blending="additive",
                    contrast_limits=(0, 255),
                    rendering="mip", multiscale=ml,
                )
                print("Probability map stored. Adjust post-process and hit Apply.",
                      flush=True)
                if after:
                    after()

            def _err(exc):
                print(f"Predict error: {exc}", flush=True)
                if after:
                    after()

            mode = f"gated by clover (dilation={dilation_um} µm)" if gate_clover else "full volume"
            print(f"Running RF predict ({mode}) ...", flush=True)
            w = _run()
            w.returned.connect(_done)
            w.errored.connect(_err)
            w.start()

        if gate_clover and state["clover_state"]["labels"] is None:
            print("gate_by_clover: clover mask not computed — running Otsu first ...",
                  flush=True)

            @thread_worker
            def _run_otsu():
                mask = otsu_mask(clover_raw)
                labels, sizes = label_and_sizes(mask)
                suggested = auto_min_size_otsu(sizes)
                filtered = filter_by_size(labels, sizes, suggested)
                return labels, sizes, suggested, filtered

            def _otsu_done(result):
                labels, sizes, suggested, filtered = result
                state["clover_state"]["labels"] = labels
                state["clover_state"]["sizes"] = sizes
                state["clover_state"]["method"] = "Otsu"
                clover_mask_layer = state.get("clover_mask_layer")
                if clover_mask_layer is not None:
                    clover_mask_layer.data = pyramid_views(filtered, levels)
                print(f"Clover Otsu done: {sizes.size} components, "
                      f"min_size={suggested}.", flush=True)
                _start_rf()

            def _otsu_err(exc):
                print(f"Clover Otsu error: {exc}", flush=True)
                if after:
                    after()

            w = _run_otsu()
            w.returned.connect(_otsu_done)
            w.errored.connect(_otsu_err)
            w.start()
        else:
            _start_rf()

    predict_w.changed.connect(_predict)

    # ------------------------------------------------------------------ #
    # Section 3b — Save / load the probability map                        #
    # ------------------------------------------------------------------ #
    # RF predict is the slow step and its output ("Adhesion proba") lives only
    # in memory, so a viewer restart otherwise means recomputing it. Persist the
    # proba volume to a .npy and reload it to resume straight at post-processing.

    proba_path_w = FileEdit(value="adhesion_proba.npy", label="proba_path",
                            mode="w")
    save_proba_w = PushButton(text="Save proba")
    load_proba_w = PushButton(text="Load proba")

    def _save_proba(*_):
        proba = state.get("adhesion_proba")
        if proba is None:
            print("No probability map to save. Run RF predict first.", flush=True)
            return
        path = str(proba_path_w.value)
        np.save(path, proba)
        print(f"Saved probability map -> {path} "
              f"({proba.nbytes / 1e6:.1f} MB)", flush=True)

    def _load_proba(*_, after=None, visible=True):
        path = str(proba_path_w.value)
        levels = state.get("levels", 1)
        scale = state.get("scale", (1.0, 1.0, 1.0))

        @thread_worker
        def _run():
            return np.load(path)

        def _done(proba):
            clover_raw = state.get("clover_raw")
            if clover_raw is not None and proba.shape != clover_raw.shape:
                print(f"Warning: proba shape {proba.shape} != current stack "
                      f"{clover_raw.shape}. Load the matching image first.",
                      flush=True)
            state["adhesion_proba"] = proba
            name = "Adhesion proba"
            ml = levels > 1
            if name in [lay.name for lay in viewer.layers]:
                viewer.layers.remove(name)
            lyr = viewer.add_image(
                pyramid_views(proba, levels), name=name, scale=scale,
                colormap="magma", blending="additive",
                contrast_limits=(0, 255), rendering="mip", multiscale=ml,
            )
            lyr.visible = visible
            print("Probability map loaded. Adjust post-process and hit Apply.",
                  flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Load proba error: {exc}", flush=True)
            if after:
                after()

        print(f"Loading probability map from {path} ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    save_proba_w.changed.connect(_save_proba)
    load_proba_w.changed.connect(_load_proba)
    state["load_proba"] = _load_proba

    # ------------------------------------------------------------------ #
    # Section 4 — Post-process (fast; re-run freely)                      #
    # ------------------------------------------------------------------ #

    prob_threshold_w = FloatSlider(value=0.5, min=0.0, max=1.0,
                                   label="prob_threshold")
    gate_by_dapi_w = CheckBox(value=False, label="gate_by_dapi")
    d_dapi_um_w = FloatSlider(value=15.0, min=0.0, max=50.0,
                               label="d_dapi_um")
    d_dapi_um_w.visible = False

    gate_by_mesh_w = CheckBox(value=False, label="gate_by_mesh")
    d_mesh_um_w = FloatSlider(value=10.0, min=0.0, max=50.0,
                               label="d_mesh_um")
    d_mesh_um_w.visible = False

    min_size_vox_w = Slider(value=100, min=0, max=50000, label="min_size_vox")
    # Default standard settings: threshold 0.5, min_size 100, split ON — gives
    # adhesion axis ratios centred in the expected 2-3 range for this cell type.
    split_merged_w = CheckBox(value=True, label="split_merged")
    split_min_size_w = Slider(value=1000, min=0, max=10000,
                               label="split_min_size_vox")
    split_min_size_w.visible = True

    apply_w = PushButton(text="Apply post-process")

    def _on_gate_dapi(val):  d_dapi_um_w.visible = bool(val)
    def _on_gate_mesh(val):  d_mesh_um_w.visible = bool(val)
    def _on_split(val):      split_min_size_w.visible = bool(val)

    gate_by_dapi_w.changed.connect(_on_gate_dapi)
    gate_by_mesh_w.changed.connect(_on_gate_mesh)
    split_merged_w.changed.connect(_on_split)

    def _apply(*_, after=None):
        proba_vol = state.get("adhesion_proba")
        if proba_vol is None:
            print("No probability map. Run RF predict first.", flush=True)
            if after:
                after()
            return

        thresh = float(prob_threshold_w.value)
        min_sz = int(min_size_vox_w.value)
        split = bool(split_merged_w.value)
        split_min = int(split_min_size_w.value)
        gate_dapi = bool(gate_by_dapi_w.value)
        gate_mesh = bool(gate_by_mesh_w.value)
        d_dapi = float(d_dapi_um_w.value)
        d_mesh = float(d_mesh_um_w.value)
        voxel_um = state.get("scale", (1.0, 1.0, 1.0))
        levels = state.get("levels", 1)
        scale = voxel_um
        state.setdefault("run_settings", {})["adhesion_postprocess"] = {
            "prob_threshold": thresh, "min_size_vox": min_sz,
            "split_merged": split, "split_min_size_vox": split_min,
            "gate_by_dapi": gate_dapi, "d_dapi_um": d_dapi,
            "gate_by_mesh": gate_mesh, "d_mesh_um": d_mesh}

        @thread_worker
        def _run():
            # Unscale uint8 → float32 [0, 1] for thresholding / gate masking.
            masked = proba_vol.astype(np.float32) * (1.0 / 255.0)

            if gate_dapi:
                dapi_labels = (state.get("dapi_state") or {}).get("labels")
                if dapi_labels is not None:
                    from scipy.ndimage import distance_transform_edt
                    sampling = tuple(float(v) for v in voxel_um) \
                               if hasattr(voxel_um, "__len__") \
                               else (float(voxel_um),) * 3
                    dist = distance_transform_edt(dapi_labels == 0,
                                                  sampling=sampling)
                    masked[dist > d_dapi] = 0.0
                    print(f"DAPI gate: kept {float((dist <= d_dapi).mean())*100:.1f}% "
                          f"of voxels", flush=True)
                else:
                    print("DAPI gate requested but no DAPI labels — skipping.",
                          flush=True)

            if gate_mesh:
                mc = state.get("mesh_cache") or {}
                if mc.get("offset_centroids") is not None and \
                        len(mc["offset_centroids"]):
                    from ..cell_mesh import signed_distance_to_mesh_vertices
                    s = np.asarray(voxel_um, dtype=np.float64) \
                        if hasattr(voxel_um, "__len__") \
                        else np.full(3, float(voxel_um))
                    idx = np.argwhere(masked > 0)
                    if len(idx):
                        sd = signed_distance_to_mesh_vertices(
                            idx.astype(np.float64) * s[None, :],
                            mc["offset_centroids"], mc["outward"],
                            tree=mc.get("kdtree"),
                        )
                        drop = idx[~((sd > -d_mesh) & (sd <= d_mesh))]
                        if len(drop):
                            masked[drop[:, 0], drop[:, 1], drop[:, 2]] = 0.0
                else:
                    print("Mesh gate requested but no cell mesh — skipping.",
                          flush=True)

            binary = (masked >= thresh).astype(np.uint8)
            labels = instance_labels_3d(binary, min_sz, split, split_min)
            n_inst = int(labels.max())
            if n_inst > 0:
                szs = np.bincount(labels.ravel())[1:]
                print(f"Adhesion instances: {n_inst}  "
                      f"[mean={szs.mean():.0f} min={szs.min()} "
                      f"max={szs.max()} vox]", flush=True)
            else:
                print("No instances — lower prob_threshold or min_size_vox.",
                      flush=True)
            return labels

        def _done(labels):
            name = "Adhesion instances"
            if name in [lay.name for lay in viewer.layers]:
                viewer.layers.remove(name)
            # Single-scale Labels: napari's 3D multiscale-Labels support is
            # flaky and can add a layer with no GPU visual, which then breaks
            # the canvas layer-reorder on the next visibility toggle.
            viewer.add_labels(labels, name=name, scale=scale)
            if after:
                after()

        def _err(exc):
            print(f"Apply error: {exc}", flush=True)
            if after:
                after()

        print("Applying post-process ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    apply_w.changed.connect(_apply)

    # Save / load the "Adhesion instances" label layer (skip re-running Apply).
    inst_path_w = FileEdit(value="adhesion_instances.npz", label="inst_path",
                           mode="w")
    save_inst_w = PushButton(text="Save instances")
    load_inst_w = PushButton(text="Load instances")
    INST_LAYER = "Adhesion instances"

    def _save_inst(*_):
        if INST_LAYER not in [l.name for l in viewer.layers]:
            print("No 'Adhesion instances' layer — run Apply post-process first.",
                  flush=True)
            return
        lyr = viewer.layers[INST_LAYER]
        d = lyr.data
        vol = np.asarray(d[0] if getattr(lyr, "multiscale", False) else d)
        vol = vol.astype(np.int32)
        path = str(inst_path_w.value)
        np.savez_compressed(path, instances=vol)
        print(f"Saved instances -> {path} ({int(vol.max())} adhesions).",
              flush=True)

    def _load_inst(*_, after=None, visible=True):
        path = str(inst_path_w.value)
        levels = state.get("levels", 1)
        scale = state.get("scale", (1.0, 1.0, 1.0))

        @thread_worker
        def _run():
            return np.load(path)["instances"].astype(np.int32)

        def _done(vol):
            if INST_LAYER in [l.name for l in viewer.layers]:
                viewer.layers.remove(INST_LAYER)
            # Single-scale Labels; add visible so the GPU visual is built, then
            # hide (a hidden multiscale Labels layer can end up with no visual).
            lyr = viewer.add_labels(vol, name=INST_LAYER, scale=scale)
            lyr.visible = visible
            print(f"Loaded instances <- {path} ({int(vol.max())} adhesions).",
                  flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Load instances error: {exc}", flush=True)
            if after:
                after()

        print(f"Loading instances from {path} ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    save_inst_w.changed.connect(_save_inst)
    load_inst_w.changed.connect(_load_inst)
    state["load_instances"] = _load_inst

    # Chainable "Run all" steps: RF predict then save proba; post-process (with
    # split_merged on) then save instances; each calls `after` when done.
    state["run_rf_predict"] = lambda after=None: _predict(
        after=(lambda: (_save_proba(), after() if after else None)))

    def _run_rf_postprocess(after=None):
        split_merged_w.value = True
        _apply(after=(lambda: (_save_inst(), after() if after else None)))
    state["run_rf_postprocess"] = _run_rf_postprocess

    # ------------------------------------------------------------------ #
    # Section 4b — Clover/Ruby ratio inside adhesions                      #
    # ------------------------------------------------------------------ #
    # Derived layer: clover / ruby per voxel, but only inside the "Adhesion
    # instances" mask (0 outside, and 0 where ruby == 0). Also reports a
    # per-adhesion mean ratio, averaged over that adhesion's ruby > 0 voxels,
    # stored in state["clover_ruby_ratio"] for downstream use (e.g. Measure).

    ratio_btn = PushButton(text="Clover/Ruby ratio (in adhesions)")

    def _ratio(*_):
        clover_native = state.get("clover_raw_native")
        ruby_native = state.get("ruby_raw_native")
        if clover_native is None or ruby_native is None:
            print("Clover/Ruby: need native clover and ruby channels — reload "
                  "the image (native channels are captured at load).", flush=True)
            return
        if "Adhesion instances" not in [lay.name for lay in viewer.layers]:
            print("Clover/Ruby: run RF predict + Apply post-process first.",
                  flush=True)
            return

        lyr = viewer.layers["Adhesion instances"]
        adh = lyr.data
        labels_iso = adh[0] if getattr(lyr, "multiscale", False) else adh
        levels = state.get("levels", 1)
        iso_scale = state.get("scale", (1.0, 1.0, 1.0))
        raw_voxel = tuple(state.get("raw_voxel_zyx_um", iso_scale))

        @thread_worker
        def _run():
            labels_iso_a = np.asarray(labels_iso).astype(np.int32)
            clover = clover_native.astype(np.float32)
            ruby = ruby_native.astype(np.float32)
            if clover.shape[1:] != labels_iso_a.shape[1:]:
                raise ValueError(
                    f"native clover XY {clover.shape[1:]} != instances XY "
                    f"{labels_iso_a.shape[1:]} — reload the image and re-Apply."
                )
            # Map the isotropic adhesion labels DOWN to the native Z planes
            # (nearest neighbour): native plane n sits at depth n*raw_z, matching
            # isotropic index round(n * raw_z / iso_z). Labels are nearest, never
            # interpolated. XY is unchanged by isotropic resampling.
            Zn = clover.shape[0]
            Ziso = labels_iso_a.shape[0]
            zf = float(raw_voxel[0]) / float(iso_scale[0])
            iso_idx = np.clip(np.rint(np.arange(Zn) * zf).astype(np.int64),
                              0, Ziso - 1)
            labels = labels_iso_a[iso_idx]           # (Zn, Y, X)

            # Per-voxel ratio on NATIVE Z, only inside adhesions and where
            # ruby > 0 (0 elsewhere, matching the chosen zero policy).
            inside = (labels > 0) & (ruby > 0)
            ratio = np.zeros(labels.shape, dtype=np.float32)
            ratio[inside] = clover[inside] / ruby[inside]

            # Per-adhesion mean over that adhesion's valid (ruby > 0) native voxels.
            lv = labels[inside]
            rv = ratio[inside]
            n_lab = int(labels.max())
            sums = np.bincount(lv, weights=rv, minlength=n_lab + 1)
            counts = np.bincount(lv, minlength=n_lab + 1)
            valid = counts[1:] > 0
            ids = np.arange(1, n_lab + 1)[valid]
            mean_ratio = (sums[1:][valid] / counts[1:][valid]).astype(np.float32)
            n_vox = counts[1:][valid].astype(np.int64)
            return ratio, ids, mean_ratio, n_vox

        def _done(result):
            ratio, ids, mean_ratio, n_vox = result
            state["clover_ruby_ratio"] = {
                "label": ids, "ratio_mean": mean_ratio, "n_vox": n_vox,
            }
            name = "Clover/Ruby (adhesions)"
            ml = levels > 1
            nz = ratio[ratio > 0]
            vmax = float(np.percentile(nz, 99)) if nz.size else 1.0
            if name in [lay.name for lay in viewer.layers]:
                viewer.layers.remove(name)
            # Native-Z layer: display at the native voxel size so it stays
            # physically aligned with the isotropic 3D scene.
            viewer.add_image(
                pyramid_views(ratio, levels), name=name, scale=raw_voxel,
                colormap="viridis", blending="additive", rendering="mip",
                contrast_limits=(0.0, max(vmax, 1e-6)), multiscale=ml,
            )
            if mean_ratio.size:
                print(f"Clover/Ruby ratio (native Z): {ids.size} adhesions  "
                      f"mean={mean_ratio.mean():.3f} ± {mean_ratio.std():.3f}  "
                      f"[min={mean_ratio.min():.3f} max={mean_ratio.max():.3f}]  "
                      f"(per-adhesion mean over ruby>0 native voxels)", flush=True)
            else:
                print("Clover/Ruby ratio: no valid voxels (ruby == 0 everywhere "
                      "inside adhesions).", flush=True)

        def _err(exc):
            print(f"Clover/Ruby error: {exc}", flush=True)

        print("Computing clover/ruby ratio inside adhesions (native Z) ...",
              flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    ratio_btn.changed.connect(_ratio)

    # Save / load the spatial ratio layer so it can be reloaded without
    # re-running post-process + the ratio computation.
    ratio_path_w = FileEdit(value="clover_ruby_ratio.npy", label="ratio_path",
                            mode="w")
    save_ratio_w = PushButton(text="Save ratio layer")
    load_ratio_w = PushButton(text="Load ratio layer")
    RATIO_LAYER = "Clover/Ruby (adhesions)"

    def _save_ratio(*_):
        if RATIO_LAYER not in [lay.name for lay in viewer.layers]:
            print("No 'Clover/Ruby (adhesions)' layer to save — compute the "
                  "ratio first.", flush=True)
            return
        lyr = viewer.layers[RATIO_LAYER]
        d = lyr.data
        vol = np.asarray(d[0] if getattr(lyr, "multiscale", False) else d)
        vol = vol.astype(np.float32, copy=False)
        path = str(ratio_path_w.value)
        np.save(path, vol)
        print(f"Saved ratio layer -> {path} ({vol.nbytes / 1e6:.1f} MB)",
              flush=True)

    def _load_ratio(*_, after=None, visible=True):
        path = str(ratio_path_w.value)
        levels = state.get("levels", 1)
        # Ratio layer is on the native-Z grid, so reload at the native voxel size.
        scale = tuple(state.get("raw_voxel_zyx_um",
                                state.get("scale", (1.0, 1.0, 1.0))))

        @thread_worker
        def _run():
            return np.load(path)

        def _done(vol):
            ml = levels > 1
            nz = vol[vol > 0]
            vmax = float(np.percentile(nz, 99)) if nz.size else 1.0
            if RATIO_LAYER in [lay.name for lay in viewer.layers]:
                viewer.layers.remove(RATIO_LAYER)
            lyr = viewer.add_image(
                pyramid_views(vol, levels), name=RATIO_LAYER, scale=scale,
                colormap="viridis", blending="additive", rendering="mip",
                contrast_limits=(0.0, max(vmax, 1e-6)), multiscale=ml,
            )
            lyr.visible = visible
            print(f"Loaded ratio layer <- {path}.", flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Load ratio error: {exc}", flush=True)
            if after:
                after()

        print(f"Loading ratio layer from {path} ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    save_ratio_w.changed.connect(_save_ratio)
    load_ratio_w.changed.connect(_load_ratio)
    state["load_ratio"] = _load_ratio

    # ------------------------------------------------------------------ #
    # Section 4c — ImageJ-style ratio (whole volume, for comparison)       #
    # ------------------------------------------------------------------ #
    # Reproduces the old ImageJ RatioStacksAveIntensity workflow so it can be
    # compared with the adhesion-segmented ratio above. Whole-volume
    # clover/ruby, kept only where the ratio >= a threshold (ImageJ used ~2.0),
    # scaled x1000 to match ImageJ values, plus an average-intensity Z-projection
    # that excludes NaN (like ImageJ's "NaN Background" + "Z Project"). This uses
    # NO adhesion mask — the ROI is the ratio threshold, exactly as in ImageJ.

    ij_thresh_w = FloatSlider(value=2.0, min=0.0, max=10.0, step=0.1,
                              label="ij_ratio_threshold")
    imagej_btn = PushButton(text="ImageJ-style ratio + Z-avg")

    def _imagej_ratio(*_):
        clover_native = state.get("clover_raw_native")
        ruby_native = state.get("ruby_raw_native")
        if clover_native is None or ruby_native is None:
            print("ImageJ ratio: need native clover and ruby channels — reload "
                  "the image.", flush=True)
            return
        thr = float(ij_thresh_w.value)
        # Native-Z grid (matches ImageJ; no isotropic resampling).
        scale = tuple(state.get("raw_voxel_zyx_um",
                                state.get("scale", (1.0, 1.0, 1.0))))

        @thread_worker
        def _run():
            import warnings
            clover = clover_native.astype(np.float32)
            ruby = ruby_native.astype(np.float32)
            pos = ruby > 0
            ratio = np.zeros(clover.shape, dtype=np.float32)
            ratio[pos] = clover[pos] / ruby[pos]
            keep = pos & (ratio >= thr)
            # x1000 to match ImageJ; NaN outside the kept region so the Z-average
            # ignores it (numpy nanmean == ImageJ NaN-Background average).
            vol = np.where(keep, ratio * 1000.0, np.nan).astype(np.float32)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                zavg = np.nanmean(vol, axis=0).astype(np.float32)  # (Y, X)
            return vol, zavg

        def _done(result):
            vol, zavg = result
            sc = tuple(float(s) for s in scale)
            fv = vol[np.isfinite(vol)]
            lo = thr * 1000.0
            hi = float(np.percentile(fv, 99)) if fv.size else lo + 1.0
            clim = (lo, max(hi, lo + 1.0))
            for nm in ("Clover/Ruby ImageJ (3D)", "Clover/Ruby ImageJ (Z-avg)"):
                if nm in [l.name for l in viewer.layers]:
                    viewer.layers.remove(nm)
            viewer.add_image(vol, name="Clover/Ruby ImageJ (3D)", scale=sc,
                             colormap="inferno", blending="additive",
                             rendering="mip", contrast_limits=clim)
            viewer.add_image(zavg, name="Clover/Ruby ImageJ (Z-avg)",
                             scale=(sc[1], sc[2]), colormap="inferno",
                             contrast_limits=clim)
            n_keep = int(np.isfinite(vol).sum())
            print(f"ImageJ-style ratio: kept clover/ruby >= {thr:.2f} "
                  f"(x1000 >= {lo:.0f}), {n_keep:,} voxels; Z-avg {zavg.shape}. "
                  f"Values are (clover/ruby)*1000 to match ImageJ.", flush=True)

        def _err(exc):
            print(f"ImageJ ratio error: {exc}", flush=True)

        print("Computing ImageJ-style ratio + Z-average ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    imagej_btn.changed.connect(_imagej_ratio)

    for _w, _fn in ((proba_path_w, "adhesion_proba.npy"),
                    (ratio_path_w, "clover_ruby_ratio.npy"),
                    (inst_path_w, "adhesion_instances.npz")):
        state.setdefault("save_path_widgets", []).append((_w, _fn))

    # ------------------------------------------------------------------ #
    # Section 5 — Save / Load                                              #
    # ------------------------------------------------------------------ #

    model_path_w = FileEdit(value="adhesion_rf.joblib", label="model_path",
                             mode="w")
    save_w = PushButton(text="Save model")
    load_w = PushButton(text="Load model")

    def _save(*_):
        rf: AdhesionRF | None = state.get("adhesion_rf")
        if rf is None or not rf.is_trained:
            print("No trained RF to save.", flush=True)
            return
        rf.save(str(model_path_w.value))

    def _load(*_):
        try:
            rf = AdhesionRF.load(str(model_path_w.value))
            state["adhesion_rf"] = rf
        except Exception as exc:
            print(f"Load error: {exc}", flush=True)

    save_w.changed.connect(_save)
    load_w.changed.connect(_load)

    # ------------------------------------------------------------------ #
    # Assemble                                                             #
    # ------------------------------------------------------------------ #

    return Container(
        widgets=[
            hint_w, create_label_w,
            n_estimators_w, train_w,
            gate_by_clover_w, clover_dilation_um_w, predict_w,
            proba_path_w, save_proba_w, load_proba_w,
            prob_threshold_w,
            gate_by_dapi_w, d_dapi_um_w,
            gate_by_mesh_w, d_mesh_um_w,
            min_size_vox_w,
            split_merged_w, split_min_size_w,
            apply_w,
            inst_path_w, save_inst_w, load_inst_w,
            ratio_btn,
            ratio_path_w, save_ratio_w, load_ratio_w,
            ij_thresh_w, imagej_btn,
            model_path_w, save_w, load_w,
        ],
        labels=True,
    )
