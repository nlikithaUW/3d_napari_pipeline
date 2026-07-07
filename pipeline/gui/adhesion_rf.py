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

    def _predict(*_):
        rf: AdhesionRF | None = state.get("adhesion_rf")
        if rf is None or not rf.is_trained:
            print("No trained RF. Train or load a model first.", flush=True)
            return
        clover_raw = state.get("clover_raw")
        if clover_raw is None:
            print("Load a stack first.", flush=True)
            return

        voxel_um = state.get("scale", (1.0, 1.0, 1.0))
        vol_f32 = clover_raw.astype(np.float32)
        levels = state.get("levels", 1)
        scale = voxel_um
        gate_clover = bool(gate_by_clover_w.value)
        dilation_um = float(clover_dilation_um_w.value)

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
                        print(f"Clover gate: {len(gate_indices):,} / {mask.size:,} "
                              f"voxels ({pct:.1f}%)", flush=True)
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

            def _err(exc):
                print(f"Predict error: {exc}", flush=True)

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

            w = _run_otsu()
            w.returned.connect(_otsu_done)
            w.errored.connect(_otsu_err)
            w.start()
        else:
            _start_rf()

    predict_w.changed.connect(_predict)

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
    split_merged_w = CheckBox(value=False, label="split_merged")
    split_min_size_w = Slider(value=1000, min=0, max=10000,
                               label="split_min_size_vox")
    split_min_size_w.visible = False

    apply_w = PushButton(text="Apply post-process")

    def _on_gate_dapi(val):  d_dapi_um_w.visible = bool(val)
    def _on_gate_mesh(val):  d_mesh_um_w.visible = bool(val)
    def _on_split(val):      split_min_size_w.visible = bool(val)

    gate_by_dapi_w.changed.connect(_on_gate_dapi)
    gate_by_mesh_w.changed.connect(_on_gate_mesh)
    split_merged_w.changed.connect(_on_split)

    def _apply(*_):
        proba_vol = state.get("adhesion_proba")
        if proba_vol is None:
            print("No probability map. Run RF predict first.", flush=True)
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
            lyr_data = pyramid_views(labels, levels)
            ml = levels > 1
            if name in [lay.name for lay in viewer.layers]:
                viewer.layers.remove(name)
            viewer.add_labels(lyr_data, name=name, scale=scale,
                              multiscale=ml)

        def _err(exc):
            print(f"Apply error: {exc}", flush=True)

        print("Applying post-process ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    apply_w.changed.connect(_apply)

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
            prob_threshold_w,
            gate_by_dapi_w, d_dapi_um_w,
            gate_by_mesh_w, d_mesh_um_w,
            min_size_vox_w,
            split_merged_w, split_min_size_w,
            apply_w,
            model_path_w, save_w, load_w,
        ],
        labels=True,
    )
