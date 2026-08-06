"""VE-Cadherin boundary RF: pixel-classify endothelial junctions.

Same idea as the adhesion RF, but the classifier is trained on the VE-Cadherin
volume to label junction pixels instead of adhesions. A random-forest pixel
classifier learns the multiscale texture of the junctional network, which is far
more robust to the intracellular VE-cad speckle than any intensity threshold —
so the predicted probability map is a clean, continuous junction network that
the endothelial-cell counter can actually segment.

Workflow (mirrors the Adhesion RF tab):
  1. Create label layer  -> paint junction (1) and background (2) on VE-cad.
  2. Train RF            -> learns junction vs background on the VE-cad volume.
  3. Run RF predict      -> clean junction probability map
                            (state["cadherin_junction_proba"]).
  4. Save the model (.joblib) to reuse across images / headless.

The probability map is what the cell-count step (or the headless
cadherin_rf_batch script) consumes instead of a thresholded VE-cad channel.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from magicgui.widgets import Container, FileEdit, PushButton, SpinBox
from napari.qt.threading import thread_worker

from ..rf import AdhesionRF

LABELS_LAYER = "VE-cad boundary labels"
PROBA_LAYER = "VE-cad junction proba"


def build(state, viewer) -> Container:
    create_label_w = PushButton(text="Create label layer")
    n_estimators_w = SpinBox(value=200, min=10, max=1000, step=10,
                             label="n_estimators")
    train_w = PushButton(text="Train RF")
    predict_w = PushButton(text="Run RF predict")
    model_path_w = FileEdit(value="cadherin_rf.joblib", label="model_path",
                            mode="w")
    save_model_w = PushButton(text="Save model")
    load_model_w = PushButton(text="Load model")
    proba_path_w = FileEdit(value="cadherin_junction_proba.npy",
                            label="proba_path", mode="w")
    save_proba_w = PushButton(text="Save proba")
    load_proba_w = PushButton(text="Load proba")

    # ---- Create the annotation layer over the VE-cad volume ----
    def _create_labels(*_):
        cadh = state.get("cadh_raw")
        if cadh is None:
            print("Load a stack with a VE-Cadherin channel first.", flush=True)
            return
        scale = state.get("scale", (1.0, 1.0, 1.0))
        if LABELS_LAYER in [l.name for l in viewer.layers]:
            viewer.layers.remove(LABELS_LAYER)
        ann = np.zeros(np.asarray(cadh).shape, dtype=np.uint8)
        viewer.add_labels(ann, name=LABELS_LAYER, scale=scale)
        print(f"'{LABELS_LAYER}' created. Select it, then paint junctions with "
              "label 1 and background with label 2.", flush=True)

    # ---- Train ----
    def _train(*_, after=None):
        cadh = state.get("cadh_raw")
        if cadh is None:
            print("Load a stack first.", flush=True)
            if after:
                after()
            return
        if LABELS_LAYER not in [l.name for l in viewer.layers]:
            print(f"No '{LABELS_LAYER}' layer — create it and paint first.",
                  flush=True)
            if after:
                after()
            return
        ann_layer = viewer.layers[LABELS_LAYER]
        ann_vol = (np.asarray(ann_layer.data[0]) if getattr(
            ann_layer, "multiscale", False) else np.asarray(ann_layer.data)
        ).astype(np.uint8)
        if not np.any(ann_vol > 0):
            print("No annotations found. Paint junction (1) and background (2) "
                  "first.", flush=True)
            if after:
                after()
            return
        n_est = int(n_estimators_w.value)
        voxel_um = state.get("scale", (1.0, 1.0, 1.0))
        vol_f32 = np.asarray(cadh).astype(np.float32)

        @thread_worker
        def _run():
            rf = AdhesionRF(n_estimators=n_est)
            rf.train(vol_f32, ann_vol, voxel_um)
            return rf

        def _done(rf):
            state["cadherin_rf"] = rf
            print("VE-cad RF training complete.", flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Training error: {exc}", flush=True)
            if after:
                after()

        print(f"Training VE-cad RF (n_estimators={n_est}) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    # ---- Predict ----
    def _predict(*_, after=None):
        cadh = state.get("cadh_raw")
        rf = state.get("cadherin_rf")
        if cadh is None or rf is None:
            print("Need a loaded stack and a trained/loaded VE-cad RF.",
                  flush=True)
            if after:
                after()
            return
        voxel_um = state.get("scale", (1.0, 1.0, 1.0))
        vol_f32 = np.asarray(cadh).astype(np.float32)

        @thread_worker
        def _run():
            return rf.predict_proba(vol_f32, voxel_um, gate_indices=None)

        def _done(proba):
            state["cadherin_junction_proba"] = proba
            scale = state.get("scale", (1.0, 1.0, 1.0))
            if PROBA_LAYER in [l.name for l in viewer.layers]:
                viewer.layers.remove(PROBA_LAYER)
            viewer.add_image(proba, name=PROBA_LAYER, scale=scale,
                             colormap="magma", blending="additive")
            print(f"VE-cad junction proba computed (max={float(proba.max()):.2f}).",
                  flush=True)
            if after:
                after()

        def _err(exc):
            print(f"Predict error: {exc}", flush=True)
            if after:
                after()

        print("Running VE-cad RF predict (full volume) ...", flush=True)
        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    # ---- Save / load model + proba ----
    def _save_model(*_):
        rf = state.get("cadherin_rf")
        if rf is None:
            print("No trained VE-cad RF to save.", flush=True)
            return
        rf.save(str(model_path_w.value))

    def _load_model(*_):
        try:
            state["cadherin_rf"] = AdhesionRF.load(str(model_path_w.value))
        except Exception as exc:
            print(f"Load model error: {exc}", flush=True)

    def _save_proba(*_):
        proba = state.get("cadherin_junction_proba")
        if proba is None:
            print("No junction proba to save. Run predict first.", flush=True)
            return
        np.save(str(proba_path_w.value), proba)
        print(f"Junction proba saved to {proba_path_w.value}", flush=True)

    def _load_proba(*_):
        try:
            proba = np.load(str(proba_path_w.value))
        except Exception as exc:
            print(f"Load proba error: {exc}", flush=True)
            return
        state["cadherin_junction_proba"] = proba
        scale = state.get("scale", (1.0, 1.0, 1.0))
        if PROBA_LAYER in [l.name for l in viewer.layers]:
            viewer.layers.remove(PROBA_LAYER)
        viewer.add_image(proba, name=PROBA_LAYER, scale=scale, colormap="magma",
                         blending="additive")
        print(f"Junction proba loaded from {proba_path_w.value}", flush=True)

    create_label_w.changed.connect(_create_labels)
    train_w.changed.connect(_train)
    predict_w.changed.connect(_predict)
    save_model_w.changed.connect(_save_model)
    load_model_w.changed.connect(_load_model)
    save_proba_w.changed.connect(_save_proba)
    load_proba_w.changed.connect(_load_proba)

    # Headless / Run-all hooks.
    state["run_cadherin_rf_train"] = lambda after=None: _train(after=after)
    state["run_cadherin_rf_predict"] = lambda after=None: _predict(
        after=(lambda: (_save_proba(), after() if after else None)))
    state["load_cadherin_proba"] = _load_proba

    return Container(widgets=[
        create_label_w, n_estimators_w, train_w, predict_w,
        model_path_w, save_model_w, load_model_w,
        proba_path_w, save_proba_w, load_proba_w,
    ])
