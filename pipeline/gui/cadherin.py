"""Unified Cadherin tab.

Per run:
  local threshold (Niblack / Sauvola / Bernsen) ANDed with a global intensity
  floor (median + k·MAD, auto-seeded from threshold_triangle on load)
  → labeled components → per-component centroid + PCA normal.

Sliders then AND a size filter with a band filter (|signed distance to the
current DAPI cell mesh| ≤ band_um); the band filter is skipped if no DAPI mesh
exists yet. Centroids/normals feed the Cell mesh "Refine with cadherin" step;
the binary mask feeds the Vessel mesh tab.
"""

from __future__ import annotations

import time

import numpy as np
from magicgui import magicgui
from magicgui.widgets import Container, FileEdit, PushButton
from napari.qt.threading import thread_worker

from ..cell_mesh import label_centroids_normals, signed_distance_to_mesh_vertices
from .. import saved_io as sio
from ..denoise import (
    auto_min_size_otsu,
    bernsen_mask,
    label_and_sizes,
    niblack_mask,
    robust_floor_value,
    sauvola_mask,
    suggest_floor_k,
)
from .widgets import (Histogram, call_all, debounce, pyramid_views,
                      skip_if_missing, with_after)

METHOD_NIBLACK = "niblack"
METHOD_SAUVOLA = "sauvola"
METHOD_BERNSEN = "bernsen"
METHODS = [METHOD_NIBLACK, METHOD_SAUVOLA, METHOD_BERNSEN]

# Floor for the auto-seeded despeckle size filter. `auto_min_size_otsu` returns
# ~5 on the WT stack, which leaves ~99% of (speckle) components in the mask and
# wrecks the downstream vessel alpha shape; ~12 drops the speckle. Fully tunable
# on the slider afterwards.
# 25, not 12, since 2026-09-02 — the operator's working value, and it is a
# **measurement**, not a preference: at the new mask (floor_k 0, window_xy 101,
# niblack_k 2) `V3_..._09step_001` keeps 17,297 components at 12 and 3,794 at
# 25, and the extra speckle changes the auto alpha from 9.77 to 6.82 um, which
# changes the solid, which loses the vessel's junction entirely — the graph
# reads 1 edge / 0 junctions instead of 3 edges / 1 junction. The Y is ground
# truth from the experimenter. min_size, alpha and topology are coupled now
# that alpha is auto, so this floor is load-bearing rather than cosmetic.
MIN_SIZE_FLOOR = 25


# --- what the cadherin labels npz holds --------------------------------------- #
# Two keys: `labels` and `sizes`. That is the whole live half — `filter_by_size`
# in `batch_mesh` is what consumes them, and `min_size` moves the auto alpha and
# with it the topology (P51), so the labeling itself is load-bearing.
#
# The per-component PCA `centroids` / `normals` / `label_ids` are **not** saved
# any more (P55). Their only consumers are the band filter and the Cell mesh
# tab's "Refine with cadherin", both of which need a per-cell DAPI mesh that no
# headless step builds — there is no cell-mesh entry in `batch.py::_RUN_STEPS`.
# They stay lazy in `state` for the interactive button, and are recomputed by
# `_ensure_centroids` on demand, which is exactly what happens after a fresh
# threshold anyway. Nothing was reading them back from disk.


def _lazy_in(d, key):
    """A field that is no longer saved: absent, size-0, or (in files written
    before P55) a pickled `None` all mean the same thing — not computed."""
    if key not in d.files:
        return None
    a = d[key]
    if a.dtype == object or a.ndim == 0 or a.size == 0:
        return None
    return a


def build(state, viewer) -> Container:
    hist = Histogram("Cadherin component sizes",
                     "log10(size in voxels)", "deeppink")
    state.setdefault("hists", []).append(hist)

    no_mesh_notified = {"done": False}

    def _redraw(min_size: int):
        sizes = state["cadh_state"]["sizes"]
        hist.update(sizes, threshold=min_size if min_size > 0 else None,
                    threshold_label=f"min_size={min_size}")

    def _ensure_centroids():
        """Per-component centroids/normals, computed on first use and cached.

        This is the expensive half of the Cadherin tab and it is **not** on the
        main path. `label_centroids_normals` is a Python loop over every
        component doing `argwhere` + covariance + `eigh` each, and the tab calls
        it with `min_voxels=1`, so it runs on every 1-voxel speck. Measured on
        the reference image, against a constant ~3.7 s for `label_and_sizes`:

            floor_k 5.0 -> 170,236 components ->  4.98 s
            floor_k 4.0 -> 375,467 components ->  9.92 s
            floor_k 3.0 -> 1,335,402 components -> 34.45 s

        Component count explodes as the floor drops (7.8x for 3.4x the voxels),
        so this dominated every run and made a permissive `floor_k` unusable.

        Nothing on the main path reads the result: the vessel mesh consumes the
        binary *mask*, not the labels. The only consumer in the codebase is
        `gui/cell_mesh.py`'s optional "Refine with cadherin" button, which is
        not in `batch.py`'s `_RUN_STEPS` — so on every headless run this was
        computed and thrown away.

        Kept memoised rather than deleted because the band filter needs it too.
        """
        cs = state["cadh_state"]
        if cs.get("centroids") is not None:
            return cs["centroids"], cs["normals"], cs["label_ids"]
        labels = cs.get("labels")
        if labels is None:
            return None, None, None
        n = cs["sizes"].size
        print(f"Cadherin: computing centroids/normals for {n:,} components "
              f"(first use) ...", flush=True)
        t0 = time.perf_counter()
        c, nm, lids = label_centroids_normals(labels, state["scale"],
                                              min_voxels=1)
        cs["centroids"], cs["normals"], cs["label_ids"] = c, nm, lids
        print(f"Cadherin: centroids DONE in {time.perf_counter() - t0:.1f} s "
              f"({len(c):,} of {n:,} components).", flush=True)
        return c, nm, lids

    def _cadh_vertices():
        """(centroids, normals, keep_mask), materialising all three on demand.

        The single entry point for consumers outside this tab, so they never
        have to know the fields are lazy or how `keep_mask` is derived.
        """
        c, nm, _ = _ensure_centroids()
        if c is None:
            return None, None, None
        cs = state["cadh_state"]
        if cs.get("keep_mask") is None:
            _compute_mask(int(filter_controls.min_size.value),
                          float(filter_controls.band_um.value))
        return cs["centroids"], cs["normals"], cs["keep_mask"]

    state["cadh_vertices"] = _cadh_vertices

    def _compute_mask(min_size: int, band_um: float):
        cs = state["cadh_state"]
        labels = cs["labels"]
        sizes = cs["sizes"]
        layer = state.get("cadh_mask_layer")
        if labels is None or layer is None:
            return None
        n_labels = sizes.size
        size_pass = sizes >= max(int(min_size), 0)
        keep_label = np.zeros(n_labels + 1, dtype=bool)
        keep_label[1:] = size_pass

        mesh = state.get("mesh_cache") or {}
        kdtree = mesh.get("kdtree")
        offset_centroids = mesh.get("offset_centroids")
        outward = mesh.get("outward")
        # Only materialise centroids if the band filter will actually use them.
        # With band_um == 0 (the default) this stays lazy.
        want_band = (band_um > 0 and kdtree is not None
                     and offset_centroids is not None and outward is not None
                     and len(offset_centroids) > 0)
        centroids, lids = cs.get("centroids"), cs.get("label_ids")
        if want_band:
            centroids, _n, lids = _ensure_centroids()
        use_band = want_band and centroids is not None and len(centroids) > 0

        if use_band:
            sd = signed_distance_to_mesh_vertices(
                centroids, offset_centroids, outward, tree=kdtree
            )
            band_pass = np.abs(sd) <= float(band_um)
            surviving = keep_label[lids] & band_pass
            keep_label = np.zeros(n_labels + 1, dtype=bool)
            keep_label[lids] = surviving
            keep_mask = surviving
        else:
            if band_um > 0 and not no_mesh_notified["done"]:
                print("Cadherin: no DAPI mesh yet — band filter skipped "
                      "(build the cell mesh first).", flush=True)
                no_mesh_notified["done"] = True
            # `lids` is None while centroids are unmaterialised. Leave
            # `keep_mask` as None rather than forcing the expensive computation
            # here — `_cadh_vertices` fills it in when a consumer actually asks.
            if lids is None:
                keep_mask = None
            else:
                keep_mask = (size_pass[lids - 1] if lids.size
                             else np.zeros(0, dtype=bool))

        cs["keep_mask"] = keep_mask
        return keep_label[labels].astype(np.uint8)

    @magicgui(
        auto_call=True,
        min_size={"widget_type": "Slider", "min": 0, "max": 5000},
        band_um={"widget_type": "FloatSlider", "min": 0.0, "max": 50.0, "step": 0.5},
    )
    def filter_controls(min_size: int = MIN_SIZE_FLOOR, band_um: float = 0.0):
        _apply(min_size, band_um)

    @debounce(150)
    def _apply(min_size: int, band_um: float):
        layer = state.get("cadh_mask_layer")
        if layer is None or state["cadh_state"]["labels"] is None:
            return
        mask = _compute_mask(min_size, band_um)
        if mask is None:
            return
        layer.data = pyramid_views(mask, state.get("levels", 1))
        layer.visible = True
        _redraw(min_size)

    @magicgui(
        call_button="Apply threshold + label",
        method={"widget_type": "ComboBox", "choices": METHODS},
        # Max 201, not 101. The shipped default is 101 and a slider ceiling is
        # a widget limit, not a physical one — CLAUDE.md's rule is that a
        # hand-tuned value sitting exactly at one may have been truncated
        # rather than converged, so the ceiling must stay clear of the default.
        window_xy={"widget_type": "Slider", "min": 5, "max": 201, "step": 1},
        niblack_k={"widget_type": "FloatSlider", "min": -1.0, "max": 4.0, "step": 0.05},
        sauvola_k={"widget_type": "FloatSlider", "min": 0.0, "max": 1.0, "step": 0.01},
        bernsen_contrast={"widget_type": "FloatSlider", "min": 0.0, "max": 1.0, "step": 0.01},
        # 0 disables the global floor entirely and is now the default; the
        # slider therefore has to reach it, which it does.
        floor_k={"widget_type": "FloatSlider", "min": 0.0, "max": 10.0, "step": 0.1},
    )
    def threshold_controls(
        method: str = METHOD_NIBLACK,
        # 101, not 53 — the operator's working value, adopted as the default on
        # 2026-09-02. A wider window makes the local mean a better estimate of
        # background under a large, dim, sparse structure like the far wall,
        # where a 53-voxel window is partly inside the wall itself.
        window_xy: int = 101,
        niblack_k: float = 2.0,
        sauvola_k: float = 0.2,
        bernsen_contrast: float = 0.15,
        # 0.0, not 1.0 — **no global intensity floor**, adopted 2026-09-02.
        # The floor is what [P28] measures deleting 59-85% of the dim wall: it
        # is a single global cut, so it removes the far wall wherever the stain
        # is weak regardless of how well that wall stands out locally. With
        # Niblack live the local term is doing the work and the floor is a
        # backstop that costs more than it earns. **This changes every mask.**
        floor_k: float = 0.0,
    ):
        # This body is the "Apply threshold + label" button action. Visibility
        # is driven separately off the `method` widget's `changed` signal, so
        # switching method re-lays-out the panel WITHOUT launching a build.
        _run_cadherin()

    # Which knob each method actually reads. `floor_k` is deliberately absent:
    # it is ANDed in for EVERY method (see `_run`), so it is always live and
    # always shown.
    _METHOD_KNOBS = {
        METHOD_NIBLACK: ("window_xy", "niblack_k"),
        METHOD_SAUVOLA: ("window_xy", "sauvola_k"),
        METHOD_BERNSEN: ("window_xy", "bernsen_contrast"),
    }

    def _sync_method_visibility(*_):
        live = _METHOD_KNOBS.get(threshold_controls["method"].value, ())
        for name in ("window_xy", "niblack_k", "sauvola_k", "bernsen_contrast"):
            threshold_controls[name].visible = name in live

    threshold_controls["method"].changed.connect(_sync_method_visibility)
    _sync_method_visibility()

    def _run_cadherin(after=None):
        cadh_raw = state.get("cadh_raw")
        if cadh_raw is None:
            print("Load a stack with a 'substrate' channel first.", flush=True)
            call_all(after)
            return
        scale = state["scale"]
        method = threshold_controls["method"].value
        window_xy = int(threshold_controls["window_xy"].value)
        niblack_k = float(threshold_controls["niblack_k"].value)
        sauvola_k = float(threshold_controls["sauvola_k"].value)
        bernsen_contrast = float(threshold_controls["bernsen_contrast"].value)
        floor_k = float(threshold_controls["floor_k"].value)
        state.setdefault("run_settings", {})["cadherin"] = {
            "method": method, "window_xy": window_xy, "niblack_k": niblack_k,
            "sauvola_k": sauvola_k, "bernsen_contrast": bernsen_contrast,
            "floor_k": floor_k}

        @thread_worker
        def _run():
            vol = cadh_raw
            if method == METHOD_NIBLACK:
                mask = niblack_mask(vol, window_xy, niblack_k, scale)
            elif method == METHOD_SAUVOLA:
                mask = sauvola_mask(vol, window_xy, sauvola_k, scale)
            else:
                mask = bernsen_mask(vol, window_xy, bernsen_contrast, scale)
            floor = robust_floor_value(vol, floor_k)
            if floor is not None:
                mask = (mask.astype(bool) & (vol > floor)).astype(np.uint8)
            labels, sizes = label_and_sizes(mask)
            # Centroids/normals are NOT computed here — see `_ensure_centroids`.
            suggested = auto_min_size_otsu(sizes)
            return labels, sizes, suggested

        def _done(result):
            labels, sizes, suggested = result
            cs = state["cadh_state"]
            cs["labels"] = labels
            cs["sizes"] = sizes
            cs["centroids"] = None      # lazy, see `_ensure_centroids`
            cs["normals"] = None
            cs["label_ids"] = None
            cs["keep_mask"] = None
            no_mesh_notified["done"] = False
            slider_max = int(sizes.max()) if sizes.size else 5000
            filter_controls.min_size.max = slider_max
            filter_controls.min_size.value = max(suggested, MIN_SIZE_FLOOR)
            # Refresh the mask layer synchronously — the auto-called filter is
            # debounced 150 ms, but the "Run all" chain's next step (vessel mesh)
            # reads this layer immediately, so it must be current *now*.
            _layer = state.get("cadh_mask_layer")
            _m = _compute_mask(int(filter_controls.min_size.value),
                               float(filter_controls.band_um.value))
            if _layer is not None and _m is not None:
                _layer.data = pyramid_views(_m, state.get("levels", 1))
                _layer.visible = True
            _redraw(int(filter_controls.min_size.value))
            filter_controls()
            # Printed HERE, not when the worker returned. The threshold itself
            # is only part of the wait: `_compute_mask` + `pyramid_views` +
            # `_redraw` run after it on the GUI thread and take seconds on a
            # full stack, so a "done" printed on return arrives while the
            # viewer is still updating and cannot be trusted to mean the image
            # is final. Carries the voxel count and the elapsed time because
            # "is it still going?" is the actual question being asked.
            vox = int(_m.sum()) if _m is not None else 0
            print(
                f"Cadherin: {method} DONE in {time.perf_counter() - _t0:.1f} s "
                f"— {sizes.size:,} components, "
                f"{vox:,} voxels in the displayed mask. "
                f"Suggested min_size = {suggested} "
                f"(applied {int(filter_controls.min_size.value)}).",
                flush=True,
            )

        def _err(exc):
            print(f"Cadherin error: {exc}", flush=True)

        floor = robust_floor_value(cadh_raw, floor_k)
        floor_str = f"{floor:.1f}" if floor is not None else "off"
        print(f"Cadherin: {method} running (window_xy={window_xy}, "
              f"floor_k={floor_k} → floor={floor_str}) ...", flush=True)
        _t0 = time.perf_counter()
        w = _run()
        w.returned.connect(with_after(_done, after, "Cadherin"))
        w.errored.connect(with_after(_err, after, "Cadherin"))
        w.start()

    def _seed_floor():
        """Set `floor_k` from `suggest_floor_k` — **explicit, not automatic.**

        This used to run from the loader on every image load, so it silently
        overwrote whatever `floor_k` was set to: the tab's default never
        survived to the build, and a value typed for the previous image was
        replaced without a word. It is the P26 hazard in miniature — the
        shipped default was not what ran (measured 2026-08-24: the default said
        1.0, the build logged `floor_k=4.595517`).

        Auto-seeding made sense while the floor was the *primary* gate, which
        is what it was for as long as the local threshold was Sauvola and
        therefore inert (P36). With Niblack live the floor is a weak backstop
        and its value is a deliberate choice, so it is no longer imposed.
        """
        raw = state.get("cadh_raw")
        if raw is None:
            return
        threshold_controls.floor_k.value = min(
            suggest_floor_k(raw), threshold_controls.floor_k.max
        )
        print(f"Cadherin: floor_k seeded to "
              f"{threshold_controls.floor_k.value:.2f} from suggest_floor_k.",
              flush=True)

    seed_floor_w = PushButton(text="Suggest floor_k")
    seed_floor_w.changed.connect(lambda *_: _seed_floor())

    # ---- Save / load the cadherin components (mask + labels for re-tuning) ----
    cadh_path_w = FileEdit(value=sio.filename(sio.CADHERIN_LABELS),
                           label="cadh_path",
                           mode="w")
    save_cadh_w = PushButton(text="Save VE-Cadherin mask")
    load_cadh_w = PushButton(text="Load VE-Cadherin mask")

    def _save_cadh(*_):
        cs = state.get("cadh_state") or {}
        if cs.get("labels") is None:
            print("No cadherin components to save — run 'Apply threshold + "
                  "label' first.", flush=True)
            return
        path = str(cadh_path_w.value)
        # Two keys only — see the module header. The three PCA fields used to be
        # written here as `np.asarray(None)`, which is a 0-d *object* array:
        # numpy pickles it on save and then refuses to unpickle it on load, so
        # every file saved without them was unreadable (P55).
        sio.save_cadherin_labels(path, cs["labels"], cs["sizes"])
        print(f"Saved VE-Cadherin components -> {path} "
              f"({np.asarray(cs['sizes']).size} components).", flush=True)

    def _load_cadh(*_, after=None, visible=True):
        path = str(cadh_path_w.value)
        if skip_if_missing(path, after, "Load cadherin",
                           "VE-Cadherin components file"):
            return

        @thread_worker
        def _run():
            # `load_cadherin_labels` is the one loader that sets
            # allow_pickle, for the pre-P55 files that stored the three PCA
            # fields as pickled `None`; it hands them back as None.
            d = sio.load_cadherin_labels(path)
            return (d.labels.astype(np.int32), d.sizes,
                    d.centroids, d.normals, d.label_ids)

        def _done(result):
            labels, sizes, centroids, normals, lids = result
            cadh_raw = state.get("cadh_raw")
            if cadh_raw is not None and labels.shape != cadh_raw.shape:
                print(f"Warning: loaded cadherin {labels.shape} != current "
                      f"{cadh_raw.shape}; load the matching image first.",
                      flush=True)
            cs = state["cadh_state"]
            cs["labels"] = labels
            cs["sizes"] = sizes
            cs["centroids"] = centroids
            cs["normals"] = normals
            cs["label_ids"] = lids
            cs["keep_mask"] = None
            no_mesh_notified["done"] = False
            suggested = auto_min_size_otsu(sizes)
            filter_controls.min_size.max = int(sizes.max()) if sizes.size else 5000
            filter_controls.min_size.value = max(suggested, MIN_SIZE_FLOOR)
            filter_controls()  # rebuild the "VE-Cadherin mask" layer
            if not visible and state.get("cadh_mask_layer") is not None:
                state["cadh_mask_layer"].visible = False
            print(f"Loaded VE-Cadherin components <- {path}: {sizes.size} "
                  f"components.", flush=True)

        def _err(exc):
            print(f"Load cadherin error: {exc}", flush=True)

        print(f"Loading VE-Cadherin components from {path} ...", flush=True)
        w = _run()
        w.returned.connect(with_after(_done, after, "Load cadherin"))
        w.errored.connect(with_after(_err, after, "Load cadherin"))
        w.start()

    save_cadh_w.changed.connect(_save_cadh)
    load_cadh_w.changed.connect(_load_cadh)
    state["load_cadh"] = _load_cadh

    # Chainable "Run all" step: force Sauvola, compute, save, then `after`.
    def _run_cadherin_step(after=None):
        # **Does NOT override `method` any more.** This forced METHOD_SAUVOLA,
        # so every saved artifact was built with Sauvola while the tab's
        # interactive default was Niblack — you tuned one method and shipped
        # another. P36 then measured that Sauvola is *inert* here (its global
        # normaliser makes the threshold fall below background, so `sauvola AND
        # floor` is bit-identical to `floor` alone), which made the drift a
        # correctness problem rather than an inconsistency: the shipped mask had
        # no local adaptivity at all. The batch path now runs whatever the tab
        # is set to, which is the `build_mesh` doctrine applied here.
        _run_cadherin(after=lambda: call_all(_save_cadh, after,
                                             label="run_cadherin"))
    state["run_cadherin"] = _run_cadherin_step
    state.setdefault("save_path_widgets", []).append(
        (cadh_path_w, sio.filename(sio.CADHERIN_LABELS)))

    return (Container(widgets=[threshold_controls, seed_floor_w,
                               filter_controls,
                              cadh_path_w, save_cadh_w, load_cadh_w],
                     labels=False),
            hist.canvas)
