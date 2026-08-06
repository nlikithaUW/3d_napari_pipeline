"""Headless batch runner: run the full pipeline on every image in a folder.

Drives the *same* pipeline code the GUI "Run all" uses, but with the napari
window hidden — so a whole folder of images is processed unattended, with no
rendering to slow it down. Each image gets its own ``<stem>_saved`` folder with
all intermediates, the measurements CSV, and the settings JSON.

Usage (in the image_analysis conda env, from the 3d_napari_pipeline folder):

    python -m pipeline.batch --dir "C:\\Users\\Nliki\\Desktop\\All Images" \
        --model "adhesion_rf.joblib"

Options:
    --dir      folder of OME-TIFFs (same channel layout as the GUI expects)
    --model    path to the trained RF model (.joblib)
    --params   params.json (default: 3d_images_folder/params.json)
    --lut      force LUT (default: the bundled Clover/mRuby2 LUT)
    --pattern  glob for images inside --dir (default: *.tif)
    --levels   pyramid levels (default 1 = fastest; no multiscale needed headless)

Run it in its own terminal; it prints progress per image/step. Smoke-test on one
image first and confirm the CSV matches a GUI "Run all" export.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Run Qt/napari truly headless BEFORE importing them: an offscreen platform has
# no real GL window, so layer draws can't crash the process (the show=False path
# still creates a hidden GL window that renders). Set only if the user hasn't
# already chosen a platform. Also force a non-interactive matplotlib backend.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")

import napari
from qtpy.QtCore import QTimer
from qtpy.QtWidgets import QApplication

from .gui import (
    adhesion_rf, cadherin, cell_mesh, clover, display, loader, measure,
    nuclei, roi, vessel_mesh,
)
from .rf import AdhesionRF

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PARAMS = ROOT / "3d_images_folder" / "params.json"
DEFAULT_LUT = ROOT / "Fret_LUT" / "Clover_mRuby2_GGSGGS7_force_new.txt"

# Same order the GUI "Run all" uses, EXCEPT the final Measure + export step,
# which is intentionally left out: it crashes headless (a GL/Qt draw in the
# measure tab). The batch produces every intermediate (proba, instances, mesh,
# curvature, DAPI, cadherin, fibroblast count) in each <stem>_saved folder;
# run Measure + export yourself per image in the GUI (Load all saved → Measure
# → export), which works fine there.
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
]


def _fresh_state() -> dict:
    return {
        "dapi_raw": None, "cadh_raw": None, "clover_raw": None,
        "force_vol": None, "adhesion_rf": None, "adhesion_proba": None,
        "scale": (1.0, 1.0, 1.0),
        "dapi_mask_layer": None, "cadh_mask_layer": None,
        "clover_mask_layer": None,
        "dapi_state": {"labels": None, "sizes": None},
        "cadh_state": {"labels": None, "sizes": None, "centroids": None,
                       "normals": None, "label_ids": None, "keep_mask": None},
        "clover_state": {"labels": None, "sizes": None, "method": None},
        "nuc_cache": {"labels_id": None, "centroids": None, "normals": None,
                      "lids": None},
        "mesh_cache": {"offset_centroids": None, "outward": None,
                       "inward_offset_um": 0.0, "kdtree": None},
        "vessel_mesh_cache": {"verts": None, "faces": None},
        "raw_layers": [], "hists": [],
        "headless": True,   # measure skips display-only plotting in batch mode
    }


def _free_gpu():
    try:
        import cupy as cp
        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass


def main() -> None:
    ap = argparse.ArgumentParser(description="Headless batch pipeline runner.")
    ap.add_argument("--dir", required=True, type=Path,
                    help="folder of OME-TIFFs to process")
    ap.add_argument("--model", required=True, type=Path,
                    help="trained RF model (.joblib)")
    ap.add_argument("--params", type=Path, default=DEFAULT_PARAMS)
    ap.add_argument("--lut", type=Path, default=DEFAULT_LUT)
    ap.add_argument("--pattern", default="*.tif")
    ap.add_argument("--levels", type=int, default=1)
    ap.add_argument("--skip-existing", action="store_true",
                    help="skip images that already have a measurements CSV in "
                         "their <stem>_saved folder")
    ap.add_argument("--timeout", type=float, default=0.0,
                    help="per-image timeout in seconds; if an image takes longer "
                         "(hang), it is skipped and the batch continues. "
                         "0 = no timeout (default).")
    args = ap.parse_args()

    images = sorted(Path(args.dir).glob(args.pattern))
    if args.skip_existing:
        kept = []
        for tif in images:
            saved = tif.parent / f"{tif.stem}_saved"
            # batch's final artifact is the instance labels (measure runs later
            # in the GUI); skip images that already have them.
            if (saved / "adhesion_instances.npz").exists():
                print(f"  skip (already processed): {tif.name}", flush=True)
            else:
                kept.append(tif)
        images = kept
    if not images:
        print(f"No images to process in {args.dir} "
              f"(pattern={args.pattern}, skip_existing={args.skip_existing})",
              flush=True)
        sys.exit(1)
    print(f"Batch: {len(images)} image(s) in {args.dir}\n"
          f"  model={args.model}\n  params={args.params}\n  lut={args.lut}\n",
          flush=True)

    # Hidden viewer keeps a GL context (so surface/labels layers work) without
    # popping a window; no logpanel install so prints go to this terminal.
    viewer = napari.Viewer(ndisplay=3, show=False)
    state = _fresh_state()

    # Build every tab once — this registers state["run_*"], state["load_*"] and
    # the per-image save-path widgets. Order matches the GUI.
    loader.build(state, viewer, args.dir, args.params, args.lut)
    display.build(state, viewer)
    nuclei.build(state, viewer)
    cadherin.build(state, viewer)
    cell_mesh.build(state, viewer)
    vessel_mesh.build(state, viewer)
    clover.build(state, viewer)
    adhesion_rf.build(state, viewer)
    roi.build(state, viewer)
    # NOTE: measure.build is intentionally NOT called — Measure isn't run in the
    # batch (it crashes headless), and skipping its build avoids creating the
    # matplotlib Qt canvas at all. Run Measure per image in the GUI afterward.

    try:
        state["adhesion_rf"] = AdhesionRF.load(str(args.model))
    except Exception as exc:
        print(f"Could not load RF model {args.model}: {exc}", flush=True)
        sys.exit(1)

    def _finish():
        print("\nBatch complete. Intermediates saved for every image "
              "(proba, instances, mesh, curvatures, DAPI, cadherin, fibroblast). "
              "Run Measure + export per image in the GUI (Load all saved → "
              "Measure → export) to produce the CSVs.", flush=True)
        _free_gpu()
        try:
            viewer.close()
        except Exception:
            pass
        app = QApplication.instance()
        if app is not None:
            app.quit()

    # `cur` tracks the image currently being processed; every "advance to next
    # image" goes through _advance(i), which is idempotent — so an image both
    # finishing AND its watchdog firing only advances once. This guarantees one
    # bad/hung image can never stop the batch or double-skip.
    cur = {"i": -1}
    wd = {"timer": None}

    def _disarm():
        if wd["timer"] is not None:
            try:
                wd["timer"].stop()
            except Exception:
                pass
            wd["timer"] = None

    def _free_image_memory():
        # Drop the finished image's big arrays + all layers BEFORE the next image
        # loads, so peak RAM is ~one image's worth instead of two. Keeps the RF
        # model, widgets, run_*/load_* functions, scale, headless flag, etc.
        import gc
        try:
            for layer in list(viewer.layers):
                viewer.layers.remove(layer)
        except Exception:
            pass
        for k in ("dapi_raw", "cadh_raw", "clover_raw", "ruby_raw",
                  "clover_raw_native", "ruby_raw_native", "force_vol",
                  "adhesion_proba", "clover_ruby_ratio", "vessel_centerline",
                  "measure_result", "dapi_mask_layer", "cadh_mask_layer",
                  "clover_mask_layer"):
            state[k] = None
        state["dapi_state"] = {"labels": None, "sizes": None}
        state["cadh_state"] = {"labels": None, "sizes": None, "centroids": None,
                               "normals": None, "label_ids": None,
                               "keep_mask": None}
        state["clover_state"] = {"labels": None, "sizes": None, "method": None}
        state["vessel_mesh_cache"] = {"verts": None, "faces": None}
        state["mesh_cache"] = {"offset_centroids": None, "outward": None,
                               "inward_offset_um": 0.0, "kdtree": None}
        state["raw_layers"] = []
        gc.collect()
        _free_gpu()

    def _advance(i: int, reason: str = ""):
        if cur["i"] != i:        # already moved past image i — ignore
            return
        cur["i"] = -2
        _disarm()
        if reason:
            print(f"  {images[i].name}: {reason}", flush=True)
        _free_image_memory()
        process(i + 1)

    def process(i: int):
        if i >= len(images):
            _finish()
            return
        cur["i"] = i
        tif = images[i]
        print(f"\n===== [{i + 1}/{len(images)}] {tif.name} =====", flush=True)
        state["run_settings"] = {}          # clean settings log per image
        state["roi_polygons_yx_um"] = None  # never gate the batch by an ROI
        state["roi_yx_um"] = None
        if args.timeout > 0:                # per-image hang guard
            _disarm()
            t = QTimer()
            t.setSingleShot(True)
            t.timeout.connect(
                lambda ci=i: _advance(ci, f"TIMED OUT after {args.timeout:.0f}s "
                                          f"— skipped"))
            t.start(int(args.timeout * 1000))
            wd["timer"] = t
        try:
            loader._load_image(state, viewer, tif, args.params, args.lut,
                               int(args.levels),
                               after=lambda: _after_load(i))
        except Exception as exc:
            _advance(i, f"load raised: {exc}; skipped")

    def _after_load(i: int):
        try:
            if state.get("force_vol") is None:
                _advance(i, "load failed; skipped")
                return
            steps = [(label, state[key]) for label, key in _RUN_STEPS
                     if state.get(key)]

            def _next(j: int):
                if j >= len(steps):
                    _advance(i, "done.")
                    return
                label, fn = steps[j]
                print(f"  step {j + 1}/{len(steps)} — {label} ...", flush=True)
                try:
                    fn(after=lambda: _next(j + 1))
                except Exception as exc:
                    print(f"  step '{label}' failed: {exc}; continuing.",
                          flush=True)
                    _next(j + 1)

            _next(0)
        except Exception as exc:
            _advance(i, f"unexpected error ({exc}); skipped")

    # Kick off once the Qt event loop is running, then run the loop.
    QTimer.singleShot(0, lambda: process(0))
    napari.run()


if __name__ == "__main__":
    main()
