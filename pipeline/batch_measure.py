"""Headless Measure + export (no napari) for the batch.

For each image whose `<stem>_saved` folder has adhesion instances but no settings
JSON, this loads the raw stack (to recompute the force volume + native
clover/ruby channels, exactly as the loader does) and the saved intermediates
(instances, vessel mesh, DAPI labels, surface curvature, centerline, fibroblast
count), runs the identical measurement (`measure_core.compute_measurements`), and
writes the timestamped CSV + settings JSON — the same outputs the GUI Measure tab
would. Pure Python: no napari, no Qt, no GL, so it can't hit the headless crash.

Usage (image_analysis env, from the 3d_napari_pipeline folder):

    python -m pipeline.batch_measure --dir "C:\\...\\images\\WT"

Smoke-test on one image first and confirm the CSV matches a GUI export.
"""

from __future__ import annotations

import argparse
import csv
import faulthandler
import gc
import json
from datetime import datetime
from pathlib import Path

import numpy as np

# Dump a Python-level traceback if the process dies from a *native* fault
# (access violation / segfault) — those can't be caught by try/except, so
# without this the smoke test just returns silently to the prompt. This makes
# the crashing C call (e.g. a cupy FRET kernel) show up in the terminal.
faulthandler.enable()

from .fret import compute_fret_volume, load_force_lut, resample_fret_outputs
from .io import load_params, load_stack, to_isotropic
from .measure_core import compute_measurements
from .rf import instance_labels_3d

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PARAMS = ROOT / "3d_images_folder" / "params.json"
DEFAULT_LUT = ROOT / "Fret_LUT" / "Clover_mRuby2_GGSGGS7_force_new.txt"

COLS = ["label", "ratio_clover_ruby", "dist_mesh_um", "path_curv_xy",
        "dist_branch_um", "curv_mean", "curv_k1", "curv_k2", "curv_gaussian",
        "dist_nuc_um", "force_mean_pN", "force_std_pN", "force_total_pN",
        "volume_um3", "axis_major_um", "axis_minor_um", "eccentricity", "nuc_id"]

# The post-process settings the batch / batch_postprocess applied (documented in
# the settings JSON so a headless-measured export is still backtrackable).
STD_POSTPROCESS = {"prob_threshold": 0.5, "min_size_vox": 100,
                   "split_merged": True, "split_min_size_vox": 1000,
                   "note": "applied by pipeline.batch / pipeline.batch_postprocess"}


def _p(msg: str) -> None:
    print(f"      . {msg}", flush=True)


def _load_raw(tif: Path, params_path: Path, lut_path: Path):
    """Recompute force_vol + native clover/ruby + scale exactly like the loader."""
    _p("load_params")
    params = load_params(params_path)
    _p(f"load_stack ({tif.name})")
    stack_raw = load_stack(tif, params)
    role_to_c = {role: stack_raw.channels.index(name)
                 for role, name in stack_raw.roles.items()}
    _p("load_force_lut")
    lut_eff, lut_force = load_force_lut(str(lut_path))
    _p("compute_fret_volume (GPU)")
    fret = compute_fret_volume(stack_raw.data, role_to_c, params["fret"],
                               lut_eff, lut_force)
    _p("copy native clover/ruby")
    clover_native = (stack_raw.data[:, role_to_c["clover"]].copy()
                     if "clover" in stack_raw.roles else None)
    ruby_native = (stack_raw.data[:, role_to_c["ruby"]].copy()
                   if "ruby" in stack_raw.roles else None)
    raw_voxel = stack_raw.voxel_zyx_um
    _p("to_isotropic")
    stack = to_isotropic(stack_raw)
    _p("resample_fret_outputs")
    fret = resample_fret_outputs(fret, raw_voxel, stack.voxel_zyx_um[0])
    del stack_raw
    _p("raw load done")
    return (fret["force_vol"], clover_native, ruby_native, raw_voxel,
            tuple(stack.voxel_zyx_um))


def _standard_instances(saved: Path):
    """Recompute adhesion instances from the saved probability map at the
    standard post-process. The on-disk adhesion_instances.npz is often STALE:
    the GUI re-post-processes to standard in memory at measure time but never
    rewrites the npz, so many folders still hold an old, non-standard labeling.
    Recomputing here guarantees every headless CSV is on the exact documented
    segmentation (and reproduces the GUI's standard result). Also rewrites the
    npz so GUI + headless agree from now on."""
    proba_path = saved / "adhesion_proba.npy"
    if not proba_path.exists():
        _p("no adhesion_proba.npy — falling back to saved instances (STALE?)")
        return np.load(saved / "adhesion_instances.npz")["instances"]
    T = int(round(STD_POSTPROCESS["prob_threshold"] * 255))
    _p(f"recompute standard instances from proba (T={T}, min="
       f"{STD_POSTPROCESS['min_size_vox']}, split ON, "
       f"split_min={STD_POSTPROCESS['split_min_size_vox']})")
    proba = np.load(proba_path)
    binary = (proba >= T).astype(np.uint8)
    del proba
    labels = instance_labels_3d(binary, STD_POSTPROCESS["min_size_vox"], True,
                                STD_POSTPROCESS["split_min_size_vox"])
    del binary
    np.savez_compressed(saved / "adhesion_instances.npz",
                        instances=labels.astype(np.int32))
    _p(f"standardized instances: {int(labels.max())} adhesions "
       f"(npz rewritten)")
    return labels


def _measure_one(tif: Path, saved: Path, params_path: Path, lut_path: Path,
                 use_saved: bool = False):
    force_vol, clover_native, ruby_native, raw_voxel, scale = _load_raw(
        tif, params_path, lut_path)
    if use_saved:
        _p("using saved adhesion_instances.npz as-is (--use-saved-instances)")
        labels = np.load(saved / "adhesion_instances.npz")["instances"]
    else:
        labels = _standard_instances(saved)
    vessel_verts = np.load(saved / "vessel_mesh.npz")["verts"]
    dapi_labels = np.load(saved / "dapi_labels.npz")["labels"]

    curvature = None
    if (saved / "vessel_curvature.npz").exists():
        cv = np.load(saved / "vessel_curvature.npz")
        curvature = {k: cv[k] for k in ("mean", "k1", "k2", "gaussian")
                     if k in cv.files}
    centerline = None
    if (saved / "vessel_path_curvature.npz").exists():
        pc = np.load(saved / "vessel_path_curvature.npz")
        centerline = {"pts_yx_um": pc["pts_yx_um"], "kappa": pc["kappa"],
                      "branch_yx_um": pc["branch_yx_um"]}
    fib = {}
    if (saved / "fibroblast_count.npz").exists():
        fc = np.load(saved / "fibroblast_count.npz")
        fib = {k: (int(fc[k]) if k != "threshold" else float(fc[k]))
               for k in fc.files}

    _p("compute_measurements")
    result = compute_measurements(labels, force_vol, vessel_verts, dapi_labels,
                                  scale, clover_native, ruby_native, raw_voxel,
                                  curvature, centerline)
    if result is None:
        print("    no adhesions; nothing written.", flush=True)
        return

    # Build columns in the GUI export order. Write with the stdlib csv module —
    # NOT pandas: pandas pulls in pyarrow, whose native DLL access-violates in
    # this env (crashes the process with no catchable exception).
    metric_cols = [c for c in COLS if c in result]
    cols = {c: np.asarray(result[c]).ravel() for c in metric_cols}
    n = len(next(iter(cols.values()))) if cols else 0
    header = ["image"] + metric_cols + ["n_fibroblasts", "n_nuclei_total"]
    n_fib = fib.get("fibroblast", "")
    n_tot = fib.get("total", "")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv_path = saved / f"{tif.stem}_adhesion_measures_{ts}.csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in range(n):
            row = [tif.name]
            for c in metric_cols:
                v = cols[c][r]
                row.append("" if (isinstance(v, float) and np.isnan(v)) else v)
            row.extend([n_fib, n_tot])
            w.writerow(row)

    meta = {"image": tif.name, "exported": ts, "csv": csv_path.name,
            "n_adhesions": int(n), "scale_um_zyx": list(scale),
            "roi_active": False, "headless_measure": True,
            "fibroblast_count": fib or None,
            "settings": {"adhesion_postprocess": STD_POSTPROCESS}}
    with open(saved / f"{tif.stem}_adhesion_measures_{ts}_settings.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    print(f"    -> {n} adhesions; {csv_path.name} + settings JSON written.",
          flush=True)
    del force_vol, clover_native, ruby_native, labels, vessel_verts, dapi_labels
    gc.collect()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", required=True, type=Path)
    ap.add_argument("--params", type=Path, default=DEFAULT_PARAMS)
    ap.add_argument("--lut", type=Path, default=DEFAULT_LUT)
    ap.add_argument("--redo-json", action="store_true",
                    help="also measure folders that already have a settings JSON")
    ap.add_argument("--only", default=None,
                    help="only process _saved folders whose name contains this "
                         "substring (use for a single-image smoke test)")
    ap.add_argument("--use-saved-instances", action="store_true",
                    help="measure the on-disk adhesion_instances.npz as-is "
                         "instead of recomputing standard instances from the "
                         "probability map. NOT recommended: many folders hold a "
                         "stale, non-standard labeling.")
    args = ap.parse_args()

    targets = []
    for saved in sorted(Path(args.dir).glob("*_saved")):
        if args.only and args.only not in saved.name:
            continue
        if not (saved / "adhesion_instances.npz").exists():
            continue
        if not args.redo_json and list(saved.glob("*_settings.json")):
            continue
        tif = Path(args.dir) / f"{saved.name[:-len('_saved')]}.tif"
        if not tif.exists():
            print(f"  skip ({saved.name}): image {tif.name} not found.", flush=True)
            continue
        targets.append((tif, saved))

    if not targets:
        print("Nothing to measure (need instances present and, unless "
              "--redo-json, no settings JSON).", flush=True)
        return

    print(f"Measuring {len(targets)} image(s):", flush=True)
    for i, (tif, saved) in enumerate(targets):
        print(f"  [{i + 1}/{len(targets)}] {tif.name} ...", flush=True)
        try:
            _measure_one(tif, saved, args.params, args.lut,
                         use_saved=args.use_saved_instances)
        except Exception as exc:
            import traceback
            print(f"    FAILED: {exc}\n{traceback.format_exc()}", flush=True)
        gc.collect()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
