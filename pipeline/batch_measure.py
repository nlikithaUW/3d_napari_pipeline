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

Check what a folder would produce before spending the compute on it:

    python -m pipeline.batch_measure --dir "..." --check-only

Exit codes: 0 all measured, 1 at least one image FAILED, 2 (with --strict) a CSV
was written but some column is silently blank.

Smoke-test on one image first and confirm the CSV matches a GUI export.
"""

from __future__ import annotations

import argparse
import csv
import faulthandler
import gc
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

# Dump a Python-level traceback if the process dies from a *native* fault
# (access violation / segfault) — those can't be caught by try/except, so
# without this the smoke test just returns silently to the prompt. This makes
# the crashing C call (e.g. a cupy FRET kernel) show up in the terminal.
faulthandler.enable()

from .fret import compute_fret_volume, load_force_lut, resample_fret_outputs
from .io import fret_configured, load_params, load_stack, to_isotropic
from .measure_core import compute_measurements
from .rf import instance_labels_3d
from .saved_check import describe_saved

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PARAMS = ROOT / "3d_images_folder" / "params.json"
DEFAULT_LUT = ROOT / "Fret_LUT" / "Clover_mRuby2_GGSGGS7_force_new.txt"

COLS = ["label", "ratio_clover_ruby", "dist_mesh_um",
        "centroid_z_um", "centroid_y_um", "centroid_x_um", "path_curv_xy",
        "dist_branch_um", "curv_mean", "curv_k1", "curv_k2", "curv_gaussian",
        "dist_nuc_um", "force_mean_pN", "force_std_pN", "force_total_pN",
        "volume_um3", "axis_major_um", "axis_minor_um", "eccentricity", "nuc_id",
        "nuc_is_fibroblast"]

# The post-process settings the batch / batch_postprocess applied (documented in
# the settings JSON so a headless-measured export is still backtrackable).
STD_POSTPROCESS = {"prob_threshold": 0.5, "min_size_vox": 100,
                   "split_merged": True, "split_min_size_vox": 1000,
                   "note": "applied by pipeline.batch / pipeline.batch_postprocess"}


def _p(msg: str) -> None:
    print(f"      . {msg}", flush=True)


def _load_raw(tif: Path, params_path: Path, lut_path: Path):
    """Recompute force_vol + native clover/ruby + scale exactly like the loader.

    ``force_vol`` is None when FRET isn't configured for this dataset — the
    same gate the GUI loader applies (see ``gui.loader.fret_configured``).
    """
    _p("load_params")
    params = load_params(params_path)
    _p(f"load_stack ({tif.name})")
    stack_raw = load_stack(tif, params)
    role_to_c = {role: stack_raw.channels.index(name)
                 for role, name in stack_raw.roles.items()}
    fret_ok, why = fret_configured(params, stack_raw.roles)
    if fret_ok:
        _p(f"load_force_lut ({why})")
        lut_eff, lut_force = load_force_lut(str(lut_path))
        _p("compute_fret_volume (GPU)")
        fret = compute_fret_volume(stack_raw.data, role_to_c, params["fret"],
                                   lut_eff, lut_force)
    else:
        _p(f"FRET skipped — {why}; force_*_pN columns will be blank")
        fret = None
    _p("copy native clover/ruby")
    clover_native = (stack_raw.data[:, role_to_c["clover"]].copy()
                     if "clover" in stack_raw.roles else None)
    ruby_native = (stack_raw.data[:, role_to_c["ruby"]].copy()
                   if "ruby" in stack_raw.roles else None)
    raw_voxel = stack_raw.voxel_zyx_um
    _p("to_isotropic")
    stack = to_isotropic(stack_raw)
    if fret is not None:
        _p("resample_fret_outputs")
        fret = resample_fret_outputs(fret, raw_voxel, stack.voxel_zyx_um[0])
    del stack_raw
    _p("raw load done")
    return (None if fret is None else fret["force_vol"],
            clover_native, ruby_native, raw_voxel, tuple(stack.voxel_zyx_um))


def _standard_instances(saved: Path, rewrite: bool = False):
    """Recompute adhesion instances from the saved probability map at the
    standard post-process. The on-disk adhesion_instances.npz can be STALE: the
    GUI re-post-processes to standard in memory at measure time but never
    rewrites the npz, so an old, non-standard labeling can persist. Recomputing
    here guarantees every headless CSV is on the exact documented segmentation
    (and reproduces the GUI's standard result).

    **Rewriting the npz is opt-in** (`--rewrite-instances`). A command named
    "measure" silently replacing a saved segmentation is the second half of
    `docs/ISSUES.md` P12. Without the flag the recomputed labels are compared
    against the file and a difference is reported, not applied. For an image
    freshly produced by `pipeline.batch` this is a no-op either way:
    `run_rf_postprocess` forces `split_merged=True` and the widget defaults are
    byte-identical to STD_POSTPROCESS, so the file already holds this labeling.
    """
    proba_path = saved / "adhesion_proba.npy"
    if not proba_path.exists():
        _p("WARN: no adhesion_proba.npy — measuring the saved instances as-is; "
           "they may be a stale, non-standard labeling.")
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
    n_new = int(labels.max())
    if rewrite:
        np.savez_compressed(saved / "adhesion_instances.npz",
                            instances=labels.astype(np.int32))
        _p(f"standardized instances: {n_new} adhesions (npz rewritten, "
           f"--rewrite-instances)")
    else:
        inst_path = saved / "adhesion_instances.npz"
        n_old = None
        if inst_path.exists():
            try:
                n_old = int(np.load(inst_path)["instances"].max())
            except Exception as exc:
                _p(f"WARN: could not read {inst_path.name} to compare ({exc})")
        if n_old is not None and n_old != n_new:
            _p(f"WARN: the saved adhesion_instances.npz holds {n_old} "
               f"adhesions but the standard post-process gives {n_new}. "
               f"Measuring the recomputed {n_new}; the file is left alone. "
               f"Pass --rewrite-instances to update it.")
        else:
            _p(f"standardized instances: {n_new} adhesions (npz unchanged)")
    return labels


def _measure_one(tif: Path, saved: Path, params_path: Path, lut_path: Path,
                 use_saved: bool = False, rewrite_instances: bool = False):
    """Measure one image. Returns ``(status, notes)``.

    ``status`` is "ok", "degraded" (a CSV was written but some column is
    silently blank) or "empty" (no adhesions). Exceptions propagate to ``main``,
    which records them as FAILED.
    """
    force_vol, clover_native, ruby_native, raw_voxel, scale = _load_raw(
        tif, params_path, lut_path)
    if use_saved:
        _p("using saved adhesion_instances.npz as-is (--use-saved-instances)")
        labels = np.load(saved / "adhesion_instances.npz")["instances"]
    else:
        labels = _standard_instances(saved, rewrite=rewrite_instances)
    vessel_verts = np.load(saved / "vessel_mesh.npz")["verts"]

    # Nucleus map: the watershed-split nuclei when they exist, else the unsplit
    # DAPI components (P04). The fallback is the backward-compatible path for a
    # _saved folder produced before nuclei_split.npz existed — it reproduces
    # that folder's previous numbers exactly and leaves nuc_is_fibroblast blank,
    # which is why it says so loudly rather than quietly doing the old thing.
    nuc_ids = nuc_is_fib = None
    nuc_split = saved / "nuclei_split.npz"
    if nuc_split.exists():
        ns = np.load(nuc_split)
        dapi_labels = ns["labels"]
        nuc_ids, nuc_is_fib = ns["ids"], ns["is_fib"]
        nucleus_map = "nuclei_split.npz"
        _p(f"nuclei: {np.asarray(nuc_ids).size} split nuclei, "
           f"{int(np.asarray(nuc_is_fib).sum())} VnTs+ "
           f"({nucleus_map})")
    else:
        dapi_labels = np.load(saved / "dapi_labels.npz")["labels"]
        nucleus_map = "dapi_labels.npz (unsplit fallback)"
        _p(f"WARN: no nuclei_split.npz — falling back to the UNSPLIT DAPI "
           f"components ({int(dapi_labels.max())} of them) for nuc_id and "
           f"dist_nuc_um, and nuc_is_fibroblast will be blank. These are DAPI "
           f"connected components, not nuclei: most are sub-nuclear debris, "
           f"and nuc_id is not joinable to the fibroblast count. Re-run "
           f"'Count fibroblasts' for this image to fix it.")

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
        # `int(fc[k])` raises TypeError on any non-scalar member, which would
        # turn a slightly-different file into a hard failure of the whole image.
        # Take the 0-d members and report the rest rather than dying on them.
        for k in fc.files:
            v = fc[k]
            if np.ndim(v) != 0:
                _p(f"WARN: fibroblast_count.npz['{k}'] has shape "
                   f"{np.shape(v)}, expected a scalar — ignored.")
                continue
            fib[k] = float(v) if k == "threshold" else int(v)

    _p("compute_measurements")
    result = compute_measurements(labels, force_vol, vessel_verts, dapi_labels,
                                  scale, clover_native, ruby_native, raw_voxel,
                                  curvature, centerline, nuc_ids=nuc_ids,
                                  nuc_is_fib=nuc_is_fib)
    if result is None:
        print("    no adhesions; nothing written.", flush=True)
        return "empty", []
    notes = list(result.get("notes", []))
    for n in notes:
        _p(n)

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
                # np.float32 is NOT a subclass of float, so an isinstance-only
                # check would leave the float32 columns (ratio, curvature)
                # writing the literal string "nan" while float64 ones blanked.
                row.append("" if (np.issubdtype(type(v), np.floating)
                                  and np.isnan(v)) else v)
            row.extend([n_fib, n_tot])
            w.writerow(row)

    meta = {"image": tif.name, "exported": ts, "csv": csv_path.name,
            "n_adhesions": int(n), "scale_um_zyx": list(scale),
            "roi_active": False, "headless_measure": True,
            "fibroblast_count": fib or None,
            # Which nucleus segmentation nuc_id / dist_nuc_um refer to. The two
            # are NOT comparable: the unsplit map's ids index DAPI components,
            # the split map's index nuclei (P04).
            "nucleus_map": nucleus_map,
            "notes": notes,
            "settings": {"adhesion_postprocess": STD_POSTPROCESS}}
    with open(saved / f"{tif.stem}_adhesion_measures_{ts}_settings.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    print(f"    -> {n} adhesions; {csv_path.name} + settings JSON written.",
          flush=True)
    del force_vol, clover_native, ruby_native, labels, vessel_verts, dapi_labels
    gc.collect()
    return ("degraded" if notes else "ok"), notes


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
                         "probability map. NOT recommended: a folder may hold a "
                         "stale, non-standard labeling.")
    ap.add_argument("--rewrite-instances", action="store_true",
                    help="write the recomputed standard instances back to "
                         "adhesion_instances.npz. Off by default: measuring "
                         "should not modify a saved segmentation.")
    ap.add_argument("--check-only", action="store_true",
                    help="describe every *_saved folder and exit without "
                         "measuring. The pre-flight for an unattended run: it "
                         "reports which columns would come out blank, before "
                         "the compute rather than after.")
    ap.add_argument("--strict", action="store_true",
                    help="exit 2 if any image is degraded (a CSV written with "
                         "a silently blank column), not just on failure.")
    args = ap.parse_args()

    if args.check_only:
        folders = [s for s in sorted(Path(args.dir).glob("*_saved"))
                   if not args.only or args.only in s.name]
        if not folders:
            print(f"No *_saved folders in {args.dir}.", flush=True)
            sys.exit(1)
        print(f"Pre-flight over {len(folders)} saved folder(s) in "
              f"{args.dir}:\n", flush=True)
        n_warn = 0
        for saved in folders:
            lines = describe_saved(saved)
            n_warn += sum(1 for ln in lines if "WARN:" in ln)
            print("\n".join(lines) + "\n", flush=True)
        print(f"{len(folders)} folder(s), {n_warn} warning(s). Nothing "
              f"measured (--check-only).", flush=True)
        sys.exit(2 if (n_warn and args.strict) else 0)

    targets = []
    for saved in sorted(Path(args.dir).glob("*_saved")):
        if args.only and args.only not in saved.name:
            continue
        if not (saved / "adhesion_instances.npz").exists():
            print(f"  skip ({saved.name}): no adhesion_instances.npz — this "
                  f"image has not been through pipeline.batch.", flush=True)
            continue
        if not args.redo_json and list(saved.glob("*_settings.json")):
            print(f"  skip ({saved.name}): already measured (a settings JSON "
                  f"exists). Pass --redo-json to measure it again.", flush=True)
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

    # Per-image outcome, so an unattended run ends with a verdict rather than a
    # bare "Done." after a screen of scrolled-past output.
    outcome: dict[str, list] = {"ok": [], "degraded": [], "failed": [],
                                "empty": []}
    degraded_notes: list[tuple[str, list[str]]] = []

    print(f"Measuring {len(targets)} image(s):", flush=True)
    for i, (tif, saved) in enumerate(targets):
        print(f"  [{i + 1}/{len(targets)}] {tif.name} ...", flush=True)
        try:
            status, notes = _measure_one(
                tif, saved, args.params, args.lut,
                use_saved=args.use_saved_instances,
                rewrite_instances=args.rewrite_instances)
            outcome[status].append(tif.name)
            if notes:
                degraded_notes.append((tif.name, notes))
        except Exception as exc:
            import traceback
            print(f"    FAILED: {exc}\n{traceback.format_exc()}", flush=True)
            outcome["failed"].append(tif.name)
        gc.collect()

    print(f"\n{len(outcome['ok'])} ok, {len(outcome['degraded'])} degraded, "
          f"{len(outcome['failed'])} failed, {len(outcome['empty'])} with no "
          f"adhesions.", flush=True)
    for name, notes in degraded_notes:
        print(f"\n  {name}:", flush=True)
        for n in notes:
            print(f"    {n}", flush=True)
    for name in outcome["failed"]:
        print(f"\n  FAILED: {name}", flush=True)

    if outcome["failed"]:
        sys.exit(1)
    # A degradation is a warning, not a failure: a legitimately FRET-less
    # dataset (this one) must not fail the run. --strict is for the case where
    # a complete CSV is the requirement.
    if args.strict and degraded_notes:
        print("\n--strict: exiting 2 because some columns are blank.",
              flush=True)
        sys.exit(2)


if __name__ == "__main__":
    main()
