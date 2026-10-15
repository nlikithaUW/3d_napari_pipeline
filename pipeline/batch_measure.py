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
from .io import fret_configured, load_stack, resolve_params, to_isotropic
from .measure_core import (
    CSV_COLUMNS, NUCLEUS_MAP_SPLIT, NUCLEUS_MAP_UNSPLIT,
    compute_measurements, export_meta,
)
from .rf import SPLIT_MIN_DISTANCE_UM, instance_labels_3d, split_distance_vox
from .saved_check import describe_saved, fibroblast_count_notes
from . import saved_io as sio

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PARAMS = ROOT / "3d_images_folder" / "params.json"
DEFAULT_LUT = ROOT / "Fret_LUT" / "Clover_mRuby2_GGSGGS7_force_new.txt"


# The post-process settings the batch / batch_postprocess applied (documented in
# the settings JSON so a headless-measured export is still backtrackable).
STD_POSTPROCESS = {"prob_threshold": 0.5, "min_size_vox": 100,
                   "split_merged": True, "split_min_size_vox": 1000,
                   "split_min_distance_um": SPLIT_MIN_DISTANCE_UM,
                   "note": "applied by pipeline.batch / pipeline.batch_postprocess"}


def _p(msg: str) -> None:
    print(f"      . {msg}", flush=True)


def _load_raw(tif: Path, params_path: Path, lut_path: Path,
              fret_enabled: bool = True):
    """Recompute force_vol + native clover/ruby + scale exactly like the loader.

    ``force_vol`` is None when FRET isn't configured for this dataset — the
    same gate the GUI loader applies (see ``gui.loader.fret_configured``).
    """
    params, params_src, params_note = resolve_params(tif, params_path)
    _p(params_note)
    _p(f"load_stack ({tif.name})")
    stack_raw = load_stack(tif, params)
    role_to_c = {role: stack_raw.channels.index(name)
                 for role, name in stack_raw.roles.items()}
    fret_ok, why = fret_configured(params, stack_raw.roles, fret_enabled)
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
    # For the per-adhesion corrected FRET (measure_core.adhesion_fret), exactly
    # as the GUI loader keeps it: only when FRET is on for this image.
    fret_native = (stack_raw.data[:, role_to_c["fret"]].copy()
                   if fret_ok else None)
    fret_cfg = ((dict(params["fret"]), (lut_eff, lut_force))
                if fret_ok else (None, None))
    raw_voxel = stack_raw.voxel_zyx_um
    _p("to_isotropic")
    stack = to_isotropic(stack_raw)
    if fret is not None:
        _p("resample_fret_outputs")
        fret = resample_fret_outputs(fret, raw_voxel, stack.voxel_zyx_um[0])
    del stack_raw
    _p("raw load done")
    params_used = {"source": str(params_src), "fret": params.get("fret"),
                   "fret_enabled": bool(fret_enabled)}
    return (None if fret is None else fret["force_vol"],
            clover_native, ruby_native, raw_voxel, tuple(stack.voxel_zyx_um),
            params_used, fret_native, fret_cfg)


def _standard_instances(saved: Path, voxel_um, rewrite: bool = False):
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
    proba_path = sio.path(saved, sio.ADHESION_PROBA)
    if not proba_path.exists():
        _p(f"WARN: no {proba_path.name} — measuring the saved instances "
           f"as-is; they may be a stale, non-standard labeling.")
        return sio.load_instances(sio.path(saved, sio.ADHESION_INSTANCES))
    T = int(round(STD_POSTPROCESS["prob_threshold"] * 255))
    split_vox = split_distance_vox(STD_POSTPROCESS["split_min_distance_um"], voxel_um)
    _p(f"recompute standard instances from proba (T={T}, min="
       f"{STD_POSTPROCESS['min_size_vox']}, split ON, "
       f"split_min={STD_POSTPROCESS['split_min_size_vox']}, seed spacing "
       f"{STD_POSTPROCESS['split_min_distance_um']:g} um = {split_vox} vox)")
    proba = sio.load_proba(proba_path)
    binary = (proba >= T).astype(np.uint8)
    del proba
    labels = instance_labels_3d(binary, STD_POSTPROCESS["min_size_vox"], True,
                                STD_POSTPROCESS["split_min_size_vox"], split_vox)
    del binary
    n_new = int(labels.max())
    if rewrite:
        sio.save_instances(sio.path(saved, sio.ADHESION_INSTANCES), labels)
        _p(f"standardized instances: {n_new} adhesions (npz rewritten, "
           f"--rewrite-instances)")
    else:
        inst_path = sio.path(saved, sio.ADHESION_INSTANCES)
        n_old = None
        if inst_path.exists():
            try:
                n_old = int(sio.load_instances(inst_path).max())
            except Exception as exc:
                _p(f"WARN: could not read {inst_path.name} to compare ({exc})")
        if n_old is not None and n_old != n_new:
            _p(f"WARN: the saved {inst_path.name} holds {n_old} "
               f"adhesions but the standard post-process gives {n_new}. "
               f"Measuring the recomputed {n_new}; the file is left alone. "
               f"Pass --rewrite-instances to update it.")
        else:
            _p(f"standardized instances: {n_new} adhesions (npz unchanged)")
    return labels


def _measure_one(tif: Path, saved: Path, params_path: Path, lut_path: Path,
                 use_saved: bool = False, rewrite_instances: bool = False,
                 fret_enabled: bool = True):
    """Measure one image. Returns ``(status, notes)``.

    ``status`` is "ok", "degraded" (a CSV was written but some column is
    silently blank) or "empty" (no adhesions). Exceptions propagate to ``main``,
    which records them as FAILED.
    """
    (force_vol, clover_native, ruby_native, raw_voxel, scale, params_used,
     fret_native, fret_cfg) = _load_raw(tif, params_path, lut_path, fret_enabled)
    if use_saved:
        _p(f"using saved {sio.filename(sio.ADHESION_INSTANCES)} as-is "
           f"(--use-saved-instances)")
        labels = sio.load_instances(sio.path(saved, sio.ADHESION_INSTANCES))
    else:
        labels = _standard_instances(saved, scale, rewrite=rewrite_instances)
    # Faces are what restrict dist_mesh_um to the mesh surface (P09a). A
    # folder whose npz predates them still measures, to every vertex, and
    # warns rather than silently reading short — so the `None` is passed
    # through rather than substituted, because `None` is what triggers that
    # warning. `batch_mesh` substitutes an empty array instead; both are
    # right for their caller, and unifying them would silence P09a here.
    _mesh = sio.load_mesh(sio.path(saved, sio.VESSEL_MESH))
    vessel_verts, vessel_faces = _mesh.verts, _mesh.faces

    # Nucleus map: the watershed-split nuclei when they exist, else the unsplit
    # DAPI components (P04). The fallback is the backward-compatible path for a
    # _saved folder produced before nuclei_split.npz existed — it reproduces
    # that folder's previous numbers exactly and leaves nuc_is_fibroblast blank,
    # which is why it says so loudly rather than quietly doing the old thing.
    nuc_ids = nuc_is_fib = None
    nuc_split = sio.path(saved, sio.NUCLEI_SPLIT)
    if nuc_split.exists():
        ns = sio.load_nuclei_split(nuc_split)
        dapi_labels = ns.labels
        nuc_ids, nuc_is_fib = ns.ids, ns.is_fib
        nucleus_map = NUCLEUS_MAP_SPLIT
        _p(f"nuclei: {np.asarray(nuc_ids).size} split nuclei, "
           f"{int(np.asarray(nuc_is_fib).sum())} VnTs+ "
           f"({nucleus_map})")
    else:
        dapi_labels = sio.load_dapi_labels(
            sio.path(saved, sio.DAPI_LABELS)).labels
        nucleus_map = NUCLEUS_MAP_UNSPLIT
        _p(f"WARN: no {nuc_split.name} — falling back to the UNSPLIT DAPI "
           f"components ({int(dapi_labels.max())} of them) for nuc_id and "
           f"dist_nuc_um, and nuc_is_fibroblast will be blank. These are DAPI "
           f"connected components, not nuclei: most are sub-nuclear debris, "
           f"and nuc_id is not joinable to the fibroblast count. Re-run "
           f"'Count fibroblasts' for this image to fix it.")

    curvature = None
    if sio.path(saved, sio.VESSEL_CURVATURE).exists():
        curvature = sio.load_curvature(sio.path(saved, sio.VESSEL_CURVATURE))
    graph = None
    if sio.path(saved, sio.VESSEL_GRAPH).exists():
        graph = sio.load_graph(sio.path(saved, sio.VESSEL_GRAPH))
    # The ragged wall-shear cloud — one row per wall point, not per centerline
    # point, which is why it is a separate file. Absent unless the graph was
    # built with shear on; `wall_columns` then leaves the tau_* columns blank
    # and says so, rather than the header changing between runs.
    wall = None
    if sio.path(saved, sio.VESSEL_SHEAR).exists():
        wall = sio.load_shear(sio.path(saved, sio.VESSEL_SHEAR))
    fib = {}
    if sio.path(saved, sio.FIBROBLAST_COUNT).exists():
        fib, skipped = sio.load_fibroblast_count(
            sio.path(saved, sio.FIBROBLAST_COUNT))
        for n in fibroblast_count_notes(skipped):
            _p(n)

    _p("compute_measurements")
    result = compute_measurements(labels, force_vol, vessel_verts, dapi_labels,
                                  scale, clover_native, ruby_native, raw_voxel,
                                  curvature, graph, nuc_ids=nuc_ids,
                                  nuc_is_fib=nuc_is_fib,
                                  vessel_faces=vessel_faces, wall=wall,
                                  fret_native=fret_native,
                                  fret_params=fret_cfg[0], fret_lut=fret_cfg[1])
    if result is None:
        print("    no adhesions; nothing written.", flush=True)
        return "empty", []
    notes = list(result.get("notes", []))
    for n in notes:
        _p(n)

    # Build columns in the GUI export order. Write with the stdlib csv module —
    # NOT pandas: pandas pulls in pyarrow, whose native DLL access-violates in
    # this env (crashes the process with no catchable exception).
    metric_cols = [c for c in CSV_COLUMNS if c in result]
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

    meta = export_meta(
        image=tif.name, exported=ts, csv_name=csv_path.name, n_adhesions=n,
        scale_um_zyx=scale, roi_active=False, fibroblast_count=fib or None,
        nucleus_map=nucleus_map, notes=notes,
        settings={"adhesion_postprocess": STD_POSTPROCESS,
                  "params": params_used}, headless=True)
    with open(saved / f"{tif.stem}_adhesion_measures_{ts}_settings.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    print(f"    -> {n} adhesions; {csv_path.name} + settings JSON written.",
          flush=True)
    del force_vol, clover_native, ruby_native, labels, vessel_verts
    del vessel_faces, dapi_labels
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
                    help="measure the on-disk instances npz as-is "
                         "instead of recomputing standard instances from the "
                         "probability map. NOT recommended: a folder may hold a "
                         "stale, non-standard labeling.")
    ap.add_argument("--rewrite-instances", action="store_true",
                    help="write the recomputed standard instances back to "
                         "the instances npz. Off by default: measuring "
                         "should not modify a saved segmentation.")
    ap.add_argument("--check-only", action="store_true",
                    help="describe every *_saved folder and exit without "
                         "measuring. The pre-flight for an unattended run: it "
                         "reports which columns would come out blank, before "
                         "the compute rather than after.")
    ap.add_argument("--no-fret", action="store_true",
                    help="skip FRET even where the experiment's params.json "
                         "has a fret block; force_*_pN export blank (the Load "
                         "tab's 'compute FRET' checkbox, unticked)")
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
        if not sio.path(saved, sio.ADHESION_INSTANCES).exists():
            print(f"  skip ({saved.name}): no "
                  f"{sio.filename(sio.ADHESION_INSTANCES)} — this "
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
                rewrite_instances=args.rewrite_instances,
                fret_enabled=not args.no_fret)
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
