"""Consistency checks over the saved intermediates — one copy, both paths.

Every function here is pure: numpy in, ``list[str]`` of human-readable notes
out. No Qt, no napari, no ``state`` dict, no printing. Callers print at their own
indent, which is why nothing here knows what indent it is at.

**Why this is a module and not inline checks.** The GUI grew nine of these checks
inline (a stale curvature file, a mismatched label volume, a missing channel) and
the headless path grew none — which is exactly the asymmetry `docs/ISSUES.md`
[P12] names: an unattended run produces a CSV that looks fine with four blank
columns and no message. Writing them inline a second time would recreate the
defect one release later. A shared module means a GUI warning and a batch warning
are *the same string*, and adding a check adds it to both.

The notes describe **silent degradations** — cases where the pipeline produces a
CSV that is structurally valid but has columns that are quietly all-NaN. They are
warnings, not errors: a legitimately FRET-less dataset must not fail a run.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

CURV_KEYS = ("mean", "k1", "k2", "gaussian")


# --------------------------------------------------------------------------- #
# Per-degradation checks. Each answers: "which columns come out blank, and why?"
# --------------------------------------------------------------------------- #

def curvature_notes(curvature, nverts: int) -> list[str]:
    """Notes for the four ``curv_*`` columns against a mesh of ``nverts``.

    Mirrors the gates in ``measure_core``: a missing file, a length mismatch, or
    an all-NaN ``k1`` blanks all four; the per-key check blanks one column while
    the others survive, which is the case that is hardest to spot in a CSV.
    """
    if curvature is None:
        return ["WARN: no surface curvature loaded — curv_mean, curv_k1, "
                "curv_k2 and curv_gaussian will all be blank."]

    k1 = curvature.get("k1")
    k1 = None if k1 is None else np.asarray(k1)
    if k1 is None:
        return ["WARN: the curvature file has no 'k1' array — curv_mean, "
                "curv_k1, curv_k2 and curv_gaussian will all be blank."]

    if k1.shape[0] != int(nverts):
        # Curvature is per-vertex, so it is only meaningful against the mesh it
        # was computed on. This is the stale-pairing case P12 was opened for.
        return [f"WARN: curvature has {k1.shape[0]} values but the mesh has "
                f"{int(nverts)} vertices — this file is from a different mesh. "
                f"curv_mean, curv_k1, curv_k2 and curv_gaussian will be blank; "
                f"recompute the curvature (or reload the matching mesh)."]

    if not bool(np.isfinite(k1).any()):
        return ["WARN: surface curvature is entirely non-finite (no surface "
                "vertices) — curv_mean, curv_k1, curv_k2 and curv_gaussian "
                "will all be blank."]

    missing = []
    for key in CURV_KEYS:
        cv = curvature.get(key)
        cv = None if cv is None else np.asarray(cv)
        if cv is None or cv.shape[0] != int(nverts):
            missing.append(f"curv_{key}")
    if missing:
        return [f"WARN: {', '.join(missing)} will be blank — missing from the "
                f"curvature file or the wrong length for a {int(nverts)}-vertex "
                f"mesh. The other curvature columns are unaffected."]
    return []


def centerline_notes(centerline) -> list[str]:
    """Notes for ``path_curv_xy`` and ``dist_branch_um``."""
    if centerline is None or not len(centerline.get("pts_yx_um", [])):
        return ["WARN: no vessel centerline — path_curv_xy and dist_branch_um "
                "will both be blank."]
    if not len(np.asarray(centerline.get("branch_yx_um", []))):
        return ["WARN: the centerline has no branch nodes — dist_branch_um "
                "will be blank. Expected on a non-branched vessel (a single "
                "tunnel or a cone); note that a bare `dist_branch_um > X` "
                "filter then drops every row rather than none."]
    return []


def native_notes(clover_native, ruby_native, labels_shape) -> list[str]:
    """Notes for ``ratio_clover_ruby`` — the project's primary output.

    The shape gate in ``measure_core`` is silent, so a mismatch here blanks the
    tension column for *every* adhesion with nothing in the log.
    """
    missing = [n for n, v in (("clover", clover_native), ("ruby", ruby_native))
               if v is None]
    if missing:
        return [f"WARN: no native {' and '.join(missing)} channel — "
                f"ratio_clover_ruby, the primary tension output, will be blank "
                f"for every adhesion."]

    c_xy = tuple(np.asarray(clover_native).shape[1:])
    l_xy = tuple(labels_shape[1:])
    if c_xy != l_xy:
        return [f"WARN: native clover XY {c_xy} != adhesion labels XY {l_xy} — "
                f"ratio_clover_ruby, the primary tension output, will be blank "
                f"for every adhesion. The saved instances and the raw image do "
                f"not match; check they are the same acquisition."]

    r_shape = tuple(np.asarray(ruby_native).shape)
    if r_shape != tuple(np.asarray(clover_native).shape):
        return [f"WARN: native ruby {r_shape} != native clover "
                f"{tuple(np.asarray(clover_native).shape)} — ratio_clover_ruby "
                f"cannot be computed from mismatched channels."]
    return []


def ratio_dropout_note(n_nan: int, n_total: int) -> list[str]:
    """Note for per-adhesion ratio dropout (``r_sum <= 0`` → NaN).

    An adhesion thinner than one native z-step can fall entirely between sampled
    planes (`docs/ISSUES.md` P15, `measure_core.py`). That is a size-dependent
    and therefore depth-dependent dropout, so a bare count makes an otherwise
    invisible depth effect visible.
    """
    if n_nan <= 0 or n_total <= 0:
        return []
    return [f"WARN: ratio_clover_ruby is blank for {n_nan}/{n_total} adhesions "
            f"({100.0 * n_nan / n_total:.1f}%) whose ruby is entirely at or "
            f"below background. This dropout is size- and therefore "
            f"depth-dependent."]


def force_notes(force_vol) -> list[str]:
    """Note for the three ``force_*_pN`` columns."""
    if force_vol is None:
        return ["WARN: no force volume (FRET not configured) — force_mean_pN, "
                "force_std_pN and force_total_pN will be blank. Expected on "
                "this dataset; the uncalibrated mClover/mRuby ratio is the "
                "tension proxy."]
    return []


def nuclei_notes(nuc_ids, nuc_is_fib, nuc_labels_max: int) -> list[str]:
    """Notes for ``nuc_is_fibroblast`` — the P04 join."""
    if nuc_ids is None or nuc_is_fib is None:
        return ["WARN: no per-nucleus fibroblast classification available — "
                "nuc_is_fibroblast will be blank. Re-run 'Count fibroblasts' "
                "to write nuclei_split.npz for this image."]
    ids = np.asarray(nuc_ids)
    is_fib = np.asarray(nuc_is_fib)
    if ids.size != is_fib.size:
        return [f"WARN: nucleus ids ({ids.size}) and classes ({is_fib.size}) "
                f"have different lengths — nuc_is_fibroblast will be blank. "
                f"nuclei_split.npz is inconsistent; recompute it."]
    if ids.size and int(ids.max()) > int(nuc_labels_max):
        return [f"WARN: nucleus ids run to {int(ids.max())} but the nucleus "
                f"label volume only reaches {int(nuc_labels_max)} — "
                f"nuc_is_fibroblast will be blank. The ids and the label "
                f"volume are from different runs."]
    return []


def mesh_notes(vessel_verts) -> list[str]:
    """Notes for ``dist_mesh_um`` and everything read off the mesh."""
    v = np.asarray(vessel_verts)
    if v.size == 0 or v.shape[0] == 0:
        return ["WARN: the vessel mesh has no vertices — dist_mesh_um and all "
                "four curv_* columns are undefined."]
    if v.ndim != 2 or v.shape[1] != 3:
        return [f"WARN: vessel mesh vertices have shape {v.shape}, expected "
                f"(N, 3) — dist_mesh_um cannot be computed."]
    return []


def volume_shape_notes(name: str, shape, ref_name: str, ref_shape) -> list[str]:
    """Generic 'these two volumes should be the same image' check."""
    a, b = tuple(shape), tuple(ref_shape)
    if a != b:
        return [f"WARN: {name} is {a} but {ref_name} is {b} — they are not the "
                f"same image, or one was saved at a different resampling."]
    return []


# --------------------------------------------------------------------------- #
# Pre-flight description of a <stem>_saved folder
# --------------------------------------------------------------------------- #

# Files described by presence + size only: reading them means decompressing
# hundreds of MB, and nothing here needs their contents.
_BULK = ("adhesion_proba.npy", "adhesion_instances.npz", "dapi_labels.npz",
         "cadherin_labels.npz")


def _size(p: Path) -> str:
    n = p.stat().st_size
    return f"{n / 1e9:.2f} GB" if n >= 1e9 else f"{n / 1e6:.1f} MB"


def _member_shapes(path: Path) -> dict:
    """{key: (shape, dtype)} for an npz, WITHOUT decompressing the payloads.

    An npz is a zip, so ``d[key]`` inflates the whole member. Each member is a
    .npy whose header is the first ~100 bytes, and a deflate stream can be read
    that far and abandoned — so this is O(number of arrays), not O(bytes).
    """
    out = {}
    try:
        with np.load(path) as d:
            zf = d.zip
            for key in d.files:
                try:
                    with zf.open(f"{key}.npy") as fh:
                        version = np.lib.format.read_magic(fh)
                        if version == (1, 0):
                            shape, _, dt = np.lib.format.read_array_header_1_0(fh)
                        else:
                            shape, _, dt = np.lib.format.read_array_header_2_0(fh)
                    out[key] = (tuple(shape), np.dtype(dt))
                except Exception:
                    out[key] = None
    except Exception:
        return {}
    return out


def _scalars(path: Path, keys) -> dict:
    """Read a few small members by name; silently skip what is not there."""
    out = {}
    try:
        with np.load(path) as d:
            for k in keys:
                if k in d.files:
                    out[k] = d[k]
    except Exception:
        pass
    return out


def describe_saved(path) -> list[str]:
    """One block of report lines for a ``<stem>_saved`` folder.

    The pre-flight for an unattended run: it answers "will measuring this folder
    produce a full CSV?" before hours of compute, not after. Cheap by
    construction — only small arrays are ever decompressed.
    """
    saved = Path(path)
    lines = [f"{saved.name}"]
    if not saved.is_dir():
        return lines + ["  FAILED: not a directory"]

    present = {p.name: p for p in saved.iterdir() if p.is_file()}

    # --- the bulk artifacts: presence + size only ---------------------------
    for name in _BULK:
        p = present.get(name)
        lines.append(f"  {name:<26} {'—  MISSING' if p is None else _size(p)}")

    # --- vessel mesh + curvature: the pair P12 is about ---------------------
    nverts = None
    mp = present.get("vessel_mesh.npz")
    if mp is None:
        lines.append(f"  {'vessel_mesh.npz':<26} —  MISSING")
    else:
        shapes = _member_shapes(mp)
        vs = shapes.get("verts")
        nverts = vs[0][0] if vs else None
        nfaces = shapes.get("faces")
        lines.append(f"  {'vessel_mesh.npz':<26} {nverts} verts, "
                     f"{nfaces[0][0] if nfaces else '?'} faces")

    cp = present.get("vessel_curvature.npz")
    if cp is None:
        lines.append(f"  {'vessel_curvature.npz':<26} —  MISSING")
    else:
        shapes = _member_shapes(cp)
        k1 = shapes.get("k1")
        nk = k1[0][0] if k1 else None
        lines.append(f"  {'vessel_curvature.npz':<26} {nk} values "
                     f"({', '.join(sorted(shapes))})")
        if nverts is not None and nk is not None:
            # Stand-ins: only the LENGTH of these matters to the per-key check,
            # and zeros keep it deterministic. k1 is the one whose values are
            # tested (all-NaN => no surface vertices), so k1 is read for real.
            curv = {k: np.zeros(nk, dtype=np.float32) for k in shapes
                    if k in CURV_KEYS}
            # Length is what the check needs; finiteness needs the real k1, and
            # k1 is one small array, so read just that one.
            if nk == nverts:
                real = _scalars(cp, ("k1",))
                if "k1" in real:
                    curv["k1"] = real["k1"]
            for n in curvature_notes(curv, nverts):
                lines.append(f"  {n}")

    pc = present.get("vessel_path_curvature.npz")
    if pc is None:
        lines.append(f"  {'vessel_path_curvature.npz':<26} —  MISSING")
        lines.append("  " + centerline_notes(None)[0])
    else:
        d = _scalars(pc, ("pts_yx_um", "kappa", "branch_yx_um"))
        npts = len(d.get("pts_yx_um", []))
        nbr = len(d.get("branch_yx_um", []))
        lines.append(f"  {'vessel_path_curvature.npz':<26} {npts} centerline "
                     f"points, {nbr} branch node(s)")
        for n in centerline_notes({"pts_yx_um": d.get("pts_yx_um", []),
                                   "branch_yx_um": d.get("branch_yx_um", [])}):
            lines.append(f"  {n}")

    # --- nuclei: the P04 join -----------------------------------------------
    np_ = present.get("nuclei_split.npz")
    if np_ is None:
        lines.append("  WARN: no nuclei_split.npz — nuc_id falls back to the "
                     "unsplit DAPI components and nuc_is_fibroblast will be "
                     "blank. Re-run 'Count fibroblasts' for this image.")
    else:
        d = _scalars(np_, ("ids", "is_fib", "min_size_vox", "threshold"))
        ids = np.asarray(d.get("ids", []))
        is_fib = np.asarray(d.get("is_fib", []))
        shapes = _member_shapes(np_)
        lab = shapes.get("labels")
        lines.append(f"  {'nuclei_split.npz':<26} {ids.size} nuclei, "
                     f"{int(is_fib.sum()) if is_fib.size else 0} VnTs+, "
                     f"labels {lab[0] if lab else '?'} "
                     f"{lab[1] if lab else ''}")
        lab_max = int(ids.max()) if ids.size else 0
        for n in nuclei_notes(ids, is_fib, lab_max):
            lines.append(f"  {n}")

    fp = present.get("fibroblast_count.npz")
    if fp is None:
        lines.append(f"  {'fibroblast_count.npz':<26} —  MISSING "
                     f"(n_fibroblasts / n_nuclei_total will be blank)")
    else:
        d = _scalars(fp, ("fibroblast", "total"))
        lines.append(f"  {'fibroblast_count.npz':<26} "
                     f"{int(d.get('fibroblast', 0))} of "
                     f"{int(d.get('total', 0))} nuclei")

    # --- what has already been measured -------------------------------------
    csvs = sorted(p.name for p in saved.glob("*_adhesion_measures_*.csv"))
    jsons = sorted(saved.glob("*_settings.json"))
    lines.append(f"  {'measured':<26} {len(csvs)} CSV(s), {len(jsons)} "
                 f"settings JSON"
                 + (f"; latest {csvs[-1]}" if csvs else ""))

    if "adhesion_instances.npz" not in present:
        lines.append("  WARN: no adhesion_instances.npz — this folder cannot "
                     "be measured at all; re-run the batch for this image.")
    return lines
