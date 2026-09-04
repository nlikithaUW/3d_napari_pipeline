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

import json
from pathlib import Path

import numpy as np

from .saved_io import BULK as _BULK
from .saved_io import CURV_KEYS
from . import saved_io as sio

# `CURV_KEYS` and `_BULK` now live in `saved_io`, the one owner of what a
# `_saved` folder holds; they are re-exported here because this module's own
# callers already import them from it.


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


def graph_notes(graph) -> list[str]:
    """Notes for the columns the 3D centerline graph feeds.

    Replaces the old ``centerline_notes``. The columns are the same two —
    ``path_curv_xy`` and ``dist_branch_um`` — but they are now read off an
    ordered graph rather than a raster-order point set, so the failure modes
    worth warning about are different ones.
    """
    if graph is None or not len(np.asarray(graph.get("points_zyx_um", []))):
        return ["WARN: no vessel centerline graph — path_curv_xy, "
                "dist_branch_um and every caliber column will be blank."]
    out = []
    # Counted from the node arrays, not read from a scalar. The saved
    # `vessel_graph.npz` holds arrays only — the scalars live in
    # vessel_recon.json — so a `graph.get("n_junctions", 0)` default of 0
    # silently reported "no junction" on a branched vessel loaded from disk.
    kind = np.asarray(graph.get("node_kind", []))
    ek = np.asarray(graph.get("node_endpoint_kind", []))
    n_junc = (int((kind == "junction").sum()) if kind.size
              else int(graph.get("n_junctions", 0)))
    n_term = (int((ek == "terminus").sum()) if ek.size
              else int(graph.get("n_terminus_endpoints", 0)))
    if n_junc == 0:
        out.append(
            "WARN: the graph has no junction — dist_branch_um will be blank. "
            "Expected on a non-branched vessel (three of the six images); note "
            "that a bare `dist_branch_um > X` filter then drops every row "
            "rather than none (docs/ISSUES.md P08(d)).")
    if n_term:
        out.append(
            f"WARN: {n_term} centerline endpoint(s) stop inside the volume "
            f"rather than at a field-of-view opening. No image in this dataset "
            f"is expected to contain a real vessel terminus "
            f"(docs/DIRECTION.md §5.1), so this usually means the lumen was "
            f"reconstructed short — check before trusting the caliber near it.")
    sec = np.asarray(graph.get("area_um2", []), float)
    if sec.size:
        finite = np.isfinite(sec)
        if not finite.any():
            out.append("WARN: no cross-section was measurable anywhere on the "
                       "centerline — every caliber column is blank.")
        else:
            clip = np.asarray(graph.get("section_clipped", []), bool)
            if clip.size == sec.size and finite.sum():
                frac = float(clip[finite].mean())
                if frac > 0.25:
                    out.append(
                        f"WARN: {frac:.0%} of cross-sections ran into the edge "
                        f"of the sampling disc or of the volume, so their "
                        f"width and area are lower bounds rather than "
                        f"measurements. Raise --section-span, or treat this "
                        f"image's caliber as censored.")
    return out


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
                f"to write {sio.filename(sio.NUCLEI_SPLIT)} for this image."]
    ids = np.asarray(nuc_ids)
    is_fib = np.asarray(nuc_is_fib)
    if ids.size != is_fib.size:
        return [f"WARN: nucleus ids ({ids.size}) and classes ({is_fib.size}) "
                f"have different lengths — nuc_is_fibroblast will be blank. "
                f"{sio.filename(sio.NUCLEI_SPLIT)} is inconsistent; "
                f"recompute it."]
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


def surface_verts_for_distance(vessel_verts, vessel_faces):
    """The vertex set ``dist_mesh_um`` should measure to, plus any notes.

    Returns ``(verts, notes)``. ``verts`` is the subset of ``vessel_verts``
    referenced by at least one face — the mesh *surface* — or all of them when
    the faces are unavailable or unusable.

    **Why this is not just ``vessel_verts``** (`docs/ISSUES.md` P09(a)). The
    default alpha-3D method returns the boundary of a *solid*, so many vertices
    lie strictly inside it and belong to no face. Measuring an adhesion's
    distance to the nearest such vertex measures to a point buried in the wall,
    not to the wall. Measured 2026-08-12 across all six images, interior
    vertices are 18.8% of the hand-tuned reference mesh but **73–85%** of the
    batched ones, where they capture 39–75% of all nearest-vertex lookups.

    This deliberately uses ``np.unique(faces)`` rather than the ``surf_mask``
    the curvature block derives from ``np.isfinite(k1)``. That mask means
    "curvature was computable here" and is absent whenever the curvature file is
    missing or stale — a degradation P12 reports — so reusing it would silently
    couple ``dist_mesh_um`` to a file it has nothing to do with. Faces ship with
    the mesh and are always present.
    """
    v = np.asarray(vessel_verts)
    if vessel_faces is None:
        return v, ["WARN: no mesh faces available, so dist_mesh_um is measured "
                   "to every mesh vertex including interior ones "
                   "(docs/ISSUES.md P09a) — it will read short."]
    f = np.asarray(vessel_faces)
    if f.size == 0 or f.ndim != 2 or f.shape[1] != 3:
        return v, [f"WARN: mesh faces have shape {f.shape}, expected (M, 3) — "
                   f"dist_mesh_um falls back to every vertex, including "
                   f"interior ones (docs/ISSUES.md P09a)."]
    surf = np.unique(f)
    if surf.size == 0 or int(surf.max()) >= v.shape[0]:
        return v, [f"WARN: mesh faces index up to vertex {int(surf.max())} but "
                   f"there are {v.shape[0]} vertices — faces and verts are not "
                   f"from the same mesh; dist_mesh_um falls back to every "
                   f"vertex."]
    return v[surf], []


def provenance_notes(has_mesh: bool, has_recon: bool) -> list[str]:
    """Note for a vessel mesh with no recorded reconstruction parameters.

    This is the regression guard for `docs/ISSUES.md` [P26]: `pipeline.batch`
    drives the GUI's mesh step at whatever the widgets default to and records
    nothing, so the five batched meshes are unattributable and the reference
    mesh's parameters survive only in one 2026-08-06 settings JSON. A mesh
    written by `pipeline.batch_mesh` carries a `vessel_recon.json`; one without
    it cannot say what produced it, and that should not pass unnoticed.
    """
    if has_mesh and not has_recon:
        return [f"WARN: {sio.filename(sio.VESSEL_MESH)} has no "
                f"{sio.filename(sio.VESSEL_RECON)} beside it — it "
                "was built by `pipeline.batch` at the widget defaults or by "
                "hand, and its parameters are unknown (docs/ISSUES.md P26). "
                "Rebuild it with `python -m pipeline.batch_mesh` to record "
                "them."]
    return []


def recon_summary(doc: dict) -> list[str]:
    """A `vessel_recon.json` as two report lines. Pure: dict in, strings out.

    Only the parameters the method actually read are shown on the mesh line —
    the ones it ignored live under `mesh_inert` and are named as such, because
    listing them as though they shaped the mesh is what made two of [P26]'s five
    "hand-tuned vs default" rows meaningless.
    """
    mesh = doc.get("mesh") or {}
    method = mesh.get("method", "?")
    knobs = " ".join(f"{k}={v}" for k, v in mesh.items() if k != "method")
    ms = doc.get("min_size") or {}
    min_size = ms.get("resolved", "?")
    mode = ms.get("mode")
    min_str = f"min_size={min_size}" + (f" ({mode})" if mode not in (None, str(min_size)) else "")
    first = f"mesh {method} {knobs} | {min_str}".rstrip()

    curv = doc.get("surface_curvature") or {}
    curv_str = (" ".join(f"{k}={v}" for k, v in curv.items()) if curv
                else "none recorded")
    counts = ""
    for name, stats in (doc.get("outputs") or {}).items():
        if "n_verts" in stats:
            counts = (f" | {stats['n_verts']:,} verts "
                      f"({stats.get('n_surface_verts', '?'):,} surface), "
                      f"{stats.get('n_faces', 0):,} faces")
            break
    inert = doc.get("mesh_inert") or {}
    inert_str = f" | {len(inert)} inert knob(s) not shown" if inert else ""
    tail = f" | {doc.get('written', '?')} {doc.get('git_sha') or ''}".rstrip()
    second = f"curv {curv_str}{counts}{inert_str}{tail}"

    out = [first, second]
    carried = doc.get("carried_over") or {}
    if carried.get("blocks"):
        out.append(f"carried over from an earlier run: "
                   f"{', '.join(carried['blocks'])}")
    if doc.get("status") and doc["status"] != "ok":
        out.append(f"status: {doc['status']}")
    return out


def fibroblast_count_notes(skipped) -> list[str]:
    """Notes for members of the fibroblast-count npz that are not scalars.

    `saved_io.load_fibroblast_count` takes the 0-d members and reports the rest
    rather than dying on them: `int(d[k])` raises TypeError on any array member,
    which would turn a slightly-different file into a hard failure of a whole
    image. The wording lives here so the GUI message and the headless one are
    the same string rather than two that drift apart (P12).
    """
    name = sio.filename(sio.FIBROBLAST_COUNT)
    return [f"WARN: {name}['{k}'] has shape {tuple(shape)}, expected a "
            f"scalar — ignored." for k, shape in skipped]


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
    mesh_name = sio.filename(sio.VESSEL_MESH)
    mp = present.get(mesh_name)
    if mp is None:
        lines.append(f"  {mesh_name:<26} —  MISSING")
    else:
        shapes = _member_shapes(mp)
        vs = shapes.get("verts")
        nverts = vs[0][0] if vs else None
        nfaces = shapes.get("faces")
        lines.append(f"  {mesh_name:<26} {nverts} verts, "
                     f"{nfaces[0][0] if nfaces else '?'} faces")

    curv_name = sio.filename(sio.VESSEL_CURVATURE)
    cp = present.get(curv_name)
    if cp is None:
        lines.append(f"  {curv_name:<26} —  MISSING")
    else:
        shapes = _member_shapes(cp)
        k1 = shapes.get("k1")
        nk = k1[0][0] if k1 else None
        lines.append(f"  {curv_name:<26} {nk} values "
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

    graph_name = sio.filename(sio.VESSEL_GRAPH)
    gp = present.get(graph_name)
    if gp is None:
        lines.append(f"  {graph_name:<26} —  MISSING")
        lines.append("  " + graph_notes(None)[0])
        if present.get(sio.VESSEL_PATH_CURVATURE_RETIRED) is not None:
            lines.append(
                f"  NOTE: {sio.VESSEL_PATH_CURVATURE_RETIRED} is present but "
                "is no longer "
                "read — the in-plane centerline was replaced by the 3D lumen "
                "graph (docs/DIRECTION.md §8B). Re-run `batch_mesh --steps "
                "mesh,graph` with a Poisson method.")
    else:
        d = _scalars(gp, ("points_zyx_um", "edge_len_um", "node_kind",
                          "node_endpoint_kind", "area_um2",
                          "section_clipped", "caliber_um"))
        npts = len(d.get("points_zyx_um", []))
        kind = np.asarray(d.get("node_kind", []))
        njunc = int((kind == "junction").sum()) if kind.size else 0
        ek = np.asarray(d.get("node_endpoint_kind", []))
        ncrop = int((ek == "crop").sum()) if ek.size else 0
        total = float(np.sum(d.get("edge_len_um", [0])))
        lines.append(f"  {graph_name:<26} {npts} points, "
                     f"{len(d.get('edge_len_um', []))} edge(s), {njunc} "
                     f"junction(s), {ncrop} cropped end(s), {total:.0f} um")
        for n in graph_notes(d):
            lines.append(f"  {n}")

    # --- where the mesh came from (P26) --------------------------------------
    tagged = sorted(n for n in present if sio.tag_of(n, sio.VESSEL_MESH))
    if tagged:
        bits = []
        for name in tagged:
            shapes = _member_shapes(present[name])
            vs = shapes.get("verts")
            tag = sio.tag_of(name, sio.VESSEL_MESH)
            bits.append(f"{tag} ({vs[0][0]:,} verts)" if vs else tag)
        lines.append(f"  {'tagged mesh variants':<26} {', '.join(bits)}")

    for name in sorted(n for n in present
                       if sio.tag_of(n, sio.VESSEL_RECON) is not None):
        try:
            doc = json.loads(present[name].read_text())
        except Exception as exc:
            lines.append(f"  {name:<26} unreadable ({exc})")
            continue
        summary = recon_summary(doc)
        lines.append(f"  {name:<26} {summary[0]}")
        lines.extend(f"  {'':<26} {extra}" for extra in summary[1:])
    for n in provenance_notes(mesh_name in present,
                              sio.filename(sio.VESSEL_RECON) in present):
        lines.append(f"  {n}")

    # --- nuclei: the P04 join -----------------------------------------------
    split_name = sio.filename(sio.NUCLEI_SPLIT)
    np_ = present.get(split_name)
    if np_ is None:
        lines.append(f"  WARN: no {split_name} — nuc_id falls back to the "
                     "unsplit DAPI components and nuc_is_fibroblast will be "
                     "blank. Re-run 'Count fibroblasts' for this image.")
    else:
        d = _scalars(np_, ("ids", "is_fib", "min_size_vox", "threshold"))
        ids = np.asarray(d.get("ids", []))
        is_fib = np.asarray(d.get("is_fib", []))
        shapes = _member_shapes(np_)
        lab = shapes.get("labels")
        lines.append(f"  {split_name:<26} {ids.size} nuclei, "
                     f"{int(is_fib.sum()) if is_fib.size else 0} VnTs+, "
                     f"labels {lab[0] if lab else '?'} "
                     f"{lab[1] if lab else ''}")
        lab_max = int(ids.max()) if ids.size else 0
        for n in nuclei_notes(ids, is_fib, lab_max):
            lines.append(f"  {n}")

    fib_name = sio.filename(sio.FIBROBLAST_COUNT)
    fp = present.get(fib_name)
    if fp is None:
        lines.append(f"  {fib_name:<26} —  MISSING "
                     f"(n_fibroblasts / n_nuclei_total will be blank)")
    else:
        d = _scalars(fp, ("fibroblast", "total"))
        lines.append(f"  {fib_name:<26} "
                     f"{int(d.get('fibroblast', 0))} of "
                     f"{int(d.get('total', 0))} nuclei")

    # --- what has already been measured -------------------------------------
    csvs = sorted(p.name for p in saved.glob("*_adhesion_measures_*.csv"))
    jsons = sorted(saved.glob("*_settings.json"))
    lines.append(f"  {'measured':<26} {len(csvs)} CSV(s), {len(jsons)} "
                 f"settings JSON"
                 + (f"; latest {csvs[-1]}" if csvs else ""))

    if sio.filename(sio.ADHESION_INSTANCES) not in present:
        lines.append(f"  WARN: no {sio.filename(sio.ADHESION_INSTANCES)} — "
                     f"this folder cannot "
                     "be measured at all; re-run the batch for this image.")
    return lines
