"""Headless vessel-mesh rebuild (no napari) — mesh, curvature, centerline, scores.

**This is an instrument, not a result.** It makes the reconstruction parameters
explicit, recorded and cheap to vary; it does not tell you which ones are right.
`docs/ISSUES.md` [P20], [P21], [P24] and [I03] are all still open after this
lands, and [P26] still stands: there is no validated parameter set, the
hand-tuned reference mesh is best-available rather than a baseline, and no
`curv_*` or `dist_mesh_um` becomes comparable between images just because it was
rebuilt here. What changes is that a mesh now says what produced it.

Why it can be headless at all: the vessel mesh is built from the **cadherin
mask**, which is saved as `cadherin_labels.npz` for every image, so rebuilding it
needs no RF predict (the expensive step) and — as it turns out — no image pixels
either:

- the isotropic **scale** comes from the OME header alone (`min(voxel_zyx_um)`,
  the target `io.to_isotropic` picks), which also gives the expected volume shape
  and so doubles as a check that the mask belongs to this image;
- the **mask** is `denoise.filter_by_size` over the saved labels, byte-identical
  to what the Cadherin tab's size slider produces;
- `cell_mesh.mesh_curvature` needs only the mesh, so a curvature-only run reads
  neither the header nor the mask.

The build itself is `vessel_mesh.build_mesh` — the *same* function the napari
tab calls, so the interactive and unattended paths cannot drift (the doctrine
`saved_check.py` was written under: one implementation, two callers).

Usage (image_analysis env, from the 3d_napari_pipeline folder):

    # what would run, resolved, without running it
    python -m pipeline.batch_mesh --dir "...\\images\\WT" --dry-run

    # rebuild one image at a recorded parameter set, as a tagged sibling
    python -m pipeline.batch_mesh --dir "...\\images\\WT" --only 04step_002 \\
        --from-settings "...\\..._20260806_145720_settings.json" --tag ref

    # P24's plateau test across every image, in one command
    python -m pipeline.batch_mesh --dir "...\\images\\WT" --steps curvature \\
        --curv-radius-um 5,10,15,20,30,40 --tag r{radius}

    # score the mesh already on disk — builds nothing, writes only the scores
    python -m pipeline.batch_mesh --dir "...\\images\\WT" --steps score

`--steps score` is the fourth step and the one that makes the other three
answerable: it reports scale stability, topology and coverage for a mesh
(`pipeline/mesh_score.py`) into the `scores` block of its `vessel_recon.json`.
It lives here rather than in a standalone script because coverage needs the
cadherin mask, and this is the only headless place that loads it — a separate
script would have to duplicate the mask loading and the OME shape check.
**Three numbers, deliberately never combined into one:** weighting them needs
weights, and an invented weight is the hidden constant this project has already
paid for twice. A selection rule belongs to whoever runs a sweep.

Exit codes: 0 all ok, 1 at least one image FAILED, 2 (with --strict) something
came out degraded.

Out of scope on purpose: promoting a tagged variant (a file copy does it),
parallelism (sweeps parallelise at the shell), and the Cadherin tab's `band_um`
filter (no batched mesh has ever used it — `batch.py::_RUN_STEPS` has no
cell-mesh step, so `band_um` has always been 0.0).
"""

from __future__ import annotations

import argparse
import faulthandler
import gc
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import tifffile

# Same reason as batch_measure: a native fault (cupy kernel, Qhull) otherwise
# returns silently to the prompt with nothing to read.
faulthandler.enable()

from .cell_mesh import mesh_curvature, suggest_alpha_um, tet_circumradii
from . import lumen_graph
from .denoise import auto_min_size_otsu, filter_by_size, local_contrast_score
from . import io as _io
from .io import _parse_ome_voxel_zyx
from .mesh_score import (
    DEFAULT_RADII_UM,
    coverage as score_coverage,
    format_scores,
    score_mesh,
    support as score_support,
)
from .saved_check import (
    _member_shapes,
    graph_notes,
    curvature_notes,
    describe_saved,
    mesh_notes,
)
from . import vessel_poisson
from .vessel_mesh import (
    ALL_PARAMS,
    METHOD_ALPHA3D,
    METHOD_POISSON,
    METHODS,
    ORIENTATIONS,
    build_mesh,
    live_params,
    mask_to_points_um,
    voxel_subsample,
)

ROOT = Path(__file__).resolve().parent.parent

STEPS = ("mesh", "curvature", "graph", "score")
# `score` is opt-in. It builds nothing, and its scale-stability sweep costs more
# than the whole rebuild does (6.6 s against 2.1 s on the reference mesh, and it
# grows with vertex count times fit radius), so defaulting it on would silently
# multiply the cost of every existing rebuild.
DEFAULT_STEPS = ("mesh", "curvature", "graph")

# Defaults are the **current widget defaults**, not the reference image's
# settings JSON. P26's 2026-08-12 correction demotes that mesh from baseline to
# best-available, so hard-coding its values here would re-enshrine a set that was
# never validated — and two of them (`close_um`, `r_smooth_bins`) were inert
# under alpha 3D anyway. Pass --from-settings to re-run it deliberately.
MESH_DEFAULTS = {
    "method": METHOD_ALPHA3D,
    "orientation": "both (OR)",
    "n_bins": 180,
    "r_smooth_bins": 3,
    "alpha2d_um": 8.0,
    # 5.0, not 1.0 — the value every mesh on disk was built at. At 1.0 the
    # cloud is ~25x denser, which is how a routine build reaches 330k points
    # and 2.2M tets: far outside `voxel_subsample`'s stated 1e4-5e4 range, and
    # the regime where the circumradius distribution degenerates (P35).
    "subsample_um": 5.0,
    # 15.0, not 10.0: every mesh on disk was built at 15.0, and `batch_mesh`
    # runs at these defaults, so a mismatch is the [P26] hazard — the shipped
    # default silently producing something other than the recorded result. It
    # is also the smallest value satisfying the >=3x subsample_um rule checked
    # in `mesh_notes` at the recorded subsample_um of 5.0.
    "alpha3d_um": 15.0,
    "fill_lumen": True,
    "close_um": 1.0,
    "min_component_um": 50.0,
    "poisson_grid_um": vessel_poisson.DEFAULT_GRID_UM,
    "poisson_smooth_um": vessel_poisson.DEFAULT_SMOOTH_UM,
    "poisson_iso_pct": vessel_poisson.DEFAULT_ISO_PCT,
    # The grid the alpha solid is rasterized onto for
    # `METHOD_ALPHA_POISSON`. Defaults to the Poisson grid so the solve
    # sees the solid at its own resolution and `occupancy_grid` reduces to
    # the identity (k = 1); a finer value hands Poisson a fractional
    # occupancy instead, which is a real choice and therefore a knob.
    "alpha_raster_um": vessel_poisson.DEFAULT_GRID_UM,
}
# --- the cadherin field (`docs/ISSUES.md` [P28]) ---------------------------
# **Not a mesh parameter and not a method.** `build_mesh` gains no argument for
# this, so `vessel_mesh`'s METHODS / LIVE_PARAMS / ALL_PARAMS are untouched and
# CLAUDE.md's four-place rule does not fire; what changes is the *input field*
# the chosen method is handed. That is also why these live in their own block
# rather than in MESH_DEFAULTS: `mesh_all` is splatted into `build_mesh`, which
# would reject them.
#
# `mask` is the default and the pixel-free rebuild contract is unchanged under
# it — no image data is decoded. `score` decodes the substrate channel and hands
# the method `denoise.local_contrast_score`, rectified at `score_k_lo`.
FIELD_MASK, FIELD_SCORE = "mask", "score"
CADHERIN_FIELDS = (FIELD_MASK, FIELD_SCORE)
# `score_k_lo` is set from a measurement, and the measurement came out worse
# than [P28] expected. `docs/scratch/p28/rays.py` (2026-08-21, five measurable
# images) found the dim gap wall at a Niblack z-score of 2.74-3.23 (p10
# 1.74-2.17) and the **no-wall controls at 3.23-3.64** — i.e. at this window the
# z-score does *not* separate real dim wall from lumen noise, so there is no
# value that keeps one and drops the other. The default is therefore the level
# that keeps ~90% of the measured gap wall (its p10), stated as what it is: a
# retention target, not a separator. What it costs in false positives is the
# thing the four-arm experiment measures.
# The default is the **lowest value the reconstruction tolerates**, not the
# retention target. `docs/scratch/p28/klo_sweep.py` measured that at k_lo = 2.0
# (the gap wall's p10) the field covers 6.2x the binary mask's voxels, the
# "wall" swallows the lumen and `lumen_core` refuses to orient — 0 vertices. At
# 2.5 it reconstructs (2.8x the voxels); at 3.0 it keeps only 1.3x, which is
# roughly the gap wall's *median* and so throws half of it away again.
P28_K_LO = 2.5
P28_K_LO_PROVENANCE = (
    "docs/scratch/p28/rays.py (2026-08-21, five measurable images: gap-wall "
    "z-score median 2.74-3.23, p10 1.74-2.17; NOT separated from the no-wall "
    "controls at 3.23-3.64) and docs/scratch/p28/klo_sweep.py (2.0 floods the "
    "lumen and yields no surface; 2.5 is the lowest that reconstructs). "
    "See docs/ISSUES.md P28")
FIELD_DEFAULTS = {"cadherin_field": FIELD_MASK, "score_window_xy": 50,
                  "score_k_lo": P28_K_LO}
CURV_DEFAULTS = {"radius_um": 5.0, "smooth_iter": 10}
# The 3D centerline graph (`lumen_graph`). Every one of these is a **ratio or a
# window, not an absolute length**, which is the point: the old centerline's
# `open_um` / `prune_um` were absolute constants and [P23] is the record of them
# amputating exactly the tapering, high-shear regions the program cares about.
GRAPH_DEFAULTS = {
    "spur_ratio": lumen_graph.DEFAULT_SPUR_RATIO,
    "crop_ratio": lumen_graph.DEFAULT_CROP_RATIO,
    "crop_margin_um": lumen_graph.DEFAULT_CROP_MARGIN_UM,
    "kappa_window_um": lumen_graph.DEFAULT_KAPPA_WINDOW_UM,
    "section_every_um": lumen_graph.DEFAULT_SECTION_EVERY_UM,
    "section_span": lumen_graph.DEFAULT_SECTION_SPAN,
}
# The radii P24's plateau table was measured at, so a `--steps score` run is
# directly comparable to the evidence already on record.
SCORE_DEFAULTS = {"radii_um": ",".join(f"{r:g}" for r in DEFAULT_RADII_UM),
                  "max_verts": None, "max_query": None}
# The despeckle cut the Cadherin tab applies (`gui/cadherin.py:42,191`):
# max(auto_min_size_otsu(sizes), 12). P19 measured the auto value at 3 on the
# reference image, so the constant is what is always in force — `auto-raw` runs
# without the floor, which is the experiment P19 asks for and nobody has run.
MIN_SIZE_FLOOR = 12


def _p(msg: str) -> None:
    print(f"      . {msg}", flush=True)


def _git_sha() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=str(ROOT), capture_output=True, text=True,
                             timeout=10)
        return out.stdout.strip() or None
    except Exception:
        return None


def _jsonable(obj):
    """Numpy scalars/arrays -> plain Python, recursively, for the recon JSON.

    ``json.dumps(default=str)`` would otherwise turn a `numpy.float64` into the
    *string* ``'0.123'``, which reads fine and compares wrong — a recon file is
    meant to be diffable against the next run. ``notes`` is dropped: it is
    already carried in the report's own notes list, and two copies drift.
    """
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items() if k != "notes"}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def _mtime(p: Path) -> str:
    return datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")


# --------------------------------------------------------------------------- #
# Scale from the OME header — no pixels decoded
# --------------------------------------------------------------------------- #

def _nanmed(a):
    """Median ignoring NaN, or None when everything is NaN — the sections are
    cut every few µm and are NaN in between by design, so a bare median is
    always NaN."""
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else None


def _header_scale(tif: Path) -> tuple[tuple, tuple, tuple]:
    """``(voxel_raw_zyx, scale_iso_zyx, expected_iso_shape)`` from the header.

    `io.to_isotropic` resamples to ``target = min(voxel_zyx_um)`` and
    `io.resample_volume` delegates to `ndimage.zoom`, whose output shape is
    ``round(n · src/target)`` per axis. So the header alone gives both the scale
    every saved artifact is on and the shape the cadherin labels must have —
    which is why the shape check below is free and catches a mask from the wrong
    image (or from a different resampling) before anything is built.
    """
    with tifffile.TiffFile(str(tif)) as tf:
        series = tf.series[0]
        if series.axes != "ZCYX":
            raise ValueError(f"{tif.name}: expected ZCYX axes, got {series.axes}")
        raw_shape = dict(zip(series.axes, series.shape))
        xml = tf.ome_metadata or ""
    voxel = _parse_ome_voxel_zyx(xml)
    if any(not np.isfinite(v) or v <= 0 for v in voxel):
        raise ValueError(f"{tif.name}: invalid OME voxel size {voxel} "
                         f"(PhysicalSizeZ/Y/X)")
    target = float(min(voxel))
    src = (raw_shape["Z"], raw_shape["Y"], raw_shape["X"])
    expected = tuple(int(round(n * v / target)) for n, v in zip(src, voxel))
    return tuple(float(v) for v in voxel), (target, target, target), expected


def _substrate_volume(tif: Path, expected, target_um: float) -> np.ndarray:
    """The substrate channel on the isotropic grid the saved mask is on.

    `io.to_isotropic` would resample all five channels; `resample_volume` does
    the one, with the same target and the same order, so the result is what
    `to_isotropic` would have produced for this channel. The shape is checked
    against the OME header's expectation for the same reason the mask is —
    a field that is not on the mask's grid is not comparable to it.

    Only `--cadherin-field score` calls this. Under the default the rebuild
    still decodes no image pixels.
    """
    params = _io.load_params(ROOT / "3d_images_folder" / "params.json")
    stack = _io.load_stack(tif, params)
    ci = stack.channels.index(stack.roles["substrate"])
    vol = _io.resample_volume(stack.data[:, ci], stack.voxel_zyx_um,
                              float(target_um), order=1)
    if tuple(vol.shape) != tuple(expected):
        raise ValueError(f"{tif.name}: the substrate channel resampled to "
                         f"{tuple(vol.shape)}, but the OME header implies "
                         f"{tuple(expected)}")
    return np.ascontiguousarray(vol)


# --------------------------------------------------------------------------- #
# Parameter resolution: built-in defaults < --from-settings < explicit flags
# --------------------------------------------------------------------------- #

def _seed_from(path: Path) -> tuple[dict, dict, dict, str | None]:
    """Read mesh / curvature / path / min_size seeds out of a JSON.

    Accepts either a measure settings JSON (``settings.vessel_mesh`` etc., what
    the GUI Measure tab writes) or a ``vessel_recon.json`` written by this tool.
    Unknown keys are ignored, so a partial file seeds what it has.
    """
    with open(path) as fh:
        d = json.load(fh)
    settings = d.get("settings", d)
    mesh = dict(settings.get("vessel_mesh") or d.get("mesh") or {})
    # A recon JSON splits the knobs the method actually read from the ones it
    # ignored; merge them back so a recon file round-trips even if the method
    # is then changed on the command line.
    mesh.update(d.get("mesh_inert") or {})
    curv = dict(settings.get("surface_curvature") or d.get("surface_curvature")
                or {})
    graphp = dict(settings.get("centerline_graph")
                  or d.get("centerline_graph") or {})
    ms = d.get("min_size") or {}
    min_size = str(ms.get("mode")) if isinstance(ms, dict) and ms.get("mode") \
        else None
    return mesh, curv, graphp, min_size


def _pick(explicit, seed, default, key):
    """explicit flag > --from-settings > built-in default, plus its source."""
    if explicit is not None:
        return explicit, "flag"
    if key in seed and seed[key] is not None:
        return seed[key], "from-settings"
    return default, "default"


def _floats(text: str, what: str) -> list[float]:
    try:
        return [float(x) for x in str(text).split(",") if x.strip() != ""]
    except ValueError:
        raise SystemExit(f"--{what}: expected a number or comma-separated "
                         f"numbers, got {text!r}")


def _ints(text: str, what: str) -> list[int]:
    return [int(round(v)) for v in _floats(text, what)]


def resolve_params(args) -> dict:
    """Everything the run needs, with where each value came from."""
    seed_mesh, seed_curv, seed_graph, seed_min = {}, {}, {}, None
    if args.from_settings:
        src = Path(args.from_settings)
        if not src.exists():
            raise SystemExit(f"--from-settings: {src} does not exist.")
        seed_mesh, seed_curv, seed_graph, seed_min = _seed_from(src)

    mesh, source = {}, {}
    explicit = {
        "method": args.method, "orientation": args.orientation,
        "n_bins": args.n_bins, "r_smooth_bins": args.r_smooth_bins,
        "alpha2d_um": args.alpha2d_um, "subsample_um": args.subsample_um,
        "alpha3d_um": args.alpha3d_um, "fill_lumen": args.fill_lumen,
        "close_um": args.close_um, "min_component_um": args.min_component_um,
        "poisson_grid_um": args.poisson_grid_um,
        "poisson_smooth_um": args.poisson_smooth_um,
        "poisson_iso_pct": args.poisson_iso_pct,
        "alpha_raster_um": args.alpha_raster_um,
    }
    for key in ALL_PARAMS:
        mesh[key], source[f"mesh.{key}"] = _pick(
            explicit[key], seed_mesh, MESH_DEFAULTS[key], key)

    curv = {}
    curv["radius_um"], source["curv.radius_um"] = _pick(
        args.curv_radius_um, seed_curv, CURV_DEFAULTS["radius_um"], "radius_um")
    curv["smooth_iter"], source["curv.smooth_iter"] = _pick(
        args.curv_smooth_iter, seed_curv, CURV_DEFAULTS["smooth_iter"],
        "smooth_iter")
    radii = _floats(curv["radius_um"], "curv-radius-um") \
        if isinstance(curv["radius_um"], str) else [float(curv["radius_um"])]
    smooths = _ints(curv["smooth_iter"], "curv-smooth-iter") \
        if isinstance(curv["smooth_iter"], str) else [int(curv["smooth_iter"])]
    combos = [(r, s) for r in radii for s in smooths]

    score = {}
    score["radii_um"], source["score.radii_um"] = _pick(
        args.score_radii, {}, SCORE_DEFAULTS["radii_um"], "radii_um")
    score["radii_um"] = _floats(score["radii_um"], "score-radii") \
        if isinstance(score["radii_um"], str) else [float(score["radii_um"])]
    score["max_verts"], source["score.max_verts"] = _pick(
        args.score_max_verts, {}, SCORE_DEFAULTS["max_verts"], "max_verts")
    score["max_query"], source["score.max_query"] = _pick(
        args.score_max_query, {}, SCORE_DEFAULTS["max_query"], "max_query")
    # The scores fit curvature with the same smoothing the curvature step uses,
    # so a score and a `vessel_curvature.npz` from the same run describe the same
    # surface. One knob, not two that can disagree.
    score["smooth_iter"] = smooths[0]

    field = {}
    explicit_field = {"cadherin_field": args.cadherin_field,
                      "score_window_xy": args.score_window_xy,
                      "score_k_lo": args.score_k_lo}
    for key, default in FIELD_DEFAULTS.items():
        field[key], source[f"field.{key}"] = _pick(
            explicit_field[key], {}, default, key)

    graph = {}
    explicit_graph = {
        "spur_ratio": args.spur_ratio, "crop_ratio": args.crop_ratio,
        "crop_margin_um": args.crop_margin_um,
        "kappa_window_um": args.kappa_window_um,
        "section_every_um": args.section_every_um,
        "section_span": args.section_span}
    for key, default in GRAPH_DEFAULTS.items():
        graph[key], source[f"graph.{key}"] = _pick(
            explicit_graph[key], seed_graph, default, key)

    if args.min_size is not None:
        min_size, source["min_size"] = str(args.min_size), "flag"
    elif seed_min:
        min_size, source["min_size"] = str(seed_min), "from-settings"
    else:
        min_size, source["min_size"] = str(MIN_SIZE_FLOOR), "default"
    if min_size not in ("auto", "auto-raw"):
        try:
            int(min_size)
        except ValueError:
            raise SystemExit(f"--min-size: expected an integer, 'auto' or "
                             f"'auto-raw', got {min_size!r}")

    # Parameters this method never reads. Recording them as though they shaped
    # the mesh is the P26 failure: the reference JSON lists close_um=5.0 and
    # r_smooth_bins=12 for an alpha-3D mesh that read neither, so two of that
    # entry's five "hand-tuned vs default" rows are meaningless.
    live = set(live_params(mesh["method"]))
    inert = {k: v for k, v in mesh.items()
             if k != "method" and k not in live}
    live_mesh = {k: v for k, v in mesh.items() if k == "method" or k in live}

    return {"mesh": live_mesh, "mesh_inert": inert, "mesh_all": mesh,
            "curv_combos": combos, "graph": graph, "score": score,
            "field": field, "min_size": min_size, "source": source}


def parameter_warnings(params: dict) -> list[str]:
    """Warnings about the resolved set itself, before any image is touched."""
    notes = []
    mesh = params["mesh_all"]
    method = mesh["method"]
    src = params["source"]
    # Only knobs that were actually *moved* are worth naming. A --from-settings
    # file supplies all ten, most of them at the default, and warning about every
    # inert one would bury the two that matter (P26's `close_um` and
    # `r_smooth_bins`, hand-tuned on a mesh that read neither).
    set_inert = sorted(
        k for k, v in params["mesh_inert"].items()
        if src.get(f"mesh.{k}") == "flag" or v != MESH_DEFAULTS[k])
    if set_inert:
        vals = ", ".join(f"{k}={mesh[k]}" for k in set_inert)
        notes.append(f"WARN: {vals} — ignored by this method ({method}); "
                     f"recorded under 'mesh_inert', not as a mesh parameter.")
    fld = params["field"]
    if fld["cadherin_field"] == FIELD_SCORE and method != METHOD_POISSON:
        notes.append(
            f"NOTE: --cadherin-field score with method '{method}', which "
            f"interpolates and so takes `score > 0` — i.e. a Niblack threshold "
            f"at k = score_k_lo = {fld['score_k_lo']} instead of the saved "
            f"mask. That is the intended fourth arm of docs/ISSUES.md P28's "
            f"experiment, not a mistake, but it is a *looser mask*, not a "
            f"continuous field: only Poisson reads the values.")
    if fld["cadherin_field"] == FIELD_MASK and (
            src.get("field.score_window_xy") == "flag"
            or src.get("field.score_k_lo") == "flag"):
        notes.append(
            "WARN: --score-window-xy / --score-k-lo were given but "
            "--cadherin-field is 'mask', so neither is read.")
    if method == METHOD_ALPHA3D:
        a3, sub = float(mesh["alpha3d_um"]), float(mesh["subsample_um"])
        if sub > 0 and a3 < 3.0 * sub:
            notes.append(
                f"WARN: alpha3d_um={a3} is below 3x subsample_um={sub}. Alpha "
                f"has to exceed the point spacing simply to connect neighbours "
                f"(docs/ISSUES.md P20), and this is the commonest route to "
                f"'N verts, 0 faces'.")
    return notes


def format_params(params: dict) -> list[str]:
    """The resolved parameters as report lines, with each value's source."""
    src = params["source"]
    out = ["mesh:"]
    for k, v in params["mesh"].items():
        out.append(f"    {k:<18} {v}   [{src.get(f'mesh.{k}', 'default')}]")
    for k, v in params["mesh_inert"].items():
        out.append(f"    {k:<18} {v}   [{src.get(f'mesh.{k}', 'default')}, "
                   f"INERT for this method]")
    out.append("surface curvature:")
    combos = params["curv_combos"]
    radii = ", ".join(str(r) for r in sorted({r for r, _ in combos}))
    smooths = ", ".join(str(s) for s in sorted({s for _, s in combos}))
    out.append(f"    {'radius_um':<18} {radii}   "
               f"[{src.get('curv.radius_um', 'default')}]")
    out.append(f"    {'smooth_iter':<18} {smooths}   "
               f"[{src.get('curv.smooth_iter', 'default')}]")
    if len(combos) > 1:
        out.append(f"    -> {len(combos)} curvature run(s) per image")
    out.append("cadherin field:")
    fld = params["field"]
    out.append(f"    {'cadherin_field':<18} {fld['cadherin_field']}   "
               f"[{src.get('field.cadherin_field', 'default')}]")
    if fld["cadherin_field"] == FIELD_SCORE:
        for k in ("score_window_xy", "score_k_lo"):
            out.append(f"    {k:<18} {fld[k]}   "
                       f"[{src.get(f'field.{k}', 'default')}]")
    else:
        out.append("    (the pixel-free rebuild: no image data is decoded)")
    out.append("centerline graph:")
    for k, v in params["graph"].items():
        out.append(f"    {k:<18} {v}   [{src.get(f'graph.{k}', 'default')}]")
    out.append("mesh scores:")
    sc = params["score"]
    out.append(f"    {'radii_um':<18} "
               f"{', '.join(f'{r:g}' for r in sc['radii_um'])}   "
               f"[{src.get('score.radii_um', 'default')}]")
    for k in ("max_verts", "max_query"):
        out.append(f"    {k:<18} {sc[k] if sc[k] is not None else 'all'}   "
                   f"[{src.get(f'score.{k}', 'default')}]")
    out.append(f"    {'smooth_iter':<18} {sc['smooth_iter']}   "
               f"[follows curv.smooth_iter]")
    out.append(f"min_size             {params['min_size']}   "
               f"[{src.get('min_size', 'default')}]")
    return out


# --------------------------------------------------------------------------- #
# The mask
# --------------------------------------------------------------------------- #

def resolve_min_size(mode: str, sizes: np.ndarray) -> tuple[int, dict]:
    """``(min_size, how)``. `how` is recorded verbatim in the recon JSON."""
    auto = int(auto_min_size_otsu(sizes))
    if mode == "auto":
        val = max(auto, MIN_SIZE_FLOOR)
        how = {"mode": "auto", "auto_otsu": auto, "floor": MIN_SIZE_FLOOR,
               "note": "max(otsu, 12) — what gui/cadherin.py applies"}
    elif mode == "auto-raw":
        val = auto
        how = {"mode": "auto-raw", "auto_otsu": auto, "floor": None,
               "note": "the Otsu cut with no floor (docs/ISSUES.md P19)"}
    else:
        val = int(mode)
        how = {"mode": mode, "auto_otsu": auto, "floor": None,
               "note": "given explicitly"}
    how["resolved"] = int(val)
    return int(val), how


# --------------------------------------------------------------------------- #
# Per-image run
# --------------------------------------------------------------------------- #

def _out(saved: Path, base: str, tag: str | None, ext: str) -> Path:
    """Tagged siblings live **beside** the canonical files, not in a subfolder:
    `describe_saved` iterates files and `batch_measure` builds literal
    filenames, so anything one directory down is invisible to both.
    """
    return saved / (f"{base}__{tag}{ext}" if tag else f"{base}{ext}")


def _existing_mesh_note(path: Path) -> str:
    """What is on disk right now, for the 'I re-ran it and nothing changed' case.

    A failure that says only "failed" leaves the previous mesh looking like the
    new one — the exact trap the GUI's debounce race sets, where the build ran
    on an empty mask, `_save_mesh` wrote nothing, and the stale file survived.
    """
    if not path.exists():
        return f"no {path.name} was on disk, and none was written"
    shapes = _member_shapes(path)
    vs, fs = shapes.get("verts"), shapes.get("faces")
    if vs and fs:
        return (f"left the existing {path.name} in place ({vs[0][0]:,} verts, "
                f"{fs[0][0]:,} faces, written {_mtime(path)})")
    return f"left the existing {path.name} in place (written {_mtime(path)})"


def run_one(tif: Path, saved: Path, params: dict, steps: tuple,
            tag_template: str | None, dry_run: bool = False,
            suggest_alpha: bool = False) -> tuple[str, list[str]]:
    """Rebuild one image. Returns ``(status, notes)``.

    Everything is computed in memory and written only once every requested step
    has succeeded, so a run can never leave a mesh from one parameter set beside
    a curvature file from another — which is precisely [P12], where
    `measure_core` then blanks all four `curv_*` with nothing in the log.
    """
    notes: list[str] = []
    t0 = time.perf_counter()
    timings: dict[str, float] = {}
    recon: dict = {}
    # `score` needs the mask too: coverage is measured from every retained
    # cadherin voxel to the surface, and this is the only headless place the
    # mask is loaded. A standalone scoring script would have to duplicate the
    # mask-loading and the OME shape check, which is the second copy the
    # project rule forbids.
    need_mask = ("mesh" in steps) or ("graph" in steps) or ("score" in steps) \
        or suggest_alpha
    combos = params["curv_combos"]
    single_tag = (tag_template.format(radius=combos[0][0], smooth=combos[0][1])
                  if tag_template else None)
    mesh_out = _out(saved, "vessel_mesh", single_tag, ".npz")

    # ---- scale + shape, from the OME header (~0.3 s, no pixels) -------------
    scale = None
    if need_mask:
        t = time.perf_counter()
        voxel_raw, scale, expected = _header_scale(tif)
        timings["header_s"] = time.perf_counter() - t
        _p(f"OME voxel {voxel_raw} um -> isotropic {scale[0]:.4f} um; "
           f"expected volume {expected}")

    # ---- the cadherin mask --------------------------------------------------
    mask = None
    labels_meta = None
    if need_mask:
        cadh = saved / "cadherin_labels.npz"
        if not cadh.exists():
            return "failed", [f"no cadherin_labels.npz in {saved.name} — this "
                              f"image has not been through the Cadherin step; "
                              f"{_existing_mesh_note(mesh_out)}."]
        t = time.perf_counter()
        with np.load(cadh) as d:
            labels = d["labels"]
            sizes = d["sizes"]
        labels_meta = {"file": cadh.name, "mtime": _mtime(cadh),
                       "bytes": cadh.stat().st_size,
                       "shape": [int(n) for n in labels.shape],
                       "expected_shape": [int(n) for n in expected]}
        if tuple(labels.shape) != tuple(expected):
            return "failed", [
                f"cadherin_labels.npz is {tuple(labels.shape)} but "
                f"{tif.name}'s OME header implies {tuple(expected)} after "
                f"to_isotropic — the mask is not from this image (or was saved "
                f"at a different resampling). {_existing_mesh_note(mesh_out)}."]
        min_size, how = resolve_min_size(params["min_size"], sizes)
        mask = filter_by_size(labels, sizes, min_size)
        del labels
        n_kept = int((sizes >= min_size).sum())
        n_vox = int(mask.sum())
        how.update({"n_components": int(sizes.size), "n_kept": n_kept,
                    "mask_voxels": n_vox})
        timings["mask_s"] = time.perf_counter() - t
        _p(f"min_size={min_size} ({how['mode']}; otsu={how['auto_otsu']}) -> "
           f"{n_kept:,}/{sizes.size:,} components, {n_vox:,} mask voxels")
        del sizes
        recon["min_size"] = how
        recon["cadherin_labels"] = labels_meta
        recon["scale_um_zyx"] = list(scale)
        recon["scale_source"] = (
            f"OME PhysicalSizeZ/Y/X of {tif.name} = {list(voxel_raw)}; "
            f"isotropic target = min = {scale[0]}")
        if n_vox == 0:
            return "failed", [f"the cadherin mask is empty at min_size="
                              f"{min_size} — nothing to mesh. "
                              f"{_existing_mesh_note(mesh_out)}."]

    if suggest_alpha:
        sub = float(params["mesh_all"]["subsample_um"])
        pts = voxel_subsample(mask_to_points_um(mask > 0, scale), sub)
        R = tet_circumradii(pts)
        if R.size == 0:
            _p(f"only {len(pts)} points — need >=4 to tetrahedralize.")
            return "ok", []
        fin = R[np.isfinite(R)]
        print(f"    suggest-alpha: {len(pts):,} points at subsample_um={sub}, "
              f"{R.size:,} tets, R {fin.min():.2f}-{fin.max():.2f} um, "
              f"suggested alpha3d_um = {suggest_alpha_um(R):.2f}", flush=True)
        return "ok", []

    if dry_run:
        return "ok", []

    # ---- the cadherin field the method is handed ---------------------------
    # `mask` stays the size-filtered binary mask throughout: it is the pinned
    # yardstick every recorded number in ISSUES.md was measured against, and the
    # centerline and the scored point set both read it. `mesh_input` is what
    # build_mesh gets, which is the only thing --cadherin-field changes.
    mesh_input = mask
    field_pts = None
    fld = params["field"]
    if "mesh" in steps and need_mask and fld["cadherin_field"] == FIELD_SCORE:
        t = time.perf_counter()
        vol = _substrate_volume(tif, expected, float(scale[0]))
        mesh_input = local_contrast_score(vol, int(fld["score_window_xy"]),
                                          scale, float(fld["score_k_lo"]))
        del vol
        gc.collect()
        n_field = int((mesh_input > 0).sum())
        timings["field_s"] = time.perf_counter() - t
        recon["cadherin_field"] = {
            "field": FIELD_SCORE,
            "window_xy": int(fld["score_window_xy"]),
            "k_lo": float(fld["score_k_lo"]),
            "k_lo_provenance": (P28_K_LO_PROVENANCE
                                if params["source"].get("field.score_k_lo")
                                != "flag"
                                else "given explicitly with --score-k-lo"),
            "source": f"substrate channel of {tif.name}, resampled to "
                      f"{scale[0]:.4f} um",
            "n_voxels_nonzero": n_field,
            "n_voxels_binary_mask": int(mask.sum()),
            "mean_score_where_nonzero": (
                float(mesh_input[mesh_input > 0].mean()) if n_field else 0.0),
            "note": ("the global intensity floor is NOT applied to the score "
                     "(P28: the floor, not the local threshold, is what "
                     "deletes the dim gap wall)"),
        }
        field_pts = mask_to_points_um(mesh_input > 0, scale)
        _p(f"cadherin field: score window_xy="
           f"{fld['score_window_xy']} k_lo={fld['score_k_lo']} -> "
           f"{n_field:,} nonzero voxels against {int(mask.sum()):,} in the "
           f"binary mask  [{timings['field_s']:.1f} s]")
        notes.append(
            f"NOTE: built on the P28 local-contrast field, not the binary "
            f"mask ({n_field:,} nonzero voxels vs {int(mask.sum()):,}). The "
            f"scores below are still measured against the **binary** mask, so "
            f"they stay comparable to every recorded number.")
    elif need_mask:
        recon["cadherin_field"] = {
            "field": FIELD_MASK,
            "note": "the size-filtered binary mask; no image pixels decoded"}

    # ---- mesh ---------------------------------------------------------------
    verts = faces = None
    # Per-vertex normals, when the method produces them (Poisson does). They are
    # not decoration: `mesh_curvature` reads the global orientation off the mesh
    # otherwise, and cannot on a surface that is open by construction, so the
    # sign of every kappa2 would be a coin flip (P01(B)).
    normals = None
    mesh_extras: dict = {}
    if "mesh" in steps:
        t = time.perf_counter()
        extras = mesh_extras
        kind, verts, faces, _vals = build_mesh(mesh_input, scale,
                                               extras=extras,
                                               **params["mesh_all"])
        normals = extras.get("normals")
        timings["mesh_s"] = time.perf_counter() - t
        if extras.get("info"):
            recon["mesh_diagnostics"] = _jsonable(extras["info"])
            for n in extras["info"].get("notes", []):
                notes.append(n)
        if kind == "empty" or verts is None or len(verts) == 0:
            return "failed", [f"the mesh has no vertices. "
                              f"{_existing_mesh_note(mesh_out)}."]
        if faces is None or len(faces) == 0:
            m = params["mesh_all"]["method"]
            if m == METHOD_ALPHA3D:
                hint = "raise alpha3d_um (P20: it must exceed ~3x subsample_um)"
            elif m == METHOD_POISSON:
                hint = ("lower poisson_iso_pct or raise poisson_smooth_um, or "
                        "lower min_component_um")
            else:
                hint = ("loosen the fill / orientation, or lower "
                        "min_component_um")
            return "failed", [f"{len(verts):,} vertices but 0 faces — {hint}. "
                              f"{_existing_mesh_note(mesh_out)}."]
        n_surf = int(np.unique(faces).size)
        _p(f"mesh: {len(verts):,} verts ({n_surf:,} on the surface), "
           f"{len(faces):,} faces  [{timings['mesh_s']:.1f} s]")
        if n_surf < 6:
            return "failed", [f"only {n_surf} surface vertices — "
                              f"mesh_curvature needs >=6 and would return "
                              f"all-NaN. {_existing_mesh_note(mesh_out)}."]
        for n in mesh_notes(verts):
            notes.append(n)
    elif "curvature" in steps or "score" in steps:
        # Curvature needs a mesh but not the header or the mask; scoring needs
        # both but still builds nothing. Prefer the tagged sibling when one
        # exists, else the canonical file — which is what makes
        # `--steps curvature --tag r{radius}` a sweep over the mesh that is
        # actually in place, and `--steps score --tag defaults` a score of the
        # variant on disk rather than of a fresh rebuild.
        src = mesh_out if mesh_out.exists() else saved / "vessel_mesh.npz"
        what = "curvature" if "curvature" in steps else "score"
        if not src.exists():
            return "failed", [f"no {src.name} to compute {what} from — run "
                              f"--steps mesh first."]
        with np.load(src) as d:
            verts = d["verts"]
            faces = d["faces"] if "faces" in d else np.zeros((0, 3), np.int64)
            # An additive key, so a mesh written before this method existed
            # still loads (the same pattern `z_um` uses in the path-curvature
            # npz). Absent means "derive the normals", which is right for every
            # closed-ish method.
            normals = d["normals"] if "normals" in d.files else None
        recon["mesh_input"] = src.name
        _p(f"{what} source: {src.name} ({len(verts):,} verts, "
           f"{len(faces):,} faces)")

    # ---- surface curvature (possibly a sweep) -------------------------------
    curv_runs: list[tuple[float, int, dict, str | None]] = []
    if "curvature" in steps:
        for radius, smooth in combos:
            t = time.perf_counter()
            curv = mesh_curvature(verts, faces, float(radius),
                                  smooth_iter=int(smooth), normals=normals)
            timings.setdefault("curvature_s", 0.0)
            timings["curvature_s"] += time.perf_counter() - t
            k1 = np.asarray(curv["k1"])
            k2 = np.asarray(curv["k2"])
            n_fin = int(np.isfinite(k1).sum())
            k2_med = (float(np.median(k2[np.isfinite(k2)]))
                      if np.isfinite(k2).any() else None)
            tg = (tag_template.format(radius=radius, smooth=smooth)
                  if tag_template else None)
            curv_runs.append((float(radius), int(smooth), curv, tg))
            _p(f"curvature radius_um={radius} smooth_iter={smooth}: "
               f"{n_fin:,} finite k1"
               + (f", k2 median {k2_med:.4f} 1/um" if k2_med is not None
                  else ", k2 all non-finite")
               + (f"  -> tag '{tg}'" if tg else ""))
            for n in curvature_notes(curv, len(verts)):
                notes.append(n)

    # ---- the 3D centerline graph -------------------------------------------
    # Built on the Poisson lumen, not on a max projection of the mask — see
    # `lumen_graph`. It therefore needs a Poisson-family mesh: the lumen field
    # comes back through `build_mesh`'s `extras`, and the interpolating methods
    # never solve for one.
    cl = None
    if "graph" in steps:
        lum = (mesh_extras or {}).get("lumen") or {}
        if lum.get("chi") is None:
            notes.append(
                f"WARN: --steps graph needs a reconstruction that solves for a "
                f"lumen field, and '{params['mesh']['method']}' does not. Use "
                f"'{METHOD_POISSON}' or '{METHOD_ALPHA_POISSON}'; no "
                f"centerline was written.")
        else:
            t = time.perf_counter()
            cl = lumen_graph.build_graph(
                lum["chi"], lum["iso_level"], lum["h_um"], lum["origin_um"],
                core=lum.get("core"),
                domain_shape_um=np.asarray(mask.shape) * float(scale[0]),
                **params["graph"])
            timings["graph_s"] = time.perf_counter() - t
            _p(f"centerline graph: {len(cl['points_zyx_um']):,} points, "
               f"{len(cl['edge_len_um'])} edge(s), {cl['n_junctions']} "
               f"junction(s), {cl['n_endpoints']} endpoint(s) "
               f"({cl['n_crop_endpoints']} cropped), "
               f"{cl['total_length_um']:.0f} um  [{timings['graph_s']:.1f} s]")
            for n in cl["notes"]:
                notes.append(n)
            for n in graph_notes(cl):
                notes.append(n)

    # ---- mesh quality scores ------------------------------------------------
    # Three numbers, never one — see `mesh_score`. This step writes no npz: the
    # scores describe a mesh that already exists, so they belong in that mesh's
    # provenance file and nowhere else.
    scores = None
    findings: list[str] = []
    if "score" in steps:
        t = time.perf_counter()
        sc = params["score"]
        # **The yardstick is pinned to the binary mask, in every arm.** Under
        # --cadherin-field score the *mesh* is built from the field, but the
        # point set it is scored against stays the size-filtered binary mask —
        # otherwise the four-vector would stop being comparable to every number
        # recorded in docs/ISSUES.md, and the P28 experiment would be measuring
        # two changes at once. The field-derived coverage is computed below and
        # labelled as the secondary number it is.
        scores = score_mesh(verts, faces, mask_to_points_um(mask > 0, scale),
                            normals=normals,
                            radii_um=sc["radii_um"],
                            smooth_iter=int(sc["smooth_iter"]),
                            max_verts=sc["max_verts"],
                            max_query=sc["max_query"])
        timings["score_s"] = time.perf_counter() - t
        for line in format_scores(scores):
            _p(line)
        # A score's warnings are findings about the **mesh**, not degradations
        # of this run, so they are reported and recorded but kept out of
        # `status`. The mesh being a genus-125 sponge is the scorer working, not
        # the run failing — and since every candidate in a sweep will produce
        # findings like these, folding them into `degraded` would make --strict
        # fire on all of them and so tell you nothing.
        findings.extend(n for n in scores["notes"] if n.startswith("WARN:"))
        if field_pts is not None and len(field_pts):
            # Secondary, and never mixed into the vector above: coverage and
            # support of the same surface against the *field's* points rather
            # than the binary mask's. It answers "is the surface where the
            # continuous field says wall is", which is a different question
            # from the one every recorded number answers.
            sc2 = params["score"]
            scores["coverage_vs_field"] = score_coverage(
                verts, faces, field_pts, max_query=sc2["max_query"])
            scores["support_vs_field"] = score_support(verts, faces, field_pts)
            scores["coverage_vs_field"]["note"] = (
                "SECONDARY: measured against the P28 local-contrast field's "
                "nonzero voxels, not the binary mask. Not comparable to the "
                "coverage numbers recorded in docs/ISSUES.md.")
            scores["support_vs_field"]["note"] = \
                scores["coverage_vs_field"]["note"]
            cvf = (scores["coverage_vs_field"].get("covered_fraction_um")
                   or {})
            cvd = (scores["coverage_vs_field"].get("distance_um") or {})
            svd = (scores["support_vs_field"].get("distance_um") or {})
            frac = "  ".join(f"<={t} um {100 * (v or 0):.1f}%"
                             for t, v in cvf.items())
            _p(f"    [SECONDARY, vs the field's {len(field_pts):,} points] "
               f"coverage median {cvd.get('median') or float('nan'):.2f} um  "
               f"{frac}")
            _p(f"    [SECONDARY, vs the field] support median "
               f"{svd.get('median') or float('nan'):.2f} um")
        for n in scores["notes"]:
            if not n.startswith("WARN:"):
                _p(n)

    # ---- write, now that every requested step has produced something --------
    # Outputs are grouped by the tag they were written under, so a curvature
    # sweep gets one self-describing recon JSON per tag rather than one that
    # claims the first radius and silently owns the others' files.
    groups: dict[str | None, dict] = {}

    def _group(tg):
        return groups.setdefault(tg, {"outputs": {}, "curv": None,
                                      "scores": None})

    if "score" in steps:
        # Scoring produces no file, so it would otherwise leave `groups` empty
        # and write no recon JSON at all on a `--steps score` run.
        _group(single_tag)["scores"] = scores

    if "mesh" in steps:
        arrays = {"verts": np.asarray(verts), "faces": np.asarray(faces)}
        if normals is not None:
            arrays["normals"] = np.asarray(normals)
        np.savez(mesh_out, **arrays)
        _group(single_tag)["outputs"][mesh_out.name] = {
            "n_verts": int(len(verts)),
            # Both figures, always. P26's own correction turned on this
            # distinction: the total-vertex count overstated the curvature
            # mismatch as 23-30x when the surface figure is 5x, because
            # curvature is fitted on `np.unique(faces)` and nothing else.
            "n_surface_verts": int(np.unique(faces).size),
            "n_faces": int(len(faces))}
        _p(f"wrote {mesh_out.name}")
        recon["mesh_input"] = mesh_out.name

    for radius, smooth, curv, tg in curv_runs:
        cp = _out(saved, "vessel_curvature", tg, ".npz")
        np.savez_compressed(cp, **{k: np.asarray(curv[k])
                                   for k in ("mean", "k1", "k2", "gaussian")})
        k1 = np.asarray(curv["k1"])
        k2 = np.asarray(curv["k2"])
        g = _group(tg)
        g["outputs"][cp.name] = {
            "radius_um": radius, "smooth_iter": smooth,
            "n_values": int(k1.size),
            "n_finite_k1": int(np.isfinite(k1).sum()),
            "k2_median": (float(np.median(k2[np.isfinite(k2)]))
                          if np.isfinite(k2).any() else None)}
        g["curv"] = {"radius_um": radius, "smooth_iter": smooth}
        _p(f"wrote {cp.name}")

    if cl is not None:
        pp = _out(saved, "vessel_graph", single_tag, ".npz")
        # Only the arrays; the notes and the nested info dicts go into
        # vessel_recon.json, which is where provenance already lives.
        np.savez_compressed(pp, **{k: v for k, v in cl.items()
                                   if isinstance(v, np.ndarray)})
        _group(single_tag)["outputs"][pp.name] = {
            "n_points": int(len(cl["points_zyx_um"])),
            "n_edges": int(len(cl["edge_len_um"])),
            "n_junctions": int(cl["n_junctions"]),
            "n_endpoints": int(cl["n_endpoints"]),
            "n_crop_endpoints": int(cl["n_crop_endpoints"]),
            "n_terminus_endpoints": int(cl["n_terminus_endpoints"]),
            "total_length_um": float(cl["total_length_um"]),
            "lumen_volume_um3": float(cl["lumen_volume_um3"]),
            "caliber_um_median": _nanmed(cl["caliber_um"]),
            "section_area_um2_median": _nanmed(cl["area_um2"]),
            "section_major_um_median": _nanmed(cl["major_um"]),
            "section_minor_um_median": _nanmed(cl["minor_um"]),
            "section_aspect_median": _nanmed(cl["aspect"]),
            "face_openings": cl.get("face_openings", {})}
        _p(f"wrote {pp.name}")

    # A mesh rebuilt without its curvature leaves a per-vertex file of the wrong
    # length beside it. `measure_core` blanks all four curv_* columns for that
    # and says nothing (P12), so delete it and say so loudly rather than let a
    # later run look successful.
    if "mesh" in steps and "curvature" not in steps:
        stale = _out(saved, "vessel_curvature", single_tag, ".npz")
        if stale.exists():
            stale.unlink()
            notes.append(f"WARN: deleted {stale.name} — it was computed on the "
                         f"previous mesh and its per-vertex length no longer "
                         f"matches. Re-run with --steps curvature.")
            _p(notes[-1])

    # ---- provenance ---------------------------------------------------------
    # `status` reads only `notes` — the degradations of this run's own outputs.
    # `findings` are reported and recorded alongside them but do not degrade it.
    status = "degraded" if notes else "ok"
    reported = notes + findings
    timings["total_s"] = time.perf_counter() - t0
    for tg, g in groups.items():
        _write_recon(saved, tg, tif, params, steps, recon, g["outputs"],
                     g["curv"], g["scores"], timings, status, reported)
    del mask, verts, faces, curv_runs, cl, groups, scores
    gc.collect()
    return status, reported


def _write_recon(saved: Path, tag: str | None, tif: Path, params: dict,
                 steps: tuple, recon: dict, outputs: dict,
                 curv_block: dict | None, scores: dict | None, timings: dict,
                 status: str, notes: list[str]) -> None:
    """Write `vessel_recon.json` beside the outputs.

    **The name must not end in `_settings.json`**: `batch_measure.py:351`,
    `batch_postprocess.py:50` and `saved_check.py` all glob that pattern to
    decide whether a folder has already been measured, and a recon file matching
    it would make every folder look measured.

    A run that did not touch a step keeps that step's recorded block rather than
    dropping it — a curvature-only sweep must not erase the mesh's provenance —
    and says which blocks it carried over.
    """
    path = _out(saved, "vessel_recon", tag, ".json")
    doc: dict = {}
    carried_from = None
    if path.exists():
        try:
            doc = json.loads(path.read_text())
            carried_from = path.name
        except Exception:
            doc = {}
    elif tag and (saved / "vessel_recon.json").exists():
        # A tagged curvature sweep over the in-place mesh: inherit that mesh's
        # provenance so the tagged file is self-contained.
        try:
            doc = json.loads((saved / "vessel_recon.json").read_text())
            carried_from = "vessel_recon.json"
        except Exception:
            doc = {}
    carried = [k for k in ("mesh", "mesh_inert", "min_size", "cadherin_labels",
                           "scale_um_zyx", "scale_source", "surface_curvature",
                           "centerline_graph", "scores") if k in doc]

    # `outputs` is keyed by filename, and a run only knows about the files it
    # wrote — so it must be **merged**, not replaced. Replacing it means a
    # `--steps curvature` run in place erases the mesh's recorded vertex counts,
    # and a `--steps score` run (which writes no file at all) erases every
    # entry, leaving a recon JSON that describes a mesh it no longer records the
    # size of. Entries whose file has since gone are dropped, so the block
    # self-heals rather than accumulating ghosts.
    #
    # **Only from this tag's own previous file.** The inheritance branch above
    # seeds a new tagged doc from `vessel_recon.json` so it is self-contained
    # about the *parameters* — but that file's `outputs` describe the canonical
    # mesh, which is a different mesh. Merging those in made a tagged recon
    # claim the canonical file's vertex counts, and `recon_summary` then
    # reported the wrong mesh's size (seen: a 166,721-vertex sub1.0 mesh
    # summarised as 8,966 verts).
    if carried_from == path.name:
        merged = {k: v for k, v in (doc.get("outputs") or {}).items()
                  if (saved / k).exists()}
        merged.update(outputs)
        outputs = merged

    doc.update({
        "tool": "pipeline.batch_mesh",
        "git_sha": _git_sha(),
        "written": datetime.now().strftime("%Y%m%d_%H%M%S"),
        "image": tif.name,
        "tag": tag,
        "steps": list(steps),
        "status": status,
        "notes": notes,
        "timings_s": {k: round(v, 3) for k, v in timings.items()},
        "outputs": outputs,
    })
    doc.update(recon)
    carried = [c for c in carried if c not in recon]
    if "mesh" in steps:
        doc["mesh"] = params["mesh"]
        doc["mesh_inert"] = params["mesh_inert"]
        carried = [c for c in carried if c not in ("mesh", "mesh_inert")]
    if curv_block is not None:
        doc["surface_curvature"] = curv_block
        carried = [c for c in carried if c != "surface_curvature"]
    elif "mesh" in steps:
        # The mesh moved and no curvature was computed for it, so any recorded
        # curvature parameters describe the file this run just deleted.
        doc.pop("surface_curvature", None)
        carried = [c for c in carried if c != "surface_curvature"]
    if "graph" in steps:
        doc["centerline_graph"] = params["graph"]
        carried = [c for c in carried if c != "centerline_graph"]
    if scores is not None:
        doc["scores"] = scores
        carried = [c for c in carried if c != "scores"]
    elif "mesh" in steps:
        # Same reasoning as the curvature block above: the mesh moved, so any
        # recorded scores describe the surface this run just replaced. A stale
        # score is worse than no score — it is the one number a later sweep
        # would rank on.
        doc.pop("scores", None)
        carried = [c for c in carried if c != "scores"]
    if carried:
        doc["carried_over"] = {"blocks": carried, "from": carried_from}
    else:
        doc.pop("carried_over", None)

    path.write_text(json.dumps(doc, indent=2, default=str))
    _p(f"wrote {path.name}")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def _targets(directory: Path, only: str | None) -> list[tuple[Path, Path]]:
    out = []
    for saved in sorted(directory.glob("*_saved")):
        if only and only not in saved.name:
            continue
        tif = directory / f"{saved.name[:-len('_saved')]}.tif"
        out.append((tif, saved))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, type=Path)
    ap.add_argument("--only", default=None,
                    help="only process _saved folders whose name contains this "
                         "substring")
    ap.add_argument("--steps", default=",".join(DEFAULT_STEPS),
                    help=f"comma-separated subset of {','.join(STEPS)} "
                         f"(default: {','.join(DEFAULT_STEPS)}). --steps "
                         f"curvature needs neither the image header nor the "
                         f"mask; --steps score builds nothing and writes only "
                         f"the 'scores' block of vessel_recon.json — it is "
                         f"opt-in because it costs more than the rebuild.")
    ap.add_argument("--tag", default=None,
                    help="write vessel_mesh__TAG.npz siblings instead of "
                         "overwriting in place. '{radius}' and '{smooth}' are "
                         "substituted per curvature run.")
    ap.add_argument("--from-settings", default=None,
                    help="seed every parameter from a measure settings JSON or "
                         "a vessel_recon.json; explicit flags still override.")

    g = ap.add_argument_group("mask")
    g.add_argument("--min-size", default=None,
                   help="cadherin despeckle cut: an integer, 'auto' "
                        "(= max(otsu, 12), what the GUI applies) or 'auto-raw' "
                        "(the Otsu cut with no floor, docs/ISSUES.md P19). "
                        "Default 12.")

    g = ap.add_argument_group(
        "mesh (defaults are the current widget defaults; no upper bounds — two "
        "of the reference image's values sit exactly at their slider ceilings)")
    g.add_argument("--method", default=None, choices=METHODS)
    g.add_argument("--orientation", default=None, choices=ORIENTATIONS)
    g.add_argument("--n-bins", type=int, default=None)
    g.add_argument("--r-smooth-bins", type=int, default=None)
    g.add_argument("--alpha2d-um", type=float, default=None)
    g.add_argument("--subsample-um", type=float, default=None)
    g.add_argument("--alpha3d-um", type=float, default=None)
    g.add_argument("--fill-lumen", dest="fill_lumen", action="store_true",
                   default=None)
    g.add_argument("--no-fill-lumen", dest="fill_lumen", action="store_false")
    g.add_argument("--close-um", type=float, default=None)
    g.add_argument("--min-component-um", type=float, default=None)
    # No upper bound on any of these: a slider ceiling is a widget limit, not a
    # physical one (CLAUDE.md, "headless equivalents of a slider").
    g.add_argument("--poisson-grid-um", type=float, default=None)
    g.add_argument("--poisson-smooth-um", type=float, default=None)
    g.add_argument("--poisson-iso-pct", type=float, default=None)
    g.add_argument("--alpha-raster-um", type=float, default=None,
                   help="grid the alpha solid is rasterized onto for "
                        "the 'alpha 3D solid -> poisson' method "
                        "(default: the poisson grid). No upper bound.")

    g = ap.add_argument_group(
        "cadherin field (docs/ISSUES.md P28) — what the method is handed, not "
        "which method runs")
    g.add_argument("--cadherin-field", default=None, choices=CADHERIN_FIELDS,
                   help=f"'{FIELD_MASK}' (default) is the size-filtered binary "
                        f"mask and decodes no image pixels. '{FIELD_SCORE}' "
                        f"decodes the substrate channel and builds "
                        f"denoise.local_contrast_score, rectified at "
                        f"--score-k-lo, so dim-but-real wall survives with a "
                        f"small positive value. Only the Poisson method reads "
                        f"the values; the interpolating methods take "
                        f"'score > 0'.")
    # No upper bound on either, deliberately (CLAUDE.md): a slider ceiling is a
    # widget limit, not a physical one.
    g.add_argument("--score-window-xy", type=int, default=None,
                   help="xy window for the local mean/std, in pixels (the "
                        "Cadherin tab's window_xy; z is scaled for anisotropy).")
    g.add_argument("--score-k-lo", type=float, default=None,
                   help="the z-score the field is hinged to zero below. "
                        "Measured, not chosen — see P28 step 1.")

    g = ap.add_argument_group("surface curvature (comma-separated = a sweep)")
    g.add_argument("--curv-radius-um", default=None)
    g.add_argument("--curv-smooth-iter", default=None)

    g = ap.add_argument_group(
        "mesh scores (--steps score; three numbers reported as a vector and "
        "never combined into one — docs/ISSUES.md I03)")
    g.add_argument("--score-radii", default=None,
                   help=f"comma-separated radius_um values for the scale-"
                        f"stability sweep. Default "
                        f"{SCORE_DEFAULTS['radii_um']} — the radii P24's "
                        f"plateau table was measured at.")
    g.add_argument("--score-max-verts", type=int, default=None,
                   help="cap the surface vertices the curvature sweep is "
                        "fitted over. Off by default: it is a cost control, "
                        "not a free approximation, and a decimated score is "
                        "comparable only to another at the same cap and seed.")
    g.add_argument("--score-max-query", type=int, default=None,
                   help="cap the cadherin voxels coverage is queried with. Off "
                        "by default; the quantiles are robust and converge "
                        "well before the full mask.")

    g = ap.add_argument_group(
        "3D centerline graph (--steps graph; built on the Poisson lumen, so it "
        "needs a Poisson-family --method)")
    g.add_argument("--spur-ratio", type=float, default=None,
                   help="a leaf branch shorter than this many local inscribed "
                        "radii is a spur. A ratio, not a length, so it adapts "
                        "to a taper (docs/ISSUES.md P23).")
    g.add_argument("--crop-ratio", type=float, default=None,
                   help="an endpoint within this many local radii of a "
                        "field-of-view opening is a crop, not a vessel "
                        "terminus.")
    g.add_argument("--crop-margin-um", type=float, default=None,
                   help="absolute floor under --crop-ratio, in um.")
    g.add_argument("--kappa-window-um", type=float, default=None,
                   help="arc-length window for the 3D path curvature fit.")
    g.add_argument("--section-every-um", type=float, default=None,
                   help="cut a perpendicular cross-section every this many um "
                        "of arc length. Between cuts the section columns are "
                        "NaN, never interpolated.")
    g.add_argument("--section-span", type=float, default=None,
                   help="sample the cut plane out to this many local inscribed "
                        "radii. The inscribed radius is the *minor* semi-axis, "
                        "so this must exceed the aspect ratio or every section "
                        "is clipped.")

    g = ap.add_argument_group("modes")
    g.add_argument("--dry-run", action="store_true",
                   help="resolve and print the parameters (per image: the "
                        "scale, the shape check and 'auto' min_size are read "
                        "from the header and the labels file) and exit.")
    g.add_argument("--check-only", action="store_true",
                   help="describe every *_saved folder and exit, including "
                        "which meshes have no recorded parameters.")
    g.add_argument("--suggest-alpha", action="store_true",
                   help="one Delaunay per image: print a suggested alpha3d_um "
                        "for the current subsample_um, write nothing.")
    g.add_argument("--strict", action="store_true",
                   help="exit 2 if any image came out degraded.")
    args = ap.parse_args()

    directory = Path(args.dir)
    if not directory.is_dir():
        raise SystemExit(f"--dir: {directory} is not a directory.")

    if args.check_only:
        folders = [s for s in sorted(directory.glob("*_saved"))
                   if not args.only or args.only in s.name]
        if not folders:
            raise SystemExit(f"No *_saved folders in {directory}.")
        print(f"Pre-flight over {len(folders)} saved folder(s) in "
              f"{directory}:\n", flush=True)
        n_warn = 0
        for saved in folders:
            lines = describe_saved(saved)
            n_warn += sum(1 for ln in lines if "WARN:" in ln)
            print("\n".join(lines) + "\n", flush=True)
        print(f"{len(folders)} folder(s), {n_warn} warning(s). Nothing built "
              f"(--check-only).", flush=True)
        sys.exit(2 if (n_warn and args.strict) else 0)

    steps = tuple(s.strip() for s in args.steps.split(",") if s.strip())
    bad = [s for s in steps if s not in STEPS]
    if bad or not steps:
        raise SystemExit(f"--steps: expected a subset of {','.join(STEPS)}, "
                         f"got {args.steps!r}")
    steps = tuple(s for s in STEPS if s in steps)     # canonical order

    params = resolve_params(args)
    combos = params["curv_combos"]
    tag = args.tag
    if len(combos) > 1:
        if "mesh" in steps or "graph" in steps or "score" in steps:
            raise SystemExit(
                "a curvature sweep varies only the curvature: run "
                "'--steps mesh,path' once, then sweep with '--steps curvature "
                "--curv-radius-um ... --tag r{radius}'. (Scoring already "
                "sweeps radius internally — see --score-radii — so it would "
                "recompute the same scores once per tag.)")
        if not tag or not ("{radius}" in tag or "{smooth}" in tag):
            raise SystemExit(
                f"{len(combos)} curvature runs need a --tag containing "
                f"{{radius}} and/or {{smooth}}, or they would overwrite each "
                f"other.")
        resolved = [tag.format(radius=r, smooth=s) for r, s in combos]
        if len(set(resolved)) != len(resolved):
            raise SystemExit(f"--tag {tag!r} does not distinguish the "
                             f"{len(combos)} runs: {resolved}")

    print("Resolved parameters:", flush=True)
    for line in format_params(params):
        print(f"  {line}", flush=True)
    warn = parameter_warnings(params)
    for w in warn:
        print(f"  {w}", flush=True)
    print(f"  steps                {','.join(steps)}", flush=True)
    if steps == ("score",):
        # Scoring writes no array file at all, so "overwrites vessel_*.npz"
        # would be a false warning about the one run that cannot lose data.
        where = ("vessel_recon.json" if not tag
                 else f"vessel_recon__{tag}.json")
        print(f"  output               scores only, into {where} "
              f"(no mesh is built and no .npz is written)", flush=True)
    else:
        print(f"  output               "
              f"{'in place (overwrites vessel_*.npz)' if not tag else f'siblings tagged {tag!r}'}",
              flush=True)
    if "mesh" in steps and "curvature" not in steps and not args.dry_run:
        print("  NOTE: rebuilding the mesh without its curvature — any "
              "existing vessel_curvature.npz for that mesh will be DELETED "
              "(it would be the wrong length; docs/ISSUES.md P12).", flush=True)

    targets = _targets(directory, args.only)
    if not targets:
        raise SystemExit(f"No *_saved folders in {directory}"
                         + (f" matching {args.only!r}." if args.only else "."))
    missing = [t for t, _ in targets if not t.exists()]
    if missing and ("mesh" in steps or "graph" in steps):
        for t in missing:
            print(f"  skip: image {t.name} not found (its header is where the "
                  f"scale comes from).", flush=True)
        targets = [(t, s) for t, s in targets if t.exists()]

    verb = ("Resolving" if args.dry_run else
            "Suggesting alpha for" if args.suggest_alpha else "Rebuilding")
    print(f"\n{verb} {len(targets)} image(s):", flush=True)
    outcome: dict[str, list] = {"ok": [], "degraded": [], "failed": []}
    all_notes: list[tuple[str, list[str]]] = []
    for i, (tif, saved) in enumerate(targets):
        print(f"  [{i + 1}/{len(targets)}] {saved.name} ...", flush=True)
        try:
            status, notes = run_one(tif, saved, params, steps, tag,
                                    dry_run=args.dry_run,
                                    suggest_alpha=args.suggest_alpha)
        except Exception as exc:
            import traceback
            print(f"    FAILED: {exc}\n{traceback.format_exc()}", flush=True)
            outcome["failed"].append(saved.name)
            continue
        outcome[status].append(saved.name)
        if status == "failed":
            print(f"    FAILED: {notes[0] if notes else 'unknown'}", flush=True)
        if notes:
            all_notes.append((saved.name, notes))
        gc.collect()

    print(f"\n{len(outcome['ok'])} ok, {len(outcome['degraded'])} degraded, "
          f"{len(outcome['failed'])} failed.", flush=True)
    for name, notes in all_notes:
        print(f"\n  {name}:", flush=True)
        for n in notes:
            print(f"    {n}", flush=True)

    if outcome["failed"]:
        sys.exit(1)
    if args.strict and outcome["degraded"]:
        print("\n--strict: exiting 2 because some images came out degraded.",
              flush=True)
        sys.exit(2)


if __name__ == "__main__":
    main()
