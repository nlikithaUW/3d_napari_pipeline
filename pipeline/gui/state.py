"""The shared session dict and the pipeline step table.

These are the two things `view.py` (the napari app) and `batch.py` (the headless
runner) both need, and each used to carry its own copy of: the initial `state`
dict was written out three times and the run-step list twice. That is the drift
CLAUDE.md's working agreement warns about — a step added to the GUI's "Run all"
but not to the batch's list is silently skipped headless, with no error.

No Qt and no napari here: this is data plus one parser.
"""

from __future__ import annotations

# --------------------------------------------------------------------------- #
# The session dict
# --------------------------------------------------------------------------- #
# The nested caches are built by their own small factories because the per-image
# reset in `batch.py` re-initialises them individually — it must NOT re-run
# `new_state`, which would also drop the loaded RF model, the scale and the
# tabs' histogram widgets.


def _dapi_state() -> dict:
    return {"labels": None, "sizes": None, "nuclei": None, "nuclei_sizes": None,
            "nuclei_ids": None, "nuclei_is_fib": None,
            "nuclei_min_size_vox": None}


def _cadh_state() -> dict:
    return {"labels": None, "sizes": None, "centroids": None, "normals": None,
            "label_ids": None, "keep_mask": None}


def _clover_state() -> dict:
    return {"labels": None, "sizes": None, "method": None}


def _nuc_cache() -> dict:
    return {"labels_id": None, "centroids": None, "normals": None, "lids": None}


def _mesh_cache() -> dict:
    return {"offset_centroids": None, "outward": None, "inward_offset_um": 0.0,
            "kdtree": None}


def new_state(*, headless: bool = False) -> dict:
    """A fresh session dict, as every tab's ``build`` expects to find it.

    ``headless=True`` is the batch runner's only difference: the Measure tab
    reads it to skip the display-only plotting.
    """
    state: dict = {
        "dapi_raw": None, "cadh_raw": None, "clover_raw": None,
        "force_vol": None, "adhesion_rf": None, "adhesion_proba": None,
        "scale": (1.0, 1.0, 1.0),
        "dapi_mask_layer": None, "cadh_mask_layer": None,
        "clover_mask_layer": None,
        "dapi_state": _dapi_state(),
        "cadh_state": _cadh_state(),
        "clover_state": _clover_state(),
        "nuc_cache": _nuc_cache(),
        "mesh_cache": _mesh_cache(),
        "vessel_mesh_cache": {"verts": None, "faces": None},
        "raw_layers": [],
        "hists": [],
    }
    if headless:
        state["headless"] = True
    return state


# Everything derived from one image, cleared before the next. Deliberately NOT
# everything in `new_state`: the RF model, the tabs' histogram widgets and the
# save-path registries outlive an image and re-creating them would break the
# batch. Belt and braces with `batch.py::_fail` — nothing derived from one image
# may survive into the next (`docs/ISSUES.md` P38).
PER_IMAGE_KEYS = (
    "dapi_raw", "cadh_raw", "clover_raw", "ruby_raw",
    "clover_raw_native", "ruby_raw_native", "force_vol",
    "adhesion_proba", "clover_ruby_ratio", "vessel_graph",
    "vessel_lumen_field", "vessel_mesh_cache",
    "measure_result", "adhesion_orient",
    "dapi_mask_layer", "cadh_mask_layer",
    "clover_mask_layer", "fibroblast_count",
)


def reset_image_state(state: dict) -> None:
    """Clear everything derived from the current image, in place."""
    for k in PER_IMAGE_KEYS:
        state[k] = None
    state["dapi_state"] = _dapi_state()
    state["cadh_state"] = _cadh_state()
    state["clover_state"] = _clover_state()
    state["vessel_mesh_cache"] = {"verts": None, "faces": None}
    state["mesh_cache"] = _mesh_cache()
    state["raw_layers"] = []


# --------------------------------------------------------------------------- #
# The pipeline step table
# --------------------------------------------------------------------------- #
# (short name, human label, state key). The short name is what `batch --steps`
# takes; it is stable and appears in scripts and in `docs/ISSUES.md`, so
# renaming the human label does not break a recorded invocation. Order is the
# pipeline's own and is the only order either driver runs them in.

RUN_STEPS = [
    ("dapi", "DAPI + components", "run_dapi"),
    ("fibroblasts", "Fibroblast count", "run_fibroblasts"),
    ("cadherin", "VE-Cadherin threshold", "run_cadherin"),
    ("mesh", "Vessel mesh", "run_vessel_mesh"),
    ("curvature", "Surface curvature", "run_vessel_curvature"),
    ("graph", "3D centerline graph", "run_vessel_path_curvature"),
    ("clover", "Clover threshold", "run_clover"),
    ("rf", "RF predict", "run_rf_predict"),
    ("postprocess", "RF post-process (split-merged)", "run_rf_postprocess"),
    ("measure", "Measure + export CSV", "run_measure"),
]

# Measure is intentionally absent from the batch: it crashes headless (a GL/Qt
# draw in the measure tab). The batch produces every intermediate; run
# `pipeline.batch_measure`, or the GUI's Load all saved -> Measure -> export,
# to get the CSVs.
BATCH_STEPS = [s for s in RUN_STEPS if s[0] != "measure"]
STEP_NAMES = [s[0] for s in BATCH_STEPS]


def resolve_steps(spec: str) -> list[str]:
    """Parse a `--steps` spec into short names, in pipeline order.

    Accepts "all", a comma-separated list, or a "first:last" range. Order is
    always the pipeline's own, never the order given — a later step consuming
    an earlier one's output cannot be run first by typo.

    **Every step still reads the tab defaults**, so re-running one is only
    reproducible if the defaults have not moved since (`docs/ISSUES.md` P26).
    """
    spec = (spec or "all").strip()
    if spec in ("", "all"):
        return list(STEP_NAMES)
    if ":" in spec:
        lo, _, hi = spec.partition(":")
        lo, hi = lo.strip() or STEP_NAMES[0], hi.strip() or STEP_NAMES[-1]
        for n in (lo, hi):
            if n not in STEP_NAMES:
                raise ValueError(f"unknown step {n!r}; expected one of "
                                 f"{', '.join(STEP_NAMES)}")
        i, j = STEP_NAMES.index(lo), STEP_NAMES.index(hi)
        if i > j:
            raise ValueError(f"empty range {spec!r}: {lo} comes after {hi}")
        return STEP_NAMES[i:j + 1]
    want = {s.strip() for s in spec.split(",") if s.strip()}
    bad = want - set(STEP_NAMES)
    if bad:
        raise ValueError(f"unknown step(s) {', '.join(sorted(bad))}; expected "
                         f"one of {', '.join(STEP_NAMES)}")
    return [n for n in STEP_NAMES if n in want]
