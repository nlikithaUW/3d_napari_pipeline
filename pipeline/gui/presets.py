"""Settings presets: every tab's tunable controls, saved and re-applied by name.

`batch.py` runs "at widget defaults" and a re-run only reproduces an earlier
one if those defaults have not moved ([P26], [I04]) — and a dataset at another
voxel size (the 20x V04 round) needs different values in several tabs at once.
A preset is a JSON of ``{"<Tab>/<control>": value}`` for every value control
in the tabs; applying it sets them all, so "Run all" (GUI) and `pipeline.batch
--preset` (headless) run with exactly the recorded values.

Controls are **discovered, not registered**: `collect` walks each tab's magicgui
widget tree (into `collapsible_group` containers too), so a control added to a
tab later is in every preset without touching this file. Excluded: buttons,
read-only labels, the per-image path fields in ``state["save_path_widgets"]``
(they point into one image's ``_saved`` folder), and the Load tab's
folder/experiment/image pickers. Keys are ``"<Tab>/<name>"``; a name repeated
within a tab gets ``#2``, ``#3`` in widget order.

`settings_from_export` maps an exported ``*_settings.json`` (its ``settings``
sections) back onto controls by name, so an earlier export's parameters can be
restored. Best effort: every key that does not land on exactly one control is
reported, never guessed.
"""

from __future__ import annotations

import json
from pathlib import Path

from magicgui.widgets import Container, FileEdit, Label, PushButton
from magicgui.widgets.bases import ValueWidget

PRESET_DIR_NAME = "presets"

# Load-tab pickers that describe *which image*, not *how to process it*, and
# the preset controls themselves.
_LOAD_EXCLUDE = {"project_dir", "experiment", "image", "settings preset",
                 "new preset name", "export settings JSON"}

# Which tab each section of an export's settings JSON belongs to.
_EXPORT_SECTION_TAB = {
    "cadherin": "Cadherin", "clover": "Clover", "fibroblast": "Nuclei",
    "vessel_mesh": "Vessel mesh", "surface_curvature": "Vessel mesh",
    "centerline_graph": "Vessel mesh", "adhesion_postprocess": "Adhesion RF",
    "rf_predict": "Adhesion RF",
}


# Export keys whose control has a different name.
_EXPORT_ALIASES = {
    ("surface_curvature", "radius_um"): "curvature_radius_um",
    ("surface_curvature", "smooth_iter"): "curvature_smooth_iter",
}


def _name(w) -> str:
    # A CheckBox built with text= (e.g. "compute FRET") has no name or label.
    return (getattr(w, "name", "") or getattr(w, "label", "")
            or getattr(w, "text", "") or "").strip()


def _walk(obj):
    """Yield value widgets under a tab part (Container, FunctionGui,
    collapsible group, or a list of them)."""
    if isinstance(obj, (list, tuple)):
        for o in obj:
            yield from _walk(o)
        return
    inner = getattr(obj, "_mg_container", None)          # collapsible_group
    if inner is not None:
        yield from _walk(inner)
        return
    if isinstance(obj, FileEdit):                        # a Container with a value
        yield obj
        return
    if isinstance(obj, Container):
        for w in obj:
            yield from _walk(w)
        return
    if isinstance(obj, ValueWidget) and not isinstance(obj, (PushButton, Label)):
        yield obj


def collect(tabs: dict, state: dict) -> dict:
    """``{"<Tab>/<name>": widget}`` for every preset-able control in `tabs`
    (``{tab title: build() result}``)."""
    per_image = {id(w) for w, _ in state.get("save_path_widgets", [])}
    out: dict = {}
    for tab, parts in tabs.items():
        seen: dict = {}
        for w in _walk(parts):
            nm = _name(w)
            if not nm or id(w) in per_image:
                continue
            if tab == "Load" and nm in _LOAD_EXCLUDE:
                continue
            seen[nm] = seen.get(nm, 0) + 1
            key = f"{tab}/{nm}" + (f"#{seen[nm]}" if seen[nm] > 1 else "")
            out[key] = w
    return out


def _jsonable(v):
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, tuple):
        return list(v)
    return v


def snapshot(tabs: dict, state: dict) -> dict:
    """Current value of every preset-able control."""
    return {k: _jsonable(w.value) for k, w in collect(tabs, state).items()}


def apply(values: dict, tabs: dict, state: dict) -> tuple[int, list[str]]:
    """Set controls from `values`. Returns ``(n_set, problems)`` — a key with
    no control, or a value the control refuses (outside a slider's range, not
    among a dropdown's choices), is reported and skipped, never coerced."""
    widgets = collect(tabs, state)
    n, problems = 0, []
    for key, val in values.items():
        w = widgets.get(key)
        if w is None:
            problems.append(f"{key}: no such control (preset from another version?)")
            continue
        try:
            if w.value != val:
                w.value = val
            n += 1
        except Exception as exc:  # noqa: BLE001 - report every refusal, keep going
            problems.append(f"{key}: could not set {val!r} ({exc})")
    return n, problems


def preset_dir(project_dir) -> Path:
    return Path(project_dir) / PRESET_DIR_NAME


def list_presets(project_dir) -> list[str]:
    d = preset_dir(project_dir)
    return sorted(p.stem for p in d.glob("*.json")) if d.is_dir() else []


def save_preset(name: str, tabs: dict, state: dict, project_dir) -> Path:
    name = name.strip()
    if not name or any(c in name for c in '\\/:*?"<>|'):
        raise ValueError(f"invalid preset name {name!r}")
    d = preset_dir(project_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.json"
    path.write_text(json.dumps({"preset": name, "values": snapshot(tabs, state)},
                               indent=2))
    return path


def load_preset(name: str, project_dir) -> dict:
    data = json.loads((preset_dir(project_dir) / f"{name}.json").read_text())
    return data["values"]


def settings_from_export(path, tabs: dict, state: dict) -> tuple[dict, list[str]]:
    """Map an exported ``*_settings.json`` onto preset keys. Returns
    ``(values, unmapped)``; `values` can be passed to `apply`."""
    meta = json.loads(Path(path).read_text())
    settings = meta.get("settings", meta)
    keys = collect(tabs, state)
    values, unmapped = {}, []
    for section, body in settings.items():
        tab = _EXPORT_SECTION_TAB.get(section)
        if tab is None or not isinstance(body, dict):
            continue
        for k, v in body.items():
            if k in ("inert", "note"):
                continue
            ctl = _EXPORT_ALIASES.get((section, k), k)
            matches = [key for key in keys if key == f"{tab}/{ctl}"]
            if len(matches) == 1 and f"{tab}/{ctl}#2" not in keys:
                values[matches[0]] = v
            else:
                unmapped.append(f"{section}.{k}")
    return values, unmapped
