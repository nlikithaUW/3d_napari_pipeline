"""Convert tab: turn the .nd2 files in an experiment folder into OME-TIFFs.

The GUI front end of `Scripts/nd2_to_ometiff.py`; both call
`pipeline/convert.py`, so the .tif written here and the one the script writes
are the same file. Workflow:

  1. Copy the .nd2 files into ``images/<experiment>/`` and pick that folder.
  2. Scan       -> every .nd2 is header-checked: OK, or why it will be skipped
                   (unreadable, empty, multi-position, channels that differ
                   from the mapping's).
  3. Map        -> for each channel in params.json, pick the ND2 channel (by
                   index) it comes from. The ND2's own names are shown to help;
                   the OME-TIFF gets the params.json name, which is what Load
                   matches against.
  4. Convert    -> writes ``<stem>.tif`` beside each selected .nd2, then
                   refreshes the Load tab's image list.

The mapping is saved beside the data as ``images/<experiment>/nd2_channels.json``
— the params.json ``index/name/role`` format plus the ND2 channel names it was
built from, which is what the per-file guard checks against. Unlike the script,
this tab never rewrites params.json.
"""

from __future__ import annotations

import json
from pathlib import Path

from magicgui.widgets import CheckBox, ComboBox, Label, PushButton, Select
from napari.qt.threading import thread_worker

from .. import convert as cv
from ..io import EXPERIMENT_PARAMS, load_params
from .widgets import collapsible_group

MAPPING_FILE = "nd2_channels.json"
NONE = "(none)"

# First-guess ND2 channel for each role, by substring of the ND2 channel name.
# Only a pre-selection: the user confirms every dropdown before converting.
_GUESS = {
    "nuclei": ("dapi", "hoechst"),
    "clover": ("clover", "gfp"),
    "ruby": ("ruby", "rfp"),
    "fret": ("fret",),
    "substrate": ("cadherin", "ve-cad", "647"),
}


def _choice(i: int, name: str, em) -> str:
    return f"{i}: {name}" + (f" ({em:.0f} nm)" if em else "")


def build(state, viewer, project_dir, params_path) -> list:
    images_root = Path(project_dir) / "images"
    params = load_params(params_path)
    # One dropdown per params.json channel, in params.json (= OME) order.
    out_channels = [{"name": c["name"], "role": c["role"]} for c in params["channels"]]

    def _experiments(*_):
        # A callable, not a list: magicgui re-evaluates it on every
        # reset_choices() (a Container does one when it adopts the widget), so
        # the folder list is always the disk's, never a stale snapshot.
        return sorted(p.name for p in images_root.iterdir() if p.is_dir()) \
            if images_root.exists() else []

    exp_w = ComboBox(choices=_experiments, label="experiment")
    # `reference`: the ND2 channel names the mapping refers to, by index.
    # `file_choices` / `map_opts` back the list widgets below — through
    # callables, because a plain list assigned to `.choices` is discarded by the
    # next reset_choices(), emptying the widget with no message.
    scan = {"infos": {}, "reference": None, "emission": [],
            "file_choices": [], "map_opts": [NONE]}
    refresh_w = PushButton(text="Refresh folders")
    scan_w = PushButton(text="Scan .nd2 files")
    files_w = Select(choices=lambda *_: scan["file_choices"], label="files")
    scan_status_w = Label(value="Pick an experiment folder and Scan.")
    map_ws = [ComboBox(choices=lambda *_: scan["map_opts"],
                       label=f"{c['name']} ({c['role']})")
              for c in out_channels]
    save_map_w = PushButton(text="Save mapping")
    overwrite_w = CheckBox(value=False, text="overwrite existing .tif")
    convert_w = PushButton(text="Convert selected")

    def _folder() -> Path | None:
        return images_root / exp_w.value if exp_w.value else None

    def _set_mapping_choices(reference, emission, saved):
        opts = [NONE] + [_choice(i, n, emission[i] if i < len(emission) else None)
                         for i, n in enumerate(reference)]
        by_role = {c["role"]: c["index"] for c in (saved or {}).get("channels", [])}
        scan["map_opts"] = opts
        for w, c in zip(map_ws, out_channels):
            w.reset_choices()
            idx = by_role.get(c["role"])
            if idx is None and not saved:
                idx = next((i for i, n in enumerate(reference)
                            if any(g in n.lower() for g in _GUESS.get(c["role"], ()))),
                           None)
            w.value = opts[idx + 1] if idx is not None and idx < len(reference) else NONE

    def _scan(*_):
        folder = _folder()
        if folder is None or not folder.is_dir():
            scan_status_w.value = "No experiment folder selected."
            return
        nd2s = sorted(folder.glob("*.nd2"))
        if not nd2s:
            scan["file_choices"] = []
            files_w.reset_choices()
            scan_status_w.value = f"No .nd2 files in images/{folder.name}."
            return
        # The OME names come from the global params.json; an experiment's own
        # params.json (io.resolve_params) must use the same ones or Load fails.
        exp_params = folder / EXPERIMENT_PARAMS
        if exp_params.is_file():
            try:
                exp_names = [c["name"] for c in load_params(exp_params)["channels"]]
            except Exception as exc:  # noqa: BLE001
                exp_names = None
                print(f"Convert: could not read {exp_params.name} ({exc}).", flush=True)
            if exp_names is not None and sorted(exp_names) != sorted(
                    c["name"] for c in out_channels):
                print(f"Convert: WARNING - images/{folder.name}/{EXPERIMENT_PARAMS} "
                      f"names channels {exp_names}, but converted files get "
                      f"{[c['name'] for c in out_channels]} (global params.json). "
                      "Load will not find them; make the names agree.", flush=True)
        infos = {p.name: cv.inspect_nd2(p) for p in nd2s}
        saved = None
        mfile = folder / MAPPING_FILE
        if mfile.exists():
            try:
                saved = json.loads(mfile.read_text())
            except Exception as exc:  # noqa: BLE001
                print(f"Convert: could not read {mfile.name} ({exc}); ignoring it.",
                      flush=True)
        ok = [i for i in infos.values() if i["ok"]]
        if saved and saved.get("nd2_channels"):
            reference = list(saved["nd2_channels"])
            source = f"saved mapping ({MAPPING_FILE})"
        elif ok:
            reference = ok[0]["channels"]
            source = f"first readable file ({ok[0]['path'].name})"
        else:
            reference = []
            source = "none"
        donor = next((i for i in ok if i["channels"] == reference), None)
        scan.update(infos=infos, reference=reference,
                    emission=donor["emission_nm"] if donor else [])

        choices, selectable = [], []
        n_skip = 0
        for name, info in infos.items():
            reason = info["reason"] or cv.channel_guard(reference, info["channels"])
            if reason:
                n_skip += 1
                choices.append((f"{name}  -  SKIP: {reason}", name))
                continue
            z = info["sizes"].get("Z")
            dz = info["voxel_zyx_um"][0]
            tag = "  [tif exists]" if (folder / name).with_suffix(".tif").exists() else ""
            choices.append((f"{name}  -  OK  Z={z}, z-step {dz:g} um{tag}", name))
            selectable.append(name)
        scan["file_choices"] = choices
        files_w.reset_choices()
        files_w.value = selectable
        _set_mapping_choices(reference, scan["emission"], saved)
        scan_status_w.value = (f"{len(selectable)} convertible, {n_skip} skipped. "
                               f"Channels from {source}.")
        print(f"Convert: scanned images/{folder.name} - {len(selectable)} "
              f"convertible, {n_skip} skipped. ND2 channels: {reference}", flush=True)
        for name, info in infos.items():
            reason = info["reason"] or cv.channel_guard(reference, info["channels"])
            if reason:
                print(f"  SKIP {name}: {reason}", flush=True)

    def _mapping() -> dict | None:
        """The current dropdowns as a params.json-style dict, or None if unusable."""
        if not scan["reference"]:
            print("Convert: Scan a folder first.", flush=True)
            return None
        chans = []
        for w, c in zip(map_ws, out_channels):
            if w.value == NONE:
                continue
            chans.append({"index": int(str(w.value).split(":", 1)[0]),
                          "name": c["name"], "role": c["role"]})
        if not chans:
            print("Convert: map at least one channel.", flush=True)
            return None
        idx = [c["index"] for c in chans]
        if len(set(idx)) != len(idx):
            print("Convert: the same ND2 channel is mapped twice - fix the "
                  "dropdowns.", flush=True)
            return None
        return {"nd2_channels": list(scan["reference"]), "channels": chans}

    def _save_mapping(*_):
        mapping = _mapping()
        folder = _folder()
        if mapping is None or folder is None:
            return
        (folder / MAPPING_FILE).write_text(json.dumps(mapping, indent=2))
        print(f"Convert: mapping saved -> images/{folder.name}/{MAPPING_FILE}",
              flush=True)

    def _convert(*_):
        mapping = _mapping()
        folder = _folder()
        if mapping is None or folder is None:
            return
        names = list(files_w.value)
        if not names:
            print("Convert: no files selected.", flush=True)
            return
        try:
            keep_idx, out = cv.keep_channel_plan(mapping, mapping["nd2_channels"])
        except ValueError as exc:
            print(f"Convert: {exc}", flush=True)
            return
        out_names = [c["name"] for c in out]
        reference = mapping["nd2_channels"]
        overwrite = bool(overwrite_w.value)
        _save_mapping()
        print(f"Convert: {len(names)} file(s), ND2 indices {keep_idx} -> "
              f"{out_names}", flush=True)

        @thread_worker
        def _run():
            done, skipped = 0, 0
            for k, name in enumerate(names, 1):
                src = folder / name
                dst = src.with_suffix(".tif")
                # Re-check at convert time: the folder may have changed since Scan.
                info = cv.inspect_nd2(src)
                reason = info["reason"] or cv.channel_guard(reference, info["channels"])
                if reason is None and dst.exists() and not overwrite:
                    reason = "tif already exists (tick 'overwrite' to replace)"
                if reason:
                    print(f"  [{k}/{len(names)}] SKIP {name}: {reason}", flush=True)
                    skipped += 1
                    continue
                print(f"  [{k}/{len(names)}] converting {name} ...", flush=True)
                cv.convert_one(src, keep_idx, out_names)
                print(f"      wrote {dst.name} ({dst.stat().st_size / 1e6:.0f} MB)",
                      flush=True)
                done += 1
            return done, skipped

        def _done(result):
            done, skipped = result
            print(f"Convert: finished - {done} written, {skipped} skipped.",
                  flush=True)
            refresh = state.get("refresh_image_list")
            if refresh is not None:
                refresh()
            _scan()

        def _err(exc):
            print(f"Convert: error - {type(exc).__name__}: {exc}", flush=True)

        w = _run()
        w.returned.connect(_done)
        w.errored.connect(_err)
        w.start()

    refresh_w.changed.connect(lambda *_: exp_w.reset_choices())
    scan_w.changed.connect(_scan)
    save_map_w.changed.connect(_save_mapping)
    convert_w.changed.connect(_convert)

    g1 = collapsible_group("1. Files", [exp_w, refresh_w, scan_w, scan_status_w,
                                        files_w], expanded=True)
    g2 = collapsible_group("2. Channel mapping (ND2 index -> params.json channel)",
                           map_ws + [save_map_w], expanded=True)
    g3 = collapsible_group("3. Convert", [overwrite_w, convert_w], expanded=True)
    return [g1, g2, g3]
