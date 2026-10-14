"""
Convert ND2 z-stacks to OME-TIFF, dropping junk channels.

Reads each .nd2 in 3d_images_folder/images/WT/, drops the channels
not listed as useful in params.json, and writes a sibling .tif (OME-TIFF,
uint16, zstd compressed) carrying pixel size and channel names in metadata.
ND2 files are never modified; existing .tif siblings are overwritten.

The conversion itself lives in `pipeline/convert.py`, shared with the napari
Convert tab. Before converting, every file is checked there: one that cannot be
read, is multi-position, or whose ND2 channel names differ from the first
convertible file's is skipped and reported instead of stopping the run or being
written with the wrong channel names.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline.convert import (  # noqa: E402
    channel_guard, convert_one, inspect_nd2, keep_channel_plan,
)

DATA_ROOT = ROOT / "3d_images_folder"
PARAMS_PATH = DATA_ROOT / "params.json"


def load_params() -> dict:
    with PARAMS_PATH.open() as f:
        return json.load(f)


def update_params(params: dict) -> dict:
    """Drop per-channel `index` fields; OME ordering is now authoritative."""
    new_channels = [{"name": c["name"], "role": c["role"]} for c in params["channels"]]
    params["channels"] = new_channels
    return params


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true", help="List actions without writing")
    ap.add_argument("--subdir", default=None, help="Only convert this subdirectory of images/")
    args = ap.parse_args()

    params = load_params()
    images_root = DATA_ROOT / "images"
    subdirs = [images_root / args.subdir] if args.subdir else sorted(p for p in images_root.iterdir() if p.is_dir())

    nd2_files: list[Path] = []
    for d in subdirs:
        nd2_files.extend(sorted(d.glob("*.nd2")))
    if not nd2_files:
        print("No ND2 files found.")
        return

    # Header-check every file first; the channel plan comes from the first
    # convertible one and every other file must carry the same ND2 channels.
    infos = [inspect_nd2(f) for f in nd2_files]
    usable = [i for i in infos if i["ok"]]
    if not usable:
        for i in infos:
            print(f"SKIP {i['path'].relative_to(DATA_ROOT)}: {i['reason']}")
        print("No convertible ND2 files found.")
        sys.exit(1)
    nd2_channel_names = usable[0]["channels"]
    keep_idx, out_channels = keep_channel_plan(params, nd2_channel_names)
    out_names = [c["name"] for c in out_channels]
    print(f"ND2 channels: {nd2_channel_names}")
    print(f"Keeping indices {keep_idx} -> {out_names}")

    converted, skipped = 0, 0
    for info in infos:
        f = info["path"]
        reason = info["reason"] or channel_guard(nd2_channel_names, info["channels"])
        if reason:
            print(f"SKIP {f.relative_to(DATA_ROOT)}: {reason}")
            skipped += 1
            continue
        out = f.with_suffix(".tif")
        if args.dry_run:
            print(f"DRY: {f} -> {out}")
            continue
        print(f"Converting {f.relative_to(DATA_ROOT)} ...", flush=True)
        convert_one(f, keep_idx, out_names)
        print(f"  wrote {out.relative_to(DATA_ROOT)} ({out.stat().st_size/1e6:.1f} MB)")
        converted += 1

    if not args.dry_run and converted:
        params = update_params(params)
        with PARAMS_PATH.open("w") as f:
            json.dump(params, f, indent=2)
        print(f"Updated {PARAMS_PATH.relative_to(ROOT)} (dropped per-channel index fields)")
    if skipped:
        print(f"{skipped} file(s) skipped - see SKIP lines above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
