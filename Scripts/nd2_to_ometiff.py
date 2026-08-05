"""
Convert ND2 z-stacks to OME-TIFF, dropping junk channels.

Reads each .nd2 in 3d_images_folder/images/WT/, drops the channels
not listed as useful in params.json, and writes a sibling .tif (OME-TIFF,
uint16, zstd compressed) carrying pixel size and channel names in metadata.
ND2 files are never modified; existing .tif siblings are overwritten.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import nd2
import numpy as np
import tifffile

ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = ROOT / "3d_images_folder"
PARAMS_PATH = DATA_ROOT / "params.json"


def load_params() -> dict:
    with PARAMS_PATH.open() as f:
        return json.load(f)


def keep_channel_plan(params: dict, nd2_channel_names: list[str]) -> tuple[list[int], list[dict]]:
    """
    Decide which ND2 channel indices to keep and produce the post-convert
    channel list (name + role, in OME order). Matches params.json entries to
    ND2 channels by the explicit `index` field.
    """
    nd2_lower = [n.lower() for n in nd2_channel_names]
    keep_idx: list[int] = []
    out_channels: list[dict] = []
    for entry in params["channels"]:
        idx = entry["index"]
        if idx >= len(nd2_channel_names):
            raise ValueError(f"params channel index {idx} out of range for ND2 ({len(nd2_channel_names)} channels)")
        if nd2_lower[idx] == "junk":
            raise ValueError(f"params lists channel index {idx} but ND2 marks it junk")
        keep_idx.append(idx)
        out_channels.append({"name": entry["name"], "role": entry["role"]})
    return keep_idx, out_channels


def convert_one(nd2_path: Path, keep_idx: list[int], channel_names: list[str]) -> Path:
    out_path = nd2_path.with_suffix(".tif")
    with nd2.ND2File(str(nd2_path)) as nd:
        sizes = nd.sizes
        if set(sizes.keys()) - {"Z", "C", "Y", "X"}:
            raise ValueError(f"Unexpected axes in {nd2_path.name}: {sizes}")
        vox = nd.voxel_size()  # (x, y, z) microns
        arr = nd.asarray()  # (Z, C, Y, X) uint16

    if arr.ndim != 4:
        raise ValueError(f"Expected 4D ZCYX, got shape {arr.shape} for {nd2_path.name}")
    arr = arr[:, keep_idx, :, :]
    assert arr.dtype == np.uint16, f"unexpected dtype {arr.dtype}"

    ome_metadata = {
        "axes": "ZCYX",
        "PhysicalSizeX": float(vox.x),
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": float(vox.y),
        "PhysicalSizeYUnit": "µm",
        "PhysicalSizeZ": float(vox.z),
        "PhysicalSizeZUnit": "µm",
        "Channel": {"Name": channel_names},
    }

    tifffile.imwrite(
        str(out_path),
        arr,
        photometric="minisblack",
        compression="zstd",
        compressionargs={"level": 5},
        ome=True,
        metadata=ome_metadata,
        bigtiff=True,
    )
    return out_path


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

    # Resolve channel plan from the first ND2 (assume consistent acquisition).
    with nd2.ND2File(str(nd2_files[0])) as nd:
        nd2_channel_names = [c.channel.name for c in nd.metadata.channels]
    keep_idx, out_channels = keep_channel_plan(params, nd2_channel_names)
    out_names = [c["name"] for c in out_channels]
    print(f"ND2 channels: {nd2_channel_names}")
    print(f"Keeping indices {keep_idx} -> {out_names}")

    for f in nd2_files:
        out = f.with_suffix(".tif")
        if args.dry_run:
            print(f"DRY: {f} -> {out}")
            continue
        print(f"Converting {f.relative_to(DATA_ROOT)} ...", flush=True)
        convert_one(f, keep_idx, out_names)
        print(f"  wrote {out.relative_to(DATA_ROOT)} ({out.stat().st_size/1e6:.1f} MB)")

    if not args.dry_run:
        params = update_params(params)
        with PARAMS_PATH.open("w") as f:
            json.dump(params, f, indent=2)
        print(f"Updated {PARAMS_PATH.relative_to(ROOT)} (dropped per-channel index fields)")


if __name__ == "__main__":
    main()
