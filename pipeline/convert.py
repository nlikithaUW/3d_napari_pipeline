"""ND2 -> OME-TIFF conversion: the one copy both converters call.

`Scripts/nd2_to_ometiff.py` (headless) and the napari Convert tab
(`gui/convert.py`) both go through this module, the same doctrine as
`saved_check.py` and `vessel_mesh.build_mesh`: two front ends, one
implementation, so the file the GUI writes and the file the script writes
cannot differ.

Channel selection is **by ND2 index**, exactly as the script always did
(`keep_channel_plan`); the OME-TIFF carries the `params.json` names, which is
what `io.load_stack` matches against. What this module adds is a **guard**:
`inspect_nd2` says whether a file can be converted at all (unreadable, empty,
multi-position, non-uint16), and `channel_guard` refuses a file whose ND2
channel names differ from the ones the mapping was built from. Before it, the
plan was taken from the first file and applied to every other one unchecked —
a file acquired with its channels in another order was written with the wrong
names and no error.

No Qt, no napari. Never prints; callers report.
"""

from __future__ import annotations

from pathlib import Path

import nd2
import numpy as np
import tifffile

# The only axes `convert_one` (and downstream `io.load_stack`) handle. A
# multi-position (P) or time-series (T) file is refused, not split.
SUPPORTED_AXES = {"Z", "C", "Y", "X"}


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


def inspect_nd2(nd2_path: str | Path) -> dict:
    """Read an ND2's header (no pixels) and say whether `convert_one` can take it.

    Returns a dict with ``ok`` (bool), ``reason`` (``None`` when ok, else why
    not), ``sizes``, ``channels`` (ND2 channel names, in index order),
    ``emission_nm`` and ``voxel_zyx_um``. Never raises: an unreadable file is a
    result, not an exception, so one bad file in a folder cannot stop the rest.
    """
    p = Path(nd2_path)
    info = {"path": p, "ok": False, "reason": None, "sizes": {}, "channels": [],
            "emission_nm": [], "voxel_zyx_um": None}
    try:
        if p.stat().st_size == 0:
            info["reason"] = "empty file (0 bytes)"
            return info
        with nd2.ND2File(str(p)) as f:
            sizes = dict(f.sizes)
            info["sizes"] = sizes
            chans = f.metadata.channels or []
            info["channels"] = [c.channel.name for c in chans]
            info["emission_nm"] = [c.channel.emissionLambdaNm for c in chans]
            v = f.voxel_size()
            info["voxel_zyx_um"] = (float(v.z), float(v.y), float(v.x))
            dtype = np.dtype(f.dtype)
    except Exception as exc:  # noqa: BLE001 - any failure to parse is "unreadable"
        info["reason"] = f"unreadable ({type(exc).__name__})"
        return info

    if "P" in sizes:
        info["reason"] = f"multi-position (P={sizes['P']}) - not supported"
    elif set(sizes) != SUPPORTED_AXES:
        info["reason"] = f"unsupported axes {sizes} (need exactly Z, C, Y, X)"
    elif dtype != np.uint16:
        info["reason"] = f"dtype {dtype} (need uint16)"
    elif len(info["channels"]) != sizes["C"]:
        info["reason"] = (f"{len(info['channels'])} channel names for "
                          f"C={sizes['C']} channels")
    else:
        info["ok"] = True
    return info


def channel_guard(reference: list[str], actual: list[str]) -> str | None:
    """``None`` if a file's ND2 channel names match the mapping's, else why not.

    Exact, ordered comparison of the full name list: the plan picks channels by
    index, so a file whose channel at any index differs — reordered, renamed,
    added or missing — would be written with the wrong names.
    """
    if list(actual) == list(reference):
        return None
    return f"ND2 channels {list(actual)} differ from the mapping's {list(reference)}"
