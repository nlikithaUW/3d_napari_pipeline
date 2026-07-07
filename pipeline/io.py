"""Load OME-TIFFs and address channels by role from params.json."""

from __future__ import annotations

import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import NamedTuple

import numpy as np
import tifffile

from . import gpu


class Stack(NamedTuple):
    data: np.ndarray            # (Z, C, Y, X) uint16
    channels: list[str]         # channel names in C order (from OME XML)
    roles: dict[str, str]       # role -> name (e.g. {"clover": "mClover"})
    voxel_zyx_um: tuple[float, float, float]

    def by_role(self, role: str) -> np.ndarray:
        name = self.roles[role]
        return self.data[:, self.channels.index(name)]


def _parse_ome_channel_names(xml: str) -> list[str]:
    return re.findall(r'<Channel[^/>]*\sName="([^"]+)"', xml)


def _parse_ome_voxel_zyx(xml: str) -> tuple[float, float, float]:
    def grab(attr: str) -> float:
        m = re.search(rf'{attr}="([0-9.eE+-]+)"', xml)
        return float(m.group(1)) if m else float("nan")
    return (grab("PhysicalSizeZ"), grab("PhysicalSizeY"), grab("PhysicalSizeX"))


def load_params(params_path: str | Path) -> dict:
    with open(params_path) as f:
        return json.load(f)


def load_stack(tif_path: str | Path, params: dict) -> Stack:
    """Read an OME-TIFF stack and return data + channel mapping built from params."""
    with tifffile.TiffFile(str(tif_path)) as tf:
        series = tf.series[0]
        if series.axes != "ZCYX":
            raise ValueError(f"{tif_path}: expected ZCYX axes, got {series.axes}")
        data = series.asarray()
        xml = tf.ome_metadata or ""

    channels = _parse_ome_channel_names(xml)
    if len(channels) != data.shape[1]:
        raise ValueError(f"{tif_path}: OME channel count {len(channels)} != data C {data.shape[1]}")
    voxel = _parse_ome_voxel_zyx(xml)

    roles = {entry["role"]: entry["name"] for entry in params["channels"]}
    missing = [r for r, n in roles.items() if n not in channels]
    if missing:
        raise ValueError(f"{tif_path}: roles {missing} not found in channels {channels}")

    return Stack(data=data, channels=channels, roles=roles, voxel_zyx_um=voxel)


def resample_volume(
    vol: np.ndarray,
    src_zyx_um: tuple[float, float, float],
    target_um: float,
    order: int = 1,
) -> np.ndarray:
    """Resample a 3D (Z,Y,X) volume to an isotropic voxel size in µm."""
    zf = (src_zyx_um[0] / target_um,
          src_zyx_um[1] / target_um,
          src_zyx_um[2] / target_um)
    if all(abs(f - 1.0) < 1e-9 for f in zf):
        return vol
    return gpu.zoom(vol, zf, order=order, mode="reflect")


def to_isotropic(stack: Stack, target_um: float | None = None) -> Stack:
    """Return a Stack resampled to isotropic voxels.

    Default target = min(voxel_zyx_um) from the source metadata, so the
    smallest existing axis sets the grid (preserves in-plane resolution).
    Each file's own OME metadata drives the resample — files with
    different z spacing each end up on their own isotropic grid.
    """
    src = stack.voxel_zyx_um
    if any(not np.isfinite(v) or v <= 0 for v in src):
        raise ValueError(
            f"Cannot resample: invalid voxel size {src}. "
            f"Check OME metadata (PhysicalSizeZ/Y/X)."
        )
    if target_um is None:
        target_um = float(min(src))

    Z, C, Y, X = stack.data.shape
    # Resample channel 0 first to discover output shape, then run the rest
    # in parallel — scipy.ndimage.zoom releases the GIL.
    first = resample_volume(stack.data[:, 0], src, target_um, order=1)
    out = np.empty((first.shape[0], C, first.shape[1], first.shape[2]),
                   dtype=stack.data.dtype)
    out[:, 0] = first.astype(stack.data.dtype, copy=False)
    del first

    def _do(c: int) -> None:
        rv = resample_volume(stack.data[:, c], src, target_um, order=1)
        out[:, c] = rv.astype(stack.data.dtype, copy=False)

    if C > 1:
        if gpu.gpu_enabled():
            # The GPU is already internally parallel; concurrent host threads
            # would just contend for the device and risk OOM. Go sequential.
            for c in range(1, C):
                _do(c)
        else:
            with ThreadPoolExecutor(max_workers=min(C - 1, 8)) as ex:
                list(ex.map(_do, range(1, C)))

    return Stack(
        data=out,
        channels=stack.channels,
        roles=stack.roles,
        voxel_zyx_um=(target_um, target_um, target_um),
    )
