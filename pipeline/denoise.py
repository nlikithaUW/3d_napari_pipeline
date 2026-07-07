"""Display-only 3D denoise helpers for the napari viewer.

Anisotropy: every kernel takes its xy size/sigma in xy pixels and the z
component is rescaled so the kernel is roughly isotropic in microns.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu, threshold_triangle

from . import gpu


def _xy_um(voxel_zyx_um: tuple[float, float, float]) -> float:
    _, y, x = voxel_zyx_um
    return 0.5 * (y + x)


def aniso_sigma(sigma_xy: float, voxel_zyx_um) -> tuple[float, float, float]:
    z_um = voxel_zyx_um[0]
    sigma_z = sigma_xy * _xy_um(voxel_zyx_um) / z_um
    return (sigma_z, sigma_xy, sigma_xy)


def aniso_size(size_xy: int, voxel_zyx_um) -> tuple[int, int, int]:
    z_um = voxel_zyx_um[0]
    size_z = max(1, int(round(size_xy * _xy_um(voxel_zyx_um) / z_um)))
    if size_z % 2 == 0:
        size_z += 1
    return (size_z, size_xy, size_xy)


def _mean_std(v: np.ndarray, win) -> tuple[np.ndarray, np.ndarray]:
    """Local mean and std over the `win` footprint (float32).

    Delegates to `gpu.local_mean_std`, which keeps the whole computation on the
    GPU (one transfer in) when available and falls back to scipy's O(1)
    ``uniform_filter`` otherwise. Returns (mean, std) in original intensity units.
    """
    return gpu.local_mean_std(v, win)


def otsu_mask(vol: np.ndarray) -> np.ndarray:
    thr = threshold_otsu(vol)
    return (vol > thr).astype(np.uint8)


def label_and_sizes(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Label 3D connected components and return (labels, sizes_per_component).

    sizes is 0-indexed by component (so size of label k is sizes[k-1]).
    """
    labels, _ = ndi.label(mask.astype(bool))
    counts = np.bincount(labels.ravel())
    sizes = counts[1:] if counts.size > 1 else np.array([], dtype=np.int64)
    return labels, sizes


def auto_min_size_otsu(sizes: np.ndarray) -> int:
    """Otsu cutoff on log10(component sizes). Returns 0 if too few components."""
    if sizes.size < 4:
        return 0
    log_sizes = np.log10(sizes.astype(np.float64))
    log_thr = threshold_otsu(log_sizes)
    return int(np.ceil(10 ** log_thr))


def filter_by_size(labels: np.ndarray, sizes: np.ndarray, min_size: int) -> np.ndarray:
    keep = np.zeros(sizes.size + 1, dtype=bool)
    keep[1:] = sizes >= max(int(min_size), 0)
    return keep[labels].astype(np.uint8)


def robust_background(vol, subsample=8):
    """Global (median, MAD·1.4826) of intensities, from a strided subsample
    for speed. Background is uniform enough that subsampling is safe."""
    v = np.asarray(vol).ravel()[::max(1, subsample)].astype(np.float32)
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med)) * 1.4826)
    return med, mad


def robust_floor_value(vol, k):
    """Absolute intensity floor = median + k·MAD. Returns None if k <= 0
    (floor disabled)."""
    if k <= 0:
        return None
    med, mad = robust_background(vol)
    return med + float(k) * mad


def suggest_floor_k(vol, k_max=10.0):
    """Auto seed for the floor_k slider: the k that reproduces
    threshold_triangle, i.e. (triangle - median)/MAD, clamped to [0, k_max]."""
    v = np.asarray(vol).ravel()[::8].astype(np.float32)
    med, mad = robust_background(vol)
    if mad <= 0:
        return 0.0
    t = float(threshold_triangle(v))
    return float(min(max((t - med) / mad, 0.0), k_max))


def niblack_mask(vol: np.ndarray, window_xy: int, k: float, voxel_zyx_um) -> np.ndarray:
    """Local adaptive threshold: foreground = vol > local_mean + k * local_std.

    k is in units of local std-devs. skimage's threshold_niblack returns
    m - k*s (built for dark-text binarization); we negate k so this follows
    the "+ k*s" convention that matches fluorescence intuition (higher k →
    stricter → less foreground).
    """
    if window_xy % 2 == 0:
        window_xy += 1
    win = aniso_size(window_xy, voxel_zyx_um)
    m, s = _mean_std(vol, win)
    thr = m + float(k) * s
    return (vol > thr).astype(np.uint8)


def sauvola_mask(vol: np.ndarray, window_xy: int, k: float, voxel_zyx_um) -> np.ndarray:
    """Sauvola local threshold: foreground = vol > m·(1 + k·(s/R − 1)).

    A Niblack variant that normalizes the std term by the local dynamic range
    R, so it is more robust to noisy/uneven background. Here R is set to half
    the volume's intensity range. Note the k sign is opposite to our Niblack:
    higher k → lower threshold → MORE foreground (so lower k for less noise).
    """
    if window_xy % 2 == 0:
        window_xy += 1
    win = aniso_size(window_xy, voxel_zyx_um)
    v = vol.astype(np.float32, copy=False)
    r = max(0.5 * float(v.max() - v.min()), 1e-6)
    m, s = _mean_std(v, win)
    thr = m * (1.0 + float(k) * (s / r - 1.0))
    return (vol > thr).astype(np.uint8)


def bernsen_mask(vol: np.ndarray, window_xy: int, contrast_frac: float, voxel_zyx_um) -> np.ndarray:
    """Bernsen local threshold with a contrast guard.

    Per voxel over a local window: contrast = max − min, midgray = (max+min)/2.
    Where local contrast ≥ `contrast_frac`·(global range), foreground = vol >
    midgray; in low-contrast (flat) windows the whole neighborhood is assigned
    via the global Otsu threshold (midgray > otsu). The guard avoids
    hallucinating walls in flat background — useful for speckly cadherin.
    `contrast_frac` is a fraction of the global intensity range (0–1), so the
    knob is data-scale independent: higher → more regions treated as flat.
    """
    if window_xy % 2 == 0:
        window_xy += 1
    win = aniso_size(window_xy, voxel_zyx_um)
    v = vol.astype(np.float32, copy=False)
    lo = gpu.minimum_filter(v, size=win)
    hi = gpu.maximum_filter(v, size=win)
    contrast = hi - lo
    midgray = 0.5 * (hi + lo)
    contrast_thresh = float(contrast_frac) * float(v.max() - v.min())
    global_thr = threshold_otsu(v) if v.max() > v.min() else float(v.max())
    fg = np.where(contrast >= contrast_thresh, v > midgray, midgray > global_thr)
    return fg.astype(np.uint8)
