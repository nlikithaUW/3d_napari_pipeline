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

    Implemented as a comparison against `_local_zscore` rather than against
    ``m + k·s`` so that it and `local_contrast_score` are the same statistic by
    construction. The two forms are algebraically identical and differ in
    float32 only at the boundary; no saved artifact is affected either way,
    because every cadherin mask in this project **up to 2026-08-24** was built
    with Sauvola: `_run_cadherin_step` forced it, overriding the tab's Niblack
    default, for both `batch.py` and the GUI's "Run all". That override is gone
    (P36 measured Sauvola as *inert* here), so masks built after that date use
    whatever the tab is set to.
    """
    z = _local_zscore(vol, window_xy, voxel_zyx_um)
    return (z > np.float32(k)).astype(np.uint8)


# The floor on the local std. Two jobs: it stops a divide-by-zero where the
# neighbourhood is perfectly flat (saturated or clipped background), and it caps
# how loud a single count of variation is allowed to be. In camera counts a
# standard deviation below this is quantisation, not signal.
_STD_FLOOR = 1e-3


def _local_zscore(vol: np.ndarray, window_xy: int, voxel_zyx_um) -> np.ndarray:
    """``(vol - m) / s`` over the anisotropy-corrected `window_xy` footprint.

    The single implementation of the Niblack statistic. `niblack_mask` and
    `local_contrast_score` both call it, so the threshold and the continuous
    field cannot be the same quantity computed two slightly different ways —
    they were, and the two disagreed on ~1e-6 of voxels purely from evaluating
    ``vol > m + k·s`` against ``(vol - m)/s > k`` in float32.
    """
    if window_xy % 2 == 0:
        window_xy += 1
    win = aniso_size(window_xy, voxel_zyx_um)
    m, s = _mean_std(vol, win)
    return (vol.astype(np.float32, copy=False) - m) / np.maximum(s, _STD_FLOOR)


def local_contrast_score(vol: np.ndarray, window_xy: int, voxel_zyx_um,
                         k_lo: float = 0.0) -> np.ndarray:
    """Niblack's statistic as a *continuous, rectified* field, not a mask.

    ``max((vol - m) / s - k_lo, 0)`` — how far each voxel stands above its own
    neighbourhood, in units of local variation, hinged to exactly zero below
    ``k_lo``. Returns float32.

    This is the same ``m``/``s`` `niblack_mask` computes, and it holds the
    identity that pins the two together::

        niblack_mask(vol, w, k) == (local_contrast_score(vol, w, vx, 0) > k)

    for ``k >= 0`` (see `docs/scratch/p28/score_check.py`). One implementation,
    two callers — the no-drift rule in `CLAUDE.md`.

    **Why continuous.** `docs/ISSUES.md` [P28] measured that the "gaps" in the
    cadherin wall mostly are not holes: the raw channel carries wall at the same
    radius, dimmer, and the threshold deletes it. A score keeps dim-but-real
    wall alive with a small positive value instead, and stays depth-corrected
    because ``m`` and ``s`` adapt with depth.

    **Why rectified.** `vessel_poisson` reads ``occ > 0`` as "this cell has
    wall". A raw z-score is nonzero nearly everywhere, so an unrectified field
    would put tens of millions of cells into `sheet_normals` and make
    `lumen_core`'s occupancy-weighted wall percentile meaningless. ``k_lo``
    hinges the field to zero exactly where the data says there is nothing —
    it is **measured**, not picked: [P28] step 1 records the z-score the gap
    wall sits at against the no-wall controls, and `batch_mesh` records which
    measurement set the value it used.

    **The global intensity floor is deliberately not applied here.** [P28]'s
    attribution found the floor, not the local threshold, is what deletes the
    dim wall; folding it back in would re-delete exactly what this exists to
    keep. A caller that wants it applies it itself.
    """
    z = _local_zscore(vol, window_xy, voxel_zyx_um)
    return np.maximum(z - np.float32(k_lo), np.float32(0.0)).astype(np.float32,
                                                                   copy=False)


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
    # Hand `vol` to `_mean_std` in its NATIVE dtype, not cast to float32 first.
    # `gpu.local_mean_std` transfers the volume as-is and casts on the device
    # precisely so a uint16 stack crosses PCIe at half the bytes — casting here
    # doubled the transfer and defeated that. On a 324x1024x1024 stack it is
    # 680 MB vs 1.36 GB, and the larger one exceeds the pinned-memory pool, so
    # the copy silently degraded to a synchronous transfer (and to a full CPU
    # fallback beyond that). Bit-identical either way: the first thing
    # `local_mean_std` does is `.astype(float32)`.
    r = max(0.5 * float(vol.max() - vol.min()), 1e-6)
    m, s = _mean_std(vol, win)
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


# ---------------------------------------------------------------------------
# Threshold + label + size-suggest, as a single dispatch point.
# ---------------------------------------------------------------------------
#
# This exists for the same reason `vessel_mesh.build_mesh` does: two callers
# needed the identical computation and each had written it out.
# `gui/clover.py::_compute` and `gui/adhesion_rf.py`'s gate fallback both did
# `otsu_mask` -> `label_and_sizes` -> `auto_min_size_otsu`, so the clover mask
# the user tuned in one tab and the one the RF gated on could differ without
# anything saying so. Both now call `build_mask`; do not write the third copy.

METHOD_OTSU = "otsu"
METHOD_NIBLACK = "niblack"
METHOD_SAUVOLA = "sauvola"
METHOD_BERNSEN = "bernsen"
MASK_METHODS = (METHOD_OTSU, METHOD_NIBLACK, METHOD_SAUVOLA, METHOD_BERNSEN)

# Which knobs each method actually reads. The rest are inert for that method
# and must be recorded as inert rather than logged as if they had an effect —
# the same contract as `vessel_mesh.LIVE_PARAMS`.
MASK_LIVE_PARAMS = {
    METHOD_OTSU: (),
    METHOD_NIBLACK: ("window_xy", "niblack_k"),
    METHOD_SAUVOLA: ("window_xy", "sauvola_k"),
    METHOD_BERNSEN: ("window_xy", "bernsen_contrast"),
}


def build_mask(vol: np.ndarray, voxel_zyx_um, method: str = METHOD_OTSU,
               window_xy: int = 50, niblack_k: float = 3.0,
               sauvola_k: float = 0.2, bernsen_contrast: float = 0.15,
               min_size: int | None = None):
    """Threshold `vol`, label components, and suggest a size cutoff.

    Returns ``(labels, sizes, suggested_min_size, mask)`` — always four, so a
    caller needing only three cannot drift from one needing four.

    `min_size=None` filters at `suggested_min_size`; pass a number to override.
    The returned `mask` is size-filtered, `labels`/`sizes` are not, so a caller
    can re-filter at a different cutoff without re-thresholding.

    **`method` is validated.** An unknown string raises rather than silently
    falling through to Otsu — [P30](../docs/ISSUES.md#p30) records that exact
    bug in `build_mesh` and there is no reason to reproduce it here.
    """
    if method not in MASK_METHODS:
        raise ValueError(
            f"unknown method {method!r}; expected one of {MASK_METHODS}")
    if method == METHOD_OTSU:
        mask = otsu_mask(vol)
    elif method == METHOD_NIBLACK:
        mask = niblack_mask(vol, int(window_xy), float(niblack_k), voxel_zyx_um)
    elif method == METHOD_SAUVOLA:
        mask = sauvola_mask(vol, int(window_xy), float(sauvola_k), voxel_zyx_um)
    else:
        mask = bernsen_mask(vol, int(window_xy), float(bernsen_contrast),
                            voxel_zyx_um)
    labels, sizes = label_and_sizes(mask)
    suggested = auto_min_size_otsu(sizes)
    cutoff = suggested if min_size is None else int(min_size)
    return labels, sizes, suggested, filter_by_size(labels, sizes, cutoff)
