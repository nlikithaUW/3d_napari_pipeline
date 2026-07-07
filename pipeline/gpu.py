"""Optional GPU shim for windowed ndimage ops, with a clean CPU fallback.

Probes for cupy once at import. If a working CUDA device is present, the
wrappers below run on the GPU (``cp.asarray`` -> ``cupyx.scipy.ndimage`` ->
``cp.asnumpy``); otherwise they transparently call ``scipy.ndimage``. An
out-of-memory error on a single call degrades to CPU for that call rather than
crashing, so a too-large volume is handled gracefully.

cupy is entirely optional: a GPU-less env (the import fails) just uses scipy.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as _sndi

try:
    import cupy as cp
    import cupyx.scipy.ndimage as cndi

    cp.zeros(1)            # force context init; any failure -> disable
    HAS_GPU = True
except Exception:
    cp = cndi = None
    HAS_GPU = False


def gpu_enabled() -> bool:
    return HAS_GPU


def _oom():
    """The cupy OOM exception type, or a never-raised sentinel on CPU."""
    if HAS_GPU:
        return cp.cuda.memory.OutOfMemoryError
    return ()


def _run(gpu_op, cpu_op, arr: np.ndarray, **kw) -> np.ndarray:
    """Run ``gpu_op`` on the GPU if available, else ``cpu_op`` on the CPU.

    Falls back to CPU for this single call on out-of-memory.
    """
    if HAS_GPU:
        try:
            return cp.asnumpy(gpu_op(cp.asarray(arr), **kw))
        except _oom():
            if cp is not None:
                cp.get_default_memory_pool().free_all_blocks()
    return cpu_op(arr, **kw)


def _take(arr, axis: int, start: int, stop: int):
    idx = [slice(None)] * arr.ndim
    idx[axis] = slice(start, stop)
    return arr[tuple(idx)]


def _box_mean_gpu(a, sizes):
    """Separable O(1) box mean on the GPU via per-axis prefix sums (cumsum).

    Independent of window size: a wide window costs the same as a narrow one,
    unlike cupy's dense ``uniform_filter`` (O(window) per voxel). Accumulates
    in float32 — fast on the GPU (float64 is ~64× slower on a 3090), and
    precise because the cumsum runs one short axis at a time (length ≤ a few
    thousand, not the whole-volume integral image that loses precision) over
    mean-centred input. ``symmetric`` padding reproduces scipy's
    ``mode="reflect"`` boundary exactly (odd windows only, which is all we use).
    """
    out = a.astype(cp.float32, copy=True)
    for axis, w in enumerate(sizes):
        w = int(w)
        if w <= 1:
            continue
        r = w // 2  # w is odd, so w == 2r + 1
        pad = [(0, 0)] * out.ndim
        pad[axis] = (r, r)
        P = cp.pad(out, pad, mode="symmetric")
        S = cp.cumsum(P, axis=axis)
        zshape = list(S.shape)
        zshape[axis] = 1
        S0 = cp.concatenate([cp.zeros(zshape, dtype=S.dtype), S], axis=axis)
        n = out.shape[axis]
        out = (_take(S0, axis, w, w + n) - _take(S0, axis, 0, n)) / w
    return out


def uniform_filter(arr: np.ndarray, size, **kw) -> np.ndarray:
    """Box-mean filter. GPU path is window-size-independent (prefix sums);
    CPU path is scipy's already-O(1) ``uniform_filter``. Drop-in for the
    subset of kwargs we use (``size``, default ``mode="reflect"``).
    """
    if HAS_GPU:
        try:
            a = cp.asarray(arr)
            sizes = tuple(size) if isinstance(size, (tuple, list)) else (int(size),) * a.ndim
            out = _box_mean_gpu(a, sizes)
            return cp.asnumpy(out).astype(arr.dtype, copy=False)
        except _oom():
            cp.get_default_memory_pool().free_all_blocks()
    return _sndi.uniform_filter(arr, size=size, **kw)


def local_mean_std(v: np.ndarray, win) -> tuple[np.ndarray, np.ndarray]:
    """Local mean and std over the ``win`` footprint, as float32 host arrays.

    The whole computation (mean-centre, both box means, ``m2 − m²``, sqrt) runs
    on the device when a GPU is available, so the volume crosses the PCIe bus
    only once — transferred as uint16 (half the bytes of float32) and cast on
    the device. This matters far more than the filter's flop count: the windowed
    math is cheap, host<->device copies are not. Mean-centres before the
    variance step to keep float32 precise. CPU fallback uses scipy's O(1)
    ``uniform_filter``.
    """
    sizes = tuple(win) if isinstance(win, (tuple, list)) else (int(win),) * v.ndim
    if HAS_GPU:
        try:
            d = cp.asarray(v).astype(cp.float32)   # single H2D (uint16 if v is)
            mu = d.mean()                           # 0-d cupy float32 scalar
            d -= mu
            m = _box_mean_gpu(d, sizes)
            d *= d                                  # vs² in place; d is free after
            m2 = _box_mean_gpu(d, sizes)
            del d
            var = m2 - m * m
            cp.maximum(var, 0.0, out=var)
            return cp.asnumpy(m + mu), cp.asnumpy(cp.sqrt(var))
        except _oom():
            cp.get_default_memory_pool().free_all_blocks()
    vf = v.astype(np.float32, copy=False)
    mu = float(vf.mean())
    vs = vf - mu
    m = _sndi.uniform_filter(vs, size=sizes)
    m2 = _sndi.uniform_filter(vs * vs, size=sizes)
    var = m2 - m * m
    np.maximum(var, 0.0, out=var)
    return m + mu, np.sqrt(var)


def median_filter(arr: np.ndarray, **kw) -> np.ndarray:
    return _run(cndi.median_filter if HAS_GPU else None,
               _sndi.median_filter, arr, **kw)


def minimum_filter(arr: np.ndarray, **kw) -> np.ndarray:
    return _run(cndi.minimum_filter if HAS_GPU else None,
               _sndi.minimum_filter, arr, **kw)


def maximum_filter(arr: np.ndarray, **kw) -> np.ndarray:
    return _run(cndi.maximum_filter if HAS_GPU else None,
               _sndi.maximum_filter, arr, **kw)


def zoom(arr: np.ndarray, *args, **kw) -> np.ndarray:
    gpu_op = (lambda a, **k: cndi.zoom(a, *args, **k)) if HAS_GPU else None
    cpu_op = lambda a, **k: _sndi.zoom(a, *args, **k)
    return _run(gpu_op, cpu_op, arr, **kw)


def gaussian_filter(arr: np.ndarray, sigma, order=0, **kw) -> np.ndarray:
    gpu_op = (lambda a, **k: cndi.gaussian_filter(a, sigma, order=order, **k)) if HAS_GPU else None
    cpu_op = lambda a, **k: _sndi.gaussian_filter(a, sigma, order=order, **k)
    return _run(gpu_op, cpu_op, arr, **kw)
