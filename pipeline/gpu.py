"""Optional GPU shim for windowed ndimage ops, with a clean CPU fallback.

Probes for cupy once at import. If a working CUDA device is present, the
wrappers below run on the GPU (``cp.asarray`` -> ``cupyx.scipy.ndimage`` ->
``cp.asnumpy``); otherwise they transparently call ``scipy.ndimage``. An
out-of-memory error on a single call degrades to CPU for that call (with a
printed note) rather than crashing, so a too-large volume is handled gracefully.
Any *other* CUDA failure propagates — see ``_oom`` / ``_is_oom``.

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


# CUDA runtime status codes (driver_types.h). cupy raises a single
# ``CUDARuntimeError`` carrying the code in ``.status`` but does not re-export
# the names (checked against cupy 14.1), so they are spelled out here.
cudaErrorMemoryAllocation = 2
cudaErrorAlreadyMapped = 208


def _oom():
    """cupy exception types that *may* mean "the device ran out of room" — the
    only condition worth retrying on the CPU.

    Two classes qualify. ``OutOfMemoryError`` is the memory-pool OOM and always
    counts. ``CUDARuntimeError`` is the base class for *every* CUDA runtime
    failure — an illegal memory access and a kernel launch failure raise it just
    as a pinned-host-alloc OOM does — so catching it alone would silently demote
    real bugs to a slow CPU path. ``_is_oom`` filters it down by status code.

    Empty (never-raised) on a CPU-only env.
    """
    if HAS_GPU:
        return (cp.cuda.memory.OutOfMemoryError,
                cp.cuda.runtime.CUDARuntimeError)
    return ()


def _is_oom(exc) -> bool:
    """True only for genuine resource exhaustion, so everything else re-raises.

    ``cudaErrorMemoryAllocation`` is device/pinned-host alloc failure.
    ``cudaErrorAlreadyMapped`` is the pinned-transfer pool being exhausted by
    concurrent host threads (see ``fret.resample_fret_outputs``). Both are
    retryable on the CPU; an illegal access or launch failure is not.
    """
    if isinstance(exc, cp.cuda.memory.OutOfMemoryError):
        return True
    status = getattr(exc, "status", None)
    return status in (cudaErrorMemoryAllocation, cudaErrorAlreadyMapped)


def free_pools() -> None:
    """Drain the cupy memory pools, tolerating a sticky CUDA error.

    Two callers, for two different reasons.

    ``_fallback`` calls it from inside an exception handler. After a sticky
    CUDA error every subsequent CUDA call — including these — raises, so an
    unguarded drain replaces the original exception with a traceback pointing
    at the cleanup. Hence the bare ``except``.

    The GUI calls it at the end of every step (``gui/widgets.py::with_after``).
    cupy's default pool keeps every block it has ever allocated for the life of
    the process, so a step that filtered a volume holds that VRAM forever —
    against a napari canvas that has to allocate GL textures for the *next*
    image out of the same device. A failed GL allocation is not reported to
    Python; it surfaces as an ``access violation`` in a later ``glDrawArrays``
    ([P48](docs/ISSUES.md#p48)). Draining is per *step*, never per op: within a
    step the pool's reuse is the point of having one.
    """
    try:
        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass


def device_memory_mb():
    """(free, total) device memory in MB, or ``None`` without a GPU.

    Reported at image load so that a later GL fault has a number next to it
    rather than a guess — see ``free_pools``.
    """
    if not HAS_GPU:
        return None
    try:
        free_b, total_b = cp.cuda.runtime.memGetInfo()
    except Exception:
        return None
    return free_b / 1e6, total_b / 1e6


def _fallback(exc, what: str) -> None:
    """Re-raise unless ``exc`` is a retryable OOM; otherwise drain and report.

    The report matters: a silent fallback turns "the GPU is out of memory" into
    "why did this run take twenty minutes", with nothing on screen to explain it.
    """
    if not _is_oom(exc):
        raise exc
    free_pools()
    print(f"  gpu: {what} fell back to CPU ({type(exc).__name__}: {exc})",
          flush=True)


def _run(gpu_op, cpu_op, arr: np.ndarray, name: str = "", **kw) -> np.ndarray:
    """Run ``gpu_op`` on the GPU if available, else ``cpu_op`` on the CPU.

    Falls back to CPU for this single call on out-of-memory; any other CUDA
    failure propagates. ``name`` labels the fallback message — pass it when
    ``cpu_op`` is a lambda and so has no useful ``__name__``.
    """
    if HAS_GPU:
        try:
            return cp.asnumpy(gpu_op(cp.asarray(arr), **kw))
        except _oom() as exc:
            _fallback(exc, name or getattr(cpu_op, "__name__", "op"))
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
        except _oom() as exc:
            _fallback(exc, "uniform_filter")
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
        except _oom() as exc:
            _fallback(exc, "local_mean_std")
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
    return _run(gpu_op, cpu_op, arr, name="zoom", **kw)


def gaussian_filter(arr: np.ndarray, sigma, order=0, **kw) -> np.ndarray:
    gpu_op = (lambda a, **k: cndi.gaussian_filter(a, sigma, order=order, **k)) if HAS_GPU else None
    cpu_op = lambda a, **k: _sndi.gaussian_filter(a, sigma, order=order, **k)
    return _run(gpu_op, cpu_op, arr, name="gaussian_filter", **kw)
