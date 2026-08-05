"""Shared GUI helpers: histogram canvas, debounce decorator, strided pyramid."""

from __future__ import annotations

from functools import wraps

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from qtpy.QtCore import QTimer


def pyramid_views(vol: np.ndarray, levels: int):
    """Return a strided-view multiscale pyramid. No copies."""
    if levels <= 1:
        return vol
    out = [vol]
    step = 2
    for _ in range(levels - 1):
        out.append(vol[::step, ::step, ::step])
        step *= 2
    return out


def finite_surface(verts, faces, values):
    """Sanitize a (verts, faces, values) surface so vispy never faults the GL
    draw. Three failure modes each cause an ``OSError: access violation`` inside
    ``glDrawArrays`` because vispy reads past a GL buffer:

    1. A face index outside the vertex array (mismatched verts/faces — e.g. a
       stale cache, or faces built against a different vertex set).
    2. A non-finite vertex coordinate.
    3. A per-vertex ``values`` array that is non-finite or the wrong length.

    Guards all three: drops out-of-range and non-finite-touching triangles,
    replaces stray vertices with the finite centroid, and makes ``values``
    finite and length-consistent. Returns a (verts, faces, values) tuple.
    """
    verts = np.asarray(verts)
    n = len(verts)
    faces = np.asarray(faces)

    # 1) Drop faces whose indices fall outside the vertex array.
    if len(faces):
        inrange = (faces >= 0).all(axis=1) & (faces < n).all(axis=1)
        if not inrange.all():
            faces = faces[inrange]

    # 2) Replace non-finite vertex coords; drop the triangles that touch them.
    bad = ~np.isfinite(verts).all(axis=1)
    if bad.any():
        if len(faces):
            faces = faces[~bad[faces].any(axis=1)]
        verts = verts.copy()
        verts[bad] = verts[~bad].mean(0) if (~bad).any() else 0.0

    # 3) Make per-vertex values finite and length-consistent (else vispy reads
    #    past the value buffer). Fill non-finite entries with the finite median.
    if values is not None:
        values = np.asarray(values)
        if values.shape[0] != n:
            values = np.zeros(n, dtype=np.float32)
        elif not np.isfinite(values).all():
            fin = np.isfinite(values)
            fill = float(np.median(values[fin])) if fin.any() else 0.0
            values = np.where(fin, values, fill)

    return verts, faces, values


def finite_vectors(vectors):
    """Keep only (N, 2, 3) napari vectors with finite, non-zero-length

    displacement. vispy builds each arrow's frame by normalizing the
    displacement (`vectors[:, 1]`), so a zero-length or non-finite vector
    divides by zero and crashes the GL draw.
    """
    vectors = np.asarray(vectors)
    if len(vectors) == 0:
        return vectors
    finite = np.isfinite(vectors).all(axis=(1, 2))
    length = np.linalg.norm(vectors[:, 1], axis=1)
    return vectors[finite & (length > 1e-6)]


def percentile_limits(vol, lo_pct: float = 1.0, hi_pct: float = 99.9,
                      subsample: int = 4):
    """Robust contrast window (lo, hi) from intensity percentiles.

    napari's auto-contrast takes min/max of a sparse 3-plane sample, which
    collapses onto the noise floor for dark-background / sparse-bright data
    (thin cadherin arcs, nuclei) and magnifies streaks. Percentiles are
    stable: `hi_pct` ignores the rare bright outliers, `lo_pct` clips the
    floor. Subsamples for speed; pass a multiscale list and level 0 is used.
    """
    # A multiscale pyramid arrives as a list/tuple or napari `MultiScaleData`
    # (a Sequence, not a list); both index level 0 = full res. Plain ndarrays
    # are used as-is — `np.asarray` on MultiScaleData would pick the coarsest.
    if not isinstance(vol, np.ndarray):
        vol = vol[0]
    vol = np.asarray(vol)
    s = vol[::subsample, ::subsample, ::subsample] if vol.ndim == 3 else vol
    lo = float(np.percentile(s, lo_pct))
    hi = float(np.percentile(s, hi_pct))
    if hi <= lo:
        hi = lo + 1.0
    return lo, hi


def debounce(ms: int):
    """Coalesce rapid magicgui auto_call invocations into one call after `ms`."""
    def deco(fn):
        timer = QTimer()
        timer.setSingleShot(True)
        pending = {}

        def _fire():
            args, kwargs = pending.pop("args", ((), {}))
            fn(*args, **kwargs)

        timer.timeout.connect(_fire)

        @wraps(fn)
        def wrapper(*args, **kwargs):
            pending["args"] = (args, kwargs)
            timer.start(ms)

        wrapper._timer = timer  # keep reference alive
        return wrapper
    return deco


class Histogram:
    """Small matplotlib canvas with one redraw method for log10 size hists."""

    def __init__(self, title: str, xlabel: str, color: str = "steelblue"):
        self._title = title
        self._xlabel = xlabel
        self._color = color
        self.fig = Figure(figsize=(4, 2.4), constrained_layout=True)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvas(self.fig)
        self.clear()

    def clear(self) -> None:
        self.update(np.array([]))

    def update(
        self,
        values: np.ndarray,
        *,
        threshold: float | None = None,
        threshold_label: str | None = None,
        log: bool = True,
    ) -> None:
        ax = self.ax
        ax.clear()
        ax.set_xlabel(self._xlabel)
        ax.set_ylabel("count")
        v = np.asarray(values) if values is not None else np.array([])
        if log:
            v = v[np.isfinite(v) & (v > 0)]
            v = np.log10(v) if v.size else v
        if v.size:
            ax.hist(v, bins=60, color=self._color)
            if threshold is not None and threshold > 0:
                x = np.log10(threshold) if log else threshold
                ax.axvline(x, color="crimson", lw=2, label=threshold_label or "")
                if threshold_label:
                    ax.legend(loc="upper right", fontsize=8)
        ax.set_title(self._title)
        self.canvas.draw_idle()
