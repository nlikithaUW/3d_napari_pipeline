"""Shared GUI helpers: histogram canvas, debounce decorator, strided pyramid,
the collapsible-group builder, the missing-file skip guard, and the
``after``-hand-off wrapper that keeps the "Run all" chain alive."""

from __future__ import annotations

import traceback
from functools import wraps
from pathlib import Path

import numpy as np
from magicgui.widgets import Container
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from qtpy.QtCore import QTimer
# superqt is a hard napari dependency (0.7.8 in both envs), so this adds
# nothing to install by hand in the Windows env.
from superqt import QCollapsible

from .. import gpu


def collapsible_group(title: str, widgets, *, extra=None,
                      expanded: bool = False) -> QCollapsible:
    """One pipeline step, one collapsible group.

    The layout rule inside a group is the Vessel Mesh tab's and now every
    tab's: inputs above the button, outputs and display-only controls below it,
    save/load last.

    Two things here are load-bearing and both were learned the hard way:

    * ``setDuration(0)`` — the collapse animation resizes badly inside the
      `QScrollArea` [P52] that `view.py::_wrap` puts around every tab.
    * ``coll._mg_container = inner`` — `addWidget` reparents the *native*
      QWidget, but nothing keeps the Python `Container` alive. When it is
      collected, the widgets underneath are deleted out from under the closures
      that still write to them, and the tab dies with "wrapped C/C++ object of
      type QLabel has been deleted" at the first status update.

    Lived as a nested function in `vessel_mesh.py` until a second tab needed
    it. Extracted rather than copied — the reference-holding line above is
    exactly the kind of non-obvious fix that a copy loses.
    """
    coll = QCollapsible(title)
    coll.setDuration(0)
    inner = Container(widgets=widgets, labels=True)
    coll.addWidget(inner.native)
    coll._mg_container = inner
    if extra is not None:
        coll.addWidget(extra)
    if expanded:
        coll.expand(animate=False)
    else:
        coll.collapse(animate=False)
    return coll


def call_all(*fns, label: str = "step") -> None:
    """Call each non-None callable in order; a failure is reported, not raised.

    Used for the "save the result, then continue the chain" hand-offs. If the
    save throws, the continuation must still run — otherwise the pipeline stops
    at that step with nothing to show for it.
    """
    for fn in fns:
        if fn is None:
            continue
        try:
            fn()
        except Exception:
            print(f"{label}: {getattr(fn, '__name__', 'callback')} failed —\n"
                  f"{traceback.format_exc()}", flush=True)


def skip_if_missing(path, after=None, label: str = "step",
                   what: str = "file", hint: str = "") -> bool:
    """True if ``path`` is absent — print a skip, hand off, and let the caller
    return without starting any work.

    A feature that was never saved is not an error, and "Load all saved" says
    so in its own docstring: a loader whose file is missing is meant to be
    skipped. Two reasons this has to happen *before* the load rather than in
    the failure handler:

    * A loader that reads inside a ``thread_worker`` cannot report a missing
      file quietly. ``create_worker`` connects its own ``reraise`` slot to
      ``errored`` (superqt ``_qthreading``) whenever the worker is built
      without a ``connect=`` dict — as every loader here is, since they connect
      their handlers afterwards — so an ``Errno 2`` arrives as a full traceback
      and a napari error popup even though the tab's own ``_err`` also runs and
      the chain advances. The user reads that as a dead pipeline.
    * Missing and broken must not print the same. Only a genuinely absent file
      short-circuits here; a truncated, unreadable or mismatched one still goes
      down the loader's error path untouched.

    Lived inline in `nuclei.py::_load_nuc`, the one loader of twelve that had
    it. Extracted rather than copied into the other eleven.
    """
    if Path(path).exists():
        return False
    msg = f"{label}: no {what} at {path} — skipped."
    print(f"{msg} {hint}" if hint else msg, flush=True)
    call_all(after, label=label)
    return True


def with_after(handler, after=None, label: str = "step"):
    """Wrap a worker completion handler so the ``after`` hand-off always fires.

    Qt swallows exceptions raised inside a signal slot. Without this, anything
    that throws part-way through ``handler`` skips the trailing ``after()`` and
    the "Run all" chain stops dead — no error, no timeout, just a pipeline that
    never reaches the next step. Runs ``handler``, prints any traceback it
    raises, then calls ``after`` from a ``finally``.

    Handlers wrapped here must NOT call ``after`` themselves.

    Every worker completion in the app goes through here, which makes it the one
    place that can drain the cupy memory pools once per *step* — see
    ``gpu.free_pools``. It runs before ``after`` so a chained next step, and any
    napari layer the handler just added, face a device with the room back.
    """
    @wraps(handler)
    def _wrapped(*args, **kwargs):
        try:
            return handler(*args, **kwargs)
        except Exception:
            print(f"{label}: handler failed —\n{traceback.format_exc()}",
                  flush=True)
        finally:
            gpu.free_pools()
            call_all(after, label=label)
    return _wrapped


def after_on_error(handler, after=None, label: str = "step"):
    """Like ``with_after``, but fires ``after`` *only* when ``handler`` raises.

    For a handler whose success path launches a further worker that owns the
    hand-off (e.g. auto-running the Clover threshold, then RF predict). Calling ``after``
    unconditionally there would advance the chain the moment the second worker
    *started*, not when it finished. Calling it never would stall the chain if
    the handler died before launching that worker — hence this middle ground.
    """
    @wraps(handler)
    def _wrapped(*args, **kwargs):
        try:
            return handler(*args, **kwargs)
        except Exception:
            print(f"{label}: handler failed —\n{traceback.format_exc()}",
                  flush=True)
            call_all(after, label=label)
    return _wrapped


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


def vector_mask(vectors, values=None):
    """Which (N, 2, 3) napari vectors are safe to draw — as a boolean mask.

    vispy builds each arrow's frame by normalizing the displacement
    (`vectors[:, 1]`), so a zero-length or non-finite vector divides by zero and
    crashes the GL draw.

    ``values`` is an optional per-vector colour feature, and passing it here is
    not a convenience: napari casts the mapped colour to uint8, so a **NaN
    feature value raises** ("cannot convert float NaN to integer") deep inside
    the colour swatch. Filtering the vectors and the colours in two separate
    steps is worse than not filtering at all — the survivors get each other's
    colours, silently. One mask, both arrays.
    """
    vectors = np.asarray(vectors)
    if len(vectors) == 0:
        return np.zeros(0, bool)
    ok = np.isfinite(vectors).all(axis=(1, 2))
    ok &= np.linalg.norm(vectors[:, 1], axis=1) > 1e-6
    if values is not None:
        ok &= np.isfinite(np.asarray(values, float))
    return ok


def finite_vectors(vectors):
    """Keep only the (N, 2, 3) napari vectors that are safe to draw.

    Thin wrapper on `vector_mask` — use that directly when the vectors carry a
    colour feature that has to be filtered with them.
    """
    vectors = np.asarray(vectors)
    if len(vectors) == 0:
        return vectors
    return vectors[vector_mask(vectors)]


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

    # Minimum canvas size, in px. Each tab stacks its controls and this canvas
    # in one QVBoxLayout (`view.py::_wrap`), so a tall control stack squeezes
    # the canvas — with no minimum it collapses to a few px and matplotlib
    # prints "constrained_layout not applied because axes sizes collapsed to
    # zero" on *every* redraw, which is once per slider move via `debounce`.
    # The floor is the decorations (title, both axis labels, two sets of ticks)
    # plus room for the bars, at the reduced font sizes set below.
    MIN_W, MIN_H = 180, 170

    def __init__(self, title: str, xlabel: str, color: str = "steelblue"):
        self._title = title
        self._xlabel = xlabel
        self._color = color
        self.fig = Figure(figsize=(4, 2.4), constrained_layout=True)
        self.ax = self.fig.add_subplot(111)
        # Smaller decorations for the same reason: this canvas lives in a dock
        # panel, not on a page, and the default 10 pt labels are most of the
        # height budget.
        self.ax.tick_params(labelsize=7)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setMinimumSize(self.MIN_W, self.MIN_H)
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
        ax.tick_params(labelsize=7)          # `clear` resets tick params
        ax.set_xlabel(self._xlabel, fontsize=8)
        ax.set_ylabel("count", fontsize=8)
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
        ax.set_title(self._title, fontsize=9)
        self.canvas.draw_idle()
