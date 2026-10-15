"""Cooperative stop: the Log panel's Stop button, checked by the long loops.

A Python thread cannot be killed safely, so a running step stops only where it
looks: `check()` raises `Cancelled` if a stop was requested. The checks sit in
the slow loops — RF predict (per z-chunk and per scale), RF training's feature
extraction (per scale), the adhesion split (per component), annotation seeding
(per stage) and ND2 conversion (per file) — and "Run all" / "Load all saved"
look between steps. A step that is stopped raises before it writes anything, so
no partial artifact is left on disk to be mistaken for a result.

One process-wide flag. It is set by `request()` and cleared by `clear()` at the
start of every action that can be stopped, so a Stop pressed while nothing is
running cannot cancel the next thing you start. No Qt; headless runners never
set it (they stop with Ctrl+C).
"""

from __future__ import annotations

import threading

_EVENT = threading.Event()


class Cancelled(Exception):
    """Raised by `check` when the user pressed Stop."""


#: What a GUI worker returns instead of raising `Cancelled`: an exception
#: escaping a napari `thread_worker` is re-raised by superqt as a traceback and
#: an error popup (P62), so a deliberate Stop would look like a crash. Workers
#: catch `Cancelled`, return this, and their `_done` reports the stop.
STOPPED = object()


def request() -> None:
    _EVENT.set()


def clear() -> None:
    _EVENT.clear()


def requested() -> bool:
    return _EVENT.is_set()


def check(where: str = "") -> None:
    if _EVENT.is_set():
        raise Cancelled(f"stopped by user{f' during {where}' if where else ''}")
