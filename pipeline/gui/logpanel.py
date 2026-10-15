"""In-napari log panel.

Mirrors ``stdout``/``stderr`` into a dock widget so all the pipeline's
``print(...)`` progress messages are visible inside napari, not just in the
terminal that launched it. The original streams are still written to (teed), so
the terminal keeps its output too.

Writes come from worker threads (``@thread_worker``), so the stream emits a Qt
signal and the GUI slot does the actual text insertion — Qt queues the
cross-thread delivery, which keeps all widget access on the GUI thread.
"""

from __future__ import annotations

import sys

from qtpy.QtCore import QObject, Signal, Slot

from .. import cancel
from qtpy.QtGui import QTextCursor
from qtpy.QtWidgets import (
    QHBoxLayout,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


class _EmitStream(QObject):
    """A text stream that tees writes to ``original`` and emits them as a
    signal for the GUI to display.

    Implements enough of the file protocol that libraries probing the stream
    behave. ``isatty``/``fileno``/``encoding`` in particular: progress bars
    (tqdm, joblib) and anything calling ``os.get_terminal_size`` check those
    first, and a stream missing them either raises or makes the bar decide it is
    writing to a pipe and dump one line per update into the panel.
    """

    text = Signal(str)

    def __init__(self, original):
        super().__init__()
        self._original = original

    def write(self, s):
        if self._original is not None:
            try:
                self._original.write(s)
            except Exception:
                pass
        if s:
            # Safe to emit from any thread; delivery to the GUI slot is queued.
            self.text.emit(s)
        return len(s) if s else 0

    def writelines(self, lines):
        for line in lines:
            self.write(line)

    def flush(self):
        if self._original is not None:
            try:
                self._original.flush()
            except Exception:
                pass

    def isatty(self) -> bool:
        # False even when the launching terminal is a tty: the panel is the
        # primary destination and it cannot interpret carriage returns or ANSI
        # cursor moves, so callers should use their non-interactive output.
        return False

    def fileno(self):
        # Delegate so C-level writers can still find a real descriptor; raise
        # the standard error (not AttributeError) when there isn't one, which is
        # what callers are written to catch.
        if self._original is not None:
            try:
                return self._original.fileno()
            except Exception:
                pass
        raise OSError("no file descriptor behind the log panel stream")

    @property
    def encoding(self) -> str:
        return getattr(self._original, "encoding", None) or "utf-8"

    @property
    def errors(self):
        return getattr(self._original, "errors", None)

    def writable(self) -> bool:
        return True

    def readable(self) -> bool:
        return False

    def seekable(self) -> bool:
        return False

    @property
    def closed(self) -> bool:
        return False

    def close(self) -> None:
        # Never close the real stdout/stderr underneath us.
        self.flush()


class LogConsole(QWidget):
    """Read-only text panel with a Clear button; appends incoming log text."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._text = QPlainTextEdit()
        self._text.setReadOnly(True)
        self._text.setMaximumBlockCount(10000)  # cap memory on long sessions
        self._text.setLineWrapMode(QPlainTextEdit.NoWrap)

        clear = QPushButton("Clear")
        clear.clicked.connect(self._text.clear)
        # Stop lives here, not on a tab, so it is reachable whatever tab is
        # showing. Cooperative (pipeline/cancel.py): "Run all" stops between
        # steps, and the long loops (RF predict/train, splitting, seeding,
        # Convert) stop at their next check without writing a result.
        stop = QPushButton("■ Stop")
        stop.setToolTip("Stop the running step at its next check point and "
                        "don't start the next step of Run all / Load all.")
        stop.clicked.connect(self._stop)
        top = QHBoxLayout()
        top.addWidget(stop)
        top.addWidget(clear)
        top.addStretch(1)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.addLayout(top)
        lay.addWidget(self._text)

    @Slot()
    def _stop(self) -> None:
        cancel.request()
        print("■ Stop requested — the running step stops at its next check "
              "point (RF predict: within ~one z-chunk); Run all will not start "
              "another step.", flush=True)

    @Slot(str)
    def append(self, s: str) -> None:
        self._text.moveCursor(QTextCursor.End)
        self._text.insertPlainText(s)
        self._text.moveCursor(QTextCursor.End)


def install(viewer) -> LogConsole:
    """Redirect stdout/stderr into a new bottom 'Log' dock panel.

    Returns the console widget. The original streams are preserved (teed), so
    the launching terminal keeps showing output.
    """
    console = LogConsole()

    out = _EmitStream(sys.stdout)
    err = _EmitStream(sys.stderr)
    # Parent the streams to the console so they stay alive as long as it does.
    out.setParent(console)
    err.setParent(console)
    out.text.connect(console.append)
    err.text.connect(console.append)
    sys.stdout = out
    sys.stderr = err

    viewer.window.add_dock_widget(console, area="bottom", name="Log")
    print("Log panel ready — pipeline progress will appear here.", flush=True)
    return console


def uninstall() -> None:
    """Restore the original ``stdout``/``stderr``.

    Without this the redirect outlives the viewer: closing napari leaves the
    interpreter (or a subsequent headless run in the same process) writing into
    a deleted Qt widget. Only unwinds streams this module installed, so calling
    it twice, or when ``install`` never ran, is harmless.
    """
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name, None)
        if isinstance(stream, _EmitStream):
            setattr(sys, name, stream._original)
