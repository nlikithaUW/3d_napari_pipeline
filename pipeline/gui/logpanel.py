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
    signal for the GUI to display."""

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

    def flush(self):
        if self._original is not None:
            try:
                self._original.flush()
            except Exception:
                pass


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
        top = QHBoxLayout()
        top.addWidget(clear)
        top.addStretch(1)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.addLayout(top)
        lay.addWidget(self._text)

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
