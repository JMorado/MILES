"""Showing unexpected errors in a dialog, so that a bug does not silently break the window."""

import sys
import traceback
from collections.abc import Callable
from types import TracebackType

from PySide6.QtWidgets import QApplication, QMessageBox, QWidget


class ErrorReporter:
    """Exception hook that prints the traceback and shows it in a dialog.

    Qt calls Python code from its event loop; an exception there (in a slot,
    a paint event, a timer) would otherwise only be printed to the terminal,
    which most users never see, while the window may be left half updated.
    Installed as ``sys.excepthook``, the reporter shows what happened and
    keeps the application running, so the user can save their work.

    Parameters
    ----------
    window : callable
        Returns the window to centre the dialog on (or None), looked up when
        an error happens.
    """

    def __init__(self, window: Callable[[], QWidget | None]) -> None:
        self.window = window
        self.shown = 0   # dialogs shown so far, for tests
        self._showing = False
        self._seen: set[tuple[str, str, str]] = set()   # errors already shown, so a repeating one shows once

    def __call__(self, kind: type[BaseException], value: BaseException, tb: TracebackType | None) -> None:
        """Print the error and show it, unless a dialog is open or the same error was shown before."""
        details = "".join(traceback.format_exception(kind, value, tb))
        print(details, file=sys.stderr, end="")
        frames = traceback.extract_tb(tb)
        where = f"{frames[-1].filename}:{frames[-1].lineno}" if frames else ""
        signature = (kind.__name__, str(value), where)
        if (issubclass(kind, KeyboardInterrupt) or self._showing or signature in self._seen
                or QApplication.instance() is None):
            return   # errors raised while a dialog is open, or again (each repaint), are only printed
        self._seen.add(signature)
        self._showing = True
        try:
            box = QMessageBox(self.window())
            box.setIcon(QMessageBox.Icon.Warning)
            box.setWindowTitle("Something went wrong")
            box.setText(f"{kind.__name__}: {value}")
            box.setInformativeText("The window is still usable, but save your work to be safe. "
                                   "The details below help to find the cause.")
            box.setDetailedText(details)
            self.shown += 1
            box.exec()
        finally:
            self._showing = False


def install(window: Callable[[], QWidget | None]) -> ErrorReporter:
    """Show unexpected errors in a dialog from now on; returns the hook installed."""
    reporter = ErrorReporter(window)
    sys.excepthook = reporter
    return reporter
