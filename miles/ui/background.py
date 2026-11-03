"""A panel's slow work on the application's worker threads, reported through signals."""

import weakref
from collections.abc import Callable
from typing import Any

from PySide6.QtCore import QObject, Signal

from miles.ui.tasks import executor


class Background(QObject):
    """Runs a function on a worker thread and reports back through signals.

    The function must not touch Qt widgets; its result arrives on the
    interface thread through :attr:`finished`, errors through :attr:`failed`.
    A new :meth:`run` replaces a request still running, whose result is then
    dropped, as is any result that arrives after the object was deleted (its
    panel was closed, the application is quitting). The work runs on the
    application's executor (:mod:`miles.ui.tasks`), which holds this object
    only weakly.
    """

    finished = Signal(object)
    failed = Signal(str)

    def run(self, func: Callable[..., Any], *args: Any) -> None:
        """Call ``func(*args)`` on a worker thread."""
        me = weakref.ref(self)

        def done(result: Any) -> None:
            owner = me()
            if owner is not None:
                owner.finished.emit(result)

        def failed(message: str) -> None:
            owner = me()
            if owner is not None:
                owner.failed.emit(message)

        executor().submit(func, *args, done=done, failed=failed, owner=self, key=("background", id(self)))
