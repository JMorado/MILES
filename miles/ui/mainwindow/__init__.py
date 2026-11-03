"""The main window of MILES.

:class:`MainWindow` puts its parts together, one module each: the tabs
(:mod:`~miles.ui.mainwindow.tabs`), files
(:mod:`~miles.ui.mainwindow.files`) and dialogs
(:mod:`~miles.ui.mainwindow.dialogs`); see
:class:`~miles.ui.mainwindow.window.MainWindow`.
"""

from miles.ui.mainwindow.files import RECENT_LIMIT
from miles.ui.mainwindow.window import APP_NAME, MainWindow

__all__ = ["APP_NAME", "RECENT_LIMIT", "MainWindow"]
