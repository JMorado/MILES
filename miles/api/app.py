"""The application as scripts and plugins see it: :class:`App`, and the running one."""

import os
from collections.abc import Callable
from enum import StrEnum
from typing import TYPE_CHECKING, Any, TypeVar

from miles.core.molecule import Molecule
from miles.files.common import JobInfo

if TYPE_CHECKING:
    from PySide6.QtWidgets import QWidget

    from miles.plugins import PluginInfo
    from miles.render.theme import Theme
    from miles.ui.mainwindow import MainWindow
    from miles.ui.plot_export import PlotSource
    from miles.ui.tasks import Task

from miles.api.document import DocumentView, document_view

T = TypeVar("T")

_current: "App | None" = None


class Severity(StrEnum):
    """How much a message matters (see :meth:`App.message`)."""

    INFO = "info"         # a short note at the top of the tab
    WARNING = "warning"   # a note that stays longer
    ERROR = "error"       # a dialog the user has to close

SECONDS = {Severity.INFO: 2.4, Severity.WARNING: 4.0}
"""How long notes of each severity stay, unless the caller says."""

def current_app() -> "App":
    """Return the running application.

    Raises
    ------
    RuntimeError
        If no MILES window has been created.
    """
    if _current is None:
        raise RuntimeError("MILES is not running")
    return _current


class App:
    """The MILES window, as seen from Python.

    Parameters
    ----------
    window : MainWindow
        The main window to drive.
    """

    def __init__(self, window: "MainWindow") -> None:
        global _current
        self._window = window
        _current = self

    @property
    def window(self) -> "MainWindow":
        """The Qt main window, for experiments in the console.

        It is internal: it may change in any release, so plugins use the
        methods of this class instead (and :attr:`widget` to parent dialogs).
        """
        return self._window

    @property
    def widget(self) -> "QWidget":
        """The main window as a widget, to parent a plugin's dialogs, menus and other Qt objects.

        Everything else goes through the methods of this class and of
        :class:`DocumentView`.
        """
        return self._window

    @property
    def documents(self) -> list["DocumentView"]:
        """Every open structure, in tab order."""
        return [document_view(self._window, canvas) for canvas in self._window.canvases]

    @property
    def current(self) -> "DocumentView | None":
        """The structure in the current tab, if any."""
        canvas = self._window.current_canvas()
        return document_view(self._window, canvas) if canvas is not None else None

    def open(self, path: str) -> "DocumentView":
        """Open a file in a new tab.

        Parameters
        ----------
        path : str
            Any file MILES reads (structures, trajectories, cube and
            Molden files, sessions).

        Returns
        -------
        DocumentView
            The new tab.

        Raises
        ------
        FileNotFoundError
            If there is no such file.
        FileFormatError
            If the file cannot be read (its format is unknown or it is
            damaged); unlike **File > Open**, no dialog is shown.
        """
        from miles.files.formats import read_file

        if not os.path.exists(path):
            raise FileNotFoundError(path)
        mol, info = read_file(path)
        self._window.remember_dir(path)
        self._window.open_structure(mol, info, path)
        return self._current_or_raise()

    def new(self, molecule: Molecule | None = None, title: str = "") -> "DocumentView":
        """Show a molecule in a new tab (an empty one without a molecule).

        Parameters
        ----------
        molecule : Molecule, optional
            Structure to show, e.g. from :func:`miles.from_smiles <miles.core.chem.from_smiles>`.
        title : str, optional
            Name of the tab.

        Returns
        -------
        DocumentView
            The new tab.
        """
        if molecule is None:
            self._window.new_document()
        else:
            self._window.open_structure(molecule, JobInfo(title=title or "molecule"))
        return self._current_or_raise()

    def commands(self) -> list[str]:
        """Every menu command, as ``"Menu > Item"`` (the names :meth:`run` takes)."""
        from miles.ui.command_palette import menu_entries

        return [entry.path for entry in menu_entries(self._window.menuBar())]

    def run(self, name: str) -> None:
        """Run a menu command, as if it were chosen from the menu.

        Parameters
        ----------
        name : str
            The command's full path (``"Edit > Clean Geometry"``) or, when
            that is unique, just its title (``"Clean Geometry"``); case and a
            trailing ``"..."`` do not matter.

        Raises
        ------
        KeyError
            If no command or more than one has that name.
        RuntimeError
            If the command is disabled at the moment.
        """
        from miles.ui.command_palette import menu_entries

        wanted = name.lower().removesuffix("...").strip()
        entries = menu_entries(self._window.menuBar())
        found = [e for e in entries if e.path.lower() == wanted]
        if not found:
            found = [e for e in entries if e.path.lower().rsplit(" > ", 1)[-1] == wanted]
        if not found:
            raise KeyError(f"No command named {name!r}")
        if len(found) > 1:
            raise KeyError(f"{name!r} could be: " + ", ".join(e.path for e in found))
        manager = getattr(self._window, "plugins", None)
        action = found[0].action
        if not action.isEnabled() or (manager is not None and not manager.available(action)):
            raise RuntimeError(f"{found[0].path} is not available now")
        action.trigger()

    def message(self, text: str, seconds: float | None = None, severity: str = Severity.INFO,
                title: str = "") -> None:
        """Tell the user something: a note at the top of the current tab, or a dialog for an error.

        Parameters
        ----------
        text : str
            What to say.
        seconds : float, optional
            How long a note stays (default: by severity).
        severity : str, optional
            A :class:`Severity`: ``"info"`` (the default) and ``"warning"``
            show a note, ``"warning"`` for longer; ``"error"`` shows a dialog
            that the user closes, e.g. when a file cannot be read.
        title : str, optional
            Title of the dialog of an error (default: "MILES").
        """
        level = Severity(severity)
        if level == Severity.ERROR:
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(self._window, title or "MILES", text)
            return
        canvas = self._window.current_canvas()
        if canvas is not None:
            canvas.toast(text, int((seconds if seconds is not None else SECONDS[level]) * 1000))

    def warn(self, title: str, text: str) -> None:
        """Show an error in a dialog, e.g. when a file cannot be written (an ``"error"`` :meth:`message`)."""
        self.message(text, severity=Severity.ERROR, title=title)

    def export_plot(self, source: "PlotSource", kind: str) -> None:
        """Copy a plot's data, or ask where to save it as CSV or its picture as PNG or SVG.

        This is what the right-click menus of the built-in plots do.

        Parameters
        ----------
        source : PlotSource
            The panel or widget with the plot
            (:class:`miles.plugin.widgets.PlotSource <miles.ui.plot_export.PlotSource>`).
        kind : str
            :data:`~miles.ui.plot_export.COPY`, :data:`~miles.ui.plot_export.DATA`,
            :data:`~miles.ui.plot_export.CURVE` or :data:`~miles.ui.plot_export.IMAGE`.
        """
        self._window.export_plot(source, kind)

    @property
    def theme(self) -> "Theme":
        """Colours of the interface, for plugin dialogs and panels (a :class:`~miles.render.theme.Theme`)."""
        return self._window.theme

    @property
    def plugins(self) -> list["PluginInfo"]:
        """The plugins found at start-up, with their state and errors."""
        return list(self._window.plugins.plugins)

    def run_in_background(self, job: Callable[..., Any], *args: Any, message: str = "Working",
                          done: Callable[[Any], Any] | None = None,
                          failed: Callable[[str], Any] | None = None, cancellable: bool = False) -> "Task":
        """Run slow work on a worker thread while the window stays responsive.

        A bar at the top right says ``message`` with a Cancel link until the
        work ends. Give ``job`` copies of what it needs (e.g.
        ``doc.molecule.copy()``): it must not touch the window or the
        documents. Apply the result in ``done``, which runs on the interface
        thread; check there that the document still is as it was
        (``doc.revision``).

        Parameters
        ----------
        job : callable
            Called as ``job(*args)`` on the worker thread.
        *args
            Its arguments.
        message : str, optional
            What the bar says, e.g. "Searching for conformers".
        done : callable, optional
            ``done(result)``, unless the task was cancelled.
        failed : callable, optional
            ``failed(message)`` when ``job`` raised (default: a message at
            the top of the canvas).
        cancellable : bool, optional
            Call ``job(*args, cancel=token)``: Cancel sets ``token.cancelled``
            (a :class:`~miles.ui.tasks.CancelToken`), so a job that checks
            it between steps stops early. Otherwise Cancel only drops the
            result and the job runs to its end.

        Returns
        -------
        Task
            The running task (:mod:`miles.ui.tasks`); ``cancel()``
            drops its result and sets the token.
        """
        def report(text: str) -> None:
            self.message(f"{message} failed: {text}", 3.2)

        return self._window.stage.tasks.start(message, job, *args, done=done,
                                              failed=failed if failed is not None else report,
                                              pass_token=cancellable)

    def wait(self, timeout: float = 60.0) -> None:
        """Wait, keeping the window responsive, until the background tasks have ended.

        Raises
        ------
        TimeoutError
            If they are still running after ``timeout`` seconds.
        """
        import time

        from PySide6.QtWidgets import QApplication

        end = time.monotonic() + timeout
        while self._window.stage.tasks.tasks:
            if time.monotonic() > end:
                raise TimeoutError("The background tasks did not end in time")
            QApplication.processEvents()
            time.sleep(0.005)

    def wait_until_closed(self) -> None:
        """Keep the window working until the user closes it, then return.

        Scripts started with plain ``python`` end with this call: without
        it the script would end, and the window close, at once. IPython and
        Jupyter do not need it, since they keep the window working between
        commands (see :func:`miles.launch`). :kbd:`Ctrl C` in the
        terminal stops waiting and raises :class:`KeyboardInterrupt`; the
        window stays open.
        """
        import signal
        import threading

        from PySide6.QtCore import QEventLoop, QTimer

        loop = QEventLoop()
        timer = QTimer()
        timer.timeout.connect(lambda: loop.quit() if not self._window.isVisible() else None)
        timer.start(50)
        interrupted: list[bool] = []

        def interrupt(*_: Any) -> None:
            interrupted.append(True)
            loop.quit()

        main = threading.current_thread() is threading.main_thread()
        previous = signal.signal(signal.SIGINT, interrupt) if main else None
        try:
            if self._window.isVisible():
                loop.exec()
        finally:
            timer.stop()
            if main:   # a handler installed outside Python reads back as None
                signal.signal(signal.SIGINT, previous if previous is not None else signal.default_int_handler)
        if interrupted:
            raise KeyboardInterrupt

    def look_up(self, text: str) -> None:
        """Look a text up in the fragment window, as if typed in its search box.

        SMILES is built at once; other text goes to a structure source
        (:meth:`~miles.plugins.PluginContext.add_structure_source`) in
        the background. A whole structure, such as a PDB entry, opens in a
        new tab; a molecule becomes the fragment to place.
        """
        self._window.look_up(text)

    def fetch(self, text: str) -> "DocumentView":
        """Build SMILES or fetch a code or name from a structure source, and open it in a new tab.

        The fetch runs on a worker thread while the window stays responsive;
        this returns when it is done.

        Parameters
        ----------
        text : str
            SMILES, or what a source takes, e.g. a PDB code or a compound name.

        Returns
        -------
        DocumentView
            The new tab.

        Raises
        ------
        FileFormatError
            If the text cannot be built or found.
        """
        import threading

        from PySide6.QtCore import QEventLoop, QTimer

        from miles.ui.lookup import look_up

        result: dict[str, Any] = {}

        def work() -> None:
            try:
                result["value"] = look_up(text, self._window.structure_sources)
            except Exception as exc:  # reported in the calling thread
                result["error"] = exc

        thread = threading.Thread(target=work, daemon=True)
        thread.start()
        loop = QEventLoop()
        timer = QTimer()
        timer.timeout.connect(lambda: loop.quit() if not thread.is_alive() else None)
        timer.start(20)
        loop.exec()
        timer.stop()
        if "error" in result:
            raise result["error"]
        mol, info, _ = result["value"]
        self._window.open_structure(mol, info if info is not None else JobInfo(title=text.strip()))
        return self._current_or_raise()

    def set_tool(self, tool: str) -> None:
        """Choose the mouse tool: ``"select"``, ``"build"`` or ``"delete"``."""
        from miles.ui.builder import Tool

        self._window.builder.set_tool(Tool(tool))

    def ask_save_path(self, title: str, name: str, filters: str = "All files (*)") -> str | None:
        """Ask where to save a file, starting in the folder used last.

        Parameters
        ----------
        title : str
            Dialog title.
        name : str
            Suggested file name.
        filters : str, optional
            Qt file filters, e.g. ``"SVG drawing (*.svg);;PNG image (*.png)"``.

        Returns
        -------
        str or None
            The chosen path, or None if the dialog was cancelled.
        """
        from PySide6.QtWidgets import QFileDialog

        start = os.path.join(self._window.last_dir(), name)
        path, _ = QFileDialog.getSaveFileName(self._window, title, start, filters)
        if not path:
            return None
        self._window.remember_dir(path)
        return path

    def ask_open_path(self, title: str, filters: str = "All files (*)") -> str | None:
        """Ask for a file to read, starting in the folder used last (None if cancelled)."""
        from PySide6.QtWidgets import QFileDialog

        path, _ = QFileDialog.getOpenFileName(self._window, title, self._window.last_dir(), filters)
        if not path:
            return None
        self._window.remember_dir(path)
        return path

    def ask_folder(self, title: str) -> str | None:
        """Ask for a folder to write files into, starting in the folder used last (None if cancelled)."""
        from PySide6.QtWidgets import QFileDialog

        path = QFileDialog.getExistingDirectory(self._window, title, self._window.last_dir())
        if not path:
            return None
        self._window.remember_dir(path, folder=True)
        return path

    def _current_or_raise(self) -> "DocumentView":
        view = self.current
        if view is None:
            raise RuntimeError("No document is open")
        return view

    def __repr__(self) -> str:
        return f"<MILES app, {len(self._window.canvases)} documents>"
