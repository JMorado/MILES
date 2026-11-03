"""Opening the MILES window from IPython, Jupyter or a Python script.

:func:`launch`, also reachable as :func:`miles.launch`, opens the
window in the running Python process and returns its
:class:`~miles.api.App`, the object the Python console calls ``app``.
A window only responds while Qt's event loop runs, and who runs it depends
on the session:

IPython and Jupyter
    :func:`launch` turns on IPython's Qt integration, as ``%gui qt`` does.
    The kernel (Jupyter, JupyterLab, qtconsole) or the terminal prompt
    then runs the loop whenever it waits for the next cell or line, so the
    window can be used between cells and every call on ``app`` shows at
    once.
Plain Python
    The plain ``python`` prompt cannot run the loop between lines, so the
    window does not respond there; :func:`launch` says so when called at
    that prompt. A script calls
    :meth:`~miles.api.App.wait_until_closed` at the end, which keeps
    the window working until the user closes it.

Errors raised inside the window, in a menu command or a panel, show in
the error dialog of MILES and their traceback is printed, as in the
``miles`` command. In IPython the dialog becomes ``sys.excepthook``
when the cell that called :func:`launch` is over, and it handles the
errors raised while IPython waits for input; the traceback goes to the
output of the cell that ran last. While a cell runs, IPython puts its own
hook in place, so an error in a command the cell starts is shown as a
normal error in that cell's output.
"""

import os
import sys
from typing import TYPE_CHECKING, Any

from miles.render.theme import Appearance

if TYPE_CHECKING:
    from PySide6.QtCore import QSettings

    from miles.api import App
    from miles.ui.mainwindow import MainWindow

PLAIN_PROMPT_NOTE = (
    "The MILES window only responds while Qt's event loop runs, which the plain Python prompt "
    "does not do between lines. Use IPython or Jupyter, where launch() keeps the window working, or "
    "call app.wait_until_closed() to use the window until you close it.\n"
)

_launched: "MainWindow | None" = None   # the window made by launch(), deleted once closed


def launch(*paths: str, plugins: bool = True, appearance: str | None = None,
           settings: "QSettings | None" = None) -> "App":
    """Open the MILES window (or bring back the open one) and return its ``app``.

    Parameters
    ----------
    *paths : str
        Files to open, each in a tab of its own.
    plugins : bool, optional
        Load third-party plugins as well as the built-in ones, as the
        ``miles`` command does without ``--no-plugins``.
    appearance : {"system", "light", "dark"}, optional
        Theme of the window, for this window only (View > Appearance
        remembers its choice; this does not). By default the saved choice
        is used.
    settings : QSettings, optional
        Where the window keeps its preferences, such as recent files and
        panel places (default: the user's MILES settings). Mostly
        useful for tests.

    Returns
    -------
    App
        The window, as scripts and the console see it.

    Raises
    ------
    FileNotFoundError
        If one of ``paths`` does not exist (checked before anything opens).
    ValueError
        If ``appearance`` is not one of the names above.
    RuntimeError
        If IPython already runs the event loop of another toolkit (tk, wx,
        the macOS backend of matplotlib...), since only one can run.

    Notes
    -----
    While a window is open, calling :func:`launch` again returns the same
    ``app``, opens the new paths in it and brings it to the front;
    ``plugins`` and ``settings`` only apply when a window is made. Once the
    window is closed, the next call makes a new one.
    """
    from PySide6.QtWidgets import QApplication

    missing = [path for path in paths if not os.path.exists(path)]
    if missing:
        raise FileNotFoundError(missing[0])
    theme = Appearance(appearance) if appearance is not None else None
    shell = _ipython()
    if shell is not None:
        _check_event_loop(shell)

    app = running_app()
    if app is None:
        from miles.ui.app import application, build_window

        qt = application()
        qt.setQuitOnLastWindowClosed(False)   # closing the window must not end the kernel's Qt loop
        _forget_closed_window()
        window = build_window(plugins=plugins, settings=settings)
        _remember(window)
        _report_errors(shell)
        app = window.api
    if theme is not None:
        _set_appearance(app.window, theme)
    for path in paths:
        app.open(path)
    app.window.raise_()
    app.window.activateWindow()
    if shell is not None:
        _start_event_loop(shell)
    elif _at_plain_prompt():
        sys.stderr.write(PLAIN_PROMPT_NOTE)
    QApplication.processEvents()   # draw the window now, whoever runs the loop later
    return app


def running_app() -> "App | None":
    """Return the ``app`` of the open MILES window, or None if no window is open.

    This is the window of the ``miles`` command as well as one made by
    :func:`launch`; a closed window does not count.
    """
    import shiboken6

    from miles.api import current_app

    try:
        app = current_app()
    except RuntimeError:
        return None
    window = app.window
    if not shiboken6.isValid(window) or not window.isVisible():
        return None
    return app


def _remember(window: "MainWindow") -> None:
    global _launched
    _launched = window


def _forget_closed_window() -> None:
    """Delete the window made by an earlier :func:`launch` once it has been closed."""
    global _launched
    import shiboken6

    if _launched is not None and shiboken6.isValid(_launched) and not _launched.isVisible():
        _launched.deleteLater()
    _launched = None


def _ipython() -> Any:
    """Return the running IPython shell (terminal or kernel), or None in plain Python."""
    module = sys.modules.get("IPython")
    if module is None:
        return None
    try:
        return module.get_ipython()
    except Exception:  # an IPython too old or half imported counts as none
        return None


def _active_loop(shell: Any) -> str:
    return str(getattr(shell, "active_eventloop", None) or "")


def _check_event_loop(shell: Any) -> None:
    """Refuse to start when IPython runs another toolkit's event loop."""
    loop = _active_loop(shell)
    if loop and not loop.startswith("qt"):
        raise RuntimeError(f"IPython is running the {loop!r} event loop and only one can run, but the "
                           "MILES window needs Qt's. Start a new session and call launch() before "
                           "using other windows, or use Qt for those too (for matplotlib: %matplotlib qt).")


def _start_event_loop(shell: Any) -> None:
    """Have IPython run Qt's event loop while it waits for input, as ``%gui qt`` does."""
    if not _active_loop(shell).startswith("qt"):   # "qt", or "qt6" in the terminal
        shell.enable_gui("qt")
    _keep_kernel_awake(shell)


_kernel_poll: Any = None   # the timer of _keep_kernel_awake, made once per process

KERNEL_POLL_MS = 50
"""How often a Jupyter kernel's Qt loop looks for requests it was not woken for."""


def _keep_kernel_awake(shell: Any) -> None:
    """Make sure a Jupyter kernel running Qt's loop gets to every request.

    ipykernel runs cells on an asyncio loop and Qt's loop inside it,
    leaving Qt's loop only when a socket notifier reports a new message.
    A request read by the asyncio side just before the kernel enters Qt's
    loop is then stuck: it waits for the previous cell's handler to finish,
    which cannot happen while Qt's loop blocks asyncio, and the cell waits
    until some other message arrives (this happens with a bare ``%gui qt``
    too). This timer gives the asyncio side a turn at a fixed interval, as
    ipykernel itself does on Windows: while the kernel waits in Qt's loop,
    it takes a request the notifier missed and leaves the loop, which the
    kernel enters again a millisecond later. It uses a few ipykernel
    internals, so anything missing simply turns it off.
    """
    global _kernel_poll
    kernel = getattr(shell, "kernel", None)   # only kernels have one; the terminal reads stdin directly
    if kernel is None or _kernel_poll is not None:
        return
    try:
        from ipykernel.eventloops import get_shell_stream
        from PySide6.QtCore import QTimer
        from PySide6.QtWidgets import QApplication

        stream = get_shell_stream(kernel)
    except Exception:  # another kernel, or an ipykernel without these parts
        return

    def poll() -> None:
        try:
            notifier = getattr(kernel, "_qt_notifier", None)
            if notifier is None or not notifier.isEnabled():
                return   # the kernel is not waiting in the Qt loop
            stream.flush(limit=1)              # a request the notifier missed
            notifier.setEnabled(False)         # leave the loop as ipykernel does
            kernel.app.qt_event_loop.quit()
        except Exception:  # never let the workaround itself break the window
            _kernel_poll.stop()

    _kernel_poll = QTimer(QApplication.instance())
    _kernel_poll.timeout.connect(poll)
    _kernel_poll.start(KERNEL_POLL_MS)


def _at_plain_prompt() -> bool:
    """Whether this is the interactive prompt of plain Python (``python``, ``python -i``)."""
    return hasattr(sys, "ps1") or bool(sys.flags.interactive)


def _report_errors(shell: Any) -> None:
    """Show errors raised inside the window in the error dialog.

    IPython swaps its own hook in while a cell runs and puts the previous
    one back when the code ends, so in IPython the dialog is installed
    once the cell that called :func:`launch` is over (its
    ``post_run_cell`` event).
    """
    from miles.ui import errors

    def window() -> "MainWindow | None":
        app = running_app()
        return app.window if app is not None else None

    if shell is None:
        errors.install(window)
        return

    def after_cell(*_: Any) -> None:
        shell.events.unregister("post_run_cell", after_cell)
        errors.install(window)

    shell.events.register("post_run_cell", after_cell)


def _set_appearance(window: "MainWindow", appearance: Appearance) -> None:
    """Use a theme without storing it as the user's choice."""
    window.appearance = str(appearance)
    window.apply_theme()
    window.menu_actions.check("appearance", str(appearance))
