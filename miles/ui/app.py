"""Application entry point: the ``miles`` command.

:func:`application` and :func:`build_window` are shared with
:func:`miles.launch`, which opens the same window from IPython,
Jupyter or a script.
"""

import argparse
import os
import signal
import sys
from typing import TYPE_CHECKING

from miles import __version__

if TYPE_CHECKING:
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QApplication

    from miles.ui.mainwindow import MainWindow


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Command-line options."""
    parser = argparse.ArgumentParser(
        prog="miles",
        description="Minimal, fast molecule builder and viewer.",
        epilog='examples:  miles   |   miles ligand.xyz   |   miles 4hhb   |   '
               'miles caffeine "CC(=O)O"')
    parser.add_argument("items", nargs="*", metavar="ITEM",
                        help="files to open, PDB codes (e.g. 1crn), SMILES strings or molecule names")
    parser.add_argument("--software", action="store_true", help="draw with QPainter instead of OpenGL")
    parser.add_argument("--script", metavar="FILE",
                        help="run a Python script in the console once the window is open")
    parser.add_argument("-c", dest="code", metavar="CODE",
                        help="run Python code in the console once the window is open (after --script)")
    parser.add_argument("--no-plugins", action="store_true",
                        help="start without third-party plugins (from packages or the plugin folder)")
    parser.add_argument("--install-app", nargs="?", const="", metavar="FOLDER",
                        help="add MILES to the Dock, the Start menu or the application menu, under its own "
                             "name and icon, then exit (FOLDER replaces the usual place)")
    parser.add_argument("--version", action="version", version=f"miles {__version__}")
    return parser.parse_args(argv)


def application() -> "QApplication":
    """Return the Qt application, creating it if there is none yet.

    Asks for the OpenGL format the canvas draws with and sets the names
    under which Qt finds the folders of MILES (the plugin folder, for
    example). An application made by someone else, such as IPython's Qt
    integration, is used as it is, apart from the names.

    Raises
    ------
    RuntimeError
        If the running application cannot show widgets (a
        ``QCoreApplication`` or ``QGuiApplication``).
    """
    from PySide6.QtWidgets import QApplication

    from miles.ui import desktop, icons
    from miles.ui.canvas import configure_gl
    from miles.ui.mainwindow import APP_NAME

    configure_gl()   # before the application, or at least before the first canvas
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv[:1])
    elif not isinstance(app, QApplication):
        raise RuntimeError(f"A {type(app).__name__} is running, which cannot show windows; "
                           "MILES needs a QApplication")
    app.setApplicationName(APP_NAME)
    app.setOrganizationName("miles")
    # the application's icon, not the window's, is what the macOS Dock and the taskbars show
    app.setWindowIcon(icons.app_icon())
    desktop.identify_process(app)
    return app


def build_window(plugins: bool = True, settings: "QSettings | None" = None) -> "MainWindow":
    """Make the main window with its plugins and show it.

    Parameters
    ----------
    plugins : bool, optional
        Load third-party plugins too (built-in ones always load).
    settings : QSettings, optional
        Where preferences are kept (default: the user's settings).
    """
    from miles.ui.mainwindow import MainWindow

    window = MainWindow(settings, plugins=plugins)
    window.show()
    return window


def main(argv: list[str] | None = None) -> int:
    """Start the application and open the items given on the command line."""
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.install_app is not None:
        from pathlib import Path

        from miles.ui import desktop

        for made in desktop.install(Path(args.install_app).expanduser() if args.install_app else None):
            print(f"Made {made}")
        return 0

    from PySide6.QtCore import QTimer

    from miles.ui import errors
    from miles.ui.canvas import SOFTWARE_ENV
    from miles.ui.mainwindow import MainWindow

    if args.software:
        os.environ[SOFTWARE_ENV] = "1"
    app = application()
    windows: list[MainWindow] = []
    errors.install(lambda: windows[0] if windows else None)   # also covers errors while the window is built
    window = build_window(plugins=not args.no_plugins)
    windows.append(window)
    # Ctrl+C in the terminal quits. Python only sees the signal while it runs code,
    # so a timer wakes it regularly; the handler stops the event loop instead of
    # raising KeyboardInterrupt inside whatever Qt callback happens to be running.
    signal.signal(signal.SIGINT, lambda *_: app.exit(130))
    wake = QTimer()
    wake.timeout.connect(lambda: None)
    wake.start(200)
    for item in args.items:
        window.open_item(item)
    if args.script:   # once the event loop runs, so the window and the files are ready
        QTimer.singleShot(0, lambda: window.run_script(args.script))
    if args.code:
        QTimer.singleShot(0, lambda: window.run_code(args.code))
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
