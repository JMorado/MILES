"""Testing plugins: pytest fixtures and helpers for a MILES window that runs without a screen.

A plugin's tests load the fixtures with one line in the ``conftest.py`` at
the top of the plugin's folder::

    pytest_plugins = ["miles.testing"]

and then drive a real window through :class:`~miles.plugin.App`::

    import my_plugin
    from miles.testing import document_from_smiles, load_plugin

    def test_command(miles_app):
        info = load_plugin(miles_app, my_plugin)
        doc = document_from_smiles(miles_app, "CCO", "ethanol")
        miles_app.run("Tools > My Command")
        assert not info.error

Importing this module prepares the environment before Qt starts: no
screen is needed (``QT_QPA_PLATFORM=offscreen``), the canvas draws with
QPainter, and the plugin and fragment folders and the autosave folder of
daily use are replaced by empty temporary ones. Each window has settings of
its own, so tests neither read nor change the user's. Only the built-in
plugins load, besides the ones a test adds with :func:`load_plugin`.

What it provides:

The window
    :func:`miles_app` (a fixture), :func:`load_plugin`,
    :func:`document_from_smiles`, :func:`document_from_text`,
    :func:`wait_until`.
The plugin's life
    :func:`switch_off`, :func:`switch_on` and :func:`reload_plugin`, as
    Help > Plugins does them, and :func:`plugin_setting`, a setting of the
    plugin by name.
What the user sees
    :func:`captured_messages` (notes and errors, without dialogs) and
    :func:`drawings` (what the plugin's layers draw for a tab).

Like :mod:`miles.plugin.widgets`, this module is provisional: it follows
the window and may change between releases, which may need small changes
to a plugin's tests. It needs pytest, which MILES itself does not.
"""

import atexit
import os
import shutil
import tempfile
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

_TEMPORARY = tempfile.mkdtemp(prefix="miles-testing-")
atexit.register(shutil.rmtree, _TEMPORARY, True)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MILES_SOFTWARE", "1")
for _variable, _folder in (("MILES_PLUGIN_DIR", "plugins"), ("MILES_FRAGMENT_DIR", "fragments"),
                           ("MILES_AUTOSAVE_DIR", "autosave")):
    os.environ.setdefault(_variable, os.path.join(_TEMPORARY, _folder))

import pytest  # noqa: E402
from PySide6.QtCore import QSettings  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from miles.api import App, DocumentView, Severity  # noqa: E402
from miles.core.chem import from_smiles  # noqa: E402
from miles.files.formats import parse_text  # noqa: E402
from miles.plugins import DrawContext, PluginInfo, PluginState, Setting  # noqa: E402

__all__ = ["Message", "captured_messages", "document_from_smiles", "document_from_text", "drawings",
           "load_plugin", "miles_app", "plugin_setting", "qt_application", "reload_plugin", "switch_off",
           "switch_on", "wait_until"]


@pytest.fixture(scope="session")
def qt_application() -> QApplication:
    """Return the Qt application, made once for all the tests."""
    instance = QApplication.instance()
    return instance if isinstance(instance, QApplication) else QApplication([])


@pytest.fixture
def miles_app(qt_application: QApplication, tmp_path: Path) -> Iterator[App]:
    """Yield a new MILES window with an empty tab, closed after the test without asking to save."""
    from miles.ui.mainwindow import MainWindow

    window = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat), plugins=False)
    window.resize(1000, 700)
    window.show()
    qt_application.processEvents()
    yield window.api
    for canvas in window.canvases:
        canvas.doc.modified = False
    window.close()
    window.deleteLater()
    qt_application.processEvents()


def load_plugin(app: App, module: ModuleType) -> PluginInfo:
    """Register a plugin's module in the window, as MILES does at start-up.

    Parameters
    ----------
    app : App
        The window, e.g. the :func:`miles_app` fixture.
    module : module
        The plugin, with its ``register(context)`` function.

    Returns
    -------
    PluginInfo
        The plugin; its ``error`` holds the traceback of any later failure
        of its callbacks, so a test can end with ``assert not info.error``.

    Raises
    ------
    AssertionError
        If the plugin failed to register, with its traceback.
    """
    info = app.window.plugins.add_module(module)
    if info.state != PluginState.LOADED:
        raise AssertionError(f"The plugin {info.name} failed to register:\n{info.error}")
    return info


def document_from_smiles(app: App, smiles: str, title: str = "") -> DocumentView:
    """Open a new tab with a structure built from SMILES, with 3D coordinates and hydrogens."""
    molecule, _ = from_smiles(smiles)
    return app.new(molecule, title or smiles)


def document_from_text(app: App, text: str, title: str = "structure") -> DocumentView:
    """Open a new tab with a structure from text, in any format MILES reads when pasted (XYZ, PDB ...)."""
    return app.new(parse_text(text), title)


def wait_until(condition: Callable[[], object], timeout: float = 5.0) -> None:
    """Keep the window working until ``condition()`` is true, e.g. for work in the background.

    Raises
    ------
    TimeoutError
        If it is still false after ``timeout`` seconds.
    """
    app = QApplication.instance()
    end = time.monotonic() + timeout
    while not condition():
        if time.monotonic() > end:
            raise TimeoutError("The condition did not become true in time")
        if app is not None:
            app.processEvents()
        time.sleep(0.005)


def _manager(app: App) -> Any:
    return app.window.plugins


def _loaded(info: PluginInfo) -> PluginInfo:
    if info.state != PluginState.LOADED:
        raise AssertionError(f"The plugin {info.name} failed to register:\n{info.error}")
    return info


def switch_off(app: App, info: PluginInfo) -> None:
    """Switch a plugin off, as its check box in Help > Plugins does: everything it added goes."""
    _manager(app).set_enabled(info.key, False)


def switch_on(app: App, info: PluginInfo) -> PluginInfo:
    """Switch a plugin on again; returns it, loaded.

    Raises
    ------
    AssertionError
        If it failed to register, with its traceback.
    """
    _manager(app).set_enabled(info.key, True)
    return _loaded(info)


def reload_plugin(app: App, info: PluginInfo) -> PluginInfo:
    """Unload a plugin and load it again, as Reload in Help > Plugins does; returns it, loaded.

    A module that was imported from a file is imported again with its own
    submodules, so the plugin's module objects a test held before are
    replaced.

    Raises
    ------
    AssertionError
        If it failed to register, with its traceback.
    """
    _manager(app).reload(info)
    return _loaded(info)


def plugin_setting(app: App, info: PluginInfo, name: str) -> Setting:
    """Return a setting the plugin added with ``context.add_setting``, to read or ``set`` its value.

    Raises
    ------
    KeyError
        If the plugin has no setting of that name.
    """
    for setting in _manager(app).settings.get(info.key, []):
        if setting.name == name:
            return setting
    raise KeyError(f"{info.name} has no setting {name!r}")


@dataclass(frozen=True)
class Message:
    """Something a plugin told the user with ``app.message`` or ``app.warn``."""

    text: str
    severity: Severity
    title: str = ""


@contextmanager
def captured_messages(app: App) -> Iterator[list[Message]]:
    """Record what the window tells the user while the block runs.

    Notes still appear on the tab; errors are recorded instead of shown in a
    dialog, which would wait for a click::

        with captured_messages(miles_app) as messages:
            miles_app.run("Tools > My Command")
        assert messages[-1].severity == Severity.WARNING
    """
    found: list[Message] = []
    shown = app.message

    def record(text: str, seconds: float | None = None, severity: str = Severity.INFO,
               title: str = "") -> None:
        level = Severity(severity)
        found.append(Message(text, level, title))
        if level != Severity.ERROR:
            shown(text, seconds, level, title)

    app.message = record   # type: ignore[method-assign]
    try:
        yield found
    finally:
        del app.message


def drawings(app: App, info: PluginInfo, document: DocumentView | None = None,
             export: bool = False) -> list[Any]:
    """Return what a plugin's layers draw for a tab now: its annotations, scene items and status text.

    Parameters
    ----------
    app : App
        The window.
    info : PluginInfo
        The plugin, from :func:`load_plugin`.
    document : DocumentView, optional
        The tab (default: the current one).
    export : bool, optional
        Draw as for an exported picture, where slow drawings are worked out
        at once.

    Returns
    -------
    list
        The results that are not None, in the order the layers were added.
    """
    document = document or app.current
    if document is None:
        return []
    canvas = document.canvas
    found = []
    for handle in _manager(app)._layers:
        if handle.info is not info:
            continue
        context = None if handle.kind == "status" else DrawContext(app.theme, export=export)
        result = handle.layer.function(canvas, context)
        if result is not None:
            found.append(result)
    return found
