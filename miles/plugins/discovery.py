"""Finding plugins: the built-in list, entry points of installed packages and the plugin folder.

Also what is recorded about each plugin found (:class:`PluginInfo`) and the
constants of the plugin interface.
"""

import functools
import importlib
import importlib.metadata
import importlib.util
import os
import sys
import traceback
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from types import ModuleType
from typing import Any

from PySide6.QtCore import QStandardPaths

PACKAGE_DIR = Path(__file__).resolve().parent
"""Folder of the plugin machinery, whose frames are left out of plugin tracebacks."""

PLUGIN_API = 1
"""Version of the plugin interface. A plugin that sets ``REQUIRES_API`` to a
higher number fails to load with a message saying so, instead of failing
later in a way that is hard to understand."""

ENTRY_POINT_GROUP = "miles.plugins"
"""Entry point group that installed packages use to declare plugins."""

BUILTIN = (
    "miles.builtin_plugins.drawing",
    "miles.builtin_plugins.movie",
    "miles.builtin_plugins.structure_lookup",
    "miles.builtin_plugins.select_atoms",
    "miles.builtin_plugins.vibrations",
    "miles.builtin_plugins.sequence",
    "miles.builtin_plugins.conformers",
    "miles.builtin_plugins.scan",
    "miles.builtin_plugins.align_frames",
    "miles.builtin_plugins.superpose",
    "miles.builtin_plugins.prepare_protein",
    "miles.builtin_plugins.binding_pocket",
    "miles.builtin_plugins.ramachandran",
    "miles.builtin_plugins.properties",
    "miles.builtin_plugins.molecular_surface",
    "miles.builtin_plugins.charges",
    "miles.builtin_plugins.hydrogen_bonds",
    "miles.builtin_plugins.mutate_residue",
    "miles.builtin_plugins.grids",
)
"""Modules of the built-in plugins, loaded even when entry points are not installed.

The order matters where commands are placed before the same menu item
(the first one loaded ends up highest) and for the panels, which stack on
the right of the window in loading order.
"""

FOLDER_ENV = "MILES_PLUGIN_DIR"
"""Environment variable naming another plugin folder (used by the tests)."""

NO_PLUGINS_ENV = "MILES_NO_PLUGINS"
"""Environment variable that, when set, skips third-party plugins like ``--no-plugins``."""


class PluginState(StrEnum):
    """What happened to a plugin at start-up."""

    LOADED = "loaded"
    FAILED = "failed"
    OFF = "switched off"
    SKIPPED = "not loaded (--no-plugins)"


@dataclass
class PluginInfo:
    """A plugin found at start-up, and what became of it.

    Attributes
    ----------
    key : str
        Identifies the plugin in the settings: its module name for
        the built-in plugins, ``entry:<name>`` for entry points and
        ``file:<file name>`` for the plugin folder.
    name : str
        ``NAME`` of the module, or a name derived from where it was found.
    source : str
        Where it came from: "built in", the package that declares it, or
        "plugin folder".
    builtin : bool
        Whether it is one of the built-in plugins.
    version, description : str
        ``VERSION`` and ``DESCRIPTION`` of the module, if any.
    state : PluginState
        Loaded, failed, switched off or skipped.
    error : str
        Traceback of the last failure (import, registration or callback).
    warnings : str
        Notes that do not stop the plugin, one per line, such as a shortcut
        another command already has.
    commands : list of str
        Menu paths of the commands it added.
    extras : list of str
        What else it added, in a few words each: panels, drawing layers,
        status text, right-click items, structure sources, fragments.
    path : str
        File or package folder of a plugin from the plugin folder.
    """

    key: str
    name: str
    source: str
    builtin: bool = False
    version: str = ""
    description: str = ""
    state: PluginState = PluginState.SKIPPED
    error: str = ""
    warnings: str = ""
    commands: list[str] = field(default_factory=list)
    extras: list[str] = field(default_factory=list)
    path: str = ""


@dataclass
class Candidate:
    """A plugin found, with the function that imports it.

    Attributes
    ----------
    info : PluginInfo
        What is known about it so far.
    load : callable
        Returns the plugin's module, or its ``register`` function for an
        entry point that names the function itself.
    """

    info: PluginInfo
    load: Callable[[], Any]


def plugin_folder() -> Path:
    """Folder searched for plugin files (it may not exist yet).

    Returns
    -------
    pathlib.Path
        ``$MILES_PLUGIN_DIR`` if set, otherwise ``plugins`` in the
        application data folder: ``~/Library/Application Support/miles/
        MILES/plugins`` on macOS, ``%APPDATA%/miles/MILES/plugins``
        on Windows and ``~/.local/share/miles/MILES/plugins`` on Linux.
    """
    override = os.environ.get(FOLDER_ENV)
    if override:
        return Path(override)
    base = QStandardPaths.writableLocation(QStandardPaths.StandardLocation.AppDataLocation)
    return Path(base or Path.home() / ".miles") / "plugins"


def _file_loader(path: Path) -> Callable[[], ModuleType]:
    def load() -> ModuleType:
        name = f"miles_user_plugins.{path.stem}"
        forget_modules(name)   # so a package's own modules are read again, not taken from the last load
        target = path / "__init__.py" if path.is_dir() else path
        locations = [str(path)] if path.is_dir() else None
        spec = importlib.util.spec_from_file_location(name, target, submodule_search_locations=locations)
        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        try:
            spec.loader.exec_module(module)
        except BaseException:
            sys.modules.pop(name, None)
            raise
        return module
    return load


def forget_modules(name: str) -> None:
    """Drop a module and its submodules from :data:`sys.modules`, so the next import reads them again."""
    for loaded in [m for m in sys.modules if m == name or m.startswith(name + ".")]:
        del sys.modules[loaded]


def discover(third_party: bool = True) -> list[Candidate]:
    """Find the plugins to load, in loading order.

    Parameters
    ----------
    third_party : bool, optional
        Also look at entry points of other packages and at the plugin folder.

    Returns
    -------
    list
        The built-in plugins, then entry points sorted by name, then the
        plugin folder's files and packages sorted by name.
    """
    found = [Candidate(PluginInfo(name, name.rsplit(".", 1)[-1].replace("_", " ").title(), "built in",
                                   builtin=True), functools.partial(importlib.import_module, name))
             for name in BUILTIN]
    if not third_party:
        return found
    try:
        points = sorted(importlib.metadata.entry_points(group=ENTRY_POINT_GROUP), key=lambda ep: ep.name)
    except Exception:   # a broken installation must not stop the application
        points = []
    for point in points:
        if point.value.split(":")[0] in BUILTIN:
            continue   # the built-in entry points: loaded above, from the list
        dist = getattr(point, "dist", None)
        source = f"package {dist.name} {dist.version}" if dist is not None else f"entry point {point.value}"
        found.append(Candidate(PluginInfo(f"entry:{point.name}", point.name, source), point.load))
    folder = plugin_folder()
    if folder.is_dir():
        for path in sorted(folder.iterdir()):
            usable = (path.suffix == ".py" or (path.is_dir() and (path / "__init__.py").exists()))
            if usable and not path.name.startswith(("_", ".")):
                info = PluginInfo(f"file:{path.name}", path.stem, "plugin folder", path=str(path))
                found.append(Candidate(info, _file_loader(path)))
    return found


def plugin_traceback() -> str:
    """Traceback of the exception being handled, without the frames of the plugin loader."""
    kind, value, tb = sys.exc_info()
    if kind is None or value is None:
        return ""
    frames = [frame for frame in traceback.extract_tb(tb)
              if Path(frame.filename).resolve().parent != PACKAGE_DIR
              and not frame.filename.startswith("<frozen ")]
    head = ["Traceback (most recent call last):\n", *traceback.format_list(frames)] if frames else []
    return "".join(head + traceback.format_exception_only(kind, value))
