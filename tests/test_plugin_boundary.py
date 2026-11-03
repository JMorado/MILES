"""Where code lives, and what plugin code may import.

MILES keeps code in three places:

* :mod:`miles.core`, with the rest of the application (files, render, ui,
  api): the app and its foundations. Every core module is used by code that
  is not a plugin, and so is every name of :mod:`miles`, apart from a few
  listed in ``SCRIPT_NAMES`` that scripts need.
* ``miles.builtin_plugins.<name>``: one feature each, its interface and its
  algorithms.
* ``miles.builtin_plugins.shared``: helpers that several built-in plugins
  use and the rest of the application does not.

Plugin code (the built-ins, the shared package and the example plugins)
imports only the two public modules, :mod:`miles` and :mod:`miles.plugin`
(with its provisional submodules), its own package, the shared package (built-ins only), the standard library
and what MILES depends on (NumPy, PySide6 for widgets, Open Babel). Tests of
a plugin may also use :mod:`miles.testing`. Plugin code never imports another
plugin, and never reaches for the window, a canvas or underscore names of
objects that are not its own.
"""

import ast
import sys
from dataclasses import dataclass
from pathlib import Path

import pytest

import miles

PACKAGE = Path(miles.__file__).resolve().parent   # the source tree, or the installed package
ROOT = PACKAGE.parent
BUILTIN = PACKAGE / "builtin_plugins"
EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "plugins"
SHARED = "miles.builtin_plugins.shared"
OTHERS = {"numpy", "PySide6", "openbabel", "__future__"}   # what MILES itself depends on
INTERNALS = {"canvas", "window"}

UNUSED_CORE_MODULES: dict[str, str] = {}
"""Core modules allowed to have no user outside plugins, with the reason (none today)."""

SCRIPT_NAMES: dict[str, str] = {
    "launch": "starts the application from a script or the console",
    "write_cube": "the writer that matches the cube reader; one format's code stays together",
    "residue_ids": "the residue grouping of the selection language, offered as a function",
    "align_frames": "trajectory geometry for scripts; Tools > Align Frames calls it",
    "conformers": "a conformer search for scripts; Tools > Generate Conformers calls it",
    "amino_acid_codes": "the amino acids of the fragment library the builder uses",
    "amino_acid_template": "the amino acids of the fragment library the builder uses",
    "is_polymer": "the per-atom form of has_polymer, which the application uses",
    "is_water": "the per-atom form of water_mask, which the application uses",
    "dipole_moment": "a property of the charges the application computes; the Charges plugin draws it",
}
"""Names of :mod:`miles` that only scripts and plugins use, with the reason they still belong there."""


@dataclass(frozen=True)
class Plugin:
    """Plugin code: its files, the package its own modules live under, and whether it is built in.

    ``base`` is the folder that dotted module names start from, so that a
    relative import can be read as the absolute one it stands for.
    """

    name: str
    own: str
    files: tuple[Path, ...]
    builtin: bool
    base: Path = ROOT


def _plugins() -> list[Plugin]:
    found = []
    for path in sorted(BUILTIN.iterdir()):
        if path.name.startswith(("_", ".")):
            continue
        if path.is_dir():
            files = tuple(sorted(path.rglob("*.py")))
        elif path.suffix == ".py":
            files = (path,)
        else:
            continue
        found.append(Plugin(f"built-in {path.stem}", f"miles.builtin_plugins.{path.stem}", files, True))
    for path in sorted(EXAMPLES.glob("*.py")):
        found.append(Plugin(f"example {path.stem}", path.stem, (path,), False, EXAMPLES))
    for package in sorted(EXAMPLES.glob("*/*/__init__.py")):   # installable examples such as the template
        folder = package.parent
        files = tuple(sorted(folder.parent.rglob("*.py")))
        found.append(Plugin(f"example {folder.parent.name}", folder.name, files, False, folder.parent))
    return found


def _within(module: str, prefix: str) -> bool:
    return module == prefix or module.startswith(prefix + ".")


def _allowed(module: str, plugin: Plugin, path: Path) -> bool:
    top = module.split(".")[0]
    if top in sys.stdlib_module_names or top in OTHERS:
        return True
    if module == "miles" or _within(module, "miles.plugin") or _within(module, plugin.own):
        return True
    if plugin.builtin and _within(module, SHARED):
        return True
    return module == "miles.testing" and path.name.startswith(("test_", "conftest"))


def _private(name: str) -> bool:
    return name.startswith("_") and not (name.startswith("__") and name.endswith("__"))


def _reaches_past(node: ast.Attribute) -> bool:
    """Whether an attribute is the window or a canvas, or an underscore name of an object not its own."""
    own = isinstance(node.value, ast.Name) and node.value.id in ("self", "cls")
    return node.attr in INTERNALS or (_private(node.attr) and not own)


def _absolute(node: ast.ImportFrom, path: Path, base: Path) -> str:
    """The module a ``from ... import`` names, with a relative one read from the file's package."""
    if node.level == 0:
        return node.module or ""
    package = list(path.relative_to(base).with_suffix("").parts[:-1])
    if node.level > 1:
        package = package[:1 - node.level]
    return ".".join([*package, *([node.module] if node.module else [])])


def problems(plugin: Plugin) -> list[str]:
    """Imports and attribute uses of plugin code that go past what it may use."""
    found = []
    for path in plugin.files:
        where = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path.name
        for node in ast.walk(ast.parse(path.read_text(), str(path))):
            if isinstance(node, ast.Import):
                found += [f"{where}:{node.lineno} imports {alias.name}" for alias in node.names
                          if not _allowed(alias.name, plugin, path)]
            elif isinstance(node, ast.ImportFrom):
                module = _absolute(node, path, plugin.base)
                if node.module is None:   # from . import name: each name is a module
                    names = [f"{module}.{alias.name}" if module else alias.name for alias in node.names]
                    found += [f"{where}:{node.lineno} imports {name}" for name in names
                              if not _allowed(name, plugin, path)]
                    continue
                if not _allowed(module, plugin, path):
                    found.append(f"{where}:{node.lineno} imports from {module}")
                elif module == "miles":
                    found += [f"{where}:{node.lineno} imports {alias.name} from miles, which is not public"
                              for alias in node.names if alias.name not in miles.__all__]
                found += [f"{where}:{node.lineno} imports {alias.name} from {module}"
                          for alias in node.names if _private(alias.name) and module.startswith("miles")
                          and not _within(module, plugin.own)]
            elif isinstance(node, ast.Attribute) and _reaches_past(node):
                found.append(f"{where}:{node.lineno} uses .{node.attr}")
    return found


def _imported_modules(path: Path) -> set[str]:
    modules = set()
    for node in ast.walk(ast.parse(path.read_text(), str(path))):
        if isinstance(node, ast.Import):
            modules |= {alias.name for alias in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            modules.add(node.module)
            modules |= {f"{node.module}.{alias.name}" for alias in node.names}   # from miles.core import chem
    return modules


def _application_files() -> list[Path]:
    """Every module of MILES that is not plugin code."""
    return [path for path in sorted(PACKAGE.rglob("*.py")) if not path.is_relative_to(BUILTIN)]


FRONTS = {PACKAGE / "__init__.py", *(PACKAGE / "plugin").glob("*.py")}
"""The public modules: they only pass names on to scripts and plugins."""


def _names_used(path: Path) -> set[str]:
    """Names a module refers to: plain names, imported names and attributes of a module or object."""
    names = set()
    for node in ast.walk(ast.parse(path.read_text(), str(path))):
        if isinstance(node, ast.Name):
            names.add(node.id)
        elif isinstance(node, ast.ImportFrom):
            names |= {alias.name for alias in node.names}
        elif isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) \
                and node.value.id not in ("self", "cls"):
            names.add(node.attr)   # internal.set_angle, perception.rebond
    return names


PLUGINS = _plugins()


def test_every_plugin_is_checked() -> None:
    names = {plugin.name for plugin in PLUGINS}
    assert {"built-in vibrations", "built-in charges", "built-in shared"} <= names
    assert {"example selection_info", "example template"} <= names


@pytest.mark.parametrize("plugin", PLUGINS, ids=[p.name for p in PLUGINS])
def test_plugin_uses_only_the_public_surface(plugin: Plugin) -> None:
    assert problems(plugin) == []


def test_the_check_catches_reaches_into_internals(tmp_path: Path) -> None:
    source = tmp_path / "bad.py"
    source.write_text("from miles.ui.canvas import Canvas\nfrom miles.core.molecule import _Table\n"
                      "import miles.render.scene\nfrom miles import core\n"
                      "from miles.builtin_plugins.scan import relaxed_scan\n\n"
                      "def run(app, doc):\n    app.window.close()\n"
                      "    doc.canvas.redraw()\n    doc._canvas.redraw()\n")
    found = problems(Plugin("bad", "miles.builtin_plugins.bad", (source,), True, tmp_path))
    assert len(found) == 9


def test_the_check_reads_relative_imports(tmp_path: Path) -> None:
    package = tmp_path / "miles" / "builtin_plugins" / "bad"
    package.mkdir(parents=True)
    source = package / "__init__.py"
    source.write_text("from . import helpers\nfrom .helpers import tidy\nfrom ..scan import relaxed_scan\n"
                      "from .. import charges\nfrom ...core.molecule import Molecule\n"
                      "from ...plugin.widgets import Table\n")
    found = problems(Plugin("bad", "miles.builtin_plugins.bad", (source,), True, tmp_path))
    assert found == [f"{source.name}:3 imports from miles.builtin_plugins.scan",
                     f"{source.name}:4 imports miles.builtin_plugins.charges",
                     f"{source.name}:5 imports from miles.core.molecule"]


def test_the_application_never_imports_plugin_code() -> None:
    found = [f"{path.relative_to(ROOT)} imports {module}" for path in _application_files()
             for module in _imported_modules(path) if _within(module, "miles.builtin_plugins")]
    assert found == []


def test_every_core_module_is_used_by_the_application() -> None:
    """A core module that only plugins use belongs to its plugin, or to the shared package."""
    core = {path.stem for path in (PACKAGE / "core").glob("*.py") if path.stem != "__init__"}
    users: dict[str, set[str]] = {name: set() for name in core}
    for path in _application_files():
        if path in FRONTS:
            continue
        user = path.stem if path.parent == PACKAGE / "core" else "app"
        for module in _imported_modules(path):
            parts = module.split(".")
            if len(parts) >= 3 and parts[:2] == ["miles", "core"] and parts[2] in core and parts[2] != user:
                users[parts[2]].add(user)
    used = {name for name, found in users.items() if "app" in found}
    while more := {name for name in core - used if users[name] & used}:
        used |= more
    assert core - used == set(UNUSED_CORE_MODULES)


def test_every_name_of_miles_is_used_by_the_application() -> None:
    """A name of :mod:`miles` that only plugins use is plugin code, unless scripts need it (listed above)."""
    used = set().union(*(_names_used(path) for path in _application_files() if path not in FRONTS))
    unused = {name for name in miles.__all__ if name not in used and not name.startswith("__")}
    assert unused == set(SCRIPT_NAMES)
