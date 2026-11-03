"""The public modules: ``miles``, ``miles.plugin`` with its submodules, and ``miles.testing``."""

import re
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

import miles
import miles.plugin
import miles.plugin.render
import miles.plugin.widgets
import miles.testing

SUBMODULES = (miles.plugin.widgets, miles.plugin.render)
"""The provisional submodules of :mod:`miles.plugin`."""

DOCS = Path(__file__).resolve().parent.parent / "docs" / "api"


def _documented(prefix: str) -> set[str]:
    """Names documented with a directive for ``prefix.<name>`` in the API reference."""
    if not DOCS.is_dir():
        pytest.skip("the documentation is not next to the tests (an installed package is tested)")
    text = "\n".join(path.read_text() for path in DOCS.glob("*.rst"))
    pattern = rf"^\.\. (?:auto)?(?:function|class|data|attribute|exception):: {re.escape(prefix)}\.(\w+)\s*$"
    return set(re.findall(pattern, text, re.MULTILINE))


def test_import_miles_loads_neither_qt_nor_open_babel() -> None:
    code = textwrap.dedent("""
        import sys
        import miles
        mol = miles.Molecule()
        heavy = [name for name in ("PySide6", "openbabel") if name in sys.modules]
        print(",".join(heavy))
    """)
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert result.stdout.strip() == ""


def test_every_public_name_resolves() -> None:
    for name in miles.__all__:
        assert getattr(miles, name) is not None, name
    for module in (miles.plugin, *SUBMODULES, miles.testing):
        for name in module.__all__:
            assert hasattr(module, name), f"{module.__name__}.{name}"
    assert set(dir(miles)) >= set(miles.__all__)


def test_the_library_comes_from_the_core_never_from_plugins() -> None:
    for name in miles.__all__:
        value = getattr(miles, name)
        module = getattr(value, "__module__", None) or miles._LAZY.get(name, "miles")
        assert not module.startswith("miles.builtin_plugins"), name


def test_each_name_has_one_home() -> None:
    modules = (miles, miles.plugin, *SUBMODULES)
    for k, first in enumerate(modules):
        for second in modules[k + 1:]:
            assert set(first.__all__) & set(second.__all__) == set(), (first.__name__, second.__name__)


def test_public_names_are_the_names_of_what_they_stand_for() -> None:
    """A public name is the object's own name, so the reference and tracebacks agree (no aliases)."""
    for module in (miles, miles.plugin, *SUBMODULES):
        for name in set(module.__all__) - {"__version__"}:
            value = getattr(module, name)
            own = getattr(value, "__name__", name)
            assert own == name or not callable(value), f"{module.__name__}.{name} is {own}"


def test_every_public_name_is_documented() -> None:
    assert set(miles.__all__) - {"__version__"} <= _documented("miles")
    assert set(miles.plugin.__all__) <= _documented("miles.plugin")
    for module in (*SUBMODULES, miles.testing):
        assert set(module.__all__) - {"miles_app", "qt_application"} <= _documented(module.__name__), module
