"""The Python namespace of the console and of scripts run inside MILES.

Scripts and console commands see the same names:

``app``
    The running application (:class:`miles.api.App`).
``doc``, ``mol``
    The current tab (:class:`miles.api.DocumentView`) and its
    molecule. They follow the current tab before every command, unless you
    have bound the names to something else yourself.
``np``, ``chem``, ``building``, ``geometry``, ``Molecule``
    NumPy and the most used parts of :mod:`miles.core`.
Helpers
    ``load``, ``select``, ``count``, ``define``, ``show``, ``fit``, ``turn``,
    ``png``, ``run``, ``dist``, ``angle`` and ``dihedral``: short commands
    for the current tab (see :class:`Helpers`). ``select``, ``count`` and
    ``define`` take expressions of the selection language
    (:mod:`miles.core.selection`), such as ``"byres within 5 of ligand"``.

The helpers only call :mod:`miles.api`; anything they do can be
written out with ``app`` and ``doc`` as well. Code runs on the interface
thread, so a long computation keeps the window waiting until it ends; the
modules of :mod:`miles.core` also work in a plain Python session,
without the window, for heavy batch work.
"""

import builtins
from collections.abc import Iterable
from typing import Any

import numpy as np

import miles
from miles.api import App, DocumentView
from miles.core import geometry
from miles.core.molecule import Molecule, Representation

BANNER = (
    "Python in MILES. Ready: app (the window), doc (the current tab), mol (its molecule),\n"
    "np, miles (the chemistry library) and Molecule. Helpers: load, select, count, define, show, fit,\n"
    "turn, png, run, dist, angle, dihedral. Try select(\"protein and name CA\") or help(select).\n"
)

_BOUND = "__miles_bound__"   # doc and mol as last set, to tell whether the user rebound them


class Helpers:
    """Short commands for the current tab, in the spirit of PyMOL's command line.

    Parameters
    ----------
    app : App
        The application the helpers act on.
    """

    def __init__(self, app: App) -> None:
        self.app = app

    def _doc(self) -> DocumentView:
        doc = self.app.current
        if doc is None:
            raise RuntimeError("No document is open")
        return doc

    def load(self, path: str) -> DocumentView:
        """Open a file in a new tab and return it."""
        return self.app.open(path)

    def select(self, what: str | Iterable[int] | None = None, add: bool = False) -> list[int]:
        """Select atoms with an expression (a string) or by 0-based indices; return the selection.

        ``select()`` clears the selection; ``select("byres within 5 of ligand")``
        selects the residues around ligands; ``select('smarts "[OX2H]"')``
        selects hydroxyl oxygens; ``select([0, 1, 2], add=True)`` adds three atoms.
        """
        doc = self._doc()
        if what is None:
            return doc.select([], add=add)
        return doc.select(what, add=add)

    def count(self, expression: str) -> int:
        """Return how many atoms an expression picks, e.g. ``count("water within 3.5 of ligand")``."""
        return self._doc().count(expression)

    def define(self, name: str, expression: str | None) -> None:
        """Name an expression for later ones: ``define("pocket", "byres within 5 of ligand")``."""
        self._doc().define(name, expression)

    def show(self, mode: str) -> None:
        """Draw the current tab as ``"ball"``, ``"tube"``, ``"wire"``, ``"spacefill"`` or ``"cartoon"``."""
        names = [r.value for r in Representation if r not in (Representation.DEFAULT, Representation.HIDDEN)]
        if mode not in names:
            raise ValueError(f"Unknown style {mode!r}; use one of: {', '.join(names)}")
        self._doc().style(mode=mode)

    def fit(self) -> None:
        """Zoom and centre so the whole structure is in view."""
        self._doc().fit()

    def turn(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> None:
        """Turn the view about the screen axes, in degrees."""
        self._doc().turn(x, y, z)

    def png(self, path: str, scale: float = 2.0, background: str = "white") -> None:
        """Save a picture of the view (``background``: white, transparent or theme)."""
        self._doc().export_image(path, scale, background)

    def run(self, command: str) -> None:
        """Run a menu command by name, e.g. ``run("Clean Geometry")``."""
        self.app.run(command)

    def _coords(self, *atoms: int) -> list[np.ndarray]:
        mol = self._doc().molecule
        n = len(mol)
        if any(not 0 <= int(i) < n for i in atoms):
            raise ValueError(f"Atom indices must be between 0 and {n - 1}")
        return [mol.coords[int(i)] for i in atoms]

    def dist(self, i: int, j: int) -> float:
        """Distance in angstrom between atoms ``i`` and ``j`` (0-based) of the current tab."""
        return geometry.distance(*self._coords(i, j))

    def angle(self, i: int, j: int, k: int) -> float:
        """Angle i-j-k in degrees."""
        return geometry.angle(*self._coords(i, j, k))

    def dihedral(self, i: int, j: int, k: int, m: int) -> float:
        """Dihedral angle i-j-k-m in degrees."""
        return geometry.dihedral(*self._coords(i, j, k, m))


def _help(*args: Any, **kwargs: Any) -> None:
    """Show help on an object, as :func:`help` does (without arguments: what the console offers)."""
    if not args and not kwargs:
        print(BANNER, end="")
        return
    builtins.help(*args, **kwargs)


def _no_exit(*args: Any) -> None:
    print("The console stays open; close the window to quit MILES.")


def _no_input(*args: Any) -> str:
    raise RuntimeError("input() is not available here; set variables instead")


def namespace(app: App) -> dict[str, Any]:
    """Names ready in the console and in scripts.

    Parameters
    ----------
    app : App
        The running application.

    Returns
    -------
    dict
        A namespace to run code in, with ``doc`` and ``mol`` for the current
        tab (see :func:`refresh`).
    """
    helpers = Helpers(app)
    names: dict[str, Any] = {
        "__name__": "__console__", "__doc__": None,
        "app": app, "np": np, "miles": miles, "Molecule": Molecule,
        "help": _help, "exit": _no_exit, "quit": _no_exit, "input": _no_input,
        _BOUND: (None, None),
    }
    for name in ("load", "select", "count", "define", "show", "fit", "turn", "png", "run", "dist", "angle",
                 "dihedral"):
        names[name] = getattr(helpers, name)
    refresh(names, app)
    return names


def refresh(names: dict[str, Any], app: App) -> None:
    """Point ``doc`` and ``mol`` at the current tab, unless the user bound them to something else.

    Parameters
    ----------
    names : dict
        A namespace made by :func:`namespace`.
    app : App
        The running application.
    """
    old_doc, old_mol = names.get(_BOUND, (None, None))
    doc = app.current
    mol = doc.molecule if doc is not None else None
    if "doc" not in names or names["doc"] is old_doc:
        names["doc"] = doc
    if "mol" not in names or names["mol"] is old_mol:
        names["mol"] = mol
    names[_BOUND] = (names["doc"] if names["doc"] is doc else old_doc,
                     names["mol"] if names["mol"] is mol else old_mol)


def snapshot(app: App) -> dict[int, tuple[int, int]]:
    """Revision and a fingerprint of every open structure, to notice changes made behind the API's back.

    Returns
    -------
    dict
        For each document (by id), its revision and a hash of its atoms,
        bonds and coordinates.
    """
    return {id(view.canvas.doc): (view.canvas.doc.revision, _fingerprint(view.molecule))
            for view in app.documents}


def _fingerprint(mol: Molecule) -> int:
    return hash((len(mol), mol.numbers.tobytes(), mol.bond_pairs.tobytes(), mol.bond_orders.tobytes(),
                 mol.coords.tobytes()))


def settle(app: App, before: dict[int, tuple[int, int]]) -> None:
    """Redraw after code ran: structures changed without an edit are announced, so the view follows.

    Changes made directly (``mol.coords = mol.coords + 1``) are not undoable steps;
    :meth:`miles.api.DocumentView.edit` makes them so.

    Parameters
    ----------
    app : App
        The running application.
    before : dict
        What :func:`snapshot` returned before the code ran.
    """
    from miles.ui.document import Change

    for view in app.documents:
        doc = view.canvas.doc
        old = before.get(id(doc))
        if old is not None and old[0] == doc.revision and old[1] != _fingerprint(view.molecule):
            doc.notify(Change.STRUCTURE)
        view.canvas.redraw()
