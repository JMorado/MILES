"""Fragment library: element fragments, rings, functional groups and biomolecules.

A :class:`Fragment` is a small molecule with an attachment point: a terminal
atom (normally a hydrogen) that is discarded when the fragment is bonded to
a structure. Rings and groups are built on first use and cached; amino
acids, nucleobases and ligands are read from ``data/fragments.json``.
"""

import json
import math
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum
from functools import cache
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np

from miles.core.building import AtomType, add_hydrogens, attach_fragment
from miles.core.chem import from_smiles
from miles.core.elements import bond_length, element
from miles.core.forcefield import clean_geometry
from miles.core.geometry import TEMPLATES, Geometry
from miles.core.molecule import Atom, Molecule, hydrogen_atom
from miles.core.perception import infer_atom_types

DATA_FILE = Path(__file__).with_name("data") / "fragments.json"


class Category(StrEnum):
    """Where a fragment comes from."""

    ELEMENT = "element"
    RING = "ring"
    GROUP = "group"
    BIO = "bio"
    LIGAND = "ligand"
    SMILES = "smiles"


@dataclass
class Fragment:
    """A molecule ready to be placed or attached.

    Attributes
    ----------
    key : str
        Unique identifier, e.g. ``"ring:benzene"``.
    name : str
        Name shown to the user.
    category : str
        A :class:`Category` value.
    mol : Molecule
        The fragment's atoms.
    attach : int or None
        Terminal atom replaced by the bond when attaching (None for bare atoms).
    symbol : str
        Element symbol, for element fragments.
    kind : AtomType or None
        Bonding type, for element fragments.
    """

    key: str
    name: str
    category: str
    mol: Molecule
    attach: int | None = None
    symbol: str = ""
    kind: AtomType | None = None

    @property
    def anchor(self) -> int:
        """Atom that bonds to the structure (the neighbour of ``attach``)."""
        return 0 if self.attach is None else self.mol.neighbors(self.attach)[0]

    def copy(self) -> "Fragment":
        """Return a copy with its own molecule."""
        return Fragment(self.key, self.name, self.category, self.mol.copy(), self.attach, self.symbol,
                        self.kind)


def molecule_fragment(mol: Molecule, name: str, key: str, attach: int | None = None,
                      category: str = Category.SMILES) -> Fragment:
    """Wrap a molecule as a fragment.

    Parameters
    ----------
    mol : Molecule
        Fragment atoms.
    name : str
        Name shown to the user.
    key : str
        Identifier, prefixed with the category.
    attach : int, optional
        Attachment atom; defaults to the first terminal hydrogen.
    category : str, optional
        A :class:`Category` value.

    Returns
    -------
    Fragment
        The fragment.
    """
    if attach is None:
        attach = next((i for i, a in enumerate(mol.atoms) if a.symbol == "H" and mol.degree(i) == 1), None)
    return Fragment(f"{category}:{key}", name, category, mol, attach)


def smiles_fragment(smiles: str, name: str, category: str = Category.SMILES) -> Fragment:
    """Build a fragment from SMILES (``*`` marks the attachment point).

    Raises
    ------
    ChemError
        If the string is not valid SMILES.
    """
    mol, attach = from_smiles(smiles)
    return molecule_fragment(mol, name, smiles, attach, category)


# -- element fragments --------------------------------------------------------------

_TRANSITION_METAL_TYPES = [
    AtomType("Octahedral", Geometry.OCTAHEDRAL, 6, 6),
    AtomType("Square planar", Geometry.SQUARE_PLANAR, 4, 4),
    AtomType("Tetrahedral", Geometry.TETRAHEDRAL, 4, 4),
    AtomType("Trigonal bipyramidal", Geometry.TBP, 5, 5),
    AtomType("Linear", Geometry.LINEAR, 2, 2),
]

_ELEMENT_TYPES = {
    "C": [
        AtomType("Tetrahedral", Geometry.TETRAHEDRAL, 4, 4),
        AtomType("Trigonal", Geometry.TRIGONAL, 3, 4),
        AtomType("Linear", Geometry.LINEAR, 2, 4),
    ],
    "N": [
        AtomType("Pyramidal", Geometry.TETRAHEDRAL, 3, 3),
        AtomType("Trigonal planar", Geometry.TRIGONAL, 3, 3),
        AtomType("Bent (double bond)", Geometry.TRIGONAL, 2, 3),
        AtomType("Linear (triple bond)", Geometry.LINEAR, 1, 3),
        AtomType("Tetrahedral (cation)", Geometry.TETRAHEDRAL, 4, 4),
    ],
    "O": [
        AtomType("Bent", Geometry.TETRAHEDRAL, 2, 2),
        AtomType("Double bond", Geometry.TRIGONAL, 1, 2),
        AtomType("Pyramidal (cation)", Geometry.TETRAHEDRAL, 3, 3),
    ],
    "S": [
        AtomType("Bent", Geometry.TETRAHEDRAL, 2, 2),
        AtomType("Double bond", Geometry.TRIGONAL, 1, 2),
        AtomType("Pyramidal (sulfoxide)", Geometry.TETRAHEDRAL, 3, 4),
        AtomType("Tetrahedral (sulfone)", Geometry.TETRAHEDRAL, 4, 6),
        AtomType("Octahedral", Geometry.OCTAHEDRAL, 6, 6),
    ],
    "P": [
        AtomType("Pyramidal", Geometry.TETRAHEDRAL, 3, 3),
        AtomType("Tetrahedral (P=O)", Geometry.TETRAHEDRAL, 4, 5),
        AtomType("Trigonal bipyramidal", Geometry.TBP, 5, 5),
    ],
    "B": [
        AtomType("Trigonal", Geometry.TRIGONAL, 3, 3),
        AtomType("Tetrahedral (anion)", Geometry.TETRAHEDRAL, 4, 4),
    ],
}
_ELEMENT_TYPES["Se"] = _ELEMENT_TYPES["Te"] = _ELEMENT_TYPES["S"][:3]
_ELEMENT_TYPES["As"] = _ELEMENT_TYPES["Sb"] = _ELEMENT_TYPES["P"]

_MAIN_GROUP_TYPES = {
    1: [AtomType("Monovalent", Geometry.LINEAR, 1, 1)],
    17: [AtomType("Monovalent", Geometry.LINEAR, 1, 1)],
    2: [AtomType("Linear", Geometry.LINEAR, 2, 2)],
    13: [AtomType("Trigonal", Geometry.TRIGONAL, 3, 3), AtomType("Tetrahedral", Geometry.TETRAHEDRAL, 4, 4)],
    14: [AtomType("Tetrahedral", Geometry.TETRAHEDRAL, 4, 4)],
}


def element_types(symbol: str) -> list[AtomType]:
    """Return the bonding types offered for an element; the bare atom comes last.

    Parameters
    ----------
    symbol : str
        Element symbol.

    Returns
    -------
    list of AtomType
        Types in menu order (the first is the default).
    """
    if symbol == "H":
        return [AtomType("Atom", Geometry.LINEAR, 1, 1, bare=True)]
    e = element(symbol)
    if symbol in _ELEMENT_TYPES:
        types = list(_ELEMENT_TYPES[symbol])
    elif e.column == 18:
        types = []
    elif e.column in _MAIN_GROUP_TYPES and not e.is_f_block:
        types = list(_MAIN_GROUP_TYPES[e.column])
    else:
        types = list(_TRANSITION_METAL_TYPES)
    first = types[0] if types else AtomType("Atom", None, 0, 0)
    return [*types, AtomType("Atom", first.geom, first.nsites, first.valence, bare=True)]


def element_fragment(symbol: str, kind: AtomType) -> Fragment:
    """Build an atom with hydrogens on all its sites (or a bare atom).

    Parameters
    ----------
    symbol : str
        Element symbol.
    kind : AtomType
        Bonding type, one of :func:`element_types`.

    Returns
    -------
    Fragment
        Element fragment; it attaches through its first hydrogen.
    """
    mol = Molecule()
    mol.add_atom(Atom(symbol, kind.geom, kind.nsites, kind.valence), np.zeros(3))
    attach = None
    if not kind.bare and kind.nsites and kind.geom is not None:
        h_len = bond_length(symbol, "H")
        for direction in TEMPLATES[kind.geom][:kind.nsites]:
            mol.add_bond(0, mol.add_atom(hydrogen_atom(), direction * h_len), 1.0)
        attach = 1
    name = element(symbol).name
    label = f"{name} atom" if kind.bare else f"{name} {kind.label.lower()}"
    return Fragment(f"el:{symbol}:{kind.label}", label, Category.ELEMENT, mol, attach,
                    symbol=symbol, kind=kind)


def default_element_fragment(symbol: str) -> Fragment:
    """Return the element fragment of the default bonding type."""
    return element_fragment(symbol, element_types(symbol)[0])


def _element(symbol: str, label: str) -> Fragment:
    """Element fragment by type label, e.g. ``("O", "Double bond")``."""
    kind = next((k for k in element_types(symbol) if k.label == label), None)
    if kind is None:
        raise KeyError(f"{symbol} {label}")
    return element_fragment(symbol, kind)


# -- rings --------------------------------------------------------------------------

class _Spec(NamedTuple):
    """Element and bonding pattern of a ring atom."""

    symbol: str
    geom: Geometry
    nsites: int
    valence: float


_AROMATIC_C = _Spec("C", Geometry.TRIGONAL, 3, 4)
_SP3_C = _Spec("C", Geometry.TETRAHEDRAL, 4, 4)
_PYRIDINE_N = _Spec("N", Geometry.TRIGONAL, 2, 3)
_PYRROLE_N = _Spec("N", Geometry.TRIGONAL, 3, 3)
_FURAN_O = _Spec("O", Geometry.TRIGONAL, 2, 2)
_THIOPHENE_S = _Spec("S", Geometry.TRIGONAL, 2, 2)
_AROMATIC_BOND = 1.40

Built = tuple[Molecule, int]
"""A built fragment molecule and its attachment atom."""


def _polygon(n: int, side: float) -> np.ndarray:
    """Regular planar polygon with one vertex on the +y axis."""
    radius = side / (2 * math.sin(math.pi / n))
    angles = [2 * math.pi * k / n + math.pi / 2 for k in range(n)]
    return np.array([[radius * math.cos(a), radius * math.sin(a), 0.0] for a in angles])


def _envelope() -> np.ndarray:
    """Cyclopentane envelope: a pentagon with one atom lifted out of the plane."""
    pts = _polygon(5, 1.54)
    pts[0, 2] = 0.55
    return pts


def _chair(n: int = 6, bond: float = 1.53, theta: float = 111.0) -> np.ndarray:
    """Cyclohexane chair with the given bond length and C-C-C angle.

    The ring radius ``r`` and pucker ``h`` solve ``r**2 + 4 h**2 = bond**2``
    and ``cos(theta) = (4 h**2 - r**2 / 2) / bond**2``.
    """
    c = math.cos(math.radians(theta))
    r2 = (bond ** 2 * (1 - c)) / 1.5
    h = math.sqrt(max(bond ** 2 - r2, 0.0)) / 2
    r = math.sqrt(r2)
    return np.array([[r * math.cos(math.pi / 3 * k), r * math.sin(math.pi / 3 * k), h if k % 2 == 0 else -h]
                     for k in range(n)])


def _finish(mol: Molecule, relax: bool = True) -> Built:
    """Add hydrogens, relax, and attach through the first hydrogen."""
    add_hydrogens(mol)
    if relax:
        clean_geometry(mol)
    return mol, next(i for i, a in enumerate(mol.atoms) if a.symbol == "H")


def _ring(coords: np.ndarray, specs: list[_Spec], orders: list[float]) -> Built:
    """Monocycle from atom positions, atom specs and the orders of bonds k-(k+1)."""
    mol = Molecule()
    for position, spec in zip(coords, specs, strict=True):
        mol.add_atom(Atom(spec.symbol, spec.geom, spec.nsites, spec.valence), position)
    for k, order in enumerate(orders):
        mol.add_bond(k, (k + 1) % len(specs), order)
    return _finish(mol)


def _six(specs: list[_Spec]) -> Built:
    """Aromatic six-membered ring."""
    return _ring(_polygon(6, 1.39), specs, [1.5] * 6)


def _five(side: float, first: _Spec, *rest: _Spec) -> Built:
    """Five-membered ring with two double bonds opposite ``first``."""
    return _ring(_polygon(5, side), [first, *rest], [1, 2, 1, 2, 1])


def _honeycomb(centers: list[tuple[float, float]]) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Vertices and edges of fused hexagons centred at the given points."""
    side = _AROMATIC_BOND
    points: list[np.ndarray] = []
    for cx, cy in centers:
        for m in range(6):
            a = math.radians(30 + 60 * m)
            p = np.array([cx + side * math.cos(a), cy + side * math.sin(a), 0.0])
            if all(np.linalg.norm(p - q) >= 0.1 for q in points):
                points.append(p)
    coords = np.array(points)
    edges = [(i, j) for i in range(len(points)) for j in range(i + 1, len(points))
             if abs(np.linalg.norm(coords[i] - coords[j]) - side) < 0.05]
    return coords, edges


def _polycyclic(centers: list[tuple[float, float]]) -> Built:
    """Fused benzenoid hydrocarbon."""
    coords, edges = _honeycomb(centers)
    mol = Molecule()
    for p in coords:
        mol.add_atom(Atom("C", Geometry.TRIGONAL, 3, 4), p)
    for i, j in edges:
        mol.add_bond(i, j, 1.5)
    return _finish(mol)


_STEP = math.sqrt(3) * _AROMATIC_BOND  # distance between centres of fused hexagons
_C, _N = _AROMATIC_C, _PYRIDINE_N

_RINGS: dict[str, tuple[str, Callable[[], Built]]] = {
    "benzene": ("Benzene", lambda: _six([_C] * 6)),
    "pyridine": ("Pyridine", lambda: _six([_N] + [_C] * 5)),
    "pyrimidine": ("Pyrimidine", lambda: _six([_N, _C, _N, _C, _C, _C])),
    "pyrazine": ("Pyrazine", lambda: _six([_N, _C, _C, _N, _C, _C])),
    "naphthalene": ("Naphthalene", lambda: _polycyclic([(0, 0), (_STEP, 0)])),
    "anthracene": ("Anthracene", lambda: _polycyclic([(0, 0), (_STEP, 0), (2 * _STEP, 0)])),
    "phenanthrene": ("Phenanthrene", lambda: _polycyclic([(0, 0), (_STEP, 0), (1.5 * _STEP, 1.5 * 1.40)])),
    "cyclopropane": ("Cyclopropane", lambda: _ring(_polygon(3, 1.51), [_SP3_C] * 3, [1] * 3)),
    "cyclobutane": ("Cyclobutane", lambda: _ring(_polygon(4, 1.55), [_SP3_C] * 4, [1] * 4)),
    "cyclopentane": ("Cyclopentane", lambda: _ring(_envelope(), [_SP3_C] * 5, [1] * 5)),
    "cyclohexane": ("Cyclohexane (chair)", lambda: _ring(_chair(), [_SP3_C] * 6, [1] * 6)),
    "cyclopentadiene": ("Cyclopentadiene", lambda: _five(1.45, _SP3_C, _C, _C, _C, _C)),
    "pyrrole": ("Pyrrole", lambda: _five(1.39, _PYRROLE_N, _C, _C, _C, _C)),
    "furan": ("Furan", lambda: _five(1.39, _FURAN_O, _C, _C, _C, _C)),
    "thiophene": ("Thiophene", lambda: _five(1.45, _THIOPHENE_S, _C, _C, _C, _C)),
    "imidazole": ("Imidazole", lambda: _five(1.37, _PYRROLE_N, _C, _N, _C, _C)),
}


# -- functional groups ----------------------------------------------------------------

class _GroupBuilder:
    """Grows a group from element fragments; hydrogen 1 stays the attachment point."""

    def __init__(self, symbol: str, label: str) -> None:
        frag = _element(symbol, label)
        assert frag.attach is not None
        self.mol = frag.mol
        self.attach: int = frag.attach

    def _free_hydrogen(self, on: int) -> int:
        return next(h for h in self.mol.hydrogens_on(on) if h != self.attach)

    def grow(self, on: int, symbol: str, label: str) -> int:
        """Replace a hydrogen on atom ``on`` by an element fragment; return the new atom."""
        piece = _element(symbol, label)
        return self.graft(on, piece)

    def graft(self, on: int, frag: Fragment) -> int:
        """Replace a hydrogen on atom ``on`` by a fragment; return its anchor atom."""
        target = self._free_hydrogen(on)
        attach_fragment(self.mol, target, frag.mol, frag.attach)
        return target

    def done(self) -> Built:
        """Relax the group and return it with its attachment atom."""
        clean_geometry(self.mol)
        return self.mol, self.attach


Step = tuple[int, str, str]
"""Growth step: (index into the atoms grown so far, element, type label); 0 is the root."""

_METHYL: Step = (0, "C", "Tetrahedral")
_CARBONYL_O: Step = (0, "O", "Double bond")

_RECIPES: dict[str, tuple[str, str, str, list[Step]]] = {
    "methyl": ("Methyl", "C", "Tetrahedral", []),
    "ethyl": ("Ethyl", "C", "Tetrahedral", [_METHYL]),
    "isopropyl": ("Isopropyl", "C", "Tetrahedral", [_METHYL] * 2),
    "tbutyl": ("tert-Butyl", "C", "Tetrahedral", [_METHYL] * 3),
    "vinyl": ("Vinyl", "C", "Trigonal", [(0, "C", "Trigonal")]),
    "ethynyl": ("Ethynyl", "C", "Linear", [(0, "C", "Linear")]),
    "hydroxyl": ("Hydroxyl", "O", "Bent", []),
    "methoxy": ("Methoxy", "O", "Bent", [_METHYL]),
    "thiol": ("Thiol", "S", "Bent", []),
    "amino": ("Amino", "N", "Pyramidal", []),
    "dimethylamino": ("Dimethylamino", "N", "Pyramidal", [_METHYL] * 2),
    "cyano": ("Cyano", "C", "Linear", [(0, "N", "Linear (triple bond)")]),
    "formyl": ("Formyl", "C", "Trigonal", [_CARBONYL_O]),
    "acetyl": ("Acetyl", "C", "Trigonal", [_CARBONYL_O, _METHYL]),
    "carboxyl": ("Carboxyl", "C", "Trigonal", [_CARBONYL_O, (0, "O", "Bent")]),
    "ester": ("Methyl ester", "C", "Trigonal", [_CARBONYL_O, (0, "O", "Bent"), (2, "C", "Tetrahedral")]),
    "amide": ("Amide", "C", "Trigonal", [_CARBONYL_O, (0, "N", "Trigonal planar")]),
    "cf3": ("Trifluoromethyl", "C", "Tetrahedral", [(0, "F", "Monovalent")] * 3),
    "sulfonic": ("Sulfonic acid", "S", "Tetrahedral (sulfone)", [_CARBONYL_O, _CARBONYL_O, (0, "O", "Bent")]),
    "tms": ("Trimethylsilyl", "Si", "Tetrahedral", [_METHYL] * 3),
}


def _from_recipe(symbol: str, label: str, steps: list[Step]) -> Built:
    builder = _GroupBuilder(symbol, label)
    grown = [0]
    for parent, step_symbol, step_label in steps:
        grown.append(builder.grow(grown[parent], step_symbol, step_label))
    return builder.done()


def _phenyl() -> Built:
    frag = ring_fragment("benzene")
    assert frag.attach is not None
    return frag.mol, frag.attach


def _benzyl() -> Built:
    builder = _GroupBuilder("C", "Tetrahedral")
    builder.graft(0, ring_fragment("benzene"))
    return builder.done()


def _nitro() -> Built:
    """Nitro group with two aromatic-order N-O bonds."""
    builder = _GroupBuilder("N", "Trigonal planar")
    builder.mol.atoms[0].valence = 4
    oxygens = [builder.grow(0, "O", "Double bond") for _ in range(2)]
    for o in oxygens:
        builder.mol.add_bond(0, o, 1.5)
    return builder.done()


_SPECIAL_GROUPS: dict[str, tuple[str, Callable[[], Built]]] = {
    "phenyl": ("Phenyl", _phenyl),
    "benzyl": ("Benzyl", _benzyl),
    "nitro": ("Nitro", _nitro),
}

_GROUP_ORDER = ["methyl", "ethyl", "isopropyl", "tbutyl", "vinyl", "ethynyl", "phenyl", "benzyl",
                "hydroxyl", "methoxy", "thiol", "amino", "dimethylamino", "nitro", "cyano", "formyl",
                "acetyl", "carboxyl", "ester", "amide", "cf3", "sulfonic", "tms"]


def _group(key: str) -> tuple[str, Callable[[], Built]]:
    if key in _SPECIAL_GROUPS:
        return _SPECIAL_GROUPS[key]
    name, symbol, label, steps = _RECIPES[key]
    return name, lambda: _from_recipe(symbol, label, steps)


# -- public access ----------------------------------------------------------------------

@cache
def _built_ring(key: str) -> Fragment:
    name, build = _RINGS[key]
    mol, attach = build()
    return Fragment(f"ring:{key}", name, Category.RING, mol, attach)


@cache
def _built_group(key: str) -> Fragment:
    name, build = _group(key)
    mol, attach = build()
    return Fragment(f"group:{key}", name, Category.GROUP, mol, attach)


def ring_fragment(key: str) -> Fragment:
    """Return a copy of a ring fragment, e.g. ``"benzene"``."""
    return _built_ring(key).copy()


def group_fragment(key: str) -> Fragment:
    """Return a copy of a functional group fragment, e.g. ``"ester"``."""
    return _built_group(key).copy()


def all_ring_fragments() -> list[Fragment]:
    """Return every ring fragment, in menu order."""
    return [ring_fragment(key) for key in _RINGS]


def all_group_fragments() -> list[Fragment]:
    """Return every functional group fragment, in menu order."""
    return [group_fragment(key) for key in _GROUP_ORDER]


# -- biomolecules and ligands (pre-built, shipped as data) -------------------------------

@cache
def _library_data() -> dict[str, list[dict[str, Any]]]:
    with open(DATA_FILE, encoding="utf-8") as fh:
        return json.load(fh)


def _from_entry(entry: dict[str, Any], category: str) -> Fragment:
    """Fragment from a data-file entry (``attach`` is the hydrogen that replaced ``*``)."""
    mol = Molecule.from_atoms([Atom(s) for s in entry["symbols"]], np.array(entry["coords"], dtype=float))
    for i, j, order in entry["bonds"]:
        mol.add_bond(i, j, order)
    for atom, name in zip(mol.atoms, entry.get("names", []), strict=False):
        atom.name = name
    infer_atom_types(mol, only_missing=False)
    return Fragment(f"{category}:{entry['key']}", entry["name"], category, mol, entry["attach"])


def library_fragments(category: str) -> list[Fragment]:
    """Return the bundled fragments of a category.

    Parameters
    ----------
    category : str
        ``"bio"`` (amino acids and nucleobases) or ``"ligand"``.

    Returns
    -------
    list of Fragment
        Fragments in file order.
    """
    return [_from_entry(entry, category) for entry in _library_data()[category]]


def read_library(path: str | Path) -> dict[str, list[Fragment]]:
    """Read a fragment library in the format of the bundled ``data/fragments.json``.

    The file maps a category to a list of entries. Each entry has a
    ``key``, a ``name``, ``symbols``, ``coords`` (angstrom), ``bonds`` as
    ``[i, j, order]`` and ``attach``, the index of the hydrogen replaced
    when the fragment is attached (or null); ``names`` (PDB atom names) is
    optional.

    Parameters
    ----------
    path : str or pathlib.Path
        JSON file.

    Returns
    -------
    dict of str to list of Fragment
        Fragments by category, in file order.

    Raises
    ------
    OSError
        If the file cannot be read.
    ValueError
        If it is not such a library.
    """
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict) or not all(isinstance(v, list) for v in data.values()):
        raise ValueError(f"{Path(path).name} is not a fragment library: it must map categories to lists")
    try:
        return {str(category): [_from_entry(entry, str(category)) for entry in entries]
                for category, entries in data.items()}
    except (KeyError, TypeError, IndexError) as exc:
        raise ValueError(f"{Path(path).name} has a fragment without {exc}") from exc


def amino_acid_template(code: str) -> Fragment:
    """Return a free amino acid with PDB atom names, by three-letter code (e.g. ``"TRP"``)."""
    entry = next((e for e in _library_data()[Category.BIO] if e.get("code") == code.upper()), None)
    if entry is None:
        raise KeyError(code)
    return _from_entry(entry, Category.BIO)


def amino_acid_codes() -> list[str]:
    """Return the three-letter codes of the bundled amino acids."""
    return [e["code"] for e in _library_data()[Category.BIO] if "code" in e]
