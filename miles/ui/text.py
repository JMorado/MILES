"""Text shown in the interface: units, atom names and short descriptions."""

import numpy as np

from miles.core.elements import element
from miles.core.geometry import GEOMETRY_NAMES
from miles.core.molecule import AtomLike, Molecule

DOT = "\u00b7"
DASH = "\u2013"
BULLET = "\u2022"
ANGSTROM = "\u00c5"
DEGREE = "\u00b0"
COMMAND = "\u2318"
SEPARATOR = f"  {DOT}  "

BOND_NAMES = {1.0: "single", 1.5: "aromatic", 2.0: "double", 3.0: "triple"}


def plural(n: int, word: str) -> str:
    """``"1 atom"``, ``"2 atoms"``."""
    return f"{n} {word}{'' if n == 1 else 's'}"


def subscript_digits(text: str) -> str:
    """Write the digits of a formula as subscripts, e.g. ``C6H6``."""
    return "".join(chr(0x2080 + int(ch)) if ch.isdigit() else ch for ch in text)


def atom_label(mol: Molecule, i: int) -> str:
    """Element symbol and 1-based number, e.g. ``C12``."""
    return f"{mol.atoms[i].symbol}{i + 1}"


def residue_label(atom: AtomLike) -> str:
    """Residue name, number and chain, e.g. ``ALA 12 A``."""
    return f"{atom.resname} {atom.resseq}{' ' + atom.chain if atom.chain else ''}"


def coordinates(position: np.ndarray, decimals: int = 3) -> str:
    """Cartesian coordinates separated by wide spaces."""
    return "  ".join(f"{v:.{decimals}f}" for v in position)


def describe_atom(mol: Molecule, i: int) -> str:
    """One line about an atom for the status bar: name, element, geometry and position."""
    a = mol.atoms[i]
    where = coordinates(mol.coords[i])
    if a.resname:
        return f"{residue_label(a)}  {a.name or a.symbol}  {element(a.symbol).name}   {where}"
    geometry = f"  {GEOMETRY_NAMES.get(a.geom, '')}" if a.geom and a.symbol != "H" else ""
    return f"{atom_label(mol, i)}  {element(a.symbol).name}{geometry}   {where}"


def describe_bond(mol: Molecule, i: int, j: int) -> str:
    """One line about a bond: its atoms, order and length."""
    length = float(np.linalg.norm(mol.coords[i] - mol.coords[j]))
    order = BOND_NAMES.get(mol.bond_order(i, j), "")
    return f"{atom_label(mol, i)}-{atom_label(mol, j)}  {order}  {length:.3f} {ANGSTROM}"
