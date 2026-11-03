"""Editing internal coordinates (bond lengths, angles, dihedrals) and rigid moves.

Changing an internal coordinate moves the part of the molecule on one side of
it. When both ends are joined through a ring only the end atom (with its
terminal neighbours) moves.
"""

import math
from collections.abc import Iterable
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

from miles.core.geometry import angle, any_perpendicular, dihedral, rotation_matrix, unit_vector
from miles.core.molecule import Molecule


class MoveMode(StrEnum):
    """Which part of the molecule moves when an internal coordinate changes."""

    LAST = "last"    # the group around the last atom
    FIRST = "first"  # the group around the first atom
    BOTH = "both"    # both groups, half each
    ATOM = "atom"    # only the last atom and its terminal neighbours


def _atom_only(mol: Molecule, i: int, exclude: Iterable[int]) -> set[int]:
    return {i, *mol.terminal_neighbors(i)} - set(exclude)


def _side(mol: Molecule, pivot: int, start: int, fallback: set[int] | None = None) -> set[int]:
    """Atoms moving with ``start`` when it turns about or slides away from ``pivot``."""
    group = mol.component(start, blocked={pivot})
    if any(n in group for n in mol.neighbors(pivot) if n != start):
        return _atom_only(mol, start, [pivot]) if fallback is None else set(fallback)
    return group


def translate_atoms(mol: Molecule, indices: Iterable[int], vector: ArrayLike) -> None:
    """Shift atoms by a vector."""
    idx = list(indices)
    if idx:
        mol.set_positions(idx, mol.coords[idx] + np.asarray(vector, dtype=float))


def rotate_atoms(mol: Molecule, indices: Iterable[int], rotation: np.ndarray,
                 center: ArrayLike | None = None) -> None:
    """Rotate atoms about a point.

    Parameters
    ----------
    mol : Molecule
        Structure.
    indices : iterable of int
        Atoms to rotate.
    rotation : numpy.ndarray
        Rotation matrix, shape ``(3, 3)``.
    center : array_like, optional
        Fixed point (default: centroid of the atoms).
    """
    idx = list(indices)
    if not idx:
        return
    c = mol.coords[idx].mean(axis=0) if center is None else np.asarray(center, dtype=float)
    mol.set_positions(idx, (mol.coords[idx] - c) @ np.asarray(rotation).T + c)


def _rotate(mol: Molecule, group: set[int], axis: np.ndarray, center: np.ndarray, phi: float) -> None:
    if group and abs(phi) > 1e-12:
        rotate_atoms(mol, group, rotation_matrix(axis, phi), center)


def set_bond_length(mol: Molecule, i: int, j: int, value: float, mode: MoveMode = MoveMode.LAST) -> None:
    """Set the distance between atoms ``i`` and ``j``.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i, j : int
        The two atoms (bonded or not).
    value : float
        New distance in angstrom.
    mode : MoveMode, optional
        Which side moves.
    """
    ri, rj = mol.coords[i], mol.coords[j]
    u = unit_vector(rj - ri)
    if not np.any(u):
        u = np.array([1.0, 0.0, 0.0])
    delta = value - float(np.linalg.norm(rj - ri))
    if mode == MoveMode.ATOM:
        translate_atoms(mol, _atom_only(mol, j, [i]), u * delta)
    elif mode == MoveMode.FIRST:
        translate_atoms(mol, _side(mol, j, i), -u * delta)
    elif mode == MoveMode.BOTH:
        last, first = _side(mol, i, j), _side(mol, j, i)
        translate_atoms(mol, last, u * delta / 2)
        translate_atoms(mol, first, -u * delta / 2)
    else:
        translate_atoms(mol, _side(mol, i, j), u * delta)


def set_angle(mol: Molecule, i: int, j: int, k: int, value: float, mode: MoveMode = MoveMode.LAST) -> None:
    """Set the angle i-j-k (vertex ``j``) by turning one side about the angle's normal.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i, j, k : int
        The three atoms; ``j`` is the vertex.
    value : float
        New angle in degrees.
    mode : MoveMode, optional
        Which side moves.
    """
    ri, rj, rk = mol.coords[i], mol.coords[j], mol.coords[k]
    normal = unit_vector(np.cross(ri - rj, rk - rj))
    if np.linalg.norm(normal) < 1e-8:
        normal = any_perpendicular(rk - rj)
    delta = math.radians(value - angle(ri, rj, rk))
    if mode == MoveMode.ATOM:
        _rotate(mol, _atom_only(mol, k, [j]), normal, rj, delta)
    elif mode == MoveMode.FIRST:
        _rotate(mol, _side(mol, j, i), normal, rj, -delta)
    elif mode == MoveMode.BOTH:
        _rotate(mol, _side(mol, j, k), normal, rj, delta / 2)
        _rotate(mol, _side(mol, j, i), normal, rj, -delta / 2)
    else:
        _rotate(mol, _side(mol, j, k), normal, rj, delta)


def set_dihedral(mol: Molecule, i: int, j: int, k: int, l: int, value: float,  # noqa: E741
                 mode: MoveMode = MoveMode.LAST) -> None:
    """Set the dihedral i-j-k-l by turning one side about the j-k bond.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i, j, k, l : int
        The four atoms.
    value : float
        New dihedral in degrees.
    mode : MoveMode, optional
        Which side turns.
    """
    rj, rk = mol.coords[j], mol.coords[k]
    axis = unit_vector(rk - rj)
    current = dihedral(mol.coords[i], rj, rk, mol.coords[l])
    delta = math.radians((value - current + 180.0) % 360.0 - 180.0)
    if mode == MoveMode.ATOM:
        _rotate(mol, _atom_only(mol, l, [k]), axis, rk, delta)
        return
    last = _side(mol, j, k, fallback=_atom_only(mol, l, [k]))
    first = _side(mol, k, j, fallback=_atom_only(mol, i, [j]))
    if mode == MoveMode.FIRST:
        _rotate(mol, first, axis, rj, -delta)
    elif mode == MoveMode.BOTH:
        _rotate(mol, last, axis, rk, delta / 2)
        _rotate(mol, first, axis, rj, -delta / 2)
    else:
        _rotate(mol, last, axis, rk, delta)
