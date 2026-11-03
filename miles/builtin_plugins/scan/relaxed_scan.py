"""Relaxed scans: step a distance, angle or dihedral and relax the rest with a force field.

At each step the coordinate is set by moving the smaller side of the
molecule (see :func:`miles.set_dihedral <miles.core.internal.set_dihedral>`), then everything is
minimised with Open Babel's MMFF94 (UFF when MMFF94 has no parameters)
while a stiff restraint holds the coordinate at its target. Each step
starts from the geometry of the one before, so the scan follows a
continuous path.
"""

from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
from openbabel import openbabel as ob

from miles import (
    FORCE_FIELDS,
    ChemError,
    Molecule,
    angle,
    dihedral,
    openbabel_coordinates,
    set_angle,
    set_bond_length,
    set_dihedral,
    to_openbabel,
)

RESTRAINT = 1e6          # restraint force constant; with Open Babel's default dihedrals drift by degrees
DESCENT_STEPS = 300      # steepest descent first: robust far from the minimum
GRADIENT_STEPS = 500     # then conjugate gradients
CONVERGENCE = 1e-8


class Coordinate(StrEnum):
    """Kind of internal coordinate, set by the number of atoms."""

    DISTANCE = "distance"   # 2 atoms, angstrom
    ANGLE = "angle"         # 3 atoms, degrees
    DIHEDRAL = "dihedral"   # 4 atoms, degrees


KINDS = {2: Coordinate.DISTANCE, 3: Coordinate.ANGLE, 4: Coordinate.DIHEDRAL}


@dataclass
class ScanResult:
    """Geometries and energies along a scan.

    Attributes
    ----------
    kind : Coordinate
        What was scanned.
    atoms : tuple of int
        The atoms defining the coordinate.
    values : numpy.ndarray
        Coordinate measured on each relaxed geometry, shape ``(n_points,)``;
        angstrom for distances, degrees otherwise.
    energies : numpy.ndarray
        Force-field energy of each geometry in kcal/mol, shape ``(n_points,)``.
    frames : numpy.ndarray
        Coordinates, shape ``(n_points, n_atoms, 3)``.
    force_field : str
        Open Babel force field used, e.g. ``"MMFF94"``.
    """

    kind: Coordinate
    atoms: tuple[int, ...]
    values: np.ndarray
    energies: np.ndarray
    frames: np.ndarray
    force_field: str


def coordinate_kind(atoms: Sequence[int]) -> Coordinate:
    """Kind of coordinate defined by a list of atoms.

    Parameters
    ----------
    atoms : sequence of int
        Two, three or four distinct atom indices.

    Returns
    -------
    Coordinate
        Distance, angle or dihedral.

    Raises
    ------
    ValueError
        If there are not 2-4 distinct atoms.
    """
    if len(atoms) not in KINDS or len(set(atoms)) != len(atoms):
        raise ValueError("Choose 2, 3 or 4 different atoms: a distance, an angle or a dihedral")
    return KINDS[len(atoms)]


def measure(coords: np.ndarray, atoms: Sequence[int]) -> float:
    """Measure the coordinate defined by ``atoms``.

    Parameters
    ----------
    coords : numpy.ndarray
        Atom positions, shape ``(n_atoms, 3)``.
    atoms : sequence of int
        Two, three or four atoms.

    Returns
    -------
    float
        Distance in angstrom, or angle or dihedral in degrees (dihedrals in
        the range -180 to 180).
    """
    points = [np.asarray(coords[i], dtype=float) for i in atoms]
    match coordinate_kind(atoms):
        case Coordinate.DISTANCE:
            return float(np.linalg.norm(points[1] - points[0]))
        case Coordinate.ANGLE:
            return float(angle(*points))
        case _:
            return float(dihedral(*points))


def _set(mol: Molecule, atoms: Sequence[int], value: float) -> None:
    match coordinate_kind(atoms):
        case Coordinate.DISTANCE:
            set_bond_length(mol, atoms[0], atoms[1], value)
        case Coordinate.ANGLE:
            set_angle(mol, atoms[0], atoms[1], atoms[2], value)
        case _:
            set_dihedral(mol, atoms[0], atoms[1], atoms[2], atoms[3], value)


def _restraint(atoms: Sequence[int], value: float) -> "ob.OBFFConstraints":
    constraints = ob.OBFFConstraints()
    constraints.SetFactor(RESTRAINT)
    ids = [int(i) + 1 for i in atoms]   # Open Babel counts atoms from 1
    match len(ids):
        case 2:
            constraints.AddDistanceConstraint(ids[0], ids[1], value)
        case 3:
            constraints.AddAngleConstraint(ids[0], ids[1], ids[2], value)
        case _:
            constraints.AddTorsionConstraint(ids[0], ids[1], ids[2], ids[3], value)
    return constraints


def _force_field(obmol: "ob.OBMol", names: Sequence[str]) -> tuple[str, "ob.OBForceField"]:
    for name in names:
        ff = ob.OBForceField.FindForceField(name)
        if ff is not None and ff.Setup(obmol):
            return name, ff
    raise ChemError("No force field parameters for this structure")


def scan(mol: Molecule, atoms: Sequence[int], start: float, stop: float, steps: int,
         force_fields: Sequence[str] = FORCE_FIELDS) -> ScanResult:
    """Relaxed scan of a distance, angle or dihedral.

    Parameters
    ----------
    mol : Molecule
        Structure (not modified).
    atoms : sequence of int
        Two atoms for a distance, three for an angle (vertex in the middle)
        or four for a dihedral.
    start, stop : float
        First and last value; angstrom for distances, degrees otherwise.
    steps : int
        Number of steps; the scan has ``steps + 1`` points.
    force_fields : sequence of str, optional
        Open Babel force fields to try, in order.

    Returns
    -------
    ScanResult
        One relaxed geometry and energy per point.

    Raises
    ------
    ValueError
        If the atoms do not define a coordinate, ``steps`` is less than 1,
        or a distance or angle is out of range.
    ChemError
        If no force field has parameters for the structure.
    """
    kind = coordinate_kind(atoms)
    if any(not 0 <= i < len(mol) for i in atoms):
        raise ValueError("Atom index out of range")
    if steps < 1:
        raise ValueError("A scan needs at least one step")
    if kind == Coordinate.DISTANCE and min(start, stop) <= 0.0:
        raise ValueError("Distances must be positive")
    if kind == Coordinate.ANGLE and not (0.0 < start < 180.0 and 0.0 < stop < 180.0):
        raise ValueError("Angles must lie between 0 and 180 degrees")
    work = mol.copy()
    obmol = to_openbabel(work)
    name, ff = _force_field(obmol, force_fields)
    targets = np.linspace(start, stop, steps + 1)
    values, energies, frames = [], [], []
    for target in targets:
        _set(work, atoms, float(target))
        for atom, position in zip(ob.OBMolAtomIter(obmol), work.coords, strict=True):
            atom.SetVector(*(float(x) for x in position))
        ff.Setup(obmol, _restraint(atoms, float(target)))
        ff.SteepestDescent(DESCENT_STEPS, CONVERGENCE)
        ff.ConjugateGradients(GRADIENT_STEPS, CONVERGENCE)
        ff.GetCoordinates(obmol)
        ff.Setup(obmol, ob.OBFFConstraints())   # the energy without the restraint
        energies.append(ff.Energy())
        work.coords = openbabel_coordinates(obmol)
        values.append(measure(work.coords, atoms))
        frames.append(work.coords.copy())
    return ScanResult(kind, tuple(int(i) for i in atoms), np.array(values), np.array(energies),
                      np.array(frames), name)
