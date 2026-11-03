"""Backbone dihedrals of proteins (phi and psi) and approximate Ramachandran regions.

Residues are linked through their peptide bonds: the carbonyl carbon of one
residue must lie within :data:`PEPTIDE_BOND` of the amide nitrogen of the
next. Residue numbers are not trusted, so gaps and chain breaks leave the
dihedrals undefined (``None``).

The regions are simple ellipses in the (phi, psi) plane that follow the
well-known shapes of the plot: right-handed helix, beta strand and
polyproline II, and left-handed helix. They are approximate, meant for a
quick look at a structure, and do not reproduce the contours of validation
programs.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from typing import NamedTuple

import numpy as np

from miles import AMINO_ACIDS, Molecule, residues
from miles.plugin.widgets import Table, format_number

PEPTIDE_BOND = 2.0   # angstrom; longest C-N distance still taken as a peptide bond


class ResidueKind(StrEnum):
    """Residue classes with different allowed backbone conformations."""

    GENERAL = "general"
    GLYCINE = "glycine"
    PROLINE = "proline"
    PRE_PROLINE = "pre-proline"   # the residue before a proline


class Region(StrEnum):
    """Where a residue falls on the Ramachandran plot."""

    FAVOURED = "favoured"
    ALLOWED = "allowed"
    OUTLIER = "outlier"


class Ellipse(NamedTuple):
    """An ellipse in the (phi, psi) plane, in degrees, wrapping around at +-180.

    Attributes
    ----------
    phi, psi : float
        Centre.
    half_phi, half_psi : float
        Semi-axes along phi and psi.
    """

    phi: float
    psi: float
    half_phi: float
    half_psi: float

    def contains(self, phi: float, psi: float) -> bool:
        """Whether a point lies inside, taking the periodicity of the angles into account."""
        dphi = (phi - self.phi + 180.0) % 360.0 - 180.0
        dpsi = (psi - self.psi + 180.0) % 360.0 - 180.0
        return (dphi / self.half_phi) ** 2 + (dpsi / self.half_psi) ** 2 <= 1.0


# Right-handed helix: a core and the bridge region towards phi -120, psi 0.
_ALPHA_R = (Ellipse(-65, -40, 30, 32), Ellipse(-100, -10, 35, 30))
_BETA = Ellipse(-105, 145, 68, 60)   # beta strand and polyproline II
_ALPHA_L = Ellipse(60, 35, 20, 30)
_ALPHA_R_WIDE = Ellipse(-90, -20, 75, 65)
_BETA_WIDE = Ellipse(-105, 150, 95, 85)
_ALPHA_L_WIDE = Ellipse(65, 30, 35, 55)
_BETWEEN = Ellipse(-110, 50, 55, 40)   # joins the helix and strand regions on the left

REGIONS: dict[ResidueKind, tuple[tuple[Ellipse, ...], tuple[Ellipse, ...]]] = {
    ResidueKind.GENERAL: ((*_ALPHA_R, _BETA, _ALPHA_L), (_ALPHA_R_WIDE, _BETA_WIDE, _BETWEEN, _ALPHA_L_WIDE)),
    # glycine has no side chain, so its plot is symmetric through the origin
    ResidueKind.GLYCINE: ((Ellipse(-75, -30, 40, 40), Ellipse(75, 30, 40, 40),
                           Ellipse(-85, 175, 35, 30), Ellipse(85, -175, 35, 30), Ellipse(180, 180, 30, 30)),
                          (Ellipse(-80, -20, 65, 70), Ellipse(80, 20, 65, 70),
                           Ellipse(-120, 170, 75, 70), Ellipse(120, -170, 75, 70))),
    ResidueKind.PROLINE: ((Ellipse(-65, 145, 20, 35), Ellipse(-65, -30, 20, 30)),
                          (Ellipse(-70, 140, 35, 60), Ellipse(-75, -20, 35, 55))),
    ResidueKind.PRE_PROLINE: ((_BETA, Ellipse(-65, -40, 25, 28)),
                              (_BETA_WIDE, Ellipse(-80, -30, 50, 45), Ellipse(60, 40, 18, 22))),
}
"""Favoured and allowed ellipses of each residue class (approximate)."""


def classify(phi: float, psi: float, kind: ResidueKind = ResidueKind.GENERAL) -> Region:
    """Region of a backbone conformation.

    Parameters
    ----------
    phi, psi : float
        Backbone dihedrals in degrees.
    kind : ResidueKind, optional
        Residue class, which selects the regions.

    Returns
    -------
    Region
        Favoured, allowed or outlier.
    """
    favoured, allowed = REGIONS[kind]
    if any(e.contains(phi, psi) for e in favoured):
        return Region.FAVOURED
    if any(e.contains(phi, psi) for e in allowed):
        return Region.ALLOWED
    return Region.OUTLIER


@dataclass(frozen=True)
class BackboneResidue:
    """An amino acid residue with the atoms that define its phi and psi.

    Attributes
    ----------
    chain : str
        Chain identifier.
    resseq : int
        Residue number from the file.
    resname : str
        Three-letter residue name.
    atoms : tuple of int
        All atoms of the residue.
    kind : ResidueKind
        Class used for the regions.
    phi_atoms, psi_atoms : tuple of int or None
        C(i-1), N, CA, C and N, CA, C, N(i+1), or None at a chain end or break.
    """

    chain: str
    resseq: int
    resname: str
    atoms: tuple[int, ...]
    kind: ResidueKind
    phi_atoms: tuple[int, int, int, int] | None
    psi_atoms: tuple[int, int, int, int] | None


@dataclass(frozen=True)
class ResidueAngles:
    """Backbone dihedrals of one residue.

    Attributes
    ----------
    residue : BackboneResidue
        The residue.
    phi, psi : float or None
        Dihedrals in degrees; None where the neighbouring residue is missing.
    """

    residue: BackboneResidue
    phi: float | None
    psi: float | None

    @property
    def region(self) -> Region | None:
        """Region on the plot, or None when phi or psi is undefined."""
        if self.phi is None or self.psi is None:
            return None
        return classify(self.phi, self.psi, self.residue.kind)


def _backbone_atoms(mol: Molecule, atoms: list[int]) -> dict[str, int]:
    """Return the indices of N, CA and C in a residue (the first atom of each name)."""
    found: dict[str, int] = {}
    for i in atoms:
        name = mol.atoms[i].name
        if name in ("N", "CA", "C") and name not in found:
            found[name] = i
    return found


def _linked(mol: Molecule, c: int, n: int) -> bool:
    """Whether a carbonyl carbon and an amide nitrogen form a peptide bond."""
    if (min(c, n), max(c, n)) in mol.bonds:
        return True
    return float(np.linalg.norm(mol.coords[c] - mol.coords[n])) < PEPTIDE_BOND


def backbone_residues(mol: Molecule) -> list[BackboneResidue]:
    """Amino acid residues with N, CA and C, and the atoms of their phi and psi.

    Neighbours are found through peptide bonds: first the residue before or
    after in file order, then any bonded residue, so missing residues and
    chain breaks give undefined dihedrals.

    Parameters
    ----------
    mol : Molecule
        Structure with residue data (atom names, residue names and numbers).

    Returns
    -------
    list of BackboneResidue
        One entry per amino acid residue that has N, CA and C, in file order.
    """
    found = []
    for atoms in residues(mol):
        first = mol.atoms[atoms[0]]
        if first.resname not in AMINO_ACIDS:
            continue
        backbone = _backbone_atoms(mol, atoms)
        if len(backbone) == 3:
            found.append((first, atoms, backbone))
    owner = {i: k for k, (_, atoms, _) in enumerate(found) for i in atoms}

    def neighbour(k: int, step: int) -> int | None:
        """Residue joined to residue ``k`` before (step -1) or after (step +1) it."""
        own = found[k][2]
        j = k + step
        if 0 <= j < len(found):
            other = found[j][2]
            c, n = (other["C"], own["N"]) if step < 0 else (own["C"], other["N"])
            if _linked(mol, c, n):
                return j
        start = own["N"] if step < 0 else own["C"]
        partner = "C" if step < 0 else "N"
        for i in mol.neighbors(start):
            j = owner.get(i, -1)
            if j not in (-1, k) and found[j][2][partner] == i:
                return j
        return None

    result = []
    for k, (first, atoms, own) in enumerate(found):
        before, after = neighbour(k, -1), neighbour(k, +1)
        if first.resname == "GLY":
            kind = ResidueKind.GLYCINE
        elif first.resname == "PRO":
            kind = ResidueKind.PROLINE
        elif after is not None and found[after][0].resname == "PRO":
            kind = ResidueKind.PRE_PROLINE
        else:
            kind = ResidueKind.GENERAL
        phi = (found[before][2]["C"], own["N"], own["CA"], own["C"]) if before is not None else None
        psi = (own["N"], own["CA"], own["C"], found[after][2]["N"]) if after is not None else None
        result.append(BackboneResidue(first.chain, first.resseq, first.resname, tuple(atoms), kind, phi, psi))
    return result


def _dihedrals(coords: np.ndarray, quads: np.ndarray) -> np.ndarray:
    """Dihedral angles in degrees of many atom quadruples at once."""
    p, q, r, s = (coords[quads[:, k]] for k in range(4))
    b1 = r - q
    b1 /= np.linalg.norm(b1, axis=1, keepdims=True)
    v = (p - q) - np.sum((p - q) * b1, axis=1, keepdims=True) * b1
    w = (s - r) - np.sum((s - r) * b1, axis=1, keepdims=True) * b1
    return np.degrees(np.arctan2(np.sum(np.cross(b1, v) * w, axis=1), np.sum(v * w, axis=1)))


def backbone_angles(found: list[BackboneResidue], coords: np.ndarray) -> np.ndarray:
    """Phi and psi of residues for one set of coordinates.

    Parameters
    ----------
    found : list of BackboneResidue
        Residues from :func:`backbone_residues`.
    coords : numpy.ndarray
        Coordinates of the whole structure, shape ``(n_atoms, 3)``.

    Returns
    -------
    numpy.ndarray
        Shape ``(n_residues, 2)``: phi and psi in degrees, NaN where undefined.
    """
    angles = np.full((len(found), 2), np.nan)
    for column, field in enumerate(("phi_atoms", "psi_atoms")):
        rows = [k for k, r in enumerate(found) if getattr(r, field) is not None]
        if rows:
            quads = np.array([getattr(found[k], field) for k in rows])
            angles[rows, column] = _dihedrals(np.asarray(coords, dtype=float), quads)
    return angles


def ramachandran(mol: Molecule) -> list[ResidueAngles]:
    """Phi and psi of every amino acid residue of a structure.

    Parameters
    ----------
    mol : Molecule
        Structure with residue data.

    Returns
    -------
    list of ResidueAngles
        One entry per residue with N, CA and C atoms, in file order.
    """
    found = backbone_residues(mol)
    angles = backbone_angles(found, mol.coords)
    return [ResidueAngles(r, None if math.isnan(phi) else float(phi), None if math.isnan(psi) else float(psi))
            for r, (phi, psi) in zip(found, angles, strict=True)]


def ramachandran_table(residues: Sequence[BackboneResidue], angles: np.ndarray) -> Table:
    """Backbone dihedrals and region of every residue of a Ramachandran plot.

    Parameters
    ----------
    residues : sequence of BackboneResidue
        Residues, as the plot lists them.
    angles : numpy.ndarray
        Phi and psi in degrees, shape ``(n_residues, 2)``, NaN where undefined.

    Returns
    -------
    Table
        Chain, residue name and number, phi, psi and region (favoured,
        allowed or outlier; empty when an angle is missing).
    """
    rows = []
    for residue, (phi, psi) in zip(residues, np.asarray(angles, dtype=float), strict=True):
        missing = math.isnan(phi) or math.isnan(psi)
        region = "" if missing else str(classify(float(phi), float(psi), residue.kind))
        rows.append([residue.chain, residue.resname, str(residue.resseq), format_number(float(phi), 2),
                     format_number(float(psi), 2), region])
    return Table(["Chain", "Residue", "Number", "Phi (deg)", "Psi (deg)", "Region"], rows)
