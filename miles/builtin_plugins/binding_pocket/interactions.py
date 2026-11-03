"""Interactions between a ligand and its binding pocket: hydrogen bonds, salt bridges and pi stacking.

The criteria are the usual geometric ones: hydrogen bonds as
:class:`~miles.builtin_plugins.shared.hbonds.HydrogenBondFinder` finds them
(the definition of View > Hydrogen Bonds), salt bridges from the distance
between the centres of oppositely charged groups, and pi stacking from the
distance and angle between aromatic rings.
"""

import math
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum

import numpy as np

from miles import AMINO_ACIDS, Geometry, Molecule, near_mask, residue_ids, water_mask
from miles.builtin_plugins.shared.hbonds import HBOND_ANGLE, HBOND_MAX, HydrogenBondFinder

SALT_BRIDGE_MAX = 4.0    # between group centres
STACKING_MAX = 5.5       # parallel rings, between centroids
T_STACKING_MAX = 6.0     # perpendicular rings
PARALLEL_ANGLE = 30.0    # largest angle between ring normals for parallel stacking
T_ANGLE = 60.0           # smallest angle between ring normals for T-shaped stacking
RING_FLATNESS = 0.1      # largest out-of-plane deviation of an aromatic ring, angstrom

Group = tuple[int, ...]


class InteractionKind(StrEnum):
    """Kinds of non-covalent interaction."""

    HYDROGEN_BOND = "hydrogen bond"
    SALT_BRIDGE = "salt bridge"
    PI_STACKING = "pi stacking"


@dataclass(frozen=True, slots=True)
class Interaction:
    """One interaction between the ligand and its pocket.

    Attributes
    ----------
    kind : InteractionKind
        What kind of interaction.
    first : tuple of int
        Ligand atoms whose centroid is one end: one atom, a charged group or a ring.
    second : tuple of int
        Pocket atoms whose centroid is the other end.
    distance : float
        Distance between the two centroids when the interaction was found, in angstrom.
    """

    kind: InteractionKind
    first: Group
    second: Group
    distance: float


# -- ligand and pocket -----------------------------------------------------------------

def ligand_atoms(mol: Molecule, i: int) -> list[int]:
    """Atoms of the ligand containing atom ``i``.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i : int
        Any atom of the ligand.

    Returns
    -------
    list of int
        The atom's residue if it has residue data, otherwise its connected molecule.
    """
    ids = residue_ids(mol)
    return np.flatnonzero(ids == ids[i]).tolist()


def pocket_atoms(mol: Molecule, ligand: Sequence[int], cutoff: float = 5.0) -> list[int]:
    """Atoms of every residue with an atom close to a ligand.

    Parameters
    ----------
    mol : Molecule
        Structure.
    ligand : sequence of int
        Ligand atoms.
    cutoff : float, optional
        Distance from any ligand atom, in angstrom.

    Returns
    -------
    list of int
        Whole residues (or whole molecules, for atoms without residue data),
        without the ligand and without water.
    """
    lig = sorted({int(i) for i in ligand})
    if not lig:
        return []
    near = near_mask(mol.coords, mol.coords[lig], cutoff)
    near[lig] = False
    near &= ~water_mask(mol)
    ids = residue_ids(mol)
    chosen = np.isin(ids, np.unique(ids[near]))
    chosen[lig] = False
    return np.flatnonzero(chosen).tolist()


# -- hydrogen bonds ----------------------------------------------------------------------

def _hydrogen_bonds(mol: Molecule, ligand: list[int], pocket: list[int], max_distance: float,
                    min_angle: float) -> list[Interaction]:
    """Hydrogen bonds between ligand and pocket, as View > Hydrogen Bonds finds them."""
    inside = set(ligand)
    found = []
    for bond in HydrogenBondFinder(mol, ligand + pocket).find(mol.coords, max_distance, min_angle):
        if (bond.donor in inside) == (bond.acceptor in inside):
            continue   # both in the ligand, or both in the pocket
        first, second = (bond.donor, bond.acceptor) if bond.donor in inside else (bond.acceptor, bond.donor)
        found.append(Interaction(InteractionKind.HYDROGEN_BOND, (first,), (second,), bond.distance))
    return found

# -- salt bridges ------------------------------------------------------------------------

def _terminal_oxygens(mol: Molecule, center: int) -> list[int]:
    """Oxygens bonded only to ``center`` and carrying no hydrogen."""
    return [o for o in mol.neighbors(center)
            if mol.atoms[o].symbol == "O" and mol.degree(o) == 1]


def _topology_groups(mol: Molecule, atoms: list[int]) -> tuple[list[Group], list[Group]]:
    """Charged groups recognised from bonds: carboxylates, phosphates, sulfates, ammonium, amidines."""
    positive: list[Group] = []
    negative: list[Group] = []
    ring = mol.ring_bonds()
    for i in atoms:
        symbol = mol.atoms[i].symbol
        oxygens = [o for o in mol.neighbors(i) if mol.atoms[o].symbol == "O"]
        terminal = _terminal_oxygens(mol, i)
        carboxylate = symbol == "C" and len(oxygens) == 2 and len(terminal) == 2
        oxyanion = symbol in ("P", "S") and len(oxygens) >= 3 and bool(terminal)
        ammonium = symbol == "N" and mol.degree(i) == 4 and all(mol.bond_order(i, j) == 1.0
                                                                 for j in mol.neighbors(i))
        if carboxylate or oxyanion:
            negative.append(tuple(terminal))
        elif ammonium:
            positive.append((i,))
        elif symbol == "C" and not oxygens:
            nitrogens = [n for n in mol.neighbors(i) if mol.atoms[n].symbol == "N"]
            double = any(mol.bond_order(i, n) == 2.0 for n in nitrogens)
            if len(nitrogens) >= 2 and double and not any(mol.key(i, n) in ring for n in nitrogens):
                positive.append(tuple(nitrogens))
    return positive, negative


def _charged_groups(mol: Molecule, atoms: list[int], histidine: bool) -> tuple[list[Group], list[Group]]:
    """Positive and negative groups among some atoms: by name in amino acids, by bonds elsewhere."""
    residues: dict[tuple[str, int, str], dict[str, int]] = {}
    others = []
    for i in atoms:
        a = mol.atoms[i]
        if a.resname in AMINO_ACIDS:
            residues.setdefault((a.chain, a.resseq, a.resname), {})[a.name] = i
        else:
            others.append(i)
    positive, negative = _topology_groups(mol, others)
    for (_, _, resname), names in residues.items():

        def group(*wanted: str, names: dict[str, int] = names) -> Group:
            return tuple(names[n] for n in wanted if n in names)

        if resname == "LYS":
            positive.append(group("NZ"))
        elif resname == "ARG":
            positive.append(group("NE", "NH1", "NH2"))
        elif resname in ("HIS", "HIP", "HSP", "HID", "HIE", "HSD", "HSE"):
            ring_n = group("ND1", "NE2")
            protonated = resname in ("HIP", "HSP") or (
                len(ring_n) == 2 and all(mol.hydrogens_on(n) for n in ring_n))
            if histidine or protonated:
                positive.append(ring_n)
        elif resname == "ASP":
            negative.append(group("OD1", "OD2"))
        elif resname == "GLU":
            negative.append(group("OE1", "OE2"))
        if "OXT" in names:
            negative.append(group("O", "OXT"))
    return [g for g in positive if g], [g for g in negative if g]


def _salt_bridges(mol: Molecule, ligand: list[int], pocket: list[int], histidine: bool) -> list[Interaction]:
    lig_pos, lig_neg = _charged_groups(mol, ligand, histidine)
    poc_pos, poc_neg = _charged_groups(mol, pocket, histidine)
    found = []
    for mine, theirs in ((lig_pos, poc_neg), (lig_neg, poc_pos)):
        for g in mine:
            for h in theirs:
                d = float(np.linalg.norm(mol.centroid(g) - mol.centroid(h)))
                if d < SALT_BRIDGE_MAX:
                    found.append(Interaction(InteractionKind.SALT_BRIDGE, g, h, d))
    return found


# -- pi stacking ---------------------------------------------------------------------------

def aromatic_rings(mol: Molecule, atoms: Sequence[int]) -> list[Group]:
    """Flat five- and six-membered rings of trigonal or aromatic atoms among some atoms.

    Parameters
    ----------
    mol : Molecule
        Structure.
    atoms : sequence of int
        Atoms to look in; a ring must lie entirely within them.

    Returns
    -------
    list of tuple of int
        Each ring's atoms, sorted.
    """
    allowed = set(atoms)
    ring = mol.ring_bonds()
    adjacent = {i: [j for j in mol.neighbors(i) if j in allowed and mol.key(i, j) in ring] for i in allowed}
    cycles: set[Group] = set()
    for start in allowed:
        stack: list[tuple[int, Group]] = [(start, (start,))]
        while stack:
            v, path = stack.pop()
            for w in adjacent[v]:
                if w == start and len(path) in (5, 6):
                    cycles.add(tuple(sorted(path)))
                elif w > start and w not in path and len(path) < 6:
                    stack.append((w, (*path, w)))
    rings = []
    for cycle in sorted(cycles):
        trigonal = all(mol.atoms[i].geom == Geometry.TRIGONAL for i in cycle)
        aromatic = all(mol.bond_order(i, j) >= 1.5 for i in cycle for j in adjacent[i] if j in cycle)
        if (trigonal or aromatic) and _flatness(mol.coords[list(cycle)]) < RING_FLATNESS:
            rings.append(cycle)
    return rings


def _flatness(points: np.ndarray) -> float:
    centered = points - points.mean(axis=0)
    return float(np.abs(centered @ _normal(points)).max())


def _normal(points: np.ndarray) -> np.ndarray:
    return np.linalg.svd(points - points.mean(axis=0))[2][2]


def _pi_stacking(mol: Molecule, ligand: list[int], pocket: list[int]) -> list[Interaction]:
    found = []
    pocket_rings = [(r, mol.centroid(r), _normal(mol.coords[list(r)])) for r in aromatic_rings(mol, pocket)]
    for ring in aromatic_rings(mol, ligand):
        center, normal = mol.centroid(ring), _normal(mol.coords[list(ring)])
        for other, other_center, other_normal in pocket_rings:
            d = float(np.linalg.norm(center - other_center))
            tilt = math.degrees(math.acos(min(1.0, abs(float(normal @ other_normal)))))
            if (tilt < PARALLEL_ANGLE and d < STACKING_MAX) or (tilt >= T_ANGLE and d < T_STACKING_MAX):
                found.append(Interaction(InteractionKind.PI_STACKING, ring, other, d))
    return found


# -- all together ----------------------------------------------------------------------------

def find_interactions(mol: Molecule, ligand: Sequence[int], pocket: Sequence[int],
                      histidine: bool = False, max_distance: float = HBOND_MAX,
                      min_angle: float = HBOND_ANGLE) -> list[Interaction]:
    """Hydrogen bonds, salt bridges and pi stacking between a ligand and its pocket.

    Parameters
    ----------
    mol : Molecule
        Structure.
    ligand : sequence of int
        Ligand atoms.
    pocket : sequence of int
        Pocket atoms, e.g. from :func:`pocket_atoms`.
    histidine : bool, optional
        Count every histidine as positively charged; by default only those
        with both ring nitrogens protonated (or named HIP or HSP) count.
    max_distance, min_angle : float, optional
        Hydrogen bond cutoffs: largest donor-acceptor distance in angstrom
        and smallest D-H...A angle in degrees (the settings of View >
        Hydrogen Bonds).

    Returns
    -------
    list of Interaction
        Each interaction once, ligand side first, sorted by kind and distance.
    """
    lig = sorted({int(i) for i in ligand})
    poc = sorted({int(i) for i in pocket} - set(lig))
    found = (_hydrogen_bonds(mol, lig, poc, max_distance, min_angle) + _salt_bridges(mol, lig, poc, histidine)
             + _pi_stacking(mol, lig, poc))
    unique = {(f.kind, f.first, f.second): f for f in found}
    order = list(InteractionKind)
    return sorted(unique.values(), key=lambda f: (order.index(f.kind), f.distance))
