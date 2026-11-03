"""Chemical perception: bonds from distances, bond orders and atom types.

Used for structures whose files carry no (or incomplete) bonding
information, and by the Rebond command.
"""

import itertools
import math

import numpy as np

from miles.core.elements import BY_Z, COVALENT_RADII, REFERENCE_ORDERS, element, max_valence, standard_valence
from miles.core.geometry import Geometry, angle, geometry_for_coordination
from miles.core.molecule import Bond, Molecule
from miles.core.neighbors import pairs_within

BOND_TOLERANCE = 1.2
"""Atoms closer than this times the sum of covalent radii are bonded."""


def perceive_bonds(mol: Molecule, factor: float = BOND_TOLERANCE) -> None:
    """Connect atoms closer than ``factor`` times the sum of their covalent radii.

    Parameters
    ----------
    mol : Molecule
        Structure to bond; existing orders are kept for pairs that stay bonded.
    factor : float, optional
        Tolerance on the covalent-radius sum.
    """
    if len(mol) < 2:
        mol.set_bonds({})
        return
    radii = COVALENT_RADII[mol.numbers]
    is_h = mol.numbers == 1
    pairs = pairs_within(mol.coords, min(4.0, factor * 2 * float(radii.max())))
    old = dict(mol.bonds)
    bonds: dict[Bond, float] = {}
    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        d = np.linalg.norm(mol.coords[i] - mol.coords[j], axis=1)
        limit = np.where(is_h[i] & is_h[j], 0.85, factor * (radii[i] + radii[j]))
        keep = (d < limit) & (d > 0.4)
        for a, b in zip(i[keep].tolist(), j[keep].tolist(), strict=True):
            bonds[(a, b)] = old.get((a, b), 1.0)
    mol.set_bonds(bonds)


def assign_bond_orders(mol: Molecule) -> None:
    """Guess multiple bonds from bond lengths, then repair valence excesses.

    Each bond takes the order whose reference length is closest. "Aromatic"
    lengths outside rings (e.g. peptide C-N) become single, and when an atom
    exceeds its maximum valence the least certain bonds are reduced first,
    so a carbonyl keeps its C=O.

    Parameters
    ----------
    mol : Molecule
        Structure whose bond orders are rewritten in place.
    """
    stretch: dict[Bond, float] = {}
    for i, j in list(mol.bonds):
        order, ratio = _order_from_length(mol, i, j)
        mol.add_bond(i, j, order)
        if ratio is not None:
            stretch[(i, j)] = ratio
    ring = mol.ring_bonds()
    for key, order in list(mol.bonds.items()):
        if order == 1.5 and key not in ring:
            mol.add_bond(*key, 1.0)
    for i in range(len(mol)):
        _reduce_excess_valence(mol, i, stretch)


def _order_from_length(mol: Molecule, i: int, j: int) -> tuple[float, float | None]:
    """Bond order whose reference length is closest, and the length relative to it."""
    a, b = mol.atoms[i].symbol, mol.atoms[j].symbol
    ref = REFERENCE_ORDERS.get((a, b) if a <= b else (b, a))
    if ref is None:
        return 1.0, None
    d = float(np.linalg.norm(mol.coords[i] - mol.coords[j]))
    candidates = [(abs(d - length), order, length)
                  for length, order in zip(ref, (1.0, 1.5, 2.0, 3.0), strict=True) if length is not None]
    _, order, length = min(candidates)
    return order, d / length


def _reduce_excess_valence(mol: Molecule, i: int, stretch: dict[Bond, float]) -> None:
    """Lower multiple bonds on atom ``i`` until it is within its maximum valence.

    Aromatic bonds go first, then the most stretched multiple bond.
    """
    limit = max_valence(mol.atoms[i].symbol)

    def doubt(j: int) -> tuple[bool, float]:
        return mol.bond_order(i, j) == 1.5, stretch.get(mol.key(i, j), 1.0)

    for _ in range(10):
        if mol.order_sum(i) <= limit + 1e-6:
            return
        multiple = [j for j in mol.neighbors(i) if mol.bond_order(i, j) > 1.0]
        if not multiple:
            return
        j = max(multiple, key=doubt)
        order = mol.bond_order(i, j)
        mol.add_bond(i, j, 1.0 if order == 1.5 else order - 1.0)


def infer_atom_type(mol: Molecule, i: int, complete: bool = False) -> None:
    """Fill in the geometry, site count and valence of an atom from its bonds.

    Parameters
    ----------
    mol : Molecule
        Structure containing the atom.
    i : int
        Atom index.
    complete : bool, optional
        If true, assume the atom needs its standard valence, so that missing
        hydrogens can be added afterwards; otherwise the atom is taken as
        saturated as it is.
    """
    a = mol.atoms[i]
    a.inferred = True
    nbrs = mol.neighbors(i)
    deg = len(nbrs)
    used = mol.order_sum(i)
    if a.symbol == "H":
        a.geom, a.nsites, a.valence = Geometry.LINEAR, max(1, deg), max(1.0, used)
        return
    std = standard_valence(a.symbol)
    if complete and used < std:
        a.valence = float(std)
        a.nsites = deg + math.floor(std - used + 1e-6)
    else:
        a.valence = max(used, 1.0) if deg else float(std)
        a.nsites = deg if deg else std
    a.geom = infer_geometry(mol, i, a.nsites)


def infer_geometry(mol: Molecule, i: int, n: int) -> Geometry | None:
    """Geometry of atom ``i`` with ``n`` bonds, from bond orders and bond angles.

    Parameters
    ----------
    mol : Molecule
        Structure containing the atom.
    i : int
        Atom index.
    n : int
        Number of bonds the atom will have, including any still to be added.

    Returns
    -------
    Geometry or None
        Linear, trigonal, tetrahedral and so on; planar for atoms next to a
        pi system (amide and aniline nitrogens).
    """
    nbrs = mol.neighbors(i)
    orders = [mol.bond_order(i, j) for j in nbrs]
    conjugated = any(o >= 1.5 for o in orders)
    if mol.atoms[i].symbol == "N" and not conjugated and n in (2, 3):
        # amide, aniline and guanidine nitrogens are planar: next to a pi system
        conjugated = any(mol.bond_order(nb, x) >= 1.5 for nb in nbrs for x in mol.neighbors(nb) if x != i)
    center = mol.coords[i]
    pairs = list(itertools.combinations(nbrs, 2))
    if n <= 1:
        return Geometry.LINEAR
    if n == 2:
        straight = len(nbrs) == 2 and angle(mol.coords[nbrs[0]], center, mol.coords[nbrs[1]]) > 165
        if straight or any(o >= 3 for o in orders) or sum(o >= 2 for o in orders) == 2:
            return Geometry.LINEAR
        return Geometry.TRIGONAL if conjugated else Geometry.TETRAHEDRAL
    if n == 3:
        planar = len(nbrs) == 3 and sum(angle(mol.coords[p], center, mol.coords[q]) for p, q in pairs) > 350
        return Geometry.TRIGONAL if conjugated or planar else Geometry.TETRAHEDRAL
    if n == 4:
        if len(nbrs) == 4 and element(mol.atoms[i].symbol).is_metal:
            trans_pairs = sum(angle(mol.coords[p], center, mol.coords[q]) > 150 for p, q in pairs)
            if trans_pairs >= 2:
                return Geometry.SQUARE_PLANAR
        return Geometry.TETRAHEDRAL
    return geometry_for_coordination(n)


_STANDARD_VALENCE = np.array([standard_valence(e.symbol) for e in BY_Z])
_IS_METAL = np.array([e.is_metal for e in BY_Z])


def infer_atom_types(mol: Molecule, only_missing: bool = True) -> None:
    """Infer the bonding pattern of every atom, as :func:`infer_atom_type` does for one.

    Works on whole columns and the bond arrays, so typing a protein of
    100,000 atoms takes a fraction of a second.

    Parameters
    ----------
    mol : Molecule
        Structure to type.
    only_missing : bool, optional
        Skip atoms that already have a geometry, site count and valence.
    """
    n = len(mol)
    if not n:
        return
    if only_missing:
        todo = np.flatnonzero(np.fromiter((g is None for g in mol.column("geom").tolist()), bool, n)
                              | (mol.column("nsites") < 0) | np.isnan(mol.column("valence")))
    else:
        todo = np.arange(n)
    if not len(todo):
        return
    numbers = mol.numbers
    pairs, orders = mol.bond_pairs, mol.bond_orders
    a, b = pairs[:, 0], pairs[:, 1]

    def per_atom(weights: np.ndarray | None = None, bonds: np.ndarray | None = None) -> np.ndarray:
        """Sum over the bonds of each atom (counts, or weights), optionally over some bonds only."""
        ia, ib = (a, b) if bonds is None else (a[bonds], b[bonds])
        w = None if weights is None else (weights if bonds is None else weights[bonds])
        return np.bincount(ia, w, n) + np.bincount(ib, w, n)

    degree = per_atom().astype(np.int64)
    used = per_atom(orders)
    hydrogen = numbers == 1
    std = _STANDARD_VALENCE[numbers]
    valence = np.where(hydrogen, np.maximum(1.0, used), np.where(degree > 0, np.maximum(used, 1.0), std))
    nsites = np.where(hydrogen, np.maximum(1, degree), np.where(degree > 0, degree, std)).astype(np.int64)

    # what decides the geometry: bond orders, conjugation and bond angles
    pi = per_atom(bonds=orders >= 1.5) > 0
    neighbour_pi = (np.bincount(a[pi[b]], minlength=n) + np.bincount(b[pi[a]], minlength=n)) > 0
    nitrogen = numbers == 7
    # amide, aniline and guanidine nitrogens are planar: next to a pi system
    conjugated = pi | (nitrogen & np.isin(nsites, (2, 3)) & neighbour_pi)
    triple = per_atom(bonds=orders >= 3) > 0
    doubles = per_atom(bonds=orders >= 2)
    start, nbr = mol.adjacency()
    coords = mol.coords
    straight, planar, trans = (np.zeros(n, dtype=bool) for _ in range(3))
    for count, k_pairs in ((2, ((0, 1),)), (3, ((0, 1), (0, 2), (1, 2))),
                           (4, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)))):
        atoms = np.flatnonzero((degree == count) & (nsites == count))
        if not len(atoms):
            continue
        ends = nbr[start[atoms][:, None] + np.arange(count)]
        arms = coords[ends] - coords[atoms][:, None, :]
        lengths = np.linalg.norm(arms, axis=2, keepdims=True)
        arms = np.divide(arms, lengths, out=np.zeros_like(arms), where=lengths > 0)
        cosines = np.stack([np.einsum("ij,ij->i", arms[:, p], arms[:, q]) for p, q in k_pairs], axis=1)
        degrees = np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))
        if count == 2:
            straight[atoms] = degrees[:, 0] > 165
        elif count == 3:
            planar[atoms] = degrees.sum(axis=1) > 350
        else:
            trans[atoms] = (degrees > 150).sum(axis=1) >= 2

    # 0 linear, 1 trigonal, 2 tetrahedral, 3 square planar; -1 from the site count alone
    code = np.full(n, -1)
    code[nsites <= 1] = 0
    two = nsites == 2
    code[two] = np.where(straight | triple | (doubles == 2), 0, np.where(conjugated, 1, 2))[two]
    three = nsites == 3
    code[three] = np.where(conjugated | planar, 1, 2)[three]
    four = nsites == 4
    code[four] = np.where(_IS_METAL[numbers] & (degree == 4) & trans, 3, 2)[four]
    code[hydrogen] = 0
    shapes = np.array([Geometry.LINEAR, Geometry.TRIGONAL, Geometry.TETRAHEDRAL, Geometry.SQUARE_PLANAR,
                       None], dtype=object)
    geometry = shapes[code]   # code -1 picks the last entry, None
    for k in np.flatnonzero(code < 0).tolist():
        geometry[k] = geometry_for_coordination(int(nsites[k]))

    mol.set_column("inferred", True, todo)
    mol.set_column("geom", geometry[todo], todo)
    mol.set_column("nsites", nsites[todo], todo)
    mol.set_column("valence", valence[todo], todo)


def rebond(mol: Molecule) -> None:
    """Recompute bonds, bond orders and atom types from the coordinates alone."""
    perceive_bonds(mol)
    assign_bond_orders(mol)
    infer_atom_types(mol, only_missing=False)
