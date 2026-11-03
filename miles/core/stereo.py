"""Stereochemistry: CIP descriptors (R/S and E/Z) from 3D coordinates, and changing them.

Priorities follow rule 1a of the Cahn-Ingold-Prelog sequence rules: the
substituent with the higher atomic number ranks higher, and ties are broken
by exploring the hierarchical digraph sphere by sphere. Multiple bonds add
duplicate atoms, and a ring closure ends in a duplicate of the atom it
returns to. The label of a centre then comes from the 3D geometry, so it is
correct for whatever structure is on screen, built or read from a file.

Approximations
--------------
* Only rule 1a is used. Isotopes (rule 2) and the stereo rules (3 to 5) are
  not, so centres whose branches differ only in their own stereochemistry,
  such as pseudoasymmetric centres, are not labelled.
* Within a sphere, nodes are explored in the order of their parent's rank
  and then of their own atomic number; equal-looking nodes are not ranked
  further by what lies beyond them.
* An atom with aromatic bonds (order 1.5) gets one duplicate, of its
  aromatic neighbour with the highest atomic number, which matches a Kekule
  structure for the usual aromatic rings.
* Exploration stops after :data:`MAX_SPHERES` spheres; branches that are
  still equal then are treated as equal.
"""

import functools
import math
from collections.abc import Iterator
from dataclasses import dataclass
from enum import StrEnum
from itertools import pairwise

import numpy as np

from miles.core.elements import element
from miles.core.geometry import dihedral, rotation_matrix, unit_vector
from miles.core.molecule import Molecule

MAX_SPHERES = 12
"""Deepest sphere of the digraph explored when comparing two branches."""

MIN_STEREO_RING = 8
"""Double bonds in smaller rings cannot be E/Z isomers and are not labelled."""

LONE_PAIR_ELEMENTS = {7, 15}
"""Double-bond ends with one substituent and a lone pair (imines, oximes, phosphaalkenes)."""

Signature = tuple[tuple[int, ...], ...]


class Chirality(StrEnum):
    """CIP descriptor of a tetrahedral stereocentre."""

    R = "R"
    S = "S"


class BondStereo(StrEnum):
    """CIP descriptor of a double bond."""

    E = "E"
    Z = "Z"


@dataclass(frozen=True, slots=True)
class _Node:
    """A node of the hierarchical digraph: an atom reached along a path."""

    atom: int
    parent: int
    path: frozenset[int]
    duplicate: bool


class _Graph:
    """Atomic numbers and neighbours with bond orders, prepared once per molecule."""

    def __init__(self, mol: Molecule) -> None:
        n = len(mol)
        self.z = mol.numbers.tolist()
        self.adj: list[list[tuple[int, float]]] = [[] for _ in range(n)]
        for (i, j), order in mol.bonds.items():
            self.adj[i].append((j, order))
            self.adj[j].append((i, order))
        self.duplicates = [self._duplicates(i) for i in range(n)]

    def _duplicates(self, i: int) -> list[int]:
        """Atoms duplicated at atom ``i`` because of its multiple bonds."""
        result = []
        aromatic = []
        for j, order in self.adj[i]:
            if order == 1.5:
                aromatic.append(j)
            elif order >= 2.0:
                result.extend([j] * (round(order) - 1))
        if aromatic:
            result.append(max(aromatic, key=lambda j: self.z[j]))
        return result

    def children(self, node: _Node) -> list[_Node]:
        """Nodes one sphere further out, highest atomic number first."""
        if node.duplicate:
            return []
        atom, path = node.atom, node.path
        kids = []
        for j, _ in self.adj[atom]:
            if j == node.parent:
                continue
            if j in path:
                kids.append(_Node(j, atom, path, True))   # ring closure
            else:
                kids.append(_Node(j, atom, path | {j}, False))
        kids.extend(_Node(j, atom, path, True) for j in self.duplicates[atom])
        kids.sort(key=lambda k: -self.z[k.atom])   # stable: real atoms before duplicates
        return kids

    def spheres(self, root: int, start: int) -> Iterator[Signature]:
        """Atomic numbers of each sphere of the branch from ``root`` through ``start``.

        The first signature is the start atom itself; each later one lists,
        for every node of the previous sphere in hierarchical order, the
        atomic numbers of its children in decreasing order.
        """
        yield ((self.z[start],),)
        frontier = [_Node(start, root, frozenset((root, start)), False)]
        for _ in range(MAX_SPHERES):
            signature = []
            following = []
            for node in frontier:
                kids = self.children(node)
                signature.append(tuple(self.z[k.atom] for k in kids))
                following.extend(kids)
            if not following:
                return
            yield tuple(signature)
            frontier = following

    def rank(self, root: int, starts: list[int]) -> list[int] | None:
        """Substituents of ``root`` from highest to lowest priority, or None if two are equal."""
        branches = {s: _Branch(self.spheres(root, s)) for s in starts}

        def compare(a: int, b: int) -> int:
            return branches[a].compare(branches[b])

        ordered = sorted(starts, key=functools.cmp_to_key(compare), reverse=True)
        if any(compare(a, b) == 0 for a, b in pairwise(ordered)):
            return None
        return ordered


class _Branch:
    """Sphere signatures of one branch, computed only as far as comparisons need."""

    def __init__(self, spheres: Iterator[Signature]) -> None:
        self._spheres = spheres
        self._known: list[Signature] = []

    def sphere(self, k: int) -> Signature:
        """Signature of sphere ``k`` (empty once the branch has ended)."""
        while len(self._known) <= k:
            following = next(self._spheres, None)
            if following is None:
                return ()
            self._known.append(following)
        return self._known[k]

    def compare(self, other: "_Branch") -> int:
        """Return 1 if this branch ranks higher, -1 if lower, 0 if equal."""
        for k in range(MAX_SPHERES + 1):
            a, b = self.sphere(k), other.sphere(k)
            if a != b:
                return 1 if a > b else -1
            if not a:
                return 0
        return 0


# -- labelling ------------------------------------------------------------------------

def stereocenters(mol: Molecule) -> dict[int, Chirality]:
    """CIP label of every tetrahedral stereocentre, from the 3D geometry.

    A stereocentre is a non-metal atom with four single bonds whose four
    branches all have different priorities.

    Parameters
    ----------
    mol : Molecule
        Structure with explicit hydrogens.

    Returns
    -------
    dict of int to Chirality
        Label of each stereocentre by atom index.
    """
    graph = _Graph(mol)
    labels: dict[int, Chirality] = {}
    for i, neighbors in enumerate(graph.adj):
        if len(neighbors) != 4 or any(order != 1.0 for _, order in neighbors):
            continue
        if sum(graph.z[j] == 1 for j, _ in neighbors) > 1 or element(mol.atoms[i].symbol).is_metal:
            continue
        ranked = graph.rank(i, [j for j, _ in neighbors])
        if ranked is not None:
            labels[i] = _chirality(mol.coords, i, ranked)
    return labels


def _chirality(coords: np.ndarray, center: int, ranked: list[int]) -> Chirality:
    """R or S from the signed volume of the four substituents, highest priority first."""
    a, b, c, d = (coords[k] - coords[center] for k in ranked)
    volume = float(np.dot(a - d, np.cross(b - d, c - d)))
    return Chirality.R if volume < 0 else Chirality.S


def double_bond_stereo(mol: Molecule) -> dict[tuple[int, int], BondStereo]:
    """E/Z label of every stereogenic double bond, from the 3D geometry.

    A double bond is labelled when it is not in a ring of fewer than
    :data:`MIN_STEREO_RING` atoms and each end has two substituents of
    different priority (or one substituent and a lone pair, on N or P).

    Parameters
    ----------
    mol : Molecule
        Structure with explicit hydrogens.

    Returns
    -------
    dict of tuple of int to BondStereo
        Label of each stereogenic double bond, keyed by ``Molecule.key(i, j)``.
    """
    graph = _Graph(mol)
    labels: dict[tuple[int, int], BondStereo] = {}
    for (i, j), order in mol.bonds.items():
        if order != 2.0:
            continue
        top_i = _top_substituent(graph, i, j)
        top_j = _top_substituent(graph, j, i) if top_i is not None else None
        if top_i is None or top_j is None or _ring_size(graph, i, j) < MIN_STEREO_RING:
            continue
        angle = dihedral(mol.coords[top_i], mol.coords[i], mol.coords[j], mol.coords[top_j])
        labels[(i, j)] = BondStereo.Z if abs(angle) < 90.0 else BondStereo.E
    return labels


def _top_substituent(graph: _Graph, end: int, other: int) -> int | None:
    """Highest-priority substituent at one end of a double bond, or None if the end is not stereogenic."""
    rest = [(k, order) for k, order in graph.adj[end] if k != other]
    if not 1 <= len(rest) <= 2 or any(order >= 2.0 for _, order in rest):
        return None
    if len(rest) == 1:
        return rest[0][0] if graph.z[end] in LONE_PAIR_ELEMENTS else None
    ranked = graph.rank(end, [k for k, _ in rest])
    return None if ranked is None else ranked[0]


def _ring_size(graph: _Graph, i: int, j: int) -> int:
    """Atoms in the smallest ring through bond i-j (:data:`MIN_STEREO_RING` if none is smaller)."""
    seen = {i}
    frontier = [i]
    for depth in range(1, MIN_STEREO_RING - 1):
        following = []
        for a in frontier:
            for b, _ in graph.adj[a]:
                if a == i and b == j:
                    continue
                if b == j:
                    return depth + 1
                if b not in seen:
                    seen.add(b)
                    following.append(b)
        frontier = following
    return MIN_STEREO_RING


# -- changing stereochemistry ------------------------------------------------------------

def invert_center(mol: Molecule, i: int) -> bool:
    """Invert a stereocentre in place by swapping two of its substituents.

    The two smallest substituents that are not part of a ring through the
    centre are turned 180 degrees about the bisector of their bonds, which
    swaps their places and keeps their own shape and configuration.

    Parameters
    ----------
    mol : Molecule
        Structure, changed in place.
    i : int
        Centre with four bonds.

    Returns
    -------
    bool
        False if the centre has no two substituents that can be swapped
        (for example when all four are in rings).
    """
    neighbors = mol.neighbors(i)
    if len(neighbors) != 4:
        return False
    others = set(neighbors)
    branches = []
    for n in neighbors:
        group = mol.component(n, blocked={i})
        if not group & (others - {n}):
            branches.append((len(group), n, group))
    if len(branches) < 2:
        return False
    branches.sort(key=lambda branch: branch[0])
    (_, a, group_a), (_, b, group_b) = branches[:2]
    center = mol.coords[i].copy()
    axis = unit_vector(unit_vector(mol.coords[a] - center) + unit_vector(mol.coords[b] - center))
    if not np.any(axis):
        return False
    idx = sorted(group_a | group_b)
    mol.set_positions(idx, (mol.coords[idx] - center) @ rotation_matrix(axis, math.pi).T + center)
    return True


def flip_double_bond(mol: Molecule, i: int, j: int) -> bool:
    """Swap E and Z in place by turning one side of a double bond half a turn.

    The smaller side is rotated 180 degrees about the bond axis.

    Parameters
    ----------
    mol : Molecule
        Structure, changed in place.
    i, j : int
        Atoms of the double bond.

    Returns
    -------
    bool
        False if the bond is not a double bond or lies in a ring.
    """
    if mol.bond_order(i, j) != 2.0:
        return False
    side_j = mol.component(j, blocked={i})
    if any(n in side_j for n in mol.neighbors(i) if n != j):   # the two sides meet in a ring
        return False
    side_i = mol.component(i, blocked={j})
    group, pivot, end = (side_j, i, j) if len(side_j) <= len(side_i) else (side_i, j, i)
    origin = mol.coords[end].copy()
    turn = rotation_matrix(mol.coords[end] - mol.coords[pivot], math.pi)
    idx = sorted(group)
    mol.set_positions(idx, (mol.coords[idx] - origin) @ turn.T + origin)
    return True
