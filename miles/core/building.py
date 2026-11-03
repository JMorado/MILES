"""Building operations on a molecule.

Placing and attaching fragments, adding hydrogens and valence, changing
elements, deleting atoms and closing rings. Every function edits the
:class:`~miles.core.molecule.Molecule` it is given. Attaching a
fragment re-uses the index of the hydrogen it replaces, so the indices of
existing atoms stay valid.
"""

import itertools
import math
from collections.abc import Callable, Iterable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from miles.core.elements import bond_length, element
from miles.core.geometry import (
    TEMPLATES,
    Geometry,
    free_sites,
    geometry_for_coordination,
    pairwise_distances,
    rotation_between,
    rotation_matrix,
    signed_angle,
    unit_vector,
)
from miles.core.molecule import Atom, AtomLike, Molecule, hydrogen_atom
from miles.core.neighbors import PointGrid
from miles.core.perception import infer_atom_type

SCAN_ANGLES = np.radians(np.arange(0.0, 360.0, 10.0))
CLASH_REACH = 8.0   # angstrom; atoms further away barely change the d**-6 clash score
"""Rotations tried about a new bond when looking for the least crowded position."""


@dataclass(frozen=True, slots=True)
class AtomType:
    """How an element fragment bonds.

    Attributes
    ----------
    label : str
        Name shown to the user, e.g. ``"Tetrahedral"``.
    geom : str or None
        A :class:`~miles.core.geometry.Geometry` value.
    nsites : int
        Number of bonds.
    valence : float
        Bond-order sum.
    bare : bool
        Place the atom without hydrogens.
    """

    label: str
    geom: str | None
    nsites: int
    valence: float
    bare: bool = False


# -- helpers ------------------------------------------------------------------------

def extra_valence(mol: Molecule, i: int, exclude: int | None = None) -> float:
    """Bond order still available at atom ``i`` beyond one single bond.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i : int
        Atom index.
    exclude : int, optional
        Neighbour whose bond is about to be replaced (not counted).

    Returns
    -------
    float
        Spare valence; 0 for atoms without a known valence.
    """
    a = mol.atoms[i]
    if a.valence is None:
        return 0.0
    return a.valence - mol.order_sum(i, () if exclude is None else (exclude,)) - 1.0


def clash_score(points: np.ndarray, coords: np.ndarray) -> float:
    """Soft repulsion between new points and existing atoms (sum of ``d**-6``)."""
    if len(points) == 0 or len(coords) == 0:
        return 0.0
    d = np.maximum(pairwise_distances(points, coords), 0.5)
    return float((d ** -6).sum())


def _clash_scores(options: np.ndarray, coords: np.ndarray) -> np.ndarray:
    """:func:`clash_score` of many sets of points at once, shape ``(k, count, 3)`` to ``(k,)``."""
    if options.size == 0 or len(coords) == 0:
        return np.zeros(len(options))
    diff = options[:, :, None, :] - coords[None, None, :, :]
    d = np.maximum(np.sqrt(np.einsum("kpmi,kpmi->kpm", diff, diff)), 0.5)
    return (d ** -6).sum(axis=(1, 2))


def _rotations(axis: np.ndarray, angles: Iterable[float]) -> np.ndarray:
    """Rotation matrices about one axis for many angles, as :func:`rotation_matrix` makes them."""
    x, y, z = unit_vector(axis)
    a = np.asarray(list(angles), dtype=float)
    c, s = np.cos(a), np.sin(a)
    t = 1.0 - c
    return np.stack([np.stack([c + x * x * t, x * y * t - z * s, x * z * t + y * s], axis=-1),
                     np.stack([y * x * t + z * s, c + y * y * t, y * z * t - x * s], axis=-1),
                     np.stack([z * x * t - y * s, z * y * t + x * s, c + z * z * t], axis=-1)], axis=1)


def _least[T](options: Iterable[T], score: Callable[[T], float]) -> T | None:
    """First option with the lowest score (later ones must win by more than 1e-12)."""
    best, best_score = None, math.inf
    for option in options:
        s = score(option)
        if s < best_score - 1e-12:
            best, best_score = option, s
    return best


def _others(mol: Molecule, skip: Iterable[int]) -> np.ndarray:
    keep = np.ones(len(mol), dtype=bool)
    keep[list(skip)] = False
    return mol.coords[keep]


def _is_trigonal(atom: AtomLike) -> bool:
    return atom.geom == Geometry.TRIGONAL


def _plane_angles(axis: np.ndarray, ref_a: np.ndarray, ref_b: np.ndarray) -> list[float] | None:
    """Rotations about ``axis`` that bring the plane (axis, ref_b) onto (axis, ref_a)."""
    na, nb = unit_vector(np.cross(axis, ref_a)), unit_vector(np.cross(axis, ref_b))
    if np.linalg.norm(na) < 1e-6 or np.linalg.norm(nb) < 1e-6:
        return None
    phi = signed_angle(nb, na, axis)
    return [phi, phi + math.pi]


# -- fragments ------------------------------------------------------------------------

def place_fragment(mol: Molecule, frag: Molecule, position: ArrayLike, rotation: np.ndarray | None = None,
                   origin: ArrayLike | None = None) -> list[int]:
    """Add a free-standing copy of a fragment.

    Parameters
    ----------
    mol : Molecule
        Structure to add to.
    frag : Molecule
        Fragment to copy.
    position : array_like
        Where the fragment's ``origin`` ends up.
    rotation : numpy.ndarray, optional
        Rotation applied to the fragment about its origin.
    origin : array_like, optional
        Point of the fragment placed at ``position`` (default: its centroid).

    Returns
    -------
    list of int
        Indices of the new atoms.
    """
    if len(frag) == 0:
        return []
    R = np.eye(3) if rotation is None else rotation
    o = frag.centroid() if origin is None else np.asarray(origin, dtype=float)
    piece = frag.copy()
    piece.coords = (piece.coords - o) @ R.T + np.asarray(position, dtype=float)
    offset = mol.extend(piece)
    return list(range(offset, offset + len(piece)))


def attach_fragment(mol: Molecule, h: int, frag: Molecule, attach: int | None) -> list[int]:
    """Replace atom ``h`` (normally a hydrogen) by a fragment.

    The fragment bonds through the neighbour of its ``attach`` atom, which is
    discarded. The bond order follows from the spare valence on both sides,
    the bond length from the elements and order, and the rotation about the
    new bond is the least crowded one (conjugated trigonal atoms are kept
    coplanar).

    Parameters
    ----------
    mol : Molecule
        Structure to build on.
    h : int
        Atom to replace.
    frag : Molecule
        Fragment to attach.
    attach : int or None
        Fragment atom standing in for the bond (None for single-atom fragments).

    Returns
    -------
    list of int
        Indices of the fragment atoms in ``mol``; the anchor takes index ``h``.
    """
    nbrs = mol.neighbors(h)
    if not nbrs:
        origin = frag.coords[0] if attach is None else frag.coords[frag.neighbors(attach)[0]]
        pos = mol.coords[h].copy()
        mol.remove_atoms([h])
        return place_fragment(mol, frag, pos, origin=origin)

    x = nbrs[0]
    if attach is None:
        anchor = 0
        spare = (frag.atoms[0].valence or 1) - frag.order_sum(0) - 1
    else:
        anchor = frag.neighbors(attach)[0]
        spare = extra_valence(frag, anchor, exclude=attach)
    order = 1.0 + math.floor(max(0.0, min(extra_valence(mol, x, exclude=h), spare, 2.0)) + 1e-6)

    xpos = mol.coords[x]
    u = unit_vector(mol.coords[h] - xpos)
    apos = xpos + u * bond_length(mol.atoms[x].symbol, frag.atoms[anchor].symbol, order)
    P = frag.coords - frag.coords[anchor]
    if attach is not None:
        P = P @ rotation_between(frag.coords[attach] - frag.coords[anchor], -u).T
    P = P + apos

    keep = [i for i in range(len(frag)) if i != attach]
    movers = [i for i in keep if i != anchor]
    if movers:
        candidates: Iterable[float] | None = None
        if _is_trigonal(mol.atoms[x]) and _is_trigonal(frag.atoms[anchor]):
            ref_x = [n for n in mol.neighbors(x) if n != h]
            ref_a = [n for n in frag.neighbors(anchor) if n != attach]
            if ref_x and ref_a:
                candidates = _plane_angles(u, mol.coords[ref_x[0]] - xpos, P[ref_a[0]] - apos)
        env = _others(mol, [h])

        def crowding(phi: float) -> float:
            return clash_score((P[movers] - apos) @ rotation_matrix(u, phi).T + apos, env)

        best = _least(SCAN_ANGLES if candidates is None else candidates, crowding)
        P[movers] = (P[movers] - apos) @ rotation_matrix(u, float(best or 0.0)).T + apos

    mapping = {anchor: h}
    mol.atoms[h] = frag.atoms[anchor]
    mol.set_positions([h], P[anchor])
    mol.add_bond(x, h, order)
    for i in movers:
        mapping[i] = mol.add_atom(frag.atoms[i], P[i])
    for (a, b), o in frag.bonds.items():
        if attach not in (a, b):
            mol.add_bond(mapping[a], mapping[b], o)
    return [mapping[i] for i in keep]


# -- hydrogens and valence ---------------------------------------------------------------

def site_directions(mol: Molecule, i: int, neighbors: list[int], geom: str | None, count: int,
                    ignore: Iterable[int] = (), around: PointGrid | None = None) -> list[np.ndarray]:
    """Directions for new bonds on an atom, avoiding nearby atoms.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i : int
        Atom receiving the bonds.
    neighbors : list of int
        Bonds the atom keeps.
    geom : str or None
        Target geometry (chosen from the coordination number when None).
    count : int
        Number of directions wanted.
    ignore : iterable of int, optional
        Atoms not counted as obstacles (e.g. hydrogens being replaced).
    around : PointGrid, optional
        The atoms' positions binned for quick searches, with the reach of
        :data:`CLASH_REACH`; makes many calls on a large structure fast.
        Its first ``len(mol)`` points must be the atoms, in order.

    Returns
    -------
    list of numpy.ndarray
        Unit vectors, at most ``count`` of them.
    """
    if count <= 0:
        return []
    center = mol.coords[i]
    dirs = [d for d in (unit_vector(mol.coords[n] - center) for n in neighbors) if np.linalg.norm(d) > 0]
    if geom is None or geom not in TEMPLATES or len(dirs) >= len(TEMPLATES[geom]):
        geom = geometry_for_coordination(len(dirs) + count) or Geometry.OCTAHEDRAL
    free = free_sites(geom, dirs)
    if not free:
        return []
    count = min(count, len(free))
    skip = {i, *ignore}
    if around is None:
        env = _others(mol, skip)
        env = env[((env - center) ** 2).sum(axis=1) < CLASH_REACH ** 2]
    else:
        near = around.near(center)
        near = near[~np.isin(near, list(skip))]
        env = around.position(near)
        env = env[((env - center) ** 2).sum(axis=1) < CLASH_REACH ** 2]
    h_len = bond_length(mol.atoms[i].symbol, "H")

    rotations = np.eye(3)[None]
    if len(dirs) == 1:
        n = neighbors[0]
        planar = None
        if geom == Geometry.TRIGONAL and _is_trigonal(mol.atoms[n]):
            refs = [m for m in mol.neighbors(n) if m != i]
            if refs:
                planar = _plane_angles(dirs[0], mol.coords[refs[0]] - mol.coords[n], free[0])
        rotations = _rotations(dirs[0], SCAN_ANGLES if planar is None else planar)

    # every rotation with every choice of sites, scored at once; the first best wins, as with _least
    combos = list(itertools.combinations(range(len(free)), count))
    turned = np.einsum("rij,fj->rfi", rotations, np.array(free))
    options = turned[:, np.array(combos)].reshape(-1, count, 3)
    scores = _clash_scores(center + options * h_len, env)
    best = _least(range(len(options)), lambda k: float(scores[k]))
    return [] if best is None else list(options[best])


def hydrogen_for(mol: Molecule, parent: int) -> Atom:
    """Make a hydrogen in the same residue as its parent and drawn the same way.

    Parameters
    ----------
    mol : Molecule
        Structure.
    parent : int
        Atom the hydrogen will bond to.

    Returns
    -------
    Atom
        Hydrogen record in the parent's residue (so a document shows it like its parent).
    """
    p = mol.atoms[parent]
    h = hydrogen_atom()
    h.resname, h.resseq, h.icode, h.chain, h.hetero = p.resname, p.resseq, p.icode, p.chain, p.hetero
    if p.resname:
        h.name = "H" + p.name[1:] if len(p.name) > 1 else "H"
    return h


def missing_hydrogens(mol: Molecule, i: int) -> int:
    """Return the number of hydrogens atom ``i`` lacks for its intended bonding."""
    a = mol.atoms[i]
    if a.symbol == "H" or a.nsites is None or a.valence is None:
        return 0
    free_bonds = a.nsites - mol.degree(i)
    free_valence = math.floor(a.valence - mol.order_sum(i) + 1e-6)
    return max(0, min(free_bonds, free_valence))


def _add_hydrogen(mol: Molecule, i: int, direction: np.ndarray) -> int:
    h = mol.add_atom(hydrogen_for(mol, i), mol.coords[i] + direction * bond_length(mol.atoms[i].symbol, "H"))
    mol.add_bond(i, h, 1.0)
    return h


def add_hydrogens(mol: Molecule, indices: Iterable[int] | None = None) -> int:
    """Fill open valences with hydrogens.

    Monatomic ions (``ion`` in the selection language: metal and halide
    atoms on their own) get none, and neither do metal atoms whose bonding
    was read from a file rather than set in the builder: their usual
    valence says nothing about how many bonds they carry.

    Parameters
    ----------
    mol : Molecule
        Structure.
    indices : iterable of int, optional
        Atoms to complete (default: all). Atoms typed from files fall back
        to their standard valence.

    Returns
    -------
    int
        Number of hydrogens added.
    """
    from miles.core.selection import find  # the selection language defines "ion" once

    targets = list(range(len(mol)) if indices is None else indices)
    if not targets:
        return 0
    ions = set(find(mol, "ion").tolist())
    around = PointGrid(mol.coords, CLASH_REACH)
    added = 0
    for i in targets:
        atom = mol.atoms[i]
        if atom.symbol == "H" or i in ions:
            continue
        guessed = atom.nsites is None or atom.inferred
        if guessed and element(atom.symbol).is_metal:
            continue
        if guessed:
            infer_atom_type(mol, i, complete=True)
        count = missing_hydrogens(mol, i)
        for d in site_directions(mol, i, mol.neighbors(i), atom.geom, count, around=around):
            h = _add_hydrogen(mol, i, d)
            around.add(mol.coords[h])
            added += 1
    return added


def add_valence(mol: Molecule, i: int) -> int | None:
    """Add one hydrogen to atom ``i``, raising its coordination when it is full.

    Returns
    -------
    int or None
        Index of the new hydrogen, or None if no position is left.
    """
    a = mol.atoms[i]
    if a.nsites is None:
        infer_atom_type(mol, i)
    if missing_hydrogens(mol, i) <= 0:
        geom = geometry_for_coordination(mol.degree(i) + 1)
        if geom is None:
            return None
        a.nsites, a.geom = mol.degree(i) + 1, geom
        a.valence = max(a.valence or 0.0, mol.order_sum(i) + 1.0)
    dirs = site_directions(mol, i, mol.neighbors(i), a.geom, 1)
    return _add_hydrogen(mol, i, dirs[0]) if dirs else None


def mutate_atom(mol: Molecule, i: int, symbol: str, kind: AtomType) -> None:
    """Turn atom ``i`` into another element or bonding type and rebuild its hydrogens.

    Bonds to heavy neighbours are kept but capped at what the new type allows.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i : int
        Atom to change.
    symbol : str
        New element.
    kind : AtomType
        New bonding type.
    """
    old_h = mol.hydrogens_on(i)
    heavy = [n for n in mol.neighbors(i) if n not in old_h]
    a = mol.atoms[i]
    a.symbol, a.geom, a.nsites, a.valence, a.inferred = symbol, kind.geom, kind.nsites, kind.valence, False

    max_order = 1.0 + max(0.0, kind.valence - kind.nsites)
    for n in heavy:
        if mol.bond_order(i, n) > max_order:
            mol.add_bond(i, n, max_order)
    while heavy and sum(mol.bond_order(i, n) for n in heavy) > kind.valence + 1e-6:
        n = max(heavy, key=lambda m: mol.bond_order(i, m))
        if mol.bond_order(i, n) <= 1.0:
            break
        mol.add_bond(i, n, max(1.0, math.ceil(mol.bond_order(i, n)) - 1.0))

    count = 0
    if not kind.bare:
        used = sum(mol.bond_order(i, n) for n in heavy)
        count = max(0, min(kind.nsites - len(heavy), math.floor(kind.valence - used + 1e-6)))
    dirs = site_directions(mol, i, heavy, kind.geom, count, ignore=old_h)
    h_len = bond_length(symbol, "H")
    for k, d in enumerate(dirs):
        if k < len(old_h):
            mol.set_positions([old_h[k]], mol.coords[i] + d * h_len)
        else:
            _add_hydrogen(mol, i, d)
    if len(old_h) > len(dirs):
        mol.remove_atoms(old_h[len(dirs):])


# -- deleting and joining ------------------------------------------------------------------

def delete_atoms(mol: Molecule, indices: Iterable[int]) -> list[int]:
    """Delete atoms together with the terminal hydrogens of deleted heavy atoms.

    Returns
    -------
    list of int
        Old-to-new index map (-1 for deleted atoms).
    """
    doomed: set[int] = set()
    for i in indices:
        doomed.add(i)
        if mol.atoms[i].symbol != "H":
            doomed.update(mol.hydrogens_on(i))
    return mol.remove_atoms(doomed)


def delete_atom(mol: Molecule, i: int) -> list[int]:
    """Delete one atom and its terminal hydrogens (see :func:`delete_atoms`)."""
    return delete_atoms(mol, [i])


def join_atoms(mol: Molecule, i: int, j: int) -> tuple[int, int] | None:
    """Bond two atoms, for example to close a ring.

    A hydrogen given as an end stands for its parent atom and is removed.
    A saturated end gives up the hydrogen pointing most towards the other end.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i, j : int
        Atoms (or hydrogens on them) to join.

    Returns
    -------
    tuple of int or None
        New indices of the bonded pair, or None if nothing was done.
    """
    ends, doomed = [], set()
    for a in (i, j):
        if mol.atoms[a].symbol == "H" and mol.degree(a) == 1:
            doomed.add(a)
            a = mol.neighbors(a)[0]
        ends.append(a)
    a, b = ends
    if a == b or mol.bond_order(a, b) > 0:
        return None
    for x, y in ((a, b), (b, a)):
        atom = mol.atoms[x]
        hs = [h for h in mol.hydrogens_on(x) if h not in doomed]
        removed = sum(1 for h in doomed if mol.bond_order(x, h) > 0)
        saturated = atom.valence is not None and mol.order_sum(x) - removed >= atom.valence - 1e-6
        if saturated and hs and removed == 0:
            doomed.add(_most_aligned(mol, x, hs, mol.coords[y] - mol.coords[x]))
    mol.add_bond(a, b, 1.0)
    mapping = mol.remove_atoms(doomed)
    return mapping[a], mapping[b]


def _most_aligned(mol: Molecule, i: int, candidates: list[int], direction: np.ndarray) -> int:
    """Neighbour of atom ``i`` whose bond points most nearly along ``direction``."""
    toward = unit_vector(direction)
    return max(candidates, key=lambda h: float(np.dot(unit_vector(mol.coords[h] - mol.coords[i]), toward)))


def best_hydrogen_towards(mol: Molecule, i: int, direction: ArrayLike,
                          project: Callable[[np.ndarray], np.ndarray]) -> int | None:
    """Hydrogen on atom ``i`` whose bond points most nearly along a direction.

    Parameters
    ----------
    mol : Molecule
        Structure.
    i : int
        Atom (if it is itself a hydrogen, it is returned).
    direction : array_like
        Direction to match, in the space ``project`` maps to.
    project : callable
        Maps a 3D position to the comparison space (e.g. screen coordinates).

    Returns
    -------
    int or None
        Index of the hydrogen, or None if the atom has none.
    """
    if mol.atoms[i].symbol == "H":
        return i
    hs = mol.hydrogens_on(i)
    if not hs:
        return None
    d = unit_vector(np.asarray(direction, dtype=float))
    origin = project(mol.coords[i])
    return max(hs, key=lambda h: float(np.dot(unit_vector(project(mol.coords[h]) - origin), d)))
