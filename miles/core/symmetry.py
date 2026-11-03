"""Point group detection and symmetrisation.

The search works on candidate symmetry elements: rotation axes and mirror
normals built from the principal axes of inertia and from shells of
equivalent atoms. Each candidate operation is tested against the structure
with a distance tolerance, the elements found are classified into a
Schoenflies symbol, and the full group is generated from a few exact
generators.
"""

import itertools
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

from miles.core.elements import element
from miles.core.geometry import rotation_matrix, unit_vector

INF = "\u221e"
"""Infinity sign used in the names of linear point groups."""

MAX_AXIS_ORDER = 8

Axis = tuple[np.ndarray, int]
"""A rotation axis (unit vector) and its order."""


@dataclass
class PointGroup:
    """A point group and its operations.

    Attributes
    ----------
    name : str
        Schoenflies symbol, e.g. ``"C2v"``; linear groups are written with :data:`INF`.
    operations : list of numpy.ndarray
        The ``(3, 3)`` matrices of the group, about ``center``.
    center : numpy.ndarray
        Centre of mass of the structure.
    """

    name: str
    operations: list[np.ndarray] = field(default_factory=list)
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))

    @property
    def order(self) -> int:
        """Number of operations."""
        return len(self.operations)


class _Rotor(StrEnum):
    """Shape of the inertia tensor, which decides where to look for axes."""

    SPHERICAL = "spherical"    # three equal moments: cubic or icosahedral
    SYMMETRIC = "symmetric"    # two equal moments: one unique main axis
    ASYMMETRIC = "asymmetric"  # no equal moments: axes along principal axes


def _reflection(normal: np.ndarray) -> np.ndarray:
    n = unit_vector(normal)
    return np.eye(3) - 2.0 * np.outer(n, n)


def _improper(axis: np.ndarray, n: int) -> np.ndarray:
    return _reflection(axis) @ rotation_matrix(axis, 2 * math.pi / n)


def _rotation(axis: np.ndarray, order: int) -> np.ndarray:
    return rotation_matrix(axis, 2 * math.pi / order)


def _dedupe(vectors: Iterable[np.ndarray], cos_tol: float = 0.995) -> list[np.ndarray]:
    """Normalise directions and drop zero vectors and (anti)parallel repeats."""
    out: list[np.ndarray] = []
    for v in vectors:
        n = np.linalg.norm(v)
        if n < 1e-6:
            continue
        v = v / n
        if all(abs(float(v @ u)) <= cos_tol for u in out):
            out.append(v)
    return out


class _Analyzer:
    """A structure centred on its centre of mass, with fast symmetry tests."""

    def __init__(self, symbols: Sequence[str], coords: np.ndarray, tol: float) -> None:
        self.tol = tol
        self.masses = np.array([element(s).mass or 1.0 for s in symbols])
        self.center = (coords * self.masses[:, None]).sum(axis=0) / self.masses.sum()
        self.x = coords - self.center
        self.radius = np.linalg.norm(self.x, axis=1)
        by_element: dict[str, list[int]] = {}
        for i, s in enumerate(symbols):
            by_element.setdefault(s, []).append(i)
        self.group_idx = [np.array(v) for v in by_element.values()]
        # Probe atoms (rarest element, farthest from the centre) reject most
        # candidate operations before the full test.
        rare = min(self.group_idx, key=len)
        far = rare[np.argsort(-self.radius[rare])[:2]]
        self.probes = [(int(i), rare) for i in far]

    def is_symmetry(self, M: np.ndarray) -> bool:
        """Whether ``M`` maps every atom onto an atom of the same element."""
        for i, members in self.probes:
            if np.min(np.linalg.norm(self.x[members] - M @ self.x[i], axis=1)) > self.tol:
                return False
        y = self.x @ M.T
        for idx in self.group_idx:
            diff = y[idx][:, None, :] - self.x[idx][None, :, :]
            if np.any(np.sqrt((diff * diff).sum(axis=-1)).min(axis=1) > self.tol):
                return False
        return True

    def permutation(self, M: np.ndarray) -> np.ndarray:
        """Index of the atom each atom is sent to by ``M``."""
        y = self.x @ M.T
        perm = np.arange(len(self.x))
        for idx in self.group_idx:
            diff = y[idx][:, None, :] - self.x[idx][None, :, :]
            perm[idx] = idx[(diff * diff).sum(axis=-1).argmin(axis=1)]
        return perm

    def shells(self) -> list[list[int]]:
        """Return groups of same-element atoms at the same distance from the centre."""
        result = []
        for idx in self.group_idx:
            r = self.radius[idx]
            order = np.argsort(r)
            current = [int(idx[order[0]])]
            for a in order[1:]:
                if abs(r[a] - self.radius[current[-1]]) >= self.tol:
                    result.append(current)
                    current = []
                current.append(int(idx[a]))
            result.append(current)
        return [s for s in result if self.radius[s[0]] > self.tol]

    def candidate_axes(self, principal: np.ndarray, rotor: _Rotor,
                       unique: np.ndarray | None = None) -> list[np.ndarray]:
        """Directions that may be rotation axes."""
        candidates = [principal[:, k] for k in range(3)]
        if rotor == _Rotor.ASYMMETRIC:
            return _dedupe(candidates)
        for shell in self.shells():
            pts = self.x[shell]
            candidates.extend(pts)
            if len(shell) <= 40:
                for a, b in itertools.combinations(range(len(shell)), 2):
                    candidates.append(pts[a] + pts[b])
                    candidates.append(np.cross(pts[a], pts[b]))
            if rotor == _Rotor.SPHERICAL and len(shell) <= 14:
                for a, b, c in itertools.combinations(range(len(shell)), 3):
                    candidates.append(pts[a] + pts[b] + pts[c])
        if rotor == _Rotor.SYMMETRIC and unique is not None:
            candidates = [unique, *(v - (v @ unique) * unique for v in candidates)]
        return _dedupe(candidates)

    def candidate_normals(self, principal: np.ndarray, rotor: _Rotor, axes: list[Axis]) -> list[np.ndarray]:
        """Directions that may be mirror-plane normals."""
        candidates = [principal[:, k] for k in range(3)]
        if rotor != _Rotor.ASYMMETRIC:
            for shell in self.shells():
                pts = self.x[shell]
                if len(shell) <= 60:
                    for a, b in itertools.combinations(range(len(shell)), 2):
                        candidates.append(pts[a] - pts[b])
                for axis, _ in axes:
                    candidates.extend(np.cross(axis, p) for p in pts[:30])
            candidates.extend(axis for axis, _ in axes)
        return _dedupe(candidates)


def _parallel(a: np.ndarray, b: np.ndarray) -> bool:
    return abs(float(a @ b)) > 0.99


def _perpendicular(a: np.ndarray, b: np.ndarray) -> bool:
    return abs(float(a @ b)) < 0.05


def _classify(axes: list[Axis], mirrors: list[np.ndarray], inversion: bool,
              impropers: list[Axis]) -> tuple[str, np.ndarray | None, np.ndarray | None]:
    """Name the point group from the symmetry elements found.

    Returns
    -------
    name : str
        Schoenflies symbol.
    main : numpy.ndarray or None
        Main axis (None for groups without one, or with several).
    second : numpy.ndarray or None
        A perpendicular C2 axis (D groups) or mirror normal (Cnv groups).
    """
    if not axes:
        return ("Cs" if mirrors else "Ci" if inversion else "C1"), None, None
    if sum(order >= 3 for _, order in axes) >= 2:
        orders = {order for _, order in axes}
        if 5 in orders:
            name = "Ih" if inversion else "I"
        elif 4 in orders:
            name = "Oh" if inversion else "O"
        else:
            name = "Th" if inversion else "Td" if mirrors else "T"
        return name, None, None

    n = max(order for _, order in axes)
    tops = [a for a, order in axes if order == n]
    main = next((a for a in tops if any(_parallel(a, s) and m == 2 * n for s, m in impropers)), tops[0])
    perp = [a for a, order in axes if order >= 2 and _perpendicular(a, main)]
    sigma_h = any(_parallel(m, main) for m in mirrors)
    sigma_v = [m for m in mirrors if _perpendicular(m, main)]
    s2n = any(_parallel(a, main) and m == 2 * n for a, m in impropers)
    if perp:
        return f"D{n}{'h' if sigma_h else 'd' if sigma_v else ''}", main, perp[0]
    if sigma_h:
        return f"C{n}h", main, None
    if sigma_v:
        return f"C{n}v", main, sigma_v[0]
    if s2n:
        return f"S{2 * n}", main, None
    return f"C{n}", main, None


def _close_group(generators: list[np.ndarray], tol: float = 1e-3, limit: int = 240) -> list[np.ndarray]:
    """All products of the generators, up to ``limit`` operations."""
    ops = [np.eye(3)]

    def known(M: np.ndarray) -> bool:
        return any(np.abs(M - op).max() < tol for op in ops)

    for g in generators:
        if not known(g):
            ops.append(np.asarray(g, dtype=float))
    changed = True
    while changed and len(ops) < limit:
        changed = False
        for a in list(ops):
            for b in list(ops):
                c = a @ b
                if not known(c):
                    ops.append(c)
                    changed = True
                    if len(ops) >= limit:
                        break
    return ops


def _is_linear(x: np.ndarray, tol: float) -> bool:
    d = unit_vector(x[np.argmax(np.linalg.norm(x, axis=1))])
    if not np.any(d):
        return True
    return bool(np.linalg.norm(x - np.outer(x @ d, d), axis=1).max() < tol)


def _merge_axes(axes: list[Axis]) -> list[Axis]:
    """Drop repeated axes, keeping the highest order of each."""
    merged: list[Axis] = []
    for a, order in sorted(axes, key=lambda t: -t[1]):
        if all(abs(float(a @ b)) <= 0.995 for b, _ in merged):
            merged.append((a, order))
    return merged


def _inertia_tensor(masses: np.ndarray, x: np.ndarray) -> np.ndarray:
    r2 = (x * x).sum(axis=1)
    return (masses * r2).sum() * np.eye(3) - (masses[:, None] * x).T @ x


def find_point_group(symbols: Sequence[str], coords: ArrayLike, tol: float = 0.1) -> PointGroup:
    """Detect the point group of a structure.

    Parameters
    ----------
    symbols : sequence of str
        Element symbols.
    coords : array_like
        Coordinates, shape ``(n, 3)``.
    tol : float, optional
        Largest displacement in angstrom at which two atoms still count as
        symmetry-equivalent.

    Returns
    -------
    PointGroup
        The group, its operations and the centre of mass.
    """
    xyz = np.asarray(coords, dtype=float)
    n = len(symbols)
    if n == 0:
        return PointGroup("C1", [np.eye(3)])
    if n == 1:
        return PointGroup("Kh", [np.eye(3)], xyz[0].copy())

    an = _Analyzer(symbols, xyz, tol)
    if n == 2 or np.linalg.matrix_rank(an.x, tol=tol) <= 1 or _is_linear(an.x, tol):
        if an.is_symmetry(-np.eye(3)):
            return PointGroup(f"D{INF}h", [np.eye(3), -np.eye(3)], an.center)
        return PointGroup(f"C{INF}v", [np.eye(3)], an.center)

    moments, principal = np.linalg.eigh(_inertia_tensor(an.masses, an.x))
    scale = max(moments[-1], 1e-6)
    rel_tol = 0.03 + 0.2 * tol
    low_pair = abs(moments[0] - moments[1]) / scale < rel_tol
    high_pair = abs(moments[1] - moments[2]) / scale < rel_tol
    unique = None
    if low_pair and high_pair:
        rotor = _Rotor.SPHERICAL
    elif low_pair or high_pair:
        rotor = _Rotor.SYMMETRIC
        unique = principal[:, 2] if low_pair else principal[:, 0]
    else:
        rotor = _Rotor.ASYMMETRIC

    axes = []
    for a in an.candidate_axes(principal, rotor, unique):
        orders = [k for k in range(2, MAX_AXIS_ORDER + 1) if an.is_symmetry(_rotation(a, k))]
        if orders:
            axes.append((a, orders[-1]))
    axes = _merge_axes(axes)
    mirrors = [m for m in an.candidate_normals(principal, rotor, axes) if an.is_symmetry(_reflection(m))]
    inversion = an.is_symmetry(-np.eye(3))
    impropers = [(a, 2 * order) for a, order in axes if an.is_symmetry(_improper(a, 2 * order))]

    name, main, second = _classify(axes, mirrors, inversion, impropers)
    ops = _group_operations(name, main, second, axes, mirrors, inversion)
    return PointGroup(name, ops, an.center)


def _in_plane(v: np.ndarray, normal: np.ndarray) -> np.ndarray:
    return unit_vector(v - (v @ normal) * normal)


def _group_operations(name: str, main: np.ndarray | None, second: np.ndarray | None, axes: list[Axis],
                      mirrors: list[np.ndarray], inversion: bool) -> list[np.ndarray]:
    """Generate the full set of operations from a few exact generators."""
    if name == "C1":
        return [np.eye(3)]
    if name == "Cs":
        return _close_group([_reflection(mirrors[0])])
    if name == "Ci":
        return _close_group([-np.eye(3)])
    if main is None:
        # Cubic and icosahedral groups: use every element found, with a loose
        # tolerance because the generators are only approximately exact.
        gens = [_rotation(a, order) for a, order in axes]
        gens += [_reflection(m) for m in mirrors[:3]]
        if inversion:
            gens.append(-np.eye(3))
        return _close_group(gens, tol=0.05, limit=120)

    n = int("".join(ch for ch in name[1:] if ch.isdigit()))
    if name.startswith("S"):
        return _close_group([_improper(main, n)])
    gens = [_rotation(main, n)]
    if name.startswith("D"):
        assert second is not None
        gens.append(rotation_matrix(_in_plane(second, main), math.pi))
        if name.endswith("h"):
            gens.append(_reflection(main))
        elif name.endswith("d"):
            sigma = next(m for m in mirrors if _perpendicular(m, main))
            gens.append(_reflection(_in_plane(sigma, main)))
    elif name.endswith("v"):
        assert second is not None
        gens.append(_reflection(_in_plane(second, main)))
    elif name.endswith("h"):
        gens.append(_reflection(main))
    return _close_group(gens)


def symmetrize(symbols: Sequence[str], coords: ArrayLike, tol: float = 0.1,
               passes: int = 3) -> tuple[np.ndarray, str]:
    """Move atoms onto exact symmetry positions of the detected point group.

    Each pass averages every atom over the images of its partners under all
    operations of the group.

    Parameters
    ----------
    symbols : sequence of str
        Element symbols.
    coords : array_like
        Coordinates, shape ``(n, 3)``.
    tol : float, optional
        Tolerance for detecting the group, in angstrom.
    passes : int, optional
        Maximum number of averaging passes.

    Returns
    -------
    coords : numpy.ndarray
        Symmetrised coordinates.
    name : str
        Point group of the result.
    """
    xyz = np.asarray(coords, dtype=float).copy()
    pg = find_point_group(symbols, xyz, tol)
    if len(symbols) < 2 or pg.order <= 1:
        return xyz, pg.name
    for _ in range(passes):
        pg = find_point_group(symbols, xyz, tol)
        if pg.order <= 1:
            break
        an = _Analyzer(symbols, xyz, max(tol, 1e-3))
        total = sum((an.x[an.permutation(M)] @ M for M in pg.operations), np.zeros_like(an.x))
        xyz = total / pg.order + an.center
    return xyz, find_point_group(symbols, xyz, tol).name
