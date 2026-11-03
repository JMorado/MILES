"""Vector maths, internal coordinates and ideal coordination geometries.

Angles passed to and returned from the public measurement functions are in
degrees; rotation helpers take radians.
"""

import itertools
import math
from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

Vector = np.ndarray


class Geometry(StrEnum):
    """Ideal arrangement of the bonds around an atom (VSEPR shapes)."""

    LINEAR = "linear"
    TRIGONAL = "trigonal"
    TETRAHEDRAL = "tetrahedral"
    SQUARE_PLANAR = "square_planar"
    TBP = "tbp"
    OCTAHEDRAL = "octahedral"

    @property
    def label(self) -> str:
        """Human-readable name, e.g. ``"trigonal bipyramidal"``."""
        return GEOMETRY_NAMES[self]


def unit_vector(v: ArrayLike) -> Vector:
    """Return the unit vector along ``v`` (the zero vector stays zero).

    Parameters
    ----------
    v : array_like
        Vector of shape ``(3,)``.

    Returns
    -------
    numpy.ndarray
        Normalised copy of ``v``.
    """
    a = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(a))
    return np.zeros(3) if n < 1e-12 else a / n


def any_perpendicular(v: ArrayLike) -> Vector:
    """Return a unit vector perpendicular to ``v``."""
    u = unit_vector(v)
    trial = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    return unit_vector(np.cross(u, trial))


def rotation_matrix(axis: ArrayLike, angle: float) -> np.ndarray:
    """Right-handed rotation about an axis (Rodrigues' formula).

    Parameters
    ----------
    axis : array_like
        Rotation axis, need not be normalised.
    angle : float
        Rotation angle in radians.

    Returns
    -------
    numpy.ndarray
        Rotation matrix, shape ``(3, 3)``.
    """
    x, y, z = unit_vector(axis)
    c, s = math.cos(angle), math.sin(angle)
    t = 1.0 - c
    return np.array([
        [c + x * x * t, x * y * t - z * s, x * z * t + y * s],
        [y * x * t + z * s, c + y * y * t, y * z * t - x * s],
        [z * x * t - y * s, z * y * t + x * s, c + z * z * t],
    ])


def rotation_between(a: ArrayLike, b: ArrayLike) -> np.ndarray:
    """Smallest rotation taking direction ``a`` onto direction ``b``.

    Parameters
    ----------
    a, b : array_like
        Directions, need not be normalised.

    Returns
    -------
    numpy.ndarray
        Rotation matrix ``R`` with ``R @ unit_vector(a) == unit_vector(b)``.
    """
    ua, ub = unit_vector(a), unit_vector(b)
    c = float(np.dot(ua, ub))
    if c > 1.0 - 1e-10:
        return np.eye(3)
    if c < -1.0 + 1e-10:
        return rotation_matrix(any_perpendicular(ua), math.pi)
    return rotation_matrix(np.cross(ua, ub), math.acos(max(-1.0, min(1.0, c))))


def signed_angle(u: ArrayLike, v: ArrayLike, axis: ArrayLike) -> float:
    """Angle in radians that turns ``u`` onto ``v`` about ``axis``.

    Both vectors are assumed perpendicular to the axis.
    """
    ua, va = np.asarray(u, dtype=float), np.asarray(v, dtype=float)
    return math.atan2(float(np.dot(np.cross(ua, va), unit_vector(axis))), float(np.dot(ua, va)))


def distance(p: ArrayLike, q: ArrayLike) -> float:
    """Distance between two points."""
    return float(np.linalg.norm(np.asarray(p) - np.asarray(q)))


def angle(p: ArrayLike, q: ArrayLike, r: ArrayLike) -> float:
    """Angle p-q-r in degrees, with ``q`` at the vertex."""
    pa, qa, ra = (np.asarray(x, dtype=float) for x in (p, q, r))
    c = float(np.dot(unit_vector(pa - qa), unit_vector(ra - qa)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def dihedral(p: ArrayLike, q: ArrayLike, r: ArrayLike, s: ArrayLike) -> float:
    """Dihedral angle p-q-r-s in degrees, in the range (-180, 180]."""
    pa, qa, ra, sa = (np.asarray(x, dtype=float) for x in (p, q, r, s))
    b1 = unit_vector(ra - qa)
    v = (pa - qa) - np.dot(pa - qa, b1) * b1
    w = (sa - ra) - np.dot(sa - ra, b1) * b1
    return math.degrees(math.atan2(float(np.dot(np.cross(b1, v), w)), float(np.dot(v, w))))


def kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Proper rotation that best maps points ``P`` onto ``Q`` (Kabsch algorithm).

    Parameters
    ----------
    P, Q : numpy.ndarray
        Corresponding points, shape ``(n, 3)``. They are not centred here.

    Returns
    -------
    numpy.ndarray
        Rotation ``R`` minimising ``sum |R p_i - q_i|^2``, with ``det(R) = +1``.
    """
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = 1.0 if np.linalg.det(Vt.T @ U.T) >= 0 else -1.0
    return Vt.T @ np.diag([1.0, 1.0, d]) @ U.T


def place_atom(a: Vector, b: Vector, c: Vector | None, r: float, theta: float, phi: float) -> Vector:
    """Position of a new atom from internal coordinates (Z-matrix style).

    Parameters
    ----------
    a, b : numpy.ndarray
        Bonded atom and angle reference.
    c : numpy.ndarray or None
        Dihedral reference; any perpendicular direction is used when None.
    r : float
        Distance to ``a``.
    theta : float
        Angle new-a-b in degrees.
    phi : float
        Dihedral new-a-b-c in degrees.

    Returns
    -------
    numpy.ndarray
        Cartesian position of the new atom.
    """
    th, ph = math.radians(theta), math.radians(phi)
    bc = unit_vector(a - b)
    ref = any_perpendicular(bc) if c is None else b - c
    n = unit_vector(np.cross(ref, bc))
    if np.linalg.norm(n) < 1e-8:
        n = any_perpendicular(bc)
    frame = np.column_stack([bc, np.cross(n, bc), n])
    local = np.array([-r * math.cos(th), r * math.sin(th) * math.cos(ph), r * math.sin(th) * math.sin(ph)])
    return a + frame @ local


def pairwise_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Distance matrix between two point sets, shape ``(len(a), len(b))``."""
    diff = a[:, None, :] - b[None, :, :]
    return np.sqrt((diff * diff).sum(axis=-1))


_S3 = math.sqrt(3.0) / 2.0
_T1 = 2.0 * math.sqrt(2.0) / 3.0
_T2 = math.sqrt(2.0) / 3.0
_T3 = math.sqrt(6.0) / 3.0

TEMPLATES: dict[str, np.ndarray] = {
    Geometry.LINEAR: np.array([[1, 0, 0], [-1, 0, 0]], dtype=float),
    Geometry.TRIGONAL: np.array([[1, 0, 0], [-0.5, _S3, 0], [-0.5, -_S3, 0]]),
    Geometry.TETRAHEDRAL: np.array([[0, 0, 1], [_T1, 0, -1 / 3], [-_T2, _T3, -1 / 3], [-_T2, -_T3, -1 / 3]]),
    Geometry.SQUARE_PLANAR: np.array([[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]], dtype=float),
    Geometry.TBP: np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0], [-0.5, _S3, 0], [-0.5, -_S3, 0]]),
    Geometry.OCTAHEDRAL: np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]],
                                  dtype=float),
}
"""Unit vectors of the bond directions of each geometry."""

GEOMETRY_NAMES: dict[str, str] = {
    Geometry.LINEAR: "linear",
    Geometry.TRIGONAL: "trigonal planar",
    Geometry.TETRAHEDRAL: "tetrahedral",
    Geometry.SQUARE_PLANAR: "square planar",
    Geometry.TBP: "trigonal bipyramidal",
    Geometry.OCTAHEDRAL: "octahedral",
}

IDEAL_ANGLES: dict[str, tuple[float, ...]] = {
    Geometry.LINEAR: (180.0,),
    Geometry.TRIGONAL: (120.0,),
    Geometry.TETRAHEDRAL: (109.4712,),
    Geometry.SQUARE_PLANAR: (90.0, 180.0),
    Geometry.TBP: (90.0, 120.0, 180.0),
    Geometry.OCTAHEDRAL: (90.0, 180.0),
}
"""Ideal angles in degrees between two bonds of each geometry."""


def geometry_for_coordination(n: int) -> Geometry | None:
    """Default geometry for an atom with ``n`` bonds (None above six)."""
    return {1: Geometry.LINEAR, 2: Geometry.LINEAR, 3: Geometry.TRIGONAL, 4: Geometry.TETRAHEDRAL,
            5: Geometry.TBP, 6: Geometry.OCTAHEDRAL}.get(n)


def free_sites(geom: str, neighbor_dirs: list[Vector]) -> list[Vector]:
    """Template directions left free once existing bonds are matched onto a geometry.

    The existing bond directions are fitted onto the template by trying every
    assignment of template vectors and keeping the Kabsch rotation with the
    smallest error.

    Parameters
    ----------
    geom : str
        A :class:`Geometry` value.
    neighbor_dirs : list of numpy.ndarray
        Unit vectors of the bonds the atom already has.

    Returns
    -------
    list of numpy.ndarray
        Unit vectors of the unused positions. With a single existing bond the
        result can still be rotated freely about that bond; callers choose the
        rotation.
    """
    template = TEMPLATES[geom]
    k = len(neighbor_dirs)
    if k == 0:
        return [t.copy() for t in template]
    if k >= len(template):
        return []
    nb = np.array(neighbor_dirs, dtype=float)
    if k == 1:
        R = rotation_between(template[0], nb[0])
        return [R @ t for t in template[1:]]
    # every assignment of template vectors to the bonds, fitted at once (Kabsch on stacked matrices)
    combos = list(itertools.permutations(range(len(template)), k))
    P = template[np.array(combos)]                          # (assignments, k, 3)
    U, _, Vt = np.linalg.svd(np.einsum("aki,kj->aij", P, nb))
    V, Ut = np.swapaxes(Vt, 1, 2), np.swapaxes(U, 1, 2)
    flip = np.where(np.linalg.det(V @ Ut) >= 0, 1.0, -1.0)
    diagonal = np.ones((len(combos), 3))
    diagonal[:, 2] = flip
    R = V @ (diagonal[:, :, None] * Ut)   # V diag(1, 1, flip) U^T for each assignment
    errors = ((np.einsum("aki,aji->akj", P, R) - nb) ** 2).sum(axis=(1, 2))
    best = 0
    for a in range(1, len(combos)):   # the first assignment that is better by more than 1e-9 wins
        if errors[a] < errors[best] - 1e-9:
            best = a
    return [R[best] @ template[i] for i in range(len(template)) if i not in combos[best]]
