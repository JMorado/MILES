r"""Electrostatic potential of point charges, for colouring molecular surfaces.

The potential at a point :math:`\mathbf r` of charges :math:`q_i` at
:math:`\mathbf r_i` is

.. math::

   \phi(\mathbf r) = k \sum_i \frac{q_i}{\varepsilon(d_i)\, d_i},
   \qquad d_i = |\mathbf r - \mathbf r_i|,

in kcal/(mol e), with :math:`k` = 332.06 kcal angstrom/(mol e^2). Two
dielectric models are offered:

``VACUUM`` (:math:`\varepsilon = 1`)
    The potential of the charges in the gas phase, the quantity that
    quantum chemistry programs map onto a density surface. Used for small
    molecules.
``DISTANCE`` (:math:`\varepsilon = 4 d`)
    A distance-dependent dielectric that stands in for the screening by
    water and the protein itself. In vacuum the potential of thousands of
    charges is dominated by the net charge and far-away groups; with
    :math:`\varepsilon = 4d` it reflects the groups near each part of the
    surface. It is the default of Coulombic surface colouring in UCSF
    Chimera and ChimeraX and is used for proteins and nucleic acids.

Following the same convention, the potential of a surface is evaluated
1.4 angstrom outside it along the normals, where the centre of a
water-sized probe touching the surface sits: what an approaching molecule
feels, rather than a value dominated by the nearest atom's charge.

Sums are evaluated with matrix products in chunks. For large surfaces the
potential is computed exactly at one point per small cell of the surface
and spread to the other points by Gaussian-weighted averaging, which keeps
a protein of several thousand atoms well under a second.
"""

from enum import StrEnum

import numpy as np
from numpy.typing import ArrayLike

from miles.core.molecule import Molecule
from miles.core.neighbors import cross_pairs
from miles.core.protein import has_polymer
from miles.core.surface import SurfaceMesh
from miles.core.units import COULOMB

PROBE_OFFSET = 1.4
"""Distance in angstrom outside the surface at which its potential is evaluated."""

MIN_DISTANCE = 0.5
"""Distances below this (angstrom) are raised to it, so a point on an atom gives a finite value."""

CHUNK_ELEMENTS = 8_000_000
"""Point-charge pairs evaluated at once."""

EXACT_LIMIT = 3e8
"""Above this many point-charge pairs, the potential is computed at fewer points and interpolated."""



class Dielectric(StrEnum):
    """How the charges are screened."""

    VACUUM = "vacuum"       # epsilon = 1
    DISTANCE = "distance"   # epsilon = 4 d


def default_dielectric(mol: Molecule) -> Dielectric:
    """``DISTANCE`` for structures with protein or nucleic acid chains, ``VACUUM`` otherwise."""
    return Dielectric.DISTANCE if has_polymer(mol) else Dielectric.VACUUM


def electrostatic_potential(points: ArrayLike, coords: ArrayLike, charges: ArrayLike,
                            dielectric: Dielectric = Dielectric.VACUUM) -> np.ndarray:
    """Potential of point charges at some points, summed over every charge.

    Parameters
    ----------
    points : array_like
        Where to evaluate, in angstrom, shape ``(n, 3)``.
    coords : array_like
        Charge positions in angstrom, shape ``(m, 3)``.
    charges : array_like
        Charges in elementary charges, shape ``(m,)``.
    dielectric : Dielectric, optional
        Screening model.

    Returns
    -------
    numpy.ndarray
        Potential in kcal/(mol e), shape ``(n,)``.
    """
    p = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    c = np.asarray(coords, dtype=np.float64).reshape(-1, 3)
    q = np.asarray(charges, dtype=np.float64).reshape(-1)
    keep = q != 0.0
    c, q = c[keep], q[keep]
    result = np.zeros(len(p))
    if not len(p) or not len(c):
        return result
    centre = c.mean(axis=0)            # float32 sums stay accurate near the origin
    c32 = (c - centre).astype(np.float32)
    p32 = (p - centre).astype(np.float32)
    q32 = q.astype(np.float32)
    cc = (c32 * c32).sum(axis=1)
    minimum = np.float32(MIN_DISTANCE ** 2)
    step = max(1, CHUNK_ELEMENTS // len(c32))
    for start in range(0, len(p32), step):
        block = p32[start:start + step]
        d2 = block @ (-2.0 * c32.T)
        d2 += (block * block).sum(axis=1)[:, None]
        d2 += cc[None, :]
        np.maximum(d2, minimum, out=d2)
        if dielectric == Dielectric.DISTANCE:
            np.reciprocal(d2, out=d2)
            d2 *= np.float32(0.25)        # epsilon = 4 d gives q / (4 d^2)
        else:
            np.sqrt(d2, out=d2)
            np.reciprocal(d2, out=d2)
        result[start:start + step] = d2 @ q32
    return result * COULOMB


def surface_potential(mol: Molecule, surface: SurfaceMesh, charges: ArrayLike,
                      dielectric: Dielectric | None = None, atoms: ArrayLike | None = None) -> np.ndarray:
    """Potential at each vertex of a surface, evaluated :data:`PROBE_OFFSET` outside it.

    Parameters
    ----------
    mol : Molecule
        Structure carrying the charges.
    surface : SurfaceMesh
        Surface with outward normals.
    charges : array_like
        Charge of every atom of ``mol``, shape ``(n_atoms,)``.
    dielectric : Dielectric, optional
        Screening model (default: :func:`default_dielectric`).
    atoms : array_like of int, optional
        Atoms whose charges count (default: all), e.g. the atoms the
        surface encloses.

    Returns
    -------
    numpy.ndarray
        Potential in kcal/(mol e), one value per vertex.
    """
    if surface.empty:
        return np.zeros(len(surface.vertices))
    model = default_dielectric(mol) if dielectric is None else dielectric
    q = np.asarray(charges, dtype=np.float64).reshape(-1)
    index = np.arange(len(mol)) if atoms is None else np.asarray(atoms, dtype=int).reshape(-1)
    coords, q = mol.coords[index], q[index]
    points = surface.vertices.astype(np.float64) + PROBE_OFFSET * surface.normals.astype(np.float64)
    count = int(np.count_nonzero(q))
    if len(points) * count <= EXACT_LIMIT:
        return electrostatic_potential(points, coords, q, model)
    return _interpolated(points, coords, q, model, max(1, int(EXACT_LIMIT // max(count, 1))))


def _interpolated(points: np.ndarray, coords: np.ndarray, charges: np.ndarray, dielectric: Dielectric,
                  target: int) -> np.ndarray:
    """Exact potential at about ``target`` representative points, spread to all points."""
    size = 1.0
    while True:
        cells = np.floor(points / size).astype(np.int64)
        _, owner, counts = np.unique(cells, axis=0, return_inverse=True, return_counts=True)
        owner = owner.reshape(-1)
        if len(counts) <= target:
            break
        size *= 1.25
    representatives = np.stack([np.bincount(owner, points[:, d]) / counts for d in range(3)], axis=1)
    exact = electrostatic_potential(representatives, coords, charges, dielectric)
    i, j = cross_pairs(points, representatives, 1.5 * size)
    weights = np.exp(-((points[i] - representatives[j]) ** 2).sum(axis=1) / (0.5 * size * size))
    total = np.bincount(i, weights, minlength=len(points))
    result = np.bincount(i, weights * exact[j], minlength=len(points))
    lonely = total < 1e-12
    result[~lonely] /= total[~lonely]
    result[lonely] = exact[owner[lonely]]
    return result
