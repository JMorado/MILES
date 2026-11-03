"""Volumetric data on a grid (orbitals, densities) and its isosurfaces.

Quantum chemistry programs write values on a regular grid, typically in
Gaussian cube files: molecular orbitals, which change sign, and electron
densities, which do not. An isosurface joins the points where the values
equal a chosen level. For orbitals it is drawn at ``+level`` and ``-level``,
giving the two lobes of opposite phase. A second grid can colour an
isosurface, for example the electrostatic potential on a density surface.
"""

import itertools
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from miles.core.surface import SurfaceMesh, surface_nets

ORBITAL_LEVEL = 0.02
"""Default isovalue of orbitals (atomic units), as most viewers use."""

DENSITY_LEVEL = 0.002
"""Default isovalue of densities (atomic units), close to the van der Waals surface."""

SIGNED_FRACTION = 0.05
"""Values of both signs count when each reaches this fraction of the largest magnitude."""

RANGE_PERCENTILE = 95.0
"""A colour scale suggested for mapped values covers this percentage of them."""


@dataclass
class Volume:
    """Values on a regular grid.

    Attributes
    ----------
    name : str
        Description, e.g. the comment line of a cube file.
    origin : numpy.ndarray
        Position in angstrom of grid point ``(0, 0, 0)``, shape ``(3,)``.
    axes : numpy.ndarray
        Step vectors in angstrom, one row per grid axis, shape ``(3, 3)``.
    values : numpy.ndarray
        Grid values, shape ``(nx, ny, nz)``, float32; point ``(i, j, k)`` lies
        at ``origin + i * axes[0] + j * axes[1] + k * axes[2]``.
    both_signs : bool or None
        True for densities that take both signs by nature, such as a spin
        density or a difference of densities, even where one sign is small;
        None leaves :attr:`signed` to the values.
    """

    name: str
    origin: np.ndarray
    axes: np.ndarray
    values: np.ndarray
    both_signs: bool | None = None

    @property
    def signed(self) -> bool:
        """Whether the values take both signs, as an orbital's or a spin density's do.

        Set by :attr:`both_signs` when given; otherwise both signs must each
        reach :data:`SIGNED_FRACTION` of the largest magnitude.
        """
        if self.both_signs is not None:
            return self.both_signs
        peak = float(np.abs(self.values).max()) if self.values.size else 0.0
        if peak == 0.0:
            return False
        threshold = SIGNED_FRACTION * peak
        return bool(self.values.max() > threshold and self.values.min() < -threshold)

    def suggested_level(self) -> float:
        """Default isovalue for drawing.

        :data:`ORBITAL_LEVEL` for signed data and :data:`DENSITY_LEVEL`
        otherwise (also for densities marked with :attr:`both_signs`), the
        values most viewers use for data in atomic units.
        Grids whose values are all much smaller (other units, diffuse
        orbitals) get a quarter of their largest magnitude instead, so the
        surface never disappears.

        Returns
        -------
        float
            A positive isovalue (0 for an all-zero grid).
        """
        peak = float(np.abs(self.values).max()) if self.values.size else 0.0
        typical = ORBITAL_LEVEL if self.signed and not self.both_signs else DENSITY_LEVEL
        return min(typical, 0.25 * peak)

    def sample(self, points: ArrayLike) -> np.ndarray:
        """Values at arbitrary positions, by trilinear interpolation.

        Parameters
        ----------
        points : array_like
            Positions in angstrom, shape ``(n, 3)``.

        Returns
        -------
        numpy.ndarray
            Values, shape ``(n,)``. Points outside the grid take the value at
            the nearest edge.
        """
        p = np.asarray(points, dtype=np.float64).reshape(-1, 3)
        shape = np.array(self.values.shape)
        frac = np.clip((p - self.origin) @ np.linalg.inv(self.axes), 0, shape - 1)
        low = np.minimum(np.floor(frac).astype(int), np.maximum(shape - 2, 0))
        high = np.minimum(low + 1, shape - 1)
        t = frac - low
        result = np.zeros(len(p))
        for corner in itertools.product((0, 1), repeat=3):
            index = tuple(np.where(c, high[:, d], low[:, d]) for d, c in enumerate(corner))
            weight = np.prod([t[:, d] if c else 1 - t[:, d] for d, c in enumerate(corner)], axis=0)
            result += weight * self.values[index]
        return result


def isosurface(volume: Volume, level: float) -> SurfaceMesh:
    """Surface where the values of a volume equal ``level``.

    Parameters
    ----------
    volume : Volume
        The data.
    level : float
        Iso-value. For ``level >= 0`` the enclosed region is where the values
        exceed it; for ``level < 0`` it is where they fall below it (a
        negative lobe).

    Returns
    -------
    SurfaceMesh
        Closed surface with normals pointing out of the enclosed region; the
        grid is padded so lobes that reach its edge are closed too. Vertices
        belong to no atom, so ``atoms`` is -1.
    """
    sign = 1.0 if level >= 0 else -1.0
    values = sign * np.asarray(volume.values, dtype=np.float32)
    target = sign * level
    outside = min(float(values.min()), target) - max(abs(target), 1e-6)
    padded = np.pad(values, 1, mode="constant", constant_values=outside)
    axes = np.asarray(volume.axes, dtype=np.float64)
    origin = np.asarray(volume.origin, dtype=np.float64) - axes.sum(axis=0)
    mesh = surface_nets(padded, target, origin, axes)
    if mesh is None:
        return SurfaceMesh()
    vertices, normals, faces = mesh
    return SurfaceMesh(vertices=vertices.astype(np.float32), normals=normals.astype(np.float32),
                       faces=faces.astype(np.uint32), atoms=np.full(len(vertices), -1, dtype=int))


def signed_isosurfaces(volume: Volume, level: float) -> tuple[SurfaceMesh, SurfaceMesh]:
    """Positive and negative lobes of signed data, at ``+level`` and ``-level``.

    Parameters
    ----------
    volume : Volume
        The data, e.g. an orbital.
    level : float
        Iso-value; its sign is ignored.

    Returns
    -------
    tuple of SurfaceMesh
        The positive lobe and the negative lobe, each with normals pointing
        out of the lobe.
    """
    magnitude = abs(level)
    return isosurface(volume, magnitude), isosurface(volume, -magnitude)


def mapped_isosurface(volume: Volume, level: float, other: Volume) -> tuple[SurfaceMesh, np.ndarray]:
    """Isosurface of one volume with the values of another at its vertices.

    Parameters
    ----------
    volume : Volume
        Data whose isosurface is drawn, e.g. the electron density.
    level : float
        Iso-value; the surface encloses the values above ``abs(level)``.
    other : Volume
        Data that colours the surface, e.g. the electrostatic potential.
        It may use a different grid, as long as both cover the molecule.

    Returns
    -------
    tuple of (SurfaceMesh, numpy.ndarray)
        The surface and the values of ``other`` at its vertices.
    """
    surface = isosurface(volume, abs(level))
    if surface.empty:
        return surface, np.zeros(0)
    return surface, other.sample(surface.vertices)


def symmetric_range(values: ArrayLike) -> float:
    """Half-width of a colour scale centred on zero that suits these values.

    Parameters
    ----------
    values : array_like
        Values to be coloured, e.g. from :func:`mapped_isosurface`.

    Returns
    -------
    float
        The :data:`RANGE_PERCENTILE` percentile of the magnitudes, so a few
        extreme values do not wash out the rest; 1 when there are no values.
    """
    magnitudes = np.abs(np.asarray(values, dtype=np.float64))
    if magnitudes.size == 0:
        return 1.0
    return float(np.percentile(magnitudes, RANGE_PERCENTILE)) or float(magnitudes.max()) or 1.0
