"""Colouring molecular surfaces and preparing them for drawing."""

from dataclasses import dataclass

import numpy as np

from miles.core.display import AtomDisplay, representations
from miles.core.molecule import Molecule, Representation
from miles.core.potential import surface_potential
from miles.core.surface import SurfaceColoring, SurfaceMesh, atom_property, molecular_surface
from miles.core.volume import Volume, isosurface, mapped_isosurface, signed_isosurfaces, symmetric_range
from miles.render.cartoon import Mesh
from miles.render.theme import Theme

SEE_THROUGH_ALPHA = 0.45

# Diverging scales: value -1 at the first colour, 0 at the middle, +1 at the last.
HYDROPHOBICITY_COLORS = np.array([[0.22, 0.56, 0.84], [0.96, 0.96, 0.96], [0.90, 0.62, 0.25]], np.float32)
CHARGE_COLORS = np.array([[0.88, 0.25, 0.25], [0.96, 0.96, 0.96], [0.25, 0.42, 0.90]], np.float32)
MAPPED_COLORS = CHARGE_COLORS   # grids colouring a surface, e.g. a potential: red negative, blue positive
POSITIVE_LOBE = np.array([0.29, 0.53, 0.91], np.float32)
NEGATIVE_LOBE = np.array([0.91, 0.51, 0.23], np.float32)
DENSITY_COLOR = np.array([0.45, 0.66, 0.92], np.float32)
UNIFORM_LIGHT = np.array([0.82, 0.84, 0.88], np.float32)
UNIFORM_DARK = np.array([0.60, 0.63, 0.68], np.float32)


@dataclass(frozen=True)
class SurfaceRequest:
    """What surface a view shows.

    Attributes
    ----------
    atoms : tuple of int or None
        Atoms enclosed (None: every atom that is drawn, so hidden water is left out).
    coloring : SurfaceColoring
        Colouring scheme.
    see_through : bool
        Draw the surface partly transparent.
    """

    atoms: tuple[int, ...] | None = None
    coloring: SurfaceColoring = SurfaceColoring.ELEMENT
    see_through: bool = False

    @property
    def alpha(self) -> float:
        """Opacity of the surface."""
        return SEE_THROUGH_ALPHA if self.see_through else 1.0


def _diverging(values: np.ndarray, scale: np.ndarray) -> np.ndarray:
    v = np.clip(values, -1.0, 1.0)[:, None]
    low = scale[1] + (scale[0] - scale[1]) * np.clip(-v, 0.0, 1.0)
    high = scale[1] + (scale[2] - scale[1]) * np.clip(v, 0.0, 1.0)
    return np.where(v < 0, low, high).astype(np.float32)


def surface_colors(mol: Molecule, surface: SurfaceMesh, coloring: SurfaceColoring,
                   theme: Theme) -> np.ndarray:
    """Colour of each surface vertex, from the atom it belongs to.

    Parameters
    ----------
    mol : Molecule
        Structure the surface encloses.
    surface : SurfaceMesh
        The surface.
    coloring : SurfaceColoring
        Element colours, hydrophobicity (blue polar to orange apolar), charge
        (red negative to blue positive) or one neutral colour. The
        electrostatic potential is not a property of atoms; it is coloured
        by :func:`surface_mesh`, and here gets the neutral colour.
    theme : Theme
        Supplies element colours and the neutral colour's brightness.

    Returns
    -------
    numpy.ndarray
        RGB colours, shape ``(n_vertices, 3)``.
    """
    n = len(surface.vertices)
    if coloring in (SurfaceColoring.UNIFORM, SurfaceColoring.POTENTIAL):
        return np.tile(UNIFORM_DARK if theme.dark else UNIFORM_LIGHT, (n, 1))
    if coloring == SurfaceColoring.ELEMENT:
        return theme.element_colors()[mol.numbers][surface.atoms]
    scale = HYDROPHOBICITY_COLORS if coloring == SurfaceColoring.HYDROPHOBICITY else CHARGE_COLORS
    return _diverging(atom_property(mol, coloring)[surface.atoms], scale)


@dataclass(frozen=True)
class VolumeRequest:
    """Which isosurface of a document's volumetric data a view shows.

    Attributes
    ----------
    index : int
        Volume in ``JobInfo.volumes``.
    level : float
        Isovalue; signed data (orbitals) show the lobes at +level and -level.
    see_through : bool
        Draw the isosurface partly transparent.
    mapped : int or None
        Another volume in ``JobInfo.volumes`` whose values colour the
        isosurface (red negative, white zero, blue positive).
    mapped_range : float
        Values of the mapped volume at ``-mapped_range`` and ``+mapped_range``
        get the full colours; 0 picks a range from the values on the surface.
    """

    index: int
    level: float
    see_through: bool = True
    mapped: int | None = None
    mapped_range: float = 0.0

    @property
    def alpha(self) -> float:
        """Opacity of the isosurface."""
        return SEE_THROUGH_ALPHA + 0.15 if self.see_through else 1.0


def _merge(parts: list[tuple[SurfaceMesh, np.ndarray]]) -> Mesh:
    """One drawable mesh from several surfaces, each in one colour."""
    parts = [(m, c) for m, c in parts if not m.empty]
    if not parts:
        return Mesh()
    offsets = np.cumsum([0] + [len(m.vertices) for m, _ in parts[:-1]])
    vertices = np.concatenate([m.vertices for m, _ in parts]).astype(np.float32)
    colors = np.concatenate([np.tile(c, (len(m.vertices), 1)) for m, c in parts]).astype(np.float32)
    return Mesh(vertices=vertices, normals=np.concatenate([m.normals for m, _ in parts]).astype(np.float32),
                colors=colors, base_colors=colors,
                indices=np.concatenate([m.faces + o for (m, _), o in zip(parts, offsets, strict=True)]))


def volume_mesh(volume: Volume, level: float, mapped: Volume | None = None,
                mapped_range: float = 0.0) -> Mesh:
    """Isosurface of volumetric data as a drawable mesh.

    Signed data such as orbitals give two lobes, blue at +level and orange
    at -level; densities give one surface. With ``mapped``, the surface at
    +level is coloured by that volume's values instead (see
    :class:`VolumeRequest`).
    """
    if mapped is not None:
        surface, values = mapped_isosurface(volume, level, mapped)
        if surface.empty:
            return Mesh()
        colors = _diverging(values / (mapped_range or symmetric_range(values)), MAPPED_COLORS)
        return Mesh(vertices=surface.vertices, normals=surface.normals, colors=colors, base_colors=colors,
                    indices=surface.faces)
    if volume.signed:
        positive, negative = signed_isosurfaces(volume, abs(level))
        return _merge([(positive, POSITIVE_LOBE), (negative, NEGATIVE_LOBE)])
    return _merge([(isosurface(volume, level), DENSITY_COLOR)])


def surface_mesh(mol: Molecule, request: SurfaceRequest, theme: Theme,
                 charges: np.ndarray | None = None, display: AtomDisplay | None = None) -> Mesh:
    """Compute and colour the surface a view asks for, as a drawable mesh.

    Parameters
    ----------
    mol : Molecule
        Structure.
    request : SurfaceRequest
        Atoms to wrap and colouring.
    theme : Theme
        Colours.
    charges : numpy.ndarray, optional
        Charge of every atom, for the ``POTENTIAL`` colouring; the charges
        of the atoms the surface encloses make the potential (see
        :mod:`miles.core.potential`). Without them the surface gets the
        neutral colour.
    display : AtomDisplay, optional
        Representations and focus-mode hiding; when the request names no
        atoms, the surface wraps the atoms that are not hidden.

    Returns
    -------
    Mesh
        The surface; for the potential, ``value_range`` holds the potential
        in kcal/(mol e) that gets the full red and blue.
    """
    atoms = request.atoms
    if atoms is None:
        reps = representations(mol, display, Representation.DEFAULT)
        atoms = tuple(np.flatnonzero(reps != Representation.HIDDEN).tolist())
    surface = molecular_surface(mol, atoms)
    if surface.empty:
        return Mesh()
    value_range = 0.0
    if request.coloring == SurfaceColoring.POTENTIAL and charges is not None:
        values = surface_potential(mol, surface, charges, atoms=atoms)
        value_range = symmetric_range(values)
        colors = _diverging(values / value_range, CHARGE_COLORS)
    else:
        colors = surface_colors(mol, surface, request.coloring, theme)
    return Mesh(vertices=surface.vertices, normals=surface.normals, colors=colors, base_colors=colors,
                indices=surface.faces, value_range=value_range)
