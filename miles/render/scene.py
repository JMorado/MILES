"""What to draw: display style, per-atom representations and the scene arrays.

A :class:`Scene` is a set of flat float32 arrays that both renderers consume:
the OpenGL renderer uploads them as instance buffers, the QPainter fallback
walks them back to front.
"""

import functools
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import IntEnum, StrEnum
from typing import Any

import numpy as np

from miles.core.display import AtomDisplay, representations
from miles.core.elements import COVALENT_RADII, VDW_RADII, element
from miles.core.molecule import Molecule, Representation
from miles.core.protein import POLYMER_RESIDUES
from miles.render.cartoon import CartoonColor, Mesh, cartoon_mesh, cartoon_residues
from miles.render.surface import CHARGE_COLORS
from miles.render.theme import Theme

REPRESENTATION_NAMES: dict[Representation, str] = {
    Representation.BALL: "Ball and Stick",
    Representation.TUBE: "Tube",
    Representation.WIRE: "Wire",
    Representation.SPACEFILL: "Space Fill",
    Representation.CARTOON: "Cartoon",
    Representation.HIDDEN: "Hidden",
}
"""Menu names of the representations, in menu order."""

CHARGE_COLOR_LIMIT = 0.5
"""Smallest charge (elementary charges) drawn in full red or blue."""

TUBE_RADIUS = 0.14
WIRE_RADIUS = 0.05
WIRE_WIDTH_PX = 1.4
BOND_RADIUS = 0.075


class Highlight(IntEnum):
    """Highlight drawn around an atom; the value indexes the renderer's colour table."""

    NONE = 0
    SELECTED = 1
    HOVER = 2
    ATTACH = 3
    DANGER = 4


class LabelMode(StrEnum):
    """Text drawn next to each atom."""

    NONE = "none"
    SYMBOL = "symbol"
    NUMBER = "number"
    BOTH = "both"


class CylinderKind(IntEnum):
    """Last column of the cylinder array."""

    SOLID = 0
    DASHED = 1
    TRACE = 2  # cartoon backbone in the QPainter fallback, drawn without outlines


@dataclass
class Style:
    """Display options of one view.

    Attributes
    ----------
    mode : Representation
        Representation of atoms that have none of their own.
    show_hydrogens : bool
        Draw hydrogen atoms.
    labels : LabelMode
        Text next to the atoms.
    depth_fade : bool
        Fade distant atoms into the background.
    cartoon_color : CartoonColor
        Colouring of protein and nucleic acid cartoons.
    stereo : bool
        Label stereocentres R/S and double bonds E/Z.
    """

    mode: Representation = Representation.BALL
    show_hydrogens: bool = True
    labels: LabelMode = LabelMode.NONE
    depth_fade: bool = True
    cartoon_color: CartoonColor = CartoonColor.SECONDARY_STRUCTURE
    stereo: bool = False


def charge_colors(charges: np.ndarray) -> np.ndarray:
    """Atom colours for partial charges: red negative, white neutral, blue positive.

    The scale is symmetric about zero and reaches full colour at the largest
    charge in the molecule, or at :data:`CHARGE_COLOR_LIMIT` if that is larger,
    so nearly non-polar molecules stay pale.

    Parameters
    ----------
    charges : numpy.ndarray
        Charge of each atom in elementary charges.

    Returns
    -------
    numpy.ndarray
        RGB colours, shape ``(n_atoms, 3)``.
    """
    q = np.asarray(charges, dtype=float)
    limit = max(float(np.abs(q).max(initial=0.0)), CHARGE_COLOR_LIMIT)
    v = np.clip(q / limit, -1.0, 1.0)[:, None]
    low = CHARGE_COLORS[1] + (CHARGE_COLORS[0] - CHARGE_COLORS[1]) * np.clip(-v, 0.0, 1.0)
    high = CHARGE_COLORS[1] + (CHARGE_COLORS[2] - CHARGE_COLORS[1]) * np.clip(v, 0.0, 1.0)
    return np.where(v < 0, low, high).astype(np.float32)


@functools.cache
def atom_radius(symbol: str, rep: str) -> float:
    """Sphere radius in angstrom of an element drawn in a representation."""
    if rep == Representation.SPACEFILL:
        return element(symbol).vdw_radius
    if rep == Representation.TUBE:
        return TUBE_RADIUS
    if rep == Representation.WIRE:
        return WIRE_RADIUS
    return 0.16 + 0.26 * element(symbol).covalent_radius


def atom_reps(mol: Molecule, style: Style, display: AtomDisplay | None = None) -> np.ndarray:
    """Representation of every atom (its own, or the view default), as an array of strings.

    Parameters
    ----------
    mol : Molecule
        Structure.
    style : Style
        Display options; ``style.mode`` is the default representation.
    display : AtomDisplay, optional
        Each atom's own representation and focus-mode hiding; hidden atoms
        count as ``HIDDEN``. Without it every atom takes the view's style.

    Returns
    -------
    numpy.ndarray
        One :class:`Representation` value per atom. Atoms outside a polymer
        residue cannot be part of a cartoon and are drawn as tubes instead.
    """
    reps = representations(mol, display, style.mode)
    cartoon = reps == Representation.CARTOON
    if cartoon.any():
        reps[cartoon & ~mol.mask("resname", POLYMER_RESIDUES)] = str(Representation.TUBE)
    return reps


def visible_atoms(mol: Molecule, style: Style, reps: np.ndarray | None = None,
                  display: AtomDisplay | None = None) -> np.ndarray:
    """Boolean mask of the atoms drawn as spheres (not hidden, not in a cartoon)."""
    reps = atom_reps(mol, style, display) if reps is None else np.asarray(reps, dtype=object)
    shown = (reps != Representation.HIDDEN) & (reps != Representation.CARTOON)
    if not style.show_hydrogens:
        shown &= mol.numbers != 1
    return np.asarray(shown, dtype=bool)


def sphere_radii(numbers: np.ndarray, reps: np.ndarray) -> np.ndarray:
    """Drawing radius of each atom in angstrom; wire atoms get minus their width in pixels."""
    return np.where(reps == Representation.SPACEFILL, VDW_RADII[numbers],
                    np.where(reps == Representation.TUBE, TUBE_RADIUS,
                             np.where(reps == Representation.WIRE, -WIRE_WIDTH_PX,
                                      0.16 + 0.26 * COVALENT_RADII[numbers])))


def _bond_lines(order: float) -> list[tuple[float, float, CylinderKind]]:
    """Offset (angstrom), radius and kind of each cylinder drawn for a bond."""
    if order == 1.5:
        return [(-0.05, 0.062, CylinderKind.SOLID), (0.125, 0.042, CylinderKind.DASHED)]
    if order == 2.0:
        return [(-0.085, 0.055, CylinderKind.SOLID), (0.085, 0.055, CylinderKind.SOLID)]
    if order >= 3.0:
        return [(0.0, 0.045, CylinderKind.SOLID), (-0.13, 0.045, CylinderKind.SOLID),
                (0.13, 0.045, CylinderKind.SOLID)]
    if 0 < order < 1.0:
        return [(0.0, 0.05, CylinderKind.DASHED)]
    return [(0.0, BOND_RADIUS, CylinderKind.SOLID)]


@dataclass
class BondLayout:
    """The bond cylinders of a scene apart from their positions.

    Everything a cylinder row takes from the bonds, representations and
    colours is filled in here once per layout; :func:`bond_cylinders` only
    adds the positions, so a trajectory frame or a drag step costs a few
    array operations whatever the number of multiple bonds.

    Attributes
    ----------
    rows : numpy.ndarray
        Cylinder rows with radius, colours and kind set and the ends zero,
        shape ``(m, 14)``: plain bonds first, then the lines of each
        multiple bond in bond order.
    trims : numpy.ndarray
        Radii of the end atoms of each row, shape ``(m, 2)``.
    start, end : numpy.ndarray
        Atom index at each end of each row.
    shift : numpy.ndarray
        How far each row is moved off the bond axis, in angstrom.
    spread : numpy.ndarray
        Row of :attr:`a_atom` (one per spread bond) whose offset direction
        moves each cylinder row, or -1 for rows on the axis.
    a_atom, b_atom : numpy.ndarray
        The two atoms of each spread bond.
    candidates : numpy.ndarray
        For each spread bond, pairs ``(atom, neighbour)`` whose plane may
        orient the lines, in order of preference, shape ``(s, k, 2)``,
        padded with -1.
    scores : numpy.ndarray
        Preference of each candidate (higher first, the earlier one on a
        tie), -1 for padding.
    """

    rows: np.ndarray
    trims: np.ndarray
    start: np.ndarray
    end: np.ndarray
    shift: np.ndarray
    spread: np.ndarray
    a_atom: np.ndarray
    b_atom: np.ndarray
    candidates: np.ndarray
    scores: np.ndarray


def bond_layout(mol: Molecule, reps: np.ndarray, shown: np.ndarray, radii: np.ndarray,
                colors: np.ndarray) -> BondLayout:
    """Cylinder rows of the bonds between drawn atoms, without their positions (see :class:`BondLayout`)."""
    pairs, orders = mol.bond_pairs, mol.bond_orders
    rep = reps
    solid = shown & (rep != Representation.SPACEFILL)
    keep = solid[pairs[:, 0]] & solid[pairs[:, 1]] if len(pairs) else np.zeros(0, bool)
    pairs, orders = pairs[keep], orders[keep]
    i0, i1 = pairs[:, 0], pairs[:, 1]
    wire = (rep[i0] == Representation.WIRE) | (rep[i1] == Representation.WIRE)
    tube = (rep[i0] == Representation.TUBE) & (rep[i1] == Representation.TUBE)
    balls = (rep[i0] == Representation.BALL) & (rep[i1] == Representation.BALL)
    radius = np.where(wire, -WIRE_WIDTH_PX, np.where(tube, TUBE_RADIUS, BOND_RADIUS)).astype(np.float32)
    simple = ~balls | (orders == 1.0)

    k = np.nonzero(simple)[0]
    start, end = [i0[k]], [i1[k]]
    shift: list[np.ndarray] = [np.zeros(len(k))]
    spread: list[np.ndarray] = [np.full(len(k), -1)]
    radius_rows: list[np.ndarray] = [radius[k]]
    kinds: list[np.ndarray] = [np.zeros(len(k))]
    trims = [np.stack([np.abs(radii[i0[k]]), np.abs(radii[i1[k]])], axis=1)]

    # each multiple bond: its lines, and the neighbours whose plane may spread them
    symbols = mol.symbols
    a_atoms, b_atoms, candidate_lists = [], [], []
    for b in np.nonzero(~simple)[0].tolist():
        i, j, order = int(i0[b]), int(i1[b]), float(orders[b])
        lines = _bond_lines(order)
        spreads = len(lines) > 1 or lines[0][0] != 0.0
        if spreads:
            found = []
            for a, other in ((i, j), (j, i)):
                for n in mol.neighbors(a):
                    if n == other:
                        continue
                    if order == 1.5:
                        score = 2 if mol.bond_order(a, n) == 1.5 else 0
                    else:
                        score = 1 if symbols[n] != "H" else 0
                    found.append((a, n, score))
            candidate_lists.append(found)
            a_atoms.append(i)
            b_atoms.append(j)
        row = len(a_atoms) - 1 if spreads else -1
        start.append(np.full(len(lines), i))
        end.append(np.full(len(lines), j))
        shift.append(np.array([line[0] for line in lines]))
        spread.append(np.full(len(lines), row))
        radius_rows.append(np.array([line[1] for line in lines]))
        kinds.append(np.array([float(line[2]) for line in lines]))
        trims.append(np.tile(np.array([radii[i], radii[j]], np.float32), (len(lines), 1)))

    start_all, end_all = np.concatenate(start).astype(int), np.concatenate(end).astype(int)
    rows = np.zeros((len(start_all), 14), np.float32)
    rows[:, 6] = np.concatenate(radius_rows)
    rows[:, 7:10] = colors[start_all]
    rows[:, 10:13] = colors[end_all]
    rows[:, 13] = np.concatenate(kinds)
    width = max((len(found) for found in candidate_lists), default=0)
    candidates = np.full((len(candidate_lists), width, 2), -1, dtype=int)
    scores = np.full((len(candidate_lists), width), -1, dtype=int)
    for s, found in enumerate(candidate_lists):
        for c, (a, n, score) in enumerate(found):
            candidates[s, c] = (a, n)
            scores[s, c] = score
    return BondLayout(rows, np.concatenate(trims).astype(np.float32), start_all, end_all,
                      np.concatenate(shift), np.concatenate(spread).astype(int),
                      np.array(a_atoms, dtype=int), np.array(b_atoms, dtype=int), candidates, scores)


def _offset_directions(coords: np.ndarray, bonds: BondLayout) -> np.ndarray:
    """Direction in which the lines of each spread bond are moved off its axis.

    The lines lie in the plane of a neighbouring atom, preferring heavy atoms
    and, for aromatic bonds, the next ring bond, so the dashed line falls
    inside the ring. A neighbour on the bond axis does not count; with none
    left, any direction perpendicular to the bond is used.
    """
    n = len(bonds.a_atom)
    if n == 0:
        return np.zeros((0, 3))
    axis = coords[bonds.b_atom] - coords[bonds.a_atom]
    length = np.linalg.norm(axis, axis=1)
    short = length <= 1e-9
    axis = np.where(short[:, None], [1.0, 0.0, 0.0], axis / np.where(short, 1.0, length)[:, None])
    result = np.zeros((n, 3))
    best = np.full(n, -1)
    if bonds.candidates.shape[1]:
        a, k = bonds.candidates[..., 0], bonds.candidates[..., 1]
        v = coords[np.maximum(k, 0)] - coords[np.maximum(a, 0)]
        v = v - np.einsum("sci,si->sc", v, axis)[..., None] * axis[:, None, :]
        norm = np.linalg.norm(v, axis=2)
        usable = (bonds.scores >= 0) & (norm >= 1e-6)
        score = np.where(usable, bonds.scores, -1)
        best = np.argmax(score, axis=1)   # the first candidate with the highest preference
        rows = np.arange(n)
        found = score[rows, best] >= 0
        result[found] = v[rows[found], best[found]] / norm[rows[found], best[found]][:, None]
        best = np.where(found, best, -1)
    missing = best < 0
    if missing.any():
        trial = np.where((np.abs(axis[missing, 2]) < 0.9)[:, None], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
        perp = np.cross(axis[missing], trial)
        result[missing] = perp / np.linalg.norm(perp, axis=1)[:, None]
    return result


def bond_cylinders(coords: np.ndarray, bonds: BondLayout) -> np.ndarray:
    """Cylinder rows of a :class:`BondLayout` for a set of positions."""
    rows = bonds.rows.copy()
    if not len(rows):
        return rows
    offsets = np.zeros((len(rows), 3))
    moved = bonds.spread >= 0
    if moved.any():
        directions = _offset_directions(coords, bonds)
        offsets[moved] = directions[bonds.spread[moved]] * bonds.shift[moved][:, None]
    rows[:, 0:3] = coords[bonds.start] + offsets
    rows[:, 3:6] = coords[bonds.end] + offsets
    return rows


@dataclass
class Scene:
    """Geometry ready for drawing.

    Attributes
    ----------
    spheres : numpy.ndarray
        One row per drawn atom, shape ``(n, 8)``: centre, radius (negative:
        a width in pixels), colour and :class:`Highlight`.
    cylinders : numpy.ndarray
        One row per drawn bond line, shape ``(m, 14)``: both ends, radius,
        the colour of each half and :class:`CylinderKind`.
    sphere_atoms : numpy.ndarray
        Atom index of each sphere row.
    cylinder_trim : numpy.ndarray
        Radii of the two end atoms of each cylinder, shape ``(m, 2)``, so the
        QPainter fallback can hide the part inside the spheres.
    mesh : Mesh or None
        Cartoon surface.
    arrows : numpy.ndarray
        Arrows such as the dipole moment, shape ``(k, 10)``: tail, tip, shaft
        radius and colour (see :mod:`miles.render.arrows`).
    """

    spheres: np.ndarray = field(default_factory=lambda: np.zeros((0, 8), np.float32))
    cylinders: np.ndarray = field(default_factory=lambda: np.zeros((0, 14), np.float32))
    sphere_atoms: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    cylinder_trim: np.ndarray = field(default_factory=lambda: np.zeros((0, 2), np.float32))
    mesh: Mesh | None = None
    arrows: np.ndarray = field(default_factory=lambda: np.zeros((0, 10), np.float32))

    @property
    def empty(self) -> bool:
        """Whether there is nothing to draw."""
        return (len(self.spheres) == 0 and len(self.cylinders) == 0 and len(self.arrows) == 0
                and (self.mesh is None or self.mesh.empty))

    def flags(self, highlights: Mapping[int, Highlight], n_atoms: int) -> np.ndarray:
        """Highlight of each sphere row (the last column of :attr:`spheres`) for some atoms.

        Parameters
        ----------
        highlights : mapping of int to Highlight
            Highlight per atom index; atoms that are not drawn are skipped.
        n_atoms : int
            Number of atoms in the molecule.

        Returns
        -------
        numpy.ndarray
            One float per sphere row, 0 for :attr:`Highlight.NONE`.
        """
        flags = np.zeros(len(self.spheres), np.float32)
        if not highlights or not len(self.spheres):
            return flags
        atoms = np.fromiter(highlights.keys(), dtype=np.int64, count=len(highlights))
        values = np.fromiter(highlights.values(), dtype=np.float32, count=len(highlights))
        inside = (atoms >= 0) & (atoms < n_atoms)
        row = np.full(n_atoms, -1)
        row[self.sphere_atoms] = np.arange(len(self.sphere_atoms))
        rows = row[atoms[inside]]
        drawn = rows >= 0
        flags[rows[drawn]] = values[inside][drawn]
        return flags

    def highlighted(self, highlights: Mapping[int, Highlight], n_atoms: int) -> np.ndarray:
        """Copy of the sphere array with highlights set for some atoms.

        Parameters
        ----------
        highlights : mapping of int to Highlight
            Highlight per atom index; atoms that are not drawn are skipped.
        n_atoms : int
            Number of atoms in the molecule.

        Returns
        -------
        numpy.ndarray
            The sphere array (the original when there is nothing to change).
        """
        if not highlights or not len(self.spheres):
            return self.spheres
        spheres = self.spheres.copy()
        spheres[:, 7] = self.flags(highlights, n_atoms)
        return spheres


@dataclass
class SceneLayout:
    """What a scene takes from the atoms, bonds, display and theme, but not from the positions.

    Kept while atoms only move (trajectory frames, drags), so a new frame
    needs only the geometry worked out again.
    """

    reps: np.ndarray
    shown: np.ndarray
    radii: np.ndarray
    colors: np.ndarray
    cartoon: np.ndarray
    residues: list[Any] | None = None
    bonds: BondLayout | None = None


def scene_layout(mol: Molecule, style: Style, theme: Theme, colors: np.ndarray | None = None,
                 display: AtomDisplay | None = None) -> SceneLayout:
    """Work out the representations, visibility, radii, colours, cartoon residues and bond layout of a scene.

    See :func:`build_scene`, which takes the result as ``layout``.
    """
    reps = atom_reps(mol, style, display)
    shown = visible_atoms(mol, style, reps)
    radii = sphere_radii(mol.numbers, reps).astype(np.float32)
    if colors is not None and len(colors) == len(mol):
        atom_colors = np.asarray(colors, dtype=np.float32)
    else:
        atom_colors = theme.element_colors()[mol.numbers]
    cartoon = np.asarray(reps == Representation.CARTOON, dtype=bool)
    residues = cartoon_residues(mol, cartoon) if cartoon.any() else None
    bonds = bond_layout(mol, reps, shown, radii, atom_colors) if len(mol.bond_pairs) else None
    return SceneLayout(reps, shown, radii, atom_colors, cartoon, residues, bonds)


def build_scene(mol: Molecule, style: Style, theme: Theme, highlights: Mapping[int, Highlight] | None = None,
                with_mesh: bool = True, colors: np.ndarray | None = None,
                display: AtomDisplay | None = None, layout: SceneLayout | None = None) -> Scene:
    """Spheres, bond cylinders and cartoon mesh for a molecule.

    Parameters
    ----------
    mol : Molecule
        Structure to draw.
    style : Style
        Display options.
    theme : Theme
        Supplies the element colours.
    highlights : mapping of int to Highlight, optional
        Highlighted atoms.
    with_mesh : bool, optional
        Build the cartoon mesh when some atoms are drawn as cartoon.
    colors : numpy.ndarray, optional
        RGB colour of each atom, shape ``(n_atoms, 3)``, in place of the
        element colours (e.g. :func:`charge_colors`).
    display : AtomDisplay, optional
        Each atom's own representation and focus-mode hiding.
    layout : SceneLayout, optional
        :func:`scene_layout` of the same atoms, display and style when it is
        at hand; only the positions are then read from ``mol``.

    Returns
    -------
    Scene
        The scene.
    """
    n = len(mol)
    if n == 0:
        return Scene()
    coords = mol.coords
    if layout is None:
        layout = scene_layout(mol, style, theme, colors, display)
    shown, radii, colors = layout.shown, layout.radii, layout.colors

    idx = np.nonzero(shown)[0]
    spheres = np.zeros((len(idx), 8), np.float32)
    spheres[:, 0:3] = coords[idx]
    spheres[:, 3] = radii[idx]
    spheres[:, 4:7] = colors[idx]
    scene = Scene(spheres=spheres, sphere_atoms=idx)
    if highlights:
        scene.spheres = scene.highlighted(highlights, n)
    if with_mesh and layout.residues is not None:
        scene.mesh = cartoon_mesh(mol, layout.cartoon, style.cartoon_color, theme, layout.residues)
    if layout.bonds is not None and len(layout.bonds.rows):
        scene.cylinders = bond_cylinders(coords, layout.bonds)
        scene.cylinder_trim = layout.bonds.trims
    return scene


def trace_cylinders(mesh: Mesh | None) -> np.ndarray:
    """Cartoon backbone as short thick cylinders, for the QPainter fallback."""
    if mesh is None or len(mesh.centers) < 2:
        return np.zeros((0, 14), np.float32)
    c, col = mesh.centers, mesh.center_colors
    ok = np.linalg.norm(c[1:] - c[:-1], axis=1) < 3.0
    rows = np.zeros((int(ok.sum()), 14), np.float32)
    rows[:, 0:3] = c[:-1][ok]
    rows[:, 3:6] = c[1:][ok]
    rows[:, 6] = 0.45
    rows[:, 7:10] = col[:-1][ok]
    rows[:, 10:13] = col[1:][ok]
    rows[:, 13] = CylinderKind.TRACE
    return rows
