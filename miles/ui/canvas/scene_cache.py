"""Scene arrays for a canvas, rebuilt only when their inputs change.

Building a scene walks every atom and bond, so it is done once per edit or
style change. Highlights and the cartoon tint change far more often (every
mouse move): the highlights are a small array of their own, one value per
sphere, so a hover change sends a few bytes to the GPU instead of every
sphere again.
"""

import hashlib
from collections.abc import Hashable, Mapping
from dataclasses import dataclass, field

import numpy as np

from miles.core.display import AtomDisplay
from miles.core.molecule import Molecule
from miles.render.arrows import COLUMNS, arrow_mesh
from miles.render.camera import Camera, DepthRange, depth_range
from miles.render.cartoon import Mesh, tint_selection
from miles.render.layers import MeshItem
from miles.render.picking import Picker
from miles.render.scene import Highlight, Scene, SceneLayout, Style, build_scene, scene_layout
from miles.render.theme import Theme

DEPTH_MARGIN = 3.0


@dataclass
class Frame:
    """Arrays to draw now, and which of them changed since the previous frame.

    :attr:`spheres` are the scene's spheres without highlights; the
    highlight of each sphere is in :attr:`flags`, so a renderer can send
    only the highlights that changed.
    """

    spheres: np.ndarray
    cylinders: np.ndarray
    mesh: Mesh | None
    mesh_colors: np.ndarray | None
    spheres_changed: bool
    cylinders_changed: bool
    mesh_changed: bool
    arrows: np.ndarray = field(default_factory=lambda: np.zeros((0, COLUMNS), np.float32))
    arrow_mesh: Mesh | None = None
    arrows_changed: bool = False
    meshes: list[MeshItem] = field(default_factory=list)
    flags: np.ndarray = field(default_factory=lambda: np.zeros(0, np.float32))

    def layers(self) -> list[tuple[Mesh | None, float]]:
        """Return the plugins' meshes with their opacity (for the QPainter renderer)."""
        return [(item.mesh, item.alpha) for item in self.meshes]

    def highlighted_spheres(self) -> np.ndarray:
        """Return the spheres with their highlights in the last column."""
        if not self.flags.any():
            return self.spheres
        spheres = self.spheres.copy()
        spheres[:, 7] = self.flags
        return spheres

    def scene(self, base: Scene) -> Scene:
        """Return the frame as a :class:`Scene` (for the QPainter renderer)."""
        return Scene(spheres=self.highlighted_spheres(), cylinders=self.cylinders,
                     sphere_atoms=base.sphere_atoms, cylinder_trim=base.cylinder_trim, mesh=self.mesh,
                     arrows=self.arrows)


class SceneCache:
    """Cached scene, picker and bounds of one canvas."""

    def __init__(self) -> None:
        self.base = Scene()
        self.center = np.zeros(3)
        self.radius = 1.0
        self._base_key: Hashable = None
        self._rows = np.zeros(0, int)   # sphere row of each atom, -1 for atoms not drawn as spheres
        self._mesh_key: np.ndarray | None = None
        self._mesh_colors: np.ndarray | None = None
        self._picker: Picker | None = None
        self._picker_key: Hashable = None
        self._arrows = np.zeros((0, COLUMNS), np.float32)
        self._arrow_mesh: Mesh | None = None
        self._colors: tuple[np.ndarray, str] | None = None   # last atom colours and their digest
        self._layout: SceneLayout | None = None
        self._layout_key: Hashable = None

    def invalidate(self) -> None:
        """Force a rebuild on the next frame (e.g. after a theme change)."""
        self._base_key = None
        self._layout_key = None

    def frame(self, mol: Molecule, revision: int, style: Style, theme: Theme,
              highlights: Mapping[int, Highlight], selection: list[int],
              atom_colors: np.ndarray | None = None, arrows: np.ndarray | None = None,
              meshes: list[MeshItem] | None = None, display: AtomDisplay | None = None,
              layout_revision: int | None = None) -> Frame:
        """Arrays for the next frame.

        Parameters
        ----------
        mol : Molecule
            Structure on screen.
        revision : int
            Document revision; a new value means the structure changed.
        style : Style
            Display options.
        theme : Theme
            Colours.
        highlights : mapping of int to Highlight
            Highlights besides the selection, such as the atom under the
            cursor; a selected atom keeps :attr:`Highlight.SELECTED`.
        selection : list of int
            Selected atoms, highlighted and tinted in the cartoon.
        atom_colors : numpy.ndarray, optional
            RGB colour of each atom in place of the element colours.
        arrows : numpy.ndarray, optional
            Arrow rows (:mod:`miles.render.arrows`).
        meshes : list of MeshItem, optional
            Meshes of plugins, such as surfaces and isosurfaces.
        display : AtomDisplay, optional
            Each atom's own representation and focus-mode hiding; a change
            of it comes with a new revision.
        layout_revision : int, optional
            Document revision that changes with everything but positions;
            while it stays, a new ``revision`` only moves the atoms, and the
            representations, colours and cartoon residues are kept.

        Returns
        -------
        Frame
            The arrays and change flags.
        """
        layout_key = (id(mol), revision if layout_revision is None else layout_revision, style.mode,
                      style.show_hydrogens, theme.name, style.cartoon_color, self._colors_key(atom_colors))
        key = (*layout_key, revision)
        base_changed = key != self._base_key
        if base_changed:
            self._base_key = key
            if layout_key != self._layout_key or self._layout is None:
                self._layout_key = layout_key
                self._layout = scene_layout(mol, style, theme, atom_colors, display)
            self.base = build_scene(mol, style, theme, layout=self._layout)
            self._rows = np.full(len(mol), -1)
            self._rows[self.base.sphere_atoms] = np.arange(len(self.base.sphere_atoms))
            if len(mol):
                c = mol.coords
                self.center = (c.min(axis=0) + c.max(axis=0)) / 2
                self.radius = float(np.linalg.norm(c - self.center, axis=1).max()) + DEPTH_MARGIN

        selected = np.asarray(selection, dtype=np.int64)
        flags = self._flags(selected, highlights)

        mesh = self.base.mesh
        mesh_key = np.unique(selected) if mesh is not None else None
        mesh_changed = base_changed or not _same(mesh_key, self._mesh_key)
        if mesh_changed:
            self._mesh_key = mesh_key
            self._mesh_colors = None
            if mesh is not None and not mesh.empty:
                self._mesh_colors = tint_selection(mesh, set(selection), theme.accent)

        arrows_changed = self._update_arrows(arrows)
        return Frame(self.base.spheres, self.base.cylinders, mesh, self._mesh_colors,
                     base_changed, base_changed, mesh_changed,
                     self._arrows, self._arrow_mesh, arrows_changed, list(meshes or []), flags)

    def _flags(self, selected: np.ndarray, highlights: Mapping[int, Highlight]) -> np.ndarray:
        """Highlight of each sphere row: the selection, then other highlights on atoms not selected."""
        flags = np.zeros(len(self.base.spheres), np.float32)
        if not len(flags):
            return flags
        flags[self._sphere_rows(selected)[0]] = Highlight.SELECTED
        if highlights:
            atoms = np.fromiter(highlights.keys(), np.int64, len(highlights))
            values = np.fromiter(highlights.values(), np.float32, len(highlights))
            rows, drawn = self._sphere_rows(atoms)
            values = values[drawn]
            free = flags[rows] == Highlight.NONE   # a selected atom stays selected
            flags[rows[free]] = values[free]
        return flags

    def _sphere_rows(self, atoms: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Sphere rows of the atoms drawn as spheres, and which of the atoms those are."""
        rows = np.full(len(atoms), -1)
        inside = (atoms >= 0) & (atoms < len(self._rows))
        rows[inside] = self._rows[atoms[inside]]
        drawn = rows >= 0
        return rows[drawn], drawn

    def _colors_key(self, colors: np.ndarray | None) -> str | None:
        """Digest of the atom colours, worked out again only when another array is passed."""
        if colors is None:
            return None
        if self._colors is None or self._colors[0] is not colors:
            data = np.ascontiguousarray(colors, dtype=np.float32)
            self._colors = (colors, hashlib.blake2b(data.tobytes(), digest_size=16).hexdigest())
        return self._colors[1]

    def _update_arrows(self, arrows: np.ndarray | None) -> bool:
        """Rebuild the arrow mesh if the arrows changed; return whether they did."""
        rows = np.zeros((0, COLUMNS), np.float32) if arrows is None else np.asarray(arrows, np.float32)
        if rows.shape == self._arrows.shape and np.array_equal(rows, self._arrows):
            return False
        self._arrows = rows
        self._arrow_mesh = arrow_mesh(rows)
        return True

    def picker(self, mol: Molecule, revision: int, style: Style,
               display: AtomDisplay | None = None) -> Picker:
        """Hit tester for the current structure, display and style.

        It reads positions at each call, so ``revision`` may be a revision
        that does not change when atoms only move.
        """
        key = (id(mol), revision, style.mode, style.show_hydrogens)
        if self._picker is None or key != self._picker_key:
            self._picker = Picker(mol, style, display)
            self._picker_key = key
        return self._picker

    def depth(self, camera: Camera) -> DepthRange:
        """Depth range of the structure in the current view."""
        return depth_range(camera, self.center, self.radius)


def _same(a: np.ndarray | None, b: np.ndarray | None) -> bool:
    """Whether two optional arrays hold the same values."""
    if a is None or b is None:
        return a is b
    return a.shape == b.shape and bool(np.array_equal(a, b))
