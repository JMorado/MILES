"""Finding the atom or bond under the cursor."""

from dataclasses import dataclass

import numpy as np

from miles.core.display import AtomDisplay
from miles.core.molecule import Molecule, Representation
from miles.render.camera import Camera
from miles.render.scene import BOND_RADIUS, TUBE_RADIUS, Style, atom_reps, sphere_radii, visible_atoms

WIRE_PICK_PX = 6.0
CARTOON_PICK_RADIUS = 1.1
MIN_PICK_PX = 5.0


@dataclass(frozen=True, slots=True)
class AtomHit:
    """The cursor is over atom ``index``."""

    index: int


@dataclass(frozen=True, slots=True)
class BondHit:
    """The cursor is over the bond between atoms ``i`` and ``j``."""

    i: int
    j: int


Hit = AtomHit | BondHit


def first_atom(hit: Hit) -> int:
    """Return the atom of an atom hit, or the first atom of a bond hit."""
    return hit.index if isinstance(hit, AtomHit) else hit.i


class Picker:
    """Hit testing for one molecule and style.

    The per-atom radii and bond list are computed once; build a new picker
    when the structure or the style changes. Coordinates are read at each
    call, so moving atoms needs no rebuild.

    Parameters
    ----------
    mol : Molecule
        Structure on screen.
    style : Style
        Display options (they decide sizes and which atoms are drawn).
    display : AtomDisplay, optional
        Each atom's own representation and focus-mode hiding.
    """

    def __init__(self, mol: Molecule, style: Style, display: AtomDisplay | None = None) -> None:
        self.mol = mol
        self.style = style
        reps = atom_reps(mol, style, display)
        self.radii, self.pickable = self._pick_radii(reps)
        self.drawn = visible_atoms(mol, style, reps)
        self.pairs = np.array(mol.bond_pairs, dtype=int).reshape(-1, 2)

    def _pick_radii(self, reps: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Pick radius per atom (angstrom; negative: pixels) and which atoms can be picked.

        In a cartoon only the CA (or P) atom of each residue can be picked.
        """
        mol = self.mol
        numbers = mol.numbers
        cartoon = reps == Representation.CARTOON
        guides = cartoon & mol.mask("name", {"CA", "P"})
        ok = (reps != Representation.HIDDEN) & ~cartoon
        if not self.style.show_hydrogens:
            hydrogen = numbers == 1
            ok &= ~hydrogen
            guides &= ~hydrogen
        radii = np.where(reps == Representation.WIRE, -WIRE_PICK_PX, sphere_radii(numbers, reps))
        radii = np.where(guides, CARTOON_PICK_RADIUS, np.where(ok, radii, 0.0)).astype(float)
        return radii, np.asarray(ok | guides, dtype=bool)

    def pick(self, camera: Camera, width: float, height: float, x: float, y: float) -> Hit | None:
        """Front-most atom or bond at a screen position.

        Parameters
        ----------
        camera : Camera
            Current view.
        width, height : float
            Size of the view in pixels.
        x, y : float
            Cursor position in pixels.

        Returns
        -------
        AtomHit, BondHit or None
            What is under the cursor; atoms win over bonds.
        """
        if not len(self.mol):
            return None
        view = camera.to_view(self.mol.coords)
        sx, sy = camera.to_screen(view, width, height)
        r = np.where(self.radii < 0, -self.radii, np.maximum(MIN_PICK_PX, self.radii * camera.scale))
        hits = np.nonzero((np.hypot(sx - x, sy - y) <= r + 2) & self.pickable)[0]
        if len(hits):
            return AtomHit(int(hits[np.argmax(view[hits, 2])]))
        if not len(self.pairs):
            return None
        i0, i1 = self.pairs[:, 0], self.pairs[:, 1]
        ax, ay, bx, by = sx[i0], sy[i0], sx[i1], sy[i1]
        dx, dy = bx - ax, by - ay
        t = np.clip(((x - ax) * dx + (y - ay) * dy) / np.maximum(dx * dx + dy * dy, 1e-9), 0.0, 1.0)
        distance = np.hypot(ax + t * dx - x, ay + t * dy - y)
        bond_radius = TUBE_RADIUS if self.style.mode == Representation.TUBE else BOND_RADIUS
        tolerance = max(MIN_PICK_PX, bond_radius * camera.scale + 2)
        near = np.nonzero(self.drawn[i0] & self.drawn[i1] & (distance <= tolerance))[0]
        if not len(near):
            return None
        z = view[i0[near], 2] + t[near] * (view[i1[near], 2] - view[i0[near], 2])
        k = near[np.argmax(z)]
        return BondHit(int(i0[k]), int(i1[k]))
