"""Commands of one canvas: editing, selecting, representations, focus mode, measurements and the clipboard.

The menus and the right-click menu call these on the current canvas
(:mod:`miles.ui.menus`); the tools call some of them too.
"""

from typing import TYPE_CHECKING

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from miles.core import building, perception, stereo
from miles.core.chem import ChemError, to_smiles
from miles.core.forcefield import CleanMode, clean_geometry, clean_selection
from miles.core.fragments import Category
from miles.core.molecule import Molecule, Representation
from miles.core.neighbors import cross_pairs
from miles.core.protein import Group, group_atoms
from miles.files.common import FileFormatError
from miles.render.camera import nice_rotation
from miles.ui.builder import Tool
from miles.ui.clipboard import copy_molecule, pasted_molecule
from miles.ui.document import Change
from miles.ui.text import plural

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

GROUP_NAMES = {Group.POLYMER: "protein", Group.LIGANDS: "ligands", Group.WATER: "water"}

PASTE_CLEARANCE = 1.2   # angstrom; a pasted piece closer than this to the structure is moved aside

BACKGROUND_ATOMS = 2000   # Add Hydrogens and Clean Geometry on larger structures run on a worker thread

CHANGED_MEANWHILE = "The structure changed while this ran; try again"

SMILES_SHOWN = 60       # characters of a copied SMILES string the message shows


class Commands:
    """What the menus do to the document and view of one canvas.

    Parameters
    ----------
    canvas : Canvas
        The canvas the commands act on.
    """

    def __init__(self, canvas: "Canvas") -> None:
        self.canvas = canvas

    # -- history ------------------------------------------------------------------------------

    def undo(self) -> None:
        """Undo the last edit."""
        self.canvas.doc.undo()

    def redo(self) -> None:
        """Redo the last undone edit."""
        self.canvas.doc.redo()

    # -- selecting ----------------------------------------------------------------------------

    def select_all(self) -> None:
        """Select every atom."""
        canvas = self.canvas
        if len(canvas.mol):
            canvas.builder.set_tool(Tool.SELECT)
            canvas.set_selection(list(range(len(canvas.mol))))

    def select_group(self, group: Group, around: int | None = None) -> None:
        """Select a residue, chain, the polymer, ligands, water or everything."""
        canvas = self.canvas
        atoms = group_atoms(canvas.mol, group, around)
        canvas.builder.set_tool(Tool.SELECT)
        canvas.set_selection(atoms)
        if not atoms:
            canvas.toast(f"No {GROUP_NAMES.get(group, group)} in this structure")

    # -- what is shown ------------------------------------------------------------------------

    def set_representation(self, rep: Representation, atoms: list[int] | None = None) -> None:
        """Draw some atoms (default: the selection, else everything) in a representation."""
        chosen = self.canvas.selected_or_all() if atoms is None else atoms
        self.canvas.doc.set_representation(rep, chosen)

    def set_display_mode(self, mode: Representation) -> None:
        """Draw the selection in a representation, or everything (clearing per-atom choices)."""
        canvas = self.canvas
        if canvas.selection:
            self.set_representation(mode)
            return
        if (canvas.doc.display.rep != Representation.DEFAULT).any():
            self.set_representation(Representation.DEFAULT, list(range(len(canvas.mol))))
        canvas.set_style(mode=mode)

    def hidden_count(self) -> int:
        """Count the atoms hidden by focus mode."""
        return int(self.canvas.doc.display.hidden.sum())

    def _set_hidden(self, hidden: set[int]) -> None:
        mask = np.zeros(len(self.canvas.mol), dtype=bool)
        mask[list(hidden)] = True
        self.canvas.doc.set_hidden(mask)

    def show_only_selection(self) -> None:
        """Hide everything except the selection and zoom to it."""
        canvas = self.canvas
        if not canvas.selection:
            canvas.toast("Select the atoms to focus on first")
            return
        keep = set(canvas.selection)
        self._set_hidden({i for i in range(len(canvas.mol)) if i not in keep})
        canvas.view.fit(animate=True)

    def hide_selection(self) -> None:
        """Hide the selected atoms."""
        canvas = self.canvas
        if not canvas.selection:
            canvas.toast("Select the atoms to hide first")
            return
        hidden = set(np.flatnonzero(canvas.doc.display.hidden).tolist()) | set(canvas.selection)
        canvas.set_selection([])
        self._set_hidden(hidden)

    def show_everything(self) -> None:
        """Show the atoms hidden by focus mode again."""
        if self.hidden_count():
            self._set_hidden(set())
            self.canvas.view.fit(animate=True)

    # -- measurements -------------------------------------------------------------------------

    def toggle_measurement(self, atoms: list[int]) -> bool:
        """Pin or unpin the distance, angle or dihedral of 2-4 atoms; return whether it is now pinned."""
        key = tuple(atoms)
        if not 2 <= len(key) <= 4:
            return False
        pinned = key in self.canvas.mol.measurements or key[::-1] in self.canvas.mol.measurements

        def apply(mol: Molecule) -> None:
            if pinned:
                mol.measurements = [m for m in mol.measurements if m not in (key, key[::-1])]
            else:
                mol.measurements.append(key)

        self.canvas.doc.edit(apply, kind=Change.STYLE)
        return not pinned

    def clear_measurements(self) -> None:
        """Unpin every measurement."""
        if self.canvas.mol.measurements:
            self.canvas.doc.edit(lambda mol: setattr(mol, "measurements", []), kind=Change.STYLE)

    # -- editing ------------------------------------------------------------------------------

    def invert_center(self, atom: int) -> None:
        """Invert a stereocentre (R to S)."""
        if not self.canvas.doc.edit(stereo.invert_center, atom, kind=Change.COORDS):
            self.canvas.toast("This centre cannot be inverted without breaking a ring")

    def flip_double_bond(self, i: int, j: int) -> None:
        """Swap E and Z at a double bond."""
        if not self.canvas.doc.edit(stereo.flip_double_bond, i, j, kind=Change.COORDS):
            self.canvas.toast("Double bonds in rings cannot be flipped")

    def delete_atoms(self, atoms: list[int]) -> None:
        """Delete atoms with their hydrogens, clearing the selection."""
        canvas = self.canvas
        canvas.selection = []
        canvas.doc.edit(building.delete_atoms, atoms)
        canvas.selectionChanged.emit()

    def delete_selection(self) -> None:
        """Delete the selected atoms."""
        self.delete_atoms(self.canvas.selection)

    def _in_background(self) -> bool:
        """Whether a slow command runs on a worker thread (large structures, in a window)."""
        return self.canvas.tasks is not None and len(self.canvas.mol) >= BACKGROUND_ATOMS

    def add_hydrogens(self) -> None:
        """Add hydrogens to the selection (keeping them selected), or to everything.

        On large structures the hydrogens are worked out on a copy on a
        worker thread and added as one undoable step when it is done.
        """
        canvas = self.canvas
        doc = canvas.doc
        selection = list(canvas.selection)
        if not self._in_background():
            start = len(canvas.mol)
            n = doc.edit(building.add_hydrogens, selection or None)
            self._hydrogens_added(selection, start, n)
            return
        revision = doc.revision

        def done(result: tuple[Molecule, int]) -> None:
            mol, n = result
            if doc.revision != revision:
                canvas.toast(CHANGED_MEANWHILE)
                return
            start = len(doc.mol)
            doc.checkpoint()
            doc.mol = mol
            doc.notify(Change.STRUCTURE)
            self._hydrogens_added(selection, start, n)

        assert canvas.tasks is not None
        canvas.tasks.start("Adding hydrogens", _add_hydrogens, doc.mol.copy(), selection or None, done=done)

    def _hydrogens_added(self, selection: list[int], start: int, n: int) -> None:
        canvas = self.canvas
        if selection:
            canvas.set_selection(selection + list(range(start, len(canvas.mol))))
        canvas.toast(f"Added {plural(n, 'hydrogen')}{' to the selection' if selection else ''}")

    def clean_geometry(self) -> None:
        """Relax the selection with its surroundings fixed, or the whole structure.

        Experimental structures are only tidied, not rebuilt. On large
        structures the work runs on a copy on a worker thread and the new
        coordinates are applied as one undoable step when it is done.
        """
        canvas = self.canvas
        doc = canvas.doc
        if len(canvas.mol) < 2:
            return
        atoms = list(canvas.selection)
        if not self._in_background():
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            try:
                mode = doc.edit(clean_selection, atoms, kind=Change.COORDS) if atoms else \
                    doc.edit(clean_geometry, kind=Change.COORDS)
            finally:
                QApplication.restoreOverrideCursor()
            self._cleaned(atoms, mode)
            return
        revision = doc.revision

        def done(result: tuple[np.ndarray, CleanMode | None]) -> None:
            coords, mode = result
            if doc.revision != revision:
                canvas.toast(CHANGED_MEANWHILE)
                return
            doc.edit(_set_coordinates, coords, kind=Change.COORDS)
            self._cleaned(atoms, mode)

        assert canvas.tasks is not None
        canvas.tasks.start("Cleaning the geometry", _clean, doc.mol.copy(), atoms, done=done)

    def _cleaned(self, atoms: list[int], mode: CleanMode | None) -> None:
        canvas = self.canvas
        if atoms:
            canvas.toast(f"Cleaned {plural(len(atoms), 'selected atom')}")
        elif mode == CleanMode.GENTLE:
            canvas.toast("Tidied, keeping the experimental structure; select atoms to rebuild them", 2600)
        else:
            canvas.toast("Cleaned with MMFF94" if mode == CleanMode.MMFF else "Cleaned")

    def rebond(self) -> None:
        """Recompute bonds from the coordinates."""
        self.canvas.doc.edit(perception.rebond)

    def move_to_origin(self) -> None:
        """Move the centre of mass to the origin and centre the view."""
        canvas = self.canvas
        if not len(canvas.mol):
            return

        def apply(mol: Molecule) -> None:
            mol.coords = mol.coords - mol.center_of_mass()

        canvas.doc.edit(apply, kind=Change.COORDS)
        canvas.view.center()

    # -- fragments ----------------------------------------------------------------------------

    def placement_rotation(self) -> np.ndarray:
        """Orientation for a new fragment: as in the preview, relative to the view."""
        builder = self.canvas.builder
        frag = builder.fragment
        R = builder.preview_rotation
        if R is None:
            R = nice_rotation(frag.mol.coords, tilt=frag.category != Category.ELEMENT)
        return self.canvas.camera.rotation.T @ R

    def place_current_fragment(self) -> None:
        """Place the current fragment in an empty document, or explain how to place it."""
        canvas = self.canvas
        frag = canvas.builder.fragment
        if len(canvas.mol):
            canvas.toast(f"Click to place {frag.name.lower()}, or click a hydrogen to attach it")
            return
        canvas.camera.pan[:] = 0
        canvas.doc.edit(building.place_fragment, frag.mol, canvas.camera.center.copy(),
                        self.placement_rotation())
        canvas.view.fit(animate=True)

    # -- clipboard ----------------------------------------------------------------------------

    def copy_atoms(self, cut: bool = False) -> int:
        """Copy (or cut) the selection, or everything, to the clipboard; returns the atom count."""
        atoms = self.canvas.selected_or_all()
        if atoms:
            copy_molecule(self.canvas.mol.subset(atoms))
            if cut:
                self.delete_atoms(atoms)
        return len(atoms)

    def copy(self, cut: bool = False) -> None:
        """Copy (or cut) the selection or everything, and say so."""
        n = self.copy_atoms(cut=cut)
        if n:
            self.canvas.toast(f"{'Cut' if cut else 'Copied'} {plural(n, 'atom')}")

    def paste(self) -> None:
        """Insert the clipboard contents next to the structure."""
        canvas = self.canvas
        try:
            piece = pasted_molecule()
        except (FileFormatError, ValueError) as exc:
            canvas.toast(str(exc).split("\n")[0][:80], 2600)
            return
        if piece is None or not len(piece):
            canvas.toast("Nothing to paste")
            return
        self.paste_molecule(piece)
        canvas.toast(f"Pasted {len(piece)} atoms")

    def duplicate(self) -> None:
        """Paste a copy of the selection or everything."""
        canvas = self.canvas
        if len(canvas.mol):
            self.paste_molecule(canvas.mol.subset(canvas.selected_or_all()))

    def copy_smiles(self) -> None:
        """Put the SMILES string of the selection or everything on the clipboard."""
        canvas = self.canvas
        if not len(canvas.mol):
            return
        try:
            smiles = to_smiles(canvas.mol.subset(canvas.selected_with_hydrogens()))
        except ChemError as exc:
            canvas.toast(str(exc))
            return
        QApplication.clipboard().setText(smiles)
        shown = smiles if len(smiles) < SMILES_SHOWN else smiles[:SMILES_SHOWN - 3] + "..."
        canvas.toast(f"Copied {shown}", 2400)

    def paste_molecule(self, piece: Molecule) -> None:
        """Insert a molecule next to the structure (not on top of it) and select it."""
        canvas = self.canvas
        if not len(piece):
            return
        piece = piece.copy()
        mol = canvas.mol
        if not len(mol):
            piece.coords = piece.coords + (canvas.camera.center - piece.coords.mean(axis=0))
        elif len(cross_pairs(piece.coords, mol.coords, PASTE_CLEARANCE)[0]):   # it would overlap
            right = canvas.camera.rotation[0]
            gap = (mol.coords @ right).max() - (piece.coords @ right).min() + 2.5
            piece.coords = piece.coords + right * gap
        start = len(mol)
        canvas.doc.edit(Molecule.extend, piece)
        canvas.builder.set_tool(Tool.SELECT)
        canvas.set_selection(list(range(start, start + len(piece))))
        canvas.view.ensure_visible()


def _add_hydrogens(mol: Molecule, atoms: list[int] | None) -> tuple[Molecule, int]:
    """Add hydrogens to a copy of a structure (on a worker thread); returns it and the number added."""
    return mol, building.add_hydrogens(mol, atoms)


def _clean(mol: Molecule, atoms: list[int]) -> tuple[np.ndarray, CleanMode | None]:
    """Clean a copy of a structure (on a worker thread); returns its coordinates and how it was cleaned."""
    if atoms:
        clean_selection(mol, atoms)
        return mol.coords.copy(), None
    mode = clean_geometry(mol)
    return mol.coords.copy(), mode


def _set_coordinates(mol: Molecule, coords: np.ndarray) -> None:
    mol.coords = coords
