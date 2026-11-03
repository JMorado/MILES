"""What clicks and drags on the canvas do with each tool.

Each tool is a small stateless object; the canvas asks the active one how to
highlight the atom under the cursor, whether a press starts a tool-specific
drag, and what a click or double click does.
"""

from dataclasses import dataclass
from enum import Enum, auto
from typing import TYPE_CHECKING

import numpy as np
from PySide6.QtCore import QPointF, Qt

from miles.core import building
from miles.core.forcefield import clean_atoms
from miles.core.fragments import Category
from miles.core.molecule import Molecule
from miles.render.picking import AtomHit, BondHit, Hit, first_atom
from miles.render.scene import Highlight
from miles.ui.builder import Tool
from miles.ui.text import atom_label

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

NEXT_ORDER = {1.0: 2.0, 2.0: 3.0, 3.0: 1.0, 1.5: 1.0}
"""Clicking a bond with the build tool cycles single, double, triple."""

ORDER_TOASTS = {1.0: "Single bond", 2.0: "Double bond", 3.0: "Triple bond"}


class Drag(Enum):
    """What dragging the mouse does."""

    ROTATE = auto()       # turn the view
    SPIN = auto()         # turn the view about the line of sight
    PAN = auto()
    ZOOM = auto()         # right button: up zooms in
    SCALE = auto()        # scale mouse mode: right or up zooms in
    DRAW = auto()         # pull a new bond out of an atom
    MOVE_ATOMS = auto()
    TURN_ATOMS = auto()


@dataclass
class Gesture:
    """A mouse press and the drag that may follow it."""

    press: QPointF
    last: QPointF
    button: Qt.MouseButton
    drag: Drag | None = None
    moved: bool = False
    atoms: list[int] | None = None        # atoms moved by MOVE_ATOMS and TURN_ATOMS
    checkpointed: bool = False            # undo step recorded for the move
    draw_from: int | None = None          # DRAW: start atom
    draw_to: int | None = None            # DRAW: atom under the cursor, if any


class CanvasTool:
    """Behaviour of a tool; the defaults only inspect."""

    def highlight(self, canvas: "Canvas", hit: Hit) -> Highlight:
        """Highlight for the atoms under the cursor."""
        return Highlight.HOVER

    def press(self, canvas: "Canvas", gesture: Gesture, hit: Hit | None) -> bool:
        """Start a tool-specific drag from a plain left press; False rotates the view instead."""
        return False

    def click(self, canvas: "Canvas", hit: Hit | None, pos: QPointF, additive: bool) -> None:
        """Handle a click that did not move."""

    def double_click(self, canvas: "Canvas", hit: Hit | None) -> bool:
        """Handle a double click; return True to treat it as a new press instead."""
        if hit is None and len(canvas.mol):
            canvas.view.fit(animate=True)
        return False


class SelectTool(CanvasTool):
    """Click selects; the order of the selected atoms defines bonds, angles and dihedrals."""

    def press(self, canvas: "Canvas", gesture: Gesture, hit: Hit | None) -> bool:
        """Dragging a selected atom moves the whole selection."""
        if isinstance(hit, AtomHit) and hit.index in canvas.selection:
            gesture.drag = Drag.MOVE_ATOMS
            gesture.atoms = list(canvas.selection)
            return True
        return False

    def click(self, canvas: "Canvas", hit: Hit | None, pos: QPointF, additive: bool) -> None:
        """Select an atom or a bond's two atoms; with Shift, add or remove them."""
        selection = list(canvas.selection)
        if isinstance(hit, BondHit):
            pair = [hit.i, hit.j]
            selection = selection + [a for a in pair if a not in selection] if additive else pair
        elif isinstance(hit, AtomHit):
            if not additive:
                selection = [hit.index]
            elif hit.index in selection:
                selection.remove(hit.index)
            else:
                selection.append(hit.index)
        elif not additive:
            selection = []
        if selection != canvas.selection:
            canvas.set_selection(selection)

    def double_click(self, canvas: "Canvas", hit: Hit | None) -> bool:
        """Select the residue, then the chain, then everything."""
        if hit is None:
            return super().double_click(canvas, hit)
        mol = canvas.mol
        start = first_atom(hit)
        current = set(canvas.selection)
        residue, chain = set(mol.residue_atoms(start)), set(mol.chain_atoms(start))
        if current == residue and chain != residue:
            canvas.set_selection(sorted(chain))
        elif current in (residue, chain) and len(current) < len(mol):
            canvas.set_selection(list(range(len(mol))))
        else:
            canvas.set_selection(sorted(residue))
        return False


class BuildTool(CanvasTool):
    """Places, attaches and grows the current fragment."""

    def highlight(self, canvas: "Canvas", hit: Hit) -> Highlight:
        """Hydrogens show where the fragment will attach."""
        if isinstance(hit, AtomHit) and canvas.mol.atoms[hit.index].symbol == "H":
            return Highlight.ATTACH
        return Highlight.HOVER

    def press(self, canvas: "Canvas", gesture: Gesture, hit: Hit | None) -> bool:
        """Dragging from an atom grows the fragment from it."""
        if isinstance(hit, AtomHit):
            gesture.drag = Drag.DRAW
            gesture.draw_from = hit.index
            return True
        return False

    def click(self, canvas: "Canvas", hit: Hit | None, pos: QPointF, additive: bool) -> None:
        """Place on empty space, attach to a hydrogen, change an element or cycle a bond order."""
        frag = canvas.builder.fragment
        mol = canvas.mol
        canvas.placed_last = hit is None
        if hit is None:
            if len(mol) == 0:
                canvas.camera.pan[:] = 0
                position = canvas.camera.center.copy()
            else:
                position = canvas.camera.unproject(pos.x(), pos.y(), canvas.width(), canvas.height())
            origin = frag.mol.coords[frag.anchor] if frag.category == Category.ELEMENT else None
            canvas.doc.edit(building.place_fragment, frag.mol, position, canvas.commands.placement_rotation(),
                            origin)
        elif isinstance(hit, BondHit):
            order = NEXT_ORDER.get(mol.bond_order(hit.i, hit.j), 1.0)
            canvas.doc.edit(lambda m: m.add_bond(hit.i, hit.j, order))
            canvas.toast(ORDER_TOASTS[order], 900)
            return
        elif _can_attach(mol, hit.index, frag.category):
            canvas.doc.edit(building.attach_fragment, hit.index, frag.mol, frag.attach)
        elif frag.kind is not None:
            canvas.doc.edit(building.mutate_atom, hit.index, frag.symbol, frag.kind)
        else:
            canvas.toast("Click a hydrogen to attach the fragment")
            return
        canvas.view.ensure_visible()

    def double_click(self, canvas: "Canvas", hit: Hit | None) -> bool:
        """Build again on a quick second click, unless the first placed a free fragment."""
        return not canvas.placed_last


def _can_attach(mol: Molecule, i: int, category: str) -> bool:
    """Whether a fragment attaches at atom ``i``: hydrogens, or terminal atoms for non-element fragments."""
    return mol.atoms[i].symbol == "H" or (category != Category.ELEMENT and mol.degree(i) <= 1)


class DeleteTool(CanvasTool):
    """Deletes atoms (with their hydrogens) and bonds."""

    def highlight(self, canvas: "Canvas", hit: Hit) -> Highlight:
        """Show what will be deleted."""
        return Highlight.DANGER

    def click(self, canvas: "Canvas", hit: Hit | None, pos: QPointF, additive: bool) -> None:
        """Delete the atom or bond under the cursor."""
        if isinstance(hit, AtomHit):
            canvas.doc.edit(building.delete_atom, hit.index)
        elif isinstance(hit, BondHit):
            canvas.doc.edit(lambda m: m.remove_bond(hit.i, hit.j))


TOOLS: dict[Tool, CanvasTool] = {
    Tool.SELECT: SelectTool(),
    Tool.BUILD: BuildTool(),
    Tool.DELETE: DeleteTool(),
}


def finish_draw(canvas: "Canvas", gesture: Gesture, pos: QPointF) -> None:
    """End a build drag: bond to the atom under the cursor, or grow towards the cursor."""
    start = gesture.draw_from
    mol = canvas.mol
    if start is None or start >= len(mol):
        return
    if gesture.draw_to is not None:
        canvas.doc.edit(_join_and_relax, start, gesture.draw_to)
        return
    w, h = canvas.width(), canvas.height()

    def project(point: np.ndarray) -> np.ndarray:
        return canvas.camera.project(point, w, h)

    direction = np.array([pos.x(), pos.y()]) - project(mol.coords[start])
    target = building.best_hydrogen_towards(mol, start, direction, project)
    if target is None:
        canvas.toast(f"No free position on {atom_label(mol, start)}")
        return
    frag = canvas.builder.fragment
    canvas.doc.edit(building.attach_fragment, target, frag.mol, frag.attach)
    canvas.view.ensure_visible()


def _join_and_relax(mol: Molecule, i: int, j: int) -> tuple[int, int] | None:
    """Bond two atoms and relax the molecule they now share (closing a ring folds it)."""
    pair = building.join_atoms(mol, i, j)
    if pair is not None:
        clean_atoms(mol, mol.component(pair[0]), fold=True)
    return pair
