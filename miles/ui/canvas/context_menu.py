"""Right-click menu of the canvas."""

from collections.abc import Callable
from functools import partial
from typing import TYPE_CHECKING

from PySide6.QtCore import QPointF
from PySide6.QtGui import QAction
from PySide6.QtWidgets import QMenu

from miles.core import building
from miles.core.elements import element
from miles.core.fragments import default_element_fragment
from miles.core.molecule import Molecule
from miles.core.protein import Group
from miles.render.picking import AtomHit, BondHit, Hit, first_atom
from miles.render.scene import REPRESENTATION_NAMES
from miles.ui.canvas.extensions import MenuItem, MenuTarget
from miles.ui.text import atom_label

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

COMMON_ELEMENTS = ["H", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I"]
BOND_ORDERS = [("Single", 1.0), ("Double", 2.0), ("Triple", 3.0), ("Aromatic", 1.5)]
GROUP_ITEMS = [(Group.POLYMER, "Protein and Nucleic Acids"), (Group.LIGANDS, "Ligands"),
               (Group.WATER, "Water"), (Group.ALL, "Everything")]


def _add(menu: QMenu, text: str, slot: Callable[[], object]) -> QAction:
    """Add an item that calls ``slot()``; Qt's checked argument is not passed on."""
    action = menu.addAction(text)
    action.triggered.connect(lambda: slot())
    return action


def _heading(menu: QMenu, text: str) -> None:
    menu.addAction(text).setEnabled(False)


def _plugin_items(menu: QMenu, canvas: "Canvas", section: str, target: MenuTarget) -> None:
    """Add the items plugins offer for this place in the menu."""
    for function in canvas.extensions.menu_items[section]:
        for item in function(canvas, target):
            _add_item(menu, item)


def _add_item(menu: QMenu, item: MenuItem) -> None:
    if item.items:
        submenu = menu.addMenu(item.title)
        submenu.setEnabled(item.enabled)
        for child in item.items:
            _add_item(submenu, child)
        return
    action = menu.addAction(item.title)
    action.setEnabled(item.enabled)
    if item.action is not None:
        run = item.action
        action.triggered.connect(lambda: run())


def show_context_menu(canvas: "Canvas", pos: QPointF, hit: Hit | None) -> None:
    """Open the menu for the atom or bond under the cursor (or the whole structure)."""
    menu = build_context_menu(canvas, hit)
    menu.exec(canvas.mapToGlobal(pos.toPoint()))
    menu.deleteLater()   # a child of the canvas; its submenus go with it


def build_context_menu(canvas: "Canvas", hit: Hit | None) -> QMenu:
    """Return the menu for the atom or bond under the cursor (or the whole structure), not yet shown."""
    menu = QMenu(canvas)
    target = MenuTarget(None if hit is None else first_atom(hit),
                        (hit.i, hit.j) if isinstance(hit, BondHit) else None)
    if isinstance(hit, AtomHit):
        _atom_items(menu, canvas, hit.index, target)
    elif isinstance(hit, BondHit):
        _bond_items(menu, canvas, hit.i, hit.j, target)
    _selection_items(menu, canvas, target)
    menu.addSeparator()
    _add(menu, "Fit to Window", partial(canvas.view.fit, animate=True))
    return menu


def _atom_items(menu: QMenu, canvas: "Canvas", i: int, target: MenuTarget) -> None:
    doc = canvas.doc
    mol = canvas.mol
    _heading(menu, f"{atom_label(mol, i)}  {element(mol.atoms[i].symbol).name}")
    change = menu.addMenu("Change Element")
    for symbol in COMMON_ELEMENTS:
        kind = default_element_fragment(symbol).kind
        if kind is not None:
            mutate = partial(doc.edit, building.mutate_atom, i, symbol, kind)
            _add(change, f"{symbol}    {element(symbol).name}", mutate)
    _add(menu, "Add Hydrogens Here", partial(doc.edit, building.add_hydrogens, [i]))
    _add(menu, "Add Valence", partial(doc.edit, building.add_valence, i))
    _add(menu, "Delete Atom", partial(doc.edit, building.delete_atom, i))
    _add(menu, "Rotate Around This Fragment", lambda: canvas.view.center_on(canvas.mol.component(i)))
    centers, _ = canvas.stereo_labels()
    if i in centers:
        _add(menu, f"Invert Stereocentre ({centers[i]})", partial(canvas.commands.invert_center, i))
    _plugin_items(menu, canvas, "atom", target)
    menu.addSeparator()


def _bond_items(menu: QMenu, canvas: "Canvas", i: int, j: int, target: MenuTarget) -> None:
    doc = canvas.doc
    _heading(menu, f"Bond {atom_label(canvas.mol, i)}-{atom_label(canvas.mol, j)}")
    current = canvas.mol.bond_order(i, j)
    for name, order in BOND_ORDERS:
        action = _add(menu, name, partial(doc.edit, Molecule.add_bond, i, j, order))
        action.setCheckable(True)
        action.setChecked(current == order)
    _add(menu, "Delete Bond", partial(doc.edit, Molecule.remove_bond, i, j))
    _, bonds = canvas.stereo_labels()
    label = bonds.get(canvas.mol.key(i, j))
    if label is not None:
        _add(menu, f"Flip Double Bond ({label})", partial(canvas.commands.flip_double_bond, i, j))
    _plugin_items(menu, canvas, "bond", target)
    menu.addSeparator()


def _selection_items(menu: QMenu, canvas: "Canvas", target: MenuTarget) -> None:
    atom = target.atom
    select = menu.addMenu("Select")
    if atom is not None:
        _add(select, "Residue", partial(canvas.commands.select_group, Group.RESIDUE, atom))
        _add(select, "Chain", partial(canvas.commands.select_group, Group.CHAIN, atom))
        select.addSeparator()
    for group, name in GROUP_ITEMS:
        _add(select, name, partial(canvas.commands.select_group, group))

    selected = bool(canvas.selection)
    show = menu.addMenu("Show Selection As" if selected else "Show Everything As")
    for rep, name in REPRESENTATION_NAMES.items():
        _add(show, name, partial(canvas.commands.set_representation, rep))
    menu.addSeparator()
    if selected:
        _add(menu, "Show Only Selection", canvas.commands.show_only_selection)
        _add(menu, "Hide Selection", canvas.commands.hide_selection)
    if canvas.commands.hidden_count():
        _add(menu, "Show Everything", canvas.commands.show_everything)
    if 2 <= len(canvas.selection) <= 4:
        _add(menu, "Pin Measurement", partial(canvas.commands.toggle_measurement, list(canvas.selection)))
    if canvas.mol.measurements:
        _add(menu, "Clear Measurements", canvas.commands.clear_measurements)
    _plugin_items(menu, canvas, "selection", target)
    menu.addSeparator()
    if selected:
        _add(menu, f"Delete Selection ({len(canvas.selection)})", canvas.commands.delete_selection)
    _add(menu, "Add Hydrogens to Selection" if selected else "Add Hydrogens", canvas.commands.add_hydrogens)
    _add(menu, "Clean Selection" if selected else "Clean Geometry", canvas.commands.clean_geometry)
