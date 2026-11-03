"""Atom list and keyboard shortcut dialogs."""

from PySide6.QtCore import Qt
from PySide6.QtGui import QFontDatabase, QKeySequence, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QHeaderView,
    QLabel,
    QMenuBar,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
)

from miles.core.elements import ELEMENTS, normalize_symbol
from miles.core.molecule import Molecule
from miles.render.theme import Theme
from miles.ui.controls import dialog_stylesheet
from miles.ui.document import Change, Document
from miles.ui.text import atom_label


def muted(widget: QLabel) -> QLabel:
    """Give a label the secondary text colour (of the theme in themed dialogs, else of the system)."""
    widget.setForegroundRole(QPalette.ColorRole.PlaceholderText)
    widget.setObjectName("Muted")
    return widget


class AtomListDialog(QDialog):
    """Editable table of element symbols and Cartesian coordinates."""

    HEADERS = ["Symbol", "X", "Y", "Z", "Bonded to"]
    SYMBOL, X, Y, Z, BONDED = range(5)

    def __init__(self, doc: Document, parent: QWidget | None = None, theme: Theme | None = None) -> None:
        super().__init__(parent)
        self.doc = doc
        self.setWindowTitle(f"Atoms - {doc.name}")
        self.resize(560, 480)
        layout = QVBoxLayout(self)
        self.table = QTableWidget(0, len(self.HEADERS))
        self.table.setHorizontalHeaderLabels(self.HEADERS)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        layout.addWidget(self.table)
        layout.addWidget(muted(QLabel("Edit a symbol or coordinate (angstrom) and press Return.")))
        self._editing = False
        self.table.itemChanged.connect(self._item_changed)
        doc.changed.connect(self.refresh)
        if theme is not None:
            self.set_theme(theme)
        self.refresh()

    def set_theme(self, theme: Theme) -> None:
        """Draw the dialog in the app's theme, which may differ from the system's."""
        self.setStyleSheet(dialog_stylesheet(theme))

    def refresh(self) -> None:
        """Fill the table from the molecule."""
        if self._editing:
            return
        mol = self.doc.mol
        self.table.blockSignals(True)
        self.table.setRowCount(len(mol))
        for i, a in enumerate(mol.atoms):
            values = [a.symbol, *(f"{v:.6f}" for v in mol.coords[i]),
                      " ".join(atom_label(mol, j) for j in mol.neighbors(i))]
            for col, text in enumerate(values):
                item = QTableWidgetItem(text)
                if col == self.BONDED:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                elif col != self.SYMBOL:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(i, col, item)
        self.table.setVerticalHeaderLabels([str(i + 1) for i in range(len(mol))])
        self.table.blockSignals(False)

    def _item_changed(self, item: QTableWidgetItem) -> None:
        row, col = item.row(), item.column()
        if row >= len(self.doc.mol):
            return
        self._editing = True
        try:
            if col == self.SYMBOL:
                symbol = normalize_symbol(item.text())
                if symbol in ELEMENTS:
                    self.doc.edit(_set_symbol, row, symbol)
            elif col in (self.X, self.Y, self.Z):
                self.doc.edit(_set_coordinate, row, col - self.X, float(item.text()), kind=Change.COORDS)
        except ValueError:
            pass
        finally:
            self._editing = False
            self.refresh()


def _set_symbol(mol: Molecule, i: int, symbol: str) -> None:
    mol.atoms[i].symbol = symbol


def _set_coordinate(mol: Molecule, i: int, axis: int, value: float) -> None:
    position = mol.coords[i].copy()
    position[axis] = value
    mol.set_positions([i], position)


SHORTCUTS = [
    ("Selecting", [
        ("Click an atom", "select it and show its details"),
        ("Shift click", "add or remove atoms; 2, 3 or 4 atoms show a bond, angle or dihedral"),
        ("Drag or scroll the dial", "change the value; Option-scroll anywhere does the same"),
        ("Drag a selected atom", "move the selection; Option-drag rotates it"),
        ("Double click", "select the residue, then the chain, then everything"),
        ("Delete or Backspace", "delete the selection, or the atom under the cursor"),
        ("Esc", "clear the selection, then go back to the select tool"),
    ]),
    ("Building", [
        ("Click empty space", "place the current fragment"),
        ("Click a hydrogen", "replace it with the fragment"),
        ("Drag from an atom", "grow the fragment in that direction"),
        ("Drag onto another atom", "bond them (closes rings)"),
        ("Click a heavy atom", "change its element (element fragments)"),
        ("Click a bond", "cycle single, double, triple"),
        ("Type a symbol: c, n, o, cl ...", "switch element; 1 - 9 picks its bonding type"),
        ("Right click", "select ligands, water, chains; show selection as cartoon, sticks ..."),
    ]),
    ("View", [
        ("Drag", "rotate"),
        ("Shift drag, middle drag", "move"),
        ("Scroll, pinch, right drag", "zoom"),
        ("Cmd drag, two-finger rotate", "spin about the screen axis"),
        ("R, T, S", "drag rotates, translates or scales (with the build tool, S is sulfur)"),
        ("Arrow keys, + and -", "turn and zoom"),
        (", and .", "previous and next frame of a trajectory or conformer search"),
        ("Space", "play or pause the frames"),
        ("Up and Down, Tab", "in the Python console: recall commands, complete names"),
    ]),
]
"""Mouse gestures and keys that are not menu commands; the dialog lists those from the menus."""


def menu_shortcuts(bar: QMenuBar) -> list[tuple[str, list[tuple[str, str]]]]:
    """Return the menu commands that have a shortcut, by top-level menu, in menu order.

    Parameters
    ----------
    bar : QMenuBar
        The window's menu bar, with the commands of the plugins loaded.

    Returns
    -------
    list of (str, list of (str, str))
        Menu title and its ``(shortcut, command)`` rows; the shortcut as the
        platform writes it, the command with its submenu path.
    """
    from miles.ui.command_palette import menu_entries

    sections: dict[str, list[tuple[str, str]]] = {}
    for entry in menu_entries(bar):
        sequence = entry.action.shortcut()
        if sequence.isEmpty():
            continue
        menu, _, path = entry.path.partition(" > ")
        key = sequence.toString(QKeySequence.SequenceFormat.NativeText)
        sections.setdefault(menu, []).append((key, path))
    return list(sections.items())


class ShortcutsDialog(QDialog):
    """Reference of mouse and keyboard shortcuts.

    Parameters
    ----------
    parent : QWidget, optional
        The main window.
    theme : Theme, optional
        Colours of the app.
    bar : QMenuBar, optional
        Menu bar whose commands with shortcuts are listed after the gestures,
        so the list matches the plugins loaded.
    """

    def __init__(self, parent: QWidget | None = None, theme: Theme | None = None,
                 bar: QMenuBar | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Shortcuts")
        self.resize(620, 640)
        if theme is not None:   # the app's theme, which may differ from the system's
            self.setStyleSheet(dialog_stylesheet(theme))
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        html = ["<div style='margin:14px 18px'>"]
        # one table per list, so the keys of all its sections share a column
        for sections in (SHORTCUTS, menu_shortcuts(bar) if bar is not None else []):
            html.append("<table cellspacing='5'>")
            for title, rows in sections:
                html.append(f"<tr><td colspan='2' style='padding-top:14px; font-weight:600'>"
                            f"{title}</td></tr>")
                html.extend(f"<tr><td style='padding-right:18px'><b>{key}</b></td><td>{what}</td></tr>"
                            for key, what in rows)
            html.append("</table>")
        html.append("</div>")
        text = QTextBrowser()
        text.setFrameShape(QTextBrowser.Shape.NoFrame)
        text.setHtml("".join(html))
        layout.addWidget(text)
