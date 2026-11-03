"""Bottom panel for the selection: atom details, or an editable bond, angle or dihedral."""

from PySide6.QtCore import QSize, Signal
from PySide6.QtWidgets import (
    QButtonGroup,
    QFrame,
    QHBoxLayout,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from miles.core import internal
from miles.core.elements import element
from miles.core.geometry import GEOMETRY_NAMES, angle, dihedral, distance
from miles.core.internal import MoveMode
from miles.core.molecule import Molecule, Representation
from miles.render.theme import Theme
from miles.ui import icons
from miles.ui.builder import Tool
from miles.ui.canvas import Canvas
from miles.ui.controls import Dial, ScrubField, button, label
from miles.ui.document import Change
from miles.ui.text import (
    ANGSTROM,
    DASH,
    DEGREE,
    DOT,
    SEPARATOR,
    atom_label,
    plural,
    residue_label,
    subscript_digits,
)

BOND_ORDERS = [(1.0, "Single"), (2.0, "Double"), (3.0, "Triple"), (1.5, "Aromatic"), (0.0, "No bond")]
MOVES = [(MoveMode.LAST, "moves {last} side"), (MoveMode.FIRST, "moves {first} side"),
         (MoveMode.BOTH, "moves both sides"), (MoveMode.ATOM, "moves {last} only")]
SHOW_AS = [(Representation.BALL, "Ball"), (Representation.TUBE, "Tube"), (Representation.WIRE, "Wire"),
           (Representation.SPACEFILL, "Space"), (Representation.CARTOON, "Cartoon"),
           (Representation.HIDDEN, "Hide")]
MAX_NEIGHBORS_SHOWN = 6

# Field range, field step per pixel, decimals; dial value per pixel, minor and major ticks, decimals.
DISTANCE_SCALE = ((0.3, 12.0, 0.004, 4), (0.002, 0.01, 0.1, 1))
ANGLE_SCALE = ((0.0, 180.0, 0.25, 2), (0.2, 1.0, 10.0, 0))
DIHEDRAL_SCALE = ((-180.0, 180.0, 0.5, 2), (0.25, 1.0, 10.0, 0))


class SelectionPanel(QFrame):
    """Floating panel shown with the select tool.

    One atom shows its details, 2-4 atoms the bond length, angle or dihedral
    they define (editable with the field and the dial), more atoms their
    formula. Larger selections can be shown in another representation.
    """

    relayout = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Panel")
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.canvas: Canvas | None = None
        self._move = 0
        self._checkpointed = False   # undo step recorded for the current edit
        self._editing = False        # ignore our own change notifications

        outer = QVBoxLayout(self)
        outer.setContentsMargins(14, 8, 10, 6)
        outer.setSpacing(2)
        row = QHBoxLayout()
        row.setSpacing(10)
        outer.addLayout(row)
        self.dial = Dial()
        self.dial.changed.connect(self._dialed)
        outer.addWidget(self.dial)
        self.show_as = self._show_as_row()
        outer.addWidget(self.show_as)

        self.atoms = label("", "Muted")
        self.field = ScrubField()
        self.field.edited.connect(self._edited)
        self.unit = label("", "Muted")
        self.info = label()
        self.orders, self.order_buttons = self._order_row()
        self.move_button = button("Link", tip="Which part of the molecule moves")
        self.move_button.clicked.connect(self._cycle_move)
        self.pin = button("Order", tip="Keep this measurement on screen and in pictures", checkable=True)
        self.pin.setFixedSize(30, 26)
        self.pin.setIconSize(QSize(18, 18))
        self.pin.clicked.connect(self._pin_clicked)
        for widget in (self.atoms, self.field, self.unit, self.info, self.orders, self.move_button, self.pin):
            row.addWidget(widget)

    def _show_as_row(self) -> QWidget:
        widget = QWidget()
        lay = QHBoxLayout(widget)
        lay.setContentsMargins(0, 2, 0, 0)
        lay.setSpacing(2)
        lay.addWidget(label("Show as", "Muted"))
        for rep, text in SHOW_AS:
            b = button("Seg", text)
            b.clicked.connect(lambda _=False, r=rep: self._show_as(r))
            lay.addWidget(b)
        lay.addStretch(1)
        return widget

    def _order_row(self) -> tuple[QWidget, dict[float, QToolButton]]:
        widget = QWidget()
        lay = QHBoxLayout(widget)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(1)
        group = QButtonGroup(self)
        buttons = {}
        for n, (order, tip) in enumerate(BOND_ORDERS):
            b = button("Order", tip=tip, checkable=True)
            b.setFixedSize(30, 26)
            b.setIconSize(QSize(20, 20))
            group.addButton(b, n)
            buttons[order] = b
            lay.addWidget(b)
        group.idClicked.connect(self._order_clicked)
        return widget, buttons

    def set_theme(self, theme: Theme) -> None:
        """Recolour the dial and the bond icons."""
        self.dial.theme = theme
        for order, b in self.order_buttons.items():
            b.setIcon(icons.bond_order_icon(order, theme.text))
        self.pin.setIcon(icons.icon("pin", theme.text))

    def set_canvas(self, canvas: Canvas) -> None:
        """Follow another document's canvas."""
        if canvas is not self.canvas:
            self._checkpointed = False
        self.canvas = canvas

    def nudge(self, steps: float) -> None:
        """Step the value by dial ticks (Option-scroll on the canvas)."""
        if self.isVisible() and self.dial.isVisible():
            self.dial.step(steps)

    def selection_changed(self) -> None:
        """Start a new undo step and show the new selection."""
        self._checkpointed = False
        self.refresh()

    def doc_changed(self, kind: str) -> None:
        """Show the new values after an edit made elsewhere."""
        if self._editing:
            return
        if kind == Change.STRUCTURE:
            self._checkpointed = False
        self.refresh()

    def _selection(self) -> list[int]:
        if self.canvas is None:
            return []
        n = len(self.canvas.mol)
        return [i for i in self.canvas.selection if i < n]

    # -- display ---------------------------------------------------------------------------------

    def refresh(self) -> None:
        """Update for the current canvas, selection and tool (hiding when there is nothing to show)."""
        before = (self.isHidden(), self.size())
        canvas = self.canvas
        selection = self._selection()
        if canvas is None or canvas.builder.tool != Tool.SELECT or not selection:
            self.hide()
            return
        mol = canvas.mol
        n = len(selection)
        editable = 2 <= n <= 4
        names = [atom_label(mol, i) for i in selection]
        self.atoms.setText(f" {DASH} ".join(names) if n <= 4 else f"{n} atoms")
        for widget in (self.field, self.unit, self.move_button, self.dial, self.pin):
            widget.setVisible(editable)
        key = tuple(selection)
        self.pin.setChecked(key in mol.measurements or key[::-1] in mol.measurements)
        self.info.setVisible(not editable)
        self.show_as.setVisible(not editable)
        self.orders.setVisible(n == 2)
        if n == 1:
            self.info.setText(self._atom_details(mol, selection[0]))
        elif n > 4:
            self.info.setText(self._group_details(mol, selection))
        else:
            self._show_measurement(mol, selection)
            first, last = names[0], names[-1]
            self.move_button.setText(MOVES[self._move][1].format(first=first, last=last))
        self.adjustSize()
        self.show()
        self.raise_()
        if (self.isHidden(), self.size()) != before:   # placed again only when it appeared or resized
            self.relayout.emit()

    @staticmethod
    def _atom_details(mol: Molecule, i: int) -> str:
        a = mol.atoms[i]
        parts = [element(a.symbol).name]
        if a.resname:
            parts.insert(0, f"{residue_label(a)}  {a.name}")
        if a.symbol != "H" and a.geom:
            parts.append(GEOMETRY_NAMES.get(a.geom, a.geom))
        neighbors = mol.neighbors(i)
        if neighbors:
            shown = ", ".join(atom_label(mol, j) for j in neighbors[:MAX_NEIGHBORS_SHOWN])
            parts.append(f"bonded to {shown}{', ...' if len(neighbors) > MAX_NEIGHBORS_SHOWN else ''}")
        else:
            parts.append("not bonded")
        x, y, z = mol.coords[i]
        return SEPARATOR.join(parts) + f"      {x:.4f}   {y:.4f}   {z:.4f}"

    @staticmethod
    def _group_details(mol: Molecule, atoms: list[int]) -> str:
        sub = mol.subset(atoms)
        residues = {(a.chain, a.resseq, a.resname) for a in sub.atoms if a.resname}
        extra = f"   {DOT}   {plural(len(residues), 'residue')}" if residues else ""
        return subscript_digits(sub.formula()) + extra

    def _show_measurement(self, mol: Molecule, selection: list[int]) -> None:
        c = mol.coords[selection]
        if len(selection) == 2:
            value, scale, unit, wrap = distance(c[0], c[1]), DISTANCE_SCALE, ANGSTROM, False
            order = mol.bond_order(selection[0], selection[1])
            for o, b in self.order_buttons.items():
                b.setChecked(o == order)
        elif len(selection) == 3:
            value, scale, unit, wrap = angle(*c), ANGLE_SCALE, DEGREE, False
        else:
            value, scale, unit, wrap = dihedral(*c), DIHEDRAL_SCALE, DEGREE, True
        (lo, hi, step, decimals), (per_px, minor, major, dial_decimals) = scale
        self.field.configure(lo, hi, step, decimals, wrap=wrap)
        self.field.set_value(value)
        self.dial.configure(lo, hi, per_px, minor, major, dial_decimals, wrap=wrap)
        self.dial.set_value(value)
        self.unit.setText(unit)

    # -- editing ---------------------------------------------------------------------------------

    def _show_as(self, rep: Representation) -> None:
        selection = self._selection()
        if self.canvas is not None and selection:
            self.canvas.commands.set_representation(rep, selection)

    def _cycle_move(self) -> None:
        self._move = (self._move + 1) % len(MOVES)
        self.refresh()

    def _dialed(self, value: float) -> None:
        self.field.set_value(value)
        self._edited(value)

    def _edited(self, value: float) -> None:
        if abs(self.dial.value - value) > 1e-9:
            self.dial.set_value(value)
        selection = self._selection()
        if self.canvas is None or not 2 <= len(selection) <= 4:
            return
        doc = self.canvas.doc
        if not self._checkpointed:
            doc.checkpoint()
            self._checkpointed = True
        mode = MOVES[self._move][0]
        match selection:
            case [a, b]:
                internal.set_bond_length(doc.mol, a, b, value, mode)
            case [a, b, c]:
                internal.set_angle(doc.mol, a, b, c, value, mode)
            case [a, b, c, d]:
                internal.set_dihedral(doc.mol, a, b, c, d, value, mode)
        self._editing = True
        try:
            doc.notify(Change.COORDS)
        finally:
            self._editing = False

    def _pin_clicked(self) -> None:
        selection = self._selection()
        if self.canvas is not None and 2 <= len(selection) <= 4:
            self.pin.setChecked(self.canvas.commands.toggle_measurement(selection))

    def _order_clicked(self, index: int) -> None:
        selection = self._selection()
        if self.canvas is None or len(selection) != 2:
            return
        i, j = selection
        order = BOND_ORDERS[index][0]
        if order == 0:
            self.canvas.doc.edit(Molecule.remove_bond, i, j, kind=Change.BONDS)
        else:
            self.canvas.doc.edit(Molecule.add_bond, i, j, order, kind=Change.BONDS)
        self._checkpointed = False
