"""Popover with the point group and the Symmetrize button."""

from PySide6.QtCore import QPoint, Qt
from PySide6.QtWidgets import QFrame, QPushButton, QVBoxLayout, QWidget

from miles.core.molecule import Molecule
from miles.core.symmetry import find_point_group, symmetrize
from miles.ui.controls import Segmented, label
from miles.ui.document import Change, Document
from miles.ui.text import SEPARATOR, plural

TOLERANCES = [("Tight", 0.01), ("Normal", 0.1), ("Loose", 0.3)]
EXACT_TOLERANCE = 0.002


class SymmetryPopover(QWidget):
    """Shows the point group at a chosen tolerance and snaps the structure to it."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent, Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setFixedWidth(240)
        self.doc: Document | None = None
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        frame = QFrame()
        frame.setObjectName("Popover")
        outer.addWidget(frame)
        lay = QVBoxLayout(frame)
        lay.setContentsMargins(16, 14, 16, 16)
        lay.setSpacing(6)
        lay.addWidget(label("Point group", "Muted"))
        self.group = label("C1", "Big")
        lay.addWidget(self.group)
        self.detail = label("", "Muted")
        lay.addWidget(self.detail)
        lay.addSpacing(6)
        lay.addWidget(label("Tolerance", "Muted"))
        self.tolerance = Segmented([name for name, _ in TOLERANCES])
        self.tolerance.set_index(1)
        self.tolerance.changed.connect(lambda _: self.refresh())
        lay.addWidget(self.tolerance)
        lay.addSpacing(8)
        self.apply = QPushButton("Symmetrize")
        self.apply.setObjectName("Primary")
        self.apply.setCursor(Qt.CursorShape.PointingHandCursor)
        self.apply.clicked.connect(self._symmetrize)
        lay.addWidget(self.apply)

    def _tolerance(self) -> float:
        return TOLERANCES[self.tolerance.index()][1]

    def refresh(self) -> None:
        """Recompute the point group."""
        if self.doc is None or not len(self.doc.mol):
            self.group.setText("-")
            self.detail.setText("No atoms")
            self.apply.setEnabled(False)
            return
        mol = self.doc.mol
        group = find_point_group(mol.symbols, mol.coords, self._tolerance())
        exact = find_point_group(mol.symbols, mol.coords, EXACT_TOLERANCE).name
        self.group.setText(group.name)
        extra = "" if exact == group.name else f"{SEPARATOR}exact: {exact}"
        self.detail.setText(plural(group.order, "operation") + extra)
        self.apply.setEnabled(group.order > 1 and exact != group.name)

    def _symmetrize(self) -> None:
        if self.doc is None:
            return
        coords, _ = symmetrize(self.doc.mol.symbols, self.doc.mol.coords, self._tolerance())

        def apply(mol: Molecule) -> None:
            mol.coords = coords

        self.doc.edit(apply, kind=Change.COORDS)
        self.refresh()

    def show_for(self, doc: Document, pos: QPoint) -> None:
        """Open for a document at a screen position."""
        self.doc = doc
        self.refresh()
        self.adjustSize()
        self.move(pos)
        self.show()
