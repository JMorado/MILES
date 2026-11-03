"""Floating panel listing the properties of the structure or the selection; clicking a value copies it."""

from typing import NamedTuple, override

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QFontMetrics, QMouseEvent
from PySide6.QtWidgets import QApplication, QGridLayout, QHBoxLayout, QLabel, QVBoxLayout, QWidget

from miles.builtin_plugins.properties.descriptors import Properties, molecular_properties
from miles.plugin import (
    SEPARATOR,
    BackgroundResults,
    Change,
    DocumentPanel,
    DocumentView,
    plural,
    subscript_digits,
)
from miles.plugin.widgets import button, label, rule

ANGSTROM2 = "\u00c5\u00b2"
VALUE_WIDTH = 230      # pixels; longer values are elided in the middle
ASYNC_ATOMS = 5000     # larger structures are measured on a worker thread
UPDATE_DELAY_MS = 200  # recompute this long after the last selection or structure change


class Row(NamedTuple):
    """One line of the panel."""

    name: str
    text: str        # as shown
    copy: str        # as copied
    tip: str = ""


def _signed(charge: int) -> str:
    return f"+{charge}" if charge > 0 else str(charge)


def property_rows(props: Properties) -> list[list[Row]]:
    """Rows to show, in sections (composition, identifiers, drug-likeness); empty sections are left out.

    Parameters
    ----------
    props : Properties
        Computed properties.

    Returns
    -------
    list of list of Row
        Sections of rows.
    """
    composition = [
        Row("Formula", subscript_digits(props.formula), props.formula),
        Row("Molecular weight", f"{props.weight:.2f} g/mol", f"{props.weight:.2f} g/mol",
            "From standard atomic weights"),
        Row("Exact mass", f"{props.exact_mass:.5f} Da", f"{props.exact_mass:.5f} Da",
            "Monoisotopic: the most abundant isotope of each element"),
        Row("Atoms", f"{props.atoms} ({props.heavy_atoms} heavy)", str(props.atoms)),
    ]
    if props.charge is not None:
        composition.append(Row("Net charge", _signed(props.charge), _signed(props.charge),
                               "Sum of the formal charges worked out from the bonding"))
    identifiers = [Row(name, value, value) for name, value in
                   (("SMILES", props.smiles), ("InChI", props.inchi), ("InChIKey", props.inchikey)) if value]
    drug: list[Row] = []
    if props.logp is not None and props.tpsa is not None and props.molar_refractivity is not None:
        drug = [
            Row("logP", f"{props.logp:.2f}", f"{props.logp:.2f}",
                "Octanol-water partition coefficient (Wildman-Crippen)"),
            Row("TPSA", f"{props.tpsa:.1f} {ANGSTROM2}", f"{props.tpsa:.1f}",
                "Topological polar surface area (Ertl)"),
            Row("Molar refractivity", f"{props.molar_refractivity:.2f}", f"{props.molar_refractivity:.2f}",
                "Wildman-Crippen"),
            Row("H-bond donors", str(props.donors), str(props.donors), "N and O atoms with hydrogen"),
            Row("H-bond acceptors", str(props.acceptors), str(props.acceptors), "N and O atoms"),
            Row("Rotatable bonds", str(props.rotatable_bonds), str(props.rotatable_bonds)),
        ]
        violations = props.lipinski_violations or []
        text = ("none" if not violations
                else f"{plural(len(violations), 'violation')}: {', '.join(violations)}")
        drug.append(Row("Rule of five", text, text,
                        "Lipinski: weight up to 500, logP up to 5, up to 5 donors and 10 acceptors"))
    return [section for section in (composition, identifiers, drug) if section]


def property_text(props: Properties) -> str:
    """All rows as plain text, one ``name: value`` per line, for the clipboard."""
    return "\n".join(f"{row.name}: {row.copy}" for section in property_rows(props) for row in section)


def property_note(props: Properties) -> str:
    """Why some properties are missing, or an empty string."""
    if props.hydrogens_missing:
        return ("Hydrogens seem to be missing. Add them (Edit > Add Hydrogens) to see the charge, "
                "identifiers and descriptors.")
    if not props.small_molecule:
        return "Identifiers and drug-likeness are shown for small molecules; select a ligand to see them."
    return ""


class CopyValue(QLabel):
    """A value that copies itself when clicked."""

    clicked = Signal()

    def __init__(self, row: Row, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.row = row
        metrics = QFontMetrics(self.font())
        self.setText(metrics.elidedText(row.text, Qt.TextElideMode.ElideMiddle, VALUE_WIDTH))
        tip = row.text if self.text() != row.text else ""
        self.setToolTip("\n".join(t for t in (tip, row.tip, "Click to copy") if t))
        self.setCursor(Qt.CursorShape.PointingHandCursor)

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if e.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()


class PropertiesPanel(DocumentPanel):
    """Composition, identifiers and drug-likeness descriptors of the structure or the selection.

    The panel follows the current tab, the selection and structure edits,
    a short moment after the last change so that editing stays smooth.
    Moving atoms and playing frames do not change these properties and
    are ignored.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedWidth(400)
        self.props: Properties | None = None
        self._results: BackgroundResults | None = None   # properties of large structures
        self._pending = QTimer(self)
        self._pending.setSingleShot(True)
        self._pending.setInterval(UPDATE_DELAY_MS)
        self._pending.timeout.connect(self.refresh)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 12, 12, 14)
        lay.setSpacing(6)
        top = QHBoxLayout()
        top.addWidget(self.chevron)
        top.addWidget(label("Properties", "Title"))
        top.addStretch(1)
        self.copy_all = button("Link", "Copy All", "Copy every value as text")
        self.copy_all.clicked.connect(self._copy_all)
        top.addWidget(self.copy_all)
        close = button("Link", "Close")
        close.clicked.connect(self.close_panel)
        top.addWidget(close)
        lay.addLayout(top)
        self.body_layout.setSpacing(6)
        self.subtitle = label("", "Muted")
        self.body_layout.addWidget(self.subtitle)
        self.rows = QVBoxLayout()
        self.rows.setSpacing(6)
        self.body_layout.addLayout(self.rows)
        self.note = label("", "Muted")
        self.note.setWordWrap(True)
        self.body_layout.addWidget(self.note)
        self.status = label("", "Muted")
        self.status.hide()
        self.body_layout.addWidget(self.status)
        lay.addWidget(self.body)

    def fit(self) -> None:
        """Take the height the contents need at the panel's width (wrapped notes included)."""
        layout = self.layout()
        if layout is None:
            return
        layout.activate()
        height = (layout.totalHeightForWidth(self.width()) if layout.hasHeightForWidth()
                  else layout.totalSizeHint().height())
        self.resize(self.width(), height)
        layout.setGeometry(self.rect())   # place the rows now, even when the size did not change

    @override
    def closing(self) -> None:
        """Forget a pending update when the user closes the panel."""
        self._pending.stop()

    @override
    def selection_changed(self) -> None:
        """Describe the new selection, after a short pause."""
        if self.wanted:
            self._pending.start()

    @override
    def document_changed(self, kind: str) -> None:
        """Follow edits that change the composition (not moves or frames), after a short pause."""
        if self.wanted and kind in (Change.STRUCTURE, Change.BONDS):
            self._pending.start()

    @override
    def has_content(self, document: DocumentView) -> bool:
        """Whether the document has atoms."""
        return len(document.molecule) > 0

    @override
    def show_contents(self, document: DocumentView) -> None:
        """Show the properties of the selection, or of everything."""
        doc = document
        atoms = doc.selection_with_hydrogens()
        selected = doc.selection
        mol = doc.molecule.subset(atoms)
        what = "Selection" if selected else doc.name
        added = len(atoms) - len(selected) if selected else 0
        extra = f" (with {plural(added, 'bonded hydrogen')})" if added else ""
        self.subtitle.setText(f"{what}{SEPARATOR}{plural(len(mol), 'atom')}{extra}")
        self.status.hide()
        if len(mol) <= ASYNC_ATOMS:
            self._fill(molecular_properties(mol))
            return
        if self._results is None and self.context is not None:
            self._results = self.context.background_results(self._arrived, "Computing the properties",
                                                            name="properties")
        if self._results is None:
            return
        key = (doc.layout_revision, tuple(atoms))
        props = self._results.get(doc, key, molecular_properties, mol)
        if self._results.has(doc, key) and props is not None:
            self._fill(props)
        else:
            self._show_message("Computing ...")

    def _arrived(self, document: DocumentView) -> None:
        if document == self.document:
            self.refresh()

    def _clear(self) -> None:
        while self.rows.count():
            item = self.rows.takeAt(0)
            widget = item.widget() if item is not None else None
            if widget is not None:
                widget.hide()   # gone from view now, deleted once control returns to the event loop
                widget.deleteLater()

    def _show_message(self, text: str) -> None:
        self._clear()
        self.props = None
        self.note.setText(text)
        self.note.setVisible(True)
        self.copy_all.setEnabled(False)
        self.fit()
        self.relayout.emit()

    def _fill(self, props: Properties) -> None:
        self._clear()
        self.props = props
        box = QWidget()
        grid = QGridLayout(box)   # one grid, so the values line up across the sections
        grid.setContentsMargins(0, 2, 0, 2)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(4)
        line = 0
        for n, section in enumerate(property_rows(props)):
            if n:
                grid.addWidget(rule(), line, 0, 1, 2)
                line += 1
            for row in section:
                grid.addWidget(label(row.name, "Muted"), line, 0)
                value = CopyValue(row)
                value.clicked.connect(lambda r=row: self._copy(r.name, r.copy))
                grid.addWidget(value, line, 1)
                line += 1
        grid.setColumnStretch(1, 1)
        self.rows.addWidget(box)
        box.show()   # now, so that fit() measures the rows (Qt would show it only later)
        note = property_note(props)
        self.note.setText(note)
        self.note.setVisible(bool(note))
        self.copy_all.setEnabled(True)
        self.fit()
        self.relayout.emit()

    def _copy(self, name: str, text: str) -> None:
        QApplication.clipboard().setText(text)
        self.status.setText(f"Copied {name}")
        self.status.show()
        self.fit()

    def _copy_all(self) -> None:
        if self.props is not None:
            self._copy("all properties", property_text(self.props))
