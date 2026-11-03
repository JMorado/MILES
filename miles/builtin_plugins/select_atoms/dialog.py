"""Selecting atoms with an expression of the selection language, with a live count and history."""

from typing import override

from PySide6.QtCore import QEvent, QObject, Qt, QTimer
from PySide6.QtGui import QFontDatabase, QKeyEvent
from PySide6.QtWidgets import QCheckBox, QComboBox, QHBoxLayout, QLineEdit

from miles import SelectionError, residue_ids
from miles.builtin_plugins.select_atoms.examples import FUNCTIONAL_GROUPS
from miles.plugin import ANGSTROM, App, DocumentView, PluginDialog, plural
from miles.plugin.widgets import label, theme_palette

EVALUATE_DELAY_MS = 120   # work the expression out this long after the last keystroke
HISTORY_LENGTH = 30

EXAMPLES: list[tuple[str, str]] = [
    ("Protein", "protein"),
    ("Ligands", "ligand"),
    (f"Pocket: residues within 5 {ANGSTROM} of ligands", "byres within 5 of ligand"),
    (f"Water within 3.5 {ANGSTROM} of ligands", "water within 3.5 of ligand"),
    ("Side chains lining the pocket", "sidechain and byres within 4 of ligand"),
    ("Protein backbone", "backbone"),
    ("Alpha carbons", "name CA"),
    ("Charged residues", "resname ASP GLU LYS ARG"),
    ("Polar hydrogens", "hydrogen and bonded to element N O"),
    ("Metal ions", "ion and metal"),
    ("Aromatic atoms", "aromatic"),
    ("Atoms not drawn", "hidden"),
]
"""Expressions offered under Examples: (description, expression)."""

FUNCTIONAL_GROUP_EXAMPLES: list[tuple[str, str]] = [(name, f'smarts "{pattern}"')
                                                    for name, pattern in FUNCTIONAL_GROUPS]
"""Common SMARTS patterns, offered under Examples after the expressions."""

history: list[str] = []
"""Expressions selected while MILES runs, the most recent first."""


def remember(expression: str) -> None:
    """Put an expression at the top of the history."""
    if expression in history:
        history.remove(expression)
    history.insert(0, expression)
    del history[HISTORY_LENGTH:]


class SelectionDialog(PluginDialog):
    """Expression field with examples and a live count; :attr:`atoms` holds the atoms to select.

    While the expression is valid, the atoms it picks are selected in the
    window so they can be seen; Cancel puts the old selection back. The Up
    and Down keys go through earlier expressions.

    Parameters
    ----------
    app : App
        The application, for its window and theme.
    doc : DocumentView
        Tab whose atoms are selected.
    """

    def __init__(self, app: App, doc: DocumentView) -> None:
        super().__init__(app, "Select")
        self.doc = doc
        self.theme = app.theme
        self.original = doc.selection
        self.found: list[int] = []
        self._recall = -1   # place in the history while going through it with the arrow keys
        self._typed = ""
        self._pending = QTimer(self)
        self._pending.setSingleShot(True)
        self._pending.setInterval(EVALUATE_DELAY_MS)
        self._pending.timeout.connect(self.evaluate)
        lay = self.body
        lay.addWidget(label("Select atoms where", "Title"))
        row = QHBoxLayout()
        row.setSpacing(8)
        self.field = QLineEdit()
        self.field.setObjectName("Field")
        self.field.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.field.setMinimumWidth(380)
        self.field.setPlaceholderText("e.g. byres within 5 of ligand")
        self.field.textEdited.connect(self._edited)
        self.field.installEventFilter(self)
        row.addWidget(self.field, 1)
        self.examples = QComboBox()
        self.examples.addItem("Examples")
        groups = [[(f"{name}  ({expression})", name) for name, expression in doc.selections.items()],
                  EXAMPLES, FUNCTIONAL_GROUP_EXAMPLES]
        for group in (g for g in groups if g):
            if self.examples.count() > 1:
                self.examples.insertSeparator(self.examples.count())
            for title, expression in group:
                self.examples.addItem(title, expression)
        # the button stays small; the list opens as wide as its longest line
        self.examples.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.examples.setMinimumContentsLength(len("Examples") + 1)
        view = self.examples.view()
        view.setMinimumWidth(view.sizeHintForColumn(0) + 32)
        self.examples.activated.connect(self._example_chosen)
        theme_palette(self.examples, self.theme)   # also themes the open list
        row.addWidget(self.examples)
        lay.addLayout(row)
        self.message = label("", "Muted")
        self.message.setWordWrap(True)
        lay.addWidget(self.message)
        self.add = QCheckBox("Add to the current selection")
        self.add.toggled.connect(lambda _: self._preview())
        lay.addWidget(self.add)
        self.select = self.buttons("Select")
        self.evaluate()

    @property
    def atoms(self) -> list[int]:
        """Atoms to select: those the expression picks, with the old selection when adding to it."""
        if self.add.isChecked():
            return sorted({*self.original, *self.found})
        return list(self.found)

    def set_expression(self, text: str) -> None:
        """Put an expression in the field and work it out at once."""
        self.field.setText(text)
        self.evaluate()

    def _edited(self, _: str) -> None:
        self._recall = -1
        self._pending.start()

    def _example_chosen(self, index: int) -> None:
        expression = self.examples.itemData(index)
        if isinstance(expression, str):
            self.set_expression(expression)
        self.examples.setCurrentIndex(0)
        self.field.setFocus()

    def evaluate(self) -> None:
        """Work the expression out now and say what it picks or what is wrong with it."""
        self._pending.stop()
        text = self.field.text().strip()
        error = False
        if not text:
            self.found = []
            note = "Type an expression such as protein and name CA, or choose one under Examples."
        else:
            try:
                self.found = self.doc.find(text, self.original)
            except SelectionError as exc:
                self.found, error = [], True
                note = str(exc)
            else:
                note = self._summary()
        self.message.setText(note)
        color = self.theme.danger if error else self.theme.muted
        self.message.setStyleSheet(f"color: {color.name()};")
        self.select.setEnabled(bool(self.found))
        self._preview()

    def _summary(self) -> str:
        if not self.found:
            return "No atom matches"
        mol = self.doc.molecule
        atoms = plural(len(self.found), "atom")
        if (mol.column("resname")[self.found] != "").any():
            residues = len(set(residue_ids(mol)[self.found].tolist()))
            return f"{atoms} in {plural(residues, 'residue')}"
        return atoms

    def _preview(self) -> None:
        self.doc.selection = self.atoms if self.found else self.original

    @override
    def accept(self) -> None:
        if not self.found:
            return
        remember(self.field.text().strip())
        super().accept()

    @override
    def reject(self) -> None:
        self.doc.selection = self.original
        super().reject()

    @override
    def eventFilter(self, watched: QObject, event: QEvent) -> bool:
        if watched is self.field and event.type() == QEvent.Type.KeyPress:
            key = event.key() if isinstance(event, QKeyEvent) else None
            if key in (Qt.Key.Key_Up, Qt.Key.Key_Down) and history:
                self._step_history(1 if key == Qt.Key.Key_Up else -1)
                return True
        return super().eventFilter(watched, event)

    def _step_history(self, step: int) -> None:
        if self._recall == -1:
            self._typed = self.field.text()
        self._recall = max(-1, min(len(history) - 1, self._recall + step))
        self.set_expression(history[self._recall] if self._recall >= 0 else self._typed)
