"""Combining two grids into their difference or sum."""

from PySide6.QtWidgets import QComboBox, QGridLayout, QLineEdit

from miles.builtin_plugins.grids.combine import Combination
from miles.plugin import App, PluginDialog
from miles.plugin.widgets import Segmented, label, theme_palette

OPERATIONS = [(Combination.DIFFERENCE, "A - B"), (Combination.SUM, "A + B")]


class CombineDialog(PluginDialog):
    """Two grids, an operation and a name for the result.

    Read :attr:`first`, :attr:`second`, :attr:`combination` and :attr:`name`
    after :meth:`run` accepts.

    Parameters
    ----------
    app : App
        The application.
    names : list of str
        Names of the document's grids, in order.
    first : int, optional
        Grid chosen as A at the start (the one on screen); B starts at the next.
    """

    def __init__(self, app: App, names: list[str], first: int = 0) -> None:
        super().__init__(app, "Combine Grids")
        self._edited = False
        grid = QGridLayout()
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(8)
        self.a, self.b = QComboBox(), QComboBox()
        for row, (text, combo) in enumerate((("A", self.a), ("B", self.b))):
            combo.addItems(names)
            theme_palette(combo, app.theme)
            combo.currentIndexChanged.connect(self._suggest_name)
            grid.addWidget(label(text), row, 0)
            grid.addWidget(combo, row, 1)
        self.operation = Segmented([text for _, text in OPERATIONS])
        self.operation.changed.connect(self._suggest_name)
        grid.addWidget(label("Make"), 2, 0)
        grid.addWidget(self.operation, 2, 1)
        self.name_field = QLineEdit()
        self.name_field.setObjectName("Field")
        self.name_field.textEdited.connect(self._name_edited)
        grid.addWidget(label("Name"), 3, 0)
        grid.addWidget(self.name_field, 3, 1)
        self.body.addLayout(grid)
        note = label("B is interpolated at the points of A when the grids differ.", "Muted")
        note.setWordWrap(True)
        self.body.addWidget(note)
        self.ok = self.buttons("Combine")
        self.a.setCurrentIndex(max(0, min(first, len(names) - 1)))
        self.b.setCurrentIndex((self.a.currentIndex() + 1) % max(len(names), 1))
        self._suggest_name()

    @property
    def first(self) -> int:
        """Index of grid A."""
        return self.a.currentIndex()

    @property
    def second(self) -> int:
        """Index of grid B."""
        return self.b.currentIndex()

    @property
    def combination(self) -> Combination:
        """The operation chosen."""
        return OPERATIONS[self.operation.index()][0]

    @property
    def name(self) -> str:
        """Name for the new grid."""
        return self.name_field.text().strip()

    def _name_edited(self, _text: str) -> None:
        self._edited = True

    def _suggest_name(self, *_args: object) -> None:
        """Name the result after the grids and the operation, until the user types a name."""
        if not self._edited:
            symbol = OPERATIONS[self.operation.index()][1][2]
            self.name_field.setText(f"{self.a.currentText()} {symbol} {self.b.currentText()}")
