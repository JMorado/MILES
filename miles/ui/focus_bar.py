"""Small bar shown while focus mode hides part of the structure."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFrame, QHBoxLayout, QWidget

from miles.ui.canvas import Canvas
from miles.ui.controls import button, label


class FocusBar(QFrame):
    """Says how many atoms are shown and offers to show everything again."""

    relayout = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Panel")
        self.canvas: Canvas | None = None
        lay = QHBoxLayout(self)
        lay.setContentsMargins(14, 5, 6, 5)
        lay.setSpacing(8)
        self.text = label("", "Muted")
        lay.addWidget(self.text)
        self.show_all = button("Link", "Show Everything")
        self.show_all.clicked.connect(self._show_everything)
        lay.addWidget(self.show_all)
        self.hide()

    def set_canvas(self, canvas: Canvas | None) -> None:
        """Follow another canvas."""
        self.canvas = canvas
        self.refresh()

    def refresh(self) -> None:
        """Show or hide the bar for the current structure (the stage places it again only if that changed)."""
        before = (self.isHidden(), self.size())
        canvas = self.canvas
        hidden = canvas.commands.hidden_count() if canvas is not None else 0
        if canvas is None or not hidden:
            self.hide()
        else:
            total = len(canvas.mol)
            self.text.setText(f"Showing {total - hidden:,} of {total:,} atoms")
            self.adjustSize()
            self.show()
        if (self.isHidden(), self.size()) != before:
            self.relayout.emit()

    def _show_everything(self) -> None:
        if self.canvas is not None:
            self.canvas.commands.show_everything()
