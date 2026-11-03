"""Small bar shown while commands work in the background."""

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QFrame, QHBoxLayout, QWidget

from miles.ui.controls import button, label
from miles.ui.tasks import TaskRunner


class BusyBar(QFrame):
    """Says what is running in the background and offers to cancel it.

    Parameters
    ----------
    runner : TaskRunner
        The window's tasks.
    parent : QWidget, optional
        The stage.
    """

    relayout = Signal()

    def __init__(self, runner: TaskRunner, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Panel")
        self.runner = runner
        lay = QHBoxLayout(self)
        lay.setContentsMargins(14, 5, 6, 5)
        lay.setSpacing(8)
        self.text = label("", "Muted")
        lay.addWidget(self.text)
        self.cancel = button("Link", "Cancel")
        self.cancel.clicked.connect(runner.cancel_last)
        lay.addWidget(self.cancel)
        runner.changed.connect(self.refresh)
        self.hide()

    def refresh(self) -> None:
        """Show the task started last, and how many others run, or hide."""
        tasks = self.runner.tasks
        if tasks:
            others = f" (and {len(tasks) - 1} more)" if len(tasks) > 1 else ""
            self.text.setText(f"{tasks[-1].message} ...{others}")
            self.adjustSize()
        self.setVisible(bool(tasks))
        self.relayout.emit()
