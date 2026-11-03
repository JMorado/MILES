"""The dialog that offers documents autosaved before MILES stopped unexpectedly."""

import datetime
import os
from typing import override

from PySide6.QtCore import QDateTime, QLocale, Qt
from PySide6.QtGui import QShowEvent
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from miles.render.theme import Theme
from miles.ui.autosave import Leftover
from miles.ui.controls import dialog_stylesheet, label
from miles.ui.text import plural

DISCARD = 2
"""Result of the dialog when the user chose Discard (Restore accepts it, closing it rejects it)."""
MAX_LIST_HEIGHT = 300


def when(seconds: float, now: datetime.datetime | None = None) -> str:
    """Time of an autosave as it reads in the list: "Today 14:32", "Yesterday 09:05" or a date."""
    moment = datetime.datetime.fromtimestamp(seconds)
    today = (now or datetime.datetime.now()).date()
    clock = moment.strftime("%H:%M")
    if moment.date() == today:
        return f"Today {clock}"
    if moment.date() == today - datetime.timedelta(days=1):
        return f"Yesterday {clock}"
    date = QDateTime.fromSecsSinceEpoch(int(seconds))
    return QLocale.system().toString(date, QLocale.FormatType.ShortFormat)


class RecoveryDialog(QDialog):
    """List of autosaved documents with Restore and Discard.

    Closing the dialog without choosing keeps the documents for the next start.

    Parameters
    ----------
    leftovers : list of Leftover
        Documents found, newest first.
    theme : Theme
        Colours.
    parent : QWidget, optional
        Window to show it over.
    """

    def __init__(self, leftovers: list[Leftover], theme: Theme, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.leftovers = leftovers
        self.setWindowTitle("Restore Unsaved Documents")
        self.setStyleSheet(dialog_stylesheet(theme))
        lay = QVBoxLayout(self)
        lay.setContentsMargins(20, 18, 20, 16)
        lay.setSpacing(10)
        lay.addWidget(label(f"Restore {plural(len(leftovers), 'unsaved document')}?", "Title"))
        text = label("MILES stopped before the changes below were saved. Restore them to carry on "
                     "where you left off, or discard them. Closing this window keeps them for the next "
                     "start.")
        text.setWordWrap(True)
        lay.addWidget(text)
        self.table = QTableWidget(len(leftovers), 2)
        self.table.setHorizontalHeaderLabels(["Document", "Autosaved"])
        self.table.verticalHeader().hide()
        self.table.setShowGrid(False)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.table.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        header.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        for row, leftover in enumerate(leftovers):
            name = QTableWidgetItem(leftover.name)
            folder = os.path.dirname(leftover.path) if leftover.path else ""
            name.setToolTip(leftover.path or "A new document that was never saved")
            if folder:
                name.setText(f"{leftover.name}   ({os.path.basename(folder) or folder})")
            self.table.setItem(row, 0, name)
            self.table.setItem(row, 1, QTableWidgetItem(when(leftover.time)))
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        lay.addWidget(self.table)
        buttons = QHBoxLayout()
        self.discard_button = QPushButton("Discard")
        self.discard_button.setToolTip("Delete the autosaved copies; the files on disk stay as they are")
        self.discard_button.clicked.connect(lambda: self.done(DISCARD))
        buttons.addWidget(self.discard_button)
        buttons.addStretch(1)
        self.restore_button = QPushButton("Restore")
        self.restore_button.setObjectName("Primary")
        self.restore_button.setDefault(True)
        self.restore_button.clicked.connect(self.accept)
        buttons.addWidget(self.restore_button)
        lay.addLayout(buttons)
        self.resize(520, self.sizeHint().height())

    @override
    def showEvent(self, event: QShowEvent) -> None:
        super().showEvent(event)
        self.fit_list()

    def fit_list(self) -> None:
        """Make the list as tall as its rows (up to a limit), now that the style has set their height."""
        rows = sum(self.table.rowHeight(row) for row in range(self.table.rowCount()))
        height = self.table.horizontalHeader().height() + rows + 2 * self.table.frameWidth()
        self.table.setFixedHeight(min(height, MAX_LIST_HEIGHT))
        self.resize(self.width(), self.sizeHint().height())
