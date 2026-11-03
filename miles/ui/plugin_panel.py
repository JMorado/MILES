"""Floating panel holding the contents a plugin provides."""

from PySide6.QtWidgets import QHBoxLayout, QVBoxLayout, QWidget

from miles.ui.controls import button, label
from miles.ui.floating import FloatingPanel


class PluginPanel(FloatingPanel):
    """A titled, movable and foldable panel around a plugin's widget.

    The plugin's widget goes in :attr:`body_layout`; the panel starts
    hidden and :meth:`~miles.ui.floating.FloatingPanel.toggle` (run by the
    plugin's menu command) shows or hides it. Its ``opened`` and ``closed``
    signals say when it appears and disappears.

    Parameters
    ----------
    title : str
        Shown in the panel's top row.
    parent : QWidget, optional
        The stage.
    """

    def __init__(self, title: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.title = title
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(6)
        top = QHBoxLayout()
        top.addWidget(self.chevron)
        top.addWidget(label(title, "Title"))
        top.addStretch(1)
        close = button("Link", "Close")
        close.clicked.connect(self.close_panel)
        top.addWidget(close)
        lay.addLayout(top)
        self.body_layout.setSpacing(6)
        lay.addWidget(self.body)
        self.hide()
