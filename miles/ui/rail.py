"""The vertical tool strip floating over the canvas."""

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtGui import QColor, QIcon
from PySide6.QtWidgets import QFrame, QToolButton, QVBoxLayout, QWidget

from miles.core.fragments import Category
from miles.render.painter import render_pixmap
from miles.render.theme import Theme
from miles.ui import icons
from miles.ui.builder import Builder, Tool
from miles.ui.controls import rule
from miles.ui.text import COMMAND

DARK_TEXT = QColor(25, 27, 31)
LIGHT_TEXT = QColor(255, 255, 255)


class Rail(QFrame):
    """Select, build and delete tools, then clean and symmetry."""

    fragmentRequested = Signal()
    toolRequested = Signal(object)
    cleanRequested = Signal()
    symmetryRequested = Signal()

    def __init__(self, builder: Builder, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Panel")
        self.builder = builder
        self.theme: Theme | None = None
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(6, 6, 6, 6)
        self._layout.setSpacing(4)

        self.select = self._button(f"Select: click an atom to inspect it, Shift-click to add  ({COMMAND} 1)",
                                   checkable=True)
        self.select.clicked.connect(lambda: self.toolRequested.emit(Tool.SELECT))
        self.build = self._button("Build: current fragment (click to choose)", checkable=True)
        self.build.setIconSize(QSize(26, 26))
        self.build.clicked.connect(self._build_clicked)
        self.delete = self._button(f"Delete atoms and bonds  ({COMMAND} 3)", checkable=True)
        self.delete.clicked.connect(lambda: self.toolRequested.emit(Tool.DELETE))
        self._layout.addWidget(rule())
        self.clean = self._button(f"Clean geometry  ({COMMAND} L)")
        self.clean.clicked.connect(self.cleanRequested.emit)
        self.symmetry = self._button(f"Point group and symmetrize  ({COMMAND} G)")
        self.symmetry.clicked.connect(self.symmetryRequested.emit)
        builder.toolChanged.connect(self.sync)
        builder.fragmentChanged.connect(self.sync)

    def _button(self, tip: str, checkable: bool = False) -> QToolButton:
        b = QToolButton()
        b.setObjectName("RailButton")
        b.setFixedSize(38, 38)
        b.setIconSize(QSize(22, 22))
        b.setToolTip(tip)
        b.setCheckable(checkable)
        b.setCursor(Qt.CursorShape.PointingHandCursor)
        self._layout.addWidget(b)
        return b

    def _build_clicked(self) -> None:
        self.builder.set_tool(Tool.BUILD)
        self.sync()
        self.fragmentRequested.emit()

    def set_theme(self, theme: Theme) -> None:
        """Redraw the icons in the theme's colours."""
        self.theme = theme
        for b, name in ((self.select, "pointer"), (self.delete, "delete"), (self.clean, "clean"),
                        (self.symmetry, "symmetry")):
            b.setIcon(icons.icon(name, theme.text))
        self.sync()

    def sync(self) -> None:
        """Show the current fragment on the build button and check the active tool."""
        if self.theme is None:
            return
        frag = self.builder.fragment
        if frag.category == Category.ELEMENT:
            fill = self.theme.atom_color(frag.symbol)
            text = DARK_TEXT if fill.lightness() > 150 else LIGHT_TEXT
            self.build.setIcon(icons.element_chip(frag.symbol, fill, text))
        else:
            self.build.setIcon(QIcon(render_pixmap(frag.mol, 26, 26, self.theme, fill=0.98)))
        self.build.setToolTip(f"Build with {frag.name.lower()}  ({COMMAND} 2, click to choose a fragment)")
        tool = self.builder.tool
        self.build.setChecked(tool == Tool.BUILD)
        self.select.setChecked(tool == Tool.SELECT)
        self.delete.setChecked(tool == Tool.DELETE)
