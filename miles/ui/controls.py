"""Small controls and the stylesheet shared by the floating panels."""

import math
from pathlib import Path
from typing import override

from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QColor, QKeyEvent, QMouseEvent, QPainter, QPaintEvent, QPalette, QWheelEvent
from PySide6.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QToolButton,
    QWidget,
)

from miles.render.theme import Theme

FINE = 0.2
"""Drag speed factor while Option is held."""

CHECK_IMAGE = (Path(__file__).parent / "images" / "check.svg").as_posix()
"""White tick drawn on checked check boxes in dialogs."""


def _rgba(c: QColor, alpha: float) -> str:
    return f"rgba({c.red()}, {c.green()}, {c.blue()}, {alpha:.3f})"


def stylesheet(t: Theme) -> str:
    """Qt stylesheet for panels, popovers and their controls in a theme."""
    border = _rgba(t.panel_border, t.panel_border.alphaF())
    return f"""
    #Panel {{ background: {_rgba(t.panel, 0.94)}; border: 1px solid {border}; border-radius: 13px; }}
    #Popover {{ background: {t.panel.name()}; border: 1px solid {border}; border-radius: 14px; }}
    #Panel QLabel, #Popover QLabel {{ color: {t.text.name()}; background: transparent; }}
    #Panel QComboBox {{ color: {t.text.name()}; }}
    QLabel#Muted {{ color: {t.muted.name()}; }}
    QLabel#Big {{ color: {t.text.name()}; font-size: 26px; font-weight: 600; }}
    QLabel#Title {{ color: {t.text.name()}; font-weight: 600; }}
    QFrame#Rule {{ background: {_rgba(t.text, 0.10)}; border: none; }}
    QToolButton#RailButton {{ border: none; border-radius: 9px; background: transparent; }}
    QToolButton#RailButton:hover {{ background: {_rgba(t.text, 0.07)}; }}
    QToolButton#RailButton:checked {{ background: {_rgba(t.accent, 0.16)}; }}
    QToolButton#Seg {{ border: none; border-radius: 7px; padding: 4px 11px; color: {t.muted.name()};
                       background: transparent; font-weight: 500; }}
    QToolButton#Seg:hover {{ color: {t.text.name()}; }}
    QToolButton#Seg:checked {{ background: {_rgba(t.text, 0.08)}; color: {t.text.name()}; }}
    QToolButton#Cell {{ border: none; border-radius: 5px; color: {t.text.name()}; background: transparent;
                        font-size: 11px; }}
    QToolButton#Cell:hover {{ background: {_rgba(t.text, 0.09)}; }}
    QToolButton#Cell:checked {{ background: {t.accent.name()}; color: white; }}
    QToolButton#Tile {{ border: none; border-radius: 10px; color: {t.text.name()}; background: transparent;
                        padding: 4px 2px; font-size: 11px; }}
    QToolButton#Tile:hover {{ background: {_rgba(t.text, 0.06)}; }}
    QToolButton#Tile:checked {{ background: {_rgba(t.accent, 0.14)}; }}
    QToolButton#Link {{ border: none; color: {t.muted.name()}; background: transparent; padding: 2px 6px;
                        border-radius: 6px; }}
    QToolButton#Link:hover {{ color: {t.text.name()}; background: {_rgba(t.text, 0.07)}; }}
    QToolButton#Order {{ border: none; border-radius: 6px; background: transparent; padding: 2px; }}
    QToolButton#Order:hover {{ background: {_rgba(t.text, 0.07)}; }}
    QToolButton#Order:checked {{ background: {_rgba(t.accent, 0.18)}; }}
    QLineEdit#Scrub {{ border: 1px solid transparent; border-radius: 8px; background: {_rgba(t.text, 0.06)};
                       color: {t.text.name()}; padding: 3px 8px; font-size: 15px; font-weight: 600;
                       selection-background-color: {t.accent.name()}; }}
    QLineEdit#Scrub:focus {{ border: 1px solid {t.accent.name()}; background: {t.panel.name()}; }}
    QScrollArea#Tiles, QWidget#TilesInner {{ background: transparent; border: none; }}
    QLineEdit#Search {{ border: 1px solid {_rgba(t.text, 0.14)}; border-radius: 9px; padding: 6px 10px;
                        background: {_rgba(t.text, 0.03)}; color: {t.text.name()};
                        selection-background-color: {t.accent.name()}; font-size: 13px; }}
    QLineEdit#Search:focus {{ border: 1px solid {t.accent.name()}; }}
    QPushButton#Primary {{ background: {t.accent.name()}; color: white; border: none; border-radius: 8px;
                           padding: 6px 14px; font-weight: 600; }}
    QPushButton#Primary:hover {{ background: {t.accent.lighter(110).name()}; }}
    QPushButton#Primary:disabled {{ background: {_rgba(t.accent, 0.35)}; }}
    QSlider#Timeline::groove:horizontal {{ height: 4px; border-radius: 2px;
                                          background: {_rgba(t.text, 0.12)}; }}
    QSlider#Timeline::sub-page:horizontal {{ border-radius: 2px; background: {t.accent.name()}; }}
    QSlider#Timeline::handle:horizontal {{ width: 14px; margin: -6px 0; border-radius: 7px;
                                          background: {t.panel.name()};
                                          border: 1px solid {_rgba(t.text, 0.25)}; }}
    QListWidget#Palette {{ border: none; background: transparent; outline: none; }}
    QListWidget#Modes {{ border: none; background: transparent; outline: none; color: {t.text.name()};
                         font-family: Menlo, Consolas, monospace; font-size: 11px; }}
    QListWidget#Modes::item {{ padding: 3px 6px; border-radius: 6px; }}
    QListWidget#Modes::item:hover {{ background: {_rgba(t.text, 0.06)}; }}
    QListWidget#Modes::item:selected {{ background: {_rgba(t.accent, 0.18)}; color: {t.text.name()}; }}
    QTabBar {{ background: {t.background.name()}; }}
    QTabBar::tab {{ background: {_rgba(t.text, 0.04)}; color: {t.muted.name()}; padding: 6px 18px;
                    border: none; border-right: 1px solid {_rgba(t.text, 0.08)}; }}
    QTabBar::tab:selected {{ background: {t.background.name()}; color: {t.text.name()}; }}
    """


def dialog_stylesheet(t: Theme) -> str:
    """Stylesheet for dialogs drawn in the app's theme rather than the system's.

    Plain buttons, number fields and text fields get theme colours too, so
    they stay readable when the app is dark and the system light, or the
    other way round.
    """
    field = _rgba(t.text, 0.18)
    return stylesheet(t) + f"""
    QDialog {{ background: {t.panel.name()}; }}
    QDialog QLabel, QDialog QCheckBox {{ color: {t.text.name()}; }}
    QDialog QPushButton {{ border: 1px solid {field}; border-radius: 8px; padding: 5px 16px;
                           background: {t.panel.name()}; color: {t.text.name()}; }}
    QDialog QPushButton:hover {{ background: {_rgba(t.text, 0.06)}; }}
    QDialog QCheckBox {{ spacing: 8px; }}
    QDialog QCheckBox::indicator, QDialog QAbstractItemView::indicator {{
        width: 13px; height: 13px; border: 1px solid {_rgba(t.text, 0.35)};
        border-radius: 4px; background: {t.panel.name()}; }}
    QDialog QCheckBox::indicator:checked, QDialog QAbstractItemView::indicator:checked {{
        border: 1px solid {t.accent.name()}; background: {t.accent.name()}; image: url({CHECK_IMAGE}); }}
    QDialog QAbstractSpinBox, QDialog QPlainTextEdit, QDialog QLineEdit#Field {{
        border: 1px solid {field}; border-radius: 6px; padding: 3px 6px; background: {t.panel.name()};
        color: {t.text.name()}; selection-background-color: {t.accent.name()}; }}
    QDialog QTableWidget {{ background: {t.panel.name()}; alternate-background-color: {_rgba(t.text, 0.04)};
                            color: {t.text.name()}; border: 1px solid {field}; border-radius: 6px;
                            gridline-color: {field}; selection-background-color: {_rgba(t.accent, 0.28)};
                            selection-color: {t.text.name()}; }}
    QDialog QHeaderView::section {{ background: {t.panel.name()}; color: {t.muted.name()}; border: none;
                                    border-bottom: 1px solid {field}; padding: 4px 6px; }}
    QDialog QTableCornerButton::section {{ background: {t.panel.name()}; border: none; }}
    QDialog QTextBrowser {{ background: {t.panel.name()}; color: {t.text.name()}; border: none; }}
    """


def theme_palette(widget: QWidget, t: Theme) -> None:
    """Give a natively drawn widget, such as a combo box, the theme's text and background colours.

    A stylesheet would replace the native look; the palette keeps it and
    only changes the colours, so the widget stays readable when the app and
    system themes differ.
    """
    palette = widget.palette()
    roles = QPalette.ColorRole
    for role in (roles.ButtonText, roles.Text, roles.WindowText):
        palette.setColor(role, t.text)
    for role in (roles.Button, roles.Base, roles.Window):
        palette.setColor(role, t.panel)
    widget.setPalette(palette)
    if isinstance(widget, QComboBox):
        widget.view().setPalette(palette)   # the open list is a widget of its own


def rule(vertical: bool = False) -> QFrame:
    """Thin separator line."""
    line = QFrame()
    line.setObjectName("Rule")
    if vertical:
        line.setFixedWidth(1)
    else:
        line.setFixedHeight(1)
    return line


def label(text: str = "", name: str = "") -> QLabel:
    """Label with a stylesheet name such as ``"Muted"`` or ``"Title"``."""
    widget = QLabel(text)
    if name:
        widget.setObjectName(name)
    return widget


def button(name: str, text: str = "", tip: str = "", checkable: bool = False) -> QToolButton:
    """Flat tool button styled by its stylesheet name."""
    b = QToolButton()
    b.setObjectName(name)
    b.setText(text)
    b.setToolTip(tip)
    b.setCheckable(checkable)
    b.setCursor(Qt.CursorShape.PointingHandCursor)
    return b


def _wrap(value: float, lo: float, hi: float) -> float:
    return (value - lo) % (hi - lo) + lo


class Segmented(QWidget):
    """Row of mutually exclusive text buttons."""

    changed = Signal(int)

    def __init__(self, labels: list[str], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(2)
        self.group = QButtonGroup(self)
        self.group.setExclusive(True)
        for n, text in enumerate(labels):
            b = button("Seg", text, checkable=True)
            self.group.addButton(b, n)
            lay.addWidget(b)
        self.set_index(0)
        self.group.idClicked.connect(self.changed.emit)

    def add(self, text: str) -> None:
        """Add a button at the end."""
        b = button("Seg", text, checkable=True)
        self.group.addButton(b, len(self.group.buttons()))
        layout = self.layout()
        if layout is not None:
            layout.addWidget(b)

    def remove_last(self) -> None:
        """Remove the last button."""
        buttons = self.group.buttons()
        if buttons:
            last = max(buttons, key=self.group.id)
            self.group.removeButton(last)
            last.deleteLater()

    def index(self) -> int:
        """Index of the checked button."""
        return self.group.checkedId()

    def set_index(self, n: int) -> None:
        """Check button ``n`` without emitting :attr:`changed`."""
        b = self.group.button(n)
        if b is not None:
            b.setChecked(True)


class ScrubField(QLineEdit):
    """Number field: drag sideways to scrub, click to type, arrows and wheel to step."""

    edited = Signal(float)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Scrub")
        self.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.setFixedWidth(104)
        self.setCursor(Qt.CursorShape.SizeHorCursor)
        self.value = 0.0
        self.step = 0.01
        self.lo, self.hi = 0.0, 1.0
        self.decimals = 3
        self.wrap = False
        self._press: QPointF | None = None
        self._last = QPointF()
        self._scrubbed = False
        self.returnPressed.connect(self._commit)
        self.editingFinished.connect(self._show)

    def configure(self, lo: float, hi: float, step: float, decimals: int, wrap: bool = False) -> None:
        """Set the range, the change per dragged pixel and the number of decimals shown."""
        self.lo, self.hi, self.step, self.decimals, self.wrap = lo, hi, step, decimals, wrap

    def set_value(self, value: float) -> None:
        """Show a value without emitting :attr:`edited` (unless the user is typing)."""
        self.value = value
        if not self.hasFocus():
            self._show()

    def _show(self) -> None:
        self.setText(f"{self.value:.{self.decimals}f}")

    def _emit(self, value: float) -> None:
        self.value = _wrap(value, self.lo, self.hi) if self.wrap else max(self.lo, min(self.hi, value))
        self._show()
        self.edited.emit(self.value)

    def _commit(self) -> None:
        try:
            value = float(self.text().replace(",", "."))
        except ValueError:
            self._show()
            return
        self._emit(value)
        self.selectAll()

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if self.hasFocus():
            super().mousePressEvent(e)
            return
        self._press = self._last = e.position()
        self._scrubbed = False

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        if self._press is None:
            super().mouseMoveEvent(e)
            return
        if not self._scrubbed and abs(e.position().x() - self._press.x()) < 3:
            return
        dx = e.position().x() - self._last.x()
        self._scrubbed = True
        self._last = e.position()
        fine = FINE if e.modifiers() & Qt.KeyboardModifier.AltModifier else 1.0
        self._emit(self.value + dx * self.step * fine)

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        if self._press is None:
            super().mouseReleaseEvent(e)
            return
        self._press = None
        if not self._scrubbed:
            self.setFocus()
            self.selectAll()

    @override
    def keyPressEvent(self, e: QKeyEvent) -> None:
        if e.key() in (Qt.Key.Key_Up, Qt.Key.Key_Down):
            big = 10 if e.modifiers() & Qt.KeyboardModifier.ShiftModifier else 1
            sign = 1 if e.key() == Qt.Key.Key_Up else -1
            self._emit(self.value + sign * self.step * 10 * big)
            self.selectAll()
        elif e.key() == Qt.Key.Key_Escape:
            self._show()
            self.clearFocus()
        else:
            super().keyPressEvent(e)

    @override
    def wheelEvent(self, e: QWheelEvent) -> None:
        steps = e.angleDelta().y() / 120.0 or e.pixelDelta().y() / 30.0
        self._emit(self.value + steps * self.step * 10)


class Dial(QWidget):
    """Horizontal ruler that is dragged or scrolled to change a value, like a tape measure."""

    changed = Signal(float)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(34)
        self.setMinimumWidth(380)
        self.setCursor(Qt.CursorShape.SizeHorCursor)
        self.theme: Theme | None = None
        self.value = 0.0
        self.lo, self.hi = 0.0, 1.0
        self.per_px = 0.01
        self.minor, self.major = 0.1, 1.0
        self.decimals = 1
        self.wrap = False
        self._last: float | None = None

    def configure(self, lo: float, hi: float, per_px: float, minor: float, major: float, decimals: int,
                  wrap: bool = False) -> None:
        """Set the range, value per pixel, tick spacing and label decimals."""
        self.lo, self.hi, self.per_px = lo, hi, per_px
        self.minor, self.major, self.decimals, self.wrap = minor, major, decimals, wrap
        self.update()

    def set_value(self, value: float) -> None:
        """Show a value without emitting :attr:`changed`."""
        self.value = value
        self.update()

    def step(self, ticks: float) -> None:
        """Move by a number of minor ticks."""
        self._move(self.value + ticks * self.minor)

    def _move(self, value: float) -> None:
        value = _wrap(value, self.lo, self.hi) if self.wrap else max(self.lo, min(self.hi, value))
        if value != self.value:
            self.value = value
            self.update()
            self.changed.emit(value)

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        self._last = e.position().x()

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        if self._last is None:
            return
        dx = e.position().x() - self._last
        self._last = e.position().x()
        fine = FINE if e.modifiers() & Qt.KeyboardModifier.AltModifier else 1.0
        self._move(self.value - dx * self.per_px * fine)

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        self._last = None

    @override
    def wheelEvent(self, e: QWheelEvent) -> None:
        if not e.pixelDelta().isNull():
            d = e.pixelDelta().x() or -e.pixelDelta().y()
            self._move(self.value - d * self.per_px)
        else:
            self.step((e.angleDelta().y() or e.angleDelta().x()) / 120.0)
        e.accept()

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        if self.theme is None:
            return
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        w = self.width()
        cx = w / 2
        font = p.font()
        font.setPixelSize(9)
        p.setFont(font)
        span = cx * self.per_px
        first = math.floor((self.value - span) / self.minor)
        last = math.ceil((self.value + span) / self.minor)
        for k in range(first, last + 1):
            v = k * self.minor
            if not self.wrap and not self.lo - 1e-9 <= v <= self.hi + 1e-9:
                continue
            x = cx + (v - self.value) / self.per_px
            fade = max(0.0, 1.0 - abs(x - cx) / cx) ** 0.7
            major = abs(v / self.major - round(v / self.major)) < 1e-6
            color = QColor(self.theme.text)
            color.setAlphaF((0.55 if major else 0.25) * fade)
            p.setPen(color)
            p.drawLine(int(x), 4 if major else 8, int(x), 14)
            if major and 20 <= x <= w - 20:   # a label cut by the edge is left out
                shown = _wrap(v, self.lo, self.hi) if self.wrap else v
                text = f"{shown:.{self.decimals}f}"
                p.drawText(int(x) - 20, 17, 40, 14, Qt.AlignmentFlag.AlignCenter, text)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(self.theme.accent)
        p.drawRoundedRect(int(cx) - 1, 1, 3, 16, 1.5, 1.5)
        p.end()
