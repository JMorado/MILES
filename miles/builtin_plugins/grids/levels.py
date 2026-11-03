"""The orbital energy diagram of the grid panel: levels to click, zoom and pan, with the HOMO-LUMO gap."""

import math
from dataclasses import dataclass
from typing import override

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QFontMetricsF,
    QMouseEvent,
    QPainter,
    QPaintEvent,
    QPen,
    QWheelEvent,
)
from PySide6.QtWidgets import QToolTip, QWidget

from miles import HARTREE_EV, OrbitalSet, Spin
from miles.plugin import (
    Theme,
)

DEGENERATE_EV = 0.01   # levels closer than this are drawn stacked, as one degenerate set
MAX_STACK = 5          # at most this many ticks share one place; zoom in to separate more
FRONTIER = 8           # the diagram first shows this many levels below the HOMO and above the LUMO
FRONTIER_EV = 20.0     # ... but no further than this from them, so core levels do not squeeze the view
MIN_SPAN_EV = 0.3      # narrowest energy window when zooming in


@dataclass
class _Level:
    """One orbital as drawn in the energy diagram."""

    index: int
    energy: float       # eV
    filled: float       # fraction of the full occupation, 0 to 1
    row: int
    slot: int = 0       # position within a stack of levels too close to draw side by side
    slots: int = 1      # size of that stack


class OrbitalLevels(QWidget):
    """Orbital energies as a strip of ticks, energy increasing to the right.

    Occupied orbitals are filled, virtual ones drawn as outlines, and levels
    too close to tell apart at the current zoom (degenerate ones at any zoom)
    are stacked. Unrestricted calculations get an alpha row
    above a beta row. The wheel zooms, dragging pans, a double click returns
    to the frontier orbitals, and clicking a tick emits :attr:`orbitalClicked`.
    """

    orbitalClicked = Signal(int)

    ROW = 28        # label band and tick band of one row, in pixels
    LABELS = 11     # height of the label band above the ticks
    AXIS = 13       # height of the energy axis below the rows
    TICK = 5        # width of a tick
    PICK = 6        # how far from a tick a click still counts, in pixels

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumWidth(380)
        self.setMouseTracking(True)
        self.theme: Theme | None = None
        self.orbitals: OrbitalSet | None = None
        self.current: int | None = None
        self.levels: list[_Level] = []
        self.rows: list[Spin] = []
        self.lo, self.hi = -10.0, 5.0
        self._limits = (-10.0, 5.0)
        self._press: float | None = None
        self._dragged = False
        self._hover: int | None = None
        self._layout_key: tuple[float, float, int, int] | None = None
        self.setFixedHeight(2 + self.ROW + self.AXIS)

    def set_orbitals(self, orbitals: OrbitalSet | None) -> None:
        """Show these orbitals, starting from the frontier region (no change for the same set)."""
        if orbitals is self.orbitals:
            return
        self.orbitals = orbitals
        self.levels = []
        if orbitals is None:
            self.update()
            return
        unrestricted = orbitals.unrestricted
        self.rows = [Spin.ALPHA, Spin.BETA] if unrestricted else [Spin.NONE]
        full = 1.0 if unrestricted else 2.0
        for row, spin in enumerate(self.rows):
            members = [k for k, o in enumerate(orbitals.orbitals)
                       if not unrestricted or (o.spin or Spin.ALPHA) == spin]
            members.sort(key=lambda k: orbitals.orbitals[k].energy)
            for k in members:
                o = orbitals.orbitals[k]
                filled = min(1.0, max(0.0, o.occupation / full))
                self.levels.append(_Level(k, o.energy * HARTREE_EV, filled, row))
        self._layout_key = None
        self.setFixedHeight(2 + self.ROW * len(self.rows) + self.AXIS)
        energies = [lv.energy for lv in self.levels]
        margin = max(0.5, 0.04 * (max(energies) - min(energies)))
        self._limits = (min(energies) - margin, max(energies) + margin)
        self.reset_view()

    def _ensure_layout(self) -> None:
        """Stack the ticks that would overlap at the current zoom; redone only when the view changes."""
        key = (self.lo, self.hi, self.width(), len(self.levels))
        if key == self._layout_key:
            return
        self._layout_key = key
        for row in range(len(self.rows)):
            group: list[_Level] = []
            for level in (lv for lv in self.levels if lv.row == row):   # sorted by energy within a row
                if group and (len(group) == MAX_STACK or not self._touching(group[-1], level)):
                    self._stack(group)
                    group = []
                group.append(level)
            self._stack(group)

    def _touching(self, a: _Level, b: _Level) -> bool:
        """Whether two neighbouring levels are degenerate or their ticks would overlap on screen."""
        close = b.energy - a.energy <= DEGENERATE_EV
        return close or self._x(b.energy) - self._x(a.energy) < self.TICK + 1

    @staticmethod
    def _stack(group: list[_Level]) -> None:
        for slot, level in enumerate(group):
            level.slot, level.slots = slot, len(group)

    def frontier(self, row: int) -> tuple[_Level | None, _Level | None]:
        """Highest occupied and lowest virtual level of a row."""
        levels = [lv for lv in self.levels if lv.row == row]
        occupied = [lv for lv in levels if lv.filled * 2 > OrbitalSet.OCCUPIED]
        virtual = [lv for lv in levels if lv.filled * 2 <= OrbitalSet.OCCUPIED]
        homo = max(occupied, key=lambda lv: (lv.energy, lv.index)) if occupied else None
        lumo = min(virtual, key=lambda lv: (lv.energy, lv.index)) if virtual else None
        return homo, lumo

    def reset_view(self) -> None:
        """Show about :data:`FRONTIER` levels on each side of the HOMO-LUMO gap."""
        if not self.levels:
            return
        energies = sorted(lv.energy for lv in self.levels if lv.row == 0)
        homo, lumo = self.frontier(0)
        below = [e for e in energies if homo is not None and e <= homo.energy]
        above = [e for e in energies if lumo is not None and e >= lumo.energy]
        lo = below[-FRONTIER - 1] if len(below) > FRONTIER else energies[0]
        hi = above[FRONTIER] if len(above) > FRONTIER else energies[-1]
        if homo is not None:
            lo = max(lo, homo.energy - FRONTIER_EV)
        if lumo is not None:
            hi = min(hi, lumo.energy + FRONTIER_EV)
        margin = max(0.4, 0.06 * (hi - lo))
        self.lo, self.hi = lo - margin, hi + margin
        self._clamp()
        self.update()

    def set_current(self, index: int | None) -> None:
        """Highlight the orbital on screen."""
        if index != self.current:
            self.current = index
            self.update()

    def _plot(self) -> tuple[float, float]:
        gutter = 18.0 if len(self.rows) > 1 else 0.0
        return gutter + 6.0, self.width() - 6.0

    def _x(self, energy: float) -> float:
        left, right = self._plot()
        return left + (energy - self.lo) / (self.hi - self.lo) * (right - left)

    def _energy(self, x: float) -> float:
        left, right = self._plot()
        return self.lo + (x - left) / max(1.0, right - left) * (self.hi - self.lo)

    def _tick_rect(self, level: _Level) -> QRectF:
        self._ensure_layout()
        top = 2 + level.row * self.ROW + self.LABELS
        band = self.ROW - self.LABELS - 2
        height = max(2.0, (band - (level.slots - 1)) / level.slots)
        x = self._x(level.energy)
        return QRectF(x - self.TICK / 2, top + level.slot * (height + 1), self.TICK, height)

    def tick_position(self, index: int) -> QPointF | None:
        """Centre of an orbital's tick in widget coordinates, or None when it is off the strip."""
        left, right = self._plot()
        for level in self.levels:
            if level.index == index:
                rect = self._tick_rect(level)
                return rect.center() if left <= rect.center().x() <= right else None
        return None

    def level_at(self, pos: QPointF) -> int | None:
        """Orbital whose tick is under (or next to) a point."""
        row = int((pos.y() - 2) // self.ROW)
        left, right = self._plot()
        if not 0 <= row < len(self.rows) or not left - self.PICK <= pos.x() <= right + self.PICK:
            return None
        best: tuple[float, float, int] | None = None
        for level in self.levels:
            if level.row != row:
                continue
            rect = self._tick_rect(level)
            dx = abs(rect.center().x() - pos.x())
            if dx > self.PICK:
                continue
            dy = 0.0 if rect.top() <= pos.y() <= rect.bottom() else abs(rect.center().y() - pos.y())
            if best is None or (dx, dy) < best[:2]:
                best = (dx, dy, level.index)
        return best[2] if best is not None else None

    def tooltip(self, index: int) -> str:
        """Text shown when hovering over an orbital."""
        assert self.orbitals is not None
        orbital = self.orbitals.orbitals[index]
        tag = self.orbitals.label(index)
        head = f"MO {index + 1}" + (f"  {tag}" if tag else "")
        return f"{head}\n{orbital.energy * HARTREE_EV:.2f} eV, occupation {orbital.occupation:g}"

    def _clamp(self) -> None:
        low, high = self._limits
        span = min(max(self.hi - self.lo, MIN_SPAN_EV), high - low)
        lo = min(max(self.lo, low), high - span)
        self.lo, self.hi = lo, lo + span

    def _zoom(self, steps: float, x: float) -> None:
        anchor = self._energy(x)
        factor = 0.85 ** steps
        self.lo = anchor - (anchor - self.lo) * factor
        self.hi = anchor + (self.hi - anchor) * factor
        self._clamp()
        self.update()

    def _pan(self, dx: float) -> None:
        left, right = self._plot()
        shift = dx / max(1.0, right - left) * (self.hi - self.lo)
        self.lo -= shift
        self.hi -= shift
        self._clamp()
        self.update()

    @override
    def wheelEvent(self, e: QWheelEvent) -> None:
        pixels = e.pixelDelta()
        if not pixels.isNull() and abs(pixels.x()) > abs(pixels.y()):
            self._pan(pixels.x())
        elif not pixels.isNull():
            self._zoom(pixels.y() / 40.0, e.position().x())
        else:
            self._zoom(e.angleDelta().y() / 120.0, e.position().x())
        e.accept()

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if e.button() == Qt.MouseButton.LeftButton:
            self._press = e.position().x()
            self._dragged = False

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        x = e.position().x()
        if self._press is not None and e.buttons() & Qt.MouseButton.LeftButton:
            if self._dragged or abs(x - self._press) > 3:
                self._dragged = True
                self.setCursor(Qt.CursorShape.ClosedHandCursor)
                self._pan(x - self._press)
                self._press = x
            return
        hit = self.level_at(e.position())
        self.setCursor(Qt.CursorShape.PointingHandCursor if hit is not None else Qt.CursorShape.ArrowCursor)
        if hit != self._hover:
            self._hover = hit
            if hit is None:
                QToolTip.hideText()
            else:
                QToolTip.showText(e.globalPosition().toPoint(), self.tooltip(hit), self)

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        if e.button() != Qt.MouseButton.LeftButton or self._press is None:
            return
        dragged = self._dragged
        self._press = None
        self._dragged = False
        self.unsetCursor()
        if not dragged:
            hit = self.level_at(e.position())
            if hit is not None:
                self.orbitalClicked.emit(hit)

    @override
    def mouseDoubleClickEvent(self, e: QMouseEvent) -> None:
        if self.level_at(e.position()) is None:
            self.reset_view()

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        if self.theme is None or not self.levels:
            return
        theme = self.theme
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        font = QFont(self.font())
        font.setPixelSize(9)
        p.setFont(font)
        left, right = self._plot()
        faint = QColor(theme.text)
        faint.setAlphaF(0.12)
        ink = QColor(theme.text)
        ink.setAlphaF(0.72)
        outline = QColor(theme.muted)
        axis_y = 2 + self.ROW * len(self.rows)
        p.setPen(QPen(faint, 1))
        p.drawLine(QPointF(left, axis_y), QPointF(right, axis_y))
        self._paint_axis(p, left, right, axis_y)
        if len(self.rows) > 1:
            p.setPen(theme.muted)
            for row, symbol in enumerate(("\u03b1", "\u03b2")):
                top = 2 + row * self.ROW + self.LABELS
                cell = QRectF(0, top, 18, self.ROW - self.LABELS - 2)
                p.drawText(cell, Qt.AlignmentFlag.AlignCenter, symbol)
        p.setClipRect(QRectF(left - self.TICK, 0, right - left + 2 * self.TICK, self.height()))
        for level in self.levels:
            rect = self._tick_rect(level)
            if rect.right() < left - self.TICK or rect.left() > right + self.TICK:
                continue
            current = level.index == self.current
            color = QColor(theme.accent) if current else ink
            p.setPen(QPen(color if current or level.filled > 0 else outline, 1))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(rect.adjusted(0.5, 0.5, -0.5, -0.5))
            if level.filled > 0:
                filled = QRectF(rect)
                filled.setTop(rect.bottom() - rect.height() * level.filled)
                p.fillRect(filled, color)
        p.setClipping(False)
        for row in range(len(self.rows)):
            self._paint_frontier(p, row)
        p.end()

    def _paint_axis(self, p: QPainter, left: float, right: float, y: float) -> None:
        assert self.theme is not None
        span = self.hi - self.lo
        raw = span / 6
        step = 10 ** math.floor(math.log10(raw))
        step *= next(m for m in (1, 2, 5, 10) if raw <= m * step)
        p.setPen(self.theme.muted)
        metrics = QFontMetricsF(p.font())
        room = right - metrics.horizontalAdvance("eV") - 6   # labels stop short of the unit
        value = math.ceil(self.lo / step) * step
        while value <= self.hi:
            x = self._x(value)
            if left <= x <= right:
                p.drawLine(QPointF(x, y), QPointF(x, y + 3))
                text = f"{value:g}"
                half = metrics.horizontalAdvance(text) / 2
                if x - half >= 0 and x + half <= room:   # a label cut by the edge is left out
                    p.drawText(QRectF(x - 24, y + 2, 48, self.AXIS - 2), Qt.AlignmentFlag.AlignCenter, text)
            value += step
        p.drawText(QRectF(right - 40, y + 2, 40, self.AXIS - 2),
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, "eV")

    def _paint_frontier(self, p: QPainter, row: int) -> None:
        """HOMO and LUMO labels centred above their ticks, and the gap between them."""
        assert self.theme is not None
        homo, lumo = self.frontier(row)
        left, right = self._plot()
        top = 2 + row * self.ROW
        metrics = p.fontMetrics()
        labels: list[tuple[str, float, float]] = []   # text, left edge, width
        for name, level in (("HOMO", homo), ("LUMO", lumo)):
            if level is not None and left <= self._x(level.energy) <= right:
                width = float(metrics.horizontalAdvance(name) + 2)
                labels.append((name, self._x(level.energy) - width / 2, width))
        if len(labels) == 2:   # push them apart if they collide
            (a, xa, wa), (b, xb, wb) = labels
            overlap = xa + wa + 4 - xb
            if overlap > 0:
                labels = [(a, xa - overlap / 2, wa), (b, xb + overlap / 2, wb)]
        p.setPen(self.theme.muted)
        bottom = Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignBottom
        for name, x, width in labels:
            p.drawText(QRectF(x, top, width, self.LABELS), bottom, name)
        if homo is None or lumo is None:
            return
        text = f"gap {lumo.energy - homo.energy:.2f} eV"
        width = metrics.horizontalAdvance(text) + 8
        start = (self._x(homo.energy) + self._x(lumo.energy)) / 2 - width / 2
        clash = any(x < start + width and start < x + w for _, x, w in labels)
        if clash or start < left or start + width > right:
            end = max([x + w for _, x, w in labels] + [left])   # no room between them: after the labels
            start = end + 6 if end + 6 + width <= right else right - width
        p.drawText(QRectF(start, top, width, self.LABELS), bottom, text)

