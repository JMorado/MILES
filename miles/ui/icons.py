"""Thin line icons drawn with QPainter in the theme's text colour."""

import math
from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QIcon, QPainter, QPainterPath, QPen, QPixmap, QPolygonF
from PySide6.QtSvg import QSvgRenderer

GRID = 24.0
"""Icons are drawn on a 24 x 24 grid and scaled to the pixmap size."""


def _canvas(size: int) -> tuple[QPixmap, QPainter]:
    pm = QPixmap(size, size)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setRenderHint(QPainter.RenderHint.TextAntialiasing)
    p.scale(size / GRID, size / GRID)
    return pm, p


def _pen(color: QColor, width: float = 1.6) -> QPen:
    pen = QPen(color, width)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    return pen


def _hexagon(cx: float, cy: float, r: float, start: float = 90) -> QPolygonF:
    return QPolygonF([QPointF(cx + r * math.cos(math.radians(start + 60 * k)),
                              cy + r * math.sin(math.radians(start + 60 * k))) for k in range(6)])


def _sparkle(p: QPainter, x: float, y: float, r: float) -> None:
    path = QPainterPath()
    path.moveTo(x, y - r)
    path.quadTo(x + r * 0.12, y - r * 0.12, x + r, y)
    path.quadTo(x + r * 0.12, y + r * 0.12, x, y + r)
    path.quadTo(x - r * 0.12, y + r * 0.12, x - r, y)
    path.quadTo(x - r * 0.12, y - r * 0.12, x, y - r)
    p.drawPath(path)


def _draw_pointer(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c, 1.6))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPolygon(QPolygonF([QPointF(6.5, 3.5), QPointF(6.5, 18.5), QPointF(10.3, 14.9),
                             QPointF(13.0, 20.6), QPointF(15.6, 19.4), QPointF(12.9, 13.8),
                             QPointF(18.2, 13.8)]))


def _draw_delete(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c))
    p.drawLine(QPointF(4, 6.5), QPointF(20, 6.5))
    p.drawLine(QPointF(9.5, 6.5), QPointF(10, 4))
    p.drawLine(QPointF(10, 4), QPointF(14, 4))
    p.drawLine(QPointF(14, 4), QPointF(14.5, 6.5))
    path = QPainterPath()
    path.moveTo(6, 6.5)
    path.lineTo(7, 19.5)
    path.quadTo(7.2, 20.5, 8.2, 20.5)
    path.lineTo(15.8, 20.5)
    path.quadTo(16.8, 20.5, 17, 19.5)
    path.lineTo(18, 6.5)
    p.drawPath(path)
    p.drawLine(QPointF(10.2, 10), QPointF(10.4, 17))
    p.drawLine(QPointF(13.8, 10), QPointF(13.6, 17))


def _draw_clean(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c, 1.5))
    p.setBrush(Qt.BrushStyle.NoBrush)
    _sparkle(p, 10, 13, 7.5)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(c)
    _sparkle(p, 18.5, 5.5, 3.2)
    _sparkle(p, 19, 18.5, 2.2)


def _draw_symmetry(p: QPainter, c: QColor) -> None:
    p.setPen(QPen(c, 1.4, Qt.PenStyle.DashLine, Qt.PenCapStyle.RoundCap))
    p.drawLine(QPointF(12, 2.5), QPointF(12, 21.5))
    p.setPen(_pen(c))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPolygon(QPolygonF([QPointF(9.5, 6), QPointF(9.5, 18), QPointF(3, 15)]))
    p.drawPolygon(QPolygonF([QPointF(14.5, 6), QPointF(14.5, 18), QPointF(21, 15)]))


def _draw_ring(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPolygon(_hexagon(12, 12, 8.5))
    p.drawEllipse(QPointF(12, 12), 4.2, 4.2)


def _draw_group(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c))
    p.drawLine(QPointF(4, 18), QPointF(10, 13))
    p.drawLine(QPointF(10, 13), QPointF(10, 5.5))
    p.drawLine(QPointF(10, 13), QPointF(17, 16.5))
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(c)
    for x, y, r in ((4, 18, 2.0), (10, 13, 2.3), (10, 5.5, 2.0), (17, 16.5, 2.0)):
        p.drawEllipse(QPointF(x, y), r, r)


def _draw_atom(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawEllipse(QPointF(12, 12), 7.5, 7.5)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(c)
    p.drawEllipse(QPointF(12, 12), 2.2, 2.2)


def _draw_play(p: QPainter, c: QColor) -> None:
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(c)
    p.drawPolygon(QPolygonF([QPointF(8, 5.5), QPointF(18.5, 12), QPointF(8, 18.5)]))


def _draw_pause(p: QPainter, c: QColor) -> None:
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(c)
    p.drawRoundedRect(QRectF(7, 5.5, 3.6, 13), 1.2, 1.2)
    p.drawRoundedRect(QRectF(13.4, 5.5, 3.6, 13), 1.2, 1.2)


def _draw_pin(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c, 1.5))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawLine(QPointF(12, 14.5), QPointF(12, 21))
    path = QPainterPath()
    path.moveTo(8.5, 3.5)
    path.lineTo(15.5, 3.5)
    path.lineTo(14.5, 9)
    path.lineTo(17.5, 14.5)
    path.lineTo(6.5, 14.5)
    path.lineTo(9.5, 9)
    path.closeSubpath()
    p.drawPath(path)


def _draw_plot(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c, 1.5))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPolyline(QPolygonF([QPointF(4, 4), QPointF(4, 20), QPointF(20, 20)]))
    p.drawPolyline(QPolygonF([QPointF(7, 15), QPointF(10.5, 9), QPointF(14, 13), QPointF(19, 6)]))


def _draw_collapse(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c, 1.8))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPolyline(QPolygonF([QPointF(7, 10), QPointF(12, 15), QPointF(17, 10)]))   # pointing down: open


def _draw_expand(p: QPainter, c: QColor) -> None:
    p.setPen(_pen(c, 1.8))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPolyline(QPolygonF([QPointF(10, 7), QPointF(15, 12), QPointF(10, 17)]))   # pointing right: folded


_DRAWERS: dict[str, Callable[[QPainter, QColor], None]] = {
    "pointer": _draw_pointer,
    "delete": _draw_delete,
    "clean": _draw_clean,
    "symmetry": _draw_symmetry,
    "ring": _draw_ring,
    "group": _draw_group,
    "atom": _draw_atom,
    "play": _draw_play,
    "pause": _draw_pause,
    "pin": _draw_pin,
    "plot": _draw_plot,
    "collapse": _draw_collapse,
    "expand": _draw_expand,
}


def icon(name: str, color: QColor, size: int = 48) -> QIcon:
    """Draw a line icon by name (see ``_DRAWERS``: pointer, delete, play, collapse, expand and so on)."""
    pm, p = _canvas(size)
    _DRAWERS[name](p, color)
    p.end()
    return QIcon(pm)


def element_chip(symbol: str, fill: QColor, text: QColor, size: int = 48) -> QIcon:
    """Draw a filled circle with an element symbol, for the current element."""
    pm, p = _canvas(size)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(fill)
    p.drawEllipse(QPointF(12, 12), 10.5, 10.5)
    font = QFont()
    font.setPixelSize(11 if len(symbol) == 1 else 9)
    font.setWeight(QFont.Weight.DemiBold)
    p.setFont(font)
    p.setPen(text)
    p.drawText(QRectF(0, 0.3, GRID, GRID), Qt.AlignmentFlag.AlignCenter, symbol)
    p.end()
    return QIcon(pm)


def bond_order_icon(order: float, color: QColor, size: int = 40) -> QIcon:
    """Lines standing for a bond order (dotted for no bond, dashed second line for aromatic)."""
    pm, p = _canvas(size)
    if order == 0:
        p.setPen(QPen(color, 1.5, Qt.PenStyle.DotLine, Qt.PenCapStyle.RoundCap))
        p.drawLine(QPointF(4, 12), QPointF(20, 12))
    elif order == 1.5:
        p.setPen(_pen(color, 1.7))
        p.drawLine(QPointF(4, 9.5), QPointF(20, 9.5))
        p.setPen(QPen(color, 1.7, Qt.PenStyle.DashLine, Qt.PenCapStyle.FlatCap))
        p.drawLine(QPointF(4, 14.5), QPointF(20, 14.5))
    else:
        p.setPen(_pen(color, 1.7))
        for y in {1: [12.0], 2: [9.5, 14.5], 3: [7.5, 12.0, 16.5]}[int(order)]:
            p.drawLine(QPointF(4, y), QPointF(20, y))
    p.end()
    return QIcon(pm)


APP_ICON = Path(__file__).parent / "images" / "miles.svg"
"""The MILES logo: a ball-and-stick M whose last stroke is a trumpet bell."""

APP_ICON_SIZES = (16, 24, 32, 48, 64, 128, 256, 512, 1024)


def app_icon() -> QIcon:
    """Application icon, drawn from :data:`APP_ICON` at every size a window, Dock or taskbar asks for."""
    renderer = QSvgRenderer(str(APP_ICON))
    result = QIcon()
    for size in APP_ICON_SIZES:
        pm = QPixmap(size, size)
        pm.fill(Qt.GlobalColor.transparent)
        p = QPainter(pm)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        renderer.render(p, QRectF(0, 0, size, size))
        p.end()
        result.addPixmap(pm)
    return result
