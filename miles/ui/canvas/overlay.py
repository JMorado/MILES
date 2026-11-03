"""Text and marks drawn over the molecule: labels, selection numbers, status line and toasts."""

import math
from itertools import pairwise
from typing import TYPE_CHECKING

import numpy as np
from PySide6.QtCore import QLineF, QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPainterPath, QPen, QPixmap

from miles.core.fragments import Category
from miles.core.geometry import angle, dihedral, distance
from miles.core.molecule import Molecule
from miles.render.camera import Camera
from miles.render.layers import Annotations, Arc
from miles.render.picking import AtomHit, BondHit
from miles.render.scene import LabelMode
from miles.render.theme import Theme
from miles.ui.builder import Tool
from miles.ui.canvas.tools import Drag
from miles.ui.text import ANGSTROM, DEGREE, SEPARATOR, describe_atom, describe_bond

if TYPE_CHECKING:
    from miles.ui.canvas.canvas import Canvas

MARGIN = 18
TOAST_WIDTH = 620.0   # longer messages wrap onto more lines

LABEL_LIMIT = 1000
"""Most atom labels drawn at once; with more atoms in view, those nearest the viewer are labelled."""

MOVING_LABEL_LIMIT = 200
"""While the view turns or moves, atom labels are left out if more than this many would be drawn."""

OUTLINED_LABELS = 200
"""Up to this many atom labels are drawn as outlines; more are copied from cached pictures of the texts."""

LABEL_MARGIN = 24.0   # pixels outside the window within which an atom still gets its label
VIEW_DRAGS = frozenset({Drag.ROTATE, Drag.SPIN, Drag.PAN, Drag.ZOOM, Drag.SCALE})
TEXT_CACHE_SIZE = 4096
LABEL_PICTURES = 1500   # cached label pictures, a little more than LABEL_LIMIT
HALO_PAD = 3.0          # pixels around a label picture for its halo
WHITE = QColor(255, 255, 255)


def paint_overlay(p: QPainter, canvas: "Canvas") -> None:
    """Draw everything that sits on top of the molecule."""
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    p.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    w, h = canvas.width(), canvas.height()
    mol = canvas.mol
    if len(mol):
        view = canvas.camera.to_view(mol.coords)
        sx, sy = canvas.camera.to_screen(view, w, h)
        if canvas.display.labels != LabelMode.NONE:
            _paint_labels(p, canvas, sx, sy, view[:, 2])
        paint_annotations(p, canvas, canvas.theme)
        if canvas.selection:
            _paint_selection_order(p, canvas, sx, sy)
        g = canvas.gesture
        if g is not None and g.drag == Drag.DRAW and g.moved and g.draw_from is not None \
                and g.draw_from < len(mol):
            _paint_rubber_band(p, canvas, sx, sy)
    else:
        _paint_empty_hint(p, canvas)

    p.setFont(_font(canvas, 11.5))
    p.setPen(canvas.theme.muted)
    if canvas.info_text:
        p.drawText(QPointF(MARGIN, h - 16), canvas.info_text)
    hover = _hover_text(canvas)
    if hover:
        fm = QFontMetricsF(p.font())
        hover = fm.elidedText(hover, Qt.TextElideMode.ElideLeft, w * 0.55)
        p.drawText(QPointF(w - MARGIN - fm.horizontalAdvance(hover), h - 16), hover)
    if canvas.toast_text:
        _paint_toast(p, canvas)


def _font(canvas: "Canvas", size: float, bold: bool = False) -> QFont:
    font = QFont(canvas.font())
    font.setPointSizeF(size)
    font.setBold(bold)
    return font


def _hover_text(canvas: "Canvas") -> str:
    hit = canvas.hover
    if isinstance(hit, BondHit):
        return describe_bond(canvas.mol, hit.i, hit.j)
    if isinstance(hit, AtomHit):
        return describe_atom(canvas.mol, hit.index)
    return ""


def paint_annotations(p: QPainter, canvas: "Canvas", theme: Theme,
                      annotations: list[Annotations] | None = None) -> None:
    """Draw the measurements, stereo labels and the lines and labels of plugins.

    Used on screen and for exported pictures. The plugins' lines go under
    the measurements and their labels on top of everything.

    Parameters
    ----------
    p : QPainter
        Painter on the canvas or the picture.
    canvas : Canvas
        Canvas whose view is drawn.
    theme : Theme
        Colours.
    annotations : list of Annotations, optional
        What plugins draw (default: the canvas's, for the screen).
    """
    mol, camera = canvas.mol, canvas.camera
    w, h = canvas.width(), canvas.height()
    layers = canvas.drawing.annotations(theme) if annotations is None else annotations
    for layer in layers:
        _paint_layer_lines(p, canvas, theme, layer)
    if mol.measurements:
        paint_measurements(p, mol, camera, w, h, theme, canvas.font(), canvas.doc.display.hidden)
    if canvas.display.stereo:
        _paint_stereo(p, canvas, theme)
    for layer in layers:
        _paint_layer_arcs(p, canvas, theme, layer)
    for layer in layers:
        _paint_layer_labels(p, canvas, theme, layer)


def _paint_layer_lines(p: QPainter, canvas: "Canvas", theme: Theme, layer: Annotations) -> None:
    """Lines of a plugin layer, with the pen each asks for."""
    lines = canvas.drawing.resolve_lines(layer)
    if not lines:
        return
    camera = canvas.camera
    ends = np.array([[a, b] for _, a, b in lines]).reshape(-1, 3)
    sx, sy = camera.to_screen(camera.to_view(ends), canvas.width(), canvas.height())
    p.setBrush(Qt.BrushStyle.NoBrush)
    current = None
    for k, (line, _, _) in enumerate(lines):
        style = (line.color, line.width, line.dashes, line.round_caps)
        if style != current:
            current = style
            cap = Qt.PenCapStyle.RoundCap if line.round_caps else Qt.PenCapStyle.FlatCap
            pen = QPen(QColor(line.color) if line.color else QColor(theme.text), line.width,
                       Qt.PenStyle.SolidLine, cap)
            if line.dashes:
                pen.setDashPattern(list(line.dashes))
            p.setPen(pen)
        p.drawLine(QLineF(sx[2 * k], sy[2 * k], sx[2 * k + 1], sy[2 * k + 1]))


def _arc_geometry(canvas: "Canvas", arc: Arc) -> tuple[float, float, float, float] | None:
    """Screen centre and the start and span angles (degrees, Qt's convention) of an arc, or None."""
    points = [canvas.drawing.anchor_point(a, arc.drawn_only) for a in (arc.center, arc.start, arc.end)]
    if any(point is None for point in points):
        return None
    camera = canvas.camera
    sx, sy = camera.to_screen(camera.to_view(np.array(points)), canvas.width(), canvas.height())
    first = math.degrees(math.atan2(-(sy[1] - sy[0]), sx[1] - sx[0]))   # Qt angles run anticlockwise
    second = math.degrees(math.atan2(-(sy[2] - sy[0]), sx[2] - sx[0]))
    span = (second - first + 180.0) % 360.0 - 180.0   # the shorter way round
    return float(sx[0]), float(sy[0]), first, span


def _paint_layer_arcs(p: QPainter, canvas: "Canvas", theme: Theme, layer: Annotations) -> None:
    """Arcs of a plugin layer, centred on their vertex."""
    for arc in layer.arcs:
        geometry = _arc_geometry(canvas, arc)
        if geometry is None:
            continue
        x, y, start, span = geometry
        pen = QPen(QColor(arc.color) if arc.color else QColor(theme.text), arc.width, Qt.PenStyle.SolidLine,
                   Qt.PenCapStyle.RoundCap)
        if arc.dashes:
            pen.setDashPattern(list(arc.dashes))
        p.setPen(pen)
        p.setBrush(Qt.BrushStyle.NoBrush)
        r = arc.radius
        p.drawArc(QRectF(x - r, y - r, 2 * r, 2 * r), round(start * 16), round(span * 16))


def annotation_tip(canvas: "Canvas", x: float, y: float) -> str:
    """Tooltip of the annotation under a point of the canvas, or an empty string.

    Labels count within their text, lines within 4 pixels and arcs within
    4 pixels of their circle.
    """
    camera = canvas.camera
    w, h = canvas.width(), canvas.height()
    for layer in reversed(canvas.drawing.annotations()):
        for label in layer.labels:
            at = canvas.drawing.resolve_label(label) if label.tip else None
            if at is None:
                continue
            sx, sy = camera.to_screen(camera.to_view(at[None, :]), w, h)
            fm = QFontMetricsF(_font(canvas, label.size, bold=label.bold))
            half = fm.horizontalAdvance(label.text) / 2 + 3
            near = abs(y - sy[0] - label.offset[1]) <= fm.height() / 2 + 2
            if near and abs(x - sx[0] - label.offset[0]) <= half:
                return label.tip
        for line, a, b in canvas.drawing.resolve_lines(layer):
            if not line.tip:
                continue
            sx, sy = camera.to_screen(camera.to_view(np.array([a, b])), w, h)
            if _segment_distance(x, y, sx[0], sy[0], sx[1], sy[1]) <= 4.0:
                return line.tip
        for arc in layer.arcs:
            geometry = _arc_geometry(canvas, arc) if arc.tip else None
            if geometry is None:
                continue
            cx, cy, start, span = geometry
            angle = math.degrees(math.atan2(-(y - cy), x - cx))
            offset = (angle - start) % 360.0 if span >= 0 else (start - angle) % 360.0
            if abs(math.hypot(x - cx, y - cy) - arc.radius) <= 4.0 and offset <= abs(span):
                return arc.tip
    return ""


def _segment_distance(x: float, y: float, x0: float, y0: float, x1: float, y1: float) -> float:
    """Distance in pixels from a point to a line segment."""
    dx, dy = x1 - x0, y1 - y0
    length = dx * dx + dy * dy
    t = 0.0 if length == 0 else max(0.0, min(1.0, ((x - x0) * dx + (y - y0) * dy) / length))
    return math.hypot(x - x0 - t * dx, y - y0 - t * dy)


def _paint_layer_labels(p: QPainter, canvas: "Canvas", theme: Theme, layer: Annotations) -> None:
    """Labels of a plugin layer, each with a halo in the background colour."""
    camera = canvas.camera
    w, h = canvas.width(), canvas.height()
    for label in layer.labels:
        at = canvas.drawing.resolve_label(label)
        if at is None:
            continue
        points = [at]
        away = canvas.drawing.anchor_point(label.away_from) if label.away_from is not None else None
        if away is not None:
            points.append(away)
        sx, sy = camera.to_screen(camera.to_view(np.array(points)), w, h)
        x, y = float(sx[0]) + label.offset[0], float(sy[0]) + label.offset[1]
        font = _font(canvas, label.size, bold=label.bold)
        font.setItalic(label.italic)
        if away is not None:
            dx, dy = sx[0] - sx[1], sy[0] - sy[1]
            span = math.hypot(dx, dy)
            ux, uy = (dx / span, dy / span) if span > 1e-6 else (0.0, -1.0)
            # the near edge of the text, not its centre, sits ``gap`` past the anchor,
            # so the text never covers what it points at (an arrow tip, for example)
            box = QFontMetricsF(font).boundingRect(label.text)
            reach = label.gap + abs(ux) * box.width() / 2 + abs(uy) * box.height() / 2
            x, y = x + ux * reach, y + uy * reach
        halo = QColor(theme.background)
        halo.setAlphaF(label.halo)
        _outlined_text(p, x, y, label.text, QColor(label.color) if label.color else QColor(theme.text), halo,
                       font)







def _paint_stereo(p: QPainter, canvas: "Canvas", theme: Theme) -> None:
    mol, camera = canvas.mol, canvas.camera
    centers, bonds = canvas.stereo_labels()
    if not centers and not bonds:
        return
    font = _font(canvas, 10.5, bold=True)
    font.setItalic(True)
    halo = QColor(theme.background)
    halo.setAlphaF(0.9)
    sx, sy = camera.to_screen(camera.to_view(mol.coords), canvas.width(), canvas.height())
    hidden = canvas.doc.display.hidden
    for i, chirality in centers.items():
        if not hidden[i]:
            _outlined_text(p, sx[i] + 12, sy[i] + 12, str(chirality), theme.accent, halo, font)
    for (i, j), bond in bonds.items():
        if not (hidden[i] or hidden[j]):
            x, y = (sx[i] + sx[j]) / 2, (sy[i] + sy[j]) / 2 - 12
            _outlined_text(p, x, y, str(bond), theme.accent, halo, font)


def paint_measurements(p: QPainter, mol: Molecule, camera: Camera, width: float, height: float, theme: Theme,
                       base_font: QFont, hidden: np.ndarray | None = None) -> None:
    """Draw the pinned distances, angles and dihedrals of a molecule (also used for exported pictures).

    Measurements of atoms in the boolean mask ``hidden`` (focus mode) are left out.
    """
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    font = QFont(base_font)
    font.setPointSizeF(10.5)
    font.setBold(True)
    halo = QColor(theme.background)
    halo.setAlphaF(0.9)
    line = QColor(theme.text)
    line.setAlphaF(0.65)
    for atoms in mol.measurements:
        if any(i >= len(mol) or (hidden is not None and hidden[i]) for i in atoms):
            continue
        c = mol.coords[list(atoms)]
        sx, sy = camera.to_screen(camera.to_view(c), width, height)
        points = [QPointF(x, y) for x, y in zip(sx, sy, strict=True)]
        p.setPen(QPen(line, 1.4, Qt.PenStyle.DashLine, Qt.PenCapStyle.RoundCap))
        for a, b in pairwise(points):
            p.drawLine(a, b)
        if len(atoms) == 2:
            text, at = f"{distance(c[0], c[1]):.2f} {ANGSTROM}", (points[0] + points[1]) / 2
        elif len(atoms) == 3:
            text, at = f"{angle(c[0], c[1], c[2]):.1f}{DEGREE}", (points[0] + points[1] + points[2]) / 3
        else:
            text, at = f"{dihedral(c[0], c[1], c[2], c[3]):.1f}{DEGREE}", (points[1] + points[2]) / 2
        _outlined_text(p, at.x(), at.y() - 10, text, theme.text, halo, font)


_text_paths: dict[tuple[str, str], tuple[QPainterPath, float, float]] = {}


def _text_path(text: str, font: QFont) -> tuple[QPainterPath, float, float]:
    """Outline of a text on its baseline from the origin, and the shift that centres it on a point.

    The outline is kept for the next time the text is drawn.
    """
    key = (text, font.key())
    found = _text_paths.get(key)
    if found is None:
        if len(_text_paths) >= TEXT_CACHE_SIZE:
            _text_paths.clear()
        fm = QFontMetricsF(font)
        path = QPainterPath()
        path.addText(QPointF(0.0, 0.0), font, text)
        path.setFillRule(Qt.FillRule.WindingFill)   # glyphs like "+" overlap themselves
        found = (path, -fm.horizontalAdvance(text) / 2, (fm.ascent() - fm.descent()) / 2)
        _text_paths[key] = found
    return found


def _outlined_text(p: QPainter, x: float, y: float, text: str, fill: QColor, halo: QColor,
                   font: QFont) -> None:
    """Text centred on (x, y) with a halo that keeps it readable over atoms (none when fully transparent)."""
    outline, dx, dy = _text_path(text, font)
    path = outline.translated(x + dx, y + dy)
    if halo.alpha():
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.setPen(QPen(halo, 3.2, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap, Qt.PenJoinStyle.RoundJoin))
        p.drawPath(path)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(fill)
    p.drawPath(path)


def _view_moving(canvas: "Canvas") -> bool:
    """Whether the view is being turned, moved or zoomed by a drag, or is animating."""
    g = canvas.gesture
    turning = g is not None and g.moved and g.drag in VIEW_DRAGS
    return turning or canvas.view.animating or canvas.dragging()


def labelled_atoms(canvas: "Canvas", sx: np.ndarray, sy: np.ndarray, depth: np.ndarray) -> np.ndarray:
    """Atoms that get a label now, in index order.

    Only atoms drawn and inside the window (with :data:`LABEL_MARGIN`)
    count. Of those, at most :data:`LABEL_LIMIT` are labelled, the ones
    nearest the viewer; while the view moves, none are when more than
    :data:`MOVING_LABEL_LIMIT` would be.
    """
    w, h = canvas.width(), canvas.height()
    inside = (canvas.picker().pickable & (sx > -LABEL_MARGIN) & (sx < w + LABEL_MARGIN)
              & (sy > -LABEL_MARGIN) & (sy < h + LABEL_MARGIN))
    atoms = np.flatnonzero(inside)
    if len(atoms) > MOVING_LABEL_LIMIT and _view_moving(canvas):
        return atoms[:0]
    if len(atoms) > LABEL_LIMIT:
        nearest = np.argpartition(-depth[atoms], LABEL_LIMIT - 1)[:LABEL_LIMIT]
        atoms = np.sort(atoms[nearest])
    return atoms


_label_pictures: dict[tuple[str, str, int, int, float], QPixmap] = {}


def _label_picture(text: str, font: QFont, fill: QColor, halo: QColor, dpr: float) -> QPixmap:
    """Picture of an outlined text, kept for the next time it is drawn."""
    key = (text, font.key(), fill.rgba(), halo.rgba(), dpr)
    picture = _label_pictures.get(key)
    if picture is None:
        if len(_label_pictures) >= LABEL_PICTURES:
            _label_pictures.clear()
        path, dx, dy = _text_path(text, font)
        box = path.boundingRect().translated(dx, dy).adjusted(-HALO_PAD, -HALO_PAD, HALO_PAD, HALO_PAD)
        picture = QPixmap(max(1, math.ceil(box.width() * dpr)), max(1, math.ceil(box.height() * dpr)))
        picture.setDevicePixelRatio(dpr)
        picture.fill(Qt.GlobalColor.transparent)
        q = QPainter(picture)
        q.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        _outlined_text(q, -box.left(), -box.top(), text, fill, halo, font)
        q.end()
        _label_pictures[key] = picture
    return picture


def _paint_labels(p: QPainter, canvas: "Canvas", sx: np.ndarray, sy: np.ndarray,
                  depth: np.ndarray) -> None:
    mol = canvas.mol
    font = _font(canvas, 10.5, bold=True)
    halo = QColor(canvas.theme.background)
    halo.setAlphaF(0.85)
    fill = canvas.theme.text
    mode = canvas.display.labels
    symbols = mol.symbols
    atoms = labelled_atoms(canvas, sx, sy, depth)
    copy = len(atoms) > OUTLINED_LABELS   # many labels: copy cached pictures instead of filling outlines
    dpr = p.device().devicePixelRatioF() if copy else 1.0
    for i in atoms.tolist():
        if mode == LabelMode.SYMBOL:
            text = symbols[i]
        elif mode == LabelMode.NUMBER:
            text = str(i + 1)
        else:
            text = f"{symbols[i]}{i + 1}"
        if copy:
            picture = _label_picture(text, font, fill, halo, dpr)
            path, dx, dy = _text_path(text, font)
            box = path.boundingRect()
            p.drawPixmap(QPointF(sx[i] + dx + box.left() - HALO_PAD, sy[i] + dy + box.top() - HALO_PAD),
                         picture)
        else:
            _outlined_text(p, sx[i], sy[i], text, fill, halo, font)


def _paint_selection_order(p: QPainter, canvas: "Canvas", sx: np.ndarray, sy: np.ndarray) -> None:
    """Dashed path through 2-4 selected atoms with their numbers, as used by the selection panel."""
    selection = [i for i in canvas.selection if i < len(canvas.mol)]
    if not 2 <= len(selection) <= 4:
        return
    accent = canvas.theme.accent
    p.setPen(QPen(accent, 1.5, Qt.PenStyle.DashLine))
    for a, b in pairwise(selection):
        p.drawLine(QPointF(sx[a], sy[a]), QPointF(sx[b], sy[b]))
    p.setFont(_font(canvas, 9.5, bold=True))
    r = 8.0
    for n, i in enumerate(selection, start=1):
        center = QPointF(sx[i] + 13, sy[i] - 13)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(accent)
        p.drawEllipse(center, r, r)
        p.setPen(WHITE)
        p.drawText(QRectF(center.x() - r, center.y() - r, 2 * r, 2 * r), Qt.AlignmentFlag.AlignCenter, str(n))


def _paint_rubber_band(p: QPainter, canvas: "Canvas", sx: np.ndarray, sy: np.ndarray) -> None:
    """Line from the start atom of a build drag to the cursor or the target atom."""
    g = canvas.gesture
    assert g is not None and g.draw_from is not None
    start = QPointF(sx[g.draw_from], sy[g.draw_from])
    end = QPointF(sx[g.draw_to], sy[g.draw_to]) if g.draw_to is not None else g.last
    color = canvas.theme.attach if g.draw_to is not None else canvas.theme.accent
    p.setPen(QPen(color, 2.0, Qt.PenStyle.DashLine, Qt.PenCapStyle.RoundCap))
    p.drawLine(start, end)
    if g.draw_to is None:
        frag = canvas.builder.fragment
        name = frag.symbol if frag.category == Category.ELEMENT else frag.name
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(color)
        p.drawEllipse(end, 4.0, 4.0)
        _outlined_text(p, end.x() + 6 + len(name) * 3.5, end.y() - 12, name, canvas.theme.text,
                       QColor(canvas.theme.background), _font(canvas, 11, bold=True))


def _paint_empty_hint(p: QPainter, canvas: "Canvas") -> None:
    w, h = canvas.width(), canvas.height()
    if canvas.builder.tool == Tool.BUILD:
        headline = f"Click anywhere to place {canvas.builder.fragment.name.lower()}"
    else:
        headline = "Press Cmd 2 or type an element symbol to start building"
    p.setPen(canvas.theme.muted)
    p.setFont(_font(canvas, 15))
    p.drawText(QRectF(0, h / 2 - 34, w, 24), Qt.AlignmentFlag.AlignCenter, headline)
    p.setFont(_font(canvas, 12))
    p.drawText(QRectF(0, h / 2 - 6, w, 20), Qt.AlignmentFlag.AlignCenter,
               f"Type an element symbol to switch{SEPARATOR}1-9 picks the bonding type")


def _paint_toast(p: QPainter, canvas: "Canvas") -> None:
    theme = canvas.theme
    p.setFont(_font(canvas, 12.5))
    # a long message (a file's notes) wraps instead of running off the window
    flags = Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextWordWrap
    limit = max(120.0, min(canvas.width() - 48.0, TOAST_WIDTH))
    text = QFontMetricsF(p.font()).boundingRect(QRectF(0, 0, limit, 10_000), flags, canvas.toast_text)
    width, height = text.width() + 28, max(30.0, text.height() + 12)
    rect = QRectF((canvas.width() - width) / 2, canvas.toast_top, width, height)
    background = QColor(theme.panel)
    background.setAlphaF(0.96)
    p.setPen(QPen(theme.panel_border, 1))
    p.setBrush(background)
    p.drawRoundedRect(rect, 15, 15)
    p.setPen(theme.text)
    p.drawText(rect, flags, canvas.toast_text)
