"""Strip along the top of the window with the sequence of the current structure."""

from dataclasses import dataclass
from typing import override

from PySide6.QtCore import QEvent, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QHelpEvent,
    QMouseEvent,
    QPainter,
    QPaintEvent,
    QResizeEvent,
    QWheelEvent,
)
from PySide6.QtWidgets import QToolTip, QWidget

from miles.core.protein import ONE_LETTER, residues
from miles.render.cartoon import SS_COLORS
from miles.render.theme import Theme
from miles.ui.builder import Tool
from miles.ui.canvas import Canvas
from miles.ui.text import residue_label

CELL = 12.0        # width of one residue in pixels
CHAIN_GAP = 22.0   # space before each chain, holding its name
HEIGHT = 30
PADDING = 10.0


@dataclass(frozen=True, slots=True)
class _Residue:
    atoms: tuple[int, ...]
    letter: str
    ss: str
    chain: str
    x: float   # left edge in content coordinates


class SequenceStrip(QWidget):
    """One cell per residue, coloured by secondary structure; click or drag to select residues."""

    relayout = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(HEIGHT)
        self.setMouseTracking(True)
        self.canvas: Canvas | None = None
        self.theme: Theme | None = None
        self._residues: list[_Residue] = []
        self._chain_labels: list[tuple[str, float]] = []
        self._width = 0.0
        self._offset = 0.0
        self._key: tuple[int, int] | None = None
        self._anchor: int | None = None
        self._base_selection: list[int] = []
        self.hide()

    # -- data -----------------------------------------------------------------------------------

    def set_canvas(self, canvas: Canvas | None) -> None:
        """Follow another canvas."""
        self.canvas = canvas
        self._offset = 0.0
        self.refresh()

    def refresh(self) -> None:
        """Rebuild the residue list if the structure changed (not when atoms only moved), and repaint."""
        canvas = self.canvas
        key = (id(canvas.mol), canvas.doc.layout_revision) if canvas is not None else None
        if key != self._key:
            self._key = key
            self._layout_residues()
        visible = bool(self._residues)
        if visible == self.isHidden():   # isHidden: right before the window is shown, too
            self.setVisible(visible)
            self.relayout.emit()
        self.update()

    def _layout_residues(self) -> None:
        self._residues, self._chain_labels = [], []
        if self.canvas is None:
            return
        mol = self.canvas.mol
        x = PADDING
        chain = None
        for atoms in residues(mol):
            a = mol.atoms[atoms[0]]
            if a.chain != chain:
                chain = a.chain
                self._chain_labels.append((chain or "-", x))
                x += CHAIN_GAP
            ss = next((mol.atoms[i].ss for i in atoms if mol.atoms[i].ss), "")
            self._residues.append(_Residue(tuple(atoms), ONE_LETTER.get(a.resname, "X"), ss, a.chain, x))
            x += CELL
        self._width = x + PADDING
        self._clamp_offset()

    def _clamp_offset(self) -> None:
        self._offset = max(0.0, min(self._offset, self._width - self.width()))

    def _residue_at(self, x: float) -> int | None:
        content_x = x + self._offset
        for k, r in enumerate(self._residues):
            if r.x <= content_x < r.x + CELL:
                return k
        return None

    # -- drawing --------------------------------------------------------------------------------

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        if self.theme is None or self.canvas is None:
            return
        theme = self.theme
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        background = QColor(theme.panel)
        background.setAlphaF(0.94)
        p.setPen(theme.panel_border)
        p.setBrush(background)
        p.drawRoundedRect(QRectF(0.5, 0.5, self.width() - 1, self.height() - 1), 10, 10)
        p.setClipRect(QRectF(4, 0, self.width() - 8, self.height()))

        font = QFont(self.font())
        font.setPixelSize(10)
        font.setBold(True)
        p.setFont(font)
        selected = set(self.canvas.selection)
        cell_top, cell_height = 7.0, self.height() - 14.0
        for name, x in self._chain_labels:
            p.setPen(theme.muted)
            rect = QRectF(x - self._offset, 0, CHAIN_GAP - 4, self.height())
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, name)
        for r in self._residues:
            x = r.x - self._offset
            if x + CELL < 0 or x > self.width():
                continue
            rect = QRectF(x + 0.5, cell_top, CELL - 1, cell_height)
            if selected.intersection(r.atoms):
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(theme.accent)
                p.drawRoundedRect(rect, 3, 3)
                p.setPen(QColor(255, 255, 255))
            else:
                color = QColor(SS_COLORS.get(r.ss, SS_COLORS[""]))
                color.setAlphaF(0.35 if r.ss else 0.12)
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(color)
                p.drawRoundedRect(rect, 3, 3)
                p.setPen(theme.text)
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, r.letter)
        p.end()

    # -- input ----------------------------------------------------------------------------------

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        k = self._residue_at(e.position().x())
        if k is None or self.canvas is None:
            return
        adding = Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier
        additive = bool(e.modifiers() & adding)
        self._base_selection = list(self.canvas.selection) if additive else []
        self._anchor = k
        self._select_range(k, k)

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        if self._anchor is None:
            return
        k = self._residue_at(e.position().x())
        if k is not None:
            self._select_range(self._anchor, k)

    @override
    def mouseReleaseEvent(self, e: QMouseEvent) -> None:
        self._anchor = None

    def _select_range(self, a: int, b: int) -> None:
        if self.canvas is None:
            return
        lo, hi = min(a, b), max(a, b)
        atoms = [i for r in self._residues[lo:hi + 1] for i in r.atoms]
        self.canvas.builder.set_tool(Tool.SELECT)
        self.canvas.set_selection(self._base_selection + [i for i in atoms if i not in self._base_selection])
        self.update()

    @override
    def wheelEvent(self, e: QWheelEvent) -> None:
        delta = e.pixelDelta().x() or e.pixelDelta().y() or (e.angleDelta().y() or e.angleDelta().x()) / 4
        self._offset -= delta
        self._clamp_offset()
        self.update()

    @override
    def event(self, e: QEvent) -> bool:
        if e.type() == QEvent.Type.ToolTip and isinstance(e, QHelpEvent) and self.canvas is not None:
            k = self._residue_at(QPointF(e.pos()).x())
            if k is not None:
                atom = self.canvas.mol.atoms[self._residues[k].atoms[0]]
                QToolTip.showText(e.globalPos(), residue_label(atom))
            else:
                QToolTip.hideText()
            return True
        return super().event(e)

    @override
    def resizeEvent(self, event: QResizeEvent) -> None:
        self._clamp_offset()
        super().resizeEvent(event)
