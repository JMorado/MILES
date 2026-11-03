"""Ramachandran plot: backbone phi and psi of every residue, linked to the 3D selection."""

import math
import os
from collections.abc import Callable
from typing import override

import numpy as np
from PySide6.QtCore import QEvent, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QHelpEvent,
    QMouseEvent,
    QPainter,
    QPainterPath,
    QPaintEvent,
    QPen,
    QPolygonF,
)
from PySide6.QtWidgets import QHBoxLayout, QToolTip, QVBoxLayout, QWidget

from miles.builtin_plugins.ramachandran.angles import (
    REGIONS,
    BackboneResidue,
    Ellipse,
    Region,
    ResidueKind,
    backbone_angles,
    backbone_residues,
    classify,
    ramachandran_table,
)
from miles.plugin import DEGREE, SEPARATOR, Change, DocumentPanel, DocumentView, Theme, plural
from miles.plugin.widgets import COPY, DATA, IMAGE, Table, button, label, show_plot_menu

PHI, PSI = "\u03c6", "\u03c8"
ALLOWED_COLOR = QColor("#E0902A")
TICKS = (-180, -90, 0, 90, 180)
PICK_RADIUS = 7.0   # pixels
KIND_NAMES = {ResidueKind.GLYCINE: "glycine", ResidueKind.PROLINE: "proline",
              ResidueKind.PRE_PROLINE: "before proline"}


class RamachandranPlot(QWidget):
    """Phi/psi scatter plot over shaded favoured and allowed regions.

    Glycines are drawn as triangles and prolines as squares; the colour
    gives the region. Click a point to select its residue, Shift-click to
    add it to the selection.
    """

    residueClicked = Signal(int, bool)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedSize(312, 318)
        self.setMouseTracking(True)
        self.theme: Theme | None = None
        self.residues: list[BackboneResidue] = []
        self.angles = np.zeros((0, 2))
        self.regions: list[Region | None] = []
        self.selected: set[int] = set()
        self.hover: int | None = None

    def set_data(self, found: list[BackboneResidue], angles: np.ndarray) -> None:
        """Show these residues with their phi and psi (NaN where undefined)."""
        self.residues = found
        self.angles = angles
        self.regions = [None if np.isnan(phi) or np.isnan(psi)
                        else classify(float(phi), float(psi), r.kind)
                        for r, (phi, psi) in zip(found, angles, strict=True)]
        if self.hover is not None and self.hover >= len(found):
            self.hover = None
        self.update()

    def _plot(self) -> QRectF:
        return QRectF(44, 6, 260, 260)

    def _point(self, phi: float, psi: float) -> QPointF:
        plot = self._plot()
        return QPointF(plot.left() + (phi + 180) / 360 * plot.width(),
                       plot.top() + (180 - psi) / 360 * plot.height())

    def _region_path(self, ellipses: tuple[Ellipse, ...]) -> QPainterPath:
        """Union of ellipses, repeated across the +-180 boundaries."""
        plot = self._plot()
        sx, sy = plot.width() / 360, plot.height() / 360
        path = QPainterPath()
        path.setFillRule(Qt.FillRule.WindingFill)
        for e in ellipses:
            for dphi in (-360, 0, 360):
                for dpsi in (-360, 0, 360):
                    center = self._point(e.phi + dphi, e.psi + dpsi)
                    path.addEllipse(center, e.half_phi * sx, e.half_psi * sy)
        return path

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        if self.theme is None:
            return
        p = QPainter(self)
        self.paint(p, self.theme)
        p.end()

    def paint(self, p: QPainter, theme: Theme, for_export: bool = False) -> None:
        """Draw the plot; an exported picture leaves out the selection and hover marks."""
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        plot = self._plot()
        font = QFont(self.font())
        font.setPixelSize(9)
        p.setFont(font)
        p.save()
        p.setClipRect(plot)
        favoured, allowed = REGIONS[ResidueKind.GENERAL]
        p.setPen(Qt.PenStyle.NoPen)
        for ellipses, alpha in ((allowed, 0.10), (favoured, 0.22)):
            fill = QColor(theme.accent)
            fill.setAlphaF(alpha)
            p.setBrush(fill)
            p.drawPath(self._region_path(ellipses))
        p.restore()
        grid = QColor(theme.text)
        grid.setAlphaF(0.12)
        p.setPen(QPen(grid, 1, Qt.PenStyle.DashLine))
        for t in (-90, 0, 90):
            p.drawLine(self._point(t, -180), self._point(t, 180))
            p.drawLine(self._point(-180, t), self._point(180, t))
        frame = QColor(theme.text)
        frame.setAlphaF(0.3)
        p.setPen(QPen(frame, 1))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(plot)
        p.setPen(theme.muted)
        for t in TICKS:
            x = self._point(t, 0).x()
            p.drawText(QRectF(x - 16, plot.bottom() + 2, 32, 12), Qt.AlignmentFlag.AlignCenter, str(t))
            y = self._point(0, t).y()
            p.drawText(QRectF(14, y - 6, plot.left() - 18, 12),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, str(t))
        title = QFont(font)
        title.setPixelSize(12)
        p.setFont(title)
        center = Qt.AlignmentFlag.AlignCenter
        p.drawText(QRectF(plot.left(), plot.bottom() + 12, plot.width(), 16), center, PHI)
        p.drawText(QRectF(0, plot.center().y() - 8, 14, 16), center, PSI)
        self._paint_points(p, theme, marks=not for_export)
        p.setFont(font)
        self._paint_legend(p, theme, plot.bottom() + 34)

    def _color(self, region: Region, theme: Theme) -> QColor:
        if region == Region.FAVOURED:
            return QColor(theme.text)
        if region == Region.ALLOWED:
            return QColor(ALLOWED_COLOR)
        return QColor(theme.danger)

    def _marker(self, p: QPainter, center: QPointF, kind: ResidueKind, size: float) -> None:
        if kind == ResidueKind.GLYCINE:
            h = size * 1.15
            x, y = center.x(), center.y()
            low = y + h * 0.8
            p.drawPolygon(QPolygonF([QPointF(x, y - h), QPointF(x + h, low), QPointF(x - h, low)]))
        elif kind == ResidueKind.PROLINE:
            p.drawRect(QRectF(center.x() - size, center.y() - size, 2 * size, 2 * size))
        else:
            p.drawEllipse(center, size, size)

    def _paint_points(self, p: QPainter, theme: Theme, marks: bool = True) -> None:
        selected = self.selected if marks else set()
        dim = bool(selected)
        # outliers and selected residues last, so they stay on top
        order = sorted(range(len(self.residues)),
                       key=lambda k: (k in selected, self.regions[k] != Region.FAVOURED))
        for k in order:
            region = self.regions[k]
            if region is None:
                continue
            phi, psi = self.angles[k]
            center = self._point(float(phi), float(psi))
            color = self._color(region, theme)
            if dim and k not in selected:
                color.setAlphaF(0.45)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(color)
            self._marker(p, center, self.residues[k].kind, 3.0)
            if k in selected:
                p.setPen(QPen(theme.accent, 1.6))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawEllipse(center, 6.0, 6.0)
        if marks and self.hover is not None and self.regions[self.hover] is not None:
            phi, psi = self.angles[self.hover]
            p.setPen(QPen(theme.hover, 1.6))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawEllipse(self._point(float(phi), float(psi)), 7.0, 7.0)

    def _paint_legend(self, p: QPainter, theme: Theme, y: float) -> None:
        x = self._plot().left()
        entries = [(Region.FAVOURED, ResidueKind.GENERAL, "Favoured"),
                   (Region.ALLOWED, ResidueKind.GENERAL, "Allowed"),
                   (Region.OUTLIER, ResidueKind.GENERAL, "Outlier"),
                   (None, ResidueKind.GLYCINE, "Gly"), (None, ResidueKind.PROLINE, "Pro")]
        for region, kind, text in entries:
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(self._color(region, theme) if region is not None else theme.muted)
            self._marker(p, QPointF(x + 3, y + 6), kind, 3.0)
            p.setPen(theme.muted)
            width = p.fontMetrics().horizontalAdvance(text)
            p.drawText(QRectF(x + 10, y, width + 2, 12), Qt.AlignmentFlag.AlignVCenter, text)
            x += width + 24

    def residue_at(self, pos: QPointF) -> int | None:
        """Residue whose point is under a position, or None."""
        best, nearest = None, PICK_RADIUS
        for k, region in enumerate(self.regions):
            if region is None:
                continue
            phi, psi = self.angles[k]
            point = self._point(float(phi), float(psi))
            d = math.hypot(point.x() - pos.x(), point.y() - pos.y())
            if d < nearest:
                best, nearest = k, d
        return best

    def describe(self, k: int) -> str:
        """Tooltip text of a residue: name, angles and region."""
        r = self.residues[k]
        phi, psi = self.angles[k]
        text = f"{r.resname} {r.resseq}{' ' + r.chain if r.chain else ''}"
        text += f"{SEPARATOR}{PHI} {phi:.0f}{DEGREE}  {PSI} {psi:.0f}{DEGREE}"
        region = self.regions[k]
        if region is not None:
            text += f"{SEPARATOR}{region}"
        if r.kind in KIND_NAMES:
            text += f" ({KIND_NAMES[r.kind]})"
        return text

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        k = self.residue_at(e.position())
        if k != self.hover:
            self.hover = k
            self.setCursor(Qt.CursorShape.PointingHandCursor if k is not None else Qt.CursorShape.ArrowCursor)
            self.update()

    @override
    def leaveEvent(self, event: QEvent) -> None:
        if self.hover is not None:
            self.hover = None
            self.update()

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if e.button() != Qt.MouseButton.LeftButton:
            return
        k = self.residue_at(e.position())
        if k is not None:
            adding = Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier
            self.residueClicked.emit(k, bool(e.modifiers() & adding))

    @override
    def event(self, e: QEvent) -> bool:
        if e.type() == QEvent.Type.ToolTip and isinstance(e, QHelpEvent):
            k = self.residue_at(QPointF(e.pos()))
            if k is not None:
                QToolTip.showText(e.globalPos(), self.describe(k), self)
            else:
                QToolTip.hideText()
            return True
        return super().event(e)


class RamachandranPanel(DocumentPanel):
    """Floating panel with the Ramachandran plot of the current document.

    The plot follows edits, trajectory frames and the 3D selection; its
    right-click menu copies or exports the angles and the picture.
    """

    unavailable = "No protein backbone to plot"

    exportRequested = Signal(str)   # one of the kinds in miles.ui.plot_export
    toolRequested = Signal(str)     # the mouse tool to use, before residues are picked in the plot

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._found: list[BackboneResidue] | None = None
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(4)
        top = QHBoxLayout()
        top.addWidget(self.chevron)
        top.addWidget(label("Ramachandran", "Title"))
        top.addStretch(1)
        close = button("Link", "Close")
        close.clicked.connect(self.close_panel)
        top.addWidget(close)
        lay.addLayout(top)
        self.body_layout.setSpacing(4)
        self.summary = label("", "Muted")
        self.body_layout.addWidget(self.summary)
        self.plot = RamachandranPlot()
        self.plot.residueClicked.connect(self._clicked)
        self.plot.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.plot.customContextMenuRequested.connect(
            lambda at: show_plot_menu(self.plot, at, self.menu_entries()))
        self.body_layout.addWidget(self.plot)
        lay.addWidget(self.body)

    @override
    def set_theme(self, theme: Theme) -> None:
        """Recolour the plot."""
        self.plot.theme = theme
        self.plot.update()
        self.set_icon_color(theme.muted)

    def menu_entries(self) -> list[tuple[str, Callable[[], None]]]:
        """Entries of the plot's right-click menu."""
        return [("Copy Data", lambda: self.exportRequested.emit(COPY)),
                ("Export Data...", lambda: self.exportRequested.emit(DATA)),
                ("Export Image...", lambda: self.exportRequested.emit(IMAGE))]

    def export_name(self) -> str:
        """File name to suggest, e.g. ``1ubq_ramachandran``."""
        name = os.path.splitext(self.document.name)[0] if self.document is not None else ""
        return f"{name or 'protein'}_ramachandran"

    def export_table(self, kind: str) -> Table | None:
        """Tabulate chain, residue, phi, psi and region of every residue plotted."""
        if not self.plot.residues:
            return None
        return ramachandran_table(self.plot.residues, self.plot.angles)

    def plot_widget(self) -> RamachandranPlot:
        """Return the plot, for exported pictures."""
        return self.plot

    @override
    def set_document(self, document: DocumentView | None) -> None:
        """Follow another tab."""
        self._found = None
        super().set_document(document)

    def _residues(self) -> list[BackboneResidue]:
        if self._found is None:
            self._found = backbone_residues(self.document.molecule) if self.document is not None else []
        return self._found

    def _has_angles(self) -> bool:
        return any(r.phi_atoms is not None or r.psi_atoms is not None for r in self._residues())

    @override
    def has_content(self, document: DocumentView) -> bool:
        """Whether the document has a protein backbone with angles to plot."""
        return document == self.document and self._has_angles()

    @override
    def document_changed(self, kind: str) -> None:
        """Recompute after edits and frame changes of the current document."""
        if kind in (Change.STRUCTURE, Change.BONDS, Change.FRAMES):
            self._found = None
        if self.wanted:
            self.refresh()

    @override
    def show_contents(self, document: DocumentView) -> None:
        """Compute the angles and plot them."""
        found = self._residues()
        self.plot.set_data(found, backbone_angles(found, document.molecule.coords))
        counts = {region: self.plot.regions.count(region) for region in Region}
        total = sum(counts.values())
        share = 100 * counts[Region.FAVOURED] / total if total else 0.0
        self.summary.setText(f"{plural(total, 'residue')}{SEPARATOR}{share:.0f}% favoured"
                             f"{SEPARATOR}{plural(counts[Region.OUTLIER], 'outlier')}")
        self.selection_changed()

    @override
    def selection_changed(self) -> None:
        """Ring the points of residues with selected atoms."""
        if self.document is None or not self.wanted:
            return
        chosen = set(self.document.selection)
        self.plot.selected = {k for k, r in enumerate(self._residues()) if chosen.intersection(r.atoms)}
        self.plot.update()

    def _clicked(self, k: int, additive: bool) -> None:
        doc = self.document
        if doc is None or k >= len(self._residues()):
            return
        atoms = list(self._residues()[k].atoms)
        base = doc.selection if additive else []
        self.toolRequested.emit("select")
        doc.select(base + [i for i in atoms if i not in base])
