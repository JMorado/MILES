"""Vibrational modes from frequency calculations: IR spectrum, mode list and animation."""

import os
from collections.abc import Callable
from typing import Any, override

import numpy as np
from PySide6.QtCore import QPoint, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QMouseEvent, QPainter, QPainterPath, QPaintEvent, QPen
from PySide6.QtWidgets import QHBoxLayout, QListWidget, QListWidgetItem, QVBoxLayout, QWidget

from miles import VibrationalMode
from miles.builtin_plugins.vibrations.modes import ir_spectrum, mode_arrows, mode_frames
from miles.builtin_plugins.vibrations.tables import mode_table, spectrum_table
from miles.plugin import Change, DocumentPanel, DocumentView, Theme
from miles.plugin.widgets import COPY, CURVE, DATA, IMAGE, Table, button, label, show_plot_menu

SPECTRUM_WIDTH = 20.0   # cm-1, full width at half maximum of each band
MODES_HEIGHT = 170      # pixels; the mode list is shorter when it has fewer rows
LOW, HIGH = 400.0, 4000.0
TICK = 500


def frequency_text(mode: VibrationalMode) -> str:
    """Wavenumber with its unit; imaginary frequencies get an ``i``."""
    if mode.frequency < 0:
        return f"{-mode.frequency:.0f}i cm\u207b\u00b9"
    return f"{mode.frequency:.0f} cm\u207b\u00b9"


class Spectrum(QWidget):
    """IR spectrum drawn as a broadened curve with a stick per mode; click a stick to choose it."""

    modeClicked = Signal(int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedSize(300, 120)
        self.theme: Theme | None = None
        self.modes: list[VibrationalMode] = []
        self.selected: int | None = None
        self._high = HIGH
        self._curve = np.zeros(0)
        self._grid = np.zeros(0)

    def set_modes(self, modes: list[VibrationalMode]) -> None:
        """Show the spectrum of these modes."""
        self.modes = modes
        top = max([HIGH] + [m.frequency + 200 for m in modes])
        self._high = float(np.ceil(top / TICK) * TICK)
        self._grid = np.linspace(LOW, self._high, 600)
        self._curve = ir_spectrum(modes, self._grid, SPECTRUM_WIDTH)
        self.update()

    @property
    def high(self) -> float:
        """Highest wavenumber on the axis, cm-1."""
        return self._high

    def _plot(self) -> QRectF:
        return QRectF(8, 6, self.width() - 16, self.height() - 26)

    def _x(self, wavenumber: float) -> float:
        plot = self._plot()   # highest wavenumbers on the left, as IR spectra are usually drawn
        return plot.left() + (self._high - wavenumber) / (self._high - LOW) * plot.width()

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        if self.theme is None:
            return
        p = QPainter(self)
        self.paint(p, self.theme)
        p.end()

    def axis_title(self) -> str:
        """Title under the axis of an exported picture."""
        return "Wavenumber (cm\u207b\u00b9)"

    def paint(self, p: QPainter, theme: Theme, for_export: bool = False) -> None:
        """Draw the spectrum; an exported picture leaves out the mark of the chosen mode."""
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        selected = None if for_export else self.selected
        plot = self._plot()
        font = QFont(self.font())
        font.setPixelSize(9)
        p.setFont(font)
        axis = QColor(theme.text)
        axis.setAlphaF(0.25)
        p.setPen(QPen(axis, 1))
        p.drawLine(QPointF(plot.left(), plot.bottom()), QPointF(plot.right(), plot.bottom()))
        p.setPen(theme.muted)
        metrics = QFontMetricsF(font)
        last_right = -1e9
        for w in range(int(self._high), int(LOW) - 1, -TICK):
            text = str(w)
            width = metrics.horizontalAdvance(text) + 4   # a little room for the glyphs' edges
            left = min(max(self._x(w) - width / 2, 1.0), self.width() - width - 1)   # end labels stay whole
            if left < last_right + 6:
                continue   # pushed inwards at an end, it would run into its neighbour
            p.drawText(QRectF(left, plot.bottom() + 3, width, 14), Qt.AlignmentFlag.AlignCenter, text)
            last_right = left + width
        intensities = [m.intensity or 0.0 for m in self.modes]
        top = max(intensities) or 1.0
        for k, mode in enumerate(self.modes):
            if mode.frequency < LOW:
                continue
            x = self._x(mode.frequency)
            height = plot.height() * (mode.intensity or 0.0) / top if mode.intensity else 6.0
            color = QColor(theme.accent if k == selected else theme.muted)
            p.setPen(QPen(color, 2.0 if k == selected else 1.0))
            p.drawLine(QPointF(x, plot.bottom()), QPointF(x, plot.bottom() - max(height, 3.0)))
        if self._curve.any():
            path = QPainterPath()
            for n, (wavenumber, value) in enumerate(zip(self._grid, self._curve, strict=True)):
                point = QPointF(self._x(float(wavenumber)), plot.bottom() - float(value) * plot.height())
                if n:
                    path.lineTo(point)
                else:
                    path.moveTo(point)
            p.setPen(QPen(theme.text, 1.3))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(path)

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if e.button() != Qt.MouseButton.LeftButton:
            return
        x = e.position().x()
        candidates = [(abs(self._x(m.frequency) - x), k)
                      for k, m in enumerate(self.modes) if m.frequency >= LOW]
        if candidates:
            nearest = min(candidates)
            if nearest[0] < 8:
                self.modeClicked.emit(nearest[1])


class VibrationsPanel(DocumentPanel):
    """Floating panel listing the vibrational modes of the current document.

    Choosing a mode animates it through the timeline; Stop restores the
    geometry and frames the document had before. Arrows shows how each atom
    moves in the chosen mode, with or without the animation: the panel keeps
    the displacements of each tab (:meth:`arrows_for`) and emits
    ``arrowsChanged`` for whoever draws them. It opens by itself for a
    structure with modes.
    """

    refresh_on = frozenset({Change.STRUCTURE, Change.FRAMES, Change.RESULTS})
    unavailable = "No vibrational modes; open the output of a frequency calculation"

    playRequested = Signal()
    arrowsChanged = Signal()
    exportRequested = Signal(str)   # one of the kinds in miles.ui.plot_export

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.wanted = True   # it opens by itself for a structure with modes
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(6)
        top = QHBoxLayout()
        top.addWidget(self.chevron)
        top.addWidget(label("Vibrations", "Title"))
        top.addStretch(1)
        self.arrows = button("Seg", "Arrows", "Show how each atom moves in the chosen mode", checkable=True)
        self.arrows.clicked.connect(self._show_arrows)
        top.addWidget(self.arrows)
        self.stop = button("Link", "Stop")
        self.stop.clicked.connect(self.stop_animation)
        top.addWidget(self.stop)
        close = button("Link", "Close")
        close.clicked.connect(self.close_panel)
        top.addWidget(close)
        lay.addLayout(top)
        self.spectrum = Spectrum()
        self.spectrum.modeClicked.connect(self.animate)
        self.spectrum.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.spectrum.customContextMenuRequested.connect(lambda at: self._menu(self.spectrum, at))
        self.body_layout.setSpacing(6)
        self.body_layout.addWidget(self.spectrum)
        self.modes = QListWidget()
        self.modes.setObjectName("Modes")
        self.modes.setFixedHeight(MODES_HEIGHT)
        self.modes.currentRowChanged.connect(self._row_chosen)
        self.modes.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.modes.customContextMenuRequested.connect(lambda at: self._menu(self.modes, at))
        self.body_layout.addWidget(self.modes)
        lay.addWidget(self.body)

    def _kept(self, doc: DocumentView) -> dict[str, Any]:
        """Return what the panel keeps for a tab, dropped when the tab closes.

        It holds the geometry an animation moves around (``"saved"``), the
        mode last chosen (``"mode"``) and the displacement arrows shown
        (``"arrows"``).
        """
        return self.store(doc)

    def arrows_for(self, doc: DocumentView) -> np.ndarray | None:
        """Return the displacement arrows shown in a tab, one row per atom, or None."""
        arrows: np.ndarray | None = self._kept(doc).get("arrows") if self.context is not None else None
        return arrows

    def menu_entries(self) -> list[tuple[str, Callable[[], None]]]:
        """Entries of the right-click menu of the spectrum and the mode list."""
        return [("Copy Data", lambda: self.exportRequested.emit(COPY)),
                ("Export Data...", lambda: self.exportRequested.emit(DATA)),
                ("Export Spectrum Curve...", lambda: self.exportRequested.emit(CURVE)),
                ("Export Image...", lambda: self.exportRequested.emit(IMAGE))]

    def _menu(self, widget: QWidget, at: QPoint) -> None:
        if self.document is not None and self.document.info.vibrations:
            show_plot_menu(widget, at, self.menu_entries())

    def export_name(self) -> str:
        """File name to suggest, e.g. ``water_ir``."""
        name = os.path.splitext(self.document.name)[0] if self.document is not None else ""
        return f"{name or 'spectrum'}_ir"

    def export_table(self, kind: str) -> Table | None:
        """Tabulate the modes, or for :data:`CURVE` the broadened spectrum."""
        modes = self.document.info.vibrations if self.document is not None else []
        if not modes:
            return None
        if kind == CURVE:
            return spectrum_table(modes, LOW, self.spectrum.high, 1.0, SPECTRUM_WIDTH)
        return mode_table(modes)

    def plot_widget(self) -> Spectrum:
        """Return the spectrum, for exported pictures."""
        return self.spectrum

    @override
    def set_theme(self, theme: Theme) -> None:
        """Recolour the spectrum."""
        self.spectrum.theme = theme
        self.spectrum.update()
        self.set_icon_color(theme.muted)

    @override
    def has_content(self, document: DocumentView) -> bool:
        """Whether the document has vibrational modes."""
        return bool(document.info.vibrations)

    @override
    def show_contents(self, document: DocumentView) -> None:
        """Fill the spectrum and the list for the document."""
        modes = document.info.vibrations
        self.spectrum.set_modes(modes)
        self.modes.blockSignals(True)
        self.modes.clear()
        for k, mode in enumerate(modes, start=1):
            intensity = f"   {mode.intensity:.0f} km/mol" if mode.intensity is not None else ""
            self.modes.addItem(QListWidgetItem(f"{k:3d}   {frequency_text(mode)}{intensity}"))
        rows = sum(self.modes.sizeHintForRow(k) for k in range(self.modes.count()))
        frame = 2 * self.modes.frameWidth() + 4
        self.modes.setFixedHeight(min(MODES_HEIGHT, rows + frame))
        doc = document
        chosen = self._kept(doc).get("mode")
        if chosen is not None and chosen < len(modes):
            self.modes.setCurrentRow(chosen)
        self.modes.blockSignals(False)
        arrows = self.arrows_for(doc) is not None
        self.arrows.setChecked(arrows)
        self.spectrum.selected = chosen if self._animating(doc) or arrows else None
        self.stop.setVisible(self._animating(doc))

    @override
    def closing(self) -> None:
        """Stop any animation and hide the arrows when the user closes the panel."""
        self.arrows.setChecked(False)
        self._show_arrows()
        self.stop_animation()

    def _animating(self, doc: DocumentView) -> bool:
        kept = self._kept(doc)
        if "saved" in kept and not doc.animating:   # an edit or undo ended it
            kept.pop("saved")
        return "saved" in kept

    def _row_chosen(self, row: int) -> None:
        if row >= 0:
            self.animate(row)

    def animate(self, index: int) -> None:
        """Play the vibration of mode ``index`` in a loop."""
        if self.document is None:
            return
        doc = self.document
        modes = doc.info.vibrations
        if not 0 <= index < len(modes):
            return
        kept = self._kept(doc)
        if not self._animating(doc):   # the geometry the modes move around (the document's own)
            kept["saved"] = doc.molecule.coords.copy()
        base = kept["saved"]
        if base.shape != doc.molecule.coords.shape:
            kept.pop("saved")
            return
        kept["mode"] = index
        frames = mode_frames(base, modes[index])
        doc.show_frames(frames, [frequency_text(modes[index])] * len(frames), loop=True)
        self.spectrum.selected = index
        self.spectrum.update()
        self.modes.blockSignals(True)
        self.modes.setCurrentRow(index)
        self.modes.blockSignals(False)
        self.stop.setVisible(True)
        self._show_arrows()
        self.playRequested.emit()

    def stop_animation(self) -> None:
        """Return to the geometry and frames from before the animation."""
        if self.document is None or "saved" not in self._kept(self.document):
            return
        doc = self.document
        self._kept(doc).pop("saved")
        doc.show_frames(None)   # the document's own frames and geometry come back
        self.spectrum.selected = self._kept(doc).get("mode") if self.arrows.isChecked() else None
        self.spectrum.update()
        self.stop.setVisible(False)

    def _show_arrows(self) -> None:
        """Draw or hide the displacement arrows of the chosen mode, following the Arrows button."""
        doc = self.document
        if doc is None:
            return
        modes = doc.info.vibrations
        kept = self._kept(doc)
        if not self.arrows.isChecked() or not modes:
            if kept.pop("arrows", None) is not None:
                self.arrowsChanged.emit()
            return
        index = kept.get("mode")
        if index is None or index >= len(modes):   # nothing chosen yet: take the highlighted row
            index = min(max(self.modes.currentRow(), 0), len(modes) - 1)
            kept["mode"] = index
            self.modes.blockSignals(True)
            self.modes.setCurrentRow(index)
            self.modes.blockSignals(False)
        self.spectrum.selected = index
        self.spectrum.update()
        kept["arrows"] = mode_arrows(modes[index])
        self.arrowsChanged.emit()
