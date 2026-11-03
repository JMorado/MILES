"""Timeline for structures with several frames: trajectories, optimisations, conformers, vibrations."""

import os
import time
import weakref
from collections.abc import Callable
from dataclasses import dataclass
from typing import override

import numpy as np
from PySide6.QtCore import QPointF, QRectF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QMouseEvent, QPainter, QPainterPath, QPaintEvent, QPen
from PySide6.QtWidgets import QComboBox, QFrame, QGridLayout, QSlider, QToolButton, QWidget

from miles.core.tables import Table, frame_table
from miles.core.trajectory import measurement_series, rmsd_atoms, rmsd_series, unwrap_degrees
from miles.render.theme import Theme
from miles.ui import icons
from miles.ui.controls import Segmented, label
from miles.ui.document import Change, Document
from miles.ui.plot_export import COPY, DATA, show_plot_menu
from miles.ui.text import ANGSTROM, DEGREE, SEPARATOR, atom_label

TICK_MS = 16          # redraw about 60 times a second while playing
BASE_RATE = 10.0      # frames per second at 1x
SPEEDS = [("0.25x", 0.25), ("0.5x", 0.5), ("1x", 1.0), ("2x", 2.0), ("4x", 4.0)]
DELTA = "\u0394"
PLOT_HEIGHT = 46
HANDLE_RADIUS = 7     # half the width of the slider handle, so the plot lines up with it
DOTS_UP_TO = 60       # mark each frame on the curve when there are this few
ENERGY, RMSD = "energy", "rmsd"
KCAL = " kcal/mol"
TEXT_WIDTH = 150      # least width of the frame text


def relative_energy_text(energies: list[float], frame: int) -> str:
    """Energy of a frame above the lowest, e.g. ``"dE 1.23 kcal/mol"`` (with a Greek Delta)."""
    return f"{DELTA}E {energies[frame] - min(energies):.2f}{KCAL}"


@dataclass(frozen=True)
class Quantity:
    """Something the timeline can plot along the frames.

    Attributes
    ----------
    key : str
        ``"energy"``, ``"rmsd"`` or the atoms of a pinned measurement, e.g. ``"3,4,5"``.
    name : str
        Text in the chooser, e.g. ``"Dihedral C1-C2-C3-N4"``.
    short : str
        Text before the value in the frame text, e.g. ``"C1-C2-C3-N4"``.
    unit : str
        Unit written after values, with its leading space if it needs one.
    decimals : int
        Decimals shown in the frame text.
    """

    key: str
    name: str
    short: str
    unit: str
    decimals: int

    @property
    def atoms(self) -> list[int]:
        """Atoms of a measurement (empty for energy and RMSD)."""
        return [] if self.key in (ENERGY, RMSD) else [int(i) for i in self.key.split(",")]


def quantities(doc: Document) -> list[Quantity]:
    """List what can be plotted for a document's frames: energy, pinned measurements, RMSD."""
    found = []
    if doc.frame_energies is not None:
        found.append(Quantity(ENERGY, "Energy", f"{DELTA}E", KCAL, 2))
    mol = doc.mol
    for atoms in mol.measurements:
        if not 2 <= len(atoms) <= 4 or any(i >= len(mol) for i in atoms):
            continue
        names = "-".join(atom_label(mol, i) for i in atoms)
        kind = {2: "Distance", 3: "Angle", 4: "Dihedral"}[len(atoms)]
        unit, decimals = (f" {ANGSTROM}", 2) if len(atoms) == 2 else (DEGREE, 1)
        found.append(Quantity(",".join(map(str, atoms)), f"{kind} {names}", names, unit, decimals))
    found.append(Quantity(RMSD, "RMSD", "RMSD", f" {ANGSTROM}", 2))
    return found


class FramePlot(QWidget):
    """One value per frame as a curve, with the shown position marked.

    Energies are drawn relative to the lowest, from zero up; other
    quantities between their smallest and largest values. Clicking or
    dragging along the plot emits :attr:`seek` with the nearest frame.
    """

    seek = Signal(int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(PLOT_HEIGHT)
        self.setToolTip("Click or drag to go to a frame")
        self.theme: Theme | None = None
        self.values = np.zeros(0)
        self.caption = ""
        self.from_zero = True
        self.position = 0.0
        self._last = -1

    def set_values(self, values: np.ndarray, caption: str, from_zero: bool) -> None:
        """Plot these values.

        Parameters
        ----------
        values : numpy.ndarray
            One value per frame.
        caption : str
            Text at the top left, e.g. the range and unit.
        from_zero : bool
            Draw from zero up (relative energies) rather than between the
            smallest and largest value.
        """
        self.values = np.asarray(values, dtype=float)
        self.caption = caption
        self.from_zero = from_zero
        self.update()

    def set_energies(self, energies: list[float]) -> None:
        """Plot energies (kcal/mol) relative to the lowest."""
        values = np.asarray(energies, dtype=float)
        values = values - values.min() if values.size else values
        top = float(values.max()) if values.size else 0.0
        self.set_values(values, f"{top:.1f}{KCAL}", from_zero=True)

    def set_position(self, position: float) -> None:
        """Mark a frame, or a point between two frames during playback."""
        self.position = position
        self.update()

    def _plot(self) -> QRectF:
        return QRectF(HANDLE_RADIUS, 14, self.width() - 2 * HANDLE_RADIUS, self.height() - 17)

    def _x(self, frame: float) -> float:
        plot = self._plot()
        return plot.left() + frame / max(len(self.values) - 1, 1) * plot.width()

    def _range(self) -> tuple[float, float]:
        low = 0.0 if self.from_zero else float(self.values.min())
        high = float(self.values.max())
        if high - low < 1e-9:   # flat: draw the line across the middle
            return low - 1.0, high + 1.0
        return low, high

    def _y(self, value: float, low: float, high: float) -> float:
        plot = self._plot()
        return plot.bottom() - (value - low) / (high - low) * plot.height()

    @override
    def paintEvent(self, event: QPaintEvent) -> None:
        if self.theme is None or len(self.values) < 2:
            return
        theme = self.theme
        values = self.values
        low, high = self._range()
        plot = self._plot()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        font = QFont(self.font())
        font.setPixelSize(9)
        p.setFont(font)
        p.setPen(theme.muted)
        corner = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        p.drawText(QRectF(plot.left(), 0, plot.width(), 12), corner, self.caption)
        grid = QColor(theme.text)
        grid.setAlphaF(0.12)
        p.setPen(QPen(grid, 1, Qt.PenStyle.DashLine))
        p.drawLine(QPointF(plot.left(), plot.top()), QPointF(plot.right(), plot.top()))
        grid.setAlphaF(0.25)
        p.setPen(QPen(grid, 1))
        p.drawLine(QPointF(plot.left(), plot.bottom()), QPointF(plot.right(), plot.bottom()))
        points = [QPointF(self._x(k), self._y(float(v), low, high)) for k, v in enumerate(values)]
        curve = QPainterPath(points[0])
        for point in points[1:]:
            curve.lineTo(point)
        area = QPainterPath(curve)
        area.lineTo(QPointF(points[-1].x(), plot.bottom()))
        area.lineTo(QPointF(points[0].x(), plot.bottom()))
        area.closeSubpath()
        fill = QColor(theme.accent)
        fill.setAlphaF(0.12)
        p.fillPath(area, fill)
        p.setPen(QPen(theme.text, 1.3))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawPath(curve)
        if len(points) <= DOTS_UP_TO:
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(theme.text)
            for point in points:
                p.drawEllipse(point, 1.6, 1.6)
        position = min(max(self.position, 0.0), len(values) - 1.0)
        k = min(int(position), len(values) - 2)
        t = position - k
        value = float(values[k]) * (1.0 - t) + float(values[k + 1]) * t
        x = self._x(position)
        p.setPen(QPen(theme.accent, 1.5))
        p.drawLine(QPointF(x, plot.top()), QPointF(x, plot.bottom()))
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(theme.accent)
        p.drawEllipse(QPointF(x, self._y(value, low, high)), 3.5, 3.5)
        p.end()

    def _frame_at(self, x: float) -> int:
        plot = self._plot()
        fraction = (x - plot.left()) / max(plot.width(), 1.0)
        return round(min(max(fraction, 0.0), 1.0) * (len(self.values) - 1))

    @override
    def mousePressEvent(self, e: QMouseEvent) -> None:
        if len(self.values) > 1 and e.button() == Qt.MouseButton.LeftButton:
            self._last = self._frame_at(e.position().x())
            self.seek.emit(self._last)

    @override
    def mouseMoveEvent(self, e: QMouseEvent) -> None:
        if len(self.values) > 1 and e.buttons() & Qt.MouseButton.LeftButton:
            frame = self._frame_at(e.position().x())
            if frame != self._last:
                self._last = frame
                self.seek.emit(frame)


class Timeline(QFrame):
    """Play button, frame slider, speed and frame label, shown when the document has frames.

    Playing advances a fractional frame position at the chosen speed and
    blends neighbouring frames, so motion stays smooth even when frames are
    far apart. A plot above the slider shows one quantity per frame: the
    energy relative to the lowest, a pinned distance, angle or dihedral, or
    the RMSD from the first frame. It appears by itself when there are
    energies or pinned measurements; the plot button shows or hides it.
    """

    relayout = Signal()
    exportRequested = Signal(str)   # one of the kinds in miles.ui.plot_export

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("Panel")
        self.doc: Document | None = None
        self.theme: Theme | None = None
        self._position = 0.0
        self._clock = 0.0
        self._quantities: list[Quantity] = []
        self._series: dict[str, np.ndarray] = {}    # raw values per quantity key, for self.doc
        self._plotted: tuple[str, int] | None = None
        self._shown: weakref.WeakKeyDictionary[Document, bool] = weakref.WeakKeyDictionary()
        self._chosen: weakref.WeakKeyDictionary[Document, str] = weakref.WeakKeyDictionary()
        lay = QGridLayout(self)
        lay.setContentsMargins(8, 4, 12, 4)
        lay.setHorizontalSpacing(10)
        lay.setVerticalSpacing(0)
        self.plot = FramePlot()
        self.plot.seek.connect(self._seek)
        self.plot.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.plot.customContextMenuRequested.connect(
            lambda at: show_plot_menu(self.plot, at, self.menu_entries()))
        self.plot.hide()
        lay.addWidget(self.plot, 0, 1)
        self.choice = QComboBox()
        self.choice.setToolTip("Quantity plotted along the frames")
        self.choice.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        self.choice.currentIndexChanged.connect(self._chose)
        self.choice.hide()
        lay.addWidget(self.choice, 0, 2, 1, 3, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.play_button = QToolButton()
        self.play_button.setObjectName("Order")
        self.play_button.setFixedSize(30, 28)
        self.play_button.setIconSize(QSize(20, 20))
        self.play_button.setToolTip("Play or pause (Space)")
        self.play_button.clicked.connect(self.toggle_playing)
        lay.addWidget(self.play_button, 1, 0)
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setObjectName("Timeline")
        self.slider.setMinimumWidth(280)
        self.slider.valueChanged.connect(self._slid)
        lay.addWidget(self.slider, 1, 1)
        self.text = label("", "Muted")
        self.text.setMinimumWidth(TEXT_WIDTH)
        lay.addWidget(self.text, 1, 2)
        self.speed = Segmented([name for name, _ in SPEEDS])
        self.speed.set_index(2)
        self.speed.setToolTip("Playback speed")
        lay.addWidget(self.speed, 1, 3)
        self.plot_button = QToolButton()
        self.plot_button.setObjectName("Order")
        self.plot_button.setCheckable(True)
        self.plot_button.setFixedSize(30, 28)
        self.plot_button.setIconSize(QSize(18, 18))
        self.plot_button.setToolTip("Plot energy, pinned measurements or RMSD along the frames")
        self.plot_button.clicked.connect(self._toggle_plot)
        lay.addWidget(self.plot_button, 1, 4)
        lay.setColumnStretch(1, 1)
        self._timer = QTimer(self)
        self._timer.setInterval(TICK_MS)
        self._timer.timeout.connect(self._advance)
        self.hide()

    @property
    def rate(self) -> float:
        """Frames per second at the chosen speed."""
        return BASE_RATE * SPEEDS[self.speed.index()][1]

    @property
    def playing(self) -> bool:
        """Whether the frames are playing."""
        return self._timer.isActive()

    @property
    def quantity(self) -> Quantity | None:
        """The quantity chosen for the plot (whether or not the plot is shown)."""
        row = self.choice.currentIndex()
        return self._quantities[row] if 0 <= row < len(self._quantities) else None

    def menu_entries(self) -> list[tuple[str, Callable[[], None]]]:
        """Entries of the plot's right-click menu."""
        return [("Copy Data", lambda: self.exportRequested.emit(COPY)),
                ("Export Data...", lambda: self.exportRequested.emit(DATA))]

    def export_name(self) -> str:
        """File name to suggest, e.g. ``alanine_dipeptide_md_frames``."""
        name = os.path.splitext(self.doc.name)[0] if self.doc is not None else ""
        return f"{name or 'trajectory'}_frames"

    def export_table(self, kind: str) -> Table | None:
        """Tabulate every quantity of the plot (energies, pinned measurements, RMSD) for all frames."""
        doc = self.doc
        if doc is None or doc.frames is None or doc.frame_count < 2:
            return None
        measurements = [(q.name, q.atoms) for q in quantities(doc) if q.atoms]
        return frame_table(doc.frames, doc.mol.symbols, doc.frame_energies, doc.frame_labels, measurements)

    def plot_widget(self) -> None:
        """Return nothing: the timeline plot is not saved as a picture."""

    def set_theme(self, theme: Theme) -> None:
        """Recolour the buttons and the plot."""
        self.theme = theme
        self.plot.theme = theme
        self.plot.update()
        self.plot_button.setIcon(icons.icon("plot", theme.text))
        self._update_icon()

    def set_document(self, doc: Document | None) -> None:
        """Follow another document."""
        if doc is not self.doc:
            self._timer.stop()
            self._series = {}
            self._plotted = None
            self.text.setMinimumWidth(TEXT_WIDTH)
        self.doc = doc
        self.refresh()

    def doc_changed(self, kind: str) -> None:
        """Follow frame changes and edits made elsewhere (keyboard, conformer search, measurements)."""
        if kind not in (Change.FRAME, Change.STYLE):   # coordinates or frames changed: measure again
            self._series = {}
            self._plotted = None
        self.refresh()

    def plot_visible(self) -> bool:
        """Whether the plot should show for the current document."""
        doc = self.doc
        if doc is None or doc.frame_count < 2:
            return False
        explicit = self._shown.get(doc)
        if explicit is not None:
            return explicit
        return doc.frame_energies is not None or any(2 <= len(m) <= 4 for m in doc.mol.measurements)

    def refresh(self) -> None:
        """Show the current frame, or hide when there is no trajectory."""
        doc = self.doc
        count = doc.frame_count if doc is not None else 0
        plotted = False
        if doc is None or count < 2:
            self._timer.stop()
            visible = False
        else:
            self.slider.blockSignals(True)
            self.slider.setRange(0, count - 1)
            self.slider.setValue(doc.frame)
            self.slider.blockSignals(False)
            plotted = self.plot_visible()
            self._fill_choice(doc)
            text = f"Frame {doc.frame + 1} of {count}"
            if doc.frame_labels and doc.frame < len(doc.frame_labels) and doc.frame_labels[doc.frame]:
                text += SEPARATOR + doc.frame_labels[doc.frame]
            if doc.frame_energies is not None:
                text += SEPARATOR + relative_energy_text(doc.frame_energies, doc.frame)
            quantity = self.quantity
            if plotted and quantity is not None:
                self._update_plot(doc, quantity)
                if quantity.key != ENERGY:
                    value = float(self._values(doc, quantity)[doc.frame])
                    text += f"{SEPARATOR}{quantity.short} {value:.{quantity.decimals}f}{quantity.unit}"
                self.plot.set_position(self._position if self.playing else float(doc.frame))
            self.text.setText(text)
            wide = self.text.sizeHint().width()
            if wide > self.text.minimumWidth():   # only grow, so the slider does not jitter while playing
                self.text.setMinimumWidth(wide)
            visible = True
        self.plot_button.setChecked(plotted)
        changed = False
        if plotted != self.plot.isVisibleTo(self):
            self.plot.setVisible(plotted)
            changed = True
        self.choice.setVisible(plotted and len(self._quantities) > 1)
        if visible == self.isHidden():   # isHidden: right before the window is shown, too
            self.setVisible(visible)
            changed = True
        if changed:
            self.adjustSize()
            self.relayout.emit()
        self._update_icon()

    def _fill_choice(self, doc: Document) -> None:
        """List what can be plotted, keeping the user's choice (or the first quantity)."""
        found = quantities(doc)
        if [q.key for q in found] != [q.key for q in self._quantities] or self.choice.count() != len(found):
            self._quantities = found
            self.choice.blockSignals(True)
            self.choice.clear()
            self.choice.addItems([q.name for q in found])
            self.choice.blockSignals(False)
        keys = [q.key for q in found]
        wanted = self._chosen.get(doc)
        row = keys.index(wanted) if wanted in keys else 0
        if self.choice.currentIndex() != row:
            self.choice.blockSignals(True)
            self.choice.setCurrentIndex(row)
            self.choice.blockSignals(False)

    def _values(self, doc: Document, quantity: Quantity) -> np.ndarray:
        """Raw values of a quantity in every frame, computed once per set of frames."""
        values = self._series.get(quantity.key)
        if values is None:
            assert doc.frames is not None
            if quantity.key == ENERGY:
                assert doc.frame_energies is not None
                values = np.asarray(doc.frame_energies, dtype=float)
            elif quantity.key == RMSD:
                values = rmsd_series(doc.frames, rmsd_atoms(doc.mol.symbols))
            else:
                values = measurement_series(doc.frames, quantity.atoms)
            self._series[quantity.key] = values
        return values

    def _update_plot(self, doc: Document, quantity: Quantity) -> None:
        values = self._values(doc, quantity)
        if self._plotted == (quantity.key, id(values)):
            return
        self._plotted = (quantity.key, id(values))
        if quantity.key == ENERGY:
            self.plot.set_energies(list(values))
            return
        if len(quantity.atoms) == 4:
            values = unwrap_degrees(values)   # a dihedral turning past 180 degrees stays one line
        d = quantity.decimals
        caption = f"{values.min():.{d}f} to {values.max():.{d}f}{quantity.unit}"
        if len(self._quantities) == 1:   # no chooser to name it
            caption = f"{quantity.name}  {caption}"
        self.plot.set_values(values, caption, from_zero=False)

    def _chose(self, row: int) -> None:
        if self.doc is not None and 0 <= row < len(self._quantities):
            self._chosen[self.doc] = self._quantities[row].key
            self.text.setMinimumWidth(TEXT_WIDTH)
            self.refresh()

    def _toggle_plot(self, checked: bool) -> None:
        if self.doc is not None:
            self._shown[self.doc] = checked
            self.text.setMinimumWidth(TEXT_WIDTH)
            self.refresh()

    def play(self) -> None:
        """Start playing the frames in a loop (no change if already playing)."""
        if not self.playing and self.doc is not None and self.doc.frame_count > 1:
            self._position = float(self.doc.frame)
            self._clock = time.monotonic()
            self._timer.start()
            self._update_icon()

    def pause(self) -> None:
        """Stop on the nearest whole frame."""
        if self.playing:
            self._timer.stop()
            if self.doc is not None:
                self.doc.set_frame(round(self._position))
            self._update_icon()

    def toggle_playing(self) -> None:
        """Start or stop playing."""
        if self.playing:
            self.pause()
        else:
            self.play()

    def step(self, delta: int = 1) -> None:
        """Show the next (or previous) frame, wrapping around."""
        if self.doc is not None and self.doc.frame_count > 1:
            self.doc.set_frame(self.doc.frame + delta)

    def _advance(self) -> None:
        if self.doc is None or self.doc.frame_count < 2:
            self._timer.stop()
            return
        now = time.monotonic()
        self._position = (self._position + (now - self._clock) * self.rate) % self.doc.frame_count
        self._clock = now
        self.doc.set_position(self._position)

    def _slid(self, value: int) -> None:
        self._seek(value)

    def _seek(self, frame: int) -> None:
        if self.doc is not None:
            self._position = float(frame)
            self._clock = time.monotonic()
            self.doc.set_frame(frame)

    def _update_icon(self) -> None:
        if self.theme is not None:
            name = "pause" if self.playing else "play"
            self.play_button.setIcon(icons.icon(name, self.theme.text))
