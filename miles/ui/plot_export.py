"""Copying and exporting the data and pictures of the plots in the panels."""

from collections.abc import Callable, Sequence
from typing import Protocol

from PySide6.QtCore import QByteArray, QPoint, QRect, QSize, Qt
from PySide6.QtGui import QColor, QImage, QImageWriter, QPainter
from PySide6.QtSvg import QSvgGenerator
from PySide6.QtWidgets import QMenu, QWidget

from miles.core.tables import Table
from miles.render.theme import LIGHT, Theme

COPY, DATA, CURVE, IMAGE = "copy", "data", "curve", "image"
"""What a plot's menu asks for: copy its table, save it as CSV, save the spectrum curve, save a picture."""

IMAGE_SCALE = 3.0   # pixels per point of the plot in exported PNG files
TITLE_HEIGHT = 18   # points added under an exported plot for its axis title


class PaintablePlot(Protocol):
    """A plot widget that can paint itself onto any painter."""

    def size(self) -> QSize:
        """Size in points."""
        ...

    def paint(self, p: QPainter, theme: Theme, for_export: bool = False) -> None:
        """Draw the plot; ``for_export`` leaves out hover and selection marks."""
        ...


class PlotSource(Protocol):
    """A panel whose plot can be copied and exported."""

    def export_name(self) -> str:
        """File name, without extension, to suggest when saving."""
        ...

    def export_table(self, kind: str) -> Table | None:
        """Table for :data:`COPY`, :data:`DATA` or :data:`CURVE`, or None when there is nothing."""
        ...

    def plot_widget(self) -> PaintablePlot | None:
        """Return the plot to save as a picture, or None when the panel has none."""
        ...


def show_plot_menu(parent: QWidget, at: QPoint, entries: Sequence[tuple[str, Callable[[], None]]]) -> None:
    """Open a context menu of plot actions at a point of ``parent``.

    Parameters
    ----------
    parent : QWidget
        Widget that was right-clicked.
    at : QPoint
        Position in ``parent``.
    entries : sequence of (str, callable)
        Menu text and action of each entry.
    """
    menu = QMenu(parent)
    for text, action in entries:
        menu.addAction(text).triggered.connect(action)
    menu.exec(parent.mapToGlobal(at))
    menu.deleteLater()


def save_plot_image(plot: PaintablePlot, path: str, theme: Theme = LIGHT) -> None:
    """Save a plot as an SVG drawing or a PNG picture, on white.

    Parameters
    ----------
    plot : PaintablePlot
        The plot.
    path : str
        File to write; ``.svg`` gives a drawing, anything else a PNG at
        :data:`IMAGE_SCALE` pixels per point.
    theme : Theme, optional
        Colours of the lines and labels (a light theme, for white paper).

    Raises
    ------
    OSError
        If the file cannot be written.
    """
    plot_size = plot.size()
    title_of = getattr(plot, "axis_title", None)
    title = title_of() if callable(title_of) else ""
    size = QSize(plot_size.width(), plot_size.height() + (TITLE_HEIGHT if title else 0))

    def draw(p: QPainter) -> None:
        plot.paint(p, theme, for_export=True)
        if title:
            font = p.font()
            font.setPixelSize(10)
            p.setFont(font)
            p.setPen(theme.muted)
            below = QRect(0, plot_size.height(), size.width(), TITLE_HEIGHT)
            p.drawText(below, Qt.AlignmentFlag.AlignCenter, title)

    if path.lower().endswith(".svg"):
        generator = QSvgGenerator()
        generator.setFileName(path)
        generator.setSize(size)
        generator.setViewBox(QRect(QPoint(0, 0), size))
        p = QPainter(generator)
        p.fillRect(QRect(QPoint(0, 0), size), QColor(Qt.GlobalColor.white))
        draw(p)
        p.end()
        return
    image = QImage(round(size.width() * IMAGE_SCALE), round(size.height() * IMAGE_SCALE),
                   QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(QColor(Qt.GlobalColor.white))
    p = QPainter(image)
    p.scale(IMAGE_SCALE, IMAGE_SCALE)
    draw(p)
    p.end()
    writer = QImageWriter(path, QByteArray(b"png"))   # PNG whatever the name ends with
    if not writer.write(image):
        raise OSError(f"Could not write {path}: {writer.errorString()}")
