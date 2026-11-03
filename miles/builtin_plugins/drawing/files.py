"""2D structure drawings as images: PNG rendering, files and the clipboard."""

from PySide6.QtCore import QByteArray, QMimeData, QRectF
from PySide6.QtGui import QColor, QImage, QPainter
from PySide6.QtSvg import QSvgRenderer

SVG_MIME = "image/svg+xml"
PNG_SCALE = 4.0   # pixels per SVG unit; a typical drug molecule comes out 1000 to 1500 pixels wide


def svg_image(svg: str, scale: float = PNG_SCALE, background: QColor | None = None) -> QImage:
    """Render an SVG drawing to an image.

    Parameters
    ----------
    svg : str
        SVG document, e.g. from :func:`~miles.builtin_plugins.drawing.depiction.depiction_svg`.
    scale : float, optional
        Pixels per SVG unit.
    background : QColor, optional
        Fill colour; None keeps the background transparent.

    Returns
    -------
    QImage
        The picture (empty if the SVG cannot be read).
    """
    renderer = QSvgRenderer(QByteArray(svg.encode()))
    size = renderer.defaultSize()
    if not renderer.isValid() or size.isEmpty():
        return QImage()
    image = QImage(max(1, round(size.width() * scale)), max(1, round(size.height() * scale)),
                   QImage.Format.Format_ARGB32_Premultiplied)
    image.fill(background if background is not None else QColor(0, 0, 0, 0))
    painter = QPainter(image)
    painter.setRenderHints(QPainter.RenderHint.Antialiasing | QPainter.RenderHint.TextAntialiasing)
    renderer.render(painter, QRectF(0, 0, image.width(), image.height()))
    painter.end()
    return image


def write_drawing(svg: str, path: str) -> None:
    """Save a drawing as SVG, or as a PNG on white when the path ends in ``.png``.

    Raises
    ------
    OSError
        If the file cannot be written.
    """
    if path.lower().endswith(".png"):
        if not svg_image(svg, background=QColor("white")).save(path):   # format from the suffix
            raise OSError(f"Could not write {path}")
        return
    with open(path, "w", encoding="utf-8") as f:
        f.write(svg)


def drawing_mime(svg: str) -> QMimeData:
    """Clipboard contents for a drawing: a PNG image on white, and the SVG for programs that take it."""
    mime = QMimeData()
    mime.setImageData(svg_image(svg, background=QColor("white")))
    mime.setData(SVG_MIME, QByteArray(svg.encode()))
    return mime
