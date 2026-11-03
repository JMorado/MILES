"""A picture of a tab's view that notebooks show inline: :class:`Picture`."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from PySide6.QtGui import QImage



class Picture:
    """A picture of a view, made by :meth:`DocumentView.image`.

    IPython and Jupyter show it inline. It has :attr:`DENSITY` pixels for
    each point of its shown size, so it is sharp on high-resolution
    screens.

    Parameters
    ----------
    image : QImage
        The rendered picture.
    """

    DENSITY = 2.0
    """Pixels per point of the shown size."""

    def __init__(self, image: "QImage") -> None:
        self.image = image

    @property
    def width(self) -> int:
        """Width in pixels."""
        return self.image.width()

    @property
    def height(self) -> int:
        """Height in pixels."""
        return self.image.height()

    def png(self) -> bytes:
        """Return the picture as PNG data."""
        from PySide6.QtCore import QBuffer, QByteArray, QIODevice
        from PySide6.QtGui import QImageWriter

        data = QByteArray()
        buffer = QBuffer(data)
        buffer.open(QIODevice.OpenModeFlag.WriteOnly)
        QImageWriter(buffer, QByteArray(b"PNG")).write(self.image)
        buffer.close()
        return bytes(data.data())

    def save(self, path: str) -> None:
        """Write the picture to a file (PNG, or the format the extension names).

        Raises
        ------
        OSError
            If the file cannot be written.
        """
        if not self.image.save(path):
            raise OSError(f"Could not write {path}")

    def _repr_png_(self) -> tuple[bytes, dict[str, int]]:
        size = {"width": round(self.width / self.DENSITY), "height": round(self.height / self.DENSITY)}
        return self.png(), size

    def __repr__(self) -> str:
        return f"<picture, {self.width} x {self.height} pixels>"
