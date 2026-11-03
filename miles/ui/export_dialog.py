"""Publication picture export: size, background, ambient occlusion and cropping, with a preview."""

import os

import numpy as np
from PySide6.QtCore import QRect, QSettings, Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPixmap
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QFileDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from miles.render.theme import Background, Theme
from miles.ui.canvas import Canvas
from miles.ui.controls import Segmented, dialog_stylesheet
from miles.ui.dialogs import muted

SCALES = [("1x", 1.0), ("2x", 2.0), ("4x", 4.0)]
BACKGROUNDS = [("Theme", Background.THEME), ("White", Background.WHITE),
               ("Transparent", Background.TRANSPARENT)]
PREVIEW_SIZE = (600, 380)
CHECKER = 10   # checkerboard square behind transparent previews, in pixels


def crop_image(image: QImage, transparent: bool, margin: float = 0.04) -> QImage:
    """Trim the empty border around the molecule, keeping a small margin.

    Parameters
    ----------
    image : QImage
        Rendered picture.
    transparent : bool
        Find the molecule by alpha instead of by difference from the corner colour.
    margin : float, optional
        Margin as a fraction of the molecule's larger dimension.

    Returns
    -------
    QImage
        The cropped picture (unchanged if it is empty).
    """
    rgba = image.convertToFormat(QImage.Format.Format_RGBA8888)
    w, h = rgba.width(), rgba.height()
    raw = np.frombuffer(rgba.constBits(), np.uint8, count=rgba.sizeInBytes())
    pixels = raw.reshape(h, rgba.bytesPerLine())[:, :w * 4].reshape(h, w, 4)
    if transparent:
        mask = pixels[..., 3] > 8
    else:
        background = pixels[0, 0, :3].astype(int)
        mask = np.abs(pixels[..., :3].astype(int) - background).sum(axis=2) > 12
    ys, xs = np.nonzero(mask)
    if not len(xs):
        return image
    left, right, top, bottom = int(xs.min()), int(xs.max()), int(ys.min()), int(ys.max())
    pad = int(margin * max(right - left, bottom - top)) + 2
    x0, y0 = max(0, left - pad), max(0, top - pad)
    x1, y1 = min(w, right + pad + 1), min(h, bottom + pad + 1)
    return image.copy(QRect(x0, y0, x1 - x0, y1 - y0))


class ExportDialog(QDialog):
    """Options for a publication picture of the current view, with a live preview."""

    def __init__(self, canvas: Canvas, theme: Theme, settings: QSettings,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.canvas = canvas
        self.settings = settings
        self.saved_path = ""
        self.setWindowTitle("Export Image")
        self.setStyleSheet(dialog_stylesheet(theme) + f"""
            QPushButton#Plain {{ border: none; color: {theme.muted.name()}; padding: 6px 12px; }}
            QPushButton#Plain:hover {{ color: {theme.text.name()}; }}""")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(18, 16, 18, 16)
        lay.setSpacing(10)
        self.preview = QLabel()
        self.preview.setFixedSize(*PREVIEW_SIZE)
        self.preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
        lay.addWidget(self.preview)

        self.scale = Segmented([name for name, _ in SCALES])
        self.scale.set_index(int(str(settings.value("export/scale", 1))))
        self.dims = muted(QLabel())
        self.background = Segmented([name for name, _ in BACKGROUNDS])
        self.background.set_index(int(str(settings.value("export/background", 1))))
        self.ao = QCheckBox("Ambient occlusion")
        self.ao.setChecked(str(settings.value("export/ao", "true")) == "true")
        self.crop = QCheckBox("Crop to molecule")
        self.crop.setChecked(str(settings.value("export/crop", "true")) == "true")
        grid = QGridLayout()
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(8)
        grid.addWidget(muted(QLabel("Size")), 0, 0)
        grid.addWidget(self.scale, 0, 1)
        grid.addWidget(self.dims, 0, 2)
        grid.addWidget(muted(QLabel("Background")), 1, 0)
        grid.addWidget(self.background, 1, 1, 1, 2)
        grid.addWidget(muted(QLabel("Options")), 2, 0)
        grid.addWidget(self.ao, 2, 1)
        grid.addWidget(self.crop, 2, 2)
        grid.setColumnStretch(3, 1)
        lay.addLayout(grid)

        buttons = QHBoxLayout()
        buttons.addStretch(1)
        cancel = QPushButton("Cancel")
        cancel.setObjectName("Plain")
        cancel.clicked.connect(self.reject)
        save = QPushButton("Save PNG...")
        save.setObjectName("Primary")
        save.setDefault(True)
        save.clicked.connect(self._save)
        buttons.addWidget(cancel)
        buttons.addWidget(save)
        lay.addLayout(buttons)

        self.scale.changed.connect(self._refresh)
        self.background.changed.connect(self._refresh)
        self.ao.toggled.connect(self._refresh)
        self.crop.toggled.connect(self._refresh)
        self._refresh()

    @property
    def _scale(self) -> float:
        return SCALES[self.scale.index()][1]

    @property
    def _background(self) -> Background:
        return BACKGROUNDS[self.background.index()][1]

    def _render(self, scale: float) -> QImage:
        image = self.canvas.render_export(scale, self._background, self.ao.isChecked())
        if self.crop.isChecked():
            return crop_image(image, self._background == Background.TRANSPARENT)
        return image

    def _refresh(self) -> None:
        image = self._render(1.0)
        self.dims.setText(f"{int(image.width() * self._scale)} x {int(image.height() * self._scale)} px")
        box = self.preview.size()
        canvas = QPixmap(box)
        canvas.fill(Qt.GlobalColor.transparent)
        p = QPainter(canvas)
        p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        pm = QPixmap.fromImage(image).scaled(box, Qt.AspectRatioMode.KeepAspectRatio,
                                             Qt.TransformationMode.SmoothTransformation)
        x, y = (box.width() - pm.width()) // 2, (box.height() - pm.height()) // 2
        if self._background == Background.TRANSPARENT:
            for i in range(0, pm.width(), CHECKER):
                for j in range(0, pm.height(), CHECKER):
                    shade = 236 if (i // CHECKER + j // CHECKER) % 2 else 250
                    p.fillRect(x + i, y + j, min(CHECKER, pm.width() - i), min(CHECKER, pm.height() - j),
                               QColor(shade, shade, shade))
        p.drawPixmap(x, y, pm)
        p.end()
        self.preview.setPixmap(canvas)

    def _save(self) -> None:
        self.settings.setValue("export/scale", self.scale.index())
        self.settings.setValue("export/background", self.background.index())
        self.settings.setValue("export/ao", "true" if self.ao.isChecked() else "false")
        self.settings.setValue("export/crop", "true" if self.crop.isChecked() else "false")
        doc = self.canvas.doc
        folder = str(self.settings.value("last_dir", os.path.expanduser("~")))
        base = os.path.splitext(doc.path)[0] if doc.path else os.path.join(folder, "molecule")
        path, _ = QFileDialog.getSaveFileName(self, "Export Image", base + ".png", "PNG image (*.png)")
        if not path:
            return
        if not path.lower().endswith(".png"):
            path += ".png"
        if self._render(self._scale).save(path):
            self.saved_path = path
            self.accept()
