"""The application icon is drawn from the logo at every size."""

from PySide6.QtCore import QSize
from PySide6.QtWidgets import QApplication

from miles.ui import icons


def test_logo_is_packaged_next_to_the_code() -> None:
    assert icons.APP_ICON.is_file()


def test_icon_has_every_size_drawn(app: QApplication) -> None:
    icon = icons.app_icon()
    assert {s.width() for s in icon.availableSizes()} == set(icons.APP_ICON_SIZES)
    for size in (16, 64, 512):
        image = icon.pixmap(QSize(size, size)).toImage()
        assert image.width() == size
        assert image.pixelColor(0, 0).alpha() == 0                  # transparent outside the rounded plate
        plate = image.pixelColor(size // 2, int(size * 0.72))       # navy plate below the M's middle atom
        assert plate.alpha() == 255 and plate.blue() > plate.red()


def test_window_uses_the_logo(app: QApplication, tmp_path) -> None:
    from PySide6.QtCore import QSettings

    from miles.ui.mainwindow import MainWindow

    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    try:
        assert {s.width() for s in win.windowIcon().availableSizes()} == set(icons.APP_ICON_SIZES)
    finally:
        for canvas in win.canvases:
            canvas.doc.modified = False
        win.close()


def test_application_has_the_logo_for_the_dock(app: QApplication) -> None:
    from miles.ui.app import application

    assert application() is app
    assert {s.width() for s in app.windowIcon().availableSizes()} == set(icons.APP_ICON_SIZES)
