"""Floating panels: folding to the top row and moving them by dragging."""

from pathlib import Path

import pytest
from PySide6.QtCore import QPoint, QPointF, QSettings, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QApplication

from miles.ui.mainwindow import MainWindow

DATA = Path(__file__).parent / "data"

def grid_panel(window):
    """The grid panel of the Orbitals and Grids plugin."""
    return window.stage.floating_panels()["grid"]


def isosurface(canvas):
    """The isosurface mesh drawn on a canvas, or None."""
    items = canvas.drawing.scene_items()
    return items.meshes[0].mesh if items.meshes else None



@pytest.fixture
def window(app, tmp_path):
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    win.resize(1100, 760)
    win.show()
    app.processEvents()
    yield win
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


def mouse(widget, kind, at, buttons=Qt.MouseButton.LeftButton):
    """Send a mouse event at a point of ``widget``."""
    button = Qt.MouseButton.NoButton if kind == QMouseEvent.Type.MouseMove else Qt.MouseButton.LeftButton
    event = QMouseEvent(kind, QPointF(at), widget.mapToGlobal(QPointF(at)), button, buttons,
                        Qt.KeyboardModifier.NoModifier)
    QApplication.sendEvent(widget, event)


def test_panels_fold_to_their_top_row(window):
    window.open_path(str(DATA / "water_freq.log"))
    panel = window.stage.floating_panels()["vibrations"]
    assert panel.isVisible() and panel.body.isVisible()
    open_height = panel.height()
    panel.chevron.click()
    assert panel.collapsed and not panel.body.isVisible()
    assert panel.height() < open_height / 3
    panel.chevron.click()
    assert not panel.collapsed and panel.body.isVisible()
    assert panel.height() == open_height


def test_panels_move_by_dragging_and_stay_inside_the_window(window):
    window.open_path(str(DATA / "water_freq.log"))
    stage, panel = window.stage, window.stage.floating_panels()["vibrations"]
    start = panel.pos()
    grip = QPoint(6, 6)   # the panel's margin, not a control
    mouse(panel, QMouseEvent.Type.MouseButtonPress, grip)
    mouse(panel, QMouseEvent.Type.MouseMove, grip + QPoint(-300, 120))
    mouse(panel, QMouseEvent.Type.MouseButtonRelease, grip + QPoint(-300, 120), Qt.MouseButton.NoButton)
    assert panel.pos() == start + QPoint(-300, 120) and panel.placed == panel.pos()
    stage.layout_overlays()   # the stage keeps the chosen place
    assert panel.pos() == start + QPoint(-300, 120)

    mouse(panel, QMouseEvent.Type.MouseButtonPress, grip)
    mouse(panel, QMouseEvent.Type.MouseMove, grip + QPoint(5000, 5000))   # far outside: kept in the window
    mouse(panel, QMouseEvent.Type.MouseButtonRelease, grip, Qt.MouseButton.NoButton)
    assert panel.geometry().right() < stage.width() and panel.geometry().bottom() < stage.height()

    QApplication.sendEvent(panel, QMouseEvent(QMouseEvent.Type.MouseButtonDblClick, QPointF(grip),
                                              panel.mapToGlobal(QPointF(grip)), Qt.MouseButton.LeftButton,
                                              Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier))
    assert panel.placed is None and panel.pos() == start   # back to its usual place


def test_grid_panel_folds_to_its_controls(window):
    window.open_path(str(DATA / "h2_sto3g.molden"))
    panel = grid_panel(window)
    assert panel.dial.isVisible()
    panel.chevron.click()
    assert not panel.dial.isVisible() and panel.choice.isVisible() and panel.level.isVisible()


def test_panel_places_and_folds_are_remembered(app, tmp_path):
    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    first = MainWindow(settings)
    first.resize(1100, 760)
    first.show()
    first.open_path(str(DATA / "water_freq.log"))
    panel = first.stage.floating_panels()["vibrations"]
    grip = QPoint(6, 6)
    mouse(panel, QMouseEvent.Type.MouseButtonPress, grip)
    mouse(panel, QMouseEvent.Type.MouseMove, grip + QPoint(-250, 80))
    mouse(panel, QMouseEvent.Type.MouseButtonRelease, grip, Qt.MouseButton.NoButton)
    moved = panel.pos()
    first.stage.floating_panels()["properties"].chevron.click()   # folded while hidden: remembered too
    for canvas in first.canvases:
        canvas.doc.modified = False
    first.close()

    second = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    second.resize(1100, 760)
    second.show()
    second.open_path(str(DATA / "water_freq.log"))
    panels = second.stage.floating_panels()
    assert panels["vibrations"].placed == moved and panels["vibrations"].pos() == moved
    assert panels["properties"].collapsed and not panels["vibrations"].collapsed
    for canvas in second.canvases:
        canvas.doc.modified = False
    second.close()


def test_fitting_keeps_the_structure_clear_of_bottom_panels(window):
    from miles.core.chem import from_smiles
    from miles.files.common import JobInfo

    window.show_console(True)
    mol, _ = from_smiles("c1ccc2ccccc2c1")
    window.open_structure(mol, JobInfo(title="naphthalene"))
    canvas = window.current_canvas()
    canvas.fit()
    QApplication.processEvents()
    _, sy = canvas.camera.to_screen(canvas.camera.to_view(canvas.mol.coords), canvas.width(), canvas.height())
    console_top = window.stage.console.geometry().top()
    assert canvas.view_insets[1] > 0
    assert sy.max() < console_top
