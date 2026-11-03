"""Arrows showing how each atom moves in a vibrational mode."""

from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QSettings

from miles.builtin_plugins.vibrations.modes import mode_arrows
from miles.core.vibrations import VibrationalMode
from miles.render.theme import Background
from miles.ui.mainwindow import MainWindow

DATA = Path(__file__).parent / "data"


@pytest.fixture
def window(app, tmp_path):
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    win.resize(1000, 700)
    win.show()
    app.processEvents()
    yield win
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


def test_mode_arrows_are_scaled_and_small_ones_dropped():
    displacements = np.array([[0.0, 0.0, 0.5], [0.0, 0.2, 0.0], [0.01, 0.0, 0.0]])
    arrows = mode_arrows(VibrationalMode(1000.0, 1.0, displacements), length=1.0, threshold=0.05)
    assert np.allclose(np.linalg.norm(arrows, axis=1), [1.0, 0.4, 0.0])
    assert np.allclose(arrows[0], [0.0, 0.0, 1.0])
    still = mode_arrows(VibrationalMode(10.0, None, np.zeros((2, 3))))
    assert still.shape == (2, 3) and not still.any()


def test_arrows_toggle_with_and_without_animation(window):
    window.open_path(str(DATA / "water_freq.log"))
    canvas = window.current_canvas()
    panel = window.stage.floating_panels()["vibrations"]
    doc = window.api.current
    assert panel.isVisible() and not panel.arrows.isChecked()
    panel.arrows.click()   # nothing chosen yet: the first mode
    assert panel.arrows_for(doc).shape == (3, 3)
    assert np.isclose(np.linalg.norm(panel.arrows_for(doc), axis=1).max(), 1.1)
    assert panel.spectrum.selected == 0
    canvas.grab()   # draws the arrows in the overlay

    panel.animate(2)
    expected = mode_arrows(canvas.doc.info.vibrations[2])
    assert np.allclose(panel.arrows_for(doc), expected)
    canvas.grab()
    panel.stop_animation()
    assert panel.arrows_for(doc) is not None and panel.spectrum.selected == 2   # still shown at rest

    panel.arrows.click()
    assert panel.arrows_for(doc) is None and len(canvas.drawing.scene_items().arrows) == 0


def test_arrows_belong_to_their_document(window):
    window.open_path(str(DATA / "water_freq.log"))
    first = window.current_canvas()
    first_view = window.api.current
    panel = window.stage.floating_panels()["vibrations"]
    panel.arrows.click()
    window.new_document()
    assert panel.arrows_for(window.api.current) is None
    window.open_path(str(DATA / "water_freq.log"))
    assert not panel.arrows.isChecked()
    window.tabbar.setCurrentIndex(window.canvases.index(first))
    assert panel.arrows.isChecked()
    panel.close_panel()
    assert panel.arrows_for(first_view) is None


def test_arrows_appear_in_exported_pictures(window):
    window.open_path(str(DATA / "water_freq.log"))
    canvas = window.current_canvas()
    plain = canvas.render_export(1.0, Background.WHITE, False)
    window.stage.floating_panels()["vibrations"].arrows.click()
    with_arrows = canvas.render_export(1.0, Background.WHITE, False)
    assert plain.size() == with_arrows.size()
    assert plain != with_arrows
