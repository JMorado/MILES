"""The grid panel: orbitals from Molden files and surfaces coloured by a second grid."""

import time
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QPointF, QSettings, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication, QMessageBox

from miles.api import document_view
from miles.builtin_plugins import grids
from miles.ui.mainwindow import MainWindow

DATA = Path(__file__).parent / "data"
BOHR = 0.529177210903
WATER = [(8, (0.0, 0.0, 0.117)), (1, (0.0, 0.757, -0.469)), (1, (0.0, -0.757, -0.469))]

def grid_panel(window):
    """The grid panel of the Orbitals and Grids plugin."""
    return window.stage.floating_panels()["grid"]


def isosurface(canvas):
    """The isosurface mesh drawn on a canvas, or None."""
    items = canvas.drawing.scene_items()
    return items.meshes[0].mesh if items.meshes else None


def shown_grid(window, canvas):
    """The grid a canvas shows as an isosurface."""
    return grid_panel(window).volume_shown(document_view(window, canvas))[0]



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


def wait_until(condition, timeout=10.0):
    """Process events until ``condition()`` holds; worker threads need the GIL, so sleep between rounds."""
    end = time.time() + timeout
    while not condition():
        assert time.time() < end, "timed out"
        QApplication.processEvents()
        time.sleep(0.01)


def write_cube(path, function, atoms=WATER, n=24):
    """Cube file (bohr) of ``function(points)`` on a 6 angstrom box around the origin."""
    origin = np.full(3, -3.0 / BOHR)
    step = 6.0 / BOHR / (n - 1)
    points = np.indices((n, n, n)).reshape(3, -1).T * step + origin
    values = function(points * BOHR)
    lines = ["test job", path.stem, f"{len(atoms):5d}" + "".join(f"{v:12.6f}" for v in origin)]
    lines += [f"{n:5d}" + "".join(f"{step if a == b else 0.0:12.6f}" for b in range(3)) for a in range(3)]
    lines += [f"{z:5d}{float(z):12.6f}" + "".join(f"{c / BOHR:12.6f}" for c in xyz) for z, xyz in atoms]
    rows = values.reshape(n * n, n)
    lines += ["".join(f"{v:13.5E}" for v in row[k:k + 6]) for row in rows for k in range(0, n, 6)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def density(points):
    return np.exp(-(points ** 2).sum(axis=1))


def potential(points):
    return 0.05 * points[:, 2]


def test_molden_orbitals_are_listed_and_computed_once(window):
    window.open_path(str(DATA / "h2_sto3g.molden"))
    canvas = window.current_canvas()
    panel = grid_panel(window)
    orbitals = canvas.doc.info.orbitals
    assert orbitals is not None and panel.isVisible()
    doc = document_view(window, canvas)
    wait_until(lambda: grid_panel(window).request(doc) is not None)   # the HOMO is shown on opening
    assert shown_grid(window, canvas).name == orbitals.name(orbitals.homo())
    assert panel.choice.count() == len(orbitals.orbitals) + 1   # and the electron density
    assert "HOMO" in panel.choice.currentText()

    lumo = next(k for k in range(len(orbitals.orbitals)) if orbitals.label(k) == "LUMO")
    panel.choice.setCurrentIndex(panel._entries.index(("orbital", lumo)))
    wait_until(lambda: shown_grid(window, canvas).name == orbitals.name(lumo))
    assert len(canvas.doc.info.volumes) == 2
    assert isosurface(canvas) is not None and not isosurface(canvas).empty

    panel.choice.setCurrentIndex(panel._entries.index(("orbital", orbitals.homo())))
    assert shown_grid(window, canvas).name == orbitals.name(orbitals.homo())   # cached
    assert len(canvas.doc.info.volumes) == 2


def test_added_grid_colours_the_density(window, tmp_path):
    write_cube(tmp_path / "density.cube", density)
    write_cube(tmp_path / "potential.cube", potential)
    window.open_path(str(tmp_path / "density.cube"))
    canvas = window.current_canvas()
    panel = grid_panel(window)
    assert canvas.doc.info.volumes[0].name == "density"
    assert grids.add_grid(window.api, grid_panel(window), str(tmp_path / "potential.cube"))
    request = grid_panel(window).request(document_view(window, canvas))
    assert request is not None and request.index == 0 and request.mapped == 1
    assert request.mapped_range > 0 and not request.see_through
    assert panel.mapping_row.isVisible() and panel.bar.isVisible()
    assert not panel.levels.isVisible()   # no orbitals in a cube file
    colors = isosurface(canvas).colors
    assert np.ptp(colors[:, 0]) > 0.3   # red on one side, blue on the other

    panel.mapped.setCurrentIndex(0)   # None
    assert grid_panel(window).request(document_view(window, canvas)).mapped is None
    panel.choice.setCurrentIndex(1)   # the potential itself: signed, so it cannot be coloured
    request = grid_panel(window).request(document_view(window, canvas))
    assert request.index == 1 and request.mapped is None
    assert not panel.mapped.isVisible() and not panel.add_grid.isVisible()
    assert panel.save_grid.isVisible()   # the row stays for Save Grid


def test_grid_of_another_structure_is_refused(window, tmp_path, monkeypatch):
    warnings = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *args: warnings.append(args[2]))
    write_cube(tmp_path / "density.cube", density)
    other = [(8, (0.0, 0.0, 0.0)), (1, (0.0, 0.9, -0.5)), (1, (0.0, -0.9, -0.5))]
    write_cube(tmp_path / "other.cube", potential, atoms=other)
    window.open_path(str(tmp_path / "density.cube"))
    assert not grids.add_grid(window.api, grid_panel(window), str(tmp_path / "other.cube"))
    assert len(window.current_canvas().doc.info.volumes) == 1
    assert "different structure" in warnings[0]


def test_orbital_diagram_chooses_by_clicking(window):
    window.open_path(str(DATA / "h2_sto3g.molden"))
    canvas = window.current_canvas()
    panel = grid_panel(window)
    orbitals = canvas.doc.info.orbitals
    wait_until(lambda: grid_panel(window).request(document_view(window, canvas)) is not None)
    levels = panel.levels
    assert levels.isVisible() and levels.current == orbitals.homo()
    lumo = next(k for k in range(len(orbitals.orbitals)) if orbitals.label(k) == "LUMO")
    homo_level, lumo_level = levels.frontier(0)
    assert (homo_level.index, lumo_level.index) == (orbitals.homo(), lumo)
    assert "LUMO" in levels.tooltip(lumo) and "eV" in levels.tooltip(lumo)

    point = levels.tick_position(lumo)
    assert point is not None and levels.level_at(point) == lumo
    assert levels.level_at(QPointF(point.x(), -50)) is None
    QTest.mouseClick(levels, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, point.toPoint())
    wait_until(lambda: shown_grid(window, canvas).name == orbitals.name(lumo))
    assert levels.current == lumo
    assert "LUMO" in panel.choice.currentText()


def test_orbital_diagram_zooms_and_resets(window):
    window.open_path(str(DATA / "h2_sto3g.molden"))
    levels = grid_panel(window).levels
    start = (levels.lo, levels.hi)
    levels._zoom(3, levels.width() / 2)
    assert levels.hi - levels.lo < start[1] - start[0]
    levels._pan(40)
    levels.reset_view()
    assert (levels.lo, levels.hi) == pytest.approx(start)


def test_orbital_diagram_starts_on_the_valence_levels(window):
    examples = Path(__file__).parent.parent / "examples" / "qm" / "formaldehyde.molden"
    window.open_path(str(examples))
    levels = grid_panel(window).levels
    homo, lumo = levels.frontier(0)
    assert min(lv.energy for lv in levels.levels) < -200   # oxygen and carbon 1s
    assert levels.lo > homo.energy - 25 and levels.hi < lumo.energy + 25


@pytest.mark.parametrize("appearance", ["light", "dark"])
def test_panel_menus_keep_the_theme_after_being_restyled(app, tmp_path, appearance):
    """A theme change while the panel was hidden used to be undone when Qt restyled its menus later."""
    from PySide6.QtGui import QPalette

    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    settings.setValue("appearance", "dark" if appearance == "light" else "light")
    win = MainWindow(settings)   # starts in the other theme, with the panel hidden
    win.show()
    app.processEvents()
    win.set_appearance(appearance)
    for combo in (grid_panel(win).choice, grid_panel(win).mapped):
        combo.style().unpolish(combo)   # what showing the panel does on macOS
        combo.style().polish(combo)
        assert combo.palette().color(QPalette.ColorRole.ButtonText).name() == win.theme.text.name()
    win.close()


def test_orbital_levels_that_would_overlap_are_stacked(window):
    examples = Path(__file__).parent.parent / "examples" / "qm" / "nitrogen_dioxide.molden"
    window.open_path(str(examples))
    levels = grid_panel(window).levels

    def overlapping():
        found = []
        for a in levels.levels:
            for b in levels.levels:
                if a.index < b.index and a.row == b.row:
                    ra, rb = levels._tick_rect(a), levels._tick_rect(b)
                    if ra.intersects(rb):
                        found.append((a.index, b.index))
        return found

    assert overlapping() == []
    for lv in levels.levels:
        levels._tick_rect(lv)   # lay out
    pairs = [(a, b) for a in levels.levels for b in levels.levels
             if a.row == b.row and 0.05 < b.energy - a.energy < 0.5 and a.slots > 1 and b.slots > 1
             and levels.lo <= a.energy <= levels.hi]
    assert pairs                              # close, non-degenerate levels share a place
    a, b = pairs[0]
    levels._zoom(30, levels._x((a.energy + b.energy) / 2))   # zoom in on them
    assert overlapping() == []
    levels._tick_rect(a)
    assert abs(levels._x(b.energy) - levels._x(a.energy)) > levels.TICK + 1
    assert levels._tick_rect(a).top() == levels._tick_rect(b).top() or a.slots == 1   # side by side now


def test_sessions_reopen_the_orbital_and_isovalue_shown(window, tmp_path):
    window.open_path(str(DATA / "h2_sto3g.molden"))
    canvas = window.current_canvas()
    panel = grid_panel(window)
    orbitals = canvas.doc.info.orbitals
    wait_until(lambda: grid_panel(window).request(document_view(window, canvas)) is not None)
    lumo = next(k for k in range(len(orbitals.orbitals)) if orbitals.label(k) == "LUMO")
    panel.choose_orbital(lumo)
    wait_until(lambda: shown_grid(window, canvas).name == orbitals.name(lumo))
    panel.level.set_value(0.037)
    panel._update()
    path = tmp_path / "h2.miles"
    canvas.doc.path = str(path)
    assert window.save_document(canvas.doc)

    window.open_path(str(path))
    reopened = window.current_canvas()
    assert reopened is not canvas and grid_panel(window).request(document_view(window, reopened)) is not None
    shown = shown_grid(window, reopened)
    assert shown.name == orbitals.name(lumo)
    assert grid_panel(window).request(document_view(window, reopened)).level == pytest.approx(0.037)
    for _ in range(20):   # no HOMO is computed on top of what the session shows
        QApplication.processEvents()
    assert shown_grid(window, reopened).name == orbitals.name(lumo)
    assert "LUMO" in panel.choice.currentText()


def test_sessions_keep_a_surface_coloured_by_another_grid(window, tmp_path):
    write_cube(tmp_path / "density.cube", density)
    write_cube(tmp_path / "potential.cube", potential)
    window.open_path(str(tmp_path / "density.cube"))
    assert grids.add_grid(window.api, grid_panel(window), str(tmp_path / "potential.cube"))
    canvas = window.current_canvas()
    before = grid_panel(window).request(document_view(window, canvas))
    path = tmp_path / "mapped.miles"
    canvas.doc.path = str(path)
    assert window.save_document(canvas.doc)
    window.open_path(str(path))
    reopened = window.current_canvas()
    assert grid_panel(window).request(document_view(window, reopened)) == before
    assert grid_panel(window).mapping_row.isVisible()


def test_densities_are_offered_with_the_orbitals(window):
    from miles.core.orbitals import DENSITY_NAME, SPIN_NAME

    examples = Path(__file__).parent.parent / "examples" / "qm"
    window.open_path(str(examples / "nitrogen_dioxide.molden"))
    canvas = window.current_canvas()
    panel = grid_panel(window)
    doc = document_view(window, canvas)
    wait_until(lambda: grid_panel(window).request(doc) is not None)   # the HOMO first
    assert ("density", 0) in panel._entries and ("spin", 0) in panel._entries
    panel.choice.setCurrentIndex(panel._entries.index(("spin", 0)))
    wait_until(lambda: shown_grid(window, canvas).name == SPIN_NAME)
    spin = shown_grid(window, canvas)
    assert spin.signed and panel.sign.isVisible()
    assert grid_panel(window).request(doc).level == pytest.approx(spin.suggested_level())
    panel.choice.setCurrentIndex(panel._entries.index(("density", 0)))
    wait_until(lambda: shown_grid(window, canvas).name == DENSITY_NAME)
    count = len(canvas.doc.info.volumes)
    panel.choice.setCurrentIndex(panel._entries.index(("spin", 0)))   # computed once, found again
    assert shown_grid(window, canvas).name == SPIN_NAME
    assert len(canvas.doc.info.volumes) == count

    window.open_path(str(DATA / "h2_sto3g.molden"))   # closed shell: no spin density
    assert ("spin", 0) not in grid_panel(window)._entries


def test_combining_two_grids(window, tmp_path, monkeypatch):
    from miles.builtin_plugins.grids.combine_dialog import CombineDialog

    write_cube(tmp_path / "density.cube", density)
    write_cube(tmp_path / "shifted.cube", lambda p: density(p + np.array([0.3, 0.0, 0.0])))
    window.open_path(str(tmp_path / "density.cube"))
    panel = grid_panel(window)
    assert not panel.combine_link.isVisible()   # one grid: nothing to combine
    assert grids.add_grid(window.api, grid_panel(window), str(tmp_path / "shifted.cube"))
    assert panel.combine_link.isVisible()
    canvas = window.current_canvas()

    zero = panel.combine(0, 0)
    assert np.abs(canvas.doc.info.volumes[zero].values).max() == 0.0

    panel.choice.setCurrentIndex(panel._entries.index(("volume", 0)))   # A is the grid on screen
    monkeypatch.setattr(CombineDialog, "exec", lambda self: 1)
    panel.combine_dialog()
    difference = canvas.doc.info.volumes[-1]
    assert difference.name == "density - shifted"
    assert grid_panel(window).request(document_view(window, canvas)).index == len(canvas.doc.info.volumes) - 1
    assert difference.signed and difference.values.max() > 0 > difference.values.min()
    assert not isosurface(canvas).empty
