"""Writing grids as cube files, and saving the grid on screen from the grid panel."""

from pathlib import Path

import numpy as np
import pytest
from openbabel import openbabel as ob
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QFileDialog

from miles.api import document_view
from miles.builtin_plugins import grids
from miles.builtin_plugins.grids import cube_file_name
from miles.core.chem import from_smiles
from miles.core.volume import Volume
from miles.files.cube import cube_text, write_cube
from miles.files.formats import read_file
from miles.ui.mainwindow import MainWindow

AXES = {
    "orthogonal": np.eye(3) * 0.25,
    "skewed": np.array([[0.2, 0.0, 0.0], [0.05, 0.25, 0.0], [0.0, 0.03, 0.3]]),
}

def grid_panel(window):
    """The grid panel of the Orbitals and Grids plugin."""
    return window.stage.floating_panels()["grid"]


def isosurface(canvas):
    """The isosurface mesh drawn on a canvas, or None."""
    items = canvas.drawing.scene_items()
    return items.meshes[0].mesh if items.meshes else None



@pytest.fixture
def water():
    mol, _ = from_smiles("O")
    return mol


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


def random_grid(axes, name="MO 5 (HOMO)"):
    values = np.random.default_rng(3).normal(size=(7, 9, 11)).astype(np.float32) * 0.1
    return Volume(name, np.array([-2.0, -2.5, -1.5]), axes, values)


@pytest.mark.parametrize("kind", sorted(AXES))
def test_cube_round_trip(tmp_path, water, kind):
    volume = random_grid(AXES[kind])
    path = tmp_path / "grid.cube"
    write_cube(str(path), water, volume, "water")
    mol, info = read_file(str(path))
    assert mol.symbols == water.symbols
    assert np.allclose(mol.coords, water.coords, atol=1e-5)
    (grid,) = info.volumes
    assert grid.name == volume.name
    assert grid.values.shape == volume.values.shape
    assert np.allclose(grid.values, volume.values, rtol=1e-5, atol=1e-6)
    assert np.allclose(grid.origin, volume.origin, atol=1e-5)
    assert np.allclose(grid.axes, volume.axes, atol=1e-5)


def test_open_babel_reads_the_file(tmp_path, water):
    path = tmp_path / "grid.cube"
    write_cube(str(path), water, random_grid(AXES["orthogonal"]))
    conv = ob.OBConversion()
    conv.SetInFormat("cube")
    obmol = ob.OBMol()
    assert conv.ReadFile(obmol, str(path))
    assert obmol.NumAtoms() == 3
    assert any(d.GetDataType() == ob.GridData for d in obmol.GetData())


def test_layout_follows_cubegen(water):
    volume = random_grid(AXES["orthogonal"], name="density\nwith a line break")
    lines = cube_text(water, volume, "").splitlines()
    assert lines[1] == "density with a line break"   # one comment line, whatever the name holds
    assert lines[2].split()[0] == "3"
    assert [int(line.split()[0]) for line in lines[3:6]] == [7, 9, 11]
    assert lines[6].split()[:2] == ["8", "8.000000"]
    data = lines[9:]
    assert len(data) == 7 * 9 * 2   # 11 values per row: a line of six and a line of five
    assert all(len(line.split()) in (5, 6) for line in data)
    assert cube_text(water, volume, "").isascii()


def test_cube_text_needs_a_three_dimensional_grid(water):
    flat = Volume("flat", np.zeros(3), np.eye(3), np.zeros((4, 4), dtype=np.float32))
    with pytest.raises(ValueError):
        cube_text(water, flat)


def test_suggested_file_name():
    assert cube_file_name(random_grid(np.eye(3), "MO 12 (alpha HOMO)")) == "MO_12_alpha_HOMO.cube"
    assert cube_file_name(random_grid(np.eye(3), "")) == "grid.cube"


def test_save_grid_from_the_panel(window, tmp_path, water, monkeypatch):
    source = tmp_path / "density.cube"
    write_cube(str(source), water, random_grid(AXES["orthogonal"], "Electron density"))
    window.open_path(str(source))
    panel = grid_panel(window)
    assert panel.isVisible() and panel.save_grid.isVisible()
    target = tmp_path / "saved.cube"
    offered = []

    def save_name(parent, caption, start, filters):
        offered.append(Path(start).name)
        return str(target), filters

    monkeypatch.setattr(QFileDialog, "getSaveFileName", save_name)
    panel.save_grid.click()
    assert offered == ["Electron_density.cube"]
    _, info = read_file(str(target))
    assert info.volumes[0].name == "Electron density"
    assert np.allclose(info.volumes[0].values, panel.selected_volume().values, atol=1e-6)


def test_computed_orbital_can_be_saved(window, tmp_path):
    import time

    from PySide6.QtWidgets import QApplication

    window.open_path(str(Path(__file__).parent / "data" / "h2_sto3g.molden"))
    canvas = window.current_canvas()
    end = time.time() + 10
    doc = document_view(window, canvas)
    while grid_panel(window).request(doc) is None:   # the HOMO is computed on opening
        assert time.time() < end
        QApplication.processEvents()
        time.sleep(0.01)
    panel = grid_panel(window)
    assert panel.save_grid.isVisible()
    path = tmp_path / "homo.cube"
    assert grids.save_grid(window.api, grid_panel(window), str(path))
    _, info = read_file(str(path))
    shown = grid_panel(window).volume_shown(document_view(window, canvas))[0]
    assert info.volumes[0].name == shown.name
    assert np.allclose(info.volumes[0].values, shown.values, rtol=1e-5, atol=1e-9)
