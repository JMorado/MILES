"""The scripting interface: documents, selection, edits, view and commands."""

from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QSettings

from miles.api import App, current_app
from miles.core.building import add_hydrogens
from miles.core.chem import from_smiles
from miles.ui.mainwindow import MainWindow

DATA = Path(__file__).parent / "data"


@pytest.fixture
def app_api(app, tmp_path):
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    win.resize(900, 650)
    win.show()
    app.processEvents()
    yield win.api
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


def test_the_window_offers_the_api(app_api):
    assert isinstance(app_api, App) and current_app() is app_api
    assert app_api.current is not None and len(app_api.documents) == 1


def test_new_documents_selection_and_undoable_edits(app_api):
    mol, _ = from_smiles("CCO")
    doc = app_api.new(mol, "ethanol")
    assert doc == app_api.current and doc.name == "ethanol"
    oxygen = doc.select(smarts="[OX2H]")
    assert [doc.molecule.atoms[i].symbol for i in oxygen] == ["O"]
    assert doc.select([0], add=True) == sorted([*oxygen, 0])
    with pytest.raises(ValueError):
        doc.select([999])

    before = len(doc.molecule)
    doc.edit(lambda m: m.remove_atoms([i for i, a in enumerate(m.atoms) if a.symbol == "H"]))
    assert len(doc.molecule) < before
    doc.undo()
    assert len(doc.molecule) == before
    with doc.editing() as m:
        m.coords = m.coords + np.array([0.0, 0.0, 1.0])
    doc.undo()
    assert doc.edit(add_hydrogens) == 0   # already complete


def test_view_style_turn_and_export(app_api, tmp_path):
    mol, _ = from_smiles("c1ccccc1")
    doc = app_api.new(mol, "benzene")
    doc.style(mode="tube", labels="symbol", hydrogen_bonds=True)   # a display option and a plugin's toggle
    assert str(doc.display.mode) == "tube" and app_api.window.plugins.toggles["hydrogen_bonds"].value(doc)
    with pytest.raises(AttributeError):
        doc.style(colour="red")
    before = doc.canvas.camera.rotation.copy()
    doc.turn(y=90)
    assert not np.allclose(doc.canvas.camera.rotation, before)
    doc.fit()
    picture = tmp_path / "benzene.png"
    doc.export_image(str(picture), scale=1.0)
    assert picture.stat().st_size > 1000
    doc.save(str(tmp_path / "benzene.miles"))
    assert doc.path.endswith("benzene.miles")


def test_files_frames_and_commands(app_api):
    doc = app_api.open(str(DATA / "water_freq.log"))
    assert doc.name == "water_freq.log" and doc.info.vibrations
    with pytest.raises(FileNotFoundError):
        app_api.open("no such file.xyz")
    commands = app_api.commands()
    assert "Edit > Clean Geometry" in commands or any(c.endswith("> Clean Geometry") for c in commands)
    mol, _ = from_smiles("CC")
    ethane = app_api.new(mol, "ethane")
    ethane.edit(lambda m: m.set_positions([1], m.coords[1] + [0.5, 0, 0]))
    app_api.run("clean geometry...")   # case and dots do not matter
    with pytest.raises(KeyError):
        app_api.run("No Such Command")
