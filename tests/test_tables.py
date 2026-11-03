"""Tables of plot data, and copying and exporting them from the panels."""

import csv
import io
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QSettings
from PySide6.QtGui import QImage
from PySide6.QtWidgets import QApplication, QFileDialog

from miles.builtin_plugins.ramachandran.angles import backbone_angles, backbone_residues, ramachandran_table
from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
from miles.builtin_plugins.vibrations.tables import mode_table, spectrum_table
from miles.core.tables import Table, frame_table
from miles.core.trajectory import measurement_series
from miles.core.vibrations import VibrationalMode
from miles.files.common import JobInfo
from miles.ui.mainwindow import MainWindow
from miles.ui.plot_export import COPY, CURVE, DATA, IMAGE

DATA_DIR = Path(__file__).parent / "data"
EXAMPLES = Path(__file__).parent.parent / "examples"


def rows_of(text: str) -> list[list[str]]:
    return list(csv.reader(io.StringIO(text)))


def test_table_text_formats():
    table = Table(["Name", "Value (unit)"], [["a, b", "1.5"], ["c", ""]])
    assert table.tsv() == "Name\tValue (unit)\na, b\t1.5\nc\t\n"
    assert rows_of(table.csv()) == [["Name", "Value (unit)"], ["a, b", "1.5"], ["c", ""]]   # commas quoted


def test_mode_and_spectrum_tables():
    still = np.zeros((2, 3))
    modes = [VibrationalMode(-120.0, None, still), VibrationalMode(1700.0, 250.0, still),
             VibrationalMode(3000.0, 50.0, still)]
    table = mode_table(modes)
    assert table.header == ["Mode", "Wavenumber (cm-1)", "IR intensity (km/mol)"]
    assert table.rows[0] == ["1", "-120.00", ""]    # imaginary, no intensity
    assert table.rows[1] == ["2", "1700.00", "250.0000"]
    curve = spectrum_table(modes, 400.0, 4000.0, 1.0, 20.0)
    assert curve.header == ["Wavenumber (cm-1)", "Absorbance (relative)"]
    assert len(curve.rows) == 3601 and curve.rows[0][0] == "4000.0" and curve.rows[-1][0] == "400.0"
    peak = max(curve.rows, key=lambda row: float(row[1]))
    assert peak == ["1700.0", "1.000000"]


def test_frame_table_has_every_quantity_for_every_frame():
    rng = np.random.default_rng(1)
    frames = rng.normal(size=(5, 4, 3)) + np.arange(4)[None, :, None] * 1.5
    table = frame_table(frames, ["C", "C", "O", "H"], energies=[3.0, 1.0, 2.0, 5.0, 1.5],
                        labels=["a", "b", "", "d", "e"],
                        measurements=[("Distance C1-C2", [0, 1]), ("Dihedral C1-C2-O3-H4", [0, 1, 2, 3])])
    assert table.header == ["Frame", "Label", "Energy (kcal/mol)", "Relative energy (kcal/mol)",
                            "Distance C1-C2 (angstrom)", "Dihedral C1-C2-O3-H4 (deg)",
                            "RMSD from frame 1 (angstrom)"]
    assert len(table.rows) == 5
    assert table.rows[1][:4] == ["2", "b", "1.0000", "0.0000"]
    distances = measurement_series(frames, [0, 1])
    assert [float(row[4]) for row in table.rows] == pytest.approx(distances, abs=1e-4)
    assert float(table.rows[0][6]) == pytest.approx(0.0, abs=1e-8)
    plain = frame_table(frames, ["C", "C", "O", "H"])
    assert plain.header == ["Frame", "RMSD from frame 1 (angstrom)"]


def test_ramachandran_table():
    mol = build_sequence("AGPAA", SequenceKind.HELIX)
    found = backbone_residues(mol)
    table = ramachandran_table(found, backbone_angles(found, mol.coords))
    assert table.header == ["Chain", "Residue", "Number", "Phi (deg)", "Psi (deg)", "Region"]
    assert [row[1] for row in table.rows] == ["ALA", "GLY", "PRO", "ALA", "ALA"]
    assert table.rows[0][3] == "" and table.rows[0][5] == ""   # no phi at the N-terminus
    assert float(table.rows[2][3]) == pytest.approx(-57, abs=15)   # helix
    assert table.rows[2][5] in ("favoured", "allowed", "outlier")


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


def save_to(monkeypatch, path):
    monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *args, **kwargs: (str(path), ""))


def test_vibrations_data_and_picture(window, tmp_path, monkeypatch):
    window.open_path(str(DATA_DIR / "water_freq.log"))
    panel = window.stage.floating_panels()["vibrations"]
    assert [text for text, _ in panel.menu_entries()] == [
        "Copy Data", "Export Data...", "Export Spectrum Curve...", "Export Image..."]
    panel.exportRequested.emit(COPY)
    copied = QApplication.clipboard().text().splitlines()
    assert copied[0] == "Mode\tWavenumber (cm-1)\tIR intensity (km/mol)" and len(copied) == 4
    save_to(monkeypatch, tmp_path / "modes.csv")
    panel.exportRequested.emit(DATA)
    assert len(rows_of((tmp_path / "modes.csv").read_text())) == 4
    save_to(monkeypatch, tmp_path / "curve.csv")
    panel.exportRequested.emit(CURVE)
    assert rows_of((tmp_path / "curve.csv").read_text())[0] == ["Wavenumber (cm-1)", "Absorbance (relative)"]
    for name in ("spectrum.png", "spectrum.svg"):
        save_to(monkeypatch, tmp_path / name)
        panel.exportRequested.emit(IMAGE)
        assert (tmp_path / name).stat().st_size > 1000
    image = QImage(str(tmp_path / "spectrum.png"))
    assert image.width() == 3 * panel.spectrum.width()
    assert (tmp_path / "spectrum.svg").read_text().lstrip().startswith("<?xml")


def test_timeline_data(window, tmp_path, monkeypatch):
    window.open_path(str(EXAMPLES / "alanine_dipeptide_md.xyz"))
    canvas = window.current_canvas()
    canvas.commands.toggle_measurement([1, 3, 4, 6])
    timeline = window.stage.timeline
    assert [text for text, _ in timeline.menu_entries()] == ["Copy Data", "Export Data..."]
    save_to(monkeypatch, tmp_path / "frames.csv")
    timeline.exportRequested.emit(DATA)
    rows = rows_of((tmp_path / "frames.csv").read_text())
    assert len(rows) == 201 and rows[0][0] == "Frame" and rows[0][-1] == "RMSD from frame 1 (angstrom)"
    assert any(title.startswith("Dihedral") and title.endswith("(deg)") for title in rows[0])


def test_ramachandran_data_and_picture(window, tmp_path, monkeypatch):
    window.open_structure(build_sequence("AEAAAKEAAAKA", SequenceKind.HELIX), JobInfo(title="helix"))
    window.api.run("Ramachandran Plot")
    panel = window.stage.floating_panels()["ramachandran"]
    assert panel.isVisible()
    panel.exportRequested.emit(COPY)
    copied = QApplication.clipboard().text().splitlines()
    assert copied[0].startswith("Chain\tResidue\tNumber") and len(copied) == 13
    save_to(monkeypatch, tmp_path / "rama.png")
    panel.exportRequested.emit(IMAGE)
    assert not QImage(str(tmp_path / "rama.png")).isNull()
