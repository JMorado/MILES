"""Quantities along trajectories and their plot in the timeline."""

import numpy as np
import pytest
from PySide6.QtCore import QPoint, QSettings, Qt
from PySide6.QtTest import QTest

from miles.builtin_plugins.superpose.matching import superpose
from miles.core.geometry import angle, dihedral, distance, rotation_matrix
from miles.core.trajectory import measurement_series, rmsd_atoms, rmsd_series, unwrap_degrees
from miles.files.common import JobInfo
from miles.ui.mainwindow import MainWindow

RNG = np.random.default_rng(3)


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


@pytest.mark.parametrize(("atoms", "measure"),
                         [((0, 3), distance), ((1, 2, 5), angle), ((0, 1, 2, 3), dihedral)])
def test_measurement_series_matches_single_measurements(atoms, measure):
    frames = RNG.normal(size=(25, 7, 3))
    expected = [measure(*frame[list(atoms)]) for frame in frames]
    assert np.allclose(measurement_series(frames, atoms), expected, atol=1e-10)


def test_measurement_needs_two_to_four_atoms():
    with pytest.raises(ValueError):
        measurement_series(np.zeros((2, 5, 3)), [0])


def test_unwrapped_dihedrals_have_no_jumps():
    turning = np.arange(0.0, 720.0, 15.0)                     # two full turns
    wrapped = (turning + 180.0) % 360.0 - 180.0
    unwrapped = unwrap_degrees(wrapped)
    assert np.abs(np.diff(unwrapped)).max() == pytest.approx(15.0)
    assert -180.0 <= unwrapped.mean() <= 180.0
    near_180 = unwrap_degrees([170.0, 179.0, -175.0, -160.0])  # crossing 180 is not a jump
    assert np.abs(np.diff(near_180)).max() < 20.0
    assert unwrap_degrees([]).size == 0


def test_rmsd_ignores_rotation_and_translation():
    first = RNG.normal(size=(10, 3))
    moved = first @ rotation_matrix([1.0, 2.0, 0.5], 1.1).T + [3.0, -2.0, 1.0]
    bent = first.copy()
    bent[0] += [0.6, 0.0, 0.0]
    rmsd = rmsd_series(np.array([first, moved, bent]))
    assert rmsd[0] == pytest.approx(0.0, abs=1e-7)
    assert rmsd[1] == pytest.approx(0.0, abs=1e-7)
    assert rmsd[2] > 0.1


def test_rmsd_matches_pairwise_superposition():
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CC(=O)NC")
    frames = np.array([mol.coords + RNG.normal(scale=0.2, size=mol.coords.shape) for _ in range(6)])
    atoms = rmsd_atoms(mol.symbols)
    assert atoms == [i for i, a in enumerate(mol.atoms) if a.symbol != "H"]
    values = rmsd_series(frames, atoms, reference=2)
    for k, frame in enumerate(frames):
        mobile, target = mol.copy(), mol.copy()
        mobile.coords, target.coords = frame.copy(), frames[2].copy()
        expected = superpose(mobile, target, [(i, i) for i in atoms]).rmsd
        assert values[k] == pytest.approx(expected, abs=1e-6)


def test_rmsd_uses_all_atoms_without_enough_heavy_atoms():
    assert rmsd_atoms(["O", "H", "H"]) == [0, 1, 2]


def test_timeline_plots_measurements_and_rmsd(window):
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CCCC")
    window.open_structure(mol, JobInfo(title="butane"))
    canvas = window.current_canvas()
    doc = canvas.doc
    frames = np.array([mol.coords + RNG.normal(scale=0.05, size=mol.coords.shape) for _ in range(8)])
    doc.set_frames(frames)
    timeline = window.stage.timeline
    assert timeline.isVisible() and not timeline.plot.isVisible()   # nothing to plot yet

    canvas.commands.toggle_measurement([0, 3])
    assert timeline.plot.isVisible()
    assert timeline.quantity is not None and timeline.quantity.key == "0,3"
    expected = distance(doc.frames[doc.frame][0], doc.frames[doc.frame][3])
    assert f"C1-C4 {expected:.2f} \u00c5" in timeline.text.text()
    assert np.allclose(timeline.plot.values, measurement_series(doc.frames, [0, 3]))

    timeline.choice.setCurrentIndex([q.key for q in timeline._quantities].index("rmsd"))
    assert timeline.quantity.key == "rmsd"
    assert "RMSD" in timeline.text.text()
    QTest.mouseClick(timeline.plot, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(timeline.plot.width() - 2, timeline.plot.height() // 2))
    assert doc.frame == 7
    timeline.plot.grab()   # paints without errors

    timeline.plot_button.click()   # hide, then show again with the same choice
    assert not timeline.plot.isVisible()
    timeline.plot_button.click()
    assert timeline.plot.isVisible() and timeline.quantity.key == "rmsd"

    canvas.commands.clear_measurements()
    assert [q.key for q in timeline._quantities] == ["rmsd"]
    assert not timeline.choice.isVisible()   # a single quantity needs no chooser


def test_energy_comes_first_and_dihedral_is_unwrapped(window):
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CCCC")
    window.open_structure(mol, JobInfo(title="butane"))
    canvas = window.current_canvas()
    doc = canvas.doc
    frames = np.array([mol.coords.copy() for _ in range(5)])
    doc.set_frames(frames, energies=[2.0, 1.0, 0.0, 3.0, 1.5])
    timeline = window.stage.timeline
    canvas.commands.toggle_measurement([0, 1, 2, 3])
    assert [q.key for q in timeline._quantities] == ["energy", "0,1,2,3", "rmsd"]
    assert timeline.quantity.key == "energy" and timeline.choice.isVisible()
    timeline.choice.setCurrentIndex(1)
    assert "\u00b0" in timeline.text.text() and "\u0394E" in timeline.text.text()
    timeline.plot.grab()   # a flat line paints without errors
