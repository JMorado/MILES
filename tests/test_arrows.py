"""Arrows in 3D: their mesh, and the dipole and vibration arrows of the canvas."""

import numpy as np
import pytest
from PySide6.QtCore import QSettings

from miles.core.charges import dipole_moment
from miles.core.chem import from_smiles
from miles.files.common import JobInfo
from miles.render.arrows import COLUMNS, HEAD_LENGTH, arrow_mesh, arrow_parts, arrow_rows
from miles.render.theme import Background
from miles.ui.mainwindow import MainWindow


def edge_counts(mesh):
    """How many triangles share each edge, with vertices at the same place merged."""
    ids = {}
    merged = np.array([ids.setdefault(tuple(np.round(v, 5)), len(ids)) for v in mesh.vertices])
    counts = {}
    for a, b, c in merged[mesh.indices.reshape(-1, 3)]:
        for edge in ((a, b), (b, c), (c, a)):
            key = (min(edge), max(edge))
            counts[key] = counts.get(key, 0) + 1
    return counts


def test_arrow_mesh_is_closed_and_points_to_the_tip():
    rows = arrow_rows([[0, 0, 0], [1, 1, 1]], [[0, 0, 2], [1, 3, 1]], 0.08, (0.2, 0.4, 0.9))
    mesh = arrow_mesh(rows)
    assert mesh is not None
    assert set(edge_counts(mesh).values()) == {2}   # watertight, so the outline pass works
    assert np.allclose(np.linalg.norm(mesh.normals, axis=1), 1.0, atol=1e-5)
    triangles = mesh.indices.reshape(-1, 3)
    a, b, c = (mesh.vertices[triangles[:, k]] for k in range(3))
    facing = np.einsum("ij,ij->i", np.cross(b - a, c - a), mesh.normals[triangles].mean(axis=1))
    assert (facing > 0).all()   # wound outwards
    for tip in ([0, 0, 2], [1, 3, 1]):
        assert np.linalg.norm(mesh.vertices - tip, axis=1).min() < 1e-6
    assert mesh.vertices[:, 2].max() == pytest.approx(2.0, abs=1e-6)   # nothing beyond the tip
    assert np.allclose(mesh.colors[0], (0.2, 0.4, 0.9))


def test_short_arrows_get_half_their_length_as_head_and_empty_ones_nothing():
    rows = arrow_rows([[0, 0, 0], [5, 5, 5]], [[0, 0, 0.2], [5, 5, 5]], 0.1, (1, 0, 0))
    _, bases, _, head_radius = arrow_parts(rows)
    assert np.allclose(bases[0], [0, 0, 0.1])           # head limited to half the arrow
    assert head_radius[0] > 0.1
    long = arrow_parts(arrow_rows([[0, 0, 0]], [[0, 0, 5]], 0.1, (1, 0, 0)))
    assert np.allclose(long[1][0], [0, 0, 5 - HEAD_LENGTH * 0.1])
    assert arrow_mesh(rows[1:]) is None                 # zero length
    assert arrow_mesh(np.zeros((0, COLUMNS))) is None


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


def test_dipole_arrow_spans_the_molecule(window):
    mol, _ = from_smiles("Nc1ccc(cc1)C(=O)O")
    window.open_structure(mol, JobInfo(title="4-aminobenzoic acid"))
    canvas = window.current_canvas()
    assert len(canvas.drawing.scene_items().arrows) == 0
    window.api.current.style(dipole=True)
    arrows = canvas.drawing.scene_items().arrows
    assert len(arrows) == 1
    tail, tip = arrows[0, 0:3].astype(float), arrows[0, 3:6].astype(float)
    centre = canvas.mol.center_of_mass()
    reach = np.linalg.norm(canvas.mol.coords - centre, axis=1).max()
    assert np.linalg.norm(tip - centre) > reach and np.linalg.norm(tail - centre) > reach   # ends outside
    direction = (tip - tail) / np.linalg.norm(tip - tail)
    moment = dipole_moment(canvas.mol, canvas.doc.charges())
    assert np.linalg.norm(moment) > 0.01
    assert np.allclose(direction, moment / np.linalg.norm(moment), atol=1e-5)
    frame = canvas.drawing.frame()
    assert frame.arrow_mesh is not None and not frame.arrow_mesh.empty
    assert not canvas.render_export(1.0, Background.WHITE, False).isNull()


def test_mode_arrows_start_at_the_moving_atoms(window):
    from pathlib import Path

    window.open_path(str(Path(__file__).parent / "data" / "water_freq.log"))
    canvas = window.current_canvas()
    panel = window.stage.floating_panels()["vibrations"]
    panel.animate(0)
    panel.stop_animation()
    panel.arrows.click()
    arrows = canvas.drawing.scene_items().arrows
    vectors = panel.arrows_for(window.api.current)
    moving = np.linalg.norm(vectors, axis=1) > 0
    assert len(arrows) == moving.sum()
    assert np.allclose(arrows[:, 0:3], canvas.mol.coords[moving])
    assert np.allclose(arrows[:, 3:6] - arrows[:, 0:3], vectors[moving], atol=1e-5)
    canvas.grab()   # the QPainter renderer draws shafts and cones
    panel.arrows.click()
    assert len(canvas.drawing.scene_items().arrows) == 0
