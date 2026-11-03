"""Scene layout reuse, highlight flags, label limits and the per-layer line cache of the canvas."""

import numpy as np
import pytest
from PySide6.QtCore import QCoreApplication, QEvent, QPointF, QSettings, QTimer
from PySide6.QtWidgets import QApplication, QMenu

from miles.core.chem import from_smiles
from miles.core.molecule import Molecule
from miles.files.common import JobInfo
from miles.render.layers import Annotations, Line
from miles.render.picking import AtomHit
from miles.render.scene import CylinderKind, Highlight, Style, build_scene, scene_layout
from miles.render.theme import LIGHT
from miles.ui.canvas import overlay
from miles.ui.canvas.context_menu import show_context_menu
from miles.ui.mainwindow import MainWindow


@pytest.fixture
def window(app, tmp_path):
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat), plugins=False)
    win.resize(1000, 700)
    win.show()
    app.processEvents()
    yield win
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


def test_moved_atoms_with_a_kept_layout_give_the_same_cylinders_as_a_new_scene():
    mol, _ = from_smiles("C=CC=Cc1ccc2ccccc2c1C#N")
    style = Style()
    layout = scene_layout(mol, style, LIGHT)
    rng = np.random.default_rng(1)
    mol.coords = mol.coords + rng.normal(scale=0.2, size=mol.coords.shape)
    kept = build_scene(mol, style, LIGHT, layout=layout)
    fresh = build_scene(mol, style, LIGHT)
    assert kept.cylinders.shape == fresh.cylinders.shape
    assert np.allclose(kept.cylinders, fresh.cylinders, atol=1e-6)
    assert np.array_equal(kept.cylinder_trim, fresh.cylinder_trim)


def test_the_dashed_line_of_an_aromatic_bond_lies_inside_the_ring():
    mol, _ = from_smiles("c1ccccc1")
    scene = build_scene(mol, Style(), LIGHT)
    ring = mol.coords[mol.numbers == 6].mean(axis=0)
    rows = scene.cylinders
    dashed = rows[rows[:, 13] == CylinderKind.DASHED]
    solid = rows[(rows[:, 13] == CylinderKind.SOLID) & (rows[:, 6] < 0.07)]   # aromatic bonds' solid lines
    assert len(dashed) == 6 and len(solid) == 6
    middle = (dashed[:, 0:3] + dashed[:, 3:6]) / 2
    other = (solid[:, 0:3] + solid[:, 3:6]) / 2
    assert (np.linalg.norm(middle - ring, axis=1).max() < np.linalg.norm(other - ring, axis=1).min())


def test_a_neighbour_on_the_bond_axis_does_not_orient_a_double_bond():
    # H-C=C=C: the middle carbon's other neighbour is on the axis of the C=C bond
    mol = Molecule()
    for symbol, x in (("C", 0.0), ("C", 1.31), ("C", 2.62)):
        mol.add_atom(symbol, [x, 0.0, 0.0])
    mol.add_atom("H", [-0.6, 0.9, 0.0])
    mol.add_bond(0, 1, 2.0)
    mol.add_bond(1, 2, 2.0)
    mol.add_bond(0, 3, 1.0)
    rows = build_scene(mol, Style(), LIGHT).cylinders
    pair = rows[np.isclose(rows[:, 0], 0.0, atol=0.2) & np.isclose(rows[:, 3], 1.31, atol=0.2)]
    assert len(pair) == 2
    offset = pair[0, 0:3] - pair[1, 0:3]
    assert np.all(np.isfinite(rows))
    assert abs(offset[0]) < 1e-6 and np.linalg.norm(offset) > 0.1   # spread across the axis, not along it


def test_a_selected_atom_stays_selected_under_the_cursor(window):
    mol, _ = from_smiles("CCO")
    window.open_structure(mol, JobInfo(title="ethanol"))
    canvas = window.current_canvas()
    canvas.set_selection([0])
    canvas.hover = AtomHit(0)
    frame = canvas.drawing.frame()
    rows = {int(a): r for r, a in enumerate(canvas.scenes.base.sphere_atoms)}
    assert frame.flags[rows[0]] == Highlight.SELECTED
    canvas.hover = AtomHit(1)
    frame = canvas.drawing.frame()
    assert frame.flags[rows[0]] == Highlight.SELECTED
    assert frame.flags[rows[1]] not in (Highlight.NONE, Highlight.SELECTED)
    assert np.count_nonzero(frame.flags) == 2
    assert frame.scene(canvas.scenes.base).spheres[rows[1], 7] == frame.flags[rows[1]]


def test_labels_only_inside_the_window_and_at_most_the_nearest(window, monkeypatch):
    mol, _ = from_smiles("CCCCCCCCCC")
    window.open_structure(mol, JobInfo(title="decane"))
    canvas = window.current_canvas()
    n = len(canvas.mol)
    sx = np.linspace(-200, canvas.width() + 200, n)
    sy = np.full(n, canvas.height() / 2)
    depth = np.arange(n, dtype=float)
    inside = overlay.labelled_atoms(canvas, sx, sy, depth)
    assert len(inside) and np.all((sx[inside] > -overlay.LABEL_MARGIN)
                                  & (sx[inside] < canvas.width() + overlay.LABEL_MARGIN))
    assert len(inside) < n
    monkeypatch.setattr(overlay, "LABEL_LIMIT", 3)
    nearest = overlay.labelled_atoms(canvas, np.full(n, 10.0), sy, depth)
    assert nearest.tolist() == sorted(np.argsort(-depth)[:3].tolist())


def test_each_annotation_layer_keeps_its_resolved_lines(window):
    mol, _ = from_smiles("CCO")
    window.open_structure(mol, JobInfo(title="ethanol"))
    canvas = window.current_canvas()
    first = Annotations(lines=[Line(0, 1)])
    second = Annotations(lines=[Line(1, 2)])
    a = canvas.drawing.resolve_lines(first)
    b = canvas.drawing.resolve_lines(second)
    assert canvas.drawing.resolve_lines(first) is a   # still cached after the other layer was resolved
    assert canvas.drawing.resolve_lines(second) is b


def test_a_context_menu_is_deleted_after_it_closes(window, app):
    mol, _ = from_smiles("CCO")
    window.open_structure(mol, JobInfo(title="ethanol"))
    canvas = window.current_canvas()

    def close() -> None:
        popup = QApplication.activePopupWidget()
        if popup is not None:
            popup.close()

    for _ in range(3):
        QTimer.singleShot(20, close)
        show_context_menu(canvas, QPointF(40, 40), None)
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    app.processEvents()
    assert canvas.findChildren(QMenu) == []


def test_a_pasted_copy_is_moved_aside_when_it_would_overlap(window):
    mol, _ = from_smiles("c1ccccc1")
    window.open_structure(mol, JobInfo(title="benzene"))
    canvas = window.current_canvas()
    before = canvas.mol.coords.copy()
    canvas.commands.paste_molecule(canvas.mol.copy())
    pasted = canvas.mol.coords[len(before):]
    gaps = np.linalg.norm(pasted[:, None, :] - before[None, :, :], axis=-1)
    assert gaps.min() > 1.2
