"""Surfaces and isosurfaces drawn by the QPainter renderer."""

import time
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QSettings
from PySide6.QtGui import QColor, QImage, QPainter
from PySide6.QtWidgets import QApplication

from miles.api import document_view
from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
from miles.core.chem import from_smiles
from miles.core.surface import molecular_surface
from miles.files.common import JobInfo
from miles.render.camera import Camera
from miles.render.cartoon import Mesh
from miles.render.painter import paint_scene
from miles.render.scene import Style, build_scene
from miles.render.theme import LIGHT, Background
from miles.render.triangles import cluster_vertices, screen_triangles, simplified
from miles.ui.mainwindow import MainWindow

SIZE = 300

def grid_panel(window):
    """The grid panel of the Orbitals and Grids plugin."""
    return window.stage.floating_panels()["grid"]


def isosurface(canvas):
    """The isosurface mesh drawn on a canvas, or None."""
    items = canvas.drawing.scene_items()
    return items.meshes[0].mesh if items.meshes else None



def surface_of(mol, color=(0.2, 0.4, 0.9)):
    s = molecular_surface(mol, tuple(range(len(mol))))
    colors = np.tile(np.array(color, np.float32), (len(s.vertices), 1))
    return Mesh(vertices=s.vertices, normals=s.normals, colors=colors, base_colors=colors, indices=s.faces)


def picture(mol, meshes=()):
    camera = Camera()
    camera.center = mol.coords.mean(axis=0)
    camera.scale = 28.0
    image = QImage(SIZE, SIZE, QImage.Format.Format_ARGB32_Premultiplied)
    p = QPainter(image)
    paint_scene(p, SIZE, SIZE, build_scene(mol, Style(), LIGHT), camera, LIGHT, Style(), meshes=meshes)
    p.end()
    return image, camera


def pixel(image, x, y):
    c = QColor(image.pixel(int(x), int(y)))
    return np.array([c.redF(), c.greenF(), c.blueF()])


def test_large_meshes_are_simplified_once():
    mol = build_sequence("ACDEFGHIKLMNPQRSTVWY" * 2, SequenceKind.HELIX)
    mesh = surface_of(mol)
    small = simplified(mesh, 4000)
    assert small.indices.size // 3 <= 4000 < mesh.indices.size // 3
    assert simplified(mesh, 4000) is small                       # cached
    assert simplified(small, 4000) is small                      # small meshes pass through
    lo, hi = mesh.vertices.min(axis=0), mesh.vertices.max(axis=0)
    assert np.all(small.vertices >= lo - 1e-3) and np.all(small.vertices <= hi + 1e-3)
    assert np.allclose(np.linalg.norm(small.normals, axis=1), 1.0, atol=1e-4)


def test_clustering_merges_close_vertices_and_drops_collapsed_triangles():
    vertices = np.array([[0, 0, 0], [0.01, 0, 0], [1, 0, 0], [0, 1, 0]], np.float32)
    normals = np.tile(np.array([0, 0, 1], np.float32), (4, 1))
    colors = np.ones((4, 3), np.float32)
    mesh = Mesh(vertices=vertices, normals=normals, colors=colors, base_colors=colors,
                indices=np.array([0, 1, 2, 1, 2, 3, 0, 2, 3], np.uint32))
    merged = cluster_vertices(mesh, 0.5)
    assert len(merged.vertices) == 3
    assert merged.indices.size // 3 == 1   # one triangle collapsed, the two others became the same


def test_surface_hides_the_atoms_and_see_through_shows_them():
    mol, _ = from_smiles("O")
    plain, camera = picture(mol)
    x, y = camera.to_screen(camera.to_view(mol.coords[[0]]), SIZE, SIZE)
    oxygen = pixel(plain, x[0], y[0])
    mesh = surface_of(mol, (0.2, 0.4, 0.9))
    opaque, _ = picture(mol, [(mesh, 1.0)])
    see_through, _ = picture(mol, [(mesh, 0.45)])
    covered = pixel(opaque, x[0], y[0])
    assert covered[2] > covered[0]                     # blue surface in front of the red oxygen
    shown = pixel(see_through, x[0], y[0])
    assert shown[0] > covered[0] + 0.1                 # some of the oxygen's red comes through
    assert oxygen[0] > shown[0]                        # but tinted by the surface
    assert np.abs(pixel(opaque, 2, 2) - pixel(plain, 2, 2)).max() < 1e-6   # background untouched


def test_back_faces_are_not_drawn():
    mol, _ = from_smiles("C")
    camera = Camera()
    camera.center = mol.coords.mean(axis=0)
    camera.scale = 28.0
    triangles = screen_triangles(surface_of(mol), camera, SIZE, SIZE)
    total = surface_of(mol).indices.size // 3
    assert 0.3 * total < len(triangles.depth) < 0.7 * total


def test_painting_a_protein_surface_is_quick():
    mol = build_sequence("ACDEFGHIKLMNPQRSTVWY" * 5, SequenceKind.HELIX)
    mesh = surface_of(mol)
    picture(mol, [(mesh, 1.0)])   # simplifies and caches
    start = time.perf_counter()
    picture(mol, [(mesh, 1.0)])
    assert time.perf_counter() - start < 1.0


@pytest.fixture
def window(app, tmp_path):
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    win.resize(900, 650)
    win.show()
    app.processEvents()
    yield win
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


def test_software_canvas_draws_surfaces_and_isosurfaces(window):
    mol, _ = from_smiles("CC(=O)O")
    window.open_structure(mol, JobInfo(title="acetic acid"))
    canvas = window.current_canvas()
    assert type(canvas.surface).__name__ == "PainterSurface"
    before = canvas.grab().toImage()
    window.api.run("View > Surface > Show Surface")
    canvas.grab()        # asks for the surface, which is worked out in the background
    window.api.wait()
    after = canvas.grab().toImage()
    assert after != before
    export = canvas.render_export(1.0, Background.WHITE, False)
    assert not export.isNull()

    window.open_path(str(Path(__file__).parent / "data" / "h2_sto3g.molden"))
    canvas = window.current_canvas()
    for _ in range(200):   # the HOMO is computed in the background
        if grid_panel(window).request(document_view(window, canvas)) is not None:
            break
        time.sleep(0.01)
        QApplication.processEvents()
    assert isosurface(canvas) is not None
    with_orbital = canvas.grab().toImage()
    grid_panel(window).visible.click()
    assert canvas.grab().toImage() != with_orbital
