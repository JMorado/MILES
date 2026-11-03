"""Smoke tests of the interface, run offscreen with the QPainter canvas (see conftest.py)."""

import numpy as np
import pytest
from PySide6.QtCore import QPoint, QPointF, QSettings, Qt, QTimer
from PySide6.QtGui import QMouseEvent
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QApplication

from miles.api import document_view
from miles.core.molecule import Representation
from miles.core.protein import Group
from miles.files.common import JobInfo
from miles.files.formats import parse_text
from miles.render.picking import AtomHit
from miles.render.scene import LabelMode
from miles.render.theme import Background
from miles.ui.builder import MouseMode, Tool
from miles.ui.canvas import Canvas, PainterSurface
from miles.ui.export_dialog import crop_image
from miles.ui.mainwindow import MainWindow

PEPTIDE = """\
ATOM      1  N   ALA A   1      -0.677  -1.230  -0.491  1.00  0.00           N
ATOM      2  CA  ALA A   1      -0.001   0.064  -0.491  1.00  0.00           C
ATOM      3  C   ALA A   1       1.499  -0.110  -0.491  1.00  0.00           C
ATOM      4  O   ALA A   1       2.030  -1.227  -0.502  1.00  0.00           O
ATOM      5  CB  ALA A   1      -0.509   0.856   0.727  1.00  0.00           C
ATOM      6  N   ALA A   2       2.250   0.989  -0.479  1.00  0.00           N
ATOM      7  CA  ALA A   2       3.700   0.900  -0.480  1.00  0.00           C
ATOM      8  C   ALA A   2       4.300   2.300  -0.470  1.00  0.00           C
ATOM      9  O   ALA A   2       3.600   3.300  -0.460  1.00  0.00           O
ATOM     10  CB  ALA A   2       4.200   0.100  -1.680  1.00  0.00           C
ATOM     11  N   ALA A   3       5.630   2.380  -0.470  1.00  0.00           N
ATOM     12  CA  ALA A   3       6.330   3.660  -0.460  1.00  0.00           C
ATOM     13  C   ALA A   3       7.830   3.480  -0.460  1.00  0.00           C
ATOM     14  O   ALA A   3       8.360   2.360  -0.470  1.00  0.00           O
ATOM     15  CB  ALA A   3       5.830   4.460   0.760  1.00  0.00           C
HETATM   16  O   HOH A 101      10.000   0.000   0.000  1.00  0.00           O
END
"""

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
    win.deleteLater()
    app.processEvents()


def plugin_menu_item(canvas, title, atom=None, section="atom"):
    """A right-click menu item that a plugin offers for an atom, found by its title (in submenus too)."""
    from miles.ui.canvas.extensions import MenuTarget

    def walk(items):
        for item in items:
            if item.title == title:
                return item
            found = walk(item.items)
            if found is not None:
                return found
        return None

    for function in canvas.extensions.menu_items[section]:
        found = walk(function(canvas, MenuTarget(atom)))
        if found is not None:
            return found
    raise LookupError(title)


def menu_action(window, path):
    """The menu item at a path such as ``"View > Charges > From File"``."""
    from miles.ui.command_palette import menu_entries

    return next(e.action for e in menu_entries(window.menuBar()) if e.path == path)


def canvas_of(window) -> Canvas:
    canvas = window.current_canvas()
    assert isinstance(canvas, Canvas) and isinstance(canvas.surface, PainterSurface)
    return canvas


def screen(canvas, atom: int) -> QPoint:
    x, y = canvas.camera.project(canvas.mol.coords[atom], canvas.width(), canvas.height())
    return QPoint(round(x), round(y))


def click(canvas, pos: QPoint, modifiers=Qt.KeyboardModifier.NoModifier) -> None:
    QTest.mouseClick(canvas, Qt.MouseButton.LeftButton, modifiers, pos)
    QApplication.processEvents()


def drag(canvas, start: QPoint, end: QPoint) -> None:
    def event(kind, pos, buttons):
        p = QPointF(pos)
        return QMouseEvent(kind, p, p, Qt.MouseButton.LeftButton, buttons, Qt.KeyboardModifier.NoModifier)

    canvas.mousePressEvent(event(QMouseEvent.Type.MouseButtonPress, start, Qt.MouseButton.LeftButton))
    for t in np.linspace(0.2, 1.0, 5):
        point = start + (end - start) * float(t)
        canvas.mouseMoveEvent(event(QMouseEvent.Type.MouseMove, point, Qt.MouseButton.LeftButton))
    canvas.mouseReleaseEvent(event(QMouseEvent.Type.MouseButtonRelease, end, Qt.MouseButton.NoButton))
    QApplication.processEvents()


def build_ethane(window) -> Canvas:
    canvas = canvas_of(window)
    click(canvas, QPoint(canvas.width() // 2, canvas.height() // 2))
    hydrogen = next(i for i, a in enumerate(canvas.mol.atoms) if a.symbol == "H")
    click(canvas, screen(canvas, hydrogen))
    return canvas


def test_build_by_clicking(window):
    canvas = build_ethane(window)
    assert canvas.mol.formula() == "C2H6"
    assert canvas.doc.modified
    window.current_canvas().commands.undo()
    assert canvas.mol.formula() == "CH4"
    window.current_canvas().commands.redo()
    assert canvas.mol.formula() == "C2H6"


def test_typing_an_element_changes_the_fragment(window):
    canvas = canvas_of(window)
    QTest.keyClicks(canvas, "n")
    assert window.builder.fragment.symbol == "N"
    QTest.keyClicks(canvas, "c")
    QTest.keyClicks(canvas, "l")
    assert window.builder.fragment.symbol == "Cl"
    QTest.keyClicks(canvas, "2")
    assert window.builder.fragment.symbol == "Cl"


def test_drag_from_an_atom_grows_a_chain(window):
    canvas = build_ethane(window)
    carbon = 0
    start = screen(canvas, carbon)
    drag(canvas, start, start + QPoint(-120, -40))
    assert canvas.mol.formula() == "C3H8"


def test_select_and_edit_a_bond_length(window):
    canvas = build_ethane(window)
    window.builder.set_tool(Tool.SELECT)
    carbons = [i for i, a in enumerate(canvas.mol.atoms) if a.symbol == "C"]
    click(canvas, screen(canvas, carbons[0]))
    click(canvas, screen(canvas, carbons[1]), Qt.KeyboardModifier.ShiftModifier)
    assert canvas.selection == carbons
    panel = window.stage.panel
    assert panel.isVisible() and panel.field.isVisible()
    panel.field.edited.emit(1.7)
    assert np.linalg.norm(canvas.mol.coords[carbons[0]] - canvas.mol.coords[carbons[1]]) == pytest.approx(1.7)
    panel.nudge(5)
    assert panel.dial.value == pytest.approx(1.75)
    window.current_canvas().commands.undo()
    restored = np.linalg.norm(canvas.mol.coords[carbons[0]] - canvas.mol.coords[carbons[1]])
    assert restored == pytest.approx(1.54, abs=0.05)


def test_delete_tool_and_key(window):
    canvas = build_ethane(window)
    window.builder.set_tool(Tool.DELETE)
    hydrogen = next(i for i, a in enumerate(canvas.mol.atoms) if a.symbol == "H")
    click(canvas, screen(canvas, hydrogen))
    assert canvas.mol.formula() == "C2H5"
    window.current_canvas().commands.select_all()
    QTest.keyClick(canvas, Qt.Key.Key_Backspace)
    assert len(canvas.mol) == 0


def test_hover_highlight_depends_on_tool(window):
    canvas = build_ethane(window)
    hydrogen = next(i for i, a in enumerate(canvas.mol.atoms) if a.symbol == "H")
    canvas.hover = AtomHit(hydrogen)
    assert canvas.drawing.highlights()[hydrogen].name == "ATTACH"
    window.builder.set_tool(Tool.DELETE)
    assert canvas.drawing.highlights()[hydrogen].name == "DANGER"


def test_copy_paste_duplicate_and_smiles(window):
    canvas = build_ethane(window)
    window.current_canvas().commands.copy()
    window.current_canvas().commands.paste()
    assert canvas.mol.formula() == "C4H12"
    window.current_canvas().commands.duplicate()
    assert len(canvas.selection) == 8
    assert canvas.mol.formula() == "C6H18"
    window.current_canvas().commands.copy_smiles()
    assert QApplication.clipboard().text() == "CC"


def test_view_menu_actions(window):
    canvas = build_ethane(window)
    for rep in (Representation.TUBE, Representation.WIRE, Representation.SPACEFILL, Representation.BALL):
        window.current_canvas().commands.set_display_mode(rep)
        assert canvas.display.mode == rep
    for labels in LabelMode:
        window.current_canvas().set_style(labels=labels)
        canvas.repaint()
    window.current_canvas().set_style(show_hydrogens=False)
    assert not canvas.display.show_hydrogens
    window.current_canvas().set_style(depth_fade=False)
    for axis in "xyz":
        window.current_canvas().view.look_along(axis)
    window.current_canvas().view.fit(True)
    window.current_canvas().view.reset()
    window.set_appearance("dark")
    assert canvas.theme.dark
    canvas.repaint()


def test_menus_are_connected(window):
    titles = [action.text() for action in window.menuBar().actions()]
    assert titles == ["File", "Edit", "Tools", "View", "Help"]
    assert window.menu_actions.choices["tool"][Tool.BUILD].isChecked()
    window.menu_actions.choices["tool"][Tool.SELECT].trigger()
    assert window.builder.tool == Tool.SELECT


def test_edit_commands(window):
    canvas = build_ethane(window)
    window.current_canvas().commands.clean_geometry()
    window.current_canvas().commands.add_hydrogens()
    window.current_canvas().commands.rebond()
    window.current_canvas().commands.move_to_origin()
    assert np.allclose(canvas.mol.center_of_mass(), 0, atol=1e-6)
    window.dialogs.show_atom_list()
    window.dialogs.show_symmetry()
    assert window.symmetry_popover.group.text() == "D3d"
    window.symmetry_popover.close()


def test_fragment_popover_tabs_and_choice(window):
    popover = window.fragment_popover
    window.dialogs.show_fragments()
    for index in range(popover.pages.count()):
        popover._show_tab(index)
    button, ring = next((b, f) for b, f in popover._tiles if f.name == "Benzene")
    button.click()
    assert window.builder.fragment is ring
    canvas = canvas_of(window)
    click(canvas, QPoint(canvas.width() // 2, canvas.height() // 2))
    assert canvas.mol.formula() == "C6H6"


def test_context_menu_opens_and_closes(window):
    canvas = build_ethane(window)
    QTimer.singleShot(50, lambda: QApplication.activePopupWidget().close())
    QTest.mouseClick(canvas, Qt.MouseButton.RightButton, Qt.KeyboardModifier.NoModifier, screen(canvas, 0))


def test_export_picture(window):
    canvas = build_ethane(window)
    image = canvas.render_export(2.0, Background.TRANSPARENT, ambient_occlusion=False)
    assert image.width() == 2 * canvas.width()
    cropped = crop_image(image, transparent=True)
    assert 0 < cropped.width() < image.width()


def test_protein_opens_as_cartoon(window):
    mol = parse_text(PEPTIDE)
    window.open_structure(mol, JobInfo(title="peptide"))
    canvas = canvas_of(window)
    assert window.builder.tool == Tool.SELECT
    assert canvas.doc.display.rep[1] == Representation.CARTOON
    canvas.repaint()
    canvas.commands.select_group(Group.WATER)
    assert [canvas.mol.atoms[i].resname for i in canvas.selection] == ["HOH"]
    canvas.commands.select_group(Group.LIGANDS)
    assert canvas.selection == []
    canvas.commands.select_group(Group.RESIDUE, 1)
    canvas.commands.set_representation(Representation.BALL)
    canvas.repaint()


def test_mouse_mode_keys(window):
    canvas = build_ethane(window)
    builder = window.builder
    empty = QPoint(40, canvas.height() - 60)
    QTest.keyClicks(canvas, "t")
    assert builder.mouse_mode == MouseMode.TRANSLATE
    rotation = canvas.camera.rotation.copy()
    drag(canvas, empty, empty + QPoint(80, 0))
    assert np.allclose(canvas.camera.rotation, rotation)
    assert canvas.camera.pan[0] > 0
    window.builder.set_tool(Tool.SELECT)
    QTest.keyClicks(canvas, "s")
    assert builder.mouse_mode == MouseMode.SCALE
    scale = canvas.camera.scale
    drag(canvas, empty, empty + QPoint(80, 0))
    assert canvas.camera.scale > scale
    QTest.keyClicks(canvas, "r")
    assert builder.mouse_mode == MouseMode.ROTATE
    assert window.menu_actions.choices["mouse_mode"][MouseMode.ROTATE].isChecked()


def test_mode_letters_still_type_elements(window):
    canvas = canvas_of(window)
    builder = window.builder
    assert builder.tool == Tool.BUILD
    QTest.keyClicks(canvas, "s")
    assert builder.fragment.symbol == "S"
    assert builder.mouse_mode == MouseMode.ROTATE
    window.builder.set_tool(Tool.SELECT)
    QTest.keyClicks(canvas, "si")
    assert builder.fragment.symbol == "Si"
    assert builder.mouse_mode == MouseMode.ROTATE
    QTest.keyClicks(canvas, "br")
    assert builder.fragment.symbol == "Br"
    assert builder.mouse_mode == MouseMode.ROTATE


def test_focus_mode_hides_and_restores(window):
    canvas = build_ethane(window)
    carbons = [i for i, a in enumerate(canvas.mol.atoms) if a.symbol == "C"]
    canvas.set_selection(carbons)
    window.current_canvas().commands.show_only_selection()
    assert canvas.commands.hidden_count() == 6
    assert window.stage.focus.isVisible()
    canvas.repaint()
    window.current_canvas().commands.show_everything()
    assert canvas.commands.hidden_count() == 0
    assert not window.stage.focus.isVisible()


def test_pinned_measurement_and_export(window):
    canvas = build_ethane(window)
    window.builder.set_tool(Tool.SELECT)
    canvas.set_selection([0, 1])
    window.stage.panel.pin.click()
    assert canvas.mol.measurements == [(0, 1)]
    image = canvas.render_export(1.0, Background.WHITE, ambient_occlusion=False)
    assert not image.isNull()
    window.current_canvas().commands.clear_measurements()
    assert canvas.mol.measurements == []


def test_conformers_fill_the_timeline(window):
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CCCCO")
    window.open_structure(mol, JobInfo(title="butanol"))
    canvas = canvas_of(window)
    window.api.run("Generate Conformers")
    window.api.wait()   # it works in the background
    assert canvas.doc.frame_count >= 2
    timeline = window.stage.timeline
    assert timeline.isVisible()
    timeline.step()
    assert canvas.doc.frame == 1
    QTest.keyClick(canvas, Qt.Key.Key_Comma)
    assert canvas.doc.frame == 0


def test_session_keeps_view_and_styles(window, tmp_path):
    canvas = build_ethane(window)
    canvas.camera.turn(0.4, 0.2)
    canvas.set_style(labels=LabelMode.SYMBOL, stereo=True)
    rotation = canvas.camera.rotation.copy()
    path = tmp_path / "ethane.miles"
    canvas.doc.path = str(path)
    assert window.save_document(canvas.doc)
    window.open_path(str(path))
    reopened = canvas_of(window)
    assert reopened is not canvas
    assert np.allclose(reopened.camera.rotation, rotation)
    assert reopened.display.labels == LabelMode.SYMBOL and reopened.display.stereo


def surface_of(canvas):
    """The molecular surface drawn on a canvas (the only mesh in these tests), or None.

    The surface is worked out on a worker thread once the canvas is drawn,
    so this draws the canvas, waits for the work and draws it again.
    """
    canvas.drawing.frame()
    canvas.window().api.wait()
    meshes = canvas.drawing.frame().meshes
    return meshes[0].mesh if meshes else None


def test_surface_is_computed_for_the_view(window):
    canvas = build_ethane(window)
    window.api.run("View > Surface > Show Surface")
    surface = surface_of(canvas)
    assert surface is not None and not surface.empty
    window.api.run("View > Surface > Hydrophobicity")
    assert menu_action(window, "View > Surface > Hydrophobicity").isChecked()
    window.api.run("View > Surface > See-Through")
    assert canvas.drawing.frame().meshes[0].alpha < 1
    window.api.run("View > Surface > Hide Surface")
    assert surface_of(canvas) is None


def test_stereo_labels_and_inversion(window):
    from miles.core.chem import from_smiles
    from miles.core.stereo import Chirality

    mol, _ = from_smiles("N[C@@H](C)C(=O)O")
    window.open_structure(mol, JobInfo(title="alanine"))
    canvas = canvas_of(window)
    window.current_canvas().set_style(stereo=True)
    centers, _ = canvas.stereo_labels()
    (center, label), = centers.items()
    assert label == Chirality.S
    canvas.repaint()
    canvas.commands.invert_center(center)
    assert canvas.stereo_labels()[0][center] == Chirality.R


def test_mutation_and_sequence_strip(window):
    mol = parse_text(PEPTIDE)
    window.open_structure(mol, JobInfo(title="peptide"))
    canvas = canvas_of(window)
    strip = window.stage.strip
    assert strip.isVisible()
    assert [r.letter for r in strip._residues] == ["A", "A", "A"]
    middle_ca = next(i for i, a in enumerate(canvas.mol.atoms) if a.name == "CA" and a.resseq == 2)
    plugin_menu_item(canvas, "PHE", middle_ca).action()
    assert [r.letter for r in strip._residues] == ["A", "F", "A"]
    x = int(strip._residues[0].x + 4)
    QTest.mouseClick(strip, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(x, 15))
    assert {canvas.mol.atoms[i].resseq for i in canvas.selection} == {1}


def test_superpose_two_tabs(window):
    from miles.builtin_plugins.superpose.dialog import superpose_documents
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CC(=O)Oc1ccccc1C(=O)O")
    target = window.api.new(mol.copy(), "target")
    moved = mol.copy()
    moved.coords = moved.coords @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float).T + 5.0
    mobile = window.api.new(moved, "mobile")
    result = superpose_documents(mobile, target, add_copy=True)
    assert result.rmsd < 1e-3
    assert len(target.molecule) == 2 * len(mol)


def test_turntable_frames(window, tmp_path):
    from miles.builtin_plugins.movie.dialog import MovieKind, render_frames

    canvas = build_ethane(window)
    rotation = canvas.camera.rotation.copy()
    written = render_frames(window.api.current, MovieKind.TURNTABLE, 4, tmp_path, 1.0, "white", False)
    assert written == 4 and len(list(tmp_path.glob("frame_*.png"))) == 4
    assert np.allclose(canvas.camera.rotation, rotation)


def test_pocket_view(window):
    ligand = (
        "HETATM   20  C1  ACT B 201       0.000   0.000   3.500  1.00  0.00           C\n"
        "HETATM   21  O1  ACT B 201       0.000   1.100   4.100  1.00  0.00           O\n"
        "HETATM   22  O2  ACT B 201       0.000  -1.100   4.100  1.00  0.00           O\n"
        "HETATM   23  C2  ACT B 201       0.000   0.000   2.000  1.00  0.00           C\n"
    )
    mol = parse_text(PEPTIDE.replace("END\n", ligand + "END\n"))
    window.open_structure(mol, JobInfo(title="complex"))
    canvas = canvas_of(window)
    c1 = next(i for i, a in enumerate(canvas.mol.atoms) if a.resname == "ACT")
    assert plugin_menu_item(canvas, "Mutate ALA 1", 0) is not None
    plugin_menu_item(canvas, "Show Binding Pocket", c1).action()
    alanine = canvas.mol.mask("resname", {"ALA"})
    assert (canvas.doc.display.rep[alanine] == Representation.TUBE).any()
    assert canvas.toast_text.startswith("Pocket of 2 residues")
    canvas.repaint()


def write_cube(path, n=24):
    """Water with an antisymmetric blob around the oxygen, in cube format (bohr)."""
    bohr = 1 / 0.529177
    atoms = [(8, (0.0, 0.0, 0.117)), (1, (0.0, 0.757, -0.469)), (1, (0.0, -0.757, -0.469))]
    origin = np.full(3, -3.0 * bohr)
    step = 6.0 * bohr / (n - 1)
    points = np.indices((n, n, n)).reshape(3, -1).T * step + origin
    centre = np.array(atoms[0][1]) * bohr
    values = (points[:, 2] - centre[2]) * np.exp(-((points - centre) ** 2).sum(axis=1) / 4.0)
    lines = ["test cube", "p-like function", f"{len(atoms):5d}" + "".join(f"{v:12.6f}" for v in origin)]
    lines += [f"{n:5d}" + "".join(f"{step if a == b else 0.0:12.6f}" for b in range(3)) for a in range(3)]
    lines += [f"{z:5d}{float(z):12.6f}" + "".join(f"{c * bohr:12.6f}" for c in xyz) for z, xyz in atoms]
    rows = values.reshape(n * n, n)
    lines += ["".join(f"{v:13.5E}" for v in row[k:k + 6]) for row in rows for k in range(0, n, 6)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_command_palette_runs_a_command(window):
    canvas = build_ethane(window)
    window.dialogs.show_command_palette()
    palette = window.palette_popover
    palette.search.setText("stereo lab")
    assert palette.results.item(0).text() == "View > Stereo Labels"
    QTest.keyClick(palette.search, Qt.Key.Key_Return)
    assert canvas.display.stereo


def test_frames_blend_during_playback(window):
    canvas = build_ethane(window)
    doc = canvas.doc
    start = canvas.mol.coords.copy()
    frames = np.array([start, start + 1.0, start + 2.0])
    doc.set_frames(frames)
    doc.set_position(0.5)
    assert np.allclose(canvas.mol.coords, start + 0.5)
    doc.set_position(2.5)   # a trajectory does not blend its last frame into the first
    assert np.allclose(canvas.mol.coords, start + 2.0)
    doc.set_frames(frames, loop=True)
    doc.set_position(2.5)
    assert np.allclose(canvas.mol.coords, start + 1.0)
    timeline = window.stage.timeline
    timeline.speed.set_index(4)
    assert timeline.rate > 10
    timeline.play()
    assert timeline.playing
    timeline.pause()
    assert not timeline.playing


def test_quantum_output_starts_on_its_last_frame(window):
    from miles.ui.document import Document

    mol = parse_text(PEPTIDE)
    frames = np.array([mol.coords - 1.0, mol.coords.copy()])
    doc = Document(mol, JobInfo(frames=frames))
    assert doc.frame == 1


def test_isosurface_panel(window, tmp_path):
    path = tmp_path / "orbital.cube"
    write_cube(path)
    window.open_path(str(path))
    canvas = canvas_of(window)
    panel = grid_panel(window)
    assert panel.isVisible()
    assert grid_panel(window).request(document_view(window, canvas)) is not None
    assert isosurface(canvas) is not None and not isosurface(canvas).empty
    panel.dial.step(5)
    panel._update()
    assert grid_panel(window).request(document_view(window, canvas)).level > 0
    panel.visible.click()
    assert grid_panel(window).request(document_view(window, canvas)) is None


def test_vibrations_panel_animates_and_restores(window):
    from pathlib import Path

    window.open_path(str(Path(__file__).parent / "data" / "water_freq.log"))
    canvas = canvas_of(window)
    panel = window.stage.floating_panels()["vibrations"]
    assert panel.isVisible() and panel.modes.count() == 3
    before = canvas.mol.coords.copy()
    panel.animate(0)
    assert canvas.doc.frame_count > 1 and canvas.doc.frames_loop
    panel.stop_animation()
    assert canvas.doc.frame_count == 0
    assert np.allclose(canvas.mol.coords, before)


def test_build_from_sequence_dialog(window):
    from miles.builtin_plugins.sequence.dialog import SequenceDialog

    dialog = SequenceDialog(window.api)
    dialog.sequence.setPlainText("ACDB")
    dialog.accept()   # B is no amino acid: the dialog stays open and says so
    assert dialog.result() == 0 and dialog.sequence_text == "" and "B" in dialog.message.text()
    dialog.sequence.setPlainText(">peptide\nACDE")
    dialog.accept()
    assert dialog.result() == 1 and dialog.sequence_text == "ACDE" and dialog.title == "Protein Helix ACDE"
    dialog.deleteLater()


def test_sequence_superpose_and_prepare_commands(window, monkeypatch):
    from miles.builtin_plugins.prepare_protein.dialog import PrepareDialog
    from miles.builtin_plugins.sequence.dialog import SequenceDialog
    from miles.builtin_plugins.superpose.dialog import SuperposeDialog

    def build(dialog):
        dialog.sequence.setPlainText(">peptide\nACDE")
        dialog.accept()
        return dialog.result()

    monkeypatch.setattr(SequenceDialog, "exec", build)
    monkeypatch.setattr(SuperposeDialog, "exec", lambda dialog: 1)
    monkeypatch.setattr(PrepareDialog, "exec", lambda dialog: 1)
    api = window.api
    api.run("Build from Sequence...")
    api.wait()   # it is built in the background
    first = api.current
    assert len({a.resseq for a in first.molecule.atoms}) == 4 and len(api.documents) == 1
    assert first.name == "Protein Helix ACDE" and window.stage.strip.isVisible()
    api.run("Superpose...")   # nothing to superpose onto: a message, nothing changes
    assert len(first.molecule) == len(api.documents[0].molecule)
    api.run("Build from Sequence...")
    api.wait()
    second = api.current
    with second.editing() as mol:
        mol.coords = mol.coords + 5.0
    api.run("Superpose...")   # onto the first tab, with a copy added there
    assert np.allclose(second.molecule.coords, first.molecule.coords[:len(second.molecule)], atol=1e-3)
    assert len(first.molecule) == 2 * len(second.molecule)
    second.undo()
    assert not np.allclose(second.molecule.coords[0], first.molecule.coords[0])

    peptide = api.new(parse_text(PEPTIDE), "peptide")
    count = len(peptide.molecule)
    api.run("Prepare Protein...")
    api.wait()   # it works in the background
    assert len(peptide.molecule) > count   # hydrogens added for pH 7.4
    peptide.undo()
    assert len(peptide.molecule) == count
    api.run("Generate Conformers")   # residues: no conformer search
    assert peptide.frames is None


def test_prepare_and_pocket_commands(window):
    from miles.builtin_plugins.prepare_protein.preparation import prepare

    ligand = (
        "HETATM   20  C1  ACT B 201       0.000   0.000   3.500  1.00  0.00           C\n"
        "HETATM   21  O1  ACT B 201       0.000   1.100   4.100  1.00  0.00           O\n"
        "HETATM   22  O2  ACT B 201       0.000  -1.100   4.100  1.00  0.00           O\n"
    )
    window.open_structure(parse_text(PEPTIDE.replace("END\n", ligand + "END\n")), JobInfo(title="complex"))
    canvas = canvas_of(window)
    report = canvas.doc.edit(prepare, waters=True)
    assert report.hydrogens_added > 0 and report.waters_removed == 1
    window.api.run("Show Binding Pocket")
    alanine = canvas.mol.mask("resname", {"ALA"})
    assert (canvas.doc.display.rep[alanine] == Representation.TUBE).any()


def test_ramachandran_panel(window):
    from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
    from miles.core.internal import set_dihedral

    mol = build_sequence("AGPAAA", SequenceKind.HELIX)
    ala = [i for i, a in enumerate(mol.atoms) if a.resseq == 4]
    n, ca, c = (next(i for i in ala if mol.atoms[i].name == name) for name in ("N", "CA", "C"))
    previous_c = next(i for i, a in enumerate(mol.atoms) if a.resseq == 3 and a.name == "C")
    set_dihedral(mol, previous_c, n, ca, c, -120.0)   # one point apart from the others
    # The helix puts every other residue on the same point, so the residue clicked first gets its own too
    gly = {mol.atoms[i].name: i for i, a in enumerate(mol.atoms) if a.resseq == 2}
    first_c = next(i for i, a in enumerate(mol.atoms) if a.resseq == 1 and a.name == "C")
    set_dihedral(mol, first_c, gly["N"], gly["CA"], gly["C"], -150.0)
    window.open_structure(mol, JobInfo(title="helix"))
    canvas = canvas_of(window)
    panel = window.stage.floating_panels()["ramachandran"]
    assert not panel.isVisible()
    window.api.run("Ramachandran Plot")
    assert panel.isVisible()
    plot = panel.plot
    plotted = [k for k, region in enumerate(plot.regions) if region is not None]
    assert len(plotted) == 4
    first = plotted[0]
    second = next(k for k in plotted if plot.residues[k].resseq == 4)
    assert plot.angles[second][0] == pytest.approx(-120.0, abs=0.5)
    plot.repaint()

    def press(k, modifiers=Qt.KeyboardModifier.NoModifier):
        phi, psi = plot.angles[k]
        QTest.mouseClick(plot, Qt.MouseButton.LeftButton, modifiers, plot._point(phi, psi).toPoint())

    press(first)
    assert set(canvas.selection) == set(plot.residues[first].atoms)
    assert plot.selected == {first}
    press(second, Qt.KeyboardModifier.ShiftModifier)
    assert plot.selected == {first, second}

    moved = canvas.mol.copy()   # the plot follows the frame shown
    set_dihedral(moved, *plot.residues[second].phi_atoms, -80.0)
    canvas.doc.set_frames(np.array([canvas.mol.coords, moved.coords]))
    canvas.doc.set_frame(1)
    assert plot.angles[second][0] == pytest.approx(-80.0, abs=0.5)
    panel.close_panel()
    assert not panel.isVisible()


def test_ramachandran_needs_a_protein(window):
    window.open_item("CCO")
    window.api.run("Ramachandran Plot")
    assert not window.stage.floating_panels()["ramachandran"].isVisible()


def test_energy_plot_follows_the_frames(window):
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CCCC")
    window.open_structure(mol, JobInfo(title="butane"))
    canvas = canvas_of(window)
    doc = canvas.doc
    frames = np.array([mol.coords + 0.01 * k for k in range(5)])
    doc.set_frames(frames, energies=[3.0, 1.0, 0.0, 2.0, 4.0])
    timeline = window.stage.timeline
    plot = timeline.plot
    assert timeline.isVisible() and plot.isVisible()
    assert "\u0394E 3.00 kcal/mol" in timeline.text.text()   # relative to the lowest frame
    QTest.mouseClick(plot, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier,
                     QPoint(plot.width() - 2, plot.height() // 2))
    assert doc.frame == 4
    assert "4.00 kcal/mol" in timeline.text.text()
    plot.grab()   # paints without errors
    doc.set_frames(frames)
    assert timeline.isVisible() and not plot.isVisible()


def test_scan_coordinate_command(window, monkeypatch):
    from miles.builtin_plugins.scan.dialog import ScanDialog
    from miles.core.chem import from_smiles

    def ask(dialog):
        dialog.steps.setValue(6)
        return 1

    monkeypatch.setattr(ScanDialog, "exec", ask)
    mol, _ = from_smiles("CCCC")
    window.open_structure(mol, JobInfo(title="butane"))
    canvas = canvas_of(window)
    before = canvas.mol.coords.copy()
    canvas.set_selection([0])
    window.api.run("Scan Coordinate...")   # one atom is not a coordinate
    assert canvas.doc.frame_count == 0
    canvas.set_selection([0, 1, 2, 3])
    window.api.run("Scan Coordinate...")
    window.api.wait()   # it scans in the background
    doc = canvas.doc
    assert doc.frame_count == 7 and doc.frame_energies is not None
    assert doc.frame_labels is not None and doc.frame_labels[0].endswith("\u00b0")
    assert window.stage.timeline.plot.isVisible()
    window.current_canvas().commands.undo()
    assert np.allclose(canvas.mol.coords, before)


def test_charges_color_labels_and_dipole(window, tmp_path):
    from pathlib import Path

    from miles.core.charges import ChargeModel

    window.open_path(str(Path(__file__).parent / "data" / "water_charges.log"))
    canvas = canvas_of(window)
    assert canvas.doc.charge_model == ChargeModel.FILE
    assert menu_action(window, "View > Charges > From File").isEnabled()
    plain = canvas.drawing.frame().spheres[:, 4:7].copy()
    for command in ("Color by Charge", "Charge Labels", "Dipole Moment"):
        window.api.run(f"View > Charges > {command}")
    colored = canvas.drawing.frame().spheres[:, 4:7]
    assert not np.allclose(plain, colored)
    arrow = canvas.drawing.scene_items().arrows[0]
    assert np.allclose((arrow[3:6] - arrow[0:3]) / np.linalg.norm(arrow[3:6] - arrow[0:3]), [0, 0, -1])
    assert "dipole 2.09 D" in canvas.info_text   # the dipole the calculation reported
    labels = canvas.drawing.annotations()[0].labels
    assert len(labels) == 4 and labels[-1].text == "2.09 D"   # three charges and the dipole
    canvas.repaint()
    image = canvas.render_export(1.0, Background.WHITE, False)
    assert not image.isNull()
    window.api.run("View > Charges > MMFF94")
    charges = canvas.doc.charges()
    assert charges is not None and abs(charges.sum()) < 1e-6
    assert np.array_equal(canvas.mol.column("charge"), charges)   # kept as a property column
    assert "dipole 2.09 D" not in canvas.info_text
    path = tmp_path / "water.miles"
    canvas.doc.path = str(path)
    assert window.save_document(canvas.doc)
    window.open_path(str(path))
    reopened = canvas_of(window)
    assert reopened is not canvas
    assert reopened.doc.charge_model == ChargeModel.MMFF94
    toggles = window.plugins.toggles
    current = window.api.current
    assert all(toggles[name].value(current) for name in ("charge_colors", "charge_labels", "dipole"))
    assert menu_action(window, "View > Charges > Color by Charge").isChecked()
    window.new_document()
    assert not menu_action(window, "View > Charges > From File").isEnabled()
    assert not menu_action(window, "View > Charges > Color by Charge").isChecked()


def test_surface_coloured_by_electrostatic_potential(window):
    from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("[NH3+]CC(=O)[O-]")
    window.open_structure(mol, JobInfo(title="glycine"))
    canvas = canvas_of(window)
    window.api.run("View > Surface > Electrostatic Potential")
    window.api.run("View > Surface > Show Surface")
    surface = surface_of(canvas)
    assert surface is not None and surface.value_range > 0
    assert "Gasteiger charges" in canvas.toast_text and "kcal/(mol e)" in canvas.toast_text
    window.api.run("View > Charges > MMFF94")
    assert surface_of(canvas) is not surface         # new charges, new colours
    assert "MMFF94 charges" in canvas.toast_text
    window.api.run("View > Surface > Element Colors")
    assert surface_of(canvas).value_range == 0

    window.open_structure(build_sequence("KAEAAKD", SequenceKind.HELIX), JobInfo(title="peptide"))
    canvas = canvas_of(window)
    window.api.run("View > Surface > Electrostatic Potential")
    window.api.run("View > Surface > Show Surface")
    assert surface_of(canvas).value_range > 0
    assert "ionisable groups" in canvas.toast_text


def test_reader_notes_are_shown_when_a_file_opens(window, app):
    mol = parse_text("3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n")
    note = "The molecular orbitals were left out (no known normalisation convention fits)."
    window.open_structure(mol, JobInfo(title="water", notes=[note]), "water.molden")
    canvas = window.current_canvas()
    assert canvas.toast_text.startswith("Opened water.molden")
    assert note in canvas.toast_text
    canvas.grab()   # the long message wraps; drawing it must not fail
    app.processEvents()


def test_a_session_keeps_its_own_representations_of_a_protein(window, tmp_path):
    window.open_structure(parse_text(PEPTIDE), JobInfo(title="peptide"))
    canvas = canvas_of(window)
    assert (canvas.doc.display.rep == Representation.CARTOON).any()   # proteins open as a cartoon
    canvas.commands.set_representation(Representation.TUBE, list(range(len(canvas.mol))))
    path = tmp_path / "peptide.miles"
    canvas.doc.path = str(path)
    assert window.save_document(canvas.doc)
    window.open_path(str(path))
    reopened = canvas_of(window)
    assert reopened is not canvas and (reopened.doc.display.rep == Representation.TUBE).all()
