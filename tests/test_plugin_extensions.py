"""What plugins can add besides commands: drawing, status text, menus, sources, fragments and state."""

import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QCoreApplication, QEvent, QSettings
from PySide6.QtWidgets import QApplication, QMenu, QTextBrowser

from miles import plugins
from miles.core.chem import from_smiles
from miles.core.fragments import DATA_FILE
from miles.core.molecule import Representation
from miles.core.volume import Volume
from miles.files import session
from miles.files.common import FileFormatError, JobInfo
from miles.render.picking import AtomHit
from miles.render.theme import Background
from miles.ui.canvas.context_menu import build_context_menu
from miles.ui.canvas.overlay import annotation_tip
from miles.ui.lookup import search_hint
from miles.ui.mainwindow import MainWindow

EXAMPLES = Path(__file__).parent.parent / "examples" / "plugins"

DRAWER = '''
"""Draws on every tab; used by the plugin tests."""
import numpy as np
from miles.plugins import Annotations, Arc, Label, Line, MenuItem, MeshItem, SceneItems
from miles.render.arrows import arrow_rows
from miles.render.cartoon import Mesh
from miles.ui.plugin_panel import PluginPanel

calls = []
GREEN = np.array([0.1, 0.9, 0.1], np.float32)
MESH = Mesh(vertices=np.array([[0, 0, 0], [3, 0, 0], [0, 3, 0]], np.float32),
            normals=np.array([[0, 0, 1]] * 3, np.float32), colors=np.array([[0.9, 0.5, 0.1]] * 3, np.float32),
            indices=np.array([[0, 1, 2]], np.uint32))

def register(context):
    def drawn(doc, drawing):
        calls.append(("annotations", drawing.export))
        return Annotations(lines=[Line(0, 1, "#FF0000", dashes=(4, 2)),
                                  Line(0, [9.0, 9.0, 9.0], drawn_only=True), Line(1, 2, drawn_only=True)],
                           labels=[Label("first", 0, tip="the first atom"),
                                   Label("hidden", 1, drawn_only=True),
                                   Label("point", [0.0, 0.0, 0.0], away_from=0)],
                           arcs=[Arc(0, 1, 2, "#00AA00", tip="an angle")])

    def items(doc, drawing):
        mol = doc.molecule
        if not len(mol):
            return None
        return SceneItems(arrows=arrow_rows([mol.coords[0]], [mol.coords[0] + [0, 0, 2]], 0.1, (1, 0, 0)),
                          meshes=[MeshItem(MESH, 0.5)], atom_colors=np.tile(GREEN, (len(mol), 1)))

    def menu(doc, target):
        if target.atom is None:
            return []
        return [MenuItem(f"Atom {target.atom}", lambda: calls.append(("menu", target.atom))),
                MenuItem("More", items=[MenuItem("Inner", lambda: calls.append(("inner",)))])]

    context.add_annotations(drawn)
    context.add_scene_items(items)
    context.add_status(lambda doc: f"{len(doc.molecule)} atoms counted")
    context.add_context_items(menu, "atom")
    def whole(doc, target):
        return [MenuItem("Whole Structure", lambda: calls.append(("all",)))]

    context.add_context_items(whole, "selection")
    context.add_menu("View > Drawer", after="Surface")
    context.add_command("View > Drawer > First", lambda app: None)
    context.add_separator("View > Drawer")
    context.add_command("View > Drawer > Second", lambda app: None)
    bottom = PluginPanel("Bottom")
    context.add_floating_panel("drawer-bottom", bottom, place="bottom")
'''

BROKEN_LAYER = '''
NAME = "Broken Layer"
calls = []

def register(context):
    def fails(doc, drawing):
        calls.append(1)
        return 1 / 0
    context.add_annotations(fails)
    context.add_status(lambda doc: ["not", "text"])
'''

SOURCES = '''
from miles.core.chem import from_smiles
from miles.files.common import FileFormatError, JobInfo
from miles.plugins import StructureSource

def fetch_code(text):
    if text.upper() == "TOY0":
        raise FileFormatError("TOY0 is not in the toy box")
    return from_smiles("CCO")[0], JobInfo(title=text.upper())

def register(context):
    context.add_structure_source(StructureSource(
        "Toys", fetch_code, matches=lambda text: text.upper().startswith("TOY"), what="toy code",
        example="toy1", message=lambda text: f"Fetching {text.upper()} from the toy box ..."))
    context.add_structure_source(StructureSource(
        "Names", lambda text: (from_smiles("c1ccccc1")[0], JobInfo()), what="name", example="benzene",
        whole=False))
'''


@pytest.fixture
def make_window(app, tmp_path):
    windows = []
    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)

    def make(**kwargs):
        win = MainWindow(settings, **kwargs)
        win.resize(900, 650)
        win.show()
        app.processEvents()
        windows.append(win)
        return win

    yield make
    for win in windows:
        for canvas in win.canvases:
            canvas.doc.modified = False
        win.close()
        win.deleteLater()   # so timers of this test's windows do not fire in the next test
    app.processEvents()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)   # processEvents leaves these


@pytest.fixture
def folder(tmp_path, monkeypatch):
    path = tmp_path / "plugins"
    path.mkdir()
    monkeypatch.setenv(plugins.FOLDER_ENV, str(path))
    monkeypatch.delenv(plugins.NO_PLUGINS_ENV, raising=False)
    return path


@pytest.fixture
def errors(monkeypatch):
    caught = []
    monkeypatch.setattr(sys, "excepthook", lambda kind, value, tb: caught.append(value))
    return caught


def wait_until(condition, timeout=5.0):
    end = time.time() + timeout
    while not condition():
        assert time.time() < end, "timed out"
        QApplication.processEvents()
        time.sleep(0.01)


def titles(menu: QMenu) -> list[str]:
    return [a.text() for a in menu.actions() if not a.isSeparator()]


def test_annotations_scene_items_status_and_menus(make_window, folder, errors):
    (folder / "drawer.py").write_text(DRAWER)
    win = make_window()
    calls = sys.modules["miles_user_plugins.drawer"].calls
    doc = win.api.new(from_smiles("O")[0], "water")   # O, H, H
    canvas = doc.canvas

    # annotations: anchors, drawn_only, labels pushed away from an anchor
    (layer,) = canvas.drawing.annotations()
    first = canvas.drawing.resolve_lines(layer)[0][0]
    assert (first.color, first.dashes) == ("#FF0000", (4, 2))
    assert len(canvas.drawing.resolve_lines(layer)) == 3
    doc.edit_style(lambda mol: setattr(mol.atoms[2], "hidden", True))
    (layer,) = canvas.drawing.annotations()
    assert len(canvas.drawing.resolve_lines(layer)) == 2   # the line to a hidden atom is left out
    assert all(canvas.drawing.resolve_label(layer.labels[k]) is not None for k in (0, 2))
    doc.edit_style(lambda mol: setattr(mol.atoms[1], "hidden", True))
    (layer,) = canvas.drawing.annotations()
    assert canvas.drawing.resolve_label(layer.labels[1]) is None
    calls.clear()
    canvas.drawing.annotations()
    assert calls == []   # kept until the document, the display or the theme changes
    canvas.repaint()

    # scene items: atom colours, arrows and a see-through mesh, on screen and in pictures
    frame = canvas.drawing.frame()
    assert np.allclose(frame.spheres[:, 4:7], [0.1, 0.9, 0.1])
    assert len(frame.arrows) == 1 and frame.arrow_mesh is not None
    assert [item.alpha for item in frame.meshes] == [0.5] and len(frame.layers()) == 1
    picture = canvas.render_export(1.0, Background.WHITE, False)
    assert not picture.isNull() and ("annotations", True) in calls

    # tooltips of labels and arcs
    doc.edit_style(lambda mol: [setattr(atom, "hidden", False) for atom in mol.atoms])
    (layer,) = canvas.drawing.annotations()
    camera = canvas.camera
    x, y = camera.to_screen(camera.to_view(canvas.mol.coords[:1]), canvas.width(), canvas.height())
    assert annotation_tip(canvas, float(x[0]), float(y[0])) == "the first atom"
    assert annotation_tip(canvas, 1.0, 1.0) == ""
    picture_with_arc = canvas.grab()
    assert not picture_with_arc.isNull()

    # status text
    canvas.update_status()
    assert canvas.info_text.endswith("3 atoms counted")

    # right-click menu items for an atom and for the selection
    items = {a.text(): a for a in build_context_menu(canvas, AtomHit(0)).actions()}
    assert "Atom 0" in items and titles(items["More"].menu()) == ["Inner"] and "Whole Structure" in items
    assert list(items).index("Atom 0") < list(items).index("Select")   # atom items before the selection's
    items["Atom 0"].trigger()
    items["More"].menu().actions()[0].trigger()
    assert ("menu", 0) in calls and ("inner",) in calls

    # a submenu placed after View > Surface, with a separator
    view = next(a.menu() for a in win.menuBar().actions() if a.text() == "View")
    names = titles(view)
    assert names[names.index("Surface") + 1] == "Drawer"
    drawer = next(a.menu() for a in view.actions() if a.text() == "Drawer")
    assert [a.isSeparator() for a in drawer.actions()] == [False, True, False]

    # a panel centred at the bottom, which fitting the view keeps clear of
    panel = win.stage.floating_panels()["drawer-bottom"]
    panel.toggle()
    assert abs(panel.x() + panel.width() / 2 - win.stage.width() / 2) <= 1
    assert panel.y() + panel.height() <= win.stage.height() - 40
    assert canvas.view_insets[1] >= panel.height()
    assert errors == []


def test_a_failing_layer_is_reported_once_and_skipped(make_window, folder, errors):
    (folder / "broken_layer.py").write_text(BROKEN_LAYER)
    win = make_window()
    calls = sys.modules["miles_user_plugins.broken_layer"].calls
    doc = win.api.new(from_smiles("CC")[0], "ethane")
    for _ in range(3):
        doc.turn(y=10)
        doc.canvas.repaint()
        doc.canvas.drawing.annotations()
    doc.canvas.update_status()
    assert len(calls) == 1   # not called again after it failed
    doc.edit(lambda mol: None)
    doc.canvas.drawing.annotations()
    assert len(calls) == 2   # but tried again after the next edit
    assert {type(e) for e in errors} == {ZeroDivisionError, TypeError}   # status text must be text
    info = next(i for i in win.plugins.plugins if i.key == "file:broken_layer.py")
    assert "ZeroDivisionError" in info.error or "TypeError" in info.error


def test_structure_sources(make_window, folder, errors):
    (folder / "sources.py").write_text(SOURCES)
    first = make_window()
    first.plugins.set_enabled("miles.builtin_plugins.structure_lookup", False)   # no network in tests
    win = make_window()
    api = win.api
    assert "File > Open from PDB" not in api.commands()
    assert search_hint(win.structure_sources) == "SMILES, toy code or name, e.g. CC(=O)O, toy1, benzene"

    doc = api.fetch("toy1")   # a source with a pattern, whole structures
    assert doc.name == "TOY1" and len(doc.molecule) == 9
    assert len(api.fetch("CCO").molecule) == 9   # SMILES is built by MILES itself
    assert len(api.fetch("benzene").molecule) == 12   # the fallback source
    with pytest.raises(FileFormatError, match="toy box"):
        api.fetch("toy0")

    count = len(api.documents)
    win.open_item("toy2")   # the command line, in the background
    wait_until(lambda: len(api.documents) == count + 1)
    assert api.current.name == "TOY2"
    popover = win.fragment_popover
    win.look_up("toy3")   # the fragment search shows what it is doing
    assert popover.hint.text() == "Fetching TOY3 from the toy box ..."
    wait_until(lambda: api.current.name == "TOY3")
    win.look_up("benzene")   # a fragment to place, not a tab
    wait_until(lambda: win.builder.fragment.name == "benzene")
    popover.close()
    assert errors == []


def test_fragment_libraries(make_window, folder, tmp_path, monkeypatch, capsys):
    bundled = json.loads(DATA_FILE.read_text())
    glycine = dict(bundled["bio"][0], key="gly2", name="Glycine Again")
    library = {"Sugars": [dict(glycine, key="sugar", name="Not a Sugar")], "bio": [glycine]}
    (folder / "sugars.json").write_text(json.dumps(library))
    (folder / "library.py").write_text(
        "from pathlib import Path\n\ndef register(context):\n"
        "    context.add_fragments(Path(__file__).with_name('sugars.json'))\n")
    user = tmp_path / "fragments"
    user.mkdir()
    (user / "mine.json").write_text(json.dumps({"Mine": [dict(glycine, key="mine", name="My Fragment")]}))
    (user / "broken.json").write_text("[1, 2")
    monkeypatch.setenv("MILES_FRAGMENT_DIR", str(user))
    win = make_window()
    assert "broken.json" in capsys.readouterr().err   # reported, not fatal
    popover = win.fragment_popover
    tabs = [b.text() for b in popover.tabs.group.buttons()]
    assert tabs[-2:] == ["Mine", "Sugars"]
    info = next(i for i in win.plugins.plugins if i.key == "file:library.py")
    assert info.extras == ["fragments from sugars.json"]
    expected = (("Sugars", "Not a Sugar"), ("Biomolecules", "Glycine Again"), ("Mine", "My Fragment"))
    for title, name in expected:
        popover._show_tab(tabs.index(title))
        assert name in [f.name for _, f in popover._tiles]


def test_sessions_keep_plugin_state_and_old_sessions_open(make_window, tmp_path):
    win = make_window()
    mol, _ = from_smiles("O")
    grid = Volume("blob", np.full(3, -2.0), np.eye(3) * 0.5, np.random.default_rng(1).random((9, 9, 9)))
    old_view = {"style": {"mode": "tube", "hydrogen_bonds": True, "charge_colors": True, "dipole": False},
                "volume": {"index": 0, "level": 0.4, "see_through": False}}
    data = session.to_dict(mol, JobInfo(title="old", view=old_view, volumes=[grid]))
    path = tmp_path / "old.miles"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        json.dump(data, fh)
    win.open_path(str(path))
    doc = win.api.current
    toggles = win.plugins.toggles
    assert str(doc.display.mode) == Representation.TUBE
    assert toggles["hydrogen_bonds"].value(doc) and toggles["charge_colors"].value(doc)
    assert not toggles["dipole"].value(doc)
    request = win.stage.floating_panels()["grid"].request(doc)
    assert request is not None and request.level == pytest.approx(0.4) and not request.see_through
    assert doc.canvas.drawing.scene_items().meshes   # the isosurface is drawn

    toggles["dipole"].set(doc, True)
    again = tmp_path / "new.miles"
    doc.save(str(again))
    with gzip.open(again, "rt") as fh:
        saved = json.load(fh)["view"]
    assert "hydrogen_bonds" not in saved["style"] and "volume" not in saved
    state = saved["plugins"]
    assert state["miles.builtin_plugins.charges"]["dipole"] is True
    assert state["miles.builtin_plugins.grids"]["volume"]["level"] == pytest.approx(0.4)
    win.open_path(str(again))
    reopened = win.api.current
    assert reopened != doc and toggles["dipole"].value(reopened)
    assert toggles["hydrogen_bonds"].value(reopened)
    assert win.stage.floating_panels()["grid"].request(reopened).level == pytest.approx(0.4)


def test_switching_off_drawing_plugins(make_window, folder, tmp_path):
    first = make_window()
    for name in ("hydrogen_bonds", "grids", "vibrations", "charges"):
        first.plugins.set_enabled(f"miles.builtin_plugins.{name}", False)
    win = make_window()
    commands = win.api.commands()
    assert "View > Hydrogen Bonds" not in commands and "Tools > Vibrations" not in commands
    assert not any(path.startswith("View > Charges") for path in commands)
    assert not {"grid", "vibrations"} & set(win.stage.floating_panels())
    doc = win.api.open(str(Path(__file__).parent / "data" / "water_freq.log"))
    assert doc.canvas.drawing.scene_items().meshes == [] and doc.canvas.drawing.annotations() == []
    with pytest.raises(AttributeError):
        doc.style(hydrogen_bonds=True)
    doc.canvas.grab()


def test_principal_axes_example(make_window, monkeypatch, errors):
    monkeypatch.setenv(plugins.FOLDER_ENV, str(EXAMPLES))
    win = make_window()
    info = next(i for i in win.plugins.plugins if i.key == "file:principal_axes.py")
    assert info.state == plugins.PluginState.LOADED, info.error
    doc = win.api.new(from_smiles("ClC(Cl)Cl")[0], "chloroform")
    canvas = doc.canvas
    assert len(canvas.drawing.scene_items().arrows) == 0
    win.api.run("View > Principal Axes")
    arrows = canvas.drawing.scene_items().arrows
    assert len(arrows) == 3
    directions = arrows[:, 3:6] - arrows[:, 0:3]
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    assert np.allclose(directions @ directions.T, np.eye(3), atol=1e-5)   # orthogonal axes
    assert len(canvas.drawing.annotations()[0].labels) == 3
    assert canvas.info_text.endswith("symmetric top")
    assert not canvas.render_export(1.0, Background.WHITE, False).isNull()
    assert errors == []


SETTINGS = '''
"""Has settings; used by the plugin tests."""
REQUIRES_API = 1
seen = []

def register(context):
    context.setting = context.add_setting("cutoff", "Cutoff", 3.5, minimum=1.0, maximum=5.0, step=0.5,
                                          unit="angstrom", changed=seen.append)
    context.add_setting("label", "Label", "x")
    context.add_setting("count", "Count", 2, minimum=1, maximum=9)
    context.add_setting("enabled", "Enabled", True)
    context.add_command("Tools > Report Cutoff", lambda app: seen.append(("cutoff", context.setting.value)))
'''


def test_settings_are_edited_in_the_settings_window(make_window, folder):
    from PySide6.QtWidgets import QCheckBox, QDoubleSpinBox, QLineEdit, QSpinBox

    (folder / "settings.py").write_text(SETTINGS)
    (folder / "too_new.py").write_text("REQUIRES_API = 99\ndef register(context):\n    pass\n")
    win = make_window()
    seen = sys.modules["miles_user_plugins.settings"].seen
    too_new = next(i for i in win.plugins.plugins if i.key == "file:too_new.py")
    assert too_new.state == plugins.PluginState.FAILED
    assert "version 99 of the plugin interface" in too_new.error

    win.dialogs.show_settings("file:settings.py")
    dialog = win.dialogs.settings_dialog
    form = dialog.form

    def field(kind, title):
        return next(f for f in dialog.findChildren(kind) if form.labelForField(f) is not None
                    and form.labelForField(f).text() == title)

    cutoff = field(QDoubleSpinBox, "Cutoff:")
    assert cutoff.value() == 3.5 and cutoff.suffix() == " angstrom" and cutoff.maximum() == 5.0
    cutoff.setValue(4.0)
    assert seen == [4.0]
    assert field(QLineEdit, "Label:").text() == "x"
    assert field(QSpinBox, "Count:").value() == 2
    assert next(b for b in dialog.findChildren(QCheckBox) if b.text() == "Enabled").isChecked()
    dialog.close()
    win.api.run("Report Cutoff")
    assert ("cutoff", 4.0) in seen
    again = make_window()   # kept with the settings
    again.api.run("Report Cutoff")
    assert sys.modules["miles_user_plugins.settings"].seen == [("cutoff", 4.0)]


def test_a_changed_plugin_file_is_reloaded(make_window, folder):
    from miles.ui.plugins_dialog import PluginsDialog

    path = folder / "greeter.py"
    path.write_text("def register(context):\n    context.add_command('Tools > Greet', lambda app: None)\n")
    win = make_window()
    assert "Tools > Greet" in win.api.commands()
    path.write_text("def register(context):\n"
                    "    context.add_command('Tools > Greet Twice', lambda app: None)\n")
    dialog = PluginsDialog(win.plugins, win.theme, win)
    dialog.table.selectRow([info.key for info in win.plugins.plugins].index("file:greeter.py"))
    dialog.reload_selected()
    dialog.close()
    commands = win.api.commands()
    assert "Tools > Greet Twice" in commands and "Tools > Greet" not in commands

    win.plugins.reload(next(i for i in win.plugins.plugins if i.key == "miles.builtin_plugins.charges"))
    assert win.api.commands().count("View > Charges > Dipole Moment") == 1   # nothing left twice
    view = next(a.menu() for a in win.menuBar().actions() if a.text() == "View")
    assert [a.text() for a in view.actions()].count("Charges") == 1


def test_background_tasks_can_be_cancelled(make_window, folder):
    import threading

    win = make_window()
    api = win.api
    bar = win.stage.busy
    release = threading.Event()
    results, failures = [], []
    task = api.run_in_background(lambda: release.wait(5) and 42, message="Counting", done=results.append)
    assert bar.isVisible() and bar.text.text() == "Counting ..."
    release.set()
    api.wait()
    assert results == [42] and not bar.isVisible() and not task.running

    release.clear()
    task = api.run_in_background(lambda: release.wait(5) and 7, message="Counting", done=results.append)
    bar.cancel.click()
    assert task.cancelled and not bar.isVisible()
    release.set()
    time.sleep(0.05)
    QApplication.processEvents()
    assert results == [42]   # a cancelled task's result is dropped

    api.run_in_background(lambda: 1 / 0, message="Dividing", failed=failures.append)
    api.wait()
    assert failures and "division by zero" in failures[0]


def test_shortcuts_dialog_lists_the_menus(make_window, folder):
    from miles.ui.dialogs import ShortcutsDialog, menu_shortcuts

    win = make_window()
    sections = dict(menu_shortcuts(win.menuBar()))
    assert "Ramachandran Plot" in [command for _, command in sections["Tools"]]
    assert any(command == "Charges > Color by Charge" for _, command in sections.get("View", [])) is False
    win.plugins.set_enabled("miles.builtin_plugins.ramachandran", False)
    sections = dict(menu_shortcuts(win.menuBar()))
    assert "Ramachandran Plot" not in [command for _, command in sections["Tools"]]
    text = ShortcutsDialog(win, win.theme, win.menuBar()).findChild(QTextBrowser).toPlainText()
    assert "Clean Geometry" in text and "Ramachandran" not in text


def test_surface_is_not_rebuilt_while_atoms_are_dragged(make_window):
    win = make_window()
    doc = win.api.new(from_smiles("CCO")[0], "ethanol")
    win.api.run("View > Surface > Show Surface")
    canvas = doc.canvas
    canvas.drawing.scene_items()   # asks for the surface, worked out in the background
    win.api.wait()
    before = canvas.drawing.scene_items().meshes[0].mesh
    doc.edit(lambda mol: setattr(mol, "coords", mol.coords + 0.5))
    during = canvas.drawing.scene_items(dragging=True).meshes[0].mesh
    assert during is before and not win.stage.tasks.tasks   # kept while dragging, nothing started
    assert canvas.drawing.scene_items().meshes[0].mesh is before   # the drag ended: worked out again
    win.api.wait()
    assert canvas.drawing.scene_items().meshes[0].mesh is not before
