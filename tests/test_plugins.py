"""Plugins: discovery, what they can add, failures, switching off, and the built-in ones."""

import importlib.metadata
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QLabel, QMenu

from miles import plugins
from miles.core.chem import from_smiles
from miles.files.common import FileFormatError
from miles.files.formats import REGISTRY, read_file, write_filter
from miles.plugin import Change
from miles.plugins import BUILTIN, PluginState
from miles.ui.command_palette import menu_entries
from miles.ui.mainwindow import MainWindow
from miles.ui.plugins_dialog import PluginsDialog, plugin_details

EXAMPLES = Path(__file__).parent.parent / "examples" / "plugins"
ALIGN = "miles.builtin_plugins.align_frames"
MOVED = {   # features moved out of the core, with the menu item each command follows
    "Tools > Build from Sequence...": "Fragments and SMILES...",
    "Tools > Generate Conformers": "Build from Sequence...",
    "Tools > Scan Coordinate...": "Generate Conformers",
    "Tools > Superpose...": "Align Frames",
    "Tools > Properties": "Superpose...",
    "Tools > Prepare Protein...": "Vibrations",   # after the separator that follows Vibrations
    "Tools > Ramachandran Plot": "Show Binding Pocket",
}

COUNTER = '''
"""Counts atoms; used by the plugin tests."""
NAME = "Counter"
VERSION = "0.3"
events = []

def count(app):
    events.append(("count", len(app.current.molecule)))

def read_toy(path):
    from miles.core.molecule import Atom, Molecule
    from miles.files.common import JobInfo
    symbols = open(path).read().split()
    coords = [[1.5 * i, 0, 0] for i in range(len(symbols))]
    return Molecule.from_atoms([Atom(s) for s in symbols], coords), JobInfo()

def register(context):
    from PySide6.QtWidgets import QLabel
    context.add_command("Tools > Counting > Count Atoms", count, shortcut="Ctrl+Alt+Shift+K")
    context.add_command("Say Hello", lambda app: events.append(("hello",)))
    context.add_command("Tools > Count If Selected", count, enabled=lambda app: bool(app.current.selection))
    context.add_command("Edit > Clean Twice", count, shortcut="Ctrl+L")   # taken by Clean Geometry
    context.add_file_reader("toy", read_toy)
    context.add_panel("Counter Panel", lambda app: QLabel("counting"))
    context.on_document_opened(lambda doc: events.append(("opened", doc.name)))
    context.on_selection_changed(lambda doc: events.append(("selection", tuple(doc.selection))))
    context.on_structure_changed(lambda doc: events.append(("structure", len(doc.molecule))))
    context.on_current_changed(lambda doc: events.append(("current", doc.name)))
'''

WATCHER = '''
"""Records the theme, frame and change events; used by the plugin tests."""
events = []

def register(context):
    from miles.ui.floating import FloatingPanel
    context.on_theme_changed(lambda theme: events.append(("theme", theme.name)))
    context.on_frame_changed(lambda doc: events.append(("frame", doc.frame)))
    context.on_document_changed(lambda doc, change: events.append(("changed", change)))
    context.add_command("Tools > Before Vibrations", lambda app: None, before="Vibrations")
    context.add_command("Tools > After Wins", lambda app: None, after="Vibrations",
                        before="Command Palette...")
    context.add_command("Tools > Neither", lambda app: None, after="Missing", before="Also Missing")
    context.add_floating_panel("watcher", FloatingPanel())
    try:
        context.add_floating_panel("properties", FloatingPanel())
    except ValueError as exc:
        events.append(("refused", str(exc)))
'''

BROKEN = {
    "a_import_error.py": "raise ImportError('missing dependency xyz')\n",
    "b_no_register.py": "NAME = 'Nothing'\n",
    "c_register_error.py": "def register(context):\n    raise RuntimeError('cannot register')\n",
    "d_command_error.py": textwrap.dedent('''
        NAME = "Faulty"
        def register(context):
            context.add_command("Tools > Faulty Command", lambda app: 1 / 0)
        '''),
    "_private.py": "raise RuntimeError('files starting with _ are not plugins')\n",
}


@pytest.fixture(autouse=True)
def clean_registry():
    readers, writers = dict(REGISTRY._readers), dict(REGISTRY._writers)
    yield
    REGISTRY._readers, REGISTRY._writers = readers, writers


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


def find(win, key):
    return next(info for info in win.plugins.plugins if info.key == key)


def titles(menu: QMenu) -> list[str]:
    return [a.text() for a in menu.actions() if not a.isSeparator()]


def top_menu(win, title):
    return next(a.menu() for a in win.menuBar().actions() if a.text() == title)


def test_builtin_plugins_keep_their_menu_places(make_window, folder):
    win = make_window()
    for module in BUILTIN:
        info = find(win, module)
        assert info.builtin and info.state == PluginState.LOADED and info.version
        assert info.commands or info.extras
    places = [("File", "Export Movie...", "Export 2D Drawing..."), ("Edit", "Copy SMILES", "Copy 2D Drawing"),
              ("Edit", "Select All", "Select..."), ("Tools", "Scan Coordinate...", "Align Frames"),
              ("Tools", "Properties", "Vibrations"), ("Tools", "Prepare Protein...", "Show Binding Pocket")]
    places += [("Tools", before, path.split(" > ")[1]) for path, before in MOVED.items()]
    for menu, before, item in places:
        items = titles(top_menu(win, menu))
        assert items[items.index(before) + 1] == item
    commands = win.api.commands()
    for path in ("File > Export 2D Drawing", "Edit > Copy 2D Drawing", "Edit > Select",
                 "Tools > Align Frames", *(path.removesuffix("...") for path in MOVED)):
        assert path in commands
        assert path in [e.path for e in menu_entries(win.menuBar())]   # the command palette lists it
    actions = {e.path: e.action for e in menu_entries(win.menuBar())}
    assert actions["Tools > Properties"].shortcut().toString() == "Ctrl+I"
    assert actions["Tools > Ramachandran Plot"].shortcut().toString() == "Ctrl+Alt+R"
    assert {"ramachandran", "properties"} <= set(win.stage.floating_panels())   # remembered as before


def test_builtin_plugins_are_listed_as_built_in(make_window, folder):
    win = make_window()
    dialog = PluginsDialog(win.plugins, win.theme, win)
    rows = {dialog.table.item(row, 1).text(): row for row in range(dialog.table.rowCount())}
    for name in ("Ramachandran Plot", "Properties", "Scan Coordinate", "Generate Conformers", "Superpose",
                 "Build from Sequence", "Prepare Protein"):
        row = rows[name]
        assert dialog.table.item(row, 3).text() == "built in"
        assert dialog.table.item(row, 4).text() == "loaded"
        assert dialog.table.item(row, 0).checkState() == Qt.CheckState.Checked
    dialog.close()


def test_builtin_entry_points_are_declared_but_not_loaded_twice(make_window, folder):
    declared = {ep.value for ep in importlib.metadata.entry_points(group=plugins.ENTRY_POINT_GROUP)}
    assert set(BUILTIN) <= declared
    win = make_window()
    assert [info.key for info in win.plugins.plugins] == list(BUILTIN)
    assert win.api.commands().count("Tools > Align Frames") == 1


def test_folder_plugin_adds_commands_panel_reader_and_hooks(make_window, folder, tmp_path, errors):
    (folder / "counter.py").write_text(COUNTER, encoding="utf-8")
    win = make_window()
    info = find(win, "file:counter.py")
    assert info.state == PluginState.LOADED and (info.name, info.version) == ("Counter", "0.3")
    assert info.description == "Counts atoms; used by the plugin tests."
    events = sys.modules["miles_user_plugins.counter"].events
    api = win.api

    # commands: in menus (a new submenu, a new top-level menu before Help), the palette and App.run
    assert "Tools > Counting > Count Atoms" in api.commands() and "Plugins > Say Hello" in api.commands()
    assert [a.text() for a in win.menuBar().actions()][-2:] == ["Plugins", "Help"]
    doc = api.new(from_smiles("CCO")[0], "ethanol")
    api.run("Count Atoms")
    api.run("Say Hello")
    assert ("count", 9) in events and ("hello",) in events
    actions = {e.path: e.action for e in menu_entries(win.menuBar())}
    assert actions["Tools > Counting > Count Atoms"].shortcut().toString() == "Ctrl+Alt+Shift+K"
    assert "already used" in info.warnings and not info.error   # Ctrl+L belongs to Clean Geometry
    assert actions["Edit > Clean Twice"].shortcut().isEmpty()

    # a command that is not available cannot be run
    events.clear()
    with pytest.raises(RuntimeError, match="not available"):
        api.run("Count If Selected")
    assert events == []
    doc.select([0])
    api.run("Count If Selected")
    assert ("count", 9) in events

    # hooks
    events.clear()
    doc.select([1, 2])
    doc.edit(lambda m: m.remove_atoms([8]))
    assert ("selection", (1, 2)) in events and ("structure", 8) in events
    toy = tmp_path / "three.toy"
    toy.write_text("C O C", encoding="utf-8")
    win.open_path(str(toy))
    assert ("opened", "three.toy") in events and ("current", "three.toy") in events
    assert [a.symbol for a in api.current.molecule.atoms] == ["C", "O", "C"]
    events.clear()
    doc.activate()
    assert events == [("current", "ethanol")]

    # the panel: hidden at first, shown and hidden by its View command, arranged like the others
    panel = win.stage.extra_panels["plugin-file_counter.py-counter_panel"]
    assert not panel.isVisible()
    api.run("View > Counter Panel")
    assert panel.isVisible() and panel.parent() is win.stage
    api.run("View > Counter Panel")
    assert not panel.isVisible()
    assert errors == []


def test_theme_frame_and_change_events(make_window, folder, errors):
    (folder / "watcher.py").write_text(WATCHER, encoding="utf-8")
    win = make_window()
    assert find(win, "file:watcher.py").state == PluginState.LOADED
    events = sys.modules["miles_user_plugins.watcher"].events
    assert events[0] == ("theme", win.theme.name)   # called at once, with the current theme
    assert events[1][0] == "refused" and "properties" in events[1][1]
    assert "watcher" in win.stage.floating_panels()

    events.clear()
    win.set_appearance("dark")
    win.set_appearance("light")
    assert events == [("theme", "dark"), ("theme", "light")]

    doc = win.api.new(from_smiles("CCCC")[0], "butane")
    events.clear()
    frames = np.array([doc.molecule.coords + 0.1 * k for k in range(3)])
    doc.set_frames(frames)   # replaced frames, then the first shown, as one undoable step
    assert ("frame", 0) in events and ("changed", "frames") in events and ("changed", "coords") in events
    events.clear()
    doc.frame = 2
    assert events == [("frame", 2), ("changed", "frame")]
    events.clear()
    doc.select([0])
    doc.canvas.commands.hide_selection()
    doc.edit(lambda m: m.remove_atoms([13]))
    changes = [change for kind, change in events if kind == "changed"]
    assert changes == ["style", "frames", "structure"]   # removing an atom drops the frames
    assert [kind for kind, _ in events].count("frame") == 1
    assert all(type(change) is Change for change in changes)   # the enum, which compares equal to its text

    tools = titles(top_menu(win, "Tools"))
    assert tools[tools.index("Before Vibrations") + 1] == "Vibrations"
    assert tools[tools.index("Vibrations") + 1] == "After Wins"
    assert tools[-1] == "Neither"
    assert errors == []


def test_failures_are_recorded_and_reported_once(make_window, folder, errors):
    for name, source in BROKEN.items():
        (folder / name).write_text(source, encoding="utf-8")
    (folder / "counter.py").write_text(COUNTER, encoding="utf-8")
    win = make_window()
    states = {info.key: info.state for info in win.plugins.plugins}
    assert states["file:a_import_error.py"] == PluginState.FAILED
    assert states["file:b_no_register.py"] == PluginState.FAILED
    assert states["file:c_register_error.py"] == PluginState.FAILED
    assert states["file:d_command_error.py"] == PluginState.LOADED
    assert states["file:counter.py"] == PluginState.LOADED   # later plugins still load
    assert "file:_private.py" not in states
    assert "missing dependency xyz" in find(win, "file:a_import_error.py").error
    assert "no register(context)" in find(win, "file:b_no_register.py").error
    assert "cannot register" in find(win, "file:c_register_error.py").error

    win.api.run("Faulty Command")   # the error goes to the error dialog; the window carries on
    assert len(errors) == 1 and isinstance(errors[0], ZeroDivisionError)
    assert "ZeroDivisionError" in find(win, "file:d_command_error.py").error
    win.api.run("Clean Geometry")

    dialog = PluginsDialog(win.plugins, win.theme, win)
    assert dialog.table.rowCount() == len(win.plugins.plugins)
    assert "missing dependency xyz" in dialog.details.toPlainText()   # the first failure is selected
    assert "Commands:\n  Tools > Faulty Command" in plugin_details(find(win, "file:d_command_error.py"))
    dialog.close()


def test_a_failing_file_reader_gives_a_message(make_window, folder, tmp_path, monkeypatch):
    (folder / "bad_reader.py").write_text(
        "def register(context):\n    context.add_file_reader('bad', lambda path: {}['missing'])\n",
        encoding="utf-8")
    win = make_window()
    path = tmp_path / "x.bad"
    path.write_text("", encoding="utf-8")
    with pytest.raises(FileFormatError, match="bad_reader plugin could not handle"):
        read_file(str(path))
    warnings = []
    monkeypatch.setattr("PySide6.QtWidgets.QMessageBox.warning", lambda *args: warnings.append(args))
    win.open_path(str(path))
    assert len(warnings) == 1 and len(win.canvases) == 1


def test_entry_point_plugins(make_window, folder, tmp_path, monkeypatch):
    package = tmp_path / "site"
    package.mkdir()
    (package / "mb_extra.py").write_text(textwrap.dedent('''
        """Adds a greeting."""
        VERSION = "2.1"
        def register(context):
            context.add_command("Tools > Greet", lambda app: app.message("Hello"))
        '''), encoding="utf-8")
    monkeypatch.syspath_prepend(str(package))
    group = plugins.ENTRY_POINT_GROUP
    points = [importlib.metadata.EntryPoint("extra", "mb_extra", group),
              importlib.metadata.EntryPoint("extra_function", "mb_extra:register", group),
              importlib.metadata.EntryPoint("broken", "mb_missing_module", group),
              importlib.metadata.EntryPoint("align_frames", ALIGN, group)]
    monkeypatch.setattr(importlib.metadata, "entry_points", lambda group: points)
    win = make_window()
    extra = find(win, "entry:extra")
    assert extra.state == PluginState.LOADED and extra.version == "2.1"
    assert extra.description == "Adds a greeting."
    assert extra.source == "entry point mb_extra"
    assert find(win, "entry:extra_function").state == PluginState.LOADED   # names register directly
    assert find(win, "entry:broken").state == PluginState.FAILED
    assert win.api.commands().count("Tools > Greet") == 2
    assert [i.key for i in win.plugins.plugins].count(ALIGN) == 1


def test_no_plugins_skips_third_party_ones(make_window, folder, monkeypatch):
    (folder / "counter.py").write_text(COUNTER, encoding="utf-8")
    win = make_window(plugins=False)
    assert find(win, "file:counter.py").state == PluginState.SKIPPED
    assert "Tools > Counting > Count Atoms" not in win.api.commands()
    assert "Tools > Align Frames" in win.api.commands()   # the built-in plugins still load
    assert all(path.removesuffix("...") in win.api.commands() for path in MOVED)
    assert all(find(win, module).state == PluginState.LOADED for module in BUILTIN)
    monkeypatch.setenv(plugins.NO_PLUGINS_ENV, "1")
    assert find(make_window(), "file:counter.py").state == PluginState.SKIPPED


def test_switching_plugins_off_and_on_takes_effect_at_once(make_window, folder):
    (folder / "counter.py").write_text(COUNTER, encoding="utf-8")
    win = make_window()
    dialog = PluginsDialog(win.plugins, win.theme, win)
    rows = [info.key for info in win.plugins.plugins]
    dialog.table.item(rows.index("file:counter.py"), 0).setCheckState(Qt.CheckState.Unchecked)
    dialog.table.item(rows.index(ALIGN), 0).setCheckState(Qt.CheckState.Unchecked)
    assert dialog.table.item(rows.index(ALIGN), 4).text() == "switched off"
    dialog.close()
    assert win.plugins.disabled == {"file:counter.py", ALIGN}
    commands = win.api.commands()
    assert "Tools > Align Frames" not in commands and "Plugins > Say Hello" not in commands
    assert [a.text() for a in win.menuBar().actions()][-2:] == ["View", "Help"]   # its own menu went too
    assert "plugin-file_counter.py-counter_panel" not in win.stage.floating_panels()
    win.plugins.set_enabled("file:counter.py", True)
    assert "Plugins > Say Hello" in win.api.commands()
    win.plugins.set_enabled("file:counter.py", False)

    again = make_window()
    assert find(again, "file:counter.py").state == PluginState.OFF
    assert find(again, ALIGN).state == PluginState.OFF
    assert "Tools > Align Frames" not in again.api.commands()
    again.plugins.set_enabled(ALIGN, True)
    assert "Tools > Align Frames" in make_window().api.commands()


def test_switching_a_moved_feature_off_removes_its_commands_and_panel(make_window, folder):
    win = make_window()
    win.api.run("Ramachandran Plot")
    for module in ("ramachandran", "scan", "properties"):
        win.plugins.set_enabled(f"miles.builtin_plugins.{module}", False)

    for window in (win, make_window()):   # gone at once, and after a restart
        commands = window.api.commands()
        for path in ("Tools > Ramachandran Plot", "Tools > Scan Coordinate", "Tools > Properties"):
            assert path not in commands
        assert not {"ramachandran", "properties"} & set(window.stage.floating_panels())
        tools = titles(top_menu(window, "Tools"))
        assert tools[tools.index("Generate Conformers") + 1] == "Align Frames"   # the others keep places
        assert tools[tools.index("Superpose...") + 1] == "Vibrations"
        shortcuts = {e.action.shortcut().toString() for e in menu_entries(window.menuBar())}
        assert not {"Ctrl+I", "Ctrl+Alt+R"} & shortcuts
        window.api.new(from_smiles("CCO")[0], "ethanol")   # the core works without them
        window.api.run("Clean Geometry")
    win.plugins.set_enabled("miles.builtin_plugins.ramachandran", True)   # back at once
    assert "Tools > Ramachandran Plot" in win.api.commands()
    assert "ramachandran" in win.stage.floating_panels()


def test_example_plugin(make_window, tmp_path, monkeypatch, errors):
    monkeypatch.setenv(plugins.FOLDER_ENV, str(EXAMPLES))
    win = make_window()
    info = find(win, "file:selection_info.py")
    assert info.state == PluginState.LOADED, info.error
    api = win.api
    doc = api.new(from_smiles("CCO")[0], "ethanol")
    panel = win.stage.extra_panels["plugin-file_selection_info.py-selection_info"]
    label = panel.body.findChild(QLabel)
    assert "All 9 atoms" in label.text() and "C2H6O" in label.text()
    api.run("Select Heavy Atoms")
    assert len(doc.selection) == 3 and "3 selected atoms" in label.text()
    assert "C2O" in label.text()

    api.run("Move Selection to Origin")
    heavy = doc.molecule.subset(doc.selection)
    assert np.allclose(heavy.center_of_mass(), 0.0, atol=1e-9)
    doc.undo()
    assert not np.allclose(doc.molecule.subset(doc.selection).center_of_mass(), 0.0)

    path = tmp_path / "ethanol.atoms"
    assert "ATOMS (*.atoms)" in write_filter()
    doc.save(str(path))
    assert path.read_text(encoding="utf-8").count("\n") >= 9
    win.open_path(str(path))
    assert api.current.molecule.formula() == "C2H6O" and len(api.current.molecule.bonds) == 8
    assert errors == []


def test_plugin_details_show_warnings_apart_from_errors():
    from miles.plugins import PluginInfo, PluginState

    info = PluginInfo(key="file:demo.py", name="Demo", source="plugin folder", state=PluginState.LOADED,
                      warnings="Shortcut Ctrl+K of Demo is already used; the command has none.\n")
    text = plugin_details(info)
    assert "Warnings:\n  Shortcut Ctrl+K of Demo is already used" in text
    assert "Error:" not in text
