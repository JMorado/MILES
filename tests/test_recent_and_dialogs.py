"""File > Open Recent, and the remaining dialogs drawn in the app's theme."""

from pathlib import Path

import pytest
from PySide6.QtCore import QSettings
from PySide6.QtWidgets import QComboBox, QMessageBox, QTextBrowser

from miles.builtin_plugins.superpose.dialog import SuperposeDialog
from miles.render.theme import DARK, LIGHT
from miles.ui.dialogs import AtomListDialog, ShortcutsDialog
from miles.ui.document import Document
from miles.ui.mainwindow import RECENT_LIMIT, MainWindow

EXAMPLE = Path(__file__).parent.parent / "examples" / "alanine_dipeptide_md.xyz"


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


def write_xyz(path):
    path.write_text("3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n", encoding="utf-8")
    return str(path)


def menu_texts(window):
    return [a.text() for a in window.menu_actions.menus["recent"].actions() if not a.isSeparator()]


def test_recent_files_are_listed_most_recent_first(window, tmp_path):
    assert window.files.recent_files() == []
    assert menu_texts(window) == ["No Recent Files", "Clear Menu"]
    first = write_xyz(tmp_path / "first.xyz")
    second = write_xyz(tmp_path / "second.xyz")
    window.open_path(first)
    window.open_path(second)
    window.open_path(first)
    assert window.files.recent_files() == [first, second]
    assert menu_texts(window) == ["first.xyz", "second.xyz", "Clear Menu"]


def test_saved_files_are_remembered_and_the_list_is_capped(window, tmp_path):
    for k in range(RECENT_LIMIT + 3):
        window.open_path(write_xyz(tmp_path / f"m{k}.xyz"))
    assert len(window.files.recent_files()) == RECENT_LIMIT
    doc = window.current_doc()
    doc.path = str(tmp_path / "saved.xyz")
    assert window.save_document(doc)
    assert window.files.recent_files()[0] == str(tmp_path / "saved.xyz")


def test_recent_entry_opens_and_missing_files_are_dropped(window, tmp_path, monkeypatch):
    warnings = []
    monkeypatch.setattr(QMessageBox, "warning", lambda *args: warnings.append(args[1]))
    kept = write_xyz(tmp_path / "kept.xyz")
    gone = write_xyz(tmp_path / "gone.xyz")
    window.open_path(kept)
    window.open_path(gone)
    Path(gone).unlink()
    tabs = len(window.canvases)
    recent = window.menu_actions.menus["recent"]
    recent.actions()[1].trigger()   # kept.xyz
    assert len(window.canvases) == tabs + 1
    window.files.open_recent(gone)
    assert warnings == ["File not found"]
    assert window.files.recent_files() == [kept]


def test_clear_menu(window, tmp_path):
    window.open_path(write_xyz(tmp_path / "a.xyz"))
    window.files.clear_recent()
    assert window.files.recent_files() == []
    assert menu_texts(window) == ["No Recent Files", "Clear Menu"]


def test_recent_files_survive_a_restart(app, tmp_path):
    settings = str(tmp_path / "settings.ini")
    first = MainWindow(QSettings(settings, QSettings.Format.IniFormat))
    path = write_xyz(tmp_path / "w.xyz")
    first.open_path(path)
    first.settings.sync()
    for canvas in first.canvases:
        canvas.doc.modified = False
    first.close()
    second = MainWindow(QSettings(settings, QSettings.Format.IniFormat))
    assert second.files.recent_files() == [path]
    assert menu_texts(second)[0] == "w.xyz"
    second.close()


@pytest.mark.parametrize("theme", [LIGHT, DARK])
def test_remaining_dialogs_use_the_app_theme(app, theme):
    doc = Document()
    atoms = AtomListDialog(doc, None, theme)
    assert theme.panel.name() in atoms.styleSheet()
    shortcuts = ShortcutsDialog(None, theme)
    assert theme.panel.name() in shortcuts.styleSheet()
    assert shortcuts.findChild(QTextBrowser) is not None


@pytest.mark.parametrize("appearance", ["light", "dark"])
def test_superpose_dialog_uses_the_app_theme(window, appearance):
    window.set_appearance(appearance)
    first = window.api.new(window.api.current.molecule.copy() if window.api.current else None)
    second = window.api.new()
    superpose = SuperposeDialog(window.api, first, [second])
    theme = window.theme
    assert theme.panel.name() in superpose.styleSheet()
    combo = superpose.findChild(QComboBox)
    assert combo.palette().color(combo.palette().ColorRole.ButtonText) == theme.text
    superpose.deleteLater()


def test_open_atom_lists_follow_the_appearance(window):
    window.open_path(str(EXAMPLE))
    window.dialogs.show_atom_list()
    dialog = window.dialogs._atom_lists[id(window.current_doc())]
    window.set_appearance("dark")
    assert DARK.panel.name() in dialog.styleSheet()
    window.set_appearance("light")
    assert LIGHT.panel.name() in dialog.styleSheet()
    dialog.close()
