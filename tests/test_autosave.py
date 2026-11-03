"""Autosave, recovery after a crash, and the Settings window."""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
from PySide6.QtCore import QSettings
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import QCheckBox, QComboBox, QDialog, QDoubleSpinBox

from miles import plugins
from miles.core.molecule import Representation
from miles.render.scene import LabelMode
from miles.ui import autosave
from miles.ui.mainwindow import MainWindow
from miles.ui.recovery_dialog import DISCARD, RecoveryDialog, when
from miles.ui.settings_dialog import FIELD_WIDTH, SettingsDialog

TRAJECTORY = Path(__file__).resolve().parent.parent / "examples" / "alanine_dipeptide_md.xyz"


@pytest.fixture
def folder(tmp_path, monkeypatch):
    path = tmp_path / "autosave"
    monkeypatch.setenv(autosave.FOLDER_ENV, str(path))
    return path


@pytest.fixture
def make_window(app, tmp_path, folder):
    windows = []
    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)

    def make():
        win = MainWindow(settings, plugins=False)
        win.resize(900, 650)
        win.show()
        app.processEvents()
        windows.append(win)
        return win

    yield make
    for win in windows:
        if win.dialogs.recovery_dialog is not None:
            win.dialogs.recovery_dialog.reject()
        for canvas in win.canvases:
            canvas.doc.modified = False
        win.close()
        win.deleteLater()
    app.processEvents()


def instance_folders(folder: Path) -> list[Path]:
    return sorted(p for p in folder.iterdir() if p.is_dir()) if folder.is_dir() else []


def edited_document(win: MainWindow):
    """Open the example trajectory, show frame 5 and hide three atoms (an edit)."""
    win.open_path(str(TRAJECTORY))
    canvas = win.current_canvas()
    canvas.doc.set_frame(5)
    canvas.set_selection([0, 1, 2])
    canvas.commands.hide_selection()
    canvas.camera.scale *= 1.7
    assert canvas.doc.modified
    return canvas


def simulate_crash(win: MainWindow) -> None:
    """Leave the window's copies behind as a crash of another process would: unlocked, never deleted."""
    saver = win.autosave
    saver._lock.unlock()
    saver._lock = None
    saver._closed = True
    saver.timer.stop()
    gone = saver.folder.with_name("999999999-" + saver.folder.name.partition("-")[2])   # no such process
    saver.folder.rename(gone)
    saver.folder = gone


def test_only_modified_documents_changed_since_their_last_copy_are_written(make_window, folder):
    win = make_window()
    assert win.autosave.save_now() == 0   # the empty untitled tab is not modified
    assert not folder.exists()            # nothing is created until there is something to keep
    canvas = edited_document(win)
    assert win.autosave.save_now() == 1
    win.autosave.wait()
    [own] = instance_folders(folder)
    assert own.name.startswith(f"{os.getpid()}-")
    assert sorted(p.suffix for p in own.iterdir() if p.name != autosave.LOCK_NAME) == [".json", ".miles"]
    assert win.autosave.save_now() == 0   # unchanged since
    canvas.doc.set_frame(7)               # showing another frame is not an edit
    assert win.autosave.save_now() == 0
    canvas.set_selection([3])
    canvas.commands.hide_selection()
    assert win.autosave.save_now() == 1
    win.autosave.wait()


def test_saving_closing_and_quitting_delete_the_copies(make_window, folder, tmp_path, monkeypatch):
    win = make_window()
    canvas = edited_document(win)
    win.new_document()
    other = win.current_canvas()
    other.commands.paste_molecule(canvas.mol.subset([3, 4, 5]))
    assert win.autosave.save_now() == 2
    win.autosave.wait()
    [own] = instance_folders(folder)
    assert len(list(own.glob("*.miles"))) == 2

    canvas.doc.path = str(tmp_path / "saved.xyz")   # not over the example file
    assert win.save_document(canvas.doc)
    win.autosave.wait()
    assert len(list(own.glob("*.miles"))) == 1               # saved: its copy is gone

    monkeypatch.setattr(win.files, "confirm_discard", lambda doc: True)   # "Don't Save"
    win.close_tab(win.canvases.index(other))
    win.autosave.wait()
    assert not list(own.glob("*.miles"))

    canvas.set_selection([4])
    canvas.commands.hide_selection()
    assert win.autosave.save_now() == 1
    win.close()   # quitting normally
    assert not instance_folders(folder)


def test_a_crash_is_offered_once_and_restores_atoms_frames_and_view(make_window, folder):
    crashed = make_window()
    canvas = edited_document(crashed)
    crashed.autosave.save_now()
    crashed.autosave.wait()
    simulate_crash(crashed)

    win = make_window()   # the next start
    dialog = win.dialogs.recovery_dialog
    assert isinstance(dialog, RecoveryDialog)
    assert [leftover.name for leftover in dialog.leftovers] == [TRAJECTORY.name]
    assert dialog.table.item(0, 1).text().startswith("Today")
    assert make_window().dialogs.recovery_dialog is None   # a third window does not offer them as well

    dialog.accept()
    restored = win.current_canvas()
    doc = restored.doc
    assert doc.path == str(TRAJECTORY) and doc.modified and doc.name == TRAJECTORY.name
    assert doc.mol.symbols == canvas.mol.symbols
    assert np.allclose(doc.mol.coords, canvas.mol.coords, atol=1e-4)
    assert dict(doc.mol.bonds) == dict(canvas.mol.bonds)
    assert np.array_equal(doc.display.hidden, canvas.doc.display.hidden) and doc.display.hidden.any()
    assert doc.frame_count == canvas.doc.frame_count and doc.frame == 5
    assert np.allclose(doc.frames, canvas.doc.frames, atol=1e-4)
    assert restored.camera.scale == pytest.approx(canvas.camera.scale, rel=1e-6)
    assert np.allclose(restored.camera.rotation, canvas.camera.rotation)
    win.autosave.wait()
    assert not crashed.autosave.folder.exists()          # the old copies are gone ...
    assert list(win.autosave.folder.glob("*.miles"))       # ... and the restored document has its own


def test_copies_of_this_process_are_never_offered(make_window, folder):
    """A window dropped without closing leaves its copies, but the process still runs: no crash."""
    dropped = make_window()
    edited_document(dropped)
    dropped.autosave.save_now()
    dropped.autosave.wait()
    dropped.autosave._lock.unlock()   # as when the window is deleted without being closed
    assert make_window().dialogs.recovery_dialog is None
    assert list(dropped.autosave.folder.glob("*.miles"))


def test_discard_deletes_the_copies_and_closing_keeps_them(make_window, folder):
    crashed = make_window()
    edited_document(crashed)
    crashed.autosave.save_now()
    crashed.autosave.wait()
    simulate_crash(crashed)

    later = make_window()
    later.dialogs.recovery_dialog.reject()   # closed without choosing
    assert crashed.autosave.folder.exists()

    win = make_window()
    assert win.dialogs.recovery_dialog is not None
    win.dialogs.recovery_dialog.done(DISCARD)
    assert not crashed.autosave.folder.exists()
    assert len(win.canvases) == 1 and not len(win.canvases[0].mol)
    assert make_window().dialogs.recovery_dialog is None


def test_a_running_miles_keeps_its_copies_until_it_stops(make_window, folder, tmp_path, app):
    """A second MILES does not take the copies of one still running, only once it is killed."""
    script = tmp_path / "other.py"
    script.write_text(textwrap.dedent(f"""
        import sys
        from PySide6.QtCore import QSettings
        from PySide6.QtWidgets import QApplication
        app = QApplication([])
        from miles.ui.mainwindow import MainWindow
        win = MainWindow(QSettings(sys.argv[1], QSettings.Format.IniFormat), plugins=False)
        win.open_path({str(TRAJECTORY)!r})
        canvas = win.current_canvas()
        canvas.set_selection([0])
        canvas.commands.hide_selection()
        win.autosave.save_now()
        win.autosave.wait()
        print("ready", flush=True)
        sys.stdin.readline()
    """))
    env = dict(os.environ, **{autosave.FOLDER_ENV: str(folder), plugins.NO_PLUGINS_ENV: "1"})
    other = subprocess.Popen([sys.executable, str(script), str(tmp_path / "other.ini")], env=env,
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    try:
        for line in other.stdout:   # anything printed before is not the signal
            if line.strip() == "ready":
                break
        win = make_window()
        assert win.dialogs.recovery_dialog is None
        assert win.autosave.find_leftovers() == []
    finally:
        other.kill()   # a crash: no clean quit
        other.wait()
    [leftover] = win.autosave.find_leftovers()
    assert leftover.path == str(TRAJECTORY)
    assert win.files.restore_documents([leftover]) == 1
    assert win.current_canvas().doc.display.hidden[0]
    assert not [p for p in instance_folders(folder) if p != win.autosave.folder]


def test_unreadable_copies_are_kept_for_another_try(make_window, folder, monkeypatch):
    crashed = make_window()
    edited_document(crashed)
    crashed.autosave.save_now()
    crashed.autosave.wait()
    simulate_crash(crashed)
    [copy] = crashed.autosave.folder.glob("*.miles")
    copy.write_bytes(b"not a session")
    win = make_window()
    leftovers = win.dialogs.recovery_dialog.leftovers
    shown = []
    monkeypatch.setattr("PySide6.QtWidgets.QMessageBox.warning", lambda *args: shown.append(args))
    win.dialogs.recovery_dialog.accept()
    assert shown and leftovers[0].file.exists()


def test_empty_folders_are_cleared_only_once_they_are_old(make_window, folder):
    from PySide6.QtCore import QLockFile

    win = make_window()
    young, old = folder / "999999998-young", folder / "999999999-old"
    young.mkdir(parents=True)
    old.mkdir()
    os.utime(old, (0, 0))
    assert win.autosave.find_leftovers() == []
    assert young.exists() and not old.exists()   # the young one may be a window just starting

    taken = QLockFile(str(win.autosave.folder / autosave.LOCK_NAME))   # as if another window got there first
    win.autosave.folder.mkdir(parents=True)
    assert taken.tryLock(0)
    first = win.autosave.folder
    edited_document(win)
    assert win.autosave.save_now() == 1
    win.autosave.wait()
    assert win.autosave.folder != first and list(win.autosave.folder.glob("*.miles"))
    taken.unlock()


def test_switching_autosave_off(make_window, folder):
    win = make_window()
    edited_document(win)
    win.autosave.save_now()
    win.autosave.wait()
    win.autosave.set_interval(0)
    win.autosave.wait()
    assert not win.autosave.timer.isActive()
    assert not list(win.autosave.folder.glob("*.miles"))
    assert win.autosave.save_now() == 0
    win.autosave.set_interval(5)
    assert win.autosave.timer.interval() == 5 * 60_000
    assert win.settings.value(autosave.INTERVAL_KEY) in (5, "5")


def test_autosave_times_read_naturally():
    import datetime

    now = datetime.datetime(2026, 11, 3, 15, 0)
    assert when(datetime.datetime(2026, 11, 3, 14, 32).timestamp(), now) == "Today 14:32"
    assert when(datetime.datetime(2026, 11, 2, 9, 5).timestamp(), now) == "Yesterday 09:05"
    assert "14" in when(datetime.datetime(2026, 10, 30, 14, 0).timestamp(), now)
    assert autosave.interval_text(1) == "Every minute" and autosave.interval_text(0) == "Off"


# -- Settings ----------------------------------------------------------------------------------

def settings_action(win: MainWindow) -> QAction:
    edit = next(a.menu() for a in win.menuBar().actions() if a.text() == "Edit")
    return next(a for a in edit.actions() if a.text() == "Settings...")


def combo(dialog: SettingsDialog, title: str) -> QComboBox:
    def titled(widget):
        name = dialog.form.labelForField(widget)
        return name is not None and name.text() == title

    return next(c for c in dialog.findChildren(QComboBox) if titled(c))


def test_settings_is_in_the_application_menu_on_macos(make_window):
    action = settings_action(make_window())
    assert action.menuRole() == QAction.MenuRole.PreferencesRole
    assert action.shortcut() == QKeySequence(QKeySequence.StandardKey.Preferences)


def test_settings_apply_at_once_and_are_kept(make_window, app):
    win = make_window()
    win.dialogs.show_settings()
    dialog = win.dialogs.settings_dialog
    assert dialog.isVisible()

    combo(dialog, "Appearance:").setCurrentIndex(2)   # Dark
    assert win.appearance == "dark" and win.menu_actions.choices["appearance"]["dark"].isChecked()
    win.set_appearance("light")                       # from the View menu
    assert combo(dialog, "Appearance:").currentText() == "Light"

    combo(dialog, "Autosave:").setCurrentText("Off")
    assert win.autosave.interval == 0 and not win.autosave.timer.isActive()

    before = win.current_canvas()
    combo(dialog, "Style:").setCurrentText("Wire")
    combo(dialog, "Labels:").setCurrentText("Element and Number")
    box = next(b for b in dialog.findChildren(QCheckBox) if b.text() == "Show hydrogens")
    box.setChecked(False)
    assert before.display.mode == Representation.BALL   # the open document keeps its style
    after = win.new_document()
    assert after.display.mode == Representation.WIRE
    assert after.display.labels == LabelMode.BOTH and not after.display.show_hydrogens

    dialog.close()
    again = make_window()   # remembered
    assert again.autosave.interval == 0 and again.appearance == "light"
    assert again.current_canvas().display.mode == Representation.WIRE
    again.dialogs.show_settings()
    dialog = again.dialogs.settings_dialog
    assert combo(dialog, "Style:").currentText() == "Wire"
    widths = {(c.minimumWidth(), c.maximumWidth()) for c in dialog.findChildren(QComboBox)}
    assert widths == {(FIELD_WIDTH, FIELD_WIDTH)}   # one column of equal fields


SETTINGS = '''
"""Has settings; used by the settings tests."""
seen = []

def register(context):
    context.setting = context.add_setting("cutoff", "Cutoff", 3.5, minimum=1.0, maximum=5.0, step=0.5,
                                          unit="angstrom", changed=seen.append)
    context.add_setting("enabled", "Enabled", True)
'''


def test_plugin_settings_are_in_settings(app, tmp_path, monkeypatch, folder):
    from miles.ui.plugins_dialog import PluginsDialog

    plugin_folder = tmp_path / "plugins"
    plugin_folder.mkdir()
    (plugin_folder / "cutoffs.py").write_text(SETTINGS)
    monkeypatch.setenv(plugins.FOLDER_ENV, str(plugin_folder))
    monkeypatch.delenv(plugins.NO_PLUGINS_ENV, raising=False)
    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    try:
        listing = PluginsDialog(win.plugins, win.theme, win)
        keys = [info.key for info in win.plugins.plugins]
        listing.table.selectRow(keys.index("file:cutoffs.py"))
        assert listing.settings_button.isEnabled()
        asked = []
        listing.settingsRequested.connect(asked.append)
        listing.open_settings()
        assert asked == ["file:cutoffs.py"] and listing.result() == QDialog.DialogCode.Accepted
        listing.table.selectRow(keys.index("miles.builtin_plugins.drawing"))
        assert not listing.settings_button.isEnabled()

        win.dialogs.show_settings("file:cutoffs.py")
        dialog = win.dialogs.settings_dialog
        assert "file:cutoffs.py" in dialog.sections
        assert "miles.builtin_plugins.hydrogen_bonds" in dialog.sections   # built-ins too
        field = next(f for f in dialog.findChildren(QDoubleSpinBox)
                     if f.suffix() == " angstrom" and dialog.form.labelForField(f).text() == "Cutoff:")
        field.setValue(4.0)
        seen = sys.modules["miles_user_plugins.cutoffs"].seen
        assert seen == [4.0]
        win.plugins.set_enabled("file:cutoffs.py", False)
        win.dialogs.show_settings()
        assert "file:cutoffs.py" not in win.dialogs.settings_dialog.sections
    finally:
        for canvas in win.canvases:
            canvas.doc.modified = False
        win.close()
        win.deleteLater()
        app.processEvents()
