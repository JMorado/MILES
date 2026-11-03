"""The Python console: running code, its namespace, history, completion and scripts."""

import sys

import pytest
from PySide6.QtCore import QSettings, Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QFileDialog

from miles.core.chem import from_smiles
from miles.files.common import JobInfo
from miles.ui.app import parse_args
from miles.ui.mainwindow import MainWindow


def make_window(app, path):
    win = MainWindow(QSettings(str(path), QSettings.Format.IniFormat))
    win.resize(1000, 700)
    win.show()
    app.processEvents()
    return win


def close(win):
    for canvas in win.canvases:
        canvas.doc.modified = False
    win.close()


@pytest.fixture
def window(app, tmp_path):
    win = make_window(app, tmp_path / "settings.ini")
    win.show_console(True)
    yield win
    close(win)


def run(win, *lines):
    console = win.stage.console
    for line in lines:
        console.execute(line)
    return console.text()


def test_statements_expressions_and_blocks(window):
    text = run(window, "x = 6 * 7", "x", "def double(a):", "    return 2 * a", "", "double(x)")
    assert ">>> x\n42\n" in text and "... " in text and "\n84\n" in text
    assert window.stage.console.prompt.text() == ">>>"


def test_output_and_errors_stay_in_the_console(window, capsys):
    stdout = sys.stdout
    text = run(window, "print('hello from the console')", "1 / 0", "undefined_name", "x = (")
    assert "hello from the console" in text
    assert "ZeroDivisionError" in text and "NameError" in text
    assert sys.stdout is stdout   # redirected only while a command runs
    assert capsys.readouterr().out == ""
    assert window.stage.console.prompt.text() == "..."   # the open parenthesis waits for more
    run(window, ")")
    assert window.stage.console.prompt.text() == ">>>"


def test_namespace_follows_the_current_tab(window):
    ethanol, _ = from_smiles("CCO")
    water, _ = from_smiles("O")
    window.open_structure(ethanol, JobInfo(title="ethanol"))
    window.open_structure(water, JobInfo(title="water"))
    assert run(window, "print(doc.name, len(mol))").rstrip().endswith("water 3")
    window.tabbar.setCurrentIndex(0)
    assert run(window, "print(doc.name, len(mol))").rstrip().endswith("ethanol 9")
    run(window, "mol = 'mine'")   # bound by the user: left alone
    window.tabbar.setCurrentIndex(1)
    assert run(window, "print(mol, doc.name)").rstrip().endswith("mine water")


def test_helpers_and_direct_changes_redraw(window):
    mol, _ = from_smiles("CCO")
    window.open_structure(mol, JobInfo(title="ethanol"))
    canvas = window.current_canvas()
    run(window, "select('smarts \"[OX2H]\"')", "show('tube')")
    assert [canvas.mol.atoms[i].symbol for i in canvas.selection] == ["O"]
    assert run(window, "print(count('element C'))").rstrip().endswith("2")
    run(window, "define('carbons', 'element C')", "select('carbons or hydrogen and bonded to element O')")
    assert sorted(canvas.mol.atoms[i].symbol for i in canvas.selection) == ["C", "C", "H"]
    assert "Unknown word" in run(window, "select('protien')")
    assert str(canvas.display.mode) == "tube"
    assert "1.5" in run(window, "print(round(dist(0, 1), 1))")
    revision = canvas.doc.revision
    run(window, "mol.coords = mol.coords + 0.5")   # changed without an edit: the view still follows
    assert canvas.doc.revision > revision
    assert "Unknown style" in run(window, "show('sticks')")


def test_exit_help_and_input_are_safe(window):
    text = run(window, "exit()", "help()", "help(select)", "input()")
    assert "stays open" in text and "Helpers:" in text and "expression" in text and "RuntimeError" in text
    assert window.isVisible()


def test_completion(window):
    console = window.stage.console
    console.input.setText("app.curr")
    console.complete()
    assert console.input.text() == "app.current"
    console.input.setText("dihe")
    QTest.keyClick(console.input, Qt.Key.Key_Tab)
    assert console.input.text().startswith("dihedral")


def test_history_is_kept_between_sessions(app, tmp_path):
    first = make_window(app, tmp_path / "settings.ini")
    first.show_console(True)
    run(first, "a = 1", "b = 2", "b = 2")   # a repeat is kept once
    console = first.stage.console
    console.set_console_height(300)
    console.save_height()
    close(first)

    second = make_window(app, tmp_path / "settings.ini")
    console = second.stage.console
    assert console.isVisible()   # it was open when the first window closed
    assert console.input.history[-2:] == ["a = 1", "b = 2"]
    assert console.wanted_height == 300
    QTest.keyClick(console.input, Qt.Key.Key_Up)
    assert console.input.text() == "b = 2"
    QTest.keyClick(console.input, Qt.Key.Key_Up)
    assert console.input.text() == "a = 1"
    QTest.keyClick(console.input, Qt.Key.Key_Down)
    QTest.keyClick(console.input, Qt.Key.Key_Down)
    assert console.input.text() == ""
    close(second)


def test_menu_toggle_close_link_and_clear(window):
    toggle = window.menu_actions.toggles["python_console"]
    console = window.stage.console
    assert toggle.isChecked() and console.isVisible()
    run(window, "print('x')")
    console.clear()
    assert console.text() == ""
    console.close_console()
    assert not console.isVisible() and not toggle.isChecked()
    toggle.trigger()
    assert console.isVisible()


def test_scripts_from_the_menu_and_the_command_line(window, tmp_path, monkeypatch):
    script = tmp_path / "make.py"
    script.write_text("from miles.core.chem import from_smiles\n"
                      "made, _ = from_smiles('CCN')\n"
                      "app.new(made, 'ethylamine')\n"
                      "print('script ran as', __name__)\n", encoding="utf-8")
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *args, **kwargs: (str(script), ""))
    window.run_script_dialog()
    assert window.current_canvas().doc.name == "ethylamine"
    text = window.stage.console.text()
    assert "script ran as __main__" in text
    assert "made" in window.stage.console.names   # what a script defines stays available
    assert window.run_code("print(len(made))")
    assert window.stage.console.text().rstrip().endswith("10")   # C2H7N

    reported = []
    monkeypatch.setattr(sys, "excepthook", lambda *exc: reported.append(exc[0]))
    broken = tmp_path / "broken.py"
    broken.write_text("raise ValueError('bad script')\n", encoding="utf-8")
    assert not window.run_script(str(broken))
    assert "ValueError: bad script" in window.stage.console.text()
    assert reported == [ValueError]   # also handed to the error dialog


def test_command_line_options():
    args = parse_args(["--script", "setup.py", "-c", "print(1)", "water.xyz"])
    assert args.script == "setup.py" and args.code == "print(1)" and args.items == ["water.xyz"]


def test_console_docks_along_the_bottom(window):
    stage, console = window.stage, window.stage.console
    stage.layout_overlays()
    assert console.width() > stage.width() * 0.8
    assert console.geometry().bottom() < stage.height() and console.geometry().bottom() > stage.height() - 60
    before = console.height()
    console.set_console_height(before + 80)
    stage.layout_overlays()
    assert console.height() == before + 80


def test_running_code_does_not_open_the_console_at_the_next_start(app, tmp_path):
    settings = tmp_path / "settings.ini"
    QSettings(str(settings), QSettings.Format.IniFormat).setValue("console/open", "false")
    win = make_window(app, settings)
    assert win.run_code("x = 1")
    assert not win.stage.console.isHidden()
    assert win.run_script(str(tmp_path / "missing.py")) is False
    assert str(win.settings.value("console/open")) == "false"
    close(win)

    again = make_window(app, settings)
    assert again.stage.console.isHidden()
    toggle = again.menu_actions.toggles["python_console"]
    toggle.trigger()   # View > Python Console: the user's choice is kept
    assert str(again.settings.value("console/open")) == "true"
    toggle.trigger()
    assert str(again.settings.value("console/open")) == "false"
    again.show_console(True)
    again.stage.console.close_console()   # the Close link
    assert str(again.settings.value("console/open")) == "false"
    close(again)


def test_console_open_at_start_is_docked_along_the_bottom(app, tmp_path):
    # the layout runs before the window is shown, when no child counts as visible yet;
    # the console must still get its place without waiting for a resize
    settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    settings.setValue("console/open", "true")
    settings.sync()
    win = make_window(app, tmp_path / "settings.ini")
    try:
        console, stage = win.stage.console, win.stage
        assert console.isVisible()
        assert console.width() > stage.width() * 0.9                          # across the window
        assert console.geometry().bottom() > stage.height() - 60               # at the bottom edge
        assert console.y() > stage.height() / 2
    finally:
        close(win)
