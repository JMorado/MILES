"""The main window's parts: the tab model, menu items that name what they call, and closing tabs."""

import gc
import weakref

import pytest
from PySide6.QtCore import QEvent, QSettings
from PySide6.QtWidgets import QApplication

from miles.core.chem import from_smiles
from miles.ui.canvas import Canvas
from miles.ui.mainwindow import MainWindow
from miles.ui.menus import MENUS, ChoiceItem, Command, Dynamic, Item, Submenu, ToggleItem, resolve


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


def flush(app: QApplication) -> None:
    """Let deleted widgets go, as the event loop would."""
    for _ in range(3):
        app.processEvents()
        app.sendPostedEvents(None, QEvent.Type.DeferredDelete.value)
    gc.collect()


def test_the_current_document_is_announced_once_per_change(window):
    seen: list[Canvas] = []
    window.currentDocumentChanged.connect(seen.append)
    first = window.current_canvas()
    second = window.new_document()
    assert seen == [second]
    window.tabs.set_current(first)
    window.tabs.set_current(first)   # already in front: nothing to announce
    assert seen == [second, first]
    seen.clear()
    window.close_canvas(first)
    assert seen == [second] and window.canvases == [second]


def test_tab_order_follows_the_tab_bar(window):
    first = window.current_canvas()
    second = window.new_document()
    third = window.new_document()
    assert window.canvases == [first, second, third]
    window.tabbar.moveTab(2, 0)   # dragging a tab
    assert window.canvases == [third, first, second]
    assert window.current_canvas() is third and window.stage.stack.currentWidget() is third
    window.close_tab(1)
    assert window.canvases == [third, second]


def test_a_closed_tab_lets_its_canvas_go(app, window):
    window.api.new(from_smiles("CCO")[0], "ethanol")
    closed = weakref.ref(window.current_canvas())
    window.current_canvas().doc.modified = False
    window.close_current_tab()
    flush(app)
    assert closed() is None   # with its scene arrays: closing big tabs gives the memory back


def test_closing_the_last_tab_opens_an_empty_one(window):
    only = window.current_canvas()
    window.close_canvas(only)
    assert len(window.canvases) == 1 and window.current_canvas() is not only


def _items(items: list[Item]) -> list[Item]:
    found: list[Item] = []
    for item in items:
        if isinstance(item, Submenu):
            found.extend(_items(item.items))
        else:
            found.append(item)
    return found


def test_every_menu_item_names_something_that_exists(window):
    canvas = window.current_canvas()
    for item in _items([entry for menu in MENUS for entry in menu.items]):
        if not isinstance(item, Command | ToggleItem | ChoiceItem | Dynamic):
            continue
        head, _, rest = item.slot.partition(".")
        owner, path = (canvas, rest) if head in ("canvas", "structure") else (window, item.slot)
        target = owner
        for name in path.split("."):
            assert hasattr(target, name), f"{item.slot} names nothing"
            target = getattr(target, name)
        assert callable(target), item.slot
        assert callable(resolve(window, item.slot))


def test_canvas_commands_act_on_the_tab_in_front(window):
    first = window.current_canvas()
    window.api.new(from_smiles("CCO")[0], "ethanol")
    second = window.current_canvas()
    window.api.run("Edit > Select All")
    assert second.selection == list(range(len(second.mol))) and first.selection == []
    window.tabs.set_current(first)
    before = first.toast_text
    window.api.run("View > Show Only Selection")   # an empty structure: nothing happens, no message
    assert first.toast_text == before


def _heavy_atoms_only(smiles: str):
    mol = from_smiles(smiles)[0]
    mol.remove_atoms([i for i, a in enumerate(mol.atoms) if a.symbol == "H"])
    return mol


def test_large_structures_get_hydrogens_on_a_worker(window, monkeypatch):
    from miles.ui.canvas import commands

    monkeypatch.setattr(commands, "BACKGROUND_ATOMS", 1)   # every structure counts as large
    doc = window.api.new(_heavy_atoms_only("CCO"), "ethanol")
    canvas = window.current_canvas()
    heavy = len(canvas.mol)
    canvas.commands.add_hydrogens()
    assert len(canvas.mol) == heavy and window.stage.tasks.tasks   # still running: the window stays usable
    window.api.wait()
    assert len(canvas.mol) == 9 and "Added 6 hydrogens" in canvas.toast_text
    doc.undo()   # one undoable step
    assert len(canvas.mol) == heavy


def test_a_background_result_is_dropped_when_the_structure_changed_meanwhile(window, monkeypatch):
    from miles.ui.canvas import commands

    monkeypatch.setattr(commands, "BACKGROUND_ATOMS", 1)
    window.api.new(_heavy_atoms_only("CCO"), "ethanol")
    canvas = window.current_canvas()
    canvas.commands.add_hydrogens()
    canvas.doc.edit(lambda mol: mol.remove_atoms([2]))   # an edit while the hydrogens are worked out
    window.api.wait()
    assert len(canvas.mol) == 2 and canvas.toast_text == commands.CHANGED_MEANWHILE


def test_large_structures_are_cleaned_on_a_worker(window, monkeypatch):
    from miles.ui.canvas import commands

    monkeypatch.setattr(commands, "BACKGROUND_ATOMS", 1)
    doc = window.api.new(from_smiles("CCO")[0], "ethanol")
    canvas = window.current_canvas()
    before = canvas.mol.coords.copy()
    canvas.doc.edit(lambda mol: setattr(mol, "coords", mol.coords * 1.3))   # stretched
    canvas.commands.clean_geometry()
    window.api.wait()
    assert canvas.toast_text.startswith(("Cleaned", "Tidied"))
    coords = canvas.mol.coords
    lengths = [float(((coords[i] - coords[j]) ** 2).sum() ** 0.5) for i, j in canvas.mol.bonds]
    assert all(0.9 < d < 1.6 for d in lengths)   # back to sensible bond lengths
    doc.undo()
    assert abs(canvas.mol.coords - before * 1.3).max() < 1e-9
