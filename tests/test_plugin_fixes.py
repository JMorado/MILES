"""Regressions of the plugin interface and the built-in plugins, each from the probe that found it."""

import importlib
import sys
import textwrap
import types

import numpy as np
import pytest
from PySide6.QtCore import QEvent, QSettings
from PySide6.QtGui import QAction, QActionGroup
from PySide6.QtWidgets import QApplication, QDialog

from miles import parse_text, water_mask
from miles.builtin_plugins.binding_pocket.interactions import (
    InteractionKind,
    find_interactions,
    ligand_atoms,
    pocket_atoms,
)
from miles.builtin_plugins.prepare_protein.preparation import add_hydrogens_at_ph, prepare, water_count
from miles.builtin_plugins.shared.hbonds import HydrogenBondFinder
from miles.plugin import PluginDialog, SceneItems, Severity
from miles.plugin.render import SurfaceColoring, SurfaceRequest, surface_mesh
from miles.plugins import FOLDER_ENV, NO_PLUGINS_ENV, PluginState
from miles.testing import (
    captured_messages,
    document_from_smiles,
    document_from_text,
    drawings,
    load_plugin,
    plugin_setting,
    reload_plugin,
    switch_off,
    switch_on,
)
from miles.ui.mainwindow import MainWindow

SURFACE = "miles.builtin_plugins.molecular_surface"


def pdb_line(serial, name, resname, resseq, xyz, element, hetero=False, icode=" "):
    record = "HETATM" if hetero else "ATOM"
    label = f" {name:<3}" if len(name) < 4 else name
    x, y, z = xyz
    return (f"{record:<6}{serial:>5} {label} {resname:>3} A{resseq:>4}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}")


def plugin_module(name, source):
    module = types.ModuleType(name)
    exec(textwrap.dedent(source), module.__dict__)
    return module


def builtin(app, key):
    return next(info for info in app.window.plugins.plugins if info.key == key)


def deleted_later():
    QApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete.value)
    QApplication.processEvents()


# -- the plugin interface ----------------------------------------------------------------------

def test_drawings_that_follow_the_selection_are_drawn_again_when_it_changes(miles_app):
    module = plugin_module("follows_selection", """
        from miles.plugin import Annotations, Label

        def register(context):
            def labels(doc, drawing):
                return Annotations(labels=[Label(str(i + 1), i) for i in doc.selection])
            context.add_annotations(labels)
    """)
    load_plugin(miles_app, module)
    doc = document_from_smiles(miles_app, "CCO", "ethanol")
    canvas = doc.canvas
    count = lambda: sum(len(a.labels) for a in canvas.drawing.annotations())   # noqa: E731
    doc.select([0])
    assert count() == 1
    doc.select([0, 1, 2])   # no update() from the plugin: the selection change is enough
    assert count() == 3


def test_reload_keeps_shortcuts_and_notes_are_not_errors(miles_app):
    module = plugin_module("shortcuts", """
        def register(context):
            context.add_command("Tools > Count Things", lambda app: None, shortcut="Ctrl+Alt+Shift+J")
            context.add_command("Tools > Clean Again", lambda app: None, shortcut="Ctrl+L")
    """)
    info = load_plugin(miles_app, module)
    assert "already used" in info.warnings and not info.error and info.state == PluginState.LOADED
    for _ in range(2):
        info = reload_plugin(miles_app, info)
        actions = {a.text(): a for a in miles_app.widget.findChildren(QAction)}
        assert actions["Count Things"].shortcut().toString() == "Ctrl+Alt+Shift+J"
        assert info.warnings.count("already used") == 1 and not info.error


def test_reload_imports_a_package_plugin_and_its_own_modules_again(app, tmp_path, monkeypatch):
    folder = tmp_path / "plugins"
    package = folder / "valued"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("from .helper import VALUE\n\n"
                                         "def register(context):\n"
                                         "    context.info.description = str(VALUE)\n", encoding="utf-8")
    (package / "helper.py").write_text("VALUE = 1\n", encoding="utf-8")
    monkeypatch.setenv(FOLDER_ENV, str(folder))
    monkeypatch.delenv(NO_PLUGINS_ENV, raising=False)
    window = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    try:
        info = next(i for i in window.plugins.plugins if i.key == "file:valued")
        assert info.description == "1"
        # another size: no stale bytecode
        (package / "helper.py").write_text("VALUE = 22\n", encoding="utf-8")
        window.plugins.reload(info)
        assert info.description == "22" and info.state == PluginState.LOADED
    finally:
        window.close()


def test_reload_imports_a_module_plugin_again(miles_app, tmp_path, monkeypatch):
    (tmp_path / "counted").mkdir()
    (tmp_path / "counted" / "__init__.py").write_text("from .helper import VALUE\n\n"
                                                     "def register(context):\n"
                                                     "    context.info.description = str(VALUE)\n",
                                                     encoding="utf-8")
    (tmp_path / "counted" / "helper.py").write_text("VALUE = 1\n", encoding="utf-8")
    monkeypatch.syspath_prepend(str(tmp_path))
    info = load_plugin(miles_app, importlib.import_module("counted"))
    (tmp_path / "counted" / "helper.py").write_text("VALUE = 22\n", encoding="utf-8")
    importlib.invalidate_caches()
    assert reload_plugin(miles_app, info).description == "22"
    sys.modules.pop("counted", None)
    sys.modules.pop("counted.helper", None)


def test_choices_leave_no_action_groups_behind(miles_app):
    count = lambda: len(miles_app.widget.findChildren(QActionGroup))   # noqa: E731
    before = count()
    for key in ("miles.builtin_plugins.charges", SURFACE):
        for _ in range(3):
            reload_plugin(miles_app, builtin(miles_app, key))
    assert count() == before   # gone from the window at once, not only once deleted later
    miles_app.run("View > Charges > Gasteiger")
    miles_app.run("View > Surface > Hydrophobicity")


def test_dialogs_are_deleted_once_closed(miles_app, monkeypatch):
    monkeypatch.setattr(PluginDialog, "exec", lambda dialog: QDialog.DialogCode.Rejected)
    document_from_smiles(miles_app, "CCO", "ethanol")
    deleted_later()
    before = len(miles_app.widget.findChildren(QDialog))
    for _ in range(3):
        for command in ("Edit > Select...", "Build from Sequence...", "File > Open from PDB...",
                        "File > Export Movie..."):
            miles_app.run(command)
    assert len(miles_app.widget.findChildren(QDialog)) == before


def test_a_command_that_is_not_available_cannot_be_run(miles_app):
    module = plugin_module("unavailable", """
        events = []

        def register(context):
            context.add_command("Tools > Never Now", lambda app: events.append(1), enabled=lambda app: False)
    """)
    load_plugin(miles_app, module)
    with pytest.raises(RuntimeError, match="not available"):
        miles_app.run("Tools > Never Now")
    assert module.events == []


def test_colors_in_any_form_and_testing_helpers(miles_app):
    module = plugin_module("colors_and_settings", """
        from PySide6.QtGui import QColor
        from miles.plugin import Annotations, Label, Line

        def register(context):
            size = context.add_setting("size", "Size", 3.0, minimum=1.0, maximum=9.0)
            def drawn(doc, drawing):
                return Annotations(lines=[Line(0, 1, (255, 0, 0)), Line(1, 2, QColor("#00ff00")),
                                          Line(0, 2, (0.0, 0.0, 1.0))],
                                   labels=[Label(f"{size.value:g}", 0, color="#123456")])
            context.add_annotations(drawn)
            context.add_command("Tools > Complain", lambda app: app.message("Too big", severity="warning"))
            context.add_command("Tools > Fail", lambda app: app.warn("Fail", "It broke"))
    """)
    info = load_plugin(miles_app, module)
    document_from_text(miles_app, "3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n", "water")
    (annotations,) = drawings(miles_app, info)
    assert [line.color for line in annotations.lines] == ["#ff0000", "#00ff00", "#0000ff"]
    assert annotations.labels[0].text == "3" and annotations.labels[0].color == "#123456"
    plugin_setting(miles_app, info, "size").set(5.0)
    assert drawings(miles_app, info)[0].labels[0].text == "5"
    with captured_messages(miles_app) as messages:
        miles_app.run("Tools > Complain")
        miles_app.run("Tools > Fail")   # an error: recorded, no dialog waits for a click
    assert [(m.text, m.severity) for m in messages] == [("Too big", Severity.WARNING),
                                                        ("It broke", Severity.ERROR)]
    switch_off(miles_app, info)
    assert info.state == PluginState.OFF and "Tools > Complain" not in miles_app.commands()
    switch_on(miles_app, info)
    assert "Tools > Complain" in miles_app.commands()


# -- the built-in plugins ----------------------------------------------------------------------

def test_hydrogen_bonds_agree_between_the_layer_and_the_pocket():
    nitrogen, oxygen = np.array([0.0, 0.0, 0.0]), np.array([-1.0, -2.7, 0.0])
    hydrogen = oxygen + 0.96 * (nitrogen - oxygen) / np.linalg.norm(nitrogen - oxygen)
    mol = parse_text("\n".join([
        pdb_line(1, "N", "ALA", 1, (0, 0, 0), "N"), pdb_line(2, "H", "ALA", 1, (-0.5, 0.87, 0), "H"),
        pdb_line(3, "CA", "ALA", 1, (1.45, 0, 0), "C"), pdb_line(4, "C", "ALA", 1, (2.0, 1.4, 0), "C"),
        pdb_line(5, "O", "ALA", 1, (1.3, 2.4, 0), "O"), pdb_line(6, "CB", "ALA", 1, (2.0, -0.8, 1.2), "C"),
        pdb_line(7, "O1", "MOH", 50, tuple(oxygen), "O", True),
        pdb_line(8, "H1", "MOH", 50, tuple(hydrogen), "H", True),
        pdb_line(9, "C1", "MOH", 50, (-2.3, -3.3, 0), "C", True), "END"]) + "\n")
    ligand = ligand_atoms(mol, 6)
    pocket = pocket_atoms(mol, ligand)
    found = {frozenset((b.donor, b.acceptor)) for b in HydrogenBondFinder(mol).find(mol.coords)}
    crossing = {pair for pair in found if len(pair & set(ligand)) == 1}
    pocket_bonds = {frozenset((*f.first, *f.second)) for f in find_interactions(mol, ligand, pocket)
                    if f.kind == InteractionKind.HYDROGEN_BOND}
    assert pocket_bonds == crossing == set()   # an amide nitrogen accepts no hydrogen bond


def test_prepare_protein_counts_water_as_molecules():
    rows = [pdb_line(1, "N", "GLY", 1, (0, 0, 0), "N"), pdb_line(2, "CA", "GLY", 1, (1.45, 0, 0), "C"),
            pdb_line(3, "C", "GLY", 1, (2.0, 1.4, 0), "C"), pdb_line(4, "O", "GLY", 1, (1.3, 2.4, 0), "O")]
    for k in range(2):
        x = 10 + 5 * k
        rows += [pdb_line(5 + 3 * k, "O", "HOH", 100 + k, (x, 0, 0), "O", True),
                 pdb_line(6 + 3 * k, "H1", "HOH", 100 + k, (x + 0.96, 0, 0), "H", True),
                 pdb_line(7 + 3 * k, "H2", "HOH", 100 + k, (x - 0.24, 0.93, 0), "H", True)]
    mol = parse_text("\n".join(rows) + "\nEND\n")
    assert water_count(mol) == 2 and water_mask(mol).sum() == 6
    assert prepare(mol, complete=False, hydrogens=False, waters=True).waters_removed == 2


def test_added_hydrogens_keep_the_insertion_code_of_their_residue():
    mol = parse_text("\n".join([
        pdb_line(1, "N", "GLY", 52, (0, 0, 0), "N", icode="A"),
        pdb_line(2, "CA", "GLY", 52, (1.45, 0, 0), "C", icode="A"),
        pdb_line(3, "C", "GLY", 52, (2.0, 1.4, 0), "C", icode="A"),
        pdb_line(4, "O", "GLY", 52, (1.3, 2.4, 0), "O", icode="A")]) + "\nEND\n")
    add_hydrogens_at_ph(mol, 7.4)
    hydrogens = [a for a in mol.atoms if a.symbol == "H"]
    assert hydrogens and {(a.resseq, a.icode) for a in hydrogens} == {(52, "A")}


def test_grids_added_by_the_user_are_undoable_and_derived_ones_are_not_edits(miles_app):
    from miles import Volume

    doc = document_from_smiles(miles_app, "O", "water")
    grid = Volume("grid", np.zeros(3), np.eye(3) * 0.5, np.zeros((4, 4, 4)))
    canvas_doc = doc.canvas.doc
    steps, modified = len(canvas_doc._undo), canvas_doc.modified
    doc.keep_volume(grid)   # e.g. an orbital worked out for display: no undo step, no unsaved change
    assert len(doc.results.volumes) == 1
    assert len(canvas_doc._undo) == steps and canvas_doc.modified == modified
    doc.add_volume(grid)    # a grid the user adds: one undoable step
    assert len(doc.results.volumes) == 2 and canvas_doc.modified
    doc.undo()
    assert len(doc.results.volumes) == 1


def surface_items(app, doc, export=False):
    found = [d for d in drawings(app, builtin(app, SURFACE), doc, export) if isinstance(d, SceneItems)]
    return found[0].meshes[0].mesh if found else None


def test_the_surface_is_worked_out_off_the_interface_thread_and_once(miles_app):
    doc = document_from_smiles(miles_app, "CCO", "ethanol")
    miles_app.run("View > Surface > Show Surface")
    assert surface_items(miles_app, doc) is None and miles_app.window.stage.tasks.tasks   # started, not drawn
    miles_app.wait()
    mesh = surface_items(miles_app, doc)
    assert mesh is not None and not mesh.empty
    for _ in range(3):   # drawn again without a change: the same mesh, no new work
        assert surface_items(miles_app, doc) is mesh and not miles_app.window.stage.tasks.tasks


def test_exported_frames_get_the_surface_of_what_they_show(miles_app):
    doc = document_from_smiles(miles_app, "CCO", "ethanol")
    base = doc.molecule.coords.copy()
    moved = base.copy()
    moved[0] += 1.5
    doc.set_frames(np.array([base, moved]))
    doc.frame = 0
    miles_app.run("View > Surface > Show Surface")
    surface_items(miles_app, doc)
    miles_app.wait()
    first = surface_items(miles_app, doc)
    doc.frame = 1   # as the movie export does, with no time for the worker
    exported = surface_items(miles_app, doc, export=True)
    assert exported is not first
    assert not np.allclose(exported.vertices.mean(axis=0), first.vertices.mean(axis=0))
    expected = surface_mesh(doc.molecule, SurfaceRequest(None, SurfaceColoring.ELEMENT), miles_app.theme)
    assert np.allclose(exported.vertices, expected.vertices)


def test_the_surface_wraps_the_same_atoms_after_others_are_deleted(miles_app):
    doc = document_from_smiles(miles_app, "CCO", "ethanol")
    doc.select([1, 2])
    miles_app.run("View > Surface > Show Surface")
    doc.edit(lambda mol: mol.remove_atoms([0]))   # the wrapped atoms are now 0 and 1
    after = surface_items(miles_app, doc, export=True)
    expected = surface_mesh(doc.molecule, SurfaceRequest((0, 1), SurfaceColoring.ELEMENT), miles_app.theme)
    assert np.allclose(after.vertices, expected.vertices)
