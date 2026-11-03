"""The selection language: reading expressions, working them out, and where the window uses them."""

import numpy as np
import pytest
from PySide6.QtCore import QEvent, QSettings, Qt
from PySide6.QtGui import QKeyEvent
from PySide6.QtWidgets import QDialog

from miles.builtin_plugins.select_atoms import dialog as select_dialog
from miles.builtin_plugins.select_atoms.dialog import SelectionDialog
from miles.core.chem import from_smiles
from miles.core.document import Document
from miles.core.jobinfo import JobInfo
from miles.core.molecule import Molecule
from miles.core.neighbors import near_mask
from miles.core.selection import (
    KEYWORDS,
    Context,
    SelectionError,
    check_name,
    evaluate,
    find,
    parse,
    residue_ids,
)
from miles.files import session
from miles.ui.command_palette import menu_entries
from miles.ui.mainwindow import MainWindow

# symbol, name, resname, resseq, chain, hetero, x
ATOMS = [
    ("N", "N", "ALA", 1, "A", False, 0.0),     # 0
    ("H", "H", "ALA", 1, "A", False, -1.0),    # 1
    ("C", "CA", "ALA", 1, "A", False, 1.5),    # 2
    ("C", "CB", "ALA", 1, "A", False, 1.5),    # 3 (moved off the axis below)
    ("C", "C", "ALA", 1, "A", False, 3.0),     # 4
    ("O", "O", "ALA", 1, "A", False, 3.0),     # 5 (moved off the axis below)
    ("N", "N", "HIS", 2, "A", False, 4.5),     # 6
    ("C", "CA", "HIS", 2, "A", False, 6.0),    # 7
    ("C", "CB", "HIS", 2, "A", False, 6.0),    # 8 (moved off the axis below)
    ("N", "ND1", "HIS", 2, "A", False, 7.0),   # 9
    ("C", "C", "HIS", 2, "A", False, 7.5),     # 10
    ("O", "O", "HIS", 2, "A", False, 8.5),     # 11
    ("C", "C1", "LIG", 100, "A", True, 10.0),  # 12
    ("O", "O1", "LIG", 100, "A", True, 11.2),  # 13
    ("N", "N1", "LIG", 100, "A", True, 12.4),  # 14
    ("O", "O", "HOH", 200, "A", True, 13.0),   # 15
    ("H", "H1", "HOH", 200, "A", True, 13.9),  # 16
    ("H", "H2", "HOH", 200, "A", True, 12.7),  # 17 (moved off the axis below)
    ("Zn", "ZN", "ZN", 300, "A", True, 30.0),  # 18
    ("C", "CA", "GLY", 1, "B", False, 50.0),   # 19
    ("C", "C", "", 0, "", False, 70.0),        # 20 a molecule without residue data
    ("O", "O", "", 0, "", False, 71.2),        # 21
]
BONDS = [(0, 1), (0, 2), (2, 3), (2, 4), (4, 5), (4, 6), (6, 7), (7, 8), (8, 9), (7, 10), (10, 11),
         (12, 13), (12, 14), (15, 16), (15, 17), (20, 21)]


def structure() -> Molecule:
    coords = np.zeros((len(ATOMS), 3))
    coords[:, 0] = [a[6] for a in ATOMS]
    coords[3, 1] = coords[5, 1] = coords[8, 1] = coords[17, 1] = 1.0
    mol = Molecule.from_columns(coords, symbol=[a[0] for a in ATOMS], name=[a[1] for a in ATOMS],
                                resname=[a[2] for a in ATOMS], resseq=[a[3] for a in ATOMS],
                                chain=[a[4] for a in ATOMS], hetero=[a[5] for a in ATOMS],
                                ss=["H"] * 6 + [""] * 16)
    for i, j in BONDS:
        mol.add_bond(i, j)
    return mol


@pytest.fixture
def mol():
    return structure()


def pick(mol, text, context=None):
    return find(mol, text, context).tolist()


# -- groups -------------------------------------------------------------------------------------

@pytest.mark.parametrize(("text", "expected"), [
    ("all", list(range(22))),
    ("none", []),
    ("protein", [*range(12), 19]),
    ("nucleic", []),
    ("water", [15, 16, 17]),
    ("ion", [18]),
    ("metal", [18]),
    ("ligand", [12, 13, 14, 20, 21]),
    ("hetero", list(range(12, 19))),
    ("hydrogen", [1, 16, 17]),
    ("backbone", [0, 2, 4, 5, 6, 7, 10, 11, 19]),
    ("sidechain", [3, 8, 9]),
    ("polar", [0, 1, 5, 6, 9, 11, 13, 14, 15, 16, 17, 21]),
    ("helix", list(range(6))),
    ("loop", [6, 7, 8, 9, 10, 11, 19]),
    ("sheet", []),
])
def test_groups(mol, text, expected):
    assert pick(mol, text) == expected


def test_apolar_is_the_rest_of_polar(mol):
    assert set(pick(mol, "apolar")) == set(range(22)) - set(pick(mol, "polar"))


def test_aromatic_atoms_come_from_aromatic_bonds():
    benzene, _ = from_smiles("c1ccccc1O")
    assert pick(benzene, "aromatic") == [i for i, s in enumerate(benzene.symbols) if s == "C"]


# -- values -------------------------------------------------------------------------------------

@pytest.mark.parametrize(("text", "expected"), [
    ("name CA", [2, 7, 19]),
    ("name ca", [2, 7, 19]),                      # capitals of what was typed match too
    ("name CA CB", [2, 3, 7, 8, 19]),
    ("name C?", [2, 3, 7, 8, 12, 19]),
    ("name c*", [2, 3, 4, 7, 8, 10, 12, 19, 20]),
    ('name "C*"', []),                            # quoted: no wildcards
    ("resname his", list(range(6, 12))),
    ("resn HIS LIG", list(range(6, 15))),
    ("resid 2", list(range(6, 12))),
    ("resid 1-2", [*range(12), 19]),
    ("resi 100 to 300", list(range(12, 19))),
    ("resid 300 100", [12, 13, 14, 18]),
    ("chain B", [19]),
    ("chain b", []),                              # chain identifiers are matched exactly
    ("index 1-3 22", [0, 1, 2, 21]),
    ("element N O", [0, 5, 6, 9, 11, 13, 14, 15, 21]),
    ("elem zn", [18]),
    ("element 30", [18]),
    ('element "Zn"', [18]),
])
def test_values(mol, text, expected):
    assert pick(mol, text) == expected


@pytest.mark.parametrize(("text", "expected"), [
    ("x > 60", [20, 21]),
    ("x >= 70 and x<=70", [20]),
    ("y == 1", [3, 5, 8, 17]),
    ("y != 0 and protein", [3, 5, 8]),
    ("resid > 99", list(range(12, 19))),
    ("index <= 2", [0, 1]),
    ("mass > 60", [18]),
])
def test_comparisons(mol, text, expected):
    assert pick(mol, text) == expected


def test_comparisons_of_property_columns_and_context_columns(mol):
    mol.set_property("b_factor", np.arange(22, dtype=float) * 10)
    assert pick(mol, "b_factor >= 200") == [20, 21]
    assert pick(mol, "beta < 10") == [0]
    mol.set_property("depth", np.where(np.arange(22) < 3, 1.0, np.nan))
    assert pick(mol, "depth > 0") == [0, 1, 2]
    charges = np.zeros(22)
    charges[[5, 11]] = -0.5
    assert pick(mol, "charge < -0.3", Context(columns={"charge": lambda: charges})) == [5, 11]
    with pytest.raises(SelectionError, match="no charge values"):
        pick(mol, "charge < 0")
    with pytest.raises(SelectionError, match="no charge values"):
        pick(mol, "charge < 0", Context(columns={"charge": lambda: None}))
    with pytest.raises(SelectionError, match="Did you mean 'depth'"):
        pick(mol, "dept > 1")
    mol.set_property("label", ["x"] * 22, dtype=object, default="")
    with pytest.raises(SelectionError, match="holds text"):
        pick(mol, "label > 1")


# -- distances and structure -------------------------------------------------------------------

def test_within_and_around(mol):
    ligand = [12, 13, 14]
    assert pick(mol, "within 1.2 of resname LIG") == [*ligand, 15, 17]
    assert pick(mol, "around 1.2 of resname LIG") == [15, 17]
    assert pick(mol, "within 2.0 of resname LIG") == [11, *ligand, 15, 16, 17]
    assert pick(mol, "water within 1.2 of ligand") == [15, 17]      # "and" understood
    assert pick(mol, "within 1.2 of resname LIG and hydrogen") == [17]
    assert pick(mol, "within 0 of index 1") == [0]             # the atoms themselves
    assert pick(mol, "within 5 of none") == []


def test_within_counts_the_cutoff_itself(mol):
    assert 0 in pick(mol, "within 1.5 of name CA and resid 1")    # N is exactly 1.5 from CA


def test_expanding_to_residues_chains_and_neighbours(mol):
    assert pick(mol, "byres name CB and resname ALA") == list(range(6))
    assert pick(mol, "byres name ND1") == list(range(6, 12))
    assert pick(mol, "byres index 21") == [20, 21]            # no residue data: the whole molecule
    assert pick(mol, "bychain resid 2") == list(range(19))
    assert pick(mol, "bychain chain B") == [19]
    assert pick(mol, "bonded to name CA and resid 1") == [0, 3, 4]
    assert pick(mol, "bonded to (element O and water)") == [16, 17]
    assert pick(mol, "byres within 1.2 of resname LIG") == [12, 13, 14, 15, 16, 17]


def test_residue_ids_group_residues_and_loose_molecules(mol):
    ids = residue_ids(mol)
    assert len(set(ids[:6])) == 1 and ids[0] != ids[6]
    assert ids[20] == ids[21] and len(set(ids.tolist())) == 7   # ALA, HIS, LIG, HOH, ZN, GLY, loose


def test_smarts_and_the_view(mol):
    ethanol, _ = from_smiles("CCO")
    assert [ethanol.symbols[i] for i in pick(ethanol, 'smarts "[OX2H]"')] == ["O"]
    assert pick(ethanol, "smarts [#8]") == pick(ethanol, "element O")
    hidden = np.zeros(22, dtype=bool)
    hidden[15:18] = True
    view = Context(selected=[0, 2, 40], hidden=hidden)
    assert pick(mol, "selected", view) == [0, 2]
    assert pick(mol, "sele or water and visible", view) == [0, 2]
    assert pick(mol, "hidden", view) == [15, 16, 17]
    assert pick(mol, "visible", Context()) == list(range(22))
    assert pick(mol, "selected") == []


# -- combining ----------------------------------------------------------------------------------

@pytest.mark.parametrize(("text", "expected"), [
    ("not protein and not water and not ion", [12, 13, 14, 20, 21]),
    ("resname LIG or water and hydrogen", [12, 13, 14, 16, 17]),    # "and" before "or"
    ("(resname LIG or water) and hydrogen", [16, 17]),
    ("!protein & !hetero | chain B", [19, 20, 21]),
    ("NOT Protein AND Element O", [13, 15, 21]),                      # keywords in any case
    ("not not water", [15, 16, 17]),
    ("byres name CB and chain A", list(range(12))),                   # byres takes the next term only
    ("byres (name CB and resid 1)", list(range(6))),
])
def test_combinations_and_precedence(mol, text, expected):
    assert pick(mol, text) == expected


def test_named_selections(mol):
    names = {"pocket": "byres within 1.2 of resname LIG", "lig": "resname LIG", "wet": "pocket and water"}
    context = Context(names=names)
    assert pick(mol, "pocket and not lig", context) == [15, 16, 17]
    assert pick(mol, "WET", context) == [15, 16, 17]                # names in any case
    with pytest.raises(SelectionError, match="refers to itself"):
        pick(mol, "loop1", Context(names={"loop1": "loop2", "loop2": "loop1 or water"}))
    with pytest.raises(SelectionError, match="In the selection 'bad'"):
        pick(mol, "bad", Context(names={"bad": "protien"}))


def test_names_must_be_plain_words():
    check_name("pocket_2")
    for bad in ("", "2nd", "my pocket", "within", "Protein"):
        with pytest.raises(ValueError):
            check_name(bad)
    assert {"and", "or", "not", "within", "byres", "protein", "name"} <= KEYWORDS


# -- errors -------------------------------------------------------------------------------------

@pytest.mark.parametrize(("text", "message", "position"), [
    ("", "Type a selection", 0),
    ("protien", "Did you mean 'protein'", 0),
    ("within of ligand", "Expected a distance", 7),
    ("within 5 ligand", "Expected 'of'", 9),
    ("within -1 of water", "cannot be negative", 7),
    ("(protein", "Missing ')'", 8),
    ("protein)", "Unexpected ')'", 7),
    ("resname", "Expected a value after 'resname'", 7),
    ("resid abc", "not a residue number", 6),
    ("element Xx", "Unknown element 'Xx'", 8),
    ("protein and", "Expected a selection after 'and'", 11),
    ("chain A protein", "Expected 'and' or 'or' before 'protein'", 8),
    ("name 'O5", "Missing closing '", 5),
    ('smarts "[C"', "Not a valid SMARTS", 7),
    ("x >", "Expected a number after 'x >'", 3),
    ("and water", "Expected a selection before 'and'", 0),
    ("bonded water", "Expected 'to'", 7),
    ("name CA >", "Unexpected '>'", 8),
    ("water and chian A", "Did you mean 'chain'", 10),
])
def test_errors_say_what_and_where(mol, text, message, position):
    with pytest.raises(SelectionError) as caught:
        evaluate(mol, text)
    assert message in caught.value.message
    assert caught.value.position == position
    assert caught.value.pointer().splitlines()[-1] == " " * position + "^"


def test_parsed_trees_are_reused():
    assert parse("protein and name CA") is parse("protein and name CA")


@pytest.mark.parametrize("centres", [5, 400])
def test_grid_search_agrees_with_all_distances(centres):
    rng = np.random.default_rng(centres)
    points, centers = rng.random((300, 3)) * 20, rng.random((centres, 3)) * 20
    distance = np.linalg.norm(points[:, None] - centers[None], axis=2)
    assert np.array_equal(near_mask(points, centers, 2.5), (distance <= 2.5).any(axis=1))


# -- documents, the API, the console and the window ----------------------------------------------

def test_documents_keep_named_selections_in_sessions(mol):
    doc = Document(mol)
    doc.define("pocket", "byres within 1.2 of resname LIG")
    assert doc.modified
    assert doc.find("pocket and hydrogen").tolist() == [16, 17]
    with pytest.raises(SelectionError):
        doc.define("broken", "within of")
    with pytest.raises(ValueError):
        doc.define("and", "water")
    reopened = Document(*session.from_dict(session.to_dict(doc.mol, doc.session_info())))
    assert reopened.selections == {"pocket": "byres within 1.2 of resname LIG"}
    reopened.define("pocket", None)
    assert reopened.selections == {}
    doc.set_representation("hidden", [12])
    doc.set_hidden(np.arange(22) == 13)
    assert doc.find("hidden").tolist() == [12, 13]


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


def test_document_view_selects_with_expressions(window):
    window.open_structure(structure(), JobInfo(title="complex"))
    doc = window.api.current
    assert doc.select("water") == [15, 16, 17]
    assert doc.select("ion", add=True) == [15, 16, 17, 18]
    assert doc.select(atoms=[0]) == [0]
    assert doc.select([1], add=True, atoms=[2]) == [0, 1, 2]
    assert doc.find("selected and element N") == [0]
    assert doc.find("selected", selected=[5]) == [5]
    assert doc.count("protein") == 13
    doc.define("pocket", "byres within 1.2 of ligand")
    assert doc.selections == {"pocket": "byres within 1.2 of ligand"}
    assert doc.select("pocket and not ligand") == [15, 16, 17]
    with pytest.raises(ValueError, match="Unknown word"):
        doc.select("pockt")


def test_select_dialog_counts_previews_and_selects(window, monkeypatch):
    window.open_structure(structure(), JobInfo(title="complex"))
    canvas = window.current_canvas()
    canvas.set_selection([0])
    select_dialog.history.clear()
    seen = {}

    def run(dialog):
        dialog.set_expression("protein and name CA")
        seen["message"] = dialog.message.text()
        seen["preview"] = list(canvas.selection)
        dialog.set_expression("protien")
        seen["error"] = dialog.message.text()
        seen["enabled"] = dialog.select.isEnabled()
        seen["restored"] = list(canvas.selection)
        dialog.set_expression("water")
        dialog.add.setChecked(True)
        dialog.accept()
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(SelectionDialog, "exec", run)
    window.api.run("Edit > Select...")
    assert seen["message"] == "3 atoms in 3 residues"
    assert seen["preview"] == [2, 7, 19]
    assert "Did you mean 'protein'" in seen["error"] and not seen["enabled"]
    assert seen["restored"] == [0]
    assert sorted(canvas.selection) == [0, 15, 16, 17]
    assert select_dialog.history == ["water"]
    actions = {e.path: e.action for e in menu_entries(window.menuBar())}
    assert actions["Edit > Select"].shortcut().toString() == "Ctrl+F"


def test_select_dialog_cancel_examples_and_history(window, app):
    window.open_structure(structure(), JobInfo(title="complex"))
    canvas = window.current_canvas()
    canvas.set_selection([21])
    select_dialog.history[:] = ["ion", "water"]
    doc = window.api.current
    doc.define("pocket", "byres within 1.2 of ligand")
    dialog = SelectionDialog(window.api, doc)
    assert "pocket" in dialog.examples.itemText(1)
    dialog._example_chosen(1)
    assert dialog.field.text() == "pocket" and dialog.atoms == [12, 13, 14, 15, 16, 17, 20, 21]
    up = QKeyEvent(QEvent.Type.KeyPress, Qt.Key.Key_Up, Qt.KeyboardModifier.NoModifier)
    down = QKeyEvent(QEvent.Type.KeyPress, Qt.Key.Key_Down, Qt.KeyboardModifier.NoModifier)
    dialog.eventFilter(dialog.field, up)
    assert dialog.field.text() == "ion"
    dialog.eventFilter(dialog.field, up)
    assert dialog.field.text() == "water"
    dialog.eventFilter(dialog.field, down)
    dialog.eventFilter(dialog.field, down)
    assert dialog.field.text() == "pocket"     # back to what was there
    dialog._example_chosen(dialog.examples.findData('smarts "[CX3]=[OX1]"'))
    assert dialog.field.text().startswith("smarts")
    dialog.set_expression("ligand")
    assert canvas.selection == [12, 13, 14, 20, 21]
    dialog.reject()
    assert canvas.selection == [21]
