"""The document's history: frames, playback, animations, results, selection and failed edits under undo."""

import gzip
import json

import numpy as np
import pytest
from PySide6.QtCore import QSettings

from miles.core.document import Change, Document
from miles.core.frames import Trajectory
from miles.core.jobinfo import JobInfo
from miles.core.molecule import Molecule
from miles.core.vibrations import VibrationalMode
from miles.core.volume import Volume
from miles.files import session
from miles.files.common import FileFormatError

H2_FRAMES = np.array([[[0, 0, 0], [0.74, 0, 0]],
                      [[0, 0, 0], [1.00, 0, 0]],
                      [[0, 0, 0], [2.00, 0, 0]]], dtype=float)


def hydrogen(coords=None) -> Molecule:
    return Molecule.from_columns(np.zeros((2, 3)) if coords is None else coords, symbol=["H", "H"])


def chain(n: int = 4) -> Molecule:
    """A row of carbons, bonded in order."""
    mol = Molecule.from_columns([[1.5 * k, 0.0, 0.0] for k in range(n)], symbol=["C"] * n)
    for k in range(n - 1):
        mol.add_bond(k, k + 1, 1)
    return mol


def recorder(doc: Document) -> list[Change]:
    events: list[Change] = []
    doc.listeners.append(lambda event: events.append(event.kind))
    return events


# -- frames and undo ----------------------------------------------------------------------------

def test_undo_after_moving_to_another_frame_restores_the_edited_frame_only():
    doc = Document(hydrogen())
    doc.set_frames(H2_FRAMES)
    doc.set_frame(1)
    doc.edit(lambda m: m.set_positions([1], [[1.5, 0, 0]]), kind=Change.COORDS)
    assert doc.frames[1][1][0] == 1.5          # the edit is kept in the frame shown
    doc.set_frame(2)
    doc.undo()
    assert doc.frames[1][1][0] == 1.0
    assert doc.frames[2][1][0] == 2.0           # the frame shown at undo time is untouched
    assert doc.frame == 1 and doc.mol.coords[1][0] == 1.0


def test_undoing_a_deletion_brings_the_frames_back():
    doc = Document(hydrogen())
    doc.set_frames(H2_FRAMES, labels=["a", "b", "c"], energies=[1.0, 2.0, 3.0])
    doc.edit(lambda m: m.remove_atoms([1]))
    assert doc.frames is None
    doc.undo()
    assert doc.frames is not None and doc.frames.shape == (3, 2, 3)
    assert doc.frame_labels == ["a", "b", "c"] and doc.frame_energies == [1.0, 2.0, 3.0]
    doc.redo()
    assert doc.frames is None


def test_an_edit_during_playback_keeps_the_frames_as_they_were():
    doc = Document(hydrogen())
    doc.set_frames(H2_FRAMES)
    doc.set_position(0.5)                          # a blend of frames 0 and 1 is shown
    assert doc.mol.coords[1][0] == pytest.approx(0.87)
    doc.set_representation("tube")
    assert doc.frames[0][1][0] == 0.74 and doc.frames[1][1][0] == 1.0
    doc.undo()
    assert doc.frames[0][1][0] == 0.74
    doc.set_frame(0)
    assert doc.mol.coords[1][0] == 0.74


def test_an_animation_never_reaches_the_history():
    orig = np.array([[0, 0, 0], [0.74, 0, 0]], dtype=float)
    doc = Document(hydrogen(orig))
    doc.set_frames(H2_FRAMES)                      # the document's own frames
    anim = np.array([orig + np.array([[0, 0, 0], [d, 0, 0]]) for d in (0.0, 0.3, 0.0, -0.3)])
    doc.show_animation(anim, loop=True)
    assert doc.animating and doc.frame_count == 4 and doc.frames_loop
    doc.set_position(1.0)                          # displaced geometry shown
    assert doc.mol.coords[1][0] == pytest.approx(1.04)
    doc.set_representation("tube")                 # restyle during the animation
    doc.stop_animation()
    assert not doc.animating and doc.frame_count == 3
    assert doc.mol.coords[1][0] == pytest.approx(0.74)
    doc.undo()
    assert doc.mol.coords[1][0] == pytest.approx(0.74)
    assert doc.frames[0][1][0] == 0.74             # the document's frames never saw the animation


def test_an_edit_that_moves_atoms_ends_an_animation():
    doc = Document(hydrogen(np.array([[0, 0, 0], [0.74, 0, 0]], dtype=float)))
    doc.show_animation(np.array([doc.mol.coords + np.array([[0, 0, 0], [d, 0, 0]]) for d in (0.0, 0.3)]))
    doc.set_frame(1)
    events = recorder(doc)
    doc.edit(lambda m: m.set_positions([0], [[0.1, 0, 0]]), kind=Change.COORDS)
    assert not doc.animating and Change.FRAMES in events
    assert doc.mol.coords[1][0] == pytest.approx(0.74)   # the edit started from the document's geometry


def test_saving_during_playback_saves_the_frames_own_coordinates():
    doc = Document(hydrogen())
    doc.set_frames(H2_FRAMES)
    doc.set_position(0.5)
    info = doc.session_info()
    assert doc.mol.coords[1][0] == 0.74 and info.frames[0][1][0] == 0.74


def test_empty_frames_remove_the_trajectory():
    doc = Document(hydrogen())
    doc.set_frames(np.zeros((0, 2, 3)))
    assert doc.frame_count == 0 and doc.frames is None
    doc.edit(lambda m: m.set_positions([1], [[0.7, 0, 0]]), kind=Change.COORDS)
    doc.edit(lambda m: m.set_positions([1], [[0.8, 0, 0]]), kind=Change.COORDS)
    assert doc.mol.coords[1][0] == 0.8


def test_frames_that_do_not_fit_the_atoms_are_refused():
    doc = Document(hydrogen())
    with pytest.raises(ValueError):
        doc.set_frames(np.zeros((2, 3, 3)))


def test_a_trajectory_changes_one_frame_without_touching_the_old_one():
    first = Trajectory.from_array(H2_FRAMES)
    assert first is not None
    second = first.with_frame(1, [[0, 0, 0], [9, 0, 0]])
    assert first.coords(1)[1][0] == 1.0 and second.coords(1)[1][0] == 9.0
    assert second.frames[0] is first.frames[0]     # the other frames are shared
    with pytest.raises(ValueError):
        first.frames[0][0, 0] = 1.0                # read-only


# -- what an edit announces ------------------------------------------------------------------------

def test_a_failing_edit_leaves_nothing_behind():
    doc = Document(chain(3))
    doc.edit(lambda m: m.set_positions([0], [[0, 1, 0]]), kind=Change.COORDS)
    doc.undo()                                     # something to redo
    events = recorder(doc)

    def half_then_fail(mol: Molecule) -> None:
        mol.add_atom("C", [9.0, 0.0, 0.0])
        raise ValueError("boom")

    with pytest.raises(ValueError):
        doc.edit(half_then_fail)
    assert len(doc.mol) == 3 and not doc.can_undo() and doc.can_redo()
    assert events == [Change.STRUCTURE]            # the window learns the molecule was put back


def test_a_change_of_atoms_counts_as_structure_whatever_the_editor_says():
    doc = Document(chain(3))
    layout = doc.layout_revision
    events = recorder(doc)
    doc.edit(lambda m: m.add_atom("C", [9.0, 0.0, 0.0]), kind=Change.COORDS)
    assert events == [Change.STRUCTURE] and doc.layout_revision > layout


def test_undo_announces_what_the_step_changed():
    doc = Document(chain(3))
    doc.edit(lambda m: m.set_positions([0], [[0, 1, 0]]), kind=Change.COORDS)
    layout = doc.layout_revision
    events = recorder(doc)
    doc.undo()
    doc.redo()
    assert events == [Change.COORDS, Change.COORDS] and doc.layout_revision == layout


def test_named_selections_are_undoable_changes():
    doc = Document(chain(3))
    events = recorder(doc)
    doc.define("ends", "index 1 3")
    assert doc.selections == {"ends": "index 1 3"} and events == [Change.NAMES] and doc.modified
    doc.undo()
    assert doc.selections == {}
    doc.redo()
    assert doc.find("ends").tolist() == [0, 2]


def test_the_charge_model_is_an_undoable_change():
    doc = Document(chain(3))
    events = recorder(doc)
    doc.set_charge_model("mmff94")
    assert str(doc.charge_model) == "mmff94" and events == [Change.PROPERTIES]
    doc.undo()
    assert str(doc.charge_model) == "gasteiger"


def test_a_document_leaves_what_the_reader_returned_alone():
    frames = H2_FRAMES.copy()
    info = JobInfo(title="h2", frames=frames, frame_labels=["a", "b", "c"])
    first = Document(hydrogen(H2_FRAMES[0]), info)
    second = Document(hydrogen(H2_FRAMES[0]), info)
    assert first.frame_count == second.frame_count == 3 and info.frames is frames
    frames[0, 1, 0] = 5.0                          # the reader's array is theirs
    assert first.frames[0][1][0] == 0.74


# -- results ---------------------------------------------------------------------------------------

def water_with_results() -> tuple[Molecule, JobInfo]:
    mol = Molecule.from_columns([[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]], symbol=["O", "H", "H"])
    mol.add_bond(0, 1, 1)
    mol.add_bond(0, 2, 1)
    mode = VibrationalMode(1600.0, 70.0, np.zeros((3, 3)))
    grid = Volume("density", np.zeros(3), np.eye(3), np.ones((2, 2, 2), dtype=np.float32))
    info = JobInfo(vibrations=[mode], volumes=[grid], charges=np.array([-0.8, 0.4, 0.4]),
                   dipole=np.array([0.0, 1.8, 0.0]))
    return mol, info


def test_results_are_stale_after_the_structure_changes_and_back_after_undo():
    doc = Document(*water_with_results())
    assert doc.info.vibrations and doc.file_charges_apply() and not doc.stale_results()
    events = recorder(doc)
    doc.edit(lambda m: m.remove_atoms([2]))
    assert doc.stale_results() == ["grids", "normal modes", "charges", "dipole moment"]
    assert not doc.info.vibrations and not doc.info.volumes and doc.info.dipole is None
    assert not doc.file_charges_apply() and events[-1] == Change.RESULTS
    doc.undo()
    assert not doc.stale_results() and doc.info.vibrations and doc.file_charges_apply()


def test_moving_atoms_keeps_the_results():
    doc = Document(*water_with_results())
    doc.edit(lambda m: m.set_positions([1], [[1.0, 0, 0]]), kind=Change.COORDS)
    assert not doc.stale_results() and doc.info.vibrations


def test_stale_results_are_not_saved():
    doc = Document(*water_with_results())
    doc.edit(lambda m: m.set_column("symbol", "S", [0]))
    info = doc.session_info()
    assert info.charges is None and not info.vibrations and not info.volumes and info.dipole is None
    reopened = Document(doc.mol, info)
    assert str(reopened.charge_model) == "gasteiger" and not reopened.file_charges_apply()


def test_grids_are_added_and_removed_as_undoable_steps():
    doc = Document(*water_with_results())
    grid = Volume("potential", np.zeros(3), np.eye(3), np.zeros((2, 2, 2), dtype=np.float32))
    events = recorder(doc)
    assert doc.add_volume(grid) == 1 and [v.name for v in doc.info.volumes] == ["density", "potential"]
    assert events == [Change.RESULTS]
    doc.remove_volume(0)
    assert [v.name for v in doc.info.volumes] == ["potential"]
    doc.undo()
    doc.undo()
    assert [v.name for v in doc.info.volumes] == ["density"]


# -- the selection ---------------------------------------------------------------------------------

def test_the_selection_stays_on_its_atoms_through_edits_and_undo():
    doc = Document(chain(4))
    doc.select([2, 3])
    first = [int(u) for u in doc.mol.uids[[2, 3]]]
    doc.edit(lambda m: m.remove_atoms([0]))       # another atom goes
    assert doc.selected == [1, 2] and [int(u) for u in doc.mol.uids[doc.selected]] == first
    doc.select([0])
    doc.undo()                                     # the selection made since stays on its atom
    assert doc.selected == [1]
    doc.select([2, 3])
    doc.edit(lambda m: m.remove_atoms([2]))       # a selected atom goes
    assert doc.selected == [2]
    doc.edit(lambda m: m.remove_atoms([2]))       # and the other one
    assert doc.selected == []
    doc.undo()                                     # nothing selected is left: the step's selection comes back
    assert doc.selected == [2]
    doc.undo()                                     # that atom is still there, so the selection stays on it
    assert doc.selected == [3]


def test_the_word_selected_means_the_documents_selection():
    doc = Document(chain(4))
    doc.select([3, 1])
    assert doc.selected == [3, 1]                  # in the order chosen
    assert doc.find("selected").tolist() == [1, 3]
    with pytest.raises(IndexError):
        doc.select([9])


def test_delete_after_undo_removes_the_atom_selected(app, tmp_path):
    from miles.core.chem import from_smiles
    from miles.ui.mainwindow import MainWindow

    win = MainWindow(QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat))
    try:
        doc = win.api.new(from_smiles("CCO")[0], "ethanol")
        canvas = doc.canvas
        canvas.commands.delete_atoms([0])                   # the methyl carbon and its hydrogens
        oxygen = doc.molecule.symbols.index("O")
        doc.select([oxygen])
        win.current_canvas().commands.undo()
        assert [doc.molecule.symbols[i] for i in doc.selection] == ["O"]
        canvas.commands.delete_selection()
        assert "O" not in doc.molecule.symbols and doc.molecule.symbols.count("C") == 2
    finally:
        for c in win.canvases:
            c.doc.modified = False
        win.close()


# -- sessions --------------------------------------------------------------------------------------

def test_a_session_keeps_everything_exactly(tmp_path):
    mol, info = water_with_results()
    mol.add_property("tag", np.int64, -1, [7, -1, 3])
    doc = Document(mol, info)
    doc.set_frames(np.array([mol.coords, mol.coords + 0.123456789]), labels=["a", "b"])
    doc.set_frame(1)
    doc.define("hs", "element H")
    doc.select([2, 1])
    path = str(tmp_path / "w.miles")
    plugin_state = {"plugins": {"p": {"x": np.float32(0.5), "n": np.int64(3)}}}   # NumPy values
    session.write(path, doc.mol, doc.session_info(plugin_state))
    again = Document(*session.read(path))
    assert np.array_equal(again.mol.coords, doc.mol.coords)          # float64, nothing rounded
    assert np.array_equal(again.frames, doc.frames) and again.frame_labels == ["a", "b"] and again.frame == 1
    assert again.selections == {"hs": "element H"} and again.selected == [2, 1]
    assert again.mol.column("tag").tolist() == [7, -1, 3]
    assert again.info.vibrations and again.file_charges_apply()
    assert again.info.view["plugins"]["p"] == {"x": 0.5, "n": 3}


def test_a_failed_save_keeps_the_previous_file(tmp_path):
    doc = Document(chain(2))
    path = tmp_path / "c.miles"
    session.write(str(path), doc.mol, doc.session_info())
    before = path.read_bytes()
    with pytest.raises(FileFormatError):
        session.write(str(path), doc.mol, doc.session_info({"plugins": {"p": {"bad": object()}}}))
    assert path.read_bytes() == before and session.read(str(path))[0].symbols == ["C", "C"]
    assert [p.name for p in tmp_path.iterdir()] == ["c.miles"]     # no temporary file left behind


def test_a_damaged_session_is_a_format_error(tmp_path):
    path = tmp_path / "d.miles"
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        json.dump({"format": session.FORMAT, "version": 2, "symbols": ["C"]}, fh)
    with pytest.raises(FileFormatError):
        session.read(str(path))


def test_every_atom_column_is_saved(tmp_path):
    from miles.core.molecule import ATOM_COLUMNS

    mol = chain(2)
    mol.set_column("name", ["CA", "CB"])
    mol.set_column("resseq", [5, 6])
    path = str(tmp_path / "c.miles")
    session.write(path, mol, JobInfo())
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        saved = json.load(fh)["atoms"]
    assert set(saved) == set(ATOM_COLUMNS) - set(session.SKIPPED_COLUMNS)
    again, _ = session.read(path)
    assert again.column("name").tolist() == ["CA", "CB"] and again.column("resseq").tolist() == [5, 6]
