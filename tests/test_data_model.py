"""The column model: shared storage, read-only arrays, properties, display state and undo."""

import tracemalloc

import numpy as np
import pytest

from miles.core.display import AtomDisplay
from miles.core.document import CHARGE, Change, Document
from miles.core.fragments import ring_fragment
from miles.core.molecule import Atom, Molecule
from miles.files import session
from miles.files.common import JobInfo
from miles.files.formats import parse_text

WATER = "3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"

PDB = """\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 21.50           N
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  0.50 17.25           C
ATOM      3  C   ALA A   1      13.161   5.902  -5.134  1.00 12.00           C
ATOM      4  O   ALA A   1      13.727   5.395  -6.104  1.00 30.75           O
HETATM    5  O   HOH A 101      20.000   5.000  -5.000  1.00 45.00           O
END
"""


def big(n: int) -> Molecule:
    coords = np.random.default_rng(0).random((n, 3)) * 50
    return Molecule.from_columns(coords, symbol=["C"] * n, name=["CA"] * n, resname=["ALA"] * n,
                                 resseq=np.arange(n) // 10, chain=["A"] * n)


def test_copies_share_storage_until_one_of_them_changes():
    mol = parse_text(WATER)
    copy = mol.copy()
    assert np.shares_memory(copy.coords, mol.coords)
    mol.set_positions([0], [1.0, 2.0, 3.0])
    assert np.allclose(copy.coords[0], 0.0) and np.allclose(mol.coords[0], [1, 2, 3])
    mol.atoms[1].symbol = "F"
    mol.add_atom("Cl", [5.0, 0.0, 0.0])
    assert copy.symbols == ["O", "H", "H"] and mol.symbols == ["O", "F", "H", "Cl"]


def test_arrays_handed_out_are_read_only():
    mol = parse_text(WATER)
    for array in (mol.coords, mol.numbers, mol.column("resname"), mol.bond_pairs, mol.bond_orders):
        with pytest.raises(ValueError):
            array[0] = array[0]
    mol.coords = mol.coords + 1.0   # whole columns are replaced, not written into
    assert np.allclose(mol.coords[0], 1.0)


def test_an_undo_step_costs_memory_only_for_what_changed():
    doc = Document(big(100_000))
    tracemalloc.start()
    before = tracemalloc.get_traced_memory()[0]
    doc.checkpoint()
    snapshot = tracemalloc.get_traced_memory()[0] - before
    doc.mol.set_positions([0], [0.0, 0.0, 0.0])   # copies the coordinates, nothing else
    after_edit = tracemalloc.get_traced_memory()[0] - before
    tracemalloc.stop()
    assert snapshot < 50_000
    coordinates = 100_000 * 3 * 8
    assert coordinates <= after_edit < coordinates + 200_000


def test_the_topology_token_follows_atoms_elements_and_bonds():
    mol = parse_text(WATER)
    token = mol.topology
    mol.coords = mol.coords * 1.1
    mol.atoms[0].name = "OW"
    assert mol.topology == token and mol.copy().topology == token
    mol.add_bond(0, 1, 2.0)
    assert mol.topology != token
    token = mol.topology
    mol.atoms[2].symbol = "D"
    assert mol.topology != token


def test_properties_follow_atoms_and_are_saved_in_sessions():
    mol = ring_fragment("benzene").mol
    mol.set_property("mass_tag", np.arange(len(mol), dtype=float))
    mol.add_property("label", object, "", ["a"] * len(mol))
    mol.remove_atoms([0])
    assert mol.column("mass_tag")[0] == 1.0
    start = mol.add_atom("C", [9.0, 0.0, 0.0])
    assert np.isnan(mol.column("mass_tag")[start]) and mol.column("label")[start] == ""
    mol.set_property(CHARGE, np.zeros(len(mol)))
    mol.property_sources[CHARGE] = "worked out"
    data = session.to_dict(mol, JobInfo())
    assert set(data["properties"]) == {"mass_tag", "label"}   # worked-out ones are left out
    again, _ = session.from_dict(data)
    assert np.array_equal(again.column("mass_tag"), mol.column("mass_tag"), equal_nan=True)
    assert again.column("label").tolist() == mol.column("label").tolist()


def test_pdb_files_give_occupancy_and_b_factors(tmp_path):
    from miles.files.formats import read_file

    path = tmp_path / "ala.pdb"
    path.write_text(PDB, encoding="utf-8")
    mol, _ = read_file(str(path))
    assert mol.has_property("b_factor") and mol.has_property("occupancy")
    by_name = dict(zip(mol.column("name").tolist(), mol.column("b_factor").tolist(), strict=True))
    assert by_name["CA"] == 17.25 and by_name["N"] == 21.5
    water = np.flatnonzero(mol.mask("resname", {"HOH"}))
    assert mol.column("b_factor")[water].tolist() == [45.0]
    assert mol.column("occupancy")[mol.column("name") == "CA"].tolist() == [0.5]


def test_the_display_stays_with_its_atoms():
    mol = ring_fragment("benzene").mol
    display = AtomDisplay(mol.uids).with_rep("tube", [3]).with_hidden(np.arange(len(mol)) == 5)
    mol.remove_atoms([0, 1])
    aligned = display.aligned(mol)
    assert aligned.rep[1] == "tube" and aligned.hidden[3] and aligned.rep.tolist().count("tube") == 1


def test_new_hydrogens_take_the_representation_of_their_atom():
    mol = Molecule()
    carbon = mol.add_atom(Atom("C", "tetrahedral", 4, 4), [0.0, 0.0, 0.0])
    display = AtomDisplay(mol.uids).with_rep("spacefill", [carbon])
    h = mol.add_atom("H", [1.09, 0.0, 0.0])
    mol.add_bond(carbon, h)
    other = mol.add_atom("N", [3.0, 0.0, 0.0])
    aligned = display.aligned(mol)
    assert aligned.rep[h] == "spacefill" and aligned.rep[other] == ""


def test_charges_are_a_column_kept_until_the_bonds_change(monkeypatch):
    import miles.core.document as document

    computed = []
    real = document.partial_charges
    monkeypatch.setattr(document, "partial_charges", lambda *args: computed.append(1) or real(*args))
    doc = Document(parse_text(WATER))
    doc.set_charge_model("gasteiger")
    charges = doc.charges()
    assert charges is not None and np.array_equal(doc.mol.column(CHARGE), charges)
    doc.edit(lambda mol: mol.set_positions([1], [0.9, 0.1, 0.0]))
    assert np.array_equal(doc.charges(), charges) and len(computed) == 1   # moving atoms keeps them
    doc.edit(lambda mol: mol.remove_atoms([2]))
    assert doc.mol.property_sources.get(CHARGE) != (doc.mol.topology, "gasteiger")
    fewer = doc.charges()
    assert fewer is not None and len(fewer) == 2
    doc.undo()
    assert len(doc.charges()) == 3   # the snapshot's charges still fit its atoms


def test_change_events_say_what_changed():
    doc = Document(parse_text(WATER))
    events = []
    doc.listeners.append(events.append)
    doc.checkpoint()
    doc.mol.set_positions([1], [1.0, 0.0, 0.0])
    doc.notify(Change.COORDS, atoms=[1])
    doc.set_charge_model("mmff94")
    assert events[0].kind == Change.COORDS and events[0].atoms.tolist() == [1]
    assert events[1].kind == Change.PROPERTIES and events[1].properties == (CHARGE,)


def test_frames_move_the_atoms_but_keep_the_scene_layout(app):
    from miles.render.scene import Style
    from miles.render.theme import LIGHT
    from miles.ui.canvas.scene_cache import SceneCache

    mol = ring_fragment("benzene").mol
    doc = Document(mol, JobInfo(frames=np.stack([mol.coords, mol.coords + 1.0])))
    cache = SceneCache()

    def frame():
        return cache.frame(doc.mol, doc.revision, Style(), LIGHT, {}, [], display=doc.display,
                           layout_revision=doc.layout_revision)

    first = frame().spheres[:, :3].copy()
    layout = cache._layout
    doc.set_frame(1)
    moved = frame().spheres[:, :3]
    assert cache._layout is layout and np.allclose(moved, first + 1.0)
    doc.edit(lambda m: m.remove_atoms([0]))
    frame()
    assert cache._layout is not layout and len(cache._layout.reps) == len(doc.mol)
