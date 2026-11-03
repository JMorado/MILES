"""Hydrogen bonds within a structure, and their display."""

import time

import numpy as np
import pytest
from PySide6.QtCore import QSettings

from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence
from miles.builtin_plugins.shared.hbonds import HydrogenBondFinder, hydrogen_bonds
from miles.core.chem import from_smiles
from miles.core.molecule import Molecule
from miles.files.common import JobInfo
from miles.render.theme import Background
from miles.ui.canvas.view_state import view_state
from miles.ui.mainwindow import MainWindow


def water_dimer(oo=2.9, tilt=0.0):
    """Two waters; the first donates along x to the second, its hydrogen tilted by ``tilt`` degrees."""
    mol = Molecule()
    t = np.radians(tilt)
    o1 = mol.add_atom("O", [0.0, 0.0, 0.0])
    h1 = mol.add_atom("H", [0.96 * np.cos(t), 0.96 * np.sin(t), 0.0])
    h2 = mol.add_atom("H", [-0.24, -0.93, 0.0])
    o2 = mol.add_atom("O", [oo, 0.0, 0.0])
    h3 = mol.add_atom("H", [oo + 0.3, 0.9, 0.0])
    h4 = mol.add_atom("H", [oo + 0.3, -0.45, 0.78])
    for a, b in ((o1, h1), (o1, h2), (o2, h3), (o2, h4)):
        mol.add_bond(a, b)
    return mol


def without_hydrogens(mol):
    return mol.subset([i for i, a in enumerate(mol.atoms) if a.symbol != "H"])


def test_water_dimer_has_one_bond():
    (bond,) = hydrogen_bonds(water_dimer())
    assert (bond.donor, bond.hydrogen, bond.acceptor) == (0, 1, 3)
    assert bond.distance == pytest.approx(2.9)


@pytest.mark.parametrize(("oo", "tilt"), [(2.2, 0.0), (3.8, 0.0), (2.9, 75.0)])
def test_no_bond_when_too_close_too_far_or_bent(oo, tilt):
    assert hydrogen_bonds(water_dimer(oo, tilt)) == []


def test_atoms_left_out_take_no_part():
    assert hydrogen_bonds(water_dimer(), atoms=[0, 1, 2]) == []


def test_atoms_close_in_the_bond_graph_are_skipped():
    mol, _ = from_smiles("OCCO")   # the hydroxyls are 1-4 related, whatever the conformation
    assert hydrogen_bonds(mol) == []


def backbone_offsets(mol):
    found = []
    for b in hydrogen_bonds(mol):
        donor, acceptor = mol.atoms[b.donor], mol.atoms[b.acceptor]
        if donor.name == "N" and acceptor.name == "O":
            found.append(donor.resseq - acceptor.resseq)
    return found


@pytest.mark.parametrize("strip", [False, True])
def test_alpha_helix_bonds_are_i_to_i_plus_4(strip):
    mol = build_sequence("AEAAAKEAAAKAAEAAAKA", SequenceKind.HELIX)
    offsets = backbone_offsets(without_hydrogens(mol) if strip else mol)
    assert len(offsets) == 15   # every carbonyl but the last four
    assert set(offsets) == {4}


@pytest.mark.parametrize("strip", [False, True])
def test_watson_crick_pairs(strip):
    mol = build_sequence("CGCGAATTCGCG", SequenceKind.DNA)
    if strip:
        mol = without_hydrogens(mol)
    pairs = {(mol.atoms[b.donor].name, mol.atoms[b.acceptor].name) for b in hydrogen_bonds(mol)
             if mol.atoms[b.donor].chain != mol.atoms[b.acceptor].chain}
    count = sum(mol.atoms[b.donor].chain != mol.atoms[b.acceptor].chain for b in hydrogen_bonds(mol))
    assert count == 8 * 3 + 4 * 2   # eight G-C pairs with three bonds, four A-T pairs with two
    assert pairs == {("N4", "O6"), ("N1", "N3"), ("N2", "O2"), ("N6", "O4"), ("N3", "N1")}


def add_water_facing(mol, atom, distance=3.0):
    """Put a water with one hydrogen pointing at ``atom``, along its free direction."""
    centre = mol.coords[atom]
    out = sum(centre - mol.coords[j] for j in mol.neighbors(atom))
    if np.linalg.norm(out) < 1e-3:
        out = np.cross(*(mol.coords[j] - centre for j in mol.neighbors(atom)[:2]))
    out = out / np.linalg.norm(out)
    o = mol.add_atom("O", centre + distance * out)
    h = mol.add_atom("H", centre + (distance - 0.96) * out)
    h2 = mol.add_atom("H", mol.coords[o] + 0.96 * np.array([0.0, 0.0, 1.0]))
    mol.add_bond(o, h)
    mol.add_bond(o, h2)
    return o


@pytest.mark.parametrize(("smiles", "accepts"), [("CN(C)C", True), ("c1ccncc1", True), ("CNC(C)=O", False),
                                                  ("Nc1ccccc1", False), ("C[N+](C)(C)C", False)])
def test_which_nitrogens_accept(smiles, accepts):
    mol, _ = from_smiles(smiles)
    n = next(i for i, a in enumerate(mol.atoms) if a.symbol == "N")
    water = add_water_facing(mol, n)
    found = [(b.donor, b.acceptor) for b in hydrogen_bonds(mol)]
    assert ((water, n) in found) == accepts


def test_large_structure_is_fast():
    mol = build_sequence("AEAAAKEAAAKA" * 30, SequenceKind.HELIX)   # about 4600 atoms
    start = time.perf_counter()
    finder = HydrogenBondFinder(mol)
    setup = time.perf_counter() - start
    start = time.perf_counter()
    found = finder.find(mol.coords)
    search = time.perf_counter() - start
    assert len(found) > 300
    assert setup < 0.5 and search < 0.1


# -- display ------------------------------------------------------------------------------------

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


def show_bonds(window, on=True, **style):
    """Switch the hydrogen bonds of the current tab on or off, with other display options."""
    window.api.current.style(hydrogen_bonds=on, **style)


def bond_lines(canvas, export=False):
    """The lines the Hydrogen Bonds plugin draws (the only annotations in these tests)."""
    layers = canvas.drawing.annotations(export=export)
    return layers[0].lines if layers else []


def ends(canvas):
    return [[line.start, line.end] for line in bond_lines(canvas)]


def points(canvas):
    """End points in angstrom of the lines drawn, shape (n, 2, 3)."""
    layers = canvas.drawing.annotations()
    resolved = canvas.drawing.resolve_lines(layers[0]) if layers else []
    return np.array([[a, b] for _, a, b in resolved]).reshape(-1, 2, 3)


def test_lines_start_at_drawn_hydrogens(window):
    window.open_structure(water_dimer(), JobInfo(title="water dimer"))
    canvas = window.current_canvas()
    show_bonds(window)
    assert ends(canvas) == [[1, 3]]   # from the hydrogen to the acceptor
    canvas.set_style(show_hydrogens=False)
    assert ends(canvas) == [[0, 3]]   # from the donor when hydrogens are hidden
    assert "1 hydrogen bond" in canvas.info_text


def test_toggle_draws_hydrogen_bonds_and_follows_frames(window):
    mol = build_sequence("AEAAAKEAAAKA", SequenceKind.HELIX)
    window.open_structure(mol, JobInfo(title="helix"))
    canvas = window.current_canvas()
    toggle = window.plugins.toggles["hydrogen_bonds"]
    window.api.run("View > Hydrogen Bonds")
    assert toggle.value(window.api.current) and toggle.action.isChecked()
    bonds = bond_lines(canvas)
    assert len(bonds) == len(hydrogen_bonds(canvas.mol)) >= 8   # the backbone's i -> i+4 bonds and more
    assert f"{len(bonds)} hydrogen bonds" in canvas.info_text
    canvas.grab()
    image = canvas.render_export(1.0, Background.WHITE, False)
    assert not image.isNull()
    assert view_state(canvas)["plugins"]["miles.builtin_plugins.hydrogen_bonds"]["hydrogen_bonds"]

    far = canvas.mol.coords.copy()
    first_half = np.array([a.resseq <= 6 for a in canvas.mol.atoms])
    far[first_half] -= 20.0   # pull the helix apart in a second frame
    canvas.doc.set_frames(np.array([canvas.mol.coords.copy(), far]))
    canvas.doc.set_frame(1)
    assert len(bond_lines(canvas, export=True)) < len(bonds)
    canvas.doc.set_frame(0)
    assert len(bond_lines(canvas, export=True)) == len(bonds)

    window.api.run("View > Hydrogen Bonds")
    assert not toggle.value(window.api.current) and "hydrogen bond" not in canvas.info_text


def follow_view_style(canvas):
    """Drop the per-atom representations that proteins and nucleic acids open with."""
    canvas.doc.set_representation("")


def test_lines_meet_the_ribbon_in_cartoon_view(window):
    from miles.core.molecule import Representation

    mol = build_sequence("AEAAAKEAAAKAAEAAAKA", SequenceKind.HELIX)
    window.open_structure(mol, JobInfo(title="helix"))
    canvas = window.current_canvas()
    follow_view_style(canvas)
    show_bonds(window, mode=Representation.BALL)
    canvas.drawing.frame()
    ball = points(canvas)
    atoms = canvas.mol.coords
    assert np.allclose(ball, atoms[np.array(ends(canvas))])   # drawn atoms: the lines end at them

    canvas.set_style(mode=Representation.CARTOON)
    canvas.drawing.frame()   # builds the cartoon mesh
    mesh = canvas.scenes.base.mesh
    start = time.perf_counter()
    lines = points(canvas)
    assert time.perf_counter() - start < 0.2
    assert len(lines) >= 12   # the backbone's i -> i+4 bonds stay visible
    vertices = mesh.vertices.astype(float)
    for point in lines.reshape(-1, 3):   # every end lies on the ribbon surface
        assert np.linalg.norm(vertices - point, axis=1).min() < 1e-6
    backbone = atoms[[i for i, a in enumerate(canvas.mol.atoms) if a.name in ("N", "O")]]
    for point in lines.reshape(-1, 3):   # and not at the hidden backbone atoms
        assert np.linalg.norm(backbone - point, axis=1).min() > 0.05


def test_base_pairs_inside_a_cartoon_are_left_out(window):
    from miles.core.molecule import Representation

    mol = build_sequence("CGCGAATTCGCG", SequenceKind.DNA)
    window.open_structure(mol, JobInfo(title="DNA"))
    canvas = window.current_canvas()
    follow_view_style(canvas)
    show_bonds(window, mode=Representation.TUBE)
    canvas.drawing.frame()
    assert len(points(canvas)) == 32   # Watson-Crick pairs between drawn bases
    canvas.set_style(mode=Representation.CARTOON)
    canvas.drawing.frame()
    assert len(points(canvas)) < 32    # bases are not drawn by the cartoon


def test_status_says_how_many_lines_a_cartoon_leaves_out(window):
    from miles.core.molecule import Representation

    mol = build_sequence("CGCGAATTCGCG", SequenceKind.DNA)
    window.open_structure(mol, JobInfo(title="DNA"))
    canvas = window.current_canvas()
    follow_view_style(canvas)
    show_bonds(window, mode=Representation.TUBE)
    assert canvas.info_text.endswith("32 hydrogen bonds")
    canvas.set_style(mode=Representation.CARTOON)
    assert "32 hydrogen bonds (" in canvas.info_text and " drawn)" in canvas.info_text
