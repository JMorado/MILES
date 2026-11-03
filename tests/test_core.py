import itertools

import numpy as np
import pytest

from miles.core import chem, internal, perception
from miles.core.elements import ELEMENTS, bond_length, element, normalize_symbol
from miles.core.forcefield import CleanMode, clean_geometry
from miles.core.fragments import (
    Category,
    all_group_fragments,
    all_ring_fragments,
    molecule_fragment,
    ring_fragment,
)
from miles.core.geometry import distance
from miles.core.neighbors import near_mask, pairs_within

# -- elements -------------------------------------------------------------------------

@pytest.mark.parametrize(("symbol", "row", "column"), [
    ("H", 1, 1), ("He", 1, 18), ("C", 2, 14), ("Fe", 4, 8), ("La", 8, 3), ("Lu", 6, 3),
    ("Ac", 9, 3), ("Lr", 7, 3), ("Rn", 6, 18),
])
def test_periodic_table_layout(symbol, row, column):
    assert (element(symbol).row, element(symbol).column) == (row, column)


@pytest.mark.parametrize(("text", "symbol"), [
    ("c", "C"), ("CL", "Cl"), ("C12", "C"), ("Fe(2+)", "Fe"), ("26", "Fe"), ("x1", "X"), ("Q", None),
    ("", None),
])
def test_normalize_symbol(text, symbol):
    assert normalize_symbol(text) == symbol


def test_bond_lengths_shrink_with_order():
    assert bond_length("C", "C", 1) > bond_length("C", "C", 1.5) > bond_length("C", "C", 2)
    assert bond_length("N", "C") == bond_length("C", "N")
    radii = ELEMENTS["Fe"].covalent_radius + ELEMENTS["Cl"].covalent_radius
    assert bond_length("Fe", "Cl") == pytest.approx(radii)


# -- neighbour search -----------------------------------------------------------------

def test_grid_search_matches_brute_force():
    rng = np.random.default_rng(1)
    points = rng.uniform(-6, 6, size=(300, 3))
    cutoff = 1.7
    found = {tuple(p) for p in pairs_within(points, cutoff).tolist()}
    expected = {(i, j) for i, j in itertools.combinations(range(len(points)), 2)
                if np.linalg.norm(points[i] - points[j]) < cutoff}
    assert found == expected

    centers = points[:5]
    mask = near_mask(points, centers, 2.5)
    brute = np.linalg.norm(points[:, None] - centers[None], axis=-1).min(axis=1) < 2.5
    assert np.array_equal(mask, brute)


def test_grid_search_handles_empty_input():
    assert pairs_within(np.zeros((0, 3)), 2.0).shape == (0, 2)
    assert not near_mask(np.zeros((3, 3)), np.zeros((0, 3)), 2.0).any()


# -- Open Babel bridge ------------------------------------------------------------------

def test_smiles_round_trip():
    mol, attach = chem.from_smiles("CCO")
    assert attach is None
    assert mol.formula() == "C2H6O"
    assert chem.to_smiles(mol) == "CCO"


def test_smiles_attachment_point_becomes_a_hydrogen():
    mol, attach = chem.from_smiles("*C(=O)OC")
    assert attach is not None and mol.atoms[attach].symbol == "H"
    frag = molecule_fragment(mol, "ester", "*C(=O)OC", attach)
    assert frag.category == Category.SMILES and frag.anchor == mol.neighbors(attach)[0]


def test_invalid_smiles_raises():
    assert not chem.is_smiles("C C")
    with pytest.raises(chem.ChemError):
        chem.from_smiles("not smiles")


def test_obmol_conversion_keeps_aromatic_bonds():
    benzene = ring_fragment("benzene").mol
    back = chem.from_openbabel(chem.to_openbabel(benzene))
    assert back.formula() == benzene.formula()
    assert sorted(back.bonds.values()) == sorted(benzene.bonds.values())


# -- perception -------------------------------------------------------------------------

def test_rebond_recovers_bonds_and_orders():
    mol, _ = chem.from_smiles("C=CC(=O)O")
    expected = dict(mol.bonds)
    mol.set_bonds({})
    perception.rebond(mol)
    assert mol.bonds == expected
    carbonyl_c = next(i for i, a in enumerate(mol.atoms) if a.symbol == "C" and mol.order_sum(i) == 4
                      and any(mol.atoms[n].symbol == "O" for n in mol.neighbors(i)))
    assert mol.atoms[carbonyl_c].geom == "trigonal"


# -- internal coordinates ---------------------------------------------------------------

def _butane():
    mol, _ = chem.from_smiles("CCCC")
    c = [i for i, a in enumerate(mol.atoms) if a.symbol == "C"]
    return mol, c


@pytest.mark.parametrize("mode", list(internal.MoveMode))
def test_bond_length_move_modes(mode):
    mol, (_, c2, c3, _) = _butane()
    before = mol.coords.copy()
    internal.set_bond_length(mol, c2, c3, 1.80, mode)
    assert distance(mol.coords[c2], mol.coords[c3]) == pytest.approx(1.80)
    moved = {i for i in range(len(mol)) if not np.allclose(mol.coords[i], before[i])}
    first_side = mol.component(c2, blocked={c3})
    last_side = mol.component(c3, blocked={c2})
    if mode == internal.MoveMode.LAST:
        assert moved == last_side
    elif mode == internal.MoveMode.FIRST:
        assert moved == first_side
    elif mode == internal.MoveMode.BOTH:
        assert moved == first_side | last_side
    else:
        assert moved == {c3, *mol.terminal_neighbors(c3)}


def test_dihedral_in_a_ring_moves_only_the_end_atom():
    mol = ring_fragment("cyclohexane").mol
    ring = [i for i, a in enumerate(mol.atoms) if a.symbol == "C"]
    a, b, c, d = ring[0], ring[1], ring[2], ring[3]
    path = [a, b, c, d]
    for x, y in itertools.pairwise(path):
        assert mol.bond_order(x, y)
    before = mol.coords.copy()
    internal.set_dihedral(mol, a, b, c, d, 20.0)
    moved = {i for i in range(len(mol)) if not np.allclose(mol.coords[i], before[i])}
    assert moved and moved <= {d, *mol.terminal_neighbors(d)}


# -- library and clean -------------------------------------------------------------------

def test_library_sizes_and_categories():
    rings, groups = all_ring_fragments(), all_group_fragments()
    assert len(rings) == 16 and len(groups) == 23
    assert {f.category for f in rings} == {Category.RING}
    assert {f.category for f in groups} == {Category.GROUP}
    assert len({f.key for f in rings + groups}) == len(rings) + len(groups)


def test_clean_mode_depends_on_residue_data():
    mol, _ = chem.from_smiles("CCO")
    assert clean_geometry(mol) == CleanMode.MMFF
    mol.atoms[0].symbol = "Si"
    assert clean_geometry(mol) == CleanMode.FULL
    mol.atoms[0].symbol = "C"
    mol.atoms[0].resname = "LIG"
    assert clean_geometry(mol) == CleanMode.GENTLE
