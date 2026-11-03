import numpy as np
import pytest

from miles.core import building, internal, perception
from miles.core.elements import bond_length
from miles.core.forcefield import ForceField, clean_geometry
from miles.core.fragments import (
    all_group_fragments,
    all_ring_fragments,
    default_element_fragment,
    element_fragment,
    element_types,
    group_fragment,
    ring_fragment,
)
from miles.core.geometry import Geometry, angle, dihedral, distance, place_atom
from miles.core.molecule import Molecule


def kind(symbol, prefix):
    return next(k for k in element_types(symbol) if k.label.startswith(prefix))


def frag(symbol, prefix=None):
    if prefix is None:
        return default_element_fragment(symbol)
    return element_fragment(symbol, kind(symbol, prefix))


def placed(f):
    m = Molecule()
    building.place_fragment(m, f.mol, [0, 0, 0], origin=f.mol.coords[f.anchor])
    return m


def test_place_atom_matches_internal_coordinates():
    a, b, c = np.zeros(3), np.array([1.5, 0, 0]), np.array([2.0, 1.4, 0.3])
    for phi in (60.0, -120.0, 180.0):
        d = place_atom(a, b, c, 1.1, 109.5, phi)
        assert distance(d, a) == pytest.approx(1.1)
        assert angle(d, a, b) == pytest.approx(109.5)
        assert dihedral(d, a, b, c) == pytest.approx(phi)


def test_ethane_is_staggered():
    m = placed(frag("C"))
    building.attach_fragment(m, 1, frag("C").mol, 1)
    assert m.formula() == "C2H6"
    assert m.bond_order(0, 1) == 1.0
    h0, h1 = m.hydrogens_on(0)[0], m.hydrogens_on(1)
    torsions = sorted(abs(dihedral(m.coords[h0], m.coords[0], m.coords[1], m.coords[h])) for h in h1)
    assert torsions == pytest.approx([60.0, 60.0, 180.0], abs=1.0)


def test_trigonal_carbons_make_planar_double_bond():
    m = placed(frag("C", "Trigonal"))
    building.attach_fragment(m, 1, frag("C", "Trigonal").mol, 1)
    assert m.formula() == "C2H4"
    assert m.bond_order(0, 1) == 2.0
    assert distance(m.coords[0], m.coords[1]) == pytest.approx(1.34)
    for a in m.hydrogens_on(0):
        for b in m.hydrogens_on(1):
            assert abs(dihedral(m.coords[a], m.coords[0], m.coords[1], m.coords[b])) % 180 == \
                pytest.approx(0.0, abs=1e-6)


def test_linear_fragments_make_triple_bonds():
    m = placed(frag("C", "Linear"))
    building.attach_fragment(m, 1, frag("N", "Linear").mol, 1)
    assert m.formula() == "CHN"
    assert m.bond_order(0, 1) == 3.0


def test_attach_keeps_existing_indices():
    m = placed(frag("C"))
    before = m.coords[:1].copy()
    new = building.attach_fragment(m, 2, ring_fragment("benzene").mol, ring_fragment("benzene").attach)
    assert 2 in new
    assert m.atoms[2].symbol == "C"
    assert np.allclose(m.coords[:1], before)
    assert m.formula() == "C7H8"


def test_mutate_carbon_to_oxygen_gives_water():
    m = placed(frag("C"))
    building.mutate_atom(m, 0, "O", kind("O", "Bent"))
    assert m.formula() == "H2O"
    assert distance(m.coords[0], m.coords[1]) == pytest.approx(0.96)


def test_mutate_ring_carbon_to_pyridine_nitrogen():
    m = ring_fragment("benzene").mol
    building.mutate_atom(m, 0, "N", kind("N", "Bent"))
    assert m.formula() == "C5H5N"


def test_add_valence_and_delete():
    m = placed(frag("C"))
    building.add_valence(m, 0)
    assert m.formula() == "CH5"
    assert m.atoms[0].geom == "tbp"
    building.delete_atom(m, 0)
    assert len(m) == 0


def test_add_hydrogens_completes_heavy_atom_skeleton():
    m = ring_fragment("benzene").mol
    m.remove_atoms([i for i, a in enumerate(m.atoms) if a.symbol == "H"])
    for a in m.atoms:
        a.geom = a.nsites = a.valence = None
    perception.infer_atom_types(m)
    assert building.add_hydrogens(m) == 6
    assert m.formula() == "C6H6"


@pytest.mark.parametrize("value", [1.2, 1.9])
def test_set_bond_length_moves_one_side(value):
    m = placed(frag("C"))
    building.attach_fragment(m, 1, frag("C").mol, 1)
    internal.set_bond_length(m, 0, 1, value)
    assert distance(m.coords[0], m.coords[1]) == pytest.approx(value)
    h = m.hydrogens_on(1)[0]
    assert distance(m.coords[1], m.coords[h]) == pytest.approx(1.09)


def test_set_bond_length_in_ring():
    m = ring_fragment("benzene").mol
    internal.set_bond_length(m, 0, 1, 1.5)
    assert distance(m.coords[0], m.coords[1]) == pytest.approx(1.5)


def test_set_angle_and_dihedral():
    m = placed(frag("C"))
    building.attach_fragment(m, 1, frag("C").mol, 1)
    h0, h1 = m.hydrogens_on(0)[0], m.hydrogens_on(1)[0]
    for mode in (internal.MoveMode.LAST, internal.MoveMode.FIRST, internal.MoveMode.BOTH):
        internal.set_dihedral(m, h0, 0, 1, h1, 33.0, mode)
        assert dihedral(m.coords[h0], m.coords[0], m.coords[1], m.coords[h1]) == pytest.approx(33.0)
        internal.set_angle(m, h0, 0, 1, 115.0, mode)
        assert angle(m.coords[h0], m.coords[0], m.coords[1]) == pytest.approx(115.0)


def test_fragment_library_builds_with_ideal_bonds():
    # Clean uses MMFF94, which shortens single bonds next to trigonal atoms
    # (conjugation) by up to about 0.13 A relative to the table values.
    for f in all_ring_fragments() + all_group_fragments():
        m = f.mol
        assert f.attach is not None
        assert m.atoms[f.attach].symbol == "H"
        for (i, j), order in m.bonds.items():
            ideal = bond_length(m.atoms[i].symbol, m.atoms[j].symbol, order)
            conjugated = Geometry.TRIGONAL in (m.atoms[i].geom, m.atoms[j].geom)
            tolerance = 0.15 if conjugated and order == 1.0 else 0.06
            assert abs(distance(m.coords[i], m.coords[j]) - ideal) < tolerance, f.name


@pytest.mark.parametrize("key,formula", [
    ("naphthalene", "C10H8"), ("pyrrole", "C4H5N"), ("imidazole", "C3H4N2"),
    ("cyclohexane", "C6H12"), ("cyclopentadiene", "C5H6"),
])
def test_ring_formulas(key, formula):
    assert ring_fragment(key).mol.formula() == formula


@pytest.mark.parametrize("key,formula", [
    ("nitro", "HNO2"), ("ester", "C2H4O2"), ("sulfonic", "H2O3S"), ("tms", "C3H10Si"),
])
def test_group_formulas(key, formula):
    assert group_fragment(key).mol.formula() == formula


def test_forcefield_gradient_matches_finite_differences():
    m = group_fragment("ester").mol
    ff = ForceField(m)
    x = m.coords + np.random.default_rng(0).normal(0, 0.1, m.coords.shape)
    _, g = ff.energy_gradient(x)
    num = np.zeros_like(x)
    h = 1e-6
    for i in range(len(x)):
        for k in range(3):
            xp, xm = x.copy(), x.copy()
            xp[i, k] += h
            xm[i, k] -= h
            num[i, k] = (ff.energy_gradient(xp)[0] - ff.energy_gradient(xm)[0]) / (2 * h)
    assert np.abs(num - g).max() < 1e-5


def test_clean_restores_bond_lengths():
    m = ring_fragment("cyclohexane").mol
    m.coords = m.coords + np.random.default_rng(0).normal(0, 0.15, m.coords.shape)
    clean_geometry(m)
    for (i, j), order in m.bonds.items():
        ideal = bond_length(m.atoms[i].symbol, m.atoms[j].symbol, order)
        assert distance(m.coords[i], m.coords[j]) == pytest.approx(ideal, abs=0.03)


def test_perceive_bonds_and_orders():
    m = ring_fragment("benzene").mol.copy()
    m.set_bonds({})
    perception.perceive_bonds(m)
    perception.assign_bond_orders(m)
    assert len(m.bonds) == 12
    assert sorted(set(m.bonds.values())) == [1.0, 1.5]
