import numpy as np
import pytest

from miles.core import building, internal
from miles.core.fragments import (
    default_element_fragment,
    element_fragment,
    element_types,
    group_fragment,
    ring_fragment,
)
from miles.core.molecule import Molecule
from miles.core.symmetry import INF, find_point_group, symmetrize


def frag(symbol, prefix=None):
    if prefix is None:
        return default_element_fragment(symbol)
    return element_fragment(symbol, next(k for k in element_types(symbol) if k.label.startswith(prefix)))


def placed(f):
    m = Molecule()
    building.place_fragment(m, f.mol, [0, 0, 0], origin=f.mol.coords[f.anchor])
    return m


def halogen():
    return next(k for k in element_types("F") if k.bare)


def ethane():
    m = placed(frag("C"))
    building.attach_fragment(m, 1, frag("C").mol, 1)
    return m


def allene():
    m = placed(frag("C", "Linear"))
    building.attach_fragment(m, 1, frag("C", "Trigonal").mol, 1)
    building.attach_fragment(m, 2, frag("C", "Trigonal").mol, 1)
    return m


def sf6():
    m = placed(frag("S", "Octahedral"))
    for h in m.hydrogens_on(0):
        building.mutate_atom(m, h, "F", halogen())
    return m


def co2():
    m = placed(frag("C", "Linear"))
    building.attach_fragment(m, 1, frag("O", "Double").mol, 1)
    building.attach_fragment(m, 2, frag("O", "Double").mol, 1)
    return m


def hydrogen_peroxide():
    m = placed(frag("O"))
    building.attach_fragment(m, 1, frag("O").mol, 1)
    a, b = m.hydrogens_on(0)[0], m.hydrogens_on(1)[0]
    internal.set_dihedral(m, a, 0, 1, b, 115.0)
    return m


CASES = [
    (lambda: placed(frag("C")), "Td"),
    (ethane, "D3d"),
    (lambda: ring_fragment("benzene").mol, "D6h"),
    (lambda: placed(frag("O")), "C2v"),
    (lambda: placed(frag("N")), "C3v"),
    (allene, "D2d"),
    (sf6, "Oh"),
    (co2, f"D{INF}h"),
    (lambda: ring_fragment("cyclohexane").mol, "D3d"),
    (lambda: ring_fragment("naphthalene").mol, "D2h"),
    (lambda: ring_fragment("cyclopropane").mol, "D3h"),
    (hydrogen_peroxide, "C2"),
    (lambda: ring_fragment("pyridine").mol, "C2v"),
    (lambda: placed(frag("P", "Trigonal bip")), "D3h"),
    (lambda: placed(frag("Pt", "Square")), "D4h"),
    (lambda: group_fragment("tbutyl").mol, "C3v"),
]


@pytest.mark.parametrize("make,expected", CASES)
def test_point_groups(make, expected):
    m = make()
    assert find_point_group(m.symbols, m.coords, 0.1).name == expected


def test_symmetrize_noisy_structure():
    m = ring_fragment("benzene").mol
    noisy = m.coords + np.random.default_rng(1).normal(0, 0.04, m.coords.shape)
    assert find_point_group(m.symbols, noisy, 0.01).name == "C1"
    coords, name = symmetrize(m.symbols, noisy, 0.15)
    assert name == "D6h"
    assert find_point_group(m.symbols, coords, 0.001).name == "D6h"
