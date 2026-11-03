import numpy as np
import pytest

from miles.core import building
from miles.core.chem import from_smiles, to_smiles
from miles.core.forcefield import clean_atoms
from miles.core.fragments import default_element_fragment, library_fragments
from miles.core.geometry import distance
from miles.core.molecule import Molecule


def chain(n: int) -> Molecule:
    c = default_element_fragment("C")
    m = Molecule()
    building.place_fragment(m, c.mol, [0, 0, 0], origin=c.mol.coords[0])
    cur = 0
    for _ in range(n - 1):
        h = max(m.hydrogens_on(cur), key=lambda i: np.linalg.norm(m.coords[i]))
        building.attach_fragment(m, h, c.mol, c.attach)
        cur = h
    return m


def test_ring_closure_removes_hydrogens_and_relaxes():
    m = chain(6)
    assert m.formula() == "C6H14"
    ends = [i for i, a in enumerate(m.atoms)
            if a.symbol == "C" and sum(m.atoms[n].symbol == "C" for n in m.neighbors(i)) == 1]
    pair = building.join_atoms(m, *ends)
    assert m.formula() == "C6H12"
    clean_atoms(m, m.component(pair[0]), fold=True)
    assert distance(m.coords[pair[0]], m.coords[pair[1]]) == pytest.approx(1.54, abs=0.05)


def test_join_through_hydrogens():
    m = chain(5)
    ends = [i for i, a in enumerate(m.atoms)
            if a.symbol == "C" and sum(m.atoms[n].symbol == "C" for n in m.neighbors(i)) == 1]
    h1, h2 = m.hydrogens_on(ends[0])[0], m.hydrogens_on(ends[1])[0]
    building.join_atoms(m, h1, h2)
    assert m.formula() == "C5H10"


def test_best_hydrogen_follows_direction():
    m = Molecule()
    c = default_element_fragment("C")
    building.place_fragment(m, c.mol, [0, 0, 0], origin=c.mol.coords[0])
    for h in m.hydrogens_on(0):
        d = m.coords[h] - m.coords[0]
        assert building.best_hydrogen_towards(m, 0, d, lambda p: p) == h


def test_subset_keeps_internal_bonds():
    m = chain(3)
    heavy = [i for i, a in enumerate(m.atoms) if a.symbol == "C"]
    sub = m.subset(heavy)
    assert sub.formula() == "C3"
    assert len(sub.bonds) == 2


def test_ring_bonds():
    m = chain(4)
    assert not m.ring_bonds()


def test_smiles_with_attachment_point():
    mol, attach = from_smiles("*C(=O)OC")
    assert mol.atoms[attach].symbol == "H"
    assert mol.atoms[mol.neighbors(attach)[0]].symbol == "C"
    assert mol.formula() == "C2H4O2"


def test_smiles_round_trip():
    mol, _ = from_smiles("c1ccccc1O")
    assert to_smiles(mol) == "Oc1ccccc1"


def test_biomolecule_and_ligand_libraries():
    bio = library_fragments("bio")
    assert len(bio) == 25
    assert all(f.attach is not None for f in bio)
    ligands = {f.name: f for f in library_fragments("ligand")}
    assert ligands["Triphenylphosphine"].mol.atoms[ligands["Triphenylphosphine"].anchor].symbol == "P"
