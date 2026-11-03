import numpy as np
import pytest

from miles.builtin_plugins.superpose.matching import apply, match_atoms, superpose
from miles.core.chem import from_smiles
from miles.core.geometry import rotation_matrix
from miles.core.molecule import Atom, Molecule

ROTATION = rotation_matrix([0.3, -1.0, 0.5], 1.1)
SHIFT = np.array([4.0, -2.5, 7.0])


def moved(mol: Molecule) -> Molecule:
    copy = mol.copy()
    copy.coords = copy.coords @ ROTATION.T + SHIFT
    return copy


def permuted(mol: Molecule, order: list[int]) -> Molecule:
    position = {old: new for new, old in enumerate(order)}
    out = Molecule.from_atoms([mol.atoms[i] for i in order], mol.coords[order])
    for (i, j), bond in mol.bonds.items():
        out.add_bond(position[i], position[j], bond)
    return out


def ca_trace(names: list[str]) -> Molecule:
    t = np.arange(len(names))
    coords = np.stack([2.3 * np.cos(1.75 * t), 2.3 * np.sin(1.75 * t), 1.5 * t], axis=1)
    atoms = [Atom("C", name="CA", resname=name, resseq=k + 1, chain="A") for k, name in enumerate(names)]
    return Molecule.from_atoms(atoms, coords)


def test_rigid_copy_overlays_exactly():
    mol = from_smiles("CC(=O)Nc1ccc(O)cc1")[0]
    target = moved(mol)
    fit = superpose(mol, target)
    assert fit.method == "all atoms in order"
    assert fit.rmsd < 1e-6
    assert np.allclose(fit.rotation, ROTATION, atol=1e-6)
    apply(mol, fit)
    assert np.allclose(mol.coords, target.coords, atol=1e-6)


def test_renumbered_copy_is_matched_by_graph():
    mol = from_smiles("CC(=O)Nc1ccc(O)cc1")[0]
    order = list(np.random.default_rng(3).permutation(len(mol)))
    target = moved(permuted(mol, order))
    pairs, method = match_atoms(mol, target)
    assert method == "heavy atoms by element and graph"
    assert len(pairs) == sum(a.symbol != "H" for a in mol.atoms)
    assert all(mol.atoms[m].symbol == target.atoms[t].symbol for m, t in pairs)
    assert superpose(mol, target).rmsd < 1e-6


def test_symmetric_molecule_is_matched():
    mol = from_smiles("c1ccccc1")[0]
    order = list(np.random.default_rng(5).permutation(len(mol)))
    target = moved(permuted(mol, order))
    assert superpose(mol, target).rmsd < 1e-6


def test_proteins_are_paired_by_sequence():
    names = ["MET", "ALA", "GLY", "SER", "TRP", "LYS", "ASP", "LEU"]
    target = ca_trace(names)
    mobile = moved(target)
    mobile.remove_atoms([4])
    fit = superpose(mobile, target)
    assert fit.method == "CA atoms by sequence"
    assert len(fit.pairs) == 7
    assert all(mobile.atoms[m].resname == target.atoms[t].resname for m, t in fit.pairs)
    assert fit.rmsd < 1e-6


def test_different_molecules_cannot_be_matched():
    with pytest.raises(ValueError, match="cannot be"):
        superpose(from_smiles("c1ccccc1")[0], from_smiles("CCO")[0])


def test_too_few_pairs():
    mol = from_smiles("CCO")[0]
    with pytest.raises(ValueError, match="three"):
        superpose(mol, mol, pairs=[(0, 0), (1, 1)])
