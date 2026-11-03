"""Partial charges, formal charges and the dipole moment."""

from pathlib import Path

import numpy as np
import pytest

from miles.core.charges import ChargeModel, dipole_moment, formal_charges, partial_charges
from miles.core.chem import ChemError, from_smiles
from miles.files import session
from miles.files.formats import read_file

DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize("model", [ChargeModel.GASTEIGER, ChargeModel.MMFF94])
@pytest.mark.parametrize("smiles", ["O", "CCO", "c1ccncc1"])
def test_neutral_molecules_sum_to_zero(smiles, model):
    mol, _ = from_smiles(smiles)
    q = partial_charges(mol, model)
    assert q.shape == (len(mol),)
    assert abs(q.sum()) < 1e-6
    for i, symbol in enumerate(mol.symbols):
        if symbol in ("O", "N"):
            assert q[i] < 0
        if symbol == "H" and mol.atoms[mol.neighbors(i)[0]].symbol == "O":
            assert q[i] > 0


@pytest.mark.parametrize("model", [ChargeModel.GASTEIGER, ChargeModel.MMFF94])
@pytest.mark.parametrize(("smiles", "net"),
                         [("CC(=O)[O-]", -1), ("C[N+](C)(C)C", 1), ("[NH3+]CC(=O)[O-]", 0)])
def test_ions_sum_to_their_charge(smiles, net, model):
    mol, _ = from_smiles(smiles)
    assert formal_charges(mol).sum() == net
    assert partial_charges(mol, model).sum() == pytest.approx(net, abs=1e-6)


def test_formal_charges_from_bonding():
    mol, _ = from_smiles("C[N+](=O)[O-]")
    charges = formal_charges(mol)
    assert sorted(charges.tolist()) == [-1, 0, 0, 0, 0, 0, 1]
    assert charges[mol.symbols.index("N")] == 1


def test_water_dipole():
    mol, _ = from_smiles("O")
    dipole = dipole_moment(mol, partial_charges(mol, ChargeModel.GASTEIGER))
    assert 0.5 < np.linalg.norm(dipole) < 3.0
    oxygen = mol.symbols.index("O")
    hydrogens = mol.coords[[i for i, s in enumerate(mol.symbols) if s == "H"]].mean(axis=0)
    assert dipole @ (hydrogens - mol.coords[oxygen]) > 0   # points from O (negative) towards the hydrogens


def test_symmetric_molecule_has_no_dipole():
    mol, _ = from_smiles("c1ccccc1")
    assert np.linalg.norm(dipole_moment(mol, partial_charges(mol, ChargeModel.MMFF94))) < 0.05


def test_unsupported_elements_are_reported():
    mol, _ = from_smiles("C[Se]C")
    with pytest.raises(ChemError, match="Se"):
        partial_charges(mol, ChargeModel.GASTEIGER)
    with pytest.raises(ValueError):
        partial_charges(mol, ChargeModel.FILE)


def test_charges_and_dipole_from_a_calculation(tmp_path):
    mol, info = read_file(str(DATA / "water_charges.log"))
    assert info.charges is not None and np.allclose(info.charges, [-0.792, 0.396, 0.396])
    assert info.dipole is not None and np.allclose(info.dipole, [0.0, 0.0, -2.092])
    path = tmp_path / "water.miles"
    session.write(str(path), mol, info)
    _, again = session.read(str(path))
    assert again.charges is not None and np.allclose(again.charges, info.charges)
    assert again.dipole is not None and np.allclose(again.dipole, info.dipole)


def test_files_without_charges():
    _, info = read_file(str(DATA / "water_freq.log"))
    assert info.charges is None and info.dipole is None
