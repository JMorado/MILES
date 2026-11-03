"""Selecting atoms by SMARTS pattern."""

import pytest

from miles.builtin_plugins.select_atoms.examples import FUNCTIONAL_GROUPS
from miles.core.chem import from_smiles
from miles.core.patterns import PatternError, compile_pattern, match_pattern


def heavy_only(mol):
    return mol.subset([i for i, s in enumerate(mol.symbols) if s != "H"])


def symbols(mol, atoms):
    return sorted(mol.symbols[i] for i in atoms)


def test_hydroxyl_with_and_without_hydrogens():
    mol, _ = from_smiles("CCO")
    assert symbols(mol, match_pattern(mol, "[OX2H]").atoms) == ["O"]
    bare = heavy_only(mol)   # no hydrogens: they are implied, as for most PDB files
    assert match_pattern(bare, "[OX2H]").atoms == [2]


def test_charged_groups_follow_from_the_bonds():
    mol, _ = from_smiles("[NH3+]CC(=O)[O-]")
    assert symbols(mol, match_pattern(mol, "[NX4+]").atoms) == ["N"]
    assert symbols(mol, match_pattern(mol, "[O-]").atoms) == ["O"]


def test_aromatic_atoms_in_kekule_and_aromatic_input():
    for smiles in ("c1ccccc1O", "C1=CC=CC=C1O"):
        mol, _ = from_smiles(smiles)
        assert symbols(mol, match_pattern(mol, "a").atoms) == ["C"] * 6


def test_matches_are_listed_once_and_map_to_our_atoms():
    mol, _ = from_smiles("OC(=O)CCC(=O)O")   # succinic acid
    found = match_pattern(mol, "[CX3](=O)[OX2H1]")
    assert len(found.matches) == 2
    for match in found.matches:
        assert [mol.symbols[i] for i in match] == ["C", "O", "O"]
    assert len(found.atoms) == 6


@pytest.mark.parametrize(("smiles", "preset", "count"), [
    ("CC(=O)Oc1ccccc1C(=O)O", "Carboxylic acid or carboxylate", 1),
    ("CC(=O)Oc1ccccc1C(=O)O", "Hydroxyl (alcohol, phenol)", 0),   # the acid OH is not a hydroxyl group
    ("CC(=O)Nc1ccc(O)cc1", "Amide", 1),
    ("CC(=O)Nc1ccc(O)cc1", "Hydroxyl (alcohol, phenol)", 1),
    ("CCN(CC)CC", "Amine nitrogen", 1),
    ("CCN(CC)CC", "Hydrogen-bond acceptors", 1),
    ("Clc1ccc(Br)cc1", "Halogens", 2),
    ("c1ccncc1", "Hydrogen-bond acceptors", 1),
])
def test_presets(smiles, preset, count):
    mol, _ = from_smiles(smiles)
    smarts = dict(FUNCTIONAL_GROUPS)[preset]
    assert len(match_pattern(mol, smarts).matches) == count


def test_every_preset_is_valid():
    for _, smarts in FUNCTIONAL_GROUPS:
        compile_pattern(smarts)


@pytest.mark.parametrize("text", ["[C", "", "   ", "C(("])
def test_invalid_patterns_raise(text):
    mol, _ = from_smiles("CCO")
    with pytest.raises(PatternError):
        match_pattern(mol, text)
