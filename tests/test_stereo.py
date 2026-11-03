import time

import numpy as np
import pytest

from miles.core.chem import from_smiles
from miles.core.molecule import Molecule
from miles.core.stereo import (
    BondStereo,
    Chirality,
    double_bond_stereo,
    flip_double_bond,
    invert_center,
    stereocenters,
)

# Expected labels are worked out by hand from the SMILES. For an atom written
# X[C@H](A)B the neighbours in order are X, H, A, B; "@" means A, H's
# successors, run anticlockwise when viewed from X. Swapping two neighbours
# swaps @ and @@, and viewing with H away from the eye reverses the sense.


def molecule(smiles: str) -> Molecule:
    return from_smiles(smiles)[0]


def only_center(mol: Molecule) -> Chirality:
    labels = stereocenters(mol)
    assert len(labels) == 1
    return next(iter(labels.values()))


def bond_lengths(mol: Molecule) -> np.ndarray:
    return np.array([np.linalg.norm(mol.coords[i] - mol.coords[j]) for i, j in mol.bonds])


def test_alanine_enantiomers():
    # N[C@@H](C)C(=O)O: [N, H, Me, COOH] @@ -> [H, N, COOH, Me] @@, so with H
    # away N -> COOH -> Me (priorities 1, 2, 3) runs anticlockwise: S.
    assert only_center(molecule("N[C@@H](C)C(=O)O")) == Chirality.S
    assert only_center(molecule("N[C@H](C)C(=O)O")) == Chirality.R


def test_bromochlorofluoromethane():
    # F[C@H](Cl)Br: [F, H, Cl, Br] @ -> [H, F, Cl, Br] @@; with H away
    # F -> Cl -> Br runs anticlockwise, so Br -> Cl -> F runs clockwise: R.
    assert only_center(molecule("F[C@H](Cl)Br")) == Chirality.R
    assert only_center(molecule("F[C@@H](Cl)Br")) == Chirality.S


def test_glyceraldehyde_ranks_by_first_point_of_difference():
    # OH (O) > CHO (O, O, H) > CH2OH (O, H, H) > H; the SMILES is the R form.
    assert only_center(molecule("O=C[C@H](O)CO")) == Chirality.R


def test_butene_e_and_z():
    for smiles, expected in (("C/C=C/C", BondStereo.E), ("C/C=C\\C", BondStereo.Z)):
        labels = double_bond_stereo(molecule(smiles))
        assert list(labels.values()) == [expected]


def test_halogen_priority_decides_e_z():
    # Br outranks Cl, CH3 outranks H; Br/C=C/C puts Br and CH3 on opposite sides.
    assert list(double_bond_stereo(molecule("Br/C(Cl)=C/C")).values()) == [BondStereo.E]
    assert list(double_bond_stereo(molecule("Br/C(Cl)=C\\C")).values()) == [BondStereo.Z]


def test_priority_is_not_size():
    # CH2OH (O, H, H) outranks the larger isopropyl (C, C, H).
    assert list(double_bond_stereo(molecule("OC/C(=C/C)C(C)C")).values()) == [BondStereo.E]
    assert list(double_bond_stereo(molecule("OC/C(=C\\C)C(C)C")).values()) == [BondStereo.Z]


def test_symmetric_double_bonds_and_centres_are_not_labelled():
    assert double_bond_stereo(molecule("CC(C)=CC")) == {}
    assert double_bond_stereo(molecule("C1=CCCCC1")) == {}
    assert stereocenters(molecule("CC(C)(Cl)Br")) == {}
    assert stereocenters(molecule("CCC(C)CC")) == {}


def test_invert_center_swaps_r_and_s_and_keeps_bonds():
    mol = molecule("O=C[C@H](O)CO")
    before = bond_lengths(mol)
    center = next(iter(stereocenters(mol)))
    assert invert_center(mol, center)
    assert stereocenters(mol) == {center: Chirality.S}
    assert np.allclose(bond_lengths(mol), before)


def test_flip_double_bond_swaps_e_and_z():
    mol = molecule("C/C=C/C")
    before = bond_lengths(mol)
    (i, j), _ = next(iter(double_bond_stereo(mol).items()))
    assert flip_double_bond(mol, i, j)
    assert double_bond_stereo(mol) == {(i, j): BondStereo.Z}
    assert np.allclose(bond_lengths(mol), before)


def test_flip_refuses_ring_and_single_bonds():
    ring = molecule("C1=CCCCC1")
    i, j = next(key for key, order in ring.bonds.items() if order == 2.0)
    assert not flip_double_bond(ring, i, j)
    chain = molecule("CCCC")
    assert not flip_double_bond(chain, *next(iter(chain.bonds)))


def test_ring_stereocentre_of_limonene():
    # Ring carbon C4: isopropenyl (C, C, C) > the ring CH2 towards the double
    # bond > the other ring CH2 > H. [C3, H, C5, C7] @@ -> [H, C3, C5, C7] @,
    # so with H away C3 -> C5 -> C7 runs clockwise, i.e. C7 -> C3 -> C5: R.
    mol = molecule("CC1=CC[C@@H](CC1)C(=C)C")
    center = only_center(mol)
    assert center == Chirality.R
    i = next(iter(stereocenters(mol)))
    ring_atoms = mol.component(next(n for n in mol.neighbors(i) if mol.degree(n) == 4 and
                                    mol.atoms[n].symbol == "C"), blocked={i})
    before = mol.coords[sorted(ring_atoms)].copy()
    assert invert_center(mol, i)
    assert stereocenters(mol) == {i: Chirality.S}
    assert np.allclose(mol.coords[sorted(ring_atoms)], before)


def test_centre_with_all_branches_in_rings_cannot_be_inverted():
    spiro = molecule("C1CCC2(C1)CCCC2")
    center = next(i for i, a in enumerate(spiro.atoms) if a.symbol == "C" and spiro.degree(i) == 4
                  and not spiro.hydrogens_on(i))
    assert not invert_center(spiro, center)


def test_peptide_alpha_carbons_and_speed():
    residues = "AAGAAAAAGAAAAAAAAAAA"
    smiles = "".join("NCC(=O)" if r == "G" else "N[C@@H](C)C(=O)" for r in residues) + "O"
    peptide = molecule(smiles)
    alpha = [i for i, a in enumerate(peptide.atoms) if a.symbol == "C" and len(peptide.hydrogens_on(i)) == 1
             and any(peptide.atoms[n].symbol == "N" for n in peptide.neighbors(i))]
    labels = stereocenters(peptide)
    assert len(alpha) == residues.count("A")
    assert {labels[i] for i in alpha} == {Chirality.S}
    assert set(labels) == set(alpha)

    protein = Molecule()
    while len(protein) < 5000:
        piece = peptide.copy()
        piece.coords = piece.coords + np.array([40.0 * (len(protein) // len(peptide)), 0.0, 0.0])
        protein.extend(piece)
    start = time.perf_counter()
    labels = stereocenters(protein)
    double_bond_stereo(protein)
    elapsed = time.perf_counter() - start
    assert len(labels) == residues.count("A") * (len(protein) // len(peptide))
    assert elapsed < 1.5


@pytest.mark.parametrize("smiles", ["C[C@H](O)[C@@H](C)Cl", "O[C@H]1CCCC[C@@H]1C"])
def test_invert_each_centre_independently(smiles):
    mol = molecule(smiles)
    labels = stereocenters(mol)
    for center, label in labels.items():
        work = mol.copy()
        assert invert_center(work, center)
        changed = stereocenters(work)
        assert changed[center] != label
        assert all(changed[k] == v for k, v in labels.items() if k != center)
