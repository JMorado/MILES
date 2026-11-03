"""Backbone dihedrals and Ramachandran regions."""

import numpy as np
import pytest

from miles.builtin_plugins.ramachandran.angles import (
    Region,
    ResidueKind,
    backbone_residues,
    classify,
    ramachandran,
)
from miles.builtin_plugins.sequence.biopolymer import (
    HELIX_PHI_PSI,
    STRAND_PHI_PSI,
    SequenceKind,
    build_sequence,
)


@pytest.mark.parametrize(("kind", "expected"), [(SequenceKind.HELIX, HELIX_PHI_PSI),
                                                (SequenceKind.STRAND, STRAND_PHI_PSI)])
def test_built_peptides_have_the_requested_angles(kind, expected):
    angles = ramachandran(build_sequence("AAAAAA", kind))
    assert len(angles) == 6
    assert angles[0].phi is None and angles[0].psi is not None
    assert angles[-1].psi is None and angles[-1].phi is not None
    for residue in angles[1:-1]:
        assert residue.phi == pytest.approx(expected[0], abs=1.0)
        assert residue.psi == pytest.approx(expected[1], abs=1.0)
        assert residue.region == Region.FAVOURED
    assert angles[0].region is None


def test_residue_kinds():
    kinds = [r.kind for r in backbone_residues(build_sequence("AGAPA", SequenceKind.HELIX))]
    assert kinds == [ResidueKind.GENERAL, ResidueKind.GLYCINE, ResidueKind.PRE_PROLINE,
                     ResidueKind.PROLINE, ResidueKind.GENERAL]


def test_neighbours_come_from_peptide_bonds_not_numbers():
    mol = build_sequence("AAAAAA", SequenceKind.HELIX)
    for atom in mol.atoms:   # a numbering gap does not break the chain
        if atom.resseq > 3:
            atom.resseq += 20
    assert all(r.phi is not None for r in ramachandran(mol)[1:])
    c3 = next(i for i, a in enumerate(mol.atoms) if a.name == "C" and a.resseq == 3)
    n4 = next(i for i, a in enumerate(mol.atoms) if a.name == "N" and a.resseq == 24)
    mol.remove_bond(c3, n4)
    later = [i for i, a in enumerate(mol.atoms) if a.resseq > 3]
    mol.set_positions(later, mol.coords[later] + np.array([12.0, 0.0, 0.0]))
    angles = ramachandran(mol)
    assert angles[2].psi is None and angles[3].phi is None
    assert angles[1].psi is not None and angles[4].phi is not None


def test_regions():
    assert classify(-63, -42) == Region.FAVOURED
    assert classify(-120, 130) == Region.FAVOURED
    assert classify(60, -120) == Region.OUTLIER
    assert classify(70, 40) != Region.OUTLIER
    assert classify(-120, -175) != Region.OUTLIER   # across the +-180 boundary
    assert classify(80, 30, ResidueKind.GLYCINE) == Region.FAVOURED
    assert classify(60, 140, ResidueKind.PROLINE) == Region.OUTLIER
    assert classify(-65, 145, ResidueKind.PROLINE) == Region.FAVOURED


def test_structures_without_protein():
    from miles.core.chem import from_smiles

    mol, _ = from_smiles("CCO")
    assert ramachandran(mol) == []
