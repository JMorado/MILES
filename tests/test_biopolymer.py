import numpy as np
import pytest

from miles.builtin_plugins.sequence.biopolymer import SequenceKind, build_sequence, clean_sequence
from miles.core.geometry import dihedral
from miles.core.neighbors import pairs_within

ALL_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def alpha_carbons(mol):
    return mol.coords[[i for i, a in enumerate(mol.atoms) if a.name == "CA"]]


def clashes(mol):
    """Non-bonded atoms that do not share a neighbour, closer than 1.5 A (H) or 2.2 A."""
    near = set(mol.bonds)
    for i, j in mol.bonds:
        near.update((min(i, k), max(i, k)) for k in mol.neighbors(j))
        near.update((min(j, k), max(j, k)) for k in mol.neighbors(i))
    found = []
    for i, j in pairs_within(mol.coords, 2.2).tolist():
        limit = 1.5 if "H" in (mol.atoms[i].symbol, mol.atoms[j].symbol) else 2.2
        if (i, j) not in near and np.linalg.norm(mol.coords[i] - mol.coords[j]) < limit:
            found.append((i, j))
    return found


def backbone_torsions(mol, k):
    idx = {(a.resseq, a.name): i for i, a in enumerate(mol.atoms)}
    phi = dihedral(*mol.coords[[idx[(k - 1, "C")], idx[(k, "N")], idx[(k, "CA")], idx[(k, "C")]]])
    psi = dihedral(*mol.coords[[idx[(k, "N")], idx[(k, "CA")], idx[(k, "C")], idx[(k + 1, "N")]]])
    omega = dihedral(*mol.coords[[idx[(k, "CA")], idx[(k, "C")], idx[(k + 1, "N")], idx[(k + 1, "CA")]]])
    return phi, psi, omega


def test_clean_sequence_strips_headers_and_checks_letters():
    assert clean_sequence(">sp|P1|test\nacd ef\n12 gh-*", SequenceKind.HELIX) == "ACDEFGH"
    assert clean_sequence("acgt", SequenceKind.DNA) == "ACGT"
    with pytest.raises(ValueError, match="U"):
        clean_sequence("ACGU", SequenceKind.DNA)
    with pytest.raises(ValueError, match="B, J"):
        clean_sequence("ABJ", SequenceKind.STRAND)
    with pytest.raises(ValueError, match="empty"):
        clean_sequence(" > header only\n", SequenceKind.RNA)
    with pytest.raises(ValueError, match="at most"):
        clean_sequence("A" * 5000, SequenceKind.HELIX)


def test_helix_of_all_amino_acids():
    mol = build_sequence(ALL_AMINO_ACIDS, SequenceKind.HELIX)
    residues = [a.resname for a in mol.atoms if a.name == "CA"]
    assert len(residues) == 20 and residues[0] == "ALA" and residues[-1] == "TYR"
    assert {a.chain for a in mol.atoms} == {"A"}
    assert {a.ss for a in mol.atoms} == {"H"}
    ca = alpha_carbons(mol)
    spacing = np.linalg.norm(ca[3:] - ca[:-3], axis=1)
    assert np.all((spacing > 5.0) & (spacing < 5.5))
    phi, psi, omega = backbone_torsions(mol, 10)
    assert phi == pytest.approx(-57, abs=1) and psi == pytest.approx(-47, abs=1)
    assert abs(omega) == pytest.approx(180, abs=1)
    assert clashes(mol) == []
    assert mol.formula() == "C107H159N29O30S2"   # neutral 20-mer with NH2 and COOH termini


def test_strand_of_all_amino_acids():
    mol = build_sequence(ALL_AMINO_ACIDS, SequenceKind.STRAND)
    assert {a.ss for a in mol.atoms} == {"E"}
    ca = alpha_carbons(mol)
    spacing = np.linalg.norm(ca[2:] - ca[:-2], axis=1)
    assert np.all((spacing > 6.4) & (spacing < 7.0))
    phi, psi, _ = backbone_torsions(mol, 12)
    assert phi == pytest.approx(-120, abs=1) and psi == pytest.approx(130, abs=1)
    assert clashes(mol) == []


def test_proline_and_glycine_backbones():
    mol = build_sequence("GPGP", SequenceKind.HELIX)
    prolines = [i for i, a in enumerate(mol.atoms) if a.resname == "PRO" and a.name == "N"]
    for n in prolines:
        bonded = {mol.atoms[j].name for j in mol.neighbors(n)}
        assert "CD" in bonded and "H" not in bonded
    glycine_ha = [a.name for a in mol.atoms if a.resname == "GLY" and a.name.startswith("HA")]
    assert sorted(glycine_ha) == ["HA2", "HA2", "HA3", "HA3"]


@pytest.mark.parametrize(("kind", "sequence", "chains", "residues", "names"), [
    (SequenceKind.DNA, "ACGT", {"A", "B"}, 8, {"DA", "DC", "DG", "DT"}),
    (SequenceKind.DNA_SINGLE, "ACGT", {"A"}, 4, {"DA", "DC", "DG", "DT"}),
    (SequenceKind.RNA, "ACGU", {"A"}, 4, {"A", "C", "G", "U"}),
])
def test_nucleic_acids(kind, sequence, chains, residues, names):
    mol = build_sequence(sequence, kind)
    assert {a.chain for a in mol.atoms} == chains
    assert len({(a.chain, a.resseq) for a in mol.atoms}) == residues
    assert {a.resname for a in mol.atoms} == names
    assert any(a.symbol == "H" for a in mol.atoms)
    assert clashes(mol) == []
    links = [np.linalg.norm(mol.coords[i] - mol.coords[j]) for i, j in mol.bonds
             if {mol.atoms[i].name, mol.atoms[j].name} == {"O3'", "P"}]
    assert links and all(1.5 < d < 1.7 for d in links)
