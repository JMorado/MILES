import numpy as np
import pytest

from miles.builtin_plugins.shared.mutation import mutate_residue
from miles.core.building import add_hydrogens
from miles.core.geometry import pairwise_distances
from miles.core.molecule import Molecule
from miles.files.formats import parse_text

PEPTIDE = """\
ATOM      1  N   ALA A   1      -0.677  -1.230  -0.491  1.00  0.00           N
ATOM      2  CA  ALA A   1      -0.001   0.064  -0.491  1.00  0.00           C
ATOM      3  C   ALA A   1       1.499  -0.110  -0.491  1.00  0.00           C
ATOM      4  O   ALA A   1       2.030  -1.227  -0.502  1.00  0.00           O
ATOM      5  CB  ALA A   1      -0.509   0.856   0.727  1.00  0.00           C
ATOM      6  N   ALA A   2       2.250   0.989  -0.479  1.00  0.00           N
ATOM      7  CA  ALA A   2       3.700   0.900  -0.480  1.00  0.00           C
ATOM      8  C   ALA A   2       4.300   2.300  -0.470  1.00  0.00           C
ATOM      9  O   ALA A   2       3.600   3.300  -0.460  1.00  0.00           O
ATOM     10  CB  ALA A   2       4.200   0.100  -1.680  1.00  0.00           C
ATOM     11  N   ALA A   3       5.630   2.380  -0.470  1.00  0.00           N
ATOM     12  CA  ALA A   3       6.330   3.660  -0.460  1.00  0.00           C
ATOM     13  C   ALA A   3       7.830   3.480  -0.460  1.00  0.00           C
ATOM     14  O   ALA A   3       8.360   2.360  -0.470  1.00  0.00           O
ATOM     15  CB  ALA A   3       5.830   4.460   0.760  1.00  0.00           C
HETATM   16  O   HOH A 101      10.000   0.000   0.000  1.00  0.00           O
END
"""

TRP_SIDE_CHAIN = {"CB", "CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2"}


def peptide() -> Molecule:
    return parse_text(PEPTIDE)


def residue(mol: Molecule, resseq: int) -> dict[str, int]:
    return {a.name: i for i, a in enumerate(mol.atoms) if a.resseq == resseq and a.resname != "HOH"}


def backbone(mol: Molecule, resseq: int) -> dict[str, np.ndarray]:
    names = residue(mol, resseq)
    return {n: mol.coords[names[n]].copy() for n in ("N", "CA", "C", "O")}


def closest_non_bonded(mol: Molecule) -> float:
    d = pairwise_distances(mol.coords, mol.coords)
    np.fill_diagonal(d, np.inf)
    for i, j in mol.bonds:
        d[i, j] = d[j, i] = np.inf
    for i in range(len(mol)):          # 1-3 pairs share a neighbour
        for a in mol.neighbors(i):
            for b in mol.neighbors(i):
                d[a, b] = np.inf
    return float(d.min())


def test_alanine_to_tryptophan_in_the_middle():
    mol = peptide()
    before = backbone(mol, 2)
    added = mutate_residue(mol, residue(mol, 2)["CA"], "TRP")
    names = residue(mol, 2)
    assert {mol.atoms[i].name for i in added} == TRP_SIDE_CHAIN
    assert set(names) == TRP_SIDE_CHAIN | {"N", "CA", "C", "O"}
    assert all(mol.atoms[i].resname == "TRP" for i in names.values())
    assert mol.bond_order(names["CA"], names["CB"]) == 1.0
    for name, position in before.items():
        assert np.allclose(mol.coords[names[name]], position)
    assert mol.formula() == "C17N4O4"
    assert closest_non_bonded(mol) > 1.0


@pytest.mark.parametrize("resseq", [1, 3])
def test_mutating_a_terminal_residue(resseq):
    mol = peptide()
    before = backbone(mol, resseq)
    mutate_residue(mol, residue(mol, resseq)["CA"], "PHE")
    names = residue(mol, resseq)
    assert {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"} <= names.keys()
    for name, position in before.items():
        assert np.allclose(mol.coords[names[name]], position)
    assert closest_non_bonded(mol) > 1.0


def test_glycine_has_no_side_chain():
    mol = peptide()
    assert mutate_residue(mol, residue(mol, 2)["CA"], "GLY") == []
    assert set(residue(mol, 2)) == {"N", "CA", "C", "O"}
    assert mol.formula() == "C8N3O4"


def test_proline_closes_its_ring_on_the_backbone_nitrogen():
    mol = peptide()
    mutate_residue(mol, residue(mol, 2)["CA"], "PRO")
    names = residue(mol, 2)
    assert mol.bond_order(names["N"], names["CD"]) == 1.0
    assert mol.bond_order(names["CB"], names["CG"]) == 1.0
    assert sorted(mol.atoms[j].name for j in mol.neighbors(names["N"])) == ["C", "CA", "CD"]


def test_hydrogens_follow_the_residue():
    mol = peptide()
    add_hydrogens(mol)
    mutate_residue(mol, residue(mol, 2)["CA"], "PHE")
    names = [mol.atoms[i] for i in range(len(mol)) if mol.atoms[i].resseq == 2]
    assert sum(a.symbol == "H" and a.name.startswith("HD") for a in names) == 2
    ca = residue(mol, 2)["CA"]
    assert [mol.atoms[h].name for h in mol.hydrogens_on(ca)] == ["HA"]
    assert closest_non_bonded(mol) > 1.0

    mutate_residue(mol, ca, "GLY")
    ca = residue(mol, 2)["CA"]
    assert sorted(mol.atoms[h].name for h in mol.hydrogens_on(ca)) == ["HA2", "HA3"]

    mutate_residue(mol, ca, "PRO")
    n = residue(mol, 2)["N"]
    assert mol.hydrogens_on(n) == []


def test_only_amino_acids_can_be_mutated():
    mol = peptide()
    water = next(i for i, a in enumerate(mol.atoms) if a.resname == "HOH")
    with pytest.raises(ValueError):
        mutate_residue(mol, water, "TRP")
