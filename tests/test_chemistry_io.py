"""Formal charges, isotopes, residues and multi-structure files through reading and writing."""

import itertools
import math

import numpy as np
import pytest
from openbabel import pybel

from miles.core.building import add_hydrogens
from miles.core.charges import formal_charges
from miles.core.chem import from_smiles, to_smiles
from miles.core.geometry import TEMPLATES, free_sites, kabsch, unit_vector
from miles.core.molecule import Molecule
from miles.core.neighbors import PointGrid, cross_pairs, near_mask, pairs_within
from miles.core.perception import infer_atom_type, infer_atom_types
from miles.core.selection import find, residue_ids
from miles.files.formats import read_file, write_file
from miles.files.openbabel import read_text
from miles.files.pdb import pdb_text

PEPTIDE = """\
ATOM      1  N   GLY B 100      10.000  10.000  10.000  1.00 11.00           N
ATOM      2  CA  GLY B 100      11.450  10.000  10.000  1.00 12.00           C
ATOM      3  C   GLY B 100      12.000  11.400  10.000  1.00 13.00           C
ATOM      4  O   GLY B 100      11.300  12.400  10.000  1.00 14.00           O
ATOM      5  N   ALA B 101      13.330  11.500  10.000  1.00 15.00           N
ATOM      6  CA  ALA B 101      14.000  12.800  10.000  1.00 16.00           C
ATOM      7  C   ALA B 101      15.500  12.700  10.000  1.00 17.00           C
ATOM      8  O   ALA B 101      16.100  11.600  10.000  1.00 18.00           O
ATOM      9  CB  ALA B 101      13.500  14.000  10.500  1.00 19.00           C
HETATM   10  O   HOH W 501      20.000  20.000  20.000  1.00 30.00           O
END
"""

INSERTION = """\
ATOM      1  N   GLY A  52      10.000  10.000  10.000  1.00 11.00           N
ATOM      2  CA  GLY A  52      11.450  10.000  10.000  1.00 12.00           C
ATOM      3  C   GLY A  52      12.000  11.400  10.000  1.00 13.00           C
ATOM      4  O   GLY A  52      11.300  12.400  10.000  1.00 14.00           O
ATOM      5  N   GLY A  52A     13.330  11.500  10.000  0.50 55.00           N
ATOM      6  CA AGLY A  52A     14.000  12.800  10.000  0.50 56.00           C
ATOM      7  CA BGLY A  52A     14.100  12.900  10.100  0.50 57.00           C
ATOM      8  C   GLY A  52A     15.500  12.700  10.000  0.50 58.00           C
ATOM      9  O   GLY A  52A     16.100  11.600  10.000  0.50 59.00           O
END
"""

IONS = """\
HETATM    1 ZN    ZN A 301      10.000  10.000  10.000  1.00 20.00          ZN2+
HETATM    2 NA    NA A 302      15.000  10.000  10.000  1.00 20.00          NA
HETATM    3 CL    CL A 303      20.000  10.000  10.000  1.00 20.00          CL
HETATM    4 CA    CA A 304      25.000  10.000  10.000  1.00 20.00          CA
HETATM    5  O   HOH A 401      30.000  10.000  10.000  1.00 20.00           O
END
"""


def _columns(mol: Molecule, *names: str) -> list[tuple]:
    return list(zip(*(mol.column(name).tolist() for name in names), strict=True))


def _sdf(tmp_path, smiles: str):
    path = tmp_path / "in.sdf"
    m = pybel.readstring("smi", smiles)
    m.addh()
    m.make3D()
    m.write("sdf", str(path), overwrite=True)
    return path


# -- formal charges and isotopes ------------------------------------------------------------------

@pytest.mark.parametrize(("smiles", "net"), [("CC(=O)[O-]", -1), ("C[NH3+]", 1)])
def test_formal_charges_survive_writing_an_sdf_file(tmp_path, smiles, net):
    mol, info = read_file(str(_sdf(tmp_path, smiles)))
    assert formal_charges(mol).sum() == net
    assert sorted(set(mol.column("formal_charge").tolist())) == sorted({0, net})   # stated by the file
    out = tmp_path / "out.sdf"
    write_file(str(out), mol, info)
    assert any(line.startswith("M  CHG") for line in out.read_text(encoding="utf-8").splitlines())
    back, _ = read_file(str(out))
    assert formal_charges(back).sum() == net


SMILES = ["[Na+].[Cl-]", "C[N+](C)(C)C", "CC(=O)[O-]", "C[13CH3]", "[2H]C([2H])([2H])O"]


@pytest.mark.parametrize("smiles", SMILES)
def test_smiles_keep_charges_and_isotopes(smiles):
    mol, _ = from_smiles(smiles)
    expected = pybel.readstring("smi", smiles).write("can").split()[0]
    assert to_smiles(mol) == expected


def test_isotopes_are_read_from_smiles():
    mol, _ = from_smiles("C[13CH3]")
    assert sorted(mol.column("isotope").tolist()) == [0] * 7 + [13]


def test_charges_are_worked_out_only_where_the_bonds_can_tell(tmp_path):
    built, _ = from_smiles("C[NH3+]")
    built.set_column("formal_charge", None)   # forget the stated charges: they are worked out again
    assert formal_charges(built).sum() == 1
    serine = tmp_path / "ser.pdb"
    serine.write_text("""\
ATOM      1  N   SER A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  SER A   1      11.450  10.000  10.000  1.00  0.00           C
ATOM      3  C   SER A   1      12.000  11.400  10.000  1.00  0.00           C
ATOM      4  O   SER A   1      11.300  12.400  10.000  1.00  0.00           O
ATOM      5  CB  SER A   1      12.000   9.000   9.000  1.00  0.00           C
ATOM      6  OG  SER A   1      13.400   9.000   9.000  1.00  0.00           O
HETATM    7 NA    NA A 301      20.000  20.000  20.000  1.00  0.00          NA
END
""", encoding="utf-8")
    mol, _ = read_file(str(serine))
    charges = formal_charges(mol).tolist()
    assert charges[:6] == [0] * 6        # no hydrogens: a hydroxyl oxygen is not taken for an anion
    assert charges[6] == 1               # an isolated sodium is a cation all the same


def test_pdb_charge_columns_are_stated_charges(tmp_path):
    path = tmp_path / "ions.pdb"
    path.write_text(IONS, encoding="utf-8")
    mol, _ = read_file(str(path))
    assert mol.atoms[0].formal_charge == 2       # "2+" in columns 79-80
    assert mol.atoms[1].formal_charge is None    # blank: not stated


# -- PDB files --------------------------------------------------------------------------------------

def test_writing_pdb_keeps_chains_residue_numbers_and_b_factors(tmp_path):
    path = tmp_path / "in.pdb"
    path.write_text(PEPTIDE, encoding="utf-8")
    mol, info = read_file(str(path))
    out = tmp_path / "out.pdb"
    write_file(str(out), mol, info)
    back, _ = read_file(str(out))
    names = ("chain", "resseq", "resname", "name", "hetero")
    assert _columns(back, *names) == _columns(mol, *names)
    assert back.column("b_factor").tolist() == [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 30.0]
    assert np.allclose(back.coords, mol.coords)


def test_insertion_codes_and_alternate_locations(tmp_path):
    path = tmp_path / "ins.pdb"
    path.write_text(INSERTION, encoding="utf-8")
    mol, info = read_file(str(path))
    assert mol.column("icode").tolist() == [""] * 4 + ["A"] * 5
    assert mol.column("altloc").tolist() == [""] * 5 + ["A", "B", "", ""]
    assert mol.column("b_factor").tolist() == [11.0, 12.0, 13.0, 14.0, 55.0, 56.0, 57.0, 58.0, 59.0]
    assert mol.column("occupancy").tolist() == [1.0] * 4 + [0.5] * 5
    # 52 and 52A are two residues everywhere residues are grouped
    assert mol.residue_atoms(0) == [0, 1, 2, 3]
    assert mol.residue_key(4) == ("A", 52, "A", "GLY")
    assert len(set(residue_ids(mol).tolist())) == 2
    assert find(mol, "byres index 5").tolist() == [4, 5, 6, 7, 8]
    out = tmp_path / "out.pdb"
    write_file(str(out), mol, info)
    back, _ = read_file(str(out))
    assert _columns(back, "resseq", "icode", "altloc") == _columns(mol, "resseq", "icode", "altloc")
    assert back.column("b_factor").tolist() == mol.column("b_factor").tolist()


def test_pdb_ligand_bonds_and_ions_are_written(tmp_path):
    ligand, _ = from_smiles("CC(=O)O")
    ligand.set_column("resname", "ACT")
    ligand.set_column("hetero", True)
    text = pdb_text(ligand)
    assert "CONECT" in text
    back = read_text(text, "pdb")
    assert sorted(back.bond_orders.tolist()) == sorted(ligand.bond_orders.tolist())   # C=O kept by repeating
    path = tmp_path / "ions.pdb"
    path.write_text(IONS, encoding="utf-8")
    ions, _ = read_file(str(path))
    assert [line[76:80] for line in pdb_text(ions).splitlines() if line.startswith("HETATM")][:2] == ["ZN2+",
                                                                                                    "NA1+"]


def test_built_molecules_get_residue_and_atom_names_in_pdb_files():
    mol, _ = from_smiles("O")
    lines = [line for line in pdb_text(mol).splitlines() if line.startswith("HETATM")]
    assert [line[12:16] for line in lines] == [" O1 ", " H1 ", " H2 "]
    assert {line[17:26] for line in lines} == {"UNL     1"}


# -- files with several structures -----------------------------------------------------------------

def test_different_molecules_of_one_file_are_not_frames(tmp_path):
    xyz = tmp_path / "two.xyz"
    xyz.write_text("3\nwater\nO 0 0 0\nH 0.957 0 0\nH -0.24 0.927 0\n"
                   "3\nhydrogen cyanide\nH 0 0 -1.066\nC 0 0 0\nN 0 0 1.156\n", encoding="utf-8")
    mol, info = read_file(str(xyz))
    assert mol.symbols == ["O", "H", "H"] and info.frames is None
    assert any("1 more structure" in note for note in info.notes)
    library = tmp_path / "lib.sdf"
    out = pybel.Outputfile("sdf", str(library), overwrite=True)
    for smiles in ("CCO", "COC"):   # isomers: the same atoms, other bonds
        m = pybel.readstring("smi", smiles)
        m.addh()
        m.make3D()
        out.write(m)
    out.close()
    mol, info = read_file(str(library))
    assert info.frames is None and info.notes


def test_frames_of_one_molecule_are_still_read(tmp_path):
    xyz = tmp_path / "traj.xyz"
    steps = (f"3\nstep {k}\nO 0 0 0\nH {0.95 + 0.01 * k} 0 0\nH -0.24 0.93 0\n" for k in range(3))
    xyz.write_text("".join(steps), encoding="utf-8")
    _, info = read_file(str(xyz))
    assert info.frames is not None and info.frames.shape == (3, 3, 3)
    assert not info.notes


# -- hydrogens ----------------------------------------------------------------------------------------

def test_ions_and_metals_from_files_get_no_hydrogens(tmp_path):
    path = tmp_path / "ions.pdb"
    path.write_text(IONS, encoding="utf-8")
    mol, _ = read_file(str(path))
    assert add_hydrogens(mol) == 2
    assert [len(mol.hydrogens_on(i)) for i in range(5)] == [0, 0, 0, 0, 2]
    heme_iron = Molecule()
    heme_iron.add_atom("Fe", (0, 0, 0))
    for k, direction in enumerate(np.eye(3).tolist() + (-np.eye(3)).tolist()[:2]):
        heme_iron.add_atom("N", np.array(direction) * 2.0)
        heme_iron.add_bond(0, k + 1)
    infer_atom_types(heme_iron, only_missing=False)   # typed as if read from a file
    add_hydrogens(heme_iron, [0])
    assert heme_iron.hydrogens_on(0) == []


def test_hydrogens_are_placed_the_same_with_the_grid():
    mol, _ = from_smiles("CC(=O)Nc1ccc(O)cc1")
    heavy = mol.subset([i for i in range(len(mol)) if mol.atoms[i].symbol != "H"])
    for i in range(len(heavy)):
        infer_atom_type(heavy, i, complete=True)
    grid = PointGrid(heavy.coords, 8.0)
    from miles.core.building import missing_hydrogens, site_directions
    for i in range(len(heavy)):
        args = (heavy, i, heavy.neighbors(i), heavy.atoms[i].geom, missing_hydrogens(heavy, i))
        plain, gridded = site_directions(*args), site_directions(*args, around=grid)
        assert len(plain) == len(gridded)
        assert all(np.allclose(a, b) for a, b in zip(plain, gridded, strict=True))


def test_fitting_template_sites_matches_a_fit_per_assignment():
    rng = np.random.default_rng(3)
    for geom in ("tetrahedral", "trigonal", "octahedral", "tbp", "square_planar"):
        template = TEMPLATES[geom]
        for k in (2, 3):
            if k >= len(template):
                continue
            dirs = [unit_vector(v) for v in rng.normal(size=(k, 3))]
            best_err, best = math.inf, None
            for combo in itertools.permutations(range(len(template)), k):
                P = template[list(combo)]
                R = kabsch(P, np.array(dirs))
                err = float(((P @ R.T - np.array(dirs)) ** 2).sum())
                if err < best_err - 1e-9:
                    best_err, best = err, [R @ template[i] for i in range(len(template)) if i not in combo]
            assert np.allclose(np.array(free_sites(geom, dirs)), np.array(best))


# -- the molecule model ---------------------------------------------------------------------------------

def test_typing_all_atoms_at_once_matches_typing_each(tmp_path):
    path = tmp_path / "ions.pdb"
    path.write_text(IONS + PEPTIDE, encoding="utf-8")
    smiles = ("CC(=O)Oc1ccccc1C(=O)O", "C=C=C", "O=C=O", "CC#N", "NC(=N)N", "C1CC1")
    cases = [from_smiles(s)[0] for s in smiles]
    cases.append(read_file(str(path))[0])
    for mol in cases:
        each, together = mol.copy(), mol.copy()
        for i in range(len(each)):
            infer_atom_type(each, i)
        infer_atom_types(together, only_missing=False)
        for name in ("geom", "nsites", "valence", "inferred"):
            assert together.column(name).tolist() == each.column(name).tolist()


def test_columns_handed_out_stay_as_they_were():
    mol, _ = from_smiles("CCO")
    before = mol.coords
    kept = before.copy()
    mol.set_positions([0], [9.0, 9.0, 9.0])
    assert np.array_equal(before, kept)
    names = mol.column("name")
    mol.set_column("name", "X", [1])
    assert names.tolist() == [""] * len(mol)


def test_atom_identifiers_cannot_be_set():
    mol, _ = from_smiles("CO")
    with pytest.raises(KeyError):
        mol.atoms[0].set("uid", int(mol.uids[1]))
    with pytest.raises(KeyError):
        mol.set_column("uid", 5)


def test_bonds_from_arrays_and_neighbour_arrays():
    mol = Molecule.from_columns(np.arange(15.0).reshape(5, 3), symbol=["C"] * 5)
    mol.set_bond_arrays([[1, 0], [1, 2], [0, 1], [3, 3], [2, 4]], [1.0, 2.0, 1.5, 1.0, 1.0])
    # each pair where it first appears, with its last order; the self pair is dropped
    assert mol.bonds.items() == [((0, 1), 1.5), ((1, 2), 2.0), ((2, 4), 1.0)]
    start, neighbours = mol.adjacency()
    assert [neighbours[start[i]:start[i + 1]].tolist() for i in range(5)] == [[1], [0, 2], [1, 4], [], [2]]
    copy = mol.copy()
    copy.add_bond(3, 4)
    assert copy.neighbors(4) == [2, 3] and mol.neighbors(4) == [2]
    assert copy.bond_order(3, 4) == 1.0 and mol.bond_order(3, 4) == 0.0


def test_neighbour_searches_include_the_cutoff():
    a, b = np.zeros((1, 3)), np.array([[2.5, 0.0, 0.0]])
    assert cross_pairs(a, b, 2.5)[0].tolist() == [0]
    assert pairs_within(np.vstack([a, b]), 2.5).tolist() == [[0, 1]]
    assert near_mask(b, a, 2.5).tolist() == [True]
    grid = PointGrid(a, 2.5)
    grid.add(b[0])
    assert grid.near(b[0]).tolist() == [0, 1]
